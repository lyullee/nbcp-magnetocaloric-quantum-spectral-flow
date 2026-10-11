"""Matrix-free spin-1/2 Hamiltonian used by the FTLM/TPQ stage.

The operator accepts arbitrary real Cartesian 3x3 pair-coupling matrices.
Parity reduction is optional and is rejected automatically if a coupling has
an odd number of transverse spin operators.  This matters because full 3D
dipolar tensors can contain xz/yz terms even though the 2D model preserves
spin-up parity.
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator


@dataclass(frozen=True)
class SpinBond:
    i: int
    j: int
    coupling_K: np.ndarray


def _local_coefficient(axis: int, bit: np.ndarray):
    if axis == 0:  # Sx
        return .5
    if axis == 1:  # Sy in basis |down>, |up>
        return .5j*(2*bit-1)
    if axis == 2:  # Sz
        return bit-.5
    raise ValueError(axis)


class SpinHamiltonianOperator(LinearOperator):
    def __init__(self, n: int, bonds, longitudinal_field_K: float = 0., parity=None,
                 precompute_target_limit_bytes: int = 256 << 20):
        self.n=int(n)
        self.bonds=list(bonds)
        self.parity=parity
        if parity not in (None,0,1):
            raise ValueError('parity must be None, 0, or 1')
        for b in self.bonds:
            if b.i==b.j or not (0<=b.i<n and 0<=b.j<n):
                raise ValueError('bond endpoints must be distinct valid sites')
            j=np.asarray(b.coupling_K)
            if j.shape!=(3,3) or np.max(np.abs(j-j.conj().T))>1e-12:
                raise ValueError('coupling matrix must be 3x3 Hermitian')
            if parity is not None:
                for a in range(3):
                    for c in range(3):
                        if abs(j[a,c])>1e-15 and ((a<2)+(c<2))%2:
                            raise ValueError('requested parity is broken by xz/yz pair coupling')
        all_states=np.arange(1<<n,dtype=np.uint64)
        if parity is None:
            self.states=all_states
            self.state_index=None
        else:
            # np.bit_count is not available on every supported NumPy build.
            keep=np.fromiter((int(s).bit_count()%2==parity for s in all_states),bool,len(all_states))
            self.states=all_states[keep]
            self.state_index=np.full(1<<n,-1,dtype=np.int32)
            self.state_index[self.states]=np.arange(len(self.states),dtype=np.int32)
        self.bits=np.vstack([((self.states>>i)&1).astype(np.int8) for i in range(n)])
        diag=-float(longitudinal_field_K)*np.sum(self.bits-.5,axis=0)
        # Sum separate physical contributions on the same unordered pair.
        # This changes no matrix element and substantially reduces M3 work,
        # where exchange and Ewald tensors coexist on nearest-neighbour pairs.
        combined={}
        self.off_terms=[]
        for bond in self.bonds:
            i,j_site=bond.i,bond.j
            tensor=np.asarray(bond.coupling_K,dtype=complex)
            if i>j_site:
                i,j_site=j_site,i
                tensor=tensor.T
            combined[(i,j_site)]=combined.get((i,j_site),np.zeros((3,3),complex))+tensor
        self.off_actions=[]
        target_bytes_per_action=len(self.states)*np.dtype(np.int32).itemsize
        allow_precompute=target_bytes_per_action*len(combined)<=precompute_target_limit_bytes
        for (i,j_site),j in combined.items():
            diag += np.real_if_close(j[2,2]*(self.bits[i]-.5)*(self.bits[j_site]-.5)).real
            grouped={}
            for a in range(3):
                for c in range(3):
                    if (a,c)==(2,2) or abs(j[a,c])<1e-15:
                        continue
                    self.off_terms.append((i,j_site,a,c,j[a,c]))
                    flip=((1<<i) if a<2 else 0)^((1<<j_site) if c<2 else 0)
                    grouped.setdefault(flip,[]).append((a,c,j[a,c]))
            for flip,components in grouped.items():
                lookup=np.empty(4,dtype=np.complex128)
                for local in range(4):
                    bit_i=local&1
                    bit_j=(local>>1)&1
                    lookup[local]=sum(value*_local_coefficient(a,bit_i)*
                                      _local_coefficient(c,bit_j)
                                      for a,c,value in components)
                target=None
                if allow_precompute:
                    target_state=self.states^np.uint64(flip)
                    mapped=target_state if self.state_index is None else self.state_index[target_state]
                    target=np.asarray(mapped,dtype=np.int32)
                self.off_actions.append((i,j_site,flip,lookup,target))
        self.diagonal=diag.astype(float)
        self.precomputed_target_bytes=int(sum(0 if action[4] is None else action[4].nbytes
                                              for action in self.off_actions))
        super().__init__(dtype=np.dtype(np.complex128),shape=(len(self.states),len(self.states)))

    def _matvec(self,x):
        x=np.asarray(x,dtype=complex)
        y=self.diagonal*x
        for i,j,flip,lookup,target in self.off_actions:
            if target is None:
                target_state=self.states^np.uint64(flip)
                target=(target_state if self.state_index is None else self.state_index[target_state])
            coefficient=lookup[self.bits[i]+2*self.bits[j]]
            # XOR is a permutation for each operator term, so target indices
            # are unique and advanced-index += is safe here.
            y[target] += coefficient*x
        return y

    def _rmatvec(self,x):
        return self._matvec(x)


def weak_model_bonds(graph, model):
    """Reproduce the three models in screen_weak_terms.py as Cartesian bonds."""
    from screen_weak_terms import (DELTA_XY_OVER_ZZ, JPP_EXTRA_MEV, JZZ_MEV,
                                   KD_MEV, GXY, GZ, MEV_TO_K,
                                   directional_edges)
    jxy=JZZ_MEV*DELTA_XY_OVER_ZZ
    jzz=JZZ_MEV
    jpp=0. if model=='XXZ_2026_parameters' else JPP_EXTRA_MEV
    if model=='XXZ_plus_nn_dipoles_plus_additional_Jpp':
        jxy -= KD_MEV*GXY**2/2
        jzz += KD_MEV*GZ**2
        jpp -= 3*KD_MEV*GXY**2/4
    bonds=[]
    for i,j,theta in directional_edges(graph):
        phi=2*theta
        coupling=np.diag([jxy,jxy,jzz]).astype(float)
        coupling[0,0]+=2*jpp*np.cos(phi)
        coupling[1,1]-=2*jpp*np.cos(phi)
        coupling[0,1]+=-2*jpp*np.sin(phi)
        coupling[1,0]+=-2*jpp*np.sin(phi)
        bonds.append(SpinBond(i,j,coupling*MEV_TO_K))
    return bonds


def benchmark_n15():
    from screen_model import lattice
    graph=lattice((4,1),(1,4))
    bonds=weak_model_bonds(graph,'XXZ_plus_nn_dipoles_plus_additional_Jpp')
    op=SpinHamiltonianOperator(15,bonds,parity=0)
    rng=np.random.default_rng(20261011)
    x=rng.normal(size=op.shape[0])+1j*rng.normal(size=op.shape[0])
    x/=np.linalg.norm(x)
    op@x  # warm-up
    start=time.perf_counter()
    repeats=5
    for _ in range(repeats):op@x
    elapsed=(time.perf_counter()-start)/repeats
    return {'N':15,'parity':0,'dimension':op.shape[0],'bonds':len(bonds),
            'off_terms':len(op.off_terms),'mean_matvec_seconds':elapsed,
            'persistent_array_bytes':int(op.states.nbytes+op.bits.nbytes+op.diagonal.nbytes+
                                         (0 if op.state_index is None else op.state_index.nbytes)+
                                         op.precomputed_target_bytes)}


if __name__=='__main__':
    import json
    print(json.dumps(benchmark_n15(),indent=2))
