import numpy as np

from ewald_cluster import EwaldSettings, cluster_cell, ewald_dipole_bonds, m3_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import lattice


def main():
    graph=lattice((3,0),(0,3))
    cell=cluster_cell(graph)
    assert abs(np.linalg.det(cell))>1.
    settings=[EwaldSettings(alpha_per_angstrom=a) for a in (.16,.22,.30)]
    bond_sets=[ewald_dipole_bonds(graph,s) for s in settings]
    assert all(len(b)==graph['N']*(graph['N']-1)//2 for b in bond_sets)
    # Ewald splitting parameter must not change any finite-cell tensor.
    for left,right in zip(bond_sets[0],bond_sets[1]):
        assert (left.i,left.j)==(right.i,right.j)
        assert np.max(np.abs(left.coupling_K-right.coupling_K))<2e-10
    for left,right in zip(bond_sets[2],bond_sets[1]):
        assert np.max(np.abs(left.coupling_K-right.coupling_K))<2e-10
    # One-layer AA symmetry cancels xz and yz, so spin-up parity is valid.
    assert max(abs(b.coupling_K[0,2]) for b in bond_sets[1])<1e-12
    assert max(abs(b.coupling_K[1,2]) for b in bond_sets[1])<1e-12
    operator=SpinHamiltonianOperator(graph['N'],m3_bonds(graph),parity=0)
    rng=np.random.default_rng(91)
    x=rng.normal(size=operator.shape[0])+1j*rng.normal(size=operator.shape[0])
    y=rng.normal(size=operator.shape[0])+1j*rng.normal(size=operator.shape[0])
    assert abs(np.vdot(x,operator@y)-np.vdot(operator@x,y))<2e-9
    print('finite-cluster Ewald mapping validation passed')


if __name__=='__main__':
    main()
