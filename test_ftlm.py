"""Regression and thermodynamic-identity tests for ftlm.py."""
import numpy as np

from ftlm import build_sector_rules, thermodynamics
from matrix_free_hamiltonian import SpinBond, SpinHamiltonianOperator


def ring_bonds(n):
    tensor=np.diag([.71,.71,1.13])
    return [SpinBond(i,(i+1)%n,tensor) for i in range(n)]


def exact_rows(operators,temperatures,n):
    values=[]
    for operator in operators:
        eye=np.eye(operator.shape[0],dtype=complex)
        dense=np.column_stack([operator@eye[:,i] for i in range(operator.shape[0])])
        values.extend(np.linalg.eigvalsh(dense))
    values=np.asarray(values)
    reference=values.min()
    result=[]
    for temperature in temperatures:
        w=np.exp(-(values-reference)/temperature)
        z=w.sum()
        e=np.dot(w,values)/z
        e2=np.dot(w,values**2)/z
        result.append((np.log(z)+(e-reference)/temperature,
                       (e2-e**2)/temperature**2))
    return result


def main():
    n=6
    operators=[SpinHamiltonianOperator(n,ring_bonds(n),.37,parity=p) for p in (0,1)]
    # Treat all but one state in each sector exactly.  The remaining
    # one-dimensional complement makes the trace estimator deterministic and
    # tests exact-state projection, sector weighting, and thermodynamic sums.
    sectors=[build_sector_rules(op,random_vectors=2,krylov_dim=32,
                                exact_low_states=op.shape[0]-1,seed=7+p)
             for p,op in enumerate(operators)]
    temperatures=np.array([.09,.2,1.,20.])
    estimated=thermodynamics(sectors,temperatures,n)
    exact=exact_rows(operators,temperatures,n)
    for row,(entropy,capacity) in zip(estimated,exact):
        assert abs(row['entropy_kB_total']-entropy)<2e-10
        assert abs(row['heat_capacity_kB_total']-capacity)<2e-9
    assert abs(estimated[-1]['entropy_kB_per_spin']-np.log(2))<2e-3
    # Direct thermodynamic identity at a temperature away from sharp finite
    # size features: dU/dT = C (k_B=1).
    t=.4
    dt=1e-5
    fine=thermodynamics(sectors,[t-dt,t,t+dt],n)
    derivative=(fine[2]['energy_K_total']-fine[0]['energy_K_total'])/(2*dt)
    assert abs(derivative-fine[1]['heat_capacity_kB_total'])<2e-7
    print('low-energy-corrected FTLM validation passed')


if __name__=='__main__':
    main()
