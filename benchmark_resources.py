"""Measured matrix-free resource table and conservative N=18 run plan."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from ewald_cluster import m3_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator, weak_model_bonds
from screen_model import lattice


ROOT=Path(__file__).resolve().parent


def persistent_bytes(operator):
    return int(operator.states.nbytes+operator.bits.nbytes+
               operator.diagonal.nbytes+
               (0 if operator.state_index is None else operator.state_index.nbytes)+
               operator.precomputed_target_bytes)


def measure(n,graph,model,repeats=3):
    if model=='M2':
        bonds=weak_model_bonds(graph,'XXZ_plus_nn_dipoles_plus_additional_Jpp')
    elif model=='M3':
        bonds=m3_bonds(graph)
    else:
        raise ValueError(model)
    operator=SpinHamiltonianOperator(n,bonds,parity=0)
    rng=np.random.default_rng(20261011+n+(0 if model=='M2' else 100))
    vector=rng.normal(size=operator.shape[0])+1j*rng.normal(size=operator.shape[0])
    vector/=np.linalg.norm(vector)
    operator@vector
    start=time.perf_counter()
    for _ in range(repeats):
        operator@vector
    seconds=(time.perf_counter()-start)/repeats
    return {'N':n,'model':model,'parity_dimension':operator.shape[0],
            'spin_bond_objects':len(bonds),'off_diagonal_terms':len(operator.off_terms),
            'grouped_off_diagonal_actions':len(operator.off_actions),
            'mean_matvec_seconds':seconds,
            'persistent_operator_bytes':persistent_bytes(operator),
            'krylov_50_complex_bytes':int(50*operator.shape[0]*16)}


def main():
    graphs={15:lattice((4,1),(1,4)),18:lattice((3,-3),(3,3))}
    rows=[measure(n,graphs[n],model) for n in (15,18) for model in ('M2','M3')]
    baseline=json.loads((ROOT/'results'/'ftlm_validation_benchmark.json').read_text(encoding='utf-8'))
    measured_n15_ftlm=baseline['benchmark_n15']['elapsed_seconds']
    base_matvec=next(r['mean_matvec_seconds'] for r in rows if r['N']==15 and r['model']=='M2')
    estimates=[]
    for row in rows:
        estimates.append({'N':row['N'],'model':row['model'],
            'rough_seconds_per_field_from_N15_M2_ratio':
                measured_n15_ftlm*row['mean_matvec_seconds']/base_matvec,
            'status':'planning estimate; eigsh convergence and reorthogonalization can change runtime'})
    payload={'measured':rows,
             'reference_ftlm':{'N':15,'model':'M2','random_vectors_per_sector':4,
                 'krylov_dimension':50,'exact_low_states_per_sector':8,
                 'elapsed_seconds':measured_n15_ftlm},
             'planning_estimates':estimates,
             'n18_execution_plan':[
                 'Run M0/M2/M3 at B=1.10, 1.20, 1.30 T first.',
                 'Use R=4, Krylov=40, k_low=8 for pilot; compare the central point with Krylov=50 and k_low=12.',
                 'Continue to a 0.05 T grid only if the pilot has a stable entropy ordering.',
                 'Refine to 0.025 T only within the provisional W95 region.',
                 'Do not use N=18 as evidence for qz=1/3; it is a thermodynamic sensitivity cluster.'
             ]}
    path=ROOT/'results'/'classical_solver_resource_table.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
