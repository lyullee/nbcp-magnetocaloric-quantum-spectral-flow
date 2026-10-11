"""Validate material FTLM at N=9, then benchmark one N=15 condition."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ftlm import run_ftlm
from matrix_free_hamiltonian import SpinHamiltonianOperator, weak_model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ, matrices, spectrum


ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
MODEL='XXZ_plus_nn_dipoles_plus_additional_Jpp'
TEMPERATURES=np.array([.09,.10,.12,.15,.20,.25,.30,2.0])


def material_operators(graph,field):
    bonds=weak_model_bonds(graph,MODEL)
    zeeman=GZ*MU_B_OVER_KB*field
    return [SpinHamiltonianOperator(graph['N'],bonds,zeeman,parity=p) for p in (0,1)]


def exact_entropy_rows(graph,field):
    values=spectrum(matrices(graph,MODEL),field)
    reference=values.min()
    rows=[]
    for temperature in TEMPERATURES:
        w=np.exp(-(values-reference)/temperature)
        z=w.sum()
        e=np.dot(w,values)/z
        entropy=(np.log(z)+(e-reference)/temperature)/graph['N']
        rows.append({'temperature_K':float(temperature),
                     'entropy_kB_per_spin':float(entropy)})
    return rows


def validate_n9(field):
    graph=lattice((3,0),(0,3))
    exact=exact_entropy_rows(graph,field)
    trials=[]
    for seed in (20261011,20261012,20261013):
        result=run_ftlm(material_operators(graph,field),TEMPERATURES,graph['N'],
                        random_vectors=6,krylov_dim=50,
                        exact_low_states=12,seed=seed)
        errors=[]
        for got,want in zip(result['rows'],exact):
            errors.append(got['entropy_kB_per_spin']-want['entropy_kB_per_spin'])
        trials.append({'seed':seed,'elapsed_seconds':result['elapsed_seconds'],
                       'entropy_errors_kB_per_spin':errors,
                       'max_abs_entropy_error_kB_per_spin':max(map(abs,errors))})
    return {'N':9,'field_T':field,'model':MODEL,'temperatures_K':TEMPERATURES.tolist(),
            'exact_rows':exact,'trials':trials,
            'max_abs_error_all_trials':max(t['max_abs_entropy_error_kB_per_spin'] for t in trials)}


def benchmark_n15(field):
    graph=lattice((4,1),(1,4))
    result=run_ftlm(material_operators(graph,field),TEMPERATURES,graph['N'],
                    random_vectors=4,krylov_dim=50,
                    exact_low_states=8,seed=20261011)
    result.update({'N':15,'field_T':field,'model':MODEL,
                   'cluster_translations':graph['translations']})
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--field',type=float,default=1.2)
    parser.add_argument('--skip-n15',action='store_true')
    args=parser.parse_args()
    payload={'validation_n9':validate_n9(args.field)}
    print(json.dumps({'validation_n9':payload['validation_n9']['trials']},indent=2),flush=True)
    if not args.skip_n15:
        payload['benchmark_n15']=benchmark_n15(args.field)
        print(json.dumps({'benchmark_n15':payload['benchmark_n15']},indent=2),flush=True)
    OUT.mkdir(exist_ok=True)
    path=OUT/'ftlm_validation_benchmark.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(path)


if __name__=='__main__':
    main()
