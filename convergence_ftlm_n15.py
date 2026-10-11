"""Numerical convergence check at the N=15 pilot point."""
from __future__ import annotations

import json
from pathlib import Path

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT=Path(__file__).resolve().parent
TEMPERATURES=[.09,.10,.12,.15,.20,.25,.30,2.0]


def main():
    graph=lattice((4,1),(1,4))
    field=1.2
    settings=[
        {'label':'pilot','random_vectors':4,'krylov_dim':40,'exact_low_states':8,'seed':20261011},
        {'label':'krylov50','random_vectors':4,'krylov_dim':50,'exact_low_states':8,'seed':20261011},
        {'label':'low12','random_vectors':4,'krylov_dim':50,'exact_low_states':12,'seed':20261011},
        {'label':'random8','random_vectors':8,'krylov_dim':50,'exact_low_states':8,'seed':20261011},
    ]
    runs=[]
    for setting in settings:
        result=run_ftlm(model_operators(graph,'M2',field),TEMPERATURES,15,
                        random_vectors=setting['random_vectors'],
                        krylov_dim=setting['krylov_dim'],
                        exact_low_states=setting['exact_low_states'],
                        seed=setting['seed'])
        runs.append({'setting':setting,'elapsed_seconds':result['elapsed_seconds'],
                     'rows':result['rows']})
        print(setting['label'],result['elapsed_seconds'],flush=True)
    reference=runs[2]
    comparisons=[]
    for run in runs:
        differences=[row['entropy_kB_per_spin']-ref['entropy_kB_per_spin']
                     for row,ref in zip(run['rows'],reference['rows'])]
        comparisons.append({'label':run['setting']['label'],
                            'entropy_difference_vs_low12':differences,
                            'max_abs_difference':max(map(abs,differences))})
    payload={'N':15,'model':'M2','field_T':field,'temperatures_K':TEMPERATURES,
             'runs':runs,'comparisons':comparisons}
    path=ROOT/'results'/'ftlm_n15_convergence.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps({'comparisons':comparisons},indent=2))


if __name__=='__main__':
    main()
