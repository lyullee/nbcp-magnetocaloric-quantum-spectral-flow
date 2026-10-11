"""Resumable N=18 three-field size-crosscheck for Gate 1."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import MODEL_NAMES, model_operators
from screen_model import lattice


ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'results'/'n18_pilot_raw'
FIELDS=[1.1,1.2,1.3]
TEMPERATURES=[.09,.10,.12,.15,.20,.25,.30]
MODELS=['M0','M2','M3']
SETTINGS={'random_vectors':2,'krylov_dim':30,'exact_low_states':8,'seed':20261011}


def calculate(graph,model,field):
    path=CACHE/f'{model}_B{field:.3f}T.json'
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    start=time.perf_counter()
    result=run_ftlm(model_operators(graph,model,field),TEMPERATURES,18,**SETTINGS)
    result.update({'model':model,'model_name':MODEL_NAMES[model],
                   'field_T':field,'boundary':'tin_foil',
                   'wall_seconds_including_operator_build':time.perf_counter()-start,
                   'scope':'N=18 pilot; sparse field grid and low stochastic count'})
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f"saved {model} B={field:.3f} T, wall={result['wall_seconds_including_operator_build']:.1f} s",flush=True)
    return result


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    graph=lattice((3,-3),(3,3))
    data={(model,field):calculate(graph,model,field)
          for model in MODELS for field in FIELDS}
    rows=[]
    best=[]
    for model in MODELS:
        for ti,temperature in enumerate(TEMPERATURES):
            candidates=[]
            for field in FIELDS:
                entropy=data[(model,field)]['rows'][ti]['entropy_kB_per_spin']
                row={'N':18,'model':model,'temperature_K':temperature,
                     'field_T':field,'entropy_kB_per_spin':entropy}
                rows.append(row); candidates.append(row)
            winner=max(candidates,key=lambda row:row['entropy_kB_per_spin'])
            best.append({'model':model,'temperature_K':temperature,
                         'B_max_entropy_sparse_T':winner['field_T'],
                         'S_max_kB_per_spin':winner['entropy_kB_per_spin'],
                         'edge_warning':winner['field_T'] in (FIELDS[0],FIELDS[-1])})
    comparisons=[]
    for temperature in TEMPERATURES:
        base=next(row for row in best if row['model']=='M0' and row['temperature_K']==temperature)
        for model in MODELS[1:]:
            full=next(row for row in best if row['model']==model and row['temperature_K']==temperature)
            comparisons.append({'temperature_K':temperature,'comparison':f'M0_to_{model}',
                                'delta_B_max_entropy_sparse_T':
                                    full['B_max_entropy_sparse_T']-base['B_max_entropy_sparse_T'],
                                'M0_B_T':base['B_max_entropy_sparse_T'],
                                'full_model_B_T':full['B_max_entropy_sparse_T'],
                                'edge_warning':base['edge_warning'] or full['edge_warning']})
    payload={'scope':'N=18 sparse triage; not final B* or W95',
             'cluster_translations':graph['translations'],'fields_T':FIELDS,
             'temperatures_K':TEMPERATURES,'settings':SETTINGS,
             'rows':rows,'best':best,'comparisons':comparisons}
    path=ROOT/'results'/'n18_gate1_pilot.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps({'comparisons':comparisons},indent=2))


if __name__=='__main__':
    main()
