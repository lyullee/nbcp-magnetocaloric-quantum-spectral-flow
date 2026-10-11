"""N=18 compact-supercell shape check at the two decisive fields."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'results'/'n18_compact_shape_raw'
FIELDS=[1.15,1.20]
MODELS=['M0','M3']
TEMPERATURES=[.09,.10,.12,.15]
SETTINGS={'random_vectors':2,'krylov_dim':30,'exact_low_states':8,'seed':20261011}
TRANSLATIONS=((4,1),(-2,4))


def calculate(graph,model,field):
    path=CACHE/f'{model}_B{field:.3f}T.json'
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    start=time.perf_counter()
    result=run_ftlm(model_operators(graph,model,field),TEMPERATURES,18,**SETTINGS)
    result.update({'model':model,'field_T':field,'boundary':'tin_foil',
                   'cluster_translations':TRANSLATIONS,
                   'wall_seconds_including_operator_build':time.perf_counter()-start})
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f"saved compact {model} B={field:.3f} T, wall={result['wall_seconds_including_operator_build']:.1f} s",flush=True)
    return result


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    graph=lattice(*TRANSLATIONS)
    data={(model,field):calculate(graph,model,field)
          for model in MODELS for field in FIELDS}
    comparison=[]
    for ti,temperature in enumerate(TEMPERATURES):
        row={'temperature_K':temperature}
        for model in MODELS:
            ranking=sorted([(field,data[(model,field)]['rows'][ti]['entropy_kB_per_spin'])
                            for field in FIELDS],key=lambda x:x[1],reverse=True)
            row[f'{model}_ranking']=ranking
        comparison.append(row)
    payload={'scope':'N=18 compact-shape two-field pilot','settings':SETTINGS,
             'cluster_translations':TRANSLATIONS,'comparison':comparison}
    path=ROOT/'results'/'n18_shape_pilot.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
