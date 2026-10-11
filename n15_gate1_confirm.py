"""High-accuracy N=15 check of the apparent 90 mK shift sign."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'results'/'n15_confirm_raw'
FIELDS=[1.1,1.15,1.2]
MODELS=['M0','M3']
TEMPERATURES=[.09,.10,.12,.15]
SETTINGS={'random_vectors':8,'krylov_dim':50,'exact_low_states':16,'seed':20261011}


def calculate(graph,model,field):
    path=CACHE/f'{model}_B{field:.3f}T.json'
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    start=time.perf_counter()
    result=run_ftlm(model_operators(graph,model,field),TEMPERATURES,15,**SETTINGS)
    result.update({'model':model,'field_T':field,'boundary':'tin_foil',
                   'wall_seconds_including_operator_build':time.perf_counter()-start})
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f"saved {model} B={field:.3f} T, wall={result['wall_seconds_including_operator_build']:.1f} s",flush=True)
    return result


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    graph=lattice((4,1),(1,4))
    data={(model,field):calculate(graph,model,field)
          for model in MODELS for field in FIELDS}
    rows=[]
    comparison=[]
    for ti,temperature in enumerate(TEMPERATURES):
        winners={}
        for model in MODELS:
            values=[]
            for field in FIELDS:
                entropy=data[(model,field)]['rows'][ti]['entropy_kB_per_spin']
                values.append((field,entropy))
                rows.append({'N':15,'model':model,'temperature_K':temperature,
                             'field_T':field,'entropy_kB_per_spin':entropy})
            winners[model]=max(values,key=lambda item:item[1])
        comparison.append({'temperature_K':temperature,
                           'M0_B_T':winners['M0'][0],'M3_B_T':winners['M3'][0],
                           'delta_B_T':winners['M3'][0]-winners['M0'][0],
                           'M0_ranking':sorted([(r['field_T'],r['entropy_kB_per_spin']) for r in rows
                                               if r['temperature_K']==temperature and r['model']=='M0'],
                                              key=lambda x:x[1],reverse=True),
                           'M3_ranking':sorted([(r['field_T'],r['entropy_kB_per_spin']) for r in rows
                                               if r['temperature_K']==temperature and r['model']=='M3'],
                                              key=lambda x:x[1],reverse=True)})
    payload={'scope':'N=15 high-accuracy three-field confirmation',
             'settings':SETTINGS,'rows':rows,'comparison':comparison}
    path=ROOT/'results'/'n15_gate1_confirmation.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps({'comparison':comparison},indent=2))


if __name__=='__main__':
    main()
