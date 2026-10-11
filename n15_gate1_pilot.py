"""Resumable N=15 M0-vs-M3 coarse Gate-1 field pilot."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from ftlm import run_ftlm
from material_models import MODEL_NAMES, model_operators
from screen_model import VOLUME_FACTOR_MJ_PER_CM3, lattice


ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'results'/'n15_pilot_raw'
FIELDS=[.9,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5]
TEMPERATURES=[.09,.10,.12,.15,.20,.25,.30]
MODELS=['M0','M1','M2','M3']
SETTINGS={'random_vectors':4,'krylov_dim':40,'exact_low_states':8,'seed':20261011}


def cache_path(model,field):
    return CACHE/f'{model}_B{field:.3f}T.json'


def calculate(graph,model,field,temperatures):
    path=cache_path(model,field)
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    wall_start=time.perf_counter()
    result=run_ftlm(model_operators(graph,model,field),temperatures,15,**SETTINGS)
    result.update({'model':model,'model_name':MODEL_NAMES[model],
                   'field_T':field,'boundary':'tin_foil',
                   'wall_seconds_including_operator_build':time.perf_counter()-wall_start})
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f"saved {model} B={field:.3f} T, {result['elapsed_seconds']:.1f} s",flush=True)
    return result


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    graph=lattice((4,1),(1,4))
    data={}
    for model in MODELS:
        data[(model,4.)]=calculate(graph,model,4.,[2.])
        for field in FIELDS:
            data[(model,field)]=calculate(graph,model,field,TEMPERATURES)
    profiles=[]
    summary=[]
    for model in MODELS:
        hot=data[(model,4.)]['rows'][0]['entropy_kB_per_spin']
        for ti,temperature in enumerate(TEMPERATURES):
            values=[]
            for field in FIELDS:
                entropy=data[(model,field)]['rows'][ti]['entropy_kB_per_spin']
                q=temperature*(entropy-hot)*VOLUME_FACTOR_MJ_PER_CM3
                row={'N':15,'model':model,'temperature_K':temperature,
                     'field_T':field,'entropy_kB_per_spin':entropy,
                     'hot_entropy_2K_4T_kB_per_spin':hot,
                     'Q_entropy_budget_mJ_per_cm3':q}
                profiles.append(row)
                values.append(row)
            best=max(values,key=lambda row:row['Q_entropy_budget_mJ_per_cm3'])
            qmax=best['Q_entropy_budget_mJ_per_cm3']
            accepted=[row['field_T'] for row in values
                      if row['Q_entropy_budget_mJ_per_cm3']>=.95*qmax]
            summary.append({'N':15,'model':model,'temperature_K':temperature,
                            'B_star_coarse_T':best['field_T'],'Qmax_coarse_mJ_per_cm3':qmax,
                            'W95_sampled_min_T':min(accepted),'W95_sampled_max_T':max(accepted),
                            'grid_warning':'coarse sampled bounds; not interpolated final W95'})
    comparisons=[]
    for temperature in TEMPERATURES:
        base=next(r for r in summary if r['model']=='M0' and r['temperature_K']==temperature)
        for model in MODELS[1:]:
            full=next(r for r in summary if r['model']==model and r['temperature_K']==temperature)
            model_at_m0=next(r['Q_entropy_budget_mJ_per_cm3'] for r in profiles
                             if r['model']==model and r['temperature_K']==temperature and
                             r['field_T']==base['B_star_coarse_T'])
            regret=(full['Qmax_coarse_mJ_per_cm3']-model_at_m0)/full['Qmax_coarse_mJ_per_cm3']*100
            comparisons.append({'temperature_K':temperature,'comparison':f'M0_to_{model}',
                                'delta_B_star_coarse_T':full['B_star_coarse_T']-base['B_star_coarse_T'],
                                'M0_B_star_T':base['B_star_coarse_T'],
                                'full_model_B_star_T':full['B_star_coarse_T'],
                                'M0_field_regret_under_full_model_percent':regret})
    payload={'scope':'N=15 0.05 T FTLM pilot; finite cluster, not bulk conclusion',
             'fields_T':FIELDS,'temperatures_K':TEMPERATURES,'settings':SETTINGS,
             'profiles':profiles,'summary':summary,'comparisons':comparisons}
    path=ROOT/'results'/'n15_gate1_pilot.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps({'comparisons':comparisons},indent=2))


if __name__=='__main__':
    main()
