"""N=9/12 M0-vs-M3 sequence to diagnose odd/even finite-size shifts."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'results'/'small_size_sequence_raw'
FIELDS=[1.10,1.15,1.20]
TEMPERATURES=[.09,.10,.12,.15]
MODELS=['M0','M3']
GRAPHS={9:lattice((3,0),(0,3)),12:lattice((2,2),(-2,4))}
SETTINGS={
    9:{'random_vectors':1,'krylov_dim':32,'exact_low_states':255,'seed':20261011},
    12:{'random_vectors':4,'krylov_dim':50,'exact_low_states':64,'seed':20261011},
}


def calculate(n,model,field):
    path=CACHE/f'N{n}_{model}_B{field:.3f}T.json'
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    start=time.perf_counter()
    result=run_ftlm(model_operators(GRAPHS[n],model,field),TEMPERATURES,n,**SETTINGS[n])
    result.update({'N':n,'model':model,'field_T':field,'boundary':'tin_foil',
                   'cluster_translations':GRAPHS[n]['translations'],
                   'wall_seconds_including_operator_build':time.perf_counter()-start})
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f"saved N={n} {model} B={field:.2f}, wall={result['wall_seconds_including_operator_build']:.1f}s",flush=True)
    return result


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    data={(n,model,field):calculate(n,model,field)
          for n in (9,12) for model in MODELS for field in FIELDS}
    comparisons=[]
    for n in (9,12):
        for ti,temperature in enumerate(TEMPERATURES):
            winners={}
            rankings={}
            for model in MODELS:
                ranking=sorted([(field,data[(n,model,field)]['rows'][ti]['entropy_kB_per_spin'])
                                for field in FIELDS],key=lambda x:x[1],reverse=True)
                rankings[model]=ranking
                winners[model]=ranking[0][0]
            comparisons.append({'N':n,'temperature_K':temperature,
                                'M0_B_T':winners['M0'],'M3_B_T':winners['M3'],
                                'delta_B_T':winners['M3']-winners['M0'],
                                'M0_ranking':rankings['M0'],'M3_ranking':rankings['M3']})
    payload={'scope':'small-size parity diagnostic; three-field window only',
             'fields_T':FIELDS,'settings':SETTINGS,'comparisons':comparisons}
    path=ROOT/'results'/'small_size_sequence.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
