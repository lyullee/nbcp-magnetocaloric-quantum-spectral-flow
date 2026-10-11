"""Diagnose N=15/N=18 field maxima and low-energy coverage."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parent
N18_RAW=ROOT/'results'/'n18_pilot_raw'
MODELS=['M0','M2','M3']
FIELDS_BY_MODEL={'M0':[1.1,1.15,1.2,1.25,1.3],
                 'M2':[1.1,1.2,1.3],
                 'M3':[1.1,1.15,1.2,1.25,1.3]}
TEMPERATURES=[.09,.10,.12,.15,.20,.25,.30]


def low_only_entropy(low_by_sector,temperature,n):
    energies=np.asarray([e for sector in low_by_sector for e in sector],float)
    reference=energies.min()
    weights=np.exp(-(energies-reference)/temperature)
    mean=np.dot(weights,energies)/weights.sum()
    return float((np.log(weights.sum())+(mean-reference)/temperature)/n)


def main():
    raw={}
    diagnostics=[]
    for model in MODELS:
        for field in FIELDS_BY_MODEL[model]:
            data=json.loads((N18_RAW/f'{model}_B{field:.3f}T.json').read_text(encoding='utf-8'))
            raw[(model,field)]=data
            energies=np.sort(np.asarray([e for sector in data['low_energies_K_by_sector'] for e in sector]))
            entry={'model':model,'field_T':field,'ground_energy_K':float(energies[0]),
                   'first_gap_K':float(energies[1]-energies[0]),
                   'states_within_0p1K':int(np.sum(energies-energies[0]<=.1)),
                   'wall_seconds':data['wall_seconds_including_operator_build'],
                   'temperature_rows':[]}
            for row in data['rows']:
                low=low_only_entropy(data['low_energies_K_by_sector'],row['temperature_K'],18)
                entry['temperature_rows'].append({
                    'temperature_K':row['temperature_K'],
                    'full_entropy_kB_per_spin':row['entropy_kB_per_spin'],
                    'low16_entropy_kB_per_spin':low,
                    'complement_increment_kB_per_spin':row['entropy_kB_per_spin']-low})
            diagnostics.append(entry)
    margins=[]
    for model in MODELS:
        for ti,temperature in enumerate(TEMPERATURES):
            values=[(field,raw[(model,field)]['rows'][ti]['entropy_kB_per_spin'])
                    for field in FIELDS_BY_MODEL[model]]
            values.sort(key=lambda item:item[1],reverse=True)
            margins.append({'N':18,'model':model,'temperature_K':temperature,
                            'winner_field_T':values[0][0],
                            'winner_minus_runnerup_kB_per_spin':values[0][1]-values[1][1],
                            'ranking':values})
    n15=json.loads((ROOT/'results'/'n15_gate1_pilot.json').read_text(encoding='utf-8'))
    size_comparison=[]
    for model in MODELS:
        for temperature in TEMPERATURES:
            n15_values=[(row['field_T'],row['entropy_kB_per_spin']) for row in n15['profiles']
                        if row['model']==model and row['temperature_K']==temperature and
                        row['field_T'] in FIELDS_BY_MODEL[model]]
            n18_values=[(field,raw[(model,field)]['rows'][TEMPERATURES.index(temperature)]['entropy_kB_per_spin'])
                        for field in FIELDS_BY_MODEL[model]]
            n15_winner=max(n15_values,key=lambda item:item[1])[0]
            n18_winner=max(n18_values,key=lambda item:item[1])[0]
            size_comparison.append({'model':model,'temperature_K':temperature,
                                    'N15_winner_on_common_grid_T':n15_winner,
                                    'N18_winner_on_common_grid_T':n18_winner,
                                    'same_winner':n15_winner==n18_winner})
    payload={'diagnostics':diagnostics,'field_margins':margins,
             'size_comparison_common_grid':size_comparison,
             'interpretation':[
                 'The low16/full difference measures omitted-state contribution, not stochastic error.',
                 'A winner at 1.1 or 1.3 T is edge-censored on this sparse grid.',
                 'No Gate-1 decision is allowed until 1.15 and 1.25 T or an equivalent interpolation is checked.'
             ]}
    path=ROOT/'results'/'n15_n18_size_diagnostic.json'
    path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
    compact={'low_temperature_margins':[row for row in margins if row['temperature_K']<=.12],
             'size_comparison':[row for row in size_comparison if row['temperature_K']<=.12]}
    print(json.dumps(compact,indent=2))


if __name__=='__main__':
    main()
