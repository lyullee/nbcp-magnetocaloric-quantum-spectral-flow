"""Focused N=12 endpoint-field decision screen for the representative terms."""
import json
import time
from pathlib import Path

from screen_model import OUT, VOLUME_FACTOR_MJ_PER_CM3, lattice, save_csv
from screen_weak_terms import spectrum, entropy
from broad_feasibility_screen import (TEMPERATURES, FIELDS, components,
                                      assemble, public_scale_by_temperature)


CACHE = OUT/'n12_decision_spectra'
CACHE.mkdir(exist_ok=True)


def cached_spectrum(blocks, label, field):
    import numpy as np
    path = CACHE/f'{label}_B{field:.4f}.npz'
    if path.exists():
        return np.load(path)['eigenvalues']
    e = spectrum(blocks, field)
    np.savez_compressed(path, eigenvalues=e)
    return e


def main():
    start = time.perf_counter()
    graph = lattice((2,2),(-2,4))
    parts = components(graph)
    models = {'base_XXZ': assemble(parts,0.,0.),
              'representative_Jpp_plus_NN_dipole': assemble(parts,1.,1.)}
    raw = {}
    rows = []
    for label, blocks in models.items():
        hot = entropy(cached_spectrum(blocks,label,4.),2.,12)
        raw[label] = {'hot':hot, 'cold':{}}
        for field in FIELDS:
            e = cached_spectrum(blocks,label,field)
            raw[label]['cold'][field] = {}
            for t in TEMPERATURES:
                s = entropy(e,t,12)
                q = t*(s-hot)*VOLUME_FACTOR_MJ_PER_CM3
                raw[label]['cold'][field][t] = {'S':s,'Q':q}
                rows.append({'N':12,'model':label,'B_end_T':field,'T_cold_K':t,
                             'S_cold_kB_per_spin':s,'S_hot_kB_per_spin':hot,
                             'Q_model_mJ_per_cm3':q})
        print(f'{label} completed; elapsed={time.perf_counter()-start:.1f}s',flush=True)
    save_csv('n12_field_decision_profiles.csv',rows)
    public = public_scale_by_temperature()
    decisions=[]
    for t in TEMPERATURES:
        base=max(FIELDS,key=lambda f:raw['base_XXZ']['cold'][f][t]['Q'])
        full=max(FIELDS,key=lambda f:raw['representative_Jpp_plus_NN_dipole']['cold'][f][t]['Q'])
        qb=raw['base_XXZ']['cold'][base][t]['Q']
        qf=raw['representative_Jpp_plus_NN_dipole']['cold'][full][t]['Q']
        decisions.append({'N':12,'T_cold_K':t,
            'base_optimal_B_on_grid_T':base,
            'perturbed_optimal_B_on_grid_T':full,
            'field_shift_T':full-base,
            'base_optimal_Q_model_mJ_per_cm3':qb,
            'perturbed_optimal_Q_model_mJ_per_cm3':qf,
            'optimized_delta_Q_model_mJ_per_cm3':qf-qb,
            'abs_optimized_delta_Q_as_percent_of_public_best_Q':100*abs(qf-qb)/public[t],
            'optimized_effect_exceeds_10_percent_screen':abs(qf-qb)>=.1*public[t]})
    save_csv('reoptimized_field_decisions_N12.csv',decisions)
    result={'field_grid_T':FIELDS,'temperatures_K':TEMPERATURES,
      'representative_scale_decisions':decisions,
      'max_reoptimized_effect':max(decisions,key=lambda r:r['abs_optimized_delta_Q_as_percent_of_public_best_Q']),
      'scope':['Classical exact diagonalization, not a QPU calculation.',
               'Field optimum is restricted to the stated coarse grid.',
               'N=12 is exact for this cluster but is not a bulk convergence result.',
               'Same model exclusions as broad_feasibility_screen.json.'],
      'elapsed_seconds':time.perf_counter()-start}
    (OUT/'n12_decision_screen.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
