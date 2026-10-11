"""Fine N=12 field scan around the entropy-budget maximum."""
import json
import time
import numpy as np

from screen_model import OUT, VOLUME_FACTOR_MJ_PER_CM3, lattice, save_csv
from screen_weak_terms import entropy
from broad_feasibility_screen import TEMPERATURES, components, assemble, public_scale_by_temperature
from n12_decision_screen import cached_spectrum

FIELDS = [round(float(x),3) for x in np.arange(.95,1.4001,.025)]


def main():
    start=time.perf_counter()
    parts=components(lattice((2,2),(-2,4)))
    models={'base_XXZ':assemble(parts,0.,0.),
            'representative_Jpp_plus_NN_dipole':assemble(parts,1.,1.)}
    values={}
    rows=[]
    for label,blocks in models.items():
        hot=entropy(cached_spectrum(blocks,label,4.),2.,12)
        values[label]={}
        for field in FIELDS:
            e=cached_spectrum(blocks,label,field)
            values[label][field]={}
            for t in TEMPERATURES:
                s=entropy(e,t,12)
                q=t*(s-hot)*VOLUME_FACTOR_MJ_PER_CM3
                values[label][field][t]=q
                rows.append({'N':12,'model':label,'B_end_T':field,'T_cold_K':t,
                             'Q_model_mJ_per_cm3':q})
        print(f'{label} fine grid complete; elapsed={time.perf_counter()-start:.1f}s',flush=True)
    save_csv('n12_fine_field_profiles.csv',rows)
    public=public_scale_by_temperature()
    decisions=[]
    for t in TEMPERATURES:
        item={'N':12,'T_cold_K':t}
        for label,prefix in [('base_XXZ','base'),('representative_Jpp_plus_NN_dipole','perturbed')]:
            best=max(FIELDS,key=lambda f:values[label][f][t])
            qmax=values[label][best][t]
            in95=[f for f in FIELDS if values[label][f][t]>=.95*qmax]
            item.update({f'{prefix}_optimal_B_on_fine_grid_T':best,
                         f'{prefix}_optimal_Q_model_mJ_per_cm3':qmax,
                         f'{prefix}_95pct_field_window_min_T':min(in95),
                         f'{prefix}_95pct_field_window_max_T':max(in95),
                         f'{prefix}_95pct_field_window_width_T':max(in95)-min(in95)})
        dq=item['perturbed_optimal_Q_model_mJ_per_cm3']-item['base_optimal_Q_model_mJ_per_cm3']
        item.update({'field_shift_T':item['perturbed_optimal_B_on_fine_grid_T']-item['base_optimal_B_on_fine_grid_T'],
                     'optimized_delta_Q_model_mJ_per_cm3':dq,
                     'abs_optimized_delta_Q_as_percent_of_public_best_Q':100*abs(dq)/public[t]})
        decisions.append(item)
    save_csv('n12_fine_field_decisions.csv',decisions)
    result={'field_grid_T':FIELDS,'decisions':decisions,
            'scope':['Classical exact diagonalization; no QPU.',
                     'Grid endpoints may truncate a 95% window; inspect min/max before using width.',
                     'Finite N=12 cluster; not a bulk prediction.'],
            'elapsed_seconds':time.perf_counter()-start}
    (OUT/'n12_fine_field_screen.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
