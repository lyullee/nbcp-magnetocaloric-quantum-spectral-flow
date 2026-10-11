"""Numerical convergence and reproducibility audit for bulk stability."""
import json
from pathlib import Path

import numpy as np

from bulk_polarized_stability import EwaldSettings, locate_critical_mode


def main():
    settings=[
        EwaldSettings(.18,42.,8.), EwaldSettings(.22,36.,8.),
        EwaldSettings(.22,42.,6.), EwaldSettings(.22,42.,8.),
        EwaldSettings(.22,42.,10.), EwaldSettings(.22,50.,8.),
        EwaldSettings(.28,42.,8.),
    ]
    convergence=[]
    for s in settings:
        r=locate_critical_mode('M4',s)
        convergence.append({**s.__dict__,
                            'critical_field_T':r['critical_field_T'],
                            'critical_wavevector_rlu':r['critical_wavevector_rlu']})
        print(json.dumps(convergence[-1]),flush=True)
    seed_rows=[]
    base=EwaldSettings()
    for seed in (1,260711,999983):
        r=locate_critical_mode('M4',base,seed=seed)
        seed_rows.append({'seed':seed,'critical_field_T':r['critical_field_T'],
                          'critical_wavevector_rlu':r['critical_wavevector_rlu']})
    fields=np.array([x['critical_field_T'] for x in convergence])
    seed_fields=np.array([x['critical_field_T'] for x in seed_rows])
    result={
        'model':'M4',
        'convergence_rows':convergence,
        'optimizer_seed_rows':seed_rows,
        'ewald_field_spread_T':float(np.ptp(fields)),
        'optimizer_seed_field_spread_T':float(np.ptp(seed_fields)),
        'base_field_T':convergence[3]['critical_field_T'],
        'reported_reference_T':1.86,
        'base_minus_reported_T':convergence[3]['critical_field_T']-1.86,
    }
    path=Path('results/bulk_stability_convergence.json')
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
