"""Compile and plot public-QMagen N=9 XTRG convergence."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results';FIGURES=ROOT/'figures'


def main():
    runs=[]
    for path in RESULTS.glob('qmagen_n9_xtrg_D*.json'):
        data=json.loads(path.read_text())
        cold=data['cold_validation_row']
        runs.append({'D':data['bond_dimension'],'elapsed_seconds':data['elapsed_seconds'],
                     'temperature_K':cold['temperature_K'],
                     'entropy_error_kB_per_spin':cold['entropy_error_kB_per_spin'],
                     'energy_error_scaled_per_spin':cold['energy_error_scaled_per_spin'],
                     'logZ_error_per_spin':cold['logZ_error_per_spin'],
                     'file':path.name})
    runs.sort(key=lambda r:r['D'])
    errors=np.abs([r['entropy_error_kB_per_spin'] for r in runs])
    if not np.all(np.diff(errors)<=1e-12):
        raise AssertionError('absolute cold entropy error is not monotone in D')
    if runs[-1]['D']<64 or errors[-1]>=3e-4:
        raise AssertionError('D64 cold validation target was not met')
    d48=next(r for r in runs if r['D']==48)
    d64=next(r for r in runs if r['D']==64)
    payload={
        'engine':'public QMagen commit dcf17db27682d4c505155682108c612bdd11274c',
        'geometry':'N=9 3x3, OBC x / PBC y',
        'condition':{'model':'M0 2026 parameters','field_T':1.50,
                     'cold_temperature_K':runs[-1]['temperature_K']},
        'octave':'11.3.0 with project-local signed libcurl shim',
        'compatibility_patch':'project shadow svdT.m changes the second find(...) to find(...,1); upstream clone unmodified',
        'runs':runs,
        'D48_to_D64_entropy_change_kB_per_spin':abs(d64['entropy_error_kB_per_spin']-d48['entropy_error_kB_per_spin']),
        'D64_absolute_entropy_error_kB_per_spin':errors[-1],
        'decision':'PASS_N9_ENGINE_REGRESSION; proceed to staggered-grid heat-capacity validation, then cylinder pilot',
    }
    (RESULTS/'qmagen_n9_xtrg_convergence.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')

    d=np.array([r['D'] for r in runs]);times=np.array([r['elapsed_seconds'] for r in runs])
    fig,axes=plt.subplots(1,2,figsize=(9.6,4.0))
    axes[0].semilogy(d,errors,'o-',label='|entropy error|')
    axes[0].semilogy(d,np.abs([r['logZ_error_per_spin'] for r in runs]),'s--',label='|logZ error|')
    axes[0].set(xlabel='Maximum bond dimension D',ylabel='Absolute error per spin',
                title='QMagen N=9 at 86.8 mK')
    axes[0].grid(alpha=.25);axes[0].legend(frameon=False)
    axes[1].plot(d,times,'o-',color='#f58518')
    axes[1].set(xlabel='Maximum bond dimension D',ylabel='Wall time (s)',
                title='Local Octave runtime')
    axes[1].grid(alpha=.25)
    fig.tight_layout();fig.savefig(FIGURES/'qmagen_n9_xtrg_validation.png',dpi=220)
    plt.close(fig)
    print(json.dumps(payload,indent=2))


if __name__=='__main__':main()
