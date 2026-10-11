"""First-order bridge from bulk critical fields to the public XTRG window.

The script deliberately brackets, rather than asserts, a finite-temperature
M4 result.  It transports the published M0 XTRG curve using two common local
similarity assumptions: constant reduced field B/Bc and constant distance
B-Bc.  Agreement of the two gives a decision-scale diagnostic.  A genuine
M4 finite-T tensor-network calculation is still required for a final claim.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from screen_model import MU_B_OVER_KB

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
FIGURES=ROOT/'figures'


def transform_field(source_b, source_hc, target_hc, method):
    source_b=np.asarray(source_b,dtype=float)
    if method=='reduced_field':
        return source_b*target_hc/source_hc
    if method=='critical_offset':
        return source_b+(target_hc-source_hc)
    raise ValueError(method)


def inverse_field(target_b, source_hc, target_hc, method):
    target_b=np.asarray(target_b,dtype=float)
    if method=='reduced_field':
        return target_b*source_hc/target_hc
    if method=='critical_offset':
        return target_b-(target_hc-source_hc)
    raise ValueError(method)


def main():
    bulk=json.loads((RESULTS/'bulk_polarized_stability.json').read_text(encoding='utf-8'))
    rows=bulk['rows']
    hc={}
    for r in rows:
        hc[(r['model'],r['boundary'])]=r['critical_field_T']
    # Gao et al. 2022 source Hamiltonian: Jxy=0.88 K, Jz=1.48 K, g=4.89.
    source_hc=(3*1.48+1.5*.88)/(4.89*MU_B_OVER_KB)
    profile=np.genfromtxt(RESULTS/'npj2022_xtrg_entropy_reconstructed.csv',
                          delimiter=',',names=True,encoding='utf-8-sig')
    profile=profile[profile['anchor_N']==12]
    methods=('reduced_field','critical_offset')
    targets={
        'M0':hc[('M0','tin_foil')],
        'M1':hc[('M1','tin_foil')],
        'M2':hc[('M2','tin_foil')],
        'M3_tin_foil':hc[('M3','tin_foil')],
        'M4_tin_foil':hc[('M4','tin_foil')],
        'M3_sphere_vacuum':hc[('M3','sphere_vacuum')],
        'M4_sphere_vacuum':hc[('M4','sphere_vacuum')],
    }
    metric_rows=[]
    regret_rows=[]
    robust_rows=[]
    curve_rows=[]
    for temperature in (.09,.10,.12,.15):
        p=np.sort(profile[np.isclose(profile['temperature_K'],temperature)],
                  order='field_T')
        bs=p['field_T']; q=p['Q_proxy_mJ_per_cm3']
        k=int(np.argmax(q)); qmax=float(q[k]); bstar=float(bs[k])
        accepted=bs[q>=.95*qmax]
        source_w=(float(accepted.min()),float(accepted.max()))
        for method in methods:
            for name,target_hc in targets.items():
                mapped=transform_field(np.array([bstar,*source_w]),source_hc,
                                       target_hc,method)
                metric_rows.append({
                    'temperature_K':temperature,'model_scenario':name,
                    'mapping_method':method,'source_Hc_T':source_hc,
                    'target_Hc_T':target_hc,'B_star_T':float(mapped[0]),
                    'W95_min_T':float(mapped[1]),'W95_max_T':float(mapped[2]),
                    'Qmax_mJ_per_cm3_inherited':qmax,
                })
                mapped_fields=transform_field(bs,source_hc,target_hc,method)
                for b,qq in zip(mapped_fields,q):
                    if b>=1.2:  # bridge is intended only for near-saturation window
                        curve_rows.append({'temperature_K':temperature,
                                           'model_scenario':name,
                                           'mapping_method':method,
                                           'field_T':float(b),
                                           'Q_proxy_mJ_per_cm3':float(qq)})
            # Consequence of operating M4 at the field selected by M0.
            m0_hc=targets['M0']
            m0_b=float(transform_field(bstar,source_hc,m0_hc,method))
            for boundary in ('tin_foil','sphere_vacuum'):
                m4_hc=targets[f'M4_{boundary}']
                source_eval=float(inverse_field(m0_b,source_hc,m4_hc,method))
                q_at=float(np.interp(source_eval,bs,q))
                regret=100*(qmax-q_at)/qmax
                regret_rows.append({
                    'temperature_K':temperature,'mapping_method':method,
                    'M4_boundary':boundary,'M0_selected_field_T':m0_b,
                    'M4_optimal_field_T':float(transform_field(bstar,source_hc,m4_hc,method)),
                    'delta_B_star_T':float(transform_field(bstar,source_hc,m4_hc,method)-m0_b),
                    'M4_Q_at_M0_field_mJ_per_cm3':q_at,
                    'M4_Qmax_mJ_per_cm3_inherited':qmax,
                    'regret_percent':regret,
                })

        # Minimax field over the four explicit M4 bridge scenarios.
        scenarios=[(method,boundary,targets[f'M4_{boundary}'])
                   for method in methods for boundary in ('tin_foil','sphere_vacuum')]
        grid=np.linspace(1.35,1.80,9001)
        regrets=[]
        for method,boundary,target_hc in scenarios:
            source_eval=inverse_field(grid,source_hc,target_hc,method)
            q_eval=np.interp(source_eval,bs,q,left=np.nan,right=np.nan)
            regrets.append((qmax-q_eval)/qmax)
        regrets=np.asarray(regrets)
        worst=np.nanmax(regrets,axis=0)
        valid=np.all(np.isfinite(regrets),axis=0)
        idx=np.where(valid)[0][np.argmin(worst[valid])]
        robust_rows.append({
            'temperature_K':temperature,
            'robust_field_T':float(grid[idx]),
            'worst_case_regret_percent':float(100*worst[idx]),
            'scenario_regret_percent':[float(100*x) for x in regrets[:,idx]],
            'scenario_order':['reduced_field/tin_foil','reduced_field/sphere_vacuum',
                              'critical_offset/tin_foil','critical_offset/sphere_vacuum'],
        })

    def write_csv(name,records):
        with (RESULTS/name).open('w',newline='',encoding='utf-8-sig') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader();writer.writerows(records)
    write_csv('bulk_bridge_metrics.csv',metric_rows)
    write_csv('bulk_bridge_regret.csv',regret_rows)
    write_csv('bulk_bridge_curves.csv',curve_rows)

    near90=[r for r in regret_rows if r['temperature_K']==.09]
    delta_range=[min(r['delta_B_star_T'] for r in near90),
                 max(r['delta_B_star_T'] for r in near90)]
    regret_range=[min(r['regret_percent'] for r in near90),
                  max(r['regret_percent'] for r in near90)]
    payload={
        'status':'decision-scale surrogate, not a finite-temperature M4 solution',
        'source_XTRG_model_Hc_T':source_hc,
        'mapping_assumptions':{
            'reduced_field':'preserve B/Bc',
            'critical_offset':'preserve B-Bc',
        },
        'bulk_Hc_T':targets,
        '90mK_M4_minus_M0_delta_B_star_bracket_T':delta_range,
        '90mK_regret_of_M0_field_under_M4_bracket_percent':regret_range,
        'predeclared_thresholds':{'delta_B_star_GO_T':.05,'regret_GO_percent':5.0},
        'threshold_screen':{
            'all_four_scenarios_delta_B_star_above_GO':bool(delta_range[0]>=.05),
            'all_four_scenarios_regret_above_GO':bool(regret_range[0]>=5.0),
            'interpretation':'candidate GO for full finite-T bulk calculation; Gate 1 remains HOLD until that calculation',
        },
        'metrics':metric_rows,'regrets':regret_rows,'robust_fields':robust_rows,
        'limitations':[
            'The source XTRG Hamiltonian and 2026 Hamiltonian differ in anisotropy and weak terms.',
            'Horizontal critical-field alignment does not calculate entropy redistribution or Qmax changes.',
            'The additive mapping is local to the near-saturation peak and is invalid around zero field.',
            'The Ewald Gamma term is sample-shape dependent; tin-foil and spherical-vacuum conventions bracket it.',
            'Final B*, W95, heat and regret require a converged finite-temperature M4 bulk solver.',
        ],
    }
    (RESULTS/'bulk_to_engineering_surrogate.json').write_text(
        json.dumps(payload,indent=2),encoding='utf-8')

    # Audit figure.
    FIGURES.mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,3,figsize=(14.4,4.2))
    hierarchy=['M0','M1','M2','M3_tin_foil','M4_tin_foil']
    vals=[targets[x] for x in hierarchy]
    axes[0].plot(range(len(vals)),vals,'o-',color='#4c78a8')
    axes[0].axhline(1.86,color='#e45756',ls='--',label='reported full model 1.86 T')
    axes[0].set_xticks(range(len(vals)),['M0','M1','M2\nNN dipole','M3\nfull dipole','M4\n+ interlayer'])
    axes[0].set(ylabel=r'$H_c$ (T)',title='Bulk model hierarchy')
    axes[0].grid(alpha=.25);axes[0].legend(frameon=False,fontsize=8)

    p=np.sort(profile[np.isclose(profile['temperature_K'],.09)],order='field_T')
    bs=p['field_T'];q=p['Q_proxy_mJ_per_cm3']
    for name,target_hc,color in [('M0',targets['M0'],'#4c78a8'),
                                 ('M4_tin_foil',targets['M4_tin_foil'],'#f58518'),
                                 ('M4_sphere_vacuum',targets['M4_sphere_vacuum'],'#54a24b')]:
        mapped=transform_field(bs,source_hc,target_hc,'reduced_field')
        mask=(mapped>=1.2)&(mapped<=1.9)&(q>=0)
        axes[1].plot(mapped[mask],q[mask],label=name.replace('_',' '),color=color)
    axes[1].set(xlabel='Field B (T)',ylabel=r'$Q$ proxy (mJ cm$^{-3}$)',
                title='90 mK critical-field bridge')
    axes[1].grid(alpha=.25);axes[1].legend(frameon=False,fontsize=8)

    labels=['red./tin','off./tin','red./sphere','off./sphere']
    ordered=[]
    for boundary in ('tin_foil','sphere_vacuum'):
        for method in methods:
            ordered.append(next(r for r in near90 if r['M4_boundary']==boundary and
                                r['mapping_method']==method))
    axes[2].bar(labels,[r['delta_B_star_T'] for r in ordered],color='#72b7b2')
    axes[2].axhline(.05,color='#e45756',ls='--',label='predeclared GO threshold')
    axes[2].set(ylabel=r'$\Delta B^*$ (T)',title='M4 minus M0 at 90 mK')
    axes[2].tick_params(axis='x',rotation=25);axes[2].grid(axis='y',alpha=.25)
    axes[2].legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES/'bulk_bridge_engineering_screen.png',dpi=220)
    plt.close(fig)
    print(json.dumps({k:v for k,v in payload.items() if k not in ('metrics','regrets')},indent=2))


if __name__=='__main__':main()
