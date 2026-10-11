"""Compile the finite-size Gate-1 evidence into a table, report, and figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from screen_model import VOLUME_FACTOR_MJ_PER_CM3


ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
FIGURES=ROOT/'figures'


def metrics(label,n,rows_m0,rows_m3,temperature,hot_entropy):
    m0=[row for row in rows_m0 if row['temperature_K']==temperature]
    m3=[row for row in rows_m3 if row['temperature_K']==temperature]
    b0=max(m0,key=lambda row:row['entropy_kB_per_spin'])['field_T']
    best=max(m3,key=lambda row:row['entropy_kB_per_spin'])
    b3=best['field_T']
    at_b0=next(row for row in m3 if row['field_T']==b0)
    denominator=best['entropy_kB_per_spin']-hot_entropy
    regret=(best['entropy_kB_per_spin']-at_b0['entropy_kB_per_spin'])/denominator*100 if denominator>0 else float('nan')
    qmax=temperature*denominator*VOLUME_FACTOR_MJ_PER_CM3
    accepted=[row['field_T'] for row in m3
              if temperature*(row['entropy_kB_per_spin']-hot_entropy)*VOLUME_FACTOR_MJ_PER_CM3>=.95*qmax]
    return {'cluster':label,'N':n,'temperature_K':temperature,'M0_B_star_sampled_T':b0,
            'M3_B_star_sampled_T':b3,'delta_B_star_sampled_T':b3-b0,
            'M0_field_regret_under_M3_percent':regret,
            'M3_Qmax_sampled_mJ_per_cm3':qmax,
            'M3_W95_sampled_min_T':min(accepted) if accepted else None,
            'M3_W95_sampled_max_T':max(accepted) if accepted else None,
            'field_grid_T':sorted(row['field_T'] for row in m3)}


def main():
    audit=json.loads((RESULTS/'screening_summary.json').read_text(encoding='utf-8'))
    hot_mj=audit['source_data_audit']['measured_entropy_at_2K_4T_mJ_per_cm3_K']
    hot=hot_mj/VOLUME_FACTOR_MJ_PER_CM3
    small=json.loads((RESULTS/'small_size_sequence.json').read_text(encoding='utf-8'))
    n15=json.loads((RESULTS/'n15_gate1_confirmation.json').read_text(encoding='utf-8'))
    n18diag=json.loads((RESULTS/'n15_n18_size_diagnostic.json').read_text(encoding='utf-8'))
    compact=json.loads((RESULTS/'n18_shape_pilot.json').read_text(encoding='utf-8'))

    datasets=[]
    for n in (9,12):
        rows=[]
        for entry in small['comparisons']:
            if entry['N']!=n: continue
            for model in ('M0','M3'):
                for field,entropy in entry[f'{model}_ranking']:
                    rows.append({'model':model,'temperature_K':entry['temperature_K'],
                                 'field_T':field,'entropy_kB_per_spin':entropy})
        datasets.append((f'N{n}',n,[r for r in rows if r['model']=='M0'],[r for r in rows if r['model']=='M3']))
    datasets.append(('N15',15,[r for r in n15['rows'] if r['model']=='M0'],
                     [r for r in n15['rows'] if r['model']=='M3']))

    primary_rows=[]
    for model in ('M0','M3'):
        for field in (1.1,1.15,1.2,1.25,1.3):
            raw=json.loads((RESULTS/'n18_pilot_raw'/f'{model}_B{field:.3f}T.json').read_text(encoding='utf-8'))
            for row in raw['rows']:
                primary_rows.append({'model':model,'temperature_K':row['temperature_K'],
                                     'field_T':field,'entropy_kB_per_spin':row['entropy_kB_per_spin']})
    datasets.append(('N18-A',18,[r for r in primary_rows if r['model']=='M0'],
                     [r for r in primary_rows if r['model']=='M3']))

    compact_rows=[]
    for entry in compact['comparison']:
        for model in ('M0','M3'):
            for field,entropy in entry[f'{model}_ranking']:
                compact_rows.append({'model':model,'temperature_K':entry['temperature_K'],
                                     'field_T':field,'entropy_kB_per_spin':entropy})
    datasets.append(('N18-B compact',18,[r for r in compact_rows if r['model']=='M0'],
                     [r for r in compact_rows if r['model']=='M3']))

    table=[]
    for label,n,m0,m3 in datasets:
        for temperature in (.09,.10):
            table.append(metrics(label,n,m0,m3,temperature,hot))

    FIGURES.mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(11,4.2))
    for temperature,marker in ((.09,'o'),(.10,'s')):
        rows=[r for r in table if r['temperature_K']==temperature]
        x=np.arange(len(rows))+( -.08 if temperature==.09 else .08)
        axes[0].scatter(x,[r['delta_B_star_sampled_T'] for r in rows],marker=marker,
                        label=f'{temperature*1000:.0f} mK',s=55)
        axes[1].scatter(x,[r['M0_field_regret_under_M3_percent'] for r in rows],marker=marker,s=55)
    labels=[r['cluster'] for r in table if r['temperature_K']==.09]
    for ax in axes:
        ax.set_xticks(np.arange(len(labels)),labels,rotation=25,ha='right')
        ax.grid(alpha=.25)
    axes[0].axhline(0,color='black',lw=.8)
    axes[0].axhline(.05,color='tab:red',ls='--',lw=.8)
    axes[0].axhline(-.05,color='tab:red',ls='--',lw=.8)
    axes[0].set_ylabel(r'$B^*_{M3}-B^*_{M0}$ (T)')
    axes[0].legend(frameon=False)
    axes[1].axhline(5,color='tab:red',ls='--',lw=.8,label='5% gate')
    axes[1].set_ylabel('M0-field regret under M3 (%)')
    axes[1].legend(frameon=False)
    fig.suptitle('Finite-cluster Gate-1 diagnostic (sampled field grids)')
    fig.tight_layout()
    figure_path=FIGURES/'gate1_finite_size_diagnostic.png'
    fig.savefig(figure_path,dpi=220,bbox_inches='tight')
    plt.close(fig)

    payload={'public_hot_entropy':{'mJ_per_cm3_K':hot_mj,'kB_per_spin':hot},
             'table':table,'gate1_status':'HOLD',
             'gate1_reason':'Individual clusters cross thresholds, but the 90 mK shift changes from +0.05 T at N12/N15 to -0.05 T in both N18 shapes.',
             'quantum_stage_authorized':False,
             'next_required_evidence':'Thermodynamic-limit or independent bulk/semi-classical calculation before any quantum algorithm work.'}
    (RESULTS/'gate1_preliminary_report.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')

    lines=['# Gate 1 예비 판정: 보류','',
           f'- 공개 실험의 2 K, 4 T 기준 엔트로피: `{hot:.6f} kB/spin` (`{hot_mj:.6f} mJ cm^-3 K^-1`)',
           '- 개별 유한계는 `|ΔB*|=0.05 T` 또는 regret 5%를 넘지만 이동 부호가 수렴하지 않는다.',
           '- 90 mK: N=9은 0, N=12/N=15는 +0.05 T, N=18 두 형상은 -0.05 T이다.',
           '- 따라서 현재 수치로 벌크 최적장 이동이나 양자알고리즘 필요성을 주장하지 않는다.',
           '- 다음 필수 단계는 큰 반고전/벌크 계산 또는 동등한 독립 열역학 방법이다.','',
           '## 판정표','',
           '|cluster|T (mK)|M0 B* (T)|M3 B* (T)|ΔB* (T)|regret (%)|sampled W95 (T)|',
           '|---|---:|---:|---:|---:|---:|---|']
    for r in table:
        lines.append(f"|{r['cluster']}|{r['temperature_K']*1000:.0f}|{r['M0_B_star_sampled_T']:.2f}|{r['M3_B_star_sampled_T']:.2f}|{r['delta_B_star_sampled_T']:+.2f}|{r['M0_field_regret_under_M3_percent']:.2f}|{r['M3_W95_sampled_min_T']}–{r['M3_W95_sampled_max_T']}|")
    lines.extend(['','모든 B*와 W95는 표시된 유한 자기장 격자에서만 계산한 값이며 벌크 외삽값이 아니다.'])
    (ROOT/'Gate1_예비판정_20261011.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
