"""Diagnose finite-cluster shell effects by critical-field alignment."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
S=.5


def cluster_critical_mode(graph,model):
    n=graph['N']
    matrix=np.zeros((2*n,2*n),dtype=float)
    for bond in model_bonds(graph,model):
        i,j=bond.i,bond.j
        tensor=np.asarray(bond.coupling_K,dtype=float)
        # Pair transverse block and the loss of longitudinal pair energy when
        # either fixed-length spin tilts away from z.
        matrix[2*i:2*i+2,2*j:2*j+2]+=tensor[:2,:2]
        matrix[2*j:2*j+2,2*i:2*i+2]+=tensor[:2,:2].T
        matrix[2*i:2*i+2,2*i:2*i+2]-=tensor[2,2]*np.eye(2)
        matrix[2*j:2*j+2,2*j:2*j+2]-=tensor[2,2]*np.eye(2)
    assert np.max(np.abs(matrix-matrix.T))<1e-11
    values,vectors=np.linalg.eigh(matrix)
    h=-S*values[0]
    return {'model':model,'critical_field_T':float(h/(GZ*MU_B_OVER_KB)),
            'minimum_transverse_eigenvalue_K':float(values[0]),
            'next_eigenvalue_K':float(values[1]),
            'soft_mode_degeneracy_1e-8':int(np.sum(np.abs(values-values[0])<1e-8))}


def aligned_rows(size, raw, hc_n, hc_bulk):
    rows=[]
    for r in raw['rows']:
        for model,key in [('M0','S_M0'),('M3','S_M3')]:
            rows.append({'N':size,'temperature_K':r['temperature_K'],
                         'model':model,'cluster_field_T':r['field_T'],
                         'aligned_bulk_field_T':r['field_T']-hc_n[model]+hc_bulk[model],
                         'entropy_kB_per_spin':r[key]})
    return rows


def main():
    shapes={
        'N9_3x3':lattice((3,0),(0,3)),
        'N12_skew':lattice((2,2),(-2,4)),
        'N15_used':lattice((4,1),(1,4)),
        'N18_used':lattice((4,1),(-2,4)),
    }
    stability=[]
    for name,graph in shapes.items():
        for model in ('M0','M3'):
            r=cluster_critical_mode(graph,model)
            r.update({'cluster':name,'N':graph['N'],
                      'translations':graph['translations']})
            stability.append(r)
            print(json.dumps(r),flush=True)
    bulk=json.loads((RESULTS/'bulk_polarized_stability.json').read_text())
    hc_bulk={m:next(r['critical_field_T'] for r in bulk['rows']
                    if r['model']==m and r['boundary']=='tin_foil')
             for m in ('M0','M3')}
    n15=json.loads((RESULTS/'n15_bulk_window_confirmation.json').read_text())
    n18=json.loads((RESULTS/'n18_bulk_window_pilot.json').read_text())
    hc_by_size={}
    for n,label in [(15,'N15_used'),(18,'N18_used')]:
        hc_by_size[n]={m:next(r['critical_field_T'] for r in stability
                              if r['cluster']==label and r['model']==m)
                       for m in ('M0','M3')}
    aligned=aligned_rows(15,n15,hc_by_size[15],hc_bulk)+aligned_rows(18,n18,hc_by_size[18],hc_bulk)
    comparison=[]
    for n in (15,18):
        for t in (.09,.10,.12):
            by_model={}
            for model in ('M0','M3'):
                z=[r for r in aligned if r['N']==n and r['temperature_K']==t and r['model']==model]
                z=sorted(z,key=lambda r:r['aligned_bulk_field_T'])
                by_model[model]=(np.array([r['aligned_bulk_field_T'] for r in z]),
                                 np.array([r['entropy_kB_per_spin'] for r in z]))
            lo=max(x[0].min() for x in by_model.values())
            hi=min(x[0].max() for x in by_model.values())
            grid=np.arange(np.ceil(lo*100)/100,np.floor(hi*100)/100+1e-9,.01)
            delta=(np.interp(grid,*by_model['M3'])-
                   np.interp(grid,*by_model['M0']))
            for b,d in zip(grid,delta):
                comparison.append({'N':n,'temperature_K':t,
                                   'aligned_bulk_field_T':float(b),
                                   'deltaS_aligned_kB_per_spin':float(d)})
    common=[]
    for t in (.09,.10,.12):
        a=[r for r in comparison if r['N']==15 and r['temperature_K']==t]
        b=[r for r in comparison if r['N']==18 and r['temperature_K']==t]
        fields=sorted(set(round(r['aligned_bulk_field_T'],8) for r in a)&
                      set(round(r['aligned_bulk_field_T'],8) for r in b))
        for field in fields:
            r15=next(r for r in a if round(r['aligned_bulk_field_T'],8)==field)
            r18=next(r for r in b if round(r['aligned_bulk_field_T'],8)==field)
            d15=r15['deltaS_aligned_kB_per_spin'];d18=r18['deltaS_aligned_kB_per_spin']
            common.append({'temperature_K':t,'aligned_bulk_field_T':field,
                           'N15_deltaS':d15,'N18_deltaS':d18,
                           'same_sign':bool(np.sign(d15)==np.sign(d18)),
                           'absolute_size_difference':abs(d15-d18)})
    payload={
        'method':'finite-cluster transverse-stability alignment to bulk Hc',
        'bulk_Hc_T':hc_bulk,'cluster_stability':stability,
        'aligned_rows':aligned,'aligned_deltaS':comparison,
        'common_size_comparison':common,
        '90mK_same_sign_fraction':float(np.mean([r['same_sign'] for r in common
                                                 if r['temperature_K']==.09])) if common else None,
        'warning':'Linear interpolation of 0.05 T entropy grids; diagnostic only.',
    }
    (RESULTS/'finite_cluster_stability_alignment.json').write_text(
        json.dumps(payload,indent=2),encoding='utf-8')
    fig,ax=plt.subplots(figsize=(6.4,4.1))
    for n,marker in [(15,'o'),(18,'s')]:
        z=[r for r in comparison if r['N']==n and r['temperature_K']==.09]
        ax.plot([r['aligned_bulk_field_T'] for r in z],
                [r['deltaS_aligned_kB_per_spin'] for r in z],
                marker=marker,label=f'N={n}')
    ax.axhline(0,color='black',lw=.8)
    ax.set(xlabel='Bulk-aligned field (T)',
           ylabel=r'$S(M3)-S(M0)$ ($k_B$/spin)',
           title='Critical-field alignment does not remove shell effect')
    ax.grid(alpha=.25);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(ROOT/'figures'/'finite_cluster_Hc_alignment.png',dpi=220)
    plt.close(fig)
    print(json.dumps({'hc_bulk':hc_bulk,'hc_by_size':hc_by_size,
                      'common_90mK':[r for r in common if r['temperature_K']==.09],
                      'same_sign_fraction':payload['90mK_same_sign_fraction']},indent=2))


if __name__=='__main__':main()
