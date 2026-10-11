"""Exact small-cluster sensitivity screen; NOT the full 3D material model.

2026 source: https://arxiv.org/html/2607.23173v1, Eq. 1,3,5 and Table 1.
Exchange and g values are a separate parameter set from screen_model.py.
Two ablations: additional bond anisotropy alone, then nearest-neighbor
in-plane dipoles plus that anisotropy. Full long-range dipoles, interlayer
couplings, incommensurate order and thermodynamic-limit convergence are absent.
These are sensitivity probes, not calibrated material predictions.
"""
import argparse
import json
import time
import numpy as np
from screen_model import ROOT, OUT, lattice, save_csv, MU_B_OVER_KB, VOLUME_FACTOR_MJ_PER_CM3

MEV_TO_K = 11.60451812155
JZZ_MEV = .1225
DELTA_XY_OVER_ZZ = .636
GXY, GZ = 4.2, 4.716
KD_MEV = .00036
JPP_EXTRA_MEV = -.0019

def directional_edges(graph):
    a, b = graph['translations']
    n = graph['N']
    def key(x,y):
        return ((b[1]*x-b[0]*y)%n,(-a[1]*x+a[0]*y)%n)
    indices = {key(x,y):i for i,(x,y) in enumerate(graph['coordinates'])}
    directions = [(1,0,0.),(0,1,np.pi/3),(1,-1,-np.pi/3)]
    result = {}
    for i,(x,y) in enumerate(graph['coordinates']):
        for dx,dy,angle in directions:
            j = indices[key(x+dx,y+dy)]
            edge = tuple(sorted((i,j)))
            assert edge not in result
            result[edge] = angle
    assert sorted(result) == graph['edges']
    return [(i,j,theta) for (i,j),theta in result.items()]

def matrices(graph, model):
    jxy = JZZ_MEV*DELTA_XY_OVER_ZZ
    jzz = JZZ_MEV
    jpp = 0. if model == 'XXZ_2026_parameters' else JPP_EXTRA_MEV
    if model == 'XXZ_plus_nn_dipoles_plus_additional_Jpp':
        # Eq. 3 rotated to each bond; only the in-plane nearest neighbors.
        jxy -= KD_MEV*GXY**2/2
        jzz += KD_MEV*GZ**2
        jpp -= 3*KD_MEV*GXY**2/4
    n = graph['N']
    edges = directional_edges(graph)
    blocks = []
    for parity in [0,1]:
        states = [s for s in range(2**n) if s.bit_count()%2 == parity]
        index = {s:p for p,s in enumerate(states)}
        h = np.zeros((len(states),len(states)),dtype=complex)
        mz = np.array([s.bit_count()-n/2 for s in states],dtype=float)
        for p,s in enumerate(states):
            for i,j,angle in edges:
                up_i,up_j = (s>>i)&1,(s>>j)&1
                same = up_i == up_j
                h[p,p] += jzz*(.25 if same else -.25)
                flipped = index[s^(1<<i)^(1<<j)]
                if not same:
                    h[flipped,p] += jxy/2
                else:
                    # gamma = exp(2i theta); S+ acts on a down (0) bit.
                    gamma = np.exp(2j*angle)
                    h[flipped,p] += jpp*(gamma if not up_i else gamma.conjugate())
        h *= MEV_TO_K
        assert np.max(np.abs(h-h.conj().T)) < 1e-12
        blocks.append((h,mz))
    return blocks

def spectrum(blocks,field):
    eigenvalues=[]
    for h,mz in blocks:
        hf=h.copy()
        hf[np.diag_indices_from(hf)] -= GZ*MU_B_OVER_KB*field*mz
        eigenvalues.extend(np.linalg.eigvalsh(hf))
    return np.array(eigenvalues)

def entropy(eigenvalues,temperature,n):
    e=eigenvalues-eigenvalues.min()
    w=np.exp(-e/temperature)
    z=w.sum()
    return float((np.log(z)+np.dot(w,e)/z/temperature)/n)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--sizes',nargs='+',type=int,default=[9])
    args=parser.parse_args()
    start=time.perf_counter()
    rows=[]
    models=['XXZ_2026_parameters','XXZ_plus_additional_Jpp','XXZ_plus_nn_dipoles_plus_additional_Jpp']
    for n in args.sizes:
        assert n in [9,12]
        graph=lattice((3,0),(0,3)) if n==9 else lattice((2,2),(-2,4))
        for model in models:
            blocks=matrices(graph,model)
            hot=spectrum(blocks,4.)
            assert len(hot)==2**n
            assert abs(entropy(hot,1e6,n)-np.log(2))<1e-10
            # At zero added anisotropy this parity calculation must reproduce
            # an independent magnetization-sector diagonalization, using the
            # 2026 exchange parameters rather than the Nature parameter set.
            if model=='XXZ_2026_parameters':
                import screen_model as sm
                saved=(sm.JXY,sm.JZ,sm.G)
                sm.JXY=JZZ_MEV*DELTA_XY_OVER_ZZ*MEV_TO_K
                sm.JZ=JZZ_MEV*MEV_TO_K
                sm.G=GZ
                # spectra() writes a cache: preserve old Nature-model cache.
                original_out=sm.OUT
                cache_out=OUT/'weak_term_reference'
                cache_out.mkdir(exist_ok=True)
                sm.OUT=cache_out
                reference=sm.spectra(graph)
                sm.OUT=original_out
                reference_hot=np.sort(reference['XXZ']-GZ*MU_B_OVER_KB*4.*reference['Mz'])
                assert np.max(np.abs(np.sort(hot)-reference_hot))<1e-10
                sm.JXY,sm.JZ,sm.G=saved
            s_hot=entropy(hot,2.,n)
            for field in [0.,1.3,1.5,1.66]:
                cold=spectrum(blocks,field)
                s_cold=entropy(cold,.1,n)
                rows.append({'N':n,'model':model,'B_end_T':field,'T_cold_K':.1,
                    'S_cold_kB_per_spin':s_cold,'S_hot_kB_per_spin':s_hot,
                    'Q_entropy_budget_mJ_per_cm3':.1*(s_cold-s_hot)*VOLUME_FACTOR_MJ_PER_CM3})
            print(f'N={n}, {model} completed, elapsed={time.perf_counter()-start:.2f}s',flush=True)
    save_csv('weak_term_sensitivity.csv',rows)
    comparisons=[]
    for n in args.sizes:
        for field in [0.,1.3,1.5,1.66]:
            relevant=[r for r in rows if r['N']==n and r['B_end_T']==field]
            for r in relevant[1:]:
                comparisons.append({'N':n,'B_end_T':field,'model':r['model'],
                    'delta_S_cold_kB_per_spin':r['S_cold_kB_per_spin']-relevant[0]['S_cold_kB_per_spin'],
                    'delta_Q_entropy_budget_mJ_per_cm3':r['Q_entropy_budget_mJ_per_cm3']-relevant[0]['Q_entropy_budget_mJ_per_cm3']})
    summary={'source':'https://arxiv.org/html/2607.23173v1','preprint':True,
      'parameters':{'Jzz_meV':JZZ_MEV,'Delta_XY_over_ZZ':DELTA_XY_OVER_ZZ,
                    'gxy':GXY,'gz':GZ,'Kd_nn_meV':KD_MEV,'additional_Jpp_meV':JPP_EXTRA_MEV},
      'scope':['Classical exact diagonalization, no quantum circuit or hardware.',
               '2D N=9/12 sensitivity ablations only. No long-range dipoles or interlayer exchange.',
               'Nearest-neighbor dipole terms and additional Jpp are counted separately as in source.',
               'Does not reproduce the 2026 full 3D model or its incommensurate phase.',
               'Negative endpoint budget flags a non-absorbing requested endpoint, not a device cooling measurement.'],
      'rows':rows,'comparisons':comparisons,'elapsed_seconds':time.perf_counter()-start}
    (OUT/'weak_term_sensitivity.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({'comparisons':comparisons,'elapsed_seconds':summary['elapsed_seconds']},indent=2))

if __name__=='__main__':
    main()
