"""Small exact benchmark for coupling-constant free-energy integration.

F1-F0 = integral_0^1 <V>_(H0+lambda V) d lambda; Delta S=-d(Delta F)/dT.
This validates a thermodynamic observable route, NOT a quantum algorithm.
The 9-spin, in-plane nearest-neighbor sensitivity proxy omits the full 3D
interactions. State preparation and actual shot budgets remain unimplemented.
"""
import json
import time
import numpy as np
import screen_weak_terms as sw
from screen_model import OUT, lattice, MU_B_OVER_KB, save_csv

def simpson(y):
    count=len(y)-1
    assert count%2==0
    w=np.ones(count+1);w[1:-1:2]=4;w[2:-1:2]=2
    return float(np.dot(w,y)/(3*count))

def thermo(spectrum,temp,n):
    energies,vdiag=spectrum
    shifted=energies-energies.min()
    weights=np.exp(-shifted/temp)
    z=weights.sum();prob=weights/z
    return {'f':float((energies.min()-temp*np.log(z))/n),
            's':float((np.log(z)+np.dot(prob,shifted)/temp)/n),
            'v':float(np.dot(prob,vdiag)/n)}

def spectra(h0,v,mz,field,lam):
    energies=[];vdiag=[]
    for base,perturb,m in zip(h0,v,mz):
        h=base+lam*perturb
        h=h.copy()
        h[np.diag_indices_from(h)]-=sw.GZ*MU_B_OVER_KB*field*m
        e,vec=np.linalg.eigh(h)
        vd=np.real(np.sum(vec.conj()*(perturb@vec),axis=0))
        energies.extend(e);vdiag.extend(vd)
    return np.asarray(energies),np.asarray(vdiag)

def main():
    start=time.perf_counter()
    n=9
    graph=lattice((3,0),(0,3))
    base=sw.matrices(graph,'XXZ_2026_parameters')
    full=sw.matrices(graph,'XXZ_plus_nn_dipoles_plus_additional_Jpp')
    h0=[pair[0] for pair in base]
    v=[full[i][0]-h0[i] for i in range(2)]
    mz=[pair[1] for pair in base]
    results=[]
    for b in [1.3,1.5,4.]:
        samples=[spectra(h0,v,mz,b,lam) for lam in np.linspace(0,1,17)]
        for t in ([2.] if b==4. else [.08,.1,.12,.15,2.]):
            values=[thermo(s,t,n) for s in samples]
            exact_delta_f=values[-1]['f']-values[0]['f']
            exact_delta_s=values[-1]['s']-values[0]['s']
            delta_f=simpson([value['v'] for value in values])
            coarse_f=simpson([values[i]['v'] for i in range(0,17,2)])
            for dt in [.001,.0025,.005]:
                above=simpson([thermo(s,t+dt,n)['v'] for s in samples])
                below=simpson([thermo(s,t-dt,n)['v'] for s in samples])
                route_s=-(above-below)/(2*dt)
                results.append({'N':n,'B_T':b,'T_K':t,'derivative_step_K':dt,
                    'exact_delta_F_K_per_spin':exact_delta_f,
                    'integration_delta_F_K_per_spin':delta_f,
                    'F_error_17_nodes_K_per_spin':delta_f-exact_delta_f,
                    'F_error_9_nodes_K_per_spin':coarse_f-exact_delta_f,
                    'exact_delta_S_kB_per_spin':exact_delta_s,
                    'correction_route_delta_S_kB_per_spin':route_s,
                    'S_error_kB_per_spin':route_s-exact_delta_s})
            print(f'B={b}, T={t}: coupling integration completed',flush=True)
    selected=[r for r in results if r['derivative_step_K']==.001]
    assert max(abs(r['F_error_17_nodes_K_per_spin']) for r in selected)<1e-6
    assert max(abs(r['S_error_kB_per_spin']) for r in selected)<1e-4
    # Conservative Pauli L1 measurement scales for the SAME expectation-value
    # error, ignoring state preparation, temperature derivatives and lambda grid.
    dxy=-sw.KD_MEV*sw.GXY**2/2*sw.MEV_TO_K
    dz=sw.KD_MEV*sw.GZ**2*sw.MEV_TO_K
    jpp=(sw.JPP_EXTRA_MEV-3*sw.KD_MEV*sw.GXY**2/4)*sw.MEV_TO_K
    correction_l1=0.
    for theta in [0.,np.pi/3,-np.pi/3]:
        c=np.cos(2*theta);s=np.sin(2*theta)
        correction_l1+=abs(dxy/4+jpp*c/2)+abs(dxy/4-jpp*c/2)+abs(jpp*s)+abs(dz/4)
    full_l1_at_1p5=3*(2*sw.JZZ_MEV*sw.DELTA_XY_OVER_ZZ+sw.JZZ_MEV)*sw.MEV_TO_K/4+sw.GZ*MU_B_OVER_KB*1.5/2
    summary={'scope':['Classical exact diagonalization, no quantum circuit or hardware.',
        'N=9 nearest-neighbor 2D sensitivity proxy; no thermodynamic-limit claim.',
        'Coupling-constant integration and entropy derivative identity benchmark only.',
        'Observation cost comparison omits Gibbs preparation, autocorrelation, derivative and integration costs. Not an end-to-end speedup.'],
        'rows':results,
        'max_F_absolute_error_at_17_nodes_K_per_spin':max(abs(r['F_error_17_nodes_K_per_spin']) for r in selected),
        'max_S_absolute_error_at_dt_1mK_kB_per_spin':max(abs(r['S_error_kB_per_spin']) for r in selected),
        'observable_l1_scales_K_per_spin':{'V':correction_l1,'H0_at_1p5T':full_l1_at_1p5,
           'squared_ratio_for_same_observable_error_only':(full_l1_at_1p5/correction_l1)**2},
        'elapsed_seconds':time.perf_counter()-start}
    save_csv('correction_route_validation.csv',results)
    (OUT/'correction_route_validation.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))

if __name__=='__main__':main()
