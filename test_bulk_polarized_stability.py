import numpy as np

from bulk_polarized_stability import (EwaldSettings, analytic_m0_field,
                                       coupling_fourier, critical_energy_at_k,
                                       locate_critical_mode)
from screen_weak_terms import (DELTA_XY_OVER_ZZ, JPP_EXTRA_MEV, JZZ_MEV,
                               MEV_TO_K)


def run():
    settings=EwaldSettings()
    m0=locate_critical_mode('M0',settings)
    assert abs(m0['critical_field_T']-analytic_m0_field())<1e-7
    # A triangular-lattice K point is symmetry-equivalent to the numerical
    # optimum and has u,v=(1/3,2/3) in this primitive reciprocal basis.
    k=np.asarray(m0['critical_wavevector_rlu'])
    candidates=np.array([[1/3,2/3],[2/3,1/3],[1/3,1/3],[2/3,2/3]])
    periodic_delta=(k[:2]-candidates+.5)%1-.5
    assert np.min(np.linalg.norm(periodic_delta,axis=1))<1e-5
    # Fourier matrices are real symmetric and reciprocal-periodic.
    probe=np.array([.213,.371,.117])
    a=coupling_fourier(probe,'M3',settings)
    b=coupling_fourier(probe+np.array([1.,-1.,2.]),'M3',settings)
    assert np.max(np.abs(a-a.T))<1e-12
    assert np.max(np.abs(a-b))<1e-9
    # Eq. (5) in the source gives an independent analytic M1 dispersion.
    u,v=.271,.593
    cosines=np.cos(2*np.pi*np.array([u,v,u-v]))
    f=cosines[0]+np.exp(-2j*np.pi/3)*cosines[1]+np.exp(2j*np.pi/3)*cosines[2]
    expected=.5*2*(3*JZZ_MEV-JZZ_MEV*DELTA_XY_OVER_ZZ*cosines.sum()+
                    2*abs(JPP_EXTRA_MEV*f))*MEV_TO_K
    actual=critical_energy_at_k(np.array([u,v,0.]),'M1',settings)
    assert abs(actual-expected)<1e-11
    print('bulk polarized stability validation passed')


if __name__=='__main__':run()
