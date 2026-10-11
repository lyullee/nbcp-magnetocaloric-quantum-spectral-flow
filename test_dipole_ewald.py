import numpy as np

from dipole_ewald import (EwaldFourierWorkspace, PREF_MICROEV_ANGSTROM3, convergence_demo,
                          ewald_fourier_tensor, ewald_tensor, free_space_tensor,
                          nbcpo_primitive_cell, reciprocal_cell,
                          spin_coupling_microeV)


def run():
    cell=nbcpo_primitive_cell()
    assert abs(np.linalg.det(cell)-5.312**2*np.sqrt(3)/2*6.983)<1e-10
    assert np.max(np.abs(cell.T@reciprocal_cell(cell)-2*np.pi*np.eye(3)))<1e-12
    a1=cell[:,0]
    d=free_space_tensor(a1)
    assert abs(np.trace(d))<1e-14
    kd=PREF_MICROEV_ANGSTROM3/5.312**3
    expected=np.diag([-2*kd*4.2**2,kd*4.2**2,kd*4.716**2])
    assert np.max(np.abs(spin_coupling_microeV(d)-expected))<1e-12
    # The paper reports Kd,nn approximately 0.36 microeV before g factors.
    assert abs(kd-.36)<.005
    check=convergence_demo()
    assert check['ewald_pair_tensor_max_alpha_spread_per_angstrom3']<2e-8
    # The two documented boundary conventions differ only by the analytic
    # spherical k=0 term.
    sc=cell@np.diag([3.,3.,3.])
    tin=ewald_tensor(a1,sc,.22,42.,8.,'tin_foil')
    sphere=ewald_tensor(a1,sc,.22,42.,8.,'sphere_vacuum')
    expected_surface=4*np.pi/(3*abs(np.linalg.det(sc)))*np.eye(3)
    assert np.max(np.abs((sphere-tin)-expected_surface))<1e-12
    # The Bravais-lattice Fourier sum at Gamma must reproduce the self-image
    # Ewald sum, and it must be periodic under reciprocal lattice vectors.
    primitive_self=ewald_tensor(np.zeros(3),cell,.22,42.,8.,'tin_foil',True)
    fourier_gamma=ewald_fourier_tensor(np.zeros(3),cell,.22,42.,8.)
    assert np.max(np.abs(primitive_self-fourier_gamma))<1e-12
    reciprocal=reciprocal_cell(cell)
    k=.213*reciprocal[:,0]+.371*reciprocal[:,1]+.117*reciprocal[:,2]
    d1=ewald_fourier_tensor(k,cell,.22,42.,8.)
    d2=ewald_fourier_tensor(k+reciprocal[:,0]-reciprocal[:,2],cell,.22,42.,8.)
    assert np.max(np.abs(d1-d2))<1e-10
    workspace=EwaldFourierWorkspace(cell,.22,42.,8.)
    assert np.max(np.abs(workspace.tensor(k)-d1))<1e-10
    assert np.max(np.abs(workspace.tensor(k+reciprocal[:,0])-d1))<1e-10
    print('dipole Ewald validation passed')


if __name__=='__main__':run()
