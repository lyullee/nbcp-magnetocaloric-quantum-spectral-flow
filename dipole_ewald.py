"""Dipolar tensors for the AA-stacked NBCP magnetic sublattice.

Coordinates and cell vectors are in angstrom.  ``ewald_tensor`` returns the
periodic dimensionless geometrical tensor in angstrom**-3 for the energy form

    m_i . D(r_ij) . m_j,

where the free-space tensor is I/r^3 - 3 rr/r^5.  ``spin_coupling_microeV``
then multiplies by mu0 mu_B^2/(4 pi) and the anisotropic g tensor, returning
the coupling matrix between spin operators in microelectronvolt.

The k=0 term of a dipolar Ewald sum is sample-shape dependent.  ``tin_foil``
omits it.  ``sphere_vacuum`` adds the spherical demagnetization term.  Every
reported M3/M4 result must state which convention was used.
"""
from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.special import erfc

MU0_OVER_4PI = 1e-7  # SI
MU_B_J_PER_T = 9.2740100783e-24
J_PER_EV = 1.602176634e-19
ANGSTROM_M = 1e-10
PREF_MICROEV_ANGSTROM3 = (
    MU0_OVER_4PI * MU_B_J_PER_T**2 / ANGSTROM_M**3 / J_PER_EV * 1e6
)


def nbcpo_primitive_cell(a_angstrom: float = 5.312,
                         c_angstrom: float = 6.983) -> np.ndarray:
    """Hexagonal primitive vectors as columns in a Cartesian matrix."""
    return np.array([[a_angstrom, .5*a_angstrom, 0.],
                     [0., np.sqrt(3)/2*a_angstrom, 0.],
                     [0., 0., c_angstrom]], dtype=float)


def reciprocal_cell(cell: np.ndarray) -> np.ndarray:
    """Reciprocal vectors (including 2 pi) as columns."""
    return 2*np.pi*np.linalg.inv(cell).T


def free_space_tensor(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    d = np.linalg.norm(r)
    if d == 0:
        raise ValueError('zero displacement has no free-space pair tensor')
    return np.eye(3)/d**3 - 3*np.outer(r, r)/d**5


def _integer_vectors(basis: np.ndarray, cutoff: float):
    # This conservative bound may enumerate extra points for skew cells; the
    # norm filter below gives the exact spherical cutoff.
    sigma_min = np.linalg.svd(basis, compute_uv=False).min()
    nmax = int(np.ceil(cutoff/sigma_min))+1
    for n in product(range(-nmax, nmax+1), repeat=3):
        yield np.asarray(n, dtype=int), basis@np.asarray(n, dtype=float)


def ewald_tensor(r: np.ndarray, cell: np.ndarray, alpha: float,
                 real_cutoff: float, reciprocal_cutoff: float,
                 boundary: str = 'tin_foil', r_is_self: bool = False) -> np.ndarray:
    """Periodic dipole tensor using a 3D Ewald sum.

    Parameters
    ----------
    r:
        Pair displacement inside the simulation cell, in angstrom.
    cell:
        Direct cell vectors as columns, in angstrom.
    alpha:
        Ewald splitting parameter in inverse angstrom.
    real_cutoff, reciprocal_cutoff:
        Spherical cutoffs in angstrom and inverse angstrom.
    boundary:
        ``tin_foil`` omits the k=0 surface term; ``sphere_vacuum`` adds
        4*pi/(3V) I.  This affects uniform magnetization.
    r_is_self:
        True computes interaction with periodic images of the same site and
        applies the standard self correction.  False requires nonzero r.
    """
    r = np.asarray(r, dtype=float)
    cell = np.asarray(cell, dtype=float)
    if cell.shape != (3, 3) or abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError('cell must be a nonsingular 3x3 column-vector matrix')
    if alpha <= 0 or real_cutoff <= 0 or reciprocal_cutoff <= 0:
        raise ValueError('alpha and cutoffs must be positive')
    if r_is_self:
        if np.linalg.norm(r) > 1e-12:
            raise ValueError('self displacement must be zero')
    elif np.linalg.norm(r) < 1e-12:
        raise ValueError('non-self displacement must be nonzero')

    sqrt_pi = np.sqrt(np.pi)
    real = np.zeros((3, 3), dtype=float)
    for n, shift in _integer_vectors(cell, real_cutoff+np.linalg.norm(r)):
        R = r+shift
        length = np.linalg.norm(R)
        if length < 1e-13 or length > real_cutoff:
            continue
        e = np.exp(-(alpha*length)**2)
        b = erfc(alpha*length)/length**3 + 2*alpha/sqrt_pi*e/length**2
        c = (3*erfc(alpha*length)/length**5 +
             2*alpha/sqrt_pi*e*(3/length**4+2*alpha**2/length**2))
        real += b*np.eye(3)-c*np.outer(R, R)

    recip_basis = reciprocal_cell(cell)
    volume = abs(np.linalg.det(cell))
    recip = np.zeros((3, 3), dtype=float)
    for n, k in _integer_vectors(recip_basis, reciprocal_cutoff):
        k2 = float(k@k)
        if k2 < 1e-24 or np.sqrt(k2) > reciprocal_cutoff:
            continue
        recip += (4*np.pi/volume)*np.exp(-k2/(4*alpha**2))*np.cos(k@r)*np.outer(k,k)/k2

    result = real+recip
    if r_is_self:
        result -= 4*alpha**3/(3*sqrt_pi)*np.eye(3)
    if boundary == 'sphere_vacuum':
        result += 4*np.pi/(3*volume)*np.eye(3)
    elif boundary != 'tin_foil':
        raise ValueError('boundary must be tin_foil or sphere_vacuum')
    return .5*(result+result.T)


def ewald_fourier_tensor(k: np.ndarray, cell: np.ndarray, alpha: float,
                         real_cutoff: float, reciprocal_cutoff: float,
                         boundary: str = 'tin_foil') -> np.ndarray:
    """Fourier transform of the periodic dipole tensor for a Bravais lattice.

    ``k`` is a Cartesian wave vector in inverse angstrom.  The returned tensor
    is the Ewald sum of ``(I/R**3 - 3 RR/R**5) exp(i k.R)`` over every nonzero
    lattice vector ``R``.  Inversion symmetry makes the result real.  At
    ``k=0`` this is identical to ``ewald_tensor(..., r_is_self=True)``.

    The Fourier transform is nonanalytic at the origin.  The documented
    surface term is consequently applied only when ``k`` is reciprocal-lattice
    equivalent to zero.  Critical modes in this project lie near K, so their
    transverse tensor is independent of this convention; the longitudinal
    uniform reference retains the physical sample-shape uncertainty.
    """
    k = np.asarray(k, dtype=float)
    cell = np.asarray(cell, dtype=float)
    if k.shape != (3,):
        raise ValueError('k must be a Cartesian three-vector')
    if cell.shape != (3, 3) or abs(np.linalg.det(cell)) < 1e-12:
        raise ValueError('cell must be a nonsingular 3x3 column-vector matrix')
    if alpha <= 0 or real_cutoff <= 0 or reciprocal_cutoff <= 0:
        raise ValueError('alpha and cutoffs must be positive')
    if boundary not in ('tin_foil', 'sphere_vacuum'):
        raise ValueError('boundary must be tin_foil or sphere_vacuum')

    sqrt_pi = np.sqrt(np.pi)
    real = np.zeros((3, 3), dtype=float)
    for _, R in _integer_vectors(cell, real_cutoff):
        length = np.linalg.norm(R)
        if length < 1e-13 or length > real_cutoff:
            continue
        e = np.exp(-(alpha*length)**2)
        b = erfc(alpha*length)/length**3 + 2*alpha/sqrt_pi*e/length**2
        c = (3*erfc(alpha*length)/length**5 +
             2*alpha/sqrt_pi*e*(3/length**4+2*alpha**2/length**2))
        real += np.cos(k@R)*(b*np.eye(3)-c*np.outer(R, R))

    recip_basis = reciprocal_cell(cell)
    volume = abs(np.linalg.det(cell))
    recip = np.zeros((3, 3), dtype=float)
    # q=k+G is cut off spherically.  The enlarged enumeration bound is needed
    # when k is close to a Brillouin-zone boundary.
    for _, G in _integer_vectors(recip_basis, reciprocal_cutoff+np.linalg.norm(k)):
        q = k+G
        q2 = float(q@q)
        if q2 < 1e-24 or np.sqrt(q2) > reciprocal_cutoff:
            continue
        recip += ((4*np.pi/volume)*np.exp(-q2/(4*alpha**2)) *
                  np.outer(q, q)/q2)

    result = real+recip-4*alpha**3/(3*sqrt_pi)*np.eye(3)
    frac = np.linalg.solve(recip_basis, k)
    at_gamma = np.max(np.abs(frac-np.rint(frac))) < 1e-10
    if boundary == 'sphere_vacuum' and at_gamma:
        result += 4*np.pi/(3*volume)*np.eye(3)
    return .5*(result+result.T)


class EwaldFourierWorkspace:
    """Vectorized reusable workspace for many Fourier Ewald evaluations."""

    def __init__(self, cell: np.ndarray, alpha: float, real_cutoff: float,
                 reciprocal_cutoff: float):
        self.cell = np.asarray(cell, dtype=float)
        self.alpha = float(alpha)
        self.real_cutoff = float(real_cutoff)
        self.reciprocal_cutoff = float(reciprocal_cutoff)
        self.recip_basis = reciprocal_cell(self.cell)
        self.volume = abs(np.linalg.det(self.cell))
        sqrt_pi = np.sqrt(np.pi)

        rv, rt = [], []
        for _, R in _integer_vectors(self.cell, self.real_cutoff):
            length = np.linalg.norm(R)
            if length < 1e-13 or length > self.real_cutoff:
                continue
            e = np.exp(-(self.alpha*length)**2)
            b = (erfc(self.alpha*length)/length**3 +
                 2*self.alpha/sqrt_pi*e/length**2)
            c = (3*erfc(self.alpha*length)/length**5 +
                 2*self.alpha/sqrt_pi*e *
                 (3/length**4+2*self.alpha**2/length**2))
            rv.append(R)
            rt.append(b*np.eye(3)-c*np.outer(R, R))
        self.real_vectors = np.asarray(rv)
        self.real_tensors = np.asarray(rt)

        # Every input is first reduced to the central reciprocal
        # parallelepiped.  This bounds the extra reciprocal vectors needed.
        max_k = .5*sum(np.linalg.norm(self.recip_basis[:, i]) for i in range(3))
        gv = []
        for _, G in _integer_vectors(self.recip_basis,
                                     self.reciprocal_cutoff+max_k):
            gv.append(G)
        self.recip_vectors = np.asarray(gv)
        self.self_term = 4*self.alpha**3/(3*sqrt_pi)*np.eye(3)

    def tensor(self, k: np.ndarray, boundary: str = 'tin_foil') -> np.ndarray:
        if boundary not in ('tin_foil', 'sphere_vacuum'):
            raise ValueError('boundary must be tin_foil or sphere_vacuum')
        frac = np.linalg.solve(self.recip_basis, np.asarray(k, dtype=float))
        reduced_frac = (frac+.5) % 1.0-.5
        reduced_k = self.recip_basis@reduced_frac
        weights = np.cos(self.real_vectors@reduced_k)
        real = np.einsum('n,nij->ij', weights, self.real_tensors)
        q = self.recip_vectors+reduced_k
        q2 = np.einsum('ni,ni->n', q, q)
        mask = (q2 > 1e-24) & (q2 <= self.reciprocal_cutoff**2)
        qm, q2m = q[mask], q2[mask]
        factors = ((4*np.pi/self.volume)*
                   np.exp(-q2m/(4*self.alpha**2))/q2m)
        recip = np.einsum('n,ni,nj->ij', factors, qm, qm)
        result = real+recip-self.self_term
        at_gamma = np.max(np.abs(reduced_frac)) < 1e-10
        if boundary == 'sphere_vacuum' and at_gamma:
            result += 4*np.pi/(3*self.volume)*np.eye(3)
        return .5*(result+result.T)


def spin_coupling_microeV(geometric_tensor: np.ndarray,
                          g=(4.2, 4.2, 4.716)) -> np.ndarray:
    gmat = np.diag(np.asarray(g, dtype=float))
    return PREF_MICROEV_ANGSTROM3*(gmat@geometric_tensor@gmat)


def convergence_demo():
    """Return reproducible unit and Ewald-splitting checks."""
    cell = nbcpo_primitive_cell()
    a1 = cell[:, 0]
    direct = spin_coupling_microeV(free_space_tensor(a1))
    kd = PREF_MICROEV_ANGSTROM3/np.linalg.norm(a1)**3
    # Use a 3x3x3 supercell for a nontrivial periodic pair check.
    supercell = cell@np.diag([3.,3.,3.])
    tensors = {}
    for alpha in [.16, .22, .30]:
        tensors[str(alpha)] = ewald_tensor(a1, supercell, alpha, 42., 8.)
    ref=tensors['0.22']
    max_spread=max(float(np.max(np.abs(t-ref))) for t in tensors.values())
    return {'prefactor_microeV_angstrom3':PREF_MICROEV_ANGSTROM3,
            'nearest_neighbor_Kd_microeV':kd,
            'nearest_neighbor_spin_tensor_microeV':direct.tolist(),
            'expected_NN_diagonal_microeV':[-2*kd*4.2**2,kd*4.2**2,kd*4.716**2],
            'ewald_alpha_values_per_angstrom':[.16,.22,.30],
            'ewald_pair_tensor_max_alpha_spread_per_angstrom3':max_spread,
            'ewald_boundary':'tin_foil'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='results/dipole_ewald_validation.json')
    args=parser.parse_args()
    result=convergence_demo()
    path=Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
