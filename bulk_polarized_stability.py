"""Thermodynamic-limit polarized-state stability for NBCP models M0--M4.

This is a zero-temperature linear-spin-wave / transverse-Hessian diagnostic.
It locates the bulk gap-closing field and wave vector.  It is not a finite-T
entropy calculation.  Its purpose is to test whether finite-cluster shifts in
the magnetocaloric window have the sign and scale expected in the bulk before
spending scarce quantum-hardware calls.

For H = 1/2 sum_(i,R) S_i^T J(R) S_(i+R) - h sum_i S_i^z, the c-polarized
state loses stability when

    h_c(k) = S [J_zz(0) - lambda_min(J_perp(k))].

For M0 this reproduces Eq. (5) of arXiv:2607.23173 exactly.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

from dipole_ewald import (EwaldFourierWorkspace, nbcpo_primitive_cell,
                          reciprocal_cell, spin_coupling_microeV)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import (DELTA_XY_OVER_ZZ, GZ, JPP_EXTRA_MEV,
                               JZZ_MEV, KD_MEV, GXY, MEV_TO_K)

MICROEV_TO_K = MEV_TO_K / 1000.0
S = 0.5
J2_MICROEV, JA_MICROEV, JB_MICROEV, DELTA_C = 2.8, 0.0, 7.9, 0.08


@dataclass(frozen=True)
class EwaldSettings:
    alpha: float = 0.22
    real_cutoff: float = 42.0
    reciprocal_cutoff: float = 8.0
    boundary: str = 'tin_foil'


def _bond_exchange_tensor(theta: float, include_jpp: bool) -> np.ndarray:
    """Intralayer exchange tensor in kelvin for one positive NN bond."""
    jxy = JZZ_MEV * DELTA_XY_OVER_ZZ * MEV_TO_K
    jzz = JZZ_MEV * MEV_TO_K
    jpp = JPP_EXTRA_MEV * MEV_TO_K if include_jpp else 0.0
    c, s = np.cos(2*theta), np.sin(2*theta)
    return np.array([[jxy+2*jpp*c, -2*jpp*s, 0.0],
                     [-2*jpp*s, jxy-2*jpp*c, 0.0],
                     [0.0, 0.0, jzz]])


def _m2_bond_exchange_tensor(theta: float) -> np.ndarray:
    """M1 plus the paper's nearest-neighbor dipolar tensor (M2)."""
    jxy = JZZ_MEV*DELTA_XY_OVER_ZZ-KD_MEV*GXY**2/2
    jzz = JZZ_MEV+KD_MEV*GZ**2
    jpp = JPP_EXTRA_MEV-3*KD_MEV*GXY**2/4
    c, s = np.cos(2*theta), np.sin(2*theta)
    return MEV_TO_K*np.array([[jxy+2*jpp*c, -2*jpp*s, 0.0],
                              [-2*jpp*s, jxy-2*jpp*c, 0.0],
                              [0.0, 0.0, jzz]])


def _exchange_fourier(frac_k: np.ndarray, model: str) -> np.ndarray:
    """Fourier exchange matrix including positive and negative bonds."""
    u, v, w = np.asarray(frac_k, dtype=float)
    include_jpp = model != 'M0'
    # Positive intralayer primitive displacements and the bond angles used in
    # the public Hamiltonian implementation.
    terms = [((1, 0, 0), 0.0),
             ((0, 1, 0), np.pi/3),
             ((1, -1, 0), -np.pi/3)]
    result = np.zeros((3, 3), dtype=float)
    for d, theta in terms:
        phase = 2*np.pi*(d[0]*u+d[1]*v+d[2]*w)
        tensor = (_m2_bond_exchange_tensor(theta) if model == 'M2' else
                  _bond_exchange_tensor(theta, include_jpp))
        result += 2*np.cos(phase)*tensor

    if model == 'M4':
        inter = [((0, 0, 1), J2_MICROEV)]
        inter += [(d, JA_MICROEV) for d in
                  ((1, 0, 1), (0, 1, 1), (1, -1, 1))]
        inter += [(d, JB_MICROEV) for d in
                  ((-1, 0, 1), (0, -1, 1), (-1, 1, 1))]
        for d, j_microev in inter:
            phase = 2*np.pi*(d[0]*u+d[1]*v+d[2]*w)
            tensor = np.diag([DELTA_C*j_microev,
                              DELTA_C*j_microev,
                              j_microev])*MICROEV_TO_K
            result += 2*np.cos(phase)*tensor
    return result


@lru_cache(maxsize=16)
def _ewald_workspace(alpha: float, real_cutoff: float,
                     reciprocal_cutoff: float) -> EwaldFourierWorkspace:
    return EwaldFourierWorkspace(nbcpo_primitive_cell(), alpha, real_cutoff,
                                 reciprocal_cutoff)


def _dipole_fourier(frac_k: np.ndarray, settings: EwaldSettings) -> np.ndarray:
    cell = nbcpo_primitive_cell()
    k = reciprocal_cell(cell) @ np.asarray(frac_k, dtype=float)
    workspace = _ewald_workspace(settings.alpha, settings.real_cutoff,
                                 settings.reciprocal_cutoff)
    geometric = workspace.tensor(k, settings.boundary)
    return spin_coupling_microeV(geometric)*MICROEV_TO_K


def coupling_fourier(frac_k: np.ndarray, model: str,
                     settings: EwaldSettings) -> np.ndarray:
    if model not in ('M0', 'M1', 'M2', 'M3', 'M4'):
        raise ValueError('bulk stability supports M0 through M4')
    result = _exchange_fourier(frac_k, model)
    if model in ('M3', 'M4'):
        result += _dipole_fourier(frac_k, settings)
    return .5*(result+result.T)


def critical_energy_at_k(frac_k: np.ndarray, model: str,
                         settings: EwaldSettings) -> float:
    j0 = coupling_fourier(np.zeros(3), model, settings)
    jk = coupling_fourier(frac_k, model, settings)
    return float(S*(j0[2, 2]-np.linalg.eigvalsh(jk[:2, :2])[0]))


def locate_critical_mode(model: str, settings: EwaldSettings,
                         seed: int = 260711) -> dict:
    dimension = 3 if model in ('M3', 'M4') else 2
    j0 = coupling_fourier(np.zeros(3), model, settings)

    def objective(x):
        frac = np.r_[x, 0.0] if dimension == 2 else np.asarray(x)
        jk = coupling_fourier(frac, model, settings)
        return -float(S*(j0[2, 2]-np.linalg.eigvalsh(jk[:2, :2])[0]))

    # Deterministic global search followed by a tighter local polishing step.
    de = differential_evolution(objective, [(0., 1.)]*dimension, seed=seed,
                                popsize=12, maxiter=80, tol=2e-8,
                                polish=False, workers=1, updating='immediate')
    local = minimize(objective, de.x, method='Nelder-Mead',
                     options={'xatol': 2e-10, 'fatol': 2e-11,
                              'maxiter': 1200})
    x = np.mod(local.x, 1.0)
    frac = np.r_[x, 0.0] if dimension == 2 else x
    h_k = critical_energy_at_k(frac, model, settings)
    bc = h_k/(GZ*MU_B_OVER_KB)
    j0 = coupling_fourier(np.zeros(3), model, settings)
    jk = coupling_fourier(frac, model, settings)
    return {
        'model': model,
        'boundary': settings.boundary,
        'critical_wavevector_rlu': frac.tolist(),
        'critical_energy_K': h_k,
        'critical_field_T': bc,
        'Jzz_Gamma_K': float(j0[2, 2]),
        'transverse_min_eigenvalue_K': float(np.linalg.eigvalsh(jk[:2, :2])[0]),
        'optimizer_success': bool(local.success),
        'optimizer_message': str(local.message),
        'global_objective': float(-de.fun),
        'local_objective': h_k,
        'ewald': settings.__dict__ if model in ('M3', 'M4') else None,
    }


def analytic_m0_field() -> float:
    jzz = JZZ_MEV*MEV_TO_K
    h = 2*S*jzz*(3+1.5*DELTA_XY_OVER_ZZ)
    return h/(GZ*MU_B_OVER_KB)


def run(output: Path) -> dict:
    rows=[]
    for boundary in ('tin_foil', 'sphere_vacuum'):
        settings=EwaldSettings(boundary=boundary)
        for model in ('M0', 'M1', 'M2', 'M3', 'M4'):
            # Boundary has no effect without long-range dipoles; avoid duplicate
            # M0/M1 rows while retaining both physical conventions for M3/M4.
            if boundary == 'sphere_vacuum' and model in ('M0', 'M1', 'M2'):
                continue
            row=locate_critical_mode(model,settings)
            print(json.dumps(row),flush=True)
            rows.append(row)
    m0=next(r for r in rows if r['model']=='M0')
    analytic=analytic_m0_field()
    result={
        'method':'thermodynamic-limit polarized-state transverse Hessian (LSWT)',
        'scope':'zero-temperature bulk critical-mode diagnostic; not finite-temperature entropy',
        'source':'arXiv:2607.23173v1 Eq. 5 and Table 1',
        'analytic_M0_field_T':analytic,
        'numeric_M0_field_T':m0['critical_field_T'],
        'M0_absolute_validation_error_T':abs(m0['critical_field_T']-analytic),
        'reported_full_model_mean_field_Hc_T':1.86,
        'rows':rows,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,
                        default=Path('results/bulk_polarized_stability.json'))
    args=parser.parse_args()
    result=run(args.output)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
