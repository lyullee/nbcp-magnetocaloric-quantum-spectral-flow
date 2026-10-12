"""One-point feasibility test for a 3-layer, 9-spin coplanar V cell."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize, minimize_scalar

from audit_field_dependent_fan_holdout import model_terms, optimize_delta, precompute_m3
from bulk_polarized_stability import (
    EwaldSettings,
    GZ,
    InterlayerParameters,
    S,
    coupling_fourier,
    locate_critical_mode_multistart,
)
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
SETTINGS = EwaldSettings(boundary="tin_foil")


def precompute_cell(parameters: InterlayerParameters) -> dict:
    positions = np.array(
        [(a, b, c) for a in range(3) for b in range(3) for c in range(3)],
        dtype=float,
    )
    kpoints = np.array(
        [(a / 3.0, b / 3.0, c / 3.0)
         for a in range(3) for b in range(3) for c in range(3)],
        dtype=float,
    )
    phases = np.exp(-2.0j * np.pi * (kpoints @ positions.T)) / len(positions)
    couplings = np.array([
        coupling_fourier(k, "M4", SETTINGS, parameters) for k in kpoints
    ])
    # Three in-plane sublattices times three explicit layers.
    site_indices = np.array([
        int(c) * 3 + ((2 * int(a) + int(b)) % 3)
        for a, b, c in positions
    ])
    return {
        "positions": positions,
        "kpoints": kpoints,
        "phases": phases,
        "couplings": couplings,
        "site_indices": site_indices,
    }


def energy(variables: np.ndarray, field_T: float, cell: dict) -> float:
    theta = np.asarray(variables[:9], dtype=float)
    phi = float(variables[9])
    unit = np.array([np.cos(phi), np.sin(phi)])
    unique = np.zeros((9, 3), dtype=float)
    unique[:, :2] = S * np.sin(theta)[:, None] * unit[None, :]
    unique[:, 2] = S * np.cos(theta)
    spins = unique[cell["site_indices"]]
    coeff = cell["phases"] @ spins
    exchange = 0.5 * np.real(np.einsum(
        "ki,kij,kj->", np.conjugate(coeff), cell["couplings"], coeff
    ))
    h = GZ * MU_B_OVER_KB * field_T
    return float(exchange - h * coeff[0, 2].real)


def energy_unrestricted(variables: np.ndarray, field_T: float, cell: dict) -> float:
    """Nine independent fixed-length spins in the 3-layer ABC cell."""
    theta = np.asarray(variables[:9], dtype=float)
    phi = np.asarray(variables[9:18], dtype=float)
    unique = np.column_stack([
        S * np.sin(theta) * np.cos(phi),
        S * np.sin(theta) * np.sin(phi),
        S * np.cos(theta),
    ])
    spins = unique[cell["site_indices"]]
    coeff = cell["phases"] @ spins
    exchange = 0.5 * np.real(np.einsum(
        "ki,kij,kj->", np.conjugate(coeff), cell["couplings"], coeff
    ))
    h = GZ * MU_B_OVER_KB * field_T
    return float(exchange - h * coeff[0, 2].real)


def exact_fan_coefficients(xi: float, max_harmonic: int = 15) -> tuple[float, dict[int, float], float]:
    """Thermodynamic-limit Fourier coefficients from the exact angle in Eq. 6."""
    phi = 2.0 * np.pi * np.arange(8192) / 8192
    cosine = np.cos(phi)
    z = S * (1.0 - 0.5 * xi * xi * cosine * cosine)
    x = S * xi * cosine * np.sqrt(np.maximum(0.0, 1.0 - 0.25 * xi * xi * cosine * cosine))
    z0 = float(np.mean(z))
    z2 = float(np.mean(z * np.cos(2.0 * phi)))
    x_coefficients = {
        n: float(np.mean(x * np.cos(n * phi)))
        for n in range(1, max_harmonic + 1, 2)
    }
    return z0, x_coefficients, z2


def exact_fan_minimum(field_T: float, parameters: InterlayerParameters) -> tuple[float, float, float]:
    """Exact-angle fan energy, retaining transverse harmonics through 15Q."""
    j0 = coupling_fourier(np.zeros(3), "M4", SETTINGS, parameters)
    candidates = []
    for delta in np.linspace(0.010, 0.030, 81):
        q = np.array([2.0 / 3.0 + 2.0 * delta,
                      1.0 / 3.0 + delta, 2.0 / 3.0])
        jq = coupling_fourier(np.mod(q, 1.0), "M4", SETTINGS, parameters)
        _, vectors = np.linalg.eigh(jq[:2, :2])
        direction = vectors[:, 0]
        transverse = {
            n: float(direction @ coupling_fourier(
                np.mod(n * q, 1.0), "M4", SETTINGS, parameters
            )[:2, :2] @ direction)
            for n in range(1, 16, 2)
        }
        j2zz = float(coupling_fourier(
            np.mod(2.0 * q, 1.0), "M4", SETTINGS, parameters
        )[2, 2])
        h = GZ * MU_B_OVER_KB * field_T

        def objective(xi: float) -> float:
            z0, xcoeff, z2 = exact_fan_coefficients(float(xi))
            return float(
                0.5 * j0[2, 2] * z0 * z0
                + sum(transverse[n] * xcoeff[n] * xcoeff[n] for n in xcoeff)
                + j2zz * z2 * z2
                - h * z0
            )

        optimum = minimize_scalar(
            objective, bounds=(0.0, 1.95), method="bounded",
            options={"xatol": 2.0e-10},
        )
        candidates.append((float(optimum.fun), float(delta), float(optimum.x)))
    return min(candidates)


def main() -> None:
    p = InterlayerParameters()
    critical = locate_critical_mode_multistart(
        "M4", SETTINGS, (260711, 260717, 260723), p
    )
    hc = float(critical["critical_field_T"])
    ratio = 0.855
    field_T = ratio * hc
    cell = precompute_cell(p)

    # Validate Fourier normalization against the polarized energy.
    polarized_variables = np.zeros(10)
    numeric_pol = energy(polarized_variables, field_T, cell)
    j0 = coupling_fourier(np.zeros(3), "M4", SETTINGS, p)
    exact_pol = 0.5 * S * S * j0[2, 2] - GZ * MU_B_OVER_KB * field_T * S

    de = differential_evolution(
        lambda x: energy(x, field_T, cell),
        [(-np.pi, np.pi)] * 10,
        seed=261012,
        popsize=22,
        maxiter=650,
        tol=1.0e-10,
        polish=False,
        workers=1,
        updating="immediate",
    )
    local = minimize(
        lambda x: energy(x, field_T, cell),
        de.x,
        method="Nelder-Mead",
        options={"xatol": 1.0e-11, "fatol": 1.0e-13, "maxiter": 15000},
    )
    v_e = float(local.fun)
    unrestricted_start = np.r_[local.x[:9], np.repeat(local.x[9], 9)]
    unrestricted_fits = []
    rng = np.random.default_rng(261013)
    for index in range(8):
        start = unrestricted_start.copy()
        if index:
            start += rng.normal(0.0, 0.08, size=18)
        unrestricted_fits.append(minimize(
            lambda x: energy_unrestricted(x, field_T, cell),
            start,
            method="L-BFGS-B",
            bounds=[(-np.pi, np.pi)] * 18,
            options={"ftol": 1.0e-15, "gtol": 1.0e-10, "maxiter": 5000},
        ))
    unrestricted = min(unrestricted_fits, key=lambda fit: fit.fun)
    unrestricted_v_e = float(unrestricted.fun)
    fan_e, delta, xi = optimize_delta(model_terms(precompute_m3(), p), field_T)
    exact_fan_e, exact_delta, exact_xi = exact_fan_minimum(field_T, p)
    result = {
        "status": "one-point 3-layer ABC-cell feasibility test",
        "field_ratio": ratio,
        "polarized_fourier_normalization_error_K": abs(numeric_pol - exact_pol),
        "fan_energy_K_per_spin": float(fan_e),
        "v_energy_K_per_spin": v_e,
        "unrestricted_abc_v_energy_K_per_spin": unrestricted_v_e,
        "unrestricted_minus_coplanar_v_microeV": float((unrestricted_v_e - v_e) / 0.011604518121550082),
        "fan_minus_v_microeV": float((fan_e - v_e) / 0.011604518121550082),
        "fan_delta_rlu": float(delta),
        "fan_xi": float(xi),
        "exact_angle_fan_energy_K_per_spin": exact_fan_e,
        "exact_angle_fan_minus_v_microeV": float((exact_fan_e - v_e) / 0.011604518121550082),
        "exact_angle_fan_delta_rlu": exact_delta,
        "exact_angle_fan_xi": exact_xi,
        "optimizer_success": bool(local.success),
        "optimized_theta_rad": local.x[:9].tolist(),
        "optimized_phi_rad": float(local.x[9]),
        "unrestricted_optimizer_success": bool(unrestricted.success),
        "unrestricted_theta_rad": unrestricted.x[:9].tolist(),
        "unrestricted_phi_rad": unrestricted.x[9:18].tolist(),
        "interpretation": (
            "The public representative transition is visually near 0.855 Hc, "
            "so a faithful V reconstruction should give a near-zero energy difference."
        ),
    }
    (ROOT / "results" / "vfan_abc_cell_feasibility.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
