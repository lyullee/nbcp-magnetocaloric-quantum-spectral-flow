"""Representative-model pilot for the commensurate-V to fan boundary.

The source paper compares the optimized incommensurate fan with a coplanar
three-sublattice V state.  This pilot reconstructs that comparison directly
from the same M4 Fourier coupling used elsewhere in the repository.  It is a
feasibility gate for a possible 69-model phase-boundary audit; it is not used
in the refrigerator decision unless it reproduces the public representative
transition near 0.86 Hc.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution, minimize

from audit_field_dependent_fan_holdout import (
    model_terms,
    optimize_delta,
    precompute_m3,
)
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
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SETTINGS = EwaldSettings(boundary="tin_foil")
Q_V = np.array([2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0])


def v_energy(variables: np.ndarray, field_T: float,
             j0: np.ndarray, jq: np.ndarray) -> float:
    """Energy per spin of an exact fixed-length coplanar 3-sublattice state."""
    theta = np.asarray(variables[:3], dtype=float)
    phi = float(variables[3])
    direction = np.array([np.cos(phi), np.sin(phi)])
    spins = np.zeros((3, 3), dtype=float)
    spins[:, :2] = S * np.sin(theta)[:, None] * direction[None, :]
    spins[:, 2] = S * np.cos(theta)
    phase = np.exp(-2.0j * np.pi * np.arange(3) / 3.0)
    c0 = spins.mean(axis=0)
    cq = (phase[:, None] * spins).mean(axis=0)
    h = GZ * MU_B_OVER_KB * field_T
    return float(
        0.5 * c0 @ j0 @ c0
        + np.real(np.conjugate(cq) @ jq @ cq)
        - h * c0[2]
    )


def optimize_v(field_T: float, j0: np.ndarray, jq: np.ndarray,
               seed: int, warm: np.ndarray | None = None) -> tuple[float, np.ndarray]:
    bounds = [(-np.pi, np.pi)] * 4
    de = differential_evolution(
        lambda x: v_energy(x, field_T, j0, jq),
        bounds,
        seed=seed,
        popsize=14,
        maxiter=220,
        tol=1.0e-10,
        polish=False,
        workers=1,
        updating="immediate",
    )
    starts = [de.x]
    if warm is not None:
        starts.append(warm)
    local = [
        minimize(
            lambda x: v_energy(x, field_T, j0, jq),
            x0,
            method="Nelder-Mead",
            options={"xatol": 1.0e-11, "fatol": 1.0e-13, "maxiter": 5000},
        )
        for x0 in starts
    ]
    best = min(local, key=lambda row: row.fun)
    return float(best.fun), np.asarray(best.x, dtype=float)


def crossing_from_grid(rows: list[dict]) -> float | None:
    for left, right in zip(rows[:-1], rows[1:]):
        y0 = left["fan_minus_v_microeV"]
        y1 = right["fan_minus_v_microeV"]
        if y0 == 0.0:
            return float(left["field_ratio"])
        if y0 * y1 < 0.0:
            x0, x1 = left["field_ratio"], right["field_ratio"]
            return float(x0 - y0 * (x1 - x0) / (y1 - y0))
    return None


def main() -> None:
    parameters = InterlayerParameters()
    critical = locate_critical_mode_multistart(
        "M4", SETTINGS, (260711, 260717, 260723), parameters
    )
    hc = float(critical["critical_field_T"])
    base = precompute_m3()
    fan_terms = model_terms(base, parameters)
    j0 = coupling_fourier(np.zeros(3), "M4", SETTINGS, parameters)
    jq = coupling_fourier(Q_V, "M4", SETTINGS, parameters)

    rows: list[dict] = []
    warm = None
    ratios = np.linspace(0.80, 0.92, 25)
    for index, ratio in enumerate(ratios):
        field_T = float(ratio * hc)
        fan_e, delta, xi = optimize_delta(fan_terms, field_T)
        v_e, warm = optimize_v(field_T, j0, jq, 261000 + index, warm)
        rows.append({
            "field_ratio": float(ratio),
            "field_T": field_T,
            "fan_energy_K_per_spin": float(fan_e),
            "v_energy_K_per_spin": v_e,
            "fan_minus_v_microeV": float((fan_e - v_e) / 0.011604518121550082),
            "fan_delta_rlu": float(delta),
            "fan_xi": float(xi),
            "v_theta_rad": warm[:3].tolist(),
            "v_plane_phi_rad": float(warm[3]),
        })
        print(f"{ratio:.3f}  dE={rows[-1]['fan_minus_v_microeV']:.6f} microeV", flush=True)

    crossing = crossing_from_grid(rows)
    published_visual = 0.855
    result = {
        "status": "representative V--fan reconstruction pilot",
        "source": "arXiv:2607.23173v1 Figure 10c-d",
        "model_critical_field_T": hc,
        "predicted_v_to_fan_ratio": crossing,
        "published_figure_visual_ratio": published_visual,
        "absolute_ratio_difference": None if crossing is None else abs(crossing - published_visual),
        "rows": rows,
        "pilot_gate": {
            "crossing_found": crossing is not None,
            "within_0p02_Hc_of_public_figure": crossing is not None and abs(crossing - published_visual) <= 0.02,
        },
    }
    result["pilot_pass"] = all(result["pilot_gate"].values())
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    (RESULTS / "vfan_transition_pilot.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.axhline(0.0, color="black", lw=1)
    ax.plot(
        [row["field_ratio"] for row in rows],
        [row["fan_minus_v_microeV"] for row in rows],
        "o-",
        color="#3b5ccc",
    )
    ax.axvline(published_visual, color="#ef6c00", ls="--", label="public figure (~0.855)")
    if crossing is not None:
        ax.axvline(crossing, color="#009688", ls=":", label=f"reconstructed ({crossing:.3f})")
    ax.set_xlabel(r"$H/H_c$")
    ax.set_ylabel(r"$E_{fan}-E_V$ ($\mu$eV/spin)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURES / "vfan_transition_pilot.png", dpi=220)
    plt.close(fig)
    print(json.dumps(result, indent=2))
    if not result["pilot_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
