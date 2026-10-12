"""Pilot reconstruction of the public field-dependent fan incommensuration.

This is a feasibility audit for a possible major manuscript extension.  It
implements the fixed-length fan ansatz of arXiv:2607.23173v1, Eqs. (6)--(7),
using the same M4 Fourier coupling already used for the bulk instability.
The source paper reports that its representative model captures the trend but
underestimates the experimental rate of change of delta(H).  Reproducing that
behaviour is a prerequisite before the whole retained M4 set can be tested.

The experimental points below are digitized from the source-vector Figure 7d
(fields 1.35--1.65 T, Hc=1.65 T).  They are used only for this pilot and carry
the approximately 0.002 r.l.u. uncertainty quoted for the near-Hc point.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

from bulk_polarized_stability import (
    EwaldSettings,
    GZ,
    InterlayerParameters,
    S,
    coupling_fourier,
    locate_critical_mode_multistart,
)
from screen_interlayer_parameter_ensemble import canonical_soft_mode
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SETTINGS = EwaldSettings(boundary="tin_foil")

# Blue data points in the source-vector Figure 7d.  The first two points are
# commensurate within the plotting resolution; the fan phase starts above the
# V--fan transition, near H/Hc=0.86.
EXPERIMENTAL_FIELD_T = np.array([1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65])
EXPERIMENTAL_DELTA_RLU = np.array([0.0, 0.0, 0.0090, 0.0138, 0.0190, 0.0226, 0.0270])
EXPERIMENTAL_HC_T = 1.65


def fan_fourier_energy_terms(delta: float, parameters: InterlayerParameters) -> dict:
    """Return the Fourier couplings needed by the paper's fan ansatz."""
    # The code's reciprocal basis represents the paper's symmetry-equivalent
    # eta branch as K=(2/3,1/3), displaced along (2,1), with qz near 2/3.
    # This is the branch selected by the independently optimized M4 soft mode.
    q = np.array([2.0 / 3.0 + 2.0 * delta,
                  1.0 / 3.0 + delta, 2.0 / 3.0])
    j0 = coupling_fourier(np.zeros(3), "M4", SETTINGS, parameters)
    jq = coupling_fourier(q, "M4", SETTINGS, parameters)
    j2q = coupling_fourier(np.mod(2.0 * q, 1.0), "M4", SETTINGS, parameters)
    j3q = coupling_fourier(np.mod(3.0 * q, 1.0), "M4", SETTINGS, parameters)
    _, eigenvectors = np.linalg.eigh(jq[:2, :2])
    transverse = eigenvectors[:, 0]
    return {
        "delta": float(delta),
        "j0_zz": float(j0[2, 2]),
        "jq_xx": float(transverse @ jq[:2, :2] @ transverse),
        "j2q_zz": float(j2q[2, 2]),
        "j3q_xx": float(transverse @ j3q[:2, :2] @ transverse),
    }


def fan_energy(xi: float, field_T: float, row: dict) -> float:
    """Classical energy per spin in kelvin for Eqs. (6)--(7).

    The amplitudes are the Fourier coefficients of the source paper's cubic
    fan expansion.  At xi=0 the expression reduces to the polarized state and
    its quadratic instability exactly matches the transverse-Hessian formula.
    """
    h = GZ * MU_B_OVER_KB * field_T
    z0 = S * (1.0 - xi * xi / 4.0)
    z2_cos = -S * xi * xi / 4.0
    x1_cos = S * (xi - 3.0 * xi**3 / 32.0)
    x3_cos = -S * xi**3 / 32.0
    return float(
        0.5 * row["j0_zz"] * z0 * z0
        + 0.25 * row["j2q_zz"] * z2_cos * z2_cos
        + 0.25 * row["jq_xx"] * x1_cos * x1_cos
        + 0.25 * row["j3q_xx"] * x3_cos * x3_cos
        - h * z0
    )


def optimize_curve(parameters: InterlayerParameters, model_hc_T: float) -> tuple[list[dict], list[dict]]:
    delta_grid = np.linspace(0.0, 0.040, 321)
    terms = [fan_fourier_energy_terms(delta, parameters) for delta in delta_grid]
    field_ratios = np.r_[np.linspace(0.84, 0.98, 15), 0.99, 0.995, 0.999]
    curve = []
    for ratio in field_ratios:
        field_T = float(ratio * model_hc_T)
        candidates = []
        for row in terms:
            optimum = minimize_scalar(
                lambda xi: fan_energy(float(xi), field_T, row),
                bounds=(0.0, 1.95), method="bounded",
                options={"xatol": 2.0e-10},
            )
            candidates.append((float(optimum.fun), float(row["delta"]), float(optimum.x)))
        energy, delta, xi = min(candidates)
        curve.append({
            "field_ratio": float(ratio),
            "field_T": field_T,
            "delta_opt_rlu": delta,
            "xi_opt": xi,
            "energy_K_per_spin": energy,
        })
    return curve, terms


def main() -> None:
    parameters = InterlayerParameters()
    mode = locate_critical_mode_multistart(
        "M4", SETTINGS, seeds=(260711, 271, 27182), interlayer=parameters
    )
    model_hc = float(mode["critical_field_T"])
    critical = canonical_soft_mode(mode["critical_wavevector_rlu"])
    curve, _ = optimize_curve(parameters, model_hc)

    ratios = np.array([row["field_ratio"] for row in curve])
    deltas = np.array([row["delta_opt_rlu"] for row in curve])
    experimental_ratio = EXPERIMENTAL_FIELD_T / EXPERIMENTAL_HC_T
    stable = experimental_ratio >= 1.45 / 1.65
    predicted = np.interp(experimental_ratio[stable], ratios, deltas)
    residual = predicted - EXPERIMENTAL_DELTA_RLU[stable]

    near_hc_delta = float(deltas[-1])
    critical_delta = float(critical["delta_inplane_rlu"])
    trend_monotone = bool(np.all(np.diff(deltas[ratios >= 0.86]) >= -1.0e-12))
    result = {
        "status": "pilot only; not yet propagated over the retained M4 set",
        "source": "arXiv:2607.23173v1 Eqs. (6)-(7) and source-vector Figure 7d",
        "representative_parameters_microeV": parameters.__dict__,
        "model_critical_field_T": model_hc,
        "critical_mode_delta_rlu": critical_delta,
        "fan_delta_at_0p999_Hc_rlu": near_hc_delta,
        "near_hc_internal_consistency_abs_rlu": abs(near_hc_delta - critical_delta),
        "monotone_delta_above_0p86_Hc": trend_monotone,
        "experimental_holdout": {
            "field_T": EXPERIMENTAL_FIELD_T.tolist(),
            "field_ratio": experimental_ratio.tolist(),
            "digitized_delta_rlu": EXPERIMENTAL_DELTA_RLU.tolist(),
            "comparison_mask_field_at_least_1p45T": stable.tolist(),
            "predicted_delta_rlu": predicted.tolist(),
            "residual_rlu": residual.tolist(),
            "mae_rlu": float(np.mean(np.abs(residual))),
            "source_near_Hc_quoted_uncertainty_rlu": 0.002,
        },
        "curve": curve,
        "pilot_gates": {
            "near_hc_matches_soft_mode_within_0p002_rlu":
                abs(near_hc_delta - critical_delta) <= 0.002,
            "delta_increases_monotonically_in_fan_window": trend_monotone,
            "reproduces_source_reported_underresponse":
                # The source says the experimental delta falls more rapidly
                # below Hc, so the mean-field curve remains above the data.
                bool(np.mean(residual[:-1]) > 0.0),
        },
        "interpretation": (
            "Passing the internal and qualitative gates only establishes that "
            "the published fan calculation can be reconstructed.  A major "
            "extension requires exact evaluation of all retained M4 models and "
            "joint constraints from delta(H), qz(H), transition fields, and the "
            "out-of-plane bandwidth."
        ),
    }
    result["pilot_pass"] = bool(all(result["pilot_gates"].values()))

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    (RESULTS / "fan_curve_holdout_pilot.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(6.7, 4.6))
    ax.plot(ratios, deltas, "o-", color="#54a24b", label="reconstructed representative M4")
    ax.plot(experimental_ratio, EXPERIMENTAL_DELTA_RLU, "o-", color="#4c5bd4",
            linewidth=1.5, markersize=7, label="public diffraction (digitized)")
    ax.axhline(critical_delta, color="black", linestyle="--", linewidth=1,
               label="independent soft-mode limit")
    ax.set(xlabel=r"reduced field $H/H_c$", ylabel=r"incommensuration $\delta$ (r.l.u.)",
           xlim=(0.81, 1.005), ylim=(-0.001, 0.031),
           title="Fan-curve reconstruction feasibility audit")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fan_curve_holdout_pilot.png", dpi=220)
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
