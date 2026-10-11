"""Set-valued interlayer screen constrained by public propagation vectors.

The 2026 diffraction paper labels its M4 interlayer parameter set as
representative and non-unique.  This script therefore evaluates a deterministic
3x5x5x5 local stress grid around that set, retains parameter combinations whose
polarized-state soft mode is compatible with the public near-saturation
propagation vector, and propagates every retained critical-field ablation to
the existing engineering decision metrics.

The grid is a designed robustness envelope, not a posterior distribution.
Fractions of retained grid points must not be interpreted as probabilities.
"""
from __future__ import annotations

import csv
import itertools
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from bulk_polarized_stability import (
    EwaldSettings,
    InterlayerParameters,
    MICROEV_TO_K,
    S,
    coupling_fourier,
    locate_critical_mode,
)
from bulk_to_engineering_surrogate import inverse_field, transform_field
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SETTINGS = EwaldSettings(boundary="tin_foil")

# Deterministic local stress envelope.  Every exchange remains below 10% of
# Jzz=122.5 microeV and Delta_c remains strongly Ising-like.  These are test
# bounds around the published representative set, not fitted confidence bounds.
GRID = {
    "j2_microev": (1.4, 2.8, 4.2),
    "ja_microev": (-2.8, -1.4, 0.0, 1.4, 2.8),
    "jb_microev": (3.95, 5.925, 7.9, 9.875, 11.85),
    "delta_c": (0.04, 0.06, 0.08, 0.10, 0.12),
}

# Public 65 mK near-saturation diffraction: delta_c=0.027(2), qz=0.328(2).
# The source model itself replaces qz by the nearest commensurate 1/3.  The
# declared deterministic bands span the quoted experimental uncertainty and
# that documented approximation.
DELTA_BAND = (0.021, 0.033)  # 3 quoted standard errors around 0.027
QZ_BAND = (0.322, 1.0 / 3.0 + 0.004)
KM_LINE_RESIDUAL_MAX = 0.003
EXPERIMENTAL_HC = (1.65, 1.70)
METHODS = ("reduced_field", "critical_offset")
TEMPERATURES = (0.09, 0.10, 0.12, 0.15)


def _periodic_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5


def canonical_soft_mode(frac_q) -> dict:
    """Return K--M incommensuration and folded qz in the code basis."""
    q = np.mod(np.asarray(frac_q, dtype=float), 1.0)
    k_points = (np.array([1 / 3, 2 / 3]), np.array([2 / 3, 1 / 3]))
    # A displacement t*d has experimental incommensuration |t|.  For example,
    # the published representative solution is d=(2,1), t about 0.025.
    directions = (np.array([2.0, 1.0]), np.array([1.0, 2.0]),
                  np.array([1.0, -1.0]))
    candidates = []
    for k in k_points:
        displacement = _periodic_delta(q[:2], k)
        for direction in directions:
            t = float(displacement @ direction / (direction @ direction))
            residual = float(np.linalg.norm(displacement - t * direction))
            candidates.append((residual, abs(t), t, k, direction))
    residual, delta, signed, k, direction = min(candidates, key=lambda x: x[0])
    return {
        "delta_inplane_rlu": delta,
        "signed_delta_inplane_rlu": signed,
        "km_line_residual_rlu": residual,
        "nearest_K_rlu": k.tolist(),
        "km_direction": direction.tolist(),
        "qz_folded_rlu": float(min(q[2], 1.0 - q[2])),
    }


def out_of_plane_bandwidth_microev(parameters: InterlayerParameters) -> float:
    """One-magnon bandwidth along (00l), matching the source convention."""
    ws = np.linspace(0.0, 1.0, 161)
    transverse = []
    for w in ws:
        matrix = coupling_fourier(np.array([0.0, 0.0, w]), "M4", SETTINGS,
                                  parameters)
        transverse.append(np.linalg.eigvalsh(matrix[:2, :2])[0])
    # The S=1/2 prefactor turns the representative result into 4.38 microeV,
    # reproducing the paper's reported 4.4 microeV value.
    return float(S * np.ptp(transverse) / MICROEV_TO_K)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def engineering_rows(retained: list[dict], m3_hc: float) -> list[dict]:
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)
    profile = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",", names=True, encoding="utf-8-sig")
    profile = profile[profile["anchor_N"] == 12]
    rows = []
    for temperature in TEMPERATURES:
        p = np.sort(profile[np.isclose(profile["temperature_K"], temperature)],
                    order="field_T")
        fields = p["field_T"]
        heat = p["Q_proxy_mJ_per_cm3"]
        source_bstar = float(fields[int(np.argmax(heat))])
        qmax = float(np.max(heat))
        for sample in retained:
            shift = sample["critical_field_shift_vs_M3_T"]
            for anchor in EXPERIMENTAL_HC:
                p2_hc = anchor
                p0_hc = anchor - shift
                for method in METHODS:
                    p2_b = float(transform_field(source_bstar, source_hc,
                                                 p2_hc, method))
                    p0_b = float(transform_field(source_bstar, source_hc,
                                                 p0_hc, method))
                    source_eval = float(inverse_field(p0_b, source_hc,
                                                      p2_hc, method))
                    q_at_p0 = float(np.interp(source_eval, fields, heat))
                    rows.append({
                        "sample_id": sample["sample_id"],
                        "temperature_K": temperature,
                        "experimental_Hc_anchor_T": anchor,
                        "mapping_method": method,
                        "critical_field_shift_vs_M3_T": shift,
                        "P0_selected_field_T": p0_b,
                        "P2_optimal_field_T": p2_b,
                        "delta_B_star_T": p2_b - p0_b,
                        "Q_at_P0_field_mJ_per_cm3": q_at_p0,
                        "Qmax_mJ_per_cm3_inherited": qmax,
                        "regret_percent": 100.0 * (qmax - q_at_p0) / qmax,
                    })
    return rows


def range_of(rows: list[dict], key: str) -> list[float]:
    values = [float(row[key]) for row in rows]
    return [min(values), max(values)]


def make_figure(all_rows: list[dict], engineering: list[dict]) -> None:
    FIGURES.mkdir(exist_ok=True)
    accepted = [row for row in all_rows if row["propagation_vector_pass"]]
    rejected = [row for row in all_rows if not row["propagation_vector_pass"]]
    e90 = [row for row in engineering if row["temperature_K"] == 0.09]

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    ax = axes[0]
    if rejected:
        ax.scatter([r["delta_inplane_rlu"] for r in rejected],
                   [r["qz_folded_rlu"] for r in rejected], s=20,
                   color="#b8b8b8", alpha=0.65, label="outside public-q band")
    scatter = ax.scatter([r["delta_inplane_rlu"] for r in accepted],
                         [r["qz_folded_rlu"] for r in accepted],
                         c=[r["critical_field_shift_vs_M3_T"] for r in accepted],
                         s=36, cmap="viridis", edgecolor="black", linewidth=0.25,
                         label="retained")
    ax.axvspan(*DELTA_BAND, color="#4c78a8", alpha=0.10)
    ax.axhspan(*QZ_BAND, color="#f58518", alpha=0.10)
    ax.errorbar([0.027], [0.328], xerr=[0.002], yerr=[0.002], fmt="o",
                color="black", label="public 65 mK")
    ax.set(xlabel=r"in-plane $\delta_c$ (r.l.u.)",
           ylabel=r"folded $q_z$ (r.l.u.)", title="Public-q constrained set")
    ax.legend(fontsize=8, loc="best")
    fig.colorbar(scatter, ax=ax, label=r"$\Delta H_c$ vs M3 (T)")

    ax = axes[1]
    ax.scatter([r["delta_B_star_T"] for r in e90],
               [r["regret_percent"] for r in e90], s=20, alpha=0.60,
               c=[r["critical_field_shift_vs_M3_T"] for r in e90],
               cmap="viridis")
    ax.axvline(0.05, color="#e45756", linestyle="--", label="0.05 T gate")
    ax.axhline(5.0, color="#e45756", linestyle=":", label="5% gate")
    ax.set(xlabel=r"90 mK $\Delta B^*$ (T)", ylabel="regret (%)",
           title="Engineering effect of retained models")
    ax.legend(fontsize=8)

    ax = axes[2]
    accepted_shift = [r["critical_field_shift_vs_M3_T"] for r in accepted]
    ax.hist(accepted_shift, bins=min(12, max(4, len(set(accepted_shift)))),
            color="#54a24b", edgecolor="white")
    ax.axvline(0.0952738, color="black", linestyle="--",
               label="published representative")
    ax.set(xlabel=r"$H_c$(M4)$-H_c$(M3) (T)", ylabel="grid count",
           title="Set-valued critical-field ablation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "interlayer_parameter_ensemble.png", dpi=220)
    plt.close(fig)


def main() -> None:
    m3 = locate_critical_mode("M3", SETTINGS)
    m3_hc = float(m3["critical_field_T"])
    keys = tuple(GRID)
    combinations = list(itertools.product(*(GRID[k] for k in keys)))
    rows = []
    for index, values in enumerate(combinations, start=1):
        p = InterlayerParameters(**dict(zip(keys, values)))
        result = locate_critical_mode("M4", SETTINGS, seed=260711 + index,
                                      interlayer=p)
        mode = canonical_soft_mode(result["critical_wavevector_rlu"])
        q_pass = (
            DELTA_BAND[0] <= mode["delta_inplane_rlu"] <= DELTA_BAND[1]
            and QZ_BAND[0] <= mode["qz_folded_rlu"] <= QZ_BAND[1]
            and mode["km_line_residual_rlu"] <= KM_LINE_RESIDUAL_MAX
        )
        row = {
            "sample_id": index,
            **asdict(p),
            "critical_field_T": float(result["critical_field_T"]),
            "critical_field_shift_vs_M3_T":
                float(result["critical_field_T"] - m3_hc),
            "q_u_rlu": float(result["critical_wavevector_rlu"][0]),
            "q_v_rlu": float(result["critical_wavevector_rlu"][1]),
            "q_w_rlu": float(result["critical_wavevector_rlu"][2]),
            "delta_inplane_rlu": mode["delta_inplane_rlu"],
            "qz_folded_rlu": mode["qz_folded_rlu"],
            "km_line_residual_rlu": mode["km_line_residual_rlu"],
            "out_of_plane_bandwidth_microeV":
                out_of_plane_bandwidth_microev(p),
            "propagation_vector_pass": bool(q_pass),
        }
        rows.append(row)
        if index % 25 == 0:
            print(f"completed {index}/{len(combinations)}", flush=True)

    retained = [row for row in rows if row["propagation_vector_pass"]]
    if not retained:
        raise RuntimeError("no stress-grid model passed the public-q constraints")
    engineering = engineering_rows(retained, m3_hc)
    e90 = [row for row in engineering if row["temperature_K"] == 0.09]

    params = list(GRID)
    sensitivity = {}
    for target in ("critical_field_shift_vs_M3_T", "qz_folded_rlu",
                   "delta_inplane_rlu", "out_of_plane_bandwidth_microeV"):
        sensitivity[target] = {}
        for param in params:
            rho, pvalue = spearmanr([r[param] for r in rows],
                                    [r[target] for r in rows])
            sensitivity[target][param] = {
                "spearman_rho_descriptive": float(rho),
                "pvalue_descriptive_only": float(pvalue),
            }

    representative = InterlayerParameters()
    representative_result = locate_critical_mode(
        "M4", SETTINGS, interlayer=representative)
    representative_mode = canonical_soft_mode(
        representative_result["critical_wavevector_rlu"])

    summary = {
        "status": "deterministic local stress ensemble; not a fitted distribution",
        "public_constraints": {
            "source": "arXiv:2607.23173v1 lines 100-103 and 275-317",
            "experimental_delta_c_rlu": {"value": 0.027, "quoted_error": 0.002},
            "experimental_qz_rlu": {"value": 0.328, "quoted_error": 0.002},
            "source_model_qz_approximation": 1.0 / 3.0,
            "declared_delta_band_rlu": list(DELTA_BAND),
            "declared_qz_band_rlu": list(QZ_BAND),
            "km_line_residual_max_rlu": KM_LINE_RESIDUAL_MAX,
        },
        "stress_grid": {key: list(value) for key, value in GRID.items()},
        "grid_size": len(rows),
        "retained_size": len(retained),
        "coverage_fraction_descriptive_not_probability": len(retained) / len(rows),
        "M3_critical_field_T": m3_hc,
        "representative_reproduction": {
            **asdict(representative),
            "critical_field_T": representative_result["critical_field_T"],
            "critical_field_shift_vs_M3_T":
                representative_result["critical_field_T"] - m3_hc,
            **representative_mode,
            "out_of_plane_bandwidth_microeV":
                out_of_plane_bandwidth_microev(representative),
            "published_bandwidth_microeV": 4.4,
        },
        "retained_ranges": {
            "critical_field_shift_vs_M3_T":
                range_of(retained, "critical_field_shift_vs_M3_T"),
            "delta_inplane_rlu": range_of(retained, "delta_inplane_rlu"),
            "qz_folded_rlu": range_of(retained, "qz_folded_rlu"),
            "out_of_plane_bandwidth_microeV":
                range_of(retained, "out_of_plane_bandwidth_microeV"),
            "90mK_delta_B_star_T": range_of(e90, "delta_B_star_T"),
            "90mK_regret_percent": range_of(e90, "regret_percent"),
        },
        "threshold_coverage": {
            "retained_models_all_90mK_cases_delta_B_at_least_0p05T":
                all(row["delta_B_star_T"] >= 0.05 for row in e90),
            "retained_models_all_90mK_cases_regret_at_least_5percent":
                all(row["regret_percent"] >= 5.0 for row in e90),
            "90mK_scenarios_above_both_gates":
                sum(row["delta_B_star_T"] >= 0.05 and
                    row["regret_percent"] >= 5.0 for row in e90),
            "90mK_scenario_count": len(e90),
        },
        "descriptive_sensitivity_full_grid": sensitivity,
        "interpretation_limits": [
            "The parameter grid is a local deterministic stress envelope, not a confidence or posterior distribution.",
            "The public q constraints validate the near-saturation soft mode only; they do not validate the full finite-temperature M4 entropy.",
            "The engineering propagation inherits the public M0 heat-profile shape and height and remains a surrogate.",
            "The qz band includes the source paper's documented replacement of 0.328(2) by the commensurate value 1/3.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    write_csv(RESULTS / "interlayer_parameter_ensemble.csv", rows)
    write_csv(RESULTS / "interlayer_parameter_ensemble_retained.csv", retained)
    write_csv(RESULTS / "interlayer_parameter_ensemble_engineering.csv",
              engineering)
    (RESULTS / "interlayer_parameter_ensemble_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    make_figure(rows, engineering)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
