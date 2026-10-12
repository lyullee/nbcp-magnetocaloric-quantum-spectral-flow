"""Field-dependent fan-curve holdout for every retained M4 model.

The near-saturation propagation vector was already used to retain the M4
models.  The four lower-field delta(H) points were not used in that selection
and therefore provide a genuine multi-point holdout.  This script reconstructs the
fixed-length fan ansatz of arXiv:2607.23173v1, Eqs. (6)--(7), for all 69 models
from the adaptive inverse audit.

The public source itself states that its representative mean-field model
captures the trend but not the rapid experimental evolution.  Consequently,
the relative-to-representative subset reported here is descriptive model
discrimination, not a confidence region or an experimental goodness-of-fit
acceptance test.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from bulk_polarized_stability import (
    EwaldSettings,
    GZ,
    InterlayerParameters,
    MICROEV_TO_K,
    S,
    coupling_fourier,
)
from pilot_fan_curve_holdout import (
    EXPERIMENTAL_DELTA_RLU,
    EXPERIMENTAL_FIELD_T,
    EXPERIMENTAL_HC_T,
)
from screen_interlayer_parameter_ensemble import engineering_rows
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCES = ROOT / "sources"
SETTINGS = EwaldSettings(boundary="tin_foil")
RETAINED_PATH = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"

DELTA_GRID = np.linspace(0.0, 0.040, 321)
CURVE_RATIOS = np.r_[np.linspace(0.86, 0.98, 13), 0.99, 0.995, 0.999]
HOLDOUT_MASK = ((EXPERIMENTAL_FIELD_T >= 1.45)
                & (EXPERIMENTAL_FIELD_T < EXPERIMENTAL_HC_T))
HOLDOUT_RATIOS = EXPERIMENTAL_FIELD_T[HOLDOUT_MASK] / EXPERIMENTAL_HC_T
HOLDOUT_EVALUATION_RATIOS = np.minimum(HOLDOUT_RATIOS, 0.999)
HOLDOUT_DATA = EXPERIMENTAL_DELTA_RLU[HOLDOUT_MASK]


def read_retained() -> list[dict]:
    rows = []
    with RETAINED_PATH.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            rows.append({
                **raw,
                "sample_id": int(float(raw["sample_id"])),
                "j2_microev": float(raw["j2_microev"]),
                "ja_microev": float(raw["ja_microev"]),
                "jb_microev": float(raw["jb_microev"]),
                "delta_c": float(raw["delta_c"]),
                "critical_field_T": float(raw["critical_field_T"]),
                "critical_field_shift_vs_M3_T":
                    float(raw["critical_field_shift_vs_M3_T"]),
                "delta_inplane_rlu": float(raw["delta_inplane_rlu"]),
            })
    return rows


def interlayer_fourier(frac_k: np.ndarray, p: InterlayerParameters) -> np.ndarray:
    u, v, w = np.asarray(frac_k, dtype=float)
    terms = [((0, 0, 1), p.j2_microev)]
    terms += [(d, p.ja_microev) for d in
              ((1, 0, 1), (0, 1, 1), (1, -1, 1))]
    terms += [(d, p.jb_microev) for d in
              ((-1, 0, 1), (0, -1, 1), (-1, 1, 1))]
    result = np.zeros((3, 3), dtype=float)
    for d, j_microev in terms:
        phase = 2.0 * np.pi * (d[0] * u + d[1] * v + d[2] * w)
        tensor = np.diag([
            p.delta_c * j_microev,
            p.delta_c * j_microev,
            j_microev,
        ]) * MICROEV_TO_K
        result += 2.0 * np.cos(phase) * tensor
    return result


def precompute_m3() -> list[dict]:
    base = []
    j0 = coupling_fourier(np.zeros(3), "M3", SETTINGS)
    for delta in DELTA_GRID:
        q = np.array([2.0 / 3.0 + 2.0 * delta,
                      1.0 / 3.0 + delta, 2.0 / 3.0])
        base.append({
            "delta": float(delta),
            "q": q,
            "j0": j0,
            "jq": coupling_fourier(np.mod(q, 1.0), "M3", SETTINGS),
            "j2q": coupling_fourier(np.mod(2.0 * q, 1.0), "M3", SETTINGS),
            "j3q": coupling_fourier(np.mod(3.0 * q, 1.0), "M3", SETTINGS),
        })
    return base


def model_terms(base: list[dict], p: InterlayerParameters) -> list[dict]:
    terms = []
    for row in base:
        q = row["q"]
        j0 = row["j0"] + interlayer_fourier(np.zeros(3), p)
        jq = row["jq"] + interlayer_fourier(np.mod(q, 1.0), p)
        j2q = row["j2q"] + interlayer_fourier(np.mod(2.0 * q, 1.0), p)
        j3q = row["j3q"] + interlayer_fourier(np.mod(3.0 * q, 1.0), p)
        _, eigenvectors = np.linalg.eigh(jq[:2, :2])
        transverse = eigenvectors[:, 0]
        terms.append({
            "delta": row["delta"],
            "j0_zz": float(j0[2, 2]),
            "jq_perp": float(transverse @ jq[:2, :2] @ transverse),
            "j2q_zz": float(j2q[2, 2]),
            "j3q_perp": float(transverse @ j3q[:2, :2] @ transverse),
        })
    return terms


def minimum_cubic(coefficients: tuple[float, float, float, float],
                  y_max: float = 1.95**2) -> tuple[float, float]:
    c0, c1, c2, c3 = coefficients
    candidates = [0.0, y_max]
    roots = np.roots([3.0 * c3, 2.0 * c2, c1])
    candidates.extend(float(root.real) for root in roots
                      if abs(root.imag) < 1.0e-10 and 0.0 <= root.real <= y_max)
    energies = [c0 + c1*y + c2*y*y + c3*y*y*y for y in candidates]
    index = int(np.argmin(energies))
    return float(energies[index]), float(np.sqrt(candidates[index]))


def optimize_delta(terms: list[dict], field_T: float) -> tuple[float, float, float]:
    h = GZ * MU_B_OVER_KB * field_T
    candidates = []
    for row in terms:
        j0, jq = row["j0_zz"], row["jq_perp"]
        j2, j3 = row["j2q_zz"], row["j3q_perp"]
        coefficients = (
            0.5 * S*S*j0 - h*S,
            S*S * (-j0 + jq) / 4.0 + h*S / 4.0,
            S*S * (2.0*j0 - 3.0*jq + j2) / 64.0,
            S*S * (9.0*jq + j3) / 4096.0,
        )
        energy, xi = minimum_cubic(coefficients)
        candidates.append((energy, row["delta"], xi))
    return min(candidates)


def evaluate_model(raw: dict, base: list[dict]) -> tuple[dict, list[dict]]:
    p = InterlayerParameters(
        j2_microev=raw["j2_microev"], ja_microev=raw["ja_microev"],
        jb_microev=raw["jb_microev"], delta_c=raw["delta_c"],
    )
    terms = model_terms(base, p)
    curve = []
    for ratio in CURVE_RATIOS:
        energy, delta, xi = optimize_delta(terms, ratio * raw["critical_field_T"])
        curve.append({
            "sample_id": raw["sample_id"],
            "field_ratio": float(ratio),
            "delta_opt_rlu": delta,
            "xi_opt": xi,
            "energy_K_per_spin": energy,
        })
    holdout_prediction = []
    for ratio in HOLDOUT_EVALUATION_RATIOS:
        _, delta, _ = optimize_delta(terms, ratio * raw["critical_field_T"])
        holdout_prediction.append(delta)
    holdout_prediction = np.asarray(holdout_prediction)
    residual = holdout_prediction - HOLDOUT_DATA
    near_hc_delta = curve[-1]["delta_opt_rlu"]
    summary = {
        "sample_id": raw["sample_id"],
        **asdict(p),
        "critical_field_T": raw["critical_field_T"],
        "critical_field_shift_vs_M3_T": raw["critical_field_shift_vs_M3_T"],
        "soft_mode_delta_rlu": raw["delta_inplane_rlu"],
        "fan_delta_at_0p999_Hc_rlu": near_hc_delta,
        "near_hc_internal_difference_rlu":
            abs(near_hc_delta - raw["delta_inplane_rlu"]),
        "holdout_mae_rlu": float(np.mean(np.abs(residual))),
        "holdout_rmse_rlu": float(np.sqrt(np.mean(residual**2))),
        "holdout_max_abs_rlu": float(np.max(np.abs(residual))),
        "holdout_mean_signed_residual_rlu": float(np.mean(residual)),
        "holdout_prediction_rlu": holdout_prediction.tolist(),
        "holdout_residual_rlu": residual.tolist(),
        "monotone_fan_curve": bool(np.all(np.diff(
            [row["delta_opt_rlu"] for row in curve]
        ) >= -1.0e-12)),
    }
    return summary, curve


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        serializable = []
        for row in rows:
            serializable.append({
                key: json.dumps(value) if isinstance(value, list) else value
                for key, value in row.items()
            })
        writer = csv.DictWriter(handle, fieldnames=list(serializable[0]))
        writer.writeheader()
        writer.writerows(serializable)


def main() -> None:
    retained = read_retained()
    base = precompute_m3()
    summaries, curves = [], []
    for index, raw in enumerate(retained, start=1):
        summary, curve = evaluate_model(raw, base)
        summaries.append(summary)
        curves.extend(curve)
        if index % 10 == 0:
            print(f"completed {index}/{len(retained)}", flush=True)

    representative_raw = {
        "sample_id": -1,
        **asdict(InterlayerParameters()),
        "critical_field_T": 1.8676712357447807,
        "critical_field_shift_vs_M3_T": 0.0952738089552219,
        "delta_inplane_rlu": 0.024826409699414364,
    }
    representative, representative_curve = evaluate_model(representative_raw, base)
    benchmark = representative["holdout_mae_rlu"]
    no_worse = [row for row in summaries if row["holdout_mae_rlu"] <= benchmark]
    conservative = [row for row in summaries
                    if row["holdout_mae_rlu"] <= benchmark + 0.002]
    best = min(summaries, key=lambda row: row["holdout_mae_rlu"])
    rho, pvalue = spearmanr(
        [row["holdout_mae_rlu"] for row in summaries],
        [row["critical_field_shift_vs_M3_T"] for row in summaries],
    )
    no_worse_engineering_90 = [
        row for row in engineering_rows(no_worse, 0.0)
        if np.isclose(row["temperature_K"], 0.09)
    ]
    no_worse_central_pass = [
        row for row in no_worse_engineering_90
        if row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
    ]
    no_worse_shifts = np.array([
        row["critical_field_shift_vs_M3_T"] for row in no_worse
    ])
    with (RESULTS / "thermodynamic_shape_gate_stress.csv").open(
            newline="", encoding="utf-8-sig") as handle:
        no_worse_shape_90 = [
            row for row in csv.DictReader(handle)
            if np.isclose(float(row["temperature_K"]), 0.09)
            and float(np.min(np.abs(
                no_worse_shifts - float(row["critical_field_shift_vs_M3_T"])
            ))) < 1.0e-10
        ]
    no_worse_shape_pass = [
        row for row in no_worse_shape_90
        if float(row["delta_B_star_T"]) >= 0.05
        and float(row["regret_percent"]) >= 5.0
    ]

    def shift_range(rows: list[dict]) -> list[float] | None:
        if not rows:
            return None
        values = [row["critical_field_shift_vs_M3_T"] for row in rows]
        return [float(min(values)), float(max(values))]

    all_internal = all(row["near_hc_internal_difference_rlu"] <= 0.001
                       for row in summaries)
    all_monotone = all(row["monotone_fan_curve"] for row in summaries)
    result = {
        "status": "field-dependent fan holdout; relative subsets are descriptive, not confidence sets",
        "source": "arXiv:2607.23173v1 Eqs. (6)-(7) and source-vector Figure 7d",
        "retained_model_count": len(summaries),
        "holdout_fields_T": EXPERIMENTAL_FIELD_T[HOLDOUT_MASK].tolist(),
        "holdout_field_ratios": HOLDOUT_RATIOS.tolist(),
        "holdout_digitized_delta_rlu": HOLDOUT_DATA.tolist(),
        "representative_model": representative,
        "best_model": best,
        "holdout_mae_range_rlu": [
            float(min(row["holdout_mae_rlu"] for row in summaries)),
            float(max(row["holdout_mae_rlu"] for row in summaries)),
        ],
        "internal_validation": {
            "all_69_near_hc_fan_vs_soft_mode_within_0p001_rlu": all_internal,
            "maximum_near_hc_difference_rlu": float(max(
                row["near_hc_internal_difference_rlu"] for row in summaries
            )),
            "all_69_monotone_in_declared_fan_window": all_monotone,
        },
        "relative_to_published_representative": {
            "representative_holdout_mae_rlu": benchmark,
            "models_no_worse_count": len(no_worse),
            "models_no_worse_shift_range_T": shift_range(no_worse),
            "models_no_worse_90mK_central_both_gates": len(no_worse_central_pass),
            "models_no_worse_90mK_central_case_count": len(no_worse_engineering_90),
            "models_no_worse_90mK_shape_both_gates": len(no_worse_shape_pass),
            "models_no_worse_90mK_shape_case_count": len(no_worse_shape_90),
            "models_within_extra_0p002_rlu_count": len(conservative),
            "models_within_extra_0p002_rlu_shift_range_T": shift_range(conservative),
            "warning": (
                "These are relative performance subsets.  The representative "
                "model itself has a documented systematic mismatch and is not "
                "an experimental acceptance threshold."
            ),
        },
        "descriptive_association": {
            "spearman_holdout_mae_vs_critical_field_shift": float(rho),
            "pvalue_descriptive_only": float(pvalue),
        },
        "major_extension_decision": {
            "pilot_and_full_evaluation_feasible": bool(all_internal and all_monotone),
            "can_claim_quantitative_curve_fit": False,
            "can_claim_multi_point_model_discrimination": True,
            "next_required_component": (
                "add qz(H) and phase-boundary energies before using this holdout "
                "to replace the current 69-model engineering envelope"
            ),
        },
        "interpretation_limits": [
            "The source paper reports that the representative mean-field curve evolves more slowly than experiment.",
            "The digitized points are not a substitute for raw reduced diffraction data.",
            "The fan ansatz is evaluated only in the source-declared fan window; it does not determine the lower V--fan transition.",
            "The relative subsets must not be called posterior probabilities or confidence intervals.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    SOURCES.mkdir(exist_ok=True)
    write_csv(RESULTS / "fan_curve_holdout_models.csv", summaries)
    write_csv(RESULTS / "fan_curve_holdout_curves.csv", curves)
    write_csv(RESULTS / "fan_curve_holdout_no_worse_models.csv", no_worse)
    (RESULTS / "fan_curve_holdout_summary.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    provenance = {
        "source_article": "https://arxiv.org/abs/2607.23173",
        "source_eprint": "https://export.arxiv.org/e-print/2607.23173",
        "source_version": "v1, 2026-07-25",
        "source_elements": ["Eqs. (6)-(7)", "Figure 7d", "Table I"],
        "digitized_fields_T": EXPERIMENTAL_FIELD_T.tolist(),
        "digitized_delta_rlu": EXPERIMENTAL_DELTA_RLU.tolist(),
        "note": (
            "Values were read from the vector figure.  The 1.45--1.60 T fan "
            "points are scored as holdout; the 1.65 T point is reserved for "
            "the independent near-Hc consistency check."
        ),
    }
    (SOURCES / "fan_curve_holdout_provenance_20261012.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    ax = axes[0]
    shift_by_id = {row["sample_id"]: row["critical_field_shift_vs_M3_T"]
                   for row in summaries}
    for sample_id in sorted(shift_by_id):
        rows = [row for row in curves if row["sample_id"] == sample_id]
        scatter = ax.scatter(
            [row["field_ratio"] for row in rows],
            [row["delta_opt_rlu"] for row in rows],
            c=[shift_by_id[sample_id]] * len(rows), vmin=0.03, vmax=0.205,
            cmap="viridis", s=5, alpha=0.45,
        )
    ax.plot(EXPERIMENTAL_FIELD_T / EXPERIMENTAL_HC_T,
            EXPERIMENTAL_DELTA_RLU, "o-", color="#d62728", linewidth=1.6,
            markersize=6, label="public diffraction")
    ax.set(xlabel=r"reduced field $H/H_c$", ylabel=r"$\delta$ (r.l.u.)",
           title="Field-dependent holdout across 69 M4 models",
           xlim=(0.81, 1.005), ylim=(-0.001, 0.033))
    ax.legend(fontsize=8)
    fig.colorbar(scatter, ax=ax, label=r"M4--M3 $\Delta H_c$ (T)")

    ax = axes[1]
    ax.scatter([row["critical_field_shift_vs_M3_T"] for row in summaries],
               [row["holdout_mae_rlu"] for row in summaries],
               color="#4c78a8", s=30, alpha=0.75)
    ax.axhline(benchmark, color="black", linestyle="--",
               label="published-representative MAE")
    ax.scatter([best["critical_field_shift_vs_M3_T"]],
               [best["holdout_mae_rlu"]], marker="*", s=150,
               color="#f2cf5b", edgecolor="black", linewidth=0.5,
               label="best retained model")
    ax.set(xlabel=r"M4--M3 $\Delta H_c$ (T)",
           ylabel=r"four-point holdout MAE in $\delta$ (r.l.u.)",
           title="Holdout fit versus engineering field shift",
           ylim=(0.00410, 0.00540))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "fan_curve_holdout_audit.png", dpi=220)
    fig.savefig(FIGURES / "fan_curve_holdout_audit.pdf")
    fig.savefig(FIGURES / "fan_curve_holdout_audit.svg")
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
