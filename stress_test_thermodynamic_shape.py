"""Stress-test the engineering decision against finite-T profile deformation.

The existing critical-field bridge shifts a public XTRG heat profile without
changing its shape.  This script adds a declared, data-scaled deformation
envelope: independent left/right width changes, a residual peak displacement,
and a heat-scale change.  It does not claim to solve the finite-temperature M4
Hamiltonian.  Its purpose is to test whether the design decision survives a
broader family than two rigid horizontal mappings.
"""
from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bulk_to_engineering_surrogate import transform_field
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

TEMPERATURES = (0.09, 0.10, 0.12, 0.15)
EXPERIMENTAL_HC = (1.65, 1.70)
METHODS = ("reduced_field", "critical_offset")
# The original size-only bounds were +/-0.02 T and 0.8--1.2 in width.
# An explicit 3 x 2 x 2 M4 finite-temperature audit showed a reproducible
# +0.10 T differential optimum shift and a 1.46--2.05 total W95 expansion
# for the strongest retained model.  The wider coordinates below contain the
# full 12-model coverage audit plus the three-seed common-random verification.
PEAK_SHIFT_T = (-0.12, 0.0, 0.12)
LEFT_WIDTH_SCALE = (0.5, 1.0, 1.7)
RIGHT_WIDTH_SCALE = (0.5, 1.0, 2.5)
HEAT_SCALE = (0.9, 1.0, 1.1)


def deformed_profile(fields, heat, source_hc, target_hc, method,
                     peak_shift_t, left_scale, right_scale, heat_scale):
    mapped = transform_field(fields, source_hc, target_hc, method)
    source_peak_index = int(np.argmax(heat))
    center = float(mapped[source_peak_index])
    delta = mapped - center
    warped = center + peak_shift_t + np.where(
        delta < 0.0, left_scale * delta, right_scale * delta)
    return warped, heat_scale * heat


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def profile_at_temperature(profile, temperature):
    p = np.sort(profile[np.isclose(profile["temperature_K"], temperature)],
                order="field_T")
    return p["field_T"], p["Q_proxy_mJ_per_cm3"]


def main() -> None:
    microscopic = json.loads((
        RESULTS / "m4_finite_temperature_cluster_summary.json"
    ).read_text(encoding="utf-8"))
    microscopic_verification = json.loads((
        RESULTS / "m4_finite_temperature_cluster_verification.json"
    ).read_text(encoding="utf-8"))
    adaptive_path = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"
    ensemble_path = (
        adaptive_path if adaptive_path.exists()
        else RESULTS / "interlayer_parameter_ensemble_retained.csv"
    )
    ensemble = list(csv.DictReader(
        ensemble_path.open(encoding="utf-8-sig")))
    if not ensemble:
        raise RuntimeError("the public-q constrained interlayer set is empty")
    shifts = sorted({float(row["critical_field_shift_vs_M3_T"])
                     for row in ensemble})
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)
    profile = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",", names=True, encoding="utf-8-sig")
    profile = profile[profile["anchor_N"] == 12]

    deformation_grid = list(itertools.product(
        PEAK_SHIFT_T, LEFT_WIDTH_SCALE, RIGHT_WIDTH_SCALE, HEAT_SCALE))
    gate_rows = []
    robust_rows = []
    curve_rows_90 = []
    envelope90 = None

    for temperature in TEMPERATURES:
        fields, heat = profile_at_temperature(profile, temperature)
        source_bstar = float(fields[int(np.argmax(heat))])
        qmax_source = float(np.max(heat))

        # Test the P2-P0 engineering gate.  The profile deformation is applied
        # to P2 as an explicit differential M4 uncertainty; P0 remains the
        # public reference transported to its counterfactual critical field.
        for shift, anchor, method, deformation in itertools.product(
                shifts, EXPERIMENTAL_HC, METHODS, deformation_grid):
            peak_shift, left_scale, right_scale, heat_scale = deformation
            p0_hc = anchor - shift
            p0_field = float(transform_field(source_bstar, source_hc,
                                             p0_hc, method))
            warped_fields, warped_heat = deformed_profile(
                fields, heat, source_hc, anchor, method, peak_shift,
                left_scale, right_scale, heat_scale)
            p2_index = int(np.argmax(warped_heat))
            p2_field = float(warped_fields[p2_index])
            q_at_p0 = float(np.interp(p0_field, warped_fields, warped_heat,
                                      left=np.nan, right=np.nan))
            qmax = float(np.max(warped_heat))
            gate_rows.append({
                "temperature_K": temperature,
                "critical_field_shift_vs_M3_T": shift,
                "experimental_Hc_anchor_T": anchor,
                "mapping_method": method,
                "residual_peak_shift_T": peak_shift,
                "left_width_scale": left_scale,
                "right_width_scale": right_scale,
                "heat_scale": heat_scale,
                "P0_selected_field_T": p0_field,
                "P2_optimal_field_T": p2_field,
                "delta_B_star_T": p2_field - p0_field,
                "Q_at_P0_field_mJ_per_cm3": q_at_p0,
                "Qmax_mJ_per_cm3": qmax,
                "regret_percent": 100.0 * (qmax - q_at_p0) / qmax,
            })

        # Actual-material minimax decision across experimental Hc, mapping,
        # and all deformation coordinates.  Interlayer shift is absent here
        # because measured Hc already contains the real interactions.
        scenarios = []
        grid = np.linspace(1.05, 1.75, 14001)
        for anchor, method, deformation in itertools.product(
                EXPERIMENTAL_HC, METHODS, deformation_grid):
            peak_shift, left_scale, right_scale, heat_scale = deformation
            warped_fields, warped_heat = deformed_profile(
                fields, heat, source_hc, anchor, method, peak_shift,
                left_scale, right_scale, heat_scale)
            values = np.interp(grid, warped_fields, warped_heat,
                               left=np.nan, right=np.nan)
            qmax = float(np.max(warped_heat))
            scenarios.append({
                "anchor": anchor,
                "method": method,
                "peak_shift": peak_shift,
                "left_scale": left_scale,
                "right_scale": right_scale,
                "heat_scale": heat_scale,
                "values": values,
                "regret": (qmax - values) / qmax,
            })
            if temperature == 0.09 and heat_scale == 1.0:
                for b, q in zip(warped_fields, warped_heat):
                    if 1.05 <= b <= 1.75:
                        curve_rows_90.append({
                            "experimental_Hc_anchor_T": anchor,
                            "mapping_method": method,
                            "residual_peak_shift_T": peak_shift,
                            "left_width_scale": left_scale,
                            "right_width_scale": right_scale,
                            "field_T": float(b),
                            "Q_mJ_per_cm3": float(q),
                        })
        regrets = np.asarray([s["regret"] for s in scenarios])
        values = np.asarray([s["values"] for s in scenarios])
        if temperature == 0.09:
            envelope90 = {
                "field_T": grid.copy(),
                "minimum_Q_mJ_per_cm3": np.nanmin(values, axis=0),
                "maximum_Q_mJ_per_cm3": np.nanmax(values, axis=0),
            }
        valid = np.all(np.isfinite(regrets), axis=0)
        worst_regret = np.nanmax(regrets, axis=0)
        valid_indices = np.where(valid)[0]
        best_index = valid_indices[np.argmin(worst_regret[valid])]
        guaranteed = valid & (worst_regret <= 0.05)
        robust_field = float(grid[best_index])
        worst_q = float(np.nanmin(values[:, best_index]))
        required_mj = 0.10e-3 * 60.0 * 60.0
        robust_rows.append({
            "temperature_K": temperature,
            "scenario_count": len(scenarios),
            "robust_field_T": robust_field,
            "worst_case_regret_percent":
                float(100.0 * worst_regret[best_index]),
            "guaranteed_W95_exists": bool(np.any(guaranteed)),
            "guaranteed_W95_min_T":
                float(grid[guaranteed].min()) if np.any(guaranteed) else None,
            "guaranteed_W95_max_T":
                float(grid[guaranteed].max()) if np.any(guaranteed) else None,
            "worst_case_Q_at_robust_field_mJ_per_cm3": worst_q,
            "conditional_volume_cm3_for_0p10uW_60min": required_mj / worst_q,
        })

    gate90 = [row for row in gate_rows if row["temperature_K"] == 0.09]
    finite90 = [row for row in gate90 if np.isfinite(row["regret_percent"])]
    def passed(row):
        return row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0

    by_shift = []
    for shift in shifts:
        subset = [row for row in finite90
                  if np.isclose(row["critical_field_shift_vs_M3_T"], shift)]
        by_shift.append({
            "critical_field_shift_vs_M3_T": shift,
            "scenario_count": len(subset),
            "pass_fraction_descriptive_not_probability":
                sum(passed(row) for row in subset) / len(subset),
            "minimum_delta_B_star_T":
                min(row["delta_B_star_T"] for row in subset),
            "minimum_regret_percent":
                min(row["regret_percent"] for row in subset),
            "all_shape_scenarios_pass": all(passed(row) for row in subset),
        })
    central = [row for row in finite90
               if row["residual_peak_shift_T"] == 0.0
               and row["left_width_scale"] == 1.0
               and row["right_width_scale"] == 1.0
               and row["heat_scale"] == 1.0]
    shifts_all_pass = [
        row["critical_field_shift_vs_M3_T"] for row in by_shift
        if row["all_shape_scenarios_pass"]]
    summary = {
        "status": "finite-cluster-calibrated shape-deformation stress test; not a converged bulk M4 solution",
        "interlayer_ensemble_source": ensemble_path.name,
        "interlayer_model_count": len(ensemble),
        "basis_for_bounds": {
            "N9_vs_N12_90mK_Qmax_ratio": 0.4971378441231503 / 0.4593937417224646,
            "N9_vs_N12_90mK_W95_width_ratio": 0.08 / 0.07,
            "heat_scale": list(HEAT_SCALE),
            "left_width_scale": list(LEFT_WIDTH_SCALE),
            "right_width_scale": list(RIGHT_WIDTH_SCALE),
            "residual_peak_shift_T": list(PEAK_SHIFT_T),
            "previous_size_only_bounds": {
                "heat_scale": [0.9, 1.0, 1.1],
                "left_width_scale": [0.8, 1.0, 1.2],
                "right_width_scale": [0.8, 1.0, 1.2],
                "residual_peak_shift_T": [-0.02, 0.0, 0.02],
            },
            "microscopic_cluster_status": microscopic["status"],
            "microscopic_90mK_ranges": microscopic["90mK_differential_ranges"],
            "common_random_verification_status": microscopic_verification["status"],
            "common_random_verification_summary": microscopic_verification["summary"],
            "note": "The heat scale remains 0.9--1.1 because all 12 coverage models and all three common-random repeats remain inside it. Width and residual-shift coordinates are enlarged to contain the N=12 M4 sensitivity and its repeated worst-model audit.",
        },
        "public_q_constrained_interlayer_shift_range_T":
            [min(shifts), max(shifts)],
        "90mK_gate_stress": {
            "scenario_count": len(gate90),
            "finite_scenario_count": len(finite90),
            "delta_B_star_T_range": [
                min(row["delta_B_star_T"] for row in finite90),
                max(row["delta_B_star_T"] for row in finite90),
            ],
            "regret_percent_range": [
                min(row["regret_percent"] for row in finite90),
                max(row["regret_percent"] for row in finite90),
            ],
            "all_finite_scenarios_pass_delta_B_gate":
                all(row["delta_B_star_T"] >= 0.05 for row in finite90),
            "all_finite_scenarios_pass_regret_gate":
                all(row["regret_percent"] >= 5.0 for row in finite90),
            "fraction_passing_both_descriptive_not_probability":
                sum(row["delta_B_star_T"] >= 0.05 and
                    row["regret_percent"] >= 5.0 for row in finite90)
                / len(finite90),
            "by_interlayer_shift": by_shift,
            "minimum_interlayer_shift_with_all_shape_scenarios_passing_T":
                min(shifts_all_pass) if shifts_all_pass else None,
            "any_interlayer_model_all_shape_scenarios_passing":
                bool(shifts_all_pass),
            "central_shape_all_interlayer_models_pass_both":
                all(passed(row) for row in central),
            "central_shape_minimum_delta_B_star_T":
                min(row["delta_B_star_T"] for row in central),
            "central_shape_minimum_regret_percent":
                min(row["regret_percent"] for row in central),
        },
        "actual_material_robust_decisions": robust_rows,
        "limitations": [
            "The deformation coordinates are deterministic stress bounds, not statistical confidence intervals.",
            "The microscopic calibration is a 3 x 2 x 2 finite-cluster sensitivity with a short-circumference seam, not a converged bulk M4 thermodynamic solution.",
            "The envelope tests peak displacement, asymmetric width and heat-scale changes but cannot represent every entropy redistribution allowed by M4.",
            "The P2-P0 gate treats the deformation as differential uncertainty on P2; the actual-material recommendation uses measured Hc and does not require a counterfactual P0 model.",
            "All heat and volume values remain magnetic-only conditional budgets.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    write_csv(RESULTS / "thermodynamic_shape_outer_bound_gate.csv", gate_rows)
    write_csv(RESULTS / "thermodynamic_shape_outer_bound_90mK_curves.csv", curve_rows_90)
    write_csv(RESULTS / "thermodynamic_shape_outer_bound_decisions.csv", robust_rows)
    (RESULTS / "thermodynamic_shape_outer_bound_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    # Compact audit figure.
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    for row in curve_rows_90[::max(1, len(curve_rows_90) // 120)]:
        pass  # rows are pointwise; envelope is plotted below from binned data
    axes[0].fill_between(envelope90["field_T"],
                         envelope90["minimum_Q_mJ_per_cm3"],
                         envelope90["maximum_Q_mJ_per_cm3"],
                         color="#4c78a8", alpha=0.35,
                         label="shape envelope")
    axes[0].axvline(robust_rows[0]["robust_field_T"], color="black",
                    linestyle="--", label="minimax field")
    axes[0].set(xlabel="Field (T)", ylabel=r"$Q$ (mJ cm$^{-3}$)",
                title="90 mK deformation envelope")
    axes[0].legend(fontsize=8)

    axes[1].scatter([r["delta_B_star_T"] for r in finite90],
                    [r["regret_percent"] for r in finite90], s=8,
                    alpha=0.25, color="#f58518")
    axes[1].axvline(0.05, color="#e45756", linestyle="--")
    axes[1].axhline(5.0, color="#e45756", linestyle=":")
    axes[1].set(xlabel=r"90 mK $\Delta B^*$ (T)", ylabel="regret (%)",
                title="Gate under shape stress")

    axes[2].plot([1000 * r["temperature_K"] for r in robust_rows],
                 [r["robust_field_T"] for r in robust_rows], "o-",
                 color="#54a24b", label="robust field")
    ax2 = axes[2].twinx()
    ax2.plot([1000 * r["temperature_K"] for r in robust_rows],
             [r["worst_case_regret_percent"] for r in robust_rows], "s--",
             color="#e45756", label="worst regret")
    axes[2].set(xlabel="Temperature (mK)", ylabel="Robust field (T)",
                title="Actual-material minimax decision")
    ax2.set_ylabel("Worst regret (%)")
    lines = axes[2].get_lines() + ax2.get_lines()
    axes[2].legend(lines, [line.get_label() for line in lines], fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "thermodynamic_shape_outer_bound.png", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
