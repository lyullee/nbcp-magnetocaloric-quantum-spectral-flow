"""Propagate correlated shape templates measured in the N=12 M4 audit.

The independent outer-bound script combines every coordinate extreme and is
kept as an adversarial destruction test.  Here each peak residual, heat ratio,
and pair of left/right widths stays tied to the same microscopic M4 model.
That preserves the correlations actually observed in the 12-model coverage
set while retaining both experimental critical-field anchors and both field
mapping rules.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from bulk_to_engineering_surrogate import transform_field
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
CURVES = RESULTS / "m4_finite_temperature_cluster_curves.csv"
MODELS = RESULTS / "m4_finite_temperature_cluster_models.csv"
VERIFICATION = RESULTS / "m4_finite_temperature_cluster_verification.json"
OUTER = RESULTS / "thermodynamic_shape_outer_bound_summary.json"
ENSEMBLE = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"

TEMPERATURES = (0.09, 0.10, 0.12, 0.15)
EXPERIMENTAL_HC = (1.65, 1.70)
METHODS = ("reduced_field", "critical_offset")
BASE_SMOOTHING_MJ_PER_CM3 = (0.005, 0.0075, 0.010, 0.0125, 0.015)
WIDTH_CONTRACTION_FLOOR = 0.8
PARAMETERS = ("j2_microev", "ja_microev", "jb_microev", "delta_c")


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def smooth_metrics(rows: list[dict], smoothing: float) -> np.ndarray:
    ordered = sorted(rows, key=lambda row: float(row["field_T"]))
    by_field = {float(row["field_T"]): float(row["Q_mJ_per_cm3"])
                for row in ordered}
    field = np.array(sorted(by_field))
    heat = np.array([by_field[value] for value in field])
    mask = (field >= 0.85) & (field <= 1.55)
    field = field[mask]
    heat = heat[mask]
    spline = UnivariateSpline(
        field, heat, k=3, s=len(field) * smoothing**2)
    dense = np.linspace(float(field.min()), float(field.max()), 14001)
    values = spline(dense)
    peak = int(np.argmax(values))
    accepted = values >= 0.95 * values[peak]
    low = peak
    high = peak
    while low > 0 and accepted[low - 1]:
        low -= 1
    while high + 1 < len(dense) and accepted[high + 1]:
        high += 1
    return np.array([
        dense[peak], values[peak], dense[peak] - dense[low],
        dense[high] - dense[peak]], dtype=float)


def build_templates(curves: list[dict], models: list[dict],
                    verification: dict, temperature: float) -> list[dict]:
    metrics = {}
    labels = sorted({row["model_label"] for row in curves})
    scale = temperature / 0.09
    for label in labels:
        selected = [row for row in curves
                    if row["model_label"] == label
                    and np.isclose(float(row["temperature_K"]), temperature)]
        repeated = np.array([
            smooth_metrics(selected, smoothing * scale)
            for smoothing in BASE_SMOOTHING_MJ_PER_CM3])
        metrics[label] = np.median(repeated, axis=0)

    baseline = metrics["M3"]
    templates = []
    for model in models:
        if model["model_label"] == "M3":
            continue
        label = model["model_label"]
        value = metrics[label]
        template = {
            "temperature_K": temperature,
            "model_label": label,
            "sample_id": model["sample_id"],
            **{key: float(model[key]) for key in PARAMETERS},
            "bulk_critical_field_shift_T":
                float(model["critical_field_shift_vs_M3_T"]),
            "cluster_B_star_shift_T": float(value[0] - baseline[0]),
            "residual_peak_shift_after_bulk_alignment_T":
                float(value[0] - baseline[0]
                      - float(model["critical_field_shift_vs_M3_T"])),
            "heat_scale_to_M3": float(value[1] / baseline[1]),
            "left_width_scale_to_M3": float(max(
                WIDTH_CONTRACTION_FLOOR, value[2] / baseline[2])),
            "right_width_scale_to_M3": float(max(
                WIDTH_CONTRACTION_FLOOR, value[3] / baseline[3])),
            "smoothing_ensemble_count": len(BASE_SMOOTHING_MJ_PER_CM3),
        }
        if (np.isclose(temperature, 0.09)
                and label == verification["worst_primary_model"]):
            summary = verification["summary"]
            template.update({
                "cluster_B_star_shift_T":
                    summary["M4_minus_M3_B_star_T"]["mean"],
                "residual_peak_shift_after_bulk_alignment_T":
                    summary["M4_minus_M3_B_star_T"]["mean"]
                    - float(model["critical_field_shift_vs_M3_T"]),
                "heat_scale_to_M3":
                    summary["M4_to_M3_Qmax_ratio"]["mean"],
                "left_width_scale_to_M3":
                    summary["M4_to_M3_W95_left_ratio"]["mean"],
                "right_width_scale_to_M3":
                    summary["M4_to_M3_W95_right_ratio"]["mean"],
                "worst_model_common_random_replacement": True,
            })
        else:
            template["worst_model_common_random_replacement"] = False
        templates.append(template)
    return templates


def nearest_template(row: dict, templates: list[dict], spans: dict) -> dict:
    point = np.array([float(row[key]) / spans[key] for key in PARAMETERS])
    distances = []
    for template in templates:
        candidate = np.array([
            float(template[key]) / spans[key] for key in PARAMETERS])
        distances.append(float(np.linalg.norm(point - candidate)))
    return templates[int(np.argmin(distances))]


def warp(fields: np.ndarray, heat: np.ndarray, source_hc: float,
         target_hc: float, method: str, template: dict):
    mapped = transform_field(fields, source_hc, target_hc, method)
    peak = int(np.argmax(heat))
    center = float(mapped[peak])
    delta = mapped - center
    warped = (
        center + template["residual_peak_shift_after_bulk_alignment_T"]
        + np.where(delta < 0.0,
                   template["left_width_scale_to_M3"] * delta,
                   template["right_width_scale_to_M3"] * delta))
    return warped, template["heat_scale_to_M3"] * heat


def main() -> None:
    curves = read_csv(CURVES)
    models = read_csv(MODELS)
    verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
    outer = json.loads(OUTER.read_text(encoding="utf-8"))
    ensemble = read_csv(ENSEMBLE)
    values = {key: np.array([float(row[key]) for row in ensemble])
              for key in PARAMETERS}
    spans = {key: max(float(np.ptp(array)), 1.0e-12)
             for key, array in values.items()}
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)
    source = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",", names=True, encoding="utf-8-sig")
    source = source[source["anchor_N"] == 12]

    all_templates = []
    robust_rows = []
    gate_rows = []
    actual_scenario_rows = []
    envelope90 = None
    for temperature in TEMPERATURES:
        templates = build_templates(curves, models, verification, temperature)
        all_templates.extend(templates)
        profile = np.sort(source[np.isclose(source["temperature_K"], temperature)],
                          order="field_T")
        fields = profile["field_T"]
        heat = profile["Q_proxy_mJ_per_cm3"]
        source_peak = float(fields[int(np.argmax(heat))])

        scenarios = []
        grid = np.linspace(1.05, 1.75, 14001)
        nominal = {
            "model_label": "M3-nominal",
            "residual_peak_shift_after_bulk_alignment_T": 0.0,
            "heat_scale_to_M3": 1.0,
            "left_width_scale_to_M3": 1.0,
            "right_width_scale_to_M3": 1.0,
        }
        for anchor in EXPERIMENTAL_HC:
            for method in METHODS:
                for template in [nominal, *templates]:
                    warped_field, warped_heat = warp(
                        fields, heat, source_hc, anchor, method, template)
                    interpolated = np.interp(
                        grid, warped_field, warped_heat, left=np.nan, right=np.nan)
                    qmax = float(np.max(warped_heat))
                    peak_field = float(warped_field[int(np.argmax(warped_heat))])
                    actual_scenario_rows.append({
                        "temperature_K": temperature,
                        "experimental_Hc_anchor_T": anchor,
                        "mapping_method": method,
                        "template_model": template["model_label"],
                        "peak_field_T": peak_field,
                        "Qmax_mJ_per_cm3": qmax,
                        "residual_peak_shift_after_bulk_alignment_T":
                            template["residual_peak_shift_after_bulk_alignment_T"],
                        "heat_scale_to_M3": template["heat_scale_to_M3"],
                        "left_width_scale_to_M3": template["left_width_scale_to_M3"],
                        "right_width_scale_to_M3": template["right_width_scale_to_M3"],
                    })
                    scenarios.append({
                        "anchor": anchor,
                        "method": method,
                        "template": template,
                        "values": interpolated,
                        "regret": (qmax - interpolated) / qmax,
                    })
        regret = np.asarray([scenario["regret"] for scenario in scenarios])
        q_values = np.asarray([scenario["values"] for scenario in scenarios])
        valid = np.all(np.isfinite(regret), axis=0)
        worst = np.nanmax(regret, axis=0)
        indices = np.where(valid)[0]
        best = indices[int(np.argmin(worst[valid]))]
        guaranteed = valid & (worst <= 0.05)
        worst_q = float(np.nanmin(q_values[:, best]))
        robust_rows.append({
            "temperature_K": temperature,
            "scenario_count": len(scenarios),
            "robust_field_T": float(grid[best]),
            "worst_case_regret_percent": float(100.0 * worst[best]),
            "guaranteed_W95_exists": bool(np.any(guaranteed)),
            "guaranteed_W95_min_T": (
                float(grid[guaranteed].min()) if np.any(guaranteed) else None),
            "guaranteed_W95_max_T": (
                float(grid[guaranteed].max()) if np.any(guaranteed) else None),
            "worst_case_Q_at_robust_field_mJ_per_cm3": worst_q,
            "conditional_volume_cm3_for_0p10uW_60min": 0.36 / worst_q,
        })
        if np.isclose(temperature, 0.09):
            envelope90 = {
                "field_T": grid,
                "minimum_Q_mJ_per_cm3": np.nanmin(q_values, axis=0),
                "maximum_Q_mJ_per_cm3": np.nanmax(q_values, axis=0),
            }

        # Each retained model receives its nearest correlated coverage
        # template.  This avoids the impossible Cartesian combinations in the
        # independent outer-bound audit.
        for row in ensemble:
            template = nearest_template(row, templates, spans)
            shift = float(row["critical_field_shift_vs_M3_T"])
            for anchor in EXPERIMENTAL_HC:
                for method in METHODS:
                    p0_hc = anchor - shift
                    p0_field = float(transform_field(
                        source_peak, source_hc, p0_hc, method))
                    warped_field, warped_heat = warp(
                        fields, heat, source_hc, anchor, method, template)
                    p2_peak = float(warped_field[int(np.argmax(warped_heat))])
                    at_p0 = float(np.interp(
                        p0_field, warped_field, warped_heat,
                        left=np.nan, right=np.nan))
                    qmax = float(np.max(warped_heat))
                    gate_rows.append({
                        "temperature_K": temperature,
                        "sample_id": row["sample_id"],
                        "assigned_template": template["model_label"],
                        "template_distance_normalized": float(np.linalg.norm(
                            np.array([float(row[key]) / spans[key]
                                      for key in PARAMETERS])
                            - np.array([float(template[key]) / spans[key]
                                        for key in PARAMETERS]))),
                        "critical_field_shift_vs_M3_T": shift,
                        "experimental_Hc_anchor_T": anchor,
                        "mapping_method": method,
                        "P0_selected_field_T": p0_field,
                        "P2_optimal_field_T": p2_peak,
                        "delta_B_star_T": p2_peak - p0_field,
                        "Q_at_P0_field_mJ_per_cm3": at_p0,
                        "Qmax_mJ_per_cm3": qmax,
                        "regret_percent": 100.0 * (qmax - at_p0) / qmax,
                    })

    gate90 = [row for row in gate_rows
              if np.isclose(row["temperature_K"], 0.09)]
    actual90 = [row for row in actual_scenario_rows
                if np.isclose(row["temperature_K"], 0.09)]
    peak_support = [min(row["peak_field_T"] for row in actual90),
                    max(row["peak_field_T"] for row in actual90)]
    sweep_start = np.floor(peak_support[0] / 0.004) * 0.004
    sweep_stop = np.ceil(peak_support[1] / 0.004) * 0.004
    sweep_fields = np.arange(sweep_start, sweep_stop + 0.002, 0.004)
    sweep_rows = [{
        "point_index": index + 1,
        "field_T": float(value),
        "purpose": "localize the 90 mK peak within 0.002 T over the finite-cluster-calibrated support",
    } for index, value in enumerate(sweep_fields)]
    finite90 = [row for row in gate90 if np.isfinite(row["regret_percent"])]
    pass_both = [row for row in finite90
                 if row["delta_B_star_T"] >= 0.05
                 and row["regret_percent"] >= 5.0]
    payload = {
        "status": "correlated microscopic-template stress test; finite-cluster calibrated, not a converged bulk M4 solution",
        "template_definition": {
            "coverage_model_count": len([row for row in models
                                          if row["model_label"] != "M3"]),
            "smoothing_values_at_90mK_mJ_per_cm3":
                list(BASE_SMOOTHING_MJ_PER_CM3),
            "smoothing_scales_linearly_with_temperature": True,
            "width_contraction_floor": WIDTH_CONTRACTION_FLOOR,
            "worst_90mK_model_replaced_by_common_random_verification":
                verification["worst_primary_model"],
            "coordinates_remain_correlated_within_each_model": True,
        },
        "90mK_template_ranges": {
            key: [min(float(row[key]) for row in all_templates
                      if np.isclose(row["temperature_K"], 0.09)),
                  max(float(row[key]) for row in all_templates
                      if np.isclose(row["temperature_K"], 0.09))]
            for key in (
                "residual_peak_shift_after_bulk_alignment_T",
                "heat_scale_to_M3",
                "left_width_scale_to_M3",
                "right_width_scale_to_M3")
        },
        "90mK_gate_stress": {
            "scenario_count": len(gate90),
            "finite_scenario_count": len(finite90),
            "delta_B_star_T_range": [
                min(row["delta_B_star_T"] for row in finite90),
                max(row["delta_B_star_T"] for row in finite90)],
            "regret_percent_range": [
                min(row["regret_percent"] for row in finite90),
                max(row["regret_percent"] for row in finite90)],
            "fraction_passing_both_descriptive_not_probability":
                len(pass_both) / len(finite90),
            "all_finite_scenarios_pass_delta_B_gate": all(
                row["delta_B_star_T"] >= 0.05 for row in finite90),
            "all_finite_scenarios_pass_regret_gate": all(
                row["regret_percent"] >= 5.0 for row in finite90),
        },
        "90mK_actual_material_peak_support_T": peak_support,
        "90mK_no_fit_field_sweep": {
            "start_T": float(sweep_fields[0]),
            "stop_T": float(sweep_fields[-1]),
            "spacing_T": 0.004,
            "point_count": len(sweep_rows),
            "maximum_distance_to_nearest_grid_point_T": 0.002,
            "scope": "peak-localization grid only; it does not by itself validate profile width or heat capacity",
        },
        "actual_material_robust_decisions": robust_rows,
        "independent_outer_bound": {
            "status": outer["status"],
            "90mK_robust_decision": outer["actual_material_robust_decisions"][0],
            "interpretation": "The outer bound independently combines all coordinate extremes and is retained as a destruction test, not the calibrated design ensemble.",
        },
        "limitations": [
            "The templates inherit the N=12 short-circumference seam and FTLM trace uncertainty.",
            "The width-contraction floor retains the prior size-based lower bound because apparent stronger contractions were not stable under common-random verification.",
            "Nearest-template assignment covers the 69 retained models but is not interpolation of a converged 4-D thermodynamic surface.",
            "The reported heat and volume remain magnetic-only conditional budgets.",
        ],
    }
    write_csv(RESULTS / "thermodynamic_shape_templates.csv", all_templates)
    write_csv(RESULTS / "thermodynamic_shape_gate_stress.csv", gate_rows)
    write_csv(RESULTS / "thermodynamic_shape_actual_material_scenarios.csv",
              actual_scenario_rows)
    write_csv(RESULTS / "recommended_90mK_microscopic_field_sweep.csv",
              sweep_rows)
    write_csv(RESULTS / "thermodynamic_shape_robust_decisions.csv", robust_rows)
    (RESULTS / "thermodynamic_shape_stress_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    # Canonical engineering figure now shows the correlated template envelope.
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    axes[0].fill_between(
        envelope90["field_T"], envelope90["minimum_Q_mJ_per_cm3"],
        envelope90["maximum_Q_mJ_per_cm3"], color="#4c78a8", alpha=0.35,
        label="correlated M4 templates")
    axes[0].axvline(robust_rows[0]["robust_field_T"], color="black",
                    ls="--", label="minimax field")
    axes[0].set(xlabel="Field (T)", ylabel=r"$Q$ (mJ cm$^{-3}$)",
                title="90 mK microscopic-template envelope")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].scatter(
        [row["delta_B_star_T"] for row in finite90],
        [row["regret_percent"] for row in finite90], s=10, alpha=0.35,
        color="#f58518")
    axes[1].axvline(0.05, color="#e45756", ls="--")
    axes[1].axhline(5.0, color="#e45756", ls=":")
    axes[1].set(xlabel=r"90 mK $\Delta B^*$ (T)", ylabel="regret (%)",
                title="69-model assigned-template gate")

    axes[2].plot(
        [1000 * row["temperature_K"] for row in robust_rows],
        [row["robust_field_T"] for row in robust_rows], "o-",
        color="#54a24b", label="robust field")
    secondary = axes[2].twinx()
    secondary.plot(
        [1000 * row["temperature_K"] for row in robust_rows],
        [row["worst_case_regret_percent"] for row in robust_rows], "s--",
        color="#e45756", label="worst regret")
    axes[2].set(xlabel="Temperature (mK)", ylabel="Robust field (T)",
                title="Calibrated minimax decision")
    secondary.set_ylabel("Worst regret (%)")
    lines = axes[2].get_lines() + secondary.get_lines()
    axes[2].legend(lines, [line.get_label() for line in lines],
                   frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / "thermodynamic_shape_stress.png", dpi=220)
    plt.close(fig)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
