"""Adversarial 90 mK design stress from the N=12 -> N=18 size escalation.

The N=18 calculation is not treated as a bulk prediction.  Its correlated
peak, heat and local-width deformation is appended to the existing N=12
scenario set as an explicit size-transfer stress.  The calculation asks
whether the previous operating field and small measurement panel survive.
"""
from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from optimize_measurement_value_of_information import (
    evaluate_panel,
    evaluate_masks,
)
from screen_model import MU_B_OVER_KB
from stress_test_microscopic_shape_templates import (
    EXPERIMENTAL_HC,
    METHODS,
    build_templates,
    read_csv,
    warp,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCE = RESULTS / "npj2022_xtrg_entropy_reconstructed.csv"
CURVES = RESULTS / "m4_finite_temperature_cluster_curves.csv"
MODELS = RESULTS / "m4_finite_temperature_cluster_models.csv"
N12_VERIFY = RESULTS / "m4_finite_temperature_cluster_verification.json"
N18 = RESULTS / "n18_stacked_m4_convergence_summary.json"
SUMMARY = RESULTS / "cluster_size_transfer_stress_summary.json"
SCENARIOS = RESULTS / "cluster_size_transfer_stress_scenarios.csv"
PAIRS = RESULTS / "cluster_size_transfer_stress_pair_search.csv"
FIGURE = FIGURES / "cluster_size_transfer_stress"

DECISION_FIELDS = np.round(np.arange(1.05, 1.85001, 0.0005), 10)
ERROR_FRACTION = 0.025


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_profiles():
    curves = read_csv(CURVES)
    models = read_csv(MODELS)
    verification = json.loads(N12_VERIFY.read_text(encoding="utf-8"))
    n18 = json.loads(N18.read_text(encoding="utf-8"))
    templates = build_templates(curves, models, verification, 0.09)
    nominal = {
        "model_label": "M3-nominal",
        "residual_peak_shift_after_bulk_alignment_T": 0.0,
        "heat_scale_to_M3": 1.0,
        "left_width_scale_to_M3": 1.0,
        "right_width_scale_to_M3": 1.0,
        "template_origin": "nominal",
    }
    for template in templates:
        template["template_origin"] = "N12"

    strongest = next(
        row for row in models if row["model_label"] == "M4-01")
    comp = n18["comparison"]
    n18_template = {
        "model_label": "M4-01-N18-size-stress",
        "residual_peak_shift_after_bulk_alignment_T": (
            comp["N18_quadratic_M4_minus_M3_peak_field_T"]
            - float(strongest["critical_field_shift_vs_M3_T"])
        ),
        "heat_scale_to_M3": comp["N18_M4_to_M3_peak_heat_ratio"],
        "left_width_scale_to_M3": comp["N18_M4_to_M3_W95_left_ratio"],
        "right_width_scale_to_M3": comp["N18_M4_to_M3_W95_right_ratio"],
        "template_origin": "N18_size_stress",
    }
    all_templates = [nominal, *templates, n18_template]

    source = np.genfromtxt(
        SOURCE, delimiter=",", names=True, encoding="utf-8-sig")
    profile = np.sort(
        source[(source["anchor_N"] == 12)
               & np.isclose(source["temperature_K"], 0.09)],
        order="field_T")
    fields = profile["field_T"]
    heat = profile["Q_proxy_mJ_per_cm3"]
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)

    records = []
    for anchor in EXPERIMENTAL_HC:
        for method in METHODS:
            for template in all_templates:
                warped_field, warped_heat = warp(
                    fields, heat, source_hc, anchor, method, template)
                qmax = float(np.max(warped_heat))
                peak = float(warped_field[int(np.argmax(warped_heat))])
                decision_heat = np.interp(
                    DECISION_FIELDS, warped_field, warped_heat,
                    left=np.nan, right=np.nan)
                records.append({
                    "scenario_id": f"Hc={anchor:.2f}|{method}|{template['model_label']}",
                    "anchor_T": anchor,
                    "mapping": method,
                    "template": template["model_label"],
                    "template_origin": template["template_origin"],
                    "peak_field_T": peak,
                    "Qmax_mJ_per_cm3": qmax,
                    "residual_peak_shift_T": template["residual_peak_shift_after_bulk_alignment_T"],
                    "heat_scale": template["heat_scale_to_M3"],
                    "left_width_scale": template["left_width_scale_to_M3"],
                    "right_width_scale": template["right_width_scale_to_M3"],
                    "warped_field": warped_field,
                    "warped_heat": warped_heat,
                    "decision_heat": decision_heat,
                })
    return records, n18_template


def robust_baseline(regret):
    worst = np.nanmax(regret, axis=0)
    valid = np.all(np.isfinite(regret), axis=0)
    index = np.where(valid)[0][int(np.argmin(worst[valid]))]
    return {
        "field_T": float(DECISION_FIELDS[index]),
        "worst_regret_percent": float(100 * worst[index]),
    }


def main():
    FIGURES.mkdir(exist_ok=True)
    records, n18_template = build_profiles()
    qmax = np.array([row["Qmax_mJ_per_cm3"] for row in records])
    qdecision = np.array([row["decision_heat"] for row in records])
    regret = (qmax[:, None] - qdecision) / qmax[:, None]
    augmented = robust_baseline(regret)

    base_indices = np.array([
        index for index, row in enumerate(records)
        if row["template_origin"] != "N18_size_stress"
    ])
    base = robust_baseline(regret[base_indices])

    peak_support = [
        min(row["peak_field_T"] for row in records),
        max(row["peak_field_T"] for row in records),
    ]
    start = np.floor(peak_support[0] / 0.004) * 0.004
    stop = np.ceil(peak_support[1] / 0.004) * 0.004
    measurement_fields = np.round(np.arange(start, stop + 0.002, 0.004), 10)
    values = np.array([
        [np.interp(field, row["warped_field"], row["warped_heat"], left=np.nan, right=np.nan)
         for row in records]
        for field in measurement_fields
    ])
    if not np.all(np.isfinite(values)):
        raise RuntimeError("expanded measurement grid is outside at least one warped profile")
    q_scales = np.full(len(measurement_fields), float(np.median(qmax)))

    best_single = None
    for index, field in enumerate(measurement_fields):
        result = evaluate_panel(values, q_scales, regret, [index], ERROR_FRACTION)
        metric = (result["worst_regret_fraction"], result["maximum_remaining_scenarios"], index)
        if best_single is None or metric < best_single[0]:
            best_single = (metric, index, result)

    cache = {}
    best_pair = None
    pair_rows = []
    for first, second in itertools.combinations(range(len(measurement_fields)), 2):
        result = evaluate_panel(
            values, q_scales, regret, [first, second], ERROR_FRACTION, cache)
        metric = (
            result["worst_regret_fraction"],
            result["maximum_remaining_scenarios"],
            result["mean_remaining_scenarios"],
            first, second,
        )
        if best_pair is None or metric < best_pair[0]:
            best_pair = (metric, first, second, result)
        pair_rows.append({
            "first_field_T": measurement_fields[first],
            "second_field_T": measurement_fields[second],
            "worst_regret_percent": 100 * result["worst_regret_fraction"],
            "maximum_remaining_scenarios": result["maximum_remaining_scenarios"],
            "mean_remaining_scenarios": result["mean_remaining_scenarios"],
        })

    full = evaluate_panel(
        values, q_scales, regret, list(range(len(measurement_fields))),
        ERROR_FRACTION)
    former_pair_indices = [
        int(np.where(np.isclose(measurement_fields, value))[0][0])
        for value in (1.296, 1.420)
    ]
    former_pair = evaluate_panel(
        values, q_scales, regret, former_pair_indices, ERROR_FRACTION)
    _, single_index, single_result = best_single
    _, pair_first, pair_second, pair_result = best_pair

    scenario_rows = [{key: value for key, value in row.items()
                      if key not in ("warped_field", "warped_heat", "decision_heat")}
                     for row in records]
    payload = {
        "status": "FAIL: the previous N=12 robust field does not transfer; its two-point panel remains sufficient at the declared error but is no longer optimal",
        "scope": "adversarial 90 mK cluster-size sensitivity, not a bulk prediction or confidence interval",
        "scenario_count": len(records),
        "scenario_definition": "52 N=12-calibrated scenarios plus four correlated N=18 size-stress scenarios",
        "N18_template": n18_template,
        "peak_support_T": peak_support,
        "baseline_N12_only": base,
        "augmented_N12_plus_N18": augmented,
        "measurement_grid": {
            "start_T": float(measurement_fields[0]),
            "stop_T": float(measurement_fields[-1]),
            "spacing_T": 0.004,
            "count": len(measurement_fields),
        },
        "per_reading_error_fraction_of_median_Qmax": ERROR_FRACTION,
        "best_single": {
            "field_T": float(measurement_fields[single_index]),
            "worst_regret_percent": 100 * single_result["worst_regret_fraction"],
            "maximum_remaining_scenarios": single_result["maximum_remaining_scenarios"],
        },
        "exact_best_pair": {
            "search_space_pairs": len(pair_rows),
            "fields_T": [float(measurement_fields[pair_first]), float(measurement_fields[pair_second])],
            "worst_regret_percent": 100 * pair_result["worst_regret_fraction"],
            "maximum_remaining_scenarios": pair_result["maximum_remaining_scenarios"],
        },
        "former_N12_optimized_pair_retest": {
            "fields_T": [1.296, 1.420],
            "worst_regret_percent": 100 * former_pair["worst_regret_fraction"],
            "maximum_remaining_scenarios": former_pair["maximum_remaining_scenarios"],
            "still_below_5percent_gate": bool(
                100 * former_pair["worst_regret_fraction"] < 5.0),
        },
        "full_sweep": {
            "worst_regret_percent": 100 * full["worst_regret_fraction"],
            "maximum_remaining_scenarios": full["maximum_remaining_scenarios"],
        },
        "engineering_decision": {
            "former_1p366T_field": "rejected as a final set point; retained only as an N=12-conditional example",
            "former_1p296_1p420T_pair": "still closes the 5% regret gate at the declared 2.5% reading bound, but loses exact optimality to 1.324/1.404 T",
            "required_next_step": "measure the expanded 90 mK heat-profile field support or obtain a bulk-converged M4 thermodynamic surface",
        },
    }
    write_csv(SCENARIOS, scenario_rows)
    write_csv(PAIRS, pair_rows)
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3), constrained_layout=True)
    for row in records:
        color = "#d34a33" if row["template_origin"] == "N18_size_stress" else "#4c78a8"
        alpha = 0.9 if row["template_origin"] == "N18_size_stress" else 0.16
        lw = 1.8 if row["template_origin"] == "N18_size_stress" else 0.8
        mask = ((row["warped_field"] >= 1.20)
                & (row["warped_field"] <= 1.75))
        axes[0].plot(row["warped_field"][mask], row["warped_heat"][mask],
                     color=color, alpha=alpha, lw=lw)
    axes[0].axvline(base["field_T"], color="#4c78a8", ls="--", label="N=12 minimax")
    axes[0].axvline(augmented["field_T"], color="#d34a33", ls="--", label="size-stressed minimax")
    axes[0].set(xlabel="Applied-field coordinate (T)", ylabel=r"90 mK heat proxy (mJ cm$^{-3}$)",
                title="(a) N=18 adds a high-field size stress", xlim=(1.20, 1.75), ylim=(0.0, 0.52))
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    top = sorted(pair_rows, key=lambda row: row["worst_regret_percent"])[:80]
    scatter = axes[1].scatter(
        [row["first_field_T"] for row in top],
        [row["second_field_T"] for row in top],
        c=[row["worst_regret_percent"] for row in top], cmap="viridis_r", s=28)
    axes[1].scatter([measurement_fields[pair_first]], [measurement_fields[pair_second]],
                    marker="*", s=160, color="#d34a33", edgecolor="black")
    axes[1].set(xlabel="First field (T)", ylabel="Second field (T)",
                title="(b) Exact two-point search after size stress")
    fig.colorbar(scatter, ax=axes[1], label="Worst regret (%)")
    axes[1].grid(alpha=0.2)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(f"{FIGURE}.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)
    print(json.dumps({
        "base": base,
        "augmented": augmented,
        "best_pair": payload["exact_best_pair"],
        "full_sweep": payload["full_sweep"],
    }, indent=2))


if __name__ == "__main__":
    main()
