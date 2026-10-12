"""Propagate the preregistered N=18 M4 cover into the 90 mK design audit.

The N=18 curves remain finite-cluster stresses.  Each microscopic model keeps
its peak residual, heat ratio, and left/right width ratios correlated.  Four
experimental-anchor/mapping combinations are generated per N=18 M4 model and
are appended to the 52 established N=12 scenarios without probabilities.
"""
from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import optimize_measurement_value_of_information as measurement_voi
from optimize_measurement_value_of_information import evaluate_panel
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
N18 = RESULTS / "n18_multimodel_cover_summary.json"
OLD_STRESS = RESULTS / "cluster_size_transfer_stress_summary.json"
SUMMARY = RESULTS / "n18_multimodel_design_stress_summary.json"
SCENARIOS = RESULTS / "n18_multimodel_design_stress_scenarios.csv"
PAIRS = RESULTS / "n18_multimodel_design_pair_search.csv"
FIGURE = FIGURES / "n18_multimodel_design_stress"

ERROR_FRACTION = 0.025
DENSITY_G_PER_CM3 = 4.206
HOLD_ENERGY_MJ = 0.10e-3 * 3600.0  # 0.10 microwatt for 60 min


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def n18_templates(data: dict) -> list[dict]:
    templates = []
    for row in data["comparisons"]:
        templates.append({
            "model_label": f"{row['model']}-N18-cover",
            "source_model": row["model"],
            "sample_id": row["sample_id"],
            "residual_peak_shift_after_bulk_alignment_T": (
                row["N18_optimum_shift_T"] - row["bulk_shift_T"]),
            "heat_scale_to_M3": row["N18_peak_heat_ratio"],
            "left_width_scale_to_M3": row["N18_W95_left_ratio"],
            "right_width_scale_to_M3": row["N18_W95_right_ratio"],
            "template_origin": "N18_multimodel_size_stress",
        })
    return templates


def build_profiles() -> tuple[list[dict], list[dict]]:
    curves = read_csv(CURVES)
    models = read_csv(MODELS)
    verification = json.loads(N12_VERIFY.read_text(encoding="utf-8"))
    n18 = json.loads(N18.read_text(encoding="utf-8"))
    templates12 = build_templates(curves, models, verification, 0.09)
    for template in templates12:
        template["template_origin"] = "N12"
    templates18 = n18_templates(n18)
    nominal = {
        "model_label": "M3-nominal",
        "residual_peak_shift_after_bulk_alignment_T": 0.0,
        "heat_scale_to_M3": 1.0,
        "left_width_scale_to_M3": 1.0,
        "right_width_scale_to_M3": 1.0,
        "template_origin": "nominal",
    }

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
            for template in [nominal, *templates12, *templates18]:
                warped_field, warped_heat = warp(
                    fields, heat, source_hc, anchor, method, template)
                qmax = float(np.max(warped_heat))
                peak = float(warped_field[int(np.argmax(warped_heat))])
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
                })
    return records, templates18


def robust_baseline(regret: np.ndarray, fields: np.ndarray) -> dict:
    worst = np.nanmax(regret, axis=0)
    valid = np.all(np.isfinite(regret), axis=0)
    index = np.where(valid)[0][int(np.argmin(worst[valid]))]
    return {
        "field_T": float(fields[index]),
        "worst_regret_percent": float(100 * worst[index]),
    }


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    records, templates18 = build_profiles()
    common_low = max(float(np.min(row["warped_field"])) for row in records)
    common_high = min(float(np.max(row["warped_field"])) for row in records)
    start = np.ceil(common_low / 0.0005) * 0.0005
    stop = np.floor(common_high / 0.0005) * 0.0005
    decision_fields = np.round(np.arange(start, stop + 0.00025, 0.0005), 10)
    # ``evaluate_panel`` reports the conditional robust field through its source
    # module's decision grid.  The multimodel support can exceed the original
    # N=12 grid, so bind that reporting coordinate to the actual regret columns.
    measurement_voi.DECISION_FIELDS = decision_fields

    qmax = np.array([row["Qmax_mJ_per_cm3"] for row in records])
    qdecision = np.array([
        np.interp(decision_fields, row["warped_field"], row["warped_heat"],
                  left=np.nan, right=np.nan)
        for row in records
    ])
    regret = (qmax[:, None] - qdecision) / qmax[:, None]
    augmented = robust_baseline(regret, decision_fields)
    base_indices = np.array([
        index for index, row in enumerate(records)
        if row["template_origin"] != "N18_multimodel_size_stress"
    ])
    base = robust_baseline(regret[base_indices], decision_fields)

    peak_support = [
        min(row["peak_field_T"] for row in records),
        max(row["peak_field_T"] for row in records),
    ]
    measure_start = np.floor(peak_support[0] / 0.004) * 0.004
    measure_stop = np.ceil(peak_support[1] / 0.004) * 0.004
    measurement_fields = np.round(
        np.arange(measure_start, measure_stop + 0.002, 0.004), 10)
    values = np.array([
        [np.interp(field, row["warped_field"], row["warped_heat"],
                   left=np.nan, right=np.nan) for row in records]
        for field in measurement_fields
    ])
    if not np.all(np.isfinite(values)):
        raise RuntimeError("measurement grid lies outside a warped profile")
    q_scales = np.full(len(measurement_fields), float(np.median(qmax)))

    best_single = None
    for index in range(len(measurement_fields)):
        result = evaluate_panel(values, q_scales, regret, [index], ERROR_FRACTION)
        metric = (result["worst_regret_fraction"],
                  result["maximum_remaining_scenarios"], index)
        if best_single is None or metric < best_single[0]:
            best_single = (metric, index, result)

    cache = {}
    best_pair = None
    pair_rows = []
    for first, second in itertools.combinations(range(len(measurement_fields)), 2):
        result = evaluate_panel(
            values, q_scales, regret, [first, second], ERROR_FRACTION, cache)
        metric = (result["worst_regret_fraction"],
                  result["maximum_remaining_scenarios"],
                  result["mean_remaining_scenarios"], first, second)
        if best_pair is None or metric < best_pair[0]:
            best_pair = (metric, first, second, result)
        pair_rows.append({
            "first_field_T": float(measurement_fields[first]),
            "second_field_T": float(measurement_fields[second]),
            "worst_regret_percent": 100 * result["worst_regret_fraction"],
            "maximum_remaining_scenarios": result["maximum_remaining_scenarios"],
            "mean_remaining_scenarios": result["mean_remaining_scenarios"],
        })

    full = evaluate_panel(
        values, q_scales, regret, list(range(len(measurement_fields))),
        ERROR_FRACTION)
    former_values = (1.296, 1.420)
    former_indices = [
        int(np.where(np.isclose(measurement_fields, value))[0][0])
        for value in former_values
    ]
    former_pair = evaluate_panel(
        values, q_scales, regret, former_indices, ERROR_FRACTION)
    _, single_index, single_result = best_single
    _, pair_first, pair_second, pair_result = best_pair

    min_peak_heat = float(np.min(qmax))
    conditional_mass = HOLD_ENERGY_MJ / min_peak_heat * DENSITY_G_PER_CM3
    old = json.loads(OLD_STRESS.read_text(encoding="utf-8"))
    scenario_rows = [
        {key: value for key, value in row.items()
         if key not in ("warped_field", "warped_heat")}
        for row in records
    ]
    payload = {
        "status": "preregistered four-model N=18 cover propagated into the cooling-design stress",
        "scope": "adversarial finite-cluster size sensitivity, not a bulk prediction or confidence interval",
        "scenario_count": len(records),
        "scenario_definition": (
            "52 N=12-calibrated scenarios plus 16 correlated scenarios from "
            "four N=18 M4 models, two experimental anchors, and two mappings"
        ),
        "N18_templates": templates18,
        "peak_support_T": peak_support,
        "decision_grid": {
            "start_T": float(decision_fields[0]),
            "stop_T": float(decision_fields[-1]),
            "spacing_T": 0.0005,
            "count": len(decision_fields),
        },
        "baseline_N12_only": base,
        "augmented_N12_plus_multimodel_N18": augmented,
        "previous_single_N18_stress": {
            "field_T": old["augmented_N12_plus_N18"]["field_T"],
            "worst_regret_percent": old["augmented_N12_plus_N18"]["worst_regret_percent"],
        },
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
            "fields_T": [float(measurement_fields[pair_first]),
                         float(measurement_fields[pair_second])],
            "worst_regret_percent": 100 * pair_result["worst_regret_fraction"],
            "maximum_remaining_scenarios": pair_result["maximum_remaining_scenarios"],
        },
        "former_N12_optimized_pair_retest": {
            "fields_T": list(former_values),
            "worst_regret_percent": 100 * former_pair["worst_regret_fraction"],
            "maximum_remaining_scenarios": former_pair["maximum_remaining_scenarios"],
            "still_below_5percent_gate": bool(
                100 * former_pair["worst_regret_fraction"] < 5.0),
        },
        "full_sweep": {
            "worst_regret_percent": 100 * full["worst_regret_fraction"],
            "maximum_remaining_scenarios": full["maximum_remaining_scenarios"],
        },
        "conditional_90mK_sizing": {
            "load_uW": 0.10,
            "hold_minutes": 60.0,
            "required_energy_mJ": HOLD_ENERGY_MJ,
            "minimum_scenario_peak_heat_mJ_per_cm3": min_peak_heat,
            "density_g_per_cm3": DENSITY_G_PER_CM3,
            "magnetic_material_mass_g": conditional_mass,
            "status": "finite-cluster stress value; not a package acceptance mass",
        },
    }
    write_csv(SCENARIOS, scenario_rows)
    write_csv(PAIRS, pair_rows)
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), constrained_layout=True)
    for row in records:
        is_n18 = row["template_origin"] == "N18_multimodel_size_stress"
        mask = ((row["warped_field"] >= measurement_fields[0] - 0.02)
                & (row["warped_field"] <= measurement_fields[-1] + 0.02))
        axes[0].plot(row["warped_field"][mask], row["warped_heat"][mask],
                     color="#d34a33" if is_n18 else "#4c78a8",
                     alpha=0.48 if is_n18 else 0.10,
                     lw=1.1 if is_n18 else 0.7)
    axes[0].axvline(base["field_T"], color="#4c78a8", ls="--", label="N=12 minimax")
    axes[0].axvline(augmented["field_T"], color="#d34a33", ls="--",
                    label="multimodel size-stressed minimax")
    axes[0].set(xlabel="Applied-field coordinate (T)",
                ylabel=r"90 mK heat proxy (mJ cm$^{-3}$)",
                title="(a) Four-model N=18 design stress")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    top = sorted(pair_rows, key=lambda row: row["worst_regret_percent"])[:120]
    scatter = axes[1].scatter(
        [row["first_field_T"] for row in top],
        [row["second_field_T"] for row in top],
        c=[row["worst_regret_percent"] for row in top],
        cmap="viridis_r", s=27)
    axes[1].scatter([measurement_fields[pair_first]],
                    [measurement_fields[pair_second]], marker="*", s=160,
                    color="#d34a33", edgecolor="black")
    axes[1].set(xlabel="First field (T)", ylabel="Second field (T)",
                title="(b) Exact two-point search")
    axes[1].grid(alpha=0.2)
    fig.colorbar(scatter, ax=axes[1], label="Worst regret (%)")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(f"{FIGURE}.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)
    print(json.dumps({
        "scenarios": len(records), "base": base, "augmented": augmented,
        "best_pair": payload["exact_best_pair"],
        "conditional_mass_g": conditional_mass,
    }, indent=2))


if __name__ == "__main__":
    main()
