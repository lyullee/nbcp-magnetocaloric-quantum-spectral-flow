"""Robust value-of-information design for the correlated M4 heat scenarios.

The existing calculation recommends a dense 90 mK field sweep that localizes
the heat-profile peak.  This audit asks a different engineering question:
before that full sweep is run, which small non-adaptive panel of heat
measurements most reduces the worst-case operating-field regret?

No probability distribution is assigned to the 52 scenarios.  A declared
per-reading error bound is expressed as a fraction of the median scenario
peak heat at that temperature.  Two scenarios are conservatively treated as
indistinguishable when their error intervals overlap.  After an observation,
the operating field is re-optimized against every scenario that remains
indistinguishable.  The reported value is the worst such post-measurement
minimax regret across every possible true scenario.

This is a model-conditional experiment-design calculation.  It does not
validate the N=12 templates or turn the stress ensemble into a posterior.
"""

from __future__ import annotations

import csv
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from screen_model import MU_B_OVER_KB
from stress_test_microscopic_shape_templates import (
    EXPERIMENTAL_HC,
    METHODS,
    TEMPERATURES,
    build_templates,
    read_csv,
    warp,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
CURVES = RESULTS / "m4_finite_temperature_cluster_curves.csv"
MODELS = RESULTS / "m4_finite_temperature_cluster_models.csv"
VERIFICATION = RESULTS / "m4_finite_temperature_cluster_verification.json"
SOURCE = RESULTS / "npj2022_xtrg_entropy_reconstructed.csv"

SUMMARY = RESULTS / "measurement_value_of_information_summary.json"
SEQUENCE = RESULTS / "measurement_value_of_information_sequence.csv"
SCENARIOS = RESULTS / "measurement_value_of_information_scenarios.csv"
PRECISION = RESULTS / "measurement_value_of_information_precision.csv"
FIGURE_STEM = FIGURES / "measurement_value_of_information"

MEASUREMENT_FIELDS = np.round(np.arange(1.292, 1.456 + 0.002, 0.004), 10)
DECISION_FIELDS = np.round(np.arange(1.05, 1.75 + 0.00025, 0.0005), 10)
ERROR_FRACTIONS = (0.01, 0.025, 0.05, 0.10)
MAX_PANEL_SIZE = 6
PRECISION_SCAN = np.round(np.arange(0.005, 0.06001, 0.0005), 10)


@dataclass(frozen=True)
class Candidate:
    temperature_K: float
    field_T: float


def scenario_key(anchor: float, method: str, template: str) -> str:
    return f"Hc={anchor:.2f}|{method}|{template}"


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_scenario_profiles():
    curves = read_csv(CURVES)
    models = read_csv(MODELS)
    verification = json.loads(VERIFICATION.read_text(encoding="utf-8"))
    source = np.genfromtxt(
        SOURCE, delimiter=",", names=True, encoding="utf-8-sig")
    source = source[source["anchor_N"] == 12]
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)

    profiles: dict[float, dict[str, dict]] = {}
    metadata: dict[str, dict] = {}
    for temperature in TEMPERATURES:
        templates = build_templates(curves, models, verification, temperature)
        nominal = {
            "model_label": "M3-nominal",
            "residual_peak_shift_after_bulk_alignment_T": 0.0,
            "heat_scale_to_M3": 1.0,
            "left_width_scale_to_M3": 1.0,
            "right_width_scale_to_M3": 1.0,
        }
        profile = np.sort(
            source[np.isclose(source["temperature_K"], temperature)],
            order="field_T")
        fields = profile["field_T"]
        heat = profile["Q_proxy_mJ_per_cm3"]
        profiles[float(temperature)] = {}
        for anchor in EXPERIMENTAL_HC:
            for method in METHODS:
                for template in [nominal, *templates]:
                    key = scenario_key(anchor, method, template["model_label"])
                    warped_field, warped_heat = warp(
                        fields, heat, source_hc, anchor, method, template)
                    measurement_heat = np.interp(
                        MEASUREMENT_FIELDS, warped_field, warped_heat,
                        left=np.nan, right=np.nan)
                    decision_heat = np.interp(
                        DECISION_FIELDS, warped_field, warped_heat,
                        left=np.nan, right=np.nan)
                    qmax = float(np.max(warped_heat))
                    profiles[float(temperature)][key] = {
                        "measurement_heat": measurement_heat,
                        "decision_heat": decision_heat,
                        "qmax": qmax,
                        "peak_field_T": float(
                            warped_field[int(np.argmax(warped_heat))]),
                    }
                    metadata.setdefault(key, {
                        "scenario_id": key,
                        "experimental_Hc_anchor_T": anchor,
                        "mapping_method": method,
                        "template_model": template["model_label"],
                    })
    keys = sorted(metadata)
    for temperature in TEMPERATURES:
        if set(profiles[float(temperature)]) != set(keys):
            raise RuntimeError("scenario identity mismatch across temperatures")
    return keys, metadata, profiles


def robust_decision(regret: np.ndarray, mask: np.ndarray,
                    cache: dict[bytes, tuple[float, float]]) -> tuple[float, float]:
    packed = np.packbits(mask).tobytes()
    if packed not in cache:
        worst = np.max(regret[mask], axis=0)
        finite = np.isfinite(worst)
        if not np.any(finite):
            cache[packed] = (float("inf"), float("nan"))
        else:
            indices = np.where(finite)[0]
            best = indices[int(np.argmin(worst[finite]))]
            cache[packed] = (float(worst[best]), float(DECISION_FIELDS[best]))
    return cache[packed]


def evaluate_masks(regret: np.ndarray, masks: np.ndarray,
                   cache: dict[bytes, tuple[float, float]]) -> dict:
    bounds = []
    fields = []
    sizes = []
    for mask in masks:
        bound, field = robust_decision(regret, mask, cache)
        bounds.append(bound)
        fields.append(field)
        sizes.append(int(np.sum(mask)))
    worst_index = int(np.argmax(bounds))
    return {
        "worst_regret_fraction": float(bounds[worst_index]),
        "median_regret_fraction": float(np.median(bounds)),
        "worst_true_scenario_index": worst_index,
        "decision_field_for_worst_true_T": float(fields[worst_index]),
        "maximum_remaining_scenarios": max(sizes),
        "mean_remaining_scenarios": float(np.mean(sizes)),
    }


def greedy_panel(values: np.ndarray, regret: np.ndarray,
                 candidates: list[Candidate], q_scales: np.ndarray,
                 error_fraction: float, allowed: np.ndarray) -> list[dict]:
    scenario_count = values.shape[1]
    # If each reading has a +/- error bound, prediction intervals overlap
    # when their central values differ by at most twice that bound.
    threshold = 2.0 * error_fraction * q_scales
    compatibility = np.empty(
        (len(candidates), scenario_count, scenario_count), dtype=bool)
    for index in range(len(candidates)):
        delta = np.abs(values[index, :, None] - values[index, None, :])
        compatibility[index] = delta <= threshold[index] + 1.0e-15

    current = np.ones((scenario_count, scenario_count), dtype=bool)
    selected: list[int] = []
    rows: list[dict] = []
    cache: dict[bytes, tuple[float, float]] = {}
    baseline = evaluate_masks(regret, current, cache)
    rows.append({
        "panel_size": 0,
        "selected_temperature_K": None,
        "selected_field_T": None,
        **baseline,
    })
    for step in range(1, MAX_PANEL_SIZE + 1):
        best_index = None
        best_metric = None
        best_result = None
        for index in np.where(allowed)[0]:
            if index in selected:
                continue
            trial = current & compatibility[index]
            result = evaluate_masks(regret, trial, cache)
            metric = (
                result["worst_regret_fraction"],
                result["maximum_remaining_scenarios"],
                result["mean_remaining_scenarios"],
                float(candidates[index].temperature_K),
                float(candidates[index].field_T),
            )
            if best_metric is None or metric < best_metric:
                best_metric = metric
                best_index = int(index)
                best_result = result
        if best_index is None or best_result is None:
            raise RuntimeError("no measurement candidate remains")
        selected.append(best_index)
        current &= compatibility[best_index]
        candidate = candidates[best_index]
        rows.append({
            "panel_size": step,
            "selected_temperature_K": candidate.temperature_K,
            "selected_field_T": candidate.field_T,
            **best_result,
        })
    return rows


def compatibility_masks(values: np.ndarray, q_scales: np.ndarray,
                        selected: list[int], error_fraction: float) -> np.ndarray:
    scenario_count = values.shape[1]
    masks = np.ones((scenario_count, scenario_count), dtype=bool)
    for index in selected:
        difference = np.abs(
            values[index, :, None] - values[index, None, :])
        masks &= difference <= (
            2.0 * error_fraction * q_scales[index] + 1.0e-15)
    return masks


def evaluate_panel(values: np.ndarray, q_scales: np.ndarray,
                   regret: np.ndarray, selected: list[int],
                   error_fraction: float,
                   cache: dict[bytes, tuple[float, float]] | None = None) -> dict:
    return evaluate_masks(
        regret,
        compatibility_masks(values, q_scales, selected, error_fraction),
        {} if cache is None else cache,
    )


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    keys, metadata, profiles = build_scenario_profiles()
    scenario_count = len(keys)
    candidates = [
        Candidate(float(temperature), float(field))
        for temperature in TEMPERATURES for field in MEASUREMENT_FIELDS
    ]
    values = np.empty((len(candidates), scenario_count), dtype=float)
    q_scales = np.empty(len(candidates), dtype=float)
    for index, candidate in enumerate(candidates):
        tprofiles = profiles[candidate.temperature_K]
        field_index = int(np.where(
            np.isclose(MEASUREMENT_FIELDS, candidate.field_T))[0][0])
        values[index] = [
            tprofiles[key]["measurement_heat"][field_index] for key in keys]
        q_scales[index] = float(np.median([
            tprofiles[key]["qmax"] for key in keys]))
    if not np.all(np.isfinite(values)):
        raise RuntimeError("measurement grid extends outside a scenario profile")

    p90 = profiles[0.09]
    q90 = np.asarray([p90[key]["decision_heat"] for key in keys])
    qmax90 = np.asarray([p90[key]["qmax"] for key in keys])
    regret = (qmax90[:, None] - q90) / qmax90[:, None]

    scenario_rows = []
    for index, key in enumerate(keys):
        scenario_rows.append({
            **metadata[key],
            "peak_field_90mK_T": p90[key]["peak_field_T"],
            "Qmax_90mK_mJ_per_cm3": p90[key]["qmax"],
        })

    sequence_rows = []
    strategies: dict[str, dict[str, list[dict]]] = {}
    temperatures = np.array([candidate.temperature_K for candidate in candidates])
    allowed_sets = {
        "90mK_only": np.isclose(temperatures, 0.09),
        "multi_temperature": np.ones(len(candidates), dtype=bool),
    }
    for error_fraction in ERROR_FRACTIONS:
        label = f"{100 * error_fraction:g}%"
        strategies[label] = {}
        for strategy, allowed in allowed_sets.items():
            rows = greedy_panel(
                values, regret, candidates, q_scales,
                error_fraction, allowed)
            strategies[label][strategy] = rows
            selected_so_far = []
            for row in rows:
                if row["panel_size"]:
                    selected_so_far.append({
                        "temperature_K": row["selected_temperature_K"],
                        "field_T": row["selected_field_T"],
                    })
                sequence_rows.append({
                    "per_reading_error_fraction_of_median_Qmax": error_fraction,
                    "strategy": strategy,
                    **row,
                    "worst_regret_percent":
                        100.0 * row["worst_regret_fraction"],
                    "median_regret_percent":
                        100.0 * row["median_regret_fraction"],
                    "selected_panel_json": json.dumps(selected_so_far),
                })

    # Evaluate the complete 42-point 90 mK sweep using the same conservative
    # interval-overlap rule for each declared measurement precision.
    full_sweep = {}
    indices90 = np.where(allowed_sets["90mK_only"])[0]
    for error_fraction in ERROR_FRACTIONS:
        masks = np.ones((scenario_count, scenario_count), dtype=bool)
        for index in indices90:
            threshold = 2.0 * error_fraction * q_scales[index]
            delta = np.abs(values[index, :, None] - values[index, None, :])
            masks &= delta <= threshold + 1.0e-15
        result = evaluate_masks(regret, masks, {})
        full_sweep[f"{100 * error_fraction:g}%"] = {
            **result,
            "worst_regret_percent": 100.0 * result["worst_regret_fraction"],
            "median_regret_percent": 100.0 * result["median_regret_fraction"],
        }

    # A full enumeration of the 42 choose 2 pairs closes the only important
    # optimality gap in the short-panel recommendation.  Larger panels remain
    # explicitly labelled greedy.
    pair_error_fraction = 0.025
    pair_cache: dict[bytes, tuple[float, float]] = {}
    best_pair = None
    for first, second in itertools.combinations(indices90, 2):
        result = evaluate_panel(
            values, q_scales, regret, [int(first), int(second)],
            pair_error_fraction, pair_cache)
        metric = (
            result["worst_regret_fraction"],
            result["maximum_remaining_scenarios"],
            result["mean_remaining_scenarios"],
            int(first), int(second),
        )
        if best_pair is None or metric < best_pair[0]:
            best_pair = (metric, result)
    if best_pair is None:
        raise RuntimeError("two-point enumeration produced no result")
    exact_pair_indices = [best_pair[0][-2], best_pair[0][-1]]
    exact_pair = {
        "search_space": "all 42 choose 2 = 861 nonadaptive 90 mK panels",
        "per_reading_error_fraction_of_median_Qmax": pair_error_fraction,
        "points": [
            {
                "temperature_K": candidates[index].temperature_K,
                "field_T": candidates[index].field_T,
            }
            for index in exact_pair_indices
        ],
        **best_pair[1],
        "worst_regret_percent":
            100.0 * best_pair[1]["worst_regret_fraction"],
    }

    # Resolution requirement: scan declared per-reading bounds for the best
    # single 90 mK point, the exact 2.5%-optimized pair, and all 42 points.
    # The fixed pair curve answers how tolerant the recommended pair itself is;
    # it does not claim global pair optimality away from 2.5%.
    precision_rows = []
    for error_fraction in PRECISION_SCAN:
        single_results = []
        for index in indices90:
            result = evaluate_panel(
                values, q_scales, regret, [int(index)],
                float(error_fraction))
            single_results.append((
                result["worst_regret_fraction"], int(index), result))
        _, best_single_index, best_single_result = min(
            single_results,
            key=lambda item: (
                item[0], item[2]["maximum_remaining_scenarios"],
                item[2]["mean_remaining_scenarios"], item[1]))
        pair_result = evaluate_panel(
            values, q_scales, regret, exact_pair_indices,
            float(error_fraction))
        full_result = evaluate_panel(
            values, q_scales, regret, [int(index) for index in indices90],
            float(error_fraction))
        precision_rows.append({
            "per_reading_error_fraction_of_median_Qmax":
                float(error_fraction),
            "best_single_field_T": candidates[best_single_index].field_T,
            "best_single_worst_regret_percent":
                100.0 * best_single_result["worst_regret_fraction"],
            "fixed_pair_fields_T": json.dumps([
                candidates[index].field_T for index in exact_pair_indices]),
            "fixed_pair_worst_regret_percent":
                100.0 * pair_result["worst_regret_fraction"],
            "full_42_worst_regret_percent":
                100.0 * full_result["worst_regret_fraction"],
        })

    def maximum_closing_error(column: str) -> float | None:
        passing = [
            row["per_reading_error_fraction_of_median_Qmax"]
            for row in precision_rows if row[column] < 5.0]
        return max(passing) if passing else None

    precision_requirements = {
        "scan_increment_fraction": 0.0005,
        "best_single_90mK_max_error_fraction_closing_5pct_regret":
            maximum_closing_error("best_single_worst_regret_percent"),
        "fixed_exact_pair_max_error_fraction_closing_5pct_regret":
            maximum_closing_error("fixed_pair_worst_regret_percent"),
        "full_42_point_max_error_fraction_closing_5pct_regret":
            maximum_closing_error("full_42_worst_regret_percent"),
        "interpretation": "grid-resolved sufficient bounds within the declared 0.05 percentage-point scan; values are not experimental confidence limits",
    }

    baseline = strategies["1%"]["90mK_only"][0]
    summary = {
        "status": "PASS",
        "scope": "model-conditional robust value-of-information; no scenario probabilities",
        "scenario_count": scenario_count,
        "scenario_definition": "13 correlated templates including nominal x 2 measured Hc anchors x 2 transport mappings",
        "target_decision": "90 mK operating field",
        "baseline": {
            "worst_regret_percent":
                100.0 * baseline["worst_regret_fraction"],
            "robust_field_T": baseline["decision_field_for_worst_true_T"],
        },
        "measurement_candidates": {
            "temperatures_K": list(TEMPERATURES),
            "field_start_T": float(MEASUREMENT_FIELDS[0]),
            "field_stop_T": float(MEASUREMENT_FIELDS[-1]),
            "field_spacing_T": 0.004,
            "fields_per_temperature": len(MEASUREMENT_FIELDS),
            "candidate_count": len(candidates),
        },
        "error_model": {
            "fractions_of_temperature_specific_median_scenario_Qmax":
                list(ERROR_FRACTIONS),
            "interpretation": "each reading has a symmetric absolute bound; scenarios remain indistinguishable when their prediction intervals overlap",
            "not_interpreted_as": [
                "experimental confidence interval",
                "posterior probability",
                "validation of the N=12 templates",
            ],
        },
        "greedy_nonadaptive_panels": strategies,
        "exact_90mK_two_point_panel_at_2p5pct": exact_pair,
        "full_42_point_90mK_sweep": full_sweep,
        "precision_requirements": precision_requirements,
        "limitations": [
            "The panel is optimized only within the finite-cluster-calibrated 52-scenario set.",
            "Greedy selection is deterministic but is not a proof of the globally optimal multi-point panel.",
            "The declared reading bounds are resolution scenarios because public NBCP heat-profile uncertainty is unavailable.",
            "A real experiment must also control field calibration, thermal equilibration, addenda subtraction, and entropy integration.",
        ],
    }
    write_csv(SEQUENCE, sequence_rows)
    write_csv(SCENARIOS, scenario_rows)
    write_csv(PRECISION, precision_rows)
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.2), constrained_layout=True)
    colors = {"1%": "#4c78a8", "2.5%": "#54a24b",
              "5%": "#f58518", "10%": "#e45756"}
    for label, color in colors.items():
        for strategy, linestyle in (("90mK_only", "--"),
                                    ("multi_temperature", "-")):
            rows = strategies[label][strategy]
            axes[0].plot(
                [row["panel_size"] for row in rows],
                [100 * row["worst_regret_fraction"] for row in rows],
                marker="o", ms=3.5, color=color, ls=linestyle,
                label=f"{label}, {'multi-T' if strategy == 'multi_temperature' else '90 mK'}")
    axes[0].axhline(5.0, color="black", lw=0.8, ls=":")
    axes[0].set(
        xlabel="Measurements in greedy panel",
        ylabel="Worst post-measurement regret (%)",
        title="Robust value of heat measurements",
        xticks=range(MAX_PANEL_SIZE + 1),
    )
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=7, ncol=2)

    axes[1].plot(
        [100 * row["per_reading_error_fraction_of_median_Qmax"]
         for row in precision_rows],
        [row["best_single_worst_regret_percent"] for row in precision_rows],
        color="#4c78a8", label="best single 90 mK point")
    axes[1].plot(
        [100 * row["per_reading_error_fraction_of_median_Qmax"]
         for row in precision_rows],
        [row["fixed_pair_worst_regret_percent"] for row in precision_rows],
        color="#54a24b", label="1.296/1.420 T pair")
    axes[1].plot(
        [100 * row["per_reading_error_fraction_of_median_Qmax"]
         for row in precision_rows],
        [row["full_42_worst_regret_percent"] for row in precision_rows],
        color="#f58518", label="full 42-point sweep")
    axes[1].axhline(5.0, color="black", lw=0.8, ls=":")
    axes[1].set(
        xlabel=r"Per-reading bound (% of median $Q_{max}$)",
        ylabel="Worst post-measurement regret (%)",
        title="Required heat resolution",
    )
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False, fontsize=7)

    example = strategies["2.5%"]["multi_temperature"][1:]
    axes[2].scatter(
        [row["selected_field_T"] for row in example],
        [1000 * row["selected_temperature_K"] for row in example],
        c=range(1, len(example) + 1), cmap="viridis", s=75,
        edgecolor="black", linewidth=0.4)
    for row in example:
        axes[2].annotate(
            str(row["panel_size"]),
            (row["selected_field_T"], 1000 * row["selected_temperature_K"]),
            ha="center", va="center", color="white", fontsize=7)
    axes[2].set(
        xlabel="Field (T)", ylabel="Temperature (mK)",
        title="2.5% robust greedy panel (selection order)",
        xlim=(MEASUREMENT_FIELDS[0] - 0.01, MEASUREMENT_FIELDS[-1] + 0.01),
        yticks=[90, 100, 120, 150],
    )
    axes[2].grid(alpha=0.2)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURE_STEM.with_suffix(f".{suffix}"),
                    dpi=300 if suffix == "png" else None)
    plt.close(fig)
    print(json.dumps({
        "summary": str(SUMMARY.relative_to(ROOT)),
        "baseline_worst_regret_percent":
            summary["baseline"]["worst_regret_percent"],
        "multi_temperature_2p5pct": [
            {
                "panel_size": row["panel_size"],
                "temperature_K": row["selected_temperature_K"],
                "field_T": row["selected_field_T"],
                "worst_regret_percent": 100 * row["worst_regret_fraction"],
                "max_remaining": row["maximum_remaining_scenarios"],
            }
            for row in strategies["2.5%"]["multi_temperature"]
        ],
        "full_sweep_2p5pct": full_sweep["2.5%"],
        "exact_two_point_2p5pct": exact_pair,
        "precision_requirements": precision_requirements,
    }, indent=2))


if __name__ == "__main__":
    main()
