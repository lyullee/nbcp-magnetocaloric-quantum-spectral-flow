"""Independent arithmetic checks for the robust measurement-design audit."""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np

from optimize_measurement_value_of_information import (
    DECISION_FIELDS,
    ERROR_FRACTIONS,
    FIGURE_STEM,
    MEASUREMENT_FIELDS,
    PRECISION,
    SUMMARY,
    TEMPERATURES,
    build_scenario_profiles,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "measurement_value_of_information_verification.json"


def panel_result(values: np.ndarray, q_scales: np.ndarray,
                 regret: np.ndarray, selected: list[int],
                 error_fraction: float) -> dict:
    n = values.shape[1]
    masks = np.ones((n, n), dtype=bool)
    for index in selected:
        difference = np.abs(values[index, :, None] - values[index, None, :])
        masks &= difference <= 2.0 * error_fraction * q_scales[index] + 1e-15
    bounds = []
    fields = []
    sizes = []
    for mask in masks:
        worst = np.max(regret[mask], axis=0)
        finite = np.isfinite(worst)
        indices = np.where(finite)[0]
        best = indices[int(np.argmin(worst[finite]))]
        bounds.append(float(worst[best]))
        fields.append(float(DECISION_FIELDS[best]))
        sizes.append(int(np.sum(mask)))
    worst_true = int(np.argmax(bounds))
    return {
        "worst_regret_fraction": bounds[worst_true],
        "median_regret_fraction": float(np.median(bounds)),
        "decision_field_for_worst_true_T": fields[worst_true],
        "maximum_remaining_scenarios": max(sizes),
        "mean_remaining_scenarios": float(np.mean(sizes)),
    }


def close(a, b, atol=2e-12):
    return bool(np.isclose(float(a), float(b), atol=atol, rtol=0.0))


def main() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    keys, _, profiles = build_scenario_profiles()
    candidates = [
        (float(temperature), float(field))
        for temperature in TEMPERATURES for field in MEASUREMENT_FIELDS
    ]
    candidate_index = {coordinate: index
                       for index, coordinate in enumerate(candidates)}
    values = np.empty((len(candidates), len(keys)), dtype=float)
    q_scales = np.empty(len(candidates), dtype=float)
    for index, (temperature, field) in enumerate(candidates):
        field_index = int(np.where(np.isclose(MEASUREMENT_FIELDS, field))[0][0])
        values[index] = [
            profiles[temperature][key]["measurement_heat"][field_index]
            for key in keys]
        q_scales[index] = float(np.median([
            profiles[temperature][key]["qmax"] for key in keys]))

    q90 = np.asarray([
        profiles[0.09][key]["decision_heat"] for key in keys])
    qmax90 = np.asarray([profiles[0.09][key]["qmax"] for key in keys])
    regret = (qmax90[:, None] - q90) / qmax90[:, None]

    recomputed = {}
    sequence_matches = True
    monotone = True
    unique_points = True
    for error_fraction in ERROR_FRACTIONS:
        label = f"{100 * error_fraction:g}%"
        recomputed[label] = {}
        for strategy in ("90mK_only", "multi_temperature"):
            selected = []
            previous = float("inf")
            rows = summary["greedy_nonadaptive_panels"][label][strategy]
            for row in rows:
                if int(row["panel_size"]) > 0:
                    coordinate = (
                        float(row["selected_temperature_K"]),
                        float(row["selected_field_T"]),
                    )
                    index = candidate_index[coordinate]
                    if index in selected:
                        unique_points = False
                    selected.append(index)
                result = panel_result(
                    values, q_scales, regret, selected, error_fraction)
                for key in (
                    "worst_regret_fraction",
                    "median_regret_fraction",
                    "decision_field_for_worst_true_T",
                    "maximum_remaining_scenarios",
                    "mean_remaining_scenarios",
                ):
                    if not close(result[key], row[key], atol=5e-12):
                        sequence_matches = False
                monotone &= result["worst_regret_fraction"] <= previous + 1e-14
                previous = result["worst_regret_fraction"]
            recomputed[label][strategy] = {
                "selected_indices": selected,
                "final": panel_result(
                    values, q_scales, regret, selected, error_fraction),
            }

        full_indices = [candidate_index[(0.09, float(field))]
                        for field in MEASUREMENT_FIELDS]
        full = panel_result(
            values, q_scales, regret, full_indices, error_fraction)
        reported_full = summary["full_42_point_90mK_sweep"][label]
        recomputed[label]["full_42_point_90mK_sweep"] = full
        for key in (
            "worst_regret_fraction",
            "median_regret_fraction",
            "decision_field_for_worst_true_T",
            "maximum_remaining_scenarios",
            "mean_remaining_scenarios",
        ):
            if not close(full[key], reported_full[key], atol=5e-12):
                sequence_matches = False

    baseline = panel_result(values, q_scales, regret, [], 0.01)
    one_1pct = summary["greedy_nonadaptive_panels"]["1%"][
        "90mK_only"][1]
    two_2p5 = summary["greedy_nonadaptive_panels"]["2.5%"][
        "90mK_only"][2]
    full_5 = summary["full_42_point_90mK_sweep"]["5%"]
    multi5_2p5 = summary["greedy_nonadaptive_panels"]["2.5%"][
        "multi_temperature"][5]
    full_2p5 = summary["full_42_point_90mK_sweep"]["2.5%"]
    indices90 = [candidate_index[(0.09, float(field))]
                 for field in MEASUREMENT_FIELDS]
    exact_pair = summary["exact_90mK_two_point_panel_at_2p5pct"]
    reported_pair_indices = [candidate_index[(
        float(point["temperature_K"]), float(point["field_T"]))]
        for point in exact_pair["points"]]
    pair_recomputed = panel_result(
        values, q_scales, regret, reported_pair_indices, 0.025)
    pair_exact_match = close(
        100 * pair_recomputed["worst_regret_fraction"],
        exact_pair["worst_regret_percent"])

    # Exhaust all 861 pairs independently of the production selector.
    exhaustive_best = None
    for first, second in itertools.combinations(indices90, 2):
        result = panel_result(
            values, q_scales, regret, [first, second], 0.025)
        metric = (
            result["worst_regret_fraction"],
            result["maximum_remaining_scenarios"],
            result["mean_remaining_scenarios"],
            first, second,
        )
        if exhaustive_best is None or metric < exhaustive_best[0]:
            exhaustive_best = (metric, result)
    pair_is_global_42_choose_2 = (
        exhaustive_best is not None
        and list(exhaustive_best[0][-2:]) == reported_pair_indices
        and close(exhaustive_best[1]["worst_regret_fraction"],
                  pair_recomputed["worst_regret_fraction"])
    )

    with PRECISION.open(newline="", encoding="utf-8-sig") as handle:
        precision_rows = [
            {key: float(value) if key != "fixed_pair_fields_T" else value
             for key, value in row.items()}
            for row in csv.DictReader(handle)]
    precision = summary["precision_requirements"]
    threshold_columns = {
        "best_single_90mK_max_error_fraction_closing_5pct_regret":
            "best_single_worst_regret_percent",
        "fixed_exact_pair_max_error_fraction_closing_5pct_regret":
            "fixed_pair_worst_regret_percent",
        "full_42_point_max_error_fraction_closing_5pct_regret":
            "full_42_worst_regret_percent",
    }
    precision_thresholds_match = len(precision_rows) == 111
    for summary_key, column in threshold_columns.items():
        passing = [
            row["per_reading_error_fraction_of_median_Qmax"]
            for row in precision_rows if row[column] < 5.0]
        recovered = max(passing) if passing else None
        precision_thresholds_match &= close(recovered, precision[summary_key])
    precision_curves_monotone = all(
        np.all(np.diff([row[column] for row in precision_rows]) >= -1e-12)
        for column in threshold_columns.values())
    checks = {
        "status_pass": summary["status"] == "PASS",
        "scenario_count_52": len(keys) == summary["scenario_count"] == 52,
        "candidate_grid_4_times_42": (
            len(candidates) == 168 and len(MEASUREMENT_FIELDS) == 42),
        "all_predictions_finite": bool(np.all(np.isfinite(values))),
        "baseline_matches_primary_robust_design": (
            close(100 * baseline["worst_regret_fraction"], 9.423091560842781)
            and close(baseline["decision_field_for_worst_true_T"], 1.366)),
        "every_reported_panel_recomputed": sequence_matches,
        "panel_regret_is_monotone": bool(monotone),
        "panel_points_are_unique": unique_points,
        "one_reading_at_1pct_closes_5pct_gate":
            100 * float(one_1pct["worst_regret_fraction"]) < 5.0,
        "two_readings_at_2p5pct_close_5pct_gate":
            100 * float(two_2p5["worst_regret_fraction"]) < 5.0,
        "full_42_at_5pct_does_not_close_gate":
            float(full_5["worst_regret_percent"]) > 5.0,
        "five_point_multiT_beats_full_42_at_2p5pct": (
            100 * float(multi5_2p5["worst_regret_fraction"])
            < float(full_2p5["worst_regret_percent"])),
        "exact_pair_recomputed": pair_exact_match,
        "exact_pair_is_global_over_861_panels": pair_is_global_42_choose_2,
        "exact_pair_is_1p296_and_1p420T": [
            candidates[index][1] for index in reported_pair_indices
        ] == [1.296, 1.42],
        "precision_thresholds_match_scan": precision_thresholds_match,
        "precision_curves_are_monotone": precision_curves_monotone,
        "figures_exist": all(
            FIGURE_STEM.with_suffix(f".{suffix}").is_file()
            for suffix in ("png", "pdf", "svg")),
    }
    payload = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "failed": [key for key, value in checks.items() if not value],
        "recomputed": recomputed,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "checks": checks,
        "failed": payload["failed"],
    }, indent=2))
    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
