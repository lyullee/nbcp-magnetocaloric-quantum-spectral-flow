"""Common-random-number verification of the strongest N=12 M4 deformation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from audit_m4_finite_temperature_cluster import (
    HOT_FIELD_T,
    PARAMETERS,
    RESULTS,
    add_interlayer_bonds,
    build_common_bonds,
    entropy_rows,
)
from bulk_polarized_stability import InterlayerParameters
from screen_model import VOLUME_FACTOR_MJ_PER_CM3


ROOT = Path(__file__).resolve().parent
MODELS = RESULTS / "m4_finite_temperature_cluster_models.csv"
COMPARISONS = RESULTS / "m4_finite_temperature_cluster_comparisons.csv"
OUTPUT = RESULTS / "m4_finite_temperature_cluster_verification.json"
FIELDS = np.round(np.arange(1.08, 1.3801, 0.01), 10)
SEEDS = (20261013, 20261014, 20261015)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def metrics(rows: list[dict]) -> dict:
    hot = next(float(row["entropy_kB_per_spin"]) for row in rows
               if np.isclose(row["field_T"], HOT_FIELD_T)
               and np.isclose(row["temperature_K"], 2.0))
    cold = sorted((row for row in rows
                   if np.isclose(row["temperature_K"], 0.09)
                   and row["field_T"] < 2.0), key=lambda row: row["field_T"])
    x = np.asarray([row["field_T"] for row in cold])
    y = np.asarray([
        0.09 * (row["entropy_kB_per_spin"] - hot)
        * VOLUME_FACTOR_MJ_PER_CM3 for row in cold])
    dense = np.linspace(float(x.min()), float(x.max()), 6001)
    smooth = PchipInterpolator(x, y)(dense)
    peak = int(np.argmax(smooth))
    threshold = 0.95 * smooth[peak]
    valid = smooth >= threshold
    low = peak
    high = peak
    while low > 0 and valid[low - 1]:
        low -= 1
    while high + 1 < len(dense) and valid[high + 1]:
        high += 1
    return {
        "B_star_T": float(dense[peak]),
        "Qmax_mJ_per_cm3": float(smooth[peak]),
        "W95_min_T": float(dense[low]),
        "W95_max_T": float(dense[high]),
        "W95_total_T": float(dense[high] - dense[low]),
        "W95_left_T": float(dense[peak] - dense[low]),
        "W95_right_T": float(dense[high] - dense[peak]),
    }


def main() -> None:
    models = read_csv(MODELS)
    comparisons = read_csv(COMPARISONS)
    comp90 = [row for row in comparisons
              if np.isclose(float(row["temperature_K"]), 0.09)]
    worst = min(comp90, key=lambda row: float(row["Qmax_ratio_to_M3"]))
    worst_label = worst["model_label"]
    worst_meta = next(row for row in models if row["model_label"] == worst_label)

    graph, common = build_common_bonds()
    parameter_sets = {
        "M3": InterlayerParameters(0.0, 0.0, 0.0, 0.08),
        worst_label: InterlayerParameters(**{
            key: float(worst_meta[key]) for key in PARAMETERS}),
    }
    trials = []
    for seed in SEEDS:
        per_model = {}
        for label, parameters in parameter_sets.items():
            bonds = add_interlayer_bonds(graph, common, parameters)
            rows, elapsed = entropy_rows(
                graph, bonds, np.r_[FIELDS, HOT_FIELD_T], seed,
                random_vectors=8, krylov_dim=60, exact_low_states=12,
                field_seed_stride=0)
            per_model[label] = {**metrics(rows), "elapsed_seconds": elapsed}
            print(json.dumps({"seed": seed, "model": label,
                              **per_model[label]}), flush=True)
        base = per_model["M3"]
        candidate = per_model[worst_label]
        trials.append({
            "seed": seed,
            "models": per_model,
            "M4_minus_M3_B_star_T":
                candidate["B_star_T"] - base["B_star_T"],
            "M4_to_M3_Qmax_ratio":
                candidate["Qmax_mJ_per_cm3"] / base["Qmax_mJ_per_cm3"],
            "M4_to_M3_W95_total_ratio":
                candidate["W95_total_T"] / base["W95_total_T"],
            "M4_to_M3_W95_left_ratio":
                candidate["W95_left_T"] / base["W95_left_T"],
            "M4_to_M3_W95_right_ratio":
                candidate["W95_right_T"] / base["W95_right_T"],
        })

    keys = (
        "M4_minus_M3_B_star_T",
        "M4_to_M3_Qmax_ratio",
        "M4_to_M3_W95_total_ratio",
        "M4_to_M3_W95_left_ratio",
        "M4_to_M3_W95_right_ratio",
    )
    summary = {}
    for key in keys:
        values = np.asarray([trial[key] for trial in trials])
        summary[key] = {
            "values": values.tolist(),
            "mean": float(np.mean(values)),
            "range": [float(np.min(values)), float(np.max(values))],
        }
    payload = {
        "status": "PASS" if (
            summary["M4_to_M3_Qmax_ratio"]["range"][0] >= 0.9
            and summary["M4_to_M3_Qmax_ratio"]["range"][1] <= 1.1
            and summary["M4_to_M3_W95_total_ratio"]["range"][0] > 1.2
        ) else "REVIEW",
        "purpose": "separate the strongest microscopic width deformation from field-to-field FTLM trace noise",
        "worst_primary_model": worst_label,
        "worst_primary_model_parameters": {
            key: float(worst_meta[key]) for key in PARAMETERS},
        "fields_T": FIELDS.tolist(),
        "common_random_number_seeds": list(SEEDS),
        "trials": trials,
        "summary": summary,
        "interpretation": [
            "The same random trace is reused across every field within a seed, reducing artificial point-to-point roughness.",
            "A total-width ratio above 1.2 in every repeat demonstrates that the previous 0.8-1.2 width stress was too narrow for this finite-cluster M4 sensitivity.",
            "This remains a finite-cluster sensitivity result, not a converged bulk confidence interval.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"status": payload["status"],
                      "summary": summary}, indent=2), flush=True)


if __name__ == "__main__":
    main()
