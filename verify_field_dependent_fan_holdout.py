"""Independent consistency checks for the field-dependent fan holdout."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def close(actual: float, expected: float, tolerance: float = 1.0e-12) -> bool:
    return abs(actual - expected) <= tolerance


def main() -> None:
    models = rows(RESULTS / "fan_curve_holdout_models.csv")
    curves = rows(RESULTS / "fan_curve_holdout_curves.csv")
    subset = rows(RESULTS / "fan_curve_holdout_no_worse_models.csv")
    summary = json.loads(
        (RESULTS / "fan_curve_holdout_summary.json").read_text(encoding="utf-8")
    )
    benchmark = float(summary["representative_model"]["holdout_mae_rlu"])

    recomputed_subset = [row for row in models
                         if float(row["holdout_mae_rlu"]) <= benchmark]
    shifts = [float(row["critical_field_shift_vs_M3_T"])
              for row in recomputed_subset]
    best = min(models, key=lambda row: float(row["holdout_mae_rlu"]))
    rho, pvalue = spearmanr(
        [float(row["holdout_mae_rlu"]) for row in models],
        [float(row["critical_field_shift_vs_M3_T"]) for row in models],
    )
    by_id: dict[int, list[dict]] = {}
    for row in curves:
        by_id.setdefault(int(row["sample_id"]), []).append(row)
    monotone = all(
        np.all(np.diff([
            float(row["delta_opt_rlu"])
            for row in sorted(group, key=lambda item: float(item["field_ratio"]))
        ]) >= -1.0e-12)
        for group in by_id.values()
    )

    checks = {
        "model_count_69": len(models) == 69,
        "curve_count_69_times_16": len(curves) == 69 * 16,
        "each_model_has_16_curve_points":
            len(by_id) == 69 and all(len(group) == 16 for group in by_id.values()),
        "all_curves_monotone": monotone,
        "near_hc_internal_max_below_0p001":
            max(float(row["near_hc_internal_difference_rlu"])
                for row in models) <= 0.001,
        "best_sample_matches_summary":
            int(best["sample_id"]) == int(summary["best_model"]["sample_id"]),
        "best_mae_matches_summary": close(
            float(best["holdout_mae_rlu"]),
            float(summary["best_model"]["holdout_mae_rlu"]),
        ),
        "relative_subset_count_matches":
            len(recomputed_subset)
            == int(summary["relative_to_published_representative"]["models_no_worse_count"])
            == len(subset),
        "relative_subset_shift_min_matches": close(
            min(shifts),
            float(summary["relative_to_published_representative"]
                  ["models_no_worse_shift_range_T"][0]),
        ),
        "relative_subset_shift_max_matches": close(
            max(shifts),
            float(summary["relative_to_published_representative"]
                  ["models_no_worse_shift_range_T"][1]),
        ),
        "spearman_matches": close(
            float(rho),
            float(summary["descriptive_association"]
                  ["spearman_holdout_mae_vs_critical_field_shift"]),
        ),
        "spearman_pvalue_matches": close(
            float(pvalue),
            float(summary["descriptive_association"]["pvalue_descriptive_only"]),
            tolerance=1.0e-38,
        ),
        "provenance_exists":
            (ROOT / "sources" / "fan_curve_holdout_provenance_20261012.json").exists(),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "recomputed": {
            "model_count": len(models),
            "curve_count": len(curves),
            "best_sample_id": int(best["sample_id"]),
            "best_mae_rlu": float(best["holdout_mae_rlu"]),
            "relative_subset_count": len(recomputed_subset),
            "relative_subset_shift_range_T": [min(shifts), max(shifts)],
            "spearman_rho": float(rho),
        },
    }
    (RESULTS / "fan_curve_holdout_verification.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
