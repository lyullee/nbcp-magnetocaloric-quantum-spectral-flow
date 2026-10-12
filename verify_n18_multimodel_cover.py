"""Independent invariants for the preregistered N=18 model-cover result."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main() -> None:
    data = json.loads(
        (RESULTS / "n18_multimodel_cover_summary.json").read_text(encoding="utf-8"))
    amendment = json.loads(
        (RESULTS / "n18_multimodel_convergence_protocol_amendment.json").read_text(
            encoding="utf-8"))
    candidate_gap = float(amendment["amended_candidate_gap_kB_per_spin"])
    errors = []
    expected = ["M3", "M4-00", "M4-02", "M4-03", "M4-01"]
    if list(data["models"]) != expected:
        errors.append("compiled model order differs from the preregistered cover")
    if len(data["comparisons"]) != 4:
        errors.append("the M4 comparison cover is incomplete")
    baseline = data["models"]["M3"]
    for row in data["comparisons"]:
        item = data["models"][row["model"]]
        shift = item["quadratic_peak_field_T"] - baseline["quadratic_peak_field_T"]
        ratio = (item["converged_peak_heat_mJ_per_cm3"]
                 / baseline["converged_peak_heat_mJ_per_cm3"])
        if not math.isclose(row["N18_optimum_shift_T"], shift, abs_tol=1e-12):
            errors.append(f"{row['model']} peak-shift mismatch")
        if not math.isclose(row["N18_peak_heat_ratio"], ratio, rel_tol=1e-12):
            errors.append(f"{row['model']} heat-ratio mismatch")
        expected_left = (item["quadratic_W95_left_T"]
                         / baseline["quadratic_W95_left_T"])
        expected_right = (item["quadratic_W95_right_T"]
                          / baseline["quadratic_W95_right_T"])
        expected_total = (item["quadratic_W95_total_T"]
                          / baseline["quadratic_W95_total_T"])
        if not math.isclose(row["N18_W95_left_ratio"], expected_left,
                            rel_tol=1e-12):
            errors.append(f"{row['model']} left-width ratio mismatch")
        if not math.isclose(row["N18_W95_right_ratio"], expected_right,
                            rel_tol=1e-12):
            errors.append(f"{row['model']} right-width ratio mismatch")
        if not math.isclose(row["N18_W95_total_ratio"], expected_total,
                            rel_tol=1e-12):
            errors.append(f"{row['model']} total-width ratio mismatch")
    for model in expected:
        item = data["models"][model]
        fields = item["dense_fields_T"]
        if len(fields) != 51 or not math.isclose(fields[0], 1.10, abs_tol=1e-12) \
                or not math.isclose(fields[-1], 2.10, abs_tol=1e-12):
            errors.append(f"{model} dense-grid coverage changed")
        if any(item[key] <= 0.0 for key in (
                "quadratic_W95_left_T", "quadratic_W95_right_T",
                "quadratic_W95_total_T")):
            errors.append(f"{model} has a nonpositive local W95 width")
        entropy = np.asarray(item["dense_entropy_kB_per_spin"], dtype=float)
        expected_candidates = [
            float(fields[index]) for index in range(1, len(fields) - 1)
            if entropy[index] >= entropy[index - 1]
            and entropy[index] >= entropy[index + 1]
            and entropy[index] >= float(entropy.max()) - candidate_gap
        ]
        actual_candidates = [
            float(row["grid_center_T"])
            for row in item["L32_candidate_peaks"]
        ]
        if actual_candidates != expected_candidates:
            errors.append(
                f"{model} L32 candidates {actual_candidates} do not match "
                f"dense-grid rule {expected_candidates}")
        for candidate in item["L32_candidate_peaks"]:
            center = float(candidate["grid_center_T"])
            expected_neighbors = [round(center - 0.02, 2), center,
                                  round(center + 0.02, 2)]
            if candidate["neighbor_fields_T"] != expected_neighbors:
                errors.append(f"{model} candidate {center} has wrong neighbors")
            if len(candidate["L32_entropy_kB_per_spin"]) != 3:
                errors.append(f"{model} candidate {center} lacks three L32 values")
        selected = max(
            item["L32_candidate_peaks"],
            key=lambda row: row["local_peak_entropy_kB_per_spin"])
        if not math.isclose(item["grid_peak_field_T"],
                            selected["grid_center_T"], abs_tol=1e-12):
            errors.append(f"{model} selected center is not the largest L32 candidate")
    checks = data["numerical_checks"]
    family = data["model_family_tests"]
    expected_counts = {
        model: len(data["models"][model]["L32_candidate_peaks"])
        for model in expected
    }
    expected_multiple = [
        model for model in expected if expected_counts[model] > 1
    ]
    expected_reselected = [
        model for model in expected
        if not math.isclose(
            data["models"][model]["dense_L12_peak_field_T"],
            data["models"][model]["grid_peak_field_T"], abs_tol=1e-12)
    ]
    if family["L32_candidate_peak_count_by_model"] != expected_counts:
        errors.append("model-family candidate counts do not match model records")
    if family["multiple_close_peak_models"] != expected_multiple:
        errors.append("multiple-close-peak model list is inconsistent")
    if family["L32_reselected_peak_branch_models"] != expected_reselected:
        errors.append("L32 branch-reselection list is inconsistent")
    if not checks["all_dense_peaks_interior"]:
        errors.append("at least one dense-grid peak lies on a scan boundary")
    if not checks["all_close_dense_local_maxima_checked_at_L32"]:
        errors.append("close dense-grid local maxima were not all checked at L32")
    if checks["cold_cutoff_crosscheck_max_abs_kB_per_spin"] > 1e-3:
        errors.append("cold cutoff convergence exceeds 1e-3 kB/spin")
    if checks["hot_last_cutoff_change_max_abs_kB_per_spin"] > 5e-3:
        errors.append("hot cutoff convergence exceeds 5e-3 kB/spin")
    if data["engineering_decision"]["finite_size_transfer_gate"] != "FAIL":
        errors.append("finite-size transfer gate status changed")
    if "not thermodynamic-limit" not in data["scope"]:
        errors.append("scope does not explicitly prevent a bulk overclaim")
    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "verified": {
            "models": expected,
            "major_change_model_count": family["major_change_model_count"],
            "multiple_close_peak_models": expected_multiple,
            "L32_reselected_peak_branch_models": expected_reselected,
            "cold_cutoff_max_abs_kB_per_spin": checks["cold_cutoff_crosscheck_max_abs_kB_per_spin"],
            "hot_last_cutoff_change_max_abs_kB_per_spin": checks["hot_last_cutoff_change_max_abs_kB_per_spin"],
        },
    }
    (RESULTS / "n18_multimodel_cover_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
