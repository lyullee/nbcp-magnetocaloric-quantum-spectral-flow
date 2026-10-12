"""Verify that the N=18 cover was selected from pre-N=18 evidence."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
PARAMETERS = ("j2_microev", "ja_microev", "jb_microev", "delta_c")


def load_json(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    plan = load_json("n18_multimodel_preregistered_plan.json")
    amendment = load_json("n18_multimodel_plan_amendment_after_coarse.json")
    convergence = load_json("n18_multimodel_convergence_protocol.json")
    convergence_amendment = load_json(
        "n18_multimodel_convergence_protocol_amendment.json")
    n12_verify = load_json("m4_finite_temperature_cluster_verification.json")
    with (RESULTS / "interlayer_adaptive_inverse_uq_retained.csv").open(
            newline="", encoding="utf-8-sig") as handle:
        retained = list(csv.DictReader(handle))
    with (RESULTS / "m4_finite_temperature_cluster_comparisons.csv").open(
            newline="", encoding="utf-8-sig") as handle:
        n12 = list(csv.DictReader(handle))

    errors = []
    retained_by_id = {row["sample_id"]: row for row in retained}
    shifts = np.array(
        [float(row["critical_field_shift_vs_M3_T"]) for row in retained])
    expected_ids = {
        "M4-00": "761", "M4-02": "188", "M4-03": "738", "M4-01": "773"
    }
    if list(plan["models"]) != ["M3", "M4-00", "M4-02", "M4-03", "M4-01"]:
        errors.append("plan model order is not the declared ordered cover")

    for model, sample_id in expected_ids.items():
        item = plan["models"][model]
        if item["sample_id"] != sample_id:
            errors.append(f"{model} sample id mismatch")
            continue
        source = retained_by_id.get(sample_id)
        if source is None:
            errors.append(f"{model} sample is absent from the retained set")
            continue
        if not close(item["bulk_shift_T"], source["critical_field_shift_vs_M3_T"]):
            errors.append(f"{model} bulk shift mismatch")
        for key in PARAMETERS:
            if not close(item["parameters"][key], source[key]):
                errors.append(f"{model} {key} mismatch")

    if not close(plan["models"]["M4-00"]["bulk_shift_T"], shifts.min()):
        errors.append("M4-00 is not the minimum retained bulk-shift model")
    if not close(plan["models"]["M4-03"]["bulk_shift_T"], np.median(shifts)):
        errors.append("M4-03 is not the median retained bulk-shift model")
    if not close(plan["models"]["M4-01"]["bulk_shift_T"], shifts.max()):
        errors.append("M4-01 is not the maximum retained bulk-shift model")

    n12_90 = {
        row["sample_id"]: row for row in n12
        if close(row["temperature_K"], 0.09)
    }
    for model in ("M4-00", "M4-02", "M4-03"):
        item = plan["models"][model]
        source = n12_90.get(item["sample_id"])
        if source is None:
            errors.append(f"{model} has no N=12 90 mK source row")
            continue
        if not close(item["n12_shift_T"], source["cluster_M4_minus_M3_B_star_T"]):
            errors.append(f"{model} N=12 shift mismatch")
        if not close(item["n12_heat_ratio"], source["Qmax_ratio_to_M3"]):
            errors.append(f"{model} N=12 heat-ratio mismatch")

    strongest = plan["models"]["M4-01"]
    if not close(strongest["n12_shift_T"],
                 n12_verify["summary"]["M4_minus_M3_B_star_T"]["mean"]):
        errors.append("M4-01 common-random N=12 shift mismatch")
    if not close(strongest["n12_heat_ratio"],
                 n12_verify["summary"]["M4_to_M3_Qmax_ratio"]["mean"]):
        errors.append("M4-01 common-random N=12 heat-ratio mismatch")

    grid = np.asarray(amendment["amended_refine_rule"]["fields_T"], dtype=float)
    if len(grid) != 51 or not close(grid[0], 1.10) or not close(grid[-1], 2.10):
        errors.append("amended dense grid endpoints or count changed")
    if not np.allclose(np.diff(grid), 0.02, atol=1e-12, rtol=0.0):
        errors.append("amended dense grid is not uniformly spaced by 0.02 T")
    if any(len(rows) < 2 for rows in
           amendment["coarse_local_maxima_above_75_percent"].values()):
        errors.append("coarse evidence does not support the multimodality trigger")
    interpretation = plan["interpretation_rule"].lower()
    if ("thermodynamic-limit" not in interpretation
            or not any(guard in interpretation for guard in ("not", "do not"))):
        errors.append("plan does not guard against a bulk overclaim")
    if "0.01 kB/spin" not in convergence["candidate_rule"]:
        errors.append("close-peak L32 candidate rule changed")
    if "+/-0.02 T" not in convergence["selection_rule"]:
        errors.append("L32 candidate-neighbor rule changed")
    if not close(convergence_amendment["original_candidate_gap_kB_per_spin"], 0.01):
        errors.append("candidate amendment does not preserve the original gap")
    if not close(convergence_amendment["amended_candidate_gap_kB_per_spin"], 0.02):
        errors.append("candidate amendment gap is not the declared conservative value")
    observed_corrections = [
        row["absolute_change_kB_per_spin"]
        for row in convergence_amendment["observed_precompletion_evidence"]
        ["M3_existing_state_count_corrections"]
    ]
    if max(observed_corrections, default=0.0) <= 0.01:
        errors.append("state-count evidence does not justify widening the candidate gap")
    if "only adds candidates" not in convergence_amendment["bias_guard"]:
        errors.append("candidate amendment lacks its one-sided bias guard")

    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "verified": {
            "retained_model_count": len(retained),
            "ordered_cover": list(plan["models"]),
            "bulk_shift_anchors_T": {
                "minimum": float(shifts.min()),
                "median": float(np.median(shifts)),
                "maximum": float(shifts.max()),
            },
            "dense_grid": {
                "count": len(grid), "minimum_T": float(grid[0]),
                "maximum_T": float(grid[-1]), "step_T": float(grid[1] - grid[0]),
            },
            "N12_metrics_checked_before_new_N18_cover": True,
            "close_peak_L32_selection_rule_checked": True,
            "candidate_gap_amendment_kB_per_spin": 0.02,
        },
    }
    (RESULTS / "n18_multimodel_plan_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
