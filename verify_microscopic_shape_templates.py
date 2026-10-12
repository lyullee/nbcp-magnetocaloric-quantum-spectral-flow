"""Regression checks for the finite-cluster-calibrated shape templates."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return abs(a - b) <= tol


def main() -> None:
    summary = json.loads(
        (RESULTS / "thermodynamic_shape_stress_summary.json").read_text(
            encoding="utf-8"))
    verification = json.loads(
        (RESULTS / "m4_finite_temperature_cluster_verification.json").read_text(
            encoding="utf-8"))
    outer = json.loads(
        (RESULTS / "thermodynamic_shape_outer_bound_summary.json").read_text(
            encoding="utf-8"))
    with (RESULTS / "thermodynamic_shape_templates.csv").open(
            encoding="utf-8-sig") as handle:
        templates = list(csv.DictReader(handle))
    with (RESULTS / "thermodynamic_shape_gate_stress.csv").open(
            encoding="utf-8-sig") as handle:
        gate = list(csv.DictReader(handle))
    with (RESULTS / "thermodynamic_shape_robust_decisions.csv").open(
            encoding="utf-8-sig") as handle:
        decisions = list(csv.DictReader(handle))

    gate90 = [row for row in gate
              if close(float(row["temperature_K"]), 0.09)]
    pass_both = [row for row in gate90
                 if float(row["delta_B_star_T"]) >= 0.05
                 and float(row["regret_percent"]) >= 5.0]
    decision90 = next(row for row in decisions
                      if close(float(row["temperature_K"]), 0.09))
    width_ratio = verification["summary"]["M4_to_M3_W95_total_ratio"]["range"]
    heat_ratio = verification["summary"]["M4_to_M3_Qmax_ratio"]["range"]

    checks = {
        "cluster_verification_passes": verification["status"] == "PASS",
        "twelve_models_at_four_temperatures": len(templates) == 48,
        "gate_grid_is_69_by_2_by_2": len(gate90) == 276,
        "correlated_gate_pass_count": len(pass_both) == 56,
        "verified_peak_shift_is_0p10T": close(
            verification["summary"]["M4_minus_M3_B_star_T"]["mean"], 0.1),
        "verified_heat_scale_remains_inside_prior_band":
            heat_ratio[0] >= 0.9 and heat_ratio[1] <= 1.1,
        "verified_width_exceeds_prior_band": width_ratio[0] > 1.2,
        "calibrated_90mK_field": close(
            float(decision90["robust_field_T"]), 1.366),
        "calibrated_90mK_regret": close(
            float(decision90["worst_case_regret_percent"]),
            9.423091560842781),
        "calibrated_90mK_heat": close(
            float(decision90["worst_case_Q_at_robust_field_mJ_per_cm3"]),
            0.382495876337409),
        "calibrated_90mK_volume": close(
            float(decision90["conditional_volume_cm3_for_0p10uW_60min"]),
            0.9411866173491376),
        "outer_bound_is_separate":
            outer["actual_material_robust_decisions"][0]["scenario_count"] == 324
            and summary["actual_material_robust_decisions"][0]["scenario_count"] == 52,
        "figures_exist": all(path.exists() and path.stat().st_size > 0 for path in (
            FIGURES / "m4_finite_temperature_cluster_audit.png",
            FIGURES / "thermodynamic_shape_stress.png",
            FIGURES / "thermodynamic_shape_outer_bound.png")),
    }
    failed = [name for name, passed in checks.items() if not passed]
    report = {
        "status": "PASS" if not failed else "FAIL",
        "checks": checks,
        "failed": failed,
        "verified_width_ratio_range": width_ratio,
        "verified_heat_ratio_range": heat_ratio,
        "90mK_correlated_gate_pass_count": len(pass_both),
    }
    (RESULTS / "microscopic_shape_template_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
