"""Regression checks for the P0--P2 and 90--300 mK engineering decision."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def main():
    inter = load("interlayer_engineering_scenarios.json")
    exp = load("experimental_hc_engineering_anchor.json")
    direct = load("engineering_case_refinement.json")
    with (RESULTS / "robust_conditional_design_map.csv").open(
            encoding="utf-8-sig") as f:
        design = list(csv.DictReader(f))

    checks = {}
    shift = inter["P2_minus_P0_Hc_shift_T_range"]
    checks["p2_shift_stable_across_boundaries"] = (
        abs(shift[1] - shift[0]) < 1e-12 and 0.09 < shift[0] < 0.10)
    db = inter["90mK_P2_minus_P0_delta_B_star_T_range"]
    regret = inter["90mK_P2_regret_at_P0_field_percent_range"]
    checks["p2_predeclared_gate_pass"] = (
        db[0] >= 0.05 and regret[0] >= 5.0 and
        inter["threshold_screen"][
            "P2_all_90mK_cases_above_delta_B_threshold"] and
        inter["threshold_screen"][
            "P2_all_90mK_cases_above_regret_threshold"])
    checks["p1_does_not_overclaim"] = not inter["threshold_screen"][
        "P1_all_90mK_cases_above_delta_B_threshold"]
    checks["robust_temperature_scope"] = (
        [row["temperature_K"] for row in inter["robust_fields"]] ==
        [0.09, 0.10, 0.12, 0.15] and
        all(row["worst_case_Q_mJ_per_cm3"] > 0
            for row in inter["robust_fields"]))
    checks["conditional_design_grid_complete"] = (
        len(design) == 4 * 3 * 2 and
        {float(row["assumed_total_cold_load_uW"]) for row in design} ==
        {0.05, 0.10, 0.50} and
        {float(row["required_hold_minutes"]) for row in design} == {30, 60})
    checks["design_energy_identity"] = all(
        abs(float(row["required_heat_absorption_mJ"]) -
            float(row["assumed_total_cold_load_uW"]) * 1e-3 *
            float(row["required_hold_minutes"]) * 60) < 1e-12
        for row in design)
    checks["direct_public_data_90_to_300mK"] = (
        [row["T_cold_K"] for row in direct["rankings"]] ==
        [0.09, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30])
    robust90 = next(row for row in exp["robust_fields"]
                    if row["temperature_K"] == 0.09)
    checks["primary_experimental_robust_field"] = (
        abs(robust90["robust_field_T"] - 1.41085) < 1e-8 and
        robust90["worst_case_regret_percent"] < 5.0)
    checks["scope_limits_present"] = (
        len(inter["limitations"]) >= 5 and len(direct["scope"]) >= 5)

    failed = [name for name, passed in checks.items() if not passed]
    payload = {"checks": checks, "failed": failed,
               "status": "PASS" if not failed else "FAIL"}
    (RESULTS / "engineering_decision_verification.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
