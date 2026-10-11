"""Verify state-preparation, live preflight and Open Plan QPU evidence."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def main():
    n6 = load("spectral_flow_circuit_unit_N6.json")
    n9 = load("spectral_flow_stateprep_N9_B1.45T.json")
    n9_scan = load("spectral_flow_stateprep_N9_field_scan.json")
    fake6 = load("fake_backend_noise_audit_N6.json")
    fake9 = load("fake_backend_noise_audit_N9.json")
    live = load("open_plan_preflight_N6.json")
    trex = load("open_plan_qpu_pilot_N6_TREX.json")
    zne = load("open_plan_qpu_pilot_N6_ZNE.json")
    summary = load("open_plan_qpu_pilot_N6_summary.json")

    checks = {}
    checks["n6_circuit_unit_pass"] = (
        n6["status"] == "PASS" and
        abs(n6["combined_error_kB_per_spin"]) <=
        n6["target_error_kB_per_spin"] and
        min(row["statevector_fidelity"] for row in n6["circuits"])
        > 1 - 1e-12)

    representatives = [row for node in n9["nodes"]
                       for row in node["representatives"]]
    checks["n9_stateprep_pass"] = (
        len(representatives) == 9 and
        min(row["qiskit_fidelity"] for row in representatives) >= 0.9993 and
        all(row["passes_target"] and row["sign_matches_exact"]
            for row in n9["thermodynamics"]))
    checks["n9_three_field_stateprep_pass"] = (
        n9_scan["status"] == "PASS" and
        n9_scan["fields_T"] == [1.4, 1.45, 1.5] and
        len(n9_scan["thermodynamic_rows"]) == 9 and
        all(row["passes_target"] and row["sign_matches_exact"]
            for row in n9_scan["thermodynamic_rows"]))
    checks["n9_fake_raw_stop"] = (
        fake9["decision"] == "STOP_RAW_QPU_CONTINUE_COMPRESSION_AND_MITIGATION")
    checks["n6_fake_conditional_go"] = (
        fake6["decision"] == "CONDITIONAL_MITIGATED_PILOT_CANDIDATE")

    checks["live_open_only_guard"] = (
        live["account"]["plans_preference"] == ["open"] and
        not live["account"]["credential_or_instance_identifier_recorded"] and
        live["decision"] == "GO_MINIMAL_TREX_PILOT")
    checks["live_compile_gate"] = (
        live["decision_basis"]["maximum_two_qubit_gate_count"] <= 84 and
        live["decision_basis"]["minimum_calibration_product_all_gates"] >= 0.80 and
        live["decision_basis"]["maximum_duration_over_median_T1"] <= 0.05)

    checks["two_distinct_completed_jobs"] = (
        trex["status"] == zne["status"] == "DONE" and
        trex["job_id"] != zne["job_id"] and
        trex["backend"] == zne["backend"] == "ibm_fez")
    checks["same_qpu_target"] = (
        abs(trex["exact_V_expectation_K"] -
            zne["exact_V_expectation_K"]) < 1e-12 and
        trex["condition"] == zne["condition"])
    checks["open_guards_on_jobs"] = (
        trex["open_plan_only_guard_passed"] and
        zne["open_plan_only_guard_passed"] and
        not trex["usage_after"]["usage_limit_reached"] and
        not zne["usage_after"]["usage_limit_reached"])

    trex_used = (trex["usage_after"]["usage_consumed_seconds"] -
                 trex["usage_before"]["usage_consumed_seconds"])
    zne_used = (zne["usage_after"]["usage_consumed_seconds"] -
                zne["usage_before"]["usage_consumed_seconds"])
    checks["usage_accounting"] = (
        trex_used == 19 and zne_used == 52 and
        summary["total_free_qpu_seconds_used"] == 71 and
        summary["free_qpu_seconds_remaining_after_pilot"] == 412)
    checks["honest_quantitative_hold"] = (
        not trex["within_two_reported_standard_errors"] and
        zne["within_two_reported_standard_errors"] and
        summary["decision"] ==
        "PASS_EXECUTION_PATH_HOLD_QUANTITATIVE_ACCURACY")

    public_text = "\n".join((
        (RESULTS / "open_plan_preflight_N6.json").read_text(encoding="utf-8"),
        (RESULTS / "open_plan_qpu_pilot_N6_TREX.json").read_text(
            encoding="utf-8"),
        (RESULTS / "open_plan_qpu_pilot_N6_ZNE.json").read_text(
            encoding="utf-8"),
    )).lower()
    checks["no_secret_or_crn_in_results"] = (
        '"token"' not in public_text and "crn:v1" not in public_text)

    failed = [name for name, passed in checks.items() if not passed]
    payload = {"checks": checks, "failed": failed,
               "status": "PASS" if not failed else "FAIL"}
    output = RESULTS / "quantum_stateprep_qpu_verification.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
