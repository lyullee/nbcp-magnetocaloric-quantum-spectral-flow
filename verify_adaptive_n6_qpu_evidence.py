"""Verify frozen adaptive N=6 live-QPU evidence and manuscript claims."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def main():
    local = load("adaptive_compressed_stateprep_N6.json")
    preflight = load("adaptive_n6_live_preflight.json")
    trex = load("adaptive_n6_qpu_complete_TWIRLED_TREX.json")
    zne = load("adaptive_n6_qpu_complete_TWIRLED_TREX_ZNE.json")
    summary = load("adaptive_n6_qpu_complete_summary.json")
    checks = {}
    checks["local_complete_and_accurate"] = (
        local["passes_entropy_target"] and
        local["minimum_fidelity"] >= 0.9999 and
        local["maximum_absolute_V_error_K"] <= 0.001 and
        len([row for node in local["nodes"] for row in node["states"]]) == 9)
    checks["preflight_open_only_go"] = (
        preflight["account"]["plans_preference"] == ["open"] and
        not preflight["account"]["credential_or_instance_identifier_recorded"] and
        preflight["decision"] == "GO_COMPLETE_TWIRLED_TREX_PILOT")
    checks["paired_jobs_complete"] = (
        trex["status"] == zne["status"] == "DONE" and
        trex["job_id"] != zne["job_id"] and
        trex["backend"] == zne["backend"] == "ibm_kingston" and
        len(trex["observations"]) == len(zne["observations"]) == 9)
    checks["same_problem_and_precision"] = (
        trex["condition"] == zne["condition"] and
        trex["target_precision_K"] == zne["target_precision_K"] == 0.02 and
        trex["mitigation"]["ZNE"] is False and zne["mitigation"]["ZNE"] is True)
    checks["open_plan_guards_and_usage"] = (
        trex["open_plan_only_guard_passed"] and
        zne["open_plan_only_guard_passed"] and
        not trex["usage_after"]["usage_limit_reached"] and
        not zne["usage_after"]["usage_limit_reached"] and
        trex["usage_after"]["usage_consumed_seconds"] -
        trex["usage_before"]["usage_consumed_seconds"] == 38 and
        zne["usage_after"]["usage_consumed_seconds"] -
        zne["usage_before"]["usage_consumed_seconds"] == 113)
    ttherm = trex["thermodynamics"]
    ztherm = zne["thermodynamics"]
    checks["sign_resolved"] = (
        ttherm["hardware_sign_matches_exact"] and
        ztherm["hardware_sign_matches_exact"] and
        ttherm["reported_uncertainty_propagation"]["does_not_cover_zero"] and
        ztherm["reported_uncertainty_propagation"]["does_not_cover_zero"])
    checks["honest_magnitude_hold"] = (
        not ttherm["reported_uncertainty_propagation"]["covers_exact_full"] and
        not ztherm["reported_uncertainty_propagation"]["covers_exact_full"] and
        summary["decision"] ==
        "PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE" and
        not summary["accuracy"]["quantitative_accuracy_pass"])
    checks["zne_moves_toward_exact"] = (
        0 < summary["accuracy"]["ZNE_bias_reduction_fraction"] < 1 and
        abs(summary["accuracy"]["ZNE_bias_reduction_fraction"] -
            0.2580371729045793) < 1e-12)
    checks["figure_and_table_exist"] = (
        (ROOT / "figures" / "adaptive_n6_qpu_complete.png").exists() and
        (ROOT / "figures" / "adaptive_n6_qpu_complete.pdf").exists() and
        (RESULTS / "adaptive_n6_qpu_observations.csv").exists())
    manuscript = (ROOT / "manuscript" /
                  "Manuscript_Cryogenics_20261012.md").read_text(
                      encoding="utf-8")
    checks["manuscript_claims_frozen"] = all(item in manuscript for item in (
        "−0.021679", "−0.026822", "−0.041611", "25.8%",
        "db6et9clf4us73c582f0", "db6eurqmb58s738bmisg",
        "fails the quantitative-magnitude gate"))
    public_text = "\n".join((
        (RESULTS / "adaptive_n6_live_preflight.json").read_text(
            encoding="utf-8"),
        (RESULTS / "adaptive_n6_qpu_complete_TWIRLED_TREX.json").read_text(
            encoding="utf-8"),
        (RESULTS / "adaptive_n6_qpu_complete_TWIRLED_TREX_ZNE.json").read_text(
            encoding="utf-8"),
    )).lower()
    checks["no_secret_or_service_crn"] = (
        '"token"' not in public_text and "crn:v1" not in public_text and
        '"instance"' not in public_text)
    failed = [name for name, passed in checks.items() if not passed]
    payload = {"checks": checks, "failed": failed,
               "status": "PASS" if not failed else "FAIL"}
    output = RESULTS / "adaptive_n6_qpu_complete_verification.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
