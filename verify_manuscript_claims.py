"""Verify that frozen manuscript numbers match authoritative result files."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
MANUSCRIPT = ROOT / "manuscript" / "Manuscript_draft_20261012.md"
OUTPUT = RESULTS / "manuscript_claim_verification.json"


def close(value: float, target: float, atol: float = 5e-6) -> bool:
    return abs(float(value) - float(target)) <= atol


def main() -> None:
    bulk = json.loads((RESULTS / "bulk_polarized_stability.json").read_text(encoding="utf-8"))
    inter = json.loads((RESULTS / "interlayer_engineering_scenarios.json").read_text(encoding="utf-8"))
    exp = json.loads((RESULTS / "experimental_hc_engineering_anchor.json").read_text(encoding="utf-8"))
    eng = json.loads((RESULTS / "engineering_case_refinement.json").read_text(encoding="utf-8"))
    state = json.loads((RESULTS / "spectral_flow_stateprep_N9_field_scan.json").read_text(encoding="utf-8"))
    qpu = json.loads((RESULTS / "open_plan_qpu_pilot_N6_summary.json").read_text(encoding="utf-8"))
    text = MANUSCRIPT.read_text(encoding="utf-8")

    rows = {(r["model"], r["boundary"]): r for r in bulk["rows"]}
    rankings = {int(round(r["T_cold_K"] * 1000)): r for r in eng["rankings"]}
    checks = {
        "manuscript_exists": MANUSCRIPT.exists(),
        "m0_field": close(rows[("M0", "tin_foil")]["critical_field_T"], 1.7743602),
        "m1_field": close(rows[("M1", "tin_foil")]["critical_field_T"], 1.7751632),
        "m2_field": close(rows[("M2", "tin_foil")]["critical_field_T"], 1.8568211),
        "m3_tin_field": close(rows[("M3", "tin_foil")]["critical_field_T"], 1.7723974),
        "m3_sphere_field": close(rows[("M3", "sphere_vacuum")]["critical_field_T"], 1.8260773),
        "m4_tin_field": close(rows[("M4", "tin_foil")]["critical_field_T"], 1.8676712),
        "m4_sphere_field": close(rows[("M4", "sphere_vacuum")]["critical_field_T"], 1.9213511),
        "p2_delta_gate": min(inter["90mK_P2_minus_P0_delta_B_star_T_range"]) >= 0.05,
        "p2_regret_gate": min(inter["90mK_P2_regret_at_P0_field_percent_range"]) >= 5.0,
        "experimental_robust_field": close(exp["robust_fields"][0]["robust_field_T"], 1.41085, 1e-8),
        "candidate_schedule": (
            rankings[90]["best_candidate_B_T"] == 1.5
            and rankings[120]["best_candidate_B_T"] == 1.5
            and rankings[150]["best_candidate_B_T"] == 1.3
            and rankings[200]["best_candidate_B_T"] == 1.3
            and rankings[250]["best_candidate_B_T"] == 0.0
            and rankings[300]["best_candidate_B_T"] == 0.0
        ),
        "all_stateprep_targets_pass": state["status"] == "PASS" and all(
            r["all_targets_pass"] and r["all_signs_match"] for r in state["field_summaries"]
        ),
        "two_qpu_jobs": len(qpu["observations"]) == 2,
        "qpu_ids_present": {j["job_id"] for j in qpu["observations"]}
        == {"db5r8184qg6s73c3sl0g", "db5r8nimb58s738atca0"},
        "qpu_usage_71s": close(qpu["total_free_qpu_seconds_used"], 71.0, 1e-12),
        "scope_no_quantum_advantage": "not a quantum-advantage demonstration" in text,
        "scope_qpu_quantitative_hold": "preclude a quantitative hardware claim" in text,
        "authors_and_orcids": all(
            token in text
            for token in [
                "Ugwiyeon Lee", "0009-0008-8976-5363",
                "Seong-Yun Bae", "0009-0003-8731-428X",
                "Hakrae Kim", "0009-0002-5048-3568",
                "Jungsu Lee", "0009-0006-9347-453X",
            ]
        ),
    }

    figure_links = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    checks["all_linked_figures_exist"] = bool(figure_links) and all(
        (MANUSCRIPT.parent / link).resolve().exists() for link in figure_links
    )
    checks["minimum_reference_count"] = text.count("doi:") >= 9
    failed = [name for name, passed in checks.items() if not passed]
    payload = {
        "checks": checks,
        "failed": failed,
        "linked_figures": figure_links,
        "status": "PASS" if not failed else "FAIL",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
