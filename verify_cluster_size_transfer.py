"""Independent decision checks for the N=18 cluster-size transfer stress."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main():
    data = json.loads((RESULTS / "cluster_size_transfer_stress_summary.json").read_text(encoding="utf-8"))
    errors = []
    if data["scenario_count"] != 56:
        errors.append("expected 52 base plus four N=18 stress scenarios")
    if data["peak_support_T"][1] < 1.60:
        errors.append("N=18 stress did not expand the high-field peak support")
    if abs(data["augmented_N12_plus_N18"]["field_T"] - data["baseline_N12_only"]["field_T"]) < 0.02:
        errors.append("reported final-setpoint rejection lacks a material field change")
    if data["exact_best_pair"]["search_space_pairs"] != (
            data["measurement_grid"]["count"] * (data["measurement_grid"]["count"] - 1) // 2):
        errors.append("pair search is not exhaustive")
    if not data["former_N12_optimized_pair_retest"]["still_below_5percent_gate"]:
        errors.append("former pair unexpectedly fails the declared 5% gate")
    if data["exact_best_pair"]["worst_regret_percent"] > data["former_N12_optimized_pair_retest"]["worst_regret_percent"]:
        errors.append("reported exact optimum is worse than the former pair")
    if "not a bulk prediction" not in data["scope"]:
        errors.append("scope overclaims the N=18 stress")
    report = {"status": "PASS" if not errors else "FAIL", "errors": errors}
    (RESULTS / "cluster_size_transfer_stress_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
