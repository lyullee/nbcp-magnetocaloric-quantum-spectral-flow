"""Independent invariants for the N=18 stacked M4 size-escalation audit."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main():
    data = json.loads((RESULTS / "n18_stacked_m4_convergence_summary.json").read_text(encoding="utf-8"))
    errors = []
    comp = data["comparison"]
    m3 = data["models"]["M3"]
    m4 = data["models"]["M4-01"]

    expected_shift = m4["quadratic_peak_field_T"] - m3["quadratic_peak_field_T"]
    if not math.isclose(comp["N18_quadratic_M4_minus_M3_peak_field_T"], expected_shift, abs_tol=1e-12):
        errors.append("quadratic peak shift mismatch")
    expected_ratio = m4["converged_peak_heat_mJ_per_cm3"] / m3["converged_peak_heat_mJ_per_cm3"]
    if not math.isclose(comp["N18_M4_to_M3_peak_heat_ratio"], expected_ratio, rel_tol=1e-12):
        errors.append("heat ratio mismatch")
    if expected_shift < 0.35:
        errors.append("N=18 peak shift no longer shows the observed size escalation")
    if expected_shift / comp["N12_common_random_M4_minus_M3_peak_field_T"] < 3.0:
        errors.append("N=18/N=12 shift amplification is below the major-effect threshold")
    if data["numerical_checks"]["cold_cutoff_crosscheck_max_abs_kB_per_spin"] > 1e-3:
        errors.append("cold cutoff convergence exceeds 1e-3 kB/spin")
    if data["numerical_checks"]["hot_last_cutoff_change_max_abs_kB_per_spin"] > 5e-3:
        errors.append("hot cutoff convergence exceeds 5e-3 kB/spin")
    if data["numerical_checks"]["odd_xz_yz_coupling_max_to_all_ratio"] > 1e-12:
        errors.append("parity-breaking Ewald terms are not negligible")
    if data["engineering_decision"]["finite_size_transfer_gate"] != "FAIL":
        errors.append("finite-size transfer must fail for these results")
    if "not thermodynamic-limit" not in data["scope"]:
        errors.append("scope does not prevent a bulk-convergence overclaim")

    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "verified": {
            "N18_peak_shift_T": expected_shift,
            "N18_peak_heat_ratio": expected_ratio,
            "N18_over_N12_shift_amplification": expected_shift / comp["N12_common_random_M4_minus_M3_peak_field_T"],
            "cold_cutoff_max_abs_kB_per_spin": data["numerical_checks"]["cold_cutoff_crosscheck_max_abs_kB_per_spin"],
            "hot_last_cutoff_change_max_abs_kB_per_spin": data["numerical_checks"]["hot_last_cutoff_change_max_abs_kB_per_spin"],
        },
    }
    (RESULTS / "n18_stacked_m4_convergence_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
