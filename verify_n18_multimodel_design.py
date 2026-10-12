"""Independent checks for the N=18 multimodel cooling-design propagation."""
from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main() -> None:
    data = json.loads(
        (RESULTS / "n18_multimodel_design_stress_summary.json").read_text(
            encoding="utf-8"))
    errors = []
    if data["scenario_count"] != 68:
        errors.append("expected 52 N=12 plus 16 N=18 multimodel scenarios")
    if len(data["N18_templates"]) != 4:
        errors.append("expected four correlated N=18 templates")
    models = [row["source_model"] for row in data["N18_templates"]]
    if models != ["M4-00", "M4-02", "M4-03", "M4-01"]:
        errors.append("N=18 template order differs from the preregistered cover")
    if not math.isclose(data["baseline_N12_only"]["field_T"], 1.366,
                        abs_tol=5.1e-4):
        errors.append("N=12 baseline no longer reproduces the established field")
    if not math.isclose(data["baseline_N12_only"]["worst_regret_percent"],
                        9.423091560842781, rel_tol=1e-10):
        errors.append("N=12 baseline regret changed")
    count = data["measurement_grid"]["count"]
    if data["exact_best_pair"]["search_space_pairs"] != count * (count - 1) // 2:
        errors.append("two-point search is not exhaustive")
    if (data["exact_best_pair"]["worst_regret_percent"]
            > data["former_N12_optimized_pair_retest"]["worst_regret_percent"]):
        errors.append("reported optimum is worse than the former panel")
    sizing = data["conditional_90mK_sizing"]
    expected_mass = (sizing["required_energy_mJ"]
                     / sizing["minimum_scenario_peak_heat_mJ_per_cm3"]
                     * sizing["density_g_per_cm3"])
    if not math.isclose(sizing["magnetic_material_mass_g"], expected_mass,
                        rel_tol=1e-12):
        errors.append("conditional mass identity mismatch")
    if "not a bulk prediction" not in data["scope"]:
        errors.append("scope overclaims the N=18 cover")
    report = {
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "verified": {
            "scenario_count": data["scenario_count"],
            "N18_models": models,
            "pair_count": data["exact_best_pair"]["search_space_pairs"],
            "conditional_mass_g": sizing["magnetic_material_mass_g"],
        },
    }
    (RESULTS / "n18_multimodel_design_stress_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
