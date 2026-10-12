from pathlib import Path
import json
import re


ROOT = Path(__file__).resolve().parent
MAIN = ROOT / "manuscript" / "Manuscript_Cryogenics_20261012.md"
SI = ROOT / "manuscript" / "Supplementary_Information_Cryogenics_20261012.md"
HIGHLIGHTS = ROOT / "manuscript" / "Cryogenics_Highlights_20261012.txt"
COVER_LETTER = ROOT / "manuscript" / "Cryogenics_Cover_Letter_Draft_20261012.md"


main = MAIN.read_text(encoding="utf-8")
si = SI.read_text(encoding="utf-8")
highlights = [line.strip() for line in HIGHLIGHTS.read_text(encoding="utf-8").splitlines() if line.strip()]
cover_letter = COVER_LETTER.read_text(encoding="utf-8")
ins = json.loads((ROOT / "results" / "public_ins_interlayer_bandwidth_summary.json").read_text(encoding="utf-8"))
ins_verification = json.loads((ROOT / "results" / "public_ins_interlayer_bandwidth_verification.json").read_text(encoding="utf-8"))
measurement_voi = json.loads((ROOT / "results" / "measurement_value_of_information_summary.json").read_text(encoding="utf-8"))
measurement_voi_verification = json.loads((ROOT / "results" / "measurement_value_of_information_verification.json").read_text(encoding="utf-8"))
n18 = json.loads((ROOT / "results" / "n18_multimodel_cover_summary.json").read_text(encoding="utf-8"))
n18_verification = json.loads((ROOT / "results" / "n18_multimodel_cover_verification.json").read_text(encoding="utf-8"))
size_stress = json.loads((ROOT / "results" / "n18_multimodel_design_stress_summary.json").read_text(encoding="utf-8"))
size_stress_verification = json.loads((ROOT / "results" / "n18_multimodel_design_stress_verification.json").read_text(encoding="utf-8"))
integration = json.loads((ROOT / "results" / "n18_multimodel_manuscript_integration.json").read_text(encoding="utf-8"))
demag = json.loads((ROOT / "results" / "sample_demagnetization_envelope_summary.json").read_text(encoding="utf-8"))
demag_verification = json.loads((ROOT / "results" / "sample_demagnetization_envelope_verification.json").read_text(encoding="utf-8"))
abstract = main.split("## Abstract", 1)[1].split("## I. Introduction", 1)[0]
abstract_words = len(re.findall(r"\b[\wµ±.–-]+\b", abstract))

reference_block = main.split("## References", 1)[1]
reference_numbers = [int(value) for value in re.findall(r"(?m)^(\d+)\. ", reference_block)]
body_citations = set()
for citation in re.findall(r"\[([0-9,–-]+)\]", main.split("## References", 1)[0]):
    for part in citation.split(","):
        if "–" in part:
            first, last = [int(value) for value in part.split("–")]
            body_citations.update(range(first, last + 1))
        elif "-" in part:
            first, last = [int(value) for value in part.split("-")]
            body_citations.update(range(first, last + 1))
        else:
            body_citations.add(int(part))
linked_figures = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", main)
linked_si_figures = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", si)
comparisons = n18["comparisons"]
shift_range = (
    min(row["N18_optimum_shift_T"] for row in comparisons),
    max(row["N18_optimum_shift_T"] for row in comparisons),
)
heat_range = (
    min(row["N18_peak_heat_ratio"] for row in comparisons),
    max(row["N18_peak_heat_ratio"] for row in comparisons),
)
shift_token = f"{shift_range[0]:.3f}–{shift_range[1]:.3f}"
heat_token = f"{heat_range[0]:.3f}–{heat_range[1]:.3f}"
augmented = size_stress["augmented_N12_plus_multimodel_N18"]
grid = size_stress["measurement_grid"]
best_pair = size_stress["exact_best_pair"]
grid_token = f"{grid['start_T']:.3f}–{grid['stop_T']:.3f} T"
pair_tokens = [f"{value:.3f}" for value in best_pair["fields_T"]]

checks = {
    "main_exists": MAIN.exists(),
    "si_exists": SI.exists(),
    "cryogenics_title": "magnetic-refrigerator stage" in main.splitlines()[0],
    "target_note": "Cryogenics-targeted" in main,
    "abstract_at_most_250_words": abstract_words <= 250,
    "engineering_abstract": all(token in abstract for token in [shift_token, heat_token, f"{augmented['field_T']:.4f} T", *pair_tokens, "0.136 T", "2.72"]),
    "decision_targeted_quantum_method": "### E. Decision-targeted spectral-flow quantum-hardware benchmark" in main,
    "complete_qpu_result_section": "### I. The complete six-qubit hardware benchmark resolves sign but not magnitude" in main,
    "honest_qpu_discussion": all(token in main for token in ["passes the sign-resolution gate", "fails the quantitative-magnitude gate", "No QPU value enters"]),
    "quantum_references_present": all(token in reference_block for token in ["IBM Quantum", "Quantum finite-temperature Lanczos", "Quantum computation of finite-temperature"]),
    "references_contiguous": reference_numbers == list(range(1, len(reference_numbers) + 1)),
    "twenty_four_references": len(reference_numbers) == 24,
    "every_reference_cited": body_citations == set(reference_numbers),
    "eleven_main_figures": len(linked_figures) == 11,
    "four_si_figures": len(linked_si_figures) == 4,
    "all_main_figures_exist": all((MAIN.parent / path).resolve().exists() for path in linked_figures),
    "all_si_figures_exist": all((SI.parent / path).resolve().exists() for path in linked_si_figures),
    "si_old_quantum_sections_removed": all(token not in si for token in ["Quantum-algorithm screening", "QPU preflight and jobs", "open_plan_qpu_pilot"]),
    "si_complete_qpu_benchmark_present": all(token in si for token in ["S7.13 Complete adaptive six-qubit QPU benchmark", "db6et9clf4us73c582f0", "db6eurqmb58s738bmisg", "PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE"]),
    "qpu_evidence_synchronized": all(token in main for token in ["IBM Kingston", "−0.021679", "−0.026822", "−0.041611", "25.8%"]),
    "si_sections_renumbered": "## S8. Verification suite" in si and "## S9. Scope table for reviewers" in si,
    "field_measurement_specification": grid_token in main and f"{grid['spacing_T']:.3f} T spacing" in main and f"{grid['count']}-point" in main,
    "n12_measurement_provenance_verified": measurement_voi["status"] == "PASS" and measurement_voi_verification["status"] == "PASS",
    "n18_size_transfer_gate": (
        n18_verification["status"] == "PASS"
        and n18["engineering_decision"]["finite_size_transfer_gate"] == "FAIL"
        and len(comparisons) == 4
        and n18["model_family_tests"]["major_change_model_count"] >= 1
        and all(token in main for token in [shift_token, heat_token, "pre-convergence", "four preregistered N=18 M4"])
    ),
    "cluster_size_stressed_measurement_design": (
        size_stress_verification["status"] == "PASS"
        and size_stress["scenario_count"] == 68
        and best_pair["search_space_pairs"] == grid["count"] * (grid["count"] - 1) // 2
        and all(token in main for token in [
            f"{augmented['field_T']:.4f} T",
            f"{augmented['worst_regret_percent']:.2f}",
            f"{best_pair['search_space_pairs']:,}",
            *pair_tokens,
            f"{best_pair['worst_regret_percent']:.3f}",
        ])
    ),
    "no_stale_single_model_n18_design_claims": all(token not in main for token in [
        "1.4795 T",
        "32.32%",
        "4,753 two-point panels",
        "1.324 and 1.404 T",
        "98-field grid",
    ]),
    "multimodel_measurement_paragraph_synchronized": all(token in main for token in [
        f"All {best_pair['search_space_pairs']:,} two-point panels",
        f"{grid['count']}-field grid",
    ]),
    "multimodel_manuscript_integration": integration["status"] == "PASS",
    "sample_demagnetization_package_gate": (
        demag_verification["status"] == "PASS"
        and all(token in main for token in ["0.136 T", "2.72", "4.68", "package shape"])
    ),
    "second_material_external_validation": all(token in main for token in ["Gd2B2MoO9", "13.8 mK", "27.6 mK"]),
    "independent_device_validation": all(token in main for token in ["YbNi1.6Sn", "15.22 g", "0.254 microWatt", "0.592%"]),
    "adaptive_inverse_audit": all(token in main for token in ["823", "69", "0.0334", "0.2021", "224 of 276"]),
    "field_dependent_fan_holdout": all(token in main for token in ["four-point", "rho=0.940", "0.1073 T", "descriptive ranking"]),
    "independent_public_ins_bandwidth_and_linewidth_audit": (
        ins_verification["status"] == "PASS"
        and ins["screen"]["central_pass_count"] == 30
        and ins["screen"]["robust_all_1sigma_corners_pass_count"] == 30
        and ins["linewidth_twin_profile_audit"]["profile_models_no_worse_than_source_count"] == 69
        and ins["linewidth_twin_profile_audit"]["equal_twin_models_no_worse_than_source_count"] == 66
        and all(token in main for token in ["18.048", "Thirty of 69", "0.14757", "rho=0.882", "66 of 69", "21.36", "26.74", "all 69", "twin mixture"])
    ),
    "stage_load_scope": "sensitivity" in main and "finished-device performance" in main,
    "authors_present": all(name in main for name in ["Ugwiyeon Lee", "Seong-Yun Bae", "Hakrae Kim", "Jungsu Lee"]),
    "five_highlights": len(highlights) == 5,
    "highlights_at_most_85_characters": all(len(line) <= 85 for line in highlights),
    "cover_letter_matches_title": main.splitlines()[0].removeprefix("# ") in cover_letter,
    "cover_letter_does_not_claim_coauthor_approval": "will replace this paragraph" in cover_letter,
}

failed = [name for name, passed in checks.items() if not passed]
payload = {
    "status": "PASS" if not failed else "FAIL",
    "abstract_words": abstract_words,
    "main_words": len(re.findall(r"\b[\wµ±.–-]+\b", main)),
    "si_words": len(re.findall(r"\b[\wµ±.–-]+\b", si)),
    "reference_count": len(reference_numbers),
    "highlight_character_counts": [len(line) for line in highlights],
    "checks": checks,
    "failed": failed,
}
(ROOT / "results" / "cryogenics_submission_verification.json").write_text(
    json.dumps(payload, indent=2), encoding="utf-8"
)
print(payload)
raise SystemExit(1 if failed else 0)
