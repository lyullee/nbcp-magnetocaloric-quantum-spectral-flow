"""Verify that frozen manuscript numbers match authoritative result files."""
from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
MANUSCRIPT = ROOT / "manuscript" / "Manuscript_draft_20261012.md"
CRYOGENICS_MANUSCRIPT = ROOT / "manuscript" / "Manuscript_Cryogenics_20261012.md"
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
    ensemble = json.loads((RESULTS / "interlayer_parameter_ensemble_summary.json").read_text(encoding="utf-8"))
    shape = json.loads((RESULTS / "thermodynamic_shape_stress_summary.json").read_text(encoding="utf-8"))
    microscopic_shape = json.loads((RESULTS / "microscopic_shape_template_verification.json").read_text(encoding="utf-8"))
    measurement_voi = json.loads((RESULTS / "measurement_value_of_information_summary.json").read_text(encoding="utf-8"))
    measurement_voi_verification = json.loads((RESULTS / "measurement_value_of_information_verification.json").read_text(encoding="utf-8"))
    measurement = json.loads((RESULTS / "peak_measurement_target_summary.json").read_text(encoding="utf-8"))
    design = json.loads((RESULTS / "refrigerator_design_chart_summary.json").read_text(encoding="utf-8"))
    boundary = json.loads((RESULTS / "interlayer_boundary_audit_summary.json").read_text(encoding="utf-8"))
    refrigerant = json.loads((RESULTS / "refrigerant_experimental_benchmark.json").read_text(encoding="utf-8"))
    refrigerant_common = json.loads((RESULTS / "refrigerant_common_basis_summary.json").read_text(encoding="utf-8"))
    refrigerant_common_verification = json.loads((RESULTS / "refrigerant_common_basis_verification.json").read_text(encoding="utf-8"))
    stage_load = json.loads((RESULTS / "adr_stage_load_envelope_summary.json").read_text(encoding="utf-8"))
    stage_load_verification = json.loads((RESULTS / "adr_stage_load_envelope_verification.json").read_text(encoding="utf-8"))
    transient_cycle = json.loads((RESULTS / "adr_transient_cycle_summary.json").read_text(encoding="utf-8"))
    transient_cycle_verification = json.loads((RESULTS / "adr_transient_cycle_verification.json").read_text(encoding="utf-8"))
    second_material = json.loads((RESULTS / "gd2b2moo9_external_validation.json").read_text(encoding="utf-8"))
    second_material_verification = json.loads((RESULTS / "gd2b2moo9_external_validation_verification.json").read_text(encoding="utf-8"))
    stage_external = json.loads((RESULTS / "ybni16sn_stage_heat_balance_summary.json").read_text(encoding="utf-8"))
    stage_external_verification = json.loads((RESULTS / "ybni16sn_stage_heat_balance_verification.json").read_text(encoding="utf-8"))
    sobol = json.loads((RESULTS / "interlayer_expanded_sobol_summary.json").read_text(encoding="utf-8"))
    sobol_shape = json.loads((RESULTS / "interlayer_expanded_sobol_shape_stress_summary.json").read_text(encoding="utf-8"))
    local_modes = json.loads((RESULTS / "interlayer_local_mode_verification.json").read_text(encoding="utf-8"))
    sobol_modes = json.loads((RESULTS / "interlayer_expanded_sobol_mode_verification.json").read_text(encoding="utf-8"))
    adaptive = json.loads((RESULTS / "interlayer_adaptive_inverse_uq_summary.json").read_text(encoding="utf-8"))
    adaptive_verification = json.loads((RESULTS / "interlayer_adaptive_inverse_uq_verification.json").read_text(encoding="utf-8"))
    fan = json.loads((RESULTS / "fan_curve_holdout_summary.json").read_text(encoding="utf-8"))
    fan_verification = json.loads((RESULTS / "fan_curve_holdout_verification.json").read_text(encoding="utf-8"))
    ins = json.loads((RESULTS / "public_ins_interlayer_bandwidth_summary.json").read_text(encoding="utf-8"))
    ins_verification = json.loads((RESULTS / "public_ins_interlayer_bandwidth_verification.json").read_text(encoding="utf-8"))
    n18 = json.loads((RESULTS / "n18_multimodel_cover_summary.json").read_text(encoding="utf-8"))
    n18_verification = json.loads((RESULTS / "n18_multimodel_cover_verification.json").read_text(encoding="utf-8"))
    size_stress = json.loads((RESULTS / "n18_multimodel_design_stress_summary.json").read_text(encoding="utf-8"))
    size_stress_verification = json.loads((RESULTS / "n18_multimodel_design_stress_verification.json").read_text(encoding="utf-8"))
    multimodel_integration = json.loads((RESULTS / "n18_multimodel_manuscript_integration.json").read_text(encoding="utf-8"))
    demag = json.loads((RESULTS / "sample_demagnetization_envelope_summary.json").read_text(encoding="utf-8"))
    demag_verification = json.loads((RESULTS / "sample_demagnetization_envelope_verification.json").read_text(encoding="utf-8"))
    text = MANUSCRIPT.read_text(encoding="utf-8")
    cryogenics_text = CRYOGENICS_MANUSCRIPT.read_text(encoding="utf-8")
    main_texts = (text, cryogenics_text)

    rows = {(r["model"], r["boundary"]): r for r in bulk["rows"]}
    rankings = {int(round(r["T_cold_K"] * 1000)): r for r in eng["rankings"]}
    n18_comparisons = n18["comparisons"]
    n18_shift_range = (
        min(row["N18_optimum_shift_T"] for row in n18_comparisons),
        max(row["N18_optimum_shift_T"] for row in n18_comparisons),
    )
    n18_heat_range = (
        min(row["N18_peak_heat_ratio"] for row in n18_comparisons),
        max(row["N18_peak_heat_ratio"] for row in n18_comparisons),
    )
    shift_token = f"{n18_shift_range[0]:.3f}–{n18_shift_range[1]:.3f}"
    heat_token = f"{n18_heat_range[0]:.3f}–{n18_heat_range[1]:.3f}"
    augmented = size_stress["augmented_N12_plus_multimodel_N18"]
    measurement_grid = size_stress["measurement_grid"]
    exact_pair = size_stress["exact_best_pair"]
    pair_tokens = [f"{value:.3f}" for value in exact_pair["fields_T"]]
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
        "ensemble_375_to_12": (
            ensemble["grid_size"] == 375
            and ensemble["retained_size"] == 12
            and ensemble["threshold_coverage"]["90mK_scenarios_above_both_gates"] == 42
            and ensemble["threshold_coverage"]["90mK_scenario_count"] == 48
        ),
        "boundary_audit_18_to_0": (
            boundary["outward_candidates"] == 18
            and boundary["outward_candidates_retained"] == 0
        ),
        "constraint_attribution": (
            boundary["constraint_attribution_counts"]["delta_only"] == 338
            and boundary["constraint_attribution_counts"]["qz_only"] == 12
            and boundary["constraint_attribution_counts"]["km_only"] == 375
            and boundary["constraint_attribution_counts"]["all_three"] == 12
        ),
        "expanded_sobol_256_to_4": (
            sobol["point_count"] == 256
            and sobol["retained_count"] == 4
            and sobol["retained_outside_original_box"] == 2
            and close(sobol["retained_ranges"]["critical_field_shift_vs_M3_T"][0], 0.07259719075801518, 1e-12)
        ),
        "expanded_central_16_of_16": (
            sobol["central_shape_threshold_coverage"]["scenario_count"] == 16
            and sobol["central_shape_threshold_coverage"]["both_gates_pass_count"] == 16
        ),
        "expanded_shape_1170_of_1296": (
            sobol_shape["scenario_count"] == 1296
            and sobol_shape["pass_count"] == 1170
            and close(sobol_shape["pass_fraction_descriptive_not_probability"], 0.9027777777777778, 1e-12)
        ),
        "local_modes_independently_verified": (
            local_modes["status"] == "PASS"
            and local_modes["retained_model_count"] == 12
            and not local_modes["models_reclassified_outside_public_q_filter"]
        ),
        "expanded_modes_independently_verified": (
            sobol_modes["status"] == "PASS"
            and sobol_modes["retained_model_count"] == 4
            and not sobol_modes["models_reclassified_outside_public_q_filter"]
        ),
        "adaptive_inverse_audit_823_to_69": (
            adaptive["design"]["initial_exact_count"] == 631
            and adaptive["design"]["adaptive_exact_count"] == 192
            and adaptive["design"]["total_exact_count"] == 823
            and adaptive["retained"]["total_count"] == 69
            and adaptive["retained"]["adaptive_new_count"] == 53
            and close(adaptive["retained"]["critical_field_shift_range_T"][0], 0.033403908987329034, 1e-12)
            and close(adaptive["retained"]["critical_field_shift_range_T"][1], 0.20210652482768032, 1e-12)
            and adaptive_verification["status"] == "PASS"
            and all(token in text for token in ["823", "69", "0.0334", "0.2021", "224 of 276"])
        ),
        "field_dependent_fan_holdout": (
            fan_verification["status"] == "PASS"
            and fan["retained_model_count"] == 69
            and close(fan["holdout_mae_range_rlu"][0], 0.0042, 1e-12)
            and close(fan["holdout_mae_range_rlu"][1], 0.00530625, 1e-12)
            and fan["relative_to_published_representative"]["models_no_worse_count"] == 23
            and close(
                fan["relative_to_published_representative"]["models_no_worse_shift_range_T"][1],
                0.10730382425956964,
                1e-12,
            )
            and fan["relative_to_published_representative"]["models_no_worse_90mK_central_both_gates"] == 40
            and fan["relative_to_published_representative"]["models_no_worse_90mK_central_case_count"] == 92
            and fan["relative_to_published_representative"]["models_no_worse_90mK_shape_both_gates"] == 0
            and fan["relative_to_published_representative"]["models_no_worse_90mK_shape_case_count"] == 92
            and close(
                fan["descriptive_association"]["spearman_holdout_mae_vs_critical_field_shift"],
                0.9403208649957064,
                1e-12,
            )
            and all(token in text for token in ["four-point", "0.00420", "0.00531", "0.940", "23 models"])
        ),
        "independent_public_ins_bandwidth_and_linewidth_audit": (
            ins_verification["status"] == "PASS"
            and ins["screen"]["input_model_count"] == 69
            and ins["screen"]["central_pass_count"] == 30
            and ins["screen"]["robust_all_1sigma_corners_pass_count"] == 30
            and close(ins["source"]["declared_sensitivity_threshold_microeV"], 18.048, 1e-9)
            and close(ins["screen"]["central_pass_shift_range_T"][1], 0.14757413094718363, 1e-12)
            and close(ins["screen"]["spearman_bandwidth_vs_critical_field_shift"], 0.8816952868103765, 1e-12)
            and ins["engineering_consequence"]["central_sensitivity_subset_90mK"]["both_gate_pass_count"] == 68
            and ins["engineering_consequence"]["central_sensitivity_subset_shape_stress_90mK"]["both_gate_pass_count"] == 4
            and ins["linewidth_twin_profile_audit"]["profile_models_no_worse_than_source_count"] == 69
            and ins["linewidth_twin_profile_audit"]["equal_twin_models_no_worse_than_source_count"] == 66
            and close(ins["linewidth_twin_profile_audit"]["profile_chi_square_range"][0], 21.36114259654802, 1e-10)
            and close(ins["linewidth_twin_profile_audit"]["profile_chi_square_range"][1], 26.74053699872553, 1e-10)
            and all(token in text for token in ["18.048", "Thirty of 69", "0.14757", "0.882", "66 of 69", "21.36", "26.74", "all 69", "twin mixture"])
        ),
        "public_time_trace_benchmark": (
            close(refrigerant["time_trace_rows_safe_for_comparison"][0]["below_1K_observed_span_min"], 95.95399, 1e-5)
            and refrigerant["caption_consistency_check"]["labels_consistent_with_article_statement"] is False
        ),
        "common_basis_refrigerant_comparison": (
            refrigerant_common["common_basis_materials"] == ["NBCP", "FAA", "CPA"]
            and refrigerant_common["common_volume_range_cm3"] == [0.395, 0.421]
            and close(
                refrigerant_common["nbcp_ratios_vs_best_hydrate_comparator"]["heat_absorption_ratio_vs_best_hydrate"],
                19.0 / 3.1,
                1e-12,
            )
            and close(
                refrigerant_common["nbcp_ratios_vs_best_hydrate_comparator"]["efficiency_ratio_vs_best_hydrate"],
                23.0 / 11.0,
                1e-12,
            )
            and close(
                refrigerant_common["nbcp_ratios_vs_best_hydrate_comparator"]["below_1K_span_ratio_vs_best_hydrate"],
                95.95399 / 32.156500000040964,
                1e-12,
            )
            and refrigerant_common_verification["status"] == "PASS"
            and all(token in text for token in ["6.13", "2.09", "2.98"])
        ),
        "steady_stage_load_envelope": (
            stage_load["rows"] == 240
            and stage_load_verification["status"] == "PASS"
            and close(stage_load["canonical_90mK_payload_0p10uW_hold_60min"][0]["conditional_NBCP_mass_g"], 3.9582202162443876, 1e-12)
            and close(stage_load["canonical_90mK_payload_0p10uW_hold_60min"][1]["conditional_NBCP_mass_g"], 11.874660648733164, 1e-12)
            and close(stage_load["canonical_90mK_payload_0p10uW_hold_60min"][2]["conditional_NBCP_mass_g"], 42.352956313814936, 1e-12)
            and all(token in text for token in ["11.87", "42.35"])
        ),
        "transient_cycle_envelope": (
            transient_cycle["rows"] == 2688
            and transient_cycle["allowance_boundary_rows"] == 35
            and transient_cycle_verification["status"] == "PASS"
            and close(
                transient_cycle["canonical_90mK_0p10uW_60min_zero_parasitic"][1]["conditional_NBCP_mass_g"],
                4.398022462493764,
                1e-12,
            )
            and close(
                transient_cycle["canonical_90mK_0p10uW_60min_zero_parasitic"][2]["conditional_NBCP_mass_g"],
                5.277626954992517,
                1e-12,
            )
            and close(
                transient_cycle["canonical_90mK_0p10uW_60min_zero_parasitic"][3]["conditional_NBCP_mass_g"],
                7.916440432488776,
                1e-12,
            )
            and close(
                transient_cycle["canonical_60min_cycle_duties"][0]["single_stage_duty_for_60min_hold"],
                3600.0 / 3682.0,
                1e-12,
            )
            and close(
                transient_cycle["canonical_60min_cycle_duties"][1]["single_stage_duty_for_60min_hold"],
                3600.0 / 3990.0,
                1e-12,
            )
            and all(token in text for token in ["4.40", "5.28", "7.92", "97.8%", "90.2%"])
        ),
        "second_material_external_validation": (
            second_material["point_count"] == 12
            and close(second_material["mean_absolute_error_K"], 0.013801210383537404, 1e-12)
            and close(second_material["maximum_absolute_error_K"], 0.02763626068590963, 1e-12)
            and second_material_verification["status"] == "PASS"
            and all(token in text for token in ["Gd2B2MoO9", "13.8 mK", "27.6 mK"])
        ),
        "independent_device_heat_balance_validation": (
            stage_external_verification["status"] == "PASS"
            and stage_external["all_validation_gates_pass"]
            and close(stage_external["fit"]["effective_mass_g"], 15.221559624354356, 1e-10)
            and close(stage_external["fit"]["parasitic_load_uW"], 0.25399304280614954, 1e-10)
            and close(stage_external["holdout_metrics"]["nrmse_vs_full_interval_enthalpy"], 0.0059175222693995345, 1e-12)
            and all(token in text for token in ["YbNi1.6Sn", "15.22 g", "0.254 microWatt", "0.592%"])
        ),
        "microscopic_shape_gate_fraction": close(
            shape["90mK_gate_stress"]["fraction_passing_both_descriptive_not_probability"],
            56 / 276,
            1e-12,
        ) and shape["90mK_gate_stress"]["scenario_count"] == 276
          and microscopic_shape["status"] == "PASS"
          and all(token in text for token in ["56", "276", "1.46", "2.05"]),
        "n12_microscopic_shape_provenance": (
            close(shape["90mK_gate_stress"]["delta_B_star_T_range"][1], 0.1, 1e-12)
            and close(shape["90mK_actual_material_peak_support_T"][0], 1.2942982936816978, 1e-12)
            and shape["90mK_no_fit_field_sweep"]["point_count"] == 42
            and all(token in text for token in ["1.366 T", "1.296", "1.420", "pre-convergence"])
        ),
        "n12_robust_field_provenance": close(
            shape["actual_material_robust_decisions"][0]["robust_field_T"],
            1.366,
            1e-12,
        ),
        "measurement_peak_target": close(
            shape["90mK_no_fit_field_sweep"]["maximum_distance_to_nearest_grid_point_T"],
            0.002, 1e-12),
        "n12_measurement_voi_provenance": (
            measurement_voi["status"] == "PASS"
            and measurement_voi_verification["status"] == "PASS"
            and measurement_voi["scenario_count"] == 52
        ),
        "n18_size_transfer_gate": (
            n18_verification["status"] == "PASS"
            and n18["engineering_decision"]["finite_size_transfer_gate"] == "FAIL"
            and len(n18_comparisons) == 4
            and n18["model_family_tests"]["major_change_model_count"] >= 1
            and n18["numerical_checks"]["all_dense_peaks_interior"]
            and n18["numerical_checks"]["all_close_dense_local_maxima_checked_at_L32"]
            and all(token in text for token in [shift_token, heat_token, "four preregistered N=18 M4", "pre-convergence"])
        ),
        "cluster_size_stressed_measurement_design": (
            size_stress_verification["status"] == "PASS"
            and size_stress["scenario_count"] == 68
            and exact_pair["search_space_pairs"] == measurement_grid["count"] * (measurement_grid["count"] - 1) // 2
            and all(token in text for token in [
                f"{measurement_grid['start_T']:.3f}–{measurement_grid['stop_T']:.3f} T",
                f"{augmented['field_T']:.4f} T",
                f"{augmented['worst_regret_percent']:.2f}",
                f"{exact_pair['search_space_pairs']:,}",
                *pair_tokens,
                f"{exact_pair['worst_regret_percent']:.3f}",
            ])
        ),
        "no_stale_single_model_n18_design_claims": all(
            token not in manuscript_text
            for manuscript_text in main_texts
            for token in (
                "1.4795 T",
                "32.32%",
                "4,753 two-point panels",
                "1.324 and 1.404 T",
                "98-field grid",
            )
        ),
        "multimodel_measurement_paragraph_synchronized": all(
            all(token in manuscript_text for token in (
                f"All {exact_pair['search_space_pairs']:,} two-point panels",
                f"{measurement_grid['count']}-field grid",
            ))
            for manuscript_text in main_texts
        ),
        "multimodel_manuscript_integration": multimodel_integration["status"] == "PASS",
        "sample_demagnetization_package_gate": (
            demag_verification["status"] == "PASS"
            and close(demag["public_B_parallel_c_specific_heat_sample"]["saturation_correction_bound_T"], 0.13643428662152593, 1e-12)
            and close(demag["shape_acceptance_requirements"]["0.020_T"]["minimum_aspect_c_over_inplane_radius"], 2.721881300941775, 1e-12)
            and close(demag["shape_acceptance_requirements"]["0.010_T"]["minimum_aspect_c_over_inplane_radius"], 4.678689436980013, 1e-12)
            and all(token in text for token in ["0.136 T", "2.72", "4.68", "package shape"])
        ),
        "design_density": close(design["density_derivation"]["density_g_cm3"], 4.205564, 1e-5),
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
        "scope_qpu_quantitative_hold": "not quantitative accuracy" in text,
        "source_label_exclusion_stated": "do not use the sheet-a/b cooling-curve labels" in text,
        "scope_q_alone_not_guarantee": "propagation vector alone cannot guarantee" in text,
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
    abstract = text.split("## Abstract", 1)[1].split("## I. Introduction", 1)[0]
    checks["abstract_at_most_250_words"] = len(re.findall(r"\b[\wµ±.–-]+\b", abstract)) <= 250
    checks["twelve_main_figures"] = len(figure_links) == 12
    checks["all_linked_figures_exist"] = bool(figure_links) and all(
        (MANUSCRIPT.parent / link).resolve().exists() for link in figure_links
    )
    supplement = ROOT / "manuscript" / "Supplementary_Information_draft_20261012.md"
    supplement_text = supplement.read_text(encoding="utf-8")
    supplement_figure_links = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", supplement_text)
    checks["all_supplementary_figures_exist"] = bool(supplement_figure_links) and all(
        (supplement.parent / link).resolve().exists() for link in supplement_figure_links
    )
    checks["minimum_reference_count"] = text.count("doi:") >= 16
    failed = [name for name, passed in checks.items() if not passed]
    payload = {
        "checks": checks,
        "failed": failed,
        "linked_figures": figure_links,
        "supplementary_linked_figures": supplement_figure_links,
        "status": "PASS" if not failed else "FAIL",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
