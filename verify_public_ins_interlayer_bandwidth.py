"""Independent arithmetic checks for the public-INS bandwidth audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RAW = ROOT / "sources" / "raw" / "data_archive_NBCPO" / "data_archive_NBCPO"


def read_csv(name: str) -> list[dict]:
    with (RESULTS / name).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def close(a: float, b: float, tolerance: float = 1e-12) -> bool:
    return abs(a-b) <= tolerance


def main() -> None:
    summary = json.loads((RESULTS / "public_ins_interlayer_bandwidth_summary.json")
                         .read_text(encoding="utf-8"))
    models = read_csv("public_ins_interlayer_bandwidth_models.csv")
    central_subset = read_csv("public_ins_interlayer_bandwidth_central_subset.csv")
    robust_subset = read_csv("public_ins_interlayer_bandwidth_robust_subset.csv")
    curves = read_csv("public_ins_interlayer_bandwidth_curves.csv")

    source = np.loadtxt(RAW / "fig17b.xy")
    linewidth_data = np.loadtxt(RAW / "fig17c.xye")
    linewidth_model = np.loadtxt(RAW / "fig17c.xy")
    threshold = float(np.max(np.ptp(source[:, 1:], axis=0))*1000.0)
    predicted = np.interp(linewidth_data[:, 0], linewidth_model[:, 0],
                          linewidth_model[:, 1])
    residual = predicted-linewidth_data[:, 1]
    chi2 = float(np.sum((residual/linewidth_data[:, 2])**2))

    bandwidth = np.asarray([
        float(row["transverse_max_domain_bandwidth_microeV"])
        for row in models
    ])
    shifts = np.asarray([
        float(row["critical_field_shift_vs_M3_T"]) for row in models
    ])
    central_ids = {
        int(float(row["sample_id"])) for row in models
        if row["central_pass"].lower() == "true"
    }
    robust_ids = {
        int(float(row["sample_id"])) for row in models
        if row["robust_pass_all_1sigma_corners"].lower() == "true"
    }
    csv_central_ids = {int(float(row["sample_id"])) for row in central_subset}
    csv_robust_ids = {int(float(row["sample_id"])) for row in robust_subset}
    rho, pvalue = spearmanr(bandwidth, shifts)
    points_by_model_domain: dict[tuple[int, int], int] = {}
    for row in curves:
        key = (int(float(row["sample_id"])), int(float(row["domain"])))
        points_by_model_domain[key] = points_by_model_domain.get(key, 0)+1

    # Independently reconstruct the calibrated unresolved-twin linewidth proxy
    # and the per-model nuisance profiles from the stored energy curves.
    wrapped_model_l = (linewidth_model[:, 0]+0.5) % 1.0-0.5
    source_sep_model = np.abs(
        np.interp(wrapped_model_l, source[:, 0], source[:, 1])
        - np.interp(wrapped_model_l, source[:, 0], source[:, 2]))
    resolution = summary["source"]["resolution_fwhm_meV"]
    equal_coefficient = float(np.sum(
        source_sep_model**2*(linewidth_model[:, 1]**2-resolution**2)
    )/np.sum(source_sep_model**4))
    source_proxy = np.sqrt(resolution**2
                           + equal_coefficient*source_sep_model**2)
    source_proxy_rmse = float(np.sqrt(np.mean(
        (source_proxy-linewidth_model[:, 1])**2)))

    curve_arrays: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    for key in points_by_model_domain:
        selected = [row for row in curves
                    if (int(float(row["sample_id"])),
                        int(float(row["domain"]))) == key]
        selected.sort(key=lambda row: float(row["l_rlu"]))
        curve_arrays[key] = (
            np.asarray([float(row["l_rlu"]) for row in selected]),
            np.asarray([float(row["energy_microeV"]) for row in selected]))
    recomputed_profile = {}
    recomputed_equal = {}
    for row in models:
        sample_id = int(float(row["sample_id"]))
        l1, e1 = curve_arrays[(sample_id, 1)]
        l2, e2 = curve_arrays[(sample_id, 2)]
        separation = np.abs(
            np.interp(linewidth_data[:, 0], l1, e1)
            - np.interp(linewidth_data[:, 0], l2, e2))/1000.0

        def score(coefficient: float) -> float:
            estimate = np.sqrt(resolution**2+coefficient*separation**2)
            return float(np.sum(
                ((estimate-linewidth_data[:, 1])/linewidth_data[:, 2])**2))

        fitted = minimize_scalar(score, bounds=(0.0, equal_coefficient),
                                 method="bounded",
                                 options={"xatol": 1e-12})
        recomputed_profile[sample_id] = float(fitted.fun)
        recomputed_equal[sample_id] = score(equal_coefficient)
    stored_profile = {
        int(float(row["sample_id"])): float(row["linewidth_profile_chi2"])
        for row in models}
    stored_equal = {
        int(float(row["sample_id"])): float(row["linewidth_equal_twin_chi2"])
        for row in models}
    profile_values = np.asarray(list(recomputed_profile.values()))
    equal_values = np.asarray(list(recomputed_equal.values()))

    central_shifts = [
        float(row["critical_field_shift_vs_M3_T"]) for row in models
        if int(float(row["sample_id"])) in central_ids
    ]
    robust_shifts = [
        float(row["critical_field_shift_vs_M3_T"]) for row in models
        if int(float(row["sample_id"])) in robust_ids
    ]
    checks = {
        "source_threshold_18p048_microeV": close(threshold, 18.048, 1e-9),
        "source_threshold_matches_summary": close(
            threshold,
            summary["source"]["declared_sensitivity_threshold_microeV"]),
        "source_linewidth_chi2_matches": close(
            chi2,
            summary["source"]["linewidth_model_vs_data"]["chi_square_no_refit"]),
        "source_linewidth_mae_matches": close(
            float(np.mean(np.abs(residual))),
            summary["source"]["linewidth_model_vs_data"]["mae_meV"]),
        "linewidth_proxy_reproduces_source_curve": (
            source_proxy_rmse < 1.5e-4
            and close(source_proxy_rmse,
                      summary["source"]["linewidth_proxy_reproduction"]
                      ["source_curve_rmse_meV"], 1e-12)
        ),
        "model_count_69": len(models) == 69,
        "curve_grid_complete": (
            len(curves) == 69*2*401
            and len(points_by_model_domain) == 69*2
            and all(count == 401 for count in points_by_model_domain.values())
        ),
        "central_pass_count_30": len(central_ids) == len(csv_central_ids) == 30,
        "central_subset_ids_match": central_ids == csv_central_ids,
        "robust_pass_count_30": len(robust_ids) == len(csv_robust_ids) == 30,
        "robust_subset_ids_match": robust_ids == csv_robust_ids,
        "central_and_robust_subsets_identical": central_ids == robust_ids,
        "central_shift_min_matches": close(
            min(central_shifts),
            summary["screen"]["central_pass_shift_range_T"][0]),
        "central_shift_max_matches": close(
            max(central_shifts),
            summary["screen"]["central_pass_shift_range_T"][1]),
        "robust_shift_min_matches": close(
            min(robust_shifts),
            summary["screen"]["robust_all_1sigma_corners_shift_range_T"][0]),
        "robust_shift_max_matches": close(
            max(robust_shifts),
            summary["screen"]["robust_all_1sigma_corners_shift_range_T"][1]),
        "bandwidth_min_matches": close(
            float(np.min(bandwidth)),
            summary["screen"]["central_bandwidth_range_microeV"][0]),
        "bandwidth_max_matches": close(
            float(np.max(bandwidth)),
            summary["screen"]["central_bandwidth_range_microeV"][1]),
        "spearman_matches": close(
            float(rho),
            summary["screen"]["spearman_bandwidth_vs_critical_field_shift"]),
        "spearman_pvalue_matches": close(
            float(pvalue),
            summary["screen"]["spearman_pvalue_descriptive_only"], 1e-34),
        "profile_chi2_recomputed_for_all_models": (
            set(stored_profile) == set(recomputed_profile)
            and all(close(stored_profile[key], recomputed_profile[key], 1e-8)
                    for key in stored_profile)
        ),
        "equal_twin_chi2_recomputed_for_all_models": (
            set(stored_equal) == set(recomputed_equal)
            and all(close(stored_equal[key], recomputed_equal[key], 1e-8)
                    for key in stored_equal)
        ),
        "all_69_profile_no_worse_than_source": (
            int(np.sum(profile_values <= chi2)) == 69
            == summary["linewidth_twin_profile_audit"]
            ["profile_models_no_worse_than_source_count"]
        ),
        "equal_twin_66_no_worse_than_source": (
            int(np.sum(equal_values <= chi2)) == 66
            == summary["linewidth_twin_profile_audit"]
            ["equal_twin_models_no_worse_than_source_count"]
        ),
        "profile_chi2_range_matches": (
            close(float(np.min(profile_values)),
                  summary["linewidth_twin_profile_audit"]
                  ["profile_chi_square_range"][0], 1e-8)
            and close(float(np.max(profile_values)),
                      summary["linewidth_twin_profile_audit"]
                      ["profile_chi_square_range"][1], 1e-8)
        ),
        "representative_axial_crosscheck_near_4p4": abs(
            summary["implementation_cross_check"]
            ["representative_M4_axial_Gamma_bandwidth_microeV"]-4.4) < 0.02,
        "microscopic_template_count_30_times_4": (
            summary["engineering_consequence"]
            ["central_sensitivity_subset_shape_stress_90mK"]
            ["scenario_count"] == 30*4
        ),
        "provenance_exists": (
            ROOT / "sources" /
            "public_ins_interlayer_bandwidth_provenance_20261012.json"
        ).is_file(),
        "figure_triplet_exists": all(
            (ROOT / "figures" /
             f"public_ins_interlayer_bandwidth_audit.{suffix}").is_file()
            for suffix in ("png", "pdf", "svg")
        ),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "recomputed": {
            "source_threshold_microeV": threshold,
            "source_linewidth_chi_square_no_refit": chi2,
            "model_count": len(models),
            "central_pass_count": len(central_ids),
            "robust_pass_count": len(robust_ids),
            "central_shift_range_T": [min(central_shifts), max(central_shifts)],
            "bandwidth_range_microeV": [float(np.min(bandwidth)),
                                          float(np.max(bandwidth))],
            "spearman_rho": float(rho),
            "linewidth_proxy_source_curve_rmse_meV": source_proxy_rmse,
            "profile_models_no_worse_than_source_count": int(np.sum(
                profile_values <= chi2)),
            "profile_chi_square_range": [float(np.min(profile_values)),
                                           float(np.max(profile_values))],
        },
    }
    (RESULTS / "public_ins_interlayer_bandwidth_verification.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
