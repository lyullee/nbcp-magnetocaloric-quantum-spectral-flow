"""Independent INS bandwidth audit of the retained M4 interlayer models.

The diffraction-constrained M4 ensemble was selected without using the
transverse-field inelastic-neutron-scattering (INS) data published with
Phys. Rev. B 112, 104413 (2025).  That paper's public Fig. 17 data contain a
3.5 T scan along the interlayer direction and a deliberately simple Jb=3 ueV
calculation.  The authors call that value an *indicative* upper bound because
the linewidth modulation is only approximate and the calculated energy
modulation is not resolved in the data.

This audit therefore uses the 18.048 ueV branch bandwidth in the public
calculation as an external sensitivity threshold, not as a likelihood or a
confidence limit.  It asks which of the 69 diffraction-compatible M4 models
predict no more interlayer modulation on the same high-field path.  The
dominant XXZ and g-factor one-sigma errors are propagated through all eight
corners.  Every reported M4 value is an explicit LSWT calculation; no
surrogate is used.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, spearmanr

from bulk_polarized_stability import (
    EwaldSettings,
    GZ,
    InterlayerParameters,
    MICROEV_TO_K,
    S,
    coupling_fourier,
)
from screen_interlayer_parameter_ensemble import engineering_rows
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GXY, JZZ_MEV, DELTA_XY_OVER_ZZ, MEV_TO_K


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCES = ROOT / "sources"
RAW = SOURCES / "raw" / "data_archive_NBCPO" / "data_archive_NBCPO"
RETAINED = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"
SHAPE = RESULTS / "thermodynamic_shape_gate_stress.csv"

SOURCE_DISPERSION = RAW / "fig17b.xy"
SOURCE_LINEWIDTH_DATA = RAW / "fig17c.xye"
SOURCE_LINEWIDTH_MODEL = RAW / "fig17c.xy"
SOURCE_ARCHIVE = SOURCES / "raw" / "data_archive_NBCPO.zip"

SETTINGS = EwaldSettings(boundary="tin_foil")
FIELD_T = 3.5
L_GRID = np.linspace(-1.0, 1.0, 401)
# The 2025 paper uses (h,k)=(-2/3,1/3) and in-plane bonds h,k,h+k.
# The bulk code uses bonds u,v,u-v, so (u,v)=(h,-k) modulo one.
K_PATHS = ((1.0 / 3.0, 2.0 / 3.0),
           (2.0 / 3.0, 1.0 / 3.0))  # structural-twin/domain counterpart

SOURCE_CENTRAL = {"jz_meV": 0.1225, "jxy_meV": 0.0779, "gab": 4.200}
SOURCE_SIGMA = {"jz_meV": 0.0010, "jxy_meV": 0.0007, "gab": 0.007}
SOURCE_JB_MEV = 0.003
SOURCE_RESOLUTION_FWHM_MEV = 0.027
MU_B_MEV_PER_T = MU_B_OVER_KB / MEV_TO_K


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def read_retained() -> list[dict]:
    rows = []
    with RETAINED.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            row = dict(raw)
            for key in ("j2_microev", "ja_microev", "jb_microev", "delta_c",
                        "critical_field_T", "critical_field_shift_vs_M3_T",
                        "out_of_plane_bandwidth_microeV"):
                row[key] = float(row[key])
            row["sample_id"] = int(float(row["sample_id"]))
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        serial = [{key: json.dumps(value) if isinstance(value, (list, dict)) else value
                   for key, value in row.items()} for row in rows]
        writer = csv.DictWriter(handle, fieldnames=list(serial[0]))
        writer.writeheader()
        writer.writerows(serial)


def source_energy_mev(l_values: np.ndarray, jz_meV: float, jxy_meV: float,
                      gab: float, twin: bool = False) -> np.ndarray:
    """Reproduce Eqs. (A13-A14) of the 2025 source for Fig. 17(b)."""
    h, k = -2.0 / 3.0, 1.0 / 3.0
    l = -np.asarray(l_values) if twin else np.asarray(l_values)
    gamma = (np.cos(2*np.pi*h) + np.cos(2*np.pi*k)
             + np.cos(2*np.pi*(h+k)))
    gamma_b = (np.cos(2*np.pi*(h-l)) + np.cos(2*np.pi*(k-l))
               + np.cos(2*np.pi*(-h-k-l)))
    a = (gab * MU_B_MEV_PER_T * FIELD_T - 3*jxy_meV
         + 0.5*(jz_meV+jxy_meV)*gamma
         + SOURCE_JB_MEV*(-3.0+gamma_b))
    b = 0.5*(jz_meV-jxy_meV)*gamma
    return 2.0*S*np.sqrt(np.maximum(0.0, a*a-b*b))


def uncertainty_corners() -> list[dict]:
    corners = []
    for signs in itertools.product((-1, 1), repeat=3):
        row = {}
        for key, sign in zip(("jz_meV", "jxy_meV", "gab"), signs):
            row[key] = SOURCE_CENTRAL[key] + sign*SOURCE_SIGMA[key]
        row["label"] = "_".join(f"{key}{'+' if sign > 0 else '-'}1s"
                                 for key, sign in zip(
                                     ("Jz", "Jxy", "gab"), signs))
        corners.append(row)
    return corners


def precompute_m3_paths() -> dict:
    paths = {}
    for path_index, (u, v) in enumerate(K_PATHS):
        paths[path_index] = np.asarray([
            coupling_fourier(np.array([u, v, l]), "M3", SETTINGS)
            for l in L_GRID
        ])
    paths["gamma"] = np.asarray([
        coupling_fourier(np.array([0.0, 0.0, l]), "M3", SETTINGS)
        for l in L_GRID
    ])
    paths["j0"] = coupling_fourier(np.zeros(3), "M3", SETTINGS)
    return paths


def interlayer_fourier(frac_k: np.ndarray,
                       parameters: InterlayerParameters) -> np.ndarray:
    u, v, w = np.asarray(frac_k, dtype=float)
    terms = [((0, 0, 1), parameters.j2_microev)]
    terms += [(d, parameters.ja_microev) for d in
              ((1, 0, 1), (0, 1, 1), (1, -1, 1))]
    terms += [(d, parameters.jb_microev) for d in
              ((-1, 0, 1), (0, -1, 1), (-1, 1, 1))]
    result = np.zeros((3, 3), dtype=float)
    for d, value in terms:
        phase = 2*np.pi*(d[0]*u+d[1]*v+d[2]*w)
        result += 2*np.cos(phase)*np.diag([
            parameters.delta_c*value,
            parameters.delta_c*value,
            value,
        ])*MICROEV_TO_K
    return result


def intralayer_correction(frac_k: np.ndarray, jz_mev: float,
                          jxy_mev: float) -> np.ndarray:
    """Shift the bulk model's dominant XXZ terms to a source-fit corner."""
    u, v, _ = np.asarray(frac_k, dtype=float)
    gamma2 = 2.0*(np.cos(2*np.pi*u) + np.cos(2*np.pi*v)
                  + np.cos(2*np.pi*(u-v)))
    delta_jz = (jz_mev-JZZ_MEV)*MEV_TO_K
    delta_jxy = (jxy_mev-JZZ_MEV*DELTA_XY_OVER_ZZ)*MEV_TO_K
    return gamma2*np.diag([delta_jxy, delta_jxy, delta_jz])


def transverse_energy_microev(jk: np.ndarray, j0: np.ndarray,
                              gab: float) -> np.ndarray:
    """LSWT mode for polarization along y, returned in microelectronvolt."""
    h = gab*MU_B_OVER_KB*FIELD_T
    ax = h-S*j0[1, 1]+S*jk[:, 0, 0]
    az = h-S*j0[1, 1]+S*jk[:, 2, 2]
    cross = S*jk[:, 0, 2]
    energy_k = np.sqrt(np.maximum(0.0, ax*az-cross*cross))
    return energy_k/MICROEV_TO_K


def axial_energy_microev(jk: np.ndarray, j0: np.ndarray) -> np.ndarray:
    """LSWT mode for polarization along z, returned in microelectronvolt."""
    h = GZ*MU_B_OVER_KB*FIELD_T
    a = h-S*j0[2, 2]+0.5*S*(jk[:, 0, 0]+jk[:, 1, 1])
    b2 = (0.5*S*(jk[:, 0, 0]-jk[:, 1, 1]))**2 + (S*jk[:, 0, 1])**2
    return np.sqrt(np.maximum(0.0, a*a-b2))/MICROEV_TO_K


def model_curves(raw: dict, base: dict, corner: dict) -> list[np.ndarray]:
    p = InterlayerParameters(
        j2_microev=raw["j2_microev"],
        ja_microev=raw["ja_microev"],
        jb_microev=raw["jb_microev"],
        delta_c=raw["delta_c"],
    )
    j0 = (base["j0"] + interlayer_fourier(np.zeros(3), p)
          + intralayer_correction(np.zeros(3), corner["jz_meV"],
                                  corner["jxy_meV"]))
    curves = []
    for path_index, (u, v) in enumerate(K_PATHS):
        jk = base[path_index] + np.asarray([
            interlayer_fourier(np.array([u, v, l]), p)
            + intralayer_correction(np.array([u, v, l]),
                                    corner["jz_meV"], corner["jxy_meV"])
            for l in L_GRID
        ])
        curves.append(transverse_energy_microev(jk, j0, corner["gab"]))
    return curves


def source_threshold(corner: dict) -> float:
    branches = [source_energy_mev(L_GRID, corner["jz_meV"],
                                  corner["jxy_meV"], corner["gab"], twin)
                for twin in (False, True)]
    return 1000.0*max(float(np.ptp(branch)) for branch in branches)


def range_of(rows: list[dict], key: str) -> list[float]:
    values = [float(row[key]) for row in rows]
    return [min(values), max(values)]


def periodic_source_separation(l_values: np.ndarray,
                               source_dispersion: np.ndarray) -> np.ndarray:
    """Twin-branch separation with the source's unit-period l convention."""
    wrapped = (np.asarray(l_values, dtype=float)+0.5) % 1.0-0.5
    first = np.interp(wrapped, source_dispersion[:, 0],
                      source_dispersion[:, 1])
    second = np.interp(wrapped, source_dispersion[:, 0],
                       source_dispersion[:, 2])
    return np.abs(first-second)


def linewidth_from_separation(separation_mev: np.ndarray,
                              mixture_coefficient: float) -> np.ndarray:
    """Calibrated single-peak FWHM proxy for two unresolved twin branches."""
    return np.sqrt(SOURCE_RESOLUTION_FWHM_MEV**2
                   + mixture_coefficient*np.asarray(separation_mev)**2)


def minor_domain_fraction_proxy(coefficient: float,
                                equal_twin_coefficient: float) -> float:
    """Map c/c_equal=4p(1-p) to the smaller effective domain fraction."""
    ratio = float(np.clip(coefficient/equal_twin_coefficient, 0.0, 1.0))
    return 0.5*(1.0-np.sqrt(1.0-ratio))


def engineering_summary(subset: list[dict]) -> dict:
    rows = [row for row in engineering_rows(subset, 0.0)
            if np.isclose(row["temperature_K"], 0.09)]
    return {
        "scenario_count": len(rows),
        "both_gate_pass_count": sum(
            row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
            for row in rows),
        "delta_B_star_range_T": range_of(rows, "delta_B_star_T"),
        "regret_range_percent": range_of(rows, "regret_percent"),
    }


def shape_summary(subset: list[dict]) -> dict:
    selected_shifts = np.asarray([
        row["critical_field_shift_vs_M3_T"] for row in subset
    ])
    rows = []
    with SHAPE.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            if not np.isclose(float(raw["temperature_K"]), 0.09):
                continue
            shift = float(raw["critical_field_shift_vs_M3_T"])
            if np.any(np.isclose(shift, selected_shifts, rtol=0.0, atol=1e-12)):
                rows.append(raw)
    passing = sum(float(row["delta_B_star_T"]) >= 0.05
                  and float(row["regret_percent"]) >= 5.0 for row in rows)
    return {
        "scenario_count": len(rows),
        "both_gate_pass_count": passing,
        "fraction_passing_both_descriptive_not_probability": passing/len(rows),
    }


def main() -> None:
    for path in (SOURCE_ARCHIVE, SOURCE_DISPERSION, SOURCE_LINEWIDTH_DATA,
                 SOURCE_LINEWIDTH_MODEL, RETAINED, SHAPE):
        if not path.is_file():
            raise FileNotFoundError(path)

    retained = read_retained()
    if len(retained) != 69:
        raise RuntimeError(f"expected 69 retained models, found {len(retained)}")

    source_dispersion = np.loadtxt(SOURCE_DISPERSION)
    source_branch_bandwidths = np.ptp(source_dispersion[:, 1:], axis=0)*1000.0
    public_threshold = float(np.max(source_branch_bandwidths))
    source_formula = np.c_[
        source_energy_mev(source_dispersion[:, 0], **SOURCE_CENTRAL,
                          twin=False),
        source_energy_mev(source_dispersion[:, 0], **SOURCE_CENTRAL,
                          twin=True),
    ]
    direct_error = float(np.max(np.abs(source_formula-source_dispersion[:, 1:])))
    swapped_error = float(np.max(np.abs(source_formula-source_dispersion[:, [2, 1]])))
    source_formula_error = min(direct_error, swapped_error)

    linewidth_data = np.loadtxt(SOURCE_LINEWIDTH_DATA)
    linewidth_model = np.loadtxt(SOURCE_LINEWIDTH_MODEL)
    linewidth_prediction = np.interp(
        linewidth_data[:, 0], linewidth_model[:, 0], linewidth_model[:, 1])
    linewidth_residual = linewidth_prediction-linewidth_data[:, 1]
    linewidth_chi2 = float(np.sum((linewidth_residual/linewidth_data[:, 2])**2))

    # The source generated its Fig. 17(c) curve by convolving two unresolved
    # structural-twin branches with the instrumental resolution and fitting one
    # apparent peak.  The following variance-form proxy reproduces that public
    # curve, while making the unknown effective twin mixture explicit.
    source_model_separation = periodic_source_separation(
        linewidth_model[:, 0], source_dispersion)
    source_model_excess = (linewidth_model[:, 1]**2
                           - SOURCE_RESOLUTION_FWHM_MEV**2)
    equal_twin_coefficient = float(np.sum(
        source_model_separation**2*source_model_excess
    )/np.sum(source_model_separation**4))
    linewidth_proxy_source_curve = linewidth_from_separation(
        source_model_separation, equal_twin_coefficient)
    linewidth_proxy_rmse = float(np.sqrt(np.mean(
        (linewidth_proxy_source_curve-linewidth_model[:, 1])**2)))
    linewidth_proxy_max_error = float(np.max(np.abs(
        linewidth_proxy_source_curve-linewidth_model[:, 1])))
    source_data_separation = periodic_source_separation(
        linewidth_data[:, 0], source_dispersion)
    source_proxy_at_data = linewidth_from_separation(
        source_data_separation, equal_twin_coefficient)
    source_proxy_chi2 = float(np.sum(
        ((source_proxy_at_data-linewidth_data[:, 1])
         / linewidth_data[:, 2])**2))

    base = precompute_m3_paths()
    central = {**SOURCE_CENTRAL, "label": "central"}
    corners = uncertainty_corners()
    all_settings = [central, *corners]
    thresholds = {row["label"]: source_threshold(row) for row in all_settings}

    model_rows = []
    curve_rows = []
    for raw in retained:
        evaluations = []
        central_curves = None
        for setting in all_settings:
            curves = model_curves(raw, base, setting)
            bandwidths = [float(np.ptp(curve)) for curve in curves]
            threshold = thresholds[setting["label"]]
            evaluations.append({
                "label": setting["label"],
                "max_domain_bandwidth_microeV": max(bandwidths),
                "source_reference_bandwidth_microeV": threshold,
                "passes_reference": max(bandwidths) <= threshold,
            })
            if setting["label"] == "central":
                central_curves = curves
                central_bandwidths = bandwidths
        assert central_curves is not None
        model_rows.append({
            "sample_id": raw["sample_id"],
            "j2_microev": raw["j2_microev"],
            "ja_microev": raw["ja_microev"],
            "jb_microev": raw["jb_microev"],
            "delta_c": raw["delta_c"],
            "critical_field_shift_vs_M3_T": raw["critical_field_shift_vs_M3_T"],
            "axial_gamma_bandwidth_microeV_existing":
                raw["out_of_plane_bandwidth_microeV"],
            "transverse_domain1_bandwidth_microeV": central_bandwidths[0],
            "transverse_domain2_bandwidth_microeV": central_bandwidths[1],
            "transverse_max_domain_bandwidth_microeV": max(central_bandwidths),
            "source_reference_bandwidth_microeV": public_threshold,
            "central_pass": max(central_bandwidths) <= public_threshold,
            "one_sigma_corner_min_bandwidth_microeV": min(
                item["max_domain_bandwidth_microeV"] for item in evaluations[1:]),
            "one_sigma_corner_max_bandwidth_microeV": max(
                item["max_domain_bandwidth_microeV"] for item in evaluations[1:]),
            "matched_corner_pass_count": sum(
                item["passes_reference"] for item in evaluations[1:]),
            "robust_pass_all_1sigma_corners": all(
                item["passes_reference"] for item in evaluations[1:]),
            "corner_evaluations": evaluations,
        })
        for domain, curve in enumerate(central_curves, start=1):
            for l_value, energy in zip(L_GRID, curve):
                curve_rows.append({
                    "sample_id": raw["sample_id"],
                    "domain": domain,
                    "l_rlu": float(l_value),
                    "energy_microeV": float(energy),
                })

    # Directly compare every M4 pair of structural-domain branches with the 11
    # public linewidth points.  A fixed-coefficient score uses the source's
    # equal-twin calculation.  A second score profiles the effective twin
    # mixture from a single-domain limit (c=0) to that equal-twin limit.  This
    # is deliberately a nuisance audit, not a new confidence construction.
    profile_linewidth_curves = []
    for row in model_rows:
        model_curve_rows = [
            item for item in curve_rows if item["sample_id"] == row["sample_id"]
        ]
        domain_curves = []
        for domain in (1, 2):
            selected = [item for item in model_curve_rows
                        if item["domain"] == domain]
            selected.sort(key=lambda item: item["l_rlu"])
            domain_curves.append(np.asarray(
                [item["energy_microeV"] for item in selected]))
        separation_at_data_mev = np.abs(
            np.interp(linewidth_data[:, 0], L_GRID, domain_curves[0])
            - np.interp(linewidth_data[:, 0], L_GRID, domain_curves[1])
        )/1000.0

        def objective(coefficient: float) -> float:
            predicted = linewidth_from_separation(
                separation_at_data_mev, coefficient)
            return float(np.sum(
                ((predicted-linewidth_data[:, 1])
                 / linewidth_data[:, 2])**2))

        profiled = minimize_scalar(
            objective, bounds=(0.0, equal_twin_coefficient), method="bounded",
            options={"xatol": 1e-12})
        fixed_chi2 = objective(equal_twin_coefficient)
        row["linewidth_equal_twin_chi2"] = fixed_chi2
        row["linewidth_profile_chi2"] = float(profiled.fun)
        row["linewidth_profile_coefficient"] = float(profiled.x)
        row["linewidth_effective_minor_domain_fraction_proxy"] = (
            minor_domain_fraction_proxy(profiled.x, equal_twin_coefficient))
        row["linewidth_profile_no_worse_than_source"] = (
            profiled.fun <= linewidth_chi2)

        separation_curve_mev = np.abs(domain_curves[0]-domain_curves[1])/1000.0
        profile_linewidth_curves.append(linewidth_from_separation(
            separation_curve_mev, profiled.x))

    central_subset_ids = {row["sample_id"] for row in model_rows
                          if row["central_pass"]}
    robust_subset_ids = {row["sample_id"] for row in model_rows
                         if row["robust_pass_all_1sigma_corners"]}
    central_subset = [row for row in retained if row["sample_id"] in central_subset_ids]
    robust_subset = [row for row in retained if row["sample_id"] in robust_subset_ids]
    if not central_subset or not robust_subset:
        raise RuntimeError("INS sensitivity screen removed every model")

    representative = {
        "sample_id": -1,
        **asdict(InterlayerParameters()),
        "critical_field_shift_vs_M3_T": 0.0952738089552219,
    }
    rep_gamma_j0 = base["j0"] + interlayer_fourier(
        np.zeros(3), InterlayerParameters())
    rep_gamma_jk = base["gamma"] + np.asarray([
        interlayer_fourier(np.array([0.0, 0.0, l]), InterlayerParameters())
        for l in L_GRID
    ])
    representative_axial_bandwidth = float(np.ptp(
        axial_energy_microev(rep_gamma_jk, rep_gamma_j0)))

    bandwidth = np.asarray([
        row["transverse_max_domain_bandwidth_microeV"] for row in model_rows])
    shifts = np.asarray([
        row["critical_field_shift_vs_M3_T"] for row in model_rows])
    rho, pvalue = spearmanr(bandwidth, shifts)
    profile_chi2_values = np.asarray([
        row["linewidth_profile_chi2"] for row in model_rows])
    equal_twin_chi2_values = np.asarray([
        row["linewidth_equal_twin_chi2"] for row in model_rows])
    linewidth_rho, linewidth_pvalue = spearmanr(profile_chi2_values, shifts)
    source_overdispersion = linewidth_chi2/len(linewidth_data)
    scaled_profile_threshold_95 = float(
        chi2.ppf(0.95, len(linewidth_data)-1)*source_overdispersion)
    profile_no_worse_ids = {
        row["sample_id"] for row in model_rows
        if row["linewidth_profile_no_worse_than_source"]
    }
    profile_no_worse_subset = [
        row for row in retained if row["sample_id"] in profile_no_worse_ids]

    summary = {
        "status": "independent public-INS sensitivity audit; not a confidence filter",
        "source": {
            "paper": "Phys. Rev. B 112, 104413 (2025)",
            "paper_doi": "10.1103/1pvl-kzjm",
            "data_doi": "10.5287/ora-qrd6mzojz",
            "field_T": FIELD_T,
            "experimental_path_hkl": "(-2/3, 1/3, l)",
            "resolution_fwhm_meV": SOURCE_RESOLUTION_FWHM_MEV,
            "source_simplified_Jb_microeV": 1000*SOURCE_JB_MEV,
            "source_language": "indicative upper bound",
            "public_branch_bandwidths_microeV": source_branch_bandwidths.tolist(),
            "declared_sensitivity_threshold_microeV": public_threshold,
            "source_formula_max_abs_error_meV_after_branch_permutation":
                source_formula_error,
            "linewidth_model_vs_data": {
                "point_count": len(linewidth_data),
                "chi_square_no_refit": linewidth_chi2,
                "chi_square_per_point_descriptive": linewidth_chi2/len(linewidth_data),
                "mae_meV": float(np.mean(np.abs(linewidth_residual))),
                "rmse_meV": float(np.sqrt(np.mean(linewidth_residual**2))),
                "interpretation": (
                    "The source model is an approximate upper-reference, not a "
                    "statistically calibrated acceptance model."
                ),
            },
            "linewidth_proxy_reproduction": {
                "form": "FWHM^2 = resolution^2 + c * twin_separation^2",
                "equal_twin_coefficient": equal_twin_coefficient,
                "source_curve_rmse_meV": linewidth_proxy_rmse,
                "source_curve_max_abs_error_meV": linewidth_proxy_max_error,
                "source_proxy_chi_square_at_11_data_points":
                    source_proxy_chi2,
            },
        },
        "implementation_cross_check": {
            "representative_M4_axial_Gamma_bandwidth_microeV":
                representative_axial_bandwidth,
            "published_2026_value_microeV": 4.4,
            "absolute_difference_microeV": abs(
                representative_axial_bandwidth-4.4),
        },
        "screen": {
            "input_model_count": len(retained),
            "path_mapping": "source (h,k)=(-2/3,1/3) -> code (u,v)=(1/3,2/3)",
            "structural_domain_rule": "maximum bandwidth over two K-domain paths",
            "central_pass_count": len(central_subset),
            "central_rejected_count": len(retained)-len(central_subset),
            "central_bandwidth_range_microeV": [float(np.min(bandwidth)),
                                                  float(np.max(bandwidth))],
            "central_pass_shift_range_T": range_of(
                central_subset, "critical_field_shift_vs_M3_T"),
            "robust_all_1sigma_corners_pass_count": len(robust_subset),
            "robust_all_1sigma_corners_shift_range_T": range_of(
                robust_subset, "critical_field_shift_vs_M3_T"),
            "source_threshold_range_over_1sigma_corners_microeV": [
                min(thresholds[row["label"]] for row in corners),
                max(thresholds[row["label"]] for row in corners),
            ],
            "spearman_bandwidth_vs_critical_field_shift": float(rho),
            "spearman_pvalue_descriptive_only": float(pvalue),
        },
        "linewidth_twin_profile_audit": {
            "status": (
                "direct public-linewidth nuisance audit; does not define a "
                "confidence filter"),
            "fixed_resolution_fwhm_meV": SOURCE_RESOLUTION_FWHM_MEV,
            "profiled_nuisance": (
                "effective structural-domain mixture coefficient between "
                "single-domain and source equal-twin limits"),
            "equal_twin_chi_square_range": [
                float(np.min(equal_twin_chi2_values)),
                float(np.max(equal_twin_chi2_values))],
            "equal_twin_models_no_worse_than_source_count": int(np.sum(
                equal_twin_chi2_values <= linewidth_chi2)),
            "profile_chi_square_range": [
                float(np.min(profile_chi2_values)),
                float(np.max(profile_chi2_values))],
            "profile_models_no_worse_than_source_count":
                len(profile_no_worse_subset),
            "profile_no_worse_shift_range_T": range_of(
                profile_no_worse_subset, "critical_field_shift_vs_M3_T"),
            "effective_minor_domain_fraction_proxy_range": range_of(
                model_rows,
                "linewidth_effective_minor_domain_fraction_proxy"),
            "maximum_profile_delta_chi_square_across_models": float(
                np.max(profile_chi2_values)-np.min(profile_chi2_values)),
            "spearman_profile_chi_square_vs_critical_field_shift":
                float(linewidth_rho),
            "spearman_pvalue_descriptive_only": float(linewidth_pvalue),
            "source_chi_square_per_point_overdispersion":
                source_overdispersion,
            "source_scaled_95pct_profile_diagnostic_threshold":
                scaled_profile_threshold_95,
            "models_below_source_scaled_95pct_profile_diagnostic": int(np.sum(
                profile_chi2_values <= scaled_profile_threshold_95)),
            "interpretation": (
                "The direct linewidth comparison does not support converting "
                "the 30-model bandwidth sensitivity set into a calibrated "
                "rejection. Profiling the unresolved domain mixture leaves all "
                "69 models no worse than the source reference chi-square."),
        },
        "engineering_consequence": {
            "full_69_shift_range_T": range_of(
                retained, "critical_field_shift_vs_M3_T"),
            "central_sensitivity_subset_90mK": engineering_summary(central_subset),
            "central_sensitivity_subset_shape_stress_90mK":
                shape_summary(central_subset),
            "robust_1sigma_subset_90mK": engineering_summary(robust_subset),
            "robust_1sigma_subset_shape_stress_90mK":
                shape_summary(robust_subset),
        },
        "interpretation_limits": [
            "The 18.048 microeV threshold is read from the source's Jb=3 microeV calculation; the source calls that coupling an indicative upper bound.",
            "The source linewidth model is approximate and its calculated energy modulation was not resolved experimentally.",
            "The M4 Ja/Jb/J2 geometry is more general than the source's Jb-only model, so this is an equivalent observable-bandwidth comparison, not a direct parameter bound.",
            "The spectroscopy-compatible subsets are sensitivity envelopes, not posterior or confidence sets; the authoritative engineering envelope remains all 69 models.",
            "A calibrated proxy reproduces the source linewidth curve, but profiling the unresolved structural-domain mixture leaves all 69 models no worse than the source reference; the public linewidth points therefore do not validate a 39-model rejection.",
            "Linear spin-wave theory is used only in the 3.5 T polarized phase where the source reports sharp, resolution-limited magnons.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    write_csv(RESULTS / "public_ins_interlayer_bandwidth_models.csv", model_rows)
    write_csv(RESULTS / "public_ins_interlayer_bandwidth_curves.csv", curve_rows)
    write_csv(RESULTS / "public_ins_interlayer_bandwidth_central_subset.csv",
              [row for row in model_rows if row["central_pass"]])
    write_csv(RESULTS / "public_ins_interlayer_bandwidth_robust_subset.csv",
              [row for row in model_rows if row["robust_pass_all_1sigma_corners"]])
    (RESULTS / "public_ins_interlayer_bandwidth_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    provenance = {
        "retrieved_utc_date": "2026-10-12",
        "official_record_url":
            "https://ora.ox.ac.uk/objects/uuid:b4da68a2-fdd2-42e9-80cd-1f13120b4790",
        "paper_url": "https://doi.org/10.1103/1pvl-kzjm",
        "data_doi": "https://doi.org/10.5287/ora-qrd6mzojz",
        "archive_download_url":
            "https://ora.ox.ac.uk/objects/uuid:b4da68a2-fdd2-42e9-80cd-1f13120b4790/files/rn870zs79x",
        "files": {
            str(SOURCE_ARCHIVE.relative_to(ROOT)): sha256(SOURCE_ARCHIVE),
            str(SOURCE_DISPERSION.relative_to(ROOT)): sha256(SOURCE_DISPERSION),
            str(SOURCE_LINEWIDTH_DATA.relative_to(ROOT)): sha256(SOURCE_LINEWIDTH_DATA),
            str(SOURCE_LINEWIDTH_MODEL.relative_to(ROOT)): sha256(SOURCE_LINEWIDTH_MODEL),
        },
        "redistribution": (
            "External source files retain their original terms and are excluded "
            "from the public release; URLs and hashes permit independent retrieval."
        ),
    }
    (SOURCES / "public_ins_interlayer_bandwidth_provenance_20261012.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.5))
    ax = axes[0, 0]
    ax.plot(source_dispersion[:, 0], source_dispersion[:, 1]*1000,
            color="#4c78a8", lw=1.8, label="source twin branch 1")
    ax.plot(source_dispersion[:, 0], source_dispersion[:, 2]*1000,
            color="#f58518", lw=1.8, label="source twin branch 2")
    lower = float(np.min(source_dispersion[:, 1:])*1000)
    upper = float(np.max(source_dispersion[:, 1:])*1000)
    ax.annotate("", xy=(-0.82, upper), xytext=(-0.82, lower),
                arrowprops={"arrowstyle": "<->", "color": "black"})
    ax.text(-0.79, 0.5*(lower+upper),
            f"{public_threshold:.3f} microeV", va="center", fontsize=8)
    ax.set(xlabel=r"$l$ (r.l.u.)", ylabel="energy (microeV)",
           title="(a) Public 3.5 T upper-reference calculation")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[0, 1]
    ax.errorbar(linewidth_data[:, 0], linewidth_data[:, 1]*1000,
                yerr=linewidth_data[:, 2]*1000, fmt="o", ms=4.5,
                color="black", capsize=2, label="public data")
    ax.plot(linewidth_model[:, 0], linewidth_model[:, 1]*1000,
            color="#e45756", lw=1.8, label="source Jb=3 microeV model")
    profile_linewidth_curves = np.asarray(profile_linewidth_curves)
    ax.fill_between(L_GRID,
                    np.min(profile_linewidth_curves, axis=0)*1000,
                    np.max(profile_linewidth_curves, axis=0)*1000,
                    color="#4c78a8", alpha=0.18,
                    label="69 M4, twin mixture profiled")
    ax.plot(L_GRID, np.median(profile_linewidth_curves, axis=0)*1000,
            color="#4c78a8", lw=1.1, alpha=0.9)
    ax.axhline(SOURCE_RESOLUTION_FWHM_MEV*1000, ls="--", lw=1,
               color="#777777", label="instrument resolution")
    ax.text(0.03, 0.95,
            (rf"source $\chi^2={linewidth_chi2:.2f}$; "
             rf"profile: {len(profile_no_worse_subset)}/69 no worse"),
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set(xlabel=r"$l$ (r.l.u.)", ylabel="FWHM (microeV)",
           title="(b) Direct linewidth nuisance audit")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    passed = np.asarray([row["central_pass"] for row in model_rows], dtype=bool)
    ax.scatter(shifts[~passed], bandwidth[~passed], s=32, color="#e45756",
               alpha=0.78, label="above reference")
    ax.scatter(shifts[passed], bandwidth[passed], s=34, color="#54a24b",
               edgecolor="black", linewidth=0.25, label="at/below reference")
    ax.axhline(public_threshold, color="black", ls="--", lw=1.2,
               label=f"{public_threshold:.3f} microeV")
    ax.text(0.03, 0.95, rf"Spearman $\rho={rho:.3f}$",
            transform=ax.transAxes, va="top")
    ax.set(xlabel=r"M4--M3 $\Delta H_c$ (T)",
           ylabel="maximum domain bandwidth (microeV)",
           title="(c) Independent INS stress of 69 M4 models")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    full_range = range_of(retained, "critical_field_shift_vs_M3_T")
    central_range = range_of(central_subset, "critical_field_shift_vs_M3_T")
    robust_range = range_of(robust_subset, "critical_field_shift_vs_M3_T")
    profile_range = range_of(profile_no_worse_subset,
                             "critical_field_shift_vs_M3_T")
    for y, values, color, label in [
        (2, full_range, "#9d9d9d", "diffraction: 69"),
        (1, central_range, "#54a24b",
         f"bandwidth sensitivity: {len(central_subset)}"),
        (0, profile_range, "#4c78a8",
         f"direct linewidth profile: {len(profile_no_worse_subset)}"),
    ]:
        ax.plot(values, [y, y], color=color, lw=7, solid_capstyle="round")
        ax.scatter(values, [y, y], color=color, s=40, zorder=3)
        ax.text(values[1]+0.003, y, label, va="center", fontsize=9)
    ax.set_yticks([])
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlim(0.02, max(full_range[1]+0.055, 0.22))
    ax.set(xlabel=r"M4--M3 $\Delta H_c$ range (T)",
           title="(d) Sensitivity contraction is not a rejection")
    ax.grid(axis="x", alpha=0.25)

    fig.suptitle("Independent public-INS bandwidth and linewidth audit",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"public_ins_interlayer_bandwidth_audit.{suffix}",
                    dpi=230, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
