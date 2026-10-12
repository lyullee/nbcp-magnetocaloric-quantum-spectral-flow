"""Validate a lumped ADR heat balance against an independent public device.

The Cambridge repository supplies zero-field heat capacity and warm-up curves
for the same approximately 15 g YbNi1.6Sn pill under labelled 2.5 and 5 uW
loads.  The low-temperature part of both curves calibrates one common pill mass
and one common parasitic load.  A disjoint higher-temperature interval is then
held out.  Raw third-party data are downloaded from their DOI repository and
are not redistributed by this project.
"""
from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "sources" / "raw"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
PROVENANCE_PATH = ROOT / "sources" / "ybni16sn_stage_validation_provenance_20261012.json"
ZIP_PATH = RAW / "ybni16sn_supporting_data.zip"
DATA_DIR = RAW / "ybni16sn_supporting_data_v1"

MOLAR_MASS_G_PER_MOL_YB = 173.045 + 1.6 * 58.6934 + 118.710
LOAD_CURVES = {
    2.5: "Fig04b_172mK.dat",
    5.0: "Fig04b_192mK.dat",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def ensure_source(provenance: dict) -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    expected = provenance["dataset"]["archive_sha256"]
    if not ZIP_PATH.exists():
        urllib.request.urlretrieve(provenance["dataset"]["bitstream_url"], ZIP_PATH)
    actual = sha256(ZIP_PATH)
    if actual != expected:
        raise RuntimeError(f"source archive SHA256 mismatch: {actual} != {expected}")
    needed = ["Fig02_0p0T.dat", *LOAD_CURVES.values()]
    if not all((DATA_DIR / name).exists() for name in needed):
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(ZIP_PATH) as archive:
            archive.extractall(DATA_DIR)


def read_warm_branch(path: Path) -> tuple[np.ndarray, np.ndarray, float, float]:
    data = np.genfromtxt(path, skip_header=2)
    data = data[np.all(np.isfinite(data), axis=1)]
    minimum = int(np.argmin(data[:, 1]))
    time = data[minimum:, 0]
    temperature = data[minimum:, 1]
    # Use first-passage times through a monotone envelope. This prevents the
    # RuO2 readout noise from creating multiple times for one temperature.
    monotone_temperature = np.maximum.accumulate(temperature)
    unique_temperature, first = np.unique(monotone_temperature, return_index=True)
    return (unique_temperature, time[first],
            float(data[minimum, 0]), float(data[minimum, 1]))


def molar_enthalpy_function(heat_capacity_path: Path, alpha: float,
                            low: float, high: float):
    data = np.genfromtxt(heat_capacity_path, skip_header=3)
    data = data[np.all(np.isfinite(data), axis=1)]
    dense_t = np.linspace(low, high, 20001)
    electronic = np.interp(dense_t, data[:, 0], data[:, 1])
    nuclear = alpha / dense_t**2
    total = electronic + nuclear
    cumulative = np.r_[
        0.0,
        np.cumsum(0.5 * (total[:-1] + total[1:]) * np.diff(dense_t)),
    ]

    def enthalpy(temperature):
        return np.interp(np.asarray(temperature), dense_t, cumulative)

    return enthalpy


def fit_balance(curves: dict, enthalpy, reference: float,
                calibration_upper: float) -> tuple[float, float]:
    grid = np.linspace(reference, calibration_upper, 81)
    design = []
    target = []
    for applied_load, (temperature, time, *_unused) in curves.items():
        elapsed = np.interp(grid, temperature, time) - np.interp(
            reference, temperature, time)
        # 1 uW min = 0.06 mJ.  k*DeltaH is the material enthalpy in mJ.
        design.extend(np.c_[enthalpy(grid), -0.06 * elapsed])
        target.extend(0.06 * applied_load * elapsed)
    k_mJ_per_J_per_mol, parasitic_uW = np.linalg.lstsq(
        np.asarray(design), np.asarray(target), rcond=None)[0]
    effective_mass_g = (
        k_mJ_per_J_per_mol * MOLAR_MASS_G_PER_MOL_YB / 1000.0)
    return float(effective_mass_g), float(parasitic_uW)


def error_metrics(rows: list[dict], split: str) -> dict:
    selected = [row for row in rows if row["role"] == split]
    residual = np.asarray([row["residual_mJ"] for row in selected])
    scale = max(row["predicted_material_enthalpy_mJ"] for row in rows)
    return {
        "row_count": len(selected),
        "rmse_mJ": float(np.sqrt(np.mean(residual**2))),
        "nrmse_vs_full_interval_enthalpy": float(np.sqrt(np.mean(residual**2)) / scale),
        "maximum_absolute_error_mJ": float(np.max(np.abs(residual))),
        "maximum_absolute_error_fraction_of_full_interval_enthalpy":
            float(np.max(np.abs(residual)) / scale),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(rows: list[dict], summary: dict) -> None:
    colors = {2.5: "#4c78a8", 5.0: "#e45756"}
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.25))
    for load in LOAD_CURVES:
        subset = [row for row in rows if row["applied_load_uW"] == load]
        axes[0].plot([row["elapsed_from_0p21K_min"] for row in subset],
                     [row["temperature_K"] for row in subset],
                     color=colors[load], label=f"{load:g} $\\mu$W")
        axes[1].plot([row["temperature_K"] for row in subset],
                     [row["measured_input_energy_with_fitted_parasitic_mJ"]
                      for row in subset], color=colors[load],
                     label=f"{load:g} $\\mu$W data")
        axes[2].plot([row["temperature_K"] for row in subset],
                     [100 * row["residual_fraction_of_full_interval_enthalpy"]
                      for row in subset], color=colors[load],
                     label=f"{load:g} $\\mu$W")

    first = [row for row in rows if row["applied_load_uW"] == 2.5]
    axes[1].plot([row["temperature_K"] for row in first],
                 [row["predicted_material_enthalpy_mJ"] for row in first],
                 color="black", linestyle="--", linewidth=1.6,
                 label="heat-capacity prediction")
    for ax in axes[1:]:
        ax.axvline(summary["protocol"]["calibration_upper_temperature_K"],
                   color="#888888", linestyle=":", linewidth=1.2)
    axes[2].axhline(0, color="black", linewidth=0.8)
    axes[0].set(xlabel="Time after first reaching 0.21 K (min)",
                ylabel="Stage temperature (K)",
                title="Independent device warm-up")
    axes[1].set(xlabel="Stage temperature (K)",
                ylabel="Accumulated energy from 0.21 K (mJ)",
                title="Energy curves collapse")
    axes[2].set(xlabel="Stage temperature (K)",
                ylabel="Residual / full-interval enthalpy (%)",
                title="Prediction residual")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.18)
    fig.suptitle(
        f"YbNi$_{{1.6}}$Sn validation: fitted mass "
        f"{summary['fit']['effective_mass_g']:.2f} g, parasitic load "
        f"{summary['fit']['parasitic_load_uW']:.3f} $\\mu$W",
        fontsize=11)
    fig.tight_layout()
    for extension in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"ybni16sn_stage_heat_balance_validation.{extension}",
                    dpi=240 if extension == "png" else None)
    plt.close(fig)


def main() -> None:
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))
    ensure_source(provenance)
    protocol = provenance["declared_analysis_protocol"]
    reference = float(protocol["temperature_reference_K"])
    calibration_upper = float(protocol["calibration_upper_temperature_K"])
    holdout_upper = float(protocol["holdout_upper_temperature_K"])
    alpha = float(provenance["supplement"][
        "zero_field_nuclear_heat_capacity_coefficient_J_K_per_mol_Yb"])

    curves = {
        load: read_warm_branch(DATA_DIR / filename)
        for load, filename in LOAD_CURVES.items()
    }
    enthalpy = molar_enthalpy_function(
        DATA_DIR / "Fig02_0p0T.dat", alpha, reference, holdout_upper)
    effective_mass_g, parasitic_uW = fit_balance(
        curves, enthalpy, reference, calibration_upper)
    k = effective_mass_g / MOLAR_MASS_G_PER_MOL_YB * 1000.0

    grid = np.linspace(reference, holdout_upper, 141)
    rows = []
    energy_profiles = {}
    for applied_load, (temperature, time, minimum_time, minimum_temperature) in curves.items():
        elapsed = np.interp(grid, temperature, time) - np.interp(
            reference, temperature, time)
        measured_energy = 0.06 * (applied_load + parasitic_uW) * elapsed
        predicted = k * enthalpy(grid)
        energy_profiles[applied_load] = measured_energy
        for temp, dt, observed, expected in zip(
                grid, elapsed, measured_energy, predicted, strict=True):
            residual = observed - expected
            rows.append({
                "applied_load_uW": applied_load,
                "temperature_K": float(temp),
                "role": "calibration" if temp <= calibration_upper + 1e-12 else "holdout",
                "elapsed_from_0p21K_min": float(dt),
                "measured_input_energy_with_fitted_parasitic_mJ": float(observed),
                "predicted_material_enthalpy_mJ": float(expected),
                "residual_mJ": float(residual),
                "residual_fraction_of_full_interval_enthalpy": float(
                    residual / (k * enthalpy(holdout_upper))),
                "raw_curve_minimum_time_min": minimum_time,
                "raw_curve_minimum_temperature_K": minimum_temperature,
            })

    calibration = error_metrics(rows, "calibration")
    holdout = error_metrics(rows, "holdout")
    profile_difference = energy_profiles[2.5] - energy_profiles[5.0]
    full_scale = k * enthalpy(holdout_upper)

    # Physical ablation: omit the public zero-field nuclear heat capacity and
    # repeat exactly the same calibration/holdout protocol.
    electronic_only_enthalpy = molar_enthalpy_function(
        DATA_DIR / "Fig02_0p0T.dat", 0.0, reference, holdout_upper)
    electronic_mass, electronic_parasitic = fit_balance(
        curves, electronic_only_enthalpy, reference, calibration_upper)
    electronic_k = electronic_mass / MOLAR_MASS_G_PER_MOL_YB * 1000.0
    electronic_residual = []
    electronic_scale = electronic_k * electronic_only_enthalpy(holdout_upper)
    holdout_grid = grid[grid > calibration_upper + 1e-12]
    for applied_load, (temperature, time, *_unused) in curves.items():
        elapsed = np.interp(holdout_grid, temperature, time) - np.interp(
            reference, temperature, time)
        observed = 0.06 * (applied_load + electronic_parasitic) * elapsed
        expected = electronic_k * electronic_only_enthalpy(holdout_grid)
        electronic_residual.extend(observed - expected)

    summary = {
        "status": "independent public-device heat-balance validation",
        "source_archive_sha256_verified": sha256(ZIP_PATH),
        "protocol": protocol,
        "model": {
            "identity": "(P_applied + P_parasitic) dt = (m_eff/M) integral(C_elec + alpha/T^2) dT",
            "molar_mass_g_per_mol_Yb": MOLAR_MASS_G_PER_MOL_YB,
            "zero_field_nuclear_alpha_J_K_per_mol_Yb": alpha,
            "load_curves_uW": list(LOAD_CURVES),
        },
        "fit": {
            "calibration_interval_K": [reference, calibration_upper],
            "effective_mass_g": effective_mass_g,
            "published_approximate_mass_g": 15.0,
            "relative_mass_difference_vs_15g": (effective_mass_g - 15.0) / 15.0,
            "parasitic_load_uW": parasitic_uW,
            "published_parasitic_upper_bound_uW": 0.4,
        },
        "calibration_metrics": calibration,
        "holdout_metrics": holdout,
        "cross_load_energy_collapse": {
            "rmse_mJ": float(np.sqrt(np.mean(profile_difference**2))),
            "nrmse_vs_full_interval_enthalpy": float(
                np.sqrt(np.mean(profile_difference**2)) / full_scale),
            "maximum_absolute_difference_mJ": float(np.max(np.abs(profile_difference))),
            "maximum_absolute_difference_fraction_of_full_interval_enthalpy":
                float(np.max(np.abs(profile_difference)) / full_scale),
        },
        "nuclear_heat_capacity_ablation": {
            "effective_mass_g_without_Cn": electronic_mass,
            "parasitic_load_uW_without_Cn": electronic_parasitic,
            "holdout_nrmse_without_Cn": float(
                np.sqrt(np.mean(np.asarray(electronic_residual)**2)) / electronic_scale),
            "holdout_nrmse_with_Cn": holdout["nrmse_vs_full_interval_enthalpy"],
        },
        "validation_gates": {
            "fitted_mass_within_10_percent_of_approximately_15g":
                bool(abs(effective_mass_g - 15.0) <= 1.5),
            "fitted_parasitic_within_published_less_than_0p4uW":
                bool(0.0 <= parasitic_uW < 0.4),
            "holdout_nrmse_below_2_percent":
                bool(holdout["nrmse_vs_full_interval_enthalpy"] < 0.02),
            "holdout_maximum_error_below_5_percent":
                bool(holdout["maximum_absolute_error_fraction_of_full_interval_enthalpy"] < 0.05),
            "cross_load_collapse_nrmse_below_2_percent":
                bool(np.sqrt(np.mean(profile_difference**2)) / full_scale < 0.02),
            "nuclear_term_improves_holdout_nrmse":
                bool(holdout["nrmse_vs_full_interval_enthalpy"]
                     < np.sqrt(np.mean(np.asarray(electronic_residual)**2)) / electronic_scale),
        },
        "interpretation": (
            "The independent device reproduces the lumped heat-balance identity and "
            "recovers its separately reported approximate mass and heat-leak bound. "
            "This validates thermodynamic-to-stage translation, not NBCP-specific "
            "package losses or the M4 entropy surrogate."),
    }
    summary["all_validation_gates_pass"] = all(summary["validation_gates"].values())

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    write_csv(RESULTS / "ybni16sn_stage_heat_balance_profiles.csv", rows)
    (RESULTS / "ybni16sn_stage_heat_balance_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    make_figure(rows, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
