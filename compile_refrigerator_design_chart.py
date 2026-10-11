"""Translate robust magnetic heat budgets into auditable design quantities."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
N_A = 6.02214076e23
MOLAR_MASS_G_MOL = (
    2 * 22.98976928 + 137.327 + 58.933194 + 2 * 30.973761998
    + 8 * 15.999
)
A_ANGSTROM = 5.312
C_ANGSTROM = 6.983
CELL_VOLUME_A3 = np.sqrt(3.0) / 2.0 * A_ANGSTROM**2 * C_ANGSTROM
DENSITY_G_CM3 = MOLAR_MASS_G_MOL / N_A / (CELL_VOLUME_A3 * 1e-24)
PUBLIC_SAMPLE_MASS_G = 1.65


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    robust = list(csv.DictReader(
        (RESULTS / "thermodynamic_shape_robust_decisions.csv").open(
            encoding="utf-8-sig")))
    rows = []
    sample_rows = []
    public_sample_volume = PUBLIC_SAMPLE_MASS_G / DENSITY_G_CM3
    for decision in robust:
        temperature = float(decision["temperature_K"])
        field = float(decision["robust_field_T"])
        q_worst = float(decision["worst_case_Q_at_robust_field_mJ_per_cm3"])
        w95_exists = decision["guaranteed_W95_exists"].lower() == "true"
        if w95_exists:
            wmin = float(decision["guaranteed_W95_min_T"])
            wmax = float(decision["guaranteed_W95_max_T"])
            tolerance_t = min(field - wmin, wmax - field)
        else:
            wmin = wmax = tolerance_t = float("nan")
        for load_uw in (0.05, 0.10, 0.50):
            for hold_min in (30, 60):
                required_mj = load_uw * 1e-3 * hold_min * 60.0
                volume = required_mj / q_worst
                rows.append({
                    "temperature_K": temperature,
                    "robust_field_T": field,
                    "guaranteed_W95_min_T": wmin,
                    "guaranteed_W95_max_T": wmax,
                    "two_sided_field_tolerance_mT": 1000.0 * tolerance_t,
                    "two_sided_field_tolerance_percent":
                        100.0 * tolerance_t / field,
                    "worst_case_Q_at_robust_field_mJ_per_cm3": q_worst,
                    "assumed_total_cold_load_uW": load_uw,
                    "required_hold_minutes": hold_min,
                    "required_heat_absorption_mJ": required_mj,
                    "conditional_material_volume_cm3": volume,
                    "conditional_material_mass_g": volume * DENSITY_G_CM3,
                })
        for load_uw in (0.05, 0.10, 0.50):
            hold = q_worst * public_sample_volume / (load_uw * 1e-3) / 60.0
            sample_rows.append({
                "temperature_K": temperature,
                "robust_field_T": field,
                "public_sample_mass_g": PUBLIC_SAMPLE_MASS_G,
                "lattice_derived_sample_volume_cm3": public_sample_volume,
                "assumed_total_cold_load_uW": load_uw,
                "conditional_hold_minutes": hold,
            })

    summary = {
        "status": "magnetic-only conditional design chart",
        "density_derivation": {
            "formula": "Na2BaCo(PO4)2",
            "molar_mass_g_mol": MOLAR_MASS_G_MOL,
            "lattice_a_angstrom": A_ANGSTROM,
            "lattice_c_angstrom": C_ANGSTROM,
            "assumed_formula_units_per_primitive_cell": 1,
            "primitive_cell_volume_angstrom3": CELL_VOLUME_A3,
            "density_g_cm3": DENSITY_G_CM3,
            "source": "public P-3 one-Co primitive magnetic lattice and a,c values in arXiv:2607.23173v1 lines 63, 78, 363-402",
        },
        "public_demagnetization_sample_crosscheck": {
            "source": "Xiang et al. Nature 625, 270-275 (2024), public Extended Data Fig. 6b workbook",
            "reported_NBCP_mass_g": PUBLIC_SAMPLE_MASS_G,
            "lattice_derived_volume_cm3": public_sample_volume,
        },
        "design_rows": len(rows),
        "scope": [
            "Loads are assumed total cold loads supplied by the designer, not measured parasitic loads.",
            "The heat budget is magnetic and reversible; addenda, heat switches, wiring, contact resistance and ramp loss are not separately predicted.",
            "Mass uses a lattice-derived crystallographic density and does not include packing voids, binder or thermal bus material.",
            "Field tolerance is the symmetric distance from the minimax field to the guaranteed 95%-efficiency interval under the declared shape envelope.",
        ],
    }

    write_csv(RESULTS / "refrigerator_design_chart.csv", rows)
    write_csv(RESULTS / "public_sample_conditional_hold_chart.csv", sample_rows)
    (RESULTS / "refrigerator_design_chart_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    temperatures = sorted({float(r["temperature_K"]) for r in rows})
    base = [next(r for r in rows if float(r["temperature_K"]) == t)
            for t in temperatures]
    fields = np.array([float(r["robust_field_T"]) for r in base])
    tolerance = np.array([float(r["two_sided_field_tolerance_mT"])
                          for r in base]) / 1000.0
    axes[0].errorbar(np.array(temperatures) * 1000.0, fields, yerr=tolerance,
                     fmt="o-", capsize=4, color="#4c78a8")
    missing = ~np.isfinite(tolerance)
    axes[0].scatter(np.array(temperatures)[missing] * 1000.0,
                    fields[missing], marker="x", s=55, color="#e45756",
                    label="no common $W_{95}$")
    axes[0].set(xlabel="Temperature (mK)", ylabel="Field (T)",
                title="Setpoint and guaranteed 95% tolerance")
    axes[0].legend(fontsize=8)

    for load in (0.05, 0.10, 0.50):
        picked = [next(r for r in rows
                       if float(r["temperature_K"]) == t
                       and float(r["assumed_total_cold_load_uW"]) == load
                       and int(r["required_hold_minutes"]) == 60)
                  for t in temperatures]
        axes[1].plot(np.array(temperatures) * 1000.0,
                     [float(r["conditional_material_mass_g"]) for r in picked],
                     "o-", label=rf"{load:.2f} $\mu$W")
    axes[1].set(xlabel="Temperature (mK)", ylabel="NBCP mass (g)",
                title="60 min magnetic-only mass")
    axes[1].legend(fontsize=8)

    for load in (0.05, 0.10, 0.50):
        picked = [r for r in sample_rows
                  if float(r["assumed_total_cold_load_uW"]) == load]
        axes[2].plot([1000 * float(r["temperature_K"]) for r in picked],
                     [float(r["conditional_hold_minutes"]) for r in picked],
                     "o-", label=rf"{load:.2f} $\mu$W")
    axes[2].set(xlabel="Temperature (mK)", ylabel="Hold time (min)",
                title="Public 1.65 g sample equivalent")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES / "refrigerator_design_chart.png", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
