"""Public-data P0--P2 interlayer scenarios and engineering consequences.

P0 uses the converged 2D M3 polarized-state instability.  P1 adds the older
public DFT estimate J_perp=0.05 K as a simple isotropic vertical bond.  P2 uses
the 2026 representative interlayer model already evaluated as M4.  Finite-T
curves are transported by reduced-field and critical-offset mappings; this is
a decision-scale uncertainty propagation, not a 3D finite-T solution.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bulk_to_engineering_surrogate import inverse_field, transform_field
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
JPERP_DFT_K = 0.05
P1_HC_SHIFT_T = 2 * JPERP_DFT_K / (GZ * MU_B_OVER_KB)
METHODS = ("reduced_field", "critical_offset")
EXPERIMENTAL_HC = {
    "neutron_diffraction_65mK": 1.65,
    "bulk_magnetization_80mK": 1.70,
}


def write_csv(name, rows):
    with (RESULTS / name).open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    bulk = json.loads((RESULTS / "bulk_polarized_stability.json").read_text(
        encoding="utf-8"))
    lookup = {(row["model"], row["boundary"]): row
              for row in bulk["rows"]}
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)
    profile = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",", names=True, encoding="utf-8-sig")
    profile = profile[profile["anchor_N"] == 12]

    scenarios = {}
    for boundary in ("tin_foil", "sphere_vacuum"):
        p0 = float(lookup[("M3", boundary)]["critical_field_T"])
        p2 = float(lookup[("M4", boundary)]["critical_field_T"])
        scenarios[("P0", boundary)] = p0
        scenarios[("P1", boundary)] = p0 + P1_HC_SHIFT_T
        scenarios[("P2", boundary)] = p2

    hc_rows = []
    for (scenario, boundary), hc in scenarios.items():
        hc_rows.append({
            "scenario": scenario,
            "boundary": boundary,
            "critical_field_T": hc,
            "delta_Hc_vs_P0_same_boundary_T":
                hc - scenarios[("P0", boundary)],
        })

    p2_shifts = [row["delta_Hc_vs_P0_same_boundary_T"]
                 for row in hc_rows if row["scenario"] == "P2"]
    if np.ptp(p2_shifts) > 1e-10:
        raise RuntimeError("P2-P0 shift unexpectedly depends on Ewald boundary")
    p2_shift = float(np.mean(p2_shifts))
    # The full theoretical model overestimates Hc by about 9%.  Use the public
    # experimental Hc as the P2 endpoint, and transport only the calculated
    # P0/P1 ablation shifts.  This avoids treating 1.86--1.92 T as a device
    # operating field.
    anchored = {}
    for anchor, hc_exp in EXPERIMENTAL_HC.items():
        anchored[("P2", anchor)] = hc_exp
        anchored[("P0", anchor)] = hc_exp - p2_shift
        anchored[("P1", anchor)] = hc_exp - p2_shift + P1_HC_SHIFT_T

    metric_rows = []
    regret_rows = []
    robust_rows = []
    design_rows = []
    # The critical-field transport is a local near-saturation approximation.
    # Above 150 mK the public optimum moves away from that window, so 200--300
    # mK is reported from the direct experimental candidate curves instead.
    for temperature in (0.09, 0.10, 0.12, 0.15):
        p = np.sort(profile[np.isclose(profile["temperature_K"], temperature)],
                    order="field_T")
        bs = p["field_T"]
        q = p["Q_proxy_mJ_per_cm3"]
        k = int(np.argmax(q))
        source_bstar = float(bs[k])
        qmax = float(q[k])
        accepted = bs[q >= 0.95 * qmax]

        for method in METHODS:
            for (scenario, anchor), target_hc in anchored.items():
                mapped = transform_field(
                    np.array([source_bstar, accepted.min(), accepted.max()]),
                    source_hc, target_hc, method)
                metric_rows.append({
                    "temperature_K": temperature,
                    "scenario": scenario,
                    "experimental_anchor": anchor,
                    "mapping_method": method,
                    "target_Hc_T": target_hc,
                    "B_star_T": float(mapped[0]),
                    "W95_min_T": float(mapped[1]),
                    "W95_max_T": float(mapped[2]),
                    "Qmax_mJ_per_cm3_inherited": qmax,
                })

                p0_hc = anchored[("P0", anchor)]
                p0_field = float(transform_field(
                    source_bstar, source_hc, p0_hc, method))
                source_eval = float(inverse_field(
                    p0_field, source_hc, target_hc, method))
                q_at_p0 = float(np.interp(source_eval, bs, q))
                regret_rows.append({
                    "temperature_K": temperature,
                    "scenario": scenario,
                    "experimental_anchor": anchor,
                    "mapping_method": method,
                    "P0_selected_field_T": p0_field,
                    "scenario_optimal_field_T": float(mapped[0]),
                    "delta_B_star_vs_P0_T": float(mapped[0] - p0_field),
                    "Q_at_P0_field_mJ_per_cm3": q_at_p0,
                    "Qmax_mJ_per_cm3_inherited": qmax,
                    "regret_percent": 100 * (qmax - q_at_p0) / qmax,
                })

        labels = [(method, scenario, anchor, hc)
                  for method in METHODS
                  for (scenario, anchor), hc in anchored.items()]
        grid = np.linspace(1.25, 1.65, 8001)
        regrets = []
        for method, _, _, target_hc in labels:
            source_eval = inverse_field(grid, source_hc, target_hc, method)
            q_eval = np.interp(source_eval, bs, q, left=np.nan, right=np.nan)
            regrets.append((qmax - q_eval) / qmax)
        regrets = np.asarray(regrets)
        valid = np.all(np.isfinite(regrets), axis=0)
        worst = np.nanmax(regrets, axis=0)
        valid_index = np.where(valid)[0]
        idx = valid_index[np.argmin(worst[valid])]
        robust_q = qmax * (1 - worst[idx])
        robust = {
            "temperature_K": temperature,
            "robust_field_T": float(grid[idx]),
            "worst_case_regret_percent": float(100 * worst[idx]),
            "worst_case_Q_mJ_per_cm3": float(robust_q),
            "scenario_order": [f"{m}/{s}/{a}" for m, s, a, _ in labels],
            "scenario_regret_percent": [float(100 * x)
                                         for x in regrets[:, idx]],
        }
        robust_rows.append(robust)

        for heat_load_uW in (0.05, 0.10, 0.50):
            for hold_minutes in (30, 60):
                required_mJ = heat_load_uW * 1e-3 * hold_minutes * 60
                design_rows.append({
                    "temperature_K": temperature,
                    "robust_field_T": robust["robust_field_T"],
                    "worst_case_regret_percent":
                        robust["worst_case_regret_percent"],
                    "worst_case_Q_mJ_per_cm3": robust_q,
                    "assumed_total_cold_load_uW": heat_load_uW,
                    "required_hold_minutes": hold_minutes,
                    "required_heat_absorption_mJ": required_mJ,
                    "conditional_magnetic_only_volume_cm3":
                        required_mJ / robust_q,
                    "conditional_hold_minutes_for_0p4_cm3":
                        robust_q * 0.4 / (heat_load_uW * 1e-3) / 60,
                })

    write_csv("interlayer_P0_P2_critical_fields.csv", hc_rows)
    write_csv("interlayer_P0_P2_metrics.csv", metric_rows)
    write_csv("interlayer_P0_P2_regret.csv", regret_rows)
    write_csv("robust_conditional_design_map.csv", design_rows)

    p2_90 = [row for row in regret_rows
             if row["scenario"] == "P2" and row["temperature_K"] == 0.09]
    payload = {
        "purpose": "public-data P0-P2 interlayer uncertainty and engineering propagation",
        "public_evidence": {
            "P0": "2D M3: long-range dipole plus bond-dependent exchange, no exchange between layers",
            "P1": {
                "source": "arXiv:2202.05242 source lines 2081-2083",
                "J_perp_K": JPERP_DFT_K,
                "model": "one isotropic vertical bond above and below",
                "analytic_Hc_shift_T": P1_HC_SHIFT_T,
            },
            "P2": {
                "source": "arXiv:2607.23173v1 Table 1 and lines 275-317",
                "J2_microeV": 2.8,
                "Ja_microeV": 0.0,
                "Jb_microeV": 7.9,
                "Delta_c": 0.08,
                "reported_out_of_plane_magnon_bandwidth_microeV": 4.4,
                "warning": "The authors call this a representative, non-unique parameter set.",
            },
        },
        "critical_field_scenarios": hc_rows,
        "experiment_anchored_Hc_scenarios": [
            {"scenario": scenario, "experimental_anchor": anchor,
             "critical_field_T": hc}
            for (scenario, anchor), hc in anchored.items()
        ],
        "P2_minus_P0_Hc_shift_T_range": [min(p2_shifts), max(p2_shifts)],
        "90mK_P2_minus_P0_delta_B_star_T_range": [
            min(row["delta_B_star_vs_P0_T"] for row in p2_90),
            max(row["delta_B_star_vs_P0_T"] for row in p2_90),
        ],
        "90mK_P2_regret_at_P0_field_percent_range": [
            min(row["regret_percent"] for row in p2_90),
            max(row["regret_percent"] for row in p2_90),
        ],
        "predeclared_thresholds": {
            "absolute_delta_B_star_GO_T": 0.05,
            "regret_GO_percent": 5.0,
        },
        "threshold_screen": {
            "P1_all_90mK_cases_above_delta_B_threshold": all(
                abs(row["delta_B_star_vs_P0_T"]) >= 0.05
                for row in regret_rows if row["scenario"] == "P1" and
                row["temperature_K"] == 0.09),
            "P2_all_90mK_cases_above_delta_B_threshold": all(
                abs(row["delta_B_star_vs_P0_T"]) >= 0.05 for row in p2_90),
            "P2_all_90mK_cases_above_regret_threshold": all(
                row["regret_percent"] >= 5.0 for row in p2_90),
        },
        "robust_fields": robust_rows,
        "design_map_rows": len(design_rows),
        "limitations": [
            "P1 is an explicit simple-bond interpretation of the published DFT scale, not a unique crystal model.",
            "P2 is the paper's representative interlayer set and is not uniquely fitted by experiment.",
            "The finite-temperature response is transported from a public XTRG curve by two similarity mappings; entropy redistribution and changes in Qmax are not calculated.",
            "Heat-load, hold-time and volume values are magnetic-only conditional scenarios; nuclear, lattice, heat-switch, addenda and thermal-contact loads are excluded.",
            "P2 is anchored to the measured 1.65--1.70 T critical field; only the calculated P0/P1 ablation shifts are transported into engineering coordinates.",
            "The direct experimental Hc anchor already contains all real interactions and remains the primary engineering decision evidence.",
        ],
    }
    (RESULTS / "interlayer_engineering_scenarios.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    x = np.arange(3)
    width = 0.35
    for offset, boundary, color in ((-width/2, "tin_foil", "#4c78a8"),
                                     (width/2, "sphere_vacuum", "#f58518")):
        vals = [scenarios[(name, boundary)] for name in ("P0", "P1", "P2")]
        axes[0].bar(x + offset, vals, width, label=boundary.replace("_", " "),
                    color=color)
    axes[0].set_xticks(x, ("P0", "P1", "P2"))
    axes[0].set_ylabel(r"$H_c$ (T)")
    axes[0].set_title("A. Model-scale $H_c$ (not operating field)")
    axes[0].legend(frameon=True, fontsize=8)

    temps_mk = [row["temperature_K"] * 1000 for row in robust_rows]
    axes[1].plot(temps_mk, [row["robust_field_T"] for row in robust_rows],
                 "o-", color="#54a24b", label=r"$B_{robust}$")
    ax2 = axes[1].twinx()
    ax2.plot(temps_mk,
             [row["worst_case_regret_percent"] for row in robust_rows],
             "s--", color="#e45756", label="worst regret")
    axes[1].set_xlabel("Temperature (mK)")
    axes[1].set_ylabel(r"$B_{robust}$ (T)")
    ax2.set_ylabel("Worst regret (%)")
    axes[1].set_title("B. Minimax operating point")
    lines = axes[1].lines + ax2.lines
    axes[1].legend(lines, [line.get_label() for line in lines],
                   frameon=True, fontsize=8)

    selected = [row for row in design_rows
                if row["assumed_total_cold_load_uW"] == 0.10 and
                row["required_hold_minutes"] == 60]
    axes[2].bar([row["temperature_K"] * 1000 for row in selected],
                [row["conditional_magnetic_only_volume_cm3"]
                 for row in selected], width=8, color="#72b7b2")
    axes[2].set_xlabel("Temperature (mK)")
    axes[2].set_ylabel(r"Conditional material volume (cm$^3$)")
    axes[2].set_title("C. 0.10 µW for 60 min")
    fig.suptitle("Interlayer uncertainty propagated to magnetic refrigeration",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    figure = FIGURES / "interlayer_engineering_scenarios.png"
    fig.savefig(figure, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "output": str(RESULTS / "interlayer_engineering_scenarios.json"),
        "figure": str(figure),
        "P1_Hc_shift_T": P1_HC_SHIFT_T,
        "P2_Hc_shift_T_range": payload["P2_minus_P0_Hc_shift_T_range"],
        "P2_delta_B_star_90mK_T_range":
            payload["90mK_P2_minus_P0_delta_B_star_T_range"],
        "P2_regret_90mK_percent_range":
            payload["90mK_P2_regret_at_P0_field_percent_range"],
        "robust_fields": robust_rows,
    }, indent=2))


if __name__ == "__main__":
    main()
