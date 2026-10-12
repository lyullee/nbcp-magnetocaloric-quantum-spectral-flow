"""Build a parametric transient and cycle envelope for the NBCP ADR stage.

The magnetic capacity is the already verified worst-case public-data budget.
No transient loss is known for a real NBCP package, so addenda, ramp and
switching terms are declared fractions of *gross* magnetic capacity.  This
makes the remaining device measurements and mass margins explicit without
fabricating an NBCP prototype prediction.

For total transient fraction f < 1,

    E_hold = (P_payload + P_parasitic) t_hold
    E_gross = E_hold / (1 - f)
    m = rho E_gross / q_worst

Cycle time is kept separate:

    duty_single = t_hold / (t_hold + t_recycle).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCES = ROOT / "sources"
PAYLOAD_LOADS_UW = (0.05, 0.10, 0.50)
PARASITIC_LOADS_UW = (0.0, 0.05, 0.10, 0.20, 0.50, 0.97, 1.0)
HOLD_TIMES_MIN = (30.0, 60.0)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with (RESULTS / "thermodynamic_shape_robust_decisions.csv").open(
            encoding="utf-8-sig") as handle:
        decisions = list(csv.DictReader(handle))
    design = json.loads(
        (RESULTS / "refrigerator_design_chart_summary.json").read_text(
            encoding="utf-8"))
    provenance = json.loads(
        (SOURCES / "adr_transient_cycle_provenance_20261012.json").read_text(
            encoding="utf-8"))
    density = float(design["density_derivation"]["density_g_cm3"])
    loss_scenarios = provenance["loss_scenarios"]
    recycle_coordinates = provenance["recycle_time_coordinates_s"]

    rows: list[dict] = []
    for decision in decisions:
        temperature = float(decision["temperature_K"])
        field = float(decision["robust_field_T"])
        q_worst = float(decision["worst_case_Q_at_robust_field_mJ_per_cm3"])
        for payload in PAYLOAD_LOADS_UW:
            for parasitic in PARASITIC_LOADS_UW:
                total_power = payload + parasitic
                for hold_min in HOLD_TIMES_MIN:
                    hold_s = 60.0 * hold_min
                    hold_energy = total_power * 1e-3 * hold_s
                    for scenario in loss_scenarios:
                        f_add = float(scenario["addenda_fraction_of_gross_capacity"])
                        f_ramp = float(scenario["ramp_fraction_of_gross_capacity"])
                        f_switch = float(scenario["switch_fraction_of_gross_capacity"])
                        f_total = f_add + f_ramp + f_switch
                        usable_fraction = 1.0 - f_total
                        gross_energy = hold_energy / usable_fraction
                        volume = gross_energy / q_worst
                        mass = density * volume
                        for recycle in recycle_coordinates:
                            recycle_s = float(recycle["value"])
                            duty = hold_s / (hold_s + recycle_s)
                            margin = hold_s - recycle_s
                            rows.append({
                                "temperature_K": temperature,
                                "robust_field_T": field,
                                "worst_case_Q_mJ_per_cm3": q_worst,
                                "payload_load_uW": payload,
                                "steady_parasitic_load_uW": parasitic,
                                "total_steady_load_uW": total_power,
                                "hold_minutes": hold_min,
                                "loss_scenario": scenario["id"],
                                "addenda_fraction_of_gross_capacity": f_add,
                                "ramp_fraction_of_gross_capacity": f_ramp,
                                "switch_fraction_of_gross_capacity": f_switch,
                                "total_transient_fraction": f_total,
                                "usable_capacity_fraction": usable_fraction,
                                "hold_energy_mJ": hold_energy,
                                "gross_magnetic_capacity_mJ": gross_energy,
                                "addenda_allowance_mJ": f_add * gross_energy,
                                "ramp_loss_allowance_mJ": f_ramp * gross_energy,
                                "switch_loss_allowance_mJ": f_switch * gross_energy,
                                "conditional_NBCP_volume_cm3": volume,
                                "conditional_NBCP_mass_g": mass,
                                "mass_multiplier_vs_zero_transient": 1.0 / usable_fraction,
                                "recycle_time_s": recycle_s,
                                "recycle_coordinate_kind": recycle["kind"],
                                "single_stage_duty_fraction": duty,
                                "average_payload_if_no_cooling_during_recycle_uW": payload * duty,
                                "two_stage_cadence_margin_s": margin,
                                "two_stage_cadence_feasible_without_handover_margin": int(margin >= 0),
                            })

    # Maximum loss energy and fraction that fit selected mass caps for the
    # canonical 90 mK, 0.10 uW payload, 60 min hold.
    decision90 = next(
        row for row in decisions if abs(float(row["temperature_K"]) - 0.09) < 1e-12)
    q90 = float(decision90["worst_case_Q_at_robust_field_mJ_per_cm3"])
    allowance_rows: list[dict] = []
    for parasitic in PARASITIC_LOADS_UW:
        hold_energy = (0.10 + parasitic) * 1e-3 * 3600.0
        for mass_cap in provenance["mass_cap_coordinates_g"]:
            gross_capacity = q90 * (float(mass_cap) / density)
            allowance = gross_capacity - hold_energy
            feasible = allowance >= 0.0
            allowance_rows.append({
                "temperature_K": 0.09,
                "payload_load_uW": 0.10,
                "steady_parasitic_load_uW": parasitic,
                "hold_minutes": 60.0,
                "mass_cap_g": mass_cap,
                "gross_magnetic_capacity_mJ": gross_capacity,
                "hold_energy_mJ": hold_energy,
                "maximum_transient_energy_allowance_mJ": allowance,
                "maximum_transient_fraction_of_gross_capacity":
                    (allowance / gross_capacity if feasible else ""),
                "mass_cap_feasible_before_transient_losses": int(feasible),
            })

    write_csv(RESULTS / "adr_transient_cycle_envelope.csv", rows)
    write_csv(RESULTS / "adr_transient_allowance_boundaries.csv", allowance_rows)

    canonical = []
    for scenario in loss_scenarios:
        row = next(
            r for r in rows
            if abs(float(r["temperature_K"]) - 0.09) < 1e-12
            and abs(float(r["payload_load_uW"]) - 0.10) < 1e-12
            and abs(float(r["steady_parasitic_load_uW"])) < 1e-12
            and abs(float(r["hold_minutes"]) - 60.0) < 1e-12
            and r["loss_scenario"] == scenario["id"]
            and abs(float(r["recycle_time_s"]) - 82.0) < 1e-12)
        canonical.append({key: row[key] for key in (
            "loss_scenario", "total_transient_fraction",
            "gross_magnetic_capacity_mJ", "conditional_NBCP_mass_g",
            "mass_multiplier_vs_zero_transient", "single_stage_duty_fraction")})

    recycle_duties = []
    for coord in recycle_coordinates:
        recycle_s = float(coord["value"])
        recycle_duties.append({
            "recycle_time_s": recycle_s,
            "kind": coord["kind"],
            "single_stage_duty_for_60min_hold": 3600.0 / (3600.0 + recycle_s),
            "two_stage_cadence_margin_for_60min_hold_s": 3600.0 - recycle_s,
        })

    bartlett = next(
        a for a in provenance["public_anchors"]
        if a["id"] == "bartlett_2015_fast_miniature_adr")
    summary = {
        "status": "parametric transient-loss and cycle-timing envelope",
        "equations": {
            "hold_energy_mJ": "(payload_uW + parasitic_uW) * 1e-3 * hold_s",
            "gross_capacity_mJ": "hold_energy_mJ / (1 - total_transient_fraction)",
            "required_mass_g": "density_g_cm3 * gross_capacity_mJ / q_worst_mJ_cm3",
            "single_stage_duty": "hold_s / (hold_s + recycle_s)",
            "mass_cap_loss_fraction": "1 - hold_energy_mJ / gross_capacity_at_mass_cap_mJ",
        },
        "rows": len(rows),
        "allowance_boundary_rows": len(allowance_rows),
        "temperature_K": sorted({float(r["temperature_K"]) for r in rows}),
        "payload_load_uW": list(PAYLOAD_LOADS_UW),
        "steady_parasitic_load_uW": list(PARASITIC_LOADS_UW),
        "hold_minutes": list(HOLD_TIMES_MIN),
        "loss_scenarios": loss_scenarios,
        "recycle_time_coordinates_s": recycle_coordinates,
        "canonical_90mK_0p10uW_60min_zero_parasitic": canonical,
        "canonical_60min_cycle_duties": recycle_duties,
        "bartlett_public_crosscheck": {
            "reported_duty_fraction": bartlett["paper_reported_duty_fraction"],
            "arithmetic_hold_over_hold_plus_recycle":
                bartlett["hold_time_s"] /
                (bartlett["hold_time_s"] + bartlett["complete_recycle_time_s"]),
            "absolute_difference": abs(
                bartlett["paper_reported_duty_fraction"] -
                bartlett["hold_time_s"] /
                (bartlett["hold_time_s"] + bartlett["complete_recycle_time_s"])),
            "interpretation": "The small difference is consistent with unmodelled switchover/handover time; the reported value is not overwritten.",
        },
        "public_anchors": provenance["public_anchors"],
        "required_future_measurements": provenance["required_future_measurements"],
        "scope": provenance["model_boundary"],
    }
    (RESULTS / "adr_transient_cycle_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    # Figure: physical design margins without claiming an NBCP prototype.
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.0))
    ax = axes[0, 0]
    f_grid = np.linspace(0.0, 0.60, 241)
    for parasitic, color in zip((0.0, 0.20, 0.97), ("#4c78a8", "#f58518", "#e45756")):
        hold_energy = (0.10 + parasitic) * 1e-3 * 3600.0
        mass = density * hold_energy / ((1.0 - f_grid) * q90)
        ax.plot(100 * f_grid, mass, color=color,
                label=rf"$P_{{par}}$={parasitic:.2f} $\mu$W")
    ax.set(xlabel="Transient allowance (% of gross capacity)",
           ylabel="Required NBCP mass (g)",
           title="(a) 90 mK, 0.10 µW payload, 60 min")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    recycle_min = np.linspace(0.0, 60.0, 301)
    for hold_min, color in zip(HOLD_TIMES_MIN, ("#54a24b", "#4c78a8")):
        duty = hold_min / (hold_min + recycle_min)
        ax.plot(recycle_min, 100 * duty, color=color,
                label=f"hold {hold_min:.0f} min")
    for x, label in ((82.0 / 60.0, "82 s measured"), (6.5, "6.5 min target midpoint")):
        ax.axvline(x, color="#777777", ls="--", lw=1)
        ax.text(x + 0.4, 52, label, rotation=90, va="bottom", fontsize=8)
    ax.set(xlabel="Recycle time (min)", ylabel="Single-stage duty (%)",
           title="(b) Hold–recycle timing trade-off", ylim=(45, 101))
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    mass_caps = [float(x) for x in provenance["mass_cap_coordinates_g"]]
    heat = np.full((len(mass_caps), len(PARASITIC_LOADS_UW)), np.nan)
    for i, mass_cap in enumerate(mass_caps):
        for j, parasitic in enumerate(PARASITIC_LOADS_UW):
            item = next(
                r for r in allowance_rows
                if abs(float(r["mass_cap_g"]) - mass_cap) < 1e-12
                and abs(float(r["steady_parasitic_load_uW"]) - parasitic) < 1e-12)
            if item["mass_cap_feasible_before_transient_losses"]:
                heat[i, j] = 100 * float(item["maximum_transient_fraction_of_gross_capacity"])
    image = ax.imshow(heat, aspect="auto", origin="lower", vmin=0, vmax=100,
                      cmap="viridis")
    ax.set_xticks(range(len(PARASITIC_LOADS_UW)),
                  [f"{x:g}" for x in PARASITIC_LOADS_UW])
    ax.set_yticks(range(len(mass_caps)), [f"{x:g}" for x in mass_caps])
    ax.set(xlabel="Steady parasitic load (µW)", ylabel="Mass cap (g)",
           title="(c) Maximum transient fraction within mass cap")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            label = "—" if np.isnan(heat[i, j]) else f"{heat[i, j]:.0f}%"
            ax.text(j, i, label, ha="center", va="center", fontsize=7,
                    color="white" if not np.isnan(heat[i, j]) and heat[i, j] < 55 else "black")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Allowed fraction (%)")

    ax = axes[1, 1]
    ids = [s["id"].replace("_stress", "").replace("_", "\n") for s in loss_scenarios]
    gross = []
    hold_parts = []
    add_parts = []
    ramp_parts = []
    switch_parts = []
    hold_energy = 0.10 * 1e-3 * 3600.0
    for s in loss_scenarios:
        f_total = float(s["total_transient_fraction"])
        energy = hold_energy / (1.0 - f_total)
        gross.append(energy)
        hold_parts.append(hold_energy)
        add_parts.append(float(s["addenda_fraction_of_gross_capacity"]) * energy)
        ramp_parts.append(float(s["ramp_fraction_of_gross_capacity"]) * energy)
        switch_parts.append(float(s["switch_fraction_of_gross_capacity"]) * energy)
    x = np.arange(len(ids))
    ax.bar(x, hold_parts, label="60 min hold", color="#4c78a8")
    bottom = np.array(hold_parts)
    ax.bar(x, add_parts, bottom=bottom, label="addenda allowance", color="#f2cf5b")
    bottom += np.array(add_parts)
    ax.bar(x, ramp_parts, bottom=bottom, label="ramp allowance", color="#e45756")
    bottom += np.array(ramp_parts)
    ax.bar(x, switch_parts, bottom=bottom, label="switch allowance", color="#72b7b2")
    ax.set_xticks(x, ids, fontsize=8)
    ax.set(ylabel="Gross magnetic capacity (mJ)",
           title="(d) Declared loss allocation (zero parasitic)")
    ax.legend(fontsize=8, ncol=2)

    for ax in axes.flat:
        ax.grid(alpha=0.22)
    fig.suptitle("Parametric ADR transient and cycle envelope", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"adr_transient_cycle_envelope.{suffix}", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
