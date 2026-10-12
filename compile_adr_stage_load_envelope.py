"""Separate payload and steady parasitic loads in the NBCP hold budget.

Published ADR systems provide contextual load anchors only.  Their values are
not transferred as an NBCP prior or a device prediction because their stages,
temperatures, interfaces and refrigerants differ.  The calculation exposes the
exact linear sensitivity of the existing magnetic-only lower bound:

    Q_required = (P_payload + P_parasitic) * t_hold
    V_NBCP = Q_required / q_worst(T, B_robust)

Transient addenda, regeneration heat, ramp loss, packing and contact resistance
remain outside this steady-hold envelope.
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
PUBLIC_SAMPLE_MASS_G = 1.65
PARASITIC_SCENARIOS_UW = (0.0, 0.05, 0.10, 0.20, 0.50, 0.97, 1.0, 2.0, 5.0, 13.14)
PAYLOAD_LOADS_UW = (0.05, 0.10, 0.50)
HOLD_TIMES_MIN = (30, 60)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select(rows: list[dict], **criteria: float) -> dict:
    for row in rows:
        if all(abs(float(row[key]) - value) < 1e-12
               for key, value in criteria.items()):
            return row
    raise KeyError(criteria)


def main() -> None:
    with (RESULTS / "thermodynamic_shape_robust_decisions.csv").open(
            encoding="utf-8-sig") as handle:
        robust = list(csv.DictReader(handle))
    design_summary = json.loads(
        (RESULTS / "refrigerator_design_chart_summary.json").read_text(
            encoding="utf-8"))
    provenance = json.loads(
        (SOURCES / "adr_device_anchor_provenance_20261012.json").read_text(
            encoding="utf-8"))
    density = float(design_summary["density_derivation"]["density_g_cm3"])
    sample_volume = PUBLIC_SAMPLE_MASS_G / density

    rows: list[dict] = []
    for decision in robust:
        temperature = float(decision["temperature_K"])
        field = float(decision["robust_field_T"])
        q_worst = float(decision["worst_case_Q_at_robust_field_mJ_per_cm3"])
        for payload in PAYLOAD_LOADS_UW:
            for parasitic in PARASITIC_SCENARIOS_UW:
                total = payload + parasitic
                for hold in HOLD_TIMES_MIN:
                    heat_mj = total * 1e-3 * hold * 60.0
                    volume = heat_mj / q_worst
                    rows.append({
                        "temperature_K": temperature,
                        "robust_field_T": field,
                        "worst_case_Q_at_robust_field_mJ_per_cm3": q_worst,
                        "payload_load_uW": payload,
                        "steady_parasitic_load_uW": parasitic,
                        "total_steady_load_uW": total,
                        "hold_minutes": hold,
                        "required_heat_absorption_mJ": heat_mj,
                        "conditional_NBCP_volume_cm3": volume,
                        "conditional_NBCP_mass_g": volume * density,
                        "mass_multiplier_vs_payload_only": total / payload,
                        "public_1p65g_conditional_hold_minutes":
                            q_worst * sample_volume / (total * 1e-3) / 60.0,
                    })

    anchors = provenance["anchors"]
    anchor_rows = []
    for anchor in anchors:
        anchor_rows.append({
            "id": anchor["id"],
            "system": anchor["system"],
            "refrigerant": anchor["refrigerant"],
            "operating_temperature_K": anchor["operating_temperature_K"],
            "parasitic_load_uW": anchor["parasitic_load_uW"],
            "reported_hold_or_observation_s": anchor["reported_hold_or_observation_s"],
            "source_kind": anchor["source_kind"],
            "doi_or_url": anchor.get("doi", anchor["url"]),
            "usage": "context anchor only; not transferred as an NBCP prediction",
        })

    canonical = []
    for parasitic in (0.0, 0.20, 0.97, 13.14):
        row = select(rows, temperature_K=0.09, payload_load_uW=0.10,
                     steady_parasitic_load_uW=parasitic, hold_minutes=60.0)
        canonical.append({key: row[key] for key in (
            "payload_load_uW", "steady_parasitic_load_uW",
            "total_steady_load_uW", "hold_minutes",
            "conditional_NBCP_mass_g", "mass_multiplier_vs_payload_only",
            "public_1p65g_conditional_hold_minutes")})

    tophat = next(a for a in anchors if a["id"] == "tophat_faa_adr")
    tophat_identity_mj = (tophat["parasitic_load_uW"] * 1e-3 *
                          tophat["reported_hold_or_observation_s"])
    summary = {
        "status": "steady-hold stage-load sensitivity envelope",
        "model": {
            "required_heat_mJ": "(payload_uW + steady_parasitic_uW) * 1e-3 * hold_s",
            "required_volume_cm3": "required_heat_mJ / worst_case_Q_mJ_per_cm3",
            "mass_multiplier": "1 + steady_parasitic_uW / payload_uW",
        },
        "rows": len(rows),
        "temperature_K": sorted({float(r["temperature_K"]) for r in rows}),
        "payload_load_uW": list(PAYLOAD_LOADS_UW),
        "steady_parasitic_load_uW": list(PARASITIC_SCENARIOS_UW),
        "hold_minutes": list(HOLD_TIMES_MIN),
        "canonical_90mK_payload_0p10uW_hold_60min": canonical,
        "published_context_anchors": anchor_rows,
        "source_crosscheck": {
            "tophat_reported_capacity_mJ": tophat["reported_cooling_capacity_mJ"],
            "tophat_load_times_hold_mJ": tophat_identity_mj,
            "relative_difference": abs(
                tophat_identity_mj - tophat["reported_cooling_capacity_mJ"]
            ) / tophat["reported_cooling_capacity_mJ"],
        },
        "scope": provenance["scope"] + [
            "The envelope covers an isothermal steady hold only.",
            "Addenda and ramp losses require transient heat capacities and a cycle model and are not represented by a fabricated number.",
            "The worst-case magnetic budget comes from the declared NBCP thermodynamic shape envelope.",
        ],
    }

    write_csv(RESULTS / "adr_stage_load_envelope.csv", rows)
    write_csv(RESULTS / "adr_literature_system_anchors.csv", anchor_rows)
    (RESULTS / "adr_stage_load_envelope_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.25))
    p_grid = np.linspace(0.0, 1.2, 241)
    decision90 = next(d for d in robust if abs(float(d["temperature_K"]) - 0.09) < 1e-12)
    q90 = float(decision90["worst_case_Q_at_robust_field_mJ_per_cm3"])
    for payload in PAYLOAD_LOADS_UW:
        heat = (payload + p_grid) * 1e-3 * 3600.0
        mass = heat / q90 * density
        axes[0].plot(p_grid, mass, label=rf"payload {payload:.2f} $\mu$W")
        hold = q90 * sample_volume / ((payload + p_grid) * 1e-3) / 60.0
        axes[1].plot(p_grid, hold, label=rf"payload {payload:.2f} $\mu$W")
    for value, label in ((0.20, "0.20"), (0.97, "0.97")):
        axes[0].axvline(value, color="#777777", ls="--", lw=1)
        axes[1].axvline(value, color="#777777", ls="--", lw=1)
        axes[0].text(value, axes[0].get_ylim()[1] * 0.96, label,
                     ha="center", va="top", fontsize=8, color="#555555")
    axes[0].set(xlabel=r"Steady parasitic load ($\mu$W)",
                ylabel="NBCP mass for 60 min (g)",
                title="90 mK steady-hold inventory")
    axes[1].set(xlabel=r"Steady parasitic load ($\mu$W)",
                ylabel="Conditional hold time (min)",
                title="Public 1.65 g sample equivalent", ylim=(0, 55))
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)

    colors = ["#4c78a8", "#f58518", "#e45756"]
    for anchor, color in zip(anchors, colors):
        axes[2].scatter(1000 * anchor["operating_temperature_K"],
                        anchor["parasitic_load_uW"], s=70, color=color)
        axes[2].annotate(anchor["id"].split("_")[0],
                         (1000 * anchor["operating_temperature_K"],
                          anchor["parasitic_load_uW"]),
                         xytext=(5, 4), textcoords="offset points", fontsize=8)
    axes[2].set_yscale("log")
    axes[2].set(xlabel="Published operating temperature (mK)",
                ylabel=r"Published parasitic load ($\mu$W)",
                title="System anchors (context only)")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"adr_stage_load_envelope.{suffix}", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
