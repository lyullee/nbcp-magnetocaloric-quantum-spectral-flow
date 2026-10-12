"""Independent identities and boundary checks for the transient ADR envelope."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
SOURCES = ROOT / "sources"


def close(a: float, b: float, rtol: float = 2e-11, atol: float = 2e-12) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def main() -> None:
    with (RESULTS / "adr_transient_cycle_envelope.csv").open(
            encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    with (RESULTS / "adr_transient_allowance_boundaries.csv").open(
            encoding="utf-8-sig") as handle:
        boundaries = list(csv.DictReader(handle))
    with (RESULTS / "adr_stage_load_envelope.csv").open(
            encoding="utf-8-sig") as handle:
        steady = list(csv.DictReader(handle))
    summary = json.loads(
        (RESULTS / "adr_transient_cycle_summary.json").read_text(encoding="utf-8"))
    provenance = json.loads(
        (SOURCES / "adr_transient_cycle_provenance_20261012.json").read_text(
            encoding="utf-8"))

    checks = {
        "row_count_matches_summary": len(rows) == summary["rows"] == 2688,
        "boundary_row_count_matches_summary":
            len(boundaries) == summary["allowance_boundary_rows"] == 35,
        "loss_components_sum": True,
        "hold_energy_identity": True,
        "gross_capacity_identity": True,
        "component_energy_identity": True,
        "mass_and_volume_identity": True,
        "mass_multiplier_identity": True,
        "duty_identity": True,
        "average_payload_identity": True,
        "cadence_identity": True,
        "zero_loss_reproduces_steady_envelope": True,
        "mass_cap_boundary_identity": True,
        "bartlett_source_crosscheck": True,
        "canonical_mass_multipliers": True,
    }

    # Density is recovered independently from mass / volume for each row and
    # required only to check internal consistency across the generated set.
    density_values = []
    for row in rows:
        f_add = float(row["addenda_fraction_of_gross_capacity"])
        f_ramp = float(row["ramp_fraction_of_gross_capacity"])
        f_switch = float(row["switch_fraction_of_gross_capacity"])
        f_total = float(row["total_transient_fraction"])
        usable = float(row["usable_capacity_fraction"])
        payload = float(row["payload_load_uW"])
        parasitic = float(row["steady_parasitic_load_uW"])
        hold_s = 60.0 * float(row["hold_minutes"])
        hold_energy = float(row["hold_energy_mJ"])
        gross = float(row["gross_magnetic_capacity_mJ"])
        q = float(row["worst_case_Q_mJ_per_cm3"])
        volume = float(row["conditional_NBCP_volume_cm3"])
        mass = float(row["conditional_NBCP_mass_g"])
        recycle_s = float(row["recycle_time_s"])
        duty = float(row["single_stage_duty_fraction"])
        margin = float(row["two_stage_cadence_margin_s"])
        density_values.append(mass / volume)

        checks["loss_components_sum"] &= close(f_add + f_ramp + f_switch, f_total)
        checks["loss_components_sum"] &= 0 <= f_total < 1 and close(usable, 1 - f_total)
        checks["hold_energy_identity"] &= close(
            hold_energy, (payload + parasitic) * 1e-3 * hold_s)
        checks["gross_capacity_identity"] &= close(gross, hold_energy / usable)
        checks["component_energy_identity"] &= close(
            float(row["addenda_allowance_mJ"]), f_add * gross)
        checks["component_energy_identity"] &= close(
            float(row["ramp_loss_allowance_mJ"]), f_ramp * gross)
        checks["component_energy_identity"] &= close(
            float(row["switch_loss_allowance_mJ"]), f_switch * gross)
        checks["mass_and_volume_identity"] &= close(volume, gross / q)
        checks["mass_multiplier_identity"] &= close(
            float(row["mass_multiplier_vs_zero_transient"]), 1 / usable)
        checks["duty_identity"] &= close(duty, hold_s / (hold_s + recycle_s))
        checks["average_payload_identity"] &= close(
            float(row["average_payload_if_no_cooling_during_recycle_uW"]),
            payload * duty)
        checks["cadence_identity"] &= close(margin, hold_s - recycle_s)
        checks["cadence_identity"] &= (
            int(row["two_stage_cadence_feasible_without_handover_margin"])
            == int(margin >= 0))

    checks["mass_and_volume_identity"] &= (
        max(density_values) - min(density_values) < 1e-11)

    steady_key = {
        (float(r["temperature_K"]), float(r["payload_load_uW"]),
         float(r["steady_parasitic_load_uW"]), float(r["hold_minutes"])):
            float(r["conditional_NBCP_mass_g"])
        for r in steady
    }
    zero_rows = [
        r for r in rows
        if r["loss_scenario"] == "zero_loss_bound"
        and abs(float(r["recycle_time_s"]) - 82.0) < 1e-12]
    for row in zero_rows:
        key = (float(row["temperature_K"]), float(row["payload_load_uW"]),
               float(row["steady_parasitic_load_uW"]), float(row["hold_minutes"]))
        checks["zero_loss_reproduces_steady_envelope"] &= (
            key in steady_key and close(
                float(row["conditional_NBCP_mass_g"]), steady_key[key]))

    for row in boundaries:
        gross = float(row["gross_magnetic_capacity_mJ"])
        hold = float(row["hold_energy_mJ"])
        allowance = float(row["maximum_transient_energy_allowance_mJ"])
        feasible = int(row["mass_cap_feasible_before_transient_losses"])
        checks["mass_cap_boundary_identity"] &= close(allowance, gross - hold)
        checks["mass_cap_boundary_identity"] &= feasible == int(allowance >= 0)
        if feasible:
            checks["mass_cap_boundary_identity"] &= close(
                float(row["maximum_transient_fraction_of_gross_capacity"]),
                allowance / gross)
        else:
            checks["mass_cap_boundary_identity"] &= (
                row["maximum_transient_fraction_of_gross_capacity"] == "")

    bartlett = next(
        a for a in provenance["public_anchors"]
        if a["id"] == "bartlett_2015_fast_miniature_adr")
    arithmetic = bartlett["hold_time_s"] / (
        bartlett["hold_time_s"] + bartlett["complete_recycle_time_s"])
    checks["bartlett_source_crosscheck"] &= close(
        arithmetic, bartlett["arithmetic_hold_over_hold_plus_recycle"])
    checks["bartlett_source_crosscheck"] &= abs(
        arithmetic - bartlett["paper_reported_duty_fraction"]) < 0.03
    checks["bartlett_source_crosscheck"] &= (
        bartlett["heat_switch_share_of_parasitic_load_fraction"] == 0.82
        and bartlett["conductivity_reduction_factor_within_5s"] == 1000.0)

    expected = {
        "zero_loss_bound": 1.0,
        "ten_percent_stress": 1 / 0.9,
        "quarter_capacity_stress": 1 / 0.75,
        "half_capacity_stress": 2.0,
    }
    for item in summary["canonical_90mK_0p10uW_60min_zero_parasitic"]:
        checks["canonical_mass_multipliers"] &= close(
            float(item["mass_multiplier_vs_zero_transient"]),
            expected[item["loss_scenario"]])

    passed = all(checks.values())
    report = {
        "status": "PASS" if passed else "FAIL",
        "checks": checks,
        "rows_checked": len(rows),
        "boundary_rows_checked": len(boundaries),
        "max_density_spread_g_cm3": max(density_values) - min(density_values),
    }
    (RESULTS / "adr_transient_cycle_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
