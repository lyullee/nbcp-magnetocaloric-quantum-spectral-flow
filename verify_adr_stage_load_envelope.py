"""Regression checks for the steady-hold ADR stage load envelope."""
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return abs(a - b) <= tol


def main() -> None:
    with (RESULTS / "adr_stage_load_envelope.csv").open(
            encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    summary = json.loads(
        (RESULTS / "adr_stage_load_envelope_summary.json").read_text(
            encoding="utf-8"))
    with (RESULTS / "refrigerator_design_chart.csv").open(
            encoding="utf-8-sig") as handle:
        base = list(csv.DictReader(handle))

    checks: dict[str, bool] = {}
    checks["grid_complete"] = (
        len(rows) == 4 * 3 * 10 * 2 and
        {float(r["temperature_K"]) for r in rows} == {0.09, 0.10, 0.12, 0.15} and
        {float(r["payload_load_uW"]) for r in rows} == {0.05, 0.10, 0.50} and
        {float(r["hold_minutes"]) for r in rows} == {30.0, 60.0})
    checks["energy_identity"] = all(close(
        float(r["required_heat_absorption_mJ"]),
        float(r["total_steady_load_uW"]) * 1e-3 *
        float(r["hold_minutes"]) * 60.0) for r in rows)
    checks["load_decomposition"] = all(close(
        float(r["total_steady_load_uW"]),
        float(r["payload_load_uW"]) + float(r["steady_parasitic_load_uW"]))
        for r in rows)
    checks["universal_mass_multiplier"] = all(close(
        float(r["mass_multiplier_vs_payload_only"]),
        1.0 + float(r["steady_parasitic_load_uW"]) /
        float(r["payload_load_uW"])) for r in rows)

    stage_base = next(r for r in rows
                      if close(float(r["temperature_K"]), 0.09)
                      and close(float(r["payload_load_uW"]), 0.10)
                      and close(float(r["steady_parasitic_load_uW"]), 0.0)
                      and close(float(r["hold_minutes"]), 60.0))
    chart_base = next(r for r in base
                      if close(float(r["temperature_K"]), 0.09)
                      and close(float(r["assumed_total_cold_load_uW"]), 0.10)
                      and close(float(r["required_hold_minutes"]), 60.0))
    checks["zero_parasitic_matches_existing_chart"] = close(
        float(stage_base["conditional_NBCP_mass_g"]),
        float(chart_base["conditional_material_mass_g"]))

    canonical = summary["canonical_90mK_payload_0p10uW_hold_60min"]
    expected = {0.0: 3.9582202162443876,
                0.20: 11.874660648733164,
                0.97: 42.352956313814936,
                13.14: 524.068356630757}
    checks["canonical_values"] = all(
        close(float(row["conditional_NBCP_mass_g"]),
              expected[float(row["steady_parasitic_load_uW"])])
        for row in canonical)
    checks["tophat_source_identity"] = (
        summary["source_crosscheck"]["relative_difference"] < 0.002)
    checks["context_not_prediction_scope"] = (
        len(summary["published_context_anchors"]) == 3 and
        all("context anchor only" in a["usage"]
            for a in summary["published_context_anchors"]) and
        any("not NBCP device assumptions" in line for line in summary["scope"]))

    failed = [name for name, passed in checks.items() if not passed]
    payload = {"checks": checks, "failed": failed,
               "status": "PASS" if not failed else "FAIL"}
    (RESULTS / "adr_stage_load_envelope_verification.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
