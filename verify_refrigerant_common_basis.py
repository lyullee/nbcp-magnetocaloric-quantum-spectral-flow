"""Verify the common-basis refrigerant comparison and its scope boundary."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
SUMMARY = RESULTS / "refrigerant_common_basis_summary.json"
OUTPUT = RESULTS / "refrigerant_common_basis_verification.json"


def close(value: float, target: float, atol: float = 1e-9) -> bool:
    return abs(float(value) - float(target)) <= atol


def main() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    with (RESULTS / "refrigerant_common_basis_comparison.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    by_material = {row["material"]: row for row in rows}
    ratios = summary["nbcp_ratios_vs_best_hydrate_comparator"]
    figures = [
        ROOT / "figures" / f"refrigerant_common_basis_comparison.{suffix}"
        for suffix in ("png", "pdf", "svg")
    ]
    checks = {
        "three_common_basis_materials": set(by_material) == {"NBCP", "FAA", "CPA"},
        "volumes_are_matched": summary["common_volume_range_cm3"] == [0.395, 0.421],
        "source_reported_entropy_values": all(
            close(by_material[m]["magnetic_entropy_density_mJ_K_cm3"], v)
            for m, v in {"NBCP": 56.0, "FAA": 53.0, "CPA": 42.0}.items()
        ),
        "source_reported_heat_values": all(
            close(by_material[m]["heat_absorption_mJ_cm3"], v)
            for m, v in {"NBCP": 19.0, "FAA": 3.0, "CPA": 3.1}.items()
        ),
        "source_reported_efficiency_values": all(
            close(by_material[m]["cycle_efficiency_percent"], v)
            for m, v in {"NBCP": 23.0, "FAA": 9.0, "CPA": 11.0}.items()
        ),
        "safe_time_trace_values": all(
            close(by_material[m]["below_1K_observed_span_min"], v, 1e-6)
            for m, v in {
                "NBCP": 95.95399,
                "FAA": 32.156500000040964,
                "CPA": 24.66499999997907,
            }.items()
        ),
        "nbcp_heat_ratio": close(ratios["heat_absorption_ratio_vs_best_hydrate"], 19.0 / 3.1),
        "nbcp_efficiency_ratio": close(ratios["efficiency_ratio_vs_best_hydrate"], 23.0 / 11.0),
        "nbcp_time_ratio": close(ratios["below_1K_span_ratio_vs_best_hydrate"], 95.95399 / 32.156500000040964),
        "ggg_is_context_only": (
            summary["separate_context_rows"][0]["material"] == "GGG"
            and "excluded" in summary["separate_context_rows"][0]["usage"]
        ),
        "ambiguous_curves_excluded": any(
            "sheet-a/b" in statement and "excluded" in statement
            for statement in summary["scope"]
        ),
        "not_a_device_prediction": any(
            "does not predict an ADR stage" in statement for statement in summary["scope"]
        ),
        "all_figures_exist": all(path.is_file() and path.stat().st_size > 0 for path in figures),
    }
    failed = [name for name, passed in checks.items() if not passed]
    payload = {"checks": checks, "failed": failed, "status": "PASS" if not failed else "FAIL"}
    OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
