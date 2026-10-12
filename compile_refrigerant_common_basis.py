"""Compile a definition-matched public refrigerant comparison.

The calculation merges source-reported volumetric metrics with independently
labelled Extended Data Fig. 6d time traces. It deliberately excludes the
ambiguous sheet-a/b material labels and keeps the 50 g GGG experiment outside
the approximately 0.4 cm3 common-volume ranking.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
PROVENANCE = ROOT / "sources" / "refrigerant_common_basis_provenance_20261012.json"
TIME_AUDIT = ROOT / "results" / "refrigerant_experimental_benchmark.json"
OUT_CSV = ROOT / "results" / "refrigerant_common_basis_comparison.csv"
OUT_JSON = ROOT / "results" / "refrigerant_common_basis_summary.json"
FIGURE_BASE = ROOT / "figures" / "refrigerant_common_basis_comparison"


def main() -> None:
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    time_audit = json.loads(TIME_AUDIT.read_text(encoding="utf-8"))
    time_rows = {
        row["material"]: row
        for row in time_audit["time_trace_rows_safe_for_comparison"]
    }

    rows: list[dict[str, object]] = []
    for source_row in provenance["common_basis_rows"]:
        material = source_row["material"]
        trace = time_rows[material]
        row = dict(source_row)
        row.update(
            {
                "sample_mass_to_volume_g_cm3": source_row["sample_mass_g"]
                / source_row["sample_volume_cm3"],
                "observed_minimum_in_sheet_d_K": trace["minimum_temperature_K"],
                "below_1K_observed_span_min": trace["below_1K_observed_span_min"],
                "time_trace_points": trace["time_trace_points"],
                "time_trace_source_curve": trace["source_curve"],
            }
        )
        rows.append(row)

    conventional = [row for row in rows if row["material"] != "NBCP"]
    nbcp = next(row for row in rows if row["material"] == "NBCP")
    comparison_metrics = {
        "magnetic_entropy_density_mJ_K_cm3": "entropy_density_ratio_vs_best_hydrate",
        "heat_absorption_mJ_cm3": "heat_absorption_ratio_vs_best_hydrate",
        "cycle_efficiency_percent": "efficiency_ratio_vs_best_hydrate",
        "below_1K_observed_span_min": "below_1K_span_ratio_vs_best_hydrate",
    }
    ratios: dict[str, float] = {}
    for metric, label in comparison_metrics.items():
        ratios[label] = float(nbcp[metric]) / max(float(row[metric]) for row in conventional)

    fields = list(rows[0].keys())
    OUT_CSV.parent.mkdir(exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "status": "definition-matched public refrigerant comparison",
        "common_basis_material_count": len(rows),
        "common_basis_materials": [row["material"] for row in rows],
        "common_volume_range_cm3": [
            min(float(row["sample_volume_cm3"]) for row in rows),
            max(float(row["sample_volume_cm3"]) for row in rows),
        ],
        "rows": rows,
        "nbcp_ratios_vs_best_hydrate_comparator": ratios,
        "separate_context_rows": provenance["separate_context_rows"],
        "source": provenance["primary_source"],
        "metric_definitions": provenance["metric_definitions"],
        "scope": provenance["scope"],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    materials = [str(row["material"]) for row in rows]
    colors = ["#5B2A86", "#8FA3AD", "#B7C5CC"]
    panels = [
        ("magnetic_entropy_density_mJ_K_cm3", "Magnetic entropy density", "mJ K$^{-1}$ cm$^{-3}$"),
        ("heat_absorption_mJ_cm3", "Volumetric heat absorption", "mJ cm$^{-3}$"),
        ("cycle_efficiency_percent", "Source-defined cycle efficiency", "%"),
        ("below_1K_observed_span_min", "Observed time below 1 K", "min"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.8), constrained_layout=True)
    for ax, (metric, title, ylabel) in zip(axes.flat, panels):
        values = np.array([float(row[metric]) for row in rows])
        bars = ax.bar(materials, values, color=colors, edgecolor="#333333", linewidth=0.6)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        upper = max(values) * 1.22
        ax.set_ylim(0, upper)
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + upper * 0.025,
                f"{value:.1f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )
    fig.suptitle(
        "Public common-basis comparison: NBCP and volume-matched hydrate salts",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "Entropy, heat absorption and efficiency are source-reported; time spans are derived from independently labelled sheet-d traces.",
        ha="center",
        fontsize=8.5,
    )
    FIGURE_BASE.parent.mkdir(exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURE_BASE.with_suffix(f".{suffix}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
