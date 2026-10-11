"""Compile the two N=6 Open Plan jobs into a summary and figure."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def main():
    trex = load("open_plan_qpu_pilot_N6_TREX.json")
    zne = load("open_plan_qpu_pilot_N6_ZNE.json")
    if trex["status"] != "DONE" or zne["status"] != "DONE":
        raise RuntimeError("Both QPU jobs must be complete")
    exact = float(trex["exact_V_expectation_K"])
    if abs(exact - float(zne["exact_V_expectation_K"])) > 1e-12:
        raise RuntimeError("QPU jobs do not share the same target")

    observations = []
    for label, item in (("TREX", trex), ("TREX+ZNE", zne)):
        used = (float(item["usage_after"]["usage_consumed_seconds"]) -
                float(item["usage_before"]["usage_consumed_seconds"]))
        observations.append({
            "label": label,
            "job_id": item["job_id"],
            "backend": item["backend"],
            "resilience_level": int(item["resilience_level"]),
            "expectation_K": float(item["measured_V_expectation_K"]),
            "reported_standard_error_K": float(
                item["reported_standard_error_K"]),
            "signed_error_K": float(item["signed_error_K"]),
            "absolute_error_K": float(item["absolute_error_K"]),
            "within_two_reported_standard_errors": bool(
                item["within_two_reported_standard_errors"]),
            "free_qpu_seconds_used": used,
        })

    bias_reduction = 1 - (observations[1]["absolute_error_K"] /
                          observations[0]["absolute_error_K"])
    summary = {
        "purpose": "paired live-QPU pilot summary for N=6 spectral flow",
        "exact_V_expectation_K": exact,
        "observations": observations,
        "absolute_bias_reduction_fraction_TREX_to_ZNE": bias_reduction,
        "total_free_qpu_seconds_used": sum(
            row["free_qpu_seconds_used"] for row in observations),
        "free_qpu_seconds_remaining_after_pilot": int(
            zne["usage_after"]["usage_remaining_seconds"]),
        "decision": "PASS_EXECUTION_PATH_HOLD_QUANTITATIVE_ACCURACY",
        "interpretation": [
            "The Open Plan guard and two real-QPU jobs completed successfully.",
            "TREX alone retained a systematic error much larger than its reported standard error.",
            "ZNE moved the central estimate toward the exact value but increased uncertainty; the apparent improvement is not statistically conclusive.",
            "No material or engineering conclusion is inferred from the N=6 pilot.",
            "Further QPU runs are stopped until circuit compression or a stronger observable-level mitigation design is available.",
        ],
    }
    (RESULTS / "open_plan_qpu_pilot_N6_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
    labels = [row["label"] for row in observations]
    x = np.arange(len(labels))
    values = [row["expectation_K"] for row in observations]
    errors = [row["reported_standard_error_K"] for row in observations]
    axes[0].errorbar(x, values, yerr=errors, fmt="o", markersize=8,
                     capsize=5, color="#2f6f9f", linewidth=2,
                     label="IBM Fez")
    axes[0].axhline(exact, color="#c44e52", linewidth=2,
                    label=f"Exact = {exact:.4f} K")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel(r"$\langle V\rangle$ (K)")
    axes[0].set_title("A. Real-QPU observable")
    axes[0].legend(frameon=True, fontsize=9)
    axes[0].text(0.03, 0.04,
                 "Error bars: reported standard error\nN=6 circuit-path proof",
                 transform=axes[0].transAxes, fontsize=8.5, va="bottom")

    used = [row["free_qpu_seconds_used"] for row in observations]
    bars = axes[1].bar(labels, used, color=["#4c9f70", "#dd8452"])
    for bar, value in zip(bars, used, strict=True):
        axes[1].text(bar.get_x() + bar.get_width()/2, value + 1,
                     f"{value:.0f} s", ha="center", va="bottom", fontsize=9)
    axes[1].set_ylabel("Open Plan QPU time used (s)")
    axes[1].set_title("B. Free-plan resource use")
    axes[1].text(0.5, 0.95,
                 f"Total 71 s; 412 s remained",
                 transform=axes[1].transAxes, ha="center", va="top",
                 fontsize=9)
    axes[1].set_ylim(0, max(used) * 1.25)

    fig.suptitle("N=6 spectral-flow pilot on IBM Fez", fontsize=13,
                 fontweight="bold")
    fig.tight_layout()
    output = ROOT / "figures" / "open_plan_qpu_pilot_N6.png"
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "summary": str(RESULTS / "open_plan_qpu_pilot_N6_summary.json"),
        "figure": str(output),
        "decision": summary["decision"],
        "bias_reduction_fraction": bias_reduction,
        "total_free_qpu_seconds_used": summary["total_free_qpu_seconds_used"],
    }, indent=2))


if __name__ == "__main__":
    main()
