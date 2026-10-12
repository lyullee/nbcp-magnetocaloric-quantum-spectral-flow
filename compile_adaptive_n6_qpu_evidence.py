"""Compile paired live-QPU evidence, tables, and the manuscript figure."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def main():
    local = load("adaptive_compressed_stateprep_N6.json")
    trex = load("adaptive_n6_qpu_complete_TWIRLED_TREX.json")
    zne = load("adaptive_n6_qpu_complete_TWIRLED_TREX_ZNE.json")
    if trex["backend"] != zne["backend"]:
        raise RuntimeError("paired jobs used different backends")
    if trex["condition"] != zne["condition"]:
        raise RuntimeError("paired jobs used different physical conditions")
    if len(trex["observations"]) != 9 or len(zne["observations"]) != 9:
        raise RuntimeError("expected nine observations per job")

    rows = []
    for a, b in zip(trex["observations"], zne["observations"], strict=True):
        if (a["lambda"], a["level"]) != (b["lambda"], b["level"]):
            raise RuntimeError("observation ordering mismatch")
        rows.append({
            "lambda": a["lambda"], "level": a["level"],
            "exact_V_K": a["exact_V_expectation_K"],
            "prepared_V_K": a["prepared_V_expectation_K"],
            "twirled_trex_V_K": a["hardware_V_expectation_K"],
            "twirled_trex_standard_error_K": a["reported_standard_error_K"],
            "twirled_trex_zne_V_K": b["hardware_V_expectation_K"],
            "twirled_trex_zne_standard_error_K": b["reported_standard_error_K"],
        })
    csv_path = RESULTS / "adaptive_n6_qpu_observations.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    exact = trex["thermodynamics"]["exact_full_delta_S_kB_per_spin"]
    prepared = trex["thermodynamics"]["prepared_low_rank_delta_S_kB_per_spin"]
    trex_s = trex["thermodynamics"]["hardware_low_rank_delta_S_kB_per_spin"]
    zne_s = zne["thermodynamics"]["hardware_low_rank_delta_S_kB_per_spin"]
    trex_u = trex["thermodynamics"]["reported_uncertainty_propagation"]
    zne_u = zne["thermodynamics"]["reported_uncertainty_propagation"]
    trex_bias = abs(trex_s - exact)
    zne_bias = abs(zne_s - exact)
    summary = {
        "purpose": "paired complete N=6 live-QPU spectral-flow benchmark",
        "backend": trex["backend"],
        "jobs": {"twirled_TREX": trex["job_id"],
                 "twirled_TREX_ZNE": zne["job_id"]},
        "qubits": int(local["condition"]["N"]),
        "circuits_per_job": len(rows),
        "maximum_two_qubit_gate_count": max(
            row["two_qubit_gate_count"] for row in zne["circuits"]),
        "target_precision_K": trex["target_precision_K"],
        "local_circuit": {
            "minimum_fidelity": local["minimum_fidelity"],
            "maximum_absolute_V_error_K": local["maximum_absolute_V_error_K"],
            "delta_S_error_kB_per_spin": local["delta_S_error_kB_per_spin"],
            "passes_predeclared_entropy_target": local["passes_entropy_target"],
        },
        "entropy_kB_per_spin": {
            "exact_full": exact,
            "prepared_low_rank": prepared,
            "twirled_TREX": trex_s,
            "twirled_TREX_95_interval": trex_u["95_percent_interval_kB_per_spin"],
            "twirled_TREX_ZNE": zne_s,
            "twirled_TREX_ZNE_95_interval": zne_u["95_percent_interval_kB_per_spin"],
        },
        "accuracy": {
            "twirled_TREX_absolute_bias": trex_bias,
            "twirled_TREX_ZNE_absolute_bias": zne_bias,
            "ZNE_bias_reduction_fraction": (trex_bias - zne_bias) / trex_bias,
            "ZNE_to_TREX_reported_sd_ratio": (
                zne_u["standard_deviation_kB_per_spin"] /
                trex_u["standard_deviation_kB_per_spin"]),
            "both_signs_match_exact": (
                trex["thermodynamics"]["hardware_sign_matches_exact"] and
                zne["thermodynamics"]["hardware_sign_matches_exact"]),
            "both_intervals_exclude_zero": (
                trex_u["does_not_cover_zero"] and
                zne_u["does_not_cover_zero"]),
            "quantitative_accuracy_pass": (
                trex_u["covers_exact_full"] or zne_u["covers_exact_full"]),
        },
        "open_plan_qpu_seconds": {
            "twirled_TREX": (
                trex["usage_after"]["usage_consumed_seconds"] -
                trex["usage_before"]["usage_consumed_seconds"]),
            "twirled_TREX_ZNE": (
                zne["usage_after"]["usage_consumed_seconds"] -
                zne["usage_before"]["usage_consumed_seconds"]),
        },
        "decision": "PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE",
        "scope": "N=6 open-patch hardware method benchmark only; N=18 and bulk engineering conclusions remain classical/public-data results",
    }
    (RESULTS / "adaptive_n6_qpu_complete_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    exact_v = np.asarray([row["exact_V_K"] for row in rows])
    trex_v = np.asarray([row["twirled_trex_V_K"] for row in rows])
    trex_e = np.asarray([row["twirled_trex_standard_error_K"] for row in rows])
    zne_v = np.asarray([row["twirled_trex_zne_V_K"] for row in rows])
    zne_e = np.asarray([row["twirled_trex_zne_standard_error_K"] for row in rows])
    labels = [f"{row['lambda']:.2f}/L{row['level']}" for row in rows]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10,
                         "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5),
                             gridspec_kw={"width_ratios": [1.05, 1.45, 1.05]})
    ax = axes[0]
    limit = max(exact_v.max(), (zne_v + 2 * zne_e).max()) * 1.08
    ax.plot([0, limit], [0, limit], color="0.4", linestyle="--", linewidth=1,
            label="identity")
    ax.errorbar(exact_v, trex_v, yerr=2 * trex_e, fmt="o", markersize=4,
                capsize=2, label="twirled TREX")
    ax.errorbar(exact_v, zne_v, yerr=2 * zne_e, fmt="s", markersize=4,
                capsize=2, label="twirled TREX + ZNE")
    ax.set(xlabel="Exact $\\langle V\\rangle$ (K)",
           ylabel="IBM Kingston $\\langle V\\rangle$ (K)",
           xlim=(0, limit), ylim=(0, limit), title="(a) Nine QPU observables")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    x = np.arange(len(rows))
    ax.plot(x, exact_v, color="black", marker="o", markersize=3,
            linewidth=1, label="exact")
    ax.errorbar(x - 0.08, trex_v, yerr=2 * trex_e, fmt="o", markersize=3,
                capsize=1.5, label="twirled TREX")
    ax.errorbar(x + 0.08, zne_v, yerr=2 * zne_e, fmt="s", markersize=3,
                capsize=1.5, label="twirled TREX + ZNE")
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.set_ylabel("$\\langle V\\rangle$ (K)")
    ax.set_title("(b) Spectral-flow node / level")
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper center")

    ax = axes[2]
    names = ["exact", "prepared", "TREX", "TREX\n+ ZNE"]
    values = np.asarray([exact, prepared, trex_s, zne_s])
    lower_t, upper_t = trex_u["95_percent_interval_kB_per_spin"]
    lower_z, upper_z = zne_u["95_percent_interval_kB_per_spin"]
    yerr = np.zeros((2, 4))
    yerr[:, 2] = [trex_s - lower_t, upper_t - trex_s]
    yerr[:, 3] = [zne_s - lower_z, upper_z - zne_s]
    colors = ["0.25", "0.55", "#377eb8", "#e41a1c"]
    ax.bar(np.arange(4), values, color=colors, width=0.72)
    ax.errorbar(np.arange(4), values, yerr=yerr, fmt="none", color="black",
                capsize=3, linewidth=1)
    ax.axhline(0, color="0.2", linewidth=0.8)
    ax.set_xticks(np.arange(4), names)
    ax.set_ylabel("$\\Delta S/N$ ($k_B$ spin$^{-1}$)")
    ax.set_title("(c) Reconstructed entropy shift")
    ax.text(0.03, 0.97, "sign: PASS\nmagnitude: HOLD", transform=ax.transAxes,
            va="top", ha="left", fontsize=8)

    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"adaptive_n6_qpu_complete.{suffix}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "summary": str(RESULTS / "adaptive_n6_qpu_complete_summary.json"),
        "observations": str(csv_path),
        "figure": str(FIGURES / "adaptive_n6_qpu_complete.png"),
        "decision": summary["decision"],
        "ZNE_bias_reduction_fraction": summary["accuracy"]["ZNE_bias_reduction_fraction"],
    }, indent=2))


if __name__ == "__main__":
    main()
