"""Compile the first local quantum-algorithm screening figure and summary."""
from __future__ import annotations

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
    mixing = load("qmetts_transition_audit_N9.json")
    base = load("qite_tangent_screen_N9.json")
    adaptive = load("qite_adaptive_weight3_N9.json")
    shots = load("qpu_shot_budget_N9.json")
    n6 = load("qite_qmetts_xbasis_N6.json")

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.3))
    colors = {0.09: "#4c78a8", 0.10: "#f58518", 0.12: "#54a24b"}
    for temperature in sorted(colors):
        rows = [row for row in mixing["rows"] if row["temperature_K"] == temperature]
        axes[0].plot([row["field_T"] for row in rows],
                     [row["transition_spectral_gap"] for row in rows],
                     marker="o", color=colors[temperature],
                     label=f"{temperature*1000:.0f} mK")
    axes[0].set(xlabel="Field (T)", ylabel="X-basis transition spectral gap",
                title="(a) N=9 QMETTS mixing", ylim=(0, 0.65))
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    markers = {0.5: "o", 2.0: "s"}
    for tau in markers:
        base_row = next(row for row in base["rows"]
                        if row["correlation_time_inverse_K"] == tau and
                        row["maximum_pauli_weight"] == 2)
        rows = [base_row] + [row for row in adaptive["rows"]
                             if row["correlation_time_inverse_K"] == tau]
        axes[1].plot([row.get("compiled_all_to_all_CX", 0) for row in rows],
                     [row["relative_tangent_residual"] for row in rows],
                     marker=markers[tau], label=rf"$\tau={tau:g}$ K$^{{-1}}$")
    axes[1].axhline(0.05, color="black", lw=1, ls="--", label="residual target")
    axes[1].set(xlabel="CX per QITE step (all-to-all)",
                ylabel="Relative tangent residual", title="(b) Accuracy–circuit tradeoff",
                yscale="log")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25, which="both")

    fields = [row["field_T"] for row in shots["rows"]]
    general = [row["budgets"]["general"]["optimally_allocated_total_shots"]
               for row in shots["rows"]]
    axes[2].bar(np.arange(len(fields)), general, color="#e45756", label=r"$\sigma_S=10^{-3}$")
    stress_row = next(row for row in shots["rows"] if abs(row["field_T"] - 1.45) < 1e-12)
    stress = stress_row["budgets"]["stress_sign_change"]["optimally_allocated_total_shots"]
    axes[2].bar(1, stress, color="#b279a2", width=0.5, label=r"1.45 T, $\sigma_S=2\times10^{-4}$")
    axes[2].set_xticks(np.arange(len(fields)), [f"{field:.2f}" for field in fields])
    axes[2].set(xlabel="Field (T)", ylabel="Optimistic total shots",
                title="(c) Entropy-correction measurement cost", yscale="log")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].grid(alpha=0.25, axis="y", which="both")

    fig.suptitle("Local quantum screening: viable sampling and the resource target for direct entropy estimators", y=1.03)
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    figure_path = FIGURES / "quantum_algorithm_screening.png"
    fig.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "N6_end_to_end_ideal_QITE_QMETTS": {
            "bias_K_per_spin": n6["end_to_end_bias_K_per_spin"],
            "minimum_state_fidelity": n6["statewise"]["minimum_fidelity"],
            "transition_spectral_gap": n6["transition_spectral_summary"]["spectral_gap"],
            "decision": "PASS_ideal_observable_path_only",
        },
        "N9_exact_QMETTS_sampling": {
            "spectral_gap_range": [min(row["transition_spectral_gap"] for row in mixing["rows"]),
                                   max(row["transition_spectral_gap"] for row in mixing["rows"])],
            "independent_sample_range_for_1e-4_K_per_spin": [
                min(row["independent_samples_for_1e-4_K_per_spin_SE"] for row in mixing["rows"]),
                max(row["independent_samples_for_1e-4_K_per_spin_SE"] for row in mixing["rows"])],
            "decision": "PASS_sampling_only",
        },
        "N9_QITE": {
            "weight2_correlated_residual_range": [
                min(row["relative_tangent_residual"] for row in base["rows"]
                    if row["maximum_pauli_weight"] == 2 and row["correlation_time_inverse_K"] > 0),
                max(row["relative_tangent_residual"] for row in base["rows"]
                    if row["maximum_pauli_weight"] == 2 and row["correlation_time_inverse_K"] > 0)],
            "adaptive_best_residual_range": [
                min(row["relative_tangent_residual"] for row in adaptive["rows"]),
                max(row["relative_tangent_residual"] for row in adaptive["rows"]
                    if row["selected_weight3_count"] == 256)],
            "decision": "STOP_current_QPU_depth",
        },
        "finite_difference_entropy_shots": {
            "general_range": [min(general), max(general)],
            "stress_1p45T": stress,
            "decision": "STOP_Open_Plan",
        },
        "figure": str(figure_path),
    }
    summary_path = RESULTS / "quantum_algorithm_screening_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
