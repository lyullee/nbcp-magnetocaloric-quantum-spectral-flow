"""Compile the second-stage quantum-algorithm screening figure."""
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
    direct = load("qpu_shot_budget_N9.json")
    covariance = load("covariance_entropy_budget_N9.json")
    low = load("low_energy_quantum_route_N9.json")
    flow9 = load("spectral_flow_entropy_N9.json")
    flow12 = load("spectral_flow_entropy_N12.json")
    complexity = load("low_energy_state_complexity_N9_N12.json")

    d = next(row for row in direct["rows"] if row["field_T"] == 1.45)
    c = next(row for row in covariance["rows"] if row["field_T"] == 1.45)
    l = next(row for row in low["delta_S_comparisons"]
             if row["field_T"] == 1.45 and row["temperature_K"] == 0.10)
    f9 = next(row for row in flow9["rows"]
              if row["field_T"] == 1.45 and row["temperature_K"] == 0.10
              and row["quadrature_order"] == 3)
    f12 = next(row for row in flow12["rows"]
               if row["field_T"] == 1.45 and row["temperature_K"] == 0.10)
    labels = ["Finite\ndifference", "Ideal\ncovariance", "Low spectrum\nfull H",
              "Spectral flow\nN=9", "Spectral flow\nN=12"]
    shots = [d["budgets"]["general"]["optimally_allocated_total_shots"],
             c["budgets"]["general"]["ideal_mean_dependent_covariance_observable_shots"],
             l["optimistic_exact_eigenstate_energy_measurement"]["optimal_total_shots"],
             f9["optimistic_optimal_V_measurement_shots"],
             f12["optimistic_optimal_V_measurement_shots"]]

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.3), constrained_layout=True)
    ax = axes[0, 0]
    colors = ["#9e9e9e", "#e45756", "#f2cf5b", "#4c78a8", "#54a24b"]
    bars = ax.bar(np.arange(len(shots)), shots, color=colors)
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(shots)), labels)
    ax.set_ylabel("Optimistic shots")
    ax.set_title(r"(a) B=1.45 T, T=100 mK, $\sigma_S=10^{-3}$")
    for bar, value in zip(bars, shots, strict=True):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.18,
                f"{value:.2g}", ha="center", va="bottom", fontsize=8)
    ax.grid(axis="y", which="both", alpha=0.25)

    ax = axes[0, 1]
    temperatures = np.array([0.09, 0.10, 0.12])
    for payload, label, color, marker in ((flow9, "N=9", "#4c78a8", "o"),
                                           (flow12, "N=12", "#54a24b", "s")):
        mins, maxs, centers = [], [], []
        for temperature in temperatures:
            data = [row["low_energy_rank"] for row in payload["rows"]
                    if row["temperature_K"] == float(temperature)
                    and row.get("quadrature_order", 3) == 3]
            mins.append(min(data)); maxs.append(max(data)); centers.append(np.mean(data))
        ax.plot(temperatures * 1000, centers, color=color, marker=marker,
                label=label)
        ax.fill_between(temperatures * 1000, mins, maxs, color=color, alpha=0.18)
    ax.set_xlabel("Temperature (mK)")
    ax.set_ylabel("Selected low-energy levels K")
    ax.set_title("(b) Rank required by the entropy-error rule")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    for cluster, color, marker in zip(complexity["clusters"],
                                      ("#4c78a8", "#54a24b"), ("o", "s"),
                                      strict=True):
        nodes = [row["lambda"] for row in cluster["nodes"]]
        chi = [row["worst_state_summary"][
            "maximum_MPS_bond_rank_for_99_9pct_per_cut"]
               for row in cluster["nodes"]]
        ax.plot(nodes, chi, marker=marker, color=color, label=f"N={cluster['N']}")
    ax.set_xlabel(r"Interaction path $\lambda$")
    ax.set_ylabel("Worst 99.9% per-cut MPS bond rank")
    ax.set_title("(c) Low-energy state compressibility at 1.45 T")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    temperature_colors = {0.09: "#4c78a8", 0.10: "#f58518", 0.12: "#54a24b"}
    for payload, marker, size_label in ((flow9, "o", "N=9"),
                                         (flow12, "^", "N=12")):
        for row in payload["rows"]:
            if row["temperature_K"] not in temperature_colors:
                continue
            if row.get("quadrature_order", 3) != 3:
                continue
            error_key = ("combined_low_rank_and_quadrature_error_kB_per_spin"
                         if "combined_low_rank_and_quadrature_error_kB_per_spin" in row
                         else "low_rank_plus_quadrature_error_vs_reference_kB_per_spin")
            ratio = abs(row[error_key]) / row["target_delta_S_error_kB_per_spin"]
            ax.scatter(row["field_T"], ratio, marker=marker, s=58,
                       color=temperature_colors[row["temperature_K"]],
                       edgecolor="white", linewidth=0.5)
    ax.axhline(1, color="black", linestyle="--", linewidth=1)
    ax.set_yscale("log")
    ax.set_xlabel("Field (T)")
    ax.set_ylabel("Absolute error / preregistered tolerance")
    ax.set_title("(d) All N=9/N=12 spectral-flow points pass")
    ax.grid(alpha=0.25, which="both")
    from matplotlib.lines import Line2D
    legend = [Line2D([0], [0], marker="o", color="none", markerfacecolor="gray",
                     label="N=9", markersize=7),
              Line2D([0], [0], marker="^", color="none", markerfacecolor="gray",
                     label="N=12", markersize=7)]
    legend.extend(Line2D([0], [0], marker="o", color="none",
                         markerfacecolor=color, label=f"{int(t*1000)} mK",
                         markersize=7)
                  for t, color in temperature_colors.items())
    ax.legend(handles=legend, frameon=False, ncol=2, fontsize=8)

    FIGURES.mkdir(exist_ok=True)
    output = FIGURES / "quantum_algorithm_second_screen.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)

    summary = {
        "purpose": "second-stage quantum algorithm screening summary",
        "comparison_condition": {"field_T": 1.45, "temperature_K": 0.10,
                                 "entropy_error_kB_per_spin": 0.001},
        "shot_comparison": dict(zip(labels, shots, strict=True)),
        "selected_route": "low-energy spectral flow of V=M3-M0",
        "N9_general_shot_range": [min(row["optimistic_optimal_V_measurement_shots"]
                                      for row in flow9["rows"]
                                      if row["quadrature_order"] == 3 and
                                      row["target_delta_S_error_kB_per_spin"] == 0.001),
                                  max(row["optimistic_optimal_V_measurement_shots"]
                                      for row in flow9["rows"]
                                      if row["quadrature_order"] == 3 and
                                      row["target_delta_S_error_kB_per_spin"] == 0.001)],
        "N12_general_shot_range": [min(row["optimistic_optimal_V_measurement_shots"]
                                       for row in flow12["rows"]
                                       if row["target_delta_S_error_kB_per_spin"] == 0.001),
                                   max(row["optimistic_optimal_V_measurement_shots"]
                                       for row in flow12["rows"]
                                       if row["target_delta_S_error_kB_per_spin"] == 0.001)],
        "figure": str(output),
    }
    (RESULTS / "quantum_algorithm_second_screen_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
