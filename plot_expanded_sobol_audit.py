"""Plot the expanded 256-point interlayer Sobol audit."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def main() -> None:
    all_rows = pd.read_csv(RESULTS / "interlayer_expanded_sobol.csv")
    retained = all_rows[all_rows["propagation_vector_pass"].astype(bool)].copy()
    engineering = pd.read_csv(RESULTS / "interlayer_expanded_sobol_engineering.csv")
    e90 = engineering[engineering["temperature_K"].round(6) == 0.09]
    central_pass = ((e90["delta_B_star_T"] >= 0.05) &
                    (e90["regret_percent"] >= 5.0)).sum()

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))

    ax = axes[0]
    ax.scatter(all_rows["delta_inplane_rlu"], all_rows["qz_folded_rlu"],
               s=16, color="#B8B8B8", alpha=0.55, label="outside filter")
    points = ax.scatter(retained["delta_inplane_rlu"], retained["qz_folded_rlu"],
                        c=retained["critical_field_shift_vs_M3_T"], s=58,
                        cmap="viridis", edgecolor="black", linewidth=0.45,
                        label="retained")
    ax.axvspan(0.021, 0.033, color="#4C78A8", alpha=0.09)
    ax.axhspan(0.322, 1 / 3 + 0.004, color="#F58518", alpha=0.09)
    ax.errorbar([0.027], [0.328], xerr=[0.002], yerr=[0.002], fmt="o",
                color="black", markersize=4, label="public 65 mK")
    ax.set(xlabel=r"in-plane $\delta_c$ (r.l.u.)",
           ylabel=r"folded $q_z$ (r.l.u.)",
           title="Expanded public-q filter")
    ax.legend(fontsize=7.7, loc="best")
    fig.colorbar(points, ax=ax, label=r"$\Delta H_c$ vs M3 (T)")

    ax = axes[1]
    inside = retained[retained["inside_original_parameter_box"].astype(bool)]
    outside = retained[~retained["inside_original_parameter_box"].astype(bool)]
    for data, marker, label in ((inside, "o", "inside local box"),
                                (outside, "s", "outside local box")):
        sc = ax.scatter(data["j2_microev"], data["jb_microev"],
                        c=data["critical_field_shift_vs_M3_T"], cmap="viridis",
                        vmin=retained["critical_field_shift_vs_M3_T"].min(),
                        vmax=retained["critical_field_shift_vs_M3_T"].max(),
                        s=85, marker=marker, edgecolor="black", linewidth=0.55,
                        label=label)
    for _, row in retained.iterrows():
        ax.annotate(f"{row['critical_field_shift_vs_M3_T']:.3f} T",
                    (row["j2_microev"], row["jb_microev"]),
                    xytext=(4, 5), textcoords="offset points", fontsize=7.5)
    ax.set(xlabel=r"$J_2$ ($\mu$eV)", ylabel=r"$J_b$ ($\mu$eV)",
           title="Retained simultaneous variations")
    ax.legend(fontsize=7.7)

    ax = axes[2]
    scatter = ax.scatter(e90["delta_B_star_T"], e90["regret_percent"],
                         c=e90["critical_field_shift_vs_M3_T"], cmap="viridis",
                         s=48, alpha=0.8, edgecolor="black", linewidth=0.3)
    ax.axvline(0.05, color="#E45756", linestyle="--", label="0.05 T gate")
    ax.axhline(5.0, color="#E45756", linestyle=":", label="5% gate")
    ax.set(xlabel=r"90 mK $\Delta B^*$ (T)", ylabel="regret (%)",
           title=f"Central-profile gate: {central_pass}/{len(e90)} pass")
    ax.legend(fontsize=7.7, loc="best")
    fig.colorbar(scatter, ax=ax, label=r"$\Delta H_c$ vs M3 (T)")

    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        kwargs = {"dpi": 240} if suffix == "png" else {}
        fig.savefig(FIGURES / f"interlayer_expanded_sobol_audit.{suffix}",
                    bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()
