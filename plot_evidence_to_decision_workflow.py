"""Create the evidence-to-decision workflow used in the manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"

COLORS = {
    "evidence": "#DCEAF7",
    "physics": "#E6E0F2",
    "decision": "#DDEEDB",
    "quantum": "#FCE5CD",
    "scope": "#F2F2F2",
    "edge": "#303030",
    "gate": "#B33A3A",
}


def box(ax, x, y, width, height, facecolor, title, lines, title_color="#111111"):
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.15, edgecolor=COLORS["edge"], facecolor=facecolor,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025 * width, y + height - 0.09 * height, title,
            transform=ax.transAxes, ha="left", va="top", fontsize=12.2,
            fontweight="bold", color=title_color)
    ax.text(x + 0.035 * width, y + height - 0.24 * height, "\n".join(lines),
            transform=ax.transAxes, ha="left", va="top", fontsize=9.15,
            linespacing=1.55, color="#222222")


def arrow(ax, x1, y1, x2, y2, color="#555555", style="-|>"):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=15,
        linewidth=1.6, color=color, transform=ax.transAxes,
        connectionstyle="arc3,rad=0.0",
    ))


def state_box(ax, x, y, width, height, title, line):
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.05, edgecolor=COLORS["edge"], facecolor="#FFFFFF",
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025 * width, y + height - 0.025, title,
            transform=ax.transAxes, ha="left", va="top", fontsize=11.4,
            fontweight="bold", color="#111111")
    ax.text(x + 0.025 * width, y + 0.030, line,
            transform=ax.transAxes, ha="left", va="bottom", fontsize=8.8,
            color="#222222")


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(14.2, 7.5))
    ax.set_axis_off()

    ax.text(0.5, 0.955, "Evidence-to-decision workflow and validation boundaries",
            transform=ax.transAxes, ha="center", va="center", fontsize=16,
            fontweight="bold")
    ax.text(0.5, 0.915,
            "Classical evidence establishes the engineering sensitivity before the quantum stage is attempted",
            transform=ax.transAxes, ha="center", va="center", fontsize=10.6,
            color="#444444")

    width, height, y = 0.215, 0.49, 0.34
    xs = (0.025, 0.275, 0.525, 0.775)
    box(ax, xs[0], y, width, height, COLORS["evidence"], "1  Public evidence", [
        "• entropy at 0, 1.3, 1.5, 1.66 T",
        "• 90–300 mK thermodynamics",
        "• near-saturation propagation vector",
        "• same-volume hold-time traces",
        "• measured transition: 1.65–1.70 T",
    ])
    box(ax, xs[1], y, width, height, COLORS["physics"], "2  Physics screen", [
        "• M0–M4 Hamiltonian hierarchy",
        "• bulk polarized-state instability",
        "• long-range Ewald dipolar sum",
        "• initial exact screen: 631 → 16",
        "• adaptive exact audit: +192 models",
        "• total: 823 → 69 retained models",
        "• extrema verified with five starts",
    ])
    box(ax, xs[2], y, width, height, COLORS["decision"], "3  Refrigerator decision", [
        "• predeclared 0.05 T and 5% gates",
        "• central gates: 224/276 pass",
        "• all-shape gates: 17,892/22,356",
        "• 90 mK minimax field: 1.405 T",
        "• conditional inventory: 4.03 g",
        "• missing-data target: ±2 mT",
        "• 26-point near-saturation sweep",
    ])
    box(ax, xs[3], y, width, height, COLORS["quantum"], "4  Quantum escalation", [
        "• target only the weak-term correction",
        "• measure low-energy derivatives",
        "• ideal N=9 and N=12 pass",
        "• N=6 IBM execution path works",
        "• hardware bias blocks feedback",
        "• no quantum advantage claimed",
    ])

    for first, second in zip(xs[:-1], xs[1:]):
        arrow(ax, first + width + 0.004, y + 0.56 * height,
              second - 0.004, y + 0.56 * height)

    # The red gate shows that the quantum stage is conditional on a material
    # decision effect, rather than being used to select the problem after the fact.
    arrow(ax, xs[2] + width + 0.004, y + 0.56 * height,
          xs[3] - 0.004, y + 0.56 * height, color=COLORS["gate"])

    # Explicit thermodynamic states prevent the magnetic-only inventory from
    # being mistaken for a full refrigerator-cycle prediction.
    ax.text(0.025, 0.285, "Conditional cold-hold state definition",
            transform=ax.transAxes, ha="left", va="center", fontsize=11.2,
            fontweight="bold")
    state_y, state_w, state_h = 0.115, 0.22, 0.115
    state_x = (0.025, 0.315, 0.605)
    state_titles = (
        "A  Magnetized reference",
        "B  Adiabatic endpoint",
        "C  Isothermal hold limit",
    )
    state_lines = (
        "2 K, 4 T; entropy Sₐ",
        "target T and B; entropy remains Sₐ",
        "same T and B; entropy rises to S(T,B)",
    )
    for sx, title, line in zip(state_x, state_titles, state_lines):
        state_box(ax, sx, state_y, state_w, state_h, title, line)
    arrow(ax, state_x[0] + state_w + 0.01, state_y + state_h / 2,
          state_x[1] - 0.01, state_y + state_h / 2)
    arrow(ax, state_x[1] + state_w + 0.01, state_y + state_h / 2,
          state_x[2] - 0.01, state_y + state_h / 2)
    ax.text(0.91, state_y + state_h / 2,
            "Q = T[S(T,B) − Sₐ]",
            transform=ax.transAxes, ha="center", va="center", fontsize=11.2,
            fontweight="bold", color="#2E6B3D")

    ax.add_patch(FancyBboxPatch(
        (0.025, 0.025), 0.95, 0.055,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=0.9, edgecolor="#777777", facecolor=COLORS["scope"],
        transform=ax.transAxes,
    ))
    ax.text(0.5, 0.052,
            "Current boundary: no converged 3D finite-temperature M4 entropy, no complete device parasitic model, and no quantitative QPU feedback",
            transform=ax.transAxes, ha="center", va="center", fontsize=9.3,
            color="#333333")

    fig.savefig(FIGURES / "evidence_to_decision_workflow.png", dpi=240,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "evidence_to_decision_workflow.pdf",
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / "evidence_to_decision_workflow.svg",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
