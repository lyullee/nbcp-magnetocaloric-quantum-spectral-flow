"""Check whether the weak-interaction entropy correction converges faster than S itself.

The absolute finite-cluster entropy has large shell effects near the field-driven
transition.  This script tests the controlled alternative

    S_bulk(M3) ~= S_bulk(M0) + lim_N [S_N(M3) - S_N(M0)].

No extrapolation produced here is labelled a confidence interval.  With only
N=12, 15, and 18 it is a diagnostic and an uncertainty bracket.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"


SOURCES = [
    (9, "N9", OUT / "small_size_sequence_raw", "N9_{model}_B{field:.3f}T.json"),
    (12, "N12", OUT / "small_size_sequence_raw", "N12_{model}_B{field:.3f}T.json"),
    (15, "N15", OUT / "n15_confirm_raw", "{model}_B{field:.3f}T.json"),
    (18, "N18-A", OUT / "n18_pilot_raw", "{model}_B{field:.3f}T.json"),
    (18, "N18-B", OUT / "n18_compact_shape_raw", "{model}_B{field:.3f}T.json"),
]

FIELDS = (1.10, 1.15, 1.20)
TEMPERATURES = (0.09, 0.10, 0.12)


def read_entropy(path: Path, temperature: float) -> float:
    d = json.loads(path.read_text(encoding="utf-8"))
    row = min(d["rows"], key=lambda r: abs(r["temperature_K"] - temperature))
    if abs(row["temperature_K"] - temperature) > 1e-9:
        raise ValueError(f"temperature {temperature} missing in {path}")
    return float(row["entropy_kB_per_spin"])


def collect() -> list[dict]:
    rows = []
    for n, cluster, directory, pattern in SOURCES:
        for b in FIELDS:
            paths = {
                m: directory / pattern.format(model=m, field=b) for m in ("M0", "M3")
            }
            if not all(p.exists() for p in paths.values()):
                continue
            for t in TEMPERATURES:
                s0 = read_entropy(paths["M0"], t)
                s3 = read_entropy(paths["M3"], t)
                rows.append(
                    {
                        "N": n,
                        "cluster": cluster,
                        "temperature_K": t,
                        "field_T": b,
                        "S_M0_kB_per_spin": s0,
                        "S_M3_kB_per_spin": s3,
                        "deltaS_kB_per_spin": s3 - s0,
                    }
                )
    return rows


def diagnostics(rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for r in rows:
        if r["N"] >= 12:
            grouped[(r["temperature_K"], r["field_T"])].append(r)

    result = []
    for (t, b), rs in sorted(grouped.items()):
        # Average shapes before fitting so N=18 is not double weighted.
        by_n = defaultdict(list)
        for r in rs:
            by_n[r["N"]].append(r["deltaS_kB_per_spin"])
        ns = np.array(sorted(by_n), dtype=float)
        yn = np.array([np.mean(by_n[int(n)]) for n in ns])
        intercept = slope = None
        if len(ns) >= 3:
            slope, intercept = np.polyfit(1.0 / ns, yn, 1)

        largest = [r["deltaS_kB_per_spin"] for r in rs if r["N"] >= 15]
        vals = [r["deltaS_kB_per_spin"] for r in rs]
        sign_consistent = all(v > 0 for v in vals) or all(v < 0 for v in vals)
        result.append(
            {
                "temperature_K": t,
                "field_T": b,
                "clusters": {r["cluster"]: r["deltaS_kB_per_spin"] for r in rs},
                "sign_consistent_N12_to_N18": sign_consistent,
                "all_cluster_range": [min(vals), max(vals)],
                "largest_cluster_range_N15_N18": [min(largest), max(largest)],
                "largest_cluster_midpoint": float((min(largest) + max(largest)) / 2),
                "linear_in_1_over_N_intercept": None if intercept is None else float(intercept),
                "linear_fit_is_diagnostic_not_CI": True,
            }
        )
    return result


def write_outputs(rows: list[dict], diag: list[dict]) -> None:
    csv_path = OUT / "deltaS_correction_size_sequence.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    payload = {
        "quantity": "deltaS = S(M3)-S(M0), in kB per spin",
        "purpose": "test finite-size cancellation before combining with a bulk M0 baseline",
        "interpretation_rule": {
            "usable": "sign stable and N15/N18 spread small relative to the correction",
            "not_a_confidence_interval": True,
            "quantum_stage_authorized": False,
        },
        "rows_csv": csv_path.name,
        "diagnostics": diag,
    }
    (OUT / "deltaS_correction_diagnostic.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.8), sharey=True)
    colors = {"N9": "#bdbdbd", "N12": "#4c78a8", "N15": "#f58518", "N18-A": "#54a24b", "N18-B": "#e45756"}
    markers = {"N9": "o", "N12": "s", "N15": "^", "N18-A": "D", "N18-B": "P"}
    for ax, t in zip(axes, TEMPERATURES):
        for cluster in colors:
            q = sorted(
                [r for r in rows if r["cluster"] == cluster and r["temperature_K"] == t],
                key=lambda r: r["field_T"],
            )
            if not q:
                continue
            ax.plot(
                [r["field_T"] for r in q],
                [r["deltaS_kB_per_spin"] for r in q],
                marker=markers[cluster], color=colors[cluster], lw=1.5, ms=5, label=cluster,
            )
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"T = {t*1000:.0f} mK")
        ax.set_xlabel("Field B (T)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$\Delta S_N=S_N(M3)-S_N(M0)$ ($k_B$/spin)")
    axes[-1].legend(frameon=False, fontsize=8, ncol=1)
    fig.suptitle("Finite-cluster cancellation test for the realistic-interaction correction")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "deltaS_correction_size_sequence.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    data = collect()
    d = diagnostics(data)
    write_outputs(data, d)
    for row in d:
        if row["temperature_K"] == 0.09:
            lo, hi = row["largest_cluster_range_N15_N18"]
            print(
                f"90 mK B={row['field_T']:.2f} T: "
                f"sign={row['sign_consistent_N12_to_N18']}, "
                f"N15/N18=[{lo:+.6f}, {hi:+.6f}], "
                f"1/N intercept={row['linear_in_1_over_N_intercept']:+.6f}"
            )
