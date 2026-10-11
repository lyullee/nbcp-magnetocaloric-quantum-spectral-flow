"""Validate a near-saturation magnon cutoff against the existing N=18 FTLM."""
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from magnon_sector import (
    deflated_sector_rules,
    deflated_slq_thermodynamics,
    FixedMagnonXXZOperator,
    periodic_graph_edges,
)
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import DELTA_XY_OVER_ZZ, GZ, JZZ_MEV, MEV_TO_K

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
TEMPERATURES = [0.09, 0.10, 0.12, 0.15, 0.20, 0.30]


def main():
    n, field = 18, 1.50
    max_q, random_vectors, krylov_steps, low_states = 9, 4, 60, 8
    graph = lattice((4, 1), (-2, 4))
    edges = periodic_graph_edges(graph)
    jzz = JZZ_MEV * MEV_TO_K
    jxy = jzz * DELTA_XY_OVER_ZZ
    h = GZ * MU_B_OVER_KB * field
    rng = np.random.default_rng(20261011)
    rules_by_sector = []
    sectors = []
    started = time.perf_counter()
    for q in range(max_q + 1):
        op = FixedMagnonXXZOperator(n, q, edges, jxy, jzz, h)
        sector_rules = deflated_sector_rules(op, low_states, random_vectors, krylov_steps, rng)
        ground = float(sector_rules.low_energies.min())
        rules_by_sector.append(sector_rules)
        sectors.append({
            "q_down": q,
            "dimension": op.shape[0],
            "ground_energy_K": ground,
            "minimum_SLQ_node_K": (min(float(rule.energies.min()) for rule in sector_rules.random_rules)
                                  if sector_rules.random_rules else None),
            "elapsed_seconds_cumulative": time.perf_counter() - started,
        })
        print(f"q={q} dim={op.shape[0]} E0={ground:.9f} elapsed={time.perf_counter()-started:.1f}s", flush=True)

    reference_raw = json.loads((RESULTS / "n18_bulk_window_compact_raw" / "M0_B1.500T.json").read_text())
    reference = {float(row["temperature_K"]): float(row["entropy_kB_per_spin"]) for row in reference_raw["rows"]}
    cutoffs = []
    for cutoff in range(max_q + 1):
        rows = deflated_slq_thermodynamics(rules_by_sector[:cutoff + 1], TEMPERATURES, n)
        for row in rows:
            row["reference_entropy_kB_per_spin"] = reference[row["temperature_K"]]
            row["entropy_error_kB_per_spin"] = row["entropy_kB_per_spin"] - row["reference_entropy_kB_per_spin"]
        cutoffs.append({"max_q_down": cutoff, "rows": rows})

    final_errors = [abs(row["entropy_error_kB_per_spin"]) for row in cutoffs[-1]["rows"]]
    payload = {
        "method": "fixed-magnon-sector stochastic Lanczos quadrature",
        "reference": "existing N=18 parity-space low-energy-corrected FTLM",
        "condition": {"N": n, "field_T": field, "model": "M0 2026 parameters"},
        "settings": {"max_q_down": max_q, "exact_low_states_per_sector": low_states, "random_vectors_per_sector": random_vectors, "krylov_steps": krylov_steps, "seed": 20261011},
        "sectors": sectors,
        "cutoffs": cutoffs,
        "max_abs_entropy_error_at_max_cutoff_kB_per_spin": max(final_errors),
        "elapsed_seconds": time.perf_counter() - started,
        "decision": "PASS" if max(final_errors) <= 0.003 else "FAIL_OR_NEEDS_MORE_SLQ_SAMPLES",
    }
    (RESULTS / "magnon_truncation_n18_validation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for t in TEMPERATURES:
        y = [next(row for row in item["rows"] if row["temperature_K"] == t)["entropy_kB_per_spin"] for item in cutoffs]
        ax.plot(range(max_q + 1), y, "o-", label=f"{int(t*1000)} mK")
        ax.axhline(reference[t], color="k", alpha=.08)
    ax.set(xlabel="Maximum number of down magnons K", ylabel=r"Entropy ($k_B$/spin)", title="N=18 magnon-cutoff validation at 1.50 T")
    ax.grid(alpha=.25); ax.legend(ncol=2, frameon=False)
    fig.tight_layout(); fig.savefig(FIGURES / "magnon_truncation_n18_validation.png", dpi=220); plt.close(fig)
    print(json.dumps({"decision": payload["decision"], "max_error": max(final_errors), "sectors": sectors}, indent=2))


if __name__ == "__main__":
    main()
