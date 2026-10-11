"""M0 near-saturation thermodynamics on three-sublattice periodic clusters."""
import argparse
import json
import time
from math import comb
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from magnon_sector import (
    FixedMagnonXXZOperator,
    deflated_sector_rules,
    deflated_slq_thermodynamics,
    periodic_graph_edges,
)
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import DELTA_XY_OVER_ZZ, GZ, JZZ_MEV, MEV_TO_K

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CELLS = {
    18: ((4, 1), (-2, 4), 4),
    21: ((-5, 1), (-1, -4), 5),
    24: ((-5, 1), (1, -5), 6),
    27: ((-6, 3), (-3, -3), 6),
}


def main(n: int):
    if n not in CELLS:
        raise ValueError(f"supported N: {sorted(CELLS)}")
    a, b, cutoff = CELLS[n]
    field = 1.50
    temperatures = [0.09, 0.10, 0.12]
    low_states = 8
    random_vectors = 2 if n >= 27 else 4
    krylov_steps = 50 if n >= 27 else 60
    graph = lattice(a, b)
    edges = periodic_graph_edges(graph)
    jzz = JZZ_MEV * MEV_TO_K
    jxy = jzz * DELTA_XY_OVER_ZZ
    h = GZ * MU_B_OVER_KB * field
    rng = np.random.default_rng(20261011 + n)
    sectors = []
    metadata = []
    started = time.perf_counter()
    for q in range(cutoff + 1):
        op = FixedMagnonXXZOperator(n, q, edges, jxy, jzz, h)
        rule = deflated_sector_rules(op, low_states, random_vectors, krylov_steps, rng)
        sectors.append(rule)
        metadata.append({
            "q_down": q,
            "dimension": op.shape[0],
            "ground_energy_K": float(rule.low_energies.min()),
            "elapsed_seconds_cumulative": time.perf_counter() - started,
        })
        print(f"N={n} q={q} dim={op.shape[0]} E0={rule.low_energies.min():.9f} elapsed={time.perf_counter()-started:.1f}s", flush=True)

    q_next = cutoff + 1
    op_next = FixedMagnonXXZOperator(n, q_next, edges, jxy, jzz, h)
    next_ground = float(eigsh(op_next, k=1, which="SA", tol=1e-9, return_eigenvectors=False)[0])
    rows = deflated_slq_thermodynamics(sectors, temperatures, n)
    ground = min(item["ground_energy_K"] for item in metadata)
    for row in rows:
        row["next_sector_partition_upper_log"] = (
            np.log(comb(n, q_next)) - (next_ground - ground) / row["temperature_K"]
        )
    payload = {
        "method": "fixed-magnon deflated stochastic Lanczos quadrature",
        "condition": {"N": n, "translations": [a, b], "field_T": field, "model": "M0 2026 parameters"},
        "settings": {"max_q_down": cutoff, "exact_low_states_per_sector": low_states, "random_vectors_per_sector": random_vectors, "krylov_steps": krylov_steps, "seed": 20261011 + n},
        "sectors": metadata,
        "first_omitted_sector": {"q_down": q_next, "dimension": op_next.shape[0], "ground_energy_K": next_ground},
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    out = RESULTS / f"magnon_m0_N{n}_B1.500T.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": out.name, "rows": rows, "first_omitted_sector": payload["first_omitted_sector"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    main(parser.parse_args().n)
