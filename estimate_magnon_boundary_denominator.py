"""Estimate second-order cutoff sensitivity using isolated boundary q blocks.

The isolated-block energy is a denominator proxy, not a rigorous lower bound for
the full omitted space, because q=K+2 can itself couple to q=K+4.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ
from truncated_magnon_operator import FixedMagnonCartesianOperator

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main():
    n = 24
    graph = lattice((-5, 1), (1, -5))
    bonds = model_bonds(graph, "M3")
    field_t, temperature = 1.50, 0.09
    field_k = GZ * MU_B_OVER_KB * field_t
    residual_path = RESULTS / "magnon_M3_boundary_residual_N24_B1.500T.json"
    residual_data = json.loads(residual_path.read_text(encoding="utf-8"))
    started = time.perf_counter()
    boundary = {}
    for q in (8, 9):
        op = FixedMagnonCartesianOperator(n, q, bonds, field_k)
        rng = np.random.default_rng(20261011 + q)
        v0 = rng.normal(size=op.shape[0]) + 1j * rng.normal(size=op.shape[0])
        value = float(eigsh(op, k=1, which="SA", tol=1e-8, v0=v0,
                            return_eigenvectors=False)[0].real)
        boundary[q] = {"dimension": op.shape[0], "isolated_block_ground_energy_K": value}
        print(f"q={q} dim={op.shape[0]} E0_block={value:.9f} elapsed={time.perf_counter()-started:.1f}s", flush=True)

    state_rows = []
    for block in residual_data["blocks"]:
        q = int(block["first_omitted_q"])
        e_boundary = boundary[q]["isolated_block_ground_energy_K"]
        for row in block["rows"]:
            gap = e_boundary - row["energy_K"]
            residual = row["first_omitted_block_residual_norm_K"]
            correction = None if gap <= 0 else residual * residual / gap
            state_rows.append({
                "parity": block["parity"],
                "state_index": row["state_index"],
                "energy_K": row["energy_K"],
                "boundary_q": q,
                "isolated_boundary_gap_K": gap,
                "boundary_residual_norm_K": residual,
                "second_order_magnitude_proxy_K": correction,
            })
    energies = np.array([row["energy_K"] for row in state_rows])
    weights = np.exp(-(energies - energies.min()) / temperature)
    weights /= weights.sum()
    corrections = np.array([
        np.nan if row["second_order_magnitude_proxy_K"] is None else row["second_order_magnitude_proxy_K"]
        for row in state_rows
    ])
    summary = {
        "all_isolated_gaps_positive": bool(np.all([row["isolated_boundary_gap_K"] > 0 for row in state_rows])),
        "minimum_isolated_boundary_gap_K": float(min(row["isolated_boundary_gap_K"] for row in state_rows)),
        "maximum_second_order_magnitude_proxy_K": float(np.nanmax(corrections)),
        "boltzmann_weighted_second_order_magnitude_proxy_K": float(np.nansum(weights * corrections)),
        "weighted_proxy_over_temperature": float(np.nansum(weights * corrections) / temperature),
        "weighted_proxy_over_N_temperature": float(np.nansum(weights * corrections) / (n * temperature)),
    }
    payload = {
        "method": "boundary residual squared divided by the isolated first-omitted-q block gap",
        "warning": "Diagnostic proxy only; the isolated boundary-block ground energy is not a rigorous lower bound for the full omitted space.",
        "condition": residual_data["condition"],
        "boundary_blocks": {str(k): v for k, v in boundary.items()},
        "state_rows": state_rows,
        "summary": summary,
        "elapsed_seconds": time.perf_counter() - started,
    }
    out = RESULTS / "magnon_M3_boundary_denominator_N24_B1.500T.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": out.name, "summary": summary, "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
