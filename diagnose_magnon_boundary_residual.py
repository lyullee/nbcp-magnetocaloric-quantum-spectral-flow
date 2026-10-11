"""Measure coupling from a truncated M3 eigenspace to the first omitted q block.

This is a diagnostic, not a rigorous free-energy error bound.  For a parity block
with cutoff K, M3 changes q by 0 or +/-2, so embedding an eigenvector into K+2
leaves its residual entirely in the newly added q=K+2 block.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ
from truncated_magnon_operator import TruncatedMagnonCartesianOperator

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CELLS = {
    9: ((3, 0), (0, 3)),
    24: ((-5, 1), (1, -5)),
}


def block_diagnostic(n, bonds, field_k, parity, cutoff, k, started):
    small = TruncatedMagnonCartesianOperator(n, cutoff, parity, bonds, field_k)
    count = min(k, small.shape[0] - 1)
    rng = np.random.default_rng(20261011 + parity)
    v0 = rng.normal(size=small.shape[0]) + 1j * rng.normal(size=small.shape[0])
    values, vectors = eigsh(small, k=count, which="SA", tol=1e-9, v0=v0)
    order = np.argsort(values.real)
    values, vectors = values.real[order], vectors[:, order]
    print(
        f"N={n} parity={parity} K={cutoff} dim={small.shape[0]} "
        f"E0={values[0]:.9f} eigenpairs elapsed={time.perf_counter()-started:.1f}s",
        flush=True,
    )

    expanded = TruncatedMagnonCartesianOperator(n, cutoff + 2, parity, bonds, field_k)
    if not np.array_equal(expanded.states[: small.shape[0]], small.states):
        raise RuntimeError("expanded basis does not preserve the smaller basis as a prefix")
    rows = []
    for j, energy in enumerate(values):
        embedded = np.zeros(expanded.shape[0], dtype=np.complex128)
        embedded[: small.shape[0]] = vectors[:, j]
        residual = expanded @ embedded - energy * embedded
        kept_norm = float(np.linalg.norm(residual[: small.shape[0]]))
        omitted_norm = float(np.linalg.norm(residual[small.shape[0] :]))
        rows.append(
            {
                "state_index": j,
                "energy_K": float(energy),
                "kept_space_residual_norm_K": kept_norm,
                "first_omitted_block_residual_norm_K": omitted_norm,
            }
        )
    print(
        f"N={n} parity={parity} added_q={cutoff+2} dim={expanded.shape[0]-small.shape[0]} "
        f"boundary residual range={min(r['first_omitted_block_residual_norm_K'] for r in rows):.3e}"
        f"..{max(r['first_omitted_block_residual_norm_K'] for r in rows):.3e} K "
        f"elapsed={time.perf_counter()-started:.1f}s",
        flush=True,
    )
    return {
        "parity": parity,
        "included_max_q": cutoff,
        "first_omitted_q": cutoff + 2,
        "included_dimension": small.shape[0],
        "first_omitted_dimension": expanded.shape[0] - small.shape[0],
        "rows": rows,
    }


def main(n: int, k: int):
    if n not in CELLS:
        raise ValueError(f"supported N: {sorted(CELLS)}")
    a, b = CELLS[n]
    graph = lattice(a, b)
    bonds = model_bonds(graph, "M3")
    field_t = 1.50
    temperature = 0.09
    field_k = GZ * MU_B_OVER_KB * field_t
    cutoffs = {0: 6, 1: 7} if n == 24 else {0: 2, 1: 3}
    started = time.perf_counter()
    blocks = [
        block_diagnostic(n, bonds, field_k, parity, cutoffs[parity], k, started)
        for parity in (0, 1)
    ]

    all_rows = [dict(row, parity=block["parity"]) for block in blocks for row in block["rows"]]
    energies = np.array([row["energy_K"] for row in all_rows])
    residuals = np.array([row["first_omitted_block_residual_norm_K"] for row in all_rows])
    weights = np.exp(-(energies - energies.min()) / temperature)
    weights /= weights.sum()
    diagnostic = {
        "boltzmann_weighted_mean_boundary_residual_K": float(np.dot(weights, residuals)),
        "boltzmann_weighted_rms_boundary_residual_K": float(np.sqrt(np.dot(weights, residuals**2))),
        "maximum_boundary_residual_K": float(residuals.max()),
        "weighted_rms_over_temperature": float(np.sqrt(np.dot(weights, residuals**2)) / temperature),
        "weighted_rms_over_N_temperature": float(np.sqrt(np.dot(weights, residuals**2)) / (n * temperature)),
    }
    payload = {
        "method": "Ritz residual after exact embedding from q<=K to the first omitted same-parity q=K+2 block",
        "warning": "Residual norms diagnose omitted coupling but are not a rigorous entropy or free-energy error bound.",
        "condition": {"N": n, "translations": [a, b], "model": "M3", "field_T": field_t, "temperature_K": temperature},
        "settings": {"low_states_per_parity": k, "eigsh_tolerance": 1e-9, "seed_base": 20261011},
        "blocks": blocks,
        "diagnostic": diagnostic,
        "elapsed_seconds": time.perf_counter() - started,
    }
    out = RESULTS / f"magnon_M3_boundary_residual_N{n}_B1.500T.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": out.name, "diagnostic": diagnostic, "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=24)
    parser.add_argument("--k", type=int, default=8)
    args = parser.parse_args()
    main(args.n, args.k)
