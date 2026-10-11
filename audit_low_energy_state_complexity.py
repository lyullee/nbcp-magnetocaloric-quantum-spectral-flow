"""Audit structure of the low-energy states required by spectral flow."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_spectral_flow_entropy_n12 import lowest_global
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
GL_NODES, _ = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
REQUIRED_RANK = {9: 5, 12: 10}


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def schmidt_audit(state, n):
    cuts = []
    for cut in range(1, n):
        matrix = state.reshape(1 << (n - cut), 1 << cut)
        singular = np.linalg.svd(matrix, compute_uv=False)
        weight = singular**2
        weight /= weight.sum()
        positive = weight[weight > 1e-15]
        entropy = float(-np.dot(positive, np.log(positive)))
        cumulative = np.cumsum(weight)
        cuts.append({
            "cut_after_low_order_qubits": cut,
            "entropy_nats": entropy,
            "exact_schmidt_rank_1e_12": int(np.sum(weight > 1e-12)),
            "rank_for_99pct_weight": int(np.searchsorted(cumulative, 0.99) + 1),
            "rank_for_99_9pct_weight": int(np.searchsorted(cumulative, 0.999) + 1),
            "rank_for_99_99pct_weight": int(np.searchsorted(cumulative, 0.9999) + 1),
        })
    return cuts


def state_audit(state, n, energy, level):
    probability = np.abs(state)**2
    hamming = np.fromiter((index.bit_count() for index in range(1 << n)),
                          dtype=int, count=1 << n)
    q_probability = np.bincount(hamming, weights=probability,
                                minlength=n + 1)
    cumulative_order = np.argsort(probability)[::-1]
    cumulative_basis = np.cumsum(probability[cumulative_order])
    cuts = schmidt_audit(state, n)
    return {
        "level": level,
        "energy_K": float(energy),
        "dominant_parity": int(np.argmax([
            q_probability[::2].sum(), q_probability[1::2].sum()])),
        "up_spin_number_probability": q_probability.tolist(),
        "mean_up_spin_number": float(np.dot(np.arange(n + 1), q_probability)),
        "mean_defects_from_fully_up": float(
            n - np.dot(np.arange(n + 1), q_probability)),
        "up_spin_sectors_above_1e_6_probability": [
            int(q) for q, value in enumerate(q_probability) if value > 1e-6
        ],
        "largest_computational_basis_probability": float(probability.max()),
        "basis_states_for_99pct_weight": int(
            np.searchsorted(cumulative_basis, 0.99) + 1),
        "basis_states_for_99_9pct_weight": int(
            np.searchsorted(cumulative_basis, 0.999) + 1),
        "inverse_participation_effective_basis_states": float(
            1 / np.sum(probability**2)),
        "maximum_bipartite_entropy_nats": max(row["entropy_nats"]
                                                for row in cuts),
        "maximum_MPS_bond_rank_for_99pct_per_cut": max(
            row["rank_for_99pct_weight"] for row in cuts),
        "maximum_MPS_bond_rank_for_99_9pct_per_cut": max(
            row["rank_for_99_9pct_weight"] for row in cuts),
        "maximum_MPS_bond_rank_for_99_99pct_per_cut": max(
            row["rank_for_99_99pct_weight"] for row in cuts),
        "cuts": cuts,
    }


def solve_n9(graph, node, field_k, rank):
    n = graph["N"]
    h0 = dense_hamiltonian(n, model_bonds(graph, "M0", boundary="tin_foil"),
                           field_k)
    h3 = dense_hamiltonian(n, model_bonds(graph, "M3", boundary="tin_foil"),
                           field_k)
    energies, vectors = np.linalg.eigh(h0 + float(node) * (h3 - h0))
    return energies[:rank], vectors[:, :rank]


def solve_n12(graph, node, field_k, rank, seed):
    n = graph["N"]
    bonds0 = model_bonds(graph, "M0", boundary="tin_foil")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    h0 = {p: SpinHamiltonianOperator(n, bonds0, field_k, parity=p)
          for p in (0, 1)}
    h3 = {p: SpinHamiltonianOperator(n, bonds3, field_k, parity=p)
          for p in (0, 1)}
    energies, vectors = lowest_global(h0, h3, float(node), max(rank, 8),
                                      seed, True)
    return energies[:rank], vectors[:, :rank]


def main():
    started = time.perf_counter()
    rows = []
    clusters = {
        9: lattice((3, 0), (0, 3)),
        12: lattice((2, 2), (-2, 4)),
    }
    for n, graph in clusters.items():
        field_k = GZ * MU_B_OVER_KB * FIELD_T
        rank = REQUIRED_RANK[n]
        node_rows = []
        for node_index, node in enumerate(GL_NODES):
            if n == 9:
                energies, vectors = solve_n9(graph, node, field_k, rank)
            else:
                energies, vectors = solve_n12(
                    graph, node, field_k, rank,
                    20261012 + 100 * node_index)
            states = [state_audit(vectors[:, level], n, energies[level], level)
                      for level in range(rank)]
            node_rows.append({
                "lambda": float(node),
                "required_low_energy_states": rank,
                "states": states,
                "worst_state_summary": {
                    "maximum_MPS_bond_rank_for_99_9pct_per_cut": max(
                        state["maximum_MPS_bond_rank_for_99_9pct_per_cut"]
                        for state in states),
                    "maximum_basis_states_for_99_9pct_weight": max(
                        state["basis_states_for_99_9pct_weight"]
                        for state in states),
                    "maximum_bipartite_entropy_nats": max(
                        state["maximum_bipartite_entropy_nats"]
                        for state in states),
                },
            })
            print(f"N={n} lambda={node:.6f}: "
                  f"chi999={node_rows[-1]['worst_state_summary']['maximum_MPS_bond_rank_for_99_9pct_per_cut']}",
                  flush=True)
        rows.append({
            "N": n,
            "translations": graph["translations"],
            "nodes": node_rows,
        })

    payload = {
        "purpose": "state-structure audit for the low-energy eigenstates used by spectral-flow entropy reconstruction",
        "condition": {"field_T": FIELD_T,
                      "lambda_nodes": GL_NODES.tolist()},
        "clusters": rows,
        "scope_limits": [
            "Schmidt ranks use the current site ordering and can improve or worsen under qubit reordering.",
            "Per-cut retained weight does not directly equal whole-state MPS fidelity because truncation errors accumulate across cuts.",
            "Small MPS bond dimension indicates compressibility but is not an executable gate count.",
            "Exact classical eigenvectors are used only for screening and are unavailable in a scalable quantum calculation.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "low_energy_state_complexity_N9_N12.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "elapsed_seconds": payload["elapsed_seconds"],
        "summary": [{
            "N": cluster["N"],
            "nodes": [{
                "lambda": node["lambda"],
                **node["worst_state_summary"],
            } for node in cluster["nodes"]],
        } for cluster in rows],
    }, indent=2))


if __name__ == "__main__":
    main()
