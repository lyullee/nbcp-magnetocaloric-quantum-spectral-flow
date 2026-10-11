"""Fit all symmetry-distinct N=9 states needed for spectral flow."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from fit_number_conserving_ansatz_n9 import (FidelityObjective,
                                             build_qiskit_circuit,
                                             dense_hamiltonian)
from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ, directional_edges


ROOT = Path(__file__).resolve().parent
DEFAULT_FIELD_T = 1.45
TEMPERATURES_K = (0.09, 0.10, 0.12)
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
GL_WEIGHTS = GL_WEIGHTS / 2
SECTOR_RESTARTS = 2
FULL_RESTARTS = 2
MAXITER = 700
SEED = 20261017


def entropy(energies, temperature):
    probability = np.exp(-(energies - energies[0]) / temperature)
    probability /= probability.sum()
    return float(-np.dot(probability[probability > 0],
                         np.log(probability[probability > 0])))


def optimize(initial, target, operations, starts):
    candidates = []
    for x0 in starts:
        objective = FidelityObjective(initial, target, operations)
        result = minimize(
            objective.value_gradient, x0, jac=True, method="L-BFGS-B",
            bounds=[(-np.pi, np.pi)] * len(operations),
            options={"maxiter": MAXITER, "ftol": 1e-14, "gtol": 1e-10,
                     "maxls": 40})
        candidates.append({
            "parameters": result.x,
            "fidelity": float(1 - result.fun),
            "success": bool(result.success),
            "message": str(result.message),
            "iterations": int(result.nit),
            "evaluations": objective.evaluations,
        })
    return max(candidates, key=lambda row: row["fidelity"])


def target_specification(field_t):
    source = json.loads((ROOT / "results" /
                         "spectral_flow_entropy_N9.json").read_text(
                             encoding="utf-8"))
    rows = [row for row in source["rows"]
            if abs(row["field_T"] - field_t) < 1e-9 and
            any(abs(row["temperature_K"] - temperature) < 1e-9
                for temperature in TEMPERATURES_K)]
    if not rows:
        raise ValueError(f"No frozen N=9 target rows for B={field_t} T")
    rank = max(int(row["low_energy_rank"]) for row in rows)
    targets = {}
    for temperature in TEMPERATURES_K:
        candidates = [row for row in rows
                      if abs(row["temperature_K"] - temperature) < 1e-9]
        targets[temperature] = min(
            float(row["target_delta_S_error_kB_per_spin"])
            for row in candidates)
    return rank, targets


def degenerate_groups(energies, derivatives, rank, tolerance=1e-10):
    groups = []
    start = 0
    while start < rank:
        stop = start + 1
        while (stop < rank and
               abs(energies[stop] - energies[start]) < tolerance and
               abs(derivatives[stop] - derivatives[start]) < tolerance):
            stop += 1
        groups.append(tuple(range(start, stop)))
        start = stop
    return groups


def main(field_t=DEFAULT_FIELD_T):
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    edges = [(i, j) for i, j, _ in directional_edges(graph)]
    low_energy_rank, target_errors = target_specification(field_t)
    field_k = GZ * MU_B_OVER_KB * field_t
    h0 = dense_hamiltonian(n, model_bonds(
        graph, "M0", boundary="tin_foil"), field_k)
    h3 = dense_hamiltonian(n, model_bonds(
        graph, "M3", boundary="tin_foil"), field_k)
    v = h3 - h0
    exact_e0 = np.linalg.eigvalsh(h0)
    exact_e3 = np.linalg.eigvalsh(h3)

    sector_operations = []
    for layer in range(2):
        ordered = edges if layer == 0 else list(reversed(edges))
        sector_operations.extend(("rz", site) for site in range(n))
        sector_operations.extend(("givens", i, j) for i, j in ordered)
    full_operations = list(sector_operations)
    full_operations.extend(("rz", site) for site in range(n))
    full_operations.extend(("pair", i, j) for i, j in edges)
    run_seed = SEED + int(round((field_t - DEFAULT_FIELD_T) * 1000))
    rng = np.random.default_rng(run_seed)
    backend = FakeSherbrooke()
    node_rows = []
    prepared_shifts = np.zeros(low_energy_rank)
    exact_shifts = np.zeros(low_energy_rank)
    group_patterns = []

    for node_index, (node, weight) in enumerate(zip(
            GL_NODES, GL_WEIGHTS, strict=True)):
        energies, vectors = np.linalg.eigh(h0 + float(node) * v)
        exact_v = np.real(np.sum(vectors[:, :low_energy_rank].conj() *
                                 (v @ vectors[:, :low_energy_rank]), axis=0))
        groups = degenerate_groups(energies, exact_v, low_energy_rank)
        group_patterns.append([list(group) for group in groups])
        exact_shifts += float(weight) * exact_v
        representative_rows = []
        approximate_v = np.zeros(low_energy_rank)

        for group in groups:
            representative = group[0]
            target = vectors[:, representative]
            q_probability = np.zeros(n + 1)
            for index, value in enumerate(np.abs(target)**2):
                q_probability[index.bit_count()] += value
            dominant_q = int(np.argmax(q_probability))
            sector_target = target.copy()
            sector_target[np.array([index.bit_count() != dominant_q
                                    for index in range(1 << n)])] = 0
            sector_weight = float(np.vdot(sector_target,
                                          sector_target).real)
            sector_target /= np.linalg.norm(sector_target)
            initial_index = int(np.argmax(np.abs(sector_target)))
            initial = np.zeros(1 << n, dtype=complex)
            initial[initial_index] = 1

            sector_starts = [rng.normal(scale=0.15,
                                        size=len(sector_operations))
                             for _ in range(SECTOR_RESTARTS)]
            sector_best = optimize(initial, sector_target,
                                   sector_operations, sector_starts)
            full_starts = []
            for restart in range(FULL_RESTARTS):
                x0 = np.zeros(len(full_operations))
                x0[:len(sector_operations)] = sector_best["parameters"]
                if restart:
                    x0[:len(sector_operations)] += rng.normal(
                        scale=0.015, size=len(sector_operations))
                x0[len(sector_operations):] = rng.normal(
                    scale=0.015 if restart == 0 else 0.04,
                    size=len(full_operations) - len(sector_operations))
                full_starts.append(x0)
            full_best = optimize(initial, target, full_operations, full_starts)
            circuit = build_qiskit_circuit(n, initial_index, full_operations,
                                           full_best["parameters"])
            prepared = np.asarray(Statevector.from_instruction(circuit).data)
            qiskit_fidelity = float(abs(np.vdot(target, prepared))**2)
            prepared_v = float(np.real(np.vdot(prepared, v @ prepared)))
            all_to_all = transpile(
                circuit, basis_gates=["rz", "sx", "x", "cx"],
                optimization_level=2, seed_transpiler=SEED)
            hardware = transpile(
                circuit, backend=backend, optimization_level=2,
                seed_transpiler=SEED)
            assert abs(qiskit_fidelity - full_best["fidelity"]) < 1e-8
            for level in group:
                approximate_v[level] = prepared_v
            representative_rows.append({
                "representative_level": representative,
                "multiplicity": len(group),
                "represented_levels": list(group),
                "dominant_up_spin_sector": dominant_q,
                "dominant_sector_weight": sector_weight,
                "sector_fit_fidelity": sector_best["fidelity"],
                "full_fit_fidelity": full_best["fidelity"],
                "qiskit_fidelity": qiskit_fidelity,
                "exact_V_expectation_K": float(exact_v[representative]),
                "prepared_V_expectation_K": prepared_v,
                "V_expectation_error_K": prepared_v -
                    float(exact_v[representative]),
                "initial_computational_basis_index": initial_index,
                "optimizer": {
                    "sector_iterations": sector_best["iterations"],
                    "sector_evaluations": sector_best["evaluations"],
                    "full_iterations": full_best["iterations"],
                    "full_evaluations": full_best["evaluations"],
                },
                "all_to_all": {
                    "cx_count": int(all_to_all.count_ops().get("cx", 0)),
                    "depth": int(all_to_all.depth()),
                },
                "fake_sherbrooke": {
                    "ecr_count": int(hardware.count_ops().get("ecr", 0)),
                    "depth": int(hardware.depth()),
                },
                "parameters_radians": full_best["parameters"].tolist(),
            })
            print(f"lambda={node:.6f} L={representative}: "
                  f"F={full_best['fidelity']:.9f}, "
                  f"dV={prepared_v-exact_v[representative]:+.3e} K",
                  flush=True)
        prepared_shifts += float(weight) * approximate_v
        node_rows.append({
            "lambda": float(node),
            "weight": float(weight),
            "exact_V_expectation_K_by_level": exact_v.tolist(),
            "prepared_V_expectation_K_by_level": approximate_v.tolist(),
            "representatives": representative_rows,
        })

    exact_reconstructed = exact_e0[:low_energy_rank] + exact_shifts
    circuit_reconstructed = exact_e0[:low_energy_rank] + prepared_shifts
    thermodynamic_rows = []
    for temperature in TEMPERATURES_K:
        exact_delta_s = (entropy(exact_e3, temperature) -
                         entropy(exact_e0, temperature)) / n
        exact_flow_delta_s = (entropy(exact_reconstructed, temperature) -
                              entropy(exact_e0[:low_energy_rank], temperature)) / n
        circuit_delta_s = (entropy(circuit_reconstructed, temperature) -
                           entropy(exact_e0[:low_energy_rank], temperature)) / n
        target_error = target_errors[temperature]
        thermodynamic_rows.append({
            "temperature_K": temperature,
            "target_error_kB_per_spin": target_error,
            "exact_full_delta_S_kB_per_spin": exact_delta_s,
            "exact_low_rank_flow_delta_S_kB_per_spin": exact_flow_delta_s,
            "circuit_prepared_delta_S_kB_per_spin": circuit_delta_s,
            "state_preparation_incremental_error_kB_per_spin":
                circuit_delta_s - exact_flow_delta_s,
            "total_error_kB_per_spin": circuit_delta_s - exact_delta_s,
            "sign_matches_exact": bool(np.sign(circuit_delta_s) ==
                                       np.sign(exact_delta_s)),
            "passes_target": bool(abs(circuit_delta_s - exact_delta_s) <=
                                  target_error),
        })

    payload = {
        "purpose": "end-to-end state-preparation bias test for all symmetry-distinct N=9 spectral-flow states",
        "condition": {"N": n, "translations": graph["translations"],
                      "field_T": field_t,
                      "lambda_nodes": GL_NODES.tolist(),
                      "low_energy_levels": low_energy_rank,
                      "degenerate_level_groups_by_node": group_patterns,
                      "unique_state_circuits_by_node":
                          [len(pattern) for pattern in group_patterns]},
        "ansatz": "two triangular Givens layers plus one triangular pair-mixing layer",
        "settings": {"sector_restarts": SECTOR_RESTARTS,
                     "full_restarts": FULL_RESTARTS,
                     "maxiter": MAXITER, "seed": run_seed},
        "nodes": node_rows,
        "thermodynamics": thermodynamic_rows,
        "scope_limits": [
            "Exact eigenvectors define variational targets; this tests circuit expressibility and bias, not scalable state discovery.",
            "One circuit represents each exactly degenerate doublet because its energy and V expectation are symmetry-identical.",
            "All-to-all and FakeSherbrooke resource counts exclude optimizer measurements, device noise and mitigation.",
            f"This file tests only B={field_t:.2f} T; other fields require their own end-to-end audit.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = (ROOT / "results" /
              f"spectral_flow_stateprep_N9_B{field_t:.2f}T.json")
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "minimum_fidelity": min(rep["full_fit_fidelity"]
                                for node in node_rows
                                for rep in node["representatives"]),
        "maximum_abs_V_error_K": max(abs(rep["V_expectation_error_K"])
                                     for node in node_rows
                                     for rep in node["representatives"]),
        "maximum_fake_sherbrooke_ECR": max(rep["fake_sherbrooke"]["ecr_count"]
                                           for node in node_rows
                                           for rep in node["representatives"]),
        "thermodynamics": thermodynamic_rows,
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=float, default=DEFAULT_FIELD_T,
                        choices=(1.40, 1.45, 1.50))
    main(parser.parse_args().field)
