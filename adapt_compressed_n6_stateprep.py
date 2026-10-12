"""Build adaptive compact N=6 circuits for the full spectral-flow pilot.

Starting from one physics-informed local layer, append only those second-layer
generators with the largest fidelity gradient.  The output stores every circuit
parameter so a read-only live-backend preflight and, if justified, one guarded
Open Plan job can reproduce the exact payload.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from fit_number_conserving_ansatz_n9 import (
    FidelityObjective,
    build_qiskit_circuit,
    dense_hamiltonian,
)
from qiskit_pauli_mapping import (
    open_triangular_directional_edges,
    open_weak_model_bonds,
)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ
from spectral_flow_circuit_unit_n6 import GL_NODES, GL_WEIGHTS, entropy


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.05
RANK = 3
SEED = 20261013
MAXITER = 1000
FIDELITY_STOP = 0.9999
V_ERROR_STOP_K = 1.0e-3


def layer_operations(n, edges, reverse=False):
    ordered = list(reversed(edges)) if reverse else list(edges)
    operations = [("rz", site) for site in range(n)]
    operations.extend(("givens", i, j) for i, j in ordered)
    operations.extend(("rz", site) for site in range(n))
    operations.extend(("pair", i, j) for i, j in ordered)
    return operations


def optimize(initial, target, operations, starts):
    rows = []
    for x0 in starts:
        objective = FidelityObjective(initial, target, operations)
        result = minimize(
            objective.value_gradient,
            x0,
            jac=True,
            method="L-BFGS-B",
            bounds=[(-np.pi, np.pi)] * len(operations),
            options={"maxiter": MAXITER, "ftol": 1e-14, "gtol": 1e-10,
                     "maxls": 50},
        )
        rows.append({
            "parameters": result.x,
            "fidelity": float(1 - result.fun),
            "iterations": int(result.nit),
            "evaluations": int(objective.evaluations),
            "success": bool(result.success),
        })
    return max(rows, key=lambda row: row["fidelity"])


def prepared_state(n, initial_index, operations, parameters):
    circuit = build_qiskit_circuit(n, initial_index, operations, parameters)
    state = np.asarray(Statevector.from_instruction(circuit).data)
    return circuit, state


def select_gradient_generator(initial, target, operations, parameters, pool):
    scores = []
    for index, candidate in enumerate(pool):
        trial_operations = operations + [candidate]
        trial_parameters = np.r_[parameters, 0.0]
        objective = FidelityObjective(initial, target, trial_operations)
        _, gradient = objective.value_gradient(trial_parameters)
        scores.append((abs(float(gradient[-1])), index, candidate))
    return max(scores, key=lambda row: (row[0], -row[1]))


def fit_adaptive(target, v, base_operations, candidate_pool, rng):
    n = int(np.log2(target.size))
    initial_index = int(np.argmax(np.abs(target)))
    initial = np.zeros(target.size, dtype=complex)
    initial[initial_index] = 1
    starts = [rng.normal(scale=scale, size=len(base_operations))
              for scale in (0.03, 0.10, 0.20, 0.35)]
    best = optimize(initial, target, base_operations, starts)
    operations = list(base_operations)
    parameters = best["parameters"].copy()
    pool = list(candidate_pool)
    history = []

    for step in range(len(candidate_pool) + 1):
        circuit, state = prepared_state(
            n, initial_index, operations, parameters)
        fidelity = float(abs(np.vdot(target, state)) ** 2)
        exact_v = float(np.real(np.vdot(target, v @ target)))
        estimated_v = float(np.real(np.vdot(state, v @ state)))
        history.append({
            "step": step,
            "operation_count": len(operations),
            "fidelity": fidelity,
            "V_error_K": estimated_v - exact_v,
            "added_operation": (None if step == 0 else list(operations[-1])),
        })
        if fidelity >= FIDELITY_STOP and abs(estimated_v - exact_v) <= V_ERROR_STOP_K:
            break
        if not pool:
            break
        score, selected_index, selected = select_gradient_generator(
            initial, target, operations, parameters, pool)
        operations.append(selected)
        pool.pop(selected_index)
        warm = np.r_[parameters, 0.0]
        perturb = warm + rng.normal(scale=0.01, size=warm.size)
        best = optimize(initial, target, operations, (warm, perturb))
        parameters = best["parameters"].copy()

    circuit, state = prepared_state(n, initial_index, operations, parameters)
    return {
        "initial_index": initial_index,
        "operations": operations,
        "parameters": parameters,
        "circuit": circuit,
        "state": state,
        "history": history,
    }


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges()
    n = graph["N"]
    edges = [(i, j) for i, j, _ in graph["directional_edges"]]
    base_operations = layer_operations(n, edges, reverse=False)
    candidate_pool = layer_operations(n, edges, reverse=True)
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    v = h2 - h0
    e0 = np.linalg.eigvalsh(h0)
    e2 = np.linalg.eigvalsh(h2)
    exact_delta_s = (entropy(e2) - entropy(e0)) / n
    fake = FakeSherbrooke()
    rng = np.random.default_rng(SEED)
    integrated_v = np.zeros(RANK)
    node_rows = []

    for node, weight in zip(GL_NODES, GL_WEIGHTS, strict=True):
        energies, vectors = np.linalg.eigh(h0 + float(node) * v)
        state_rows = []
        for level in range(RANK):
            target = vectors[:, level]
            fitted = fit_adaptive(
                target, v, base_operations, candidate_pool, rng)
            circuit = fitted["circuit"]
            state = fitted["state"]
            exact_v = float(np.real(np.vdot(target, v @ target)))
            prepared_v = float(np.real(np.vdot(state, v @ state)))
            integrated_v[level] += float(weight) * prepared_v
            all_to_all = transpile(
                circuit, basis_gates=["rz", "sx", "x", "cx"],
                optimization_level=3, seed_transpiler=SEED)
            hardware = transpile(
                circuit, backend=fake, optimization_level=3,
                seed_transpiler=SEED)
            row = {
                "lambda": float(node),
                "weight": float(weight),
                "level": level,
                "initial_basis_index": fitted["initial_index"],
                "operations": [list(item) for item in fitted["operations"]],
                "parameters_radians": fitted["parameters"].tolist(),
                "adaptive_history": fitted["history"],
                "fidelity": float(abs(np.vdot(target, state)) ** 2),
                "exact_V_expectation_K": exact_v,
                "prepared_V_expectation_K": prepared_v,
                "V_error_K": prepared_v - exact_v,
                "all_to_all": {
                    "cx_count": int(all_to_all.count_ops().get("cx", 0)),
                    "depth": int(all_to_all.depth()),
                },
                "fake_sherbrooke": {
                    "ecr_count": int(hardware.count_ops().get("ecr", 0)),
                    "depth": int(hardware.depth()),
                },
            }
            state_rows.append(row)
            print(f"lambda={node:.6f} level={level} "
                  f"added={len(fitted['operations'])-len(base_operations)} "
                  f"F={row['fidelity']:.8f} dV={row['V_error_K']:+.3e} "
                  f"CX={row['all_to_all']['cx_count']} "
                  f"ECR={row['fake_sherbrooke']['ecr_count']}", flush=True)
        node_rows.append({"lambda": float(node), "weight": float(weight),
                          "states": state_rows})

    reconstructed = e0[:RANK] + integrated_v
    reconstructed_delta_s = (entropy(reconstructed) -
                             entropy(e0[:RANK])) / n
    flat = [row for node in node_rows for row in node["states"]]
    payload = {
        "purpose": "adaptive compact N=6 circuits for a complete spectral-flow QPU pilot",
        "condition": {"N": n, "field_T": FIELD_T,
                      "temperature_K": 0.10,
                      "path": "M0 + lambda*(M2-M0)",
                      "low_energy_rank": RANK},
        "ansatz": "one local Givens/pair layer plus gradient-selected second-layer generators",
        "stop": {"minimum_fidelity": FIDELITY_STOP,
                 "maximum_absolute_V_error_K": V_ERROR_STOP_K},
        "exact_delta_S_kB_per_spin": exact_delta_s,
        "reconstructed_delta_S_kB_per_spin": reconstructed_delta_s,
        "delta_S_error_kB_per_spin": reconstructed_delta_s - exact_delta_s,
        "passes_entropy_target": bool(
            abs(reconstructed_delta_s - exact_delta_s) <= 0.001),
        "minimum_fidelity": min(row["fidelity"] for row in flat),
        "maximum_absolute_V_error_K": max(abs(row["V_error_K"])
                                           for row in flat),
        "maximum_all_to_all_cx": max(row["all_to_all"]["cx_count"]
                                      for row in flat),
        "maximum_fake_sherbrooke_ecr": max(
            row["fake_sherbrooke"]["ecr_count"] for row in flat),
        "maximum_added_generators": max(
            len(row["operations"]) - len(base_operations) for row in flat),
        "nodes": node_rows,
        "elapsed_seconds": time.perf_counter() - started,
        "scope_limits": [
            "Exact eigenvectors are fitting targets; the experiment audits execution, not scalable eigenstate discovery.",
            "The open N=6 patch is a circuit benchmark and does not determine the NBCP material conclusion.",
            "No IBM service is contacted by this script.",
        ],
    }
    output = ROOT / "results" / "adaptive_compressed_stateprep_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({key: payload[key] for key in (
        "passes_entropy_target", "minimum_fidelity",
        "maximum_absolute_V_error_K", "maximum_all_to_all_cx",
        "maximum_fake_sherbrooke_ecr", "maximum_added_generators",
        "delta_S_error_kB_per_spin", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
