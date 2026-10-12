"""Screen compact physics-informed N=6 state-preparation circuits.

The existing live-QPU pilot used exact amplitude loading.  This script tests
whether the same nine low-energy spectral-flow states can be represented by
short local Givens/pair-mixing circuits before any additional Open Plan time
is consumed.  It is deliberately local-only and never connects to IBM.
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
SEED = 20261012
RESTARTS = 5
MAXITER = 1200


def operation_structure(n, edges, layers):
    operations = []
    for layer in range(layers):
        ordered = edges if layer % 2 == 0 else list(reversed(edges))
        operations.extend(("rz", site) for site in range(n))
        operations.extend(("givens", i, j) for i, j in ordered)
        operations.extend(("rz", site) for site in range(n))
        operations.extend(("pair", i, j) for i, j in ordered)
    return operations


def fit_state(target, operations, rng, warm=None):
    initial_index = int(np.argmax(np.abs(target)))
    initial = np.zeros(target.size, dtype=complex)
    initial[initial_index] = 1
    candidates = []
    for restart in range(RESTARTS):
        if restart == 0 and warm is not None:
            x0 = np.zeros(len(operations))
            x0[: min(len(warm), len(x0))] = warm[: min(len(warm), len(x0))]
            if len(x0) > len(warm):
                x0[len(warm):] = rng.normal(scale=0.02,
                                             size=len(x0) - len(warm))
        else:
            x0 = rng.normal(scale=0.15 if restart else 0.05,
                            size=len(operations))
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
        candidates.append({
            "parameters": result.x,
            "fidelity": float(1 - result.fun),
            "iterations": int(result.nit),
            "evaluations": int(objective.evaluations),
            "success": bool(result.success),
            "message": str(result.message),
        })
    best = max(candidates, key=lambda item: item["fidelity"])
    circuit = build_qiskit_circuit(
        int(np.log2(target.size)), initial_index, operations,
        best["parameters"])
    state = np.asarray(Statevector.from_instruction(circuit).data)
    qiskit_fidelity = float(abs(np.vdot(target, state)) ** 2)
    if abs(qiskit_fidelity - best["fidelity"]) > 1e-8:
        raise RuntimeError("custom/Qiskit ansatz mismatch")
    return best, circuit, state, initial_index


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges()
    n = graph["N"]
    edges = [(i, j) for i, j, _ in graph["directional_edges"]]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    v = h2 - h0
    e0 = np.linalg.eigvalsh(h0)
    e2 = np.linalg.eigvalsh(h2)
    rank = 3
    exact_delta_s = (entropy(e2) - entropy(e0)) / n
    rng = np.random.default_rng(SEED)
    fake = FakeSherbrooke()
    layer_results = []
    previous = {}

    for layers in (1, 2, 3):
        operations = operation_structure(n, edges, layers)
        node_rows = []
        integrated_v = np.zeros(rank)
        minimum_fidelity = 1.0
        maximum_v_error = 0.0
        maximum_cx = 0
        maximum_ecr = 0
        for node_index, (node, weight) in enumerate(zip(
                GL_NODES, GL_WEIGHTS, strict=True)):
            energies, vectors = np.linalg.eigh(h0 + float(node) * v)
            rows = []
            for level in range(rank):
                key = (node_index, level)
                target = vectors[:, level]
                best, circuit, state, initial_index = fit_state(
                    target, operations, rng, previous.get(key))
                previous[key] = best["parameters"].copy()
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
                    "level": level,
                    "initial_basis_index": initial_index,
                    "fidelity": best["fidelity"],
                    "exact_V_expectation_K": exact_v,
                    "prepared_V_expectation_K": prepared_v,
                    "V_error_K": prepared_v - exact_v,
                    "parameter_count": len(operations),
                    "all_to_all_cx": int(all_to_all.count_ops().get("cx", 0)),
                    "all_to_all_depth": int(all_to_all.depth()),
                    "fake_sherbrooke_ecr": int(
                        hardware.count_ops().get("ecr", 0)),
                    "fake_sherbrooke_depth": int(hardware.depth()),
                    "optimizer_iterations": best["iterations"],
                    "optimizer_success": best["success"],
                }
                rows.append(row)
                minimum_fidelity = min(minimum_fidelity, best["fidelity"])
                maximum_v_error = max(maximum_v_error,
                                      abs(prepared_v - exact_v))
                maximum_cx = max(maximum_cx, row["all_to_all_cx"])
                maximum_ecr = max(maximum_ecr, row["fake_sherbrooke_ecr"])
                print(f"L={layers} lambda={node:.6f} level={level} "
                      f"F={best['fidelity']:.8f} "
                      f"dV={prepared_v-exact_v:+.3e} "
                      f"CX={row['all_to_all_cx']} ECR={row['fake_sherbrooke_ecr']}",
                      flush=True)
            node_rows.append({"lambda": float(node), "states": rows})

        reconstructed = e0[:rank] + integrated_v
        reconstructed_delta_s = (
            entropy(reconstructed) - entropy(e0[:rank])) / n
        layer_results.append({
            "layers": layers,
            "parameter_count": len(operations),
            "minimum_fidelity": minimum_fidelity,
            "maximum_V_error_K": maximum_v_error,
            "maximum_all_to_all_cx": maximum_cx,
            "maximum_fake_sherbrooke_ecr": maximum_ecr,
            "reconstructed_delta_S_kB_per_spin": reconstructed_delta_s,
            "exact_delta_S_kB_per_spin": exact_delta_s,
            "delta_S_error_kB_per_spin": reconstructed_delta_s - exact_delta_s,
            "passes_entropy_target": bool(
                abs(reconstructed_delta_s - exact_delta_s) <= 0.001),
            "nodes": node_rows,
        })

    payload = {
        "purpose": "local-only compact N=6 state-preparation screen before additional QPU use",
        "condition": {"N": n, "field_T": FIELD_T,
                      "path": "M0 + lambda*(M2-M0)",
                      "low_energy_rank": rank},
        "ansatz": "local RZ, Givens and pair-mixing layers on the open triangular graph",
        "settings": {"restarts": RESTARTS, "maxiter": MAXITER,
                     "seed": SEED},
        "layers": layer_results,
        "recommended_layer": next((
            row["layers"] for row in layer_results
            if row["passes_entropy_target"] and
            row["maximum_fake_sherbrooke_ecr"] < 129), None),
        "elapsed_seconds": time.perf_counter() - started,
        "scope_limits": [
            "Exact eigenvectors are fitting targets; this is circuit compression, not scalable state discovery.",
            "Fake-backend compilation does not establish live-hardware accuracy.",
            "No IBM service is contacted and no QPU time is consumed.",
        ],
    }
    output = ROOT / "results" / "compressed_stateprep_N6_screen.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "recommended_layer": payload["recommended_layer"],
        "layers": [{key: row[key] for key in (
            "layers", "minimum_fidelity", "maximum_V_error_K",
            "maximum_all_to_all_cx", "maximum_fake_sherbrooke_ecr",
            "delta_S_error_kB_per_spin", "passes_entropy_target")}
                   for row in layer_results],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
