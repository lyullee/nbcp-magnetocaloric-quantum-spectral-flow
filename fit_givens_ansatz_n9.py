"""Compress the N=9 central low state with a two-layer Givens ansatz."""
from __future__ import annotations

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
FIELD_T = 1.45
LAMBDA = 0.5
LAYERS = 2
RESTARTS = 5
MAXITER = 600
SEED = 20261015


def operations_for(n, edges):
    operations = []
    for layer in range(LAYERS):
        ordered = edges if layer % 2 == 0 else list(reversed(edges))
        operations.extend(("rz", site) for site in range(n))
        operations.extend(("givens", i, j) for i, j in ordered)
    return operations


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    edges = [(i, j) for i, j, _ in directional_edges(graph)]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, model_bonds(
        graph, "M0", boundary="tin_foil"), field_k)
    h3 = dense_hamiltonian(n, model_bonds(
        graph, "M3", boundary="tin_foil"), field_k)
    energies, vectors = np.linalg.eigh(h0 + LAMBDA * (h3 - h0))
    exact = vectors[:, 0]
    q_probability = np.zeros(n + 1)
    for index, value in enumerate(np.abs(exact)**2):
        q_probability[index.bit_count()] += value
    dominant_q = int(np.argmax(q_probability))
    target = exact.copy()
    target[np.array([index.bit_count() != dominant_q
                     for index in range(1 << n)])] = 0
    sector_weight = float(np.vdot(target, target).real)
    target /= np.linalg.norm(target)
    initial_index = int(np.argmax(np.abs(target)))
    initial = np.zeros(1 << n, dtype=complex)
    initial[initial_index] = 1
    operations = operations_for(n, edges)

    rng = np.random.default_rng(SEED)
    candidates = []
    for restart in range(RESTARTS):
        x0 = rng.normal(scale=0.15, size=len(operations))
        objective = FidelityObjective(initial, target, operations)
        result = minimize(
            objective.value_gradient, x0, jac=True, method="L-BFGS-B",
            bounds=[(-np.pi, np.pi)] * len(operations),
            options={"maxiter": MAXITER, "ftol": 1e-13, "gtol": 1e-9,
                     "maxls": 30})
        sector_fidelity = float(1 - result.fun)
        candidates.append({
            "parameters": result.x,
            "sector_fidelity": sector_fidelity,
            "full_state_fidelity": sector_weight * sector_fidelity,
            "success": bool(result.success),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_gradient_evaluations": objective.evaluations,
        })
        print(f"restart={restart}: sector={sector_fidelity:.12f}, "
              f"full={sector_weight*sector_fidelity:.12f}", flush=True)
    best = max(candidates, key=lambda row: row["sector_fidelity"])
    circuit = build_qiskit_circuit(n, initial_index, operations,
                                   best["parameters"])
    prepared = np.asarray(Statevector.from_instruction(circuit).data)
    qiskit_sector_fidelity = float(abs(np.vdot(target, prepared))**2)
    all_to_all = transpile(
        circuit, basis_gates=["rz", "sx", "x", "cx"],
        optimization_level=2, seed_transpiler=SEED)
    hardware = transpile(
        circuit, backend=FakeSherbrooke(), optimization_level=2,
        seed_transpiler=SEED)
    assert abs(qiskit_sector_fidelity - best["sector_fidelity"]) < 1e-8

    payload = {
        "purpose": "compressed N=9 central-state preparation for spectral flow",
        "condition": {"N": n, "translations": graph["translations"],
                      "field_T": FIELD_T, "lambda": LAMBDA,
                      "target_level": 0},
        "ansatz": "two layers of RZ plus number-conserving Givens gates on all 27 triangular bonds",
        "dominant_up_spin_sector": dominant_q,
        "dominant_sector_weight": sector_weight,
        "initial_computational_basis_index": initial_index,
        "parameter_count": len(operations),
        "two_qubit_parameterized_gate_count": sum(
            operation[0] == "givens" for operation in operations),
        "settings": {"layers": LAYERS, "restarts": RESTARTS,
                     "maxiter": MAXITER, "seed": SEED},
        "best": {
            "sector_fidelity": best["sector_fidelity"],
            "full_state_fidelity": best["full_state_fidelity"],
            "qiskit_sector_fidelity": qiskit_sector_fidelity,
            "success": best["success"],
            "message": best["message"],
            "iterations": best["iterations"],
            "function_gradient_evaluations":
                best["function_gradient_evaluations"],
            "parameters_radians": best["parameters"].tolist(),
        },
        "all_to_all": {
            "cx_count": int(all_to_all.count_ops().get("cx", 0)),
            "depth": int(all_to_all.depth()),
            "operation_count": int(all_to_all.size()),
        },
        "fake_sherbrooke": {
            "ecr_count": int(hardware.count_ops().get("ecr", 0)),
            "depth": int(hardware.depth()),
            "operation_count": int(hardware.size()),
        },
        "scope_limits": [
            "Only the ground state at lambda=0.5 is fitted; the complete algorithm needs several excited states at three nodes.",
            "The target amplitudes come from exact diagonalization, so this is an expressibility/resource test rather than scalable state discovery.",
            "The full-state fidelity ceiling is set by omitted weak mixing outside the dominant fixed-up-spin sector.",
            "Topology-aware compilation excludes hardware noise and optimization measurement cost.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "givens_ansatz_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "sector_weight": sector_weight,
        "best_sector_fidelity": best["sector_fidelity"],
        "best_full_state_fidelity": best["full_state_fidelity"],
        "all_to_all": payload["all_to_all"],
        "fake_sherbrooke": payload["fake_sherbrooke"],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
