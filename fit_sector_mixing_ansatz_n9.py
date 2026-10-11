"""Add one parity-preserving pair-mixing layer to the N=9 Givens state."""
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
RESTARTS = 5
MAXITER = 800
SEED = 20261016


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
    _, vectors = np.linalg.eigh(h0 + LAMBDA * (h3 - h0))
    target = vectors[:, 0]
    dominant = json.loads((ROOT / "results" / "givens_ansatz_N9.json").read_text(
        encoding="utf-8"))
    initial_index = dominant["initial_computational_basis_index"]
    initial = np.zeros(1 << n, dtype=complex)
    initial[initial_index] = 1

    operations = []
    for layer in range(2):
        ordered = edges if layer == 0 else list(reversed(edges))
        operations.extend(("rz", site) for site in range(n))
        operations.extend(("givens", i, j) for i, j in ordered)
    # A final local-phase plus pair-rotation layer captures the weak q±2
    # admixture created by J±± and dipolar transverse anisotropy.
    operations.extend(("rz", site) for site in range(n))
    operations.extend(("pair", i, j) for i, j in edges)
    dominant_parameters = np.asarray(
        dominant["best"]["parameters_radians"], dtype=float)
    assert len(dominant_parameters) == 72

    rng = np.random.default_rng(SEED)
    candidates = []
    for restart in range(RESTARTS):
        x0 = np.zeros(len(operations))
        x0[:len(dominant_parameters)] = dominant_parameters
        if restart:
            x0[:len(dominant_parameters)] += rng.normal(
                scale=0.02, size=len(dominant_parameters))
        x0[len(dominant_parameters):] = rng.normal(
            scale=0.01 if restart == 0 else 0.04,
            size=len(operations) - len(dominant_parameters))
        objective = FidelityObjective(initial, target, operations)
        result = minimize(
            objective.value_gradient, x0, jac=True, method="L-BFGS-B",
            bounds=[(-np.pi, np.pi)] * len(operations),
            options={"maxiter": MAXITER, "ftol": 1e-14, "gtol": 1e-10,
                     "maxls": 40})
        fidelity = float(1 - result.fun)
        candidates.append({
            "parameters": result.x,
            "fidelity": fidelity,
            "success": bool(result.success),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_gradient_evaluations": objective.evaluations,
        })
        print(f"restart={restart}: fidelity={fidelity:.12f}, "
              f"iter={result.nit}", flush=True)
    best = max(candidates, key=lambda row: row["fidelity"])
    circuit = build_qiskit_circuit(n, initial_index, operations,
                                   best["parameters"])
    prepared = np.asarray(Statevector.from_instruction(circuit).data)
    qiskit_fidelity = float(abs(np.vdot(target, prepared))**2)
    v = h3 - h0
    exact_v = float(np.real(np.vdot(target, v @ target)))
    prepared_v = float(np.real(np.vdot(prepared, v @ prepared)))
    all_to_all = transpile(
        circuit, basis_gates=["rz", "sx", "x", "cx"],
        optimization_level=2, seed_transpiler=SEED)
    hardware = transpile(
        circuit, backend=FakeSherbrooke(), optimization_level=2,
        seed_transpiler=SEED)
    assert abs(qiskit_fidelity - best["fidelity"]) < 1e-8

    payload = {
        "purpose": "recover weak fixed-magnetization-sector mixing in the N=9 spectral-flow state",
        "condition": {"N": n, "translations": graph["translations"],
                      "field_T": FIELD_T, "lambda": LAMBDA,
                      "target_level": 0},
        "ansatz": "two triangular Givens layers plus one triangular pair-mixing layer",
        "initial_computational_basis_index": initial_index,
        "parameter_count": len(operations),
        "two_qubit_parameterized_gate_count": sum(
            operation[0] in {"givens", "pair"} for operation in operations),
        "settings": {"restarts": RESTARTS, "maxiter": MAXITER,
                     "seed": SEED},
        "best": {
            "fidelity": best["fidelity"],
            "qiskit_fidelity": qiskit_fidelity,
            "success": best["success"],
            "message": best["message"],
            "iterations": best["iterations"],
            "function_gradient_evaluations":
                best["function_gradient_evaluations"],
            "parameters_radians": best["parameters"].tolist(),
        },
        "observable_check": {
            "exact_V_expectation_K": exact_v,
            "prepared_V_expectation_K": prepared_v,
            "V_expectation_error_K": prepared_v - exact_v,
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
            "This fits only the central-node ground state; all required excited states and nodes remain to be tested.",
            "Exact amplitudes define the optimization target, so the result tests ansatz expressibility rather than scalable eigenstate discovery.",
            "The topology-aware circuit excludes noise and the measurements needed to optimize parameters on hardware.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "sector_mixing_ansatz_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "best_fidelity": best["fidelity"],
        "V_expectation_error_K": prepared_v - exact_v,
        "all_to_all": payload["all_to_all"],
        "fake_sherbrooke": payload["fake_sherbrooke"],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
