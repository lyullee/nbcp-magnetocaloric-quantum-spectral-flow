"""Fit a physics-informed number-conserving circuit to one N=9 low state.

This is a state-preparation feasibility screen for the spectral-flow route.
It targets the dominant fixed-up-spin sector of the B=1.45 T, lambda=0.5
ground state using Givens mixers plus local and pair phases on triangular bonds.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import RZZGate, XXMinusYYGate, XXPlusYYGate
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ, directional_edges


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
LAMBDA = 0.5
LAYERS = (1, 2, 3)
RESTARTS = 3
MAXITER = 500
SEED = 20261012


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def operation_structure(n, edges, layers):
    operations = []
    for layer in range(layers):
        edge_order = edges if layer % 2 == 0 else list(reversed(edges))
        operations.extend(("rz", site) for site in range(n))
        operations.extend(("givens", i, j) for i, j in edge_order)
        operations.extend(("rzz", i, j) for i, j in edge_order)
    return operations


def apply_operation(state, operation, theta, inverse=False):
    theta = -theta if inverse else theta
    result = state.copy()
    kind = operation[0]
    indices = np.arange(state.size)
    if kind == "rz":
        site = operation[1]
        signed = 2 * ((indices >> site) & 1) - 1
        result *= np.exp(0.5j * theta * signed)
    elif kind == "rzz":
        i, j = operation[1:]
        zi = 1 - 2 * ((indices >> i) & 1)
        zj = 1 - 2 * ((indices >> j) & 1)
        result *= np.exp(-0.5j * theta * zi * zj)
    elif kind == "givens":
        i, j = operation[1:]
        mask = (((indices >> i) & 1) == 1) & (((indices >> j) & 1) == 0)
        a = indices[mask]
        b = a ^ (1 << i) ^ (1 << j)
        old_a = state[a]
        old_b = state[b]
        c, s = np.cos(theta), np.sin(theta)
        result[a] = c * old_a - s * old_b
        result[b] = s * old_a + c * old_b
    elif kind == "pair":
        i, j = operation[1:]
        mask = (((indices >> i) & 1) == 0) & (((indices >> j) & 1) == 0)
        a = indices[mask]
        b = a ^ (1 << i) ^ (1 << j)
        old_a = state[a]
        old_b = state[b]
        c, s = np.cos(theta), np.sin(theta)
        result[a] = c * old_a - s * old_b
        result[b] = s * old_a + c * old_b
    else:
        raise ValueError(kind)
    return result


def derivative_from_output(output, operation):
    derivative = np.zeros_like(output)
    kind = operation[0]
    indices = np.arange(output.size)
    if kind == "rz":
        site = operation[1]
        signed = 2 * ((indices >> site) & 1) - 1
        derivative = 0.5j * signed * output
    elif kind == "rzz":
        i, j = operation[1:]
        zi = 1 - 2 * ((indices >> i) & 1)
        zj = 1 - 2 * ((indices >> j) & 1)
        derivative = -0.5j * zi * zj * output
    elif kind == "givens":
        i, j = operation[1:]
        mask = (((indices >> i) & 1) == 1) & (((indices >> j) & 1) == 0)
        a = indices[mask]
        b = a ^ (1 << i) ^ (1 << j)
        derivative[a] = -output[b]
        derivative[b] = output[a]
    elif kind == "pair":
        i, j = operation[1:]
        mask = (((indices >> i) & 1) == 0) & (((indices >> j) & 1) == 0)
        a = indices[mask]
        b = a ^ (1 << i) ^ (1 << j)
        derivative[a] = -output[b]
        derivative[b] = output[a]
    else:
        raise ValueError(kind)
    return derivative


class FidelityObjective:
    def __init__(self, initial, target, operations):
        self.initial = initial
        self.target = target
        self.operations = operations
        self.evaluations = 0

    def value_gradient(self, parameters):
        self.evaluations += 1
        states = [self.initial]
        for operation, theta in zip(self.operations, parameters, strict=True):
            states.append(apply_operation(states[-1], operation, theta))
        overlap = np.vdot(self.target, states[-1])
        backward = self.target.copy()
        gradient = np.empty(len(parameters))
        for index in range(len(parameters) - 1, -1, -1):
            derivative = derivative_from_output(states[index + 1],
                                                self.operations[index])
            derivative_overlap = np.vdot(backward, derivative)
            gradient[index] = -2 * np.real(np.conj(overlap) *
                                           derivative_overlap)
            backward = apply_operation(backward, self.operations[index],
                                       parameters[index], inverse=True)
        return float(1 - abs(overlap)**2), gradient


def build_qiskit_circuit(n, initial_index, operations, parameters):
    circuit = QuantumCircuit(n)
    for site in range(n):
        if (initial_index >> site) & 1:
            circuit.x(site)
    for operation, theta in zip(operations, parameters, strict=True):
        if operation[0] == "rz":
            circuit.rz(float(theta), operation[1])
        elif operation[0] == "rzz":
            circuit.append(RZZGate(float(theta)), operation[1:])
        elif operation[0] == "givens":
            # beta=pi/2 implements the same real rotation; XXPlusYY uses
            # half of its theta argument as the subspace rotation angle.
            circuit.append(XXPlusYYGate(2 * float(theta), np.pi / 2),
                           operation[1:])
        elif operation[0] == "pair":
            circuit.append(XXMinusYYGate(2 * float(theta), np.pi / 2),
                           operation[1:])
    return circuit


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    edges = [(i, j) for i, j, _ in directional_edges(graph)]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(
        n, model_bonds(graph, "M0", boundary="tin_foil"), field_k)
    h3 = dense_hamiltonian(
        n, model_bonds(graph, "M3", boundary="tin_foil"), field_k)
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

    rng = np.random.default_rng(SEED)
    runs = []
    previous_best = None
    for layers in LAYERS:
        operations = operation_structure(n, edges, layers)
        candidates = []
        for restart in range(RESTARTS):
            if restart == 0 and previous_best is not None:
                x0 = np.zeros(len(operations))
                x0[:len(previous_best)] = previous_best
                x0[len(previous_best):] = rng.normal(
                    scale=0.03, size=len(operations) - len(previous_best))
            else:
                x0 = rng.normal(scale=0.15, size=len(operations))
            objective = FidelityObjective(initial, target, operations)
            result = minimize(
                objective.value_gradient, x0, jac=True, method="L-BFGS-B",
                bounds=[(-np.pi, np.pi)] * len(operations),
                options={"maxiter": MAXITER, "ftol": 1e-13,
                         "gtol": 1e-9, "maxls": 30})
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
            print(f"L={layers} restart={restart}: "
                  f"Fsector={sector_fidelity:.9f} Ffull={sector_weight*sector_fidelity:.9f} "
                  f"iter={result.nit}", flush=True)
        best = max(candidates, key=lambda row: row["sector_fidelity"])
        previous_best = best["parameters"].copy()
        circuit = build_qiskit_circuit(n, initial_index, operations,
                                       best["parameters"])
        simulated = np.asarray(Statevector.from_instruction(circuit).data)
        qiskit_sector_fidelity = float(abs(np.vdot(target, simulated))**2)
        all_to_all = transpile(
            circuit, basis_gates=["rz", "sx", "x", "cx"],
            optimization_level=2, seed_transpiler=SEED)
        hardware = transpile(
            circuit, backend=FakeSherbrooke(), optimization_level=2,
            seed_transpiler=SEED)
        runs.append({
            "layers": layers,
            "parameter_count": len(operations),
            "best_sector_fidelity": best["sector_fidelity"],
            "best_full_state_fidelity": best["full_state_fidelity"],
            "qiskit_sector_fidelity": qiskit_sector_fidelity,
            "optimizer": {k: v for k, v in best.items() if k != "parameters"},
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
        })
        # Validate that the custom simulator and the emitted Qiskit circuit
        # implement the same ansatz convention.
        if abs(qiskit_sector_fidelity - best["sector_fidelity"]) > 1e-8:
            raise RuntimeError("custom/Qiskit ansatz convention mismatch")

    payload = {
        "purpose": "physics-informed N=9 low-energy state-preparation screen",
        "condition": {"N": n, "translations": graph["translations"],
                      "field_T": FIELD_T, "lambda": LAMBDA,
                      "target_level": 0},
        "ansatz": "triangular-bond Givens mixers plus RZ and RZZ phases",
        "dominant_up_spin_sector": dominant_q,
        "dominant_sector_weight": sector_weight,
        "initial_computational_basis_index": initial_index,
        "edge_count": len(edges),
        "settings": {"layers": list(LAYERS), "restarts": RESTARTS,
                     "maxiter": MAXITER, "seed": SEED},
        "runs": runs,
        "scope_limits": [
            "Only the ground state at the central lambda node is fitted; spectral flow needs several excited states and all three nodes.",
            "Fidelity is optimized in the dominant fixed-up-spin sector, with full-state fidelity reduced by omitted M3 sector mixing.",
            "Exact target amplitudes are used for local fitting and do not constitute a scalable state-discovery method.",
            "FakeSherbrooke compilation is topology-aware but excludes noise and optimizer measurements.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "number_conserving_ansatz_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "dominant_sector_weight": sector_weight,
        "runs": runs,
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
