"""Optimistic local-statevector screen of one global variational-QITE step.

This uses exact state amplitudes to fit the QITE generator, so it is a lower
bound on algorithmic error and does not yet include the measurements needed to
form the QITE linear system on hardware.
"""
from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm
import qiskit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.synthesis import LieTrotter

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds,
                                  spin_hamiltonian_to_sparse_pauli)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
STEP_SIZES = (0.02, 0.05, 0.10)  # inverse-kelvin imaginary-time steps
CORRELATION_TIMES = (0.0, 0.5, 2.0)
POOL_WEIGHTS = (1, 2, 3)


def pauli_pool(n: int, maximum_weight: int):
    labels = []
    for weight in range(1, maximum_weight + 1):
        for sites in combinations(range(n), weight):
            for axes in product("XYZ", repeat=weight):
                label = ["I"] * n
                for site, axis in zip(sites, axes, strict=True):
                    label[n - 1 - site] = axis
                labels.append("".join(label))
    operator = SparsePauliOp(labels, np.ones(len(labels)))
    matrices = np.asarray([SparsePauliOp(label).to_matrix() for label in labels])
    return labels, operator, matrices


def normalized_imaginary_state(h, initial, tau):
    state = expm(-tau * h) @ initial
    return state / np.linalg.norm(state)


def exact_step(h, state, dt):
    target = expm(-dt * h) @ state
    return target / np.linalg.norm(target)


def fit_qite_generator(h, state, labels, matrices):
    energy = float(np.real(np.vdot(state, h @ state)))
    tangent = -(h @ state - energy * state)
    columns = -1j * np.einsum("pij,j->pi", matrices, state).T
    design = np.vstack((columns.real, columns.imag))
    target = np.concatenate((tangent.real, tangent.imag))
    coefficients, _, rank, singular = np.linalg.lstsq(design, target, rcond=1e-10)
    residual = design @ coefficients - target
    keep = np.abs(coefficients) > 1e-10
    if np.any(keep):
        generator = SparsePauliOp(
            [label for label, use in zip(labels, keep, strict=True) if use],
            coefficients[keep],
        ).simplify(atol=1e-10)
    else:
        generator = SparsePauliOp("I" * int(round(np.log2(len(state)))), coeffs=[0.0])
    nonzero = singular[singular > singular[0] * 1e-10] if len(singular) else singular
    condition = float(nonzero[0] / nonzero[-1]) if len(nonzero) else float("inf")
    return generator, {
        "energy_K": energy,
        "tangent_norm": float(np.linalg.norm(target)),
        "relative_tangent_residual": float(np.linalg.norm(residual) / np.linalg.norm(target)),
        "linear_system_rows": int(design.shape[0]),
        "linear_system_columns": int(design.shape[1]),
        "effective_rank": int(rank),
        "retained_pauli_terms": int(len(generator)),
        "resolved_condition_number": condition,
    }


def circuit_step(generator, state, dt):
    n = generator.num_qubits
    gate = PauliEvolutionGate(generator, time=dt, synthesis=LieTrotter(reps=1))
    circuit = QuantumCircuit(n)
    circuit.append(gate, range(n))
    evolved = Statevector(state).evolve(circuit).data
    compiled = transpile(circuit, basis_gates=["rz", "sx", "x", "cx"],
                         optimization_level=3, seed_transpiler=20261011)
    operations = {str(name): int(count) for name, count in compiled.count_ops().items()}
    return evolved, {
        "logical_depth": int(circuit.decompose(reps=3).depth()),
        "compiled_all_to_all_depth": int(compiled.depth()),
        "compiled_all_to_all_two_qubit_gates": int(operations.get("cx", 0)),
        "compiled_operation_counts": operations,
    }


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges(3, 2)
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds = open_weak_model_bonds(graph, "M2")
    exact_operator = SpinHamiltonianOperator(n, bonds, field_k)
    h = np.asarray(exact_operator @ np.eye(1 << n, dtype=complex))
    pauli_h = spin_hamiltonian_to_sparse_pauli(n, bonds, field_k)
    assert np.max(np.abs(pauli_h.to_matrix() - h)) < 1e-12

    initial = np.zeros(1 << n, dtype=complex)
    initial[0b101010] = 1.0  # non-eigenstate product configuration
    pools = {weight: pauli_pool(n, weight) for weight in POOL_WEIGHTS}
    rows = []
    for tau in CORRELATION_TIMES:
        state = normalized_imaginary_state(h, initial, tau)
        for weight in POOL_WEIGHTS:
            labels, _, matrices = pools[weight]
            generator, fit = fit_qite_generator(h, state, labels, matrices)
            for dt in STEP_SIZES:
                target = exact_step(h, state, dt)
                evolved, resources = circuit_step(generator, state, dt)
                fidelity = float(abs(np.vdot(target, evolved)) ** 2)
                row = {
                    "correlation_time_inverse_K": tau,
                    "imaginary_step_inverse_K": dt,
                    "maximum_pauli_weight": weight,
                    "pool_size": len(labels),
                    **fit,
                    **resources,
                    "one_step_infidelity": 1 - fidelity,
                    "projected_steps_to_90mK_beta_over_2": int(np.ceil((1 / (2 * 0.09)) / dt)),
                    "projected_CX_per_90mK_METTS_state": int(
                        np.ceil((1 / (2 * 0.09)) / dt) *
                        resources["compiled_all_to_all_two_qubit_gates"]),
                }
                rows.append(row)
                print(f"tau={tau:.1f}, weight<={weight}, dt={dt:.2f}: "
                      f"residual={fit['relative_tangent_residual']:.3e}, "
                      f"infidelity={1-fidelity:.3e}", flush=True)

    payload = {
        "purpose": "optimistic N=6 one-step QITE feasibility screen before sampled QMETTS",
        "runtime": {"qiskit": qiskit.__version__, "statevector": True},
        "cluster": {"N": n, "boundary": "open 3x2 triangular patch", "model": "M2"},
        "field_T": FIELD_T,
        "initial_state": "computational product state |101010>",
        "rows": rows,
        "scope_limits": [
            "Generator coefficients use exact state amplitudes rather than measured QITE matrices.",
            "The Pauli exponential uses one first-order Lie-Trotter step and optimistic all-to-all connectivity.",
            "One-step agreement does not establish a full beta/2 thermal preparation or QMETTS sampling cost.",
            "This N=6 M2 patch is a circuit-development test; N=9 periodic M3 is the scientific target.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qite_tangent_screen_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
