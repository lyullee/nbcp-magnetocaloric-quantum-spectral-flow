"""Optimistic full beta/2 QITE trajectory for the N=6 circuit test.

The QITE coefficients are obtained from exact statevectors at every step.  The
result therefore tests accumulated ansatz/Trotter error while deliberately
excluding hardware measurement and linear-system noise.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import Statevector
from qiskit.synthesis import LieTrotter

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds,
                                  spin_hamiltonian_to_sparse_pauli)
from qite_tangent_screen_n6 import circuit_step, fit_qite_generator, pauli_pool
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
TARGET_TEMPERATURE_K = 0.10
TARGET_TAU = 1 / (2 * TARGET_TEMPERATURE_K)
STEP_SIZES = (0.02, 0.05, 0.10)
POOL_WEIGHT = 2


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def apply_generator(generator, state, dt):
    circuit = QuantumCircuit(generator.num_qubits)
    circuit.append(PauliEvolutionGate(generator, time=dt,
                                      synthesis=LieTrotter(reps=1)),
                   range(generator.num_qubits))
    return Statevector(state).evolve(circuit).data


def normalized_exact(h, initial, tau):
    state = expm(-tau * h) @ initial
    return state / np.linalg.norm(state)


def hadamard_transform(state):
    output = np.asarray(state, dtype=complex).copy()
    stride = 1
    scale = 1 / np.sqrt(2)
    while stride < len(output):
        for start in range(0, len(output), 2 * stride):
            a = output[start:start + stride].copy()
            b = output[start + stride:start + 2 * stride].copy()
            output[start:start + stride] = (a + b) * scale
            output[start + stride:start + 2 * stride] = (a - b) * scale
        stride *= 2
    return output


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges(3, 2)
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds0 = open_weak_model_bonds(graph, "M0")
    bonds2 = open_weak_model_bonds(graph, "M2")
    h0 = dense_hamiltonian(n, bonds0, field_k)
    h2 = dense_hamiltonian(n, bonds2, field_k)
    v = h2 - h0
    pauli_h2 = spin_hamiltonian_to_sparse_pauli(n, bonds2, field_k)
    assert np.max(np.abs(pauli_h2.to_matrix() - h2)) < 1e-12

    z_initial = np.zeros(1 << n, dtype=complex)
    z_initial[0b101010] = 1.0
    initial_states = {
        "Z:101010": z_initial,
        "X:101010": hadamard_transform(z_initial),
    }
    labels, _, matrices = pauli_pool(n, POOL_WEIGHT)

    rows = []
    for initial_label, initial in initial_states.items():
        exact = normalized_exact(h2, initial, TARGET_TAU)
        exact_energy = float(np.real(np.vdot(exact, h2 @ exact)))
        exact_v = float(np.real(np.vdot(exact, v @ exact)))
        for nominal_dt in STEP_SIZES:
            steps = int(round(TARGET_TAU / nominal_dt))
            dt = TARGET_TAU / steps
            state = initial.copy()
            residuals = []
            retained = []
            conditions = []
            resource_samples = []
            checkpoints = {0, steps // 4, steps // 2, 3 * steps // 4, steps - 1}
            trajectory_start = time.perf_counter()
            for step in range(steps):
                generator, fit = fit_qite_generator(h2, state, labels, matrices)
                residuals.append(fit["relative_tangent_residual"])
                retained.append(fit["retained_pauli_terms"])
                conditions.append(fit["resolved_condition_number"])
                if step in checkpoints:
                    _, resources = circuit_step(generator, state, dt)
                    resource_samples.append({"step": step, **resources})
                state = apply_generator(generator, state, dt)
                state /= np.linalg.norm(state)
                if (step + 1) % max(1, steps // 5) == 0:
                    print(f"{initial_label}, dt={nominal_dt:.2f}: {step+1}/{steps}", flush=True)

            fidelity = float(abs(np.vdot(exact, state)) ** 2)
            energy = float(np.real(np.vdot(state, h2 @ state)))
            mean_v = float(np.real(np.vdot(state, v @ state)))
            cx_samples = [sample["compiled_all_to_all_two_qubit_gates"]
                          for sample in resource_samples]
            depth_samples = [sample["compiled_all_to_all_depth"]
                             for sample in resource_samples]
            rows.append({
                "initial_product_basis_state": initial_label,
                "nominal_step_inverse_K": nominal_dt,
                "actual_step_inverse_K": dt,
                "step_count": steps,
                "final_infidelity": 1 - fidelity,
                "exact_energy_K": exact_energy,
                "qite_energy_K": energy,
                "energy_error_K": energy - exact_energy,
                "exact_mean_V_K": exact_v,
                "qite_mean_V_K": mean_v,
                "mean_V_error_K": mean_v - exact_v,
                "max_relative_tangent_residual": float(max(residuals)),
                "median_relative_tangent_residual": float(np.median(residuals)),
                "median_retained_pauli_terms": float(np.median(retained)),
                "max_resolved_condition_number": float(max(conditions)),
                "resource_samples": resource_samples,
                "sampled_CX_range_per_step": [int(min(cx_samples)), int(max(cx_samples))],
                "sampled_depth_range_per_step": [int(min(depth_samples)), int(max(depth_samples))],
                "projected_CX_range_per_METTS_state": [int(steps * min(cx_samples)),
                                                        int(steps * max(cx_samples))],
                "trajectory_wall_seconds": time.perf_counter() - trajectory_start,
            })

    payload = {
        "purpose": "optimistic accumulated-QITE screen for one N=6 METTS state",
        "cluster": {"N": n, "boundary": "open 3x2 triangular patch", "model": "M2"},
        "field_T": FIELD_T,
        "temperature_K": TARGET_TEMPERATURE_K,
        "target_beta_over_2_inverse_K": TARGET_TAU,
        "initial_states": list(initial_states),
        "maximum_pauli_weight": POOL_WEIGHT,
        "pool_size": len(labels),
        "rows": rows,
        "scope_limits": [
            "Every QITE generator is fitted from the exact statevector, with no shot or hardware noise.",
            "Only one product-state trajectory is tested; QMETTS burn-in and autocorrelation are absent.",
            "CX counts assume all-to-all connectivity and sample only five steps per trajectory.",
            "The N=6 M2 patch is a circuit-development test, not the N=9 M3 scientific benchmark.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qite_full_trajectory_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": rows,
                      "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
