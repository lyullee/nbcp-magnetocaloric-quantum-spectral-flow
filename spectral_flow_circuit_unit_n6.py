"""End-to-end N=6 circuit unit test for low-energy spectral flow.

The open 3x2 patch and M0-to-M2 path are a circuit test only.  The scientific
material benchmark remains the periodic N=9/N=12 M0-to-M3 calculation.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import qiskit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds,
                                  spin_hamiltonian_to_sparse_pauli)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.05
TEMPERATURE_K = 0.10
TARGET_DELTA_S_ERROR = 0.001
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
GL_WEIGHTS = GL_WEIGHTS / 2


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def entropy(energies):
    probability = np.exp(-(energies - energies[0]) / TEMPERATURE_K)
    probability /= probability.sum()
    return float(-np.dot(probability[probability > 0],
                         np.log(probability[probability > 0])))


def entropy_gradient(energies):
    probability = np.exp(-(energies - energies[0]) / TEMPERATURE_K)
    probability /= probability.sum()
    mean = float(np.dot(probability, energies))
    return -(probability * (energies - mean)) / TEMPERATURE_K**2


def group_standard_deviations(group_matrices, vectors):
    rows = []
    for matrix in group_matrices:
        action = matrix @ vectors
        mean = np.real(np.sum(vectors.conj() * action, axis=0))
        second = np.sum(np.abs(action) ** 2, axis=0).real
        rows.append(np.sqrt(np.maximum(0.0, second - mean**2)))
    return np.asarray(rows)


def state_preparation_circuit(state):
    state = np.asarray(state, complex).copy()
    state[np.abs(state) < 1e-12] = 0
    state /= np.linalg.norm(state)
    circuit = QuantumCircuit(int(np.log2(state.size)))
    circuit.append(StatePreparation(state, normalize=True), circuit.qubits)
    return circuit.decompose(reps=15), state


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges()
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds0 = open_weak_model_bonds(graph, "M0")
    bonds2 = open_weak_model_bonds(graph, "M2")
    h0 = dense_hamiltonian(n, bonds0, field_k)
    h2 = dense_hamiltonian(n, bonds2, field_k)
    v = h2 - h0
    e0 = np.linalg.eigvalsh(h0)
    e2 = np.linalg.eigvalsh(h2)
    exact_delta_s = (entropy(e2) - entropy(e0)) / n
    rank = next(rank for rank in range(1, (1 << n) + 1)
                if (abs(entropy(e0[:rank]) - entropy(e0)) +
                    abs(entropy(e2[:rank]) - entropy(e2))) / n
                <= TARGET_DELTA_S_ERROR)

    pauli0 = spin_hamiltonian_to_sparse_pauli(n, bonds0, 0.0)
    pauli2 = spin_hamiltonian_to_sparse_pauli(n, bonds2, 0.0)
    pauli_v = (pauli2 - pauli0).simplify(atol=1e-13)
    groups = pauli_v.group_commuting(qubit_wise=True)
    group_matrices = [np.asarray(group.to_matrix()) for group in groups]
    assert np.max(np.abs(np.asarray(pauli_v.to_matrix()) - v)) < 1e-11

    shifts = np.zeros(rank)
    node_rows = []
    circuit_rows = []
    backend = FakeSherbrooke()
    for node_index, (node, weight) in enumerate(zip(
            GL_NODES, GL_WEIGHTS, strict=True)):
        energies, vectors = np.linalg.eigh(h0 + float(node) * v)
        vectors = vectors[:, :rank]
        v_mean = np.real(np.sum(vectors.conj() * (v @ vectors), axis=0))
        shifts += float(weight) * v_mean
        group_stddev = group_standard_deviations(group_matrices, vectors)
        node_rows.append({
            "lambda": float(node),
            "weight": float(weight),
            "V_expectation_K_by_level": v_mean.tolist(),
            "V_group_stddev_K_by_group_and_level": group_stddev.tolist(),
        })
        for level in range(rank):
            circuit, target = state_preparation_circuit(vectors[:, level])
            prepared = np.asarray(Statevector.from_instruction(circuit).data)
            fidelity = float(abs(np.vdot(target, prepared))**2)
            circuit_v = float(np.real(
                Statevector(prepared).expectation_value(pauli_v)))
            exact_v = float(v_mean[level])
            all_to_all = transpile(
                circuit, basis_gates=["rz", "sx", "x", "cx"],
                optimization_level=1)
            hardware = transpile(circuit, backend=backend,
                                 optimization_level=1, seed_transpiler=20261012)
            circuit_rows.append({
                "lambda": float(node),
                "level": level,
                "statevector_fidelity": fidelity,
                "circuit_V_expectation_K": circuit_v,
                "exact_V_expectation_K": exact_v,
                "V_expectation_abs_error_K": abs(circuit_v - exact_v),
                "all_to_all": {
                    "cx_count": int(all_to_all.count_ops().get("cx", 0)),
                    "depth": int(all_to_all.depth()),
                    "operation_count": int(all_to_all.size()),
                },
                "fake_sherbrooke": {
                    "ecr_count": int(hardware.count_ops().get("ecr", 0)),
                    "cx_count": int(hardware.count_ops().get("cx", 0)),
                    "depth": int(hardware.depth()),
                    "operation_count": int(hardware.size()),
                    "physical_qubits": [int(q) for q in hardware.layout.initial_layout.get_physical_bits()],
                },
            })
            print(f"lambda={node:.6f} level={level}: "
                  f"CX={circuit_rows[-1]['all_to_all']['cx_count']} "
                  f"ECR={circuit_rows[-1]['fake_sherbrooke']['ecr_count']}",
                  flush=True)

    reconstructed = e0[:rank] + shifts
    reconstructed_delta_s = (entropy(reconstructed) - entropy(e0[:rank])) / n
    exact_truncated_delta_s = (entropy(e2[:rank]) - entropy(e0[:rank])) / n
    gradient = np.abs(entropy_gradient(reconstructed)) / n
    amplitude = 0.0
    for row in node_rows:
        stddev = np.asarray(row["V_group_stddev_K_by_group_and_level"])
        amplitude += row["weight"] * float(np.sum(
            stddev[:, :rank] * gradient[None, :]))
    shots = int(np.ceil((amplitude / TARGET_DELTA_S_ERROR)**2))

    assert min(row["statevector_fidelity"] for row in circuit_rows) > 1 - 1e-12
    assert max(row["V_expectation_abs_error_K"] for row in circuit_rows) < 1e-10
    assert abs(reconstructed_delta_s - exact_delta_s) <= TARGET_DELTA_S_ERROR

    payload = {
        "purpose": "N=6 executable circuit unit test for spectral-flow entropy reconstruction",
        "scientific_scope": "circuit pipeline only; open M0-to-M2 patch is not the periodic material benchmark",
        "runtime": {"qiskit": qiskit.__version__,
                    "backend_model": backend.name},
        "condition": {"N": n, "geometry": "3x2 open triangular patch",
                      "field_T": FIELD_T, "temperature_K": TEMPERATURE_K,
                      "path": "M0 + lambda*(M2-M0)"},
        "observable": {"pauli_term_count": int(len(pauli_v)),
                       "qwc_group_count": len(groups),
                       "coefficient_l1_K": float(np.sum(np.abs(pauli_v.coeffs)))},
        "low_energy_rank": rank,
        "exact_delta_S_kB_per_spin": exact_delta_s,
        "exact_truncated_delta_S_kB_per_spin": exact_truncated_delta_s,
        "reconstructed_delta_S_kB_per_spin": reconstructed_delta_s,
        "combined_error_kB_per_spin": reconstructed_delta_s - exact_delta_s,
        "target_error_kB_per_spin": TARGET_DELTA_S_ERROR,
        "optimistic_optimal_V_measurement_shots": shots,
        "node_rows": node_rows,
        "circuits": circuit_rows,
        "status": "PASS",
        "scope_limits": [
            "StatePreparation is synthesized from exact classical eigenvectors and is an upper-bound circuit unit test, not a scalable eigensolver.",
            "FakeSherbrooke transpilation checks topology and native-gate expansion but does not simulate device noise.",
            "Shot budget excludes eigenstate discovery, optimizer repetitions and error mitigation.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "spectral_flow_circuit_unit_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "status": payload["status"],
        "rank": rank, "error": payload["combined_error_kB_per_spin"],
        "shots": shots,
        "minimum_fidelity": min(row["statevector_fidelity"] for row in circuit_rows),
        "maximum_all_to_all_CX": max(row["all_to_all"]["cx_count"] for row in circuit_rows),
        "maximum_fake_sherbrooke_ECR": max(row["fake_sherbrooke"]["ecr_count"] for row in circuit_rows),
        "maximum_fake_sherbrooke_depth": max(row["fake_sherbrooke"]["depth"] for row in circuit_rows),
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
