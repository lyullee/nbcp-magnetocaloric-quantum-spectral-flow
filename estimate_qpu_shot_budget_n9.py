"""Exact thermal-variance shot budget for the N=9 M0-to-M3 correction."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import spin_hamiltonian_to_sparse_pauli
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
CENTRAL_T_K = 0.10
DERIVATIVE_HALF_WIDTH_K = 0.005
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
GL_WEIGHTS = GL_WEIGHTS / 2
GENERAL_ENTROPY_ERROR_KB_PER_SPIN = 0.001
STRESS_ENTROPY_ERROR_KB_PER_SPIN = 0.0002
OPTIMISTIC_CX_PER_N9_METTS_STATE = 10_800


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def group_variance_amplitudes(h, group_matrices, temperature):
    energies, vectors = np.linalg.eigh(h)
    shifted = energies - energies.min()
    probability = np.exp(-shifted / temperature)
    probability /= probability.sum()
    amplitudes = []
    variances = []
    for group in group_matrices:
        action = group @ vectors
        mean_by_eigenstate = np.real(np.sum(vectors.conj() * action, axis=0))
        second_by_eigenstate = np.sum(np.abs(action) ** 2, axis=0).real
        mean = float(np.dot(probability, mean_by_eigenstate))
        second = float(np.dot(probability, second_by_eigenstate))
        variance = max(0.0, second - mean**2)
        variances.append(variance)
        amplitudes.append(np.sqrt(variance))
    return np.asarray(amplitudes), variances


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    bonds0 = model_bonds(graph, "M0", boundary="tin_foil")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    h0_zero = dense_hamiltonian(n, bonds0, 0.0)
    h3_zero = dense_hamiltonian(n, bonds3, 0.0)
    v = h3_zero - h0_zero
    pauli0 = spin_hamiltonian_to_sparse_pauli(n, bonds0, 0.0)
    pauli3 = spin_hamiltonian_to_sparse_pauli(n, bonds3, 0.0)
    pauli_v = (pauli3 - pauli0).simplify(atol=1e-13)
    groups = pauli_v.group_commuting(qubit_wise=True)
    group_matrices = [np.asarray(group.to_matrix()) for group in groups]
    assert np.max(np.abs(sum(group_matrices) - v)) < 1e-11

    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        h0 = dense_hamiltonian(n, bonds0, field_k)
        temperature_amplitude_sums = {}
        detailed = []
        for temperature in (CENTRAL_T_K - DERIVATIVE_HALF_WIDTH_K,
                            CENTRAL_T_K + DERIVATIVE_HALF_WIDTH_K):
            weighted_amplitude_sum = 0.0
            for node, weight in zip(GL_NODES, GL_WEIGHTS, strict=True):
                h = h0 + float(node) * v
                amplitudes, variances = group_variance_amplitudes(
                    h, group_matrices, temperature)
                amplitude_sum = float(np.sum(amplitudes))
                weighted_amplitude_sum += float(weight) * amplitude_sum
                detailed.append({
                    "temperature_K": temperature,
                    "lambda": float(node),
                    "quadrature_weight": float(weight),
                    "sum_group_standard_deviations_K": amplitude_sum,
                    "group_variances_K2": variances,
                })
            temperature_amplitude_sums[temperature] = weighted_amplitude_sum

        # For S=-[F(T+dT)-F(T-dT)]/(2*dT), optimal allocation over
        # temperatures, quadrature nodes and commuting groups gives this bound.
        combined_amplitude = sum(temperature_amplitude_sums.values()) / (
            2 * DERIVATIVE_HALF_WIDTH_K)
        targets = {"general": GENERAL_ENTROPY_ERROR_KB_PER_SPIN}
        if abs(field_t - 1.45) < 1e-12:
            targets["stress_sign_change"] = STRESS_ENTROPY_ERROR_KB_PER_SPIN
        budgets = {}
        for label, entropy_error in targets.items():
            total_entropy_error = n * entropy_error
            shots = int(np.ceil((combined_amplitude / total_entropy_error) ** 2))
            quadratic_queries = int(np.ceil(np.sqrt(shots)))
            budgets[label] = {
                "target_entropy_error_kB_per_spin": entropy_error,
                "optimally_allocated_total_shots": shots,
                "optimistic_total_CX_using_weight2_QITE": shots * OPTIMISTIC_CX_PER_N9_METTS_STATE,
                "quadratic_precision_query_scale_sqrt_shots": quadratic_queries,
                "optimistic_CX_if_one_QITE_prep_per_query": quadratic_queries * OPTIMISTIC_CX_PER_N9_METTS_STATE,
            }
        row = {
            "field_T": field_t,
            "central_temperature_K": CENTRAL_T_K,
            "derivative_half_width_K": DERIVATIVE_HALF_WIDTH_K,
            "qubitwise_commuting_group_count": len(groups),
            "pauli_term_count": int(len(pauli_v)),
            "combined_derivative_amplitude": combined_amplitude,
            "budgets": budgets,
            "details": detailed,
        }
        rows.append(row)
        print(f"B={field_t:.2f}: general shots={budgets['general']['optimally_allocated_total_shots']:,}", flush=True)

    payload = {
        "purpose": "optimistic shot lower bound from exact thermal variances of QWC correction groups",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "observable": "V=M3-M0",
        "pauli_term_count": int(len(pauli_v)),
        "qubitwise_commuting_group_count": len(groups),
        "quadrature": {"type": "3-point Gauss-Legendre", "nodes": GL_NODES.tolist(),
                       "weights": GL_WEIGHTS.tolist()},
        "rows": rows,
        "scope_limits": [
            "Budgets assume exact Gibbs/QMETTS state preparation and optimal shot allocation.",
            "State-preparation error, QMETTS autocorrelation, calibration, readout and routing overhead are excluded.",
            "Group variances use the exact thermal mixed state; a finite METTS ensemble can require more shots.",
            "The CX product uses the optimistic all-to-all weight-2 QITE projection and is not a hardware runtime estimate.",
            "The square-root query scale only illustrates the precision scaling of ideal amplitude estimation; it excludes controlled state preparation, inverses, reflections and fault-tolerance overhead.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qpu_shot_budget_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": [
        {"field_T": row["field_T"], "budgets": row["budgets"]} for row in rows],
        "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
