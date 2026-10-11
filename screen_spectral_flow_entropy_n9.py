"""Screen low-energy spectral-flow reconstruction of the M0-to-M3 entropy.

For ordered levels E_k(lambda), Hellmann--Feynman gives
  E_k(1)-E_k(0) = integral <k,lambda|V|k,lambda> d lambda.
Only the small correction V=M3-M0 is measured.  A classical M0 spectrum and
the reconstructed M3 low-energy spectrum then yield the entropy correction
without temperature finite differences or direct entropy measurements.
"""
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
TEMPERATURES_K = (0.09, 0.10, 0.12, 0.30)
QUADRATURE_ORDERS = (3, 5)
GENERAL_DELTA_S_ERROR = 0.001
STRESS_DELTA_S_ERROR = 0.0002


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def entropy(energies, temperature):
    probability = np.exp(-(energies - energies[0]) / temperature)
    probability /= probability.sum()
    return float(-np.dot(probability[probability > 0],
                         np.log(probability[probability > 0])))


def entropy_gradient(energies, temperature):
    probability = np.exp(-(energies - energies[0]) / temperature)
    probability /= probability.sum()
    mean = float(np.dot(probability, energies))
    return -(probability * (energies - mean)) / temperature**2


def group_standard_deviations(groups, vectors):
    result = []
    for group in groups:
        matrix = np.asarray(group.to_matrix())
        action = matrix @ vectors
        mean = np.real(np.sum(vectors.conj() * action, axis=0))
        second = np.sum(np.abs(action) ** 2, axis=0).real
        result.append(np.sqrt(np.maximum(0.0, second - mean**2)))
    return np.asarray(result)


def main():
    started = time.perf_counter()
    low_rank = json.loads((ROOT / "results" /
                           "low_energy_quantum_route_N9.json").read_text(
                               encoding="utf-8"))
    rank_lookup = {
        (row["field_T"], row["temperature_K"]):
            row["minimum_common_low_energy_rank"]
        for row in low_rank["delta_S_comparisons"]
    }

    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    bonds0 = model_bonds(graph, "M0", boundary="tin_foil")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    h0_zero = dense_hamiltonian(n, bonds0, 0.0)
    h3_zero = dense_hamiltonian(n, bonds3, 0.0)
    v = h3_zero - h0_zero
    pauli_v = (spin_hamiltonian_to_sparse_pauli(n, bonds3, 0.0) -
               spin_hamiltonian_to_sparse_pauli(n, bonds0, 0.0)).simplify(
                   atol=1e-13)
    groups = pauli_v.group_commuting(qubit_wise=True)
    assert np.max(np.abs(np.asarray(pauli_v.to_matrix()) - v)) < 1e-11

    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        h0 = dense_hamiltonian(n, bonds0, field_k)
        h3 = dense_hamiltonian(n, bonds3, field_k)
        exact_e0 = np.linalg.eigvalsh(h0)
        exact_e3 = np.linalg.eigvalsh(h3)
        max_rank = max(rank_lookup[(field_t, temperature)]
                       for temperature in TEMPERATURES_K)

        for order in QUADRATURE_ORDERS:
            nodes, weights = np.polynomial.legendre.leggauss(order)
            nodes = (nodes + 1) / 2
            weights = weights / 2
            shifts = np.zeros(max_rank)
            node_data = []
            for node, weight in zip(nodes, weights, strict=True):
                energies, vectors = np.linalg.eigh(h0 + float(node) * v)
                low_vectors = vectors[:, :max_rank]
                v_action = v @ low_vectors
                v_mean = np.real(np.sum(low_vectors.conj() * v_action, axis=0))
                shifts += float(weight) * v_mean
                group_stddev = group_standard_deviations(groups, low_vectors)
                node_data.append({
                    "lambda": float(node),
                    "weight": float(weight),
                    "V_expectation_K_by_level": v_mean.tolist(),
                    "V_group_stddev_K_by_group_and_level": group_stddev.tolist(),
                })
            reconstructed = exact_e0[:max_rank] + shifts

            for temperature in TEMPERATURES_K:
                rank = rank_lookup[(field_t, temperature)]
                target = (STRESS_DELTA_S_ERROR
                          if abs(field_t - 1.45) < 1e-12 and
                          abs(temperature - 0.09) < 1e-12
                          else GENERAL_DELTA_S_ERROR)
                exact_delta_s = (entropy(exact_e3, temperature) -
                                 entropy(exact_e0, temperature)) / n
                reconstructed_delta_s = (
                    entropy(reconstructed[:rank], temperature) -
                    entropy(exact_e0[:rank], temperature)) / n
                # Separate quadrature error from the already-certified low-rank
                # error by comparing exact and reconstructed truncated spectra.
                exact_truncated_delta_s = (
                    entropy(exact_e3[:rank], temperature) -
                    entropy(exact_e0[:rank], temperature)) / n
                quadrature_error = reconstructed_delta_s - exact_truncated_delta_s
                total_error = reconstructed_delta_s - exact_delta_s

                gradient = np.abs(entropy_gradient(
                    reconstructed[:rank], temperature)) / n
                amplitude = 0.0
                for item in node_data:
                    group_stddev = np.asarray(
                        item["V_group_stddev_K_by_group_and_level"])
                    amplitude += item["weight"] * float(
                        np.sum(group_stddev[:, :rank] * gradient[None, :]))
                shots = int(np.ceil((amplitude / target) ** 2))

                rows.append({
                    "field_T": field_t,
                    "temperature_K": temperature,
                    "quadrature_order": order,
                    "low_energy_rank": rank,
                    "target_delta_S_error_kB_per_spin": target,
                    "exact_full_delta_S_kB_per_spin": exact_delta_s,
                    "exact_truncated_delta_S_kB_per_spin": exact_truncated_delta_s,
                    "reconstructed_delta_S_kB_per_spin": reconstructed_delta_s,
                    "spectral_flow_quadrature_error_kB_per_spin": quadrature_error,
                    "combined_low_rank_and_quadrature_error_kB_per_spin": total_error,
                    "maximum_low_level_energy_error_K": float(np.max(np.abs(
                        reconstructed[:rank] - exact_e3[:rank]))),
                    "optimistic_optimal_V_measurement_shots": shots,
                    "weighted_measurement_amplitude": amplitude,
                    "node_data": node_data,
                })
                print(f"B={field_t:.2f} T={temperature:.2f} GL{order}: "
                      f"K={rank}, err={total_error:.2e}, shots={shots:,}",
                      flush=True)

    payload = {
        "purpose": "low-energy spectral-flow entropy correction using only V=M3-M0 measurements",
        "identity": "E_k(1)-E_k(0)=integral_0^1 <k,lambda|V|k,lambda> d lambda",
        "cluster": {"N": n, "translations": graph["translations"],
                    "boundary": "3x3 periodic; M3 tin-foil Ewald"},
        "observable": {
            "pauli_term_count": int(len(pauli_v)),
            "qwc_group_count": len(groups),
            "coefficient_l1_K": float(np.sum(np.abs(pauli_v.coeffs))),
        },
        "rows": rows,
        "acceptance": {
            "general_delta_S_abs_error_kB_per_spin": GENERAL_DELTA_S_ERROR,
            "stress_B1_45_T0_09_abs_error_kB_per_spin": STRESS_DELTA_S_ERROR,
            "stress_sign_must_match": True,
        },
        "scope_limits": [
            "M0 low-energy levels are treated as a classical baseline; the budget covers only the M3-M0 spectral correction.",
            "Exact eigenvectors stand in for a future VQD/QSE/eigensolver and therefore make the shot counts optimistic.",
            "Ordered-level Hellmann-Feynman integration is insensitive to isolated label exchanges in thermodynamic sums, but finite-size avoided crossings can increase state-tracking cost.",
            "Shot budgets omit eigenstate optimization, overlap measurements, routing, noise, error mitigation and autocorrelation.",
            "The same low-rank and quadrature accuracy must be retested at N=12 before any QPU run.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "spectral_flow_entropy_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "elapsed_seconds": payload["elapsed_seconds"],
        "GL3": [{
            "field_T": row["field_T"],
            "temperature_K": row["temperature_K"],
            "rank": row["low_energy_rank"],
            "error": row["combined_low_rank_and_quadrature_error_kB_per_spin"],
            "shots": row["optimistic_optimal_V_measurement_shots"],
        } for row in rows if row["quadrature_order"] == 3],
    }, indent=2))


if __name__ == "__main__":
    main()
