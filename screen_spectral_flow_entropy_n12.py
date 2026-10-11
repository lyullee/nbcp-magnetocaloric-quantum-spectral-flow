"""N=12 size-extension of the low-energy spectral-flow entropy route."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import eigsh

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import spin_hamiltonian_to_sparse_pauli
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURES_K = (0.09, 0.10, 0.12)
K_PER_PARITY = 48
GENERAL_DELTA_S_ERROR = 0.001
STRESS_DELTA_S_ERROR = 0.0002
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
GL_WEIGHTS = GL_WEIGHTS / 2


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


def binary_entropy(probability):
    if probability <= 0 or probability >= 1:
        return 0.0
    return float(-probability * np.log(probability) -
                 (1 - probability) * np.log(1 - probability))


def omitted_entropy_bound(energies, temperature, dimension):
    """Fannes-Audenaert bound from a worst-case omitted partition tail."""
    weights = np.exp(-(energies - energies[0]) / temperature)
    partial_z = float(np.sum(weights))
    omitted = dimension - len(energies)
    tail_z_upper = omitted * float(weights[-1])
    tail_mass_upper = tail_z_upper / (partial_z + tail_z_upper)
    entropy_bound = (binary_entropy(tail_mass_upper) +
                     tail_mass_upper * np.log(dimension - 1))
    return tail_mass_upper, entropy_bound


def lowest_global(h0_ops, h3_ops, lam, count_per_parity, seed, vectors):
    records = []
    for parity in (0, 1):
        operator = h0_ops[parity] + float(lam) * (
            h3_ops[parity] - h0_ops[parity])
        rng = np.random.default_rng(seed + parity)
        v0 = rng.normal(size=operator.shape[0]) + 1j * rng.normal(
            size=operator.shape[0])
        if vectors:
            values, local_vectors = eigsh(
                operator, k=count_per_parity, which="SA", tol=1e-10,
                v0=v0)
            for index, value in enumerate(values):
                full = np.zeros(1 << h0_ops[parity].n, dtype=complex)
                full[h0_ops[parity].states] = local_vectors[:, index]
                records.append((float(value), full))
        else:
            values = eigsh(operator, k=count_per_parity, which="SA",
                           tol=1e-10, v0=v0, return_eigenvectors=False)
            records.extend((float(value), None) for value in values)
    records.sort(key=lambda item: item[0])
    energies = np.asarray([item[0] for item in records])
    if not vectors:
        return energies, None
    return energies, np.column_stack([item[1] for item in records])


def group_standard_deviations(group_matrices, vectors):
    rows = []
    for matrix in group_matrices:
        action = matrix @ vectors
        mean = np.real(np.sum(vectors.conj() * action, axis=0))
        second = np.sum(np.abs(action) ** 2, axis=0).real
        rows.append(np.sqrt(np.maximum(0.0, second - mean**2)))
    return np.asarray(rows)


def main():
    started = time.perf_counter()
    graph = lattice((2, 2), (-2, 4))
    n = graph["N"]
    dimension = 1 << n
    bonds0 = model_bonds(graph, "M0", boundary="tin_foil")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    pauli_v = (spin_hamiltonian_to_sparse_pauli(n, bonds3, 0.0) -
               spin_hamiltonian_to_sparse_pauli(n, bonds0, 0.0)).simplify(
                   atol=1e-13)
    groups = pauli_v.group_commuting(qubit_wise=True)
    group_matrices = [group.to_matrix(sparse=True) for group in groups]
    rows = []

    for field_index, field_t in enumerate(FIELDS_T):
        field_k = GZ * MU_B_OVER_KB * field_t
        h0_ops = {parity: SpinHamiltonianOperator(
            n, bonds0, field_k, parity=parity) for parity in (0, 1)}
        h3_ops = {parity: SpinHamiltonianOperator(
            n, bonds3, field_k, parity=parity) for parity in (0, 1)}
        e0, _ = lowest_global(h0_ops, h3_ops, 0.0, K_PER_PARITY,
                              20261012 + 100 * field_index, False)
        e3, _ = lowest_global(h0_ops, h3_ops, 1.0, K_PER_PARITY,
                              20261112 + 100 * field_index, False)
        print(f"B={field_t:.2f}: endpoint spectra complete", flush=True)

        tail = {}
        chosen_rank = {}
        reference_entropy = {}
        for temperature in TEMPERATURES_K:
            target = (STRESS_DELTA_S_ERROR
                      if abs(field_t - 1.45) < 1e-12 and
                      abs(temperature - 0.09) < 1e-12
                      else GENERAL_DELTA_S_ERROR)
            tail[temperature] = {}
            reference_entropy[temperature] = {}
            for model, energies in (("M0", e0), ("M3", e3)):
                mass, bound = omitted_entropy_bound(
                    energies, temperature, dimension)
                tail[temperature][model] = {
                    "omitted_probability_upper_bound": mass,
                    "entropy_kB_error_upper_bound": bound,
                }
                reference_entropy[temperature][model] = entropy(
                    energies, temperature)
            tail_delta_s_bound = sum(
                tail[temperature][model]["entropy_kB_error_upper_bound"]
                for model in ("M0", "M3")) / n
            selected = None
            for rank in range(1, len(e0) + 1):
                truncation_bound = sum(abs(
                    entropy(energies[:rank], temperature) -
                    reference_entropy[temperature][model])
                    for model, energies in (("M0", e0), ("M3", e3))) / n
                if truncation_bound + tail_delta_s_bound <= target:
                    selected = rank
                    break
            if selected is None:
                raise RuntimeError(
                    f"K={len(e0)} cannot certify B={field_t}, T={temperature}")
            chosen_rank[temperature] = selected

        max_rank = max(chosen_rank.values())
        shifts = np.zeros(max_rank)
        node_data = []
        for node_index, (node, weight) in enumerate(zip(
                GL_NODES, GL_WEIGHTS, strict=True)):
            energies, vectors = lowest_global(
                h0_ops, h3_ops, float(node), K_PER_PARITY,
                20262012 + 1000 * field_index + 10 * node_index, True)
            vectors = vectors[:, :max_rank]
            v_action = sum(matrix @ vectors for matrix in group_matrices)
            v_mean = np.real(np.sum(vectors.conj() * v_action, axis=0))
            shifts += float(weight) * v_mean
            group_stddev = group_standard_deviations(group_matrices, vectors)
            node_data.append({
                "lambda": float(node),
                "weight": float(weight),
                "V_expectation_K_by_level": v_mean.tolist(),
                "V_group_stddev_K_by_group_and_level": group_stddev.tolist(),
            })
            print(f"B={field_t:.2f}: lambda={node:.6f} complete", flush=True)
        reconstructed = e0[:max_rank] + shifts

        for temperature in TEMPERATURES_K:
            rank = chosen_rank[temperature]
            target = (STRESS_DELTA_S_ERROR
                      if abs(field_t - 1.45) < 1e-12 and
                      abs(temperature - 0.09) < 1e-12
                      else GENERAL_DELTA_S_ERROR)
            exact_reference_delta_s = (
                reference_entropy[temperature]["M3"] -
                reference_entropy[temperature]["M0"]) / n
            exact_truncated_delta_s = (
                entropy(e3[:rank], temperature) -
                entropy(e0[:rank], temperature)) / n
            reconstructed_delta_s = (
                entropy(reconstructed[:rank], temperature) -
                entropy(e0[:rank], temperature)) / n
            gradient = np.abs(entropy_gradient(
                reconstructed[:rank], temperature)) / n
            amplitude = 0.0
            for item in node_data:
                stddev = np.asarray(
                    item["V_group_stddev_K_by_group_and_level"])
                amplitude += item["weight"] * float(np.sum(
                    stddev[:, :rank] * gradient[None, :]))
            shots = int(np.ceil((amplitude / target) ** 2))
            tail_delta_s_bound = sum(
                tail[temperature][model]["entropy_kB_error_upper_bound"]
                for model in ("M0", "M3")) / n
            rows.append({
                "field_T": field_t,
                "temperature_K": temperature,
                "low_energy_rank": rank,
                "target_delta_S_error_kB_per_spin": target,
                "reference_delta_S_kB_per_spin": exact_reference_delta_s,
                "reference_tail_delta_S_error_bound_kB_per_spin":
                    tail_delta_s_bound,
                "exact_truncated_delta_S_kB_per_spin":
                    exact_truncated_delta_s,
                "reconstructed_delta_S_kB_per_spin":
                    reconstructed_delta_s,
                "low_rank_plus_quadrature_error_vs_reference_kB_per_spin":
                    reconstructed_delta_s - exact_reference_delta_s,
                "maximum_selected_level_energy_error_K": float(np.max(
                    np.abs(reconstructed[:rank] - e3[:rank]))),
                "sign_matches_reference": bool(
                    np.sign(reconstructed_delta_s) ==
                    np.sign(exact_reference_delta_s)),
                "optimistic_optimal_V_measurement_shots": shots,
                "weighted_measurement_amplitude": amplitude,
                "tail_bounds": tail[temperature],
            })
            print(f"B={field_t:.2f} T={temperature:.2f}: K={rank}, "
                  f"err={rows[-1]['low_rank_plus_quadrature_error_vs_reference_kB_per_spin']:.2e}, "
                  f"shots={shots:,}", flush=True)

    payload = {
        "purpose": "N=12 size-extension of low-energy spectral-flow entropy reconstruction",
        "cluster": {"N": n, "translations": graph["translations"],
                    "boundary": "periodic; M3 tin-foil Ewald"},
        "settings": {
            "low_eigenpairs_per_parity": K_PER_PARITY,
            "quadrature": "3-point Gauss-Legendre",
            "eigsh_tolerance": 1e-10,
        },
        "observable": {
            "pauli_term_count": int(len(pauli_v)),
            "qwc_group_count": len(groups),
            "coefficient_l1_K": float(np.sum(np.abs(pauli_v.coeffs))),
        },
        "rows": rows,
        "scope_limits": [
            "The endpoint reference uses 96 low levels and a worst-case omitted-tail/Fannes-Audenaert error bound, not full N=12 diagonalization.",
            "Exact eigsh eigenvectors stand in for a future quantum eigensolver, so shot counts omit the dominant state-preparation/optimization risk.",
            "The classical M0 low-energy spectrum is treated as an available baseline.",
            "Routing, noise, mitigation, overlap penalties and optimizer repetitions are excluded.",
            "N=12 is a finite cluster and does not establish thermodynamic-limit quantum resource scaling.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "spectral_flow_entropy_N12.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "elapsed_seconds": payload["elapsed_seconds"],
        "summary": [{
            "field_T": row["field_T"],
            "temperature_K": row["temperature_K"],
            "rank": row["low_energy_rank"],
            "error": row["low_rank_plus_quadrature_error_vs_reference_kB_per_spin"],
            "tail_bound": row["reference_tail_delta_S_error_bound_kB_per_spin"],
            "shots": row["optimistic_optimal_V_measurement_shots"],
        } for row in rows],
    }, indent=2))


if __name__ == "__main__":
    main()
