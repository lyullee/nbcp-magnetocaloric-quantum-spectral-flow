"""End-to-end ideal-statevector X-basis QITE-QMETTS test for N=6.

All 64 X-basis product states are propagated with global weight-2 QITE.  Their
collapse probabilities form an approximate QMETTS transition matrix whose
stationary thermal observable is compared with exact diagonalization.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm, hadamard

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import open_triangular_directional_edges, open_weak_model_bonds
from qite_full_trajectory_n6 import apply_generator
from qite_tangent_screen_n6 import fit_qite_generator, pauli_pool
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
TEMPERATURE_K = 0.10
TAU = 1 / (2 * TEMPERATURE_K)
STEP = 0.10
STEPS = int(round(TAU / STEP))
POOL_WEIGHT = 2


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def stationary_distribution(transition):
    eigenvalues, eigenvectors = np.linalg.eig(transition)
    order = np.argsort(np.abs(eigenvalues - 1))
    stationary = np.real_if_close(eigenvectors[:, order[0]], tol=1000).real
    if stationary.sum() < 0:
        stationary *= -1
    stationary[np.abs(stationary) < 1e-14] = 0
    if np.min(stationary) < -1e-10:
        raise RuntimeError("stationary vector has material negative entries")
    stationary = np.clip(stationary, 0, None)
    stationary /= stationary.sum()
    moduli = np.sort(np.abs(eigenvalues))[::-1]
    return stationary, {
        "unit_eigenvalue_error": float(abs(eigenvalues[order[0]] - 1)),
        "second_eigenvalue_modulus": float(moduli[1]),
        "spectral_gap": float(1 - moduli[1]),
        "stationary_residual": float(np.max(np.abs(transition @ stationary - stationary))),
    }


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges(3, 2)
    n = graph["N"]
    dimension = 1 << n
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    v = h2 - h0
    propagator = expm(-TAU * h2)
    x_basis = hadamard(dimension).astype(complex) / np.sqrt(dimension)
    labels, _, matrices = pauli_pool(n, POOL_WEIGHT)

    exact_weights = np.zeros(dimension)
    exact_values = np.zeros(dimension)
    qite_values = np.zeros(dimension)
    fidelities = np.zeros(dimension)
    transition = np.zeros((dimension, dimension))
    max_tangent_residuals = np.zeros(dimension)
    retained_medians = np.zeros(dimension)
    for index in range(dimension):
        initial = x_basis[:, index]
        exact_unnormalized = propagator @ initial
        exact_weights[index] = float(np.vdot(exact_unnormalized, exact_unnormalized).real)
        exact_state = exact_unnormalized / np.sqrt(exact_weights[index])
        exact_values[index] = float(np.vdot(exact_state, v @ exact_state).real)

        state = initial.copy()
        residuals = []
        retained = []
        for _ in range(STEPS):
            generator, fit = fit_qite_generator(h2, state, labels, matrices)
            residuals.append(fit["relative_tangent_residual"])
            retained.append(fit["retained_pauli_terms"])
            state = apply_generator(generator, state, STEP)
            state /= np.linalg.norm(state)
        qite_values[index] = float(np.vdot(state, v @ state).real)
        fidelities[index] = float(abs(np.vdot(exact_state, state)) ** 2)
        amplitudes = x_basis.conj().T @ state
        transition[:, index] = np.abs(amplitudes) ** 2
        transition[:, index] /= transition[:, index].sum()
        max_tangent_residuals[index] = max(residuals)
        retained_medians[index] = np.median(retained)
        if (index + 1) % 8 == 0:
            print(f"completed X-basis states {index+1}/{dimension}", flush=True)

    exact_probability = exact_weights / exact_weights.sum()
    exact_mean = float(np.dot(exact_probability, exact_values))
    best_case_qite_mean = float(np.dot(exact_probability, qite_values))
    qite_stationary, spectral = stationary_distribution(transition)
    end_to_end_mean = float(np.dot(qite_stationary, qite_values))
    probability_l1 = float(np.sum(np.abs(qite_stationary - exact_probability)))
    per_state_error = qite_values - exact_values
    payload = {
        "purpose": "all-state ideal QITE-QMETTS X-basis benchmark",
        "cluster": {"N": n, "boundary": "open 3x2 triangular patch", "model": "M2"},
        "observable": "V=M2-M0",
        "field_T": FIELD_T,
        "temperature_K": TEMPERATURE_K,
        "beta_over_2_inverse_K": TAU,
        "imaginary_step_inverse_K": STEP,
        "step_count_per_METTS_state": STEPS,
        "maximum_pauli_weight": POOL_WEIGHT,
        "pool_size": len(labels),
        "exact_thermal_mean_V_K": exact_mean,
        "best_case_exact_weight_QITE_mean_V_K": best_case_qite_mean,
        "best_case_bias_K": best_case_qite_mean - exact_mean,
        "best_case_bias_K_per_spin": (best_case_qite_mean - exact_mean) / n,
        "QITE_transition_stationary_mean_V_K": end_to_end_mean,
        "end_to_end_bias_K": end_to_end_mean - exact_mean,
        "end_to_end_bias_K_per_spin": (end_to_end_mean - exact_mean) / n,
        "stationary_probability_L1_error": probability_l1,
        "transition_spectral_summary": spectral,
        "statewise": {
            "minimum_fidelity": float(np.min(fidelities)),
            "median_infidelity": float(np.median(1 - fidelities)),
            "maximum_abs_V_error_K": float(np.max(np.abs(per_state_error))),
            "exact_probability_weighted_RMSE_V_K": float(np.sqrt(np.dot(exact_probability, per_state_error**2))),
            "maximum_tangent_residual": float(np.max(max_tangent_residuals)),
            "median_retained_pauli_terms": float(np.median(retained_medians)),
        },
        "acceptance": {
            "provisional_V_bias_limit_K_per_spin": 1e-4,
            "passes_provisional_V_bias": abs((end_to_end_mean - exact_mean) / n) <= 1e-4,
        },
        "scope_limits": [
            "QITE coefficients use exact statevectors; measurement and linear-system noise are absent.",
            "Observable shots and QMETTS finite-chain variance are absent because the stationary transition matrix is solved exactly.",
            "This tests one temperature, field and lambda endpoint on N=6 M2; it does not yet establish delta-F or delta-S accuracy.",
            "The scientific target remains N=9 periodic M3.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qite_qmetts_xbasis_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items()
                      if key not in {"scope_limits"}}, indent=2))


if __name__ == "__main__":
    main()
