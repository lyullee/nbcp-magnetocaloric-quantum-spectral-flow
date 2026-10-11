"""Exact transition-matrix audit of N=6 QMETTS collapse bases."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm, hadamard

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import open_triangular_directional_edges, open_weak_model_bonds
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
TEMPERATURES_K = (0.09, 0.10, 0.12)


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def ensemble(propagator, observable, basis):
    thermal_columns = propagator @ basis
    weights = np.sum(np.abs(thermal_columns) ** 2, axis=0)
    observable_columns = observable @ thermal_columns
    values = np.real(np.sum(thermal_columns.conj() * observable_columns, axis=0)) / weights
    probabilities = weights / np.sum(weights)
    mean = float(np.dot(probabilities, values))
    variance = float(np.dot(probabilities, (values - mean) ** 2))
    return thermal_columns, weights, probabilities, values, mean, variance


def transition(thermal_columns, weights, target_basis):
    amplitudes = target_basis.conj().T @ thermal_columns
    matrix = np.abs(amplitudes) ** 2 / weights[np.newaxis, :]
    assert np.max(np.abs(np.sum(matrix, axis=0) - 1)) < 1e-11
    return matrix


def spectral_summary(matrix):
    eigenvalues = np.linalg.eigvals(matrix)
    moduli = np.sort(np.abs(eigenvalues))[::-1]
    multiplicity = int(np.sum(np.abs(eigenvalues - 1) < 1e-9))
    second = float(moduli[1])
    return {
        "unit_eigenvalue_multiplicity": multiplicity,
        "second_eigenvalue_modulus": second,
        "spectral_gap": float(1 - second),
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
    z_basis = np.eye(dimension, dtype=complex)
    x_basis = hadamard(dimension).astype(complex) / np.sqrt(dimension)
    odd = np.array([index.bit_count() % 2 == 1 for index in range(dimension)])

    rows = []
    for temperature in TEMPERATURES_K:
        propagator = expm(-h2 / (2 * temperature))
        z = ensemble(propagator, v, z_basis)
        x = ensemble(propagator, v, x_basis)
        z_columns, z_weights, pz, _, z_mean, z_variance = z
        x_columns, x_weights, px, _, x_mean, x_variance = x
        t_z_from_z = transition(z_columns, z_weights, z_basis)
        t_x_from_z = transition(z_columns, z_weights, x_basis)
        t_z_from_x = transition(x_columns, x_weights, z_basis)
        t_x_from_x = transition(x_columns, x_weights, x_basis)
        two_step_z = t_z_from_x @ t_x_from_z
        stationary_z_error = float(np.max(np.abs(two_step_z @ pz - pz)))
        stationary_x_error = float(np.max(np.abs(t_x_from_z @ pz - px)))
        parity_crossing = float(np.max(t_z_from_z[np.ix_(odd, ~odd)]))
        trace_gibbs = float(np.trace(propagator.conj().T @ propagator).real)
        exact_mean = float(np.trace(v @ (propagator @ propagator)).real / trace_gibbs)
        row = {
            "temperature_K": temperature,
            "exact_thermal_mean_V_K": exact_mean,
            "Z_ensemble_mean_error_K": z_mean - exact_mean,
            "X_ensemble_mean_error_K": x_mean - exact_mean,
            "Z_ensemble_variance_K2": z_variance,
            "X_ensemble_variance_K2": x_variance,
            "Z_stationary_odd_parity_weight": float(np.sum(pz[odd])),
            "Z_only_max_cross_parity_transition": parity_crossing,
            "Z_only_transition": spectral_summary(t_z_from_z),
            "X_only_transition": spectral_summary(t_x_from_x),
            "alternating_ZX_two_step_transition": spectral_summary(two_step_z),
            "alternating_ZX_Z_stationary_max_error": stationary_z_error,
            "alternating_ZX_X_stationary_max_error": stationary_x_error,
        }
        rows.append(row)
        print(f"T={temperature:.2f}: Z-only mult={row['Z_only_transition']['unit_eigenvalue_multiplicity']}, "
              f"ZX gap={row['alternating_ZX_two_step_transition']['spectral_gap']:.3e}", flush=True)

    max_ensemble_error = max(
        max(abs(row["Z_ensemble_mean_error_K"]), abs(row["X_ensemble_mean_error_K"]))
        for row in rows)
    assert max_ensemble_error < 1e-10
    assert all(row["Z_only_transition"]["unit_eigenvalue_multiplicity"] >= 2 for row in rows)
    assert all(row["alternating_ZX_two_step_transition"]["unit_eigenvalue_multiplicity"] == 1
               for row in rows)

    payload = {
        "purpose": "exact QMETTS transition audit for parity trapping and alternating-basis ergodicity",
        "cluster": {"N": n, "boundary": "open 3x2 triangular patch", "model": "M2"},
        "observable": "V=M2-M0",
        "field_T": FIELD_T,
        "rows": rows,
        "status": "PASS_ALTERNATING_ZX; STOP_Z_ONLY",
        "interpretation": [
            "Each complete Z or X METTS ensemble reproduces the exact Gibbs observable.",
            "Z-only collapse has at least two stationary sectors because it cannot cross spin parity.",
            "Alternating Z/X collapse has a unique stationary distribution in this N=6 test.",
        ],
        "scope_limits": [
            "This is an exact small-system transition audit with no QITE, shot or hardware noise.",
            "A unique N=6 stationary distribution does not guarantee fast mixing at N=9 or N=12.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qmetts_transition_audit_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "status": payload["status"],
                      "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
