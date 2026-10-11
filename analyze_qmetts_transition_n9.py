"""Exact X-basis QMETTS variance and mixing audit for periodic N=9 M3."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import hadamard

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURES_K = (0.09, 0.10, 0.12)
TARGET_V_ERROR_K_PER_SPIN = 1e-4


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    dimension = 1 << n
    x_basis = hadamard(dimension).astype(complex) / np.sqrt(dimension)
    bonds0 = model_bonds(graph, "M0", boundary="tin_foil")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    # V is field independent.
    v = dense_hamiltonian(n, bonds3, 0.0) - dense_hamiltonian(n, bonds0, 0.0)
    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        h3 = dense_hamiltonian(n, bonds3, field_k)
        energies, vectors = np.linalg.eigh(h3)
        gap = float(energies[1] - energies[0])
        for temperature in TEMPERATURES_K:
            factors = np.exp(-(energies - energies[0]) / (2 * temperature))
            propagator = (vectors * factors) @ vectors.conj().T
            thermal_columns = propagator @ x_basis
            weights = np.sum(np.abs(thermal_columns) ** 2, axis=0)
            states = thermal_columns / np.sqrt(weights)[np.newaxis, :]
            observable_values = np.real(np.sum(states.conj() * (v @ states), axis=0))
            probability = weights / weights.sum()
            mean_v = float(np.dot(probability, observable_values))
            variance_v = float(np.dot(probability, (observable_values - mean_v) ** 2))
            amplitudes = x_basis.conj().T @ states
            transition = np.abs(amplitudes) ** 2
            transition /= transition.sum(axis=0, keepdims=True)
            stationary_error = float(np.max(np.abs(transition @ probability - probability)))
            eigenvalues = np.linalg.eigvals(transition)
            moduli = np.sort(np.abs(eigenvalues))[::-1]
            unit_multiplicity = int(np.sum(np.abs(eigenvalues - 1) < 1e-9))
            independent_samples = int(np.ceil(
                variance_v / (n**2 * TARGET_V_ERROR_K_PER_SPIN**2)))
            row = {
                "field_T": field_t,
                "temperature_K": temperature,
                "spectral_gap_H_K": gap,
                "exact_mean_V_K": mean_v,
                "X_METTS_variance_V_K2": variance_v,
                "X_METTS_standard_deviation_V_K": float(np.sqrt(variance_v)),
                "independent_samples_for_1e-4_K_per_spin_SE": max(1, independent_samples),
                "transition_unit_eigenvalue_multiplicity": unit_multiplicity,
                "transition_second_eigenvalue_modulus": float(moduli[1]),
                "transition_spectral_gap": float(1 - moduli[1]),
                "stationary_distribution_max_error": stationary_error,
                "minimum_basis_probability": float(np.min(probability)),
                "maximum_basis_probability": float(np.max(probability)),
            }
            rows.append(row)
            print(f"B={field_t:.2f}, T={temperature:.2f}: gap={row['transition_spectral_gap']:.3f}, "
                  f"n_ind={row['independent_samples_for_1e-4_K_per_spin_SE']}", flush=True)

    assert all(row["transition_unit_eigenvalue_multiplicity"] == 1 for row in rows)
    assert max(row["stationary_distribution_max_error"] for row in rows) < 1e-10
    payload = {
        "purpose": "exact N=9 X-basis QMETTS mixing and observable-variance audit",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "observable": "V=M3-M0",
        "target_independent_sample_standard_error_K_per_spin": TARGET_V_ERROR_K_PER_SPIN,
        "rows": rows,
        "scope_limits": [
            "Imaginary-time propagation and transition matrices are exact; QITE and hardware errors are absent.",
            "Independent-sample counts do not include observable shot noise.",
            "N=9 is a finite-cluster algorithm benchmark, not a bulk-material prediction.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qmetts_transition_audit_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
