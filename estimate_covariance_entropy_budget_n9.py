"""Optimistic direct-covariance entropy estimator for the N=9 correction."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURE_K = 0.10
GL_NODES, GL_WEIGHTS = np.polynomial.legendre.leggauss(3)
GL_NODES = (GL_NODES + 1) / 2
GL_WEIGHTS = GL_WEIGHTS / 2
GENERAL_ERROR = 0.001
STRESS_ERROR = 0.0002


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def moments(operator, vectors, probability):
    action = operator @ vectors
    diagonal = np.real(np.sum(vectors.conj() * action, axis=0))
    second = np.sum(np.abs(action) ** 2, axis=0).real
    mean = float(np.dot(probability, diagonal))
    variance = max(0.0, float(np.dot(probability, second)) - mean**2)
    return mean, variance


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    h0_zero = dense_hamiltonian(n, model_bonds(graph, "M0"), 0.0)
    h3_zero = dense_hamiltonian(n, model_bonds(graph, "M3", boundary="tin_foil"), 0.0)
    v = h3_zero - h0_zero
    exact_targets = json.loads((ROOT / "results" / "quantum_local_exact_target_N9.json").read_text(encoding="utf-8"))
    exact_lookup = {(row["field_T"], row["temperature_K"]): row["exact_delta_S_kB_per_spin"]
                    for row in exact_targets["rows"] if row["quadrature_order"] == 5}

    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        h0 = dense_hamiltonian(n, model_bonds(graph, "M0"), field_k)
        covariance_integral = 0.0
        ideal_amplitude_integral = 0.0
        separate_amplitude_integral = 0.0
        details = []
        for node, weight in zip(GL_NODES, GL_WEIGHTS, strict=True):
            h = h0 + float(node) * v
            energies, vectors = np.linalg.eigh(h)
            probability = np.exp(-(energies - energies.min()) / TEMPERATURE_K)
            probability /= probability.sum()
            mean_h, var_h = moments(h, vectors, probability)
            mean_v, var_v = moments(v, vectors, probability)
            anticommutator = (h @ v + v @ h) / 2
            mean_a, var_a = moments(anticommutator, vectors, probability)
            covariance = mean_a - mean_h * mean_v
            # Linearized covariance fluctuation. Exact mean-dependent constants
            # do not affect its variance.
            fluctuation_operator = anticommutator - mean_v * h - mean_h * v
            _, var_fluctuation = moments(fluctuation_operator, vectors, probability)
            covariance_integral += float(weight) * covariance
            ideal_amplitude_integral += float(weight) * np.sqrt(var_fluctuation) / TEMPERATURE_K**2
            separate_amplitude_integral += float(weight) * (
                np.sqrt(var_a) + abs(mean_v) * np.sqrt(var_h) +
                abs(mean_h) * np.sqrt(var_v)) / TEMPERATURE_K**2
            details.append({
                "lambda": float(node),
                "weight": float(weight),
                "mean_H_K": mean_h,
                "mean_V_K": mean_v,
                "symmetrized_covariance_K2": covariance,
                "variance_ideal_linearized_covariance_K4": var_fluctuation,
                "variance_A_K4": var_a,
                "variance_H_K2": var_h,
                "variance_V_K2": var_v,
            })
        estimated_delta_s = -covariance_integral / (TEMPERATURE_K**2 * n)
        exact_delta_s = exact_lookup[(field_t, TEMPERATURE_K)]
        targets = {"general": GENERAL_ERROR}
        if abs(field_t - 1.45) < 1e-12:
            targets["stress_sign_change"] = STRESS_ERROR
        budgets = {}
        for label, error in targets.items():
            total_error = n * error
            ideal_shots = int(np.ceil((ideal_amplitude_integral / total_error) ** 2))
            separate_shots = int(np.ceil((separate_amplitude_integral / total_error) ** 2))
            budgets[label] = {
                "target_entropy_error_kB_per_spin": error,
                "ideal_mean_dependent_covariance_observable_shots": ideal_shots,
                "separate_A_H_V_delta_method_shots": separate_shots,
            }
        row = {
            "field_T": field_t,
            "temperature_K": TEMPERATURE_K,
            "covariance_GL3_delta_S_kB_per_spin": estimated_delta_s,
            "exact_delta_S_kB_per_spin": exact_delta_s,
            "quadrature_error_kB_per_spin": estimated_delta_s - exact_delta_s,
            "ideal_amplitude_integral": ideal_amplitude_integral,
            "separate_amplitude_integral": separate_amplitude_integral,
            "budgets": budgets,
            "details": details,
        }
        rows.append(row)
        print(f"B={field_t:.2f}: ideal covariance shots={budgets['general']['ideal_mean_dependent_covariance_observable_shots']:,}", flush=True)

    payload = {
        "purpose": "optimistic direct-covariance alternative to finite-difference entropy correction",
        "identity": "Delta S = - integral_0^1 Cov_sym(H_lambda,V)/T^2 d lambda",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "quadrature": {"type": "3-point Gauss-Legendre", "nodes": GL_NODES.tolist(),
                       "weights": GL_WEIGHTS.tolist()},
        "rows": rows,
        "scope_limits": [
            "The ideal covariance budget assumes direct projective measurement of a mean-dependent many-body observable.",
            "The separate A,H,V budget uses first-order independent-error propagation and optimal shot allocation.",
            "Neither budget includes Pauli decomposition, basis-change circuits, QITE, QMETTS correlation or hardware noise.",
            "These are optimistic screening bounds, not executable QPU resource estimates.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "covariance_entropy_budget_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": [
        {"field_T": row["field_T"], "delta_S": row["covariance_GL3_delta_S_kB_per_spin"],
         "quadrature_error": row["quadrature_error_kB_per_spin"], "budgets": row["budgets"]}
        for row in rows], "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
