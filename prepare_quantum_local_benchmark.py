"""Freeze the exact N=9 target for the local quantum-algorithm screening.

The target is the coupling-constant correction from M0 to the full, tin-foil
M3 finite-cluster Hamiltonian.  Energies are expressed as H/k_B in kelvin.
This script produces classical exact answers only; it does not emulate or run
a quantum circuit.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import time

import numpy as np

from material_models import model_operators
from screen_model import lattice


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURES_K = (0.09, 0.10, 0.12)
QUADRATURE_ORDERS = (3, 5)
N = 9


def dense_blocks(model: str, field_t: float) -> list[np.ndarray]:
    graph = lattice((3, 0), (0, 3))
    blocks = []
    for operator in model_operators(graph, model, field_t, boundary="tin_foil"):
        eye = np.eye(operator.shape[0], dtype=complex)
        matrix = np.asarray(operator @ eye)
        hermiticity_error = float(np.max(np.abs(matrix - matrix.conj().T)))
        if hermiticity_error > 1e-11:
            raise RuntimeError(f"non-Hermitian block: {hermiticity_error}")
        blocks.append((matrix + matrix.conj().T) / 2)
    return blocks


def spectrum_with_v(h0: list[np.ndarray], v: list[np.ndarray], lam: float):
    energies = []
    v_diagonal = []
    for base, correction in zip(h0, v, strict=True):
        eigvals, eigvecs = np.linalg.eigh(base + lam * correction)
        diag_v = np.real(np.einsum("ij,ij->j", eigvecs.conj(), correction @ eigvecs))
        energies.append(eigvals)
        v_diagonal.append(diag_v)
    return np.concatenate(energies), np.concatenate(v_diagonal)


def thermo(spectrum, temperature_k: float):
    energies, v_diagonal = spectrum
    minimum = float(np.min(energies))
    shifted = energies - minimum
    weights = np.exp(-shifted / temperature_k)
    partition_shifted = float(np.sum(weights))
    probability = weights / partition_shifted
    mean_shifted = float(np.dot(probability, shifted))
    return {
        "free_energy_K_per_spin": (minimum - temperature_k * np.log(partition_shifted)) / N,
        "entropy_kB_per_spin": (np.log(partition_shifted) + mean_shifted / temperature_k) / N,
        "mean_V_K_per_spin": float(np.dot(probability, v_diagonal)) / N,
    }


def gauss_legendre_nodes(order: int):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1) / 2, w / 2


def main():
    started = time.perf_counter()
    rows = []
    field_metadata = []
    for field_t in FIELDS_T:
        h0 = dense_blocks("M0", field_t)
        h1 = dense_blocks("M3", field_t)
        v = [target - base for base, target in zip(h0, h1, strict=True)]
        endpoint = {
            0.0: spectrum_with_v(h0, v, 0.0),
            1.0: spectrum_with_v(h0, v, 1.0),
        }
        quadrature = {}
        for order in QUADRATURE_ORDERS:
            nodes, weights = gauss_legendre_nodes(order)
            quadrature[order] = {
                "nodes": nodes,
                "weights": weights,
                "spectra": [spectrum_with_v(h0, v, float(node)) for node in nodes],
            }

        correction_norms = [float(np.linalg.norm(block, 2)) for block in v]
        field_metadata.append({
            "field_T": field_t,
            "parity_block_dimensions": [int(block.shape[0]) for block in h0],
            "V_spectral_norm_K_by_block": correction_norms,
        })

        for temperature_k in TEMPERATURES_K:
            base = thermo(endpoint[0.0], temperature_k)
            target = thermo(endpoint[1.0], temperature_k)
            exact_delta_f = target["free_energy_K_per_spin"] - base["free_energy_K_per_spin"]
            exact_delta_s = target["entropy_kB_per_spin"] - base["entropy_kB_per_spin"]
            for order, data in quadrature.items():
                expectation_values = [
                    thermo(spectrum, temperature_k)["mean_V_K_per_spin"]
                    for spectrum in data["spectra"]
                ]
                integrated_delta_f = float(np.dot(data["weights"], expectation_values))
                rows.append({
                    "N": N,
                    "boundary": "3x3 periodic, M3 tin-foil Ewald",
                    "field_T": field_t,
                    "temperature_K": temperature_k,
                    "quadrature_order": order,
                    "exact_delta_F_K_per_spin": exact_delta_f,
                    "integrated_delta_F_K_per_spin": integrated_delta_f,
                    "delta_F_error_K_per_spin": integrated_delta_f - exact_delta_f,
                    "exact_delta_S_kB_per_spin": exact_delta_s,
                    "M0_entropy_kB_per_spin": base["entropy_kB_per_spin"],
                    "M3_entropy_kB_per_spin": target["entropy_kB_per_spin"],
                    "lambda_nodes": [float(x) for x in data["nodes"]],
                    "mean_V_K_per_spin_at_nodes": expectation_values,
                })
        print(f"completed B={field_t:.2f} T", flush=True)

    payload = {
        "purpose": "exact target for local quantum-algorithm screening; no quantum circuit or QPU",
        "hamiltonian_path": "M0 + lambda*(M3-M0)",
        "units": "Hamiltonian divided by k_B, in kelvin",
        "fields_T": list(FIELDS_T),
        "temperatures_K": list(TEMPERATURES_K),
        "quadrature_orders": list(QUADRATURE_ORDERS),
        "field_metadata": field_metadata,
        "rows": rows,
        "acceptance_criteria": {
            "delta_S_absolute_error_kB_per_spin_max": 0.001,
            "sign_must_match": True,
            "propagated_field_uncertainty_T_max": 0.02,
            "qpu_submission": "forbidden until ideal and noisy local simulations pass",
        },
        "scope_limits": [
            "N=9 is a finite periodic cluster and is not a bulk prediction.",
            "The target validates an observable route, not quantum advantage.",
            "The experimental critical-field discrepancy is not attributed to M3 alone.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    OUT.mkdir(exist_ok=True)
    json_path = OUT / "quantum_local_exact_target_N9.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    csv_path = OUT / "quantum_local_exact_target_N9.csv"
    scalar_keys = [
        "N", "boundary", "field_T", "temperature_K", "quadrature_order",
        "exact_delta_F_K_per_spin", "integrated_delta_F_K_per_spin",
        "delta_F_error_K_per_spin", "exact_delta_S_kB_per_spin",
        "M0_entropy_kB_per_spin", "M3_entropy_kB_per_spin",
    ]
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=scalar_keys)
        writer.writeheader()
        writer.writerows([{key: row[key] for key in scalar_keys} for row in rows])

    max_errors = {
        order: max(abs(row["delta_F_error_K_per_spin"]) for row in rows
                   if row["quadrature_order"] == order)
        for order in QUADRATURE_ORDERS
    }
    print(json.dumps({
        "json": str(json_path),
        "csv": str(csv_path),
        "max_delta_F_quadrature_error": max_errors,
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
