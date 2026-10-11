"""Exact N=6 QMETTS sampling and autocorrelation benchmark.

This isolates Markov-chain sampling from QITE state-preparation error.  It
compares a Z-only collapse basis, which preserves parity, with alternating Z/X
product bases.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import open_triangular_directional_edges, open_weak_model_bonds
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
TEMPERATURES_K = (0.09, 0.10, 0.12)
BURN_IN = 500
SAMPLES = 20_000
SEED = 20261011


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def hadamard_transform(state):
    """Unitary H^tensor(n), using Qiskit's little-endian computational order."""
    output = np.asarray(state, dtype=complex).copy()
    size = len(output)
    stride = 1
    scale = 1 / np.sqrt(2)
    while stride < size:
        for start in range(0, size, 2 * stride):
            a = output[start:start + stride].copy()
            b = output[start + stride:start + 2 * stride].copy()
            output[start:start + stride] = (a + b) * scale
            output[start + stride:start + 2 * stride] = (a - b) * scale
        stride *= 2
    return output


def product_state(index, basis, dimension):
    state = np.zeros(dimension, dtype=complex)
    state[index] = 1.0
    return state if basis == "Z" else hadamard_transform(state)


def integrated_autocorrelation(values):
    values = np.asarray(values, dtype=float)
    centered = values - np.mean(values)
    n = len(values)
    fft_size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, fft_size)
    covariance = np.fft.irfft(spectrum * spectrum.conj(), fft_size)[:n]
    covariance /= np.arange(n, 0, -1)
    if covariance[0] <= 0:
        return 0.5, 1.0
    acf = covariance / covariance[0]
    stop = 1
    while stop < n and acf[stop] > 0:
        stop += 1
    tau = max(0.5, 0.5 + float(np.sum(acf[1:stop])))
    effective = min(float(n), n / (2 * tau))
    return tau, effective


def exact_thermal_mean(h, observable, temperature):
    energies, vectors = np.linalg.eigh(h)
    shifted = energies - energies.min()
    weights = np.exp(-shifted / temperature)
    diagonal = np.real(np.einsum("ij,ij->j", vectors.conj(), observable @ vectors))
    return float(np.dot(weights, diagonal) / np.sum(weights))


def run_chain(propagator, observable, strategy, rng, burn_in, samples):
    dimension = propagator.shape[0]
    current_basis = "Z"
    current_index = int(rng.integers(dimension))
    values = []
    parity = []
    total = burn_in + samples
    for iteration in range(total):
        initial = product_state(current_index, current_basis, dimension)
        thermal = propagator @ initial
        thermal /= np.linalg.norm(thermal)
        value = float(np.real(np.vdot(thermal, observable @ thermal)))
        if strategy == "Z_only":
            next_basis = "Z"
        elif strategy == "alternating_ZX":
            next_basis = "X" if current_basis == "Z" else "Z"
        else:
            raise ValueError(strategy)
        amplitudes = thermal if next_basis == "Z" else hadamard_transform(thermal)
        probabilities = np.abs(amplitudes) ** 2
        probabilities /= probabilities.sum()
        current_index = int(rng.choice(dimension, p=probabilities))
        current_basis = next_basis
        if iteration >= burn_in:
            values.append(value)
            parity.append(current_index.bit_count() % 2 if current_basis == "Z" else -1)
    values = np.asarray(values)
    tau, effective = integrated_autocorrelation(values)
    standard_error = float(np.std(values, ddof=1) / np.sqrt(effective))
    z_parities = [p for p in parity if p >= 0]
    return {
        "sample_mean_K": float(np.mean(values)),
        "sample_standard_deviation_K": float(np.std(values, ddof=1)),
        "integrated_autocorrelation_time": tau,
        "effective_sample_count": effective,
        "standard_error_K": standard_error,
        "Z_collapse_even_fraction": (float(np.mean(np.asarray(z_parities) == 0))
                                      if z_parities else None),
    }


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges(3, 2)
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    v = h2 - h0
    rows = []
    for temperature in TEMPERATURES_K:
        propagator = expm(-h2 / (2 * temperature))
        exact = exact_thermal_mean(h2, v, temperature)
        for offset, strategy in enumerate(("Z_only", "alternating_ZX")):
            rng = np.random.default_rng(SEED + int(round(temperature * 1000)) + offset)
            result = run_chain(propagator, v, strategy, rng, BURN_IN, SAMPLES)
            result.update({
                "temperature_K": temperature,
                "strategy": strategy,
                "exact_thermal_mean_V_K": exact,
                "bias_K": result["sample_mean_K"] - exact,
                "bias_per_spin_K": (result["sample_mean_K"] - exact) / n,
                "z_score": ((result["sample_mean_K"] - exact) / result["standard_error_K"]
                            if result["standard_error_K"] > 0 else None),
            })
            rows.append(result)
            print(f"T={temperature:.2f} {strategy}: bias={result['bias_K']:+.3e} K, "
                  f"ESS={result['effective_sample_count']:.0f}", flush=True)

    payload = {
        "purpose": "exact QMETTS Markov-chain and collapse-basis benchmark",
        "cluster": {"N": n, "boundary": "open 3x2 triangular patch", "model": "M2"},
        "observable": "V=M2-M0",
        "field_T": FIELD_T,
        "burn_in": BURN_IN,
        "samples_per_chain": SAMPLES,
        "rows": rows,
        "decision_rule": "alternating Z/X must remove parity trapping and agree with the exact Gibbs mean within 3 autocorrelation-adjusted standard errors",
        "scope_limits": [
            "Imaginary-time propagation is exact; QITE, shot and hardware errors are absent.",
            "The standard error uses a positive-sequence autocorrelation estimate from one chain.",
            "N=6 M2 is a sampling unit test; the scientific target remains N=9 periodic M3.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qmetts_exact_sampling_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
