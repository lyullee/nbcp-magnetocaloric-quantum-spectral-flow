"""Finite-time exact ramp screen for N=9 TPM Jarzynski at 1.45 T."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter
from scipy.special import logsumexp

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import spin_hamiltonian_to_sparse_pauli
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
TEMPERATURES_K = (0.095, 0.105)
RAMP_TIMES_INVERSE_K = (0.25, 0.50, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0)
SLICES = 24
DERIVATIVE_HALF_WIDTH_K = 0.005
GENERAL_ENTROPY_ERROR = 0.001
STRESS_ENTROPY_ERROR = 0.0002


def dense_block(n, bonds, field_k, parity):
    operator = SpinHamiltonianOperator(n, bonds, field_k, parity=parity)
    return np.asarray(operator @ np.eye(operator.shape[0], dtype=complex))


def propagate(path_eigensystems, ramp_time, dimension):
    unitary = np.eye(dimension, dtype=complex)
    dt = ramp_time / len(path_eigensystems)
    for energies, vectors in path_eigensystems:
        step = (vectors * np.exp(-1j * energies * dt)) @ vectors.conj().T
        unitary = step @ unitary
    return unitary


def work_stats(blocks, transition_blocks, temperature):
    all_e0 = np.concatenate([blocks[p][0] for p in (0, 1)])
    log_z0 = float(logsumexp(-all_e0 / temperature))
    beta = 1 / temperature
    log_terms_1 = []
    log_terms_2 = []
    survival = 0.0
    for parity in (0, 1):
        e0, e1 = blocks[parity]
        transition = transition_blocks[parity]
        probability = np.abs(transition) ** 2
        with np.errstate(divide="ignore"):
            log_probability = np.log(probability)
        log_p0 = -beta * e0 - log_z0
        work = e1[:, None] - e0[None, :]
        log_joint = log_probability + log_p0[None, :]
        log_terms_1.append((log_joint - beta * work).ravel())
        log_terms_2.append((log_joint - 2 * beta * work).ravel())
        survival += float(np.dot(np.exp(log_p0), np.diag(probability)))
    log_mean = float(logsumexp(np.concatenate(log_terms_1)))
    log_second = float(logsumexp(np.concatenate(log_terms_2)))
    relative_variance = float(np.expm1(log_second - 2 * log_mean))
    delta_f = -temperature * log_mean
    all_e1 = np.concatenate([blocks[p][1] for p in (0, 1)])
    exact_delta_f = -temperature * (
        logsumexp(-all_e1 / temperature) - log_z0)
    return delta_f, exact_delta_f, relative_variance, survival


def shot_budget(relative_variance, temperature, entropy_error, n):
    target_f_error = n * np.sqrt(2) * DERIVATIVE_HALF_WIDTH_K * entropy_error
    shots = relative_variance * (temperature / target_f_error) ** 2
    return {"log10_shots": float(np.log10(shots)),
            "shots": int(np.ceil(shots)) if shots < 1e18 else None}


def slice_resource(n, bonds0, bonds3, field_k):
    p0 = spin_hamiltonian_to_sparse_pauli(n, bonds0, field_k)
    p1 = spin_hamiltonian_to_sparse_pauli(n, bonds3, field_k)
    midpoint = (0.5 * p0 + 0.5 * p1).simplify(atol=1e-13)
    circuit = QuantumCircuit(n)
    circuit.append(PauliEvolutionGate(midpoint, time=1 / SLICES,
                                      synthesis=LieTrotter(reps=1)), range(n))
    compiled = transpile(circuit, basis_gates=["rz", "sx", "x", "cx"],
                         optimization_level=3, seed_transpiler=20261011)
    counts = {str(key): int(value) for key, value in compiled.count_ops().items()}
    return {"midpoint_pauli_terms": int(len(midpoint)),
            "all_to_all_CX_per_slice": int(counts.get("cx", 0)),
            "all_to_all_depth_per_slice": int(compiled.depth()),
            "projected_CX_for_ramp": int(SLICES * counts.get("cx", 0))}


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds0 = model_bonds(graph, "M0")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    blocks = {}
    paths = {}
    for parity in (0, 1):
        h0 = dense_block(n, bonds0, field_k, parity)
        h1 = dense_block(n, bonds3, field_k, parity)
        e0, u0 = np.linalg.eigh(h0)
        e1, u1 = np.linalg.eigh(h1)
        blocks[parity] = (e0, e1, u0, u1)
        paths[parity] = []
        for step in range(SLICES):
            lam = (step + 0.5) / SLICES
            paths[parity].append(np.linalg.eigh((1 - lam) * h0 + lam * h1))
    resource = slice_resource(n, bonds0, bonds3, field_k)

    rows = []
    for ramp_time in RAMP_TIMES_INVERSE_K:
        transitions = {}
        energy_blocks = {}
        for parity in (0, 1):
            e0, e1, u0, u1 = blocks[parity]
            unitary = propagate(paths[parity], ramp_time, len(e0))
            transitions[parity] = u1.conj().T @ unitary @ u0
            energy_blocks[parity] = (e0, e1)
        for temperature in TEMPERATURES_K:
            delta_f, exact_delta_f, relative_variance, survival = work_stats(
                energy_blocks, transitions, temperature)
            budgets = {
                "general": shot_budget(relative_variance, temperature,
                                       GENERAL_ENTROPY_ERROR, n),
                "stress_sign_change": shot_budget(relative_variance, temperature,
                                                  STRESS_ENTROPY_ERROR, n),
            }
            row = {
                "ramp_time_inverse_K": ramp_time,
                "temperature_K": temperature,
                "jarzynski_delta_F_K": delta_f,
                "exact_delta_F_K": exact_delta_f,
                "identity_error_K": delta_f - exact_delta_f,
                "relative_variance_exp_betaW": relative_variance,
                "thermal_same_rank_survival_probability": survival,
                "budgets": budgets,
            }
            rows.append(row)
            print(f"ramp={ramp_time:.2f}, T={temperature:.3f}: survival={survival:.3f}, "
                  f"log10shots={budgets['general']['log10_shots']:.2f}", flush=True)

    payload = {
        "purpose": "finite-time exact-dynamics Jarzynski shot-depth tradeoff at the sign-change field",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "field_T": FIELD_T,
        "ramp": {"path": "linear M0 to M3", "piecewise_constant_exact_slices": SLICES,
                 "times_inverse_K": list(RAMP_TIMES_INVERSE_K)},
        "optimistic_first_order_Pauli_slice_resource": resource,
        "rows": rows,
        "scope_limits": [
            "Ramp dynamics use exact exponentials of 24 midpoint Hamiltonians.",
            "The circuit resource uses one first-order Pauli-Trotter slice and does not reproduce the exact-slice error.",
            "State preparation and two projective energy measurements/phase estimation are excluded.",
            "All-to-all CX counts exclude hardware routing and fault-tolerance overhead.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "jarzynski_ramp_screen_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "resource": resource,
                      "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
