"""Exact sudden-quench Jarzynski variance screen for N=9 M0-to-M3."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURES_K = (0.095, 0.105)
DERIVATIVE_HALF_WIDTH_K = 0.005
GENERAL_ENTROPY_ERROR = 0.001
STRESS_ENTROPY_ERROR = 0.0002


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def log_partition(energies, temperature):
    return float(logsumexp(-energies / temperature))


def jarzynski_statistics(e0, u0, e1, u1, temperature):
    beta = 1 / temperature
    log_z0 = log_partition(e0, temperature)
    log_p0 = -beta * e0 - log_z0
    overlap = np.abs(u1.conj().T @ u0) ** 2
    with np.errstate(divide="ignore"):
        log_overlap = np.log(overlap)
    work = e1[:, np.newaxis] - e0[np.newaxis, :]
    log_joint = log_overlap + log_p0[np.newaxis, :]
    log_mean_y = float(logsumexp(log_joint - beta * work))
    log_mean_y2 = float(logsumexp(log_joint - 2 * beta * work))
    log_relative_second = log_mean_y2 - 2 * log_mean_y
    relative_variance = float(np.expm1(log_relative_second)) if log_relative_second < 700 else float("inf")
    delta_f = -temperature * log_mean_y
    exact_delta_f = -temperature * (log_partition(e1, temperature) - log_z0)
    return {
        "jarzynski_delta_F_K": delta_f,
        "exact_delta_F_K": exact_delta_f,
        "identity_error_K": delta_f - exact_delta_f,
        "log_relative_second_moment": log_relative_second,
        "relative_variance_exp_betaW": relative_variance,
    }


def adiabatic_statistics(e0_blocks, e1_blocks, temperature):
    """Optimistic sorted level mapping within each conserved parity block."""
    e0 = np.concatenate(e0_blocks)
    e1 = np.concatenate(e1_blocks)
    beta = 1 / temperature
    log_z0 = log_partition(e0, temperature)
    log_p0 = -beta * e0 - log_z0
    work = e1 - e0
    log_mean_y = float(logsumexp(log_p0 - beta * work))
    log_mean_y2 = float(logsumexp(log_p0 - 2 * beta * work))
    log_relative_second = log_mean_y2 - 2 * log_mean_y
    relative_variance = float(np.expm1(log_relative_second)) if log_relative_second < 700 else float("inf")
    delta_f = -temperature * log_mean_y
    exact_delta_f = -temperature * (log_partition(e1, temperature) - log_z0)
    return {
        "jarzynski_delta_F_K": delta_f,
        "exact_delta_F_K": exact_delta_f,
        "identity_error_K": delta_f - exact_delta_f,
        "log_relative_second_moment": log_relative_second,
        "relative_variance_exp_betaW": relative_variance,
    }


def shot_log10(relative_variance, temperature, target_total_F_error):
    # Delta F=-T log(mean Y): first-order variance is
    # T^2 Var(Y)/(M mean(Y)^2).
    if np.isfinite(relative_variance):
        shots = relative_variance * (temperature / target_total_F_error) ** 2
        return float(np.log10(shots)), (int(np.ceil(shots)) if shots < 1e18 else None)
    return float("inf"), None


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    bonds0 = model_bonds(graph, "M0")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        h0 = dense_hamiltonian(n, bonds0, field_k)
        h1 = dense_hamiltonian(n, bonds3, field_k)
        e0, u0 = np.linalg.eigh(h0)
        e1, u1 = np.linalg.eigh(h1)
        e0_blocks = []
        e1_blocks = []
        for parity in (0, 1):
            op0 = SpinHamiltonianOperator(n, bonds0, field_k, parity=parity)
            op1 = SpinHamiltonianOperator(n, bonds3, field_k, parity=parity)
            m0 = np.asarray(op0 @ np.eye(op0.shape[0], dtype=complex))
            m1 = np.asarray(op1 @ np.eye(op1.shape[0], dtype=complex))
            e0_blocks.append(np.linalg.eigvalsh(m0))
            e1_blocks.append(np.linalg.eigvalsh(m1))
        targets = {"general": GENERAL_ENTROPY_ERROR}
        if abs(field_t - 1.45) < 1e-12:
            targets["stress_sign_change"] = STRESS_ENTROPY_ERROR
        for temperature in TEMPERATURES_K:
            protocols = {
                "sudden_quench": jarzynski_statistics(e0, u0, e1, u1, temperature),
                "ideal_sorted_parity_adiabatic": adiabatic_statistics(e0_blocks, e1_blocks, temperature),
            }
            for protocol, stats in protocols.items():
                budgets = {}
                for label, entropy_error in targets.items():
                    # Equal independent Delta-F errors at T-dT and T+dT.
                    target_f_error_per_spin = np.sqrt(2) * DERIVATIVE_HALF_WIDTH_K * entropy_error
                    target_total_f_error = n * target_f_error_per_spin
                    log10_shots, shots = shot_log10(
                        stats["relative_variance_exp_betaW"], temperature,
                        target_total_f_error)
                    budgets[label] = {
                        "target_entropy_error_kB_per_spin": entropy_error,
                        "target_each_delta_F_error_K_per_spin": target_f_error_per_spin,
                        "log10_shots_for_this_temperature": log10_shots,
                        "shots_if_below_1e18": shots,
                    }
                rows.append({"field_T": field_t, "temperature_K": temperature,
                             "protocol": protocol, **stats, "budgets": budgets})
                print(f"B={field_t:.2f}, T={temperature:.3f}, {protocol}: log10 shots="
                      f"{budgets['general']['log10_shots_for_this_temperature']:.2f}", flush=True)

    max_identity_error = max(abs(row["identity_error_K"]) for row in rows)
    assert max_identity_error < 1e-10
    payload = {
        "purpose": "exact sudden-quench two-projective-measurement Jarzynski variance screen",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "protocols": [
            "sudden quench M0 to M3 with exact initial/final energy measurements",
            "optimistic infinitely-slow level mapping sorted within conserved parity sectors",
        ],
        "rows": rows,
        "scope_limits": [
            "The exact Jarzynski identity holds, but sampling can be dominated by exponentially rare work values.",
            "State preparation, energy measurement/phase estimation, circuit depth and hardware noise are excluded.",
            "The first-order Delta-F error formula is used to convert work variance to shots.",
            "This is the exact two-projective-measurement protocol, not the approximate METTS pseudo-work implementation.",
            "The sorted-parity adiabatic mapping is a best-case variance endpoint and has no finite-depth circuit estimate.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "jarzynski_shot_budget_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "max_identity_error_K": max_identity_error,
                      "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
