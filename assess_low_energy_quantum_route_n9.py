"""Screen a low-energy-spectrum route to N=9 NBCP thermodynamics.

Instead of estimating entropy directly on a quantum processor, this route
measures a small number of low-energy levels and reconstructs Z, F and S
classically.  Exact diagonalization is used only to decide whether that route
is scientifically viable and to produce optimistic measurement budgets.
"""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from scipy.special import gammaincc, logsumexp

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import spin_hamiltonian_to_sparse_pauli
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURES_K = (0.09, 0.10, 0.12, 0.30)
GENERAL_DELTA_S_ERROR = 0.001
STRESS_DELTA_S_ERROR = 0.0002
MASS_LEVELS = (0.95, 0.99, 0.999, 0.9999)
TAYLOR_FIDELITIES = (0.99, 0.999, 0.9999)
MAX_TAYLOR_ORDER = 2000


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def thermal_distribution(energies, temperature):
    log_weight = -(energies - energies[0]) / temperature
    probability = np.exp(log_weight - logsumexp(log_weight))
    entropy = float(-np.dot(probability[probability > 0],
                            np.log(probability[probability > 0])))
    return probability, entropy


def truncated_entropy(energies, temperature, rank):
    probability, entropy = thermal_distribution(energies[:rank], temperature)
    return probability, entropy


def mass_rank(probability, target):
    return int(np.searchsorted(np.cumsum(probability), target) + 1)


def low_rank_fidelity_ceiling(probability, ancilla_qubits):
    """Maximum squared Uhlmann fidelity for rank <= 2**ancilla_qubits."""
    return float(np.sum(probability[:min(1 << ancilla_qubits, len(probability))]))


def shifted_positive_taylor_distribution(energies, temperature, order):
    """Best-conditioned ordinary Taylor screen using Emax-H >= 0.

    exp[-beta(H-Emax)] is expanded with positive terms.  The scalar energy
    shift cancels after normalization.  gammaincc evaluates the truncated
    exponential without explicitly summing huge powers.
    """
    y = (energies[-1] - energies) / temperature
    cdf = gammaincc(order + 1, y)
    # At very small cdf scipy can underflow.  Such a component is immaterial
    # after normalization; retain -inf rather than manufacturing probability.
    with np.errstate(divide="ignore"):
        log_weight = y + np.log(cdf)
    log_weight -= logsumexp(log_weight)
    return np.exp(log_weight)


def classical_fidelity(probability, approximation):
    return float(np.sum(np.sqrt(probability * approximation)) ** 2)


def minimum_taylor_orders(energies, temperature, probability):
    remaining = set(TAYLOR_FIDELITIES)
    answer = {}
    for order in range(MAX_TAYLOR_ORDER + 1):
        approximation = shifted_positive_taylor_distribution(
            energies, temperature, order)
        fidelity = classical_fidelity(probability, approximation)
        for target in tuple(remaining):
            if fidelity >= target:
                answer[str(target)] = {
                    "minimum_order": order,
                    "achieved_fidelity": fidelity,
                }
                remaining.remove(target)
        if not remaining:
            break
    for target in remaining:
        answer[str(target)] = {
            "minimum_order": None,
            "achieved_fidelity_at_max_order": classical_fidelity(
                probability,
                shifted_positive_taylor_distribution(
                    energies, temperature, MAX_TAYLOR_ORDER)),
        }
    return answer


def qwc_group_amplitudes(groups, eigenvectors, state_count):
    """sum_g sqrt(Var_{|E_k>}(G_g)) for each low-energy eigenstate."""
    amplitudes = np.zeros(state_count)
    for group in groups:
        matrix = np.asarray(group.to_matrix())
        vectors = eigenvectors[:, :state_count]
        action = matrix @ vectors
        mean = np.real(np.sum(vectors.conj() * action, axis=0))
        second = np.sum(np.abs(action) ** 2, axis=0).real
        amplitudes += np.sqrt(np.maximum(0.0, second - mean**2))
    return amplitudes


def entropy_energy_gradient(energies, temperature):
    probability, _ = thermal_distribution(energies, temperature)
    mean_energy = float(np.dot(probability, energies))
    # dS/dE_j in k_B per kelvin when H/k_B is expressed in kelvin.
    gradient = -(probability * (energies - mean_energy)) / temperature**2
    return probability, gradient


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    model_bond_sets = {
        "M0": model_bonds(graph, "M0", boundary="tin_foil"),
        "M3": model_bonds(graph, "M3", boundary="tin_foil"),
    }
    rows = []
    comparison_rows = []

    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        solved = {}
        for model, bonds in model_bond_sets.items():
            hamiltonian = dense_hamiltonian(n, bonds, field_k)
            energies, vectors = np.linalg.eigh(hamiltonian)
            pauli = spin_hamiltonian_to_sparse_pauli(n, bonds, field_k)
            groups = pauli.group_commuting(qubit_wise=True)
            solved[model] = {
                "energies": energies,
                "vectors": vectors,
                "pauli_terms": int(len(pauli)),
                "groups": groups,
            }

            for temperature in TEMPERATURES_K:
                probability, entropy = thermal_distribution(energies, temperature)
                mass_ranks = {str(target): mass_rank(probability, target)
                              for target in MASS_LEVELS}
                ancilla_ceilings = {
                    str(q): low_rank_fidelity_ceiling(probability, q)
                    for q in range(0, n + 1)
                }
                rows.append({
                    "model": model,
                    "field_T": field_t,
                    "temperature_K": temperature,
                    "energy_min_K": float(energies[0]),
                    "energy_max_K": float(energies[-1]),
                    "spectral_width_K": float(energies[-1] - energies[0]),
                    "log10_gibbs_condition_number": float(
                        (energies[-1] - energies[0]) /
                        (temperature * np.log(10))),
                    "entropy_kB": entropy,
                    "effective_rank_exp_entropy": float(np.exp(entropy)),
                    "ground_state_probability": float(probability[0]),
                    "mass_ranks": mass_ranks,
                    "minimum_ancilla_qubits_for_99_9pct_mass": int(
                        np.ceil(np.log2(mass_ranks["0.999"]))),
                    "low_rank_fidelity_ceilings": ancilla_ceilings,
                    "positive_shift_taylor_orders": minimum_taylor_orders(
                        energies, temperature, probability),
                })

        for temperature in TEMPERATURES_K:
            target = (STRESS_DELTA_S_ERROR
                      if abs(field_t - 1.45) < 1e-12 and
                      abs(temperature - 0.09) < 1e-12
                      else GENERAL_DELTA_S_ERROR)
            exact = {}
            for model in ("M0", "M3"):
                _, exact[model] = thermal_distribution(
                    solved[model]["energies"], temperature)

            chosen_rank = None
            rank_detail = None
            for rank in range(1, (1 << n) + 1):
                truncated = {}
                for model in ("M0", "M3"):
                    _, truncated[model] = truncated_entropy(
                        solved[model]["energies"], temperature, rank)
                direct_error = ((truncated["M3"] - truncated["M0"]) -
                                (exact["M3"] - exact["M0"])) / n
                triangle_bound = (abs(truncated["M3"] - exact["M3"]) +
                                  abs(truncated["M0"] - exact["M0"])) / n
                if triangle_bound <= target:
                    chosen_rank = rank
                    rank_detail = {
                        "truncated_delta_S_kB_per_spin":
                            (truncated["M3"] - truncated["M0"]) / n,
                        "direct_delta_S_error_kB_per_spin": direct_error,
                        "triangle_error_bound_kB_per_spin": triangle_bound,
                        "individual_entropy_kB": truncated,
                    }
                    break
            assert chosen_rank is not None

            # Optimistic measurement budget: exact low-energy eigenstates,
            # QWC Hamiltonian groups, optimal allocation over states/groups.
            weighted_amplitudes = []
            per_model = {}
            for model in ("M0", "M3"):
                energies = solved[model]["energies"][:chosen_rank]
                _, gradient = entropy_energy_gradient(energies, temperature)
                group_amplitude = qwc_group_amplitudes(
                    solved[model]["groups"], solved[model]["vectors"],
                    chosen_rank)
                coefficients = np.abs(gradient) / n
                contribution = coefficients * group_amplitude
                weighted_amplitudes.extend(contribution.tolist())
                per_model[model] = {
                    "pauli_term_count": solved[model]["pauli_terms"],
                    "qwc_group_count": len(solved[model]["groups"]),
                    "entropy_gradient_abs_sum_inverse_K":
                        float(np.sum(np.abs(gradient))),
                    "group_stddev_sum_by_state_K": group_amplitude.tolist(),
                    "weighted_measurement_amplitude_sum":
                        float(np.sum(contribution)),
                }
            total_amplitude = float(np.sum(weighted_amplitudes))
            total_shots = int(np.ceil((total_amplitude / target) ** 2))
            exact_delta_s = (exact["M3"] - exact["M0"]) / n
            comparison_rows.append({
                "field_T": field_t,
                "temperature_K": temperature,
                "target_delta_S_error_kB_per_spin": target,
                "exact_delta_S_kB_per_spin": exact_delta_s,
                "minimum_common_low_energy_rank": chosen_rank,
                "minimum_ancilla_qubits_for_that_rank": int(
                    np.ceil(np.log2(chosen_rank))),
                **rank_detail,
                "optimistic_exact_eigenstate_energy_measurement": {
                    "optimal_total_shots": total_shots,
                    "weighted_measurement_amplitude_sum": total_amplitude,
                    "models": per_model,
                },
            })
            print(f"B={field_t:.2f} T={temperature:.2f}: "
                  f"K={chosen_rank}, shots={total_shots:,}", flush=True)

    payload = {
        "purpose": "screen low-energy eigenspectrum reconstruction as an alternative to direct entropy estimation",
        "cluster": {"N": n, "translations": graph["translations"],
                    "boundary": "3x3 periodic; M3 tin-foil Ewald"},
        "principle": "measure K low-energy levels; reconstruct Z, F and S classically",
        "rows": rows,
        "delta_S_comparisons": comparison_rows,
        "literature_basis": [
            "https://arxiv.org/abs/2203.12757",
            "https://arxiv.org/abs/2303.11276",
        ],
        "scope_limits": [
            "Exact diagonalization supplies the eigenstates and spectra only for local route screening.",
            "Shot budgets assume exact preparation of every selected eigenstate and omit VQD/QSE optimization, overlap penalties, routing, noise and mitigation.",
            "The entropy-error rank uses a conservative triangle bound, not cancellation between M0 and M3 truncation errors.",
            "The ordinary Taylor screen uses the most stable scalar shift Emax-H; implementing H^m still incurs rapidly growing observable cost.",
            "N=9 finite-cluster low-rank behavior does not prove the same rank suffices in the bulk or at N=12/15/18.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "low_energy_quantum_route_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "elapsed_seconds": payload["elapsed_seconds"],
        "comparisons": [{
            "field_T": row["field_T"],
            "temperature_K": row["temperature_K"],
            "rank": row["minimum_common_low_energy_rank"],
            "shots": row["optimistic_exact_eigenstate_energy_measurement"]["optimal_total_shots"],
        } for row in comparison_rows],
    }, indent=2))


if __name__ == "__main__":
    main()
