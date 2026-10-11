"""Local N=6 statevector validation of the Qiskit Hamiltonian mapping."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import qiskit
from qiskit.quantum_info import Statevector

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds,
                                  spin_hamiltonian_to_sparse_pauli)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45


def dense_matrix(operator):
    eye = np.eye(operator.shape[0], dtype=complex)
    return np.asarray(operator @ eye)


def main():
    started = time.perf_counter()
    graph = open_triangular_directional_edges(3, 2)
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    rng = np.random.default_rng(20261011)
    vector = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
    vector /= np.linalg.norm(vector)
    state = Statevector(vector)

    models = {}
    matrices = {}
    paulis = {}
    for model in ("M0", "M2"):
        bonds = open_weak_model_bonds(graph, model)
        exact = dense_matrix(SpinHamiltonianOperator(n, bonds, field_k))
        pauli = spin_hamiltonian_to_sparse_pauli(n, bonds, field_k)
        qiskit_matrix = pauli.to_matrix()
        matrix_error = float(np.max(np.abs(qiskit_matrix - exact)))
        exact_expectation = complex(np.vdot(vector, exact @ vector))
        qiskit_expectation = complex(state.expectation_value(pauli))
        expectation_error = float(abs(qiskit_expectation - exact_expectation))
        groups_global = pauli.group_commuting(qubit_wise=False)
        groups_qubitwise = pauli.group_commuting(qubit_wise=True)
        models[model] = {
            "bond_count": len(bonds),
            "pauli_term_count": int(len(pauli)),
            "global_commuting_group_count": len(groups_global),
            "qubitwise_commuting_group_count": len(groups_qubitwise),
            "dense_matrix_max_abs_error_K": matrix_error,
            "random_state_expectation_K": qiskit_expectation.real,
            "random_state_expectation_abs_error_K": expectation_error,
        }
        matrices[model] = exact
        paulis[model] = pauli

    exact_v = matrices["M2"] - matrices["M0"]
    pauli_v = (paulis["M2"] - paulis["M0"]).simplify(atol=1e-13)
    v_matrix_error = float(np.max(np.abs(pauli_v.to_matrix() - exact_v)))
    exact_v_expectation = complex(np.vdot(vector, exact_v @ vector))
    pauli_v_expectation = complex(state.expectation_value(pauli_v))
    v_expectation_error = float(abs(pauli_v_expectation - exact_v_expectation))

    tolerance = 1e-12
    assert max(row["dense_matrix_max_abs_error_K"] for row in models.values()) < tolerance
    assert max(row["random_state_expectation_abs_error_K"] for row in models.values()) < tolerance
    assert v_matrix_error < tolerance
    assert v_expectation_error < tolerance

    payload = {
        "purpose": "N=6 open-patch circuit-development unit test; not a material or bulk result",
        "runtime": {"qiskit": qiskit.__version__, "simulator": "qiskit.quantum_info.Statevector"},
        "cluster": graph,
        "field_T": FIELD_T,
        "field_K": field_k,
        "basis_mapping": {"Sx": "X/2", "Sy": "-Y/2", "Sz": "-Z/2"},
        "models": models,
        "correction_M2_minus_M0": {
            "pauli_term_count": int(len(pauli_v)),
            "global_commuting_group_count": len(pauli_v.group_commuting(qubit_wise=False)),
            "qubitwise_commuting_group_count": len(pauli_v.group_commuting(qubit_wise=True)),
            "dense_matrix_max_abs_error_K": v_matrix_error,
            "random_state_expectation_K": pauli_v_expectation.real,
            "random_state_expectation_abs_error_K": v_expectation_error,
        },
        "acceptance_tolerance_K": tolerance,
        "status": "PASS",
        "scope_limits": [
            "The open N=6 patch uses M2 nearest-neighbour corrections because periodic M3 is ill-defined at N=6.",
            "The scientific algorithm benchmark remains the N=9 periodic M0-to-M3 target.",
            "Statevector expectation agreement validates mapping only, not thermal-state preparation.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qiskit_pauli_mapping_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
