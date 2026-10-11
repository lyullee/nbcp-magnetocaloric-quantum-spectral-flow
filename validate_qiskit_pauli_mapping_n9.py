"""Validate and count the N=9 periodic M0-to-M3 Qiskit observables."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import qiskit
from qiskit.quantum_info import Statevector

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import spin_hamiltonian_to_sparse_pauli
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45


def dense_matrix(operator):
    eye = np.eye(operator.shape[0], dtype=complex)
    return np.asarray(operator @ eye)


def operator_summary(pauli):
    return {
        "pauli_term_count": int(len(pauli)),
        "global_commuting_group_count": len(pauli.group_commuting(qubit_wise=False)),
        "qubitwise_commuting_group_count": len(pauli.group_commuting(qubit_wise=True)),
        "coefficient_l1_K": float(np.sum(np.abs(pauli.coeffs))),
        "largest_coefficient_K": float(np.max(np.abs(pauli.coeffs))),
    }


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    rng = np.random.default_rng(20261011)
    vector = rng.normal(size=1 << n) + 1j * rng.normal(size=1 << n)
    vector /= np.linalg.norm(vector)
    state = Statevector(vector)

    matrices = {}
    paulis = {}
    models = {}
    for model in ("M0", "M3"):
        bonds = model_bonds(graph, model, boundary="tin_foil")
        exact = dense_matrix(SpinHamiltonianOperator(n, bonds, field_k))
        pauli = spin_hamiltonian_to_sparse_pauli(n, bonds, field_k)
        matrix_error = float(np.max(np.abs(pauli.to_matrix() - exact)))
        expected = complex(np.vdot(vector, exact @ vector))
        measured = complex(state.expectation_value(pauli))
        models[model] = {
            "raw_bond_contribution_count": len(bonds),
            **operator_summary(pauli),
            "dense_matrix_max_abs_error_K": matrix_error,
            "random_state_expectation_K": measured.real,
            "random_state_expectation_abs_error_K": float(abs(measured - expected)),
        }
        matrices[model] = exact
        paulis[model] = pauli
        print(f"mapped {model}: {len(pauli)} Pauli terms", flush=True)

    exact_v = matrices["M3"] - matrices["M0"]
    pauli_v = (paulis["M3"] - paulis["M0"]).simplify(atol=1e-13)
    expected_v = complex(np.vdot(vector, exact_v @ vector))
    measured_v = complex(state.expectation_value(pauli_v))
    correction = {
        **operator_summary(pauli_v),
        "dense_matrix_max_abs_error_K": float(np.max(np.abs(pauli_v.to_matrix() - exact_v))),
        "random_state_expectation_K": measured_v.real,
        "random_state_expectation_abs_error_K": float(abs(measured_v - expected_v)),
    }

    tolerance = 1e-11
    assert max(row["dense_matrix_max_abs_error_K"] for row in models.values()) < tolerance
    assert max(row["random_state_expectation_abs_error_K"] for row in models.values()) < tolerance
    assert correction["dense_matrix_max_abs_error_K"] < tolerance
    assert correction["random_state_expectation_abs_error_K"] < tolerance

    payload = {
        "purpose": "N=9 periodic scientific observable mapping before thermal-state simulation",
        "runtime": {"qiskit": qiskit.__version__, "simulator": "qiskit.quantum_info.Statevector"},
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "field_T": FIELD_T,
        "field_K": field_k,
        "basis_mapping": {"Sx": "X/2", "Sy": "-Y/2", "Sz": "-Z/2"},
        "models": models,
        "correction_M3_minus_M0": correction,
        "acceptance_tolerance_K": tolerance,
        "status": "PASS",
        "scope_limits": [
            "This validates the observable representation and measurement grouping, not Gibbs-state preparation.",
            "Commuting-group counts are Qiskit greedy partitions and may not be globally minimal.",
            "N=9 remains a finite-cluster algorithm benchmark rather than a bulk-material prediction.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qiskit_pauli_mapping_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
