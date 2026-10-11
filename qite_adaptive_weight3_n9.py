"""Adaptive weight-3 QITE operator selection for correlated N=9 M3 states."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qite_tangent_screen_n9 import (evolve_and_resources,
                                    hadamard_product_state,
                                    labels_up_to_weight,
                                    pauli_action_columns,
                                    spectral_state)
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
DT = 0.05
CORRELATION_TIMES = (0.5, 2.0)
SELECTION_CHECKPOINTS = (64, 128, 256)
BATCH = 32


def real_design(labels, state):
    columns = -1j * pauli_action_columns(labels, state)
    return np.vstack((columns.real, columns.imag))


def solve(design, target):
    coefficients, _, rank, singular = np.linalg.lstsq(design, target, rcond=1e-10)
    residual = design @ coefficients - target
    nonzero = singular[singular > singular[0] * 1e-10]
    return coefficients, residual, int(rank), float(nonzero[0] / nonzero[-1])


def generator_from(labels, coefficients):
    keep = np.abs(coefficients) > 1e-10
    return SparsePauliOp(
        [label for label, use in zip(labels, keep, strict=True) if use],
        coefficients[keep],
    ).simplify(atol=1e-10)


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds = model_bonds(graph, "M3", boundary="tin_foil")
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    h = np.asarray(operator @ np.eye(1 << n, dtype=complex))
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    initial = hadamard_product_state(n, 0b101010101)
    labels2 = labels_up_to_weight(n, 2)
    all_to_3 = labels_up_to_weight(n, 3)
    labels3 = all_to_3[len(labels2):]
    rows = []

    for tau in CORRELATION_TIMES:
        state = spectral_state(eigenvalues, eigenvectors, initial, tau)
        energy = float(np.vdot(state, h @ state).real)
        tangent = -(h @ state - energy * state)
        target = np.concatenate((tangent.real, tangent.imag))
        design2 = real_design(labels2, state)
        design3 = real_design(labels3, state)
        selected = []
        available = np.ones(len(labels3), dtype=bool)
        while len(selected) < max(SELECTION_CHECKPOINTS):
            current = (design2 if not selected else
                       np.column_stack((design2, design3[:, selected])))
            coefficients, residual, _, _ = solve(current, target)
            scores = np.abs(design3.T @ residual)
            scores[~available] = -np.inf
            add_count = min(BATCH, max(SELECTION_CHECKPOINTS) - len(selected))
            candidates = np.argpartition(scores, -add_count)[-add_count:]
            candidates = candidates[np.argsort(scores[candidates])[::-1]]
            selected.extend(int(index) for index in candidates)
            available[candidates] = False

            if len(selected) in SELECTION_CHECKPOINTS:
                labels = labels2 + [labels3[index] for index in selected]
                design = np.column_stack((design2, design3[:, selected]))
                coefficients, residual, rank, condition = solve(design, target)
                generator = generator_from(labels, coefficients)
                exact_next = spectral_state(eigenvalues, eigenvectors, state, DT)
                evolved, resources = evolve_and_resources(generator, state, DT)
                infidelity = 1 - float(abs(np.vdot(exact_next, evolved)) ** 2)
                row = {
                    "correlation_time_inverse_K": tau,
                    "imaginary_step_inverse_K": DT,
                    "base_weight2_pool_size": len(labels2),
                    "weight3_candidate_count": len(labels3),
                    "selected_weight3_count": len(selected),
                    "combined_pool_size": len(labels),
                    "effective_rank": rank,
                    "retained_pauli_terms": int(len(generator)),
                    "relative_tangent_residual": float(np.linalg.norm(residual) / np.linalg.norm(target)),
                    "resolved_condition_number": condition,
                    "one_step_infidelity": infidelity,
                    **resources,
                    "projected_CX_per_100mK_METTS_state": int(
                        np.ceil(5.0 / DT) * resources["compiled_all_to_all_CX"]),
                }
                rows.append(row)
                print(f"tau={tau:.1f}, +w3={len(selected)}: "
                      f"residual={row['relative_tangent_residual']:.3e}, "
                      f"infidelity={infidelity:.3e}, CX={resources['compiled_all_to_all_CX']}", flush=True)

    payload = {
        "purpose": "adaptive weight-3 QITE accuracy-resource Pareto for N=9 M3",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "field_T": FIELD_T,
        "initial_state": "X-basis product state X:101010101",
        "selection": "greedy batches ranked by absolute overlap with the current real tangent residual",
        "rows": rows,
        "scope_limits": [
            "Operator selection and coefficients use exact statevectors with no measurement noise.",
            "Selection is greedy and is not guaranteed to be the smallest possible operator set.",
            "Only one step at two correlated states is tested; full thermal trajectories remain unproven.",
            "CX counts assume all-to-all connectivity and exclude routing overhead.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qite_adaptive_weight3_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
