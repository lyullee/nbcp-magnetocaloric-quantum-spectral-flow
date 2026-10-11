"""N=9 periodic M3 one-step QITE locality and resource screen."""
from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path
import time

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.synthesis import LieTrotter

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELD_T = 1.45
DT = 0.05
CORRELATION_TIMES = (0.0, 0.5, 2.0)
POOL_WEIGHTS = (1, 2)


def labels_up_to_weight(n, maximum_weight):
    labels = []
    for weight in range(1, maximum_weight + 1):
        for sites in combinations(range(n), weight):
            for axes in product("XYZ", repeat=weight):
                label = ["I"] * n
                for site, axis in zip(sites, axes, strict=True):
                    label[n - 1 - site] = axis
                labels.append("".join(label))
    return labels


def pauli_action_columns(labels, state):
    """Return columns P|state> without constructing dense Pauli matrices."""
    dimension = len(state)
    source = np.arange(dimension, dtype=np.uint64)
    columns = np.empty((dimension, len(labels)), dtype=complex)
    for column, label in enumerate(labels):
        x_mask = 0
        z_mask = 0
        y_count = 0
        for label_position, axis in enumerate(label):
            site = len(label) - 1 - label_position
            if axis in "XY":
                x_mask |= 1 << site
            if axis in "YZ":
                z_mask |= 1 << site
            if axis == "Y":
                y_count += 1
        parity = np.fromiter(((int(index) & z_mask).bit_count() & 1 for index in source),
                             dtype=np.int8, count=dimension)
        phase = (1j ** y_count) * (1 - 2 * parity)
        target = source ^ np.uint64(x_mask)
        columns[target.astype(np.int64), column] = phase * state
    return columns


def fit_generator(h, state, labels):
    energy = float(np.vdot(state, h @ state).real)
    tangent = -(h @ state - energy * state)
    pauli_columns = pauli_action_columns(labels, state)
    columns = -1j * pauli_columns
    design = np.vstack((columns.real, columns.imag))
    target = np.concatenate((tangent.real, tangent.imag))
    coefficients, _, rank, singular = np.linalg.lstsq(design, target, rcond=1e-10)
    residual = design @ coefficients - target
    keep = np.abs(coefficients) > 1e-10
    generator = SparsePauliOp(
        [label for label, use in zip(labels, keep, strict=True) if use],
        coefficients[keep],
    ).simplify(atol=1e-10)
    nonzero = singular[singular > singular[0] * 1e-10]
    return generator, {
        "energy_K": energy,
        "relative_tangent_residual": float(np.linalg.norm(residual) / np.linalg.norm(target)),
        "effective_rank": int(rank),
        "retained_pauli_terms": int(len(generator)),
        "resolved_condition_number": float(nonzero[0] / nonzero[-1]),
    }


def evolve_and_resources(generator, state, dt):
    n = generator.num_qubits
    circuit = QuantumCircuit(n)
    circuit.append(PauliEvolutionGate(generator, time=dt,
                                      synthesis=LieTrotter(reps=1)), range(n))
    evolved = Statevector(state).evolve(circuit).data
    compiled = transpile(circuit, basis_gates=["rz", "sx", "x", "cx"],
                         optimization_level=3, seed_transpiler=20261011)
    counts = {str(key): int(value) for key, value in compiled.count_ops().items()}
    return evolved, {
        "compiled_all_to_all_depth": int(compiled.depth()),
        "compiled_all_to_all_CX": int(counts.get("cx", 0)),
        "compiled_operation_counts": counts,
    }


def hadamard_product_state(n, index):
    signs = np.ones(1 << n)
    for basis_index in range(1 << n):
        if (basis_index & index).bit_count() & 1:
            signs[basis_index] = -1
    return signs.astype(complex) / np.sqrt(1 << n)


def spectral_state(eigenvalues, eigenvectors, initial, tau):
    coefficients = eigenvectors.conj().T @ initial
    state = eigenvectors @ (np.exp(-tau * eigenvalues) * coefficients)
    return state / np.linalg.norm(state)


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
    pools = {weight: labels_up_to_weight(n, weight) for weight in POOL_WEIGHTS}

    # Validate the efficient Pauli action convention against Qiskit.
    check_labels = pools[2][::max(1, len(pools[2]) // 11)][:11]
    fast = pauli_action_columns(check_labels, initial)
    reference = np.column_stack([SparsePauliOp(label).to_matrix() @ initial
                                 for label in check_labels])
    action_error = float(np.max(np.abs(fast - reference)))
    assert action_error < 1e-12

    rows = []
    for tau in CORRELATION_TIMES:
        state = spectral_state(eigenvalues, eigenvectors, initial, tau)
        exact_next = spectral_state(eigenvalues, eigenvectors, state, DT)
        for weight in POOL_WEIGHTS:
            labels = pools[weight]
            generator, fit = fit_generator(h, state, labels)
            evolved, resources = evolve_and_resources(generator, state, DT)
            infidelity = 1 - float(abs(np.vdot(exact_next, evolved)) ** 2)
            row = {
                "correlation_time_inverse_K": tau,
                "imaginary_step_inverse_K": DT,
                "maximum_pauli_weight": weight,
                "pool_size": len(labels),
                **fit,
                **resources,
                "one_step_infidelity": infidelity,
                "projected_steps_to_100mK_beta_over_2": int(np.ceil(5.0 / DT)),
                "projected_CX_per_100mK_METTS_state": int(np.ceil(5.0 / DT) * resources["compiled_all_to_all_CX"]),
            }
            rows.append(row)
            print(f"tau={tau:.1f}, weight<={weight}: residual={fit['relative_tangent_residual']:.3e}, "
                  f"infidelity={infidelity:.3e}, CX={resources['compiled_all_to_all_CX']}", flush=True)

    payload = {
        "purpose": "N=9 periodic M3 one-step QITE locality and optimistic circuit-resource screen",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "field_T": FIELD_T,
        "initial_state": "X-basis product state X:101010101",
        "efficient_pauli_action_max_error": action_error,
        "rows": rows,
        "scope_limits": [
            "QITE coefficients use exact statevectors and exclude measurement noise.",
            "Only one imaginary-time step is tested at each pre-correlated state.",
            "CX counts assume all-to-all connectivity; hardware routing can only increase them.",
            "Projected full-path CX is a step-count product, not a compiled concatenated circuit.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "qite_tangent_screen_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
