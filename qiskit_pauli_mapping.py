"""Map the project's spin Hamiltonian convention to Qiskit Pauli operators.

The matrix-free code uses the one-site basis |down>, |up>.  Relative to
Qiskit's |0>, |1> Pauli matrices this means Sx=X/2, Sy=-Y/2 and Sz=-Z/2.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from matrix_free_hamiltonian import SpinBond
from screen_weak_terms import (DELTA_XY_OVER_ZZ, GXY, GZ, JPP_EXTRA_MEV,
                               JZZ_MEV, KD_MEV, MEV_TO_K)


PAULI_BY_AXIS = ("X", "Y", "Z")
SPIN_SIGNS = (1.0, -1.0, -1.0)


def open_triangular_directional_edges(nx: int = 3, ny: int = 2):
    """Return a rectangular open patch in triangular-lattice coordinates."""
    coordinates = [(x, y) for y in range(ny) for x in range(nx)]
    index = {coordinate: site for site, coordinate in enumerate(coordinates)}
    directions = ((1, 0, 0.0), (0, 1, np.pi / 3), (1, -1, -np.pi / 3))
    edges = []
    for i, (x, y) in enumerate(coordinates):
        for dx, dy, theta in directions:
            target = (x + dx, y + dy)
            if target in index:
                edges.append((i, index[target], float(theta)))
    return {"N": len(coordinates), "coordinates": coordinates, "directional_edges": edges}


def open_weak_model_bonds(graph, model: str):
    """M0--M2 Cartesian bonds on an open patch; intended for circuit tests."""
    jxy = JZZ_MEV * DELTA_XY_OVER_ZZ
    jzz = JZZ_MEV
    jpp = 0.0 if model == "M0" else JPP_EXTRA_MEV
    if model == "M2":
        jxy -= KD_MEV * GXY**2 / 2
        jzz += KD_MEV * GZ**2
        jpp -= 3 * KD_MEV * GXY**2 / 4
    if model not in {"M0", "M1", "M2"}:
        raise ValueError("the open circuit-test patch supports M0, M1 and M2 only")
    bonds = []
    for i, j, theta in graph["directional_edges"]:
        phase = 2 * theta
        coupling = np.diag([jxy, jxy, jzz]).astype(float)
        coupling[0, 0] += 2 * jpp * np.cos(phase)
        coupling[1, 1] -= 2 * jpp * np.cos(phase)
        coupling[0, 1] += -2 * jpp * np.sin(phase)
        coupling[1, 0] += -2 * jpp * np.sin(phase)
        bonds.append(SpinBond(i, j, coupling * MEV_TO_K))
    return bonds


def spin_hamiltonian_to_sparse_pauli(n: int, bonds, longitudinal_field_K: float = 0.0):
    """Convert bonds and -h sum(Sz) to a simplified Qiskit SparsePauliOp."""
    coefficients = defaultdict(complex)

    def add_local(operators, coefficient):
        label = ["I"] * n
        for site, pauli in operators.items():
            label[n - 1 - site] = pauli
        coefficients["".join(label)] += coefficient

    # -h*Sz = +(h/2)*Z in the |down>,|up> convention.
    for site in range(n):
        add_local({site: "Z"}, longitudinal_field_K / 2)

    for bond in bonds:
        tensor = np.asarray(bond.coupling_K, dtype=complex)
        for axis_i in range(3):
            for axis_j in range(3):
                value = tensor[axis_i, axis_j]
                if abs(value) < 1e-15:
                    continue
                coefficient = value * SPIN_SIGNS[axis_i] * SPIN_SIGNS[axis_j] / 4
                add_local({bond.i: PAULI_BY_AXIS[axis_i],
                           bond.j: PAULI_BY_AXIS[axis_j]}, coefficient)

    operator = SparsePauliOp(list(coefficients), list(coefficients.values())).simplify(atol=1e-13)
    if np.max(np.abs(operator.coeffs.imag), initial=0.0) < 1e-13:
        operator = SparsePauliOp(operator.paulis, operator.coeffs.real)
    return operator

