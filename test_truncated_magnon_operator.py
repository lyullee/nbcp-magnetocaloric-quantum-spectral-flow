import numpy as np

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ
from truncated_magnon_operator import FixedMagnonCartesianOperator, TruncatedMagnonCartesianOperator


def main():
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    bonds = model_bonds(graph, "M3")
    field = GZ * MU_B_OVER_KB * 1.5
    for parity in (0, 1):
        full = SpinHamiltonianOperator(n, bonds, field, parity=(n - parity) % 2)
        projected = TruncatedMagnonCartesianOperator(n, n, parity, bonds, field)
        rng = np.random.default_rng(100 + parity)
        x = rng.normal(size=projected.shape[0]) + 1j * rng.normal(size=projected.shape[0])
        # Existing full-space operator uses bit=1 for spin up; the projected
        # near-saturation operator uses bit=1 for a down-spin magnon.
        complemented = projected.states ^ np.uint64((1 << n) - 1)
        full_indices = full.state_index[complemented]
        x_full = np.empty_like(x)
        x_full[full_indices] = x
        y_full = (full @ x_full)[full_indices]
        assert np.max(np.abs(y_full - projected @ x)) < 1e-11
        for q in range(parity, n + 1, 2):
            fixed = FixedMagnonCartesianOperator(n, q, bonds, field)
            start = projected.offsets[q]
            stop = start + fixed.shape[0]
            z = rng.normal(size=fixed.shape[0]) + 1j * rng.normal(size=fixed.shape[0])
            embedded = np.zeros(projected.shape[0], dtype=complex)
            embedded[start:stop] = z
            assert np.max(np.abs((projected @ embedded)[start:stop] - fixed @ z)) < 1e-11
    print("truncated-magnon Cartesian operator validation passed")


if __name__ == "__main__":
    main()
