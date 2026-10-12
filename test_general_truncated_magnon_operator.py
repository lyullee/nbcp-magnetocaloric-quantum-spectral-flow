"""Exact small-space validation of the parity-breaking magnon projection."""
import numpy as np

from general_truncated_magnon_operator import GeneralTruncatedMagnonCartesianOperator
from matrix_free_hamiltonian import SpinBond, SpinHamiltonianOperator


def main():
    n = 6
    rng = np.random.default_rng(20261012)
    bonds = []
    for i in range(n):
        j = (i + 1) % n
        tensor = rng.normal(scale=0.07, size=(3, 3))
        tensor = (tensor + tensor.T) / 2
        # Keep every Cartesian component, especially the parity-breaking ones.
        bonds.append(SpinBond(i, j, tensor))
    field = 0.83
    full = SpinHamiltonianOperator(n, bonds, field, parity=None)
    projected = GeneralTruncatedMagnonCartesianOperator(n, n, bonds, field)
    assert projected.shape == full.shape

    x = rng.normal(size=projected.shape[0]) + 1j * rng.normal(size=projected.shape[0])
    complemented = projected.states ^ np.uint64((1 << n) - 1)
    x_full = np.empty_like(x)
    x_full[complemented] = x
    y_full = (full @ x_full)[complemented]
    error = float(np.max(np.abs(y_full - projected @ x)))
    if error >= 1e-11:
        raise AssertionError(f"general projected operator mismatch: {error}")

    # A cutoff projection must equal the retained principal block of the exact
    # full operator; transitions out of the retained space are discarded.
    cutoff = 3
    truncated = GeneralTruncatedMagnonCartesianOperator(n, cutoff, bonds, field)
    kept = truncated.states
    z = rng.normal(size=kept.size) + 1j * rng.normal(size=kept.size)
    embedded = np.zeros(projected.shape[0], dtype=complex)
    embedded[:kept.size] = z
    block_error = float(np.max(np.abs((projected @ embedded)[:kept.size] - truncated @ z)))
    if block_error >= 1e-11:
        raise AssertionError(f"cutoff principal-block mismatch: {block_error}")
    print({"full_error": error, "cutoff_block_error": block_error, "status": "PASS"})


if __name__ == "__main__":
    main()
