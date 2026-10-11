"""Projected low-magnon Hamiltonian for parity-preserving Cartesian bonds."""
from __future__ import annotations

from math import comb

import numpy as np
from numba import njit
from scipy.sparse.linalg import LinearOperator

from magnon_sector import _colex_rank, _fill_fixed_popcount_basis


@njit(cache=True)
def _projected_diagonal(states, pair_i, pair_j, pair_zz, n, field_energy):
    out = np.empty(states.size, dtype=np.float64)
    for p in range(states.size):
        s = states[p]
        q = 0
        x = s
        while x:
            q += int(x & np.uint64(1)); x >>= np.uint64(1)
        value = -field_energy * (0.5 * n - q)
        for b in range(pair_i.size):
            bi = int((s >> np.uint64(pair_i[b])) & np.uint64(1))
            bj = int((s >> np.uint64(pair_j[b])) & np.uint64(1))
            value += pair_zz[b] * (0.5 - bi) * (0.5 - bj)
        out[p] = value
    return out


@njit(cache=True)
def _projected_matvec(states, diagonal, pair_i, pair_j, lookup, offsets, binom, cutoff, x):
    y = diagonal.astype(np.complex128) * x
    for p in range(states.size):
        s = states[p]
        q_state = 0
        tmp = s
        while tmp:
            q_state += int(tmp & np.uint64(1)); tmp >>= np.uint64(1)
        for b in range(pair_i.size):
            i, j = pair_i[b], pair_j[b]
            bi = int((s >> np.uint64(i)) & np.uint64(1))
            bj = int((s >> np.uint64(j)) & np.uint64(1))
            coefficient = lookup[b, bi + 2 * bj]
            if coefficient == 0:
                continue
            q_change = 2 - 2 * (bi + bj)
            target_q = q_state + q_change
            if target_q < 0 or target_q > cutoff or offsets[target_q] < 0:
                continue
            target = s ^ (np.uint64(1) << np.uint64(i)) ^ (np.uint64(1) << np.uint64(j))
            target_index = offsets[target_q] + _colex_rank(target, binom)
            y[target_index] += coefficient * x[p]
    return y


@njit(cache=True)
def _fixed_projected_matvec(states, diagonal, pair_i, pair_j, lookup, binom, x):
    """Cartesian-bond matvec projected to one fixed-q block."""
    y = diagonal.astype(np.complex128) * x
    for p in range(states.size):
        s = states[p]
        for b in range(pair_i.size):
            i, j = pair_i[b], pair_j[b]
            bi = int((s >> np.uint64(i)) & np.uint64(1))
            bj = int((s >> np.uint64(j)) & np.uint64(1))
            if bi + bj != 1:
                continue
            coefficient = lookup[b, bi + 2 * bj]
            if coefficient == 0:
                continue
            target = s ^ (np.uint64(1) << np.uint64(i)) ^ (np.uint64(1) << np.uint64(j))
            target_index = _colex_rank(target, binom)
            y[target_index] += coefficient * x[p]
    return y


class TruncatedMagnonCartesianOperator(LinearOperator):
    def __init__(self, n: int, cutoff: int, parity: int, bonds, field_energy_k: float):
        self.n, self.cutoff, self.parity = int(n), int(cutoff), int(parity)
        if parity not in (0, 1):
            raise ValueError("parity must be 0 or 1")
        combined = {}
        for bond in bonds:
            i, j = int(bond.i), int(bond.j)
            tensor = np.asarray(bond.coupling_K, dtype=np.complex128)
            if i > j:
                i, j = j, i; tensor = tensor.T
            combined[(i, j)] = combined.get((i, j), np.zeros((3, 3), complex)) + tensor
        items = sorted(combined.items())
        self.pair_i = np.array([key[0] for key, _ in items], dtype=np.int16)
        self.pair_j = np.array([key[1] for key, _ in items], dtype=np.int16)
        self.pair_zz = np.array([tensor[2, 2].real for _, tensor in items], dtype=np.float64)
        lookup = np.zeros((len(items), 4), dtype=np.complex128)
        for b, (_, tensor) in enumerate(items):
            if np.max(np.abs(tensor[[0, 1], 2])) > 1e-12 or np.max(np.abs(tensor[2, [0, 1]])) > 1e-12:
                raise ValueError("single-spin-flip xz/yz terms break parity")
            for local in range(4):
                bi, bj = local & 1, (local >> 1) & 1
                cx_i, cx_j = .5, .5
                cy_i, cy_j = .5j * (1 - 2 * bi), .5j * (1 - 2 * bj)
                lookup[b, local] = (
                    tensor[0, 0] * cx_i * cx_j + tensor[0, 1] * cx_i * cy_j
                    + tensor[1, 0] * cy_i * cx_j + tensor[1, 1] * cy_i * cy_j
                )
        self.lookup = lookup
        self.offsets = np.full(cutoff + 1, -1, dtype=np.int64)
        blocks = []
        offset = 0
        for q in range(parity, cutoff + 1, 2):
            self.offsets[q] = offset
            block = _fill_fixed_popcount_basis(n, q, comb(n, q))
            blocks.append(block); offset += block.size
        self.states = np.concatenate(blocks)
        self.binom = np.zeros((n + 1, cutoff + 2), dtype=np.int64)
        for p in range(n + 1):
            for r in range(min(p, cutoff + 1) + 1):
                self.binom[p, r] = comb(p, r)
        self.diagonal = _projected_diagonal(
            self.states, self.pair_i, self.pair_j, self.pair_zz, n, float(field_energy_k)
        )
        super().__init__(dtype=np.dtype(np.complex128), shape=(self.states.size, self.states.size))

    def _matvec(self, x):
        return _projected_matvec(
            self.states, self.diagonal, self.pair_i, self.pair_j, self.lookup,
            self.offsets, self.binom, self.cutoff, np.asarray(x, dtype=np.complex128)
        )

    def _rmatvec(self, x):
        return self._matvec(x)


class FixedMagnonCartesianOperator(LinearOperator):
    """The q-conserving diagonal block of the general Cartesian Hamiltonian."""

    def __init__(self, n: int, q: int, bonds, field_energy_k: float):
        self.n, self.q = int(n), int(q)
        helper = TruncatedMagnonCartesianOperator(n, q, q % 2, bonds, field_energy_k)
        start = int(helper.offsets[q])
        self.states = helper.states[start:].copy()
        self.diagonal = helper.diagonal[start:].copy()
        self.pair_i = helper.pair_i
        self.pair_j = helper.pair_j
        self.lookup = helper.lookup
        self.binom = helper.binom
        super().__init__(dtype=np.dtype(np.complex128), shape=(self.states.size, self.states.size))

    def _matvec(self, x):
        return _fixed_projected_matvec(
            self.states, self.diagonal, self.pair_i, self.pair_j,
            self.lookup, self.binom, np.asarray(x, dtype=np.complex128)
        )

    def _rmatvec(self, x):
        return self._matvec(x)
