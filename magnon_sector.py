"""Fixed-magnon sectors for the near-saturation triangular XXZ model.

Bits mark spins flipped down from the fully polarized state.  Basis masks are
ordered by their integer value, whose fixed-popcount rank is the combinadic
colex rank sum_j C(position_j, j).
"""
from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np
from numba import njit
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import LinearOperator, eigsh


@njit(cache=True)
def _fill_fixed_popcount_basis(n: int, q: int, size: int) -> np.ndarray:
    out = np.empty(size, dtype=np.uint64)
    if q == 0:
        out[0] = 0
        return out
    x = (np.uint64(1) << np.uint64(q)) - np.uint64(1)
    limit = np.uint64(1) << np.uint64(n)
    p = 0
    while x < limit:
        out[p] = x
        p += 1
        c = x & (-x)
        r = x + c
        x = (((r ^ x) >> np.uint64(2)) // c) | r
    return out


@njit(cache=True)
def _colex_rank(mask: np.uint64, binom: np.ndarray) -> int:
    rank = np.int64(0)
    order = 1
    pos = 0
    x = mask
    while x:
        if x & np.uint64(1):
            rank += binom[pos, order]
            order += 1
        pos += 1
        x >>= np.uint64(1)
    return int(rank)


@njit(cache=True)
def _diagonal_values(basis: np.ndarray, edges: np.ndarray, n: int, q: int,
                     jz: float, field_energy: float) -> np.ndarray:
    out = np.empty(basis.size, dtype=np.float64)
    magnetization = 0.5 * n - q
    for p in range(basis.size):
        s = basis[p]
        unequal = 0
        for e in range(edges.shape[0]):
            i, j = edges[e, 0], edges[e, 1]
            unequal += int(((s >> np.uint64(i)) ^ (s >> np.uint64(j))) & np.uint64(1))
        out[p] = jz * (0.25 * edges.shape[0] - 0.5 * unequal) - field_energy * magnetization
    return out


@njit(cache=True)
def _sector_matvec(basis: np.ndarray, diagonal: np.ndarray, neighbors: np.ndarray,
                   binom: np.ndarray, q: int, hopping: float, x: np.ndarray) -> np.ndarray:
    y = diagonal * x
    for p in range(basis.size):
        s = basis[p]
        remaining = s
        while remaining:
            low = remaining & (-remaining)
            site = 0
            tmp = low
            while tmp > np.uint64(1):
                tmp >>= np.uint64(1)
                site += 1
            for a in range(neighbors.shape[1]):
                other = neighbors[site, a]
                if other < 0:
                    break
                if ((s >> np.uint64(other)) & np.uint64(1)) == 0:
                    target = s ^ low ^ (np.uint64(1) << np.uint64(other))
                    target_index = _colex_rank(target, binom)
                    y[target_index] += hopping * x[p]
            remaining ^= low
    return y


def cylinder_edges(lx: int, ly: int) -> list[tuple[int, int]]:
    """YC-like triangular cylinder matching QMagen OBCx/PBCy ordering."""
    edges: set[tuple[int, int]] = set()
    site = lambda x, y: x * ly + (y % ly)
    for x in range(lx):
        for y in range(ly):
            edges.add(tuple(sorted((site(x, y), site(x, y + 1)))))
            if x + 1 < lx:
                edges.add(tuple(sorted((site(x, y), site(x + 1, y)))))
                edges.add(tuple(sorted((site(x, y), site(x + 1, y + 1)))))
    return sorted(edges)


def periodic_graph_edges(graph: dict) -> list[tuple[int, int]]:
    return [tuple(map(int, edge)) for edge in graph["edges"]]


class FixedMagnonXXZOperator(LinearOperator):
    def __init__(self, n: int, q: int, edges, jxy_k: float, jz_k: float,
                 field_energy_k: float):
        self.n, self.q = int(n), int(q)
        self.edges = np.asarray(edges, dtype=np.int16)
        self.basis = _fill_fixed_popcount_basis(self.n, self.q, comb(self.n, self.q))
        self.binom = np.zeros((self.n + 1, self.q + 2), dtype=np.int64)
        for p in range(self.n + 1):
            for r in range(min(p, self.q + 1) + 1):
                self.binom[p, r] = comb(p, r)
        adjacency = [[] for _ in range(self.n)]
        for i, j in self.edges:
            adjacency[int(i)].append(int(j)); adjacency[int(j)].append(int(i))
        width = max(map(len, adjacency))
        self.neighbors = np.full((self.n, width), -1, dtype=np.int16)
        for i, row in enumerate(adjacency):
            self.neighbors[i, :len(row)] = sorted(row)
        self.diagonal = _diagonal_values(
            self.basis, self.edges, self.n, self.q, float(jz_k), float(field_energy_k)
        )
        self.hopping = 0.5 * float(jxy_k)
        super().__init__(dtype=np.dtype(np.float64), shape=(self.basis.size, self.basis.size))

    def _matvec(self, x):
        return _sector_matvec(
            self.basis, self.diagonal, self.neighbors, self.binom,
            self.q, self.hopping, np.asarray(x, dtype=np.float64)
        )


@dataclass
class SLQRule:
    dimension: int
    energies: np.ndarray
    weights: np.ndarray
    last_beta: float


def stochastic_lanczos_rule(operator: LinearOperator, steps: int, rng: np.random.Generator,
                            deflation_vectors: np.ndarray | None = None,
                            trace_dimension: int | None = None) -> SLQRule:
    """Three-vector stochastic Lanczos quadrature rule for a real operator."""
    n = operator.shape[0]
    q = rng.choice(np.array([-1.0, 1.0]), size=n)
    if deflation_vectors is not None and deflation_vectors.size:
        q -= deflation_vectors @ (deflation_vectors.T @ q)
    q /= np.linalg.norm(q)
    q_prev = np.zeros(n)
    beta_prev = 0.0
    alpha = []
    beta = []
    for step in range(min(steps, n)):
        z = operator @ q
        if step:
            z -= beta_prev * q_prev
        a = float(np.dot(q, z))
        z -= a * q
        if deflation_vectors is not None and deflation_vectors.size:
            z -= deflation_vectors @ (deflation_vectors.T @ z)
        # One local correction limits roundoff without storing the Krylov basis.
        z -= float(np.dot(q, z)) * q
        b = float(np.linalg.norm(z))
        alpha.append(a)
        if step + 1 == min(steps, n) or b < 1e-13:
            break
        beta.append(b)
        q_prev, q = q, z / b
        beta_prev = b
    values, vectors = eigh_tridiagonal(np.asarray(alpha), np.asarray(beta))
    return SLQRule(n if trace_dimension is None else trace_dimension,
                   values, vectors[0] ** 2, beta[-1] if beta else 0.0)


@dataclass
class DeflatedSLQSector:
    dimension: int
    low_energies: np.ndarray
    random_rules: list[SLQRule]


def deflated_sector_rules(operator: LinearOperator, low_states: int, random_vectors: int,
                          steps: int, rng: np.random.Generator) -> DeflatedSLQSector:
    n = operator.shape[0]
    if n <= max(32, low_states + 1):
        eye = np.eye(n)
        dense = np.column_stack([operator @ eye[:, i] for i in range(n)])
        return DeflatedSLQSector(n, np.linalg.eigvalsh(dense), [])
    count = min(low_states, n - 1)
    energies, vectors = eigsh(operator, k=count, which="SA", tol=1e-10)
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    complement = n - count
    rules = [stochastic_lanczos_rule(operator, steps, rng, vectors, complement)
             for _ in range(random_vectors)]
    return DeflatedSLQSector(n, energies, rules)


def deflated_slq_thermodynamics(sectors: list[DeflatedSLQSector], temperatures, spins: int):
    reference = min(float(sector.low_energies.min()) for sector in sectors)
    rows = []
    for temperature in temperatures:
        z = e1 = e2 = 0.0
        for sector in sectors:
            exact_w = np.exp(-(sector.low_energies - reference) / temperature)
            z += float(exact_w.sum())
            e1 += float(np.dot(exact_w, sector.low_energies))
            e2 += float(np.dot(exact_w, sector.low_energies ** 2))
            if sector.random_rules:
                scale = (sector.dimension - len(sector.low_energies)) / len(sector.random_rules)
                for rule in sector.random_rules:
                    w = rule.weights * np.exp(-(rule.energies - reference) / temperature)
                    z += scale * float(w.sum())
                    e1 += scale * float(np.dot(w, rule.energies))
                    e2 += scale * float(np.dot(w, rule.energies ** 2))
        energy = e1 / z
        entropy_total = np.log(z) + (energy - reference) / temperature
        rows.append({
            "temperature_K": float(temperature),
            "energy_K_total": energy,
            "entropy_kB_total": entropy_total,
            "entropy_kB_per_spin": entropy_total / spins,
            "heat_capacity_kB_total": (e2 / z - energy ** 2) / temperature ** 2,
            "shift_energy_K": reference,
        })
    return rows


def slq_thermodynamics(rules_by_sector: list[list[SLQRule]], temperatures, spins: int):
    all_min = min(float(rule.energies.min()) for rules in rules_by_sector for rule in rules)
    rows = []
    for temperature in temperatures:
        z = e1 = e2 = 0.0
        for rules in rules_by_sector:
            dim = rules[0].dimension
            scale = dim / len(rules)
            for rule in rules:
                w = rule.weights * np.exp(-(rule.energies - all_min) / temperature)
                z += scale * float(w.sum())
                e1 += scale * float(np.dot(w, rule.energies))
                e2 += scale * float(np.dot(w, rule.energies ** 2))
        energy = e1 / z
        entropy_total = np.log(z) + (energy - all_min) / temperature
        rows.append({
            "temperature_K": float(temperature),
            "energy_K_total": energy,
            "entropy_kB_total": entropy_total,
            "entropy_kB_per_spin": entropy_total / spins,
            "heat_capacity_kB_total": (e2 / z - energy ** 2) / temperature ** 2,
            "shift_energy_K": all_min,
        })
    return rows
