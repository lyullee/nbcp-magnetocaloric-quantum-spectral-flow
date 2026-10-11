from itertools import combinations

import numpy as np

from magnon_sector import FixedMagnonXXZOperator, cylinder_edges


def dense_reference(n, q, edges, jxy, jz, h):
    basis = np.array([sum(1 << i for i in c) for c in combinations(range(n), q)], dtype=np.uint64)
    basis.sort()
    index = {int(state): p for p, state in enumerate(basis)}
    out = np.zeros((basis.size, basis.size))
    for p, state_u in enumerate(basis):
        state = int(state_u)
        out[p, p] -= h * (n / 2 - q)
        for i, j in edges:
            unequal = ((state >> i) & 1) != ((state >> j) & 1)
            out[p, p] += jz * (-0.25 if unequal else 0.25)
            if unequal:
                out[index[state ^ (1 << i) ^ (1 << j)], p] += jxy / 2
    return out


def main():
    n, q = 9, 3
    edges = cylinder_edges(3, 3)
    op = FixedMagnonXXZOperator(n, q, edges, .88, 1.48, 4.89 * .67171381563 * 1.5)
    dense = dense_reference(n, q, edges, .88, 1.48, 4.89 * .67171381563 * 1.5)
    assembled = np.column_stack([op @ np.eye(op.shape[0])[:, i] for i in range(op.shape[0])])
    assert np.max(np.abs(assembled - dense)) < 1e-12
    assert np.max(np.abs(assembled - assembled.T)) < 1e-12
    print("fixed-magnon XXZ operator validation passed")


if __name__ == "__main__":
    main()
