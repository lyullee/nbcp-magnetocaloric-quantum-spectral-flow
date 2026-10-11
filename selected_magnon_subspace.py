"""Selected-M0-eigenstate subspace approximation for low-energy M3 spectra.

The approximation is validated against direct M3 eigenspectra before it may be
used at N=27. Basis states can come from M0 or from the q-conserving block of
M3 in every fixed-q sector q<=7. M3 is projected exactly into their span.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

from magnon_sector import FixedMagnonXXZOperator, periodic_graph_edges
from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import DELTA_XY_OVER_ZZ, GZ, JZZ_MEV, MEV_TO_K
from truncated_magnon_operator import FixedMagnonCartesianOperator, TruncatedMagnonCartesianOperator

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CELLS = {
    18: ((4, 1), (-2, 4)),
    21: ((-5, 1), (-1, -4)),
    24: ((-5, 1), (1, -5)),
    27: ((-6, 3), (-3, -3)),
}
TEMPERATURES = (0.09, 0.10, 0.12)


def low_eigenpairs(op, count, seed):
    n = op.shape[0]
    if n <= count:
        eye = np.eye(n)
        dense = np.column_stack([op @ eye[:, j] for j in range(n)])
        values, vectors = np.linalg.eigh(dense)
        return values, vectors
    rng = np.random.default_rng(seed)
    v0 = rng.normal(size=n)
    if np.issubdtype(op.dtype, np.complexfloating):
        v0 = v0 + 1j * rng.normal(size=n)
    request = min(count, n - 1)
    values, vectors = eigsh(op, k=request, which="SA", tol=1e-9, v0=v0)
    # Complex ARPACK can return non-orthogonal vectors inside exact or near
    # degeneracies.  Orthonormalize the returned span, then solve the small
    # Rayleigh--Ritz problem again before using it as a projection basis.
    vectors, _ = np.linalg.qr(vectors)
    hv = np.column_stack([op @ vectors[:, j] for j in range(vectors.shape[1])])
    rayleigh = vectors.conj().T @ hv
    rayleigh = (rayleigh + rayleigh.conj().T) / 2
    values, rotation = eigh(rayleigh)
    vectors = vectors @ rotation
    order = np.argsort(values.real)
    values, vectors = values.real[order], vectors[:, order]
    gram_error = np.max(np.abs(vectors.conj().T @ vectors - np.eye(vectors.shape[1])))
    if gram_error > 1e-10:
        raise RuntimeError(f"fixed-q Ritz vectors are not orthonormal: {gram_error}")
    return values, vectors


def entropy(energies, temperature, n):
    e = np.sort(np.asarray(energies, float))
    shifted = e - e[0]
    w = np.exp(-shifted / temperature)
    return float((np.log(w.sum()) + np.dot(w, shifted) / w.sum() / temperature) / n)


def reference_low16(n, model):
    if n == 18:
        data = json.loads((RESULTS / "n18_bulk_window_compact_raw" / f"{model}_B1.500T.json").read_text())
        return sorted(np.concatenate([np.asarray(x)[:8] for x in data["low_energies_K_by_sector"]]).tolist())
    data = json.loads((RESULTS / f"magnon_M3_low_energy_N{n}_B1.500T.json").read_text())
    return data["runs"][-1]["models"][model]["combined_low_energies_K"]


def build_parity_projection(n, graph, m3_bonds, field_k, parity, max_basis, basis_source, started):
    cutoff = 6 if parity == 0 else 7
    m3 = TruncatedMagnonCartesianOperator(n, cutoff, parity, m3_bonds, field_k)
    blocks = []
    column_count = 0
    for q in range(parity, cutoff + 1, 2):
        if basis_source == "M0":
            fixed = FixedMagnonXXZOperator(
                n, q, periodic_graph_edges(graph),
                JZZ_MEV * DELTA_XY_OVER_ZZ * MEV_TO_K,
                JZZ_MEV * MEV_TO_K, field_k,
            )
        elif basis_source == "M3q":
            fixed = FixedMagnonCartesianOperator(n, q, m3_bonds, field_k)
        else:
            raise ValueError("basis_source must be M0 or M3q")
        values, vectors = low_eigenpairs(fixed, max_basis, 20261011 + 100 * parity + q)
        start = int(m3.offsets[q]); stop = start + fixed.shape[0]
        fixed_basis = fixed.basis if hasattr(fixed, "basis") else fixed.states
        if not np.array_equal(m3.states[start:stop], fixed_basis):
            raise RuntimeError(f"basis mismatch at q={q}")
        blocks.append({
            "q": q, "start": start, "stop": stop, "dimension": fixed.shape[0],
            "values": values, "vectors": vectors.astype(np.complex128, copy=False),
            "column_start": column_count, "column_stop": column_count + vectors.shape[1],
        })
        column_count += vectors.shape[1]
        print(
            f"N={n} source={basis_source} parity={parity} q={q} dim={fixed.shape[0]} m={vectors.shape[1]} "
            f"E0={values[0]:.9f} elapsed={time.perf_counter()-started:.1f}s",
            flush=True,
        )

    hproj = np.zeros((column_count, column_count), dtype=np.complex128)
    x = np.zeros(m3.shape[0], dtype=np.complex128)
    for source in blocks:
        for local in range(source["vectors"].shape[1]):
            x[source["start"]:source["stop"]] = source["vectors"][:, local]
            y = m3 @ x
            x[source["start"]:source["stop"]] = 0
            col = source["column_start"] + local
            for target in blocks:
                hproj[target["column_start"]:target["column_stop"], col] = (
                    target["vectors"].conj().T @ y[target["start"]:target["stop"]]
                )
        print(
            f"N={n} parity={parity} projected source q={source['q']} "
            f"elapsed={time.perf_counter()-started:.1f}s",
            flush=True,
        )
    hermiticity = float(np.max(np.abs(hproj - hproj.conj().T)))
    hproj = (hproj + hproj.conj().T) / 2
    return m3, blocks, hproj, hermiticity


def selected_indices(blocks, basis_per_q):
    out = []
    for block in blocks:
        out.extend(range(block["column_start"], min(block["column_start"] + basis_per_q, block["column_stop"])))
    return np.asarray(out, dtype=int)


def approximate_vector(m3, blocks, selected, coefficients):
    psi = np.zeros(m3.shape[0], dtype=np.complex128)
    position = {int(full): j for j, full in enumerate(selected)}
    for block in blocks:
        columns = [c for c in range(block["column_start"], block["column_stop"]) if c in position]
        if not columns:
            continue
        local = np.asarray([c - block["column_start"] for c in columns])
        coeff = np.asarray([coefficients[position[c]] for c in columns])
        psi[block["start"]:block["stop"]] = block["vectors"][:, local] @ coeff
    return psi


def main(n, max_basis, basis_source):
    if n not in CELLS:
        raise ValueError(f"supported N: {sorted(CELLS)}")
    a, b = CELLS[n]
    graph = lattice(a, b)
    field_t = 1.50
    field_k = GZ * MU_B_OVER_KB * field_t
    bonds = model_bonds(graph, "M3")
    started = time.perf_counter()
    parity_data = []
    for parity in (0, 1):
        m3, blocks, hproj, hermiticity = build_parity_projection(
            n, graph, bonds, field_k, parity, max_basis, basis_source, started
        )
        parity_data.append((m3, blocks, hproj, hermiticity))

    reference = {model: reference_low16(n, model) for model in ("M0", "M3")} if n < 27 else None
    runs = []
    basis_values = sorted(set([4, 8, max_basis]))
    for basis_per_q in basis_values:
        combined = []
        parity_records = []
        for parity, (m3, blocks, hproj, hermiticity) in enumerate(parity_data):
            selected = selected_indices(blocks, basis_per_q)
            values, vectors = eigh(hproj[np.ix_(selected, selected)])
            take = min(8, values.size)
            low_values = values[:take].real
            residuals = []
            if basis_per_q == max_basis:
                for j, energy in enumerate(low_values):
                    psi = approximate_vector(m3, blocks, selected, vectors[:, j])
                    residuals.append(float(np.linalg.norm(m3 @ psi - energy * psi)))
            combined.extend(low_values.tolist())
            parity_records.append({
                "parity": parity,
                "selected_dimension": int(selected.size),
                "projected_hermiticity_error_K": hermiticity,
                "low_energies_K": low_values.tolist(),
                "full_space_residual_norms_K": residuals,
            })
        rows = []
        for t in TEMPERATURES:
            s_m3 = entropy(combined, t, n)
            row = {"temperature_K": t, "S_M3_selected": s_m3}
            if reference is not None:
                s_ref = entropy(reference["M3"], t, n)
                s0_ref = entropy(reference["M0"], t, n)
                row.update({
                    "S_M3_reference_low16": s_ref,
                    "M3_entropy_error_kB_per_spin": s_m3 - s_ref,
                    "deltaS_selected_vs_reference_M0": s_m3 - s0_ref,
                })
            rows.append(row)
        runs.append({"basis_states_per_q": basis_per_q, "parities": parity_records, "rows": rows})

    payload = {
        "method": "M3 projected into low fixed-q eigenstates",
        "condition": {"N": n, "translations": [a, b], "field_T": field_t, "max_q_even": 6, "max_q_odd": 7},
        "settings": {"basis_source": basis_source, "maximum_basis_states_per_q": max_basis, "eigsh_tolerance": 1e-9},
        "reference": None if reference is None else "direct M3 lowest eight states per parity",
        "runs": runs,
        "elapsed_seconds": time.perf_counter() - started,
    }
    tag = "M0basis" if basis_source == "M0" else "M3qbasis"
    out = RESULTS / f"selected_magnon_{tag}_N{n}_B1.500T.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": out.name,
        "rows": [{"m": run["basis_states_per_q"], "rows": run["rows"]} for run in runs],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--max-basis", type=int, default=12)
    parser.add_argument("--basis-source", choices=("M0", "M3q"), default="M3q")
    args = parser.parse_args()
    main(args.n, args.max_basis, args.basis_source)
