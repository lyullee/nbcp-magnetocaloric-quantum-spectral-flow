"""Small-system validation of an MPO-squaring XTRG route using quimb.

This is not yet a bulk result.  It builds the same arbitrary Cartesian
spin-1/2 Hamiltonians used by the FTLM stage as an MPO, starts from the
high-temperature approximation ``I - beta_0 H``, and repeatedly squares and
compresses the normalized density operator.  Exact diagonalization is used
only for the N=9 validation reported by the command-line entry point.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np
import quimb
import quimb.operator as qop
import quimb.tensor as qtn

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
AXES = ("sx", "sy", "sz")


def build_spin_mpo(n, bonds, longitudinal_field_K=0.0):
    """Build an MPO for the Cartesian pair Hamiltonian in kelvin units."""
    hs = qop.HilbertSpace(n)
    builder = qop.SparseOperatorBuilder(hilbert_space=hs)
    for bond in bonds:
        coupling = np.asarray(bond.coupling_K)
        for a, opa in enumerate(AXES):
            for c, opc in enumerate(AXES):
                value = coupling[a, c]
                if abs(value) > 1e-15:
                    builder.add_term(value, (opa, bond.i), (opc, bond.j))
    for i in range(n):
        if abs(longitudinal_field_K) > 0:
            builder.add_term(-longitudinal_field_K, ("sz", i))
    return builder, builder.build_mpo(method="greedy")


def _real_scalar(value, name, imag_tolerance=1e-8):
    value = complex(value)
    if abs(value.imag) > imag_tolerance * max(1.0, abs(value.real)):
        raise FloatingPointError(f"{name} has a material imaginary part: {value}")
    return float(value.real)


def exact_thermodynamics(eigenvalues, beta, n):
    e0 = float(np.min(eigenvalues))
    weights = np.exp(-beta * (eigenvalues - e0))
    z_shift = float(weights.sum())
    energy = float(np.dot(weights, eigenvalues) / z_shift)
    log_z = -beta * e0 + math.log(z_shift)
    entropy = (log_z + beta * energy) / n
    return {"logZ": log_z, "energy_K": energy, "entropy_kB_per_spin": entropy}


def xtrg_square(
    h_mpo,
    n,
    beta_final,
    squarings=20,
    max_bond=64,
    cutoff=1e-12,
    eigenvalues=None,
    positive_square=True,
):
    """Cool an MPO density operator by repeated squaring and compression."""
    beta = beta_final / (2**squarings)
    identity = qtn.MPO_identity(n, dtype="complex128")
    rho = identity + h_mpo.multiply(-beta)
    z0 = _real_scalar(rho.trace(), "initial trace")
    if z0 <= 0:
        raise FloatingPointError(f"nonpositive initial trace {z0}")
    log_z = math.log(z0)
    rho.multiply_(1.0 / z0)
    rows = []

    def measure(step, elapsed):
        energy = _real_scalar(rho.apply(h_mpo, compress=False).trace(), "energy")
        entropy = (log_z + beta * energy) / n
        row = {
            "step": step,
            "beta_K_inv": beta,
            "temperature_K": 1.0 / beta,
            "logZ": log_z,
            "energy_K": energy,
            "entropy_kB_per_spin": entropy,
            "rho_trace": _real_scalar(rho.trace(), "normalized trace"),
            "rho_max_bond": int(rho.max_bond()),
            "elapsed_s": elapsed,
        }
        if eigenvalues is not None:
            exact = exact_thermodynamics(eigenvalues, beta, n)
            row["exact"] = exact
            row["errors"] = {
                "energy_K": energy - exact["energy_K"],
                "entropy_kB_per_spin": entropy - exact["entropy_kB_per_spin"],
                "logZ": log_z - exact["logZ"],
            }
        rows.append(row)

    measure(0, 0.0)
    for step in range(1, squarings + 1):
        start = time.perf_counter()
        # Compression can introduce a small anti-Hermitian component.  Using
        # rho^dagger rho instead of rho rho removes that component before each
        # cooling step and is identical to squaring when rho is exact.
        if positive_square:
            # In quimb ``.H`` conjugates the tensors; selecting the upper
            # legs of both operators supplies the transpose needed for the
            # dense-matrix equivalent rho^dagger @ rho.
            product = rho.H.apply(
                rho,
                which_A="upper",
                which_B="upper",
                compress=True,
                max_bond=max_bond,
                cutoff=cutoff,
            )
        else:
            product = rho.apply(
                rho,
                compress=True,
                max_bond=max_bond,
                cutoff=cutoff,
            )
        ratio = _real_scalar(product.trace(), "partition ratio")
        if ratio <= 0:
            raise FloatingPointError(f"nonpositive partition ratio {ratio}")
        log_z = 2.0 * log_z + math.log(ratio)
        product.multiply_(1.0 / ratio)
        rho = product
        beta *= 2.0
        measure(step, time.perf_counter() - start)

    dense_rho = np.asarray(rho.to_dense()) if n <= 10 else None
    diagnostic = {}
    if dense_rho is not None:
        hermitian_error = float(np.max(np.abs(dense_rho - dense_rho.conj().T)))
        evals = np.linalg.eigvalsh((dense_rho + dense_rho.conj().T) / 2)
        diagnostic = {
            "final_hermitian_max_abs_error": hermitian_error,
            "final_min_density_eigenvalue": float(evals.min()),
            "final_max_density_eigenvalue": float(evals.max()),
        }
    return rows, diagnostic


def run_validation(field_T=1.50, max_bonds=(16, 32, 64), squarings=20):
    graph = lattice((3, 0), (0, 3))
    bonds = model_bonds(graph, "M0")
    field_k = GZ * MU_B_OVER_KB * field_T
    builder, h_mpo = build_spin_mpo(graph["N"], bonds, field_k)
    dense_h = np.asarray(builder.build_dense())
    mpo_dense_error = float(np.max(np.abs(np.asarray(h_mpo.to_dense()) - dense_h)))
    eigenvalues = np.linalg.eigvalsh(dense_h)
    beta_final = 1.0 / 0.09

    out = ROOT / "results" / "xtrg_quimb_n9_validation_positive_square.json"
    runs = []
    for max_bond in max_bonds:
        start = time.perf_counter()
        rows, diagnostic = xtrg_square(
            h_mpo,
            graph["N"],
            beta_final,
            squarings=squarings,
            max_bond=max_bond,
            cutoff=1e-12,
            eigenvalues=eigenvalues,
            positive_square=True,
        )
        runs.append(
            {
                "max_bond": max_bond,
                "wall_s": time.perf_counter() - start,
                "final": rows[-1],
                "diagnostic": diagnostic,
                "trajectory": rows,
            }
        )
        # Checkpoint every bond dimension so a costly later run cannot erase
        # already completed validation data.
        partial = {
            "purpose": "N=9 exact validation of local MPO-squaring XTRG route; not a bulk claim",
            "model": "M0 current 2026 parameters",
            "N": graph["N"],
            "field_T": field_T,
            "target_temperature_K": 0.09,
            "squarings": squarings,
            "initial_beta_K_inv": beta_final / (2**squarings),
            "initialization": "first-order high-temperature MPO I-beta0*H",
            "squaring": "rho_dagger_times_rho before compression",
            "hamiltonian_initial_max_bond": int(h_mpo.max_bond()),
            "hamiltonian_MPO_dense_max_abs_error": mpo_dense_error,
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "quimb": quimb.__version__,
            },
            "runs": runs,
        }
        out.write_text(json.dumps(partial, indent=2), encoding="utf-8")
        final = runs[-1]["final"]
        print(
            f"checkpoint D={max_bond:>3} wall={runs[-1]['wall_s']:.2f}s "
            f"S={final['entropy_kB_per_spin']:.10f} "
            f"dS={final['errors']['entropy_kB_per_spin']:+.3e} "
            f"dE={final['errors']['energy_K']:+.3e}",
            flush=True,
        )

    payload = {
        "purpose": "N=9 exact validation of local MPO-squaring XTRG route; not a bulk claim",
        "model": "M0 current 2026 parameters",
        "N": graph["N"],
        "field_T": field_T,
        "target_temperature_K": 0.09,
        "squarings": squarings,
        "initial_beta_K_inv": beta_final / (2**squarings),
        "initialization": "first-order high-temperature MPO I-beta0*H",
        "squaring": "rho_dagger_times_rho before compression",
        "hamiltonian_initial_max_bond": int(h_mpo.max_bond()),
        "hamiltonian_MPO_dense_max_abs_error": mpo_dense_error,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "quimb": quimb.__version__,
        },
        "runs": runs,
    }
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload, out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--field", type=float, default=1.50)
    parser.add_argument("--max-bonds", type=int, nargs="+", default=[16, 32, 64])
    parser.add_argument("--squarings", type=int, default=20)
    args = parser.parse_args()
    payload, out = run_validation(args.field, tuple(args.max_bonds), args.squarings)
    print(out)
    for run in payload["runs"]:
        final = run["final"]
        print(
            f"D={run['max_bond']:>3} wall={run['wall_s']:.2f}s "
            f"S={final['entropy_kB_per_spin']:.10f} "
            f"dS={final['errors']['entropy_kB_per_spin']:+.3e} "
            f"dE={final['errors']['energy_K']:+.3e}"
        )


if __name__ == "__main__":
    main()
