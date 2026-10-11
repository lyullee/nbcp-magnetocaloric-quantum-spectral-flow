"""Positive MPO-factor thermal squaring, validated against N=9 exact data.

The factor X represents ``exp(-beta H / 2) / sqrt(Z)`` and the physical
density operator is always formed as ``rho = X^dagger X``.  Compression can
make X non-Hermitian, but cannot make rho non-positive.  This is a numerical
feasibility test, not a bulk-material result.
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
import quimb.tensor as qtn

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ
from xtrg_quimb import build_spin_mpo, exact_thermodynamics, _real_scalar


ROOT = Path(__file__).resolve().parent


def factor_norm_squared(x_mpo):
    rho = x_mpo.H.apply(
        x_mpo,
        which_A="upper",
        which_B="upper",
        compress=False,
    )
    return _real_scalar(rho.trace(), "factor norm")


def measure_factor(x_mpo, h_mpo, log_z, beta, n, eigenvalues=None):
    rho = x_mpo.H.apply(
        x_mpo,
        which_A="upper",
        which_B="upper",
        compress=False,
    )
    trace = _real_scalar(rho.trace(), "rho trace")
    energy = _real_scalar(
        rho.apply(h_mpo, compress=False).trace(),
        "energy",
        imag_tolerance=1e-5,
    ) / trace
    entropy = (log_z + beta * energy) / n
    row = {
        "beta_K_inv": beta,
        "temperature_K": 1.0 / beta,
        "logZ": log_z,
        "energy_K": energy,
        "entropy_kB_per_spin": entropy,
        "rho_trace": trace,
        "factor_max_bond": int(x_mpo.max_bond()),
    }
    if eigenvalues is not None:
        exact = exact_thermodynamics(eigenvalues, beta, n)
        row["exact"] = exact
        row["errors"] = {
            "energy_K": energy - exact["energy_K"],
            "entropy_kB_per_spin": entropy - exact["entropy_kB_per_spin"],
            "logZ": log_z - exact["logZ"],
        }
    return row, rho


def purified_square(
    h_mpo,
    n,
    beta_final,
    squarings=18,
    max_bond=24,
    cutoff=1e-12,
    eigenvalues=None,
    compression_method="direct",
    fit_max_iterations=4,
):
    beta = beta_final / (2**squarings)
    x_mpo = qtn.MPO_identity(n, dtype="complex128") + h_mpo.multiply(-beta / 2)
    z0 = factor_norm_squared(x_mpo)
    log_z = math.log(z0)
    x_mpo.multiply_(1 / math.sqrt(z0))
    rows = [{
        "step": 0,
        "beta_K_inv": beta,
        "temperature_K": 1.0 / beta,
        "logZ": log_z,
        "rho_trace": 1.0,
        "factor_max_bond": int(x_mpo.max_bond()),
        "elapsed_s": 0.0,
    }]

    for step in range(1, squarings + 1):
        start = time.perf_counter()
        # For exact Hermitian X, X @ X and X^dagger @ X are identical.
        # The latter damps anti-Hermitian compression noise at every cooling
        # step while retaining the same beta-doubling recursion.
        lazy_product = x_mpo.H.apply(
            x_mpo,
            which_A="upper",
            which_B="upper",
            contract=False,
            compress=False,
        )
        compress_opts = {}
        if compression_method.startswith("fit"):
            compress_opts.update(
                max_iterations=fit_max_iterations,
                tol=1e-9,
                progbar=False,
            )
        product = qtn.tensor_network_1d_compress(
            lazy_product,
            max_bond=max_bond,
            cutoff=cutoff,
            method=compression_method,
            **compress_opts,
        )
        ratio = factor_norm_squared(product)
        if ratio <= 0:
            raise FloatingPointError(f"nonpositive partition ratio {ratio}")
        log_z = 2 * log_z + math.log(ratio)
        product.multiply_(1 / math.sqrt(ratio))
        x_mpo = product
        beta *= 2
        if step == squarings:
            row, _ = measure_factor(x_mpo, h_mpo, log_z, beta, n, eigenvalues)
        else:
            row = {
                "beta_K_inv": beta,
                "temperature_K": 1.0 / beta,
                "logZ": log_z,
                "rho_trace": 1.0,
                "factor_max_bond": int(x_mpo.max_bond()),
            }
        row.update(step=step, elapsed_s=time.perf_counter() - start)
        rows.append(row)

    diagnostic = {}
    if n <= 10:
        _, rho = measure_factor(x_mpo, h_mpo, log_z, beta, n, eigenvalues)
        dense = np.asarray(rho.to_dense())
        dense /= np.trace(dense)
        hermitian_error = float(np.max(np.abs(dense - dense.conj().T)))
        evals = np.linalg.eigvalsh((dense + dense.conj().T) / 2)
        diagnostic = {
            "final_hermitian_max_abs_error": hermitian_error,
            "final_min_density_eigenvalue": float(evals.min()),
            "final_max_density_eigenvalue": float(evals.max()),
        }
    return rows, diagnostic


def run(
    field_T=1.50,
    max_bonds=(16, 24),
    squarings=18,
    compression_method="direct",
):
    graph = lattice((3, 0), (0, 3))
    bonds = model_bonds(graph, "M0")
    builder, h_mpo = build_spin_mpo(
        graph["N"], bonds, GZ * MU_B_OVER_KB * field_T
    )
    dense_h = np.asarray(builder.build_dense())
    h_mpo_bond_before = int(h_mpo.max_bond())
    h_mpo.compress(cutoff=1e-13)
    h_mpo_bond_after = int(h_mpo.max_bond())
    eigenvalues = np.linalg.eigvalsh(dense_h)
    beta_final = 1 / 0.09
    payload = {
        "purpose": "N=9 exact validation of positive MPO-factor thermal squaring; not a bulk claim",
        "representation": "rho=X_dagger_X; X approximates exp(-beta H/2)/sqrt(Z)",
        "factor_update": "X_new=X_dagger_X before compression",
        "compression_method": compression_method,
        "model": "M0 current 2026 parameters",
        "N": graph["N"],
        "field_T": field_T,
        "target_temperature_K": 0.09,
        "squarings": squarings,
        "initial_beta_K_inv": beta_final / (2**squarings),
        "hamiltonian_max_bond_before_canonical_compression": h_mpo_bond_before,
        "hamiltonian_initial_max_bond": h_mpo_bond_after,
        "hamiltonian_MPO_dense_max_abs_error": float(
            np.max(np.abs(np.asarray(h_mpo.to_dense()) - dense_h))
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "quimb": quimb.__version__,
        },
        "runs": [],
    }
    safe_method = compression_method.replace("-", "_")
    out = ROOT / "results" / f"xtrg_quimb_n9_purified_{safe_method}_validation.json"
    for max_bond in max_bonds:
        start = time.perf_counter()
        rows, diagnostic = purified_square(
            h_mpo,
            graph["N"],
            beta_final,
            squarings=squarings,
            max_bond=max_bond,
            cutoff=1e-12,
            eigenvalues=eigenvalues,
            compression_method=compression_method,
        )
        result = {
            "max_bond": max_bond,
            "wall_s": time.perf_counter() - start,
            "final": rows[-1],
            "diagnostic": diagnostic,
            "trajectory": rows,
        }
        payload["runs"].append(result)
        out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        f = result["final"]
        print(
            f"checkpoint D={max_bond:>3} wall={result['wall_s']:.2f}s "
            f"S={f['entropy_kB_per_spin']:.10f} "
            f"dS={f['errors']['entropy_kB_per_spin']:+.3e} "
            f"dE={f['errors']['energy_K']:+.3e} "
            f"minrho={diagnostic['final_min_density_eigenvalue']:+.3e}",
            flush=True,
        )
    return payload, out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--field", type=float, default=1.50)
    p.add_argument("--max-bonds", type=int, nargs="+", default=[16, 24])
    p.add_argument("--squarings", type=int, default=18)
    p.add_argument(
        "--compression-method",
        choices=["direct", "dm", "zipup", "fit", "fit-zipup", "fit-projector"],
        default="direct",
    )
    a = p.parse_args()
    _, out = run(
        a.field,
        tuple(a.max_bonds),
        a.squarings,
        compression_method=a.compression_method,
    )
    print(out)


if __name__ == "__main__":
    main()
