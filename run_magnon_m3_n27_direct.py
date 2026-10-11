"""One decisive direct N=27 M3 low-energy calculation at 1.50 T.

M0 fixed-q eigenvalues are reused from the completed size scan.  M3 is solved
directly in even q<=6 and odd q<=7 parity spaces.  A checkpoint is written after
each parity block.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import eigsh

from material_models import model_bonds
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ
from truncated_magnon_operator import TruncatedMagnonCartesianOperator

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "magnon_M3_low_energy_N27_B1.500T.json"
TEMPERATURES = (0.09, 0.10, 0.12)


def entropy(energies, temperature, n=27):
    e = np.sort(np.asarray(energies, float))
    shifted = e - e[0]
    w = np.exp(-shifted / temperature)
    return float((np.log(w.sum()) + np.dot(w, shifted) / w.sum() / temperature) / n)


def write_checkpoint(payload):
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main():
    n = 27
    graph = lattice((-6, 3), (-3, -3))
    bonds_by_model = {model: model_bonds(graph, model) for model in ("M0", "M3")}
    field_t = 1.50
    field_k = GZ * MU_B_OVER_KB * field_t
    started = time.perf_counter()
    payload = {
        "method": "direct eigsh lowest 8 states in each projected M3 magnon-parity space",
        "condition": {"N": n, "translations": graph["translations"], "field_T": field_t, "boundary": "tin_foil"},
        "settings": {"max_q_even": 6, "max_q_odd": 7, "low_states_per_parity": 8, "eigsh_tolerance": 1e-9, "seed_base": 20261011},
        "M0_source": "direct parity eigsh; absolute M0 thermodynamics cross-reference results/magnon_m0_N27_B1.500T.json",
        "models": {"M0": {"parities": []}, "M3": {"parities": []}},
        "status": "running",
    }
    for model in ("M0", "M3"):
        for parity, cutoff in ((0, 6), (1, 7)):
            op = TruncatedMagnonCartesianOperator(n, cutoff, parity, bonds_by_model[model], field_k)
            rng = np.random.default_rng(20261011 + 100 * (model == "M3") + parity)
            v0 = rng.normal(size=op.shape[0]) + 1j * rng.normal(size=op.shape[0])
            values = eigsh(op, k=8, which="SA", tol=1e-9, v0=v0, return_eigenvectors=False)
            values = np.sort(values.real)
            record = {
                "down_magnon_parity": parity,
                "effective_max_q_down": cutoff,
                "dimension": op.shape[0],
                "low_energies_K": values.tolist(),
                "elapsed_seconds_cumulative": time.perf_counter() - started,
            }
            payload["models"][model]["parities"].append(record)
            payload["elapsed_seconds"] = time.perf_counter() - started
            write_checkpoint(payload)
            print(
                f"N=27 {model} parity={parity} K={cutoff} dim={op.shape[0]} "
                f"E0={values[0]:.9f} elapsed={payload['elapsed_seconds']:.1f}s checkpointed",
                flush=True,
            )

    combined = {}
    for model in ("M0", "M3"):
        energies = sorted(sum((block["low_energies_K"] for block in payload["models"][model]["parities"]), []))
        combined[model] = energies
        payload["models"][model]["combined_low_energies_K"] = energies
        payload["models"][model]["low_entropy_kB_per_spin"] = {str(t): entropy(energies, t) for t in TEMPERATURES}
    payload["rows"] = []
    for t in TEMPERATURES:
        s0 = payload["models"]["M0"]["low_entropy_kB_per_spin"][str(t)]
        s3 = payload["models"]["M3"]["low_entropy_kB_per_spin"][str(t)]
        payload["rows"].append({"temperature_K": t, "S_M0_low16": s0, "S_M3_low16": s3, "deltaS_low16": s3 - s0})
    payload["status"] = "complete"
    payload["elapsed_seconds"] = time.perf_counter() - started
    write_checkpoint(payload)
    print(json.dumps({"output": OUT.name, "rows": payload["rows"], "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
