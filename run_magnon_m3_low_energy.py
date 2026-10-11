"""Low-energy M3-M0 correction in controlled projected-magnon spaces."""
import argparse
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
CELLS = {
    21: ((-5, 1), (-1, -4), [5, 7]),
    24: ((-5, 1), (1, -5), [6, 7]),
    27: ((-6, 3), (-3, -3), [6, 7]),
}
TEMPERATURES = [0.09, 0.10, 0.12]


def low_entropy(energies, temperature, n):
    energies = np.sort(np.asarray(energies, float))
    shifted = energies - energies[0]
    w = np.exp(-shifted / temperature)
    return float((np.log(w.sum()) + np.dot(w, shifted) / w.sum() / temperature) / n)


def main(n: int):
    if n not in CELLS:
        raise ValueError(f"supported N: {sorted(CELLS)}")
    a, b, cutoffs = CELLS[n]
    graph = lattice(a, b)
    field_t = 1.50
    field_k = GZ * MU_B_OVER_KB * field_t
    started = time.perf_counter()
    runs = []
    cache = {}
    bonds = {model: model_bonds(graph, model) for model in ("M0", "M3")}
    print(f"N={n} Ewald/model bonds ready at {time.perf_counter()-started:.1f}s", flush=True)
    for cutoff in cutoffs:
        by_model = {}
        for model in ("M0", "M3"):
            energies = []
            parity_runs = []
            for parity in (0, 1):
                parity_cutoff = cutoff if cutoff % 2 == parity else cutoff - 1
                key = (model, parity, parity_cutoff)
                if key not in cache:
                    op = TruncatedMagnonCartesianOperator(n, cutoff, parity, bonds[model], field_k)
                    count = min(8, op.shape[0] - 1)
                    vals = (np.array([op.diagonal[0]]) if op.shape[0] == 1 else
                            eigsh(op, k=count, which="SA", tol=1e-9, return_eigenvectors=False))
                    vals = np.sort(np.real(vals))
                    cache[key] = {"down_magnon_parity": parity, "effective_max_q_down": parity_cutoff,
                                  "dimension": op.shape[0], "low_energies_K": vals.tolist()}
                    print(f"N={n} K={cutoff} {model} parity={parity} dim={op.shape[0]} E0={vals[0]:.9f} elapsed={time.perf_counter()-started:.1f}s", flush=True)
                else:
                    print(f"N={n} K={cutoff} {model} parity={parity} reused", flush=True)
                record = cache[key]
                energies.extend(record["low_energies_K"])
                parity_runs.append(record)
            by_model[model] = {
                "parities": parity_runs,
                "combined_low_energies_K": sorted(energies),
                "low_entropy_kB_per_spin": {str(t): low_entropy(energies, t, n) for t in TEMPERATURES},
            }
        rows = []
        for t in TEMPERATURES:
            s0 = by_model["M0"]["low_entropy_kB_per_spin"][str(t)]
            s3 = by_model["M3"]["low_entropy_kB_per_spin"][str(t)]
            rows.append({"temperature_K": t, "S_M0_low16": s0, "S_M3_low16": s3, "deltaS_low16": s3 - s0})
        runs.append({"max_q_down": cutoff, "models": by_model, "rows": rows})
    payload = {
        "method": "lowest 8 states in each magnon-parity block of the projected Cartesian Hamiltonian",
        "condition": {"N": n, "translations": [a, b], "field_T": field_t, "boundary": "tin_foil"},
        "validation_basis": "At N=18 and 90 mK, low-16 deltaS differed from full FTLM by 1.45e-4 kB/spin.",
        "runs": runs,
        "elapsed_seconds": time.perf_counter() - started,
    }
    out = RESULTS / f"magnon_M3_low_energy_N{n}_B1.500T.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": out.name, "runs": [{"K": r["max_q_down"], "rows": r["rows"]} for r in runs]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--n", type=int, required=True)
    main(parser.parse_args().n)
