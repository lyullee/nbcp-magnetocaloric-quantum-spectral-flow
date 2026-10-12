"""Targeted N=18 convergence audit for the strongest stacked M4 deformation.

The N=12 calculation establishes that the rigid critical-field transport can
fail, but its 3x2x2 seam is short.  This script uses the three-sublattice-
compatible 3x3x2 cluster and exact low-energy spectra in controlled projected
magnon spaces.  It is resumable because every (model, field, cutoff) spectrum
is stored separately.

The cold 90 mK entropy uses explicitly converged low eigenvalues.  The 2 K,
4 T reference entropy is handled separately by low-energy-corrected FTLM; it
must not be inferred from the cold low-state list.
"""
from __future__ import annotations

import argparse
import json
import pickle
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ftlm import lowest_eigenpairs, run_ftlm
from ewald_cluster import MICROEV_TO_K
from matrix_free_hamiltonian import SpinBond
from screen_model import MU_B_OVER_KB, VOLUME_FACTOR_MJ_PER_CM3
from screen_weak_terms import GZ
from truncated_magnon_operator import TruncatedMagnonCartesianOperator


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CACHE = RESULTS / "n18_stacked_m4_work"
COMMON_CACHE = RESULTS / "n18_stacked_common_bonds.pkl"
OUT = RESULTS / "n18_stacked_m4_convergence_summary.json"
TEMPERATURE = 0.09


@dataclass(frozen=True)
class InterlayerParameters:
    j2_microev: float
    ja_microev: float
    jb_microev: float
    delta_c: float


MODELS = {
    "M3": InterlayerParameters(0.0, 0.0, 0.0, 0.08),
    "M4-01": InterlayerParameters(
        4.5861779235303395,
        3.4409039821475744,
        13.63851090040989,
        0.03836957771331072,
    ),
}


def entropy_from_energies(energies, temperature=TEMPERATURE, spins=18):
    energies = np.sort(np.asarray(energies, dtype=float))
    shifted = energies - energies[0]
    weights = np.exp(-shifted / temperature)
    z = float(weights.sum())
    return float((np.log(z) + np.dot(weights, shifted) / z / temperature) / spins)


def load_common():
    if COMMON_CACHE.exists():
        graph, bonds = pickle.loads(COMMON_CACHE.read_bytes())
    else:
        # Building the Ewald cache requires the broader stability module.  It
        # is intentionally lazy so ordinary resumed runs only need SciPy's
        # linear-algebra path.
        from audit_m4_finite_temperature_cluster import build_common_bonds
        graph, bonds = build_common_bonds(3, 3, 2)
        COMMON_CACHE.write_bytes(pickle.dumps((graph, bonds)))
    return graph, bonds


def add_interlayer_bonds(graph, common, parameters):
    bonds = list(common)
    values = {
        "J2": parameters.j2_microev,
        "Ja": parameters.ja_microev,
        "Jb": parameters.jb_microev,
    }
    for bond in graph["bonds"]:
        if bond.kind == "J1":
            continue
        coupling = values[bond.kind]
        tensor = np.diag([
            parameters.delta_c * coupling,
            parameters.delta_c * coupling,
            coupling,
        ]) * MICROEV_TO_K
        bonds.append(SpinBond(bond.i, bond.j, tensor))
    return bonds


def spectrum_path(model, field, even_cutoff, odd_cutoff, states):
    return CACHE / (
        f"{model}_B{field:.4f}_K{even_cutoff}-{odd_cutoff}_L{states}.json"
    )


def cold_spectrum(graph, common, model, field, even_cutoff, odd_cutoff, states):
    CACHE.mkdir(exist_ok=True)
    path = spectrum_path(model, field, even_cutoff, odd_cutoff, states)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    bonds = add_interlayer_bonds(graph, common, MODELS[model])
    records = []
    combined = []
    started = time.perf_counter()
    for parity, cutoff in ((0, even_cutoff), (1, odd_cutoff)):
        operator = TruncatedMagnonCartesianOperator(
            graph["N"], cutoff, parity, bonds,
            GZ * MU_B_OVER_KB * float(field),
        )
        part_started = time.perf_counter()
        values, _ = lowest_eigenpairs(
            operator, states, tolerance=1e-10,
            seed=20261012 + 1000 * parity + int(round(1000 * field)),
        )
        combined.extend(values.tolist())
        records.append({
            "parity": parity,
            "cutoff": cutoff,
            "dimension": int(operator.shape[0]),
            "elapsed_seconds": time.perf_counter() - part_started,
            "energies_K": values.tolist(),
        })
    payload = {
        "model": model,
        "field_T": float(field),
        "temperature_K": TEMPERATURE,
        "even_cutoff": even_cutoff,
        "odd_cutoff": odd_cutoff,
        "low_states_per_parity": states,
        "entropy_kB_per_spin": entropy_from_energies(combined),
        "records": records,
        "elapsed_seconds": time.perf_counter() - started,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def scan(graph, common, fields, even_cutoff, odd_cutoff, states):
    rows = []
    for model in MODELS:
        for index, field in enumerate(fields):
            row = cold_spectrum(
                graph, common, model, float(field),
                even_cutoff, odd_cutoff, states,
            )
            rows.append(row)
            print(json.dumps({
                "model": model,
                "completed": index + 1,
                "total": len(fields),
                "field_T": float(field),
                "entropy_kB_per_spin": row["entropy_kB_per_spin"],
                "elapsed_seconds": row["elapsed_seconds"],
            }), flush=True)
    return rows


def hot_entropy(graph, common, model, even_cutoff, odd_cutoff):
    path = CACHE / f"{model}_hot4T_K{even_cutoff}-{odd_cutoff}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    bonds = add_interlayer_bonds(graph, common, MODELS[model])
    operators = [
        TruncatedMagnonCartesianOperator(
            graph["N"], cutoff, parity, bonds, GZ * MU_B_OVER_KB * 4.0,
        )
        for parity, cutoff in ((0, even_cutoff), (1, odd_cutoff))
    ]
    result = run_ftlm(
        operators, [2.0], graph["N"], random_vectors=4,
        krylov_dim=60, exact_low_states=12,
        seed=20261012,
    )
    payload = {
        "model": model,
        "field_T": 4.0,
        "temperature_K": 2.0,
        "even_cutoff": even_cutoff,
        "odd_cutoff": odd_cutoff,
        "entropy_kB_per_spin": result["rows"][0]["entropy_kB_per_spin"],
        "ftlm": result,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def summarize(rows, hot_rows=None):
    summary = {}
    for model in MODELS:
        selected = [row for row in rows if row["model"] == model]
        selected.sort(key=lambda row: row["field_T"])
        peak = max(selected, key=lambda row: row["entropy_kB_per_spin"])
        item = {
            "peak_field_on_grid_T": peak["field_T"],
            "peak_cold_entropy_kB_per_spin": peak["entropy_kB_per_spin"],
            "fields_T": [row["field_T"] for row in selected],
            "cold_entropy_kB_per_spin": [row["entropy_kB_per_spin"] for row in selected],
        }
        if hot_rows:
            hot = hot_rows[model]["entropy_kB_per_spin"]
            heats = [
                TEMPERATURE * (row["entropy_kB_per_spin"] - hot)
                * VOLUME_FACTOR_MJ_PER_CM3 for row in selected
            ]
            qpeak = int(np.argmax(heats))
            item.update({
                "hot_entropy_kB_per_spin": hot,
                "heat_mJ_per_cm3": heats,
                "Qmax_mJ_per_cm3": heats[qpeak],
                "Q_peak_field_on_grid_T": selected[qpeak]["field_T"],
            })
        summary[model] = item
    if all(model in summary for model in ("M3", "M4-01")):
        summary["comparison"] = {
            "M4_minus_M3_peak_field_T": (
                summary["M4-01"]["peak_field_on_grid_T"]
                - summary["M3"]["peak_field_on_grid_T"]
            ),
        }
        if hot_rows:
            summary["comparison"]["Qmax_ratio_M4_to_M3"] = (
                summary["M4-01"]["Qmax_mJ_per_cm3"]
                / summary["M3"]["Qmax_mJ_per_cm3"]
            )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("coarse", "refine", "hot"), default="coarse")
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    graph, common = load_common()
    if graph["N"] != 18 or tuple(graph["shape"]) != (3, 3, 2):
        raise RuntimeError("common-bond cache is not the 3x3x2 cluster")

    if args.stage == "coarse":
        fields = np.round(np.arange(1.10, 2.0001, 0.10), 10)
        rows = scan(graph, common, fields, 4, 5, 8)
        payload = {
            "status": "coarse low-magnon location scan; not a convergence result",
            "cluster": {"shape": [3, 3, 2], "spins": 18},
            "settings": {"even_cutoff": 4, "odd_cutoff": 5, "low_states_per_parity": 8},
            "summary": summarize(rows),
        }
        OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    elif args.stage == "refine":
        coarse = json.loads(OUT.read_text(encoding="utf-8"))
        rows = []
        for model in MODELS:
            center = coarse["summary"][model]["peak_field_on_grid_T"]
            fields = np.round(np.arange(center - 0.08, center + 0.0801, 0.02), 10)
            for index, field in enumerate(fields):
                row = cold_spectrum(graph, common, model, float(field), 6, 7, 12)
                rows.append(row)
                print(json.dumps({
                    "model": model, "completed": index + 1, "total": len(fields),
                    "field_T": float(field),
                    "entropy_kB_per_spin": row["entropy_kB_per_spin"],
                    "elapsed_seconds": row["elapsed_seconds"],
                }), flush=True)
        hot_rows = {model: hot_entropy(graph, common, model, 6, 7) for model in MODELS}
        payload = {
            "status": "targeted N=18 three-sublattice-compatible finite-temperature convergence audit",
            "cluster": {
                "shape": [3, 3, 2], "spins": 18, "boundary": "tin_foil",
                "three_sublattice_compatible": True,
            },
            "settings": {"even_cutoff": 6, "odd_cutoff": 7, "low_states_per_parity": 12},
            "summary": summarize(rows, hot_rows),
        }
        OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    else:
        payload = {
            model: {
                f"K{even}-{odd}": hot_entropy(graph, common, model, even, odd)
                for even, odd in ((4, 5), (6, 7), (8, 9))
            }
            for model in MODELS
        }
        (RESULTS / "n18_stacked_m4_hot_entropy_convergence.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
        print(json.dumps({
            model: {key: value["entropy_kB_per_spin"] for key, value in levels.items()}
            for model, levels in payload.items()
        }, indent=2))


if __name__ == "__main__":
    main()
