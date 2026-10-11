"""Dominant-state adiabatic coupling diagnostics along N=9 M0-to-M3."""
from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np

from material_models import model_bonds
from matrix_free_hamiltonian import SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, lattice
from screen_weak_terms import GZ


ROOT = Path(__file__).resolve().parent
FIELDS_T = (1.40, 1.45, 1.50)
TEMPERATURE_K = 0.10
LAMBDA_GRID = np.linspace(0, 1, 101)
COUPLING_THRESHOLD_K = 1e-9
THERMAL_STATE_PROBABILITY_THRESHOLD = 1e-8


def dense_block(n, bonds, field_k, parity):
    operator = SpinHamiltonianOperator(n, bonds, field_k, parity=parity)
    return np.asarray(operator @ np.eye(operator.shape[0], dtype=complex))


def main():
    started = time.perf_counter()
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    bonds0 = model_bonds(graph, "M0")
    bonds3 = model_bonds(graph, "M3", boundary="tin_foil")
    rows = []
    for field_t in FIELDS_T:
        field_k = GZ * MU_B_OVER_KB * field_t
        blocks = {}
        all_initial = []
        for parity in (0, 1):
            h0 = dense_block(n, bonds0, field_k, parity)
            h1 = dense_block(n, bonds3, field_k, parity)
            blocks[parity] = (h0, h1 - h0)
            all_initial.extend((float(value), parity, rank)
                               for rank, value in enumerate(np.linalg.eigvalsh(h0)))
        global_min = min(item[0] for item in all_initial)
        weights = np.array([np.exp(-(item[0] - global_min) / TEMPERATURE_K)
                            for item in all_initial])
        weights /= weights.sum()
        parity_weights = {parity: float(sum(weight for weight, item in zip(weights, all_initial)
                                            if item[1] == parity))
                          for parity in (0, 1)}
        ground_weights = {parity: float(sum(weight for weight, item in zip(weights, all_initial)
                                            if item[1] == parity and item[2] == 0))
                          for parity in (0, 1)}

        parity_results = {}
        for parity, (h0, v) in blocks.items():
            max_ratio = 0.0
            max_ratio_lambda = None
            min_coupled_gap = float("inf")
            min_coupled_gap_lambda = None
            min_raw_gap = float("inf")
            max_all_spectrum_ratio = 0.0
            max_thermal_spectrum_ratio = 0.0
            all_ratio_lambda = None
            thermal_ratio_lambda = None
            for lam in LAMBDA_GRID:
                energies, vectors = np.linalg.eigh(h0 + float(lam) * v)
                gaps = energies[1:] - energies[0]
                min_raw_gap = min(min_raw_gap, float(np.min(gaps)))
                couplings = np.abs(vectors[:, 1:].conj().T @ (v @ vectors[:, 0]))
                mask = couplings > COUPLING_THRESHOLD_K
                if np.any(mask):
                    local_min = float(np.min(gaps[mask]))
                    if local_min < min_coupled_gap:
                        min_coupled_gap = local_min
                        min_coupled_gap_lambda = float(lam)
                    ratios = couplings[mask] / gaps[mask]**2
                    local_ratio = float(np.max(ratios))
                    if local_ratio > max_ratio:
                        max_ratio = local_ratio
                        max_ratio_lambda = float(lam)
                # Audit every coupled level pair.  Columns are possible initial
                # levels and rows are coupled destination levels.
                v_eigen = np.abs(vectors.conj().T @ (v @ vectors))
                pair_gaps = np.abs(energies[:, None] - energies[None, :])
                pair_mask = ((v_eigen > COUPLING_THRESHOLD_K) &
                             (pair_gaps > 1e-12))
                ratios_all = np.zeros_like(pair_gaps)
                ratios_all[pair_mask] = v_eigen[pair_mask] / pair_gaps[pair_mask]**2
                local_all = float(np.max(ratios_all))
                if local_all > max_all_spectrum_ratio:
                    max_all_spectrum_ratio = local_all
                    all_ratio_lambda = float(lam)
                instantaneous_probability = np.exp(-(energies - energies[0]) / TEMPERATURE_K)
                instantaneous_probability /= instantaneous_probability.sum()
                initial_mask = instantaneous_probability > THERMAL_STATE_PROBABILITY_THRESHOLD
                thermal_mask = pair_mask & initial_mask[None, :]
                ratios_thermal = np.zeros_like(pair_gaps)
                ratios_thermal[thermal_mask] = (v_eigen[thermal_mask] /
                                                pair_gaps[thermal_mask]**2)
                local_thermal = float(np.max(ratios_thermal))
                if local_thermal > max_thermal_spectrum_ratio:
                    max_thermal_spectrum_ratio = local_thermal
                    thermal_ratio_lambda = float(lam)
            parity_results[str(parity)] = {
                "initial_thermal_parity_weight": parity_weights[parity],
                "initial_ground_state_probability": ground_weights[parity],
                "minimum_raw_ground_gap_K": min_raw_gap,
                "minimum_V_coupled_ground_gap_K": min_coupled_gap,
                "minimum_coupled_gap_lambda": min_coupled_gap_lambda,
                "maximum_abs_Vmn_over_gap2_inverse_K": max_ratio,
                "maximum_ratio_lambda": max_ratio_lambda,
                "indicative_ramp_time_for_1pct_amplitude_inverse_K": max_ratio / 0.01,
                "maximum_all_spectrum_abs_Vmn_over_gap2_inverse_K": max_all_spectrum_ratio,
                "maximum_all_spectrum_ratio_lambda": all_ratio_lambda,
                "maximum_thermal_spectrum_abs_Vmn_over_gap2_inverse_K": max_thermal_spectrum_ratio,
                "maximum_thermal_spectrum_ratio_lambda": thermal_ratio_lambda,
                "indicative_thermal_spectrum_ramp_time_for_1pct_inverse_K": max_thermal_spectrum_ratio / 0.01,
            }
        rows.append({
            "field_T": field_t,
            "temperature_K": TEMPERATURE_K,
            "lambda_grid_points": len(LAMBDA_GRID),
            "parities": parity_results,
            "combined_initial_ground_probability": sum(ground_weights.values()),
        })
        print(f"B={field_t:.2f}: ratios p0={parity_results['0']['maximum_abs_Vmn_over_gap2_inverse_K']:.2e}, "
              f"p1={parity_results['1']['maximum_abs_Vmn_over_gap2_inverse_K']:.2e}", flush=True)

    payload = {
        "purpose": "dominant-state adiabatic depth diagnostic for the low-variance Jarzynski endpoint",
        "cluster": {"N": n, "translations": graph["translations"], "boundary": "periodic; M3 tin-foil Ewald"},
        "path": "H(lambda)=M0+lambda*(M3-M0)",
        "rows": rows,
        "scope_limits": [
            "The ratio |<m|dH/dlambda|0>|/gap^2 is an adiabatic diagnostic, not a rigorous complete runtime bound.",
            "Only the lowest state in each conserved spin-parity sector is audited.",
            "Translation and other symmetry-protected crossings are suppressed by the coupling matrix element threshold.",
            "Excited states that contribute to high moments of exp(-beta W) can require slower ramps.",
            "The thermal-spectrum audit uses instantaneous Boltzmann probability >1e-8 and does not exactly reproduce Jarzynski rare-event importance.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "adiabatic_path_diagnostic_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(output), "rows": rows,
                      "elapsed_seconds": payload["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
