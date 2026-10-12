"""Low-discrepancy audit of a broader interlayer-parameter hyperrectangle.

The original 375-point Cartesian grid is deliberately local. This independent
Sobol design tests simultaneous parameter changes and disconnected admissible
regions at fixed computational cost. It is a deterministic space-filling
stress design, not a posterior sample and not an estimate of probability.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from bulk_polarized_stability import (
    EwaldSettings,
    InterlayerParameters,
    locate_critical_mode,
    locate_critical_mode_multistart,
)
from screen_interlayer_parameter_ensemble import (
    DELTA_BAND,
    GRID,
    KM_LINE_RESIDUAL_MAX,
    QZ_BAND,
    canonical_soft_mode,
    engineering_rows,
    out_of_plane_bandwidth_microev,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
SETTINGS = EwaldSettings(boundary="tin_foil")
PARAMETERS = ("j2_microev", "ja_microev", "jb_microev", "delta_c")
BOUNDS = {
    "j2_microev": (0.0, 5.6),
    "ja_microev": (-4.2, 4.2),
    "jb_microev": (1.975, 13.825),
    "delta_c": (0.02, 0.14),
}
N_POINTS = 256
SEED = 261012
EXTRA_OPTIMIZER_SEEDS = (271, 27182)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def in_original_box(row: dict) -> bool:
    return all(min(GRID[key]) <= float(row[key]) <= max(GRID[key]) for key in PARAMETERS)


def range_of(rows: list[dict], key: str) -> list[float] | None:
    if not rows:
        return None
    values = [float(row[key]) for row in rows]
    return [min(values), max(values)]


def main() -> None:
    sampler = qmc.Sobol(d=len(PARAMETERS), scramble=True, seed=SEED)
    unit = sampler.random_base2(m=8)
    lower = np.array([BOUNDS[key][0] for key in PARAMETERS])
    upper = np.array([BOUNDS[key][1] for key in PARAMETERS])
    points = qmc.scale(unit, lower, upper)

    m3_hc = float(locate_critical_mode("M3", SETTINGS)["critical_field_T"])
    rows: list[dict] = []
    for index, values in enumerate(points, start=1):
        params = dict(zip(PARAMETERS, map(float, values), strict=True))
        parameters = InterlayerParameters(**params)
        result = locate_critical_mode_multistart(
            "M4", SETTINGS,
            seeds=(SEED + index, *EXTRA_OPTIMIZER_SEEDS),
            interlayer=parameters,
        )
        mode = canonical_soft_mode(result["critical_wavevector_rlu"])
        delta_pass = DELTA_BAND[0] <= mode["delta_inplane_rlu"] <= DELTA_BAND[1]
        qz_pass = QZ_BAND[0] <= mode["qz_folded_rlu"] <= QZ_BAND[1]
        km_pass = mode["km_line_residual_rlu"] <= KM_LINE_RESIDUAL_MAX
        rows.append({
            "sample_id": index,
            **params,
            "critical_field_T": float(result["critical_field_T"]),
            "critical_field_shift_vs_M3_T": float(result["critical_field_T"] - m3_hc),
            "q_u_rlu": float(result["critical_wavevector_rlu"][0]),
            "q_v_rlu": float(result["critical_wavevector_rlu"][1]),
            "q_w_rlu": float(result["critical_wavevector_rlu"][2]),
            "delta_inplane_rlu": mode["delta_inplane_rlu"],
            "qz_folded_rlu": mode["qz_folded_rlu"],
            "km_line_residual_rlu": mode["km_line_residual_rlu"],
            "out_of_plane_bandwidth_microeV": out_of_plane_bandwidth_microev(parameters),
            "delta_pass": delta_pass,
            "qz_pass": qz_pass,
            "km_pass": km_pass,
            "propagation_vector_pass": delta_pass and qz_pass and km_pass,
            "inside_original_parameter_box": in_original_box(params),
            "multistart_selected_seed": result["multistart"]["selected_seed"],
            "multistart_min_field_T": result["multistart"]["critical_field_range_T"][0],
            "multistart_field_spread_T": (
                result["multistart"]["critical_field_range_T"][1]
                - result["multistart"]["critical_field_range_T"][0]
            ),
        })
        if index % 16 == 0:
            print(f"completed expanded Sobol audit {index}/{N_POINTS}", flush=True)

    retained = [row for row in rows if row["propagation_vector_pass"]]
    engineering = engineering_rows(retained, m3_hc) if retained else []
    e90 = [row for row in engineering if row["temperature_K"] == 0.09]
    outside_retained = [row for row in retained if not row["inside_original_parameter_box"]]

    summary = {
        "status": "deterministic expanded Sobol stress design; not a probability distribution",
        "seed": SEED,
        "optimizer_multistart": {
            "starts_per_point": 1 + len(EXTRA_OPTIMIZER_SEEDS),
            "point_specific_seed": "261012 + sample_id",
            "additional_seeds": list(EXTRA_OPTIMIZER_SEEDS),
            "selection_rule": "largest critical field (first mode to soften on lowering field)",
        },
        "point_count": N_POINTS,
        "bounds": {key: list(value) for key, value in BOUNDS.items()},
        "original_box": {key: [min(GRID[key]), max(GRID[key])] for key in PARAMETERS},
        "constraint_counts": {
            "delta": sum(row["delta_pass"] for row in rows),
            "qz": sum(row["qz_pass"] for row in rows),
            "km": sum(row["km_pass"] for row in rows),
            "all": len(retained),
        },
        "retained_count": len(retained),
        "retained_outside_original_box": len(outside_retained),
        "retained_ranges": {
            "critical_field_shift_vs_M3_T": range_of(retained, "critical_field_shift_vs_M3_T"),
            "delta_inplane_rlu": range_of(retained, "delta_inplane_rlu"),
            "qz_folded_rlu": range_of(retained, "qz_folded_rlu"),
            "out_of_plane_bandwidth_microeV": range_of(retained, "out_of_plane_bandwidth_microeV"),
            "90mK_delta_B_star_T": range_of(e90, "delta_B_star_T"),
            "90mK_regret_percent": range_of(e90, "regret_percent"),
        },
        "central_shape_threshold_coverage": {
            "scenario_count": len(e90),
            "both_gates_pass_count": sum(
                row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
                for row in e90
            ),
            "all_pass": bool(e90) and all(
                row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
                for row in e90
            ),
        },
        "interpretation_limits": [
            "Sobol coverage is space filling but finite and cannot prove that no disconnected admissible region exists.",
            "Retained fractions are deterministic design coverage, not probabilities.",
            "The near-saturation wavevector does not validate the full finite-temperature entropy.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    write_csv(RESULTS / "interlayer_expanded_sobol.csv", rows)
    write_csv(RESULTS / "interlayer_expanded_sobol_retained.csv", retained)
    write_csv(RESULTS / "interlayer_expanded_sobol_engineering.csv", engineering)
    (RESULTS / "interlayer_expanded_sobol_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
