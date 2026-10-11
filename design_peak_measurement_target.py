"""Turn the 90 mK robustness gap into a concrete field-sweep target.

This is a deterministic value-of-information calculation.  It asks how small
the residual uncertainty in the heat-profile peak must become before every
remaining anchor, mapping, and asymmetric-width scenario has a common W95
field.  Fractions are not interpreted as probabilities.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np

from screen_model import MU_B_OVER_KB
from stress_test_thermodynamic_shape import (
    EXPERIMENTAL_HC,
    LEFT_WIDTH_SCALE,
    METHODS,
    PEAK_SHIFT_T,
    RIGHT_WIDTH_SCALE,
    deformed_profile,
    profile_at_temperature,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
TEMPERATURE_K = 0.09
DECISION_GRID_T = np.linspace(1.25, 1.60, 7001)
UNCERTAINTY_GRID_T = np.arange(0.0, 0.0200001, 0.00025)
SWEEP_STEP_T = 0.004  # guarantees a sampled point within 2 mT


def scenario_matrix(fields, heat, source_hc, peak_uncertainty_t):
    peak_shifts = np.linspace(-peak_uncertainty_t, peak_uncertainty_t, 9)
    curves = []
    for anchor, method, shift, left, right in itertools.product(
        EXPERIMENTAL_HC,
        METHODS,
        peak_shifts,
        LEFT_WIDTH_SCALE,
        RIGHT_WIDTH_SCALE,
    ):
        warped_field, warped_heat = deformed_profile(
            fields, heat, source_hc, anchor, method, shift, left, right, 1.0
        )
        curves.append(
            np.interp(
                DECISION_GRID_T,
                warped_field,
                warped_heat,
                left=np.nan,
                right=np.nan,
            )
            / np.max(warped_heat)
        )
    return np.asarray(curves)


def minimax_row(fields, heat, source_hc, uncertainty):
    normalized = scenario_matrix(fields, heat, source_hc, uncertainty)
    valid = np.all(np.isfinite(normalized), axis=0)
    worst_regret = 1.0 - np.nanmin(normalized, axis=0)
    valid_indices = np.flatnonzero(valid)
    best_index = valid_indices[np.argmin(worst_regret[valid])]
    common_w95 = valid & (worst_regret <= 0.05)
    return {
        "peak_location_uncertainty_T": float(uncertainty),
        "scenario_count": int(normalized.shape[0]),
        "minimax_field_T": float(DECISION_GRID_T[best_index]),
        "worst_case_regret_percent": float(100.0 * worst_regret[best_index]),
        "common_W95_exists": bool(np.any(common_w95)),
        "common_W95_min_T": (
            float(DECISION_GRID_T[common_w95].min())
            if np.any(common_w95)
            else None
        ),
        "common_W95_max_T": (
            float(DECISION_GRID_T[common_w95].max())
            if np.any(common_w95)
            else None
        ),
    }


def main():
    profile = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",",
        names=True,
        encoding="utf-8-sig",
    )
    profile = profile[profile["anchor_N"] == 12]
    fields, heat = profile_at_temperature(profile, TEMPERATURE_K)
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)

    rows = [
        minimax_row(fields, heat, source_hc, uncertainty)
        for uncertainty in UNCERTAINTY_GRID_T
    ]
    feasible = [row for row in rows if row["common_W95_exists"]]
    max_uncertainty = max(row["peak_location_uncertainty_T"] for row in feasible)

    # Bound where the peak may lie under the original +/-20 mT envelope.
    peaks = []
    for anchor, method, shift, left, right in itertools.product(
        EXPERIMENTAL_HC,
        METHODS,
        PEAK_SHIFT_T,
        LEFT_WIDTH_SCALE,
        RIGHT_WIDTH_SCALE,
    ):
        warped_field, warped_heat = deformed_profile(
            fields, heat, source_hc, anchor, method, shift, left, right, 1.0
        )
        peaks.append(float(warped_field[int(np.argmax(warped_heat))]))
    peak_min, peak_max = min(peaks), max(peaks)
    sweep_min = math.floor(peak_min / SWEEP_STEP_T) * SWEEP_STEP_T
    sweep_max = math.ceil(peak_max / SWEEP_STEP_T) * SWEEP_STEP_T
    sweep_fields = np.arange(sweep_min, sweep_max + SWEEP_STEP_T / 2, SWEEP_STEP_T)

    with (RESULTS / "peak_location_uncertainty_scan.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (RESULTS / "recommended_90mK_field_sweep.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence", "field_T"])
        writer.writeheader()
        writer.writerows(
            {"sequence": i + 1, "field_T": float(field)}
            for i, field in enumerate(sweep_fields)
        )

    payload = {
        "status": "deterministic measurement-design target; not a statistical power calculation",
        "temperature_K": TEMPERATURE_K,
        "remaining_scenario_coordinates": {
            "experimental_Hc_T": list(EXPERIMENTAL_HC),
            "mapping_method": list(METHODS),
            "left_width_scale": list(LEFT_WIDTH_SCALE),
            "right_width_scale": list(RIGHT_WIDTH_SCALE),
        },
        "largest_symmetric_peak_uncertainty_with_common_W95_T": max_uncertainty,
        "worst_case_peak_support_under_original_envelope_T": [peak_min, peak_max],
        "recommended_no_fit_sweep": {
            "start_T": float(sweep_fields[0]),
            "stop_T": float(sweep_fields[-1]),
            "step_T": SWEEP_STEP_T,
            "point_count": int(len(sweep_fields)),
            "maximum_grid_localization_error_T": SWEEP_STEP_T / 2,
        },
        "interpretation": (
            "At 90 mK, localizing the heat-profile peak to +/-2 mT or better "
            "restores a common W95 field across the retained anchor, mapping, "
            "and asymmetric-width envelope.  A 4 mT no-fit sweep across the "
            "full predicted peak support meets that discretization target."
        ),
        "limitations": [
            "The sweep spacing controls discretization error only; instrument field calibration and heat/entropy uncertainty must also be below the target.",
            "A smooth peak fit could use fewer points, but its guarantee would depend on the chosen line shape.",
            "The target applies to the declared deformation envelope and does not replace a full ADR cycle measurement.",
        ],
    }
    (RESULTS / "peak_measurement_target_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
