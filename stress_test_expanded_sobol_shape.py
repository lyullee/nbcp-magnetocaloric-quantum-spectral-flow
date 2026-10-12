"""Apply the declared 90 mK shape envelope to expanded Sobol models."""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np

from bulk_to_engineering_surrogate import transform_field
from screen_model import MU_B_OVER_KB
from stress_test_thermodynamic_shape import (
    EXPERIMENTAL_HC,
    HEAT_SCALE,
    LEFT_WIDTH_SCALE,
    METHODS,
    PEAK_SHIFT_T,
    RIGHT_WIDTH_SCALE,
    deformed_profile,
    profile_at_temperature,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main() -> None:
    with (RESULTS / "interlayer_expanded_sobol_retained.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        retained = list(csv.DictReader(handle))
    shifts = sorted({float(row["critical_field_shift_vs_M3_T"]) for row in retained})
    source_hc = (3 * 1.48 + 1.5 * 0.88) / (4.89 * MU_B_OVER_KB)
    profile = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",", names=True, encoding="utf-8-sig"
    )
    profile = profile[profile["anchor_N"] == 12]
    fields, heat = profile_at_temperature(profile, 0.09)
    source_bstar = float(fields[int(np.argmax(heat))])
    deformations = list(itertools.product(
        PEAK_SHIFT_T, LEFT_WIDTH_SCALE, RIGHT_WIDTH_SCALE, HEAT_SCALE
    ))

    rows: list[dict] = []
    for shift, anchor, method, deformation in itertools.product(
        shifts, EXPERIMENTAL_HC, METHODS, deformations
    ):
        peak_shift, left_scale, right_scale, heat_scale = deformation
        p0_field = float(transform_field(
            source_bstar, source_hc, anchor - shift, method
        ))
        warped_fields, warped_heat = deformed_profile(
            fields, heat, source_hc, anchor, method, peak_shift,
            left_scale, right_scale, heat_scale
        )
        p2_field = float(warped_fields[int(np.argmax(warped_heat))])
        q_at_p0 = float(np.interp(p0_field, warped_fields, warped_heat))
        qmax = float(np.max(warped_heat))
        rows.append({
            "critical_field_shift_vs_M3_T": shift,
            "experimental_Hc_anchor_T": anchor,
            "mapping_method": method,
            "residual_peak_shift_T": peak_shift,
            "left_width_scale": left_scale,
            "right_width_scale": right_scale,
            "heat_scale": heat_scale,
            "delta_B_star_T": p2_field - p0_field,
            "regret_percent": 100.0 * (qmax - q_at_p0) / qmax,
        })

    def passes(row: dict) -> bool:
        return row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0

    by_shift = []
    for shift in shifts:
        subset = [row for row in rows if np.isclose(row["critical_field_shift_vs_M3_T"], shift)]
        by_shift.append({
            "critical_field_shift_vs_M3_T": shift,
            "scenario_count": len(subset),
            "pass_count": sum(passes(row) for row in subset),
            "pass_fraction_descriptive_not_probability": sum(passes(row) for row in subset) / len(subset),
            "minimum_delta_B_star_T": min(row["delta_B_star_T"] for row in subset),
            "minimum_regret_percent": min(row["regret_percent"] for row in subset),
            "all_shape_scenarios_pass": all(passes(row) for row in subset),
        })
    central = [
        row for row in rows
        if row["residual_peak_shift_T"] == 0.0
        and row["left_width_scale"] == 1.0
        and row["right_width_scale"] == 1.0
        and row["heat_scale"] == 1.0
    ]
    summary = {
        "status": "expanded Sobol models under declared 90 mK shape stress",
        "model_count": len(shifts),
        "scenario_count": len(rows),
        "pass_count": sum(passes(row) for row in rows),
        "pass_fraction_descriptive_not_probability": sum(passes(row) for row in rows) / len(rows),
        "central_scenario_count": len(central),
        "central_pass_count": sum(passes(row) for row in central),
        "central_all_pass": all(passes(row) for row in central),
        "by_interlayer_shift": by_shift,
        "minimum_shift_with_all_shape_scenarios_passing_T": min(
            row["critical_field_shift_vs_M3_T"] for row in by_shift
            if row["all_shape_scenarios_pass"]
        ),
        "limitations": [
            "Sobol and deformation fractions are deterministic coverage, not probabilities.",
            "The transported profile remains a surrogate rather than an M4 finite-temperature solution.",
        ],
    }

    with (RESULTS / "interlayer_expanded_sobol_shape_stress.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (RESULTS / "interlayer_expanded_sobol_shape_stress_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
