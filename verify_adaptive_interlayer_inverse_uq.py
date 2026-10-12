"""Independent checks for the adaptive interlayer inverse-UQ audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from adaptive_interlayer_inverse_uq import (
    DELTA_BAND,
    KM_LINE_RESIDUAL_MAX,
    OBSERVABLES,
    PARAMETERS,
    QZ_BAND,
    SETTINGS,
    exact_pass,
)
from bulk_polarized_stability import (
    InterlayerParameters,
    locate_critical_mode,
    locate_critical_mode_multistart,
)
from screen_interlayer_parameter_ensemble import canonical_soft_mode, engineering_rows


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SUMMARY_PATH = RESULTS / "interlayer_adaptive_inverse_uq_summary.json"
ALL_PATH = RESULTS / "interlayer_adaptive_inverse_uq.csv"
RETAINED_PATH = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"
OUTPUT_PATH = RESULTS / "interlayer_adaptive_inverse_uq_verification.json"
INDEPENDENT_SEEDS = (314159, 271828, 161803, 141421, 173205)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return str(value).strip().lower() == "true"


def obs_array(rows: list[dict]) -> np.ndarray:
    return np.array([[float(row[key]) for key in OBSERVABLES] for row in rows])


def source_rows() -> list[dict]:
    rows: list[dict] = []
    for name in (
        "interlayer_parameter_ensemble.csv",
        "interlayer_expanded_sobol.csv",
    ):
        rows.extend(read_csv(RESULTS / name))
    return rows


def max_source_difference(combined: list[dict], source: list[dict]) -> float:
    keys = (*PARAMETERS, "critical_field_T", *OBSERVABLES)
    differences = []
    for current, original in zip(combined[:len(source)], source, strict=True):
        differences.extend(
            abs(float(current[key]) - float(original[key])) for key in keys
        )
    return float(max(differences, default=0.0))


def independent_extreme_check(model: dict, m3_hc: float) -> dict:
    parameters = InterlayerParameters(**{
        key: float(model[key]) for key in PARAMETERS
    })
    result = locate_critical_mode_multistart(
        "M4", SETTINGS, seeds=INDEPENDENT_SEEDS, interlayer=parameters
    )
    mode = canonical_soft_mode(result["critical_wavevector_rlu"])
    recalculated = {
        "critical_field_shift_vs_M3_T": float(result["critical_field_T"] - m3_hc),
        "delta_inplane_rlu": mode["delta_inplane_rlu"],
        "qz_folded_rlu": mode["qz_folded_rlu"],
        "km_line_residual_rlu": mode["km_line_residual_rlu"],
    }
    differences = {
        key: abs(recalculated[key] - float(model[key])) for key in recalculated
    }
    return {
        "independent_seed_count": len(INDEPENDENT_SEEDS),
        "recalculated": recalculated,
        "absolute_differences": differences,
        "field_and_mode_match": (
            differences["critical_field_shift_vs_M3_T"] < 1.0e-8
            and differences["delta_inplane_rlu"] < 1.0e-7
            and differences["qz_folded_rlu"] < 1.0e-7
            and differences["km_line_residual_rlu"] < 1.0e-7
        ),
        "still_public_q_feasible": bool(exact_pass(np.array([[
            recalculated["delta_inplane_rlu"],
            recalculated["qz_folded_rlu"],
            recalculated["km_line_residual_rlu"],
            recalculated["critical_field_shift_vs_M3_T"],
        ]]))[0]),
        "multistart_field_spread_T": float(
            result["multistart"]["critical_field_range_T"][1]
            - result["multistart"]["critical_field_range_T"][0]
        ),
    }


def main() -> None:
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    rows = read_csv(ALL_PATH)
    retained = read_csv(RETAINED_PATH)
    source = source_rows()

    labels_from_observables = exact_pass(obs_array(rows))
    labels_in_file = np.array([
        as_bool(row["propagation_vector_pass"]) for row in rows
    ])
    shift = np.array([
        float(row["critical_field_shift_vs_M3_T"]) for row in retained
    ])
    minimum_model = summary["retained"]["minimum_shift_model"]
    maximum_model = summary["retained"]["maximum_shift_model"]
    m3_hc = float(locate_critical_mode("M3", SETTINGS)["critical_field_T"])
    minimum_check = independent_extreme_check(minimum_model, m3_hc)
    maximum_check = independent_extreme_check(maximum_model, m3_hc)

    typed_retained: list[dict] = []
    for row in retained:
        typed = {key: float(value) if key not in (
            "source_design", "acquisition_role", "propagation_vector_pass"
        ) else value for key, value in row.items()}
        typed["sample_id"] = int(float(row["sample_id"]))
        typed_retained.append(typed)
    e90 = [
        row for row in engineering_rows(typed_retained, m3_hc)
        if row["temperature_K"] == 0.09
    ]
    central_count = sum(
        row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
        for row in e90
    )

    checks = {
        "total_row_count": len(rows) == summary["design"]["total_exact_count"] == 823,
        "adaptive_row_count": sum(int(float(row["adaptive_round"])) > 0 for row in rows) == 192,
        "retained_row_count": len(retained) == summary["retained"]["total_count"] == 69,
        "all_labels_recomputed": bool(np.array_equal(labels_from_observables, labels_in_file)),
        "initial_sources_unchanged": max_source_difference(rows, source) < 1.0e-12,
        "minimum_matches_summary": abs(float(np.min(shift)) - summary["retained"]["critical_field_shift_range_T"][0]) < 1.0e-12,
        "maximum_matches_summary": abs(float(np.max(shift)) - summary["retained"]["critical_field_shift_range_T"][1]) < 1.0e-12,
        "minimum_extreme_independent_multistart": minimum_check["field_and_mode_match"] and minimum_check["still_public_q_feasible"],
        "maximum_extreme_independent_multistart": maximum_check["field_and_mode_match"] and maximum_check["still_public_q_feasible"],
        "central_engineering_gate_recomputed": (
            len(e90) == summary["central_90mK_engineering_gate"]["scenario_count"]
            and central_count == summary["central_90mK_engineering_gate"]["both_gates_pass_count"]
            and not summary["central_90mK_engineering_gate"]["all_pass"]
        ),
        "constraints_match_declared": (
            summary["public_q_constraints"]["delta_inplane_rlu"] == list(DELTA_BAND)
            and summary["public_q_constraints"]["qz_folded_rlu"] == list(QZ_BAND)
            and summary["public_q_constraints"]["km_line_residual_max_rlu"] == KM_LINE_RESIDUAL_MAX
        ),
        "figure_outputs_exist": all(
            (FIGURES / f"interlayer_adaptive_inverse_uq.{suffix}").exists()
            for suffix in ("png", "pdf", "svg")
        ),
    }
    output = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "rows_checked": len(rows),
        "retained_rows_checked": len(retained),
        "independent_minimum_check": minimum_check,
        "independent_maximum_check": maximum_check,
        "recomputed_central_90mK": {
            "scenario_count": len(e90),
            "both_gates_pass_count": central_count,
        },
    }
    OUTPUT_PATH.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2), flush=True)
    if output["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
