"""Verify expanded-Sobol retained modes against optimizer seed changes."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from bulk_polarized_stability import EwaldSettings, InterlayerParameters, locate_critical_mode
from screen_interlayer_parameter_ensemble import (
    DELTA_BAND,
    KM_LINE_RESIDUAL_MAX,
    QZ_BAND,
    canonical_soft_mode,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
SETTINGS = EwaldSettings(boundary="tin_foil")
FIXED_SEEDS = (271, 2718, 27182, 271828)
FIELD_CLUSTER_TOLERANCE_T = 2e-7


def main() -> None:
    with (RESULTS / "interlayer_expanded_sobol_retained.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        retained = list(csv.DictReader(handle))

    rows: list[dict] = []
    model_checks: list[dict] = []
    for model in retained:
        params = InterlayerParameters(
            j2_microev=float(model["j2_microev"]),
            ja_microev=float(model["ja_microev"]),
            jb_microev=float(model["jb_microev"]),
            delta_c=float(model["delta_c"]),
        )
        seeds = (261012 + int(model["sample_id"]), *FIXED_SEEDS)
        repeats = []
        for seed in seeds:
            result = locate_critical_mode("M4", SETTINGS, seed=seed, interlayer=params)
            mode = canonical_soft_mode(result["critical_wavevector_rlu"])
            row = {
                "sample_id": int(model["sample_id"]),
                "seed": seed,
                "critical_field_T": float(result["critical_field_T"]),
                "delta_inplane_rlu": mode["delta_inplane_rlu"],
                "qz_folded_rlu": mode["qz_folded_rlu"],
                "km_line_residual_rlu": mode["km_line_residual_rlu"],
                "optimizer_success": bool(result["optimizer_success"]),
            }
            repeats.append(row)
            rows.append(row)
        best_field = max(row["critical_field_T"] for row in repeats)
        selected = [row for row in repeats
                    if best_field - row["critical_field_T"] <= FIELD_CLUSTER_TOLERANCE_T]
        best = selected[0]
        public_q_pass = (
            DELTA_BAND[0] <= best["delta_inplane_rlu"] <= DELTA_BAND[1]
            and QZ_BAND[0] <= best["qz_folded_rlu"] <= QZ_BAND[1]
            and best["km_line_residual_rlu"] <= KM_LINE_RESIDUAL_MAX
        )
        model_checks.append({
            "sample_id": int(model["sample_id"]),
            "critical_field_full_spread_T": float(np.ptp(
                [row["critical_field_T"] for row in repeats]
            )),
            "selected_global_field_T": best_field,
            "selected_mode_repeat_count": len(selected),
            "selected_delta_inplane_rlu": best["delta_inplane_rlu"],
            "selected_qz_folded_rlu": best["qz_folded_rlu"],
            "selected_public_q_pass": public_q_pass,
            "original_single_start_public_q_pass": True,
            "all_optimizers_successful": all(row["optimizer_success"] for row in repeats),
            "verification_pass": len(selected) >= 2
                                 and all(row["optimizer_success"] for row in repeats),
        })

    with (RESULTS / "interlayer_expanded_sobol_mode_repeats.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "status": "PASS" if all(row["verification_pass"] for row in model_checks) else "FAIL",
        "retained_model_count": len(retained),
        "fixed_seeds_plus_original_point_seed": list(FIXED_SEEDS),
        "field_cluster_tolerance_T": FIELD_CLUSTER_TOLERANCE_T,
        "models_reclassified_outside_public_q_filter": [
            row["sample_id"] for row in model_checks
            if not row["selected_public_q_pass"]
        ],
        "model_checks": model_checks,
    }
    (RESULTS / "interlayer_expanded_sobol_mode_verification.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if summary["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
