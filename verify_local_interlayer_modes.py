"""Verify local-grid retained modes against independent optimizer starts."""

from __future__ import annotations

import csv
import json
from pathlib import Path

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
    with (RESULTS / "interlayer_parameter_ensemble_retained.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        retained = list(csv.DictReader(handle))
    checks = []
    rows = []
    for model in retained:
        sample_id = int(model["sample_id"])
        params = InterlayerParameters(
            j2_microev=float(model["j2_microev"]),
            ja_microev=float(model["ja_microev"]),
            jb_microev=float(model["jb_microev"]),
            delta_c=float(model["delta_c"]),
        )
        repeats = []
        for seed in (260711 + sample_id, *FIXED_SEEDS):
            result = locate_critical_mode("M4", SETTINGS, seed=seed, interlayer=params)
            mode = canonical_soft_mode(result["critical_wavevector_rlu"])
            repeat = {
                "sample_id": sample_id,
                "seed": seed,
                "critical_field_T": float(result["critical_field_T"]),
                "delta_inplane_rlu": mode["delta_inplane_rlu"],
                "qz_folded_rlu": mode["qz_folded_rlu"],
                "km_line_residual_rlu": mode["km_line_residual_rlu"],
                "optimizer_success": bool(result["optimizer_success"]),
            }
            repeats.append(repeat)
            rows.append(repeat)
        best_field = max(row["critical_field_T"] for row in repeats)
        selected = [row for row in repeats
                    if best_field - row["critical_field_T"] <= FIELD_CLUSTER_TOLERANCE_T]
        best = selected[0]
        public_q_pass = (
            DELTA_BAND[0] <= best["delta_inplane_rlu"] <= DELTA_BAND[1]
            and QZ_BAND[0] <= best["qz_folded_rlu"] <= QZ_BAND[1]
            and best["km_line_residual_rlu"] <= KM_LINE_RESIDUAL_MAX
        )
        checks.append({
            "sample_id": sample_id,
            "selected_global_field_T": best_field,
            "selected_mode_repeat_count": len(selected),
            "selected_delta_inplane_rlu": best["delta_inplane_rlu"],
            "selected_qz_folded_rlu": best["qz_folded_rlu"],
            "selected_public_q_pass": public_q_pass,
            "verification_pass": len(selected) >= 2
                                 and all(row["optimizer_success"] for row in repeats),
        })

    with (RESULTS / "interlayer_local_mode_repeats.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "PASS" if all(row["verification_pass"] for row in checks) else "FAIL",
        "retained_model_count": len(retained),
        "models_reclassified_outside_public_q_filter": [
            row["sample_id"] for row in checks if not row["selected_public_q_pass"]
        ],
        "model_checks": checks,
    }
    (RESULTS / "interlayer_local_mode_verification.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if summary["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
