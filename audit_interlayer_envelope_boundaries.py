"""Audit whether the retained interlayer set terminates at grid boundaries.

This is an adaptive boundary check, not a second rectangular parameter scan.
For each model retained by the 375-point local grid, every parameter that lies
on an original grid edge is moved one additional step outward while the other
parameters are held fixed.  The calculation therefore answers the narrow but
important question: does the public-q-compatible set visibly continue beyond
the declared local envelope?
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

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
    out_of_plane_bandwidth_microev,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
SETTINGS = EwaldSettings(boundary="tin_foil")

OUTWARD_VALUES = {
    "j2_microev": {min(GRID["j2_microev"]): 0.0, max(GRID["j2_microev"]): 5.6},
    "ja_microev": {min(GRID["ja_microev"]): -4.2, max(GRID["ja_microev"]): 4.2},
    "jb_microev": {min(GRID["jb_microev"]): 1.975, max(GRID["jb_microev"]): 13.825},
    "delta_c": {min(GRID["delta_c"]): 0.02, max(GRID["delta_c"]): 0.14},
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def passes(mode: dict[str, float]) -> tuple[bool, bool, bool, bool]:
    delta_pass = DELTA_BAND[0] <= mode["delta_inplane_rlu"] <= DELTA_BAND[1]
    qz_pass = QZ_BAND[0] <= mode["qz_folded_rlu"] <= QZ_BAND[1]
    km_pass = mode["km_line_residual_rlu"] <= KM_LINE_RESIDUAL_MAX
    return delta_pass, qz_pass, km_pass, delta_pass and qz_pass and km_pass


def main() -> None:
    all_rows = read_csv(RESULTS / "interlayer_parameter_ensemble.csv")
    retained = read_csv(RESULTS / "interlayer_parameter_ensemble_retained.csv")

    condition_counts = {
        "delta_only": 0,
        "qz_only": 0,
        "km_only": 0,
        "delta_and_qz": 0,
        "delta_and_km": 0,
        "qz_and_km": 0,
        "all_three": 0,
    }
    for row in all_rows:
        delta_pass = DELTA_BAND[0] <= float(row["delta_inplane_rlu"]) <= DELTA_BAND[1]
        qz_pass = QZ_BAND[0] <= float(row["qz_folded_rlu"]) <= QZ_BAND[1]
        km_pass = float(row["km_line_residual_rlu"]) <= KM_LINE_RESIDUAL_MAX
        condition_counts["delta_only"] += int(delta_pass)
        condition_counts["qz_only"] += int(qz_pass)
        condition_counts["km_only"] += int(km_pass)
        condition_counts["delta_and_qz"] += int(delta_pass and qz_pass)
        condition_counts["delta_and_km"] += int(delta_pass and km_pass)
        condition_counts["qz_and_km"] += int(qz_pass and km_pass)
        condition_counts["all_three"] += int(delta_pass and qz_pass and km_pass)

    m3_hc = float(locate_critical_mode("M3", SETTINGS)["critical_field_T"])
    candidates: dict[tuple[float, float, float, float], dict] = {}
    keys = tuple(GRID)
    for source in retained:
        base = {key: float(source[key]) for key in keys}
        for key in keys:
            replacement = OUTWARD_VALUES[key].get(base[key])
            if replacement is None:
                continue
            values = dict(base)
            values[key] = replacement
            token = tuple(values[name] for name in keys)
            candidates.setdefault(token, {
                "source_sample_ids": [],
                "outward_parameter": key,
                **values,
            })["source_sample_ids"].append(int(source["sample_id"]))

    audited: list[dict] = []
    for index, candidate in enumerate(candidates.values(), start=1):
        parameters = InterlayerParameters(**{key: candidate[key] for key in keys})
        result = locate_critical_mode_multistart(
            "M4", SETTINGS,
            seeds=(261012 + index, 271, 27182),
            interlayer=parameters,
        )
        mode = canonical_soft_mode(result["critical_wavevector_rlu"])
        delta_pass, qz_pass, km_pass, all_pass = passes(mode)
        audited.append({
            "audit_id": index,
            "source_sample_ids": ";".join(map(str, candidate["source_sample_ids"])),
            "outward_parameter": candidate["outward_parameter"],
            **{key: candidate[key] for key in keys},
            "critical_field_T": float(result["critical_field_T"]),
            "critical_field_shift_vs_M3_T": float(result["critical_field_T"] - m3_hc),
            "delta_inplane_rlu": mode["delta_inplane_rlu"],
            "qz_folded_rlu": mode["qz_folded_rlu"],
            "km_line_residual_rlu": mode["km_line_residual_rlu"],
            "out_of_plane_bandwidth_microeV": out_of_plane_bandwidth_microev(parameters),
            "delta_pass": delta_pass,
            "qz_pass": qz_pass,
            "km_pass": km_pass,
            "all_constraints_pass": all_pass,
        })
        print(f"completed boundary audit {index}/{len(candidates)}", flush=True)

    csv_path = RESULTS / "interlayer_boundary_audit.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audited[0]))
        writer.writeheader()
        writer.writerows(audited)

    passed = [row for row in audited if row["all_constraints_pass"]]
    summary = {
        "status": "one-step outward boundary audit; not a full expanded grid",
        "original_grid_size": len(all_rows),
        "original_retained_size": len(retained),
        "constraint_attribution_counts": condition_counts,
        "outward_candidates": len(audited),
        "outward_candidates_retained": len(passed),
        "retained_outward_parameter_counts": {
            key: sum(row["all_constraints_pass"] and row["outward_parameter"] == key for row in audited)
            for key in keys
        },
        "retained_outward_shift_range_T": (
            [
                min(row["critical_field_shift_vs_M3_T"] for row in passed),
                max(row["critical_field_shift_vs_M3_T"] for row in passed),
            ]
            if passed else None
        ),
        "audit_bounds": OUTWARD_VALUES,
        "interpretation": (
            "At least one public-q-compatible model continues outside the original local envelope; "
            "the seven-model set is therefore a local sensitivity set, not an exhaustive admissible set."
            if passed else
            "No one-step outward neighbor retained all public-q constraints; the retained set closes at the tested local boundary."
        ),
        "limitations": [
            "Only one coordinate is extended at a time.",
            "A disconnected admissible region or a combination of multiple outward moves is not tested.",
            "Constraint counts describe this designed grid and are not probabilities or parameter confidence intervals.",
        ],
    }
    (RESULTS / "interlayer_boundary_audit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
