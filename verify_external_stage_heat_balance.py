"""Independent checks for the YbNi1.6Sn stage heat-balance validation."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
RAW_ZIP = ROOT / "sources" / "raw" / "ybni16sn_supporting_data.zip"
PROVENANCE = ROOT / "sources" / "ybni16sn_stage_validation_provenance_20261012.json"


def close(a: float, b: float, rtol: float = 2e-10, atol: float = 2e-11) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def main() -> None:
    summary = json.loads((RESULTS / "ybni16sn_stage_heat_balance_summary.json").read_text(
        encoding="utf-8"))
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    with (RESULTS / "ybni16sn_stage_heat_balance_profiles.csv").open(
            encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    mass = float(summary["fit"]["effective_mass_g"])
    molar_mass = float(summary["model"]["molar_mass_g_per_mol_Yb"])
    k = mass / molar_mass * 1000.0
    parasitic = float(summary["fit"]["parasitic_load_uW"])
    identity = True
    for row in rows:
        load = float(row["applied_load_uW"])
        elapsed = float(row["elapsed_from_0p21K_min"])
        observed = float(row["measured_input_energy_with_fitted_parasitic_mJ"])
        predicted = float(row["predicted_material_enthalpy_mJ"])
        residual = float(row["residual_mJ"])
        identity &= close(observed, 0.06 * (load + parasitic) * elapsed)
        identity &= close(residual, observed - predicted)
        identity &= predicted >= -1e-12
        identity &= k > 0

    holdout_rows = [r for r in rows if r["role"] == "holdout"]
    residual = np.asarray([float(r["residual_mJ"]) for r in holdout_rows])
    scale = max(float(r["predicted_material_enthalpy_mJ"]) for r in rows)
    independent_holdout_nrmse = float(np.sqrt(np.mean(residual**2)) / scale)
    checks = {
        "source_archive_hash_matches_provenance":
            file_hash(RAW_ZIP) == provenance["dataset"]["archive_sha256"],
        "profile_row_count": len(rows) == 282,
        "two_loads_present": {float(r["applied_load_uW"]) for r in rows} == {2.5, 5.0},
        "energy_identities": identity,
        "independent_holdout_nrmse_matches_summary": close(
            independent_holdout_nrmse,
            float(summary["holdout_metrics"]["nrmse_vs_full_interval_enthalpy"])),
        "effective_mass_matches_approximate_public_mass": abs(mass - 15.0) <= 1.5,
        "parasitic_load_respects_public_upper_bound": 0 <= parasitic < 0.4,
        "holdout_nrmse_below_declared_gate": independent_holdout_nrmse < 0.02,
        "all_primary_gates_pass": bool(summary["all_validation_gates_pass"]),
        "figure_outputs_exist": all(
            (ROOT / "figures" / f"ybni16sn_stage_heat_balance_validation.{ext}").exists()
            for ext in ("png", "pdf", "svg")),
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "rows_checked": len(rows),
        "independent_holdout_nrmse": independent_holdout_nrmse,
    }
    (RESULTS / "ybni16sn_stage_heat_balance_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
