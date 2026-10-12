"""Verify identities and frozen values for the Gd2B2MoO9 public-data audit."""

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "results" / "gd2b2moo9_external_validation.json"
VERIFICATION = ROOT / "results" / "gd2b2moo9_external_validation_verification.json"

payload = json.loads(RESULT.read_text(encoding="utf-8"))
points = payload["points"]
errors = np.asarray([point["error_indirect_minus_direct_K"] for point in points])
absolute = np.asarray([point["absolute_error_K"] for point in points])
checks = {
    "twelve_points": len(points) == 12,
    "four_fields": sorted({point["field_T"] for point in points}) == [0.5, 1.0, 2.0, 5.0],
    "error_identity": all(abs(point["indirect_interpolated_delta_Tad_K"] - point["direct_delta_Tad_K"] - point["error_indirect_minus_direct_K"]) < 1e-12 for point in points),
    "absolute_error_identity": np.allclose(absolute, np.abs(errors), rtol=0, atol=1e-12),
    "mae_identity": abs(payload["mean_absolute_error_K"] - float(np.mean(np.abs(errors)))) < 1e-12,
    "rmse_identity": abs(payload["root_mean_square_error_K"] - float(np.sqrt(np.mean(errors**2)))) < 1e-12,
    "maximum_identity": abs(payload["maximum_absolute_error_K"] - float(np.max(np.abs(errors)))) < 1e-12,
    "frozen_mae": abs(payload["mean_absolute_error_K"] - 0.013801210383537404) < 1e-12,
    "frozen_maximum": abs(payload["maximum_absolute_error_K"] - 0.02763626068590963) < 1e-12,
    "declared_source_hash": payload["workbook_sha256_observed"] == payload["source"]["workbook_sha256"],
    "figures_exist": all((ROOT / "figures" / f"gd2b2moo9_external_validation.{suffix}").exists() for suffix in ["png", "pdf", "svg"]),
    "scope_not_nbcp_surrogate_validation": any("does not validate the NBCP" in sentence for sentence in payload["interpretation"]),
}
failed = [name for name, passed in checks.items() if not passed]
verification = {"status": "PASS" if not failed else "FAIL", "checks": checks, "failed": failed}
VERIFICATION.write_text(json.dumps(verification, indent=2), encoding="utf-8")
print(json.dumps(verification, indent=2))
raise SystemExit(1 if failed else 0)
