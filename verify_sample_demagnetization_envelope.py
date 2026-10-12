"""Independent invariants for the NBCP sample-demagnetization envelope."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def main():
    data = json.loads((RESULTS / "sample_demagnetization_envelope_summary.json").read_text(encoding="utf-8"))
    errors = []
    plate = data["public_B_parallel_c_specific_heat_sample"]
    if not (0.10 < plate["saturation_correction_bound_T"] < 0.16):
        errors.append("public plate correction is outside the expected 0.10-0.16 T scale")
    aspects = data["shape_acceptance_requirements"]
    if not (aspects["0.050_T"]["minimum_aspect_c_over_inplane_radius"]
            < aspects["0.020_T"]["minimum_aspect_c_over_inplane_radius"]
            < aspects["0.010_T"]["minimum_aspect_c_over_inplane_radius"]):
        errors.append("aspect requirements are not monotonic")
    if data["sphere_saturation_correction_T"] <= 0.05:
        errors.append("sphere correction unexpectedly below the prior 0.05 T gate")
    if data["unknown_public_1p65g_column"]["dimensions_publicly_available"]:
        errors.append("unknown sample-column dimensions were treated as public")
    report = {"status": "PASS" if not errors else "FAIL", "errors": errors}
    (RESULTS / "sample_demagnetization_envelope_verification.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
