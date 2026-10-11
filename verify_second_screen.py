import json
from pathlib import Path
R=Path("results")
for path in R.glob("*.json"):
    json.loads(path.read_text(encoding="utf-8"))
for name in ("spectral_flow_entropy_N9.json","spectral_flow_entropy_N12.json"):
    data=json.loads((R/name).read_text(encoding="utf-8"))
    for row in data["rows"]:
        if row.get("quadrature_order",3)!=3: continue
        key=("combined_low_rank_and_quadrature_error_kB_per_spin" if "combined_low_rank_and_quadrature_error_kB_per_spin" in row else "low_rank_plus_quadrature_error_vs_reference_kB_per_spin")
        assert abs(row[key]) <= row["target_delta_S_error_kB_per_spin"]
        if "sign_matches_reference" in row: assert row["sign_matches_reference"]
        else: assert (row["exact_full_delta_S_kB_per_spin"]>0)==(row["reconstructed_delta_S_kB_per_spin"]>0)
print("second-screen JSON and acceptance checks: PASS")
