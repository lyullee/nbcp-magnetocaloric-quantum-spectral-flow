"""Regression checks for the frozen N=9 quantum-algorithm target."""
import json
from pathlib import Path


def main():
    path = Path(__file__).resolve().parent / "results" / "quantum_local_exact_target_N9.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = [row for row in data["rows"] if row["quadrature_order"] == 5]
    assert len(rows) == 9
    assert max(abs(row["delta_F_error_K_per_spin"]) for row in rows) < 1e-12
    by_key = {(row["field_T"], row["temperature_K"]): row for row in rows}
    assert by_key[(1.40, 0.09)]["exact_delta_S_kB_per_spin"] < 0
    assert by_key[(1.45, 0.09)]["exact_delta_S_kB_per_spin"] > 0
    assert by_key[(1.50, 0.09)]["exact_delta_S_kB_per_spin"] > 0
    assert abs(by_key[(1.45, 0.09)]["exact_delta_S_kB_per_spin"] -
               0.0006669278613436191) < 1e-14
    print("quantum local exact target regression passed")


if __name__ == "__main__":
    main()
