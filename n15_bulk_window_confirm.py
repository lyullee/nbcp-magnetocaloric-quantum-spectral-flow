"""High-accuracy N=15 M3-M0 correction in the XTRG-indicated bulk window."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "results" / "n15_bulk_window_raw"
FIELDS = [1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70]
MODELS = ["M0", "M3"]
TEMPERATURES = [0.09, 0.10, 0.12, 0.15, 0.20, 0.30]
SETTINGS = {"random_vectors": 8, "krylov_dim": 50, "exact_low_states": 16, "seed": 20261011}


def calculate(graph, model, field):
    path = CACHE / f"{model}_B{field:.3f}T.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    start = time.perf_counter()
    result = run_ftlm(model_operators(graph, model, field), TEMPERATURES, 15, **SETTINGS)
    result.update(
        {
            "N": 15,
            "model": model,
            "field_T": field,
            "boundary": "tin_foil",
            "cluster_translations": [[4, 1], [1, 4]],
            "wall_seconds_including_operator_build": time.perf_counter() - start,
        }
    )
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"saved {model} B={field:.3f} T, wall={result['wall_seconds_including_operator_build']:.1f} s", flush=True)
    return result


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    graph = lattice((4, 1), (1, 4))
    data = {(m, b): calculate(graph, m, b) for m in MODELS for b in FIELDS}
    rows, comparison = [], []
    for ti, t in enumerate(TEMPERATURES):
        for b in FIELDS:
            s0 = data[("M0", b)]["rows"][ti]["entropy_kB_per_spin"]
            s3 = data[("M3", b)]["rows"][ti]["entropy_kB_per_spin"]
            rows.append({"N": 15, "temperature_K": t, "field_T": b, "S_M0": s0, "S_M3": s3, "deltaS": s3 - s0})
        q = [r for r in rows if r["temperature_K"] == t]
        comparison.append(
            {
                "temperature_K": t,
                "M0_B_T": max(q, key=lambda r: r["S_M0"])["field_T"],
                "M3_B_T": max(q, key=lambda r: r["S_M3"])["field_T"],
                "deltaS_by_field": {str(r["field_T"]): r["deltaS"] for r in q},
            }
        )
    payload = {
        "scope": "N=15 correction scan in XTRG-indicated 1.4--1.7 T bulk window",
        "reason_for_window": "reconstructed published XTRG M0 entropy peaks near 1.50 T; dimensionless parameter mapping places the current M0 peak near 1.60 T",
        "settings": SETTINGS,
        "rows": rows,
        "comparison": comparison,
    }
    (ROOT / "results" / "n15_bulk_window_confirmation.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"comparison": comparison}, indent=2))


if __name__ == "__main__":
    main()
