"""Targeted N=18 correction at the two XTRG-implied bulk field centers."""
from __future__ import annotations

import json
from pathlib import Path
import time

from ftlm import run_ftlm
from material_models import model_operators
from screen_model import lattice


ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "results" / "n18_bulk_window_compact_raw"
FIELDS = [1.50, 1.55, 1.60]
MODELS = ["M0", "M3"]
TEMPERATURES = [0.09, 0.10, 0.12, 0.15, 0.20, 0.30]
SETTINGS = {"random_vectors": 2, "krylov_dim": 30, "exact_low_states": 8, "seed": 20261011}


def calculate(graph, model, field):
    path = CACHE / f"{model}_B{field:.3f}T.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    start = time.perf_counter()
    result = run_ftlm(model_operators(graph, model, field), TEMPERATURES, 18, **SETTINGS)
    result.update(
        {
            "N": 18,
            "model": model,
            "field_T": field,
            "boundary": "tin_foil",
            "cluster_translations": [[4, 1], [-2, 4]],
            "wall_seconds_including_operator_build": time.perf_counter() - start,
        }
    )
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"saved {model} B={field:.3f} T, wall={result['wall_seconds_including_operator_build']:.1f} s", flush=True)
    return result


def low_spectrum_entropy(data, temperature):
    import numpy as np
    e = np.concatenate([np.asarray(x) for x in data["low_energies_K_by_sector"]])
    e = e - e.min()
    w = np.exp(-e / temperature)
    return float((np.log(w.sum()) + np.dot(w, e) / w.sum() / temperature) / 18)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    graph = lattice((4, 1), (-2, 4))
    data = {(m, b): calculate(graph, m, b) for m in MODELS for b in FIELDS}
    rows = []
    for ti, t in enumerate(TEMPERATURES):
        for b in FIELDS:
            a, c = data[("M0", b)], data[("M3", b)]
            s0, s3 = a["rows"][ti]["entropy_kB_per_spin"], c["rows"][ti]["entropy_kB_per_spin"]
            low0, low3 = low_spectrum_entropy(a, t), low_spectrum_entropy(c, t)
            rows.append(
                {
                    "N": 18,
                    "temperature_K": t,
                    "field_T": b,
                    "S_M0": s0,
                    "S_M3": s3,
                    "deltaS": s3 - s0,
                    "deltaS_low16_only": low3 - low0,
                    "deltaS_full_minus_low16": (s3 - s0) - (low3 - low0),
                }
            )
    payload = {
        "scope": "compact N=18 targeted correction in the XTRG-indicated bulk window",
        "settings": SETTINGS,
        "rows": rows,
        "warning": "pilot stochastic count; compare full FTLM with the deterministic low-16 result before interpreting",
    }
    (ROOT / "results" / "n18_bulk_window_pilot.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"rows_90mK": [r for r in rows if r["temperature_K"] == .09]}, indent=2))


if __name__ == "__main__":
    main()
