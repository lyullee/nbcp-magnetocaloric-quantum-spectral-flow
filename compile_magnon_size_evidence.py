"""Compile absolute M0 and low-energy M3-M0 size evidence at 1.50 T."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def low_entropy(energies, n, t=.09):
    e = np.sort(np.asarray(energies, float)); w = np.exp(-(e - e[0]) / t)
    return float((np.log(w.sum()) + np.dot(w, e - e[0]) / w.sum() / t) / n)


def raw_delta(n, folder):
    values = {}
    for model in ("M0", "M3"):
        data = json.loads((RESULTS / folder / f"{model}_B1.500T.json").read_text())
        energies = np.concatenate([np.asarray(block)[:8] for block in data["low_energies_K_by_sector"]])
        values[model] = low_entropy(energies, n)
    return {"N": n, "S_M0_low16": values["M0"], "S_M3_low16": values["M3"],
            "deltaS_low16": values["M3"] - values["M0"], "source": folder}


def main():
    m0 = [{"N": 18, "entropy_kB_per_spin": 0.07112174042109666,
           "source": "full parity-space low-energy-corrected FTLM"}]
    for n in (21, 24, 27):
        path = RESULTS / f"magnon_m0_N{n}_B1.500T.json"
        if path.exists():
            data = json.loads(path.read_text())
            row = next(row for row in data["rows"] if np.isclose(row["temperature_K"], .09))
            m0.append({"N": n, "entropy_kB_per_spin": row["entropy_kB_per_spin"],
                       "source": path.name, "max_q_down": data["settings"]["max_q_down"],
                       "first_omitted_sector_upper_log": row["next_sector_partition_upper_log"]})

    corrections = [raw_delta(15, "n15_bulk_window_raw"),
                   raw_delta(18, "n18_bulk_window_compact_raw")]
    for n in (21, 24, 27):
        path = RESULTS / f"magnon_M3_low_energy_N{n}_B1.500T.json"
        if path.exists():
            data = json.loads(path.read_text())
            if "runs" in data:
                run = data["runs"][-1]
                row = next(row for row in run["rows"] if np.isclose(row["temperature_K"], .09))
                max_q = run["max_q_down"]
            else:
                row = next(row for row in data["rows"] if np.isclose(row["temperature_K"], .09))
                max_q = max(block["effective_max_q_down"] for block in data["models"]["M3"]["parities"])
            corrections.append({"N": n, **row, "max_q_down": max_q, "source": path.name})

    coverage = []
    by_n_correction = {row["N"]: row for row in corrections}
    for row in m0:
        if row["N"] in by_n_correction:
            low = by_n_correction[row["N"]]["S_M0_low16"]
            coverage.append({
                "N": row["N"],
                "M0_absolute_entropy_kB_per_spin": row["entropy_kB_per_spin"],
                "M0_low16_entropy_kB_per_spin": low,
                "missing_from_low16_kB_per_spin": row["entropy_kB_per_spin"] - low,
            })

    payload = {
        "condition": {"temperature_K": .09, "field_T": 1.50, "boundary": "tin_foil"},
        "M0_absolute_size_sequence": m0,
        "M3_minus_M0_low16_size_sequence": corrections,
        "M0_low16_coverage": coverage,
        "interpretation": "Low-16 deltaS is validated at N=18, but its absolute-state coverage worsens with N; N=27 is a sign/magnitude diagnostic, not a converged bulk correction or a B*/regret result.",
    }
    (RESULTS / "magnon_size_evidence_B1.500T.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].plot([r["N"] for r in m0], [r["entropy_kB_per_spin"] for r in m0], "o-")
    axes[0].set(xlabel="Cluster size N", ylabel=r"M0 entropy ($k_B$/spin)", title="90 mK, 1.50 T absolute M0")
    axes[1].axhline(0, color="k", lw=1)
    axes[1].plot([r["N"] for r in corrections], [r["deltaS_low16"] for r in corrections], "o-")
    axes[1].set(xlabel="Cluster size N", ylabel=r"$S_{M3}-S_{M0}$ ($k_B$/spin)", title="Low-energy realistic-interaction correction")
    axes[2].plot([r["N"] for r in coverage], [r["missing_from_low16_kB_per_spin"] for r in coverage], "o-", color="#e45756")
    axes[2].set(xlabel="Cluster size N", ylabel=r"M0 entropy omitted by low 16 ($k_B$/spin)", title="Low-energy coverage warning")
    for ax in axes: ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(FIGURES / "magnon_size_evidence_B1.500T.png", dpi=220); plt.close(fig)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__": main()
