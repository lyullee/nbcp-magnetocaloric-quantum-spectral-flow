"""Compile the N=9 state-preparation audits at 1.40, 1.45 and 1.50 T."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def load(field):
    return json.loads((RESULTS /
                       f"spectral_flow_stateprep_N9_B{field:.2f}T.json").read_text(
                           encoding="utf-8"))


def main():
    fields = (1.40, 1.45, 1.50)
    rows = []
    summaries = []
    for field in fields:
        data = load(field)
        representatives = [rep for node in data["nodes"]
                           for rep in node["representatives"]]
        thermo = data["thermodynamics"]
        rows.extend({"field_T": field, **row} for row in thermo)
        summaries.append({
            "field_T": field,
            "low_energy_levels": int(data["condition"]["low_energy_levels"]),
            "distinct_state_circuits": len(representatives),
            "minimum_qiskit_fidelity": min(
                rep["qiskit_fidelity"] for rep in representatives),
            "maximum_abs_V_error_K": max(
                abs(rep["V_expectation_error_K"]) for rep in representatives),
            "maximum_all_to_all_CX": max(
                rep["all_to_all"]["cx_count"] for rep in representatives),
            "maximum_fake_sherbrooke_ECR": max(
                rep["fake_sherbrooke"]["ecr_count"]
                for rep in representatives),
            "maximum_absolute_total_delta_S_error_kB_per_spin": max(
                abs(row["total_error_kB_per_spin"]) for row in thermo),
            "maximum_error_to_target_ratio": max(
                abs(row["total_error_kB_per_spin"]) /
                row["target_error_kB_per_spin"] for row in thermo),
            "all_targets_pass": all(row["passes_target"] for row in thermo),
            "all_signs_match": all(row["sign_matches_exact"] for row in thermo),
        })

    payload = {
        "purpose": "N=9 executable state-preparation field scan across both entropy-change signs and the zero crossing",
        "fields_T": list(fields),
        "temperatures_K": [0.09, 0.10, 0.12],
        "field_summaries": summaries,
        "thermodynamic_rows": rows,
        "status": "PASS" if all(row["all_targets_pass"] and
                                 row["all_signs_match"]
                                 for row in summaries) else "FAIL",
        "interpretation": [
            "The same physics-informed circuit family passes at negative, near-zero and positive entropy-change fields.",
            "The 1.40 T rank-7 states include a harder excited doublet with lower full-state fidelity, but the V observable and final entropy remain within the predeclared tolerance.",
            "This validates circuit expressibility and observable bias using exact classical targets; it is not a scalable eigensolver or a noisy-device accuracy claim.",
        ],
    }
    (RESULTS / "spectral_flow_stateprep_N9_field_scan.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    colors = {1.40: "#4c78a8", 1.45: "#f58518", 1.50: "#54a24b"}
    for field in fields:
        chosen = [row for row in rows if row["field_T"] == field]
        t = np.asarray([row["temperature_K"] * 1000 for row in chosen])
        exact = np.asarray([row["exact_full_delta_S_kB_per_spin"]
                            for row in chosen])
        circuit = np.asarray([row["circuit_prepared_delta_S_kB_per_spin"]
                              for row in chosen])
        axes[0].plot(t, exact, "o-", color=colors[field],
                     label=f"{field:.2f} T exact")
        axes[0].plot(t, circuit, "x--", color=colors[field], alpha=0.85,
                     label=f"{field:.2f} T circuit")
        ratio = [abs(row["total_error_kB_per_spin"]) /
                 row["target_error_kB_per_spin"] for row in chosen]
        axes[1].plot(t, ratio, "o-", color=colors[field],
                     label=f"{field:.2f} T")
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set(xlabel="Temperature (mK)",
                ylabel=r"$\Delta S$ ($k_B$/spin)",
                title="A. Exact and circuit-prepared entropy")
    axes[0].legend(frameon=True, fontsize=7, ncol=2)
    axes[1].axhline(1, color="#c44e52", linestyle="--",
                    label="acceptance limit")
    axes[1].set(xlabel="Temperature (mK)",
                ylabel="|total error| / target",
                title="B. Predeclared accuracy test")
    axes[1].legend(frameon=True, fontsize=8)

    x = np.arange(len(fields))
    fidelity_loss = [1 - row["minimum_qiskit_fidelity"]
                     for row in summaries]
    axes[2].bar(x, fidelity_loss, color=[colors[field] for field in fields])
    axes[2].set_yscale("log")
    axes[2].set_xticks(x, [f"{field:.2f}" for field in fields])
    axes[2].set(xlabel="Field (T)", ylabel="1 - minimum fidelity",
                title="C. Hardest prepared state")
    fig.suptitle("N=9 spectral-flow state preparation across field",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    figure = ROOT / "figures" / "spectral_flow_stateprep_N9_field_scan.png"
    fig.savefig(figure, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "summary": str(RESULTS /
                       "spectral_flow_stateprep_N9_field_scan.json"),
        "figure": str(figure),
        "status": payload["status"],
        "field_summaries": summaries,
    }, indent=2))


if __name__ == "__main__":
    main()
