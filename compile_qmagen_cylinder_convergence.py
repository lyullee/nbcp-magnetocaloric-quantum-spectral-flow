"""Compile public-QMagen YC6-width cylinder convergence without overstating pilots."""
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def public_reference(temperature_k: float, field_t: float = 1.50) -> float:
    data = np.genfromtxt(
        RESULTS / "npj2022_xtrg_entropy_reconstructed.csv",
        delimiter=",",
        names=True,
        encoding="utf-8-sig",
    )
    q = data[(data["anchor_N"] == 12) & np.isclose(data["field_T"], field_t)]
    order = np.argsort(q["temperature_K"])
    return float(np.interp(temperature_k, q["temperature_K"][order], q["entropy_kB_per_spin"][order]))


def main() -> None:
    runs = []
    pattern = re.compile(r"qmagen_cylinder_9x6_D(\d+)\.json$")
    for path in RESULTS.glob("qmagen_cylinder_9x6_D*.json"):
        if not pattern.fullmatch(path.name):
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        rows = data["rows"]
        cold = data["cold_row"]
        ref = public_reference(cold["temperature_K"], data["field_T"])
        s = np.array([row["entropy_kB_per_spin"] for row in rows], dtype=float)
        # Rows run from hot to cold; positive differences are unphysical upturns on cooling.
        upturns = np.diff(s)
        runs.append(
            {
                "D": int(data["bond_dimension"]),
                "elapsed_seconds": float(data["elapsed_seconds"]),
                "temperature_K": float(cold["temperature_K"]),
                "entropy_kB_per_spin": float(cold["entropy_kB_per_spin"]),
                "public_reconstructed_reference_kB_per_spin": ref,
                "entropy_error_kB_per_spin": float(cold["entropy_kB_per_spin"] - ref),
                "max_entropy_upturn_on_cooling_kB_per_spin": float(max(0.0, np.max(upturns))),
                "file": path.name,
                "rows": rows,
            }
        )
    if not runs:
        raise RuntimeError("No 9x6 cylinder results found")
    runs.sort(key=lambda item: item["D"])
    for i, run in enumerate(runs):
        run["change_from_previous_D_kB_per_spin"] = (
            None if i == 0 else run["entropy_kB_per_spin"] - runs[i - 1]["entropy_kB_per_spin"]
        )

    last = runs[-1]
    converged = (
        len(runs) >= 2
        and last["D"] >= 64
        and abs(last["entropy_error_kB_per_spin"]) <= 0.005
        and abs(last["change_from_previous_D_kB_per_spin"]) <= 0.002
        and last["max_entropy_upturn_on_cooling_kB_per_spin"] <= 0.002
    )
    payload = {
        "engine": "public QMagen commit dcf17db27682d4c505155682108c612bdd11274c",
        "geometry": "9x6 OBC x / PBC y (YC6 width, N=54)",
        "condition": {"model": "Gao et al. 2022 M0", "field_T": 1.50},
        "reference": {
            "method": "reconstructed YC6x9 D=400 C/T, integrated and N=12-anchored at 1 K",
            "validation_RMSE_J_mol_K": 0.033549922615027414,
            "note": "Reference is a digitized reconstruction, not an author-supplied numeric table.",
        },
        "predeclared_convergence_check": {
            "minimum_D": 64,
            "absolute_reference_error_limit_kB_per_spin": 0.005,
            "successive_D_change_limit_kB_per_spin": 0.002,
            "entropy_upturn_limit_kB_per_spin": 0.002,
        },
        "runs": runs,
        "decision": "M0_CYLINDER_CONVERGED" if converged else "PILOT_ONLY_M0_CYLINDER_NOT_YET_CONVERGED",
    }
    (RESULTS / "qmagen_cylinder_convergence.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for run in runs:
        rows = run["rows"]
        t = np.array([row["temperature_K"] for row in rows])
        s = np.array([row["entropy_kB_per_spin"] for row in rows])
        mask = t <= 1.5
        axes[0].semilogx(t[mask], s[mask], "o-", label=f"D={run['D']}")
    axes[0].axhline(runs[-1]["public_reconstructed_reference_kB_per_spin"], color="k", ls="--", lw=1, label="public proxy at 90.3 mK")
    axes[0].set(xlabel="Temperature (K)", ylabel=r"Entropy ($k_B$/spin)", title="YC6-width M0 physicality check")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    d = np.array([run["D"] for run in runs])
    cold_s = np.array([run["entropy_kB_per_spin"] for run in runs])
    axes[1].plot(d, cold_s, "o-", label="QMagen at 90.332 mK")
    axes[1].axhline(runs[-1]["public_reconstructed_reference_kB_per_spin"], color="k", ls="--", lw=1, label="public reconstructed proxy")
    axes[1].set(xlabel="Maximum bond dimension D", ylabel=r"Entropy ($k_B$/spin)", title="Cold entropy convergence")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "qmagen_cylinder_convergence.png", dpi=220)
    plt.close(fig)
    print(json.dumps({"decision": payload["decision"], "runs": [{k: v for k, v in run.items() if k != "rows"} for run in runs]}, indent=2))


if __name__ == "__main__":
    main()
