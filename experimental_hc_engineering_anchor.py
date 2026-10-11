"""Decision-scale engineering bridge anchored to public experimental Hc values.

The 2026 public source reports Hc=1.65 T from 65 mK neutron diffraction and
about 1.70 T from 80 mK bulk magnetization/thermodynamic anomalies.  These
experimental anchors are transported through the same two similarity mappings
as the model-only bridge.  This remains a surrogate, not a finite-T M4 solve.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bulk_to_engineering_surrogate import inverse_field, transform_field
from screen_model import MU_B_OVER_KB

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def main():
    bulk = json.loads((RESULTS / "bulk_polarized_stability.json").read_text())
    m0_hc = next(r["critical_field_T"] for r in bulk["rows"] if r["model"] == "M0" and r["boundary"] == "tin_foil")
    source_hc = (3 * 1.48 + 1.5 * .88) / (4.89 * MU_B_OVER_KB)
    experimental = {
        "neutron_diffraction_65mK": 1.65,
        "bulk_magnetization_80mK": 1.70,
    }
    profile = np.genfromtxt(RESULTS / "npj2022_xtrg_entropy_reconstructed.csv", delimiter=",", names=True, encoding="utf-8-sig")
    profile = profile[profile["anchor_N"] == 12]
    records = []
    robust = []
    for temperature in (.09, .10, .12):
        p = np.sort(profile[np.isclose(profile["temperature_K"], temperature)], order="field_T")
        bs, q = p["field_T"], p["Q_proxy_mJ_per_cm3"]
        k = int(np.argmax(q)); source_bstar = float(bs[k]); qmax = float(q[k])
        accepted = bs[q >= .95 * qmax]
        for method in ("reduced_field", "critical_offset"):
            m0_bstar = float(transform_field(source_bstar, source_hc, m0_hc, method))
            for label, target_hc in experimental.items():
                mapped = transform_field(np.array([source_bstar, accepted.min(), accepted.max()]), source_hc, target_hc, method)
                source_eval = float(inverse_field(m0_bstar, source_hc, target_hc, method))
                q_at_m0 = float(np.interp(source_eval, bs, q))
                records.append({
                    "temperature_K": temperature,
                    "experimental_anchor": label,
                    "mapping_method": method,
                    "source_XTRG_Hc_T": source_hc,
                    "M0_2026_Hc_T": m0_hc,
                    "experimental_Hc_T": target_hc,
                    "M0_selected_field_T": m0_bstar,
                    "experimental_anchor_B_star_T": float(mapped[0]),
                    "delta_B_star_vs_M0_T": float(mapped[0] - m0_bstar),
                    "W95_min_T": float(mapped[1]),
                    "W95_max_T": float(mapped[2]),
                    "Qmax_mJ_per_cm3_inherited": qmax,
                    "Q_at_M0_field_mJ_per_cm3": q_at_m0,
                    "regret_percent": float(100 * (qmax - q_at_m0) / qmax),
                })
        scenarios = [(method, label, hc) for method in ("reduced_field", "critical_offset") for label, hc in experimental.items()]
        grid = np.linspace(1.30, 1.60, 6001)
        regrets = []
        for method, _, hc in scenarios:
            source_eval = inverse_field(grid, source_hc, hc, method)
            q_eval = np.interp(source_eval, bs, q, left=np.nan, right=np.nan)
            regrets.append((qmax - q_eval) / qmax)
        regrets = np.asarray(regrets)
        valid = np.all(np.isfinite(regrets), axis=0)
        worst = np.nanmax(regrets, axis=0)
        idx = np.where(valid)[0][np.argmin(worst[valid])]
        robust.append({
            "temperature_K": temperature,
            "robust_field_T": float(grid[idx]),
            "worst_case_regret_percent": float(100 * worst[idx]),
            "scenario_order": [f"{m}/{l}" for m, l, _ in scenarios],
            "scenario_regret_percent": [float(100 * x) for x in regrets[:, idx]],
        })

    near90 = [r for r in records if r["temperature_K"] == .09]
    deltas = [r["delta_B_star_vs_M0_T"] for r in near90]
    regrets90 = [r["regret_percent"] for r in near90]
    payload = {
        "status": "experimental-Hc decision-scale surrogate; not a finite-temperature M4 solution",
        "public_evidence": [
            {"source": "arXiv:2607.23173v1 main.tex lines 80, 96, 101", "condition": "neutron diffraction, 0.065 K", "Hc_T": 1.65},
            {"source": "arXiv:2607.23173v1 main.tex lines 264 and 341", "condition": "bulk magnetization anomaly, 0.08 K", "Hc_T": 1.70},
            {"source": "arXiv:2607.23173v1 main.tex lines 293 and 341", "condition": "full-model mean field", "Hc_T": 1.86, "warning": "authors report about 9% overestimate versus experiment"},
        ],
        "mapping_assumptions": {"reduced_field": "preserve B/Bc", "critical_offset": "preserve B-Bc"},
        "90mK_delta_B_star_vs_M0_T_range": [float(min(deltas)), float(max(deltas))],
        "90mK_absolute_delta_B_star_min_T": float(min(abs(x) for x in deltas)),
        "90mK_regret_of_M0_field_percent_range": [float(min(regrets90)), float(max(regrets90))],
        "predeclared_thresholds": {"absolute_delta_B_star_GO_T": .05, "regret_GO_percent": 5.0},
        "threshold_screen": {
            "all_four_scenarios_absolute_delta_B_star_above_GO": bool(min(abs(x) for x in deltas) >= .05),
            "all_four_scenarios_regret_above_GO": bool(min(regrets90) >= 5.0),
        },
        "records": records,
        "robust_fields": robust,
        "limitations": [
            "The experimental Hc anchors and source XTRG profile come from different public parameterizations.",
            "Horizontal transport does not calculate entropy redistribution or changes in Qmax.",
            "The 1.65 and 1.70 T values bracket probe and temperature dependence; they are not fit parameters.",
            "A converged finite-temperature realistic-model solver is still required for a final microscopic attribution.",
        ],
    }
    (RESULTS / "experimental_hc_engineering_anchor.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    with (RESULTS / "experimental_hc_engineering_anchor.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    labels = [f"{r['experimental_anchor'].split('_')[0]}\n{r['mapping_method'].replace('_',' ')}" for r in near90]
    axes[0].bar(labels, deltas, color="#4c78a8")
    axes[0].axhline(-.05, color="#e45756", ls="--", label="|shift| GO threshold")
    axes[0].set(ylabel=r"$B^*_{experiment}-B^*_{M0}$ (T)", title="90 mK experimental-Hc anchor")
    axes[0].tick_params(axis="x", rotation=20); axes[0].grid(axis="y", alpha=.25); axes[0].legend(frameon=False, fontsize=8)
    axes[1].bar(labels, regrets90, color="#72b7b2")
    axes[1].axhline(5, color="#e45756", ls="--", label="regret GO threshold")
    axes[1].set(ylabel="Regret at M0-selected field (%)", title="Cost of ignoring experimental Hc")
    axes[1].tick_params(axis="x", rotation=20); axes[1].grid(axis="y", alpha=.25); axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIGURES / "experimental_hc_engineering_anchor.png", dpi=220); plt.close(fig)
    print(json.dumps({k: v for k, v in payload.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
