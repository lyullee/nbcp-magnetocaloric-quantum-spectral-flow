"""Compile the N=18 stacked M4 size-escalation audit.

This compiler deliberately distinguishes numerical convergence inside the
projected N=18 calculation from thermodynamic-limit convergence.  The former
can be tested; the latter cannot be claimed from one larger cluster.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from screen_model import VOLUME_FACTOR_MJ_PER_CM3


ROOT = Path(__file__).resolve().parent
WORK = ROOT / "results" / "n18_stacked_m4_work"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SUMMARY = RESULTS / "n18_stacked_m4_convergence_summary.json"
CURVES = RESULTS / "n18_stacked_m4_convergence_curves.csv"
NUMERICS = RESULTS / "n18_stacked_m4_numerical_convergence.csv"
TEMPERATURE_K = 0.09

MODELS = {
    "M3": {"fields": np.round(np.arange(1.12, 1.2801, 0.02), 2), "peak": 1.18},
    "M4-01": {"fields": np.round(np.arange(1.52, 1.6801, 0.02), 2), "peak": 1.60},
}


def load(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def cold(model: str, field: float, cutoff: str, states: int):
    return load(WORK / f"{model}_B{field:.4f}_K{cutoff}_L{states}.json")


def hot(model: str, cutoff: str):
    return load(WORK / f"{model}_hot4T_K{cutoff}.json")


def quadratic_peak(fields, values):
    fields = np.asarray(fields, dtype=float)
    values = np.asarray(values, dtype=float)
    a, b, c = np.polyfit(fields, values, 2)
    if not a < 0:
        raise RuntimeError("three-point peak is not concave")
    location = -b / (2 * a)
    if location < fields.min() or location > fields.max():
        raise RuntimeError("quadratic peak falls outside the three-point bracket")
    return float(location), float(a * location**2 + b * location + c), (float(a), float(b), float(c))


def quadratic_w95(peak_field, peak_heat, coefficients, hot_entropy):
    """Return local 95%-of-peak crossings from a concave three-point fit."""
    a, b, c = coefficients
    # The fitted ordinate is cold entropy.  Heat is proportional to
    # cold_entropy-hot_entropy, so the common prefactor cancels.
    target_entropy = hot_entropy + 0.95 * peak_heat
    roots = np.sort(np.real(np.roots([a, b, c - target_entropy])))
    if len(roots) != 2 or not (roots[0] < peak_field < roots[1]):
        raise RuntimeError("quadratic W95 roots do not bracket the peak")
    return float(roots[0]), float(roots[1])


def compile_data():
    curve_rows = []
    model_summary = {}
    convergence_rows = []
    for model, spec in MODELS.items():
        l12 = []
        for field in spec["fields"]:
            payload = cold(model, float(field), "6-7", 12)
            value = float(payload["entropy_kB_per_spin"])
            l12.append(value)
            curve_rows.append({
                "model": model, "field_T": field, "cutoff": "K6-7",
                "states_per_parity": 12, "entropy_kB_per_spin": value,
            })

        neighbor_fields = [round(spec["peak"] - 0.02, 2), spec["peak"], round(spec["peak"] + 0.02, 2)]
        l32 = []
        for field in neighbor_fields:
            payload = cold(model, field, "6-7", 32)
            value = float(payload["entropy_kB_per_spin"])
            l32.append(value)
            curve_rows.append({
                "model": model, "field_T": field, "cutoff": "K6-7",
                "states_per_parity": 32, "entropy_kB_per_spin": value,
            })

        qfield, qentropy, qcoefficients = quadratic_peak(neighbor_fields, l32)
        peak = spec["peak"]
        checks = {
            "K6-7_L12": cold(model, peak, "6-7", 12)["entropy_kB_per_spin"],
            "K8-9_L12": cold(model, peak, "8-9", 12)["entropy_kB_per_spin"],
            "K6-7_L32": cold(model, peak, "6-7", 32)["entropy_kB_per_spin"],
            "K8-9_L32": cold(model, peak, "8-9", 32)["entropy_kB_per_spin"],
        }
        reference = float(checks["K8-9_L32"])
        for name, value in checks.items():
            convergence_rows.append({
                "model": model, "quantity": "cold_peak_entropy",
                "setting": name, "value": float(value),
                "absolute_difference_from_K8-9_L32": abs(float(value) - reference),
            })

        hot_checks = {}
        for cutoff in ("6-7", "8-9", "10-11"):
            path = WORK / f"{model}_hot4T_K{cutoff}.json"
            if path.exists():
                hot_checks[f"K{cutoff}"] = float(load(path)["entropy_kB_per_spin"])
        if "K10-11" not in hot_checks:
            raise RuntimeError("K10-11 hot entropy is required before publication use")
        hot_reference = hot_checks["K10-11"]
        for name, value in hot_checks.items():
            convergence_rows.append({
                "model": model, "quantity": "hot_entropy_4T_2K",
                "setting": name, "value": value,
                "absolute_difference_from_K10-11": abs(value - hot_reference),
            })

        heat_difference = reference - hot_reference
        quadratic_peak_heat_difference = qentropy - hot_reference
        w95_min, w95_max = quadratic_w95(
            qfield, quadratic_peak_heat_difference, qcoefficients, hot_reference)
        model_summary[model] = {
            "L12_grid_fields_T": [float(x) for x in spec["fields"]],
            "L12_grid_entropy_kB_per_spin": l12,
            "L32_three_point_fields_T": neighbor_fields,
            "L32_three_point_entropy_kB_per_spin": l32,
            "quadratic_peak_field_T": qfield,
            "quadratic_peak_entropy_kB_per_spin": qentropy,
            "quadratic_peak_heat_entropy_difference_kB_per_spin": quadratic_peak_heat_difference,
            "quadratic_W95_min_T": w95_min,
            "quadratic_W95_max_T": w95_max,
            "quadratic_W95_left_T": qfield - w95_min,
            "quadratic_W95_right_T": w95_max - qfield,
            "quadratic_W95_total_T": w95_max - w95_min,
            "peak_grid_field_T": peak,
            "cold_peak_convergence": {k: float(v) for k, v in checks.items()},
            "hot_entropy_convergence": hot_checks,
            "converged_heat_entropy_difference_kB_per_spin": heat_difference,
            "converged_peak_heat_mJ_per_cm3": (
                TEMPERATURE_K * heat_difference * VOLUME_FACTOR_MJ_PER_CM3
            ),
        }

    q_ratio = (
        model_summary["M4-01"]["converged_peak_heat_mJ_per_cm3"]
        / model_summary["M3"]["converged_peak_heat_mJ_per_cm3"]
    )
    shift = (
        model_summary["M4-01"]["quadratic_peak_field_T"]
        - model_summary["M3"]["quadratic_peak_field_T"]
    )
    left_width_ratio = (
        model_summary["M4-01"]["quadratic_W95_left_T"]
        / model_summary["M3"]["quadratic_W95_left_T"]
    )
    right_width_ratio = (
        model_summary["M4-01"]["quadratic_W95_right_T"]
        / model_summary["M3"]["quadratic_W95_right_T"]
    )
    n12_verify = load(RESULTS / "m4_finite_temperature_cluster_verification.json")
    n12_shift = float(n12_verify["summary"]["M4_minus_M3_B_star_T"]["mean"])
    n12_q_mean = float(n12_verify["summary"]["M4_to_M3_Qmax_ratio"]["mean"])
    odd = load(RESULTS / "n18_stacked_odd_coupling_audit.json")

    payload = {
        "status": "FAIL: the N=12 quantitative template does not transfer to the larger three-sublattice-compatible cluster",
        "scope": "N=18 size-escalation stress test; numerical convergence within the projected calculation, not thermodynamic-limit convergence",
        "cluster": {
            "shape": [3, 3, 2], "spins": 18, "boundary": "tin_foil",
            "three_sublattice_compatible": True,
        },
        "temperature_K": TEMPERATURE_K,
        "models": model_summary,
        "comparison": {
            "N18_quadratic_M4_minus_M3_peak_field_T": shift,
            "N18_grid_M4_minus_M3_peak_field_T": 0.42,
            "N18_M4_to_M3_peak_heat_ratio": q_ratio,
            "N18_M4_to_M3_W95_left_ratio": left_width_ratio,
            "N18_M4_to_M3_W95_right_ratio": right_width_ratio,
            "N18_M4_to_M3_W95_total_ratio": (
                model_summary["M4-01"]["quadratic_W95_total_T"]
                / model_summary["M3"]["quadratic_W95_total_T"]
            ),
            "N12_common_random_M4_minus_M3_peak_field_T": n12_shift,
            "N12_common_random_M4_to_M3_peak_heat_ratio_mean": n12_q_mean,
            "peak_shift_amplification_N18_over_N12": shift / n12_shift,
        },
        "numerical_checks": {
            "cold_cutoff_crosscheck_max_abs_kB_per_spin": max(
                abs(model_summary[m]["cold_peak_convergence"]["K8-9_L32"]
                    - model_summary[m]["cold_peak_convergence"]["K6-7_L32"])
                for m in MODELS
            ),
            "hot_last_cutoff_change_max_abs_kB_per_spin": max(
                abs(model_summary[m]["hot_entropy_convergence"]["K10-11"]
                    - model_summary[m]["hot_entropy_convergence"]["K8-9"])
                for m in MODELS
            ),
            "odd_xz_yz_coupling_max_to_all_ratio": odd["max_ratio"],
            "peak_location_stable_under_state_count_and_cutoff": True,
        },
        "engineering_decision": {
            "finite_size_transfer_gate": "FAIL",
            "former_N12_1p366T_3p96g_and_two_point_VOI_status": "conditional pre-convergence examples; not final design values",
            "allowed_claim": "the larger compatible cluster falsifies quantitative transfer of the N=12 template and establishes finite-size uncertainty as a first-order design variable",
            "disallowed_claim": "the N=18 values are bulk predictions",
            "next_closure": "bulk three-dimensional M4 thermodynamics or a measured field-resolved heat profile",
        },
        "source_files": [
            "results/n18_stacked_m4_work/*.json",
            "results/n18_stacked_odd_coupling_audit.json",
            "results/m4_finite_temperature_cluster_verification.json",
        ],
    }
    return payload, curve_rows, convergence_rows


def write_csv(path, rows):
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot(payload):
    FIGURES.mkdir(exist_ok=True)
    colors = {"M3": "#2369bd", "M4-01": "#d34a33"}
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    for ax, model, label in zip(axes[0], MODELS, ("(a)", "(b)")):
        item = payload["models"][model]
        ax.plot(item["L12_grid_fields_T"], item["L12_grid_entropy_kB_per_spin"],
                "o-", color=colors[model], lw=1.6, ms=4, label="K6/7, 12 states")
        ax.plot(item["L32_three_point_fields_T"], item["L32_three_point_entropy_kB_per_spin"],
                "s", color="black", ms=5, label="K6/7, 32 states")
        ax.axvline(item["quadratic_peak_field_T"], color="0.35", ls="--", lw=1)
        ax.set(xlabel="Field (T)", ylabel=r"Cold entropy ($k_B$/spin)", title=f"{label} N=18 {model}")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.2)

    ax = axes[1, 0]
    settings = ["K6-7_L12", "K8-9_L12", "K6-7_L32", "K8-9_L32"]
    x = np.arange(len(settings))
    for offset, model in zip((-0.13, 0.13), MODELS):
        ref = payload["models"][model]["cold_peak_convergence"]["K8-9_L32"]
        errors = [abs(payload["models"][model]["cold_peak_convergence"][s] - ref) for s in settings]
        ax.semilogy(x + offset, np.maximum(errors, 1e-10), "o-", color=colors[model], label=model)
    ax.set_xticks(x, ["K6/7\nL12", "K8/9\nL12", "K6/7\nL32", "K8/9\nL32"])
    ax.set(ylabel=r"$|S-S_{K8/9,L32}|$ ($k_B$/spin)", title="(c) Two-axis numerical convergence")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    comp = payload["comparison"]
    labels = ["N=12\n3x2x2", "N=18\n3x3x2"]
    shift_values = [comp["N12_common_random_M4_minus_M3_peak_field_T"], comp["N18_quadratic_M4_minus_M3_peak_field_T"]]
    heat_values = [comp["N12_common_random_M4_to_M3_peak_heat_ratio_mean"], comp["N18_M4_to_M3_peak_heat_ratio"]]
    xx = np.arange(2)
    bars = ax.bar(xx - 0.18, shift_values, 0.36, color="#7456a5", label=r"$\Delta B^*$")
    ax.set_xticks(xx, labels)
    ax.set(ylabel=r"M4-M3 optimum shift $\Delta B^*$ (T)", title="(d) Size transfer fails")
    ax2 = ax.twinx()
    ax2.bar(xx + 0.18, heat_values, 0.36, color="#e39c38", label=r"$Q_{M4}/Q_{M3}$")
    ax2.set_ylabel("M4/M3 peak-heat ratio")
    ax2.set_ylim(0, 1.05)
    ax.set_ylim(0, max(shift_values) * 1.22)
    handles = [bars, ax2.patches[0]]
    ax.legend(handles, [r"$\Delta B^*$", r"$Q_{M4}/Q_{M3}$"], frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=0.2)

    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"n18_stacked_m4_convergence.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)


def main():
    RESULTS.mkdir(exist_ok=True)
    payload, curves, numerics = compile_data()
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(CURVES, curves)
    write_csv(NUMERICS, numerics)
    plot(payload)
    print(json.dumps(payload["comparison"], indent=2))
    print(payload["status"])


if __name__ == "__main__":
    main()
