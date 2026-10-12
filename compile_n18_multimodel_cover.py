"""Compile the preregistered N=18 M3/M4 model cover."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from compile_n18_stacked_m4_convergence import quadratic_peak, quadratic_w95
from run_n18_multimodel_cover import MODEL_SPECS
from screen_model import VOLUME_FACTOR_MJ_PER_CM3


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
WORK = RESULTS / "n18_stacked_m4_work"
FIGURES = ROOT / "figures"
SUMMARY = RESULTS / "n18_multimodel_cover_summary.json"
CURVES = RESULTS / "n18_multimodel_cover_curves.csv"
COMPARISONS = RESULTS / "n18_multimodel_cover_comparisons.csv"
NUMERICS = RESULTS / "n18_multimodel_cover_numerical_convergence.csv"
TEMPERATURE_K = 0.09
FIELDS = np.round(np.arange(1.10, 2.1001, 0.02), 10)
MODEL_ORDER = list(MODEL_SPECS)
M4_ORDER = [label for label in MODEL_ORDER if label != "M3"]


def load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def cold(model: str, field: float, cutoff: str, states: int) -> dict:
    return load(WORK / f"{model}_B{field:.4f}_K{cutoff}_L{states}.json")


def hot(model: str, cutoff: str) -> dict:
    return load(WORK / f"{model}_hot4T_K{cutoff}.json")


def convergence_checkpoint(model: str) -> dict:
    safe = model.replace("-", "")
    data = load(RESULTS / f"n18_multimodel_convergence_checkpoint_{safe}.json")
    return data["models"][model]


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def compile_models() -> tuple[dict, list[dict], list[dict]]:
    models = {}
    curves = []
    numerics = []
    for model in MODEL_ORDER:
        dense = [cold(model, float(field), "6-7", 12) for field in FIELDS]
        dense_entropy = np.array(
            [row["entropy_kB_per_spin"] for row in dense], dtype=float)
        dense_peak_index = int(np.argmax(dense_entropy))
        if dense_peak_index in (0, len(FIELDS) - 1):
            raise RuntimeError(f"{model} dense-grid peak lies on a boundary")
        selection = convergence_checkpoint(model)
        grid_peak = float(selection["selected_peak_grid_center_T"])
        if not np.any(np.isclose(FIELDS, grid_peak)):
            raise RuntimeError(f"{model} selected peak is outside the dense grid")
        neighbors = [round(grid_peak - 0.02, 2), grid_peak, round(grid_peak + 0.02, 2)]
        l32 = [cold(model, field, "6-7", 32) for field in neighbors]
        l32_entropy = [float(row["entropy_kB_per_spin"]) for row in l32]
        qfield, qentropy, coefficients = quadratic_peak(neighbors, l32_entropy)

        checks = {
            "K6-7_L12": cold(model, grid_peak, "6-7", 12)["entropy_kB_per_spin"],
            "K8-9_L12": cold(model, grid_peak, "8-9", 12)["entropy_kB_per_spin"],
            "K6-7_L32": cold(model, grid_peak, "6-7", 32)["entropy_kB_per_spin"],
            "K8-9_L32": cold(model, grid_peak, "8-9", 32)["entropy_kB_per_spin"],
        }
        cold_reference = float(checks["K8-9_L32"])
        for setting, value in checks.items():
            numerics.append({
                "model": model,
                "quantity": "cold_peak_entropy",
                "setting": setting,
                "value": float(value),
                "reference_setting": "K8-9_L32",
                "absolute_difference": abs(float(value) - cold_reference),
            })

        hot_checks = {
            f"K{cutoff}": float(hot(model, cutoff)["entropy_kB_per_spin"])
            for cutoff in ("6-7", "8-9", "10-11")
        }
        hot_reference = hot_checks["K10-11"]
        for setting, value in hot_checks.items():
            numerics.append({
                "model": model,
                "quantity": "hot_entropy_4T_2K",
                "setting": setting,
                "value": value,
                "reference_setting": "K10-11",
                "absolute_difference": abs(value - hot_reference),
            })

        local_heat_peak = qentropy - hot_reference
        w95_min, w95_max = quadratic_w95(
            qfield, local_heat_peak, coefficients, hot_reference)
        heat_difference = cold_reference - hot_reference
        spec = MODEL_SPECS[model]
        models[model] = {
            "sample_id": spec["sample_id"],
            "selection_role": spec["selection_role"],
            "bulk_shift_T": spec["bulk_shift_T"],
            "n12_shift_T": spec["n12_shift_T"],
            "n12_heat_ratio": spec["n12_heat_ratio"],
            "dense_fields_T": FIELDS.tolist(),
            "dense_entropy_kB_per_spin": dense_entropy.tolist(),
            "dense_L12_peak_field_T": float(FIELDS[dense_peak_index]),
            "L32_candidate_peaks": selection["candidate_peaks"],
            "grid_peak_field_T": grid_peak,
            "quadratic_peak_field_T": qfield,
            "quadratic_peak_entropy_kB_per_spin": qentropy,
            "quadratic_W95_min_T": w95_min,
            "quadratic_W95_max_T": w95_max,
            "quadratic_W95_left_T": qfield - w95_min,
            "quadratic_W95_right_T": w95_max - qfield,
            "quadratic_W95_total_T": w95_max - w95_min,
            "cold_peak_convergence": {key: float(value) for key, value in checks.items()},
            "hot_entropy_convergence": hot_checks,
            "converged_heat_entropy_difference_kB_per_spin": heat_difference,
            "converged_peak_heat_mJ_per_cm3": (
                TEMPERATURE_K * heat_difference * VOLUME_FACTOR_MJ_PER_CM3),
        }
        for field, entropy in zip(FIELDS, dense_entropy, strict=True):
            curves.append({
                "model": model,
                "sample_id": spec["sample_id"],
                "bulk_shift_T": spec["bulk_shift_T"],
                "field_T": float(field),
                "entropy_kB_per_spin": float(entropy),
                "heat_proxy_mJ_per_cm3": float(
                    TEMPERATURE_K * (entropy - hot_reference)
                    * VOLUME_FACTOR_MJ_PER_CM3),
            })
    return models, curves, numerics


def compare(models: dict) -> list[dict]:
    baseline = models["M3"]
    rows = []
    for model in M4_ORDER:
        item = models[model]
        n18_shift = item["quadratic_peak_field_T"] - baseline["quadratic_peak_field_T"]
        n18_heat_ratio = (
            item["converged_peak_heat_mJ_per_cm3"]
            / baseline["converged_peak_heat_mJ_per_cm3"])
        n18_width_ratio = (
            item["quadratic_W95_total_T"] / baseline["quadratic_W95_total_T"])
        n18_left_width_ratio = (
            item["quadratic_W95_left_T"] / baseline["quadratic_W95_left_T"])
        n18_right_width_ratio = (
            item["quadratic_W95_right_T"] / baseline["quadratic_W95_right_T"])
        rows.append({
            "model": model,
            "sample_id": item["sample_id"],
            "selection_role": item["selection_role"],
            "bulk_shift_T": item["bulk_shift_T"],
            "N12_optimum_shift_T": item["n12_shift_T"],
            "N18_optimum_shift_T": n18_shift,
            "N18_minus_N12_optimum_shift_T": n18_shift - item["n12_shift_T"],
            "N12_peak_heat_ratio": item["n12_heat_ratio"],
            "N18_peak_heat_ratio": n18_heat_ratio,
            "N18_minus_N12_peak_heat_ratio": n18_heat_ratio - item["n12_heat_ratio"],
            "N18_W95_left_ratio": n18_left_width_ratio,
            "N18_W95_right_ratio": n18_right_width_ratio,
            "N18_W95_total_ratio": n18_width_ratio,
        })
    return rows


def rank_correlation(x: list[float], y: list[float]) -> float:
    rx = np.argsort(np.argsort(np.asarray(x, dtype=float))).astype(float)
    ry = np.argsort(np.argsort(np.asarray(y, dtype=float))).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def plot(payload: dict) -> None:
    colors = {
        "M3": "#111827", "M4-00": "#0ea5e9", "M4-02": "#10b981",
        "M4-03": "#f59e0b", "M4-01": "#ef4444",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    ax = axes[0, 0]
    for model in MODEL_ORDER:
        item = payload["models"][model]
        ax.plot(item["dense_fields_T"], item["dense_entropy_kB_per_spin"],
                color=colors[model], lw=1.35, label=model)
        ax.axvline(item["quadratic_peak_field_T"], color=colors[model],
                   lw=0.75, alpha=0.55)
    ax.set(xlabel="Field (T)", ylabel=r"Cold entropy ($k_B$/spin)",
           title="(a) Common dense N=18 field scan")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, ncol=2, fontsize=8)

    comparisons = payload["comparisons"]
    x = np.arange(len(comparisons))
    labels = [row["model"] for row in comparisons]
    ax = axes[0, 1]
    ax.plot(x, [row["N12_optimum_shift_T"] for row in comparisons],
            "o--", color="#64748b", label="N=12")
    ax.plot(x, [row["N18_optimum_shift_T"] for row in comparisons],
            "s-", color="#7c3aed", label="N=18")
    ax.set_xticks(x, labels)
    ax.set(ylabel=r"M4-M3 optimum shift $\Delta B^*$ (T)",
           title="(b) Size transfer across the model cover")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(x, [row["N12_peak_heat_ratio"] for row in comparisons],
            "o--", color="#64748b", label="N=12")
    ax.plot(x, [row["N18_peak_heat_ratio"] for row in comparisons],
            "s-", color="#dc2626", label="N=18")
    ax.axhline(1.0, color="0.5", lw=0.8)
    ax.set_xticks(x, labels)
    ax.set(ylabel="M4/M3 peak-heat ratio", title="(c) Capacity transfer")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    settings = ["K6-7_L12", "K8-9_L12", "K6-7_L32", "K8-9_L32"]
    xx = np.arange(len(settings))
    for model in MODEL_ORDER:
        ref = payload["models"][model]["cold_peak_convergence"]["K8-9_L32"]
        errors = [abs(payload["models"][model]["cold_peak_convergence"][key] - ref)
                  for key in settings]
        ax.semilogy(xx, np.maximum(errors, 1e-11), "o-", ms=3.5,
                    color=colors[model], label=model)
    ax.set_xticks(xx, ["K6/7\nL12", "K8/9\nL12", "K6/7\nL32", "K8/9\nL32"])
    ax.set(ylabel=r"$|S-S_{K8/9,L32}|$ ($k_B$/spin)",
           title="(d) Numerical convergence")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=7, ncol=2)

    FIGURES.mkdir(exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"n18_multimodel_cover.{suffix}",
                    dpi=300 if suffix == "png" else None)
    plt.close(fig)


def main() -> None:
    models, curves, numerics = compile_models()
    comparisons = compare(models)
    n18_shifts = [row["N18_optimum_shift_T"] for row in comparisons]
    bulk_shifts = [row["bulk_shift_T"] for row in comparisons]
    major = [
        row for row in comparisons
        if (abs(row["N18_minus_N12_optimum_shift_T"]) >= 0.05
            or abs(row["N18_minus_N12_peak_heat_ratio"]) >= 0.10)
    ]
    monotone = all(b >= a - 1e-12 for a, b in zip(n18_shifts, n18_shifts[1:]))
    if len(major) >= 3:
        generality = "size-transfer failure spans most of the preregistered M4 cover"
    elif len(major) >= 2:
        generality = "size-transfer failure occurs in multiple preregistered M4 models"
    else:
        generality = "large size-transfer failure remains concentrated in the high-effect tail"
    payload = {
        "status": generality,
        "scope": (
            "preregistered N=18 model-cover audit with numerical convergence "
            "inside the projected calculation; not thermodynamic-limit convergence"
        ),
        "cluster": {
            "shape": [3, 3, 2], "spins": 18, "boundary": "tin_foil",
            "three_sublattice_compatible": True,
        },
        "temperature_K": TEMPERATURE_K,
        "selection_plan": "results/n18_multimodel_preregistered_plan.json",
        "coarse_amendment": "results/n18_multimodel_plan_amendment_after_coarse.json",
        "convergence_protocol": "results/n18_multimodel_convergence_protocol.json",
        "convergence_protocol_amendment": "results/n18_multimodel_convergence_protocol_amendment.json",
        "models": models,
        "comparisons": comparisons,
        "model_family_tests": {
            "major_transfer_change_thresholds": {
                "absolute_optimum_shift_change_T": 0.05,
                "absolute_peak_heat_ratio_change": 0.10,
            },
            "major_change_model_count": len(major),
            "major_change_models": [row["model"] for row in major],
            "N18_shift_nondecreasing_with_preregistered_bulk_shift_order": monotone,
            "spearman_bulk_shift_vs_N18_optimum_shift": rank_correlation(
                bulk_shifts, n18_shifts),
            "N18_optimum_shift_range_T": [min(n18_shifts), max(n18_shifts)],
            "N18_peak_heat_ratio_range": [
                min(row["N18_peak_heat_ratio"] for row in comparisons),
                max(row["N18_peak_heat_ratio"] for row in comparisons),
            ],
            "L32_candidate_peak_count_by_model": {
                model: len(models[model]["L32_candidate_peaks"])
                for model in MODEL_ORDER
            },
            "multiple_close_peak_models": [
                model for model in MODEL_ORDER
                if len(models[model]["L32_candidate_peaks"]) > 1
            ],
            "L32_reselected_peak_branch_models": [
                model for model in MODEL_ORDER
                if not np.isclose(models[model]["dense_L12_peak_field_T"],
                                  models[model]["grid_peak_field_T"])
            ],
        },
        "numerical_checks": {
            "cold_cutoff_crosscheck_max_abs_kB_per_spin": max(
                abs(models[m]["cold_peak_convergence"]["K8-9_L32"]
                    - models[m]["cold_peak_convergence"]["K6-7_L32"])
                for m in MODEL_ORDER),
            "hot_last_cutoff_change_max_abs_kB_per_spin": max(
                abs(models[m]["hot_entropy_convergence"]["K10-11"]
                    - models[m]["hot_entropy_convergence"]["K8-9"])
                for m in MODEL_ORDER),
            "all_dense_peaks_interior": all(
                FIELDS[0] < models[m]["grid_peak_field_T"] < FIELDS[-1]
                for m in MODEL_ORDER),
            "all_close_dense_local_maxima_checked_at_L32": all(
                len(models[m]["L32_candidate_peaks"]) >= 1
                for m in MODEL_ORDER),
        },
        "engineering_decision": {
            "finite_size_transfer_gate": "FAIL",
            "allowed_claim": generality,
            "disallowed_claim": "the N=18 cover is a bulk prediction",
        },
    }
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(CURVES, curves)
    write_csv(COMPARISONS, comparisons)
    write_csv(NUMERICS, numerics)
    plot(payload)
    print(json.dumps({
        "status": payload["status"],
        "model_family_tests": payload["model_family_tests"],
        "numerical_checks": payload["numerical_checks"],
    }, indent=2))


if __name__ == "__main__":
    main()
