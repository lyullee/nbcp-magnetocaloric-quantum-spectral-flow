"""Adaptive set-valued inverse audit for the M4 interlayer model.

The existing 375-point Cartesian design and 256-point expanded Sobol design
are exact evaluations but finite.  This script uses their exact observables to
choose additional points near the public-q acceptance boundary, predicted
critical-field extrema, and poorly covered parts of the declared 4-D box.

The radial-basis surrogate is an acquisition device only.  Every reported
point, feasible label, and critical-field extremum is recomputed with the
original multistart M4 calculation.  The result is an adaptive stress audit,
not a posterior distribution or a mathematical proof over a continuum.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.spatial import cKDTree
from scipy.stats import qmc

from bulk_polarized_stability import (
    EwaldSettings,
    InterlayerParameters,
    locate_critical_mode,
    locate_critical_mode_multistart,
)
from screen_interlayer_expanded_sobol import BOUNDS, PARAMETERS
from screen_interlayer_parameter_ensemble import (
    DELTA_BAND,
    GRID,
    KM_LINE_RESIDUAL_MAX,
    QZ_BAND,
    canonical_soft_mode,
    engineering_rows,
    out_of_plane_bandwidth_microev,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SETTINGS = EwaldSettings(boundary="tin_foil")

INITIAL_FILES = (
    ("local_cartesian", RESULTS / "interlayer_parameter_ensemble.csv"),
    ("expanded_sobol", RESULTS / "interlayer_expanded_sobol.csv"),
)
OBSERVABLES = (
    "delta_inplane_rlu",
    "qz_folded_rlu",
    "km_line_residual_rlu",
    "critical_field_shift_vs_M3_T",
)
ACTIVE_ROUNDS = 3
POINTS_PER_ROUND = 64
POOL_POWER = 15
POOL_SEED = 271013
OPTIMIZER_SEEDS = (271, 27182)
RBF_NEIGHBORS = 96
RBF_SMOOTHING = 1.0e-8


def read_rows(path: Path, source: str) -> list[dict]:
    rows: list[dict] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            row = {key: float(raw[key]) for key in PARAMETERS}
            row.update({key: float(raw[key]) for key in OBSERVABLES})
            row.update({
                "critical_field_T": float(raw["critical_field_T"]),
                "q_u_rlu": float(raw["q_u_rlu"]),
                "q_v_rlu": float(raw["q_v_rlu"]),
                "q_w_rlu": float(raw["q_w_rlu"]),
                "out_of_plane_bandwidth_microeV": float(
                    raw["out_of_plane_bandwidth_microeV"]
                ),
                "propagation_vector_pass": str(
                    raw["propagation_vector_pass"]
                ).lower() == "true",
                "multistart_selected_seed": int(float(
                    raw["multistart_selected_seed"]
                )),
                "multistart_min_field_T": float(
                    raw["multistart_min_field_T"]
                ),
                "multistart_field_spread_T": float(
                    raw["multistart_field_spread_T"]
                ),
                "source_design": source,
                "adaptive_round": 0,
                "acquisition_role": "initial_exact_design",
            })
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def arrays(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    x = np.array([[float(row[key]) for key in PARAMETERS] for row in rows])
    y = np.array([[float(row[key]) for key in OBSERVABLES] for row in rows])
    return x, y


def bounds_arrays() -> tuple[np.ndarray, np.ndarray]:
    lower = np.array([BOUNDS[key][0] for key in PARAMETERS], dtype=float)
    upper = np.array([BOUNDS[key][1] for key in PARAMETERS], dtype=float)
    return lower, upper


def normalize(x: np.ndarray) -> np.ndarray:
    lower, upper = bounds_arrays()
    return (np.asarray(x, dtype=float) - lower) / (upper - lower)


def constraint_violation(observables: np.ndarray) -> np.ndarray:
    """Return <=0 for points satisfying every public-q constraint."""
    obs = np.atleast_2d(np.asarray(observables, dtype=float))
    delta = obs[:, 0]
    qz = obs[:, 1]
    residual = obs[:, 2]
    delta_width = DELTA_BAND[1] - DELTA_BAND[0]
    qz_width = QZ_BAND[1] - QZ_BAND[0]
    components = np.column_stack((
        (DELTA_BAND[0] - delta) / delta_width,
        (delta - DELTA_BAND[1]) / delta_width,
        (QZ_BAND[0] - qz) / qz_width,
        (qz - QZ_BAND[1]) / qz_width,
        (residual - KM_LINE_RESIDUAL_MAX) / KM_LINE_RESIDUAL_MAX,
    ))
    return np.max(components, axis=1)


def exact_pass(obs: np.ndarray) -> np.ndarray:
    return constraint_violation(obs) <= 0.0


def fit_surrogate(rows: list[dict]) -> RBFInterpolator:
    x, y = arrays(rows)
    return RBFInterpolator(
        normalize(x), y,
        neighbors=min(RBF_NEIGHBORS, len(rows)),
        kernel="thin_plate_spline",
        smoothing=RBF_SMOOTHING,
    )


def select_diverse(
    pool_x: np.ndarray,
    evaluated_x: np.ndarray,
    ordered_groups: list[tuple[str, np.ndarray, int]],
    target: int,
) -> tuple[np.ndarray, list[str]]:
    pool_n = normalize(pool_x)
    evaluated_n = normalize(evaluated_x)
    distances, _ = cKDTree(evaluated_n).query(pool_n, k=1)
    selected: list[int] = []
    roles: list[str] = []

    def eligible(index: int) -> bool:
        if distances[index] < 0.006:
            return False
        if selected:
            separation = np.linalg.norm(pool_n[index] - pool_n[selected], axis=1)
            if float(np.min(separation)) < 0.024:
                return False
        return True

    for role, indices, quota in ordered_groups:
        added = 0
        for index in map(int, indices):
            if index in selected or not eligible(index):
                continue
            selected.append(index)
            roles.append(role)
            added += 1
            if added >= quota or len(selected) >= target:
                break
        if len(selected) >= target:
            break

    if len(selected) < target:
        for index in np.argsort(-distances):
            index = int(index)
            if index in selected or not eligible(index):
                continue
            selected.append(index)
            roles.append("space_filling_fallback")
            if len(selected) >= target:
                break

    if len(selected) != target:
        raise RuntimeError(f"selected only {len(selected)} of {target} active points")
    return np.asarray(selected, dtype=int), roles


def acquire(rows: list[dict], round_index: int) -> tuple[np.ndarray, list[str], dict]:
    model = fit_surrogate(rows)
    sampler = qmc.Sobol(
        d=len(PARAMETERS), scramble=True, seed=POOL_SEED + round_index
    )
    unit = sampler.random_base2(m=POOL_POWER)
    lower, upper = bounds_arrays()
    pool = qmc.scale(unit, lower, upper)
    predicted = model(unit)
    violation = constraint_violation(predicted)
    predicted_feasible = violation <= 0.0

    evaluated_x, _ = arrays(rows)
    nearest_distance, _ = cKDTree(normalize(evaluated_x)).query(unit, k=1)
    feasible_indices = np.flatnonzero(predicted_feasible)
    boundary_order = np.argsort(np.abs(violation) - 0.12 * nearest_distance)
    if feasible_indices.size:
        min_shift_order = feasible_indices[
            np.argsort(predicted[feasible_indices, 3])
        ]
        max_shift_order = feasible_indices[
            np.argsort(-predicted[feasible_indices, 3])
        ]
        feasible_distance_order = feasible_indices[
            np.argsort(-nearest_distance[feasible_indices])
        ]
    else:
        min_shift_order = boundary_order
        max_shift_order = boundary_order
        feasible_distance_order = boundary_order
    exploration_order = np.argsort(-nearest_distance)

    groups = [
        ("predicted_feasible_min_shift", min_shift_order, 12),
        ("predicted_feasible_max_shift", max_shift_order, 12),
        ("public_q_boundary", boundary_order, 24),
        ("predicted_feasible_coverage", feasible_distance_order, 8),
        ("space_filling", exploration_order, 8),
    ]
    selected, roles = select_diverse(
        pool, evaluated_x, groups, POINTS_PER_ROUND
    )
    diagnostics = {
        "pool_count": int(len(pool)),
        "predicted_feasible_pool_count": int(np.sum(predicted_feasible)),
        "selected_nearest_existing_distance_range_normalized": [
            float(np.min(nearest_distance[selected])),
            float(np.max(nearest_distance[selected])),
        ],
        "selected_predicted_violation_range": [
            float(np.min(violation[selected])),
            float(np.max(violation[selected])),
        ],
    }
    return pool[selected], roles, diagnostics


def evaluate_point(
    values: np.ndarray,
    role: str,
    round_index: int,
    sequence_index: int,
    m3_hc: float,
) -> dict:
    params = dict(zip(PARAMETERS, map(float, values), strict=True))
    interlayer = InterlayerParameters(**params)
    result = locate_critical_mode_multistart(
        "M4", SETTINGS,
        seeds=(POOL_SEED + 1000 * round_index + sequence_index, *OPTIMIZER_SEEDS),
        interlayer=interlayer,
    )
    mode = canonical_soft_mode(result["critical_wavevector_rlu"])
    obs = np.array([[
        mode["delta_inplane_rlu"],
        mode["qz_folded_rlu"],
        mode["km_line_residual_rlu"],
        float(result["critical_field_T"] - m3_hc),
    ]])
    return {
        **params,
        "critical_field_T": float(result["critical_field_T"]),
        "critical_field_shift_vs_M3_T": float(result["critical_field_T"] - m3_hc),
        "q_u_rlu": float(result["critical_wavevector_rlu"][0]),
        "q_v_rlu": float(result["critical_wavevector_rlu"][1]),
        "q_w_rlu": float(result["critical_wavevector_rlu"][2]),
        "delta_inplane_rlu": mode["delta_inplane_rlu"],
        "qz_folded_rlu": mode["qz_folded_rlu"],
        "km_line_residual_rlu": mode["km_line_residual_rlu"],
        "out_of_plane_bandwidth_microeV": out_of_plane_bandwidth_microev(interlayer),
        "propagation_vector_pass": bool(exact_pass(obs)[0]),
        "multistart_selected_seed": int(result["multistart"]["selected_seed"]),
        "multistart_min_field_T": float(
            result["multistart"]["critical_field_range_T"][0]
        ),
        "multistart_field_spread_T": float(
            result["multistart"]["critical_field_range_T"][1]
            - result["multistart"]["critical_field_range_T"][0]
        ),
        "source_design": "adaptive_exact",
        "adaptive_round": round_index,
        "acquisition_role": role,
    }


def range_of(rows: list[dict], key: str) -> list[float] | None:
    if not rows:
        return None
    values = [float(row[key]) for row in rows]
    return [min(values), max(values)]


def inside_original_box(row: dict) -> bool:
    return all(
        min(GRID[key]) <= float(row[key]) <= max(GRID[key])
        for key in PARAMETERS
    )


def cross_validation(initial_rows: list[dict]) -> dict:
    """Deterministic diagnostic; exact calculations remain authoritative."""
    retained_indices = [
        i for i, row in enumerate(initial_rows) if row["propagation_vector_pass"]
    ]
    nonretained_indices = [
        i for i, row in enumerate(initial_rows) if not row["propagation_vector_pass"]
    ]
    holdout_indices = retained_indices[::2] + nonretained_indices[::11]
    holdout = set(holdout_indices)
    train = [row for i, row in enumerate(initial_rows) if i not in holdout]
    test = [row for i, row in enumerate(initial_rows) if i in holdout]
    model = fit_surrogate(train)
    x_test, y_test = arrays(test)
    predicted = model(normalize(x_test))
    truth = exact_pass(y_test)
    pred_label = exact_pass(predicted)
    mae = np.mean(np.abs(predicted - y_test), axis=0)
    return {
        "role": "acquisition-surrogate diagnostic only",
        "train_count": len(train),
        "holdout_count": len(test),
        "holdout_feasible_count": int(np.sum(truth)),
        "observable_mae": {
            key: float(value) for key, value in zip(OBSERVABLES, mae, strict=True)
        },
        "classification": {
            "true_positive": int(np.sum(truth & pred_label)),
            "false_positive": int(np.sum(~truth & pred_label)),
            "false_negative": int(np.sum(truth & ~pred_label)),
            "true_negative": int(np.sum(~truth & ~pred_label)),
        },
        "interpretation": (
            "This split diagnoses acquisition quality only. Every final label and "
            "extremum comes from an exact multistart M4 evaluation."
        ),
    }


def make_figure(rows: list[dict], convergence: list[dict]) -> None:
    x, _ = arrays(rows)
    xn = normalize(x)
    centered = xn - np.mean(xn, axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    projected = centered @ vt[:2].T

    initial_rejected = np.array([
        row["adaptive_round"] == 0 and not row["propagation_vector_pass"]
        for row in rows
    ])
    initial_retained = np.array([
        row["adaptive_round"] == 0 and row["propagation_vector_pass"]
        for row in rows
    ])
    active_rejected = np.array([
        row["adaptive_round"] > 0 and not row["propagation_vector_pass"]
        for row in rows
    ])
    active_retained = np.array([
        row["adaptive_round"] > 0 and row["propagation_vector_pass"]
        for row in rows
    ])

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.3))
    ax = axes[0]
    ax.scatter(projected[initial_rejected, 0], projected[initial_rejected, 1],
               s=9, c="#c7c7c7", alpha=0.35, label="initial rejected")
    ax.scatter(projected[active_rejected, 0], projected[active_rejected, 1],
               s=18, c="#e45756", alpha=0.40, label="adaptive rejected")
    ax.scatter(projected[initial_retained, 0], projected[initial_retained, 1],
               s=38, c="#4c78a8", edgecolor="black", linewidth=0.3,
               label="initial retained")
    ax.scatter(projected[active_retained, 0], projected[active_retained, 1],
               s=48, c="#f2cf5b", edgecolor="black", linewidth=0.4,
               label="adaptive retained")
    ax.set(xlabel="parameter-space principal coordinate 1",
           ylabel="principal coordinate 2", title="Adaptive 4-D coverage")
    ax.legend(fontsize=7, loc="best")

    ax = axes[1]
    retained = [row for row in rows if row["propagation_vector_pass"]]
    scatter = ax.scatter(
        [row["qz_folded_rlu"] for row in rows],
        [row["critical_field_shift_vs_M3_T"] for row in rows],
        c=[row["adaptive_round"] for row in rows], s=11, alpha=0.35,
        cmap="viridis",
    )
    ax.scatter(
        [row["qz_folded_rlu"] for row in retained],
        [row["critical_field_shift_vs_M3_T"] for row in retained],
        facecolors="none", edgecolors="black", s=44, linewidth=0.7,
    )
    ax.axvspan(*QZ_BAND, color="#54a24b", alpha=0.12,
               label="public qz band")
    ax.set(xlabel=r"folded $q_z$ (r.l.u.)",
           ylabel=r"$H_c$(M4)$-$H_c$(M3) (T)",
           title="Exact retained extrema")
    ax.legend(fontsize=8)
    fig.colorbar(scatter, ax=ax, label="adaptive round")

    ax = axes[2]
    rounds = [item["round"] for item in convergence]
    minima = [item["retained_shift_range_T"][0] for item in convergence]
    maxima = [item["retained_shift_range_T"][1] for item in convergence]
    ax.plot(rounds, minima, "o-", label="minimum retained shift")
    ax.plot(rounds, maxima, "s-", label="maximum retained shift")
    ax.set(xlabel="adaptive round", ylabel="retained shift (T)",
           title="Extremum convergence", xticks=rounds)
    ax.legend(fontsize=8)
    ax2 = ax.twinx()
    ax2.plot(rounds, [item["retained_count"] for item in convergence],
             "^-", color="#e45756", alpha=0.65, label="retained count")
    ax2.set_ylabel("retained exact models")

    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        options = {"dpi": 240} if suffix == "png" else {}
        fig.savefig(FIGURES / f"interlayer_adaptive_inverse_uq.{suffix}",
                    **options)
    plt.close(fig)


def main() -> None:
    initial_rows: list[dict] = []
    for source, path in INITIAL_FILES:
        initial_rows.extend(read_rows(path, source))
    rows = list(initial_rows)
    m3_hc = float(locate_critical_mode("M3", SETTINGS)["critical_field_T"])
    convergence: list[dict] = []
    acquisition_diagnostics: list[dict] = []

    def record_convergence(round_index: int) -> None:
        retained = [row for row in rows if row["propagation_vector_pass"]]
        convergence.append({
            "round": round_index,
            "evaluated_count": len(rows),
            "retained_count": len(retained),
            "retained_shift_range_T": range_of(
                retained, "critical_field_shift_vs_M3_T"
            ),
            "retained_qz_range_rlu": range_of(retained, "qz_folded_rlu"),
        })

    record_convergence(0)
    for round_index in range(1, ACTIVE_ROUNDS + 1):
        selected, roles, diagnostics = acquire(rows, round_index)
        diagnostics["round"] = round_index
        acquisition_diagnostics.append(diagnostics)
        new_rows = []
        for index, (values, role) in enumerate(zip(selected, roles, strict=True), start=1):
            new_rows.append(evaluate_point(
                values, role, round_index, index, m3_hc
            ))
            if index % 8 == 0:
                print(
                    f"adaptive exact evaluation round {round_index}: "
                    f"{index}/{POINTS_PER_ROUND}", flush=True
                )
        rows.extend(new_rows)
        record_convergence(round_index)
        print(json.dumps(convergence[-1], indent=2), flush=True)

    # engineering_rows uses sample_id only as a row identifier.  Reassign a
    # globally unique deterministic id because the two initial designs reuse
    # their own local ids and newly acquired points have no legacy id.
    for sample_id, row in enumerate(rows, start=1):
        row["sample_id"] = sample_id

    retained = [row for row in rows if row["propagation_vector_pass"]]
    adaptive = [row for row in rows if row["adaptive_round"] > 0]
    adaptive_retained = [row for row in adaptive if row["propagation_vector_pass"]]
    retained_sorted = sorted(
        retained, key=lambda row: row["critical_field_shift_vs_M3_T"]
    )
    e90 = [
        row for row in engineering_rows(retained, m3_hc)
        if row["temperature_K"] == 0.09
    ]
    last_change = {
        "minimum_shift_change_T": float(
            convergence[-1]["retained_shift_range_T"][0]
            - convergence[-2]["retained_shift_range_T"][0]
        ),
        "maximum_shift_change_T": float(
            convergence[-1]["retained_shift_range_T"][1]
            - convergence[-2]["retained_shift_range_T"][1]
        ),
    }
    summary = {
        "status": "adaptive exact set-valued inverse stress audit",
        "scope": (
            "Surrogate-guided acquisition with exact multistart M4 decisions; "
            "not a posterior and not a proof of continuous global coverage."
        ),
        "declared_parameter_bounds": {
            key: list(BOUNDS[key]) for key in PARAMETERS
        },
        "public_q_constraints": {
            "delta_inplane_rlu": list(DELTA_BAND),
            "qz_folded_rlu": list(QZ_BAND),
            "km_line_residual_max_rlu": KM_LINE_RESIDUAL_MAX,
        },
        "design": {
            "initial_exact_count": len(initial_rows),
            "adaptive_rounds": ACTIVE_ROUNDS,
            "exact_points_per_round": POINTS_PER_ROUND,
            "adaptive_exact_count": len(adaptive),
            "total_exact_count": len(rows),
            "pool_count_per_round": 2 ** POOL_POWER,
            "acquisition_roles": sorted(set(row["acquisition_role"] for row in adaptive)),
            "multistart_count_per_new_point": 1 + len(OPTIMIZER_SEEDS),
        },
        "cross_validation": cross_validation(initial_rows),
        "acquisition_diagnostics": acquisition_diagnostics,
        "convergence": convergence,
        "last_round_extremum_change": last_change,
        "retained": {
            "initial_count": sum(
                row["propagation_vector_pass"] for row in initial_rows
            ),
            "adaptive_new_count": len(adaptive_retained),
            "total_count": len(retained),
            "outside_original_box_count": sum(
                not inside_original_box(row) for row in retained
            ),
            "critical_field_shift_range_T": range_of(
                retained, "critical_field_shift_vs_M3_T"
            ),
            "delta_range_rlu": range_of(retained, "delta_inplane_rlu"),
            "qz_range_rlu": range_of(retained, "qz_folded_rlu"),
            "out_of_plane_bandwidth_range_microeV": range_of(
                retained, "out_of_plane_bandwidth_microeV"
            ),
            "minimum_shift_model": retained_sorted[0],
            "maximum_shift_model": retained_sorted[-1],
        },
        "central_90mK_engineering_gate": {
            "scenario_count": len(e90),
            "both_gates_pass_count": sum(
                row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
                for row in e90
            ),
            "all_pass": bool(e90) and all(
                row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
                for row in e90
            ),
            "delta_B_star_range_T": range_of(e90, "delta_B_star_T"),
            "regret_range_percent": range_of(e90, "regret_percent"),
        },
        "interpretation_limits": [
            "The acquisition surrogate selects calculations but never supplies a reported physical value.",
            "A finite adaptive audit cannot prove that no disconnected admissible island exists in the continuum.",
            "Retained fractions are deterministic coverage counts and are not probabilities.",
            "The near-saturation propagation vector does not validate the finite-temperature entropy surrogate.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    write_csv(RESULTS / "interlayer_adaptive_inverse_uq.csv", rows)
    write_csv(RESULTS / "interlayer_adaptive_inverse_uq_retained.csv", retained)
    write_csv(RESULTS / "interlayer_adaptive_inverse_uq_new.csv", adaptive)
    (RESULTS / "interlayer_adaptive_inverse_uq_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    make_figure(rows, convergence)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
