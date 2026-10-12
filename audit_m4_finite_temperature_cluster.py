"""Finite-temperature M4 sensitivity on an explicit 3-D periodic cluster.

The engineering surrogate transports a public 2-D XTRG heat profile by its
critical field and then applies a declared shape-stress envelope.  This audit
asks a narrower but more microscopic question: on the same finite 3-D
cluster, how much does adding each diffraction-compatible interlayer model
change the near-saturation entropy/heat profile relative to M3?

The 3 x 2 x 2 cluster is deliberately described as a sensitivity cluster.  It
contains both in-plane and interlayer bonds and the complete periodic Ewald
dipole tensor, but its short y circumference is not a converged representation
of the incommensurate thermodynamic limit.  The calculation can therefore
test the rigid-shift/shape assumption, not replace bulk finite-T XTRG.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bulk_polarized_stability import InterlayerParameters, _bond_exchange_tensor
from dipole_ewald import ewald_tensor, nbcpo_primitive_cell, spin_coupling_microeV
from ewald_cluster import MICROEV_TO_K
from ftlm import run_ftlm
from matrix_free_hamiltonian import SpinBond, SpinHamiltonianOperator
from screen_model import MU_B_OVER_KB, VOLUME_FACTOR_MJ_PER_CM3
from screen_weak_terms import GZ
from stacked_lattice import stacked_triangular_multigraph


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RETAINED = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"
SUMMARY = RESULTS / "m4_finite_temperature_cluster_summary.json"
CURVES = RESULTS / "m4_finite_temperature_cluster_curves.csv"
MODELS = RESULTS / "m4_finite_temperature_cluster_models.csv"
FIGURE_STEM = FIGURES / "m4_finite_temperature_cluster_audit"

PARAMETERS = ("j2_microev", "ja_microev", "jb_microev", "delta_c")
TEMPERATURES = np.array([0.09, 0.10, 0.12, 0.15, 2.0])
COARSE_FIELDS = np.round(np.arange(0.60, 2.1001, 0.05), 10)
PILOT_FIELDS = np.round(np.arange(0.60, 1.8001, 0.10), 10)
HOT_FIELD_T = 4.0
PRIMARY_SEED = 20261012


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_retained() -> list[dict]:
    with RETAINED.open(newline="", encoding="utf-8-sig") as handle:
        rows = []
        for raw in csv.DictReader(handle):
            row = {key: float(raw[key]) for key in PARAMETERS}
            row.update({
                "critical_field_shift_vs_M3_T":
                    float(raw["critical_field_shift_vs_M3_T"]),
                "critical_field_T": float(raw["critical_field_T"]),
                "sample_id": str(raw["sample_id"]),
                "source_design": str(raw["source_design"]),
            })
            rows.append(row)
    return rows


def select_cover(rows: list[dict], target: int) -> list[dict]:
    """Deterministic farthest-point cover in parameters plus bulk shift."""
    keys = (*PARAMETERS, "critical_field_shift_vs_M3_T")
    values = np.array([[row[key] for key in keys] for row in rows], dtype=float)
    span = np.ptp(values, axis=0)
    span[span == 0.0] = 1.0
    normalized = (values - values.min(axis=0)) / span
    shifts = values[:, -1]
    target = max(3, min(target, len(rows)))
    selected = [int(np.argmin(shifts)), int(np.argmax(shifts))]

    default = np.array([2.8, 0.0, 7.9, 0.08], dtype=float)
    pspan = span[:4]
    nearest_default = int(np.argmin(np.linalg.norm(
        (values[:, :4] - default) / pspan, axis=1)))
    if nearest_default not in selected:
        selected.append(nearest_default)

    while len(selected) < target:
        distances = np.min(
            np.linalg.norm(normalized[:, None, :] - normalized[selected][None, :, :], axis=2),
            axis=1,
        )
        distances[selected] = -1.0
        selected.append(int(np.argmax(distances)))
    return [rows[index] for index in selected[:target]]


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_common_bonds(lx: int = 3, ly: int = 2, lz: int = 2):
    graph = stacked_triangular_multigraph(lx, ly, lz)
    primitive = nbcpo_primitive_cell()
    cell = primitive @ np.diag([lx, ly, lz])
    angle = {(1, 0): 0.0, (0, 1): np.pi / 3, (1, -1): -np.pi / 3}
    bonds: list[SpinBond] = []
    for bond in graph["bonds"]:
        if bond.kind != "J1":
            continue
        theta = angle[bond.displacement_primitive[:2]]
        bonds.append(SpinBond(
            bond.i, bond.j, _bond_exchange_tensor(theta, include_jpp=True)))

    # The pair tensor includes all periodic images.  Same-site image terms are
    # a constant for spin 1/2 and cannot affect entropy or the field optimum.
    coordinates = np.asarray(graph["sites"], dtype=float)
    for i in range(graph["N"]):
        for j in range(i + 1, graph["N"]):
            displacement = primitive @ (coordinates[j] - coordinates[i])
            geometric = ewald_tensor(
                displacement, cell, alpha=0.22, real_cutoff=42.0,
                reciprocal_cutoff=8.0, boundary="tin_foil")
            coupling = spin_coupling_microeV(geometric) * MICROEV_TO_K
            bonds.append(SpinBond(i, j, coupling))
    return graph, bonds


def add_interlayer_bonds(graph: dict, common: list[SpinBond],
                         parameters: InterlayerParameters) -> list[SpinBond]:
    bonds = list(common)
    values = {
        "J2": parameters.j2_microev,
        "Ja": parameters.ja_microev,
        "Jb": parameters.jb_microev,
    }
    for bond in graph["bonds"]:
        if bond.kind == "J1":
            continue
        coupling = values[bond.kind]
        tensor = np.diag([
            parameters.delta_c * coupling,
            parameters.delta_c * coupling,
            coupling,
        ]) * MICROEV_TO_K
        bonds.append(SpinBond(bond.i, bond.j, tensor))
    return bonds


def entropy_rows(graph: dict, bonds: list[SpinBond], fields: np.ndarray,
                 seed: int, random_vectors: int, krylov_dim: int,
                 exact_low_states: int,
                 field_seed_stride: int = 7919) -> tuple[list[dict], float]:
    rows = []
    elapsed = 0.0
    for index, field in enumerate(fields):
        operator = SpinHamiltonianOperator(
            graph["N"], bonds, GZ * MU_B_OVER_KB * float(field), parity=None)
        result = run_ftlm(
            [operator], TEMPERATURES, graph["N"],
            random_vectors=random_vectors, krylov_dim=krylov_dim,
            exact_low_states=exact_low_states,
            seed=seed + field_seed_stride * index)
        elapsed += result["elapsed_seconds"]
        for row in result["rows"]:
            rows.append({"field_T": float(field), **row})
    return rows, elapsed


def extract_entropy(rows: list[dict], temperature: float, field: float) -> float:
    selected = [row for row in rows
                if np.isclose(row["temperature_K"], temperature)
                and np.isclose(row["field_T"], field)]
    if len(selected) != 1:
        raise RuntimeError("entropy coordinate missing or duplicated")
    return float(selected[0]["entropy_kB_per_spin"])


def refine_fields(rows: list[dict], hot_entropy: float) -> np.ndarray:
    cold = [row for row in rows if np.isclose(row["temperature_K"], 0.09)
            and row["field_T"] < 3.0]
    q = np.array([
        0.09 * (row["entropy_kB_per_spin"] - hot_entropy)
        * VOLUME_FACTOR_MJ_PER_CM3 for row in cold])
    center = float(cold[int(np.argmax(q))]["field_T"])
    return np.round(np.arange(center - 0.04, center + 0.0401, 0.01), 10)


def profile_metrics(rows: list[dict], label: str, model_meta: dict) -> tuple[list[dict], list[dict]]:
    hot_entropy = extract_entropy(rows, 2.0, HOT_FIELD_T)
    curve_rows = []
    metrics = []
    for temperature in TEMPERATURES[:-1]:
        selected = sorted(
            [row for row in rows if np.isclose(row["temperature_K"], temperature)
             and row["field_T"] < 3.0], key=lambda row: row["field_T"])
        # Remove a duplicate field if a refinement coordinate coincides with
        # the coarse grid.
        by_field = {float(row["field_T"]): row for row in selected}
        selected = [by_field[field] for field in sorted(by_field)]
        fields = np.array([row["field_T"] for row in selected])
        entropy = np.array([row["entropy_kB_per_spin"] for row in selected])
        heat = temperature * (entropy - hot_entropy) * VOLUME_FACTOR_MJ_PER_CM3
        peak = int(np.argmax(heat))
        qmax = float(heat[peak])
        accepted = fields[heat >= 0.95 * qmax]
        metrics.append({
            "model_label": label,
            **model_meta,
            "temperature_K": float(temperature),
            "hot_entropy_kB_per_spin": hot_entropy,
            "B_star_T": float(fields[peak]),
            "Qmax_mJ_per_cm3": qmax,
            "W95_min_T": float(accepted.min()),
            "W95_max_T": float(accepted.max()),
            "W95_left_width_T": float(fields[peak] - accepted.min()),
            "W95_right_width_T": float(accepted.max() - fields[peak]),
        })
        for field, s, q in zip(fields, entropy, heat):
            curve_rows.append({
                "model_label": label,
                **model_meta,
                "temperature_K": float(temperature),
                "field_T": float(field),
                "entropy_kB_per_spin": float(s),
                "Q_mJ_per_cm3": float(q),
                "hot_entropy_kB_per_spin": hot_entropy,
            })
    return curve_rows, metrics


def compare_metrics(all_metrics: list[dict]) -> list[dict]:
    baseline = {
        row["temperature_K"]: row for row in all_metrics
        if row["model_label"] == "M3"
    }
    comparisons = []
    for row in all_metrics:
        if row["model_label"] == "M3":
            continue
        base = baseline[row["temperature_K"]]
        left_scale = (row["W95_left_width_T"] / base["W95_left_width_T"]
                      if base["W95_left_width_T"] > 0 else None)
        right_scale = (row["W95_right_width_T"] / base["W95_right_width_T"]
                       if base["W95_right_width_T"] > 0 else None)
        comparisons.append({
            **row,
            "cluster_M4_minus_M3_B_star_T": row["B_star_T"] - base["B_star_T"],
            "Qmax_ratio_to_M3": row["Qmax_mJ_per_cm3"] / base["Qmax_mJ_per_cm3"],
            "left_W95_scale_to_M3": left_scale,
            "right_W95_scale_to_M3": right_scale,
        })
    return comparisons


def plot(curve_rows: list[dict], comparisons: list[dict], selected_models: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.5), constrained_layout=True)
    ax = axes[0, 0]
    baseline = [row for row in curve_rows if row["model_label"] == "M3"
                and np.isclose(row["temperature_K"], 0.09)]
    ax.plot([row["field_T"] for row in baseline],
            [row["Q_mJ_per_cm3"] for row in baseline],
            color="black", lw=1.8, label="M3")
    for model in selected_models:
        label = f"M4-{model['cover_index']:02d}"
        rows = [row for row in curve_rows if row["model_label"] == label
                and np.isclose(row["temperature_K"], 0.09)]
        ax.plot([row["field_T"] for row in rows],
                [row["Q_mJ_per_cm3"] for row in rows],
                color="#3b82f6", alpha=0.28, lw=0.9)
    ax.set(xlabel="Field (T)", ylabel=r"$Q$ proxy (mJ cm$^{-3}$)",
           title="(a) N=12 microscopic 90 mK profiles")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)

    comp90 = [row for row in comparisons if np.isclose(row["temperature_K"], 0.09)]
    ax = axes[0, 1]
    ax.scatter([row["critical_field_shift_vs_M3_T"] for row in comp90],
               [row["cluster_M4_minus_M3_B_star_T"] for row in comp90],
               c=[row["delta_c"] for row in comp90], cmap="viridis", s=38)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set(xlabel=r"Bulk $Delta H_c$ (T)", ylabel=r"Cluster $Delta B^*$ (T)",
           title="(b) Field-shift transfer")
    ax.grid(alpha=0.2)

    ax = axes[1, 0]
    ax.scatter([row["critical_field_shift_vs_M3_T"] for row in comp90],
               [row["Qmax_ratio_to_M3"] for row in comp90],
               color="#ef4444", s=38)
    ax.axhspan(0.9, 1.1, color="#f59e0b", alpha=0.14,
               label="declared heat-scale stress")
    ax.axhline(1, color="0.4", lw=0.8)
    ax.set(xlabel=r"Bulk $Delta H_c$ (T)", ylabel=r"$Q_{max}$(M4)/$Q_{max}$(M3)",
           title="(c) Peak-capacity deformation")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    left = [row["left_W95_scale_to_M3"] for row in comp90
            if row["left_W95_scale_to_M3"] is not None]
    right = [row["right_W95_scale_to_M3"] for row in comp90
             if row["right_W95_scale_to_M3"] is not None]
    x = np.arange(len(comp90))
    ax.scatter(x[:len(left)], left, marker="o", color="#10b981", label="left")
    ax.scatter(x[:len(right)], right, marker="s", color="#8b5cf6", label="right")
    ax.axhspan(0.8, 1.2, color="#f59e0b", alpha=0.14,
               label="declared width stress")
    ax.axhline(1, color="0.4", lw=0.8)
    ax.set(xlabel="Coverage model index", ylabel="W95 width / M3 width",
           title="(d) Profile-width deformation")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURE_STEM.with_suffix(f".{suffix}"),
                    dpi=300 if suffix == "png" else None)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--cover-count", type=int, default=12)
    parser.add_argument("--random-vectors", type=int, default=8)
    parser.add_argument("--krylov-dim", type=int, default=60)
    parser.add_argument("--exact-low-states", type=int, default=12)
    args = parser.parse_args()

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    retained = read_retained()
    selected = select_cover(retained, 1 if args.pilot else args.cover_count)
    for index, row in enumerate(selected):
        row["cover_index"] = index
    model_rows = [{"model_label": "M3", "cover_index": -1,
                   **{key: 0.0 for key in PARAMETERS},
                   "critical_field_shift_vs_M3_T": 0.0,
                   "critical_field_T": None, "sample_id": "M3",
                   "source_design": "baseline"}, *[
        {"model_label": f"M4-{row['cover_index']:02d}", **row}
        for row in selected]]
    write_csv(MODELS, model_rows)

    build_start = time.perf_counter()
    graph, common = build_common_bonds()
    build_seconds = time.perf_counter() - build_start
    all_rows: list[dict] = []
    all_metrics: list[dict] = []
    timings = []
    fields = PILOT_FIELDS if args.pilot else COARSE_FIELDS

    for model_index, meta in enumerate(model_rows):
        if meta["model_label"] == "M3":
            parameters = InterlayerParameters(0.0, 0.0, 0.0, 0.08)
        else:
            parameters = InterlayerParameters(**{
                key: float(meta[key]) for key in PARAMETERS})
        bonds = add_interlayer_bonds(graph, common, parameters)
        initial_fields = np.unique(np.r_[fields, HOT_FIELD_T])
        rows, elapsed = entropy_rows(
            graph, bonds, initial_fields, PRIMARY_SEED,
            args.random_vectors, args.krylov_dim, args.exact_low_states)
        hot = extract_entropy(rows, 2.0, HOT_FIELD_T)
        if not args.pilot:
            refine = np.setdiff1d(refine_fields(rows, hot), initial_fields)
            if refine.size:
                added, refine_elapsed = entropy_rows(
                    graph, bonds, refine, PRIMARY_SEED,
                    args.random_vectors, args.krylov_dim, args.exact_low_states)
                rows.extend(added)
                elapsed += refine_elapsed
        curve_rows, metrics = profile_metrics(rows, meta["model_label"], {
            key: meta.get(key) for key in
            ("cover_index", *PARAMETERS, "critical_field_shift_vs_M3_T",
             "critical_field_T", "sample_id", "source_design")})
        all_rows.extend(curve_rows)
        all_metrics.extend(metrics)
        timings.append({"model_label": meta["model_label"],
                        "elapsed_seconds": elapsed})
        print(json.dumps({
            "completed": model_index + 1,
            "total": len(model_rows),
            "model": meta["model_label"],
            "elapsed_seconds": elapsed,
            "B_star_90mK_T": metrics[0]["B_star_T"],
            "Qmax_90mK_mJ_per_cm3": metrics[0]["Qmax_mJ_per_cm3"],
        }), flush=True)

    comparisons = compare_metrics(all_metrics)
    comp90 = [row for row in comparisons if np.isclose(row["temperature_K"], 0.09)]
    write_csv(CURVES, all_rows)
    write_csv(RESULTS / "m4_finite_temperature_cluster_metrics.csv", all_metrics)
    write_csv(RESULTS / "m4_finite_temperature_cluster_comparisons.csv", comparisons)
    plot(all_rows, comparisons, selected)

    def finite_values(key: str) -> list[float]:
        return [float(row[key]) for row in comp90 if row[key] is not None]

    def finite_range(key: str) -> list[float | None]:
        values = finite_values(key)
        return [min(values), max(values)] if values else [None, None]

    payload = {
        "status": "finite 3-D cluster thermodynamic sensitivity; not a converged bulk M4 solution",
        "retained_source": {
            "path": str(RETAINED.relative_to(ROOT)),
            "sha256": sha256(RETAINED),
            "retained_count": len(retained),
        },
        "cluster": {
            "shape": list(graph["shape"]),
            "spins": graph["N"],
            "boundary": "tin_foil",
            "common_bond_count": len(common),
            "three_sublattice_compatible": False,
            "role": "near-saturation differential sensitivity cluster",
        },
        "selection": {
            "method": "farthest-point coverage in four interlayer parameters plus bulk critical-field shift, seeded by both shift extrema and the nearest representative model",
            "cover_count": len(selected),
            "models": model_rows,
        },
        "ftlm": {
            "primary_seed": PRIMARY_SEED,
            "random_vectors": args.random_vectors,
            "krylov_dimension": args.krylov_dim,
            "exact_low_states": args.exact_low_states,
            "temperatures_K": TEMPERATURES.tolist(),
            "coarse_fields_T": fields.tolist(),
            "hot_state": {"temperature_K": 2.0, "field_T": HOT_FIELD_T},
            "build_seconds": build_seconds,
            "timings": timings,
        },
        "90mK_differential_ranges": {
            "cluster_M4_minus_M3_B_star_T":
                finite_range("cluster_M4_minus_M3_B_star_T"),
            "Qmax_ratio_to_M3": finite_range("Qmax_ratio_to_M3"),
            "left_W95_scale_to_M3": finite_range("left_W95_scale_to_M3"),
            "right_W95_scale_to_M3": finite_range("right_W95_scale_to_M3"),
        },
        "declared_shape_stress_contains_cluster_cover": {
            "heat_scale_0p9_to_1p1": all(
                0.9 <= value <= 1.1 for value in finite_values("Qmax_ratio_to_M3")),
            "left_width_0p8_to_1p2": all(
                0.8 <= value <= 1.2 for value in finite_values("left_W95_scale_to_M3")),
            "right_width_0p8_to_1p2": all(
                0.8 <= value <= 1.2 for value in finite_values("right_W95_scale_to_M3")),
        },
        "metrics": all_metrics,
        "comparisons": comparisons,
        "limitations": [
            "The 3 x 2 x 2 cluster contains a short-circumference seam and is not a converged representation of the incommensurate bulk state.",
            "FTLM has stochastic trace error; common settings and seed make M4-M3 differences correlated but do not create a confidence interval.",
            "This audit tests differential profile deformation near saturation. It does not replace the measured entropy surface or a bulk 3-D finite-temperature solver.",
        ],
    }
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "summary": str(SUMMARY),
        "ranges_90mK": payload["90mK_differential_ranges"],
        "stress_contains_cover": payload["declared_shape_stress_contains_cluster_cover"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
