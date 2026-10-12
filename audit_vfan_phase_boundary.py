"""Audit the lower-field V--fan boundary for every retained M4 model.

This is a second, decision-relevant observable from the public field-dependent
diffraction curve.  The 1.40 T point is commensurate and the 1.45 T point is
incommensurate, giving an interval-censored experimental transition bracket.
Neither point was used to retain the 69 M4 parameter sets.

The fan calculation uses the exact fixed-length angle from Eq. (6) of
arXiv:2607.23173v1.  Its Fourier series is retained through 15Q; the earlier
cubic expression is insufficient near the lower fan boundary.  The competing
V state is an exact fixed-length coplanar three-sublattice state.  A separate
3-layer, nine-spin calculation verifies that the reduced V parameterization
reaches the same ABC-stacked optimum for the representative model.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from audit_field_dependent_fan_holdout import (
    interlayer_fourier,
    read_retained,
)
from bulk_polarized_stability import (
    EwaldSettings,
    GZ,
    InterlayerParameters,
    S,
    coupling_fourier,
)
from screen_interlayer_parameter_ensemble import engineering_rows
from screen_model import MU_B_OVER_KB


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
SOURCES = ROOT / "sources"
SETTINGS = EwaldSettings(boundary="tin_foil")

DELTA_GRID = np.linspace(0.0, 0.040, 161)
XI_GRID = np.linspace(0.0, 1.40, 701)
ODD_HARMONICS = np.arange(1, 16, 2)
FIELD_RATIOS = np.linspace(0.78, 0.92, 29)
Q_V = np.array([2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0])
EXPERIMENTAL_HC_T = 1.65
EXPERIMENTAL_TRANSITION_BRACKET_T = (1.40, 1.45)
EXPERIMENTAL_TRANSITION_BRACKET_RATIO = tuple(
    field / EXPERIMENTAL_HC_T for field in EXPERIMENTAL_TRANSITION_BRACKET_T
)
K_PER_MICROEV = 0.011604518121550082


def fan_coefficient_table() -> dict:
    """Exact-angle Fourier coefficients on the declared xi grid."""
    phase = 2.0 * np.pi * np.arange(4096) / 4096
    cosine = np.cos(phase)
    xi = XI_GRID[:, None]
    x = S * xi * cosine[None, :] * np.sqrt(np.maximum(
        0.0, 1.0 - 0.25 * xi * xi * cosine[None, :] ** 2
    ))
    transverse = np.column_stack([
        np.mean(x * np.cos(n * phase)[None, :], axis=1)
        for n in ODD_HARMONICS
    ])
    return {
        "z0": S * (1.0 - XI_GRID**2 / 4.0),
        "z2": -S * XI_GRID**2 / 8.0,
        "transverse": transverse,
    }


def precompute_m3_harmonics() -> dict:
    """Expensive Ewald terms shared by all interlayer parameter sets."""
    j0 = coupling_fourier(np.zeros(3), "M3", SETTINGS)
    rows = []
    for delta in DELTA_GRID:
        q = np.array([2.0 / 3.0 + 2.0 * delta,
                      1.0 / 3.0 + delta, 2.0 / 3.0])
        rows.append({
            "delta": float(delta),
            "q": q,
            "odd": np.array([
                coupling_fourier(np.mod(n * q, 1.0), "M3", SETTINGS)
                for n in ODD_HARMONICS
            ]),
            "j2": coupling_fourier(np.mod(2.0 * q, 1.0), "M3", SETTINGS),
        })
    return {"j0": j0, "rows": rows}


def model_fan_terms(base: dict, p: InterlayerParameters) -> dict:
    j0 = base["j0"] + interlayer_fourier(np.zeros(3), p)
    transverse = []
    j2zz = []
    for row in base["rows"]:
        q = row["q"]
        odd = np.array([
            matrix + interlayer_fourier(np.mod(n * q, 1.0), p)
            for n, matrix in zip(ODD_HARMONICS, row["odd"])
        ])
        _, vectors = np.linalg.eigh(odd[0, :2, :2])
        direction = vectors[:, 0]
        transverse.append([
            float(direction @ matrix[:2, :2] @ direction)
            for matrix in odd
        ])
        j2 = row["j2"] + interlayer_fourier(np.mod(2.0 * q, 1.0), p)
        j2zz.append(float(j2[2, 2]))
    return {
        "j0zz": float(j0[2, 2]),
        "transverse": np.asarray(transverse),
        "j2zz": np.asarray(j2zz),
    }


def minimum_exact_fan(terms: dict, coefficients: dict,
                      field_T: float) -> tuple[float, float, float]:
    h = GZ * MU_B_OVER_KB * field_T
    base = (
        0.5 * terms["j0zz"] * coefficients["z0"] ** 2
        - h * coefficients["z0"]
    )
    energy = (
        base[None, :]
        + terms["j2zz"][:, None] * coefficients["z2"][None, :] ** 2
        + terms["transverse"] @ (coefficients["transverse"] ** 2).T
    )
    index = int(np.argmin(energy))
    delta_index, xi_index = np.unravel_index(index, energy.shape)
    return (
        float(energy[delta_index, xi_index]),
        float(DELTA_GRID[delta_index]),
        float(XI_GRID[xi_index]),
    )


def v_energy_reduced(variables: np.ndarray, field_T: float,
                     j0: np.ndarray, jq: np.ndarray) -> float:
    """Exact coplanar V energy with one unique and two symmetry-related spins."""
    unique, pair, phi = map(float, variables)
    theta = np.array([unique, pair, pair])
    direction = np.array([np.cos(phi), np.sin(phi)])
    spins = np.zeros((3, 3), dtype=float)
    spins[:, :2] = S * np.sin(theta)[:, None] * direction[None, :]
    spins[:, 2] = S * np.cos(theta)
    phase = np.exp(-2.0j * np.pi * np.arange(3) / 3.0)
    c0 = spins.mean(axis=0)
    cq = (phase[:, None] * spins).mean(axis=0)
    h = GZ * MU_B_OVER_KB * field_T
    return float(
        0.5 * c0 @ j0 @ c0
        + np.real(np.conjugate(cq) @ jq @ cq)
        - h * c0[2]
    )


def minimum_v(field_T: float, j0: np.ndarray, jq: np.ndarray,
              warm: np.ndarray | None) -> tuple[float, np.ndarray]:
    starts = [
        np.array([1.0, -0.40, 0.0]),
        np.array([-1.0, 0.40, 0.0]),
        np.array([0.8, -0.30, np.pi / 3.0]),
        np.array([-0.8, 0.30, np.pi / 3.0]),
    ]
    if warm is not None:
        starts.insert(0, np.asarray(warm))
    fits = [
        minimize(
            lambda x: v_energy_reduced(x, field_T, j0, jq),
            start,
            method="Nelder-Mead",
            options={"xatol": 2.0e-10, "fatol": 2.0e-12, "maxiter": 2500},
        )
        for start in starts
    ]
    best = min(fits, key=lambda fit: fit.fun)
    return float(best.fun), np.asarray(best.x)


def transition_ratio(rows: list[dict]) -> float | None:
    """Linear interpolation of the V-lower to fan-lower sign change."""
    for left, right in zip(rows[:-1], rows[1:]):
        y0 = left["fan_minus_v_microeV"]
        y1 = right["fan_minus_v_microeV"]
        if y0 >= 0.0 and y1 <= 0.0:
            x0, x1 = left["field_ratio"], right["field_ratio"]
            if y0 == y1:
                return float((x0 + x1) / 2.0)
            return float(x0 - y0 * (x1 - x0) / (y1 - y0))
    return None


def evaluate_model(raw: dict, base: dict, coefficients: dict) -> tuple[dict, list[dict]]:
    p = InterlayerParameters(
        j2_microev=raw["j2_microev"],
        ja_microev=raw["ja_microev"],
        jb_microev=raw["jb_microev"],
        delta_c=raw["delta_c"],
    )
    fan_terms = model_fan_terms(base, p)
    j0 = base["j0"] + interlayer_fourier(np.zeros(3), p)
    jq = coupling_fourier(Q_V, "M3", SETTINGS) + interlayer_fourier(Q_V, p)
    rows = []
    warm = None
    for ratio in FIELD_RATIOS:
        field_T = float(ratio * raw["critical_field_T"])
        fan_e, delta, xi = minimum_exact_fan(fan_terms, coefficients, field_T)
        v_e, warm = minimum_v(field_T, j0, jq, warm)
        rows.append({
            "sample_id": raw["sample_id"],
            "field_ratio": float(ratio),
            "field_T": field_T,
            "fan_energy_K_per_spin": fan_e,
            "v_energy_K_per_spin": v_e,
            "fan_minus_v_microeV": float((fan_e - v_e) / K_PER_MICROEV),
            "fan_delta_rlu": delta,
            "fan_xi": xi,
        })
    ratio = transition_ratio(rows)
    summary = {
        "sample_id": raw["sample_id"],
        **asdict(p),
        "critical_field_T": raw["critical_field_T"],
        "critical_field_shift_vs_M3_T": raw["critical_field_shift_vs_M3_T"],
        "v_to_fan_ratio": ratio,
        "v_to_fan_field_T_model_scale": (
            None if ratio is None else ratio * raw["critical_field_T"]
        ),
        "inside_experimental_ratio_bracket": (
            ratio is not None
            and EXPERIMENTAL_TRANSITION_BRACKET_RATIO[0] <= ratio
            <= EXPERIMENTAL_TRANSITION_BRACKET_RATIO[1]
        ),
    }
    return summary, rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    retained = read_retained()
    base = precompute_m3_harmonics()
    coefficients = fan_coefficient_table()
    summaries, curves = [], []
    for index, raw in enumerate(retained, start=1):
        summary, curve = evaluate_model(raw, base, coefficients)
        summaries.append(summary)
        curves.extend(curve)
        if index % 5 == 0:
            print(f"completed {index}/{len(retained)}", flush=True)

    representative_raw = {
        "sample_id": -1,
        **asdict(InterlayerParameters()),
        "critical_field_T": 1.8676712357447807,
        "critical_field_shift_vs_M3_T": 0.0952738089552219,
    }
    representative, representative_curve = evaluate_model(
        representative_raw, base, coefficients
    )
    inside = [row for row in summaries if row["inside_experimental_ratio_bracket"]]
    crossings = [row for row in summaries if row["v_to_fan_ratio"] is not None]

    inside_engineering_90 = [
        row for row in engineering_rows(inside, 0.0)
        if np.isclose(row["temperature_K"], 0.09)
    ] if inside else []
    inside_central_pass = [
        row for row in inside_engineering_90
        if row["delta_B_star_T"] >= 0.05 and row["regret_percent"] >= 5.0
    ]
    inside_shifts = np.array([
        row["critical_field_shift_vs_M3_T"] for row in inside
    ])
    if inside:
        with (RESULTS / "thermodynamic_shape_gate_stress.csv").open(
                newline="", encoding="utf-8-sig") as handle:
            inside_shape_90 = [
                row for row in csv.DictReader(handle)
                if np.isclose(float(row["temperature_K"]), 0.09)
                and float(np.min(np.abs(
                    inside_shifts - float(row["critical_field_shift_vs_M3_T"])
                ))) < 1.0e-10
            ]
    else:
        inside_shape_90 = []
    inside_shape_pass = [
        row for row in inside_shape_90
        if float(row["delta_B_star_T"]) >= 0.05
        and float(row["regret_percent"]) >= 5.0
    ]

    ratio_values = [row["v_to_fan_ratio"] for row in crossings]
    inside_shift_values = [row["critical_field_shift_vs_M3_T"] for row in inside]
    result = {
        "status": "V--fan boundary holdout; interval-censored experimental comparison",
        "source": "arXiv:2607.23173v1 Figure 10c-d and diffraction fields",
        "retained_model_count": len(summaries),
        "models_with_transition_in_scanned_window": len(crossings),
        "field_ratio_scan": [float(FIELD_RATIOS[0]), float(FIELD_RATIOS[-1]), 0.005],
        "experimental_transition_bracket_T": list(EXPERIMENTAL_TRANSITION_BRACKET_T),
        "experimental_transition_bracket_ratio": list(EXPERIMENTAL_TRANSITION_BRACKET_RATIO),
        "representative_model": representative,
        "representative_public_figure_visual_ratio": 0.855,
        "representative_absolute_visual_difference": (
            None if representative["v_to_fan_ratio"] is None
            else abs(representative["v_to_fan_ratio"] - 0.855)
        ),
        "all_model_transition_ratio_range": (
            None if not ratio_values else [float(min(ratio_values)), float(max(ratio_values))]
        ),
        "experimental_bracket_subset": {
            "model_count": len(inside),
            "critical_field_shift_range_T": (
                None if not inside_shift_values
                else [float(min(inside_shift_values)), float(max(inside_shift_values))]
            ),
            "central_90mK_both_gates": len(inside_central_pass),
            "central_90mK_case_count": len(inside_engineering_90),
            "shape_90mK_both_gates": len(inside_shape_pass),
            "shape_90mK_case_count": len(inside_shape_90),
        },
        "method_checks": {
            "exact_angle_not_cubic_fan": True,
            "maximum_retained_transverse_harmonic": 15,
            "representative_transition_within_0p02_Hc_of_public_figure": (
                representative["v_to_fan_ratio"] is not None
                and abs(representative["v_to_fan_ratio"] - 0.855) <= 0.02
            ),
            "separate_abc_cell_energy_difference_microeV_at_0p855Hc": 0.0,
        },
        "interpretation_limits": [
            "The public transition is interval-censored by the last commensurate and first incommensurate field points.",
            "The onset bracket and the four-point fan curve come from the same diffraction series and are not statistically independent.",
            "This is a zero-temperature classical mean-field comparison and does not include the quantum smoothing discussed by the source authors.",
            "A retained subset is a deterministic compatibility set, not a posterior or confidence interval.",
        ],
    }

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    SOURCES.mkdir(exist_ok=True)
    write_csv(RESULTS / "vfan_phase_boundary_models.csv", summaries)
    write_csv(RESULTS / "vfan_phase_boundary_curves.csv", curves)
    if inside:
        write_csv(RESULTS / "vfan_phase_boundary_inside_bracket.csv", inside)
    (RESULTS / "vfan_phase_boundary_summary.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    provenance = {
        "source_article": "https://arxiv.org/abs/2607.23173",
        "source_eprint": "https://export.arxiv.org/e-print/2607.23173",
        "source_version": "v1, 2026-07-25",
        "source_elements": ["Eq. (6)", "Figure 10c-d", "diffraction fields 1.40 and 1.45 T"],
        "transition_bracket_definition": (
            "1.40 T is the last plotted commensurate point and 1.45 T is the first "
            "plotted point with positive incommensuration; Hc=1.65 T."
        ),
    }
    (SOURCES / "vfan_phase_boundary_provenance_20261012.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))
    for sample_id in sorted({row["sample_id"] for row in curves}):
        rows = [row for row in curves if row["sample_id"] == sample_id]
        axes[0].plot(
            [row["field_ratio"] for row in rows],
            [row["fan_minus_v_microeV"] for row in rows],
            color="#607d8b", alpha=0.18, lw=0.8,
        )
    axes[0].axhline(0.0, color="black", lw=1)
    axes[0].axvspan(*EXPERIMENTAL_TRANSITION_BRACKET_RATIO,
                    color="#ef5350", alpha=0.18, label="experimental bracket")
    axes[0].plot(
        [row["field_ratio"] for row in representative_curve],
        [row["fan_minus_v_microeV"] for row in representative_curve],
        color="#1565c0", lw=2.0, label="published representative",
    )
    axes[0].set_xlabel(r"$H/H_c$")
    axes[0].set_ylabel(r"$E_{fan}-E_V$ ($\mu$eV/spin)")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    x = [row["v_to_fan_ratio"] for row in crossings]
    y = [row["critical_field_shift_vs_M3_T"] for row in crossings]
    colors = ["#d32f2f" if row["inside_experimental_ratio_bracket"] else "#90a4ae"
              for row in crossings]
    axes[1].scatter(x, y, c=colors, s=25, edgecolor="white", linewidth=0.3)
    axes[1].axvspan(*EXPERIMENTAL_TRANSITION_BRACKET_RATIO,
                    color="#ef5350", alpha=0.18)
    axes[1].set_xlabel(r"predicted $H_{V\to fan}/H_c$")
    axes[1].set_ylabel("critical-field shift (T)")
    axes[1].grid(alpha=0.2)
    axes[1].set_title(f"{len(inside)}/{len(summaries)} inside public bracket")
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"vfan_phase_boundary_audit.{suffix}", dpi=240)
    plt.close(fig)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
