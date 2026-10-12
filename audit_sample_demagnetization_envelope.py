"""Convert public NBCP sample geometry into an applied/internal-field bound.

The calculation uses an equal-volume spheroid as a transparent engineering
shape model.  Saturation magnetization gives a conservative correction
bound: B_internal = B_applied - N_z * mu0*M.  It does not infer the unknown
dimensions of the published 1.65 g sample column.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

MU0 = 4.0 * math.pi * 1e-7
MU_B = 9.2740100783e-24
G_C = 4.79  # public B||c magnetometry value used with the geometry source
A_M = 5.312e-10
C_M = 6.983e-10
CELL_VOLUME_M3 = math.sqrt(3.0) * A_M**2 * C_M / 2.0
MU0_MS_T = MU0 * G_C * MU_B / (2.0 * CELL_VOLUME_M3)


def nz_spheroid(aspect_c_over_a):
    """Demagnetizing factor along c for a spheroid with semiaxes (a,a,c)."""
    r = np.asarray(aspect_c_over_a, dtype=float)
    out = np.empty_like(r)
    sphere = np.isclose(r, 1.0, atol=1e-10)
    oblate = r < 1.0 - 1e-10
    prolate = r > 1.0 + 1e-10
    out[sphere] = 1.0 / 3.0
    xi = np.sqrt(1.0 / r[oblate] ** 2 - 1.0)
    out[oblate] = (1.0 + xi**2) / xi**3 * (xi - np.arctan(xi))
    eccentricity = np.sqrt(1.0 - 1.0 / r[prolate] ** 2)
    out[prolate] = ((1.0 - eccentricity**2) / (2.0 * eccentricity**3)
                    * (np.log((1.0 + eccentricity) / (1.0 - eccentricity))
                       - 2.0 * eccentricity))
    return out


def aspect_for_bound(bound_T):
    low, high = 1.0, 100.0
    for _ in range(100):
        middle = 0.5 * (low + high)
        correction = float(nz_spheroid([middle])[0] * MU0_MS_T)
        if correction > bound_T:
            low = middle
        else:
            high = middle
    return high


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    aspects = np.geomspace(0.05, 20.0, 401)
    nz = nz_spheroid(aspects)
    correction = nz * MU0_MS_T
    rows = [{
        "aspect_c_over_inplane_radius": float(r),
        "demagnetizing_factor_Nz": float(n),
        "saturation_field_correction_bound_T": float(d),
    } for r, n, d in zip(aspects, nz, correction)]
    with (RESULTS / "sample_demagnetization_envelope.csv").open(
            "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # Public B||c heat-capacity crystal: field along the 0.18 mm thickness.
    # Equal-area in-plane semiaxis sqrt((1.87/2)*(1.16/2)).
    plate_aspect = (0.18 / 2.0) / math.sqrt((1.87 / 2.0) * (1.16 / 2.0))
    plate_nz = float(nz_spheroid([plate_aspect])[0])
    plate_bound = plate_nz * MU0_MS_T
    targets = {
        "0.010_T": aspect_for_bound(0.010),
        "0.020_T": aspect_for_bound(0.020),
        "0.050_T": aspect_for_bound(0.050),
    }
    payload = {
        "status": "PASS: sample shape is a first-order field-coordinate requirement",
        "model": "equal-volume spheroid, field parallel to c axis",
        "relation": "B_internal = B_applied - N_z * mu0*M; saturation M gives the reported upper bound",
        "material_constants": {
            "g_c": G_C,
            "cell_a_A": A_M / 1e-10,
            "cell_c_A": C_M / 1e-10,
            "one_effective_spin_half_per_cell": True,
            "mu0_Ms_T": MU0_MS_T,
        },
        "public_B_parallel_c_specific_heat_sample": {
            "dimensions_mm": [1.87, 1.16, 0.18],
            "field_axis_dimension_mm": 0.18,
            "equivalent_spheroid_aspect_c_over_a": plate_aspect,
            "equivalent_spheroid_Nz": plate_nz,
            "saturation_correction_bound_T": plate_bound,
            "fraction_of_1p70T_percent": 100.0 * plate_bound / 1.70,
            "source_reported_upper_modification": "<8% for B parallel c",
            "interpretation": "the simple spheroid reproduces the published correction scale; it is not an exact rectangular-prism correction",
        },
        "shape_acceptance_requirements": {
            key: {
                "minimum_aspect_c_over_inplane_radius": value,
                "meaning": f"conservative saturation correction <= {float(key[:5]):.3f} T",
            }
            for key, value in targets.items()
        },
        "sphere_saturation_correction_T": MU0_MS_T / 3.0,
        "unknown_public_1p65g_column": {
            "mass_g": 1.65,
            "volume_cm3_from_public_comparison": 0.395,
            "dimensions_publicly_available": False,
            "decision": "do not assign an exact applied-field set point until the packed-column dimensions, orientation, and magnetization curve are specified",
        },
        "engineering_conclusion": {
            "0p02T_previous_stopping_band_survives_shape_unknown": False,
            "reason": "a sphere alone gives a correction above 0.05 T and the public plate-scale bound is about 0.13 T",
            "actionable_specification": "align the long refrigerant dimension with c and meet the reported aspect ratio requirement, or solve the measured M(B) self-consistently for the actual package geometry",
        },
        "sources": [
            {
                "citation": "Possible itinerant excitations and quantum spin state transitions in Na2BaCo(PO4)2, Nature Communications 11, 3996 (2020)",
                "doi": "10.1038/s41467-020-18041-3",
                "use": "1.87 x 1.16 x 0.18 mm^3 B||c heat-capacity sample and reported <8% demagnetization modification",
            },
            {
                "citation": "Giant magnetocaloric effect in spin supersolid candidate Na2BaCo(PO4)2, Nature 625, 270-275 (2024)",
                "doi": "10.1038/s41586-023-06885-w",
                "use": "1.65 g sample-column mass and 0.395 cm^3 comparison volume; exact column dimensions are not reported",
            },
        ],
        "limitations": [
            "The spheroid is a package-design envelope rather than an exact rectangular-prism solution.",
            "Saturation magnetization is conservative; a field-dependent correction requires public M(B,T) on the same sample and geometry.",
            "Porosity, grains, binders, and local fields in a packed pill require an effective-medium or measured correction.",
        ],
    }
    (RESULTS / "sample_demagnetization_envelope_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    ax.semilogx(aspects, 1e3 * correction, color="#3569a8", lw=2)
    ax.scatter([plate_aspect], [1e3 * plate_bound], color="#d34a33", zorder=4,
               label="public B||c heat-capacity plate")
    for bound in (0.01, 0.02, 0.05):
        aspect = aspect_for_bound(bound)
        ax.axhline(1e3 * bound, color="0.55", lw=0.8, ls=":")
        ax.axvline(aspect, color="0.55", lw=0.8, ls="--")
        ax.text(aspect * 1.03, 1e3 * bound * 1.08,
                f"c/a >= {aspect:.2f} for {1e3*bound:.0f} mT", fontsize=8)
    ax.set(xlabel="Refrigerant aspect ratio c/a (field along c)",
           ylabel="Saturation demagnetization bound (mT)",
           title="Package shape converts internal-field physics to applied-field control")
    ax.grid(alpha=0.2, which="both")
    ax.legend(frameon=False)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"sample_demagnetization_envelope.{suffix}",
                    dpi=300 if suffix == "png" else None)
    plt.close(fig)
    print(json.dumps({
        "mu0_Ms_T": MU0_MS_T,
        "public_plate_bound_T": plate_bound,
        "minimum_aspects": targets,
    }, indent=2))


if __name__ == "__main__":
    main()
