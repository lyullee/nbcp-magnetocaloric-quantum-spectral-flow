"""Recover the published classical-MC C_m(T,B) vectors from Fig. 4d.

This is vector extraction, not manual pixel digitisation.  The PDF contains one
119-segment path for each 120-point classical-MC curve.  Axis coordinates are recovered
from the vector frame.  The output retains both PDF coordinates and SI data so
that every transform can be audited.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import fitz
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
PDF = ROOT / "sources" / "raw" / "arxiv_2202.05242_source" / "Fig" / "Fig4.pdf"
OUTDIR = ROOT / "data" / "derived"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Exact axis frame found in the vector drawing.  Fig. 4d spans T=0..0.30 K and
# C_m=0..7 J mol^-1 K^-1 (see the printed ticks and caption).
X_LEFT, X_RIGHT = 642.755, 1076.756
Y_TOP, Y_BOTTOM = 459.530, 786.031

# Draw order in the PDF and field labels in the printed legend.
FIELD_BY_PATH_ORDER = [0.00, 1.33, 1.07, 0.80, 0.53, 0.27, 0.13, 1.60]


def data_paths(page):
    paths = []
    for d in page.get_drawings():
        r = d["rect"]
        if (
            len(d["items"]) == 119
            and abs(r.x0 - 655.432) < 0.05
            and abs(r.x1 - 1032.292) < 0.05
            and d.get("color") is not None
        ):
            paths.append(d)
    if len(paths) != 8:
        raise RuntimeError(f"expected 8 classical-MC paths; found {len(paths)}")
    return paths


def path_points(path):
    items = path["items"]
    points = [items[0][1]] + [item[2] for item in items]
    if len(points) != 120:
        raise RuntimeError("curve does not contain 120 points")
    return points


def map_point(p):
    t = (p.x - X_LEFT) / (X_RIGHT - X_LEFT) * 0.30
    c = (Y_BOTTOM - p.y) / (Y_BOTTOM - Y_TOP) * 7.0
    return t, c


def main():
    doc = fitz.open(PDF)
    page = doc[0]
    rows = []
    curves = {}
    for field, path in zip(FIELD_BY_PATH_ORDER, data_paths(page)):
        points = path_points(path)
        curve = []
        for i, p in enumerate(points):
            t, c = map_point(p)
            row = {
                "field_T": field,
                "point_index": i,
                "temperature_K": t,
                "Cm_J_mol_K": c,
                "pdf_x_pt": p.x,
                "pdf_y_pt": p.y,
            }
            rows.append(row)
            curve.append(row)
        curves[str(field)] = curve

    csv_path = OUTDIR / "npj2022_fig4d_classical_mc_Cm_vector.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    sha = hashlib.sha256(PDF.read_bytes()).hexdigest().upper()
    meta = {
        "source_article": "Gao et al., npj Quantum Materials 7, 89 (2022)",
        "source_url": "https://doi.org/10.1038/s41535-022-00500-3",
        "source_arxiv": "https://arxiv.org/abs/2202.05242",
        "source_pdf_relative_path": str(PDF.relative_to(ROOT)),
        "source_pdf_sha256": sha,
        "extraction": "PyMuPDF vector path extraction; no raster/pixel tracing",
        "method": "classical Monte Carlo on up to 48x48 spins; this is not XTRG quantum data",
        "figure_panel": "Fig. 4d",
        "axis_map": {
            "pdf_frame_points": [X_LEFT, Y_TOP, X_RIGHT, Y_BOTTOM],
            "temperature_K": [0.0, 0.30],
            "Cm_J_mol_K": [0.0, 7.0],
        },
        "fields_T_in_pdf_draw_order": FIELD_BY_PATH_ORDER,
        "points_per_curve": 120,
        "output_csv": str(csv_path.relative_to(ROOT)),
        "model_in_source": {
            "Jxy_K": 0.88,
            "Jz_K": 1.48,
            "J_Gamma_and_J_PD": "negligible / fixed to zero",
            "g_c": 4.89,
        },
        "current_M0_model": {
            "Jxy_K_approx": 0.1225 * 0.636 * 11.60451812,
            "Jz_K_approx": 0.1225 * 11.60451812,
            "g_c": 4.716,
        },
        "limitations": [
            "The source curves begin near 8.76 mK; entropy below the first point requires an explicit low-T assumption.",
            "Only 1.07 T and 1.33 T bracket the current 1.10--1.25 T window.",
            "The source and current M0 parameters differ by several percent and must not be silently merged.",
        ],
    }
    (OUTDIR / "npj2022_fig4d_classical_mc_Cm_vector_metadata.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for field in sorted(FIELD_BY_PATH_ORDER):
        q = curves[str(field)]
        ax.plot([r["temperature_K"] for r in q], [r["Cm_J_mol_K"] for r in q], label=f"{field:.2f} T")
    ax.set(xlim=(0, 0.30), ylim=(0, 7), xlabel="Temperature (K)", ylabel=r"$C_m$ (J mol$^{-1}$ K$^{-1}$)")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ax.set_title("Vector-recovered classical-MC curves from npj QM (2022), Fig. 4d")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "npj2022_fig4d_classical_mc_vector_recovery.png", dpi=220)
    plt.close(fig)

    print(f"wrote {len(rows)} points; T=[{min(r['temperature_K'] for r in rows):.6f}, {max(r['temperature_K'] for r in rows):.6f}] K")
    print(f"source sha256 {sha}")


if __name__ == "__main__":
    main()
