"""Recover the XTRG C_m/T heatmap in Fig. 4a from the embedded PDF images.

The panel and its color bar are embedded as lossless indexed images.  Vector
axes provide the coordinate calibration: B=0..2 T and a logarithmic T axis
with the 1 K and 0.1 K ticks at known PDF coordinates.  RGB values are mapped
to the separately embedded 0..15 J mol^-1 K^-2 color bar.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import fitz
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
PDF = ROOT / "sources" / "raw" / "arxiv_2202.05242_source" / "Fig" / "Fig4.pdf"
OUTDIR = ROOT / "data" / "derived"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Vector-axis calibration in PDF points.
X_B0, X_B2 = 79.335, 464.835
Y_T1, Y_T01 = 34.335, 315.521


def rgb_pixmap(doc, xref):
    p = fitz.Pixmap(doc, xref)
    p = fitz.Pixmap(fitz.csRGB, p)
    return np.frombuffer(p.samples, dtype=np.uint8).reshape(p.height, p.width, 3)


def page_to_pixel(x, y, rect, width, height):
    u = (x - rect.x0) / rect.width * (width - 1)
    v = (y - rect.y0) / rect.height * (height - 1)
    return u, v


def bilinear(img, u, v):
    u = float(np.clip(u, 0, img.shape[1] - 1))
    v = float(np.clip(v, 0, img.shape[0] - 1))
    x0, y0 = int(np.floor(u)), int(np.floor(v))
    x1, y1 = min(x0 + 1, img.shape[1] - 1), min(y0 + 1, img.shape[0] - 1)
    a, b = u - x0, v - y0
    return (
        (1 - a) * (1 - b) * img[y0, x0]
        + a * (1 - b) * img[y0, x1]
        + (1 - a) * b * img[y1, x0]
        + a * b * img[y1, x1]
    )


def main():
    doc = fitz.open(PDF)
    page = doc[0]
    images = page.get_images(full=True)
    panel = next(x for x in images if x[2] > 2000 and x[3] > 1500)
    bar = next(x for x in images if x[2] == 32 and x[3] > 600)
    panel_img = rgb_pixmap(doc, panel[0])
    bar_img = rgb_pixmap(doc, bar[0])
    panel_rect = page.get_image_rects(panel[0])[0]

    # Median over the bar width suppresses subpixel edge artefacts. Top=15,
    # bottom=0 as printed beside the bar.
    bar_rgb = np.median(bar_img[:, 4:-4, :], axis=1)
    bar_values = np.linspace(15.0, 0.0, len(bar_rgb))

    def sample(b, t):
        x = X_B0 + (b / 2.0) * (X_B2 - X_B0)
        y = Y_T1 + (-np.log10(t)) * (Y_T01 - Y_T1)
        u, v = page_to_pixel(x, y, panel_rect, panel_img.shape[1], panel_img.shape[0])
        rgb = bilinear(panel_img, u, v)
        distance = np.sum((bar_rgb - rgb) ** 2, axis=1)
        k = int(np.argmin(distance))
        return float(bar_values[k]), rgb, float(np.sqrt(distance[k]))

    fields = np.linspace(0.0, 2.0, 201)
    temperatures = np.geomspace(0.07, 1.0, 181)
    rows = []
    for t in temperatures:
        for b in fields:
            value, rgb, residual = sample(float(b), float(t))
            rows.append(
                {
                    "field_T": b,
                    "temperature_K": t,
                    "Cm_over_T_J_mol_K2": value,
                    "rgb_R": rgb[0],
                    "rgb_G": rgb[1],
                    "rgb_B": rgb[2],
                    "colorbar_RGB_residual": residual,
                }
            )

    path = OUTDIR / "npj2022_fig4a_xtrg_Cm_over_T_heatmap.csv"
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    meta = {
        "source_article": "Gao et al., npj Quantum Materials 7, 89 (2022)",
        "source_url": "https://doi.org/10.1038/s41535-022-00500-3",
        "source_pdf_sha256": hashlib.sha256(PDF.read_bytes()).hexdigest().upper(),
        "panel": "Fig. 4a",
        "quantity": "XTRG magnetic specific heat divided by temperature",
        "units": "J mol^-1 K^-2",
        "xtrg_geometry": "YC6x9, bond dimension up to D=400 (article Methods)",
        "axis_calibration_pdf_points": {"B0_x": X_B0, "B2_x": X_B2, "T1_y": Y_T1, "T0.1_y": Y_T01},
        "embedded_images": {
            "panel_xref": panel[0], "panel_shape": list(panel_img.shape),
            "colorbar_xref": bar[0], "colorbar_shape": list(bar_img.shape),
        },
        "grid": {"field_T": [0.0, 2.0, 0.01], "temperature_K_logspace": [0.07, 1.0, 181]},
        "limitations": [
            "Values are recovered from the lossless plotted heatmap and color scale, not the authors' numerical table.",
            "Color quantisation is about 15/653 = 0.023 J mol^-1 K^-2 before interpolation effects.",
            "The paper reports XTRG truncation errors below about 100 mK and finite-cylinder effects at low temperature.",
        ],
        "output_csv": str(path.relative_to(ROOT)),
    }
    (OUTDIR / "npj2022_fig4a_xtrg_Cm_over_T_heatmap_metadata.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    z = np.array([r["Cm_over_T_J_mol_K2"] for r in rows]).reshape(len(temperatures), len(fields))
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.pcolormesh(fields, temperatures, z, shading="auto", vmin=0, vmax=15, cmap="viridis")
    ax.set_yscale("log")
    ax.set(xlabel="Field B (T)", ylabel="Temperature (K)", title="Recovered XTRG heatmap, npj QM (2022) Fig. 4a")
    fig.colorbar(im, ax=ax, label=r"$C_m/T$ (J mol$^{-1}$ K$^{-2}$)")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "npj2022_fig4a_xtrg_heatmap_recovered.png", dpi=220)
    plt.close(fig)

    residuals = np.array([r["colorbar_RGB_residual"] for r in rows])
    print(f"wrote {len(rows)} heatmap samples; median RGB residual={np.median(residuals):.3f}, p99={np.quantile(residuals,.99):.3f}")


if __name__ == "__main__":
    main()
