"""Externally validate public thermodynamic processing against direct Gd2B2MoO9 data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parent
RAW_ZIP = ROOT / "sources" / "raw" / "gd2b2moo9_source_data.zip"
WORKBOOK = ROOT / "sources" / "raw" / "gd2b2moo9_source_data" / "Data-all.xlsx"
PROVENANCE = ROOT / "sources" / "gd2b2moo9_external_validation_provenance_20261012.json"
RESULT = ROOT / "results" / "gd2b2moo9_external_validation.json"
FIGURE_STEM = ROOT / "figures" / "gd2b2moo9_external_validation"

FIELD_COLUMNS = {
    0.5: {"indirect_temperature": 19, "indirect_delta": 20, "direct_temperature": 38, "direct_delta": 39},
    1.0: {"indirect_temperature": 21, "indirect_delta": 22, "direct_temperature": 36, "direct_delta": 37},
    2.0: {"indirect_temperature": 25, "indirect_delta": 26, "direct_temperature": 34, "direct_delta": 35},
    5.0: {"indirect_temperature": 27, "indirect_delta": 28, "direct_temperature": 32, "direct_delta": 33},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_public_workbook() -> None:
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    RAW_ZIP.parent.mkdir(parents=True, exist_ok=True)
    request = Request(provenance["source_data_url"], headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, timeout=60) as response, RAW_ZIP.open("wb") as output:
        output.write(response.read())
    WORKBOOK.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(RAW_ZIP) as archive:
        archive.extract("Data-all.xlsx", WORKBOOK.parent)


def numeric_column(rows: list[tuple], one_based_column: int, first_row: int = 4) -> np.ndarray:
    values = []
    index = one_based_column - 1
    for row in rows[first_row - 1 :]:
        value = row[index] if index < len(row) else None
        if isinstance(value, (int, float)):
            values.append(float(value))
    return np.asarray(values, dtype=float)


def load_points() -> list[dict]:
    workbook = load_workbook(WORKBOOK, read_only=True, data_only=True)
    rows = list(workbook["Fig. 6 Direct measurement"].iter_rows(values_only=True))
    points = []
    for field, columns in FIELD_COLUMNS.items():
        indirect_temperature = numeric_column(rows, columns["indirect_temperature"])
        indirect_delta = numeric_column(rows, columns["indirect_delta"])
        count = min(indirect_temperature.size, indirect_delta.size)
        indirect_temperature = indirect_temperature[:count]
        indirect_delta = indirect_delta[:count]
        order = np.argsort(indirect_temperature)
        indirect_temperature = indirect_temperature[order]
        indirect_delta = indirect_delta[order]

        direct_temperature = numeric_column(rows, columns["direct_temperature"])
        direct_delta = numeric_column(rows, columns["direct_delta"])
        count = min(direct_temperature.size, direct_delta.size)
        direct_temperature = direct_temperature[:count]
        direct_delta = direct_delta[:count]
        predicted = np.interp(direct_temperature, indirect_temperature, indirect_delta)
        for temperature, measured, estimate in zip(direct_temperature, direct_delta, predicted):
            points.append(
                {
                    "field_T": field,
                    "starting_temperature_K": float(temperature),
                    "direct_delta_Tad_K": float(measured),
                    "indirect_interpolated_delta_Tad_K": float(estimate),
                    "error_indirect_minus_direct_K": float(estimate - measured),
                    "absolute_error_K": float(abs(estimate - measured)),
                }
            )
    return points


def summarize(points: list[dict]) -> dict:
    errors = np.asarray([point["error_indirect_minus_direct_K"] for point in points])
    measured = np.asarray([point["direct_delta_Tad_K"] for point in points])
    predicted = np.asarray([point["indirect_interpolated_delta_Tad_K"] for point in points])
    per_field = []
    for field in FIELD_COLUMNS:
        selected = [point for point in points if point["field_T"] == field]
        field_errors = np.asarray([point["error_indirect_minus_direct_K"] for point in selected])
        per_field.append(
            {
                "field_T": field,
                "point_count": len(selected),
                "mean_absolute_error_K": float(np.mean(np.abs(field_errors))),
                "maximum_absolute_error_K": float(np.max(np.abs(field_errors))),
                "mean_bias_K": float(np.mean(field_errors)),
            }
        )
    return {
        "status": "public second-material external thermodynamic validation",
        "source": json.loads(PROVENANCE.read_text(encoding="utf-8")),
        "workbook_sha256_observed": sha256(WORKBOOK),
        "point_count": len(points),
        "field_T": list(FIELD_COLUMNS),
        "mean_absolute_error_K": float(np.mean(np.abs(errors))),
        "root_mean_square_error_K": float(np.sqrt(np.mean(errors**2))),
        "maximum_absolute_error_K": float(np.max(np.abs(errors))),
        "mean_bias_K": float(np.mean(errors)),
        "pearson_correlation": float(np.corrcoef(measured, predicted)[0, 1]),
        "per_field": per_field,
        "points": points,
        "interpretation": [
            "Heat-capacity-derived Delta T_ad reproduces the source's twelve direct comparison points to a mean absolute error of about 0.014 K.",
            "This supports the thermodynamic-processing step on a second frustrated refrigerant with direct demagnetization data.",
            "It does not validate the NBCP critical-field transport surrogate, microscopic Hamiltonian, or package loss model.",
        ],
    }


def plot(payload: dict) -> None:
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), constrained_layout=True)
    colors = {0.5: "#3b82f6", 1.0: "#10b981", 2.0: "#f59e0b", 5.0: "#dc2626"}
    markers = {0.5: "o", 1.0: "s", 2.0: "^", 5.0: "D"}
    for field in FIELD_COLUMNS:
        points = [point for point in payload["points"] if point["field_T"] == field]
        direct = [point["direct_delta_Tad_K"] for point in points]
        predicted = [point["indirect_interpolated_delta_Tad_K"] for point in points]
        axes[0].scatter(
            direct,
            predicted,
            s=42,
            marker=markers[field],
            color=colors[field],
            edgecolor="white",
            linewidth=0.6,
            label=f"{field:g} T",
            zorder=3,
        )
        axes[1].scatter(
            [field] * len(points),
            [point["error_indirect_minus_direct_K"] * 1000 for point in points],
            s=42,
            marker=markers[field],
            color=colors[field],
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
    bounds = (-0.14, 2.82)
    axes[0].plot(bounds, bounds, color="0.25", lw=1.1, ls="--", label="identity")
    axes[0].set(xlim=bounds, ylim=bounds, xlabel="Direct $\\Delta T_{ad}$ (K)", ylabel="Indirect $\\Delta T_{ad}$ (K)")
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].legend(frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0.10, 1.0))
    axes[0].text(0.03, 0.96, "a", transform=axes[0].transAxes, fontweight="bold", va="top")
    axes[0].text(
        0.98,
        0.05,
        f"MAE = {payload['mean_absolute_error_K'] * 1000:.1f} mK\nmax = {payload['maximum_absolute_error_K'] * 1000:.1f} mK",
        transform=axes[0].transAxes,
        ha="right",
        va="bottom",
    )
    axes[1].axhline(0, color="0.25", lw=1.0, ls="--")
    axes[1].set(xlim=(0.2, 5.3), xticks=[0.5, 1, 2, 5], xticklabels=["0.5", "1", "2", "5"], xlabel="Field change (T)", ylabel="Indirect − direct (mK)")
    axes[1].set_ylim(-34, 16)
    axes[1].text(0.03, 0.96, "b", transform=axes[1].transAxes, fontweight="bold", va="top")
    axes[1].grid(axis="y", color="0.9", linewidth=0.7)
    fig.suptitle("External thermodynamic check with public Gd$_2$B$_2$MoO$_9$ data", fontsize=10.5)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURE_STEM.with_suffix(f".{suffix}"), dpi=300 if suffix == "png" else None)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="download the publisher source-data ZIP if absent")
    args = parser.parse_args()
    if not WORKBOOK.exists():
        if not args.download:
            raise FileNotFoundError("Public workbook missing; rerun with --download")
        fetch_public_workbook()
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    if sha256(WORKBOOK) != provenance["workbook_sha256"]:
        raise RuntimeError("Public workbook SHA-256 does not match declared provenance")
    payload = summarize(load_points())
    RESULT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    plot(payload)
    print(json.dumps({key: payload[key] for key in ["point_count", "mean_absolute_error_K", "root_mean_square_error_K", "maximum_absolute_error_K", "mean_bias_K", "pearson_correlation"]}, indent=2))


if __name__ == "__main__":
    main()
