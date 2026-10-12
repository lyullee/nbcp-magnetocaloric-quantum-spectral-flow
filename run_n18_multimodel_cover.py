"""Run the preregistered N=18 interlayer-model cover.

This runner extends the already validated M3/M4-01 calculation to three
additional M4 Hamiltonians chosen *before* inspecting their N=18 results:

* M4-00: minimum retained bulk critical-field shift;
* M4-02: nominal/local-Cartesian representative;
* M4-03: median retained bulk shift with a distinct coupling composition;
* M4-01: maximum retained bulk critical-field shift (already calculated).

All raw spectra use the same projected-magnon settings and cache format as
``audit_n18_stacked_m4_convergence.py``.  Stages are resumable and never
overwrite the publication summary produced by the compiler.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path

import numpy as np

import audit_n18_stacked_m4_convergence as base


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
PLAN = RESULTS / "n18_multimodel_preregistered_plan.json"
PILOT = RESULTS / "n18_multimodel_coarse_checkpoint.json"
REFINE = RESULTS / "n18_multimodel_refine_checkpoint.json"
CONVERGENCE = RESULTS / "n18_multimodel_convergence_checkpoint.json"
AMENDMENT = RESULTS / "n18_multimodel_plan_amendment_after_coarse.json"
CONVERGENCE_PROTOCOL = RESULTS / "n18_multimodel_convergence_protocol.json"
CONVERGENCE_AMENDMENT = RESULTS / "n18_multimodel_convergence_protocol_amendment.json"
CANDIDATE_GAP_KB_PER_SPIN = 0.02


MODEL_SPECS = {
    "M3": {
        "parameters": base.InterlayerParameters(0.0, 0.0, 0.0, 0.08),
        "sample_id": "M3",
        "selection_role": "common no-interlayer baseline",
        "bulk_shift_T": 0.0,
        "n12_shift_T": 0.0,
        "n12_heat_ratio": 1.0,
    },
    "M4-00": {
        "parameters": base.InterlayerParameters(
            0.008018729090690612, -2.456390716508031,
            5.570390689605847, 0.12898443758487702),
        "sample_id": "761",
        "selection_role": "minimum retained bulk critical-field shift",
        "bulk_shift_T": 0.033403908987329034,
        "n12_shift_T": 0.0,
        "n12_heat_ratio": 1.0127213300003153,
    },
    "M4-02": {
        "parameters": base.InterlayerParameters(2.8, 0.0, 7.9, 0.08),
        "sample_id": "188",
        "selection_role": "nominal local-Cartesian representative",
        "bulk_shift_T": 0.09527380895522186,
        "n12_shift_T": 0.040000000000000036,
        "n12_heat_ratio": 0.9943352822895036,
    },
    "M4-03": {
        "parameters": base.InterlayerParameters(
            0.27222977876663207, -2.8672347117215393,
            13.435032871365546, 0.04629844818264246),
        "sample_id": "738",
        "selection_role": (
            "median retained bulk shift with a coupling composition distinct "
            "from the nominal and maximum-shift anchors"
        ),
        "bulk_shift_T": 0.11555865530097242,
        "n12_shift_T": 0.050000000000000044,
        "n12_heat_ratio": 0.9607978595680514,
    },
    "M4-01": {
        "parameters": base.InterlayerParameters(
            4.5861779235303395, 3.4409039821475744,
            13.63851090040989, 0.03836957771331072),
        "sample_id": "773",
        "selection_role": "maximum retained bulk critical-field shift",
        "bulk_shift_T": 0.20210652482768032,
        "n12_shift_T": 0.10000000000000009,
        "n12_heat_ratio": 0.9192300000169081,
    },
}


def configure_models() -> None:
    base.MODELS = {
        label: spec["parameters"] for label, spec in MODEL_SPECS.items()
    }


def write_plan() -> None:
    if PLAN.exists():
        return
    retained = RESULTS / "interlayer_adaptive_inverse_uq_retained.csv"
    n12 = RESULTS / "m4_finite_temperature_cluster_comparisons.csv"
    payload = {
        "status": "frozen before new N=18 M4-00/M4-02/M4-03 calculations",
        "cluster": {
            "shape": [3, 3, 2],
            "spins": 18,
            "boundary": "tin_foil",
            "three_sublattice_compatible": True,
        },
        "selection_rule": (
            "M3 plus minimum retained bulk shift, nominal representative, "
            "median retained bulk shift with distinct coupling composition, "
            "and maximum retained bulk shift"
        ),
        "models": {
            label: {
                key: value for key, value in spec.items() if key != "parameters"
            } | {
                "parameters": {
                    "j2_microev": spec["parameters"].j2_microev,
                    "ja_microev": spec["parameters"].ja_microev,
                    "jb_microev": spec["parameters"].jb_microev,
                    "delta_c": spec["parameters"].delta_c,
                }
            }
            for label, spec in MODEL_SPECS.items()
        },
        "numerical_protocol": {
            "coarse": {
                "fields_T": np.round(np.arange(1.10, 2.0001, 0.10), 10).tolist(),
                "even_odd_cutoffs": [4, 5],
                "states_per_parity": 8,
            },
            "refine": {
                "field_window_T": "coarse peak +/- 0.08 in 0.02 T steps",
                "even_odd_cutoffs": [6, 7],
                "states_per_parity": 12,
            },
            "peak_and_cutoff": {
                "three_point_window_T": "refine peak +/- 0.02",
                "states_per_parity": [12, 32],
                "even_odd_cutoffs": [[6, 7], [8, 9]],
            },
            "hot_reference": {
                "field_T": 4.0,
                "temperature_K": 2.0,
                "even_odd_cutoffs": [[6, 7], [8, 9], [10, 11]],
            },
        },
        "source_files": [
            str(retained.relative_to(ROOT)), str(n12.relative_to(ROOT))
        ],
        "interpretation_rule": (
            "Report the complete ordered cover irrespective of monotonicity; "
            "do not promote N=18 to a thermodynamic-limit prediction."
        ),
    }
    PLAN.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def selected_labels(text: str | None) -> list[str]:
    if not text:
        return list(MODEL_SPECS)
    labels = [item.strip() for item in text.split(",") if item.strip()]
    unknown = [label for label in labels if label not in MODEL_SPECS]
    if unknown:
        raise ValueError(f"unknown model labels: {unknown}")
    return labels


def checkpoint(path: Path, stage: str, rows: list[dict], labels: list[str]) -> None:
    path.write_text(json.dumps({
        "status": f"completed {stage} checkpoint",
        "models": labels,
        "rows": rows,
        "summary": base.summarize(rows),
    }, indent=2), encoding="utf-8")


def stage_path(path: Path, labels: list[str]) -> Path:
    """Keep subset-worker checkpoints separate while sharing raw caches."""
    if labels == list(MODEL_SPECS):
        return path
    suffix = "_".join(label.replace("-", "") for label in labels)
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")


def write_coarse_amendment(pilot_path: Path) -> None:
    if AMENDMENT.exists():
        return
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    candidates = {}
    for model in pilot["models"]:
        selected = sorted(
            [row for row in pilot["rows"] if row["model"] == model],
            key=lambda row: row["field_T"])
        values = np.array([row["entropy_kB_per_spin"] for row in selected])
        threshold = 0.75 * float(values.max())
        local = []
        for index in range(1, len(selected) - 1):
            if (values[index] >= values[index - 1]
                    and values[index] >= values[index + 1]
                    and values[index] >= threshold):
                local.append({
                    "field_T": selected[index]["field_T"],
                    "entropy_kB_per_spin": selected[index]["entropy_kB_per_spin"],
                })
        candidates[model] = local
    AMENDMENT.write_text(json.dumps({
        "status": "rule amended after the declared low-accuracy coarse pilot and before any new K6/K7 dense-cover calculation",
        "trigger": (
            "The K4/K5 pilot showed multiple separated local maxima within "
            "75% of the largest entropy for several models. A single +/-0.08 T "
            "refinement window could therefore select the wrong global peak."
        ),
        "amended_refine_rule": {
            "fields_T": np.round(np.arange(1.10, 2.1001, 0.02), 10).tolist(),
            "even_odd_cutoffs": [6, 7],
            "states_per_parity": 12,
            "selection": "global maximum on the common dense field grid",
        },
        "coarse_local_maxima_above_75_percent": candidates,
        "interpretation_rule_unchanged": (
            "Report the complete ordered cover irrespective of monotonicity; "
            "do not promote N=18 to a thermodynamic-limit prediction."
        ),
    }, indent=2), encoding="utf-8")


def write_convergence_protocol() -> None:
    if CONVERGENCE_PROTOCOL.exists():
        return
    CONVERGENCE_PROTOCOL.write_text(json.dumps({
        "status": "frozen before completion of the new K6/K7 dense-cover calculations",
        "trigger": (
            "The declared K4/K5 pilot showed separated peaks. Increasing the "
            "low-state count could reorder close K6/K7 dense-grid maxima."
        ),
        "candidate_rule": (
            "Retain every interior K6/K7-L12 local maximum whose entropy is "
            "within 0.01 kB/spin of the dense-grid maximum."
        ),
        "selection_rule": (
            "Evaluate K6/K7-L32 at each candidate and its +/-0.02 T neighbors; "
            "select the candidate with the largest local quadratic peak entropy."
        ),
        "final_crosscheck": (
            "At the selected grid center, compare K6/K7 and K8/K9 with both "
            "12 and 32 low states per parity."
        ),
        "scope": (
            "This selects the global peak inside the declared N=18 projection; "
            "it does not establish thermodynamic-limit convergence."
        ),
    }, indent=2), encoding="utf-8")


def write_convergence_amendment() -> None:
    """Conservatively widen candidate inclusion before the dense scan finishes."""
    if CONVERGENCE_AMENDMENT.exists():
        return
    m3_corrections = []
    for field in (1.16, 1.18, 1.20):
        l12_path = base.spectrum_path("M3", field, 6, 7, 12)
        l32_path = base.spectrum_path("M3", field, 6, 7, 32)
        if l12_path.exists() and l32_path.exists():
            l12 = json.loads(l12_path.read_text(encoding="utf-8"))
            l32 = json.loads(l32_path.read_text(encoding="utf-8"))
            m3_corrections.append({
                "field_T": field,
                "L12_entropy_kB_per_spin": l12["entropy_kB_per_spin"],
                "L32_entropy_kB_per_spin": l32["entropy_kB_per_spin"],
                "absolute_change_kB_per_spin": abs(
                    l32["entropy_kB_per_spin"] - l12["entropy_kB_per_spin"]),
            })
    partial_m403 = {}
    for field in (1.30, 1.54):
        path = base.spectrum_path("M4-03", field, 6, 7, 12)
        if path.exists():
            partial_m403[f"{field:.2f}_T_entropy_kB_per_spin"] = json.loads(
                path.read_text(encoding="utf-8"))["entropy_kB_per_spin"]
    if len(partial_m403) == 2:
        values = list(partial_m403.values())
        partial_m403["absolute_gap_kB_per_spin"] = abs(values[1] - values[0])
    CONVERGENCE_AMENDMENT.write_text(json.dumps({
        "status": (
            "frozen before completion of the five K6/K7 dense scans; conservative "
            "candidate expansion only"
        ),
        "created_local": datetime.now().astimezone().isoformat(),
        "trigger": (
            "Existing M3 L12-to-L32 entropy corrections reach approximately 0.0102 "
            "kB/spin, comparable to the original 0.01 candidate gap. During the "
            "still-incomplete common-grid scan, the M4-03 local peaks at 1.30 and "
            "1.54 T differed by approximately 0.01035 kB/spin. The original rule "
            "could therefore discard a branch within the observed state-count error."
        ),
        "original_candidate_gap_kB_per_spin": 0.01,
        "amended_candidate_gap_kB_per_spin": CANDIDATE_GAP_KB_PER_SPIN,
        "observed_precompletion_evidence": {
            "M3_existing_state_count_corrections": m3_corrections,
            "M4-03_partial_dense_peak_gap": partial_m403,
        },
        "unchanged": [
            "ordered model cover", "1.10-2.10 T common grid", "0.02 T spacing",
            "L32 three-point candidate calculation", "largest local quadratic peak selection",
            "reported metrics", "thermodynamic-limit scope guard",
        ],
        "bias_guard": (
            "The amendment only adds candidates and was recorded before any of the "
            "five 51-point dense scans completed."
        ),
    }, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage", choices=("coarse", "refine", "convergence", "all"),
        default="all")
    parser.add_argument(
        "--models", help="comma-separated subset; defaults to the frozen cover")
    args = parser.parse_args()

    configure_models()
    write_plan()
    write_convergence_protocol()
    write_convergence_amendment()
    labels = selected_labels(args.models)
    base.MODELS = {label: MODEL_SPECS[label]["parameters"] for label in labels}
    pilot_path = stage_path(PILOT, labels)
    refine_path = stage_path(REFINE, labels)
    convergence_path = stage_path(CONVERGENCE, labels)
    graph, common = base.load_common()
    if graph["N"] != 18 or tuple(graph["shape"]) != (3, 3, 2):
        raise RuntimeError("common-bond cache is not the 3x3x2 N=18 cluster")

    def coarse() -> list[dict]:
        rows = []
        fields = np.round(np.arange(1.10, 2.0001, 0.10), 10)
        for model in labels:
            for index, field in enumerate(fields):
                row = base.cold_spectrum(
                    graph, common, model, float(field), 4, 5, 8)
                rows.append(row)
                print(json.dumps({
                    "stage": "coarse", "model": model,
                    "completed": index + 1, "total": len(fields),
                    "field_T": float(field),
                    "entropy_kB_per_spin": row["entropy_kB_per_spin"],
                    "elapsed_seconds": row["elapsed_seconds"],
                }), flush=True)
        checkpoint(pilot_path, "coarse", rows, labels)
        if labels == list(MODEL_SPECS):
            write_coarse_amendment(pilot_path)
        return rows

    def refine() -> list[dict]:
        # All subset workers use the same completed full-cover coarse pilot.
        source_pilot = PILOT if PILOT.exists() else pilot_path
        if not source_pilot.exists():
            coarse()
        rows = []
        fields = np.round(np.arange(1.10, 2.1001, 0.02), 10)
        for model in labels:
            for index, field in enumerate(fields):
                row = base.cold_spectrum(
                    graph, common, model, float(field), 6, 7, 12)
                rows.append(row)
                print(json.dumps({
                    "stage": "refine", "model": model,
                    "completed": index + 1, "total": len(fields),
                    "field_T": float(field),
                    "entropy_kB_per_spin": row["entropy_kB_per_spin"],
                    "elapsed_seconds": row["elapsed_seconds"],
                }), flush=True)
        checkpoint(refine_path, "dense-refine", rows, labels)
        return rows

    def convergence() -> None:
        if not refine_path.exists():
            refine()
        refined = json.loads(refine_path.read_text(encoding="utf-8"))
        output = {}
        for model in labels:
            dense = sorted(
                [row for row in refined["rows"] if row["model"] == model],
                key=lambda row: row["field_T"])
            fields = np.array([row["field_T"] for row in dense], dtype=float)
            entropy = np.array(
                [row["entropy_kB_per_spin"] for row in dense], dtype=float)
            local_indices = [
                index for index in range(1, len(dense) - 1)
                if entropy[index] >= entropy[index - 1]
                and entropy[index] >= entropy[index + 1]
                and entropy[index] >= float(entropy.max()) - CANDIDATE_GAP_KB_PER_SPIN
            ]
            if not local_indices:
                local_indices = [int(np.argmax(entropy))]
            if any(index in (0, len(dense) - 1) for index in local_indices):
                raise RuntimeError(f"{model} candidate peak lies on dense-grid boundary")

            cold_rows = []
            candidates = []
            for candidate_index in local_indices:
                center = float(fields[candidate_index])
                neighbors = [round(center - 0.02, 2), center,
                             round(center + 0.02, 2)]
                candidate_rows = [
                    base.cold_spectrum(
                        graph, common, model, field, 6, 7, 32)
                    for field in neighbors
                ]
                cold_rows.extend(candidate_rows)
                values = np.array(
                    [row["entropy_kB_per_spin"] for row in candidate_rows],
                    dtype=float)
                a, b, c = np.polyfit(neighbors, values, 2)
                location = -b / (2 * a) if a < 0 else center
                if a < 0 and neighbors[0] <= location <= neighbors[-1]:
                    peak_entropy = float(a * location**2 + b * location + c)
                    peak_field = float(location)
                else:
                    best = int(np.argmax(values))
                    peak_entropy = float(values[best])
                    peak_field = float(neighbors[best])
                candidates.append({
                    "grid_center_T": center,
                    "neighbor_fields_T": neighbors,
                    "L32_entropy_kB_per_spin": values.tolist(),
                    "local_peak_field_T": peak_field,
                    "local_peak_entropy_kB_per_spin": peak_entropy,
                    "quadratic_concave_and_bracketed": bool(
                        a < 0 and neighbors[0] <= location <= neighbors[-1]),
                })
            selected = max(
                candidates, key=lambda row: row["local_peak_entropy_kB_per_spin"])
            center = selected["grid_center_T"]
            cold_rows.append(base.cold_spectrum(
                graph, common, model, center, 8, 9, 12))
            cold_rows.append(base.cold_spectrum(
                graph, common, model, center, 8, 9, 32))
            hot_rows = {
                f"K{even}-{odd}": base.hot_entropy(
                    graph, common, model, even, odd)
                for even, odd in ((6, 7), (8, 9), (10, 11))
            }
            output[model] = {
                "dense_L12_peak_T": float(fields[int(np.argmax(entropy))]),
                "candidate_entropy_gap_kB_per_spin": CANDIDATE_GAP_KB_PER_SPIN,
                "candidate_peaks": candidates,
                "selected_peak_grid_center_T": center,
                "cold_files": [
                    base.spectrum_path(
                        model, row["field_T"], row["even_cutoff"],
                        row["odd_cutoff"], row["low_states_per_parity"]
                    ).name for row in cold_rows
                ],
                "hot_entropy_kB_per_spin": {
                    key: value["entropy_kB_per_spin"]
                    for key, value in hot_rows.items()
                },
            }
            print(json.dumps({
                "stage": "convergence", "model": model,
                "dense_L12_peak_T": output[model]["dense_L12_peak_T"],
                "selected_L32_center_T": center,
                "candidate_count": len(candidates),
                "hot_entropy_kB_per_spin": output[model]["hot_entropy_kB_per_spin"],
            }), flush=True)
        convergence_path.write_text(json.dumps({
            "status": "completed preregistered N=18 multimodel convergence runs",
            "models": output,
        }, indent=2), encoding="utf-8")

    if args.stage == "coarse":
        coarse()
    elif args.stage == "refine":
        refine()
    elif args.stage == "convergence":
        convergence()
    else:
        coarse()
        refine()
        convergence()


if __name__ == "__main__":
    main()
