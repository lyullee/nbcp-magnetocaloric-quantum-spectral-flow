"""Calibration-product noise diagnostic for compiled N=9 state circuits.

This is not a density-matrix simulation.  It multiplies FakeSherbrooke's
per-instruction error complements along each compiled circuit and compares the
scheduled duration with T1/T2.  The result is a conservative gate-screening
diagnostic before investing in a full noise/mitigation study.
"""
from __future__ import annotations

import json
from pathlib import Path
import time
import warnings

import numpy as np
from qiskit import transpile
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from fit_number_conserving_ansatz_n9 import build_qiskit_circuit
from screen_model import lattice
from screen_weak_terms import directional_edges


ROOT = Path(__file__).resolve().parent
SEED = 20261017


def operations(n, edges):
    result = []
    for layer in range(2):
        ordered = edges if layer == 0 else list(reversed(edges))
        result.extend(("rz", site) for site in range(n))
        result.extend(("givens", i, j) for i, j in ordered)
    result.extend(("rz", site) for site in range(n))
    result.extend(("pair", i, j) for i, j in edges)
    return result


def circuit_error_product(circuit, backend):
    log_survival = 0.0
    log_two_qubit = 0.0
    expected_errors = 0.0
    expected_two_qubit_errors = 0.0
    used_qubits = set()
    counted = 0
    counted_two = 0
    missing = []
    for instruction in circuit.data:
        name = instruction.operation.name
        if name in {"barrier", "delay"}:
            continue
        qargs = tuple(circuit.find_bit(qubit).index
                      for qubit in instruction.qubits)
        used_qubits.update(qargs)
        properties = backend.target[name].get(qargs)
        error = None if properties is None else properties.error
        if error is None:
            missing.append((name, qargs))
            continue
        error = min(max(float(error), 0.0), 1.0 - 1e-15)
        log_survival += np.log1p(-error)
        expected_errors += error
        counted += 1
        if len(qargs) == 2:
            log_two_qubit += np.log1p(-error)
            expected_two_qubit_errors += error
            counted_two += 1
    return {
        "calibration_product_all_gates": float(np.exp(log_survival)),
        "calibration_product_two_qubit_gates": float(np.exp(log_two_qubit)),
        "sum_reported_gate_errors": expected_errors,
        "sum_reported_two_qubit_errors": expected_two_qubit_errors,
        "counted_gates": counted,
        "counted_two_qubit_gates": counted_two,
        "used_physical_qubits": sorted(used_qubits),
        "missing_error_entries": [
            {"instruction": name, "qubits": list(qargs)}
            for name, qargs in missing
        ],
    }


def main():
    started = time.perf_counter()
    source = json.loads((ROOT / "results" /
                         "spectral_flow_stateprep_N9_B1.45T.json").read_text(
                             encoding="utf-8"))
    graph = lattice((3, 0), (0, 3))
    n = graph["N"]
    edges = [(i, j) for i, j, _ in directional_edges(graph)]
    ansatz_operations = operations(n, edges)
    backend = FakeSherbrooke()
    rows = []
    for node in source["nodes"]:
        for representative in node["representatives"]:
            circuit = build_qiskit_circuit(
                n, representative["initial_computational_basis_index"],
                ansatz_operations,
                np.asarray(representative["parameters_radians"], dtype=float))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                compiled = transpile(
                    circuit, backend=backend, optimization_level=2,
                    seed_transpiler=SEED, scheduling_method="alap")
                duration_dt = int(compiled.duration)
            audit = circuit_error_product(compiled, backend)
            qubit_properties = [backend.qubit_properties(qubit)
                                for qubit in audit["used_physical_qubits"]]
            t1 = np.asarray([item.t1 for item in qubit_properties])
            t2 = np.asarray([item.t2 for item in qubit_properties])
            duration_seconds = duration_dt * backend.dt
            row = {
                "lambda": node["lambda"],
                "representative_level":
                    representative["representative_level"],
                "compiled_ecr_count": int(
                    compiled.count_ops().get("ecr", 0)),
                "compiled_depth": int(compiled.depth()),
                "scheduled_duration_dt": duration_dt,
                "scheduled_duration_seconds": duration_seconds,
                "median_used_T1_seconds": float(np.median(t1)),
                "median_used_T2_seconds": float(np.median(t2)),
                "duration_over_median_T1": float(
                    duration_seconds / np.median(t1)),
                "duration_over_median_T2": float(
                    duration_seconds / np.median(t2)),
                **audit,
            }
            rows.append(row)
            print(f"lambda={node['lambda']:.6f} L={representative['representative_level']}: "
                  f"survival={audit['calibration_product_all_gates']:.3e}, "
                  f"duration/T1={row['duration_over_median_T1']:.2f}",
                  flush=True)

    payload = {
        "purpose": "topology/calibration noise pre-screen for N=9 spectral-flow state preparation",
        "backend_model": backend.name,
        "source": "results/spectral_flow_stateprep_N9_B1.45T.json",
        "rows": rows,
        "decision": "STOP_RAW_QPU_CONTINUE_COMPRESSION_AND_MITIGATION",
        "decision_basis": {
            "minimum_calibration_product_all_gates": min(
                row["calibration_product_all_gates"] for row in rows),
            "maximum_scheduled_duration_over_median_T1": max(
                row["duration_over_median_T1"] for row in rows),
            "maximum_scheduled_duration_over_median_T2": max(
                row["duration_over_median_T2"] for row in rows),
        },
        "scope_limits": [
            "The product of reported gate fidelities is a heuristic probability of no calibrated gate error, not circuit fidelity.",
            "Correlated, coherent, leakage and time-dependent errors are not modeled.",
            "No measurement/readout circuits are included.",
            "FakeSherbrooke is a fixed historical calibration model and does not represent the live Open Plan backend.",
            "The result supports rejecting an unmitigated raw run; it does not predict the outcome after circuit compression or error mitigation.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "fake_backend_noise_audit_N9.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "decision": payload["decision"],
        **payload["decision_basis"],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
