"""Read-only IBM Open Plan preflight for adaptive N=6 spectral flow.

This script authenticates with the saved Open-only account, compiles all nine
frozen circuits against every operational free-plan QPU, and records only
sanitized usage and calibration diagnostics.  It never submits a job.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import time
import warnings

import numpy as np
from qiskit import transpile
from qiskit_ibm_runtime import QiskitRuntimeService

from audit_fake_backend_noise_n9 import circuit_error_product
from fit_number_conserving_ansatz_n9 import build_qiskit_circuit
from qiskit_pauli_mapping import (
    open_triangular_directional_edges,
    open_weak_model_bonds,
    spin_hamiltonian_to_sparse_pauli,
)


ROOT = Path(__file__).resolve().parent
ACCOUNT_NAME = "psa_open"
SEED = 20261013


def safe_usage(raw):
    keep = ("usage_consumed_seconds", "usage_consumed_ru",
            "usage_limit_seconds", "usage_limit_ru", "usage_limit_reached",
            "usage_remaining_seconds", "usage_period")
    return {key: raw.get(key) for key in keep if key in raw}


def all_circuits(source):
    n = int(source["condition"]["N"])
    for node in source["nodes"]:
        for row in node["states"]:
            operations = [tuple(item) for item in row["operations"]]
            parameters = np.asarray(row["parameters_radians"], dtype=float)
            circuit = build_qiskit_circuit(
                n, int(row["initial_basis_index"]), operations, parameters)
            circuit.metadata = {"lambda": float(row["lambda"]),
                                "level": int(row["level"])}
            yield row, circuit


def compile_backend(backend, source):
    rows = []
    for source_row, circuit in all_circuits(source):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            compiled = transpile(
                circuit, backend=backend, optimization_level=3,
                seed_transpiler=SEED, scheduling_method="alap")
        audit = circuit_error_product(compiled, backend)
        used = audit["used_physical_qubits"]
        properties = [backend.qubit_properties(qubit) for qubit in used]
        duration_seconds = int(compiled.duration) * float(backend.dt)
        median_t1 = float(np.median([item.t1 for item in properties]))
        median_t2 = float(np.median([item.t2 for item in properties]))
        two_qubit_count = sum(
            int(count) for gate, count in compiled.count_ops().items()
            if gate in {"cx", "cz", "ecr", "iswap", "rzz"})
        rows.append({
            "lambda": float(source_row["lambda"]),
            "level": int(source_row["level"]),
            "fidelity_local": float(source_row["fidelity"]),
            "V_error_local_K": float(source_row["V_error_K"]),
            "depth": int(compiled.depth()),
            "operation_count": int(compiled.size()),
            "operation_counts": {str(name): int(count)
                                 for name, count in compiled.count_ops().items()},
            "two_qubit_gate_count": two_qubit_count,
            "scheduled_duration_seconds": duration_seconds,
            "duration_over_median_T1": duration_seconds / median_t1,
            "duration_over_median_T2": duration_seconds / median_t2,
            **audit,
        })
    return rows


def main():
    started = time.perf_counter()
    accounts = QiskitRuntimeService.saved_accounts()
    if ACCOUNT_NAME not in accounts:
        raise RuntimeError(f"Saved account {ACCOUNT_NAME!r} was not found")
    account = accounts[ACCOUNT_NAME]
    preferences = list(account.get("plans_preference") or [])
    if preferences != ["open"]:
        raise RuntimeError(
            f"Refusing preflight: account is not Open-only: {preferences!r}")
    service = QiskitRuntimeService(name=ACCOUNT_NAME)
    usage = safe_usage(service.usage())
    source = json.loads((ROOT / "results" /
                         "adaptive_compressed_stateprep_N6.json").read_text(
                             encoding="utf-8"))
    graph = open_triangular_directional_edges()
    pauli0 = spin_hamiltonian_to_sparse_pauli(
        graph["N"], open_weak_model_bonds(graph, "M0"), 0.0)
    pauli2 = spin_hamiltonian_to_sparse_pauli(
        graph["N"], open_weak_model_bonds(graph, "M2"), 0.0)
    observable = (pauli2 - pauli0).simplify(atol=1e-13)

    backend_rows = []
    for backend in service.backends(operational=True, simulator=False):
        status = backend.status()
        compiled = compile_backend(backend, source)
        backend_rows.append({
            "name": backend.name,
            "num_qubits": int(backend.num_qubits),
            "pending_jobs": int(status.pending_jobs),
            "status_message": str(status.status_msg),
            "circuits": compiled,
            "summary": {
                "maximum_two_qubit_gate_count": max(
                    row["two_qubit_gate_count"] for row in compiled),
                "total_two_qubit_gate_count": sum(
                    row["two_qubit_gate_count"] for row in compiled),
                "minimum_calibration_product_all_gates": min(
                    row["calibration_product_all_gates"] for row in compiled),
                "maximum_duration_over_median_T1": max(
                    row["duration_over_median_T1"] for row in compiled),
                "maximum_duration_over_median_T2": max(
                    row["duration_over_median_T2"] for row in compiled),
            },
        })
    if not backend_rows:
        raise RuntimeError("No operational Open Plan QPU is visible")
    backend_rows.sort(key=lambda row: (
        -row["summary"]["minimum_calibration_product_all_gates"],
        row["summary"]["total_two_qubit_gate_count"],
        row["pending_jobs"], row["name"]))
    chosen = backend_rows[0]
    basis_count = len(observable.group_commuting(qubit_wise=True))
    remaining = usage.get("usage_remaining_seconds")
    gates_ok = (
        chosen["summary"]["maximum_two_qubit_gate_count"] <= 100 and
        chosen["summary"]["minimum_calibration_product_all_gates"] >= 0.65 and
        chosen["summary"]["maximum_duration_over_median_T1"] <= 0.10)
    time_ok = remaining is not None and float(remaining) >= 180
    decision = ("GO_COMPLETE_TWIRLED_TREX_PILOT"
                if gates_ok and time_ok else
                "HOLD_NO_LIVE_SUBMISSION")
    payload = {
        "purpose": "read-only live-backend preflight for complete adaptive N=6 spectral flow",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "account": {"saved_name": ACCOUNT_NAME,
                    "plans_preference": preferences,
                    "credential_or_instance_identifier_recorded": False},
        "open_plan_usage": usage,
        "source": "results/adaptive_compressed_stateprep_N6.json",
        "observable": {"name": "V=M2-M0",
                       "pauli_term_count": int(len(observable)),
                       "qwc_group_count": basis_count},
        "backends": backend_rows,
        "chosen_backend": chosen["name"],
        "chosen_summary": chosen["summary"],
        "decision": decision,
        "decision_basis": {"gates_ok": gates_ok, "time_ok": time_ok,
                           "minimum_remaining_free_seconds": 180,
                           "maximum_two_qubit_gates_per_circuit": 100,
                           "minimum_calibration_product": 0.65,
                           "maximum_duration_over_median_T1": 0.10},
        "proposed_job": {
            "pub_count": 9,
            "resilience_level": 1,
            "gate_twirling": True,
            "dynamical_decoupling": "XpXm",
            "target_precision_K": 0.02,
            "job_submitted_by_this_script": False,
            "scientific_claim": "complete N=6 entropy-correction reconstruction with hardware uncertainty; no bulk-material claim"},
        "scope_limits": [
            "Calibration products are heuristics derived from the current backend snapshot, not circuit fidelities.",
            "The N=6 open patch is a method benchmark and does not establish NBCP bulk thermodynamics.",
            "No Runtime job is submitted by this script.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "adaptive_n6_live_preflight.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "usage": usage,
        "chosen_backend": chosen["name"],
        "chosen_summary": chosen["summary"],
        "decision": decision, "elapsed_seconds": payload["elapsed_seconds"]
    }, indent=2))


if __name__ == "__main__":
    main()
