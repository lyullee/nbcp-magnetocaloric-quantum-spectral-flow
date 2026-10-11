"""Read-only live-backend preflight for the N=6 spectral-flow QPU pilot.

The script never submits a job.  It verifies that the saved credential selects
only an IBM Open Plan instance, records the remaining free allowance, compiles
the nine N=6 state-preparation circuits against the least-busy live backend,
and writes a sanitized JSON record without credentials or the service CRN.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import time
import warnings

import numpy as np
import qiskit
from qiskit import transpile
import qiskit_ibm_runtime
from qiskit_ibm_runtime import QiskitRuntimeService

from audit_fake_backend_noise_n9 import circuit_error_product
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ
from spectral_flow_circuit_unit_n6 import (FIELD_T, GL_NODES,
                                           state_preparation_circuit)


ROOT = Path(__file__).resolve().parent
ACCOUNT_NAME = "psa_open"
SEED = 20261012


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def safe_usage(raw):
    """Keep usage facts while omitting account/instance identifiers."""
    keep = ("usage_consumed_seconds", "usage_consumed_ru", "usage_limit_seconds",
            "usage_limit_ru", "usage_limit_reached", "usage_remaining_seconds",
            "usage_period")
    return {key: raw.get(key) for key in keep if key in raw}


def backend_row(backend):
    status = backend.status()
    return {
        "name": backend.name,
        "num_qubits": int(backend.num_qubits),
        "operational": bool(status.operational),
        "pending_jobs": int(status.pending_jobs),
        "status_message": str(status.status_msg),
    }


def main():
    started = time.perf_counter()
    accounts = QiskitRuntimeService.saved_accounts()
    if ACCOUNT_NAME not in accounts:
        raise RuntimeError(f"Saved account {ACCOUNT_NAME!r} was not found")
    account = accounts[ACCOUNT_NAME]
    preferences = list(account.get("plans_preference") or [])
    if preferences != ["open"]:
        raise RuntimeError(
            "Refusing preflight because the saved account is not restricted "
            f"to the Open Plan: {preferences!r}")

    service = QiskitRuntimeService(name=ACCOUNT_NAME)
    usage = safe_usage(service.usage())
    backends = service.backends(operational=True, simulator=False)
    if not backends:
        raise RuntimeError("No operational Open Plan QPU is visible")
    backend_rows = sorted((backend_row(item) for item in backends),
                          key=lambda row: (row["pending_jobs"], row["name"]))
    chosen_name = backend_rows[0]["name"]
    backend = service.backend(chosen_name)

    source = json.loads((ROOT / "results" /
                         "spectral_flow_circuit_unit_N6.json").read_text(
                             encoding="utf-8"))
    rank = int(source["low_energy_rank"])
    graph = open_triangular_directional_edges()
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    v = h2 - h0

    rows = []
    for node in GL_NODES:
        _, vectors = np.linalg.eigh(h0 + float(node) * v)
        for level in range(rank):
            circuit, _ = state_preparation_circuit(vectors[:, level])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                compiled = transpile(
                    circuit, backend=backend, optimization_level=3,
                    seed_transpiler=SEED, scheduling_method="alap")
            audit = circuit_error_product(compiled, backend)
            duration_dt = int(compiled.duration)
            duration_seconds = duration_dt * float(backend.dt)
            properties = [backend.qubit_properties(qubit)
                          for qubit in audit["used_physical_qubits"]]
            median_t1 = float(np.median([item.t1 for item in properties]))
            median_t2 = float(np.median([item.t2 for item in properties]))
            two_qubit_count = sum(
                count for gate, count in compiled.count_ops().items()
                if gate in {"cx", "cz", "ecr", "iswap", "rzz"})
            row = {
                "lambda": float(node),
                "level": int(level),
                "depth": int(compiled.depth()),
                "operation_count": int(compiled.size()),
                "two_qubit_gate_count": int(two_qubit_count),
                "operation_counts": {
                    str(name): int(count)
                    for name, count in compiled.count_ops().items()
                },
                "scheduled_duration_seconds": duration_seconds,
                "duration_over_median_T1": duration_seconds / median_t1,
                "duration_over_median_T2": duration_seconds / median_t2,
                **audit,
            }
            rows.append(row)
            print(f"{chosen_name} lambda={node:.6f} L={level}: "
                  f"2q={two_qubit_count}, "
                  f"cal-product={audit['calibration_product_all_gates']:.3f}, "
                  f"duration/T1={row['duration_over_median_T1']:.2f}",
                  flush=True)

    min_product = min(row["calibration_product_all_gates"] for row in rows)
    max_t1_ratio = max(row["duration_over_median_T1"] for row in rows)
    decision = ("GO_MINIMAL_TREX_PILOT"
                if min_product >= 0.20 and max_t1_ratio <= 0.50
                else "HOLD_COMPRESS_BEFORE_LIVE_RUN")
    payload = {
        "purpose": "read-only live Open Plan preflight for the N=6 QPU pilot",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "qiskit": qiskit.__version__,
            "qiskit_ibm_runtime": qiskit_ibm_runtime.__version__,
        },
        "account": {
            "saved_name": ACCOUNT_NAME,
            "channel": account.get("channel"),
            "region": account.get("region"),
            "plans_preference": preferences,
            "credential_or_instance_identifier_recorded": False,
        },
        "open_plan_usage": usage,
        "visible_operational_backends": backend_rows,
        "chosen_backend": chosen_name,
        "source": "results/spectral_flow_circuit_unit_N6.json",
        "compiled_circuits": rows,
        "decision": decision,
        "decision_basis": {
            "minimum_calibration_product_all_gates": min_product,
            "maximum_duration_over_median_T1": max_t1_ratio,
            "maximum_duration_over_median_T2": max(
                row["duration_over_median_T2"] for row in rows),
            "maximum_two_qubit_gate_count": max(
                row["two_qubit_gate_count"] for row in rows),
        },
        "proposed_minimal_pilot": {
            "scope": "central Gauss-Legendre node, ground state, V expectation",
            "lambda": 0.5,
            "level": 0,
            "resilience_level": 1,
            "mitigation": "TREX readout mitigation; no ZNE in first run",
            "target_precision": 0.05,
            "job_submitted_by_this_script": False,
            "scientific_claim": "execution-path demonstration only; no material conclusion",
        },
        "scope_limits": [
            "Calibration products use the backend snapshot visible during preflight and are not circuit fidelities.",
            "The selected backend and calibration can change before submission.",
            "The returned usage period is recorded verbatim and must be refreshed immediately before a job.",
            "The script does not submit a job and cannot consume QPU time.",
            "N=6 is a circuit proof only; material conclusions require the N=9/N=12 simulator evidence.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "open_plan_preflight_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "chosen_backend": chosen_name,
        "usage": usage,
        "decision": decision,
        **payload["decision_basis"],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
