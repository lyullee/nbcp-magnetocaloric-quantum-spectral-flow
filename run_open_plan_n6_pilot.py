"""Validate or submit the minimal N=6 IBM Open Plan QPU pilot.

Default execution is a dry run.  Passing ``--submit`` is the only path that
creates a workload.  The script refuses to continue unless the saved account
is restricted to the Open Plan and free QPU allowance remains.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime.executor_estimator import Estimator

from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds,
                                  spin_hamiltonian_to_sparse_pauli)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ
from spectral_flow_circuit_unit_n6 import FIELD_T, state_preparation_circuit


ROOT = Path(__file__).resolve().parent
ACCOUNT_NAME = "psa_open"
LAMBDA = 0.5
LEVEL = 0
TARGET_PRECISION = 0.01
SEED = 20261012


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def sanitized_usage(raw):
    keep = ("usage_consumed_seconds", "usage_consumed_ru", "usage_limit_seconds",
            "usage_limit_ru", "usage_limit_reached", "usage_remaining_seconds",
            "usage_period")
    return {key: raw.get(key) for key in keep if key in raw}


def scalar(value):
    array = np.asarray(value)
    return float(array.reshape(-1)[0])


def make_problem():
    graph = open_triangular_directional_edges()
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    bonds0 = open_weak_model_bonds(graph, "M0")
    bonds2 = open_weak_model_bonds(graph, "M2")
    h0 = dense_hamiltonian(n, bonds0, field_k)
    h2 = dense_hamiltonian(n, bonds2, field_k)
    v = h2 - h0
    _, vectors = np.linalg.eigh(h0 + LAMBDA * v)
    circuit, target = state_preparation_circuit(vectors[:, LEVEL])
    pauli0 = spin_hamiltonian_to_sparse_pauli(n, bonds0, 0.0)
    pauli2 = spin_hamiltonian_to_sparse_pauli(n, bonds2, 0.0)
    observable = (pauli2 - pauli0).simplify(atol=1e-13)
    exact = float(np.real(Statevector(target).expectation_value(observable)))
    return circuit, observable, exact


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true",
                        help="submit one mitigated Estimator job")
    parser.add_argument("--resilience-level", type=int, choices=(1, 2),
                        default=1,
                        help="1=TREX; 2=TREX plus gate twirling and ZNE")
    args = parser.parse_args()
    started = time.perf_counter()

    accounts = QiskitRuntimeService.saved_accounts()
    if ACCOUNT_NAME not in accounts:
        raise RuntimeError(f"Saved account {ACCOUNT_NAME!r} was not found")
    account = accounts[ACCOUNT_NAME]
    preferences = list(account.get("plans_preference") or [])
    if preferences != ["open"]:
        raise RuntimeError("Submission is allowed only through an Open-only account")

    service = QiskitRuntimeService(name=ACCOUNT_NAME)
    before = sanitized_usage(service.usage())
    if before.get("usage_limit_reached"):
        raise RuntimeError("Open Plan QPU allowance is exhausted")
    remaining = before.get("usage_remaining_seconds")
    if remaining is not None and float(remaining) < 30:
        raise RuntimeError(f"Only {remaining} free QPU seconds remain; refusing pilot")

    preflight = json.loads((ROOT / "results" /
                            "open_plan_preflight_N6.json").read_text(
                                encoding="utf-8"))
    if preflight["decision"] != "GO_MINIMAL_TREX_PILOT":
        raise RuntimeError("Live-backend preflight did not pass")
    backend = service.backend(preflight["chosen_backend"])
    if not backend.status().operational:
        raise RuntimeError(f"Backend {backend.name} is not operational")

    circuit, observable, exact = make_problem()
    isa_circuit = transpile(
        circuit, backend=backend, optimization_level=3,
        seed_transpiler=SEED, scheduling_method="alap")
    isa_observable = observable.apply_layout(isa_circuit.layout)
    resilience_level = int(args.resilience_level)
    mitigation = ("TREX readout mitigation" if resilience_level == 1 else
                  "TREX + gate twirling + zero-noise extrapolation")
    options = {
        "resilience_level": resilience_level,
        "default_precision": TARGET_PRECISION,
        "environment": {
            "job_tags": ["nbcp", "spectral-flow", "n6-pilot", "open-plan",
                         f"resilience-{resilience_level}"]
        },
    }
    submitted_utc = datetime.now(timezone.utc).isoformat()
    base = {
        "purpose": "minimal live-QPU execution-path proof for N=6 spectral flow",
        "created_utc": submitted_utc,
        "submission_requested": bool(args.submit),
        "open_plan_only_guard_passed": True,
        "account_name": ACCOUNT_NAME,
        "backend": backend.name,
        "usage_before": before,
        "condition": {"N": 6, "field_T": FIELD_T, "lambda": LAMBDA,
                      "level": LEVEL, "path": "M0 + lambda*(M2-M0)"},
        "observable": "V=M2-M0",
        "exact_V_expectation_K": exact,
        "target_precision_K": TARGET_PRECISION,
        "resilience_level": resilience_level,
        "mitigation": mitigation,
        "compiled": {
            "depth": int(isa_circuit.depth()),
            "size": int(isa_circuit.size()),
            "operation_counts": {
                str(name): int(count)
                for name, count in isa_circuit.count_ops().items()
            },
        },
        "scientific_scope": (
            "QPU execution-path demonstration only; N=6 does not support the "
            "NBCP material or engineering conclusion"),
    }

    suffix = "TREX" if resilience_level == 1 else "ZNE"
    filename = (f"open_plan_qpu_pilot_N6_{suffix}.json" if args.submit else
                f"open_plan_qpu_pilot_N6_{suffix}_dry_run.json")
    output = ROOT / "results" / filename
    if not args.submit:
        base.update({
            "status": "LOCAL_SUBMISSION_PAYLOAD_VALIDATED",
            "job_submitted": False,
            "validation_note": (
                "Circuit transpilation, observable layout, Open-only account "
                "guard and Estimator options were validated locally. No "
                "Runtime job was created."),
            "elapsed_seconds": time.perf_counter() - started,
        })
        output.write_text(json.dumps(base, indent=2), encoding="utf-8")
        print(json.dumps(base, indent=2))
        return

    estimator = Estimator(mode=backend, options=options)
    job = estimator.run([(isa_circuit, isa_observable)],
                        precision=TARGET_PRECISION)
    job_id = job.job_id()
    base.update({"job_submitted": True, "job_id": job_id,
                 "status_after_submission": str(job.status())})
    output.write_text(json.dumps(base, indent=2), encoding="utf-8")
    print(json.dumps({"job_id": job_id, "backend": backend.name,
                      "status": str(job.status())}, indent=2), flush=True)

    result = job.result()
    pub = result[0]
    measured = scalar(pub.data.evs)
    standard_error = scalar(pub.data.stds)
    metadata = dict(pub.metadata)
    # Metadata is expected to contain only numerical execution information.
    base.update({
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "status": "DONE",
        "measured_V_expectation_K": measured,
        "reported_standard_error_K": standard_error,
        "signed_error_K": measured - exact,
        "absolute_error_K": abs(measured - exact),
        "within_two_reported_standard_errors": (
            abs(measured - exact) <= 2 * standard_error),
        "result_metadata": metadata,
        "usage_after": sanitized_usage(service.usage()),
        "elapsed_seconds": time.perf_counter() - started,
    })
    output.write_text(json.dumps(base, indent=2, default=str), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "job_id": job_id,
        "measured_V_expectation_K": measured,
        "exact_V_expectation_K": exact,
        "reported_standard_error_K": standard_error,
        "absolute_error_K": abs(measured - exact),
        "usage_before": base["usage_before"],
        "usage_after": base["usage_after"],
        "elapsed_seconds": base["elapsed_seconds"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
