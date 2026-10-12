"""Submit or dry-run the complete adaptive N=6 Open Plan QPU pilot.

Nine frozen circuits estimate all three low levels at all three spectral-flow
quadrature nodes.  The default is a dry run.  ``--submit`` is guarded by an
Open-only saved account, a passing live preflight, and remaining free time.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np
from qiskit import transpile
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime.executor_estimator import Estimator

from fit_number_conserving_ansatz_n9 import (
    build_qiskit_circuit,
    dense_hamiltonian,
)
from qiskit_pauli_mapping import (
    open_triangular_directional_edges,
    open_weak_model_bonds,
    spin_hamiltonian_to_sparse_pauli,
)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ
from spectral_flow_circuit_unit_n6 import entropy


ROOT = Path(__file__).resolve().parent
ACCOUNT_NAME = "psa_open"
TARGET_PRECISION_K = 0.02
SEED = 20261013
MINIMUM_FREE_SECONDS = 180


def safe_usage(raw):
    keep = ("usage_consumed_seconds", "usage_consumed_ru",
            "usage_limit_seconds", "usage_limit_ru", "usage_limit_reached",
            "usage_remaining_seconds", "usage_period")
    return {key: raw.get(key) for key in keep if key in raw}


def scalar(value):
    return float(np.asarray(value).reshape(-1)[0])


def build_problems(source, backend, observable):
    problems = []
    for node in source["nodes"]:
        for row in node["states"]:
            operations = [tuple(item) for item in row["operations"]]
            parameters = np.asarray(row["parameters_radians"], dtype=float)
            circuit = build_qiskit_circuit(
                source["condition"]["N"], int(row["initial_basis_index"]),
                operations, parameters)
            circuit.metadata = {"lambda": float(row["lambda"]),
                                "level": int(row["level"])}
            isa_circuit = transpile(
                circuit, backend=backend, optimization_level=3,
                seed_transpiler=SEED, scheduling_method="alap")
            isa_observable = observable.apply_layout(isa_circuit.layout)
            two_qubit_count = sum(
                int(count) for gate, count in isa_circuit.count_ops().items()
                if gate in {"cx", "cz", "ecr", "iswap", "rzz"})
            problems.append({
                "lambda": float(row["lambda"]),
                "weight": float(row["weight"]),
                "level": int(row["level"]),
                "exact_V_expectation_K": float(row["exact_V_expectation_K"]),
                "prepared_V_expectation_K": float(
                    row["prepared_V_expectation_K"]),
                "local_stateprep_V_error_K": float(row["V_error_K"]),
                "compiled_depth": int(isa_circuit.depth()),
                "compiled_size": int(isa_circuit.size()),
                "two_qubit_gate_count": two_qubit_count,
                "operation_counts": {str(name): int(count) for name, count
                                     in isa_circuit.count_ops().items()},
                "circuit": isa_circuit,
                "observable": isa_observable,
            })
    return problems


def thermodynamic_summary(problems, measured, standard_errors):
    graph = open_triangular_directional_edges()
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * float(
        json.loads((ROOT / "results" /
                    "adaptive_compressed_stateprep_N6.json").read_text(
                        encoding="utf-8"))["condition"]["field_T"])
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    e0 = np.linalg.eigvalsh(h0)
    e2 = np.linalg.eigvalsh(h2)
    rank = 3
    shifts = np.zeros(rank)
    exact_shifts = np.zeros(rank)
    prepared_shifts = np.zeros(rank)
    for problem, value in zip(problems, measured, strict=True):
        level = problem["level"]
        weight = problem["weight"]
        shifts[level] += weight * value
        exact_shifts[level] += weight * problem["exact_V_expectation_K"]
        prepared_shifts[level] += weight * problem["prepared_V_expectation_K"]
    exact_full = (entropy(e2) - entropy(e0)) / n
    exact_low_rank = (entropy(e0[:rank] + exact_shifts) -
                      entropy(e0[:rank])) / n
    prepared_low_rank = (entropy(e0[:rank] + prepared_shifts) -
                         entropy(e0[:rank])) / n
    hardware_low_rank = (entropy(e0[:rank] + shifts) -
                         entropy(e0[:rank])) / n

    rng = np.random.default_rng(SEED)
    samples = np.empty(100000)
    measured_array = np.asarray(measured)
    error_array = np.asarray(standard_errors)
    for index in range(samples.size):
        draw = rng.normal(measured_array, error_array)
        draw_shifts = np.zeros(rank)
        for problem, value in zip(problems, draw, strict=True):
            draw_shifts[problem["level"]] += problem["weight"] * value
        samples[index] = (entropy(e0[:rank] + draw_shifts) -
                          entropy(e0[:rank])) / n
    lower, upper = np.quantile(samples, (0.025, 0.975))
    return {
        "exact_full_delta_S_kB_per_spin": float(exact_full),
        "exact_low_rank_flow_delta_S_kB_per_spin": float(exact_low_rank),
        "prepared_low_rank_delta_S_kB_per_spin": float(prepared_low_rank),
        "hardware_low_rank_delta_S_kB_per_spin": float(hardware_low_rank),
        "hardware_minus_exact_full_kB_per_spin": float(
            hardware_low_rank - exact_full),
        "hardware_sign_matches_exact": bool(
            np.sign(hardware_low_rank) == np.sign(exact_full)),
        "reported_uncertainty_propagation": {
            "method": "100000 independent normal draws from Estimator-reported standard errors",
            "seed": SEED,
            "mean_kB_per_spin": float(np.mean(samples)),
            "standard_deviation_kB_per_spin": float(np.std(samples, ddof=1)),
            "95_percent_interval_kB_per_spin": [float(lower), float(upper)],
            "covers_exact_full": bool(lower <= exact_full <= upper),
            "does_not_cover_zero": bool(not (lower <= 0 <= upper)),
            "limitation": "propagates reported statistical errors only; residual hardware bias is assessed against the known N=6 exact benchmark",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--resilience-level", type=int, choices=(1, 2),
                        default=1,
                        help="1=twirled TREX; 2=twirled TREX plus ZNE")
    args = parser.parse_args()
    started = time.perf_counter()
    accounts = QiskitRuntimeService.saved_accounts()
    if ACCOUNT_NAME not in accounts:
        raise RuntimeError(f"Saved account {ACCOUNT_NAME!r} was not found")
    preferences = list(accounts[ACCOUNT_NAME].get("plans_preference") or [])
    if preferences != ["open"]:
        raise RuntimeError(f"Refusing non-Open account: {preferences!r}")
    service = QiskitRuntimeService(name=ACCOUNT_NAME)
    before = safe_usage(service.usage())
    if before.get("usage_limit_reached"):
        raise RuntimeError("Open Plan usage limit is already reached")
    if float(before.get("usage_remaining_seconds") or 0) < MINIMUM_FREE_SECONDS:
        raise RuntimeError("Insufficient free QPU allowance for guarded job")
    preflight = json.loads((ROOT / "results" /
                            "adaptive_n6_live_preflight.json").read_text(
                                encoding="utf-8"))
    if preflight["decision"] != "GO_COMPLETE_TWIRLED_TREX_PILOT":
        raise RuntimeError("Frozen live preflight did not pass")
    backend = service.backend(preflight["chosen_backend"])
    if not backend.status().operational:
        raise RuntimeError(f"Backend {backend.name} is not operational")
    source = json.loads((ROOT / "results" /
                         "adaptive_compressed_stateprep_N6.json").read_text(
                             encoding="utf-8"))
    graph = open_triangular_directional_edges()
    pauli0 = spin_hamiltonian_to_sparse_pauli(
        graph["N"], open_weak_model_bonds(graph, "M0"), 0.0)
    pauli2 = spin_hamiltonian_to_sparse_pauli(
        graph["N"], open_weak_model_bonds(graph, "M2"), 0.0)
    observable = (pauli2 - pauli0).simplify(atol=1e-13)
    problems = build_problems(source, backend, observable)
    max_2q = max(problem["two_qubit_gate_count"] for problem in problems)
    if max_2q > 110:
        raise RuntimeError(f"Live recompile exceeded gate guard: {max_2q}")

    created = datetime.now(timezone.utc).isoformat()
    resilience_level = int(args.resilience_level)
    mitigation_label = "TWIRLED_TREX" if resilience_level == 1 else "TWIRLED_TREX_ZNE"
    base = {
        "purpose": "complete live-QPU N=6 spectral-flow entropy reconstruction",
        "created_utc": created,
        "submission_requested": bool(args.submit),
        "open_plan_only_guard_passed": True,
        "account_name": ACCOUNT_NAME,
        "backend": backend.name,
        "usage_before": before,
        "condition": source["condition"],
        "observable": {"name": "V=M2-M0", "pauli_term_count": int(len(observable)),
                       "qwc_group_count": len(observable.group_commuting(qubit_wise=True))},
        "target_precision_K": TARGET_PRECISION_K,
        "mitigation": {
            "resilience_level": resilience_level,
            "TREX_measurement_mitigation": True,
            "gate_twirling": True,
            "dynamical_decoupling": "XpXm",
            "ZNE": resilience_level == 2,
        },
        "circuits": [{key: value for key, value in problem.items()
                      if key not in {"circuit", "observable"}}
                     for problem in problems],
        "scope": "complete N=6 method benchmark; no N=18 or bulk-material inference",
    }
    output = ROOT / "results" / (
        f"adaptive_n6_qpu_complete_{mitigation_label}.json" if args.submit else
        f"adaptive_n6_qpu_complete_{mitigation_label}_dry_run.json")
    if not args.submit:
        base.update({"status": "LOCAL_PAYLOAD_VALIDATED", "job_submitted": False,
                     "elapsed_seconds": time.perf_counter() - started})
        output.write_text(json.dumps(base, indent=2), encoding="utf-8")
        print(json.dumps({"output": str(output), "status": base["status"],
                          "backend": backend.name, "pub_count": len(problems),
                          "maximum_two_qubit_gate_count": max_2q,
                          "usage_before": before}, indent=2))
        return

    options = {
        "resilience_level": resilience_level,
        "default_precision": TARGET_PRECISION_K,
        "dynamical_decoupling": {"enable": True,
                                 "sequence_type": "XpXm"},
        "twirling": {"enable_gates": True, "enable_measure": True},
        "environment": {"job_tags": ["nbcp", "spectral-flow", "n6-complete",
                                      "open-plan", mitigation_label.lower()]},
    }
    estimator = Estimator(mode=backend, options=options)
    pubs = [(problem["circuit"], problem["observable"])
            for problem in problems]
    job = estimator.run(pubs, precision=TARGET_PRECISION_K)
    base.update({"job_submitted": True, "job_id": job.job_id(),
                 "status_after_submission": str(job.status())})
    output.write_text(json.dumps(base, indent=2), encoding="utf-8")
    print(json.dumps({"job_id": job.job_id(), "backend": backend.name,
                      "pub_count": len(pubs), "status": str(job.status())},
                     indent=2), flush=True)
    result = job.result()
    measured = [scalar(pub.data.evs) for pub in result]
    standard_errors = [scalar(pub.data.stds) for pub in result]
    observations = []
    for problem, value, error, pub in zip(
            problems, measured, standard_errors, result, strict=True):
        observations.append({
            "lambda": problem["lambda"], "weight": problem["weight"],
            "level": problem["level"],
            "exact_V_expectation_K": problem["exact_V_expectation_K"],
            "prepared_V_expectation_K": problem["prepared_V_expectation_K"],
            "hardware_V_expectation_K": value,
            "reported_standard_error_K": error,
            "hardware_minus_exact_K": value - problem["exact_V_expectation_K"],
            "hardware_minus_prepared_K": value - problem["prepared_V_expectation_K"],
            "metadata": dict(pub.metadata),
        })
    base.update({
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "status": "DONE",
        "observations": observations,
        "thermodynamics": thermodynamic_summary(
            problems, measured, standard_errors),
        "usage_after": safe_usage(service.usage()),
        "elapsed_seconds": time.perf_counter() - started,
    })
    output.write_text(json.dumps(base, indent=2, default=str), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "job_id": base["job_id"],
        "thermodynamics": base["thermodynamics"],
        "usage_before": base["usage_before"],
        "usage_after": base["usage_after"],
        "elapsed_seconds": base["elapsed_seconds"],
    }, indent=2, default=str))


if __name__ == "__main__":
    main()
