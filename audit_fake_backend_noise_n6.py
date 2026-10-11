"""FakeSherbrooke calibration-product audit for the N=6 circuit unit test."""
from __future__ import annotations

import json
from pathlib import Path
import time
import warnings

import numpy as np
from qiskit import transpile
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke

from audit_fake_backend_noise_n9 import circuit_error_product
from matrix_free_hamiltonian import SpinHamiltonianOperator
from qiskit_pauli_mapping import (open_triangular_directional_edges,
                                  open_weak_model_bonds)
from screen_model import MU_B_OVER_KB
from screen_weak_terms import GZ
from spectral_flow_circuit_unit_n6 import (FIELD_T, GL_NODES,
                                           state_preparation_circuit)


ROOT = Path(__file__).resolve().parent
SEED = 20261012


def dense_hamiltonian(n, bonds, field_k):
    operator = SpinHamiltonianOperator(n, bonds, field_k)
    return np.asarray(operator @ np.eye(1 << n, dtype=complex))


def main():
    started = time.perf_counter()
    source = json.loads((ROOT / "results" /
                         "spectral_flow_circuit_unit_N6.json").read_text(
                             encoding="utf-8"))
    rank = source["low_energy_rank"]
    graph = open_triangular_directional_edges()
    n = graph["N"]
    field_k = GZ * MU_B_OVER_KB * FIELD_T
    h0 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M0"), field_k)
    h2 = dense_hamiltonian(n, open_weak_model_bonds(graph, "M2"), field_k)
    backend = FakeSherbrooke()
    rows = []
    for node in GL_NODES:
        _, vectors = np.linalg.eigh(h0 + float(node) * (h2 - h0))
        for level in range(rank):
            circuit, _ = state_preparation_circuit(vectors[:, level])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                compiled = transpile(
                    circuit, backend=backend, optimization_level=1,
                    seed_transpiler=SEED, scheduling_method="alap")
                duration_dt = int(compiled.duration)
            audit = circuit_error_product(compiled, backend)
            properties = [backend.qubit_properties(qubit)
                          for qubit in audit["used_physical_qubits"]]
            duration_seconds = duration_dt * backend.dt
            median_t1 = float(np.median([item.t1 for item in properties]))
            median_t2 = float(np.median([item.t2 for item in properties]))
            rows.append({
                "lambda": float(node), "level": level,
                "compiled_ecr_count": int(
                    compiled.count_ops().get("ecr", 0)),
                "compiled_depth": int(compiled.depth()),
                "scheduled_duration_seconds": duration_seconds,
                "duration_over_median_T1": duration_seconds / median_t1,
                "duration_over_median_T2": duration_seconds / median_t2,
                **audit,
            })
            print(f"lambda={node:.6f} L={level}: "
                  f"survival={audit['calibration_product_all_gates']:.3f}, "
                  f"duration/T1={duration_seconds/median_t1:.2f}", flush=True)

    payload = {
        "purpose": "noise pre-screen for the N=6 free-QPU circuit pilot",
        "backend_model": backend.name,
        "source": "results/spectral_flow_circuit_unit_N6.json",
        "rows": rows,
        "decision": "CONDITIONAL_MITIGATED_PILOT_CANDIDATE",
        "decision_basis": {
            "minimum_calibration_product_all_gates": min(
                row["calibration_product_all_gates"] for row in rows),
            "maximum_duration_over_median_T1": max(
                row["duration_over_median_T1"] for row in rows),
            "maximum_duration_over_median_T2": max(
                row["duration_over_median_T2"] for row in rows),
        },
        "scope_limits": [
            "The calibration product is a heuristic rather than circuit fidelity.",
            "FakeSherbrooke is historical and must be replaced with a live-backend check immediately before any run.",
            "Readout, measurement grouping, zero-noise extrapolation and optimizer costs are excluded.",
            "N=6 is a circuit proof only and cannot support the material conclusion by itself.",
        ],
        "elapsed_seconds": time.perf_counter() - started,
    }
    output = ROOT / "results" / "fake_backend_noise_audit_N6.json"
    output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output), "decision": payload["decision"],
        **payload["decision_basis"],
        "elapsed_seconds": payload["elapsed_seconds"],
    }, indent=2))


if __name__ == "__main__":
    main()
