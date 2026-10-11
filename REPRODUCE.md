# Reproduction guide

## 1. Environment

The frozen environment used Python 3.14.8 on Windows 11 with the versions pinned in `requirements.txt`. Octave 11.3.0 was used for independent checks; its local runtime is excluded from version control.

## 2. Obtain public inputs

Download the three Nature source-data workbooks and verify their hashes:

```text
python fetch_sources.py
```

Download and verify the public arXiv and Oxford archives:

```text
python fetch_public_inputs.py --download --extract
```

Third-party QMagen checkouts are optional because their frozen derived results are included. To obtain them at the recorded commits:

```text
python fetch_public_inputs.py --with-qmagen
```

Official URLs, byte counts, commits, and SHA-256 hashes are recorded in `source_manifest.json` and `sources/input_manifest_20261011.json`.

## 3. Fast audit

The fast audit uses frozen results and normally completes in seconds:

```text
python verify_engineering_decision.py
python verify_quantum_stateprep_qpu.py
python verify_second_screen.py
python test_quantum_local_benchmark.py
python test_bulk_polarized_stability.py
python test_interlayer_robustness.py
python verify_manuscript_claims.py
python make_bulk_manifest.py
```

Expected status is `PASS` for every verifier. The final manifest records hashes and runtime versions.

## 4. Principal classical calculations

Run in this order when regenerating the principal evidence:

```text
python refine_engineering_case.py
python bulk_polarized_stability.py
python validate_bulk_stability_convergence.py
python bulk_to_engineering_surrogate.py
python finite_cluster_stability_alignment.py
python experimental_hc_engineering_anchor.py
python interlayer_engineering_scenarios.py
python screen_interlayer_parameter_ensemble.py
python stress_test_thermodynamic_shape.py
python design_peak_measurement_target.py
python compile_refrigerator_design_chart.py
python verify_engineering_decision.py
```

N=15/N=18 FTLM and the larger low-magnon paths are more expensive. Their individual settings, seeds, and outputs are frozen in the corresponding scripts and JSON files. Start with `benchmark_ftlm.py`, then follow task order T004–T010K in `연구실행_체크리스트_20261011.csv`.

## 5. Quantum simulation

The non-hardware path is:

```text
python prepare_quantum_local_benchmark.py
python screen_spectral_flow_entropy_n9.py
python screen_spectral_flow_entropy_n12.py
python spectral_flow_circuit_unit_n6.py
python fit_spectral_flow_states_n9.py --field 1.40
python fit_spectral_flow_states_n9.py --field 1.45
python fit_spectral_flow_states_n9.py --field 1.50
python compile_n9_stateprep_field_scan.py
python verify_quantum_stateprep_qpu.py
```

The field values accepted by `fit_spectral_flow_states_n9.py` are 1.40, 1.45, and 1.50 T. The compiled scan must pass all nine field-temperature cases.

## 6. IBM hardware safeguard

No hardware job is required to reproduce the paper conclusions. Frozen sanitized records contain both completed job IDs. The submission script:

- refuses any account whose plan preference is not exactly `open`;
- checks remaining free usage;
- defaults to dry run;
- requires explicit `--submit`;
- records no token or cloud resource name.

Do not resubmit merely to duplicate the archived pilot. The recorded stopping rule requires prior circuit compression and observable-level noisy evidence.

## 7. Interpretation boundaries

- P2 is a public-propagation-vector-constrained set of seven models retained from a deterministic 375-member local grid; it is not a posterior distribution.
- The critical-field transport is a decision-scale surrogate, not a finite-temperature three-dimensional M4 solution.
- The finite-temperature deformation coordinates and their pass fractions are deterministic stress coverage, not probabilities or confidence intervals.
- The 4 mT measurement sweep controls peak-grid localization only; calibration and thermodynamic measurement uncertainty remain experimental requirements.
- The 90–300 mK `W95` result is a set of measured field candidates.
- Conditional volumes include magnetic entropy only.
- N=9 state-preparation targets use exact classical eigenvectors.
- The N=6 QPU pilot does not establish quantitative accuracy or quantum advantage.
