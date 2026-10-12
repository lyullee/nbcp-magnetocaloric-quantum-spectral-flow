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

The independent YbNi1.6Sn device audit downloads its 277,226-byte public archive directly from Cambridge DOI 10.17863/CAM.106770 when absent, verifies SHA-256 `DC42ADBDEA8159967B0A09F35010F90F04C0A9EE4EEFAD41AE47D61D0342997C`, and leaves the CC BY-NC-ND raw files under the ignored `sources/raw/` directory.

## 3. Fast audit

The fast audit uses frozen results and normally completes in seconds:

```text
python verify_engineering_decision.py
python verify_quantum_stateprep_qpu.py
python compile_adaptive_n6_qpu_evidence.py
python verify_adaptive_n6_qpu_evidence.py
python verify_second_screen.py
python test_quantum_local_benchmark.py
python test_bulk_polarized_stability.py
python test_interlayer_robustness.py
python audit_interlayer_envelope_boundaries.py
python screen_interlayer_expanded_sobol.py
python verify_local_interlayer_modes.py
python verify_expanded_sobol_modes.py
python stress_test_expanded_sobol_shape.py
python extract_refrigerant_benchmark.py
python compile_refrigerant_common_basis.py
python verify_refrigerant_common_basis.py
python compile_adr_stage_load_envelope.py
python verify_adr_stage_load_envelope.py
python compile_adr_transient_cycle_envelope.py
python verify_adr_transient_cycle_envelope.py
python audit_second_material_public_validation.py
python verify_second_material_validation.py
python audit_external_stage_heat_balance.py
python verify_external_stage_heat_balance.py
python verify_measurement_value_of_information.py
python test_general_truncated_magnon_operator.py
python verify_n18_multimodel_plan.py
python verify_n18_stacked_m4_convergence.py
python verify_n18_multimodel_cover.py
python verify_n18_multimodel_design.py
python integrate_n18_multimodel_manuscript.py
python verify_cluster_size_transfer.py
python verify_sample_demagnetization_envelope.py
python verify_manuscript_claims.py
python prepare_cryogenics_submission.py
python verify_cryogenics_submission.py
python build_cryogenics_docx.py
python verify_cryogenics_docx.py
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
python audit_interlayer_envelope_boundaries.py
python screen_interlayer_expanded_sobol.py
python verify_local_interlayer_modes.py
python verify_expanded_sobol_modes.py
python stress_test_expanded_sobol_shape.py
python plot_expanded_sobol_audit.py
python adaptive_interlayer_inverse_uq.py
python verify_adaptive_interlayer_inverse_uq.py
python audit_field_dependent_fan_holdout.py
python verify_field_dependent_fan_holdout.py
python audit_public_ins_interlayer_bandwidth.py
python verify_public_ins_interlayer_bandwidth.py
python audit_m4_finite_temperature_cluster.py
python verify_m4_finite_temperature_cluster.py
python test_general_truncated_magnon_operator.py
python compile_n18_stacked_m4_convergence.py
python verify_n18_stacked_m4_convergence.py
python verify_n18_multimodel_plan.py
python run_n18_multimodel_cover.py --stage all
python compile_n18_multimodel_cover.py
python verify_n18_multimodel_cover.py
python stress_test_n18_multimodel_design.py
python verify_n18_multimodel_design.py
python integrate_n18_multimodel_manuscript.py
python stress_test_thermodynamic_shape.py
python stress_test_microscopic_shape_templates.py
python verify_microscopic_shape_templates.py
python optimize_measurement_value_of_information.py
python verify_measurement_value_of_information.py
python stress_test_cluster_size_transfer.py
python verify_cluster_size_transfer.py
python audit_sample_demagnetization_envelope.py
python verify_sample_demagnetization_envelope.py
python extract_refrigerant_benchmark.py
python compile_refrigerant_common_basis.py
python verify_refrigerant_common_basis.py
python compile_refrigerator_design_chart.py
python compile_adr_stage_load_envelope.py
python verify_adr_stage_load_envelope.py
python compile_adr_transient_cycle_envelope.py
python verify_adr_transient_cycle_envelope.py
python audit_second_material_public_validation.py
python verify_second_material_validation.py
python audit_external_stage_heat_balance.py
python verify_external_stage_heat_balance.py
python plot_evidence_to_decision_workflow.py
python verify_engineering_decision.py
```

N=15/N=18 FTLM and the larger low-magnon paths are more expensive. Their individual settings, seeds, and outputs are frozen in the corresponding scripts and JSON files. The N=18 model cover fixes M3 plus minimum-shift, nominal, median-shift/distinct-coupling, and maximum-shift M4 anchors before the new calculations. Its low-accuracy pilot revealed separated competing peaks, so `results/n18_multimodel_plan_amendment_after_coarse.json` expands the high-accuracy calculation to the common 1.10–2.10 T grid at 0.02 T spacing before any new K6/K7 cover result is inspected. Before that dense scan completed, `results/n18_multimodel_convergence_protocol.json` fixed a close-peak L32 rule. The original 0.01 kB/spin candidate gap was conservatively widened to 0.02 in `results/n18_multimodel_convergence_protocol_amendment.json` before any dense worker completed, because existing M3 L12-to-L32 corrections reached 0.01023 kB/spin. Every included local maximum receives an L32 three-point calculation before the peak is selected. Start with `benchmark_ftlm.py`, then follow task order T004–T010K in `연구실행_체크리스트_20261011.csv`.

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
python adapt_compressed_n6_stateprep.py
python compile_adaptive_n6_qpu_evidence.py
python verify_adaptive_n6_qpu_evidence.py
```

The field values accepted by `fit_spectral_flow_states_n9.py` are 1.40, 1.45, and 1.50 T. The compiled scan must pass all nine field-temperature cases. The adaptive N=6 preparation must reconstruct all nine low-energy states and pass the frozen entropy-error target before the archived hardware evidence is interpreted.

## 6. IBM hardware safeguard

No new hardware job is required to reproduce the paper conclusions. Frozen sanitized records contain the complete paired nine-circuit runs on IBM Kingston: twirled TREX job `db6et9clf4us73c582f0` and twirled TREX+ZNE job `db6eurqmb58s738bmisg`. They used 38 and 113 seconds of Open Plan QPU time, respectively. Against the exact entropy correction of -0.041611 kB/spin, the jobs returned -0.021679 and -0.026822 kB/spin. Both resolve the negative sign; neither establishes quantitative magnitude. ZNE reduces absolute bias by 25.8% while increasing the propagated reported-error standard deviation by a factor of 4.60.

The hardware script `run_adaptive_n6_qpu.py`:

- refuses any account whose plan preference is not exactly `open`;
- checks remaining free usage;
- defaults to dry run;
- requires explicit `--submit`;
- records no token or cloud resource name.

Do not resubmit merely to duplicate the archived benchmark. The frozen stopping rule is `PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE`; no hardware value enters the bulk material or refrigerator calculations.

## 7. Interpretation boundaries

- The local and Sobol designs provide 631 exact M4 evaluations and initially retain 16 models. Three adaptive rounds add 192 exact evaluations and retain 69 models in total, spanning 0.03340–0.20211 T. The acquisition surrogate only chooses points; every reported value is an exact multistart calculation. The result is not a posterior distribution or a continuous-space exhaustion proof.
- The four lower-field fan points at 1.45–1.60 T were not used for retention. Their deterministic holdout error correlates with the critical-field shift at Spearman rho=0.940. The 23 models no worse than the published representative are a relative ranking only; the full 69-model envelope remains authoritative until raw delta(H), qz(H), and phase-boundary information are fitted jointly.
- The public 3.5 T INS interlayer scan was not used for diffraction retention or fan ranking. Its archived Jb=3 microelectronvolt calculation defines an 18.048 microelectronvolt indicative bandwidth reference. Thirty of 69 M4 models remain below it at the central fit and all matched one-standard-deviation exchange/g-factor corners; their shift range is 0.03340–0.14757 T. The direct linewidth audit calibrates `FWHM^2 = resolution^2 + c*twin_separation^2` to the archived calculation (0.000141 meV RMSE) and profiles `c` between the single-domain and equal-twin limits. All 69 M4 models then have chi-square no larger than the source curve, so the 30-model contraction is a sensitivity scenario rather than a calibrated rejection and the full 69-model range remains authoritative.
- The one-coordinate boundary audit closes the tested grid faces but does not exclude disconnected parameter regions or simultaneous multi-parameter extensions.
- Extended Data Fig. 6 sheet-a/b labels fail an identifiable article-caption consistency check and are excluded from new cross-material calculations; sheet-d time traces are used only for context.
- The N=12 calculation rejects rigid-profile transfer. The preregistered four-model N=18 cover gives optimum shifts of 0.293–0.425 T and peak-heat ratios of 0.537–0.852; 4 of four models cross the major-transfer threshold. Common-grid, close-peak, low-state, and cutoff checks pass inside N=18, which remains short of a thermodynamic-limit solution.
- The former independent shape-coordinate box is retained only as an adversarial destruction test. The primary design keeps peak shift, heat scale, and left/right widths correlated within each microscopic coverage model.
- The correlated-template pass fraction is deterministic coverage, not a probability or confidence interval.
- The size-stressed 1.292–1.680 T, 4 mT sweep controls peak-grid localization only; calibration, profile width, demagnetization, and thermodynamic measurement uncertainty remain experimental requirements.
- The exact 1.292/1.716 T panel is globally optimal among all 5,671 two-point panels on the 107-field 90 mK grid at a ±2.5% heat bound. Error fractions are deterministic resolution scenarios, and the design is conditional on 52 N=12 profiles plus 16 correlated four-model N=18 stresses.
- The equal-volume spheroid uses saturation magnetization as a conservative demagnetization bound. Exact package correction requires dimensions, orientation, packing fraction, and measured M(B,T).
- The 90–300 mK `W95` result is a set of measured field candidates.
- Conditional volumes include magnetic entropy only.
- The YbNi1.6Sn result validates the lumped heat-capacity-to-hold-time identity on an independent device; it does not transfer that device's addenda, switch, contact, or ramp losses to NBCP.
- N=9 state-preparation targets use exact classical eigenvectors.
- The complete N=6 QPU benchmark resolves the sign of the selected entropy correction but does not establish quantitative magnitude, bulk-material accuracy, or quantum advantage.
