# NBCP magnetocaloric model risk and quantum spectral flow

This repository contains a public-data study of how weak interactions in the frustrated quantum refrigerant Na2BaCo(PO4)2 change sub-kelvin magnetic-refrigeration decisions, followed by a resource-gated quantum-algorithm validation.

The study asks whether long-range dipolar coupling, bond-dependent exchange, and interlayer exchange move the optimal terminal field enough to matter in a refrigerator. Quantum circuits are tested only after the classical engineering gate is passed.

## Main findings

- A nearest-neighbor dipolar truncation shifts the polarized critical field by about 0.0825 T, but the full Ewald sum nearly cancels that shift under the tin-foil convention. Local dipoles are therefore a poor bulk proxy.
- A deterministic 375-member interlayer grid is constrained by the public near-saturation propagation vector. Seven models remain, with critical-field ablations of 0.06832–0.16365 T relative to the no-interlayer model.
- Under the central heat profile, all seven retained models exceed the 0.05 T field-shift and 5% regret thresholds. Across a broader 2,268-case shape envelope, 88.1% pass both; every declared deformation passes once the ablation reaches 0.10850 T.
- The actual-material 90 mK minimax field is 1.4052 T with 9.21% worst regret. No common 95% field exists under the full shape envelope, but constraining the peak position to ±0.002 T restores one. A no-fit target is a 26-point sweep from 1.376 to 1.476 T at 0.004 T spacing.
- At 90 mK, a 0.10 microWatt magnetic-only load held for 60 min requires 0.959 cm3 or 4.03 g in the worst declared scenario.
- Direct public entropy curves select 1.5 T at 90–120 mK, 1.3 T at 150–200 mK, and 0 T at 250–300 mK among the measured candidate fields.
- A low-energy spectral-flow estimator passes all N=9 and N=12 ideal entropy tests. Physics-informed N=9 state-preparation circuits pass at 1.40, 1.45, and 1.50 T.
- Two N=6 IBM Fez jobs used 71 seconds of free Open Plan QPU time. They prove the execution path, while remaining bias and ZNE variance prevent a quantitative hardware claim.

No quantum advantage is claimed. The material and engineering conclusions use classical calculations and public experimental data. The hardware run is a scoped circuit-path pilot.

## Manuscript

- [Main manuscript](manuscript/Manuscript_draft_20261012.md)
- [Supplementary Information](manuscript/Supplementary_Information_draft_20261012.md)
- [References](manuscript/references.bib)
- [Coauthor review checklist](manuscript/COAUTHOR_REVIEW_CHECKLIST.md)
- [Completion audit](연구목표_완료감사_20261012.md)

## Repository map

- `sources/`: input manifests and magnetic-lattice provenance. Large public originals are downloaded locally and are not redistributed.
- `data/`: compact digitized or derived inputs used by the calculations.
- `results/`: machine-readable JSON and CSV outputs, including stopping decisions.
- `figures/`: generated figures.
- `manuscript/`: article draft, supplement, references, and author metadata.
- `octave/`: independent Octave cross-check code. The 3.7 GB local Octave runtime is excluded.
- top-level Python files: calculation, compilation, validation, and audit scripts.

## Quick verification

Python 3.14 was used for the frozen run. Create an isolated environment and install:

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

Run the fast frozen-result checks:

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

The commands do not submit a quantum job. `run_open_plan_n6_pilot.py` defaults to a local dry run and requires the explicit `--submit` flag plus an Open-only account guard.

See [REPRODUCE.md](REPRODUCE.md) for public-input download and longer calculation paths.

## Data scope

The repository preserves hashes and official URLs for every external input. Publisher workbooks, article source archives, and third-party software checkouts are not relicensed here. `fetch_sources.py` and `fetch_public_inputs.py` download or verify them from their public locations.

The reported cooling heat is a conditional magnetic entropy budget. It excludes nuclear and lattice heat capacity, addenda, heat switches, wiring, contact resistance, magnetic-ramp loss, and thermal leakage.

## License and citation

Original code in this repository is released under the MIT License. Original text, tables, figures, and derived data are released under CC BY 4.0 unless a file identifies a different upstream license. External source files retain their original terms.

Citation metadata are provided in `CITATION.cff`, and Zenodo deposit metadata are provided in `.zenodo.json`. The archival DOI will be inserted after coauthor confirmation and the author-approved GitHub release are connected to Zenodo.
