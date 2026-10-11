# Release notes: v0.1.0

This is the first reproducible research release for the NBCP magnetocaloric model-risk and quantum spectral-flow study.

## Included evidence

- Public-input manifests, magnetic-lattice provenance, and SHA-256 hashes.
- M0--M4 finite-size, polarized-state, long-range dipolar, and interlayer-ablation results.
- Model-risk engineering outputs for optimal field, 95% operating set, conditional magnetic heat, and regret.
- Ideal N=6/9/12 quantum-algorithm screens and N=9 field-dependent state-preparation results.
- Two IBM Fez Open Plan pilot records: TREX job `db5r8184qg6s73c3sl0g` and TREX+ZNE job `db5r8nimb58s738atca0`, totaling 71 seconds of QPU usage.
- Main manuscript draft, Supplementary Information, references, and machine-checked claim audit.

## Frozen conclusions

- The representative interlayer ablation changes the inferred 90 mK optimum by 0.0815--0.0953 T and produces 8.80--10.00% design regret, exceeding the preregistered 0.05 T or 5% engineering gate.
- Direct public entropy data select 1.5 T at 90--120 mK, 1.3 T at 150--200 mK, and 0 T at 250--300 mK among the measured candidate fields.
- Ideal spectral-flow and state-preparation tests pass at the declared small-system sizes.
- The hardware runs demonstrate an execution path. Residual bias and ZNE variance keep the quantitative QPU result on hold.

## Scope

No quantum advantage is claimed. The cooling-heat values are conditional magnetic entropy budgets rather than complete refrigerator-cycle predictions. The P2 interlayer network is representative and non-unique, so microscopic attribution remains on hold.

## Verification

Follow `REPRODUCE.md`. The fast claim checks require no QPU submission. External publisher archives and third-party repositories are fetched from their frozen public URLs and retain their original licenses.
