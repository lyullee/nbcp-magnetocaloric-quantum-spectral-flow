# Release notes: v0.2.0

This release completes the Cryogenics-targeted finite-size, quantum-hardware, and package-shape audit.

## New evidence

- A preregistered four-model N=18 M4 cover rejects universal quantitative transfer from N=12. The N=18 optimum shifts span 0.293–0.425 T and peak-heat ratios span 0.537–0.852; all four models cross the declared major-transfer threshold.
- Four N=18 templates expand the measurement design to 68 correlated scenarios. Exhaustive search over all 5,671 two-point panels selects 1.292 and 1.716 T, leaving 4.246% worst regret at a ±2.5% heat bound.
- A package-demagnetization envelope gives a conservative 0.136 T correction for the published plate geometry and requires c/a of at least 2.72 to bound the correction below 0.02 T.
- Two complete IBM Kingston Open Plan jobs each execute all nine adaptive six-qubit circuits. Twirled TREX job `db6et9clf4us73c582f0` uses 38 s; twirled TREX+ZNE job `db6eurqmb58s738bmisg` uses 113 s.
- Both hardware runs resolve the exact negative entropy-correction sign. ZNE reduces absolute bias by 25.8%, from 0.019932 to 0.014789 kB/spin, while increasing the propagated reported-error standard deviation by a factor of 4.60. The quantitative magnitude remains on hold.
- The release includes the 22-page Word manuscript, 19-page Word Supplementary Information, 11 main figures, 4 supplementary figures, 24 cited references, and machine-readable verification records.

## Scope

No quantum advantage is claimed, and no QPU value enters the bulk material or refrigerator calculations. The QPU result is a measured method-accuracy gate. N=18 remains a finite cluster rather than a bulk solution. Cooling-heat values are conditional magnetic entropy budgets, and final package performance still requires measured heat profiles, geometry, magnetization, addenda, contacts, switching, ramp loss, and warm-side regeneration.

## Verification

`verify_adaptive_n6_qpu_evidence.py`, `verify_cryogenics_submission.py`, and `verify_cryogenics_docx.py` all return `PASS`. The frozen QPU stopping decision is `PASS_SIGN_RESOLUTION_HOLD_QUANTITATIVE_MAGNITUDE`. Follow `REPRODUCE.md`; no new QPU submission is required.
