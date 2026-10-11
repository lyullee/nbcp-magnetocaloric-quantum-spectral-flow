# Supplementary Information

## Propagation-vector-constrained sub-kelvin refrigerator design with a frustrated quantum magnet and a quantum spectral-flow route

Ugwiyeon Lee, Seong-Yun Bae, Hakrae Kim, and Jungsu Lee

## S1. Reproducibility policy

The study was run as a sequence of gates. Later calculations were allowed only after a recorded earlier gate passed. Thresholds were fixed before the P0–P2 engineering calculation:

- engineering GO if \(|\Delta B^*|\ge0.05\) T or regret ≥5%;
- engineering STOP if \(|\Delta B^*|<0.02\) T and regret <2%;
- ordinary quantum entropy target: absolute error ≤0.001 \(k_B\) per spin;
- sign-change quantum target: absolute error ≤0.0002 \(k_B\) per spin;
- inferred field uncertainty target: \(\sigma_B\le0.02\) T.

The work log is `연구실행_체크리스트_20261011.csv`. The final hash inventory is `results/bulk_path_manifest_20261011.json`. Each negative route was retained rather than deleted so that algorithm selection can be audited.

The interlayer parameter-grid and finite-temperature shape-envelope analyses were added as reviewer-oriented robustness tests after the original P0–P2 result. Their coordinate values were frozen before each Cartesian screen was run. They are deterministic sensitivity bounds, not registered confidence intervals or posterior samples.

## S2. Public sources and provenance

### S2.1 Experimental thermodynamics

The following source-data workbooks associated with Xiang *et al.* were downloaded from the publisher and preserved unchanged:

| Local file | Use | SHA-256 |
|---|---|---|
| `sources/figure_2.xlsx` | demagnetization and heat-capacity source data | `6c0ff63711b67815e368a229a2ec53a60a338b7728f4b5cf1eb7657907c42c27` |
| `sources/extended_figure_4.xlsx` | entropy and correction curves | `ee95f89f71d58c0ef06d16ed14afa6b5261a1c0018334f45c9bf902a16897805` |
| `sources/extended_figure_6.xlsx` | refrigeration comparisons and hold-time source data | `26a31b0d22ff37f808763679b9aab3423d6eb6b68fb727d3e6cc02dfa068dbe7` |

The source workbooks contain central-value curves, not a complete set of original synchronized instrument logs. No confidence interval was invented where the public file did not supply one.

### S2.2 Magnetic structure and Hamiltonian

The 2026 arXiv source archive was preserved with SHA-256 `8e6b98abf97340684470416ab4f3229993548c10677cf66b3a047eaa17881b13`. The magnetic-lattice input is documented in `sources/crystal_input_provenance_20261012.json`. A public full atomic CIF was not found in the preserved archives. A full CIF is unnecessary for the frozen spin model because the Ewald and exchange calculations require only the Co magnetic sublattice, lattice vectors, bond displacements, and g tensor.

## S3. Hamiltonian conventions

### S3.1 Nearest-neighbor exchange

The M0 Hamiltonian is

\[
H_0=J_{zz}\sum_{\langle ij\rangle}\left[S_i^zS_j^z+\Delta(S_i^xS_j^x+S_i^yS_j^y)\right]-h\sum_iS_i^z,
\quad h=g_{zz}\mu_BB.
\]

The parameters are \(J_{zz}=0.1225\) meV, \(\Delta=0.636\), and \(g_{zz}=4.716\). The M1 correction is

\[
H_{\pm\pm}=\sum_{\langle ij\rangle}\left(J_{\pm\pm}\gamma_{ij}S_i^+S_j^++J_{\pm\pm}^*\gamma_{ij}^*S_i^-S_j^-\right),
\]

with \(J_{\pm\pm}=-1.9\) microelectronvolt and \(\gamma_{ij}=1,e^{i2\pi/3},e^{-i2\pi/3}\) for the three positive triangular bonds. \(J_{z\pm}\) is zero in the representative parameterization used here.

### S3.2 Dipolar interaction

For magnetic moment \({\bf m}_i=\mu_B{\bf g}{\bf S}_i\),

\[
H_{\rm dip}=\sum_{i<j}\frac{\mu_0}{4\pi r_{ij}^3}
\left[{\bf m}_i\cdot{\bf m}_j-3({\bf m}_i\cdot\hat{\bf r}_{ij})({\bf m}_j\cdot\hat{\bf r}_{ij})\right].
\]

M2 retains only the nearest-neighbor tensor. M3 and M4 use the periodic Ewald tensor with \(\alpha=0.22\) Å−1, real-space cutoff 42 Å, and reciprocal cutoff 8 Å−1 in the frozen production calculation. Tin-foil omits the \(k=0\) surface term. Spherical-vacuum adds \(4\pi I/(3V)\). Convergence checks are stored in `results/bulk_stability_convergence.json`.

### S3.3 Interlayer paths

The AA-stacked triangular convention is implemented in `stacked_lattice.py`. Each site has one vertical \(J_2\) path and two symmetry-inequivalent sets of three diagonal paths, \(J_a\) and \(J_b\), to the next layer. Reversing the crystallographic \(c\) direction exchanges the labels \(J_a\) and \(J_b\) without changing the spectrum.

P1 is a sensitivity model with one isotropic vertical bond above and below and \(J_\perp=0.05\) K. The published P2 seed uses \(J_2=2.8\) microelectronvolt, \(J_a=0\), \(J_b=7.9\) microelectronvolt, and \(\Delta_c=0.08\).

The set-valued P2 screen uses:

- \(J_2=\{1.4,2.8,4.2\}\) microelectronvolt;
- \(J_a=\{-2.8,-1.4,0,1.4,2.8\}\) microelectronvolt;
- \(J_b=\{3.95,5.925,7.9,9.875,11.85\}\) microelectronvolt;
- \(\Delta_c=\{0.04,0.06,0.08,0.10,0.12\}\).

This gives 375 deterministic combinations. All exchange magnitudes remain below 10% of \(J_{zz}\), and every \(\Delta_c\) remains strongly Ising-like. The ranges are a local robustness envelope around the published representative set, not confidence limits.

## S4. Public entropy processing and engineering metrics

The workbooks were read sheet by sheet without changing the raw files. Duplicate points were collapsed only when their entropy values agreed within \(10^{-9}\). Temperature interpolation is piecewise linear and restricted to the published range. The direct candidate table contains no interpolation between magnetic fields.

The hot-state entropy is obtained from the 4 T curve at 2 K. For each cold endpoint,

\[
Q=T(S_{\rm cold}-S_{\rm hot}).
\]

Negative \(Q\) means that the endpoint cannot absorb heat relative to the selected hot-state entropy under this endpoint accounting. The 90–300 mK direct results are stored in:

- `results/temperature_candidate_comparison.csv`;
- `results/temperature_candidate_ranking.csv`;
- `results/engineering_case_refinement.json`.

The reconstructed continuous reference profile comes from top-down integration of public XTRG \(C_m/T\) data. Its field dependence at the 1 K anchor is supplied by exact N=9/N=12 entropy differences. This hybrid reconstruction is used only as a decision-scale reference and is validated independently at zero field.

## S5. Finite-size and low-magnon diagnostics

### S5.1 Matrix-free validation

The matrix-free implementation was checked against explicitly constructed Pauli matrices at N=9 and N=12. Maximum elementwise discrepancies were below numerical roundoff. Spin conventions, bond multiplicities, and periodic cells are unit-tested in `test_matrix_free_hamiltonian.py`, `test_stacked_lattice.py`, and `test_dipole_ewald.py`.

### S5.2 FTLM and TPQ tests

FTLM follows Jaklič and Prelovšek. Random trace vectors are propagated through Krylov spaces and combined with exact low-energy states. TPQ calculations provide an independent typical-state check. The following identities or limits were verified:

1. \(S/N\rightarrow\ln2\) at high temperature;
2. \(C=dU/dT\) within the finite-difference tolerance;
3. exact diagonalization agreement at N=9 and N=12;
4. random-vector convergence at N=15;
5. deterministic low-energy replacement at N=18.

At 90 mK and 1.50 T, the N=18 full estimator and the lowest 16 deterministic levels differ by \(1.45\times10^{-4}\;k_B\) per spin for the M3–M0 correction. This is smaller than the N=15/N=18 disagreement and confirms that shell effects dominate stochastic error.

### S5.3 Fixed-magnon truncation

Near saturation, the Hilbert space was partitioned by the number of down spins for M0 and by even/odd parity for M3. N=18 through N=27 calculations progressively included low-magnon sectors. For N=24, the first omitted parity-compatible sectors have boundary residuals up to 0.056 K, but the minimum isolated energy denominator is 2.5816 K. A statewise \(r^2/\Delta E\) diagnostic gives a Boltzmann-weighted shift of 0.586 mK and a 90 mK entropy influence of \(8.72\times10^{-5}\;k_B\) per spin.

The size sequence of the M3–M0 correction at 90 mK alternates in sign: N=15 positive, N=18 negative, N=21 positive, N=24 negative, and N=27 negative but closer to zero. This is why the finite clusters are used as a diagnostic rather than a bulk extrapolation.

## S6. Polarized-state stability

For a fully polarized state, the quadratic transverse fluctuation matrix is constructed in momentum space and its lowest eigenvalue is minimized over the Brillouin zone. The critical field is where the lowest mode reaches zero. The M0 analytic result is reproduced to \(4.44\times10^{-16}\) T.

The frozen values are:

| Model | Boundary | \(H_c\) (T) | Critical wavevector (r.l.u.) |
|---|---|---:|---|
| M0 | tin foil | 1.7743602 | (1/3, 2/3, 0) |
| M1 | tin foil | 1.7751632 | (0.34050, 0.65950, 0) |
| M2 | tin foil | 1.8568211 | (0.73061, 0.36530, 0) |
| M3 | tin foil | 1.7723974 | (0.28357, 0.64179, 0.50000) |
| M4 | tin foil | 1.8676712 | (0.71632, 0.35816, 0.66571) |
| M3 | sphere vacuum | 1.8260773 | same finite-k minimum |
| M4 | sphere vacuum | 1.9213511 | same finite-k minimum |

The surface term changes the uniform contribution and absolute field but leaves the P2–P0 ablation at 0.0952738 T.

### S6.1 Public-propagation-vector screen

Near saturation, Ref. [4] reports \(\delta_c=0.027(2)\) along a K–M boundary and \(q_z=0.328(2)\). The same source approximates \(q_z\) by \(1/3\) in its interlayer model. Before the grid calculation, the acceptance bands were fixed to \(\delta_c=0.021\)–0.033, folded \(q_z=0.322\)–0.3373, and a K–M line residual no larger than 0.003 reciprocal-lattice units. Seven of 375 models pass:

| \(J_2\) | \(J_a\) | \(J_b\) | \(\Delta_c\) | \(\Delta H_c\) vs M3 (T) | \(\delta_c\) | folded \(q_z\) | \((00l)\) bandwidth (microelectronvolt) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.4 | 0.0 | 5.925 | 0.10 | 0.068321 | 0.024893 | 0.329831 | 3.979 |
| 1.4 | 2.8 | 11.85 | 0.04 | 0.163650 | 0.025245 | 0.327820 | 3.772 |
| 2.8 | −2.8 | 11.85 | 0.06 | 0.108503 | 0.024450 | 0.332444 | 3.738 |
| 2.8 | 0.0 | 5.925 | 0.12 | 0.073709 | 0.024754 | 0.325444 | 5.082 |
| 2.8 | 0.0 | 7.90 | 0.08 | 0.095274 | 0.024826 | 0.334289 | 4.384 |
| 4.2 | −1.4 | 11.85 | 0.06 | 0.128791 | 0.024600 | 0.329369 | 4.410 |
| 4.2 | 1.4 | 7.90 | 0.08 | 0.115488 | 0.025037 | 0.328267 | 5.280 |

Exchange entries are in microelectronvolt. The representative set reproduces the paper's 4.4 microelectronvolt bandwidth to 0.016 microelectronvolt. Over the full designed grid, the critical-field ablation correlates most strongly with \(J_b\) (descriptive Spearman \(\rho=0.824\)); this is a sensitivity diagnostic, not a fitted causal coefficient.

## S7. Critical-field transport and robust optimization

Let \(Q_0(B,T)\) be the public reference curve with source critical field \(H_{c,0}\). For a target critical field \(H_c\), the two mappings are:

\[
Q_{\rm red}(B,T;H_c)=Q_0\left(B\frac{H_{c,0}}{H_c},T\right),
\]

\[
Q_{\rm off}(B,T;H_c)=Q_0\left(B-H_c+H_{c,0},T\right).
\]

The mappings bracket two simple assumptions about how the near-critical profile moves. They do not alter \(Q_{\max}\), which is why the manuscript calls the inherited heat scale a surrogate.

The 2026 full-model mean-field field overestimates the measured field. P2 is therefore set to either 1.65 or 1.70 T and the calculated shifts are subtracted to obtain P0 and P1. Under all four anchor/mapping combinations, P2–P0 passes both engineering thresholds.

The robust-field grid spacing is 0.00005 T in the final minimax calculation. The identity \(E=P t\) was checked for every row of the 24-condition load/hold map.

### S7.1 Finite-temperature shape envelope

The rigid mappings are augmented by four deformation coordinates. Relative to the mapped peak, the low- and high-field sides are independently rescaled by \(s_-,s_+\in\{0.8,1.0,1.2\}\), the peak is displaced by \(\Delta B_{\rm peak}\in\{-0.02,0,0.02\}\) T, and the entire heat curve is scaled by \(a_Q\in\{0.9,1.0,1.1\}\). The 10% heat and 20% width ranges cover the N=9/N=12 reconstruction differences at 90 mK. The residual peak bound equals twice the heatmap field pixel and twice the observed N=9/N=12 peak-location difference at 120 mK.

Combining seven retained interlayer shifts, two measured \(H_c\) anchors, two mappings, and 81 deformations gives 2,268 gate cases at each temperature. At 90 mK, 1,998 cases pass both the 0.05 T and 5% gates. The pass count is coverage of a Cartesian stress design and has no probability interpretation.

| \(\Delta H_c\) (T) | Cases | Fraction passing both | Minimum \(\Delta B^*\) (T) | Minimum regret (%) |
|---:|---:|---:|---:|---:|
| 0.068321 | 324 | 0.556 | 0.03844 | 2.24 |
| 0.073709 | 324 | 0.667 | 0.04305 | 2.68 |
| 0.095274 | 324 | 0.944 | 0.06150 | 4.39 |
| 0.108503 | 324 | 1.000 | 0.07281 | 5.47 |
| 0.115488 | 324 | 1.000 | 0.07879 | 6.22 |
| 0.128791 | 324 | 1.000 | 0.09017 | 7.99 |
| 0.163650 | 324 | 1.000 | 0.11998 | 10.98 |

The smallest retained shift that passes every declared deformation is therefore 0.108503 T. The central deformation \((\Delta B_{\rm peak},s_-,s_+,a_Q)=(0,1,1,1)\) passes for every retained model.

### S7.2 Density and design conversion

Using \(a=5.312\) Å, \(c=6.983\) Å, one formula unit per \(P\bar 3\) primitive cell, and a formula mass of 432.1793 g mol\(^{-1}\), the crystallographic density is 4.2056 g cm\(^{-3}\). The reported 1.65 g NBCP sample in the public Extended Data Fig. 6b workbook therefore corresponds to 0.3923 cm\(^3\). Required magnetic mass is \(m=\rho P_{\rm load}t_{\rm hold}/Q\). Packing and nonmagnetic components are not included.

### S7.3 Minimum near-saturation field sweep

The full shape envelope has no common 90 mK \(W_{95}\), so we varied the residual symmetric uncertainty in peak position while retaining both measured critical-field anchors, both transport mappings, and all nine independent left/right width combinations. The largest tested uncertainty with a common \(W_{95}\) was ±0.002 T. At ±0.002 T, the minimax field is 1.4106 T, the worst regret is 4.97%, and the common interval is 1.41035–1.41070 T. At ±0.00225 T the common interval disappears on the 0.00005 T decision grid.

Under the original ±0.02 T displacement envelope, the possible 90 mK peak positions span 1.376405–1.474155 T. Sampling 1.376–1.476 T every 0.004 T gives 26 points and guarantees a grid point within 0.002 T of any peak in that support without fitting a line shape. This is only a discretization target. Field calibration, thermal equilibration, entropy reconstruction, and measurement noise must be controlled separately. The scan and field list are stored in `results/peak_location_uncertainty_scan.csv` and `results/recommended_90mK_field_sweep.csv`.

## S8. Quantum-algorithm screening

The following routes were evaluated under a common final entropy-error target:

| Route | Main failure or result | Decision |
|---|---|---|
| direct temperature finite difference | cancellation amplifies measurement variance; up to \(3.87\times10^9\) shots | STOP |
| covariance/thermodynamic integration | too many Hamiltonian covariance terms at required precision | STOP for Open Plan |
| Jarzynski ramp | exponential-work variance and ramp depth | STOP |
| variational Gibbs state | low-temperature conditioning and auxiliary resources | STOP |
| QITE–QMETTS | N=9 adaptive three-body tangent reduces residual but reaches up to 1057 CX per step | STOP on current QPU |
| low-energy spectral flow | passes N=9/N=12 ideal accuracy and admits observable-specific circuits | SELECT |

The screening summary and Pareto evidence are in `results/quantum_algorithm_screening_summary.json` and `results/quantum_algorithm_second_screen_summary.json`.

## S9. Spectral-flow error accounting

For each field and temperature, total entropy error includes:

1. low-energy truncation and certified thermal tail;
2. Gauss–Legendre integration error;
3. state-preparation bias in \(\langle V\rangle\);
4. propagated energy error in the partition function;
5. sampling target where a shot budget is reported.

The strict sign-change point uses a smaller absolute tolerance because the physical correction is near zero. All N=9 and N=12 exact-state tests preserve the sign. The N=9 circuit-prepared scan also preserves the sign in all nine field-temperature combinations.

At 1.40 T, the hardest degenerate pair has fidelity 0.97636. Its \(V\) expectation and final entropy remain within target. At 1.45 T, the minimum fidelity is 0.99935 and the maximum entropy error is \(1.51\times10^{-4}\;k_B\) per spin. At 1.50 T, the minimum fidelity is 0.99806 and maximum entropy error is \(5.67\times10^{-4}\;k_B\) per spin.

## S10. QPU preflight and jobs

The submission script performs the following guards before any job:

1. load only the named saved account;
2. require `plans_preference=['open']` exactly;
3. require at least 30 remaining QPU seconds;
4. select an operational non-simulator backend;
5. transpile every circuit before submission;
6. record only sanitized backend, job, usage, circuit, observable, and error data.

The N=6 live preflight on IBM Fez gave a maximum of 84 CZ gates, minimum calibration-product heuristic 0.8049, and maximum circuit duration divided by median T1 of 0.03986. The two completed jobs were:

| Mitigation | Job ID | Estimate (K) | Reported standard error (K) | Absolute error (K) | QPU seconds |
|---|---|---:|---:|---:|---:|
| TREX | `db5r8184qg6s73c3sl0g` | 0.0209552 | 0.0018476 | 0.0277818 | 19 |
| TREX + twirling + ZNE | `db5r8nimb58s738atca0` | 0.0311438 | 0.0158933 | 0.0175932 | 52 |

The exact value is 0.04873698 K. TREX is inconsistent with the exact value at two reported standard errors. The ZNE interval contains the exact value within two standard errors, but the uncertainty is too large for a thermodynamic claim. No additional QPU execution is justified until the N=9 ECR count is reduced and a stronger noisy preflight predicts target-level accuracy.

![N=6 Open Plan hardware pilot](../figures/open_plan_qpu_pilot_N6.png)

**Figure S1.** (a) IBM Fez estimates compared with the exact N=6 observable. TREX alone leaves a large systematic discrepancy; ZNE reduces the central discrepancy but increases variance. (b) QPU time. The two jobs used 71 free seconds. This figure is supplementary because neither hardware estimate enters the material or refrigerator result.

## S11. Verification suite

The frozen verification commands are:

```text
python verify_engineering_decision.py
python verify_quantum_stateprep_qpu.py
python verify_second_screen.py
python test_quantum_local_benchmark.py
python test_bulk_polarized_stability.py
python test_interlayer_robustness.py
python screen_interlayer_parameter_ensemble.py
python stress_test_thermodynamic_shape.py
python compile_refrigerator_design_chart.py
python make_bulk_manifest.py
```

All passed on Python 3.14.8 with NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2 under Windows 11. Independent Octave 11.3.0 checks reproduced the M0–M4 polarized-state calculation and the N=9 exact thermodynamic regression. The public QMagen low-bond-dimension N=54 route was stopped because it failed convergence and log-partition-function convexity tests.

## S12. Scope table for reviewers

| Statement | Supported? | Evidence boundary |
|---|---|---|
| every public-\(q\)-constrained interlayer model shifts the central-profile optimum by more than 0.05 T | yes | seven retained models, two mappings, and two measured \(H_c\) anchors |
| every retained model passes both gates under every declared shape deformation | no | 88.1% of Cartesian cases pass; all-shape guarantee begins at \(\Delta H_c=0.108503\) T |
| P2 is the unique microscopic Hamiltonian of NBCP | no | source explicitly calls it representative and non-unique; seven stress-grid models are retained |
| full dipolar summation is needed before judging nearest-neighbor dipoles | yes | M2 versus M3 bulk stability and Ewald convergence |
| the calculated volumes predict a finished refrigerator | no | magnetic-only reversible budgets |
| spectral flow is more sample-efficient than direct finite-difference entropy for the tested ideal instances | yes | common error target and shot estimates |
| a quantum advantage has been demonstrated | no | classical targets and N=6 noisy pilot |
| IBM hardware gives a quantitatively reliable NBCP correction | no | N=6 bias and variance; N=9 not run |
