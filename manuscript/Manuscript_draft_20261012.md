# Propagation-vector-constrained sub-kelvin refrigerator design with a frustrated quantum magnet and a quantum spectral-flow route

**Ugwiyeon Lee**<sup>1,*</sup> (ORCID: 0009-0008-8976-5363), **Seong-Yun Bae**<sup>1</sup> (0009-0003-8731-428X), **Hakrae Kim**<sup>1</sup> (0009-0002-5048-3568), and **Jungsu Lee**<sup>1</sup> (0009-0006-9347-453X)

<sup>1</sup>Korea Gas Safety Corporation, Republic of Korea
<sup>*</sup>Corresponding author: lyullee@kgs.or.kr

Draft prepared in a *Physical Review Applied*-style structure. Numerical values and scope statements are frozen to the reproducibility manifest dated 12 October 2026.

## Abstract

Frustrated quantum magnets can provide large entropy changes below 1 K, but a refrigerator is controlled through a small set of engineering decisions: terminal field, field tolerance, and refrigerant inventory. We quantify how uncertain weak interactions change those decisions for the spin-1/2 triangular-lattice refrigerant Na2BaCo(PO4)2 using only public thermodynamic and diffraction data. Small-cluster finite-temperature calculations show alternating shell effects and cannot determine the bulk entropy correction. A thermodynamic-limit polarized-state calculation instead shows that a nearest-neighbor dipolar truncation produces a 0.0825 T critical-field shift that is almost cancelled by the long-range Ewald sum under tin-foil conditions. We then replace a single published interlayer model by a deterministic 375-member stress grid. Public near-saturation propagation-vector constraints retain seven models. Their M4–M3 critical-field ablations span 0.0683–0.1636 T; under the central transported heat profile, all retained models move the 90 mK optimum by 0.0584–0.1636 T and cause 5.13–20.44% regret, crossing predeclared 0.05 T and 5% engineering thresholds. A broader 2,268-case finite-temperature shape envelope, covering residual peak displacement, asymmetric width, and heat-scale changes, preserves both thresholds in 88.1% of cases and in every case once the critical-field ablation exceeds 0.1085 T. For the measured 1.65–1.70 T transition interval, the full envelope gives a 90 mK minimax field of 1.405 T and 9.21% worst-case regret; it also shows that no field is guaranteed to remain within 95% of optimum without denser near-saturation thermodynamic data. A deterministic measurement-design calculation identifies the missing specification: the 90 mK heat-profile peak must be localized to within ±0.002 T. Direct public entropy curves independently select 1.5 T at 90–120 mK, 1.3 T at 150–200 mK, and 0 T at 250–300 mK. Finally, a low-energy quantum spectral-flow estimator passes ideal N=9 and N=12 error tests, while two N=6 IBM Fez executions retain too much bias for thermodynamic use. The study converts microscopic model ambiguity into a bounded refrigerator decision, a minimum measurement target, and a hardware-aware quantum route.

## I. Introduction

Magnetic refrigeration below 1 K depends on entropy stored in a magnetic material. In an adiabatic-demagnetization refrigerator (ADR), the useful material properties must ultimately be converted into stage field, refrigerant mass or volume, and hold-time requirements [12–15]. Frustration and quantum fluctuations can retain unusually large entropy at low temperature, which makes frustrated magnets promising refrigerants. Na2BaCo(PO4)2 (NBCP) is a particularly strong example. Public measurements report a giant magnetocaloric response and cooling below 100 mK, together with an easy-axis triangular-lattice spin-1/2 description [1–4]. Recent spectroscopy further resolves field-dependent Kramers-doublet excitations and a transverse-field evolution from a continuum to sharp magnons [16,17].

The material question and the engineering question are different. A microscopic study asks which spin Hamiltonian produces the observed phases. A refrigerator designer must choose a terminal field, determine how accurately that field must be reached, and size a material inventory for a heat load. A model can reproduce broad thermodynamics and still move the narrow low-temperature optimum enough to change these decisions. This risk is greatest near a critical field, where a small energy-scale error can translate into a large loss at a fixed operating point.

The commonly used nearest-neighbor XXZ Hamiltonian captures much of the NBCP phenomenology [1–4]. Recent neutron diffraction, however, found a field-tuned incommensurate fan phase immediately below saturation and showed that dipolar coupling, bond-dependent exchange, and frustrated interlayer paths are required to describe the full sequence of ordered states [4]. These terms are small relative to the main exchange, but their effect on a refrigerator operating field has not been separated from their role in explaining magnetic order.

This work asks four linked questions.

1. Does a realistic interaction change the optimal terminal field by at least 0.05 T or cause at least 5% loss when a simpler-model field is used?
2. Does the conclusion survive finite-size diagnostics and uncertainty in the long-range dipolar boundary convention?
3. What field choices and magnetic entropy budgets follow directly from public measurements between 90 and 300 mK?
4. If the decision depends on a small free-energy correction, can a quantum algorithm estimate that correction with a lower sampling burden than direct entropy evaluation, and can the circuit path be demonstrated on free-access hardware?

We answer these questions with a staged calculation. Classical calculations determine whether the engineering effect is large enough to justify any quantum work. Only after the predeclared engineering threshold is passed do we test quantum algorithms. This order prevents a small hardware demonstration from becoming the reason for selecting the scientific problem.

Our main physical finding is that weak terms cannot be ranked from their local magnitude alone. The nearest-neighbor part of the dipolar interaction creates a large critical-field shift, but the full long-range sum can cancel it. Public propagation-vector data then reduce a 375-member interlayer stress grid to seven admissible near-saturation models and expose a quantitative boundary between a guaranteed design change and a data-limited one. Our main quantum result is similarly bounded: a low-energy spectral-flow estimator solves the ideal sampling problem for N=9 and N=12, while current N=6 hardware remains too biased for quantitative thermodynamics.

## II. Decision problem and model hierarchy

### A. Engineering quantities

We use the public 2 K, 4 T state as the hot reference. Its measured magnetic entropy density is

\[
S_{\rm hot}=6.6877263\;{\rm mJ\,cm^{-3}\,K^{-1}}.
\]

For a cold endpoint at temperature \(T\) and field \(B\), the conditional magnetic entropy budget is

\[
Q(B,T)=T\,[S(B,T)-S_{\rm hot}].
\]

This is the reversible magnetic contribution available during a cold hold. It is not cooling power, coefficient of performance, or a complete cycle simulation. The optimal measured or calculated field is

\[
B^*(T)=\arg\max_B Q(B,T),
\]

and the 95% operating set is

\[
W_{95}(T)=\{B:Q(B,T)\ge 0.95Q_{\max}(T)\}.
\]

For continuous reconstructed profiles, \(W_{95}\) is an interval. For direct public data at only 0, 1.3, 1.5, and 1.66 T, it is a set of measured candidates and must not be interpreted as a continuous field range.

If a field selected from model \(a\) is used under model or scenario \(b\), the design regret is

\[
R_{a\rightarrow b}(T)=100\frac{Q_b(B_b^*,T)-Q_b(B_a^*,T)}{Q_b(B_b^*,T)}\;\%.
\]

We fixed the decision thresholds before the final calculations: an interaction is engineering-relevant if \(|\Delta B^*|\ge0.05\) T or regret is at least 5%. A shift below 0.02 T and regret below 2% support stopping the quantum stage.

### B. M0–M4 hierarchy

The baseline model is the spin-1/2 nearest-neighbor easy-axis XXZ Hamiltonian for a field along the crystallographic \(c\) axis,

\[
H_0=J_{zz}\sum_{\langle ij\rangle}\left[S_i^zS_j^z+\Delta(S_i^xS_j^x+S_i^yS_j^y)\right]-g_{zz}\mu_BB\sum_iS_i^z.
\]

We use \(J_{zz}=0.1225\) meV, \(\Delta=0.636\), and \(g_{zz}=4.716\) from the public interaction analysis [4]. The model hierarchy is:

- **M0:** two-dimensional nearest-neighbor XXZ;
- **M1:** M0 plus the reported bond-dependent \(J_{\pm\pm}=-1.9\) microelectronvolt term;
- **M2:** M1 plus only the nearest-neighbor dipolar tensor;
- **M3:** M1 plus the full long-range dipolar interaction within and between AA-stacked layers, evaluated by a three-dimensional Ewald sum;
- **M4:** M3 plus the symmetry-allowed \(J_2,J_a,J_b\) interlayer exchange paths.

The magnetic lattice uses space group \(P\bar{3}\), \(a=5.312\) Å, \(c=6.983\) Å, and one Co magnetic site per primitive cell. The complete equivalent magnetic-lattice input and hashes are provided in `sources/crystal_input_provenance_20261012.json`.

For the interlayer uncertainty analysis, P0 denotes M3 without interlayer exchange. P1 adds an explicit simple vertical bond \(J_\perp=0.05\) K based on a published density-functional energy scale. The published representative P2 seed uses \(J_2=2.8\) microelectronvolt, \(J_a=0\), \(J_b=7.9\) microelectronvolt, and \(\Delta_c=0.08\). Because Ref. [4] explicitly calls this set non-unique, the final P2 analysis uses a deterministic local stress grid rather than treating the seed as a fitted Hamiltonian or probability distribution.

## III. Methods

### A. Public data and entropy reconstruction

All experimental inputs are public source-data workbooks associated with Ref. [1], public article source files, or a public diffraction data archive. The preserved raw files are hashed in the input manifests. We interpolate only in temperature within a published field curve. Direct 90–300 mK rankings do not interpolate across magnetic field.

For the continuous near-saturation diagnostic, we reconstructed the published XTRG entropy surface by integrating the recovered \(C_m/T\) field-temperature map downward from its 1 K entropy anchor. The zero-field reconstructed entropy agrees with the independently published vector curve with a root-mean-square error of 0.0335 J mol−1 K−1. We use the reconstructed surface as the M0 bulk reference profile; it is not relabeled as an M4 solution.

### B. Finite-size thermodynamics

We implemented matrix-free Hamiltonian-vector products and verified them against exact N=9 and N=12 matrices. N=15 and N=18 thermodynamics were computed with finite-temperature Lanczos sampling and thermal-pure-state checks [5,6]. Internal tests included the high-temperature \(\ln2\) entropy limit, \(dU/dT=C\), exact small-system agreement, random-vector convergence, and deterministic low-energy replacement.

The cluster calculations serve two roles. First, they determine whether a weak-term correction is visible at the target temperatures. Second, disagreement between sizes acts as a stopping signal against an unjustified bulk claim. Fixed-magnon-sector calculations through N=27 were added near saturation to separate truncation error from shell effects.

### C. Bulk polarized-state stability and dipolar summation

For each model, the zero-temperature polarized-state critical field was obtained from the minimum eigenvalue of the transverse Hessian over the three-dimensional Brillouin zone. The M0 implementation reproduces its analytic field to \(4.4\times10^{-16}\) T. Long-range dipoles were summed by Ewald decomposition. Because the \(k=0\) dipolar term depends on the dielectric or magnetic environment at the macroscopic boundary, conditionally convergent lattice sums require an explicit surface convention [18]. We therefore report both tin-foil and spherical-vacuum conventions. Real- and reciprocal-space cutoffs and the Ewald splitting parameter were converged independently.

The stability calculation is a bulk critical-mode diagnostic, not a finite-temperature entropy solver. To connect it to refrigeration, we horizontally transport the public M0 entropy profile by two decision-scale mappings: fixed reduced field \(B/H_c\) and fixed critical-field offset \(B-H_c\). Reporting both prevents a single similarity assumption from determining the conclusion.

The published full-model mean-field critical field is 1.86 T, whereas neutron diffraction and bulk probes place the measured transition at 1.65–1.70 T [4]. Consequently, the primary P0–P2 engineering analysis anchors P2 to the measured interval and transports only the calculated ablation shifts to P0 and P1. This removes the known absolute mean-field bias from the recommended operating field.

We test the non-uniqueness of P2 with a deterministic \(3\times5\times5\times5\) grid: \(J_2=1.4,2.8,4.2\) microelectronvolt; \(J_a=-2.8\) to \(2.8\) microelectronvolt; \(J_b=3.95\) to \(11.85\) microelectronvolt; and \(\Delta_c=0.04\) to 0.12. This is a local robustness envelope around the published representative set, not a statistical sample. For every point we calculate \(H_c\), the soft-mode wavevector, and the one-magnon bandwidth along \((00l)\). A model is retained only if its soft mode lies on the K–M line within 0.003 reciprocal-lattice units, its in-plane incommensuration lies in 0.021–0.033, and its folded \(q_z\) lies in 0.322–0.3373. These bands span the public 65 mK values \(\delta_c=0.027(2)\) and \(q_z=0.328(2)\) and the source model's documented replacement \(q_z\rightarrow1/3\).

### D. Robust operating fields and conditional sizing

For a scenario set \(\Theta\), the robust field minimizes the worst regret,

\[
B_{\rm robust}=\arg\min_B\max_{\theta\in\Theta}R(B;\theta).
\]

We distinguish three scenario sets. The central set contains the measured 1.65 and 1.70 T critical fields and the two transport mappings. The model set adds every public-\(q\)-constrained interlayer Hamiltonian. The finite-temperature shape set additionally applies a residual peak displacement of \(-0.02,0,+0.02\) T, independent left and right width factors of 0.8, 1.0, and 1.2, and heat-scale factors of 0.9, 1.0, and 1.1. The heat and width bounds cover the N=9 versus N=12 reconstruction spread; the peak bound is twice the recovered heatmap field pixel and twice the observed anchor-dependent peak shift at 120 mK. These coordinates are deterministic stress bounds and carry no probabilities.

Material volume for a prescribed magnetic-only hold is

\[
V=\frac{P_{\rm load}t_{\rm hold}}{Q},
\]

with consistent conversion of microWatts, seconds, and millijoules. These volumes exclude nuclear and lattice heat capacities, sample holder, wiring, heat switches, contact resistance, and irreversible ramp losses.

Mass is obtained from a crystallographic density of 4.206 g cm\(^{-3}\), calculated from the public \(P\bar 3\) primitive cell, one formula unit per cell, and the formula mass of Na\(_2\)BaCo(PO\(_4\))\(_2\). As a cross-check, the public 1.65 g demagnetization sample corresponds to 0.392 cm\(^3\).

### E. Quantum spectral-flow estimator

The quantum target is not the full entropy surface. It is the low-energy correction between a reference and a realistic Hamiltonian,

\[
H(\lambda)=H_{\rm M0}+\lambda V,\qquad V=H_{\rm M3}-H_{\rm M0}.
\]

For an instantaneous eigenstate \(|k,\lambda\rangle\), the Hellmann–Feynman relation gives

\[
\frac{dE_k}{d\lambda}=\langle k,\lambda|V|k,\lambda\rangle.
\]

We prepare only the thermally relevant low-energy states, measure \(V\) at three Gauss–Legendre nodes, integrate the derivative to recover the corrected levels, and calculate the low-rank partition function and entropy classically. Degenerate states with identical energy and derivative are grouped. Tail bounds are obtained from the omitted Boltzmann weight and a Fannes–Audenaert entropy bound.

This observable-specific route was selected after direct finite-difference entropy, covariance integration, Jarzynski estimation, variational Gibbs preparation, and QITE–QMETTS alternatives were evaluated under the same error target. At the N=9 sign-change point, direct finite differences required approximately \(3.87\times10^9\) shots. The selected low-energy route reduced the ideal estimate by orders of magnitude, although state search and noisy-hardware overhead remain.

N=9 physics-informed circuits use number-conserving Givens rotations plus a small sector-mixing layer required by \(J_{\pm\pm}\). The optimization target is an exact classical eigenvector, so this test measures representability and observable bias; it does not establish a scalable eigensolver. N=12 uses ideal low-energy spectral data and certified tail bounds.

### F. Hardware protocol

All hardware submissions were guarded to accept only an IBM account whose plan preference was exactly `open`. N=9 raw hardware execution was stopped after a fake-backend audit predicted excessive depth. A smaller N=6 open 3×2 fragment was compiled live to IBM Fez and executed twice: once with twirled readout error extinction (TREX), and once with TREX, gate twirling, and zero-noise extrapolation (ZNE). IBM documents an Open Plan allowance of up to 10 QPU minutes per rolling 28-day window and describes the applied mitigation levels [10,11]. Job identifiers and QPU seconds are frozen; credentials and cloud resource names are excluded.

## IV. Results

### A. Small clusters diagnose rather than remove the bulk problem

The N=15 and N=18 calculations detect a low-temperature correction, but they disagree in sign. At 90 mK and 1.50 T, \(S_{\rm M3}-S_{\rm M0}\) is +0.00664 \(k_B\) per spin for N=15 and −0.00331 \(k_B\) per spin for N=18. Aligning the cluster critical fields to the bulk instability does not remove the sign reversal. N=21, N=24, and N=27 low-magnon calculations continue an alternating size sequence whose magnitude decreases toward zero.

For N=18, a deterministic 16-level calculation agrees with the full finite-temperature estimate within \(2.42\times10^{-4}\;k_B\) per spin, showing that random trace noise is not responsible. N=24 boundary-residual and energy-denominator tests constrain the first omitted sectors to an estimated 90 mK entropy influence of \(8.72\times10^{-5}\;k_B\) per spin. The dominant uncertainty is therefore the finite-size shell sequence rather than the magnon cutoff.

This negative result sets the scope of the paper: finite clusters cannot determine the bulk sign of the local M3 entropy correction. We retain the public XTRG surface for absolute finite-temperature thermodynamics and use the independently converged bulk instability only to test whether the field scale changes enough to matter.

### B. A nearest-neighbor dipolar model gives the wrong hierarchy

The bulk critical fields are summarized in Table I.

**Table I. Thermodynamic-limit polarized-state critical fields.** Values are model-scale diagnostics, not direct operating fields.

| Model | Added physics | Tin-foil \(H_c\) (T) | Spherical-vacuum \(H_c\) (T) |
|---|---|---:|---:|
| M0 | nearest-neighbor XXZ | 1.77436 | 1.77436 |
| M1 | + \(J_{\pm\pm}\) | 1.77516 | 1.77516 |
| M2 | + nearest-neighbor dipole | 1.85682 | 1.85682 |
| M3 | + full long-range dipole | 1.77240 | 1.82608 |
| M4 | + representative interlayer exchange | 1.86767 | 1.92135 |

The \(J_{\pm\pm}\) contribution to the critical field is 0.00080 T. The nearest-neighbor dipolar term appears much larger, increasing \(H_c\) by 0.08246 T. That conclusion does not survive the long-range sum: M3 returns to within −0.00196 T of M0 under tin-foil conditions and lies +0.05172 T above M0 for a spherical-vacuum boundary. Thus a nearest-neighbor dipolar truncation can assign both the wrong magnitude and, relative to other weak terms, the wrong design priority. The macroscopic dipolar convention is itself a reportable source of model risk.

The P2 interlayer network shifts M4 relative to M3 by 0.09527 T under either boundary convention. This invariance occurs because the same dipolar surface term is present on both sides of the ablation. The simple P1 vertical bond gives only 0.03157 T. The comparison shows that the topology and anisotropy of the frustrated interlayer network, not the existence of an arbitrary small interlayer energy scale, drive the threshold-crossing result.

![Bulk model hierarchy and decision-scale bridge](../figures/bulk_bridge_engineering_screen.png)

**Figure 1.** (a) Critical fields across M0–M4 for the tin-foil Ewald convention. M2 demonstrates why nearest-neighbor dipoles are a poor proxy for the long-range interaction. (b) The public 90 mK reference profile transported by the calculated critical-field shifts; the vertical scale is inherited and is not a direct M4 heat calculation. (c) All four representative-M4 transport and boundary scenarios exceed the predeclared 0.05 T threshold.

### C. Public diffraction converts P2 into a constrained model set

The representative P2 calculation is reproduced independently: \(\delta_c=0.02483\), folded \(q_z=0.33429\), an out-of-plane one-magnon bandwidth of 4.384 microelectronvolt, and an M4–M3 critical-field ablation of 0.09527 T. The source reports 4.4 microelectronvolt for the same bandwidth.

Only seven of the 375 deterministic stress-grid points satisfy the public near-saturation propagation-vector constraints. Their \(\delta_c\) values span 0.02445–0.02525, folded \(q_z\) spans 0.32544–0.33429, and out-of-plane bandwidth spans 3.738–5.280 microelectronvolt. Most parameter combinations are therefore rejected by an observable that was not used in the refrigerator calculation.

The retained models produce M4–M3 critical-field ablations of 0.06832–0.16365 T. With the central transported heat shape, every retained model crosses both engineering gates at 90 mK: \(\Delta B^*=0.05844\)–0.16365 T and regret \(=5.13\)–20.44%. The published representative set is no longer the logical basis of the result; it lies inside a public-data-constrained set.

![Public propagation-vector constrained interlayer ensemble](../figures/interlayer_parameter_ensemble.png)

**Figure 2.** (a) Soft modes of all 375 grid points. Colored points satisfy the public in-plane and out-of-plane propagation-vector bands; gray points fail at least one condition. (b) The 90 mK decision effects of the seven retained models under two measured critical-field anchors and two transport mappings. Every central-shape case lies above both engineering gates. (c) The retained critical-field ablations; the dashed line marks the published representative set. Grid fractions are descriptive coverage, not probabilities.

### D. Finite-temperature shape stress identifies the guarantee boundary

The rigid-profile result does not survive every allowed deformation. Across 2,268 combinations of retained interlayer shifts, measured critical-field anchors, transport mappings, residual peak positions, asymmetric widths, and heat scales, 88.1% cross both engineering gates. All central-shape cases pass, but the full envelope reaches \(\Delta B^*=0.03844\) T and 2.24% regret in the least separated case.

The failures are localized rather than arbitrary. Both gates are guaranteed over all 324 shape cases for each retained model whose M4–M3 critical-field ablation is at least 0.10850 T. The two smallest retained ablations, 0.06832 and 0.07371 T, pass only 55.6% and 66.7% of their shape cases. The published representative ablation of 0.09527 T passes 94.4%. The dominant adverse combination is a \(-0.02\) T residual peak displacement with a broadened low-field side. This gives a concrete validation target: resolving the differential peak displacement or establishing an ablation above 0.1085 T closes the decision.

We next reduced the symmetric peak-location uncertainty while retaining both critical-field anchors, both transport mappings, and every asymmetric-width combination. A common 95% operating field first appears when the residual peak uncertainty is no larger than ±0.002 T. Under the original envelope, the possible 90 mK peak spans 1.3764–1.4742 T. A line-shape-free measurement plan is therefore a 26-point sweep from 1.376 to 1.476 T in 0.004 T increments; the nearest sampled field is then at most 0.002 T from the peak. This is a field-grid requirement, not a statistical power calculation, and calibration and heat-measurement errors must also fit within the target.

For the actual material, measured \(H_c\) already includes every real interaction. Minimizing worst regret across the 324 measured-\(H_c\) shape scenarios at each temperature gives Table II.

**Table II. Shape-robust actual-material decisions.** Heat and volume are magnetic-only worst cases at the minimax field.

| Temperature (mK) | Robust field (T) | Worst regret (%) | Guaranteed \(W_{95}\) (T) | Worst \(Q\) (mJ cm\(^{-3}\)) | Volume for 0.10 microWatt, 60 min (cm\(^3\)) |
|---:|---:|---:|---|---:|---:|
| 90 | 1.4052 | 9.21 | none | 0.3754 | 0.959 |
| 100 | 1.4030 | 7.98 | none | 0.4868 | 0.740 |
| 120 | 1.3970 | 6.17 | none | 0.7388 | 0.487 |
| 150 | 1.3915 | 4.31 | 1.3836–1.3942 | 1.1688 | 0.308 |

At 90–120 mK no single field remains within 95% of every deformed optimum. This is a result of the deliberately broad shape envelope, not a claim that operation is impossible. The minimax field remains tightly located near 1.40 T, while a guaranteed high-efficiency tolerance requires denser near-saturation thermodynamic data.

![Thermodynamic shape stress and minimax decision](../figures/thermodynamic_shape_stress.png)

**Figure 3.** (a) The 90 mK heat-profile envelope and minimax field. (b) P2–P0 field displacement and regret across the shape stress; points below either dashed gate identify the unresolved corner. (c) Actual-material minimax field and worst regret from 90 to 150 mK.

### E. The result translates into field and inventory specifications

The lattice-derived density is 4.206 g cm\(^{-3}\). Under the full shape envelope, a 0.10 microWatt total cold load held for 60 min requires 4.03, 3.11, 2.05, and 1.30 g of ideal magnetic material at 90, 100, 120, and 150 mK, respectively. The public 1.65 g sample corresponds to 0.392 cm\(^3\); at 0.10 microWatt its conditional magnetic-only holds are 24.5, 31.8, 48.3, and 76.4 min at those temperatures.

These numbers specify magnetic inventory rather than finished-device performance. Packing voids, binder, a thermal bus, addenda, heat switches, wiring, contact resistance, leakage, and irreversible ramp loss must be added to the total load and volume. The design chart is therefore directly scalable by a measured or specified total cold load without presenting generic parasitic-load assumptions as NBCP data.

![Refrigerator field and material design chart](../figures/refrigerator_design_chart.png)

**Figure 4.** (a) Minimax field and guaranteed \(W_{95}\) tolerance; crosses mark temperatures with no common 95% interval under the full envelope. (b) Required NBCP mass for a 60 min magnetic-only hold. (c) Conditional hold time of a volume equivalent to the public 1.65 g sample.

### F. Direct public data give a temperature-dependent field schedule

Table III uses only published central-value entropy curves and compares the measured fields without field interpolation.

**Table III. Direct 90–300 mK candidate-field decisions.** \(W_{95}\) is a candidate set.

| Temperature (mK) | Best candidate (T) | Candidates in \(W_{95}\) (T) | \(Q_{\max}\) (mJ cm−3) |
|---:|---:|---|---:|
| 90 | 1.5 | 1.5 | 0.5386 |
| 100 | 1.5 | 1.5 | 0.7223 |
| 120 | 1.5 | 1.5, 1.3 | 1.0315 |
| 150 | 1.3 | 1.3, 1.5 | 1.5856 |
| 200 | 1.3 | 1.3, 0 | 2.6007 |
| 250 | 0 | 0, 1.3 | 3.9135 |
| 300 | 0 | 0 | 5.4379 |

The best measured candidate changes near 128 mK from 1.5 to 1.3 T and near 215 mK from 1.3 to 0 T when the published temperature curves are linearly interpolated within each field. At 100 mK and 1.5 T, the magnetic entropy budget is 0.7223 mJ cm−3. A 0.10 microWatt load for 60 min therefore requires 0.498 cm3 of magnetic material under the stated idealized accounting; 0.4 cm3 provides 48.2 min. These values are transparent conditional budgets rather than a device forecast.

### G. Low-energy spectral flow passes the ideal quantum gate

The selected spectral-flow estimator passes all nine N=9 field-temperature cases at 1.40, 1.45, and 1.50 T and 90, 100, and 120 mK. N=12 also passes all nine cases with 96 low-energy levels and certified tail bounds. At N=9, 1.45 T, and 100 mK, a direct finite-difference entropy estimate requires approximately \(1.55\times10^8\) shots, whereas the ideal spectral-flow budget is approximately \(3.83\times10^4\) shots for the ordinary accuracy target. At the strict 1.45 T, 90 mK sign-change point, the conservative budget is \(1.42\times10^6\) shots for a \(2\times10^{-4}\;k_B\) per-spin target.

The N=9 physics-informed state preparation results are listed in Table IV. All final entropy errors remain below their predeclared targets and all signs are correct, even though the hardest 1.40 T doublet has a lower full-state fidelity. This shows why observable bias must be checked directly rather than using fidelity as the only circuit criterion.

**Table IV. N=9 state-preparation and observable tests.**

| Field (T) | Low levels | Distinct circuits | Minimum fidelity | Maximum \(|\Delta S|\) error (\(k_B\)/spin) | Maximum error/target |
|---:|---:|---:|---:|---:|---:|
| 1.40 | 7 | 12 | 0.97636 | 0.000615 | 0.615 |
| 1.45 | 5 | 9 | 0.99935 | 0.000151 | 0.350 |
| 1.50 | 5 | 9 | 0.99806 | 0.000567 | 0.567 |

![N=9 state preparation across magnetic field](../figures/spectral_flow_stateprep_N9_field_scan.png)

**Figure 5.** (a) Exact and circuit-prepared M3–M0 entropy corrections. (b) Absolute total error divided by the predeclared target; every value remains below one. (c) Infidelity of the hardest state at each field. The lowest fidelity does not correspond to failure of the target observable.

### H. IBM Fez proves the execution path, not quantitative accuracy

The live preflight selected IBM Fez. Across the N=6 circuits, the maximum transpiled two-qubit count was 84 CZ, the minimum product of instruction success probabilities was 0.805, and the maximum circuit-duration-to-median-T1 ratio was 0.040. These diagnostics supported a minimal N=6 pilot but not N=9 execution.

For the hardware point \(B=1.05\) T and \(\lambda=0.5\), the exact \(\langle V\rangle\) is 0.04873698 K. TREX returned 0.020955±0.001848 K, with an absolute error of 0.027782 K. TREX+ZNE returned 0.031144±0.015893 K, with an absolute error of 0.017593 K. The ZNE central bias is 36.7% smaller, but its uncertainty is much larger; the improvement is not statistically established. The two jobs consumed 19 and 52 Open Plan QPU seconds, respectively.

The job identifiers, mitigation comparison, and hardware figure are reported in the Supplementary Information because the QPU values do not enter the refrigerator decision.

## V. Discussion

### A. What is established for magnetic-refrigerator design

Three conclusions are robust within the available public evidence.

First, an interaction that is weak in energy can matter because the refrigerator operates near a narrow critical-field feature. Public propagation-vector data retain seven distinct interlayer models, and every one crosses both decision thresholds under the central heat-profile transport after the known mean-field bias is removed by experimental anchoring.

Second, a local dipolar truncation is not a safe engineering surrogate. M2 suggests a large upward shift, while M3 shows that the rest of the lattice can cancel it. Long-range electrostatic or magnetic interactions must be summed with an explicit boundary convention before ranking weak terms.

Third, the appropriate deliverable is a robust field, its regret, and the boundary at which the conclusion stops being guaranteed. The full shape envelope gives a 90 mK minimax field of 1.405 T but no common \(W_{95}\). An M4–M3 ablation of at least 0.1085 T guarantees both engineering gates over the declared envelope. Below that value, a near-saturation heat-profile measurement constraining the residual peak position to ±0.002 T is more valuable than another small-cluster calculation. The proposed 1.376–1.476 T sweep converts that statement into a directly testable measurement specification.

### B. What is established for quantum computing

The quantum contribution is problem-specific. The desired output is a small correction to a few thermally occupied levels, not a full thermal density matrix or a complete wavefunction. This lets spectral flow replace a noisy entropy derivative with expectation values of the weak-term operator. In ideal simulation, that replacement produces a decisive sampling advantage over the direct finite-difference estimator for the tested instances.

The result is not a quantum-advantage demonstration. Exact classical eigenvectors supplied the N=9 state-preparation targets, and the N=12 test used ideal low-energy spectral data. The IBM experiment used only N=6 and did not feed into the material conclusion. The hardware result identifies the next bottleneck quantitatively: N=9 circuits reach up to 374 ECR gates on a representative heavy-hex backend, while N=6 already shows a large gate-noise bias. Symmetry tapering, shallower problem-informed preparation, and observable-specific mitigation are required before another QPU run is justified.

### C. Limitations

The largest physical limitation is the absence of a converged three-dimensional finite-temperature M4 solution. The expanded shape envelope tests peak displacement, asymmetric width, and heat-scale changes but cannot span every entropy redistribution allowed by M4. We therefore report both the central-profile model-set result and the full-envelope guarantee boundary rather than attributing a unique heat curve to any retained Hamiltonian.

The public experimental curves provide only four terminal fields. The direct 90–300 mK \(W_{95}\) result is a measured-candidate set. The absence of a common continuous \(W_{95}\) at 90–120 mK under the shape envelope is a value-of-information result, not an experimental confidence statement. Experimental confidence intervals are unavailable for several source curves, so we do not fabricate statistical error bars or interpret stress-grid coverage as probability.

The heat budgets include only magnetic entropy. A complete refrigerator model must add nuclear and phonon heat capacities, addenda, heat switches, wiring, contact resistance, magnetic-ramp losses, and thermal leakage. Conditional volume calculations are lower-bound design scenarios.

Finally, Ewald boundary conventions bracket a macroscopic shape effect but do not replace a demagnetization calculation for a specified sample geometry. A future device calculation should couple the microscopic susceptibility to its actual shape.

## VI. Conclusions

This study converts a quantum-magnet Hamiltonian question into a falsifiable refrigerator decision. In NBCP, \(J_{\pm\pm}\) alone has little effect on the polarized critical field. A nearest-neighbor dipolar model appears important but is qualitatively misleading because the long-range sum cancels much of its shift. More importantly, the interlayer conclusion no longer rests on one representative Hamiltonian: public propagation-vector data retain seven of 375 stress-grid models, with critical-field ablations of 0.0683–0.1636 T.

Under the central heat profile, all seven retained models pass the predeclared field-shift and regret gates. Under a broader finite-temperature shape envelope, 88.1% of cases pass both, and every shape case passes when the critical-field ablation is at least 0.1085 T. This boundary defines the next measurement: at 90 mK, a 1.376–1.476 T sweep in 0.004 T increments can localize the peak to ±0.002 T without assuming a fitted line shape, which restores a common 95% operating field within the declared envelope.

Public measurements support a temperature-dependent candidate schedule of 1.5 T at 90–120 mK, 1.3 T at 150–200 mK, and 0 T at 250–300 mK. The shape-robust 90 mK minimax field is 1.405 T with 9.21% worst regret, while the lack of a common \(W_{95}\) prevents a guaranteed narrow field-tolerance claim. A 0.10 microWatt, 60 min magnetic-only hold requires 0.959 cm\(^3\), or 4.03 g, at the 90 mK worst-case heat budget.

The same weak-term correction has a quantum-computing structure: only a few low-energy derivatives are needed. Spectral-flow simulations and N=9 state-preparation tests pass their predeclared accuracy targets, but N=6 IBM hardware retains too much bias for quantitative thermodynamics. Present hardware therefore supports a reproducible execution proof, while the refrigerator conclusion remains classical and public-data based.

The broader design lesson is that microscopic realism should be judged by the decision it changes. For frustrated refrigerants, long-range summation and interlayer topology can matter more than the apparent magnitude of a local weak interaction, and a quantum algorithm becomes relevant only after that decision sensitivity has been demonstrated classically.

## Data and code availability

All experimental inputs are public source-data files or public article archives. Raw public files are preserved unchanged with SHA-256 hashes. Derived CSV and JSON files, source code, unit tests, figures, stopping decisions, IBM job identifiers, and a complete file manifest are included in the associated repository. IBM credentials and cloud resource names are excluded. A versioned GitHub release and Zenodo DOI will be added after author review.

## Author contributions

Ugwiyeon Lee: Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Visualization, Writing—original draft, Writing—review and editing, Project administration.
Seong-Yun Bae: Writing—review and editing.
Hakrae Kim: Writing—review and editing.
Jungsu Lee: Writing—review and editing.

The review-only roles are provisional until each coauthor confirms the manuscript and contribution statement.

## Competing interests

The authors declare no competing interests. This statement requires confirmation by all authors before submission.

## Acknowledgments

The authors acknowledge the providers of the public experimental source data and the IBM Quantum Open Plan used for the two hardware pilot jobs. No statement here implies endorsement by those providers.

## References

1. J. Xiang *et al.*, “Giant magnetocaloric effect in spin supersolid candidate Na2BaCo(PO4)2,” *Nature* **625**, 270–275 (2024), doi:10.1038/s41586-023-06885-w.
2. Y. Gao *et al.*, “Spin supersolidity in nearly ideal easy-axis triangular quantum antiferromagnet Na2BaCo(PO4)2,” *npj Quantum Materials* **7** (2022), doi:10.1038/s41535-022-00500-3.
3. J. Sheng *et al.*, “Two-dimensional quantum universality in the spin-1/2 triangular-lattice quantum antiferromagnet Na2BaCo(PO4)2,” *PNAS* **119**, e2211193119 (2022), doi:10.1073/pnas.2211193119.
4. D. Flavián *et al.*, “Field-tuned incommensurate fan phase in an Ising-like triangular antiferromagnet,” arXiv:2607.23173v1 (2026), doi:10.48550/arXiv.2607.23173.
5. J. Jaklič and P. Prelovšek, “Lanczos method for the calculation of finite-temperature quantities in correlated systems,” *Physical Review B* **49**, 5065–5068 (1994), doi:10.1103/PhysRevB.49.5065.
6. S. Sugiura and A. Shimizu, “Thermal Pure Quantum States at Finite Temperature,” *Physical Review Letters* **108**, 240401 (2012), doi:10.1103/PhysRevLett.108.240401.
7. S. R. White, “Minimally Entangled Typical Quantum States at Finite Temperature,” *Physical Review Letters* **102**, 190601 (2009), doi:10.1103/PhysRevLett.102.190601.
8. S.-N. Sun *et al.*, “Quantum Computation of Finite-Temperature Static and Dynamical Properties of Spin Systems Using Quantum Imaginary Time Evolution,” *PRX Quantum* **2**, 010317 (2021), doi:10.1103/PRXQuantum.2.010317.
9. G. Gentinetta *et al.*, “Quantum finite-temperature Lanczos method,” *Physical Review Research* **8**, 023285 (2026), doi:10.1103/f3q5-cm1n.
10. IBM Quantum, “Plans overview,” https://quantum.cloud.ibm.com/docs/en/guides/plans-overview (accessed 12 October 2026).
11. IBM Quantum, “Error mitigation and suppression techniques,” https://quantum.cloud.ibm.com/docs/en/guides/error-mitigation-and-suppression-techniques (accessed 12 October 2026).
12. P. J. Shirron, “Applications of the magnetocaloric effect in single-stage, multi-stage and continuous adiabatic demagnetization refrigerators,” *Cryogenics* **62**, 130–139 (2014), doi:10.1016/j.cryogenics.2014.03.014.
13. P. J. Shirron, “Optimization strategies for single-stage, multi-stage and continuous ADRs,” *Cryogenics* **62**, 140–149 (2014), doi:10.1016/j.cryogenics.2014.03.012.
14. P. Wikus, E. Canavan, S. T. Heine, K. Matsumoto, and T. Numazawa, “Magnetocaloric materials and the optimization of cooling power density,” *Cryogenics* **62**, 150–162 (2014), doi:10.1016/j.cryogenics.2014.04.005.
15. P. Wikus, G. Burghart, and E. Figueroa-Feliciano, “Optimum operating regimes of common paramagnetic refrigerants,” *Cryogenics* **51**, 555–558 (2011), doi:10.1016/j.cryogenics.2011.07.001.
16. T. I. Popescu *et al.*, “Zeeman Split Kramers Doublets in Spin-Supersolid Candidate Na2BaCo(PO4)2,” *Physical Review Letters* **134**, 136703 (2025), doi:10.1103/PhysRevLett.134.136703.
17. L. Woodland, R. Okuma, J. R. Stewart, C. Balz, and R. Coldea, “From continuum excitations to sharp magnons via transverse magnetic field in the spin-1/2 Ising-like triangular lattice antiferromagnet Na2BaCo(PO4)2,” *Physical Review B* **112**, 104413 (2025), doi:10.1103/1pvl-kzjm.
18. S. W. de Leeuw, J. W. Perram, and E. R. Smith, “Simulation of electrostatic systems in periodic boundary conditions. I. Lattice sums and dielectric constants,” *Proceedings of the Royal Society of London A* **373**, 27–56 (1980), doi:10.1098/rspa.1980.0135.
