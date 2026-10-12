"""Integrate frozen adaptive N=6 live-QPU evidence into the Cryogenics draft."""
from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANUSCRIPT = ROOT / "manuscript" / "Manuscript_Cryogenics_20261012.md"


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new, 1)


def main():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    if "adaptive six-qubit hardware benchmark" in text:
        print("QPU manuscript integration already present")
        return

    text = replace_once(
        text,
        "# Finite-size and package-shape gates for a sub-kelvin Na2BaCo(PO4)2 magnetic-refrigerator stage",
        "# Finite-size, quantum-hardware accuracy, and package-shape gates for a sub-kelvin Na2BaCo(PO4)2 magnetic-refrigerator stage",
        "title")

    old_abstract = (
        "Selecting a sub-kelvin magnetic refrigerant requires proving that a microscopic heat profile transfers across numerical size and macroscopic sample shape. Using public data for Na2BaCo(PO4)2, we audit both gates. An adaptive search retains 69 of 823 exact interlayer models with a 0.0334–0.2021 T critical-field-shift range. We then recalculate a preregistered minimum, nominal, median-distinct, and maximum-effect M4 cover, plus the M3 baseline, on a three-sublattice-compatible 18-spin cluster. The four N=18 optimum shifts span 0.293–0.425 T and their peak-heat ratios span 0.537–0.852; 4 of four models cross the declared major size-transfer threshold. The common grid exposes multiple close peaks in 5 of five Hamiltonians; close-peak, low-state, and projected-magnon-cutoff checks meet the declared local tolerances inside the N=18 projection, while N=18 is not treated as bulk. Propagating all four correlated N=18 templates expands the design set to 68 scenarios and gives 34.55% unmeasured worst regret at 1.4890 T. Exhaustive search over 5,671 two-point panels selects 1.292 and 1.716 T, leaving 4.246% regret at a ±2.5% heat bound. A separate demagnetization envelope gives a 0.136 T plate correction and requires c/a at least 2.72 for a 0.02 T bound. The defensible outputs are therefore a falsifiable heat-measurement plan and a package-geometry specification.")
    new_abstract = (
        "Selecting a sub-kelvin magnetic refrigerant requires testing whether a microscopic heat profile transfers across numerical size, whether a quantum processor resolves the required correction, and whether internal field maps to package field. Using public data for Na2BaCo(PO4)2, we audit these gates. An adaptive search retains 69 of 823 exact interlayer models with a 0.0334–0.2021 T critical-field-shift range. A preregistered four-model cover on a three-sublattice-compatible 18-spin cluster gives optimum shifts of 0.293–0.425 T and peak-heat ratios of 0.537–0.852; all four models cross the declared major size-transfer threshold, while N=18 is not treated as bulk. An adaptive six-qubit hardware benchmark then reconstructs a complete three-level spectral-flow entropy correction from nine circuits on IBM Kingston. Twirled readout mitigation and zero-noise extrapolation both recover the exact negative sign and exclude zero, but the best hardware value, −0.02682 kB per spin, remains biased relative to the exact −0.04161; hardware magnitude therefore fails the quantitative gate. Propagating the correlated N=18 stresses gives 34.55% unmeasured worst regret at 1.4890 T. Exhaustive search over 5,671 two-point panels selects 1.292 and 1.716 T, leaving 4.246% regret at a ±2.5% heat bound. A demagnetization envelope gives a 0.136 T plate correction and requires c/a at least 2.72 for a 0.02 T bound. The defensible outputs are a falsifiable heat-measurement plan, a measured quantum-hardware accuracy gate, and a package-geometry specification.")
    text = replace_once(text, old_abstract, new_abstract, "abstract")

    old_questions = (
        "We therefore ask three engineering questions. First, can a realistic weak interaction shift the optimal terminal field by at least 0.05 T or cause at least 5% loss when a simpler-model field is used? Second, does that conclusion survive finite-size diagnostics, long-range dipolar boundary conventions, and non-unique interlayer parameters constrained by public diffraction? Third, what robust field, material inventory, loss allowance, minimum heat-measurement panel, and measurement precision follow from the public evidence between 90 and 300 mK?")
    new_questions = (
        "We therefore ask four engineering questions. First, can a realistic weak interaction shift the optimal terminal field by at least 0.05 T or cause at least 5% loss when a simpler-model field is used? Second, does that conclusion survive finite-size diagnostics, long-range dipolar boundary conventions, and non-unique interlayer parameters constrained by public diffraction? Third, can a finite-temperature correction selected by those decision gates be reconstructed on present quantum hardware with the correct sign and usable magnitude? Fourth, what robust field, material inventory, loss allowance, minimum heat-measurement panel, and measurement precision follow from the public evidence between 90 and 300 mK?")
    text = replace_once(text, old_questions, new_questions, "questions")

    old_flow = (
        "The calculation deliberately separates what is measured, what is modelled, and what remains a device input. Public thermodynamic data supply the absolute magnetic entropy. Bulk stability calculations supply interaction-induced field shifts. Local-grid, expanded Sobol, and adaptive boundary/extremum searches test whether public propagation vectors identify a unique engineering correction. A lower-field fan-curve holdout then tests those retained models against field dependence not used for retention. Finally, steady-load and transient-loss envelopes translate the magnetic budget into stage mass, duty, and measurements needed for cycle closure.")
    new_flow = (
        "The calculation deliberately separates what is measured, what is modelled, what is executed on quantum hardware, and what remains a device input. Public thermodynamic data supply the absolute magnetic entropy. Bulk stability calculations supply interaction-induced field shifts. Local-grid, expanded Sobol, and adaptive boundary/extremum searches test whether public propagation vectors identify a unique engineering correction. A lower-field fan-curve holdout then tests those retained models against field dependence not used for retention. A decision-targeted spectral-flow benchmark tests a complete low-energy correction on six physical qubits without assigning that small cluster a bulk role. Finally, steady-load and transient-loss envelopes translate the magnetic budget into stage mass, duty, and measurements needed for cycle closure.")
    text = replace_once(text, old_flow, new_flow, "workflow")

    old_central_end = (
        "This rejects universal quantitative transfer and reclassifies N=12-derived field and inventory values as pre-convergence examples. An independent demagnetization envelope converts applied/internal-field ambiguity into a package aspect-ratio requirement.")
    new_central_end = (
        "This rejects universal quantitative transfer and reclassifies N=12-derived field and inventory values as pre-convergence examples. On IBM Kingston, paired nine-circuit jobs recover the negative sign of a six-spin spectral-flow entropy correction, while zero-noise extrapolation reduces absolute bias by 25.8% but does not pass the quantitative-magnitude gate. An independent demagnetization envelope converts applied/internal-field ambiguity into a package aspect-ratio requirement.")
    text = replace_once(text, old_central_end, new_central_end, "central result")

    methods_anchor = (
        "To connect the internal microscopic field to the applied package field, we use the ellipsoidal relation \\(B_{\\rm int}=B_{\\rm app}-N_z\\mu_0M\\). The public \\(B\\parallel c\\) value \\(g_c=4.79\\), one effective spin-1/2 per unit cell, and the public lattice constants give \\(\\mu_0M_s=0.16357\\) T. An equal-volume spheroid is used as a reproducible shape envelope. Saturation magnetization gives a conservative upper bound; an exact package correction requires its dimensions, orientation, packing fraction, and measured \\(M(B,T)\\). The public 1.87 by 1.16 by 0.18 mm heat-capacity plate provides an external scale for this approximation [20].")
    qpu_methods = methods_anchor + (
        "\n\n### E. Decision-targeted spectral-flow quantum-hardware benchmark\n\n"
        "The finite-size audit first establishes that weak interactions change a refrigerator decision before quantum hardware is invoked. We then benchmark one complete low-energy correction on an open \\(3\\times2\\), \\(N=6\\) triangular patch. This patch is a circuit benchmark rather than a material or bulk approximation. The path \\(H(\\lambda)=H_{\\rm M0}+\\lambda V\\), with \\(V=H_{\\rm M2}-H_{\\rm M0}\\), converts the finite-temperature problem into low-energy expectation values through the Hellmann–Feynman relation \\(dE_k/d\\lambda=\\langle k,\\lambda|V|k,\\lambda\\rangle\\). Three-point Gauss–Legendre integration reconstructs the three lowest levels; those levels meet the predeclared \\(10^{-3} k_B\\) per-spin entropy-error target at 100 mK. This spectral-flow route avoids direct entropy measurement and follows finite-temperature quantum-algorithm proposals that reconstruct thermal information from controlled state or subspace calculations [21,22].\n\n"
        "Each state starts from one local layer of \\(R_z\\), number-conserving Givens, and pair-mixing gates on the patch graph. Second-layer generators are appended in descending fidelity-gradient order until fidelity exceeds 0.9999 and the absolute \\(V\\)-expectation error is below 0.001 K. Exact eigenvectors are fitting targets, so this procedure tests circuit representation and hardware execution rather than scalable eigenstate discovery. Across nine node/level circuits, the minimum local fidelity is 0.999916, the maximum local \\(V\\) error is 0.000900 K, and the reconstructed entropy error is \\(-9.53\\times10^{-5} k_B\\) per spin.\n\n"
        "The frozen circuits were compiled to IBM Kingston with at most 96 CZ gates. Two Open Plan Estimator jobs used a 0.02 K target precision, XpXm dynamical decoupling, gate twirling, and twirled readout-error extinction; the second additionally used zero-noise extrapolation [23,24]. Each job contains all nine circuits. Estimator-reported standard errors are propagated through the nonlinear entropy reconstruction using 100,000 independent normal draws with a fixed seed. The resulting interval represents reported sampling uncertainty only. Residual hardware bias is evaluated separately against the known exact six-spin result.")
    text = replace_once(text, methods_anchor, qpu_methods, "QPU methods")

    results_anchor = (
        "**Figure 10.** (a) N=12 scenario envelope and 16 correlated stresses derived from four preregistered N=18 M4 templates, with the baseline and augmented minimax fields. (b) The best 120 results from the exhaustive 5,671-panel search; the star marks 1.292/1.716 T. The N=18 curves are adversarial finite-cluster stresses, not bulk predictions.")
    qpu_results = results_anchor + (
        "\n\n### I. The complete six-qubit hardware benchmark resolves sign but not magnitude\n\n"
        "The exact full-spectrum M2–M0 entropy correction for the open six-spin patch is −0.041611 \\(k_B\\) per spin at 100 mK. The adaptive prepared circuits give −0.041707 \\(k_B\\) per spin, confirming that local representation error is below the declared \\(10^{-3} k_B\\) per-spin target before hardware noise.\n\n"
        "The complete twirled-TREX job on IBM Kingston gives −0.021679 \\(k_B\\) per spin, with a reported-error 95% interval of [−0.024627, −0.018727]. Adding zero-noise extrapolation moves the result toward the exact value, to −0.026822 \\(k_B\\) per spin with interval [−0.040300, −0.013193]. Both intervals exclude zero and both central values recover the exact negative sign. Zero-noise extrapolation reduces absolute bias from 0.019932 to 0.014789 \\(k_B\\) per spin, or 25.8%, while increasing the propagated standard deviation by a factor of 4.60. The exact value remains outside both reported-error intervals. The hardware therefore passes the sign-resolution gate and fails the quantitative-magnitude gate. The paired jobs used 38 and 113 s of free QPU time; their identifiers are `db6et9clf4us73c582f0` and `db6eurqmb58s738bmisg`.\n\n"
        "![Complete adaptive N=6 QPU benchmark](../figures/adaptive_n6_qpu_complete.png)\n\n"
        "**Figure 11.** (a) Nine IBM Kingston expectation values against exact values. (b) The three spectral-flow nodes and three low levels; pointwise error bars show two Estimator-reported standard errors. (c) Exact, locally prepared, twirled-TREX, and twirled-TREX-plus-ZNE entropy corrections; hardware bars show propagated 95% reported-error intervals. Shared calibration and residual device bias are not contained in those intervals, so accuracy is judged against the exact six-spin benchmark.")
    text = replace_once(text, results_anchor, qpu_results, "QPU results")

    text = replace_once(
        text,
        "Three conclusions are robust within the available public evidence.",
        "Four conclusions are robust within the declared evidence and hardware tests.",
        "discussion count")

    discussion_anchor = (
        "The multimodel value-of-information calculation turns that failure into a revised experiment. The 68-scenario peak support is 1.2943–1.7136 T, and the unmeasured minimax decision has 34.551% worst regret. At a ±2.5% heat bound, exhaustive enumeration of 5,671 pairs selects 1.292 and 1.716 T and leaves 4.246% regret. The full 107-point, 1.292–1.716 T sweep remains the model-falsification stage. Measurement precision and finite-size transfer are therefore coupled design variables.")
    discussion_new = discussion_anchor + (
        "\n\nFourth, the quantum-hardware result has a decision-scale interpretation rather than a generic execution claim. All nine low-energy expectation values are measured and propagated to the final entropy correction. Present hardware resolves the sign, which distinguishes a correction from zero, but its residual bias is about twenty times the local state-preparation entropy error and prevents quantitative use in the refrigerator field or inventory calculation. The paired mitigation result identifies circuit and device error, rather than sampling alone, as the next quantum-computing bottleneck.")
    text = replace_once(text, discussion_anchor, discussion_new, "QPU discussion")

    table_anchor = (
        "| Does the stage heat balance predict an independent physical device? | two-load YbNi1.6Sn calibration and disjoint temperature holdout [19] | recovers 15.22 g and 0.254 microWatt; holdout NRMSE 0.592% | validates the lumped identity, not NBCP-specific package losses |")
    table_new = table_anchor + (
        "\n| Can present quantum hardware recover the selected low-energy correction? | nine adaptive six-qubit circuits on IBM Kingston, paired twirled-TREX and ZNE jobs | sign resolved and separated from zero; best magnitude −0.02682 versus exact −0.04161 \\(k_B\\) per spin | do not use hardware magnitude for bulk field or inventory; reduce circuit/device bias and demonstrate scalable state discovery |")
    text = replace_once(text, table_anchor, table_new, "QPU scope row")

    limitations_anchor = (
        "The heat budgets include only magnetic entropy. The stage and transient envelopes validate the scaling identities and expose the required device inputs, but the failed size-transfer gate removes the N=12 heat value as a final inventory basis. The YbNi1.6Sn holdout validates the lumped heat-capacity-to-hold-time identity; it does not infer NBCP losses. A complete prediction still requires a bulk-converged or measured NBCP heat profile together with nuclear and phonon heat capacities, package addenda, heat-switch conduction, wiring, contact resistance, ramp-rate-dependent loss, and warm-side regeneration heat.")
    limitations_new = (
        "The QPU benchmark uses six qubits, an open patch, and exact classical eigenvectors as circuit-fitting targets. It therefore demonstrates a complete hardware measurement-and-propagation path, not quantum advantage, scalable eigenstate discovery, or an N=18 calculation. Reported Estimator errors also omit correlated calibration drift and residual device bias; the exact benchmark shows that this omitted bias remains decisive even after ZNE. No QPU value enters the bulk scenario envelope, robust field, or material inventory.\n\n" + limitations_anchor)
    text = replace_once(text, limitations_anchor, limitations_new, "QPU limitation")

    conclusion_anchor = (
        "An independent public-device audit closes the next link. Two YbNi1.6Sn warm-up curves calibrate to 15.22 g effective mass and 0.254 microWatt parasitic load, consistent with separate reports of approximately 15 g and less than 0.4 microWatt. The frozen balance predicts the disjoint 0.25–0.28 K interval with 0.592% normalized error. This supports the stage energy identity on an actual device without transferring its package losses to NBCP.")
    conclusion_new = conclusion_anchor + (
        "\n\nThe quantum-hardware audit closes a different link. Two complete nine-circuit IBM Kingston jobs recover the negative sign of the six-spin spectral-flow entropy correction, and ZNE reduces absolute bias by 25.8%. The exact magnitude remains outside the reported-error interval, so current hardware is retained as a sign-resolving method benchmark and is not used as a quantitative refrigerator input.")
    text = replace_once(text, conclusion_anchor, conclusion_new, "QPU conclusion")

    text = replace_once(
        text,
        "The broader design lesson is that microscopic realism should be judged by the decision it changes. For frustrated refrigerants, long-range summation and interlayer topology can matter more than the apparent magnitude of a local weak interaction, and a quantum algorithm becomes relevant only after that decision sensitivity has been demonstrated classically.",
        "The broader design lesson is that microscopic realism and quantum computation should be judged by the decision they change. For frustrated refrigerants, long-range summation and interlayer topology can matter more than the apparent magnitude of a local weak interaction. Classical sensitivity analysis identifies the observable that matters; a complete hardware benchmark then determines whether present QPU accuracy is sufficient to use it.",
        "final lesson")

    data_anchor = (
        "All experimental inputs are public source-data files or public article archives. Raw public files are preserved unchanged with SHA-256 hashes. Derived CSV and JSON files, source code, unit tests, figures, stopping decisions, and a complete file manifest are included in the associated repository. A versioned GitHub release and Zenodo DOI will be added after author review.")
    data_new = (
        "All experimental inputs are public source-data files or public article archives. Raw public files are preserved unchanged with SHA-256 hashes. Derived CSV and JSON files, source code, unit tests, figures, stopping decisions, and a complete file manifest are included in the associated repository. Sanitized QPU records include circuit resources, mitigation settings, usage accounting, the IBM job identifiers, and returned values; account tokens and service identifiers are excluded. A versioned GitHub release and Zenodo DOI will be added after author review.")
    text = replace_once(text, data_anchor, data_new, "data availability")

    reference_anchor = (
        "20. H. Li *et al.*, “Possible itinerant excitations and quantum spin state transitions in the effective spin-1/2 triangular-lattice antiferromagnet Na2BaCo(PO4)2,” *Nature Communications* **11**, 3996 (2020), doi:10.1038/s41467-020-18041-3.")
    reference_new = reference_anchor + (
        "\n21. S.-N. Sun, M. Motta, R. N. Tazhigulov, A. T. K. Tan, G. K.-L. Chan, and A. J. Minnich, “Quantum computation of finite-temperature static and dynamical properties of spin systems using quantum imaginary time evolution,” *PRX Quantum* **2**, 010317 (2021), doi:10.1103/PRXQuantum.2.010317.\n"
        "22. G. Gentinetta, F. Metz, W. Kirby, and G. Carleo, “Quantum finite-temperature Lanczos method,” *Physical Review Research* **8**, 023285 (2026), doi:10.1103/f3q5-cm1n.\n"
        "23. IBM Quantum, “Plans overview,” https://quantum.cloud.ibm.com/docs/en/guides/plans-overview (accessed 12 October 2026).\n"
        "24. IBM Quantum, “Error mitigation and suppression techniques,” https://quantum.cloud.ibm.com/docs/en/guides/error-mitigation-and-suppression-techniques (accessed 12 October 2026).")
    text = replace_once(text, reference_anchor, reference_new, "references")

    MANUSCRIPT.write_text(text, encoding="utf-8")
    print(f"updated {MANUSCRIPT}")


if __name__ == "__main__":
    main()
