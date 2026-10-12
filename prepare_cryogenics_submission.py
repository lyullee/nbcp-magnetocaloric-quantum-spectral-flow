from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
MASTER = ROOT / "manuscript" / "Manuscript_draft_20261012.md"
MASTER_SI = ROOT / "manuscript" / "Supplementary_Information_draft_20261012.md"
OUTPUT = ROOT / "manuscript" / "Manuscript_Cryogenics_20261012.md"
OUTPUT_SI = ROOT / "manuscript" / "Supplementary_Information_Cryogenics_20261012.md"


TITLE = "Finite-size and package-shape gates for a sub-kelvin Na2BaCo(PO4)2 magnetic-refrigerator stage"

ABSTRACT = r"""Selecting a sub-kelvin magnetic refrigerant requires proving that a microscopic heat profile transfers across numerical size and macroscopic sample shape. Using public data for Na2BaCo(PO4)2, we audit both gates. An adaptive search retains 69 of 823 exact interlayer models with a 0.0334–0.2021 T critical-field-shift range; unused fan and neutron data constrain sensitivity without identifying one model. Direct finite-temperature calculations on a 12-spin stacked cluster give a 0.100 T optimum shift and a mean peak-heat ratio of 0.919 for the strongest model. A three-sublattice-compatible 18-spin escalation changes these to 0.425 T and 0.537. Increasing both the low-state count and projected-magnon cutoff leaves the N=18 peak location stable, so the N=12 quantitative template fails a declared size-transfer gate even though N=18 is not a bulk solution. Adding this correlated size stress moves the 90 mK minimax field from 1.366 to 1.4795 T and raises unmeasured worst regret from 9.42% to 32.32%; exhaustive search over 4,753 panels changes the best two fields to 1.324 and 1.404 T at a ±2.5% heat bound. A separate demagnetization envelope gives a 0.136 T saturation correction for a published plate geometry and requires an equivalent c-axis aspect ratio of at least 2.72 to bound it below 0.02 T. The defensible outputs are therefore measurement and package-geometry specifications; the former field and mass values remain pre-convergence examples. External thermodynamic and device audits give 13.8 mK error and recover 15.22 g with 0.59% holdout error."""

INTRODUCTION = r"""Magnetic refrigeration below 1 K depends on entropy stored in a magnetic material. In an adiabatic-demagnetization refrigerator (ADR), useful material properties must ultimately be converted into stage field, refrigerant mass or volume, hold time, and recycle requirements [7–10]. Frustration and quantum fluctuations can retain unusually large entropy at low temperature, making frustrated magnets promising refrigerants. Na2BaCo(PO4)2 (NBCP) is a strong example. Public measurements report a giant magnetocaloric response and cooling below 100 mK, together with an easy-axis triangular-lattice spin-1/2 description [1–4]. Recent spectroscopy further resolves field-dependent Kramers-doublet excitations and a transverse-field evolution from a continuum to sharp magnons [11,12].

The public experiments already establish a material-screening premise. In same-volume time traces, NBCP remains below 1 K for an observed span of 96.0 min, compared with 32.2 min for ferric ammonium alum and 24.7 min for chromium potassium alum. The source article further reports a volumetric heat absorption of approximately 19 mJ cm(^{-3}) and a cycle efficiency of approximately 23%, compared with 3.0–3.1 mJ cm(^{-3}) and 9–11% for the two hydrate salts [1,17]. Relative to the better hydrate comparator for each metric, these values correspond to factors of 6.13 in volumetric heat absorption, 2.09 in source-defined efficiency, and 2.98 in observed time below 1 K. These results show that NBCP merits refrigerator development, but they do not determine a robust near-critical terminal field, the useful field tolerance, or the inventory required for a specified stage load.

The material question and the engineering question are different. A microscopic study asks which spin Hamiltonian produces the observed phases. A refrigerator designer must choose a terminal field, determine how accurately that field must be reached, size the refrigerant, and allocate capacity to steady and transient losses. A model can reproduce broad thermodynamics while moving a narrow low-temperature optimum enough to change those decisions. This risk is greatest near a critical field, where a small error in an interaction energy can produce a substantial loss at a fixed operating point.

The commonly used nearest-neighbor XXZ Hamiltonian captures much of the NBCP phenomenology [1–4]. Recent neutron diffraction, however, found a field-tuned incommensurate fan phase immediately below saturation and showed that dipolar coupling, bond-dependent exchange, and frustrated interlayer paths are required to describe the ordered-state sequence [4]. These terms are small relative to the main exchange, but their consequences for a refrigerator terminal field have not been separated from their role in explaining magnetic order.

We therefore ask three engineering questions. First, can a realistic weak interaction shift the optimal terminal field by at least 0.05 T or cause at least 5% loss when a simpler-model field is used? Second, does that conclusion survive finite-size diagnostics, long-range dipolar boundary conventions, and non-unique interlayer parameters constrained by public diffraction? Third, what robust field, material inventory, loss allowance, minimum heat-measurement panel, and measurement precision follow from the public evidence between 90 and 300 mK?

The calculation deliberately separates what is measured, what is modelled, and what remains a device input. Public thermodynamic data supply the absolute magnetic entropy. Bulk stability calculations supply interaction-induced field shifts. Local-grid, expanded Sobol, and adaptive boundary/extremum searches test whether public propagation vectors identify a unique engineering correction. A lower-field fan-curve holdout then tests those retained models against field dependence not used for retention. Finally, steady-load and transient-loss envelopes translate the magnetic budget into stage mass, duty, and measurements needed for cycle closure.

The central result is that neither local interaction magnitude nor one finite cluster supplies a device set point. A nearest-neighbor dipolar truncation creates a large critical-field shift, while the full long-range sum can cancel it. Public diffraction narrows the interlayer Hamiltonian but does not determine its finite-temperature cooling consequence. Direct N=12 thermodynamics first rejects rigid critical-field transport. A three-sublattice-compatible N=18 escalation then changes the strongest-model optimum shift from 0.100 to 0.425 T and its peak-heat ratio from 0.919 to 0.537, despite convergence within the projected calculation. This fails quantitative size transfer and reclassifies the N=12-derived field and inventory as pre-convergence examples. An independent demagnetization envelope converts the remaining applied/internal-field ambiguity into a package aspect-ratio requirement."""

DISCUSSION = r"""### A. Implications for sub-kelvin refrigerator design

Three conclusions are robust within the available public evidence. First, a local interaction ranking is not a safe engineering surrogate. The nearest-neighbor dipolar model suggests a large upward critical-field shift, while the full Ewald sum shows that the rest of the lattice can cancel it. Long-range magnetic interactions must therefore be summed with an explicit boundary convention before weak terms are assigned a design priority.

Second, public propagation-vector data constrain the interlayer Hamiltonian strongly but do not guarantee a unique refrigerator correction. The local and Sobol designs initially retain 16 models. Three adaptive rounds find 53 more and expand the critical-field-shift range to 0.0334–0.2021 T. Independent five-start calculations reproduce both extrema. Only 224 of 276 central cases cross both gates. The audit therefore gives a set-valued decision range; it does not identify one Hamiltonian or prove that every allowed heat-profile deformation crosses the engineering thresholds.

The lower-field fan curve provides a stronger unused observable. Across the 69 retained models, four-point holdout error rises with critical-field shift at Spearman rho=0.940. The 23 models no worse than the published representative span only 0.0334–0.1073 T, but the representative model itself misses the rapid experimental evolution. This subset is therefore a descriptive ranking rather than a calibrated confidence set. Raw delta(H), qz(H), and phase-boundary data must be fitted jointly before the full engineering range can be narrowed.

The public transverse-field INS archive supplies a second unused observable with a different physical origin. Thirty of 69 models remain at or below the source's 18.048 microelectronvolt bandwidth reference, and membership is unchanged across all matched one-standard-deviation exchange and \(g\)-factor corners. The sensitivity-envelope upper shift falls by 27% from 0.2021 to 0.1476 T, while its shape-stress gate fraction falls to 55.0%. A more direct audit changes the interpretation: a calibrated proxy reproduces the archived linewidth calculation with 0.000141 meV root-mean-square error, but profiling the unresolved twin mixture leaves all 69 M4 models with chi-square no larger than the source curve's 28.35. The INS data therefore demonstrate non-identifiability rather than a calibrated 39-model rejection. The full 0.0334–0.2021 T engineering range remains authoritative.

Third, a useful actual-material decision is possible without making that microscopic attribution. Combining the measured transition interval with 324 public-profile shape scenarios gives a 90 mK minimax field of 1.405 T. The absence of a common 95% interval at 90–120 mK identifies the missing evidence rather than making the design unusable. Localizing the near-saturation heat-profile peak to ±0.002 T is the most direct closure measurement. A 26-point sweep from 1.376 to 1.476 T at 0.004 T spacing implements that requirement without assuming a line shape.

The thermodynamic-processing step also survives an external material check. For 12 public Gd2B2MoO9 comparison points, heat-capacity-derived adiabatic temperature changes reproduce direct measurements with 13.8 mK mean absolute and 27.6 mK maximum error. This result supports the equilibrium processing used to connect field-resolved thermodynamics to cooling response. It does not validate NBCP's microscopic field-shift surrogate or stage package.

The next translation step survives an independent physical-stage check. Using two public YbNi1.6Sn warm-up curves, a calibration below 0.25 K recovers 15.22 g effective mass and 0.254 microWatt parasitic load, consistent with the separately reported approximate 15 g mass and less-than-0.4 microWatt heat leak. With both quantities frozen, the 0.25–0.28 K holdout has 0.592% normalized root-mean-square error. This validates the lumped heat-capacity-to-hold-time identity on another material and device; it does not assign YbNi1.6Sn package losses to NBCP.

### B. From magnetic capacity to stage requirements

The public material comparison and the stage calculation answer different questions. The approximately volume-matched comparison establishes that NBCP has unusually favorable material-level heat absorption and hold behavior. The robust calculation then determines where that capacity should be used and how model uncertainty affects the result. At 90 mK, the reversible lower bound for a 0.10 microWatt, 60 min hold is 4.03 g. This value is a transparent magnetic-capacity baseline, not a packaged-refrigerator prediction.

The separated envelopes show which device quantities dominate the next decision. A steady parasitic load of 0.20 microWatt triples the required mass to 12.10 g, while the 0.97 microWatt literature coordinate raises it to 43.15 g. At zero specified parasitic load, reserving 10%, 25%, or 50% of gross magnetic capacity for transient losses gives 4.48, 5.38, or 8.07 g. The 82 s measured recycle coordinate and 6.5 min target midpoint imply 97.8% and 90.2% single-stage duty for a 60 min hold. Because those timing and load coordinates come from other architectures and refrigerants, they define sensitivities and acceptance thresholds rather than NBCP performance.

**Table IV. Engineering use envelope.** Each output is paired with the decision it can support and the evidence required for a stronger claim.

| Question | Current evidence | Usable output | Boundary or next action |
|---|---|---|---|
| Is NBCP worth refrigerator development? | public same-volume cooling and time traces [1] | yes; material-level screening is established | no new material-ranking claim is needed |
| Which 90 mK field is robust to the declared uncertainty? | measured critical-field interval and 324 actual-material shape scenarios | 1.405 T minimax planning field | no common 95% interval; do not treat it as a calibrated device set point |
| What measurement closes the field ambiguity? | peak-location value-of-information scan | 1.376–1.476 T at 0.004 T spacing | also control field calibration, equilibration, and heat uncertainty |
| How much NBCP is required? | worst declared magnetic entropy budget, crystallographic density, separated steady loads, and transient envelope | 4.03 g reversible bound; 4.48, 5.38, and 8.07 g at 10%, 25%, and 50% transient allocations | measure stage parasitics, transient energy, packing/contact efficiency, and recycle sequence |
| Is one interlayer Hamiltonian uniquely identified? | 69 of 823 exact models satisfy public-(q); fan error correlates with shift at rho=0.940; an indicative INS bandwidth rule retains 30, but direct linewidth/twin profiling leaves all 69 no worse than the source curve | no; propagate the 0.0334–0.2021 T range | jointly constrain raw delta(H), qz(H), phase boundaries, twin fractions, the resolution kernel, polarized interlayer dispersion, and finite-temperature heat |
| Does thermodynamic processing track direct cooling in another material? | 12 public Gd2B2MoO9 indirect/direct points [18] | 13.8 mK mean absolute and 27.6 mK maximum error | validates processing only; it does not transfer a Hamiltonian or package model |
| Does the stage heat balance predict an independent device? | two-load YbNi1.6Sn calibration and disjoint temperature holdout [19] | recovers 15.22 g and 0.254 microWatt; holdout NRMSE 0.592% | validates the lumped identity only; NBCP package losses remain inputs |

### C. Evidence limits

The largest material-model limitation is the absence of a converged three-dimensional finite-temperature M4 solution. The expanded shape envelope tests peak displacement, asymmetric width, and heat-scale changes but cannot span every entropy redistribution allowed by M4. We therefore report the central-profile model-set result and the smallest all-shape-passing shift among evaluated models, rather than assigning a unique heat curve or a continuous guarantee boundary to a retained Hamiltonian.

The public experimental curves provide only four terminal fields. The direct 90–300 mK 95% result is consequently a measured-candidate set. The absence of a common continuous 95% interval at 90–120 mK under the shape envelope is a value-of-information result, not an experimental confidence statement. Several source curves lack experimental covariance, so deterministic stress-grid coverage is not interpreted as probability.

The interlayer search remains finite. The adaptive audit improves coverage from 631 to 823 exact M4 evaluations and targets the public-(q) boundary and shift extrema, but its last-round minimum still changes by 0.00228 T. It does not establish a posterior distribution or mathematically exhaust the continuous parameter box. The source-data audit also finds that the labels of Extended Data Fig. 6 sheets a/b do not reproduce the article's identifiable 94 mK NBCP point. Those curves are excluded from new cross-material calculations; only the independently labelled sheet-d time traces are used for context.

The fan holdout uses four values digitized from a vector figure because a reduced numerical table is not public. Every retained mean-field curve evolves too slowly at the lowest holdout field, so the relative ranking cannot be converted into posterior weights. Simultaneous V–fan and qz(H) calculations against raw diffraction data are needed before the 23-model relative subset can replace the full envelope.

The INS audit compares M4 bandwidths with the source's simplified \(J_b=3\) microelectronvolt calculation. The source calls this value indicative: its linewidth curve gives \(\chi^2=28.35\) for 11 public points without refitting, and its predicted energy modulation was not resolved. The M4 \(J_2,J_a,J_b\) topology is also more general than the source model. A calibrated linewidth proxy reproduces the source curve, but it profiles only an effective twin-mixture coefficient while fixing the published 0.027 meV resolution. Its all-model non-rejection diagnoses insufficient public information; a calibrated reduction requires the full resolution function, independently constrained twin populations, reduced counts, and a joint multi-branch fit.

The lumped heat-capacity-to-hold-time identity is externally supported by the YbNi1.6Sn device holdout, but the NBCP thermal accounting remains conditional. The stage envelope separates specified steady parasitics and declared transient allowances without inferring NBCP package losses. A complete cycle prediction requires nuclear and phonon heat capacities, package addenda, time-resolved heat-switch conduction, wiring, contact resistance, ramp-rate-dependent loss, packing and thermal-boundary efficiency, and warm-side heat rejection. Likewise, the Ewald boundary bracket does not replace a demagnetization calculation for a specified sample geometry."""

CONCLUSIONS = r"""This study converts microscopic uncertainty in a frustrated quantum magnet into a falsifiable sub-kelvin refrigerator design. For NBCP, the nearest-neighbor dipolar truncation gives a qualitatively misleading interaction ranking because the long-range sum cancels most of its apparent critical-field shift. Starting from local and Sobol designs, an adaptive inverse audit brings the total to 823 exact M4 evaluations and retains 69 public-(q)-compatible models. Their 0.0334–0.2021 T critical-field shifts are wider than either initial design, and five independent starts reproduce both extrema. A four-point lower-field fan holdout adds an unused test: its error rises with critical-field shift at rho=0.940, while no retained mean-field model quantitatively reproduces the rapid measured evolution. Independent 3.5 T INS adds an orthogonal test: an indicative 18.048 microelectronvolt bandwidth rule retains 30 models with an upper shift of 0.1476 T, but direct linewidth/twin profiling leaves all 69 no worse than the source curve. The public observables therefore constrain sensitivity without identifying a unique interlayer correction or supporting a calibrated contraction of the full design range.

The resulting engineering boundary is conditional but actionable. Of 276 central cases, 224 cross both predeclared gates; under the full shape envelope, 17,892 of 22,356 cases do so. Compatible low-shift models reduce the field displacement to 0.00857 T and regret to 0.17% under the most adverse profile. At 90 mK, the actual-material minimax terminal field remains 1.405 T with 9.21% worst regret. A 1.376–1.476 T sweep at 0.004 T spacing would localize the heat-profile peak to ±0.002 T and recover a common 95% operating field within the declared envelope.

An external check on 12 public Gd2B2MoO9 points reproduces direct adiabatic temperature changes with 13.8 mK mean absolute error. It validates the equilibrium thermodynamic-processing step on a second frustrated refrigerant without transferring that material's microscopic model or device behavior to NBCP.

An independent public YbNi1.6Sn device test closes the next link in the chain. A low-temperature calibration recovers 15.22 g effective mass and 0.254 microWatt parasitic load against separate reports of approximately 15 g and less than 0.4 microWatt. The frozen model predicts the disjoint 0.25–0.28 K interval with 0.592% normalized error. The governing stage heat balance is therefore supported by an actual device, while NBCP-specific package losses remain to be measured.

For a 0.10 microWatt, 60 min hold, the 90 mK reversible magnetic bound is 0.959 cm3 or 4.03 g. The mass scales exactly with total steady load: the 0.20 and 0.97 microWatt literature coordinates produce 12.10 and 43.15 g before transient losses. Reserving 10%, 25%, and 50% of gross capacity for transient loss raises the zero-parasitic result to 4.48, 5.38, and 8.07 g. These values define package acceptance criteria. Measured stage parasitics, addenda, switch conduction, ramp loss, recycle sequence, contact efficiency, and warm-side heat rejection are the remaining quantities required to replace the conditional envelope with a device prediction.

The broader design lesson is that microscopic realism should be judged by the refrigerator decision it changes. Public diffraction and thermodynamics can support a robust planning field and inventory before a unique Hamiltonian or a complete device exists, provided that model uncertainty, stage losses, and missing measurements are kept explicit."""


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(re.escape(start) + r".*?(?=" + re.escape(end) + r")", re.S)
    updated, count = pattern.subn(lambda _: start + "\n\n" + replacement.rstrip() + "\n\n", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not replace section {start!r}")
    return updated


def remove_between(text: str, start: str, end: str) -> str:
    pattern = re.compile(re.escape(start) + r".*?(?=" + re.escape(end) + r")", re.S)
    updated, count = pattern.subn("", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not remove section {start!r}")
    return updated


def prepare_main() -> None:
    text = MASTER.read_text(encoding="utf-8")
    lines = text.splitlines()
    lines[0] = f"# {TITLE}"
    lines[7] = "Cryogenics-targeted working draft. Numerical values and scope statements are frozen to the reproducibility manifest dated 12 October 2026."
    text = "\n".join(lines) + "\n"

    text = replace_between(text, "## Abstract", "## I. Introduction", ABSTRACT)
    text = replace_between(text, "## I. Introduction", "## II. Decision problem and model hierarchy", INTRODUCTION)
    text = remove_between(text, "### E. Quantum spectral-flow estimator", "## IV. Results")
    text = remove_between(text, "### I. Low-energy spectral flow passes the ideal quantum gate", "## V. Discussion")
    text = remove_between(text, "### B. What is established for quantum computing", "### C. Limitations")
    text = text.replace("### C. Limitations", "### B. Limitations")
    text = re.sub(
        r"The same weak-term correction has a quantum-computing structure:.*?remains classical and public-data based\.\n\n",
        "", text, flags=re.S)
    text = re.sub(
        r"Figure 12 summarizes the evidence chain.*?\*\*Table V\. Engineering use envelope\.\*\*",
        "**Table V. Engineering use envelope.**", text, flags=re.S)
    text = re.sub(r"\| Can present QPU results modify the refrigerator design\?.*?\n", "", text)
    text = text.replace(
        "A shift below 0.02 T and regret below 2% support stopping the quantum stage.",
        "A shift below 0.02 T together with regret below 2% defines the lower stopping band for this engineering screen.",
    )

    text = text.replace(
        "All experimental inputs are public source-data files or public article archives. Raw public files are preserved unchanged with SHA-256 hashes. Derived CSV and JSON files, source code, unit tests, figures, stopping decisions, IBM job identifiers, and a complete file manifest are included in the associated repository. IBM credentials and cloud resource names are excluded. A versioned GitHub release and Zenodo DOI will be added after author review.",
        "All experimental inputs are public source-data files or public article archives. Raw public files are preserved unchanged with SHA-256 hashes. Derived CSV and JSON files, source code, unit tests, figures, stopping decisions, and a complete file manifest are included in the associated repository. A versioned GitHub release and Zenodo DOI will be added after author review.",
    )
    text = text.replace(
        "The authors acknowledge the providers of the public experimental source data and the IBM Quantum Open Plan used for the two hardware pilot jobs. No statement here implies endorsement by those providers.",
        "The authors acknowledge the providers of the public experimental source data. No statement here implies endorsement by those providers.",
    )

    body, refs = text.split("## References", 1)
    citation_replacements = {
        "[12–15]": "[7–10]",
        "[1,22]": "[1,17]",
        "[16,17]": "[11,12]",
        "[19–21]": "[14–16]",
        "[19,20]": "[14,15]",
        "[18]": "[13]",
        "[19]": "[14]",
        "[20]": "[15]",
        "[21]": "[16]",
        "[23]": "[18]",
        "[24]": "[19]",
        "[25]": "[20]",
    }
    for old, new in citation_replacements.items():
        body = body.replace(old, new)

    ref_lines = [line for line in refs.strip().splitlines() if line.strip()]
    kept = []
    for line in ref_lines:
        match = re.match(r"(\d+)\.\s+(.*)", line)
        if not match:
            continue
        number = int(match.group(1))
        if 7 <= number <= 11:
            continue
        new_number = number if number <= 6 else number - 5
        kept.append(f"{new_number}. {match.group(2)}")
    text = body.rstrip() + "\n\n## References\n\n" + "\n".join(kept) + "\n"
    OUTPUT.write_text(text, encoding="utf-8")


def prepare_si() -> None:
    text = MASTER_SI.read_text(encoding="utf-8")
    text = text.replace(
        "## Finite-size and package-shape gates for sub-kelvin refrigerator design with a frustrated quantum magnet",
        f"## {TITLE}",
    )
    text = re.sub(r"- ordinary quantum entropy target:.*?\n", "", text)
    text = re.sub(r"- sign-change quantum target:.*?\n", "", text)
    text = remove_between(text, "## S8. Quantum-algorithm screening", "## S11. Verification suite")
    text = text.replace("## S11. Verification suite", "## S8. Verification suite")
    text = text.replace("## S12. Scope table for reviewers", "## S9. Scope table for reviewers")
    text = text.replace("python verify_quantum_stateprep_qpu.py\n", "")
    text = text.replace("python verify_second_screen.py\n", "")
    text = text.replace("python test_quantum_local_benchmark.py\n", "")
    text = text.replace("python plot_evidence_to_decision_workflow.py\n", "")
    text = text.replace(
        "All passed on Python 3.14.8 with NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2 under Windows 11. Independent Octave 11.3.0 checks reproduced the M0–M4 polarized-state calculation and the N=9 exact thermodynamic regression. The public QMagen low-bond-dimension N=54 route was stopped because it failed convergence and log-partition-function convexity tests.",
        "All listed checks passed on Python 3.14.8 with NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2 under Windows 11. Independent Octave 11.3.0 calculations reproduced the M0–M4 polarized-state result.",
    )
    text = re.sub(r"\| spectral flow is more sample-efficient.*?\n", "", text)
    text = re.sub(r"\| a quantum advantage has been demonstrated.*?\n", "", text)
    text = re.sub(r"\| IBM hardware gives a quantitatively reliable.*?\n", "", text)
    OUTPUT_SI.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    prepare_main()
    prepare_si()
    multimodel_summary = ROOT / "results" / "n18_multimodel_cover_summary.json"
    multimodel_design = ROOT / "results" / "n18_multimodel_design_stress_summary.json"
    if multimodel_summary.exists() and multimodel_design.exists():
        from integrate_n18_multimodel_manuscript import main as integrate_multimodel
        integrate_multimodel()
    from integrate_adaptive_qpu_manuscript import main as integrate_qpu
    integrate_qpu()
    written = OUTPUT.read_text(encoding="utf-8")
    abstract = written.split("## Abstract", 1)[1].split("## I. Introduction", 1)[0]
    abstract_words = len(re.findall(r"\b[\wµ±.–-]+\b", abstract))
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"wrote {OUTPUT_SI.relative_to(ROOT)}")
    print(f"abstract_words={abstract_words}")
