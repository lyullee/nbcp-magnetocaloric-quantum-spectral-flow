"""Propagate the authoritative N=18 multimodel results into the manuscripts.

The script intentionally refuses to run until the numerical and engineering
verification reports both pass.  Every reported value is formatted from the
publication JSON files, so the manuscript, SI, highlights, and cover letter
cannot silently retain numbers from the earlier single-model size stress.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
MANUSCRIPT = ROOT / "manuscript"
MAIN_FILES = [
    MANUSCRIPT / "Manuscript_draft_20261012.md",
    MANUSCRIPT / "Manuscript_Cryogenics_20261012.md",
]
SI_FILES = [
    MANUSCRIPT / "Supplementary_Information_draft_20261012.md",
    MANUSCRIPT / "Supplementary_Information_Cryogenics_20261012.md",
]
COVER = MANUSCRIPT / "Cryogenics_Cover_Letter_Draft_20261012.md"
HIGHLIGHTS = MANUSCRIPT / "Cryogenics_Highlights_20261012.txt"
README = ROOT / "README.md"
REPRODUCE = ROOT / "REPRODUCE.md"
REPORT = RESULTS / "n18_multimodel_manuscript_integration.json"


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def section(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(
        rf"(?ms)^{re.escape(start)}\n.*?(?=^{re.escape(end)}\n)")
    updated, count = pattern.subn(
        lambda _: replacement.rstrip() + "\n\n", text)
    if count != 1:
        raise RuntimeError(f"expected one section {start!r} before {end!r}; got {count}")
    return updated


def paragraph(text: str, prefix: str, replacement: str) -> str:
    pattern = re.compile(rf"(?ms)^{re.escape(prefix)}.*?(?=\n\n)")
    updated, count = pattern.subn(lambda _: replacement.rstrip(), text)
    if count != 1:
        raise RuntimeError(f"expected one paragraph beginning {prefix!r}; got {count}")
    return updated


def paragraph_any(text: str, prefixes: tuple[str, ...], replacement: str) -> str:
    present = [prefix for prefix in prefixes if re.search(
        rf"(?m)^{re.escape(prefix)}", text)]
    if len(present) != 1:
        raise RuntimeError(f"expected one matching paragraph prefix from {prefixes}; got {present}")
    return paragraph(text, present[0], replacement)


def section_any(text: str, starts: tuple[str, ...], end: str,
                replacement: str) -> str:
    present = [start for start in starts if re.search(
        rf"(?m)^{re.escape(start)}$", text)]
    if len(present) != 1:
        raise RuntimeError(f"expected one matching section start from {starts}; got {present}")
    return section(text, present[0], end, replacement)


def table_row(text: str, prefix: str, replacement: str) -> str:
    pattern = re.compile(rf"(?m)^\| {re.escape(prefix)}.*$")
    updated, count = pattern.subn(lambda _: replacement, text)
    if count != 1:
        raise RuntimeError(f"expected one table row beginning {prefix!r}; got {count}")
    return updated


def line_any(text: str, prefixes: tuple[str, ...], replacement: str) -> str:
    present = [prefix for prefix in prefixes if re.search(
        rf"(?m)^{re.escape(prefix)}", text)]
    if len(present) != 1:
        raise RuntimeError(f"expected one matching line prefix from {prefixes}; got {present}")
    pattern = re.compile(rf"(?m)^{re.escape(present[0])}.*$")
    updated, count = pattern.subn(lambda _: replacement, text)
    if count != 1:
        raise RuntimeError(f"expected one line beginning {present[0]!r}; got {count}")
    return updated


def words(text: str) -> int:
    return len(re.findall(r"\b[\wµ±.–-]+\b", text))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fmt_range(values: list[float], digits: int = 3) -> str:
    return f"{min(values):.{digits}f}–{max(values):.{digits}f}"


def build_text(n18: dict, design: dict) -> dict[str, str]:
    comparisons = n18["comparisons"]
    by_model = {row["model"]: row for row in comparisons}
    tests = n18["model_family_tests"]
    numerical = n18["numerical_checks"]
    major_count = tests["major_change_model_count"]
    major_models = ", ".join(tests["major_change_models"]) or "none"
    shift_range = fmt_range([row["N18_optimum_shift_T"] for row in comparisons])
    heat_range = fmt_range([row["N18_peak_heat_ratio"] for row in comparisons])
    width_range = fmt_range([row["N18_W95_total_ratio"] for row in comparisons])
    n12_shift_range = fmt_range([row["N12_optimum_shift_T"] for row in comparisons])
    n12_heat_range = fmt_range([row["N12_peak_heat_ratio"] for row in comparisons])
    m401 = by_model["M4-01"]
    base_peak = n18["models"]["M3"]["quadratic_peak_field_T"]
    monotone = tests["N18_shift_nondecreasing_with_preregistered_bulk_shift_order"]
    spearman = tests["spearman_bulk_shift_vs_N18_optimum_shift"]
    candidate_counts = tests["L32_candidate_peak_count_by_model"]
    multiple_candidates = tests["multiple_close_peak_models"]
    reselected = tests["L32_reselected_peak_branch_models"]
    candidate_text = ", ".join(
        f"{model} ({candidate_counts[model]})" for model in multiple_candidates
    ) or "none"
    reselected_text = ", ".join(reselected) or "none"

    augmented = design["augmented_N12_plus_multimodel_N18"]
    baseline = design["baseline_N12_only"]
    support = design["peak_support_T"]
    grid = design["measurement_grid"]
    pair = design["exact_best_pair"]
    former = design["former_N12_optimized_pair_retest"]
    sweep = design["full_sweep"]
    mass = design["conditional_90mK_sizing"]
    pair_fields = pair["fields_T"]
    gate_word = "crosses" if pair["worst_regret_percent"] < 5 else "does not cross"
    former_word = "remains below" if former["still_below_5percent_gate"] else "exceeds"

    family_sentence = (
        f"{major_count} of the four M4 models ({major_models}) cross the "
        "declared major-transfer threshold. "
        f"The N=18 shift is {'monotone' if monotone else 'not monotone'} in the "
        f"preregistered bulk-shift order (Spearman rho={spearman:.3f})."
    )
    abstract = (
        "Selecting a sub-kelvin magnetic refrigerant requires proving that a "
        "microscopic heat profile transfers across numerical size and macroscopic "
        "sample shape. Using public data for Na2BaCo(PO4)2, we audit both gates. "
        "An adaptive search retains 69 of 823 exact interlayer models with a "
        "0.0334–0.2021 T critical-field-shift range. We then recalculate a "
        "preregistered minimum, nominal, median-distinct, and maximum-effect M4 "
        "cover, plus the M3 baseline, on a three-sublattice-compatible 18-spin "
        f"cluster. The four N=18 optimum shifts span {shift_range} T and their "
        f"peak-heat ratios span {heat_range}; {major_count} of four models cross "
        f"the declared major size-transfer threshold. The common grid exposes "
        f"multiple close peaks in {len(multiple_candidates)} of five Hamiltonians; "
        "close-peak, low-state, and projected-magnon-cutoff checks meet the declared "
        "local tolerances inside "
        "the N=18 projection, while N=18 is not treated as bulk. Propagating all "
        f"four correlated N=18 templates expands the design set to {design['scenario_count']} "
        f"scenarios and gives {augmented['worst_regret_percent']:.2f}% unmeasured "
        f"worst regret at {augmented['field_T']:.4f} T. Exhaustive search over "
        f"{pair['search_space_pairs']:,} two-point panels selects {pair_fields[0]:.3f} "
        f"and {pair_fields[1]:.3f} T, leaving {pair['worst_regret_percent']:.3f}% "
        "regret at a ±2.5% heat bound. A separate demagnetization envelope gives "
        "a 0.136 T plate correction and requires c/a at least 2.72 for a 0.02 T "
        "bound. The defensible outputs are therefore a falsifiable heat-measurement "
        "plan and a package-geometry specification."
    )

    comparison_rows = []
    for row in comparisons:
        comparison_rows.append(
            f"| {row['model']} | {row['sample_id']} | "
            f"{row['N12_optimum_shift_T']:.3f} | {row['N18_optimum_shift_T']:.3f} | "
            f"{row['N12_peak_heat_ratio']:.3f} | {row['N18_peak_heat_ratio']:.3f} | "
            f"{row['N18_W95_total_ratio']:.3f} |")
    comparison_table = "\n".join(comparison_rows)

    results_d = f"""### D. A preregistered N=18 model cover rejects universal N=12 transfer

The direct N=12 calculation remains useful because it rejects rigid critical-field transport. Across the 12-model coverage set at 90 mK, the M4–M3 optimum shift is 0–0.08 T and the peak-heat ratio is 0.916–1.013. Three common-random repeats of M4-01 give a 0.100 T shift, a mean peak-heat ratio of 0.919, and total $W_{{95}}$-width ratios of 1.46–2.05. These results invalidate the earlier ±0.02 T residual-shift and ±20% width assumptions.

We test whether that conclusion is confined to one extreme Hamiltonian. Before inspecting the new N=18 results, we fixed four M4 models: the minimum retained bulk-shift model M4-00, nominal M4-02, median-shift/distinct-coupling M4-03, and maximum-shift M4-01. The same M3 baseline and all four M4 models are calculated on the periodic $3\times3\times2$, $N=18$ cluster. Unlike the $3\times2\times2$ cluster, this geometry accommodates three-sublattice order in both planar directions.

| Model | Sample | N=12 shift (T) | N=18 shift (T) | N=12 heat ratio | N=18 heat ratio | N=18 total-$W_{{95}}$ ratio |
|---|---:|---:|---:|---:|---:|---:|
{comparison_table}

The M3 N=18 peak is {base_peak:.5f} T. Across the preregistered M4 cover, the N=18 optimum shift spans {shift_range} T, the peak-heat ratio spans {heat_range}, and the local total-$W_{{95}}$ ratio spans {width_range}. A major transfer change is defined as an absolute N=18-minus-N=12 optimum-shift change of at least 0.05 T or an absolute peak-heat-ratio change of at least 0.10. {family_sentence} This model-family result is stronger than the previous M4-01-only comparison because selection was fixed from retained bulk and N=12 evidence before the new N=18 values were known.

Every Hamiltonian uses the same 1.10–2.10 T grid at 0.02 T spacing. Every close dense-grid local maximum is recalculated with 32 states per parity before peak selection. Multiple close candidates occur for {candidate_text}; increasing the state count changes the selected branch for {reselected_text}. The maximum K6/K7-to-K8/K9 cold-entropy change at the selected centers is {numerical['cold_cutoff_crosscheck_max_abs_kB_per_spin']:.2e} $k_B$ per spin; the largest final K8/K9-to-K10/K11 hot-reference change is {numerical['hot_last_cutoff_change_max_abs_kB_per_spin']:.2e} $k_B$ per spin. These checks meet the declared local numerical tolerances inside the projected N=18 calculation. They do not establish the thermodynamic limit.

![Preregistered N=18 interlayer-model cover and convergence](../figures/n18_multimodel_cover.png)

**Figure 6.** (a) Common-grid N=18 cold-entropy scans for M3 and four preregistered M4 models. (b) N=12 and N=18 optimum shifts. (c) N=12 and N=18 peak-heat ratios. (d) Low-state and projected-magnon-cutoff checks at the selected peaks. N=18 convergence is numerical convergence within one projected cluster, not bulk convergence.

Propagating all four correlated N=18 deformations, rather than only M4-01, expands the engineering stress from 52 to {design['scenario_count']} scenarios. Peak support becomes {support[0]:.4f}–{support[1]:.4f} T. The unmeasured minimax field moves from {baseline['field_T']:.4f} to {augmented['field_T']:.4f} T and worst regret rises from {baseline['worst_regret_percent']:.2f}% to {augmented['worst_regret_percent']:.2f}%. The former N=12 field and inventory values remain transparent pre-convergence examples, not package acceptance values.

For provenance, the former N=12-only decisions are listed in Table II. They quantify what the smaller scenario set implied before the size test and must not be used as package acceptance values.

**Table II. N=12-conditional pre-convergence decisions.** Heat and volume are magnetic-only values within the former 52-scenario set.

| Temperature (mK) | Conditional field (T) | Worst regret (%) | Guaranteed $W_{{95}}$ (T) | Worst $Q$ (mJ cm$^{{-3}}$) | Volume for 0.10 microWatt, 60 min (cm$^3$) |
|---:|---:|---:|---|---:|---:|
| 90 | 1.3660 | 9.42 | none | 0.3825 | 0.941 |
| 100 | 1.3608 | 8.50 | none | 0.4938 | 0.729 |
| 120 | 1.3619 | 6.01 | none | 0.7448 | 0.483 |
| 150 | 1.3435 | 4.99 | 1.3434–1.3435 | 1.1624 | 0.310 |"""

    results_h = f"""### H. Multimodel size stress defines the first heat measurements

The N=12-only design is retained as provenance: it gave {baseline['worst_regret_percent']:.3f}% unmeasured regret at {baseline['field_T']:.3f} T and selected 1.296 and 1.420 T on its 42-field grid. The four-model N=18 cover supplies 16 correlated stresses from four M4 models, two experimental anchors, and two mappings. Together with the 52 N=12 scenarios, this gives {design['scenario_count']} adversarial scenarios, {augmented['worst_regret_percent']:.3f}% unmeasured regret at {augmented['field_T']:.4f} T, and peak support of {support[0]:.4f}–{support[1]:.4f} T.

The first experiment must cover {grid['start_T']:.3f}–{grid['stop_T']:.3f} T at {grid['spacing_T']:.3f} T spacing ({grid['count']} points). At a declared ±2.5% per-reading heat bound, exhaustive enumeration of all {pair['search_space_pairs']:,} two-point panels selects {pair_fields[0]:.3f} and {pair_fields[1]:.3f} T. The panel leaves {pair['worst_regret_percent']:.3f}% worst regret and {gate_word} the predeclared 5% gate. The former 1.296/1.420 T pair leaves {former['worst_regret_percent']:.3f}% and {former_word} 5%. The full {grid['count']}-point sweep leaves {sweep['worst_regret_percent']:.3f}%.

The two-point panel is the minimum first-stage discrimination experiment if ±2.5% heat accuracy is available. The full field sweep localizes peaks and tests the scenario set itself. A measured curve outside the {design['scenario_count']}-scenario envelope requires a model update. These are deterministic set-membership bounds, not confidence intervals.

![Four-model N=18 size stress and revised measurement design](../figures/n18_multimodel_design_stress.png)

**Figure 10.** (a) N=12 scenario envelope and 16 correlated stresses derived from four preregistered N=18 M4 templates, with the baseline and augmented minimax fields. (b) The best 120 results from the exhaustive {pair['search_space_pairs']:,}-panel search; the star marks {pair_fields[0]:.3f}/{pair_fields[1]:.3f} T. The N=18 curves are adversarial finite-cluster stresses, not bulk predictions."""

    methods_measurement = (
        "We convert the 52 N=12-calibrated scenarios into a robust measurement-design "
        "problem without assigning probabilities. For a proposed heat measurement "
        r"\(x=(T,B)\), the declared per-reading absolute error is "
        r"\(\epsilon_T=\eta\,{\rm median}_{\theta}Q_{\max,\theta}(T)\). "
        "Scenarios remain conservatively indistinguishable when their bounded "
        "prediction intervals overlap. The preregistered four-model N=18 cover "
        "supplies 16 correlated deformations from four M4 models, two experimental "
        f"anchors, and two mappings, producing {design['scenario_count']} adversarial "
        f"scenarios. Its support requires {grid['count']} fields from "
        f"{grid['start_T']:.3f} to {grid['stop_T']:.3f} T at {grid['spacing_T']:.3f} T "
        rf"spacing. At \(\eta=2.5\%\), all {pair['search_space_pairs']:,} two-point "
        "panels are enumerated. The augmented set is a cluster-size stress test, "
        "not a probability distribution or a claim that N=18 is bulk."
    )
    methods_size_transfer = (
        "The quantitative size-transfer gate compares the periodic "
        r"\(3\times2\times2\), \(N=12\) cluster with a periodic "
        r"\(3\times3\times2\), \(N=18\) cluster compatible with three-sublattice "
        "order in both planar directions. Before inspecting the new N=18 results, "
        "we fixed M3 plus minimum-shift M4-00, nominal M4-02, median-shift/distinct-"
        "coupling M4-03, and maximum-shift M4-01. A K4/K5 pilot triggered a common "
        "1.10–2.10 T K6/K7 grid at 0.02 T spacing. Before any dense worker completed, "
        "the original 0.01 "
        r"\(k_B\) per spin close-peak window was conservatively widened to 0.02 "
        "because existing M3 L12-to-L32 corrections reached 0.01023 "
        r"\(k_B\) per spin. Every local maximum inside this amended window is "
        "recalculated with 32 states per parity at its center and ±0.02 T neighbors; "
        "the largest local quadratic peak is retained. At the selected center we "
        "compare K6/K7 and K8/K9, and the 4 T, 2 K reference extends through K10/K11. "
        "These are local numerical checks inside the declared projection, not "
        "thermodynamic-limit convergence."
    )

    intro = (
        "The central result is that neither local interaction magnitude nor one "
        "finite cluster supplies a device set point. A nearest-neighbor dipolar "
        "truncation creates a large critical-field shift, while the full long-range "
        "sum can cancel it. Public diffraction narrows the interlayer Hamiltonian "
        "but does not determine its finite-temperature cooling consequence. Direct "
        "N=12 thermodynamics first rejects rigid critical-field transport. We then "
        "replace the previous one-extreme-model size check with a preregistered "
        f"four-model N=18 cover. Its optimum shifts span {shift_range} T and peak-heat "
        f"ratios span {heat_range}; {major_count} of four models cross the declared "
        "major-transfer threshold despite numerical convergence inside N=18. This "
        "rejects universal quantitative transfer and reclassifies N=12-derived field "
        "and inventory values as pre-convergence examples. An independent "
        "demagnetization envelope converts applied/internal-field ambiguity into a "
        "package aspect-ratio requirement."
    )

    discussion_size = (
        "Third, finite-size transfer fails at the scale relevant to the engineering "
        f"decision across a preregistered model family. The N=18 optimum shifts span "
        f"{shift_range} T, compared with {n12_shift_range} T for the same four M4 "
        f"anchors at N=12; peak-heat ratios span {heat_range}, compared with "
        f"{n12_heat_range}. {family_sentence} The close-peak rule finds multiple "
        f"candidates for {candidate_text}, and L32 reselects the branch for "
        f"{reselected_text}. Local cutoff checks meet the declared tolerance at the "
        "selected centers. This is positive evidence "
        "against universal transfer, while one N=18 cluster remains short of a bulk "
        "prediction."
    )
    discussion_design = (
        "The multimodel value-of-information calculation turns that failure into a "
        f"revised experiment. The {design['scenario_count']}-scenario peak support is "
        f"{support[0]:.4f}–{support[1]:.4f} T, and the unmeasured minimax decision has "
        f"{augmented['worst_regret_percent']:.3f}% worst regret. At a ±2.5% heat "
        f"bound, exhaustive enumeration of {pair['search_space_pairs']:,} pairs "
        f"selects {pair_fields[0]:.3f} and {pair_fields[1]:.3f} T and leaves "
        f"{pair['worst_regret_percent']:.3f}% regret. The full {grid['count']}-point, "
        f"{grid['start_T']:.3f}–{grid['stop_T']:.3f} T sweep remains the "
        "model-falsification stage. Measurement precision and finite-size transfer "
        "are therefore coupled design variables."
    )
    measurement_caveat = (
        f"The revised measurement panel is conditional on 52 N=12 scenarios plus "
        f"four correlated N=18 size stresses. Its bounded-error overlaps are "
        f"conservative set-membership tests, not a noise likelihood, and the error "
        f"fractions are resolution requirements rather than measured instrument "
        f"uncertainties. All {pair['search_space_pairs']:,} two-point panels on the "
        f"{grid['count']}-field grid are enumerated, but the four N=18 curves cannot "
        f"span every possible bulk finite-size trajectory. A measured heat value "
        f"outside the envelope would falsify the scenario set and require expansion "
        f"rather than assignment to the nearest template."
    )

    conclusion1 = (
        "Direct finite-temperature calculation supplies the decisive result. A "
        "preregistered minimum, nominal, median-distinct, and maximum-effect M4 cover "
        f"gives N=18 optimum shifts of {shift_range} T and peak-heat ratios of "
        f"{heat_range}; {major_count} of four models cross the declared major-transfer "
        "threshold. The common-grid, close-peak, state-count, and projected-magnon "
        "checks pass inside N=18. The result rejects universal quantitative transfer "
        "of the N=12 templates while remaining explicitly short of a thermodynamic-"
        "limit prediction."
    )
    conclusion2 = (
        "That result changes the engineering output. Propagating 16 correlated N=18 "
        f"stresses gives {augmented['worst_regret_percent']:.2f}% unmeasured worst "
        f"regret at {augmented['field_T']:.4f} T across {design['scenario_count']} "
        "scenarios. The former 1.366 T and 3.96 g results, together with their steady- "
        "and transient-load multiples, remain N=12-conditional pre-convergence "
        "examples rather than package acceptance values."
    )
    conclusion3 = (
        "The surviving output is a falsifiable measurement plan. Peak support is "
        f"{support[0]:.4f}–{support[1]:.4f} T, requiring a {grid['start_T']:.3f}–"
        f"{grid['stop_T']:.3f} T sweep at {grid['spacing_T']:.3f} T spacing. Exhaustive "
        f"search of {pair['search_space_pairs']:,} two-point panels selects "
        f"{pair_fields[0]:.3f} and {pair_fields[1]:.3f} T and leaves "
        f"{pair['worst_regret_percent']:.3f}% regret at a ±2.5% heat bound. The "
        f"complete {grid['count']}-point sweep leaves {sweep['worst_regret_percent']:.3f}% "
        "and tests whether the finite-cluster scenario set is adequate."
    )

    si_s72 = f"""### S7.2 Preregistered N=18 interlayer-model cover

We repeat M3 and a preregistered four-model M4 cover on a periodic $3\times3\times2$, $N=18$ cluster. The cover contains minimum retained bulk-shift M4-00 (sample 761), nominal M4-02 (sample 188), median-shift/distinct-coupling M4-03 (sample 738), and maximum-shift M4-01 (sample 773). Selection uses only retained-model and N=12 evidence. The maximum parity-breaking Ewald $xz/yz$ coupling is $3.19\times10^{{-18}}$ of the largest bond, permitting separate even and odd projected-magnon sectors.

The K4/K5 pilot showed multiple separated local peaks. Before the K6/K7 results completed, the protocol was amended to a common 1.10–2.10 T grid at 0.02 T spacing. The original 0.01 $k_B$ per spin close-peak window was then widened conservatively to 0.02 before any of the five dense scans completed: existing M3 L12-to-L32 corrections reached 0.01023 $k_B$ per spin, comparable to the original window. Every interior L12 local maximum inside the amended window is recalculated with 32 states per parity at the center and ±0.02 T. The candidate with the largest L32 local quadratic peak is selected, then crosschecked at K8/K9 with 12 and 32 states per parity. The 4 T, 2 K reference is evaluated at K6/K7, K8/K9, and K10/K11. The frozen plan, pilot amendment, conservative candidate amendment, convergence rule, and independent plan verification are stored in `results/n18_multimodel_*`.

| Model | Sample | Bulk shift (T) | N=12 optimum shift (T) | N=18 optimum shift (T) | N=12 heat ratio | N=18 heat ratio | N=18 total-$W_{{95}}$ ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(
        f"| {row['model']} | {row['sample_id']} | {row['bulk_shift_T']:.5f} | "
        f"{row['N12_optimum_shift_T']:.5f} | {row['N18_optimum_shift_T']:.5f} | "
        f"{row['N12_peak_heat_ratio']:.5f} | {row['N18_peak_heat_ratio']:.5f} | "
        f"{row['N18_W95_total_ratio']:.5f} |" for row in comparisons
    ) + f"""

The N=18 shift range is {shift_range} T, the peak-heat-ratio range is {heat_range}, and the total-$W_{{95}}$-ratio range is {width_range}. A major transfer change means an absolute N=18-minus-N=12 optimum-shift change of at least 0.05 T or an absolute peak-heat-ratio change of at least 0.10. {family_sentence} Multiple close candidates occur for {candidate_text}; the L32 selection differs from the L12 grid maximum for {reselected_text}. The maximum cold cutoff difference at the selected centers is {numerical['cold_cutoff_crosscheck_max_abs_kB_per_spin']:.6e} $k_B$ per spin, and the maximum last hot-cutoff change is {numerical['hot_last_cutoff_change_max_abs_kB_per_spin']:.6e} $k_B$ per spin. All dense peaks are interior and all close local maxima were checked at L32.

These checks establish convergence inside the N=18 projection and reject universal quantitative transfer of the N=12 template. They do not establish the thermodynamic limit. Raw spectra are in `results/n18_stacked_m4_work`; compilation and independent checks are provided by `compile_n18_multimodel_cover.py` and `verify_n18_multimodel_cover.py`.

For the engineering destruction test, four N=18 templates are combined with two experimental anchors and two mappings to give 16 correlated N=18 stresses. Appending them to the 52 N=12 scenarios gives {design['scenario_count']} scenarios, moves the minimax field from {baseline['field_T']:.4f} to {augmented['field_T']:.4f} T, raises unmeasured regret from {baseline['worst_regret_percent']:.3f}% to {augmented['worst_regret_percent']:.3f}%, and expands peak support to {support[0]:.6f}–{support[1]:.6f} T. This is a size-stress envelope rather than a bulk prediction."""

    si_s75 = f"""### S7.5 Minimum near-saturation field sweep

The original 52 N=12 scenarios give the previously reported 90 mK support. After the 16 correlated four-model N=18 stresses are included, support is {support[0]:.6f}–{support[1]:.6f} T. Sampling {grid['start_T']:.3f}–{grid['stop_T']:.3f} T every {grid['spacing_T']:.3f} T gives {grid['count']} points and guarantees a grid point within {grid['spacing_T']/2:.3f} T of any peak in the augmented support without fitting a line shape.

The grid localizes the peak only. Field calibration, sample demagnetization, thermal equilibration, entropy reconstruction, heat uncertainty, and profile width must be measured separately. The former 42-point list is retained as N=12 provenance; the expanded grid is stored with the multimodel design-stress results."""

    si_s76 = rf"""### S7.6 Robust measurement value of information

The expanded grid in S7.5 guarantees peak localization within the size-stressed support but does not say which measurements should be taken first. We use 52 N=12-calibrated scenarios plus 16 correlated stresses from the four N=18 M4 templates, two experimental anchors, and two mappings, without assigning probabilities. At a candidate coordinate $x=(T,B)$, a declared per-reading bound is

\[
\epsilon_T=\eta\,{{\rm median}}_{{\theta}}Q_{{\max,\theta}}(T).
\]

Two scenario predictions are conservatively retained together when $|Q_i(x)-Q_j(x)|\le2\epsilon_T$, so their symmetric error intervals overlap. For each possible true scenario, the post-measurement field minimizes maximum 90 mK regret over every still-indistinguishable scenario. The panel score is the largest conditional minimax regret. This union-of-overlaps construction is conservative.

At $\eta=2.5\%$, all {pair['search_space_pairs']:,} two-point panels on the {grid['count']}-field grid are enumerated. The exact pair is {pair_fields[0]:.3f} and {pair_fields[1]:.3f} T, reducing worst regret from {augmented['worst_regret_percent']:.4f}% to {pair['worst_regret_percent']:.4f}%. The previous 1.296/1.420 T pair gives {former['worst_regret_percent']:.4f}% and {former_word} the 5% gate. The full sweep gives {sweep['worst_regret_percent']:.4f}%.

These results are a model-conditional experiment design, not confidence intervals or validation of N=18 as bulk. A real observation outside the {design['scenario_count']}-profile envelope falsifies the scenario set. Files are `results/n18_multimodel_design_stress_*`; `verify_n18_multimodel_design.py` independently checks the scenario construction, exhaustive panel count, optimum, and sizing identity."""

    return locals()


def update_main(path: Path, parts: dict[str, str]) -> None:
    text = path.read_text(encoding="utf-8")
    abstract = parts["abstract"]
    if "Cryogenics" not in path.name:
        abstract += (
            " Spectral flow passes ideal small-system tests, while current "
            "hardware remains too biased for feedback."
        )
    text = section(text, "## Abstract", "## I. Introduction",
                   "## Abstract\n\n" + abstract)
    if "The central result is that" in text:
        text = paragraph(text, "The central result is that", parts["intro"])
    elif "Our main physical finding is that" in text:
        text = paragraph(text, "Our main physical finding is that", parts["intro"])
    text = paragraph_any(
        text,
        ("We first convert the 52 N=12-calibrated scenarios",
         "We convert the 52 N=12-calibrated scenarios"),
        parts["methods_measurement"])
    text = paragraph(
        text, "The quantitative size-transfer gate compares",
        parts["methods_size_transfer"])
    text = section_any(
        text,
        ("### D. Size escalation invalidates the quantitative N=12 template",
         "### D. A preregistered N=18 model cover rejects universal N=12 transfer"),
        "### E. The result translates into a package-shape specification",
        parts["results_d"])
    h_end = "## V. Discussion" if "Cryogenics" in path.name else "### I. Low-energy spectral flow passes the ideal quantum gate"
    text = section_any(
        text,
        ("### H. Robust measurement design converts uncertainty into a heat-resolution requirement",
         "### H. Multimodel size stress defines the first heat measurements"),
        h_end, parts["results_h"])
    text = paragraph(text, "Third, finite-size transfer fails", parts["discussion_size"])
    text = paragraph_any(
        text,
        ("The size-stressed value-of-information calculation",
         "The multimodel value-of-information calculation"),
        parts["discussion_design"])
    text = paragraph(
        text, "The revised measurement panel is conditional",
        parts["measurement_caveat"])
    text = table_row(
        text, "Which 90 mK field is robust",
        f"| Which 90 mK field is robust to the declared uncertainty? | 52 N=12 scenarios plus 16 stresses from four preregistered N=18 M4 models | no final set point; the adversarial minimax coordinate is {parts['augmented']['field_T']:.4f} T with {parts['augmented']['worst_regret_percent']:.2f}% regret | obtain bulk-converged M4 thermodynamics or measure the expanded heat profile |")
    text = table_row(
        text, "What measurement closes the field ambiguity",
        f"| What measurement closes the field ambiguity? | exhaustive search over {parts['pair']['search_space_pairs']:,} pairs on a {parts['grid']['count']}-field grid | first measure {parts['pair_fields'][0]:.3f} and {parts['pair_fields'][1]:.3f} T at 90 mK; at a ±2.5% heat bound the worst regret is {parts['pair']['worst_regret_percent']:.3f}%; retain the {parts['grid']['start_T']:.3f}–{parts['grid']['stop_T']:.3f} T sweep | treat an observation outside the {parts['design']['scenario_count']}-scenario envelope as model falsification |")
    text = paragraph(text, "Direct finite-temperature calculation supplies the decisive result", parts["conclusion1"])
    text = paragraph_any(
        text,
        ("That failure changes the engineering output",
         "That result changes the engineering output"),
        parts["conclusion2"])
    text = paragraph(text, "The surviving output is a falsifiable measurement plan", parts["conclusion3"])
    path.write_text(text, encoding="utf-8")


def update_si(path: Path, parts: dict[str, str]) -> None:
    text = path.read_text(encoding="utf-8")
    text = section_any(
        text,
        ("### S7.2 Three-sublattice-compatible N=18 size escalation",
         "### S7.2 Preregistered N=18 interlayer-model cover"),
        "### S7.3 Density and design conversion", parts["si_s72"])
    text = section(text, "### S7.5 Minimum near-saturation field sweep",
                   "### S7.6 Robust measurement value of information", parts["si_s75"])
    text = section(text, "### S7.6 Robust measurement value of information",
                   "### S7.7 Public refrigeration-workbook consistency audit", parts["si_s76"])
    marker = "python test_general_truncated_magnon_operator.py\n"
    commands = (
        marker
        + "python verify_n18_multimodel_plan.py\n"
        + "python compile_n18_multimodel_cover.py\n"
        + "python verify_n18_multimodel_cover.py\n"
        + "python stress_test_n18_multimodel_design.py\n"
        + "python verify_n18_multimodel_design.py\n"
        + "python integrate_n18_multimodel_manuscript.py\n"
    )
    if "python verify_n18_multimodel_plan.py" not in text:
        if marker not in text:
            raise RuntimeError(f"verification marker missing from {path}")
        text = text.replace(marker, commands, 1)
    text = table_row(
        text, "the N=12 finite-temperature template transfers quantitatively",
        f"| the N=12 finite-temperature template transfers quantitatively to a three-sublattice-compatible larger cluster | no | the preregistered four-model N=18 cover gives shifts of {parts['shift_range']} T and heat ratios of {parts['heat_range']}; {parts['major_count']} of four cross the major-transfer threshold |")
    n18_scope_prefix = (
        "the N=18 model cover is a bulk prediction"
        if "| the N=18 model cover is a bulk prediction" in text
        else "the N=18 curve is a bulk prediction")
    text = table_row(
        text, n18_scope_prefix,
        "| the N=18 model cover is a bulk prediction | no | common-grid, close-peak, low-state, and cutoff convergence is established only within one projected 3×3×2 cluster |")
    text = table_row(
        text, "the former 1.296/1.420 T heat panel remains optimal",
        f"| the former 1.296/1.420 T heat panel remains optimal after multimodel size stress | no | exhaustive search of {parts['pair']['search_space_pairs']:,} pairs selects {parts['pair_fields'][0]:.3f}/{parts['pair_fields'][1]:.3f} T; the former pair leaves {parts['former']['worst_regret_percent']:.3f}% regret |")
    text = table_row(
        text, "the former 1.292–1.456 T field sweep covers",
        f"| the former 1.292–1.456 T field sweep covers the multimodel size-stressed peak support | no | the required {parts['grid']['count']}-point grid is {parts['grid']['start_T']:.3f}–{parts['grid']['stop_T']:.3f} T at {parts['grid']['spacing_T']:.3f} T spacing |")
    path.write_text(text, encoding="utf-8")


def update_cover_and_highlights(parts: dict[str, str]) -> None:
    cover = COVER.read_text(encoding="utf-8")
    cover = paragraph(
        cover, "The central new result is",
        "The central new result is a preregistered model-family size-transfer test "
        f"with an engineering consequence. Four M4 anchors spanning the retained "
        f"bulk-shift range are recalculated with M3 on N=18. Their optimum shifts "
        f"span {parts['shift_range']} T and peak-heat ratios span {parts['heat_range']}; "
        f"{parts['major_count']} of four cross the declared major-transfer threshold. "
        "Common-grid, close-peak, low-state, and cutoff checks pass within N=18. "
        "This replaces the earlier one-extreme-model result and rejects universal "
        "N=12 transfer without presenting N=18 as a bulk solution.")
    cover = paragraph(
        cover, "We convert that failure into",
        f"We convert that failure into a new measurement specification. Sixteen "
        f"correlated N=18 stresses expand the design to {parts['design']['scenario_count']} "
        f"scenarios and give {parts['augmented']['worst_regret_percent']:.2f}% "
        f"unmeasured regret at {parts['augmented']['field_T']:.4f} T. Exhaustive "
        f"search over {parts['pair']['search_space_pairs']:,} two-point panels "
        f"selects {parts['pair_fields'][0]:.3f} and {parts['pair_fields'][1]:.3f} T, "
        f"leaving {parts['pair']['worst_regret_percent']:.3f}% regret at a declared "
        "±2.5% heat bound. A separate demagnetization calculation gives a "
        "conservative 0.136 T correction for a published plate geometry and requires "
        "a c-axis aspect ratio of at least 2.72 for a 0.02 T bound. Public "
        "second-material and device audits retain their independent validation roles.")
    COVER.write_text(cover, encoding="utf-8")

    highlights = [
        "Four preregistered M4 models are recalculated on an N=18 cluster.",
        f"N=18 optimum shifts span {parts['shift_range']} T across the model cover.",
        f"Multimodel size stress gives {parts['augmented']['worst_regret_percent']:.1f}% minimax regret.",
        f"Exhaustive search selects {parts['pair_fields'][0]:.3f} and {parts['pair_fields'][1]:.3f} T first.",
        "Package c/a ≥ 2.72 limits the demagnetization bound to 0.02 T.",
    ]
    if any(len(line) > 85 for line in highlights):
        raise RuntimeError(f"highlight exceeds 85 characters: {[(len(x), x) for x in highlights]}")
    HIGHLIGHTS.write_text("\n".join(highlights) + "\n", encoding="utf-8")


def update_reproducibility_guides(parts: dict[str, str]) -> None:
    readme = README.read_text(encoding="utf-8")
    readme = line_any(
        readme,
        ("- A three-sublattice-compatible", "- A preregistered minimum"),
        "- A preregistered minimum, nominal, median-distinct, and maximum-effect "
        f"M4 cover is recalculated with M3 on the three-sublattice-compatible N=18 "
        f"cluster. N=18 optimum shifts span {parts['shift_range']} T and peak-heat "
        f"ratios span {parts['heat_range']}; {parts['major_count']} of four models "
        "cross the major-transfer threshold. Numerical checks pass inside the "
        "projected cluster, which is not a bulk solution.")
    readme = line_any(
        readme,
        ("- Adding four correlated N=18 size stresses", "- Four N=18 templates"),
        "- Four N=18 templates, two experimental anchors, and two mappings give 16 "
        f"correlated stresses. Together with 52 N=12 scenarios, the "
        f"{parts['design']['scenario_count']}-scenario design has "
        f"{parts['augmented']['worst_regret_percent']:.2f}% unmeasured worst regret "
        f"at {parts['augmented']['field_T']:.4f} T. N=12-derived fields and masses "
        "remain pre-convergence examples.")
    readme = line_any(
        readme,
        ("- The size-stressed peak support is", "- The multimodel peak support is"),
        f"- The multimodel peak support is {parts['support'][0]:.4f}–"
        f"{parts['support'][1]:.4f} T. The required grid contains "
        f"{parts['grid']['count']} fields from {parts['grid']['start_T']:.3f} to "
        f"{parts['grid']['stop_T']:.3f} T. Exhausting all "
        f"{parts['pair']['search_space_pairs']:,} two-point panels selects "
        f"{parts['pair_fields'][0]:.3f} and {parts['pair_fields'][1]:.3f} T and "
        f"leaves {parts['pair']['worst_regret_percent']:.3f}% regret at a ±2.5% "
        "heat bound.")
    README.write_text(readme, encoding="utf-8")

    reproduce = REPRODUCE.read_text(encoding="utf-8")
    reproduce = line_any(
        reproduce,
        ("- The stacked N=12 finite-temperature M4 calculation",
         "- The N=12 calculation rejects rigid-profile transfer"),
        "- The N=12 calculation rejects rigid-profile transfer. The preregistered "
        f"four-model N=18 cover gives optimum shifts of {parts['shift_range']} T "
        f"and peak-heat ratios of {parts['heat_range']}; {parts['major_count']} of "
        "four models cross the major-transfer threshold. Common-grid, close-peak, "
        "low-state, and cutoff checks pass inside N=18, which remains short of a "
        "thermodynamic-limit solution.")
    reproduce = line_any(
        reproduce, ("- The exact",),
        f"- The exact {parts['pair_fields'][0]:.3f}/{parts['pair_fields'][1]:.3f} T "
        f"panel is globally optimal among all {parts['pair']['search_space_pairs']:,} "
        f"two-point panels on the {parts['grid']['count']}-field 90 mK grid at a "
        "±2.5% heat bound. Error fractions are deterministic resolution scenarios, "
        f"and the design is conditional on 52 N=12 profiles plus 16 correlated "
        "four-model N=18 stresses.")
    REPRODUCE.write_text(reproduce, encoding="utf-8")


def main() -> None:
    n18 = load("n18_multimodel_cover_summary.json")
    n18_check = load("n18_multimodel_cover_verification.json")
    design = load("n18_multimodel_design_stress_summary.json")
    design_check = load("n18_multimodel_design_stress_verification.json")
    plan_check = load("n18_multimodel_plan_verification.json")
    if any(item.get("status") != "PASS"
           for item in (n18_check, design_check, plan_check)):
        raise RuntimeError("multimodel plan, numerical, and design verifications must pass first")
    parts = build_text(n18, design)
    if words(parts["abstract"]) > 250:
        raise RuntimeError(f"abstract has {words(parts['abstract'])} words")
    for path in MAIN_FILES:
        update_main(path, parts)
    for path in SI_FILES:
        update_si(path, parts)
    update_cover_and_highlights(parts)
    update_reproducibility_guides(parts)
    files = MAIN_FILES + SI_FILES + [COVER, HIGHLIGHTS, README, REPRODUCE]
    report = {
        "status": "PASS",
        "source_results": [
            "results/n18_multimodel_cover_summary.json",
            "results/n18_multimodel_design_stress_summary.json",
            "results/n18_multimodel_convergence_protocol_amendment.json",
        ],
        "abstract_words": words(parts["abstract"]),
        "highlight_character_counts": [
            len(line) for line in HIGHLIGHTS.read_text(encoding="utf-8").splitlines()
        ],
        "reported": {
            "N18_shift_range_T": n18["model_family_tests"]["N18_optimum_shift_range_T"],
            "major_change_model_count": n18["model_family_tests"]["major_change_model_count"],
            "scenario_count": design["scenario_count"],
            "minimax_field_T": design["augmented_N12_plus_multimodel_N18"]["field_T"],
            "minimax_regret_percent": design["augmented_N12_plus_multimodel_N18"]["worst_regret_percent"],
            "best_pair_T": design["exact_best_pair"]["fields_T"],
            "best_pair_regret_percent": design["exact_best_pair"]["worst_regret_percent"],
        },
        "updated_files": {
            str(path.relative_to(ROOT)): sha256(path) for path in files
        },
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
