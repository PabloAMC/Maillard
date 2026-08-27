"""Guard tests pinning the HONEST validation headline numbers to reproducible evidence.

ADDED 2026-08-27 (Wave J2 of the audit-remediation campaign).

WHY THIS FILE EXISTS
--------------------
The validation badge in this repo has drifted back toward flattery before. Nine remediation
waves established a set of deliberately unflattering headline numbers -- 0/6, 0/5, 0/14 --
and every one of them lives, today, only in English prose in README.md and AUDIT.md. Prose
is not enforced. Nothing stops a future change from quietly improving a number, and nothing
stops the docs and the model from silently diverging.

Each test below RECOMPUTES a published number from evidence and fails if it moves. Where a
number is a literal it carries a dated comment naming the artifact or code path it came from.

TWO RULES FOR MAINTAINERS
-------------------------
1. If a test here fails because the model got BETTER, that is good news -- but you must
   re-pin the constant AND update the prose in the same change, with a dated comment saying
   what moved and why. Never widen a bound to make this file green.
2. Never relax a pin to accommodate a regression. A red test here means the published
   headline is now false; the fix is the model or the prose, not the assertion.

WHERE THE NUMBERS COME FROM
---------------------------
Most are recomputed LIVE from tracked benchmark files rather than read from
``results/validation/*.json``. That is deliberate: ``.gitignore`` excludes
``results/validation/*`` wholesale, so ``benchmark_summary.json`` and
``validation_overview.json`` -- the artifacts behind "14 benchmarks", "0/6" and "0/14" --
DO NOT EXIST in a fresh clone. A guard that read them would skip itself into uselessness on
CI, which is the self-excusing-skip pattern this wave exists to remove. The full panel
re-evaluates in ~1.3 s, so live recomputation is cheap and strictly stronger.
``prediction_uncertainty.json`` and ``external_validation_report.json`` ARE force-added and
tracked, so those two are read from disk.
"""

import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (  # noqa: E402
    evaluate_benchmark,
    evaluate_benchmark_payload,
    get_benchmark_files,
    load_benchmark,
    summarize_evaluation,
)

README = ROOT / "README.md"
AUDIT = ROOT / "AUDIT.md"
PREDICTION_UNCERTAINTY = ROOT / "results" / "validation" / "prediction_uncertainty.json"
EXTERNAL_VALIDATION = ROOT / "results" / "validation" / "external_validation_report.json"


# --------------------------------------------------------------------------------------
# Shared panel evaluation. Module-scoped so the 14 benchmarks are run once, not once per
# test. ~1.3 s total on the documented conda environment (measured 2026-08-27).
# --------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def panel():
    summaries = []
    for bench_file in get_benchmark_files():
        evaluation = evaluate_benchmark(bench_file)
        bench = load_benchmark(bench_file)
        summaries.append(
            summarize_evaluation(evaluation, protein_type=bench.get("protein_type", "free"))
        )
    return summaries


def _doc_text(path: Path) -> str:
    assert path.exists(), f"{path.name} is missing; the headline claims have no home"
    return path.read_text(encoding="utf-8")


def _assert_quoted(text: str, token: str, where: str, what: str) -> None:
    assert token in text, (
        f"{where} no longer quotes {token!r} for {what}. Either the number changed and this "
        f"guard was not updated with it, or the prose was edited away from the evidence. "
        f"Docs and evidence must not diverge silently."
    )


# --------------------------------------------------------------------------------------
# 1. Panel size and strict-readiness
# --------------------------------------------------------------------------------------


def test_calibration_panel_is_14_benchmarks_and_none_is_strict_ready(panel):
    """PANEL 14 · STRICT-READY 0/14.

    Pinned 2026-08-27 from a live evaluation of the tracked benchmark files, cross-checked
    against results/validation/validation_overview.json (benchmark_count 14,
    strict_ready_count 0) at the time of writing. The panel was 16 before Wave I quarantined
    two fabricated benchmarks (spi_hvp_xylose_120C_PMC9905368,
    wheat_gluten_hvp_xylose_120C_PMC9905368) -- if this count returns to 16, check that the
    quarantine has not been reverted before re-pinning.
    """
    assert len(panel) == 14, (
        f"Calibration panel is {len(panel)} benchmarks, not the published 14. "
        f"Adding or removing a benchmark changes every downstream headline; re-pin this "
        f"guard and the README/AUDIT tables in the same change."
    )

    strict_ready = [s.benchmark_id for s in panel if s.strict_ready]
    assert strict_ready == [], (
        f"The published claim is 0/14 strict-ready. These are now strict-ready: "
        f"{strict_ready}. If that is real, say so in README.md -- the 'no high tier' "
        f"statement in 'When to trust the predictions' depends on this being zero."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "**14**", "README.md", "the panel size")
    _assert_quoted(readme, "0/14", "README.md", "the strict-ready count")


# --------------------------------------------------------------------------------------
# 2. The 0/6 predictive headline -- the single most important number in the repo
# --------------------------------------------------------------------------------------


def test_zero_of_six_predictive_benchmarks_are_free_of_blocking_gaps(panel):
    """PREDICTIVE 0/6 · fit-recovery 1/4 · internal-synthetic 4/4.

    This is the number the whole audit turns on. The panel scores 5/14 "pass" overall, and
    a 7/16 version of that number was once published as evidence of accuracy. Split by what
    each benchmark can actually support (``BenchmarkSummary.evidence_role``), EVERY passing
    benchmark sits in a non-evidence bucket: 1 is a fit-recovery scoring constants that were
    back-solved from it, and 4 are internal synthetic reproducibility rows. Of the 6
    benchmarks that could constitute predictive evidence, ZERO pass.

    RE-PINNED 2026-08-27 (Wave M): fit-recovery passes 3 -> 1 and the aggregate 7 -> 5.
    CAUSE: the Wave K/M content correction of the two Pratap-Singh benchmarks (see their
    `content_correction_note`). Molecules 2021, 26, 4104 Table 1 (PMC8271896) gives hexanal
    = 1138.00 ppb (pea) and 1621.71 ppb (soy), not the 260 / 380 this repo carried, and
    reports hexanol as n.d. rather than 80 / 120 ppb. Their observability factors were
    back-solved from the erroneous values, so with the paper's real numbers the two lanes
    stop recovering and miss by 4.37x and 4.27x -- exactly the size of the transcription
    error. This is the strongest possible demonstration of why fit-recovery is not evidence:
    a recovery that no longer recovers was never measuring the model. The factors were
    deliberately NOT refitted (that is an owner decision, and refitting them in the same
    pass as the Wave N chemistry change would make the agreement unattributable).

    Pinned from a live panel evaluation; matches
    results/validation/benchmark_summary.json at the time of writing.

    The bucket totals (6/4/4) are pinned as well as the pass counts. Without them, moving a
    benchmark from `predictive` to `fit_recovery` would improve the headline while changing
    nothing about the model -- which is precisely the laundering route this guard blocks.
    """
    by_role = {}
    for summary in panel:
        role = summary.evidence_role
        by_role.setdefault(role, []).append(summary)

    assert sorted(by_role) == ["fit_recovery", "internal_synthetic", "predictive"], (
        f"Unexpected evidence_role buckets: {sorted(by_role)}"
    )

    totals = {role: len(rows) for role, rows in by_role.items()}
    assert totals == {"predictive": 6, "fit_recovery": 4, "internal_synthetic": 4}, (
        f"Evidence-role split moved to {totals}, published as 6/4/4. Reclassifying a "
        f"benchmark changes the denominator of the headline claim -- justify it in AUDIT.md."
    )

    def passing(role):
        return sorted(s.benchmark_id for s in by_role[role] if s.overall_status == "pass")

    predictive_passes = passing("predictive")
    assert predictive_passes == [], (
        f"The published headline is 0/6 PREDICTIVE benchmarks without blocking gaps. These "
        f"now pass: {predictive_passes}. Before celebrating, confirm the benchmark is still "
        f"genuinely predictive (its constants were not fitted to it) -- then re-pin here and "
        f"correct README.md and AUDIT.md together."
    )

    assert len(passing("fit_recovery")) == 1, (
        f"fit-recovery passes moved to {len(passing('fit_recovery'))}/4 (published 1/4 "
        f"since the 2026-08-27 Wave K/M Pratap-Singh content correction). The single "
        f"survivor is pea_isolate_uht_140C_Trikusuma2019, whose source has NOT yet been "
        f"content-verified -- if it rises again, check whether a corrected benchmark was "
        f"reverted rather than a model improved."
    )
    assert len(passing("internal_synthetic")) == 4, (
        f"internal-synthetic passes moved to {len(passing('internal_synthetic'))}/4 "
        f"(published 4/4)"
    )

    total_passes = sum(1 for s in panel if s.overall_status == "pass")
    assert total_passes == 5, (
        f"Aggregate panel passes moved to {total_passes}/14 (published 5/14 since the "
        f"2026-08-27 Wave K/M content correction; 7/14 before it, explicitly "
        f"labelled do-not-quote because every pass is in a non-evidence bucket)."
    )

    for doc, path in (("README.md", README), ("AUDIT.md", AUDIT)):
        _assert_quoted(_doc_text(path), "0/6 predictive", doc, "the predictive headline")


# --------------------------------------------------------------------------------------
# 3. Honest external-literature coverage
# --------------------------------------------------------------------------------------


def test_honest_external_literature_coverage_is_1_of_3_with_fitted_rows_excluded():
    """EXTERNAL LITERATURE 1/3 evaluable · 4 not evaluable · 2 fitted rows excluded BOTH sides.

    Source artifact: results/validation/prediction_uncertainty.json (TRACKED / force-added),
    key ``summary.honest_literature_coverage``. Read 2026-08-27.

    The load-bearing detail is ``excluded_fitted_rows_that_would_have_been_hits == 2``: both
    fitted rows WOULD have been counted as literature hits under the old accounting, which
    is how a 3/5 could be reported for a model with one real external hit. They are removed
    from numerator and denominator alike, not merely from the numerator.

    ``median_ci_width_log10`` is pinned too, because 1/3 is only meaningful alongside the
    width of the interval that produced it: 0.856 dex is a factor of ~7.2 from end to end.
    A future change that widened intervals until coverage improved would raise this number
    and get caught here.
    """
    payload = json.loads(PREDICTION_UNCERTAINTY.read_text(encoding="utf-8"))
    summary = payload["summary"]
    coverage = summary["honest_literature_coverage"]

    assert summary["benchmark_count"] == 11, (
        f"MC panel is {summary['benchmark_count']} benchmarks, published as 11"
    )
    assert summary["matched_compound_count"] == 35, (
        f"MC panel matched rows moved to {summary['matched_compound_count']}, published as 35"
    )

    assert (coverage["hits"], coverage["total"]) == (1, 3), (
        f"Honest external-literature coverage moved to "
        f"{coverage['hits']}/{coverage['total']}, published as 1/3."
    )
    assert coverage["not_evaluable"] == 4, (
        f"not_evaluable moved to {coverage['not_evaluable']}, published as 4. Rows that "
        f"cannot be evaluated must stay visible; silently dropping them inflates the rate."
    )
    assert coverage["excluded_fitted_rows"] == 2
    assert coverage["excluded_fitted_rows_that_would_have_been_hits"] == 2, (
        "Both fitted rows would have counted as literature hits under the old accounting. "
        "If this drops to 0, check that fitted rows are still being detected at all."
    )
    assert coverage["median_ci_width_log10"] == pytest.approx(0.8558, abs=5e-4), (
        f"Median CI width moved to {coverage['median_ci_width_log10']:.4f} dex from the "
        f"published 0.856. A coverage rate that improves while this widens is the interval "
        f"getting looser, not the model getting better."
    )


# --------------------------------------------------------------------------------------
# 4. The external hold-out
# --------------------------------------------------------------------------------------


def test_holdout_scores_1_of_5_genuine_extrapolations_at_the_pre_widening_prior():
    """HOLD-OUT 1/5 genuine extrapolations pre-widening · median 15.31x · worst 2474x.

    Source artifact: results/validation/external_validation_report.json (TRACKED /
    force-added). Read 2026-08-27.

    RE-PINNED 2026-08-27 (Wave M): 0/5 -> 1/5 and median 32.79x -> 15.31x. CAUSE: the
    Wave K/M CONTENT correction of the Li 2026 hold-out bundle, NOT a model change -- the
    predictions are byte-identical. Two of its four points had been transcribed from
    adjacent table rows: 2-pentylfuran 221.5 was the paper's *Maltol* row (true 5625.80 ppb)
    and nonanal 29.42 was its *Decanal* row (true 72.66 ppb), both verified against Europe
    PMC fullTextXML (PMC12984281) and both labelled `reported_point_value` while wrong. With
    the paper's own numbers the 2-pentylfuran point goes 49.8x over -> 1.96x over (that is
    the new hit) and nonanal 673x -> 273x. The extreme process-state misses are untouched:
    roasted pea 2474x and HME 1-hexanol 1117x are unchanged, which is why `max_fold_error`
    below is unchanged. Read this as "one reference was wrong", not "the model improved".

    Three things are pinned and each blocks a different way of flattering this number:

    * ``pre_widening_coverage.genuine_extrapolation_hits == 0``. The uncalibrated matrix
      ln-sigma was raised 2.0 -> 2.86 on 2026-08-26. At the SHIPPED sigma the same
      predictions score 2/5. Pinning the pre-widening figure keeps "the interval got wider"
      from being read as "the model got better".
    * The genuine-extrapolation / in-panel-rescoring split (5 and 3). Only bundles at a
      process state absent from the calibration panel test transfer at all; the other 3
      re-run an in-panel anchor at its own conditions. Reclassifying a genuine
      extrapolation as a rescoring would improve the rate without touching the model.
    * ``median_accuracy_fold`` and ``max_fold_error``, which NO prior can change. These are
      the honest error magnitudes and they are the numbers to quote.
    """
    payload = json.loads(EXTERNAL_VALIDATION.read_text(encoding="utf-8"))
    summary = payload["summary"]

    pre = summary["pre_widening_coverage"]
    assert pre["genuine_extrapolation_hits"] == 1, (
        f"Hold-out now covers {pre['genuine_extrapolation_hits']}/"
        f"{pre['genuine_extrapolation_total']} genuine extrapolations at the pre-widening "
        f"prior (published 1/5 since the 2026-08-27 Wave K/M Li 2026 correction; 0/5 before "
        f"it). Verify this is a model change and not a prior change or another reference "
        f"correction before re-pinning: matrix_sigma here must still be "
        f"{pre['matrix_sigma']}."
    )
    assert pre["genuine_extrapolation_total"] == 5
    assert pre["matrix_sigma"] == pytest.approx(2.0), (
        "pre_widening_coverage is no longer scored at the pre-widening sigma of 2.0, so it "
        "is no longer the comparison the README claims it is."
    )

    split = summary["holdout_kind_split"]
    assert split["genuine_extrapolation"]["total"] == 5, (
        f"Genuine extrapolations moved to {split['genuine_extrapolation']['total']} "
        f"(published 5). Reclassifying one as in-panel rescoring improves the headline "
        f"without improving the model."
    )
    assert split["in_panel_rescoring"]["total"] == 3
    assert summary["matched_compound_count"] == 8, (
        f"Hold-out matched points moved to {summary['matched_compound_count']}, published "
        f"as 8 -- of which only 4 are measurements at all."
    )

    assert summary["median_accuracy_fold"] == pytest.approx(15.31, abs=0.01), (
        f"Median hold-out fold error is {summary['median_accuracy_fold']:.2f}x, published "
        f"as 15.31x (32.79x before the 2026-08-27 Wave K/M Li 2026 reference correction)."
    )
    assert summary["max_fold_error"] == pytest.approx(2474.4, abs=0.1), (
        f"Worst hold-out fold error is {summary['max_fold_error']:.1f}x, published as 2474x."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "1/5", "README.md", "the genuine-extrapolation hold-out result")
    _assert_quoted(readme, "2474", "README.md", "the worst hold-out fold error")


# --------------------------------------------------------------------------------------
# 5. The pentose >> hexose ordering claim, and how little of it is structural
# --------------------------------------------------------------------------------------


def test_pentose_hexose_mft_ordering_is_3_39x_not_the_retired_8_98x_or_15_8x():
    """PENTOSE >> HEXOSE 3.39x (ribose 370.3 ppb vs glucose 109.3 ppb at matched conditions).

    RE-PINNED 2026-08-27 (Wave N -- MFT ROUTE CORRECTION). Was 8.98x (ribose 981.3), and
    15.8x before that. CAUSE: the norfuraneol -> MFT step was retired on isotope evidence
    (Cerny & Davidek 2003, 10.1021/jf026123f: authentic norfuraneol spiked into a
    [13C5]ribose/cysteine system does NOT label the MFT; 2004, 10.1021/jf035265m:
    1,4-dideoxypento-2,3-diulose positionally confirmed with [1-13C]ribose). The pentose
    limb now runs Deoxyosone_Reduction (28.0) + Thiol_Addition_Pentodiulose (28.60) instead
    of a single Thiol_Addition_Norfuraneol at the Wave-H-fitted 26.85. Only the pentose side
    moved: glucose is unchanged at 109.3 ppb, which is what identifies the cause.

    This guard exists because the claim is the LAST surviving quantitative boast in the
    README, which makes it exactly the number most likely to drift back up. It
    is pinned to a tolerance of a few percent rather than a floor, in both directions: a
    floor assertion (``>= 3.0``, as in test_pentose_hexose_sulfur_ordering.py) cannot
    detect the ratio silently climbing back toward 8.98x or 15.8x. Note how close 3.39x now
    sits to that 3.0 floor -- one more honest correction of this size would break the
    ordering claim outright, and that would be a finding, not a regression to be tuned away.

    Read it with the README's own caveat, which this test cannot assert cheaply: only
    ~1.13x of the separation is structural. Setting ``thiol_addition_pentodiulose`` (28.60)
    equal to ``thiol_addition_hexose`` (29.65) -- i.e. zeroing the 1.05 kcal/mol gap --
    collapses the ratio to 1.13x, so ~3x of it is carried by an ESTIMATED, UNCONSTRAINED
    barrier rather than by the mechanism. What improved is the KIND of unconstrained: 28.60
    is the un-fitted sulfur-addition class value, where the retired 26.85 was fitted through
    a route the isotope evidence contradicts. The ordering agrees with Hofmann & Schieberle
    1998 (10.1021/jf9705983); the reason the model reproduces it is weaker than the
    agreement.
    """
    from tests.scientific.test_pentose_hexose_sulfur_ordering import (
        _mft_predicted_ppb,
        _payload,
    )

    pentose = evaluate_benchmark_payload(
        _payload("guard_cys_ribose_ordinal", "D-Ribose", ["2-methyl-3-furanthiol"])
    )
    hexose = evaluate_benchmark_payload(
        _payload("guard_cys_glucose_ordinal", "D-Glucose", ["2-methyl-3-furanthiol"])
    )

    ribose_ppb = _mft_predicted_ppb(pentose)
    glucose_ppb = _mft_predicted_ppb(hexose)
    ratio = ribose_ppb / glucose_ppb

    assert ribose_ppb == pytest.approx(370.3, rel=0.01), (
        f"Ribose MFT moved to {ribose_ppb:.1f} ppb from the published 370.3 "
        f"(981.3 before the 2026-08-27 Wave N route correction)"
    )
    assert glucose_ppb == pytest.approx(109.3, rel=0.01), (
        f"Glucose MFT moved to {glucose_ppb:.1f} ppb from the published 109.3. This side "
        f"was NOT touched by Wave N -- the hexose limb still runs the demoted one-step lump "
        f"-- so a move here has a different cause than a move in the ribose value."
    )
    assert ratio == pytest.approx(3.39, rel=0.01), (
        f"Pentose/hexose MFT ratio is {ratio:.2f}x, published as 3.39x. If it went UP, do "
        f"not report it as improved sugar discrimination without first checking how much of "
        f"the change is structural -- ~3x of this ratio rides on `thiol_addition_pentodiulose` "
        f"= 28.60, an ESTIMATED and explicitly UNCONSTRAINED barrier, not on the mechanism."
    )

    for doc, path in (("README.md", README), ("AUDIT.md", AUDIT)):
        _assert_quoted(_doc_text(path), "3.39", doc, "the pentose/hexose ordering margin")


# --------------------------------------------------------------------------------------
# 6. The no_verifiable_source census
# --------------------------------------------------------------------------------------


def _iter_records(node):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from _iter_records(value)
    elif isinstance(node, list):
        for value in node:
            yield from _iter_records(value)


def _contains_number(node) -> bool:
    if isinstance(node, bool):
        return False
    if isinstance(node, (int, float)):
        return True
    if isinstance(node, dict):
        return any(_contains_number(v) for v in node.values())
    if isinstance(node, list):
        return any(_contains_number(v) for v in node)
    return False


def _no_verifiable_source_census():
    """Return {relative_path: (flagged_count, flagged_records_containing_a_number)}.

    Scope and definitions, established by measurement 2026-08-27 and chosen because they
    reproduce README.md's published triple (102 / 80 / 62) exactly:

        flagged  = any object, at any depth, with source_status == "no_verifiable_source"
        numeric  = a flagged record containing at least one numeric value anywhere inside it
        scope    = every parseable .json / .yml / .yaml under data/

    No generator emits this census, so before this test the README's counts were
    unreproducible assertions. They are now recomputed on every run.
    """
    import yaml  # noqa: PLC0415 - only this census needs the YAML surface

    census = {}
    for pattern in ("*.json", "*.yml", "*.yaml"):
        for path in sorted(ROOT.glob(f"data/**/{pattern}")):
            try:
                text = path.read_text(encoding="utf-8")
                payload = json.loads(text) if path.suffix == ".json" else yaml.safe_load(text)
            except Exception:  # unparseable / non-record file; the citation gate owns those
                continue
            flagged = [
                record
                for record in _iter_records(payload)
                if isinstance(record.get("source_status"), str)
                and record["source_status"].strip() == "no_verifiable_source"
            ]
            if flagged:
                census[str(path.relative_to(ROOT))] = (
                    len(flagged),
                    sum(1 for record in flagged if _contains_number(record)),
                )
    return census


def test_no_verifiable_source_census_is_102_records_80_numeric_62_reaching_runtime():
    """no_verifiable_source: 102 flagged · 80 carrying numbers · 62 of those reaching runtime.

    Pinned 2026-08-27 against README.md's "On literature provenance" note. The three
    numbers are pinned separately because they fail for different reasons and a maintainer
    needs to know which one moved.

    The 84 -> 102 step earlier the same day is the reason this guard exists at all. It was
    not a regression: `data/qm/` had been hidden from git by `.gitignore`, so every sweep
    this repo ever ran silently excluded it. Tracking it exposed 18 further records. That is
    exactly the failure mode a prose-only number cannot protect against -- the scope of the
    measurement changed and the published figure did not -- so the SPLIT is pinned below,
    not just the total.

    The runtime figure stays at 62 because the 18 data/qm records are read only by a loader
    test and by the skip-heavy Phase 3 authority lane; none reaches the model. 62 is the
    count that actually matters: an unverifiable citation attached to a number the runtime
    consumes is a fabricated parameter, not merely a bad footnote.

    All three counts are expected to FALL as anchors get verified. When they do, re-pin here
    and in README.md in the same change. They must never rise silently -- a rise means new
    unverifiable numbers entered the registries.
    """
    census = _no_verifiable_source_census()

    total = sum(flagged for flagged, _ in census.values())
    numeric = sum(with_number for _, with_number in census.values())
    runtime = sum(
        with_number
        for path, (_, with_number) in census.items()
        if path.startswith("data/lit/")
    )

    assert total == 102, (
        f"Repo-wide no_verifiable_source count is {total}, published as 102. "
        f"Per file: { {k: v[0] for k, v in census.items()} }"
    )
    assert numeric == 80, (
        f"{numeric} flagged records carry numeric payloads, published as 80."
    )
    assert runtime == 62, (
        f"{runtime} flagged records with numeric payloads sit in data/lit and therefore "
        f"reach the runtime, published as 62. If a data/qm record becomes runtime-reachable "
        f"this number must rise and the README's 'none of those 18 reach the model' claim "
        f"becomes false."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "102 records", "README.md", "the no_verifiable_source census")
    _assert_quoted(readme, "80 carry numeric payloads", "README.md", "the numeric subset")
    _assert_quoted(readme, "62", "README.md", "the runtime-consumed subset")


def test_the_data_qm_records_that_moved_the_census_from_84_to_102_are_still_counted():
    """Pins the 84 (data/lit) + 18 (data/qm) = 102 split, not just the total.

    Separate from the census test above on purpose. A total of 102 can be preserved while
    `data/qm/` silently drops back out of scope and 18 new data/lit records appear -- the
    aggregate would look stable and the honesty gain of 2026-08-27 would be quietly undone.
    This test makes the specific files that caused the 84 -> 102 correction load-bearing.
    """
    census = _no_verifiable_source_census()
    outside_lit = {
        path: flagged
        for path, (flagged, _) in census.items()
        if not path.startswith("data/lit/")
    }

    assert outside_lit == {
        # Measured 2026-08-27: the 18 records exposed when data/qm stopped being gitignored.
        "data/qm/phase33_barrier_benchmarks.json": 9,
        "data/qm/phase35_double_hybrid_benchmarks.json": 9,
    }, (
        f"The no_verifiable_source records OUTSIDE data/lit changed to {outside_lit}. These "
        f"18 are what moved the published census from 84 to 102; if they disappear, check "
        f"whether data/qm was re-hidden from git rather than whether the records were fixed."
    )

    in_lit = sum(
        flagged for path, (flagged, _) in census.items() if path.startswith("data/lit/")
    )
    assert in_lit == 84, (
        f"data/lit carries {in_lit} flagged records, measured as 84 on 2026-08-27 (the "
        f"figure the README published before data/qm was brought into scope)."
    )
