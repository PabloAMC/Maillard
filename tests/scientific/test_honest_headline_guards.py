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


def test_calibration_panel_is_17_benchmarks_and_none_is_strict_ready(panel):
    """PANEL 17 · STRICT-READY 0/17.

    Pinned 2026-08-27 from a live evaluation of the tracked benchmark files, cross-checked
    against results/validation/validation_overview.json (benchmark_count 14,
    strict_ready_count 0) at the time of writing. The panel was 16 before Wave I quarantined
    two fabricated benchmarks (spi_hvp_xylose_120C_PMC9905368,
    wheat_gluten_hvp_xylose_120C_PMC9905368) -- if this count returns to 16, check that the
    quarantine has not been reverted before re-pinning.

    RE-PINNED 2026-08-28 (Wave W): 14 -> 17, strict-ready 0/14 -> 0/17. CAUSE: three
    absolute, isotope-dilution literature anchors were added from the full text of Hofmann &
    Schieberle 1998 (10.1021/jf9705983), obtained by interlibrary loan --
    ``hofmann1998_{ribose,glucose,fructose}_cysteine_145C_20min_pH5``, the pH-5.0 aqueous
    rows of its Table 1. THIS IS THE WAVE THAT ENDED "THE SULFUR BRANCH HAS ZERO ABSOLUTE
    LITERATURE ANCHORS": the branch now has three, from the very paper whose misattributed
    values Wave S2b traced to a repo-internal derivation and Wave S2c retired.
    STRICT-READY DID NOT MOVE, and that is the point worth reading twice. All three new rows
    FAIL, and badly: ribose 12.27x on FFT, glucose 29.58x, fructose 14.46x, with a panel-side
    mean |log10| of 0.92 dex across the six new comparisons. Gaining a real anchor made the
    panel look worse, not better, because the previous anchor was a number the repository had
    written for itself. A panel that gets worse when it is given real data is a panel that was
    not previously being tested.
    """
    assert len(panel) == 17, (
        f"Calibration panel is {len(panel)} benchmarks, not the published 17. "
        f"Adding or removing a benchmark changes every downstream headline; re-pin this "
        f"guard and the README/AUDIT tables in the same change."
    )

    strict_ready = [s.benchmark_id for s in panel if s.strict_ready]
    assert strict_ready == [], (
        f"The published claim is 0/17 strict-ready. These are now strict-ready: "
        f"{strict_ready}. If that is real, say so in README.md -- the 'no high tier' "
        f"statement in 'When to trust the predictions' depends on this being zero."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "**17**", "README.md", "the panel size")
    _assert_quoted(readme, "0/17", "README.md", "the strict-ready count")


# --------------------------------------------------------------------------------------
# 2. The 0/6 predictive headline -- the single most important number in the repo
# --------------------------------------------------------------------------------------


def test_zero_of_nine_predictive_benchmarks_are_free_of_blocking_gaps(panel):
    """PREDICTIVE 0/9 · fit-recovery 0/4 · internal-synthetic 4/4.

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

    2026-08-27 (Wave O refit to content-corrected anchors, owner-approved) -- THE NUMBERS
    ABOVE ARE UNCHANGED, AND THAT NEEDED CHECKING RATHER THAN ASSUMING. The owner decision
    was taken: the two ambient hexanal factors were refitted against the VERIFIED values
    (one shared scale of 4.317249x; record
    results/validation/matrix_observability_refit_pratap_singh.json), so both Pratap-Singh
    rows stopped missing by 4.37x/4.27x and their ``overall_status`` moved
    ``scale-gap`` -> ``pass-no-ranking``. They are NOT counted as passes here, because
    ``overall_status == "pass"`` is the predicate this guard and the published headline
    both use, and ``pass-no-ranking`` is a different (weaker) predicate. So fit-recovery
    passes stay 1/4 and the aggregate stays 5/14.
    That is a real and defensible outcome, but "the guard did not move" is exactly how a
    laundering route hides, so the ``pass-no-ranking`` count is now pinned explicitly below
    (0 -> 2) and the two rows are named. Nobody reading this file can now see 1/4 and be
    unaware that two rows improved.
    AND NOTE THE PRICE, which is not in this test: the same refit made the external hold-out
    WORSE (median 15.31x -> 42.62x, coverage 5/8 -> 4/8). See
    ``test_holdout_scores_1_of_5_genuine_extrapolations_at_the_pre_widening_prior``. The
    in-panel rows that improved are fit-recovery; the out-of-sample rows that degraded are
    the evidence.

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
    assert totals == {"predictive": 9, "fit_recovery": 4, "internal_synthetic": 4}, (
        f"Evidence-role split moved to {totals}, published as 9/4/4. Reclassifying a "
        f"benchmark changes the denominator of the headline claim -- justify it in AUDIT.md."
    )

    def passing(role):
        return sorted(s.benchmark_id for s in by_role[role] if s.overall_status == "pass")

    predictive_passes = passing("predictive")
    assert predictive_passes == [], (
        f"The published headline is 0/9 PREDICTIVE benchmarks without blocking gaps. These "
        f"now pass: {predictive_passes}. Before celebrating, confirm the benchmark is still "
        f"genuinely predictive (its constants were not fitted to it) -- then re-pin here and "
        f"correct README.md and AUDIT.md together."
    )

    # RE-PINNED 2026-08-27 (Wave P item 4): 1 -> 0. THE LAST FIT-RECOVERY PASS IS GONE.
    # CAUSE, and it is the cleanest demonstration this campaign has produced of why fit
    # recovery is not evidence: `pea_isolate_uht_140C_Trikusuma2019` was the single
    # survivor, and it passed because THREE observability factors had been back-solved so
    # it reproduced its own measured 782 / 163 / 24 ppb. Wave P corrected the SUBSTRATE
    # nonanal is cleaved from -- oleate, not linoleate (Miyazaki 2023,
    # 10.1093/bbb/zbac189: nonanal appears in neither linoleate hydroperoxide isomer's
    # product list; `LipidProfile.oleic_acid_pct` had been dead code) -- and the nonanal
    # row immediately fell to 2.2727x under, which is EXACTLY 1 / (22.0 / 50.0), the pea
    # oleic/linoleic ratio. A constant that tracks a substrate error to five significant
    # figures was never measuring the model.
    # The factor was deliberately NOT refitted: doing so would re-absorb the correction
    # into the same constant and make the fix invisible, and Trikusuma 2019 is in any case
    # still the last content-unverified pillar of the matrix lane (Wave O [P] item 5), so
    # there is no verified anchor to refit against. See the dated note on the record in
    # src/matrix_calibration_registry.py.
    assert len(passing("fit_recovery")) == 0, (
        f"fit-recovery passes moved to {len(passing('fit_recovery'))}/4 (published 0/4 "
        f"since the 2026-08-27 Wave P oleate substrate correction). A row returning here "
        f"means either a genuine improvement or a back-solved constant being refitted to "
        f"its own benchmark -- check which before re-pinning."
    )
    assert len(passing("internal_synthetic")) == 4, (
        f"internal-synthetic passes moved to {len(passing('internal_synthetic'))}/4 "
        f"(published 4/4)"
    )

    # 2026-08-27 (Wave O): pinned so a fit-recovery improvement cannot land silently.
    # `pass-no-ranking` means "coverage and scale both pass, but the benchmark declares no
    # ranked target list to check", i.e. a weaker pass than `pass`. It went 0 -> 2 with the
    # refit, and both rows are fit-recovery.
    no_ranking = sorted(s.benchmark_id for s in panel if s.overall_status == "pass-no-ranking")
    assert no_ranking == [
        "pea_isolate_40C_PratapSingh2021",
        "soy_isolate_40C_PratapSingh2021",
    ], (
        f"`pass-no-ranking` rows moved to {no_ranking} (published as the two Pratap-Singh "
        f"fit-recovery rows since the 2026-08-27 Wave O refit). A row appearing here is a "
        f"status improvement that the 0/6 and 5/14 headlines do NOT show; say so in "
        f"README.md and AUDIT.md rather than only re-pinning."
    )
    for benchmark_id in no_ranking:
        role = next(s.evidence_role for s in panel if s.benchmark_id == benchmark_id)
        assert role == "fit_recovery", (
            f"{benchmark_id} improved to pass-no-ranking while classified `{role}`. A "
            f"PREDICTIVE benchmark reaching this status is a genuine result and must be "
            f"promoted deliberately, not absorbed by this guard."
        )

    total_passes = sum(1 for s in panel if s.overall_status == "pass")
    # RE-PINNED 2026-08-27 (Wave P): 5 -> 4, same cause as the fit-recovery line above.
    # EVERY remaining pass is now an internal synthetic reproducibility row, i.e. the model
    # agreeing with its own frozen output. There is no longer a single benchmark in the
    # panel that passes on anything but its own snapshot.
    assert total_passes == 4, (
        f"Aggregate panel passes moved to {total_passes}/14 (published 4/14 since the "
        f"2026-08-27 Wave P oleate substrate correction; 5/14 after Wave K/M, 7/14 before "
        f"that -- explicitly labelled do-not-quote because every pass is in a non-evidence "
        f"bucket)."
    )
    assert passing("internal_synthetic") == sorted(
        s.benchmark_id for s in panel if s.overall_status == "pass"
    ), (
        "a pass appeared outside the internal-synthetic bucket; that is a real result and "
        "must be promoted deliberately, in README.md and AUDIT.md, not absorbed here."
    )

    # 2026-08-27 (Wave P): the PUBLISHED MARKDOWN USES A DIFFERENT PREDICATE, and the two
    # have quietly disagreed since Wave O. `src/presentation.py::_is_pass` counts
    # `pass-no-ranking` and `partial-pass` as passes, so benchmark_summary.md prints
    # "0/6 + 2/4 + 4/4 = 6/14" while this guard, on `overall_status == "pass"`, sees
    # 0/6 + 0/4 + 4/4 = 4/14. Both are defensible; publishing them without saying which is
    # which is not. Pinned here so the divergence cannot widen unnoticed. [P] reconcile.
    lenient_passes = sum(
        1 for s in panel
        if s.overall_status in {"pass", "pass-no-ranking", "partial-pass"}
    )
    assert lenient_passes == 6, (
        f"the lenient (presentation-layer) pass count moved to {lenient_passes}/14, "
        f"published 6/14. It is the number benchmark_summary.md prints; the strict count "
        f"asserted above is the number this guard and the README headline use."
    )

    for doc, path in (("README.md", README), ("AUDIT.md", AUDIT)):
        _assert_quoted(_doc_text(path), "0/6 predictive", doc, "the predictive headline")


# --------------------------------------------------------------------------------------
# 3. Honest external-literature coverage
# --------------------------------------------------------------------------------------


def test_honest_external_literature_coverage_is_4_of_9_on_intervals_too_wide_to_mean_much():
    """EXTERNAL LITERATURE 4/9 evaluable · 4 not evaluable · 2 fitted rows excluded BOTH sides.
    MEDIAN 90% CI WIDTH 2.627 dex -- A FACTOR OF ~424 FROM END TO END. READ THE TWO TOGETHER.

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

    assert summary["benchmark_count"] == 14, (
        f"MC panel is {summary['benchmark_count']} benchmarks, published as 14"
    )
    assert summary["matched_compound_count"] == 41, (
        f"MC panel matched rows moved to {summary['matched_compound_count']}, published as 41"
    )

    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR).
    # 1/3 -> 0/3. THIS IS A DEGRADATION AND IT IS NOT RELAXED HERE: the assertion is still
    # two-sided and exact, so the number cannot drift back in either direction unnoticed.
    # CAUSE, identified per row: the ONE literature hit was
    # `resconi_2023_pbma_beef_identity_benchmark` / furfural, whose `inside_ci` flipped
    # True -> False. Furfural is produced by `Enolisation_1_2`, which Wave S1b connected to
    # the acid-favoured enolisation route-selection term for the first time; its p50 rose
    # 2504.50 -> 3462.43 ppb while its 90% CI NARROWED (1.4048 -> 0.7460 dex), and the
    # measured value fell outside the narrower interval. A narrower interval that now misses
    # every literature row is the worst of both, and that is the honest reading.
    # The COUNTS that did not move -- benchmark_count 11, matched rows 35,
    # not_evaluable 4, excluded_fitted_rows 2 (both would-have-been hits) -- are what
    # identify this as one row leaving its interval rather than an accounting change.
    # RE-PINNED 2026-08-28 (Wave W): 0/3 -> 4/9, and THE INSTRUCTION IN THE OLD MESSAGE WAS
    # FOLLOWED BEFORE RE-PINNING. It said: "If this rises, verify it is the model getting
    # closer and not the interval getting wider -- check median_ci_width_log10 in the same
    # breath." IT IS THE INTERVAL GETTING WIDER. The median 90% CI went 0.746 -> 2.627 dex,
    # i.e. 3.5x wider, from a factor of 5.6 end-to-end to a factor of ~424.
    # CAUSE: six new external-literature rows arrived -- the three Hofmann & Schieberle 1998
    # panel anchors, two compounds each -- and the sulfur branch's barrier priors are so
    # loose that four of the six "cover" their measurement with intervals 2.90, 2.95, 3.04
    # and 3.74 dex wide. An interval spanning three orders of magnitude will contain almost
    # any measurement; those four hits are close to vacuous and are NOT evidence the model
    # improved. The two rows that MISS (fructose FFT and MFT) miss despite intervals 1.44 and
    # 2.63 dex wide, which is the more informative half of the result.
    # NOTHING WAS WIDENED TO ACHIEVE THIS. No prior, threshold or interval was touched in
    # Wave W; the population changed, not the model. The headline is re-pinned with the width
    # assertion below deliberately kept two-sided so the vacuity stays visible.
    assert (coverage["hits"], coverage["total"]) == (4, 9), (
        f"Honest external-literature coverage moved to "
        f"{coverage['hits']}/{coverage['total']}, published as 4/9 (0/3 before Wave W). "
        f"If this rises, verify it is the model getting closer and not the interval getting "
        f"wider -- check `median_ci_width_log10` in the same breath."
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
    # RE-PINNED 2026-08-27 (Wave S1): 0.8495 -> 0.9463 dex, i.e. the external-literature
    # interval got 1.25x WIDER while its coverage stayed at 1/3. CAUSE: the flux propagator
    # became additive over parallel channels, so a compound reached by several routes now
    # samples the barrier uncertainty of ALL of them instead of only its fastest route --
    # more routes, more spread. THIS IS A HONEST WIDENING AND IT IS NOT AN IMPROVEMENT:
    # coverage did not rise with it, so the model is not paying for the extra width. The
    # COUNTS again did not move -- hits 1/3, not_evaluable 4, excluded_fitted_rows 2,
    # benchmark_count 11, matched rows 35 -- which is what identifies this as an
    # interval-width change. Companion widths moved the OTHER way: fitted_row
    # 2.2767 -> 2.2083, internal_synthetic 3.6929 -> 3.5612, because on those rows the
    # extra channels concentrate the allocation rather than spreading it.
    # (History: 0.8558 under Wave O, 0.8495 under Wave P, 0.9463 under Wave S1.)
    # RE-PINNED 2026-08-27 (Wave S1b): 0.9463 -> 0.7460 dex. The interval NARROWED by 1.27x
    # and coverage fell 1/3 -> 0/3 at the same time. Read those two together: this is not
    # the reassuring case (tighter intervals, same coverage) -- it is the interval tightening
    # around a point prediction that moved further from the measurement. Companion widths
    # moved the other way again: fitted_row 2.2083 -> 2.9573, internal_synthetic
    # 3.5612 -> 4.1780.
    # RE-PINNED 2026-08-28 (Wave W): 0.7460 -> 2.6269 dex. See the block above the hits
    # assertion. This is the case the old comment warned about in as many words -- coverage
    # improving while the interval widens -- and it is being reported as that, not as a win.
    assert coverage["median_ci_width_log10"] == pytest.approx(2.6269, abs=5e-4), (
        f"Median CI width moved to {coverage['median_ci_width_log10']:.4f} dex from the "
        f"published 2.627. A coverage rate that improves while this widens is the interval "
        f"getting looser, not the model getting better."
    )
    # The vacuity is asserted, not just described: 4/9 on a ~424x interval must never be
    # quoted as though it were 4/9 on a tight one.
    assert coverage["median_ci_width_log10"] > 2.0, (
        "The median external-literature interval has narrowed below 2 dex. If coverage is "
        "still 4/9 or better at that width, that IS a real improvement -- re-pin this guard "
        "and say so loudly in README.md, because the current prose says the opposite."
    )


# --------------------------------------------------------------------------------------
# 4. The external hold-out
# --------------------------------------------------------------------------------------


def test_holdout_scores_1_of_5_genuine_extrapolations_at_the_pre_widening_prior():
    """HOLD-OUT 1/5 genuine extrapolations pre-widening · median 42.62x · worst 2474x.

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

    RE-PINNED AGAIN 2026-08-27 (Wave O refit to content-corrected anchors, owner-approved),
    median 15.31x -> 42.62x. THE HOLD-OUT GOT WORSE, and this is the wave's headline number.
    The ambient hexanal observability factors were refitted against the CONTENT-VERIFIED
    Pratap-Singh anchors (one shared scale, 4.317249x; record
    results/validation/matrix_observability_refit_pratap_singh.json). Per point:
        Bi 2020 raw pea hexanal   5.37x  ->   1.24x   (IMPROVED)
        Liu 2023 hexanal          4.52x  ->  19.50x   (WORSE)
        Li 2026 hexanal          21.58x  ->  93.15x   (WORSE, via the propagated soy
                                                       heated-matrix baseline)
        the other five points are byte-identical.
    Why a refit to VERIFIED values makes the hold-out worse, stated so it is not read as a
    bug: the pea ambient lane carries two mutually contradictory external measurements --
    Bi 2020 at 1260 ppb and Liu 2023's band midpoint at 51.96 ppb, a 24x spread at nominally
    identical conditions -- and the erroneous 260 ppb the old constants reproduced sat almost
    exactly at their geometric mean (sqrt(1260 * 51.96) = 255.9). Moving onto the verified
    anchor (1138 ppb, which agrees with Bi to 1.11x) necessarily moves off Liu. No single
    observability factor can satisfy both; the disagreement is in the literature, not in the
    fit. What the refit bought is that the constant is now anchored to a number that exists.

    NOTE what did NOT move, because it is the discriminating evidence: the pre-widening
    ``genuine_extrapolation_hits`` is still 1/5 and ``max_fold_error`` is still 2474.4. The
    refit touched one lane, not the model's transfer behaviour.

    2026-08-27 (Wave P item 4) — TWO POINTS IMPROVED AND EVERY PINNED NUMBER BELOW HELD.
    This is the only untuned hold-out movement the campaign has produced, so it is worth
    stating exactly what it was and what it was not. Nonanal is the C9 fragment of the
    OLEATE double bond, not a linoleate product (Miyazaki 2023, 10.1093/bbb/zbac189, read
    in full text: nonanal is in neither linoleate hydroperoxide isomer's product list;
    Hung, Katrib & Martin 2005, 10.1021/jp0500900, "1-nonanal (30 +/- 3% carbon yield)"
    from oleate cleavage). The model computed the whole pool from `linoleic_acid_pct` and
    `LipidProfile.oleic_acid_pct` was dead code. Correcting the SUBSTRATE, with no constant
    fitted and no factor refitted:
        Li 2026 HME nonanal      272.63x -> 118.31x   (IMPROVED, x0.434 = soy oleic/linoleic)
        Liu 2023 PPI nonanal      10.86x ->   4.78x   (IMPROVED, x0.440 = pea oleic/linoleic)
        the other six points are byte-identical.
    And yet: ``median_accuracy_fold`` is UNCHANGED at 42.6159x, ``ci_coverage_hits`` is
    UNCHANGED at 4/8, ``max_fold_error`` is UNCHANGED at 2474.4, and the pre-widening
    genuine-extrapolation count is UNCHANGED at 1/5 — because the median sits between two
    points that did not move and the two that did were already outside the interval and
    remain so. A real, mechanistically-motivated correction improved two of eight points by
    2.3x each and moved the headline by exactly nothing. That is what an honest hold-out
    looks like, and it is the reason the headline is quoted rather than the per-point table.
    Both nonanal points are STILL over-predicted (118x and 4.8x), and the model still treats
    oleate as being as oxidisable as linoleate, which biases them high by roughly another
    order of magnitude — see `src.lipid_oxidation.MARKER_HYDROPEROXIDE_POOL`.

    2026-08-27 (Wave R) — THE LIU HOLD-OUT WAS SCORED AGAINST FABRICATED REFERENCE VALUES.
    The paragraph above quotes "Liu 2023 PPI nonanal 10.86x -> 4.78x (IMPROVED)". Both of
    those folds were computed against a number that does not exist in the cited source. The
    primary document (Yaozheng Liu, "Flavor Chemistry of Pea Proteins", NC State MS thesis
    2021; published as Liu, Cadwallader & Drake 2023, Food Chem. 406:134998) was retrieved
    and read in full. Its Table 2.7 reports, across nine quantified commercial pea proteins,
    hexanal 2445-52454 ug/L and nonanal 0.188-3.42 ug/L of the rehydrated 10%-solids slurry.
    The repo carried hexanal 15-180 ppb (50-300x LOW) and nonanal 5-50 ppb (6-266x HIGH).
    Neither band, and neither of the OAV pairs attached to them, matches any row of any table
    in the thesis. Correcting them (no prediction moved):
        Liu 2023 PPI hexanal   meas 51.96 -> 11320 ppb    fold 19.50x -> 11.17x  (IMPROVED)
        Liu 2023 PPI nonanal   meas 15.81 -> 0.8018 ppb   fold  4.78x -> 94.22x  (WORSE)
    and the headline moved 42.62x -> 93.68x with shipped-sigma coverage 4/8 -> 3/8. The
    nonanal row is now the sharpest lipid-lane over-prediction the repo has against a
    directly-quantified reference: 75.5 ppb predicted against a band whose TOP is 3.42 ppb.
    Wave P's oleate-substrate fix is the partial mitigation already landed — it took this
    same point from 214x to 94x — and it was not enough. NOTE the DI-water calibration
    caveat travelling with the corrected values: Liu's curves were built in deionized water,
    not in the protein matrix, so protein binding is uncorrected and 0.188-3.42 is a LOWER
    bound. The gap is therefore, if anything, understated here.

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

    # RE-PINNED 2026-08-27 (Wave R): 42.62 -> 93.68. NOT a model change -- no constant, prior
    # or factor moved. The Liu 2023 PPI hold-out's two reference values were found to match
    # NOTHING in their cited source and were corrected against the primary document (Yaozheng
    # Liu, "Flavor Chemistry of Pea Proteins", NC State MS thesis 2021, Table 2.7; published
    # as Food Chem. 406:134998). MEASURED values moved, predictions did not:
    #     hexanal  51.96 -> 11320 ppb (band 15-180 -> 2445-52454)   fold 19.50x -> 11.17x
    #     nonanal  15.81 ->  0.8018 ppb (band 5-50 -> 0.188-3.42)   fold  4.78x -> 94.22x
    # The nonanal point is the one that hurts, and it is reported plainly: the model predicts
    # 75.5 ppb against a real inter-lot band whose TOP is 3.42 ppb. The other six points are
    # byte-identical. Wave P's oleate-substrate correction already took this point from 214x
    # to 94x before the reference was fixed, i.e. the mitigation landed first and was not
    # enough. If this number IMPROVES, check it was not bought by reverting either the Wave O
    # refit to the verified anchors or the Wave R reference correction.
    assert summary["median_accuracy_fold"] == pytest.approx(93.68, abs=0.01), (
        f"Median hold-out fold error is {summary['median_accuracy_fold']:.2f}x, published "
        f"as 93.68x since the 2026-08-27 Wave R Liu reference correction (42.62x from the "
        f"Wave O observability refit, 15.31x before it, 32.79x before the Wave K/M Li 2026 "
        f"reference correction). If this IMPROVED, check that it was not bought by reverting "
        f"the refit to the verified anchors or the Liu correction to the verified table."
    )
    # Shipped-sigma coverage: 5/8 -> 4/8 at the Wave O refit, then 4/8 -> 3/8 at the Wave R
    # Liu correction (the nonanal point left the interval; the hexanal point stayed inside
    # it). Pinned here so the headline cannot be quoted from a stale run.
    assert summary["ci_coverage_hits"] == 3
    assert summary["ci_coverage_rate"] == pytest.approx(0.375)
    assert summary["max_fold_error"] == pytest.approx(2474.4, abs=0.1), (
        f"Worst hold-out fold error is {summary['max_fold_error']:.1f}x, published as 2474x."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "1/5", "README.md", "the genuine-extrapolation hold-out result")
    _assert_quoted(readme, "2474", "README.md", "the worst hold-out fold error")


# --------------------------------------------------------------------------------------
# 5. The pentose >> hexose ordering claim, and how little of it is structural
# --------------------------------------------------------------------------------------


def test_pentose_hexose_mft_ordering_is_8_26x_not_the_retired_18_27x_7_78x_6_15x_3_39x_8_98x_or_15_8x():
    """PENTOSE >> HEXOSE 8.26x (ribose 169.1 ppb vs glucose 20.5 ppb at matched conditions).

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

    # RE-PINNED 2026-08-27 (Wave S1 -- THE ADDITIVE FLUX PROPAGATOR). NO BARRIER MOVED.
    # The ratio went UP again, 6.15x -> 7.78x, and this test's own warning applies with
    # more force than ever: DO NOT REPORT THAT AS IMPROVED SUGAR DISCRIMINATION. The
    # pentose limb is reached by more parallel routes than the hexose limb, and a
    # propagator that sums parallel channels rewards exactly that. Re-measured
    # decomposition, in-process, by setting `thiol_addition_pentodiulose` equal to
    # `thiol_addition_hexose` (29.65):
    #     shipped (26.35 vs 29.65)   ribose 824.72 / glucose 105.95 = 7.7838x
    #     equalised                  ribose 332.35 / glucose 105.95 = 3.1368x
    # so 3.14x is now STRUCTURAL and the remaining ~2.5x rides on the 3.30 kcal/mol gap
    # between a FITTED barrier and an UNCONSTRAINED LEGACY FIT. History of the split:
    # 1.13x of 3.39x under Wave N, 2.31x of 6.15x under Wave P, 3.14x of 7.78x now. The
    # structural share HAS grown, and for a defensible reason -- richer topology now
    # reaches the number instead of being discarded by the propagator -- but the gap
    # between two barriers still carries a third of the claim, and the hexose limb still
    # runs the demoted one-step lump.
    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR).
    # NO BARRIER MOVED. 7.78x -> 18.27x, and THIS IS EMPHATICALLY NOT IMPROVED SUGAR
    # DISCRIMINATION -- it is the denominator collapsing. BOTH absolute numbers FELL:
    # ribose 824.7 -> 374.0 ppb, glucose 106.0 -> 20.5 ppb. CAUSE: at aw 0.98 the
    # water-activity correction now reaches the sulfur families, and the hexose limb runs
    # `Thiol_Addition_Hexose_Legacy_Shortcut`, a lumped step that sheds THREE waters, while
    # the pentose limb runs `Deoxyosone_Reduction` + `Thiol_Addition_Pentodiulose`. The
    # hexose route is penalised harder, so the ratio rose while the chemistry got no better.
    # Re-measured decomposition, in-process, setting `thiol_addition_pentodiulose` equal to
    # `thiol_addition_hexose` (29.65):
    #     shipped (26.35 vs 29.65)   ribose 374.03 / glucose 20.47 = 18.2744x
    #     equalised                  ribose  87.31 / glucose 20.47 =  4.2659x
    # History of the structural split: 1.13x of 3.39x (Wave N), 2.31x of 6.15x (Wave P),
    # 3.14x of 7.78x (Wave S1), 4.27x of 18.27x (Wave S1b).
    # RE-PINNED 2026-08-27 (Wave S2c -- THE HOFMANN ANCHOR RETIREMENT). 18.27x -> 8.2607x,
    # ribose 374.0 -> 169.1 ppb, GLUCOSE UNCHANGED at 20.5 ppb. THE RATIO FELL BY MORE THAN
    # HALF AND IS PINNED LOWER; nothing was tuned to hold it up.
    # CAUSE: `thiol_addition_pentodiulose` was REVERTED 26.35 -> 28.60 (the un-fitted Wave N
    # class value) because the Wave P refit that produced 26.35 had exactly one fit target,
    # `cys_ribose_140C_Hofmann1998`, and Wave S2b showed that benchmark's MFT 342 / FFT 200 ppb
    # are a REPO-INTERNAL DERIVATION -- interior points of two invented mol % bands in
    # data/benchmarks/maillard_validation_benchmarks.md section 1.3, committed in the SAME
    # commit as the benchmark JSON. Only the PENTOSE limb runs that barrier, which is exactly
    # what identifies the cause: glucose does not move at all.
    # THE CLAIM IS SMALLER BUT ITS EVIDENCE IS BETTER, and both halves must be reported:
    #     shipped (28.60 vs 29.65)   ribose 169.08 / glucose 20.47 =  8.2607x
    #     equalised (both 29.65)     ribose  87.31 / glucose 20.47 =  4.2659x
    # so the STRUCTURAL share is unchanged at 4.27x while the total fell, i.e. the fraction of
    # the ordering carried by mechanism rather than by a barrier gap went from 23% (4.27 of
    # 18.27) to 52% (4.27 of 8.26). The residual gap is now the 1.05 kcal/mol between an
    # ESTIMATED class value and an unconstrained legacy fit, not the 3.30 kcal/mol between a
    # FITTED barrier and an unconstrained legacy fit -- and no part of the claim now traces to
    # a number this repository invented. The hexose limb still runs the demoted one-step lump.
    # STILL BELOW THE 8.98x IT SAT AT AFTER WAVE N, and the 3.0x floor in
    # test_pentose_hexose_sulfur_ordering.py is now 2.75x away rather than 6x away.
    assert ribose_ppb == pytest.approx(169.1, rel=0.01), (
        f"Ribose MFT moved to {ribose_ppb:.1f} ppb from the published 169.1 (374.0 under "
        f"Wave S1b, 824.7 under Wave S1, "
        f"686.8 after the Wave P refit, 370.3 after the Wave N route correction, "
        f"981.3 before it)"
    )
    assert glucose_ppb == pytest.approx(20.5, rel=0.01), (
        f"Glucose MFT moved to {glucose_ppb:.1f} ppb from the published 20.5 (106.0 under "
        f"Wave S1, 111.6 "
        f"under Wave P, 109.3 before it). The hexose limb still runs the demoted one-step "
        f"lump; it FELL under Wave S1 because the additive propagator gives the parallel "
        f"sulfur channels a larger share of the same fixed volatile budget, so a LARGE "
        f"move here has a different cause than a move in the ribose value. Wave S2c's "
        f"barrier revert left it untouched, which is what identified the cause."
    )
    assert ratio == pytest.approx(8.26, rel=0.01), (
        f"Pentose/hexose MFT ratio is {ratio:.2f}x, published as 8.26x. If it moved, do "
        f"not report it as changed sugar discrimination without first re-measuring the "
        f"structural share -- only 4.27x of this ratio survives setting "
        f"`thiol_addition_pentodiulose` equal to `thiol_addition_hexose`; the rest is the "
        f"gap between an estimated class value and an unconstrained legacy fit. And if it "
        f"ROSE, check first whether a constant was refitted against "
        f"cys_ribose_140C_Hofmann1998, whose values are not measurements."
    )

    for doc, path in (("README.md", README), ("AUDIT.md", AUDIT)):
        _assert_quoted(_doc_text(path), "8.26", doc, "the pentose/hexose ordering margin")


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


def test_no_verifiable_source_census_is_120_records_98_numeric_80_reaching_runtime():
    """no_verifiable_source: 120 flagged · 98 carrying numbers · 80 of those reaching runtime.

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

    RE-PINNED 2026-08-27 (Wave T3): 102/80/62 -> 120/98/80. This is the SECOND rise in one
    day and, like the first, it is a labelling change, not a data change. Eighteen records
    that were already shipping unverifiable numbers were finally labelled:

      * 15 in `data/lit/protein_source_registry.json` -- the file-level record plus all 14
        protein-source profiles. That file has always described itself as "Mocked values for
        14 protein sources based on Report 06 requirements", but the sentence sat in a JSON
        field nothing read while the numbers underneath it drove `matrix_uncertainty_factor`
        and the meaty-potential score at prediction time (Wave T1 finding T1-01).
      * 2 in `data/lit/retention_reference_payloads.json` -- the two `runtime_surrogate`
        blocks whose `log_slope = 0.235` is exactly ln(1.60)/2, back-solved from an invented
        "~55-65%" band in an in-repo brief (T1-02).
      * 1 in `data/lit/computational_priors.json` -- `ref41_ppi_sulfur_volatile_binding_v1`,
        which cited reference *number* 41 inside an LLM research dump (T1-08).

    NO VALUE WAS ADDED, CHANGED OR INVENTED IN THAT WAVE. All 18 are in data/lit and all
    carry numbers, so the runtime figure moves with the total: 62 -> 80. The honest reading
    is that 80 was always the true runtime figure and 62 was an undercount.

    All three counts are expected to FALL as anchors get verified. When they do, re-pin here
    and in README.md in the same change. A rise is only acceptable when it is a LABELLING
    correction of numbers that were already shipping, and the re-pin must say which records
    moved and why -- a rise with no such account means new unverifiable numbers entered the
    registries.
    """
    census = _no_verifiable_source_census()

    total = sum(flagged for flagged, _ in census.values())
    numeric = sum(with_number for _, with_number in census.values())
    runtime = sum(
        with_number
        for path, (_, with_number) in census.items()
        if path.startswith("data/lit/")
    )

    assert total == 120, (
        f"Repo-wide no_verifiable_source count is {total}, published as 120. "
        f"Per file: { {k: v[0] for k, v in census.items()} }"
    )
    assert numeric == 98, (
        f"{numeric} flagged records carry numeric payloads, published as 98."
    )
    assert runtime == 80, (
        f"{runtime} flagged records with numeric payloads sit in data/lit and therefore "
        f"reach the runtime, published as 80. If a data/qm record becomes runtime-reachable "
        f"this number must rise and the README's 'none of those 18 reach the model' claim "
        f"becomes false."
    )

    readme = _doc_text(README)
    _assert_quoted(readme, "120 records", "README.md", "the no_verifiable_source census")
    _assert_quoted(readme, "98 carry numeric payloads", "README.md", "the numeric subset")
    _assert_quoted(readme, "80 of those are", "README.md", "the runtime-consumed subset")


def test_the_data_qm_records_that_moved_the_census_from_84_to_102_are_still_counted():  # noqa: N802
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
    # RE-PINNED 2026-08-27 (Wave T3): 84 -> 102. The data/qm side of the split is UNCHANGED
    # at 18, which is what this test exists to protect. The data/lit side rose because Wave
    # T3 labelled 18 already-shipping records (15 protein_source_registry, 2
    # retention_reference_payloads runtime_surrogate blocks, 1 ref41) -- see the census test
    # above for the full account. Total is now 102 + 18 = 120.
    assert in_lit == 102, (
        f"data/lit carries {in_lit} flagged records, measured as 102 on 2026-08-27 after "
        f"Wave T3's labelling pass (it was 84 before that pass, and 84 is still the figure "
        f"the README published before data/qm was brought into scope)."
    )
