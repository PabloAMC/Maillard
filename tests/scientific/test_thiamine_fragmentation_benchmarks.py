import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark  # noqa: E402


def test_bolton_1994_thiamine_fragmentation_benchmark_is_executable_and_under_predicts():
    # Replaces test_hofmann_1996_...: that benchmark was deleted on 2026-08-26 after source
    # recovery confirmed its DOI (10.1021/jf960062o) was never registered with any DOI agency
    # and its reported MFT value 3.14 ppb is numerically pi. It was rebuilt from the verified
    # Bolton et al. 1994 chapter (10.1021/bk-1994-0543.ch022) as
    # data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json, which IS in the panel.
    #
    # The executability + coverage assertions carry over unchanged: they exercise the thiamine
    # fragmentation wiring (that thiamine + cysteine + sugar can reach MFT at all).
    #
    # The numeric assertion is INVERTED relative to the old test, and deliberately so. Against
    # a real reference value the model under-predicts MFT here by a large factor.
    # The old test's "ratio <= 1.5" band passed only because both the reference value and the
    # tolerance had been chosen to make it pass. The bound below is a drift guard on a KNOWN
    # FAILURE: it pins the gap loosely so that neither a silent regression nor a silent
    # "improvement" goes unnoticed, and it fails loudly if the gap ever closes -- at which
    # point this test should be rewritten as a passing one, not loosened.
    #
    # 2026-08-27 (projection budget retune + Boltzmann de-duplication, audit remediation
    # Part 1): the gap moved 162.8x -> 173.3x (predicted 0.0799 -> 0.0750 ppb against the
    # same 13 ppb target). This is a side effect of the allocation change, NOT a tuning
    # target -- nothing in the retune was aimed at this benchmark.
    #
    # 2026-08-27 (Wave G1 chemistry rebuild + Wave H sulfur refit): the gap moved again,
    # 173.3x -> 744.9x (predicted 0.0750 -> 0.01745 ppb). Cause, in order of size: G1
    # deleted the fabricated lipid radical chain and rebuilt the thiamine cascade, and the
    # MFT route moved off the fabricated one-step shortcut onto the real norfuraneol route
    # (van den Ouweland & Peer 1975), which cost the whole sulfur branch 5-40x of absolute
    # yield; Wave H then refit `thiol_addition_norfuraneol` against Hofmann1998 alone, which
    # gave a little of that back on the ribose lanes but essentially nothing here. Again a
    # side effect: nothing in either wave was aimed at Bolton, and no contract tolerance in
    # data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json was touched, so this stays a
    # loud failure in every regenerated artifact. The guard below is RE-PINNED to the new
    # magnitude rather than left at a bound so loose it would have stopped detecting drift.
    evaluation = evaluate_benchmark(
        ROOT / "data" / "benchmarks" / "thiamine_cys_glucose_120C_Bolton1994.json"
    )

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "2-Methyl-3-furanthiol (MFT)"
    assert comparison.measured_ppb == 13.0
    assert 0.0 < comparison.predicted_ppb < 13.0, (
        "The Bolton 1994 anchor is a known under-prediction. If predicted MFT has reached or "
        f"exposed the 13 ppb target (got {comparison.predicted_ppb:.4f} ppb), the Family 03 "
        "thiamine gap has closed and this test must be rewritten as a pass, not relaxed."
    )
    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR). NO
    # CONSTANT MOVED. The gap moved 748.0x -> 6730.9x (predicted 0.01738 -> 0.001931 ppb
    # against the same 13 ppb target). CAUSE, ATTRIBUTED BY MEASUREMENT rather than assumed:
    # the three fixes were run one at a time (by emptying the other two family sets at
    # runtime, no source edit) and this row is carried ALMOST ENTIRELY by the water-activity
    # fix -- none 0.01738, route-pH-only 0.01659, pyrazine-only 0.01744, aw-only 0.002055,
    # all three 0.001931 ppb. MECHANISM: `_water_activity_correction` used to reach 3 of the
    # 29 emitted families and its dehydration branch keyed on the substring "furfural",
    # which matches none of them, so it was dead. Membership is now decided by MEASURED net
    # water stoichiometry, and at this benchmark's aw 0.98 every water-shedding step takes
    # 1.3 - 0.98 = 0.32, i.e. +0.89 kcal/mol of effective barrier. The thiamine -> MFT lane's
    # terminal `Furan_Ring_Aromatisation` (net +1 water) takes that penalty while the
    # `Additive_Thermal_Degradation` steps above it do NOT -- that family is the only
    # stoichiometrically non-uniform one (+2/0/-1/-2 across its steps) and is deliberately
    # excluded, because one family-level factor cannot honestly represent it. MFT is a tiny
    # share of a large fixed budget here, so its competitors absorb what it loses and the
    # fold error amplifies. NOTHING WAS CLAWED BACK and no contract tolerance in
    # data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json was touched, so this stays a
    # loud failure in every regenerated artifact -- it is now the WORST quantitative point in
    # the whole panel, having overtaken the CML row at 1203.7x.
    assert 3000.0 < comparison.ratio < 14000.0, (
        f"Known Family 03 gap is {comparison.ratio:.1f}x; the pinned post-Wave-S1b value is "
        "6730.9x (744.9x post-G1/Wave-H, 748.0x post-Wave-S1). A move outside [3000, 14000] "
        "means something other than the 2026-08-27 approved changes has moved the thiamine "
        "lane -- investigate before re-pinning, and if the gap has genuinely closed, rewrite "
        "this test as a pass."
    )


def test_cerny_2008_thiamine_fragmentation_benchmark_is_executable_and_under_predicts():
    # RE-PINNED AGAIN 2026-08-27 (Wave H), and RENAMED a second time: this benchmark has
    # flipped from over- to UNDER-prediction. Its own test docstring said the rule for this
    # situation -- "re-derive, not relax" -- so here is the derivation rather than a widened
    # band.
    #   (3) 2026-08-27, Wave G1 chemistry rebuild: MFT stopped being made by the fabricated
    #       one-step 3-deoxyosone + 2 H2S shortcut and now comes from the accepted
    #       1-deoxyosone -> norfuraneol -> MFT route (van den Ouweland & Peer 1975), the
    #       fabricated lipid radical chain (93% of the emitted network) was deleted, and the
    #       thiamine cascade was rebalanced. Sulfur absolute yields fell 5-40x branch-wide.
    #       Cerny went from 3.419x OVER to under-predicting.
    #   (4) 2026-08-27, Wave H: `thiol_addition_norfuraneol` refit 28.60 -> 26.85 against
    #       cys_ribose_140C_Hofmann1998 ONLY (the sole surviving literature constraint on the
    #       sulfur branch). Cerny was NOT a fit target and must not become one -- see the
    #       standing caveat below. Final: predicted 0.7730 ppb vs 2.47 ppb = 3.195x UNDER.
    #       CORRECTED 2026-08-27 (Wave S2c): the parenthesis above is FALSE and is kept
    #       because it records what the repo believed. cys_ribose_140C_Hofmann1998 was never
    #       a literature constraint -- Wave S2b traced its MFT 342 / FFT 200 ppb to
    #       data/benchmarks/maillard_validation_benchmarks.md section 1.3, an
    #       abstract-reconstructed table committed in the SAME commit as the benchmark JSON,
    #       with both values interior points of two invented mol % bands. THE SULFUR BRANCH
    #       HAS ZERO ABSOLUTE LITERATURE ANCHORS. Note the compounding: this benchmark's own
    #       2.47 ppb reference is ALSO unverified (see the standing caveat below), so the
    #       Family 03 lane is now measured against one unverified value using a barrier that
    #       was fitted to a fabricated one and has since been reverted.
    # The assertion below is inverted to match the measured outcome and pinned to the
    # observed magnitude in both directions. Neither this benchmark's contract tolerance nor
    # any other was touched, so the failure remains visible in the regenerated artifacts.
    #
    # Earlier history, kept because it is the record of how this expectation moved:
    # this test originally asserted `ratio <= 1.5`; that band failed under two separately
    # owner-approved changes:
    #   (1) 2026-08-26, thiol_addition barrier re-centred 28.85 -> 28.60 kcal/mol (the
    #       Hofmann1998-only optimum, after the Mottram1994/Farmer1999 fitting targets were
    #       quarantined): ratio 1.05 -> 1.641.
    #   (2) 2026-08-27, projection budget retune + Boltzmann de-duplication (audit
    #       remediation Part 1): ratio 1.641 -> 3.419. Removing the double-applied
    #       Boltzmann raises MFT's share of the volatile budget across the whole
    #       free-precursor panel; Cerny over-predicts as a direct consequence.
    # The expectation is UPDATED, not widened away: the benchmark's own contract tolerance in
    # data/benchmarks/thiamine_cys_xylose_145C_Cerny2008.json is untouched, so the panel
    # summary still reports this benchmark as `scale-gap` and the failure stays visible in
    # every regenerated artifact. This test pins the observed magnitude as a drift guard
    # in both directions.
    # Standing caveat: this benchmark carries source_metadata.VALUES_NEED_RE_EXAMINATION --
    # the repaired DOI (10.1021/jf801762c) is a paper about 5-hydroxy-3-mercapto-2-pentanone,
    # not MFT, so the 2.47 ppb reference is not yet verified against the real source and must
    # not be used to constrain the sulfur barriers.
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "thiamine_cys_xylose_145C_Cerny2008.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert evaluation.reference_signal_origin == "reference_volatiles"
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "2-Methyl-3-furanthiol (MFT)"
    assert 0.0 < comparison.predicted_ppb < comparison.measured_ppb, (
        "Cerny 2008 is currently an UNDER-prediction (it was an over-prediction before the "
        f"2026-08-27 Wave G1 chemistry rebuild). Got {comparison.predicted_ppb:.4f} ppb vs "
        f"{comparison.measured_ppb} ppb. If the model has swung back to over-predicting, the "
        "Family 03 sulfur allocation has changed character again and this test must be "
        "re-derived, not relaxed."
    )
    # RE-PINNED 2026-08-27 (Wave S1b -- THE pH / WATER-ACTIVITY ROUTING REPAIR). NO
    # CONSTANT MOVED. 2.787x -> 23.406x (predicted 0.7730 -> 0.1055 ppb vs 2.47 ppb). SAME
    # CAUSE AND SAME LANE as the Bolton row above -- this is the second thiamine benchmark
    # and it moves for the identical reason, which is itself the evidence that the movement
    # is the aw routing repair and not something benchmark-specific. The direction did not
    # change: still an UNDER-prediction, asserted separately above. Not clawed back.
    # RE-RECORDED 2026-08-27 (Wave S2c -- THE HOFMANN ANCHOR RETIREMENT). 23.406x -> 25.741x
    # (predicted 0.1055 -> 0.0959 ppb vs 2.47 ppb). CAUSE: `thiol_addition_pentodiulose`
    # REVERTED 26.35 -> 28.60 because the Wave P refit that produced 26.35 was fitted against
    # cys_ribose_140C_Hofmann1998 alone, whose values are not measurements. The xylose lane
    # runs the same pentodiulose step as the ribose lane, which is why this row moves at all.
    # WORSE, AND NOT CLAWED BACK. THE BAND IS DELIBERATELY NOT WIDENED: [12, 46] was set
    # around 23.406x and 25.741x sits comfortably inside it, so the guard keeps exactly the
    # sensitivity it had. Only the recorded value below is updated.
    assert 12.0 < comparison.ratio < 46.0, (
        f"Cerny 2008 MFT fold error is {comparison.ratio:.3f}x; the pinned post-Wave-S2c "
        "value is 25.741x (under); it was 23.406x post-Wave-S1b, 3.195x post-G1/Wave-H and "
        "2.787x post-Wave-S1. A "
        "move outside [12, 46] means something other than the 2026-08-26/27 approved changes "
        "has moved the sulfur branch -- investigate before re-pinning."
    )