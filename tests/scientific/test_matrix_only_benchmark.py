import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import (
    evaluate_benchmark,
    get_matrix_ranking_contract,
    snapshot_all_benchmark_targets,
    snapshot_benchmark_targets,
    summarize_evaluation,
)


MATRIX_ONLY_BENCHMARKS = {
    "pea_iso": ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json",
    "soy_iso": ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json",
}
TRIKUSUMA_BENCHMARK = ROOT / "data" / "benchmarks" / "pea_isolate_uht_140C_Trikusuma2019.json"


@pytest.mark.parametrize("bench_file", MATRIX_ONLY_BENCHMARKS.values())
def test_matrix_only_benchmark_is_executable_with_full_coverage(bench_file):

    evaluation = evaluate_benchmark(bench_file)

    assert evaluation.supported is True
    assert evaluation.coverage == 1.0
    # RE-PINNED 2026-08-27 (Wave M): 3 -> 2 comparisons. CAUSE: the Wave K/M content
    # correction removed the hexanol row from both Pratap-Singh benchmarks. Molecules 2021,
    # 26, 4104 Table 1 (Europe PMC, PMC8271896) reports hexanol as *n.d.* for pea and for
    # soy, and the paper's text states pea proteins "contained no alcohol compounds"; soy's
    # entire alcohol fraction is 40 +/- 9 ppb of 1-octen-3-ol. The 80 / 120 ppb rows this
    # repo scored against had no source. See each file's `content_correction_note`.
    assert len(evaluation.comparisons) == 2


# RE-DERIVED 2026-08-27 (Wave M). CAUSE: the Wave K/M content correction of both
# Pratap-Singh benchmarks against Molecules 2021, 26, 4104 Table 1 (Europe PMC, PMC8271896).
# Three things changed and each is recorded here rather than folded away:
#
#   1. The hexanol row is GONE (paper: n.d. in both matrices), so the three-way ordering
#      assertion has only two members left.
#   2. The measured hexanal is 1138.00 ppb (pea) and 1621.71 ppb (soy), not 260 / 380. The
#      observability factors were back-solved from 260 / 380 and are UNCHANGED (refitting
#      them is an owner science decision -- AUDIT.md, Round 3), so the lane now misses by
#      4.366x and 4.269x, i.e. exactly the size of the transcription error.
#   3. The pea ORDERING claim inverted. With the paper's real values hexanal (1138) outranks
#      2-pentylfuran (638), but the model still predicts 2-pentylfuran (638.3) above hexanal
#      (260.6). That is a genuine, newly visible ranking failure -- previously masked because
#      the reference agreed with the constant fitted to it -- and it is asserted here as a
#      failure rather than removed. Soy's ordering (2-pentylfuran > hexanal) is unchanged and
#      still correct.
#
# NOTHING WAS RELAXED. The `hexanal` limits below were <=1.25 / <=1.10 pass-bands around an
# algebraic recovery; they are now two-sided pins on an honest miss, which is a STRICTER
# contract (drift in either direction fails). The 2-pentylfuran limits are untouched.
#
# RE-DERIVED AGAIN 2026-08-27 (Wave O refit to content-corrected anchors, owner-approved).
# The owner decision item 2 above flags has been taken. The ambient hexanal observability
# factors were refitted against the VERIFIED 1138.00 / 1621.71 ppb using ONE shared scale of
# 4.317249x (scripts/generators/refit_matrix_observability_pratap_singh.py ->
# results/validation/matrix_observability_refit_pratap_singh.{json,md}). Consequences pinned
# below:
#   * hexanal ratios 4.366x / 4.269x -> a COMMON 1.0113x. Not 1.0000x: one shared scale
#     against two rows leaves a residual, and that residual is the only informative number
#     in the fit -- it says the two corrected anchors agree to 1.1%, i.e. the transcription
#     error was a common absolute-scale error and the pea-vs-soy release ratio survived it.
#   * item 3's ordering inversion is RESOLVED for pea: the model now predicts hexanal
#     (1125.3) above 2-pentylfuran (638.3), matching the paper. Read as fit recovery -- the
#     winning compound is the one whose constant was solved from its own measurement.
#   * item 2's "under-prediction" direction now SPLITS: pea under, soy over, because a single
#     shared scale must bracket the two per-lane required scales (4.36606 and 4.26899). That
#     bracketing is asserted below and is the fingerprint of a one-parameter fit.
# Still nothing relaxed: same tolerances, same strict-gate exclusion, and both benchmarks
# stay `fit_recovery` in the evidence-role split.
@pytest.mark.parametrize(
    ("protein_type", "bench_file", "expected_limits", "expects_hexanal_first"),
    [
        (
            "pea_iso",
            MATRIX_ONLY_BENCHMARKS["pea_iso"],
            {"hexanal": 1.0113, "2-pentylfuran": 1.2},
            True,
        ),
        (
            "soy_iso",
            MATRIX_ONLY_BENCHMARKS["soy_iso"],
            {"hexanal": 1.0113, "2-pentylfuran": 1.05},
            False,
        ),
    ],
)
def test_matrix_only_benchmark_preserves_measured_ordering_without_entering_strict_gate(
    protein_type,
    bench_file,
    expected_limits,
    expects_hexanal_first,
):

    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type=protein_type)
    ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}
    predicted = {comparison.compound: comparison.predicted_ppb for comparison in evaluation.comparisons}
    measured = {comparison.compound: comparison.measured_ppb for comparison in evaluation.comparisons}

    assert set(predicted) == {"hexanal", "2-pentylfuran"}

    # The measured ordering the paper reports.
    assert (measured["hexanal"] > measured["2-pentylfuran"]) is expects_hexanal_first

    # 2026-08-27 (Wave O): the ordering the model produces now AGREES with the paper on both
    # lanes -- pea puts hexanal first, soy puts 2-pentylfuran first. Asserted as agreement
    # with the measured ordering rather than as a fixed compound order, so the test states
    # the property instead of a coincidence. It is fit recovery on the pea side: the compound
    # that moved to the top is the one whose observability factor was solved from its own
    # measured value.
    assert (predicted["hexanal"] > predicted["2-pentylfuran"]) is expects_hexanal_first, (
        "the predicted ordering no longer matches the measured ordering; before re-pinning, "
        "check whether the Wave O refit was reverted."
    )

    assert ratios["2-pentylfuran"] <= expected_limits["2-pentylfuran"]
    assert ratios["hexanal"] == pytest.approx(expected_limits["hexanal"], rel=0.01), (
        f"{protein_type} hexanal ratio is {ratios['hexanal']:.4f}x, pinned at "
        f"{expected_limits['hexanal']}x -- the residual left by the 2026-08-27 Wave O "
        "one-parameter refit against the content-verified anchors. 1.0000 would mean someone "
        "re-added per-lane freedom; 4.37x/4.27x would mean the refit was reverted."
    )
    # 2026-08-27 (Wave O): direction is now lane-specific and is the fingerprint of the
    # one-parameter fit. The shared scale 4.317249x sits between the pea lane's required
    # 4.36606x and the soy lane's 4.26899x, so it MUST under-shoot pea and over-shoot soy.
    if protein_type == "pea_iso":
        assert predicted["hexanal"] < measured["hexanal"], (
            "pea hexanal must sit just UNDER its anchor under a shared-scale fit."
        )
    else:
        assert predicted["hexanal"] > measured["hexanal"], (
            "soy hexanal must sit just OVER its anchor under a shared-scale fit."
        )
    assert summary.strict_ready is False


def test_matrix_only_benchmark_is_deliberately_excluded_from_target_snapshots():
    bench_file = MATRIX_ONLY_BENCHMARKS["pea_iso"]
    soy_bench_file = MATRIX_ONLY_BENCHMARKS["soy_iso"]

    assert snapshot_benchmark_targets(bench_file) == []
    assert snapshot_benchmark_targets(soy_bench_file) == []

    rows = snapshot_all_benchmark_targets([
        bench_file,
        soy_bench_file,
        ROOT / "data" / "benchmarks" / "cys_ribose_140C_Hofmann1998.json",
    ])

    assert rows
    assert {row.benchmark_id for row in rows} == {"cys_ribose_140C_Hofmann1998"}


def test_matrix_only_benchmark_exposes_ranking_contract_and_calibration_metadata():
    bench_file = MATRIX_ONLY_BENCHMARKS["soy_iso"]

    evaluation = evaluate_benchmark(bench_file)
    summary = summarize_evaluation(evaluation, protein_type="soy_iso")
    contract = get_matrix_ranking_contract(bench_file)

    assert contract["process_state"] == "ambient_slurry"
    assert contract["calibration_mode"] == "compound_specific_headspace"
    # RE-PINNED 2026-08-27 (Wave M): the "hexanol" target is gone. CAUSE: the Wave K/M
    # content correction -- Molecules 2021, 26, 4104 Table 1 reports hexanol as n.d. for soy
    # (and the paper's whole soy alcohol fraction is 40 +/- 9 ppb of 1-octen-3-ol, i.e. 3x
    # LESS than the 120 ppb this file had attributed to a compound the paper says was not
    # detected). The soy ORDER (2-pentylfuran > hexanal) is unchanged and still correct.
    assert [item["name"] for item in contract["observable_targets"]] == [
        "2-pentylfuran",
        "hexanal",
    ]

    hexanal_meta = evaluation.projection_metadata["Hexanal"]
    # RE-PINNED 2026-08-27 (Wave I): `literature_anchored` -> `fitted_to_benchmark`.
    # This is a relabel, not a value change -- the factor is still 1.0. The pea ambient lane
    # is the REFERENCE lane, and the base marker yields it multiplies
    # (benchmark_validation.MATRIX_BENCHMARK_BASE_MARKER_YIELDS) are this very benchmark's
    # own measured ppb divided by a single common scale: 0.205 x 1268.3 = 260 ppb,
    # 0.502 x 1270.9 = 638, 0.063 x 1269.8 = 80. So the Pearson 1.000 / max_ratio 1.002 this
    # benchmark scores is the lane reproducing the numbers it was built from. Calling that
    # `literature_anchored` and then counting the agreement as validation is the circularity
    # the cold-start red team found; `fitted_to_benchmark` says what it is.
    # See src/matrix_calibration_registry.py (the arithmetic is written out there).
    #
    # 2026-08-27 (Wave M) -- and the label is now doing even more work than it was. Wave K
    # read the paper: the 260 / 380 hexanal and the 80 / 120 hexanol those factors were
    # solved from are NOT in it (true hexanal 1138.00 / 1621.71; hexanol n.d.). So this lane
    # is fitted to values with no source, and the max_ratio it now scores is 4.366x rather
    # than 1.002x. `fitted_to_benchmark` remains the correct label; the factors were NOT
    # refitted (owner decision -- AUDIT.md, Round 3).
    #
    # 2026-08-27 (Wave O refit to content-corrected anchors, owner-approved): the factors HAVE
    # now been refitted, against the verified values this time, and the label is UNCHANGED --
    # which is exactly the point. A constant solved from a benchmark is fit-recovery whether
    # the value it was solved from was right or wrong; refitting improved the arithmetic and
    # changed nothing about the evidence. The assertion below is therefore identical before
    # and after the refit, and that stability is deliberate.
    assert hexanal_meta["calibration_evidence_strength"] == "fitted_to_benchmark"
    assert hexanal_meta["calibration_fallback_mode"] == "compound_specific"
    assert hexanal_meta["process_state"] == "ambient_slurry"
    assert hexanal_meta["evidence_state"] == "externally_benchmarked"
    assert hexanal_meta["target_class"] == "adverse_lipid_markers"
    assert summary.ranking_contract_status == "pass"
    assert summary.adverse_markers == ["2-pentylfuran", "hexanal"]


def test_trikusuma_heated_pea_matrix_fit_is_recovered_to_within_1_05x():
    """FIT RECOVERY, NOT PREDICTIVE ACCURACY. The <=1.05x below is expected BECAUSE fitted.

    RENAMED 2026-08-27 (Wave J2, red-team finding: deceptive test names). The previous name
    was ``test_trikusuma_heated_pea_matrix_benchmark_is_quantitatively_supported``, which
    reads as "the model quantitatively predicts this literature benchmark". It does not. The
    heated-matrix projection constants this test scores against were BACK-FITTED to this very
    benchmark -- see ``results/validation/projection_constant_refit.json`` and the
    ``fitted_to_benchmark`` calibration_evidence_strength asserted a few tests above, where
    the arithmetic (measured ppb / a single common scale) is written out in full.

    So a tight ratio here is arithmetic, not agreement: the lane is reproducing the numbers
    it was solved from. Recovering a fit to within 5% is a REGRESSION CHECK on the fitting
    machinery -- it catches an optimiser, registry or unit change that stops the constants
    reproducing their own target -- and it is worth having for exactly that. It is NOT
    evidence of predictive accuracy and must never be counted as a validation hit. The panel
    accounting already excludes it on both sides of the ratio: this benchmark sits in the
    ``fit_recovery`` evidence_role bucket, and the honest headline is 0/6 on the PREDICTIVE
    bucket (see ``tests/scientific/test_honest_headline_guards.py``).

    Note the last assertion, which is the honest half of this test and stays: even at 1.02x
    recovery, ``strict_ready`` is False.
    """
    evaluation = evaluate_benchmark(TRIKUSUMA_BENCHMARK)
    summary = summarize_evaluation(evaluation, protein_type="pea_iso")
    predicted = {comparison.compound: comparison.predicted_ppb for comparison in evaluation.comparisons}
    ratios = {comparison.compound: comparison.ratio for comparison in evaluation.comparisons}

    assert predicted["hexanal"] > predicted["2-pentylfuran"] > predicted["nonanal"]
    assert ratios["hexanal"] <= 1.05
    assert ratios["2-pentylfuran"] <= 1.05
    assert ratios["nonanal"] <= 1.05
    assert summary.process_state == "heated_matrix"
    assert summary.ranking_contract_status == "pass"
    assert summary.overall_status == "pass"
    assert summary.strict_ready is False