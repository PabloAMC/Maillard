from pathlib import Path

from src.benchmark_validation import evaluate_benchmark
from src.safety import predict_furosine


ROOT = Path(__file__).resolve().parents[2]


def _comparison_map(evaluation):
    return {comparison.compound: comparison for comparison in evaluation.comparisons}


def test_hvp_xylose_benchmarks_are_executable_and_rank_wheat_above_spi():
    spi_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "spi_hvp_xylose_120C_PMC9905368.json")
    wheat_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "wheat_gluten_hvp_xylose_120C_PMC9905368.json")

    assert spi_eval.supported, spi_eval.reason
    assert wheat_eval.supported, wheat_eval.reason
    assert spi_eval.coverage == 1.0
    assert wheat_eval.coverage == 1.0

    spi_mft = _comparison_map(spi_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    wheat_mft = _comparison_map(wheat_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    assert wheat_mft > spi_mft

    # RE-PINNED 2026-08-27 (cause: projection budget retune + Boltzmann de-duplication,
    # audit remediation Part 1). These two benchmarks were inside 1.5x before the retune
    # (SPI 2.465 -> was 1.5-band at an earlier calibration, wheat 2.128) and are now at
    # 5.165x and 5.170x respectively. The driver is the Boltzmann de-duplication: the
    # projection layer was applying exp(-dspan/RT) twice -- once explicitly as
    # `span_activity` at a 0.65*RT window, and again inside the pathway flux it
    # multiplied by -- for a net selectivity evaluated at T/2.19. Removing the duplicate
    # flattens the allocation, which lifts MFT's share and drops methional's on exactly
    # these hydrolysate systems (methional 1.93 -> 0.36 ppb on SPI, 3.61 -> 0.67 on
    # wheat; MFT 0.44 -> 0.61 and 0.72 -> 1.36).
    # The ORDERING claim this test is named for (wheat > SPI) is unaffected and still
    # asserted above. The magnitude expectation is updated rather than deleted, and the
    # benchmarks' own contract tolerances are untouched -- both now report `ranking-gap`
    # in the panel summary and both LOST strict-ready status (6/16 -> 4/16), so the
    # regression stays loud in every regenerated artifact.
    # RE-PINNED AGAIN 2026-08-27 (Wave H). Causes, in order of size:
    #   (1) Wave G1 chemistry rebuild -- MFT moved off the fabricated one-step
    #       3-deoxyosone + 2 H2S shortcut onto the accepted 1-deoxyosone -> norfuraneol ->
    #       MFT route, and the fabricated lipid radical chain (93% of the emitted network)
    #       was deleted. Sulfur absolute yields fell 5-40x branch-wide; these two lanes
    #       went from ~5x to 40-95x under.
    #   (2) Wave H sulfur refit -- `thiol_addition_norfuraneol` 28.60 -> 26.85, fitted
    #       against cys_ribose_140C_Hofmann1998 ONLY. Gave back ~1.9x of MFT here as a
    #       side effect; these benchmarks were NOT fit targets.
    #   (3) Wave H observability re-derivation -- the Methional `base_factor` was
    #       re-derived 0.0045 -> 0.05623 against THESE two literature benchmarks (the same
    #       two it was originally fitted to, with the barriers frozen; see
    #       results/validation/hydrolysate_observability_rederivation.md). Methional went
    #       16.2x/12.4x under to 1.30x under / 1.01x, so methional is no longer the worst
    #       row on either lane -- FFT is.
    # The MFT and FFT observability factors were deliberately NOT re-derived: their
    # unconstrained optima are 3.5x and 8.6x ABOVE the physical maximum of 1.0 for a
    # surviving fraction, i.e. the observability layer is saturated and CANNOT explain
    # this residual. The gap is real and is reported rather than absorbed.
    # The ORDERING claim this test is named for (wheat > SPI) is unaffected and asserted
    # above. Contract tolerances are untouched, so both benchmarks stay loud failures.
    spi_worst = max(comparison.ratio for comparison in spi_eval.comparisons)
    wheat_worst = max(comparison.ratio for comparison in wheat_eval.comparisons)
    assert 60.0 <= spi_worst <= 140.0, (
        f"SPI hydrolysate worst ratio is {spi_worst:.3f}x; the post-G1/Wave-H value is "
        "94.617x (FFT). A move outside [60, 140] means something beyond the 2026-08-27 "
        "approved chemistry changes has shifted this lane -- investigate before re-pinning."
    )
    assert 30.0 <= wheat_worst <= 70.0, (
        f"Wheat gluten hydrolysate worst ratio is {wheat_worst:.3f}x; the post-G1/Wave-H "
        "value is 46.8x (FFT). A move outside [30, 70] means something beyond the "
        "2026-08-27 approved chemistry changes has shifted this lane -- investigate "
        "before re-pinning."
    )


def test_fast_acrylamide_extrusion_benchmark_under_predicts_after_sme_desaturation():
    # RENAMED + RE-PINNED 2026-08-27 (cause: SME response desaturation, Wave G4 / audit
    # item 3.1). This test used to assert `ratio <= 1.35` and the benchmark read 1.04x.
    # That agreement was read off a SATURATED surface: src/extrusion.py's temperature
    # offset hit its 40 C cap at 250 kJ/kg, which is exactly this benchmark's SME, so the
    # "model" contributing here was a constant. With the offset replaced by an explicit
    # melt energy balance the offset at these conditions is 18.04 C, the effective
    # temperature falls 170.0 C -> 148.0 C, and on the Knol 2009 lumped Ea of 129 kJ/mol
    # that is a factor of 6.2 -- the whole of the observed move to 6.42x under.
    # The SME model is deliberately NOT bent back and the benchmark's own contract
    # (1.5x / 0.2 dex) is deliberately UNTOUCHED, so this benchmark FAILS honestly in the
    # panel. See the model_change_note in the benchmark JSON.
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "acrylamide_spi_extrusion_130C_ACSRef3.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    comparison = evaluation.comparisons[0]
    assert comparison.compound == "acrylamide"
    assert comparison.predicted_ppb < comparison.measured_ppb, (
        "This benchmark is a known UNDER-prediction since the SME desaturation. If it has "
        f"swung to over-predicting (got {comparison.predicted_ppb:.4f} ppb vs "
        f"{comparison.measured_ppb} ppb), re-derive this expectation rather than relaxing it."
    )
    assert 4.5 <= comparison.ratio <= 9.0, (
        f"Acrylamide extrusion ratio is {comparison.ratio:.2f}x; the post-desaturation value "
        "is 6.42x. A move outside [4.5, 9.0] means something other than the known "
        "2026-08-27 SME recalibration has shifted this lane -- investigate before re-pinning. "
        "Closing this gap needs a calorimetric anchor for the melt heat capacity / retained "
        "fraction, or SME-resolved acrylamide data, NOT a wider contract."
    )


def test_cml_cel_and_furosine_benchmarks_execute_with_expected_directionality():
    age_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "cml_cel_commercial_pbma_Foods2023.json")
    furosine_eval = evaluate_benchmark(ROOT / "data" / "benchmarks" / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json")

    assert age_eval.supported, age_eval.reason
    assert age_eval.coverage == 1.0
    assert furosine_eval.supported, furosine_eval.reason
    assert furosine_eval.coverage == 1.0

    age_by_compound = _comparison_map(age_eval)
    assert age_by_compound["Nε-(Carboxyethyl)lysine (CEL)"].predicted_ppb > age_by_compound["Nε-(Carboxymethyl)lysine (CML)"].predicted_ppb

    peak = predict_furosine(140.0, 20.0, lysine_mM=35.0, reducing_sugar_mM=35.0, protein_type="pea_iso", water_activity=0.55)
    severe = predict_furosine(165.0, 20.0, lysine_mM=35.0, reducing_sugar_mM=35.0, protein_type="pea_iso", water_activity=0.55)
    assert peak > severe

    # RE-PINNED 2026-08-27 (cause: safety unit-collision repair, Wave G2). The
    # benchmark's conc_ppb held 87.0 for an anchor that is 8.7 mg per 100 g
    # protein; on the declared 20% protein basis the true value is 17400 ppb, so
    # the old `<= 2.0` was measuring a 1000x unit error against itself. The
    # predictor is now visibly ~200x low and the benchmark FAILS its own
    # unchanged 2.0x contract in the panel, which is the intended, honest state
    # (see the unit_correction_note in the benchmark JSON). The directionality
    # claim this test is named for is asserted above and is unaffected.
    worst_ratio = max(comparison.ratio for comparison in furosine_eval.comparisons)
    assert 150.0 <= worst_ratio <= 260.0, (
        f"Furosine worst ratio is {worst_ratio:.1f}x; the post-unit-correction value is 201x. "
        "A move outside [150, 260] means something other than the known 2026-08-27 unit "
        "correction has shifted the furosine lane -- investigate before re-pinning. Closing "
        "this gap requires calibrating predict_furosine against real ug/kg data, NOT loosening "
        "the benchmark contract."
    )


def test_resconi_identity_gap_subset_benchmark_executes_on_supported_marker_subset():
    evaluation = evaluate_benchmark(ROOT / "data" / "benchmarks" / "resconi_2023_pbma_beef_identity_benchmark.json")

    assert evaluation.supported, evaluation.reason
    assert evaluation.coverage == 1.0
    assert len(evaluation.comparisons) == 1

    comparison = evaluation.comparisons[0]
    assert comparison.compound == "furfural"
    assert comparison.matched_name is not None
    # RE-PINNED 2026-08-27 (cause: projection budget retune + Boltzmann de-duplication,
    # audit remediation Part 1). Furfural here went 891 ppb -> 3313 ppb against a
    # 1040 ppb reference, i.e. ratio 1.167 -> 3.186. Furfural is the species that gains
    # most from the de-duplication: its own path activity was pinned at 1.0 and its
    # output was pure normalisation against a softmax denominator, so flattening the
    # selectivity hands it a much larger share of the budget. This is the same defect
    # that produced the spurious furfural temperature maximum at ~145-150 C, which is
    # now gone (the curve is monotone from 110 to 190 C). The benchmark's own tolerance
    # is untouched -- it now reports `scale-gap` in the panel summary.
    assert 2.5 <= comparison.ratio <= 4.0, (
        f"Resconi furfural ratio is {comparison.ratio:.3f}x; post-retune value is 3.186x. "
        "A move outside [2.5, 4.0] means something beyond the 2026-08-27 retune has "
        "shifted the carbohydrate-pyrolysis lane -- investigate before re-pinning."
    )