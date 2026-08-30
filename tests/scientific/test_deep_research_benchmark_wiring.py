from pathlib import Path

import json

from src.benchmark_validation import evaluate_benchmark, evaluate_benchmark_payload
from src.safety import predict_furosine


ROOT = Path(__file__).resolve().parents[2]


def _comparison_map(evaluation):
    return {comparison.compound: comparison for comparison in evaluation.comparisons}


def _hydrolysate_payload(benchmark_id: str, protein_source: str, cysteine_mM: float,
                         methionine_mM: float, glycine_mM: float) -> dict:
    """A SYNTHETIC hydrolysate + xylose formulation. No reference values are asserted.

    2026-08-27 (Wave I). The `measured_volatiles` entries below are placeholders that exist
    only so the evaluator emits comparison rows; the test scores model-vs-model, never
    model-vs-any-number.
    """
    return {
        "benchmark_id": benchmark_id,
        "protein_source": protein_source,
        "precursors": {
            f"{'Soy protein isolate' if protein_source == 'soy_isolate' else 'Wheat gluten'} hydrolysate": {
                "concentration_mM": 50.0
            },
            "L-Cysteine": {"concentration_mM": cysteine_mM},
            "L-Methionine": {"concentration_mM": methionine_mM},
            "Glycine": {"concentration_mM": glycine_mM},
            "D-Xylose": {"concentration_mM": 50.0},
        },
        "conditions": {"temp_C": 120, "ph": 6.0, "water_activity": 0.92, "time_min": 30},
        "metadata": {
            "tier": "DIAGNOSTIC",
            "family": "alternative_protein_matrix_scope",
            "execution_path": "free_precursor",
            "notes": (
                "Synthetic hydrolysate formulation for a model-internal ordering check. The "
                "conc_ppb placeholders are NOT reference values and nothing scores against "
                "them."
            ),
        },
        "validation_contract": {
            "scale_thresholds": {"max_ratio": 1e9, "mean_abs_log10_error": 9.0}
        },
        "measured_volatiles": {
            "2-Furfurylthiol (FFT)": {"conc_ppb": 1.0, "uncertainty_pct": 100},
            "2-Methyl-3-furanthiol (MFT)": {"conc_ppb": 1.0, "uncertainty_pct": 100},
            "Methional": {"conc_ppb": 1.0, "uncertainty_pct": 100},
        },
    }


def test_hydrolysate_lane_is_executable_and_ranks_wheat_above_spi():
    """RESTRUCTURED 2026-08-27 (Wave I) — was
    `test_hvp_xylose_benchmarks_are_executable_and_rank_wheat_above_spi`.

    THE OLD TEST SCORED AGAINST TWO FABRICATED BENCHMARKS. It loaded
    `spi_hvp_xylose_120C_PMC9905368.json` and `wheat_gluten_hvp_xylose_120C_PMC9905368.json`
    and pinned worst-case fold errors against their `conc_ppb` values. Those values have no
    possible source: the cited paper (10.1007/s10068-022-01194-w) reacts protein
    hydrolysates with glucose and fructose at pH 7.5 for 90 min, reports only RELATIVE GC-MS
    peak areas, and never mentions 2-furfurylthiol or 2-methyl-3-furanthiol at all. Both
    files are now in `data/benchmarks/quarantined/`.

    The magnitude assertions are therefore DELETED, not re-pinned. Re-pinning them would
    have preserved a number whose only meaning was agreement with an invented measurement,
    and each previous wave's careful re-derivation of those bounds (94.6x, 46.8x, and the
    5.165x before them) was, in hindsight, precision about nothing.

    What is KEPT is the one claim that never depended on those values: the lane is
    executable, and a wheat-gluten hydrolysate carrying more cysteine outranks a soy-isolate
    hydrolysate carrying less, at matched conditions. That is a model-internal ordering,
    asserted here on explicitly SYNTHETIC formulations.

    HONESTY NOTE, and it matters: this ordering now has **no literature support in this
    repository**. The two quarantined files were its only anchor. It is retained as a
    structural regression guard -- it would catch the lane silently losing its
    protein-source discrimination -- and NOT as evidence. Any README or report claim that
    'the wheat >> soy hydrolysate ranking holds' as a validated result is withdrawn.
    """
    spi_eval = evaluate_benchmark_payload(
        _hydrolysate_payload("synthetic_spi_hvp_xylose", "soy_isolate", 3.0, 4.0, 8.0)
    )
    wheat_eval = evaluate_benchmark_payload(
        _hydrolysate_payload("synthetic_wheat_gluten_hvp_xylose", "wheat_gluten", 5.0, 6.0, 10.0)
    )

    assert spi_eval.supported, spi_eval.reason
    assert wheat_eval.supported, wheat_eval.reason
    assert spi_eval.coverage == 1.0
    assert wheat_eval.coverage == 1.0

    spi_mft = _comparison_map(spi_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    wheat_mft = _comparison_map(wheat_eval)["2-Methyl-3-furanthiol (MFT)"].predicted_ppb
    assert spi_mft > 0.0 and wheat_mft > 0.0, "degenerate comparison: a lane predicts zero"
    assert wheat_mft > spi_mft, (
        f"wheat-gluten hydrolysate MFT ({wheat_mft:.4g} ppb) must exceed soy-isolate "
        f"hydrolysate MFT ({spi_mft:.4g} ppb) at matched conditions with more cysteine; "
        "the lane has lost its protein-source discrimination."
    )


def test_quarantined_hvp_xylose_benchmarks_stay_out_of_the_panel():
    """2026-08-27 (Wave I): the quarantine must not silently reverse.

    A fabricated benchmark that reappears in `data/benchmarks/` is back in every headline
    count, because the loader's glob is non-recursive. Assert both halves: gone from the
    panel, still present in quarantine for forensics.
    """
    panel = ROOT / "data" / "benchmarks"
    for name in (
        "spi_hvp_xylose_120C_PMC9905368.json",
        "wheat_gluten_hvp_xylose_120C_PMC9905368.json",
    ):
        assert not (panel / name).exists(), (
            f"{name} is back in the panel; its values have no possible source "
            "(10.1007/s10068-022-01194-w reports relative peak areas for glucose/fructose "
            "at pH 7.5 and never mentions FFT or MFT)."
        )
        quarantined = panel / "quarantined" / name
        assert quarantined.exists(), f"{name} is missing from quarantine"
        payload = json.loads(quarantined.read_text())
        assert payload["source_metadata"]["quarantined"] is True


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
    # RE-PINNED AGAIN 2026-08-27 (Wave M). Band [4.5, 9.0] around 6.42x -> [11.0, 20.0]
    # around 15.39x. Cause: the Wave K/M CONTENT correction of the reference, not a model
    # change -- the prediction is unchanged at 9.7484 ppb. See the benchmark's
    # `content_correction_note`: Ma et al. 2024 (Int J Mol Sci 25:8668, PMC11354377) reports
    # acrylamide ONLY in Figure 2D, whose 130 C bar reads ~150 ug/kg; the 62.62 ppb this file
    # carried appears nowhere in the paper and no derivation from the companion Fu 2023
    # reproduces it. The reference is now a declared `figure_readoff` at 150 ppb with a 20%
    # uncertainty, so the honest miss WIDENS 6.42x -> 15.39x. The band width is unchanged in
    # relative terms (~[0.71x, 1.30x] of the pinned value); it is NOT widened.
    assert 11.0 <= comparison.ratio <= 20.0, (
        f"Acrylamide extrusion ratio is {comparison.ratio:.2f}x; the value after the "
        "2026-08-27 Wave K/M reference correction is 15.39x. A move outside [11.0, 20.0] "
        "means something other than the known SME recalibration and reference correction "
        "has shifted this lane -- investigate before re-pinning. "
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
    # RE-PINNED AGAIN 2026-08-27 (Wave M). Band [2.5, 4.0] around 3.186x -> [3.6, 5.7]
    # around 4.561x. Cause: the Wave K/M CONTENT correction of the reference, not a model
    # change -- the prediction is unchanged at 3261.91 ppb. See the `content_correction_note`
    # on the furfural row: 1040.0 +/- 5% was the mean of only TWO of the paper's three PBMA
    # products (987.41 and 1093.54), silently dropping Impossible Burger at 64.71 -- an
    # undeclared cherry-pick inside a 5% tolerance. The honest 3-product mean is 715.22 with
    # SD 565.85 (79%), so the ratio widens 3.186x -> 4.561x. Band width unchanged in relative
    # terms (~[0.79x, 1.25x] of the pinned value); it is NOT widened.
    assert 3.6 <= comparison.ratio <= 5.7, (
        f"Resconi furfural ratio is {comparison.ratio:.3f}x; the value after the 2026-08-27 "
        "Wave K/M reference correction is 4.561x. A move outside [3.6, 5.7] means something "
        "beyond the 2026-08-27 retune and reference correction has shifted the "
        "carbohydrate-pyrolysis lane -- investigate before re-pinning."
    )