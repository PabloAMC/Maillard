"""Wave G2 safety-layer regression tests (2026-08-27).

One test per fix in the un-audited-module-hunt TIER 1 list: the unit contract,
the [0, 1] risk score, dose-dependent mitigation, pH plumbing and shape, the
molar_ratios unit declaration, flag hygiene, and analyte-context splitting.
"""

import json
import warnings
from pathlib import Path

import pytest

from src import safety
from src.literature_runtime import (
    MOLAR_RATIO_CANONICAL_UNIT,
    _assess_concentration_unit,
    _build_protein_damage_markers_lane,
    _sum_effective_molar_ratios,
    _token_matches,
)
from src.safety import (
    ACRYLAMIDE_ACTION_PPB,
    ACRYLAMIDE_ATTENTION_PPB,
    FUROSINE_ANCHOR_PROTEIN_FRACTION,
    MG_PER_KG_TO_PPB,
    UNITS,
    build_safety_reference_context,
    evaluate_formulation_safety,
    get_safety_reference_entries,
    mg_per_100g_protein_to_ppb,
    predict_acrylamide,
    predict_cel,
    predict_cml,
    predict_furosine,
)


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = ROOT / "data" / "benchmarks"


# --------------------------------------------------------------------------
# Fix 1: the unit collision
# --------------------------------------------------------------------------


def test_units_contract_declares_every_predictor_output_once():
    for key in (
        "predict_acrylamide.acrylamide_ppb",
        "predict_cml",
        "predict_cel",
        "predict_furosine",
        "evaluate_formulation_safety.risk",
        "precursor_concentration",
    ):
        assert key in UNITS, f"{key} has no declared unit"
    # All four analyte predictors must agree on one unit; the collision was two
    # of them being read as mg/kg elsewhere.
    assert (
        UNITS["predict_cml"]
        == UNITS["predict_cel"]
        == UNITS["predict_furosine"]
        == UNITS["predict_acrylamide.acrylamide_ppb"]
        == "ppb (ug/kg food)"
    )
    assert UNITS["evaluate_formulation_safety.risk"] == "dimensionless [0, 1]"


def test_mg_per_100g_protein_conversion_is_explicit_and_documented():
    # 8.7 mg/100 g protein at 20% protein -> 17.4 mg/kg -> 17400 ppb.
    assert FUROSINE_ANCHOR_PROTEIN_FRACTION == pytest.approx(0.20)
    assert mg_per_100g_protein_to_ppb(8.7) == pytest.approx(17400.0)
    assert mg_per_100g_protein_to_ppb(8.7, protein_fraction=1.0) == pytest.approx(87000.0)
    assert 32.0 * MG_PER_KG_TO_PPB == pytest.approx(32000.0)


def test_corrected_benchmarks_hold_true_ppb_and_declare_the_expected_miss():
    age = json.loads((BENCHMARK_DIR / "cml_cel_commercial_pbma_Foods2023.json").read_text())
    furosine = json.loads(
        (BENCHMARK_DIR / "furosine_extrusion_crossover_140C_RamirezJimenez2000.json").read_text()
    )

    assert age["measured_volatiles"]["Nε-(Carboxymethyl)lysine (CML)"]["conc_ppb"] == 32000.0
    assert age["measured_volatiles"]["Nε-(Carboxyethyl)lysine (CEL)"]["conc_ppb"] == 55000.0
    assert furosine["measured_volatiles"]["furosine"]["conc_ppb"] == 17400.0

    for bench in (age, furosine):
        note = bench["unit_correction_note"]
        assert note["date"] == "2026-08-27"
        assert note["status"] == "corrected_expected_to_fail"
        assert note["expected_miss"]
        # Contracts must NOT have been loosened to hide the honest miss.
        assert bench["validation_contract"]["scale_thresholds"]["max_ratio"] <= 2.0

    # The miss magnitudes recorded in the notes must still be the real ones.
    cml = predict_cml(45.0, 45.0, 150.0, 20.0, water_activity=0.45)
    cel = predict_cel(45.0, 45.0, 150.0, 20.0, water_activity=0.45)
    fur = predict_furosine(
        140.0, 20.0, lysine_mM=35.0, reducing_sugar_mM=35.0, protein_type="pea_iso", water_activity=0.55
    )
    assert 1000.0 < 32000.0 / cml < 1400.0
    assert 900.0 < 55000.0 / cel < 1300.0
    assert 150.0 < 17400.0 / fur < 260.0


def test_runtime_damage_anchors_are_converted_to_ppb_and_furosine_ratio_is_not_pinned():
    def _lane(temperature, duration):
        return _build_protein_damage_markers_lane(
            normalized_sugars=["glucose"],
            normalized_amino=["lysine"],
            normalized_additives=[],
            protein_type="pea_iso",
            process_state="extrusion_structured",
            temperature_celsius=temperature,
            time_minutes=duration,
            water_activity=0.45,
            effective_molar_ratios={"D-Glucose": 45.0, "L-Lysine": 45.0},
            pH=6.0,
            concentration_unit="mM",
        )

    mild = _lane(130.0, 10.0)
    severe = _lane(165.0, 30.0)

    assert mild["cml_anchor_ppb"] == pytest.approx(32000.0)
    assert mild["cel_anchor_ppb"] == pytest.approx(55000.0)
    assert mild["furosine_anchor_ppb"] == pytest.approx(17400.0)
    assert mild["damage_marker_units"] == "ppb (ug/kg food)"

    # Pre-fix the furosine ratio was min()-capped at 2.0 for every heated
    # formulation; it must now be a real, varying, sub-unity quantity.
    assert mild["furosine_ratio_vs_crossover"] != pytest.approx(2.0)
    assert severe["furosine_ratio_vs_crossover"] != pytest.approx(2.0)
    assert mild["furosine_ratio_vs_crossover"] != pytest.approx(
        severe["furosine_ratio_vs_crossover"]
    )
    assert 0.0 < mild["furosine_ratio_vs_crossover"] < 1.0


# --------------------------------------------------------------------------
# Fix 2: a real [0, 1] risk score
# --------------------------------------------------------------------------


def test_risk_score_stays_in_unit_interval_across_a_1e9_acrylamide_span():
    scores = []
    for asparagine in (1.0e-4, 1.0e-2, 1.0, 100.0, 1.0e4, 1.0e5):
        score, _flagged = evaluate_formulation_safety(
            {"L-Asparagine": asparagine, "D-Glucose": 50.0}, 180.0, 15.0, 6.5
        )
        assert 0.0 <= score <= 1.0, f"score {score} outside [0, 1] at asn={asparagine}"
        scores.append(score)
    # Monotone in the precursor loading, and the span is actually resolved
    # (the old log-compressed aggregate moved by ~0.9 over 1e9).
    assert scores == sorted(scores)
    assert scores[-1] > scores[0]
    assert scores[-1] == pytest.approx(1.0)


def test_risk_score_is_zero_only_when_nothing_was_assessed():
    score, flagged = evaluate_formulation_safety({"Glycine": 10.0, "D-Glucose": 10.0}, 160.0, 20.0, 7.0)
    assert score == 0.0
    assert "Acrylamide" not in flagged


def test_flag_threshold_is_band_derived_not_the_dead_100_ppb_literal():
    assert ACRYLAMIDE_ATTENTION_PPB == pytest.approx(31.81)
    assert ACRYLAMIDE_ACTION_PPB == pytest.approx(748.0)
    # Both bounds must come from the payload, not from a literal in the module.
    band = safety._reference_band_ppb(analyte="acrylamide")
    assert band["attention_ppb"] == ACRYLAMIDE_ATTENTION_PPB
    assert band["action_ppb"] == ACRYLAMIDE_ACTION_PPB
    assert "foods_2023_pbma_acrylamide_ages" in band["source_ids"]
    # The dead literal and its false regulatory justification are gone.
    source = (ROOT / "src" / "safety.py").read_text()
    assert "aa_ppb < 100.0" not in source
    assert "no such category exists in Reg. (EU)" in source


def test_safety_result_flagged_is_consumed_and_band_relative():
    low = predict_acrylamide(0.05, 1.0, 140.0, 5.0, 6.0)
    high = predict_acrylamide(50.0, 50.0, 180.0, 20.0, 6.5)
    assert low.acrylamide_ppb > 0.0
    assert low.flagged is False
    assert high.flagged is True
    _score, flagged = evaluate_formulation_safety(
        {"L-Asparagine": 50.0, "D-Glucose": 50.0}, 180.0, 20.0, 6.5
    )
    assert "Acrylamide above reference band" in flagged


# --------------------------------------------------------------------------
# Fix 3: dose-dependent mitigation
# --------------------------------------------------------------------------


def test_mitigation_is_dose_dependent_and_monotone():
    def _risk(dose):
        return evaluate_formulation_safety(
            {"L-Asparagine": 1.0, "D-Glucose": 5.0, "L-Cysteine": dose, "Glycine": dose},
            180.0,
            15.0,
            6.5,
        )[0]

    trace = _risk(0.01)
    small = _risk(1.0)
    large = _risk(25.0)
    none = evaluate_formulation_safety(
        {"L-Asparagine": 1.0, "D-Glucose": 5.0}, 180.0, 15.0, 6.5
    )[0]

    assert none >= trace >= small >= large
    assert trace == pytest.approx(none, rel=0.02), "a trace must not earn the full discount"
    assert large < small


def test_mitigation_dose_response_shape_saturates():
    adjustment_trace = safety._resolve_acrylamide_runtime_adjustment(
        {"L-Cysteine": 0.05, "Glycine": 0.05}, {}
    )
    adjustment_full = safety._resolve_acrylamide_runtime_adjustment(
        {"L-Cysteine": 50.0, "Glycine": 50.0}, {}
    )
    assert adjustment_trace["mitigation_fraction"] < 0.05
    assert adjustment_full["mitigation_fraction"] == pytest.approx(safety.MITIGATION_CEILING_BOTH, rel=1e-6)
    # Never above the declared ceiling, whatever the dose.
    huge = safety._resolve_acrylamide_runtime_adjustment({"L-Cysteine": 1.0e6, "Glycine": 1.0e6}, {})
    assert huge["mitigation_fraction"] <= safety.MITIGATION_CEILING_BOTH + 1e-9


# --------------------------------------------------------------------------
# Fix 4: pH plumbing and a bounded pH response
# --------------------------------------------------------------------------


def test_runtime_safety_lane_uses_the_formulation_ph():
    def _lane(ph):
        return _build_protein_damage_markers_lane(
            normalized_sugars=["glucose"],
            normalized_amino=["asparagine"],
            normalized_additives=[],
            protein_type="pea_iso",
            process_state="extrusion_structured",
            temperature_celsius=150.0,
            time_minutes=20.0,
            water_activity=0.45,
            effective_molar_ratios={"D-Glucose": 40.0, "L-Asparagine": 40.0},
            pH=ph,
            concentration_unit="mM",
        )

    acidic = _lane(5.0)
    neutral = _lane(7.0)
    assert acidic["acrylamide_pH_used"] == pytest.approx(5.0)
    assert neutral["acrylamide_pH_used"] == pytest.approx(7.0)
    assert neutral["predicted_acrylamide_ppb"] > acidic["predicted_acrylamide_ppb"]
    # The old hard pin made this comparison tautological.
    assert acidic["predicted_acrylamide_ppb"] != pytest.approx(
        neutral["predicted_acrylamide_ppb"]
    )


def test_acrylamide_ph_response_peaks_and_declines_instead_of_growing_without_bound():
    levels = {
        ph: predict_acrylamide(10.0, 10.0, 150.0, 20.0, ph).acrylamide_ppb
        for ph in (4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0)
    }
    assert levels[5.0] < levels[6.0] < levels[7.0] < levels[8.0]
    assert levels[9.0] < levels[8.0]
    assert levels[11.0] < levels[9.0]
    peak_ph = max(levels, key=levels.get)
    assert 7.0 <= peak_ph <= 9.0
    # pH 6 behaviour is preserved exactly (the decline term is normalised there).
    assert safety._acrylamide_ph_factor(6.0) == pytest.approx(1.0 / (1.0 + 10 ** (8.8 - 6.0)))
    assert safety.ACRYLAMIDE_PH_REFERENCE_ID == "de_vleeschouwer_2006_acrylamide_aqueous"


# --------------------------------------------------------------------------
# Fix 5: molar_ratios unit declaration
# --------------------------------------------------------------------------


def test_undeclared_small_molar_ratios_raise_a_loud_warning_without_rescaling():
    ratios = {"L-Cysteine": 1.0, "D-Ribose": 1.0}
    with pytest.warns(RuntimeWarning, match="NO declared"):
        assessment = _assess_concentration_unit(ratios, context="test")
    assert assessment["looks_like_unitless_ratios"] is True
    assert assessment["assumed_concentration_unit"] == MOLAR_RATIO_CANONICAL_UNIT
    # Values must be untouched.
    assert ratios == {"L-Cysteine": 1.0, "D-Ribose": 1.0}
    assert _sum_effective_molar_ratios(ratios, tokens=["cysteine"]) == pytest.approx(1.0)


def test_declared_concentration_unit_passes_through_and_silences_the_heuristic():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assessment = _assess_concentration_unit(
            {"L-Cysteine": 1.0, "concentration_unit": "mM"}, context="test"
        )
    assert assessment["declared_concentration_unit"] == "mM"
    assert assessment["looks_like_unitless_ratios"] is False


def test_non_mM_declared_unit_warns_and_is_never_silently_rescaled():
    with pytest.warns(RuntimeWarning, match="concentration_unit"):
        _assess_concentration_unit({"L-Cysteine": 0.001, "concentration_unit": "M"}, context="test")
    # Metadata keys never leak into the numeric sum.
    assert _sum_effective_molar_ratios(
        {"L-Cysteine": 0.001, "concentration_unit": "M"}, tokens=["cysteine"]
    ) == pytest.approx(0.001)


def test_large_undeclared_ratios_are_not_warned_about():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assessment = _assess_concentration_unit({"D-Glucose": 45.0, "L-Lysine": 45.0}, context="test")
    assert assessment["looks_like_unitless_ratios"] is False


# --------------------------------------------------------------------------
# Fix 6: flag hygiene
# --------------------------------------------------------------------------


def test_damage_flags_require_actual_thermal_exposure():
    ambient_modifiers = {
        "__extrusion_process__": {
            "moisture_regime": "hme",
            "water_activity": 0.55,
            "total_damage_load": {"furosine_mg_per_kg": 18.0, "lal_mg_per_kg": 68.0},
        }
    }
    score, flagged = evaluate_formulation_safety({"L-Lysine": 45.0}, 25.0, 0.0, 6.5, modifiers=ambient_modifiers)
    assert flagged == []
    assert score == 0.0

    heated = dict(ambient_modifiers)
    _score, heated_flagged = evaluate_formulation_safety(
        {"L-Lysine": 45.0}, 150.0, 20.0, 6.5, modifiers=heated
    )
    assert "Furosine" in heated_flagged
    assert "LAL" in heated_flagged


def test_asnase_is_not_parsed_as_asparagine():
    # safety.py precursor parsing
    score_enzyme, flagged_enzyme = evaluate_formulation_safety(
        {"asnase": 20.0, "D-Glucose": 40.0}, 180.0, 15.0, 6.5
    )
    assert flagged_enzyme == []
    assert score_enzyme == 0.0
    # literature_runtime molar-ratio parsing
    assert _sum_effective_molar_ratios({"asnase": 20.0}, tokens=["asparagine", "asn"]) == 0.0
    assert _sum_effective_molar_ratios({"L-Asparagine": 20.0}, tokens=["asparagine", "asn"]) == pytest.approx(20.0)


def test_word_boundary_matching_covers_the_imp_culture_ferment_traps():
    assert _token_matches("imp", "imp") is True
    assert _token_matches("important supplier note", "imp") is False
    assert _token_matches("imported yeast extract", "imp") is False
    assert _token_matches("yeast fermentation", "ferment") is True
    assert _token_matches("fermented koji paste", "ferment") is True
    assert _token_matches("starter culture", "culture") is True
    assert _token_matches("agriculture residue", "culture") is False
    assert _token_matches("reducing sugars", "sugar") is True


def test_acrylamide_absence_surfaces_as_not_assessed_rather_than_a_clean_zero():
    result = predict_acrylamide(0.0, 40.0, 180.0, 20.0, 6.5)
    assert result.assessed is False
    assert result.acrylamide_ppb == 0.0
    assert "not assessed" in result.description.lower()

    lane = _build_protein_damage_markers_lane(
        normalized_sugars=["glucose"],
        normalized_amino=["lysine"],
        normalized_additives=[],
        protein_type="pea_iso",
        process_state="extrusion_structured",
        temperature_celsius=150.0,
        time_minutes=20.0,
        water_activity=0.45,
        effective_molar_ratios={"D-Glucose": 45.0, "L-Lysine": 45.0},
        pH=6.2,
        concentration_unit="mM",
    )
    assert lane["acrylamide_assessment_status"] == "not_assessed"
    assert lane["acrylamide_not_assessed_reason"] == "no_asparagine_in_formulation"
    assert lane["predicted_acrylamide_ppb"] == 0.0


def test_non_numeric_acrylamide_modifier_does_not_crash():
    # Regression: `ea_mod = v` used to hand a string/dict straight into the
    # Arrhenius term and raise TypeError.
    with pytest.warns(RuntimeWarning, match="non-numeric"):
        score, flagged = evaluate_formulation_safety(
            {"L-Asparagine": 20.0, "D-Glucose": 20.0},
            180.0,
            15.0,
            6.5,
            modifiers={"acrylamide_policy": "aggressive"},
        )
    assert 0.0 <= score <= 1.0
    assert "Acrylamide" in flagged

    score_dict, _flagged = evaluate_formulation_safety(
        {"L-Asparagine": 20.0, "D-Glucose": 20.0},
        180.0,
        15.0,
        6.5,
        modifiers={"__runtime_context__": {"additives": [], "interventions": []}},
    )
    assert 0.0 <= score_dict <= 1.0


def test_sucrose_is_not_counted_as_a_reducing_sugar():
    sucrose_score, sucrose_flagged = evaluate_formulation_safety(
        {"L-Asparagine": 20.0, "Sucrose": 40.0}, 180.0, 15.0, 6.5
    )
    glucose_score, glucose_flagged = evaluate_formulation_safety(
        {"L-Asparagine": 20.0, "D-Glucose": 40.0}, 180.0, 15.0, 6.5
    )
    assert sucrose_flagged == []
    assert sucrose_score == 0.0
    assert "Acrylamide" in glucose_flagged
    assert glucose_score > 0.0


# --------------------------------------------------------------------------
# Fix 7: analyte-context splitting
# --------------------------------------------------------------------------


def test_underscore_composite_analytes_are_visible():
    cml_ids = {entry["id"] for entry in get_safety_reference_entries(analyte="cml", visibility="all")}
    cel_ids = {entry["id"] for entry in get_safety_reference_entries(analyte="cel", visibility="all")}
    furosine_ids = {
        entry["id"] for entry in get_safety_reference_entries(analyte="furosine", visibility="all")
    }
    # "CML_CEL", "cml_cel", "cml_cel_acrylamide" and "furosine_cml" were all
    # invisible while the splitter only handled "+" and "/".
    assert "poulsen_2023_pbma_cml_cel" in cml_ids
    assert "yu_2017_cml_cel_meat_review" in cel_ids
    assert "ma_2024_pbma_extrusion_damage" in cml_ids
    assert "charissou_2007_cookie_cml_furosine" in furosine_ids


def test_reference_context_declares_excluded_evidence():
    context = build_safety_reference_context(analyte="acrylamide")
    assert "excluded_entries" in context
    excluded = context["excluded_entries"]
    # The legacy `category`-schema entries whose coarse category ("ages",
    # "protein_damage_markers") cannot resolve to "acrylamide" are named rather
    # than silently dropped.
    assert len(excluded) >= 1
    excluded_ids = {row["id"] for row in excluded}
    assert "foods_2023_cml_cel_proxy_benchmark" in excluded_ids
    assert all(
        row["reason"] in {"no_analyte_or_category_field", "category_only_entry_not_matched"}
        for row in excluded
    )
    assert str(len(excluded)) in context["exclusion_note"]
    assert "EXCLUDED" in context["exclusion_note"]
    # An entry whose category DOES resolve is included, not excluded.
    assert "acs_ref3_spi_acrylamide_fast_kinetics" not in excluded_ids


def test_category_only_entries_are_matchable_through_the_category_fallback():
    furosine_ids = {
        entry["id"] for entry in get_safety_reference_entries(analyte="furosine", visibility="all")
    }
    assert "krause_2003_furosine_hydrolysis_yields" in furosine_ids
