"""
tests/unit/test_kinetic_core_b4.py

Focused unit tests for the MATRIX / OAV OUTPUT LAYER (Build Wave B4).
Deliberately narrow and fast, in B1/B2/B3's style:

  * THRESHOLD-DOMAIN SELECTION: a threshold is returned only for a matrix that
    has one, and NOTHING is ever borrowed across matrices;
  * INTERVAL-NOT-POINT enforcement for absolutes: an AbsoluteConcentration has
    no float coercion, and OAV comes back as an interval;
  * the MFT DIMER's potency in OAV -- the 15.6x that changes a conclusion;
  * RESIDUAL DECOMPOSITION SUMS: terms + residual == measured, exactly, in
    decades, and the reversible term cannot outrun its evidence ceiling
    unflagged;
  * the 32-of-47 NEGATIVE GATE and the pH-3 adduct gate, enforced;
  * the covalent term is INERT: it contributes exactly zero;
  * the OLD LANE is not reproduced: no protein_source_registry value reaches
    this layer;
  * the HOLD-OUT FIREWALL: the exposed Hong and Brewer literals appear in no
    B4 runtime or fit source file, and the frozen hold-out JSON is not read by
    any of them.

Run with:
  ./scripts/docker_maillard.sh run "pytest tests/unit/test_kinetic_core_b*.py -q"
"""

from __future__ import annotations

import io
import json
import math
import re
import tokenize
from pathlib import Path

import pytest

from src.kinetic_core.matrix_oav import (
    MATRIX_THRESHOLDS,
    SEALED_OR_REFUSED_MATRICES,
    WATER_THRESHOLDS,
    AbsoluteConcentration,
    NoMeasuredThreshold,
    absolute_concentration,
    compare_formulations,
    covalent_channel_state,
    decompose_residual,
    fit_class_binding_constants,
    fit_unsaturation_penalty,
    layer_metadata,
    oav_table,
    odour_activity,
    predict_matrix_shift,
    select_threshold,
    select_threshold_verbose,
)
from src.kinetic_core.parameters_matrix import (
    ADDUCT_POSITIVE_CLASSES,
    COMPOUND_STRUCTURE,
    COVALENT_CEILING,
    HOLDOUT_SEALED_BINDING,
    MATRIX_LOADING,
    NO_ADDUCT_GATE_COMPOUNDS,
    PH_ADDUCT_GATE_BELOW,
    REVERSIBLE_BINDING,
    REVERSIBLE_LOG_SHIFT_CEILING,
    assert_no_dft_matrix,
    assert_no_mocked_protein_source,
    binding_constant_for,
    matrix_registry_metadata,
)

REPO = Path(__file__).resolve().parents[2]


# =========================================================================
# 1. THRESHOLD-DOMAIN SELECTION
# =========================================================================

def test_threshold_selection_returns_a_measured_record_in_its_own_domain():
    record = select_threshold("hexanal", "gelatin_3pct", temperature_c=22.0)
    assert not isinstance(record, NoMeasuredThreshold)
    assert record.matrix == "gelatin_3pct"
    assert record.value_ug_per_l == 58.0
    assert record.thermal_step_after_dosing is False


def test_threshold_is_never_borrowed_across_matrices():
    """A matrix with no measured value gets an explicit state, not a number."""
    result = select_threshold("hexanal", "soy_paste_hong")
    assert isinstance(result, NoMeasuredThreshold)
    assert "GATING HOLD-OUT" in result.reason
    # and the state is not silently a zero or a None that would divide
    with pytest.raises(AttributeError):
        _ = result.value_ug_per_l


def test_every_sealed_matrix_returns_its_governance_reason():
    for matrix in SEALED_OR_REFUSED_MATRICES:
        result = select_threshold("hexanal", matrix)
        assert isinstance(result, NoMeasuredThreshold), matrix
        assert result.reason == SEALED_OR_REFUSED_MATRICES[matrix]


def test_no_hong_or_brewer_or_leksrisompong_threshold_is_tabulated():
    """The hold-out matrices contribute ZERO table entries."""
    tabulated = {r.matrix for r in MATRIX_THRESHOLDS}
    for forbidden in ("soy_paste_hong", "cooked_beef", "caseinate_1pct",
                      "emulsion_10pct_fat", "emulsion_20pct_fat", "soybean_oil",
                      "milk_tian"):
        assert forbidden not in tabulated


def test_a_compound_with_no_threshold_anywhere_is_an_explicit_gap():
    # 2,3-dimethylpyrazine: the corpus has 2,5- and 2,6-DMP and NOT this one.
    result = select_threshold("2_3_dimethylpyrazine", "water")
    assert isinstance(result, NoMeasuredThreshold)


def test_temperature_resolved_domain_picks_the_nearest_measured_point():
    record, diagnostics = select_threshold_verbose(
        "hexanal", "gelatin_3pct", temperature_c=35.0)
    assert record.temperature_c == 37.0
    assert diagnostics["temperature_gap_c"] == pytest.approx(2.0)
    assert diagnostics["borrowed"] is False


def test_multiple_water_sources_report_their_spread_rather_than_averaging():
    _, diagnostics = select_threshold_verbose("hexanal", "water")
    assert diagnostics["n_candidates"] == 2          # Guadagni 4.5 and Xin-recovered 5.0
    assert diagnostics["spread_x"] == pytest.approx(5.0 / 4.5)


def test_every_threshold_record_is_an_input_never_a_target():
    for record in WATER_THRESHOLDS + MATRIX_THRESHOLDS:
        assert record.role == "reference_input"


# =========================================================================
# 2. INTERVAL-NOT-POINT ENFORCEMENT FOR ABSOLUTES
# =========================================================================

def test_absolute_concentration_has_no_float_coercion():
    absolute = absolute_concentration(400.0)
    assert isinstance(absolute, AbsoluteConcentration)
    with pytest.raises(TypeError):
        float(absolute)


def test_absolute_interval_carries_the_measured_reliability_band():
    absolute = absolute_concentration(400.0, via_partition=True)
    # dispersion component: sqrt(23) two-sided
    assert absolute.components["hs_spme_same_sample_dispersion"] == pytest.approx(
        math.log10(math.sqrt(23.0)))
    assert absolute.components["air_water_partition_constant"] == pytest.approx(0.5)
    assert absolute.lo_ug_per_l < 400.0 < absolute.hi_ug_per_l
    # the band is wide, and it is supposed to be
    assert absolute.band_x > 40.0


def test_a_bare_float_handed_to_oav_is_wrapped_not_divided_as_a_point():
    result = odour_activity("MFT", 1.59, "water")
    assert result.oav_lo is not None and result.oav_hi is not None
    assert result.oav_lo < result.oav_point < result.oav_hi
    assert result.concentration is not None
    assert isinstance(result.concentration, AbsoluteConcentration)


def test_oav_without_a_threshold_is_a_state_not_a_number():
    result = odour_activity("hexanal", 400.0, "soy_paste_hong")
    assert result.threshold_state == "no_measured_threshold_for_this_matrix"
    assert result.oav_point is None
    assert result.diagnostics["borrowed_from_another_matrix"] is False


def test_absolute_interval_is_wider_when_it_passed_through_a_partition_step():
    with_partition = absolute_concentration(400.0, via_partition=True)
    without = absolute_concentration(400.0, via_partition=False)
    assert with_partition.halfwidth_decades > without.halfwidth_decades


# =========================================================================
# 3. THE DIMER'S POTENCY IN OAV -- the 15.6x that changes a conclusion
# =========================================================================

def test_dimer_potency_ratio_is_15_6x():
    mft = select_threshold("MFT", "water").value_ug_per_l
    dimer = select_threshold("MFTD", "water").value_ug_per_l
    assert mft / dimer == pytest.approx(15.625)


def test_dimer_oav_matches_the_monomer_at_under_ten_percent_of_the_mass():
    """
    The structural fact the whole thiol-consumption module hangs on: only
    6.5-9.6 % of MFT-equivalents sit in the dimer, but the dimer is 15.6x more
    potent, so its OAV MATCHES the monomer's. Dimerisation is NOT aroma loss.
    """
    mft_ug, dimer_ug = 1.0, 0.08          # 8 % of the mass, in the measured band
    table = oav_table({"MFT": mft_ug, "MFTD": dimer_ug}, matrix="water")
    ratio = table["dimer_over_monomer_OAV"]
    assert ratio > 1.0, "the dimer's OAV must not be reported as a pure loss"
    assert table["dimer_potency_ratio_over_monomer"] == pytest.approx(15.625)


def test_oav_table_tracks_the_dimer_as_its_own_species():
    table = oav_table({"MFT": 1.0, "MFTD": 0.08, "MMFT": 0.4}, matrix="water")
    assert "MFTD" in table["per_species"]
    assert table["per_species"]["MFTD"]["OAV_point"] is not None
    # MMFT has no threshold in the corpus -> an explicit gap, not a zero
    assert table["per_species"]["MMFT"]["OAV_point"] is None
    assert table["n_without_threshold"] == 1


# =========================================================================
# 4. RESIDUAL DECOMPOSITION SUMS
# =========================================================================

def test_residual_decomposition_sums_exactly_in_decades():
    prediction = predict_matrix_shift("hexanal", "gelatin_3pct")
    decomposition = decompose_residual(prediction, measured_ratio=12.888888)
    total = sum(decomposition.per_term_decades.values()) + decomposition.residual_decades
    assert total == pytest.approx(decomposition.measured_decades, abs=1e-12)
    assert decomposition.explained_decades == pytest.approx(
        sum(decomposition.per_term_decades.values()), abs=1e-12)


def test_residual_x_is_the_antilog_of_the_residual_decades():
    prediction = predict_matrix_shift("hexanal", "gelatin_3pct")
    decomposition = decompose_residual(prediction, 12.888888)
    assert decomposition.residual_x == pytest.approx(
        10.0 ** decomposition.residual_decades)
    # and the three factors times the residual reproduce the measurement
    product = decomposition.residual_x
    for value in decomposition.per_term_decades.values():
        product *= 10.0 ** value
    assert product == pytest.approx(12.888888, rel=1e-9)


def test_the_covalent_term_is_inert_it_contributes_exactly_zero():
    for compound in ("hexanal", "t_2_hexenal", "2_methylbutanal"):
        prediction = predict_matrix_shift(compound, "gelatin_3pct")
        assert prediction.terms["covalent_ceiling"]["contribution_to_point_prediction"] == 0.0
        decomposition = decompose_residual(prediction, 10.0)
        assert decomposition.per_term_decades["covalent_ceiling"] == 0.0


def test_reversible_term_cannot_outrun_its_evidence_ceiling_unflagged():
    """
    Amendment 6 ruling 2 caps reversible binding at ~25 % of an observed
    log-shift. A decomposition in which it claims more must FLAG, because a term
    that outruns its own evidence is a bug in the layer, not a discovery.
    """
    prediction = predict_matrix_shift("hexanal", "gelatin_3pct")
    # a deliberately small "measured" shift, so the term claims most of the log
    decomposition = decompose_residual(prediction, measured_ratio=1.4)
    share = (decomposition.per_term_decades["reversible_binding"]
             / decomposition.measured_decades)
    assert share > REVERSIBLE_LOG_SHIFT_CEILING
    assert any("EVIDENCE CEILING" in f for f in decomposition.flags)


def test_a_compound_with_no_term_reports_the_whole_shift_as_residual():
    prediction = predict_matrix_shift("butyric_acid", "soy_paste_hong")
    assert prediction.state == "no_binding_constant_for_class"
    assert prediction.predicted_ratio == 1.0
    decomposition = decompose_residual(prediction, 2000.0)
    assert decomposition.explained_decades == pytest.approx(0.0)
    assert decomposition.residual_decades == pytest.approx(math.log10(2000.0))
    assert any("NO TERM AVAILABLE" in f for f in decomposition.flags)
    assert any("SIGN NOT PREDICTED" in f for f in decomposition.flags)


def test_exactly_1x_is_reported_as_the_absence_of_a_prediction_not_a_flat_one():
    prediction = predict_matrix_shift("4_ethylphenol", "soy_paste_hong")
    decomposition = decompose_residual(prediction, 30.0)
    assert decomposition.sign_predicted == "no_prediction"
    assert decomposition.sign_correct is False


# =========================================================================
# 5. THE NEGATIVE GATE AND THE pH GATE
# =========================================================================

def test_the_negative_gate_has_exactly_32_members():
    assert len(NO_ADDUCT_GATE_COMPOUNDS) == 32
    assert len(set(NO_ADDUCT_GATE_COMPOUNDS)) == 32


def test_named_members_of_the_negative_gate_get_no_covalent_channel():
    for compound in ("4_vinylphenol", "butyric_acid", "dimethyl_disulfide",
                     "2_heptanone", "2_nonanone", "furaneol"):
        state = covalent_channel_state(compound)["state"]
        assert state.startswith("blocked_by_negative_gate"), (compound, state)


def test_gated_classes_reach_compounds_the_panel_did_not_contain():
    """2,3-dimethylpyrazine was never dosed; 2,5-dimethylpyrazine was, and is
    one of the 32. The gate reaches the untested isomer by CLASS."""
    assert covalent_channel_state("2_3_dimethylpyrazine")["state"] == (
        "blocked_by_negative_gate_class")
    assert covalent_channel_state("4_ethylphenol")["state"].startswith(
        "blocked_by_negative_gate")


def test_alkanals_and_thiols_do_get_a_structurally_allowed_channel():
    for compound in ("hexanal", "2_methylbutanal", "t_2_hexenal", "FFT"):
        assert covalent_channel_state(compound)["state"] == "structurally_allowed"
        assert COMPOUND_STRUCTURE[compound].binding_class in ADDUCT_POSITIVE_CLASSES


def test_the_ph_3_adduct_gate_abolishes_the_channel():
    assert covalent_channel_state("hexanal", ph=7.0)["state"] == "structurally_allowed"
    assert covalent_channel_state("hexanal", ph=PH_ADDUCT_GATE_BELOW)["state"] == (
        "blocked_by_pH_gate")
    assert covalent_channel_state("hexanal", ph=2.5)["state"] == "blocked_by_pH_gate"


def test_the_intermediate_ph_band_is_reported_uncertain_not_interpolated():
    state = covalent_channel_state("hexanal", ph=4.5)
    assert state["state"] == "structurally_allowed"
    assert "uncertain rather than interpolated" in state["ph_state"]


def test_the_dmds_source_contradiction_rides_on_every_dmds_row():
    state = covalent_channel_state("dimethyl_disulfide")
    assert "source_contradiction" in state
    assert "Table 2" in state["source_contradiction"]


def test_a_generic_protein_binding_loss_term_is_impossible_by_construction():
    """The 32-of-47 gate exists to falsify exactly that. Verify no more than a
    minority of the structural registry can even receive a covalent channel."""
    allowed = sum(1 for key in COMPOUND_STRUCTURE
                  if covalent_channel_state(key)["state"] == "structurally_allowed")
    assert allowed < len(COMPOUND_STRUCTURE) / 2


# =========================================================================
# 6. THE FIT, AND WHAT IT REFUSES TO DO
# =========================================================================

def test_aldehyde_binding_constants_are_never_pooled_across_methods():
    headspace = binding_constant_for("hexanal", "headspace")
    dialysis = binding_constant_for("hexanal", "dialysis")
    assert headspace is not None and dialysis is not None
    assert headspace.method == "static_headspace_partition"
    assert dialysis.method == "gel_filtration"
    assert headspace.value != dialysis.value
    # and the fitted class constant uses the headspace family only
    classes = fit_class_binding_constants()
    assert classes["n_alkanal"]["method_family"] == "headspace"
    assert classes["n_alkanal"]["k_g_l_per_g"] == pytest.approx(headspace.value)


def test_the_quarantined_alkenal_constant_never_becomes_a_binding_class():
    classes = fit_class_binding_constants()
    assert "alkenal" not in classes


def test_measured_enhancements_survive_the_fit_as_negative_constants():
    """Any matrix model that can only increase thresholds is refuted by two
    independent papers. The layer must be able to emit a ratio below 1."""
    classes = fit_class_binding_constants()
    assert classes["furanone"]["k_g_l_per_g"] < 0
    assert classes["lactone"]["k_g_l_per_g"] < 0
    prediction = predict_matrix_shift("furaneol", "caseinate_1pct")
    assert prediction.predicted_ratio < 1.0


def test_enhancement_is_not_extrapolated_to_another_loading():
    prediction = predict_matrix_shift("furaneol", "soy_paste_hong")
    assert prediction.terms["reversible_binding"]["state"] == (
        "measured_enhancement_extrapolation_refused")
    assert any("MEASURED ENHANCEMENT" in w for w in prediction.warnings)


def test_the_unsaturation_penalty_excludes_the_holdout_beef_observation():
    penalty = fit_unsaturation_penalty()
    assert penalty["n_fit_rows"] == 2
    assert 2.01 not in penalty["observations"].values()
    assert penalty["excluded"]["unsat_penalty_beef"] == 2.01
    assert penalty["penalty_x"] == pytest.approx(math.sqrt(2.81 * 4.95), rel=1e-9)


def test_the_unsaturation_penalty_needs_a_carbonyl_not_just_a_double_bond():
    """4-vinyl phenol is conjugated and has NO carbonyl -> no penalty."""
    assert COMPOUND_STRUCTURE["4_vinylphenol"].alpha_beta_unsaturated_carbonyl is False
    prediction = predict_matrix_shift("4_vinylphenol", "soy_paste_hong")
    assert prediction.terms["alpha_beta_unsaturation"]["factor_x"] == 1.0
    assert COMPOUND_STRUCTURE["t_2_hexenal"].alpha_beta_unsaturated_carbonyl is True


def test_a_declared_assumption_loading_propagates_a_band_and_a_warning():
    prediction = predict_matrix_shift("hexanal", "soy_paste_hong")
    assert prediction.predicted_lo < prediction.predicted_ratio < prediction.predicted_hi
    assert any("DECLARED ASSUMPTION" in w for w in prediction.warnings)
    assert any("EXTRAPOLATION WARNING" in w for w in prediction.warnings)


def test_holdout_binding_constants_are_registered_but_have_no_values():
    assert HOLDOUT_SEALED_BINDING
    for key, reason in HOLDOUT_SEALED_BINDING.items():
        assert "HOLD-OUT" in reason
    carried = {p.key for p in REVERSIBLE_BINDING}
    assert carried.isdisjoint(HOLDOUT_SEALED_BINDING)


def test_ratios_lead_and_unresolved_ratios_say_so():
    comparison = compare_formulations(
        {"MFT": 10.0, "hexanal": 100.0}, {"MFT": 1.0, "hexanal": 90.0}, "A", "B")
    rows = {r["compound"]: r for r in comparison["rows"]}
    # 100/90 = 1.11x sits inside the sqrt(23) = 4.80x same-sample band
    assert rows["hexanal"]["within_reliability_band"] is True
    assert "NOT RESOLVED" in rows["hexanal"]["note"]
    # 10x clears it
    assert rows["MFT"]["within_reliability_band"] is False
    assert rows["MFT"]["ratio_a_over_b"] == pytest.approx(10.0)
    assert rows["MFT"]["direction"] == "higher_in_A"
    assert comparison["n_resolved"] == 1
    assert comparison["reliability_band_x"] == pytest.approx(math.sqrt(23.0))


# =========================================================================
# 7. POLICY: NO DFT, NO MOCKED PROTEIN REGISTRY, NO OLD LANE
# =========================================================================

def test_no_dft_and_no_mocked_protein_source():
    assert_no_dft_matrix()
    assert_no_mocked_protein_source()


def test_no_b4_source_file_references_the_unsourced_protein_registry():
    for relative in B4_SOURCE_FILES:
        text = (REPO / relative).read_text()
        assert "protein_source_registry" not in text or "no_verifiable_source" in text, (
            f"{relative} names the unsourced registry outside a policy note")
        for forbidden in ("meaty_potential_multiplier", "off_note_penalty",
                          "lox_activity_flag", "hydrolysate_observability_bias",
                          "methoxypyrazine_ceiling", "matrix_uncertainty_factor"):
            assert forbidden not in text, f"{relative} reproduces {forbidden}"


def test_b4_does_not_import_the_old_lane():
    import src.kinetic_core.matrix_oav as module
    import src.kinetic_core.parameters_matrix as registry
    for namespace in (vars(module), vars(registry)):
        for name, value in namespace.items():
            module_name = getattr(value, "__module__", "")
            assert module_name not in ("src.matrix_correction", "src.headspace"), (
                f"{name} came from the old lane")


def test_layer_metadata_names_what_it_refuses():
    refused = matrix_registry_metadata()["refused_by_policy"]
    joined = " ".join(refused)
    assert "general matrix correction factor" in joined
    assert "log P" in joined
    assert "protein_source_registry" in joined
    assert "HOLD-OUT" in joined
    assert "DFT" in joined
    assert layer_metadata()["primary_output"].startswith("formulation")


def test_covalent_ceiling_carries_an_unmeasured_ea_rather_than_a_guess():
    assert COVALENT_CEILING.ea_kj_mol is None
    assert COVALENT_CEILING.k2_m_per_s_at_20c == 2.5e-5
    assert COVALENT_CEILING.reversible_share_headspace_timescale >= 0.98
    report = matrix_registry_metadata()["covalent_ceiling"]
    assert "UNMEASURED" in report["ea_status"]


# =========================================================================
# 8. THE HOLD-OUT FIREWALL
# =========================================================================

B4_SOURCE_FILES = (
    "src/kinetic_core/parameters_matrix.py",
    "src/kinetic_core/matrix_oav.py",
    "scripts/generators/generate_kinetic_core_b4_fit.py",
)

#: Literals from hold-out datasets that this wave was EXPOSED to during its
#: instructed reading of k2 and k4b (disclosed in the fit report). None of them
#: may appear in executable code in any B4 runtime or fit file.
#:
#: Deliberately NOT listed, and the reason matters: bare "5.0", "3.0" and "1.5"
#: are too generic to discriminate -- they match ordinary constants and a guard
#: that fires on those trains the reader to ignore it.
HOLDOUT_LITERALS = (
    # Brewer 1995 beef thresholds (ppm) and their computed ratios
    "2.67", "5.87", "7.87", "4.20", "0.47",
    "223", "1304", "1_304", "2623", "2_623", "1400", "6714", "6_714",
    # Brewer/gelatin derived ratios
    "101", "65", "72", "38.5",
    # Hong 2020 aggregate statistics exposed via k4b
    "2035", "2_035", "29.1", "132.5", "277", "1669.6", "1_669.6", "12.6",
    "0.155", "0.668", "1.206", "0.052", "0.070", "24942",
    # Leksrisompong BETs (not fit-eligible)
    "40.8", "44.9", "546", "294", "21.8", "329", "99.5", "1550", "1_550",
    "90.8", "46.4", "56.4", "67.8", "66.9",
    # Barallat-Perez lupin / mucin constants (Module 6 HOLD-OUT)
    "0.0515", "5.15e-2", "0.317", "3.17e-1", "0.0471", "4.71e-2",
    "0.64", "6.4e-1", "2.35", "1.36",
)


def _executable_code(path: Path) -> str:
    """Strip comments, docstrings and string literals; return only code."""
    out = []
    with open(path, "rb") as handle:
        for token in tokenize.tokenize(handle.readline):
            if token.type in (tokenize.COMMENT, tokenize.STRING):
                continue
            out.append(token.string)
    return " ".join(out)


def test_holdout_firewall_no_exposed_literal_reaches_executable_code():
    for relative in B4_SOURCE_FILES:
        code = _executable_code(REPO / relative)
        numbers = set(re.findall(r"\d+\.?\d*(?:[eE][+-]?\d+)?", code))
        for literal in HOLDOUT_LITERALS:
            assert literal not in numbers, (
                f"{relative}: hold-out literal {literal!r} appears in executable code")


def test_holdout_firewall_no_runtime_or_fit_file_reads_the_frozen_values():
    frozen_name = "hong2020_paired_thresholds"
    for relative in B4_SOURCE_FILES:
        # Executable code only: the fit report NAMES the file in prose, in the
        # disclosure section, which is exactly where it should be named.
        assert frozen_name not in _executable_code(REPO / relative), (
            f"{relative} reads the frozen hold-out file; only the scorer may")
    scorer = REPO / "scripts/generators/generate_kinetic_core_b4_holdout.py"
    assert frozen_name in scorer.read_text()


def test_holdout_firewall_no_b4_file_touches_the_frozen_bundles():
    for relative in B4_SOURCE_FILES + (
            "scripts/generators/generate_kinetic_core_b4_holdout.py",):
        text = (REPO / relative).read_text()
        for line in text.splitlines():
            if "data/benchmarks/external_validation" in line:
                assert ("NEVER" in line or "never" in line or "#" in line
                        or line.strip().startswith(('"', "'", "*", "-"))), (
                    f"{relative}: live reference to the frozen bundles: {line!r}")


def test_the_scorer_contains_no_optimiser():
    # Executable code only: the module docstring says the words in order to
    # promise their absence, which is the point.
    code = _executable_code(REPO / "scripts/generators/generate_kinetic_core_b4_holdout.py")
    for forbidden in ("least_squares", "minimize", "curve_fit", "differential_evolution",
                      "scipy"):
        assert forbidden not in code


def test_frozen_predictions_exist_and_predate_the_score():
    path = REPO / "results/validation/kinetic_core_b4_frozen_predictions.json"
    if not path.exists():
        pytest.skip("run generate_kinetic_core_b4_fit.py first")
    frozen = json.loads(path.read_text())
    assert len(frozen["predictions"]) == 10
    assert frozen["pre_registered_expectations"]["expected_failures"]
    # the pre-registration must name the inversion compound as an expected miss
    joined = json.dumps(frozen["pre_registered_expectations"])
    assert "ethyl_4_methylpentanoate" in joined
    assert "SIGN WRONG" in joined


# =========================================================================
# 9. PINNED STRUCTURAL REGRESSIONS
# =========================================================================

def test_pinned_fit_row_reproduction_meynier_hexanal():
    """The layer must reproduce the FIT row it was built from, to 4 digits."""
    prediction = predict_matrix_shift("hexanal", "skim_milk")
    assert prediction.predicted_ratio == pytest.approx(1.390, abs=1e-3)


def test_pinned_soy_hexanal_prediction_is_small_and_says_so():
    """
    The pinned honesty regression: the layer predicts a single-digit shift for
    hexanal in a dense protein paste, because a single-digit shift is all
    reversible binding can produce. Amendment 6 already computed that this
    under-explains a real matrix shift by ~0.75 of the log.
    """
    prediction = predict_matrix_shift("hexanal", "soy_paste_hong")
    assert 1.5 < prediction.predicted_ratio < 5.0
    assert prediction.terms["covalent_ceiling"]["state"] == "structurally_allowed"
    assert "UNMEASURED" in prediction.terms["covalent_ceiling"][
        "process_temperature_report"]


def test_the_two_branched_butanals_are_indistinguishable_and_that_is_recorded():
    two = predict_matrix_shift("2_methylbutanal", "soy_paste_hong")
    three = predict_matrix_shift("3_methylbutanal", "soy_paste_hong")
    assert two.predicted_ratio == three.predicted_ratio
    assert two.state == "predicted_with_chain_length_surrogate"
    assert "Branch position is unmeasured" in (
        two.terms["reversible_binding"]["surrogate"])


def test_matrix_loadings_are_printed_composition_or_declared_assumption():
    for key, loading in MATRIX_LOADING.items():
        assert loading.evidence_class in ("measured_ratio", "declared_assumption")
        if loading.evidence_class == "declared_assumption":
            assert loading.protein_lo_g_per_l < loading.protein_hi_g_per_l, key
