"""FROZEN-WAVE REGRESSION RECORD (labelled 2026-09-03, test audit). The wave generator this file
tests is frozen (scripts/generators/WAVES.md); these tests fail only if the frozen report, the
network or the parameter tables change. They are the contract of a finished wave, not live checks
of new behaviour.


tests/unit/test_kinetic_core_b3.py

Focused unit tests for the ACRYLAMIDE / SAFETY MODULE (Build Wave B3).
Deliberately narrow and fast, in B1's and B2's style:

  * carbon / nitrogen / sulfur balance, at import and per reaction;
  * conservation of all three elements through an integration;
  * THE ELIMINATION-EXISTS PROPERTY: acrylamide has real sinks, its trajectory
    has an interior maximum, and it is NOT linear in time -- the specific
    defect that made the old FAST lane ~40x under-responsive to time;
  * competition is mass action, not a multiplier: the effect scales with the
    competitor's concentration and vanishes as it goes to zero;
  * the low-confidence FRUCTOSE reporting: no fructose-specific parameter
    exists, and the hold-out scorer reports an interval rather than a value;
  * the fabricated 129 kJ/mol barrier cannot enter;
  * the hold-out firewall: the declared hold-out literals appear in no runtime
    or fit source file;
  * a pinned structural regression at Claeys' condition.

Run with:
  ./scripts/docker_maillard.sh run "pytest tests/unit/test_kinetic_core_b*.py -q"
"""

from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import operative_acrylamide_parameters
from src.kinetic_core.network import BALANCED_ELEMENTS, Reaction
from src.kinetic_core.acrylamide import (
    ACRYLAMIDE_REACTIONS,
    ACRYLAMIDE_SINK_REACTIONS,
    ACRYLAMIDE_SOURCE_REACTIONS,
    FULL_ACRYLAMIDE_REACTIONS,
    OUT_OF_SCOPE,
    acrylamide_flux_budget,
    apparent_activation_energy,
    apparent_lumped_constants,
    describe_acrylamide,
    integrate_acrylamide,
    validate_acrylamide_balance,
)
from src.kinetic_core.parameters_acrylamide import (
    ASN_SCHIFF_AMADORI_SPLIT,
    CROSS_LAB_CONFLICTS,
    DELIBERATE_UNDERFITS,
    EA_E2_MEASURED_KJ_MOL,
    FITTED_ACRYLAMIDE_BOUNDS_LOG10K,
    FITTED_ACRYLAMIDE_EA_KEYS,
    FITTED_ACRYLAMIDE_KEYS,
    HOLDOUT_EXPOSURE_DISCLOSURE,
    MEASURED_ACRYLAMIDE,
    REFUSED_PARAMETERS,
    T_REF_A_K,
    AcrylamideParameter,
    acrylamide_placeholders,
    acrylamide_registry_metadata,
    assert_no_dft_acrylamide,
    assert_no_fabricated_barrier,
    per_molar_to_per_mmol,
    with_fitted_acrylamide,
)
from src.kinetic_core.species_acrylamide import (
    ACRYLAMIDE_BY_KEY,
    ACRYLAMIDE_INDEX,
    ACRYLAMIDE_SPECIES,
    ACRYLAMIDE_STATE,
    ACRYLAMIDE_TERMINAL_POOLS,
    COMPETITOR_KEYS,
    acrylamide_ppb,
    ppb_to_mmol_per_litre,
    total_element_acrylamide,
)
from src.kinetic_core.species_sulfur import SULFUR_STATE

ROOT = Path(__file__).resolve().parents[2]
CELSIUS = 273.15

B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

#: A parameter set that is NOT the fit. Every fitted step gets the same neutral
#: value, so these tests exercise STRUCTURE and never the fit's numbers: nothing
#: here breaks when the fit is re-run, and nothing here can quietly become a
#: second, undeclared copy of the fitted values.
NEUTRAL_LOG10K = {key: -2.0 for key in FITTED_ACRYLAMIDE_KEYS}
NEUTRAL_EA = {key: 120.0 for key in FITTED_ACRYLAMIDE_EA_KEYS}


@pytest.fixture(scope="module")
def parameters():
    return operative_acrylamide_parameters(B1_FITTED, NEUTRAL_LOG10K, NEUTRAL_EA)


# ---------------------------------------------------------------------------
# Stoichiometry and the three-element balance
# ---------------------------------------------------------------------------


def test_every_reaction_balances_all_three_elements():
    validate_acrylamide_balance()
    for reaction in FULL_ACRYLAMIDE_REACTIONS:
        for element in BALANCED_ELEMENTS:
            left, right = reaction.atom_balance(element, ACRYLAMIDE_BY_KEY)
            assert left == right, f"{reaction.key}: {element} {left} -> {right}"


def test_a_deliberately_unbalanced_step_is_refused():
    """A balance checker that cannot fail is not a balance checker."""
    planted = Reaction("planted", {"Asn": 1}, {"ACR": 1}, None, "")
    with pytest.raises(ValueError, match="does not balance"):
        validate_acrylamide_balance([planted])


def test_the_cysteine_adduct_carries_its_sulfur():
    """The one step where the sulfur element is a real constraint."""
    adduct = next(r for r in ACRYLAMIDE_REACTIONS if r.key == "a_acr_cys")
    left, right = adduct.atom_balance("sulfur", ACRYLAMIDE_BY_KEY)
    assert left == right == 1


def test_sulfur_state_is_a_prefix_of_the_acrylamide_state():
    assert ACRYLAMIDE_STATE[: len(SULFUR_STATE)] == SULFUR_STATE
    assert len(ACRYLAMIDE_STATE) == len(SULFUR_STATE) + len(ACRYLAMIDE_SPECIES)
    # every index B1 and B2 published still resolves to the same species
    assert ACRYLAMIDE_INDEX["Glc"] == 0
    assert ACRYLAMIDE_INDEX["Cys"] < ACRYLAMIDE_INDEX["Asn"]


def test_the_cysteine_pool_is_shared_with_b2_and_not_redeclared():
    """The whole reason the acrylamide block sits on top of the sulfur block."""
    keys = [s.key for s in ACRYLAMIDE_SPECIES]
    assert "Cys" not in keys, (
        "a private cysteine would let B2 and B3 each spend the same pool"
    )
    assert "Cys" in COMPETITOR_KEYS
    assert any("Cys" in r.reactants for r in ACRYLAMIDE_REACTIONS)


def test_every_acrylamide_intermediate_has_formation_and_consumption():
    formed, consumed = set(), set()
    for reaction in FULL_ACRYLAMIDE_REACTIONS:
        formed.update(reaction.products)
        consumed.update(reaction.reactants)
    for species in ACRYLAMIDE_SPECIES:
        if species.role == "reactant":
            assert species.key in consumed, f"{species.key} is never consumed"
            continue
        if species.key in ACRYLAMIDE_TERMINAL_POOLS:
            assert species.key not in consumed
            continue
        assert species.key in formed, f"{species.key} is never formed"
        assert species.key in consumed, f"{species.key} is never consumed"


# ---------------------------------------------------------------------------
# Conservation through an integration
# ---------------------------------------------------------------------------


def test_all_three_elements_are_conserved(parameters):
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS,
        {"Asn": 10.0, "Glc": 10.0, "Cys": 10.0, "Lys": 5.0},
        np.linspace(0.0, 45.0, 21), water_activity=1.0,
    )
    for element in BALANCED_ELEMENTS:
        closure = run.element_closure(element)
        assert closure["max_relative_drift"] < 1e-8, (element, closure)


def test_no_negative_concentrations(parameters):
    run = integrate_acrylamide(
        parameters, 200.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0},
        np.linspace(0.0, 45.0, 21), water_activity=1.0,
    )
    assert run.concentrations.min() >= 0.0


# ---------------------------------------------------------------------------
# THE ELIMINATION-EXISTS PROPERTY
# ---------------------------------------------------------------------------


def test_acrylamide_has_named_sinks_and_sources():
    assert ACRYLAMIDE_SOURCE_REACTIONS == ("a_int1_acr",)
    assert "a_acr_dp" in ACRYLAMIDE_SINK_REACTIONS
    assert "a_acr_cys" in ACRYLAMIDE_SINK_REACTIONS
    assert len(ACRYLAMIDE_SINK_REACTIONS) >= 3


def test_fed_acrylamide_decays(parameters):
    """The most direct statement of the property: a fed pool must go away."""
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS, {"ACR": 1.0},
        np.array([0.0, 30.0]), water_activity=1.0,
    )
    assert run.final("ACR") < 1.0


def test_the_acrylamide_trajectory_has_an_interior_maximum(parameters):
    """
    THE defect this module removes. A formation term with no sink accumulates
    monotonically; Claeys measures a peak at 40-50 min at 160 C. The peak must
    be a consequence of the network, so it is asserted structurally: the
    maximum is strictly inside the window, not at its end.
    """
    times = np.linspace(0.0, 240.0, 121)
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0}, times,
        water_activity=1.0,
    )
    series = run.series("ACR")
    peak = int(np.argmax(series))
    assert 0 < peak < len(times) - 1, "acrylamide never turns over"
    assert series[-1] < series[peak]


def test_the_sinks_flatten_the_response_to_time(parameters):
    """
    The quantitative form of the same property, and the one that speaks to the
    measured defect. The old FAST lane was ~40x UNDER-responsive to time
    because it had no sink: a sinkless pool grows for as long as its precursors
    last. The test compares the SAME network with its acrylamide sinks on and
    off, so it asserts the effect of the sinks rather than an absolute number
    that would depend on the neutral parameter values.
    """
    log10k = dict(NEUTRAL_LOG10K)
    for key in ("k_acr_dp", "k_acr_gln", "k_acr_lys", "k_acr_ala"):
        log10k[key] = FITTED_ACRYLAMIDE_BOUNDS_LOG10K[key][0]
    sinkless = operative_acrylamide_parameters(B1_FITTED, log10k, NEUTRAL_EA)

    def growth(params):
        def at(minutes):
            run = integrate_acrylamide(
                params, 160.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0},
                np.array([0.0, minutes]), water_activity=1.0,
            )
            return run.final("ACR")

        return at(80.0) / at(10.0)

    assert growth(parameters) < growth(sinkless), (
        "the acrylamide sinks must flatten the response to bake time"
    )


def test_removing_the_elimination_channels_removes_the_maximum():
    """
    The control for the two tests above: with every sink turned off, the same
    network accumulates monotonically. This is what the old lane did.
    """
    log10k = dict(NEUTRAL_LOG10K)
    for key in ("k_acr_dp", "k_acr_gln", "k_acr_lys", "k_acr_ala"):
        log10k[key] = FITTED_ACRYLAMIDE_BOUNDS_LOG10K[key][0]
    parameters = operative_acrylamide_parameters(B1_FITTED, log10k, NEUTRAL_EA)
    times = np.linspace(0.0, 240.0, 121)
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0}, times,
        water_activity=1.0,
    )
    series = run.series("ACR")
    assert np.all(np.diff(series) >= -1e-12), (
        "with no sink the pool must be monotone; if it is not, some other step "
        "is quietly consuming acrylamide and the sink list is incomplete"
    )
    assert int(np.argmax(series)) == len(times) - 1


# ---------------------------------------------------------------------------
# COMPETITION IS MASS ACTION, NOT A MULTIPLIER
# ---------------------------------------------------------------------------


def test_there_is_no_competition_multiplier_anywhere():
    for key, parameter in MEASURED_ACRYLAMIDE.items():
        assert parameter.order in (1, 2)
        assert parameter.unit in ("1/min", "L/(mmol*min)")
    for key in FITTED_ACRYLAMIDE_KEYS:
        assert key.startswith("k_"), (
            f"{key} is not a rate constant; the registry holds rate constants "
            f"on balanced reactions and nothing else"
        )
    metadata = acrylamide_registry_metadata(dict(MEASURED_ACRYLAMIDE))
    assert metadata["competition_is_a_multiplier"] is False


def test_a_competitor_lowers_acrylamide_and_the_effect_scales_with_its_dose(parameters):
    def peak(competitor_mmol_l):
        initial = {"Asn": 10.0, "Glc": 10.0}
        if competitor_mmol_l:
            initial["Lys"] = competitor_mmol_l
        run = integrate_acrylamide(
            parameters, 160.0 + CELSIUS, initial, np.linspace(0.0, 60.0, 31),
            water_activity=1.0,
        )
        return float(np.max(run.series("ACR")))

    none, some, more = peak(0.0), peak(5.0), peak(20.0)
    assert some < none, "a competitor must lower acrylamide"
    assert more < some, "and the effect must grow with the dose"
    # and it must vanish continuously, which a flag or a fixed factor cannot do
    assert peak(1e-9) == pytest.approx(none, rel=1e-6)


def test_the_competition_channels_can_be_switched_off_by_the_data():
    """Bounds must ALLOW ~zero, so the fit may reject a channel."""
    for key in ("k_gln_glc", "k_lys_glc", "k_ala_glc",
                "k_acr_gln", "k_acr_lys", "k_acr_ala"):
        low, _ = FITTED_ACRYLAMIDE_BOUNDS_LOG10K[key]
        # at a 10 mmol/L competitor loading this is < 1e-7 per minute
        assert (10.0 ** low) * 10.0 < 1e-7


def test_the_cysteine_competitor_has_no_fitted_parameter():
    """
    Claeys' cysteine row is the panel's one parameter-free prediction, and that
    is only true if BOTH of cysteine's channels are measured.
    """
    assert "k_acr_cys" in MEASURED_ACRYLAMIDE
    assert "k_cys_glc" in MEASURED_ACRYLAMIDE
    assert "k_cys_sink" in MEASURED_ACRYLAMIDE
    for key in ("k_acr_cys", "k_cys_glc", "k_cys_sink"):
        assert key not in FITTED_ACRYLAMIDE_KEYS


# ---------------------------------------------------------------------------
# The apparent-constant extractor
# ---------------------------------------------------------------------------


def test_the_apparent_elimination_constant_is_exact_on_a_first_order_system():
    """
    The recipe that lets a mass-action model be compared with a published
    first-order constant must be EXACT when the system really is first order.
    Fed acrylamide with no competitor decays by one first-order channel, so the
    extracted constant must equal that channel's rate constant.
    """
    log10k = dict(NEUTRAL_LOG10K)
    log10k["k_acr_dp"] = math.log10(0.037)
    parameters = operative_acrylamide_parameters(B1_FITTED, log10k, NEUTRAL_EA)
    observed = apparent_lumped_constants(
        parameters, T_REF_A_K, {"ACR": 1.0}, 30.0, n_points=241,
        water_activity=1.0,
    )
    assert observed["k_E_app_per_min"] == pytest.approx(0.037, rel=2e-4)


def test_the_apparent_activation_energy_recovers_a_planted_barrier():
    ea = apparent_activation_energy(1.0, 400.0, math.exp(1.0), 450.0)
    expected = 8.314462618e-3 * 1.0 / (1.0 / 400.0 - 1.0 / 450.0)
    assert ea == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# Units -- the place a factor of 1000 hides
# ---------------------------------------------------------------------------


def test_the_per_molar_to_per_millimolar_conversion_is_applied_once():
    assert per_molar_to_per_mmol(49.36) == pytest.approx(0.04936)
    assert MEASURED_ACRYLAMIDE["k_acr_cys"].k_ref == pytest.approx(0.04936)
    assert MEASURED_ACRYLAMIDE["k_asn_glc"].k_ref == pytest.approx(1.70e-3)
    # and a second-order constant must carry the millimolar unit
    assert MEASURED_ACRYLAMIDE["k_acr_cys"].unit == "L/(mmol*min)"


def test_the_ppb_conversion_round_trips_and_matches_a_hand_check():
    assert ppb_to_mmol_per_litre(acrylamide_ppb(0.037)) == pytest.approx(0.037)
    # 1 mmol/L of acrylamide is 71.08 mg/L = 71 080 000 ug/kg at unit density
    assert acrylamide_ppb(1.0) == pytest.approx(71.08e3 * 1000.0 / 1000.0 * 1.0)


# ---------------------------------------------------------------------------
# Policy: no DFT, no fabricated barrier, no silent defaults
# ---------------------------------------------------------------------------


def test_no_dft():
    assert_no_dft_acrylamide()
    assert_no_dft_acrylamide(dict(MEASURED_ACRYLAMIDE))


def test_the_fabricated_129_barrier_cannot_enter():
    assert_no_fabricated_barrier()
    planted = dict(MEASURED_ACRYLAMIDE)
    planted["k_planted"] = AcrylamideParameter(
        key="k_planted", transformation="planted", k_ref=1.0, ea_kj_mol=129.0,
        unit="1/min", order=1, evidence_class="measured_rate",
        source_anchor="planted", dossier_anchor="planted", conditions="planted",
        ph_of_measurement=6.0, aw_of_measurement=1.0, t_ref_k=T_REF_A_K,
        temperature_range_c=(120.0, 200.0), rate_transfer="not_licensed",
    )
    with pytest.raises(ValueError, match="129"):
        assert_no_fabricated_barrier(planted)


def test_the_refusal_list_names_the_fabrications_and_their_locations():
    quantities = " ".join(r["quantity"] for r in REFUSED_PARAMETERS)
    assert "activation energy" in quantities
    assert "_acrylamide_ph_factor" in quantities
    for row in REFUSED_PARAMETERS:
        assert row["verdict"].startswith("REFUSE")
        assert row["why"] and row["dossier_anchor"]


def test_a_populated_parameter_without_a_barrier_is_refused():
    with pytest.raises(ValueError, match="activation energy"):
        AcrylamideParameter(
            key="k_bad", transformation="bad", k_ref=1.0, ea_kj_mol=None,
            unit="1/min", order=1, evidence_class="measured_rate",
            source_anchor="", dossier_anchor="", conditions="",
            ph_of_measurement=None, aw_of_measurement=None, t_ref_k=T_REF_A_K,
            temperature_range_c=(120.0, 200.0), rate_transfer="not_licensed",
        )


def test_unpopulated_parameters_refuse_to_integrate():
    from src.kinetic_core import operative_parameters

    parameters = dict(operative_parameters(B1_FITTED))
    parameters.update(MEASURED_ACRYLAMIDE)
    parameters.update(acrylamide_placeholders())
    with pytest.raises(ValueError, match="unpopulated"):
        integrate_acrylamide(
            parameters, 433.15, {"Asn": 10.0, "Glc": 10.0}, np.array([0.0, 1.0])
        )


def test_the_scavenging_barrier_is_measured_and_cannot_be_overridden():
    populated = with_fitted_acrylamide(NEUTRAL_LOG10K, NEUTRAL_EA)
    for key in ("k_acr_gln", "k_acr_lys", "k_acr_ala"):
        assert populated[key].ea_kj_mol == EA_E2_MEASURED_KJ_MOL
    assert EA_E2_MEASURED_KJ_MOL != NEUTRAL_EA["Ea_acr_dp"]


def test_every_parameter_carries_its_provenance_including_water_activity():
    for key, parameter in MEASURED_ACRYLAMIDE.items():
        assert parameter.source_anchor and parameter.dossier_anchor
        assert parameter.conditions
        assert parameter.rate_transfer
        assert parameter.t_ref_k == T_REF_A_K
        assert hasattr(parameter, "aw_of_measurement")
        assert parameter.evidence_class == "measured_rate"


def test_the_water_activity_gap_is_reported_on_every_run(parameters):
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS, {"Asn": 42.0, "Glc": 42.0},
        np.array([0.0, 1.0]), water_activity=0.35,
    )
    warnings = " ".join(run.metadata["extrapolation_warnings"])
    assert "a_w" in warnings
    assert run.metadata["aw_term_present"] is False


def test_there_is_no_ph_term():
    metadata = acrylamide_registry_metadata(dict(MEASURED_ACRYLAMIDE))
    assert metadata["ph_term_present"] is False
    assert "_acrylamide_ph_factor" in metadata["ph_term_absent_because"]


def test_the_schiff_split_is_inert_and_says_so():
    """
    The pinned Schiff/Amadori ratio is reused from B1's GLYCINE fit. That is
    only defensible while the Schiff base has exactly one fate, so the claim is
    tested rather than asserted in prose.
    """
    assert ASN_SCHIFF_AMADORI_SPLIT["rate_transfer"] == "not_licensed"
    sinks = [r for r in ACRYLAMIDE_REACTIONS if "SBA" in r.reactants]
    assert len(sinks) == 1 and sinks[0].products == {"INT1": 1}


def test_the_deliberate_underfit_is_recorded_before_the_fit():
    assert DELIBERATE_UNDERFITS
    row = DELIBERATE_UNDERFITS[0]
    assert "glutamine" in row["row"].lower()
    assert "hold-out" in row["why_no_promotion_TERM_IS_ADDED"].lower()


def test_the_cross_lab_conflicts_are_carried_not_averaged():
    quantities = [r["quantity"] for r in CROSS_LAB_CONFLICTS]
    assert any("elimination activation energy" in q.lower() for q in quantities)
    assert any("formation activation energy" in q.lower() for q in quantities)
    for row in CROSS_LAB_CONFLICTS:
        assert len(row["values"]) >= 3
        assert row["verdict"]


def test_out_of_scope_lanes_name_what_they_strand():
    assert len(OUT_OF_SCOPE) >= 3
    for row in OUT_OF_SCOPE:
        assert row["lane"] and row["what_is_stranded"] and row["why"]
    lanes = " ".join(r["lane"] for r in OUT_OF_SCOPE)
    assert "sucrose" in lanes


# ---------------------------------------------------------------------------
# LOW-CONFIDENCE FRUCTOSE REPORTING
# ---------------------------------------------------------------------------

HOLDOUT_SCRIPT = ROOT / "scripts/generators/generate_kinetic_core_b3_holdout.py"


def _load_holdout_module():
    spec = importlib.util.spec_from_file_location("b3_holdout", HOLDOUT_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_no_fructose_specific_parameter_exists():
    """
    The declaration holds fructose out as the SUGAR-TRANSFER test. A fitted
    fructose constant would make the hold-out unfalsifiable, so there must not
    be one -- fructose may reach acrylamide only through the trunk's measured
    isomerisation.
    """
    keys = set(MEASURED_ACRYLAMIDE) | set(FITTED_ACRYLAMIDE_KEYS)
    for key in keys:
        assert "fru" not in key.lower(), key
    for reaction in ACRYLAMIDE_REACTIONS:
        assert "Fru" not in reaction.reactants, (
            f"{reaction.key} consumes fructose directly; the sugar-transfer "
            f"hold-out requires that fructose reach acrylamide only through "
            f"the trunk's measured isomerisation"
        )


def test_the_fructose_prediction_is_reported_as_a_wide_interval(parameters):
    """
    'Low confidence' has to be a number, not an adjective. The Monte-Carlo band
    over the PUBLISHED parameter uncertainties must be wide enough that a point
    value would be misleading.
    """
    module = _load_holdout_module()

    def peak(params):
        return module.observe(
            params, {"Asn": 100.0, "Fru": 100.0}, 160.0, 0.92, n_points=31
        )["peak_acrylamide_ppb"]

    band = module.interval(peak, parameters, n_draws=12, seed=7)
    assert band["n_draws"] >= 8
    assert band["high"] > band["low"] > 0.0
    assert band["spread_factor"] > 2.0, (
        "a fructose prediction with a narrow band would be claiming a "
        "confidence the corpus does not support"
    )



def test_the_exposure_disclosure_is_carried_in_the_runtime_registry():
    for field in ("what_was_seen", "what_was_not_seen",
                  "what_was_done_about_it", "why_it_is_recorded"):
        assert HOLDOUT_EXPOSURE_DISCLOSURE[field]
    assert "fructose" in HOLDOUT_EXPOSURE_DISCLOSURE["what_was_seen"].lower()


# ---------------------------------------------------------------------------
# THE HOLD-OUT FIREWALL
# ---------------------------------------------------------------------------

#: Numeric values of declared Module 3 HOLD-OUT rows. None may appear as a
#: NUMBER in the runtime modules or in the fit script. They are allowed in the
#: hold-out scorer, which is where they belong.
#:
#: WHY "AS A NUMBER" AND NOT "ANYWHERE IN THE FILE": a hold-out value is
#: dangerous when it is USED -- as a target, a bound, an initialisation or a
#: constant. It is not dangerous when a docstring explains why the row is held
#: out; forbidding that would mean the module could not document its own
#: firewall. The stripper removes docstrings, string literals and comments and
#: searches only what is left, which is executable code.
HOLDOUT_NUMBERS = (
    # De Vleeschouwer 2009 I, the FRUCTOSE column
    "7.40", "9.48", "13.62", "4.08", "122.0", "66.5", "108.3", "149.1", "87.7",
    # De Vleeschouwer 2009 I, the SUCROSE column
    "3.56", "7.29", "96.2", "108.9", "328.2", "180.7", "109.4", "140.6",
    # De Vleeschouwer 2009 II, the GLUTAMINE column
    "8.05", "124.1", "92.4", "35.9",
    # inventory sec. B5.5, the glutamine promotion percentages
    "267", "322",
    # DELIBERATELY NOT LISTED, and the reason matters. Knol 2010's Table 2 has
    # no transcribed values to guard. Knol 2009's band (9.3e3-2.6e4) and the
    # sucrose k_HY / k_I (0.47, 0.50) are too generic to discriminate: "0.47"
    # matches a dozen ordinary constants and a guard that fires on those trains
    # the reader to ignore it. Note also that DV-I sucrose's Ea_I is 113.2,
    # which is NUMERICALLY IDENTICAL to the FIT glucose column's Ea_E and is
    # therefore un-guardable without banning a fit value.
)

RUNTIME_FILES = (
    "src/kinetic_core/acrylamide.py",
    "src/kinetic_core/parameters_acrylamide.py",
    "src/kinetic_core/species_acrylamide.py",
)
FIT_FILE = "scripts/generators/generate_kinetic_core_b3_fit.py"


from tests.support import executable_code, strip_prose  # noqa: E402


@pytest.mark.parametrize("relative", RUNTIME_FILES + (FIT_FILE,))
def test_no_holdout_value_reaches_the_runtime_or_the_fit(relative):
    code = executable_code(ROOT / relative)
    for literal in HOLDOUT_NUMBERS:
        pattern = r"(?<![\w.])" + re.escape(literal) + r"(?![\w.])"
        assert re.search(pattern, code) is None, (
            f"{relative} uses the hold-out value {literal} as a NUMBER in "
            f"executable code. That is a firewall breach."
        )


def test_the_firewall_stripper_actually_finds_a_planted_breach():
    """A firewall test that cannot fail is not a firewall test."""
    planted = 'x = "the note says 8.05"\ny = 8.05  # a breach\n'
    code = strip_prose(planted)
    assert re.search(r"(?<![\w.])8\.05(?![\w.])", code) is not None
    assert code.count("8.05") == 1, "the string literal must have been stripped"


def test_the_fit_script_reads_no_holdout_file():
    code = strip_prose((ROOT / FIT_FILE).read_text())
    assert "external_validation" not in code
    assert "mp_holdout" not in code



# ---------------------------------------------------------------------------
# PINNED REGRESSION at the Claeys condition
# ---------------------------------------------------------------------------


def test_pinned_claeys_condition_regression(parameters):
    """
    A pinned STRUCTURAL regression at Claeys' control condition (10 mmol/L
    asparagine + 10 mmol/L glucose, 160 C, aqueous pH 6).

    It pins the NEUTRAL-parameter network's behaviour, not the fit's numbers,
    so it does not have to be rewritten every time the fit is re-run -- and so
    it cannot become a second, undeclared copy of the fitted values. What it
    asserts is the set of qualitative facts the corpus is unambiguous about:

      * acrylamide is formed at all, and in ppb rather than in percent;
      * it goes through a maximum inside the heating window;
      * asparagine is consumed by BOTH of its measured fates;
      * the yield is a small fraction of the asparagine charged (Claeys' own
        derived plateau is 2886 ppb from 10 mmol/L, i.e. 0.4 mol %);
      * all three elements close.
    """
    times = np.linspace(0.0, 240.0, 121)
    run = integrate_acrylamide(
        parameters, 160.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0}, times,
        water_activity=1.0,
    )
    assert run.peak_acrylamide_ppb() > 0.0
    assert 0.0 < run.molar_yield_on_asparagine() < 0.5, (
        "a yield above 50 mol % would mean the Int1 partition is missing"
    )
    assert 0.0 < run.time_of_peak_min() < 240.0
    assert run.final("Asn") < 10.0
    assert run.final("Asp") > 0.0, "the deamidation lane must run"
    for element in BALANCED_ELEMENTS:
        closure = run.element_closure(element)
        # this pot contains no sulfur at all, so the SULFUR invariant is
        # 0 -> 0 and its RELATIVE drift is undefined. The absolute form is the
        # one that is meaningful for every element.
        initial = closure[f"initial_{element}_mmol_L"]
        assert closure["max_abs_drift_mmol_L"] < 1e-8 * max(initial, 1.0), (
            element, closure
        )


def test_the_flux_budget_accounts_for_every_reaction(parameters):
    budget = acrylamide_flux_budget(
        parameters, 160.0 + CELSIUS, {"Asn": 10.0, "Glc": 10.0, "Cys": 10.0},
        45.0, n_points=41,
    )
    assert set(budget) == {r.key for r in FULL_ACRYLAMIDE_REACTIONS}
    assert budget["a_acr_cys"] > 0.0, "the measured scavenging channel is dead"
    assert all(v >= -1e-15 for v in budget.values())


def test_network_shape_is_pinned():
    described = describe_acrylamide()
    # B2.1 added four species to the sulfur block (TTCA, BND_F, PRB,
    # PROT_SS) and B2.2 added a fifth (ACID, the titratable organic-acid pool),
    # all of which the acrylamide state vector inherits unchanged. The
    # acrylamide lane composes ZERO sulfur reactions, asserted two lines below,
    # so a wider state vector changes no acrylamide prediction.
    # 56 through B2.2; B2.3 adds ONE terminal accounting pool (CBX) to the
    # SULFUR block this lane inherits, and no species any rate law can see.
    # B7 adds SEVEN more, all inherited and none of them acrylamide's: five on
    # the trunk (the furanic channel) and two on the sulfur block. 57 through
    # B6.
    assert described["n_species"] == 64
    # ... and ELEVEN more trunk REACTIONS, which this lane DOES compose,
    # because the furanic channel hangs on the trunk rather than in a lane of
    # its own. That is deliberate and it is what lets a glucose/alanine pot --
    # which has no glycine and therefore no Amadori compound -- still make HMF
    # and DMHF through the AMINE-FREE routes of Kocadagli's melt.
    assert described["trunk_reactions"] == 26
    assert described["acrylamide_reactions"] == 16
    assert described["n_reactions"] == 42
    assert described["sulfur_reactions_composed_in"] == 0
    assert described["reference_temperature_K"] == T_REF_A_K
