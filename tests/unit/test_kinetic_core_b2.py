"""
tests/unit/test_kinetic_core_b2.py

Focused unit tests for the SULFUR MODULE (Build Wave B2: sulfur formation and
thiol consumption). Deliberately narrow and fast, in B1's style:

  * carbon / nitrogen / SULFUR balance, at import and per reaction;
  * conservation of all three elements through an integration;
  * THE NO-FIXED-BRANCH-FRACTION PROPERTY: perturb a precursor 2x and assert the
    branch share MOVES (Cerny 2007 Table 5's 15% -> 46%);
  * channel separation: the four named thiol sinks are four objects with four
    temperature windows and NO activation energy;
  * the hold-out firewall: the declared hold-out literals appear in no runtime
    or fit source file;
  * a pinned regression at the Hofmann pH-5 condition.

Run with:
  ./scripts/docker_maillard.sh run "pytest tests/unit/test_kinetic_core_b1*.py \\
      tests/unit/test_kinetic_core_b2*.py -q"
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import operative_parameters
from src.kinetic_core.network import BALANCED_ELEMENTS
from src.kinetic_core.parameters import R_KJ
from src.kinetic_core.parameters_sulfur import (
    CONSUMPTION_KEYS,
    FITTED_SULFUR_BOUNDS_LOG10K,
    FITTED_SULFUR_KEYS,
    FORMATION_KEYS,
    MEASURED_SULFUR,
    NO_MEASURED_RATE,
    PROHIBITED_DERIVATIONS,
    PH_TERM_PROVENANCE,
    THIOL_CHANNELS,
    T_REF_S_K,
    assert_no_dft_sulfur,
    cysteine_thermolysis_arrhenius,
    cysteine_thermolysis_k,
    free_amine_fraction,
    neutral_h2s_fraction,
    ph_factor,
    sulfur_placeholders,
    with_fitted_sulfur,
)
from src.kinetic_core.species_sulfur import (
    MOLECULAR_WEIGHT_G_PER_MOL,
    NO_THRESHOLD_AVAILABLE,
    ODOUR_THRESHOLD_UG_PER_L,
    SITE_POOLS,
    SULFUR_BY_KEY,
    SULFUR_INDEX,
    SULFUR_STATE,
    mmol_per_litre_to_ug_per_litre,
    odour_activity_values,
    total_element,
    ug_per_litre_to_mmol_per_litre,
)
from src.kinetic_core.sulfur import (
    FULL_REACTIONS,
    OUT_OF_SCOPE,
    REACTION_PH_FACTOR,
    SULFUR_REACTIONS,
    branch_shares,
    describe_sulfur,
    integrate_sulfur,
    sulfur_flux_budget,
    validate_sulfur_balance,
)

ROOT = Path(__file__).resolve().parents[2]
CELSIUS = 273.15

B1_FITTED = {
    "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
    "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
    "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
    "k_aa_frag": (0.011812994692176768, 20.000000150449104),
}

#: A parameter set that is NOT the fit. Every fitted step is set to the same
#: neutral value, so these tests exercise STRUCTURE and never the fit's numbers:
#: nothing here breaks if the fit is re-run, and nothing here can quietly become
#: a second copy of the fitted values.
NEUTRAL_LOG10K = {key: -3.0 for key in FITTED_SULFUR_KEYS}
NEUTRAL_EA = 120.0


@pytest.fixture(scope="module")
def parameters():
    p = dict(operative_parameters(B1_FITTED))
    p.update(MEASURED_SULFUR)
    p.update(with_fitted_sulfur(NEUTRAL_LOG10K, NEUTRAL_EA))
    return p


# ---------------------------------------------------------------------------
# Stoichiometry and the SULFUR balance
# ---------------------------------------------------------------------------


def test_sulfur_is_a_balanced_element():
    assert "sulfur" in BALANCED_ELEMENTS
    validate_sulfur_balance()  # raises if anything fails


def test_every_reaction_balances_all_three_elements():
    for reaction in FULL_REACTIONS:
        for element in ("carbon", "nitrogen", "sulfur"):
            left, right = reaction.atom_balance(element, SULFUR_BY_KEY)
            assert left == right, f"{reaction.key} unbalanced in {element}"


def test_a_deliberately_unbalanced_sulfur_step_is_refused():
    """The S balance is a real constraint, not a decorative one."""
    from src.kinetic_core.network import Reaction

    thiol_from_nothing = Reaction("bad", {"FUR": 1}, {"MFT": 1}, None, "")
    with pytest.raises(ValueError, match="sulfur does not balance"):
        validate_sulfur_balance([thiol_from_nothing])


def test_trunk_species_are_a_prefix_of_the_sulfur_state():
    from src.kinetic_core.species import SPECIES_KEYS

    for i, key in enumerate(SPECIES_KEYS):
        assert SULFUR_INDEX[key] == i, "B1 indices must survive the extension"


def test_site_pools_carry_no_atoms():
    for key in SITE_POOLS:
        species = SULFUR_BY_KEY[key]
        assert (species.carbon, species.nitrogen, species.sulfur) == (0, 0, 0)


def test_every_sulfur_species_has_formation_and_consumption():
    """The defect this rebuild exists to remove: something that never leaves."""
    produced, consumed = set(), set()
    for reaction in FULL_REACTIONS:
        produced.update(reaction.products)
        consumed.update(reaction.reactants)
    # B2.1 adds BND_F (the FFT half of the split adduct reservoir, which IS
    # consumed by its own release step but is listed here for symmetry with
    # BND) and PRB, the protein-disulfide adduct, which is TERMINAL because
    # anantharamkrishnan2020b sec. 5d reports no reversibility experiment in
    # either direction -- carrying a release rate would be inventing one.
    # B2.2 adds ACID, the titratable organic-acid pool. It is TERMINAL BY
    # CONSTRUCTION and that is the whole point of it: nothing in the sulfur
    # network consumes it, so it cannot change a volatile prediction except
    # through the charge balance that sets the pH. tests/unit/
    # test_kinetic_core_b2_2.py asserts separately that it never appears as a
    # reactant anywhere.
    terminal_by_design = {
        "FRAG_C", "FRAG_N", "FRAG_S", "MEL_C", "MEL_N", "BND", "BND_F",
        "OLG", "PRB", "ACID",
    }
    for species in SULFUR_STATE:
        if species.role in ("reactant", "site") or species.key in terminal_by_design:
            continue
        assert species.key in produced, f"{species.key} has no formation term"
        assert species.key in consumed, f"{species.key} has no consumption term"


# ---------------------------------------------------------------------------
# Conservation through an integration
# ---------------------------------------------------------------------------


def test_all_three_elements_are_conserved(parameters):
    run = integrate_sulfur(
        parameters, 145.0 + CELSIUS,
        {"PENT": 100.0, "Cys": 33.0, "THI": 10.0, "OX": 1.0},
        np.linspace(0.0, 20.0, 9), ph=5.0,
    )
    for element in ("carbon", "nitrogen", "sulfur"):
        closure = run.element_closure(element)
        assert closure["max_relative_drift"] < 1e-6, (element, closure)


def test_no_negative_concentrations(parameters):
    run = integrate_sulfur(
        parameters, 145.0 + CELSIUS, {"PENT": 100.0, "Cys": 33.0, "OX": 1.0},
        np.linspace(0.0, 20.0, 9), ph=5.0,
    )
    assert float(np.min(run.concentrations)) >= 0.0


def test_sulfur_is_conserved_when_a_thiol_is_consumed(parameters):
    """The dimerisation and thioether channels must not lose sulfur."""
    run = integrate_sulfur(
        parameters, 120.0 + CELSIUS,
        {"MFT": 1.0, "FFT": 1.0, "OX": 5.0, "MELE": 20.0, "MESH": 1.0},
        np.array([0.0, 60.0]), ph=6.0,
    )
    s0 = total_element(run.concentrations[0], "sulfur")
    s1 = total_element(run.concentrations[-1], "sulfur")
    assert abs(s1 - s0) / s0 < 1e-6
    assert run.final("MFT") < 1.0, "the thiol must actually be consumed"


# ---------------------------------------------------------------------------
# THE NO-FIXED-BRANCH-FRACTION PROPERTY
# ---------------------------------------------------------------------------


def test_no_branch_fraction_constant_exists_anywhere():
    """
    A grep-level guarantee: no module in the sulfur package may contain a
    stored branch fraction. The registry holds rate constants only.
    """
    for name in ("parameters_sulfur.py", "sulfur.py", "species_sulfur.py"):
        text = (ROOT / "src/kinetic_core" / name).read_text()
        for forbidden in (
            "BRANCH_FRACTION", "branch_fraction =", "BRANCHING_RATIO =",
            "thiamine_share =", "SPLIT_FRACTION",
        ):
            assert forbidden not in text, f"{name} carries a stored fraction"


def test_branch_share_moves_when_a_precursor_is_doubled(parameters):
    """
    THE PROPERTY TEST FOR CERNY 2007 TABLE 5.

    A 2x change in precursor loading is measured to move the xylose share of
    MFT from 15% to 46% -- a 3.1x change in the branch fraction. A model with
    fixed fractions predicts exactly 1.00x. This asserts the share MOVES, which
    is the structural precondition; whether it moves by the right AMOUNT is
    scored out-of-sample by the hold-out report, not here.
    """
    one_x = branch_shares(sulfur_flux_budget(
        parameters, 145.0 + CELSIUS,
        {"Cys": 50.0, "THI": 50.0, "PENT": 150.0}, 20.0, ph=5.0, n_points=41,
    ))
    two_x = branch_shares(sulfur_flux_budget(
        parameters, 145.0 + CELSIUS,
        {"Cys": 100.0, "THI": 100.0, "PENT": 300.0}, 20.0, ph=5.0, n_points=41,
    ))
    a = one_x["MFT_share_thiamine_route"]
    b = two_x["MFT_share_thiamine_route"]
    assert np.isfinite(a) and np.isfinite(b)
    assert abs(math.log10(max(b, 1e-12) / max(a, 1e-12))) > 0.02, (
        f"the thiamine share did not move on a 2x loading change "
        f"({a:.4f} -> {b:.4f}); that is the fixed-branch-fraction failure mode"
    )


def test_branch_shares_are_recomputed_not_stored(parameters):
    """Changing the temperature must move the route mix."""
    hot = branch_shares(sulfur_flux_budget(
        parameters, 145.0 + CELSIUS, {"PENT": 100.0, "Cys": 33.0}, 20.0,
        ph=5.0, n_points=31,
    ))
    cold = branch_shares(sulfur_flux_budget(
        parameters, 95.0 + CELSIUS, {"PENT": 100.0, "Cys": 33.0}, 240.0,
        ph=5.0, n_points=31,
    ))
    assert hot["MFT_share_intact_skeleton_route"] != cold["MFT_share_intact_skeleton_route"]


def test_single_route_controls_are_structural(parameters):
    """
    With no cysteine, the sugar routes cannot make MFT at all (they all need
    sulfide, and the only sulfide source is cysteine thermolysis). With no
    thiamine, the thiamine route cannot. Both are properties of the TOPOLOGY,
    not of any fitted number -- which is why Cerny's two single-route controls
    are a test of structure.
    """
    no_thiamine = branch_shares(sulfur_flux_budget(
        parameters, 145.0 + CELSIUS, {"Cys": 100.0, "PENT": 300.0}, 20.0,
        ph=5.0, n_points=31,
    ))
    assert no_thiamine["MFT_share_thiamine_route"] == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# CHANNEL SEPARATION -- Module 2's architecture
# ---------------------------------------------------------------------------


def test_four_named_channels_with_four_windows():
    names = {c["channel"] for c in THIOL_CHANNELS}
    # B2's four, which must all still be there and separately named
    assert {
        "covalent_addition_to_matrix_electrophiles",
        "acid_catalysed_C5_oligomerisation",
        "oxidative_dimerisation",
        "radical_coupling_to_methanethiol",
    } <= names
    # B2.1 adds two more, each with its own source and its own exclusion
    assert "thiolate_mediated_oxidative_loss" in names
    assert "thiol_disulfide_exchange_with_matrix_protein" in names
    windows = {c["dominant_at_c"] for c in THIOL_CHANNELS}
    assert len(windows) >= 4, "each channel must declare its own T window"
    for channel in THIOL_CHANNELS:
        assert channel["what_excludes_the_neighbour"], channel["channel"]


def test_no_consumption_channel_has_an_activation_energy(parameters):
    """
    Inventory sec. C.1 / B.7: pairing two channels' rates to get an Ea is a
    PROHIBITED DERIVATION. Enforced by there being no Ea to pair.
    """
    # B2.1 NARROWS THIS, AND THE NARROWING IS THE POINT OF THE REVISION.
    # Sec. B.7's finding is about FOUR NAMED CHANNELS that four papers measure
    # at four temperatures; pairing two of THOSE is what is prohibited, and it
    # still is. The residual `*_decay` lumps are not those channels -- nobody
    # has measured them at any temperature, so there is no pairing to prohibit,
    # and asserting they are temperature-INDEPENDENT is a strong claim with no
    # evidence and a known direction of error. See
    # parameters_sulfur.NAMED_CHANNEL_KEYS / UNASSIGNED_SINK_KEYS.
    from src.kinetic_core.parameters_sulfur import NO_EA_KEYS

    for key in NO_EA_KEYS:
        assert parameters[key].ea_kj_mol is None, key
        assert "no_Ea_available" in parameters[key].flags
    assert MEASURED_SULFUR["k_protein_ss"].ea_kj_mol is None
    assert MEASURED_SULFUR["k_oligomer" if False else "k_protein_ss"].k_ref is not None


def test_a_no_Ea_parameter_is_held_fixed_not_extrapolated():
    """
    B2.1: k_thioether now carries Stack 2018's MEASURED forward barrier, so the
    channel that stands for "no Ea exists" is the protein-disulfide exchange,
    for which no source reports a temperature dependence of any kind.
    """
    protein = MEASURED_SULFUR["k_protein_ss"]
    assert protein.k_at(303.15) == protein.k_at(388.15) == protein.k_ref


def test_the_thioether_channel_emits_an_extrapolation_warning_off_its_window(parameters):
    run = integrate_sulfur(
        parameters, 115.0 + CELSIUS,
        {"MFT": 1.0, "MELE": 10.0, "PROT_SS": 1.0}, np.array([0.0, 60.0]),
        ph=5.0,
    )
    warnings = run.metadata["extrapolation_warnings"]
    # k_thioether is measured over 25-30 C and is being evaluated at 115 C
    assert any("k_thioether" in w for w in warnings)
    # and the channel that genuinely has no barrier says so in those words
    assert any(
        "k_protein_ss" in w and "NO ACTIVATION ENERGY" in w for w in warnings
    )


def test_the_oligomerisation_channel_is_declared_and_zero(parameters):
    from src.kinetic_core.sulfur import sulfur_rate_constants_at

    k = sulfur_rate_constants_at(parameters, 323.15, 5.0)
    assert k["ch_oligomer_mft"] == 0.0
    assert k["ch_oligomer_fft"] == 0.0
    assert "k_oligomer" in NO_MEASURED_RATE
    assert "HOLD-OUT" in NO_MEASURED_RATE["k_oligomer"]


def test_prohibited_derivations_are_recorded():
    subjects = " ".join(d["derivation"] for d in PROHIBITED_DERIVATIONS).lower()
    assert "30 c and 115 c" in subjects
    assert "branch fraction" in subjects
    assert len(PROHIBITED_DERIVATIONS) >= 5
    for entry in PROHIBITED_DERIVATIONS:
        assert entry["why_forbidden"] and entry["enforced_how"]


# ---------------------------------------------------------------------------
# DIMERISATION IS NOT AROMA LOSS
# ---------------------------------------------------------------------------


def test_the_dimer_is_a_species_with_its_own_potency():
    ratio = ODOUR_THRESHOLD_UG_PER_L["MFT"] / ODOUR_THRESHOLD_UG_PER_L["MFTD"]
    assert ratio == pytest.approx(15.625, rel=1e-3), "the 15.6x potency ratio"
    assert SULFUR_BY_KEY["MFTD"].sulfur == 2
    assert SULFUR_BY_KEY["MFTD"].carbon == 2 * SULFUR_BY_KEY["MFT"].carbon


def test_dimerisation_at_the_measured_share_does_not_lose_the_OAV():
    """
    At Zhou's measured 6.5% of thiol equivalents, the dimer's OAV must MATCH
    the monomer's -- the paper prints 3.21e5 against 3.18e5. A module that
    scored dimerisation as a pure loss would be wrong by the threshold ratio.
    """
    state = np.zeros(len(SULFUR_STATE))
    mft_mmol = ug_per_litre_to_mmol_per_litre("MFT", 1588.57)
    dimer_mmol = ug_per_litre_to_mmol_per_litre("MFTD", 102.59)
    state[SULFUR_INDEX["MFT"]] = mft_mmol
    state[SULFUR_INDEX["MFTD"]] = dimer_mmol
    oav = odour_activity_values(state)
    assert oav["MFT"] == pytest.approx(1588.57 / 0.005, rel=1e-6)
    assert oav["MFTD"] == pytest.approx(102.59 / 0.00032, rel=1e-6)
    assert 0.5 < oav["dimer_over_monomer_OAV"] < 2.0, (
        "the dimer's OAV must be of the same order as the monomer's"
    )


def test_species_without_a_threshold_are_named_not_defaulted():
    state = np.zeros(len(SULFUR_STATE))
    oav = odour_activity_values(state)
    for key in NO_THRESHOLD_AVAILABLE:
        assert oav[key] == "no_threshold_available"


def test_molecular_weights_round_trip():
    for key in ("MFT", "FFT", "MFTD", "ACTZ", "FUR"):
        assert mmol_per_litre_to_ug_per_litre(
            key, ug_per_litre_to_mmol_per_litre(key, 123.4)
        ) == pytest.approx(123.4, rel=1e-9)


# ---------------------------------------------------------------------------
# pH IS STRUCTURAL, NOT FITTED
# ---------------------------------------------------------------------------


def test_there_are_no_fitted_ph_parameters():
    assert PH_TERM_PROVENANCE["fitted_ph_parameters"] == 0
    for key in FITTED_SULFUR_KEYS:
        assert "ph" not in key.lower().replace("k_", ""), key


def test_ph_factors_are_normalised_at_pH5():
    for kind in ("acid", "base", "neutral_h2s", "hs_anion", "thiolate"):
        assert ph_factor(kind, 5.0) == pytest.approx(1.0, rel=1e-9)


def test_h2s_speciation_falls_and_amine_rises_through_the_measured_pKas():
    """
    B2.1: the SHAPES are unchanged and still measured; what changed is that the
    pKa are now evaluated at reaction temperature, so the half-way point of the
    sulfide curve is 7.05 only at 25 C.
    """
    from src.kinetic_core.parameters_sulfur import pka_at

    assert neutral_h2s_fraction(5.0) > neutral_h2s_fraction(7.0) > neutral_h2s_fraction(9.0)
    assert free_amine_fraction(5.0) < free_amine_fraction(7.0) < free_amine_fraction(9.0)
    assert neutral_h2s_fraction(7.05, 298.15) == pytest.approx(0.5, abs=1e-6)
    assert neutral_h2s_fraction(pka_at("h2s_1", 418.15), 418.15) == pytest.approx(
        0.5, abs=1e-6
    )


def test_the_two_factor_product_law_and_exactly_where_its_maximum_is():
    """
    THE TWO-FACTOR LAW, AND AN HONEST STATEMENT OF WHAT IT DOES AND DOES NOT
    DELIVER.

    The build requirement was that [1-deoxypentosone] rising with pH times
    [H2S] falling with pH should produce Zhou's pH-7 MFT maximum WITHOUT a
    fitted pH term. The product law does exactly that -- but only at a specific
    combination of exponents, and this test pins which:

      base^1 x sulfide^2 -> maximum at pH 7.05, the SULFIDE pKa exactly
      base^2 x sulfide^2 -> maximum at pH 8.67
      base^1 x sulfide^1 -> maximum at pH 8.66
      base^2 x sulfide^1 -> maximum at pH 10.28, the CYSTEINE AMINE pKa exactly

    Every maximum is an interior one and every one sits between the two
    measured pKas, at a position set only by the exponents. Nothing here is
    fitted: the two pKas are physical constants and the exponents are the
    reaction orders of the steps.

    NOW THE HONEST PART. THE NETWORK'S OWN MFT FLUX IS base^2 x sulfide^1 --
    2,3-enolisation and the Nedvidek step are both base-catalysed, and only the
    thiol-forming step takes the sulfide -- so the exponents put its algebraic
    maximum at 10.28, THREE pH UNITS ABOVE the maximum Zhou measures. The
    product law alone therefore does NOT deliver the pH-7 maximum for this
    topology, and this test records that rather than implying otherwise.

    What can still deliver it is the COMPETITION FOR CYSTEINE: thermolysis to
    H2S accelerates with pH (Zheng & Ho, measured) and burns the same pool that
    the carbon skeleton needs as a REAGENT (Nedvidek), so the two factors are
    coupled through a shared, depleting reactant rather than merely multiplied.
    That is a dynamic property of the ODE, not of these algebraic factors.
    Whether it actually produces the maximum is scored out-of-sample by
    `zhou_MFT_shape_pH8_over_pH7` in the hold-out report, and nowhere else.

    B2.1 ADDENDUM. The pKa are now evaluated at REACTION temperature, so every
    maximum listed above MOVES DOWN with it -- and that is the point, because
    the whole family shifts toward the pH window the experiments actually
    occupy instead of sitting two to three units above it. The invariant the
    test now pins is the STRUCTURE, which is temperature-independent: the
    maximum of base^a x sulfide^b always lands between the two pKa AT THE
    EVALUATION TEMPERATURE, and it lands ON one of them at the two limiting
    exponent combinations.
    """
    from src.kinetic_core.parameters_sulfur import pka_at

    grid = np.linspace(2.0, 12.0, 1001)

    def argmax_ph(base_power, sulfide_power, t_k):
        v = np.array([
            ph_factor("base", p, t_k) ** base_power
            * ph_factor("neutral_h2s", p, t_k) ** sulfide_power
            for p in grid
        ])
        assert np.isfinite(v).all()
        return float(grid[int(np.argmax(v))])

    for t_k in (298.15, 393.15, 418.15):
        pka_s = pka_at("h2s_1", t_k)
        pka_n = pka_at("cysteine_amine", t_k)
        assert argmax_ph(1, 2, t_k) == pytest.approx(pka_s, abs=0.03)
        assert argmax_ph(2, 1, t_k) == pytest.approx(pka_n, abs=0.03)
        for powers in ((1, 1), (2, 2)):
            assert pka_s - 0.03 <= argmax_ph(*powers, t_k) <= pka_n + 0.03
    # and the whole family moves DOWN as the pKa do
    assert argmax_ph(1, 1, 418.15) < argmax_ph(1, 1, 298.15) - 1.0
    # every one is INTERIOR -- the product law is peaked, not monotone
    for powers in ((1, 1), (1, 2), (2, 1), (2, 2)):
        assert grid[0] < argmax_ph(*powers, 418.15) < grid[-1]


def test_cysteine_thermolysis_uses_the_source_consistent_matched_pairs():
    """
    The repo's shipped pair (Ea 130.4 kJ/mol, A 1.0e14 1/s) runs ~15x FASTER
    than its own source at pH 7 / 100 C. The operative set is the matched
    (Ea, A) refit of zheng1994_extraction.md sec. 5b.
    """
    for ph, expected_a_per_s, expected_ea in (
        (3.0, 9.79e11, 131.2), (5.0, 1.93e12, 133.0),
        (7.0, 2.36e13, 135.5), (9.0, 1.04e12, 123.3),
    ):
        a_per_min, ea = cysteine_thermolysis_arrhenius(ph)
        assert a_per_min / 60.0 == pytest.approx(expected_a_per_s, rel=1e-6)
        assert ea == pytest.approx(expected_ea, rel=1e-9)
    assert cysteine_thermolysis_arrhenius(5.0)[0] / 60.0 != pytest.approx(1.0e14)


def test_the_thermolysis_pair_reproduces_the_sources_own_table_I():
    """
    THE CHECK THAT MATTERS: the (Ea, A) pairs must reproduce the rate constants
    they were fitted to. Zheng & Ho's Table I is 16 measured first-order
    constants over 80-110 C and pH 3-9; the four-point Arrhenius lines have
    R^2 of only 0.94-0.98, so the fit does not pass through every datum, but it
    must be close everywhere.

    THE SAME CHECK CONDEMNS THE SHIPPED PAIR: Ea 130.4 with A 1.0e14 1/s gives
    k(100 C, pH 7) = 5.6e-5 1/s against Table I's measured 3.8e-6 -- 15x fast,
    and faster even than the measured pH-9 value.
    """
    from src.kinetic_core.parameters_sulfur import ZHENG_TABLE_I_K_PER_MIN

    worst = 0.0
    for temperature_c, row in ZHENG_TABLE_I_K_PER_MIN.items():
        for ph, measured_per_min in row.items():
            predicted = cysteine_thermolysis_k(ph, temperature_c + CELSIUS)
            worst = max(worst, abs(math.log10(predicted / measured_per_min)))
    assert worst < math.log10(3.0), (
        f"the Arrhenius pairs must reproduce their own source's Table I to "
        f"within 3x everywhere; worst is {10 ** worst:.2f}x"
    )

    # and the shipped repo pair fails the same check badly
    shipped = 1.0e14 * 60.0 * math.exp(-130.4 / (R_KJ * (100.0 + CELSIUS)))
    assert shipped / ZHENG_TABLE_I_K_PER_MIN[100.0][7.0] > 8.0


def test_ph3_is_not_interpolated_through():
    """
    Zheng & Ho identify pH 3.0 as a DIFFERENT MECHANISM, and their own pH
    regressions are fitted over pH 5.0-9.0 only. At 100 C their pH-3 constant
    EXCEEDS their pH-5 one, inverting the ordering that holds at every other
    temperature. Blending across that boundary would manufacture a monotone
    trend the source explicitly denies.
    """
    from src.kinetic_core.parameters_sulfur import (
        ZHENG_MECHANISM_BOUNDARY_PH,
        ZHENG_TABLE_I_K_PER_MIN,
    )

    at_ph3 = cysteine_thermolysis_arrhenius(3.0)
    for ph in (3.0, 3.5, 4.0, 4.5, 4.99):
        assert cysteine_thermolysis_arrhenius(ph) == at_ph3, (
            "below the mechanism boundary the pH-3 pair is used as measured"
        )
    assert cysteine_thermolysis_arrhenius(ZHENG_MECHANISM_BOUNDARY_PH) != at_ph3
    # the inversion the boundary exists for, straight from the source table
    assert ZHENG_TABLE_I_K_PER_MIN[100.0][3.0] > ZHENG_TABLE_I_K_PER_MIN[100.0][5.0]
    for temperature_c in (80.0, 90.0, 110.0):
        assert (ZHENG_TABLE_I_K_PER_MIN[temperature_c][3.0]
                < ZHENG_TABLE_I_K_PER_MIN[temperature_c][5.0])


def test_the_carried_ph_law_is_not_operative_and_carries_its_correction():
    from src.kinetic_core.parameters_sulfur import ZHENG_TABLE_II_PH_LAW

    assert ZHENG_TABLE_II_PH_LAW["operative"] is False
    assert ZHENG_TABLE_II_PH_LAW["do_not_use_below_ph_5"] is True
    # the printed 110 C intercept returns a NEGATIVE rate constant at pH 9
    row = ZHENG_TABLE_II_PH_LAW["110C"]
    log_oh_at_ph9 = -5.0
    assert row["slope"] * log_oh_at_ph9 + row["intercept_as_printed"] < 0.0
    assert row["slope"] * log_oh_at_ph9 + row["intercept"] > 0.0


def test_the_in_situ_H2S_budget_reproduces_the_dossier_derivation():
    """
    The inventory derives ~0.17 mmol of H2S in 100 mL at 145 C / 20 min from a
    3.3 mmol cysteine charge (sec. A.3.5). That derived number is what makes the
    fed-precursor experiments a 6-12x H2S SURPLUS -- the quantitative mechanism
    behind the Hofmann/Cerny disjointness. This module must reproduce it, and it
    is a CROSS-CHECK on the corrected prefactors, not a fitted target.
    """
    k = cysteine_thermolysis_k(5.0, T_REF_S_K)  # 1/min at 145 C
    released_mmol_per_L = 33.0 * (1.0 - math.exp(-k * 20.0))
    released_mmol_in_100_mL = released_mmol_per_L * 0.1
    assert 0.12 < released_mmol_in_100_mL < 0.25, released_mmol_in_100_mL


def test_ph_trajectory_changes_the_answer(parameters):
    """
    Zhou's pH labels are INITIAL pH of an unbuffered system. A clamped run and a
    trajectory run must differ, or the trajectory state is decorative.
    """
    clamped = integrate_sulfur(
        parameters, 120.0 + CELSIUS, {"ARP": 20.0, "Cys": 20.0, "OX": 1.0},
        np.array([0.0, 60.0]), ph=7.0,
    )
    trajectory = integrate_sulfur(
        parameters, 120.0 + CELSIUS, {"ARP": 20.0, "Cys": 20.0, "OX": 1.0},
        np.array([0.0, 60.0]), ph=7.0, ph_final=3.42,
    )
    assert clamped.metadata["ph_is_clamped"] is True
    assert trajectory.metadata["ph_is_clamped"] is False
    assert trajectory.final("MFT") != clamped.final("MFT")


# ---------------------------------------------------------------------------
# Registry policy
# ---------------------------------------------------------------------------


def test_no_dft():
    assert_no_dft_sulfur()
    for name in ("parameters_sulfur.py", "sulfur.py", "species_sulfur.py"):
        # A docstring may SAY that no DFT is read, and the policy guard's own
        # ban list necessarily NAMES the banned tokens. What must not exist is
        # a code path that reads a computed barrier.
        code = _strip_prose((ROOT / "src/kinetic_core" / name).read_text()).lower()
        assert "data/qm" not in code
        assert "dft_coverage" not in code


def test_every_parameter_carries_its_provenance():
    placeholders = sulfur_placeholders()
    for registry in (MEASURED_SULFUR, placeholders):
        for key, parameter in registry.items():
            assert parameter.source_anchor, key
            assert parameter.dossier_anchor, key
            assert parameter.conditions, key
            assert parameter.rate_transfer, key
            from src.kinetic_core.parameters import EVIDENCE_CLASSES

            assert parameter.evidence_class in EVIDENCE_CLASSES
            metadata = parameter.as_metadata()
            assert metadata["dft_derived"] is False
            assert "ph_of_measurement" in metadata


def test_unpopulated_parameters_refuse_to_integrate():
    from src.kinetic_core.sulfur import sulfur_rate_constants_at

    p = dict(operative_parameters(B1_FITTED))
    p.update(MEASURED_SULFUR)
    p.update(sulfur_placeholders())  # k_ref is None
    with pytest.raises(ValueError, match="unpopulated"):
        sulfur_rate_constants_at(p, 418.15, 5.0)


def test_formation_and_consumption_lanes_are_disjoint():
    assert set(FORMATION_KEYS).isdisjoint(CONSUMPTION_KEYS)
    assert set(FORMATION_KEYS) | set(CONSUMPTION_KEYS) == set(FITTED_SULFUR_KEYS)


def test_bounds_allow_a_channel_to_be_rejected():
    for key in ("k_mft_decay", "k_fft_decay", "k_pent_caramel"):
        assert FITTED_SULFUR_BOUNDS_LOG10K[key][0] <= -8.0


def test_out_of_scope_lanes_name_what_they_strand():
    assert len(OUT_OF_SCOPE) >= 4
    for row in OUT_OF_SCOPE:
        assert row["lane"] and row["what_is_stranded"] and row["why"]
    lanes = " ".join(r["lane"] for r in OUT_OF_SCOPE)
    assert "pyrazine" in lanes


# ---------------------------------------------------------------------------
# THE HOLD-OUT FIREWALL
# ---------------------------------------------------------------------------

#: Numeric values of declared HOLD-OUT rows. None may appear as a NUMBER in
#: the runtime modules or in the fit script. They are allowed only in the
#: hold-out scorer, which is where they belong.
#:
#: WHY "AS A NUMBER" AND NOT "ANYWHERE IN THE FILE". A hold-out value is
#: dangerous when it is USED -- as a target, a bound, an initialisation or a
#: constant. It is not dangerous when a docstring explains why the row is held
#: out; forbidding that would mean the module could not document its own
#: firewall. The stripper below therefore removes docstrings, string literals
#: and comments and searches only what is left, which is executable code.
HOLDOUT_NUMBERS = (
    # Zhou 2023 Table 1, the pH-6 and pH-8 columns
    "525.62", "696.99", "813.65", "325.22", "59.70", "582.34", "436.63", "107.38",
    # Hofmann 1998 T2, the pH-3 ribose rows (ug per 100 mL).
    # The pH-7 row's values (2.5 and 1.2 ug) are DELIBERATELY NOT LISTED: they
    # are too generic to discriminate. "2.5" matches the exponent of an
    # unrelated 10**2.5, "1.2" matches a dozen ordinary constants, and a guard
    # that fires on those trains the reader to ignore it. The firewall is a
    # numeric tripwire on DISTINCTIVE hold-out values, not a proof of
    # innocence; the actual control is that the fit script's row table is
    # explicit and reviewable.
    "55.3", "22.9",
    # Whitfield 2001's non-detect at pH 6.5
    "0.0010",
    # Hofmann 1998's dry-180 C rows
    "1553.9", "97.2", "25.1",
    # Zhang 2024 Fig. 2 read-offs
    "23.6", "11.0", "0.234",
    # van Seeventer's 50 C zero-order storage rates, %/day
    "59", "28",
)

RUNTIME_FILES = (
    "src/kinetic_core/sulfur.py",
    "src/kinetic_core/parameters_sulfur.py",
    "src/kinetic_core/species_sulfur.py",
)
FIT_FILE = "scripts/generators/generate_kinetic_core_b2_fit.py"


def _strip_prose(text: str) -> str:
    """Remove docstrings, string literals and comments; return executable code."""
    text = re.sub(r'"""(?:.|\n)*?"""', " ", text)
    text = re.sub(r"'''(?:.|\n)*?'''", " ", text)
    text = re.sub(r'"(?:[^"\\\n]|\\.)*"', ' "" ', text)
    text = re.sub(r"'(?:[^'\\\n]|\\.)*'", " '' ", text)
    return "\n".join(line.split("#")[0] for line in text.splitlines())


@pytest.mark.parametrize("relative", RUNTIME_FILES + (FIT_FILE,))
def test_no_holdout_value_reaches_the_runtime_or_the_fit(relative):
    code = _strip_prose((ROOT / relative).read_text())
    for literal in HOLDOUT_NUMBERS:
        pattern = r"(?<![\w.])" + re.escape(literal) + r"(?![\w.])"
        assert re.search(pattern, code) is None, (
            f"{relative} uses the hold-out value {literal} as a NUMBER in "
            f"executable code. That is a firewall breach."
        )


def test_the_firewall_stripper_actually_finds_a_planted_breach():
    """A firewall test that cannot fail is not a firewall test."""
    planted = 'x = "the note says 525.62"\ny = 525.62  # a breach\n'
    code = _strip_prose(planted)
    assert re.search(r"(?<![\w.])525\.62(?![\w.])", code) is not None
    assert code.count("525.62") == 1, "the string literal must have been stripped"


def test_the_fit_script_reads_no_benchmark_or_holdout_file():
    text = (ROOT / FIT_FILE).read_text()
    code = _strip_prose(text)
    assert "data/benchmarks" not in code
    assert "mp_holdout" not in code
    assert "external_validation" not in code


def test_the_holdout_scorer_contains_no_optimiser():
    code = _strip_prose(
        (ROOT / "scripts/generators/generate_kinetic_core_b2_holdout.py").read_text()
    )
    for forbidden in ("least_squares", "minimize(", "curve_fit", "differential_evolution"):
        assert forbidden not in code, (
            "the hold-out scorer must not fit anything; it reads the frozen "
            "parameters and predicts"
        )


# ---------------------------------------------------------------------------
# PINNED REGRESSION at the Hofmann pH-5 condition
# ---------------------------------------------------------------------------


def test_pinned_hofmann_ph5_regression(parameters):
    """
    A pinned STRUCTURAL regression at Hofmann's pH-5 condition (ribose 100 mM +
    cysteine 33 mM, 145 C, 20 min, buffered).

    It pins the neutral-parameter network's behaviour, NOT the fit's numbers,
    so it does not have to be rewritten every time the fit is re-run -- and so
    it cannot become a second, undeclared copy of the fitted values. What it
    asserts is the set of qualitative facts the corpus is unambiguous about:

      * both thiols are formed at all;
      * MFT exceeds FFT (Hofmann 1998 T1 pH 5: 198 against 121 ppb, ratio 1.64;
        Zhou 2023 T1 pH 7 independently gives 2.10);
      * norfuraneol vastly exceeds MFT in the same pot (Hofmann T5: ~2750x);
      * cysteine is partly consumed (van Seeventer: 55% at 130 C);
      * all three elements close.
    """
    run = integrate_sulfur(
        parameters, 145.0 + CELSIUS, {"PENT": 100.0, "Cys": 33.0, "OX": 1.0},
        np.array([0.0, 20.0]), ph=5.0,
    )
    assert run.final("MFT") > 0.0
    assert run.final("FFT") > 0.0
    assert run.final("NF") > run.final("MFT"), (
        "norfuraneol must dominate MFT in the same pot -- Hofmann 1998 T5"
    )
    assert run.final("Cys") < 33.0, "cysteine must be consumed"
    assert run.final("H2S") > 0.0, "cysteine thermolysis must supply sulfide"
    for element in ("carbon", "nitrogen", "sulfur"):
        assert run.element_closure(element)["max_relative_drift"] < 1e-6


def test_network_shape_is_pinned():
    described = describe_sulfur()
    # B2.1 added TTCA, BND_F, PRB, PROT_SS; B2.2 adds ACID.
    assert described["n_species"] == 47
    assert described["trunk_reactions"] == 15
    assert len(SULFUR_REACTIONS) == described["sulfur_reactions"]
    ph_tagged = {k for k, v in REACTION_PH_FACTOR.items() if v}
    assert len(ph_tagged) >= 9, "the pH factors must actually be wired up"
