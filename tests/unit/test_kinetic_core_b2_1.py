"""FROZEN-WAVE REGRESSION RECORD (labelled 2026-09-03, test audit). The wave generator this file
tests is frozen (scripts/generators/WAVES.md); these tests fail only if the frozen report, the
network or the parameter tables change. They are the contract of a finished wave, not live checks
of new behaviour.


tests/unit/test_kinetic_core_b2_1.py

BUILD WAVE B2.1 -- the properties the sulfur revision must not lose.

B2's own hold-out report named one dominant defect (a structural pH slope of
~1 decade per pH unit) and several smaller ones. These tests pin the fixes as
PROPERTIES rather than as numbers, so that a later wave cannot quietly undo
them, and they pin the two firewalls this wave depends on.

They are unit tests: nothing here integrates a fitted network or reads a
results/ artefact, so the whole file runs in seconds.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import sulfur as S
from src.kinetic_core.species_sulfur import (
    SITE_POOLS,
    SULFUR_INDEX,
    TERMINAL_POOLS,
)

ROOT = Path(__file__).resolve().parents[2]
T145 = ps.T_REF_S_K
T121 = 394.15
T120 = 393.15


# ===========================================================================
# 1. THE pH TERM -- the defect B2 named, and the two structural fixes
# ===========================================================================


def test_pka_fall_with_temperature_and_by_the_right_amount():
    """
    Ionisation of a neutral acid is endothermic, so every pKa FALLS as
    temperature rises. B2 used 25 C values at 115-145 C and called it a
    declared approximation; it moves the cysteine amine by more than two units.
    """
    for name in ps.PKA_25C:
        cold = ps.pka_at(name, 298.15)
        hot = ps.pka_at(name, T145)
        assert cold == pytest.approx(ps.PKA_25C[name])
        assert hot < cold, f"{name}: pKa must fall with temperature"
    # The two that matter move INTO the operating window rather than sitting
    # three units outside it, which is what flattens the slope.
    assert 5.5 < ps.pka_at("h2s_1", T145) < 6.5
    assert 7.5 < ps.pka_at("cysteine_amine", T145) < 8.7


def test_base_lane_slope_is_shallower_than_b2s():
    """
    THE DEFECT, MEASURED. B2's base factor used the 25 C cysteine pKa of 10.28,
    which is three units above the operating window, so the free-amine fraction
    was in its asymptotic region and moved a full decade per pH unit. At
    reaction temperature the midpoint is inside the window and the same span
    costs measurably less.
    """
    def b2_slope(ph_lo, ph_hi):
        f = lambda p: 1.0 / (1.0 + 10.0 ** (10.28 - p))
        return math.log10(f(ph_hi) / f(ph_lo)) / (ph_hi - ph_lo)

    def b2_1_slope(ph_lo, ph_hi):
        return math.log10(
            ps.ph_factor("base", ph_hi, T120) / ps.ph_factor("base", ph_lo, T120)
        ) / (ph_hi - ph_lo)

    assert b2_slope(5.0, 8.0) == pytest.approx(1.0, abs=0.01)
    assert b2_1_slope(5.0, 8.0) < 0.98
    assert b2_1_slope(5.0, 8.0) > 0.0  # still monotone, still base-catalysed


def test_no_ph_sensitive_product_has_only_one_route():
    """
    THE SECOND HALF OF THE FIX, AND THE STRUCTURAL ONE. If a product has one
    route and that route is fully rate-limited by a catalyst that moves a
    decade per pH unit, the product moves a decade per pH unit BY CONSTRUCTION.
    Every pH-sensitive product must therefore be reachable by at least two
    steps carrying DIFFERENT pH factors.
    """
    producers: dict = {}
    for reaction in S.SULFUR_REACTIONS:
        kind = S.REACTION_PH_FACTOR.get(reaction.key, "")
        for product in reaction.products:
            producers.setdefault(product, set()).add(kind)

    for product in ("TDP", "DPO", "MFT", "FFT"):
        kinds = producers[product]
        assert len(kinds) > 1, (
            f"{product} is produced only through pH factor(s) {kinds}. A single "
            f"fully-catalysed route makes its pH response a decade per unit by "
            f"construction -- the defect B2's hold-out report named."
        )


def test_the_sulfide_nucleophile_has_both_protonation_states():
    """A nucleophile with two protonation states adds through both."""
    kinds = {
        S.REACTION_PH_FACTOR.get(r.key, "")
        for r in S.SULFUR_REACTIONS
        if "H2S" in r.reactants and r.key.startswith("r_")
    }
    assert "neutral_h2s" in kinds
    assert "hs_anion" in kinds


def test_the_two_sulfide_factors_move_in_opposite_directions():
    lo = ps.ph_factor("neutral_h2s", 8.0, T120) / ps.ph_factor("neutral_h2s", 5.0, T120)
    hi = ps.ph_factor("hs_anion", 8.0, T120) / ps.ph_factor("hs_anion", 5.0, T120)
    assert lo < 1.0 < hi


def test_every_ph_factor_is_one_at_the_normalisation_point_at_every_temperature():
    """
    k_ref always means "k at pH 5". If the normalisation itself carried a
    temperature dependence, the fitted constants would silently absorb one.
    """
    for kind in ("", "acid", "base", "neutral_h2s", "hs_anion", "thiolate"):
        for t_k in (323.15, 373.15, 394.15, 418.15):
            assert ps.ph_factor(kind, ps.PH_NORMALISATION, t_k) == pytest.approx(1.0)


def test_there_are_still_zero_fitted_ph_parameters():
    """
    The whole point of the revision is that the slope was fixed WITHOUT buying
    a pH knob. No fitted key may be a pH exponent, and the provenance block
    must still say zero.
    """
    assert ps.PH_TERM_PROVENANCE["fitted_ph_parameters"] == 0
    banned = ("exponent", "ph_slope", "ph_shape", "ph_power", "n_ph")
    for key in ps.FITTED_SULFUR_KEYS:
        for token in banned:
            assert token not in key.lower(), f"{key} looks like a fitted pH knob"


def test_thiol_consumption_now_carries_a_ph_dependence():
    """
    B2 had NO pH dependence on any consumption channel, so it had to load the
    whole observed FFT-versus-pH response onto formation. Kumazawa measures
    the consumption half directly.
    """
    ph_sensitive_sinks = [
        r.key for r in S.SULFUR_REACTIONS
        if S.REACTION_PH_FACTOR.get(r.key, "") == "thiolate"
    ]
    assert len(ph_sensitive_sinks) >= 4
    # and it must rise with pH, which is what a thiolate-mediated sink does
    assert ps.ph_factor("thiolate", 7.0, T121) > ps.ph_factor("thiolate", 4.0, T121)


# ===========================================================================
# 2. TEMPERATURE STRUCTURE
# ===========================================================================


def test_measured_activation_energies_are_not_fittable():
    """
    Kang's 55.1 kJ/mol is a MEASUREMENT. The fit may move the rate of
    r_cys_thermal; it may not move that barrier.
    """
    assert ps.MEASURED_EA_OVERRIDES["k_cys_thermal"] == pytest.approx(55.1)
    for lumped in (20.0, 120.0, 250.0):
        built = ps.with_fitted_sulfur({"k_cys_thermal": -3.0}, lumped)
        assert built["k_cys_thermal"].ea_kj_mol == pytest.approx(55.1)


def test_no_new_global_sulfur_ea_was_added():
    """
    The brief forbids adding a fitted global sulfur Ea. B2 already had exactly
    one lumped formation Ea and B2.1 must still have exactly one.
    """
    assert isinstance(ps.LUMPED_FORMATION_EA_BOUNDS, tuple)
    assert len(ps.LUMPED_FORMATION_EA_BOUNDS) == 2
    # STILL exactly ONE lumped, FITTED sulfur Ea. That is what the brief
    # forbids adding to, and B8 did not add to it.
    #
    # UPDATED BY B8, DELIBERATELY. This assertion used to read
    # `len(ps.MEASURED_EA_OVERRIDES) == 1`. B8 (FIT_HOLDOUT_DECLARATION.md
    # Amendment 17 clause 5) adds FOUR entries on TWO measurements: Zhang 2026
    # k17's 122.2 kJ/mol for thiol -> disulfide on both dimer channels, and its
    # k16's 85.7 for Cys-Amadori -> alpha-dicarbonyl on both amine-catalysed
    # fed-Amadori enolisations. A MEASURED barrier the fit cannot move is the
    # opposite of a new fitted global Ea, so the test's purpose is served by
    # counting the FITTED ones and pinning the measured ones by name.
    assert len(ps.MEASURED_EA_OVERRIDES) == 5
    assert set(ps.MEASURED_EA_OVERRIDES) == {
        "k_cys_thermal", "k_dimer_mft", "k_dimer_fft", "k_arp_dpo", "k_arp_tdp",
    }
    # every override must name the source it came from
    assert set(ps.MEASURED_EA_OVERRIDE_ANCHORS) == set(ps.MEASURED_EA_OVERRIDES)
    for key, anchor in ps.MEASURED_EA_OVERRIDE_ANCHORS.items():
        assert anchor.strip(), key


def test_policy_2_still_holds_for_the_named_channels():
    """
    The prohibited derivation is an Arrhenius line through the FOUR NAMED
    SINKS that four papers measure at four temperatures. No code path can pair
    them, and that prohibition is untouched.

    UPDATED BY B8. The loop used to run over every NAMED_CHANNEL_KEY. Two of
    those -- `k_dimer_mft` and `k_dimer_fft` -- now carry Zhang 2026 k17's
    MEASURED 122.2 kJ/mol, which is ONE paper's ONE step at THREE temperatures
    (R^2 0.971) and therefore not the cross-paper pairing policy 2 exists to
    prevent. The test now runs over NO_EA_KEYS, which is what "carries no
    activation energy" actually means, and separately pins that the two dimer
    channels carry the MEASUREMENT rather than the lumped fit.
    """
    built = ps.with_fitted_sulfur({k: -3.0 for k in ps.FITTED_SULFUR_KEYS}, 90.0)
    for key in ps.NO_EA_KEYS:
        assert built[key].ea_kj_mol is None, f"{key} must carry no Ea"
    assert set(ps.NO_EA_KEYS) == {"k_mmft"}
    for key in ("k_dimer_mft", "k_dimer_fft"):
        assert built[key].ea_kj_mol == pytest.approx(122.2), key
        assert built[key].ea_kj_mol != pytest.approx(90.0), key
    assert ps.MEASURED_SULFUR["k_thioether"].channel.startswith("covalent_addition")
    assert ps.MEASURED_SULFUR["k_protein_ss"].ea_kj_mol is None
    assert ps.oligomerisation_rate() == 0.0


def test_the_unassigned_sinks_are_no_longer_temperature_independent():
    """
    B2 pinned every fitted decay lump to 145 C and evaluated it unchanged at
    50 C, which is not conservatism -- it is a strong claim with a known
    direction of error, and it is most of why B2 lost FFT 18x too fast in an
    80 C brew. They now share the lumped Ea.
    """
    built = ps.with_fitted_sulfur({k: -1.0 for k in ps.UNASSIGNED_SINK_KEYS}, 90.0)
    for key in ps.UNASSIGNED_SINK_KEYS:
        assert built[key].ea_kj_mol == pytest.approx(90.0), key
        # and the constant must actually fall at low temperature
        assert built[key].k_at(323.15) < built[key].k_at(T145) / 100.0


# ===========================================================================
# 3. THE STACK 2018 CORRECTION
# ===========================================================================


def test_stack_equilibrium_is_the_measured_pair_not_b2s_range():
    """
    B2 used K = 10^2.5 = 316 M^-1 on Stack's authority. Stack measures both
    constants: K = 496/88 = 5.64 M^-1 at 19.4 C. B2 was 56x too large.
    """
    k19 = ps.THIOL_QUINONE_EQUILIBRIUM_USED_L_PER_MMOL
    assert k19 == pytest.approx(5.636e-3, rel=1e-3)
    assert k19 < 0.3162 / 40.0  # emphatically not B2's number


def test_the_adduct_equilibrium_falls_with_temperature():
    """
    The measured van 't Hoff enthalpy is NEGATIVE, so heating pushes the
    equilibrium back toward free thiol. B2's K had no temperature dependence
    at all, which is why it could not express Hofmann's 80 C brew being slower
    than his own 30 C models.
    """
    assert ps.STACK_DELTA_H_ADDUCT_KJ_MOL < 0
    k30 = ps.thiol_adduct_equilibrium_l_per_mmol(303.15)
    k80 = ps.thiol_adduct_equilibrium_l_per_mmol(353.15)
    assert k80 < k30 / 3.0


def test_the_release_rate_is_derived_and_not_a_parameter():
    """k_reverse = k_forward(T) / K(T). It is in no fitted key list."""
    assert "k_thioether_release" not in ps.FITTED_SULFUR_KEYS
    release = [r for r in S.SULFUR_REACTIONS if r.key.startswith("ch_thioether_release")]
    assert len(release) == 2
    for reaction in release:
        assert reaction.parameter_key is None


def test_each_thiol_releases_the_thiol_that_was_bound():
    """
    B2 carried ONE lumped C5 adduct for both thiols and recorded the
    consequence as a limitation: a pot containing both slowly converted FFT
    into MFT through the shared reservoir. The reservoirs are now separate.
    """
    by_key = {r.key: r for r in S.SULFUR_REACTIONS}
    assert set(by_key["ch_thioether_mft"].products) == {"BND"}
    assert set(by_key["ch_thioether_fft"].products) == {"BND_F"}
    assert "MFT" in by_key["ch_thioether_release"].products
    assert "FFT" in by_key["ch_thioether_release_fft"].products


# ===========================================================================
# 4. THE PROTEIN-DISULFIDE CHANNEL
# ===========================================================================


def test_protein_channel_is_structurally_present_and_rate_zero_without_protein():
    """
    The declaration requires the channel; the corpus supplies no rate and no
    FIT row. It is therefore present with a BOUNDED constant, and it is inert
    in a protein-free pot by MASS ACTION -- not by a flag that could be flipped.
    """
    by_key = {r.key: r for r in S.SULFUR_REACTIONS}
    for key in ("ch_protein_ss_fft", "ch_protein_ss_mft"):
        assert "PROT_SS" in by_key[key].reactants
    assert "PROT_SS" in SITE_POOLS
    assert "k_protein_ss" not in ps.FITTED_SULFUR_KEYS
    assert ps.MEASURED_SULFUR["k_protein_ss"].evidence_class == (
        "bounded_from_a_timescale_bracket"
    )
    lo, hi = ps.PROTEIN_DISULFIDE_EXCHANGE_BRACKET_L_PER_MMOL_MIN
    assert lo <= ps.PROTEIN_DISULFIDE_EXCHANGE_USED_L_PER_MMOL_MIN <= hi


def test_protein_channel_carries_zero_flux_in_a_protein_free_pot():
    import numpy as np

    from src.kinetic_core.sulfur import (
        FULL_REACTION_KEYS,
        sulfur_reaction_rates,
    )
    from src.kinetic_core.species_sulfur import initial_sulfur_state

    state = initial_sulfur_state({"MFT": 1.0, "FFT": 1.0, "PENT": 10.0})
    k = {key: 1.0 for key in FULL_REACTION_KEYS}
    rates = sulfur_reaction_rates(state, k)
    index = {key: i for i, key in enumerate(FULL_REACTION_KEYS)}
    for key in ("ch_protein_ss_fft", "ch_protein_ss_mft"):
        assert rates[index[key]] == 0.0
    # and it fires when protein IS present
    state2 = initial_sulfur_state({"MFT": 1.0, "FFT": 1.0, "PROT_SS": 1.09})
    rates2 = sulfur_reaction_rates(state2, k)
    assert rates2[index["ch_protein_ss_fft"]] > 0.0


def test_protein_bound_thiol_is_terminal_and_says_why():
    assert "PRB" in TERMINAL_POOLS
    assert "PRB" in SULFUR_INDEX


# ===========================================================================
# 5. STRUCTURAL INVARIANTS INHERITED FROM B2 -- none may be lost
# ===========================================================================


def test_three_element_balance_still_holds_over_the_extended_network():
    S.validate_sulfur_balance()  # raises if any step fails


def test_still_no_branch_fraction_constant_anywhere():
    registry_text = (ROOT / "src/kinetic_core/parameters_sulfur.py").read_text()
    assert "branch_fraction" not in registry_text.replace(
        "branch_fraction_policy", ""
    ).replace("branch_fractions_present", "")
    assert ps.sulfur_registry_metadata({})["branch_fractions_present"] is False


def test_no_dft_anywhere():
    ps.assert_no_dft_sulfur()
    ps.assert_no_dft_sulfur(ps.MEASURED_SULFUR)


# ===========================================================================
# 6. THE FIREWALL -- values this wave unavoidably SAW must not have entered
# ===========================================================================
# DISCLOSURE, recorded here as well as in the report: the Kang 140 C hold-out
# values are PRINTED in kang2026_SI_extraction.md, and the build brief directed
# this wave to read that dossier. They were therefore seen before the fit ran.
# What follows is the mechanical check that they were not used.

#: Literals that belong ONLY to the hold-out scorer.
_HOLDOUT_LITERALS = (
    # Kang 2026 SI Table S4, the 140 C column
    "5.907", "11.439", "60.400", "12.757", "73.157",
    # Kang Fig. S4's 140 C rung
    "62.6", "8.195",
    # Sun 2019's pH-9 column
    "0.4 / 0.5 / 0.5",
    # everything B2 already held out
    "525.62", "325.22", "50.07", "582.34", "436.63",
    "696.99", "813.65", "59.70", "553.0", "229.0",
)

_FIT_SIDE_FILES = (
    "src/kinetic_core/sulfur.py",
    "src/kinetic_core/parameters_sulfur.py",
    "src/kinetic_core/species_sulfur.py",
    "scripts/generators/generate_kinetic_core_b2_1_fit.py",
)


def test_no_holdout_literal_appears_on_the_fit_side():
    offences = []
    for relative in _FIT_SIDE_FILES:
        text = (ROOT / relative).read_text()
        # the fit script is allowed to NAME the hold-out in prose, so strip
        # comment-only prose lines that explicitly say the value is held out
        for literal in _HOLDOUT_LITERALS:
            for line_no, line in enumerate(text.splitlines(), 1):
                if literal not in line:
                    continue
                if "HOLD-OUT" in line.upper() or "IS NOT HERE" in line.upper():
                    continue
                offences.append(f"{relative}:{line_no}: {literal!r} in {line.strip()[:90]}")
    assert not offences, "hold-out literals leaked onto the fit side:\n" + "\n".join(offences)


def test_the_frozen_external_validation_bundles_were_never_opened():
    """
    The 21 frozen bundles under data/benchmarks/external_validation/ are
    untouchable. Naming them in a disclosure is fine and is in fact required;
    READING one is not. This checks that no occurrence of the path sits on a
    line that also does file I/O.
    """
    io_tokens = ("open(", "read_text", "read_bytes", "json.load", "glob(", "iterdir")
    for relative in _FIT_SIDE_FILES + (
        "scripts/generators/generate_kinetic_core_b2_1_holdout.py",
    ):
        for line_no, line in enumerate(
            (ROOT / relative).read_text().splitlines(), 1
        ):
            if "external_validation" not in line:
                continue
            assert not any(token in line for token in io_tokens), (
                f"{relative}:{line_no} reads a frozen bundle: {line.strip()[:100]}"
            )


def test_the_fit_never_integrates_a_holdout_condition():
    """
    This catches a hold-out condition being smuggled in as a 'system' that then
    quietly constrains a parameter. The conditions are checked PER SYSTEM, not
    by grepping the file, because 140 C and pH 4.5 are perfectly legitimate FIT
    conditions elsewhere in the panel (Whitfield 1999's fed-norfuraneol rows).
    """
    import scripts.generators.generate_kinetic_core_b2_1_fit as fit

    for name, spec in fit.SYSTEMS.items():
        t_c = float(spec["t_c"])
        ph = float(spec["ph"])
        # Kang's ladder: 100 and 120 C are FIT, 140 C is the gating hold-out
        if name.startswith("kang"):
            assert t_c in (100.0, 120.0), f"{name} at {t_c} C is the hold-out rung"
        # Zhou's columns: pH 7 is FIT, pH 6 and pH 8 are held out
        if name.startswith("zhou"):
            assert ph == 7.0, f"{name} at pH {ph} is a held-out column"
        # Sun 2019's pH-9 column, Whitfield 2001's pH 6.5 collapse
        assert ph != 9.0, f"{name} sits on Sun 2019's held-out pH-9 column"
        assert ph != 6.5, f"{name} sits on Whitfield 2001's held-out pH"
        # Hofmann's own held-out pH rows are pH 3 and pH 7 in the AQUEOUS
        # free-sugar pots; the fed-precursor T10 pH ladder at 3/5/7 IS declared
        # FIT, so the exclusion is by system and not by pH alone.
        if name.startswith("hofmann_") and "pH5" not in name:
            raise AssertionError(f"unexpected Hofmann fit system {name}")


def test_branch_shares_count_both_halves_of_every_split_route():
    """
    A ROUTE-SHARE BUG THIS WAVE FOUND IN ITSELF, pinned so it cannot return.
    Splitting a step into two parallel elementary steps (neutral H2S and HS-)
    silently breaks any share that names one of them: the intact-skeleton share
    would be under-counted, and worse, it would then MOVE WITH pH for a reason
    that has nothing to do with the route mix -- exactly the artefact the
    no-fixed-fraction architecture exists to avoid.
    """
    flux = {key: 0.0 for key in S.FULL_REACTION_KEYS}
    flux["r_ddp_mft"] = 1.0
    flux["r_ddp_mft_hs"] = 3.0
    flux["r_hmp_mft"] = 4.0
    shares_a = S.branch_shares(flux)
    assert shares_a["MFT_share_intact_skeleton_route"] == pytest.approx(0.5)
    assert shares_a["MFT_share_thiamine_route"] == pytest.approx(0.5)

    flux2 = dict(flux)
    flux2["r_fur_fft"] = 1.0
    flux2["r_fur_fft_hs"] = 1.0
    flux2["r_fur_decay"] = 2.0
    assert S.branch_shares(flux2)["FFT_share_of_furfural_flux"] == pytest.approx(0.5)
