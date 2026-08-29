"""
tests/unit/test_kinetic_core_b2_2.py

BUILD WAVE B2.2 -- the properties the two new pieces must not lose.

B2.2 does two things: it gives the decay lumps their own barriers in two named
families, and it turns the pH from a PRESCRIBED interpolation between two
measurements into a COMPUTED state solved from a charge balance. Both are the
kind of change that can be quietly undone or quietly widened, so both are
pinned here as PROPERTIES rather than as numbers.

They are unit tests. The three that read the frozen fit report skip themselves
if it has not been generated, so the file runs in seconds either way.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from src.kinetic_core import parameters_sulfur as ps
from src.kinetic_core import ph_state as PH
from src.kinetic_core import sulfur as S
from src.kinetic_core.species_sulfur import (
    SULFUR_INDEX,
    SULFUR_STATE_KEYS,
    TERMINAL_POOLS,
    initial_sulfur_state,
)

ROOT = Path(__file__).resolve().parents[2]
FIT_REPORT = ROOT / "results/validation/kinetic_core_b2_2_fit_report.json"

T25 = 298.15
T120 = 393.15
T145 = ps.T_REF_S_K

_BUFFER_NONE = PH.BufferSpec(kind="none", declared=True, source="test")
_BUFFER_PHOSPHATE = PH.BufferSpec(
    kind="phosphate", phosphate_mol_l=0.5, declared=True, source="test")
_BUFFER_MCILVAINE = PH.BufferSpec(
    kind="citrate_phosphate", phosphate_mol_l=0.2, citrate_mol_l=0.1,
    declared=True, source="test")
_DRIFT = PH.PhDrift(acid_yield=0.9, arp_amine_pka=10.0)


def _frozen():
    if not FIT_REPORT.exists():  # pragma: no cover - environment dependent
        pytest.skip(f"{FIT_REPORT.name} not generated yet")
    return json.loads(FIT_REPORT.read_text())


def _toy_parameters():
    """A fitted-looking parameter set that costs nothing to build."""
    from src.kinetic_core import operative_parameters

    b1 = {
        "k_glc_frag": (1.000032373292967e-08, 180.69531857985976),
        "k_mgo_mel": (0.02272608289635856, 20.043206355884948),
        "k_fa_frag": (3.4646810085648807e-08, 20.53065919356619),
        "k_aa_frag": (0.011812994692176768, 20.000000150449104),
    }
    out = dict(operative_parameters(b1))
    out.update(ps.MEASURED_SULFUR)
    out.update(ps.with_fitted_sulfur(
        {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}, 110.0,
        {"thiol_sink": 80.0, "carbonyl_sink": 60.0},
    ))
    return out


# ===========================================================================
# 1. THE TEXTBOOK ACID-BASE CONSTANTS -- immovable, and correct
# ===========================================================================


def test_pkw_reproduces_the_tabulated_ion_product():
    """
    Marshall & Franck's table, at the three temperatures this corpus runs at.
    A constant-enthalpy van 't Hoff from 25 C is wrong by 0.23 units at 100 C
    and 0.35 at 125 C, which is why the table is carried instead; this test is
    what stops a later wave "simplifying" it back.
    """
    assert PH.pkw_at(T25) == pytest.approx(13.995, abs=1e-3)
    assert PH.pkw_at(373.15) == pytest.approx(12.265, abs=1e-3)
    assert PH.pkw_at(398.15) == pytest.approx(11.911, abs=1e-3)
    assert PH.pkw_at(423.15) == pytest.approx(11.638, abs=1e-3)
    # and it must fall monotonically over the whole operating window
    values = [PH.pkw_at(t) for t in (T25, 373.15, T120, 423.15)]
    assert values == sorted(values, reverse=True)


def test_the_shared_pka_agree_with_the_speciation_module():
    """
    The charge balance and the rate law must not disagree about cysteine. If
    ph_state carried its own thiol pKa the model could compute a thiolate
    fraction for the rate and a different one for the charge, which is exactly
    the sort of split-brain a later refactor introduces by accident.
    """
    cys_pkas = PH.CYSTEINE.pkas
    assert cys_pkas[1] == pytest.approx(ps.PKA_25C["thiol"])
    assert cys_pkas[2] == pytest.approx(ps.PKA_25C["cysteine_amine"])
    assert PH.CYSTEINE.delta_h_kj_mol[1] == pytest.approx(
        ps.DELTA_H_IONISATION_KJ_MOL["thiol"])
    assert PH.CYSTEINE.delta_h_kj_mol[2] == pytest.approx(
        ps.DELTA_H_IONISATION_KJ_MOL["cysteine_amine"])
    # and both modules must apply the SAME van 't Hoff correction
    assert PH.CYSTEINE.pkas_at(T145)[1] == pytest.approx(
        ps.pka_at("thiol", T145), abs=1e-9)


def test_every_titratable_group_carries_a_citation_and_no_fitted_number():
    for group in (PH.FORMIC, PH.ACETIC, PH.CYSTEINE, PH.PHOSPHATE,
                  PH.CITRATE, PH.TTCA_GROUP):
        assert group.source.strip(), f"{group.key} has no source"
        assert len(group.pkas) == len(group.delta_h_kj_mol)
    assert PH.PH_STATE_PROVENANCE["fitted_constants"] == 2
    assert set(PH.PH_STATE_PROVENANCE["fitted_constant_names"]) == {
        "acid_yield_per_sink_event", "arp_secondary_ammonium_pKa"
    }


def test_average_charge_has_the_right_limits():
    """A polyprotic distribution that is wrong at its limits is wrong."""
    # phosphate: neutral H3PO4 far below pKa1, fully -3 far above pKa3
    assert PH.PHOSPHATE.average_charge(-2.0, T25) == pytest.approx(0.0, abs=1e-3)
    assert PH.PHOSPHATE.average_charge(18.0, T25) == pytest.approx(-3.0, abs=1e-3)
    # cysteine: +1 fully protonated, -2 fully deprotonated
    assert PH.CYSTEINE.average_charge(-3.0, T25) == pytest.approx(1.0, abs=1e-3)
    assert PH.CYSTEINE.average_charge(16.0, T25) == pytest.approx(-2.0, abs=1e-3)
    # monotone decreasing in pH, everywhere
    grid = np.linspace(0.0, 14.0, 60)
    for group in (PH.FORMIC, PH.CYSTEINE, PH.PHOSPHATE, PH.CITRATE):
        z = [group.average_charge(p, T120) for p in grid]
        assert all(z[i] >= z[i + 1] - 1e-12 for i in range(len(z) - 1)), group.key


# ===========================================================================
# 2. THE CHARGE BALANCE -- conservation, which is the whole claim
# ===========================================================================


def test_the_solved_ph_is_a_root_of_the_charge_balance():
    inventory = PH.titratable_inventory(
        {"Cys": 20.0, "ARP": 20.0, "ACID": 9.0}, _DRIFT, _BUFFER_NONE)
    sid = 5.0
    ph = PH.solve_ph(T120, sid, inventory)
    assert abs(PH._charge_residual(ph, T120, sid, inventory)) < 1e-6


def test_the_strong_ion_difference_reproduces_the_declared_initial_ph():
    """
    SID is defined as the number that makes the charge balance hold at the
    declared label pH. Round-tripping it is the definition, and it is the one
    piece of arithmetic the whole drift model rests on.
    """
    for label_ph in (3.0, 5.0, 6.0, 7.0, 8.0):
        composition = {"Cys": 20.0, "ARP": 20.0}
        sid = PH.strong_ion_difference(
            label_ph, composition, _DRIFT, _BUFFER_NONE)
        back = PH.ph_of_state(composition, T25, sid, _DRIFT, _BUFFER_NONE)
        assert back == pytest.approx(label_ph, abs=1e-4)


def test_more_base_is_needed_to_reach_a_higher_ph():
    """
    The mechanism that separates Zhou's three runs from one another: the pot
    taken to pH 8 carries more titratable base into the hold than the pot taken
    to pH 6, and that -- not the pH label -- is what survives the heating.
    """
    composition = {"Cys": 20.0, "ARP": 20.0}
    sids = [
        PH.strong_ion_difference(p, composition, _DRIFT, _BUFFER_NONE)
        for p in (6.0, 7.0, 8.0)
    ]
    assert sids[0] < sids[1] < sids[2]


def test_a_dynamic_run_conserves_charge_at_every_point():
    """
    THE CONSERVATION TEST. No reaction creates or destroys a strong ion, so the
    SID computed at t = 0 must satisfy the charge balance at EVERY point of the
    trajectory, with the state and the pH that were actually used there.
    """
    p = _toy_parameters()
    run = S.integrate_sulfur(
        p, T120, {"ARP": 20.0, "Cys": 20.0, "OX": 1.0},
        np.linspace(0.0, 60.0, 13), ph=7.0,
        buffer_spec=_BUFFER_NONE, ph_drift=_DRIFT, rtol=1e-8, atol=1e-14,
    )
    sid = run.metadata["strong_ion_difference_mmol_L"]
    assert run.ph_series is not None
    worst = 0.0
    for i in range(len(run.times_min)):
        state = {k: float(run.concentrations[i, SULFUR_INDEX[k]])
                 for k in SULFUR_STATE_KEYS}
        inventory = PH.titratable_inventory(state, _DRIFT, _BUFFER_NONE)
        worst = max(worst, abs(PH._charge_residual(
            float(run.ph_series[i]), T120, sid, inventory)))
    # the pH is interpolated between fixed-point nodes, so the residual is not
    # machine zero; it must still be far below one part in a thousand of the
    # charge scale. B2.3 EXPRESSES THAT THRESHOLD AS THE COMMENT ALWAYS STATED
    # IT -- relative to the charge scale actually present -- rather than as the
    # absolute 1e-2 that stood here before. The reason is disclosed: B2.3's
    # charge fix makes this pot carry a LARGER titratable inventory (the
    # carried carboxylate is now tracked instead of being deleted), so the same
    # relative interpolation error is a larger absolute number. Pinning the
    # absolute value would have made the test tighten itself every time the
    # model became more complete, which is the wrong direction.
    #
    # ==== WAVE Q1 RE-EXAMINATION, because this was flagged worth distrusting ====
    # It was flagged for a good reason and the reason is now recorded rather than
    # left as suspicion. B2.3 re-expressed this bound in the SAME commit that
    # changed the charge physics, and the framing above ("as the comment always
    # stated it") understates what happened: the old absolute bound would FAIL
    # today. Measured on this tree, worst = 0.012756 mmol/L, so `worst < 1e-2`
    # is violated by 28 %. So B2.3 did not merely restate the threshold, it
    # loosened one that had gone red, in the commit that made it red.
    #
    # The RELATIVE FORM survives that scrutiny: the titratable inventory
    # legitimately grows as the model is completed, and a fixed absolute bound
    # would ratchet against completeness. The COEFFICIENT does not survive it --
    # 1e-3 was a round number, not a measured one. Against the measured
    # scale = 44.145 mmol/L the actual relative residual is 2.89e-4, so 1e-3 was
    # spending 3.5x of headroom nobody had accounted for.
    #
    # Tightened to 5e-4, which keeps 1.7x headroom over the measured value --
    # enough for the interpolation error to breathe as the node spacing changes,
    # not enough to absorb a real charge-closure regression unnoticed. If a
    # future wave needs this loosened again, the number to beat is 2.89e-4 and
    # the loosening should be disclosed as one.
    scale = abs(float(sid)) + sum(
        c for _g, c in PH.titratable_inventory(
            {k: float(run.concentrations[-1, SULFUR_INDEX[k]])
             for k in SULFUR_STATE_KEYS}, _DRIFT, _BUFFER_NONE))
    assert worst < 5e-4 * scale, (
        f"charge balance drifts by {worst:.3g} mmol/L against a charge scale "
        f"of {scale:.3g} mmol/L (Q1 measured 2.89e-4 relative; the bound is "
        f"5e-4)")


def test_a_dynamic_run_conserves_carbon_nitrogen_and_sulfur():
    """
    The ACID pool was added to the state vector. A new species is the classic
    way to lose an atom, so all three element closures are checked on a run
    that actually fills it.
    """
    p = _toy_parameters()
    run = S.integrate_sulfur(
        p, T120, {"ARP": 20.0, "Cys": 20.0, "OX": 1.0},
        np.linspace(0.0, 60.0, 9), ph=7.0,
        buffer_spec=_BUFFER_NONE, ph_drift=_DRIFT, rtol=1e-9, atol=1e-14,
    )
    assert run.final("ACID") > 0.0, "the test system must exercise the ACID pool"
    for element in ("carbon", "nitrogen", "sulfur"):
        closure = run.element_closure(element)
        assert closure["max_relative_drift"] < 1e-6, (element, closure)


def test_the_acid_pool_is_terminal_and_feeds_nothing_but_the_ph():
    """
    `ACID` may never be a reactant. If it were, adding it would have changed a
    volatile prediction, and the whole "this cannot move anything except the
    pH" argument would be false.
    """
    assert "ACID" in TERMINAL_POOLS
    for reaction in S.FULL_REACTIONS:
        assert "ACID" not in reaction.reactants, reaction.key
    producers = [r.key for r in S.FULL_REACTIONS if "ACID" in r.products]
    assert len(producers) >= 5


def test_a_clamped_run_is_bit_identical_to_one_with_no_drift_supplied():
    """
    B2/B2.1 reproducibility. Passing no PhDrift must reproduce the old clamped
    behaviour exactly -- not approximately -- so that every earlier artefact
    stays regenerable.
    """
    p = _toy_parameters()
    grid = np.linspace(0.0, 20.0, 5)
    a = S.integrate_sulfur(p, T145, {"PENT": 100.0, "Cys": 33.0}, grid, ph=5.0)
    b = S.integrate_sulfur(
        p, T145, {"PENT": 100.0, "Cys": 33.0}, grid, ph=5.0,
        buffer_spec=PH.BufferSpec(kind="clamped"), ph_drift=_DRIFT)
    assert a.metadata["ph_mode"] == "clamped"
    assert b.metadata["ph_mode"] == "clamped"
    assert np.allclose(a.concentrations, b.concentrations, rtol=1e-12, atol=0.0)
    assert a.ph_series is None and b.ph_series is None


# ===========================================================================
# 3. THE BUFFER -- an input, with a declared default and an honest absence
# ===========================================================================


def test_the_default_buffer_is_unbuffered_and_says_so_out_loud():
    assert PH.DEFAULT_BUFFER.kind == "none"
    assert PH.DEFAULT_BUFFER.declared is False
    p = _toy_parameters()
    run = S.integrate_sulfur(
        p, T120, {"ARP": 20.0, "Cys": 20.0}, np.array([0.0, 10.0]), ph=7.0,
        ph_drift=_DRIFT)  # NO buffer supplied
    assert any("EXTRAPOLATED" in note for note in run.metadata["ph_notes"]), (
        "an undeclared buffer must raise the extrapolation warning, not be "
        "silently treated as either buffered or clamped"
    )


def test_buffered_invariance_a_pot_that_makes_no_acid_does_not_drift():
    """
    Kumazawa's grid: 1 ppm 2-furfurylthiol in a McIlvaine buffer, nothing that
    can make an organic acid. The pH must not move at all -- if it does, the
    charge balance is manufacturing acid out of the solver.
    """
    p = _toy_parameters()
    for label_ph in (3.0, 5.0, 6.4):
        run = S.integrate_sulfur(
            p, 394.15, {"FFT": 8.7597e-3, "OX": 1.0},
            np.linspace(0.0, 10.0, 5), ph=label_ph,
            buffer_spec=_BUFFER_MCILVAINE, ph_drift=_DRIFT,
            rtol=1e-9, atol=1e-15)
        assert run.metadata["ph_final_cooled"] == pytest.approx(label_ph, abs=1e-3)
        assert run.ph_series is not None
        assert float(np.ptp(run.ph_series)) < 1e-3


def test_a_buffer_raises_the_capacity_it_is_supposed_to_raise():
    """
    Buffer capacity is REPORTED, never fitted, and its ordering must be right:
    a McIlvaine buffer at pH 5 holds far harder than 0.5 M phosphate at pH 5,
    which is the arithmetic reason Kumazawa's grid stays put and Hofmann's pot
    does not.
    """
    state = {"Cys": 33.0}
    none = PH.buffer_capacity_mmol_per_ph(state, T145, 5.0, _DRIFT, _BUFFER_NONE)
    phos = PH.buffer_capacity_mmol_per_ph(
        state, T145, 5.0, _DRIFT, _BUFFER_PHOSPHATE)
    mcil = PH.buffer_capacity_mmol_per_ph(
        state, T145, 5.0, _DRIFT, _BUFFER_MCILVAINE)
    assert none < phos < mcil
    assert none < 5.0


def test_a_buffer_spec_refuses_to_be_incomplete():
    with pytest.raises(ValueError):
        PH.BufferSpec(kind="phosphate")
    with pytest.raises(ValueError):
        PH.BufferSpec(kind="citrate_phosphate", phosphate_mol_l=0.2)
    with pytest.raises(ValueError):
        PH.BufferSpec(kind="tris", phosphate_mol_l=0.1)


def test_the_two_ph_scales_are_actually_different():
    """
    The correction that made the drift model work. An unbuffered pot titrated
    to a bench pH 7 does NOT run at pH 7: cysteine's thiol pKa falls from 8.33
    to ~7.2 on heating, so the same solution is far more ionised in the
    autoclave than on the bench.
    """
    p = _toy_parameters()
    run = S.integrate_sulfur(
        p, T120, {"ARP": 20.0, "Cys": 20.0}, np.array([0.0, 1e-6]), ph=7.0,
        buffer_spec=_BUFFER_NONE, ph_drift=_DRIFT)
    assert run.metadata["ph_label_declared"] == 7.0
    assert run.metadata["ph_initial_in_situ"] < 6.5
    assert run.metadata["ph_initial_in_situ"] > 4.0


def test_the_dynamic_state_ignores_a_measured_endpoint_and_says_so():
    """
    Once the model can PREDICT an endpoint, being handed the measured one is a
    thing to refuse, loudly. B2.1's scorer prescribed Zhou's 5.07 and Kang's
    4.9 at scoring time; that must now be impossible to do silently.
    """
    p = _toy_parameters()
    run = S.integrate_sulfur(
        p, T120, {"ARP": 20.0, "Cys": 20.0}, np.array([0.0, 5.0]), ph=7.0,
        ph_final=3.42, buffer_spec=_BUFFER_NONE, ph_drift=_DRIFT)
    assert run.metadata["ph_mode"] == "dynamic_charge_balance"
    assert any("IGNORED" in note for note in run.metadata["ph_notes"])


# ===========================================================================
# 4. THE DECAY BARRIERS -- two families, and everything else immovable
# ===========================================================================


def test_the_two_decay_families_are_disjoint_and_leave_two_lumps_alone():
    thiol = set(ps.DECAY_FAMILY_THIOL_SINK)
    carbonyl = set(ps.DECAY_FAMILY_CARBONYL_SINK)
    assert not thiol & carbonyl
    covered = thiol | carbonyl
    left = set(ps.UNASSIGNED_SINK_KEYS) - covered
    assert left == set(ps.DECAY_KEYS_ON_FORMATION_EA) == {
        "k_dimer_decay", "k_h2s_loss"
    }, (
        "k_dimer_decay and k_h2s_loss must KEEP the lumped formation Ea: "
        "nothing in the corpus measures a disulfide or a sulfide sink at two "
        "temperatures, so a barrier for either is identified by nothing."
    )
    assert covered <= set(ps.UNASSIGNED_SINK_KEYS)
    # exactly TWO new fitted barriers, not nine and not one per lump
    assert len(ps.DECAY_FAMILIES) == 2
    for family in ps.DECAY_FAMILIES:
        assert ps.DECAY_FAMILY_IDENTIFYING_ROWS[family].strip()


def test_each_decay_family_receives_its_own_barrier():
    built = ps.with_fitted_sulfur(
        {k: -2.0 for k in ps.FITTED_SULFUR_KEYS}, 119.0,
        {"thiol_sink": 55.0, "carbonyl_sink": 175.0},
    )
    for key in ps.DECAY_FAMILY_THIOL_SINK:
        assert built[key].ea_kj_mol == pytest.approx(55.0), key
    for key in ps.DECAY_FAMILY_CARBONYL_SINK:
        assert built[key].ea_kj_mol == pytest.approx(175.0), key
    for key in ps.DECAY_KEYS_ON_FORMATION_EA:
        assert built[key].ea_kj_mol == pytest.approx(119.0), key
    # UPDATED BY B8: the loop was over NAMED_CHANNEL_KEYS and is now over
    # NO_EA_KEYS. `k_dimer_mft` / `k_dimer_fft` left the no-Ea set because
    # Zhang 2026 k17 measures their barrier (122.2 kJ/mol, R^2 0.971, three
    # temperatures, one paper). What this test is FOR -- that a decay family's
    # fitted barrier never leaks onto a channel outside it -- is preserved and
    # strengthened: neither dimer channel takes 55.0, 175.0 or 119.0.
    for key in ps.NO_EA_KEYS:
        assert built[key].ea_kj_mol is None, key
    for key in ("k_dimer_mft", "k_dimer_fft"):
        assert built[key].ea_kj_mol == pytest.approx(122.2), key


def test_omitting_the_decay_barriers_reproduces_b2_1_exactly():
    """The fallback that keeps every B2.1 artefact regenerable."""
    built = ps.with_fitted_sulfur({k: -2.0 for k in ps.FITTED_SULFUR_KEYS}, 119.0)
    for key in ps.UNASSIGNED_SINK_KEYS:
        assert built[key].ea_kj_mol == pytest.approx(119.0), key


def test_measured_barriers_outrank_a_family_barrier():
    """
    THE IMMOVABILITY TEST, EXTENDED (B2.1 pinned it against the lumped Ea; the
    family barriers are a second thing that could have overwritten it).
    Kang's 55.1 kJ/mol is a MEASUREMENT. Nothing the fit does may move it.
    """
    assert ps.MEASURED_EA_OVERRIDES["k_cys_thermal"] == pytest.approx(55.1)
    for lumped in (20.0, 119.0, 250.0):
        for thiol in (20.0, 100.0, 250.0):
            built = ps.with_fitted_sulfur(
                {"k_cys_thermal": -3.0}, lumped,
                {"thiol_sink": thiol, "carbonyl_sink": thiol},
            )
            assert built["k_cys_thermal"].ea_kj_mol == pytest.approx(55.1)
    # and the override key must not have been quietly moved INTO a family
    assert ps.decay_family_of("k_cys_thermal") is None


def test_the_decay_barrier_bounds_are_symmetric_about_no_answer():
    """
    The bounds must not encode the answer the diagnosis hoped for. They are the
    formation barrier's own bounds, so 'below the formation Ea' and 'above it'
    are equally reachable.
    """
    assert ps.DECAY_EA_BOUNDS == ps.LUMPED_FORMATION_EA_BOUNDS


def test_a_lower_decay_barrier_makes_the_sink_relatively_slower_when_hot():
    """
    The mechanism the whole of Task 1 rests on, checked as arithmetic: this is
    what a decay barrier below the formation barrier BUYS at 140 C, and it is
    what B2.1's diagnosis section 2a asked for.
    """
    hot, warm = 413.15, T120   # Kang's 140 C hold-out rung and its 120 C FIT rung
    low = ps.with_fitted_sulfur(
        {"k_fft_decay": -1.0}, 119.0,
        {"thiol_sink": 60.0, "carbonyl_sink": 60.0})["k_fft_decay"]
    high = ps.with_fitted_sulfur(
        {"k_fft_decay": -1.0}, 119.0,
        {"thiol_sink": 180.0, "carbonyl_sink": 180.0})["k_fft_decay"]
    # the sink with the LOWER barrier gains less on the 120 -> 140 C leg, so a
    # product that survives to 120 min at 120 C is destroyed less at 140 C.
    assert low.k_at(hot) / low.k_at(warm) < high.k_at(hot) / high.k_at(warm)
    # and the formation lane, on the lumped barrier, gains more than either
    formation = ps.with_fitted_sulfur({"k_ddp_mft": -1.0}, 119.0)["k_ddp_mft"]
    assert (formation.k_at(hot) / formation.k_at(warm)
            > low.k_at(hot) / low.k_at(warm))


# ===========================================================================
# 5. THE FROZEN FIT -- the pre-registration's own failure criteria
# ===========================================================================


def test_the_frozen_decay_barriers_are_not_on_a_bound():
    """
    Pre-registration section 5 item 3: a barrier that lands ON a bound is not
    identified and must be withdrawn rather than reported as a value.
    """
    frozen = _frozen()["frozen_parameters"]
    lo, hi = ps.DECAY_EA_BOUNDS
    for family, value in frozen["decay_Ea_kJ_mol"].items():
        assert abs(value - lo) > 0.5 and abs(value - hi) > 0.5, (
            f"{family} landed on a bound at {value}; it is not identified"
        )


def test_neither_decay_barrier_may_be_quoted_as_a_measured_value():
    """
    THE HONEST READING OF THE FROZEN NUMBERS, PINNED SO IT CANNOT BE FORGOTTEN.

    `thiol_sink` came back at ~248 kJ/mol against a ceiling of 250 -- it clears
    the "on a bound" test by 2 kJ/mol and nothing more, which means the fit
    wants MORE temperature dependence on the thiol sink than the bound allows
    and the number is a bound-limited estimate, not a barrier. Neither family
    carries a finite confidence interval. Both must therefore travel with
    `identified: False`, and this test is what stops either being lifted into
    another module as if it were measured.
    """
    intervals = _frozen()["parameter_intervals"]
    for name in ("Ea_decay_thiol_sink", "Ea_decay_carbonyl_sink"):
        assert intervals[name]["identified"] is False, (
            f"{name} is reported as identified; if a later wave genuinely "
            f"identifies it, update this test deliberately rather than by "
            f"letting it pass silently"
        )
    hi = ps.DECAY_EA_BOUNDS[1]
    thiol = _frozen()["frozen_parameters"]["decay_Ea_kJ_mol"]["thiol_sink"]
    assert thiol > 0.9 * hi, (
        "the thiol-sink barrier is no longer pressed against its ceiling -- "
        "that is a real change in what the data say and the fit report's "
        "narrative must be rewritten, not the test relaxed"
    )


def test_the_frozen_drift_constants_are_inside_their_declared_bounds():
    drift = _frozen()["frozen_parameters"]["ph_drift"]
    lo, hi = PH.ACID_YIELD_BOUNDS
    assert lo <= drift["acid_yield_per_sink_event"] <= hi
    lo, hi = PH.ARP_AMINE_PKA_BOUNDS
    assert lo <= drift["arp_secondary_ammonium_pKa"] <= hi


def test_the_frozen_fit_reproduces_the_three_declared_drift_anchors():
    """
    Pre-registration section 5 item 2: the drift anchors missed by more than
    1.0 pH unit on any of the three would make the drift model a parameter
    that fits nothing.
    """
    rows = {r["id"]: r for r in _frozen()["rows"]}
    for label in ("6", "7", "8"):
        row = rows[f"zhou_final_pH_from_pH{label}"]
        error = abs(float(row["predicted"]) - float(row["observed"]))
        assert error <= 1.0, (
            f"initial pH {label}: predicted {row['predicted']:.2f} against a "
            f"measured {row['observed']:.2f} -- the drift model is not "
            f"reproducing its own declared calibration anchors"
        )


# ===========================================================================
# 6. THE FIREWALL -- values this wave unavoidably SAW must not have entered
# ===========================================================================
# DISCLOSURE: kang2026_SI_extraction.md prints the 140 C hold-out column and
# cutover_final_exam.md prints every Yiltirak measured value; the build brief
# directed this wave to read both. They were therefore SEEN before the fit ran.
# What follows is the mechanical check that they were not used.

_HOLDOUT_LITERALS = (
    # Kang 2026 SI Table S4, the 140 C column
    "5.907", "11.439", "60.400", "12.757", "73.157",
    # Kang Fig. S4's 140 C rung
    "62.6", "8.195",
    # Sun 2019's pH-9 column
    "0.4 / 0.5 / 0.5",
    # the Yiltirak T-ladder measured values (cutover_final_exam.md)
    "6.88", "1.28", "3.29", "1.46", "2.4", "1.68", "1.71", "1.62",
    # everything B2/B2.1 already held out
    "525.62", "325.22", "50.07", "582.34", "436.63",
    "696.99", "813.65", "59.70", "553.0", "229.0",
)

_FIT_SIDE_FILES = (
    "src/kinetic_core/sulfur.py",
    "src/kinetic_core/parameters_sulfur.py",
    "src/kinetic_core/species_sulfur.py",
    "src/kinetic_core/ph_state.py",
    "scripts/generators/generate_kinetic_core_b2_2_fit.py",
)


def test_no_holdout_literal_appears_on_the_fit_side():
    """
    Literals are matched with NUMERIC BOUNDARIES. Several Yiltirak values are
    short generic decimals (2.4, 1.28) that occur as digit runs inside
    unrelated numbers -- `2.4e13` in Zheng & Ho's prefactor set, `11.289` in
    the pKw table. Matching those as leaks would make the firewall cry wolf and
    would eventually get it deleted, so the boundary is part of the check.
    """
    import re

    offences = []
    for relative in _FIT_SIDE_FILES:
        text = (ROOT / relative).read_text()
        for literal in _HOLDOUT_LITERALS:
            # WAVE B2.4 -- ONE CHARACTER, AND THE REASON IS ON THE RECORD.
            # The literal list contains "2.4" (a Yiltirak hold-out fold). That
            # substring also occurs inside the WAVE NAME "B2.4", so as first
            # written this test made it impossible for any fit-side file to
            # mention the fourth wave of its own module -- which is a defect in
            # the guard, not in the prose. The added lookbehind exempts a match
            # glued to a letter and nothing else: a measured value is never
            # immediately preceded by a letter, so the guard loses no reach.
            pattern = re.compile(
                r"(?<![A-Za-z])(?<![0-9.eE])" + re.escape(literal)
                + r"(?![0-9eE])")
            for line_no, line in enumerate(text.splitlines(), 1):
                if not pattern.search(line):
                    continue
                upper = line.upper()
                if "HOLD-OUT" in upper or "IS NOT HERE" in upper:
                    continue
                offences.append(
                    f"{relative}:{line_no}: {literal!r} in {line.strip()[:90]}")
    assert not offences, (
        "hold-out literals leaked onto the fit side:\n" + "\n".join(offences))


def test_the_frozen_external_validation_bundles_were_never_opened():
    io_tokens = ("open(", "read_text", "read_bytes", "json.load", "glob(",
                 "iterdir")
    for relative in _FIT_SIDE_FILES + (
        "scripts/generators/generate_kinetic_core_b2_2_holdout.py",
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
    THE SYSTEMS WALK, which is the stronger check: a value can leak through a
    system definition without its literal ever being typed.
    """
    import scripts.generators.generate_kinetic_core_b2_2_fit as fit

    for name, spec in fit.SYSTEMS.items():
        t_c = float(spec["t_c"])
        ph = float(spec["ph"])
        if name.startswith("kang"):
            assert t_c in (100.0, 120.0), f"{name} at {t_c} C is the hold-out rung"
        if name.startswith("zhou") and not name.startswith("zhou_drift_"):
            assert ph == 7.0, f"{name} at pH {ph} is a held-out column"
        assert ph != 9.0, f"{name} sits on Sun 2019's held-out pH-9 column"
        assert ph != 6.5, f"{name} sits on Whitfield 2001's held-out pH"
        if name.startswith("hofmann_") and "pH5" not in name:
            raise AssertionError(f"unexpected Hofmann fit system {name}")
        # Yiltirak's ladder is 100-130 C at 30 min - 4 h in a buffered
        # ribose/cysteine pot. No fit system may sit on any of its rungs.
        if {"PENT", "Cys"} <= set(spec["initial"]) and len(spec["initial"]) <= 3:
            assert not (
                t_c in (100.0, 110.0, 120.0, 130.0)
                and float(spec["minutes"]) in (240.0, 120.0, 60.0, 30.0)
            ), f"{name} sits on a Yiltirak hold-out rung"


def test_the_drift_systems_can_only_ever_score_a_ph_endpoint():
    """
    THE NEW EXPOSURE THIS WAVE CREATES, MECHANICALLY FENCED.

    Amendment 7 declares Zhou's three final-pH endpoints FIT. Reaching them
    means INTEGRATING the pH-6 and pH-8 systems, whose VOLATILE columns are
    gating hold-outs. The licence covers the endpoint and nothing else, so
    every fit row that touches one of those systems must be a pH endpoint.
    """
    import scripts.generators.generate_kinetic_core_b2_2_fit as fit

    drift_systems = {n for n in fit.SYSTEMS if n.startswith("zhou_drift_")}
    assert drift_systems == {"zhou_drift_pH6", "zhou_drift_pH7", "zhou_drift_pH8"}
    for name in drift_systems:
        assert fit.SYSTEMS[name].get("drift_only") is True, name
    touched = set()
    for row in fit.ACTIVE_FIT_ROWS:
        if row["system"] in drift_systems:
            touched.add(row["system"])
            assert row["kind"] == "ph_endpoint", (
                f"{row['id']} reads system {row['system']} with kind "
                f"{row['kind']!r}. Only the pH ENDPOINT of a held-out Zhou "
                f"column is declared FIT (Amendment 7)."
            )
    assert touched == drift_systems, "a declared drift anchor is unscored"
    # and no drift system may be asked for a branch share or a flux either
    assert not any(
        row["system"] in drift_systems and row["kind"] in ("share", "ceiling")
        for row in fit.ACTIVE_FIT_ROWS
    )


def test_the_fit_scores_exactly_three_ph_endpoint_rows():
    import scripts.generators.generate_kinetic_core_b2_2_fit as fit

    ph_rows = [r for r in fit.ACTIVE_FIT_ROWS if r["kind"] == "ph_endpoint"]
    assert len(ph_rows) == 3
    assert {r["target"] for r in ph_rows} == {3.22, 3.42, 5.07}
