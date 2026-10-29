"""src/kinetic_core/observations.py: each operator against a hand calculation (sub-model 1, 2026-10-29)."""
from __future__ import annotations

import pytest

from src.kinetic_core import observations as obs
from src.kinetic_core import ph_state
from src.kinetic_core.ph_state import BufferSpec, PhDrift

DRIFT = PhDrift(acid_yield=0.4, arp_amine_pka=8.0)
PHOSPHATE = BufferSpec(kind="phosphate", phosphate_mol_l=0.1, declared=True)


def test_base_to_hold_counts_a_strong_enough_acid_one_for_one():
    # 5 mM formic acid (pKa 3.75) at pH 7 is 99.9 % dissociated: a pH-stat must add 5 mM NaOH to hold the pH
    # (to ~1 %, the rest the ionic strength the formate adds, which moves phosphate's activity coefficients)
    uptake = obs.base_to_hold(7.0, {"FA": 5.0}, {"FA": 0.0}, DRIFT, PHOSPHATE, 298.15)
    assert uptake == pytest.approx(5.0 * (1.0 - 1.0 / (1.0 + 10 ** (7.0 - 3.75))), rel=0.02)


def test_base_to_hold_counts_only_the_dissociated_share_below_the_pka():
    # at pH 3.75 (formic's pKa) half of the added formic acid is dissociated
    uptake = obs.base_to_hold(3.75, {"FA": 4.0}, {"FA": 0.0}, DRIFT, BufferSpec(kind="none", declared=True), 298.15)
    assert uptake == pytest.approx(2.0, rel=0.1)


def test_the_cooled_ph_returns_the_label_at_the_label_sid():
    state = {"Gly": 200.0, "Glc": 200.0}
    sid = ph_state.strong_ion_difference(6.0, state, DRIFT, PHOSPHATE)
    assert obs.ph_cooled(state, sid, DRIFT, PHOSPHATE) == pytest.approx(6.0, abs=1e-6)


def test_a_dilution_of_one_is_the_sample_itself_and_a_buffer_moves_little_on_dilution():
    state = {"Cys": 400.0, "PENT": 400.0}
    sid = ph_state.strong_ion_difference(5.5, state, DRIFT, BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True))
    spec = BufferSpec(kind="phosphate", phosphate_mol_l=0.2, declared=True)
    assert obs.ph_diluted(state, sid, DRIFT, spec, 1.0) == pytest.approx(obs.ph_cooled(state, sid, DRIFT, spec), abs=1e-9)
    # a buffered sample keeps its ratio of acid to base on dilution: the pH moves only through activity, by tenths
    assert abs(obs.ph_diluted(state, sid, DRIFT, spec, 6.0) - 5.5) < 0.3
    with pytest.raises(ValueError):
        obs.ph_diluted(state, sid, DRIFT, spec, 0.5)


def test_a_diluted_buffer_divides_every_molarity():
    spec = obs.diluted(BufferSpec(kind="citrate_phosphate", phosphate_mol_l=0.2, citrate_mol_l=0.1, declared=True), 4.0)
    assert (spec.phosphate_mol_l, spec.citrate_mol_l, spec.kind) == (0.05, 0.025, "citrate_phosphate")


def test_formic_plus_acetic_counts_the_lumped_pool_at_the_calibrated_yield():
    assert obs.formic_plus_acetic({"FA": 1.0, "AA": 2.0, "ACID": 10.0}, DRIFT) == pytest.approx(7.0)


def test_amounts_convert_to_the_printed_unit():
    assert obs.amount_in("umol_L", "MGO", 0.5) == pytest.approx(500.0)
    assert obs.amount_in("ug_L", "MFT", 1e-3) == pytest.approx(114.16)     # 1 umol/L of MFT is 114.16 ug/L
    assert obs.amount_in("ug_kg", "FUR", 1.0) == pytest.approx(96080.0)
    assert obs.amount_in("percent_of_unheated", "FFT", 0.25, 1.0) == pytest.approx(25.0)
    with pytest.raises(ValueError):
        obs.amount_in("percent_of_unheated", "FFT", 0.25, 0.0)
