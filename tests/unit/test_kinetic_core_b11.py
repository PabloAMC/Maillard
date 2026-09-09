"""Wave B11 (2026-09-07): oxygen as a two-pool input on the sulfur lane -- the STRUCTURE, inert
until a B11 report ships its consumers."""
from __future__ import annotations

import pytest

from src.kinetic_core import engine, panel, sulfur, species_sulfur as ss
from src.kinetic_core.engine import SULFUR, core_parameters, oxygen_reservoir_units, predict
from src.kinetic_core.parameters_sulfur import (
    MEASURED_SULFUR, OX_RESERVOIR_DEFAULT_UNITS, OX_SAT_MMOL_L, OXYGEN_KEYS, oxygen_parameters,
)
from src import data_paths


def test_the_structure_is_in_the_network_and_balances():
    keys = set(sulfur.FULL_REACTION_KEYS)
    assert {"ox_supply", "ch_cys_ox", "ch_red_ox_dpo", "ch_red_ox_nf"} <= keys
    assert "OXR" in ss.SULFUR_INDEX and "OXV" in ss.SULFUR_INDEX
    assert "OXR" in ss.SITE_POOLS and "OXV" in ss.SITE_POOLS
    dimer = next(r for r in sulfur.FULL_REACTIONS if r.key == "ch_dimer_mft")
    assert dimer.products.get("OXV") == 1
    sulfur.validate_sulfur_balance()          # atoms; raises on a bad step


def test_the_defaults_are_inert_so_every_earlier_wave_reproduces():
    for key in OXYGEN_KEYS:
        assert key in MEASURED_SULFUR
    assert MEASURED_SULFUR["k_cys_ox"].k_ref == 0.0 and MEASURED_SULFUR["k_red_ox"].k_ref == 0.0
    assert MEASURED_SULFUR["k_ox_supply"].k_ref == 0.0          # switched on with the consumers
    assert oxygen_parameters(k_cys_ox=1e-3)["k_ox_supply"].k_ref > 0.0
    built = core_parameters(SULFUR)
    frozen = engine.frozen_parameters(SULFUR)
    if "oxygen" not in frozen:
        assert built["k_cys_ox"].k_ref == 0.0 and built["k_red_ox"].k_ref == 0.0


def test_a_hofmann_pot_reproduces_the_pre_b11_number_when_the_consumers_are_zero():
    bench = panel.load_bundle(data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json")
    spec = panel.core_spec(bench, use_buffer=True)
    assert spec.process.vessel is not None and spec.process.vessel.fill_mL == 100.0
    decl = engine.declare_envelope(spec, list(panel.bundle_targets(bench)))
    operative = core_parameters(SULFUR)
    operative = {**operative, **oxygen_parameters(0.0, 0.0)}
    with_reservoir, meta = engine._integrate_program(SULFUR, operative, dict(decl.mapped_precursors), spec.process)
    no_reservoir, _ = engine._integrate_program(SULFUR, operative, {**decl.mapped_precursors, "OXR": 0.0}, spec.process)
    for key in ("MFT", "FFT", "OX"):
        assert with_reservoir[key] == pytest.approx(no_reservoir[key], rel=1e-9), key
    assert with_reservoir["OX"] == pytest.approx(1.0, rel=1e-3)   # only the dimers' own draw, not refilled
    # inert consumers: the reservoir changes no rate, so the vessel is not an extrapolation flag
    assert not any("OXYGEN RESERVOIR (B11)" in w for w in decl.warnings)


def test_the_reservoir_arithmetic_and_defaults():
    bench = panel.load_bundle(data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_ribose_cysteine_buffer_100C_4h_Yiltirak2026.json")
    spec = panel.core_spec(bench, use_buffer=True)
    units, basis = oxygen_reservoir_units(spec.process)
    assert units == pytest.approx(49.6 / OX_SAT_MMOL_L, rel=0.02)     # 0.149 mmol over 3 mL
    assert "headspace" in basis
    hof = panel.core_spec(panel.load_bundle(data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json"))
    # Hofmann's own vessel: 100 mL headspace of air (0.871 mmol) + dissolved, over 0.1 L, / 0.3
    assert oxygen_reservoir_units(hof.process)[0] == pytest.approx((0.8706 + 0.027) / 0.1 / OX_SAT_MMOL_L, rel=0.01)
    assert OX_RESERVOIR_DEFAULT_UNITS == pytest.approx(8.7 / OX_SAT_MMOL_L)
    bare = engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(145.0, 20.0))
    units, basis = oxygen_reservoir_units(bare)
    assert units == OX_RESERVOIR_DEFAULT_UNITS and "default" in basis


def test_with_consumers_on_the_reservoir_depletes_and_the_thiols_respond():
    bench = panel.load_bundle(data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json")
    spec = panel.core_spec(bench, use_buffer=True)
    decl = engine.declare_envelope(spec, list(panel.bundle_targets(bench)))
    base = core_parameters(SULFUR)
    # The thiolate factor is 1 at bench pH 5 but the in-situ pH of a 0.5 M phosphate pot at
    # 145 C sits lower, so a large constant is used here to prove the MECHANICS (depletion,
    # then the thiols' response), not a plausible rate: 33 mM cysteine at 5 per unit per
    # minute drains the 30-unit reservoir inside the 20 min cook.
    fast = {**base, **oxygen_parameters(k_cys_ox=50.0, k_red_ox=0.0)}
    state, _ = engine._integrate_program(SULFUR, fast, dict(decl.mapped_precursors), spec.process)
    assert state["OXR"] < 1.0 and state["OX"] < 0.6          # the reservoir is spent (probe: OX 0.52)
    slow = {**base, **oxygen_parameters(k_cys_ox=1e-5, k_red_ox=0.0)}
    state2, _ = engine._integrate_program(SULFUR, slow, dict(decl.mapped_precursors), spec.process)
    assert state2["OX"] == pytest.approx(1.0, abs=0.05)
    # the thiols respond -- and in the direction the chemistry sets: autoxidation removes the
    # sulfur SOURCE (cysteine 0.04 mM left of 33), so the thiols fall with heavy consumption
    assert state["Cys"] < 0.1 * state2["Cys"]
    assert state["FFT"] < state2["FFT"]


def test_a_draw_can_move_the_consumers_and_the_reservoir():
    run0 = predict(_pot(), ["2-furfurylthiol"])
    draw = engine.CoreDraw(maillard={"oxygen": {"k_cys_ox": 0.1, "k_red_ox": 0.0}}, oxygen_reservoir_scale=0.1)
    run1 = predict(_pot(), ["2-furfurylthiol"], draw=draw)
    assert run1.require("2-furfurylthiol") != run0.require("2-furfurylthiol")


def _pot():
    return engine.FormulationSpec(name="pot", precursors={"D-Ribose": 100.0, "L-Cysteine": 33.0},
                                  process=engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(145.0, 20.0), ph=5.0))
