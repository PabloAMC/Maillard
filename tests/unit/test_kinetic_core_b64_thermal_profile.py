"""
tests/unit/test_kinetic_core_b64_thermal_profile.py

WAVE B64 -- the thermal history as the source records it (2026-10-28).

A pot whose source records a heat-up is a piecewise-linear PROFILE (``ThermalProgram.from_profile``).
The single-temperature integrators read its staircase; the sulfur lane integrates the profile
itself, with its rate constants at T(t) on its pH nodes. These tests check that (i) a program built
from holds runs exactly as before, (ii) a flat profile is the isothermal run, and (iii) a ramp
reproduces the exact Arrhenius integral of a first-order step and agrees with a fine staircase.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from src.kinetic_core import engine as E
from src.kinetic_core import ph_state as PH
from src.kinetic_core import sulfur as S

T145 = 418.15
_PHOSPHATE = PH.BufferSpec(kind="phosphate", phosphate_mol_l=0.5, declared=True, source="test")
_CLAMPED = PH.BufferSpec(kind="clamped")
#: the Hofmann 1995 reference pot's recorded heating: 20 -> 145 C within 20 min
_RAMP_K = (np.array([0.0, 20.0]), np.array([20.0, 145.0]) + 273.15)


@pytest.fixture(scope="module")
def parameters():
    return E.core_parameters(E.SULFUR)


def _state(run):
    return {k: float(run.concentrations[-1, S.SULFUR_INDEX[k]]) for k in S.SULFUR_STATE_KEYS}


def test_staircase_steps_and_duration():
    program = E.ThermalProgram.from_profile(((0.0, 20.0), (20.0, 145.0), (40.0, 145.0)))
    assert program.total_minutes == pytest.approx(40.0)
    assert program.peak_temperature_c == 145.0 and program.min_temperature_c == 20.0
    temps = [c for _, c in program.segments]
    assert max(abs(b - a) for a, b in zip(temps, temps[1:])) <= E.RAMP_STEP_K + 1e-9
    # the ramp's staircase sits at mid temperatures; the hold is one step at 145 C
    assert temps[0] == pytest.approx(20.0 + 125.0 / 63 / 2)
    assert program.segments[-1] == (20.0, 145.0)
    assert program.temperature_c_at([0.0, 10.0, 30.0]).tolist() == pytest.approx([20.0, 82.5, 145.0])
    assert E.ThermalProgram.isothermal(145.0, 20.0).profile is None


def test_flat_profile_is_the_isothermal_run(parameters):
    charge = {"PENT": 100.0, "Cys": 33.0}
    common = dict(ph=5.0, buffer_spec=_PHOSPHATE, ph_drift=E.core_ph_drift(), rtol=1e-8, atol=1e-14)
    plain = S.integrate_sulfur(parameters, T145, charge, [0.0, 20.0], **common)
    flat = S.integrate_sulfur(parameters, T145, charge, [0.0, 20.0], **common,
                              temperature_profile=(np.array([0.0, 20.0]), np.array([T145, T145])))
    assert np.array_equal(plain.concentrations, flat.concentrations)


def test_ramp_reproduces_the_arrhenius_integral(parameters):
    """Furfural alone, pH clamped: its loss is first order, so FUR/FUR0 = exp(-int k(T(t)) dt)."""
    first_order = [r for r in S.FULL_REACTIONS if dict(r.reactants) == {"FUR": 1}]
    assert first_order, "the network has a first-order furfural sink"

    def k_total(temperature_k):
        k = S.sulfur_rate_constants_at(parameters, float(temperature_k), 5.0)
        return sum(float(k[r.key]) for r in first_order)

    t = np.linspace(0.0, 20.0, 4001)
    exponent = np.trapezoid([k_total(T) for T in np.interp(t, *_RAMP_K)], t)
    run = S.integrate_sulfur(parameters, T145, {"FUR": 10.0}, [0.0, 20.0], ph=5.0,
                             buffer_spec=_CLAMPED, rtol=1e-10, atol=1e-14,
                             temperature_profile=_RAMP_K)
    assert _state(run)["FUR"] / 10.0 == pytest.approx(math.exp(-exponent), rel=2e-3)
    # and it is NOT the isothermal answer, which the ramp exposes far less
    hot = S.integrate_sulfur(parameters, T145, {"FUR": 10.0}, [0.0, 20.0], ph=5.0,
                             buffer_spec=_CLAMPED, rtol=1e-10, atol=1e-14)
    assert _state(hot)["FUR"] < _state(run)["FUR"]


def test_ramp_agrees_with_a_fine_staircase(parameters):
    """A full pentose + cysteine pot, pH clamped: one profile run against 0.25 K isothermal steps."""
    charge = {"PENT": 100.0, "Cys": 33.0}
    run = S.integrate_sulfur(parameters, T145, charge, [0.0, 20.0], ph=5.0, buffer_spec=_CLAMPED,
                             rtol=1e-9, atol=1e-14, temperature_profile=_RAMP_K)
    stairs = E.ThermalProgram.from_profile(((0.0, 20.0), (20.0, 145.0)), max_step_k=0.25)
    state = dict(charge)
    for duration, temperature_c in stairs.segments:
        step = S.integrate_sulfur(parameters, temperature_c + 273.15, state, [0.0, duration], ph=5.0,
                                  buffer_spec=_CLAMPED, rtol=1e-9, atol=1e-14)
        state = _state(step)
    ramp = _state(run)
    for key in ("MFT", "FFT", "NF", "FUR", "Cys", "PENT"):
        if state[key] > 1e-9:
            assert ramp[key] == pytest.approx(state[key], rel=1e-2), key


def test_engine_routes_a_profile_through_one_sulfur_integration(parameters):
    program = E.ThermalProgram.from_profile(((0.0, 20.0), (20.0, 145.0)))
    process = E.ProcessSpec(thermal=program, ph=5.0, buffer=_PHOSPHATE)
    final, metadata = E._integrate_program(E.SULFUR, parameters, {"PENT": 100.0, "Cys": 33.0}, process)
    assert len(metadata["segments"]) == 1 and metadata["segments"][0]["temperature_C"] == 145.0
    hot = E.ProcessSpec(thermal=E.ThermalProgram.isothermal(145.0, 20.0), ph=5.0, buffer=_PHOSPHATE)
    final_hot, _ = E._integrate_program(E.SULFUR, parameters, {"PENT": 100.0, "Cys": 33.0}, hot)
    assert final["Cys"] > final_hot["Cys"]


def test_heating_block_conventions():
    from src.kinetic_core.thermal_history import HeatingUnstated, profile_points

    assert profile_points(None, 145.0, 20.0) is None
    assert profile_points({"time_convention": "isothermal"}, 145.0, 20.0) is None
    assert profile_points({"time_convention": "ramp_only"}, 145.0, 20.0) == ((0.0, 20.0), (20.0, 145.0))
    assert profile_points({"time_convention": "ramp_plus_hold", "ramp_min": 10, "hold_min": 20},
                          130.0, 20.0) == ((0.0, 20.0), (10.0, 130.0), (30.0, 130.0))
    assert profile_points({"time_convention": "hold_after_reaching", "ramp_min": 3}, 120.0, 60.0) \
        == ((0.0, 20.0), (3.0, 120.0), (63.0, 120.0))
    assert profile_points({"time_convention": "includes_come_up", "ramp_min": 3, "cool_min": 2,
                           "cool_to_C": 40}, 120.0, 60.0) \
        == ((0.0, 20.0), (3.0, 120.0), (60.0, 120.0), (62.0, 40.0))
    with pytest.raises(HeatingUnstated):
        profile_points({"time_convention": "hold_after_reaching"}, 120.0, 60.0)


def test_panel_reads_heating_only_when_asked():
    from src.kinetic_core import panel

    bench = {"conditions": {"temp_C": 145.0, "time_min": 20.0, "heating": {"time_convention": "ramp_only"}}}
    saved = panel.THERMAL_HISTORY
    try:
        assert panel.bundle_thermal_program(bench).profile is None
        panel.THERMAL_HISTORY = "as_recorded"
        assert panel.bundle_thermal_program(bench).profile == ((0.0, 20.0), (20.0, 145.0))
    finally:
        panel.THERMAL_HISTORY = saved


def test_fast_charge_residual_is_bit_identical():
    """C5: the bisection's precomputed residual equals ``_charge_residual`` exactly."""
    rng = np.random.default_rng(64)
    state = {"Cys": 33.0, "PENT": 100.0, "ACID": 12.0, "H2S": 0.4, "TTCA": 5.0}
    inventory = PH.titratable_inventory(state, E.core_ph_drift(), _PHOSPHATE)
    for temperature_k in (298.15, 373.15, 418.15):
        for ionic in (0.0, 0.31):
            fast = PH._residual_at_fixed_conditions(temperature_k, -40.0, inventory, ionic)
            for ph in rng.uniform(-1.0, 15.0, 50):
                assert fast(ph) == PH._charge_residual(ph, temperature_k, -40.0, inventory, ionic)


def test_inline_ph_coupling_matches_the_converged_fixed_point(parameters):
    """C5: the inline charge balance is the limit of B60's node fixed point (81 nodes, rtol 1e-8)."""
    charge = {"TTCA": 10.0}
    common = dict(ph=7.0, buffer_spec=PH.DEFAULT_BUFFER, ph_drift=E.core_ph_drift(), ph_iterations=4,
                  rtol=1e-8, atol=1e-14)
    saved = S.PH_COUPLING
    try:
        S.PH_COUPLING = "fixed_point_nodes"
        nodes = S.integrate_sulfur(parameters, 373.15, charge, [0.0, 120.0], ph_nodes=81, **common)
        S.PH_COUPLING = "inline"
        inline = S.integrate_sulfur(parameters, 373.15, charge, [0.0, 120.0], **common)
    finally:
        S.PH_COUPLING = saved
    assert inline.metadata["ph_coupling"] == "inline" and inline.metadata["ph_fixed_point_sweeps"] == 0
    assert inline.metadata["ph_final_cooled"] == pytest.approx(nodes.metadata["ph_final_cooled"], abs=2e-3)
    for key in ("MFT", "FFT", "H2S", "Cys", "TTCA", "FUR", "NF"):
        a, b = nodes.final(key), inline.final(key)
        if max(a, b) > 1e-6:
            assert abs(math.log10((b + 1e-9) / (a + 1e-9))) < 5e-3, key


def test_compiled_charge_balance_agrees_with_the_python_solve():
    """C5: ph_kernel (numba where installed) and the Python warm solve find the same pH and I."""
    rng = np.random.default_rng(5)
    for _ in range(40):
        state = {"Cys": rng.uniform(0, 60), "ACID": rng.uniform(0, 30), "TTCA": rng.uniform(0, 20),
                 "ARP": rng.uniform(0, 20), "FA": rng.uniform(0, 5), "CBX": rng.uniform(0, 5)}
        buffer = _PHOSPHATE if rng.uniform() < 0.5 else PH.DEFAULT_BUFFER
        inventory = PH.titratable_inventory(state, E.core_ph_drift(), buffer)
        for temperature_k in (298.15, 373.15, 418.15):
            sid = float(rng.uniform(-80.0, 40.0))
            guess = float(rng.uniform(2.0, 9.0))
            a = PH.solve_ph_and_ionic_strength_warm(temperature_k, sid, inventory, guess, 0.1)
            b = PH.solve_ph_and_ionic_strength_compiled(temperature_k, sid, inventory, guess, 0.1)
            assert b[0] == pytest.approx(a[0], abs=1e-8) and b[1] == pytest.approx(a[1], rel=1e-7, abs=1e-12)


def test_b59a_strecker_sulfur_split_is_conservative(parameters):
    """B59a V1/V2: off, the new channel is zero; on, the two channels sum to the old Strecker constant."""
    import dataclasses

    from src.kinetic_core import parameters_sulfur as ps

    off = S.sulfur_rate_constants_at(parameters, T145, 6.0)
    assert off["r_dpo_ddp_h2s"] == 0.0
    saved = ps.STRECKER_SULFUR
    try:
        for share in saved.alternatives:
            ps.STRECKER_SULFUR = dataclasses.replace(saved, enabled=True, h2s_share=share)
            on = S.sulfur_rate_constants_at(parameters, T145, 6.0)
            assert on["r_dpo_ddp"] + on["r_dpo_ddp_h2s"] == pytest.approx(off["r_dpo_ddp"], rel=1e-12)
            assert on["r_dpo_ddp_h2s"] == pytest.approx(share * off["r_dpo_ddp"], rel=1e-12)
            run = S.integrate_sulfur(parameters, T145, {"DPO": 10.0, "Cys": 10.0}, [0.0, 20.0], ph=6.0,
                                     buffer_spec=_CLAMPED, rtol=1e-9, atol=1e-14)
            for element in ("carbon", "nitrogen", "sulfur"):
                assert run.metadata[f"{element}_closure"]["max_relative_drift"] < 1e-9, element
            assert run.final("H2S") > 0.0
    finally:
        ps.STRECKER_SULFUR = saved
