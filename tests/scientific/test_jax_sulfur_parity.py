"""The JAX sulfur lane (``src/kinetic_core_jax/sulfur_lane.py``) against the NumPy engine's sulfur lane, LIVE: every
engine number is computed in this process from the code on disk, so a change to ``src/kinetic_core/sulfur.py`` (or to
the ``parameters_sulfur`` / ``ph_state`` logic it uses) that the JAX lane does not follow fails here -- update
src/kinetic_core_jax/sulfur_lane.py in the same PR (docs/guides/A2_MODULE_CONTRACT.md sec. 6.2.2: the one-writer
rule's replacement while M1b(i) and M1c still edit the lane).

Skipped unless jax and diffrax import (the core env); CI runs it in the ``maillard-jax`` env (job "JAX sulfur lane
parity"). Scientific tier: ~10 - 20 min on one machine (25 panel pots, 6 draws, 4 mode pots, one reverse-mode
gradient), too slow for the unit tier.

THE TOLERANCE, from the solvers' own tolerances (``build_jax_sulfur_lane_module``'s docstring): both sides at the
engine's shipped rtol 1e-8 / atol 1e-14; local error per step <= rtol |y| + atol, no growth in a dissipative
mass-action network, at most 1e4 steps per run: <= 1e-4 relative per side, 8.7e-5 dex for the pair, so 1e-4 dex on
every species above atol / rtol = 1e-6 mmol/L at every pH node, 2e-10 mmol/L absolute below, 1e-4 on the pH (a log10
quantity), and the fixed point's sweep count EQUAL.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")

import jax.numpy as jnp  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"
UPDATE = "update src/kinetic_core_jax/sulfur_lane.py in the same PR"


@pytest.fixture(scope="module")
def g():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_sulfur_lane_module as gen

    return gen


@pytest.fixture(scope="module")
def h(g):
    return g.Harness()


@pytest.fixture(scope="module")
def panel(g, h):
    return [h.prepare(p) for p in g.panel_pots()]


def _check(g, h, pot, run, eng, what, theta=None, drift=None):
    """``Harness.matched`` (a sweep-count tie at the fixed point's tolerance is rerun at the engine's count) then the
    tolerance checks."""
    tol = {"rtol": g.FIT_TOL[0], "atol": g.FIT_TOL[1]}
    _, c = h.matched(pot, run, eng, theta, drift=drift, **tol)
    bad = g.within(c, **tol)
    assert not bad, f"{what}: {bad} -- {UPDATE}"
    return c


def test_the_jax_lane_has_the_engines_reactions_species_and_parameter_keys(h):
    from src.kinetic_core import engine, sulfur
    from src.kinetic_core.species_sulfur import SULFUR_STATE_KEYS
    from src.kinetic_core_jax import sulfur_structure as ss

    net = h.lane.network
    assert tuple(net.reactions) == tuple(sulfur.FULL_REACTION_KEYS), f"the reaction set differs -- {UPDATE}"
    assert tuple(net.species) == tuple(SULFUR_STATE_KEYS), f"the species set differs -- {UPDATE}"
    assert np.array_equal(net.stoichiometry, sulfur.SULFUR_STOICHIOMETRY), f"the stoichiometry differs -- {UPDATE}"
    generic = {r.parameter_key for r in sulfur.FULL_REACTIONS if r.key not in ss.HAND_PORTED}
    behind = {k.partition(":")[2] for k in h.lane.module.params if ":" in k}
    assert behind == generic, f"parameter keys: JAX {sorted(behind ^ generic)} differ -- {UPDATE}"
    assert set(engine.core_parameters(engine.SULFUR)) >= generic
    assert ss.classify()["custom"] == set(ss.HAND_PORTED), f"the hand-ported set differs -- {UPDATE}"
    assert h.lane.settings.coupling == "fixed_point_nodes" and h.lane.network.ph_mode == "inline"


def test_rate_constants_equal_the_engine(g, h):
    from src.kinetic_core import engine

    cases = [(engine.core_parameters(engine.SULFUR), engine.core_ph_drift())] + [g.draw(s)[:2] for s in g.DRAW_SEEDS[:1]]
    for params, drift in cases:
        out = g.rate_probe(h, params, drift)
        assert out["values"] > 30_000
        assert out["worst_relative"] < 1e-12, f"{out} -- {UPDATE}"
        assert out["nonzero_where_engine_zero"] == 0, f"{out} -- {UPDATE}"


def test_the_charge_balance_equals_ph_state(g, h):
    out = g.stage_balance(h)
    assert out["sid"] < 1e-9 and out["ionic_rel"] < 1e-8, out
    assert out["ph"] < 1e-8 and out["cool"] < 1e-8, out      # the engine's bisection stops at 1e-9 pH


def test_every_panel_pot_the_sulfur_lane_answers_matches_the_engine(g, h, panel):
    from src.kinetic_core import engine

    assert len(panel) >= 20, [p["id"] for p in panel]
    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    for pot in panel:
        run, _ = h.solve(pot)
        eng = g.engine_run(h.lane, pot, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        _check(g, h, pot, run, eng, pot["id"])
        assert not bool(np.any(np.asarray(run.failed)))
        e = h.lane.network.element_matrix
        np.testing.assert_allclose(e @ np.asarray(run.ys)[-1, -1], e @ pot["jax"]["y0"], rtol=1e-7, atol=1e-10)


def test_the_pot_prelude_is_the_engines(g, h, panel):
    """The JAX pot's state equals what engine._integrate_program integrates (OX, PROT_SS, OXR, OXV; MET / PRO / PYRL
    dropped): checked through the final state of the engine's own path on three pots."""
    from src.kinetic_core import engine

    params = engine.core_parameters(engine.SULFUR)
    for pot in panel[:: max(1, len(panel) // 3)][:3]:
        run, _ = h.solve(pot)
        final, _ = engine._integrate_program(engine.SULFUR, params, dict(pot["initial"]), pot["process"])
        ef = np.array([float(final.get(k, 0.0)) for k in h.lane.species])
        jf = np.asarray(run.ys)[-1, -1]
        m = np.maximum(ef, jf) >= g.FIT_TOL[1] / g.FIT_TOL[0]
        assert np.max(np.abs(np.log10(jf[m]) - np.log10(ef[m]))) < g.DEX_TOL, f"{pot['id']} -- {UPDATE}"


def test_random_parameter_draws_match_the_engine(g, h, panel):
    pots = [p for p in panel if p["id"] in g.DRAW_POTS]
    assert len(pots) == len(g.DRAW_POTS)
    for seed in g.DRAW_SEEDS:
        params, drift, coords = g.draw(seed)
        assert coords, "the draw moved nothing"
        theta = h.sl.theta_from(h.lane, params, drift)
        for pot in pots:
            run, _ = h.solve(pot, theta)
            eng = g.engine_run(h.lane, pot, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
            _check(g, h, pot, run, eng, f"draw {seed}, {pot['id']}", theta, drift)


def test_profile_staircase_clamped_and_inline_match_the_engine(g, h, panel):
    from src.kinetic_core import config, engine

    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    for pot in [h.prepare(p) for p in g.synthetic_pots()]:
        run, _ = h.solve(pot)
        _check(g, h, pot, run, g.engine_run(h.lane, pot, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1]), pot["id"])
    cfg = config.Config.from_globals().replace({"ph.coupling": "inline"})
    inline = g.Harness(cfg)
    for pot in [inline.prepare(p) for p in g.panel_pots() if p["id"] == g.DRAW_POTS[0]]:
        run, _ = inline.solve(pot)
        with config.applied(cfg):
            eng = g.engine_run(inline.lane, pot, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        c = g.compare(run, eng, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        c["jax_sweeps"] = c["engine_sweeps"]   # no fixed point under "inline"
        assert not g.within(c, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1]), f"inline {pot['id']}: {c} -- {UPDATE}"


def test_jax_grad_through_the_lane_is_finite_and_matches_central_differences(g, h, panel):
    """Reverse mode (the fixed point as a scan) on five coordinates, two of them the drift's (they enter only through
    the charge balance), at rtol 1e-10: against the port's own central differences (truncation O(step^2) plus solver
    noise rtol / step, both < 1e-5 here) and against forward mode (the same derivative another way)."""
    pot = next(p for p in panel if p["id"] == g.GRADIENT_POT)
    out = g.gradient_check(h, pot)
    assert out["finite"]
    grad, fd = np.asarray(out["grad"]), np.asarray(out["central_fd"])
    assert np.count_nonzero(np.abs(fd) > 1e-6) >= 3, out
    assert np.all(np.asarray(out["rel_grad_vs_fd"]) < 1e-4), out
    assert np.all(np.asarray(out["rel_grad_vs_jacfwd"]) < 1e-8), out


def test_unported_switches_and_pots_are_refused(g, h):
    from src.kinetic_core import config, engine, ph_state
    from src.kinetic_core_jax import sulfur_lane as sl
    from src.kinetic_core_jax.lane import UnsupportedConfiguration

    base = config.Config.from_globals()
    rivals = {
        "sulfur.ttca_formation": dataclasses.replace(base["sulfur.ttca_formation"], enabled=True),
        "sulfur.cys_amadori_route": dataclasses.replace(base["sulfur.cys_amadori_route"], enabled=True),
        "sulfur.strecker_sulfur": dataclasses.replace(base["sulfur.strecker_sulfur"], enabled=True),
        "ph.trunk_terms": True,
        "ph.trunk_solutes_in_balance": True,
        "ph.b51": dataclasses.replace(base["ph.b51"], ionic_strength=False),
        "sulfur.ddp_formation_ph": dataclasses.replace(base["sulfur.ddp_formation_ph"], law="none"),
    }
    for key, value in rivals.items():
        with pytest.raises(UnsupportedConfiguration, match=key.replace(".", r"\.")):
            sl.build(base.replace({key: value}))
    process = engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(120.0, 10.0), ph=6.0,
                                 ph_drift=ph_state.PhDrift(acid_yield=0.1, arp_amine_pka=8.0))
    with pytest.raises(UnsupportedConfiguration, match="ph_drift"):
        sl.make_pot(h.lane, {"Cys": 1.0}, process)
    with pytest.raises(UnsupportedConfiguration):
        sl.make_pot(h.lane, {"NOT_A_SPECIES": 1.0},
                    engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(120.0, 10.0), ph=6.0))
    import src.kinetic_core_jax.network as network

    with pytest.raises(UnsupportedConfiguration):
        network.integrate(h.lane.network, h.lane.theta0(), {"y0": jnp.zeros(h.lane.network.n_species),
                                                            "segments": np.array([[1.0, 100.0]])})
