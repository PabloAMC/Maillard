"""The one JAX network with the sulfur lane in it (A2 contract v2.3), LIVE against the engine: the composed lane
(``sulfur_lane.build_composed``: three modules, the pH solved by ``network.integrate``) on two panel pots at the
sulfur parity test's bound, and ``jax.grad`` through it with the default adjoint (no ``reverse`` flag) against central
differences. The full record (five pots, the fit driver, the density bridge, timings) is
results/validation/jax_network_join.md, written by scripts/generators/build_jax_network_join.py.

Skipped unless jax and diffrax import; CI runs it in the ``maillard-jax`` env beside the sulfur parity test.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"


@pytest.fixture(scope="module")
def join():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_network_join as gen

    return gen, gen.Join()


def test_the_composed_network_matches_the_engine(join):
    from src.kinetic_core import engine

    gen, j = join
    g = gen.g
    params, drift = engine.core_parameters(engine.SULFUR), engine.core_ph_drift()
    for pid in ("thiamine_cys_glucose_120C_Bolton1994", "mp_holdout_arp_glutathione_120C_30min_pH75_Zhou2023c"):
        pot = j.pot(j.composed, pid)
        run, _ = j.composed.solve(pot)
        eng = g.engine_run(j.composed.lane, pot, params, drift, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        _, c = j.composed.matched(pot, run, eng)
        bad = g.within(c, rtol=g.FIT_TOL[0], atol=g.FIT_TOL[1])
        assert not bad, f"{pid}: {bad}"


def test_an_unconverged_fixed_point_is_a_failed_pot_and_sweep_counts_are_logged(join):
    """The engine marks a run whose last fixed-point residual exceeds PH_FIXED_POINT_FLAG as not converged; the
    network reports it as failed, the fit driver scores it as failed, and logs every pot's sweep counts."""
    import dataclasses

    from src.kinetic_core_jax import fit, network

    gen, j = join
    lane = j.one.lane
    pot = j.pot(j.one, "thiamine_cys_glucose_120C_Bolton1994")["jax"]
    strict = network.compose([dataclasses.replace(lane.module, ph_model=dataclasses.replace(
        lane.module.ph_model, flag=1e-9))])
    arrays, mode = network.pot_arrays(lane.network, pot), network.pot_mode(lane.network, pot)
    ok = jax.jit(lambda th, a: network.integrate(lane.network, th, a, mode=mode))(lane.theta0(), arrays)
    bad = jax.jit(lambda th, a: network.integrate(strict, th, a, mode=mode))(strict.theta0(), arrays)
    assert not bool(np.any(np.asarray(ok.failed))) and not bool(np.any(np.asarray(ok.ph.unconverged)))
    assert bool(np.all(np.asarray(bad.ph.unconverged))) and bool(np.all(np.asarray(bad.failed)))
    assert not bool(np.any(np.asarray(bad.ph.failed)))          # the integrator itself succeeded
    key = ("sulfur", "acid_yield")
    obs = (fit.Observation(pot=0, row=network.segment_end_rows(lane.network, pot)[-1], observable="mmol_per_L:MFT",
                           value=-3.0, sigma=0.1),)
    for net, failed in ((lane.network, []), (strict, [0])):
        comp = fit.Compiled(fit.JaxProblem(net=net, keys=(key,), kinds=("acid_yield",), lower=np.array([0.0]),
                                           upper=np.array([1.0]), pots=(pot,), observations=obs, batch="loop"))
        x = np.array([float(net.theta0()["sulfur"]["acid_yield"])])
        r, got = comp.residual(x)
        assert got == failed
        assert comp.sweep_log[0] == [tuple(int(s) for s in np.asarray(ok.ph.sweeps))]
        assert comp.sweep_changes() == {}


def test_jax_grad_through_the_composed_network_needs_no_flag(join):
    gen, j = join
    out = gen.stage_gradient(j)
    assert out["finite"]
    fd = np.asarray(out["central_fd"])
    assert np.count_nonzero(np.abs(fd) > 1e-6) >= 3, out
    # A central difference's error is absolute (solver tolerance / step), the same on every coordinate, so a component
    # ~1000x below the largest cannot be held to 1e-4 of itself (k_hmf_cys: |fd| 7e-4 against 0.98; its relative error
    # moved 3e-5 -> 2e-4 between runners). Hold each to 1e-4 of max(|fd_i|, 1 % of the largest); the exact check
    # against forward mode below stays per-component at 1e-8.
    gv = np.asarray(out["grad"])
    floor = np.maximum(np.abs(fd), 1e-2 * np.max(np.abs(fd)))
    assert np.all(np.abs(gv - fd) / floor < 1e-4), out
    assert np.all(np.asarray(out["rel_grad_vs_jacfwd"]) < 1e-8), out
