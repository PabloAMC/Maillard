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


def test_jax_grad_through_the_composed_network_needs_no_flag(join):
    gen, j = join
    out = gen.stage_gradient(j)
    assert out["finite"]
    fd = np.asarray(out["central_fd"])
    assert np.count_nonzero(np.abs(fd) > 1e-6) >= 3, out
    assert np.all(np.asarray(out["rel_grad_vs_fd"]) < 1e-4), out
    assert np.all(np.asarray(out["rel_grad_vs_jacfwd"]) < 1e-8), out
