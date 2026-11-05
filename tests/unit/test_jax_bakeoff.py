"""src/kinetic_core_jax (the JAX bake-off port of the sulfur lane, STATUS_AND_PLAN W4 step H2; 2026-11-04).

Skipped unless jax and diffrax import, so the core env's CI never runs (or breaks on) them; run them in the separate
``maillard-jax`` env (environment-jax.yml). They pin what results/validation/jax_bakeoff_sulfur_lane.md claims:

* every rate constant equals the engine's ``sulfur_rate_constants_at`` (R3En configuration) to 1e-12 relative on a
  (pH, label pH, I) grid, at the fit's optimum and at a perturbed point;
* the charge balance gives the engine's (pH, I), hot and cooled, and its SID;
* the pH solve's implicit-function derivative equals a finite difference;
* the port refuses an engine switch it does not implement;
* (slow) a 60 min integration matches the engine's tight-tolerance trajectory, and forward-mode gradients equal
  central differences of the port.
"""
from __future__ import annotations

import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")

import jax.numpy as jnp  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"
REFERENCE = ROOT / "results" / "validation" / "jax_bakeoff_sulfur_lane_reference.json"


@pytest.fixture(scope="module")
def setup():
    """The port built in the R3En engine state, and the engine's own numbers computed in that state, then every
    switch and script global restored."""
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_bakeoff as bk
    import fit_joint_deoxyosone_thiol as joint
    from src.kinetic_core import config, engine, ph_state, sulfur

    reference = json.loads(REFERENCE.read_text())
    with pytest.MonkeyPatch.context() as mp:
        for name in ("DDP_LAW_STATIC", "TDP_EXIT_STATIC", "STRECKER_SHARE", "EA_STATIC"):
            mp.setattr(joint, name, getattr(joint, name))
        with config.applied(config.Config.from_globals()):
            lane, pot, job = bk.setup_pot(reference)
            params = engine.core_parameters(engine.SULFUR)
            for key, k_ref in job["k_ref"].items():
                params[key] = dataclasses.replace(params[key], k_ref=float(k_ref))
            t_k = job["t_c"] + 273.15
            grid = [(ph, phl, i) for ph in (4.1, 5.0, 6.7, 8.8) for phl in (4.6, 5.2, 7.1) for i in (0.0, 0.3)]
            k_engine = [sulfur.rate_vector(sulfur.sulfur_rate_constants_at(params, t_k, ph, ph_label=phl,
                                                                           ionic_strength_m=i)) for ph, phl, i in grid]
            import build_submodel_probe as probe

            drift, spec = probe._drift(job), probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"])
            sid = reference["trajectory_tight_tol"]["sid_mmol_L"]
            states = reference["trajectory_tight_tol"]["states"][::9]
            ph_engine = []
            for row in states:
                state = {k: float(row[sulfur.SULFUR_INDEX[k]]) for k in sulfur.SULFUR_STATE_KEYS}
                ph_engine.append((ph_state.ph_and_ionic_strength_of_state(state, t_k, sid, drift, spec),
                                  ph_state.ph_and_ionic_strength_of_state(state, 298.15, sid, drift, spec)))
            yield {"bk": bk, "reference": reference, "lane": lane, "pot": pot, "job": job, "grid": grid,
                   "k_engine": k_engine, "states": states, "ph_engine": ph_engine, "sid": sid}


def test_rate_constants_equal_the_engine(setup):
    from src.kinetic_core_jax import model

    lane, pot = setup["lane"], setup["pot"]
    theta = jnp.asarray(setup["reference"]["theta"])
    f = jax.jit(lambda ph, phl, i: model.rate_constants(theta, pot.data["t_k"], ph, phl, i, lane))
    for (ph, phl, i), k_ref in zip(setup["grid"], setup["k_engine"]):
        k = np.asarray(f(ph, phl, i))
        np.testing.assert_allclose(k, k_ref, rtol=1e-12, atol=0.0, err_msg=f"pH {ph}, label {phl}, I {i}")


def test_rate_constants_at_a_perturbed_point(setup):
    """The reference file's perturbed-theta probe (written by the core env)."""
    from src.kinetic_core_jax import model

    lane, pot = setup["lane"], setup["pot"]
    probe = setup["reference"]["rate_constant_probe"]["theta_perturbed"]
    theta = jnp.asarray(probe["theta"])
    f = jax.jit(lambda ph, phl, i: model.rate_constants(theta, pot.data["t_k"], ph, phl, i, lane))
    for (ph, phl, i), k_ref in list(zip(probe["grid_ph_phlabel_I"], probe["k"]))[::7]:
        np.testing.assert_allclose(np.asarray(f(ph, phl, i)), k_ref, rtol=1e-12, atol=0.0)


def test_charge_balance_equals_the_engine(setup):
    pot = setup["pot"]
    assert float(pot.data["sid"]) == pytest.approx(setup["sid"], rel=1e-12)
    ph_fn = jax.jit(pot.ph_state)
    for row, ((ph_hot, i_hot), (ph_cool, i_cool)) in zip(setup["states"], setup["ph_engine"]):
        hot, ionic, cool, ionic_cool = map(float, ph_fn(jnp.asarray(row)))
        # the engine's bisection stops at 1e-9 pH and its I fixed point at 1e-9 M; the port's Newton is exact
        assert hot == pytest.approx(ph_hot, abs=2e-9)
        assert cool == pytest.approx(ph_cool, abs=2e-9)
        assert ionic == pytest.approx(i_hot, rel=1e-8)
        assert ionic_cool == pytest.approx(i_cool, rel=1e-8)


def test_ph_derivative_is_the_implicit_function(setup):
    """d pH / d inventory by the custom JVP equals a central difference of the forward solve."""
    from src.kinetic_core_jax import model

    pot, lane = setup["pot"], setup["lane"]
    y = jnp.asarray(setup["states"][3])
    conc = model.inventory(y, lane, pot.data["buffer"])
    solve = pot.solve_ph
    t_k, sid = pot.data["t_k"], pot.data["sid"]
    jac = np.asarray(jax.jacfwd(lambda c: solve(t_k, sid, c))(conc))
    rev = np.asarray(jax.jacrev(lambda c: solve(t_k, sid, c))(conc))
    np.testing.assert_allclose(rev, jac, rtol=1e-10, atol=1e-14)
    for g in range(len(lane.group_names)):
        h = 1e-4 * max(float(conc[g]), 1.0)
        up = np.asarray(solve(t_k, sid, conc.at[g].add(h)))
        dn = np.asarray(solve(t_k, sid, conc.at[g].add(-h)))
        np.testing.assert_allclose(jac[:, g], (up - dn) / (2 * h), rtol=1e-5, atol=1e-12,
                                   err_msg=lane.group_names[g])


def test_the_port_refuses_a_switch_it_does_not_implement(setup):
    from src.kinetic_core import config, engine
    from src.kinetic_core import parameters_sulfur as ps
    from src.kinetic_core_jax import lane as lane_mod

    import build_submodel_probe as probe

    job = setup["job"]
    with config.applied(config.Config.from_globals().replace({"ph.coupling": "inline"})):
        cfg = config.Config.from_globals()
        with config.applied(cfg.replace({"sulfur.tdp_own_decay": ps.TdpOwnDecay(enabled=True)})):
            with pytest.raises(lane_mod.UnsupportedConfiguration):
                lane_mod.build(engine.core_parameters(engine.SULFUR), setup["reference"]["keys"],
                               probe.buffer_of({"buffer": job["buffer"]}, job["label_ph"]), probe._drift(job))


@pytest.mark.slow
def test_a_short_integration_matches_the_engine(setup):
    pot, reference = setup["pot"], setup["reference"]
    tight = reference["trajectory_tight_tol"]
    times = [t for t in tight["times"] if t <= 60.0]
    ys = np.asarray(jax.jit(lambda th: pot.solve(th, times, rtol=1e-10, atol=1e-16).ys)(jnp.asarray(reference["theta"])))
    ref = np.asarray(tight["states"])[: len(times)]
    for species in ("MFT", "FFT", "Cys", "NF", "TTCA", "H2S"):
        i = reference["species"].index(species)
        for r, t in enumerate(times):
            if t >= 1.0:
                assert abs(math.log10(ys[r, i]) - math.log10(ref[r, i])) < 1e-6, (species, t)


@pytest.mark.slow
def test_forward_mode_gradients_equal_central_differences(setup):
    import diffrax

    pot, reference = setup["pot"], setup["reference"]
    species = reference["species"]
    idx = [species.index("MFT"), species.index("FFT")]
    times = [0.0, 30.0, 60.0]
    theta = jnp.asarray(reference["theta"])

    def obs(th, adjoint=None):
        ys = pot.solve(th, times, rtol=1e-11, atol=1e-17, adjoint=adjoint).ys
        return jnp.log10(ys[1:, idx]).ravel()

    jac = np.asarray(jax.jit(jax.jacfwd(lambda th: obs(th, diffrax.ForwardMode())))(theta))
    f = jax.jit(obs)
    keys = reference["keys"]
    for c in (keys.index("offset_k_dpo_ddp"), keys.index("offset_k_nf_mft"), keys.index("offset_fft_formation"),
              keys.index("mgo_share_nf")):
        h = 1e-4
        e = jnp.zeros(len(keys)).at[c].set(h)
        fd = (np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h)
        np.testing.assert_allclose(jac[:, c], fd, rtol=1e-4, atol=1e-6, err_msg=keys[c])
