"""src/kinetic_core_jax: the A2 network infrastructure (contract v2.2) -- the forced-input driver (``forcing.py``), the
fit driver (``fit.py``), the declared-competition and titratable-table hooks (``network.compose``) and the bridge to
M1.1's predictive densities (``densities.py``).

Skipped unless jax and diffrax import, so the core env's CI never runs (or breaks on) them; run them in the separate
``maillard-jax`` env (environment-jax.yml). They pin what results/validation/jax_network_infra.md claims:

* a forced constant equals pH mode "prescribed"; a 'previous' pH + temperature course equals the engine's staircase
  run; a forced species equals the network with its row removed and the forced-exchange ledger closes the elements;
  a linear course gets a knot where it crosses a module's declared breakpoint; every refusal fires;
* the JAX fit driver reaches the engine's polished optimum on the lipid toy (live) and the trunk-lane toy (the
  artifact), and refuses an unconverged fit;
* ``densities.Predictor`` equals single solves and ``log_likelihoods`` equals ``gaussian_log_density``.
"""
from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"
REFERENCE = ROOT / "results" / "validation" / "jax_network_infra_reference.json"
ARTIFACT = ROOT / "results" / "validation" / "jax_network_infra.json"
TIGHT = {"rtol": 1e-10, "atol": 1e-16}


@pytest.fixture(scope="module")
def env():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_network_infra as bk
    from src.kinetic_core_jax import densities, fit, forcing, lipid, network, trunk_lane

    return {"bk": bk, "reference": json.loads(REFERENCE.read_text()), "trunk": bk.trunk_net(), "lipid": bk.lipid_net(),
            "fz": forcing, "fit": fit, "network": network, "trunk_lane": trunk_lane, "lipid_mod": lipid,
            "densities": densities}


def _forced(env, net, F, fp):
    return jax.jit(lambda th, p: env["fz"].integrate_forced(net, th, p, F, **TIGHT))(net.theta0(), fp.pot)


# --------------------------------------------------------------------------- forcing
def test_forcing_refusals(env):
    from src.kinetic_core_jax.lane import UnsupportedConfiguration

    fz, net, lnet = env["fz"], env["trunk"], env["lipid"]
    with pytest.raises(ValueError):
        fz.Forcing(inputs=("ph_label",), ph_scale="label_25C")            # no temperature course
    with pytest.raises(ValueError):
        fz.Forcing(inputs=("t_c", "ph_label"))                            # undeclared pH scale
    with pytest.raises(UnsupportedConfiguration):
        fz.Forcing(inputs=("t_c", "ph_label"), ph_scale="in_situ")
    with pytest.raises(ValueError):
        fz.Forcing(inputs=("t_c", "ionic_m"))
    with pytest.raises(ValueError):                                       # the lipid lane reads no pH
        fz.forced_pot(lnet, fz.Forcing(inputs=("t_c", "ph_label"), ph_scale="cooled"), np.zeros(lnet.n_species),
                      {"t_c": ([0, 10], [100, 100]), "ph_label": ([0, 10], [6, 6])})
    with pytest.raises(UnsupportedConfiguration):
        fz.forced_pot(lnet, fz.Forcing(inputs=("t_c",), species=("Gly",)), np.zeros(lnet.n_species),
                      {"t_c": ([0, 10], [100, 100]), "Gly": ([0, 10], [1, 1])})
    F = fz.Forcing(inputs=("t_c",))
    y0 = np.zeros(net.n_species)
    with pytest.raises(ValueError):                                       # knots must increase
        fz.forced_pot(net, F, y0, {"t_c": ([0, 10, 10], [100, 110, 120])}, scalars={"ph_label": 6.0, "aw": 0.9})
    with pytest.raises(ValueError):                                       # a scalar the network reads is missing
        fz.forced_pot(net, F, y0, {"t_c": ([0, 10], [100, 110])}, scalars={"ph_label": 6.0})
    with pytest.raises(ValueError):                                       # NaN is a scalar's None, never a course
        fz.forced_pot(net, fz.Forcing(inputs=("t_c", "aw")), y0, {"t_c": ([0, 10], [100, 110]),
                      "aw": ([0, 10], [0.9, float("nan")])}, scalars={"ph_label": 6.0})


def test_acrylamide_forced_course_keeps_the_cf106_refusal(env):
    from src.kinetic_core_jax import acrylamide, trunk
    from src.kinetic_core_jax.lane import UnsupportedConfiguration

    fz = env["fz"]
    net = env["network"].compose([trunk.build(), acrylamide.build()])
    with pytest.raises(UnsupportedConfiguration):
        fz.forced_pot(net, fz.Forcing(inputs=("t_c",)), np.zeros(net.n_species), {"t_c": ([0, 10], [60, 100])},
                      scalars={"ph_label": 6.0, "aw": 0.9})


def test_a_linear_course_gets_a_knot_at_each_declared_breakpoint(env):
    fz, net, lnet = env["fz"], env["trunk"], env["lipid"]
    assert fz.breakpoints(net, "ph_label") == (7.0,)
    assert fz.breakpoints(lnet, "t_c") == (100.0,)
    F = fz.Forcing(inputs=("t_c", "ph_label"), ph_scale="label_25C")
    fp = fz.forced_pot(net, F, np.zeros(net.n_species), {"t_c": ([0, 40], [120, 120]), "ph_label": ([0, 40], [7.5, 5.5])},
                       scalars={"aw": 0.9}, save_times=[20.0, 40.0])
    assert fp.crossings == [("ph_label", 7.0, 10.0)]
    assert np.allclose(fp.knots, [0.0, 10.0, 20.0, 40.0])
    assert np.allclose(fp.pot["start"][:, 1], [7.5, 7.0, 6.5]) and np.allclose(fp.pot["stop"][:, 1], [7.0, 6.5, 5.5])
    assert fp.rows == [1, 2]
    lfp = fz.forced_pot(lnet, fz.Forcing(inputs=("t_c",)), np.zeros(lnet.n_species), {"t_c": ([0, 20], [60, 160])})
    assert lfp.crossings == [("t_c", 100.0, 8.0)]
    # a staircase crosses only at its knots: no knot is added
    S = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("ph_label", "previous"),), ph_scale="cooled")
    sp = fz.forced_pot(net, S, np.zeros(net.n_species), {"t_c": ([0, 40], [120, 120]),
                                                          "ph_label": ([0, 20, 40], [7.5, 5.5, 5.5])},
                       scalars={"aw": 0.9})
    assert sp.crossings == [] and np.allclose(sp.pot["start"][:, 1], sp.pot["stop"][:, 1])
    assert fz.ph_mode(net, F) == "prescribed (trajectory)" and fz.ph_mode(net, None) == "prescribed"


def test_forced_constant_equals_the_prescribed_mode(env):
    fz, net, network = env["fz"], env["trunk"], env["network"]
    segs = [(6.0, 120.0), (14.0, 140.0)]
    pot = env["trunk_lane"].make_pot(net, {"Glc": 100.0, "Gly": 50.0, "MGO": 5.0}, segs, ph=6.4, water_activity=0.9)
    plain = np.asarray(jax.jit(lambda th, p: network.integrate(net, th, p, **TIGHT).ys)(net.theta0(), pot))
    F = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"),), ph_scale="label_25C")
    fp = fz.forced_pot(net, F, pot["y0"], {"t_c": fz.program_course(segs), "ph_label": ([0, 20], [6.4, 6.4])},
                       scalars={"aw": 0.9}, save_times=[6.0, 20.0])
    ys = np.asarray(_forced(env, net, F, fp).ys)[fp.rows]
    assert env["bk"].max_dex(ys, plain, 100.0) < 1e-10


def test_a_staircase_course_equals_the_engines_staircase_run(env):
    fz, net, st = env["fz"], env["trunk"], env["reference"]["forced_staircase"]
    segs = [tuple(s) for s in st["segments"]]
    knots = np.concatenate([[0.0], np.cumsum([d for d, _ in segs])])
    F = fz.Forcing(inputs=("t_c", "ph_label"), interpolation=(("t_c", "previous"), ("ph_label", "previous")),
                   ph_scale="label_25C")
    pot = env["trunk_lane"].make_pot(net, st["initial"], segs, ph=st["ph"][0], water_activity=st["water_activity"])
    fp = fz.forced_pot(net, F, pot["y0"], {"t_c": fz.program_course(segs), "ph_label": (knots, st["ph"] + st["ph"][-1:])},
                       scalars={"aw": st["water_activity"]}, save_times=knots[1:])
    ys = np.asarray(_forced(env, net, F, fp).ys)[fp.rows]
    eng = np.asarray(st["states"])[:, [st["engine_species"].index(s) for s in net.species]]
    assert env["bk"].max_dex(ys, eng, 100.0) < 1e-6


def test_a_forced_species_leaves_the_state_and_its_ledger_closes(env):
    fz, net, network = env["fz"], env["trunk"], env["network"]
    pot = env["trunk_lane"].make_pot(net, {"Glc": 100.0, "Gly": 50.0}, [(20.0, 130.0)], ph=6.8, water_activity=None)
    stoich = np.array(net.stoichiometry)
    stoich[net.index["Gly"], :] = 0.0
    net_wo = dataclasses.replace(net, stoichiometry=stoich)
    y_wo = np.asarray(jax.jit(lambda th, p: network.integrate(net_wo, th, p, **TIGHT).ys)(net.theta0(), pot))
    F = fz.Forcing(inputs=("t_c",), species=("Gly",))
    scal = {"ph_label": 6.8, "aw": float("nan")}
    fp = fz.forced_pot(net, F, pot["y0"], {"t_c": ([0, 20], [130, 130]), "Gly": ([0, 20], [50, 50])}, scalars=scal)
    run = _forced(env, net, F, fp)
    assert env["bk"].max_dex(np.asarray(run.ys), y_wo, 100.0) < 1e-9
    fp2 = fz.forced_pot(net, F, pot["y0"], {"t_c": ([0, 20], [130, 130]), "Gly": ([0, 5, 20], [50, 45, 30])},
                        scalars=scal, save_times=[5.0, 20.0])
    run2 = _forced(env, net, F, fp2)
    assert np.allclose(np.asarray(run2.ys)[:, net.index["Gly"]], [45.0, 30.0])
    assert float(np.max(np.abs(fz.closure(net, F, fp2.pot["y0"], run2)))) < 1e-10
    assert float(np.asarray(run2.exchange)[-1, 0]) < 0.0          # the network draws glycine from the clamp


# --------------------------------------------------------------------------- compose hooks
def test_declared_competition_and_the_titratable_table(env):
    from src.kinetic_core_jax import trunk, trunk_lane
    from src.kinetic_core_jax.lane import UnsupportedConfiguration

    network = env["network"]
    net = env["trunk"]
    undeclared = set(net.undeclared_competition)
    assert ("Glc", "trunk_lane") in undeclared and ("MGO", "trunk_lane") in undeclared
    assert all(owner == "trunk" for s, (owner, *_) in net.competition.items() if s in ("Glc", "MGO", "Gly"))
    mods = [trunk.build(lane="trunk"), trunk_lane.build(), trunk_lane.build_fed_3deoxy()]
    with pytest.raises(ValueError, match="undeclared shared consumption"):
        network.compose(mods, competition="refuse")
    declared = [mods[0]] + [dataclasses.replace(m, competes=tuple(sorted({s for s, c in undeclared if c == m.name})))
                            for m in mods[1:]]
    ok = network.compose(declared, competition="refuse")
    assert ok.undeclared_competition == () and tuple(ok.reactions) == tuple(net.reactions)
    with pytest.raises(ValueError, match="no module owns"):
        network.compose([mods[0], dataclasses.replace(mods[1], competes=("NOT_A_SPECIES",)), mods[2]])
    t = network.Titratable("Gly", (2.34, 9.6), (1, 0, -1), "test")
    with_t = network.compose([dataclasses.replace(mods[0], titratables=(t,)), mods[1], mods[2]])
    assert with_t.titratables == (t,)
    with pytest.raises(ValueError, match="declared twice"):
        network.compose([dataclasses.replace(mods[0], titratables=(t,)), dataclasses.replace(mods[1], titratables=(t,)),
                         mods[2]])
    with pytest.raises(ValueError, match="charge states"):
        network.compose([dataclasses.replace(mods[0], titratables=(network.Titratable("Gly", (2.3,), (1,)),)),
                         mods[1], mods[2]])
    with pytest.raises(UnsupportedConfiguration):
        network.charge_closure(with_t)


# --------------------------------------------------------------------------- the fit driver
def test_the_fit_driver_reaches_the_engines_polished_optimum_on_the_lipid_toy(env):
    fit, toy = env["fit"], env["reference"]["lipid_toy"]
    problem = env["bk"].toy_problem(toy, "lipid")
    result = fit.fit(problem, toy["start"], budget=200, require_converged="raise")
    eng = np.asarray(list(toy["polished"]["theta"].values()))
    sigma = fit.laplace(problem, result)["sigma"]
    assert np.max(np.abs(result.theta - eng) / sigma) < 1e-3
    assert abs(result.cost - toy["polished"]["cost"]) / toy["polished"]["cost"] < 1e-6


def test_the_fit_artifact_records_equivalence_on_both_toys():
    if not ARTIFACT.exists():
        pytest.skip("results/validation/jax_network_infra.json not written yet")
    fits = json.loads(ARTIFACT.read_text())["fit"]
    for lane in ("trunk", "lipid"):
        for batch in ("vmap", "loop"):
            row = fits[lane][batch]
            assert row["status"] > 0 and not row["failed_pots"]
            assert row["max_delta_in_sigma"] < 1e-2, (lane, batch, row["max_delta_in_sigma"])
            assert abs(row["cost_relative_difference"]) < 1e-5, (lane, batch, row["cost_relative_difference"])
        assert fits[lane]["fd_check_at_optimum"]["max_relative_column_error"] < 1e-4


def test_an_unconverged_fit_is_refused(env):
    from src.kinetic_core.fitting import FitNotConverged

    fit, toy = env["fit"], env["reference"]["lipid_toy"]
    problem = env["bk"].toy_problem(toy, "lipid")
    with pytest.raises(FitNotConverged):
        fit.fit(problem, toy["start"], budget=1, require_converged="raise")
    flagged = fit.fit(problem, toy["start"], budget=1, require_converged="flag")
    assert flagged.converged is False


def test_problem_kinds_and_bounds(env):
    fit, net = env["fit"], env["trunk"]
    assert fit.kind_of("log10_k:k_schiff") == "log10k" and fit.kind_of("ea:k_schiff") == "Ea"
    assert fit.kind_of("k_ref:k_go_met") == "linear" and fit.kind_of("ph_exponent") == "ph_slope"
    with pytest.raises(AssertionError, match="bounded below by 0"):
        fit.JaxProblem(net=net, keys=(("trunk_lane", "k_ref:k_go_met"),), kinds=("linear",), lower=np.array([-1.0]),
                       upper=np.array([1.0]), pots=(), observations=())
    with pytest.raises(AssertionError, match="unknown coordinate kinds"):
        fit.JaxProblem(net=net, keys=(("trunk", "ph_exponent"),), kinds=("slope",), lower=np.array([0.0]),
                       upper=np.array([1.0]), pots=(), observations=())


# --------------------------------------------------------------------------- the density bridge
def test_predictor_equals_single_solves_and_the_likelihoods_equal_gaussian_log_density(env):
    from src.kinetic_core import predictive_densities as pdm

    densities, network, fit = env["densities"], env["network"], env["fit"]
    toy = env["reference"]["lipid_toy"]
    problem = env["bk"].toy_problem(toy, "lipid")
    net = problem.net
    obs = [o for o in problem.observations if o.pot == 1]
    rng = np.random.default_rng(3)
    draws = np.asarray(toy["start"]) + rng.normal(size=(5, 4)) * np.asarray([0.1, 3.0, 0.02, 0.05])
    draws[:, 2:] = np.clip(draws[:, 2:], 0.0, 0.9)
    pred = densities.Predictor(net, [problem.pots[1]], obs, problem.keys, chunk=4)
    values, failed = pred(draws)
    assert not failed.any() and values.shape == (5, len(obs))
    table = {k: f for m in net.modules for k, f in m.observables.items()}
    base = np.asarray(network.theta_vector(net, net.theta0()))
    where = {c: i for i, c in enumerate(net.coordinates)}
    for x, v in zip(draws, values):
        vec = base.copy()
        for k, xi in zip(problem.keys, x):
            vec[where[k]] = xi
        ys = np.asarray(network.integrate(net, network.theta_from_vector(net, vec), problem.pots[1], rtol=1e-8,
                                          atol=1e-14).ys)
        ref = [np.log10(float(table[o.observable](ys[o.row], net.index))) for o in obs]
        assert np.max(np.abs(np.asarray(ref) - v)) < 1e-10
    y = np.asarray([o.value for o in obs])
    cov = pdm.series_covariance([o.sigma for o in obs], "absolute", omega_dex=0.1)
    ll = densities.log_likelihoods(values, y, cov)
    assert np.allclose(ll, [pdm.gaussian_log_density(y, v, cov) for v in values], rtol=0, atol=1e-9)
    bad = values.copy()
    bad[0, 0] = np.nan
    with pytest.raises(pdm.MeasurementModelRefusal):
        densities.log_likelihoods(bad, y, cov)
    assert densities.log_likelihoods(bad, y, cov, on_failed="impossible")[0] == -np.inf
    g = densities.log_integrand(pred, y, cov, lambda th: np.zeros(len(th)))
    assert np.allclose(g(draws), ll)
