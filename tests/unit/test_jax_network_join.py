"""The network's own pH solve and the composed sulfur lane (A2 contract v2.3; ``src/kinetic_core_jax/ph_coupling.py``,
``sulfur_lane.build_composed``): structure only, no integration (the integrations are
tests/scientific/test_jax_network_join_parity.py and results/validation/jax_network_join.md).

Skipped unless jax and diffrax import (the core env); run in the ``maillard-jax`` env.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")


@pytest.fixture(scope="module")
def lanes():
    from src.kinetic_core_jax import sulfur_lane as sl

    return sl.build(), sl.build_composed()


def test_the_composed_lane_is_the_one_module_lane(lanes):
    one, comp = lanes
    a, b = one.network, comp.network
    assert [m.name for m in comp.modules] == ["sulfur_trunk", "sulfur_fed_3deoxy", "sulfur"]
    assert tuple(a.reactions) == tuple(b.reactions) and tuple(a.species) == tuple(b.species)
    assert np.array_equal(a.stoichiometry, b.stoichiometry)
    assert b.undeclared_competition == () and b.competition       # declared, and refused if not
    assert sorted(k for m in comp.modules for k in m.params) == sorted(one.module.params)
    assert b.ph_mode == "inline" and b.ph_model == a.ph_model


def test_the_merged_table_is_the_engines_inventory_and_closes(lanes):
    """The ONE closure check over the merged table equals the engine's two tables: the per-species centres are
    ``TITRATABLE_CENTRES`` and each step's movement is ``ph_state.centre_delta``; the lane's ``CENTRE_LEDGER``
    declares exactly those steps."""
    from src.kinetic_core import ph_state, sulfur
    from src.kinetic_core_jax import network

    for lane in lanes:
        net = lane.network
        report = network.charge_closure(net, ph=6.0)
        assert report["closes"], report["problems"]
        engine = {r.key: {k: v for k, v in ph_state.centre_delta(r.reactants, r.products).items() if v}
                  for r in sulfur.FULL_REACTIONS}
        engine = {k: v for k, v in engine.items() if v}
        assert set(report["moved"]) == set(engine) == set(sulfur.CENTRE_LEDGER)
        for r, d in engine.items():
            assert all(abs(report["moved"][r].get(k, 0.0) - v) < 1e-9 for k, v in d.items()), r
        sc = net.charge.structural_weights().T @ net.charge.centres.T
        for s, c in ph_state.TITRATABLE_CENTRES.items():
            if s in net.index:
                for k, v in c.items():
                    assert abs(sc[net.index[s], net.charge.kinds.index(k)] - v) < 1e-12, (s, k)
        assert sorted(net.charge.buffer_rows) and net.charge.yields and net.charge.pka_slots


def test_a_ledger_that_drifts_from_the_table_is_refused(lanes):
    from src.kinetic_core_jax import network

    _, comp = lanes
    mods = list(comp.modules)
    trunk = mods[0]
    key = sorted(trunk.centre_ledger)[0]
    for ledger, what in (({k: v for k, v in trunk.centre_ledger.items() if k != key}, "not declared"),
                         ({**trunk.centre_ledger, key: {**trunk.centre_ledger[key], "carboxyl": 7}}, "declares"),
                         ({**trunk.centre_ledger, "r_glc_fru": {"carboxyl": 1, "basis": "stale"}}, "stale"),
                         ({**trunk.centre_ledger, key: {**trunk.centre_ledger[key], "sulfhydryl": 1}},
                          "no pool of the table carries")):
        broken = [dataclasses.replace(trunk, centre_ledger=ledger)] + mods[1:]
        with pytest.raises(ValueError, match="charge closure fails"):          # compose refuses it ...
            network.compose(broken, layout=comp.network.species, competition="refuse")
        net = network.compose(broken, layout=comp.network.species, competition="refuse", check_closure=False)
        with pytest.raises(ValueError, match=what):                            # ... and says why
            network.charge_closure(net, strict=True)


def test_refusals_of_the_v2_3_contract(lanes):
    from src.kinetic_core_jax import network, trunk
    from src.kinetic_core_jax.lane import UnsupportedConfiguration

    one, comp = lanes
    mods = list(comp.modules)
    sulfur_mod = mods[2]
    # a second, different PhModel
    other = dataclasses.replace(sulfur_mod.ph_model, tolerance=2e-3)
    with pytest.raises(ValueError, match="different PhModels"):
        network.compose([dataclasses.replace(mods[0], ph_model=other)] + mods[1:])
    # a yield coordinate the module does not own
    pools = list(sulfur_mod.titratables)
    bad = [t._replace(yield_coordinate="not_a_coordinate") if t.yield_coordinate else t for t in pools]
    with pytest.raises(ValueError, match="not sulfur's"):
        network.compose(mods[:2] + [dataclasses.replace(sulfur_mod, titratables=tuple(bad))])
    # a pool on a species the module neither owns nor requires
    stray = [dataclasses.replace(mods[0], titratables=(network.Titratable("x", (4.0,), (0, -1),
                                                                          weights={"Cys": 1.0}),))]
    with pytest.raises(ValueError, match="neither owns nor requires"):
        network.compose(stray + mods[1:])
    # the declared-pH trunk beside the inline-pH sulfur lane: two pH conventions
    with pytest.raises(UnsupportedConfiguration, match="two things"):
        network.compose([trunk.build(lane="trunk"), dataclasses.replace(one.module, species=(), elements={},
                                                                       reactions=(), stoichiometry={},
                                                                       reactants={}, titratables=())])
    # an a_w reader, and a pH-free module's pools, beside the inline pH (not carried / A2b)
    with pytest.raises(UnsupportedConfiguration, match="a_w"):
        network.compose([dataclasses.replace(mods[0], pot_inputs=("aw",))] + mods[1:])
    lipid_like = dataclasses.replace(mods[0], needs_ph=False, titratables=(network.Titratable("Gly", (2.3, 9.6),
                                                                                              (1, 0, -1)),))
    with pytest.raises(UnsupportedConfiguration, match="solve no pH"):
        network.compose([lipid_like] + mods[1:], check_closure=False)
    # a pH-free pot is refused by the network's pH solve, not run without its pH
    with pytest.raises(ValueError, match="inline-pH pot"):
        network.integrate(one.network, one.theta0(), {"y0": np.zeros(one.network.n_species),
                                                      "segments": np.array([[1.0, 100.0]])})


def test_rows_outside_a_pot_are_refused_by_the_fit_driver_and_the_density_bridge(lanes):
    from src.kinetic_core import engine
    from src.kinetic_core_jax import densities, fit, network
    from src.kinetic_core_jax import sulfur_lane as sl

    one, _ = lanes
    net = one.network
    pot = sl.make_pot(one, {"Cys": 10.0, "PENT": 10.0},
                      engine.ProcessSpec(thermal=engine.ThermalProgram.isothermal(120.0, 10.0), ph=6.0))
    n = network.n_rows(net, pot)
    key = ("sulfur", "acid_yield")
    for row, ok in ((n - 1, True), (n, False), (-1, False)):
        obs = (fit.Observation(pot=0, row=row, observable="mmol_per_L:MFT", value=-3.0, sigma=0.1),)
        build = [lambda: fit.Compiled(fit.JaxProblem(net=net, keys=(key,), kinds=("acid_yield",),
                                                     lower=np.array([0.0]), upper=np.array([1.0]), pots=(pot,),
                                                     observations=obs, batch="loop")),
                 lambda: densities.Predictor(net, [pot], obs, [key])]
        for make in build:
            if ok:
                make()
            else:
                with pytest.raises(ValueError, match="outside their pot"):
                    make()


def test_sweep_changes_names_the_pots_whose_count_moved():
    from types import SimpleNamespace

    from src.kinetic_core_jax import fit

    log = SimpleNamespace(sweep_log={0: [(7,), (7,), (8,), (7,)], 1: [(6,), (6,)], 2: []})
    assert fit.Compiled.sweep_changes(log) == {0: [(7,), (8,)]}


def test_the_debug_root_check_fires_only_on_an_unconverged_root(lanes, monkeypatch):
    import jax.numpy as jnp

    from src.kinetic_core_jax import ph_coupling

    one, _ = lanes
    net = one.network
    theta = one.theta0()
    y = np.zeros(net.n_species)
    y[net.index["Cys"]], y[net.index["CBX"]] = 20.0, 5.0
    buffer = jnp.zeros(net.charge.n_pools)
    monkeypatch.setattr(ph_coupling, "DEBUG_ROOT", True)
    solver = ph_coupling.solver_of(net)
    ph, ionic, cool = solver.state(jnp.asarray(y), 393.15, 3.0, theta, buffer)
    assert np.isfinite(float(ph)) and 0.0 < float(ionic)
    monkeypatch.setattr(ph_coupling, "NEWTON_STEPS", 0)               # the bisection alone: 1e-3 pH, not a root
    with pytest.raises(Exception, match="did not converge"):
        ph_coupling.solver_of(net).state(jnp.asarray(y), 393.15, 3.0, theta, buffer)


def test_a_profile_pot_merges_near_coincident_knots_for_the_step_controller(lanes):
    from src.kinetic_core_jax import ph_coupling

    one, _ = lanes
    net = one.network
    pot = ph_coupling.pot_from_program(net, np.zeros(net.n_species), 6.0, np.zeros(net.charge.n_pools),
                                       profile=((0.0, 25.0), (5.0, 140.0), (5.0 + 1e-14, 140.0), (20.0, 140.0)))
    knots, jumps = np.asarray(pot["knot_s"])[0], np.asarray(pot["knot_jump_s"])[0]
    assert np.all(np.diff(knots) > 0.0) and np.min(np.diff(jumps)[np.diff(jumps) > 0]) >= ph_coupling.JUMP_MERGE


def test_pot_rows_and_arrays(lanes):
    from src.kinetic_core import engine
    from src.kinetic_core_jax import network
    from src.kinetic_core_jax import sulfur_lane as sl

    one, _ = lanes
    process = engine.ProcessSpec(thermal=engine.ThermalProgram(((5.0, 100.0), (15.0, 140.0))), ph=6.0)
    pot = sl.make_pot(one, {"Cys": 10.0, "PENT": 10.0}, process)
    rows = network.segment_end_rows(one.network, pot)
    n = np.shape(pot["node_s"])[1]
    assert rows == [n - 1, 2 * n - 1]
    assert set(network.pot_arrays(one.network, pot)) == set(sl.ph_coupling.POT_KEYS)
    assert network.pot_mode(one.network, pot) == "dynamic"
