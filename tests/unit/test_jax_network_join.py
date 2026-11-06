"""The network's own pH solve and the composed sulfur lane (A2 contract v2.3; ``src/kinetic_core_jax/ph_coupling.py``,
``sulfur_lane.build_composed``): structure only, no integration (the integrations are
tests/scientific/test_jax_network_join.py and results/validation/jax_network_join.md).

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
                         ({**trunk.centre_ledger, "r_glc_fru": {"carboxyl": 1, "basis": "stale"}}, "stale")):
        net = network.compose([dataclasses.replace(trunk, centre_ledger=ledger)] + mods[1:],
                              layout=comp.network.species, competition="refuse")
        with pytest.raises(ValueError, match=what):
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
    # a pH-free pot is refused by the network's pH solve, not run without its pH
    with pytest.raises(ValueError, match="inline-pH pot"):
        network.integrate(one.network, one.theta0(), {"y0": np.zeros(one.network.n_species),
                                                      "segments": np.array([[1.0, 100.0]])})


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
