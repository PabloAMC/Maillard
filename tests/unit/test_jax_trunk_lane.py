"""src/kinetic_core_jax: the trunk LANE (``network.TRUNK_REACTIONS``, 55 steps) as three A2 contract v2.1 modules --
``trunk.build(lane="trunk")``, ``trunk_lane.build`` (the six trunk-only blocks) and ``trunk_lane.build_fed_3deoxy``.

Skipped unless jax and diffrax import, so the core env's CI never runs (or breaks on) them; run them in the separate
``maillard-jax`` env (environment-jax.yml). They pin what results/validation/jax_trunk_lane_module.md claims:

* every rate constant equals the engine's (``network.rate_constants_at`` after ``trunk_conditions.apply``) to 1e-12
  relative, on the reference probe and live, with the declared pH read in pH mode "prescribed";
* the panel pot the engine routes to the trunk lane, at its point and both band corners, four synthetic isothermal
  pots (one per block it leaves idle) and a staircase profile at the shipped point and with B22 / B24 on match the
  engine's states after every segment to 1e-6 dex at tight tolerance; the elements close;
* the forward-mode Jacobian equals central differences; the panel predictions at rtol 1e-6 equal engine.predict's to
  1e-4 dex;
* every switch the lane reads is ported or refused.
"""
from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")

import jax.numpy as jnp  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "generators"
REFERENCE = ROOT / "results" / "validation" / "jax_trunk_lane_module_reference.json"


@pytest.fixture(scope="module")
def setup():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_trunk_lane_module as bk
    import build_jax_trunk_lane_reference as refgen
    from src.kinetic_core_jax import network, trunk, trunk_lane

    reference, net = bk.setup()
    return {"bk": bk, "refgen": refgen, "reference": reference, "net": net, "trunk": trunk,
            "trunk_lane": trunk_lane, "network": network}


def _tight(net, network):
    return jax.jit(lambda th, pot: network.integrate(net, th, pot, rtol=1e-10, atol=1e-16).ys)


def test_the_composed_lane_is_the_engines_trunk_lane(setup):
    from src.kinetic_core import network as engine_network
    from src.kinetic_core.species import SPECIES_KEYS

    net = setup["net"]
    assert tuple(net.reactions) == engine_network.TRUNK_REACTION_KEYS and len(net.reactions) == 55
    assert sorted(net.species) == sorted(SPECIES_KEYS)
    assert [m.name for m in net.modules] == ["trunk", "trunk_lane", "fed_3deoxy"]
    assert [len(m.reactions) for m in net.modules] == [26, 26, 3]
    assert net.ph_mode == "prescribed"


def test_rate_constants_equal_the_engine_on_the_reference_probe(setup):
    out = setup["bk"].stage_rates(setup["reference"], setup["net"])
    assert out["n_rate_values"] > 100_000
    assert out["worst_relative"] < 1e-12, out


def test_rate_constants_equal_the_engine_live(setup):
    """Against the engine's functions in this process, at points and a parameter set the reference does not hold,
    with B22 / B24 on through the engine's own override blocks."""
    from src.kinetic_core import config, engine, network as engine_network, trunk_conditions

    bk, refgen, network = setup["bk"], setup["refgen"], setup["network"]
    spec = {"methionine_b22": True, "proline_b24": True}
    params = engine.core_parameters(engine.TRUNK, frozen=refgen.frozen_for(spec))
    params["k_glyc"] = dataclasses.replace(params["k_glyc"], k_ref=float(params["k_glyc"].k_ref) * 3.0, ea_kj_mol=80.0)
    params["k_da_sink"] = dataclasses.replace(params["k_da_sink"], k_ref=2e-4)
    cond = refgen.full_conditions({"aw_scale": 0.9, "ph_exponent": 0.45, "pyrazine_ph_slope_above_knot": 0.12,
                                   "pyrazine_ph_slope_below_knot": 0.71})
    net = bk.network_for(spec)
    theta = bk.theta_of({"parameters": refgen.parameter_table(params), "conditions": cond, "config_spec": spec}, net)
    with config.applied(config.Config.from_globals()):
        for t_c, ph, aw in ((151.3, 5.3, 0.77), (97.0, 8.6, 0.30), (205.0, None, 0.95), (133.0, 6.8, None),
                            (118.0, 7.0, 0.6), (140.0, 6.95, 0.99)):
            applied, _ = trunk_conditions.apply(params, SimpleNamespace(ph=ph, water_activity=aw, atmosphere=None),
                                                **refgen.apply_kwargs(cond))
            expect = engine_network.rate_constants_at(applied, t_c + 273.15)
            ctx = network.RateCtx(t_c=jnp.asarray(t_c), t_k=jnp.asarray(t_c + 273.15),
                                  ph_label=jnp.asarray(setup["trunk"].nan_or(ph)),
                                  aw=jnp.asarray(setup["trunk"].nan_or(aw)))
            k = np.asarray(net.rate_constants(theta, ctx))
            e = np.array([expect[r] for r in net.reactions])
            np.testing.assert_allclose(k, e, rtol=1e-12, atol=0.0, err_msg=str((t_c, ph, aw)))


@pytest.mark.parametrize("group", ["panel", "synthetic isothermal"])
def test_trajectories_match_the_engine_on_the_isothermal_pots(setup, group):
    bk, reference, network = setup["bk"], setup["reference"], setup["network"]
    recs = [r for r in reference["pots"] if bk.group_of(r) == group]
    assert len(recs) == (3 if group == "panel" else 4)
    cols = bk.engine_columns_of_species(reference, setup["net"])
    for rec in recs:
        net = bk.network_for(rec["config_spec"])
        pot = bk.pot_of(rec, net)
        ys = np.asarray(_tight(net, network)(bk.theta_of(rec, net), pot))
        c = bk.compare(ys, rec, pot["y0"], net, cols)
        assert c["max_dex"] < 1e-6, (rec["id"], rec["variant"], c)
        assert c["max_abs_below_floor_share_of_y0"] < 1e-9, (rec["id"], rec["variant"], c)
        e0 = net.element_matrix @ pot["y0"]
        np.testing.assert_allclose(net.element_matrix @ ys[-1], e0, rtol=1e-9, atol=1e-12)


def test_a_staircase_profile_with_every_block_matches_the_engine(setup):
    bk, reference, network = setup["bk"], setup["reference"], setup["network"]
    recs = [r for r in reference["pots"] if bk.group_of(r) == "profile"]
    assert {r["variant"] for r in recs} == {"point", "switches_on"}
    cols = bk.engine_columns_of_species(reference, setup["net"])
    for rec in recs:
        assert len(rec["refined_segments"]) > 50
        net = bk.network_for(rec["config_spec"])
        pot = bk.pot_of(rec, net)
        theta = bk.theta_of(rec, net)
        ys = np.asarray(_tight(net, network)(theta, pot))
        c = bk.compare(ys, rec, pot["y0"], net, cols)
        assert c["max_dex"] < 1e-6, (rec["variant"], c)
        # every step with a non-zero constant carries flux mid-programme
        k = np.asarray(net.rate_constants(theta, network.RateCtx(t_c=160.0, t_k=433.15, ph_label=5.5, aw=0.8)))
        rates = k * np.asarray(network.mass_action(jnp.asarray(ys[len(ys) // 2]), net))
        idle = [r for r, kk, v in zip(net.reactions, k, rates) if kk > 0.0 and not v > 0.0]
        assert not idle, idle
    on = next(r for r in recs if r["variant"] == "switches_on")
    final = dict(zip(reference["engine_species"], on["states"][-1]))
    assert final["MTAL"] > 0.0 and final["AP"] > 0.0      # the B22 / B24 steps run when switched on


def test_exact_jacobian_equals_central_differences(setup):
    bk, reference, network = setup["bk"], setup["reference"], setup["network"]
    pot_id = "mp_holdout_glucose_only_autoclave_121C_Steinhagen2021"
    fd = reference["jacobian_central_fd"][pot_id]
    rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
    net = bk.network_for(fd["config_spec"])
    obs_idx = jnp.asarray([net.index[s] for s in fd["species"]])
    rows = jnp.asarray(fd["rows"])

    def observe(ys):
        return jnp.log10(ys[rows][:, obs_idx].T.ravel())

    pot = bk.pot_of(rec, net)
    theta = network.theta_vector(net, bk.theta_of(fd, net))
    _, jac = network.jacobian(net, observe, rtol=1e-10, atol=1e-16)(theta, pot)
    jac = np.asarray(jac)
    keys = [k for _, k in net.coordinates]
    cols = [keys.index(c) for c in fd["columns"]]
    errors = bk.column_errors(jac[:, cols], np.asarray(fd["jacobian"]))
    assert errors["n_relative_columns"] >= 20
    assert errors["worst"] < 2e-4, errors
    assert errors["worst_absolute"] < 1e-6, errors
    f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), pot,
                                                    rtol=1e-12, atol=1e-20).ys))
    h = fd["step"]
    for name in ("log10_k:k_glc_g", "ea:k_g_go", "ph_exponent", "ph_exponent_exit:k_tdg_fa", "log10_k:k_tdg_ddg",
                 "ea:k_ddg_dgal"):
        c = keys.index(name)
        e = jnp.zeros(len(theta)).at[c].set(h)
        col = (np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h)
        np.testing.assert_allclose(jac[:, c], col, rtol=1e-5, atol=1e-8, err_msg=name)


def test_the_jacobian_in_a_linear_k_ref_at_zero_equals_a_one_sided_difference(setup):
    """An unshipped step (k_go_met = 0, B22) has a LINEAR coordinate; its derivative at zero is the physical one-sided
    one. The difference's step is scaled to the column: a fixed 1e-5 L/(mmol min) is ~30x glycine's own constant."""
    trunk_lane, network, net = setup["trunk_lane"], setup["network"], setup["net"]
    pot = trunk_lane.make_pot(net, {"GO": 5.0, "MGO": 5.0, "Gly": 5.0, "MET": 5.0}, [(5.0, 120.0), (5.0, 120.0)],
                              ph=7.5, water_activity=None)
    idx = jnp.asarray([net.index[s] for s in ("AKG", "PZ", "MPZ")])

    def observe(ys):
        return jnp.log10(ys[:, idx].ravel())

    theta = network.theta_vector(net, net.theta0())
    _, jac = network.jacobian(net, observe, rtol=1e-10, atol=1e-16)(theta, pot)
    jac = np.asarray(jac)
    keys = [k for _, k in net.coordinates]
    c = keys.index("k_ref:k_go_met")
    assert float(theta[c]) == 0.0 and np.max(np.abs(jac[:, c])) > 1.0
    f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), pot,
                                                    rtol=1e-12, atol=1e-20).ys))
    h = 1e-5 / float(np.max(np.abs(jac[:, c])))
    e = jnp.zeros(len(theta)).at[c].set(h)
    col = (-3.0 * np.asarray(f(theta)) + 4.0 * np.asarray(f(theta + e)) - np.asarray(f(theta + 2 * e))) / (2 * h)
    np.testing.assert_allclose(jac[:, c], col, rtol=1e-6, atol=1e-6 * np.max(np.abs(col)))


def test_panel_predictions_at_fit_tolerance_equal_the_engine(setup):
    out = setup["bk"].stage_residuals(setup["reference"], setup["net"])
    for row in out["rows"]:
        assert row["n_rows"] == 6
        assert row["max_abs_diff_dex"] < 1e-4, row


def test_build_phase_switches_are_ported_or_refused(setup):
    from src.kinetic_core import config, engine
    from src.kinetic_core.parameters_methionine import FROZEN_B22

    trunk_lane, network = setup["trunk_lane"], setup["network"]
    shipped = config.Config.shipped()
    off = trunk_lane.build(shipped)
    assert "k_ref:k_go_met" in off.params and off.params["k_ref:k_go_met"] == 0.0
    on = trunk_lane.build(shipped.replace({"methionine.shipped_b22": True}))
    expect = engine.core_parameters(engine.TRUNK, frozen={"methionine": dict(FROZEN_B22)})["k_go_met"]
    assert on.params["log10_k:k_go_met"] == pytest.approx(np.log10(expect.k_ref), abs=1e-14)
    assert "log10_k:k_pyrl_ap" in trunk_lane.build(shipped.replace({"proline.shipped_b24": True})).params
    for key in ("methionine.shipped_b22b", "proline.shipped_b24b"):
        with pytest.raises(network.UnsupportedConfiguration, match="empty"):
            trunk_lane.build(shipped.replace({key: True}))
    # dicarbonyl.shipped_b39 moves no number
    a, b = (trunk_lane.build_fed_3deoxy(shipped.replace({"dicarbonyl.shipped_b39": v})) for v in (True, False))
    assert a.params == b.params


def test_rate_phase_switches_are_ported_or_refused(setup):
    from src.kinetic_core import config

    trunk_lane, trunk, network = setup["trunk_lane"], setup["trunk"], setup["network"]
    shipped = config.Config.shipped()
    m = dict(shipped["dicarbonyl.aqueous_glyoxal_parameters"])
    m["k_ama_g"] = dataclasses.replace(m["k_ama_g"], k_ref=float(m["k_ama_g"].k_ref) * 2.0)
    moved = trunk_lane.build(shipped.replace({"dicarbonyl.aqueous_glyoxal_parameters": m}))
    base = trunk_lane.build(shipped)
    assert moved.params["log10_k:k_ama_g"] - base.params["log10_k:k_ama_g"] == pytest.approx(np.log10(2.0), abs=1e-14)
    exit_map = dict(shipped["trunk.three_deoxy_exit_ph"])
    for step in ("k_go_sink", "k_ddg_tdg"):
        named = shipped.replace({"trunk.three_deoxy_exit_ph": dict(exit_map, **{step: (0.5, (0.4, 0.6), "x")})})
        with pytest.raises(network.UnsupportedConfiguration, match=step):
            trunk_lane.compose_lane(named)
    # the exit term off is ported (a structurally different trunk module)
    off = trunk_lane.compose_lane(shipped.replace({"trunk.three_deoxy_exit_ph_term": False}))
    assert not any(k.startswith("ph_exponent_exit:") for _, k in off.coordinates)
    assert trunk.build(shipped, lane="trunk").pot_inputs == ("ph_label", "aw")


def test_the_pot_folds_amines_and_refuses_what_the_engine_refuses(setup):
    trunk_lane, network, net = setup["trunk_lane"], setup["network"], setup["net"]
    pot = trunk_lane.make_pot(net, {"Glc": 10.0, "Gly": 1.0, "MET": 2.0, "PRO": 3.0}, [(10.0, 120.0)], ph=6.0,
                              water_activity=None)
    assert pot["y0"][net.index["Gly"]] == 6.0 and pot["y0"][net.index["MET"]] == 2.0
    assert pot["y0"][net.index["PRO"]] == 3.0
    lysp = trunk_lane.make_pot(net, {"Glc": 10.0}, [(10.0, 120.0)], ph=6.0, water_activity=None, lysp=4.0)
    assert lysp["y0"][net.index["LYSP"]] == 4.0
    for atmosphere in ("argon", "air_cu"):
        with pytest.raises(network.UnsupportedConfiguration, match="atmosphere"):
            trunk_lane.make_pot(net, {"Glc": 10.0}, [(10.0, 120.0)], ph=6.0, water_activity=None,
                                atmosphere=atmosphere)
    with pytest.raises(network.UnsupportedConfiguration, match="Asn"):
        trunk_lane.make_pot(net, {"Glc": 10.0, "Asn": 1.0}, [(10.0, 120.0)], ph=6.0, water_activity=None)
    with pytest.raises(ValueError, match="lacks"):
        network.integrate(net, net.theta0(), {"y0": pot["y0"], "segments": pot["segments"]})


def test_compose_refuses_the_blocks_without_the_trunk(setup):
    network, net = setup["network"], setup["net"]
    trunk_m, lane_m, fed_m = net.modules
    with pytest.raises(ValueError, match="requires"):
        network.compose([lane_m, fed_m])
    bad = dict(lane_m.stoichiometry)
    bad["r_go_ak"] = {"GO": -1.0, "Gly": -1.0, "AKG": 1.0, "FRAG_C": 1.0}
    with pytest.raises(ValueError, match="does not balance"):
        network.compose([trunk_m, dataclasses.replace(lane_m, stoichiometry=bad), fed_m])


def test_the_panel_pot_has_one_programme_under_both_thermal_histories(setup):
    histories = setup["reference"]["panel_thermal_history"]
    assert histories and all(h["isothermal"] == h["as_recorded"] for h in histories.values())
    assert json.dumps(setup["reference"]["blocks"]["FED_3DEOXY_REACTIONS"]) == json.dumps(
        list(setup["net"].modules[2].reactions))
