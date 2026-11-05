"""src/kinetic_core_jax: the trunk (``network.REACTIONS``) and the acrylamide block as A2 contract v2 modules, composed
into the engine's acrylamide lane (trunk.py, acrylamide.py).

Skipped unless jax and diffrax import, so the core env's CI never runs (or breaks on) them; run them in the separate
``maillard-jax`` env (environment-jax.yml). They pin what results/validation/jax_trunk_acrylamide_module.md claims:

* every rate constant equals the engine's (``acrylamide_rate_constants_at`` after ``acrylamide_conditions.apply``;
  ``network.rate_constants_at`` after ``trunk_conditions.apply`` for the trunk-lane terms) to 1e-12 relative, on the
  reference probe and live;
* every acrylamide panel pot, at its point and both band corners, and a staircase profile, match the engine's states
  after every segment to 1e-6 dex at tight tolerance; the elements close;
* the forward-mode Jacobian equals central differences; the panel predictions at rtol 1e-6 equal engine.predict's to
  1e-4 dex;
* every switch the lane reads is ported or refused: CF-106's low-temperature refusal, an uncarried species, a 3-DG
  exit term outside the module, a non-air atmosphere on the trunk lane; compose refuses acrylamide without the trunk
  and acrylamide with another owner of cysteine.
"""
from __future__ import annotations

import dataclasses
import json
import math
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
REFERENCE = ROOT / "results" / "validation" / "jax_trunk_acrylamide_module_reference.json"


@pytest.fixture(scope="module")
def setup():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_trunk_acrylamide_module as bk
    from src.kinetic_core_jax import acrylamide, network, trunk

    reference, net = bk.setup()
    return {"bk": bk, "reference": reference, "net": net, "trunk": trunk, "acrylamide": acrylamide,
            "network": network}


def test_rate_constants_equal_the_engine_on_the_reference_probe(setup):
    out = setup["bk"].stage_rates(setup["reference"], setup["net"])
    assert out["acrylamide_lane"]["n_rate_values"] > 50_000
    assert out["acrylamide_lane"]["worst_relative"] < 1e-12, out["acrylamide_lane"]
    assert out["trunk_lane_terms"]["n_rate_values"] > 50_000
    assert out["trunk_lane_terms"]["worst_relative"] < 1e-12, out["trunk_lane_terms"]


def test_rate_constants_equal_the_engine_live(setup):
    """Against the engine's functions in this process, at a point and a draw the reference does not hold."""
    from src.kinetic_core import acrylamide as lane
    from src.kinetic_core import acrylamide_conditions as ac
    from src.kinetic_core import engine

    bk, net, network = setup["bk"], setup["net"], setup["network"]
    params = engine.core_parameters(engine.ACRYLAMIDE)
    params["k_lys_glc"] = dataclasses.replace(params["k_lys_glc"], k_ref=0.004, ea_kj_mol=44.0)
    draw = {"aw_scale": 0.77, "aw_elimination_scale": 1.1, "ph_exponent_formation": 0.26,
            "ph_exponent_elimination": 0.155, "ea_int1": 101.0, "ea_asn_glc": 98.0}
    import build_jax_trunk_acrylamide_reference as refgen

    theta = bk.theta_of(refgen.effective_theta(params, draw), net)
    for t_c, ph, aw in ((151.3, 5.3, 0.77), (97.0, 8.6, 0.30), (205.0, None, 0.95), (133.0, 6.8, None)):
        applied, _ = ac.apply(params, SimpleNamespace(ph=ph, water_activity=aw), **refgen.apply_kwargs(draw))
        expect = lane.acrylamide_rate_constants_at(applied, t_c + 273.15)
        ctx = network.RateCtx(t_c=jnp.asarray(t_c), t_k=jnp.asarray(t_c + 273.15),
                              ph_label=jnp.asarray(setup["trunk"].nan_or(ph)),
                              aw=jnp.asarray(setup["trunk"].nan_or(aw)))
        k = np.asarray(net.rate_constants(theta, ctx))
        e = np.array([expect[r] for r in net.reactions])
        np.testing.assert_allclose(k, e, rtol=1e-12, atol=0.0)


def test_trajectories_match_the_engine_on_the_acrylamide_panel_pots(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    recs = [r for r in reference["pots"] if r["panel"] != "synthetic"]
    assert len({r["id"] for r in recs}) == 6 and len(recs) == 18
    cols = bk.engine_columns_of_species(reference, net)
    s = max(len(r["refined_segments"]) for r in recs)
    run = jax.jit(lambda th, pot: network.integrate(net, th, pot, rtol=1e-10, atol=1e-16).ys)
    for rec in recs:
        pot = bk.pot_of(rec, net, s)
        ys = np.asarray(run(bk.theta_of(rec["theta"], net), pot))
        c = bk.compare(ys, rec, pot["y0"], net, cols)
        assert c["max_dex"] < 1e-6, (rec["id"], rec["variant"], c)
        assert c["max_abs_below_floor_share_of_y0"] < 1e-9, (rec["id"], rec["variant"], c)
        e0 = net.element_matrix @ pot["y0"]
        np.testing.assert_allclose(net.element_matrix @ ys[-1], e0, rtol=1e-9, atol=1e-12)


def test_a_staircase_profile_with_every_step_matches_the_engine(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    rec = next(r for r in reference["pots"] if r["panel"] == "synthetic")
    assert len(rec["refined_segments"]) > 100
    pot = bk.pot_of(rec, net)
    ys = np.asarray(jax.jit(lambda th, p: network.integrate(net, th, p, rtol=1e-10, atol=1e-16).ys)(
        bk.theta_of(rec["theta"], net), pot))
    c = bk.compare(ys, rec, pot["y0"], net, bk.engine_columns_of_species(reference, net))
    assert c["max_dex"] < 1e-6, c
    # every step carries flux in this pot (the profile exists to exercise all 42)
    k = np.asarray(net.rate_constants(net.theta0(), network.RateCtx(t_c=180.0, t_k=453.15, ph_label=5.0, aw=0.8)))
    rates = k * np.asarray(network.mass_action(jnp.asarray(ys[len(ys) // 2]), net))
    assert np.all(rates > 0.0), [r for r, v in zip(net.reactions, rates) if not v > 0.0]


def test_exact_jacobian_equals_central_differences(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    pot_id = "mp_holdout_glucose_asparagine_180C_30min_Chang2021"
    fd = reference["jacobian_central_fd"][pot_id]
    rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
    obs_idx = jnp.asarray([net.index[s] for s in reference["observable_species"]])
    rows = jnp.asarray(fd["rows"])

    def observe(ys):
        return jnp.log10(ys[rows][:, obs_idx].T.ravel())

    pot = bk.pot_of(rec, net)
    theta = network.theta_vector(net, bk.theta_of(fd["theta"], net))
    _, jac = network.jacobian(net, observe, rtol=1e-10, atol=1e-16)(theta, pot)
    jac = np.asarray(jac)
    keys = [k for _, k in net.coordinates]
    cols = [keys.index(c) for c in fd["columns"]]
    errors = bk.column_errors(jac[:, cols], np.asarray(fd["jacobian"]))
    assert errors["worst"] < 1e-4, errors
    # against the port's own central differences on a few columns
    f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), pot,
                                                    rtol=1e-12, atol=1e-20).ys))
    h = fd["step"]
    for name in ("log10_k:k_int1_acr", "ea:k_acr_dp", "ph_exponent_formation", "log10_k:k_glc_tdg"):
        c = keys.index(name)
        e = jnp.zeros(len(theta)).at[c].set(h)
        col = (np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h)
        np.testing.assert_allclose(jac[:, c], col, rtol=1e-5, atol=1e-8, err_msg=name)


def test_panel_predictions_at_fit_tolerance_equal_the_engine(setup):
    out = setup["bk"].stage_residuals(setup["reference"], setup["net"])
    for row in out["rows"]:
        assert row["n_rows"] == 8
        assert row["max_abs_diff_dex"] < 1e-4, row


def test_cf106_refusal_and_uncarried_species_are_refused(setup):
    acrylamide, network, net = setup["acrylamide"], setup["network"], setup["net"]
    for peak in (61.0, 100.0, 119.9):
        with pytest.raises(network.UnsupportedConfiguration, match="CF-106"):
            acrylamide.make_pot(net, {"Glc": 10.0, "Asn": 10.0}, [(5.0, 40.0), (30.0, peak)], ph=6.0,
                                water_activity=None)
    for peak in (60.0, 120.0, 180.0):
        acrylamide.make_pot(net, {"Glc": 10.0, "Asn": 10.0}, [(30.0, peak)], ph=6.0, water_activity=None)
    with pytest.raises(network.UnsupportedConfiguration, match="Rib"):
        acrylamide.make_pot(net, {"Glc": 10.0, "Rib": 1.0}, [(30.0, 160.0)], ph=6.0, water_activity=None)
    assert net.ph_mode == "prescribed"
    pot = acrylamide.make_pot(net, {"Glc": 10.0, "Asn": 10.0}, [(30.0, 160.0)], ph=6.0, water_activity=None)
    with pytest.raises(ValueError, match="lacks"):
        network.integrate(net, net.theta0(), {"y0": pot["y0"], "segments": pot["segments"]})


def test_compose_refuses_acrylamide_without_trunk_and_a_second_owner_of_cysteine(setup):
    network, net = setup["network"], setup["net"]
    trunk_m, acr = net.modules
    with pytest.raises(ValueError, match="requires"):
        network.compose([acr])
    other = dataclasses.replace(acr, name="sulfur_stand_in", requires=("trunk",), reactions=(), stoichiometry={},
                                reactants={}, params={}, species=("Cys",), elements={"Cys": acr.elements["Cys"]})
    with pytest.raises(ValueError, match="defined by both"):
        network.compose([trunk_m, acr, other])
    bad = dict(acr.stoichiometry)
    bad["a_int1_acr"] = {"INT1": -1.0, "ACR": 1.0, "FRAG_C": 6.0, "FRAG_N": 1.0}
    with pytest.raises(ValueError, match="does not balance"):
        network.compose([trunk_m, dataclasses.replace(acr, stoichiometry=bad)])


def test_trunk_lane_switches_are_ported_or_refused(setup):
    from src.kinetic_core import config

    trunk, network = setup["trunk"], setup["network"]
    shipped = config.Config.shipped()
    on = trunk.build(shipped, lane="trunk")
    assert "ph_exponent_exit:k_tdg_fa" in on.params and on.pot_inputs == ("ph_label", "aw")
    off = trunk.build(shipped.replace({"trunk.three_deoxy_exit_ph_term": False}), lane="trunk")
    assert not any(k.startswith("ph_exponent_exit:") for k in off.params)
    exit_map = dict(shipped["trunk.three_deoxy_exit_ph"])
    outside = dict(exit_map, k_go_sink=(0.5, (0.4, 0.6), "a step of the dicarbonyl block"))
    with pytest.raises(network.UnsupportedConfiguration, match="k_go_sink"):
        trunk.build(shipped.replace({"trunk.three_deoxy_exit_ph": outside}), lane="trunk")
    with pytest.raises(ValueError):
        trunk.build(shipped, lane="sulfur")
    acr_trunk = trunk.build(shipped)
    assert acr_trunk.pot_inputs == () and not any(k in acr_trunk.params for k in ("aw_scale", "ph_exponent"))
    with pytest.raises(network.UnsupportedConfiguration, match="atmosphere"):
        trunk.check_atmosphere("argon")
    trunk.check_atmosphere("air")
    trunk.check_atmosphere(None)


def test_the_trunk_is_the_shared_prefix_of_the_other_two_trunks(setup):
    sets = setup["reference"]["trunk_sets"]
    assert sets["network.REACTIONS"] == 26 and sets["acrylamide.ACRYLAMIDE_REACTIONS"] == 16
    assert sets["network.TRUNK_REACTIONS starts with REACTIONS"]
    assert sets["sulfur.FULL_REACTIONS starts with REACTIONS + FED_3DEOXY_REACTIONS"]
    assert math.isfinite(float(len(setup["net"].reactions))) and len(setup["net"].reactions) == 42
