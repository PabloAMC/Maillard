"""src/kinetic_core_jax: the A2 contract v2 runtime (network.py) and its first module, the lipid lane (lipid.py).

Skipped unless jax and diffrax import, so the core env's CI never runs (or breaks on) them; run them in the separate
``maillard-jax`` env (environment-jax.yml). They pin what results/validation/jax_lipid_module.md claims:

* every rate constant equals the engine's (k_looh_decomp_per_min, DecadienalBreakage.k_loss_per_min, the per-cell slate
  yields, the 1981 ratios) to 1e-12 relative, on the reference probe and live;
* every lipid panel pot, at its point and both declared-band corners, and a staircase profile, match the engine's
  states after every segment to 1e-6 dex at tight tolerance; carbon closes;
* the forward-mode Jacobian equals the engine's central differences;
* compose refuses duplicate and unbalanced definitions; build refuses an unported switch; a pH-coupled network is
  refused by integrate; padding is the identity and a batch equals its single solves.
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
REFERENCE = ROOT / "results" / "validation" / "jax_lipid_module_reference.json"


@pytest.fixture(scope="module")
def setup():
    if str(GEN) not in sys.path:
        sys.path.insert(0, str(GEN))
    import build_jax_lipid_module as bk
    from src.kinetic_core_jax import lipid, network

    reference = json.loads(REFERENCE.read_text())
    net = network.compose([lipid.build()])
    return {"bk": bk, "reference": reference, "net": net, "lipid": lipid, "network": network}


def test_rate_constants_equal_the_engine_on_the_reference_probe(setup):
    out = setup["bk"].stage_rates(setup["reference"], setup["net"])
    assert out["n_rate_values"] > 4000
    assert out["worst_relative"] < 1e-12, out


def test_rate_constants_equal_the_engine_live(setup):
    """Against the engine's functions in this process, at a point the reference does not hold."""
    from src.kinetic_core import lipid as lane
    from src.kinetic_core import parameters_lipid as pl

    net, network, lipid = setup["net"], setup["network"], setup["lipid"]
    mod = net.modules[0]
    t_c, q10, donor = 151.3, 2.1, 0.37
    theta = {k: jnp.asarray(float(v)) for k, v in mod.params.items()}
    theta.update(q10=jnp.asarray(q10), donor=jnp.asarray(donor))
    k = np.asarray(mod.rate_constants(theta, network.RateCtx(t_c=jnp.asarray(t_c), t_k=jnp.asarray(t_c + 273.15))))
    k_looh = pl.k_looh_decomp_per_min(t_c, q10)
    k_brk = lane.DECADIENAL_BREAKAGE.k_loss_per_min(t_c)
    cells = lipid.cell_yields(donor)
    ix = {r: j for j, r in enumerate(net.reactions)}
    assert k[ix["lipid_decomp_LOOH_OL"]] == pytest.approx(k_looh, rel=1e-13)
    assert k[ix["lipid_slate_LOOH_9_tt_DECADIENAL"]] == pytest.approx(k_looh * cells["LOOH_9_tt"]["DECADIENAL"],
                                                                      rel=1e-13)
    assert k[ix["lipid_slate_LOOH_13_ct_HEXANAL"]] == pytest.approx(k_looh * cells["LOOH_13_ct"]["HEXANAL"],
                                                                    rel=1e-13)
    assert k[ix["lipid_decadienal_loss"]] == pytest.approx(k_brk, rel=1e-13)
    assert k[ix["lipid_decadienal_hexanal"]] == pytest.approx(
        k_brk * lane.DECADIENAL_BREAKAGE.share_hexanal, rel=1e-13)
    below = np.asarray(mod.rate_constants(theta, network.RateCtx(t_c=jnp.asarray(99.0), t_k=jnp.asarray(372.15))))
    assert below[ix["lipid_decadienal_loss"]] == 0.0   # the engine's breakage window


def test_trajectories_match_the_engine_on_the_lipid_panel_pots(setup):
    bk, net, reference = setup["bk"], setup["net"], setup["reference"]
    network = setup["network"]
    recs = [r for r in reference["pots"] if r["panel"] != "synthetic"]
    assert len({r["id"] for r in recs}) == 8 and len(recs) == 24
    s = max(len(r["refined_segments"]) for r in recs)
    run = jax.jit(lambda th, pot: network.integrate(net, th, pot, rtol=1e-10, atol=1e-16).ys)
    for rec in recs:
        pot = bk.pot_of(rec, net, s)
        ys = np.asarray(run(bk.theta_of(rec["theta"]), pot))
        c = bk.compare(ys, rec, pot["y0"], net.species)
        assert c["max_dex"] < 1e-6, (rec["id"], rec["variant"], c)
        assert c["max_abs_below_floor_share_of_y0"] < 1e-9, (rec["id"], rec["variant"], c)
        carbon = net.element_matrix[0]
        assert float(carbon @ ys[-1]) == pytest.approx(float(carbon @ pot["y0"]), rel=1e-10)


def test_a_staircase_profile_matches_the_engine(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    rec = next(r for r in reference["pots"] if r["panel"] == "synthetic")
    assert len(rec["refined_segments"]) > 100
    pot = bk.pot_of(rec, net)
    ys = np.asarray(jax.jit(lambda th, p: network.integrate(net, th, p, rtol=1e-10, atol=1e-16).ys)(
        bk.theta_of(rec["theta"]), pot))
    c = bk.compare(ys, rec, pot["y0"], net.species)
    assert c["max_dex"] < 1e-6, c


def test_exact_jacobian_equals_the_engines_central_differences(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    pot_id = "pea_isolate_uht_140C_Trikusuma2019"
    fd = reference["jacobian_central_fd"][pot_id]
    rec = next(r for r in reference["pots"] if r["id"] == pot_id and r["variant"] == "point")
    obs_idx = jnp.asarray([net.index[s] for s in reference["observable_species"]])
    rows = jnp.asarray(fd["rows"])

    def observe(ys):
        return jnp.log10(ys[rows][:, obs_idx].T.ravel())

    theta = jnp.asarray([fd["theta"][k] for _, k in net.coordinates])
    _, jac = network.jacobian(net, observe, rtol=1e-10, atol=1e-16)(theta, bk.pot_of(rec, net))
    jac = np.asarray(jac)
    errors = bk.column_errors(jac[:, bk.engine_columns(net, fd)], np.asarray(fd["jacobian"]))
    assert errors["worst"] < 1e-6, errors
    # log10_k_looh, which the engine cannot move: against the port's own central difference
    c = [k for _, k in net.coordinates].index("log10_k_looh")
    f = jax.jit(lambda v: observe(network.integrate(net, network.theta_from_vector(net, v), bk.pot_of(rec, net),
                                                    rtol=1e-12, atol=1e-18).ys))
    h = 1e-5
    e = jnp.zeros(len(theta)).at[c].set(h)
    fd_col = (np.asarray(f(theta + e)) - np.asarray(f(theta - e))) / (2 * h)
    np.testing.assert_allclose(jac[:, c], fd_col, rtol=1e-6, atol=1e-8)


def test_compose_refuses_duplicates_and_unbalanced_reactions(setup):
    network, net = setup["network"], setup["net"]
    mod = net.modules[0]
    with pytest.raises(ValueError, match="duplicate module"):
        network.compose([mod, mod])
    with pytest.raises(ValueError, match="defined by both"):
        network.compose([mod, dataclasses.replace(mod, name="lipid_copy")])
    key = "lipid_decadienal_loss"
    bad = dict(mod.stoichiometry)
    bad[key] = {"DECADIENAL": -1.0, "LIPID_FRAG_C": 9.0}
    with pytest.raises(ValueError, match="does not balance"):
        network.compose([dataclasses.replace(mod, stoichiometry=bad)])


def test_build_reads_the_config_it_is_given_and_refuses_an_unported_switch(setup):
    from src.kinetic_core import config

    lipid, network = setup["lipid"], setup["network"]
    shipped = config.Config.shipped()
    sink = shipped["lipid.covalent_sink"]
    with pytest.raises(network.UnsupportedConfiguration):
        lipid.build(shipped.replace({"lipid.covalent_sink": dataclasses.replace(sink, enabled=True)}))
    brk = shipped["lipid.decadienal_breakage"]
    off = lipid.build(shipped.replace({"lipid.decadienal_breakage": dataclasses.replace(brk, enabled=False)}))
    assert not any("decadienal" in r for r in off.reactions)
    q10 = shipped["lipid.q10_assumption"]
    other = lipid.build(shipped.replace({"lipid.q10_assumption": dataclasses.replace(q10, default=3.0)}))
    assert other.params["q10"] == 3.0 and other.config_digest != lipid.build(shipped).config_digest


def test_a_ph_coupled_network_is_refused(setup):
    network, mod = setup["network"], setup["net"].modules[0]
    net = network.compose([dataclasses.replace(mod, needs_ph=True)])
    with pytest.raises(network.UnsupportedConfiguration):
        network.integrate(net, net.theta0(), {"y0": np.zeros(net.n_species), "segments": [[1.0, 100.0]]})


def test_padding_is_the_identity_and_a_batch_equals_its_single_solves(setup):
    bk, net, reference, network = setup["bk"], setup["net"], setup["reference"], setup["network"]
    rec = next(r for r in reference["pots"] if r["id"] == "external_validation_bi_2020_roasted_pea_hexanal"
               and r["variant"] == "point")
    run = jax.jit(lambda th, p: network.integrate(net, th, p, rtol=1e-9, atol=1e-15).ys)
    short = np.asarray(run(net.theta0(), bk.pot_of(rec, net)))
    padded = np.asarray(run(net.theta0(), bk.pot_of(rec, net, len(rec["refined_segments"]) + 3)))
    np.testing.assert_allclose(padded[len(short):], np.repeat(short[-1:], 3, axis=0), rtol=1e-14, atol=0.0)
    draws = [net.theta0(), net.theta0()]
    draws[1]["lipid"]["q10"] = jnp.asarray(3.1)
    draws[1]["lipid"]["donor"] = jnp.asarray(0.2)
    batch = np.asarray(network.batched(net, pot_axis=None, rtol=1e-9, atol=1e-15)(
        network.stack_theta(draws), bk.pot_of(rec, net)))
    for i, th in enumerate(draws):
        one = np.asarray(run(th, bk.pot_of(rec, net)))
        mask = one > 1e-9 * float(np.max(one))
        assert np.max(np.abs(np.log10(batch[i][mask]) - np.log10(one[mask]))) < 1e-9
    assert math.isfinite(float(batch.sum()))
