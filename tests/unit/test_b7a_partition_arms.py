"""B7a's rival arms (2026-10-31; results/validation/kinetic_core_dpo_partition_prereg.md sec. 8, as revised after its
design review): R-S, the cysteine Amadori compound's enolisation split under a saturating base law
(CysAmadoriRoute.dpo_share_at); R-T, 3-deoxypentosone's own decay at the trunk's 3-DG -> formic acid constant
(TDP_OWN_DECAY); B7D, the coded r_tdp_mgo freed. All off by default and then bit-identical; the joint fit's install()
leaves none set on a job that does not carry it."""
from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

import pytest

from src.kinetic_core import engine, parameters as trunk, parameters_sulfur as ps, sulfur

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))
import fit_joint_deoxyosone_thiol as joint  # noqa: E402

#: (T K, in-situ pH, label pH): a spread of the states the B7a pots visit
STATES = ((373.15, 5.5, 5.5), (393.15, 7.0, 6.6), (413.15, 8.0, 7.4))


def _k(t_k, ph, ph_label):
    return sulfur.sulfur_rate_constants_at(engine.core_parameters(engine.SULFUR), t_k, ph, ph_label=ph_label)


def _odds(share):
    return (1.0 - share) / share


def test_the_split_is_dpo_share_exactly_under_the_null_law_and_saturating_otherwise():
    route = ps.CysAmadoriRoute()
    assert (route.split_law, route.split_ph_c, route.dpo_share_ref_ph) == ("none", 7.5, 7.0)
    for ph in (3.0, 5.5, 7.0, 8.3):
        assert route.dpo_share_at(ph) == route.dpo_share      # identity, no arithmetic on the share
    assert dataclasses.replace(route, dpo_share=0.123456789, split_ph_c=6.0).dpo_share_at(6.1) == 0.123456789
    law = dataclasses.replace(route, split_law="base_saturating", dpo_share=0.30, split_ph_c=7.5)
    assert law.dpo_share_at(7.0) == pytest.approx(0.30, rel=1e-14)       # dpo_share is the share at pH 7
    for ph in (4.0, 5.5, 8.0, 9.5):
        odds = _odds(0.30) * (1.0 + 10.0 ** (7.0 - 7.5)) / (1.0 + 10.0 ** (ph - 7.5))
        assert law.dpo_share_at(ph) == pytest.approx(1.0 / (1.0 + odds), rel=1e-12)
    # flat in acid (saturated), one decade of TDP : DPO odds per pH unit well above pH_c
    assert _odds(law.dpo_share_at(3.0)) / _odds(law.dpo_share_at(4.0)) == pytest.approx(1.0, abs=2e-3)
    assert _odds(law.dpo_share_at(11.0)) / _odds(law.dpo_share_at(12.0)) == pytest.approx(10.0, rel=1e-3)
    with pytest.raises(ValueError):
        dataclasses.replace(route, split_law="logit").dpo_share_at(7.0)


@pytest.mark.parametrize("t_k, ph, ph_label", STATES)
def test_every_switch_at_default_leaves_the_derived_rates_bit_identical(monkeypatch, t_k, ph, ph_label):
    monkeypatch.setattr(ps, "CYS_AMADORI_ROUTE", dataclasses.replace(ps.CYS_AMADORI_ROUTE, enabled=True))
    base = _k(t_k, ph, ph_label)
    route = ps.CYS_AMADORI_ROUTE
    decay = route.k_decay(t_k, ph)
    assert base["r_cama_dpo"] == route.dpo_share * decay
    assert base["r_cama_tdp"] == (1.0 - route.dpo_share) * decay
    assert base["r_osone_decay_tdp"] == base["r_osone_decay_dpo"]   # the shared k_osone_decay
    assert base["r_tdp_mgo"] == 0.0
    # the inert fields of either switch cannot move anything while it is off / at its null
    monkeypatch.setattr(ps, "TDP_OWN_DECAY", ps.TdpOwnDecay(enabled=False, k_ref_per_min=99.0, ea_kj_mol=1.0))
    monkeypatch.setattr(ps, "CYS_AMADORI_ROUTE", dataclasses.replace(route, split_ph_c=5.1, dpo_share_ref_ph=3.3))
    assert _k(t_k, ph, ph_label) == base


def test_the_tdp_exit_defaults_are_the_trunks_k_tdg_fa():
    k_tdg_fa = trunk.MARTINS_M4["k_tdg_fa"]
    d = ps.TdpOwnDecay()
    assert not d.enabled and not ps.TDP_OWN_DECAY.enabled
    assert (d.k_ref_per_min, d.ea_kj_mol) == (k_tdg_fa.k_ref, k_tdg_fa.ea_kj_mol)
    assert d.t_ref_c + 273.15 == pytest.approx(trunk.T_REF_K, abs=1e-12)
    for t_c in (100.0, 120.0, 140.0):
        assert d.k_at(t_c + 273.15) == pytest.approx(k_tdg_fa.k_at(t_c + 273.15), rel=1e-12)


@pytest.mark.parametrize("t_c", (100.0, 120.0, 140.0))
def test_the_tdp_exit_runs_at_its_constant_and_its_sensitivity(monkeypatch, t_c):
    t_k = t_c + 273.15
    monkeypatch.setattr(ps, "TDP_OWN_DECAY", ps.TdpOwnDecay(enabled=True))
    k = _k(t_k, 6.0, 6.0)
    assert k["r_osone_decay_tdp"] == pytest.approx(trunk.MARTINS_M4["k_tdg_fa"].k_at(t_k), rel=1e-12)
    assert k["r_osone_decay_tdp"] != k["r_osone_decay_dpo"]               # DPO keeps the shared decay
    assert _k(t_k, 8.0, 8.0)["r_osone_decay_tdp"] == k["r_osone_decay_tdp"]   # no pH term
    # Martins 2005b's sensitivity pairing: the prereg's arithmetic check, 1.0e-2 / 1.6e-2 / 2.5e-2 1/min
    monkeypatch.setattr(ps, "TDP_OWN_DECAY", ps.TdpOwnDecay(enabled=True, k_ref_per_min=1.0e-2, ea_kj_mol=29.6))
    expected = {100.0: 1.0e-2, 120.0: 1.6e-2, 140.0: 2.5e-2}[t_c]
    exact = 1.0e-2 * math.exp(-29.6 / ps.R_KJ * (1.0 / t_k - 1.0 / 373.15))
    assert _k(t_k, 6.0, 6.0)["r_osone_decay_tdp"] == pytest.approx(exact, rel=1e-12)
    assert exact == pytest.approx(expected, rel=0.03)


def test_the_split_moves_the_cama_branches_with_the_in_situ_ph(monkeypatch):
    monkeypatch.setattr(ps, "CYS_AMADORI_ROUTE", dataclasses.replace(
        ps.CYS_AMADORI_ROUTE, enabled=True, split_law="base_saturating", dpo_share=0.30, split_ph_c=6.5))
    k = _k(393.15, 8.0, 7.0)
    share = k["r_cama_dpo"] / (k["r_cama_dpo"] + k["r_cama_tdp"])
    assert share == pytest.approx(ps.CYS_AMADORI_ROUTE.dpo_share_at(8.0), rel=1e-12)   # in-situ, not label
    assert share > 0.30


def _job(arm, monkeypatch, **flags):
    monkeypatch.setattr(joint, "TDP_EXIT_STATIC", joint.tdp_exit_static(joint.ARMS[arm], **flags))
    keys, _kinds, _lo, _hi, start = joint.coordinates(False, joint.ARMS[arm])
    theta = dict(zip(keys, start))
    if "dpo_share7" in theta:
        theta.update(dpo_share7=0.45, split_ph_c=6.2)
    return dict(joint.trial([theta[k] for k in keys], keys), strecker_share=0.042)


def test_the_arms_carry_their_coordinates_and_static_switch(monkeypatch):
    keys, kinds, lower, upper, start = joint.coordinates(False, joint.ARMS["B7ST"])
    i, j = keys.index("dpo_share7"), keys.index("split_ph_c")
    assert (kinds[i], lower[i], upper[i], start[i]) == ("acid_yield", 0.02, 0.98, 0.30)
    assert (kinds[j], lower[j], upper[j], start[j]) == ("pka", 5.0, 10.0, 7.5)
    keys_d, kinds_d, lower_d, upper_d, start_d = joint.coordinates(False, joint.ARMS["B7D"])
    m = keys_d.index("log10_k_tdp_mgo")
    assert (kinds_d[m], lower_d[m], upper_d[m], start_d[m]) == ("log10_offset", -8.0, -1.0, -5.08)
    n_null = len(joint.coordinates(False, joint.ARMS["R3En"])[0])
    assert len(joint.coordinates(False, joint.ARMS["B7T"])[0]) == n_null and len(keys_d) == n_null + 1
    assert joint.ARMS["B7ST"] - joint.ARMS["R3En"] == {"split_ph", "tdp_exit"}
    assert joint.ARMS["B7D"] - joint.ARMS["R3En"] == {"tdp_mgo"}
    assert _job("B7T", monkeypatch)["tdp_exit"] == dataclasses.asdict(ps.TdpOwnDecay(enabled=True))
    assert _job("B7T", monkeypatch, k_ref_per_min=1e-2, ea_kj_mol=84.0)["tdp_exit"]["ea_kj_mol"] == 84.0
    assert "tdp_exit" not in _job("B7S", monkeypatch) and "cama_split" not in _job("R3En", monkeypatch)
    routes = _job("B7D", monkeypatch)["routes"]
    assert routes["k_tdp_mgo_100c_ph5"] == pytest.approx(10.0 ** -5.08)
    assert routes["ea_tdp_mgo_kj_mol"] == ps.PentoseMgoRoutes().ea_tdp_mgo_kj_mol == 100.0
    assert _job("R3En", monkeypatch)["routes"]["k_tdp_mgo_100c_ph5"] == 0.0
    with pytest.raises(SystemExit):
        joint.tdp_exit_static(joint.ARMS["B7S"], k_ref_per_min=1e-2)


def test_install_sets_every_switch_and_resets_them_on_a_job_without_them(monkeypatch):
    job = _job("B7ST", monkeypatch, ea_kj_mol=43.9)
    job["routes"] = dict(job["routes"], k_tdp_mgo_100c_ph5=1e-5)
    joint.install(job)
    assert ps.TDP_OWN_DECAY.enabled and ps.TDP_OWN_DECAY.ea_kj_mol == 43.9
    route = ps.CYS_AMADORI_ROUTE
    assert (route.enabled, route.split_law, route.dpo_share, route.split_ph_c) == (True, "base_saturating", 0.45, 6.2)
    assert ps.PENTOSE_MGO_ROUTES.k_tdp_mgo_100c_ph5 == 1e-5
    joint.install(_job("R3En", monkeypatch))
    assert ps.TDP_OWN_DECAY == ps.TdpOwnDecay()
    assert ps.CYS_AMADORI_ROUTE == dataclasses.replace(ps.CysAmadoriRoute(), enabled=True)
    assert ps.PENTOSE_MGO_ROUTES.k_tdp_mgo_100c_ph5 == 0.0
