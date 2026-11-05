"""Plan step M1b part (i)'s build (src/kinetic_core/ph_balance_rivals.py, scripts/generators/fit_m1b_ph_balance.py;
prereg results/validation/kinetic_core_m1b_ph_balance_pentose_steps_prereg.md v3.6): every rival off by default and
never entered then, each rival's direction on a toy pot, the B0 corners, the ARP forcing, B3's branch rule, the pH
scale guard, the measurement model's pieces and the refusals."""
from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src.kinetic_core import config, engine, observations, ph_state, sulfur  # noqa: E402
from src.kinetic_core import ph_balance_rivals as pbr  # noqa: E402
from src.kinetic_core import predictive_densities as pdens  # noqa: E402
from src.kinetic_core.ph_state import BufferSpec  # noqa: E402

NEW_KEYS = ("ph.amine_survival", "ph.acid_lag", "ph.acid_supply", "ph.arp_forcing")
WATER = BufferSpec(kind="none", declared=True, source="test water")
TIMES = np.array([0.0, 20.0, 40.0])
BREAKPOINTS = (5.0, 6.0, 6.5)


@pytest.fixture(scope="module")
def shipped():
    return engine.core_parameters(engine.SULFUR), engine.core_ph_drift()


def _run(changes, shipped, initial=None, drift=None, times=TIMES, label=7.0):
    params, d0 = shipped
    cfg = config.Config.shipped().replace({"ph.coupling": "inline", **changes})
    with config.applied(cfg):
        run = sulfur.integrate_sulfur(params, 393.15, initial or {"ARP": 5.0, "OX": 1.0}, times, ph=label,
                                      buffer_spec=WATER, ph_drift=drift or d0, rtol=1e-8, atol=1e-14)
    sid = run.metadata["strong_ion_difference_mmol_L"]
    pools = (run.metadata.get("ph_balance_rivals") or {}).get("pools") or {}
    out = []
    for i in range(len(times)):
        st = {k: float(run.series(k)[i]) for k in sulfur.SULFUR_STATE_KEYS}
        st.update(pbr.balance_pools({k: v[i] for k, v in pools.items()}))
        out.append(observations.ph_cooled(st, sid, drift or d0, WATER))
    return run, np.array(out)


# --- the shipped default ------------------------------------------------------------------------------------------

def test_the_new_keys_ship_off_and_match_the_module():
    shipped = config.Config.shipped()
    for name in NEW_KEYS:
        value = shipped[name]
        assert value.enabled is False
        assert value == getattr(pbr, config.KEY_BY_NAME[name].bindings[0][1])
    with config.applied(shipped):
        assert not pbr.any_active()


def test_the_default_never_enters_the_rival_path(shipped, monkeypatch):
    """Bit-identity by construction: with every key off, integrate_sulfur runs the shipped inline code; the rival
    integrator is not called (it would raise here)."""
    def boom(*a, **k):
        raise AssertionError("the rival path was entered with every switch off")

    monkeypatch.setattr(pbr, "integrate_inline", boom)
    run, ph = _run({}, shipped)
    assert "ph_balance_rivals" not in run.metadata
    assert np.all(np.isfinite(ph))


def test_the_inventory_without_a_rival_pool_is_unchanged():
    drift = engine.core_ph_drift()
    state = {"ARP": 3.0, "CBX": 2.0, "ACID": 1.0, "Cys": 1.0}
    base = ph_state.titratable_inventory(state, drift, WATER)
    with_zero = ph_state.titratable_inventory({**state, "NH4": 0.0, "RET_AA": 0.0}, drift, WATER)
    assert base == with_zero
    more = ph_state.titratable_inventory({**state, "NH4": 1.0}, drift, WATER)
    assert more[-1][0] is ph_state.AMMONIUM or any(g is ph_state.AMMONIUM for g, _ in more)


# --- each rival's direction on a toy pot (an unbuffered 5 mM ARP pot, 120 C) ----------------------------------------

def test_b1_amine_survival_raises_the_ph(shipped):
    _, b0 = _run({}, shipped)
    _, b1 = _run({"ph.amine_survival": pbr.AmineSurvival(True, 0.0, 1.0)}, shipped)
    assert b1[-1] > b0[-1] + 0.5


def test_b1_at_phi_zero_reproduces_b0(shipped):
    _, b0 = _run({}, shipped)
    _, b1 = _run({"ph.amine_survival": pbr.AmineSurvival(True, 0.0, 0.0)}, shipped)
    assert np.max(np.abs(b1 - b0)) <= 1e-6


def test_b2_returned_amino_acid_and_lag_make_less_acid(shipped):
    _, b0 = _run({}, shipped)
    run, ret = _run({"ph.acid_lag": pbr.AcidLag(True, 1.0, 10.0, 100.0)}, shipped)
    assert ret[-1] > b0[-1] + 0.3
    assert run.metadata["ph_balance_rivals"]["pools"]["RET_AA"][-1] > 0.0
    _, slow = _run({"ph.acid_lag": pbr.AcidLag(True, 0.0, 0.01, 100.0)}, shipped)
    assert slow[1] > b0[1]          # the acid arrives late
    _, fast = _run({"ph.acid_lag": pbr.AcidLag(True, 0.0, 1e6, 100.0)}, shipped)
    assert np.max(np.abs(fast - b0)) <= 1e-6


def test_b3_acid_supply_makes_more_acid_and_holds_dpo_total(shipped):
    d0 = shipped[1]
    no_lump = ph_state.PhDrift(acid_yield=0.0, arp_amine_pka=d0.arp_amine_pka)
    off = pbr.AcidSupply(True, 0.0, 0.0, 0.0, BREAKPOINTS, 0.1)
    on = pbr.AcidSupply(True, 0.8, 0.9, 1.0, BREAKPOINTS, 0.1)
    _, y0 = _run({"ph.acid_supply": off}, shipped, drift=no_lump)
    run, y8 = _run({"ph.acid_supply": on}, shipped, drift=no_lump)
    assert y8[-1] < y0[-1] - 0.1
    diag = run.metadata["ph_balance_rivals"]["diagnostics"]
    assert diag["dpo_total_rel_dev_max"] <= 1e-9
    assert run.metadata["ph_balance_rivals"]["pools"]["LAA"][-1] >= 0.0
    added = dataclasses.replace(on, branch_form="added_sink", ph_in_situ_breakpoints=(3.0, 4.0, 4.5))
    run_a, ya = _run({"ph.acid_supply": added}, shipped, drift=no_lump)
    assert run_a.metadata["ph_balance_rivals"]["diagnostics"]["added_sink_cap_reached"]   # y = 1 above pH 4.5
    assert np.all(np.isfinite(ya))


def test_b3_share_is_linear_in_the_in_situ_ph_and_flat_outside():
    s = pbr.AcidSupply(True, 0.2, 0.6, 0.9, (5.0, 6.0, 6.5), 0.1)
    assert s.share(4.0) == pytest.approx(0.2)
    assert s.share(5.5) == pytest.approx(0.4)
    assert s.share(6.25) == pytest.approx(0.75)
    assert s.share(8.0) == pytest.approx(0.9)


def test_the_forcing_follows_the_measured_arp(shipped):
    times = (0.0, 20.0, 40.0)
    arp = (5.0, 3.2, 2.18)
    run, _ = _run({"ph.arp_forcing": pbr.ArpForcing(True, times, arp)}, shipped)
    assert np.allclose(run.series("ARP"), arp, rtol=1e-6)
    diag = run.metadata["ph_balance_rivals"]["diagnostics"]
    assert diag["forced_loss_never_negative"] and diag["forced_unresolvable_calls"] == 0


# --- refusals --------------------------------------------------------------------------------------------------------

def test_a_rival_runs_only_inline(shipped):
    params, d0 = shipped
    cfg = config.Config.shipped().replace({"ph.amine_survival": pbr.AmineSurvival(True, 0.5, 0.5)})
    with config.applied(cfg), pytest.raises(ValueError, match="inline"):
        sulfur.integrate_sulfur(params, 393.15, {"ARP": 5.0}, TIMES, ph=7.0, buffer_spec=WATER, ph_drift=d0)


def test_undeclared_combinations_and_bad_values_are_refused(shipped):
    with pytest.raises(ValueError, match="B12"):
        _run({"ph.acid_lag": pbr.AcidLag(True, 0.5), "ph.amine_survival": pbr.AmineSurvival(True, 0.5, 0.5)}, shipped)
    with pytest.raises(ValueError, match="order"):
        pbr.AcidSupply(True, 0.7, 0.6, 0.9, BREAKPOINTS)
    with pytest.raises(ValueError, match="breakpoints"):
        pbr.AcidSupply(True, 0.1, 0.6, 0.9, ())
    with pytest.raises(ValueError, match="monotone"):
        pbr.ArpForcing(True, (0.0, 20.0), (5.0, 5.5))
    with pytest.raises(ValueError, match="fraction"):
        pbr.AmineSurvival(True, 1.5, 0.0)
    with pytest.raises(ValueError):
        config.Config.shipped().replace({"ph.acid_lag": pbr.AcidLag(True, 0.5, 1e7)})   # outside its field band


def test_a_mixed_ph_scale_comparison_is_refused():
    with pytest.raises(pdens.MeasurementModelRefusal):
        pdens.ph_difference(pdens.PhReading(4.0, "in_situ", 120.0), pdens.PhReading(4.0, "cooled", 25.0))
    assert pdens.ph_difference(pdens.PhReading(4.2, "cooled"), pdens.PhReading(4.0, "cooled")) == pytest.approx(0.2)


# --- the generator: measurement model, rule, refusals -----------------------------------------------------------------

@pytest.fixture(scope="module")
def G():
    from tests.support import wave_generator
    with wave_generator("fit_m1b_ph_balance") as module:
        yield module
    for name in ("fit_m1b_ph_balance", "fit_m1_supply_combinations", "m1_preflight_checks"):
        sys.modules.pop(name, None)


def test_the_censored_likelihood_is_the_normal_interval(G):
    data = G.UnitData("C-Mey-Cys", within=(True, False), censored=True)
    f = np.array([[0.0, 0.5]])
    dl = np.array([[0.1, 0.1]])
    s = math.sqrt(2 * 0.05 ** 2 + 0.1 ** 2)
    from scipy.stats import norm
    p0 = norm.cdf(0.2 / s) - norm.cdf(-0.2 / s)
    p1 = 1.0 - (norm.cdf((0.2 - 0.5) / s) - norm.cdf((-0.2 - 0.5) / s))
    assert G.unit_loglik(data, f, dl)[0] == pytest.approx(math.log(p0) + math.log(p1))


def test_the_rule_credits_a_clearly_better_rival_and_refuses_the_reverse(G):
    rng = np.random.default_rng(1)
    data = {"C-Mey-Cys": G.UnitData("C-Mey-Cys", within=(True,) * 5, censored=True),
            "C-Mey-Gly": G.UnitData("C-Mey-Gly", within=(True,) * 5, censored=True),
            "C-Hou": G.UnitData("C-Hou", (5.0, 4.5), (0.05, 0.05)),
            "C-Hou-unforced": G.UnitData("C-Hou-unforced", (5.0, 4.5), (0.05, 0.05)),
            "C-Zhou-ARP": G.UnitData("C-Zhou-ARP", (4.0,), (0.05,)), "C-Zhou-Cys": G.UnitData("C-Zhou-Cys", (3.3,),
                                                                                              (0.05,))}

    def ens(name, mey, hou, zhou):
        S = 300
        j = rng.normal(size=(S, 1)) * 0.01
        f = {"C-Mey-Cys": np.full((S, 5), mey) + j, "C-Mey-Gly": np.zeros((S, 5)) + j,
             "C-Hou": np.tile(hou, (S, 1)) + j, "C-Hou-unforced": np.tile(hou, (S, 1)) + j,
             "C-Zhou-ARP": np.full((S, 1), zhou) + j, "C-Zhou-Cys": np.full((S, 1), 3.3) + j}
        d = {u: {k: np.full_like(a, 0.05) for k in G.OMEGAS} for u, a in f.items()}
        return G.Ensemble(name, "primary", [{}] * S, f, d)

    good = ens("B1", 0.0, np.array([5.0, 4.5]), 4.0)
    bad = ens("B0", -0.7, np.array([6.3, 5.5]), 4.0)
    e = {"B1": good, "B0": bad}
    res = G.rule_over_forms("B1", "B0", e, data, rng, require_resolved=False)
    assert res["favoured"], res["breaking_forms"]
    rev = G.rule_over_forms("B1", "B0", {"B1": bad, "B0": good}, data, rng, require_resolved=False)
    assert not rev["favoured"]


def test_draws_respect_the_bands_and_b3_order(G):
    draws = G.draw_coordinates("B3", "primary", 50, 3)
    assert len(draws) == 50
    for d in draws:
        assert d["y5"] <= d["y6"] <= d["y_hi"]
        assert 0.20 * 10 ** -0.3 - 1e-12 <= d["y5"] <= 0.56
        assert 0.025 <= d["k_iso"] <= 0.24
        assert 7.07 <= d["arp_amine_pka"] <= 8.18
    unwidened = G.draw_coordinates("B3", "delta_t_0", 20, 3)
    assert min(d["y_hi"] for d in unwidened) >= 0.80
    sohn = G.draw_coordinates("B1", "sohn_phi_cys", 20, 3)
    assert min(d["phi_cys"] for d in sohn) >= 0.76
    logy = G.draw_coordinates("B0", "yield_loguniform", 64, 3)
    assert min(d["acid_yield"] for d in logy) < 1e-2 and max(d["acid_yield"] for d in logy) <= 1.0


def test_score_and_fit_refuse_without_a_freeze_and_records(G, monkeypatch):
    monkeypatch.setattr(G, "read_record", lambda name: {})

    def no_freeze(commit):
        raise SystemExit("--freeze-commit is required")

    problems = G.score_refusals("", freeze_check=no_freeze)
    assert any("freeze" in p for p in problems) and any("screen" in p for p in problems)
    assert any("power" in p for p in problems)
    fp = G.fit_refusals("B0", "", freeze_check=no_freeze)
    assert any("incumbent" in p for p in fp) and any("score" in p for p in fp)


def test_score_runs_once(G, monkeypatch):
    rec = {"checks": {"passed": True}, "screen": {"complete": True}, "power": {"complete": True},
           "score": {"verdicts": {}}}
    monkeypatch.setattr(G, "read_record", lambda name: rec.get(name, {}))
    problems = G.score_refusals("abc", freeze_check=lambda c: c)
    assert problems == ["the score has run (run once): its record stands"]


def test_power_and_screen_refuse_a_screen_with_failed_pots(G, monkeypatch):
    monkeypatch.setattr(G, "read_record", lambda name: {"complete": True, "refused_for_failed_pot_runs": True,
                                                        "failed_pot_runs": 3} if name == "screen" else {})
    with pytest.raises(SystemExit, match="failed pot runs"):
        G.run_power()


def test_a_failed_pot_is_counted_not_raised(G, monkeypatch):
    G._init_pots()

    def broken(*a, **k):
        raise RuntimeError("stiff corner")

    monkeypatch.setattr(G.sulfur, "integrate_sulfur", broken)
    failures = []
    jobs = G.build_tasks("B0_pinned", "primary", [{}], BREAKPOINTS, G.pots_for("primary")[:2])
    out = G.run_jobs(jobs, failures)
    assert len(failures) == 2 and all("failed" in r for r in out)
    ens = G.assemble("B0_pinned", "primary", [{}], out)
    assert ens.failed == 2
