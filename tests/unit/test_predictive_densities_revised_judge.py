"""The revised decision rule (sec. 7 item 5 and M1.1, adopted 2026-11-06, PR #102; the unresolved-replicate bound,
PR #105) in src/kinetic_core/predictive_densities.py: the declared noise and its refusals, the spent-conditioned
density, the forms, every verdict branch, the bound, the prototype's verdict (scripts/diagnostics/judge_calibration.py
``_verdicts``) on shared inputs, and planted-truth worlds: the observed statistic's rank in the null (uniform), power
against a Neyman-Pearson oracle, and false credit by the exact binomial test in df-carrying and spent-row worlds.
Synthetic only; no case is scored. The statistical worlds are marked ``slow``."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from src.kinetic_core import predictive_densities as pdn
from src.kinetic_core.predictive_densities import MeasurementModelRefusal

ROOT = Path(__file__).resolve().parents[2]
ALL_INPUTS = {k: False for k in pdn.NULL_INPUTS} | {"decisive_data": True}
PREC = 100.0                     # the declared precision requirement: ESS >= max(100, PREC), i.e. the shared floor
GATE = {"misspecified_false_credit": (0, 1000)}   # a supplied, passing power check (the synthetic worlds assume one)


# --- the declared noise (item 5) -----------------------------------------------------------------------------------

def _within(unit="u1", obs="shape", sd=0.1, df=4.0):
    return pdn.UnitNoiseModel(unit, obs, "within_run", (pdn.NoiseTerm("within_run_sd", sd, obs, df=df),
                                                       pdn.NoiseTerm("instrument_floor", 0.02, obs)))


def test_a_bound_is_refused_as_an_sd():
    with pytest.raises(MeasurementModelRefusal, match="bound"):
        pdn.NoiseTerm("discrepancy", 0.4, "ph", df=3.0, is_bound=True)


def test_a_term_from_another_observable_is_not_transferred():
    thiol_omega = pdn.NoiseTerm("discrepancy", 0.419, "thiol_concentration", df=10.0)
    with pytest.raises(MeasurementModelRefusal, match="another observable"):
        pdn.UnitNoiseModel("hou", "ph", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "ph", df=1.0), thiol_omega))


def test_an_instrument_floor_is_never_the_whole_noise():
    with pytest.raises(MeasurementModelRefusal, match="floor"):
        pdn.UnitNoiseModel("u", "ph", "absolute", (pdn.NoiseTerm("instrument_floor", 0.05, "ph"),))
    with pytest.raises(MeasurementModelRefusal, match="floor"):
        pdn.UnitNoiseModel("u", "ph", "within_run", (pdn.NoiseTerm("instrument_floor", 0.05, "ph"),))
    # with a declared within-run prior component, bounded below by the floor, the within-run value is scoreable
    m = pdn.UnitNoiseModel("u", "ph", "within_run", (pdn.NoiseTerm("instrument_floor", 0.05, "ph"),
                                                     pdn.NoiseTerm("within_run_prior", 0.03, "ph", df=3.0)))
    assert m.role == "credit"
    assert all(np.all(pp >= 0.05 ** 2 - 1e-15) for pp, _ in m.variance_nodes(3))


def test_printed_sd_reading():
    """Item 5: an SD whose replicate level is unstated is repeatability (a floor); an unstated n is taken as 2."""
    assert pdn.printed_run_sd(0.1, "ph").kind == "instrument_floor"
    t = pdn.printed_run_sd(0.1, "ph", replicate_level="independent_runs")
    assert t.kind == "run_sd" and t.df == 1.0
    assert pdn.printed_run_sd(0.1, "ph", n=3, replicate_level="independent_runs").df == 2.0


def test_roles_follow_the_declaration():
    rep = pdn.NoiseTerm("reproducibility", 0.4, "ph", df=1.0)
    disc = pdn.NoiseTerm("discrepancy", 0.1, "ph", df=5.0)
    assert pdn.UnitNoiseModel("u", "ph", "absolute", (rep,)).role == "falsification_only"
    assert pdn.UnitNoiseModel("u", "ph", "absolute", (rep, disc)).role == "credit"
    assert pdn.UnitNoiseModel("u", "ph", "cross_lab", (rep, disc)).role == "not_decisive"
    tau = pdn.NoiseTerm("between_lab", 0.9, "ph", structure="run_offset", df=8.0)
    assert pdn.UnitNoiseModel("u", "ph", "cross_lab", (rep, disc, tau)).role == "credit"
    lane = pdn.NoiseTerm("lane_within_lab", 0.419, "thiol_concentration")
    assert pdn.UnitNoiseModel("u", "ph", "absolute", (rep, disc, lane)).role == "falsification_only"
    with pytest.raises(MeasurementModelRefusal, match="cancels"):
        pdn.UnitNoiseModel("u", "ph", "within_run", (pdn.NoiseTerm("reproducibility", 0.4, "ph", "run_offset", 1.0),))


def test_fixed_noise_log_lik_is_the_gaussian_density():
    a = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])
    m = pdn.UnitNoiseModel("u", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", [0.1, 0.2, 0.3], "lvl"),
                                                    pdn.NoiseTerm("reproducibility", 0.5, "lvl", "run_offset")),
                           design=a)
    y, f = np.array([0.3, -0.2]), np.array([[0.1, 0.0], [0.0, 0.4]])
    cov = a @ np.diag([0.01, 0.04, 0.09]) @ a.T + 0.25 * np.outer(a.sum(1), a.sum(1))
    want = [pdn.gaussian_log_density(y, fi, cov) for fi in f]
    assert np.allclose(m.log_lik(y, f), want)
    assert m.log_lik(np.vstack([y, y]), f).shape == (2, 2)


def test_integrated_sd_tracks_the_student_t():
    """A per-point SD with nu df integrated over its scaled inverse chi-square is a Student-t; the K = 4 quadrature
    is within 0.15 nats of it in the body."""
    m = pdn.UnitNoiseModel("u", "x", "within_run", (pdn.NoiseTerm("within_run_sd", 1.0, "x", df=6.0),))
    for r in (0.0, 1.0, 2.0):
        assert abs(m.log_lik(np.array([r]), np.zeros((1, 1)))[0] - stats.t(6).logpdf(r)) < 0.15


def test_bootstrap_draw_redraws_the_estimate():
    m = _within(sd=0.1, df=4.0)
    rng = np.random.default_rng(0)
    sds = [m.bootstrap_draw(np.zeros(3), rng)[1].terms[0].sd for _ in range(400)]
    assert np.std(sds) > 0.01 and m.terms[1].sd == 0.02
    assert m.redraws == frozenset({"replicate_sds"})


# --- the spent-conditioned density ---------------------------------------------------------------------------------

def test_conditioned_lppd_matches_the_conjugate_gaussian():
    """theta ~ N(0, 1); y_spent, y_dec ~ N(theta, 0.5^2): p(y_dec | y_spent) is Gaussian in closed form."""
    rng = np.random.default_rng(2)
    th = rng.standard_normal(20000)
    s, ys, yd = 0.5, 0.4, 0.9
    est = pdn.conditioned_lppd(stats.norm(th, s).logpdf(yd), stats.norm(th, s).logpdf(ys), precision_ess=PREC)
    post_var = 1.0 / (1.0 + 1.0 / s ** 2)
    exact = stats.norm(post_var * ys / s ** 2, math.sqrt(post_var + s ** 2)).logpdf(yd)
    assert est.accepted and abs(est.lppd - exact) < 0.02


@pytest.mark.parametrize("bad", [None, 0.0, -5.0, float("inf")])
def test_precision_ess_has_no_default_and_zero_is_refused(bad):
    ll = np.zeros(200)
    with pytest.raises(MeasurementModelRefusal, match="precision_ess"):
        pdn.conditioned_lppd(ll, ll, precision_ess=bad)
    with pytest.raises(TypeError):
        pdn.conditioned_lppd(ll, ll)                                   # required keyword: no silent default


# --- forms ---------------------------------------------------------------------------------------------------------

def test_forms_are_classified_as_the_text_says():
    f = pdn.classify_forms([pdn.FormDeclaration("r02", "experimental_fact"), pdn.FormDeclaration("r05",
                            "experimental_fact"), pdn.FormDeclaration("unforced", "modelling_choice"),
                            pdn.FormDeclaration("omega_half", "omega_hedge")])
    assert f["fact_weights"] == {"r02": 0.5, "r05": 0.5}
    assert f["scored_forms"] == ["unforced"] and f["dropped_forms"] == ["omega_half"]
    with pytest.raises(MeasurementModelRefusal, match="equal unless sourced"):
        pdn.classify_forms([pdn.FormDeclaration("x", "experimental_fact", 1.0),
                            pdn.FormDeclaration("y", "experimental_fact", 3.0)])
    g = pdn.classify_forms([pdn.FormDeclaration("x", "experimental_fact", 1.0, "paper table 2"),
                            pdn.FormDeclaration("y", "experimental_fact", 3.0, "paper table 2")])
    assert g["fact_weights"] == {"x": 0.25, "y": 0.75} and g["scored_forms"] == ["base"]


# --- the null's tail test, keeping +-inf ------------------------------------------------------------------------

def test_tail_test_keeps_infinities_that_weighted_quantile_drops():
    x = np.r_[np.random.default_rng(0).normal(0, 1, 950), np.full(50, np.inf)]
    assert math.isfinite(pdn.weighted_quantile(x, None, 0.95)[0])       # the prototype's estimator drops them
    r = pdn.null_tail_test(x, None, 1e9, 0.05, "upper")
    assert r["threshold"] == np.inf and not r["crosses"] and abs(r["infinite_mass"] - 0.05) < 1e-12
    r = pdn.null_tail_test(-x, None, -1e9, 0.05, "lower")               # the mirror
    assert r["threshold"] == -np.inf and not r["crosses"]
    r = pdn.null_tail_test(np.r_[x[:950], 0.0, np.full(49, np.inf)], None, 1e9, 0.05, "upper")
    assert r["crosses"] and math.isfinite(r["threshold"])
    with pytest.raises(MeasurementModelRefusal, match="NaN"):
        pdn.null_tail_test(np.r_[x[:10], np.nan], None, 0.0, 0.05, "upper")


def test_tail_test_crossing_is_the_strict_tail_and_matches_its_threshold():
    rng = np.random.default_rng(3)
    x, w = rng.normal(0, 1, 2000), rng.uniform(0.2, 1.0, 2000)
    for t in rng.normal(1.6, 0.4, 200):
        r = pdn.null_tail_test(x, w, t, 0.05, "upper")
        tail = float(np.sum(w[x >= t]) / w.sum())
        assert r["crosses"] == (tail < 0.05) == (t > r["threshold"])
        assert abs(r["tail_weight"] - tail) < 1e-12 and r["tail_weight_mc_se"] > 0


# --- the verdict's branches (densities and nulls given directly) ---------------------------------------------------

def _obs(T, p_a=None, p_b=None, accepted=True, units=("u1",)):
    """A FormObservation whose T keys are ``T`` (rival b's densities 0)."""
    def de(v, ok=True):
        return pdn.DensityEstimate(float(v), "exact", ok, "" if ok else "k-hat 0.9 >= 0.7")
    la = {k: de(v, accepted if k == "full" else True) for k, v in T.items()}
    lb = {k: de(0.0) for k in T}
    pa = p_a or {u: 0.5 for u in units}
    pb = p_b or {u: 0.5 for u in units}
    return pdn.FormObservation({"a": la, "b": lb}, {"a": pa, "b": pb})


def _null(truth, stats_, n=1000, redrawn=frozenset({"decisive_data"}), outside=True, resolved=None, weights=None,
          extra=None):
    return pdn.NullDistribution(truth, {"base": {k: np.asarray(v) for k, v in stats_.items()}}, frozenset(redrawn),
                                outside, np.ones(n, bool) if resolved is None else resolved, weights, extra or {})


ONE_UNIT = {"u1": _within("u1")}
RNG = np.random.default_rng(7)
NULL_A = _null("a", {"full": RNG.normal(5.0, 1.0, 1000)})
NULL_B = _null("b", {"full": RNG.normal(-5.0, 1.0, 1000)})


def test_a_credit_and_the_mirror_credit():
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "credited" and v.credited == "a"
    d = v.as_dict()["per_form"]["base"]
    assert d["a"]["threshold"] < -3 and d["a"]["threshold_mc_se"] > 0 and d["b"]["threshold_quantile"] == 0.05
    assert v.details["null_a"]["unresolved_fraction"] == 0.0 and not v.details["bound_differs"]
    v = pdn.judge_pair({"base": _obs({"full": -8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "credited" and v.credited == "b"


def test_both_qualify_credits_neither():
    """T = 0 lies above the 95th percentile under b AND below the 5th under a: both meet their conditions."""
    v = pdn.judge_pair({"base": _obs({"full": 0.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "neither_credited_both_qualify" and v.credited is None and v.carried == ("a", "b")


def test_the_adequacy_gate_after_holm_and_its_mc_error():
    p_bad = {"u1": 0.001, "u2": 0.5, "lvl": 0.5}           # Holm-adjusted 0.003 < 0.05
    p_ok = {"u1": 0.02, "u2": 0.5, "lvl": 0.5}             # raw 0.02, Holm-adjusted 0.06 >= 0.05: adequate
    units = {"u1": _within("u1"), "u2": _within("u2"), "lvl": pdn.UnitNoiseModel(
        "lvl", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "lvl", df=1.0),))}
    T = {"full": 9.0, "-u1": 4.5, "-u2": 4.5, "u1": 4.5, "u2": 4.5}
    st = {k: RNG.normal(-5.0, 1.0, 1000) for k in T}
    nb, na = _null("b", st), _null("a", {k: RNG.normal(5.0, 1.0, 1000) for k in T})
    v = pdn.judge_pair({"base": _obs(T, p_a=p_bad, units=units)}, na, nb, units, ALL_INPUTS, **GATE)
    assert v.outcome == "undecided" and not v.details["per_form"]["base"]["a"]["adequate"]
    assert v.details["unit_roles"]["falsification_only"] == ["lvl"]
    ob = _obs(T, p_a=p_ok, units=units)
    ob = pdn.FormObservation(ob.lppd, ob.adequacy_p, {"a": {u: 0.007 for u in units}, "b": {u: 0.02 for u in units}})
    v = pdn.judge_pair({"base": ob}, na, nb, units, ALL_INPUTS, **GATE)
    assert v.outcome == "credited" and v.credited == "a"
    assert v.details["per_form"]["base"]["a"]["adequacy"]["u1"]["p_mc_se"] == 0.007


def test_a_credit_resting_on_one_unit_needs_that_unit_alone():
    units = {"u1": _within("u1"), "u2": _within("u2")}
    T = {"full": 3.0, "-u1": -2.0, "-u2": 3.0, "u1": 0.5, "u2": 0.1}
    nb = _null("b", {"full": RNG.normal(0, 1, 1000), "-u1": RNG.normal(0, 1, 1000), "-u2": RNG.normal(0, 1, 1000),
                     "u1": RNG.normal(0, 1, 1000), "u2": RNG.normal(0, 1, 1000)})
    na = _null("a", {k: RNG.normal(-10, 1, 1000) for k in T})          # b never crosses its mirror threshold
    v = pdn.judge_pair({"base": _obs(T, units=units)}, na, nb, units, ALL_INPUTS, **GATE)
    assert v.details["per_form"]["base"]["a"]["crosses"] and v.details["per_form"]["base"]["a"]["adequate"]
    loo = v.details["per_form"]["base"]["a"]["leave_one_unit_out"]
    assert loo["u1"]["rests_on"] and not loo["u1"]["alone_credits"] and not loo["u2"]["rests_on"]
    assert v.outcome == "undecided"


def test_sole_survivor_is_carried_not_credited():
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS,
                       falsified={"b": True}, incumbent="b", **GATE)
    assert v.outcome == "carried_sole_survivor" and v.credited is None and v.carried == ("a",)
    assert v.details["s3_refusal_applies"]
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, None, None, ONE_UNIT, ALL_INPUTS,
                       falsified={"a": True, "b": True})
    assert v.outcome == "both_falsified" and v.carried == ()


def test_falsification_is_read_before_not_decisive():
    """D1 (review of PR #105): with no credit unit, a falsified rival was returned as 'carried' by 'not_decisive'."""
    lvl = {"lvl": pdn.UnitNoiseModel("lvl", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "lvl", df=1.),))}
    ob = {"base": _obs({"full": 8.0}, units=("lvl",))}
    v = pdn.judge_pair(ob, NULL_A, NULL_B, lvl, ALL_INPUTS, falsified={"b": True})
    assert v.outcome == "carried_sole_survivor" and v.carried == ("a",) and "b" not in v.carried
    v = pdn.judge_pair(ob, NULL_A, NULL_B, lvl, ALL_INPUTS, falsified={"a": True, "b": True})
    assert v.outcome == "both_falsified" and v.carried == ()
    v = pdn.judge_pair(ob, NULL_A, NULL_B, lvl, ALL_INPUTS)
    assert v.outcome == "not_decisive" and v.carried == ("a", "b")


def test_an_unresolved_density_is_not_computable():
    v = pdn.judge_pair({"base": _obs({"full": 8.0}, accepted=False)}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "not_computable" and v.credited is None and v.carried == ("a", "b")
    assert any("unresolved density" in r for r in v.reasons)


@pytest.mark.parametrize("kind", ["forced_inputs", "few_replicates", "inside_ensemble"])
def test_an_incomplete_null_is_not_computable(kind):
    inputs = dict(ALL_INPUTS)
    nb = NULL_B
    if kind == "forced_inputs":
        inputs["forced_inputs"] = True                      # the case has forced inputs; the null did not redraw them
    elif kind == "few_replicates":
        nb = _null("b", {"full": RNG.normal(-5, 1, 999)}, n=999)
    else:
        nb = _null("b", {"full": RNG.normal(-5, 1, 1000)}, outside=False)
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, nb, ONE_UNIT, inputs, **GATE)
    assert v.outcome == "not_computable", v.reasons


# --- the bound on unresolved replicates (PR #105) ----------------------------------------------------------------

def _grid(n, loc=0.0):
    """n replicates at the normal's (i + 1/2) / n quantiles: a null with no sampling noise."""
    return loc + stats.norm.ppf((np.arange(n) + 0.5) / n)


def _with_unresolved(truth, finite, n_bad, fill=np.nan, **kw):
    x = np.r_[finite, np.full(n_bad, fill)]
    res = np.r_[np.ones(finite.size, bool), np.zeros(n_bad, bool)]
    return _null(truth, {"full": x}, n=x.size, resolved=res, **kw)


def test_unresolved_replicates_are_bounded_not_refused():
    """30 of 1000 replicates unresolved (stored as NaN): the pair is computable; T = 2.0 would credit a with them
    dropped (q95 of N(0, 1) ~ 1.64) but not once they sit at +inf (the threshold moves to the 98th percentile of the
    finite part, ~2.06); T = 3.0 survives the bound."""
    nb = _with_unresolved("b", _grid(970), 30)
    na = _null("a", {"full": _grid(1000, -10.0)})                     # b never qualifies (T far above its 5 %)
    v = pdn.judge_pair({"base": _obs({"full": 2.0})}, na, nb, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "undecided" and any(r.startswith("a would be credited") for r in v.reasons), v.reasons
    assert abs(v.details["null_b"]["unresolved_fraction"] - 0.03) < 1e-12 and v.details["null_b"]["unresolved_count"] == 30
    assert v.details["bound_differs"]
    v = pdn.judge_pair({"base": _obs({"full": 3.0})}, na, nb, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome == "credited" and v.credited == "a"


@pytest.mark.parametrize("n_bad, credits", [(49, True), (50, False), (80, False)])
def test_an_unresolved_fraction_reaching_alpha_cannot_credit(n_bad, credits):
    far = 1e6
    nb = _with_unresolved("b", _grid(1000 - n_bad), n_bad)
    na = _with_unresolved("a", _grid(1000 - n_bad), n_bad)
    v = pdn.judge_pair({"base": _obs({"full": far})}, NULL_A, nb, ONE_UNIT, ALL_INPUTS, **GATE)
    assert (v.credited == "a") == credits, v.reasons
    v = pdn.judge_pair({"base": _obs({"full": -far})}, na, NULL_B, ONE_UNIT, ALL_INPUTS, **GATE)
    assert (v.credited == "b") == credits, v.reasons


def test_the_unresolved_fraction_is_weighted():
    """10 unresolved of 2000, but carrying weight 12 each: 120 / 2110 = 0.057 >= 0.05, so nothing credits."""
    w = np.r_[np.ones(1990), np.full(10, 12.0)]
    nb = _with_unresolved("b", _grid(1990), 10, weights=w, extra={"weights_k_hat": 0.1})
    assert nb.n_effective >= 1000
    v = pdn.judge_pair({"base": _obs({"full": 1e6})}, NULL_A, nb, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.outcome != "credited" and abs(v.details["null_b"]["unresolved_fraction"] - 120 / 2110) < 1e-12


def test_the_bound_applies_in_both_nulls():
    """Testing a's credit, null a's unresolved replicates sit at +inf too, which makes b qualify more easily (and both
    qualifying blocks a). Grid null a: 960 finite + 40 unresolved; 49 finite replicates lie at or below T = -1.635.
    Dropped, that is 49 / 960 = 5.1 %: b does not cross, and a (far above null b) would be credited. At +inf it is
    49 / 1000 = 4.9 %: b qualifies too, and a is not credited."""
    na_drop = _null("a", {"full": _grid(960)}, n=960)
    na = _with_unresolved("a", _grid(960), 40)
    nb = _null("b", {"full": _grid(1000, -10.0)})
    obs = {"base": _obs({"full": -1.635})}
    v = pdn.judge_pair(obs, na_drop, nb, ONE_UNIT, ALL_INPUTS, min_replicates=900, **GATE)
    assert v.outcome == "credited" and v.credited == "a"
    v = pdn.judge_pair(obs, na, nb, ONE_UNIT, ALL_INPUTS, **GATE)
    assert v.credited is None and v.details["qualifies_under_bound"]["testing_a_unresolved_plus_inf"]["b"]
    assert any(r.startswith("a would be credited") for r in v.reasons), v.reasons


def test_case_inputs_must_be_stated_in_full():
    with pytest.raises(MeasurementModelRefusal, match="every listed null input"):
        pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, {"decisive_data": True})


def test_no_credit_unit_is_not_decisive():
    lvl = {"lvl": pdn.UnitNoiseModel("lvl", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "lvl", df=1.),))}
    v = pdn.judge_pair({"base": _obs({"full": 8.0}, units=("lvl",))}, NULL_A, NULL_B, lvl, ALL_INPUTS)
    assert v.outcome == "not_decisive"


def test_the_power_check_gate_is_required_for_a_credit():
    """D2 (review of PR #105): a credit needs the power check supplied AND passing."""
    obs = {"base": _obs({"full": 8.0})}
    v = pdn.judge_pair(obs, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "undecided" and v.credited is None and "not supplied" in v.reasons[0]
    assert v.details["power_check_gate"] == "not supplied" and v.details["qualifies"]["a"]
    v = pdn.judge_pair(obs, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, misspecified_false_credit=(25, 200))
    assert v.outcome == "undecided" and "Wilson lower end" in v.reasons[0]        # Wilson lower 0.086 > 0.05
    v = pdn.judge_pair(obs, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS, misspecified_false_credit=(12, 200))
    assert v.outcome == "credited"
    v = pdn.judge_pair({"base": _obs({"full": 0.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "neither_credited_both_qualify"                          # no credit: the gate is not read


def test_required_agreement_needs_every_scored_form():
    forms = [pdn.FormDeclaration("f1", "modelling_choice"), pdn.FormDeclaration("f2", "modelling_choice"),
             pdn.FormDeclaration("hedge", "omega_hedge")]
    na = pdn.NullDistribution("a", {f: {"full": RNG.normal(5, 1, 1000)} for f in ("f1", "f2")},
                              frozenset({"decisive_data"}), True)
    nb = pdn.NullDistribution("b", {f: {"full": RNG.normal(-5, 1, 1000)} for f in ("f1", "f2")},
                              frozenset({"decisive_data"}), True)
    v = pdn.judge_pair({"f1": _obs({"full": 8.0}), "f2": _obs({"full": -4.5})}, na, nb, ONE_UNIT, ALL_INPUTS, forms,
                       **GATE)
    assert v.outcome == "undecided" and not v.details["per_form"]["f2"]["a"]["crosses"]
    assert v.details["per_form"]["f1"]["a"]["qualifies"] and not v.details["qualifies"]["a"]
    v = pdn.judge_pair({"f1": _obs({"full": 8.0}), "f2": _obs({"full": 7.0})}, na, nb, ONE_UNIT, ALL_INPUTS, forms,
                       **GATE)
    assert v.outcome == "credited" and v.details["forms"]["dropped_forms"] == ["hedge"]
    with pytest.raises(MeasurementModelRefusal, match="scored forms"):
        pdn.judge_pair({"f1": _obs({"full": 8.0})}, na, nb, ONE_UNIT, ALL_INPUTS, forms)


# --- the prototype's verdict on shared inputs ----------------------------------------------------------------------

def _prototype():
    spec = importlib.util.spec_from_file_location("judge_calibration_proto",
                                                  ROOT / "scripts" / "diagnostics" / "judge_calibration.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_reproduces_the_prototype_verdict_where_the_text_agrees():
    """``_verdicts`` (the --revised prototype) and ``judge_pair`` on the same T, nulls and adequacy: identical
    credits on every draw except where both rivals qualify (the prototype credits a; the text credits neither). The
    crossing here is strict in the null's tail (the prototype's is not): they differ only when T falls between two
    adjacent order statistics at the 5 % boundary, which these draws never hit."""
    proto = _prototype()
    rng = np.random.default_rng(11)
    units = {"u1": _within("u1"), "u2": _within("u2")}
    keys = ["full", "-u1", "-u2", "u1", "u2"]
    n_cmp = n_both = 0
    for _ in range(150):
        T = {k: float(rng.normal(0, 3)) for k in keys}
        sa = {k: rng.normal(rng.normal(2, 2), 1, 1000) for k in keys}
        sb = {k: rng.normal(rng.normal(-2, 2), 1, 1000) for k in keys}
        ad = {"a": bool(rng.random() < 0.8), "b": bool(rng.random() < 0.8)}
        pa = {"u1": 0.5 if ad["a"] else 0.001, "u2": 0.5}
        pb = {"u1": 0.5 if ad["b"] else 0.001, "u2": 0.5}
        mine = pdn.judge_pair({"base": _obs(T, pa, pb, units=units)}, _null("a", sa), _null("b", sb), units,
                              ALL_INPUTS, **GATE)
        theirs = proto._verdicts({"base": T}, lambda f, k, q: pdn.weighted_quantile(sa[k], None, q),
                                 lambda f, k, q: pdn.weighted_quantile(sb[k], None, q), ["u1", "u2"], ["base"],
                                 lambda f: ad["a"], lambda f: ad["b"])
        if mine.outcome == "neither_credited_both_qualify":
            n_both += 1
            assert theirs["credit"] == "a"
            continue
        n_cmp += 1
        assert (mine.credited or "none") == theirs["credit"], (T, ad, mine.reasons)
    assert n_cmp > 100


# --- worlds from draws -------------------------------------------------------------------------------------------

# The world: two within-run shapes (credit) and one absolute level (falsification only: no discrepancy estimate);
# rival a predicts theta x, rival b theta x^1.35, theta ~ N(1, PRIOR_SD); per-point noise 0.2 at a fixed same-observable
# within-run repeat (``df_shape`` makes it a printed within-run SD with that many df instead). Synthetic only.
X1, X2 = np.array([0.2, 0.4, 0.6, 0.8, 1.0]), np.array([1.2, 1.4, 1.6, 1.8])
P_EXP, SD, PRIOR_SD, S_DRAWS, POOL = 1.35, 0.2, 0.015, 1500, 1000


def _world_noise(sd=SD, floor=0.02, df_shape=None):
    shape_term = (pdn.NoiseTerm("reproducibility", sd, "shape") if df_shape is None
                  else pdn.NoiseTerm("within_run_sd", sd, "shape", df=float(df_shape)))
    m = {u: pdn.UnitNoiseModel(u, "shape", "within_run", (shape_term, pdn.NoiseTerm("instrument_floor", floor, "shape")))
         for u in ("shape1", "shape2")}
    m["level"] = pdn.UnitNoiseModel("level", "level", "absolute",
                                    (pdn.NoiseTerm("reproducibility", 0.1, "level"),
                                     pdn.NoiseTerm("reproducibility", 0.15, "level", "run_offset", 1.0)))
    return m


def _world_draws(rival, n, rng, p_exp=P_EXP, prior_sd=PRIOR_SD):
    th = rng.normal(1.0, prior_sd, n)
    p = 1.0 if rival == "a" else p_exp
    return pdn.RivalDraws({"shape1": th[:, None] * X1 ** p, "shape2": th[:, None] * X2 ** p,
                           "level": np.column_stack([th, th])})


def _build(design, draws, pool=POOL, seed=0, **pool_kw):
    """Scoring ensembles, held-out truths and both pools, seeded as the PR #105 review's experiment (its run1.jsonl)."""
    rng = np.random.default_rng(1 + 100 * seed)
    ea, eb = {"base": {"base": draws("a", S_DRAWS, rng)}}, {"base": {"base": draws("b", S_DRAWS, rng)}}
    ta = draws("a", S_DRAWS, np.random.default_rng(2 + 100 * seed))
    tb = draws("b", S_DRAWS, np.random.default_rng(3 + 100 * seed))
    pa = pdn.build_null_pool("a", {"base": ta}, design, ea, eb, pool, np.random.default_rng(4 + 100 * seed),
                             precision_ess=PREC, **pool_kw)
    pb = pdn.build_null_pool("b", {"base": tb}, design, ea, eb, pool, np.random.default_rng(5 + 100 * seed),
                             precision_ess=PREC, **pool_kw)
    return {"design": design, "ea": ea, "eb": eb, "ta": ta, "tb": tb, "pa": pa, "pb": pb}


def _tally(w, n_sets, seed=0, ppc_reps=100, spent=False):
    """Verdicts on ``n_sets`` synthetic sets per truth (truth draw, then data at it), the review's seeds."""
    design = w["design"]
    out = {}
    for truth, td in (("a", w["ta"]), ("b", w["tb"])):
        r = np.random.default_rng((10 if truth == "a" else 11) + 100 * seed)
        tally: dict = {}
        for _ in range(n_sets):
            j = int(r.integers(td.n))
            y = {u: design.noise[u].sample(td.pred(u)[j][None, :], r)[0] for u in design.noise}
            v = pdn.judge_from_draws(y, design, w["ea"], w["eb"], w["pa"], w["pb"], r, ppc_reps=ppc_reps,
                                     precision_ess=PREC, **GATE)
            key = v.credited or v.outcome
            tally[key] = tally.get(key, 0) + 1
        out[truth] = tally
    return out


@pytest.fixture(scope="module")
def world():
    design = pdn.MeasurementDesign(_world_noise())
    rng = np.random.default_rng(1)
    ea, eb = {"base": {"base": _world_draws("a", S_DRAWS, rng)}}, {"base": {"base": _world_draws("b", S_DRAWS, rng)}}
    return {"design": design, "ea": ea, "eb": eb}


def test_unresolved_draw_densities_are_not_computable(world):
    """At a noise far below the draws' spread the joint ESS collapses: the verdict refuses rather than reads it."""
    design = pdn.MeasurementDesign(_world_noise(sd=0.0005, floor=0.0001))
    rng = np.random.default_rng(21)
    y = {u: world["ea"]["base"]["base"].pred(u)[3] + 0.001 for u in design.noise}
    obs = pdn.score_from_draws(y, design, world["ea"], world["eb"], rng, ppc_reps=50, precision_ess=PREC)
    assert not obs["base"].lppd["a"]["full"].accepted
    assert all(0 <= se < 0.08 for se in obs["base"].adequacy_mc_se["a"].values())        # the adequacy p's MC error
    pa = pdn.build_null_pool("a", {"base": _world_draws("a", 200, np.random.default_rng(26))}, design, world["ea"],
                             world["eb"], 5, rng, precision_ess=PREC)
    pb = pdn.build_null_pool("b", {"base": _world_draws("b", 200, np.random.default_rng(27))}, design, world["ea"],
                             world["eb"], 5, rng, precision_ess=PREC)
    v = pdn.judge_pair(obs, pdn.null_from_pool(pa), pdn.null_from_pool(pb),
                       {u: design.noise[u] for u in design.decisive_units()}, design.case_inputs(), **GATE)
    assert v.outcome == "not_computable" and any("unresolved density" in r for r in v.reasons)


def test_a_pool_without_forced_inputs_is_incomplete(world):
    design = pdn.MeasurementDesign(_world_noise(), forced_inputs_present=True)
    rng = np.random.default_rng(22)
    ta = _world_draws("a", 500, np.random.default_rng(23))
    tb = _world_draws("b", 500, np.random.default_rng(24))
    pa = pdn.build_null_pool("a", {"base": ta}, design, world["ea"], world["eb"], 20, rng, precision_ess=PREC)
    pb = pdn.build_null_pool("b", {"base": tb}, design, world["ea"], world["eb"], 20, rng, precision_ess=PREC)
    assert "forced_inputs" not in pa.redrawn
    y = {u: ta.pred(u)[0] for u in design.noise}
    v = pdn.judge_from_draws(y, design, world["ea"], world["eb"], pa, pb, rng, ppc_reps=50, precision_ess=PREC, **GATE)
    assert v.outcome == "not_computable"
    assert any("forced_inputs" in r for r in v.reasons) and any("effective replicates" in r for r in v.reasons)


def test_truths_inside_the_scoring_ensemble_are_flagged(world):
    design, ens = world["design"], world["ea"]["base"]["base"]
    p = pdn.build_null_pool("a", {"base": ens}, design, world["ea"], world["eb"], 3, np.random.default_rng(25),
                            precision_ess=PREC)
    assert not p.truths_outside_scoring_ensemble
    # copied row by row (new buffers, a different number of draws): still the scoring ensemble's draws
    rows = np.random.default_rng(28).choice(ens.n, 40, replace=False)
    copied = pdn.RivalDraws({u: np.array([np.array(ens.pred(u)[i]) for i in rows]) for u in ens.predictions})
    p = pdn.build_null_pool("a", {"base": copied}, design, world["ea"], world["eb"], 3, np.random.default_rng(25),
                            precision_ess=PREC)
    assert not p.truths_outside_scoring_ensemble
    # one copied row among fresh draws is enough
    fresh = _world_draws("a", 40, np.random.default_rng(29))
    mixed = pdn.RivalDraws({u: np.vstack([fresh.pred(u)[:-1], ens.pred(u)[7:8]]) for u in ens.predictions})
    assert not pdn.build_null_pool("a", {"base": mixed}, design, world["ea"], world["eb"], 3,
                                   np.random.default_rng(25), precision_ess=PREC).truths_outside_scoring_ensemble
    assert pdn.build_null_pool("a", {"base": fresh}, design, world["ea"], world["eb"], 3, np.random.default_rng(25),
                               precision_ess=PREC).truths_outside_scoring_ensemble


def test_spent_rows_condition_the_null():
    """With a spent unit the null is the pool re-weighted to the truth rival's spent-conditioned posterior."""
    noise = _world_noise()
    noise["spent"] = pdn.UnitNoiseModel("spent", "level", "absolute",
                                        (pdn.NoiseTerm("reproducibility", 0.3, "level", df=3.0),))
    design = pdn.MeasurementDesign(noise, spent_units=("spent",))
    assert design.case_inputs()["conditioning_rows"]

    def dr(rival, n, rng):
        d = _world_draws(rival, n, rng)
        return pdn.RivalDraws({**d.predictions, "spent": d.predictions["level"][:, :1]})
    rng = np.random.default_rng(31)
    ea, eb = {"base": {"base": dr("a", 600, rng)}}, {"base": {"base": dr("b", 600, rng)}}
    pool = pdn.build_null_pool("a", {"base": dr("a", 600, np.random.default_rng(32))}, design, ea, eb, 40, rng,
                               precision_ess=PREC)
    assert "conditioning_rows" in pool.redrawn
    nl = pdn.null_from_pool(pool, {"spent": np.array([1.1])})
    assert nl.weights is not None and 1.0 < nl.n_effective < 40 and "weights_k_hat" in nl.extra


def test_an_unresolved_replicate_is_re_estimated_once_with_the_retry_draws():
    """Few scoring draws and a printed SD with 2 df: some replicates' re-drawn SDs are sharp enough to leave ESS < 100.
    Each is re-scored once with an independent, larger draw set; most are rescued; the rest stay unresolved."""
    design = pdn.MeasurementDesign(_world_noise(df_shape=2.0))
    rng = np.random.default_rng(41)
    ea, eb = {"base": {"base": _world_draws("a", 150, rng)}}, {"base": {"base": _world_draws("b", 150, rng)}}
    big = ({"base": {"base": _world_draws("a", 3000, rng)}}, {"base": {"base": _world_draws("b", 3000, rng)}})
    truths = {"base": _world_draws("a", 500, np.random.default_rng(42))}
    plain = pdn.build_null_pool("a", truths, design, ea, eb, 80, np.random.default_rng(43), precision_ess=PREC)
    retried = pdn.build_null_pool("a", truths, design, ea, eb, 80, np.random.default_rng(43), precision_ess=PREC,
                                  retry_draws=big)
    assert plain.extra["n_retried"] == 0 and plain.extra["n_unresolved_after_retry"] > 0
    assert retried.extra["n_retried"] == plain.extra["n_unresolved_after_retry"]           # the same replicates
    assert retried.extra["n_rescued_by_retry"] > 0
    assert retried.extra["n_unresolved_after_retry"] == retried.extra["n_retried"] - retried.extra["n_rescued_by_retry"]
    assert pdn.null_from_pool(retried).extra["n_rescued_by_retry"] == retried.extra["n_rescued_by_retry"]
    inside = pdn.build_null_pool("a", {"base": big[0]["base"]["base"]}, design, ea, eb, 2, rng, precision_ess=PREC,
                                 retry_draws=big)
    assert not inside.truths_outside_scoring_ensemble                    # the retry ensembles are checked too


# --- statistical worlds (slow) -------------------------------------------------------------------------------------

def _binom_upper_p(k, n, rate=pdn.NOMINAL_ALPHA):
    """Exact one-sided binomial p of k or more false credits in n sets at the stated rate."""
    return float(stats.binomtest(int(k), int(n), rate, alternative="greater").pvalue)


@pytest.fixture(scope="module")
def planted():
    """The fixed-noise world (no df-carrying term enters T): pools of 1000 replicates under each truth."""
    return _build(pdn.MeasurementDesign(_world_noise()), _world_draws)


R_RANK = 400


@pytest.mark.slow
def test_observed_rank_is_uniform_in_the_null_rival_bootstrap(planted):
    """(a) Under rival r as truth (fresh truth draws from r's prior, data at them), the observed T's randomised rank
    among null r's replicates, (#{T_null < T} + U) / (n + 1), is Uniform(0, 1): KS at R = 400 per truth (1 % level;
    the noise entering T is fixed, so the null's generator is the data's exactly)."""
    design = planted["design"]
    for truth, pool, seed in (("a", planted["pa"], 51), ("b", planted["pb"], 52)):
        nl = pdn.null_from_pool(pool)
        assert nl.unresolved_fraction == 0.0
        z = nl.bounded("base", "full", 1.0)
        rng = np.random.default_rng(seed)
        u = []
        for _ in range(R_RANK):
            td = _world_draws(truth, 1, rng)
            y = {k: design.noise[k].sample(td.pred(k), rng)[0] for k in design.noise}
            ob = pdn.score_from_draws(y, design, planted["ea"], planted["eb"], rng, ppc_reps=0, precision_ess=PREC)
            assert all(e.accepted for e in ob["base"].lppd["a"].values())
            t = ob["base"].T()["full"]
            u.append((np.sum(z < t) + rng.random() * (np.sum(z == t) + 1)) / (z.size + 1))
        p = stats.kstest(u, "uniform").pvalue
        print(f"[rank uniformity] truth {truth}: KS p = {p:.3f} (R = {R_RANK})")
        assert p >= 0.01, (truth, p)


R_POWER, ORACLE_LOSS = 200, 0.10


@pytest.mark.slow
def test_power_is_near_the_neyman_pearson_oracle_in_a_near_simple_world():
    """(b) theta ~ N(1, 0.002): each rival is near a simple hypothesis, so the most powerful 5 % test of b against a
    is the likelihood ratio, with power Phi(Delta - z_0.95), Delta^2 = sum (x - x^1.35)^2 / 0.2^2 (here 0.882). The
    rule may lose at most ORACLE_LOSS = 0.10 to it (the adequacy family's Holm false alarm, <= 0.05, the strict tail
    and the thresholds' MC error): the exact one-sided binomial test of power >= oracle - 0.10 is not rejected at 5 %
    (R = 200 per truth: a shortfall of 0.10 more is detected with probability ~0.98)."""
    x = np.r_[X1, X2]
    delta = math.sqrt(np.sum((x - x ** P_EXP) ** 2) / SD ** 2)
    oracle = float(stats.norm.cdf(delta - stats.norm.ppf(0.95)))

    def draws(r, n, rng):
        return _world_draws(r, n, rng, prior_sd=0.002)
    w = _build(pdn.MeasurementDesign(_world_noise()), draws, seed=3)
    rates = _tally(w, R_POWER, seed=3)
    for truth, wrong in (("a", "b"), ("b", "a")):
        k = rates[truth].get(truth, 0)
        p = float(stats.binomtest(k, R_POWER, oracle - ORACLE_LOSS, alternative="less").pvalue)
        print(f"[power vs oracle] truth {truth}: power {k / R_POWER:.3f}, oracle {oracle:.3f}, binomial p {p:.3f}; "
              f"false credit {rates[truth].get(wrong, 0)}/{R_POWER}; {rates[truth]}")
        assert p >= 0.05, rates
        assert _binom_upper_p(rates[truth].get(wrong, 0), R_POWER) >= 0.05, rates


R_FALSE = 300


def _same_mean_world(seed):
    """(c)'s world: rivals a (theta ~ N(1, 0.03)) and b (theta ~ N(1, 0.12)) both predict theta x: the same mean,
    different spread. The shapes carry a printed within-run SD with 4 df; a spent level row (theta, a run SD with 3 df)
    conditions every density and the null (re-weighted pools of 1600)."""
    noise = _world_noise(df_shape=4.0)
    noise["spent"] = pdn.UnitNoiseModel("spent", "level", "absolute",
                                        (pdn.NoiseTerm("run_sd", 0.2, "level", df=3.0),))
    design = pdn.MeasurementDesign(noise, spent_units=("spent",))

    def draws(r, n, rng):
        th = rng.normal(1.0, 0.03 if r == "a" else 0.12, n)
        return pdn.RivalDraws({"shape1": th[:, None] * X1, "shape2": th[:, None] * X2,
                               "level": np.column_stack([th, th]), "spent": th[:, None]})
    return _build(design, draws, pool=1600, seed=seed)


@pytest.mark.slow
def test_same_mean_different_spread_with_a_spent_row_does_not_over_credit():
    """(c) The exact one-sided binomial p of the false credits is >= 0.05 against a false-credit rate of 0.05, under
    each truth (R = 300 per truth)."""
    w = _same_mean_world(seed=6)
    rates = _tally(w, R_FALSE, seed=6)
    for truth, wrong in (("a", "b"), ("b", "a")):
        k, nc = rates[truth].get(wrong, 0), rates[truth].get("not_computable", 0)
        p = _binom_upper_p(k, R_FALSE)
        print(f"[same mean, different spread] truth {truth}: false credit {k}/{R_FALSE} (binomial p {p:.3f}), "
              f"not computable {nc}; {rates[truth]}")
        assert p >= 0.05, rates
        assert nc <= 0.1 * R_FALSE, rates                        # the check is not emptied by refusals


@pytest.mark.slow
@pytest.mark.parametrize("nu", [2, 4, 9])
def test_the_bound_keeps_false_credit_within_alpha_in_df_carrying_worlds(nu):
    """(d) The review's df worlds (its seeds; a printed within-run SD with nu df on both shapes): the pools hold
    unresolved replicates (4 %, 2 %, 1 % at nu = 2, 4, 9), which before PR #105's bound made every pair "not
    computable". Now every pair is computable, and the exact one-sided binomial p of the false credits is >= 0.05
    against 0.05 (R = 300 per truth)."""
    w = _build(pdn.MeasurementDesign(_world_noise(df_shape=nu)), _world_draws)
    assert (~w["pa"].resolved).any() or (~w["pb"].resolved).any()
    rates = _tally(w, R_FALSE)
    for truth, wrong in (("a", "b"), ("b", "a")):
        k = rates[truth].get(wrong, 0)
        p = _binom_upper_p(k, R_FALSE)
        print(f"[df {nu}] truth {truth}: false credit {k}/{R_FALSE} (binomial p {p:.3f}), power "
              f"{rates[truth].get(truth, 0) / R_FALSE:.3f}, pool unresolved {np.mean(~w['p' + truth].resolved):.3f}; "
              f"{rates[truth]}")
        assert p >= 0.05, rates
        assert rates[truth].get("not_computable", 0) == 0, rates
