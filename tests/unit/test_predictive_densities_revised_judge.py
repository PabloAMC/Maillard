"""The revised decision rule (sec. 7 item 5 and M1.1, adopted 2026-11-06, PR #102) in
src/kinetic_core/predictive_densities.py: the declared noise and its refusals, the spent-conditioned density, the forms,
every verdict branch, the prototype's verdict (scripts/diagnostics/judge_calibration.py ``_verdicts``) on shared
inputs, and planted-truth worlds for power and false credit. Synthetic only; no case is scored."""
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
    est = pdn.conditioned_lppd(stats.norm(th, s).logpdf(yd), stats.norm(th, s).logpdf(ys))
    post_var = 1.0 / (1.0 + 1.0 / s ** 2)
    exact = stats.norm(post_var * ys / s ** 2, math.sqrt(post_var + s ** 2)).logpdf(yd)
    assert est.accepted and abs(est.lppd - exact) < 0.02


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


def _null(truth, stats_, n=1000, redrawn=frozenset({"decisive_data"}), outside=True):
    return pdn.NullDistribution(truth, {"base": {k: np.asarray(v) for k, v in stats_.items()}}, frozenset(redrawn),
                                outside, np.ones(n, bool))


ONE_UNIT = {"u1": _within("u1")}
RNG = np.random.default_rng(7)
NULL_A = _null("a", {"full": RNG.normal(5.0, 1.0, 1000)})
NULL_B = _null("b", {"full": RNG.normal(-5.0, 1.0, 1000)})


def test_a_credit_and_the_mirror_credit():
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "credited" and v.credited == "a"
    d = v.as_dict()["per_form"]["base"]
    assert d["a"]["threshold"] < -3 and d["a"]["threshold_mc_se"] > 0 and d["b"]["threshold_quantile"] == 0.05
    v = pdn.judge_pair({"base": _obs({"full": -8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "credited" and v.credited == "b"


def test_both_qualify_credits_neither():
    """T = 0 lies above the 95th percentile under b AND below the 5th under a: both meet their conditions."""
    v = pdn.judge_pair({"base": _obs({"full": 0.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "neither_credited_both_qualify" and v.credited is None and v.carried == ("a", "b")


def test_the_adequacy_gate_after_holm():
    p_bad = {"u1": 0.001, "u2": 0.5, "lvl": 0.5}           # Holm-adjusted 0.003 < 0.05
    p_ok = {"u1": 0.02, "u2": 0.5, "lvl": 0.5}             # raw 0.02, Holm-adjusted 0.06 >= 0.05: adequate
    units = {"u1": _within("u1"), "u2": _within("u2"), "lvl": pdn.UnitNoiseModel(
        "lvl", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "lvl", df=1.0),))}
    T = {"full": 9.0, "-u1": 4.5, "-u2": 4.5, "u1": 4.5, "u2": 4.5}
    st = {k: RNG.normal(-5.0, 1.0, 1000) for k in T}
    nb, na = _null("b", st), _null("a", {k: RNG.normal(5.0, 1.0, 1000) for k in T})
    v = pdn.judge_pair({"base": _obs(T, p_a=p_bad, units=units)}, na, nb, units, ALL_INPUTS)
    assert v.outcome == "undecided" and not v.details["per_form"]["base"]["a"]["adequate"]
    assert v.details["unit_roles"]["falsification_only"] == ["lvl"]
    v = pdn.judge_pair({"base": _obs(T, p_a=p_ok, units=units)}, na, nb, units, ALL_INPUTS)
    assert v.outcome == "credited" and v.credited == "a"


def test_a_credit_resting_on_one_unit_needs_that_unit_alone():
    units = {"u1": _within("u1"), "u2": _within("u2")}
    T = {"full": 3.0, "-u1": -2.0, "-u2": 3.0, "u1": 0.5, "u2": 0.1}
    nb = _null("b", {"full": RNG.normal(0, 1, 1000), "-u1": RNG.normal(0, 1, 1000), "-u2": RNG.normal(0, 1, 1000),
                     "u1": RNG.normal(0, 1, 1000), "u2": RNG.normal(0, 1, 1000)})
    na = _null("a", {k: RNG.normal(-10, 1, 1000) for k in T})          # b never crosses its mirror threshold
    v = pdn.judge_pair({"base": _obs(T, units=units)}, na, nb, units, ALL_INPUTS)
    assert v.details["per_form"]["base"]["a"]["crosses"] and v.details["per_form"]["base"]["a"]["adequate"]
    loo = v.details["per_form"]["base"]["a"]["leave_one_unit_out"]
    assert loo["u1"]["rests_on"] and not loo["u1"]["alone_credits"] and not loo["u2"]["rests_on"]
    assert v.outcome == "undecided"


def test_sole_survivor_is_carried_not_credited():
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS,
                       falsified={"b": True}, incumbent="b")
    assert v.outcome == "carried_sole_survivor" and v.credited is None and v.carried == ("a",)
    assert v.details["s3_refusal_applies"]
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, None, None, ONE_UNIT, ALL_INPUTS,
                       falsified={"a": True, "b": True})
    assert v.outcome == "both_falsified" and v.carried == ()


def test_an_unresolved_density_is_not_computable():
    v = pdn.judge_pair({"base": _obs({"full": 8.0}, accepted=False)}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS)
    assert v.outcome == "not_computable" and v.credited is None and v.carried == ("a", "b")
    assert any("unresolved density" in r for r in v.reasons)


@pytest.mark.parametrize("kind", ["forced_inputs", "few_replicates", "inside_ensemble", "unresolved_replicate"])
def test_an_incomplete_null_is_not_computable(kind):
    inputs = dict(ALL_INPUTS)
    nb = NULL_B
    if kind == "forced_inputs":
        inputs["forced_inputs"] = True                      # the case has forced inputs; the null did not redraw them
    elif kind == "few_replicates":
        nb = _null("b", {"full": RNG.normal(-5, 1, 999)}, n=999)
    elif kind == "inside_ensemble":
        nb = _null("b", {"full": RNG.normal(-5, 1, 1000)}, outside=False)
    else:
        res = np.ones(1000, bool)
        res[17] = False
        nb = pdn.NullDistribution("b", NULL_B.stats, NULL_B.redrawn, True, res)
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, nb, ONE_UNIT, inputs)
    assert v.outcome == "not_computable", v.reasons


def test_case_inputs_must_be_stated_in_full():
    with pytest.raises(MeasurementModelRefusal, match="every listed null input"):
        pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, {"decisive_data": True})


def test_no_credit_unit_is_not_decisive():
    lvl = {"lvl": pdn.UnitNoiseModel("lvl", "lvl", "absolute", (pdn.NoiseTerm("reproducibility", 0.4, "lvl", df=1.),))}
    v = pdn.judge_pair({"base": _obs({"full": 8.0}, units=("lvl",))}, NULL_A, NULL_B, lvl, ALL_INPUTS)
    assert v.outcome == "not_decisive"


def test_the_power_check_gate():
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS,
                       misspecified_false_credit=(25, 200))         # Wilson lower 0.086 > 0.05: cannot credit
    assert v.outcome == "undecided" and "cannot credit" in v.reasons[0]
    v = pdn.judge_pair({"base": _obs({"full": 8.0})}, NULL_A, NULL_B, ONE_UNIT, ALL_INPUTS,
                       misspecified_false_credit=(12, 200))
    assert v.outcome == "credited"


def test_required_agreement_needs_every_scored_form():
    forms = [pdn.FormDeclaration("f1", "modelling_choice"), pdn.FormDeclaration("f2", "modelling_choice"),
             pdn.FormDeclaration("hedge", "omega_hedge")]
    na = pdn.NullDistribution("a", {f: {"full": RNG.normal(5, 1, 1000)} for f in ("f1", "f2")},
                              frozenset({"decisive_data"}), True)
    nb = pdn.NullDistribution("b", {f: {"full": RNG.normal(-5, 1, 1000)} for f in ("f1", "f2")},
                              frozenset({"decisive_data"}), True)
    v = pdn.judge_pair({"f1": _obs({"full": 8.0}), "f2": _obs({"full": -4.5})}, na, nb, ONE_UNIT, ALL_INPUTS, forms)
    assert v.outcome == "undecided" and not v.details["per_form"]["f2"]["a"]["crosses"]
    assert v.details["per_form"]["f1"]["a"]["qualifies"] and not v.details["qualifies"]["a"]
    v = pdn.judge_pair({"f1": _obs({"full": 8.0}), "f2": _obs({"full": 7.0})}, na, nb, ONE_UNIT, ALL_INPUTS, forms)
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
    credits on every draw except where both rivals qualify (the prototype credits a; the text credits neither)."""
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
                              ALL_INPUTS)
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


# --- planted-truth worlds, from draws --------------------------------------------------------------------------------

# The world: two within-run shapes (credit) and one absolute level (falsification only: no discrepancy estimate);
# rival a predicts theta x, rival b theta x^1.35, theta ~ N(1, 0.015) (a weakly identified scale, so every density
# of the pools and the sets resolves with 1 500 prior draws); per-point noise 0.2 at a fixed same-observable within-run
# repeat (a df-carrying SD with few df makes some bootstrap replicates' re-drawn SDs sharp enough to leave ESS < 100:
# the pool is then incomplete and the rule refuses, as it should). Tuned once, before these tests were frozen, to put
# power in the confirmation run's informative-case range (M3 route, 0.70 - 0.90): seed-fixed, 0.81 / 0.77 power and
# 0.04 / 0.04 false credit at R = 100 per truth.
X1, X2 = np.array([0.2, 0.4, 0.6, 0.8, 1.0]), np.array([1.2, 1.4, 1.6, 1.8])
P_EXP, SD, PRIOR_SD, S_DRAWS, POOL, R_SETS = 1.35, 0.2, 0.015, 1500, 1000, 100


def _world_noise(sd=SD, floor=0.02):
    m = {u: pdn.UnitNoiseModel(u, "shape", "within_run", (pdn.NoiseTerm("reproducibility", sd, "shape"),
                                                          pdn.NoiseTerm("instrument_floor", floor, "shape")))
         for u in ("shape1", "shape2")}
    m["level"] = pdn.UnitNoiseModel("level", "level", "absolute",
                                    (pdn.NoiseTerm("reproducibility", 0.1, "level"),
                                     pdn.NoiseTerm("reproducibility", 0.15, "level", "run_offset", 1.0)))
    return m


def _world_draws(rival, n, rng, p_exp=P_EXP):
    th = rng.normal(1.0, PRIOR_SD, n)
    p = 1.0 if rival == "a" else p_exp
    return pdn.RivalDraws({"shape1": th[:, None] * X1 ** p, "shape2": th[:, None] * X2 ** p,
                           "level": np.column_stack([th, th])})


@pytest.fixture(scope="module")
def planted():
    design = pdn.MeasurementDesign(_world_noise())
    rng = np.random.default_rng(1)
    ea, eb = {"base": {"base": _world_draws("a", S_DRAWS, rng)}}, {"base": {"base": _world_draws("b", S_DRAWS, rng)}}
    ta, tb = _world_draws("a", S_DRAWS, np.random.default_rng(2)), _world_draws("b", S_DRAWS, np.random.default_rng(3))
    pa = pdn.build_null_pool("a", {"base": ta}, design, ea, eb, POOL, np.random.default_rng(4))
    pb = pdn.build_null_pool("b", {"base": tb}, design, ea, eb, POOL, np.random.default_rng(5))
    rates = {}
    for truth, td in (("a", ta), ("b", tb)):
        r = np.random.default_rng(10 if truth == "a" else 11)
        tally = {}
        for _ in range(R_SETS):
            j = int(r.integers(td.n))
            y = {u: design.noise[u].sample(td.pred(u)[j][None, :], r)[0] for u in design.noise}
            v = pdn.judge_from_draws(y, design, ea, eb, pa, pb, r, ppc_reps=100)
            key = v.credited or v.outcome
            tally[key] = tally.get(key, 0) + 1
        rates[truth] = tally
    return {"design": design, "ea": ea, "eb": eb, "pa": pa, "pb": pb, "rates": rates}


def test_planted_truth_power(planted):
    """The true rival is credited at about the rate the confirmation run reported where the data are informative
    (M3 route: 0.70 - 0.90, the estimated diagnostic): this world is built to sit there."""
    ra, rb = planted["rates"]["a"], planted["rates"]["b"]
    power_a, power_b = ra.get("a", 0) / R_SETS, rb.get("b", 0) / R_SETS
    assert 0.6 <= power_a <= 0.97 and 0.6 <= power_b <= 0.97, planted["rates"]


def test_planted_truth_false_credit_within_the_bound(planted):
    """False credit: the Wilson 95 % interval's lower end does not exceed 0.05 (the stated bound)."""
    for truth, wrong in (("a", "b"), ("b", "a")):
        k = planted["rates"][truth].get(wrong, 0)
        assert pdn.wilson_interval(k, R_SETS)[0] <= 0.05, planted["rates"]


def test_planted_pools_are_complete(planted):
    for p in (planted["pa"], planted["pb"]):
        assert p.truths_outside_scoring_ensemble and p.resolved.all()
        assert {"decisive_data", "replicate_sds"} <= p.redrawn


def test_unresolved_draw_densities_are_not_computable(planted):
    """At a noise far below the draws' spread the joint ESS collapses: the verdict refuses rather than reads it."""
    design = pdn.MeasurementDesign(_world_noise(sd=0.0005, floor=0.0001))
    rng = np.random.default_rng(21)
    y = {u: planted["ea"]["base"]["base"].pred(u)[3] + 0.001 for u in design.noise}
    obs = pdn.score_from_draws(y, design, planted["ea"], planted["eb"], rng, ppc_reps=50)
    assert not obs["base"].lppd["a"]["full"].accepted
    v = pdn.judge_pair(obs, pdn.null_from_pool(planted["pa"]), pdn.null_from_pool(planted["pb"]),
                       {u: design.noise[u] for u in design.decisive_units()}, design.case_inputs())
    assert v.outcome == "not_computable" and any("unresolved density" in r for r in v.reasons)


def test_a_pool_without_forced_inputs_is_incomplete(planted):
    design = pdn.MeasurementDesign(_world_noise(), forced_inputs_present=True)
    rng = np.random.default_rng(22)
    ta = _world_draws("a", 500, np.random.default_rng(23))
    tb = _world_draws("b", 500, np.random.default_rng(24))
    pa = pdn.build_null_pool("a", {"base": ta}, design, planted["ea"], planted["eb"], 20, rng)
    pb = pdn.build_null_pool("b", {"base": tb}, design, planted["ea"], planted["eb"], 20, rng)
    assert "forced_inputs" not in pa.redrawn
    y = {u: ta.pred(u)[0] for u in design.noise}
    v = pdn.judge_from_draws(y, design, planted["ea"], planted["eb"], pa, pb, rng, ppc_reps=50)
    assert v.outcome == "not_computable"
    assert any("forced_inputs" in r for r in v.reasons) and any("effective replicates" in r for r in v.reasons)


def test_truths_inside_the_scoring_ensemble_are_flagged(planted):
    design = planted["design"]
    p = pdn.build_null_pool("a", {"base": planted["ea"]["base"]["base"]}, design, planted["ea"], planted["eb"], 3,
                            np.random.default_rng(25))
    assert not p.truths_outside_scoring_ensemble


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
    pool = pdn.build_null_pool("a", {"base": dr("a", 600, np.random.default_rng(32))}, design, ea, eb, 40, rng)
    assert "conditioning_rows" in pool.redrawn
    nl = pdn.null_from_pool(pool, {"spent": np.array([1.1])})
    assert nl.weights is not None and 1.0 < nl.n_effective < 40 and "weights_k_hat" in nl.extra
