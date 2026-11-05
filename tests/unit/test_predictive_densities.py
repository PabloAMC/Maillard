"""src/kinetic_core/predictive_densities.py (plan step M1.1): PSIS on known targets, the censored likelihood, Holm,
the adequacy rule, the refusals, the pH discrepancy, and the regression against M1's frozen records (EX-061)."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate, stats

from src.kinetic_core import predictive_densities as pdn
from src.kinetic_core.predictive_densities import MeasurementModelRefusal

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / "results" / "validation"


# --- PSIS and importance sampling on known targets ------------------------------------------------------------------

def _is_estimate(target, proposal, n, seed, log_c=0.0, smooth=True):
    rng = np.random.default_rng(seed)
    x = proposal.rvs(size=n, random_state=rng)
    return pdn.importance_lppd(log_c + target.logpdf(x), proposal.logpdf(x), smooth=smooth)


def test_psis_khat_and_smoothed_weights_on_a_known_pareto_tail():
    """Log ratios whose ratios ARE generalised-Pareto with shape 0.5: k-hat recovers 0.5; the body's log weights are
    returned unchanged; the smoothed tail never exceeds the largest raw ratio."""
    rng = np.random.default_rng(3)
    k_true = 0.5
    w = stats.genpareto(c=k_true).rvs(size=8000, random_state=rng) + 1e-3
    lr = np.log(w)
    lw, k = pdn.psis(lr)
    assert abs(k - k_true) < 0.12
    order = np.argsort(lr)
    m = int(math.ceil(min(lr.size / 5.0, 3.0 * math.sqrt(lr.size))))
    body = order[: lr.size - m - 1]
    assert np.allclose(lw[body], lr[body], rtol=0, atol=1e-12)
    assert np.max(lw) <= np.max(lr) + 1e-12
    # a light tail (bounded ratios) gives a small, even negative, k-hat
    _, k_light = pdn.psis(np.log(rng.uniform(0.5, 1.0, size=4000)))
    assert k_light < 0.0


def test_gaussian_proposal_for_a_heavier_tailed_target_is_flagged():
    """q = N(0, 1), p = Student-t(3, scale 3): the ratios p / q are unbounded (infinite variance); over 20 fixed
    seeds the median k-hat lands above 0.7 and most estimates are refused. (A finite-sample k-hat scatters by about
    +-0.15 at S 4 000, so one seed alone does not decide the regime.)"""
    ests = [_is_estimate(stats.t(df=3, scale=3), stats.norm(0, 1), 4000, seed=s) for s in range(20)]
    assert np.median([e.k_hat for e in ests]) > pdn.K_HAT_MAX
    assert sum(not e.accepted for e in ests) >= 12
    assert all("k-hat" in e.reason for e in ests if not e.accepted)


def test_heavy_proposal_for_a_gaussian_target_converges_to_the_analytic_value():
    """q = Student-t(3, scale 1.5), integrand = c N(0.3, 0.8): bounded ratios, k-hat < 0.5, the estimate converges to
    log c and the error shrinks with S (within 3 MC SE at every S)."""
    log_c = -2.7
    errs = []
    for n, seed in ((500, 1), (4000, 2), (32000, 3)):
        est = _is_estimate(stats.norm(0.3, 0.8), stats.t(df=3, scale=1.5), n, seed, log_c=log_c)
        assert est.k_hat < 0.5 and est.accepted
        assert abs(est.lppd - log_c) < 3.0 * est.mc_se
        errs.append(est.mc_se)
    assert errs[0] > errs[1] > errs[2]
    raw = _is_estimate(stats.norm(0.3, 0.8), stats.t(df=3, scale=1.5), 32000, 3, log_c=log_c, smooth=False)
    assert abs(raw.lppd - log_c) < 3.0 * raw.mc_se


def test_ess_and_the_precision_rule_decide_acceptance():
    est = _is_estimate(stats.norm(0, 1), stats.t(df=5, scale=1.2), 300, seed=5)
    assert est.accepted and est.ess >= 100
    need = pdn.tie_band_precision_ess(sd_rep=0.05)          # 2 / (0.15 x 0.05)^2 = 35 556 draws' worth
    assert need == pytest.approx(2.0 / (0.15 * 0.05) ** 2)
    rng = np.random.default_rng(5)
    x = stats.t(df=5, scale=1.2).rvs(size=300, random_state=rng)
    tight = pdn.importance_lppd(stats.norm(0, 1).logpdf(x), stats.t(df=5, scale=1.2).logpdf(x), precision_ess=need)
    assert not tight.accepted and "ESS" in tight.reason
    assert pdn.effective_sample_size(np.zeros(50)) == pytest.approx(50.0)


def test_importance_lppd_refuses_nan_and_counts_zero_weight_draws():
    with pytest.raises(MeasurementModelRefusal):
        pdn.importance_lppd([0.0, float("nan")], [0.0, 0.0])
    li = np.array([0.0] * 9 + [-np.inf])
    est = pdn.importance_lppd(li, np.zeros(10), smooth=False)
    assert est.n_nonfinite == 1 and est.lppd == pytest.approx(math.log(0.9))


def test_multivariate_t_proposal_density_matches_scipy():
    rng = np.random.default_rng(0)
    x, lq = pdn.multivariate_t_draws([1.0, -1.0], np.array([[0.5, 0.1], [0.1, 0.3]]), 200, rng)
    assert x.shape == (200, 2)
    ref = stats.multivariate_t(loc=[1.0, -1.0], shape=[[0.5, 0.1], [0.1, 0.3]], df=5).logpdf(x)
    assert np.allclose(lq, ref)


# --- quadrature and the plug-in / linearised forms ----------------------------------------------------------------

def test_quadrature_is_exact_on_gaussian_integrands_and_refuses_three_dimensions():
    log_c = 1.3
    g1 = lambda x: log_c + stats.norm(0.4, 0.2).logpdf(x[:, 0])
    assert pdn.quadrature_lppd(g1, [0.35], [[0.05]], n_nodes=30).lppd == pytest.approx(log_c, abs=1e-8)
    cov = np.array([[0.3, 0.1], [0.1, 0.2]])
    g2 = lambda x: log_c + stats.multivariate_normal([0.1, -0.2], cov).logpdf(x)
    assert pdn.quadrature_lppd(g2, [0.0, 0.0], cov, n_nodes=30).lppd == pytest.approx(log_c, abs=1e-6)
    with pytest.raises(MeasurementModelRefusal):
        pdn.quadrature_lppd(lambda x: np.zeros(len(x)), [0, 0, 0], np.eye(3))


def test_linearised_equals_the_exact_density_when_the_prediction_is_linear():
    """f(theta) = f0 + J theta, theta ~ N(0, S): the predictive is exactly N(f0, Sigma + J S J^T)."""
    y, f0, sig = np.array([0.3]), np.array([0.0]), np.array([[0.04]])
    j, s = np.array([[2.0]]), np.array([[0.01]])
    lin = pdn.linearised_lppd(y, f0, sig, j, s).lppd
    g = lambda th: np.array([pdn.gaussian_log_density(y, f0 + j @ t, sig) for t in th]) + stats.norm(0, 0.1).logpdf(th[:, 0])
    assert pdn.quadrature_lppd(g, [0.0], s, n_nodes=40).lppd == pytest.approx(lin, abs=1e-10)


# --- likelihood pieces ----------------------------------------------------------------------------------------------

def test_dex_and_ln_densities_differ_by_n_ln_ln10():
    y, f = np.array([0.5, 0.1, -0.2]), np.array([0.4, 0.0, 0.0])
    cov_dex = pdn.series_covariance([0.05, 0.05, 0.05], "absolute", "dex", omega_dex=0.39)
    cov_ln = pdn.series_covariance(pdn.dex_to_ln([0.05, 0.05, 0.05]), "absolute", "ln", omega_dex=0.39)
    assert np.allclose(cov_ln, cov_dex * pdn.LN10 ** 2)
    l_dex = pdn.gaussian_log_density(y, f, cov_dex)
    l_ln = pdn.gaussian_log_density(pdn.dex_to_ln(y), pdn.dex_to_ln(f), cov_ln)
    assert l_ln == pytest.approx(pdn.log_density_dex_to_ln(l_dex, 3))
    assert pdn.log_density_ln_to_dex(l_ln, 3) == pytest.approx(l_dex)
    # the default is H0's ship value (the worse set, 0.419 dex); the reset set's 0.390 is the named sensitivity
    assert (pdn.OMEGA_DEX, pdn.OMEGA_DEX_RESET) == (0.419, 0.390)
    assert pdn.OMEGA_LN == pytest.approx(0.965, abs=1e-3)
    assert pdn.dex_to_ln(pdn.OMEGA_DEX_RESET) == pytest.approx(0.898, abs=1e-3)   # M1c + M2 prereg sec. 6's figure


def test_series_covariance_ratio_design_and_omega():
    """M1's D1 unit: ratios 360/30 and 720/360 share the 360 min point; omega enters each predicted point."""
    v = {30: 0.02 ** 2, 360: 0.03 ** 2, 720: 0.01 ** 2}
    sd = np.sqrt([v[30], v[360], v[720]])
    a = np.array([[-1.0, 1.0, 0.0], [0.0, -1.0, 1.0]])
    m1 = np.array([[v[360] + v[30], -v[360]], [-v[360], v[720] + v[360]]])
    assert np.allclose(pdn.series_covariance(sd, "ratio", omega_dex=None, design=a), m1)
    om = 0.419 ** 2                                          # the default omega
    assert np.allclose(pdn.series_covariance(sd, "ratio", design=a), m1 + om * np.array([[2.0, -1.0], [-1.0, 2.0]]))
    om = 0.390 ** 2                                          # the reset-set sensitivity, passed explicitly
    assert np.allclose(pdn.series_covariance(sd, "ratio", omega_dex=pdn.OMEGA_DEX_RESET, design=a),
                       m1 + om * np.array([[2.0, -1.0], [-1.0, 2.0]]))


def test_censored_likelihood_against_analytic_cases():
    # one value, below its floor: log Phi((log floor - f) / sd)
    got = pdn.censored_gaussian_log_density([np.nan], [0.2], [[0.09]], [True], [0.0])
    assert got == pytest.approx(stats.norm.logcdf((0.0 - 0.2) / 0.3))
    # independent values: the observed density plus the censored cdf
    cov = np.diag([0.04, 0.09])
    got = pdn.censored_gaussian_log_density([0.1, np.nan], [0.0, 0.5], cov, [False, True], [0.0, 0.3])
    assert got == pytest.approx(stats.norm(0.0, 0.2).logpdf(0.1) + stats.norm.logcdf((0.3 - 0.5) / 0.3))
    # correlated (a shared level offset): against numerical integration of the bivariate density
    cov = np.array([[0.05, 0.03], [0.03, 0.06]])
    y_o, f, lf = 0.25, np.array([0.1, 0.4]), 0.2
    dens = stats.multivariate_normal(f, cov)
    num, _ = integrate.quad(lambda z: dens.pdf([y_o, z]), -np.inf, lf)
    got = pdn.censored_gaussian_log_density([y_o, np.nan], f, cov, [False, True], [np.nan, lf])
    assert got == pytest.approx(math.log(num), abs=1e-7)
    # two censored values: the bivariate cdf
    got = pdn.censored_gaussian_log_density([np.nan, np.nan], f, cov, [True, True], [0.0, 0.3])
    assert got == pytest.approx(math.log(dens.cdf([0.0, 0.3])), abs=1e-5)
    with pytest.raises(MeasurementModelRefusal):
        pdn.censored_gaussian_log_density([np.nan], [0.0], [[1.0]], [True], [np.nan])


def test_assay_floor_and_the_floor_at_a_datum_refusal():
    assert pdn.assay_floor(lod=0.02, loq=0.07) == 0.02 and pdn.assay_floor(lod=0.02, loq=0.07, use="loq") == 0.07
    with pytest.raises(MeasurementModelRefusal):
        pdn.assay_floor(loq=0.07)
    # M1's primary floors: Cys 0.80 mM against the 50-min datum 0.8016, PENT 1.22 against 1.22
    with pytest.raises(MeasurementModelRefusal, match="equals the data value"):
        pdn.check_floor_not_data(0.80, [20.26, 4.01, 0.8016])
    with pytest.raises(MeasurementModelRefusal):
        pdn.check_floor_not_data(1.22, [18.52, 1.22])
    pdn.check_floor_not_data(0.020, [20.26, 4.01, 0.8016])        # Zhai 2021's Xyl LOD: accepted


def test_tau_is_refused_on_relative_value_bases():
    for basis in ("c_over_c0", "ratio", "first_point_normalised", "contrast", "share"):
        with pytest.raises(MeasurementModelRefusal, match="tau refused"):
            pdn.series_covariance([0.03, 0.03], basis, tau_dex=0.964)
    cov = pdn.series_covariance([0.03, 0.03], "absolute", omega_dex=None, tau_dex=0.964)
    assert np.allclose(cov, np.diag([0.0009, 0.0009]) + 0.964 ** 2)
    with pytest.raises(MeasurementModelRefusal):
        pdn.series_covariance([0.03], "levels")
    with pytest.raises(MeasurementModelRefusal):
        pdn.series_covariance([0.03], "absolute", basis="log2")


# --- decision helpers -------------------------------------------------------------------------------------------------

def test_holm_against_a_hand_case():
    out = pdn.holm({"a": 0.01, "b": 0.04, "c": 0.03, "d": 0.005}, alpha=0.05)
    # sorted d .005 <= .05/4, a .01 <= .05/3, c .03 > .05/2: stop; b not rejected
    assert {g: v["rejected"] for g, v in out.items()} == {"a": True, "b": False, "c": False, "d": True}
    assert {g: round(v["p_holm"], 10) for g, v in out.items()} == {"d": 0.02, "a": 0.03, "c": 0.06, "b": 0.06}
    assert pdn.holm({"x": 0.06, "y": 0.2}, alpha=0.10)["x"]["rejected"] is False   # .06 > .10 / 2


def test_adequacy_rule_three_branches():
    # |d| within the band: a tie, whatever the adequacy
    assert pdn.unit_reading(1.0, 2.0, 0.5, 0.5, [0.0], [0.0], 0.39) == "tie"
    # favours a: past the band and a adequate
    assert pdn.unit_reading(5.0, 2.0, 0.30, 0.0, [0.0], [0.1], 0.39) == "a"
    # favours b on a mean gap > omega although b is inadequate
    assert pdn.unit_reading(-5.0, 2.0, 0.0, 0.001, [0.0, 0.0], [0.1, 0.6], 0.39) == "b"
    # past the band, the winner inadequate and the means within omega: both inadequate
    assert pdn.unit_reading(5.0, 2.0, 0.01, 0.0, [0.0], [0.2], 0.39) == "both inadequate"


def test_posterior_predictive_p_is_the_chi2_tail():
    assert pdn.posterior_predictive_p([1.96], [0.0], [[1.0]]) == pytest.approx(0.05, abs=1e-3)
    assert pdn.posterior_predictive_p([1.0, 1.0], [0.0, 0.0], np.eye(2)) == pytest.approx(stats.chi2.sf(2.0, 2))


def test_tie_band_by_simulation_matches_the_linear_case():
    rng = np.random.default_rng(9)
    y, cov, w = np.array([0.1, 0.2]), np.diag([0.04, 0.01]), np.array([3.0, -2.0])
    reps = pdn.gaussian_replicates(y, cov, 4000, rng)
    tb = pdn.tie_band(lambda yy: float(w @ yy), reps)
    sd = math.sqrt(w @ cov @ w)
    assert tb.sd_rep == pytest.approx(sd, rel=0.05) and tb.band == pytest.approx(pdn.Z_90 * tb.sd_rep)
    assert tb.n_replicates == 4000


def test_required_agreement_over_forms():
    bands = {"primary": 1.0, "plug_in": 1.0, "floors_alt": 0.5}
    assert pdn.required_agreement({"primary": 3.0, "plug_in": 2.0, "floors_alt": 0.6}, bands)["side"] == "a"
    assert pdn.required_agreement({"primary": -3.0, "plug_in": -2.0, "floors_alt": -0.6}, bands)["side"] == "b"
    out = pdn.required_agreement({"primary": 3.0, "plug_in": 0.4, "floors_alt": 0.6}, bands)
    assert out["side"] == "tie" and out["breaking_forms"] == ["plug_in"]
    assert pdn.required_agreement({"primary": 3.0, "plug_in": -2.0, "floors_alt": 0.6}, bands)["side"] == "tie"
    with pytest.raises(MeasurementModelRefusal):
        pdn.required_agreement({"primary": 3.0, "extra": 1.0}, {"primary": 1.0})


# --- pH ---------------------------------------------------------------------------------------------------------------

def _toy_balance(amounts):
    """A strong acid HA and a strong base B in water at 25 C: [H+] + B = HA + Kw / [H+]; pH = -log10 [H+] (mol/L)."""
    d = amounts["B"] - amounts["HA"]
    h = (-d + math.sqrt(d * d + 4e-14)) / 2.0
    return -math.log10(h)


def test_ph_species_discrepancy_is_the_half_spread_of_the_corners():
    produced = {"HA": 2e-3, "B": 1e-3}
    out = pdn.ph_species_discrepancy(_toy_balance, produced, acids={"HA"}, per_species=True)
    assert out["omega_dex"] == 0.419                       # the default: H0's ship value
    up, dn = 10 ** 0.419, 10 ** -0.419
    ph_a = _toy_balance({"HA": 2e-3 * up, "B": 1e-3 * dn})
    ph_b = _toy_balance({"HA": 2e-3 * dn, "B": 1e-3 * up})
    assert out["delta_ph"] == pytest.approx(0.5 * abs(ph_b - ph_a))
    assert out["corners_bracket_central"] and out["ph_central"] == pytest.approx(_toy_balance(produced))
    # the corners bound every independent +-omega combination
    combos = [_toy_balance({"HA": 2e-3 * s1, "B": 1e-3 * s2}) for s1 in (up, dn) for s2 in (up, dn)]
    assert out["delta_ph"] == pytest.approx(0.5 * (max(combos) - min(combos)))
    assert set(out["per_species_half_spread"]) == {"HA", "B"}
    with pytest.raises(MeasurementModelRefusal, match="not net"):
        pdn.ph_species_discrepancy(_toy_balance, {"HA": -1e-3, "B": 1e-3}, acids={"HA"})
    with pytest.raises(MeasurementModelRefusal):
        pdn.ph_species_discrepancy(_toy_balance, produced, acids={"ACID"})


def test_mixed_ph_scale_is_refused():
    assert pdn.require_ph_scale("cooled", "cooled") == "cooled"
    with pytest.raises(MeasurementModelRefusal, match="mixed-scale"):
        pdn.require_ph_scale("cooled", "in_situ")
    with pytest.raises(MeasurementModelRefusal):
        pdn.require_ph_scale("25C")
    with pytest.raises(MeasurementModelRefusal):
        pdn.ph_difference(pdn.PhReading(5.1, "in_situ", 120.0), pdn.PhReading(5.3, "cooled", 25.0))
    assert pdn.ph_difference(pdn.PhReading(5.1, "cooled"), pdn.PhReading(5.3, "cooled")) == pytest.approx(-0.2)
    with pytest.raises(MeasurementModelRefusal):
        pdn.ph_species_discrepancy(_toy_balance, {"HA": 1e-3, "B": 1e-3}, acids={"HA"}, scale="at_T")


# --- the key regression: M1's numbers from its frozen records -----------------------------------------------------

#: M1's settings (fit_m1_supply_combinations.py): the reservoir pot, its primary floors (data points) and the
#: alternative floors (the ELSD curve's zero crossing), tau_ABS(all)
RES_POT = "zhang2025b_blank_fig2gh"
RES_FLOOR = {"Cys": 0.80, "PENT": 1.22}
RES_FLOOR_ALT = {"Cys": 0.22, "PENT": 0.22}


@pytest.fixture(scope="module")
def m1():
    from src.kinetic_core import fit_rows

    dec = json.loads((V / "kinetic_core_m1_decisive_comparison.json").read_text())
    carried = {r: json.loads((V / f"kinetic_core_m1_supply_combinations_carried_{r}.json").read_text())
               for r in ("R0", "C1")}
    _, rows = fit_rows.load("route_submodel")
    data = {s: [r for r in rows if r.pot == RES_POT and r.species == s] for s in ("Cys", "PENT")}
    return dec, carried, data


def _sigma_theta(c):
    """M1's ``identified_parts``: identified, not-at-bound coordinates with their Laplace covariance, bound-active ones
    with their profile variance (uncorrelated), in the problem's key order."""
    lap, prof = c["laplace"], c.get("profiles") or {}
    keys = list(lap["sigma"])
    ident = [i for i, k in enumerate(keys) if lap["identified"][k] and not lap["at_bound"][k]]
    bound = [i for i, k in enumerate(keys) if lap["at_bound"][k] and "variance" in prof.get(k, {})]
    cov = np.asarray(lap["covariance"], dtype=float)
    out = np.zeros((len(ident) + len(bound),) * 2)
    out[:len(ident), :len(ident)] = cov[np.ix_(ident, ident)]
    for n, i in enumerate(bound):
        out[len(ident) + n, len(ident) + n] = prof[keys[i]]["variance"]
    return out


def _old_unit(dec, carried, data, rival, species, floors, tau, reading="isothermal"):
    """One reservoir unit under M1's OLD measurement model, built with this module: the prediction floored AT the floor
    (zero Jacobian row where floored), tau on the C/C0 series (the override), no omega; plug-in and linearised."""
    lv = dec["levels"][rival][reading]
    idx = {lbl: i for i, lbl in enumerate(dec["point_labels"])}
    rs = data[species]
    pts = [idx[f"RES|{species}|{r.t_min:g}"] for r in rs]
    levels = np.asarray(lv["levels"], dtype=float)[pts]
    jac = np.asarray(lv["jac"], dtype=float)[pts]
    f, mask = pdn.legacy_floored_prediction(levels, math.log10(floors[species]))
    jac = np.where(mask[:, None], 0.0, jac)
    y = np.log10([r.value for r in rs])
    cov = pdn.series_covariance([float(r.sigma) for r in rs], "c_over_c0", "dex", omega_dex=None, tau_dex=tau,
                                legacy_tau_override=True)
    return {"plug_in": pdn.plug_in_lppd(y, f, cov).lppd,
            "linearised": pdn.linearised_lppd(y, f, cov, jac, _sigma_theta(carried[rival])).lppd,
            "y": y, "f": f, "levels": levels}


def test_m1_old_settings_reproduce_the_decisive_table_for_res_cys(m1):
    """Under M1's old settings (linearised density, data-point floors, tau on C/C0) the module reproduces the frozen
    decisive comparison's R0 vs C1 numbers for RES_Cys (isothermal) to the record's precision."""
    dec, carried, data = m1
    tau = dec["tau_ABS_all"]
    settings = {"primary": (RES_FLOOR, tau), "floors_alternative": (RES_FLOOR_ALT, tau),
                "tau_x0.5": (RES_FLOOR, 0.5 * tau), "tau_x2": (RES_FLOOR, 2.0 * tau)}
    got = {}
    for name, (floors, t) in settings.items():
        u = {r: _old_unit(dec, carried, data, r, "Cys", floors, t) for r in ("R0", "C1")}
        rec = {r: dec["scored"][r]["isothermal"][name]["RES_Cys"] for r in ("R0", "C1")}
        for r in ("R0", "C1"):
            assert np.allclose(u[r]["y"], rec[r]["y"], atol=1e-12)
            assert np.allclose(u[r]["f"], rec[r]["f"], atol=1e-12)
            assert u[r]["linearised"] == pytest.approx(rec[r]["linearised"], rel=1e-9)
            assert u[r]["plug_in"] == pytest.approx(rec[r]["plug_in"], rel=1e-9)
        got[name] = {k: u["C1"][k] - u["R0"][k] for k in ("linearised", "plug_in")}
    table = dec["per_reading"]["isothermal"]["units"]["RES_Cys"]["d"]
    assert got["primary"]["linearised"] == pytest.approx(table["linearised"], rel=1e-9)        # +638.01
    assert got["primary"]["plug_in"] == pytest.approx(table["plug_in"], rel=1e-9)              # +375.31
    assert got["floors_alternative"]["linearised"] == pytest.approx(table["floors_alternative"], rel=1e-9)  # -322.38
    assert got["tau_x0.5"]["linearised"] == pytest.approx(table["tau_x0.5"], rel=1e-9)         # +638.75
    assert got["tau_x2"]["linearised"] == pytest.approx(table["tau_x2"], rel=1e-9)             # +637.83
    # the review's deciding number: R0's linearised RES_Cys lppd under the alternative floor
    r0_alt = _old_unit(dec, carried, data, "R0", "Cys", RES_FLOOR_ALT, tau)
    assert r0_alt["linearised"] == pytest.approx(-175.4, abs=0.05)
    # and the new rules refuse both of the old choices
    with pytest.raises(MeasurementModelRefusal):
        pdn.series_covariance([0.0261] * 8, "c_over_c0", tau_dex=tau)
    with pytest.raises(MeasurementModelRefusal):
        pdn.check_floor_not_data(RES_FLOOR["Cys"], [r.value for r in data["Cys"]])
    with pytest.raises(MeasurementModelRefusal):
        pdn.check_floor_not_data(RES_FLOOR["PENT"], [r.value for r in data["PENT"]])


#: The engine scan (scripts/generators/scan_m1_r0_cysteine_thermal.py): R0's decisive levels at k_cys_thermal + delta
#: dex, delta = -1.2 ... +1.2 in steps of 0.1, the other coordinates at R0's carried M1 optimum, M1's scoring numerics.
SCAN = ROOT / "tests" / "fixtures" / "m1_r0_cysteine_thermal_scan.json"


@pytest.fixture(scope="module")
def scan():
    return json.loads(SCAN.read_text())


def test_m1_exact_density_of_r0_cys_under_the_alternative_floor_from_the_engine_scan(m1, scan):
    """The review's finding, from engine values: R0's RES_Cys predictive under the alternative floor (the old form
    otherwise: prediction floored at 0.22 mM, tau 0.964 on the series, no omega), integrated over the one identified
    coordinate that moves R0's cysteine (k_cys_thermal, Laplace sd 0.279), is about -746, not the linearised -175.
    A trapezoid on the 25 engine points (the review's construction) and, between them, a cubic spline in delta of each
    log level feeding a fine trapezoid, PSIS and adaptive Gauss-Hermite: all agree. The integrand has a kink at its
    peak (the legacy-floored 50 min prediction crosses 0.22 mM), and PSIS refuses a proposal scaled from the local
    Hessian there."""
    from scipy.interpolate import CubicSpline
    from scipy.optimize import minimize_scalar

    dec, carried, data = m1
    tau = dec["tau_ABS_all"]
    sd = carried["R0"]["laplace"]["sigma"]["k_cys_thermal"]
    assert scan["laplace_sd"] == sd and scan["failed_pot_runs"] == 0
    assert scan["carried_theta"] == carried["R0"]["carried_theta"]["k_cys_thermal"]
    assert scan["point_labels"] == dec["point_labels"] and scan["numerics"] == dec["forms"]["R0"]["numerics"]
    deltas = np.asarray(scan["deltas"])
    levels = np.asarray(scan["levels"], dtype=float)
    # the scan's centre IS the record's optimum: the engine reproduces M1's levels exactly
    assert np.array_equal(levels[int(np.argmin(np.abs(deltas)))], np.asarray(dec["levels"]["R0"]["isothermal"]["levels"]))
    idx = {lbl: i for i, lbl in enumerate(scan["point_labels"])}
    rs = data["Cys"]
    cols = [idx[f"RES|Cys|{r.t_min:g}"] for r in rs]
    y = np.log10([r.value for r in rs])
    cov = pdn.series_covariance([float(r.sigma) for r in rs], "c_over_c0", omega_dex=None, tau_dex=tau,
                                legacy_tau_override=True)
    lf = math.log10(RES_FLOOR_ALT["Cys"])

    def plug_of(lv):
        return pdn.gaussian_log_density(y, pdn.legacy_floored_prediction(lv, lf)[0], cov)

    prior = stats.norm(0.0, sd)
    on_grid = np.array([plug_of(lv[cols]) for lv in levels])
    assert on_grid[int(np.argmin(np.abs(deltas)))] == pytest.approx(-966.1, abs=0.1)    # the review's delta-0 value
    g_grid = on_grid + prior.logpdf(deltas)
    coarse = float(np.max(g_grid) + math.log(np.trapezoid(np.exp(g_grid - np.max(g_grid)), deltas)))
    assert coarse == pytest.approx(-746.0, abs=1.5)                    # the review's construction: "about -746"

    spline = CubicSpline(deltas, levels[:, cols], axis=0)

    def g(points):
        d = np.asarray(points, dtype=float)[:, 0]
        return np.array([plug_of(spline(min(max(x, deltas[0]), deltas[-1]))) for x in d]) + prior.logpdf(d)

    xs = np.linspace(deltas[0], deltas[-1], 48001)
    gs = g(xs[:, None])
    ref = float(np.max(gs) + math.log(np.trapezoid(np.exp(gs - np.max(gs)), xs)))
    assert ref == pytest.approx(-746.0, abs=1.5) and ref == pytest.approx(coarse, abs=1.5)
    mode = minimize_scalar(lambda d: -g([[d]])[0], bounds=(0.0, 1.2), method="bounded", options={"xatol": 1e-6}).x
    # the kink: the local inverse Hessian is far narrower than the integrand, and PSIS refuses that proposal
    h = 1e-3
    local = 1.0 / (-(g([[mode + h]])[0] - 2 * g([[mode]])[0] + g([[mode - h]])[0]) / h ** 2)
    assert math.sqrt(local) < 0.01
    rng = np.random.default_rng(746)
    x, lq = pdn.multivariate_t_draws([mode], [[local]], 2000, rng)
    narrow = pdn.importance_lppd(g(x), lq)
    assert not narrow.accepted and narrow.k_hat >= pdn.K_HAT_MAX
    # a proposal as wide as the integrand is accepted and agrees with the reference
    x, lq = pdn.multivariate_t_draws([mode], [[0.05 ** 2]], 2000, rng)
    est = pdn.importance_lppd(g(x), lq)
    assert est.accepted and est.k_hat < pdn.K_HAT_MAX and est.ess >= pdn.ESS_MIN
    assert est.lppd == pytest.approx(ref, abs=max(3.0 * est.mc_se, 0.1))
    # adaptive Gauss-Hermite at the same scale, the preregs' fallback
    quad = pdn.quadrature_lppd(g, [mode], [[0.05 ** 2]], n_nodes=60)
    assert quad.accepted and quad.lppd == pytest.approx(ref, abs=0.5)
    # the prior as proposal (an untuned IS) is not accepted either
    xp = rng.normal(0.0, sd, size=2000)[:, None]
    assert not pdn.importance_lppd(g(xp), prior.logpdf(xp[:, 0])).accepted
    # the linearised value it replaces
    lin = _old_unit(dec, carried, data, "R0", "Cys", RES_FLOOR_ALT, tau)["linearised"]
    assert lin - ref > 500.0
