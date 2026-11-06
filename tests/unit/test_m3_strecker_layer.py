"""M3 phase 1 (src/m3_strecker_layer.py; scripts/generators/fit_m3_strecker_layer.py): the layer's rates, the
acid-load function (V-acid), the O2 budget, the balances (V2), EX-049 at its structure limit (V1), R0-pinned from the
code (V4), the refusals (memory, run-once, freeze, compute cap, the keyed table) and the blind-order check (V-blind on
a scratch git repository). Nothing here reads hofmann2000b Figs. 2 - 3."""
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from src import m3_strecker_layer as L
from src.kinetic_core.predictive_densities import MeasurementModelRefusal

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))
import fit_m3_strecker_layer as G  # noqa: E402


@pytest.fixture(scope="module")
def data():
    return L.load_data()


def _params(data, rival="R1", form=L.PRIMARY, n=4, seed=3):
    u = L.sobol_unit(n, len(L.coords(rival, form)), seed)
    return L.derive(L.from_unit(u, rival, form), rival, form, data)


# ------------------------------------------------------------------------------------------- V-acid and pH

def test_v_acid_reproduces_the_prereg_table_within_0_02_pH():
    rows = L.v_acid_table()
    diffs = [d for r in rows for d in r.get("abs_diff", [])]
    assert diffs and max(diffs) <= 0.02, rows
    assert all(r.get("pass", True) for r in rows)


def test_buffer_ph_is_henderson_hasselbalch_with_a_floor_past_exhaustion():
    ph, ex = L.buffer_ph(np.array([0.0]), 0.5)
    assert ph[0] == pytest.approx(7.0, abs=1e-12) and not ex[0]
    r0 = 10 ** (7.0 - L.PKA2)
    hpo4 = 0.5 * r0 / (1 + r0)
    ph, ex = L.buffer_ph(np.array([0.05]), 0.5)
    assert ph[0] == pytest.approx(L.PKA2 + math.log10((hpo4 - 0.05) / (0.5 - hpo4 + 0.05)), abs=1e-12)
    ph, ex = L.buffer_ph(np.array([hpo4 + 0.01]), 0.5)
    assert ph[0] == L.PH_ACID_FLOOR and ex[0]
    ph, ex = L.buffer_ph(np.array([0.001]), L.PHOSPHATE_UNBUFFERED_M)      # the "0.5 mmol/L" reading
    assert ex[0]


def test_ph_factor_and_the_label_scale():
    assert L.ph_factor(7.0, 0.5) == pytest.approx(1.0)
    assert L.ph_factor(6.0, 0.54) == pytest.approx(10 ** -0.54)
    assert L.PH_SCALE == "label_25C"


# ------------------------------------------------------------------------------------------------ rates

def test_arrhenius_factor_at_85_C_is_the_prereg_0_363():
    assert float(L.arrhenius(np.array([75.0]), L.T85)[0]) == pytest.approx(0.363, abs=0.001)
    assert float(L.arrhenius(np.array([75.0]), L.T100)[0]) == 1.0


def test_dirichlet_stick_breaking_sums_to_one_and_has_uniform_means():
    u = L.sobol_unit(4096, 2, 1)
    x = np.column_stack(L.stick3(u[:, 0], u[:, 1]))
    assert np.allclose(x.sum(axis=1), 1.0) and np.all(x >= 0)
    assert np.allclose(x.mean(axis=0), 1 / 3, atol=0.01)


def test_rival_structures(data):
    p2 = _params(data, "R2")
    assert np.all(p2.g[2] == 0) and np.allclose(p2.g[0] + p2.g[1], 1.0)        # R2: g_B = 0
    p1 = _params(data, "R1")
    assert np.allclose(sum(p1.g), 1.0) and np.all(p1.g[2] > 0)
    p3 = _params(data, "R3")
    assert not np.allclose(p3.k_d7, p3.k_d7_glc)                               # R3: its own glucose-pot decay
    assert np.array_equal(p1.k_d7, p1.k_d7_glc)
    p4 = _params(data, "R4")
    assert p4.ph_on_fed and p4.ph_on_glc and p4.ph_on_trap and np.all((p4.n_order >= 1) & (p4.n_order <= 2))
    assert not p1.ph_on_fed and not p1.ph_on_trap
    assert _params(data, "R1", L.Form(ph_drift_trapped=True)).ph_on_trap
    names = {r: [c.name for c in L.coords(r)] for r in ("R1", "R2", "R3", "R4")}
    assert "log10_kd7_glc" in names["R3"] and "n_order" in names["R4"] and "psi" not in names["R2"]
    assert [c.name for c in L.coords("R0")] == ["cal_mgo", "cal_go"]
    assert min(c.lo for c in L.coords("R1", L.Form(kd_floor=True)) if c.name == "log10_kd7") == L.LOG10_KD_FLOOR


def test_shares_and_the_unmeasured_donors_follow_the_declared_ranges(data):
    p = _params(data, "R2", n=16)
    kap, cal = p.cal["kappa"], p.cal
    s_g, s_mgo = cal["k"]["G"] / kap["G"], cal["k"]["MGO"] / kap["MGO"]
    lo, hi = np.minimum(s_g, s_mgo), np.maximum(s_g, s_mgo)
    s_fox = p.donor["FOX"]["k"] / p.donor["FOX"]["kappa"]
    s_dg1 = p.donor["DG1"]["k"] / p.donor["DG1"]["kappa"]
    assert np.all(s_fox >= lo * (1 - 1e-12)) and np.all(s_fox <= hi * (1 + 1e-12))
    assert np.all(s_dg1 <= hi * 10 ** L.GRANVOGL_DEX * (1 + 1e-12))
    p73 = _params(data, "R2", L.Form(cf073=True), n=16)
    assert np.max(p73.donor["DG1"]["k"] / p73.donor["DG1"]["kappa"]) > np.max(s_dg1)
    assert np.all((p.donor["DG1"]["k0"] >= 0.07 - 1e-12) & (p.donor["DG1"]["k0"] <= 1.0))


def test_first_order_decay_through_the_exits_when_no_donor_reacts(data):
    """With pH off and no order term (R1), O2 = 1 and no Cu, ARP-Phe decays at exactly k_d = k_anox + k_ox."""
    p = _params(data, "R1", n=3)
    sig = np.array([0.1, 0.4, 0.9])
    r = L.solve_arp_pot(p, L.Pot("air_test", "air", 120.0, (60.0, 120.0)), sig, np.ones(3))
    for j, t in enumerate((60.0, 120.0)):
        assert np.allclose(r["ARP"][:, j], p.th["arp0"] * np.exp(-p.k_d7 * t), rtol=1e-10)


def test_argon_slows_the_decay_by_the_oxygen_share(data):
    p = _params(data, "R1", n=2)
    sig = np.array([0.3, 0.6])
    r = L.solve_arp_pot(p, L.Pot("ar", "argon", 120.0, (120.0,)), sig)
    assert np.allclose(r["ARP"][:, 0], p.th["arp0"] * np.exp(-(1 - sig) * p.k_d7 * 120.0), rtol=1e-10)


# --------------------------------------------------------------------------------------------- the O2 budget

def test_o2_capacity_matches_the_prereg_arithmetic():
    assert L.O2_UMOL_PER_ML_AIR == pytest.approx(8.56, abs=0.01)
    assert L.o2_capacity_mM(5.0) == pytest.approx((5 * 8.563 + 2.6) / 10, rel=1e-3)
    assert L.o2_capacity_mM(100.0) == pytest.approx(85.9, abs=0.1)


def test_o2_budget_is_conserved_depletes_at_5_mL_and_tends_to_unlimited(data):
    sig, cu = np.full(4, 0.5), np.full(4, 3.0)
    pot = L.POT_T2["air_cu5"]
    unl = L.solve_arp_pot(_params(data, "R1", n=4), pot, sig, cu)
    small = L.solve_arp_pot(_params(data, "R1", L.Form(headspace_ml=5.0), n=4), pot, sig, cu, ledger=True)
    g = small["ledger"]
    ok = g["o2_clipped"] == 0
    assert ok.any()
    assert np.allclose((g["o2"] + g["o2_used"])[ok], g["cap"], rtol=1e-9)
    assert np.all(small["O2"][:, 0] < 1.0)
    big = L.solve_arp_pot(_params(data, "R1", L.Form(headspace_ml=1e6), n=4), pot, sig, cu)
    assert np.allclose(big["PA"], unl["PA"], rtol=1e-3)
    assert np.all(small["PA"][:, 0] <= unl["PA"][:, 0] * (1 + 1e-9))


# ----------------------------------------------------------------------------------------------- V2 and V1

@pytest.mark.parametrize("rival", ["R1", "R2", "R4"])
def test_v2_phe_and_sugar_balance_in_every_arp_pot(data, rival):
    p = _params(data, rival, L.Form(headspace_ml=20.0, nu_o2=0.5), n=4)
    for pot in list(L.POT_T2.values()) + [L.POT_FIG1] + list(L.POT_TRAP.values()):
        g = L.solve_arp_pot(p, pot, np.full(4, 0.3), np.full(4, 3.0), ledger=True)["ledger"]
        phe = g["arp"] + g["phe"] + g["pa"] + g["pa_lost"] + g["paa"] + g["phe_sink"] + g["phe_bound"]
        assert np.allclose(phe, g["arp0"] + g["phe0"], rtol=1e-9, atol=0), pot.name
        sugar = g["donors"] + g["trapped"] + g["donor_sink"] + g["pa"] + g["pa_lost"] + g["paa"] + g["sugar_frag"]
        assert np.allclose(sugar, g["arp0"] - g["arp"], rtol=1e-9, atol=1e-12), pot.name


def test_v1_route_a_and_the_glucose_dicarbonyl_part_are_bit_exact_against_ex049(data):
    v1 = G.v1_ex049(data)
    assert v1["route_a_calibration_bit_exact"] and v1["glucose_pot_dicarbonyl_bit_exact"]


def test_v4_r0_pinned_reads_the_shipped_b18_constants():
    from src.kinetic_core.parameters_pyrazine import FROZEN_B18, pyrazine_ph_factor
    c = L.r0_pinned_constants()
    assert c["log10_k_go_ak"] == FROZEN_B18["log10_k_go_ak_100C"]
    assert c["log10_k_mgo_ak"] == FROZEN_B18["log10_k_mgo_ak_100C"]
    assert c["ph_factor_7"] == pyrazine_ph_factor(7.0)
    assert round(c["log10_k_go_ak"], 2) == -6.54 and round(c["log10_k_mgo_ak"], 2) == -7.53


# ------------------------------------------------------------------------------------- calibration and units

def test_calibration_matches_table_2_and_one_evaluation_is_finite(data):
    u = L.sobol_unit(8, len(L.coords("R1")), 5)
    ev = L.evaluate(u, "R1", L.PRIMARY, data)
    f = ev.feasible
    assert f.any() and ev.failed == 0
    meas_ex = sum(data.arp_t2["air"]) - sum(data.arp_t2["argon"])
    tot = ev.diag["t2_total"][f]
    assert np.allclose(tot[:, 1] - tot[:, 0], meas_ex, rtol=1e-4)
    assert np.allclose(tot[:, 2], sum(data.arp_t2["air_cu5"]), rtol=1e-4)
    assert np.all((ev.diag["c_cu5"][f] >= 1) & (ev.diag["c_cu5"][f] <= 25))
    assert ev.pred["U3_lin"].shape == (8, 40)
    units = L.build_units(data, ev.pred["U3_lin"][np.flatnonzero(f)[0]])
    for name in ("U1e", "U1l", "U4", "U5a", "U5b", "U3", "U5a_PA_shape"):
        ll = L.unit_loglik(units[name], L.points_for(name, ev)[f])
        assert np.all(np.isfinite(ll)), name


def test_shared_reference_covariance_and_the_run_offset_refusal(data):
    units = L.build_units(data)
    u1e = units["U1e"]
    cov = L.unit_covariance(u1e)
    om2 = L.OMEGA_DEX ** 2
    assert cov[0, 1] == pytest.approx(om2 + u1e.point_sd[3] ** 2)        # omega^2 (I + 11^T): the shared 120-min point
    assert cov[0, 0] == pytest.approx(2 * om2 + u1e.point_sd[0] ** 2 + u1e.point_sd[3] ** 2)
    with pytest.raises(MeasurementModelRefusal):
        L.unit_covariance(L.replace(u1e, run_offset_dex=0.075))
    u5 = units["U5a"]
    c5 = L.unit_covariance(u5, L.Form(run_between_points=True))
    assert c5[0, 1] == pytest.approx(L.SIGMA_RUN_DEX ** 2)
    assert c5[0, 0] - L.unit_covariance(u5)[0, 0] == pytest.approx(L.SIGMA_RUN_DEX ** 2)


def test_u3_unit_censors_at_the_axis_resolution_and_the_censored_density_is_a_probability():
    y = np.full(40, 0.05)
    y[10:15] = 0.001
    u = L.u3_unit(y)
    assert u.censored.sum() == 5
    f = np.log10(np.full((2, 40), 0.05))
    f[1, 10:15] = np.log10(0.0001)                 # a prediction far below the floor is favoured on those points
    ll = L.unit_loglik(u, f)
    assert ll[1] > ll[0]


def test_u3a_and_u3b_units_follow_the_amendment():
    y = np.linspace(0.01, 0.6, 40)
    lay = L.u3_layout()
    ua = L.u3a_unit(y, 0.06, 0.03)
    dg = [k for k, (a, s, _) in enumerate(lay) if a == "argon" and s == "DG1"]
    assert np.allclose(ua.y, np.log10(y[dg]) - (math.log10(0.06) - 2.0))
    assert ua.value_basis == "contrast" and ua.design.shape == (5, 6)
    cov = L.unit_covariance(ua)
    om2 = L.OMEGA_DEX ** 2
    # off-diagonal: the shared PA point (omega^2 + its sd^2), both pots' run offsets, the 1-DG calibration
    assert cov[0, 1] == pytest.approx(om2 + ua.point_sd[5] ** 2 + 2 * L.SIGMA_RUN_DEX ** 2 + 0.20 ** 2)
    assert L.unit_covariance(ua, L.Form(dg1_cal_sd=0.30))[0, 1] - cov[0, 1] == pytest.approx(0.30 ** 2 - 0.20 ** 2)
    ub = L.u3b_unit(y)
    g = [k for k, (a, s, _) in enumerate(lay) if a == "air_cu25" and s == "G"]
    ph = [k for k, (a, s, _) in enumerate(lay) if a == "air_cu25" and s == "PHE"]
    assert np.allclose(ub.y, np.log10(y[g]) - np.log10(y[ph]))
    cb = L.unit_covariance(ub)
    assert cb[0, 1] == pytest.approx(0.10 ** 2 + 0.05 ** 2)           # the two analytes' calibrations; no run offset
    pts = np.log10(np.tile(y, (2, 1)))
    u4 = np.zeros((2, 20))
    u4[:, 0] = math.log10(0.06)
    assert np.allclose(L.unit_loglik(ua, L.u3a_points(pts, u4)), L.unit_loglik(ua, L.u3a_points(pts, u4))[0])
    yc = y.copy()
    yc[g[0]] = 0.001
    assert L.u3b_unit(yc).censored.sum() == 1


def test_smc_recovers_a_32_dim_gaussian_where_reweighting_fails():
    out = G.v_smc_gaussian_32(np.random.default_rng(3), n=512, k=4)
    assert not out["reweighting_8192"]["resolved"]
    assert out["smc"]["resolved"] and out["pass"], out


def test_uarp_likelihood_rises_with_the_allowed_decay(data):
    y = L.arp_interval_measured(data)
    c_lo = np.zeros((1, 4)) + 0.0
    c_hi = np.zeros((1, 4)) + 2.0
    lo = L.uarp_loglik(data, c_lo, np.array([0.05]), measured=y)
    hi = L.uarp_loglik(data, c_hi, np.array([0.05]), measured=y)
    assert hi[0] > lo[0] and hi[0] <= 0.0
    # the 30 -> 60 min fall needs k >= ln(48 / 32) / 30 at face value
    assert y[0, 1] == pytest.approx(math.log(4.8 / 3.2))


def test_nulls_are_finite_on_every_unit_layout(data):
    y = np.linspace(0.01, 0.6, 40)
    for unit in ("U1e", "U4", "U5a", "U5b", "U5a_PA_shape", "U3"):
        pts = L.null_points(unit, data, y)
        assert np.all(np.isfinite(pts)), unit


# ------------------------------------------------------------------------------------------------- refusals

def test_the_keyed_trapped_table_is_never_read_by_the_layer():
    with pytest.raises(PermissionError):
        L._read_table(L.KEYED_TRAPPED_TABLE)


def test_run_once_memory_guard_and_compute_cap_refuse(tmp_path, monkeypatch):
    rec = tmp_path / "x.json"
    rec.write_text("{}")
    with pytest.raises(G.Refusal):
        G.run_once(rec)
    monkeypatch.setattr(G, "memory_free_pct", lambda: 25)
    with pytest.raises(G.Refusal):
        G.memory_guard(True)
    assert G.memory_guard(False)["free_pct_at_launch"] == 25
    b = G.Budget(hours=1e-9)
    b.t0 -= 1.0
    with pytest.raises(G.Refusal):
        b.check("test")


def test_fit_and_score_refuse_while_the_prereg_is_a_draft(tmp_path):
    assert G.freeze_commit("## 20. Results\n\n*(the freeze commit is entered here first)*") is None
    with pytest.raises(G.Refusal):
        G.freeze_check()
    assert G.main(["--fit", "--out-dir", str(tmp_path), "--no-memory-guard"]) == 2
    assert G.main(["--score", "--out-dir", str(tmp_path), "--no-memory-guard"]) == 2


def _git(cwd, *args):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


def _repo(tmp_path):
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "config", "user.email", "t@t")
    _git(tmp_path, "config", "user.name", "t")
    (tmp_path / "prereg.md").write_text("design\n## 20. Results\n")
    _git(tmp_path, "add", "prereg.md")
    _git(tmp_path, "commit", "-qm", "design")
    return _git(tmp_path, "rev-parse", "HEAD")


def test_v_blind_passes_only_when_the_freeze_precedes_the_keyed_commit_and_the_hash_holds(tmp_path):
    freeze = _repo(tmp_path)
    keyed = "keyed.yml"
    with pytest.raises(G.Refusal, match="not committed"):
        G.v_blind(freeze, keyed, cwd=tmp_path)
    (tmp_path / keyed).write_text("rows: [1]\n")
    _git(tmp_path, "add", keyed)
    _git(tmp_path, "commit", "-qm", "keyed")
    ok = G.v_blind(freeze, keyed, cwd=tmp_path)
    assert ok["passed"] and ok["keyed_commit"] != freeze
    (tmp_path / keyed).write_text("rows: [2]\n")                       # (iii): an uncommitted edit
    with pytest.raises(G.Refusal, match=r"\(iii\)"):
        G.v_blind(freeze, keyed, cwd=tmp_path)
    _git(tmp_path, "commit", "-qam", "edit keyed")                     # (iii): a later commit changed it
    with pytest.raises(G.Refusal, match=r"\(iii\)"):
        G.v_blind(freeze, keyed, cwd=tmp_path)


def test_v_blind_refuses_a_keyed_commit_that_precedes_the_freeze(tmp_path):
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "config", "user.email", "t@t")
    _git(tmp_path, "config", "user.name", "t")
    (tmp_path / "keyed.yml").write_text("rows: [1]\n")
    _git(tmp_path, "add", "keyed.yml")
    _git(tmp_path, "commit", "-qm", "keyed first")
    (tmp_path / "prereg.md").write_text("design\n")
    _git(tmp_path, "add", "prereg.md")
    _git(tmp_path, "commit", "-qm", "freeze after")
    freeze = _git(tmp_path, "rev-parse", "HEAD")
    with pytest.raises(G.Refusal, match=r"\(i\)"):
        G.v_blind(freeze, "keyed.yml", cwd=tmp_path)


def test_freeze_check_on_a_scratch_repo(tmp_path):
    first = _repo(tmp_path)
    (tmp_path / "prereg.md").write_text(f"design\n## 20. Results\nFreeze commit: {first}\n")
    _git(tmp_path, "commit", "-qam", "enter the freeze commit")
    assert G.freeze_check(tmp_path, "prereg.md")["freeze_commit"] == first
    (tmp_path / "prereg.md").write_text(f"design CHANGED\n## 20. Results\nFreeze commit: {first}\n")
    with pytest.raises(G.Refusal, match="changed"):
        G.freeze_check(tmp_path, "prereg.md")


# ---------------------------------------------------------------------------------- decisions and the fit cell

def test_binomial_decision_and_the_outcome_rows():
    d = G.binomial_decision(120, 20, 200)
    assert d["decidable"] and not d["low_power"] and d["power_ci95"][0] < 0.6 < d["power_ci95"][1]
    assert not G.binomial_decision(20, 18, 200)["decidable"]
    assert G.outcome_row("R1", "R1", "none", True)[0] == "1"
    assert G.outcome_row("R1", "R1", "R4", True)[0] == "1b"
    assert G.outcome_row("R1", "R3", "none", True)[0] == "2"
    assert G.outcome_row("R1", "none", "none", True)[0] == "3"
    assert G.outcome_row("R2", "none", "none", True)[0] == "4"
    assert G.outcome_row("none", "R1", "none", True)[0] == "5"
    assert G.contrast_credit({"reading": "a"}, {"U1e": {"reading": "tie"}}) == "a"
    assert G.contrast_credit({"reading": "a"}, {"U1e": {"reading": "b"}}) == "none"


def test_a_rival_with_no_feasible_draw_reads_one_sided(data):
    empty = G.DrawSet("R2", L.PRIMARY, 0, 100, np.zeros((0, 3)), {}, None, None, None, None, {}, {}, 0, 0.0, {})
    full = G.DrawSet("R1", L.PRIMARY, 0, 1, np.zeros((1, 3)), {}, None, None, None, None, {}, {}, 0, 0.0, {})
    v = G.unit_verdict(full, empty, "U3", (), G.Obs({}, np.zeros((1, 4))), data, L.PRIMARY,
                       np.random.default_rng(0), 2)
    assert v["reading"] == "a" and v["by_infeasibility"] and v["d"] == math.inf


def test_score_contrast_needs_agreement_adequacy_and_the_primary():
    d = {"U3": {"primary": 5.0, "cu": 4.0}, "U1e": {"primary": 0.1, "cu": 0.0}}
    b = {"U3": {"primary": 1.0, "cu": 1.0}, "U1e": {"primary": 1.0, "cu": 1.0}}
    ok = G.score_contrast(d, b, {"U3": {"a": True}}, "U3", {})
    assert ok["credit"] == "a"
    assert G.score_contrast(d, b, {"U3": {"a": False}}, "U3", {})["credit"] == "none"
    d["U3"]["cu"] = 0.5
    assert G.score_contrast(d, b, {"U3": {"a": True}}, "U3", {})["credit"] == "none"


def test_fit_cell_recovers_an_analytic_gaussian_lppd():
    """The decisive fits' cell (mode search + Student-t PSIS) on a toy: L = N(y; theta, s^2 I) on the unit square,
    y in the middle and s small, so lppd = log of the Gaussian's mass inside the square (~ 0)."""
    rng = np.random.default_rng(0)
    y, s = np.array([0.4, 0.6]), 0.05

    def ev(pts):
        r = (pts - y) / s
        ll = -0.5 * np.sum(r ** 2, axis=1) - 2 * math.log(s) - math.log(2 * math.pi)
        return ll, r

    out = G.fit_cell(ev, 2, [np.full(2, 0.5), np.array([0.2, 0.8]), np.array([0.7, 0.3])], rng, evals=60)
    assert np.allclose(out["theta_mode"], y, atol=1e-4)
    assert out["estimate"].accepted
    est = out["estimate"]
    assert abs(est.lppd) <= 3.0 * est.mc_se + 0.01, est
