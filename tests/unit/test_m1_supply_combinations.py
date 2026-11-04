"""Plan step M1's build (scripts/generators/fit_m1_supply_combinations.py; prereg
results/validation/kinetic_core_m1_supply_loss_combinations_prereg.md v3.3): the decision rule with its tie band, the
comparison rule, the assay floors, the encounter ceiling, the refusals, and the shipped Config unchanged."""
from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

from src.kinetic_core import config, engine  # noqa: E402
from src.kinetic_core import parameters_sulfur as ps  # noqa: E402

NEW_KEYS = ("fit.dimer_encounter_ceiling", "fit.mft_decay_prior")


@pytest.fixture(scope="module")
def M():
    """The generator, imported through tests.support.wave_generator: it pulls in B60 -> B47 -> ... -> B2.3, which
    splice B2.3's row table and wrap its integrator; imported bare, they stay spliced for every later test (B10, B16,
    B2.4's tests failed so, 2026-11-04). Dropped again on exit."""
    from tests.support import wave_generator
    with wave_generator("fit_m1_supply_combinations") as module:
        yield module
    for name in ("fit_m1_supply_combinations", "m1_preflight_checks"):
        sys.modules.pop(name, None)


# --- the decision rule (sec. 7) ------------------------------------------------------------------------------------

def test_block_se_with_two_units_is_the_absolute_difference(M):
    for d1, d2 in ((1.0, 3.0), (-2.5, 4.0), (0.3, 0.3), (-7.0, -1.0)):
        assert M.block_se([d1, d2]) == pytest.approx(abs(d1 - d2))
    assert M.block_se([1.0]) == float("inf")


@pytest.mark.parametrize("d1, d2", [(1.0, 3.0), (-1.0, 3.0), (-1.0, -3.0), (2.0, -0.5), (0.1, 5.0)])
def test_the_one_se_rule_on_two_units_is_a_sign_test(M, d1, d2):
    dl, se = d1 + d2, M.block_se([d1, d2])
    assert (dl > se) == (d1 > 0 and d2 > 0)
    assert (dl < -se) == (d1 < 0 and d2 < 0)


def test_unit_side_needs_one_sign_under_every_density_and_sensitivity(M):
    assert M.unit_side([1.0, 0.5, 2.0, 0.1, 3.0]) == "C1"
    assert M.unit_side([-1.0, -0.5, -2.0, -0.1, -3.0]) == "R0"
    assert M.unit_side([1.0, 0.5, -0.01, 0.1, 3.0]) == "tie"     # one sensitivity flips: a tie (v3.3)
    assert M.unit_side([1.0, 0.0, 1.0, 1.0, 1.0]) == "tie"       # an exact zero is no side


def test_the_three_outcomes(M):
    assert M.block_outcome(["C1", "C1"], ["C1", "R0"]) == M.C1_FAVOURED
    assert M.block_outcome(["C1", "C1"], ["R0", "R0"]) == M.NOT_SEPARATED     # C1 worse on the thiol shape
    assert M.block_outcome(["R0", "R0"], ["C1", "C1"]) == M.R0_FAVOURED      # R0 favoured whatever D1 says
    assert M.block_outcome(["C1", "R0"], ["C1", "C1"]) == M.NOT_SEPARATED     # conflict: not evidence against C1


def test_the_tie_band_a_reservoir_tie_is_not_separated_and_a_d1_tie_is_not_worse(M):
    assert M.block_outcome(["C1", "tie"], ["C1", "C1"]) == M.NOT_SEPARATED
    assert M.block_outcome(["R0", "tie"], ["R0", "R0"]) == M.NOT_SEPARATED
    assert M.block_outcome(["C1", "C1"], ["R0", "tie"]) == M.C1_FAVOURED
    assert M.block_outcome(["C1", "C1"], ["tie", "tie"]) == M.C1_FAVOURED


def _per_form(c1, r0):
    """{form: {setting: {unit: {linearised, plug_in, block}}}} from per-unit lppd tuples (lin, plug, floors, t05, t2)."""
    out = {}
    for name, vals in (("C1", c1), ("R0", r0)):
        out[name] = {}
        for k, setting in enumerate(("primary", "floors_alternative", "tau_x0.5", "tau_x2")):
            out[name][setting] = {}
            for u, v in vals.items():
                lin = v[0] if setting == "primary" else v[1 + k]
                out[name][setting][u] = {"linearised": lin, "plug_in": v[1] if setting == "primary" else lin,
                                         "block": "D1" if u.startswith("D1") else "reservoir"}
    return out


def test_decide_reads_the_sign_pattern_with_the_tie_band(M):
    base = {"D1_MFT": (0, 0, 0, 0, 0, 0), "D1_FFT": (0, 0, 0, 0, 0, 0)}
    c1 = {**base, "RES_Cys": (5, 5, 5, 5, 5, 5), "RES_PENT": (3, 3, 3, 3, 3, 3)}
    r0 = {**base, "RES_Cys": (1, 1, 1, 1, 1, 1), "RES_PENT": (1, 1, 1, 1, 1, 1)}
    rec = M.decide(_per_form(c1, r0))
    assert rec["outcome"] == M.C1_FAVOURED                # D1 exactly tied: ties count as not worse
    assert rec["blocks"]["reservoir"]["SE"] == pytest.approx(abs(4 - 2))
    # the plug-in density reverses one reservoir unit: that unit is a tie -> not separated
    c1b = dict(c1, RES_PENT=(3, 0, 3, 3, 3, 3))
    assert M.decide(_per_form(c1b, r0))["units"]["RES_PENT"]["side"] == "tie"
    assert M.decide(_per_form(c1b, r0))["outcome"] == M.NOT_SEPARATED


def test_the_verdict_counts_only_if_every_reading_agrees(M):
    same = {r: M.C1_FAVOURED for r in M.READINGS}
    assert M.verdict_across_readings(same) == {"verdict": M.C1_FAVOURED, "agree": True, "disagreeing_readings": []}
    mixed = dict(same, ramp_plus_hold=M.R0_FAVOURED)
    v = M.verdict_across_readings(mixed)
    assert v["verdict"] == M.NOT_SEPARATED and v["disagreeing_readings"] == ["ramp_plus_hold"]


def test_the_verdict_text_carries_the_cf063_note_and_the_credit(M):
    text = M.verdict_text(M.C1_FAVOURED, gigl=True)
    assert "CF-063" in text and "Zhang 2025b alone" in text and "not attributed" in text and "Gigl ceiling" in text
    assert "not evidence against C1" in M.verdict_text(M.NOT_SEPARATED, gigl=False)


# --- the comparison rule (sec. 4.3) --------------------------------------------------------------------------------

def _starts(*rows):
    return [{"start": s, "status": st, "cost": c} for s, st, c in rows]


def test_comparison_rule_a_is_two_converged_starts_within_1_92_nats(M):
    assert M.comparison_rule(_starts((0, 2, 10.0), (1, 1, 11.91), (2, 0, 9.0)))["a"] is True
    r = M.comparison_rule(_starts((0, 2, 10.0), (1, 1, 11.93), (2, 0, 9.0)))
    assert r["a"] is False and r["status"] == "not_compared" and r["best_start"] == 0   # status 0 never counts
    assert M.comparison_rule(_starts((0, 0, 1.0), (1, -9, 1.0)))["best_start"] is None


def test_comparison_rule_b_is_prediction_agreement_without_the_data(M):
    starts = _starts((0, 2, 10.0), (1, 2, 10.5), (2, 1, 30.0))
    tol = np.array([0.06, 0.03])
    agree = {0: np.array([1.0, 2.0]), 1: np.array([1.05, 2.02])}
    assert M.comparison_rule(starts, agree, tol)["status"] == "compared"
    disagree = {0: np.array([1.0, 2.0]), 1: np.array([1.0, 2.05])}
    r = M.comparison_rule(starts, disagree, tol)
    assert r["b"] is False and r["status"] == "compared_if_verdicts_agree"
    assert r["qualifying_starts"] == [0, 1]


# --- the assay floors (sec. 5) -------------------------------------------------------------------------------------

def test_a_prediction_below_the_floor_is_scored_at_the_floor_and_a_failed_run_is_not(M):
    assert M.floored(0.0, 1.22) == (1.22, True)
    assert M.floored(5.0, 1.22) == (5.0, False)
    v, fl = M.floored(float("nan"), 1.22)
    assert math.isnan(v) and fl is False


def test_units_floor_the_reservoir_and_zero_its_jacobian_rows(M):
    labels = M.point_labels()
    levels = np.zeros(len(labels))
    for i, lbl in enumerate(labels):
        block, s, t = lbl.split("|")
        levels[i] = math.log10(1000.0) if block == "D1" else (-30.0 if s == "PENT" else math.log10(10.0))
    jac = np.ones((len(labels), 2))
    us = M.units(levels, jac, M.Scoring(tau=0.2))
    pent = us["RES_PENT"]
    assert pent["n_floored"] == 8 and np.allclose(pent["f"], math.log10(M.RES_FLOOR_MMOL_L["PENT"]))
    assert np.all(pent["J"] == 0.0) and np.all(us["RES_Cys"]["J"] == 1.0)
    # the reservoir covariance carries the integrated level offset tau^2 11'
    cov = us["RES_Cys"]["Sigma"]
    assert cov[0, 1] == pytest.approx(0.2 ** 2)
    # D1: a bivariate unit, the shared 360 min point gives the ratios a negative covariance
    d1 = us["D1_MFT"]["Sigma"]
    assert d1[0, 1] < 0 and math.sqrt(d1[0, 0]) == pytest.approx(0.059, abs=0.002)
    alt = M.units(levels, jac, M.Scoring(tau=0.2, floors="alternative"))["RES_PENT"]
    assert np.allclose(alt["f"], math.log10(M.RES_FLOOR_ALT_MMOL_L["PENT"]))


# --- the encounter ceiling (sec. 2.5) ------------------------------------------------------------------------------

def test_the_encounter_ceiling_is_log10_9_47_at_145_c(M):
    c = ps.DimerEncounterCeiling(enabled=True)
    assert c.water_viscosity_pa_s(418.15) == pytest.approx(0.19e-3, rel=0.02)
    assert c.encounter_limit_l_per_mmol_min(418.15) == pytest.approx(2.9e9, rel=0.03)
    assert c.log10_ceiling_at(145.0) == pytest.approx(9.47, abs=0.005)
    # a cooler pot allows a higher 145 C constant (the Arrhenius factor of 122.2 kJ/mol); a hotter one a lower
    assert c.log10_ceiling_at(100.0) > c.log10_ceiling_at(145.0) > c.log10_ceiling_at(160.0)
    assert c.log10_ceiling([70.0, 100.0, 145.0]) == pytest.approx(c.log10_ceiling_at(145.0))
    assert dataclasses.replace(c, mft_cap_mmol_l=10.0).log10_ceiling_at(145.0) == pytest.approx(8.47, abs=0.005)


def test_the_mft_cap_rule(M):
    assert M.mft_cap_from_panel_max(0.05) == 1.0
    assert M.mft_cap_from_panel_max(0.1) == 1.0
    assert M.mft_cap_from_panel_max(0.2) == pytest.approx(2.0)


def test_the_prior_is_centred_on_b60s_lambda(M):
    p = ps.MftDecayPrior(enabled=True)
    assert p.centre_log10_per_min == pytest.approx(-1.391, abs=5e-4)
    assert abs(p.residual(M.prior_log10_rate(M.b60_point().x(), p))) < 1e-3
    assert M.lambda_max_log10(101.25) == pytest.approx(math.log10(0.094), abs=0.01)


# --- the refusals --------------------------------------------------------------------------------------------------

def _ok_freeze(commit):
    return "f" * 40


def _bad_freeze(commit):
    raise SystemExit("the prereg changed since the freeze commit")


def _records(status=2, failed=0):
    return {n: [{"start": 0, "status": status, "failed_pot_runs": failed}] for n in ("R0", "C1")}


def test_score_refuses_on_each_missing_condition(M):
    carried = {"R0": {"status": "compared"}, "C1": {"status": "compared"}}
    assert M.score_refusals("x", 0.3, _records(), carried, _ok_freeze) == []
    assert any("prereg changed" in p for p in M.score_refusals("x", 0.3, _records(), carried, _bad_freeze))
    assert any("status > 0" in p for p in M.score_refusals("x", 0.3, _records(status=0), carried, _ok_freeze))
    assert any("failed_pot_runs" in p for p in M.score_refusals("x", 0.3, _records(failed=3), carried, _ok_freeze))
    assert any("tau_ABS" in p for p in M.score_refusals("x", None, _records(), carried, _ok_freeze))
    assert any("finite and > 0" in p for p in M.score_refusals("x", 0.0, _records(), carried, _ok_freeze))
    assert any("no fit record" in p for p in M.score_refusals("x", 0.3, {"R0": []}, carried, _ok_freeze))
    assert any("carried" in p for p in M.score_refusals("x", 0.3, _records(), {}, _ok_freeze))


def test_fit_refuses_without_the_freeze_commit_and_without_passing_checks(M):
    with pytest.raises(SystemExit, match="freeze-commit"):
        M.check_freeze("")
    with pytest.raises(SystemExit, match="--freeze-commit"):
        M.main(["--fit", "R0", "--start", "0"])
    good = {"checks": {v: {"passed": True, "integrated_without_failure": True} for v in ("V1", "V2", "V5", "V6", "V7", "V8")}}
    assert M.fit_refusals("x", good, _ok_freeze) == []
    v7_failed = json.loads(json.dumps(good))
    v7_failed["checks"]["V7"]["passed"] = False
    assert any("V7 failed" in p for p in M.fit_refusals("x", v7_failed, _ok_freeze))
    assert any("has not been run" in p for p in M.fit_refusals("x", {"checks": {}}, _ok_freeze))
    v8_failed = json.loads(json.dumps(good))
    v8_failed["checks"]["V8"]["passed"] = False
    assert M.numerics_for_fits(v8_failed).coupling == "fixed_point_nodes"
    assert M.numerics_for_fits(good).coupling == "inline"


def test_check_freeze_reads_the_prereg_against_the_commit(M, tmp_path):
    with pytest.raises(SystemExit, match="not a commit"):
        M.check_freeze("0000000000000000000000000000000000000000")


def test_read_tau_takes_the_command_line_then_the_h0_record(M, tmp_path):
    assert M.read_tau(0.4, tmp_path / "absent.json") == 0.4
    assert M.read_tau(None, tmp_path / "absent.json") is None
    p = tmp_path / "tau.json"
    p.write_text(json.dumps({"D2": {"tau_ABS_all_used": 0.31}}))
    assert M.read_tau(None, p) == pytest.approx(0.31)


def test_the_memory_plan(M):
    assert M.memory_plan(2 * 10 ** 9, 5 * 10 ** 9) == "parallel"
    assert M.memory_plan(3 * 10 ** 9, 5 * 10 ** 9) == "sequential"
    assert M.memory_plan(6 * 10 ** 9, 5 * 10 ** 9) == "refuse"
    assert M.memory_plan(0, 5 * 10 ** 9) == "refuse"


# --- the shipped engine unchanged ------------------------------------------------------------------------------------

def test_the_shipped_config_is_unchanged_with_every_new_key_at_its_default(M):
    shipped = config.Config.shipped()
    assert config.Config.from_globals() == shipped
    for name in NEW_KEYS:
        assert shipped[name].enabled is False
    # the digest over the keys that existed before M1 is M1.0's (EX-060: configuration digest 40b746476905786c)
    blob = json.dumps([[n, config._canonical(v)] for n, v in shipped.items if n not in NEW_KEYS], sort_keys=True,
                      separators=(",", ":"))
    assert hashlib.sha256(blob.encode("utf-8")).hexdigest().startswith("40b746476905786c")
    route = shipped["sulfur.cys_amadori_route"]
    assert route.enabled is False and route.ea_form_kj_mol == 100.9
    assert config.KEY_BY_NAME["sulfur.cys_amadori_route"].bands["ea_form_kj_mol"] == (40.0, 135.0)


def test_the_kinetics_never_read_the_fit_side_keys(M):
    on = config.Config.shipped().replace({
        "fit.dimer_encounter_ceiling": ps.DimerEncounterCeiling(enabled=True, mft_cap_mmol_l=3.0),
        "fit.mft_decay_prior": ps.MftDecayPrior(enabled=True, sd_dex=0.1)})
    params = engine.core_parameters(engine.SULFUR)

    def run():
        from src.kinetic_core.sulfur import integrate_sulfur
        r = integrate_sulfur(params, 120.0 + 273.15, {"PENT": 20.0, "Cys": 20.0, "OX": 1.0}, [0.0, 5.0], ph=6.0,
                             rtol=1e-6, atol=1e-12)
        return [float(r.final(k)) for k in ("MFT", "FFT", "Cys", "PENT")]

    shipped_run = run()
    with config.applied(on):
        flipped_run = run()
    assert flipped_run == shipped_run


def test_the_rivals_coordinates_are_the_prereg_s(M):
    assert len(M.RIVALS["R0"].keys) == 17 and len(M.RIVALS["C1"].keys) == 19
    assert "k_mft_decay" in M.RIVALS["R0"].keys and "k_ttca_cys" not in M.RIVALS["R0"].keys
    assert not set(M.HELD_B60_UNIDENTIFIED) & set(M.RIVALS["C1"].keys)
    rows, jobs, _ = M.route_block()
    assert len(rows) == 296 - 6 and len(jobs) == 34 - 1      # the Kang duplicate pot dropped (sec. 4.1's check)
    assert not {r.pot for r in rows} & set(M.ROUTE_OUT_POTS)


# --- amendment A1: deterministic exploratory pot failures ------------------------------------------------------------

def _fit_rec(cost=10.0, n=1, seg=1, evals=100, pots=("b60:kang_freecys_100",)):
    r = np.array([3.0, 3.0, 2.0])   # 0.5 * r.r = 11.0
    scale = math.sqrt(cost / 11.0)
    return {"artifact": "kinetic_core_m1_supply_combinations_fit_C1_s1", "cost": cost, "evaluations": evals,
            "failed_pot_runs": n, "history": [{"failed_pot_runs_segment": seg}],
            "residual_by_row": {f"r{i}": float(v * scale) for i, v in enumerate(r)},
            "failures_sample": [{"pot": p, "reason": "exceeded 120 s", "timeout": True} for p in pots]}


def test_a_reproduced_exploratory_failure_with_a_clean_optimum_is_accepted(M, tmp_path, monkeypatch):
    monkeypatch.setattr(M.data_paths, "VALIDATION_DIR", tmp_path)
    rec = _fit_rec()
    (tmp_path / f"{rec['artifact']}_failed_run1.json").write_text(json.dumps(rec))
    ok, why = M.failure_accepted(rec)
    assert ok and "amendment A1" in why


def test_a_failure_is_refused_without_reproduction_or_at_the_optimum(M, tmp_path, monkeypatch):
    monkeypatch.setattr(M.data_paths, "VALIDATION_DIR", tmp_path)
    rec = _fit_rec()
    assert not M.failure_accepted(rec)[0]                                   # no same-settings re-run on disk
    twin = dict(rec, cost=10.5)
    (tmp_path / f"{rec['artifact']}_failed_run1.json").write_text(json.dumps(twin))
    assert not M.failure_accepted(rec)[0]                                   # the re-run differs: machine state
    (tmp_path / f"{rec['artifact']}_failed_run1.json").write_text(json.dumps(rec))
    assert not M.failure_accepted(_fit_rec(n=2, seg=1))[0]                  # a failure at the optimum itself
    bad = _fit_rec()
    bad["cost"] = 12.0                                                       # residuals do not give the cost
    assert not M.failure_accepted(bad)[0]
    assert M.failure_accepted(_fit_rec(n=0, seg=0))[0]                      # no failure: nothing to accept
