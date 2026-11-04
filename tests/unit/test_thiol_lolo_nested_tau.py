"""H0's scoring step (results/validation/kinetic_core_between_lab_tau_h0_prereg.md): the estimator's four identity
tests of sec. 2, the vectorised predictive against fitting.MeasurementModel, the St quadrature against its Gaussian
limit, and the blinding precondition (--score refuses to run without the ten pair records, each with its freeze
commit). Synthetic rows only."""
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "generators"))

import fit_thiol_lolo as lolo  # noqa: E402
import fit_thiol_lolo_pairs as pair_runner  # noqa: E402
import score_thiol_lolo_nested_tau as h0  # noqa: E402
from src.kinetic_core import fitting  # noqa: E402


def _row(paper, species, method, i, lab="lab_a", sigma=0.05, table="thiol"):
    return SimpleNamespace(ledger_id=f"{paper}_t{i % 2}", species=species, method=method, table=table, t_min=float(i),
                           sigma=sigma, series=f"{paper}|{species}|{method}", laboratory=lab, pot=f"{paper}_p{i}",
                           value=1.0, unit="mmol/L")


def _abs_rows(lab, n_groups, rows_per_group, method="sida", tag=""):
    rows = []
    for g in range(n_groups):
        for k in range(rows_per_group if not isinstance(rows_per_group, list) else rows_per_group[g]):
            rows.append(_row(f"p{tag}{lab[-1]}{g}{method[0]}", f"S{g}", method, k, lab=lab,
                             sigma=0.03 + 0.01 * ((g + k) % 4)))
    return rows


def _model_and_s2(rows, noise=0.1):
    model = lolo._model({fitting_series(r): noise for r in rows}, level=lolo.level_of)
    s2 = np.array([model.sigma(r) ** 2 for r in rows])
    return model, s2


def fitting_series(r):
    return lolo.fj.series_of(r)


def _residuals(rows, rng, tau=0.7):
    lay = h0.Layout.of_rows(rows)
    level = tau * rng.standard_normal(lay.G)
    return level[lay.gid] + 0.1 * rng.standard_normal(len(rows))


# ---- sec. 2 identity tests -------------------------------------------------------------------------------------

def test_identity_1_one_lab_one_method_equals_the_shipped_estimator():
    rng = np.random.default_rng(1)
    rows = _abs_rows("lab_a", 6, [1, 2, 3, 4, 1, 5])
    model, s2 = _model_and_s2(rows)
    r = _residuals(rows, rng)
    shipped = fitting.estimate_between_tau(model, r, rows)["sida"]
    est = h0.estimate([h0.Block("lab_a", h0.Layout.of_rows(rows), s2, r[None, :])], 1, omega2_fixed=0.0)
    assert shipped > 0.0
    assert abs(float(est["ABS"]["tau"][0]) - shipped) <= 1e-12 * max(1.0, shipped)
    assert est["ABS"]["status"] == "single-lab"   # one laboratory: flagged, the estimate itself unchanged


def test_identity_2_pooled_dbar_is_the_group_count_weighted_mean_of_the_methods():
    rng = np.random.default_rng(2)
    rows = _abs_rows("lab_a", 4, [2, 3, 1, 4], "sida") + _abs_rows("lab_a", 3, [3, 1, 2], "calibrated", tag="c")
    _, s2 = _model_and_s2(rows)
    r = _residuals(rows, rng)
    lay = h0.Layout.of_rows(rows)
    dbar = {}
    for method in ("sida", "calibrated"):
        d = []
        for g in range(lay.G):
            if lay.method_g[g] != method:
                continue
            idx = lay.gid == g
            w = 1.0 / s2[idx]
            m = float(np.sum(w * r[idx]) / np.sum(w))
            d.append(m * m - 1.0 / float(np.sum(w)))
        dbar[method] = (len(d), float(np.mean(d)))
    expected = sum(n * v for n, v in dbar.values()) / sum(n for n, _ in dbar.values())
    est = h0.estimate([h0.Block("lab_a", lay, s2, r[None, :])], 1, omega2_fixed=0.0)
    assert abs(float(est["ABS"]["tau2_raw"][0]) - expected) <= 1e-12 * max(1.0, abs(expected))


def test_identity_3_within_group_moment_has_expectation_c_omega2():
    rng = np.random.default_rng(3)
    rows = _abs_rows("lab_a", 3, [2, 3, 5])
    lay = h0.Layout.of_rows(rows)
    s2 = np.array([0.02, 0.05, 0.01, 0.03, 0.08, 0.02, 0.04, 0.06, 0.01, 0.05])
    omega2 = 0.09
    R = 100_000
    level = 3.0 * rng.standard_normal((R, lay.G))   # an arbitrary level: it must cancel
    r = level[:, lay.gid] + np.sqrt(s2 + omega2)[None, :] * rng.standard_normal((R, lay.n))
    _, e, c, sel = h0.within_moments(lay, s2, r)
    assert sel.all()
    mean, se = e.mean(axis=0), e.std(axis=0, ddof=1) / math.sqrt(R)
    assert np.all(np.abs(mean - c * omega2) <= 3.0 * se), (mean - c * omega2) / se


def test_identity_4_equal_lab_sizes_make_lab_weights_equal_group_weights():
    rng = np.random.default_rng(4)
    blocks, pooled = [], []
    for lab in ("lab_a", "lab_b", "lab_c"):
        rows = _abs_rows(lab, 4, [1, 3, 2, 2])
        _, s2 = _model_and_s2(rows)
        r = _residuals(rows, rng, tau=0.9)
        blocks.append(h0.Block(lab, h0.Layout.of_rows(rows), s2, r[None, :]))
        pooled.append((rows, r))
    rows_all = [x for rows, _ in pooled for x in rows]
    r_all = np.concatenate([r for _, r in pooled])
    model, _ = _model_and_s2(rows_all)
    shipped = fitting.estimate_between_tau(model, r_all, rows_all)["sida"]
    est = h0.estimate(blocks, 1, omega2_fixed=0.0)
    assert est["ABS"]["status"] == "normal" and shipped > 0.0
    assert abs(float(est["ABS"]["tau"][0]) - shipped) <= 1e-12 * max(1.0, shipped)


# ---- the vectorised predictive equals MeasurementModel; St reduces to the Gaussian ---------------------------

def test_vectorised_arm_matches_measurement_model_and_row_intervals():
    rng = np.random.default_rng(5)
    rows = (_abs_rows("lab_a", 3, [1, 3, 2])
            + [_row("figA", "DPO", "semi_quantitative", k, table="deoxyosone_submodel") for k in range(4)]
            + [_row("zhou2023", "MFT", "relative", k) for k in range(5)]
            + [_row("zhou2024", "FFT", "relative", 0)])
    _, s2 = _model_and_s2(rows)
    r = rng.standard_normal(len(rows))
    omega2, taus = 0.04, {"sida": 0.5, "calibrated": 0.5, "semi_figure": 0.45}
    arm = h0.gaussian_arm(h0.Layout.of_rows(rows), s2, r[None, :], omega2, taus)
    noise = {lolo.fj.series_of(x): math.sqrt(v - x.sigma ** 2 + omega2) for x, v in zip(rows, s2)}
    pm = lolo._model(noise, {"sida": 0.5, "calibrated": 0.5}, lolo.level_of)
    import dataclasses
    pm = dataclasses.replace(pm, scale_tau={**lolo.sm3.SCALE_TAU, "semi_figure": 0.45})
    assert abs(float(arm["lpd"][0]) - pm.log_density(r, rows)) <= 1e-9 * abs(pm.log_density(r, rows))
    for p, c, v in zip(lolo.row_intervals(pm, r, rows), arm["centred"][0], arm["var"][0]):
        if p["var"] is None:
            assert not np.isfinite(v)
        else:
            assert abs(p["var"] - v) <= 1e-12 and abs(p["centred"] - c) <= 1e-12


def test_st_reduces_to_the_gaussian_at_large_df():
    rng = np.random.default_rng(6)
    r, st2 = rng.standard_normal(4) * 0.5, np.array([0.01, 0.04, 0.02, 0.03])
    tau = 0.6
    W, B, A = np.sum(1 / st2), np.sum(r / st2), np.sum(r * r / st2)
    gauss = -0.5 * (A - tau ** 2 * B * B / (1 + tau ** 2 * W) + np.sum(np.log(st2)) + math.log1p(tau ** 2 * W)
                    + 4 * math.log(2 * math.pi))
    assert abs(h0.st_group_logdensity(r, st2, tau, 1e7) - gauss) < 1e-5
    from scipy.stats import norm
    x, sd = np.array([0.3, 1.0, 2.0]), np.array([0.2, 0.5, 1.0])
    for scale in (0.1, 0.8):
        got = h0.st_cdf_abs(x, sd, np.full(3, scale), 1e7, nodes=200)
        want = 2 * norm.cdf(x / np.sqrt(sd ** 2 + scale ** 2)) - 1
        assert np.allclose(got, want, atol=1e-6)


# ---- blinding: --score refuses without the ten pair records ------------------------------------------------

def _pair_record(a, b, commit="abc1234"):
    return {"artifact": f"{pair_runner.STEM}_{a}__{b}", "labs": [a, b], "freeze_commit": commit, "complete": True,
            "row_keys": {}, "held_engine_raw": {}, "train_engine_raw": {}}


def test_score_refuses_without_any_pair_record(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    monkeypatch.setattr(h0, "OUT_DIR", tmp_path / "out")  # the real primary record exists since 2026-11-04
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "refuses" in str(exc.value.code)


def test_score_refuses_with_nine_of_ten_pair_records(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    monkeypatch.setattr(h0, "OUT_DIR", tmp_path / "out")  # the real primary record exists since 2026-11-04
    for a, b in pair_runner.pairs()[:9]:
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b)))
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "missing" in str(exc.value.code)


def test_score_refuses_when_a_pair_record_has_no_freeze_commit(tmp_path, monkeypatch):
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    monkeypatch.setattr(h0, "OUT_DIR", tmp_path / "out")  # the real primary record exists since 2026-11-04
    ps = pair_runner.pairs()
    for i, (a, b) in enumerate(ps):
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b, commit="" if i == 3 else "abc1234")))
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "freeze_commit" in str(exc.value.code)


def test_score_refuses_a_pair_record_with_failed_pot_runs(tmp_path, monkeypatch):
    # amendment A2: a failed pot run (a wall-clock PotTimeout) makes the record depend on machine state
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    monkeypatch.setattr(h0, "OUT_DIR", tmp_path / "out")  # the real primary record exists since 2026-11-04
    for i, (a, b) in enumerate(pair_runner.pairs()):
        rec = _pair_record(a, b)
        rec["failed_pot_runs"] = 55 if i == 0 else 0
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(rec))
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert exc.value.code not in (0, None) and "failed pot runs" in str(exc.value.code)


def test_the_ten_complete_pair_records_pass_the_precondition(tmp_path, monkeypatch):
    for a, b in pair_runner.pairs():
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(_pair_record(a, b)))
    assert len(h0.check_pair_records(tmp_path)) == 10


def test_the_sign_flip_and_the_exact_bootstrap_are_exact():
    counts, w = h0.multisets(5)
    assert counts.shape == (126, 5) and abs(w.sum() - 1.0) < 1e-12
    assert h0.sign_patterns(5).shape == (32, 5)
    d = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    res = h0.d3_decision(d[None, :], -d[None, :])
    assert res["p_better"][0] == 1 / 32 and res["lpd_pass"][0] and not res["sign_trigger"][0]


# ---- amendment A1 (kinetic_core_thiol_lolo_nested_tau_amendment_A1.md) -----------------------------------------

@pytest.mark.parametrize("change, needle", [("complete", "not complete"), ("held_engine_raw", "missing key"),
                                            ("row_keys", "missing key"), ("train_engine_raw", "missing key")])
def test_score_refuses_an_incomplete_pair_record_or_one_without_its_residual_keys(tmp_path, change, needle):
    for i, (a, b) in enumerate(pair_runner.pairs()):
        rec = _pair_record(a, b)
        if i == 4:
            if change == "complete":
                rec["complete"] = False
            else:
                del rec[change]
        h0.pair_path(a, b, tmp_path).write_text(json.dumps(rec))
    with pytest.raises(SystemExit) as exc:
        h0.check_pair_records(tmp_path)
    assert needle in str(exc.value.code)


GUARD_LAB = h0.GUARD_RECORD_LAB


def _fake_terms(raw, train, n):
    """A deterministic stand-in for training_terms (the guard's wiring, not training_terms itself, is under test)."""
    raw = np.asarray(raw, dtype=float)
    return {"class_noise": {"sida": 0.1 + 0.01 * n}, "between_tau": {"sida": float(np.nansum(raw ** 2)) + len(train)},
            "n_training_rows_scored": int(np.isfinite(raw).sum())}


def _fake_null(kind, train, rows):
    return np.linspace(-1.0, 1.0, len(rows)) * (1.0 if kind == "typical_yield" else 2.0), 3


def _guard_world():
    data = tuple(_abs_rows(GUARD_LAB, 2, [2, 1]) + _abs_rows("lab_b", 2, [2, 2]))
    folds_rows = lolo.make_folds(data)
    train_r, held_r = folds_rows[GUARD_LAB]
    train_raw = np.linspace(0.2, 0.9, len(train_r))
    held_raw = np.array([0.1, np.nan, 0.3])
    folds = {}
    for lab, (train, held) in folds_rows.items():
        terms = {f"null_{k}": _fake_terms(_fake_null(k, train, train)[0], train, 3) for k in lolo.NULLS}
        if lab == GUARD_LAB:
            terms["engine_warm"] = _fake_terms(train_raw, train, h0.N_CONSTANTS_ENGINE)
        folds[lab] = h0.Fold(lab=lab, held=held, paired=np.ones(len(held), bool),
                             raw={"engine_warm": held_raw[None, :] if lab == GUARD_LAB else np.zeros((1, len(held)))},
                             terms=terms, leak_prone=[], record={"warm": {"fit": {"cost": 100.0}}})
    guard_rec = {"labs": [GUARD_LAB], "complete": True,
                 "row_keys": {"train": h0.row_keys(train_r), "held": {GUARD_LAB: h0.row_keys(held_r)}},
                 "train_engine_raw": {"warm": train_raw.tolist()},
                 "held_engine_raw": {"warm": {GUARD_LAB: [0.1, float("nan"), 0.305]}},
                 "warm": {"fit": {"cost": 100.05}}}
    return data, folds_rows, folds, guard_rec


def _run_guard(world):
    data, folds_rows, folds, guard_rec = world
    return h0.reproduction_guards(data, folds_rows, folds, None, guard_rec, null_fn=_fake_null, terms_fn=_fake_terms)


def test_the_reproduction_guards_pass_when_every_saved_term_is_reproduced():
    rg = _run_guard(_guard_world())
    assert rg["pass"] and rg["terms"]["n_checks"] == 2 * 2 + 1
    assert rg["terms"]["max_relative_difference"] == 0.0
    assert rg["engine_guard_4_2c"]["held_out_raw"]["same_finite_mask"]


@pytest.mark.parametrize("where", ["null", "engine"])
def test_a_perturbed_saved_term_makes_the_reproduction_guard_refuse(where):
    world = _guard_world()
    folds = world[2]
    if where == "null":
        folds["lab_b"].terms["null_no_change"]["class_noise"]["sida"] *= 1.0 + 1e-7
    else:
        folds[GUARD_LAB].terms["engine_warm"]["between_tau"]["sida"] *= 1.0 + 1e-7
    rg = _run_guard(world)
    assert not rg["pass"] and not rg["terms"]["pass"]
    failed = h0._failed_guard_checks(rg)
    assert failed == (["lab_b null_no_change"] if where == "null" else [f"{GUARD_LAB} engine_warm (guard record "
                                                                         "train_engine_raw['warm'])"])


def test_engine_guard_4_2c_refuses_on_cost_raw_or_mask_and_on_misaligned_row_keys():
    for mutate, expect in ((lambda g: g["warm"]["fit"].update(cost=101.0), "4.2(c) cost"),
                           (lambda g: g["held_engine_raw"]["warm"].update({GUARD_LAB: [0.1, float("nan"), 0.32]}),
                            "4.2(c) held_out_raw"),
                           (lambda g: g["held_engine_raw"]["warm"].update({GUARD_LAB: [0.1, 0.2, 0.3]}),
                            "4.2(c) held_out_raw"),
                           (lambda g: g["row_keys"]["train"].reverse(), "training row keys differ")):
        world = _guard_world()
        mutate(world[3])
        rg = _run_guard(world)
        assert not rg["pass"] and any(expect in x for x in h0._failed_guard_checks(rg)), expect


def test_the_freeze_check_covers_the_engine_files():
    for f in ("scripts/generators/fit_thiol_lolo.py", "scripts/generators/fit_joint_deoxyosone_thiol.py",
              "scripts/generators/fit_deoxyosone_submodel.py", "scripts/generators/build_submodel_probe.py",
              "src/kinetic_core/fitting.py", "scripts/generators/fit_thiol_lolo_pairs.py", h0.GENERATED_BY):
        assert f in h0.FROZEN_FILES


CAL = "CALIBRATED, NOT SHOWN BETTER THAN A NULL (D3 power = 0.05, engine and null errors simulated independent)"


def test_fragile_reads_the_two_sets_combined_by_precedence_per_dropped_lab():
    labs = ("l1", "l2", "l3")
    # warm alone would be 'fragile' (dropping l2 moves it to NOT CALIBRATED), but reset is INCONCLUSIVE whatever is
    # dropped, so the combined verdict never moves: not FRAGILE
    jack = {"warm": {"l1": CAL, "l2": "NOT CALIBRATED", "l3": CAL},
            "reset_where_leak_prone": {lab: "INCONCLUSIVE" for lab in labs}}
    res = h0.combined_fragility(jack, h0.combine_by_precedence([CAL, "INCONCLUSIVE"]))
    assert not res["fragile"] and set(res["combined_per_drop"].values()) == {"INCONCLUSIVE"}
    # here the combined verdict moves when l3 is dropped in both sets
    jack = {"warm": {"l1": CAL, "l2": CAL, "l3": "ADEQUATE on EX-056's 5 labs, nested"},
            "reset_where_leak_prone": {"l1": CAL, "l2": CAL, "l3": "ADEQUATE on EX-056's 5 labs, nested"}}
    res = h0.combined_fragility(jack, CAL)
    assert res["fragile"] and res["changes_when_dropped"] == ["l3"]
    # one set moving alone does not move the combination if the other set governs
    jack["reset_where_leak_prone"]["l3"] = CAL
    assert not h0.combined_fragility(jack, CAL)["fragile"]


def _rel(passes):
    return {null: {"lpd_pass": p, "sign_trigger": False} for null, p in zip(lolo.NULLS, passes)}


def test_the_relative_branch_note_fires_unless_better_than_both_nulls_in_every_set():
    both = {"warm": _rel((True, True)), "reset_where_leak_prone": _rel((True, True))}
    assert not h0.relative_note_fires("(i)", both)
    one_set = {"warm": _rel((True, True)), "reset_where_leak_prone": _rel((True, False))}
    assert h0.relative_note_fires("(i)", one_set)          # the worse set governs
    one_null = {"warm": _rel((True, False)), "reset_where_leak_prone": _rel((True, False))}
    assert h0.relative_note_fires("(i)", one_null)
    assert not h0.relative_note_fires("(ii)", one_null)     # only when D2 ships the relative branch


def test_d2_bound_reporting_at_the_grid_edges_and_when_undefined():
    top = h0.d2_bound_report(np.array([1.60]), np.array([0.5]))
    assert top["U_text"] == "U >= 1.60" and top["U_used"] == 1.60 and top["L_text"] == "L = 0.50"
    low = h0.d2_bound_report(np.array([np.nan]), np.array([np.nan]))
    assert low["U_text"] == "U <= 0.10" and low["U_used"] == 0.10 and low["U"] is None and low["L_text"] == "L: none"
    assert h0.d2_bound_report(0.45, 0.10)["L_text"] == "L <= 0.10"
    # invert: tau_obs below the 5 % quantile at every tau0 -> U undefined -> 'U <= 0.10'
    dist = np.tile(np.linspace(1.0, 2.0, 200), (len(h0.D2_GRID), 1))
    U, L = h0.invert(np.array([0.0]), dist)
    assert np.isnan(U[0]) and h0.d2_bound_report(U, L)["U_text"] == "U <= 0.10"


def _d2_set(tau, U, multi, unstable=False):
    return {"tau_ABS_all": tau, "f1": h0.d2_bound_report(U, np.nan), "multi_row_tau_ABS_all": multi,
            "unstable": unstable}


def test_d2_decision_is_nan_safe_and_order_free():
    a, b = _d2_set(0.6, np.nan, np.nan), _d2_set(0.7, 1.2, 0.5)
    for order in ((a, b), (b, a)):
        dec = h0.d2_decision({f"s{i}": x for i, x in enumerate(order)}, n5_pass=True, unstable_informative=False)
        assert dec["branch"] == "(ii)" and dec["U_used"] == 1.2 and dec["multi_row_tau_ABS_all_used"] == 0.5
    # every U undefined: each is used as 0.10, so U <= 0.8 -> branch (iii) (formerly the dead `is not None` check)
    dec = h0.d2_decision({"s0": _d2_set(0.3, np.nan, 0.2), "s1": _d2_set(0.2, np.nan, np.nan)}, True, False)
    assert dec["branch"] == "(iii)" and dec["U_used"] == 0.10
    # N5 failed: branch (i) needs the multi-row tau > 0.8; a NaN multi-row tau must not satisfy it in either order
    for order in ((_d2_set(0.9, 1.6, np.nan), _d2_set(0.95, 1.6, 0.7)), (_d2_set(0.95, 1.6, 0.7), _d2_set(0.9, 1.6, np.nan))):
        dec = h0.d2_decision({f"s{i}": x for i, x in enumerate(order)}, n5_pass=False, unstable_informative=False)
        assert dec["branch"] == "(ii)"
    with pytest.raises(SystemExit):
        h0.d2_decision({"s0": _d2_set(np.nan, 1.0, 0.5), "s1": _d2_set(0.5, 1.0, 0.5)}, True, False)


def test_score_runs_once_and_numbers_its_amendments(tmp_path, monkeypatch):
    assert h0.score_output("", tmp_path) == {"stem": h0.STEM, "kind": "primary", "amendment": ""}
    pre = h0.score_output("A1 fixes", tmp_path)
    assert pre["stem"] == h0.STEM and pre["kind"] == "pre-scoring amendment" and pre["amendment"] == "A1 fixes"
    (tmp_path / f"{h0.STEM}.json").write_text("{}")
    with pytest.raises(SystemExit) as exc:
        h0.score_output("", tmp_path)
    assert "exists" in str(exc.value.code)
    with pytest.raises(SystemExit):
        h0.score_output("   ", tmp_path)
    assert h0.score_output("why", tmp_path)["stem"] == f"{h0.STEM}_amendment_1"
    (tmp_path / f"{h0.STEM}_amendment_1.json").write_text("{}")
    (tmp_path / f"{h0.STEM}_amendment_2.md").write_text("")
    assert h0.score_output("why again", tmp_path)["stem"] == f"{h0.STEM}_amendment_3"
    # the CLI refuses before anything is loaded
    monkeypatch.setattr(h0, "OUT_DIR", tmp_path)
    monkeypatch.setattr(h0, "PAIR_DIR", tmp_path)
    with pytest.raises(SystemExit) as exc:
        h0.main(["--score"])
    assert "exists" in str(exc.value.code)


def _nest_world(converged_ab=True):
    """Outer lab_a; inner lab_b (one row only lab_a and lab_b measure: no strict null answers it) and lab_c (every
    strict engine residual NaN: emptied by the strict masks)."""
    terms = {"class_noise": {"sida": 0.1}, "between_tau": {"sida": 0.3}, "n_training_rows_scored": 1}
    rows = {lab: _abs_rows(lab, 2, [2, 2]) for lab in ("lab_a", "lab_b", "lab_c", "lab_d")}
    rows["lab_a"].append(_row("pa9s", "ONLY_AB", "sida", 0, lab="lab_a"))
    rows["lab_b"].append(_row("pb9s", "ONLY_AB", "sida", 0, lab="lab_b"))
    for lab in ("lab_b", "lab_c"):     # the nulls (trained on lab_d) can answer the other species
        for x in rows[lab]:
            if x.species != "ONLY_AB":
                rows["lab_d"].append(_row(f"pd{x.species}", x.species, "sida", 0, lab="lab_d"))
    data = tuple(x for lab in rows for x in rows[lab])
    rng = np.random.default_rng(8)
    folds = {lab: h0.Fold(lab=lab, held=tuple(r), paired=np.ones(len(r), bool),
                          raw={"engine_warm": rng.standard_normal((1, len(r)))}, terms={"engine_warm": terms},
                          leak_prone=[]) for lab, r in rows.items()}
    pairs = {}
    for inner in ("lab_b", "lab_c"):
        train, held = pair_runner.split(data, ("lab_a", inner))
        raw = {lab: rng.standard_normal((1, len(held[lab]))) for lab in held}
        if inner == "lab_c":
            raw["lab_c"][:] = np.nan
        pairs[("lab_a", inner)] = h0.Pair(labs=("lab_a", inner), held=held, train=train, raw_held={"warm": raw},
                                          raw_train={}, converged={"warm": converged_ab or inner != "lab_b"},
                                          leak_prone=[])
    return h0.Nesting(data, folds, pairs, None, pair_terms=lambda p, v, lab: terms, check_null_masks=False)


def test_a_mask_emptied_inner_lab_takes_the_minimum_data_rule_and_caps_nothing():
    nest = _nest_world()
    est = nest.calibrate("lab_a", "engine", "warm", inner_labs=("lab_a", "lab_b", "lab_c"))
    assert est["nonconverged_inner"] == []
    assert any("lab_c: no row left by the strict masks" in f and "minimum-data" in f for f in est["flags"])
    assert est["ABS"]["status"] == "single-lab" and est["labs"] == ["lab_b"]
    # the one-deep sensitivity uses the same strict mask: lab_b loses its ONLY_AB row there too, lab_c is emptied
    one = nest.calibrate("lab_a", "engine", "warm", nesting="onedeep", inner_labs=("lab_a", "lab_b", "lab_c"))
    strict_b = next(b for b in est["blocks"] if b.lab == "lab_b")
    one_b = next(b for b in one["blocks"] if b.lab == "lab_b")
    assert one_b.rows == strict_b.rows and all(x.species != "ONLY_AB" for x in one_b.rows)
    assert len(one_b.rows) == len(nest.folds["lab_b"].held) - 1 and one["labs"] == ["lab_b"]


def test_only_a_non_converged_pair_fit_caps_the_verdict():
    nest = _nest_world(converged_ab=False)
    est = nest.calibrate("lab_a", "engine", "warm", inner_labs=("lab_a", "lab_b", "lab_c"))
    assert est["nonconverged_inner"] == ["pair lab_a+lab_b warm: not converged (status <= 0)"]
    assert len(est["missing_inner"]) == 2
    d3 = {null: {"lpd_pass": True, "sign_trigger": False} for null in lolo.NULLS}
    assert h0._set_verdict(True, d3, bool(est["nonconverged_inner"]), True, 0.05) == "INCONCLUSIVE"
    assert h0._set_verdict(True, d3, False, True, 0.05).startswith("ADEQUATE")


def test_the_not_shown_better_label_names_the_independence_of_the_simulated_errors():
    d3 = {"typical_yield": {"lpd_pass": False, "sign_trigger": False}, "no_change": {"lpd_pass": True, "sign_trigger": False}}
    assert h0._set_verdict(True, d3, False, True, 0.054) == CAL
