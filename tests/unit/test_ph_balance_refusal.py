"""The S3 application of M1b part (i), 2026-11-07: outputs that depend on the falsified pH balance are refused.

Prereg sec. 21.10 (R1 - R3): a refused LEVEL carries the reason (CF-084, CF-134, the M1b prereg sec. 21, S3) and the
lift clause and stays in the arm for ratios, which the ratio rule (omega, sign) decides; compare's state vector carries
no refused species; the balance's pH is withheld against the pH-level term in any pot; an output the balance does not
touch is unchanged, bit for bit; a held pH (the measured-trajectory path) is answered as before; the rivals' points are
the screen's passing extremes; the engine reproduces the record; the generator runs on a toy pot.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

from src import data_access
from src.comparative_cli import compare_core, predict_core, spec_to_core, validate_spec
from src.kinetic_core import engine
from src.kinetic_core import ph_balance_refusal as s3

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "generators"))

MEYNIER = {"name": "cys_rib_pyro_pH5", "precursors": {"cysteine": 50.0, "ribose": 35.0}, "temp_C": 140.0,
           "time_min": 60.0, "ph": 5.0, "aw": 0.98, "buffer": {"kind": "pyrophosphate", "pyrophosphate_mol_l": 0.2}}
CLAMPED = dict(MEYNIER, name="cys_rib_clamped_pH5", buffer={"kind": "clamped"})
TRUNK = {"name": "glc_gly", "precursors": {"glucose": 100.0, "glycine": 100.0}, "temp_C": 140.0, "time_min": 30.0, "aw": 0.98,
         "ph": 6.0}
MFT = "2-methyl-3-furanthiol (MFT)"


def _core(spec):
    return spec_to_core(validate_spec(spec, label="s"))


#: unbuffered cysteine + ribose, pH 4.5, 120 C, 60 min: a grid cell of the record whose levels move with the carried
#: rivals beyond tau (results/validation/ph_balance_s3_sensitivity.md)
UNBUFFERED = {"name": "cys_rib_water_pH4.5", "precursors": {"cysteine": 50.0, "ribose": 35.0}, "temp_C": 120.0,
              "time_min": 60.0, "ph": 4.5, "aw": 0.98, "buffer": {"kind": "none"}}
FFT, FUR, ACTZ, H2S = "2-furfurylthiol (FFT)", "furfural", "2-acetylthiazole", "hydrogen sulfide"
LEVELS = [MFT, FFT, FUR, ACTZ, H2S, "bis(2-methyl-3-furyl) disulfide"]


def test_a_refused_level_carries_the_reason_and_the_lift_text_and_stays_for_ratios():
    run = engine.predict(_core(UNBUFFERED), LEVELS)
    refused = run.declaration.balance_refused_targets
    assert run.answered and refused, run.run_metadata.get("balance_s3")
    for name in refused:
        assert name in run.declaration.absolute_refused_targets
        assert name in run.concentrations_ug_per_l          # kept: a ratio between pots still needs it
    payload = predict_core(UNBUFFERED, targets=LEVELS)
    by_name = {r["target"]: r for r in payload["refused_targets"]}
    for name in refused:
        reason = " ".join(by_name[name]["reasons"])
        assert reason.startswith(s3.BALANCE_REFUSED_TAG)
        for needle in ("CF-084", "CF-134", "prereg sec. 21", "S3", "WHAT WOULD LIFT IT",
                       "union of the remaining carried", "clamped", "0.964"):
            assert needle in reason, needle
        assert name not in {r["compound"] for r in payload["rows"]}


@pytest.fixture(scope="module")
def buffered_pair():
    """One compare of two 0.2 M pyrophosphate pots (145 C / 20 min against 120 C / 60 min), shared by the tests that
    read its rows (2026-11-13: each test ran the same two arms' 23 rival points, ~5 min apiece on CI). The targets are
    the union of what those tests ask: each compound's verdict is its own spread, so asking more changes none."""
    a = {"name": "A", "precursors": {"cysteine": 50.0, "ribose": 35.0}, "temp_C": 145.0, "time_min": 20.0, "ph": 6.5,
         "aw": 0.98, "buffer": {"kind": "pyrophosphate", "pyrophosphate_mol_l": 0.2}}
    b = dict(a, name="B", temp_C=120.0, time_min=60.0)
    return compare_core(a, b, targets=[H2S, ACTZ, FFT, "methanethiol", MFT])["comparison"]


@pytest.mark.slow
def test_the_ratio_rule_decides_every_ratio_and_the_state_vector_leaks_nothing(buffered_pair):
    """A refused level is not a refused ratio: compare keeps both arms' numbers and judges each ratio on omega and the
    sign; the arms' state vectors carry no species whose level S3 refused (review, PR #108)."""
    cmp = buffered_pair
    assert cmp["comparable"] is True
    for arm in ("a", "b"):
        decl = cmp[f"declaration_{arm}"]
        run = cmp[f"run_{arm}"]
        for name in decl["balance_refused_targets"]:
            key = decl["mapped_targets"][name]
            assert key not in run["species_mmol_per_l"], (arm, name)
            assert name not in run["concentrations_ug_per_l"], (arm, name)
    for row in cmp["ratios"]["rows"]:
        if row.get("balance_refused"):
            assert row["ratio_a_over_b"] is None and row["note"].startswith(s3.BALANCE_REFUSED_TAG)
        elif row.get("direction") != "undefined":
            assert row["ratio_a_over_b"] is not None


def test_the_balance_ph_is_withheld_against_the_level_term_in_a_buffered_pot():
    """R2: Meynier's 0.2 M pyrophosphate pots move their cooled pH ~0.9 between B0 and the rivals, beyond 0.478."""
    run = engine.predict(_core(MEYNIER), [MFT, ACTZ])
    s = run.run_metadata["balance_s3"]
    assert s["ph_level_term"] == pytest.approx(0.4776, abs=1e-3)
    assert s["ph_refused"] is True and run.run_metadata["segments"][-1]["ph_cooled_end"] is None


def test_the_short_reason_names_its_basis_and_its_lift():
    text = s3.short_reason()
    assert text.startswith(s3.BALANCE_REFUSED_TAG)
    for needle in ("CF-084", "sec. 21", "S3", "lift clause", "tau", "omega"):
        assert needle in text


@pytest.mark.parametrize("spec,targets", [
    (TRUNK, ["5-HMF", "melanoidins"]),
    ({"name": "asn_glc", "precursors": {"asparagine": 20.0, "glucose": 20.0}, "temp_C": 160.0, "time_min": 10.0, "aw": 0.98,
      "ph": 6.0}, ["acrylamide"]),
])
def test_an_output_the_balance_does_not_touch_is_unchanged(spec, targets):
    core = _core(spec)
    shipped = engine.predict(core, targets)
    unrefused = engine.predict(core, targets, balance_refusal=False)
    assert shipped.answered
    assert shipped.concentrations_ug_per_l == unrefused.concentrations_ug_per_l
    assert not any(r.startswith(s3.BALANCE_REFUSED_TAG) for r in shipped.declaration.reasons)


def test_a_held_ph_pot_still_runs_the_measured_trajectory_path():
    core = _core(CLAMPED)
    shipped = engine.predict(core, [MFT, "2-acetylthiazole"])
    unrefused = engine.predict(core, [MFT, "2-acetylthiazole"], balance_refusal=False)
    assert shipped.answered
    assert shipped.concentrations_ug_per_l == unrefused.concentrations_ug_per_l
    assert all(seg["ph_mode"] != "dynamic" for seg in shipped.run_metadata["segments"])
    assert not s3.balance_solved(core, shipped.declaration.lanes)


def test_the_rivals_are_the_screens_passing_extremes():
    """R3: every point is a stored screen draw that passes, built as the screen built it, and the record is the one
    the generator derives from the stored draws now."""
    import fit_m1b_ph_balance as FM
    import build_ph_balance_s3_sensitivity as gen

    record = data_access.load_json(s3.POINTS_RECORD)
    fresh = gen.build_points()
    assert fresh["rivals"] == record["rivals"] and fresh["counts"] == record["counts"]
    points = s3.rival_points()
    for rival, pts in record["rivals"].items():
        for p in pts:
            assert p["passes"]
            assert points[p["name"]] == FM.rival_changes(rival, "primary", p["coordinates"], False,
                                                         record["davidek_breakpoints_in_situ"])
        whys = {w for p in pts for w in p["why"]}
        for coord in record["counts"][rival]["coordinates"]:
            assert f"min {coord}" in whys and f"max {coord}" in whys
    assert s3.EA_ISO_KJ_MOL == FM.EA_ISO_KJ_MOL


def test_the_rival_budget_is_deterministic_and_refuses_on_record():
    """Fix 4: no wall clock. A pot whose rival cannot be integrated is refused for a deterministic solver reason:
    LSODA exceeding the evaluation budget, after which BDF retries under a fresh budget and either exceeds it too or
    fails on step size (Chen 2002's pH 5 ribose + cysteine arms in the record)."""
    import build_ph_balance_s3_sensitivity as gen
    from src.kinetic_core import ph_balance_rivals as pbr

    assert isinstance(pbr.RIVAL_RHS_BUDGET, int) and pbr.RIVAL_RHS_BUDGET > 0
    failed = [v["failed"] for r in gen.read_detail()["recipes"] for p, v in r["points"].items()
              if v.get("failed") and p not in ("B0", "B0_inline") and "skipped" not in v["failed"]]
    deterministic = ("RIVAL_RHS_BUDGET", "Required step size is less than spacing between numbers")
    assert failed and all(any(d in f for d in deterministic) for f in failed), failed[:3]


def test_the_sign_floor_sits_above_the_coupling_noise():
    import build_ph_balance_s3_sensitivity as gen

    noise = [e["b0_vs_inline"] for e in gen.read_detail()["outputs"]
             if e.get("kind") == "absolute" and e.get("b0_vs_inline") is not None]
    assert noise and max(noise) < s3.SIGN_FLOOR_DEX


@pytest.mark.slow
def test_the_engine_reproduces_the_records_verdicts():
    """The record audits the rule the engine applies: on recipes it ran, the live verdict agrees output by output."""
    import build_ph_balance_s3_sensitivity as gen

    record = data_access.load_json(s3.SENSITIVITY_RECORD)
    assert record["summary"]["clamped_controls"]["all_independent"] is True
    by_rid = {}
    for e in gen.read_detail()["outputs"]:
        if e.get("kind") == "absolute" and e.get("dependent") is not None:
            by_rid.setdefault(e["rid"], {})[e["output"]] = e["dependent"]   # against tau (R1)
    wanted = ["grid:cys50_rib35|pyrophosphate_0.2M|pH6.5|145C_20min", "grid:thi10_rib35|none|pH5.5|120C_60min"]
    recipes = {r.rid: r for r in gen.grid_recipes() if r.rid in wanted}
    for rid in wanted:
        recipe = recipes[rid]
        targets = gen._answerable_targets(recipe.spec, recipe.targets)
        verdict = s3.assess(recipe.spec, targets)
        for output, dependent in by_rid[rid].items():
            assert (output in verdict.refused) is dependent, (rid, output)


def test_the_sensitivity_generator_on_a_toy_pot():
    import build_ph_balance_s3_sensitivity as gen

    spec = _core(dict(MEYNIER, name="toy"))
    recipe = gen.Recipe("toy:cys_rib", "toy", spec, (MFT, "furfural"), {"buffer": "pyrophosphate_0.2M"})
    toy = {"B0": {MFT: 100.0, "furfural": 50.0, "ph": 4.1}, "B0_inline": {MFT: 100.0, "furfural": 50.0, "ph": 4.1},
           "R1": {MFT: 80.0, "furfural": 10.0, "ph": 4.9}, "R2": {MFT: 90.0, "furfural": 5.0, "ph": 4.4}}
    points = {"B0": ({}, None), "B0_inline": ({"x": 1}, None), "R1": ({"x": 2}, None), "R2": ({"x": 3}, None)}
    name_of = {0: "B0", 1: "B0_inline", 2: "R1", 3: "R2"}

    def fake(spec_, targets, changes, drift):
        v = toy[name_of[changes.get("x", 0)]]
        return {"answered": True, "concentrations_ug_per_l": {t: v[t] for t in targets},
                "ph_cooled_end": v["ph"], "ph_in_situ_end": v["ph"] + 0.2, "ph_mode": "dynamic"}

    row = gen.evaluate_recipe(recipe, points=points, predict_fn=fake)
    th = gen.thresholds()
    analysis = gen.analyse([row], [], th)
    by_output = {e["output"]: e for e in analysis["outputs"]}
    assert by_output[MFT]["dependent"] is False                    # log10(100 / 80) = 0.097 < tau
    assert by_output["furfural"]["dependent"] is True              # log10(50 / 5) = 1.0 > tau = 0.964
    assert by_output["ph_cooled_end"]["dependent"] is True         # 0.8 pH > 0.478 (R2: buffered or not)
    assert by_output["ph_cooled_end"]["tolerance"] == pytest.approx(th["ph_level"])
    assert by_output["ph_cooled_end"]["dependent_sensitivity"] is True    # > 0.265
    assert by_output[MFT]["dependent_sensitivity"] is False        # 0.097 < omega
    assert gen.log_spread([1.0, 0.0]) == float("inf") and gen.log_spread([0.0, 0.0]) == 0.0


@pytest.mark.slow
def test_a_thiol_refused_on_both_grounds_keeps_h0s_interval_when_its_ratio_stands(buffered_pair):
    """Confirmation read of PR #108: an H0 thiol whose level the balance ALSO refuses must keep H0's ratio treatment
    (the omega interval, the relative caveat, the not-resolved check) whenever its ratio stands under S3."""
    cmp = buffered_pair
    checked = 0
    for row in cmp["ratios"]["rows"]:
        if row["compound"] not in (MFT, FFT):      # the H0 thiols; the shared pair also holds non-H0 targets
            continue
        if row.get("direction") == "undefined" or row.get("balance_refused"):
            continue
        assert "interval_90_ratio_a_over_b" in row and row.get("relative_caveat"), row["compound"]
        checked += 1
    assert checked >= 1


def test_the_rivals_run_in_parallel_with_the_serial_verdict(monkeypatch):
    """The 23 rival integrations go through ``parallel.imap_ordered`` (each under its own Config): the verdict with
    workers is the serial one, bit for bit (2026-11-11: the serial loop doubled CI time and timed main out)."""
    spec = _core(UNBUFFERED)
    verdicts = []
    for n in ("1", "3"):
        monkeypatch.setenv("MAILLARD_WORKERS", n)
        s3._CACHE.clear()
        verdicts.append(s3.assess(spec, LEVELS))
    s3._CACHE.clear()
    serial, forked = verdicts
    assert serial.failure is None and forked.failure is None
    assert forked.values == serial.values and forked.ph == serial.ph
    assert forked.refused == serial.refused and forked.spread_dex == serial.spread_dex
