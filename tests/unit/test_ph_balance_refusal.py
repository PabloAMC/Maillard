"""The S3 application of M1b part (i), 2026-11-07: outputs that depend on the falsified pH balance are refused.

A refused output carries the reason (CF-084, CF-134, the M1b prereg sec. 21, S3) and the lift clause; an output the
balance does not touch is unchanged, bit for bit; a pot whose pH is held (the measured-trajectory path) is answered as
before; the rule covers every output the sensitivity record marks dependent; and the sensitivity generator's
evaluation and analysis run on a toy pot.
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


#: unbuffered cysteine + ribose, pH 4.5, 120 C, 60 min: the record's grid cell where FFT (0.74 dex) and furfural (0.96)
#: move beyond omega between B0 and the carried rivals, MFT (0.21) does not
UNBUFFERED = {"name": "cys_rib_water_pH4.5", "precursors": {"cysteine": 50.0, "ribose": 35.0}, "temp_C": 120.0,
              "time_min": 60.0, "ph": 4.5, "aw": 0.98, "buffer": {"kind": "none"}}
FFT, FUR = "2-furfurylthiol (FFT)", "furfural"


def test_a_refused_output_carries_the_reason_and_the_lift_text():
    run = engine.predict(_core(UNBUFFERED), [FFT, FUR])
    assert not run.answered and not run.concentrations_ug_per_l
    with pytest.raises(engine.OutOfEnvelope):
        run.require(FFT)
    reason = next(r for r in run.declaration.reasons if r.startswith(s3.BALANCE_REFUSED_TAG))
    for needle in ("CF-084", "CF-134", "prereg sec. 21", "S3", "WHAT WOULD LIFT IT", "union of the remaining carried",
                   "clamped", FFT, FUR):
        assert needle in reason, needle
    assert "0.419" in reason
    # the CLI says it; a compare across the pH lever is not a ratio
    payload = predict_core(UNBUFFERED, targets=[FFT])
    assert payload["answered"] is False
    assert any(r.startswith(s3.BALANCE_REFUSED_TAG) for r in payload["declaration"]["reasons"])
    cmp = compare_core(UNBUFFERED, dict(UNBUFFERED, name="pH5.5", ph=5.5), targets=[FFT])
    assert cmp["comparison"]["comparable"] is False


def test_some_targets_refused_by_name_the_rest_answered():
    run = engine.predict(_core(UNBUFFERED), [MFT, FUR])
    assert run.answered and set(run.concentrations_ug_per_l) == {MFT}
    assert FUR in run.run_metadata["refused_targets"]
    assert any(w.startswith("NOT ANSWERED, BY NAME") and s3.BALANCE_REFUSED_TAG in w for w in run.declaration.warnings)
    # the balance's end pH is withheld (its spread over the points is ~2 pH, against 0.478 unbuffered)
    seg = run.run_metadata["segments"][-1]
    assert seg["ph_cooled_end"] is None and seg["ph_in_situ_end"] is None and "ph_refused" in seg


def test_meynier_pots_keep_their_thiols_and_lose_their_ph():
    """The data that falsified the balance: 0.2 M pyrophosphate damps the thiols' spread below omega, not the pH's."""
    run = engine.predict(_core(MEYNIER), [MFT, "2-acetylthiazole"])
    assert run.answered and set(run.concentrations_ug_per_l) == {MFT, "2-acetylthiazole"}
    assert run.run_metadata["balance_s3"]["ph_refused"] is True
    assert run.run_metadata["segments"][-1]["ph_cooled_end"] is None


def test_the_short_reason_names_its_basis_and_its_lift():
    text = s3.short_reason()
    assert text.startswith(s3.BALANCE_REFUSED_TAG)
    for needle in ("CF-084", "sec. 21", "S3", "lift clause"):
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


def test_the_rivals_are_the_points_the_screen_built():
    import fit_m1b_ph_balance as FM

    best = data_access.load_json(s3.FALSIFICATION_RECORD)["best_draws"]
    breakpoints = data_access.load_json(s3.SCREEN_RECORD)["davidek_breakpoints_in_situ"]
    b1 = dict(best["B1"]["coordinates"], phi_cys=1.0)
    points = s3.rival_points()
    assert points["B1_phi_cys_1"] == FM.rival_changes("B1", "primary", b1, False, breakpoints)
    assert points["B13_best"] == FM.rival_changes("B13", "primary", best["B13"]["coordinates"], False, breakpoints)
    assert s3.EA_ISO_KJ_MOL == FM.EA_ISO_KJ_MOL


def test_the_sign_floor_sits_above_the_coupling_noise():
    import build_ph_balance_s3_sensitivity as gen

    noise = [e["b0_vs_inline"] for e in gen.read_detail()["outputs"]
             if e.get("kind") == "absolute" and e.get("b0_vs_inline") is not None]
    assert noise and max(noise) < s3.SIGN_FLOOR_DEX


def test_the_engine_reproduces_the_records_verdicts():
    """The record audits the rule the engine applies: on recipes it ran, the live verdict agrees output by output."""
    import build_ph_balance_s3_sensitivity as gen

    record = data_access.load_json(s3.SENSITIVITY_RECORD)
    assert record["summary"]["clamped_controls"]["all_independent"] is True
    by_rid = {}
    for e in gen.read_detail()["outputs"]:
        if e.get("kind") == "absolute" and e.get("dependent") is not None:
            by_rid.setdefault(e["rid"], {})[e["output"]] = e["dependent"]
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
           "B1": {MFT: 80.0, "furfural": 10.0, "ph": 4.9}, "B13": {MFT: 90.0, "furfural": 12.0, "ph": 5.0}}
    points = {"B0": ({}, None), "B0_inline": ({"x": 1}, None), "B1_phi_cys_1": ({"x": 2}, None),
              "B13_best": ({"x": 3}, None)}
    name_of = {0: "B0", 1: "B0_inline", 2: "B1", 3: "B13"}

    def fake(spec_, targets, changes, drift):
        v = toy[name_of[changes.get("x", 0)]]
        return {"answered": True, "concentrations_ug_per_l": {t: v[t] for t in targets},
                "ph_cooled_end": v["ph"], "ph_in_situ_end": v["ph"] + 0.2, "ph_mode": "dynamic"}

    row = gen.evaluate_recipe(recipe, points=points, predict_fn=fake)
    analysis = gen.analyse([row], [], gen.thresholds())
    by_output = {e["output"]: e for e in analysis["outputs"]}
    assert by_output[MFT]["dependent"] is False                    # log10(100 / 80) = 0.097 < omega
    assert by_output["furfural"]["dependent"] is True              # log10(50 / 10) = 0.70 > omega
    assert by_output["ph_cooled_end"]["dependent"] is True         # 0.9 pH > 0.135 (buffered)
    assert by_output["ph_cooled_end"]["tolerance"] == pytest.approx(gen.thresholds()["ph_buffered"])
    assert gen.log_spread([1.0, 0.0]) == float("inf") and gen.log_spread([0.0, 0.0]) == 0.0
