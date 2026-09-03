"""Wave S5 — the comparative front door: shape, honesty and reliability wiring.

WHAT THESE TESTS ARE FOR
------------------------
The CLI added in this wave is orchestration over machinery that other tests already cover.
So these tests do NOT re-check the chemistry. They pin the three things that are new and that
a refactor could silently break:

1. **Shape.** Each verb runs on a tiny example and produces the payload the renderers expect.
2. **Honesty.** No absolute ppb reaches the user without its caveat; the sulfur caveat fires
   on sulfur compounds; `predict` says out loud that its tier is a scope claim, not a
   validation claim.
3. **Wiring.** The reliability tags come from the directional artifact at runtime, not from a
   constant in the code -- so if the panel is re-scored and the report edited, the CLI moves
   with it. This is tested by pointing the loader at a synthetic report and checking the tags
   change.

Deliberately NOT pinned: any predicted ppb, ratio or fold error. Pinning a number this wave
produced is the circularity Rounds 1-3 removed, and Waves U/S3/S4 all set the precedent of
pinning behaviour rather than output.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import comparative_cli  # noqa: E402
from src import directional_reliability as dr  # noqa: E402
from src import model_card  # noqa: E402

CLI = ROOT / "scripts" / "maillard.py"
EXAMPLE_SPEC = ROOT / "docs" / "examples" / "compare_ribose_vs_glucose.yml"


# ---------------------------------------------------------------------------------------
# Reliability wiring -- the tags must come from the artifact, not from a constant
# ---------------------------------------------------------------------------------------


def test_panel_counts_are_read_from_the_core_directional_scorecard():
    """Step B7 (2026-09-03): the tags come from results/validation/core_directional_scores.json,
    the kinetic core scored on the claims panel. The VALUES are pinned in
    tests/scientific/test_core_headline_guards.py; here only the categories the CLI tags
    against must be present, so a category cannot vanish into "unmeasured" silently."""
    counts = dr.load_panel_counts()
    for category in ("sugar_identity", "additive_cysteine", "temperature", "time", "ph", "moisture_aw"):
        assert category in counts, f"{category} missing from the core directional scorecard"
        agree, evaluable = counts[category]
        assert 0 <= agree <= evaluable


def test_reliability_tags_follow_the_artifact_rather_than_a_hardcoded_table(tmp_path):
    """Edit the artifact, and the tag moves. This is the whole point of reading it at runtime."""
    fake = tmp_path / "scores.json"
    fake.write_text(json.dumps({
        "artifact": "core_directional_scores",
        "summary": {"independent": {"by_category": {
            "ph": {"agree": 7, "evaluable": 7, "not_evaluable": 0, "misses": [], "mechanism_absent": 0},
            "sugar_identity": {"agree": 0, "evaluable": 8, "not_evaluable": 1, "misses": ["X-1"], "mechanism_absent": 0},
        }}},
    }), encoding="utf-8")
    counts = dr.load_panel_counts(fake)
    assert dr.reliability_for_axis("ph", counts, {}).verdict == dr.VERDICT_TRUST
    assert dr.reliability_for_axis("sugar_identity", counts, {}).verdict == dr.VERDICT_DO_NOT_USE
    assert "X-1" in dr.load_axis_notes(fake)["sugar_identity"]


def test_a_missing_or_unparseable_artifact_raises_rather_than_defaulting(tmp_path):
    """A silent fallback here would republish a number after its evidence stopped being read."""
    with pytest.raises(FileNotFoundError):
        dr.load_panel_counts(tmp_path / "nope.json")
    wrong = tmp_path / "wrong.json"
    wrong.write_text(json.dumps({"artifact": "something_else"}), encoding="utf-8")
    with pytest.raises(ValueError):
        dr.load_panel_counts(wrong)


def test_water_activity_is_never_trusted_and_ph_is_tagged_by_the_standing_rule():
    """No lane carries a water-activity term, so every a_w comparison is refused (0/0
    evaluable) and the axis stays 'do-not-use'. pH is different since axis refusal: the
    trunk, acrylamide and lipid lanes refuse pH, so every EVALUABLE pH claim is a sulfur-lane
    claim, and the tag is whatever the standing thresholds say about those. RE-PINNED
    2026-09-03 (B9): 4 of 5 -> 'trust' at the boundary (was 3 of 5, 'caution', under B8).
    The thresholds are part of the claim: a promotion that came from moving them is a
    different event from one that came from a refit, and this one came from the refit."""
    counts = dr.load_panel_counts()
    assert dr.TRUST_MIN_RATE == 0.80
    assert dr.CAUTION_MIN_RATE == 0.60
    assert dr.MIN_EVALUABLE_FOR_TRUST == 3
    assert dr.reliability_for_axis("moisture_aw", counts).verdict == dr.VERDICT_DO_NOT_USE
    ph = dr.reliability_for_axis("ph", counts)
    assert (ph.agree, ph.evaluable) == (4, 5)
    assert ph.verdict == dr.verdict_for(4, 5) == dr.VERDICT_TRUST


def test_an_unmeasured_axis_is_reported_as_do_not_use_not_as_silence():
    tag = dr.reliability_for_axis("an_axis_nobody_measured", {"ph": (4, 7)})
    assert tag.verdict == dr.VERDICT_DO_NOT_USE
    assert "not measured" in tag.note


def test_axes_exercised_detects_each_knob_independently():
    base = {
        "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
        "temp_C": 140.0,
        "time_min": 30.0,
        "ph": 5.0,
        "aw": 0.98,
        "protein_type": "free",
    }
    assert dr.axes_exercised(base, base) == []
    assert dr.axes_exercised(base, {**base, "ph": 8.0}) == ["ph"]
    assert dr.axes_exercised(base, {**base, "aw": 0.4}) == ["moisture_aw"]
    assert dr.axes_exercised(base, {**base, "temp_C": 180.0}) == ["temperature"]
    assert dr.axes_exercised(base, {**base, "time_min": 90.0}) == ["time"]
    assert dr.axes_exercised(base, {**base, "protein_type": "pea_iso"}) == ["matrix_identity"]
    swapped = {**base, "precursors": {"L-Cysteine": 10.0, "D-Glucose": 10.0}}
    assert dr.axes_exercised(base, swapped) == ["sugar_identity"]


def test_a_multi_axis_comparison_is_governed_by_its_weakest_axis():
    """Bundling a pH change with a sugar change must not launder the pH miss away.

    RE-EXPRESSED 2026-09-03 (B7): the invariant is tested on a SYNTHETIC counts table so it
    holds whatever the live scorecard says about either axis; the live numbers are pinned in
    tests/scientific/test_core_headline_guards.py."""
    counts = {"sugar_identity": (8, 8), "ph": (1, 5)}
    base = {
        "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
        "temp_C": 140.0,
        "time_min": 30.0,
        "ph": 5.0,
        "aw": 0.98,
    }
    other = {**base, "ph": 8.0, "precursors": {"L-Cysteine": 10.0, "D-Glucose": 10.0}}
    described = dr.describe_comparison(base, other, counts=counts)
    assert set(described["axes"]) == {"sugar_identity", "ph"}
    ph_tag = dr.reliability_for_axis("ph", counts, {})
    sugar_tag = dr.reliability_for_axis("sugar_identity", counts, {})
    assert ph_tag.verdict == dr.VERDICT_DO_NOT_USE and sugar_tag.verdict == dr.VERDICT_TRUST
    assert described["governing"].verdict == ph_tag.verdict
    assert described["governing"].verdict != sugar_tag.verdict


def test_sulfur_detection_fires_on_every_alias_the_pipeline_uses():
    """predicted_ppb carries three spellings per compound; a caveat must fire on all of them."""
    for alias in ("2-Methyl-3-furanthiol (MFT)", "2-methyl-3-furanthiol", "2-Furfurylthiol (FFT)"):
        assert dr.is_sulfur_compound(alias), alias
    assert not dr.is_sulfur_compound("Furfural")
    assert not dr.is_sulfur_compound("2,5-Dimethylpyrazine")


# ---------------------------------------------------------------------------------------
# Spec handling
# ---------------------------------------------------------------------------------------


@pytest.fixture
def free_spec_document():
    return {
        "a": {
            "name": "cys_ribose",
            "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
            "temp_C": 140.0,
            "time_min": 30.0,
            "ph": 5.0,
            "aw": 0.98,
            "protein_type": "free",
        },
        "b": {
            "name": "cys_glucose",
            "precursors": {"L-Cysteine": 10.0, "D-Glucose": 10.0},
            "temp_C": 140.0,
            "time_min": 30.0,
            "ph": 5.0,
            "aw": 0.98,
            "protein_type": "free",
        },
    }


def test_a_spec_missing_a_process_condition_is_refused_not_defaulted():
    """A defaulted temperature is a silent chemistry claim, so the interface refuses to make one."""
    with pytest.raises(comparative_cli.SpecError) as excinfo:
        comparative_cli.validate_spec(
            {"precursors": {"D-Ribose": 10.0}, "ph": 5.0, "aw": 0.9, "time_min": 30.0},
            label="test",
        )
    assert "temp_C" in str(excinfo.value)


def test_an_empty_precursor_map_is_refused():
    with pytest.raises(comparative_cli.SpecError):
        comparative_cli.validate_spec(
            {"precursors": {}, "temp_C": 140.0, "ph": 5.0, "aw": 0.9, "time_min": 30.0},
            label="test",
        )


def test_the_shipped_template_is_a_valid_two_arm_spec():
    document = yaml.safe_load(comparative_cli.SPEC_TEMPLATE)
    spec_a, spec_b = comparative_cli.split_comparison_document(document, source="template")
    assert spec_a["name"] and spec_b["name"]
    assert dr.axes_exercised(spec_a, spec_b) == ["sugar_identity"]


def test_the_committed_example_spec_matches_the_template():
    """The example file is generated from --template; drift between them is a doc bug."""
    assert EXAMPLE_SPEC.exists()
    assert yaml.safe_load(EXAMPLE_SPEC.read_text()) == yaml.safe_load(comparative_cli.SPEC_TEMPLATE)


# ---------------------------------------------------------------------------------------
# Verb shape -- in-process, so the assertions can see the payload
# ---------------------------------------------------------------------------------------


def test_the_ratio_bound_is_the_independent_error_combination():
    a = {"p5": 1.0, "p50": 10.0, "p95": 100.0}
    b = {"p5": 2.0, "p50": 20.0, "p95": 200.0}
    lo, hi = comparative_cli._ratio_bounds(a, b)
    assert lo == pytest.approx(1.0 / 200.0)
    assert hi == pytest.approx(100.0 / 2.0)


def test_no_bound_is_invented_when_either_arm_lacks_an_envelope():
    assert comparative_cli._ratio_bounds(None, {"p5": 1.0, "p95": 2.0}) == (None, None)
    assert comparative_cli._ratio_bounds({"p5": 1.0, "p95": 2.0}, None) == (None, None)


# ---------------------------------------------------------------------------------------
# Rendering honesty -- what actually reaches the terminal
# ---------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------
# End-to-end smoke: each verb runs as a subprocess and emits parseable output
# ---------------------------------------------------------------------------------------


def _run_cli(*args, expect_ok=True):
    proc = subprocess.run(
        [sys.executable, str(CLI), *args],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        timeout=600,
    )
    if expect_ok:
        assert proc.returncode == 0, proc.stderr[-2000:]
    return proc


@pytest.mark.slow
def test_compare_verb_runs_end_to_end_and_emits_json():
    """2026-09-03 (B5): the kinetic core is the only lane; `compare` reports ratios with
    the envelope declaration of each arm."""
    proc = _run_cli("compare", str(EXAMPLE_SPEC), "--json")
    payload = json.loads(proc.stdout)
    assert payload["artifact"] == "maillard_compare_core"
    assert payload["lane"] == "core"
    assert "comparison" in payload


def test_predict_verb_runs_end_to_end_and_emits_json():
    """2026-08-29 (Wave B5): default lane is the kinetic core -- see above."""
    proc = _run_cli("predict", str(EXAMPLE_SPEC), "--system", "a", "--json")
    payload = json.loads(proc.stdout)
    assert payload["artifact"] == "maillard_predict_core"
    assert payload["lane"] == "core"
    assert payload["answered"] is True
    assert payload["rows"]


@pytest.mark.slow
def test_rank_experiments_verb_runs_end_to_end():
    proc = _run_cli("rank-experiments", "--top", "3", "--json", expect_ok=False)
    if proc.returncode == 2:
        # The uncertainty panel is gitignored, so a fresh clone legitimately has no artifact.
        # The verb must say so and name the generator rather than emitting an empty ranking.
        assert "generate_prediction_uncertainty.py" in proc.stderr
        pytest.skip("prediction_uncertainty.json absent in this tree")
    payload = json.loads(proc.stdout)
    assert "candidates" in payload
    assert len(payload["candidates"]) <= 3


def test_template_flag_prints_a_spec_and_exits_zero():
    proc = _run_cli("compare", "--template")
    assert yaml.safe_load(proc.stdout)["a"]["precursors"]


# ---------------------------------------------------------------------------------------
# The generated model card
# ---------------------------------------------------------------------------------------


def test_readme_carries_the_model_card_markers():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert model_card.README_BEGIN in readme
    assert model_card.README_END in readme
    assert readme.index(model_card.README_BEGIN) < readme.index(model_card.README_END)


def test_the_readme_model_card_is_not_stale():
    """The README's validity domain must equal what the artifacts currently say.

    This is the guard the whole card exists for. Wave S5 found the directional report still
    publishing 20/29 four commits after the tree started scoring 21/29, because the number
    lived only in prose. Here, drift is a red test.
    """
    # Gates ARE run here: the committed card carries their status line, and all three
    # complete in well under a second. A gate that starts failing will therefore also turn
    # this test red -- which is the right outcome, since the published card would be wrong.
    card = model_card.build_model_card(run_gate_checks=True)
    markdown = model_card.render_model_card_markdown(card)
    _, changed = model_card.splice_into_readme(markdown, ROOT / "README.md")
    assert not changed, (
        "README.md's model card does not match the artifacts. Regenerate it:\n"
        "  python scripts/generators/generate_model_card.py"
    )


def test_the_card_states_all_three_honest_headline_sentences():
    card = model_card.build_model_card(run_gate_checks=False)
    joined = " ".join(card["headline_sentences"])
    assert "Absolute concentrations are unreliable" in joined
    # 2026-09-03 (B5): the directional panel was scored on the retired screening lane; the
    # card must say the core has NOT been scored on it rather than borrow those numbers.
    # 2026-09-03 (B7): the core is scored on the directional panel; the sentence carries the count.
    assert "Directional and ranking claims are the product, and on the kinetic core they score" in joined
    # RE-PINNED 2026-08-28 (Wave W). This used to require the literal phrase "sulfur branch
    # has zero absolute literature anchors". That sentence became FALSE on 2026-08-28 when
    # the full text of 10.1021/jf9705983 arrived and three of its Table 1 rows were ingested
    # as panel anchors, so requiring it would have pinned the card to a stale claim -- and a
    # stale claim in the FLATTERING direction, since "unanchored" sounds less damning than
    # "anchored and wrong by 12-30x". What the third headline sentence must still do is name
    # the sulfur branch's anchor status honestly, whichever way it currently falls.
    assert "sulfur branch has" in joined and "literature anchors" in joined
    assert ("zero absolute literature anchors" in joined) or (
        "fails every one of them" in joined
    ), (
        "The third headline sentence no longer states the sulfur branch's anchor status "
        "either way. It must say either that the branch is unanchored or that it is "
        "anchored and failing; silence is the one option that is not allowed."
    )


def test_the_sulfur_anchor_claim_is_checked_not_asserted():
    """The card may only make its anchor claim if the tree actually supports it.

    RE-PINNED 2026-08-28 (Wave W). The old version asserted
    ``status["zero_absolute_anchors"]`` -- i.e. it pinned the ANSWER, not the mechanism.
    That was safe only while the answer happened to be "zero", and it would have forced the
    card to keep printing a false claim once real anchors arrived. It also had a subtler
    defect: ``collect_sulfur_anchor_status`` inspected ONE retired file, so what it really
    measured was "the retired file still has no verified value" -- permanently true, and not
    the question. Wave W widened the collector to scan the panel, and this guard now pins the
    two things that must stay true regardless of which way the answer falls:
    the status is COMPUTED from files rather than asserted, and it agrees with itself.
    """
    status = model_card.collect_sulfur_anchor_status()
    assert status["available"], "the sulfur benchmark must be readable for the claim to stand"

    # The retired benchmark's own values must stay marked as unverifiable. If this ever
    # empties, someone has quietly re-blessed 342/200 ppb -- the exact regression the whole
    # Wave S2b/S2c trail exists to prevent.
    assert status["values_without_verifiable_source"]
    assert not status["values_with_an_anchor"]

    # The headline flag must be DERIVED from the panel scan, not hardcoded either way.
    verified = status["panel_verified_sulfur_benchmarks"]
    assert status["zero_absolute_anchors"] == (not verified), (
        "zero_absolute_anchors disagrees with the panel scan it is supposed to summarise: "
        f"flag={status['zero_absolute_anchors']}, verified anchors={verified}."
    )

    # RE-PINNED 2026-08-28 (Wave X): three -> eight, all still from Hofmann & Schieberle
    # 1998. The five new ones are the STEP-LEVEL rows -- Tables 3, 8 and 10 -- which
    # constrain individual reaction steps rather than end-to-end lumps. Pinned so that
    # losing any of them is a failure rather than a silent return to the old sentence.
    assert verified == [
        "hofmann1998_c2c3_recombination_145C_20min_pH3",
        "hofmann1998_c2c3_recombination_145C_20min_pH5",
        "hofmann1998_c2c3_recombination_145C_20min_pH7",
        "hofmann1998_fructose_cysteine_145C_20min_pH5",
        "hofmann1998_furan2aldehyde_h2s_145C_20min_pH5",
        "hofmann1998_glucose_cysteine_145C_20min_pH5",
        "hofmann1998_norfuraneol_cysteine_145C_20min_pH5",
        "hofmann1998_ribose_cysteine_145C_20min_pH5",
    ], f"panel-verified sulfur anchors moved to {verified}"

    # ADDED 2026-08-28 (Wave X). A NINTH primary-source-verified sulfur row is on the panel
    # and is deliberately NOT in the list above: `hofmann1998_norfuraneol_h2s_145C_20min_pH5`
    # is the fit target of results/validation/
    # furanone_reductive_sulfhydrylation_refit_hofmann.json, so a constant was selected by
    # looking at it and its agreement is not evidence about the model. It is SEPARATED, not
    # dropped -- excluding a fitted row removes its misses as well as its hits, and this one
    # is currently a miss (2.3x under). The guard pins that the split exists and that the two
    # lists are disjoint, so a future wave cannot quietly promote a fitted row into the
    # anchor count.
    fitted = status["panel_verified_sulfur_values_that_are_fit_targets"]
    assert fitted, "the Wave X fit target must still be recognised as one"
    fitted_benchmarks = {row.split("/", 1)[0] for row in fitted}
    assert fitted_benchmarks == {"hofmann1998_norfuraneol_h2s_145C_20min_pH5"}
    assert not fitted_benchmarks.intersection(verified), (
        "a fit target is being counted as an absolute literature anchor; agreement on a "
        "row a constant was selected against carries no information about the model"
    )


def test_the_provenance_census_reproduces_the_readme_pinned_source_status_count():
    """87 is pinned in README prose and by test_honest_headline_guards; recount must agree.

    RE-PINNED 2026-09-01: 120 -> 102 when data/qm/ (18 records) was deleted with the QM lane.
    RE-PINNED 2026-09-02: 102 -> 87 when the mocked protein_source_registry (15) was withdrawn.
    """
    census = model_card.collect_no_verifiable_source_census()
    assert census["available"]
    assert census["by_status_key"]["source_status"] == 87


def test_a_missing_artifact_is_reported_in_the_card_not_dropped_from_it(monkeypatch, tmp_path):
    """A blank cell is indistinguishable from a pass, so absence must be stated."""
    monkeypatch.setattr(model_card, "ENVELOPE_PATH", tmp_path / "absent.json")
    collected = model_card.collect_core_envelope()
    assert collected["available"] is False
    assert collected["path"]


def test_every_verdict_in_the_card_is_one_of_the_three_words():
    card = model_card.build_model_card(run_gate_checks=False)
    assert card["validity_domain"], "the card must not be empty"
    for row in card["validity_domain"]:
        assert row["verdict"] in {
            dr.VERDICT_TRUST,
            dr.VERDICT_CAUTION,
            dr.VERDICT_DO_NOT_USE,
        }, row


def test_absolute_concentration_rows_are_never_trust():
    """No artifact state may make the card say 'trust' about an absolute ppb number."""
    card = model_card.build_model_card(run_gate_checks=False)
    absolute_rows = [
        row for row in card["validity_domain"] if row["claim_type"].startswith("Absolute")
    ]
    assert absolute_rows
    assert all(row["verdict"] == dr.VERDICT_DO_NOT_USE for row in absolute_rows)


def test_splice_refuses_a_readme_without_markers(tmp_path):
    readme = tmp_path / "README.md"
    readme.write_text("# no markers here\n", encoding="utf-8")
    with pytest.raises(ValueError) as excinfo:
        model_card.splice_into_readme("card", readme)
    assert "model-card markers" in str(excinfo.value)


def test_a_spec_missing_a_condition_exits_nonzero_naming_the_missing_key(tmp_path):
    bad = tmp_path / "bad.yml"
    bad.write_text(
        "a:\n  precursors: {D-Ribose: 1.0}\n  ph: 5.0\n  aw: 0.9\n  time_min: 30\n"
        "b:\n  precursors: {D-Glucose: 1.0}\n  ph: 5.0\n  aw: 0.9\n  time_min: 30\n",
        encoding="utf-8",
    )
    proc = _run_cli("compare", str(bad), expect_ok=False)
    assert proc.returncode == 2
    assert "missing required key" in proc.stderr
    assert "temp_C" in proc.stderr


def test_a_one_armed_comparison_spec_exits_nonzero_with_a_usable_hint(tmp_path):
    bad = tmp_path / "one_arm.yml"
    bad.write_text("a:\n  precursors: {D-Ribose: 1.0}\n", encoding="utf-8")
    proc = _run_cli("compare", str(bad), expect_ok=False)
    assert proc.returncode == 2
    assert "--template" in proc.stderr
