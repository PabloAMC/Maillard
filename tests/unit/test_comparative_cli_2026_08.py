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
EXAMPLE_SPEC = ROOT / "data" / "cli_examples" / "compare_ribose_vs_glucose.yml"


# ---------------------------------------------------------------------------------------
# Reliability wiring -- the tags must come from the artifact, not from a constant
# ---------------------------------------------------------------------------------------


def test_panel_counts_are_read_from_the_committed_directional_report():
    counts = dr.load_panel_counts()
    # The categories the CLI tags against must all be present. Their VALUES are not pinned
    # here -- the report is allowed to move; what must not happen is a category vanishing and
    # the CLI silently reporting "unmeasured".
    for category in (
        "sugar_identity",
        "additive_cysteine",
        "temperature",
        "time",
        "lipid_lane",
        "matrix_identity",
        "ph",
        "moisture_aw",
    ):
        assert category in counts, f"{category} missing from the CURRENT STANDING table"
        agree, evaluable = counts[category]
        assert 0 <= agree <= evaluable


def test_reliability_tags_follow_the_artifact_rather_than_a_hardcoded_table(tmp_path):
    """Edit the report, and the tag moves. This is the whole point of reading it at runtime."""
    fake = tmp_path / "report.md"
    fake.write_text(
        "# CURRENT STANDING — synthetic\n\n"
        "| category | agree | evaluable |\n|---|---:|---:|\n"
        "| ph | 7 | 7 |\n| sugar_identity | 0 | 8 |\n",
        encoding="utf-8",
    )
    counts = dr.load_panel_counts(fake)
    assert dr.reliability_for_axis("ph", counts).verdict == dr.VERDICT_TRUST
    assert dr.reliability_for_axis("sugar_identity", counts).verdict == dr.VERDICT_DO_NOT_USE


def test_a_missing_or_unparseable_report_raises_rather_than_defaulting(tmp_path):
    """A silent fallback here would republish a number after its evidence stopped being read."""
    with pytest.raises(FileNotFoundError):
        dr.load_panel_counts(tmp_path / "nope.md")
    empty = tmp_path / "empty.md"
    empty.write_text("# something else\n", encoding="utf-8")
    with pytest.raises(ValueError):
        dr.load_panel_counts(empty)


def test_ph_and_water_activity_are_never_trusted_on_the_shipped_panel():
    """The report's own conclusion, enforced. It says the model 'licenses no pH recommendation'.

    RE-PINNED 2026-08-28 (Wave W), and the reason matters more than the new value.

    This test used to assert ``ph`` == ``do-not-use``. It no longer holds, because the
    EVIDENCE moved, not because the threshold did: Wave W added two independent pH rows
    from Mottram & Nobrega 2002 (``MOT-01``, ``MOT-02``, both quoted from the paper's
    Table 1) and the model AGREED with both, taking ``ph`` from 4/7 = 0.571 to 6/9 = 0.667.
    ``CAUTION_MIN_RATE`` is still 0.60 and ``TRUST_MIN_RATE`` is still 0.80; neither was
    touched in that wave, and the guard below asserts that they weren't.

    So the assertion is re-expressed as the thing the report actually licenses -- **pH and
    water activity are never TRUSTED** -- rather than as one specific tag. That is the
    invariant §8/§A6 depends on ("licenses no pH recommendation, no moisture
    recommendation"), it survives the panel gaining rows in either direction, and it still
    fails loudly if someone widens a threshold to promote either axis.
    """
    counts = dr.load_panel_counts()

    # The thresholds themselves are part of the claim: a "promotion" that came from moving
    # these is a different event from one that came from new measurements.
    assert dr.TRUST_MIN_RATE == 0.80
    assert dr.CAUTION_MIN_RATE == 0.60

    for axis in ("ph", "moisture_aw"):
        assert dr.reliability_for_axis(axis, counts).verdict != dr.VERDICT_TRUST, (
            f"{axis} is now tagged 'trust'. The report's licence section says the model "
            f"licenses no pH and no moisture recommendation; promoting either axis to "
            f"'trust' contradicts a shipped licence and must be done in the report first."
        )

    # Water activity has moved for nobody: 0/3, and the miss is structural (see
    # src/directional_reliability._AXIS_NOTES).
    assert dr.reliability_for_axis("moisture_aw", counts).verdict == dr.VERDICT_DO_NOT_USE

    # pH, as of Wave W, cleared the caution floor by two rows and no more. Pinned so that a
    # drop back below 0.60 is visible as a change rather than absorbed silently.
    # RE-PINNED 2026-08-28 (Wave X): 6/9 -> 6/10, and the axis is now EXACTLY ON the floor
    # (0.600 against CAUTION_MIN_RATE 0.60). NO THRESHOLD MOVED -- both are asserted above.
    # The new row is `HOX-03`: Hofmann & Schieberle 1998 Table 8 measures MFT rising 20x from
    # pH 3 to pH 7 in a system containing no amino acid and no sugar (hydroxyacetaldehyde +
    # mercapto-2-propanone), and the model predicts the IDENTICAL value at pH 3, 5 and 7,
    # because that lane's three families appear in none of the pH sets in src/conditions.py.
    # ONE MORE pH MISS TAKES THIS AXIS BACK TO do-not-use. That is not a reason to avoid
    # adding pH rows; it is the reason to add them.
    assert dr.reliability_for_axis("ph", counts).verdict == dr.VERDICT_CAUTION
    assert counts["ph"] == (6, 10)


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
    """Bundling a pH change with a sugar change must not launder the pH miss away."""
    counts = dr.load_panel_counts()
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

    # RE-EXPRESSED 2026-08-28 (Wave W). This used to assert the governing verdict was the
    # literal string "do-not-use", which was true only because `ph` happened to sit there.
    # Wave W's two new independent pH rows (MOT-01/MOT-02, both agreeing) moved `ph` to
    # `caution` on the evidence, with no threshold change -- and the old assertion would
    # then have failed for a reason that has nothing to do with what this test is for.
    # What the test is actually for is the LAUNDERING invariant: the bundled comparison must
    # inherit the WEAKER of its two axes, never the stronger one. That is asserted directly
    # below and it holds at any tag either axis may take in future.
    ph_tag = dr.reliability_for_axis("ph", counts)
    sugar_tag = dr.reliability_for_axis("sugar_identity", counts)
    assert ph_tag.rate < sugar_tag.rate, (
        "This test needs pH to be the weaker of the two axes to say anything. It is not "
        f"any more (ph {ph_tag.rate:.3f} vs sugar_identity {sugar_tag.rate:.3f}); pick a "
        "different weak axis for the bundle rather than deleting the check."
    )
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


def test_spec_to_bench_reuses_the_benchmark_input_contract(free_spec_document):
    """Precursor values must land as concentration_mM, the unit the projection budget reads."""
    spec = comparative_cli.validate_spec(free_spec_document["a"], label="a")
    bench = comparative_cli._spec_to_bench(spec)
    assert bench["precursors"]["D-Ribose"]["concentration_mM"] == 10.0
    assert bench["conditions"]["water_activity"] == 0.98
    assert bench["conditions"]["time_min"] == 30.0


# ---------------------------------------------------------------------------------------
# Verb shape -- in-process, so the assertions can see the payload
# ---------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def free_runs():
    document = yaml.safe_load(comparative_cli.SPEC_TEMPLATE)
    spec_a, spec_b = comparative_cli.split_comparison_document(document, source="template")
    return (
        comparative_cli.evaluate_system(spec_a),
        comparative_cli.evaluate_system(spec_b),
    )


def test_each_compound_appears_exactly_once_despite_the_three_alias_keys(free_runs):
    run_a, _ = free_runs
    displays = list(run_a.compounds)
    assert len(displays) == len(set(displays))
    # The canonical SMILES spelling must not leak into the user-facing name.
    assert not any(display.startswith("O=C") or display.startswith("Cc1") for display in displays)


def test_compare_payload_has_the_expected_shape(free_runs):
    payload = comparative_cli.compare_systems(*free_runs)
    assert payload["artifact"] == "maillard_compare"
    assert payload["axes_exercised"] == ["sugar_identity"]
    assert payload["rows"], "the free-precursor sulfur system must produce compounds"
    for row in payload["rows"]:
        assert set(row) >= {
            "compound",
            "ratio",
            "ratio_lo",
            "ratio_hi",
            "dominant_pathway_a",
            "reliability",
            "reliability_verdict",
        }
        assert row["reliability_verdict"] in {
            dr.VERDICT_TRUST,
            dr.VERDICT_CAUTION,
            dr.VERDICT_DO_NOT_USE,
            "n/a",
        }


def test_at_least_one_compound_carries_a_dominant_pathway(free_runs):
    """The route trace side-channel is wired; without it every pathway column would be '-'."""
    payload = comparative_cli.compare_systems(*free_runs)
    assert any(row["dominant_pathway_a"] for row in payload["rows"])


def test_the_ratio_bound_is_the_independent_error_combination():
    a = {"p5": 1.0, "p50": 10.0, "p95": 100.0}
    b = {"p5": 2.0, "p50": 20.0, "p95": 200.0}
    lo, hi = comparative_cli._ratio_bounds(a, b)
    assert lo == pytest.approx(1.0 / 200.0)
    assert hi == pytest.approx(100.0 / 2.0)


def test_no_bound_is_invented_when_either_arm_lacks_an_envelope():
    assert comparative_cli._ratio_bounds(None, {"p5": 1.0, "p95": 2.0}) == (None, None)
    assert comparative_cli._ratio_bounds({"p5": 1.0, "p95": 2.0}, None) == (None, None)


def test_predict_reports_a_range_not_a_bare_point(free_runs):
    run_a, _ = free_runs
    payload = comparative_cli.predict_system(run_a)
    assert payload["artifact"] == "maillard_predict"
    assert payload["rows"]
    assert any(row["range_available"] for row in payload["rows"])
    for row in payload["rows"]:
        assert "range_p5" in row and "range_p95" in row


def test_the_sulfur_caveat_fires_on_a_sulfur_system(free_runs):
    payload = comparative_cli.compare_systems(*free_runs)
    assert any(row["sulfur"] for row in payload["rows"])
    assert payload["caveats"]["sulfur"]
    assert "zero absolute literature anchors" in payload["caveats"]["sulfur"]


# ---------------------------------------------------------------------------------------
# Rendering honesty -- what actually reaches the terminal
# ---------------------------------------------------------------------------------------


def test_ratios_lead_and_absolutes_do_not_appear_unless_asked(free_runs):
    payload = comparative_cli.compare_systems(*free_runs)
    default = comparative_cli.render_compare_text(payload, show_absolute=False)
    assert "A/B" in default
    assert "ABSOLUTE ppb (requested with --absolute)" not in default

    with_absolute = comparative_cli.render_compare_text(payload, show_absolute=True)
    assert "ABSOLUTE ppb (requested with --absolute)" in with_absolute
    # The caveat is not optional: it ships in the same block as the numbers.
    assert "ABSOLUTE ppb ARE NOT RELIABLE" in with_absolute


def test_predict_prints_the_absolute_caveat_unconditionally(free_runs):
    run_a, _ = free_runs
    text = comparative_cli.render_predict_text(comparative_cli.predict_system(run_a))
    flat = " ".join(text.split())  # the renderer wraps, so compare on normalised whitespace
    assert "ABSOLUTE ppb ARE NOT RELIABLE" in flat
    # And it must not let a run-level "high" tier read as a validation claim.
    # RE-PINNED 2026-08-28 (Wave W): 14 -> 17 as the panel gained three Hofmann anchors.
    # The NUMERATOR is the load-bearing half and it did not move: the three new rows are
    # literature anchors and all three fail their contracts, so nothing became strict-ready.
    assert "0 of 17 benchmarks in the panel are strict-ready" in flat


def test_predict_points_the_user_at_the_comparative_verb(free_runs):
    run_a, _ = free_runs
    text = comparative_cli.render_predict_text(comparative_cli.predict_system(run_a))
    assert "maillard compare" in text


def test_the_per_axis_note_reaches_the_rendered_table():
    """A do-not-use tag must arrive with its reason, not as a bare word."""
    counts = dr.load_panel_counts()
    payload = {
        "a": {"name": "A", "spec": {}},
        "b": {"name": "B", "spec": {}},
        "axes_exercised": ["ph"],
        "governing_reliability": "do-not-use (4/7)",
        "per_axis": [
            {
                "axis": "ph",
                "score": "%d/%d" % counts["ph"],
                "verdict": dr.VERDICT_DO_NOT_USE,
                "note": dr.reliability_for_axis("ph", counts).note,
            }
        ],
        "rows": [],
        "warnings_a": [],
        "warnings_b": [],
        "caveats": {"ratio": "r", "absolute": "a", "sulfur": None},
    }
    text = comparative_cli.render_compare_text(payload)
    assert "licenses no pH recommendation" in text


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
    proc = _run_cli("compare", str(EXAMPLE_SPEC), "--json")
    payload = json.loads(proc.stdout)
    assert payload["artifact"] == "maillard_compare"
    assert payload["axes_exercised"] == ["sugar_identity"]
    assert payload["rows"]


@pytest.mark.slow
def test_predict_verb_runs_end_to_end_and_emits_json():
    proc = _run_cli("predict", str(EXAMPLE_SPEC), "--system", "a", "--json")
    payload = json.loads(proc.stdout)
    assert payload["artifact"] == "maillard_predict"
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
    assert "Directional and ranking claims are the measured product" in joined
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
    """120 is pinned in README prose and by test_honest_headline_guards; recount must agree."""
    census = model_card.collect_no_verifiable_source_census()
    assert census["available"]
    assert census["by_status_key"]["source_status"] == 120


def test_a_missing_artifact_is_reported_in_the_card_not_dropped_from_it(monkeypatch, tmp_path):
    """A blank cell is indistinguishable from a pass, so absence must be stated."""
    monkeypatch.setattr(model_card, "HOLDOUT_PATH", tmp_path / "absent.json")
    collected = model_card.collect_free_precursor_holdout()
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
