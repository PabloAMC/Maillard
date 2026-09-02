"""Step 6 (2026-09-03): bring your own measurement -- scored on the core, never calibrated.

What is pinned: the input contract (refusals, not repairs), that a measured row is scored
the way the panel scorecard scores a bundle, that an unrepresented compound is a refusal
with the engine's reason, that the record lands under results/ and never data/, and that
the record carries what a future fit wave needs.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from src import data_paths
from src.kinetic_core import user_scoring as us
from src.kinetic_core.engine import predict


@pytest.fixture(scope="module")
def document():
    return yaml.safe_load(us.TEMPLATE)


@pytest.fixture(scope="module")
def scored(document):
    return us.score_document(document)


def test_the_template_is_a_valid_document(document):
    systems = us.systems_of(document)
    assert len(systems) == 1
    assert systems[0]["measured"]["2-methyl-3-furanthiol"]["value"] == 198.0


def test_a_measurement_without_a_value_or_in_another_unit_is_refused():
    with pytest.raises(us.MeasurementSpecError):
        us.systems_of({"systems": [{"name": "x", "precursors": {"D-Ribose": 10}, "temp_C": 140, "time_min": 30,
                                    "ph": 5, "aw": 0.98, "measured": {"furfural": {"unit": "ppb"}}}]})
    with pytest.raises(us.MeasurementSpecError):
        us.systems_of({"systems": [{"name": "x", "precursors": {"D-Ribose": 10}, "temp_C": 140, "time_min": 30,
                                    "ph": 5, "aw": 0.98, "measured": {"furfural": {"value": 1.0, "unit": "mg/kg"}}}]})
    with pytest.raises(us.MeasurementSpecError):
        us.systems_of({"systems": []})


def test_a_missing_process_condition_is_refused_not_defaulted():
    with pytest.raises(Exception):
        us.systems_of({"systems": [{"name": "x", "precursors": {"D-Ribose": 10}, "temp_C": 140,
                                    "measured": {"furfural": 1.0}}]})


def test_rows_are_scored_like_the_panel_scorecard(scored):
    system = scored["systems"][0]
    assert system["scored"] == 2 and not system["refused"]
    from src.comparative_cli import spec_to_core
    spec = us.systems_of(yaml.safe_load(us.TEMPLATE))[0]
    for row in system["rows"]:
        run = predict(spec_to_core(spec), [row["compound"]])
        assert row["predicted_ppb"] == pytest.approx(run.concentrations_ug_per_l[row["compound"]])
        assert row["fold_error"] >= 1.0
        assert row["within_band"] == (row["fold_error"] <= scored["pass_band_level"])
        assert row["interval_ug_per_L"][0] <= row["interval_ug_per_L"][1]
        assert row["lane"] == "sulfur"


def test_an_unrepresented_compound_is_a_refusal_with_the_engines_reason():
    document = yaml.safe_load(us.TEMPLATE)
    document["systems"][0]["measured"]["2-pentylfuran"] = {"value": 5.0, "unit": "ppb"}
    payload = us.score_document(document)
    system = payload["systems"][0]
    assert [x["compound"] for x in system["refused"]] == ["2-pentylfuran"]
    assert "UNREPRESENTED" in system["refused"][0]["reason"]
    assert payload["summary"]["refused_rows"] == 1 and payload["summary"]["scored_rows"] == 2


def test_records_land_under_results_and_carry_what_a_wave_needs(scored, tmp_path):
    written = us.write_records(scored, out_dir=tmp_path / "user")
    names = sorted(p.name for p in written)
    assert names == ["scored.json", "user_my_ribose_cysteine_run.json"]
    record = json.loads((tmp_path / "user" / "user_my_ribose_cysteine_run.json").read_text())
    assert record["evidence_class"] == "user_measurement"
    assert record["source_metadata"]["origin"] == "user_measurement"
    assert record["measured_volatiles"]["2-furfurylthiol"]["conc_ppb"] == 121.0
    assert record["content_verification"]["quantification_class"] == "stable_isotope_dilution_gcms"
    assert "fit_target_ids" in record["next_wave_note"]
    assert not str(us.USER_RESULTS_DIR).startswith(str(data_paths.DATA_ROOT))


def test_the_default_output_directory_is_under_results_not_data():
    assert us.USER_RESULTS_DIR == data_paths.RESULTS_ROOT / "user"


def test_text_render_names_every_row_and_says_it_is_not_calibration(scored):
    text = us.render_text(scored)
    for row in scored["systems"][0]["rows"]:
        assert row["compound"] in text
    assert "Scoring, not calibration" in text
