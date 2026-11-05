"""H0 / EX-059, D2 branch (i), applied 2026-11-05: the thiol lane ships relative within-laboratory predictions only.

`predict` refuses the absolute level of every in-scope observable BY NAME and still answers a co-requested
out-of-scope target; `compare` answers the between-condition ratio with omega (and the measurement noise) in its
interval exactly once, resolves it only when that interval excludes 1, and withholds the arms' absolute levels; every
refusal names EX-059 and what would lift it; only a calibration made by `calibrate` from data rescues a level. The
scope itself is re-derived here from H0's own rows and the engine's own envelope, so it cannot drift from either.
"""
from __future__ import annotations

import json
import math

import pytest

from src import api, data_paths
from src.comparative_cli import render_compare_core_text, render_predict_core_text, spec_to_core, validate_spec
from src.kinetic_core import thiol_relative_only as h0
from src.kinetic_core.calibration import Calibration, ResponseFactor
from src.kinetic_core.engine import SULFUR, TARGET_ALIASES, declare_envelope, predict
from src.report_html import render_compare_report, render_predict_report

POT = {"name": "ribose_cysteine", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
       "temp_C": 145.0, "time_min": 20.0, "ph": 5.0, "aw": 0.98}
COOLER = dict(POT, name="ribose_cysteine_130C", temp_C=130.0)
IN_SCOPE = {"MFT": "2-methyl-3-furanthiol", "FFT": "2-furfurylthiol", "MFTD": "bis(2-methyl-3-furyl) disulfide",
            "FUR": "furfural", "NF": "norfuraneol", "MGO": "methylglyoxal"}
OUT_OF_SCOPE = "2-acetylthiazole"     # H0 did not score it (one unpaired Jiangnan row, D-all only)
SCORABLE = ("dalian_polytechnic_chen_d", "dfa_hofmann_schieberle", "jiangnan_zhang_x", "ochanomizu_murata",
            "reading_mottram")
MFT, FFT = IN_SCOPE["MFT"], IN_SCOPE["FFT"]


def _core(spec):
    return spec_to_core(validate_spec(spec, label="s"))


def _calibration(**overrides):
    fields = dict(lab="test_lab", base_wave="b60", created="2026-11-05", matrix="water",
                  response_factors={MFT: ResponseFactor(MFT, 0.5, 0.1, 4)}, overrides=(),
                  fit_records=("fit_1", "fit_2"), validate_records=(),
                  provenance={"generated_by": h0.CALIBRATE_GENERATOR})
    fields.update(overrides)
    return Calibration(**fields)


# --- the scope ----------------------------------------------------------------------------------------------------

def test_the_scope_is_what_h0_scored_and_the_lane_emits():
    """EX-056's paired rows of the five scorable laboratories carry H0's ship-rule quantities (scorer choice 19);
    the ruling covers exactly those species that the shipped sulfur lane ANSWERS for a pentose + cysteine pot, read
    from declare_envelope (review M1), not from the target-lane table."""
    scored = set()
    for lab in SCORABLE:
        record = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_thiol_lolo_{lab}.json").read_text())
        scored |= {row["species"] for row in record["rows"] if row["paired"]}
    assert scored == {"MFT", "FFT", "MFTD", "FUR", "NF", "MGO", "TDP"}
    names = {}
    for name, key in TARGET_ALIASES.items():
        names.setdefault(key, name)
    core = _core(POT)
    emitted = set()
    for key in scored:
        if key not in names:
            continue                                   # TDP: no target name at all
        declaration = declare_envelope(core, [names[key]])
        if declaration.is_answerable and SULFUR in declaration.lanes:
            emitted.add(key)
    assert h0.RELATIVE_ONLY_SPECIES == emitted == set(IN_SCOPE)
    actz = declare_envelope(core, [OUT_OF_SCOPE])
    assert actz.is_answerable and SULFUR in actz.lanes and not actz.absolute_refused_targets


def test_a_sulfur_lane_mgo_is_refused_and_a_trunk_mgo_is_answered():
    sulfur = api.predict(POT, targets=["methylglyoxal", OUT_OF_SCOPE])
    assert sulfur["declaration"]["lane"] == SULFUR and [r["compound"] for r in sulfur["rows"]] == [OUT_OF_SCOPE]
    assert [r["target"] for r in sulfur["refused_targets"]] == ["methylglyoxal"]
    trunk = api.predict(dict(POT, precursors={"D-Glucose": 10.0, "Glycine": 10.0}), targets=["methylglyoxal"])
    assert trunk["declaration"]["lane"] == "trunk" and trunk["answered"]
    assert [r["compound"] for r in trunk["rows"]] == ["methylglyoxal"] and not trunk["refused_targets"]


def test_the_ruling_is_read_from_the_frozen_record():
    r = h0.ruling()
    assert r.branch == "(i)" and r.tau_abs_dex == pytest.approx(0.963970, abs=1e-6)
    assert r.tau_abs_set == "reset_where_leak_prone" and r.upper_bound_text == "U >= 1.60"
    assert r.omega_dex == pytest.approx(0.418685, abs=1e-6) and r.omega_set == "warm"   # the larger set's omega
    assert r.class_noise_dex == pytest.approx(0.0973, abs=1e-4)                         # Reading fold, calibrated


def test_the_refusal_names_ex059_what_would_lift_it_and_the_transfer():
    text = h0.refusal_text()
    assert text.startswith(h0.ABSOLUTE_REFUSED_TAG) and "EX-059" in text
    assert "0.964" in text and "U >= 1.60" in text and "one laboratory" in text
    assert "second independent laboratory" in text and "WUR" in text and "H0's rule" in text
    assert "R3En" in text and "B60" in text and "conservative transfer" in text
    assert "EX-059" in h0.short_reason() and "EX-059" in h0.relative_caveat()
    assert "the null is better" in h0.relative_caveat() and "NOT SHOWN BETTER THAN A NULL" in h0.relative_caveat()


# --- predict ------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("key", sorted(IN_SCOPE))
def test_predict_refuses_each_in_scope_observable_by_name_and_answers_the_rest(key):
    name = IN_SCOPE[key]
    payload = api.predict(POT, targets=[name, OUT_OF_SCOPE])
    assert payload["answered"] is True
    assert [r["compound"] for r in payload["rows"]] == [OUT_OF_SCOPE]
    refused = {r["target"]: r for r in payload["refused_targets"]}
    assert set(refused) == {name} and refused[name]["absolute_refused"] is True
    assert "EX-059" in " ".join(refused[name]["reasons"])
    assert key not in (payload["oav_table"].get("per_species") or {})
    assert name in payload["declaration"]["absolute_refused_targets"]


def test_a_lone_in_scope_request_emits_no_number():
    payload = api.predict(POT, targets=[MFT])
    assert payload["answered"] is False and payload["rows"] == []
    assert [r["target"] for r in payload["refused_targets"]] == [MFT]


def test_the_renderers_print_no_in_scope_level():
    targets = [MFT, FFT, IN_SCOPE["FUR"], IN_SCOPE["MFTD"], OUT_OF_SCOPE]
    values = predict(_core(POT), targets).concentrations_ug_per_l
    levels = {c: v for c, v in values.items() if c != OUT_OF_SCOPE}
    p = api.predict(POT, targets=targets)
    c = api.compare(POT, COOLER, targets=targets)
    for text in (render_predict_core_text(p), render_predict_report(p), render_compare_core_text(c),
                 render_compare_report(c), json.dumps(p, default=str), json.dumps(c, default=str)):
        for compound, v in levels.items():
            for formatted in {repr(float(v)), f"{v:.6g}", f"{v:.4g}", f"{v:.2f}", f"{v:.3g}"}:
                if len(formatted.replace(".", "")) < 4:
                    continue                           # "15" is everywhere in a page (font sizes, dates)
                assert formatted not in text, (compound, formatted)


# --- compare ------------------------------------------------------------------------------------------------------

def test_compare_answers_the_ratio_with_omega_and_noise_once_and_withholds_the_levels():
    payload = api.compare(POT, COOLER, targets=[MFT, OUT_OF_SCOPE])
    comparison = payload["comparison"]
    assert comparison["comparable"] is True
    rows = {r["compound"]: r for r in comparison["ratios"]["rows"]}
    mft, actz = rows[MFT], rows[OUT_OF_SCOPE]
    ratio = mft["ratio_a_over_b"]
    assert ratio > 0 and math.isfinite(ratio)
    lo, hi = mft["interval_90_ratio_a_over_b"]
    r = h0.ruling()
    once = h0.Z90 * math.sqrt(2.0 * r.omega_dex ** 2 + 2.0 * r.class_noise_dex ** 2)   # prereg sec. 2 / sec. 5, score 6
    assert math.log10(hi / ratio) == pytest.approx(once, rel=1e-9)
    assert math.log10(ratio / lo) == pytest.approx(once, rel=1e-9)
    assert math.log10(hi / ratio) < h0.Z90 * 2.0 * math.hypot(r.omega_dex, r.class_noise_dex) * 0.99   # not twice
    assert mft["omega_dex_per_arm"] == r.omega_dex and "EX-059" in mft["relative_caveat"]
    assert "interval_90_ratio_a_over_b" not in actz                 # out of scope: no omega
    assert "thiol_relative_only" in payload["caveats"] and "no tighter interval" not in payload["caveats"]["ratio"]
    for arm in ("a", "b"):
        assert MFT not in comparison[f"run_{arm}"]["concentrations_ug_per_l"]
        assert not set(comparison[f"run_{arm}"]["species_mmol_per_l"]) & h0.RELATIVE_ONLY_SPECIES
        assert [x["compound"] for x in comparison[f"rows_{arm}"]] == [OUT_OF_SCOPE]
        assert "MFT" not in comparison[f"oav_table_{arm}"]["per_species"]
    assert set(comparison["absolute_refused"]) == {MFT}


def test_a_thiol_ratio_is_resolved_only_when_its_interval_excludes_one():
    comparison = api.compare(POT, COOLER, targets=[MFT, IN_SCOPE["MFTD"]])["comparison"]
    rows = {r["compound"]: r for r in comparison["ratios"]["rows"]}
    for row in rows.values():
        lo, hi = row["interval_90_ratio_a_over_b"]
        spans = lo <= 1.0 <= hi
        assert row.get("interval_spans_one", False) == spans
        if spans:
            assert row["within_reliability_band"] is True
    assert rows[MFT]["interval_spans_one"]                       # ~1.9x against a factor-10 interval
    resolved = sum(1 for r in rows.values() if r["direction"] != "undefined" and not r["within_reliability_band"])
    assert comparison["ratios"]["n_resolved"] == resolved
    assert "90 % interval spans 1" in render_compare_core_text(api.compare(POT, COOLER, targets=[MFT]))


# --- calibrations -------------------------------------------------------------------------------------------------

def test_a_calibration_made_by_calibrate_is_answered_with_a_two_row_ratio_interval():
    calibration = _calibration()
    payload = api.predict(POT, targets=[MFT, FFT], calibration=calibration)
    assert [r["compound"] for r in payload["rows"]] == [MFT]                      # the lab's own level
    assert [r["target"] for r in payload["refused_targets"]] == [FFT]             # no factor: refused
    extra = payload["run_metadata"]["calibration_extra_decades"][MFT]
    omega = h0.ruling().omega_dex
    assert extra == pytest.approx(h0.Z90 * math.sqrt(0.1 ** 2 + omega ** 2 + omega ** 2 / 4), rel=1e-12)
    compared = api.compare(POT, COOLER, targets=[MFT, FFT], calibration=calibration)["comparison"]
    assert MFT in compared["run_a"]["concentrations_ug_per_l"] and FFT not in compared["run_a"]["concentrations_ug_per_l"]
    assert set(compared["absolute_refused"]) == {FFT}


@pytest.mark.parametrize("change", [
    {"response_factors": {MFT: ResponseFactor(MFT, 0.0, 0.0, 0)}},       # the review's hand-written file
    {"response_factors": {MFT: ResponseFactor(MFT, 0.0, 0.1, 0)}},       # no rows
    {"response_factors": {MFT: ResponseFactor(MFT, 0.0, 0.0, 3)}},       # no sigma
    {"fit_records": ()},                                                 # no fit records
    {"provenance": {}},                                                  # not made by calibrate
    {"response_factors": {}},                                            # empty
])
def test_a_hand_written_or_empty_calibration_does_not_rescue_a_level(change):
    calibration = _calibration(**change)
    payload = api.predict(POT, targets=[MFT], calibration=calibration)
    assert payload["answered"] is False and payload["rows"] == [] and payload["refused_targets"][0]["absolute_refused"]
    from src.kinetic_core.user_scoring import score_document

    scored = score_document({"systems": [dict(POT, measured={MFT: {"value": 200.0, "unit": "ppb"}})]},
                            calibration=calibration)["systems"][0]
    assert scored["rows"] == [] and scored["refused"][0]["absolute_refused"]
    compared = api.compare(POT, COOLER, targets=[MFT], calibration=calibration)["comparison"]
    assert MFT not in compared["run_a"]["concentrations_ug_per_l"] and set(compared["absolute_refused"]) == {MFT}


def test_calibrate_shows_no_shipped_absolute_level_in_its_hold_out_rows():
    from src.kinetic_core import user_fit as F

    document = {"systems": [dict(POT, name="held", role="validate",
                                 measured={MFT: {"value": 200.0, "unit": "ppb"}, OUT_OF_SCOPE: {"value": 2.0, "unit": "ppb"}})]}
    records = F.records_from_document(document)
    assert F.engine_predict(records[0].spec, [MFT, OUT_OF_SCOPE], None).absolute_refused == (MFT,)
    before = F.Fitter([], records, F.engine_predict).holdout_scores(None)
    rows = {r["compound"]: r for r in before["rows"]}
    assert rows[MFT]["predicted_ppb"] is None and rows[MFT]["fold_error"] is None
    assert rows[MFT]["refused"] == h0.ABSOLUTE_REFUSED_TAG
    assert rows[OUT_OF_SCOPE]["predicted_ppb"] > 0 and before["n"] == 1


def test_the_score_verb_keeps_the_measurement_and_refuses_the_level():
    from src.kinetic_core.user_scoring import score_document

    document = {"systems": [dict(POT, measured={MFT: {"value": 200.0, "unit": "ppb"},
                                                OUT_OF_SCOPE: {"value": 5.0, "unit": "ppb"}})]}
    scored = score_document(document)["systems"][0]
    assert [r["compound"] for r in scored["rows"]] == [OUT_OF_SCOPE]
    refused = scored["refused"][0]
    assert refused["compound"] == MFT and refused["absolute_refused"] is True
    assert refused["measured_ppb"] == 200.0 and "EX-059" in refused["reason"]
