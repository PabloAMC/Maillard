"""H0 / EX-059, D2 branch (i), applied 2026-11-05: the thiol lane ships relative within-laboratory predictions only.

`predict` refuses the absolute level of every in-scope observable BY NAME and still answers a co-requested
out-of-scope target; `compare` answers the between-condition ratio with omega in its interval exactly once and
withholds the arms' absolute levels; every refusal names EX-059 and what would lift it. The scope itself is
re-derived here from H0's own rows, so it cannot drift from the record.
"""
from __future__ import annotations

import json
import math

import pytest

from src import api, data_paths
from src.kinetic_core import thiol_relative_only as h0
from src.kinetic_core.calibration import Calibration, ResponseFactor
from src.kinetic_core.engine import SULFUR, _TARGET_LANE

POT = {"name": "ribose_cysteine", "precursors": {"L-Cysteine": 10.0, "D-Ribose": 10.0},
       "temp_C": 145.0, "time_min": 20.0, "ph": 5.0, "aw": 0.98}
COOLER = dict(POT, name="ribose_cysteine_130C", temp_C=130.0)
IN_SCOPE = {"MFT": "2-methyl-3-furanthiol", "FFT": "2-furfurylthiol", "MFTD": "bis(2-methyl-3-furyl) disulfide",
            "FUR": "furfural", "NF": "norfuraneol"}
OUT_OF_SCOPE = "2-acetylthiazole"     # H0 did not score it (one unpaired Jiangnan row, D-all only)
SCORABLE = ("dalian_polytechnic_chen_d", "dfa_hofmann_schieberle", "jiangnan_zhang_x", "ochanomizu_murata",
            "reading_mottram")


def test_the_scope_is_what_h0_scored_and_the_lane_emits():
    """EX-056's paired rows of the five scorable laboratories carry H0's ship-rule quantities (scorer choice 19);
    the ruling covers exactly those species that the shipped sulfur lane answers as a target."""
    scored = set()
    for lab in SCORABLE:
        record = json.loads((data_paths.VALIDATION_DIR / f"kinetic_core_thiol_lolo_{lab}.json").read_text())
        scored |= {row["species"] for row in record["rows"] if row["paired"]}
    assert scored == {"MFT", "FFT", "MFTD", "FUR", "NF", "MGO", "TDP"}
    sulfur_targets = {k for k, lane in _TARGET_LANE.items() if lane == SULFUR}
    assert h0.RELATIVE_ONLY_SPECIES == scored & sulfur_targets == set(IN_SCOPE)
    assert "ACTZ" in sulfur_targets and "ACTZ" not in h0.RELATIVE_ONLY_SPECIES


def test_the_ruling_is_read_from_the_frozen_record():
    r = h0.ruling()
    assert r.branch == "(i)" and r.tau_abs_dex == pytest.approx(0.963970, abs=1e-6)
    assert r.tau_abs_set == "reset_where_leak_prone" and r.upper_bound_text == "U >= 1.60"
    assert r.omega_dex == pytest.approx(0.418685, abs=1e-6) and r.omega_set == "warm"   # the larger set's omega


def test_the_refusal_names_ex059_and_what_would_lift_it():
    text = h0.refusal_text()
    assert text.startswith(h0.ABSOLUTE_REFUSED_TAG) and "EX-059" in text
    assert "0.964" in text and "U >= 1.60" in text and "one laboratory" in text
    assert "second independent laboratory" in text and "WUR" in text and "H0's rule" in text
    assert "EX-059" in h0.short_reason() and "EX-059" in h0.relative_caveat()
    assert "not shown better than a null".lower() in h0.relative_caveat().lower()


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
    payload = api.predict(POT, targets=["2-methyl-3-furanthiol"])
    assert payload["answered"] is False and payload["rows"] == []
    assert [r["target"] for r in payload["refused_targets"]] == ["2-methyl-3-furanthiol"]


def test_compare_answers_the_ratio_with_omega_once_and_withholds_the_levels():
    payload = api.compare(POT, COOLER, targets=["2-methyl-3-furanthiol", OUT_OF_SCOPE])
    comparison = payload["comparison"]
    assert comparison["comparable"] is True
    rows = {r["compound"]: r for r in comparison["ratios"]["rows"]}
    mft, actz = rows["2-methyl-3-furanthiol"], rows[OUT_OF_SCOPE]
    ratio = mft["ratio_a_over_b"]
    assert ratio > 0 and math.isfinite(ratio)
    lo, hi = mft["interval_90_ratio_a_over_b"]
    omega = h0.ruling().omega_dex
    once = h0.Z90 * math.sqrt(2.0) * omega             # omega in each of the two arms, in quadrature, once
    assert math.log10(hi / ratio) == pytest.approx(once, rel=1e-9)
    assert math.log10(ratio / lo) == pytest.approx(once, rel=1e-9)
    assert math.log10(hi / ratio) < h0.Z90 * 2.0 * omega * 0.99   # not twice (linear sum per arm)
    assert mft["omega_dex_per_arm"] == omega and "EX-059" in mft["relative_caveat"]
    assert "interval_90_ratio_a_over_b" not in actz                 # out of scope: no omega
    assert "thiol_relative_only" in payload["caveats"]
    for arm in ("a", "b"):
        assert "2-methyl-3-furanthiol" not in comparison[f"run_{arm}"]["concentrations_ug_per_l"]
        assert not set(comparison[f"run_{arm}"]["species_mmol_per_l"]) & h0.RELATIVE_ONLY_SPECIES
        assert [r["compound"] for r in comparison[f"rows_{arm}"]] == [OUT_OF_SCOPE]
        assert "MFT" not in comparison[f"oav_table_{arm}"]["per_species"]
    assert set(comparison["absolute_refused"]) == {"2-methyl-3-furanthiol"}


def test_a_laboratory_calibration_with_a_response_factor_is_answered_with_omega_once():
    sigma = 0.1
    calibration = Calibration(
        lab="test_lab", base_wave="b60", created="2026-11-05", matrix="water",
        response_factors={"2-methyl-3-furanthiol": ResponseFactor("2-methyl-3-furanthiol", 0.5, sigma, 4)},
        overrides=(), fit_records=(), validate_records=(),
    )
    payload = api.predict(POT, targets=["2-methyl-3-furanthiol", "2-furfurylthiol"], calibration=calibration)
    assert [r["compound"] for r in payload["rows"]] == ["2-methyl-3-furanthiol"]     # the lab's own level
    assert [r["target"] for r in payload["refused_targets"]] == ["2-furfurylthiol"]  # no factor: refused
    extra = payload["run_metadata"]["calibration_extra_decades"]["2-methyl-3-furanthiol"]
    assert extra == pytest.approx(h0.Z90 * math.hypot(sigma, h0.ruling().omega_dex), rel=1e-12)


def test_the_score_verb_keeps_the_measurement_and_refuses_the_level():
    from src.kinetic_core.user_scoring import score_document

    document = {"systems": [dict(POT, measured={"2-methyl-3-furanthiol": {"value": 200.0, "unit": "ppb"},
                                                OUT_OF_SCOPE: {"value": 5.0, "unit": "ppb"}})]}
    scored = score_document(document)["systems"][0]
    assert [r["compound"] for r in scored["rows"]] == [OUT_OF_SCOPE]
    refused = scored["refused"][0]
    assert refused["compound"] == "2-methyl-3-furanthiol" and refused["absolute_refused"] is True
    assert refused["measured_ppb"] == 200.0 and "EX-059" in refused["reason"]
