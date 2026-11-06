"""Unit tests for scripts/generators/build_s4_band_sources.py (plan sec. 7, S4 band sourcing)."""
from __future__ import annotations

import json
import math

import pytest

from scripts.generators import build_s4_band_sources as G

CLASSES = {G.CARB: {"n": 3, "median": 100.0, "lo": 30.0, "hi": 237.0},
           G.SULF: {"n": 3, "median": 100.0, "lo": 29.675, "hi": 133.0}}


def _row(coordinate="log10_k_ref", value=-2.0, lane="sulfur", t_ref=145.0, bounds=(-10.0, 0.5), band=None):
    return {"lane": lane, "coordinate": coordinate, "value": value, "unit": "log10(1/min)",
            "reference_temperature_c": t_ref, "identification": {"bounds": list(bounds)}, "band": band,
            "sensitivity": {"species_of_its_steps": ["X"], "max_odorant_path_share": 0.5, "max_pot_flux_share": 0.1}}


def test_arrhenius_shift_sign_and_size():
    # 121 -> 145 C at 101.3 kJ/mol: k rises by ~5.9x (0.77 dex)
    assert G.arrhenius_shift_dex(101.3, 121.0, 145.0) == pytest.approx(0.770, abs=0.005)
    assert G.arrhenius_shift_dex(101.3, 145.0, 121.0) == pytest.approx(-0.770, abs=0.005)
    assert G.arrhenius_shift_dex(50.0, 145.0, 145.0) == 0.0


def test_upper_bound_uses_the_loose_edge_of_the_ea_band():
    src = {"relation": "upper", "value": 1e-3, "temperature_c": 120.0, "terms": {}}
    ev = G.evaluate_bound(src, t_ref_c=145.0, ea_band=(30.0, 237.0), between_lab=None, search=(-10, 0.5),
                          shipped=-2.0)
    loose = G.arrhenius_shift_dex(237.0, 120.0, 145.0)
    assert ev["upper"]["raw_edge"] == pytest.approx(-3.0 + loose)
    assert ev["upper"]["edge_95"] is None and ev["upper"]["binding"] is None  # no between-lab term: undefined
    assert "lower" not in ev


def test_lower_bound_uses_the_tight_end_and_the_lane_term_loosens_it():
    src = {"relation": "lower", "value": 1e-4, "temperature_c": 145.0, "terms": {"noise": 0.1}}
    ev = G.evaluate_bound(src, t_ref_c=145.0, ea_band=None, between_lab=0.964, search=(-10, 0.5), shipped=-3.0)
    sig = math.sqrt(0.1 ** 2 + 0.964 ** 2)
    assert ev["lower"]["edge_95"] == pytest.approx(-4.0 - 1.96 * sig)
    assert ev["lower"]["raw_edge"] == pytest.approx(-4.0 - 1.96 * 0.1)
    assert ev["lower"]["binding"] is True and ev["lower"]["shipped_outside_95"] is False


def test_bound_at_another_temperature_without_an_ea_band_cannot_be_carried():
    src = {"relation": "upper", "value": 1e-3, "temperature_c": 120.0, "terms": {}}
    ev = G.evaluate_bound(src, t_ref_c=145.0, ea_band=None, between_lab=0.964, search=(-10, 0.5), shipped=-2.0)
    assert "undefined" in ev["transport"]


def test_no_sources_is_refused_never_a_search_bound():
    e = G.propose("k_x.log10_k_ref", {"sources": [], "transport": G.SULF}, _row(), classes=CLASSES,
                  lane_between={"sulfur": 0.964}, dex_per_kj=0.015)
    assert e["proposal"]["route"] == "refused" and e["consequence"] == "refuse"


def test_lane_without_between_lab_term_refuses_a_raw_bound(monkeypatch):
    monkeypatch.setitem(G.SOURCES, "tmp_upper", dict(
        paper_keys=["x"], where="w", lab="a", conditions="c", relation="upper", value=1e-2, temperature_c=160.0,
        terms={}, reader="first", usable=True, why="y", conversion="z"))
    e = G.propose("k_x.log10_k_ref", {"sources": ["tmp_upper"], "transport": None},
                  _row(lane="acrylamide", t_ref=160.0, bounds=(-9, 0)), classes=CLASSES, lane_between={"sulfur": 0.964},
                  dex_per_kj=0.015)
    assert e["proposal"]["route"] == "refused"
    assert "between-lab term" in e["proposal"]["basis"]


def test_two_sided_measurement_band(monkeypatch):
    monkeypatch.setitem(G.SOURCES, "tmp_two", dict(
        paper_keys=["x"], where="w", lab="a", conditions="c", relation="two_sided", value=1e-2, temperature_c=145.0,
        terms={}, reader="second", usable=True, why="y", conversion="z"))
    e = G.propose("k_x.log10_k_ref", {"sources": ["tmp_two"], "transport": G.CARB}, _row(), classes=CLASSES,
                  lane_between={"sulfur": 0.964}, dex_per_kj=0.015)
    p = e["proposal"]
    assert p["route"] == "measurement_two_sided"
    assert p["lo"] == pytest.approx(-2.0 - 1.96 * 0.964) and p["hi"] == pytest.approx(-2.0 + 1.96 * 0.964)
    assert e["consequence"] == "widen"  # fixed today
    assert any("single laboratory" in f for f in e["flags"])
    assert any("second reader" in f for f in e["flags"])


def test_ea_class_band_and_consequences():
    legit = {"source": "b55", "distribution": "empirical_measured_class", "lo": 30.0, "hi": 237.0, "sampled": True}
    same = G.propose("k_x.ea_kj_mol", {"sources": [], "ea_class": G.CARB},
                     _row("ea_kj_mol", 64.0, band=legit), classes=CLASSES, lane_between={}, dex_per_kj=0.015)
    assert same["proposal"]["route"] == "class_distribution" and same["consequence"] == "unchanged"
    gigl = {"source": "b8", "distribution": "normal", "lo": 7.0, "hi": 102.0, "sampled": True}
    shift = G.propose("k_x.ea_kj_mol", {"sources": [], "ea_class": G.SULF},
                      _row("ea_kj_mol", 101.3, band=gigl), classes=CLASSES, lane_between={}, dex_per_kj=0.015)
    assert shift["consequence"] == "shift" and any("upgrade" in f for f in shift["flags"])
    outside = G.propose("k_x.ea_kj_mol", {"sources": [], "ea_class": G.CARB},
                        _row("ea_kj_mol", 20.0, band=gigl), classes=CLASSES, lane_between={}, dex_per_kj=0.015)
    assert any("SHIPPED barrier" in f for f in outside["flags"])
    none = G.propose("k_x.ea_kj_mol", {"sources": [], "ea_class": None},
                     _row("ea_kj_mol", 64.0), classes=CLASSES, lane_between={}, dex_per_kj=0.015)
    assert none["proposal"]["route"] == "refused"


def test_curated_tables_are_consistent():
    for cid, spec in {**G.COORDINATES, **G.EXTRA}.items():
        for s in spec.get("sources", []):
            assert s in G.SOURCES, (cid, s)
        if cid.endswith(".ea_kj_mol"):
            assert spec.get("ea_class") in (G.CARB, G.SULF), cid
    for sid, s in G.SOURCES.items():
        assert s["reader"] in ("first", "second"), sid
        assert s["paper_keys"] and s["lab"] and s["why"], sid
        assert "data/articles" not in json.dumps(s), sid  # the repo is public
        if s.get("usable"):
            assert s["relation"] in ("upper", "lower", "two_sided"), sid
            assert s["value"] > 0 and s["temperature_c"] is not None, sid
    assert len(G.COORDINATES) == 43


def test_committed_artifact_matches_its_record():
    """The committed JSON lists the audit's 43, and its markdown is the generator's rendering of it."""
    audit = json.loads(G.AUDIT.read_text(encoding="utf-8"))
    listed = sorted(e["coordinate"] for e in audit["summary"]["s4_sensitive_unidentified"])
    assert listed == sorted(G.COORDINATES)
    payload = json.loads(G.OUT.read_text(encoding="utf-8"))
    assert sorted(e["coordinate"] for e in payload["coordinates"]) == listed
    for e in payload["coordinates"]:
        if e["coordinate_kind"] == "ea_kj_mol" and e["proposal"]["route"] == "class_distribution":
            c = payload["classes"][e["proposal"]["class"]]
            assert (e["proposal"]["lo"], e["proposal"]["hi"]) == (c["lo"], c["hi"])
    assert G.OUT.with_suffix(".md").read_text(encoding="utf-8") == G.markdown(payload)
