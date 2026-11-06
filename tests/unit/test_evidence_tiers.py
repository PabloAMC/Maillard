"""Unit tests for the incumbent evidence audit reader (src/evidence_tiers.py; plan sec. 7, S1)."""
from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

from src import evidence_tiers as et


def _param(k_ref=1e-2, ea=100.0, t_ref_k=373.15):
    return SimpleNamespace(k_ref=k_ref, ea_kj_mol=ea, t_ref_k=t_ref_k)


def _dm(**over):
    dm = {
        "id": "DM-T", "constant": "k_x", "coordinate": "log10_k_ref", "lane": "trunk", "measurement": "test",
        "sources": ["src2026"], "ledger_rows": [], "ledger_note": "", "cite": "", "lab": "one lab", "system": "",
        "channels": "", "rate_law": "first order", "engine_rate_law": "first order", "concentration": "",
        "unstated": [], "relation": "same_channel",
        "measured": {"kind": "rate", "value": 1e-3, "unit": "1/min", "temperature_c": 100.0},
        "band": {"noise": 0.1, "within_lab_between_paper": None, "between_lab": 0.3, "reaction_class": None},
        "matches": {"channels": True, "rate_law": True, "concentration": True, "conditions_stated": True},
        "note": "",
    }
    dm.update(over)
    return dm


def test_arrhenius_and_ea_dex():
    assert et.arrhenius(1e-2, 100.0, 373.15, 100.0) == pytest.approx(1e-2)
    assert et.arrhenius(1e-2, None, 373.15, 145.0) == pytest.approx(1e-2)  # no barrier: held
    k145 = et.arrhenius(1e-2, 100.0, 373.15, 145.0)
    assert math.log10(k145 / 1e-2) == pytest.approx(et.ea_dex(100.0), rel=1e-9)


def test_interval_is_quadrature_and_needs_a_between_lab_term():
    sd, half, _ = et.comparison_interval_dex({"noise": 0.3, "between_lab": 0.4}, "trunk")
    assert sd == pytest.approx(0.5) and half == pytest.approx(0.98)
    sd, half, note = et.comparison_interval_dex({"noise": 0.3, "between_lab": None}, "trunk")
    assert sd is None and half is None and "undefined" in note
    sd, _, note = et.comparison_interval_dex({"noise": None, "between_lab": None}, "sulfur", fitted_sigma_dex=0.2)
    assert sd == pytest.approx(math.hypot(0.964, 0.2)) and "lane" in note


def test_local_criteria():
    pinned = et.local_unidentified_criteria(value=-2.0, lo=-10.0, hi=0.5, ci95_halfwidth=0.4, fit_says_identified=True)
    assert pinned == []
    edge = et.local_unidentified_criteria(value=20.04, lo=20.0, hi=260.0, ci95_halfwidth=61.6,
                                          fit_says_identified=None, is_barrier=True, ci_lower=-41.5)
    assert any("band edge" in c for c in edge) and any("unphysical sign" in c for c in edge)
    spans = et.local_unidentified_criteria(value=-4.0, lo=-6.0, hi=2.0, ci95_halfwidth=5.0, fit_says_identified=True)
    assert any("spans the band" in c for c in spans)
    assert any("flat" in c for c in et.local_unidentified_criteria(
        value=-4.0, lo=-10.0, hi=0.5, ci95_halfwidth=0.5, fit_says_identified=True, slice_verdict="flat"))


def test_subset_lower_bound_satisfied_is_not_a_contradiction():
    dm = _dm(relation="subset_lower_bound", measured={"kind": "lower_bound", "value": 1e-3, "unit": "1/min",
                                                       "temperature_c": 100.0})
    out = et.evaluate_measurement(dm, _param(k_ref=1e-2), "fitted_identified", "trunk")
    assert not out["violated"] and out["s2_status"].startswith("consistent")


def test_s2_targets_fitted_identified_only():
    dm = _dm()
    unid = et.evaluate_measurement(dm, _param(), "fitted_unidentified", "trunk")
    assert unid["violated"] and unid["s2_status"].startswith("not S2") and "S4" in unid["s2_status"]
    ident = et.evaluate_measurement(dm, _param(), "fitted_identified", "trunk", fitted_sigma_dex=0.1)
    assert ident["s2_status"].startswith("S2-eligible")  # 1 dex against a ~0.65 dex 95 % interval
    other = et.evaluate_measurement(dm, _param(), "analogue", "trunk")
    assert other["s2_status"].startswith("contradicted, not fitted")


def test_s2_blocks():
    blocked = et.evaluate_measurement(_dm(inferred_bound=True), _param(), "fitted_identified", "trunk")
    assert "(d)" in blocked["s2_status"]
    mismatch = et.evaluate_measurement(_dm(matches={"channels": False, "rate_law": True, "concentration": True,
                                                    "conditions_stated": True}), _param(), "fitted_identified", "trunk")
    assert "(a) channels" in mismatch["s2_status"]
    unstated = et.evaluate_measurement(_dm(depends_on_unstated=["O2"]), _param(), "fitted_identified", "trunk")
    assert unstated["interval95_halfwidth_dex"] is None and "(b)" in unstated["s2_status"]
    inside = et.evaluate_measurement(_dm(band={"noise": 1.0, "between_lab": 1.0}), _param(), "fitted_identified", "trunk")
    assert "includes zero" in inside["s2_status"]


def test_declarations_are_well_formed():
    for dm in et.DIRECT_MEASUREMENTS:
        assert dm["relation"] in ("same_channel", "superset_upper_bound", "subset_lower_bound", "lumped_total_vs_channel")
        assert dm["measured"]["kind"] in ("rate", "upper_bound", "lower_bound", "ea")
        assert set(dm["matches"]) == {"channels", "rate_law", "concentration", "conditions_stated"}
        assert dm.get("s4_route", "band") in ("band", "rival") and dm["sources"] and dm["cite"]
        if dm["relation"] == "subset_lower_bound":
            assert dm["measured"]["kind"] == "lower_bound"
        if dm["relation"] == "superset_upper_bound":
            assert dm["measured"]["kind"] == "upper_bound"
    assert len({s["id"] for s in et.STRUCTURAL_CHOICES}) == len(et.STRUCTURAL_CHOICES)
    for tier, _, sides in et.BAND_ROW_TIERS.values():
        assert tier in et.TIERS
        assert (sides is not None) == (tier == "bounded")


@pytest.fixture(scope="module")
def record():
    return et.build_record()


def test_every_shipped_coordinate_has_a_tier(record):
    from src.kinetic_core import engine
    ids = {r["id"] for r in record["rows"]}
    for lane in engine.MAILLARD_LANES:
        for key in engine.core_parameters(lane):
            assert f"{key}.log10_k_ref" in ids and f"{key}.ea_kj_mol" in ids
    assert all(r["tier"] in et.TIERS for r in record["rows"]), record["summary"]["unclassified"]
    assert len(ids) == len(record["rows"])


def test_known_readings(record):
    rows = {r["id"]: r for r in record["rows"]}
    mgo_ea = rows["k_mgo_mel.ea_kj_mol"]
    assert mgo_ea["tier"] == "fitted_unidentified"
    assert any(lc["label"] == "identified" for lc in mgo_ea["label_conflicts"])  # uncertainty.py's label, overridden
    assert rows["k_schiff.log10_k_ref"]["tier"] == "measured"
    assert rows["k_dimer_release.log10_k_ref"]["tier"] == "structural"  # inert arm (B17 not shipped)
    assert rows["k_nf_decay.log10_k_ref"]["identification"]["status"] == "unidentified"
    assert not any(r["tier"] == "fitted_identified" for r in record["rows"])  # no full profile on disk
    assert record["summary"]["s2_candidates"] == []
    assert rows["k_mgo_mel.log10_k_ref"]["contradicted_by"][0]["s2_status"].startswith("consistent")
    b3 = rows["k_acr_dp.log10_k_ref"]["fit_rows_withdrawn_or_corrected"]
    assert b3 and b3[0].startswith("B3: spi_extrusion_130C_acrylamide_ppb")
    assert {w["fit"] for w in record["withdrawn_or_corrected_fit_rows"]} >= {"B3", "B20"}


def test_search_bounds_are_never_a_band(record):
    for e in record["summary"]["s4_sensitive_unidentified"]:
        row = next(r for r in record["rows"] if r["id"] == e["coordinate"])
        if row.get("band_is_search_bound"):
            assert e["band_width_dex"] is None


def test_unidentified_lists_reconciled(record):
    rec = record["summary"]["unidentified_lists_reconciled"]
    assert rec["on_a_list_but_not_unidentified_here"] == []
    assert set(rec["unpinned_constants_yml"]) <= set(rec["sulfur_evidence_unpinned_that_matter"])
