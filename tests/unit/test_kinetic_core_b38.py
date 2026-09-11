"""
WAVE B38 (2026-09-11): the identifiability audit. Pre-registration:
results/validation/kinetic_core_b38_prereg.md. The artifact is derived from the shipped fit
reports and moves no constant; these tests hold its shape and its accounting.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "results" / "validation" / "kinetic_core_b38_identifiability.json"
EXPECTED_FITS = ("B8", "B3", "B18", "B20", "B21")


def _payload():
    return json.loads(ART.read_text())


def test_the_audit_covers_every_shipped_fit():
    fits = [f["fit"].split(" ")[0] for f in _payload()["fits"]]
    assert tuple(fits) == EXPECTED_FITS, fits


def test_verdict_counts_add_up_per_fit_and_in_total():
    p = _payload()
    total = 0
    for f in p["fits"]:
        n = sum(f["verdict_counts"].values())
        assert n == f["n_free"] == len(f["coordinates"]), f["fit"]
        assert f["not_pinned"] == f["verdict_counts"]["WEAK"] + f["verdict_counts"]["UNIDENTIFIED"], f["fit"]
        assert 0 <= f["collinear_not_insensitive"] <= f["not_pinned"], f["fit"]
        total += n
    assert total == p["summary"]["total_free_coordinates"]


def test_every_prediction_has_a_verdict_and_the_thresholds_are_the_preregistered_ones():
    p = _payload()
    assert set(p["predictions"]) == {f"P{i}_" + k for i, k in enumerate(
        ("fewer_than_40pct_pinned", "b3_three_null_and_competitors_insensitive",
         "b8_thiol_sink_at_bound_plus_one_more_barrier", "collinearity_dominant",
         "priors_mostly_prior_dominated", "b18_barriers_pinned_by_band", "nothing_moves"), start=1)}
    for v in p["predictions"].values():
        assert isinstance(v["held"], bool)
    th = p["thresholds"]
    assert th["pinned_ci95_halfwidth"]["log10k"] == 0.5 and th["pinned_ci95_halfwidth"]["ea"] == 30.0
    assert th["weak_ci95_halfwidth"]["log10k"] == 1.5 and th["null_eigenvalue_fraction"] == 1e-3


def test_a_coordinate_on_a_bound_is_never_called_pinned():
    for f in _payload()["fits"]:
        for r in f["coordinates"]:
            if r["at_bound"]:
                assert r["verdict"] == "AT_BOUND", (f["fit"], r["key"])


def test_the_prior_cross_classifies_every_matched_prior():
    c = _payload()["priors_cross"]
    assert c["prior_dominated"] + c["data_dominated"] == c["n_priors_on_fitted_coordinates"] == len(c["rows"])
    assert all(r["dominance"] in ("prior-dominated", "data-dominated") for r in c["rows"])


def test_the_artifact_declares_its_inputs_for_the_freshness_gate():
    inputs = {i["path"] for i in _payload()["provenance"]["inputs"]}
    for stem in ("kinetic_core_b38_prereg.md", "kinetic_core_b8_laplace_covariance.json", "kinetic_core_b3_fit_report.json"):
        assert any(s.endswith(stem) for s in inputs), stem
