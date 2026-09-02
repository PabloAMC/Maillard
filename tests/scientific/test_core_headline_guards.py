"""Headline guards for the KINETIC CORE's own artifacts (retirement step B4, 2026-09-03).

The legacy guards in ``test_honest_headline_guards.py`` pin numbers the legacy harness
computes. These pin the numbers the core's scorecard and envelope publish, and the README
sentences that quote them, under the same rule: a moved number must move the prose in the
same change. Every pin carries the date and the cause of its last move.

Sources:
  * ``results/validation/core_panel_scores.json`` -- tracked, AND recomputed live here
    (~15 s) so a stale artifact fails rather than being echoed.
  * ``results/validation/core_prediction_uncertainty.json`` -- tracked (n=200, ~40 min to
    regenerate; not recomputed here).
  * ``results/validation/cutover_final_exam.json`` -- tracked.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import scoring

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
AUDIT = ROOT / "AUDIT.md"
SCORES = data_paths.VALIDATION_DIR / "core_panel_scores.json"
ENVELOPE = data_paths.VALIDATION_DIR / "core_prediction_uncertainty.json"
EXAM = data_paths.VALIDATION_DIR / "cutover_final_exam.json"


def _doc_text(path: Path) -> str:
    assert path.exists(), f"{path.name} is missing; the headline claims have no home"
    return path.read_text(encoding="utf-8")


def _assert_quoted(text: str, token: str, where: str, what: str) -> None:
    assert token in text, (
        f"{where} no longer quotes {token!r} for {what}. Either the number changed and this "
        f"guard was not updated with it, or the prose was edited away from the evidence."
    )


@pytest.fixture(scope="module")
def tracked_scores():
    return json.loads(SCORES.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def live_scores():
    return scoring.score_panel()


# --------------------------------------------------------------------------------------
# 1. The scorecard artifact is what the code produces today
# --------------------------------------------------------------------------------------


def test_tracked_scorecard_is_not_stale(tracked_scores, live_scores):
    """The cutover exam drifted 4/34 -> 3/34 silently before B2 regenerated it. Not again."""
    for key in ("panel_benchmark_count", "scored_benchmark_count", "matched_compound_count",
                "within_band_count", "predictive_passes", "strict_ready", "evidence_role_totals"):
        assert tracked_scores["summary"][key] == live_scores["summary"][key], key
    for key in ("honest_literature", "out_of_sample", "in_core_fit"):
        for field in ("rows", "within_band", "benchmarks"):
            assert tracked_scores["summary"][key][field] == live_scores["summary"][key][field], (key, field)
    live = {(b["benchmark_id"], r["compound"]): r["predicted"] for b in live_scores["benchmarks"] for r in b["compounds"]}
    tracked = {(b["benchmark_id"], r["compound"]): r["predicted"] for b in tracked_scores["benchmarks"] for r in b["compounds"]}
    assert set(live) == set(tracked)
    for key, value in live.items():
        assert tracked[key] == pytest.approx(value, rel=1e-9), key


# --------------------------------------------------------------------------------------
# 2. Panel, rows, band, strict-ready
# --------------------------------------------------------------------------------------


def test_core_panel_is_40_bundles_32_answered_49_rows(tracked_scores):
    """Pinned 2026-09-03 (B3). The union panel: 19 scored trust-loop bundles (the two
    *_Internal2026 synthetic snapshots are off the panel), 17 maillard_path hold-outs, 4
    external matrix bundles. 8 bundles are refused whole (H2S / hydroxyacetaldehyde /
    mercapto-2-propanone precursors, CML/CEL/furosine targets, 2-pentylfuran)."""
    s = tracked_scores["summary"]
    assert s["panel_benchmark_count"] == 40
    assert s["scored_benchmark_count"] == 32
    assert s["matched_compound_count"] == 49
    assert s["refused_compound_count"] == 18
    assert {k: v["benchmarks"] for k, v in s["by_panel"].items()} == {
        "trust_loop": 19, "maillard_path_holdout": 17, "external_matrix": 4,
    }


def test_core_evidence_roles_are_40_predictive_and_the_legacy_split_is_kept_beside(tracked_scores):
    """Pinned 2026-09-03 (B3). The core never fitted the rows the legacy matrix factors and
    projection budget were solved from, so those five bundles are predictive on the core and
    are listed as differing from the legacy role. Moving a bundle INTO fit_recovery without a
    core fit record that names it is the laundering route this guard blocks."""
    s = tracked_scores["summary"]
    assert s["evidence_role_totals"] == {"predictive": 40}
    assert s["evidence_role_differs_from_legacy"] == [
        "cys_ribose_140C_Hofmann1998",
        "hofmann1998_norfuraneol_h2s_145C_20min_pH5",
        "pea_isolate_40C_PratapSingh2021",
        "pea_isolate_uht_140C_Trikusuma2019",
        "soy_isolate_40C_PratapSingh2021",
    ]


def test_within_3x_is_8_of_49_and_out_of_sample_4_of_40(tracked_scores):
    """Pinned 2026-09-03 (B3). honest_literature == all 49 rows (no fit-recovery, no synthetic
    on the core); out_of_sample removes the 9 rows the sulfur fit read (4 of which are hits)."""
    s = tracked_scores["summary"]
    assert (s["within_band_count"], s["matched_compound_count"]) == (8, 49)
    assert (s["honest_literature"]["within_band"], s["honest_literature"]["rows"]) == (8, 49)
    assert (s["out_of_sample"]["within_band"], s["out_of_sample"]["rows"]) == (4, 40)
    assert (s["in_core_fit"]["within_band"], s["in_core_fit"]["rows"]) == (4, 9)
    assert s["holdout_within_band"] == {"hits": 3, "total": 30}
    readme = _doc_text(README)
    _assert_quoted(readme, "8 of 49", "README.md", "the core's within-3x count")
    _assert_quoted(readme, "4 of 40", "README.md", "the core's out-of-sample count")


def test_one_core_benchmark_is_strict_ready_and_it_is_bolton_1994(tracked_scores):
    """Pinned 2026-09-03 (B3). Bolton 1994 (thiamine + cysteine + glucose, 120 C): MFT 13 ->
    17.4 ppb, 1.34x, under the bundle's OWN 3x / 0.48-dex contract; PRIMARY, free_precursor,
    not a row any core fit read. The legacy lane has 0/23. A second strict-ready bundle is
    either a real result or a contract that was loosened: check the bundle's
    validation_contract before re-pinning."""
    s = tracked_scores["summary"]
    assert s["strict_ready"] == ["thiamine_cys_glucose_120C_Bolton1994"]
    assert s["predictive_passes"] == ["thiamine_cys_glucose_120C_Bolton1994"]
    bolton = next(b for b in tracked_scores["benchmarks"] if b["benchmark_id"] == "thiamine_cys_glucose_120C_Bolton1994")
    assert bolton["scale_thresholds"]["max_ratio"] == pytest.approx(3.0)
    assert not any(r["in_core_fit"] for r in bolton["compounds"])
    readme = _doc_text(README)
    _assert_quoted(readme, "1 of 40", "README.md", "the core strict-ready count")
    _assert_quoted(readme, "thiamine_cys_glucose_120C_Bolton1994", "README.md", "the strict-ready bundle")
    _assert_quoted(_doc_text(AUDIT), "strict-ready 1/40", "AUDIT.md", "the core strict-ready count")


def test_the_sulfur_fit_read_the_xylose_holdout_row(tracked_scores):
    """Pinned 2026-09-03 (B3). Declared in src/kinetic_core/fit_targets.py; the row must stay
    flagged until it is either moved off the hold-out panel or dropped from the fit."""
    xyl = next(
        b for b in tracked_scores["benchmarks"]
        if b["benchmark_id"] == "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5"
    )
    assert xyl["panel"] == "maillard_path_holdout"
    assert xyl["fit_target_of"]["core"] == ["kinetic_core_b8_fit_report.json"]
    assert all(r["in_core_fit"] for r in xyl["compounds"])


# --------------------------------------------------------------------------------------
# 3. The envelope
# --------------------------------------------------------------------------------------


def test_core_envelope_covers_7_of_25_evaluable_literature_rows_and_6_of_23_out_of_sample():
    """Pinned 2026-09-03 (B3 re-run of the B2 envelope, n=200 seed 0). 24 rows are not
    evaluable: sulfur rows quantified by extraction, on a lane whose fit report carries no
    uncertainty. Coverage 5/20 at B2 -> 7/25 because the five legacy fit-recovery rows are
    predictive on the core. A sulfur Laplace covariance will move ALL of these."""
    payload = json.loads(ENVELOPE.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert (s["n_samples"], s["seed"]) == (200, 0)
    lit = s["honest_literature_coverage"]
    assert (lit["hits"], lit["total"], lit["not_evaluable"]) == (7, 25, 24)
    assert lit["median_ci_width_log10"] == pytest.approx(0.9407, abs=5e-4)
    oos = s["out_of_sample_literature_coverage"]
    assert (oos["hits"], oos["total"]) == (6, 23)
    assert s["unsampled_lanes"] == ["sulfur"]
    assert s["observable_multiplier_policy"]["rows_by_family"] == {
        "headspace": 2, "extraction": 37, "undeclared": 10,
    }
    readme = _doc_text(README)
    _assert_quoted(readme, "7 of 25", "README.md", "the core envelope's literature coverage")
    _assert_quoted(readme, "6 of 23", "README.md", "the core envelope's out-of-sample coverage")
    _assert_quoted(readme, "24 sulfur rows not evaluable", "README.md", "the not-evaluable count")


# --------------------------------------------------------------------------------------
# 4. The exam, regenerated 2026-09-03 after it had drifted silently
# --------------------------------------------------------------------------------------


def test_cutover_exam_is_34_answered_3_within_band_paired_24_78x_vs_10_86x():
    """Pinned 2026-09-03. The tracked artifact said 4/34 while the code produced 3/34 (four
    MFT/FFT rows moved with the B8 sulfur refit); the README said 42.23x paired (the B7
    scoring). Both now match the artifact the code writes."""
    payload = json.loads(EXAM.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert (s["core_answered"], s["core_declined"], s["core_within_band"]) == (34, 6, 3)
    assert s["old_within_band"] == 7
    assert s["paired_subset"]["core"]["median_fold_error"] == pytest.approx(24.78, abs=0.01)
    assert s["paired_subset"]["old"]["median_fold_error"] == pytest.approx(10.86, abs=0.01)
    assert s["core"]["median_fold_error"] == pytest.approx(19.08, abs=0.01)
    readme = _doc_text(README)
    _assert_quoted(readme, "**3 / 34**", "README.md", "the exam's within-band count")
    _assert_quoted(readme, "**24.78x**", "README.md", "the exam's paired core median")
    _assert_quoted(readme, "**10.86x**", "README.md", "the exam's paired old-lane median")
    _assert_quoted(readme, "**19.08x**", "README.md", "the exam's all-answered core median")
