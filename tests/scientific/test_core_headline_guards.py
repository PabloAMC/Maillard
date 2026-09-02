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
DIRECTIONAL = data_paths.CORE_DIRECTIONAL_SCORES


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
    # 2026-09-03: the xylose pH-5 bundle moved from the hold-out to the trust loop (in-fit)
    assert {k: v["benchmarks"] for k, v in s["by_panel"].items()} == {
        "trust_loop": 20, "maillard_path_holdout": 16, "external_matrix": 4,
    }


def test_core_evidence_roles_are_40_predictive_and_the_legacy_split_is_kept_beside(tracked_scores):
    """Pinned 2026-09-03 (B3). The core never fitted the rows the legacy matrix factors and
    projection budget were solved from, so those five bundles are predictive on the core and
    are listed as differing from the legacy role. Moving a bundle INTO fit_recovery without a
    core fit record that names it is the laundering route this guard blocks."""
    s = tracked_scores["summary"]
    assert s["evidence_role_totals"] == {"predictive": 40}


def test_within_3x_is_8_of_49_and_out_of_sample_4_of_40(tracked_scores):
    """Pinned 2026-09-03 (B3). honest_literature == all 49 rows (no fit-recovery, no synthetic
    on the core); out_of_sample removes the 9 rows the sulfur fit read (4 of which are hits)."""
    s = tracked_scores["summary"]
    assert (s["within_band_count"], s["matched_compound_count"]) == (8, 49)
    assert (s["honest_literature"]["within_band"], s["honest_literature"]["rows"]) == (8, 49)
    assert (s["out_of_sample"]["within_band"], s["out_of_sample"]["rows"]) == (4, 40)
    assert (s["in_core_fit"]["within_band"], s["in_core_fit"]["rows"]) == (4, 9)
    assert s["holdout_within_band"] == {"hits": 2, "total": 28}
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


def test_the_xylose_row_the_sulfur_fit_read_is_declared_in_fit(tracked_scores):
    """Pinned 2026-09-03. Found at B3 on the hold-out panel; moved to the trust loop the same
    day (owner decision) and declared in-fit by the generated fit-target record."""
    xyl = next(
        b for b in tracked_scores["benchmarks"]
        if b["benchmark_id"] == "hofmann1998_xylose_cysteine_145C_20min_pH5"
    )
    # moved 2026-09-03: the fit read it, so it lives on the trust loop, declared in-fit
    assert xyl["panel"] == "trust_loop"
    assert xyl["fit_target_of"]["core"] == ["kinetic_core_b8_fit_targets.json"]
    assert all(r["in_core_fit"] for r in xyl["compounds"])


# --------------------------------------------------------------------------------------
# 3. The envelope
# --------------------------------------------------------------------------------------


def test_core_envelope_covers_11_of_44_evaluable_literature_rows_and_8_of_35_out_of_sample():
    """Pinned 2026-09-03 (B8, n=200 seed 0). The sulfur lane is now sampled jointly from the
    Laplace covariance at its B8 optimum (18 of 23 free coordinates identified, reduced
    chi-square 1.03), so the 24 rows that were "not evaluable" at B3 are evaluable: 7/25 ->
    11/44, median width 0.94 -> 1.07 dex, 5 rows still not evaluable (zero-width or refused
    draws). Out of sample (every row the sulfur fit read removed): 8/35."""
    payload = json.loads(ENVELOPE.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert (s["n_samples"], s["seed"]) == (200, 0)
    lit = s["honest_literature_coverage"]
    assert (lit["hits"], lit["total"], lit["not_evaluable"]) == (11, 44, 5)
    assert lit["median_ci_width_log10"] == pytest.approx(1.0687, abs=5e-4)
    oos = s["out_of_sample_literature_coverage"]
    assert (oos["hits"], oos["total"]) == (8, 35)
    assert s["unsampled_lanes"] == []
    assert s["sulfur_laplace"]["identified"] == 18 and s["sulfur_laplace"]["free"] == 23
    assert s["sulfur_laplace"]["reduced_chi_square"] == pytest.approx(1.03, abs=0.01)
    assert s["observable_multiplier_policy"]["rows_by_family"] == {
        "headspace": 2, "extraction": 37, "undeclared": 10,
    }
    readme = _doc_text(README)
    _assert_quoted(readme, "11 of 44", "README.md", "the core envelope's literature coverage")
    _assert_quoted(readme, "8 of 35", "README.md", "the core envelope's out-of-sample coverage")
    _assert_quoted(readme, "18 of 23", "README.md", "the identified sulfur coordinates")


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


# --------------------------------------------------------------------------------------
# 5. The directional claims panel, on the core (step B7)
# --------------------------------------------------------------------------------------


def test_core_scores_18_of_27_independent_directional_claims():
    """Pinned 2026-09-03 (B7). 69 claims; 16 are prose-only and 22 independent claims are
    not evaluable on the core (refused arms: 2,5-dimethylpyrazine, 2-pentylfuran, nonanal from
    oleate, H2S / hydroxyacetaldehyde precursors; a zero disulfide). Of the evaluable
    independent claims the core agrees on 18 of 30 -- 15 of 22 with pH and water activity
    set aside, 3 of 8 on those two axes, 3 of the independent misses (7 over all claims) being
    identical predictions across an axis the lane carries no term for. The retired lane's 24/36 was a different engine
    on a different evaluable subset and is NOT a baseline for this number."""
    payload = json.loads(DIRECTIONAL.read_text(encoding="utf-8"))
    s = payload["summary"]
    assert payload["panel"]["claims"] == 69
    # RE-PINNED 2026-09-03 (step 5): comparisons that move an axis the lane has no term for
    # are REFUSED by the engine (water activity everywhere; pH on trunk / acrylamide / lipid),
    # so those claims are not evaluable instead of identical-prediction misses: 18/30 -> 18/27.
    assert s["headline"] == [18, 27]
    ind = s["independent"]
    assert (ind["excluding_ph_aw"]["agree"], ind["excluding_ph_aw"]["evaluable"]) == (15, 22)
    assert (ind["ph_aw"]["agree"], ind["ph_aw"]["evaluable"]) == (3, 5)
    assert ind["total"]["not_evaluable"] == 25
    assert ind["total"]["mechanism_absent"] == 0
    assert s["not_evaluable_reasons"]["refused by the engine"] >= 4
    cats = {k: (v["agree"], v["evaluable"]) for k, v in ind["by_category"].items()}
    assert cats["sugar_identity"] == (5, 9)
    assert cats["temperature"] == (5, 7)
    assert cats["ph"] == (3, 5)
    assert cats["moisture_aw"] == (0, 0)   # every a_w comparison is refused
    assert cats["additive_cysteine"] == (3, 3)
    assert cats["time"] == (2, 2)
    readme = _doc_text(README)
    _assert_quoted(readme, "18 of 27", "README.md", "the core's directional headline")
    _assert_quoted(readme, "15 of 22", "README.md", "the directional count excluding pH and water activity")
    _assert_quoted(readme, "3 of 5", "README.md", "the directional count on pH")


def test_directional_scorecard_is_not_stale():
    """The tracked artifact must be what the code scores today (~90 s)."""
    from src.kinetic_core import directional

    tracked = json.loads(DIRECTIONAL.read_text(encoding="utf-8"))
    live = directional.score_panel()
    assert live["summary"]["headline"] == tracked["summary"]["headline"]
    assert {r["claim_id"]: r["status"] for r in live["claims"]} == {
        r["claim_id"]: r["status"] for r in tracked["claims"]
    }
