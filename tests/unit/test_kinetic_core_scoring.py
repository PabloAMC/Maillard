"""Retirement step B3: the core's own panel scorecard.

What is pinned here is the CONTRACT (fields, refusal vocabulary, the contract
thresholds resolved the way the legacy harness resolved them, agreement with the
cutover exam's points) -- not the numbers. The numbers are re-pinned at B4.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from src import data_paths
from src.kinetic_core import fit_targets, panel, scoring
from src.kinetic_core.engine import predict

ROOT = Path(__file__).resolve().parents[2]
MFT = "2-Methyl-3-furanthiol (MFT)"
FFT = "2-Furfurylthiol (FFT)"

_SMALL_PANEL = (
    (data_paths.BENCHMARKS_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json", "trust_loop"),
    (data_paths.BENCHMARKS_DIR / "acrylamide_spi_extrusion_130C_ACSRef3.json", "trust_loop"),
    (data_paths.BENCHMARKS_DIR / "hofmann1998_norfuraneol_h2s_145C_20min_pH5.json", "trust_loop"),
    (data_paths.BENCHMARKS_DIR / "pea_isolate_uht_140C_Trikusuma2019.json", "trust_loop"),
    (
        data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5.json",
        "maillard_path_holdout",
    ),
    (
        data_paths.MAILLARD_PATH_HOLDOUT_DIR / "mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7.json",
        "maillard_path_holdout",
    ),
)

BENCHMARK_FIELDS = {
    "benchmark_id", "bench_file", "panel", "tier", "family", "execution_path", "protein_type",
    "conditions", "signal_origin", "evidence_role", "fit_target_of",
    "core_fit_targets", "quantification_family", "scale_thresholds", "compounds",
    "refused_compounds", "matched_compounds", "total_compounds", "coverage", "max_ratio",
    "mean_abs_log10_error", "pearson_r_log10", "ranking_status", "scale_status",
    "overall_status", "strict_gate_eligible", "strict_ready", "blocking_issues",
    "within_band_count", "lanes",
}
ROW_FIELDS = {
    "compound", "target_unit", "measured", "measured_ppb", "predicted", "fold_error",
    "abs_log10_error", "within_band", "within_contract_ratio", "interval_ug_per_L",
    "measured_within_interval", "lane", "envelope_state", "shared_with", "in_core_fit",
}
SUMMARY_FIELDS = {
    "panel_benchmark_count", "scored_benchmark_count", "matched_compound_count",
    "refused_compound_count", "within_band_count", "within_band_rate", "fold_summary",
    "evidence_role_totals", "predictive_passes",
    "strict_ready", "by_panel", "by_evidence_role", "by_lane", "honest_literature",
    "out_of_sample", "in_core_fit", "holdout_within_band",
}


@pytest.fixture(scope="module")
def small():
    return scoring.score_panel(_SMALL_PANEL)


@pytest.fixture(scope="module")
def full():
    return scoring.score_panel()


# ---------------------------------------------------------------------------
# The metadata helpers have one engine-neutral home
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "path, ratio, log10, source",
    [
        ("hofmann1998_ribose_cysteine_145C_20min_pH5.json", 1.1, 0.0414, "bundle"),
        ("acrylamide_spi_extrusion_130C_ACSRef3.json", 1.5, 0.2, "bundle"),
        ("thiamine_cys_glucose_120C_Bolton1994.json", 3.0, 0.48, "bundle"),
        ("cys_ribose_140C_Hofmann1998.json", 1.5, 0.10, "default"),  # contract RETIRED -> free default
    ],
)
def test_scale_contract_is_the_bundles_own_else_the_global_default(path, ratio, log10, source):
    b = scoring.score_benchmark(data_paths.BENCHMARKS_DIR / path, "trust_loop")
    assert b["scale_thresholds"]["max_ratio"] == pytest.approx(ratio)
    assert b["scale_thresholds"]["mean_abs_log10_error"] == pytest.approx(log10)
    assert source in b["scale_thresholds"]["source"]


# ---------------------------------------------------------------------------
# The core fit-target ledger
# ---------------------------------------------------------------------------


_bundle_path = data_paths.bundle_path  # panel-wide resolver (2026-09-03)


def _fit_target_record():
    paths = fit_targets.record_paths()
    assert [p.name for p in paths] == ["kinetic_core_b9_fit_targets.json"]
    return json.loads(paths[0].read_text())


def test_the_generated_fit_target_record_names_real_rows_of_real_bundles():
    """The sulfur fit's level rows carry benchmark_id / benchmark_compound in the generator;
    generate_core_fit_targets.py writes them in the fit-target index's vocabulary. Every
    declared (bundle, compound) must exist, and a ppb anchor must quote the bundle's value."""
    record = _fit_target_record()
    assert record["artifact"] == "kinetic_core_fit_targets"
    assert record["fit_target_ids"] and record["rows"]
    for row in record["rows"]:
        bench = panel.load_bundle(_bundle_path(row["benchmark_id"]))
        targets = panel.bundle_targets(bench)
        assert row["compound"] in targets, (row["id"], row["benchmark_id"], sorted(targets))
        if row["kind"] == "conc":
            measured = panel.measured_value(bench, row["compound"], targets[row["compound"]])
            assert f"{measured:g} ppb" in row["anchor"], (row["id"], measured, row["anchor"])


def test_sulfur_fit_leverage_is_read_from_the_record_and_is_low():
    record = _fit_target_record()
    lev = record["fit_leverage"]
    assert lev["free_parameters"] < lev["fitted_rows"]
    report = json.loads((data_paths.VALIDATION_DIR / "kinetic_core_b9_fit_report.json").read_text())
    assert lev["free_parameters"] == report["free_set"]["n_free"]
    targets = fit_targets.core_fit_targets("hofmann1998_c2c3_recombination_145C_20min_pH5")
    assert len(targets) == 1 and targets[0].leverage_class == "global_low_leverage"


def test_in_core_fit_flags_only_the_step_level_bundles_after_b9():
    # B9 (2026-09-03) removed the Hofmann Table 1 level rows from the objective
    assert not fit_targets.in_core_fit("hofmann1998_ribose_cysteine_145C_20min_pH5", MFT)
    assert not fit_targets.in_core_fit("mp_holdout_hofmann1998_xylose_cysteine_145C_20min_pH5", FFT)
    assert fit_targets.in_core_fit("hofmann1998_c2c3_recombination_145C_20min_pH5", MFT)
    assert not fit_targets.in_core_fit("mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7", MFT)
    assert not fit_targets.in_core_fit("hofmann1998_ribose_cysteine_145C_20min_pH5", "furfural")


def test_core_evidence_role_ignores_legacy_fit_declarations_but_reports_them(small):
    by_id = {b["benchmark_id"]: b for b in small["benchmarks"]}
    trik = by_id["pea_isolate_uht_140C_Trikusuma2019"]
    # the legacy matrix observability factors were fitted to this bundle; the core's were not
    assert trik["evidence_role"] == "predictive"
    assert not trik["fit_target_of"]["core"]
    hof = by_id["hofmann1998_ribose_cysteine_145C_20min_pH5"]
    assert hof["evidence_role"] == "predictive"
    assert hof["fit_target_of"]["core"] == []  # B9: the level rows left the objective
    assert not any(r["in_core_fit"] for r in hof["compounds"])


# ---------------------------------------------------------------------------
# The scorecard contract
# ---------------------------------------------------------------------------


def test_artifact_contract_fields(small):
    assert SUMMARY_FIELDS <= set(small["summary"])
    assert small["parameter_sources"] and all(len(p["sha256"]) == 64 for p in small["parameter_sources"])
    for b in small["benchmarks"]:
        assert BENCHMARK_FIELDS <= set(b), b["benchmark_id"]
        assert b["overall_status"] in scoring.OVERALL_STATUSES
        assert b["matched_compounds"] + len(b["refused_compounds"]) == b["total_compounds"]
        for r in b["compounds"]:
            assert ROW_FIELDS <= set(r)
            if r["fold_error"] is not None:
                assert r["within_band"] == (r["fold_error"] <= small["pass_band_level"])
                assert r["within_contract_ratio"] == (r["fold_error"] <= b["scale_thresholds"]["max_ratio"])
        for r in b["refused_compounds"]:
            assert r["reason"].strip()


def test_predicted_point_is_the_deterministic_prediction(small):
    for b in small["benchmarks"]:
        spec = panel.core_spec(panel.load_bundle(ROOT / b["bench_file"]))
        for r in b["compounds"]:
            if r["target_unit"] == "ppb":
                run = predict(spec, [r["compound"]])
                assert r["predicted"] == pytest.approx(run.concentrations_ug_per_l[r["compound"]])


def test_refusals_use_the_engine_vocabulary_and_a_fully_refused_bundle_is_refused(small):
    by_id = {b["benchmark_id"]: b for b in small["benchmarks"]}
    h2s = by_id["hofmann1998_norfuraneol_h2s_145C_20min_pH5"]
    assert h2s["overall_status"] == "refused" and not h2s["compounds"]
    assert any("UNMAPPED PRECURSORS" in r["reason"] for r in h2s["refused_compounds"])
    assert h2s["blocking_issues"] and not h2s["strict_ready"]


def test_status_vocabulary_matches_the_legacy_harness(small):
    for b in small["benchmarks"]:
        if not b["compounds"]:
            continue
        if b["scale_status"] == "fail":
            assert b["overall_status"] in {"scale-gap", "coverage-gap", "ranking-gap"}
            assert any("max ratio" in i for i in b["blocking_issues"])
        if b["overall_status"] in {"pass", "pass-no-ranking"}:
            assert b["scale_status"] == "pass" and b["coverage"] == 1.0
        if b["strict_ready"]:
            assert b["tier"] == "PRIMARY" and b["execution_path"] == "free_precursor"


def test_shared_hofmann_rows_carry_shared_with(small):
    by_id = {b["benchmark_id"]: b for b in small["benchmarks"]}
    rows = {r["compound"]: r for r in by_id["mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7"]["compounds"]}
    assert rows[MFT]["shared_with"] == "hofmann_ribose_pH7_MFT"
    assert small["summary"]["shared_with_holdout_panel"] >= 1


def test_coverage_splits_are_consistent(small):
    s = small["summary"]
    assert s["honest_literature"]["rows"] == s["out_of_sample"]["rows"] + s["in_core_fit"]["rows"]
    assert s["honest_literature"]["within_band"] == (
        s["out_of_sample"]["within_band"] + s["in_core_fit"]["within_band"]
    )
    # B9 (2026-09-03): no Hofmann level row is a fit row any more and the small panel
    # holds none of the six step-level bundles, so the in-fit split is empty here.
    assert s["in_core_fit"]["rows"] == 0


def test_markdown_renders(small):
    text = scoring.render_markdown(small)
    for token in ("out-of-sample", "## Benchmarks", "## Rows", "hofmann1998_ribose_cysteine_145C_20min_pH5"):
        assert token in text


# ---------------------------------------------------------------------------
# Whole panel: every bundle, and agreement with the cutover exam's points
# ---------------------------------------------------------------------------


def test_every_panel_bundle_is_scored_or_refused_with_a_reason(full):
    scored, _ = panel.panel_bundles()
    assert {b["benchmark_id"] for b in full["benchmarks"]} == {
        panel.load_bundle(p)["benchmark_id"] for p, _ in scored
    }
    for b in full["benchmarks"]:
        targets = set(panel.bundle_targets(panel.load_bundle(ROOT / b["bench_file"])))
        assert {r["compound"] for r in b["compounds"]} | {r["compound"] for r in b["refused_compounds"]} == targets
    assert set(full["summary"]["by_panel"]) == {"trust_loop", "maillard_path_holdout", "external_matrix"}


