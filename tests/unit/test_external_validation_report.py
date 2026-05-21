"""Tests for Lane A.3 external validation reporting."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.external_validation import (
    build_external_validation_report,
    get_holdout_benchmark_files,
    render_external_validation_markdown,
    write_external_validation_artifact,
)


@pytest.fixture(scope="module")
def report_payload():
    return build_external_validation_report(n_samples=8, seed=0)


def test_external_validation_report_matches_holdout_panel(report_payload):
    summary = report_payload["summary"]
    assert summary["holdout_benchmark_count"] == len(get_holdout_benchmark_files())
    assert summary["matched_compound_count"] >= 8
    assert summary["ci_coverage_rate"] is not None
    assert summary["median_abs_log10_error"] is not None
    assert summary["median_accuracy_fold"] is not None
    benchmark_ids = {row["benchmark_id"] for row in report_payload["benchmarks"]}
    assert "external_validation_liu_2023_ppi_offnote_baseline" in benchmark_ids
    assert "external_validation_li_2026_spi_wg_hme_control" in benchmark_ids


def test_external_validation_markdown_and_json_write(report_payload, tmp_path: Path):
    text = render_external_validation_markdown(report_payload)
    assert text.startswith("<!-- Auto-regenerated")
    assert "Headline external trust metric" in text
    paths = write_external_validation_artifact(report_payload, output_dir=tmp_path, basename="ext_val_report_test")
    assert paths["md"].exists()
    assert paths["json"].exists()
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert payload["summary"]["holdout_benchmark_count"] == report_payload["summary"]["holdout_benchmark_count"]


def test_lane_f_failing_compounds_flagged(report_payload):
    """Lane F (sprint 2026-05-10b): hexanal/nonanal should land in
    `external_failing_compounds` (mean |log10 error| > 1 dex)."""
    failing = {
        str(row["compound"]).strip().lower()
        for row in report_payload.get("external_failing_compounds", []) or []
    }
    assert "hexanal" in failing
    assert "nonanal" in failing
    assert report_payload.get("external_failing_threshold_dex", 0.0) > 0.0


def test_lane_f_sidecar_is_written(report_payload, tmp_path: Path):
    paths = write_external_validation_artifact(report_payload, output_dir=tmp_path)
    assert "failing_sidecar" in paths
    sidecar = paths["failing_sidecar"]
    assert sidecar.exists()
    sidecar_payload = json.loads(sidecar.read_text(encoding="utf-8"))
    failing = {
        str(row["compound"]).strip().lower()
        for row in sidecar_payload.get("external_failing_compounds", []) or []
    }
    assert "hexanal" in failing
