"""Unit tests for src.cross_validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from src.cross_validation import (
    _abs_log10_residual,
    compute_leverage,
    render_markdown,
    write_artifact,
)


def test_abs_log10_residual_handles_invalid_inputs():
    assert _abs_log10_residual(0.0, 1.0) is None
    assert _abs_log10_residual(1.0, 0.0) is None
    assert _abs_log10_residual(-1.0, 1.0) is None
    assert _abs_log10_residual(10.0, 100.0) == pytest.approx(1.0)


def test_compute_leverage_panel_and_benchmarks():
    payload = {
        "benchmarks": [
            {
                "benchmark_id": "all_in",
                "compounds": [
                    {"measured_ppb": 10.0, "predicted_p50": 10.0, "inside_ci": True},
                    {"measured_ppb": 20.0, "predicted_p50": 22.0, "inside_ci": True},
                ],
            },
            {
                "benchmark_id": "all_out",
                "compounds": [
                    {"measured_ppb": 100.0, "predicted_p50": 1.0, "inside_ci": False},
                ],
            },
        ]
    }
    report = compute_leverage(payload)
    panel = report["panel"]
    assert panel["matched_compound_count"] == 3
    assert panel["inside_ci_count"] == 2
    assert panel["coverage"] == pytest.approx(2 / 3)

    rows = {b["benchmark_id"]: b for b in report["benchmarks"]}
    # Removing all_in lowers coverage from 2/3 to 0/1, so leverage > 0.
    assert rows["all_in"]["leverage"] == pytest.approx(2 / 3 - 0.0)
    # Removing all_out raises coverage from 2/3 to 2/2, so leverage < 0.
    assert rows["all_out"]["leverage"] == pytest.approx(2 / 3 - 1.0)
    # Sorted descending by leverage so all_in comes first.
    assert report["benchmarks"][0]["benchmark_id"] == "all_in"


def test_write_artifact_round_trip(tmp_path: Path):
    payload = {
        "benchmarks": [
            {
                "benchmark_id": "b1",
                "compounds": [
                    {"measured_ppb": 1.0, "predicted_p50": 1.5, "inside_ci": True},
                ],
            }
        ]
    }
    report = compute_leverage(payload)
    md = render_markdown(report)
    assert "Leave-One-Benchmark-Out" in md
    paths = write_artifact(report, output_dir=tmp_path, basename="loo_smoke")
    assert paths["json"].exists() and paths["md"].exists()


def test_compute_leverage_handles_empty_payload():
    report = compute_leverage({"benchmarks": []})
    assert report["panel"]["matched_compound_count"] == 0
    assert report["panel"]["coverage"] is None
    assert report["benchmarks"] == []
