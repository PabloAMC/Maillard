"""Unit tests for src.uncertainty_propagation."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from src.uncertainty_propagation import (
    DEFAULT_FAMILY_PRIORS,
    ParameterPrior,
    _barrier_offsets,
    _percentile,
    default_priors,
    propagate_benchmarks,
    render_markdown,
    sample_offset_vectors,
    write_artifact,
)


def test_default_priors_cover_known_families():
    priors = default_priors()
    keys = {p.key for p in priors}
    # Every default family appears.
    assert keys == set(DEFAULT_FAMILY_PRIORS)
    for p in priors:
        assert p.sigma_kcal > 0.0


def test_sample_offset_vectors_reproducible_and_zero_mean():
    priors = [ParameterPrior(key="schiff_condensation", sigma_kcal=2.0)]
    a = sample_offset_vectors(priors, n=64, seed=42)
    b = sample_offset_vectors(priors, n=64, seed=42)
    assert a == b
    flat = [s["schiff_condensation"] for s in a]
    assert abs(sum(flat) / len(flat)) < 0.6  # near zero with sigma=2, n=64


def test_barrier_offsets_context_restores_environ():
    sentinel = "{\"original\": 1.0}"
    os.environ["BARRIER_OFFSETS"] = sentinel
    try:
        with _barrier_offsets({"schiff_condensation": 1.5}):
            assert json.loads(os.environ["BARRIER_OFFSETS"]) == {"schiff_condensation": 1.5}
        assert os.environ["BARRIER_OFFSETS"] == sentinel
    finally:
        os.environ.pop("BARRIER_OFFSETS", None)
    with _barrier_offsets({"amadori_rearrangement": -0.5}):
        assert "BARRIER_OFFSETS" in os.environ
    assert "BARRIER_OFFSETS" not in os.environ


def test_percentile_basic():
    assert _percentile([1.0, 2.0, 3.0, 4.0, 5.0], 50.0) == 3.0
    assert _percentile([1.0, 2.0], 0.0) == 1.0
    assert _percentile([1.0, 2.0], 100.0) == 2.0


@pytest.mark.scientific_regression
def test_propagate_benchmarks_smoke(tmp_path: Path):
    """Quick MC smoke test: small N, single-prior. Verify schema + that
    coverage rate is within the unit interval and CI bounds are ordered."""

    priors = [ParameterPrior(key="schiff_condensation", sigma_kcal=1.0)]
    payload = propagate_benchmarks(n_samples=12, seed=7, priors=priors)

    assert "summary" in payload and "benchmarks" in payload
    summary = payload["summary"]
    assert summary["n_samples"] == 12
    coverage = summary["ci_coverage_rate"]
    if coverage is not None:
        assert 0.0 <= coverage <= 1.0

    # At least one benchmark should produce envelopes.
    populated = [b for b in payload["benchmarks"] if b["compounds"]]
    assert populated, "Expected at least one benchmark with envelopes"

    for benchmark in populated:
        for compound in benchmark["compounds"]:
            assert compound["predicted_p5"] <= compound["predicted_p50"] <= compound["predicted_p95"]
            assert compound["predicted_std"] >= 0.0

    # Markdown render and artifact write should not raise.
    md = render_markdown(payload)
    assert "Prediction Uncertainty Envelope" in md
    paths = write_artifact(payload, output_dir=tmp_path, basename="smoke")
    assert paths["json"].exists() and paths["md"].exists()
