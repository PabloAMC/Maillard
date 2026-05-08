"""Unit tests for src.uncertainty_propagation."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from src.uncertainty_propagation import (
    DEFAULT_FAMILY_PRIORS,
    DEFAULT_OBSERVABLE_PRIORS,
    ParameterPrior,
    _barrier_offsets,
    _observable_multipliers,
    _percentile,
    default_priors,
    propagate_benchmarks,
    render_markdown,
    sample_offset_vectors,
    write_artifact,
)


def test_default_priors_cover_known_families_and_observables():
    priors = default_priors()
    barrier_keys = {p.key for p in priors if p.kind == "barrier"}
    observable_keys = {p.key for p in priors if p.kind == "observable"}
    assert barrier_keys == set(DEFAULT_FAMILY_PRIORS)
    assert observable_keys == set(DEFAULT_OBSERVABLE_PRIORS)
    for p in priors:
        assert p.sigma_kcal > 0.0


def test_default_priors_can_skip_observables():
    priors = default_priors(include_observable=False)
    assert all(p.kind == "barrier" for p in priors)


def test_sample_offset_vectors_reproducible_partitioned():
    priors = [
        ParameterPrior(key="schiff_condensation", sigma_kcal=2.0, kind="barrier"),
        ParameterPrior(key="matrix_headspace", sigma_kcal=0.3, kind="observable"),
    ]
    a = sample_offset_vectors(priors, n=64, seed=42)
    b = sample_offset_vectors(priors, n=64, seed=42)
    assert a == b
    assert all("barrier" in s and "observable" in s for s in a)
    barrier_flat = [s["barrier"]["schiff_condensation"] for s in a]
    obs_flat = [s["observable"]["matrix_headspace"] for s in a]
    # barrier offsets near zero, observable multipliers near 1.
    assert abs(sum(barrier_flat) / len(barrier_flat)) < 0.6
    assert all(0.1 <= v <= 10.0 for v in obs_flat)


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

    priors = [ParameterPrior(key="schiff_condensation", sigma_kcal=1.0, kind="barrier")]
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


def test_observable_multipliers_scale_headspace_factor():
    """Inside the context manager the matrix-headspace factor should be
    multiplied; outside, the original is restored."""
    from src.headspace import HeadspaceModel

    model = HeadspaceModel()
    baseline = float(
        model.get_matrix_benchmark_headspace_factor(
            "hexanal",
            protein_type="soy_iso",
            pH=6.5,
            temperature_celsius=100.0,
        )
    )
    with _observable_multipliers(matrix_headspace=2.0):
        scaled = float(
            model.get_matrix_benchmark_headspace_factor(
                "hexanal",
                protein_type="soy_iso",
                pH=6.5,
                temperature_celsius=100.0,
            )
        )
    restored = float(
        model.get_matrix_benchmark_headspace_factor(
            "hexanal",
            protein_type="soy_iso",
            pH=6.5,
            temperature_celsius=100.0,
        )
    )
    assert scaled == pytest.approx(baseline * 2.0)
    assert restored == pytest.approx(baseline)
