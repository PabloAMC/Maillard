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


# --- S27 Workstream B: process-state-aware (uncalibrated) matrix observable priors ---

from src.uncertainty_propagation import (  # noqa: E402
    DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS,
    _OBSERVABLE_CLAMP,
    sample_offset_vectors,
)


def _observables(priors):
    return {p.key: p.sigma_kcal for p in priors if p.kind == "observable"}


def test_calibrated_tier_is_unchanged_default():
    # The default tier must keep the original tight observable sigmas so the
    # in-panel prediction-uncertainty run is byte-identical.
    obs = _observables(default_priors())
    assert obs == dict(DEFAULT_OBSERVABLE_PRIORS)
    assert _observables(default_priors(matrix_tier="calibrated")) == dict(DEFAULT_OBSERVABLE_PRIORS)


def test_uncalibrated_tier_widens_only_observables():
    cal = default_priors(matrix_tier="calibrated")
    unc = default_priors(matrix_tier="uncalibrated")
    # Observables widen...
    assert _observables(unc) == dict(DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS)
    for key in DEFAULT_OBSERVABLE_PRIORS:
        assert _observables(unc)[key] > _observables(cal)[key]
    # ...barriers are identical between tiers.
    bar_cal = {p.key: p.sigma_kcal for p in cal if p.kind == "barrier"}
    bar_unc = {p.key: p.sigma_kcal for p in unc if p.kind == "barrier"}
    assert bar_cal == bar_unc


def test_unknown_matrix_tier_raises():
    with pytest.raises(ValueError):
        default_priors(matrix_tier="bogus")


def test_widened_clamp_is_inert_for_calibrated_sampling():
    # The clamp upper/lower bounds were widened to let the uncalibrated tier
    # express its width. This must NOT change calibrated sampling: at sigma 0.2-0.3
    # no sample should reach even the OLD bounds (0.1, 10.0) / (0.05, 1.0), so the
    # in-panel run is unaffected by the clamp change.
    samples = sample_offset_vectors(default_priors(matrix_tier="calibrated"), n=2000, seed=0)
    hs = [s["observable"]["matrix_headspace"] for s in samples]
    kaw = [s["observable"]["henry_kaw"] for s in samples]
    ret = [s["observable"]["matrix_retention"] for s in samples]
    assert min(hs) > 0.1 and max(hs) < 10.0
    assert min(kaw) > 0.1 and max(kaw) < 10.0
    assert min(ret) > 0.05 and max(ret) <= 1.0


def test_uncalibrated_sampling_actually_spans_wider_range():
    samples = sample_offset_vectors(default_priors(matrix_tier="uncalibrated"), n=2000, seed=0)
    hs = [s["observable"]["matrix_headspace"] for s in samples]
    # With ln-sigma 2.0 and the widened clamp, the headspace multiplier must reach
    # well beyond the old +/-10x band in at least one tail.
    assert max(hs) > 15.0 or min(hs) < 0.05
    # Retention stays physically bounded at <= 1.0.
    ret = [s["observable"]["matrix_retention"] for s in samples]
    assert max(ret) <= 1.0


# --- S27 followup #1: runtime process-state-aware envelope widening ---

from src.uncertainty_propagation import (  # noqa: E402
    PREDICTION_UNCERTAINTY_PATH,
    UNCALIBRATED_ENVELOPE_LOWER_RATIO,
    UNCALIBRATED_ENVELOPE_UPPER_RATIO,
    _load_prediction_interval_library,
    build_formulation_uncertainty_envelopes,
)


def _cached_library_keys(n=3):
    lib = _load_prediction_interval_library(str(PREDICTION_UNCERTAINTY_PATH))
    return list(lib.get("by_compound", {}).keys())[:n]


def test_calibrated_runtime_envelope_unchanged_by_default():
    keys = _cached_library_keys()
    if not keys:
        pytest.skip("no cached uncertainty library available")
    predicted = {k: 100.0 for k in keys}
    default = build_formulation_uncertainty_envelopes(predicted)
    explicit = build_formulation_uncertainty_envelopes(predicted, uncalibrated=False)
    assert default == explicit
    for env in default.values():
        assert env["envelope_source"] == "prediction_uncertainty"


def test_uncalibrated_runtime_envelope_widens_to_structural_floor():
    keys = _cached_library_keys()
    if not keys:
        pytest.skip("no cached uncertainty library available")
    predicted = {k: 100.0 for k in keys}
    cal = build_formulation_uncertainty_envelopes(predicted)
    unc = build_formulation_uncertainty_envelopes(predicted, uncalibrated=True)
    assert unc, "expected at least one cached-library compound to match"
    for name, env in unc.items():
        value = env["predicted_ppb"]
        # The band is never tighter than the structural-ignorance floor...
        assert env["predicted_p5"] <= value * UNCALIBRATED_ENVELOPE_LOWER_RATIO * (1 + 1e-9)
        assert env["predicted_p95"] >= value * UNCALIBRATED_ENVELOPE_UPPER_RATIO * (1 - 1e-9)
        assert env["envelope_source"] == "prediction_uncertainty_uncalibrated"
        # ...and strictly no narrower than the calibrated band in both directions.
        if name in cal:
            assert env["predicted_p5"] <= cal[name]["predicted_p5"]
            assert env["predicted_p95"] >= cal[name]["predicted_p95"]
        # p50 (the point estimate) is unchanged.
        assert env["predicted_p50"] == value


# --- S27 followup #2: provenance-tiered barrier sigma (narrow-only, exact-match) ---

def test_provenance_tiering_is_narrow_only_no_widening_on_current_file():
    # The provenance file catalogues only exploratory families (not in the MC set)
    # at surrogate tiers, so nothing should narrow AND nothing should widen: the
    # core barrier sigmas must equal the flat defaults. This guards the in-panel
    # headline against accidental barrier-sigma drift.
    barrier = {p.key: p.sigma_kcal for p in default_priors() if p.kind == "barrier"}
    assert barrier == dict(DEFAULT_FAMILY_PRIORS)


def test_provenance_narrowing_never_widens_a_core_family():
    from src.uncertainty_propagation import _apply_qm_provenance_narrowing

    priors = {
        key: ParameterPrior(key=key, sigma_kcal=sigma, source="x", kind="barrier")
        for key, sigma in DEFAULT_FAMILY_PRIORS.items()
    }
    _apply_qm_provenance_narrowing(priors)
    for key, sigma in DEFAULT_FAMILY_PRIORS.items():
        assert priors[key].sigma_kcal <= sigma  # narrow-only invariant
