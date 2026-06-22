"""S20.1 — Monte Carlo uncertainty propagation through the benchmark surface.

Reuses the BARRIER_OFFSETS environment-variable hook already understood by
`src.barrier_constants.get_barrier()` so we do not duplicate evaluator logic.

Per benchmark / matched compound we report P5 / P50 / P95 predicted ppb plus
whether the measured value lies inside the 90% CI. The headline metric is the
fraction of matched compounds across the panel whose measured value sits
inside their MC envelope: a stronger contract than the existing 1.5x-band
ratio test because it integrates the effect of the prior uncertainty on each
upstream barrier.

The propagation is intentionally cheap (single-shot, family-uncorrelated
Gaussian offsets) so it can be re-run inside CI; it is *not* a Sobol sweep
and does not attempt to capture covariance between calibrated families.
"""

from __future__ import annotations

import json
import math
import os
import random
import re
import statistics
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from src.benchmark_types import BenchmarkEvaluation
from src.benchmark_validation import (
    DEFAULT_TARGET_TAG,
    evaluate_benchmark,
    get_benchmark_files,
    get_benchmark_metadata,
    load_benchmark,
)


ROOT = Path(__file__).resolve().parents[1]
QM_PROVENANCE_PATH = ROOT / "results" / "validation" / "qm_barrier_provenance.json"
PREDICTION_UNCERTAINTY_PATH = ROOT / "results" / "validation" / "prediction_uncertainty.json"

# Default Gaussian sigma (kcal/mol) per barrier-family offset key. The keys
# match the ones already accepted by BARRIER_OFFSETS via barrier_constants.
# Sigmas come from `get_barrier()` returning ±3.5 kcal/mol for FAST_BARRIERS
# entries and ±5.0 for the heuristic default, divided by 2 so ±2σ ≈ ±band.
DEFAULT_FAMILY_PRIORS: Dict[str, float] = {
    "schiff_condensation": 1.5,
    "amadori_rearrangement": 1.5,
    "1,2-enolisation": 1.5,
    "2,3-enolisation": 1.5,
    "strecker_degradation": 1.5,
    "cysteine_thermolysis": 1.5,
    "thiol_addition": 1.5,
    "thiol_addition_hexose": 1.5,
    "thiol_oxidation": 1.5,
    "aminoketone_condensation": 1.5,
    "retro_aldol": 1.5,
    "dehydration": 1.5,
    "beta_elimination": 1.5,
    "lipid_thiazole": 1.5,
}

# Default observable-multiplier priors. `sigma` is the natural-log standard
# deviation of a lognormal multiplier (mean=1.0). 0.3 ≈ ±35% one-sigma —
# wide enough to register matrix/headspace/retention uncertainty as a
# competing source of variance, narrow enough to keep the MC stable.
DEFAULT_OBSERVABLE_PRIORS: Dict[str, float] = {
    "matrix_headspace": 0.30,
    "henry_kaw": 0.30,
    "matrix_retention": 0.20,
}

# S27 Workstream B — observable priors for UNCALIBRATED matrix process-states.
# When a matrix prediction is NOT pinned by the per-(matrix, process_state)
# calibration registry (e.g. roasting, HME extrusion, any novel formulation), the
# headspace/retention translation is genuinely uncertain by ~1-2 orders of
# magnitude: protein-volatile binding, denaturation, and volatile stripping under
# severe processing all vary widely and are unmeasured for these conditions. These
# sigmas are a STRUCTURAL IGNORANCE PRIOR stating that absence-of-calibration =>
# large uncertainty. They are deliberately set from physical reasoning, NOT fitted
# to the external hold-out values (the 4 hold-out bundles remain frozen); the
# resulting coverage is reported as-is, not tuned to a target.
DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS: Dict[str, float] = {
    "matrix_headspace": 2.0,   # ln-sigma; 90% CI ~ +/- 1.43 dex (~27x)
    "henry_kaw": 1.0,
    "matrix_retention": 0.7,
}

# Multiplicative 90% CI band implied by the dominant uncalibrated matrix sigma
# above (matrix_headspace ln-sigma 2.0). Used by the RUNTIME report path to widen
# a formulation's per-compound envelope when the formulation falls outside the
# calibration scope, so a user's novel/out-of-registry run shows honestly wide CIs
# instead of the tight in-panel band. ~[1/27, 27] around the point estimate.
_UNCALIBRATED_HEADSPACE_SIGMA = DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS["matrix_headspace"]
UNCALIBRATED_ENVELOPE_LOWER_RATIO = math.exp(-1.645 * _UNCALIBRATED_HEADSPACE_SIGMA)
UNCALIBRATED_ENVELOPE_UPPER_RATIO = math.exp(1.645 * _UNCALIBRATED_HEADSPACE_SIGMA)

# Hard clamps on the sampled multipliers — guard against absurd MC tails. The
# upper bounds are wide enough that the calibrated priors (sigma 0.2-0.3) never
# reach them in practice (so the in-panel run is unaffected), while the
# uncalibrated tier above can still express its intended width.
_OBSERVABLE_CLAMP: Dict[str, Tuple[float, float]] = {
    "matrix_headspace": (0.01, 100.0),
    "henry_kaw": (0.01, 100.0),
    "matrix_retention": (0.01, 1.0),  # retention factor itself is bounded [0.01, 1]
}


@dataclass(frozen=True)
class ParameterPrior:
    """A prior distribution over either a barrier-family offset or an
    observable multiplier.

    For `kind == "barrier"`: `key` matches the BARRIER_OFFSETS dictionary
    consumed by `src.barrier_constants.get_barrier()`; `sigma_kcal` is the
    standard deviation of an additive Gaussian offset (mean 0).

    For `kind == "observable"`: `key` is one of the observable multiplier
    names (matrix_headspace, henry_kaw, matrix_retention); `sigma_kcal`
    becomes the lognormal log-sigma of a multiplier with mean 1.0.
    """

    key: str
    sigma_kcal: float
    source: str = "barrier_constants_default"
    kind: str = "barrier"


@dataclass(frozen=True)
class CompoundEnvelope:
    benchmark_id: str
    compound: str
    measured_ppb: float
    predicted_p5: float
    predicted_p50: float
    predicted_p95: float
    predicted_mean: float
    predicted_std: float
    inside_ci: bool


@dataclass(frozen=True)
class BenchmarkEnvelope:
    benchmark_id: str
    bench_file: str
    execution_path: str
    matched_compounds: int
    coverage_rate: float
    compounds: List[CompoundEnvelope] = field(default_factory=list)


def normalize_compound_key(name: str) -> str:
    if not name:
        return ""
    stripped = re.sub(r"\([^)]*\)", " ", str(name))
    stripped = stripped.lower()
    return re.sub(r"[^a-z0-9]+", " ", stripped).strip()


@lru_cache(maxsize=4)
def _load_prediction_interval_library(path_str: str) -> Dict[str, Any]:
    try:
        payload = json.loads(Path(path_str).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"summary": {"ci_level_pct": 90}, "by_compound": {}}

    by_compound: Dict[str, List[Dict[str, float]]] = {}
    for benchmark in payload.get("benchmarks", []) or []:
        for compound in benchmark.get("compounds", []) or []:
            compound_name = str(compound.get("compound", "")).strip()
            key = normalize_compound_key(compound_name)
            if not key:
                continue
            p5 = float(compound.get("predicted_p5", 0.0) or 0.0)
            p50 = float(compound.get("predicted_p50", 0.0) or 0.0)
            p95 = float(compound.get("predicted_p95", 0.0) or 0.0)
            if p50 <= 0.0 or not all(math.isfinite(value) for value in (p5, p50, p95)):
                continue
            by_compound.setdefault(key, []).append(
                {
                    "lower_ratio": max(0.0, p5 / p50),
                    "upper_ratio": max(1.0, p95 / p50),
                }
            )

    return {
        "summary": dict(payload.get("summary", {})),
        "by_compound": by_compound,
    }


def build_formulation_uncertainty_envelopes(
    predicted_ppb: Mapping[str, float],
    *,
    prediction_path: Path | str = PREDICTION_UNCERTAINTY_PATH,
    uncalibrated: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """Map cached benchmark envelopes onto current predicted compounds.

    The report path consumes the already-materialized trust-loop artifact rather
    than rerunning Monte Carlo propagation. Each formulation compound receives a
    benchmark-derived interval only when that compound already appears in the
    uncertainty panel.

    ``uncalibrated`` (S27 followup #1): when True, the formulation falls outside
    the calibration scope (out-of-registry matrix/process_state), so the tight
    in-panel band understates the true uncertainty. Each compound's envelope is
    widened to at least the structural-ignorance band (~[1/27, 27]x), matching the
    ``matrix_tier="uncalibrated"`` prior used for the external hold-out, and the
    envelope_source is tagged accordingly. In-scope formulations (default) are
    unchanged.
    """

    library = _load_prediction_interval_library(str(Path(prediction_path)))
    ci_level_pct = int(library.get("summary", {}).get("ci_level_pct", 90) or 90)
    by_compound = library.get("by_compound", {})

    envelopes: Dict[str, Dict[str, Any]] = {}
    for compound, raw_value in predicted_ppb.items():
        try:
            predicted_value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if predicted_value <= 0.0 or not math.isfinite(predicted_value):
            continue
        matches = by_compound.get(normalize_compound_key(str(compound)), [])
        if not matches:
            continue
        lower_ratio = statistics.median(row["lower_ratio"] for row in matches)
        upper_ratio = statistics.median(row["upper_ratio"] for row in matches)
        source = "prediction_uncertainty"
        if uncalibrated:
            # Out-of-calibration: never report a band tighter than the structural
            # ignorance prior. Widen (never narrow) in both directions.
            lower_ratio = min(lower_ratio, UNCALIBRATED_ENVELOPE_LOWER_RATIO)
            upper_ratio = max(upper_ratio, UNCALIBRATED_ENVELOPE_UPPER_RATIO)
            source = "prediction_uncertainty_uncalibrated"
        envelopes[str(compound)] = {
            "compound": str(compound),
            "predicted_ppb": predicted_value,
            "predicted_p5": predicted_value * lower_ratio,
            "predicted_p50": predicted_value,
            "predicted_p95": predicted_value * upper_ratio,
            "ci_level_pct": ci_level_pct,
            "support_count": len(matches),
            "envelope_source": source,
        }
    return envelopes


@contextmanager
def _barrier_offsets(offsets: Mapping[str, float]) -> Iterator[None]:
    previous = os.environ.get("BARRIER_OFFSETS")
    try:
        os.environ["BARRIER_OFFSETS"] = json.dumps(dict(offsets))
        yield
    finally:
        if previous is None:
            os.environ.pop("BARRIER_OFFSETS", None)
        else:
            os.environ["BARRIER_OFFSETS"] = previous


@contextmanager
def _observable_multipliers(
    *,
    matrix_headspace: float = 1.0,
    henry_kaw: float = 1.0,
    matrix_retention: float = 1.0,
) -> Iterator[None]:
    """Context-managed monkey-patches for matrix-headspace, Henry, and
    matrix-retention multipliers. The wrapped functions receive the same
    arguments and return ``original_value * multiplier`` (clamped where
    the underlying function clamps). Restored on exit.

    This is **diagnostic only**: it never persists, never writes back to
    calibration tables, and is never used in production prediction paths.
    """

    # Local imports keep the propagation module importable even if these
    # heavyweight evaluator modules fail to import in some lanes.
    from src import headspace as _headspace_module
    from src import matrix_correction as _matrix_module

    original_headspace = _headspace_module.HeadspaceModel.get_matrix_benchmark_headspace_factor
    original_kaw = _headspace_module.HeadspaceModel.get_kaw_at_temp
    original_retention = _matrix_module.resolve_compound_matrix_retention

    def _wrapped_headspace(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        return float(original_headspace(self, *args, **kwargs)) * matrix_headspace

    def _wrapped_kaw(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        return float(original_kaw(self, *args, **kwargs)) * henry_kaw

    def _wrapped_retention(*args, **kwargs):  # type: ignore[no-untyped-def]
        # resolve_compound_matrix_retention already clamps to [0.01, 1.0];
        # we re-clamp after multiplication so a perturbed retention factor
        # cannot leak above the physical bound.
        value = float(original_retention(*args, **kwargs)) * matrix_retention
        return max(0.01, min(1.0, value))

    try:
        _headspace_module.HeadspaceModel.get_matrix_benchmark_headspace_factor = _wrapped_headspace  # type: ignore[assignment]
        _headspace_module.HeadspaceModel.get_kaw_at_temp = _wrapped_kaw  # type: ignore[assignment]
        _matrix_module.resolve_compound_matrix_retention = _wrapped_retention  # type: ignore[assignment]
        yield
    finally:
        _headspace_module.HeadspaceModel.get_matrix_benchmark_headspace_factor = original_headspace  # type: ignore[assignment]
        _headspace_module.HeadspaceModel.get_kaw_at_temp = original_kaw  # type: ignore[assignment]
        _matrix_module.resolve_compound_matrix_retention = original_retention  # type: ignore[assignment]


def default_priors(
    *,
    include_observable: bool = True,
    matrix_tier: str = "calibrated",
) -> List[ParameterPrior]:
    """Build the default prior set, overriding sigmas where qm_barrier_provenance
    has narrowed the bound for a specific anchor target.

    qm_barrier_provenance currently catalogues per-target literature anchors
    rather than per-family runtime offsets, so we apply its evidence by
    *narrowing* (never widening) the family sigma if the provenance states a
    bounded-calibration tier. When ``include_observable`` is True we also
    add lognormal priors over the matrix-headspace, Henry, and matrix-
    retention multipliers (S20.4).

    ``matrix_tier`` (S27 Workstream B) selects the observable-multiplier sigmas:
    ``"calibrated"`` (default) uses the tight in-registry priors and is unchanged
    from prior behaviour; ``"uncalibrated"`` uses the wide structural-ignorance
    priors for matrix predictions that are not pinned by the calibration registry
    (the external hold-out is uncalibrated by construction). Barrier priors are
    identical in both tiers — only the matrix observables widen.
    """

    priors: Dict[str, ParameterPrior] = {
        key: ParameterPrior(
            key=key,
            sigma_kcal=sigma,
            source="barrier_constants_default",
            kind="barrier",
        )
        for key, sigma in DEFAULT_FAMILY_PRIORS.items()
    }

    _apply_qm_provenance_narrowing(priors)

    result: List[ParameterPrior] = list(priors.values())
    if include_observable:
        if matrix_tier == "uncalibrated":
            observable_priors = DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS
            source = "observable_uncalibrated_matrix"
        elif matrix_tier == "calibrated":
            observable_priors = DEFAULT_OBSERVABLE_PRIORS
            source = "observable_multiplier_default"
        else:
            raise ValueError(f"Unknown matrix_tier: {matrix_tier!r} (expected 'calibrated' or 'uncalibrated')")
        for key, sigma in observable_priors.items():
            result.append(
                ParameterPrior(
                    key=key,
                    sigma_kcal=sigma,
                    source=source,
                    kind="observable",
                )
            )
    return result


# Current-tier labels in qm_barrier_provenance.json that represent a genuine
# measured/DFT anchor strong enough to NARROW a core family's MC sigma. The file's
# other labels (literature_family_surrogate, literature_derived_transfer,
# family_rule_surrogate, no_literature_anchor) are weaker than the flat default and
# must NOT widen a well-anchored core family (see the design note below).
_QM_ANCHORED_TIERS = frozenset(
    {"bounded_calibration", "selective_dft_anchor", "dft_validated", "wet_lab_anchor"}
)


def _apply_qm_provenance_narrowing(priors: Dict[str, ParameterPrior]) -> None:
    """Narrow (never widen) core-family sigmas using qm_barrier_provenance.json.
    Mutates ``priors`` in place. Silently no-ops if the file is missing or
    malformed — the wider default sigma is the safe fallback.

    Design note (S27 followup #2): this provenance file catalogues the EXPLORATORY
    families (11-16: radical quench, quinone-Cys Michael, PE Schiff/Amadori,
    lysinoalanine, ascorbic dicarbonyl), which are NOT in the MC's core-14
    ``DEFAULT_FAMILY_PRIORS`` sampling set. Their tiers are surrogate / literature-
    derived. We therefore (a) match the provenance ``active_arrhenius_key`` EXACTLY
    against a core family key — never by substring — so an exploratory surrogate
    barrier (e.g. ``quinone_cys_michael_thiol_addition_family``) can no longer
    misattribute its uncertainty onto the well-anchored core ``thiol_addition`` used
    in benchmarked Cys+ribose chemistry; and (b) only ever NARROW, and only when the
    provenance ``current_tier`` is a genuine measured/DFT anchor. With the present
    file this correctly no-ops (no core family carries an anchored tier), preserving
    the in-panel headline; it activates automatically if a future entry anchors a
    core family. The earlier code read a non-existent ``entries`` key (the schema
    uses ``targets``) and was a silent no-op — this also fixes that latent bug.
    """
    if not QM_PROVENANCE_PATH.exists():
        return
    try:
        provenance = json.loads(QM_PROVENANCE_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return
    targets = provenance.get("targets") if isinstance(provenance, Mapping) else None
    if not isinstance(targets, list):
        return
    for entry in targets:
        if not isinstance(entry, Mapping):
            continue
        target_key = str(entry.get("active_arrhenius_key", ""))
        tier = str(entry.get("current_tier", ""))
        # Exact match only — no substring cross-mapping onto core families.
        if target_key not in priors:
            continue
        if tier in _QM_ANCHORED_TIERS:
            current = priors[target_key]
            if current.sigma_kcal > 1.0:
                priors[target_key] = ParameterPrior(
                    key=target_key,
                    sigma_kcal=1.0,
                    source="qm_barrier_provenance",
                    kind="barrier",
                )


def sample_offset_vectors(
    priors: Sequence[ParameterPrior],
    *,
    n: int,
    seed: int = 0,
) -> List[Dict[str, Dict[str, float]]]:
    """Independent samples per prior. Returns a list of
    ``{"barrier": {key: kcal_offset}, "observable": {key: multiplier}}``.

    Barrier priors yield additive Gaussian offsets (mean 0, sigma=sigma_kcal).
    Observable priors yield lognormal multipliers (mean 1.0, log-sigma=sigma_kcal),
    clamped to a per-key safe range.
    """

    rng = random.Random(seed)
    samples: List[Dict[str, Dict[str, float]]] = []
    for _ in range(n):
        barriers: Dict[str, float] = {}
        observables: Dict[str, float] = {}
        for prior in priors:
            if prior.kind == "observable":
                multiplier = math.exp(rng.gauss(0.0, prior.sigma_kcal))
                lo, hi = _OBSERVABLE_CLAMP.get(prior.key, (0.01, 100.0))
                observables[prior.key] = max(lo, min(hi, multiplier))
            else:
                barriers[prior.key] = rng.gauss(0.0, prior.sigma_kcal)
        samples.append({"barrier": barriers, "observable": observables})
    return samples


def _percentile(sorted_values: Sequence[float], pct: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return float(sorted_values[0])
    idx = (len(sorted_values) - 1) * (pct / 100.0)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return float(sorted_values[lo])
    frac = idx - lo
    return float(sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * frac)


def _baseline_targets(evaluation: BenchmarkEvaluation) -> List[Tuple[str, float]]:
    """Return (compound, measured_ppb) for matched comparisons in the
    baseline (zero-offset) evaluation. We only track compounds the
    evaluator already maps to a prediction at baseline; new matches
    appearing under perturbed offsets would not have a stable identity."""

    targets: List[Tuple[str, float]] = []
    for comparison in evaluation.comparisons:
        if comparison.matched_name is None:
            continue
        try:
            measured = float(comparison.measured_ppb)
        except (TypeError, ValueError):
            continue
        if measured <= 0.0 or not math.isfinite(measured):
            continue
        targets.append((comparison.compound, measured))
    return targets


def _comparison_predicted(evaluation: BenchmarkEvaluation, compound: str) -> Optional[float]:
    for comparison in evaluation.comparisons:
        if comparison.compound == compound and comparison.matched_name is not None:
            return float(comparison.predicted_ppb)
    return None


def propagate_benchmarks(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    n_samples: int = 200,
    seed: int = 0,
    priors: Optional[Sequence[ParameterPrior]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
    execution_paths: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Run Monte Carlo barrier-offset propagation across the benchmark panel.

    Returns a JSON-serialisable payload with per-compound P5/P50/P95 envelopes
    and an aggregate 90% CI coverage rate.
    """

    if priors is None:
        priors = default_priors()
    bench_files = list(benchmark_files) if benchmark_files is not None else list(get_benchmark_files())
    allowed_execution_paths = {
        str(path).strip()
        for path in (
            execution_paths
            if execution_paths is not None
            else ("free_precursor", "matrix_precursor_augmented")
        )
        if str(path).strip()
    }

    # Pre-load benchmarks and decide which compounds to track per benchmark
    # using the baseline (zero-offset) evaluation.
    bench_index: List[Dict[str, Any]] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in allowed_execution_paths:
            continue
        baseline = evaluate_benchmark(bench_file, target_tag=target_tag)
        targets = _baseline_targets(baseline)
        if not targets:
            continue
        bench_index.append(
            {
                "bench_file": Path(bench_file),
                "benchmark_id": str(bench.get("benchmark_id") or Path(bench_file).stem),
                "execution_path": metadata.execution_path,
                "targets": targets,
            }
        )

    samples = sample_offset_vectors(priors, n=n_samples, seed=seed)

    # accumulator[(bench_id, compound)] = list[predicted_ppb]
    accumulator: Dict[Tuple[str, str], List[float]] = {}

    for sample in samples:
        observable = sample.get("observable", {}) or {}
        with _barrier_offsets(sample.get("barrier", {}) or {}), _observable_multipliers(
            matrix_headspace=float(observable.get("matrix_headspace", 1.0)),
            henry_kaw=float(observable.get("henry_kaw", 1.0)),
            matrix_retention=float(observable.get("matrix_retention", 1.0)),
        ):
            for entry in bench_index:
                evaluation = evaluate_benchmark(entry["bench_file"], target_tag=target_tag)
                for compound, _measured in entry["targets"]:
                    predicted = _comparison_predicted(evaluation, compound)
                    if predicted is None or not math.isfinite(predicted):
                        continue
                    accumulator.setdefault((entry["benchmark_id"], compound), []).append(predicted)

    benchmark_envelopes: List[BenchmarkEnvelope] = []
    coverage_hits = 0
    coverage_total = 0

    for entry in bench_index:
        compound_envelopes: List[CompoundEnvelope] = []
        matched = 0
        for compound, measured in entry["targets"]:
            samples_for_compound = accumulator.get((entry["benchmark_id"], compound), [])
            if len(samples_for_compound) < max(5, n_samples // 4):
                # Too few finite samples to form a meaningful envelope.
                continue
            sorted_values = sorted(samples_for_compound)
            p5 = _percentile(sorted_values, 5.0)
            p50 = _percentile(sorted_values, 50.0)
            p95 = _percentile(sorted_values, 95.0)
            mean_val = statistics.fmean(samples_for_compound)
            std_val = statistics.pstdev(samples_for_compound) if len(samples_for_compound) > 1 else 0.0
            inside = bool(p5 <= measured <= p95)
            compound_envelopes.append(
                CompoundEnvelope(
                    benchmark_id=entry["benchmark_id"],
                    compound=compound,
                    measured_ppb=measured,
                    predicted_p5=p5,
                    predicted_p50=p50,
                    predicted_p95=p95,
                    predicted_mean=mean_val,
                    predicted_std=std_val,
                    inside_ci=inside,
                )
            )
            matched += 1
            coverage_total += 1
            if inside:
                coverage_hits += 1
        coverage_rate = (matched / max(len(entry["targets"]), 1)) if entry["targets"] else 0.0
        benchmark_envelopes.append(
            BenchmarkEnvelope(
                benchmark_id=entry["benchmark_id"],
                bench_file=str(entry["bench_file"]),
                execution_path=entry["execution_path"],
                matched_compounds=matched,
                coverage_rate=coverage_rate,
                compounds=compound_envelopes,
            )
        )

    payload: Dict[str, Any] = {
        "summary": {
            "n_samples": n_samples,
            "seed": seed,
            "benchmark_count": len(benchmark_envelopes),
            "matched_compound_count": coverage_total,
            "ci_coverage_hits": coverage_hits,
            "ci_coverage_rate": (coverage_hits / coverage_total) if coverage_total else None,
            "ci_level_pct": 90,
        },
        "priors": [
            {
                "key": p.key,
                "sigma_kcal": p.sigma_kcal,
                "source": p.source,
                "kind": p.kind,
            }
            for p in priors
        ],
        "benchmarks": [
            {
                "benchmark_id": env.benchmark_id,
                "bench_file": env.bench_file,
                "execution_path": env.execution_path,
                "matched_compounds": env.matched_compounds,
                "coverage_rate": env.coverage_rate,
                "compounds": [
                    {
                        "compound": c.compound,
                        "measured_ppb": c.measured_ppb,
                        "predicted_p5": c.predicted_p5,
                        "predicted_p50": c.predicted_p50,
                        "predicted_p95": c.predicted_p95,
                        "predicted_mean": c.predicted_mean,
                        "predicted_std": c.predicted_std,
                        "inside_ci": c.inside_ci,
                        "ci_width_log10": (
                            math.log10(c.predicted_p95 / c.predicted_p5)
                            if c.predicted_p5 > 0.0 and c.predicted_p95 > 0.0
                            else None
                        ),
                    }
                    for c in env.compounds
                ],
            }
            for env in benchmark_envelopes
        ],
    }
    return payload


def render_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    coverage_rate = summary.get("ci_coverage_rate")
    coverage_str = f"{coverage_rate * 100:.1f}%" if coverage_rate is not None else "n/a"
    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_prediction_uncertainty.py". Manual edits will be overwritten. -->',
        "",
        "# Prediction Uncertainty Envelope",
        "",
        "_Monte Carlo propagation of barrier-family offset priors (additive Gaussian, kcal/mol) through the benchmark evaluator. CI = 90% (P5–P95)._",
        "",
        f"**Headline trust metric**: measured value lies inside 90% CI for **{summary.get('ci_coverage_hits', 0)} / {summary.get('matched_compound_count', 0)}** matched compounds (**{coverage_str}**).",
        "",
        f"Samples per benchmark: {summary.get('n_samples', 0)}; seed {summary.get('seed', 0)}; benchmarks evaluated: {summary.get('benchmark_count', 0)}.",
        "",
        "## Priors",
        "",
        "| Key | Kind | σ (kcal/mol or log) | Source |",
        "| --- | --- | --- | --- |",
    ]
    for prior in payload.get("priors", []):
        lines.append(
            f"| `{prior.get('key')}` | {prior.get('kind', 'barrier')} | "
            f"{float(prior.get('sigma_kcal', 0.0)):.2f} | {prior.get('source')} |"
        )
    lines.append("")
    lines.append("## Per-benchmark envelopes")
    lines.append("")
    for benchmark in payload.get("benchmarks", []):
        lines.append(f"### `{benchmark.get('benchmark_id')}`")
        lines.append("")
        lines.append(
            f"- Execution path: `{benchmark.get('execution_path')}`"
        )
        lines.append(
            f"- Matched compounds with envelope: {benchmark.get('matched_compounds', 0)}"
        )
        lines.append("")
        lines.append("| Compound | Measured (ppb) | P5 | P50 | P95 | log₁₀ width | Inside 90% CI |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for compound in benchmark.get("compounds", []):
            inside = "✓" if compound.get("inside_ci") else "✗"
            width = compound.get("ci_width_log10")
            width_str = f"{float(width):.2f}" if width is not None else "n/a"
            lines.append(
                f"| {compound.get('compound')} | {float(compound.get('measured_ppb', 0.0)):.3g} | "
                f"{float(compound.get('predicted_p5', 0.0)):.3g} | {float(compound.get('predicted_p50', 0.0)):.3g} | "
                f"{float(compound.get('predicted_p95', 0.0)):.3g} | {width_str} | {inside} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_artifact(
    payload: Mapping[str, Any],
    *,
    output_dir: Path | str = ROOT / "results" / "validation",
    basename: str = "prediction_uncertainty",
) -> Dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{basename}.json"
    md_path = output_dir / f"{basename}.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    return {"json": json_path, "md": md_path}
