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
from src.fit_target_index import fit_target_records_for
from src.benchmark_validation import (
    DEFAULT_TARGET_TAG,
    SYNTHETIC_BENCHMARK_ORIGINS,
    benchmark_evidence_role,
    benchmark_signal_origin,
    evaluate_benchmark,
    get_benchmark_files,
    get_benchmark_metadata,
    load_benchmark,
    matrix_source_anchor,
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
#
# Residual-based sizing (audit follow-up 2026-08-26, the derivation tasks/todo.md
# had promised): leave-lane-out transfer error over the in-panel matrix anchors
# (hold-out untouched) gives RMS ln-sigma = 2.86, 90% chi-square interval
# [1.98, 5.48] with n=6 — see scripts/generators/derive_matrix_sigma_from_residuals.py
# and results/validation/matrix_sigma_residual_derivation.md. The original 2.0 sat
# at the LOWER edge of that interval; on 2026-08-26 the owner approved raising the
# value to the residual-derived point estimate 2.86, so the tier is now sized from
# data rather than judgment.
DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS: Dict[str, float] = {
    "matrix_headspace": 2.86,  # ln-sigma from leave-lane-out residuals; 90% CI ~ +/- 2.04 dex (~110x)
    "henry_kaw": 1.0,
    "matrix_retention": 0.7,
}

# Multiplicative 90% CI band implied by the dominant uncalibrated matrix sigma
# above (auto-derived from the constant). Used by the RUNTIME report path to widen
# a formulation's per-compound envelope when the formulation falls outside the
# calibration scope, so a user's novel/out-of-registry run shows honestly wide CIs
# instead of the tight in-panel band. At ln-sigma 2.86: ~[1/110, 110].
_UNCALIBRATED_HEADSPACE_SIGMA = DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS["matrix_headspace"]
UNCALIBRATED_ENVELOPE_LOWER_RATIO = math.exp(-1.645 * _UNCALIBRATED_HEADSPACE_SIGMA)
UNCALIBRATED_ENVELOPE_UPPER_RATIO = math.exp(1.645 * _UNCALIBRATED_HEADSPACE_SIGMA)

# Hard clamps on the sampled multipliers — a guard against absurd MC tails.
#
# 2026-08-27 (Wave I) — THIS COMMENT USED TO BE FALSE, and the falsehood was load-bearing.
# It claimed the bounds were "wide enough that ... the uncalibrated tier can still express
# its intended width". At the sigma this tier shipped with since 2026-08-26 (2.86) that was
# simply not true:
#
#     advertised 90% CI  = exp(+/- 1.645 * 2.86) = +/- 110.4x
#     clamp              = +/- 100x
#     P(|multiplier| clamped) = 2 * P(Z > ln(100)/2.86) = 2 * P(Z > 1.610) = 10.7%
#
# So roughly one draw in nine was pinned to the clamp, the realised band was +-100x, not
# the +-110x quoted in the README and the reports, and the truncation fell hardest exactly
# where the uncalibrated tier is meant to be honest about not knowing. A "tail guard" that
# reshapes 10.7% of the distribution is not a guard; it is an undeclared prior.
#
# FIX: the clamp is now DERIVED FROM THE SIGMA rather than fixed, at
# `_OBSERVABLE_CLAMP_SIGMA_MULTIPLE` sigmas, so it sits far out in the tail no matter how
# the tier is sized and the sigma expresses fully. At 3.0 sigma it truncates
# 2 * P(Z > 3) = 0.27% of draws — a genuine outlier guard. The legacy fixed bounds are kept
# as a FLOOR so the calibrated tier's guard cannot become tighter than it was (its sigmas
# are 0.2-0.3, so it never came near either bound and its results are unchanged).
#
# `matrix_retention` keeps its hard physical ceiling of 1.0: a retention FRACTION above
# unity is not a wide tail, it is a different quantity.
_OBSERVABLE_CLAMP: Dict[str, Tuple[float, float]] = {
    "matrix_headspace": (0.01, 100.0),
    "henry_kaw": (0.01, 100.0),
    "matrix_retention": (0.01, 1.0),  # retention factor itself is bounded [0.01, 1]
}

# How many sigmas out the tail guard sits. 3.0 -> 0.27% of draws truncated.
_OBSERVABLE_CLAMP_SIGMA_MULTIPLE = 3.0

# Keys whose clamp is a PHYSICAL bound, not a tail guard, and must not be widened.
_PHYSICALLY_BOUNDED_OBSERVABLES = frozenset({"matrix_retention"})


def _observable_clamp_for(key: str, sigma: float) -> Tuple[float, float]:
    """The sampling clamp for one observable prior. 2026-08-27 (Wave I).

    Returns the wider of (a) the legacy fixed bound and (b) a symmetric multiplicative
    bound `_OBSERVABLE_CLAMP_SIGMA_MULTIPLE` sigmas out, so a wide tier expresses its
    full width instead of being silently truncated at a constant. Physically bounded
    keys are returned unchanged.
    """
    lo, hi = _OBSERVABLE_CLAMP.get(key, (0.01, 100.0))
    if key in _PHYSICALLY_BOUNDED_OBSERVABLES:
        return lo, hi
    sigma = abs(float(sigma))
    if sigma <= 0.0:
        return lo, hi
    span = math.exp(_OBSERVABLE_CLAMP_SIGMA_MULTIPLE * sigma)
    return min(lo, 1.0 / span), max(hi, span)


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


@lru_cache(maxsize=None)
def _routed_offset_keys(prior_key: str) -> Tuple[str, ...]:
    """Every BARRIER_OFFSETS key a prior named `prior_key` must set to actually bite.

    2026-08-27 (Wave I) — MEASURED DEFECT, not a refactor.

    `barrier_constants.get_barrier()` resolves BARRIER_OFFSETS against the NORMALISED RAW
    family label and only canonicalises afterwards, so a prior keyed on a canonical family
    name never matches the label the engine emits. Measured directly on this tree, **10 of
    the 14 keys in `DEFAULT_FAMILY_PRIORS` moved no barrier at all**:

        schiff_condensation, 1,2-enolisation, 2,3-enolisation, cysteine_thermolysis,
        thiol_addition, thiol_addition_hexose, retro_aldol, dehydration,
        beta_elimination, lipid_thiazole

    (`1,2-enolisation` could never have matched anything: its target contains a `-` that
    key normalisation has already replaced with `_`.) Only amadori_rearrangement,
    strecker_degradation, thiol_oxidation and aminoketone_condensation were live.

    So the Monte-Carlo barrier channel was ~70% inert and every published CI was narrower
    than the priors it claimed to sample. **Note the direction: this made the intervals
    too NARROW, so the coverage numbers this repo published were too LOW, not too high.**
    Fixing it widens the intervals and therefore RAISES coverage — a change that flatters
    the model, which is exactly why the before/after must be reported as an interval-width
    change and never as an accuracy improvement. Nothing about any prediction changes; only
    the width of the envelope drawn around it.

    Same root cause as the `src/family_sensitivity.py` false zeros; this reuses that
    module's resolver so the two cannot drift apart.
    """
    from src.family_sensitivity import resolve_offset_keys

    # One prior key is not a family name at all and so canonicalises to itself, leaving it
    # inert even after routing: `lipid_thiazole`. The engine emits
    # `Lipid_Thiazole_Condensation`, which canonicalises to `lipid_condensation`. Aliased
    # here rather than renamed, because the key is published in the artifact's `priors`
    # list and renaming it would break comparison against older runs.
    _PRIOR_KEY_ALIASES = {"lipid_thiazole": "lipid_condensation"}
    resolved_family = _PRIOR_KEY_ALIASES.get(prior_key, prior_key)

    # STILL INERT AFTER THE FIX, and deliberately so: `dehydration`. That canonical family
    # exists in FAST_BARRIERS, but the engine emits `Sugar_Dehydration` and
    # `Thiol_Dehydration`, which canonicalise to themselves — nothing canonicalises to
    # `dehydration`. So this prior samples a family the network does not contain. It is NOT
    # aliased here: which dehydration barrier the MC should perturb (one, the other, or
    # both) is a science decision that changes every published interval, and it belongs to
    # the owner with a stated reason, not to a lookup table in a routing helper. Pinned by
    # tests/unit/test_uncertainty_propagation.py::
    # test_every_default_family_prior_actually_moves_a_barrier, which names it explicitly
    # so it cannot quietly stay broken or quietly become live.
    return tuple(resolve_offset_keys(resolved_family, prior_key))


@contextmanager
def _barrier_offsets(offsets: Mapping[str, float]) -> Iterator[None]:
    previous = os.environ.get("BARRIER_OFFSETS")
    # 2026-08-27 (Wave I): route each prior onto every engine label it names. See
    # `_routed_offset_keys` -- 10 of 14 priors were previously exact no-ops.
    routed: Dict[str, float] = {}
    for key, value in dict(offsets).items():
        for routed_key in _routed_offset_keys(str(key)):
            routed[routed_key] = float(value)
    try:
        os.environ["BARRIER_OFFSETS"] = json.dumps(routed)
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
    matrix_sigma_override: Optional[float] = None,
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
            # 2026-08-27 (Wave I): `matrix_sigma_override` re-scores an identical run at a
            # different matrix_headspace sigma. Used by the external-validation report to
            # show coverage at the PRE-WIDENING sigma (2.0) beside the shipped one (2.86),
            # so "the model got better" and "the interval got wider" stay distinguishable.
            # It affects only the matrix_headspace observable and never the shipped
            # defaults; nothing in the production path passes it.
            if matrix_sigma_override is not None and key == "matrix_headspace":
                sigma = float(matrix_sigma_override)
                source = f"{source}_sigma_override_{sigma:g}"
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
                # 2026-08-27 (Wave I): sigma-derived tail guard. The previous fixed
                # +/-100x clamp truncated 10.7% of draws at the uncalibrated tier's
                # sigma 2.86, so the realised band was narrower than the advertised one.
                lo, hi = _observable_clamp_for(prior.key, prior.sigma_kcal)
                observables[prior.key] = max(lo, min(hi, multiplier))
            else:
                barriers[prior.key] = rng.gauss(0.0, prior.sigma_kcal)
        samples.append({"barrier": barriers, "observable": observables})
    return samples


# Envelopes narrower than this (in dex) are treated as degenerate: they either
# trivially contain a value copied from the model or trivially miss everything.
_MIN_EVALUABLE_CI_WIDTH_LOG10 = 0.01

# 2026-08-27 (Wave I): the implementation moved to `benchmark_validation` so the benchmark
# summary layer can share it without a circular import (that module cannot import this one).
# These names are kept as aliases because several scripts and tests import them by name.
_SYNTHETIC_ORIGINS = SYNTHETIC_BENCHMARK_ORIGINS


def _benchmark_signal_origin(bench_file: Path) -> str:
    """Classify a benchmark's comparator signal as literature-measured or internal/synthetic.

    Thin alias for `src.benchmark_validation.benchmark_signal_origin`, which is now the
    single definition. Kept so existing importers keep working.
    """
    return benchmark_signal_origin(bench_file)


def _benchmark_protein_type(bench_file: Path) -> str:
    """The benchmark's declared protein system ("free", "pea_iso", ...). 2026-08-27."""
    try:
        bench = json.loads(Path(bench_file).read_text())
    except (OSError, json.JSONDecodeError):
        return ""
    return str(bench.get("protein_type", "") or "")


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

    # Audit 2026-08-26 (4.2): the aggregate coverage number silently mixed
    # literature-measured rows with internal synthetic comparators, and counted
    # degenerate zero-width envelopes as ordinary pass/fail rows. Split them.
    origin_by_benchmark: Dict[str, str] = {}
    role_by_benchmark: Dict[str, str] = {}
    for env in benchmark_envelopes:
        origin_by_benchmark[env.benchmark_id] = _benchmark_signal_origin(Path(env.bench_file))
        role_by_benchmark[env.benchmark_id] = benchmark_evidence_role(
            env.benchmark_id, Path(env.bench_file)
        )

    # 2026-08-27 (Wave I). THIRD bucket: `fitted_row`.
    #
    # The 2026-08-26 split separated literature-measured rows from internal synthetic
    # comparators, which was the right first cut but left the worst case hiding inside the
    # literature bucket: a row whose own constant was FITTED to it. Those rows are
    # literature-sourced, so they landed in `external_literature` and were counted as
    # evidence -- and because a fit reproduces its target, they were the rows most likely
    # to be "inside the CI". In Wave H they were literally the only two hits in the
    # headline. A fitted row must be reported, but it must never be in the numerator or
    # the denominator of a coverage claim.
    split: Dict[str, Dict[str, Any]] = {
        "external_literature": {"hits": 0, "total": 0, "not_evaluable": 0, "ci_widths_log10": []},
        "internal_synthetic": {"hits": 0, "total": 0, "not_evaluable": 0, "ci_widths_log10": []},
        "fitted_row": {"hits": 0, "total": 0, "not_evaluable": 0, "ci_widths_log10": []},
    }
    for env in benchmark_envelopes:
        if role_by_benchmark[env.benchmark_id] == "fit_recovery":
            bucket = split["fitted_row"]
        else:
            bucket = split[origin_by_benchmark[env.benchmark_id]]
        for c in env.compounds:
            width = (
                math.log10(c.predicted_p95 / c.predicted_p5)
                if c.predicted_p5 > 0.0 and c.predicted_p95 > 0.0
                else None
            )
            if width is None or width <= _MIN_EVALUABLE_CI_WIDTH_LOG10:
                # A degenerate envelope cannot meaningfully contain (or miss) a
                # measurement — report it, but keep it out of the coverage rate.
                bucket["not_evaluable"] += 1
                continue
            bucket["total"] += 1
            bucket["ci_widths_log10"].append(width)
            if c.inside_ci:
                bucket["hits"] += 1
    for bucket in split.values():
        widths = sorted(bucket.pop("ci_widths_log10"))
        bucket["median_ci_width_log10"] = (
            _percentile(widths, 50.0) if widths else None
        )
        bucket["ci_coverage_rate"] = (
            bucket["hits"] / bucket["total"] if bucket["total"] else None
        )

    lit = split["external_literature"]
    fitted = split["fitted_row"]
    payload: Dict[str, Any] = {
        "summary": {
            "n_samples": n_samples,
            "seed": seed,
            "benchmark_count": len(benchmark_envelopes),
            "matched_compound_count": coverage_total,
            "ci_coverage_hits": coverage_hits,
            "ci_coverage_rate": (coverage_hits / coverage_total) if coverage_total else None,
            "ci_level_pct": 90,
            "signal_origin_split": split,
            # 2026-08-27 (Wave I). The one number to read. `ci_coverage_rate` above is the
            # mixed-population aggregate, kept only for continuity with older artifacts:
            # it pools genuine literature rows with the model reproducing its own frozen
            # snapshots and with rows whose constants were fitted to them.
            "honest_literature_coverage": {
                "hits": lit["hits"],
                "total": lit["total"],
                "rate": lit["ci_coverage_rate"],
                "not_evaluable": lit["not_evaluable"],
                "median_ci_width_log10": lit["median_ci_width_log10"],
                "excluded_fitted_rows": fitted["total"] + fitted["not_evaluable"],
                "excluded_fitted_rows_that_would_have_been_hits": fitted["hits"],
                "definition": (
                    "External-literature rows only, with fitted rows (constants back-solved "
                    "from the benchmark) and internal synthetic reproducibility rows removed "
                    "from BOTH numerator and denominator. Read it together with "
                    "median_ci_width_log10: a wide interval covering a measurement is a weak "
                    "claim."
                ),
            },
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
                "signal_origin": origin_by_benchmark[env.benchmark_id],
                # 2026-08-27 (Wave I): "predictive" | "fit_recovery" | "internal_synthetic".
                "evidence_role": role_by_benchmark[env.benchmark_id],
                # The disclosure flag the Wave I fit-target gate requires: True means the
                # constants under test were derived from THIS benchmark with enough freedom
                # to reproduce it row by row, so its rows are excluded from the
                # literature-coverage numerator and denominator.
                "fitted_row": role_by_benchmark[env.benchmark_id] == "fit_recovery",
                # Every live fit record naming this benchmark, INCLUDING low-leverage
                # global fits that do not trigger exclusion. A reader can then see that no
                # literature row here is strictly out-of-sample even when it still counts.
                "fit_target_of": list(fit_target_records_for(env.benchmark_id)),
                "execution_path": env.execution_path,
                # 2026-08-27 (Wave I): consumed by src.experiment_value._suggest_template,
                # which was prescribing protein-matrix DoE protocols for free-precursor
                # benchmarks because it could not see the system.
                "protein_type": _benchmark_protein_type(Path(env.bench_file)),
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


def _render_signal_origin_split(split: Mapping[str, Any]) -> str:
    """Render the honest coverage split (audit 2026-08-26).

    The aggregate headline mixes literature-measured rows with internal
    synthetic comparators (whose reference values are frozen model output) and
    with degenerate zero-width envelopes. Only the literature slice is
    validation evidence.
    """
    if not split:
        return ""
    lines = [
        "> **Coverage split — this is the headline; the aggregate is not.** The aggregate",
        "> pools three populations that support completely different claims:",
        ">",
        "> * **External literature** — a published measurement the model was not fitted to.",
        ">   Only this row is validation evidence.",
        "> * **Fitted rows** — the constants under test were back-solved from these",
        ">   benchmarks with enough freedom to reproduce them row by row, so agreement here",
        ">   is algebraic recovery. They are excluded from the literature numerator AND",
        ">   denominator and reported separately. Read their outcomes: a row the model",
        ">   still *fails* after being fitted to it is a strong negative result.",
        "> * **Internal synthetic** — the comparator is the model's own frozen output.",
        ">   Agreement means the model reproduces itself. Not evidence about chemistry.",
        ">",
        "> Zero-width (degenerate) envelopes are excluded from coverage and counted as",
        "> not-evaluable — a predicted==measured synthetic hit trivially 'contains' its own",
        "> value. Coverage is only interpretable next to its median CI width: a wide",
        "> interval makes coverage cheap.",
        ">",
        "> | Signal origin | Inside 90% CI | Not evaluable | Median CI width (dex) |",
        "> | --- | ---: | ---: | ---: |",
    ]
    labels = {
        "external_literature": "External literature (the only validation evidence)",
        "fitted_row": "Fitted rows (fit recovery — NOT evidence)",
        "internal_synthetic": "Internal synthetic (reproducibility only — NOT evidence)",
    }
    for key, label in labels.items():
        bucket = split.get(key) or {}
        hits = bucket.get("hits", 0)
        total = bucket.get("total", 0)
        rate = bucket.get("ci_coverage_rate")
        rate_str = f" ({rate * 100:.0f}%)" if rate is not None else ""
        width = bucket.get("median_ci_width_log10")
        width_str = f"{width:.2f}" if width is not None else "n/a"
        lines.append(
            f"> | {label} | {hits}/{total}{rate_str} | {bucket.get('not_evaluable', 0)} | {width_str} |"
        )
    return "\n".join(lines)


def _render_honest_headline(summary: Mapping[str, Any]) -> str:
    """The one number a reader should take away. 2026-08-27 (Wave I)."""
    honest = summary.get("honest_literature_coverage") or {}
    if not honest:
        return (
            "**Headline trust metric**: unavailable — this artifact predates the "
            "Wave I coverage split. Regenerate it."
        )
    hits = honest.get("hits", 0)
    total = honest.get("total", 0)
    rate = honest.get("rate")
    rate_str = f" (**{rate * 100:.1f}%**)" if rate is not None else ""
    width = honest.get("median_ci_width_log10")
    width_str = (
        f" Median CI width **{width:.2f} dex** (~{10 ** width:.0f}× end to end) — read the "
        "coverage with the width."
        if width is not None
        else ""
    )
    excluded = honest.get("excluded_fitted_rows", 0)
    would_have_hit = honest.get("excluded_fitted_rows_that_would_have_been_hits", 0)
    excluded_str = (
        f" **{excluded}** fitted row(s) are excluded from both numerator and denominator "
        f"({would_have_hit} of them would have counted as hits); see the split below."
        if excluded
        else ""
    )
    not_evaluable = honest.get("not_evaluable", 0)
    ne_str = f" {not_evaluable} literature row(s) are not evaluable (degenerate envelope)." if not_evaluable else ""
    return (
        f"**Headline trust metric — external literature only**: the measured value lies "
        f"inside the 90% CI for **{hits} / {total}** literature rows{rate_str}."
        f"{excluded_str}{ne_str}{width_str}"
    )


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
        # 2026-08-27 (Wave I): the aggregate used to be the headline. It pooled genuine
        # literature rows with the model reproducing its own frozen snapshots and with
        # rows whose constants had been fitted to them -- and in Wave H the only two
        # "hits" in the literature slice were fitted rows. The honest number leads now;
        # the aggregate is demoted to a labelled secondary line.
        _render_honest_headline(summary),
        "",
        _render_signal_origin_split(summary.get("signal_origin_split") or {}),
        "",
        f"_Secondary, mixed-population figure, retained only for continuity with older "
        f"reports — do not quote it: measured value lies inside 90% CI for "
        f"{summary.get('ci_coverage_hits', 0)} / {summary.get('matched_compound_count', 0)} "
        f"matched compounds ({coverage_str}), pooling literature, fitted and synthetic rows._",
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
