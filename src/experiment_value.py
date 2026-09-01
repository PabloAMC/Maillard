"""S22.1 — Experiment-value (VoI) ranker.

Consumes the Monte-Carlo envelope from `src.uncertainty_propagation`
(`results/validation/prediction_uncertainty.json`) and ranks `(benchmark,
compound)` pairs by the value of running a confirmatory experiment now.

Score components (deliberately simple — no ML, easily auditable):

* `envelope_miss`   — 0 if the measured value is inside the 90% CI, otherwise
  ``log10`` distance from the nearest CI bound.  Captures "the model is
  demonstrably wrong about this compound today".
* `ci_width_log10`  — width of the 90% CI in dex.  Captures "the model is
  also uncertain here".
* `decision_relevance` — ``log10(max(measured, P50) / ODT)`` clipped to
  [0.5, 5.0] when an odour threshold is known.  Compounds that exceed their
  threshold by orders of magnitude matter more than near-noise tails.

We do not try to invent monetary cost. The output is an ordered backlog;
`scripts/request_experiment.py` (S22.2) will turn the top-N into concrete
lab protocols.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml

from src import data_paths
from src import data_access

PREDICTION_UNCERTAINTY_PATH = data_paths.PREDICTION_UNCERTAINTY
DESIRABLE_TARGETS_PATH = data_paths.DESIRABLE_TARGETS
OFF_FLAVOUR_TARGETS_PATH = data_paths.OFF_FLAVOUR_TARGETS

# Compound classes for DoE template suggestion.
_MEATY_KEYWORDS = ("furanthiol", "furfurylthiol", "methional", "thiazole", "mft", "fft")
_OFFNOTE_KEYWORDS = ("hexanal", "nonanal", "octenal", "pentylfuran", "hexanol")
_SAFETY_KEYWORDS = ("acrylamide", "cml", "cel", "furosine", "hmf")

# Canonical short matrix families (must align with `src.matrix_experiment_intake._PROTEIN_TYPES`
# and the matrix-family registry under data/lit/matrix_family_coverage_registry.json).
# We intentionally use the short protein-type names because that is what the formulation
# pipeline already emits via `FormulationResult.confidence_metadata['scope_assessment']`.
# End-to-end filtering then needs no translation table.
MATRIX_FAMILIES: Tuple[str, ...] = ("free", "pea_iso", "soy_iso", "wheat_gluten", "myco", "unknown")

# Substring -> canonical matrix family. Order matters: more specific tokens first.
_BENCHMARK_ID_MATRIX_RULES: Tuple[Tuple[str, str], ...] = (
    ("pea_isolate", "pea_iso"),
    ("pea_iso", "pea_iso"),
    ("soy_isolate", "soy_iso"),
    ("soy_iso", "soy_iso"),
    ("spi_", "soy_iso"),
    ("wheat_gluten", "wheat_gluten"),
    ("mycoprotein", "myco"),
    ("myco_", "myco"),
    # free-precursor / model-aqueous benchmarks (match by typical precursor naming)
    ("cys_", "free"),
    ("thiamine_", "free"),
    ("acrylamide_asparagine", "free"),
    ("furosine_extrusion", "free"),
    ("cml_cel_commercial", "free"),
)


def infer_matrix_family(
    benchmark_id: str,
    benchmark: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the canonical short matrix family for a benchmark.

    Resolution order:
    1. Explicit ``protein_type`` field on the benchmark JSON when present.
    2. Substring match on ``benchmark_id`` against ``_BENCHMARK_ID_MATRIX_RULES``.
    3. Fallback ``"unknown"``.

    Returned values use the same short-name convention as
    ``src.matrix_experiment_intake._PROTEIN_TYPES`` so a formulation's
    ``protein_type`` can be passed through unchanged as a filter.
    """
    if benchmark is not None:
        explicit = benchmark.get("protein_type")
        if explicit:
            return str(explicit).strip() or "unknown"
    bid = (benchmark_id or "").lower()
    for token, family in _BENCHMARK_ID_MATRIX_RULES:
        if token in bid:
            return family
    return "unknown"


@dataclass(frozen=True)
class CompoundSpec:
    canonical_name: str
    odour_threshold_ug_per_kg: Optional[float]
    priority: Optional[str]


@dataclass(frozen=True)
class ExperimentCandidate:
    rank: int
    benchmark_id: str
    compound: str
    measured_ppb: float
    predicted_p5: float
    predicted_p50: float
    predicted_p95: float
    inside_ci: bool
    envelope_miss_log10: float
    ci_width_log10: Optional[float]
    odour_threshold_ug_per_kg: Optional[float]
    decision_relevance: float
    voi_score: float
    suggested_doe_template: str
    rationale: str
    # S23 — matrix-family attribution. Defaulted so legacy callers/fixtures
    # constructing ExperimentCandidate directly without the field keep working.
    matrix_family: str = "unknown"


# ---------------------------------------------------------------------------
# Compound lookup
# ---------------------------------------------------------------------------

def _normalise(name: str) -> str:
    """Lowercase + strip parenthetical aliases + collapse non-alnum to single
    spaces. Examples:
        '2-Methyl-3-furanthiol (MFT)' -> '2 methyl 3 furanthiol'
        '2-furfurylthiol'             -> '2 furfurylthiol'
    """
    if not name:
        return ""
    # remove parenthetical aliases
    stripped = re.sub(r"\([^)]*\)", " ", name)
    stripped = stripped.lower()
    stripped = re.sub(r"[^a-z0-9]+", " ", stripped).strip()
    return stripped


def _alias_keys(canonical: str) -> List[str]:
    """All match keys for a compound: full canonical + each parenthetical token."""
    keys = [_normalise(canonical)]
    for token in re.findall(r"\(([^)]*)\)", canonical):
        token_norm = _normalise(token)
        if token_norm and token_norm not in keys:
            keys.append(token_norm)
    return keys


def load_compound_specs(
    paths: Sequence[Path] = (DESIRABLE_TARGETS_PATH, OFF_FLAVOUR_TARGETS_PATH),
) -> Dict[str, CompoundSpec]:
    """Return {normalised_alias: CompoundSpec}.

    Multiple aliases (full name + parenthetical) point at the same spec.
    """
    specs: Dict[str, CompoundSpec] = {}
    for path in paths:
        # Strict since 2026-09-01: a missing or malformed species file used to be skipped.
        data = data_access.load_yaml(path) or {}
        for entry in data.get("compounds", []) or []:
            name = str(entry.get("name") or "").strip()
            if not name:
                continue
            spec = CompoundSpec(
                canonical_name=name,
                odour_threshold_ug_per_kg=(
                    float(entry["odour_threshold_ug_per_kg"])
                    if entry.get("odour_threshold_ug_per_kg") not in (None, "")
                    else None
                ),
                priority=entry.get("priority"),
            )
            for alias in _alias_keys(name):
                # First spec wins for each alias.
                specs.setdefault(alias, spec)
    return specs


def lookup_spec(
    compound_name: str, specs: Mapping[str, CompoundSpec]
) -> Optional[CompoundSpec]:
    key = _normalise(compound_name)
    if not key:
        return None
    if key in specs:
        return specs[key]
    # token containment fallback for noisy names like "2-methyl-3-furanthiol-derived adduct"
    for alias, spec in specs.items():
        if alias and alias in key:
            return spec
    return None


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def _envelope_miss_log10(measured: float, p5: float, p95: float) -> float:
    if measured <= 0.0 or not math.isfinite(measured):
        return 0.0
    if p5 <= measured <= p95:
        return 0.0
    if measured > p95 and p95 > 0.0:
        return float(math.log10(measured / p95))
    if measured < p5 and measured > 0.0:
        return float(math.log10(p5 / measured))
    return 0.0


def _decision_relevance(
    measured: float, p50: float, odt: Optional[float]
) -> float:
    if odt is None or odt <= 0.0 or not math.isfinite(odt):
        return 1.0
    reference = max(measured, p50, 0.0)
    if reference <= 0.0:
        return 0.5
    ratio = reference / odt
    if ratio <= 0.0:
        return 0.5
    log_ratio = math.log10(ratio)
    return float(max(0.5, min(5.0, log_ratio)))


def _suggest_template(
    compound: str,
    ci_width_log10: Optional[float],
    inside_ci: bool,
    *,
    execution_path: str = "",
    protein_type: str = "",
) -> Tuple[str, str]:
    """Return (template_key, rationale_fragment).

    2026-08-27 (Wave I) — the template now depends on the BENCHMARK'S SYSTEM, not only on
    a substring of the compound name.

    The red team found the rank-1 experiment card prescribing a protein-matrix protocol
    ("Matrix (SPI, PPI)", "Standard PBMA formulation baseline") for
    `thiamine_cys_glucose_120C_Bolton1994` -- a FREE-PRECURSOR aqueous thiamine + cysteine
    + glucose system with no protein in it at all. The cause: selection matched "furanthiol"
    against `_MEATY_KEYWORDS` and returned `blocking_benchmark_gap` unconditionally, and
    that template's whole point is to vary the protein matrix. A CRO following the card
    would have run the wrong experiment.
    """
    name = compound.lower()
    free_precursor_system = (
        str(protein_type).strip().lower() in {"free", "free_amino_acid", ""}
        and str(execution_path).strip().lower() in {"free_precursor", ""}
    ) or str(execution_path).strip().lower() == "free_precursor"

    if any(k in name for k in _SAFETY_KEYWORDS):
        return (
            "missing_absolute_anchor",
            "safety marker — needs SIDA-grade absolute anchor",
        )
    if any(k in name for k in _MEATY_KEYWORDS):
        if free_precursor_system:
            # Varying "Matrix (SPI, PPI)" is meaningless where there is no protein matrix.
            return (
                "free_precursor_sulfur_yield",
                "critical meaty odorant in a FREE-PRECURSOR system — the open question is "
                "absolute yield vs precursor dose and temperature, not matrix transfer",
            )
        return (
            "blocking_benchmark_gap",
            "critical meaty odorant — multi-factor SIDA closes precursor × matrix gap",
        )
    if any(k in name for k in _OFFNOTE_KEYWORDS):
        return (
            "missing_positive_flavor_anchor",
            "off-note marker — focused band quantitation",
        )
    if ci_width_log10 is not None and ci_width_log10 >= 1.5:
        return (
            "missing_kinetic_dataset",
            "wide envelope — time-course narrows the rate-limiting step",
        )
    return (
        "missing_absolute_anchor",
        "anchor measurement to constrain prior",
    )


def _voi_score(
    *,
    inside_ci: bool,
    envelope_miss_log10: float,
    ci_width_log10: Optional[float],
    decision_relevance: float,
) -> float:
    miss_term = 1.0 + envelope_miss_log10  # >=1 outside, 1 inside (no penalty)
    if inside_ci:
        miss_term = 0.0  # only consider width term when inside CI
    width = ci_width_log10 if ci_width_log10 is not None else 0.0
    width_term = 0.3 * width
    return float((miss_term + width_term) * decision_relevance)


def _build_rationale(
    *,
    inside_ci: bool,
    envelope_miss_log10: float,
    ci_width_log10: Optional[float],
    decision_relevance: float,
    measured: float,
    p50: float,
    odt: Optional[float],
    template_rationale: str,
) -> str:
    parts: List[str] = []
    if not inside_ci:
        parts.append(
            f"measured outside 90% CI by {envelope_miss_log10:.2f} dex"
        )
    if ci_width_log10 is not None:
        parts.append(f"CI width {ci_width_log10:.2f} dex")
    if odt is not None and odt > 0.0:
        ratio_ref = max(measured, p50, 0.0)
        if ratio_ref > 0.0:
            parts.append(
                f"≈{ratio_ref / odt:.1g}× ODT (decision_relevance={decision_relevance:.2f})"
            )
    parts.append(template_rationale)
    return "; ".join(parts)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def rank_experiments(
    payload: Mapping[str, Any],
    *,
    specs: Optional[Mapping[str, CompoundSpec]] = None,
    top_n: Optional[int] = None,
) -> List[ExperimentCandidate]:
    """Rank `(benchmark, compound)` pairs from a prediction-uncertainty payload."""
    if specs is None:
        specs = load_compound_specs()

    rows: List[Tuple[float, ExperimentCandidate]] = []
    for benchmark in payload.get("benchmarks", []) or []:
        bench_id = str(benchmark.get("benchmark_id") or "unknown")
        matrix_family = infer_matrix_family(bench_id, benchmark)
        for compound in benchmark.get("compounds", []) or []:
            measured = float(compound.get("measured_ppb", 0.0) or 0.0)
            p5 = float(compound.get("predicted_p5", 0.0) or 0.0)
            p50 = float(compound.get("predicted_p50", 0.0) or 0.0)
            p95 = float(compound.get("predicted_p95", 0.0) or 0.0)
            inside_ci = bool(compound.get("inside_ci"))
            ci_width_log10 = compound.get("ci_width_log10")
            if ci_width_log10 is not None:
                try:
                    ci_width_log10 = float(ci_width_log10)
                except (TypeError, ValueError):
                    ci_width_log10 = None
            spec = lookup_spec(str(compound.get("compound", "")), specs)
            odt = spec.odour_threshold_ug_per_kg if spec else None
            envelope_miss = _envelope_miss_log10(measured, p5, p95)
            decision_relevance = _decision_relevance(measured, p50, odt)
            template_key, template_rationale = _suggest_template(
                str(compound.get("compound", "")),
                ci_width_log10,
                inside_ci,
                # 2026-08-27 (Wave I): the benchmark's own system decides whether a
                # matrix-varying protocol makes any sense. See _suggest_template.
                execution_path=str(benchmark.get("execution_path", "") or ""),
                protein_type=str(benchmark.get("protein_type", "") or ""),
            )
            score = _voi_score(
                inside_ci=inside_ci,
                envelope_miss_log10=envelope_miss,
                ci_width_log10=ci_width_log10,
                decision_relevance=decision_relevance,
            )
            rationale = _build_rationale(
                inside_ci=inside_ci,
                envelope_miss_log10=envelope_miss,
                ci_width_log10=ci_width_log10,
                decision_relevance=decision_relevance,
                measured=measured,
                p50=p50,
                odt=odt,
                template_rationale=template_rationale,
            )
            candidate = ExperimentCandidate(
                rank=0,  # filled after sort
                benchmark_id=bench_id,
                matrix_family=matrix_family,
                compound=str(compound.get("compound", "")),
                measured_ppb=measured,
                predicted_p5=p5,
                predicted_p50=p50,
                predicted_p95=p95,
                inside_ci=inside_ci,
                envelope_miss_log10=envelope_miss,
                ci_width_log10=ci_width_log10,
                odour_threshold_ug_per_kg=odt,
                decision_relevance=decision_relevance,
                voi_score=score,
                suggested_doe_template=template_key,
                rationale=rationale,
            )
            rows.append((score, candidate))

    rows.sort(key=lambda pair: pair[0], reverse=True)
    ranked: List[ExperimentCandidate] = []
    for idx, (_score, candidate) in enumerate(rows, start=1):
        ranked.append(
            ExperimentCandidate(
                rank=idx,
                benchmark_id=candidate.benchmark_id,
                matrix_family=candidate.matrix_family,
                compound=candidate.compound,
                measured_ppb=candidate.measured_ppb,
                predicted_p5=candidate.predicted_p5,
                predicted_p50=candidate.predicted_p50,
                predicted_p95=candidate.predicted_p95,
                inside_ci=candidate.inside_ci,
                envelope_miss_log10=candidate.envelope_miss_log10,
                ci_width_log10=candidate.ci_width_log10,
                odour_threshold_ug_per_kg=candidate.odour_threshold_ug_per_kg,
                decision_relevance=candidate.decision_relevance,
                voi_score=candidate.voi_score,
                suggested_doe_template=candidate.suggested_doe_template,
                rationale=candidate.rationale,
            )
        )
    if top_n is not None:
        ranked = ranked[: max(top_n, 0)]
    return ranked


def filter_by_matrix(
    candidates: Sequence[ExperimentCandidate],
    matrix_filter: Optional[Sequence[str]],
) -> List[ExperimentCandidate]:
    """Restrict candidates to the requested matrix families and re-rank.

    ``matrix_filter`` accepts the canonical short names (``pea_iso``,
    ``soy_iso``, ``wheat_gluten``, ``myco``, ``free``). Empty / ``None``
    returns the input unchanged.
    """
    if not matrix_filter:
        return list(candidates)
    wanted = {str(m).strip().lower() for m in matrix_filter if str(m).strip()}
    if not wanted:
        return list(candidates)
    kept = [c for c in candidates if c.matrix_family.lower() in wanted]
    re_ranked: List[ExperimentCandidate] = []
    for idx, c in enumerate(kept, start=1):
        re_ranked.append(
            ExperimentCandidate(
                rank=idx,
                benchmark_id=c.benchmark_id,
                matrix_family=c.matrix_family,
                compound=c.compound,
                measured_ppb=c.measured_ppb,
                predicted_p5=c.predicted_p5,
                predicted_p50=c.predicted_p50,
                predicted_p95=c.predicted_p95,
                inside_ci=c.inside_ci,
                envelope_miss_log10=c.envelope_miss_log10,
                ci_width_log10=c.ci_width_log10,
                odour_threshold_ug_per_kg=c.odour_threshold_ug_per_kg,
                decision_relevance=c.decision_relevance,
                voi_score=c.voi_score,
                suggested_doe_template=c.suggested_doe_template,
                rationale=c.rationale,
            )
        )
    return re_ranked


def build_ranking_payload(
    candidates: Sequence[ExperimentCandidate],
    *,
    source_path: Optional[Path] = None,
    matrix_filter: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Convert ranked candidates into a JSON-serialisable dict."""
    return {
        "source": str(source_path) if source_path else None,
        "matrix_filter": list(matrix_filter) if matrix_filter else None,
        "candidate_count": len(candidates),
        "miss_count": sum(1 for c in candidates if not c.inside_ci),
        "candidates": [
            {
                "rank": c.rank,
                "benchmark_id": c.benchmark_id,
                "matrix_family": c.matrix_family,
                "compound": c.compound,
                "voi_score": c.voi_score,
                "inside_ci": c.inside_ci,
                "envelope_miss_log10": c.envelope_miss_log10,
                "ci_width_log10": c.ci_width_log10,
                "decision_relevance": c.decision_relevance,
                "odour_threshold_ug_per_kg": c.odour_threshold_ug_per_kg,
                "measured_ppb": c.measured_ppb,
                "predicted_p5": c.predicted_p5,
                "predicted_p50": c.predicted_p50,
                "predicted_p95": c.predicted_p95,
                "suggested_doe_template": c.suggested_doe_template,
                "rationale": c.rationale,
            }
            for c in candidates
        ],
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    candidates = payload.get("candidates", []) or []
    matrix_filter = payload.get("matrix_filter")
    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_experiment_value_ranking.py". Manual edits will be overwritten. -->',
        "",
        "# Experiment Value Ranking (VoI)",
        "",
        "_Ranks `(benchmark, compound)` pairs by the value of a confirmatory experiment now: combines envelope miss, CI width, and ODT-anchored decision relevance._",
        "",
        f"Total candidates: **{payload.get('candidate_count', 0)}** "
        f"(out-of-CI: **{payload.get('miss_count', 0)}**). "
        f"Source: `{payload.get('source') or 'in-memory'}`.",
        "",
    ]
    if matrix_filter:
        lines.append(
            f"Matrix filter: **{', '.join(matrix_filter)}** "
            "(pass `--matrix '<family>'` to `experiment-value-ranking` / `next-experiment` to change)."
        )
        lines.append("")
    lines.extend([
        "| Rank | VoI | Benchmark | Matrix | Compound | In CI | Miss (dex) | Width (dex) | Meas (ppb) | P50 (ppb) | DoE template | Rationale |",
        "| ---: | ---: | --- | --- | --- | :---: | ---: | ---: | ---: | ---: | --- | --- |",
    ])
    for c in candidates:
        in_ci = "✓" if c.get("inside_ci") else "✗"
        width = c.get("ci_width_log10")
        width_str = f"{float(width):.2f}" if width is not None else "n/a"
        lines.append(
            f"| {c.get('rank')} | {float(c.get('voi_score', 0.0)):.2f} | "
            f"`{c.get('benchmark_id')}` | `{c.get('matrix_family', 'unknown')}` | "
            f"{c.get('compound')} | {in_ci} | "
            f"{float(c.get('envelope_miss_log10', 0.0)):.2f} | {width_str} | "
            f"{float(c.get('measured_ppb', 0.0)):.3g} | {float(c.get('predicted_p50', 0.0)):.3g} | "
            f"`{c.get('suggested_doe_template')}` | {c.get('rationale')} |"
        )
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_artifact(
    payload: Mapping[str, Any],
    *,
    output_dir: Path | str = data_paths.VALIDATION_DIR,
    basename: str = "experiment_value_ranking",
) -> Dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{basename}.json"
    md_path = output_dir / f"{basename}.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    return {"json": json_path, "md": md_path}


def load_prediction_payload(
    path: Path | str = PREDICTION_UNCERTAINTY_PATH,
) -> Dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
