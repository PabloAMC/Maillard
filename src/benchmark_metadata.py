"""Engine-neutral facts about a benchmark bundle: metadata, signal origin, evidence role,
and the scale contract it declares.

2026-09-03 (retirement step B3). These helpers used to live in the legacy validation harness
(``src/benchmark_validation.py``, deleted at B5b) and are what the kinetic core's scorecard and
envelope read. Nothing here runs an engine. The evidence-role vocabulary (predictive /
fit_recovery / internal_synthetic) has ONE implementation, ``kinetic_core.fit_targets.core_evidence_role``;
this module supplies the signal origin it builds on.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping

from src.benchmark_types import BenchmarkMetadata
from src.validation_contract import BenchmarkThresholds, DEFAULT_VALIDATION_CONTRACT


def _infer_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    benchmark_id = bench.get("benchmark_id", "unknown")
    protein_type = bench.get("protein_type", "free")
    if protein_type != "free":
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="matrix_headspace",
            execution_path="matrix_only",
            benchmark_engine="matrix_intake_headspace",
            comparator_signal="predicted_ppb",
            cantera_role="not_authoritative",
            target_snapshot_policy="excluded",
            thermodynamic_gating_policy="not_applicable",
            notes="Matrix-only benchmark requiring a dedicated precursor-accessibility path.",
        )
    if "cys_" in benchmark_id:
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="free_aa_sulfur",
            execution_path="free_precursor",
            benchmark_engine="fast_observable",
            comparator_signal="predicted_ppb",
            cantera_role="diagnostic_reference_only",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="diagnostic_only",
            notes="Free amino-acid sulfur benchmark.",
        )
    if "acrylamide" in benchmark_id or bench.get("benchmark_type") == "safety":
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="safety",
            execution_path="free_precursor",
            benchmark_engine="fast_observable",
            comparator_signal="predicted_ppb",
            cantera_role="diagnostic_reference_only",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="diagnostic_only",
            notes="Safety-critical benchmark (e.g., Acrylamide).",
        )
    return BenchmarkMetadata(
        tier="SECONDARY",
        family="general",
        execution_path="free_precursor",
        benchmark_engine="fast_observable",
        comparator_signal="predicted_ppb",
        cantera_role="diagnostic_reference_only",
        target_snapshot_policy="included",
        thermodynamic_gating_policy="diagnostic_only",
        notes=None,
    )

def get_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    metadata = bench.get("metadata") or {}
    inferred = _infer_benchmark_metadata(bench)
    policy = DEFAULT_VALIDATION_CONTRACT.policy_for_execution_path(str(metadata.get("execution_path", inferred.execution_path)))
    return BenchmarkMetadata(
        tier=str(metadata.get("tier", inferred.tier)),
        family=str(metadata.get("family", inferred.family)),
        execution_path=policy.execution_path,
        benchmark_engine=str(metadata.get("benchmark_engine", policy.benchmark_engine)),
        comparator_signal=str(metadata.get("comparator_signal", policy.comparator_signal)),
        cantera_role=str(metadata.get("cantera_role", policy.cantera_role)),
        target_snapshot_policy=str(metadata.get("target_snapshot_policy", policy.target_snapshot_policy)),
        thermodynamic_gating_policy=str(metadata.get("thermodynamic_gating_policy", policy.thermodynamic_gating_policy)),
        notes=metadata.get("notes", inferred.notes),
    )

def resolve_scale_thresholds(
    bench: dict,
    *,
    protein_type: str,
    thresholds: BenchmarkThresholds,
) -> Dict[str, float]:
    configured = (bench.get("validation_contract") or {}).get("scale_thresholds") or {}
    return {
        "max_ratio": float(configured.get("max_ratio", thresholds.ratio_threshold_for(protein_type))),
        "mean_abs_log10_error": float(
            configured.get(
                "mean_abs_log10_error",
                thresholds.mean_abs_log10_error_threshold_for(protein_type),
            )
        ),
    }

def matrix_source_anchor(bench: Mapping[str, Any]) -> str:
    """Return the benchmark's citable external anchor, or "" if it has none.

    2026-08-27 (audit remediation, DOI-less identifier retyping): four sources in
    this repo are genuinely DOI-less (two theses, a US patent, a journal with no
    DOI registration) and their identifiers used to be stored in fields named
    ``doi`` / ``source_doi``, which is both dishonest and a citation-gate
    violation. They now carry a typed ``identifier`` + ``identifier_scheme``
    pair. A typed non-DOI identifier is exactly as much an external anchor as a
    DOI, so every consumer that tested ``source_doi`` truthiness reads this
    helper instead — otherwise retyping an identifier would silently *downgrade*
    a real external source to "unspecified origin".
    """
    doi = str(bench.get("source_doi") or "").strip()
    if doi:
        return doi
    identifier = str(bench.get("identifier") or "").strip()
    if identifier:
        return identifier
    return ""

# 2026-08-27 (Wave I). Canonical home of the signal-origin classifier. It used to live
# only in `uncertainty_propagation._benchmark_signal_origin`, but the benchmark summary
# layer needs it too and `uncertainty_propagation` imports THIS module, so keeping the
# implementation there would have meant a cycle or a second, drifting copy.
SYNTHETIC_BENCHMARK_ORIGINS = frozenset(
    {
        "synthetic_diagnostic",
        "internal_reproducibility_candidate",
        "internal_experiment",
    }
)

def benchmark_signal_origin(bench_file: Path) -> str:
    """Classify a benchmark's comparator signal as literature-measured or internal/synthetic.

    Returns ``"external_literature"`` or ``"internal_synthetic"``.
    """
    try:
        bench = json.loads(Path(bench_file).read_text())
    except (OSError, json.JSONDecodeError):
        return "internal_synthetic"
    origin = str((bench.get("source_metadata") or {}).get("origin", "")).strip().lower()
    if origin in SYNTHETIC_BENCHMARK_ORIGINS:
        return "internal_synthetic"
    # `matrix_source_anchor` also accepts the typed `identifier`/`identifier_scheme` pair
    # that DOI-less sources (theses, patents, DOI-free journals) now carry, so retyping an
    # identifier out of a `doi`-named field cannot reclassify a literature row as
    # internal/synthetic.
    if matrix_source_anchor(bench) or origin.startswith("external"):
        return "external_literature"
    return "internal_synthetic"



__all__ = [
    "BenchmarkMetadata",
    "SYNTHETIC_BENCHMARK_ORIGINS",
    "benchmark_signal_origin",
    "get_benchmark_metadata",
    "matrix_source_anchor",
    "resolve_scale_thresholds",
]
