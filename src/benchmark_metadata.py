"""Engine-neutral facts about a benchmark bundle: metadata, signal origin, evidence role,
and the scale contract it declares.

2026-09-03 (retirement step B3). These helpers used to live in ``src/benchmark_validation.py``,
the LEGACY validation harness, and were imported from there by the kinetic core's envelope
(``src/kinetic_core/uncertainty.py``). The harness is scheduled for deletion (retirement
step B5); nothing here runs an engine, so it moves to a module that will outlive it.
``benchmark_validation`` re-exports every name so its 36 test files keep working.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping

from src.benchmark_types import BenchmarkMetadata
from src.fit_target_index import is_per_row_fit_target
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

def benchmark_evidence_role(benchmark_id: str, bench_file: Path) -> str:
    """What KIND of claim this benchmark can support: see BenchmarkSummary.evidence_role.

    2026-08-27 (Wave I). The mechanical status says whether predictions and measurements
    agree; this says whether that agreement is evidence. Fit-recovery is checked FIRST and
    beats signal origin: a lane whose observable factor was solved from a literature
    benchmark is still literature-sourced, and is still not a prediction.

    Two independent sources of fit-recovery are consulted:

    * the calibration registry's own `fitted_to_benchmark` declarations (the matrix
      observability factors, one free factor per compound per lane); and
    * `src.fit_target_index`, which reads the fit records under `results/validation/` and
      classifies each by LEVERAGE. Only high-leverage fits (enough free parameters to
      reproduce their target row by row) make a benchmark fit-recovery. A low-leverage
      global fit -- two constants across two dozen rows -- does NOT, because excluding
      those rows would delete genuine failures from the count rather than expose them.
    """
    # 2026-09-03 (B5): the legacy matrix observability factors' `fitted_to_benchmark`
    # declarations are no longer consulted; the fit-target index is the one source.
    if is_per_row_fit_target(benchmark_id):
        return "fit_recovery"
    if benchmark_signal_origin(Path(bench_file)) == "internal_synthetic":
        return "internal_synthetic"
    return "predictive"


#: Backwards-compatible alias (the legacy harness called it by this name).
_resolve_scale_thresholds = resolve_scale_thresholds

__all__ = [
    "BenchmarkMetadata",
    "SYNTHETIC_BENCHMARK_ORIGINS",
    "benchmark_evidence_role",
    "benchmark_signal_origin",
    "get_benchmark_metadata",
    "matrix_source_anchor",
    "resolve_scale_thresholds",
]
