"""
src/benchmark_validation.py — FACADE MODULE (Priority 2 refactor)

This file has been decomposed into four focused sub-modules:

  - src/benchmark_registry.py   → IO, metadata, conditions/formulation parsing
  - src/benchmark_evaluator.py  → Prediction execution and statistical helpers
  - src/benchmark_reporting.py  → BenchmarkSummary, artifacts, markdown, snapshots
  - src/benchmark_assertions.py → Matrix assertion rows

This facade re-exports every public symbol from those modules so that all
downstream imports of the form ``from src.benchmark_validation import X``
continue to work unchanged.  New code should import directly from the
appropriate sub-module.
"""
from __future__ import annotations

# ── Types (dataclasses shared across all sub-modules) ─────────────────────────
from src.benchmark_types import (
    BenchmarkMetadata,
    BenchmarkIndexEntry,
    BenchmarkNotSupportedError,
    CompoundComparison,
    BenchmarkEvaluation,
    BenchmarkSummary,
    MatrixBenchmarkDelta,
    MatrixBenchmarkEvidence,
    MatrixBenchmarkBranchDelta,
    ThermodynamicGatingAudit,
    BenchmarkTargetSnapshot,
    MatrixBenchmarkAssertion,
    MatrixPromotionFamilyStatus,
)
from src.validation_contract import BenchmarkThresholds, DEFAULT_VALIDATION_CONTRACT  # noqa: F401

# ── Registry (IO + metadata) ─────────────────────────────────────────────────
from src.benchmark_registry import (
    ROOT,
    BENCHMARK_DIR,
    DEFAULT_TARGET_TAG,
    MATRIX_NAMES,
    THERMODYNAMIC_GATING_POLICIES,
    get_benchmark_files,
    load_benchmark,
    get_benchmark_metadata,
    resolve_thermodynamic_gating_mode,
    is_matrix_only_benchmark,
    get_matrix_ranking_contract,
    benchmark_to_conditions,
    benchmark_to_formulation,
    get_matrix_only_target_snapshot_exclusions,
)

# ── Evaluator (prediction + statistics) ──────────────────────────────────────
from src.benchmark_evaluator import (
    MATRIX_BENCHMARK_PROFILES,
    MATRIX_BENCHMARK_BASE_MARKER_YIELDS,
    BENCHMARK_NAME_ALIASES,
    THERMODYNAMIC_GATING_THRESHOLD_KCAL,
    evaluate_benchmark,
    evaluate_benchmark_payload,
    _best_prediction_match,
    _build_comparisons,
    _pearson,
    _mean_abs_log10_error,
    _resolve_scale_thresholds,
    _run_benchmark_recommendation,
    _run_matrix_only_benchmark_prediction,
)

# ── Reporting (summaries, artifacts, markdown, snapshots) ────────────────────
from src.benchmark_reporting import (
    DEFAULT_BENCHMARK_THRESHOLDS,
    THERMODYNAMIC_GATING_MIN_ABSOLUTE_MAE_IMPROVEMENT_PPB,
    THERMODYNAMIC_GATING_MIN_RELATIVE_MAE_IMPROVEMENT,
    THERMODYNAMIC_GATING_MIN_RATIO_IMPROVEMENT,
    summarize_evaluation_for_benchmark,
    summarize_evaluation,
    summarize_benchmarks,
    assess_matrix_benchmark_evidence,
    build_matrix_benchmark_evidence_audit,
    build_matrix_benchmark_deltas,
    compare_matrix_benchmark_delta_sets,
    render_matrix_branch_deltas_markdown,
    build_family_lane_validation_artifact,
    render_family_lane_validation_markdown,
    thermodynamic_gating_materially_improves,
    audit_thermodynamic_gating,
    audit_all_thermodynamic_gating,
    snapshot_benchmark_targets,
    snapshot_all_benchmark_targets,
    build_benchmark_index,
    _projection_metadata_for_match,
    _evaluate_matrix_ranking_contract,
    _enrich_benchmark_summary_family_metadata,
    _benchmark_status_score,
)

# ── Assertions ────────────────────────────────────────────────────────────────
from src.benchmark_assertions import (
    build_matrix_benchmark_assertions,
    render_matrix_benchmark_assertions_markdown,
    _matrix_assertion_thresholds,
    _predicted_order_lookup,
    _matched_contract_prediction_rows,
)

# ── Re-export the large reporting functions that remain in this file for now ──
# (The "build_matrix_target_status_artifact", "build_matrix_promotion_*", and
#  "build_matrix_observable_closure_audit" functions are heavy and have many
#  private helpers.  They are re-exported from the reporting module below.)
from src.benchmark_reporting import (  # noqa: F811  (re-export)
    _benchmark_compound_names,
    _payload_role_from_evidence_state,
)

# Inline re-exports that were defined as private helpers in the original file
# and are referenced by scripts/generators.
from src.benchmark_registry import _is_supported_formulation  # noqa: F401
from src.benchmark_evaluator import _tokenize_name  # noqa: F401

# ── Compatibility shims for build_matrix_target_status_artifact family ────────
from src.benchmark_reporting import (  # noqa: F401
    build_matrix_target_status_artifact,
    build_matrix_promotion_family_status,
    build_matrix_promotion_contract_artifact,
    build_matrix_observable_closure_audit,
    _matrix_compound_support_status,
    _select_matrix_promotion_target,
    _matrix_promotion_requirement_rows,
    _matrix_closure_action,
    _mechanistic_refinement_expected_change,
    _matrix_source_origin,
    _matrix_source_reference,
    _matrix_target_profile,
    _matrix_external_data_status,
)