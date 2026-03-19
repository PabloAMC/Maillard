from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
from src.projection_metadata import ProjectionMetadataMap


@dataclass(frozen=True)
class BenchmarkMetadata:
    tier: str
    family: str
    execution_path: str
    benchmark_engine: str
    comparator_signal: str
    cantera_role: str
    target_snapshot_policy: str
    thermodynamic_gating_policy: str
    notes: Optional[str] = None


@dataclass(frozen=True)
class BenchmarkIndexEntry:
    benchmark_id: str
    bench_file: Path
    tier: str
    family: str
    protein_type: str
    execution_path: str
    benchmark_engine: str
    cantera_role: str
    thermodynamic_gating_policy: str
    supported: bool
    reason: Optional[str]
    status: str
    strict_ready: bool
    process_state: Optional[str] = None
    ranking_contract_status: str = "n/a"


class BenchmarkNotSupportedError(RuntimeError):
    pass


@dataclass(frozen=True)
class CompoundComparison:
    compound: str
    measured_ppb: float
    predicted_ppb: float
    matched_name: Optional[str]
    uncertainty_pct: Optional[float]
    match_score: float = 0.0

    @property
    def ratio(self) -> float:
        smallest = min(self.measured_ppb, self.predicted_ppb)
        largest = max(self.measured_ppb, self.predicted_ppb)
        if smallest <= 0.0:
            return math.inf if largest > 0.0 else 1.0
        return largest / smallest


@dataclass(frozen=True)
class BenchmarkEvaluation:
    benchmark_id: str
    bench_file: Path
    supported: bool
    reason: Optional[str]
    predicted_ppb: Dict[str, float]
    comparisons: List[CompoundComparison]
    pearson_r: Optional[float]
    mae_ppb: Optional[float]
    projection_metadata: ProjectionMetadataMap = field(default_factory=dict)
    reference_signal_origin: str = "measured_volatiles"

    @property
    def coverage(self) -> float:
        if not self.comparisons:
            return 0.0
        matched = sum(1 for comparison in self.comparisons if comparison.matched_name is not None)
        return matched / len(self.comparisons)


@dataclass(frozen=True)
class BenchmarkSummary:
    benchmark_id: str
    bench_file: Path
    tier: str
    family: str
    execution_path: str
    benchmark_engine: str
    comparator_signal: str
    cantera_role: str
    target_snapshot_policy: str
    thermodynamic_gating_policy: str
    supported: bool
    reason: Optional[str]
    protein_type: str
    coverage: float
    matched_compounds: int
    total_compounds: int
    pearson_r: Optional[float]
    mae_ppb: Optional[float]
    max_ratio: Optional[float]
    mean_ratio: Optional[float]
    ranking_status: str
    scale_status: str
    overall_status: str
    strict_ready: bool
    blocking_issues: List[str]
    conditions: Dict[str, float]
    process_state: Optional[str] = None
    ranked_observable_targets: List[str] = field(default_factory=list)
    adverse_markers: List[str] = field(default_factory=list)
    ranking_contract_status: str = "n/a"
    calibration_mode: Optional[str] = None
    reference_signal_origin: str = "measured_volatiles"
    mean_abs_log10_error: Optional[float] = None


@dataclass(frozen=True)
class MatrixBenchmarkDelta:
    benchmark_id: str
    bench_file: Path
    protein_type: str
    execution_path: str
    process_state: Optional[str]
    reference_signal_origin: str
    ranking_contract_status: str
    compound: str
    role: str
    reference_ppb: float
    predicted_ppb: float
    abs_delta_ppb: float
    pct_delta: Optional[float]
    ratio: float
    calibration_source: str
    calibration_evidence_strength: str
    calibration_fallback_mode: str


@dataclass(frozen=True)
class MatrixBenchmarkEvidence:
    benchmark_id: str
    bench_file: Path
    protein_type: str
    execution_path: str
    process_state: Optional[str]
    reference_signal_origin: str
    source_origin: str
    source_reference: str
    target_profile: str
    external_data_status: str
    promotable: bool
    promotion_blocker: str


@dataclass(frozen=True)
class MatrixPromotionFamilyStatus:
    protein_type: str
    off_flavour_anchor_count: int
    meaty_candidate_count: int
    external_meaty_anchor_count: int
    candidate_set_ready: bool
    external_assessment_unlocked: bool
    blocker: str


@dataclass(frozen=True)
class MatrixBenchmarkBranchDelta:
    benchmark_id: str
    compound: str
    change_type: str
    current_present: bool
    baseline_present: bool
    current_execution_path: str
    baseline_execution_path: str
    current_reference_signal_origin: str
    baseline_reference_signal_origin: str
    current_source_origin: str
    baseline_source_origin: str
    current_external_data_status: str
    baseline_external_data_status: str
    current_predicted_ppb: Optional[float] = None
    baseline_predicted_ppb: Optional[float] = None
    predicted_delta_ppb: Optional[float] = None
    current_ratio: Optional[float] = None
    baseline_ratio: Optional[float] = None
    ratio_delta: Optional[float] = None


@dataclass(frozen=True)
class ThermodynamicGatingAudit:
    benchmark_id: str
    bench_file: Path
    execution_path: str
    applicable: bool
    baseline_overall_status: str
    gated_overall_status: str
    baseline_mae_ppb: Optional[float] = None
    gated_mae_ppb: Optional[float] = None
    baseline_max_ratio: Optional[float] = None
    gated_max_ratio: Optional[float] = None
    delta_mae_ppb: Optional[float] = None
    delta_max_ratio: Optional[float] = None
    material_improvement: bool = False
    recommended_policy: str = "diagnostic_only"
    notes: str = ""


@dataclass(frozen=True)
class MatrixBenchmarkAssertion:
    benchmark_id: str
    bench_file: Path
    protein_type: str
    execution_path: str
    process_state: Optional[str]
    target_profile: str
    ranking_contract_status: str
    coverage: float
    min_coverage: float
    top_k: int
    top_k_hits: int
    top_k_status: str
    adverse_order_status: str
    max_ratio: Optional[float]
    ratio_tolerance: float
    ratio_status: str
    overall_status: str
    strict_gate_blocked: bool
    blocker: str


@dataclass(frozen=True)
class BenchmarkTargetSnapshot:
    benchmark_id: str
    bench_file: Path
    target_name: str
    target_type: str
    roles: List[str]
    predicted_ppb: float
    proxy_ppb: float
    observable_ratio: float
    weighted_flux: float
    span: float
    depth: int
    volatile_class: str
    matrix_factor: float
    headspace_factor: float
    headspace_observable: bool
    headspace_class: str
    henry_kaw_25c: Optional[float] = None
    henry_source_name: Optional[str] = None
