from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class BenchmarkThresholds:
    min_matched_for_ranking: int = 3
    ranking_threshold: float = 0.85
    full_coverage_threshold: float = 1.0
    free_aa_ratio_threshold: float = 1.5
    matrix_ratio_threshold: float = 2.0
    free_aa_mean_abs_log10_error_threshold: float = 0.10
    matrix_mean_abs_log10_error_threshold: float = 0.12

    def ratio_threshold_for(self, protein_type: str) -> float:
        return self.free_aa_ratio_threshold if protein_type == "free" else self.matrix_ratio_threshold

    def mean_abs_log10_error_threshold_for(self, protein_type: str) -> float:
        return (
            self.free_aa_mean_abs_log10_error_threshold
            if protein_type == "free"
            else self.matrix_mean_abs_log10_error_threshold
        )


@dataclass(frozen=True)
class BenchmarkExecutionPolicy:
    execution_path: str
    benchmark_engine: str
    comparator_signal: str
    cantera_role: str
    target_snapshot_policy: str
    thermodynamic_gating_policy: str
    notes: str


@dataclass(frozen=True)
class ValidationContract:
    replication_meaning: str
    directional_validity: str
    quantitative_replication: str
    formulation_utility: str
    benchmark_policy: str
    strict_gate_execution_paths: tuple[str, ...] = ("free_precursor",)
    strict_gate_tiers: tuple[str, ...] = ("PRIMARY",)
    execution_policies: tuple[BenchmarkExecutionPolicy, ...] = field(default_factory=tuple)
    thresholds: BenchmarkThresholds = field(default_factory=BenchmarkThresholds)

    def is_strict_gate_eligible(self, *, tier: str, execution_path: str) -> bool:
        return tier in self.strict_gate_tiers and execution_path in self.strict_gate_execution_paths

    def policy_for_execution_path(self, execution_path: str) -> BenchmarkExecutionPolicy:
        for policy in self.execution_policies:
            if policy.execution_path == execution_path:
                return policy
        raise KeyError(f"No benchmark execution policy registered for execution_path={execution_path!r}")


DEFAULT_VALIDATION_CONTRACT = ValidationContract(
    replication_meaning=(
        "Literature replication means executing a benchmark through the intended model path and "
        "matching the reported experimental system with explicit coverage, ranking, and scale criteria."
    ),
    directional_validity=(
        "Directional validity asks whether the model preserves the experimentally observed ordering or trend, "
        "without claiming quantitative concentration accuracy."
    ),
    quantitative_replication=(
        "Quantitative replication asks whether supported benchmarks also meet the configured Pearson, ratio, and "
        "log-scale error thresholds at full coverage."
    ),
    formulation_utility=(
        "Formulation utility is downstream usefulness for ranking recipes or interventions; it may remain useful "
        "even when a benchmark is only directionally valid and not strict-ready."
    ),
    benchmark_policy=(
        "PRIMARY free-precursor benchmarks compare the FAST observable projection as the authoritative regression signal; "
        "Cantera remains a diagnostic reference lane until a benchmark-facing dual-lane contract is explicitly promoted. "
        "Matrix-only benchmarks compare the dedicated intake/headspace path and stay outside the strict gate. "
        "Matrix precursor-augmented benchmarks are allowed as reproducible matrix target-ranking harnesses, but they also remain outside the strict gate until tied to external evidence."
    ),
    execution_policies=(
        BenchmarkExecutionPolicy(
            execution_path="free_precursor",
            benchmark_engine="fast_observable",
            comparator_signal="predicted_ppb",
            cantera_role="diagnostic_reference_only",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="diagnostic_only",
            notes="FAST observable projection is the benchmark-facing contract; Cantera is used for diagnosis and calibration, not for strict-gate pass/fail. Thermodynamic gating stays diagnostic-only unless a benchmark is explicitly promoted after a materiality audit.",
        ),
        BenchmarkExecutionPolicy(
            execution_path="matrix_only",
            benchmark_engine="matrix_intake_headspace",
            comparator_signal="predicted_ppb",
            cantera_role="not_authoritative",
            target_snapshot_policy="excluded",
            thermodynamic_gating_policy="not_applicable",
            notes="Matrix-only benchmarks validate the dedicated intake/headspace path and remain outside the strict gate until a precursor-resolved benchmark model exists.",
        ),
        BenchmarkExecutionPolicy(
            execution_path="matrix_precursor_augmented",
            benchmark_engine="fast_observable_matrix",
            comparator_signal="predicted_ppb",
            cantera_role="not_authoritative",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="not_applicable",
            notes="Matrix precursor-augmented benchmarks run the FAST observable path with non-free protein accessibility enabled; they are useful for reproducible matrix ranking harnesses but remain outside the strict gate until externally benchmarked.",
        ),
    ),
)