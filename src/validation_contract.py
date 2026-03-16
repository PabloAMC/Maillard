from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class BenchmarkThresholds:
    min_matched_for_ranking: int = 3
    ranking_threshold: float = 0.85
    full_coverage_threshold: float = 1.0
    free_aa_ratio_threshold: float = 1.5
    matrix_ratio_threshold: float = 2.0

    def ratio_threshold_for(self, protein_type: str) -> float:
        return self.free_aa_ratio_threshold if protein_type == "free" else self.matrix_ratio_threshold


@dataclass(frozen=True)
class ValidationContract:
    replication_meaning: str
    directional_validity: str
    quantitative_replication: str
    formulation_utility: str
    strict_gate_execution_paths: tuple[str, ...] = ("free_precursor",)
    strict_gate_tiers: tuple[str, ...] = ("PRIMARY",)
    thresholds: BenchmarkThresholds = field(default_factory=BenchmarkThresholds)

    def is_strict_gate_eligible(self, *, tier: str, execution_path: str) -> bool:
        return tier in self.strict_gate_tiers and execution_path in self.strict_gate_execution_paths


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
        "Quantitative replication asks whether supported benchmarks also meet the configured Pearson and ratio "
        "thresholds at full coverage."
    ),
    formulation_utility=(
        "Formulation utility is downstream usefulness for ranking recipes or interventions; it may remain useful "
        "even when a benchmark is only directionally valid and not strict-ready."
    ),
)