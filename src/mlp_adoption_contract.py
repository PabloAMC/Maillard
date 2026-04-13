from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_MLP_CANDIDATE_REGISTRY = _repo_root() / "data" / "lit" / "mlp_candidate_registry.json"


OFFLINE_ACCELERATOR_ROLES = {
    "conformer_screening",
    "geom_preopt",
    "ts_initialization",
    "local_barrier_surrogate",
}


@dataclass(frozen=True)
class MLPModelCandidate:
    candidate_id: str
    model_family: str
    model_name: str
    chemistry_domain: str
    materials_first: bool
    proposed_role: str
    target_motif_families: List[str]
    benchmark_visible_gap: str
    expected_speedup: float
    likely_failure_modes: List[str]
    fallback_comparator: str
    backend_locator: Optional[str] = None
    benchmark_results_path: Optional[str] = None
    geometry_benchmark_path: Optional[str] = None
    ts_seed_benchmark_path: Optional[str] = None
    external_evidence_id: Optional[str] = None
    status: str = "candidate_shortlist"


@dataclass(frozen=True)
class MLPAdoptionDecision:
    candidate_id: str
    model_family: str
    model_name: str
    proposed_role: str
    decision: str
    benchmark_set_id: str
    coverage_ratio: float
    rank_correlation: Optional[float]
    mean_abs_error_kcal: Optional[float]
    max_abs_error_kcal: Optional[float]
    stop_reasons: List[str]
    rationale: str
    fallback_comparator: str
    benchmark_visible_gap: str
    approved_for_default: bool = False


def load_mlp_candidate_registry(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_MLP_CANDIDATE_REGISTRY
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_mlp_candidates(file_path: Optional[Path | str] = None) -> List[MLPModelCandidate]:
    payload = load_mlp_candidate_registry(file_path)
    supported_fields = {field.name for field in fields(MLPModelCandidate)}
    candidates = [
        MLPModelCandidate(**{key: value for key, value in row.items() if key in supported_fields})
        for row in payload.get("candidates", [])
    ]
    for candidate in candidates:
        if candidate.proposed_role not in OFFLINE_ACCELERATOR_ROLES:
            raise ValueError(f"Unsupported MLP role for {candidate.candidate_id}: {candidate.proposed_role}")
    return candidates


def build_adoption_note_payload(
    decisions: List[MLPAdoptionDecision],
    *,
    benchmark_set_id: str,
) -> Dict[str, Any]:
    adopted = [decision for decision in decisions if decision.decision == "adopt_offline"]
    return {
        "summary": {
            "benchmark_set_id": benchmark_set_id,
            "decision_count": len(decisions),
            "adopted_count": len(adopted),
            "default_policy": "no_default_mlp_until_reaction_benchmark_passes",
        },
        "decisions": [asdict(decision) for decision in decisions],
    }


def render_adoption_note_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# MLP Adoption Notes",
        "",
        "| Candidate | Role | Decision | Coverage | Rank Correlation | MAE (kcal/mol) | Stop Reasons | Fallback |",
        "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in payload.get("decisions", []):
        correlation = row.get("rank_correlation")
        mae = row.get("mean_abs_error_kcal")
        lines.append(
            f"| {row.get('candidate_id', 'unknown')} | {row.get('proposed_role', 'unknown')} | {row.get('decision', 'unknown')} | "
            f"{100.0 * float(row.get('coverage_ratio', 0.0)):.0f}% | "
            f"{'' if correlation is None else f'{float(correlation):.2f}'} | "
            f"{'' if mae is None else f'{float(mae):.2f}'} | "
            f"{', '.join(str(item) for item in row.get('stop_reasons', [])) or 'none'} | {row.get('fallback_comparator', 'unknown')} |"
        )
    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Benchmark set: {summary.get('benchmark_set_id', 'unknown')}",
            f"Decisions recorded: {int(summary.get('decision_count', 0))}",
            f"Adopted candidates: {int(summary.get('adopted_count', 0))}",
            f"Default policy: {summary.get('default_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"