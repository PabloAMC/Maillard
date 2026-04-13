from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src.mlp_adoption_contract import MLPModelCandidate, load_mlp_candidates


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_MLP_EXTERNAL_EVIDENCE = _repo_root() / "data" / "lit" / "mlp_external_benchmark_evidence.json"


@dataclass(frozen=True)
class ExternalMLPEvidence:
    evidence_id: str
    model_family: str
    model_name: str
    prior_strength: str
    domain_relevance: str
    benchmark_scope: str
    supported_tasks: List[str]
    supported_roles: List[str]
    unproven_roles: List[str]
    why_relevant: str
    caution: str
    source_notes: List[str]
    sources: List[str]


def load_external_mlp_evidence(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_MLP_EXTERNAL_EVIDENCE
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_external_mlp_evidence_entries(file_path: Optional[Path | str] = None) -> List[ExternalMLPEvidence]:
    payload = load_external_mlp_evidence(file_path)
    return [ExternalMLPEvidence(**row) for row in payload.get("models", [])]


def build_external_mlp_evidence_index(file_path: Optional[Path | str] = None) -> Dict[str, ExternalMLPEvidence]:
    entries = load_external_mlp_evidence_entries(file_path)
    return {entry.evidence_id: entry for entry in entries}


def _selection_priority(candidate: MLPModelCandidate, evidence: Optional[ExternalMLPEvidence]) -> str:
    if evidence is None:
        return "unknown"
    if candidate.proposed_role in evidence.supported_roles and evidence.prior_strength == "strong":
        return "high"
    if candidate.proposed_role in evidence.supported_roles and evidence.prior_strength == "moderate":
        return "medium"
    if evidence.prior_strength in {"strong", "moderate"}:
        return "watch"
    return "low"


def build_external_mlp_landscape_payload() -> Dict[str, Any]:
    candidates = load_mlp_candidates()
    evidence_index = build_external_mlp_evidence_index()
    rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        evidence = evidence_index.get(str(candidate.external_evidence_id or ""))
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "model_family": candidate.model_family,
                "model_name": candidate.model_name,
                "proposed_role": candidate.proposed_role,
                "external_evidence_id": candidate.external_evidence_id,
                "prior_strength": evidence.prior_strength if evidence else "unknown",
                "domain_relevance": evidence.domain_relevance if evidence else "unknown",
                "benchmark_scope": evidence.benchmark_scope if evidence else "unknown",
                "supported_roles": [] if evidence is None else list(evidence.supported_roles),
                "unproven_roles": [] if evidence is None else list(evidence.unproven_roles),
                "selection_priority": _selection_priority(candidate, evidence),
                "why_relevant": "" if evidence is None else evidence.why_relevant,
                "caution": "" if evidence is None else evidence.caution,
                "sources": [] if evidence is None else list(evidence.sources),
            }
        )

    return {
        "summary": {
            "candidate_count": len(rows),
            "policy": "external_benchmarks_prioritize_candidates_but_do_not_replace_local_maillard_validation",
        },
        "candidates": rows,
        "evidence": [asdict(entry) for entry in evidence_index.values()],
    }


def render_external_mlp_landscape_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# P4 External MLP Landscape",
        "",
        "| Candidate | Role | External Prior | Domain Relevance | Selection Priority | Externally Supported | Still Unproven |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("candidates", []):
        lines.append(
            f"| {row.get('candidate_id', 'unknown')} | {row.get('proposed_role', 'unknown')} | {row.get('prior_strength', 'unknown')} | "
            f"{row.get('domain_relevance', 'unknown')} | {row.get('selection_priority', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('supported_roles', [])) or 'none'} | "
            f"{', '.join(str(item) for item in row.get('unproven_roles', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Candidates reviewed: {int(summary.get('candidate_count', 0))}",
            f"Policy: {summary.get('policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"