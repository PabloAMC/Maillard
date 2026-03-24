from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_MATRIX_FAMILY_COVERAGE_REGISTRY = _repo_root() / "data" / "lit" / "matrix_family_coverage_registry.json"


def load_matrix_family_coverage_registry(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_MATRIX_FAMILY_COVERAGE_REGISTRY
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def build_matrix_family_coverage_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    payload = load_matrix_family_coverage_registry(file_path)
    families = [dict(row) for row in payload.get("families", [])]

    posture_counts: Dict[str, int] = {}
    category_counts: Dict[str, int] = {}
    for row in families:
        posture = str(row.get("runtime_posture", "unknown"))
        category = str(row.get("category", "unknown"))
        posture_counts[posture] = posture_counts.get(posture, 0) + 1
        category_counts[category] = category_counts.get(category, 0) + 1

    explicit_supported = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) in {"quantitative_core", "directional_matrix", "qualitative_intake_only"}
    ]
    indirect_only = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) == "indirect_generic_support"
    ]
    open_gaps = [
        row["matrix_family"]
        for row in families
        if str(row.get("runtime_posture", "")) == "open_gap"
    ]

    return {
        "summary": {
            "family_count": len(families),
            "posture_counts": dict(sorted(posture_counts.items())),
            "category_counts": dict(sorted(category_counts.items())),
            "explicit_supported_families": explicit_supported,
            "indirect_only_families": indirect_only,
            "open_gap_families": open_gaps,
            "policy": "matrix_family_scope_must_distinguish_explicit_support_from_generic_indirect_support",
        },
        "families": families,
    }


def render_matrix_family_coverage_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Matrix Family Coverage",
        "",
        "| Matrix Family | Category | Runtime Posture | Evidence Surface | Importance | Next Best Action |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {row.get('category', 'unknown')} | {row.get('runtime_posture', 'unknown')} | "
            f"{row.get('evidence_surface', 'unknown')} | {row.get('importance_tier', 'unknown')} | {row.get('next_best_action', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Scope Notes",
        "",
        "| Matrix Family | Supported Today | Not Supported Today |",
        "| --- | --- | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {'; '.join(str(item) for item in row.get('what_is_supported', [])) or 'none'} | "
            f"{'; '.join(str(item) for item in row.get('what_is_not_supported', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Matrix families tracked: {int(summary.get('family_count', 0))}",
            f"Explicitly supported families: {', '.join(str(item) for item in summary.get('explicit_supported_families', [])) or 'none'}",
            f"Indirect-only families: {', '.join(str(item) for item in summary.get('indirect_only_families', [])) or 'none'}",
            f"Open-gap families: {', '.join(str(item) for item in summary.get('open_gap_families', [])) or 'none'}",
            f"Policy: {summary.get('policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"