from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from src.artifact_io import load_json_mapping, repo_root, resolve_optional_path


DEFAULT_CHEMISTRY_FAMILY_SCOPE_REGISTRY = repo_root() / "data" / "lit" / "chemistry_family_scope_registry.json"


def load_chemistry_family_scope_registry(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    return load_json_mapping(resolve_optional_path(file_path, DEFAULT_CHEMISTRY_FAMILY_SCOPE_REGISTRY))


def build_chemistry_family_scope_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    payload = load_chemistry_family_scope_registry(file_path)
    families = [dict(row) for row in payload.get("families", [])]

    status_counts: Dict[str, int] = {}
    role_counts: Dict[str, int] = {}
    for row in families:
        status = str(row.get("current_status", "unknown"))
        role = str(row.get("strategic_role", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1
        role_counts[role] = role_counts.get(role, 0) + 1

    first_class = [
        row["family_id"]
        for row in families
        if str(row.get("current_status", "")) == "first_class_core"
    ]
    expansion_candidates = [
        row["family_id"]
        for row in families
        if str(row.get("current_status", "")) in {"partially_encoded_high_priority", "bounded_lane"}
    ]
    open_gaps = [
        row["family_id"]
        for row in families
        if str(row.get("current_status", "")) == "open_gap"
    ]

    recommended_next_family = ""
    for row in families:
        if str(row.get("current_status", "")) == "partially_encoded_high_priority":
            recommended_next_family = str(row.get("family_id", ""))
            break

    return {
        "summary": {
            "family_count": len(families),
            "status_counts": dict(sorted(status_counts.items())),
            "strategic_role_counts": dict(sorted(role_counts.items())),
            "first_class_families": first_class,
            "expansion_candidates": expansion_candidates,
            "open_gap_families": open_gaps,
            "recommended_next_family": recommended_next_family,
            "policy": "expand_product_scope_by_adding_benchmark_visible_family_lanes_not_by_diluting_the_core",
            "ingestion_policy": "beyond_amino_acid_sugar_use_the_same_runtime_contract_but_choose_family_specific_payload_types_benchmark_payloads_for_observable_panels_priors_for_transfer_and_gap_registries_for_non_closable_scope",
        },
        "families": families,
    }


def render_chemistry_family_scope_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Chemistry Family Scope",
        "",
        "| Family | Status | Strategic Role | Preferred Ingestion Mode | Priority | Next Best Action |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('family_id', 'unknown')} | {row.get('current_status', 'unknown')} | {row.get('strategic_role', 'unknown')} | "
            f"{row.get('preferred_ingestion_mode', 'unknown')} | {row.get('priority', 'unknown')} | {row.get('next_best_action', 'unknown')} |"
        )

    lines.extend([
        "",
        "## Why These Families Matter",
        "",
        "| Family | Why It Matters | Current Runtime Assets | Ingestion Surfaces |",
        "| --- | --- | --- | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('family_id', 'unknown')} | {row.get('why_it_matters', 'unknown')} | "
            f"{'; '.join(str(item) for item in row.get('current_runtime_assets', [])) or 'none'} | "
            f"{'; '.join(str(item) for item in row.get('ingestion_surfaces', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Chemistry families tracked: {int(summary.get('family_count', 0))}",
            f"First-class families: {', '.join(str(item) for item in summary.get('first_class_families', [])) or 'none'}",
            f"Expansion candidates: {', '.join(str(item) for item in summary.get('expansion_candidates', [])) or 'none'}",
            f"Open-gap families: {', '.join(str(item) for item in summary.get('open_gap_families', [])) or 'none'}",
            f"Recommended next family: {summary.get('recommended_next_family', 'none')}",
            f"Policy: {summary.get('policy', 'unknown')}",
            f"Ingestion policy: {summary.get('ingestion_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"