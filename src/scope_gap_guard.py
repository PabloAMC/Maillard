from __future__ import annotations

from typing import Any, Dict, Mapping

from src.matrix_family_priority_ranking import build_matrix_family_priority_ranking_artifact


def build_scope_gap_guard_artifact() -> Dict[str, Any]:
    ranking = build_matrix_family_priority_ranking_artifact()
    guard_rows = []
    for row in ranking.get("families", []):
        if str(row.get("scope_priority", "")) != "scope_gap_to_rank":
            continue
        guard_rows.append(
            {
                "matrix_family": str(row.get("matrix_family", "unknown")),
                "decision": "blocked_from_active_expansion",
                "evidence_landing_required": str(row.get("evidence_landing", "unknown")),
                "primary_blocker": str(row.get("primary_blocker", "unknown")),
                "next_best_action": str(row.get("next_best_action", "unknown")),
            }
        )

    return {
        "summary": {
            "blocked_family_count": len(guard_rows),
            "blocked_families": [row["matrix_family"] for row in guard_rows],
            "policy": "scope_gap_families_must_remain_out_of_active_expansion_until_family_specific_evidence_landing_exists",
        },
        "guards": guard_rows,
    }


def render_scope_gap_guard_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Scope Gap Guard",
        "",
        "| Matrix Family | Guard Decision | Evidence Landing Required | Primary Blocker | Next Best Action |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("guards", []):
        lines.append(
            f"| {row.get('matrix_family', 'unknown')} | {row.get('decision', 'unknown')} | {row.get('evidence_landing_required', 'unknown')} | {row.get('primary_blocker', 'unknown')} | {row.get('next_best_action', 'unknown')} |"
        )
    lines.append("")
    lines.append(f"Policy: {payload.get('summary', {}).get('policy', 'unknown')}")
    return "\n".join(lines) + "\n"