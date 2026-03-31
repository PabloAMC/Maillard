from __future__ import annotations

from typing import Any, Dict, List, Mapping

from src.matrix_family_priority_ranking import build_matrix_family_priority_ranking_artifact


def build_matrix_family_next_action_artifact() -> Dict[str, Any]:
    ranking = build_matrix_family_priority_ranking_artifact()
    families = list(ranking.get("families", []))

    bounded_candidates = [
        row for row in families if str(row.get("scope_priority", "")) == "bounded_next_candidate"
    ]
    chosen = bounded_candidates[0] if bounded_candidates else None

    decision_rows: List[Dict[str, Any]] = []
    for row in families:
        family = str(row.get("matrix_family", "unknown"))
        scope_priority = str(row.get("scope_priority", "unknown"))
        if chosen is not None and family == str(chosen.get("matrix_family", "")):
            decision = "advance_now"
            rationale = "Highest-ranked bounded next candidate with a calibration-grade landing that does not overclaim broad matrix closure."
        elif scope_priority == "active_matrix_priority":
            decision = "defer_until_primary_matrix_lane_moves"
            rationale = "Pea and soy remain the main product blocker, so these stay primary closure work rather than scope-expansion work."
        elif scope_priority == "scope_gap_to_rank":
            decision = "defer_as_scope_gap"
            rationale = "Family remains explicitly out of scope until family-specific evidence landing is available."
        else:
            decision = "hold_current_posture"
            rationale = "Current support is intentionally bounded and should not be broadened in this cycle."

        decision_rows.append(
            {
                "matrix_family": family,
                "priority_rank": int(row.get("priority_rank", 0)),
                "scope_priority": scope_priority,
                "evidence_landing": str(row.get("evidence_landing", "unknown")),
                "decision": decision,
                "next_best_action": str(row.get("next_best_action", "unknown")),
                "rationale": rationale,
            }
        )

    return {
        "summary": {
            "chosen_family": None if chosen is None else str(chosen.get("matrix_family", "unknown")),
            "chosen_evidence_landing": None if chosen is None else str(chosen.get("evidence_landing", "unknown")),
            "chosen_action": "none" if chosen is None else "advance_bounded_next_family",
            "policy": "advance_at_most_one_bounded_next_family_while_explicitly_deferring_scope_gaps_and_preserving_primary_pea_soy_closure_work",
        },
        "decisions": decision_rows,
    }


def render_matrix_family_next_action_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    lines = [
        "# Matrix Family Next Action",
        "",
        f"Chosen family: {summary.get('chosen_family', 'none')}",
        f"Chosen evidence landing: {summary.get('chosen_evidence_landing', 'none')}",
        f"Chosen action: {summary.get('chosen_action', 'none')}",
        f"Policy: {summary.get('policy', 'unknown')}",
        "",
        "| Rank | Matrix Family | Scope Priority | Evidence Landing | Decision | Next Best Action |",
        "| ---: | --- | --- | --- | --- | --- |",
    ]
    for row in payload.get("decisions", []):
        lines.append(
            f"| {int(row.get('priority_rank', 0))} | {row.get('matrix_family', 'unknown')} | {row.get('scope_priority', 'unknown')} | {row.get('evidence_landing', 'unknown')} | {row.get('decision', 'unknown')} | {row.get('next_best_action', 'unknown')} |"
        )

    lines.extend(["", "## Rationale", ""])
    for row in payload.get("decisions", []):
        lines.append(f"- {row.get('matrix_family', 'unknown')}: {row.get('rationale', 'unknown')}")
    return "\n".join(lines) + "\n"