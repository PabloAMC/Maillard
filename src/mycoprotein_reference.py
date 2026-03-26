from __future__ import annotations

from typing import Any, Dict, Mapping

from src.matrix_family_coverage import build_matrix_family_coverage_artifact
from src.matrix_family_next_action import build_matrix_family_next_action_artifact
from src.matrix_prior_registry import (
    get_accessibility_window_entry,
    get_denaturation_heuristic_entry,
    get_volatile_class_profile_entry,
    summarize_matrix_prior_bundle,
)


def build_mycoprotein_reference_artifact() -> Dict[str, Any]:
    coverage = build_matrix_family_coverage_artifact()
    next_action = build_matrix_family_next_action_artifact()
    family_row = next(
        row for row in coverage.get("families", []) if str(row.get("matrix_family", "")) == "mycoprotein"
    )
    next_row = next(
        row for row in next_action.get("decisions", []) if str(row.get("matrix_family", "")) == "mycoprotein"
    )

    accessibility = get_accessibility_window_entry("myco") or {}
    denaturation = get_denaturation_heuristic_entry("myco") or {}
    volatile_profile = get_volatile_class_profile_entry("myco") or {}

    return {
        "summary": {
            "matrix_family": "mycoprotein",
            "protein_type": "myco",
            "runtime_posture": str(family_row.get("runtime_posture", "unknown")),
            "evidence_surface": str(family_row.get("evidence_surface", "unknown")),
            "decision": str(next_row.get("decision", "unknown")),
            "next_best_action": str(next_row.get("next_best_action", family_row.get("next_best_action", "unknown"))),
            "policy": "mycoprotein_reference_is_runtime_executable_but_remains_bounded_directional_support_until_family_specific_measurements_are_landed",
        },
        "family_scope": {
            "supported_today": list(family_row.get("what_is_supported", [])),
            "not_supported_today": list(family_row.get("what_is_not_supported", [])),
            "primary_blocker": str(family_row.get("primary_blocker", "none")),
        },
        "matrix_prior_bundle": summarize_matrix_prior_bundle("myco"),
        "reference_windows": {
            "accessibility": {
                "lysine_min": accessibility.get("lysine_min"),
                "lysine_max": accessibility.get("lysine_max"),
                "cysteine_min": accessibility.get("cysteine_min"),
                "cysteine_max": accessibility.get("cysteine_max"),
                "uncertainty_posture": accessibility.get("uncertainty_posture"),
                "source": accessibility.get("source"),
            },
            "denaturation": {
                "midpoint_celsius": denaturation.get("midpoint_celsius"),
                "width_celsius": denaturation.get("width_celsius"),
                "time_gain_celsius": denaturation.get("time_gain_celsius"),
                "uncertainty_posture": denaturation.get("uncertainty_posture"),
                "source": denaturation.get("source"),
            },
            "volatile_classes": {
                "native_factors": dict(volatile_profile.get("native_factors", {})),
                "denatured_factors": dict(volatile_profile.get("denatured_factors", {})),
                "uncertainty_posture": volatile_profile.get("uncertainty_posture"),
                "source": volatile_profile.get("source"),
            },
        },
    }


def render_mycoprotein_reference_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    scope = payload.get("family_scope", {})
    windows = payload.get("reference_windows", {})
    acc = windows.get("accessibility", {})
    den = windows.get("denaturation", {})
    vol = windows.get("volatile_classes", {})
    lines = [
        "# Mycoprotein Reference",
        "",
        f"Matrix family: {summary.get('matrix_family', 'unknown')}",
        f"Runtime posture: {summary.get('runtime_posture', 'unknown')}",
        f"Evidence surface: {summary.get('evidence_surface', 'unknown')}",
        f"Decision: {summary.get('decision', 'unknown')}",
        f"Next best action: {summary.get('next_best_action', 'unknown')}",
        f"Policy: {summary.get('policy', 'unknown')}",
        "",
        "## Bounded Reference Windows",
        "",
        f"Accessibility window: lysine {acc.get('lysine_min', 'n/a')} to {acc.get('lysine_max', 'n/a')}; cysteine {acc.get('cysteine_min', 'n/a')} to {acc.get('cysteine_max', 'n/a')}",
        f"Denaturation heuristic: midpoint {den.get('midpoint_celsius', 'n/a')} C; width {den.get('width_celsius', 'n/a')} C; time gain {den.get('time_gain_celsius', 'n/a')} C",
        f"Volatile profile uncertainty posture: {vol.get('uncertainty_posture', 'unknown')}",
        "",
        "## Scope",
        "",
        f"Supported today: {', '.join(str(item) for item in scope.get('supported_today', [])) or 'none'}",
        f"Not supported today: {', '.join(str(item) for item in scope.get('not_supported_today', [])) or 'none'}",
        f"Primary blocker: {scope.get('primary_blocker', 'none')}",
    ]
    return "\n".join(lines) + "\n"