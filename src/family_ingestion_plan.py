from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


DEFAULT_FAMILY_INGESTION_PLAN = _repo_root() / "data" / "lit" / "family_ingestion_plan.json"


def load_family_ingestion_plan(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    path = Path(file_path) if file_path is not None else DEFAULT_FAMILY_INGESTION_PLAN
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def build_family_ingestion_plan_artifact(file_path: Optional[Path | str] = None) -> Dict[str, Any]:
    payload = load_family_ingestion_plan(file_path)
    families = [dict(row) for row in payload.get("families", [])]

    payload_type_counts: Dict[str, int] = {}
    posture_counts: Dict[str, int] = {}
    wave_counts: Dict[str, int] = {}
    mapped_scope_count = 0
    unmapped_scope_families = []

    families.sort(key=lambda row: (int(row.get("implementation_wave", 99)), int(row.get("order_in_wave", 99)), str(row.get("slr_family", "99"))))

    for row in families:
        posture = str(row.get("strategic_posture", "unknown"))
        posture_counts[posture] = posture_counts.get(posture, 0) + 1

        wave = str(row.get("implementation_wave", "unknown"))
        wave_counts[wave] = wave_counts.get(wave, 0) + 1

        if str(row.get("scope_family_id", "")).strip():
            mapped_scope_count += 1
        else:
            unmapped_scope_families.append(str(row.get("slr_family", "unknown")))

        for payload_type in row.get("preferred_payload_types", []):
            key = str(payload_type)
            payload_type_counts[key] = payload_type_counts.get(key, 0) + 1

    first_wave = [
        str(row.get("slr_family", "unknown"))
        for row in families
        if int(row.get("implementation_wave", 99)) == 1
    ]
    active_sequence = [str(row.get("slr_family", "unknown")) for row in families]

    return {
        "source": str((Path(file_path) if file_path is not None else DEFAULT_FAMILY_INGESTION_PLAN).resolve().relative_to(_repo_root().resolve()).as_posix()),
        "template_source": str(payload.get("template_source", "")),
        "summary": {
            "family_count": len(families),
            "payload_type_counts": dict(sorted(payload_type_counts.items())),
            "strategic_posture_counts": dict(sorted(posture_counts.items())),
            "implementation_wave_counts": dict(sorted(wave_counts.items())),
            "recommended_first_wave": first_wave,
            "active_build_sequence": active_sequence,
            "mapped_scope_family_count": mapped_scope_count,
            "unmapped_scope_families": unmapped_scope_families,
            "state_of_art_policy": "extend_the_quantitative_core_by_explicit_family_lanes_with_machine_readable_payloads_not_narrative_only_docs",
            "identifier_policy": "scope_family_id_uses_the_same_canonical_chemistry_family_id_as_payloads_and_validation_outputs",
            "axis_policy": "chemistry_family_and_matrix_family_are_separate_axes_and_should_not_share_identifier_spaces",
        },
        "families": families,
    }


def render_family_ingestion_plan_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Family Ingestion Plan",
        "",
        "| SLR | Family | Posture | Runtime Concept | Payload Types | Wave | Source File | Next Build Action |",
        "| --- | --- | --- | --- | --- | ---: | --- | --- |",
    ]
    for row in payload.get("families", []):
        next_action = "; ".join(str(item) for item in row.get("next_curation_actions", [])[:1]) or "none"
        src = str(row.get("source_slr_file", "")) or "—"
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {row.get('family_id', 'unknown')} | {row.get('strategic_posture', 'unknown')} | {row.get('runtime_concept', 'unknown')} | "
            f"{', '.join(str(item) for item in row.get('preferred_payload_types', [])) or 'none'} | {row.get('implementation_wave', 'unknown')} | {src} | {next_action} |"
        )

    lines.extend([
        "",
        "## Target Modules",
        "",
        "| SLR | Target Modules | Target Observables or State Variables |",
        "| --- | --- | --- |",
    ])
    for row in payload.get("families", []):
        lines.append(
            f"| {row.get('slr_family', 'unknown')} | {'; '.join(str(item) for item in row.get('target_runtime_modules', [])) or 'none'} | "
            f"{'; '.join(str(item) for item in row.get('target_compounds_or_state_variables', [])) or 'none'} |"
        )

    summary = payload.get("summary", {})
    lines.extend(
        [
            "",
            f"Families tracked: {int(summary.get('family_count', 0))}",
            f"Recommended first wave: {', '.join(str(item) for item in summary.get('recommended_first_wave', [])) or 'none'}",
            f"Active build sequence: {', '.join(str(item) for item in summary.get('active_build_sequence', [])) or 'none'}",
            f"Mapped scope families: {int(summary.get('mapped_scope_family_count', 0))}",
            f"Unmapped scope families: {', '.join(str(item) for item in summary.get('unmapped_scope_families', [])) or 'none'}",
            f"Policy: {summary.get('state_of_art_policy', 'unknown')}",
            f"Identifier policy: {summary.get('identifier_policy', 'unknown')}",
            f"Axis policy: {summary.get('axis_policy', 'unknown')}",
        ]
    )
    return "\n".join(lines) + "\n"