from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional


ROOT = Path(__file__).resolve().parents[1]
INTAKE_REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

LEGACY_STATUS_TO_TRIAGE_STATUS = {
    "ready_for_reference_encoding": "ready_reference",
    "ready_for_calibration_encoding": "ready_calibration",
    "ready_for_runtime_encoding": "ready_calibration",
    "ready_for_intake_encoding": "ready_benchmark",
    "ready_for_directional_prior_encoding": "ready_directional_prior",
    "reviewed_qualitative_only": "reviewed_qualitative_only",
    "encoded_runtime_artifact": "encoded",
}

READY_TRIAGE_STATUSES = {
    "ready_reference",
    "ready_calibration",
    "ready_benchmark",
    "ready_directional_prior",
}

TRIAGE_STATUS_DEFAULT_TEMPLATE_KIND = {
    "ready_reference": "safety_payload",
    "ready_calibration": "process_state_calibration",
    "ready_benchmark": "benchmark_payload",
    "ready_directional_prior": "computational_prior",
}

ARTIFACT_TYPE_TO_TEMPLATE_KIND = {
    "benchmark": "benchmark_payload",
    "process_state_calibration": "process_state_calibration",
    "computational_prior": "computational_prior",
    "directional_prior": "computational_prior",
    "safety_reference": "safety_payload",
    "flavor_reference_payload": "flavor_reference_payload",
    "retention_payload": "retention_payload",
    "structural_gap_entry": "structural_gap_entry",
}


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_intake_registry(root: Path = ROOT) -> Dict[str, Any]:
    return _load_json(root / "data" / "lit" / "benchmark_intake_registry.json")


def _iter_payload_entries(payload: Any, *, section_name: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    if isinstance(payload, list):
        for entry in payload:
            if isinstance(entry, Mapping):
                row = dict(entry)
                if section_name is not None:
                    row.setdefault("_section_name", section_name)
                yield row
        return
    if isinstance(payload, Mapping):
        for key, value in payload.items():
            if isinstance(value, list):
                next_section = str(key) if section_name is None else section_name
                yield from _iter_payload_entries(value, section_name=next_section)
            elif isinstance(value, Mapping):
                next_section = str(key) if section_name is None else section_name
                yield from _iter_payload_entries(value, section_name=next_section)


def canonical_template_kind_for_artifact(artifact_type: str) -> str:
    normalized = str(artifact_type).strip()
    return ARTIFACT_TYPE_TO_TEMPLATE_KIND.get(normalized, normalized)


def _artifact_exists(artifact: Mapping[str, Any], *, root: Path = ROOT) -> bool:
    path = root / str(artifact.get("path", ""))
    if not path.exists():
        return False

    artifact_id = str(artifact.get("artifact_id", "")).strip()
    artifact_type = str(artifact.get("artifact_type", "")).strip()
    if not artifact_id or artifact_type == "benchmark":
        return True

    payload = _load_json(path)
    if artifact_type == "intake_registry_entry":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("eligible_references", []))
    if artifact_type == "slr_incorporation_ledger":
        return any(str(entry.get("paper_id", "")) == artifact_id for entry in payload.get("entries", []))
    if artifact_type == "process_state_calibration":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("entries", []))
    if artifact_type == "safety_reference":
        return any(str(entry.get("id", "")) == artifact_id for entry in payload.get("entries", []))

    sections = artifact.get("sections", [])
    if isinstance(sections, list) and sections:
        for section_name in sections:
            for entry in _iter_payload_entries(payload.get(str(section_name), []), section_name=str(section_name)):
                if str(entry.get("id", entry.get("protein_type", ""))) == artifact_id:
                    return True
        return False

    for entry in _iter_payload_entries(payload):
        if str(entry.get("id", entry.get("protein_type", ""))) == artifact_id:
            return True
    return False


def infer_target_payload_types(entry: Mapping[str, Any]) -> List[str]:
    payload_types: List[str] = []
    for artifact in entry.get("runtime_artifacts", []) or []:
        template_kind = canonical_template_kind_for_artifact(str(artifact.get("artifact_type", "")))
        if template_kind and template_kind not in payload_types:
            payload_types.append(template_kind)
    if payload_types:
        return payload_types

    triage_status = normalize_triage_status(entry)
    fallback = TRIAGE_STATUS_DEFAULT_TEMPLATE_KIND.get(triage_status, "reference_payload")
    return [fallback]


def resolve_primary_template_kind(entry: Mapping[str, Any]) -> str:
    payload_types = infer_target_payload_types(entry)
    return payload_types[0] if payload_types else "reference_payload"


def normalize_triage_status(entry: Mapping[str, Any]) -> str:
    explicit_status = str(entry.get("triage_status", "")).strip()
    if explicit_status:
        return explicit_status
    legacy_status = str(entry.get("status", "")).strip()
    return LEGACY_STATUS_TO_TRIAGE_STATUS.get(legacy_status, legacy_status or "unknown")


def _normalize_runtime_artifacts(entry: Mapping[str, Any], *, root: Path = ROOT) -> List[Dict[str, Any]]:
    runtime_artifacts: List[Dict[str, Any]] = []
    for artifact in entry.get("runtime_artifacts", []) or []:
        artifact_row = dict(artifact)
        artifact_row["template_kind"] = canonical_template_kind_for_artifact(str(artifact_row.get("artifact_type", "")))
        artifact_row["exists"] = _artifact_exists(artifact_row, root=root)
        runtime_artifacts.append(artifact_row)
    return runtime_artifacts


def _determine_encoding_status(entry: Mapping[str, Any], runtime_artifacts_present: bool) -> str:
    explicit_status = str(entry.get("encoding_status", "")).strip()
    if explicit_status:
        return explicit_status
    return "encoded_runtime_artifact" if runtime_artifacts_present else "template_required"


def _determine_backlog_queue(*, triage_status: str, encoding_status: str, template_kind: str) -> Optional[str]:
    if encoding_status == "encoded_runtime_artifact":
        return None
    if triage_status not in READY_TRIAGE_STATUSES:
        return None
    if template_kind == "benchmark_payload":
        return "ready_benchmark"
    return "ready_runtime"


def normalize_intake_entry(entry: Mapping[str, Any], *, root: Path = ROOT) -> Dict[str, Any]:
    runtime_artifacts = _normalize_runtime_artifacts(entry, root=root)
    runtime_artifacts_present = bool(runtime_artifacts) and all(bool(item.get("exists", False)) for item in runtime_artifacts)
    triage_status = normalize_triage_status(entry)
    template_kind = resolve_primary_template_kind(entry)
    encoding_status = _determine_encoding_status(entry, runtime_artifacts_present)
    backlog_queue = _determine_backlog_queue(
        triage_status=triage_status,
        encoding_status=encoding_status,
        template_kind=template_kind,
    )
    return {
        "id": str(entry.get("id", "unknown")),
        "citation": str(entry.get("citation", "unknown")),
        "kind": str(entry.get("kind", "unknown")),
        "chemistry_family": str(entry.get("chemistry_family", "")),
        "slr_family_source": str(entry.get("slr_family_source", "")),
        "source_payload_role": str(entry.get("payload_role", "unknown")),
        "target_payload_types": infer_target_payload_types(entry),
        "matrix_family": str(entry.get("matrix_family", "unknown")),
        "legacy_status": str(entry.get("status", "unknown")),
        "triage_status": triage_status,
        "encoding_status": encoding_status,
        "template_kind": template_kind,
        "runtime_artifacts_present": runtime_artifacts_present,
        "backlog_queue": backlog_queue,
        "requires_primary_data": bool(entry.get("requires_primary_data", False)),
        "target_modules": [str(item) for item in entry.get("target_modules", []) or []],
        "repo_next_action": str(entry.get("repo_next_action", "")),
        "runtime_artifacts": runtime_artifacts,
    }


def build_intake_reference_rows(root: Path = ROOT) -> List[Dict[str, Any]]:
    intake = load_intake_registry(root)
    rows: List[Dict[str, Any]] = []
    for entry in intake.get("eligible_references", []):
        normalized = normalize_intake_entry(entry, root=root)
        if normalized["triage_status"] in READY_TRIAGE_STATUSES or normalized["encoding_status"] == "encoded_runtime_artifact":
            rows.append(normalized)

    queue_order = {None: 0, "ready_runtime": 1, "ready_benchmark": 2}
    rows.sort(
        key=lambda row: (
            0 if row["encoding_status"] == "encoded_runtime_artifact" else 1,
            queue_order.get(row.get("backlog_queue"), 9),
            row["template_kind"],
            row["id"],
        )
    )
    return rows


def build_literature_backlog_artifact(root: Path = ROOT) -> Dict[str, Any]:
    intake = load_intake_registry(root)
    reference_rows = build_intake_reference_rows(root)
    ready_runtime = [row for row in reference_rows if row.get("backlog_queue") == "ready_runtime"]
    ready_benchmark = [row for row in reference_rows if row.get("backlog_queue") == "ready_benchmark"]
    encoded_rows = [row for row in reference_rows if row.get("encoding_status") == "encoded_runtime_artifact"]
    queue_conflicts = [
        row["id"]
        for row in reference_rows
        if row.get("backlog_queue") is not None and row.get("encoding_status") == "encoded_runtime_artifact"
    ]

    wet_lab_blocked: List[Dict[str, Any]] = []
    for gap in intake.get("structural_gaps", []):
        if str(gap.get("closure_outcome", "")) != "wet_lab_only":
            continue
        wet_lab_blocked.append(
            {
                "gap_id": str(gap.get("id", "unknown")),
                "priority": str(gap.get("priority", "unknown")),
                "triage_decision": str(gap.get("triage_decision", "")),
                "closure_outcome": str(gap.get("closure_outcome", "unknown")),
                "evidence_state": str(gap.get("evidence_state", "unknown")),
                "benchmark_contract_missing": [str(item) for item in gap.get("benchmark_contract_missing", []) or []],
                "near_miss_candidates": [
                    {
                        "entry_id": str(item.get("entry_id", "unknown")),
                        "reason": str(item.get("reason", "")),
                    }
                    for item in gap.get("near_miss_candidates", []) or []
                    if isinstance(item, Mapping)
                ],
                "why": str(gap.get("why", "")),
                "backlog_queue": "wet_lab_blocked",
            }
        )

    return {
        "schema_version": "1.0",
        "source": INTAKE_REGISTRY_PATH.relative_to(ROOT).as_posix(),
        "queue_policy": "ready queues are exclusive to non-encoded intake rows; wet_lab_blocked comes only from structural gaps with closure_outcome=wet_lab_only",
        "intake_reference_rows": reference_rows,
        "ready_runtime": ready_runtime,
        "ready_benchmark": ready_benchmark,
        "wet_lab_blocked": wet_lab_blocked,
        "encoded_reference_rows": encoded_rows,
        "minimum_primary_experiment": dict(intake.get("minimum_primary_experiment", {})),
        "conflicts": {
            "encoded_and_ready_ids": sorted(queue_conflicts),
        },
        "summary": {
            "intake_reference_count": len(reference_rows),
            "encoded_reference_count": len(encoded_rows),
            "ready_runtime_count": len(ready_runtime),
            "ready_benchmark_count": len(ready_benchmark),
            "wet_lab_blocked_count": len(wet_lab_blocked),
            "queue_conflict_count": len(queue_conflicts),
        },
    }


def render_literature_backlog_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# Literature Backlog",
        "",
        f"Queue policy: {payload.get('queue_policy', 'unknown')}",
        "",
        f"Encoded references: {payload.get('summary', {}).get('encoded_reference_count', 0)}",
        f"Ready runtime: {payload.get('summary', {}).get('ready_runtime_count', 0)}",
        f"Ready benchmark: {payload.get('summary', {}).get('ready_benchmark_count', 0)}",
        f"Wet-lab blocked: {payload.get('summary', {}).get('wet_lab_blocked_count', 0)}",
        f"Queue conflicts: {payload.get('summary', {}).get('queue_conflict_count', 0)}",
        "",
        "## Ready Runtime",
        "",
        "| ID | Template | Family | Matrix | Triage | Encoding | Next Action |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]

    ready_runtime = list(payload.get("ready_runtime", []))
    if not ready_runtime:
        lines.append("| none | n/a | n/a | n/a | n/a | n/a | no runtime-only backlog remains in the curated intake registry |")
    for row in ready_runtime:
        lines.append(
            f"| {row['id']} | {row['template_kind']} | {row['chemistry_family']} | {row['matrix_family']} | {row['triage_status']} | {row['encoding_status']} | {row.get('repo_next_action', 'none')} |"
        )

    lines.extend([
        "",
        "## Ready Benchmark",
        "",
        "| ID | Template | Family | Matrix | Triage | Encoding | Next Action |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ])

    ready_benchmark = list(payload.get("ready_benchmark", []))
    if not ready_benchmark:
        lines.append("| none | n/a | n/a | n/a | n/a | n/a | no benchmark-ready backlog remains in the curated intake registry |")
    for row in ready_benchmark:
        lines.append(
            f"| {row['id']} | {row['template_kind']} | {row['chemistry_family']} | {row['matrix_family']} | {row['triage_status']} | {row['encoding_status']} | {row.get('repo_next_action', 'none')} |"
        )

    lines.extend([
        "",
        "## Wet-Lab Blocked",
        "",
        "| Gap | Priority | Decision | Near Misses | Missing Contract |",
        "| --- | --- | --- | --- | --- |",
    ])
    for row in payload.get("wet_lab_blocked", []):
        near_misses = ", ".join(item.get("entry_id", "unknown") for item in row.get("near_miss_candidates", [])) or "none"
        missing_contract = ", ".join(row.get("benchmark_contract_missing", [])[:3]) or "none"
        lines.append(
            f"| {row['gap_id']} | {row['priority']} | {row.get('triage_decision', 'unknown')} | {near_misses} | {missing_contract} |"
        )

    minimum_primary_experiment = dict(payload.get("minimum_primary_experiment", {}))
    if minimum_primary_experiment:
        lines.extend([
            "",
            "## Minimum Primary Experiment",
            "",
            f"Matrices: {', '.join(minimum_primary_experiment.get('matrices', [])) or 'none'}",
            f"Exogenous precursors: {', '.join(f'{key}={value}' for key, value in minimum_primary_experiment.get('exogenous_precursors', {}).items()) or 'none'}",
            f"Analytical panel: {', '.join(minimum_primary_experiment.get('analytical_panel', [])) or 'none'}",
            f"Companion assays: {', '.join(minimum_primary_experiment.get('companion_assays', [])) or 'none'}",
            f"Instrumentation: {minimum_primary_experiment.get('instrumentation', 'unknown')}",
            "",
            "| Matrix | pH | Temp C | Time Points min |",
            "| --- | ---: | ---: | --- |",
        ])
        for condition in minimum_primary_experiment.get("conditions", []):
            lines.append(
                f"| {condition.get('matrix', 'unknown')} | {float(condition.get('ph', 0.0)):.1f} | {int(condition.get('temp_C', 0))} | {', '.join(str(item) for item in condition.get('time_points_min', []))} |"
            )

    return "\n".join(lines) + "\n"