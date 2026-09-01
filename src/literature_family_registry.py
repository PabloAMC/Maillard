from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src import data_paths
from src import data_access
from src.family_ingestion_plan import load_family_ingestion_plan


# Module-level so tests can monkeypatch the intake registry location.
BENCHMARK_INTAKE_REGISTRY_PATH = data_paths.BENCHMARK_INTAKE_REGISTRY
PENDING_BENCHMARK_INTAKE_STATUSES = {"pending_json_payload"}

_CANONICAL_FAMILY_ALIASES = {
    "advanced_glycation_and_damage": "protein_damage_markers",
    "microbial_fermentation_modulation": "fermentation_pretreatment",
    "phospholipid_amine_maillard": "phospholipid_amine_sink",
    "lipid_oxidation_crosstalk": "lipid_oxidation_and_carbonylic_crosstalk",
    "carbohydrate_pyrolysis_caramelization": "carbohydrate_pyrolysis_and_caramelization",
}


def _load_json(path: Path) -> Dict[str, Any]:
    return data_access.load_json(path)


_FAMILY_PLAN = load_family_ingestion_plan()
_FAMILY_BY_SLR = {
    str(entry.get("slr_family", "")).zfill(2): dict(entry)
    for entry in _FAMILY_PLAN.get("families", [])
    if isinstance(entry, Mapping)
}
_FAMILY_BY_ID = {
    str(entry.get("family_id", "")).strip(): dict(entry)
    for entry in _FAMILY_PLAN.get("families", [])
    if isinstance(entry, Mapping)
}


def resolve_family_descriptor(family_ref: Optional[str]) -> Dict[str, Any]:
    normalized = str(family_ref or "").strip()
    normalized = _CANONICAL_FAMILY_ALIASES.get(normalized, normalized)
    if not normalized:
        return {}
    if normalized.isdigit():
        return dict(_FAMILY_BY_SLR.get(normalized.zfill(2), {}))
    return dict(_FAMILY_BY_ID.get(normalized, {}))


def _family_matches(row: Mapping[str, Any], family_ref: Optional[str]) -> bool:
    if family_ref is None:
        return True
    descriptor = resolve_family_descriptor(family_ref)
    target_family_id = str(descriptor.get("family_id", family_ref)).strip()
    target_slr = str(descriptor.get("slr_family", family_ref)).zfill(2) if str(family_ref).strip().isdigit() else str(descriptor.get("slr_family", "")).zfill(2)

    primary_family = str(row.get("chemistry_family", "")).strip()
    if primary_family and primary_family == target_family_id:
        return True

    primary_slr = str(row.get("slr_family_source", "")).zfill(2)
    if target_slr and primary_slr and primary_slr == target_slr:
        return True

    supporting = [str(item).strip() for item in row.get("supporting_families", []) or []]
    return target_family_id in supporting or target_slr in [item.zfill(2) for item in supporting if item.isdigit()]


def _apply_metadata_defaults(entry: Mapping[str, Any], defaults: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    row = dict(defaults or {})
    row.update(dict(entry))
    chemistry_family = str(row.get("chemistry_family", "")).strip()
    if chemistry_family:
        row["chemistry_family"] = _CANONICAL_FAMILY_ALIASES.get(chemistry_family, chemistry_family)
    supporting = []
    for source in [defaults or {}, entry]:
        values = source.get("supporting_families", []) if isinstance(source, Mapping) else []
        if isinstance(values, list):
            supporting.extend(
                _CANONICAL_FAMILY_ALIASES.get(str(item).strip(), str(item).strip())
                for item in values
                if str(item).strip()
            )
    if supporting:
        row["supporting_families"] = list(dict.fromkeys(supporting))
    observable_tags = []
    for source in [defaults or {}, entry]:
        values = source.get("observable_panel_tags", []) if isinstance(source, Mapping) else []
        if isinstance(values, list):
            observable_tags.extend(str(item).strip() for item in values if str(item).strip())
    if observable_tags:
        row["observable_panel_tags"] = list(dict.fromkeys(observable_tags))
    process_states = []
    for source in [defaults or {}, entry]:
        values = source.get("process_state_scope", []) if isinstance(source, Mapping) else []
        if isinstance(values, list):
            process_states.extend(str(item).strip() for item in values if str(item).strip())
    if process_states:
        row["process_state_scope"] = list(dict.fromkeys(process_states))
    family_descriptor = resolve_family_descriptor(str(row.get("chemistry_family", "")).strip())
    if not family_descriptor:
        slr_family_source = str(row.get("slr_family_source", row.get("slr_family", ""))).strip()
        if slr_family_source:
            family_descriptor = resolve_family_descriptor(slr_family_source)
    if family_descriptor:
        row.setdefault("chemistry_family", str(family_descriptor.get("family_id", "unknown")))
        row.setdefault("slr_family_source", str(family_descriptor.get("slr_family", "")).zfill(2))
        row["family_descriptor"] = {
            "slr_family": str(family_descriptor.get("slr_family", "")).zfill(2),
            "family_id": str(family_descriptor.get("family_id", "unknown")),
            "display_name": str(family_descriptor.get("display_name", "unknown")),
            "strategic_posture": str(family_descriptor.get("strategic_posture", "unknown")),
        }
    return row


def iter_benchmark_intake_entries(
    *,
    family: Optional[str] = None,
    include_pending: bool = False,
) -> Iterable[Dict[str, Any]]:
    payload = _load_json(BENCHMARK_INTAKE_REGISTRY_PATH)
    for entry in payload.get("eligible_references", []):
        if not isinstance(entry, Mapping):
            continue
        if not include_pending and str(entry.get("status", "")).strip() in PENDING_BENCHMARK_INTAKE_STATUSES:
            continue
        row = _apply_metadata_defaults(entry)
        if _family_matches(row, family):
            yield row


def iter_flavor_reference_entries(*, family: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    payload = _load_json(data_paths.FLAVOR_REFERENCE_PAYLOADS)
    metadata = payload.get("section_family_metadata", {}) if isinstance(payload.get("section_family_metadata", {}), Mapping) else {}
    for section_name, entries in payload.items():
        if section_name == "section_family_metadata" or not isinstance(entries, list):
            continue
        defaults = metadata.get(section_name, {}) if isinstance(metadata, Mapping) else {}
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            row = _apply_metadata_defaults(entry, defaults)
            row["_section_name"] = str(section_name)
            if _family_matches(row, family):
                yield row


def iter_retention_reference_entries(*, family: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    payload = _load_json(data_paths.RETENTION_REFERENCE_PAYLOADS)
    metadata = payload.get("section_family_metadata", {}) if isinstance(payload.get("section_family_metadata", {}), Mapping) else {}
    for section_name, section_payload in payload.items():
        if section_name == "section_family_metadata" or not isinstance(section_payload, Mapping):
            continue
        defaults = metadata.get(section_name, {}) if isinstance(metadata, Mapping) else {}
        for matrix_family, entries in section_payload.items():
            if not isinstance(entries, list):
                continue
            for entry in entries:
                if not isinstance(entry, Mapping):
                    continue
                row = _apply_metadata_defaults(entry, defaults)
                row["_section_name"] = str(section_name)
                row["_matrix_family"] = str(matrix_family)
                if _family_matches(row, family):
                    yield row


def iter_computational_prior_entries(
    *,
    family: Optional[str] = None,
    protein_type: Optional[str] = None,
) -> Iterable[Dict[str, Any]]:
    payload = _load_json(data_paths.COMPUTATIONAL_PRIORS)
    metadata = payload.get("section_family_metadata", {}) if isinstance(payload.get("section_family_metadata", {}), Mapping) else {}
    normalized_protein = str(protein_type).strip() if protein_type is not None else None
    for section_name, section_payload in payload.items():
        if section_name == "section_family_metadata":
            continue
        if isinstance(section_payload, list):
            entries = section_payload
        elif isinstance(section_payload, Mapping) and isinstance(section_payload.get("entries"), list):
            entries = section_payload.get("entries", [])
        else:
            continue
        defaults = metadata.get(section_name, {}) if isinstance(metadata, Mapping) else {}
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            row = _apply_metadata_defaults(entry, defaults)
            row["_section_name"] = str(section_name)
            if normalized_protein is not None and str(row.get("protein_type", "")).strip() not in {"", normalized_protein}:
                continue
            if _family_matches(row, family):
                yield row


def iter_process_gap_entries(*, family: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    payload = _load_json(data_paths.PROCESS_GAP_REGISTRY)
    for entry in payload.get("entries", []):
        if not isinstance(entry, Mapping):
            continue
        row = _apply_metadata_defaults(entry)
        if _family_matches(row, family):
            yield row


def iter_matrix_decision_panel_entries(*, family: Optional[str] = None) -> Iterable[Dict[str, Any]]:
    payload = _load_json(data_paths.MATRIX_DECISION_PANEL)
    metadata = payload.get("target_class_family_metadata", {}) if isinstance(payload.get("target_class_family_metadata", {}), Mapping) else {}
    entries = payload.get("entries", {})
    if not isinstance(entries, Mapping):
        return
    for canonical_name, entry in entries.items():
        if not isinstance(entry, Mapping):
            continue
        defaults = metadata.get(str(entry.get("target_class", "")), {}) if isinstance(metadata, Mapping) else {}
        row = _apply_metadata_defaults(entry, defaults)
        row["canonical_name"] = str(canonical_name)
        if _family_matches(row, family):
            yield row


def get_family_prior_entries(*, family: str, protein_type: Optional[str] = None) -> List[Dict[str, Any]]:
    return list(iter_computational_prior_entries(family=family, protein_type=protein_type))


def build_family_payload_coverage_artifact() -> Dict[str, Any]:
    payload_groups = {
        "benchmark_intake": list(iter_benchmark_intake_entries()),
        "flavor_reference_payload": list(iter_flavor_reference_entries()),
        "retention_payload": list(iter_retention_reference_entries()),
        "computational_prior": list(iter_computational_prior_entries()),
        "structural_gap_entry": list(iter_process_gap_entries()),
        "decision_panel_entry": list(iter_matrix_decision_panel_entries()),
    }

    family_rows: List[Dict[str, Any]] = []
    supported_families = 0
    for family in sorted(_FAMILY_PLAN.get("families", []), key=lambda row: (int(row.get("implementation_wave", 99)), int(row.get("order_in_wave", 99)), str(row.get("slr_family", "99")))):
        family_id = str(family.get("family_id", "unknown"))
        slr_family = str(family.get("slr_family", "")).zfill(2)
        primary_counts: Dict[str, int] = {}
        supporting_counts: Dict[str, int] = {}
        for role, rows in payload_groups.items():
            primary_counts[role] = sum(1 for row in rows if str(row.get("chemistry_family", "")).strip() == family_id or str(row.get("slr_family_source", "")).zfill(2) == slr_family)
            supporting_counts[role] = sum(1 for row in rows if family_id in [str(item).strip() for item in row.get("supporting_families", []) or []] or slr_family in [str(item).zfill(2) for item in row.get("supporting_families", []) or [] if str(item).isdigit()])
        total_primary = sum(primary_counts.values())
        total_supporting = sum(supporting_counts.values())
        if total_primary > 0 or total_supporting > 0:
            supported_families += 1
        family_rows.append(
            {
                "slr_family": slr_family,
                "family_id": family_id,
                "display_name": str(family.get("display_name", "unknown")),
                "strategic_posture": str(family.get("strategic_posture", "unknown")),
                "implementation_wave": int(family.get("implementation_wave", 99)),
                "primary_payload_counts": primary_counts,
                "supporting_payload_counts": supporting_counts,
                "total_primary_payload_count": int(total_primary),
                "total_supporting_payload_count": int(total_supporting),
                "has_runtime_support": bool(total_primary > 0 or total_supporting > 0),
            }
        )

    return {
        "source": "data/lit/* family-aware payload metadata",
        "summary": {
            "family_count": len(family_rows),
            "families_with_primary_payload_support": supported_families,
            "payload_roles": list(payload_groups.keys()),
        },
        "families": family_rows,
    }
