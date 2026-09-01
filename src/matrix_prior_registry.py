from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from src import data_paths
from src import data_access

_PRIORS_PATH = data_paths.COMPUTATIONAL_PRIORS


def _load_priors() -> dict[str, Any]:
    return data_access.load_json(_PRIORS_PATH)


COMPUTATIONAL_PRIORS_PAYLOAD = _load_priors()


def _index_by_protein(section_name: str) -> Dict[str, Dict[str, Any]]:
    indexed: Dict[str, Dict[str, Any]] = {}
    for entry in COMPUTATIONAL_PRIORS_PAYLOAD.get(section_name, []):
        protein_type = str(entry.get("protein_type", "")).strip()
        if protein_type:
            indexed[protein_type] = dict(entry)
    return indexed


_ACCESSIBILITY_WINDOWS = _index_by_protein("accessibility_windows")
_DENATURATION_HEURISTICS = _index_by_protein("denaturation_heuristics")
_VOLATILE_CLASS_PROFILES = _index_by_protein("volatile_class_profiles")
_MATRIX_CORRECTIONS = _index_by_protein("matrix_corrections")


def _normalize_string_list(values: Any) -> List[str]:
    if not isinstance(values, list):
        return []
    return [str(item).strip() for item in values if str(item).strip()]


def _row_process_state_matches(row: Dict[str, Any], process_state: Optional[str]) -> bool:
    normalized_state = str(process_state or "").strip()
    if not normalized_state:
        return True
    scope = _normalize_string_list(row.get("process_state_applicability", row.get("process_state_scope", [])))
    return not scope or normalized_state in scope


def get_accessibility_window_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _ACCESSIBILITY_WINDOWS.get(str(protein_type))


def get_denaturation_heuristic_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _DENATURATION_HEURISTICS.get(str(protein_type))


def get_volatile_class_profile_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _VOLATILE_CLASS_PROFILES.get(str(protein_type))


def get_matrix_correction_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _MATRIX_CORRECTIONS.get(str(protein_type))


def query_family_prior_entries(
    *,
    chemistry_family: Optional[str] = None,
    protein_type: Optional[str] = None,
    process_state: Optional[str] = None,
    payload_role: Optional[str] = None,
    observable_panel_tag: Optional[str] = None,
    supporting_family: Optional[str] = None,
) -> List[Dict[str, Any]]:
    normalized_family = str(chemistry_family or "").strip()
    normalized_protein = str(protein_type or "").strip()
    normalized_payload_role = str(payload_role or "").strip()
    normalized_panel_tag = str(observable_panel_tag or "").strip()
    normalized_supporting_family = str(supporting_family or "").strip()
    rows: List[Dict[str, Any]] = []
    section_family_metadata = COMPUTATIONAL_PRIORS_PAYLOAD.get("section_family_metadata", {})
    if not isinstance(section_family_metadata, dict):
        section_family_metadata = {}

    for section_name, entries in COMPUTATIONAL_PRIORS_PAYLOAD.items():
        if section_name == "section_family_metadata" or not isinstance(entries, list):
            continue
        defaults = section_family_metadata.get(section_name, {}) if isinstance(section_family_metadata.get(section_name, {}), dict) else {}
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            row = dict(defaults)
            row.update(dict(entry))
            if normalized_family and str(row.get("chemistry_family", "")).strip() != normalized_family:
                continue
            if normalized_protein and str(row.get("protein_type", "")).strip() not in {"", normalized_protein}:
                continue
            if normalized_payload_role and str(row.get("payload_role", "")).strip() != normalized_payload_role:
                continue
            if normalized_panel_tag and normalized_panel_tag not in _normalize_string_list(row.get("observable_panel_tags", [])):
                continue
            if normalized_supporting_family and normalized_supporting_family not in _normalize_string_list(row.get("supporting_families", [])):
                continue
            if not _row_process_state_matches(row, process_state):
                continue
            row["_section_name"] = section_name
            rows.append(row)
    return rows


def summarize_family_prior_bundle(
    *,
    protein_type: Optional[str] = None,
    process_state: Optional[str] = None,
    chemistry_family: Optional[str] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    bundle: Dict[str, List[Dict[str, Any]]] = {}
    for row in query_family_prior_entries(
        chemistry_family=chemistry_family,
        protein_type=protein_type,
        process_state=process_state,
    ):
        family_key = str(row.get("chemistry_family", "")).strip()
        if not family_key:
            continue
        process_state_scope = _normalize_string_list(row.get("process_state_scope", row.get("process_state_applicability", [])))
        summary_row = {
            "section_name": str(row.get("_section_name", "unknown")),
            "payload_role": str(row.get("payload_role", "unknown")),
            "source": str(row.get("source", "unknown")),
            "confidence_tier": str(row.get("confidence_tier", "unknown")),
            "observable_panel_tags": _normalize_string_list(row.get("observable_panel_tags", [])),
            "supporting_families": _normalize_string_list(row.get("supporting_families", [])),
            "process_state_scope": process_state_scope,
        }
        protein_value = str(row.get("protein_type", "")).strip()
        if protein_value:
            summary_row["protein_type"] = protein_value
        bundle.setdefault(family_key, []).append(summary_row)
    return bundle


def summarize_matrix_prior_bundle(protein_type: str) -> Dict[str, Dict[str, Any]]:
    def _summarize(entry: Optional[Dict[str, Any]], parameter: str) -> Optional[Dict[str, Any]]:
        if entry is None:
            return None
        process_state_applicability = entry.get("process_state_applicability", [])
        if not isinstance(process_state_applicability, list):
            process_state_applicability = []
        return {
            "parameter": parameter,
            "source": str(entry.get("source", "unknown")),
            "provenance_tier": str(entry.get("provenance_tier", "unknown")),
            "confidence_tier": str(entry.get("confidence_tier", "unknown")),
            "process_state_applicability": [str(item) for item in process_state_applicability],
            "notes": str(entry.get("notes", "")),
        }

    normalized = str(protein_type)
    bundle = {
        "accessibility_window": _summarize(get_accessibility_window_entry(normalized), "accessibility_window"),
        "denaturation_heuristic": _summarize(get_denaturation_heuristic_entry(normalized), "denaturation_heuristic"),
        "volatile_class_profile": _summarize(get_volatile_class_profile_entry(normalized), "volatile_class_profile"),
        "matrix_correction": _summarize(get_matrix_correction_entry(normalized), "matrix_correction"),
    }
    return {key: value for key, value in bundle.items() if value is not None}