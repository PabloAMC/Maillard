from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional


ROOT = Path(__file__).resolve().parents[1]
_PRIORS_PATH = ROOT / "data" / "lit" / "computational_priors.json"


def _load_priors() -> dict[str, Any]:
    with open(_PRIORS_PATH, "r", encoding="utf-8") as handle:
        return json.load(handle)


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


def get_accessibility_window_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _ACCESSIBILITY_WINDOWS.get(str(protein_type))


def get_denaturation_heuristic_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _DENATURATION_HEURISTICS.get(str(protein_type))


def get_volatile_class_profile_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _VOLATILE_CLASS_PROFILES.get(str(protein_type))


def get_matrix_correction_entry(protein_type: str) -> Optional[Dict[str, Any]]:
    return _MATRIX_CORRECTIONS.get(str(protein_type))


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
            "uncertainty_posture": str(entry.get("uncertainty_posture", "unknown")),
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