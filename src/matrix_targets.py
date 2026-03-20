import json
import re
from dataclasses import dataclass
from typing import Dict, Optional
from pathlib import Path

EVIDENCE_STATES = {
    "externally_benchmarked",
    "internally_benchmarked",
    "conditional_calibration",
    "safety_reference",
    "transferred_prior",
    "still_missing",
}

TARGET_CLASSES = {
    "sulfur_positives",
    "strecker_aldehydes",
    "pyrazines",
    "furans_furanones",
    "adverse_lipid_markers",
    "severity_markers",
    "safety_markers",
}

ROOT = Path(__file__).resolve().parents[1]
PANEL_PATH = ROOT / "data" / "lit" / "matrix_decision_panel.json"

@dataclass(frozen=True)
class TargetCompound:
    name: str
    target_class: str
    evidence_state: str
    display_name: str
    
    def __post_init__(self):
        if self.evidence_state not in EVIDENCE_STATES:
            raise ValueError(f"Invalid evidence state: {self.evidence_state}")
        if self.target_class not in TARGET_CLASSES:
            raise ValueError(f"Invalid target class: {self.target_class}")


def _normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _load_panel() -> tuple[Dict[str, TargetCompound], Dict[str, str]]:
    with open(PANEL_PATH, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    entries = payload.get("entries", {})
    panel: Dict[str, TargetCompound] = {}
    alias_lookup: Dict[str, str] = {}
    for canonical_name, entry in entries.items():
        normalized_name = _normalize_name(canonical_name)
        compound = TargetCompound(
            name=normalized_name,
            target_class=str(entry["target_class"]),
            evidence_state=str(entry["evidence_state"]),
            display_name=str(entry.get("display_name", canonical_name)),
        )
        panel[normalized_name] = compound
        alias_lookup[normalized_name] = normalized_name
        for alias in entry.get("aliases", []):
            alias_lookup[_normalize_name(alias)] = normalized_name
        alias_lookup[_normalize_name(compound.display_name)] = normalized_name
    return panel, alias_lookup


MATRIX_TARGET_PANEL, MATRIX_TARGET_ALIASES = _load_panel()

def export_matrix_target_panel(path: Path) -> None:
    data = {}
    for name, compound in MATRIX_TARGET_PANEL.items():
        data[name] = {
            "name": compound.name,
            "display_name": compound.display_name,
            "target_class": compound.target_class,
            "evidence_state": compound.evidence_state,
        }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)

def get_compound_evidence_state(name: str) -> str:
    n = _normalize_name(name)
    canonical = MATRIX_TARGET_ALIASES.get(n)
    if canonical and canonical in MATRIX_TARGET_PANEL:
        return MATRIX_TARGET_PANEL[canonical].evidence_state
    return "still_missing"
    
def get_compound_target_class(name: str) -> Optional[str]:
    n = _normalize_name(name)
    canonical = MATRIX_TARGET_ALIASES.get(n)
    if canonical and canonical in MATRIX_TARGET_PANEL:
        return MATRIX_TARGET_PANEL[canonical].target_class
    return None


def get_compound_panel_entry(name: str) -> Optional[dict]:
    n = _normalize_name(name)
    canonical = MATRIX_TARGET_ALIASES.get(n)
    if not canonical:
        return None
    compound = MATRIX_TARGET_PANEL.get(canonical)
    if compound is None:
        return None
    return {
        "canonical_name": canonical,
        "display_name": compound.display_name,
        "target_class": compound.target_class,
        "evidence_state": compound.evidence_state,
        "decision_panel_source": str(PANEL_PATH.relative_to(ROOT)),
    }
