import json
from dataclasses import dataclass
from typing import Dict, List, Optional

from src import data_paths

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
    "umami_support_markers",
    "pretreatment_state_markers",
    "safety_markers",
}

PANEL_ROLES = {
    "scored",
    "constrained",
    "diagnostic",
    "report_only",
}

OBSERVABLE_KINDS = {
    "volatile",
    "state_variable",
}

PANEL_PATH = data_paths.MATRIX_DECISION_PANEL

@dataclass(frozen=True)
class TargetCompound:
    name: str
    target_class: str
    evidence_state: str
    display_name: str
    chemistry_family: Optional[str] = None
    supporting_families: tuple[str, ...] = ()
    observable_panel_tags: tuple[str, ...] = ()
    panel_role: str = "diagnostic"
    observable_kind: str = "volatile"
    modeling_regimes: tuple[str, ...] = ()
    
    def __post_init__(self):
        if self.evidence_state not in EVIDENCE_STATES:
            raise ValueError(f"Invalid evidence state: {self.evidence_state}")
        if self.target_class not in TARGET_CLASSES:
            raise ValueError(f"Invalid target class: {self.target_class}")
        if self.panel_role not in PANEL_ROLES:
            raise ValueError(f"Invalid panel role: {self.panel_role}")
        if self.observable_kind not in OBSERVABLE_KINDS:
            raise ValueError(f"Invalid observable kind: {self.observable_kind}")


from src.text_utils import normalize_compound_name as _normalize_name


def _load_panel() -> tuple[Dict[str, TargetCompound], Dict[str, str]]:
    with open(PANEL_PATH, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    entries = payload.get("entries", {})
    target_class_metadata = payload.get("target_class_family_metadata", {})
    panel: Dict[str, TargetCompound] = {}
    alias_lookup: Dict[str, str] = {}
    for canonical_name, entry in entries.items():
        normalized_name = _normalize_name(canonical_name)
        defaults = target_class_metadata.get(str(entry.get("target_class", "")), {}) if isinstance(target_class_metadata, dict) else {}
        supporting_families = []
        for source in [defaults, entry]:
            values = source.get("supporting_families", []) if isinstance(source, dict) else []
            if isinstance(values, list):
                supporting_families.extend(str(item) for item in values)
        observable_panel_tags = []
        for source in [defaults, entry]:
            values = source.get("observable_panel_tags", []) if isinstance(source, dict) else []
            if isinstance(values, list):
                observable_panel_tags.extend(str(item) for item in values)
        modeling_regimes = []
        for source in [defaults, entry]:
            values = source.get("modeling_regimes", []) if isinstance(source, dict) else []
            if isinstance(values, list):
                modeling_regimes.extend(str(item) for item in values)
        compound = TargetCompound(
            name=normalized_name,
            target_class=str(entry["target_class"]),
            evidence_state=str(entry["evidence_state"]),
            display_name=str(entry.get("display_name", canonical_name)),
            chemistry_family=str(entry.get("chemistry_family", defaults.get("chemistry_family", ""))) or None,
            supporting_families=tuple(dict.fromkeys(supporting_families)),
            observable_panel_tags=tuple(dict.fromkeys(observable_panel_tags)),
            panel_role=str(entry.get("panel_role", defaults.get("panel_role", "diagnostic"))),
            observable_kind=str(entry.get("observable_kind", defaults.get("observable_kind", "volatile"))),
            modeling_regimes=tuple(dict.fromkeys(modeling_regimes)),
        )
        panel[normalized_name] = compound
        alias_lookup[normalized_name] = normalized_name
        for alias in entry.get("aliases", []):
            alias_lookup[_normalize_name(alias)] = normalized_name
        alias_lookup[_normalize_name(compound.display_name)] = normalized_name
    return panel, alias_lookup


MATRIX_TARGET_PANEL, MATRIX_TARGET_ALIASES = _load_panel()

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
        "chemistry_family": compound.chemistry_family,
        "supporting_families": list(compound.supporting_families),
        "observable_panel_tags": list(compound.observable_panel_tags),
        "panel_role": compound.panel_role,
        "observable_kind": compound.observable_kind,
        "modeling_regimes": list(compound.modeling_regimes),
        "decision_panel_source": data_paths.rel(PANEL_PATH),
    }


def iter_target_panel_entries(
    *,
    target_class: Optional[str] = None,
    chemistry_family: Optional[str] = None,
    panel_role: Optional[str] = None,
    observable_kind: Optional[str] = None,
) -> List[dict]:
    rows: List[dict] = []
    for compound in MATRIX_TARGET_PANEL.values():
        if target_class is not None and compound.target_class != target_class:
            continue
        if chemistry_family is not None and compound.chemistry_family != chemistry_family:
            continue
        if panel_role is not None and compound.panel_role != panel_role:
            continue
        if observable_kind is not None and compound.observable_kind != observable_kind:
            continue
        rows.append(
            {
                "canonical_name": compound.name,
                "display_name": compound.display_name,
                "target_class": compound.target_class,
                "evidence_state": compound.evidence_state,
                "chemistry_family": compound.chemistry_family,
                "supporting_families": list(compound.supporting_families),
                "observable_panel_tags": list(compound.observable_panel_tags),
                "panel_role": compound.panel_role,
                "observable_kind": compound.observable_kind,
                "modeling_regimes": list(compound.modeling_regimes),
            }
        )
    rows.sort(key=lambda row: str(row.get("display_name", row.get("canonical_name", "unknown"))))
    return rows
