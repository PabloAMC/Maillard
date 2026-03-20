import json
from dataclasses import dataclass
from typing import Dict, List, Optional
from pathlib import Path

EVIDENCE_STATES = {
    "externally_benchmarked",
    "internally_benchmarked",
    "transferred_prior",
    "still_missing",
}

TARGET_CLASSES = {
    "sulfur_positives",
    "strecker_aldehydes",
    "pyrazines",
    "furans_furanones",
    "adverse_lipid_markers",
}

@dataclass(frozen=True)
class TargetCompound:
    name: str
    target_class: str
    evidence_state: str
    
    def __post_init__(self):
        if self.evidence_state not in EVIDENCE_STATES:
            raise ValueError(f"Invalid evidence state: {self.evidence_state}")
        if self.target_class not in TARGET_CLASSES:
            raise ValueError(f"Invalid target class: {self.target_class}")


MATRIX_TARGET_PANEL: Dict[str, TargetCompound] = {
    # Sulfur positives
    "2methyl3furanthiol": TargetCompound("2methyl3furanthiol", "sulfur_positives", "externally_benchmarked"),
    "2furfurylthiol": TargetCompound("2furfurylthiol", "sulfur_positives", "externally_benchmarked"),
    "methional": TargetCompound("methional", "sulfur_positives", "externally_benchmarked"),
    "bis2methyl3furyldisulfide": TargetCompound("bis2methyl3furyldisulfide", "sulfur_positives", "transferred_prior"),
    
    # Strecker aldehydes
    "3methylbutanal": TargetCompound("3methylbutanal", "strecker_aldehydes", "externally_benchmarked"),
    "2methylbutanal": TargetCompound("2methylbutanal", "strecker_aldehydes", "externally_benchmarked"),
    
    # Pyrazines
    "25dimethylpyrazine": TargetCompound("25dimethylpyrazine", "pyrazines", "externally_benchmarked"),
    "23dimethylpyrazine": TargetCompound("23dimethylpyrazine", "pyrazines", "transferred_prior"),
    "2ethyl35dimethylpyrazine": TargetCompound("2ethyl35dimethylpyrazine", "pyrazines", "transferred_prior"),
    "tetramethylpyrazine": TargetCompound("tetramethylpyrazine", "pyrazines", "still_missing"),
    
    # Furans/Furanones
    "4hydroxy25dimethyl32hfuranone": TargetCompound("4hydroxy25dimethyl32hfuranone", "furans_furanones", "transferred_prior"),
    "2methyltetrahydrofuran3one": TargetCompound("2methyltetrahydrofuran3one", "furans_furanones", "transferred_prior"),
    
    # Adverse Lipid Markers (from Pratt-Singh family references and typical markers)
    "hexanal": TargetCompound("hexanal", "adverse_lipid_markers", "externally_benchmarked"),
    "nonanal": TargetCompound("nonanal", "adverse_lipid_markers", "externally_benchmarked"),
    "2pentylfuran": TargetCompound("2pentylfuran", "adverse_lipid_markers", "externally_benchmarked"),
    "1hexanol": TargetCompound("1hexanol", "adverse_lipid_markers", "externally_benchmarked"),
}

def export_matrix_target_panel(path: Path) -> None:
    data = {}
    for name, compound in MATRIX_TARGET_PANEL.items():
        data[name] = {
            "name": compound.name,
            "target_class": compound.target_class,
            "evidence_state": compound.evidence_state,
        }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)

def get_compound_evidence_state(name: str) -> str:
    from src.benchmark_validation import _normalize_name, BENCHMARK_NAME_ALIASES
    n = _normalize_name(name)
    
    if n in MATRIX_TARGET_PANEL:
        return MATRIX_TARGET_PANEL[n].evidence_state
        
    for alias, alternatives in BENCHMARK_NAME_ALIASES.items():
        if n in alternatives:
            if alias in MATRIX_TARGET_PANEL:
                return MATRIX_TARGET_PANEL[alias].evidence_state
                
    return "still_missing"
    
def get_compound_target_class(name: str) -> Optional[str]:
    from src.benchmark_validation import _normalize_name, BENCHMARK_NAME_ALIASES
    n = _normalize_name(name)
    
    if n in MATRIX_TARGET_PANEL:
        return MATRIX_TARGET_PANEL[n].target_class
        
    for alias, alternatives in BENCHMARK_NAME_ALIASES.items():
        if n in alternatives:
            if alias in MATRIX_TARGET_PANEL:
                return MATRIX_TARGET_PANEL[alias].target_class
                
    return None
