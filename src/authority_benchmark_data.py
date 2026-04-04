from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from src.artifact_io import load_json_mapping


ROOT = Path(__file__).resolve().parents[1]
PHASE33_DATASET = ROOT / "data" / "qm" / "phase33_barrier_benchmarks.json"
PHASE35_DATASET = ROOT / "data" / "qm" / "phase35_double_hybrid_benchmarks.json"
IRC_FIXTURES_ROOT = ROOT / "data" / "qm" / "irc_validation_cases"


def _load_json_file(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Benchmark dataset not found: {path}")
    payload = load_json_mapping(path)
    if not isinstance(payload, dict):
        raise ValueError(f"Benchmark dataset must be a JSON object: {path}")
    return payload


def _require_numeric(entry: Mapping[str, Any], key: str, dataset_path: Path) -> float:
    value = entry.get(key)
    if not isinstance(value, (int, float)):
        raise ValueError(f"Expected numeric '{key}' in {dataset_path}")
    return float(value)


def load_phase33_barrier_benchmarks(dataset_path: Optional[Path] = None) -> Dict[str, Dict[str, Any]]:
    path = dataset_path or PHASE33_DATASET
    payload = _load_json_file(path)
    benchmarks = payload.get("benchmarks")
    if not isinstance(benchmarks, list) or not benchmarks:
        raise ValueError(f"Phase 3.3 benchmark dataset has no benchmarks: {path}")

    by_family: Dict[str, Dict[str, Any]] = {}
    for raw_entry in benchmarks:
        if not isinstance(raw_entry, dict):
            raise ValueError(f"Invalid benchmark entry in {path}")
        family = str(raw_entry.get("family", "")).strip()
        if not family:
            raise ValueError(f"Missing family in Phase 3.3 dataset: {path}")
        literature = raw_entry.get("literature")
        if not isinstance(literature, dict):
            raise ValueError(f"Missing literature block for {family} in {path}")
        for key in ("low", "high", "best"):
            _require_numeric(literature, key, path)
        _require_numeric(raw_entry, "xtb_reference_kcal_mol", path)
        _require_numeric(raw_entry, "wb97mv_kcal_mol", path)
        by_family[family] = dict(raw_entry)
    return by_family


def load_phase35_double_hybrid_benchmarks(dataset_path: Optional[Path] = None) -> Dict[str, Dict[str, Any]]:
    path = dataset_path or PHASE35_DATASET
    payload = _load_json_file(path)
    benchmarks = payload.get("benchmarks")
    if not isinstance(benchmarks, list) or not benchmarks:
        raise ValueError(f"Phase 3.5 benchmark dataset has no benchmarks: {path}")

    by_family: Dict[str, Dict[str, Any]] = {}
    for raw_entry in benchmarks:
        if not isinstance(raw_entry, dict):
            raise ValueError(f"Invalid double-hybrid benchmark entry in {path}")
        family = str(raw_entry.get("family", "")).strip()
        if not family:
            raise ValueError(f"Missing family in Phase 3.5 dataset: {path}")
        _require_numeric(raw_entry, "wb97mv_kcal_mol", path)
        _require_numeric(raw_entry, "revdsd_pbep86_d4_kcal_mol", path)
        by_family[family] = dict(raw_entry)
    return by_family


def _read_xyz(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing XYZ fixture: {path}")
    return path.read_text(encoding="utf-8").strip() + "\n"


def load_irc_validation_cases(fixtures_root: Optional[Path] = None) -> List[Dict[str, Any]]:
    root = fixtures_root or IRC_FIXTURES_ROOT
    if not root.exists():
        raise FileNotFoundError(f"IRC fixtures root not found: {root}")

    cases: List[Dict[str, Any]] = []
    for case_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        profile_path = case_dir / "profile.json"
        if not profile_path.exists():
            continue

        profile = _load_json_file(profile_path)
        case_id = str(profile.get("case_id", case_dir.name)).strip()
        family = str(profile.get("family", "")).strip()
        energies = profile.get("energies_hartree")
        if not case_id or not family:
            raise ValueError(f"IRC case metadata incomplete in {profile_path}")
        if not isinstance(energies, list) or len(energies) != 3:
            raise ValueError(f"IRC case must contain exactly 3 energy points: {profile_path}")
        energy_values = [float(value) for value in energies]

        cases.append(
            {
                "case_id": case_id,
                "family": family,
                "reaction_label": str(profile.get("reaction_label", case_id)),
                "reactant_smiles": str(profile.get("reactant_smiles", "")).strip(),
                "product_smiles": str(profile.get("product_smiles", "")).strip(),
                "expected_bond_change": str(profile.get("expected_bond_change", "")).strip(),
                "direct_barrier_kcal_mol": float(profile.get("direct_barrier_kcal_mol", 0.0)),
                "energies": energy_values,
                "ts_xyz": _read_xyz(case_dir / "ts.xyz"),
                "backward_endpoint_xyz": _read_xyz(case_dir / "backward.xyz"),
                "forward_endpoint_xyz": _read_xyz(case_dir / "forward.xyz"),
            }
        )

    if not cases:
        raise ValueError(f"No IRC validation cases found under {root}")
    return cases