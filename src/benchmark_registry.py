"""
src/benchmark_registry.py — Benchmark IO, metadata, and condition parsing.

Extracted from src/benchmark_validation.py (Priority 2 decomposition).
Contains all functions responsible for loading benchmark JSON files from disk,
inferring/resolving their metadata, and converting them to ReactionConditions
or Formulation dicts.  No evaluation or reporting logic lives here.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from src.conditions import ReactionConditions
from src.benchmark_types import BenchmarkMetadata
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT
from src.matrix_calibration_registry import determine_matrix_process_state
from src.precursor_resolver import resolve_many

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_DIR = ROOT / "data" / "benchmarks"
DEFAULT_TARGET_TAG = "meaty"

MATRIX_NAMES = (
    "pea protein isolate",
    "soy protein isolate",
    "brown rice protein isolate",
    "pea protein",
    "soy protein",
    "mycoprotein",
)

THERMODYNAMIC_GATING_POLICIES = {
    "off",
    "on",
    "auto",
    "diagnostic_only",
    "benchmark_facing",
    "not_applicable",
}

# ---------------------------------------------------------------------------
# File I/O helpers
# ---------------------------------------------------------------------------

def get_benchmark_files(benchmark_dir: Path = BENCHMARK_DIR) -> List[Path]:
    if not benchmark_dir.exists():
        return []
    return sorted(benchmark_dir.glob("*.json"))


def load_benchmark(bench_file: Path | str) -> dict:
    bench_path = Path(bench_file)
    with open(bench_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


# ---------------------------------------------------------------------------
# Metadata resolution
# ---------------------------------------------------------------------------

def _infer_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    benchmark_id = bench.get("benchmark_id", "unknown")
    protein_type = bench.get("protein_type", "free")
    if protein_type != "free":
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="matrix_headspace",
            execution_path="matrix_only",
            benchmark_engine="matrix_intake_headspace",
            comparator_signal="predicted_ppb",
            cantera_role="not_authoritative",
            target_snapshot_policy="excluded",
            thermodynamic_gating_policy="not_applicable",
            notes="Matrix-only benchmark requiring a dedicated precursor-accessibility path.",
        )
    if "cys_" in benchmark_id:
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="free_aa_sulfur",
            execution_path="free_precursor",
            benchmark_engine="fast_observable",
            comparator_signal="predicted_ppb",
            cantera_role="diagnostic_reference_only",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="diagnostic_only",
            notes="Free amino-acid sulfur benchmark.",
        )
    if "acrylamide" in benchmark_id or bench.get("benchmark_type") == "safety":
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="safety",
            execution_path="free_precursor",
            benchmark_engine="fast_observable",
            comparator_signal="predicted_ppb",
            cantera_role="diagnostic_reference_only",
            target_snapshot_policy="included",
            thermodynamic_gating_policy="diagnostic_only",
            notes="Safety-critical benchmark (e.g., Acrylamide).",
        )
    return BenchmarkMetadata(
        tier="SECONDARY",
        family="general",
        execution_path="free_precursor",
        benchmark_engine="fast_observable",
        comparator_signal="predicted_ppb",
        cantera_role="diagnostic_reference_only",
        target_snapshot_policy="included",
        thermodynamic_gating_policy="diagnostic_only",
        notes=None,
    )


def get_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    metadata = bench.get("metadata") or {}
    inferred = _infer_benchmark_metadata(bench)
    policy = DEFAULT_VALIDATION_CONTRACT.policy_for_execution_path(
        str(metadata.get("execution_path", inferred.execution_path))
    )
    return BenchmarkMetadata(
        tier=str(metadata.get("tier", inferred.tier)),
        family=str(metadata.get("family", inferred.family)),
        execution_path=policy.execution_path,
        benchmark_engine=str(metadata.get("benchmark_engine", policy.benchmark_engine)),
        comparator_signal=str(metadata.get("comparator_signal", policy.comparator_signal)),
        cantera_role=str(metadata.get("cantera_role", policy.cantera_role)),
        target_snapshot_policy=str(metadata.get("target_snapshot_policy", policy.target_snapshot_policy)),
        thermodynamic_gating_policy=str(
            metadata.get("thermodynamic_gating_policy", policy.thermodynamic_gating_policy)
        ),
        notes=metadata.get("notes", inferred.notes),
    )


def resolve_thermodynamic_gating_mode(bench: dict, requested_mode: str = "auto") -> str:
    metadata = get_benchmark_metadata(bench)
    normalized_mode = str(requested_mode).strip().lower()
    if normalized_mode not in THERMODYNAMIC_GATING_POLICIES:
        raise ValueError(f"Unsupported thermodynamic_gating mode: {requested_mode!r}")
    if normalized_mode == "auto":
        normalized_mode = metadata.thermodynamic_gating_policy
    if normalized_mode in {"off", "diagnostic_only", "not_applicable"}:
        return "off"
    if normalized_mode == "benchmark_facing":
        return "on"
    return normalized_mode


def is_matrix_only_benchmark(bench: dict | Path | str) -> bool:
    if isinstance(bench, (str, Path)):
        bench = load_benchmark(bench)
    return get_benchmark_metadata(bench).execution_path == "matrix_only"


# ---------------------------------------------------------------------------
# Contract helpers
# ---------------------------------------------------------------------------

def get_matrix_ranking_contract(bench: dict | Path | str) -> Dict[str, Any]:
    from src.matrix_targets import get_compound_evidence_state, get_compound_target_class

    if isinstance(bench, (Path, str)):
        bench = load_benchmark(bench)
    contract = bench.get("matrix_ranking_contract") or {}
    process_metadata = bench.get("process_metadata") or {}

    raw_observable_targets = [
        item
        for item in contract.get("observable_targets", [])
        if isinstance(item, dict) and item.get("name")
    ]
    observable_targets = []
    for item in raw_observable_targets:
        name = str(item["name"])
        enriched_item = dict(item)
        if "evidence_state" not in enriched_item:
            enriched_item["evidence_state"] = get_compound_evidence_state(name)
        if "target_class" not in enriched_item:
            enriched_item["target_class"] = get_compound_target_class(name) or "unknown"
        observable_targets.append(enriched_item)

    adverse_markers = [str(item) for item in contract.get("adverse_markers", [])]
    citation_provenance = [str(item) for item in contract.get("citation_provenance", [])]
    return {
        "process_state": process_metadata.get("state"),
        "process_metadata": process_metadata,
        "observable_targets": observable_targets,
        "adverse_markers": adverse_markers,
        "citation_provenance": citation_provenance,
        "calibration_mode": contract.get("calibration_mode"),
        "notes": contract.get("notes"),
    }


# ---------------------------------------------------------------------------
# Formulation / condition converters
# ---------------------------------------------------------------------------

def benchmark_to_conditions(bench: dict) -> ReactionConditions:
    conditions_payload = bench["conditions"]
    return ReactionConditions(
        pH=conditions_payload["ph"],
        temperature_celsius=conditions_payload["temp_C"],
        water_activity=conditions_payload["water_activity"],
        protein_type=bench.get("protein_type", "free"),
        temperature_profile=conditions_payload.get("temperature_profile"),
        water_activity_profile=conditions_payload.get("water_activity_profile"),
        pH_profile=conditions_payload.get("pH_profile"),
    )


def benchmark_to_formulation(bench: dict) -> dict:
    molar_ratios = {
        name: data["concentration_mM"]
        for name, data in bench["precursors"].items()
    }

    sugars: List[str] = []
    amino_acids: List[str] = []
    lipids: List[str] = []
    skipped_matrix_precursors: List[str] = []

    for name in bench["precursors"]:
        name_lower = name.lower()
        if any(matrix_name in name_lower for matrix_name in MATRIX_NAMES):
            skipped_matrix_precursors.append(name)
            continue
        if any(
            token in name_lower
            for token in ["ribose", "glucose", "fructose", "xylose", "maltose", "sugar"]
        ):
            sugars.append(name)
        elif any(
            token in name_lower
            for token in ["hexanal", "nonanal", "lipid", "fat", "furan"]
        ):
            lipids.append(name)
        else:
            amino_acids.append(name)

    return {
        "name": bench["benchmark_id"],
        "sugars": sugars,
        "amino_acids": amino_acids,
        "lipids": lipids,
        "molar_ratios": molar_ratios,
        "ph": bench["conditions"]["ph"],
        "temp": bench["conditions"]["temp_C"],
        "aw": bench["conditions"]["water_activity"],
        "time_minutes": bench["conditions"]["time_min"],
        "protein_type": bench.get("protein_type", "free"),
        "denaturation_state": bench.get("denaturation_state", 0.5),
        "_skipped_matrix_precursors": skipped_matrix_precursors,
    }


def _is_supported_formulation(formulation: dict) -> tuple[bool, Optional[str]]:
    protein_type = str(formulation.get("protein_type", "free"))
    if protein_type != "free" and formulation.get("_skipped_matrix_precursors"):
        from src.benchmark_validation import MATRIX_BENCHMARK_PROFILES
        if protein_type in MATRIX_BENCHMARK_PROFILES:
            return True, None
        return False, f"No executable matrix-only benchmark path for protein_type={protein_type}"

    candidate_precursors = (
        formulation["sugars"] + formulation["amino_acids"] + formulation["lipids"]
    )
    if not candidate_precursors:
        skipped = ", ".join(formulation.get("_skipped_matrix_precursors", [])) or "none"
        return False, f"No resolvable free-precursor system in benchmark. Matrix-only precursors: {skipped}"

    try:
        resolve_many(candidate_precursors)
    except ValueError as exc:
        return False, str(exc)
    return True, None


def get_matrix_only_target_snapshot_exclusions(
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> List[str]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    excluded: List[str] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        if is_matrix_only_benchmark(bench):
            excluded.append(str(bench.get("benchmark_id", Path(bench_file).stem)))
    return excluded
