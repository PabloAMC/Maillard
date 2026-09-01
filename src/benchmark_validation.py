from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

from src import data_paths
from src import data_access
from src.conditions import ReactionConditions
from src.barrier_constants import effective_barrier_from_rate_constant
from src.headspace import HeadspaceModel
from src.kinetics import KineticsEngine
from src.lipid_oxidation import (
    HYDROPEROXIDE_POOL_KEYS,
    MARKER_HYDROPEROXIDE_POOL,
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    hydroperoxide_pool_key_for_marker,
    predict_hexanal_generation,
)
from src.fit_target_index import fit_target_records_for, is_per_row_fit_target
from src.matrix_calibration_registry import (
    describe_matrix_calibration,
    determine_matrix_process_state,
    is_fit_recovery_benchmark,
)
from src.pipeline import MaillardPipeline
from src.precursor_resolver import resolve_many
from src.protein_binding import (
    binding_mode_active,
    observability_factor as binding_observability_factor,
    observability_mode,
    resolve_binding_context,
)
from src.projection_metadata import make_projection_metadata_row
from src.smirks_engine import SmirksEngine
from src.validation_contract import BenchmarkThresholds, DEFAULT_VALIDATION_CONTRACT
from src.safety import predict_acrylamide, predict_cel, predict_cml, predict_furosine
from src.matrix_targets import get_compound_panel_entry
from src.literature_family_registry import iter_benchmark_intake_entries, resolve_family_descriptor
from src.literature_runtime import build_family_upstream_contract
from src.benchmark_types import (
    BenchmarkMetadata,
    BenchmarkIndexEntry,
    BenchmarkNotSupportedError,
    CompoundComparison,
    BenchmarkEvaluation,
    BenchmarkSummary,
    MatrixBenchmarkDelta,
    MatrixBenchmarkEvidence,
    MatrixBenchmarkBranchDelta,
    ThermodynamicGatingAudit,
    BenchmarkTargetSnapshot,
    MatrixBenchmarkAssertion,
    MatrixPromotionFamilyStatus,
)


BENCHMARK_DIR = data_paths.BENCHMARKS_DIR
DEFAULT_TARGET_TAG = "meaty"
MATRIX_NAMES = (
    "pea protein isolate",
    "soy protein isolate",
    "brown rice protein isolate",
    "pea protein",
    "soy protein",
    "wheat gluten",
    "mycoprotein",
)

SUPPORT_ADDITIVE_TOKENS = (
    "thiamine",
    "vitamin b1",
    "imp",
    "gmp",
    "inosinate",
    "guanylate",
    "yeast extract",
    "ascorbic acid",
    "ascorbate",
    "vitamin c",
    "lecithin",
    "phospholipid",
    "phosphatidyl",
)

# The hand-written BENCHMARK_NAME_ALIASES table that lived here until 2026-09-01 is now
# derived from data/keys/compounds.yml (src.compound_keys.match_norms); the literal is
# frozen in tests/unit/test_compound_keys.py as the floor the registry must cover.

MATRIX_BENCHMARK_PROFILES = {
    "pea_iso": {
        "lipid_profile": PEA_LIPID_PROFILE,
    },
    "soy_iso": {
        "lipid_profile": SOY_LIPID_PROFILE,
    },
}

MATRIX_BENCHMARK_BASE_MARKER_YIELDS = {
    # Pea-referenced observable yields from the Pratap-Singh 2021 ambient-slurry
    # family. Soy-vs-pea release differences now live explicitly in HeadspaceModel
    # so the matrix-only path separates oxidation intake from headspace observability.
    #
    # WHAT THESE ARE, DIMENSIONALLY (2026-08-28, Wave Y). The matrix_only lane computes
    #
    #     observable_ppb = L(pool, lane, T, t) * Y(compound) * release * cal
    #
    # so `Y` converts a pseudo-hydroperoxide load into marker ppb. `Y` is DEGENERATE with
    # `kinetics.hydroperoxide_scale` (1.0e6, a round unanchored number in
    # data/lit/lipid_oxidation_calibration.json): only the product is identified. So a
    # `Y` above 1 is NOT a unit violation -- which is exactly the asymmetry that makes the
    # relocation below correct. An OBSERVABILITY factor is a fraction of a total and
    # cannot exceed 1 (Wave S4 (b)); a marker yield against an arbitrary scale can.
    #
    # 2026-08-28 (Wave Y) -- Hexanal 0.205 -> 0.885036. THE SCALE MOVED SIDES; NO NEW
    # PARAMETER WAS INTRODUCED. Wave O fitted ONE shared scale s = 4.317249 against the
    # two CONTENT-VERIFIED Pratap-Singh anchors (pea 1138.00, soy 1621.71 ppb, Molecules
    # 2021, 26, 4104 Table 1 via PMC8271896) and wrote it into the ambient-slurry HEXANAL
    # observability factors, which took the pea reference lane's factor from its
    # by-construction 1.0 to 4.31725 -- "a fraction of the total that the measurement
    # sees", equal to 4.32. Wave S4 (b) named that as the tell that the deficit lived
    # here instead. This wave moves the same single constant to this side of the product:
    #
    #     0.205 * 4.317249 = 0.885036       and every hexanal `cal` is divided by 4.317249
    #     (src/matrix_calibration_registry.py: 4.31725 -> 1.0, 0.228776 -> 0.0529912,
    #      9.54007 -> 0.453/0.205, 2.80478 -> (0.453/0.205)*(1-0.7060))
    #
    # The CALIBRATED tier is preserved to 6 significant figures on every lane (measured;
    # every hexanal lane in the repository resolves to a compound_specific record, so
    # none reaches the class-level aldehyde anchor, which is deliberately not moved). The
    # UNCALIBRATED tier -- `_uncalibrated_prediction_ppb`, which reads this dict and never
    # reads an observability factor -- moves by the full 4.317249x, and that movement is
    # the point: it is the first time the matrix sigma derivation has been able to see
    # the correction at all.
    #
    # Derivation, corpus and identifiability:
    #   scripts/generators/rederive_matrix_marker_yields.py
    #   results/validation/matrix_marker_yield_rederivation.{json,md}
    #
    # previous_value: 0.205 (Wave O era and earlier; itself 260/1268.3, i.e. back-solved
    # from a hexanal figure Wave K proved is not in the paper).
    "Hexanal": 0.885036,
    # CONFIRMED, NOT MOVED (Wave Y). 638 ppb was verified VERBATIM against the paper, and
    # the pea-reference requirement re-derives to 0.5017897 -- 0.000182 dex from the
    # shipped value, below the 0.01 dex materiality floor Wave O used. Moving it would be
    # fitting rounding.
    "2-Pentylfuran": 0.502,
    # UNANCHORED (Wave Y, restating Wave T3). Pratap-Singh report n.d. for hexanol in BOTH
    # matrices; 0.063 is 80/1269.8 and Wave T3 traced the 80 to the repository's own
    # abstract-reconstructed brief. There is nothing to fit against, so it was NOT fitted.
    # The 1-hexanol lane is unanchored end to end -- yield AND factor. [P]
    "1-Hexanol": 0.063,
    # UNANCHORED (Wave Y). Pratap-Singh report no nonanal on either ambient lane, so the
    # only nonanal target in the repository is Trikusuma 2019 (content-unverified) and its
    # observability factor is back-solved from that very row. Nothing here is fitted.
    "Nonanal": 0.150,
}







THERMODYNAMIC_GATING_POLICIES = {
    "off",
    "on",
    "auto",
    "diagnostic_only",
    "benchmark_facing",
    "not_applicable",
}


class _CompatibilityStatus(str):
    def __new__(cls, value: str, *, aliases: tuple[str, ...] = ()):
        obj = super().__new__(cls, value)
        obj._aliases = tuple(str(alias) for alias in aliases)
        return obj

    def __eq__(self, other: object) -> bool:
        return bool(super().__eq__(other) or other in self._aliases)


def _build_runtime_benchmark_family_map() -> tuple[Dict[str, str], Dict[str, List[str]]]:
    benchmark_family_map: Dict[str, str] = {}
    payload_roles_by_benchmark: Dict[str, set[str]] = defaultdict(set)
    for entry in iter_benchmark_intake_entries():
        chemistry_family = str(entry.get("chemistry_family", "")).strip()
        if not chemistry_family:
            continue
        payload_role = str(entry.get("payload_role", "benchmark_intake")).strip() or "benchmark_intake"
        for artifact in entry.get("runtime_artifacts", []) or []:
            if str(artifact.get("artifact_type", "")).strip() != "benchmark":
                continue
            artifact_id = str(artifact.get("artifact_id", "")).strip()
            if not artifact_id:
                continue
            benchmark_family_map[artifact_id] = chemistry_family
            payload_roles_by_benchmark[artifact_id].add(payload_role)
    return benchmark_family_map, {key: sorted(value) for key, value in payload_roles_by_benchmark.items()}


_RUNTIME_BENCHMARK_FAMILY_MAP, _RUNTIME_BENCHMARK_PAYLOAD_ROLES = _build_runtime_benchmark_family_map()


DEFAULT_BENCHMARK_THRESHOLDS = DEFAULT_VALIDATION_CONTRACT.thresholds
THERMODYNAMIC_GATING_THRESHOLD_KCAL = 30.0
THERMODYNAMIC_GATING_MIN_ABSOLUTE_MAE_IMPROVEMENT_PPB = 5.0
THERMODYNAMIC_GATING_MIN_RELATIVE_MAE_IMPROVEMENT = 0.05
THERMODYNAMIC_GATING_MIN_RATIO_IMPROVEMENT = 0.05


def get_matrix_ranking_contract(bench: dict | Path | str) -> Dict[str, Any]:
    from src.matrix_targets import get_compound_evidence_state, get_compound_target_class

    if isinstance(bench, (Path, str)):
        bench = load_benchmark(bench)
    contract = bench.get("matrix_ranking_contract") or {}
    process_metadata = bench.get("process_metadata") or {}
    
    raw_observable_targets = [
        item for item in contract.get("observable_targets", []) if isinstance(item, dict) and item.get("name")
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


def _evaluate_matrix_ranking_contract(
    bench: dict,
    predicted_ppb: Dict[str, float],
) -> Dict[str, object]:
    contract = get_matrix_ranking_contract(bench)
    observable_targets = contract.get("observable_targets", [])
    if not observable_targets:
        return {
            "status": "missing_contract",
            "ranked_observable_targets": [],
            "adverse_markers": contract.get("adverse_markers", []),
        }

    expected = sorted(
        observable_targets,
        key=lambda item: int(item.get("expected_rank", 999)),
    )
    predicted_rows = []
    missing_targets = []
    for item in expected:
        matched_name, predicted_value, _score = _best_prediction_match(str(item.get("name", "")), predicted_ppb)
        if matched_name is None:
            missing_targets.append(str(item.get("name", "")))
            continue
        predicted_rows.append((str(item.get("name", "")), float(predicted_value)))

    if missing_targets:
        status = "missing_targets"
    else:
        predicted_order = [name for name, _value in sorted(predicted_rows, key=lambda row: row[1], reverse=True)]
        expected_order = [str(item.get("name", "")) for item in expected]
        status = "pass" if predicted_order == expected_order else "order_mismatch"

    return {
        "status": status,
        "ranked_observable_targets": [str(item.get("name", "")) for item in expected],
        "adverse_markers": contract.get("adverse_markers", []),
    }


def is_matrix_only_benchmark(bench: dict | Path | str) -> bool:
    if isinstance(bench, (str, Path)):
        bench = load_benchmark(bench)
    return get_benchmark_metadata(bench).execution_path == "matrix_only"


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
    policy = DEFAULT_VALIDATION_CONTRACT.policy_for_execution_path(str(metadata.get("execution_path", inferred.execution_path)))
    return BenchmarkMetadata(
        tier=str(metadata.get("tier", inferred.tier)),
        family=str(metadata.get("family", inferred.family)),
        execution_path=policy.execution_path,
        benchmark_engine=str(metadata.get("benchmark_engine", policy.benchmark_engine)),
        comparator_signal=str(metadata.get("comparator_signal", policy.comparator_signal)),
        cantera_role=str(metadata.get("cantera_role", policy.cantera_role)),
        target_snapshot_policy=str(metadata.get("target_snapshot_policy", policy.target_snapshot_policy)),
        thermodynamic_gating_policy=str(metadata.get("thermodynamic_gating_policy", policy.thermodynamic_gating_policy)),
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


def get_benchmark_files(benchmark_dir: Path = BENCHMARK_DIR) -> List[Path]:
    if not benchmark_dir.exists():
        return []
    return sorted(benchmark_dir.glob("*.json"))


def load_benchmark(bench_file: Path | str) -> dict:
    # Cached and strict (2026-09-01): this is called at ~15 sites per run.
    return data_access.load_json(Path(bench_file))


def _get_condition_water_activity(conditions: Dict[str, Any], *, required: bool = False) -> Optional[float]:
    aw = conditions.get("water_activity")
    if aw is None:
        aw = conditions.get("aw")
    if aw is None:
        if required:
            raise KeyError("Benchmark conditions must include 'water_activity' or 'aw'.")
        return None
    return float(aw)


def benchmark_to_conditions(bench: dict) -> ReactionConditions:
    conditions = bench["conditions"]
    return ReactionConditions(
        pH=conditions["ph"],
        temperature_celsius=conditions["temp_C"],
        water_activity=_get_condition_water_activity(conditions, required=True),
        protein_type=bench.get("protein_type", "free"),
        sme_kj_per_kg=bench.get("sme_kj_per_kg", conditions.get("sme_kj_per_kg", 0.0)),
        moisture_regime=bench.get("moisture_regime", conditions.get("moisture_regime")),
        sterilization_temperature_celsius=bench.get(
            "sterilization_temperature_celsius",
            conditions.get("sterilization_temperature_celsius"),
        ),
        sterilization_time_minutes=bench.get(
            "sterilization_time_minutes",
            conditions.get("sterilization_time_minutes", 0.0),
        ),
    )


def benchmark_to_formulation(bench: dict) -> dict:
    conditions = bench["conditions"]
    molar_ratios = {
        name: data["concentration_mM"]
        for name, data in bench["precursors"].items()
    }

    sugars: List[str] = []
    amino_acids: List[str] = []
    additives: List[str] = []
    lipids: List[str] = []
    skipped_matrix_precursors: List[str] = []

    for name in bench["precursors"]:
        name_lower = name.lower()
        if any(matrix_name in name_lower for matrix_name in MATRIX_NAMES):
            skipped_matrix_precursors.append(name)
            continue
        if any(token in name_lower for token in ["ribose", "glucose", "fructose", "xylose", "maltose", "sugar"]):
            sugars.append(name)
        elif any(token in name_lower for token in ["hexanal", "nonanal", "lipid", "fat", "furan"]):
            lipids.append(name)
        elif any(token in name_lower for token in SUPPORT_ADDITIVE_TOKENS):
            additives.append(name)
        else:
            amino_acids.append(name)

    formulation = {
        "name": bench["benchmark_id"],
        "sugars": sugars,
        "amino_acids": amino_acids,
        "additives": additives,
        "lipids": lipids,
        "molar_ratios": molar_ratios,
        "ph": conditions["ph"],
        "temp": conditions["temp_C"],
        "aw": _get_condition_water_activity(conditions, required=True),
        "time_minutes": conditions["time_min"],
        "protein_type": bench.get("protein_type", "free"),
        "protein_source": bench.get("protein_source"),
        "support_cues": skipped_matrix_precursors,
        "denaturation_state": bench.get("denaturation_state", 0.5),
        "_skipped_matrix_precursors": skipped_matrix_precursors,
    }
    return formulation


from src.text_utils import normalize_compound_name as _normalize_name
from src import compound_keys


def _tokenize_name(name: str) -> List[str]:
    return [token for token in re.split(r"[^a-z0-9]+", name.lower()) if token]


def _is_supported_formulation(formulation: dict) -> tuple[bool, Optional[str]]:
    protein_type = str(formulation.get("protein_type", "free"))
    if protein_type != "free" and formulation.get("_skipped_matrix_precursors"):
        if protein_type in MATRIX_BENCHMARK_PROFILES:
            return True, None
        return False, f"No executable matrix-only benchmark path for protein_type={protein_type}"

    candidate_precursors = formulation["sugars"] + formulation["amino_acids"] + formulation.get("additives", []) + formulation["lipids"]
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


def _best_prediction_match(target_name: str, predicted_ppb: Dict[str, float]) -> tuple[Optional[str], float, float]:
    target_norm = _normalize_name(target_name)
    target_tokens = {token for token in _tokenize_name(target_name) if len(token) >= 4}
    target_aliases = compound_keys.match_norms(target_name)
    best_name: Optional[str] = None
    best_score = -1.0

    for candidate_name in predicted_ppb:
        candidate_norm = _normalize_name(candidate_name)
        if not candidate_norm:
            continue
        if len(candidate_norm) < 4:
            continue

        candidate_tokens = {token for token in _tokenize_name(candidate_name) if len(token) >= 4}

        if candidate_norm == target_norm or candidate_norm in target_aliases:
            score = 1.0
        elif target_tokens and candidate_tokens and target_tokens.intersection(candidate_tokens):
            overlap = len(target_tokens.intersection(candidate_tokens))
            score = overlap / max(len(target_tokens), len(candidate_tokens))
        else:
            score = SequenceMatcher(None, target_norm, candidate_norm).ratio()

        if score > best_score:
            best_score = score
            best_name = candidate_name

    if best_score < 0.75:
        return None, 0.0, 0.0
    return best_name, predicted_ppb[best_name], best_score


def _pearson(values_a: Iterable[float], values_b: Iterable[float]) -> Optional[float]:
    data_a = list(values_a)
    data_b = list(values_b)
    if len(data_a) < 3 or len(data_a) != len(data_b):
        return None

    mean_a = sum(data_a) / len(data_a)
    mean_b = sum(data_b) / len(data_b)
    centered_a = [value - mean_a for value in data_a]
    centered_b = [value - mean_b for value in data_b]
    numerator = sum(a * b for a, b in zip(centered_a, centered_b))
    denom_a = math.sqrt(sum(a * a for a in centered_a))
    denom_b = math.sqrt(sum(b * b for b in centered_b))
    if denom_a == 0.0 or denom_b == 0.0:
        return None
    return numerator / (denom_a * denom_b)


def _mean_abs_log10_error(comparisons: Iterable[CompoundComparison]) -> Optional[float]:
    matched = [comparison for comparison in comparisons if comparison.matched_name is not None]
    errors = [
        abs(math.log10(comparison.predicted_ppb / comparison.measured_ppb))
        for comparison in matched
        if comparison.measured_ppb > 0.0 and comparison.predicted_ppb > 0.0
    ]
    if not errors:
        return None
    return sum(errors) / len(errors)


def _resolve_scale_thresholds(
    bench: dict,
    *,
    protein_type: str,
    thresholds: BenchmarkThresholds,
) -> Dict[str, float]:
    configured = (bench.get("validation_contract") or {}).get("scale_thresholds") or {}
    return {
        "max_ratio": float(configured.get("max_ratio", thresholds.ratio_threshold_for(protein_type))),
        "mean_abs_log10_error": float(
            configured.get(
                "mean_abs_log10_error",
                thresholds.mean_abs_log10_error_threshold_for(protein_type),
            )
        ),
    }


def _build_comparisons(bench: dict, predicted_ppb: Dict[str, float]) -> List[CompoundComparison]:
    comparisons: List[CompoundComparison] = []
    signal_map = bench.get("measured_volatiles") or bench.get("reference_volatiles") or {}
    for compound, measured in signal_map.items():
        matched_name, predicted_value, match_score = _best_prediction_match(compound, predicted_ppb)
        comparisons.append(
            CompoundComparison(
                compound=compound,
                measured_ppb=float(measured.get("conc_ppb", 0.0)),
                predicted_ppb=float(predicted_value),
                matched_name=matched_name,
                uncertainty_pct=measured.get("uncertainty_pct"),
                match_score=match_score,
            )
        )
    return comparisons


def _run_benchmark_recommendation(
    bench: dict,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thermodynamic_gating: str = "auto",
) -> dict:
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)
    designer = MaillardPipeline(target_tag=target_tag)
    kinetics = KineticsEngine(temperature_k=conditions.temperature_kelvin)
    gating_mode = resolve_thermodynamic_gating_mode(bench, thermodynamic_gating)

    process_state = determine_matrix_process_state(
        temperature_celsius=float(conditions.temperature_celsius),
        time_minutes=float(formulation.get("time_minutes", 60.0)),
        water_activity=conditions.water_activity,
    )
    protein_type = formulation.get("protein_type", "free")
    family_upstream_contract = build_family_upstream_contract(
        sugars=formulation.get("sugars", []),
        amino_acids=formulation.get("amino_acids", []),
        additives=formulation.get("additives", []),
        lipids=formulation.get("lipids", []),
        support_cues=formulation.get("support_cues", []),
        interventions=formulation.get("interventions", []),
        protein_type=protein_type,
        pH=formulation.get("ph", conditions.pH),
        molar_ratios=formulation.get("molar_ratios", {}),
        process_state=process_state,
        temperature_celsius=float(conditions.temperature_celsius),
        time_minutes=float(formulation.get("time_minutes", 60.0)),
        water_activity=conditions.water_activity,
    )
    effective_ratios = dict(family_upstream_contract.get("effective_molar_ratios", {})) or dict(formulation.get("molar_ratios", {}))
    for precursor_name, ratio_value in (family_upstream_contract.get("added_precursor_ratios", {}) or {}).items():
        effective_ratios.setdefault(str(precursor_name), float(ratio_value))
    added_precursors = [
        precursor_name
        for precursor_name in (family_upstream_contract.get("added_precursors", []) or [])
        if precursor_name not in formulation.get("sugars", [])
        + formulation.get("amino_acids", [])
        + formulation.get("additives", [])
        + formulation.get("lipids", [])
    ]

    names = (
        formulation["sugars"]
        + formulation["amino_acids"]
        + formulation.get("additives", [])
        + formulation.get("lipids", [])
        + added_precursors
    )
    precursors = resolve_many(names)
    steps = SmirksEngine(conditions).enumerate(precursors, max_generations=4)

    heuristic_barriers = {}
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        bar, _source, unc = designer.db.get_best_barrier(reactants, products, step.reaction_family or "unknown")
        barrier_for_rate = float(bar)
        if gating_mode == "on":
            thermo = kinetics.get_reaction_thermo(reactants, products, conditions.temperature_kelvin)
            dg = float(thermo.get("delta_g_kcal_mol", 0.0))
            if dg > THERMODYNAMIC_GATING_THRESHOLD_KCAL:
                barrier_for_rate = 99.0
            else:
                barrier_for_rate = max(barrier_for_rate, max(0.5, dg + 0.5))
        k = conditions.get_rate_constant(
            step.reaction_family or "unknown",
            ea_override_kcal=barrier_for_rate,
            reactant_labels=[species.label for species in step.reactants],
        )
        bar_eff = effective_barrier_from_rate_constant(
            k,
            conditions.temperature_kelvin,
            step.reaction_family or "unknown",
        )
        rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)

    from src.recommend import Recommender
    from src.chem_utils import canonicalize_smiles
    from src.lipid_oxidation import predict_hexanal_generation

    initial_concentrations = {}
    ratios = effective_ratios
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[canonicalize_smiles(precursor.smiles, fallback_to_original=True, strip_salts=True)] = qty

    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is not None:
        oxidation = predict_hexanal_generation(
            model["lipid_profile"],
            temp_C=float(conditions.temperature_celsius),
            time_min=float(formulation.get("time_minutes", 60.0)),
            oxygen_availability=1.0,
        )
        # SUBSTRATE CORRECTION 2026-08-27 (Wave P item 4): the load is now
        # PER MARKER, because nonanal comes off the oleate pool and the other
        # three come off the linoleate pool. See
        # `src.lipid_oxidation.MARKER_HYDROPEROXIDE_POOL` for the evidence.
        def _oxidation_load_ppb_for(compound_name: str) -> float:
            return float(oxidation[hydroperoxide_pool_key_for_marker(compound_name)]) * 1000.0

        # Convert ppb to a proxy molar concentration. The framework's recommend engine converts molar -> ppb using ppb_conversion_factor later.
        # But wait, predict_hexanal_generation gives total_hydroperoxide load which drives hexanal and nonanal.
        # Since the matrix_only path uses MATRIX_BENCHMARK_BASE_MARKER_YIELDS, and recommend.py applies its own logic,
        # we can just inject the specific SMILES for Hexanal, Nonanal directly into initial_concentrations so the recommender projects them.
        # Actually, if we inject them, recommender will project them as observables.
        # Let's inject Hexanal (CCCCCC=O), Nonanal (CCCCCCCCC=O), 1-Hexanol (CCCCCCO), 2-Pentylfuran (CCCCC1=CC=CO1)
        # We need their proxy mass to align with the yield factors.
        # Proxy ppb in the matrix_only path was: oxidation_load_ppb * yield_factor.
        # So we want initial_concentrations to be (oxidation_load_ppb * yield_factor / MW) / ppb_conversion_factor
        from src.benchmark_validation import MATRIX_BENCHMARK_BASE_MARKER_YIELDS
        from src.projection import DEFAULT_PROJECTION_STRATEGY
        for compound_name, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            smiles_map = {
                "Hexanal": "CCCCCC=O",
                "Nonanal": "CCCCCCCCC=O",
                "1-Hexanol": "CCCCCCO",
                # Audit 2026-08-26: was CCCCC1=CC=CO1 (2-butylfuran, C8); the curated
                # registry uses the correct C9 pentyl form, so canonical matching failed.
                "2-Pentylfuran": "CCCCCC1=CC=CO1",
            }
            smi = smiles_map.get(compound_name)
            if smi:
                mol = Chem.MolFromSmiles(smi)
                # Full molecular weight (implicit H included) to mirror the reverse
                # molar->ppb transform in recommend.py; the previous heavy-atom-only
                # sum inflated the injected molarity by 11-16%.
                mw = float(Descriptors.MolWt(mol)) if mol else 100.0
                target_proxy_ppb = _oxidation_load_ppb_for(compound_name) * float(yield_factor)
                # target_proxy_ppb = molar_concentration * mw * ppb_conversion_factor
                molar_conc = target_proxy_ppb / (mw * DEFAULT_PROJECTION_STRATEGY.ppb_conversion_factor)
                # initial_concentrations is consumed downstream in mM (projection.py
                # applies limiting_pool_to_molar_factor=1e-3); molar_conc above is in
                # mol/L. The missing x1000 silently suppressed the whole volatile
                # budget by 1000^(k/n) via the geometric-mean pool (audit 2026-08-26).
                conc_mM = molar_conc * 1000.0
                canon = canonicalize_smiles(smi, fallback_to_original=True, strip_salts=True)
                initial_concentrations[canon] = initial_concentrations.get(canon, 0.0) + conc_mM

    rec = Recommender()
    return rec.predict_from_steps(
        steps,
        heuristic_barriers,
        initial_concentrations,
        temperature_kelvin=conditions.temperature_kelvin,
        time_minutes=float(formulation.get("time_minutes", 60.0)),
        protein_type=protein_type,
        denaturation_state=float(bench.get("denaturation_state", 0.5)),
        protein_source=formulation.get("protein_source"),
        family_upstream_contract=family_upstream_contract,
    )


def _run_matrix_only_benchmark_prediction(bench: dict) -> dict:
    protein_type = str(bench.get("protein_type", "free"))
    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is None:
        raise BenchmarkNotSupportedError(f"No matrix-only predictor for protein_type={protein_type}")

    conditions = bench["conditions"]
    water_activity = _get_condition_water_activity(conditions)
    process_state = str((bench.get("process_metadata") or {}).get(
        "state",
        determine_matrix_process_state(
            temperature_celsius=float(conditions["temp_C"]),
            time_minutes=float(conditions["time_min"]),
            water_activity=water_activity,
        ),
    ))
    oxidation = predict_hexanal_generation(
        model["lipid_profile"],
        temp_C=float(conditions["temp_C"]),
        time_min=float(conditions["time_min"]),
        oxygen_availability=1.0,
    )
    # SUBSTRATE CORRECTION 2026-08-27 (Wave P item 4): per-marker oxidation load.
    # Nonanal is the C9 fragment of the OLEATE double bond and is now scaled off
    # `oleic_acid_pct`; hexanal / 2-pentylfuran / 1-hexanol stay on the linoleate
    # pool, which is what Miyazaki 2023's isomer-resolved product lists support.
    # See `src.lipid_oxidation.MARKER_HYDROPEROXIDE_POOL`.
    oxidation_load_ppb_by_pool = {
        key: float(oxidation[key]) * 1000.0 for key in HYDROPEROXIDE_POOL_KEYS.values()
    }
    headspace_model = HeadspaceModel()
    pH = conditions.get("ph")

    # 2026-08-27 (Wave S4): the per-lane binding context, built from
    # data/lit/binding_constants.yml -- i.e. from each lane's OWN source, never from a
    # repository guess. It is inert unless the binding observability mode is selected.
    binding_context = resolve_binding_context(bench)

    predicted_ppb: Dict[str, float] = {}
    predicted_proxy_ppb: Dict[str, float] = {}
    projection_metadata: Dict[str, Dict[str, Any]] = {}
    binding_residual_ratios: Dict[str, float] = {}
    for compound, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
        headspace_factor = headspace_model.get_matrix_benchmark_headspace_factor(
            compound,
            protein_type=protein_type,
            pH=pH,
            temperature_celsius=float(conditions["temp_C"]),
            time_minutes=float(conditions["time_min"]),
            water_activity=water_activity,
            binding_context=binding_context,
        )
        calibration = describe_matrix_calibration(
            compound,
            protein_type=protein_type,
            process_state=process_state,
        )
        panel_entry = get_compound_panel_entry(compound) or {}
        calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
        release_factor = headspace_factor / calibration_factor if calibration_factor > 0.0 else headspace_factor
        marker_load_ppb = oxidation_load_ppb_by_pool[hydroperoxide_pool_key_for_marker(compound)]
        proxy_ppb = marker_load_ppb * float(yield_factor) * release_factor
        observable_ppb = proxy_ppb * calibration_factor
        binding_row: Dict[str, Any] = {}
        if binding_mode_active():
            # 2026-08-27 (Wave S4) -- THE NO-DOUBLE-COUNT CHECK (assembled here, asserted
            # after the loop).
            #
            # `calibration_factor` (the FITTED / back-solved observability constant) is
            # divided out of `release_factor` above and multiplied back in here, so it
            # cancels exactly. In binding mode the surviving net observability must
            # therefore be the BINDING factor times the pH release factor, with NO
            # contribution from the fitted registry and none from the dynamic-retention
            # composition -- the latter matters most, because
            # `compose_dynamic_retention` routes through `resolve_compound_matrix_retention`
            # whose `volatile_retention` is documented as "fraction escaping matrix (rest
            # is bound)", i.e. it is itself an unanchored binding model.
            #
            # The test is the RATIO net / (f_free x pH), collected per compound. It cannot
            # be asserted to equal 1.0 point-by-point because
            # `src.uncertainty_propagation._observable_multipliers` legitimately wraps this
            # method with a sampled scalar during Monte-Carlo propagation. But that scalar
            # is GLOBAL to the sample, while every factor this check is guarding against
            # (registry factor, dynamic retention) is PER COMPOUND. So the invariant that
            # survives the sampler and still catches a leak is: the ratio must be the SAME
            # for every compound on the lane.
            binding = binding_observability_factor(compound, context=binding_context)
            ph_factor = headspace_model.get_matrix_ph_release_factor(
                compound, protein_type=protein_type, pH=pH
            )
            expected_net = float(binding.f_free) * float(ph_factor)
            if expected_net <= 0.0:
                raise AssertionError(
                    f"binding observability for {compound!r} on {bench.get('benchmark_id')!r} "
                    "is non-positive; f_free must be in (0, 1]."
                )
            binding_residual_ratios[compound] = float(headspace_factor) / expected_net
            binding_row = binding.to_dict()
            binding_row["binding_ph_release_factor"] = float(ph_factor)
        predicted_proxy_ppb[compound] = proxy_ppb
        predicted_ppb[compound] = observable_ppb
        projection_metadata[compound] = make_projection_metadata_row(
            compound=compound,
            proxy_ppb=proxy_ppb,
            observable_ppb=observable_ppb,
            extras={
                "matrix_factor": 1.0,
                # Wave P item 4: which fatty-acid pool this marker was cleaved
                # from. Nonanal = oleate, everything else = linoleate.
                "hydroperoxide_pool": MARKER_HYDROPEROXIDE_POOL.get(compound, "linoleate"),
                "oxidation_load_ppb": marker_load_ppb,
                "headspace_factor": release_factor,
                "total_observable_factor": headspace_factor,
                "process_state": process_state,
                "observability_mode": observability_mode(),
                **panel_entry,
                **calibration,
                **binding_row,
            },
        )

    if binding_residual_ratios:
        ratios = list(binding_residual_ratios.values())
        reference = ratios[0]
        for compound, ratio in binding_residual_ratios.items():
            if not math.isclose(ratio, reference, rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(
                    "binding-physics observability mode is active but the net matrix "
                    f"observability on {bench.get('benchmark_id')!r} is not a single global "
                    "multiple of (binding f_free x pH release factor): "
                    f"{binding_residual_ratios!r}. A COMPOUND-SPECIFIC term other than the "
                    "measured binding model is contributing -- most likely the fitted "
                    "registry factor or the dynamic-retention composition. Double counting "
                    "refused."
                )

    return {
        "targets": [],
        "metrics": {
            "matrix_model": protein_type,
            "observability_mode": observability_mode(),
            "binding_no_double_count_ratio": (
                sorted(binding_residual_ratios.values())[0] if binding_residual_ratios else None
            ),
            # Wave P item 4: the LINOLEATE load keeps the historical key, so existing
            # readers keep the meaning they had; the oleate pool that now drives
            # nonanal is reported alongside it.
            "oxidation_load_ppb": oxidation_load_ppb_by_pool["total_hydroperoxide"],
            "oxidation_load_ppb_oleate": oxidation_load_ppb_by_pool["total_hydroperoxide_oleate"],
        },
        "predicted_ppb": predicted_ppb,
        "predicted_proxy_ppb": predicted_proxy_ppb,
        "projection_metadata": projection_metadata,
        "debug_paths": {},
        # 2026-08-27 (Wave S1): mirrors `predict_from_steps`'s per-channel flux breakdown so
        # the two execution paths return the same key set. It is EMPTY here and always will
        # be: the matrix-only lane never enumerates a reaction network, so it has no channels
        # to break down -- which is also why neither Wave S1 fix reaches this path.
        "debug_channel_flux": {},
        "species_names": {compound: compound for compound in predicted_ppb},
    }


def _evaluate_loaded_benchmark(
    bench: dict,
    *,
    bench_path: Path,
    target_tag: str = DEFAULT_TARGET_TAG,
    thermodynamic_gating: str = "auto",
) -> BenchmarkEvaluation:
    formulation = benchmark_to_formulation(bench)
    supported, reason = _is_supported_formulation(formulation)

    if not supported:
        return BenchmarkEvaluation(
            benchmark_id=bench["benchmark_id"],
            bench_file=bench_path,
            supported=False,
            reason=reason,
            predicted_ppb={},
            comparisons=[],
            pearson_r=None,
            mae_ppb=None,
            projection_metadata={},
            reference_signal_origin="measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles",
        )

    metadata = get_benchmark_metadata(bench)
    if metadata.family == "safety":
        # Use the dedicated safety kinetics solver for safety-critical targets
        conditions = benchmark_to_conditions(bench)
        asn_conc = 0.0
        sugar_conc = 0.0
        lysine_conc = 0.0
        for name, data in bench["precursors"].items():
            n_low = name.lower()
            if "asparagine" in n_low or "asn" in n_low:
                asn_conc = data["concentration_mM"]
            if "lysine" in n_low or "lys" in n_low:
                lysine_conc += data["concentration_mM"]
            if any(s in n_low for s in ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar"]):
                sugar_conc += data["concentration_mM"]

        time_min = bench["conditions"].get("time_min", 20.0)
        signal_map = bench.get("measured_volatiles") or bench.get("reference_volatiles") or {}
        signal_names = {str(name).strip().lower() for name in signal_map}
        predicted_safety: Dict[str, float] = {}
        if any("acrylamide" in name for name in signal_names):
            safety_res = predict_acrylamide(
                asparagine_mM=asn_conc,
                reducing_sugar_mM=sugar_conc,
                temp_C=conditions.temperature_celsius,
                time_min=time_min,
                pH=conditions.pH,
                water_activity=conditions.water_activity,
                moisture_regime=getattr(conditions, "moisture_regime", None),
                effective_temp_c=conditions.effective_temperature_celsius,
            )
            predicted_safety["acrylamide"] = safety_res.acrylamide_ppb
        if any(name in {"cml", "n carboxymethyl lysine", "carboxymethyllysine"} or "cml" in name for name in signal_names):
            predicted_safety["Nε-(Carboxymethyl)lysine (CML)"] = predict_cml(
                lysine_mM=lysine_conc,
                reducing_sugar_mM=sugar_conc,
                temp_C=conditions.temperature_celsius,
                time_min=time_min,
                water_activity=conditions.water_activity,
                effective_temp_c=conditions.effective_temperature_celsius,
            )
        if any(name in {"cel", "n carboxyethyl lysine", "carboxyethyllysine"} or "cel" in name for name in signal_names):
            predicted_safety["Nε-(Carboxyethyl)lysine (CEL)"] = predict_cel(
                lysine_mM=lysine_conc,
                reducing_sugar_mM=sugar_conc,
                temp_C=conditions.temperature_celsius,
                time_min=time_min,
                water_activity=conditions.water_activity,
                effective_temp_c=conditions.effective_temperature_celsius,
            )
        if any("furosine" in name for name in signal_names):
            safety_context = bench.get("safety_context") or {}
            damage_protein_type = str(
                safety_context.get("protein_type")
                or bench.get("damage_protein_type")
                or bench.get("protein_type")
                or "free"
            )
            predicted_safety["furosine"] = predict_furosine(
                temp_C=conditions.temperature_celsius,
                time_min=time_min,
                lysine_mM=lysine_conc,
                reducing_sugar_mM=sugar_conc,
                protein_type=damage_protein_type,
                water_activity=conditions.water_activity,
                effective_temp_c=conditions.effective_temperature_celsius,
            )
        rec_result = {"predicted_ppb": predicted_safety, "projection_metadata": {}}
    elif metadata.execution_path == "matrix_only":
        rec_result = _run_matrix_only_benchmark_prediction(bench)
    else:
        rec_result = _run_benchmark_recommendation(
            bench,
            target_tag=target_tag,
            thermodynamic_gating=thermodynamic_gating,
        )
    comparisons = _build_comparisons(bench, rec_result["predicted_ppb"])

    matched_comparisons = [comparison for comparison in comparisons if comparison.matched_name is not None]
    measured_values = [comparison.measured_ppb for comparison in matched_comparisons]
    predicted_values = [comparison.predicted_ppb for comparison in matched_comparisons]
    pearson_r = _pearson(measured_values, predicted_values)
    mae_ppb = None
    if matched_comparisons:
        mae_ppb = sum(abs(measured - predicted) for measured, predicted in zip(measured_values, predicted_values)) / len(matched_comparisons)

    return BenchmarkEvaluation(
        benchmark_id=bench["benchmark_id"],
        bench_file=bench_path,
        supported=True,
        reason=None,
        predicted_ppb=rec_result["predicted_ppb"],
        comparisons=comparisons,
        pearson_r=pearson_r,
        mae_ppb=mae_ppb,
        projection_metadata=rec_result.get("projection_metadata", {}),
        reference_signal_origin="measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles",
    )


def evaluate_benchmark(
    bench_file: Path | str,
    target_tag: str = DEFAULT_TARGET_TAG,
    thermodynamic_gating: str = "auto",
) -> BenchmarkEvaluation:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
    return _evaluate_loaded_benchmark(
        bench,
        bench_path=bench_path,
        target_tag=target_tag,
        thermodynamic_gating=thermodynamic_gating,
    )


def evaluate_benchmark_payload(
    bench: dict,
    *,
    benchmark_id: Optional[str] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
    thermodynamic_gating: str = "auto",
) -> BenchmarkEvaluation:
    normalized = dict(bench)
    if benchmark_id and not normalized.get("benchmark_id"):
        normalized["benchmark_id"] = benchmark_id
    payload_id = str(normalized.get("benchmark_id", benchmark_id or "matrix_experiment_payload"))
    pseudo_path = data_paths.VALIDATION_DIR / f"{payload_id}.synthetic_benchmark.json"
    return _evaluate_loaded_benchmark(
        normalized,
        bench_path=pseudo_path,
        target_tag=target_tag,
        thermodynamic_gating=thermodynamic_gating,
    )


def summarize_evaluation_for_benchmark(
    evaluation: BenchmarkEvaluation,
    bench: dict,
    *,
    protein_type: str = "free",
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> BenchmarkSummary:
    matched = [comparison for comparison in evaluation.comparisons if comparison.matched_name is not None]
    ratios = [comparison.ratio for comparison in matched if math.isfinite(comparison.ratio)]
    max_ratio = max(ratios) if ratios else None
    mean_ratio = sum(ratios) / len(ratios) if ratios else None
    mean_abs_log10_error = _mean_abs_log10_error(matched)

    scale_thresholds = _resolve_scale_thresholds(
        bench,
        protein_type=protein_type,
        thresholds=thresholds,
    )
    ratio_threshold = scale_thresholds["max_ratio"]
    log_error_threshold = scale_thresholds["mean_abs_log10_error"]
    metadata = get_benchmark_metadata(bench)
    conditions = bench.get("conditions", {})
    ranking_contract = _evaluate_matrix_ranking_contract(bench, evaluation.predicted_ppb) if metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"} else {
        "status": "n/a",
        "ranked_observable_targets": [],
        "adverse_markers": [],
    }
    matrix_contract = get_matrix_ranking_contract(bench) if metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"} else {}
    process_state = matrix_contract.get("process_state")

    if not evaluation.supported:
        return BenchmarkSummary(
            benchmark_id=evaluation.benchmark_id,
            bench_file=evaluation.bench_file,
            tier=metadata.tier,
            family=metadata.family,
            execution_path=metadata.execution_path,
            benchmark_engine=metadata.benchmark_engine,
            comparator_signal=metadata.comparator_signal,
            cantera_role=metadata.cantera_role,
            target_snapshot_policy=metadata.target_snapshot_policy,
            thermodynamic_gating_policy=metadata.thermodynamic_gating_policy,
            supported=False,
            reason=evaluation.reason,
            protein_type=protein_type,
            coverage=0.0,
            matched_compounds=0,
            total_compounds=0,
            pearson_r=None,
            mae_ppb=None,
            max_ratio=None,
            mean_ratio=None,
            ranking_status="unsupported",
            scale_status="unsupported",
            overall_status="unsupported",
            strict_ready=False,
            blocking_issues=[evaluation.reason or "unsupported"],
            conditions=conditions,
            process_state=process_state,
            ranked_observable_targets=list(ranking_contract.get("ranked_observable_targets", [])),
            adverse_markers=list(ranking_contract.get("adverse_markers", [])),
            ranking_contract_status=str(ranking_contract.get("status", "n/a")),
            calibration_mode=matrix_contract.get("calibration_mode"),
            reference_signal_origin=evaluation.reference_signal_origin,
            mean_abs_log10_error=None,
        )

    if len(matched) >= thresholds.min_matched_for_ranking and evaluation.pearson_r is not None:
        ranking_status = "pass" if evaluation.pearson_r >= thresholds.ranking_threshold else "fail"
    elif len(matched) > 0:
        ranking_status = "n/a"
    else:
        ranking_status = "fail"

    if max_ratio is None:
        scale_status = "fail"
    elif max_ratio <= ratio_threshold and (
        mean_abs_log10_error is None or mean_abs_log10_error <= log_error_threshold
    ):
        scale_status = "pass"
    else:
        scale_status = "fail"

    overall_status = "pass"
    if evaluation.coverage < thresholds.full_coverage_threshold:
        overall_status = "coverage-gap"
    elif ranking_status == "fail":
        overall_status = "ranking-gap"
    elif scale_status == "fail":
        overall_status = "scale-gap"
    elif ranking_status == "n/a":
        overall_status = "pass-no-ranking"

    blocking_issues: List[str] = []
    if evaluation.coverage < thresholds.full_coverage_threshold:
        blocking_issues.append(
            f"coverage {evaluation.coverage:.1%} < {thresholds.full_coverage_threshold:.0%}"
        )
    if ranking_status == "fail":
        pearson = "n/a" if evaluation.pearson_r is None else f"{evaluation.pearson_r:.3f}"
        blocking_issues.append(
            f"ranking {pearson} < {thresholds.ranking_threshold:.2f}"
        )
    if scale_status == "fail":
        ratio_value = "n/a" if max_ratio is None else f"{max_ratio:.3f}"
        blocking_issues.append(
            f"max ratio {ratio_value} > {ratio_threshold:.2f}"
        )
        if mean_abs_log10_error is not None and mean_abs_log10_error > log_error_threshold:
            blocking_issues.append(
                f"mean |log10 ratio| {mean_abs_log10_error:.3f} > {log_error_threshold:.3f}"
            )

    strict_ready = (
        evaluation.coverage >= thresholds.full_coverage_threshold
        and ranking_status != "fail"
        and scale_status == "pass"
    )
    if not DEFAULT_VALIDATION_CONTRACT.is_strict_gate_eligible(
        tier=metadata.tier,
        execution_path=metadata.execution_path,
    ):
        strict_ready = False
        if not blocking_issues:
            blocking_issues.append("matrix-only intake path is executable but not yet in the strict release gate")

    if metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"} and ranking_contract.get("status") not in {"pass", "n/a"}:
        blocking_issues.append(f"matrix ranking contract: {ranking_contract.get('status')}")

    return BenchmarkSummary(
        benchmark_id=evaluation.benchmark_id,
        bench_file=evaluation.bench_file,
        tier=metadata.tier,
        family=metadata.family,
        execution_path=metadata.execution_path,
        benchmark_engine=metadata.benchmark_engine,
        comparator_signal=metadata.comparator_signal,
        cantera_role=metadata.cantera_role,
        target_snapshot_policy=metadata.target_snapshot_policy,
        thermodynamic_gating_policy=metadata.thermodynamic_gating_policy,
        supported=True,
        reason=None,
        protein_type=protein_type,
        coverage=evaluation.coverage,
        matched_compounds=len(matched),
        total_compounds=len(evaluation.comparisons),
        pearson_r=evaluation.pearson_r,
        mae_ppb=evaluation.mae_ppb,
        max_ratio=max_ratio,
        mean_ratio=mean_ratio,
        ranking_status=ranking_status,
        scale_status=scale_status,
        overall_status=overall_status,
        strict_ready=strict_ready,
        blocking_issues=blocking_issues,
        conditions=conditions,
        process_state=process_state,
        ranked_observable_targets=list(ranking_contract.get("ranked_observable_targets", [])),
        adverse_markers=list(ranking_contract.get("adverse_markers", [])),
        ranking_contract_status=str(ranking_contract.get("status", "n/a")),
        calibration_mode=matrix_contract.get("calibration_mode"),
        reference_signal_origin=evaluation.reference_signal_origin,
        mean_abs_log10_error=mean_abs_log10_error,
        # 2026-08-27 (Wave I): the honesty column. `overall_status` above is mechanical --
        # it says whether predictions and measurements agree. This says whether that
        # agreement is evidence about the model or a recovery of the model's own inputs.
        evidence_role=benchmark_evidence_role(evaluation.benchmark_id, evaluation.bench_file),
    )


def summarize_evaluation(
    evaluation: BenchmarkEvaluation,
    *,
    protein_type: str = "free",
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> BenchmarkSummary:
    bench = load_benchmark(evaluation.bench_file)
    return summarize_evaluation_for_benchmark(
        evaluation,
        bench,
        protein_type=protein_type,
        thresholds=thresholds,
    )


def _projection_metadata_for_match(
    evaluation: BenchmarkEvaluation,
    comparison: CompoundComparison,
) -> Dict[str, Any]:
    if comparison.matched_name and comparison.matched_name in evaluation.projection_metadata:
        return evaluation.projection_metadata[comparison.matched_name]
    normalized = _normalize_name(comparison.compound)
    for key, meta in evaluation.projection_metadata.items():
        if _normalize_name(str(key)) == normalized or _normalize_name(str(meta.get("compound", ""))) == normalized:
            return meta
    return {}


def build_matrix_benchmark_deltas(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[MatrixBenchmarkDelta]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkDelta] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue
        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation,
            protein_type=bench.get("protein_type", "free"),
        )
        role_lookup = {
            str(item.get("name", "")).strip().lower(): str(item.get("role", "target"))
            for item in get_matrix_ranking_contract(bench).get("observable_targets", [])
        }
        adverse = {str(item).strip().lower() for item in summary.adverse_markers}
        for comparison in evaluation.comparisons:
            meta = _projection_metadata_for_match(evaluation, comparison)
            compound_key = comparison.compound.strip().lower()
            role = role_lookup.get(compound_key, "adverse_marker" if compound_key in adverse else "reference")
            abs_delta = abs(float(comparison.measured_ppb) - float(comparison.predicted_ppb))
            pct_delta = None if comparison.measured_ppb <= 0.0 else abs_delta / float(comparison.measured_ppb)
            rows.append(
                MatrixBenchmarkDelta(
                    benchmark_id=evaluation.benchmark_id,
                    bench_file=evaluation.bench_file,
                    protein_type=bench.get("protein_type", "free"),
                    execution_path=metadata.execution_path,
                    process_state=summary.process_state,
                    reference_signal_origin=evaluation.reference_signal_origin,
                    ranking_contract_status=summary.ranking_contract_status,
                    compound=comparison.compound,
                    role=role,
                    reference_ppb=float(comparison.measured_ppb),
                    predicted_ppb=float(comparison.predicted_ppb),
                    abs_delta_ppb=abs_delta,
                    pct_delta=pct_delta,
                    ratio=float(comparison.ratio),
                    calibration_source=str(meta.get("calibration_source", "class_fallback")),
                    calibration_evidence_strength=str(meta.get("calibration_evidence_strength", "heuristic")),
                    calibration_fallback_mode=str(meta.get("calibration_fallback_mode", "class_level")),
                )
            )
    return rows


def matrix_source_anchor(bench: Mapping[str, Any]) -> str:
    """Return the benchmark's citable external anchor, or "" if it has none.

    2026-08-27 (audit remediation, DOI-less identifier retyping): four sources in
    this repo are genuinely DOI-less (two theses, a US patent, a journal with no
    DOI registration) and their identifiers used to be stored in fields named
    ``doi`` / ``source_doi``, which is both dishonest and a citation-gate
    violation. They now carry a typed ``identifier`` + ``identifier_scheme``
    pair. A typed non-DOI identifier is exactly as much an external anchor as a
    DOI, so every consumer that tested ``source_doi`` truthiness reads this
    helper instead — otherwise retyping an identifier would silently *downgrade*
    a real external source to "unspecified origin".
    """
    doi = str(bench.get("source_doi") or "").strip()
    if doi:
        return doi
    identifier = str(bench.get("identifier") or "").strip()
    if identifier:
        return identifier
    return ""


# 2026-08-27 (Wave I). Canonical home of the signal-origin classifier. It used to live
# only in `uncertainty_propagation._benchmark_signal_origin`, but the benchmark summary
# layer needs it too and `uncertainty_propagation` imports THIS module, so keeping the
# implementation there would have meant a cycle or a second, drifting copy.
SYNTHETIC_BENCHMARK_ORIGINS = frozenset(
    {
        "synthetic_diagnostic",
        "internal_reproducibility_candidate",
        "internal_experiment",
    }
)


def benchmark_signal_origin(bench_file: Path) -> str:
    """Classify a benchmark's comparator signal as literature-measured or internal/synthetic.

    Returns ``"external_literature"`` or ``"internal_synthetic"``.
    """
    try:
        bench = json.loads(Path(bench_file).read_text())
    except (OSError, json.JSONDecodeError):
        return "internal_synthetic"
    origin = str((bench.get("source_metadata") or {}).get("origin", "")).strip().lower()
    if origin in SYNTHETIC_BENCHMARK_ORIGINS:
        return "internal_synthetic"
    # `matrix_source_anchor` also accepts the typed `identifier`/`identifier_scheme` pair
    # that DOI-less sources (theses, patents, DOI-free journals) now carry, so retyping an
    # identifier out of a `doi`-named field cannot reclassify a literature row as
    # internal/synthetic.
    if matrix_source_anchor(bench) or origin.startswith("external"):
        return "external_literature"
    return "internal_synthetic"


def benchmark_evidence_role(benchmark_id: str, bench_file: Path) -> str:
    """What KIND of claim this benchmark can support: see BenchmarkSummary.evidence_role.

    2026-08-27 (Wave I). The mechanical status says whether predictions and measurements
    agree; this says whether that agreement is evidence. Fit-recovery is checked FIRST and
    beats signal origin: a lane whose observable factor was solved from a literature
    benchmark is still literature-sourced, and is still not a prediction.

    Two independent sources of fit-recovery are consulted:

    * the calibration registry's own `fitted_to_benchmark` declarations (the matrix
      observability factors, one free factor per compound per lane); and
    * `src.fit_target_index`, which reads the fit records under `results/validation/` and
      classifies each by LEVERAGE. Only high-leverage fits (enough free parameters to
      reproduce their target row by row) make a benchmark fit-recovery. A low-leverage
      global fit -- two constants across two dozen rows -- does NOT, because excluding
      those rows would delete genuine failures from the count rather than expose them.
    """
    if is_fit_recovery_benchmark(benchmark_id):
        return "fit_recovery"
    if is_per_row_fit_target(benchmark_id):
        return "fit_recovery"
    if benchmark_signal_origin(Path(bench_file)) == "internal_synthetic":
        return "internal_synthetic"
    return "predictive"


def _matrix_source_origin(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    origin = str(source_metadata.get("origin", "")).strip()
    if origin:
        return origin
    if matrix_source_anchor(bench):
        return "external_literature"
    return "unspecified"


def _matrix_source_reference(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    anchor = matrix_source_anchor(bench)
    if anchor:
        return anchor
    generator = str(source_metadata.get("generator", "")).strip()
    origin = str(source_metadata.get("origin", "")).strip()
    if origin and generator:
        return f"{origin}:{generator}"
    if origin:
        return origin
    if generator:
        return generator
    return "unspecified"


def _matrix_target_profile(bench: dict) -> str:
    contract = get_matrix_ranking_contract(bench)
    roles = {str(item.get("role", "")).strip().lower() for item in contract.get("observable_targets", [])}
    has_desirable = "desirable_marker" in roles
    has_adverse = bool(contract.get("adverse_markers"))
    if has_desirable and has_adverse:
        return "mixed"
    if has_desirable:
        return "meaty_positive"
    if has_adverse:
        return "adverse_only"
    return "untyped"


def _is_internal_measured_matrix_source(source_origin: str) -> bool:
    return str(source_origin).strip().lower() == "internal_experiment"


def _is_internal_candidate_support_status(support_status: str) -> bool:
    return str(support_status).strip().lower() in {
        "internal_candidate",
        "internal_measured_candidate",
        "internal_reference_candidate",
    }


def _init_matrix_support_counts() -> Dict[str, int]:
    return {
        "quantitative_closed": 0,
        "internal_candidate": 0,
        "internal_measured_candidate": 0,
        "internal_reference_candidate": 0,
        "directional_support": 0,
        "open_gap": 0,
    }


def _increment_matrix_support_counts(counts: Dict[str, int], support_status: str) -> None:
    normalized = str(support_status).strip().lower()
    if normalized not in counts:
        counts[normalized] = 0
    counts[normalized] += 1
    if normalized in {"internal_measured_candidate", "internal_reference_candidate"}:
        counts["internal_candidate"] = int(counts.get("internal_candidate", 0)) + 1


def _matrix_external_data_status(bench: dict) -> str:
    # This module is the single live implementation. A drifted duplicate family
    # (benchmark_registry/evaluator/reporting/assertions/markdown) was confirmed
    # dead code and deleted on 2026-08-27 — see tasks/audit_remediation.md.
    evidence_class = str(
        (bench.get("source_metadata") or {}).get(
            "evidence_class",
            (bench.get("metadata") or {}).get("evidence_class", bench.get("evidence_class", "calibration_candidate")),
        )
    ).strip()
    if evidence_class == "external_validation_only":
        return "external_validation_only"

    source_origin = _matrix_source_origin(bench)
    has_measured = bool(bench.get("measured_volatiles"))
    # 2026-08-27: reads the typed identifier as well as `source_doi` (see
    # matrix_source_anchor) so a DOI-less external source is not downgraded.
    if has_measured and (matrix_source_anchor(bench) or source_origin.startswith("external")):
        return "external_quantitative"
    if has_measured and source_origin.strip().lower() == "synthetic_diagnostic":
        # Audit 2026-08-26: payloads whose "measured" values are frozen model
        # output (e.g. the ProtocolPilot intakes) must surface as synthetic,
        # never as measured evidence of any grade.
        return _CompatibilityStatus(
            "synthetic_diagnostic_only",
            aliases=("quantitative_unspecified_origin",),
        )
    if has_measured and _is_internal_measured_matrix_source(source_origin):
        return _CompatibilityStatus(
            "internal_measured_quantitative",
            aliases=("quantitative_unspecified_origin",),
        )
    if has_measured:
        return "quantitative_unspecified_origin"
    if bench.get("reference_volatiles"):
        return "internal_reference_only"
    return "no_comparator_signal"


def assess_matrix_benchmark_evidence(bench: dict | Path | str) -> MatrixBenchmarkEvidence:
    if isinstance(bench, (Path, str)):
        bench_path = Path(bench)
        bench = load_benchmark(bench_path)
    else:
        bench_path = data_paths.benchmark_path(bench.get('benchmark_id', 'unknown'))

    metadata = get_benchmark_metadata(bench)
    process_state = get_matrix_ranking_contract(bench).get("process_state")
    reference_signal_origin = "measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles"
    source_origin = _matrix_source_origin(bench)
    target_profile = _matrix_target_profile(bench)
    external_data_status = _matrix_external_data_status(bench)
    promotable = (
        metadata.execution_path in {"matrix_only", "matrix_precursor_augmented"}
        and target_profile in {"meaty_positive", "mixed"}
        and external_data_status == "external_quantitative"
        and reference_signal_origin == "measured_volatiles"
    )

    if target_profile == "adverse_only":
        blocker = "benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present"
    elif external_data_status == "external_validation_only":
        blocker = "external-validation hold-out only; explicitly excluded from calibration and promotion"
    elif external_data_status == "synthetic_diagnostic_only":
        blocker = "missing external quantitative matrix evidence for meaty-positive targets; current comparator is synthetic model output (diagnostic only)"
    elif external_data_status == "internal_measured_quantitative":
        blocker = "missing external quantitative matrix evidence for meaty-positive targets; current comparator is an internal measured experiment"
    elif external_data_status == "internal_reference_only":
        blocker = "missing external quantitative matrix evidence for meaty-positive targets; current comparator is internal reference-only"
    elif external_data_status != "external_quantitative":
        blocker = "missing external quantitative matrix evidence for meaty-positive targets"
    elif reference_signal_origin != "measured_volatiles":
        blocker = "comparator signal is not wet-lab measured_volatiles"
    else:
        blocker = ""

    return MatrixBenchmarkEvidence(
        benchmark_id=str(bench.get("benchmark_id", bench_path.stem)),
        bench_file=bench_path,
        protein_type=str(bench.get("protein_type", "free")),
        execution_path=metadata.execution_path,
        process_state=process_state,
        reference_signal_origin=reference_signal_origin,
        source_origin=source_origin,
        source_reference=_matrix_source_reference(bench),
        target_profile=target_profile,
        external_data_status=external_data_status,
        promotable=promotable,
        promotion_blocker=blocker,
    )


def build_matrix_benchmark_evidence_audit(
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> List[MatrixBenchmarkEvidence]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkEvidence] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue
        rows.append(assess_matrix_benchmark_evidence(Path(bench_file)))
    return rows





def _matrix_assertion_thresholds(
    bench: dict,
    *,
    protein_type: str,
    thresholds: BenchmarkThresholds,
) -> Dict[str, float]:
    contract = bench.get("matrix_ranking_contract") or {}
    # P1: Check for either the legacy assertion_thresholds or the new validation_contract.scale_thresholds
    configured = contract.get("assertion_thresholds") or contract.get("validation_contract", {}).get("scale_thresholds") or {}
    observable_targets = get_matrix_ranking_contract(bench).get("observable_targets", [])
    return {
        "min_coverage": float(configured.get("min_coverage", thresholds.full_coverage_threshold)),
        "top_k": float(configured.get("top_k", len(observable_targets))),
        "max_ratio": float(configured.get("max_ratio", thresholds.ratio_threshold_for(protein_type))),
    }


def _matched_contract_prediction_rows(
    evaluation: BenchmarkEvaluation,
    contract_names: Iterable[str],
) -> List[tuple[str, float]]:
    rows: List[tuple[str, float]] = []
    for contract_name in contract_names:
        matched_name, predicted_value, _score = _best_prediction_match(str(contract_name), evaluation.predicted_ppb)
        if matched_name is None:
            continue
        rows.append((str(contract_name), float(predicted_value)))
    return rows


def build_matrix_benchmark_assertions(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> List[MatrixBenchmarkAssertion]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    rows: List[MatrixBenchmarkAssertion] = []
    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation,
            protein_type=bench.get("protein_type", "free"),
            thresholds=thresholds,
        )
        evidence = assess_matrix_benchmark_evidence(bench_file)
        assertion_thresholds = _matrix_assertion_thresholds(
            bench,
            protein_type=str(bench.get("protein_type", "free")),
            thresholds=thresholds,
        )
        top_k = max(0, int(assertion_thresholds["top_k"]))
        matched_ranked_targets = _matched_contract_prediction_rows(
            evaluation,
            summary.ranked_observable_targets,
        )
        ranked_targets = [name for name, _value in matched_ranked_targets[:top_k]]
        predicted_top_k = {
            name for name, _value in sorted(matched_ranked_targets, key=lambda row: row[1], reverse=True)[:top_k]
        }
        top_k_hits = sum(1 for name in ranked_targets if name in predicted_top_k)
        if top_k == 0:
            top_k_status = "n/a"
        elif top_k_hits == top_k and summary.ranking_contract_status == "pass":
            top_k_status = "pass"
        else:
            top_k_status = "fail"

        adverse_rows = _matched_contract_prediction_rows(
            evaluation,
            summary.adverse_markers,
        )
        adverse_markers = [name for name in summary.adverse_markers]
        if len(adverse_markers) <= 1:
            adverse_order_status = "n/a"
        else:
            predicted_adverse = [name for name, _value in sorted(adverse_rows, key=lambda row: row[1], reverse=True)]
            if len(predicted_adverse) != len(adverse_markers):
                adverse_order_status = "missing"
            else:
                adverse_order_status = "pass" if predicted_adverse == adverse_markers else "fail"

        min_coverage = float(assertion_thresholds["min_coverage"])
        ratio_tolerance = float(assertion_thresholds["max_ratio"])
        if summary.max_ratio is None:
            ratio_status = "fail"
        elif summary.max_ratio <= ratio_tolerance:
            ratio_status = "pass"
        else:
            ratio_status = "fail"

        overall_pass = (
            summary.coverage >= min_coverage
            and top_k_status in {"pass", "n/a"}
            and adverse_order_status in {"pass", "n/a"}
            and ratio_status == "pass"
            and summary.ranking_contract_status == "pass"
        )
        blocker = evidence.promotion_blocker or "matrix strict gate remains disabled by contract"
        rows.append(
            MatrixBenchmarkAssertion(
                benchmark_id=summary.benchmark_id,
                bench_file=summary.bench_file,
                protein_type=summary.protein_type,
                execution_path=summary.execution_path,
                process_state=summary.process_state,
                target_profile=evidence.target_profile,
                ranking_contract_status=summary.ranking_contract_status,
                coverage=summary.coverage,
                min_coverage=min_coverage,
                top_k=top_k,
                top_k_hits=top_k_hits,
                top_k_status=top_k_status,
                adverse_order_status=adverse_order_status,
                max_ratio=summary.max_ratio,
                ratio_tolerance=ratio_tolerance,
                ratio_status=ratio_status,
                overall_status="pass" if overall_pass else "fail",
                strict_gate_blocked=True,
                blocker=blocker,
            )
        )
    return rows


def render_matrix_benchmark_assertions_markdown(rows: Iterable[MatrixBenchmarkAssertion]) -> str:
    assertion_rows = list(rows)
    lines = [
        "# Matrix Benchmark Assertions",
        "",
        "| Benchmark | Protein | Path | Process State | Target Profile | Ranking Contract | Coverage | Min Coverage | Top-k | Top-k Hits | Top-k Status | Adverse Order | Max Ratio | Ratio Tol. | Ratio Status | Overall | Strict Gate Blocked | Blocker |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in assertion_rows:
        max_ratio = f"{row.max_ratio:.3f}" if row.max_ratio is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.protein_type} | {row.execution_path} | {row.process_state or 'n/a'} | {row.target_profile} | {row.ranking_contract_status} | {row.coverage:.3f} | {row.min_coverage:.3f} | {row.top_k} | {row.top_k_hits} | {row.top_k_status} | {row.adverse_order_status} | {max_ratio} | {row.ratio_tolerance:.3f} | {row.ratio_status} | {row.overall_status} | {'yes' if row.strict_gate_blocked else 'no'} | {row.blocker} |"
        )
    lines.extend([
        "",
        f"Benchmarks asserted: {len(assertion_rows)}",
        f"Assertion passes: {sum(1 for row in assertion_rows if row.overall_status == 'pass')}",
        "Strict gate remains blocked for all matrix benchmarks by contract until external evidence exists.",
    ])
    return "\n".join(lines) + "\n"


def build_matrix_promotion_family_status(
    benchmark_files: Optional[Iterable[Path | str]] = None,
) -> List[MatrixPromotionFamilyStatus]:
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    proteins = sorted({row.protein_type for row in evidence_rows})
    rows: List[MatrixPromotionFamilyStatus] = []
    for protein_type in proteins:
        subset = [row for row in evidence_rows if row.protein_type == protein_type]
        off_flavour_anchor_count = sum(
            1 for row in subset if row.target_profile == "adverse_only" and row.external_data_status == "external_quantitative"
        )
        meaty_candidate_count = sum(
            1 for row in subset if row.target_profile in {"meaty_positive", "mixed"}
        )
        external_meaty_anchor_count = sum(
            1
            for row in subset
            if row.target_profile in {"meaty_positive", "mixed"}
            and row.external_data_status == "external_quantitative"
            and row.reference_signal_origin == "measured_volatiles"
        )
        candidate_set_ready = off_flavour_anchor_count > 0 and meaty_candidate_count > 0
        external_assessment_unlocked = off_flavour_anchor_count > 0 and external_meaty_anchor_count > 0
        if off_flavour_anchor_count == 0:
            blocker = "missing external off-flavour anchor"
        elif meaty_candidate_count == 0:
            blocker = "missing meaty-positive benchmark candidate"
        elif external_meaty_anchor_count == 0:
            blocker = "missing external meaty-positive benchmark"
        else:
            blocker = "none"
        rows.append(
            MatrixPromotionFamilyStatus(
                protein_type=protein_type,
                off_flavour_anchor_count=off_flavour_anchor_count,
                meaty_candidate_count=meaty_candidate_count,
                external_meaty_anchor_count=external_meaty_anchor_count,
                candidate_set_ready=candidate_set_ready,
                external_assessment_unlocked=external_assessment_unlocked,
                blocker=blocker,
            )
        )
    return rows


def _matrix_compound_support_status(
    *,
    evidence_state: str,
    calibration_evidence_strength: str,
    reference_signal_origin: str,
    source_origin: str,
) -> str:
    evidence = str(evidence_state).strip().lower()
    strength = str(calibration_evidence_strength).strip().lower()
    signal_origin = str(reference_signal_origin).strip().lower()
    origin = str(source_origin).strip().lower()

    if evidence == "externally_benchmarked" and signal_origin == "measured_volatiles" and origin.startswith("external"):
        return "quantitative_closed"
    if signal_origin == "reference_volatiles":
        return "internal_reference_candidate"
    if origin == "synthetic_diagnostic":
        # Frozen model output carries reference-grade support at most, never
        # measured-grade (audit 2026-08-26).
        return "internal_reference_candidate"
    if evidence in {"internally_benchmarked", "conditional_calibration"} and _is_internal_measured_matrix_source(origin):
        return "internal_measured_candidate"
    if evidence in {"internally_benchmarked", "conditional_calibration"}:
        return "internal_candidate"
    if evidence in {"transferred_prior", "safety_reference"} or strength in {"class_anchored", "directional_transferred"}:
        return "directional_support"
    return "open_gap"


def build_matrix_target_status_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    benchmark_rows: List[Dict[str, Any]] = []
    support_totals = _init_matrix_support_counts()

    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"matrix_only", "matrix_precursor_augmented"}:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        evidence = assess_matrix_benchmark_evidence(bench_file)
        contract = get_matrix_ranking_contract(bench)
        adverse_markers = {str(item).strip().lower() for item in contract.get("adverse_markers", [])}
        compounds: List[Dict[str, Any]] = []
        benchmark_counts = _init_matrix_support_counts()

        for item in contract.get("observable_targets", []):
            compound_name = str(item.get("name", "")).strip()
            if not compound_name:
                continue
            meta = _projection_metadata_for_match(
                evaluation,
                CompoundComparison(
                    compound=compound_name,
                    measured_ppb=0.0,
                    predicted_ppb=0.0,
                    matched_name=None,
                    uncertainty_pct=None,
                ),
            )
            role = str(item.get("role", "adverse_marker" if compound_name.lower() in adverse_markers else "desirable_marker"))
            support_status = _matrix_compound_support_status(
                evidence_state=str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
                calibration_evidence_strength=str(meta.get("calibration_evidence_strength", "heuristic")),
                reference_signal_origin=summary.reference_signal_origin,
                source_origin=evidence.source_origin,
            )
            _increment_matrix_support_counts(benchmark_counts, support_status)
            _increment_matrix_support_counts(support_totals, support_status)
            compounds.append(
                {
                    "compound": compound_name,
                    "role": role,
                    "target_class": str(meta.get("target_class", item.get("target_class", "unknown"))),
                    "evidence_state": str(meta.get("evidence_state", item.get("evidence_state", "still_missing"))),
                    "calibration_source": str(meta.get("calibration_source", "unknown")),
                    "calibration_evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                    "support_status": support_status,
                }
            )

        external_decision_ready = (
            evidence.target_profile in {"mixed", "meaty_positive"}
            and summary.ranking_contract_status == "pass"
            and benchmark_counts["quantitative_closed"] >= 2
            and benchmark_counts["internal_candidate"] == 0
            and benchmark_counts["directional_support"] == 0
        )
        mechanistic_priority_ready = (
            evidence.target_profile in {"mixed", "meaty_positive"}
            and summary.ranking_contract_status == "pass"
            and not external_decision_ready
            and (benchmark_counts["internal_candidate"] + benchmark_counts["directional_support"]) >= 1
        )
        if evidence.target_profile not in {"mixed", "meaty_positive"}:
            promotion_blocker = "benchmark lacks meaty-positive targets"
        elif summary.ranking_contract_status != "pass":
            promotion_blocker = "ranking contract not yet passing"
        elif benchmark_counts["quantitative_closed"] < 2:
            if evidence.external_data_status == "synthetic_diagnostic_only":
                promotion_blocker = "insufficient externally measured target closure; current comparator is synthetic model output (diagnostic only)"
            elif evidence.external_data_status == "internal_measured_quantitative":
                promotion_blocker = "insufficient externally measured target closure; current comparator is internal measured only"
            elif evidence.external_data_status == "internal_reference_only":
                promotion_blocker = "insufficient externally measured target closure; current comparator is internal reference-only"
            else:
                promotion_blocker = "insufficient externally measured target closure"
        elif benchmark_counts["internal_candidate"] > 0 or benchmark_counts["directional_support"] > 0:
            if benchmark_counts["internal_measured_candidate"] > 0 and benchmark_counts["internal_reference_candidate"] == 0 and benchmark_counts["directional_support"] == 0:
                promotion_blocker = "depends on internally measured support"
            elif benchmark_counts["internal_reference_candidate"] > 0 and benchmark_counts["internal_measured_candidate"] == 0 and benchmark_counts["directional_support"] == 0:
                promotion_blocker = "depends on internal reference-only support"
            else:
                promotion_blocker = "depends on internal or transferred support"
        else:
            promotion_blocker = "none"

        if external_decision_ready:
            next_best_action = "use_for_external_decision"
        elif mechanistic_priority_ready:
            next_best_action = "prioritize_mechanistic_refinement"
        elif evidence.target_profile in {"mixed", "meaty_positive"}:
            next_best_action = "seek_external_data"
        else:
            next_best_action = "retain_as_adverse_anchor"

        benchmark_rows.append(
            {
                "benchmark_id": summary.benchmark_id,
                "bench_file": str(summary.bench_file),
                "protein_type": summary.protein_type,
                "execution_path": summary.execution_path,
                "process_state": summary.process_state,
                "target_profile": evidence.target_profile,
                "external_data_status": evidence.external_data_status,
                "reference_signal_origin": summary.reference_signal_origin,
                "source_origin": evidence.source_origin,
                "ranking_contract_status": summary.ranking_contract_status,
                "support_counts": benchmark_counts,
                "quantitative_support_ready": benchmark_counts["quantitative_closed"] > 0,
                "promotion_ready": external_decision_ready,
                "mechanistic_priority_ready": mechanistic_priority_ready,
                "promotion_blocker": promotion_blocker,
                "next_best_action": next_best_action,
                "compounds": compounds,
            }
        )

    return {
        "schema_version": "1.0",
        "description": "Matrix target support status artifact distinguishing external decision-ready support from mechanistic-priority candidates and unresolved external gaps.",
        "benchmarks": benchmark_rows,
        "summary": {
            "total_benchmarks": len(benchmark_rows),
            "quantitative_support_ready": sum(1 for row in benchmark_rows if row["quantitative_support_ready"]),
            "promotion_ready": sum(1 for row in benchmark_rows if row["promotion_ready"]),
            "mechanistic_priority_ready": sum(1 for row in benchmark_rows if row["mechanistic_priority_ready"]),
            **support_totals,
        },
    }



def _matrix_promotion_requirement_rows(
    benchmark_row: Mapping[str, Any],
    evidence_row: MatrixBenchmarkEvidence,
) -> List[Dict[str, Any]]:
    support_counts = benchmark_row.get("support_counts", {})
    requirement_rows = [
        {
            "key": "meaty_positive_targets_present",
            "label": "Target profile includes meaty-positive compounds",
            "passed": benchmark_row.get("target_profile") in {"mixed", "meaty_positive"},
            "detail": str(benchmark_row.get("target_profile", "unknown")),
        },
        {
            "key": "ranking_contract_passes",
            "label": "Ranking contract passes",
            "passed": benchmark_row.get("ranking_contract_status") == "pass",
            "detail": str(benchmark_row.get("ranking_contract_status", "unknown")),
        },
        {
            "key": "comparator_is_measured_volatiles",
            "label": "Comparator signal is wet-lab measured_volatiles",
            "passed": benchmark_row.get("reference_signal_origin") == "measured_volatiles",
            "detail": str(benchmark_row.get("reference_signal_origin", "unknown")),
        },
        {
            "key": "external_quantitative_origin",
            "label": "Source is externally quantitative",
            "passed": evidence_row.external_data_status == "external_quantitative",
            "detail": evidence_row.external_data_status,
        },
        {
            "key": "minimum_quantitative_closed_targets",
            "label": "At least two compounds are quantitatively closed",
            "passed": int(support_counts.get("quantitative_closed", 0)) >= 2,
            "detail": str(int(support_counts.get("quantitative_closed", 0))),
        },
        {
            "key": "no_internal_or_directional_dependencies",
            "label": "No internal-candidate or directional dependencies remain",
            "passed": int(support_counts.get("internal_candidate", 0)) == 0 and int(support_counts.get("directional_support", 0)) == 0,
            "detail": f"internal={int(support_counts.get('internal_candidate', 0))}; directional={int(support_counts.get('directional_support', 0))}",
        },
    ]
    return requirement_rows


def _select_matrix_promotion_target(
    benchmark_rows: Iterable[Mapping[str, Any]],
    evidence_rows: Iterable[MatrixBenchmarkEvidence],
) -> Optional[Dict[str, Any]]:
    benchmark_list = list(benchmark_rows)
    evidence_list = list(evidence_rows)
    candidates = [
        row for row in benchmark_list
        if row.get("target_profile") in {"mixed", "meaty_positive"}
        and row.get("ranking_contract_status") == "pass"
    ]
    if not candidates:
        return None

    external_anchor_counts: Dict[str, int] = defaultdict(int)
    distinct_external_states: Dict[str, set[str]] = defaultdict(set)
    for row in evidence_list:
        if row.external_data_status != "external_quantitative":
            continue
        external_anchor_counts[row.protein_type] += 1
        if row.process_state:
            distinct_external_states[row.protein_type].add(str(row.process_state))

    def rank_tuple(row: Mapping[str, Any]) -> tuple[int, int, int, int, str, str]:
        protein_type = str(row.get("protein_type", "free"))
        counts = row.get("support_counts", {})
        return (
            int(external_anchor_counts.get(protein_type, 0)),
            len(distinct_external_states.get(protein_type, set())),
            int(counts.get("quantitative_closed", 0)),
            -1 * (int(counts.get("internal_candidate", 0)) + int(counts.get("directional_support", 0)) + int(counts.get("open_gap", 0))),
            "1" if row.get("target_profile") == "mixed" else "0",
            str(row.get("benchmark_id", "unknown")),
        )

    selected = sorted(candidates, key=rank_tuple, reverse=True)[0]
    protein_type = str(selected.get("protein_type", "free"))
    rationale = []
    rationale.append(f"same_protein_external_anchor_count={int(external_anchor_counts.get(protein_type, 0))}")
    rationale.append(f"same_protein_external_process_states={len(distinct_external_states.get(protein_type, set()))}")
    rationale.append(f"quantitative_closed={int(selected.get('support_counts', {}).get('quantitative_closed', 0))}")
    rationale.append(f"internal_candidate={int(selected.get('support_counts', {}).get('internal_candidate', 0))}")
    rationale.append(f"internal_measured_candidate={int(selected.get('support_counts', {}).get('internal_measured_candidate', 0))}")
    rationale.append(f"internal_reference_candidate={int(selected.get('support_counts', {}).get('internal_reference_candidate', 0))}")
    return {
        "benchmark_id": selected.get("benchmark_id"),
        "protein_type": protein_type,
        "process_state": selected.get("process_state"),
        "target_profile": selected.get("target_profile"),
        "rationale": rationale,
        "selection_policy": "prefer_mixed_matrix_lanes_with_the_broadest_same_protein_external_anchor_span_before_mechanistic_escalation",
    }


def build_matrix_promotion_contract_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    status_payload = build_matrix_target_status_artifact(benchmark_files, target_tag=target_tag)
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    evidence_lookup = {row.benchmark_id: row for row in evidence_rows}

    benchmark_assessments: List[Dict[str, Any]] = []
    for row in status_payload.get("benchmarks", []):
        benchmark_id = str(row.get("benchmark_id", "unknown"))
        evidence = evidence_lookup.get(benchmark_id)
        if evidence is None:
            continue
        requirements = _matrix_promotion_requirement_rows(row, evidence)
        benchmark_assessments.append(
            {
                "benchmark_id": benchmark_id,
                "protein_type": row.get("protein_type"),
                "process_state": row.get("process_state"),
                "target_profile": row.get("target_profile"),
                "promotion_ready": bool(row.get("promotion_ready", False)),
                "promotion_blocker": row.get("promotion_blocker", "unknown"),
                "requirements": requirements,
            }
        )

    selected_target = _select_matrix_promotion_target(status_payload.get("benchmarks", []), evidence_rows)
    summary = status_payload.get("summary", {})
    return {
        "schema_version": "1.0",
        "description": "Explicit matrix promotion contract defining how a matrix benchmark moves from internal candidate support to external decision readiness.",
        "promotion_rule": {
            "contract_id": "matrix_external_decision_ready_v1",
            "minimum_quantitative_closed_targets": 2,
            "disallow_internal_candidate_support": True,
            "disallow_directional_support": True,
            "requires_measured_volatiles": True,
            "requires_external_quantitative_origin": True,
            "requires_mixed_or_meaty_positive_target_profile": True,
            "requires_passing_ranking_contract": True,
            "notes": [
                "External decision readiness is a benchmark-level promotion state, not a generic matrix-family claim.",
                "Internal reproducibility candidates and transferred priors can strengthen triage, but they do not by themselves unlock promotion.",
                "Mechanistic refinement stays secondary until the observable audit says the remaining blocker is no longer external evidence or transfer dependence.",
            ],
        },
        "selected_promotion_target": selected_target,
        "benchmarks": benchmark_assessments,
        "summary": {
            "benchmarks_assessed": len(benchmark_assessments),
            "promotion_ready": int(summary.get("promotion_ready", 0)),
            "mechanistic_priority_ready": int(summary.get("mechanistic_priority_ready", 0)),
        },
    }


def _matrix_closure_action(
    *,
    compound_row: Mapping[str, Any],
    benchmark_row: Mapping[str, Any],
) -> str:
    support_status = str(compound_row.get("support_status", "open_gap"))
    evidence_state = str(compound_row.get("evidence_state", "still_missing"))
    calibration_strength = str(compound_row.get("calibration_evidence_strength", "heuristic"))

    if support_status == "quantitative_closed":
        return "already_closed"
    if calibration_strength in {"literature_anchored", "conditional_literature_anchored"}:
        return "literature_anchor_available"
    if _is_internal_candidate_support_status(support_status) and calibration_strength == "heuristic":
        return "mechanistic_blocker"
    if (_is_internal_candidate_support_status(support_status) or support_status == "directional_support") and (
        calibration_strength in {"class_anchored", "directional_transferred"}
        or evidence_state in {"externally_benchmarked", "transferred_prior", "safety_reference"}
    ):
        return "class_level_transfer_acceptable"
    if benchmark_row.get("mechanistic_priority_ready") and _is_internal_candidate_support_status(support_status):
        return "mechanistic_blocker"
    return "external_data_blocker"


def _mechanistic_refinement_expected_change(compound_rows: Iterable[Mapping[str, Any]]) -> str:
    rows = list(compound_rows)
    roles = {str(row.get("role", "unknown")) for row in rows}
    if roles == {"adverse_marker"}:
        return "clarify whether the lane remains meaty-positive once named adverse-marker closure is resolved"
    if roles == {"desirable_marker"}:
        return "clarify whether named desirable-marker closure can move the lane toward external decision readiness"
    return "clarify whether named mechanistic blockers materially change benchmark readiness before broader retuning"


def build_matrix_observable_closure_audit(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    status_payload = build_matrix_target_status_artifact(benchmark_files, target_tag=target_tag)
    evidence_rows = build_matrix_benchmark_evidence_audit(benchmark_files)
    evidence_lookup = {row.benchmark_id: row for row in evidence_rows}
    selected_target = _select_matrix_promotion_target(status_payload.get("benchmarks", []), evidence_rows)

    audit_rows: List[Dict[str, Any]] = []
    action_counts: Dict[str, int] = defaultdict(int)
    mechanistic_watchlist: List[Dict[str, Any]] = []
    for benchmark_row in status_payload.get("benchmarks", []):
        if benchmark_row.get("target_profile") not in {"mixed", "meaty_positive"}:
            continue
        benchmark_id = str(benchmark_row.get("benchmark_id", "unknown"))
        evidence = evidence_lookup.get(benchmark_id)
        compound_rows: List[Dict[str, Any]] = []
        benchmark_action_counts: Dict[str, int] = defaultdict(int)
        for compound_row in benchmark_row.get("compounds", []):
            closure_action = _matrix_closure_action(compound_row=compound_row, benchmark_row=benchmark_row)
            benchmark_action_counts[closure_action] += 1
            action_counts[closure_action] += 1
            compound_rows.append(
                {
                    **compound_row,
                    "closure_action": closure_action,
                }
            )
        audit_rows.append(
            {
                "benchmark_id": benchmark_id,
                "protein_type": benchmark_row.get("protein_type"),
                "process_state": benchmark_row.get("process_state"),
                "target_profile": benchmark_row.get("target_profile"),
                "promotion_blocker": benchmark_row.get("promotion_blocker"),
                "source_origin": evidence.source_origin if evidence is not None else "unknown",
                "compounds": compound_rows,
                "closure_action_counts": dict(sorted(benchmark_action_counts.items())),
            }
        )

        if benchmark_row.get("mechanistic_priority_ready"):
            blocker_rows = [
                row for row in compound_rows
                if str(row.get("closure_action", "unknown")) == "mechanistic_blocker"
            ]
            if blocker_rows:
                mechanistic_watchlist.append(
                    {
                        "benchmark_id": benchmark_id,
                        "protein_type": benchmark_row.get("protein_type"),
                        "process_state": benchmark_row.get("process_state"),
                        "promotion_blocker": benchmark_row.get("promotion_blocker"),
                        "target_compounds": [str(row.get("compound", "unknown")) for row in blocker_rows],
                        "target_roles": sorted({str(row.get("role", "unknown")) for row in blocker_rows}),
                        "expected_decision_change": _mechanistic_refinement_expected_change(blocker_rows),
                        "allowed_scope": "named_compound_refinement_only",
                        "offline_compute_gate": "escalate_only_if_named_compounds_change_benchmark_visible_decision_readiness",
                    }
                )

    return {
        "schema_version": "1.0",
        "description": "Compound-level observable closure audit for mixed matrix lanes, labeling the next closure action needed for each decision-driving compound.",
        "selected_promotion_target": selected_target,
        "benchmarks": audit_rows,
        "mechanistic_refinement_watchlist": mechanistic_watchlist,
        "summary": {
            "benchmarks_audited": len(audit_rows),
            "mechanistic_watchlist_count": len(mechanistic_watchlist),
            "closure_action_counts": dict(sorted(action_counts.items())),
        },
    }





def compare_matrix_benchmark_delta_sets(
    current_rows: Iterable[MatrixBenchmarkDelta],
    baseline_rows: Iterable[MatrixBenchmarkDelta],
    *,
    current_evidence: Optional[Iterable[MatrixBenchmarkEvidence]] = None,
    baseline_evidence: Optional[Iterable[MatrixBenchmarkEvidence]] = None,
) -> List[MatrixBenchmarkBranchDelta]:
    current_lookup = {(row.benchmark_id, _normalize_name(row.compound)): row for row in current_rows}
    baseline_lookup = {(row.benchmark_id, _normalize_name(row.compound)): row for row in baseline_rows}
    current_evidence_lookup = {row.benchmark_id: row for row in (current_evidence or [])}
    baseline_evidence_lookup = {row.benchmark_id: row for row in (baseline_evidence or [])}

    rows: List[MatrixBenchmarkBranchDelta] = []
    for key in sorted(set(current_lookup) | set(baseline_lookup)):
        current = current_lookup.get(key)
        baseline = baseline_lookup.get(key)
        current_meta = current_evidence_lookup.get(key[0])
        baseline_meta = baseline_evidence_lookup.get(key[0])

        if current is None:
            change_type = "removed"
        elif baseline is None:
            change_type = "added"
        else:
            predicted_changed = not math.isclose(current.predicted_ppb, baseline.predicted_ppb, rel_tol=1.0e-9, abs_tol=1.0e-12)
            ratio_changed = not math.isclose(current.ratio, baseline.ratio, rel_tol=1.0e-9, abs_tol=1.0e-12)
            metadata_changed = (
                current.execution_path != baseline.execution_path
                or current.reference_signal_origin != baseline.reference_signal_origin
                or (current_meta.source_origin if current_meta else "n/a") != (baseline_meta.source_origin if baseline_meta else "n/a")
                or (current_meta.external_data_status if current_meta else "n/a") != (baseline_meta.external_data_status if baseline_meta else "n/a")
            )
            if predicted_changed or ratio_changed:
                change_type = "modified"
            elif metadata_changed:
                change_type = "metadata_changed"
            else:
                continue

        predicted_delta = None if current is None or baseline is None else current.predicted_ppb - baseline.predicted_ppb
        ratio_delta = None if current is None or baseline is None else current.ratio - baseline.ratio
        rows.append(
            MatrixBenchmarkBranchDelta(
                benchmark_id=key[0],
                compound=current.compound if current is not None else baseline.compound,
                change_type=change_type,
                current_present=current is not None,
                baseline_present=baseline is not None,
                current_execution_path=current.execution_path if current is not None else "n/a",
                baseline_execution_path=baseline.execution_path if baseline is not None else "n/a",
                current_reference_signal_origin=current.reference_signal_origin if current is not None else "n/a",
                baseline_reference_signal_origin=baseline.reference_signal_origin if baseline is not None else "n/a",
                current_source_origin=current_meta.source_origin if current_meta is not None else "n/a",
                baseline_source_origin=baseline_meta.source_origin if baseline_meta is not None else "n/a",
                current_external_data_status=current_meta.external_data_status if current_meta is not None else "n/a",
                baseline_external_data_status=baseline_meta.external_data_status if baseline_meta is not None else "n/a",
                current_predicted_ppb=current.predicted_ppb if current is not None else None,
                baseline_predicted_ppb=baseline.predicted_ppb if baseline is not None else None,
                predicted_delta_ppb=predicted_delta,
                current_ratio=current.ratio if current is not None else None,
                baseline_ratio=baseline.ratio if baseline is not None else None,
                ratio_delta=ratio_delta,
            )
        )
    return rows


def render_matrix_branch_deltas_markdown(
    rows: Iterable[MatrixBenchmarkBranchDelta],
    *,
    base_ref: str,
) -> str:
    delta_rows = list(rows)
    lines = [
        f"# Matrix Benchmark Branch Comparison vs {base_ref}",
        "",
        "| Benchmark | Compound | Change | Current Path | Base Path | Current Origin | Base Origin | Current Data Status | Base Data Status | Current Predicted ppb | Base Predicted ppb | Δ Predicted ppb |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: |",
    ]
    for row in delta_rows:
        current_predicted = f"{row.current_predicted_ppb:.3f}" if row.current_predicted_ppb is not None else "n/a"
        baseline_predicted = f"{row.baseline_predicted_ppb:.3f}" if row.baseline_predicted_ppb is not None else "n/a"
        predicted_delta = f"{row.predicted_delta_ppb:.3f}" if row.predicted_delta_ppb is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.compound} | {row.change_type} | {row.current_execution_path} | {row.baseline_execution_path} | {row.current_source_origin} | {row.baseline_source_origin} | {row.current_external_data_status} | {row.baseline_external_data_status} | {current_predicted} | {baseline_predicted} | {predicted_delta} |"
        )
    lines.extend([
        "",
        f"Changed rows: {len(delta_rows)}",
        f"Added rows: {sum(1 for row in delta_rows if row.change_type == 'added')}",
        f"Removed rows: {sum(1 for row in delta_rows if row.change_type == 'removed')}",
        f"Modified rows: {sum(1 for row in delta_rows if row.change_type == 'modified')}",
        f"Metadata-only changes: {sum(1 for row in delta_rows if row.change_type == 'metadata_changed')}",
    ])
    if delta_rows and not any(row.baseline_present for row in delta_rows):
        lines.extend([
            f"Base ref {base_ref} exposes no matrix benchmark rows under the current evaluator.",
            "All reported rows are additions from the current branch/worktree.",
        ])
    return "\n".join(lines) + "\n"





def summarize_benchmarks(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> List[BenchmarkSummary]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    summaries: List[BenchmarkSummary] = []
    for bench_file in bench_files:
        bench_path = Path(bench_file)
        bench = load_benchmark(bench_path)
        evaluation = evaluate_benchmark(bench_path, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation,
            protein_type=bench.get("protein_type", "free"),
            thresholds=thresholds,
        )
        if not evaluation.supported:
            metadata = get_benchmark_metadata(bench)
            summary = BenchmarkSummary(
                benchmark_id=summary.benchmark_id,
                bench_file=summary.bench_file,
                tier=metadata.tier,
                family=metadata.family,
                execution_path=metadata.execution_path,
                benchmark_engine=metadata.benchmark_engine,
                comparator_signal=metadata.comparator_signal,
                cantera_role=metadata.cantera_role,
                target_snapshot_policy=metadata.target_snapshot_policy,
                thermodynamic_gating_policy=metadata.thermodynamic_gating_policy,
                supported=summary.supported,
                reason=summary.reason,
                protein_type=summary.protein_type,
                coverage=summary.coverage,
                matched_compounds=summary.matched_compounds,
                total_compounds=summary.total_compounds,
                pearson_r=summary.pearson_r,
                mae_ppb=summary.mae_ppb,
                max_ratio=summary.max_ratio,
                mean_ratio=summary.mean_ratio,
                ranking_status=summary.ranking_status,
                scale_status=summary.scale_status,
                overall_status=summary.overall_status,
                strict_ready=summary.strict_ready,
                blocking_issues=summary.blocking_issues,
                conditions=bench.get("conditions", {}),
                # 2026-08-27 (Wave I): an unsupported benchmark still has an evidence role,
                # and losing it here would make an unsupported fit-recovery row read as
                # "predictive" in the summary.
                evidence_role=benchmark_evidence_role(summary.benchmark_id, bench_path),
            )
        summaries.append(_enrich_benchmark_summary_family_metadata(summary, bench))
    return summaries


def _payload_role_from_evidence_state(evidence_state: str) -> str:
    normalized = str(evidence_state).strip().lower()
    if normalized in {"externally_benchmarked", "internally_benchmarked"}:
        return "benchmark_payload"
    if normalized == "conditional_calibration":
        return "calibration_payload"
    if normalized in {"transferred_prior", "safety_reference"}:
        return "directional_prior"
    return "structural_gap_extrapolation"


def _benchmark_compound_names(bench: dict, summary: BenchmarkSummary) -> List[str]:
    names: List[str] = []
    measured = bench.get("measured_volatiles", {}) or {}
    if isinstance(measured, dict):
        names.extend(str(name) for name in measured.keys())
    reference = bench.get("reference_volatiles", {}) or {}
    if isinstance(reference, dict):
        names.extend(str(name) for name in reference.keys())
    names.extend(str(name) for name in summary.ranked_observable_targets)
    names.extend(str(name) for name in summary.adverse_markers)
    deduped: List[str] = []
    seen: set[str] = set()
    for name in names:
        normalized = _normalize_name(name)
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(str(name))
    return deduped


def _enrich_benchmark_summary_family_metadata(summary: BenchmarkSummary, bench: dict) -> BenchmarkSummary:
    chemistry_families: List[str] = []
    slr_families: List[str] = []
    family_lane_names: List[str] = []
    payload_roles: List[str] = []

    for compound_name in _benchmark_compound_names(bench, summary):
        panel_entry = get_compound_panel_entry(compound_name) or {}
        chemistry_family = str(panel_entry.get("chemistry_family", "")).strip()
        if not chemistry_family:
            continue
        descriptor = resolve_family_descriptor(chemistry_family)
        chemistry_families.append(chemistry_family)
        if descriptor:
            slr_families.append(str(descriptor.get("slr_family", "")).zfill(2))
            family_lane_names.append(str(descriptor.get("display_name", chemistry_family)))
        payload_roles.append(_payload_role_from_evidence_state(str(panel_entry.get("evidence_state", "still_missing"))))

    mapped_family = _RUNTIME_BENCHMARK_FAMILY_MAP.get(summary.benchmark_id, "")
    if mapped_family:
        descriptor = resolve_family_descriptor(mapped_family)
        chemistry_families.append(mapped_family)
        if descriptor:
            slr_families.append(str(descriptor.get("slr_family", "")).zfill(2))
            family_lane_names.append(str(descriptor.get("display_name", mapped_family)))
        payload_roles.extend(_RUNTIME_BENCHMARK_PAYLOAD_ROLES.get(summary.benchmark_id, ["benchmark_intake"]))

    return BenchmarkSummary(
        benchmark_id=summary.benchmark_id,
        bench_file=summary.bench_file,
        tier=summary.tier,
        family=summary.family,
        execution_path=summary.execution_path,
        benchmark_engine=summary.benchmark_engine,
        comparator_signal=summary.comparator_signal,
        cantera_role=summary.cantera_role,
        target_snapshot_policy=summary.target_snapshot_policy,
        thermodynamic_gating_policy=summary.thermodynamic_gating_policy,
        supported=summary.supported,
        reason=summary.reason,
        protein_type=summary.protein_type,
        coverage=summary.coverage,
        matched_compounds=summary.matched_compounds,
        total_compounds=summary.total_compounds,
        pearson_r=summary.pearson_r,
        mae_ppb=summary.mae_ppb,
        max_ratio=summary.max_ratio,
        mean_ratio=summary.mean_ratio,
        ranking_status=summary.ranking_status,
        scale_status=summary.scale_status,
        overall_status=summary.overall_status,
        strict_ready=summary.strict_ready,
        blocking_issues=list(summary.blocking_issues),
        conditions=dict(summary.conditions),
        process_state=summary.process_state,
        ranked_observable_targets=list(summary.ranked_observable_targets),
        adverse_markers=list(summary.adverse_markers),
        ranking_contract_status=summary.ranking_contract_status,
        calibration_mode=summary.calibration_mode,
        reference_signal_origin=summary.reference_signal_origin,
        mean_abs_log10_error=summary.mean_abs_log10_error,
        chemistry_families=sorted(dict.fromkeys(chemistry_families)),
        slr_families=sorted(dict.fromkeys(slr_families)),
        payload_roles=sorted(dict.fromkeys(payload_roles)),
        family_lane_names=sorted(dict.fromkeys(family_lane_names)),
        # 2026-08-27 (Wave I): this rebuild-with-families constructor is the one the
        # summary generator actually returns, so the evidence role must be carried
        # through here or every row silently reports "predictive".
        evidence_role=summary.evidence_role,
    )


def build_family_lane_validation_artifact(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> Dict[str, Any]:
    summaries = summarize_benchmarks(benchmark_files, target_tag=target_tag, thresholds=thresholds)
    enriched_summaries: List[BenchmarkSummary] = []
    for summary in summaries:
        bench = load_benchmark(summary.bench_file)
        enriched_summaries.append(_enrich_benchmark_summary_family_metadata(summary, bench))

    family_rows: List[Dict[str, Any]] = []
    family_groups: Dict[str, Dict[str, Any]] = {}
    lane_groups: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
        "execution_path": "unknown",
        "benchmark_count": 0,
        "strict_ready_count": 0,
        "supported_count": 0,
        "status_counts": defaultdict(int),
        "chemistry_families": set(),
        "payload_roles": set(),
    })

    for summary in enriched_summaries:
        lane_bucket = lane_groups[summary.execution_path]
        lane_bucket["execution_path"] = summary.execution_path
        lane_bucket["benchmark_count"] += 1
        lane_bucket["strict_ready_count"] += int(bool(summary.strict_ready))
        lane_bucket["supported_count"] += int(bool(summary.supported))
        lane_bucket["status_counts"][summary.overall_status] += 1
        lane_bucket["chemistry_families"].update(summary.chemistry_families)
        lane_bucket["payload_roles"].update(summary.payload_roles)

        for chemistry_family in summary.chemistry_families:
            descriptor = resolve_family_descriptor(chemistry_family)
            family_bucket = family_groups.setdefault(
                chemistry_family,
                {
                    "chemistry_family": chemistry_family,
                    "slr_family": str(descriptor.get("slr_family", "")).zfill(2) if descriptor else "",
                    "display_name": str(descriptor.get("display_name", chemistry_family)) if descriptor else chemistry_family,
                    "strategic_posture": str(descriptor.get("strategic_posture", "unknown")) if descriptor else "unknown",
                    "benchmark_count": 0,
                    "strict_ready_count": 0,
                    "supported_count": 0,
                    "status_counts": defaultdict(int),
                    "execution_paths": defaultdict(int),
                    "payload_roles": set(),
                    "benchmark_ids": [],
                },
            )
            family_bucket["benchmark_count"] += 1
            family_bucket["strict_ready_count"] += int(bool(summary.strict_ready))
            family_bucket["supported_count"] += int(bool(summary.supported))
            family_bucket["status_counts"][summary.overall_status] += 1
            family_bucket["execution_paths"][summary.execution_path] += 1
            family_bucket["payload_roles"].update(summary.payload_roles)
            family_bucket["benchmark_ids"].append(summary.benchmark_id)

    for chemistry_family, row in sorted(family_groups.items(), key=lambda item: (item[1]["slr_family"], item[0])):
        family_rows.append(
            {
                "chemistry_family": chemistry_family,
                "slr_family": row["slr_family"],
                "display_name": row["display_name"],
                "strategic_posture": row["strategic_posture"],
                "benchmark_count": int(row["benchmark_count"]),
                "strict_ready_count": int(row["strict_ready_count"]),
                "supported_count": int(row["supported_count"]),
                "status_counts": dict(sorted(row["status_counts"].items())),
                "execution_paths": dict(sorted(row["execution_paths"].items())),
                "payload_roles": sorted(row["payload_roles"]),
                "benchmark_ids": sorted(dict.fromkeys(row["benchmark_ids"])),
            }
        )

    lane_rows: List[Dict[str, Any]] = []
    for execution_path, row in sorted(lane_groups.items()):
        lane_rows.append(
            {
                "execution_path": execution_path,
                "benchmark_count": int(row["benchmark_count"]),
                "strict_ready_count": int(row["strict_ready_count"]),
                "supported_count": int(row["supported_count"]),
                "status_counts": dict(sorted(row["status_counts"].items())),
                "chemistry_families": sorted(row["chemistry_families"]),
                "payload_roles": sorted(row["payload_roles"]),
            }
        )

    return {
        "summary": {
            "benchmark_count": len(enriched_summaries),
            "family_count": len(family_rows),
            "lane_count": len(lane_rows),
            "tracked_execution_paths": [row["execution_path"] for row in lane_rows],
        },
        "families": family_rows,
        "lanes": lane_rows,
        "benchmarks": [
            {
                "benchmark_id": summary.benchmark_id,
                "execution_path": summary.execution_path,
                "overall_status": summary.overall_status,
                "strict_ready": bool(summary.strict_ready),
                "chemistry_families": list(summary.chemistry_families),
                "slr_families": list(summary.slr_families),
                "payload_roles": list(summary.payload_roles),
                "family_lane_names": list(summary.family_lane_names),
            }
            for summary in enriched_summaries
        ],
    }


def render_family_lane_validation_markdown(payload: Dict[str, Any]) -> str:
    lines = [
        "# Family Lane Validation",
        "",
        "| SLR | Family | Posture | Benchmarks | Strict Ready | Supported | Payload Roles | Execution Paths |",
        "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in payload.get("families", []):
        execution_paths = ", ".join(f"{key}={value}" for key, value in row.get("execution_paths", {}).items()) or "none"
        lines.append(
            f"| {row.get('slr_family', '') or 'n/a'} | {row.get('chemistry_family', 'unknown')} | {row.get('strategic_posture', 'unknown')} | {int(row.get('benchmark_count', 0))} | {int(row.get('strict_ready_count', 0))} | {int(row.get('supported_count', 0))} | {', '.join(str(item) for item in row.get('payload_roles', [])) or 'none'} | {execution_paths} |"
        )
    lines.extend([
        "",
        "## Lane Summary",
        "",
        "| Execution Path | Benchmarks | Strict Ready | Supported | Chemistry Families | Payload Roles |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ])
    for row in payload.get("lanes", []):
        lines.append(
            f"| {row.get('execution_path', 'unknown')} | {int(row.get('benchmark_count', 0))} | {int(row.get('strict_ready_count', 0))} | {int(row.get('supported_count', 0))} | {', '.join(str(item) for item in row.get('chemistry_families', [])) or 'none'} | {', '.join(str(item) for item in row.get('payload_roles', [])) or 'none'} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Benchmarks summarized: {int(summary.get('benchmark_count', 0))}",
        f"Chemistry families summarized: {int(summary.get('family_count', 0))}",
        f"Execution lanes summarized: {int(summary.get('lane_count', 0))}",
    ])
    return "\n".join(lines) + "\n"


def _benchmark_status_score(summary: BenchmarkSummary) -> int:
    ranking = {
        "unsupported": 0,
        "coverage-gap": 1,
        "ranking-gap": 2,
        "scale-gap": 3,
        "partial-pass": 4,
        "pass-no-ranking": 5,
        "pass": 5,
    }
    return ranking.get(summary.overall_status, 0)


def thermodynamic_gating_materially_improves(
    baseline: BenchmarkSummary,
    gated: BenchmarkSummary,
) -> bool:
    if not baseline.supported or not gated.supported:
        return False
    if gated.coverage + 1.0e-12 < baseline.coverage:
        return False
    if _benchmark_status_score(gated) < _benchmark_status_score(baseline):
        return False

    baseline_mae = baseline.mae_ppb if baseline.mae_ppb is not None else float("inf")
    gated_mae = gated.mae_ppb if gated.mae_ppb is not None else float("inf")
    baseline_ratio = baseline.max_ratio if baseline.max_ratio is not None else float("inf")
    gated_ratio = gated.max_ratio if gated.max_ratio is not None else float("inf")

    mae_improvement = baseline_mae - gated_mae
    ratio_improvement = baseline_ratio - gated_ratio
    mae_threshold = max(
        THERMODYNAMIC_GATING_MIN_ABSOLUTE_MAE_IMPROVEMENT_PPB,
        baseline_mae * THERMODYNAMIC_GATING_MIN_RELATIVE_MAE_IMPROVEMENT if math.isfinite(baseline_mae) else float("inf"),
    )

    return (
        (math.isfinite(mae_improvement) and mae_improvement >= mae_threshold)
        or (math.isfinite(ratio_improvement) and ratio_improvement >= THERMODYNAMIC_GATING_MIN_RATIO_IMPROVEMENT)
    )


def audit_thermodynamic_gating(
    bench_file: Path | str,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> ThermodynamicGatingAudit:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
    metadata = get_benchmark_metadata(bench)
    protein_type = bench.get("protein_type", "free")

    if metadata.execution_path != "free_precursor" or metadata.family == "safety":
        return ThermodynamicGatingAudit(
            benchmark_id=bench["benchmark_id"],
            bench_file=bench_path,
            execution_path=metadata.execution_path,
            applicable=False,
            baseline_overall_status="n/a",
            gated_overall_status="n/a",
            baseline_mae_ppb=None,
            gated_mae_ppb=None,
            baseline_max_ratio=None,
            gated_max_ratio=None,
            delta_mae_ppb=None,
            delta_max_ratio=None,
            material_improvement=False,
            recommended_policy="diagnostic_only",
            notes="Thermodynamic gating audit is currently only meaningful for non-safety free-precursor FAST benchmarks.",
        )

    baseline_eval = evaluate_benchmark(bench_path, target_tag=target_tag, thermodynamic_gating="off")
    gated_eval = evaluate_benchmark(bench_path, target_tag=target_tag, thermodynamic_gating="on")
    baseline_summary = summarize_evaluation(baseline_eval, protein_type=protein_type, thresholds=thresholds)
    gated_summary = summarize_evaluation(gated_eval, protein_type=protein_type, thresholds=thresholds)

    baseline_mae = baseline_summary.mae_ppb
    gated_mae = gated_summary.mae_ppb
    baseline_ratio = baseline_summary.max_ratio
    gated_ratio = gated_summary.max_ratio

    delta_mae = None if baseline_mae is None or gated_mae is None else baseline_mae - gated_mae
    delta_ratio = None if baseline_ratio is None or gated_ratio is None else baseline_ratio - gated_ratio
    material = thermodynamic_gating_materially_improves(baseline_summary, gated_summary)
    if material:
        notes = "Thermodynamic gating materially improves benchmark error without degrading supported coverage/status."
        policy = "benchmark_facing_candidate"
    else:
        notes = "Thermodynamic gating does not materially improve benchmark error under the current threshold and remains diagnostic-only."
        policy = "diagnostic_only"

    return ThermodynamicGatingAudit(
        benchmark_id=bench["benchmark_id"],
        bench_file=bench_path,
        execution_path=metadata.execution_path,
        applicable=True,
        baseline_overall_status=baseline_summary.overall_status,
        gated_overall_status=gated_summary.overall_status,
        baseline_mae_ppb=baseline_mae,
        gated_mae_ppb=gated_mae,
        baseline_max_ratio=baseline_ratio,
        gated_max_ratio=gated_ratio,
        delta_mae_ppb=delta_mae,
        delta_max_ratio=delta_ratio,
        material_improvement=material,
        recommended_policy=policy,
        notes=notes,
    )


def audit_all_thermodynamic_gating(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> List[ThermodynamicGatingAudit]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    return [
        audit_thermodynamic_gating(bench_file, target_tag=target_tag, thresholds=thresholds)
        for bench_file in bench_files
    ]



def snapshot_benchmark_targets(
    bench_file: Path | str,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkTargetSnapshot]:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
    if is_matrix_only_benchmark(bench):
        # Matrix-only benchmarks are benchmarked as intake/headspace checks,
        # not as FAST target-ranking snapshots.
        return []
    formulation = benchmark_to_formulation(bench)
    supported, _reason = _is_supported_formulation(formulation)
    if not supported:
        return []
    rec_result = _run_benchmark_recommendation(bench, target_tag=target_tag)

    snapshots: List[BenchmarkTargetSnapshot] = []
    for target in sorted(rec_result.get("targets", []), key=lambda row: (row.get("type", ""), -float(row.get("concentration", 0.0)), row.get("name", ""))):
        snapshots.append(
            BenchmarkTargetSnapshot(
                benchmark_id=bench["benchmark_id"],
                bench_file=bench_path,
                target_name=str(target.get("name", "")),
                target_type=str(target.get("type", "unknown")),
                roles=list(target.get("roles", [target.get("type", "unknown")])),
                target_class=str(target.get("projection", {}).get("target_class", "unknown")),
                evidence_state=str(target.get("projection", {}).get("evidence_state", "still_missing")),
                predicted_ppb=float(target.get("concentration", 0.0)),
                proxy_ppb=float(target.get("proxy_concentration", target.get("concentration", 0.0))),
                observable_ratio=float(target.get("projection", {}).get("proxy_to_observable_ratio", 1.0)),
                weighted_flux=float(target.get("weighted_flux", 0.0)),
                span=float(target.get("span", math.inf)),
                depth=int(target.get("depth", 0)),
                volatile_class=str(target.get("projection", {}).get("volatile_class", "other")),
                matrix_factor=float(target.get("projection", {}).get("matrix_factor", 1.0)),
                headspace_factor=float(target.get("projection", {}).get("headspace_factor", 1.0)),
                headspace_observable=bool(target.get("headspace_observable", True)),
                headspace_class=str(target.get("headspace_class", "assumed_observable")),
                henry_kaw_25c=float(target["henry_kaw_25c"]) if target.get("henry_kaw_25c") is not None else None,
                henry_source_name=target.get("henry_source_name"),
            )
        )
    return snapshots


def snapshot_all_benchmark_targets(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkTargetSnapshot]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    snapshots: List[BenchmarkTargetSnapshot] = []
    for bench_file in bench_files:
        snapshots.extend(snapshot_benchmark_targets(bench_file, target_tag=target_tag))
    return snapshots


def build_benchmark_index(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkIndexEntry]:
    summaries = summarize_benchmarks(benchmark_files=benchmark_files, target_tag=target_tag)
    return [
        BenchmarkIndexEntry(
            benchmark_id=summary.benchmark_id,
            bench_file=summary.bench_file,
            tier=summary.tier,
            family=summary.family,
            protein_type=summary.protein_type,
            execution_path=summary.execution_path,
            benchmark_engine=summary.benchmark_engine,
            cantera_role=summary.cantera_role,
            thermodynamic_gating_policy=summary.thermodynamic_gating_policy,
            supported=summary.supported,
            reason=summary.reason,
            status=summary.overall_status,
            strict_ready=summary.strict_ready,
            process_state=summary.process_state,
            ranking_contract_status=summary.ranking_contract_status,
        )
        for summary in summaries
    ]