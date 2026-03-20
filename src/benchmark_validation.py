from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from src.conditions import ReactionConditions
from src.barrier_constants import effective_barrier_from_rate_constant
from src.headspace import HeadspaceModel
from src.kinetics import KineticsEngine
from src.lipid_oxidation import PEA_LIPID_PROFILE, SOY_LIPID_PROFILE, predict_hexanal_generation
from src.matrix_calibration_registry import (
    describe_matrix_calibration,
    determine_matrix_process_state,
)
from src.pipeline import MaillardPipeline
from src.precursor_resolver import resolve_many
from src.projection_metadata import ProjectionMetadataMap, make_projection_metadata_row
from src.smirks_engine import SmirksEngine
from src.validation_contract import BenchmarkThresholds, DEFAULT_VALIDATION_CONTRACT
from src.safety import predict_acrylamide
from src.projection_utils import build_projection_rows
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

BENCHMARK_NAME_ALIASES = {
    "2methyl3furanthiol": {
        "2methyl3furanthiolmft",
        "2methyl3furylthiol",
        "2methylfuran3thiol",
        "mft",
    },
    "2furfurylthiol": {
        "2furfurylthiolfft",
        "2furylmethanethiol",
        "fft",
    },
    "bis2methyl3furyldisulfide": {
        "bis2methyl3furyldisulfide",
        "2methyl3furyl2methyl3furyldisulfide",
    },
    "pyrazine": {
        "25dimethylpyrazine",
        "23dimethylpyrazine",
        "26dimethylpyrazine",
        "2ethyl35dimethylpyrazine",
        "2ethylpyrazine",
        "trimethylpyrazine",
        "tetramethylpyrazine",
        "dimethylpyrazine",
    },
    "3methylbutanal": {
        "3methylbutanal",
        "isovaleraldehyde",
    },
    "2methylbutanal": {
        "2methylbutanal",
        "2methylbutyraldehyde",
    },
    "methional": {
        "3methylthiopropanal",
        "methional",
    },
    "acrylamide": {
        "acrylamide",
        "2-propenamide",
    },
}

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
    "Hexanal": 0.205,
    "2-Pentylfuran": 0.502,
    "1-Hexanol": 0.063,
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
    bench_path = Path(bench_file)
    with open(bench_path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def benchmark_to_conditions(bench: dict) -> ReactionConditions:
    return ReactionConditions(
        pH=bench["conditions"]["ph"],
        temperature_celsius=bench["conditions"]["temp_C"],
        water_activity=bench["conditions"]["water_activity"],
        protein_type=bench.get("protein_type", "free"),
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
        if any(token in name_lower for token in ["ribose", "glucose", "fructose", "xylose", "maltose", "sugar"]):
            sugars.append(name)
        elif any(token in name_lower for token in ["hexanal", "nonanal", "lipid", "fat", "furan"]):
            lipids.append(name)
        else:
            amino_acids.append(name)

    formulation = {
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
    return formulation


def _normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def _tokenize_name(name: str) -> List[str]:
    return [token for token in re.split(r"[^a-z0-9]+", name.lower()) if token]


def _is_supported_formulation(formulation: dict) -> tuple[bool, Optional[str]]:
    protein_type = str(formulation.get("protein_type", "free"))
    if protein_type != "free" and formulation.get("_skipped_matrix_precursors"):
        if protein_type in MATRIX_BENCHMARK_PROFILES:
            return True, None
        return False, f"No executable matrix-only benchmark path for protein_type={protein_type}"

    candidate_precursors = formulation["sugars"] + formulation["amino_acids"] + formulation["lipids"]
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
    target_aliases = BENCHMARK_NAME_ALIASES.get(target_norm, set())
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

    names = formulation["sugars"] + formulation["amino_acids"] + formulation.get("additives", []) + formulation.get("lipids", [])
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
        k = conditions.get_rate_constant(step.reaction_family or "unknown", ea_override_kcal=barrier_for_rate)
        bar_eff = effective_barrier_from_rate_constant(
            k,
            conditions.temperature_kelvin,
            step.reaction_family or "unknown",
        )
        rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)

    from src.recommend import Recommender, _canon
    from src.lipid_oxidation import predict_hexanal_generation

    initial_concentrations = {}
    ratios = formulation.get("molar_ratios", {})
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[_canon(precursor.smiles)] = qty

    protein_type = formulation.get("protein_type", "free")
    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is not None:
        oxidation = predict_hexanal_generation(
            model["lipid_profile"],
            temp_C=float(conditions.temperature_celsius),
            time_min=float(formulation.get("time_minutes", 60.0)),
            oxygen_availability=1.0,
        )
        oxidation_load_ppb = float(oxidation["total_hydroperoxide"]) * 1000.0
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
            smiles_map = {
                "Hexanal": "CCCCCC=O",
                "Nonanal": "CCCCCCCCC=O",
                "1-Hexanol": "CCCCCCO",
                "2-Pentylfuran": "CCCCC1=CC=CO1"
            }
            smi = smiles_map.get(compound_name)
            if smi:
                mol = Chem.MolFromSmiles(smi)
                mw = sum(atom.GetMass() for atom in mol.GetAtoms()) if mol else 100.0
                target_proxy_ppb = oxidation_load_ppb * float(yield_factor)
                # target_proxy_ppb = molar_concentration * mw * ppb_conversion_factor
                molar_conc = target_proxy_ppb / (mw * DEFAULT_PROJECTION_STRATEGY.ppb_conversion_factor)
                canon = _canon(smi)
                initial_concentrations[canon] = initial_concentrations.get(canon, 0.0) + molar_conc

    rec = Recommender()
    return rec.predict_from_steps(
        steps,
        heuristic_barriers,
        initial_concentrations,
        temperature_kelvin=conditions.temperature_kelvin,
        time_minutes=float(formulation.get("time_minutes", 60.0)),
        protein_type=protein_type,
        denaturation_state=float(bench.get("denaturation_state", 0.5)),
    )


def _run_matrix_only_benchmark_prediction(bench: dict) -> dict:
    protein_type = str(bench.get("protein_type", "free"))
    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is None:
        raise BenchmarkNotSupportedError(f"No matrix-only predictor for protein_type={protein_type}")

    conditions = bench["conditions"]
    process_state = str((bench.get("process_metadata") or {}).get(
        "state",
        determine_matrix_process_state(
            temperature_celsius=float(conditions["temp_C"]),
            time_minutes=float(conditions["time_min"]),
        ),
    ))
    oxidation = predict_hexanal_generation(
        model["lipid_profile"],
        temp_C=float(conditions["temp_C"]),
        time_min=float(conditions["time_min"]),
        oxygen_availability=1.0,
    )
    oxidation_load_ppb = float(oxidation["total_hydroperoxide"]) * 1000.0
    headspace_model = HeadspaceModel()
    pH = conditions.get("ph")

    predicted_ppb: Dict[str, float] = {}
    predicted_proxy_ppb: Dict[str, float] = {}
    projection_metadata: Dict[str, Dict[str, Any]] = {}
    for compound, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
        headspace_factor = headspace_model.get_matrix_benchmark_headspace_factor(
            compound,
            protein_type=protein_type,
            pH=pH,
            temperature_celsius=float(conditions["temp_C"]),
            time_minutes=float(conditions["time_min"]),
        )
        calibration = describe_matrix_calibration(
            compound,
            protein_type=protein_type,
            process_state=process_state,
        )
        calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
        release_factor = headspace_factor / calibration_factor if calibration_factor > 0.0 else headspace_factor
        proxy_ppb = oxidation_load_ppb * float(yield_factor) * release_factor
        observable_ppb = proxy_ppb * calibration_factor
        predicted_proxy_ppb[compound] = proxy_ppb
        predicted_ppb[compound] = observable_ppb
        projection_metadata[compound] = make_projection_metadata_row(
            compound=compound,
            proxy_ppb=proxy_ppb,
            observable_ppb=observable_ppb,
            extras={
                "matrix_factor": 1.0,
                "headspace_factor": release_factor,
                "total_observable_factor": headspace_factor,
                "process_state": process_state,
                **calibration,
            },
        )

    return {
        "targets": [],
        "metrics": {
            "matrix_model": protein_type,
            "oxidation_load_ppb": oxidation_load_ppb,
        },
        "predicted_ppb": predicted_ppb,
        "predicted_proxy_ppb": predicted_proxy_ppb,
        "projection_metadata": projection_metadata,
        "debug_paths": {},
        "species_names": {compound: compound for compound in predicted_ppb},
    }


def evaluate_benchmark(
    bench_file: Path | str,
    target_tag: str = DEFAULT_TARGET_TAG,
    thermodynamic_gating: str = "auto",
) -> BenchmarkEvaluation:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
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
        for name, data in bench["precursors"].items():
            n_low = name.lower()
            if "asparagine" in n_low or "asn" in n_low:
                asn_conc = data["concentration_mM"]
            if any(s in n_low for s in ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar"]):
                sugar_conc += data["concentration_mM"]
        
        # We assume 180 min if not specified, but Parker 2012 has it in conditions
        time_min = bench["conditions"].get("time_min", 20.0)
        safety_res = predict_acrylamide(
            asparagine_mM=asn_conc,
            reducing_sugar_mM=sugar_conc,
            temp_C=conditions.temperature_celsius,
            time_min=time_min,
            pH=conditions.pH,
        )
        rec_result = {"predicted_ppb": {"acrylamide": safety_res.acrylamide_ppb}, "projection_metadata": {}}
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


def summarize_evaluation(
    evaluation: BenchmarkEvaluation,
    *,
    protein_type: str = "free",
    thresholds: BenchmarkThresholds = DEFAULT_BENCHMARK_THRESHOLDS,
) -> BenchmarkSummary:
    matched = [comparison for comparison in evaluation.comparisons if comparison.matched_name is not None]
    ratios = [comparison.ratio for comparison in matched if math.isfinite(comparison.ratio)]
    max_ratio = max(ratios) if ratios else None
    mean_ratio = sum(ratios) / len(ratios) if ratios else None
    mean_abs_log10_error = _mean_abs_log10_error(matched)

    bench = load_benchmark(evaluation.bench_file)
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


def _matrix_source_origin(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    origin = str(source_metadata.get("origin", "")).strip()
    if origin:
        return origin
    if bench.get("source_doi"):
        return "external_literature"
    return "unspecified"


def _matrix_source_reference(bench: dict) -> str:
    source_metadata = bench.get("source_metadata") or {}
    if bench.get("source_doi"):
        return str(bench["source_doi"])
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


def _matrix_external_data_status(bench: dict) -> str:
    source_origin = _matrix_source_origin(bench)
    has_measured = bool(bench.get("measured_volatiles"))
    if has_measured and (bench.get("source_doi") or source_origin.startswith("external")):
        return "external_quantitative"
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
        bench_path = ROOT / "data" / "benchmarks" / f"{bench.get('benchmark_id', 'unknown')}.json"

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


def _predicted_order_lookup(evaluation: BenchmarkEvaluation) -> List[str]:
    normalized_rows = []
    seen: set[str] = set()
    for name, value in sorted(evaluation.predicted_ppb.items(), key=lambda item: item[1], reverse=True):
        normalized = _normalize_name(name)
        if normalized in seen:
            continue
        seen.add(normalized)
        normalized_rows.append(name)
    return normalized_rows


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
            )
        summaries.append(summary)
    return summaries


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