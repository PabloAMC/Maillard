"""
src/benchmark_evaluator.py — Benchmark prediction engine and statistical helpers.

Extracted from src/benchmark_validation.py (Priority 2 decomposition).
Contains all functions responsible for *running* a benchmark prediction and
computing comparison statistics.  No reporting, assertion, or markdown
generation lives here.
"""
from __future__ import annotations

import math
import re
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from src.benchmark_types import (
    BenchmarkMetadata,
    BenchmarkNotSupportedError,
    CompoundComparison,
    BenchmarkEvaluation,
)
from src.benchmark_registry import (
    ROOT,
    DEFAULT_TARGET_TAG,
    load_benchmark,
    get_benchmark_files,
    get_benchmark_metadata,
    benchmark_to_conditions,
    benchmark_to_formulation,
    resolve_thermodynamic_gating_mode,
    _is_supported_formulation,
)

from src.text_utils import normalize_compound_name as _normalize_name
from src.barrier_constants import effective_barrier_from_rate_constant
from src.headspace import HeadspaceModel
from src.kinetics import KineticsEngine
from src.lipid_oxidation import (
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    predict_hexanal_generation,
)
from src.matrix_calibration_registry import (
    describe_matrix_calibration,
    determine_matrix_process_state,
)
from src.pipeline import MaillardPipeline
from src.precursor_resolver import resolve_many
from src.projection import DEFAULT_PROJECTION_STRATEGY
from src.projection_metadata import ProjectionMetadataMap, make_projection_metadata_row
from src.smirks_engine import SmirksEngine
from src.validation_contract import BenchmarkThresholds, DEFAULT_VALIDATION_CONTRACT
from src.safety import predict_acrylamide
from src.projection_utils import build_projection_rows
from src.matrix_targets import get_compound_panel_entry


# Re-export profiles so consumers that used to import from benchmark_validation
# still work without changes even when they import from this module directly.
MATRIX_BENCHMARK_PROFILES = {
    "pea_iso": {"lipid_profile": PEA_LIPID_PROFILE},
    "soy_iso": {"lipid_profile": SOY_LIPID_PROFILE},
}

MATRIX_BENCHMARK_BASE_MARKER_YIELDS = {
    "Hexanal": 0.205,
    "2-Pentylfuran": 0.502,
    "1-Hexanol": 0.063,
    "Nonanal": 0.150,
}

THERMODYNAMIC_GATING_THRESHOLD_KCAL = 30.0

# ---------------------------------------------------------------------------
# Name-matching helpers
# ---------------------------------------------------------------------------

BENCHMARK_NAME_ALIASES: Dict[str, set] = {
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
    "3methylbutanal": {"3methylbutanal", "isovaleraldehyde"},
    "2methylbutanal": {"2methylbutanal", "2methylbutyraldehyde"},
    "methional": {"3methylthiopropanal", "methional"},
    "acrylamide": {"acrylamide", "2-propenamide"},
}


def _tokenize_name(name: str) -> List[str]:
    return [token for token in re.split(r"[^a-z0-9]+", name.lower()) if token]


def _best_prediction_match(
    target_name: str, predicted_ppb: Dict[str, float]
) -> tuple[Optional[str], float, float]:
    target_norm = _normalize_name(target_name)
    target_tokens = {token for token in _tokenize_name(target_name) if len(token) >= 4}
    target_aliases = BENCHMARK_NAME_ALIASES.get(target_norm, set())
    best_name: Optional[str] = None
    best_score = -1.0

    for candidate_name in predicted_ppb:
        candidate_norm = _normalize_name(candidate_name)
        if not candidate_norm or len(candidate_norm) < 4:
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


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def _pearson(values_a: Iterable[float], values_b: Iterable[float]) -> Optional[float]:
    data_a = list(values_a)
    data_b = list(values_b)
    if len(data_a) < 3 or len(data_a) != len(data_b):
        return None
    mean_a = sum(data_a) / len(data_a)
    mean_b = sum(data_b) / len(data_b)
    centered_a = [v - mean_a for v in data_a]
    centered_b = [v - mean_b for v in data_b]
    numerator = sum(a * b for a, b in zip(centered_a, centered_b))
    denom_a = math.sqrt(sum(a * a for a in centered_a))
    denom_b = math.sqrt(sum(b * b for b in centered_b))
    if denom_a == 0.0 or denom_b == 0.0:
        return None
    return numerator / (denom_a * denom_b)


def _mean_abs_log10_error(comparisons: Iterable[CompoundComparison]) -> Optional[float]:
    matched = [c for c in comparisons if c.matched_name is not None]
    errors = [
        abs(math.log10(c.predicted_ppb / c.measured_ppb))
        for c in matched
        if c.measured_ppb > 0.0 and c.predicted_ppb > 0.0
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
        "max_ratio": float(
            configured.get("max_ratio", thresholds.ratio_threshold_for(protein_type))
        ),
        "mean_abs_log10_error": float(
            configured.get(
                "mean_abs_log10_error",
                thresholds.mean_abs_log10_error_threshold_for(protein_type),
            )
        ),
    }


def _build_comparisons(
    bench: dict, predicted_ppb: Dict[str, float]
) -> List[CompoundComparison]:
    comparisons: List[CompoundComparison] = []
    signal_map = (
        bench.get("measured_volatiles") or bench.get("reference_volatiles") or {}
    )
    for compound, measured in signal_map.items():
        matched_name, predicted_value, match_score = _best_prediction_match(
            compound, predicted_ppb
        )
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


# ---------------------------------------------------------------------------
# Prediction runners
# ---------------------------------------------------------------------------

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

    names = (
        formulation["sugars"]
        + formulation["amino_acids"]
        + formulation.get("additives", [])
        + formulation.get("lipids", [])
    )
    precursors = resolve_many(names)
    steps = SmirksEngine(conditions).enumerate(precursors, max_generations=4)

    heuristic_barriers: Dict[str, tuple] = {}
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        bar, _source, unc = designer.db.get_best_barrier(
            reactants, products, step.reaction_family or "unknown"
        )
        barrier_for_rate = float(bar)
        if gating_mode == "on":
            thermo = kinetics.get_reaction_thermo(
                reactants, products, conditions.temperature_kelvin
            )
            dg = float(thermo.get("delta_g_kcal_mol", 0.0))
            if dg > THERMODYNAMIC_GATING_THRESHOLD_KCAL:
                barrier_for_rate = 99.0
            else:
                barrier_for_rate = max(barrier_for_rate, max(0.5, dg + 0.5))
        k = conditions.get_rate_constant(
            step.reaction_family or "unknown", ea_override_kcal=barrier_for_rate
        )
        bar_eff = effective_barrier_from_rate_constant(
            k, conditions.temperature_kelvin, step.reaction_family or "unknown"
        )
        rxn_key = (
            f"{'+'.join(sorted(r.smiles for r in step.reactants))}"
            f"->{'+'.join(sorted(p.smiles for p in step.products))}"
        )
        heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)

    from src.recommend import Recommender
    from src.chem_utils import canonicalize_smiles

    initial_concentrations: Dict[str, float] = {}
    ratios = formulation.get("molar_ratios", {})
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[
            canonicalize_smiles(precursor.smiles, fallback_to_original=True, strip_salts=True)
        ] = qty

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
        smiles_map = {
            "Hexanal": "CCCCCC=O",
            "Nonanal": "CCCCCCCCC=O",
            "1-Hexanol": "CCCCCCO",
            "2-Pentylfuran": "CCCCC1=CC=CO1",
        }
        try:
            from rdkit import Chem
        except ImportError:
            Chem = None  # type: ignore[assignment]

        for compound_name, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
            smi = smiles_map.get(compound_name)
            if smi and Chem is not None:
                mol = Chem.MolFromSmiles(smi)
                mw = (
                    sum(atom.GetMass() for atom in mol.GetAtoms()) if mol else 100.0
                )
                target_proxy_ppb = oxidation_load_ppb * float(yield_factor)
                molar_conc = target_proxy_ppb / (
                    mw * DEFAULT_PROJECTION_STRATEGY.ppb_conversion_factor
                )
                canon = canonicalize_smiles(smi, fallback_to_original=True, strip_salts=True)
                initial_concentrations[canon] = initial_concentrations.get(canon, 0.0) + molar_conc

    process_state = str(
        (bench.get("process_metadata") or {}).get(
            "state",
            determine_matrix_process_state(
                temperature_celsius=float(conditions.temperature_celsius),
                time_minutes=float(formulation.get("time_minutes", 60.0)),
                water_activity=conditions.water_activity,
            ),
        )
    )

    rec = Recommender()
    return rec.predict_from_steps(
        steps,
        heuristic_barriers,
        initial_concentrations,
        temperature_kelvin=conditions.temperature_kelvin,
        time_minutes=float(formulation.get("time_minutes", 60.0)),
        protein_type=protein_type,
        denaturation_state=float(bench.get("denaturation_state", 0.5)),
        process_state=process_state,
    )


def _run_matrix_only_benchmark_prediction(bench: dict) -> dict:
    protein_type = str(bench.get("protein_type", "free"))
    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is None:
        raise BenchmarkNotSupportedError(
            f"No matrix-only predictor for protein_type={protein_type}"
        )

    conditions = bench["conditions"]
    process_state = str(
        (bench.get("process_metadata") or {}).get(
            "state",
            determine_matrix_process_state(
                temperature_celsius=float(conditions["temp_C"]),
                time_minutes=float(conditions["time_min"]),
                water_activity=conditions.get("aw"),
            ),
        )
    )
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
            water_activity=conditions.get("aw"),
        )
        calibration = describe_matrix_calibration(
            compound,
            protein_type=protein_type,
            process_state=process_state,
        )
        panel_entry = get_compound_panel_entry(compound) or {}
        calibration_factor = float(calibration.get("calibration_observable_factor") or 1.0)
        release_factor = (
            headspace_factor / calibration_factor if calibration_factor > 0.0 else headspace_factor
        )
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
                **panel_entry,
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
        "species_names": {c: c for c in predicted_ppb},
    }


# ---------------------------------------------------------------------------
# Top-level evaluation entry points
# ---------------------------------------------------------------------------

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
            reference_signal_origin=(
                "measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles"
            ),
        )

    metadata = get_benchmark_metadata(bench)
    if metadata.family == "safety":
        conditions = benchmark_to_conditions(bench)
        asn_conc = 0.0
        sugar_conc = 0.0
        for name, data in bench["precursors"].items():
            n_low = name.lower()
            if "asparagine" in n_low or "asn" in n_low:
                asn_conc = data["concentration_mM"]
            if any(
                s in n_low
                for s in ["ribose", "glucose", "fructose", "maltose", "xylose", "sugar"]
            ):
                sugar_conc += data["concentration_mM"]
        time_min = bench["conditions"].get("time_min", 20.0)
        safety_res = predict_acrylamide(
            asparagine_mM=asn_conc,
            reducing_sugar_mM=sugar_conc,
            temp_C=conditions.temperature_celsius,
            time_min=time_min,
            pH=conditions.pH,
        )
        rec_result = {
            "predicted_ppb": {"acrylamide": safety_res.acrylamide_ppb},
            "projection_metadata": {},
        }
    elif metadata.execution_path == "matrix_only":
        rec_result = _run_matrix_only_benchmark_prediction(bench)
    else:
        rec_result = _run_benchmark_recommendation(
            bench,
            target_tag=target_tag,
            thermodynamic_gating=thermodynamic_gating,
        )
    comparisons = _build_comparisons(bench, rec_result["predicted_ppb"])

    matched_comparisons = [c for c in comparisons if c.matched_name is not None]
    measured_values = [c.measured_ppb for c in matched_comparisons]
    predicted_values = [c.predicted_ppb for c in matched_comparisons]
    pearson_r = _pearson(measured_values, predicted_values)
    mae_ppb = None
    if matched_comparisons:
        mae_ppb = sum(
            abs(m - p) for m, p in zip(measured_values, predicted_values)
        ) / len(matched_comparisons)

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
        reference_signal_origin=(
            "measured_volatiles" if bench.get("measured_volatiles") else "reference_volatiles"
        ),
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
    pseudo_path = ROOT / "results" / "validation" / f"{payload_id}.synthetic_benchmark.json"
    return _evaluate_loaded_benchmark(
        normalized,
        bench_path=pseudo_path,
        target_tag=target_tag,
        thermodynamic_gating=thermodynamic_gating,
    )
