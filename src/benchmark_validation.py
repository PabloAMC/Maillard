from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.conditions import ReactionConditions
from src.inverse_design import InverseDesigner
from src.precursor_resolver import resolve_many
from src.smirks_engine import SmirksEngine


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
    "2methyl3furanthiol": {"2methyl3furanthiolmft"},
    "2furfurylthiol": {"2furfurylthiolfft"},
    "bis2methyl3furyldisulfide": {"bis2methyl3furyldisulfide"},
    "pyrazine": {"25dimethylpyrazine", "23dimethylpyrazine", "dimethylpyrazine"},
}


@dataclass(frozen=True)
class BenchmarkMetadata:
    tier: str
    family: str
    execution_path: str
    notes: Optional[str] = None


@dataclass(frozen=True)
class BenchmarkIndexEntry:
    benchmark_id: str
    bench_file: Path
    tier: str
    family: str
    protein_type: str
    execution_path: str
    supported: bool
    reason: Optional[str]
    status: str
    strict_ready: bool


class BenchmarkNotSupportedError(RuntimeError):
    pass


@dataclass(frozen=True)
class CompoundComparison:
    compound: str
    measured_ppb: float
    predicted_ppb: float
    matched_name: Optional[str]
    uncertainty_pct: Optional[float]
    match_score: float = 0.0

    @property
    def ratio(self) -> float:
        smallest = min(self.measured_ppb, self.predicted_ppb)
        largest = max(self.measured_ppb, self.predicted_ppb)
        if smallest <= 0.0:
            return math.inf if largest > 0.0 else 1.0
        return largest / smallest


@dataclass(frozen=True)
class BenchmarkEvaluation:
    benchmark_id: str
    bench_file: Path
    supported: bool
    reason: Optional[str]
    predicted_ppb: Dict[str, float]
    comparisons: List[CompoundComparison]
    pearson_r: Optional[float]
    mae_ppb: Optional[float]

    @property
    def coverage(self) -> float:
        if not self.comparisons:
            return 0.0
        matched = sum(1 for comparison in self.comparisons if comparison.matched_name is not None)
        return matched / len(self.comparisons)


@dataclass(frozen=True)
class BenchmarkSummary:
    benchmark_id: str
    bench_file: Path
    tier: str
    family: str
    execution_path: str
    supported: bool
    reason: Optional[str]
    protein_type: str
    coverage: float
    matched_compounds: int
    total_compounds: int
    pearson_r: Optional[float]
    mae_ppb: Optional[float]
    max_ratio: Optional[float]
    mean_ratio: Optional[float]
    ranking_status: str
    scale_status: str
    overall_status: str
    strict_ready: bool
    blocking_issues: List[str]


@dataclass(frozen=True)
class BenchmarkThresholds:
    ranking_threshold: float = 0.85
    free_aa_ratio_threshold: float = 1.5
    matrix_ratio_threshold: float = 2.0


@dataclass(frozen=True)
class BenchmarkTargetSnapshot:
    benchmark_id: str
    bench_file: Path
    target_name: str
    target_type: str
    roles: List[str]
    predicted_ppb: float
    weighted_flux: float
    span: float
    depth: int
    headspace_observable: bool
    headspace_class: str
    henry_kaw_25c: Optional[float]
    henry_source_name: Optional[str]


DEFAULT_BENCHMARK_THRESHOLDS = BenchmarkThresholds()


def _infer_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    benchmark_id = bench.get("benchmark_id", "unknown")
    protein_type = bench.get("protein_type", "free")
    if protein_type != "free":
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="matrix_headspace",
            execution_path="matrix_only",
            notes="Matrix-only benchmark requiring a dedicated precursor-accessibility path.",
        )
    if "cys_" in benchmark_id:
        return BenchmarkMetadata(
            tier="PRIMARY",
            family="free_aa_sulfur",
            execution_path="free_precursor",
            notes="Free amino-acid sulfur benchmark.",
        )
    return BenchmarkMetadata(
        tier="SECONDARY",
        family="general",
        execution_path="free_precursor",
        notes=None,
    )


def get_benchmark_metadata(bench: dict) -> BenchmarkMetadata:
    metadata = bench.get("metadata") or {}
    inferred = _infer_benchmark_metadata(bench)
    return BenchmarkMetadata(
        tier=str(metadata.get("tier", inferred.tier)),
        family=str(metadata.get("family", inferred.family)),
        execution_path=str(metadata.get("execution_path", inferred.execution_path)),
        notes=metadata.get("notes", inferred.notes),
    )


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
    candidate_precursors = formulation["sugars"] + formulation["amino_acids"] + formulation["lipids"]
    if not candidate_precursors:
        skipped = ", ".join(formulation.get("_skipped_matrix_precursors", [])) or "none"
        return False, f"No resolvable free-precursor system in benchmark. Matrix-only precursors: {skipped}"

    try:
        resolve_many(candidate_precursors)
    except ValueError as exc:
        return False, str(exc)
    return True, None


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


def _build_comparisons(bench: dict, predicted_ppb: Dict[str, float]) -> List[CompoundComparison]:
    comparisons: List[CompoundComparison] = []
    for compound, measured in bench["measured_volatiles"].items():
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
) -> dict:
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)
    designer = InverseDesigner(target_tag=target_tag)

    names = formulation["sugars"] + formulation["amino_acids"] + formulation.get("additives", []) + formulation.get("lipids", [])
    precursors = resolve_many(names)
    steps = SmirksEngine(conditions).enumerate(precursors, max_generations=4)

    heuristic_barriers = {}
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        bar, _source, unc = designer.db.get_best_barrier(reactants, products, step.reaction_family or "unknown")
        k = conditions.get_rate_constant(step.reaction_family or "unknown", ea_override_kcal=bar)
        rt = 0.001987 * conditions.temperature_kelvin
        pre_exponential = 1e11
        bar_eff = -rt * math.log(k / pre_exponential) if k > 0 else 99.0
        rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)

    from src.recommend import Recommender, _canon

    initial_concentrations = {}
    ratios = formulation.get("molar_ratios", {})
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[_canon(precursor.smiles)] = qty

    rec = Recommender()
    return rec.predict_from_steps(
        steps,
        heuristic_barriers,
        initial_concentrations,
        temperature_kelvin=conditions.temperature_kelvin,
        time_minutes=formulation.get("time_minutes"),
        protein_type=formulation.get("protein_type", "free"),
        denaturation_state=formulation.get("denaturation_state", 0.5),
    )


def evaluate_benchmark(bench_file: Path | str, target_tag: str = DEFAULT_TARGET_TAG) -> BenchmarkEvaluation:
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
        )

    rec_result = _run_benchmark_recommendation(bench, target_tag=target_tag)
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
    ratio_threshold = (
        thresholds.free_aa_ratio_threshold
        if protein_type == "free"
        else thresholds.matrix_ratio_threshold
    )

    if not evaluation.supported:
        return BenchmarkSummary(
            benchmark_id=evaluation.benchmark_id,
            bench_file=evaluation.bench_file,
            tier="UNKNOWN",
            family="unknown",
            execution_path="unknown",
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
        )

    if len(matched) >= 3 and evaluation.pearson_r is not None:
        ranking_status = "pass" if evaluation.pearson_r >= thresholds.ranking_threshold else "fail"
    elif len(matched) > 0:
        ranking_status = "n/a"
    else:
        ranking_status = "fail"

    if max_ratio is None:
        scale_status = "fail"
    elif max_ratio <= ratio_threshold:
        scale_status = "pass"
    else:
        scale_status = "fail"

    overall_status = "pass"
    if evaluation.coverage < 1.0:
        overall_status = "coverage-gap"
    elif ranking_status == "fail":
        overall_status = "ranking-gap"
    elif scale_status == "fail":
        overall_status = "scale-gap"
    elif ranking_status == "n/a":
        overall_status = "partial-pass"

    blocking_issues: List[str] = []
    if evaluation.coverage < 1.0:
        blocking_issues.append(f"coverage {evaluation.coverage:.1%} < 100%")
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

    strict_ready = evaluation.coverage == 1.0 and ranking_status != "fail" and scale_status == "pass"

    bench = load_benchmark(evaluation.bench_file)
    metadata = get_benchmark_metadata(bench)

    return BenchmarkSummary(
        benchmark_id=evaluation.benchmark_id,
        bench_file=evaluation.bench_file,
        tier=metadata.tier,
        family=metadata.family,
        execution_path=metadata.execution_path,
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
    )


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
        metadata = get_benchmark_metadata(bench)
        evaluation = evaluate_benchmark(bench_path, target_tag=target_tag)
        summary = summarize_evaluation(
            evaluation,
            protein_type=bench.get("protein_type", "free"),
            thresholds=thresholds,
        )
        if not evaluation.supported:
            summary = BenchmarkSummary(
                benchmark_id=summary.benchmark_id,
                bench_file=summary.bench_file,
                tier=metadata.tier,
                family=metadata.family,
                execution_path=metadata.execution_path,
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
            )
        summaries.append(summary)
    return summaries


def snapshot_benchmark_targets(
    bench_file: Path | str,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> List[BenchmarkTargetSnapshot]:
    bench_path = Path(bench_file)
    bench = load_benchmark(bench_path)
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
                weighted_flux=float(target.get("weighted_flux", 0.0)),
                span=float(target.get("span", math.inf)),
                depth=int(target.get("depth", 0)),
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


def render_benchmark_targets_markdown(snapshots: Iterable[BenchmarkTargetSnapshot]) -> str:
    rows = list(snapshots)
    lines = [
        "# Benchmark Targets",
        "",
        "| Benchmark | Target | Type | Roles | Predicted ppb | Span | Depth | Headspace | Kaw 25C | Henry Name |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for row in rows:
        kaw = f"{row.henry_kaw_25c:.3e}" if row.henry_kaw_25c is not None else "n/a"
        lines.append(
            f"| {row.benchmark_id} | {row.target_name} | {row.target_type} | {', '.join(row.roles)} | {row.predicted_ppb:.3f} | {row.span:.3f} | {row.depth} | {row.headspace_class} | {kaw} | {row.henry_source_name or 'n/a'} |"
        )

    lines.extend([
        "",
        f"Target rows: {len(rows)}",
        f"Low-headspace rows: {sum(1 for row in rows if row.headspace_class == 'low_headspace')}",
    ])
    return "\n".join(lines) + "\n"


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
            supported=summary.supported,
            reason=summary.reason,
            status=summary.overall_status,
            strict_ready=summary.strict_ready,
        )
        for summary in summaries
    ]


def render_benchmark_index_markdown(entries: Iterable[BenchmarkIndexEntry]) -> str:
    rows = list(entries)
    lines = [
        "# Benchmark Index",
        "",
        "| Benchmark | Tier | Family | Protein | Execution Path | Supported | Status | Strict Ready | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for entry in rows:
        notes = entry.reason or "indexed"
        lines.append(
            f"| {entry.benchmark_id} | {entry.tier} | {entry.family} | {entry.protein_type} | {entry.execution_path} | {'yes' if entry.supported else 'no'} | {entry.status} | {'yes' if entry.strict_ready else 'no'} | {notes} |"
        )
    return "\n".join(lines) + "\n"


def render_benchmark_summary_markdown(summaries: Iterable[BenchmarkSummary]) -> str:
    rows = list(summaries)
    lines = [
        "# Benchmark Summary",
        "",
        "| Benchmark | Tier | Family | Protein | Execution Path | Status | Strict Ready | Coverage | Pearson R | Max Ratio | MAE ppb | Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for summary in rows:
        if not summary.supported:
            notes = summary.reason or "Unsupported"
            pearson = "n/a"
            max_ratio = "n/a"
            mae = "n/a"
            coverage = "0.0%"
            strict_ready = "no"
        else:
            notes = ", ".join(summary.blocking_issues) or "validated"
            pearson = f"{summary.pearson_r:.3f}" if summary.pearson_r is not None else "n/a"
            max_ratio = f"{summary.max_ratio:.3f}" if summary.max_ratio is not None else "n/a"
            mae = f"{summary.mae_ppb:.2f}" if summary.mae_ppb is not None else "n/a"
            coverage = f"{summary.coverage:.1%}"
            strict_ready = "yes" if summary.strict_ready else "no"

        lines.append(
            f"| {summary.benchmark_id} | {summary.tier} | {summary.family} | {summary.protein_type} | {summary.execution_path} | {summary.overall_status} | {strict_ready} | {coverage} | {pearson} | {max_ratio} | {mae} | {notes} |"
        )

    supported_count = sum(1 for summary in rows if summary.supported)
    pass_count = sum(1 for summary in rows if summary.overall_status in {"pass", "partial-pass"})
    strict_ready_count = sum(1 for summary in rows if summary.strict_ready)
    lines.extend([
        "",
        f"Supported benchmarks: {supported_count}/{len(rows)}",
        f"Benchmarks without blocking coverage/ranking gaps: {pass_count}/{len(rows)}",
        f"Strict-ready benchmarks: {strict_ready_count}/{len(rows)}",
    ])
    return "\n".join(lines) + "\n"