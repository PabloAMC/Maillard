from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from src.benchmark_validation import (
    DEFAULT_TARGET_TAG,
    benchmark_to_conditions,
    benchmark_to_formulation,
    evaluate_benchmark,
    get_benchmark_files,
    get_benchmark_metadata,
    load_benchmark,
    summarize_evaluation,
)
from src.dft_refinement_contract import DFTTargetCandidate, build_offline_dft_job, triage_dft_candidate
from src.precursor_resolver import resolve_many
from src.results_db import ResultsDB
from src.smirks_engine import SmirksEngine


def _benchmark_priority_weight(bench: dict, summary_status: str, reference_signal_origin: str) -> float:
    contract = bench.get("matrix_ranking_contract") or {}
    roles = {str(item.get("role", "")).strip().lower() for item in contract.get("observable_targets", [])}
    if "desirable_marker" in roles and bool(contract.get("adverse_markers")):
        base = 1.0
    elif "desirable_marker" in roles:
        base = 0.9
    else:
        base = 0.75

    gap_weight = {
        "scale-gap": 1.0,
        "ranking-gap": 0.95,
        "coverage-gap": 0.9,
        "pass-no-ranking": 0.7,
        "pass": 0.65,
        "unsupported": 0.4,
    }.get(str(summary_status), 0.6)
    signal_weight = 1.0 if str(reference_signal_origin) == "measured_volatiles" else 0.8
    return base * gap_weight * signal_weight


def _source_gap_score(source: str) -> float:
    normalized = str(source).strip().lower()
    if normalized.startswith("db:wb97m-v"):
        return 0.2
    if normalized.startswith("db:r2scan"):
        return 0.3
    if normalized.startswith("db:mace-off"):
        return 0.55
    if normalized.startswith("db:xtb"):
        return 0.7
    if normalized.startswith("heuristic"):
        return 0.95
    return 0.8


def _requires_explicit_solvent(family: str) -> bool:
    normalized = str(family).strip().lower()
    return any(token in normalized for token in ["hydrolysis", "schiff", "rearrangement", "proton"])


def _requires_irc(family: str) -> bool:
    normalized = str(family).strip().lower()
    return any(token in normalized for token in ["rearrangement", "elimination", "strecker", "thiol_addition"])


def _enumerate_benchmark_steps(bench: dict) -> List[Any]:
    formulation = benchmark_to_formulation(bench)
    names = formulation["sugars"] + formulation["amino_acids"] + formulation.get("additives", []) + formulation.get("lipids", [])
    if not names:
        return []
    precursors = resolve_many(names)
    conditions = benchmark_to_conditions(bench)
    return SmirksEngine(conditions).enumerate(precursors, max_generations=4)


def build_refinement_watchlist(
    benchmark_files: Optional[Iterable[Path | str]] = None,
    *,
    target_tag: str = DEFAULT_TARGET_TAG,
) -> Dict[str, Any]:
    bench_files = list(benchmark_files) if benchmark_files is not None else get_benchmark_files()
    db = ResultsDB()
    family_rows: Dict[str, Dict[str, Any]] = {}

    for bench_file in bench_files:
        bench = load_benchmark(bench_file)
        metadata = get_benchmark_metadata(bench)
        if metadata.execution_path not in {"free_precursor", "matrix_precursor_augmented"}:
            continue

        steps = _enumerate_benchmark_steps(bench)
        if not steps:
            continue

        evaluation = evaluate_benchmark(bench_file, target_tag=target_tag)
        summary = summarize_evaluation(evaluation, protein_type=str(bench.get("protein_type", "free")))
        bench_weight = _benchmark_priority_weight(bench, summary.overall_status, summary.reference_signal_origin)

        for step in steps:
            family = str(step.reaction_family or "unknown").strip()
            if not family or family == "unknown":
                continue
            reactants = [species.smiles for species in step.reactants]
            products = [species.smiles for species in step.products]
            _barrier, source, uncertainty = db.get_best_barrier(reactants, products, family)

            row = family_rows.setdefault(
                family,
                {
                    "reaction_family": family,
                    "occurrence_count": 0,
                    "weighted_occurrence": 0.0,
                    "benchmark_ids": set(),
                    "execution_paths": set(),
                    "protein_types": set(),
                    "reference_signal_origins": set(),
                    "barrier_sources": defaultdict(int),
                    "uncertainties": [],
                    "gap_scores": [],
                },
            )
            row["occurrence_count"] += 1
            row["weighted_occurrence"] += bench_weight
            row["benchmark_ids"].add(str(bench.get("benchmark_id", Path(bench_file).stem)))
            row["execution_paths"].add(metadata.execution_path)
            row["protein_types"].add(str(bench.get("protein_type", "free")))
            row["reference_signal_origins"].add(summary.reference_signal_origin)
            row["barrier_sources"][str(source)] += 1
            row["uncertainties"].append(float(uncertainty))
            row["gap_scores"].append(_source_gap_score(source))

    if not family_rows:
        return {"summary": {"candidate_count": 0, "run_now": 0, "watchlist": 0, "reject": 0}, "candidates": [], "offline_jobs": []}

    max_weight = max(float(row["weighted_occurrence"]) for row in family_rows.values())
    max_benchmark_span = max(len(row["benchmark_ids"]) for row in family_rows.values())
    candidates: List[Dict[str, Any]] = []
    offline_jobs: List[Dict[str, Any]] = []

    for family, row in family_rows.items():
        dominant_source = max(row["barrier_sources"].items(), key=lambda item: item[1])[0]
        sensitivity_score = float(row["weighted_occurrence"]) / max(max_weight, 1.0e-9)
        benchmark_impact_score = len(row["benchmark_ids"]) / max(max_benchmark_span, 1)
        evidence_gap_score = sum(row["gap_scores"]) / max(len(row["gap_scores"]), 1)
        benchmark_ids = sorted(str(item) for item in row["benchmark_ids"])

        candidate = DFTTargetCandidate(
            reaction_id=f"family::{family.lower()}",
            reaction_family=family,
            sensitivity_score=min(1.0, sensitivity_score),
            benchmark_impact_score=min(1.0, benchmark_impact_score),
            evidence_gap_score=min(1.0, evidence_gap_score),
            current_barrier_source=dominant_source,
            rationale=(
                f"Family appears in {len(benchmark_ids)} benchmark-visible systems ({', '.join(benchmark_ids[:4])}) "
                f"with current source {dominant_source}."
            ),
            requires_explicit_solvent=_requires_explicit_solvent(family),
            requires_irc=_requires_irc(family),
            current_confidence="heuristic" if dominant_source.lower().startswith("heuristic") else "mixed",
        )
        decision = triage_dft_candidate(candidate)
        row_payload = {
            "reaction_id": candidate.reaction_id,
            "reaction_family": family,
            "occurrence_count": int(row["occurrence_count"]),
            "benchmark_ids": benchmark_ids,
            "execution_paths": sorted(str(item) for item in row["execution_paths"]),
            "protein_types": sorted(str(item) for item in row["protein_types"]),
            "reference_signal_origins": sorted(str(item) for item in row["reference_signal_origins"]),
            "current_barrier_source": dominant_source,
            "avg_uncertainty_kcal": sum(row["uncertainties"]) / max(len(row["uncertainties"]), 1),
            "sensitivity_score": candidate.sensitivity_score,
            "benchmark_impact_score": candidate.benchmark_impact_score,
            "evidence_gap_score": candidate.evidence_gap_score,
            "decision": decision.decision,
            "priority_tier": decision.priority_tier,
            "rationale": decision.rationale,
            "surrogate_plan": decision.surrogate_plan,
            "required_outputs": list(decision.required_outputs),
            "expected_benchmark_impact": f"Can affect {len(benchmark_ids)} benchmark-visible systems; mixed/meaty matrix candidates are prioritized via weighted occurrence.",
        }
        candidates.append(row_payload)
        if decision.decision == "run_now":
            job = build_offline_dft_job(candidate, artifact_prefix=f"watchlist_{family.lower()}")
            offline_jobs.append(
                {
                    "reaction_id": job.reaction_id,
                    "priority_tier": job.priority_tier,
                    "solvent_name": job.solvent_name,
                    "temperature_k": job.temperature_k,
                    "charge": job.charge,
                    "spin": job.spin,
                    "use_explicit_solvent": job.use_explicit_solvent,
                    "n_water": job.n_water,
                    "generate_irc": job.generate_irc,
                    "expected_outputs": list(job.expected_outputs),
                    "artifact_prefix": job.artifact_prefix,
                }
            )

    candidates.sort(key=lambda row: ({"run_now": 0, "defer": 1, "reject": 2}.get(row["decision"], 3), -(0.4 * row["sensitivity_score"] + 0.4 * row["benchmark_impact_score"] + 0.2 * row["evidence_gap_score"]), row["reaction_family"]))
    return {
        "summary": {
            "candidate_count": len(candidates),
            "run_now": sum(1 for row in candidates if row["decision"] == "run_now"),
            "watchlist": sum(1 for row in candidates if row["decision"] == "defer"),
            "reject": sum(1 for row in candidates if row["decision"] == "reject"),
        },
        "candidates": candidates,
        "offline_jobs": offline_jobs,
    }


def render_refinement_watchlist_markdown(payload: Dict[str, Any]) -> str:
    lines = [
        "# Refinement Watchlist",
        "",
        "| Reaction Family | Decision | Tier | Source | Sensitivity | Benchmark Impact | Evidence Gap | Benchmarks |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in payload.get("candidates", []):
        lines.append(
            f"| {row['reaction_family']} | {row['decision']} | {row['priority_tier']} | {row['current_barrier_source']} | {row['sensitivity_score']:.2f} | {row['benchmark_impact_score']:.2f} | {row['evidence_gap_score']:.2f} | {', '.join(row['benchmark_ids'][:4])} |"
        )
    summary = payload.get("summary", {})
    lines.extend([
        "",
        f"Candidates: {int(summary.get('candidate_count', 0))}",
        f"Run-now candidates: {int(summary.get('run_now', 0))}",
        f"Deferred watchlist candidates: {int(summary.get('watchlist', 0))}",
    ])
    return "\n".join(lines) + "\n"