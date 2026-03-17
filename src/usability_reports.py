from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from src.benchmark_validation import BenchmarkSummary, summarize_benchmarks
from src.inverse_design import FormulationResult, InverseDesigner
from src.smirks_engine import ReactionConditions
from src.validation_contract import DEFAULT_VALIDATION_CONTRACT


@dataclass(frozen=True)
class ValidatedEnvelopeReport:
    target_tag: str
    strict_ready_benchmarks: List[str]
    supported_benchmarks: int
    total_benchmarks: int
    matrix_only_benchmarks: List[str]
    warnings: List[str]
    next_priorities: List[str]


def build_validated_envelope_report(target_tag: str = "meaty") -> ValidatedEnvelopeReport:
    summaries = summarize_benchmarks(target_tag=target_tag)
    strict_ready = [summary.benchmark_id for summary in summaries if summary.strict_ready]
    matrix_only = [summary.benchmark_id for summary in summaries if summary.execution_path == "matrix_only"]
    warnings = []
    if matrix_only:
        warnings.append(
            "Matrix benchmarks are executable intake/headspace checks, but remain outside the strict gate and target snapshots."
        )
    warnings.append(
        "Benchmark-facing concentrations still use the FAST observable projection; Cantera remains diagnostic-reference-only."
    )
    warnings.append(
        "Peptide-bound and intact-protein accessibility remain outside the validated precursor envelope."
    )
    warnings.append(
        "The validated plant-matrix scope is currently limited to pea/soy matrix-only systems and not yet a broad release gate."
    )
    next_priorities = [
        "Expose matrix-state and projection explainability in user-facing artifacts.",
        "Promote domain-of-applicability warnings into routine CLI/report outputs.",
        "Replace bulk matrix retention with broader compound-aware observability across plant-matrix systems.",
    ]
    return ValidatedEnvelopeReport(
        target_tag=target_tag,
        strict_ready_benchmarks=strict_ready,
        supported_benchmarks=sum(1 for summary in summaries if summary.supported),
        total_benchmarks=len(summaries),
        matrix_only_benchmarks=matrix_only,
        warnings=warnings,
        next_priorities=next_priorities,
    )


def render_validated_envelope_markdown(report: ValidatedEnvelopeReport) -> str:
    lines = [
        "# Validated Envelope",
        "",
        f"Target tag: {report.target_tag}",
        f"Supported benchmarks: {report.supported_benchmarks}/{report.total_benchmarks}",
        f"Strict-ready benchmarks: {', '.join(report.strict_ready_benchmarks) if report.strict_ready_benchmarks else 'none'}",
        f"Matrix-only executable benchmarks: {', '.join(report.matrix_only_benchmarks) if report.matrix_only_benchmarks else 'none'}",
        "",
        "## Warnings",
    ]
    for warning in report.warnings:
        lines.append(f"- {warning}")
    lines.extend([
        "",
        "## Next Priorities",
    ])
    for item in report.next_priorities:
        lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def build_formulation_explainability_payload(
    formulation: Dict[str, object],
    result: FormulationResult,
    *,
    target_tag: str,
    minimize_tag: str,
) -> Dict[str, object]:
    top_rows = []
    seen = set()
    for name, observable in sorted(result.predicted_ppb.items(), key=lambda item: item[1], reverse=True):
        if not isinstance(name, str) or name in seen:
            continue
        if name in result.projection_metadata:
            metadata = result.projection_metadata[name]
        else:
            continue
        seen.add(name)
        top_rows.append({
            "compound": name,
            "proxy_ppb": metadata.get("proxy_ppb", observable),
            "observable_ppb": metadata.get("observable_ppb", observable),
            "observable_ratio": metadata.get("proxy_to_observable_ratio", 1.0),
            "matrix_factor": metadata.get("matrix_factor", 1.0),
            "headspace_factor": metadata.get("headspace_factor", 1.0),
            "volatile_class": metadata.get("volatile_class", "other"),
        })
        if len(top_rows) >= 8:
            break

    return {
        "formulation_name": formulation.get("name", "unknown"),
        "target_tag": target_tag,
        "minimize_tag": minimize_tag,
        "protein_type": formulation.get("protein_type", "free"),
        "matrix_explainability": result.matrix_explainability,
        "scores": {
            "target_score": result.target_score,
            "off_flavour_risk": result.off_flavour_risk,
            "safety_score": result.safety_score,
            "texture_risk": result.texture_risk,
            "lysine_budget": result.lysine_budget,
            "trapping_efficiency": result.trapping_efficiency,
        },
        "top_projection_rows": top_rows,
        "detected_targets": result.detected_targets,
        "detected_minimize": result.detected_minimize,
    }


def render_formulation_explainability_markdown(payload: Dict[str, object]) -> str:
    lines = [
        "# Formulation Explainability",
        "",
        f"Formulation: {payload['formulation_name']}",
        f"Protein type: {payload['protein_type']}",
        f"Target tag: {payload['target_tag']}",
        f"Minimize tag: {payload['minimize_tag']}",
        "",
        "## Matrix State",
    ]
    matrix = payload["matrix_explainability"]
    lines.extend([
        f"- Effective denaturation state: {matrix['effective_denaturation_state']:.3f}",
        f"- Lysine accessibility: {matrix['lysine_accessibility']:.3f}",
        f"- Cysteine accessibility: {matrix['cysteine_accessibility']:.3f}",
        f"- Bulk volatile retention: {matrix['bulk_volatile_retention']:.3f}",
        f"- Denaturation source: {matrix['denaturation_source']}",
        "",
        "## Scores",
    ])
    for key, value in payload["scores"].items():
        lines.append(f"- {key}: {value:.3f}")
    lines.extend([
        "",
        "## Projection Rows",
        "",
        "| Compound | Proxy ppb | Observable ppb | Obs/Proxy | Matrix | Headspace | Class |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ])
    for row in payload["top_projection_rows"]:
        lines.append(
            f"| {row['compound']} | {row['proxy_ppb']:.3f} | {row['observable_ppb']:.3f} | {row['observable_ratio']:.3f} | {row['matrix_factor']:.3f} | {row['headspace_factor']:.3f} | {row['volatile_class']} |"
        )
    return "\n".join(lines) + "\n"