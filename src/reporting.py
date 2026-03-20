"""
src/reporting.py

Consolidated report generation for Maillard framework simulations.
Outputs machine-readable JSON and human-readable Markdown.
"""

import datetime
import hashlib
import io
import json
import platform
import shlex
import subprocess
import sys
from collections import defaultdict
from contextlib import redirect_stdout
from pathlib import Path
from typing import List, Dict, Any, Optional, Iterable

from src.pipeline import FormulationResult
from src.literature_learning_loop import build_literature_learning_loop_payload
from src.projection_metadata import ProjectionMetadataMap, normalize_projection_metadata_row
from src.safety import build_safety_reference_context
from src.literature_runtime import build_flavor_reference_policy_summary
from src.usability_reports import DomainWarning
from src.projection_utils import build_projection_rows, build_artifact_provenance
from src.presentation import (
    render_decision_summary_cli,
    render_flavor_axis_markdown,
    render_projection_rows_markdown,
    render_provenance_markdown,
)

SCHEMA_VERSION = "2026-03-18"


def _sorted_projection_metadata(result: FormulationResult) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for raw_row in (result.projection_metadata or {}).values():
        normalized = normalize_projection_metadata_row(
            raw_row,
            compound_fallback=str(raw_row.get("compound", "unknown")),
            observable_ppb_fallback=float(raw_row.get("observable_ppb", 0.0) or 0.0),
        )
        rows.append(dict(normalized))
    rows.sort(key=lambda row: float(row.get("observable_ppb", 0.0) or 0.0), reverse=True)
    return rows


def _evidence_ladder_flags(meta: Dict[str, Any]) -> Dict[str, bool]:
    evidence_state = str(meta.get("evidence_state", "")).lower()
    source = str(meta.get("calibration_source", "")).lower()
    strength = str(meta.get("calibration_evidence_strength", "")).lower()
    fallback = str(meta.get("calibration_fallback_mode", "")).lower()
    notes_blob = " ".join(
        [
            source,
            str(meta.get("calibration_notes", "")).lower(),
            str(meta.get("retention_runtime_mode", "")).lower(),
        ]
    )

    direct_anchor = evidence_state in {"externally_benchmarked", "internally_benchmarked"} or (
        strength == "literature_anchored" and fallback == "compound_specific"
    )
    transferred_prior = evidence_state == "transferred_prior" or (
        strength in {"conditional_literature_anchored", "process_state_mismatch"}
        or fallback in {"nearest_process_state", "compound_specific_process_state"}
        or "transfer" in source
        or "carryover" in source
        or "ratio" in source
    )
    computational_refinement = any(token in notes_blob for token in ["dft", "xtb", "qm", "semiempirical", "computational", "refinement"])
    mechanistic_surrogate = (
        not direct_anchor and not transferred_prior and not computational_refinement
    ) or evidence_state == "still_missing" or strength == "heuristic" or float(meta.get("melanoidin_trapping_factor", 1.0) or 1.0) < 1.0

    return {
        "direct_anchor": bool(direct_anchor),
        "transferred_prior": bool(transferred_prior),
        "mechanistic_surrogate": bool(mechanistic_surrogate),
        "computational_refinement": bool(computational_refinement),
    }


def _support_origin(meta: Dict[str, Any]) -> str:
    process_state = str(meta.get("process_state", "unknown"))
    calibration_state = str(meta.get("calibration_process_state", "unknown"))
    fallback = str(meta.get("calibration_fallback_mode", "class_level"))
    source = str(meta.get("calibration_source", "class_fallback")).lower()
    strength = str(meta.get("calibration_evidence_strength", "heuristic")).lower()
    extrusion_state = process_state in {"aqueous_pre_extrusion_model", "extrusion_structured"}

    if extrusion_state:
        if fallback == "nearest_process_state" or (calibration_state not in {"unknown", process_state}):
            return "lower_regime_transfer"
        if strength == "heuristic" or source == "class_fallback":
            return "extrusion_extrapolation"
        return "extrusion_specific_support"

    return "standard_matrix_support"


def _build_compound_evidence_ladder(result: FormulationResult, *, top_n: int = 8) -> List[Dict[str, Any]]:
    ladder_rows: List[Dict[str, Any]] = []
    for meta in _sorted_projection_metadata(result)[:top_n]:
        flags = _evidence_ladder_flags(meta)
        ladder_rows.append(
            {
                "compound": str(meta.get("compound", "unknown")),
                "observable_ppb": float(meta.get("observable_ppb", 0.0) or 0.0),
                "direct_anchor": flags["direct_anchor"],
                "transferred_prior": flags["transferred_prior"],
                "mechanistic_surrogate": flags["mechanistic_surrogate"],
                "computational_refinement": flags["computational_refinement"],
                "evidence_state": str(meta.get("evidence_state", "still_missing")),
                "target_class": str(meta.get("target_class", "unknown")),
                "decision_panel_source": str(meta.get("decision_panel_source", "")),
                "support_origin": _support_origin(meta),
                "calibration_source": str(meta.get("calibration_source", "class_fallback")),
                "calibration_evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                "calibration_fallback_mode": str(meta.get("calibration_fallback_mode", "class_level")),
            }
        )
    return ladder_rows


def _build_calibration_summary(result: FormulationResult, *, top_n: int = 5) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    for meta in _sorted_projection_metadata(result):
        source = str(meta.get("calibration_source", "class_fallback"))
        if not source:
            continue
        support_origin = _support_origin(meta)
        bucket = grouped.setdefault(
            f"{source}::{support_origin}",
            {
                "source": source,
                "support_origin": support_origin,
                "compounds": [],
                "observable_ppb_total": 0.0,
                "evidence_strength": str(meta.get("calibration_evidence_strength", "heuristic")),
                "fallback_mode": str(meta.get("calibration_fallback_mode", "class_level")),
            },
        )
        compound = str(meta.get("compound", "unknown"))
        if compound not in bucket["compounds"]:
            bucket["compounds"].append(compound)
        bucket["observable_ppb_total"] += float(meta.get("observable_ppb", 0.0) or 0.0)

    rows = sorted(grouped.values(), key=lambda row: float(row["observable_ppb_total"]), reverse=True)
    return rows[:top_n]


def _flatten_condition_values(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for nested in value.values():
            yield from _flatten_condition_values(nested)
        return
    if isinstance(value, (list, tuple, set)):
        for nested in value:
            yield from _flatten_condition_values(nested)
        return
    if value is None:
        return
    yield str(value)


def _build_benchmark_neighborhood_summary(
    result: FormulationResult,
    conditions_dict: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    confidence = result.confidence_metadata or {}
    neighborhood = str(confidence.get("benchmark_neighborhood", "unknown"))
    prediction_mode = str(confidence.get("prediction_mode", "unknown"))
    process_regime = str(confidence.get("process_regime", "unknown"))
    condition_tokens = " ".join(
        token.lower() for token in _flatten_condition_values(conditions_dict or {})
    )
    hydrolysate_proxy = any(token in condition_tokens for token in ["hydrolysate", "peptide", "glutathione"])

    if hydrolysate_proxy:
        category = "hydrolysate_proxy"
        summary = "Run depends on hydrolysate/peptide-like inputs that sit outside the primary free-precursor benchmark surface."
    elif process_regime in {"extrusion_like", "extrusion_heavy"}:
        category = "extrusion_regime"
        summary = str(confidence.get("process_regime_summary", "Run is being interpreted through an extrusion-like regime without direct benchmark closure."))
    elif neighborhood in {"primary_free_precursor", "free_precursor_partial_analogy"}:
        category = "free_system_anchor"
        summary = "Run is anchored primarily by the free-system benchmark family, with varying degrees of analogy completeness."
    elif neighborhood == "matrix_intake_only":
        category = "matrix_transfer"
        summary = "Run uses matrix intake/headspace support and transferred accessibility priors rather than direct ranking benchmarks."
    else:
        category = "structural_gap"
        summary = "Run still sits beyond the strongest benchmark neighborhood and should be treated as a structural-gap extrapolation."

    return {
        "benchmark_neighborhood": neighborhood,
        "category": category,
        "prediction_mode": prediction_mode,
        "process_regime": process_regime,
        "summary": summary,
    }


def _build_missing_data_summary(result: FormulationResult) -> Dict[str, Any]:
    items: List[str] = []
    hypothesis_only: List[str] = []
    structurally_unsupported: List[str] = []

    for meta in _sorted_projection_metadata(result)[:8]:
        compound = str(meta.get("compound", "unknown"))
        source = str(meta.get("calibration_source", "class_fallback"))
        strength = str(meta.get("calibration_evidence_strength", "heuristic"))
        if source == "class_fallback" or strength == "heuristic":
            structurally_unsupported.append(compound)
            items.append(f"{compound}: still relies on class-level fallback or heuristic calibration.")
        elif str(meta.get("calibration_fallback_mode", "")) == "nearest_process_state":
            hypothesis_only.append(compound)
            items.append(f"{compound}: uses nearest-process-state transfer rather than a direct benchmark-condition anchor.")

    for compound in result.flavor_axis_summary.get("furanone_missing", []) or []:
        hypothesis_only.append(str(compound))
        items.append(f"{compound}: mechanistically expected but still unobserved in the current prediction surface.")

    benchmark_neighborhood = str((result.confidence_metadata or {}).get("benchmark_neighborhood", "unknown"))
    if benchmark_neighborhood in {"matrix_intake_only", "exploratory_matrix", "sparse_precursor_analogy"}:
        items.append(
            "Benchmark neighborhood remains extrapolative relative to the primary free-precursor validation envelope."
        )

    process_regime = str((result.confidence_metadata or {}).get("process_regime", "unknown"))
    extrusion_panel = (result.confidence_metadata or {}).get("extrusion_observable_panel", {})
    if process_regime in {"extrusion_like", "extrusion_heavy"} and not bool(extrusion_panel.get("minimum_panel_ready", False)):
        missing_categories = [
            category.replace("_", " ")
            for category in ("meaty_positive", "off_notes", "severity_markers")
            if not extrusion_panel.get(category, {}).get("present")
        ]
        items.append(
            "Extrusion observable panel remains incomplete: missing support for " + ", ".join(missing_categories) + "."
        )

    deduped_items = list(dict.fromkeys(items))[:8]
    return {
        "items": deduped_items,
        "hypothesis_only_compounds": list(dict.fromkeys(hypothesis_only)),
        "structurally_unsupported_compounds": list(dict.fromkeys(structurally_unsupported)),
    }


def _strecker_support_marker(result: FormulationResult) -> str:
    score = float(result.strecker_balance_score)
    if score >= 0.75:
        return "strong"
    if score >= 0.4:
        return "moderate"
    return "weak"


def _sulfur_trapping_summary(result: FormulationResult) -> Dict[str, Any]:
    sulfur_rows = []
    for meta in _sorted_projection_metadata(result):
        volatile_class = str(meta.get("volatile_class", "")).lower()
        compound = str(meta.get("compound", "")).lower()
        if volatile_class == "sulfur" or any(token in compound for token in ["thiol", "sulfide", "sulfur", "methional", "thiazole", "thiophene"]):
            sulfur_rows.append(meta)

    if not sulfur_rows:
        return {
            "status": "not_applicable",
            "avg_trapping_factor": 1.0,
            "summary": "No sulfur-focused observable rows were present in this run.",
        }

    avg_trapping = sum(float(row.get("melanoidin_trapping_factor", 1.0) or 1.0) for row in sulfur_rows) / len(sulfur_rows)
    if avg_trapping < 0.55:
        status = "severe"
    elif avg_trapping < 0.85:
        status = "moderate"
    else:
        status = "mild"
    return {
        "status": status,
        "avg_trapping_factor": float(avg_trapping),
        "summary": f"Average sulfur trapping factor is {avg_trapping:.2f} across {len(sulfur_rows)} sulfur-linked observable rows.",
    }


def _build_safety_reference_summary(result: FormulationResult) -> Dict[str, Any]:
    analyte = "acrylamide" if "Acrylamide" in (result.flagged_toxics or []) or float(result.safety_score) > 0.0 else "acrylamide"
    return build_safety_reference_context(analyte=analyte)


def _build_flavor_reference_policy(result: FormulationResult) -> List[Dict[str, Any]]:
    policy_rows = result.flavor_axis_summary.get("flavor_reference_policy") if getattr(result, "flavor_axis_summary", None) else None
    if isinstance(policy_rows, list) and policy_rows:
        return [dict(row) for row in policy_rows]
    return build_flavor_reference_policy_summary()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _to_repo_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _safe_git_output(root: Path, args: List[str]) -> Optional[str]:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _build_scientific_surface(root: Path) -> Dict[str, str]:
    references = {
        "scientific_reference": root / "docs/reference/SCIENTIFIC_REFERENCE.md",
        "benchmark_summary": root / "results/validation/benchmark_summary.md",
        "validated_envelope": root / "results/validation/validated_envelope.md",
        "validation_overview": root / "results/validation/validation_overview.md",
        "matrix_decision_panel": root / "data/lit/matrix_decision_panel.json",
        "matrix_family_coverage_registry": root / "data/lit/matrix_family_coverage_registry.json",
        "benchmark_intake_registry": root / "data/lit/benchmark_intake_registry.json",
        "computational_priors": root / "data/lit/computational_priors.json",
        "slr_incorporation_matrix": root / "data/lit/slr_incorporation_matrix.json",
        "flavor_reference_payloads": root / "data/lit/flavor_reference_payloads.json",
        "process_state_calibrations": root / "data/lit/process_state_calibrations.json",
        "retention_reference_payloads": root / "data/lit/retention_reference_payloads.json",
        "process_gap_registry": root / "data/lit/process_gap_registry.json",
        "safety_reference_payloads": root / "data/lit/safety_reference_payloads.json",
        "primary_benchmark_protocol": root / "docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md",
        "primary_benchmark_contract": root / "data/protocols/ppi_spi_primary_benchmark_contract.json",
        "literature_learning_loop": root / "results/validation/literature_learning_loop.md",
        "literature_learning_loop_json": root / "results/validation/literature_learning_loop.json",
        "literature_runtime_templates": root / "results/validation/literature_runtime_templates.json",
        "matrix_target_status": root / "results/validation/matrix_target_status.md",
        "matrix_target_status_json": root / "results/validation/matrix_target_status.json",
        "matrix_family_coverage": root / "results/validation/matrix_family_coverage.md",
        "matrix_family_coverage_json": root / "results/validation/matrix_family_coverage.json",
        "refinement_watchlist": root / "results/validation/refinement_watchlist.md",
        "refinement_watchlist_json": root / "results/validation/refinement_watchlist.json",
        "offline_dft_jobs": root / "results/validation/offline_dft_jobs.json",
        "family_sensitivity": root / "results/validation/family_sensitivity.md",
        "family_sensitivity_json": root / "results/validation/family_sensitivity.json",
        "p3_global_sensitivity": root / "results/validation/p3_global_sensitivity.md",
        "p3_global_sensitivity_json": root / "results/validation/p3_global_sensitivity.json",
        "cheap_refinement_screening": root / "results/validation/cheap_refinement_screening.md",
        "cheap_refinement_screening_json": root / "results/validation/cheap_refinement_screening.json",
        "selective_dft_plan": root / "results/validation/selective_dft_plan.md",
        "selective_dft_plan_json": root / "results/validation/selective_dft_plan.json",
        "p3_offline_dft_jobs": root / "results/validation/p3_offline_dft_jobs.json",
        "refinement_impact": root / "results/validation/refinement_impact.md",
        "refinement_impact_json": root / "results/validation/refinement_impact.json",
        "refinement_surrogate_patches": root / "data/lit/refinement_surrogate_patches.json",
        "reaction_benchmark_set": root / "data/lit/reaction_benchmark_set.json",
        "mlp_candidate_registry": root / "data/lit/mlp_candidate_registry.json",
        "mlp_external_benchmark_evidence": root / "data/lit/mlp_external_benchmark_evidence.json",
        "p4_geometry_benchmark_set": root / "data/lit/p4_geometry_benchmark_set.json",
        "p4_geometry_benchmark": root / "results/validation/p4_geometry_benchmark.md",
        "p4_geometry_benchmark_json": root / "results/validation/p4_geometry_benchmark.json",
        "p4_geometry_assessment": root / "results/validation/p4_geometry_assessment.md",
        "p4_geometry_assessment_json": root / "results/validation/p4_geometry_assessment.json",
        "p4_reaction_benchmark": root / "results/validation/p4_reaction_benchmark.md",
        "p4_reaction_benchmark_json": root / "results/validation/p4_reaction_benchmark.json",
        "p4_mlp_assessment": root / "results/validation/p4_mlp_assessment.md",
        "p4_mlp_assessment_json": root / "results/validation/p4_mlp_assessment.json",
        "p4_external_mlp_landscape": root / "results/validation/p4_external_mlp_landscape.md",
        "p4_external_mlp_landscape_json": root / "results/validation/p4_external_mlp_landscape.json",
        "p4_adoption_notes": root / "results/validation/p4_adoption_notes.md",
        "p4_adoption_notes_json": root / "results/validation/p4_adoption_notes.json",
    }
    payload: Dict[str, str] = {}
    for key, path in references.items():
        if path.exists():
            payload[key] = _to_repo_relative(path, root)
    return payload


def _build_literature_evidence_summary(root: Optional[Path] = None) -> Dict[str, Any]:
    repo_root = root or _repo_root()
    intake_path = repo_root / "data" / "lit" / "benchmark_intake_registry.json"
    if not intake_path.exists():
        return {}

    with open(intake_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    eligible = payload.get("eligible_references", []) or []
    structural_gaps = payload.get("structural_gaps", []) or []
    ready_refs = [entry for entry in eligible if str(entry.get("status", "")).startswith("ready_for_")]
    no_primary_data_refs = [entry for entry in eligible if not bool(entry.get("requires_primary_data", False))]

    modules: Dict[str, int] = defaultdict(int)
    matrix_families: Dict[str, int] = defaultdict(int)
    kinds: Dict[str, int] = defaultdict(int)
    for entry in ready_refs:
        kinds[str(entry.get("kind", "unknown"))] += 1
        matrix_families[str(entry.get("matrix_family", "unknown"))] += 1
        for module in entry.get("target_modules", []) or []:
            modules[str(module)] += 1

    return {
        "source": _to_repo_relative(intake_path, repo_root),
        "eligible_reference_count": len(eligible),
        "ready_reference_count": len(ready_refs),
        "closable_without_primary_data_count": len(no_primary_data_refs),
        "structural_gap_count": len(structural_gaps),
        "ready_reference_ids": [str(entry.get("id", "unknown")) for entry in ready_refs[:8]],
        "ready_by_kind": dict(sorted(kinds.items())),
        "ready_by_module": dict(sorted(modules.items())),
        "ready_by_matrix_family": dict(sorted(matrix_families.items())),
        "structural_gap_ids": [str(entry.get("id", entry.get("gap_id", "unknown"))) for entry in structural_gaps[:8]],
    }


def _build_literature_learning_loop_summary(root: Optional[Path] = None) -> Dict[str, Any]:
    repo_root = root or _repo_root()
    payload = build_literature_learning_loop_payload(repo_root)
    return dict(payload.get("summary", {}))


def generate_report(
    result: FormulationResult, 
    warnings: List[DomainWarning], 
    conditions_dict: Dict[str, Any],
    output_dir: Optional[Path] = None,
    campaign_metadata: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Generates a consolidated report (JSON + MD) for a formulation result.
    If output_dir is None, creates a timestamped folder in 'reports/'.
    Returns the path to the report directory.
    """
    if output_dir is None:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"reports/run_{timestamp}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = build_artifact_provenance(
        artifact_kind="single_run_report",
        output_dir=output_dir,
        inputs=conditions_dict,
        campaign_metadata=campaign_metadata,
    )
    compound_evidence_ladder = _build_compound_evidence_ladder(result)
    calibration_summary = _build_calibration_summary(result)
    missing_data_summary = _build_missing_data_summary(result)
    benchmark_neighborhood_summary = _build_benchmark_neighborhood_summary(result, conditions_dict)
    safety_reference_summary = _build_safety_reference_summary(result)
    flavor_reference_policy = _build_flavor_reference_policy(result)
    literature_evidence_summary = _build_literature_evidence_summary()
    literature_learning_loop_summary = _build_literature_learning_loop_summary()
    
    # 1. Save JSON Report
    json_path = output_dir / "report.json"
    report_data = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "provenance": provenance,
        "inputs": conditions_dict,
        "results": {
            "name": result.name,
            "target_score": float(result.target_score),
            "off_flavour_risk": float(result.off_flavour_risk),
            "safety_score": float(result.safety_score),
            "lysine_budget": float(result.lysine_budget),
            "trapping_efficiency": float(result.trapping_efficiency),
            "mft_to_furfural_ratio": float(result.mft_to_furfural_ratio),
            "meaty_quality_penalty": float(result.meaty_quality_penalty),
            "strecker_balance_score": float(result.strecker_balance_score),
            "strecker_gap_penalty": float(result.strecker_gap_penalty),
            "pyrazine_propensity": float(result.pyrazine_propensity),
            "pyrazine_burden": float(result.pyrazine_burden),
            "pyrazine_penalty": float(result.pyrazine_penalty),
            "furanone_penalty": float(result.furanone_penalty),
            "flagged_toxics": result.flagged_toxics,
            "radar": {k: float(v[0]) for k, v in result.radar.items()},
            "matrix_explainability": result.matrix_explainability,
            "confidence_metadata": result.confidence_metadata,
            "compound_evidence_ladder": compound_evidence_ladder,
            "calibration_summary": calibration_summary,
            "missing_data_summary": missing_data_summary,
            "benchmark_neighborhood_summary": benchmark_neighborhood_summary,
            "safety_reference_summary": safety_reference_summary,
            "flavor_reference_policy": flavor_reference_policy,
            "literature_evidence_summary": literature_evidence_summary,
            "literature_learning_loop_summary": literature_learning_loop_summary,
            "projection_metadata": dict(result.projection_metadata),
            "flavor_axis_summary": result.flavor_axis_summary,
            "predicted_ppb": {k: float(v) for k, v in result.predicted_ppb.items()},
            "detected_targets": result.detected_targets,
            "detected_minimize": result.detected_minimize
        },
        "domain_warnings": [
            {"category": w.category, "level": w.level, "message": w.message}
            for w in warnings
        ]
    }
    
    with open(json_path, "w") as f:
        json.dump(report_data, f, indent=4, default=str)
        
    # 2. Save Markdown Report
    md_path = output_dir / "report.md"
    with open(md_path, "w") as f:
        f.write(f"# Maillard Simulation Report - {result.name}\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 1. Input Formulation & Conditions\n")
        f.write("| Parameter | Value |\n")
        f.write("| :--- | :--- |\n")
        for k, v in conditions_dict.items():
            f.write(f"| {k} | {v} |\n")
        f.write("\n")
        
        f.write("## 2. Decision Summary\n")
        f.write("```text\n")
        with io.StringIO() as buf, redirect_stdout(buf):
            render_decision_summary_cli(result, warnings)
            f.write(buf.getvalue())
        f.write("```\n\n")
        
        f.write("## 3. Detailed Results\n")
        f.write(f"- **Target Score:** {result.target_score:.2f}\n")
        f.write(f"- **Off-Flavour Risk:** {result.off_flavour_risk:.2f}\n")
        f.write(f"- **Safety Score:** {result.safety_score:.2f}\n\n")
        f.write(f"- **MFT/Furfural Ratio:** {result.mft_to_furfural_ratio:.4f}\n")
        f.write(f"- **Meaty Quality Penalty:** {result.meaty_quality_penalty:.2f}\n\n")
        f.write(f"- **Strecker Balance Score:** {result.strecker_balance_score:.2f}\n")
        f.write(f"- **Strecker Gap Penalty:** {result.strecker_gap_penalty:.2f}\n")
        f.write(f"- **Pyrazine Propensity:** {result.pyrazine_propensity:.2f}\n")
        f.write(f"- **Pyrazine Burden:** {result.pyrazine_burden:.2f}\n")
        f.write(f"- **Pyrazine Penalty:** {result.pyrazine_penalty:.2f}\n")
        f.write(f"- **Furanone Penalty:** {result.furanone_penalty:.2f}\n\n")

        if result.confidence_metadata:
            f.write("### Confidence & Support\n")
            f.write(f"- **tier:** {result.confidence_metadata.get('tier', 'unknown')}\n")
            f.write(f"- **score:** {result.confidence_metadata.get('score', 0):.1f}\n")
            f.write(f"- **benchmark_neighborhood:** {result.confidence_metadata.get('benchmark_neighborhood', 'unknown')}\n")
            f.write(f"- **decision_mode:** {result.confidence_metadata.get('decision_mode', 'directional_hypothesis')}\n")
            f.write(f"- **prediction_mode:** {result.confidence_metadata.get('prediction_mode', 'unknown')}\n")
            f.write(f"- **recommended_posture:** {result.confidence_metadata.get('recommended_posture', '')}\n")
            for factor in result.confidence_metadata.get("dominant_factors", []):
                f.write(f"- **factor:** {factor}\n")
            f.write("\n")

            calibration = result.confidence_metadata.get("calibration_diagnostics", {})
            if calibration:
                f.write("### Calibration Diagnostics\n")
                f.write(f"- **supported_envelope:** {calibration.get('supported_envelope', False)}\n")
                f.write(f"- **summary:** {calibration.get('summary', '')}\n")
                axes = calibration.get("extrapolation_axes", [])
                if axes:
                    f.write(f"- **extrapolation_axes:** {', '.join(str(axis) for axis in axes)}\n")
                f.write("\n")

            f.write("### Benchmark Neighborhood\n")
            f.write(f"- **benchmark_neighborhood:** {benchmark_neighborhood_summary.get('benchmark_neighborhood', 'unknown')}\n")
            f.write(f"- **category:** {benchmark_neighborhood_summary.get('category', 'unknown')}\n")
            f.write(f"- **prediction_mode:** {benchmark_neighborhood_summary.get('prediction_mode', 'unknown')}\n")
            f.write(f"- **summary:** {benchmark_neighborhood_summary.get('summary', '')}\n\n")

            if calibration_summary:
                f.write("### Calibration Summary\n")
                f.write("| Source | Support Origin | Evidence | Fallback | Compounds | Observable ppb |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- | ---: |\n")
                for row in calibration_summary:
                    f.write(
                        f"| {row.get('source', 'unknown')} | {row.get('support_origin', 'standard_matrix_support')} | {row.get('evidence_strength', 'unknown')} | {row.get('fallback_mode', 'unknown')} | {', '.join(str(item) for item in row.get('compounds', []))} | {float(row.get('observable_ppb_total', 0.0)):.2f} |\n"
                    )
                f.write("\n")

            extrusion_panel = result.confidence_metadata.get("extrusion_observable_panel", {})
            if result.confidence_metadata.get("process_regime") in {"extrusion_like", "extrusion_heavy"} and extrusion_panel:
                f.write("### Extrusion Observable Panel\n")
                f.write("| Category | Present | Missing | Ready |\n")
                f.write("| :--- | :--- | :--- | :---: |\n")
                for category in ("meaty_positive", "off_notes", "severity_markers"):
                    row = extrusion_panel.get(category, {})
                    f.write(
                        f"| {category.replace('_', ' ')} | {', '.join(row.get('present', [])) or '-'} | {', '.join(row.get('missing', [])) or '-'} | {'yes' if row.get('present') else '-'} |\n"
                    )
                f.write(f"\n- **minimum_panel_ready:** {extrusion_panel.get('minimum_panel_ready', False)}\n\n")

            compound_rows = result.confidence_metadata.get("compound_confidence", [])
            if compound_rows:
                f.write("### Compound Confidence\n")
                f.write("| Compound | Observable ppb | Tier | Score | Mode | Reachability | Calibration Source | Observable Assumption |\n")
                f.write("| :--- | ---: | :---: | ---: | :--- | :--- | :--- | :--- |\n")
                for row in compound_rows:
                    f.write(
                        f"| {row.get('compound', 'unknown')} | {float(row.get('observable_ppb', 0.0)):.2f} | {row.get('tier', 'unknown')} | {float(row.get('score', 0.0)):.1f} | {row.get('prediction_mode', 'unknown')} | {row.get('reachability_status', 'merely_plausible')} | {row.get('calibration_source', 'unknown')} | {row.get('observable_assumption_summary', 'unknown')} |\n"
                    )
                f.write("\n")

            aggregate_rows = result.confidence_metadata.get("aggregate_confidence", {})
            if aggregate_rows:
                f.write("### Aggregate Sensory Confidence\n")
                f.write("| Tag | Sensory Score | Supporting Compounds | Tier | Mode |\n")
                f.write("| :--- | ---: | ---: | :---: | :--- |\n")
                for tag, row in aggregate_rows.items():
                    f.write(
                        f"| {tag} | {float(row.get('score', 0.0)):.2f} | {int(row.get('support_count', 0))} | {row.get('tier', 'unknown')} | {row.get('prediction_mode', 'unknown')} |\n"
                    )
                f.write("\n")

        if compound_evidence_ladder:
            f.write("### Compound Evidence Ladder\n")
            f.write("| Compound | Class | Evidence State | Direct Anchor | Transferred Prior | Mechanistic Surrogate | Computational Refinement | Support Origin | Source |\n")
            f.write("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- |\n")
            for row in compound_evidence_ladder:
                f.write(
                    f"| {row.get('compound', 'unknown')} | {row.get('target_class', 'unknown')} | {row.get('evidence_state', 'still_missing')} | {'yes' if row.get('direct_anchor') else '-'} | {'yes' if row.get('transferred_prior') else '-'} | {'yes' if row.get('mechanistic_surrogate') else '-'} | {'yes' if row.get('computational_refinement') else '-'} | {row.get('support_origin', 'standard_matrix_support')} | {row.get('calibration_source', 'unknown')} |\n"
                )
            f.write("\n")

        f.write("### Missing Data\n")
        missing_items = missing_data_summary.get("items", [])
        if missing_items:
            for item in missing_items:
                f.write(f"- {item}\n")
        else:
            f.write("- No high-priority missing-data flags were triggered for the top reported compounds.\n")
        f.write("\n")

        safety_defaults = safety_reference_summary.get("default_entries", [])
        if safety_defaults:
            f.write("### Safety Reference Context\n")
            f.write("| Visibility | Kind | Source | Summary |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for row in safety_defaults:
                f.write(
                    f"| {row.get('report_visibility', 'default')} | {row.get('kind', 'unknown')} | {row.get('source_citation', 'unknown')} | {row.get('summary', '')} |\n"
                )
            extended_count = len(safety_reference_summary.get("extended_entries", []))
            if extended_count:
                f.write(f"\n- Extended safety provenance entries available in JSON: {extended_count}\n")
            f.write("\n")

        if flavor_reference_policy:
            f.write("### Flavor Reference Policy\n")
            f.write("| Compound | Pipeline Role | Benchmark Role | Source |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for row in flavor_reference_policy:
                f.write(
                    f"| {row.get('compound', 'unknown')} | {row.get('pipeline_role', 'reference_only')} | {row.get('benchmark_role', 'unknown')} | {row.get('source_citation', 'unknown')} |\n"
                )
            f.write("\n")

        if literature_evidence_summary:
            f.write("### Literature Evidence Summary\n")
            f.write(f"- **source:** {literature_evidence_summary.get('source', 'unknown')}\n")
            f.write(f"- **eligible_reference_count:** {literature_evidence_summary.get('eligible_reference_count', 0)}\n")
            f.write(f"- **ready_reference_count:** {literature_evidence_summary.get('ready_reference_count', 0)}\n")
            f.write(f"- **closable_without_primary_data_count:** {literature_evidence_summary.get('closable_without_primary_data_count', 0)}\n")
            f.write(f"- **structural_gap_count:** {literature_evidence_summary.get('structural_gap_count', 0)}\n")
            if literature_evidence_summary.get('ready_reference_ids'):
                f.write(f"- **ready_reference_ids:** {', '.join(str(item) for item in literature_evidence_summary.get('ready_reference_ids', []))}\n")
            if literature_evidence_summary.get('structural_gap_ids'):
                f.write(f"- **structural_gap_ids:** {', '.join(str(item) for item in literature_evidence_summary.get('structural_gap_ids', []))}\n")
            f.write("\n")

        if literature_learning_loop_summary:
            f.write("### Literature Learning Loop Summary\n")
            for key in [
                "ready_reference_count",
                "encoded_runtime_reference_count",
                "template_queue_count",
                "matrix_family_count",
                "intake_structural_gap_count",
                "process_gap_count",
            ]:
                f.write(f"- **{key}:** {literature_learning_loop_summary.get(key, 0)}\n")
            if literature_learning_loop_summary.get("matrix_prior_families"):
                f.write(f"- **matrix_prior_families:** {', '.join(str(item) for item in literature_learning_loop_summary.get('matrix_prior_families', []))}\n")
            f.write("\n")

        if result.confidence_metadata:
            sensitivity = result.confidence_metadata.get("sensitivity_summary", {})
            if sensitivity:
                f.write("### Sensitivity Summary\n")
                f.write(f"- **mode:** {sensitivity.get('mode', 'unknown')}\n")
                f.write(f"- **evaluated_perturbations:** {sensitivity.get('evaluated_perturbations', 0)}\n")
                for item in sensitivity.get("ranking_drivers", [])[:3]:
                    f.write(
                        f"- **ranking_driver:** {item.get('input', 'unknown')} | {item.get('perturbation', 'unknown')} | "
                        f"Δdecision {float(item.get('decision_delta', 0.0)):+.2f} | Δsafety {float(item.get('safety_delta', 0.0)):+.2f}\n"
                    )
                for item in sensitivity.get("safety_drivers", [])[:3]:
                    f.write(
                        f"- **safety_driver:** {item.get('input', 'unknown')} | {item.get('perturbation', 'unknown')} | "
                        f"Δsafety {float(item.get('safety_delta', 0.0)):+.2f}\n"
                    )
                f.write("\n")
        
        if result.detected_targets:
            f.write("### Predicted Desirable Targets\n")
            f.write("| Compound | Barrier |\n")
            f.write("| :--- | :--- |\n")
            for t in result.detected_targets[:10]:
                f.write(f"| {t} | - |\n")
            f.write("\n")

        if getattr(result, "projection_metadata", None):
            projection_rows = build_projection_rows(result)
            f.write(render_projection_rows_markdown(projection_rows, heading="### Projection Calibration", variant="compact"))

            if projection_rows:
                f.write("### Trust Surface\n")
                f.write(f"- **decision_mode:** {result.confidence_metadata.get('decision_mode', 'directional_hypothesis')}\n")
                f.write(f"- **benchmark_neighborhood:** {benchmark_neighborhood_summary.get('benchmark_neighborhood', 'unknown')}\n")
                f.write(f"- **extrapolation_axes:** {', '.join(str(axis) for axis in calibration.get('extrapolation_axes', [])) if calibration else 'none'}\n")
                top_row = projection_rows[0]
                f.write(f"- **top_calibration_source:** {top_row.get('calibration_source', 'unknown')}\n")
                f.write(f"- **top_observable_assumption:** {top_row.get('observable_assumption_summary', 'unknown')}\n")
                f.write(f"- **top_reachability_status:** {top_row.get('reachability_status', 'merely_plausible')}\n\n")

        if getattr(result, "flavor_axis_summary", None):
            f.write(render_flavor_axis_markdown(result.flavor_axis_summary, heading="### Flavor Axis Diagnostics", variant="detailed"))

        f.write("## 4. Analytical Metadata\n")
        f.write("### Matrix Explainability\n")
        for k, v in result.matrix_explainability.items():
            f.write(f"- **{k}:** {v}\n")
        f.write("\n")
        f.write(render_provenance_markdown(provenance))
            
    return output_dir

def generate_comparison_report(
    results: List[FormulationResult],
    conditions_list: List[Dict[str, Any]],
    warnings_list: Optional[List[List[DomainWarning]]] = None,
    output_dir: Optional[Path] = None,
    campaign_metadata: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Generates a side-by-side comparison report for multiple formulation results.
    """
    if output_dir is None:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"reports/comparison_{timestamp}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = build_artifact_provenance(
        artifact_kind="formulation_comparison",
        output_dir=output_dir,
        inputs=conditions_list,
        campaign_metadata=campaign_metadata,
    )
    
    # 1. Save JSON Comparison
    json_path = output_dir / "comparison.json"
    comparison_data = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "provenance": provenance,
        "campaign": campaign_metadata or {},
        "runs": []
    }
    
    for res, cond in zip(results, conditions_list):
        compound_evidence_ladder = _build_compound_evidence_ladder(res)
        calibration_summary = _build_calibration_summary(res)
        missing_data_summary = _build_missing_data_summary(res)
        benchmark_neighborhood_summary = _build_benchmark_neighborhood_summary(res, cond)
        sulfur_trapping_summary = _sulfur_trapping_summary(res)
        safety_reference_summary = _build_safety_reference_summary(res)
        flavor_reference_policy = _build_flavor_reference_policy(res)
        literature_evidence_summary = _build_literature_evidence_summary()
        literature_learning_loop_summary = _build_literature_learning_loop_summary()
        comparison_data["runs"].append({
            "name": res.name,
            "inputs": cond,
            "metrics": {
                "target_score": float(res.target_score),
                "off_flavour_risk": float(res.off_flavour_risk),
                "safety_score": float(res.safety_score),
                "lysine_budget": float(res.lysine_budget),
                "trapping_efficiency": float(res.trapping_efficiency),
                "mft_to_furfural_ratio": float(res.mft_to_furfural_ratio),
                "meaty_quality_penalty": float(res.meaty_quality_penalty),
                "strecker_balance_score": float(res.strecker_balance_score),
                "strecker_gap_penalty": float(res.strecker_gap_penalty),
                "pyrazine_propensity": float(res.pyrazine_propensity),
                "pyrazine_burden": float(res.pyrazine_burden),
                "pyrazine_penalty": float(res.pyrazine_penalty),
                "furanone_penalty": float(res.furanone_penalty),
                "confidence_tier": res.confidence_metadata.get("tier", "unknown"),
                "confidence_score": float(res.confidence_metadata.get("score", 0.0)),
                "benchmark_neighborhood": res.confidence_metadata.get("benchmark_neighborhood", "unknown"),
                "prediction_mode": res.confidence_metadata.get("prediction_mode", "unknown"),
                "strecker_support_marker": _strecker_support_marker(res),
            },
            "compound_evidence_ladder": compound_evidence_ladder,
            "calibration_summary": calibration_summary,
            "missing_data_summary": missing_data_summary,
            "benchmark_neighborhood_summary": benchmark_neighborhood_summary,
            "sulfur_trapping_summary": sulfur_trapping_summary,
            "safety_reference_summary": safety_reference_summary,
            "flavor_reference_policy": flavor_reference_policy,
            "literature_evidence_summary": literature_evidence_summary,
            "literature_learning_loop_summary": literature_learning_loop_summary,
            "flavor_axis_summary": res.flavor_axis_summary,
        })
    
    with open(json_path, "w") as f:
        json.dump(comparison_data, f, indent=4, default=str)
        
    # 2. Save Markdown Comparison
    md_path = output_dir / "comparison.md"
    with open(md_path, "w") as f:
        f.write("# Maillard Formulation Comparison Report\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 1. Metric Overview\n")
        f.write("| Metric | " + " | ".join([res.name for res in results]) + " |\n")
        f.write("| :--- | " + " | ".join([":---:"] * len(results)) + " |\n")
        
        f.write("| **Target Score** | " + " | ".join([f"{res.target_score:.2f}" for res in results]) + " |\n")
        f.write("| **Off-Flavour Risk** | " + " | ".join([f"{res.off_flavour_risk:.2f}" for res in results]) + " |\n")
        f.write("| **Safety Score** | " + " | ".join([f"{res.safety_score:.2f}" for res in results]) + " |\n")
        f.write("| **Lysine Budget** | " + " | ".join([f"{res.lysine_budget:.1f}%" for res in results]) + " |\n")
        f.write("| **Trapping Eff.** | " + " | ".join([f"{res.trapping_efficiency:.1f}%" for res in results]) + " |\n")
        f.write("| **MFT/Furfural Ratio** | " + " | ".join([f"{res.mft_to_furfural_ratio:.4f}" for res in results]) + " |\n")
        f.write("| **Meaty Quality Penalty** | " + " | ".join([f"{res.meaty_quality_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Strecker Balance** | " + " | ".join([f"{res.strecker_balance_score:.2f}" for res in results]) + " |\n")
        f.write("| **Strecker Penalty** | " + " | ".join([f"{res.strecker_gap_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Pyrazine Burden** | " + " | ".join([f"{res.pyrazine_burden:.2f}" for res in results]) + " |\n")
        f.write("| **Pyrazine Penalty** | " + " | ".join([f"{res.pyrazine_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Furanone Penalty** | " + " | ".join([f"{res.furanone_penalty:.2f}" for res in results]) + " |\n")
        f.write("| **Confidence** | " + " | ".join([f"{res.confidence_metadata.get('tier', 'unknown')} ({float(res.confidence_metadata.get('score', 0.0)):.0f})" for res in results]) + " |\n")
        f.write("| **Prediction Mode** | " + " | ".join([str(res.confidence_metadata.get('prediction_mode', 'unknown')) for res in results]) + " |\n")
        f.write("\n")
        
        f.write("## 2. Competitive Highlights\n")
        best_target = max(results, key=lambda x: x.target_score)
        best_safety = min(results, key=lambda x: x.safety_score)
        best_risk = min(results, key=lambda x: x.off_flavour_risk)
        
        f.write(f"- 🏆 **Highest Target Score:** {best_target.name} ({best_target.target_score:.2f})\n")
        f.write(f"- 🛡️ **Safest Formulation:** {best_safety.name} ({best_safety.safety_score:.2f})\n")
        f.write(f"- 🍃 **Lowest Off-Flavour Risk:** {best_risk.name} ({best_risk.off_flavour_risk:.2f})\n\n")

        f.write("## 3. Cross-Marker Context\n")
        f.write("| Formulation | Strecker Balance | Strecker Support | Pyrazine Burden | Sulfur Trapping | Furanone Penalty | Benchmark Neighborhood | Thiamine Pathway | Thiamine Source | Expected Furanones | Missing Furanones |\n")
        f.write("| :--- | ---: | :---: | ---: | :--- | ---: | :--- | :---: | :--- | :--- | :--- |\n")
        for res in results:
            flavor_axis = res.flavor_axis_summary or {}
            trapping = _sulfur_trapping_summary(res)
            benchmark_summary = _build_benchmark_neighborhood_summary(res)
            f.write(
                f"| {res.name} | {res.strecker_balance_score:.2f} | {_strecker_support_marker(res)} | {res.pyrazine_burden:.2f} | {trapping.get('status', 'n/a')} ({float(trapping.get('avg_trapping_factor', 1.0)):.2f}) | {res.furanone_penalty:.2f} | {benchmark_summary.get('category', 'unknown')} | {str(flavor_axis.get('thiamine_pathway_active', False))} | {str(flavor_axis.get('thiamine_availability_source', 'unknown'))} | {', '.join(str(item) for item in flavor_axis.get('furanone_expected', [])) or '-'} | {', '.join(str(item) for item in flavor_axis.get('furanone_missing', [])) or '-'} |\n"
            )
        f.write("\n")

        f.write("## 4. Calibration Contrast\n")
        f.write("| Formulation | Decision Mode | Benchmark Neighborhood | Top Calibration Source | Observable Assumption | Extrapolation Axes | Missing Data Flags | Benchmark Summary |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- | ---: | :--- |\n")
        for res, cond in zip(results, conditions_list):
            calibration_summary = _build_calibration_summary(res)
            top_calibration = calibration_summary[0] if calibration_summary else {
                "source": "class_fallback",
                "evidence_strength": "heuristic",
            }
            missing_data = _build_missing_data_summary(res)
            benchmark_summary = _build_benchmark_neighborhood_summary(res, cond)
            projection_rows = build_projection_rows(res)
            top_projection = projection_rows[0] if projection_rows else {}
            diagnostics = (res.confidence_metadata or {}).get("calibration_diagnostics", {})
            f.write(
                f"| {res.name} | {res.confidence_metadata.get('decision_mode', 'directional_hypothesis')} | {benchmark_summary.get('benchmark_neighborhood', 'unknown')} | {top_calibration.get('source', 'class_fallback')} | {top_projection.get('observable_assumption_summary', 'unknown')} | {', '.join(str(axis) for axis in diagnostics.get('extrapolation_axes', [])) or 'none'} | {len(missing_data.get('items', []))} | {benchmark_summary.get('summary', '')} |\n"
            )
        f.write("\n")

        f.write("## Trust Surface\n")
        f.write("| Formulation | Prediction Mode | Decision Mode | Top Reachability | Support Origin | Recommended Use |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")
        for res in results:
            projection_rows = build_projection_rows(res)
            top_projection = projection_rows[0] if projection_rows else {}
            f.write(
                f"| {res.name} | {res.confidence_metadata.get('prediction_mode', 'unknown')} | {res.confidence_metadata.get('decision_mode', 'directional_hypothesis')} | {top_projection.get('reachability_status', 'merely_plausible')} | {top_projection.get('support_origin', 'standard_matrix_support')} | {res.confidence_metadata.get('recommended_posture', '')} |\n"
            )
        f.write("\n")

        for index, res in enumerate(results):
            f.write(f"### {res.name}\n")
            f.write("```text\n")
            with io.StringIO() as buf, redirect_stdout(buf):
                item_warnings = warnings_list[index] if warnings_list and index < len(warnings_list) else []
                render_decision_summary_cli(res, item_warnings)
                f.write(buf.getvalue())
            f.write("```\n\n")
        f.write(render_provenance_markdown(provenance))
            
    return output_dir


def generate_campaign_report(
    campaign_spec: Dict[str, Any],
    results: List[FormulationResult],
    conditions_list: List[Dict[str, Any]],
    run_artifacts: List[Dict[str, Any]],
    warnings_list: Optional[List[List[DomainWarning]]] = None,
    output_dir: Optional[Path] = None,
) -> Path:
    campaign_metadata = campaign_spec.get("campaign", {})
    comparison_dir = generate_comparison_report(
        results=results,
        conditions_list=conditions_list,
        warnings_list=warnings_list,
        output_dir=output_dir,
        campaign_metadata=campaign_metadata,
    )

    leaderboard = []
    for result, conditions in zip(results, conditions_list):
        leaderboard.append(
            {
                "name": result.name,
                "protein_type": conditions.get("protein_type", "free"),
                "ph": float(conditions.get("ph", 0.0)),
                "temp": float(conditions.get("temp", 0.0)),
                "target_score": float(result.target_score),
                "off_flavour_risk": float(result.off_flavour_risk),
                "safety_score": float(result.safety_score),
                "confidence_tier": result.confidence_metadata.get("tier", "unknown"),
                "prediction_mode": result.confidence_metadata.get("prediction_mode", "unknown"),
            }
        )
    leaderboard.sort(key=lambda item: item["target_score"], reverse=True)

    provenance = build_artifact_provenance(
        artifact_kind="campaign_screen",
        output_dir=comparison_dir,
        inputs=campaign_spec,
        campaign_metadata=campaign_metadata,
    )

    campaign_payload = {
        "timestamp": datetime.datetime.now().isoformat(),
        "schema_version": SCHEMA_VERSION,
        "campaign": campaign_metadata,
        "shared_conditions": campaign_spec.get("shared_conditions", {}),
        "leaderboard": leaderboard,
        "run_artifacts": run_artifacts,
        "comparison_artifacts": {
            "markdown": str(comparison_dir / "comparison.md"),
            "json": str(comparison_dir / "comparison.json"),
        },
        "provenance": provenance,
    }

    json_path = comparison_dir / "campaign.json"
    with open(json_path, "w") as handle:
        json.dump(campaign_payload, handle, indent=4, default=str)

    md_path = comparison_dir / "campaign.md"
    with open(md_path, "w") as handle:
        handle.write("# Maillard Campaign Report\n\n")
        handle.write(f"**Campaign:** {campaign_metadata.get('name', 'Unnamed campaign')}\n\n")
        handle.write(f"**Objective:** {campaign_metadata.get('objective', campaign_metadata.get('description', ''))}\n\n")
        if campaign_metadata.get("audience"):
            handle.write(f"**Audience:** {campaign_metadata.get('audience')}\n\n")
        handle.write("## 1. Shared Conditions\n")
        handle.write("| Parameter | Value |\n")
        handle.write("| :--- | :--- |\n")
        for key, value in campaign_spec.get("shared_conditions", {}).items():
            handle.write(f"| {key} | {value} |\n")
        handle.write("\n")
        handle.write("## 2. Leaderboard\n")
        handle.write("| Formulation | Protein | pH | Temp C | Target | Off-flavour | Safety | Confidence | Mode |\n")
        handle.write("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :---: | :--- |\n")
        for row in leaderboard:
            handle.write(
                f"| {row['name']} | {row['protein_type']} | {row['ph']:.2f} | {row['temp']:.1f} | {row['target_score']:.2f} | {row['off_flavour_risk']:.2f} | {row['safety_score']:.2f} | {row['confidence_tier']} | {row['prediction_mode']} |\n"
            )
        handle.write("\n")
        handle.write("## 3. Generated Artifacts\n")
        handle.write(f"- comparison markdown: {comparison_dir / 'comparison.md'}\n")
        handle.write(f"- comparison json: {comparison_dir / 'comparison.json'}\n")
        handle.write(f"- campaign json: {json_path}\n")
        for artifact in run_artifacts:
            handle.write(
                f"- run artifact: {artifact.get('name', 'unknown')} -> {artifact.get('directory', '')}\n"
            )
        handle.write("\n")
        handle.write(render_provenance_markdown(provenance))

    return comparison_dir
