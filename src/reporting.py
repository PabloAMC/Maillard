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
from contextlib import redirect_stdout
from pathlib import Path
from typing import List, Dict, Any, Optional

from src.pipeline import FormulationResult
from src.projection_metadata import ProjectionMetadataMap
from src.usability_reports import DomainWarning
from src.projection_utils import build_projection_rows, build_artifact_provenance
from src.presentation import (
    render_decision_summary_cli,
    render_flavor_axis_markdown,
    render_projection_rows_markdown,
    render_provenance_markdown,
)

SCHEMA_VERSION = "2026-03-18"


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
    }
    payload: Dict[str, str] = {}
    for key, path in references.items():
        if path.exists():
            payload[key] = _to_repo_relative(path, root)
    return payload


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
            "flagged_toxics": result.flagged_toxics,
            "radar": {k: float(v[0]) for k, v in result.radar.items()},
            "matrix_explainability": result.matrix_explainability,
            "confidence_metadata": result.confidence_metadata,
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
        f.write(f"- **Pyrazine Penalty:** {result.pyrazine_penalty:.2f}\n\n")

        if result.confidence_metadata:
            f.write("### Confidence & Support\n")
            f.write(f"- **tier:** {result.confidence_metadata.get('tier', 'unknown')}\n")
            f.write(f"- **score:** {result.confidence_metadata.get('score', 0):.1f}\n")
            f.write(f"- **benchmark_neighborhood:** {result.confidence_metadata.get('benchmark_neighborhood', 'unknown')}\n")
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

            compound_rows = result.confidence_metadata.get("compound_confidence", [])
            if compound_rows:
                f.write("### Compound Confidence\n")
                f.write("| Compound | Observable ppb | Tier | Score | Mode |\n")
                f.write("| :--- | ---: | :---: | ---: | :--- |\n")
                for row in compound_rows:
                    f.write(
                        f"| {row.get('compound', 'unknown')} | {float(row.get('observable_ppb', 0.0)):.2f} | {row.get('tier', 'unknown')} | {float(row.get('score', 0.0)):.1f} | {row.get('prediction_mode', 'unknown')} |\n"
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
                "confidence_tier": res.confidence_metadata.get("tier", "unknown"),
                "confidence_score": float(res.confidence_metadata.get("score", 0.0)),
                "benchmark_neighborhood": res.confidence_metadata.get("benchmark_neighborhood", "unknown"),
                "prediction_mode": res.confidence_metadata.get("prediction_mode", "unknown"),
            },
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
        f.write("| Formulation | Strecker Balance | Pyrazine Burden | Thiamine Pathway | Expected Furanones |\n")
        f.write("| :--- | ---: | ---: | :---: | :--- |\n")
        for res in results:
            flavor_axis = res.flavor_axis_summary or {}
            f.write(
                f"| {res.name} | {res.strecker_balance_score:.2f} | {res.pyrazine_burden:.2f} | {str(flavor_axis.get('thiamine_pathway_active', False))} | {', '.join(str(item) for item in flavor_axis.get('furanone_expected', [])) or '-'} |\n"
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
