"""
src/reporting.py

Consolidated report generation for Maillard framework simulations.
Outputs machine-readable JSON and human-readable Markdown.
"""

import json
import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

from src.inverse_design import FormulationResult
from src.usability_reports import DomainWarning, render_decision_summary_cli

def generate_report(
    result: FormulationResult, 
    warnings: List[DomainWarning], 
    conditions_dict: Dict[str, Any],
    output_dir: Optional[Path] = None
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
    
    # 1. Save JSON Report
    json_path = output_dir / "report.json"
    report_data = {
        "timestamp": datetime.datetime.now().isoformat(),
        "inputs": conditions_dict,
        "results": {
            "name": result.name,
            "target_score": float(result.target_score),
            "off_flavour_risk": float(result.off_flavour_risk),
            "safety_score": float(result.safety_score),
            "lysine_budget": float(result.lysine_budget),
            "trapping_efficiency": float(result.trapping_efficiency),
            "flagged_toxics": result.flagged_toxics,
            "radar": {k: float(v[0]) for k, v in result.radar.items()},
            "matrix_explainability": result.matrix_explainability,
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
        json.dump(report_data, f, indent=4)
        
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
        # We'll use a string buffer or just capture stdout if we want exact same box
        # But better to have a clean MD version. For now, let's embed the CLI box for visual parity.
        import io
        from contextlib import redirect_stdout
        with io.StringIO() as buf, redirect_stdout(buf):
            render_decision_summary_cli(result, warnings)
            f.write(buf.getvalue())
        f.write("```\n\n")
        
        f.write("## 3. Detailed Results\n")
        f.write(f"- **Target Score:** {result.target_score:.2f}\n")
        f.write(f"- **Off-Flavour Risk:** {result.off_flavour_risk:.2f}\n")
        f.write(f"- **Safety Score:** {result.safety_score:.2f}\n\n")
        
        if result.detected_targets:
            f.write("### Predicted Desirable Targets\n")
            f.write("| Compound | Barrier |\n")
            f.write("| :--- | :--- |\n")
            # This is a bit simplified, but good for a start
            for t in result.detected_targets[:10]:
                f.write(f"| {t} | - |\n") # Optimization results don't easily expose barriers here
            f.write("\n")

        f.write("## 4. Analytical Metadata\n")
        f.write("### Matrix Explainability\n")
        for k, v in result.matrix_explainability.items():
            f.write(f"- **{k}:** {v}\n")
            
    return output_dir

def generate_comparison_report(
    results: List[FormulationResult],
    conditions_list: List[Dict[str, Any]],
    output_dir: Optional[Path] = None
) -> Path:
    """
    Generates a side-by-side comparison report for multiple formulation results.
    """
    if output_dir is None:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(f"reports/comparison_{timestamp}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Save JSON Comparison
    json_path = output_dir / "comparison.json"
    comparison_data = {
        "timestamp": datetime.datetime.now().isoformat(),
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
                "trapping_efficiency": float(res.trapping_efficiency)
            }
        })
    
    with open(json_path, "w") as f:
        json.dump(comparison_data, f, indent=4)
        
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
        f.write("\n")
        
        f.write("## 2. Competitive Highlights\n")
        # Find winners
        best_target = max(results, key=lambda x: x.target_score)
        best_safety = min(results, key=lambda x: x.safety_score)
        best_risk = min(results, key=lambda x: x.off_flavour_risk)
        
        f.write(f"- 🏆 **Highest Target Score:** {best_target.name} ({best_target.target_score:.2f})\n")
        f.write(f"- 🛡️ **Safest Formulation:** {best_safety.name} ({best_safety.safety_score:.2f})\n")
        f.write(f"- 🍃 **Lowest Off-Flavour Risk:** {best_risk.name} ({best_risk.off_flavour_risk:.2f})\n\n")

        f.write("## 3. Decision Summaries\n")
        import io
        from contextlib import redirect_stdout
        for res in results:
            f.write(f"### {res.name}\n")
            f.write("```text\n")
            with io.StringIO() as buf, redirect_stdout(buf):
                # We don't have the warnings easily here without passing them in, 
                # but for comparison we can show the boxes.
                # Assuming empty warnings for simplicity in comparison view if not provided.
                render_decision_summary_cli(res, [])
                f.write(buf.getvalue())
            f.write("```\n\n")
            
    return output_dir
