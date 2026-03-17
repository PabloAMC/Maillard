#!/usr/bin/env python3
"""
scripts/optimize_formulation.py

Runs Bayesian optimization to find the Pareto-optimal Maillard formulation 
balancing target flavor traits against off-flavours and safety risks.
"""

import sys
import argparse
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.bayesian_optimizer import FormulationOptimizer  # noqa: E402
from src.usability_reports import (
    DomainOfValidityChecker, 
    render_domain_warnings_cli, 
    render_decision_summary_cli,
    render_deep_explainability_cli
)  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description="Bayesian Formulation Optimizer")
    parser.add_argument("--sugars", type=str, required=True, help="Comma-separated list of sugars")
    parser.add_argument("--amino-acids", type=str, required=True, help="Comma-separated list of amino acids")
    parser.add_argument("--lipids", type=str, default="", help="Comma-separated list of lipids/off-flavor precursors")
    parser.add_argument("--target-tag", type=str, default="meaty", help="Target flavor profile (default: meaty)")
    parser.add_argument("--minimize-tag", type=str, default="beany", help="Flavor to suppress (default: beany)")
    parser.add_argument("--n-iterations", type=int, default=50, help="Number of Optuna trials (default: 50)")
    parser.add_argument("--risk-aversion", type=float, default=1.0, help="Penalty weight for toxic markers (default: 1.0)")
    parser.add_argument("--pre-process", type=str, choices=["none", "yeast", "protease", "both"], default="none", help="Biological pre-processing step")
    parser.add_argument("--protein-type", choices=["free", "pea_conc", "pea_iso", "soy_conc", "soy_iso", "myco"], default="free", help="Protein matrix type (default: free)")
    parser.add_argument("--denaturation-state", type=float, default=None, help="Optional denaturation state 0-1. If omitted, infer from temperature/time/pH.")
    parser.add_argument("--report", action="store_true", help="Generate consolidated JSON/Markdown report for the best trial")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the report")
    args = parser.parse_args()

    sugars = [s.strip() for s in args.sugars.split(",")] if args.sugars else []
    aas = [a.strip() for a in args.amino_acids.split(",")] if args.amino_acids else []
    lipids = [lip.strip() for lip in args.lipids.split(",")] if args.lipids else []

    print("======================================================")
    print("      Maillard Bayesian Formulation Optimizer")
    print("======================================================\n")
    print(f"Target:       Maximize {args.target_tag.upper()}")
    print(f"Constraint:   Minimize {args.minimize_tag.upper()} and SAFETY RISKS")
    print(f"Risk Weight:  {args.risk_aversion}")
    print(f"Iterations:   {args.n_iterations}")
    print("-" * 54)
    print(f"System:       {'+'.join(sugars + aas + lipids)}")
    print(f"Pre-process:  {args.pre_process}")
    print("-" * 54)
    print("Starting optimization... this may take a few minutes depending on complexity.\n")

    optimizer = FormulationOptimizer(
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
        risk_aversion=args.risk_aversion,
        protein_type=args.protein_type,
        denaturation_state=args.denaturation_state
    )
    
    study = optimizer.optimize(
        fixed_sugars=sugars,
        fixed_amino_acids=aas,
        fixed_lipids=lipids,
        n_trials=args.n_iterations
    )
    
    print("\n======================================================")
    print("Optimization Complete!")
    print(f"Best Objective Score: {study.best_value:.4f}")
    print("======================================================")
    
    best_trial = study.best_trial
    
    print("\n⭐ Optimized Parameter Set:")
    for key, value in best_trial.params.items():
        if "conc" in key:
            print(f"  {key.ljust(15)} : {value:.3f} M")
        elif key == "ph":
            print(f"  {key.ljust(15)} : {value:.2f}")
        elif key == "temp":
            print(f"  {key.ljust(15)} : {value:.1f} °C")
        elif key == "aw":
            print(f"  {key.ljust(15)} : {value:.2f}")
        else:
            if isinstance(value, (int, float)):
                print(f"  {key.ljust(15)} : {value:.1f} min")
            else:
                print(f"  {key.ljust(15)} : {value}")
            
    print("\n🔍 Predicted Outcomes:")
    print(f"  Target Score   : {best_trial.user_attrs.get('target_score', 0):.4f}")
    safe_score = best_trial.user_attrs.get('safety_score', 0)
    print(f"  Safety Risk    : {safe_score:.4f} " + ("⚠️" if safe_score > 0 else "✅"))
    print(f"  Off-Flavours   : {best_trial.user_attrs.get('off_flavour_risk', 0):.4f}")
    flagged = best_trial.user_attrs.get('flagged_toxics', [])
    if flagged:
        print(f"  Flagged Toxins : {', '.join(set(flagged))}")
    
    # --- Decision Summary ---
    print("\n🔍 Re-evaluating best formulation for detailed summary...")
    # Re-construct formulation for evaluate_single
    best_formulation = {
        "name": "Best_Trial",
        "sugars": sugars,
        "amino_acids": aas,
        "lipids": lipids,
        "molar_ratios": best_trial.params, # This might need some mapping if names differ
        "ph": best_trial.params.get('ph', 6.0),
        "temp": best_trial.params.get('temp', 150.0),
        "aw": best_trial.params.get('aw', 0.8),
        "time_minutes": best_trial.params.get('time_minutes', 60.0),
        "protein_type": args.protein_type,
        "denaturation_state": args.denaturation_state
    }
    # InverseDesigner is needed here
    from src.inverse_design import InverseDesigner
    designer = InverseDesigner(args.target_tag, args.minimize_tag)
    from src.smirks_engine import ReactionConditions
    cond = ReactionConditions(
        pH=best_formulation["ph"],
        temperature_celsius=best_formulation["temp"],
        water_activity=best_formulation["aw"],
        protein_type=args.protein_type
    )
    res = designer.evaluate_single(best_formulation, cond)
    
    checker = DomainOfValidityChecker(args.target_tag)
    warnings = checker.check(
        precursor_names=sugars + aas + lipids, 
        protein_type=args.protein_type,
        temp_c=res.matrix_explainability.get("temperature", 150.0), # or just from params
        ph=res.matrix_explainability.get("pH", 6.0)
    )
    render_decision_summary_cli(res, warnings)
    render_deep_explainability_cli(res)
    
    if args.report:
        from src.reporting import generate_report
        report_dir = generate_report(
            res, 
            warnings, 
            vars(args), 
            output_dir=Path(args.output_dir) if args.output_dir else None
        )
        print(f"📄 Best Trial Report generated in: {report_dir}")
    
    print("\n(Note: Optimization trace saved in memory. For persistent tracking, configure Optuna DB storage.)")
    print("======================================================")

if __name__ == "__main__":
    main()
