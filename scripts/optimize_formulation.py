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

def main():
    parser = argparse.ArgumentParser(description="Bayesian Formulation Optimizer")
    parser.add_argument("--sugars", type=str, required=True, help="Comma-separated list of sugars")
    parser.add_argument("--amino-acids", type=str, required=True, help="Comma-separated list of amino acids")
    parser.add_argument("--lipids", type=str, default="", help="Comma-separated list of lipids/off-flavor precursors")
    parser.add_argument("--target-tag", type=str, default="meaty", help="Target flavor profile (default: meaty)")
    parser.add_argument("--minimize-tag", type=str, default="beany", help="Flavor to suppress (default: beany)")
    parser.add_argument("--n-iterations", type=int, default=50, help="Number of Optuna trials (default: 50)")
    parser.add_argument("--risk-aversion", type=float, default=1.0, help="Penalty weight for toxic markers (default: 1.0)")
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
    print("-" * 54)
    print("Starting optimization... this may take a few minutes depending on complexity.\n")

    optimizer = FormulationOptimizer(
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
        risk_aversion=args.risk_aversion
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
            print(f"  {key.ljust(15)} : {value:.1f} min")
            
    print("\n🔍 Predicted Outcomes:")
    print(f"  Target Score   : {best_trial.user_attrs.get('target_score', 0):.4f}")
    safe_score = best_trial.user_attrs.get('safety_score', 0)
    print(f"  Safety Risk    : {safe_score:.4f} " + ("⚠️" if safe_score > 0 else "✅"))
    print(f"  Off-Flavours   : {best_trial.user_attrs.get('off_flavour_risk', 0):.4f}")
    flagged = best_trial.user_attrs.get('flagged_toxics', [])
    if flagged:
        print(f"  Flagged Toxins : {', '.join(set(flagged))}")
    
    print("\n(Note: Optimization trace saved in memory. For persistent tracking, configure Optuna DB storage.)")
    print("======================================================")

if __name__ == "__main__":
    main()
