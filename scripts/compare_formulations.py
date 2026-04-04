#!/usr/bin/env python3
"""
scripts/compare_formulations.py

Compare two or more formulations from the grid under specific conditions.
Generates a consolidated side-by-side report.
"""

import sys
import argparse
from pathlib import Path
from typing import List

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.pipeline import MaillardPipeline
from src.conditions import ReactionConditions
from src.reporting import generate_comparison_report
from src.usability_reports import DomainOfValidityChecker, build_confidence_package

def main():
    parser = argparse.ArgumentParser(description="Compare multiple Maillard formulations.")
    parser.add_argument("--names", type=str, required=True, help="Comma-separated formulation names to compare")
    parser.add_argument("--ph", type=float, default=6.0, help="Reaction pH (default 6.0)")
    parser.add_argument("--temp", type=float, default=150.0, help="Temperature in Celsius (default 150)")
    parser.add_argument("--aw", "--water-activity", type=float, default=1.0, help="Water activity (default 1.0)")
    parser.add_argument("--target-tag", type=str, default="meaty", help="Target flavor profile (default: meaty)")
    parser.add_argument("--minimize-tag", type=str, default="beany", help="Flavor to minimize (default: beany)")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the comparison report")
    
    args = parser.parse_args()
    
    names = [n.strip() for n in args.names.split(",")]
    if len(names) < 2:
        print("ERROR: Please specify at least two formulation names to compare.")
        sys.exit(1)
    
    print("======================================================")
    print("      Maillard Formulation Comparison Tool")
    print("======================================================\n")
    print(f"Comparing: {', '.join(names)}")
    print(f"Target:    {args.target_tag}")
    print(f"Minimize:  {args.minimize_tag}")
    print(f"Conditions: pH {args.ph}, {args.temp}°C, aᵥ {args.aw}")
    print("-" * 54)
    
    designer = MaillardPipeline(args.target_tag, args.minimize_tag)
    
    # Filter grid for requested names
    requested_forms = []
    for f_name in names:
        match = next((f for f in designer.grid if f.get("name") == f_name), None)
        if not match:
            print(f"WARNING: Formulation '{f_name}' not found in grid. Skipping.")
            continue
        requested_forms.append(match)
        
    if not requested_forms:
        print("ERROR: No valid formulations found to compare.")
        sys.exit(1)
        
    conditions = ReactionConditions(
        pH=args.ph,
        temperature_celsius=args.temp,
        water_activity=args.aw
    )
    
    print(f"Evaluating {len(requested_forms)} formulations...")
    results = designer.evaluate_all(conditions, grid_override=requested_forms)
    results_by_name = {result.name: result for result in results}
    ordered_results = [results_by_name[formulation["name"]] for formulation in requested_forms if formulation["name"] in results_by_name]
    checker = DomainOfValidityChecker(args.target_tag)
    warnings_list = []
    for formulation in requested_forms:
        result = results_by_name.get(formulation["name"])
        if result is None:
            continue
        precursor_names = []
        for key in ("sugars", "amino_acids", "additives", "lipids"):
            precursor_names.extend(formulation.get(key, []))
        protein_type = formulation.get("protein_type", "free")
        item_warnings = checker.check(
            precursor_names=precursor_names,
            protein_type=protein_type,
            temp_c=float(formulation.get("temp", args.temp)),
            ph=float(formulation.get("ph", args.ph)),
            aw=float(formulation.get("aw", args.aw)),
            matrix_explainability=result.matrix_explainability,
        )
        result.confidence_metadata = build_confidence_package(
            result,
            item_warnings,
            precursor_names=precursor_names,
            protein_type=protein_type,
            formulation=formulation,
            baseline_conditions=conditions,
            designer=designer,
        )
        warnings_list.append(item_warnings)
    
    # Prepare conditions_list for the report
    conditions_list = [vars(args)] * len(ordered_results)
    
    report_dir = generate_comparison_report(
        ordered_results, 
        conditions_list, 
        warnings_list=warnings_list,
        output_dir=Path(args.output_dir) if args.output_dir else None
    )
    
    print(f"\n✅ Comparison report generated in: {report_dir}")
    print(f"   - {report_dir / 'comparison.md'}")
    print(f"   - {report_dir / 'comparison.json'}")
    print("\n" + "=" * 54)

if __name__ == "__main__":
    main()
