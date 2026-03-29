#!/usr/bin/env python3
"""
scripts/scientist_quickstart.py — High-level scientist entrypoint for the Maillard Framework.
Answers the 3 core questions:
1. What should I try next?
2. Why does the model think that?
3. What is benchmarked vs. model-derived?
"""

import sys
import os
from pathlib import Path

# Add src to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

# Imports will be handled inside main to ensure sys.path is set

def main():
    if len(sys.argv) < 2:
        print("\n=== Maillard Framework Scientist Quickstart ===\n")
        print("Usage: python3 scripts/scientist_quickstart.py <formulation_name>")
        print("Example: python3 scripts/scientist_quickstart.py pea_meaty")
        print("\nAvailable formulations (from formulation_grid.yml):")
        # Quick peek at grid
        grid_path = ROOT / "data" / "formulation_grid.yml"
        if grid_path.exists():
            import yaml
            with open(grid_path) as f:
                grid = yaml.safe_load(f)
                for name in list(grid.keys())[:5]:
                    print(f"  - {name}")
        print("\n" + "="*45 + "\n")
        sys.exit(0)

    name = sys.argv[1]
    print(f"Evaluating formulation: {name} ...")
    
    try:
        from src.pipeline import MaillardPipeline
        from src.conditions import ReactionConditions
        from src.reporting import build_family_role_explanation
    except ImportError as e:
        print(f"Error: Could not import Maillard source modules. {e}")
        sys.exit(1)

    name = sys.argv[1]
    
    try:
        # 1. Initialize pipeline
        pipeline = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
        
        # 2. Find formulation in grid
        formulation = None
        for f in pipeline.grid:
            if f.get("name") == name:
                formulation = f
                break
        
        if not formulation:
            print(f"Error: Formulation '{name}' not found in data/formulation_grid.yml")
            sys.exit(1)

        # 3. Run evaluation
        print(f"Evaluating formulation: {name} ...")
        conditions = ReactionConditions()
        result = pipeline.evaluate_single(formulation, conditions)
        
        # 4. Get family explanation
        explanation = build_family_role_explanation(result)
        summary = explanation.get("summary", {})
        
        # 3. Answer the 3 Core Questions
        print(f"\n=== RESULTS FOR: {name} ===\n")
        
        # Q1: What should I try next?
        print("1. What should I try next?")
        missing = explanation.get("missing_or_transferred", [])
        if missing:
            top_missing = missing[0]
            print(f"   -> CONSIDER ACTIVATING: {top_missing['display_name']} ({top_missing['slr_family']})")
            print(f"      Reason: This family has {top_missing['primary_payload_count']} benchmark anchors but is currently inactive.")
        else:
            print("   -> FULL COVERAGE: All 16 canonical families are represented.")
            print("      Optimization: Review 'modifier' lanes in the full report to refine precursor ratios.")
        
        # Q2: Why does the model think that?
        print("\n2. Why does the model think that?")
        drivers = explanation.get("drivers", [])
        if drivers:
            print(f"   -> TOP DRIVER: {drivers[0]['display_name']} ({drivers[0]['slr_family']})")
            print(f"      Delta: {drivers[0]['target_score_delta']:+.3f} (target) | {drivers[0]['off_flavour_risk_delta']:+.3f} (off-note)")
            print(f"      Insight: {drivers[0]['summary']}")
        else:
            print("   -> No strong kinetic drivers detected. Result is dominated by matrix-scope defaults.")

        # Q3: Trust Boundary
        print("\n3. What is benchmark-backed vs. model-derived?")
        active = summary.get("active_lane_count", 0)
        total = summary.get("total_canonical_family_count", 16)
        print(f"   -> {active}/{total} canonical families are active.")
        print(f"   -> {summary.get('driver_count', 0)} lanes are directly shaping the result via benchmark-linked priors.")
        print(f"   -> (Refer to 'report.md' in the results folder for the full evidence ladder)")

        print("\n" + "="*45)
        print(f"Detailed artifacts generated in: results/{name}/")
        print("="*45 + "\n")

    except Exception as e:
        print(f"Error during evaluation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
