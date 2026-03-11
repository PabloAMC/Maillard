#!/usr/bin/env python3
"""
scripts/run_curated_screening.py — Run xTB screening on hand-curated Maillard pathways.

Feeds the 5 curated Maillard cascades (A–E) through the xTB pipeline:
  curated_pathways → PathwayRanker → ranked barriers → JSON output
"""

import sys
import json
from pathlib import Path

# Project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.reactions.curated_pathways import PATHWAYS # noqa: E402
from src.pathway_ranker import PathwayRanker  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402


def main():
    # Typical Maillard cooking conditions
    conditions = ReactionConditions(
        pH=6.0,
        temperature_celsius=150,
        water_activity=0.7,
    )
    
    print("=" * 70)
    print("Maillard Curated Pathway Screening (xTB Tier 1)")
    print(f"Conditions: pH={conditions.pH}, T={conditions.temperature_celsius}°C, a_w={conditions.water_activity}")
    print("=" * 70)
    
    # List pathways
    print(f"\n{len(PATHWAYS)} pathways loaded:")
    for name, steps in PATHWAYS.items():
        print(f"  {name}: {len(steps)} steps")
    
    # Run xTB screening
    print("\nRunning xTB screening (this may take a few minutes)...\n")
    ranker = PathwayRanker(n_cores=1, conditions=conditions)  # n_cores=1 for stability
    profiles = ranker.screen_pathways(PATHWAYS)
    
    # Print ranked results
    print("\n" + "=" * 70)
    print("RANKED PATHWAYS (fastest → slowest)")
    print("=" * 70)
    for rank, profile in enumerate(profiles, 1):
        print(f"\n#{rank}: {profile}")
        for i, step in enumerate(profile.steps):
            barrier = profile.barrier_kcal_list[i] if i < len(profile.barrier_kcal_list) else "N/A"
            delta_e = profile.deltaE_kcal_list[i] if i < len(profile.deltaE_kcal_list) else "N/A"
            print(f"    Step {i+1}: {step}")
            if isinstance(barrier, float):
                print(f"           ΔE‡ = {barrier:.1f} kcal/mol, ΔE = {delta_e:.1f} kcal/mol")
    
    # Save results to JSON
    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    output_path = results_dir / "curated_screening_results.json"
    
    results = []
    for profile in profiles:
        results.append({
            "pathway": profile.pathway_name,
            "max_barrier_kcal": profile.rate_limiting_barrier,
            "energetic_span_kcal": profile.energetic_span,
            "overall_deltaE_kcal": profile.overall_thermodynamics,
            "steps": [
                {
                    "reaction": str(step),
                    "barrier_kcal": profile.barrier_kcal_list[i],
                    "deltaE_kcal": profile.deltaE_kcal_list[i],
                }
                for i, step in enumerate(profile.steps)
            ],
        })
    
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
