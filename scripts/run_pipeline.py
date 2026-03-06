#!/usr/bin/env python3
"""
scripts/run_pipeline.py — End-to-End Generative Maillard Pipeline

The core Phase 7 orchestrator.
1. Accepts custom formulations via CLI.
2. Generates the full reaction network via SmirksEngine.
3. Screens barriers (with quick Hammond fallback by default).
4. Recommends target volatile outputs.
"""

import sys
import argparse
from pathlib import Path
from typing import List, Dict

# Add project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.conditions import ReactionConditions
from src.smirks_engine import SmirksEngine
from src.pathway_extractor import ElementaryStep
from src.xtb_screener import XTBScreener
from src.recommend import Recommender, _trunc
from src import precursor_resolver
from src.barrier_constants import get_barrier, HEME_CATALYST_FAMILIES, HEME_CATALYST_REDUCTION


def print_table(active_pathways: list):
    """Update table to include sensory descriptors."""
    # Column widths: 22, 16, 12, 30, 20
    print("    ┌" + "─"*22 + "┬" + "─"*16 + "┬" + "─"*12 + "┬" + "─"*30 + "┬" + "─"*20 + "┐")
    print("    │ PREDICTED COMPOUND   │ FORMATION TAG    │ BARRIER    │ SENSORY CHARACTER            │ TOXICITY/RISK        │")
    print("    ├" + "─"*22 + "┼" + "─"*16 + "┼" + "─"*12 + "┼" + "─"*30 + "┼" + "─"*20 + "┤")
    
    for p in active_pathways:
        target_str = p['target'].label if p['target'] else "Unknown"
        
        tag = ""
        if p['type'] == 'desirable': tag = "[✅ AROMA]"
        elif p['type'] == 'competing': tag = "[⚠️ COMPETING]"
        elif p['type'] == 'toxic': tag = "[☠️ TOXIC]"
        elif p['type'] == 'masking': tag = "[🛡️ MASKING]"
        
        barrier_str = f"{p['span']:.1f} kcal"
        sensory_str = p.get('sensory', '-')
        
        tox_str = "-"
        if p.get('toxicity'):
            meta = p['toxicity']
            tox_str = f"[{meta['priority']}] {meta['name']}"
            
        col1 = _trunc(target_str, 20)
        col2 = _trunc(tag, 14)
        col3 = _trunc(barrier_str, 10)
        col4 = _trunc(sensory_str, 28)
        col5 = _trunc(tox_str, 18)
        
        print(f"    │ {col1} │ {col2} │ {col3} │ {col4} │ {col5} │")
        
    print("    └" + "─"*22 + "┴" + "─"*16 + "┴" + "─"*12 + "┴" + "─"*30 + "┴" + "─"*20 + "┘")


def main():
    parser = argparse.ArgumentParser(description="Maillard formulation screening pipeline.")
    parser.add_argument("--sugars", type=str, default="", help="Comma-separated sugars (e.g. ribose,glucose)")
    parser.add_argument("--amino-acids", type=str, default="", help="Comma-separated amino acids (e.g. cysteine,glycine)")
    parser.add_argument("--additives", type=str, default="", help="Comma-separated additives (e.g. thiamine,glutathione)")
    parser.add_argument("--lipids", type=str, default="", help="Comma-separated lipid aldehydes (e.g. hexanal,nonanal)")
    parser.add_argument("--ph", type=float, default=6.0, help="Reaction pH (default 6.0)")
    parser.add_argument("--temp", type=float, default=150.0, help="Heating temperature in Celsius (default 150)")
    parser.add_argument("--ratios", type=str, default="", help="Optional comma-separated molar ratios (e.g. cysteine:2.0,ribose:1.0). Default 1.0.")
    parser.add_argument("--catalyst", type=str, default=None, choices=["heme"], help="Apply catalyst effect (e.g. heme)")
    parser.add_argument("--aw", "--water-activity", type=float, default=1.0, help="Water activity (default 1.0)")
    parser.add_argument("--target", type=str, default=None, help="Inverse design target sensory tag (e.g. meaty, roasted)")
    parser.add_argument("--minimize", type=str, default="beany", help="Inverse design off-flavour tag to minimize (default: beany)")
    parser.add_argument("--xtb", action="store_true", help="Run full GFN2-xTB structural optimizations (SLOW!). Defaults to fast Hammond estimating.")
    parser.add_argument("--list-precursors", action="store_true", help="List available precursors and exit")
    
    # Catch simple cases where user only wants to list
    if "--list-precursors" in sys.argv:
        print("Available Precursors:")
        for name in precursor_resolver.list_available():
            print(f"  - {name}")
        sys.exit(0)

    args = parser.parse_args()

    # 1. Setup Conditions
    conditions = ReactionConditions(
        pH=args.ph, 
        temperature_celsius=args.temp,
        water_activity=args.aw
    )

    # =================================================================
    # Inverse Design Mode Execution
    # =================================================================
    if args.target:
        print("======================================================")
        print("      Maillard Inverse Design Mode (Phase 7.5)")
        print("======================================================\n")
        print(f"Targeting Sensory Profile: '{args.target}'")
        print(f"Minimizing Risk Profile:   '{args.minimize}'")
        print(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}")
        print("-" * 60)
        
        try:
            from src.inverse_design import InverseDesigner
            designer = InverseDesigner(args.target, args.minimize)
            print(f"Evaluating {len(designer.grid)} industrial formulations against tags...")
            
            results = designer.evaluate_all(conditions)
            
            print("\n  Top Recommended Formulations:")
            print("    ┌──────────────────────────────┬──────────────┬──────────────┬─────────────┬─────────────┐")
            print("    │ FORMULATION                  │ TARGET SCORE │ RISK PENALTY │ LYS BUDGET  │ TRAP EFF    │")
            print("    ├──────────────────────────────┼──────────────┼──────────────┼─────────────┼─────────────┤")
            
            for res in results[:10]:  # Top 10
                t_score = f"{res.target_score:.1f}"
                r_score = f"{res.off_flavour_risk:.1f}"
                lys = f"{res.lysine_budget:.1f}%"
                trap = f"{res.trapping_efficiency:.1f}%"
                print(f"    │ {res.name[:28]:<28} │ {t_score:>12} │ {r_score:>12} │ {lys:>11} │ {trap:>11} │")
            
            print("    └──────────────────────────────┴──────────────┴──────────────┴─────────────┴─────────────┘")
            sys.exit(0)
            
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # =================================================================
    # Standard Forward Pipeline Execution
    # =================================================================

    # 1. Resolve Precursors
    names = []
    if args.sugars: names += args.sugars.split(",")
    if args.amino_acids: names += args.amino_acids.split(",")
    if args.additives: names += args.additives.split(",")
    if args.lipids: names += args.lipids.split(",")
    names = [n.strip() for n in names if n.strip()]

    if not names:
        print("ERROR: No precursors specified. Use --sugars, --amino-acids, --additives, --lipids OR --target.")
        sys.exit(1)

    try:
        precursors = precursor_resolver.resolve_many(names)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # Parse ratios
    ratio_dict = {}
    if args.ratios:
        for pair in args.ratios.split(","):
            if ":" in pair:
                k, v = pair.split(":")
                try:
                    ratio_dict[k.strip().lower()] = float(v.strip())
                except ValueError:
                    print(f"ERROR: Invalid ratio value '{v}' for '{k}'. Must be a float.")
                    sys.exit(1)

    # Build canonical smiles to concentration map
    from src.recommend import _canon
    initial_concentrations = {}
    for p in precursors:
        qty = ratio_dict.get(p.label.lower(), 1.0)
        initial_concentrations[_canon(p.smiles)] = qty

    # Print Forward Pipeline Header
    print("======================================================")
    print("      Maillard Generative Pipeline (Phase 7)")
    print("======================================================\n")
    
    print(f"Inputs: {', '.join(p.label for p in precursors)}")
    if args.ratios:
        print(f"Molar Ratios: {', '.join(f'{k}: {v}' for k, v in ratio_dict.items())}")
    print(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}, Catalyst: {args.catalyst or 'None'}")
    print("-" * 60)

    # 2. Enumerate Pathways (SmirksEngine)
    engine = SmirksEngine(conditions)
    print("Running rule-based generation (up to 4 generations)...")
    steps = engine.enumerate(precursors, max_generations=4)
    print(f"Generated {len(steps)} unique elementary steps in the reaction network.")

    # 3. Screen Barriers
    barriers_dict: Dict[str, float] = {}
    screener = XTBScreener()

    print(f"\nEvaluating thermodynamics for {len(steps)} steps...")
    if not args.xtb:
        print("  Using FAST mode (Hammond fallbacks without full structural optimization).")
        print("  Pass --xtb for rigorous GFN2-xTB geometries and NEB approximations.")
        
    for i, step in enumerate(steps):
        sys.stdout.write(f"\r  Evaluating step {i+1}/{len(steps)}...")
        sys.stdout.flush()
        
        step_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        
        if args.xtb:
            try:
                # Full execution
                dE, bar = screener.compute_reaction_energy(step)
            except Exception as e:
                dE, bar = 99.0, 99.0 # Extreme penalty for failed convergence
            barriers_dict[step_key] = bar
        else:
            # Calculate boundaries from shared constants
            bar = get_barrier(step.reaction_family)
                
            # Apply condition modifiers
            ph_mult = conditions.get_ph_multiplier(step.reaction_family or "")
            
            # Accelerate by dividing barrier if multiplier > 1 (higher kinetic probability)
            adjusted_bar = bar / ph_mult

            # 3a. Apply Catalyst Heme Heuristic 
            if args.catalyst == "heme" and step.reaction_family in HEME_CATALYST_FAMILIES:
                adjusted_bar -= HEME_CATALYST_REDUCTION 
                
            barriers_dict[step_key] = adjusted_bar
            
    print("\n\nScreening complete.")

    # 4. Recommend Targets
    recommender = Recommender(None)
    
    results = recommender.predict_from_steps(steps, barriers_dict, initial_concentrations)
    active_pathways = results["targets"]
    metrics = results["metrics"]
    
    print("-" * 60)
    if not active_pathways:
        print("  [!] No target compounds (desirable/off-flavour/toxic) found generated from these precursors.")
    else:
        print("  Predicted Targets (ranked by pathway bottleneck barrier):")
        print_table(active_pathways)

    # 5. Display PBMA Metrics
    if metrics["trapping_efficiency"] or metrics.get("lysine_budget_dha", 0) > 0:
        print("\n  PBMA Formulation Metrics:")
        for lipid, eff in metrics["trapping_efficiency"].items():
            bar_len = int(eff / 5)
            bar_str = "█" * bar_len + "░" * (20 - bar_len)
            print(f"    │ Lipid Trapping ({lipid}): {eff:5.1f}% | {bar_str} |")
        
        lys_budget = metrics.get("lysine_budget_dha", 0)
        if lys_budget > 0:
            bar_len = int(lys_budget / 5)
            bar_str = "█" * bar_len + "░" * (20 - bar_len)
            print(f"    │ Lysine Budget (DHA):    {lys_budget:5.1f}% | {bar_str} |")
            if lys_budget > 50.0:
                print("    ⚠️  WARNING: High Lysine consumption by DHA pathway significantly reduces aroma yield.")
        
    print("\n" + "═"*85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - FAST mode barrier estimates are purely qualitative heuristics.")
    print("    - Use --xtb for valid semi-empirical calculations.")
    print("═"*85 + "\n")

if __name__ == "__main__":
    main()
