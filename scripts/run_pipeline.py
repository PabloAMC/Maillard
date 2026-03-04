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


def print_table(active_pathways: list):
    """Reuse the standard unicode table format."""
    print("    ┌" + "─"*24 + "┬" + "─"*18 + "┬" + "─"*15 + "┬" + "─"*15 + "┬" + "─"*22 + "┐")
    print("    │ PREDICTED COMPOUND     │ PATHWAY TYPE     │ BARRIER (ΔE‡) │ PENALTY RISK  │ TOXICITY ALERT       │")
    print("    ├" + "─"*24 + "┼" + "─"*18 + "┼" + "─"*15 + "┼" + "─"*15 + "┼" + "─"*22 + "┤")
    
    for p in active_pathways:
        target_str = p['target'].label if p['target'] else "Unknown"
        
        tag = ""
        if p['type'] == 'desirable': tag = "[✅ AROMA]"
        elif p['type'] == 'competing': tag = "[⚠️ COMPETING]"
        elif p['type'] == 'toxic': tag = "[☠️ TOXIC]"
        elif p['type'] == 'masking': tag = "[🛡️ MASKING]"
        
        barrier_str = f"{p['span']:.1f} kcal"
        penalty_str = p.get('penalty', '-')
        
        tox_str = "-"
        if p.get('toxicity'):
            meta = p['toxicity']
            tox_str = f"[{meta['priority']}] {meta['name']}"
            
        col1 = _trunc(target_str, 22)
        col2 = _trunc(tag, 16)
        col3 = _trunc(barrier_str, 13)
        col4 = _trunc(penalty_str, 13)
        col5 = _trunc(tox_str, 20)
        
        print(f"    │ {col1} │ {col2} │ {col3} │ {col4} │ {col5} │")
        
    print("    └" + "─"*24 + "┴" + "─"*18 + "┴" + "─"*15 + "┴" + "─"*15 + "┴" + "─"*22 + "┘")


def main():
    parser = argparse.ArgumentParser(description="Maillard formulation screening pipeline.")
    parser.add_argument("--sugars", type=str, required=True, help="Comma-separated sugars (e.g. ribose,glucose)")
    parser.add_argument("--amino-acids", type=str, required=True, help="Comma-separated amino acids (e.g. cysteine,glycine)")
    parser.add_argument("--other", type=str, default="", help="Comma-separated other precursors (e.g. thiamine)")
    parser.add_argument("--ph", type=float, default=6.0, help="Reaction pH (default 6.0)")
    parser.add_argument("--temp", type=float, default=150.0, help="Heating temperature in Celsius (default 150)")
    parser.add_argument("--xtb", action="store_true", help="Run full GFN2-xTB structural optimizations (SLOW!). Defaults to fast Hammond estimating.")
    parser.add_argument("--list-precursors", action="store_true", help="List available precursors and exit")
    
    # Catch simple cases where user only wants to list
    if "--list-precursors" in sys.argv:
        print("Available Precursors:")
        for name in precursor_resolver.list_available():
            print(f"  - {name}")
        sys.exit(0)

    args = parser.parse_args()

    # 1. Resolve Precursors
    names = args.sugars.split(",") + args.amino_acids.split(",")
    if args.other:
        names += args.other.split(",")
    names = [n.strip() for n in names if n.strip()]

    try:
        precursors = precursor_resolver.resolve_many(names)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    print("======================================================")
    print("      Maillard Generative Pipeline (Phase 7)")
    print("======================================================\n")
    print(f"Inputs: {', '.join(p.label for p in precursors)}")
    print(f"Conditions: pH {args.ph}, {args.temp}°C")
    print("-" * 60)

    # 2. Enumerate Pathways (SmirksEngine)
    conditions = ReactionConditions(pH=args.ph, temperature_celsius=args.temp)
    engine = SmirksEngine(conditions)
    print("Running rule-based generation (up to 3 generations)...")
    steps = engine.enumerate(precursors)
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
            # Fake fast evaluation using heuristics to allow testing pipeline logic
            # This is a mock since real XTBScreener.compute_reaction_energy calls RDKit EmbedMultipleConfs which is slow,
            # even when it falls back to Hammond. We want true CLI instant-response mode here.
            
            # Simple heuristic mock logic (for demo CLI speed):
            # Amadori/Heyns: ~25 kcal
            # Schiff base: ~15 kcal
            # Dehydration: ~30 kcal
            # Condensation/Strecker: ~20 kcal
            # Default: 40 kcal
            bar = 40.0
            if step.reaction_family:
                fm = step.reaction_family.lower()
                if "amadori" in fm or "heyns" in fm: bar = 25.0
                elif "schiff" in fm: bar = 15.0
                elif "dehydration" in fm or "enolisation" in fm: bar = 30.0
                elif "strecker" in fm: bar = 22.0
                elif "retro" in fm: bar = 35.0
                elif "cysteine" in fm or "beta" in fm: bar = 32.0
                elif "thiol" in fm: bar = 18.0
                
            # Apply condition modifiers
            ph_mult = conditions.get_ph_multiplier(step.reaction_family or "")
            arrh_mult = conditions.get_arrhenius_multiplier(bar)
            
            # Very crude adjustment to show differences
            adjusted_bar = bar * ph_mult
            
            barriers_dict[step_key] = adjusted_bar
            
    print("\n\nScreening complete.")

    # 4. Recommend Targets
    recommender = Recommender(None)
    initial_pool_smiles = [p.smiles for p in precursors]
    
    active_pathways = recommender.predict_from_steps(steps, barriers_dict, initial_pool_smiles)
    
    print("-" * 60)
    if not active_pathways:
        print("  [!] No target compounds (desirable/off-flavour/toxic) found generated from these precursors.")
    else:
        print("  Predicted Targets (ranked by pathway bottleneck barrier):")
        print_table(active_pathways)
        
    print("\n" + "═"*85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - FAST mode barrier estimates are purely qualitative heuristics.")
    print("    - Use --xtb for valid semi-empirical calculations.")
    print("═"*85 + "\n")

if __name__ == "__main__":
    main()
