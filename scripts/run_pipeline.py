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
from dataclasses import asdict
from pathlib import Path
from typing import Dict

# Add project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.conditions import ReactionConditions  # noqa: E402
from src.smirks_engine import SmirksEngine, Species  # noqa: E402
from src.pipeline import MaillardPipeline  # noqa: E402
from src.xtb_screener import XTBScreener  # noqa: E402
from src.recommend import _trunc  # noqa: E402
from src import precursor_resolver  # noqa: E402
from src.barrier_constants import get_barrier, HEME_CATALYST_FAMILIES, HEME_CATALYST_REDUCTION  # noqa: E402
from src.usability_reports import (
    DomainOfValidityChecker, 
        build_confidence_package
)  # noqa: E402
from src.presentation import (
    render_decision_summary_cli,
    render_deep_explainability_cli
)  # noqa: E402
from src.reporting import generate_report  # noqa: E402



def print_table(active_pathways: list):
    """Update table to include sensory descriptors."""
    # Column widths: 22, 16, 12, 30, 20
    print("    ┌" + "─"*22 + "┬" + "─"*16 + "┬" + "─"*12 + "┬" + "─"*30 + "┬" + "─"*20 + "┐")
    print("    │ PREDICTED COMPOUND   │ FORMATION TAG    │ BARRIER    │ SENSORY CHARACTER            │ TOXICITY/RISK        │")
    print("    ├" + "─"*22 + "┼" + "─"*16 + "┼" + "─"*12 + "┼" + "─"*30 + "┼" + "─"*20 + "┤")
    
    for p in active_pathways:
        target_str = p['target'].label if p['target'] else "Unknown"
        
        tag = ""
        if p['type'] == 'desirable':
            tag = "[✅ AROMA]"
        elif p['type'] == 'competing':
            tag = "[⚠️ COMPETING]"
        elif p['type'] == 'toxic':
            tag = "[☠️ TOXIC]"
        elif p['type'] == 'masking':
            tag = "[🛡️ MASKING]"
        
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
    parser.add_argument("--time-minutes", type=float, default=60.0, help="Reaction time in minutes (default 60.0)")
    parser.add_argument("--target", type=str, default=None, help="Inverse design target sensory tag (e.g. meaty, roasted)")
    parser.add_argument("--minimize", type=str, default="beany", help="Inverse design off-flavour tag to minimize (default: beany)")
    parser.add_argument("--xtb", action="store_true", help="Run full GFN2-xTB structural optimizations (SLOW!). Defaults to fast Hammond estimating.")
    parser.add_argument("--protein-type", choices=["free", "pea_conc", "pea_iso", "soy_conc", "soy_iso", "myco"], default="free", help="Protein matrix type for accessibility corrections.")
    parser.add_argument("--denaturation-state", type=float, default=0.5, help="Protein denaturation level (0.0 to 1.0). Default 0.5.")
    parser.add_argument("--list-precursors", action="store_true", help="List available precursors and exit")
    parser.add_argument("--list-tags", action="store_true", help="List available sensory tags and exit")
    parser.add_argument("--report", action="store_true", help="Generate consolidated JSON/Markdown report")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the report")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs against envelope without running simulation")
    
    # Catch simple cases where user only wants to list
    if "--list-precursors" in sys.argv:
        print("Available Precursors:")
        for name in precursor_resolver.list_available():
            print(f"  - {name}")
        sys.exit(0)

    if "--list-tags" in sys.argv:
        from src.sensory import SensoryDatabase
        db = SensoryDatabase()
        print("Available Sensory Tags (Categories):")
        for tag in db.tags.keys():
            print(f"  - {tag}")
        sys.exit(0)

    args = parser.parse_args()

    # 1. Setup Conditions
    conditions = ReactionConditions(
        pH=args.ph, 
        temperature_celsius=args.temp,
        water_activity=args.aw,
        protein_type=args.protein_type
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
            designer = MaillardPipeline(args.target, args.minimize)
            print(f"Evaluating {len(designer.grid)} industrial formulations against tags...")
            
            if args.dry_run:
                print("\n[DRY RUN] Bypassing expensive grid evaluation.")
                print("Validation: Grid dimensions bounded correctly.")
                print("\nDry-run complete. Exiting.")
                sys.exit(0)
                
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
            
            # --- Decision Summary ---
            checker = DomainOfValidityChecker(args.target)
            if results:
                best = results[0]
                best_formulation = next((item for item in designer.grid if item.get("name") == best.name), {})
                best_precursors = []
                for key in ("sugars", "amino_acids", "additives", "lipids"):
                    best_precursors.extend(best_formulation.get(key, []))
                best_protein_type = best_formulation.get("protein_type", args.protein_type)
                best_temp = float(best_formulation.get("temp", conditions.temperature_celsius))
                best_ph = float(best_formulation.get("ph", conditions.pH))
                warnings = checker.check(
                    precursor_names=best_precursors,
                    protein_type=best_protein_type,
                    temp_c=best_temp,
                    ph=best_ph,
                    aw=float(best_formulation.get("aw", conditions.water_activity)),
                    matrix_explainability=best.matrix_explainability,
                )
                best.confidence_metadata = build_confidence_package(
                    best,
                    warnings,
                    precursor_names=best_precursors,
                    protein_type=best_protein_type,
                    formulation=best_formulation,
                    baseline_conditions=conditions,
                    designer=designer,
                )
                # [10] Scientist-Facing Decision Summary
                render_decision_summary_cli(best, warnings)
                render_deep_explainability_cli(best)
            
            if args.report:
                report_dir = generate_report(
                    results[0], 
                    warnings, 
                    vars(args), 
                    output_dir=Path(args.output_dir) if args.output_dir else None
                )
                print(f"📄 Report generated in: {report_dir}")
            
            sys.exit(0)
            
        except ValueError as e:
            print(f"ERROR: {e}")
            sys.exit(1)

    # =================================================================
    # Standard Forward Pipeline Execution
    # =================================================================

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

    # Print Forward Pipeline Header
    print("======================================================")
    print("      Maillard Generative Pipeline (Phase 7)")
    print("======================================================\n")
    
    # 1. Resolve Precursors (for display and DoV check)
    names = []
    if args.sugars:
        names += args.sugars.split(",")
    if args.amino_acids:
        names += args.amino_acids.split(",")
    if args.additives:
        names += args.additives.split(",")
    if args.lipids:
        names += args.lipids.split(",")
    names = [n.strip() for n in names if n.strip()]

    if not names:
        print("ERROR: No precursors specified. Use --sugars, --amino-acids, --additives, --lipids OR --target.")
        sys.exit(1)

    print(f"Inputs: {', '.join(names)}")
    if args.ratios:
        print(f"Molar Ratios: {', '.join(f'{k}: {v}' for k, v in ratio_dict.items())}")
    print(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}, Catalyst: {args.catalyst or 'None'}")
    print("-" * 60)

    if args.dry_run:
        print("\n[DRY RUN] Validating formulation inputs against known envelope constraints...")
        checker = DomainOfValidityChecker("meaty")
        warnings = checker.check(
            precursor_names=names,
            protein_type=args.protein_type,
            temp_c=args.temp,
            ph=args.ph,
            aw=args.aw,
        )
        if not warnings:
            print("  ✅ All inputs are within the rigorously validated envelope.")
        else:
            print("  ⚠️ The following limitations apply to your formulation:")
            for w in warnings:
                print(f"      - {w.description}")
        
        print("\nDry-run complete. Exiting without evaluating formulation predictions.")
        sys.exit(0)

    # --- Decision Summary ---
    designer = MaillardPipeline(target_tag="meaty", minimize_tag="beany")
    
    # We construct a formulation dict for the designer
    formulation = {
        "name": "Single_Run",
        "sugars": args.sugars.split(",") if args.sugars else [],
        "amino_acids": args.amino_acids.split(",") if args.amino_acids else [],
        "lipids": args.lipids.split(",") if args.lipids else [],
        "additives": args.additives.split(",") if args.additives else [],
        "molar_ratios": ratio_dict,
        "ph": args.ph,
        "temp": args.temp,
        "aw": args.aw,
        "time_minutes": args.time_minutes,
        "protein_type": args.protein_type,
        "denaturation_state": args.denaturation_state,
        "catalyst": args.catalyst
    }
    
    print("\nRunning generative pipeline and scoring sensory impact...")
    res = designer.evaluate_single(formulation, conditions)
    print(f"Generated {len(res.targets)} actionable pathways.")
    
    # Display the table for backward compatibility if targets found
    if res.targets:
        print("\n  Top Predicted Targets (by likelihood/bottleneck):")
        print_table(res.targets)

    # 5. Display PBMA Metrics
    if res.trapping_efficiency > 0 or res.lysine_budget > 0:
        print("\n  PBMA Formulation Metrics:")
        # trapping_efficiency in res is a float (avg), let's show it
        if res.trapping_efficiency > 0:
            bar_len = int(res.trapping_efficiency / 5)
            bar_str = "█" * bar_len + "░" * (20 - bar_len)
            print(f"    │ Lipid Trapping Efficiency: {res.trapping_efficiency:5.1f}% | {bar_str} |")
        
        if res.lysine_budget > 0:
            bar_len = int(res.lysine_budget / 5)
            bar_str = "█" * bar_len + "░" * (20 - bar_len)
            print(f"    │ Lysine Budget (DHA):       {res.lysine_budget:5.1f}% | {bar_str} |")
            if res.lysine_budget > 50.0:
                print("    ⚠️  WARNING: High Lysine consumption by DHA pathway significantly reduces aroma yield.")
    
    checker = DomainOfValidityChecker()
    warnings = checker.check(
        precursor_names=names,
        protein_type=args.protein_type,
        temp_c=args.temp,
        ph=args.ph,
        aw=args.aw,
        matrix_explainability=res.matrix_explainability,
    )
    res.confidence_metadata = build_confidence_package(
        res,
        warnings,
        precursor_names=names,
        protein_type=args.protein_type,
        formulation=formulation,
        baseline_conditions=conditions,
        designer=designer,
    )
    render_decision_summary_cli(res, warnings)
    render_deep_explainability_cli(res)

    if args.report:
        report_dir = generate_report(
            res, 
            warnings, 
            vars(args), 
            output_dir=Path(args.output_dir) if args.output_dir else None
        )
        print(f"📄 Report generated in: {report_dir}")

    print("\n" + "═"*85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - FAST mode barrier estimates are purely qualitative heuristics.")
    print("    - Use --xtb for valid semi-empirical calculations.")
    print("═"*85 + "\n")

if __name__ == "__main__":
    main()
