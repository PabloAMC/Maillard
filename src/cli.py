import argparse
import sys
import logging
import json
import re
from pathlib import Path
from typing import List, Optional, Dict, Any

from src.pipeline import MaillardPipeline
from src.formulation import Formulation
from src.smirks_engine import ReactionConditions
from src.bayesian_optimizer import FormulationOptimizer
from src import precursor_resolver

# Locate project root
ROOT = Path(__file__).resolve().parents[1]

# Configure logging for CLI
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

def run_pipeline(args):
    """Handler for 'run' subcommand."""
    from src.presentation import render_decision_summary_cli, render_deep_explainability_cli
    from src.usability_reports import DomainOfValidityChecker, build_confidence_package
    from src.reporting import generate_report

    # Parse ratios
    ratio_dict = {}
    if args.ratios:
        for pair in args.ratios.split(","):
            if ":" in pair:
                k, v = pair.split(":")
                ratio_dict[k.strip().lower()] = float(v.strip())

    conditions = ReactionConditions(
        pH=args.ph,
        temperature_celsius=args.temp,
        water_activity=args.aw,
        protein_type=args.protein_type
    )

    designer = MaillardPipeline(target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    
    if args.target_design and not (args.sugars or args.amino_acids):
        # Inverse Design Mode
        logger.info(f"Targeting Sensory Profile: '{args.target_tag}'")
        results = designer.evaluate_all(conditions)
        if results:
            # Print summary table of top candidates
            print("\n  Top Recommended Formulations:")
            print("    ┌──────────────────────────────┬──────────────┬──────────────┬─────────────┬─────────────┐")
            print("    │ FORMULATION                  │ TARGET SCORE │ RISK PENALTY │ LYS BUDGET  │ TRAP EFF    │")
            print("    ├──────────────────────────────┼──────────────┼──────────────┼─────────────┼─────────────┤")
            for res in results[:10]:
                t_score = f"{res.target_score:.1f}"
                r_score = f"{res.off_flavour_risk:.1f}"
                lys = f"{res.lysine_budget:.1f}%"
                trap = f"{res.trapping_efficiency:.1f}%"
                print(f"    │ {res.name[:28]:<28} │ {t_score:>12} │ {r_score:>12} │ {lys:>11} │ {trap:>11} │")
            print("    └──────────────────────────────┴──────────────┴──────────────┴─────────────┴─────────────┘")
            
            best = results[0]
            # Find the best formulation object to extract precursors for reporting
            best_form = next((f for f in designer.grid if f.name == best.name), None)
            precursor_names = []
            if best_form:
                precursor_names = best_form.sugars + best_form.amino_acids + best_form.lipids + best_form.additives
            
            # Use presentation helpers
            render_decision_summary_cli(best, [])
            render_deep_explainability_cli(best)
            
            if args.report:
                report_dir = generate_report(best, [], vars(args), output_dir=Path(args.output_dir) if args.output_dir else None)
                logger.info(f"Report generated in: {report_dir}")
    else:
        # Forward Pipeline Mode
        form = Formulation(
            name="CLI_Run",
            sugars=args.sugars.split(",") if args.sugars else [],
            amino_acids=args.amino_acids.split(",") if args.amino_acids else [],
            lipids=args.lipids.split(",") if args.lipids else [],
            additives=args.additives.split(",") if args.additives else [],
            molar_ratios=ratio_dict,
            ph=args.ph,
            temperature=args.temp,
            water_activity=args.aw,
            protein_type=args.protein_type,
            catalyst=args.catalyst
        )
        res = designer.evaluate_single(form, conditions)
        
        # Resolve all input names for reporting
        input_names = form.sugars + form.amino_acids + form.lipids + form.additives
        
        render_decision_summary_cli(res, [])
        render_deep_explainability_cli(res)
        
        if args.report:
            report_dir = generate_report(res, [], vars(args), output_dir=Path(args.output_dir) if args.output_dir else None)
            logger.info(f"Report generated in: {report_dir}")

def run_optimize(args):
    """Handler for 'optimize' subcommand."""
    optimizer = FormulationOptimizer(
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
        risk_aversion=args.risk_aversion,
        protein_type=args.protein_type
    )
    
    sugars = [s.strip() for s in args.sugars.split(",")] if args.sugars else []
    aas = [a.strip() for a in args.amino_acids.split(",")] if args.amino_acids else []
    lipids = [lip.strip() for lip in args.lipids.split(",")] if args.lipids else []

    study = optimizer.optimize(
        fixed_sugars=sugars,
        fixed_amino_acids=aas,
        fixed_lipids=lipids,
        n_trials=args.n_iterations
    )
    logger.info(f"Optimization Complete. Best Value: {study.best_value:.4f}")

def run_explain(args):
    """Handler for 'explain' subcommand."""
    from src.usability_reports import build_formulation_explainability_payload
    from src.presentation import render_formulation_explainability_markdown

    designer = MaillardPipeline(target_tag=args.target_tag, minimize_tag=args.minimize_tag)
    formulation = next((f for f in designer.grid if f.name == args.name), None)
    if formulation is None:
        logger.error(f"Unknown formulation name: {args.name}")
        sys.exit(1)

    conditions = ReactionConditions(
        pH=formulation.ph,
        temperature_celsius=formulation.temperature,
        water_activity=formulation.water_activity,
        protein_type=formulation.protein_type,
    )
    result = designer.evaluate_single(formulation, conditions)
    payload = build_formulation_explainability_payload(
        formulation.to_dict(),
        result,
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
    )
    markdown = render_formulation_explainability_markdown(payload)

    output_dir = Path(args.output_dir) if args.output_dir else ROOT / "results" / "validation"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^a-z0-9_]+", "_", formulation.name.lower()).strip("_")
    
    markdown_path = output_dir / f"{stem}_explainability.md"
    json_path = output_dir / f"{stem}_explainability.json"
    
    markdown_path.write_text(markdown, encoding="utf-8")
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(markdown)
    logger.info(f"Wrote {markdown_path}")
    logger.info(f"Wrote {json_path}")

def run_list(args):
    """Handler for 'list' subcommand."""
    if args.type == "precursors":
        print("Available Precursors:")
        for name in precursor_resolver.list_available():
            print(f"  - {name}")
    elif args.type == "tags":
        from src.sensory import SensoryDatabase
        db = SensoryDatabase()
        print("Available Sensory Tags (Categories):")
        for tag in db.tags.keys():
            print(f"  - {tag}")

def main():
    parser = argparse.ArgumentParser(prog="maillard", description="Maillard Framework Unified CLI")
    subparsers = parser.add_subparsers(dest="command", help="Subcommand to run")

    # 'run' command
    run_parser = subparsers.add_parser("run", help="Run forward or inverse pipeline")
    run_parser.add_argument("--sugars", type=str, default="", help="Comma-separated sugars")
    run_parser.add_argument("--amino-acids", type=str, default="", help="Comma-separated amino acids")
    run_parser.add_argument("--additives", type=str, default="", help="Comma-separated additives")
    run_parser.add_argument("--lipids", type=str, default="", help="Comma-separated lipids")
    run_parser.add_argument("--ph", type=float, default=6.0, help="pH (default 6.0)")
    run_parser.add_argument("--temp", type=float, default=150.0, help="Temperature Celsius (default 150)")
    run_parser.add_argument("--aw", type=float, default=0.8, help="Water activity (default 0.8)")
    run_parser.add_argument("--ratios", type=str, default="", help="Molar ratios e.g. glucose:2.0")
    run_parser.add_argument("--target-design", action="store_true", help="Enable inverse design mode")
    run_parser.add_argument("--target-tag", type=str, default="meaty", help="Target flavor tag")
    run_parser.add_argument("--minimize-tag", type=str, default="beany", help="Flavor tag to minimize")
    run_parser.add_argument("--protein-type", default="free", help="Protein matrix type")
    run_parser.add_argument("--catalyst", choices=["heme"], default=None, help="Optional catalyst")
    run_parser.add_argument("--report", action="store_true", help="Generate report")
    run_parser.add_argument("--output-dir", type=str, help="Output directory for reports")
    run_parser.set_defaults(func=run_pipeline)

    # 'optimize' command
    opt_parser = subparsers.add_parser("optimize", help="Run Bayesian optimization")
    opt_parser.add_argument("--sugars", type=str, required=True)
    opt_parser.add_argument("--amino-acids", type=str, required=True)
    opt_parser.add_argument("--lipids", type=str, default="")
    opt_parser.add_argument("--target-tag", default="meaty")
    opt_parser.add_argument("--minimize-tag", default="beany")
    opt_parser.add_argument("--n-iterations", type=int, default=20)
    opt_parser.add_argument("--risk-aversion", type=float, default=1.0)
    opt_parser.add_argument("--protein-type", default="free")
    opt_parser.set_defaults(func=run_optimize)

    # 'explain' command
    exp_parser = subparsers.add_parser("explain", help="Generate explainability report for a grid formulation")
    exp_parser.add_argument("--name", required=True, help="Name of formulation in grid")
    exp_parser.add_argument("--target-tag", default="meaty")
    exp_parser.add_argument("--minimize-tag", default="beany")
    exp_parser.add_argument("--output-dir", help="Output directory")
    exp_parser.set_defaults(func=run_explain)

    # 'list' command
    list_parser = subparsers.add_parser("list", help="List available internal datasets")
    list_parser.add_argument("type", choices=["precursors", "tags"], help="What to list")
    list_parser.set_defaults(func=run_list)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        return

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
