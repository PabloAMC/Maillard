#!/usr/bin/env python3
"""
scripts/run_pipeline.py — End-to-End Generative Maillard Pipeline

Legacy compatibility entrypoint preserved for tests and downstream callers.
The implementation keeps the historical CLI surface while delegating to the
current pipeline/reporting stack.
"""

import sys
import argparse
import csv
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.config import DEFAULTS  # noqa: E402
from src.logger import get_logger  # noqa: E402

logger = get_logger(__name__)

from src.conditions import ReactionConditions  # noqa: E402
from src.pipeline import MaillardPipeline  # noqa: E402
from src.recommend import _trunc  # noqa: E402
from src import precursor_resolver  # noqa: E402
from src.barrier_constants import HEME_CATALYST_FAMILIES, HEME_CATALYST_REDUCTION  # noqa: E402
from src.usability_reports import (  # noqa: E402
    DomainOfValidityChecker,
    build_confidence_package,
)
from src.presentation import (  # noqa: E402
    render_decision_summary_cli,
    render_deep_explainability_cli,
)
from src.reporting import generate_report  # noqa: E402


def _load_profile_csv(path: str, *, value_column: str) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "time" not in (reader.fieldnames or []) or value_column not in (reader.fieldnames or []):
            raise ValueError(f"Profile CSV must contain 'time' and '{value_column}' columns")
        for row in reader:
            rows.append((float(row["time"]), float(row[value_column])))
    return rows


def _formulation_to_dict(formulation: Any) -> Dict[str, Any]:
    if hasattr(formulation, "to_dict"):
        return formulation.to_dict()
    if isinstance(formulation, dict):
        return dict(formulation)
    return {
        "name": getattr(formulation, "name", "unknown"),
        "sugars": list(getattr(formulation, "sugars", [])),
        "amino_acids": list(getattr(formulation, "amino_acids", [])),
        "additives": list(getattr(formulation, "additives", [])),
        "lipids": list(getattr(formulation, "lipids", [])),
        "ph": getattr(formulation, "ph", None),
        "temp": getattr(formulation, "temperature", None),
        "aw": getattr(formulation, "water_activity", None),
        "time_minutes": getattr(formulation, "time_minutes", None),
        "protein_type": getattr(formulation, "protein_type", "free"),
        "denaturation_state": getattr(formulation, "denaturation_state", None),
    }


def print_table(active_pathways: List[Dict[str, Any]]) -> None:
    print("    ┌" + "─" * 22 + "┬" + "─" * 16 + "┬" + "─" * 12 + "┬" + "─" * 30 + "┬" + "─" * 20 + "┐")
    print("    │ PREDICTED COMPOUND   │ FORMATION TAG    │ BARRIER    │ SENSORY CHARACTER            │ TOXICITY/RISK        │")
    print("    ├" + "─" * 22 + "┼" + "─" * 16 + "┼" + "─" * 12 + "┼" + "─" * 30 + "┼" + "─" * 20 + "┤")

    for pathway in active_pathways:
        target_str = pathway["target"].label if pathway["target"] else "Unknown"

        tag = ""
        if pathway["type"] == "desirable":
            tag = "[✅ AROMA]"
        elif pathway["type"] == "competing":
            tag = "[⚠️ COMPETING]"
        elif pathway["type"] == "toxic":
            tag = "[☠️ TOXIC]"
        elif pathway["type"] == "masking":
            tag = "[🛡️ MASKING]"

        barrier_str = f"{pathway['span']:.1f} kcal"
        sensory_str = pathway.get("sensory", "-")

        tox_str = "-"
        if pathway.get("toxicity"):
            meta = pathway["toxicity"]
            tox_str = f"[{meta['priority']}] {meta['name']}"

        col1 = _trunc(target_str, 20)
        col2 = _trunc(tag, 14)
        col3 = _trunc(barrier_str, 10)
        col4 = _trunc(sensory_str, 28)
        col5 = _trunc(tox_str, 18)

        print(f"    │ {col1} │ {col2} │ {col3} │ {col4} │ {col5} │")

    print("    └" + "─" * 22 + "┴" + "─" * 16 + "┴" + "─" * 12 + "┴" + "─" * 30 + "┴" + "─" * 20 + "┘")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Maillard formulation screening pipeline.")
    parser.add_argument("--sugars", type=str, default="", help="Comma-separated sugars (e.g. ribose,glucose)")
    parser.add_argument("--amino-acids", type=str, default="", help="Comma-separated amino acids (e.g. cysteine,glycine)")
    parser.add_argument("--additives", type=str, default="", help="Comma-separated additives (e.g. thiamine,glutathione)")
    parser.add_argument("--lipids", type=str, default="", help="Comma-separated lipid aldehydes (e.g. hexanal,nonanal)")
    parser.add_argument("--ph", type=float, default=DEFAULTS.default_ph, help="Reaction pH (default 6.0)")
    parser.add_argument("--temp", type=float, default=DEFAULTS.default_temp_c, help="Heating temperature in Celsius (default 150)")
    parser.add_argument("--ratios", type=str, default="", help="Optional comma-separated molar ratios (e.g. cysteine:2.0,ribose:1.0). Default 1.0.")
    parser.add_argument("--catalyst", type=str, default=None, choices=["heme"], help="Apply catalyst effect (e.g. heme)")
    parser.add_argument("--aw", "--water-activity", type=float, default=DEFAULTS.default_aw, help="Water activity (default 1.0)")
    parser.add_argument("--time-minutes", type=float, default=DEFAULTS.default_time_minutes, help="Reaction time in minutes (default 60.0)")
    parser.add_argument("--target", type=str, default=None, help="Inverse design target sensory tag (e.g. meaty, roasted)")
    parser.add_argument("--minimize", type=str, default=DEFAULTS.default_minimize_tag, help="Inverse design off-flavour tag to minimize (default: beany)")
    parser.add_argument("--xtb", action="store_true", help="Run full GFN2-xTB structural optimizations (SLOW!). Defaults to fast Hammond estimating.")
    parser.add_argument("--protein-type", choices=["free", "pea_conc", "pea_iso", "soy_conc", "soy_iso", "myco"], default=DEFAULTS.default_protein_type, help="Protein matrix type for accessibility corrections.")
    parser.add_argument("--denaturation-state", type=float, default=DEFAULTS.default_denaturation_state, help="Protein denaturation level (0.0 to 1.0). Default 0.5.")
    parser.add_argument("--prediction-mode", choices=["projection", "kinetic"], default="projection", help="Prediction engine to use. Projection stays the default; kinetic runs the ODE microkinetic lane.")
    parser.add_argument("--temp-ramp-csv", type=str, default=None, help="Optional CSV with columns time,temp defining a temperature profile in minutes and Celsius.")
    parser.add_argument("--water-activity-profile-csv", type=str, default=None, help="Optional CSV with columns time,water_activity defining a dynamic water-activity profile.")
    parser.add_argument("--list-precursors", action="store_true", help="List available precursors and exit")
    parser.add_argument("--list-tags", action="store_true", help="List available sensory tags and exit")
    parser.add_argument("--report", action="store_true", help="Generate consolidated JSON/Markdown report")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the report")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs against envelope without running simulation")
    return parser.parse_args()


def output_decision_summary_and_reports(
    res: Any,
    args: argparse.Namespace,
    conditions: ReactionConditions,
    designer: MaillardPipeline,
    formulation: Dict[str, Any],
    precursor_names: List[str],
) -> None:
    checker = DomainOfValidityChecker(args.target if args.target else DEFAULTS.default_target_tag)

    warnings = checker.check(
        precursor_names=precursor_names,
        protein_type=formulation.get("protein_type", args.protein_type),
        temp_c=float(formulation.get("temp", conditions.temperature_celsius)),
        ph=float(formulation.get("ph", conditions.pH)),
        aw=float(formulation.get("aw", conditions.water_activity)),
        matrix_explainability=res.matrix_explainability,
    )
    confidence_payload = build_confidence_package(
        res,
        warnings,
        precursor_names=precursor_names,
        protein_type=formulation.get("protein_type", args.protein_type),
        formulation=formulation,
        baseline_conditions=conditions,
        designer=designer,
    )
    if res.confidence_metadata:
        confidence_payload["prediction_engine"] = res.confidence_metadata.get("prediction_engine", "path_span_projection")
        confidence_payload["kinetics"] = deepcopy(res.confidence_metadata.get("kinetics", {}))
    res.confidence_metadata = confidence_payload

    render_decision_summary_cli(res, warnings)
    render_deep_explainability_cli(res)

    if args.report:
        report_dir = generate_report(
            res,
            warnings,
            vars(args),
            output_dir=Path(args.output_dir) if args.output_dir else None,
        )
        logger.info(f"📄 Report generated in: {report_dir}")


def run_inverse_design(args: argparse.Namespace, conditions: ReactionConditions) -> None:
    print("======================================================")
    print("      Maillard Inverse Design Mode (Phase 7.5)")
    print("======================================================\n")
    logger.info(f"Targeting Sensory Profile: '{args.target}'")
    logger.info(f"Minimizing Risk Profile:   '{args.minimize}'")
    logger.info(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}")
    print("-" * 60)

    designer = MaillardPipeline(args.target, args.minimize)
    logger.info(f"Evaluating {len(designer.grid)} industrial formulations against tags...")

    if args.dry_run:
        logger.info("\n[DRY RUN] Bypassing expensive grid evaluation.")
        logger.info("Validation: Grid dimensions bounded correctly.")
        logger.info("\nDry-run complete. Exiting.")
        sys.exit(0)

    results = designer.evaluate_all(conditions)

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

    if results:
        best = results[0]
        best_formulation_obj = next((item for item in designer.grid if getattr(item, "name", None) == best.name), None)
        best_formulation = _formulation_to_dict(best_formulation_obj) if best_formulation_obj is not None else {}
        best_precursors: List[str] = []
        for key in ("sugars", "amino_acids", "additives", "lipids"):
            best_precursors.extend(best_formulation.get(key, []))

        output_decision_summary_and_reports(
            best,
            args,
            conditions,
            designer,
            best_formulation,
            best_precursors,
        )

    sys.exit(0)


def run_forward_pipeline(args: argparse.Namespace, conditions: ReactionConditions) -> None:
    ratio_dict = {}
    if args.ratios:
        for pair in args.ratios.split(","):
            if ":" in pair:
                key, value = pair.split(":")
                try:
                    ratio_dict[key.strip().lower()] = float(value.strip())
                except ValueError:
                    logger.error(f"Invalid ratio value '{value}' for '{key}'. Must be a float.")
                    sys.exit(1)

    print("======================================================")
    print("      Maillard Generative Pipeline (Phase 7)")
    print("======================================================\n")

    names: List[str] = []
    if args.sugars:
        names += args.sugars.split(",")
    if args.amino_acids:
        names += args.amino_acids.split(",")
    if args.additives:
        names += args.additives.split(",")
    if args.lipids:
        names += args.lipids.split(",")
    names = [name.strip() for name in names if name.strip()]

    if not names:
        logger.error("No precursors specified. Use --sugars, --amino-acids, --additives, --lipids OR --target.")
        sys.exit(1)

    logger.info(f"Inputs: {', '.join(names)}")
    if args.ratios:
        logger.info(f"Molar Ratios: {', '.join(f'{k}: {v}' for k, v in ratio_dict.items())}")
    logger.info(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}, Catalyst: {args.catalyst or 'None'}")
    print("-" * 60)

    if args.dry_run:
        logger.info("\n[DRY RUN] Validating formulation inputs against known envelope constraints...")
        checker = DomainOfValidityChecker(DEFAULTS.default_target_tag)
        warnings = checker.check(
            precursor_names=names,
            protein_type=args.protein_type,
            temp_c=args.temp,
            ph=args.ph,
            aw=args.aw,
        )
        if not warnings:
            logger.info("  ✅ All inputs are within the rigorously validated envelope.")
        else:
            logger.warning("  ⚠️ The following limitations apply to your formulation:")
            for warning in warnings:
                logger.warning(f"      - {warning.description}")

        logger.info("\nDry-run complete. Exiting without evaluating formulation predictions.")
        sys.exit(0)

    designer = MaillardPipeline(target_tag=DEFAULTS.default_target_tag, minimize_tag=DEFAULTS.default_minimize_tag)

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
        "catalyst": args.catalyst,
    }

    logger.info("\nRunning generative pipeline and scoring sensory impact...")
    res = designer.evaluate_single(formulation, conditions)
    logger.info(f"Generated {len(res.targets)} actionable pathways.")

    if res.targets:
        print("\n  Top Predicted Targets (by likelihood/bottleneck):")
        print_table(res.targets)

    if res.trapping_efficiency > 0 or res.lysine_budget > 0:
        print("\n  PBMA Formulation Metrics:")
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

    output_decision_summary_and_reports(res, args, conditions, designer, formulation, names)

    print("\n" + "═" * 85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - FAST mode barrier estimates are purely qualitative heuristics.")
    print("    - Use --xtb for valid semi-empirical calculations.")
    print("═" * 85 + "\n")


def main() -> None:
    args = parse_args()

    if args.list_precursors:
        print("Available Precursors:")
        for name in precursor_resolver.list_available():
            print(f"  - {name}")
        sys.exit(0)

    if args.list_tags:
        from src.sensory import SensoryDatabase

        db = SensoryDatabase()
        print("Available Sensory Tags (Categories):")
        for tag in db.tags.keys():
            print(f"  - {tag}")
        sys.exit(0)

    conditions = ReactionConditions(
        pH=args.ph,
        temperature_celsius=args.temp,
        water_activity=args.aw,
        protein_type=args.protein_type,
        prediction_mode=args.prediction_mode,
        temperature_profile=(
            _load_profile_csv(args.temp_ramp_csv, value_column="temp") if args.temp_ramp_csv else None
        ),
        water_activity_profile=(
            _load_profile_csv(args.water_activity_profile_csv, value_column="water_activity")
            if args.water_activity_profile_csv
            else None
        ),
    )

    if args.target:
        run_inverse_design(args, conditions)
    else:
        run_forward_pipeline(args, conditions)


if __name__ == "__main__":
    main()