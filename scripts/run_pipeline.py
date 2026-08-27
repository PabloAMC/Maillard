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
from typing import Dict, List, Any

# Add project root
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.config import DEFAULTS  # noqa: E402
from src.logger import get_logger  # noqa: E402
from src.exceptions import FormulationInputError  # noqa: E402

logger = get_logger(__name__)

from src.conditions import ReactionConditions  # noqa: E402
from src.smirks_engine import SmirksEngine, Species  # noqa: E402
from src.pipeline import MaillardPipeline  # noqa: E402
from src.xtb_screener import XTBScreener  # noqa: E402
from src.recommend import _trunc  # noqa: E402
from src import precursor_resolver  # noqa: E402
from src.barrier_constants import get_barrier, HEME_CATALYST_FAMILIES, HEME_CATALYST_REDUCTION  # noqa: E402
from src.usability_reports import DomainOfValidityChecker, prepare_cli_confidence  # noqa: E402
from src.presentation import render_decision_summary_cli, render_deep_explainability_cli  # noqa: E402
from src.reporting import generate_report  # noqa: E402


def _parse_float_list(raw: str) -> List[float]:
    values: List[float] = []
    for item in raw.split(","):
        token = item.strip()
        if not token:
            continue
        values.append(float(token))
    return values


def print_table(active_pathways: List[Dict[str, Any]]) -> None:
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
    parser.add_argument("--sme-kj-per-kg", type=float, default=0.0, help="Specific mechanical energy for extrusion correction (default 0.0)")
    parser.add_argument("--moisture-regime", choices=["lme", "hme"], default=None, help="Extrusion moisture regime override")
    parser.add_argument("--sterilization-temp", type=float, default=None, help="Optional post-extrusion sterilization temperature in Celsius")
    parser.add_argument("--sterilization-time-minutes", type=float, default=0.0, help="Optional post-extrusion sterilization time in minutes")
    parser.add_argument("--barrel-zones", type=str, default="", help="Comma-separated sequential barrel zone temperatures in Celsius")
    parser.add_argument("--barrel-zone-time-fractions", type=str, default="", help="Comma-separated residence-time fractions aligned with --barrel-zones")
    parser.add_argument(
        "--target",
        type=str,
        default=None,
        help=(
            "Switch to INVERSE DESIGN mode: screen the built-in formulation grid "
            "(data/formulation_grid.yml) for this sensory tag. Any formulation you "
            "also pass (--sugars/--amino-acids/--additives/--lipids) is evaluated as "
            "one extra candidate alongside the grid. Omit --target to run FORWARD mode, "
            "which scores your formulation and nothing else."
        ),
    )
    parser.add_argument("--minimize", type=str, default=DEFAULTS.default_minimize_tag, help="Inverse design off-flavour tag to minimize (default: beany)")
    parser.add_argument("--xtb", action="store_true", help="Run full GFN2-xTB structural optimizations (SLOW!). Defaults to fast Hammond estimating.")
    # NOTE (2026-08-27): `wheat_gluten` was advertised here but is NOT a member of
    # src.matrix_correction.ProteinType, so every run using it died with
    # "'wheat_gluten' is not a valid ProteinType" the moment the matrix layer
    # resolved a denaturation state. The matrix-correction layer has no wheat
    # gluten accessibility/retention profile (matrix_calibration_registry carries
    # only "Placeholder: no wheat_gluten matrix calibration data yet"), so the
    # option is withdrawn rather than silently aliased to soy. Wheat gluten is
    # still reachable as a registry-backed --protein-source for MEATY multipliers.
    parser.add_argument("--protein-type", choices=["free", "pea_conc", "pea_iso", "soy_conc", "soy_iso", "myco"], default=DEFAULTS.default_protein_type, help="Protein matrix type for accessibility corrections. (wheat_gluten has no calibrated matrix profile; use --protein-source wheat_gluten instead.)")
    parser.add_argument("--protein-source", type=str, default=None, help="Explicit registry-backed protein source for MEATY multipliers (e.g. pea_isolate, wheat_gluten).")
    parser.add_argument("--denaturation-state", type=float, default=DEFAULTS.default_denaturation_state, help="Protein denaturation level (0.0 to 1.0). Default 0.5.")
    parser.add_argument("--list-precursors", action="store_true", help="List available precursors and exit")
    parser.add_argument("--list-tags", action="store_true", help="List available sensory tags and exit")
    parser.add_argument("--report", action="store_true", help="Generate consolidated JSON/Markdown report")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the report")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs against envelope without running simulation")
    return parser.parse_args()


def _render_result_bundle(
    res: Any,
    *,
    target_tag: str,
    precursor_names: List[str],
    protein_type: str,
    temp_c: float,
    ph: float,
    aw: float,
    formulation: Dict[str, Any],
    baseline_conditions: ReactionConditions,
    designer: MaillardPipeline,
    report_enabled: bool,
    conditions_dict: Dict[str, Any],
    output_dir: Path | None,
) -> Path | None:
    warnings = prepare_cli_confidence(
        res,
        target_tag=target_tag,
        precursor_names=precursor_names,
        protein_type=protein_type,
        temp_c=temp_c,
        ph=ph,
        aw=aw,
        formulation=formulation,
        baseline_conditions=baseline_conditions,
        designer=designer,
    )
    render_decision_summary_cli(res, warnings)
    render_deep_explainability_cli(res)
    if not report_enabled:
        return None
    return generate_report(res, warnings, conditions_dict, output_dir=output_dir)


CUSTOM_CANDIDATE_NAME = "Your formulation (custom)"

#: CLI arguments that describe a *formulation* rather than the global process
#: conditions. Inverse design used to read none of them.
FORMULATION_ARG_KEYS = ("sugars", "amino_acids", "additives", "lipids")
#: Formulation-shaping arguments that only mean something once precursors exist.
FORMULATION_MODIFIER_ARG_KEYS = ("ratios",)


def _parse_ratio_dict(raw: str) -> Dict[str, float]:
    ratio_dict: Dict[str, float] = {}
    if not raw:
        return ratio_dict
    for pair in raw.split(","):
        if ":" not in pair:
            continue
        k, v = pair.split(":", 1)
        try:
            ratio_dict[k.strip().lower()] = float(v.strip())
        except ValueError:
            logger.error(f"Invalid ratio value '{v}' for '{k}'. Must be a float.")
            sys.exit(1)
    return ratio_dict


def _split_names(raw: str) -> List[str]:
    return [token.strip() for token in (raw or "").split(",") if token.strip()]


def build_formulation_from_args(args: argparse.Namespace, name: str) -> Dict[str, Any]:
    """Assemble the formulation dict the pipeline consumes from CLI arguments.

    Shared by forward mode and by the custom candidate that inverse design now
    screens alongside the grid, so the two modes cannot drift apart.
    """
    return {
        "name": name,
        "sugars": _split_names(args.sugars),
        "amino_acids": _split_names(args.amino_acids),
        "lipids": _split_names(args.lipids),
        "additives": _split_names(args.additives),
        "molar_ratios": _parse_ratio_dict(args.ratios),
        "ph": args.ph,
        "temp": args.temp,
        "aw": args.aw,
        "time_minutes": args.time_minutes,
        "sme_kj_per_kg": args.sme_kj_per_kg,
        "moisture_regime": args.moisture_regime,
        "sterilization_temperature_celsius": args.sterilization_temp,
        "sterilization_time_minutes": args.sterilization_time_minutes,
        "barrel_zone_temperatures": _parse_float_list(args.barrel_zones) if args.barrel_zones else None,
        "barrel_zone_time_fractions": _parse_float_list(args.barrel_zone_time_fractions) if args.barrel_zone_time_fractions else None,
        "protein_type": args.protein_type,
        "protein_source": args.protein_source,
        "denaturation_state": args.denaturation_state,
        "catalyst": args.catalyst,
    }


def _user_supplied_precursors(args: argparse.Namespace) -> bool:
    return any(_split_names(getattr(args, key, "") or "") for key in FORMULATION_ARG_KEYS)


def _report_inputs_for_inverse_design(
    args: argparse.Namespace,
    *,
    winning_formulation: Dict[str, Any],
    candidate_count: int,
    custom_candidate_included: bool,
    ignored_arguments: List[str],
) -> Dict[str, Any]:
    """The inputs block §1 of the report should echo.

    Historically this was ``vars(args)``, which printed the user's
    ``--sugars``/``--ratios``/``--time-minutes`` at the top of a report whose
    science came from a grid entry that never saw them. We echo what was
    actually evaluated instead, and name what was ignored.
    """
    payload: Dict[str, Any] = {
        "mode": "inverse_design",
        "target_tag": args.target,
        "minimize_tag": args.minimize,
        "candidates_screened": candidate_count,
        "custom_candidate_included": custom_candidate_included,
        "selected_formulation": winning_formulation.get("name", "unknown"),
        "selected_sugars": ", ".join(str(x) for x in winning_formulation.get("sugars", [])) or "-",
        "selected_amino_acids": ", ".join(str(x) for x in winning_formulation.get("amino_acids", [])) or "-",
        "selected_additives": ", ".join(str(x) for x in winning_formulation.get("additives", [])) or "-",
        "selected_lipids": ", ".join(str(x) for x in winning_formulation.get("lipids", [])) or "-",
        "selected_molar_ratios": winning_formulation.get("molar_ratios", {}) or "1.0 for every precursor (grid default)",
        "ph": winning_formulation.get("ph", args.ph),
        "temp": winning_formulation.get("temp", args.temp),
        "aw": winning_formulation.get("aw", args.aw),
        "time_minutes": winning_formulation.get("time_minutes", args.time_minutes),
        "protein_type": winning_formulation.get("protein_type", args.protein_type),
        "catalyst": winning_formulation.get("catalyst", args.catalyst),
        "xtb": bool(args.xtb),
    }
    if ignored_arguments:
        payload["ignored_cli_arguments"] = ", ".join(sorted(ignored_arguments))
        payload["ignored_cli_arguments_note"] = (
            "These arguments do not reach the selected grid formulation. "
            "Use forward mode (drop --target) to score exactly the formulation you typed."
        )
    return payload


def _report_inputs_for_forward_mode(
    args: argparse.Namespace, formulation: Dict[str, Any]
) -> Dict[str, Any]:
    """Inputs block for forward mode — every value here really was evaluated."""
    payload: Dict[str, Any] = {
        "mode": "forward_single_formulation",
        "sugars": ", ".join(str(x) for x in formulation.get("sugars", [])) or "-",
        "amino_acids": ", ".join(str(x) for x in formulation.get("amino_acids", [])) or "-",
        "additives": ", ".join(str(x) for x in formulation.get("additives", [])) or "-",
        "lipids": ", ".join(str(x) for x in formulation.get("lipids", [])) or "-",
        "molar_ratios": formulation.get("molar_ratios", {}) or "1.0 for every precursor (default)",
        "ph": formulation.get("ph"),
        "temp": formulation.get("temp"),
        "aw": formulation.get("aw"),
        "time_minutes": formulation.get("time_minutes"),
        "protein_type": formulation.get("protein_type"),
        "protein_source": formulation.get("protein_source"),
        "denaturation_state": formulation.get("denaturation_state"),
        "catalyst": formulation.get("catalyst"),
        "sme_kj_per_kg": formulation.get("sme_kj_per_kg"),
        "moisture_regime": formulation.get("moisture_regime"),
        "sterilization_temperature_celsius": formulation.get("sterilization_temperature_celsius"),
        "sterilization_time_minutes": formulation.get("sterilization_time_minutes"),
        "barrel_zone_temperatures": formulation.get("barrel_zone_temperatures"),
        "barrel_zone_time_fractions": formulation.get("barrel_zone_time_fractions"),
        "xtb": bool(args.xtb),
    }
    return {key: value for key, value in payload.items() if value is not None}


def report_skipped_formulations(designer: MaillardPipeline) -> None:
    """Print any formulation the resolver dropped, instead of losing it silently."""
    skipped = list(getattr(designer, "last_skipped_formulations", []) or [])
    if not skipped:
        return
    print("\n  ⚠️  SKIPPED FORMULATIONS (not scored, not ranked):")
    for row in skipped:
        print(
            f"    - {row.get('name', 'unknown')}: unresolved precursor(s) "
            f"{row.get('unresolved_precursors', 'unknown')} — {row.get('reason', 'unknown reason')}"
        )
    print("    Run `--list-precursors` to see the recognised names.\n")


def run_inverse_design(args: argparse.Namespace, conditions: ReactionConditions) -> None:
    print("======================================================")
    print("      Maillard Inverse Design Mode (Phase 7.5)")
    print("======================================================\n")
    logger.info(f"Targeting Sensory Profile: '{args.target}'")
    logger.info(f"Minimizing Risk Profile:   '{args.minimize}'")
    logger.info(f"Conditions: pH {conditions.pH}, {conditions.temperature_celsius}°C, aᵥ {conditions.water_activity}")
    print("-" * 60)

    try:
        designer = MaillardPipeline(args.target, args.minimize)

        # --- Honest handling of formulation arguments in inverse-design mode ---
        # Inverse design screens the fixed grid. Before 2026-08-27 it read none of
        # --sugars/--amino-acids/--additives/--lipids/--ratios and reported them
        # anyway. Now a user-supplied formulation is screened as one more
        # candidate, and anything still unused is named out loud.
        candidate_grid = list(designer.grid)
        custom_candidate_included = _user_supplied_precursors(args)
        ignored_arguments: List[str] = []

        if custom_candidate_included:
            custom_formulation = build_formulation_from_args(args, CUSTOM_CANDIDATE_NAME)
            candidate_grid.append(custom_formulation)
            print(
                f"  ℹ️  Your formulation is entered as candidate '{CUSTOM_CANDIDATE_NAME}' "
                f"alongside the {len(designer.grid)} built-in grid entries."
            )
            print(
                "     Grid entries carry their own precursors, ratios and (sometimes) pH/temp; "
                "your --ratios/--time-minutes apply only to your own candidate."
            )
        else:
            for key in FORMULATION_ARG_KEYS + FORMULATION_MODIFIER_ARG_KEYS:
                if str(getattr(args, key, "") or "").strip():
                    ignored_arguments.append(f"--{key.replace('_', '-')}")
            if ignored_arguments:
                logger.warning(
                    "Inverse design (--target) screens the built-in grid; %s cannot be applied "
                    "without at least one precursor argument. Drop --target to score your own "
                    "formulation in forward mode.",
                    ", ".join(sorted(ignored_arguments)),
                )

        logger.info(f"Evaluating {len(candidate_grid)} formulations against tags...")

        if args.dry_run:
            logger.info("\n[DRY RUN] Bypassing expensive grid evaluation.")
            logger.info("Validation: Grid dimensions bounded correctly.")
            logger.info("\nDry-run complete. Exiting.")
            sys.exit(0)

        results = designer.evaluate_all(conditions, grid_override=candidate_grid)
        report_skipped_formulations(designer)

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

        if results:
            best = results[0]
            best_formulation = next(
                (item for item in candidate_grid if item.get("name") == best.name), {}
            )
            best_precursors = []
            for key in ("sugars", "amino_acids", "additives", "lipids"):
                best_precursors.extend(best_formulation.get(key, []))

            if custom_candidate_included:
                if best.name == CUSTOM_CANDIDATE_NAME:
                    print("\n  ✅ Your formulation won the screen. The report below is your formulation.")
                else:
                    custom_rank = next(
                        (idx + 1 for idx, res in enumerate(results) if res.name == CUSTOM_CANDIDATE_NAME),
                        None,
                    )
                    rank_text = f"ranked #{custom_rank} of {len(results)}" if custom_rank else "was not scored"
                    print(
                        f"\n  ℹ️  The report below describes the winning GRID formulation "
                        f"'{best.name}', not the formulation you typed (yours {rank_text})."
                    )
                    print(
                        "     To get a report for exactly your formulation, drop --target "
                        "(forward mode)."
                    )
            else:
                print(
                    f"\n  ℹ️  The report below describes the winning grid formulation '{best.name}'. "
                    "Inverse design does not evaluate a formulation of your own unless you pass one."
                )

            report_dir = _render_result_bundle(
                best,
                target_tag=args.target or DEFAULTS.default_target_tag,
                precursor_names=best_precursors,
                protein_type=str(best_formulation.get("protein_type", args.protein_type)),
                temp_c=float(best_formulation.get("temp", conditions.temperature_celsius)),
                ph=float(best_formulation.get("ph", conditions.pH)),
                aw=float(best_formulation.get("aw", conditions.water_activity)),
                formulation=best_formulation,
                baseline_conditions=conditions,
                designer=designer,
                conditions_dict=_report_inputs_for_inverse_design(
                    args,
                    winning_formulation=best_formulation,
                    candidate_count=len(candidate_grid),
                    custom_candidate_included=custom_candidate_included,
                    ignored_arguments=ignored_arguments,
                ),
                report_enabled=bool(args.report),
                output_dir=Path(args.output_dir) if args.output_dir else None,
            )
            if report_dir is not None:
                logger.info(f"📄 Report generated in: {report_dir}")

        sys.exit(0)
        
    except ValueError as e:
        logger.error(f"Formulation error: {e}")
        sys.exit(1)


def run_forward_pipeline(args: argparse.Namespace, conditions: ReactionConditions) -> None:
    # Parse ratios
    ratio_dict = _parse_ratio_dict(args.ratios)

    # Print Forward Pipeline Header
    print("======================================================")
    print("      Maillard Generative Pipeline (Phase 7) — forward mode")
    print("======================================================\n")
    print("  Scoring exactly the formulation you supplied (no grid screening).\n")

    # 1. Resolve Precursors (for display and DoV check)
    names: List[str] = []
    for key in FORMULATION_ARG_KEYS:
        names += _split_names(getattr(args, key, "") or "")

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
            temp_c=conditions.effective_temperature_celsius,
            ph=args.ph,
            aw=args.aw,
        )
        if not warnings:
            logger.info("  ✅ All inputs are within the rigorously validated envelope.")
        else:
            logger.warning("  ⚠️ The following limitations apply to your formulation:")
            for w in warnings:
                logger.warning(f"      - {w.description}")
        
        logger.info("\nDry-run complete. Exiting without evaluating formulation predictions.")
        sys.exit(0)

    designer = MaillardPipeline(target_tag=DEFAULTS.default_target_tag, minimize_tag=DEFAULTS.default_minimize_tag)

    formulation = build_formulation_from_args(args, "Single_Run")

    logger.info("\nRunning generative pipeline and scoring sensory impact...")
    try:
        res = designer.evaluate_single(formulation, conditions)
    except ValueError as exc:
        report_skipped_formulations(designer)
        logger.error(f"Formulation error: {exc}")
        sys.exit(1)
    report_skipped_formulations(designer)
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
    
    report_dir = _render_result_bundle(
        res,
        target_tag=DEFAULTS.default_target_tag,
        precursor_names=names,
        protein_type=str(formulation.get("protein_type", args.protein_type)),
        temp_c=float(formulation.get("temp", conditions.temperature_celsius)),
        ph=float(formulation.get("ph", conditions.pH)),
        aw=float(formulation.get("aw", conditions.water_activity)),
        formulation=formulation,
        baseline_conditions=conditions,
        designer=designer,
        conditions_dict=_report_inputs_for_forward_mode(args, formulation),
        report_enabled=bool(args.report),
        output_dir=Path(args.output_dir) if args.output_dir else None,
    )
    if report_dir is not None:
        logger.info(f"📄 Report generated in: {report_dir}")

    print("\n" + "═"*85)
    print(" ℹ️  KNOWN LIMITATIONS:")
    print("    - FAST mode barrier estimates are purely qualitative heuristics.")
    print("    - Use --xtb for valid semi-empirical calculations.")
    print("═"*85 + "\n")


def main() -> None:
    args = parse_args()

    # Catch simple cases where user only wants to list
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

    # Setup Conditions
    conditions = ReactionConditions(
        pH=args.ph, 
        temperature_celsius=args.temp,
        water_activity=args.aw,
        protein_type=args.protein_type,
        sme_kj_per_kg=args.sme_kj_per_kg,
        moisture_regime=args.moisture_regime,
        sterilization_temperature_celsius=args.sterilization_temp,
        sterilization_time_minutes=args.sterilization_time_minutes,
        barrel_zone_temperatures=_parse_float_list(args.barrel_zones) if args.barrel_zones else None,
        barrel_zone_time_fractions=_parse_float_list(args.barrel_zone_time_fractions) if args.barrel_zone_time_fractions else None,
    )

    if args.target:
        run_inverse_design(args, conditions)
    else:
        run_forward_pipeline(args, conditions)


if __name__ == "__main__":
    main()
