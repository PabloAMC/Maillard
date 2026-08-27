#!/usr/bin/env python3
"""
scripts/optimize_formulation.py

Runs Bayesian optimization to find the Pareto-optimal Maillard formulation 
balancing target flavor traits against off-flavours and safety risks.
"""

import sys
import argparse
from pathlib import Path
from typing import Any, Dict, List

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.bayesian_optimizer import (  # noqa: E402
    CATEGORICAL_PARAM_KEYS,
    CONCENTRATION_PARAM_KEYS,
    FormulationOptimizer,
    build_molar_ratios,
    formulation_from_params,
)
from src.presentation import render_decision_summary_cli, render_deep_explainability_cli  # noqa: E402
from src.reporting import generate_report  # noqa: E402
from src.usability_reports import prepare_cli_confidence  # noqa: E402

#: A study whose objective never moved has not optimized anything — the "best"
#: trial is just whichever random draw came first. Optuna reports it as a win
#: regardless, so we detect and say so.
DEGENERATE_OBJECTIVE_TOLERANCE = 1.0e-12

_PARAM_LABELS = {
    "sugar_conc": "sugar (all)",
    "aa_conc_sulfur": "amino acids (sulfur)",
    "aa_conc_branched": "amino acids (branched)",
    "aa_conc_basic": "amino acids (basic)",
    "aa_conc_other": "amino acids (other)",
}


def describe_objective_spread(study) -> Dict[str, Any]:
    """Summarise how much the objective actually moved across the study."""
    values = [
        float(trial.value)
        for trial in study.trials
        if trial.value is not None
    ]
    if not values:
        return {"n": 0, "degenerate": True, "spread": 0.0, "best": None, "worst": None}
    spread = max(values) - min(values)
    return {
        "n": len(values),
        "degenerate": len(values) > 1 and spread <= DEGENERATE_OBJECTIVE_TOLERANCE,
        "spread": spread,
        "best": max(values),
        "worst": min(values),
    }


def print_optimized_parameters(
    params: Dict[str, Any],
    *,
    sugars: List[str],
    aas: List[str],
    lipids: List[str],
) -> None:
    """Print the winning trial, routed by parameter kind.

    The previous version formatted every value with ``"{:.3f} M"`` after a
    substring test on the key, which mislabelled the unit (values are mM) and
    crashed on the categorical string parameters.
    """
    print("\n⭐ Optimized Parameter Set:")

    print("  Concentrations (mM):")
    for key in CONCENTRATION_PARAM_KEYS:
        if key in params:
            label = _PARAM_LABELS.get(key, key)
            print(f"    {label.ljust(24)} : {float(params[key]):.3f} mM")

    resolved = build_molar_ratios(params, sugars, aas, lipids)
    if resolved:
        print("  Resolved per-precursor loading (mM):")
        for name in sorted(resolved):
            print(f"    {str(name).ljust(24)} : {float(resolved[name]):.3f} mM")

    print("  Process conditions:")
    if "ph" in params:
        print(f"    {'pH'.ljust(24)} : {float(params['ph']):.2f}")
    if "temp" in params:
        print(f"    {'temperature'.ljust(24)} : {float(params['temp']):.1f} °C")
    if "aw" in params:
        print(f"    {'water activity'.ljust(24)} : {float(params['aw']):.2f}")
    if "time_minutes" in params:
        print(f"    {'time'.ljust(24)} : {float(params['time_minutes']):.1f} min")

    categorical = {key: params[key] for key in CATEGORICAL_PARAM_KEYS if key in params}
    if categorical or "intervention_dose" in params:
        print("  Interventions:")
        for key, value in categorical.items():
            print(f"    {key.ljust(24)} : {value}")
        if "intervention_dose" in params:
            print(f"    {'intervention_dose'.ljust(24)} : {float(params['intervention_dose']):.3f}")

    known = set(CONCENTRATION_PARAM_KEYS) | set(CATEGORICAL_PARAM_KEYS) | {
        "ph", "temp", "aw", "time_minutes", "intervention_dose",
    }
    leftovers = {key: value for key, value in params.items() if key not in known}
    for key, value in sorted(leftovers.items()):
        print(f"    {key.ljust(24)} : {value}")


def main():
    parser = argparse.ArgumentParser(description="Bayesian Formulation Optimizer")
    parser.add_argument("--sugars", type=str, required=True, help="Comma-separated list of sugars")
    parser.add_argument("--amino-acids", type=str, required=True, help="Comma-separated list of amino acids")
    parser.add_argument("--lipids", type=str, default="", help="Comma-separated list of lipids/off-flavor precursors")
    parser.add_argument("--target-tag", type=str, default="meaty", help="Target flavor profile (default: meaty)")
    parser.add_argument("--minimize-tag", type=str, default="beany", help="Flavor to suppress (default: beany)")
    parser.add_argument("--n-iterations", type=int, default=50, help="Number of Optuna trials (default: 50)")
    parser.add_argument("--risk-aversion", type=float, default=1.0, help="Penalty weight for toxic markers (default: 1.0)")
    parser.add_argument("--pre-process", type=str, choices=["none", "yeast", "protease", "both"], default="none", help="Biological pre-processing step")
    # NOTE (2026-08-27): `wheat_gluten` used to be offered here but is not a
    # member of src.matrix_correction.ProteinType, so selecting it killed every
    # trial with "'wheat_gluten' is not a valid ProteinType". The matrix layer
    # has no wheat-gluten accessibility/retention profile (the calibration
    # registry only carries a placeholder), so the choice is withdrawn rather
    # than aliased to a different protein behind the user's back.
    parser.add_argument("--protein-type", choices=["free", "pea_conc", "pea_iso", "soy_conc", "soy_iso", "myco"], default="free", help="Protein matrix type (default: free). wheat_gluten has no calibrated matrix profile and is not selectable.")
    parser.add_argument("--denaturation-state", type=float, default=None, help="Optional denaturation state 0-1. If omitted, infer from temperature/time/pH.")
    parser.add_argument("--report", action="store_true", help="Generate consolidated JSON/Markdown report for the best trial")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save the report")
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
    if args.pre_process != "none":
        # Honest UX: the optimizer treats pre-processing as a searched
        # categorical (`pre_processing` in the trial params), so a fixed
        # --pre-process value has never been applied. Say so instead of
        # printing it as if it were in force.
        print(
            f"NOTE: --pre-process {args.pre_process} is NOT applied. This optimizer "
            "searches pre-processing itself (see the 'pre_processing' trial parameter "
            "in the result below). Use scripts/run_pipeline.py to pin a pre-processing step."
        )
        print("-" * 54)
    print("Starting optimization... this may take a few minutes depending on complexity.\n")

    optimizer = FormulationOptimizer(
        target_tag=args.target_tag,
        minimize_tag=args.minimize_tag,
        risk_aversion=args.risk_aversion,
        protein_type=args.protein_type,
        denaturation_state=args.denaturation_state
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

    spread = describe_objective_spread(study)
    if spread["degenerate"]:
        print(
            "\n⚠️  DEGENERATE SEARCH: the objective returned the identical value "
            f"({spread['best']:.6f}) for all {spread['n']} completed trials "
            f"(spread {spread['spread']:.3e})."
        )
        print(
            "    Nothing was optimized. The 'best' trial below is simply the first "
            "random draw, not a discovered optimum — do not report it as one."
        )
        print(
            "    Likely causes: the target tag scores zero for this precursor set, "
            "or every sampled point saturates the same bound."
        )
    elif spread["n"] > 1:
        print(
            f"\n   Objective spread across {spread['n']} trials: "
            f"{spread['worst']:.4f} → {spread['best']:.4f} (Δ {spread['spread']:.4f})."
        )

    best_trial = study.best_trial

    print_optimized_parameters(best_trial.params, sugars=sugars, aas=aas, lipids=lipids)

    print("\n🔍 Predicted Outcomes:")
    print(f"  Target Score   : {best_trial.user_attrs.get('target_score', 0):.4f}")
    safe_score = best_trial.user_attrs.get('safety_score', 0)
    print(f"  Safety Risk    : {safe_score:.4f} " + ("⚠️" if safe_score > 0 else "✅"))
    print(f"  Off-Flavours   : {best_trial.user_attrs.get('off_flavour_risk', 0):.4f}")
    flagged = best_trial.user_attrs.get('flagged_toxics', [])
    if flagged:
        print(f"  Flagged Toxins : {', '.join(set(flagged))}")
    
    # --- Decision Summary ---
    print("\n🔍 Re-evaluating best formulation for detailed summary...")
    # Rebuild the formulation the winning trial actually scored. This used to
    # pass `best_trial.params` straight in as `molar_ratios`, which (a) never
    # matched a precursor label so every concentration silently reverted to a
    # flat 1:1:1, (b) leaked the trial pH into L-Phenylalanine via substring
    # matching, and (c) crashed on float('rosemary_extract') from the
    # categorical params. `formulation_from_params` is the same mapper the
    # objective uses, so the re-evaluation reproduces the trial.
    best_formulation = formulation_from_params(
        best_trial.params,
        name="Best_Trial",
        fixed_sugars=sugars,
        fixed_amino_acids=aas,
        fixed_lipids=lipids,
        protein_type=args.protein_type,
        denaturation_state=args.denaturation_state,
    )
    # MaillardPipeline is needed here
    from src.pipeline import MaillardPipeline
    designer = MaillardPipeline(args.target_tag, args.minimize_tag)
    from src.smirks_engine import ReactionConditions
    cond = ReactionConditions(
        pH=best_formulation["ph"],
        temperature_celsius=best_formulation["temp"],
        water_activity=best_formulation["aw"],
        protein_type=args.protein_type
    )
    res = designer.evaluate_single(best_formulation, cond)


    warnings = prepare_cli_confidence(
        res,
        target_tag=args.target_tag,
        precursor_names=sugars + aas + lipids,
        protein_type=args.protein_type,
        temp_c=float(res.matrix_explainability.get("temperature_celsius", best_formulation["temp"])),
        ph=float(res.matrix_explainability.get("pH", best_formulation["ph"])),
        aw=float(best_formulation.get("aw", cond.water_activity)),
        formulation=best_formulation,
        baseline_conditions=cond,
        designer=designer,
    )
    render_decision_summary_cli(res, warnings)
    render_deep_explainability_cli(res)

    report_dir = None
    if args.report:
        # Report the formulation that was actually evaluated, not the raw CLI
        # namespace (which says nothing about the optimized concentrations).
        report_inputs: Dict[str, Any] = {
            "mode": "bayesian_optimization_best_trial",
            "target_tag": args.target_tag,
            "minimize_tag": args.minimize_tag,
            "risk_aversion": args.risk_aversion,
            "n_iterations": args.n_iterations,
            "trials_completed": spread["n"],
            "objective_spread": spread["spread"],
            "degenerate_objective": spread["degenerate"],
            "best_objective": study.best_value,
            "best_trial_number": best_trial.number,
            "sugars": ", ".join(sugars) or "-",
            "amino_acids": ", ".join(aas) or "-",
            "lipids": ", ".join(lipids) or "-",
            "molar_ratios_mM": best_formulation.get("molar_ratios", {}),
            "ph": best_formulation.get("ph"),
            "temp": best_formulation.get("temp"),
            "aw": best_formulation.get("aw"),
            "time_minutes": best_formulation.get("time_minutes"),
            "interventions": best_formulation.get("interventions", []),
            "protein_type": args.protein_type,
            "denaturation_state": args.denaturation_state,
        }
        report_dir = generate_report(
            res,
            warnings,
            report_inputs,
            output_dir=Path(args.output_dir) if args.output_dir else None,
        )
    
    if report_dir is not None:
        print(f"📄 Best Trial Report generated in: {report_dir}")
    
    print("\n(Note: Optimization trace saved in memory. For persistent tracking, configure Optuna DB storage.)")
    print("======================================================")

if __name__ == "__main__":
    main()
