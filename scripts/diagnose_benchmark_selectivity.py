#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import load_benchmark, benchmark_to_conditions, benchmark_to_formulation
from src.inverse_design import InverseDesigner
from src.precursor_resolver import resolve_many
from src.recommend import Recommender, _canon
from src.smirks_engine import SmirksEngine


def _build_rec_result(bench_path: str):
    bench = load_benchmark(bench_path)
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)

    names = formulation["sugars"] + formulation["amino_acids"] + formulation.get("additives", []) + formulation.get("lipids", [])
    precursors = resolve_many(names)

    designer = InverseDesigner(target_tag="meaty")
    engine = SmirksEngine(conditions)
    steps = engine.enumerate(precursors, max_generations=4)

    heuristic_barriers = {}
    for step in steps:
        reactants = [s.smiles for s in step.reactants]
        products = [s.smiles for s in step.products]
        bar, source, unc = designer.db.get_best_barrier(reactants, products, step.reaction_family or "unknown")
        k = conditions.get_rate_constant(step.reaction_family or "unknown", ea_override_kcal=bar)
        rt = 0.001987 * conditions.temperature_kelvin
        pre_exponential = 1e11
        bar_eff = -rt * math.log(k / pre_exponential) if k > 0 else 99.0
        rxn_key = f"{'+'.join(sorted(r.smiles for r in step.reactants))}->{'+'.join(sorted(p.smiles for p in step.products))}"
        heuristic_barriers[rxn_key] = (max(0.0, bar_eff), unc)

    initial_concentrations = {}
    ratios = formulation.get("molar_ratios", {})
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[_canon(precursor.smiles)] = qty

    rec = Recommender()
    rec_result = rec.predict_from_steps(
        steps,
        heuristic_barriers,
        initial_concentrations,
        temperature_kelvin=conditions.temperature_kelvin,
        time_minutes=formulation.get("time_minutes"),
        protein_type=formulation.get("protein_type", "free"),
        denaturation_state=formulation.get("denaturation_state", 0.5),
    )
    return bench, rec_result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lit", required=True)
    args = parser.parse_args()

    bench, rec_result = _build_rec_result(args.lit)
    debug_paths = rec_result.get("debug_paths", {})
    predicted_ppb = rec_result.get("predicted_ppb", {})
    species_names = rec_result.get("species_names", {})

    print(f"Benchmark: {bench['benchmark_id']}")
    for compound in bench["measured_volatiles"].keys():
        print(f"\nCompound: {compound}")
        matched_name = None
        matched_ppb = 0.0
        matched_path = None
        for canon, path in debug_paths.items():
            display_name = species_names.get(canon, canon)
            if compound.lower() in display_name.lower() or display_name.lower() in compound.lower():
                matched_name = display_name
                matched_ppb = predicted_ppb.get(display_name, predicted_ppb.get(canon, 0.0))
                matched_path = path
                break

        print(f"Matched: {matched_name or '-'} | Predicted ppb: {matched_ppb:.6g}")
        if not matched_path:
            print("No path trace available")
            continue
        for idx, step in enumerate(matched_path, start=1):
            print(
                f"  {idx}. {step['family']} | barrier={step['barrier']:.2f} | span={step['path_span']:.2f} | "
                f"reactants={', '.join(step['reactants'])} -> products={', '.join(step['products'])}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())