#!/usr/bin/env python3
"""Compare the current recommend implementation against the pre-refactor HEAD version.

This diagnostic focuses on the failing matrix benchmark where Hexanal rises above
its expected rank after the P2 refactor.
"""

from __future__ import annotations

import subprocess
import types
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_evaluator import (
    MATRIX_BENCHMARK_BASE_MARKER_YIELDS,
    MATRIX_BENCHMARK_PROFILES,
)
from src.benchmark_registry import benchmark_to_conditions, benchmark_to_formulation, load_benchmark
from src.barrier_constants import effective_barrier_from_rate_constant
from src.chem_utils import canonicalize_smiles
from src.kinetics import KineticsEngine
from src.lipid_oxidation import predict_hexanal_generation
from src.pipeline import MaillardPipeline
from src.precursor_resolver import resolve_many
from src.projection import DEFAULT_PROJECTION_STRATEGY
from src.smirks_engine import SmirksEngine

import src.recommend as current_recommend
BENCHMARK = ROOT / "data" / "benchmarks" / "pea_isolate_ribose_cysteine_100C_45min_Internal2026.json"


def build_benchmark_inputs() -> tuple[dict, dict, object, dict, dict]:
    bench = load_benchmark(BENCHMARK)
    formulation = benchmark_to_formulation(bench)
    conditions = benchmark_to_conditions(bench)
    designer = MaillardPipeline(target_tag="meaty")
    _kinetics = KineticsEngine(temperature_k=conditions.temperature_kelvin)

    names = (
        formulation["sugars"]
        + formulation["amino_acids"]
        + formulation.get("additives", [])
        + formulation.get("lipids", [])
    )
    precursors = resolve_many(names)
    steps = SmirksEngine(conditions).enumerate(precursors, max_generations=4)

    heuristic_barriers = {}
    for step in steps:
        reactants = [species.smiles for species in step.reactants]
        products = [species.smiles for species in step.products]
        barrier, _source, unc = designer.db.get_best_barrier(
            reactants,
            products,
            step.reaction_family or "unknown",
        )
        rate_constant = conditions.get_rate_constant(
            step.reaction_family or "unknown",
            ea_override_kcal=float(barrier),
        )
        effective_barrier = effective_barrier_from_rate_constant(
            rate_constant,
            conditions.temperature_kelvin,
            step.reaction_family or "unknown",
        )
        rxn_key = (
            f"{'+'.join(sorted(species.smiles for species in step.reactants))}"
            f"->{'+'.join(sorted(species.smiles for species in step.products))}"
        )
        heuristic_barriers[rxn_key] = (max(0.0, effective_barrier), unc)

    initial_concentrations = {}
    ratios = formulation.get("molar_ratios", {})
    for precursor in precursors:
        qty = 1.0
        for key, value in ratios.items():
            if key.lower() in precursor.label.lower() or precursor.label.lower() in key.lower():
                qty = float(value)
                break
        initial_concentrations[
            canonicalize_smiles(precursor.smiles, fallback_to_original=True, strip_salts=True)
        ] = qty

    protein_type = formulation.get("protein_type", "free")
    model = MATRIX_BENCHMARK_PROFILES.get(protein_type)
    if model is not None:
        oxidation = predict_hexanal_generation(
            model["lipid_profile"],
            temp_C=float(conditions.temperature_celsius),
            time_min=float(formulation.get("time_minutes", 60.0)),
            oxygen_availability=1.0,
        )
        oxidation_load_ppb = float(oxidation["total_hydroperoxide"]) * 1000.0
        smiles_map = {
            "Hexanal": "CCCCCC=O",
            "Nonanal": "CCCCCCCCC=O",
            "1-Hexanol": "CCCCCCO",
            "2-Pentylfuran": "CCCCC1=CC=CO1",
        }
        from rdkit import Chem

        for compound_name, yield_factor in MATRIX_BENCHMARK_BASE_MARKER_YIELDS.items():
            smiles = smiles_map.get(compound_name)
            mol = Chem.MolFromSmiles(smiles)
            mw = sum(atom.GetMass() for atom in mol.GetAtoms()) if mol else 100.0
            target_proxy_ppb = oxidation_load_ppb * float(yield_factor)
            molar_conc = target_proxy_ppb / (mw * DEFAULT_PROJECTION_STRATEGY.ppb_conversion_factor)
            canon = canonicalize_smiles(smiles, fallback_to_original=True, strip_salts=True)
            initial_concentrations[canon] = initial_concentrations.get(canon, 0.0) + molar_conc

    return bench, formulation, conditions, steps, heuristic_barriers, initial_concentrations


def load_old_recommend_module() -> types.ModuleType:
    source = subprocess.check_output(["git", "show", "HEAD:src/recommend.py"], cwd=ROOT, text=True)
    module = types.ModuleType("old_recommend")
    module.__dict__["__file__"] = str(ROOT / "src" / "recommend.py")
    sys.modules[module.__name__] = module
    exec(source, module.__dict__)
    return module


def instrument(module: types.ModuleType) -> dict:
    capture: dict = {}

    original_apply = module.apply_matrix_correction
    def wrapped_apply(*args, **kwargs):
        vol, corrected = original_apply(*args, **kwargs)
        capture["corrected_initial"] = dict(corrected)
        return vol, corrected
    module.apply_matrix_correction = wrapped_apply

    original_project = module._project_weighted_flux_to_ppb
    def wrapped_project(*args, **kwargs):
        out = original_project(*args, **kwargs)
        capture["raw_concentrations"] = dict(out)
        return out
    module._project_weighted_flux_to_ppb = wrapped_project

    original_output = module._apply_output_projection
    def wrapped_output(*args, **kwargs):
        obs, meta = original_output(*args, **kwargs)
        capture["observable_volatiles"] = dict(obs)
        capture["projection_metadata"] = dict(meta)
        return obs, meta
    module._apply_output_projection = wrapped_output

    return capture


def top_rows(mapping: dict, limit: int = 5) -> list[tuple[str, float]]:
    return [
        (key, round(value, 6))
        for key, value in sorted(mapping.items(), key=lambda item: item[1], reverse=True)[:limit]
    ]


def main() -> None:
    bench, formulation, conditions, steps, heuristic_barriers, initial_concentrations = build_benchmark_inputs()
    old_recommend = load_old_recommend_module()

    current_capture = instrument(current_recommend)
    old_capture = instrument(old_recommend)

    kwargs = {
        "temperature_kelvin": conditions.temperature_kelvin,
        "time_minutes": float(formulation.get("time_minutes", 60.0)),
        "protein_type": formulation.get("protein_type", "free"),
        "denaturation_state": float(bench.get("denaturation_state", 0.5)),
    }

    current_result = current_recommend.Recommender().predict_from_steps(
        steps, heuristic_barriers, initial_concentrations, **kwargs
    )
    old_result = old_recommend.Recommender().predict_from_steps(
        steps, heuristic_barriers, initial_concentrations, **kwargs
    )

    hexanal_canon = canonicalize_smiles("CCCCCC=O", fallback_to_original=True, strip_salts=True)

    print("CURRENT_TOP5", top_rows(current_result["predicted_ppb"]))
    print("OLD_TOP5", top_rows(old_result["predicted_ppb"]))
    print("CURRENT_HEXANAL", current_result["predicted_ppb"].get("Hexanal"))
    print("OLD_HEXANAL", old_result["predicted_ppb"].get("Hexanal"))
    print("CURRENT_CORRECTED_HEXANAL", current_capture.get("corrected_initial", {}).get(hexanal_canon))
    print("OLD_CORRECTED_HEXANAL", old_capture.get("corrected_initial", {}).get(hexanal_canon))
    print("CURRENT_RAW_TOP5", top_rows(current_capture.get("raw_concentrations", {})))
    print("OLD_RAW_TOP5", top_rows(old_capture.get("raw_concentrations", {})))
    print("CURRENT_OBS_TOP5", top_rows(current_capture.get("observable_volatiles", {})))
    print("OLD_OBS_TOP5", top_rows(old_capture.get("observable_volatiles", {})))
    print("CURRENT_HEX_META", current_capture.get("projection_metadata", {}).get(hexanal_canon))
    print("OLD_HEX_META", old_capture.get("projection_metadata", {}).get(hexanal_canon))


if __name__ == "__main__":
    main()