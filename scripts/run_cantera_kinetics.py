#!/usr/bin/env python3
"""
scripts/run_cantera_kinetics.py

Phase 12: Cantera Microkinetic Simulation CLI.
Ties together barrier input, mechanism generation, ODE integration, 
and sensory prediction.
"""

import argparse
import json
import pandas as pd
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src import data_paths  # noqa: E402
from src.config import DEFAULTS  # noqa: E402
from src.logger import get_logger  # noqa: E402
from src.exceptions import KineticsError  # noqa: E402

logger = get_logger(__name__)

from typing import Optional, List, Dict, Any
from src.cantera_export import CanteraExporter  # noqa: E402
from src.kinetics import KineticsEngine  # noqa: E402
from src.sensory import SensoryPredictor  # noqa: E402
from src.results_db import ResultsDB  # noqa: E402
from src.smirks_engine import SmirksEngine, Species  # noqa: E402
from src.barrier_constants import get_barrier  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402


NAME_MAP = {
    "O=CC(O)C(O)C(O)CO": "ribose",
    "OCC1OC(O)C(O)C1O": "ribose_ring",
    "NCC(O)=O": "glycine",
    "NC(CS)C(=O)O": "cysteine",
    "CC(C)CC(N)C(=O)O": "leucine",
    "OCC1OC(O)C(O)C1N=CC(O)=O": "amadori",
    "SCc1ccco1": "2-furfurylthiol",
    "C1=C(SC=C1)CS": "2-furfurylthiol_legacy",
    "O=Cc1ccco1": "furfural",
    "CC1=NC=C(N=C1)C": "2-5-dimethylpyrazine",
    "CC=O": "acetaldehyde",
    "Cc1cc(SSC2=C(C)OC=C2)co1": "bis(2-methyl-3-furyl) disulfide",
    "O": "water",
    "O=C=O": "CO2",
    "N": "ammonia",
    "S": "H2S",
    "[HH]": "H2"
}

PRECURSOR_LOOKUP = {
    "ribose": "O=CC(O)C(O)C(O)CO", 
    "glycine": "NCC(O)=O",
    "cysteine": "NC(CS)C(=O)O",
    "leucine": "CC(C)CC(N)C(=O)O"
}

LEGACY_SMILES_MAP = {
    "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
    "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), 
}


def build_mechanism(
    exporter: CanteraExporter,
    conditions: ReactionConditions,
    barriers: dict,
    is_db: bool,
    db: Optional[ResultsDB],
    precursors: dict,
    from_smirks: bool,
    verbose_reactions: bool,
    no_gating: bool
) -> None:
    if from_smirks:
        logger.info("Generating dynamic network from precursors using SmirksEngine...")
        engine_smirks = SmirksEngine(conditions=conditions)
        
        precursor_objs = []
        for name, conc in precursors.items():
            smi = PRECURSOR_LOOKUP.get(name.lower(), name) or name
            precursor_objs.append(Species(label=name, smiles=smi))
            
        steps = engine_smirks.enumerate(precursor_objs, max_generations=4)
        logger.info(f"SmirksEngine discovered {len(steps)} elementary steps.")
        
        count = 0
        for step in steps:
            reactants = [s.smiles for s in step.reactants]
            products = [s.smiles for s in step.products]
            
            if is_db and db:
                barrier_kcal, source, _ = db.get_best_barrier(reactants, products, step.reaction_family or "unknown")
            else:
                barrier_kcal = get_barrier(step.reaction_family or "unknown")
                source = "Heuristic"

            for s in step.reactants + step.products:
                if s.smiles in NAME_MAP:
                    exporter.add_species(s.smiles, name=NAME_MAP[s.smiles])
                else:
                    exporter.add_species(s.smiles, name=s.label)
            
            try:
                barrier_val = barrier_kcal[0] if isinstance(barrier_kcal, tuple) else barrier_kcal
                exporter.add_reaction(
                    reactants, products, barrier_val, 
                    thermo_gating=(not no_gating),
                    reaction_family=step.reaction_family,
                    conditions=conditions
                )
                if verbose_reactions:
                    react_str = " + ".join([NAME_MAP.get(s.smiles, s.label) for s in step.reactants])
                    prod_str = " + ".join([NAME_MAP.get(s.smiles, s.label) for s in step.products])
                    logger.info(f"  [{step.reaction_family}] {react_str} -> {prod_str} | Ea: {barrier_kcal:.2f} kcal/mol ({source})")
                count += 1
            except ValueError as e:
                if verbose_reactions:
                    logger.warning(f"  Skipping unbalanced reaction ({step.reaction_family}): {e}")
        
        logger.info(f"Integrated {count} valid reactions into mechanism.")

    elif is_db:
        logger.info("Discovering mechanism from database...")
    else:
        for key, (reactants, products) in LEGACY_SMILES_MAP.items():
            barrier_val = barriers.get(key)
            if barrier_val is not None:
                exporter.add_reaction(reactants, products, barrier_val)


def _plot_results(df: pd.DataFrame, temp_ramp_csv: Optional[str], track_species: Optional[List[str]], output_prefix: str) -> None:
    try:
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        for name in df.columns:
            if name not in ["time", "temperature"] and not name.endswith("_X"):
                linestyle = "-"
                alpha = 0.7
                if track_species and name in track_species:
                    linestyle = "--"
                    alpha = 1.0
                    ax1.plot(df["time"], df[name], label=f"*{name}*", linestyle=linestyle, linewidth=2.5)
                else:
                    ax1.plot(df["time"], df[name], label=name, alpha=alpha)
        
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Concentration (kmol/m3)")
        ax1.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
        ax1.grid(True, alpha=0.3)
        
        ax2 = ax1.twinx()
        ax2.plot(df["time"], df["temperature"], color="red", linestyle=":", alpha=0.5, label="Temp")
        ax2.set_ylabel("Temperature (K)", color="red")
        ax2.tick_params(axis='y', labelcolor="red")
        
        plt.title(f"Maillard Microkinetics: {'Ramp' if temp_ramp_csv else 'Isothermal'}")
        plt.tight_layout()
        plot_path = f"{output_prefix}_plot.png"
        plt.savefig(plot_path)
        logger.info(f"Plot saved to {plot_path}")
    except ImportError:
        logger.warning("Matplotlib not installed, skipping plot generation.")


def _predict_sensory_profile(df: pd.DataFrame) -> None:
    predictor = SensoryPredictor()
    final_concs = {name: df[name].iloc[-1] for name in df.columns if name not in ["time", "temperature"] and not name.endswith("_X")}
    ppm_concs = {name: conc * 1e6 for name, conc in final_concs.items()}
    
    profile = predictor.predict_profile(ppm_concs)
    print("\n====================================")
    print("PREDICTED SENSORY PROFILE (OAV)")
    print("====================================")
    notes = predictor.get_dominant_notes(profile)
    if not notes:
        logger.info("No significant volatiles detected above threshold.")
    else:
        for note, score in notes:
            print(f"{note.upper():<10} : {score:>10.2f}")
    print("====================================\n")


def run_simulation(barriers_json: str, precursors: dict, temp_c: Optional[float] = None, 
                   time_sec: float = 600.0, temp_ramp_csv: Optional[str] = None,
                   predict_sensory: bool = False, output_prefix: str = "simulation",
                   track_species: Optional[List[str]] = None,
                   verbose_reactions: bool = False,
                   no_gating: bool = False,
                   pH: float = 7.0,
                   solvent: str = "water",
                   **kwargs):
    """
    Ties together the microkinetic workflow.
    """
    is_db = barriers_json.endswith(".db")
    db = None
    barriers = {}
    
    if is_db:
        logger.info(f"Loading barriers from database: {barriers_json}")
        db = ResultsDB(db_path=barriers_json)
    else:
        logger.info(f"Loading barriers from JSON: {barriers_json}")
        with open(barriers_json, "r") as f:
            barriers = json.load(f)
    
    conditions = ReactionConditions(pH=pH, solvent_name=solvent, temperature_celsius=temp_c if temp_c else DEFAULTS.default_temp_c)
    exporter = CanteraExporter()
    
    build_mechanism(
        exporter, conditions, barriers, is_db, db, 
        precursors, kwargs.get("from_smirks", False), 
        verbose_reactions, no_gating
    )
            
    logger.info(f"Built Cantera mechanism with {len(exporter.reactions)} reactions.")
    mech_path = f"{output_prefix}_mech.yaml"
    exporter.export_yaml(mech_path)
    
    temp_profile = None
    if temp_ramp_csv:
        logger.info(f"Loading temperature ramp from {temp_ramp_csv}...")
        ramp_df = pd.read_csv(temp_ramp_csv)
        temp_profile = list(zip(ramp_df['time'], ramp_df['temp'] + 273.15))
        max_time = float(ramp_df['time'].max())
        time_sec = max_time
        logger.info(f"Ramp detected: {len(temp_profile)} points over {time_sec} s.")
    else:
        logger.info(f"Running isothermal simulation at {temp_c} C for {time_sec} s...")
    
    engine = KineticsEngine(temperature_k=(temp_c + 273.15) if temp_c else 423.15)
    
    init_state = {}
    for smiles, conc in precursors.items():
        found = False
        norm_smiles = PRECURSOR_LOOKUP.get(smiles.lower(), smiles)
        
        for k, v in exporter.species.items():
            if v["smiles"] == norm_smiles or v["name"] == smiles:
                init_state[v["name"]] = conc
                found = True
                break
        if not found:
            logger.warning(f"Warning: Precursor {smiles} not found in mechanism.")

    try:
        results = engine.simulate_network_cantera(
            mech_path, 
            init_state, 
            (0, time_sec),
            temperature_profile=temp_profile
        )
    except Exception as e:
        raise KineticsError(f"Cantera microkinetic simulation failed: {e}") from e
    
    df = pd.DataFrame(results)
    csv_path = f"{output_prefix}_results.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Saved concentration profiles to {csv_path}")

    if track_species:
        print("\n--- TRACKED SPECIES (Final Yield) ---")
        for s in track_species:
            if s in df.columns:
                val = df[s].iloc[-1]
                print(f"  {s:<20}: {val:.4e} kmol/m3")
            else:
                matches = [c for c in df.columns if s.lower() in c.lower()]
                for m in matches:
                    val = df[m].iloc[-1]
                    print(f"  {m:<20}: {val:.4e} kmol/m3 (partial match)")
        print("-------------------------------------\n")
    
    _plot_results(df, temp_ramp_csv, track_species, output_prefix)

    if predict_sensory:
        _predict_sensory_profile(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cantera microkinetic simulation for Maillard.")
    parser.add_argument("--input", "-i", type=str, default=data_paths.rel(data_paths.RESULTS_DB),
                        help="Path to the Results DB (.db) or legacy JSON.")
    parser.add_argument("--precursors", "-p", type=str, required=True,
                        help="Comma-separated precursors and molarities (e.g., 'ribose:0.1,glycine:0.1')")
    parser.add_argument("--temp", "-T", type=float, help="Temperature in Celsius (isothermal)")
    parser.add_argument("--temp-ramp", type=str, help="Path to CSV file with 'time,temp' columns for a ramp.")
    parser.add_argument("--time", "-t", type=float, default=DEFAULTS.default_cantera_time_sec, help="Total simulation time in seconds")
    parser.add_argument("--output", "-o", type=str, default=data_paths.rel(data_paths.RESULTS_ROOT / "sim_maillard"), help="Output file prefix")
    parser.add_argument("--predict-sensory", action="store_true", help="Predict sensory aroma profile")
    parser.add_argument("--from-smirks", action="store_true", help="Generate network dynamically via SmirksEngine")
    parser.add_argument("--track", type=str, help="Comma-separated species names to track even if not precursors")
    parser.add_argument("--verbose-reactions", action="store_true", help="Print verbose details of each reaction added")
    parser.add_argument("--no-gating", action="store_true", help="Disable thermodynamic gating (skip Delta G check)")
    parser.add_argument("--ph", type=float, default=DEFAULTS.default_ph, help="Reaction pH")
    parser.add_argument("--solvent", type=str, default="water", help="Solvent name (water, lipid, etc.)")
    
    args = parser.parse_args()
    
    if args.temp is None and args.temp_ramp is None:
        args.temp = DEFAULTS.default_temp_c  # Default isothermal if nothing provided
    
    precursor_dict = {}
    for p in args.precursors.split(","):
        name, conc = p.split(":")
        precursor_dict[name] = float(conc)
        
    track_list = args.track.split(",") if args.track else None
    
    run_simulation(args.input, precursor_dict, args.temp, args.time, 
                   temp_ramp_csv=args.temp_ramp,
                   predict_sensory=args.predict_sensory, output_prefix=args.output,
                   from_smirks=args.from_smirks,
                   track_species=track_list,
                   verbose_reactions=args.verbose_reactions,
                   no_gating=args.no_gating,
                   pH=args.ph,
                   solvent=args.solvent)
