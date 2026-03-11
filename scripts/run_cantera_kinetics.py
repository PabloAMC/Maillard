#!/usr/bin/env python3
"""
scripts/run_cantera_kinetics.py

Phase 12: Cantera Microkinetic Simulation CLI.
Ties together barrier input, mechanism generation, ODE integration, 
and sensory prediction.
"""

import argparse
import json
import sqlite3
import pandas as pd
from pathlib import Path
from typing import Optional, List
from src.cantera_export import CanteraExporter
from src.kinetics import KineticsEngine
from src.sensory import SensoryPredictor
from src.results_db import ResultsDB
from src.smirks_engine import SmirksEngine, Species
from src.barrier_constants import get_barrier
from src.conditions import ReactionConditions

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
    # 1. Load barriers (JSON or SQLite)
    is_db = barriers_json.endswith(".db")
    if is_db:
        print(f"Loading barriers from database: {barriers_json}")
        db = ResultsDB(db_path=barriers_json)
        barriers = {}
    else:
        print(f"Loading barriers from JSON: {barriers_json}")
        with open(barriers_json, "r") as f:
            barriers = json.load(f)
    
    # 2. Build the mechanism
    conditions = ReactionConditions(pH=pH, solvent_name=solvent, temperature_celsius=temp_c if temp_c else 150.0)
    exporter = CanteraExporter()
    
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
        "CC1=NC=C(N=C1)C": "2,5-dimethylpyrazine",
        "CC=O": "acetaldehyde",
        "O": "water",
        "O=C=O": "CO2",
        "N": "ammonia",
        "S": "H2S",
        "[HH]": "H2"
    }
    
    # Precursor normalization lookup
    lookup = {
        "ribose": "O=CC(O)C(O)C(O)CO", 
        "glycine": "NCC(O)=O",
        "cysteine": "NC(CS)C(=O)O",
        "leucine": "CC(C)CC(N)C(=O)O"
    }

    if kwargs.get("from_smirks", False):
        print(f"Generating dynamic network from precursors using SmirksEngine...")
        engine_smirks = SmirksEngine(conditions=conditions)
        
        # Convert precursors to Species objects
        precursor_objs = []
        for name, conc in precursors.items():
            smi = lookup.get(name.lower(), name)
            precursor_objs.append(Species(label=name, smiles=smi))
            
        # Discover network
        steps = engine_smirks.enumerate(precursor_objs, max_generations=3)
        print(f"SmirksEngine discovered {len(steps)} elementary steps.")
        
        count = 0
        for step in steps:
            reactants = [s.smiles for s in step.reactants]
            products = [s.smiles for s in step.products]
            
            # --- DUAL-LOOKUP BARRIER ---
            # Single source of truth: DB first, then heuristic fallback
            if is_db:
                barrier_kcal, source, _ = db.get_best_barrier(reactants, products, step.reaction_family)
            else:
                barrier_kcal = get_barrier(step.reaction_family)
                source = "Heuristic"

            # Add species names if known
            for s in step.reactants + step.products:
                if s.smiles in NAME_MAP:
                    exporter.add_species(s.smiles, name=NAME_MAP[s.smiles])
                else:
                    # Fallback label from SmirksEngine if no common name
                    exporter.add_species(s.smiles, name=s.label)
            
            try:
                exporter.add_reaction(reactants, products, barrier_kcal, 
                                      thermo_gating=(not no_gating),
                                      reaction_family=step.reaction_family,
                                      conditions=conditions)
                if verbose_reactions:
                    react_str = " + ".join([NAME_MAP.get(s.smiles, s.label) for s in step.reactants])
                    prod_str = " + ".join([NAME_MAP.get(s.smiles, s.label) for s in step.products])
                    print(f"  [{step.reaction_family}] {react_str} -> {prod_str} | Ea: {barrier_kcal:.2f} kcal/mol ({source})")
                count += 1
            except ValueError as e:
                if verbose_reactions:
                    print(f"  Skipping unbalanced reaction ({step.reaction_family}): {e}")
        
        print(f"Integrated {count} valid reactions into mechanism.")

    elif is_db:
        print("Discovering mechanism from database...")
    else:
        # Legacy hardcoded fallback for JSON
        SMILES_MAP = {
            "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
            "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), 
        }
        for key, (reactants, products) in SMILES_MAP.items():
            barrier_val = barriers.get(key)
            if barrier_val is not None:
                exporter.add_reaction(reactants, products, barrier_val)
            
    print(f"Built Cantera mechanism with {len(exporter.reactions)} reactions.")
    mech_path = f"{output_prefix}_mech.yaml"
    exporter.export_yaml(mech_path)
    
    # 3. Handle Temperature Profile
    temp_profile = None
    if temp_ramp_csv:
        print(f"Loading temperature ramp from {temp_ramp_csv}...")
        ramp_df = pd.read_csv(temp_ramp_csv)
        # Expect columns 'time' (sec) and 'temp' (Celsius)
        temp_profile = list(zip(ramp_df['time'], ramp_df['temp'] + 273.15))
        max_time = ramp_df['time'].max()
        time_sec = max_time
        print(f"Ramp detected: {len(temp_profile)} points over {time_sec} s.")
    else:
        print(f"Running isothermal simulation at {temp_c} C for {time_sec} s...")
    
    # 4. Simulate
    engine = KineticsEngine(temperature_k=(temp_c + 273.15) if temp_c else 423.15)
    
    # Precursors: species names from the mechanism (mapped from SMILES)
    init_state = {}
    for smiles, conc in precursors.items():
        found = False
        # Normalize common names to SMILES if needed
        lookup = {
            "ribose": "O=CC(O)C(O)C(O)CO", 
            "glycine": "NCC(O)=O",
            "cysteine": "NC(CS)C(=O)O",
            "leucine": "CC(C)CC(N)C(=O)O"
        }
        norm_smiles = lookup.get(smiles.lower(), smiles)
        
        for k, v in exporter.species.items():
            if v["smiles"] == norm_smiles or v["name"] == smiles:
                init_state[v["name"]] = conc
                found = True
                break
        if not found:
            print(f"Warning: Precursor {smiles} not found in mechanism.")

    results = engine.simulate_network_cantera(
        mech_path, 
        init_state, 
        (0, time_sec),
        temperature_profile=temp_profile
    )
    
    # 5. Process results
    df = pd.DataFrame(results)
    csv_path = f"{output_prefix}_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved concentration profiles to {csv_path}")

    # Track specific species
    if track_species:
        print("\n--- TRACKED SPECIES (Final Yield) ---")
        for s in track_species:
            if s in df.columns:
                val = df[s].iloc[-1]
                print(f"  {s:<20}: {val:.4e} kmol/m3")
            else:
                # Try partial match or case-insensitive
                matches = [c for c in df.columns if s.lower() in c.lower()]
                for m in matches:
                    val = df[m].iloc[-1]
                    print(f"  {m:<20}: {val:.4e} kmol/m3 (partial match)")
        print("-------------------------------------\n")
    
    # 6. Visualization
    try:
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # Plot concentrations on primary axis
        # If tracking, only plot tracked species + anything else significant?
        # Actually, let's plot all but highlight tracked ones if specified.
        for name in results.keys():
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
        
        # Plot temperature on secondary axis
        ax2 = ax1.twinx()
        ax2.plot(df["time"], df["temperature"], color="red", linestyle=":", alpha=0.5, label="Temp")
        ax2.set_ylabel("Temperature (K)", color="red")
        ax2.tick_params(axis='y', labelcolor="red")
        
        plt.title(f"Maillard Microkinetics: {'Ramp' if temp_ramp_csv else 'Isothermal'}")
        plt.tight_layout()
        plot_path = f"{output_prefix}_plot.png"
        plt.savefig(plot_path)
        print(f"Plot saved to {plot_path}")
    except ImportError:
        print("Matplotlib not installed, skipping plot generation.")

    # 7. Sensory prediction
    if predict_sensory:
        predictor = SensoryPredictor()
        # Take the final concentrations (index -1)
        final_concs = {name: df[name].iloc[-1] for name in results.keys() if name not in ["time", "temperature"] and not name.endswith("_X")}
        # Convert fractions to ppm (approximation)
        ppm_concs = {name: conc * 1e6 for name, conc in final_concs.items()}
        
        profile = predictor.predict_profile(ppm_concs)
        print("\n====================================")
        print("PREDICTED SENSORY PROFILE (OAV)")
        print("====================================")
        notes = predictor.get_dominant_notes(profile)
        if not notes:
            print("No significant volatiles detected above threshold.")
        else:
            for note, score in notes:
                print(f"{note.upper():<10} : {score:>10.2f}")
        print("====================================\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cantera microkinetic simulation for Maillard.")
    parser.add_argument("--input", "-i", type=str, default="results/maillard_results.db",
                        help="Path to the Results DB (.db) or legacy JSON.")
    parser.add_argument("--precursors", "-p", type=str, required=True,
                        help="Comma-separated precursors and molarities (e.g., 'ribose:0.1,glycine:0.1')")
    parser.add_argument("--temp", "-T", type=float, help="Temperature in Celsius (isothermal)")
    parser.add_argument("--temp-ramp", type=str, help="Path to CSV file with 'time,temp' columns for a ramp.")
    parser.add_argument("--time", "-t", type=float, default=600.0, help="Total simulation time in seconds")
    parser.add_argument("--output", "-o", type=str, default="results/sim_maillard", help="Output file prefix")
    parser.add_argument("--predict-sensory", action="store_true", help="Predict sensory aroma profile")
    parser.add_argument("--from-smirks", action="store_true", help="Generate network dynamically via SmirksEngine")
    parser.add_argument("--track", type=str, help="Comma-separated species names to track even if not precursors")
    parser.add_argument("--verbose-reactions", action="store_true", help="Print verbose details of each reaction added")
    parser.add_argument("--no-gating", action="store_true", help="Disable thermodynamic gating (skip Delta G check)")
    parser.add_argument("--ph", type=float, default=7.0, help="Reaction pH (default 7.0)")
    parser.add_argument("--solvent", type=str, default="water", help="Solvent name (water, lipid, etc.)")
    
    args = parser.parse_args()
    
    if args.temp is None and args.temp_ramp is None:
        args.temp = 150.0  # Default to isothermal 150 if nothing provided
    
    # Parse precursors
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
