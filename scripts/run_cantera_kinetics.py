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
from pathlib import Path
from src.cantera_export import CanteraExporter
from src.kinetics import KineticsEngine
from src.sensory import SensoryPredictor

def run_simulation(barriers_json: str, precursors: dict, temp_c: float, time_sec: float, 
                   predict_sensory: bool = False, output_prefix: str = "simulation"):
    """
    Ties together the microkinetic workflow.
    """
    # 1. Load barriers
    with open(barriers_json, "r") as f:
        barriers = json.load(f)
    
    # 2. Build the mechanism
    print(f"Building Cantera mechanism from {len(barriers)} available barriers...")
    exporter = CanteraExporter()
    
    # Balanced representative Maillard pathway:
    # Ribose (C5H10O5) + Glycine (C2H5NO2) -> Schiff Base (C7H13NO6) + H2O
    # Schiff Base (C7H13NO6) -> FFT (C5H6OS) + ... (Simplified)
    # For the smoke test, we'll use a single balanced step or dummy balanced species
    SMILES_MAP = {
        "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
        "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), # Still unbalanced, but let's try
    }
    
    for key, barrier in barriers.items():
        if barrier is None or key not in SMILES_MAP:
            continue
        reactants, products = SMILES_MAP[key]
        exporter.add_reaction(reactants, products, barrier)
        
    mech_path = f"{output_prefix}_mech.yaml"
    exporter.export_yaml(mech_path)
    
    # 3. Simulate
    print(f"Running isothermal simulation at {temp_c} C for {time_sec} s...")
    engine = KineticsEngine(temperature_k=temp_c + 273.15)
    
    # Precursors: species names from the mechanism (mapped from SMILES)
    # Mapping SMILES to generic S_0... names used by exporter
    spec_map = {v["smiles"]: k for k, v in exporter.species.items()}
    init_state = {}
    for smiles, conc in precursors.items():
        # Check if SMILES or name matches
        found = False
        for k, v in exporter.species.items():
            if v["smiles"] == smiles or v["name"] == smiles:
                init_state[v["name"]] = conc
                found = True
                break
        if not found:
            print(f"Warning: Precursor {smiles} not found in mechanism.")

    results = engine.simulate_network_cantera(mech_path, init_state, (0, time_sec))
    
    # 4. Process results
    df = pd.DataFrame(results)
    csv_path = f"{output_prefix}_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved concentration profiles to {csv_path}")
    
    # 5. Sensory prediction
    if predict_sensory:
        predictor = SensoryPredictor()
        # Take the final concentrations (index -1)
        final_concs = {name: df[name].iloc[-1] for name in results.keys() if name != "time"}
        # Conver fractions to ppm (hacky for demonstration)
        ppm_concs = {name: conc * 1e6 for name, conc in final_concs.items()}
        
        profile = predictor.predict_profile(ppm_concs)
        print("\n====================================")
        print("PREDICTED SENSORY PROFILE (OAV)")
        print("====================================")
        for note, score in predictor.get_dominant_notes(profile):
            print(f"{note.upper():<10} : {score:>10.2f}")
        print("====================================\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cantera microkinetic simulation for Maillard.")
    parser.add_argument("--input", "-i", type=str, default="results/dft_tier2/refinement_all.json",
                        help="Path to the JSON output from run_tier2_dft.py")
    parser.add_argument("--precursors", "-p", type=str, required=True,
                        help="Comma-separated precursors and molarities (e.g., 'ribose:0.1,glycine:0.1')")
    parser.add_argument("--temp", "-T", type=float, default=150.0, help="Temperature in Celsius (default: 150)")
    parser.add_argument("--time", "-t", type=float, default=600.0, help="Total simulation time in seconds")
    parser.add_argument("--output", "-o", type=str, default="results/sim_maillard", help="Output file prefix")
    parser.add_argument("--predict-sensory", action="store_true", help="Predict sensory aroma profile")
    
    args = parser.parse_args()
    
    # Parse precursors
    precursor_dict = {}
    for p in args.precursors.split(","):
        name, conc = p.split(":")
        precursor_dict[name] = float(conc)
        
    run_simulation(args.input, precursor_dict, args.temp, args.time, 
                   predict_sensory=args.predict_sensory, output_prefix=args.output)
