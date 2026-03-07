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

def run_simulation(barriers_json: str, precursors: dict, temp_c: Optional[float] = None, 
                   time_sec: float = 600.0, temp_ramp_csv: Optional[str] = None,
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
    
    # Representative Maillard pathways (balanced via exporter heuristic)
    SMILES_MAP = {
        "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
        "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), 
    }
    
    for key, barrier in barriers.items():
        if barrier is None or key not in SMILES_MAP:
            continue
        reactants, products = SMILES_MAP[key]
        exporter.add_reaction(reactants, products, barrier)
        
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
        lookup = {"ribose": "OCC1OC(O)C(O)C1O", "glycine": "NCC(O)=O"}
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
    
    # 6. Visualization
    try:
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # Plot concentrations on primary axis
        for name in results.keys():
            if name not in ["time", "temperature"] and not name.endswith("_X"):
                ax1.plot(df["time"], df[name], label=name)
        
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Concentration (kmol/m3)")
        ax1.legend(loc="upper left")
        ax1.grid(True, alpha=0.3)
        
        # Plot temperature on secondary axis
        ax2 = ax1.twinx()
        ax2.plot(df["time"], df["temperature"], color="red", linestyle="--", alpha=0.5, label="Temp")
        ax2.set_ylabel("Temperature (K)", color="red")
        ax2.tick_params(axis='y', labelcolor="red")
        
        plt.title(f"Maillard Microkinetics: {'Ramp' if temp_ramp_csv else 'Isothermal'}")
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
        # Convert fractions to ppm
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
    parser.add_argument("--input", "-i", type=str, default="results/dft_tier2/refinement_all.json",
                        help="Path to the JSON output from run_tier2_dft.py")
    parser.add_argument("--precursors", "-p", type=str, required=True,
                        help="Comma-separated precursors and molarities (e.g., 'ribose:0.1,glycine:0.1')")
    parser.add_argument("--temp", "-T", type=float, help="Temperature in Celsius (isothermal)")
    parser.add_argument("--temp-ramp", type=str, help="Path to CSV file with 'time,temp' columns for a ramp.")
    parser.add_argument("--time", "-t", type=float, default=600.0, help="Total simulation time in seconds")
    parser.add_argument("--output", "-o", type=str, default="results/sim_maillard", help="Output file prefix")
    parser.add_argument("--predict-sensory", action="store_true", help="Predict sensory aroma profile")
    
    args = parser.parse_args()
    
    if args.temp is None and args.temp_ramp is None:
        args.temp = 150.0  # Default to isothermal 150 if nothing provided
    
    # Parse precursors
    precursor_dict = {}
    for p in args.precursors.split(","):
        name, conc = p.split(":")
        precursor_dict[name] = float(conc)
        
    run_simulation(args.input, precursor_dict, args.temp, args.time, 
                   temp_ramp_csv=args.temp_ramp,
                   predict_sensory=args.predict_sensory, output_prefix=args.output)
