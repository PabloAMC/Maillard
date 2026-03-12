"""
scripts/compare_sim_to_lit.py

Automation script to validate microkinetic simulations against experimental literature data.
"""

import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from pathlib import Path

def run_targeted_simulation(temp_c, time_min, precursors):
    """Executes the kinetics CLI with specific parameters."""
    print(f"Running simulation for {precursors} at {temp_c}C for {time_min}min...")
    
    cmd = [
        "python", "scripts/run_cantera_kinetics.py",
        "--precursors", precursors,
        "--temp", str(temp_c),
        "--time", str(time_min * 60),
        "--output", "results/validation/sim_raw"
    ]
    
    import os
    env = os.environ.copy()
    env["PYTHONPATH"] = env.get("PYTHONPATH", "") + ":" + str(Path.cwd())
    
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"Simulation failed:\n{result.stderr}")
        return None
    
    # Load results
    return pd.read_csv("results/validation/sim_raw_results.csv")

def compare(lit_path: str):
    with open(lit_path, "r") as f:
        lit_data = json.load(f)
    
    conditions = lit_data["conditions"]
    exp_yields = lit_data["experimental_yields_ppm"]
    
    # Ensure output directory exists before simulation
    out_dir = Path("results/validation")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Run simulation
    # Defaulting to ribose/glycine 0.1M each as per typical model systems
    df = run_targeted_simulation(
        conditions["temperature_c"], 
        conditions["time_min"], 
        "ribose:0.1,glycine:0.1"
    )
    
    if df is None:
        return

    # Extract final yields (terminal row)
    final_row = df.iloc[-1]
    
    # Map species names in df to lit names (lowercase comparison)
    sim_yields = {}
    for compound in exp_yields.keys():
        # Try to find a column name that contains the compound name and ends with _X
        match = None
        for col in df.columns:
            if compound.lower() in col.lower() and col.endswith("_X"):
                match = col
                break
        
        # Fallback to any match
        if not match:
            for col in df.columns:
                if compound.lower() in col.lower():
                    match = col
                    break
        
        if match:
            # We use mole fraction scaled to ppm
            sim_yields[compound] = float(final_row[match]) * 1e6
        else:
            sim_yields[compound] = 0.0

    # Calculate metrics
    y_exp = np.array(list(exp_yields.values()))
    y_sim = np.array([sim_yields[k] for k in exp_yields.keys()])
    
    # Handle scaling: DFT barriers might be off by a constant offset,
    # so we correlate the relative abundance.
    correlation = np.corrcoef(y_exp, y_sim)[0, 1]
    mae = np.mean(np.abs(y_exp - y_sim))
    
    print("\nComparison Results:")
    print("-" * 30)
    for k in exp_yields.keys():
        print(f"{k:20} | Exp: {exp_yields[k]:8.2f} | Sim: {sim_yields[k]:8.2f}")
    
    print("-" * 30)
    print(f"Pearson R: {correlation:.4f}")
    print(f"MAE:       {mae:.2f} ppm")

    # Plotting
    Path("results/validation").mkdir(parents=True, exist_ok=True)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    x = np.arange(len(exp_yields))
    width = 0.35
    
    # 1. Absolute comparison
    ax1.bar(x - width/2, y_exp, width, label='Exp (Lit)')
    ax1.bar(x + width/2, y_sim, width, label='Sim (Framework)')
    ax1.set_ylabel('Yield / Abundance (ppm)')
    ax1.set_title(f'Absolute Yields (R = {correlation:.3f})')
    ax1.set_xticks(x)
    ax1.set_xticklabels(exp_yields.keys(), rotation=45)
    ax1.legend()
    
    # 2. Relative comparison (Normalized to 1.0)
    y_exp_norm = y_exp / np.max(y_exp) if np.max(y_exp) > 0 else y_exp
    y_sim_norm = y_sim / np.max(y_sim) if np.max(y_sim) > 0 else y_sim
    
    ax2.bar(x - width/2, y_exp_norm, width, label='Exp (Normalized)')
    ax2.bar(x + width/2, y_sim_norm, width, label='Sim (Normalized)')
    ax2.set_ylabel('Relative Abundance')
    ax2.set_title('Normalized Abundance Profiles (Selectivity)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(exp_yields.keys(), rotation=45)
    ax2.legend()
    
    plt.tight_layout()
    plot_path = "results/validation/comparison_plot.png"
    plt.savefig(plot_path)
    print(f"\nPlot saved to {plot_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lit", default="data/lit/ribose_glycine_2021.json")
    args = parser.parse_args()
    
    compare(args.lit)
