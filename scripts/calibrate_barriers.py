#!/usr/bin/env python3
"""
scripts/calibrate_barriers.py

Automates the fitting of "Barrier Offsets" per reaction family to minimize the error 
between framework simulations and experimental literature benchmarks (Mottram, Farmer).

Uses Optuna for hyperparameter optimization of activation energies.
"""

import os
import json
import optuna
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any
import subprocess

# Add project root to sys.path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.barrier_constants import FAST_BARRIERS

# Benchmarks to fit against
BENCHMARKS = [
    "data/benchmarks/cys_ribose_150C_Mottram1994.json",
    "data/benchmarks/cys_glucose_150C_Farmer1999.json"
]

def run_sim(offsets: Dict[str, float]) -> float:
    """
    Runs all benchmarks with a given set of offsets and returns total MAE.
    """
    # 1. Create temporary barrier override file
    temp_barriers = FAST_BARRIERS.copy()
    for family, offset in offsets.items():
        if family in temp_barriers:
            val, note = temp_barriers[family]
            temp_barriers[family] = (val + offset, note)
            
    # We'll write this to a temp JSON and have run_cantera_kinetics use it
    # No, run_cantera_kinetics uses src.barrier_constants. 
    # Better: Write a temp_barriers.json and monkeypatch or use a mock.
    # Actually, run_cantera_kinetics.py takes a --input file.
    # If the input is a JSON, it uses it for hardcoded mappings.
    # But for dynamic discovery, it uses get_barrier().
    
    # We'll modify run_cantera_kinetics.py slightly or use an env var to pass offsets.
    
    # For now, let's create a temporary JSON that run_cantera_kinetics can use.
    # Actually, let's just use the subprocess approach with a wrapper.
    
    total_mae = 0.0
    valid_benchmarks = 0
    
    for bench_path in BENCHMARKS:
        try:
            with open(bench_path, "r") as f:
                lit_data = json.load(f)
            
            # Simulation setup
            conditions = lit_data.get("conditions", {})
            temp_c = conditions.get("temp_C", 150.0)
            time_min = conditions.get("time_min", 30.0)
            ph = conditions.get("ph", 7.0)
            
            # Precursors
            precursors = []
            for p, d in lit_data["precursors"].items():
                # Map names to common identifiers
                name = p.replace("L-", "").replace("D-", "").lower()
                conc = d.get("concentration_mM", 100.0)
                precursors.append(f"{name}:{conc/1000.0}")
            prec_str = ",".join(precursors)
            
            env = os.environ.copy()
            env["BARRIER_OFFSETS"] = json.dumps(offsets)
            env["PYTHONPATH"] = str(ROOT)
            
            out_prefix = f"/tmp/calib_{Path(bench_path).stem}"
            cmd = [
                "python", "scripts/run_cantera_kinetics.py",
                "--precursors", prec_str,
                "--temp", str(temp_c),
                "--time", str(time_min * 60),
                "--ph", str(ph),
                "--output", out_prefix,
                "--from-smirks"
            ]
            
            res = subprocess.run(cmd, capture_output=True, text=True, env=env)
            if res.returncode != 0:
                print(f"Sim failed for {bench_path}:\n{res.stderr}")
                continue
            
            df = pd.read_csv(f"{out_prefix}_results.csv")
            print(f"Sim succeeded for {bench_path}. Found {len(df.columns)} species. First 10: {list(df.columns[:10])}...")
            final_row = df.iloc[-1]
            
            # Match using RDKit canonical SMILES
            from rdkit import Chem
            def _canon(s):
                if not s or s.lower() in ["time", "temperature"]: return ""
                m = Chem.MolFromSmiles(s)
                return Chem.MolToSmiles(m, isomericSmiles=False) if m else s

            measured = lit_data.get("measured_volatiles", {})
            # Look for common volatile SMILES in the mechanism
            NAME_TO_SMILES = {
                "furfural": "O=Cc1ccco1",
                "2-methyl-3-furanthiol": "Cc1cc(S)co1",
                "2-furfurylthiol": "SCc1ccco1",
                "bis(2-methyl-3-furyl) disulfide": "Cc1cc(SSC2=C(C)OC=C2)co1",
                "methional": "CSCCC=O",
                "2-5-dimethylpyrazine": "CC1=NC=C(N=C1)C",
                "2-3-dimethylpyrazine": "CC1=NC=CN=C1C",
                "pyrazine": "CC1=NC=C(N=C1)C", 
                "methanethiol": "CS",
                "hydrogen sulfide": "S",
                "dimethyl disulfide": "CSSC"
            }

            benchmark_mae_sum = 0.0
            benchmark_matched_compounds = 0

            for compound, data in measured.items():
                exp_ppb = data.get("conc_ppb", 0.0)
                name_clean = compound.lower().replace("-", " ").replace("_", " ").strip()
                target_smi = _canon(NAME_TO_SMILES.get(compound, ""))
                
                match = None
                # Sort columns to prefer non-_X (concentration) columns
                cols_to_check = sorted(df.columns, key=lambda x: ("_X" in x, len(x)))
                
                # First try: exact or partial name match
                for col in cols_to_check:
                    col_base = col[:-2] if col.endswith("_X") else col
                    col_clean = col_base.lower().replace("-", " ").replace("_", " ").strip()
                    if name_clean == col_clean or name_clean in col_clean or col_clean in name_clean:
                        match = col
                        break
                
                # Second try: SMILES match
                if not match and target_smi:
                    for col in cols_to_check:
                        col_base = col[:-2] if col.endswith("_X") else col
                        # Only try parsing if it looks like a SMILES
                        if any(c in col_base for c in "()=123"):
                            try:
                                col_smi = _canon(col_base)
                                if col_smi == target_smi:
                                    match = col
                                    break
                            except:
                                continue
                
                if match:
                    # moles -> ppb conversion (Assuming 1 kmol/m3 = 1 M)
                    col_base = match[:-2] if match.endswith("_X") else match
                    mol = Chem.MolFromSmiles(col_base)
                    mw = Descriptors.MolWt(mol) if mol else 100.0
                    
                    if match.endswith("_X"):
                        # Mole fraction to ppb: X * 55.5 * MW * 1e6 (55.5 mol/L is water molarity)
                        sim_ppb = float(final_row[match]) * 55.5 * mw * 1e6
                    else:
                        # Molar to ppb: M * MW * 1e6
                        sim_ppb = float(final_row[match]) * mw * 1e6
                        
                    print(f"    FOUND match: {match:25} | Sim: {sim_ppb:10.2f} ppb | Exp: {exp_ppb:10.2f} ppb")
                    total_mae += abs(exp_ppb - sim_ppb)
                    valid_benchmarks += 1
                else:
                    print(f"    NO match for {compound} (SMI: {target_smi})")
                    
        except Exception as e:
            print(f"Error in benchmark {bench_path}: {e}")
            import traceback
            traceback.print_exc()
            
    if valid_benchmarks == 0:
        return 999999.0
    
    # Return MAE in ppb, but print in ppm
    mae_ppb = total_mae / valid_benchmarks
    print(f"  --> Average MAE: {mae_ppb/1000.0:.3f} ppm")
    return mae_ppb

def objective(trial):
    # Optimize major reaction families
    offsets = {
        "schiff_condensation": trial.suggest_float("schiff", -5.0, 20.0),
        "amadori_rearrangement": trial.suggest_float("amadori", -5.0, 20.0),
        "1,2-enolisation": trial.suggest_float("enol", -5.0, 20.0),
        "strecker_degradation": trial.suggest_float("strecker", -5.0, 20.0),
        "cysteine_thermolysis": trial.suggest_float("cys", -5.0, 20.0)
    }
    return run_sim(offsets)

if __name__ == "__main__":
    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=25)
    
    print("\nCalibration Complete!")
    print(f"Best MAE: {study.best_value:.2f} ppb")
    print("Best Offsets:")
    for k, v in study.best_params.items():
        print(f"  {k}: {v:+.2f} kcal/mol")
    
    # Save results
    with open("data/lit/calibration_offsets.json", "w") as f:
        json.dump(study.best_params, f, indent=2)
