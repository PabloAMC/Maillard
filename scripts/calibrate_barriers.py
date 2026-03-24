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
from typing import Dict, List, Any, Tuple
import subprocess
import concurrent.futures

# Add project root to sys.path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.logger import get_logger
from src.exceptions import BenchmarkError

logger = get_logger(__name__)

from src.barrier_constants import FAST_BARRIERS

# Benchmarks to fit against
BENCHMARKS = [
    "data/benchmarks/cys_ribose_150C_Mottram1994.json",
    "data/benchmarks/cys_glucose_150C_Farmer1999.json"
]

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

def compute_mae_for_benchmark(bench_path: str, out_prefix: str, lit_data: dict) -> Tuple[float, int]:
    """
    Computes the Mean Absolute Error between simulation output and benchmark data.
    Returns (total_mae, matched_compounds).
    """
    df = pd.read_csv(f"{out_prefix}_results.csv")
    logger.info(f"Sim succeeded for {bench_path}. Found {len(df.columns)} species. First 10: {list(df.columns[:10])}...")
    final_row = df.iloc[-1]
    
    # Match using RDKit canonical SMILES
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        logger.warning("RDKit is not installed. Continuing without accurate molecular weights.")
        Chem = None
        Descriptors = None
    
    def _canon(s: str) -> str:
        if not s or s.lower() in ["time", "temperature"]: return ""
        if Chem:
            m = Chem.MolFromSmiles(s)
            return Chem.MolToSmiles(m, isomericSmiles=False) if m else s
        return s

    measured = lit_data.get("measured_volatiles", {})
    
    total_mae = 0.0
    matched_compounds = 0

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
                if any(c in col_base for c in "()=123"):
                    try:
                        col_smi = _canon(col_base)
                        if col_smi == target_smi:
                            match = col
                            break
                    except:
                        continue
        
        if match:
            col_base = match[:-2] if match.endswith("_X") else match
            if Chem:
                mol = Chem.MolFromSmiles(col_base)
                mw = Descriptors.MolWt(mol) if mol else 100.0
            else:
                mw = 100.0
            
            if match.endswith("_X"):
                sim_ppb = float(final_row[match]) * 55.5 * mw * 1e6
            else:
                sim_ppb = float(final_row[match]) * mw * 1e6
                
            logger.info(f"    FOUND match: {match:25} | Sim: {sim_ppb:10.2f} ppb | Exp: {exp_ppb:10.2f} ppb")
            total_mae += abs(exp_ppb - sim_ppb)
            matched_compounds += 1
        else:
            logger.debug(f"    NO match for {compound} (SMI: {target_smi})")
            
    return total_mae, matched_compounds


def __run_single_benchmark(bench_path: str, offsets: Dict[str, float]) -> Tuple[float, int]:
    try:
        with open(bench_path, "r") as f:
            lit_data = json.load(f)
    except Exception as e:
        raise BenchmarkError(f"Failed to load benchmark {bench_path}: {e}") from e
    
    conditions = lit_data.get("conditions", {})
    temp_c = conditions.get("temp_C", 150.0)
    time_min = conditions.get("time_min", 30.0)
    ph = conditions.get("ph", 7.0)
    
    precursors = []
    for p, d in lit_data["precursors"].items():
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
        logger.error(f"Sim failed for {bench_path}:\n{res.stderr}")
        return 0.0, 0
    
    return compute_mae_for_benchmark(bench_path, out_prefix, lit_data)


def run_sim(offsets: Dict[str, float]) -> float:
    """
    Runs all benchmarks with a given set of offsets and returns total MAE.
    """
    total_mae = 0.0
    valid_benchmarks = 0
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(__run_single_benchmark, bp, offsets): bp for bp in BENCHMARKS}
        for future in concurrent.futures.as_completed(futures):
            bp = futures[future]
            try:
                b_mae, b_matches = future.result()
                total_mae += b_mae
                valid_benchmarks += b_matches
            except BenchmarkError as e:
                logger.error(f"Benchmark Configuration Error for {bp}: {e}")
            except Exception as e:
                logger.error(f"Unexpected error in benchmark {bp}: {e}")
                import traceback
                logger.error(traceback.format_exc())
            
    if valid_benchmarks == 0:
        return 999999.0
    
    mae_ppb = total_mae / valid_benchmarks
    logger.info(f"  --> Average MAE: {mae_ppb/1000.0:.3f} ppm")
    return mae_ppb


def objective(trial: optuna.Trial) -> float:
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
    
    logger.info("Calibration Complete!")
    logger.info(f"Best MAE: {study.best_value:.2f} ppb")
    logger.info("Best Offsets:")
    for k, v in study.best_params.items():
        logger.info(f"  {k}: {v:+.2f} kcal/mol")
    
    with open("data/lit/calibration_offsets.json", "w") as f:
        json.dump(study.best_params, f, indent=2)
