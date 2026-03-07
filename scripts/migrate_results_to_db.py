"""
scripts/migrate_results_to_db.py

Migrates legacy JSON results to the new Structured Results Database (SQLite).
"""

import json
import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from src.results_db import ResultsDB

# Legacy SMILES mapping from run_cantera_kinetics.py
SMILES_MAP = {
    "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
    "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), 
}

def migrate(json_path: str, db_path: str):
    print(f"Migrating {json_path} to {db_path}...")
    
    if not Path(json_path).exists():
        print(f"Error: {json_path} not found.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    db = ResultsDB(db_path=db_path)
    
    count = 0
    for key, barrier in data.items():
        if key in SMILES_MAP:
            reactants, products = SMILES_MAP[key]
            db.add_barrier(
                reactants=reactants,
                products=products,
                delta_g_kcal=barrier,
                method="wB97M-V", # Assumed legacy method
                basis="def2-tzvp",
                solvation="water",
                family=key,
                converged=True
            )
            print(f"  Migrated {key}: {barrier:.2f} kcal/mol")
            count += 1
        else:
            print(f"  Warning: No SMILES mapping for {key}, skipping.")

    print(f"Done. Migrated {count} barriers.")

if __name__ == "__main__":
    migrate("results/dft_tier2/refinement_all.json", "results/maillard_results.db")
