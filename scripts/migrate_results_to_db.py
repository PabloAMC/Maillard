"""
scripts/migrate_results_to_db.py

Populates the ResultsDB with:
1. Literature heuristic barriers mapped to exact SMILES from curated_pathways.py
2. Legacy JSON DFT results (if available)
"""

import json
import sqlite3
import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from src.results_db import ResultsDB
from data.reactions.curated_pathways import PATHWAYS
from src.barrier_constants import get_barrier, DEFAULT_BARRIER

def migrate(legacy_json_path: str, db_path: str):
    db_path_obj = Path(db_path)
    db_path_obj.parent.mkdir(parents=True, exist_ok=True)
    db = ResultsDB(db_path=db_path)

    print(f"Populating database at {db_path}...")
    
    # 1. Populate Heuristic Barriers from Curated Pathways
    print("\n--- Populating Heuristic Pathways ---")
    count_heuristic = 0
    seen_reactions = set()
    
    for pathway_name, steps in PATHWAYS.items():
        for step in steps:
            # Sort SMILES to ensure uniqueness
            r_smiles = tuple(sorted([s.smiles for s in step.reactants]))
            p_smiles = tuple(sorted([s.smiles for s in step.products]))
            reaction_key = (r_smiles, p_smiles)
            
            if reaction_key in seen_reactions:
                continue
            seen_reactions.add(reaction_key)
            
            barrier_kcal = get_barrier(step.reaction_family)
            if barrier_kcal == DEFAULT_BARRIER:
                print(f"  Warning: Default barrier used for {step.reaction_family}")
            
            db.add_barrier(
                reactants=list(r_smiles),
                products=list(p_smiles),
                delta_g_kcal=barrier_kcal,
                method="literature_heuristic",
                basis=None,
                solvation="water",
                family=step.reaction_family,
                converged=True
            )
            count_heuristic += 1
            print(f"  Inserted heuristic: {step.reaction_family} ({barrier_kcal:.1f} kcal/mol)")
            
    # 2. Migrate legacy DFT JSON if available
    count_dft = 0
    if Path(legacy_json_path).exists():
        print(f"\n--- Migrating Legacy JSON ({legacy_json_path}) ---")
        SMILES_MAP = {
            "amadori": (["OCC1OC(O)C(O)C1O", "NCC(O)=O"], ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]),
            "strecker": (["OCC1OC(O)C(O)C1N=CC(O)=O"], ["C1=C(SC=C1)CS"]), 
        }
        
        with open(legacy_json_path, "r") as f:
            data = json.load(f)
            
        for key, barrier in data.items():
            if barrier is None:
                continue
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
                print(f"  Migrated DFT: {key} ({barrier:.2f} kcal/mol)")
                count_dft += 1
    
    print(f"\nDone. Populated {count_heuristic} heuristics and {count_dft} DFT barriers.")

if __name__ == "__main__":
    migrate("results/dft_tier2/refinement_all.json", "results/maillard_results.db")
