#!/usr/bin/env python3
"""
scripts/run_tier2_dft.py

Executes the Tier 2 DFT Refinement workflow for the 8 target bifurcations/rate-limiting steps.
Uses the r2SCAN-3c + wB97M-V protocol defined in `src/dft_refiner.py`.
"""

import argparse
from pathlib import Path
import json

from src.dft_refiner import DFTRefiner  # noqa: E402

TARGET_REACTIONS = {
    "amadori": "3.3a: Amadori rearrangement",
    "enolisation": "3.3b: 2,3-enolisation vs 1,2-enolisation",
    "strecker": "3.3c: Strecker decarboxylation",
    "cys_ribose": "3.3d: Cysteine + ribose -> FFT",
    "retro_aldol": "3.3e: Ribose retro-aldol -> MFT",
    "dha": "3.3f: DHA beta-elimination",
    "trapping": "3.3g: Off-flavour trapping",
    "pyrazine": "3.3h: alpha-aminoketone dimerisation"
}

def load_geometries(reaction_key: str) -> dict:
    """Loads reactant and TS geometries from the standard data/geometries/xtb_inputs/ directory."""
    # Special case for amadori which might have a legacy test file
    if reaction_key == "amadori":
        legacy_path = Path("data/geometries/amadori_test.json")
        if legacy_path.exists():
            with open(legacy_path, "r") as f:
                return json.load(f)

    base_path = Path("data/geometries/xtb_inputs") / reaction_key
    r_path = base_path / "reactant.xyz"
    ts_path = base_path / "xtbpath_ts.xyz"
    
    if not r_path.exists() or not ts_path.exists():
        raise FileNotFoundError(f"Missing geometries for {reaction_key} in {base_path}. Have you run scripts/generate_mapped_geometries.py and run_xtb.sh?")
        
    with open(r_path, "r") as f:
        r_xyz = f.read()
    with open(ts_path, "r") as f:
        ts_xyz = f.read()
        
    # Standard Maillard targets are neutral singlets (charge=0, spin=0)
    return {"reactant_xyz": r_xyz, "ts_xyz": ts_xyz, "charge": 0, "spin": 0}

def main():
    parser = argparse.ArgumentParser(description="Tier 2 DFT Refinement execution script.")
    parser.add_argument("--reaction", choices=list(TARGET_REACTIONS.keys()) + ["all"], required=True, 
                        help="Which parameterised reaction to refine.")
    parser.add_argument("--temp", type=float, default=423.15, help="Temperature in Kelvin (default 150 C)")
    parser.add_argument("--solvent", type=str, default="water", help="Implicit solvent name")
    parser.add_argument("--fast", action="store_true", help="Use minimal basis and HF for ultra-fast testing")
    parser.add_argument("--irc", action="store_true", help="Run Phase 3.4 IRC validation for the Transition State")
    
    args = parser.parse_args()
    
    refiner = DFTRefiner(solvent_name=args.solvent, temp_k=args.temp, 
                         db_path="results/maillard_results.db")
    if args.fast:
        print("[FAST MODE] Using minimal basis and methods.")
        refiner.opt_basis = 'sto-3g'
        refiner.opt_method = 'hf'
        refiner.refinement_method = 'hf'
        refiner.refinement_basis = 'sto-3g'
        
    targets = list(TARGET_REACTIONS.keys()) if args.reaction == "all" else [args.reaction]
    
    # Pre-defined SMILES for targets to ensure DB consistency
    # (In a full production engine, these would come from the SMIRKS generator)
    SMILES_DATA = {
        "amadori": {"reactants": ["OCC1OC(O)C(O)C1O", "NCC(O)=O"], "products": ["OCC1OC(O)C(O)C1N=CC(O)=O", "O"]},
        "strecker": {"reactants": ["OCC1OC(O)C(O)C1N=CC(O)=O"], "products": ["C1=C(SC=C1)CS"]}, # Simplified for demo
    }
    
    results = {}
    
    for t in targets:
        print("\n========================================================")
        print(f"Refining: {TARGET_REACTIONS[t]}")
        print("========================================================")
        
        try:
            data = load_geometries(t)
        except Exception as e:
            print(f"SKIPPED: {str(e)}")
            continue
        
        # Meta for DB
        meta = SMILES_DATA.get(t, {})
        if meta:
            meta["family"] = [t]

        try:
            barrier = refiner.calculate_barrier(
                data["reactant_xyz"],
                data["ts_xyz"],
                charge=data.get("charge", 0),
                run_irc=args.irc,
                reaction_meta=meta
            )
            print(f"SUCCESS: Delta G‡ (qh-corrected) = {barrier:.2f} kcal/mol")
            results[t] = barrier
        except Exception as e:
            print(f"FAILED: {str(e)}")
            results[t] = None
            
    # Output to results json
    out_dir = Path("results/dft_tier2")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"refinement_{args.reaction}.json"
    
    with open(out_file, "w") as f:
        json.dump(results, f, indent=4)
        
    print(f"\nResults written to {out_file}")

if __name__ == "__main__":
    main()
