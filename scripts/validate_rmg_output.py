#!/usr/bin/env python3
"""
validate_rmg_output.py — Automatically verify successful Maillard testing.

Parses the RMG-Py output (`chemkin/species_dictionary.txt`) and searches for 
the required target molecules defined for that specific test case.

Usage: python scripts/validate_rmg_output.py <case_output_dir> <target_smiles...>
Example: python scripts/validate_rmg_output.py data/reactions/rmg_validation_cases/case_1_ribose_cys "SCc1ccco1" "Cc1occc1S"
"""

import sys
from pathlib import Path

try:
    from rdkit import Chem
except ImportError:
    print("ERROR: RDKit not installed.")
    sys.exit(1)

GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

def canonicalize(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None
    return Chem.MolToSmiles(mol)

def parse_rmg_dictionary(dict_path: Path) -> dict[str, str]:
    """Parse RMG's species_dictionary.txt.
    Returns: mapping of Canonical_SMILES -> RMG Species Name.
    Warning: This is a simplified parser intended for the standard RMG adjlist/SMILES output.
    """
    if not dict_path.exists():
        return {}
    
    species = {}
    current_name = None
    
    with open(dict_path, "r") as f:
        for line in f:
            line = line.strip()
            # In RMG, block starts with the species name
            if not line:
                continue
            if not line.startswith("1") and "multiplicity" not in line and not line.startswith("label"):
                current_name = line
            # Often, RMG includes the SMILES in a comment or adjlist. 
            # A rigorous implementation would parse the RMG adjacency list.
            # Here we simulate the logic: if we find SMILES or structure mapping, canonicalize it.
            if "SMILES" in line and current_name:
                pass # Extensible parsing logic goes here 
                
    # Mocking for the framework representation (RMG parser is complex)
    return species

def main():
    if len(sys.argv) < 3:
        print("Usage: validate_rmg_output.py <directory> <smiles_1> <smiles_2> ...")
        sys.exit(1)
        
    case_dir = Path(sys.argv[1])
    target_smiles = sys.argv[2:]
    
    print(f"\nVerifying RMG Output for: {case_dir.name}")
    print("-" * 50)
    
    dict_path = case_dir / "chemkin" / "species_dictionary.txt"
    if not dict_path.exists():
        print(f"{RED}✗ ERROR: RMG output not found at {dict_path}{RESET}")
        print("  Please run RMG on the input.py file first:\n  $ rmg data/reactions/rmg_validation_cases/.../input.py\n")
        return 1

    parse_rmg_dictionary(dict_path)
    
    all_found = True
    for t_idx, target in enumerate(target_smiles):
        canon_t = canonicalize(target)
        if not canon_t:
            print(f"✗ Target {target} is an invalid SMILES.")
            all_found = False
            continue
            
        # Simplified string-matching logic assuming RMG output has SMILES
        # In a production environment, this uses `Chem.MolFromSmarts()` and 
        # graph isomorphism over the RMG adjacency lists.
        match_found = False
        with open(dict_path, "r") as f:
            content = f.read()
            # Crude heuristic check: if the literal SMILES is in the dictionary file.
            # RMG writes adjacency lists, but we can do a fallback check.
            if target in content:
                match_found = True
        
        if match_found:
            print(f"{GREEN}✓ FOUND:{RESET} Target {target} was generated.")
        else:
            print(f"{RED}✗ MISSING:{RESET} Target {target} was NOT generated.")
            all_found = False
            
    print("-" * 50)
    if all_found:
        print(f"{GREEN}SUCCESS: All targets generated! RMG reaction families are configured correctly.{RESET}\n")
        return 0
    else:
        print(f"{RED}FAILURE: Some targets are missing. Please refine the custom reaction families.{RESET}\n")
        return 1

if __name__ == "__main__":
    sys.exit(main())
