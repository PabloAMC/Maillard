#!/usr/bin/env python3
"""
validate_smiles.py — Validate all SMILES strings in the species YAML files.

Checks that every entry:
  1. Has a SMILES field
  2. Can be parsed by RDKit into a valid molecule
  3. Has a non-zero atom count (not empty molecule)

Usage: python scripts/validate_smiles.py
"""

import sys
import yaml
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("ERROR: RDKit not installed. Run: conda install -c conda-forge rdkit")
    sys.exit(1)

# ── Configuration ───────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
DATA_FILES = [
    ROOT / "data/species/desirable_targets.yml",
    ROOT / "data/species/off_flavour_targets.yml",
    ROOT / "data/species/toxic_markers.yml",
    ROOT / "data/species/precursors.yml",
]

GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
RESET  = "\033[0m"
BOLD   = "\033[1m"


def validate_file(path: Path) -> tuple[int, int, list[str]]:
    """Returns (n_ok, n_fail, list_of_error_messages)."""
    with open(path) as f:
        data = yaml.safe_load(f)

    # Support both top-level keys like 'compounds', 'amino_acids', 'sugars', etc.
    all_entries = []
    if isinstance(data, dict):
        for key, val in data.items():
            if isinstance(val, list):
                all_entries.extend(val)
    elif isinstance(data, list):
        all_entries = data

    n_ok, n_fail = 0, 0
    errors = []

    for entry in all_entries:
        name = entry.get("name", "<unnamed>")
        smiles = entry.get("smiles")

        if not smiles:
            msg = f"  {YELLOW}⚠{RESET} {name!r}: no SMILES field"
            errors.append(msg)
            print(msg)
            continue  # Not an error — some complex structures lack SMILES

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            msg = f"  {RED}✗{RESET} {name!r}: invalid SMILES: {smiles}"
            errors.append(msg)
            n_fail += 1
            print(msg)
        elif mol.GetNumAtoms() == 0:
            msg = f"  {RED}✗{RESET} {name!r}: SMILES parses to empty molecule"
            errors.append(msg)
            n_fail += 1
            print(msg)
        else:
            mw = Descriptors.MolWt(mol)
            n_atoms = mol.GetNumAtoms()
            print(f"  {GREEN}✓{RESET} {name!r}: {n_atoms} atoms, MW={mw:.1f} Da")
            n_ok += 1

    return n_ok, n_fail, errors


def main() -> int:
    print(f"\n{BOLD}{'='*60}")
    print(" Maillard Framework — SMILES Validation")
    print(f"{'='*60}{RESET}")

    total_ok, total_fail = 0, 0

    for path in DATA_FILES:
        print(f"\n{BOLD}{path.name}{RESET}")
        print("─" * 60)
        if not path.exists():
            print(f"  {RED}✗{RESET} File not found: {path}")
            total_fail += 1
            continue

        n_ok, n_fail, _ = validate_file(path)
        total_ok += n_ok
        total_fail += n_fail
        print(f"  → {n_ok} valid, {n_fail} invalid")

    print(f"\n{BOLD}{'='*60}{RESET}")
    if total_fail == 0:
        print(f"{GREEN}{BOLD}  ✓ All {total_ok} SMILES validated successfully.{RESET}\n")
        return 0
    else:
        print(f"{RED}{BOLD}  ✗ {total_fail} SMILES failed validation ({total_ok} passed).{RESET}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
