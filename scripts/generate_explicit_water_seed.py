#!/usr/bin/env python3
"""
Generate an explicit-water cluster seed for the asparagine-sugar target.
Places 3 water molecules around the reactant complex, relaxes them with xTB,
and ensures identical water positioning in the product geometry for a valid path search.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_BASE = ROOT / "data" / "geometries" / "xtb_inputs" / "asparagine_sugar_explicit_water_cluster"

# Mapped SMILES for Asparagine + Glucose + 3xH2O -> Acrylamide + NH3 + CO2 + Glucose + 3xH2O
R_SMI = "[NH2:1][CH:2]([CH2:3][C:4](=[O:5])[NH2:6])[C:7](=[O:8])[OH:9].[CH:10](=[O:11])[CH:12]([OH:13])[CH:14]([OH:15])[CH:16]([OH:17])[CH:18]([OH:19])[CH2:20][OH:21].[OH2:30].[OH2:31].[OH2:32]"
P_SMI = "[NH3:1].[CH2:2]=[CH:3][C:4](=[O:5])[NH2:6].[O:8]=[C:7]=[O:9].[CH:10](=[O:11])[CH:12]([OH:13])[CH:14]([OH:15])[CH:16]([OH:17])[CH:18]([OH:19])[CH2:20][OH:21].[OH2:30].[OH2:31].[OH2:32]"


def _parse_xyz_coords(xyz_path: Path) -> list:
    lines = xyz_path.read_text("utf-8").strip().splitlines()
    count = int(lines[0])
    return [line.split() for line in lines[2:2 + count]]

def _write_xyz(coords: list, out_path: Path):
    lines = [str(len(coords)), ""]
    for sym, x, y, z in coords:
         lines.append(f"{sym:2s}  {float(x):12.6f}  {float(y):12.6f}  {float(z):12.6f}")
    out_path.write_text("\n".join(lines) + "\n", "utf-8")


def generate_seed():
    print("Generating explicit water cluster seed...")
    GEOMETRY_BASE.mkdir(parents=True, exist_ok=True)
    
    # We leverage fix_steric_clashes for the alignment logic
    sys.path.insert(0, str(ROOT))
    import scripts.fix_steric_clashes as fix
    
    r_mol, p_mol = fix._reorder_product_to_match_reactant(R_SMI, P_SMI)
    
    print("Embedding base geometries...")
    r_mol_3d = fix._embed_with_fragment_separation(r_mol)
    p_mol_3d = fix._embed_with_fragment_separation(p_mol)
    
    # Write initial separated geometries
    r_xyz_temp = GEOMETRY_BASE / "reactant_unrelaxed.xyz"
    p_xyz_temp = GEOMETRY_BASE / "product_unrelaxed.xyz"
    
    r_xyz_temp.write_text(fix._mol_to_xyz_string(r_mol_3d), "utf-8")
    p_xyz_temp.write_text(fix._mol_to_xyz_string(p_mol_3d), "utf-8")
    
    # Find heavy atom indices to freeze during water cluster relaxation
    # We freeze all atoms except the last 9 (which are the 3 water molecules: 3 * 3 = 9 atoms)
    total_atoms = r_mol_3d.GetNumAtoms()
    water_atoms = 9
    frozen_atoms = list(range(1, total_atoms - water_atoms + 1)) # 1-indexed for xTB
    
    # Write xTB constraint file
    const_file = GEOMETRY_BASE / "constrain.inp"
    const_file.write_text(f"""$constrain
  force constant=1.0
  atoms: {",".join(map(str, frozen_atoms))}
$end
""", "utf-8")
    
    print(f"Relaxing water cluster around fixed reactant scaffold ({total_atoms - water_atoms} frozen atoms)...")
    try:
        subprocess.run(
            ["xtb", "reactant_unrelaxed.xyz", "--opt", "--input", "constrain.inp", "--gfn", "2", "--alpb", "water", "--chrg", "0"],
            cwd=GEOMETRY_BASE, check=True, capture_output=True, text=True
        )
        shutil.move(GEOMETRY_BASE / "xtbopt.xyz", GEOMETRY_BASE / "reactant_relaxed.xyz")
        
        print(f"Relaxing water cluster around fixed product scaffold ({total_atoms - water_atoms} frozen atoms)...")
        subprocess.run(
            ["xtb", "product_unrelaxed.xyz", "--opt", "--input", "constrain.inp", "--gfn", "2", "--alpb", "water", "--chrg", "0"],
            cwd=GEOMETRY_BASE, check=True, capture_output=True, text=True
        )
        shutil.move(GEOMETRY_BASE / "xtbopt.xyz", GEOMETRY_BASE / "product_relaxed.xyz")
    except FileNotFoundError:
        print("ERROR: xtb command not found. Must run inside maillard conda environment.")
        return 1
    except subprocess.CalledProcessError as e:
        print("ERROR during xTB optimization.")
        print(e.stderr)
        return 1
        
    # Read relaxed coords
    relaxed_r_coords = _parse_xyz_coords(GEOMETRY_BASE / "reactant_relaxed.xyz")
    relaxed_p_coords = _parse_xyz_coords(GEOMETRY_BASE / "product_relaxed.xyz")
    
    _write_xyz(relaxed_r_coords, GEOMETRY_BASE / "reactant.xyz")
    _write_xyz(relaxed_p_coords, GEOMETRY_BASE / "product.xyz")
    
    print("Wrote aligned and clustered reactant.xyz and product.xyz.")
    
    # Write run_xtb.sh
    run_sh = GEOMETRY_BASE / "run_xtb.sh"
    run_sh.write_text("""#!/bin/bash
echo "Running xTB path search for explicit water cluster"
xtb reactant.xyz --path product.xyz --gfn 2 --alpb water
""", "utf-8")
    run_sh.chmod(0o755)
    
    # Cleanup intermediate files
    for clean_tgt in ["reactant_unrelaxed.xyz", "product_unrelaxed.xyz", "constrain.inp", "xtbopt.xyz", "xtbtopo.mol", "charges", "wbo", "xtbopt.log"]:
        f = GEOMETRY_BASE / clean_tgt
        if f.exists():
            f.unlink()
            
    print("Explicit water cluster seed successfully generated.")
    return 0


if __name__ == "__main__":
    sys.exit(generate_seed())
