#!/usr/bin/env python3
"""
scripts/benchmark_mace_rmsd.py

Phase 10.6: VALIDATION GATE
Quantifies the out-of-distribution error of the MACE foundation model natively on
Maillard-specific geometries by calculating the Root Mean Square Deviation (RMSD)
between a MACE-optimized structure and a DFT-optimized source ground truth.

Used to evaluate when fine-tuning has achieved the target < 0.1 Å tolerance required
for SOTA architecture alignment.
"""

import argparse
import numpy as np

try:
    from ase.io import read
    from ase.build import minimize_rotation_and_translation
    from src.mlp_optimizer import MLPOptimizer
except ImportError:
    print("Warning: ASE or src.mlp_optimizer not found. Ensure mace-torch is installed.")

def calculate_rmsd(xyz_ref_string: str, xyz_test_string: str) -> float:
    """ Calculate RMSD between two XYZ coordinate sets using ASE mapping """
    import tempfile
    
    with tempfile.NamedTemporaryFile("w+", suffix=".xyz") as f1, \
         tempfile.NamedTemporaryFile("w+", suffix=".xyz") as f2:
        f1.write(xyz_ref_string)
        f1.flush()
        f2.write(xyz_test_string)
        f2.flush()
        
        atoms_ref = read(f1.name, format="xyz")
        atoms_test = read(f2.name, format="xyz")
        
    if len(atoms_ref) != len(atoms_test):
        raise ValueError("Structures have a different number of atoms!")
        
    # Align structures to minimize rigid body translation/rotation differences
    minimize_rotation_and_translation(atoms_ref, atoms_test)
    
    pos_ref = atoms_ref.get_positions()
    pos_test = atoms_test.get_positions()
    
    diff = pos_ref - pos_test
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return float(rmsd)

def benchmark_reaction_geometry(dft_xyz_path: str, model_name: str = "medium"):
    """
    Takes a confirmed DFT structure, passes it through MACE optimization,
    and returns the resulting RMSD drift.
    """
    print(f"Loading '{model_name}' MACE model...")
    optimizer = MLPOptimizer(model_name=model_name, device="cpu")
    
    with open(dft_xyz_path, "r") as f:
        dft_xyz = f.read()
        
    print(f"Running MACE {model_name} geometry relaxation...")
    mace_xyz = optimizer.optimize_geometry(dft_xyz, fmax=0.01)
    
    rmsd = calculate_rmsd(dft_xyz, mace_xyz)
    print("====================================")
    print(f"DFT to MACE RMSD Drift: {rmsd:.4f} Å")
    print("====================================")
    
    if rmsd > 0.1:
        print("[WARNING] OUT-OF-DISTRIBUTION ERROR.")
        print("MACE structure differs from DFT by more than 0.1 Angstroms.")
        print("The network must be fine-tuned via mace-torch before production use.")
    else:
        print("[SUCCESS] MACE represents the ground truth structure securely (< 0.1 Å).")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark MACE structural accuracy against DFT.")
    parser.add_argument("--dft_xyz", type=str, required=True, help="Path to the reference DFT .xyz file")
    parser.add_argument("--model", type=str, default="small", choices=["small", "medium", "large"],
                        help="The mace-mp-0 model scale to benchmark (default: small)")
    
    args = parser.parse_args()
    benchmark_reaction_geometry(args.dft_xyz, model_name=args.model)
