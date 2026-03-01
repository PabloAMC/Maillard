#!/usr/bin/env python3
"""
scripts/run_tier2_dft.py — Execute high-accuracy Skala/DFT barrier refinement.

This script demonstrates processing the 8 critical rate-limiting steps 
(defined in Phase 3.3) using the SkalaRefiner module. 
Since true DFT TS scanning (especially of large molecules) is computationally 
expensive (hours/days), this acts as the control logic for HPC deployment.
"""

from src.skala_refiner import SkalaRefiner

# We mock the initial XYZ coordinate blocks for the Amadori Rearrangement
# In a real workflow, these geometries are passed from the Tier 1 xTB outputs.

# Using Methanol + Explicit Water for a fast, closed-shell mock 
# capable of running without odd-electron spin errors, representing
# explicit solvation/Grotthuss mechanism setup.
MOCK_REACTANT_XYZ = """9

C 0.000000 0.000000 0.000000
O 1.400000 0.000000 0.000000
H 1.700000 0.900000 0.000000
H 0.300000 -0.500000 0.800000
H -0.300000 -0.500000 -0.800000
H -0.500000 0.900000 0.000000
O -2.500000 0.000000 0.000000
H -1.500000 0.000000 0.000000
H -2.800000 0.800000 0.000000
"""

# Stretched C-H...O bond for mock TS with explicit water shuttle
MOCK_TS_XYZ = """9

C 0.000000 0.000000 0.000000
O 1.400000 0.000000 0.000000
H 1.700000 0.900000 0.000000
H 0.300000 -0.500000 0.800000
H -0.300000 -0.500000 -0.800000
H -0.800000 0.500000 0.000000
O -1.800000 0.000000 0.000000
H -1.600000 -0.900000 0.000000
H -2.500000 0.400000 0.000000
"""

def refine_core_reactions():
    print("Initializing Tier 2 Skala DFT Refiner...")
    # Use standard basis set def2-SVP for real runs.
    # We pass use_skala=True (assuming module installed on HPC). 
    # It safely falls back to B3LYP if not present on the current runner.
    refiner = SkalaRefiner(basis='def2-svp', solvent_name='water', use_skala=True)
    
    tasks = {
        "3.3a": "Amadori rearrangement (Schiff base -> 1-amino-1-deoxy-2-ketose)",
        "3.3b": "2,3-enolisation vs 1,2-enolisation bifurcation point",
        "3.3c": "Strecker decarboxylation (a-dicarbonyl + amino acid)",
        "3.3d": "Cysteine + ribose -> FFT (via furfural + H2S)",
        "3.3e": "Ribose retro-aldol -> 1,4-dideoxyosone -> MFT",
        "3.3f": "DHA b-elimination (Ser -> dehydroalanine)",
        "3.3g": "Off-flavour trapping: hexanal + amino acid -> Schiff base",
        "3.3h": "a-aminoketone dimerisation -> pyrazine"
    }
    
    print("\nQueuing the following 8 critical barriers for HPC refinement:")
    for task_id, desc in tasks.items():
        print(f"  [{task_id}] {desc}")
        
    print("\n--- Dry Run: Submitting 3.3a (Amadori Refinement) ---")
    try:
        # In this mock demo run, we will just use single points because an 
        # actual DFT TS optimization takes a long time.
        r_energy = refiner.single_point(MOCK_REACTANT_XYZ).energy_hartree
        ts_energy = refiner.single_point(MOCK_TS_XYZ).energy_hartree
        
        barrier_kcal = (ts_energy - r_energy) * 627.509
        print(f"Calculated Barrier for 3.3a: {barrier_kcal:.2f} kcal/mol")
        
        if barrier_kcal > 100 or barrier_kcal < -100:
            print("Note: The mock XYZ geometries generate unphysical barriers. This is expected without true xTB pre-optimization.")
            
    except Exception as e:
        print(f"Execution skipped - PySCF/HPC requirements missing locally. Error: {e}")

if __name__ == "__main__":
    refine_core_reactions()
