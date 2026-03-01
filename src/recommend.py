#!/usr/bin/env python3
"""
src/recommend.py — Maillard Reaction Pathway Recommender

Prototype integration of the three-tier computational framework:
Tier 0: Graph traversal (mocked RMG-Py here, as RMG requires a separate environment)
Tier 1: xTB Energy Screening (GFN2-xTB)
Tier 2: DFT Refinement (PySCF + Skala)

Outputs an ordered list of viable pathways and estimates sensory profile impact
based on thermodynamic favourability under given environmental conditions.
"""

import sys
from pathlib import Path

# Add project root to python path if run directly
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.pathway_extractor import PathwayExtractor, Species, ElementaryStep
from src.xtb_screener import XTBScreener
from src.pathway_ranker import PathwayRanker
from src.conditions import ReactionConditions
from src.skala_refiner import SkalaRefiner

def main():
    print("======================================================")
    print("      Maillard Reaction Pathway Recommender")
    print("======================================================\n")
    
    # 1. Define conditions
    cond = ReactionConditions(pH=6.5, temperature_celsius=140.0, water_activity=0.7)
    print(f"[1] Setting Conditions: pH {cond.pH:.1f}, Temp {cond.temperature_celsius:.1f} °C, aw {cond.water_activity:.2f}")
    
    # 2. Define precursors
    print(f"[2] Initializing Precursors: D-Ribose + Cysteine (Meat/Savory model)")
    ribose = Species("D-Ribose", "O=C[C@H](O)[C@H](O)[C@H](O)CO")
    cysteine = Species("L-Cysteine", "N[C@@H](CS)C(=O)O")
    
    # 3. Extract Pathways (Tier 0)
    print(f"\n[Tier 0] Extracting Pathways (Graph/RMG)...")
    extractor = PathwayExtractor(rmg_output_dir="mock/path/to/rmg") 
    
    # Mocking Tier 0 output for the prototype (since RMG is a heavy external dep)
    mock_schiff = Species("SchiffBase", "OC[C@@H](O)[C@@H](O)[C@@H](O)/C=N/[C@@H](CS)C(=O)O")
    mock_water = Species("Water", "O")
    mock_arp = Species("AmadoriProduct", "OC[C@@H](O)[C@@H](O)[C@@H](O)C(=O)CN[C@@H](CS)C(=O)O")
    mock_fft = Species("2-Furfurylthiol", "SCc1ccco1")
    
    step1 = ElementaryStep(
        reactants=[ribose, cysteine], 
        products=[mock_schiff, mock_water], 
        reaction_family="Schiff Base Formation"
    )
    step2 = ElementaryStep(
        reactants=[mock_schiff], 
        products=[mock_arp], 
        reaction_family="Amadori Rearrangement"
    )
    step3 = ElementaryStep(
        reactants=[mock_arp],
        products=[mock_fft],
        reaction_family="Furan-thiol Synthesis" # simplified multiple steps
    )
    
    pathways = {
        "Savory_FFT_Pathway": [step1, step2, step3]
    }
    
    # 4. Screen with xTB (Tier 1)
    print(f"\n[Tier 1] xTB Screening ({len(pathways)} pathways)...")
    screener = XTBScreener()
    ranker = PathwayRanker(n_cores=1, conditions=cond)
    
    # In a real run, XTBScreener would process each step in pathways to get `barrier_kcal_list`.
    # PathwayRanker currently mocks evaluate_single_step if xtb is missing.
    ranked_profiles = ranker.screen_pathways(pathways)
    
    for prof in ranked_profiles:
        print(f"  -> Analyzed {prof.pathway_name}: Span = {prof.energetic_span:.1f} kcal/mol, Thermo = {prof.overall_thermodynamics:.1f} kcal/mol")
        
    top_pathway = ranked_profiles[0]
        
    # 5. Refine with Skala (Tier 2)
    print(f"\n[Tier 2] DFT/Skala Refinement of Top Candidates...")
    print(f"  Selecting rate-limiting step from {top_pathway.pathway_name} for high-accuracy DFT...")
    refiner = SkalaRefiner(basis='def2-svp', use_skala=True)
    # Mock TS refinement geometry
    mock_ts_xyz = "5\n\nC 0.0 0.0 0.0\nH 1.2 1.2 1.2\nH -0.6 -0.6 0.6\nH -0.6 0.6 -0.6\nH 0.6 -0.6 -0.6"
    mock_r_xyz = "5\n\nC 0.0 0.0 0.0\nH 0.6 0.6 0.6\nH -0.6 -0.6 0.6\nH -0.6 0.6 -0.6\nH 0.6 -0.6 -0.6"
    
    try:
        # Fast dry-run
        refiner.fallback_xc = 'lda'
        refiner.basis = 'sto-3g'
        refined_barrier = refiner.refine_barrier(mock_r_xyz, mock_ts_xyz)
        print(f"  Refined Barrier: {refined_barrier:.1f} kcal/mol (using {refiner.fallback_xc}/{refiner.basis} fallback)")
    except Exception as e:
        print(f"  Skipping refinement, PySCF/Skala setup incomplete: {e}")
        
    # 6. Output Recommendations
    print(f"\n======================================================")
    print(f" Final Recommendation: {top_pathway.pathway_name} is highly favoured.")
    print(f" Expect formation of Savory/Meat notes (FFT).")
    print(f"======================================================")

if __name__ == "__main__":
    main()
