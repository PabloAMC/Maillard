#!/usr/bin/env python3
"""
generate_rmg_inputs.py — Generate RMG-Py input files for Maillard validation cases.

This creates the `input.py` configuration files for three critical test cases:
1. D-Ribose + L-Cysteine   (Target: FFT, MFT)
2. D-Glucose + Glycine     (Target: Furfural, Amadori intermediates)
3. Ribose + Cys + Leucine  (Target: FFT, MFT, 3-methylbutanal, pyrazines)

Usage: python scripts/generate_rmg_inputs.py
"""

from pathlib import Path

ROOT = Path(__file__).parent.parent
CASES_DIR = ROOT / "data/reactions/rmg_validation_cases"

# Standard RMG database configuration for aqueous-phase Maillard
DB_TEMPLATE = """database(
    thermoLibraries=['primaryThermoLibrary', 'thermo_DFT_CCSDTF12_BAC', 'CBS_QB3_1dHR'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    # Using families relevant for generic aqueous phase and condensation
    kineticsFamilies=['Schiff_Base_Formation', 'Amadori_Rearrangement', 'Heyns_Rearrangement', 'Enolisation', 'Retro_Aldol_Fragmentation', 'Thiol_Addition', 'Strecker_Degradation', 'Aminoketone_Condensation', 'Cysteine_Degradation', 'Lipid_Schiff_Base', 'Lipid_Thiazole_Condensation', 'Beta_Elimination', 'DHA_Crosslinking', 'Sugar_Ring_Opening'],
    kineticsEstimator='rate rules',
)
"""

REACTOR_TEMPLATE = """
# Liquid-phase reactor (150°C, typical Maillard cooking/extrusion temperature)
liquidReactor(
    temperature=(423.15,'K'),
    initialConcentrations=[CONCENTRATIONS_PLACEHOLDER],
    terminationConversion=[TERMINATION_PLACEHOLDER],
    terminationTime=(10,'hr'),
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000,
)

options(
    saveEdgeSpecies=False,
    saveSimulationProfiles=True,
    generatePlots=True,
)
"""

CASES = {
    "case_1_ribose_cys": {
        "species": [
            ("ribose",   "O=CC(O)C(O)C(O)CO", "1.0"),  # Open-chain aldehyde form
            ("cysteine", "N[C@@H](CS)C(=O)O", "1.0"),
            ("water",    "O", "55.0"), # Solvent
        ],
        "termination": "'ribose': 0.9"
    },
    "case_2_glucose_gly": {
        "species": [
            ("glucose", "O=CC(O)C(O)C(O)C(O)CO", "1.0"),  # Open-chain aldehyde form
            ("glycine", "NCC(=O)O", "1.0"),
            ("water",   "O", "55.0"),
        ],
        "termination": "'glucose': 0.9"
    },
    "case_3_ribose_cys_leu": {
        "species": [
            ("ribose",   "O=CC(O)C(O)C(O)CO", "1.0"),  # Open-chain aldehyde form
            ("cysteine", "N[C@@H](CS)C(=O)O", "0.5"),
            ("leucine",  "CC(C)C[C@@H](N)C(=O)O", "0.5"),
            ("water",    "O", "55.0"),
        ],
        "termination": "'ribose': 0.9"
    }
}

def generate_species_block(species_list):
    blocks = []
    concs = []
    for name, smiles, conc in species_list:
        blocks.append(f"species(\n    label='{name}',\n    reactive=True,\n    structure=SMILES('{smiles}')\n)")
        concs.append(f"'{name}': ({conc}, 'mol/l')")
    return "\n\n".join(blocks), ", ".join(concs)

def main():
    CASES_DIR.mkdir(parents=True, exist_ok=True)
    
    for case_name, data in CASES.items():
        case_path = CASES_DIR / case_name
        case_path.mkdir(exist_ok=True)
        
        species_block, concs_str = generate_species_block(data["species"])
        reactor_block = REACTOR_TEMPLATE.replace("[CONCENTRATIONS_PLACEHOLDER]", "{" + concs_str + "}")
        reactor_block = reactor_block.replace("[TERMINATION_PLACEHOLDER]", "{" + data["termination"] + "}")
        
        input_py_content = f"{DB_TEMPLATE}\n{species_block}\n{reactor_block}"
        
        with open(case_path / "input.py", "w") as f:
            f.write(input_py_content)
            
        print(f"Generated {case_path.relative_to(ROOT)}/input.py")

if __name__ == "__main__":
    main()
