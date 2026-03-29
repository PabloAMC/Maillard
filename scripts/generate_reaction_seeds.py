import os
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path

# Reaction targets and their reactant/product SMILES (atom-mapped where possible)
# format: (name, reactant_smiles, product_smiles)
REACTIONS = [
    ("hexanal_radical_quench", "CCCCCC=O.[OH]", "CCCCCC([O])[OH]"),
    ("mft_protein_noncovalent", "CC1=C(C=CO1)S.C", "CC1=C(C=CO1)SC"), # simplified methyl-thiol proxy
    ("quinone_cys_michael", "O=C1C=CC(=O)C=C1.CS", "O=C1C(SC)=CC(=O)C=C1"),
    ("quinone_lys_schiff", "O=C1C=CC(=O)C=C1.CN", "CN=C1C=CC(=O)C=C1"),
    ("aa_ring_open_dicarbonyl", "O=C1C(=O)C(O)C(=O)O1", "O=C(CO)C(O)C(=O)C(=O)O"), # simplified DHA ring-open
    ("pe_schiff_base", "CN.CCCCCC=O", "CN=CCCCCCC"),
    ("pe_amadori", "CN.CCCCCC=O", "CN(C)CCCCCC=O"), # simplified Amadori proxy
    ("melanoidin_radical_trapping", "C=CC=O.[S]", "CSCC=O"), # simplified vinylic radical trapping
    ("lysinoalanine_crosslink", "NC(C)C(=O)O.NC(C)C(=O)O", "NC(CSC)C(=O)O"),
    ("furosine_amadori_hydrolysis", "C1=CC=C(O1)CN.N", "C1=CC=C(O1)CO.N") # simplified furosine proxy
]

def generate_xyz(name, r_smiles, p_smiles):
    print(f"Generating 3D seeds for {name}...")
    
    # Reactant
    mol_r = Chem.MolFromSmiles(r_smiles)
    mol_r = Chem.AddHs(mol_r)
    AllChem.EmbedMolecule(mol_r, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_r)
    
    # Product (simplified generation)
    mol_p = Chem.MolFromSmiles(p_smiles)
    mol_p = Chem.AddHs(mol_p)
    AllChem.EmbedMolecule(mol_p, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_p)
    
    # Save to data/geometries/xtb_inputs/{name}/
    out_dir = Path(f"data/geometries/xtb_inputs/{name}")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    Chem.MolToXYZFile(mol_r, str(out_dir / "reactant.xyz"))
    Chem.MolToXYZFile(mol_p, str(out_dir / "product.xyz"))
    print(f"  - Done: {out_dir}")

def main():
    for name, r, p in REACTIONS:
        try:
            generate_xyz(name, r, p)
        except Exception as e:
            print(f"Error generating {name}: {e}")

if __name__ == "__main__":
    main()
