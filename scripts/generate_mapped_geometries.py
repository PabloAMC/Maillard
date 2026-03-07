import sys
import types
import pathlib

# Monkeypatch pkg_resources for rxnmapper in Python 3.14
dummy_pkg = types.ModuleType("pkg_resources")
def resource_filename(package_or_requirement, resource_name):
    # We know rxnmapper looks for models/transformers/...
    # We will resolve it relative to the rxnmapper package installation directory
    import importlib.util
    spec = importlib.util.find_spec("rxnmapper")
    base_path = pathlib.Path(spec.origin).parent
    return str(base_path / resource_name)
dummy_pkg.resource_filename = resource_filename
sys.modules["pkg_resources"] = dummy_pkg

from rxnmapper import RXNMapper
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from pathlib import Path

def generate_mapped_pair(rxnmapper, name, reactant_smiles, product_smiles, output_dir):
    print(f"[{name}] Generating mapped geometries...")
    rxn_smiles = f"{reactant_smiles}>>{product_smiles}"
    
    # 1. Map heavy atoms with IBM rxnmapper
    try:
        results = rxnmapper.get_attention_guided_atom_maps([rxn_smiles])
        mapped_rxn = results[0]['mapped_rxn']
    except Exception as e:
        print(f"  [Error] rxnmapper failed: {e}")
        return False

    r_smi, p_smi = mapped_rxn.split('>>')
    from rdkit import Geometry
    
    # Generate separated reactant
    frags = r_smi.split('.')
    r_mol_separated = None
    offset = 0.0
    for frag_smi in frags:
        frag = Chem.MolFromSmiles(frag_smi)
        frag = Chem.AddHs(frag)
        AllChem.EmbedMolecule(frag, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(frag, maxIters=1000)
        
        # Translate
        conf = frag.GetConformer()
        for i in range(frag.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + offset, pos.y, pos.z))
            
        if r_mol_separated is None:
            r_mol_separated = frag
        else:
            r_mol_separated = Chem.CombineMols(r_mol_separated, frag)
            
        offset += 3.0 # 3.0 A is a safe distance that prevents overlapping but allows H-bonding

    # Relax the combined VdW complex so fragments dock naturally!
    Chem.GetSSSR(r_mol_separated)
    AllChem.MMFFOptimizeMolecule(r_mol_separated, maxIters=5000)

    r_mol = r_mol_separated
    p_mol = Chem.MolFromSmiles(p_smi)
    
    if not r_mol or not p_mol:
        print("  [Error] Failed to parse mapped SMILES.")
        return False

    # 2. Add Hydrogens to product (r_mol already has Hs and 3D coords!)
    p_mol = Chem.AddHs(p_mol)
    
    # Extract map indices for heavy atoms
    r_map_to_idx = {a.GetAtomMapNum(): a.GetIdx() for a in r_mol.GetAtoms() if a.GetAtomMapNum() > 0}
    p_map_to_idx = {a.GetAtomMapNum(): a.GetIdx() for a in p_mol.GetAtoms() if a.GetAtomMapNum() > 0}
    
    ordered_maps = [a.GetAtomMapNum() for a in r_mol.GetAtoms() if a.GetAtomMapNum() > 0]
    unmapped_p = [a.GetIdx() for a in p_mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetAtomMapNum() == 0]
    
    # 3. Map Hydrogens
    p_h_indices = {a.GetIdx() for a in p_mol.GetAtoms() if a.GetAtomicNum() == 1}
    r_h_indices = {a.GetIdx() for a in r_mol.GetAtoms() if a.GetAtomicNum() == 1}
    
    unmatched_p_h = set(p_h_indices)
    unmatched_r_h = set(r_h_indices)
    p_to_r_h = {}
    
    for m in ordered_maps:
        if m not in p_map_to_idx: continue
        r_idx = r_map_to_idx[m]
        p_idx = p_map_to_idx[m]
        
        r_hs = [n.GetIdx() for n in r_mol.GetAtomWithIdx(r_idx).GetNeighbors() if n.GetAtomicNum() == 1]
        p_hs = [n.GetIdx() for n in p_mol.GetAtomWithIdx(p_idx).GetNeighbors() if n.GetAtomicNum() == 1]
        
        matched_count = min(len(r_hs), len(p_hs))
        for i in range(matched_count):
            p_to_r_h[p_hs[i]] = r_hs[i]
            unmatched_p_h.remove(p_hs[i])
            unmatched_r_h.remove(r_hs[i])
            
    unmatched_p_h = list(unmatched_p_h)
    unmatched_r_h = list(unmatched_r_h)
    for p_h, r_h in zip(unmatched_p_h, unmatched_r_h):
        p_to_r_h[p_h] = r_h
        
    r_to_p_h = {v: k for k, v in p_to_r_h.items()}
    
    # 4. Construct final ordering for Product
    final_p_order = []
    for i in range(r_mol.GetNumAtoms()):
        a = r_mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() > 1:
            m = a.GetAtomMapNum()
            if m > 0 and m in p_map_to_idx:
                final_p_order.append(p_map_to_idx[m])
            else:
                if unmapped_p: final_p_order.append(unmapped_p.pop(0))
                else: return False
        else:
            final_p_order.append(r_to_p_h.get(i, -1))
            
    if -1 in final_p_order: return False
        
    # 5. Embed Product
    try:
        p_mol_aligned = Chem.RenumberAtoms(p_mol, final_p_order)
        # Use reactant coordinates as a starting point.
        # Since components are 6A apart, the bonded product will be highly strained.
        # MMFF94 will smoothly pull them together into the realistic bonded shape.
        conf = r_mol.GetConformer()
        p_conf = Chem.Conformer(r_mol.GetNumAtoms())
        for i in range(r_mol.GetNumAtoms()):
            p_conf.SetAtomPosition(i, conf.GetAtomPosition(i))
        p_mol_aligned.AddConformer(p_conf)
        
        AllChem.MMFFOptimizeMolecule(p_mol_aligned, maxIters=5000)
    except Exception as e:
        print(f"  [Error] Failed to embed/optimize: {e}")
        return False
        
    # 6. Save
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "reactant.xyz", "w") as f:
        f.write(Chem.MolToXYZBlock(r_mol))
    with open(output_dir / "product.xyz", "w") as f:
        f.write(Chem.MolToXYZBlock(p_mol_aligned))
    
    print(f"  [Success] Saved perfectly aligned XYZ to {output_dir}")
    return True

def main():
    base_dir = Path("data/geometries/xtb_inputs")
    rxnmapper = RXNMapper()
    
    # Must be perfectly mass-balanced for mapping to succeed
    REACTIONS = {
        "enolisation": {
            "r": "OCC(O)C(O)C(=O)CNCC(=O)O", # Amadori
            "p": "O=CC(=O)CC(O)CO.NCC(=O)O"  # 3-Deoxyosone + Glycine (balanced)
        },
        "strecker": {
            "r": "CC(=O)C=O.CC(C)CC(N)C(=O)O", # Pyruvaldehyde + Leucine
            "p": "CC(C)CC=O.CC(=O)CN.O=C=O"    # 3-Methylbutanal + Aminoacetone + CO2
        },
        "cys_ribose": {
            "r": "O=Cc1ccco1.S.[H][H]",        # Furfural + H2S + H2 (reducing equivalent!)
            "p": "SCc1ccco1.O"                 # FFT + Water
        },
        "trapping": {
            "r": "CCCCCC=O.NCC(=O)O",          # Hexanal + Glycine
            "p": "CCCCC/C=N/CC(=O)O.O"         # Schiff base + water
        },
        "dha": {
            "r": "NC(CS)C(=O)O",               # Cysteine
            "p": "C=C(N)C(=O)O.S"              # DHA + H2S
        },
        "pyrazine": {
            "r": "CC(=O)CN.CC(=O)CN",          # 2x Aminoacetone
            "p": "CC1=NCC(C)=NC1.O.O"          # Dihydropyrazine + 2x Water (balanced)
        },
        "retro_aldol": {
            "r": "OCC(O)C(O)C(O)C=O",          # Linear Ribose
            "p": "OCC(O)C=O.OCC=O"             # Glyceraldehyde + Glycolaldehyde
        }
    }

    for name, s in REACTIONS.items():
        r_dir = base_dir / name
        generate_mapped_pair(rxnmapper, name, s['r'], s['p'], r_dir)

if __name__ == "__main__":
    main()
