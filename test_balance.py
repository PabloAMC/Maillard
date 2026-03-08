from data.reactions.curated_pathways import PATHWAYS
from rdkit import Chem

def get_comp(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol: return None
    mol = Chem.AddHs(mol)
    comp = {}
    for a in mol.GetAtoms():
        s = a.GetSymbol()
        comp[s] = comp.get(s, 0) + 1
    return comp

for name, steps in PATHWAYS.items():
    for step in steps:
        r_comp = {}
        for s in step.reactants:
            c = get_comp(s.smiles)
            if not c: print(f"Failed to parse {s.smiles}"); continue
            for k,v in c.items(): r_comp[k] = r_comp.get(k, 0) + v
        p_comp = {}
        for s in step.products:
            c = get_comp(s.smiles)
            if not c: print(f"Failed to parse {s.smiles}"); continue
            for k,v in c.items(): p_comp[k] = p_comp.get(k, 0) + v
        if r_comp != p_comp:
            print("UNBALANCED:")
            print("Reactants:", [s.smiles for s in step.reactants], r_comp)
            print("Products:", [s.smiles for s in step.products], p_comp)

