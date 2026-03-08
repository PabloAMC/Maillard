from data.reactions.curated_pathways import PATHWAYS
from rdkit import Chem
def get_comp(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if not mol: return None
        mol = Chem.AddHs(mol)
        comp = {}
        for a in mol.GetAtoms():
            s = a.GetSymbol()
            comp[s] = comp.get(s, 0) + 1
        return comp
    except:
        return None

for name, steps in PATHWAYS.items():
    for i, step in enumerate(steps):
        r_comp = {}
        for s in step.reactants:
            c = get_comp(s.smiles)
            if c:
                for k,v in c.items(): r_comp[k] = r_comp.get(k, 0) + v
        p_comp = {}
        for s in step.products:
            c = get_comp(s.smiles)
            if c:
                for k,v in c.items(): p_comp[k] = p_comp.get(k, 0) + v
        if r_comp != p_comp:
            print(f"[{name}] Step {i+1} ({step.reaction_family}) is UNBALANCED:")
            
            all_elements = set(r_comp.keys()) | set(p_comp.keys())
            diffs = {}
            for e in all_elements:
                d = p_comp.get(e,0) - r_comp.get(e,0)
                if d != 0: diffs[e] = d
            
            print(f"  Missing in products: { {k:-v for k,v in diffs.items() if v < 0} }")
            print(f"  Extra in products: { {k:v for k,v in diffs.items() if v > 0} }")
