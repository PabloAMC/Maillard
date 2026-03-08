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

r1 = get_comp("NC(CS)C(=O)O") # cys
r2 = get_comp("O") # water
print(f"Reactants: Cys {r1} + Water {r2}")

p1 = get_comp("CC(=O)C=O") # pyruv
p2 = get_comp("N") # nh3
p3 = get_comp("O=C=O") # co2
p4 = get_comp("S") # h2s
print(f"Products: Pyruv {p1} + NH3 {p2} + CO2 {p3} + H2S {p4}")
