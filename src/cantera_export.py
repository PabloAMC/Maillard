"""
src/cantera_export.py

Phase 12: Cantera mechanism export.
Converts the SMILES-based reaction network and computed barriers 
into a Cantera-compatible YAML mechanism file.
"""

import yaml
from typing import Dict, List, Optional, Any
from pathlib import Path
from src.kinetics import KineticsEngine

class CanteraExporter:
    """
    Exports Maillard reaction networks to Cantera YAML format.
    """
    
    def __init__(self, name: str = "Maillard_Mechanism"):
        self.name = name
        self.species = {}  # SMILES -> Metadata
        self.reactions = []
        self.kinetics = KineticsEngine()

    def add_species(self, smiles: str, name: Optional[str] = None):
        """Append a species to the mechanism and calculate composition via RDKit."""
        if smiles not in self.species:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                comp = {}
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    comp[symbol] = comp.get(symbol, 0) + 1
            else:
                # Fallback composition when SMILES is not parseable
                comp = {"C": 1, "H": 1, "O": 1}
                
            self.species[smiles] = {
                "name": name or f"S_{len(self.species)}",
                "smiles": smiles,
                "composition": comp,
            }
        return self.species[smiles]["name"]

    def add_reaction(self, reactants: List[str], products: List[str], 
                     barrier_kcal: float, temperature_k: float = 423.15,
                     A: float = 1e13, b: float = 0.0,
                     thermo_gating: bool = True,
                     gating_threshold: float = 30.0,
                     reaction_family: Optional[str] = None,
                     conditions: Optional[Any] = None):
        """
        Add an elementary step. Enforces strict mass balance and optional thermo-gating.
        """
        # 1. Thermo-Gating (Phase 12.3)
        if thermo_gating:
            is_feasible, dg = self.kinetics.is_reaction_feasible(reactants, products, gating_threshold)
            if not is_feasible:
                # Log and skip unphysical reactions
                # print(f"  Skipping unphysical reaction (Delta G = {dg:.1f} kcal/mol): {reactants} -> {products}")
                return

        # Apply pH/Environment multipliers (Phase 16.2 fix)
        if conditions and reaction_family:
            A *= conditions.get_ph_multiplier(reaction_family)
            A *= conditions.get_water_activity_multiplier()

        r_names = [self.add_species(r) for r in reactants]
        p_names = [self.add_species(p) for p in products]
        
        # Compute total reactant composition
        r_comp = {}
        for r in reactants:
            for k, v in self.species[r]["composition"].items():
                r_comp[k] = r_comp.get(k, 0) + v
        
        p_comp = {}
        for p in products:
            for k, v in self.species[p]["composition"].items():
                p_comp[k] = p_comp.get(k, 0) + v

        # Enforce strict balance
        all_elements = set(r_comp.keys()) | set(p_comp.keys())
        diffs = {}
        for elem in all_elements:
            d = p_comp.get(elem, 0) - r_comp.get(elem, 0)
            if d != 0:
                diffs[elem] = d

        if diffs:
            raise ValueError(f"Reaction is structurally unbalanced! "
                             f"Reactants {reactants} vs Products {products}. "
                             f"Atom difference (Products - Reactants): {diffs}")

        # Support for stoichiometric coefficients in equation string
        def group_and_count(names):
            counts = {}
            for n in names:
                counts[n] = counts.get(n, 0) + 1
            return " + ".join([(f"{v} " if v > 1 else "") + k for k, v in counts.items()])

        reaction_str = f"{group_and_count(r_names)} <=> {group_and_count(p_names)}"
        
        # Scale A-factor for trimolecular reactions to reduce stiffness (Phase 16.3 fix)
        # Standard collision theory: A_tri << A_bi
        if len(reactants) >= 3:
            A = A * 1e-4

        self.reactions.append({
            "equation": reaction_str,
            "rate-constant": {
                "A": A,
                "b": b,
                "Ea": f"{barrier_kcal} kcal/mol"
            }
        })

    def export_yaml(self, output_path: str):
        """Write the Cantera YAML file."""
        elements = set()
        for s in self.species.values():
            elements.update(s["composition"].keys())
            
        data = {
            "generator": "Maillard Framework CanteraExporter",
            "cantera-version": "3.0.0",
            "date": "2026-03-08",
            "units": {"length": "m", "quantity": "kmol", "mass": "kg", "time": "s", "energy": "J"},
            "phases": [{
                "name": "maillard_phase",
                "thermo": "ideal-gas", 
                "species": [s["name"] for s in self.species.values()],
                "kinetics": "gas", 
                "reactions": "all",
                "state": {"T": 423.15, "P": 101325}
            }],
            "species": [],
            "reactions": self.reactions
        }
        
        from src.thermo import get_nasa_coefficients
        
        for s in self.species.values():
            try:
                coeffs = get_nasa_coefficients(s["smiles"])
                thermo_block = {
                    "model": "NASA7",
                    "temperature-ranges": [300.0, 1000.0],
                    "data": [coeffs]
                }
            except Exception as e:
                # Fallback to constant-cp if estimation fails
                print(f"  Warning: Thermo estimation failed for {s['name']} ({s['smiles']}): {e}")
                thermo_block = {
                    "model": "constant-cp",
                    "t0": 273.15,
                    "h0": 0.0,
                    "s0": 0.0,
                    "cp0": 75000.0
                }

            data["species"].append({
                "name": s["name"],
                "composition": s["composition"],
                "thermo": thermo_block
            })
            
        with open(output_path, "w") as f:
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)
            
        print(f"Exported Cantera mechanism to {output_path}")

if __name__ == "__main__":
    exporter = CanteraExporter()
    exporter.add_reaction(["C(C1C(C(C(C(O1)O)O)O)O)=O"], ["Product_SMILES"], 25.0)
    exporter.export_yaml("test_mech.yaml")
