"""
src/cantera_export.py

Phase 12: Cantera mechanism export.
Converts the SMILES-based reaction network and computed barriers 
into a Cantera-compatible YAML mechanism file.

R.1 Fix: Uses `ideal-condensed` phase model for liquid/solid Maillard
kinetics, with per-species molar volumes estimated from molecular weight
and Girolami's density approximation.
"""

import yaml
from typing import Dict, List, Optional, Any

# Scientific stack
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

from src.kinetics import KineticsEngine  # noqa: E402
from src.barrier_constants import get_arrhenius_params  # noqa: E402
from src.thermo import get_nasa_coefficients  # noqa: E402


def _estimate_molar_volume(smiles: str) -> float:
    """
    Estimate molar volume in m³/kmol using Girolami's atom/group contribution method.

    Girolami's method estimates liquid density from atom counts:
      density ≈ MW / (sum of atom volume increments)
    Then: V_m = MW / density (in L/mol) → convert to m³/kmol.

    For organic Maillard intermediates, typical densities are 0.8–1.3 g/cm³.
    We use a simplified version that gives reasonable estimates for small organics.

    We use a simplified version that gives reasonable estimates for small organics.

    Returns molar volume in m³/kmol (Cantera's default unit for this field).
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.1  # Fallback: ~100 mL/mol

    mol = Chem.AddHs(mol)
    mw = Descriptors.MolWt(mol)

    # Count atoms by element
    atom_counts = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        atom_counts[sym] = atom_counts.get(sym, 0) + 1

    # Girolami volume increments (cm³/mol per atom)
    # Source: Girolami, J. Chem. Educ. 1994, 71, 962-964
    vol_increments = {
        "C": 16.35, "H": 8.71, "O": 12.43, "N": 14.39,
        "S": 22.91, "Cl": 20.95, "Br": 26.21, "F": 10.48,
    }

    total_vol = 0.0
    for sym, count in atom_counts.items():
        increment = vol_increments.get(sym, 15.0)  # Default for rare atoms
        total_vol += increment * count

    if total_vol <= 0:
        total_vol = mw  # Absolute fallback

    # total_vol is in cm³/mol = mL/mol
    # Convert to m³/kmol: 1 mL/mol = 1e-3 L/mol = 1e-6 m³/mol = 1e-3 m³/kmol
    molar_volume_m3_per_kmol = total_vol * 1e-3

    return max(molar_volume_m3_per_kmol, 0.01)  # Floor at 10 mL/mol


class CanteraExporter:
    """
    Exports Maillard reaction networks to Cantera YAML format.
    """
    
    def __init__(self, name: str = "Maillard_Mechanism"):
        self.name: str = name
        self.species: Dict[str, Dict[str, Any]] = {}  # SMILES -> Metadata
        self.reactions: List[Dict[str, Any]] = []
        self.kinetics: KineticsEngine = KineticsEngine()

    def add_species(self, smiles: str, name: Optional[str] = None) -> str:
        """Append a species to the mechanism and calculate composition via RDKit."""
        if smiles not in self.species:
            
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
                "molar_volume": _estimate_molar_volume(smiles),
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
                return

        # 2. Literature Arrhenius Calibration (Phase A)
        source_quality = "heuristic"
        if reaction_family:
            lit_params = get_arrhenius_params(reaction_family)
            if lit_params:
                # lit_params is (A, Ea, source, quality)
                A_lit, barrier_lit, source_q_lit, _ = lit_params
                A = float(A_lit)
                barrier_kcal = float(barrier_lit)
                source_quality = str(source_q_lit)

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

        reaction_str = f"{group_and_count(r_names)} => {group_and_count(p_names)}"
        
        # Scale A-factor for trimolecular reactions to reduce stiffness (Phase 16.3 fix)
        if len(reactants) >= 3:
            A = float(A) * 1e-4

        self.reactions.append({
            "equation": reaction_str,
            "rate-constant": {
                "A": A,
                "b": b,
                "Ea": f"{barrier_kcal} kcal/mol"
            },
            "note": f"Source: {source_quality} ({reaction_family})"
        })

    def export_yaml(self, output_path: str):
        """Write the Cantera YAML file."""
        elements = set()
        for s in self.species.values():
            elements.update(s["composition"].keys())
            
        # Explicit typing helps some linters with the heterogeneous dict
        data: Dict[str, Any] = {
            "generator": "Maillard Framework CanteraExporter",
            "cantera-version": "3.0.0",
            "date": "2026-03-11",
            "units": {"length": "m", "quantity": "kmol", "mass": "kg", "time": "s", "energy": "J"},
            "phases": [{
                "name": "maillard_phase",
                "thermo": "ideal-condensed",
                "species": [s["name"] for s in self.species.values()],
                "kinetics": "gas",
                "reactions": "all",
                "standard-concentration-basis": "unity",
                "state": {"T": 423.15, "P": 101325}
            }],
            "species": [],
            "reactions": []
        }
        
        # 2. Format reactions as irreversible to avoid thermodynamic bias (P0 fix)
        cantera_reactions: List[Dict[str, Any]] = []
        for r in self.reactions:
            cantera_reac = dict(r)
            cantera_reac["reversible"] = False
            cantera_reactions.append(cantera_reac)
        data["reactions"] = cantera_reactions
        
        cantera_species: List[Dict[str, Any]] = []
        for s in self.species.values():
            try:
                # P3 Fix: get_nasa_coefficients now returns 14 floats (2 ranges)
                coeffs = get_nasa_coefficients(s["smiles"])
                thermo_block = {
                    "model": "NASA7",
                    "temperature-ranges": [300.0, 1000.0, 3000.0],
                    "data": [coeffs[:7], coeffs[7:]]
                }
            except Exception as e:
                # Fallback to constant-cp if estimation fails
                print(f"  Warning: Thermo estimation failed for {s['name']} ({s['smiles']}): {e}")
                thermo_block = {
                    "model": "constant-cp",
                    "t0": 298.15,
                    "h0": 0.0,
                    "s0": 0.0,
                    "cp0": 150000.0 # Higher for condensed phase
                }

            # R.1: Each species needs an equation-of-state block for ideal-condensed
            cantera_species.append({
                "name": s["name"],
                "composition": s["composition"],
                "thermo": thermo_block,
                "equation-of-state": {
                    "model": "constant-volume",
                    "molar-volume": s["molar_volume"],
                }
            })
        data["species"] = cantera_species
            
        with open(output_path, "w") as f:
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)
            
        print(f"Exported Cantera mechanism to {output_path}")

if __name__ == "__main__":
    exporter = CanteraExporter()
    exporter.add_reaction(["C(C1C(C(C(C(O1)O)O)O)O)=O"], ["Product_SMILES"], 25.0)
    exporter.export_yaml("test_mech.yaml")
