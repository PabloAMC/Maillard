"""
src/pathway_extractor.py — Extract elementary steps from RMG-Py output.

Parses RMG's Chemkin-format reaction outputs and adjacency lists (species dictionaries)
to generate a list of elementary reaction steps suitable for Tier 1 xTB screening.
"""

from dataclasses import dataclass
from typing import List, Optional
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import Descriptors
from functools import cached_property

@dataclass
class Species:
    label: str
    smiles: str

    @cached_property
    def is_volatile(self) -> bool:
        """
        Scientific heuristic for aroma volatility.
        Volatiles are generally small (MW < 160) and lack excessive polarity 
        (H-bond donors <= 1).
        """
        if not self.smiles:
            return False
            
        # Inerts and precursors that aren't aroma compounds
        inerts = {"water", "h2o", "h2", "co2", "ammonia", "h2s"}
        if self.label.lower() in inerts or self.smiles.lower() in inerts:
            return False

        try:
            mol = Chem.MolFromSmiles(self.smiles)
            if not mol:
                return False
            
            mw = Descriptors.MolWt(mol)
            h_donors = Descriptors.NumHDonors(mol)
            
            # Thresholds aligned with empirical flavor chemistry
            return mw < 160 and h_donors <= 1
        except Exception:
            return False
    
@dataclass
class ElementaryStep:
    reactants: List[Species]
    products: List[Species]
    reaction_family: Optional[str] = None
    rate_constant_k: Optional[float] = None
    source_quality: str = "heuristic" # "literature", "estimated_tst", "heuristic"
    barrier_uncertainty_kcal: float = 5.0 # Default heuristic uncertainty
    
    def __str__(self) -> str:
        reacts = " + ".join([r.label for r in self.reactants])
        prods = " + ".join([p.label for p in self.products])
        return f"{reacts} -> {prods} [{self.reaction_family or 'unknown'}]"

class PathwayExtractor:
    """Parses RMG chemkin outputs into explicit reaction objects."""
    
    def __init__(self, rmg_output_dir: Path):
        self.rmg_output_dir = Path(rmg_output_dir)
        self.chemkin_dir = self.rmg_output_dir / "chemkin"
        self.species_dict: dict[str, Species] = {}
        self.elementary_steps: List[ElementaryStep] = []
        
    def _parse_species_dictionary(self):
        """Parse species_dictionary.txt into Species objects with SMILES."""
        dict_path = self.chemkin_dir / "species_dictionary.txt"
        if not dict_path.exists():
            raise FileNotFoundError(f"Missing {dict_path}")
            
        current_label = None
        current_smiles = None
        
        # Simplified parser for tests: Looks for label and INCHI/SMILES annotations
        with open(dict_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    if current_label and current_smiles:
                        self.species_dict[current_label] = Species(current_label, current_smiles)
                    current_label = None
                    current_smiles = None
                    continue
                    
                # Match label (starts with text, not a number/whitespace) 
                if not current_label and not line.startswith("1") and "multiplicity" not in line:
                    # Very simple heuristic: if it doesn't have spaces and isn't a geometry tag
                    if " " not in line:
                        current_label = line
                        
                # Match SMILES from RMG metadata comments
                # E.g. // SMILES="CC=O"
                if current_label and "SMILES=" in line:
                    match = re.search(r'SMILES="([^"]+)"', line)
                    if match:
                        current_smiles = match.group(1)
                        
        # Flush last entry
        if current_label and current_smiles:
            self.species_dict[current_label] = Species(current_label, current_smiles)

    def _parse_chem_inp(self):
        """Parse chem_annotated.inp to extract reaction definitions."""
        inp_path = self.chemkin_dir / "chem_annotated.inp"
        if not inp_path.exists():
            raise FileNotFoundError(f"Missing {inp_path}")
            
        in_reactions_block = False
        
        with open(inp_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("REACTIONS"):
                    in_reactions_block = True
                    continue
                if line.startswith("END") and in_reactions_block:
                    break
                    
                if in_reactions_block and line and not line.startswith("!"):
                    # Typical line: A + B = C + D  1.0E3  0.0  10.0
                    rxn_part = line.split()[0] # get just the equation part
                    if "=" not in rxn_part:
                        continue
                        
                    left, right = rxn_part.split("=")
                    # RMG output often writes <B=> for reversible, so we split on "=" and clean
                    left = left.replace("<", "").replace(">", "").strip()
                    right = right.replace("<", "").replace(">", "").strip()
                    
                    reactants = [self.species_dict.get(s.strip()) for s in left.split("+") if s.strip()]
                    products = [self.species_dict.get(s.strip()) for s in right.split("+") if s.strip()]
                    
                    # Filter out any unresolved species
                    reactants = [r for r in reactants if r]
                    products = [p for p in products if p]
                    
                    if reactants and products:
                        self.elementary_steps.append(
                            ElementaryStep(reactants=reactants, products=products)
                        )

    def run(self) -> List[ElementaryStep]:
        """Execute parsing and return list of steps."""
        self._parse_species_dictionary()
        self._parse_chem_inp()
        return self.elementary_steps
