from dataclasses import dataclass, field
from typing import List, Optional
from functools import cached_property, lru_cache

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    Chem = None
    Descriptors = None


@dataclass
class Species:
    label: str
    smiles: str

    @cached_property
    def is_volatile(self) -> bool:
        """Scientific heuristic for aroma volatility.

        Volatiles are generally small (MW < 160) and lack excessive polarity
        (H-bond donors <= 1).
        """
        if not self.smiles:
            return False

        inerts = {"water", "h2o", "h2", "co2", "ammonia", "h2s"}
        if self.label.lower() in inerts or self.smiles.lower() in inerts:
            return False

        if Chem is None or Descriptors is None:
            return False

        try:
            mol = Chem.MolFromSmiles(self.smiles)
            if not mol:
                return False

            mw = Descriptors.MolWt(mol)
            h_donors = Descriptors.NumHDonors(mol)
            return mw < 160 and h_donors <= 1
        except Exception:
            return False


@dataclass
class ElementaryStep:
    reactants: List[Species]
    products: List[Species]
    reaction_family: Optional[str] = None
    rate_constant_k: Optional[float] = None
    source_quality: str = "heuristic"
    barrier_uncertainty_kcal: float = 5.0
    barrier_kcal_mol: Optional[float] = None
    reversible: bool = False
    direction: str = "forward"
    unresolved_species: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        reacts = " + ".join([r.label for r in self.reactants])
        prods = " + ".join([p.label for p in self.products])
        arrow = "<=>" if self.reversible else "->"
        suffix = f" [UNRESOLVED: {', '.join(self.unresolved_species)}]" if self.unresolved_species else ""
        return f"{reacts} {arrow} {prods} [{self.reaction_family or 'unknown'}]{suffix}"


@lru_cache(maxsize=4096)
def canonicalize_smiles(smi: str, fallback_to_original: bool = False, strip_salts: bool = False) -> Optional[str]:
    """Generates a canonical SMILES string using RDKit, stripping stereochemistry

    if necessary, and optionally isolating the largest fragment (strip_salts).
    Returns None if the SMILES is entirely invalid or RDKit is not installed,
    unless fallback_to_original is True.
    """
    if Chem is None or not smi:
        return smi if fallback_to_original else None

    try:
        m = Chem.MolFromSmiles(smi)
        if m:
            can_smi = Chem.MolToSmiles(m)
            if strip_salts and "." in can_smi:
                return max(set(can_smi.split(".")), key=len)
            return can_smi
    except Exception:
        pass
    return smi if fallback_to_original else None
