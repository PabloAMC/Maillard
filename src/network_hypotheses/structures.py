"""Species structures: data/species/structures.yml joined with the compound registry, as RDKit molecules."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from src import data_access, data_paths


def _rdkit():
    try:
        from rdkit import Chem  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without rdkit
        raise RuntimeError(
            "the hypothesis layer needs RDKit (conda-forge `rdkit`); it is optional for the engine "
            "and required here"
        ) from exc
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    return Chem


@dataclass(frozen=True)
class Structure:
    key: str                 # engine species key, or the canonical SMILES for a product not in the table
    kind: str                # "molecule" | "lump"
    smiles: Optional[str]
    canonical: Optional[str]  # canonical SMILES without stereochemistry: the identity used for matching
    inchikey: Optional[str]
    registry_id: Optional[str] = None
    note: str = ""


def canonical(smiles: str) -> Optional[str]:
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def inchikey(smiles: str) -> Optional[str]:
    Chem = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(mol) if mol is not None else None


def load(literature: bool = True) -> Dict[str, Structure]:
    """Every engine species, as a Structure; lumps have no SMILES. With ``literature`` (the default)
    the table also holds the species the rules name that the engine does not: literature
    intermediates and products (data/species/literature_structures.yml), keyed so they never collide
    with an engine key; the placement logic treats them as known compounds without an engine
    reaction, which is what "mechanism known, not modelled" means."""
    table = dict(data_access.load_yaml(data_paths.SPECIES_STRUCTURES)["species"])
    if literature:
        extra = data_access.load_yaml(data_paths.LITERATURE_STRUCTURES)["species"]
        clash = sorted(set(extra) & set(table))
        if clash:
            raise ValueError(f"literature_structures.yml re-uses engine species keys: {clash}")
        table.update(extra)
    registry = {c["id"]: c for c in data_access.load_yaml(data_paths.COMPOUND_REGISTRY)["compounds"]}
    out: Dict[str, Structure] = {}
    for key, entry in table.items():
        if entry["kind"] == "lump":
            out[key] = Structure(key, "lump", None, None, None, None, entry.get("note", ""))
            continue
        rid = entry.get("registry_id")
        smiles = entry.get("smiles") or registry[rid]["smiles"]
        out[key] = Structure(key, "molecule", smiles, canonical(smiles), inchikey(smiles), rid, entry.get("note", ""))
    return out


def registry_by_canonical() -> Dict[str, str]:
    """canonical SMILES -> registry compound id, for products that are known compounds but not engine species."""
    registry = data_access.load_yaml(data_paths.COMPOUND_REGISTRY)["compounds"]
    out = {}
    for c in registry:
        if c.get("smiles"):
            can = canonical(c["smiles"])
            if can:
                out.setdefault(can, c["id"])
    return out
