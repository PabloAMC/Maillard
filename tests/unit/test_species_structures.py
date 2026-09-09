"""data/species/structures.yml: one structure (or a declared lump) per engine species, consistent with the engine.

2026-09-08. The hypothesis layer enumerates on these structures, so a wrong SMILES would propose wrong
chemistry; the engine's own atom counts and molar masses are the check. RDKit is the registry's tool
(data/keys/compounds.yml is keyed by InChIKey computed from SMILES), so it is required here too.
"""
from __future__ import annotations

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402

from src import data_access, data_paths  # noqa: E402
from src.kinetic_core import species, species_acrylamide, species_lipid, species_sulfur  # noqa: E402


def _engine_species():
    out = {}
    for mod in (species, species_sulfur, species_acrylamide, species_lipid):
        for name in dir(mod):
            obj = getattr(mod, name)
            if isinstance(obj, tuple) and obj and all(isinstance(s, species.Species) for s in obj):
                for s in obj:
                    out.setdefault(s.key, s)
    return out


def _engine_masses():
    out = {}
    for mod in (species_sulfur, species_acrylamide, species_lipid):
        for name in dir(mod):
            if "WEIGHT" in name or "MASS" in name or "MW" in name:
                obj = getattr(mod, name)
                if isinstance(obj, dict):
                    out.update({k: v for k, v in obj.items() if isinstance(v, (int, float))})
    return out


@pytest.fixture(scope="module")
def table():
    return data_access.load_yaml(data_paths.SPECIES_STRUCTURES)["species"]


@pytest.fixture(scope="module")
def registry():
    return {c["id"]: c for c in data_access.load_yaml(data_paths.COMPOUND_REGISTRY)["compounds"]}


def _smiles(entry, registry):
    if "smiles" in entry:
        return entry["smiles"]
    return registry[entry["registry_id"]]["smiles"]


def test_every_engine_species_has_exactly_one_entry(table):
    engine = _engine_species()
    assert set(table) == set(engine), (sorted(set(engine) - set(table)), sorted(set(table) - set(engine)))


def test_entries_are_molecules_or_declared_lumps(table, registry):
    for key, entry in table.items():
        assert entry["kind"] in ("molecule", "lump"), key
        if entry["kind"] == "lump":
            assert entry.get("note"), f"{key}: a lump must say why it has no structure"
            assert "smiles" not in entry and "registry_id" not in entry, key
        else:
            assert ("smiles" in entry) != ("registry_id" in entry), f"{key}: exactly one of smiles / registry_id"
            if "registry_id" in entry:
                assert entry["registry_id"] in registry, f"{key}: unknown registry id {entry['registry_id']}"
                assert registry[entry["registry_id"]].get("smiles"), f"{key}: registry entry has no SMILES"


def test_molecule_atom_counts_match_the_engine(table, registry):
    engine = _engine_species()
    bad = []
    for key, entry in table.items():
        if entry["kind"] != "molecule":
            continue
        mol = Chem.MolFromSmiles(_smiles(entry, registry))
        assert mol is not None, f"{key}: RDKit cannot parse the SMILES"
        counts = tuple(sum(1 for a in mol.GetAtoms() if a.GetSymbol() == el) for el in ("C", "N", "S"))
        want = (engine[key].carbon, engine[key].nitrogen, getattr(engine[key], "sulfur", 0))
        if counts != want:
            bad.append((key, counts, want))
    assert not bad, bad


#: Species whose declared molar mass in the engine disagrees with the structure (recorded 2026-09-08,
#: not silently corrected here, since the engine converts mmol/L to ug/L with these numbers):
#:   SBA        Asn + Glc - 2 H2O in the engine; the Schiff base is Asn + Glc - H2O (294.26)
#:   THI        the engine reports thiamine as its hydrochloride (337.27); the cation is 265.36
#:   DPO, TDP   130.10 in the engine; C5H8O4 is 132.11
#:   LOOH_*     310.47 (the hydroxide, C19H34O3); the hydroperoxide C19H34O4 is 326.48
#: Recorded in the backlog (tasks/data_restructure_plan.md section 7); each correction moves a reported
#: concentration, so it is a change of its own with the guards re-pinned.
KNOWN_MASS_DISAGREEMENTS = {"SBA", "THI", "DPO", "TDP", "LOOH_13_ct", "LOOH_13_tt", "LOOH_9_ct", "LOOH_9_tt"}


def test_molecule_masses_match_the_engine_where_declared(table, registry):
    masses = _engine_masses()
    bad = []
    for key, entry in table.items():
        if entry["kind"] != "molecule" or key not in masses or key in KNOWN_MASS_DISAGREEMENTS:
            continue
        mol = Chem.MolFromSmiles(_smiles(entry, registry))
        if abs(Descriptors.MolWt(mol) - masses[key]) > 0.15:
            bad.append((key, round(Descriptors.MolWt(mol), 2), masses[key]))
    assert not bad, bad


def test_known_disagreements_are_still_disagreements(table, registry):
    """When the engine's masses are corrected, remove the key from KNOWN_MASS_DISAGREEMENTS."""
    masses = _engine_masses()
    for key in KNOWN_MASS_DISAGREEMENTS:
        mol = Chem.MolFromSmiles(_smiles(table[key], registry))
        assert abs(Descriptors.MolWt(mol) - masses[key]) > 0.15, f"{key}: the engine's mass now agrees; drop it from the list"
