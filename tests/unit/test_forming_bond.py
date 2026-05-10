from pathlib import Path

import src.forming_bond as forming_bond
from src.dft_geometry_preflight import repair_steric_clash


def test_infer_forming_bond_metadata_from_geometry_diff(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(forming_bond, "ROOT", tmp_path)
    target_dir = tmp_path / "data/geometries/xtb_inputs/lysinoalanine_crosslink"
    target_dir.mkdir(parents=True)
    reactant = target_dir / "reactant.xyz"
    product = target_dir / "product.xyz"
    reactant.write_text(
        "4\nreactant\nC 0.0 0.0 0.0\nN 0.0 0.0 3.0\nH 0.0 0.0 4.0\nH 0.0 0.0 1.0\n",
        encoding="utf-8",
    )
    product.write_text(
        "4\nproduct\nC 0.0 0.0 0.0\nN 0.0 0.0 1.4\nH 0.0 0.0 2.4\nH 0.0 0.0 1.1\n",
        encoding="utf-8",
    )

    payload = forming_bond.infer_forming_bond_metadata(
        {
            "surrogate_basis": ["data/rmg_extensions/families/DHA_Crosslinking/rules.py"],
        },
        reactant_relative_path="data/geometries/xtb_inputs/lysinoalanine_crosslink/reactant.xyz",
        product_relative_path="data/geometries/xtb_inputs/lysinoalanine_crosslink/product.xyz",
    )

    assert payload["available"] is True
    assert payload["atom_indices_zero_based"] == [0, 1]
    assert payload["atom_symbols"] == ["C", "N"]
    assert payload["family_recipe_source"] == "data/rmg_extensions/families/DHA_Crosslinking/groups.py"


def test_repair_steric_clash_pushes_atoms_apart():
    payload = repair_steric_clash("2\nclash\nC 0.0 0.0 0.0\nN 0.0 0.0 0.4\n")

    assert payload["repaired"] is True
    assert payload["min_interatomic_distance_angstrom"] is not None
    assert payload["min_interatomic_distance_angstrom"] >= 0.8