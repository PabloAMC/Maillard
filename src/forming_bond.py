from __future__ import annotations

import re
from math import dist
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from src.xyz_common import parse_xyz as _parse_xyz


ROOT = Path(__file__).resolve().parents[1]

_COVALENT_RADII_ANGSTROM = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "P": 1.07,
    "S": 1.05,
    "F": 0.57,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}


def _covalent_radius(symbol: str) -> float:
    return _COVALENT_RADII_ANGSTROM.get(symbol, 0.77)


def _bond_cutoff(symbol_left: str, symbol_right: str) -> float:
    scale = 1.25 if "H" in {symbol_left, symbol_right} else 1.20
    return scale * (_covalent_radius(symbol_left) + _covalent_radius(symbol_right))


def _infer_bonds(
    atoms: Sequence[str],
    coords: Sequence[Tuple[float, float, float]],
) -> List[Tuple[int, int]]:
    bonds: List[Tuple[int, int]] = []
    for left_index, left_symbol in enumerate(atoms):
        left_coord = coords[left_index]
        for right_index in range(left_index + 1, len(atoms)):
            right_symbol = atoms[right_index]
            cutoff = _bond_cutoff(left_symbol, right_symbol)
            if dist(left_coord, coords[right_index]) <= cutoff:
                bonds.append((left_index, right_index))
    return bonds


def _connected_components(atom_count: int, bonds: Sequence[Tuple[int, int]]) -> List[int]:
    adjacency = {index: set() for index in range(atom_count)}
    for left_index, right_index in bonds:
        adjacency[left_index].add(right_index)
        adjacency[right_index].add(left_index)

    component_ids = [-1] * atom_count
    component = 0
    for atom_index in range(atom_count):
        if component_ids[atom_index] >= 0:
            continue
        stack = [atom_index]
        component_ids[atom_index] = component
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if component_ids[neighbor] >= 0:
                    continue
                component_ids[neighbor] = component
                stack.append(neighbor)
        component += 1
    return component_ids


def _find_family_recipe_path(source_row: Mapping[str, Any]) -> Optional[str]:
    for basis_entry in source_row.get("surrogate_basis", []):
        basis_text = str(basis_entry).strip()
        match = re.search(r"data/rmg_extensions/families/([^/]+)/rules\.py", basis_text)
        if not match:
            continue
        family = match.group(1)
        return f"data/rmg_extensions/families/{family}/groups.py"
    return None


def _extract_form_bond_actions(recipe_relative_path: Optional[str]) -> List[Dict[str, Any]]:
    if not recipe_relative_path:
        return []
    recipe_path = ROOT / recipe_relative_path
    if not recipe_path.exists():
        return []
    content = recipe_path.read_text(encoding="utf-8")
    actions: List[Dict[str, Any]] = []
    for match in re.finditer(r"\['FORM_BOND',\s*'([^']+)',\s*([0-9]+),\s*'([^']+)'\]", content):
        actions.append(
            {
                "left_label": match.group(1),
                "order": int(match.group(2)),
                "right_label": match.group(3),
            }
        )
    return actions


def infer_forming_bond_metadata(
    source_row: Mapping[str, Any],
    *,
    reactant_relative_path: str,
    product_relative_path: str,
) -> Dict[str, Any]:
    recipe_relative_path = _find_family_recipe_path(source_row)
    recipe_actions = _extract_form_bond_actions(recipe_relative_path)
    reactant_path = ROOT / reactant_relative_path
    product_path = ROOT / product_relative_path

    if not reactant_path.exists() or not product_path.exists():
        return {
            "available": False,
            "method": "unavailable",
            "reason": "missing_geometry_pair",
            "family_recipe_source": recipe_relative_path,
            "family_recipe_form_bonds": recipe_actions,
        }

    reactant_atoms, reactant_coords = _parse_xyz(reactant_path.read_text(encoding="utf-8"))
    product_atoms, product_coords = _parse_xyz(product_path.read_text(encoding="utf-8"))
    if reactant_atoms != product_atoms:
        return {
            "available": False,
            "method": "unavailable",
            "reason": "element_sequence_mismatch",
            "family_recipe_source": recipe_relative_path,
            "family_recipe_form_bonds": recipe_actions,
        }

    reactant_bonds = set(_infer_bonds(reactant_atoms, reactant_coords))
    product_bonds = set(_infer_bonds(product_atoms, product_coords))
    reactant_components = _connected_components(len(reactant_atoms), sorted(reactant_bonds))
    new_bonds = sorted(product_bonds - reactant_bonds)

    candidates: List[Dict[str, Any]] = []
    for left_index, right_index in new_bonds:
        left_symbol = reactant_atoms[left_index]
        right_symbol = reactant_atoms[right_index]
        product_distance = dist(product_coords[left_index], product_coords[right_index])
        reactant_distance = dist(reactant_coords[left_index], reactant_coords[right_index])
        candidates.append(
            {
                "atom_indices_zero_based": [left_index, right_index],
                "atom_indices_one_based": [left_index + 1, right_index + 1],
                "atom_symbols": [left_symbol, right_symbol],
                "product_distance_angstrom": round(product_distance, 4),
                "reactant_distance_angstrom": round(reactant_distance, 4),
                "distance_drop_angstrom": round(reactant_distance - product_distance, 4),
                "interfragment_in_reactant": reactant_components[left_index] != reactant_components[right_index],
                "heavy_atom_bond": left_symbol != "H" and right_symbol != "H",
            }
        )

    selected = None
    if candidates:
        heavy_candidates = [candidate for candidate in candidates if candidate["heavy_atom_bond"]]
        interfragment_heavy = [candidate for candidate in heavy_candidates if candidate["interfragment_in_reactant"]]
        ranked = interfragment_heavy or heavy_candidates or candidates
        selected = max(
            ranked,
            key=lambda candidate: (
                int(candidate["interfragment_in_reactant"]),
                int(candidate["heavy_atom_bond"]),
                float(candidate["distance_drop_angstrom"]),
                -float(candidate["product_distance_angstrom"]),
            ),
        )

    return {
        "available": selected is not None,
        "method": "geometry_topology_diff",
        "selection_rule": "interfragment_heavy_bond" if selected and selected["interfragment_in_reactant"] else "largest_distance_drop",
        "family_recipe_source": recipe_relative_path,
        "family_recipe_form_bonds": recipe_actions,
        "all_new_bonds": candidates,
        "atom_indices_zero_based": list(selected["atom_indices_zero_based"]) if selected else [],
        "atom_indices_one_based": list(selected["atom_indices_one_based"]) if selected else [],
        "atom_symbols": list(selected["atom_symbols"]) if selected else [],
        "product_distance_angstrom": selected["product_distance_angstrom"] if selected else None,
        "reactant_distance_angstrom": selected["reactant_distance_angstrom"] if selected else None,
        "distance_drop_angstrom": selected["distance_drop_angstrom"] if selected else None,
        "interfragment_in_reactant": bool(selected["interfragment_in_reactant"]) if selected else False,
    }