#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import urllib.request
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


ROOT = Path(__file__).resolve().parents[2]
XTB_INPUT_ROOT = ROOT / "data" / "geometries" / "xtb_inputs"
CACHE_ROOT = ROOT / "results" / "computational_gap_refinement" / "pubchem_geometry_refresh"


TARGET_SPECS = {
    "hexanal_radical_quench": {
        "reactant_fragments": [
            {
                "smiles": "CCCCCC=O",
                "label": "hexanal",
                "source": {"kind": "pubchem", "cid": 6184},
            },
            {
                "smiles": "[SH]",
                "label": "thiyl_radical",
                "source": {"kind": "manual"},
            },
        ],
        "product": {
            "smiles": "CCCCCC(O)[S]",
            "label": "thiohemiacetal_radical_adduct",
            "source": {
                "kind": "pubchem_proxy",
                "cid": 157932425,
                "component_index": 0,
                "remove_sulfur_hydrogen": True,
                "fallback_smiles": "CCCCCC(O)S",
                "fallback_to_target_embedding": True,
            },
            "note": "Uses PubChem 1-sulfanylhexan-1-ol as the closest stable 3D proxy and removes the sulfur hydrogen to recover the radical atom count.",
        },
    },
    "lysinoalanine_crosslink": {
        "reactant_fragments": [
            {
                "smiles": "C=C(C(=O)O)N",
                "label": "dehydroalanine",
                "source": {"kind": "pubchem", "cid": 123991},
            },
            {
                "smiles": "NCCCCC(N)C(=O)O",
                "label": "l_lysine",
                "source": {"kind": "pubchem", "cid": 5962},
            },
        ],
        "product": {
            "smiles": "C(CCNCC(C(=O)O)N)CC(C(=O)O)N",
            "label": "lysinoalanine",
            "source": {"kind": "pubchem", "cid": 29269},
        },
    },
    "aa_ring_open_dicarbonyl": {
        "reactant_fragments": [
            {
                "smiles": "C(C(C1C(=O)C(=O)C(=O)O1)O)O",
                "label": "dehydroascorbic_acid",
                "source": {"kind": "pubchem", "cid": 440667},
            },
            {
                "smiles": "O",
                "label": "water",
                "source": {"kind": "pubchem", "cid": 962},
            },
        ],
        "product": {
            "smiles": "C(C(C(C(=O)C(=O)C(=O)O)O)O)O",
            "label": "diketogulonic_acid",
            "source": {"kind": "pubchem", "cid": 18871},
        },
    },
}


def _build_hydrogenated_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"RDKit could not embed a 3D conformer for {smiles}")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:
        pass
    return mol


def _download_pubchem_sdf(cid: int, cache_dir: Path) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"cid_{cid}.sdf"
    if path.exists():
        return path
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
    with urllib.request.urlopen(url) as response:
        path.write_bytes(response.read())
    return path


def _load_pubchem_mol(source: dict, cache_dir: Path) -> tuple[Chem.Mol, Path]:
    sdf_path: Path | None = None
    mol: Chem.Mol | None = None
    try:
        sdf_path = _download_pubchem_sdf(int(source["cid"]), cache_dir)
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
        if mol is None:
            text = sdf_path.read_text(encoding="utf-8")
            mol = Chem.MolFromMolBlock(text, removeHs=False, strictParsing=False)
    except Exception:
        fallback_smiles = source.get("fallback_smiles")
        if fallback_smiles:
            mol = _build_hydrogenated_mol(str(fallback_smiles))
        else:
            raise
    if mol is None:
        raise RuntimeError(f"Could not obtain a geometry template for PubChem CID {source['cid']}")
    component_index = source.get("component_index")
    if component_index is not None:
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        mol = fragments[int(component_index)]
    if source.get("remove_sulfur_hydrogen"):
        mol = _remove_sulfur_hydrogen(mol)
    return mol, sdf_path


def _remove_sulfur_hydrogen(mol: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(mol)
    remove_idx = None
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 1:
                remove_idx = neighbor.GetIdx()
                break
        if remove_idx is not None:
            break
    if remove_idx is None:
        raise RuntimeError("Could not find a sulfur-bound hydrogen to remove from the proxy geometry")
    rw_mol.RemoveAtom(remove_idx)
    new_mol = rw_mol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[0].strip())
    atom_lines = [line for line in lines[2:] if line.strip()]
    if len(atom_lines) != atom_count:
        raise RuntimeError(f"Unexpected XYZ length in {path}")
    symbols: list[str] = []
    coords = np.zeros((atom_count, 3), dtype=float)
    for index, line in enumerate(atom_lines):
        parts = line.split()
        symbols.append(parts[0])
        coords[index] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return symbols, coords


def _write_xyz(path: Path, symbols: list[str], coords: np.ndarray, comment: str) -> None:
    lines = [str(len(symbols)), comment]
    for symbol, row in zip(symbols, coords):
        lines.append(f"{symbol:<2} {row[0]: .12f} {row[1]: .12f} {row[2]: .12f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _conformer_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    coords = np.zeros((mol.GetNumAtoms(), 3), dtype=float)
    for index in range(mol.GetNumAtoms()):
        position = conf.GetAtomPosition(index)
        coords[index] = [position.x, position.y, position.z]
    return coords


def _renumber_template_to_target(template: Chem.Mol, target: Chem.Mol) -> Chem.Mol:
    full_match = template.GetSubstructMatch(target, useChirality=False)
    if full_match and len(full_match) == target.GetNumAtoms():
        return Chem.RenumberAtoms(template, list(full_match))

    target_heavy = Chem.RemoveHs(target)
    template_heavy = Chem.RemoveHs(template)
    heavy_match = template_heavy.GetSubstructMatch(target_heavy, useChirality=False)
    if not heavy_match or len(heavy_match) != target_heavy.GetNumAtoms():
        raise RuntimeError("Could not map the public geometry onto the local target ordering")

    target_heavy_original = [atom.GetIdx() for atom in target.GetAtoms() if atom.GetAtomicNum() > 1]
    template_heavy_original = [atom.GetIdx() for atom in template.GetAtoms() if atom.GetAtomicNum() > 1]
    heavy_map = {
        target_heavy_original[target_index]: template_heavy_original[template_index]
        for target_index, template_index in enumerate(heavy_match)
    }

    template_hydrogens_by_parent: dict[int, list[int]] = {}
    for atom in template.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        parent = next((neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1), None)
        if parent is None:
            continue
        template_hydrogens_by_parent.setdefault(parent, []).append(atom.GetIdx())

    target_hydrogens_by_parent: dict[int, list[int]] = {}
    for atom in target.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        parent = next((neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1), None)
        if parent is None:
            continue
        target_hydrogens_by_parent.setdefault(parent, []).append(atom.GetIdx())

    full_order = [-1] * target.GetNumAtoms()
    for target_index, template_index in heavy_map.items():
        full_order[target_index] = template_index

    for target_parent, target_hydrogens in target_hydrogens_by_parent.items():
        template_parent = heavy_map[target_parent]
        template_hydrogens = template_hydrogens_by_parent.get(template_parent, [])
        if len(template_hydrogens) < len(target_hydrogens):
            raise RuntimeError("Public geometry does not contain enough hydrogens for the local target")
        for target_hydrogen, template_hydrogen in zip(target_hydrogens, template_hydrogens):
            full_order[target_hydrogen] = template_hydrogen

    if any(index < 0 for index in full_order):
        raise RuntimeError("Failed to assign a full atom ordering from the public geometry")
    return Chem.RenumberAtoms(template, full_order)


def _align_coords_to_reference(mobile: np.ndarray, reference: np.ndarray, heavy_mask: np.ndarray) -> np.ndarray:
    mobile_core = mobile[heavy_mask]
    reference_core = reference[heavy_mask]
    mobile_centroid = mobile_core.mean(axis=0)
    reference_centroid = reference_core.mean(axis=0)
    mobile_centered = mobile_core - mobile_centroid
    reference_centered = reference_core - reference_centroid
    covariance = mobile_centered.T @ reference_centered
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1, :] *= -1
        rotation = right_t.T @ left.T
    return (mobile - mobile_centroid) @ rotation + reference_centroid


def _symbols_from_mol(mol: Chem.Mol) -> list[str]:
    return [atom.GetSymbol() for atom in mol.GetAtoms()]


def _translate_to_reference_centroid(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return mobile - mobile.mean(axis=0) + reference.mean(axis=0)


def _backup_current_xyz(target_id: str, path: Path, cache_dir: Path) -> str:
    backup_dir = cache_dir / target_id
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup_path = backup_dir / f"{path.stem}.previous.xyz"
    backup_path.write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    return str(backup_path.relative_to(ROOT))


def _refresh_reactant(target_id: str, spec: dict, target_dir: Path, cache_dir: Path) -> list[dict]:
    reactant_path = target_dir / "reactant.xyz"
    current_symbols, current_coords = _read_xyz(reactant_path)
    refreshed_symbols: list[str] = []
    refreshed_coords_chunks: list[np.ndarray] = []
    metadata: list[dict] = []
    cursor = 0

    pending_specs = list(spec["reactant_fragments"])
    while pending_specs:
        matched_index = None
        fragment_spec = None
        fragment_target = None
        atom_count = 0
        fragment_symbols: list[str] = []
        align_mode = "atomwise"
        for index, candidate in enumerate(pending_specs):
            candidate_target = _build_hydrogenated_mol(candidate["smiles"])
            candidate_symbols = _symbols_from_mol(candidate_target)
            candidate_count = candidate_target.GetNumAtoms()
            if current_symbols[cursor:cursor + candidate_count] == candidate_symbols:
                matched_index = index
                fragment_spec = candidate
                fragment_target = candidate_target
                atom_count = candidate_count
                fragment_symbols = candidate_symbols
                align_mode = "atomwise"
                break
            current_slice = current_symbols[cursor:cursor + candidate_count]
            if len(current_slice) == candidate_count and sorted(current_slice) == sorted(candidate_symbols):
                matched_index = index
                fragment_spec = candidate
                fragment_target = candidate_target
                atom_count = candidate_count
                fragment_symbols = candidate_symbols
                align_mode = "centroid"
                break
        if fragment_spec is None or fragment_target is None or matched_index is None:
            raise RuntimeError(f"Could not match the current fragment ordering for {target_id} at atom index {cursor}")

        fragment_reference = current_coords[cursor:cursor + atom_count]
        pending_specs.pop(matched_index)

        source = fragment_spec["source"]
        if source["kind"] == "manual":
            refreshed_symbols.extend(fragment_symbols)
            refreshed_coords_chunks.append(fragment_reference.copy())
            metadata.append({
                "label": fragment_spec["label"],
                "kind": "manual",
                "note": "Kept the manually oriented local radical fragment because no exact public 3D radical record is available.",
            })
        else:
            template_mol, sdf_path = _load_pubchem_mol(source, cache_dir)
            ordered_template = _renumber_template_to_target(template_mol, fragment_target)
            mobile_coords = _conformer_coords(ordered_template)
            if align_mode == "atomwise":
                heavy_mask = np.array([atom.GetAtomicNum() > 1 for atom in ordered_template.GetAtoms()], dtype=bool)
                aligned_coords = _align_coords_to_reference(mobile_coords, fragment_reference, heavy_mask)
            else:
                aligned_coords = _translate_to_reference_centroid(mobile_coords, fragment_reference)
            refreshed_symbols.extend(_symbols_from_mol(ordered_template))
            refreshed_coords_chunks.append(aligned_coords)
            metadata.append({
                "label": fragment_spec["label"],
                "kind": source["kind"],
                "cid": int(source["cid"]),
                "alignment_mode": align_mode,
                "sdf_path": None if sdf_path is None else str(sdf_path.relative_to(ROOT)),
            })

        cursor += atom_count

    refreshed_coords = np.vstack(refreshed_coords_chunks)
    comment = f"{target_id} reactant refreshed from public 3D templates"
    _write_xyz(reactant_path, refreshed_symbols, refreshed_coords, comment)
    return metadata


def _refresh_product(target_id: str, spec: dict, target_dir: Path, cache_dir: Path) -> dict:
    product_path = target_dir / "product.xyz"
    current_symbols, current_coords = _read_xyz(product_path)
    product_target = _build_hydrogenated_mol(spec["product"]["smiles"])
    target_symbols = _symbols_from_mol(product_target)
    align_mode = "atomwise" if current_symbols == target_symbols else "centroid"

    template_mol, sdf_path = _load_pubchem_mol(spec["product"]["source"], cache_dir)
    used_target_embedding = False
    try:
        ordered_template = _renumber_template_to_target(template_mol, product_target)
    except RuntimeError:
        if spec["product"]["source"].get("fallback_to_target_embedding"):
            ordered_template = product_target
            used_target_embedding = True
        else:
            raise
    mobile_coords = _conformer_coords(ordered_template)
    if align_mode == "atomwise":
        heavy_mask = np.array([atom.GetAtomicNum() > 1 for atom in ordered_template.GetAtoms()], dtype=bool)
        aligned_coords = _align_coords_to_reference(mobile_coords, current_coords, heavy_mask)
    else:
        aligned_coords = _translate_to_reference_centroid(mobile_coords, current_coords)
    comment = f"{target_id} product refreshed from public 3D template"
    _write_xyz(product_path, target_symbols, aligned_coords, comment)
    payload = {
        "label": spec["product"]["label"],
        "kind": spec["product"]["source"]["kind"],
        "cid": int(spec["product"]["source"]["cid"]),
        "alignment_mode": align_mode,
        "sdf_path": None if sdf_path is None else str(sdf_path.relative_to(ROOT)),
    }
    if spec["product"].get("note"):
        payload["note"] = spec["product"]["note"]
    if used_target_embedding:
        payload["note"] = (
            (payload.get("note", "") + " ").strip()
            + "Fell back to a local ETKDG embedding of the radical target because the stable sulfur-hydrogen proxy could not be mapped atom-for-atom onto the radical graph."
        ).strip()
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target",
        action="append",
        default=[],
        help="Refresh only the named target(s). Repeat the flag to select multiple entries.",
    )
    args = parser.parse_args()

    selected = set(args.target) if args.target else set(TARGET_SPECS)
    unknown = sorted(target for target in selected if target not in TARGET_SPECS)
    if unknown:
        raise SystemExit(f"Unknown target(s): {', '.join(unknown)}")

    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    report: dict[str, dict] = {}
    for target_id in sorted(selected):
        spec = TARGET_SPECS[target_id]
        target_dir = XTB_INPUT_ROOT / target_id
        if not target_dir.exists():
            raise RuntimeError(f"Target directory does not exist: {target_dir}")

        reactant_path = target_dir / "reactant.xyz"
        product_path = target_dir / "product.xyz"
        backups = {
            "reactant_previous": _backup_current_xyz(target_id, reactant_path, CACHE_ROOT),
            "product_previous": _backup_current_xyz(target_id, product_path, CACHE_ROOT),
        }
        reactant_sources = _refresh_reactant(target_id, spec, target_dir, CACHE_ROOT)
        product_source = _refresh_product(target_id, spec, target_dir, CACHE_ROOT)
        report[target_id] = {
            "target_dir": str(target_dir.relative_to(ROOT)),
            "backups": backups,
            "reactant_sources": reactant_sources,
            "product_source": product_source,
        }
        print(f"[{target_id}] refreshed reactant.xyz and product.xyz from public 3D templates")

    report_path = CACHE_ROOT / "refresh_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {report_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())