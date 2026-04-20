#!/usr/bin/env python3
"""
Regenerate reactant/product XYZ files for computational-gap targets
with correct atom-type alignment (required by xTB path search).

Root cause of previous failures:
1. RDKit embeds multi-fragment SMILES with overlapping fragments (steric clash)
2. Independent embedding of reactant and product produces different atom orderings
   (e.g., N at index 0 in reactant, O at index 0 in product), causing xTB to abort
   with "Atom type mismatch between reactant and product"

This script fixes both problems by:
1. Embedding fragments independently with multi-seed search, then combining with spacing
2. Using atom-mapped SMILES and RDKit's RenumberAtoms to ensure atom ordering matches

After running this script, re-run xTB path search:
    ./scripts/docker_maillard.sh computational-gap-xtb <target>
"""

from __future__ import annotations

import argparse
import shutil
import sys
from math import dist
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("ERROR: RDKit is required. Install via: conda install -c conda-forge rdkit")
    sys.exit(1)


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_BASE = ROOT / "data" / "geometries" / "xtb_inputs"
MIN_INTERATOMIC_DISTANCE = 0.60  # Angstroms, matches dft_geometry_preflight.py
FRAGMENT_SEPARATION_TARGET = 2.5  # Angstroms, minimum inter-fragment spacing


# ──────────────────────────────────────────────────────────────
# Target definitions with atom-mapped SMILES
# ──────────────────────────────────────────────────────────────
# Atom maps (:N) ensure that atom N in the reactant corresponds to
# atom N in the product. The product will be reordered to match
# the reactant's atom sequence.

TARGETS_TO_FIX: Dict[str, dict] = {
    "lysinoalanine_crosslink": {
        # Dehydroalanine + Lysine epsilon-amine -> Lysinoalanine crosslink
        # Michael addition: Lys-NH2 attacks DHA β-carbon (C=CH2)
        "reactant_smiles": "[CH2:1]=[C:2]([NH2:3])[C:4](=[O:5])[OH:6].[NH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH:12]([NH2:13])[C:14](=[O:15])[OH:16]",
        "product_smiles": "[CH2:1]([NH:7][CH2:8][CH2:9][CH2:10][CH2:11][CH:12]([NH2:13])[C:14](=[O:15])[OH:16])[CH:2]([NH2:3])[C:4](=[O:5])[OH:6]",
        "reactive_bonds": [(1, 7)],
        "note": "DHA + Lys Michael addition (35 atoms). Two reactant fragments.",
    },
    "pe_schiff_base": {
        # Ethanolamine + open-chain aldose -> Schiff base + water
        # Condensation: amine + aldehyde -> imine + H2O
        "reactant_smiles": "[NH2:1][CH2:2][CH2:3][OH:4].[CH:5](=[O:6])[CH:7]([OH:8])[CH:9]([OH:10])[CH:11]([OH:12])[CH2:13][OH:14]",
        "product_smiles": "[CH2:2]([N:1]=[CH:5][CH:7]([OH:8])[CH:9]([OH:10])[CH:11]([OH:12])[CH2:13][OH:14])[CH2:3][OH:4].[OH2:6]",
        "note": "Ethanolamine + aldose Schiff base + water (31 atoms). Two fragments each side.",
    },
    "quinone_cys_michael": {
        # 1,4-Benzoquinone + methanethiol -> 1,4-Michael addition product
        # Thiol S adds to C3, H migrates to C4 (C3-C4 now single bond)
        "reactant_smiles": "[O:1]=[C:2]1[CH:3]=[CH:4][C:5](=[O:6])[CH:7]=[CH:8]1.[CH3:9][SH:10]",
        "product_smiles": "[O:1]=[C:2]1[CH:3]([S:10][CH3:9])[CH2:4][C:5](=[O:6])[CH:7]=[CH:8]1",
        "note": "Benzoquinone + methanethiol 1,4-addition (18 atoms). Two reactant fragments.",
    },
}


# ──────────────────────────────────────────────────────────────
# Geometry helpers
# ──────────────────────────────────────────────────────────────

def _parse_xyz_coords(xyz_content: str) -> List[Tuple[float, float, float]]:
    """Parse coordinates from XYZ content."""
    lines = [line.rstrip() for line in xyz_content.splitlines() if line.strip()]
    atom_count = int(lines[0].strip())
    remaining = lines[1:]
    parts0 = remaining[0].strip().split()
    if len(parts0) != 4:
        remaining = remaining[1:]
    else:
        try:
            float(parts0[1]); float(parts0[2]); float(parts0[3])
        except ValueError:
            remaining = remaining[1:]
    coords = []
    for line in remaining[:atom_count]:
        parts = line.split()[:4]
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return coords


def _min_pairwise_distance(coords: List[Tuple[float, float, float]]) -> float:
    """Compute minimum interatomic distance."""
    minimum = float("inf")
    for i, left in enumerate(coords):
        for right in coords[i + 1:]:
            minimum = min(minimum, dist(left, right))
    return minimum


def _mol_min_dist(mol: Chem.Mol) -> float:
    """Compute minimum interatomic distance from an RDKit mol with conformer."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    min_d = float("inf")
    for i in range(n):
        pi = conf.GetAtomPosition(i)
        for j in range(i + 1, n):
            pj = conf.GetAtomPosition(j)
            min_d = min(min_d, dist((pi.x, pi.y, pi.z), (pj.x, pj.y, pj.z)))
    return min_d


# ──────────────────────────────────────────────────────────────
# Embedding with fragment separation and multi-seed search
# ──────────────────────────────────────────────────────────────

def _embed_single_fragment(mol: Chem.Mol, n_seeds: int = 50) -> Chem.Mol:
    """Embed a single-fragment mol, trying multiple seeds for best geometry."""
    best_mol = None
    best_min_dist = -1.0
    for seed in range(n_seeds):
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        candidate = Chem.RWMol(mol)
        if AllChem.EmbedMolecule(candidate, params) != 0:
            continue
        AllChem.MMFFOptimizeMolecule(candidate, maxIters=1000)
        md = _mol_min_dist(candidate)
        if md > best_min_dist:
            best_min_dist = md
            best_mol = candidate
        if best_min_dist >= MIN_INTERATOMIC_DISTANCE:
            break
    if best_mol is None:
        raise RuntimeError(f"EmbedMolecule failed for all {n_seeds} seeds")
    return best_mol


def _embed_with_fragment_separation(mol: Chem.Mol, reactive_bonds: List[Tuple[int, int]] = None) -> Chem.Mol:
    """Embed mol, handling multi-fragment SMILES with independent optimization + spacing."""
    mol = Chem.AddHs(mol)
    frag_indices = Chem.GetMolFrags(mol)

    if reactive_bonds:
        # Pre-Docking: We use an MMFF Distance Constraint to pull reacting atoms together to ~3.0 Å
        # This properly docks intermolecular complexes without breaking explicit valences.
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("EmbedMolecule failed for distance-constrained pre-docking")
            
        conf = mol.GetConformer()
        
        # Apply initial rigid axial spacing to avoid overlapping random embeds
        frag_indices = Chem.GetMolFrags(mol)
        x_offset = 0.0
        for frag_atom_indices in frag_indices:
            positions = [conf.GetAtomPosition(i) for i in frag_atom_indices]
            cx = sum(p.x for p in positions) / len(positions)
            for gi in frag_atom_indices:
                pos = conf.GetAtomPosition(gi)
                conf.SetAtomPosition(gi, (pos.x - cx + x_offset, pos.y, pos.z))
            x_offset += 6.0
            
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
        
        for map_a, map_b in reactive_bonds:
            idx_a, idx_b = -1, -1
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == map_a:
                    idx_a = atom.GetIdx()
                elif atom.GetAtomMapNum() == map_b:
                    idx_b = atom.GetIdx()
            if idx_a != -1 and idx_b != -1:
                ff.MMFFAddDistanceConstraint(idx_a, idx_b, False, 2.7, 3.3, 1000.0)
                
        ff.Initialize()
        ff.Minimize(maxIts=5000)
        
        return mol

    # Fallback to simple X-axis separation for non-reactive separated targets
    frag_indices = Chem.GetMolFrags(mol)
    if len(frag_indices) <= 1:
        return _embed_single_fragment(mol)

    frag_mols = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    optimized_frags = []
    for frag in frag_mols:
        frag = Chem.AddHs(frag)
        optimized_frags.append(_embed_single_fragment(frag))

    params = AllChem.ETKDGv3()
    for attempt in range(10):
        params.randomSeed = 42 + attempt
        if AllChem.EmbedMolecule(mol, params) == 0:
            break
    else:
        raise RuntimeError("EmbedMolecule failed for combined mol")

    conf = mol.GetConformer()
    x_offset = 0.0
    for frag_atom_indices, frag_mol in zip(frag_indices, optimized_frags):
        frag_conf = frag_mol.GetConformer()
        positions = [frag_conf.GetAtomPosition(i) for i in range(frag_mol.GetNumAtoms())]
        cx = sum(p.x for p in positions) / len(positions)
        cy = sum(p.y for p in positions) / len(positions)
        cz = sum(p.z for p in positions) / len(positions)

        for local_idx, global_idx in enumerate(frag_atom_indices):
            pos = frag_conf.GetAtomPosition(local_idx)
            conf.SetAtomPosition(global_idx, (pos.x - cx + x_offset, pos.y - cy, pos.z - cz))

        placed_positions = [conf.GetAtomPosition(gi) for gi in frag_atom_indices]
        max_x = max(p.x for p in placed_positions)
        x_offset = max_x + FRAGMENT_SEPARATION_TARGET

    return mol


# ──────────────────────────────────────────────────────────────
# Atom-map-based reordering
# ──────────────────────────────────────────────────────────────

def _strip_maps(mol: Chem.Mol) -> Chem.Mol:
    """Remove atom map numbers from a mol (in place)."""
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def _get_map_to_idx(mol: Chem.Mol) -> Dict[int, int]:
    """Return {atom_map_number: atom_index} for mapped atoms."""
    return {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }


def _reorder_product_to_match_reactant(
    reactant_smiles: str,
    product_smiles: str,
) -> Tuple[Chem.Mol, Chem.Mol]:
    """
    Parse atom-mapped SMILES for reactant and product. Return (reactant_mol, product_mol)
    where product_mol has been renumbered so that atom i in product has the same element
    as atom i in reactant.

    The atom maps in SMILES define the correspondence:
      - Mapped atoms are placed in the same position
      - Unmapped hydrogens are appended in order after their parent heavy atom
    """
    r_mol = Chem.MolFromSmiles(reactant_smiles)
    p_mol = Chem.MolFromSmiles(product_smiles)
    if r_mol is None:
        raise ValueError(f"Invalid reactant SMILES: {reactant_smiles}")
    if p_mol is None:
        raise ValueError(f"Invalid product SMILES: {product_smiles}")

    # Add explicit H before renumbering
    r_mol = Chem.AddHs(r_mol)
    p_mol = Chem.AddHs(p_mol)

    # Build map: atom_map_number -> reactant_idx
    r_map = _get_map_to_idx(r_mol)
    p_map = _get_map_to_idx(p_mol)

    # All mapped atoms should appear in both
    mapped_keys = sorted(set(r_map.keys()) & set(p_map.keys()))
    if not mapped_keys:
        raise ValueError("No atom maps found in SMILES — cannot align")

    # Build the reordering for the product:
    # For each mapped heavy atom, place it at the reactant's position.
    # Then, for each heavy atom, place its implicit H neighbors right after.
    n_r = r_mol.GetNumAtoms()
    n_p = p_mol.GetNumAtoms()

    if n_r != n_p:
        raise ValueError(
            f"Atom count mismatch after AddHs: reactant={n_r}, product={n_p}. "
            f"SMILES are not balanced."
        )

    # Step 1: Build heavy-atom ordering by following reactant's atom order
    # For each atom in reactant order, find its map number, look up the product atom
    used_p_indices = set()
    p_new_order = [None] * n_p  # p_new_order[new_idx] = old_p_idx

    # First pass: place mapped (heavy) atoms
    pos = 0
    r_atom_order = []  # Track reactant atom order: (r_idx, map_num_or_None)
    for r_idx in range(n_r):
        r_atom = r_mol.GetAtomWithIdx(r_idx)
        map_num = r_atom.GetAtomMapNum()
        if map_num > 0 and map_num in p_map:
            p_idx = p_map[map_num]
            r_atom_order.append((r_idx, map_num, p_idx))

    # Build the product reorder: follow reactant order, placing mapped atoms
    # and their H neighbors together
    new_order = []
    visited_p = set()

    for r_idx in range(n_r):
        r_atom = r_mol.GetAtomWithIdx(r_idx)
        map_num = r_atom.GetAtomMapNum()

        if map_num > 0 and map_num in p_map:
            p_idx = p_map[map_num]
            if p_idx not in visited_p:
                new_order.append(p_idx)
                visited_p.add(p_idx)
        elif r_atom.GetAtomicNum() == 1:
            # This is an H in reactant. Find its parent's map number.
            r_parent = None
            for nbr in r_atom.GetNeighbors():
                if nbr.GetAtomMapNum() > 0:
                    r_parent = nbr.GetAtomMapNum()
                    break
            if r_parent is not None and r_parent in p_map:
                p_parent_idx = p_map[r_parent]
                p_parent_atom = p_mol.GetAtomWithIdx(p_parent_idx)
                # Find an unvisited H neighbor of this parent in product
                for nbr in p_parent_atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 1 and nbr.GetIdx() not in visited_p:
                        new_order.append(nbr.GetIdx())
                        visited_p.add(nbr.GetIdx())
                        break
                else:
                    # No H neighbor found on this parent; find any unvisited H
                    pass

    # Append any remaining product atoms not yet placed (safety net)
    for p_idx in range(n_p):
        if p_idx not in visited_p:
            new_order.append(p_idx)
            visited_p.add(p_idx)

    if len(new_order) != n_p:
        raise RuntimeError(f"Reorder bug: expected {n_p} atoms, got {len(new_order)}")

    # Apply renumbering
    p_mol = Chem.RenumberAtoms(p_mol, new_order)

    # Strip maps before returning
    _strip_maps(r_mol)
    _strip_maps(p_mol)

    return r_mol, p_mol


# ──────────────────────────────────────────────────────────────
# XYZ serialization
# ──────────────────────────────────────────────────────────────

def _mol_to_xyz_string(mol: Chem.Mol) -> str:
    """Convert RDKit mol with conformer to XYZ format string."""
    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    lines = [str(n_atoms), ""]
    for i in range(n_atoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}")
    return "\n".join(lines) + "\n"


# ──────────────────────────────────────────────────────────────
# Per-target fix
# ──────────────────────────────────────────────────────────────

def fix_target(target_id: str, *, dry_run: bool = False) -> dict:
    """Fix a single target's reactant/product geometries with atom alignment."""
    if target_id not in TARGETS_TO_FIX:
        return {"target_id": target_id, "status": "skipped", "reason": "not in fix list"}

    target = TARGETS_TO_FIX[target_id]
    target_dir = GEOMETRY_BASE / target_id
    reactant_path = target_dir / "reactant.xyz"
    product_path = target_dir / "product.xyz"

    print(f"\n=== {target_id} ===")
    print(f"  Note: {target['note']}")

    # Check current state
    if reactant_path.exists():
        old_coords = _parse_xyz_coords(reactant_path.read_text(encoding="utf-8"))
        old_min = _min_pairwise_distance(old_coords)
        print(f"  Current reactant min distance: {old_min:.4f} Å (threshold: {MIN_INTERATOMIC_DISTANCE})")
    else:
        old_min = None
        print("  No current reactant.xyz found")

    # Step 1: Parse and align atom ordering via mapped SMILES
    print(f"  Aligning atoms via mapped SMILES...")
    r_mol, p_mol = _reorder_product_to_match_reactant(
        target["reactant_smiles"],
        target["product_smiles"],
    )

    # Verify element-sequence match
    r_elements = [r_mol.GetAtomWithIdx(i).GetSymbol() for i in range(r_mol.GetNumAtoms())]
    p_elements = [p_mol.GetAtomWithIdx(i).GetSymbol() for i in range(p_mol.GetNumAtoms())]
    element_match = r_elements == p_elements
    if not element_match:
        mismatches = [(i, r_elements[i], p_elements[i]) for i in range(min(len(r_elements), len(p_elements))) if r_elements[i] != p_elements[i]]
        print(f"  ERROR: Element sequence mismatch after alignment: {mismatches[:5]}...")
        return {"target_id": target_id, "status": "error_alignment_failed",
                "detail": f"mismatches: {mismatches[:10]}"}

    print(f"  Element sequence match: ✓ ({len(r_elements)} atoms)")

    reactive_bonds = target.get("reactive_bonds")

    # Step 2: Embed 3D coordinates
    print(f"  Embedding reactant ({len(r_elements)} atoms)...")
    r_mol_3d = _embed_with_fragment_separation(r_mol, reactive_bonds=reactive_bonds)
    r_xyz = _mol_to_xyz_string(r_mol_3d)

    print(f"  Embedding product ({len(p_elements)} atoms)...")
    # For the product, the bond is usually already formed in the SMILES so we don't strictly need reactive_bonds,
    # but passing it does no harm since the bond will just not be double-added.
    p_mol_3d = _embed_with_fragment_separation(p_mol)
    p_xyz = _mol_to_xyz_string(p_mol_3d)

    # Step 3: Validate
    new_r_coords = _parse_xyz_coords(r_xyz)
    new_p_coords = _parse_xyz_coords(p_xyz)
    new_r_min = _min_pairwise_distance(new_r_coords)
    new_p_min = _min_pairwise_distance(new_p_coords)
    r_count = len(new_r_coords)
    p_count = len(new_p_coords)

    # Re-verify element sequence from the actual XYZ
    r_xyz_elements = [line.split()[0] for line in r_xyz.strip().splitlines()[2:2 + r_count]]
    p_xyz_elements = [line.split()[0] for line in p_xyz.strip().splitlines()[2:2 + p_count]]
    xyz_element_match = r_xyz_elements == p_xyz_elements

    print(f"  New reactant: {r_count} atoms, min distance: {new_r_min:.4f} Å")
    print(f"  New product:  {p_count} atoms, min distance: {new_p_min:.4f} Å")
    print(f"  XYZ element sequence match: {'✓' if xyz_element_match else '✗ FAIL'}")

    result = {
        "target_id": target_id,
        "reactant_atoms": r_count,
        "product_atoms": p_count,
        "atom_count_match": r_count == p_count,
        "element_sequence_match": xyz_element_match,
        "reactant_min_distance": new_r_min,
        "product_min_distance": new_p_min,
        "reactant_passes_gate": new_r_min >= MIN_INTERATOMIC_DISTANCE,
        "product_passes_gate": new_p_min >= MIN_INTERATOMIC_DISTANCE,
        "old_reactant_min_distance": old_min,
    }

    if r_count != p_count:
        print(f"  ERROR: Atom count mismatch ({r_count} vs {p_count})!")
        result["status"] = "error_atom_count_mismatch"
        return result

    if not xyz_element_match:
        print(f"  ERROR: Element sequence mismatch in final XYZ!")
        if r_count == p_count:
            mismatches = [(i, r_xyz_elements[i], p_xyz_elements[i]) for i in range(r_count) if r_xyz_elements[i] != p_xyz_elements[i]]
            print(f"    Mismatches: {mismatches[:5]}")
        result["status"] = "error_element_mismatch"
        return result

    if new_r_min < MIN_INTERATOMIC_DISTANCE:
        print(f"  WARNING: Reactant steric clash ({new_r_min:.4f} < {MIN_INTERATOMIC_DISTANCE})")
        result["status"] = "warning_reactant_clash"
    elif new_p_min < MIN_INTERATOMIC_DISTANCE:
        print(f"  WARNING: Product steric clash ({new_p_min:.4f} < {MIN_INTERATOMIC_DISTANCE})")
        result["status"] = "warning_product_clash"
    else:
        result["status"] = "ok"

    if dry_run:
        print("  DRY RUN: not writing files")
        result["written"] = False
        return result

    # Backup and clean
    target_dir.mkdir(parents=True, exist_ok=True)
    for path in [reactant_path, product_path]:
        if path.exists():
            backup = path.with_suffix(".xyz.bak")
            shutil.copy2(path, backup)
            print(f"  Backed up {path.name} -> {backup.name}")

    for pattern in ["xtbpath.xyz", "xtbpath_ts.xyz", "xtbpath_*.xyz"]:
        for stale_file in target_dir.glob(pattern):
            stale_file.unlink()
            print(f"  Removed stale: {stale_file.name}")

    reactant_path.write_text(r_xyz, encoding="utf-8")
    product_path.write_text(p_xyz, encoding="utf-8")
    print(f"  Wrote {reactant_path.name} ({r_count} atoms)")
    print(f"  Wrote {product_path.name} ({p_count} atoms)")
    result["written"] = True

    return result


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fix steric clashes and atom alignment in computational-gap geometry seeds."
    )
    parser.add_argument("--target", default="all", help="Target ID or 'all'")
    parser.add_argument("--dry-run", action="store_true", help="Validate without writing files")
    args = parser.parse_args()

    targets = list(TARGETS_TO_FIX.keys()) if args.target == "all" else [args.target]

    print("=" * 60)
    print("Computational-Gap Geometry Steric Clash & Alignment Fixer")
    print("=" * 60)

    results = []
    for target_id in targets:
        result = fix_target(target_id, dry_run=args.dry_run)
        results.append(result)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for r in results:
        status = r.get("status", "unknown")
        status_icon = "✓" if status == "ok" else "✗" if "error" in status else "⚠"
        written = " [written]" if r.get("written") else " [dry-run]" if r.get("status") != "skipped" else ""
        print(f"  {status_icon} {r['target_id']}: {status}{written}")
        if r.get("atom_count_match") is False:
            print(f"    Atom count: {r.get('reactant_atoms')} vs {r.get('product_atoms')} MISMATCH")
        if r.get("element_sequence_match") is False:
            print(f"    Element sequence: MISMATCH")

    has_errors = any("error" in r.get("status", "") for r in results)
    if has_errors:
        print("\nSome targets have errors. Fix the SMILES definitions and re-run.")
        return 1

    if not args.dry_run:
        print("\nNext steps:")
        for r in results:
            if r.get("written"):
                print(f"  ./scripts/docker_maillard.sh computational-gap-xtb {r['target_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
