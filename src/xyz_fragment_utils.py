"""Utilities for splitting multi-fragment XYZ files by bond connectivity.

Used by the DFT pipeline to optimise bimolecular reactants as individual
fragments, which is thermodynamically correct for computing activation
barriers:  ΔG‡ = G(TS) − Σᵢ G(fragmentᵢ).
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from src.xyz_common import parse_xyz as _parse_xyz

# Covalent radii (Å) — sufficient subset for organic + S/N/O/P chemistry.
# Source: Cordero et al. 2008, Dalton Trans.
_COVALENT_RADII: Dict[str, float] = {
    "H": 0.31, "He": 0.28,
    "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39,
    "Si": 1.11, "Se": 1.20,
}
_DEFAULT_RADIUS = 0.77  # generic fallback

_BOND_TOLERANCE = 1.3  # factor × (r_i + r_j) for bond detection


def _build_adjacency(
    symbols: List[str],
    coords: List[Tuple[float, float, float]],
) -> List[List[int]]:
    """Build an adjacency list from covalent-radius bond detection."""
    n = len(symbols)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        ri = _COVALENT_RADII.get(symbols[i], _DEFAULT_RADIUS)
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            rj = _COVALENT_RADII.get(symbols[j], _DEFAULT_RADIUS)
            xj, yj, zj = coords[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            dist_sq = dx * dx + dy * dy + dz * dz
            threshold = _BOND_TOLERANCE * (ri + rj)
            if dist_sq < threshold * threshold:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _connected_components(adj: List[List[int]]) -> List[List[int]]:
    """BFS to find connected components.  Returns list of sorted atom-index lists."""
    n = len(adj)
    visited = [False] * n
    components: List[List[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        queue = [start]
        visited[start] = True
        component: List[int] = []
        while queue:
            node = queue.pop(0)
            component.append(node)
            for nb in adj[node]:
                if not visited[nb]:
                    visited[nb] = True
                    queue.append(nb)
        components.append(sorted(component))
    return components


def split_xyz_into_fragments(xyz_str: str) -> List[str]:
    """Split a multi-fragment XYZ string into per-fragment XYZ strings.

    If the molecule is a single connected component, returns a list
    with one element (the original XYZ).

    Returns:
        List of XYZ strings, one per connected fragment.
    """
    symbols, coords = _parse_xyz(xyz_str)
    adj = _build_adjacency(symbols, coords)
    components = _connected_components(adj)

    if len(components) <= 1:
        return [xyz_str]

    lines = xyz_str.strip().splitlines()
    # Preserve the original comment line style
    comment = lines[1] if len(lines) > 1 else ""
    atom_lines = lines[2:2 + len(symbols)]

    fragments: List[str] = []
    for comp in components:
        n = len(comp)
        frag_lines = [str(n), comment]
        for idx in comp:
            frag_lines.append(atom_lines[idx])
        fragments.append("\n".join(frag_lines) + "\n")
    return fragments


def count_unpaired_electrons(symbols: List[str]) -> int:
    """Heuristic: count unpaired electrons for a neutral radical fragment.

    Rules (for the Maillard chemistry scope):
    - [SH] / [S] radical → 1 unpaired electron
    - [OH] radical → 1 unpaired electron
    - Organic molecules with even total electrons → 0 (singlet)
    - Otherwise → 1

    This is a best-effort heuristic.  For edge cases, the manifest
    ``fragment_spins`` field should be used instead.
    """
    # Atomic numbers
    _Z = {
        "H": 1, "He": 2, "C": 6, "N": 7, "O": 8, "F": 9,
        "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53, "Si": 14, "Se": 34,
    }
    total_electrons = sum(_Z.get(s, 0) for s in symbols)
    return total_electrons % 2  # 0 → singlet, 1 → doublet


def assign_fragment_charge_spin(
    fragments_xyz: List[str],
    total_charge: int,
    total_spin: int,
    manifest_fragment_spins: List[int] | None = None,
) -> List[Tuple[int, int]]:
    """Assign (charge, spin) to each fragment.

    Strategy:
    1. If ``manifest_fragment_spins`` is provided, use it directly.
    2. Otherwise, use ``count_unpaired_electrons`` per fragment to infer
       2S values.  Charge is distributed as: charge = 0 for all fragments
       (valid for neutral bimolecular reactions — which is the entire
       Maillard scope).

    Args:
        fragments_xyz: Per-fragment XYZ strings.
        total_charge: Overall system charge.
        total_spin: Overall system 2S (e.g. 1 = doublet).
        manifest_fragment_spins: Optional explicit per-fragment 2S values.

    Returns:
        List of (charge, spin) tuples, same length as fragments_xyz.
    """
    n = len(fragments_xyz)

    if manifest_fragment_spins is not None and len(manifest_fragment_spins) == n:
        # Charge: all-zero for now (Maillard scope is neutral reactions)
        return [(0, s) for s in manifest_fragment_spins]

    # Auto-detect from electron count parity
    result: List[Tuple[int, int]] = []
    remaining_spin = total_spin
    for frag_xyz in fragments_xyz:
        syms, _ = _parse_xyz(frag_xyz)
        unpaired = count_unpaired_electrons(syms)
        frag_spin = min(unpaired, remaining_spin)
        remaining_spin -= frag_spin
        result.append((0, frag_spin))

    # Sanity: total spin should match
    assigned_total = sum(s for _, s in result)
    if assigned_total != total_spin:
        # Fallback: put all remaining spin on the smallest fragment (radical)
        result_sorted = sorted(range(n), key=lambda i: len(_parse_xyz(fragments_xyz[i])[0]))
        result = [(0, 0)] * n
        result[result_sorted[0]] = (0, total_spin)

    return result
