from __future__ import annotations

import re
from typing import List, Tuple


def is_xyz_coordinate_line(line: str) -> bool:
    parts = line.strip().split()
    if len(parts) != 4:
        return False
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", parts[0]):
        return False
    try:
        float(parts[1])
        float(parts[2])
        float(parts[3])
    except ValueError:
        return False
    return True


def parse_xyz(xyz_content: str) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Parse an XYZ string into (symbols, coords).

    Tolerant of empty lines and a comment line after the atom count.
    Uses :func:`is_xyz_coordinate_line` to skip non-coordinate lines
    (e.g. comment, blank).

    Parameters
    ----------
    xyz_content : str
        Standard XYZ format: first line = atom count, optional comment,
        then ``atom_count`` coordinate lines.

    Returns
    -------
    Tuple[List[str], List[Tuple[float, float, float]]]
        (element_symbols, coordinates) where each coordinate is an
        (x, y, z) tuple in Ångström.

    Raises
    ------
    ValueError
        If the content is empty or atom count doesn't match coordinate
        lines found.
    """
    lines = [line.rstrip() for line in xyz_content.splitlines() if line.strip()]
    if not lines:
        raise ValueError("XYZ content is empty")

    atom_count = int(lines[0].strip())
    atom_lines = [line for line in lines[1:] if is_xyz_coordinate_line(line)][:atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(
            f"XYZ atom count mismatch: expected {atom_count} coordinates, "
            f"found {len(atom_lines)}"
        )

    symbols: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    for line in atom_lines:
        parts = line.split()[:4]
        symbols.append(parts[0])
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return symbols, coords


def extract_xyz_last_frame(xyz_content: str) -> str:
    lines = xyz_content.strip().splitlines()
    if not lines:
        return xyz_content
    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        return xyz_content
    frame_size = n_atoms + 2
    if frame_size < 3 or len(lines) < frame_size:
        return xyz_content
    n_frames = len(lines) // frame_size
    if n_frames <= 1:
        return xyz_content
    for frame_index in range(n_frames - 1, -1, -1):
        start = frame_index * frame_size
        candidate = lines[start:start + frame_size]
        if len(candidate) != frame_size:
            continue
        try:
            candidate_n = int(candidate[0].strip())
        except ValueError:
            continue
        if candidate_n == n_atoms:
            return "\n".join(candidate) + "\n"


def format_xyz(
    symbols: List[str],
    coords: List[Tuple[float, float, float]],
    comment: str = "",
) -> str:
    """Format atom symbols and coordinates into an XYZ string.

    Parameters
    ----------
    symbols : List[str]
        Element symbols, one per atom.
    coords : List[Tuple[float, float, float]]
        (x, y, z) coordinates in Ångström, one per atom.
    comment : str
        Optional comment for the second line.

    Returns
    -------
    str
        Standard XYZ-format string (atom-count, comment, coordinates).
    """
    lines = [str(len(symbols)), comment]
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"{sym:<2} {x: .12f} {y: .12f} {z: .12f}")
    return "\n".join(lines) + "\n"
    return xyz_content