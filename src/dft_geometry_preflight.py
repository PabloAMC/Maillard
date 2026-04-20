from __future__ import annotations

from math import dist, sqrt
from typing import Any, Dict, List, Tuple

from src.xyz_common import parse_xyz as _parse_xyz, format_xyz as _format_xyz


DEFAULT_MIN_INTERATOMIC_DISTANCE_ANGSTROM = 0.80
DEFAULT_MIN_PAIRWISE_RMS_DELTA_ANGSTROM = 0.30  # reject TS guess ≈ reactant
DEFAULT_MAX_PAIRWISE_RMS_DELTA_ANGSTROM = 1.75
DEFAULT_MAX_PAIRWISE_DELTA_ANGSTROM = 6.00


def _min_pairwise_distance(coords: List[Tuple[float, float, float]]) -> float | None:
    if len(coords) < 2:
        return None
    minimum = float("inf")
    for index, left in enumerate(coords):
        for right in coords[index + 1 :]:
            minimum = min(minimum, dist(left, right))
    return minimum


def _pairwise_distance_delta(
    reactant_coords: List[Tuple[float, float, float]],
    ts_coords: List[Tuple[float, float, float]],
) -> Tuple[float | None, float | None]:
    if len(reactant_coords) != len(ts_coords) or len(reactant_coords) < 2:
        return None, None

    squared_error = 0.0
    sample_count = 0
    max_delta = 0.0
    for index, reactant_left in enumerate(reactant_coords):
        ts_left = ts_coords[index]
        for offset, reactant_right in enumerate(reactant_coords[index + 1 :], start=index + 1):
            reactant_distance = dist(reactant_left, reactant_right)
            ts_distance = dist(ts_left, ts_coords[offset])
            delta = abs(reactant_distance - ts_distance)
            squared_error += delta * delta
            sample_count += 1
            max_delta = max(max_delta, delta)

    if sample_count == 0:
        return None, None
    return sqrt(squared_error / sample_count), max_delta


def repair_steric_clash(
    xyz_content: str,
    *,
    min_interatomic_distance_angstrom: float = DEFAULT_MIN_INTERATOMIC_DISTANCE_ANGSTROM,
    target_distance_angstrom: float = 1.00,
    max_iterations: int = 6,
) -> Dict[str, Any]:
    atoms, coords = _parse_xyz(xyz_content)
    repaired_coords = [list(coord) for coord in coords]
    repairs: List[Dict[str, Any]] = []

    for _ in range(max_iterations):
        minimum_distance = _min_pairwise_distance([tuple(coord) for coord in repaired_coords])
        if minimum_distance is None or minimum_distance >= min_interatomic_distance_angstrom:
            break

        worst_pair = None
        worst_distance = float("inf")
        for left_index, left_coord in enumerate(repaired_coords):
            for right_index in range(left_index + 1, len(repaired_coords)):
                current_distance = dist(tuple(left_coord), tuple(repaired_coords[right_index]))
                if current_distance < worst_distance:
                    worst_distance = current_distance
                    worst_pair = (left_index, right_index)

        if worst_pair is None:
            break

        left_index, right_index = worst_pair
        left_coord = repaired_coords[left_index]
        right_coord = repaired_coords[right_index]
        axis = [right_coord[dimension] - left_coord[dimension] for dimension in range(3)]
        axis_norm = sqrt(sum(component * component for component in axis))
        if axis_norm < 1e-8:
            axis = [1.0, 0.0, 0.0]
            axis_norm = 1.0

        unit_axis = [component / axis_norm for component in axis]
        displacement = max(0.0, target_distance_angstrom - worst_distance) / 2.0
        for dimension in range(3):
            left_coord[dimension] -= unit_axis[dimension] * displacement
            right_coord[dimension] += unit_axis[dimension] * displacement

        repairs.append(
            {
                "atom_indices_zero_based": [left_index, right_index],
                "atom_indices_one_based": [left_index + 1, right_index + 1],
                "atom_symbols": [atoms[left_index], atoms[right_index]],
                "distance_before_angstrom": round(worst_distance, 4),
                "target_distance_angstrom": target_distance_angstrom,
            }
        )

    final_coords = [tuple(coord) for coord in repaired_coords]
    final_min = _min_pairwise_distance(final_coords)
    repair_complete = final_min is None or final_min >= min_interatomic_distance_angstrom
    return {
        "xyz": _format_xyz(atoms, final_coords),
        "repairs": repairs,
        "min_interatomic_distance_angstrom": final_min,
        "repaired": bool(repairs),
        "repair_complete": repair_complete,
    }


def build_dft_geometry_preflight(
    reactant_xyz: str,
    ts_xyz: str,
    *,
    min_interatomic_distance_angstrom: float = DEFAULT_MIN_INTERATOMIC_DISTANCE_ANGSTROM,
    min_pairwise_rms_delta_angstrom: float = DEFAULT_MIN_PAIRWISE_RMS_DELTA_ANGSTROM,
    max_pairwise_rms_delta_angstrom: float = DEFAULT_MAX_PAIRWISE_RMS_DELTA_ANGSTROM,
    max_pairwise_delta_angstrom: float = DEFAULT_MAX_PAIRWISE_DELTA_ANGSTROM,
) -> Dict[str, Any]:
    reactant_atoms, reactant_coords = _parse_xyz(reactant_xyz)
    ts_atoms, ts_coords = _parse_xyz(ts_xyz)

    reactant_min_distance = _min_pairwise_distance(reactant_coords)
    ts_min_distance = _min_pairwise_distance(ts_coords)
    pairwise_rms_delta, pairwise_max_delta = _pairwise_distance_delta(reactant_coords, ts_coords)

    blockers: List[str] = []
    if len(reactant_atoms) != len(ts_atoms):
        blockers.append("atom_count_mismatch")
    if reactant_atoms != ts_atoms:
        blockers.append("element_sequence_mismatch")
    if reactant_min_distance is not None and reactant_min_distance < min_interatomic_distance_angstrom:
        blockers.append("reactant_steric_clash")
    if ts_min_distance is not None and ts_min_distance < min_interatomic_distance_angstrom:
        blockers.append("ts_steric_clash")
    if len(reactant_coords) >= 3 and pairwise_rms_delta is not None and pairwise_rms_delta < min_pairwise_rms_delta_angstrom:
        blockers.append("ts_too_similar_to_reactant")
    if pairwise_rms_delta is not None and pairwise_rms_delta > max_pairwise_rms_delta_angstrom:
        blockers.append("ts_pairwise_distortion_rms")
    if pairwise_max_delta is not None and pairwise_max_delta > max_pairwise_delta_angstrom:
        blockers.append("ts_pairwise_distortion_max")

    return {
        "reactant_atom_count": len(reactant_atoms),
        "ts_atom_count": len(ts_atoms),
        "atom_count_match": len(reactant_atoms) == len(ts_atoms),
        "element_sequence_match": reactant_atoms == ts_atoms,
        "reactant_min_interatomic_distance_angstrom": reactant_min_distance,
        "ts_min_interatomic_distance_angstrom": ts_min_distance,
        "pairwise_distance_rms_delta_angstrom": pairwise_rms_delta,
        "pairwise_distance_max_delta_angstrom": pairwise_max_delta,
        "thresholds": {
            "min_interatomic_distance_angstrom": min_interatomic_distance_angstrom,
            "min_pairwise_rms_delta_angstrom": min_pairwise_rms_delta_angstrom,
            "max_pairwise_rms_delta_angstrom": max_pairwise_rms_delta_angstrom,
            "max_pairwise_delta_angstrom": max_pairwise_delta_angstrom,
        },
        "blockers": blockers,
        "quality_gate_passed": not blockers,
    }


def preflight_candidate_score(report: Dict[str, Any]) -> Tuple[float, float]:
    rms_delta = float(report.get("pairwise_distance_rms_delta_angstrom") or 0.0)
    max_delta = float(report.get("pairwise_distance_max_delta_angstrom") or 0.0)
    return rms_delta, -max_delta