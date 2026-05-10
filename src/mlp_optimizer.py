"""MLP-backed geometry optimization helpers and thin backend wrapper."""

from __future__ import annotations

import importlib
import io
import logging
import math
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from typing import Dict, List, Optional, Tuple

import numpy as np

from src.logger import get_logger
from src.xyz_common import parse_xyz

logger = get_logger(__name__)

try:
    from ase.io import read, write
    from ase.optimize import BFGS
    from mace.calculators import mace_mp
    _MLP_AVAILABLE = True
except ImportError:
    read = write = BFGS = mace_mp = None
    _MLP_AVAILABLE = False

def _parse_xyz_coordinates(xyz_string: str) -> Tuple[List[str], np.ndarray]:
    atoms, coords_tuples = parse_xyz(xyz_string)
    return atoms, np.asarray(coords_tuples, dtype=float)


def _angle_degrees(origin: np.ndarray, left: np.ndarray, right: np.ndarray) -> Optional[float]:
    left_vec = left - origin
    right_vec = right - origin
    left_norm = float(np.linalg.norm(left_vec))
    right_norm = float(np.linalg.norm(right_vec))
    if left_norm < 1e-12 or right_norm < 1e-12:
        return None
    cosine = float(np.dot(left_vec, right_vec) / (left_norm * right_norm))
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def compute_geometry_drift_metrics(reference_xyz: str, test_xyz: str) -> Dict[str, Optional[float]]:
    ref_atoms, ref_coords = _parse_xyz_coordinates(reference_xyz)
    test_atoms, test_coords = _parse_xyz_coordinates(test_xyz)
    if ref_atoms != test_atoms:
        raise ValueError("Geometry drift metrics require identical atom ordering and element sequence")

    displacements = np.linalg.norm(test_coords - ref_coords, axis=1)
    sulfur_indices = [index for index, symbol in enumerate(ref_atoms) if symbol == "S"]

    sulfur_local_rmsd: Optional[float] = None
    sulfur_neighbor_max_delta: Optional[float] = None
    sulfur_bond_max_delta: Optional[float] = None
    sulfur_angle_max_delta: Optional[float] = None

    for sulfur_index in sulfur_indices:
        neighbor_indices = [
            index
            for index in range(len(ref_atoms))
            if index != sulfur_index and float(np.linalg.norm(ref_coords[sulfur_index] - ref_coords[index])) <= 2.4
        ]
        if not neighbor_indices:
            continue

        local_indices = [sulfur_index, *neighbor_indices]
        local_rmsd = float(np.sqrt(np.mean(np.sum((test_coords[local_indices] - ref_coords[local_indices]) ** 2, axis=1))))
        sulfur_local_rmsd = local_rmsd if sulfur_local_rmsd is None else max(sulfur_local_rmsd, local_rmsd)

        distance_deltas = [
            abs(
                float(np.linalg.norm(ref_coords[sulfur_index] - ref_coords[index]))
                - float(np.linalg.norm(test_coords[sulfur_index] - test_coords[index]))
            )
            for index in neighbor_indices
        ]
        if distance_deltas:
            max_delta = max(distance_deltas)
            sulfur_neighbor_max_delta = max_delta if sulfur_neighbor_max_delta is None else max(sulfur_neighbor_max_delta, max_delta)
            sulfur_bond_max_delta = max_delta if sulfur_bond_max_delta is None else max(sulfur_bond_max_delta, max_delta)

        if len(neighbor_indices) >= 2:
            angle_deltas: List[float] = []
            for left_offset, left_index in enumerate(neighbor_indices):
                for right_index in neighbor_indices[left_offset + 1:]:
                    ref_angle = _angle_degrees(ref_coords[sulfur_index], ref_coords[left_index], ref_coords[right_index])
                    test_angle = _angle_degrees(test_coords[sulfur_index], test_coords[left_index], test_coords[right_index])
                    if ref_angle is None or test_angle is None:
                        continue
                    angle_deltas.append(abs(test_angle - ref_angle))
            if angle_deltas:
                max_angle_delta = max(angle_deltas)
                sulfur_angle_max_delta = max_angle_delta if sulfur_angle_max_delta is None else max(sulfur_angle_max_delta, max_angle_delta)

    return {
        "max_atom_displacement_angstrom": float(np.max(displacements)) if len(displacements) else None,
        "sulfur_local_rmsd_angstrom": sulfur_local_rmsd,
        "sulfur_neighbor_max_delta_angstrom": sulfur_neighbor_max_delta,
        "sulfur_bond_max_delta_angstrom": sulfur_bond_max_delta,
        "sulfur_angle_max_delta_degrees": sulfur_angle_max_delta,
    }


def _should_reject_bounded_preopt(
    metrics: Dict[str, Optional[float]],
    *,
    max_atom_drift_threshold: float = 1.0,
    sulfur_neighbor_delta_threshold: float = 0.25,
) -> bool:
    max_atom_drift = float(metrics.get("max_atom_displacement_angstrom") or 0.0)
    sulfur_neighbor_delta = float(metrics.get("sulfur_neighbor_max_delta_angstrom") or 0.0)
    sulfur_bond_delta = float(metrics.get("sulfur_bond_max_delta_angstrom") or 0.0)
    return (
        max_atom_drift > max_atom_drift_threshold
        or sulfur_neighbor_delta > sulfur_neighbor_delta_threshold
        or sulfur_bond_delta > sulfur_neighbor_delta_threshold
    )


def _resolve_backend_locator(model_family: str, model_name: str, backend_locator: Optional[str]) -> Optional[str]:
    family = str(model_family or "mace_mp").strip().lower()
    locator = str(backend_locator or "").strip()
    if family in {"mace_mp", "mace_off"}:
        if locator.startswith("builtin:"):
            built_in_name = locator.split(":", 1)[1].strip()
            return built_in_name or model_name
        return locator or None
    if family == "mace_omol":
        if not locator:
            raise ImportError("mace_omol requires an explicit backend locator")
        return locator
    return locator or None


def _read_single_atoms(xyz_string: str):
    assert read is not None
    with io.StringIO(xyz_string) as handle:
        atoms = read(handle, format="xyz")
    if isinstance(atoms, list):
        if not atoms:
            raise ValueError("XYZ did not contain any readable frames")
        return atoms[-1]
    return atoms


def _write_atoms(atoms) -> str:
    assert write is not None
    with io.StringIO() as handle:
        write(handle, atoms, format="xyz")
        return handle.getvalue()


class MLPOptimizer:
    """Wrapper for ASE-driven geometry optimization using MACE-compatible backends."""

    def __init__(
        self,
        model_name: str = "medium",
        device: str = "cpu",
        default_dtype: str = "float64",
        model_family: str = "mace_mp",
        backend_locator: Optional[str] = None,
    ):
        self.model_name = model_name
        self.device = device
        self.model_family = str(model_family or "mace_mp").strip().lower()
        self.backend_locator = backend_locator
        resolved_locator = _resolve_backend_locator(self.model_family, model_name, backend_locator)

        if self.model_family in {"mace_mp", "mace_off"}:
            if not _MLP_AVAILABLE:
                raise ImportError("MLPOptimizer dependencies (mace-torch, ase) are not installed in the current environment.")
            selected_model = resolved_locator or model_name
            assert mace_mp is not None
            with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
                self.calc = mace_mp(
                    model=selected_model,
                    dispersion=False,
                    default_dtype=default_dtype,
                    device=device,
                )
            return

        if self.model_family == "mace_omol":
            if not resolved_locator:
                raise ImportError("mace_omol requires an explicit backend locator")
            if not resolved_locator.startswith("python:"):
                raise ImportError(
                    "mace_omol backend locator is configured, but only python:module adapters are supported in-repo"
                )
            module_name = resolved_locator.split(":", 1)[1].strip()
            if not module_name:
                raise ImportError("mace_omol backend locator is invalid")
            plugin_module = importlib.import_module(module_name)
            if not hasattr(plugin_module, "build_calculator"):
                raise ImportError("Configured mace_omol backend does not expose build_calculator")
            self.calc = plugin_module.build_calculator(
                model_name=model_name,
                locator=backend_locator,
                device=device,
                default_dtype=default_dtype,
            )
            return

        raise ImportError(f"Unsupported MLP model family: {self.model_family}")

    def optimize_geometry(self, xyz_string: str, fmax: float = 0.01, max_steps: int = 500, drift_threshold: float = 1.0) -> str:
        atoms = _read_single_atoms(xyz_string)
        atoms.calc = self.calc

        assert BFGS is not None
        with io.StringIO() as buf, redirect_stdout(buf):
            opt = BFGS(atoms, logfile=None)
            opt.run(fmax=fmax, steps=max_steps)

        optimized_xyz = _write_atoms(atoms)
        metrics = compute_geometry_drift_metrics(xyz_string, optimized_xyz)
        if _should_reject_bounded_preopt(
            metrics,
            max_atom_drift_threshold=drift_threshold,
            sulfur_neighbor_delta_threshold=min(0.25, drift_threshold),
        ):
            logger.warning(
                ">>> [MLPOptimizer] WARNING: Excessive local drift detected "
                f"(max atom {float(metrics.get('max_atom_displacement_angstrom') or 0.0):.2f} Å, "
                f"sulfur neighbor {float(metrics.get('sulfur_neighbor_max_delta_angstrom') or 0.0):.2f} Å). "
                "MLP result may be chemically invalid. Reverting to original geometry."
            )
            text = xyz_string.strip()
            return f"{text}\n" if text else ""

        return optimized_xyz

    def count_imaginary_modes(self, xyz_string: str, delta: float = 0.01) -> int:
        atoms = _read_single_atoms(xyz_string)
        atoms.calc = self.calc

        from ase.vibrations import Vibrations

        with tempfile.TemporaryDirectory(prefix="mlp-vib-") as tmpdir:
            vib = Vibrations(atoms, name=f"{tmpdir}/mode", delta=delta)
            vib.run()
            try:
                frequencies = vib.get_frequencies()
            finally:
                vib.clean()

        count = 0
        for freq in frequencies:
            real_part = float(np.real(freq))
            imag_part = float(np.imag(freq))
            if abs(imag_part) > 1e-6 or real_part < -50.0:
                count += 1
        return count

    def find_ts_relaxed_scan(
        self,
        reactant_xyz: str,
        forming_bond_atoms: tuple[int, int],
        *,
        fmax: float = 0.05,
        max_steps: int = 100,
        start_distance: float | None = None,
        stop_distance: float = 1.60,
        step_size: float = 0.10,
    ) -> tuple[str, dict[str, float | int]]:
        seed_atoms = _read_single_atoms(reactant_xyz)

        left_index, right_index = forming_bond_atoms
        current_distance = float(seed_atoms.get_distance(left_index, right_index))
        scan_start = float(start_distance if start_distance is not None else max(current_distance, stop_distance + 0.5))
        if scan_start <= stop_distance:
            scan_start = stop_distance + step_size
        if scan_start <= stop_distance:
            raise ValueError(
                f"Relaxed scan cannot proceed: start ({scan_start:.2f} Å) <= stop ({stop_distance:.2f} Å). "
                f"Current distance is {current_distance:.2f} Å."
            )

        from ase.constraints import FixBondLengths

        distances = []
        cursor = scan_start
        while cursor >= stop_distance - 1e-8:
            distances.append(round(cursor, 4))
            cursor -= step_size
        if not distances:
            distances = [round(scan_start, 4)]

        best_atoms = None
        best_energy = None
        best_distance = None
        for target_distance in distances:
            atoms = seed_atoms.copy()
            atoms.calc = self.calc
            atoms.set_distance(left_index, right_index, target_distance, fix=0)
            atoms.set_constraint(FixBondLengths([(left_index, right_index)]))
            with io.StringIO() as buf, redirect_stdout(buf):
                opt = BFGS(atoms, logfile=None)
                opt.run(fmax=fmax, steps=max_steps)
            energy = float(atoms.get_potential_energy())
            if best_energy is None or energy > best_energy:
                best_atoms = atoms.copy()
                best_energy = energy
                best_distance = target_distance

        if best_atoms is None or best_energy is None or best_distance is None:
            raise RuntimeError("MLP relaxed scan failed to produce a TS seed")

        return _write_atoms(best_atoms), {
            "best_energy_ev": best_energy,
            "best_scan_distance_angstrom": best_distance,
            "scan_point_count": len(distances),
            "initial_distance_angstrom": current_distance,
        }

    def optimize_ts(self, xyz_string: str, fmax: float = 0.05, max_steps: int = 200) -> str:
        from .ts_optimizer import TSOptimizer

        ts_opt = TSOptimizer(fmax=fmax, max_steps=max_steps)
        atoms = _read_single_atoms(xyz_string)

        try:
            result = ts_opt.find_ts(atoms, self.calc)
        except (ImportError, RuntimeError) as exc:
            logging.getLogger(__name__).warning(
                "Sella TS search failed; falling back to geometry minimization: %s",
                exc,
            )
            return self.optimize_geometry(xyz_string, fmax=fmax, max_steps=max_steps)

        if isinstance(result, list):
            atoms = result[-1]
        else:
            atoms = result
        return _write_atoms(atoms)
