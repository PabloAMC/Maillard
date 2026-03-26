"""
src/mlp_optimizer.py

Phase 10: MLP-Accelerated Geometry Optimization via MACE (mace-mp-0).

Replaces costly DFT structural optimizations with rapid Machine Learning Potential
optimizations. By default, relies on the `medium` MACE model via mace_mp, acting
as an ASE Calculator.
"""

import io
import logging
import os
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from src.logger import get_logger

logger = get_logger(__name__)

try:
    from ase.io import read, write
    from ase.optimize import BFGS
    import mace.calculators as mace_calculators
    _MLP_AVAILABLE = True
except ImportError:
    read = write = BFGS = mace_calculators = None
    _MLP_AVAILABLE = False


def _resolve_backend_locator(model_family: str, model_name: str, backend_locator: str | None) -> str:
    if not backend_locator:
        if model_family == "mace_omol":
            raise ImportError(
                "MACE-OMol requires an explicit backend locator; configure backend_locator in the candidate registry"
            )
        return model_name

    if backend_locator.startswith("builtin:"):
        locator = backend_locator.split(":", 1)[1].strip()
        if not locator:
            raise ImportError(f"Invalid builtin backend locator for {model_family}: {backend_locator}")
        return locator

    if backend_locator.startswith("env:"):
        env_var = backend_locator.split(":", 1)[1].strip()
        resolved = os.environ.get(env_var)
        if not resolved:
            raise ImportError(f"Backend locator environment variable is not set: {env_var}")
        return resolved

    if backend_locator.startswith("path:"):
        locator = backend_locator.split(":", 1)[1].strip()
    else:
        locator = backend_locator

    path = Path(locator).expanduser()
    if not path.exists():
        raise ImportError(f"Backend locator path does not exist for {model_family}: {path}")
    return str(path)


def _read_xyz_atoms(xyz_string: str):
    assert read is not None
    with io.StringIO(xyz_string.strip()) as handle:
        atoms = read(handle, format="xyz")
    if isinstance(atoms, list):
        return atoms[-1]
    return atoms


def compute_geometry_drift_metrics(
    reference_xyz: str,
    test_xyz: str,
    *,
    reactive_radius_angstrom: float = 2.3,
) -> dict[str, float | None]:
    try:
        from ase.build import minimize_rotation_and_translation
    except ImportError as exc:
        raise RuntimeError("ASE is required to compute geometry drift metrics") from exc

    ref_atoms = _read_xyz_atoms(reference_xyz)
    test_atoms = _read_xyz_atoms(test_xyz)
    if len(ref_atoms) != len(test_atoms):
        raise ValueError("Geometry drift metrics require atom-count preservation")

    ref_atoms = ref_atoms.copy()
    test_atoms = test_atoms.copy()
    minimize_rotation_and_translation(ref_atoms, test_atoms)

    delta = ref_atoms.get_positions() - test_atoms.get_positions()
    squared = np.sum(delta**2, axis=1)
    displacements = np.sqrt(squared)
    symbols = list(ref_atoms.get_chemical_symbols())

    hetero_indices = [index for index, symbol in enumerate(symbols) if symbol in {"N", "O", "S"}]
    sulfur_indices = [index for index, symbol in enumerate(symbols) if symbol == "S"]

    hetero_rmsd = None
    if hetero_indices:
        hetero_delta = delta[hetero_indices]
        hetero_rmsd = float(np.sqrt(np.mean(np.sum(hetero_delta**2, axis=1))))

    sulfur_local_rmsd = None
    sulfur_neighbor_delta = None
    sulfur_angle_max_delta = None
    sulfur_bond_max_delta = None
    if sulfur_indices:
        sulfur_index = sulfur_indices[0]
        ref_positions = ref_atoms.get_positions()
        ref_distances = np.linalg.norm(ref_positions - ref_positions[sulfur_index], axis=1)
        local_indices = [index for index, dist in enumerate(ref_distances) if dist <= reactive_radius_angstrom or index == sulfur_index]
        local_delta = delta[local_indices]
        sulfur_local_rmsd = float(np.sqrt(np.mean(np.sum(local_delta**2, axis=1))))
        test_positions = test_atoms.get_positions()
        ref_local_distances = np.linalg.norm(ref_positions[local_indices] - ref_positions[sulfur_index], axis=1)
        test_local_distances = np.linalg.norm(test_positions[local_indices] - test_positions[sulfur_index], axis=1)
        sulfur_neighbor_delta = float(np.max(np.abs(ref_local_distances - test_local_distances)))
        neighbor_indices = [index for index in local_indices if index != sulfur_index]
        if neighbor_indices:
            sulfur_bond_max_delta = float(
                np.max(
                    np.abs(
                        np.linalg.norm(ref_positions[neighbor_indices] - ref_positions[sulfur_index], axis=1)
                        - np.linalg.norm(test_positions[neighbor_indices] - test_positions[sulfur_index], axis=1)
                    )
                )
            )
        if len(neighbor_indices) >= 2:
            def _angle_degrees(positions: np.ndarray, center: int, left: int, right: int) -> float:
                vec_a = positions[left] - positions[center]
                vec_b = positions[right] - positions[center]
                cosine = float(np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b)))
                cosine = max(-1.0, min(1.0, cosine))
                return float(np.degrees(np.arccos(cosine)))

            angle_deltas = []
            for offset, left in enumerate(neighbor_indices[:-1]):
                for right in neighbor_indices[offset + 1:]:
                    ref_angle = _angle_degrees(ref_positions, sulfur_index, left, right)
                    test_angle = _angle_degrees(test_positions, sulfur_index, left, right)
                    angle_deltas.append(abs(ref_angle - test_angle))
            if angle_deltas:
                sulfur_angle_max_delta = float(max(angle_deltas))

    return {
        "rmsd_angstrom": float(np.sqrt(np.mean(squared))),
        "max_atom_displacement_angstrom": float(np.max(displacements)),
        "hetero_atom_rmsd_angstrom": hetero_rmsd,
        "sulfur_local_rmsd_angstrom": sulfur_local_rmsd,
        "sulfur_neighbor_max_delta_angstrom": sulfur_neighbor_delta,
        "sulfur_bond_max_delta_angstrom": sulfur_bond_max_delta,
        "sulfur_angle_max_delta_degrees": sulfur_angle_max_delta,
    }


def _should_reject_bounded_preopt(
    metrics: dict[str, float | None],
    *,
    max_atom_drift_threshold: float,
    sulfur_neighbor_delta_threshold: float,
) -> bool:
    if float(metrics.get("max_atom_displacement_angstrom") or 0.0) > max_atom_drift_threshold:
        return True
    sulfur_neighbor_delta = metrics.get("sulfur_neighbor_max_delta_angstrom")
    if sulfur_neighbor_delta is not None and float(sulfur_neighbor_delta) > sulfur_neighbor_delta_threshold:
        return True
    return False


def _build_mace_calculator(model_family: str, model_name: str, device: str, default_dtype: str, backend_locator: str | None):
    if not _MLP_AVAILABLE:
        raise ImportError("MLPOptimizer dependencies (mace-torch, ase) are not installed in the current environment.")

    assert mace_calculators is not None
    resolved_locator = _resolve_backend_locator(model_family, model_name, backend_locator)
    if model_family == "mace_mp":
        if not hasattr(mace_calculators, "mace_mp"):
            raise ImportError("mace_mp calculator is not available in the installed MACE package")
        return mace_calculators.mace_mp(  # type: ignore[attr-defined]
            model=resolved_locator,
            dispersion=False,
            default_dtype=default_dtype,
            device=device,
        )

    if model_family == "mace_omol":
        if hasattr(mace_calculators, "MACECalculator"):
            return mace_calculators.MACECalculator(  # type: ignore[attr-defined]
                model_paths=resolved_locator,
                device=device,
                default_dtype=default_dtype,
            )
        raise ImportError(
            "MACE-OMol backend requires MACECalculator support in the installed MACE package"
        )

    raise ImportError(f"Unsupported MACE model family: {model_family}")


class MLPOptimizer:
    """
    Wrapper for ASE-driven geometric optimization using the MACE neural network potential.
    """

    def __init__(
        self,
        model_name: str = "medium",
        device: str = "cpu",
        default_dtype: str = "float64",
        model_family: str = "mace_mp",
        backend_locator: str | None = None,
    ):
        """
        Initialize the MACE ASE Calculator.
        
        Args:
            model_name: "small", "medium", or "large" (from the mace-mp-0 foundation models).
            device: 'cpu' or 'cuda'/'mps'.
            default_dtype: Float precision for the neural network.
        """
        self.model_name = model_name
        self.device = device
        self.model_family = model_family
        self.backend_locator = backend_locator
        
        # Suppress extremely verbose MACE weight-loading output
        with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
            self.calc = _build_mace_calculator(model_family, model_name, device, default_dtype, backend_locator)

    def optimize_geometry(
        self,
        xyz_string: str,
        fmax: float = 0.01,
        max_steps: int = 500,
        drift_threshold: float = 1.0,
        *,
        bounded: bool = False,
        max_atom_drift_threshold: float = 0.8,
        sulfur_neighbor_delta_threshold: float = 0.25,
    ) -> str:
        """
        Given a starting XYZ geometry, optimize it using MACE and return the 
        converged XYZ geometry string.
        
        Args:
            xyz_string: Initial structure.
            fmax: Force convergence tolerance in eV/Å.
            max_steps: Maximum optimization steps.
            drift_threshold: Max allowed RMSD displacement (Å) before rejecting the result.
            
        Returns:
            The optimized Cartesian coordinates as an XYZ string.
        """
        # 1. Convert XYZ string to ASE Atoms object
        atoms = _read_xyz_atoms(xyz_string)
        initial_positions = atoms.positions.copy()

        # 2. Attach the MACE Neural Network Potential
        atoms.calc = self.calc
        
        # 3. Optimize using BFGS algorithm
        # Suppress ASE optimization step logs to keep output clean
        assert BFGS is not None
        with io.StringIO() as buf, redirect_stdout(buf):
            opt = BFGS(atoms, logfile=None)
            opt.run(fmax=fmax, steps=max_steps)
            
        # 4. Chemical Identity Guard: Check for excessive drift
        # This prevents MACE from distorting sensitive sulfur chemistries
        displacements = np.linalg.norm(atoms.positions - initial_positions, axis=1)
        max_drift = np.max(displacements)
        if max_drift > drift_threshold:
            logger.warning(f">>> [MLPOptimizer] WARNING: Excessive drift detected ({max_drift:.2f} Å). "
                  "MLP result may be chemically invalid. Reverting to original or suggest DFT.")

        # 5. Convert back to standard XYZ formatting
        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')
            optimized_xyz = f.getvalue()

        if bounded:
            metrics = compute_geometry_drift_metrics(xyz_string, optimized_xyz)
            if _should_reject_bounded_preopt(
                metrics,
                max_atom_drift_threshold=max_atom_drift_threshold,
                sulfur_neighbor_delta_threshold=sulfur_neighbor_delta_threshold,
            ):
                logger.warning(
                    ">>> [MLPOptimizer] Bounded geom_preopt rejected due to drift metrics: %s",
                    metrics,
                )
                return xyz_string
            
        return optimized_xyz

    def optimize_ts(self, xyz_string: str, fmax: float = 0.05, max_steps: int = 200) -> str:
        """
        Perform a rigorous Eigenvector Following saddle-point search using Sella.
        """
        from .ts_optimizer import TSOptimizer
        ts_opt = TSOptimizer(fmax=fmax, max_steps=max_steps)
        
        # 1. Convert to Atoms
        assert read is not None
        with io.StringIO(xyz_string) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
            
        # 2. Run Sella + MACE
        try:
            result = ts_opt.find_ts(atoms, self.calc)
        except ImportError as exc:
            logging.getLogger(__name__).warning(
                "TS optimizer backend unavailable; falling back to geometry minimization: %s",
                exc,
            )
            return self.optimize_geometry(xyz_string, fmax=fmax, max_steps=max_steps)
        if isinstance(result, list):
            atoms = result[-1]
        else:
            atoms = result
        
        # 3. Convert back to XYZ
        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')  # type: ignore
            return f.getvalue()
