"""
src/mlp_optimizer.py

Phase 10: MLP-Accelerated Geometry Optimization via MACE (mace-mp-0).

Replaces costly DFT structural optimizations with rapid Machine Learning Potential
optimizations. By default, relies on the `medium` MACE model via mace_mp, acting
as an ASE Calculator.
"""

import io
import logging
import numpy as np
from contextlib import redirect_stdout, redirect_stderr

try:
    from ase.io import read, write
    from ase.optimize import BFGS
    from mace.calculators import mace_mp
    _MLP_AVAILABLE = True
except ImportError:
    read = write = BFGS = mace_mp = None
    _MLP_AVAILABLE = False


class MLPOptimizer:
    """
    Wrapper for ASE-driven geometric optimization using the MACE neural network potential.
    """

    def __init__(self, model_name: str = "medium", device: str = "cpu", default_dtype: str = "float64"):
        """
        Initialize the MACE ASE Calculator.
        
        Args:
            model_name: "small", "medium", or "large" (from the mace-mp-0 foundation models).
            device: 'cpu' or 'cuda'/'mps'.
            default_dtype: Float precision for the neural network.
        """
        self.model_name = model_name
        self.device = device
        
        # Suppress extremely verbose MACE weight-loading output
        if not _MLP_AVAILABLE:
            raise ImportError("MLPOptimizer dependencies (mace-torch, ase) are not installed in the current environment.")

        assert mace_mp is not None
        with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
            self.calc = mace_mp(
                model=model_name,
                dispersion=False,
                default_dtype=default_dtype,
                device=device
            )

    def optimize_geometry(self, xyz_string: str, fmax: float = 0.01, max_steps: int = 500, drift_threshold: float = 1.0) -> str:
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
        assert read is not None
        with io.StringIO(xyz_string) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
            
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
            print(f">>> [MLPOptimizer] WARNING: Excessive drift detected ({max_drift:.2f} Å). "
                  "MLP result may be chemically invalid. Reverting to original or suggest DFT.")

        # 5. Convert back to standard XYZ formatting
        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')
            optimized_xyz = f.getvalue()
            
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
