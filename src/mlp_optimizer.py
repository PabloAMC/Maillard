"""
src/mlp_optimizer.py

Phase 10: MLP-Accelerated Geometry Optimization via MACE (mace-mp-0).

Replaces costly DFT structural optimizations with rapid Machine Learning Potential
optimizations. By default, relies on the `medium` MACE model via mace_mp, acting
as an ASE Calculator.
"""

import io
import tempfile
import numpy as np
from typing import Optional
from contextlib import redirect_stdout, redirect_stderr

try:
    from ase.io import read, write
    from ase.optimize import BFGS
    from mace.calculators import mace_mp
except ImportError:
    # Graceful degradation if ML dependencies are missing
    pass


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
        with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
            self.calc = mace_mp(
                model=model_name,
                dispersion=False,
                default_dtype=default_dtype,
                device=device
            )

    def optimize_geometry(self, xyz_string: str, fmax: float = 0.01, max_steps: int = 500) -> str:
        """
        Given a starting XYZ geometry, optimize it using MACE and return the 
        converged XYZ geometry string.
        
        Args:
            xyz_string: Initial structure.
            fmax: Force convergence tolerance in eV/Å.
            max_steps: Maximum optimization steps.
            
        Returns:
            The optimized Cartesian coordinates as an XYZ string.
        """
        # 1. Convert XYZ string to ASE Atoms object
        with tempfile.NamedTemporaryFile('w', suffix='.xyz') as tmp_in:
            tmp_in.write(xyz_string)
            tmp_in.flush()
            atoms = read(tmp_in.name, format='xyz')
            
        # 2. Attach the MACE Neural Network Potential
        atoms.calc = self.calc
        
        # 3. Optimize using BFGS algorithm
        # Suppress ASE optimization step logs to keep output clean
        with io.StringIO() as buf, redirect_stdout(buf):
            opt = BFGS(atoms, logfile=None)
            opt.run(fmax=fmax, steps=max_steps)
            
        # 4. Convert back to standard XYZ formatting
        with tempfile.NamedTemporaryFile('w+', suffix='.xyz') as tmp_out:
            write(tmp_out.name, atoms, format='xyz')
            tmp_out.seek(0)
            optimized_xyz = tmp_out.read()
            
        return optimized_xyz

    def optimize_ts(self, xyz_string: str, fmax: float = 0.05, max_steps: int = 200) -> str:
        """
        Placeholder for true Eigenvector Following saddle-point search (Sella - Phase 11).
        
        For Phase 10, this emits a warning and falls back to standard minimization
        so that the algorithmic routing logic in DFTRefiner can be tested.
        """
        print(">>> [MLPOptimizer] WARNING: Sella (Phase 11) is not yet integrated. "
              "Falling back to standard BFGS minimization. THIS IS NOT A TRUE TS SEARCH.")
        return self.optimize_geometry(xyz_string, fmax=fmax, max_steps=max_steps)
