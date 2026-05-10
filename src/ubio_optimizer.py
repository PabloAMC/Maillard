"""
src/ubio_optimizer.py

Alternative Geometry Pre-relaxation backend using the UBio-MolFM-V1 
Equivariant Transformer (E2Former-V2) specifically tuned for protein 
and bio-molecular transition states.
"""

import io
import os
import sys
import logging
import tempfile
import numpy as np
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr

from src.logger import get_logger

logger = get_logger(__name__)

# Dynamically inject the UBio-MolFM source code path so we can import it
repo_path = Path("models/ubio_repo/src").resolve()
if repo_path.exists() and str(repo_path) not in sys.path:
    sys.path.insert(0, str(repo_path))

try:
    from ase.io import read, write
    from ase.optimize import BFGS
    from molfm.interface.ase.calculator.e2former_calculator import E2FormerCalculator
    _UBIO_AVAILABLE = True
except ImportError as exc:
    read = write = BFGS = E2FormerCalculator = None
    _UBIO_AVAILABLE = False
    _UBIO_ERROR = str(exc)

class UBioOptimizer:
    """
    Wrapper for ASE-driven geometric optimization using the UBio-MolFM-V1 
    neural network potential. Designed specifically for peptide and amino acid interactions.
    """

    def __init__(self, checkpoint_path: str = "models/ubio-molfm/molfm-v1-stage-3.pt", config_name: str = "models/ubio-molfm/config.yaml", device: str = "cpu"):
        self.device = device
        
        if not _UBIO_AVAILABLE:
            raise ImportError(f"UBioOptimizer dependencies (ase, molfm, e3nn, torchaudio etc.) are not installed. Details: {_UBIO_ERROR}")

        checkpoint_abs = Path(checkpoint_path).resolve()
        config_abs = Path(config_name).resolve()
        
        if not checkpoint_abs.exists():
            raise FileNotFoundError(f"UBio-MolFM checkpoint not found at {checkpoint_abs}")
        if not config_abs.exists():
            raise FileNotFoundError(f"UBio-MolFM config not found at {config_abs}")

        logger.info(f"Loading UBio-MolFM Calculator from {checkpoint_abs}...")
        
        with io.StringIO() as buf, redirect_stdout(buf), redirect_stderr(buf):
            # Currently use_compile=False for compatibility, can toggle True if PyTorch >= 2.0
            self.calc = E2FormerCalculator(
                checkpoint_path=str(checkpoint_abs),
                config_name=str(config_abs),
                head_name="omol25", # default regression head in the UBio Stage-3 model
                device=device,
                use_tf32=False, # Disable to prevent precision drift
                use_compile=False
            )
        logger.info("UBio-MolFM Calculator initialized successfully.")

    def optimize_geometry(self, xyz_string: str, fmax: float = 0.01, max_steps: int = 500, drift_threshold: float = 1.0) -> str:
        """
        Optimize using UBio-MolFM-V1 and return the converged XYZ.
        """
        assert read is not None
        with io.StringIO(xyz_string) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
            
        # Optional: Setup cell/pbc if the model strictly requires it
        atoms.set_cell([100, 100, 100]) # Large vacuum box
        atoms.pbc = [False, False, False]
        
        initial_positions = atoms.positions.copy()

        atoms.calc = self.calc
        
        assert BFGS is not None
        with io.StringIO() as buf, redirect_stdout(buf):
            opt = BFGS(atoms, logfile=None)
            opt.run(fmax=fmax, steps=max_steps)
            
        displacements = np.linalg.norm(atoms.positions - initial_positions, axis=1)
        max_drift = np.max(displacements)
        if max_drift > drift_threshold:
            logger.warning(f">>> [UBioOptimizer] WARNING: Excessive drift detected ({max_drift:.2f} Å). "
                 "Returning optimized but recommend double-checking structural integrity.")

        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')
            optimized_xyz = f.getvalue()
            
        return optimized_xyz

    def optimize_ts(self, xyz_string: str, fmax: float = 0.05, max_steps: int = 200) -> str:
        """
        Perform Saddle-Point Search using Sella paired with UBio-MolFM.
        """
        from .ts_optimizer import TSOptimizer
        ts_opt = TSOptimizer(fmax=fmax, max_steps=max_steps)
        
        assert read is not None
        with io.StringIO(xyz_string) as f:
            atoms = read(f, format='xyz')
            if isinstance(atoms, list):
                atoms = atoms[-1]
            
        atoms.set_cell([100, 100, 100])
        atoms.pbc = [False, False, False]
        
        try:
            result = ts_opt.find_ts(atoms, self.calc)
        except (ImportError, RuntimeError) as exc:
            logging.getLogger(__name__).warning("Sella TS search failed; falling back to minimisation: %s", exc)
            return self.optimize_geometry(xyz_string, fmax=fmax, max_steps=max_steps)
            
        if isinstance(result, list):
            atoms = result[-1]
        else:
            atoms = result
        
        assert write is not None
        with io.StringIO() as f:
            write(f, atoms, format='xyz')  # type: ignore
            return f.getvalue()
