"""
src/ts_optimizer.py

Phase 11: Sella Eigenvector-Following TS Search.
Provides a rigorous Transition State (TS) optimization engine using the 
Sella algorithm, compatible with any ASE calculator (MACE, PySCF, etc.).
"""

import io
from typing import Optional, Tuple
from contextlib import redirect_stdout

_SELLA_IMPORT_ERROR: Optional[Exception] = None

try:
    from ase import Atoms
except ImportError:
    Atoms = object  # type: ignore[assignment]

try:
    from sella import Sella
except Exception as exc:
    # Graceful degradation if Sella is unavailable or broken at import time.
    Sella = None
    _SELLA_IMPORT_ERROR = exc


class TSOptimizer:
    """
    Wrapper for the Sella transition state optimizer.
    """

    def __init__(self, fmax: float = 0.05, max_steps: int = 200):
        self.fmax = fmax
        self.max_steps = max_steps

    @staticmethod
    def probe_availability() -> Tuple[bool, str]:
        if Sella is None:
            if _SELLA_IMPORT_ERROR is None:
                return False, "Sella is not installed in the active environment."
            return False, f"Sella is unavailable in the active environment: {_SELLA_IMPORT_ERROR}"
        return True, ""

    def find_ts(self, atoms: Atoms, calculator, fmax: Optional[float] = None, max_steps: Optional[int] = None) -> Atoms:
        """
        Perform a transition state search starting from an initial Atoms object.
        
        Args:
            atoms: ASE Atoms object (initial TS guess).
            calculator: ASE Calculator to compute forces/energies.
            fmax: Force convergence tolerance.
            max_steps: Maximum optimization steps.
            
        Returns:
            Optimized ASE Atoms object.
        """
        if Sella is None:
            available, reason = self.probe_availability()
            if available:
                reason = "Sella is unavailable in the active environment."
            raise ImportError(reason)

        fmax = fmax or self.fmax
        max_steps = max_steps or self.max_steps

        # Attach calculator
        atoms.calc = calculator

        # Run Sella optimization
        # We suppress stdout to keep logs clean unless there is an error
        with io.StringIO() as buf, redirect_stdout(buf):
            dyn = Sella(atoms, logfile=None)  # type: ignore
            dyn.run(fmax=fmax, steps=max_steps)
            
        return atoms

    def is_converged(self, atoms: Atoms) -> bool:
        """
        Check if the atoms object reached convergence in the last Sella run.
        Directly checking forces since ASE dynamics objects don't always 
        persist the 'converged' flag in a simple way.
        """
        forces = atoms.get_forces()
        f_max = (forces**2).sum(axis=1).max()**0.5
        return f_max <= self.fmax
