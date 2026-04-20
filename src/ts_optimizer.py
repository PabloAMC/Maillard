"""
src/ts_optimizer.py

Phase 11: Sella Eigenvector-Following TS Search.
Provides a rigorous Transition State (TS) optimization engine using the 
Sella algorithm, compatible with any ASE calculator (MACE, PySCF, etc.).
"""

import time
import numpy as np
from pathlib import Path
from typing import Optional, Tuple

from .logger import get_logger

_SELLA_IMPORT_ERROR: Optional[Exception] = None

try:
    from ase import Atoms
except ImportError:
    Atoms = object  # type: ignore[assignment]

try:
    from sella import Sella
except Exception as exc:
    Sella = None
    _SELLA_IMPORT_ERROR = exc

logger = get_logger(__name__)


class TSOptimizationStall(RuntimeError):
    """Raised when TS search stalls or exceeds wall-time limit."""
    pass


class TSPlateauConverged(TSOptimizationStall):
    """Raised when the TS search has flat-lined on a numerical plateau.

    Signals the caller that further geometry optimization is futile and that
    the best geometry seen so far should be used directly for single-point /
    Hessian post-processing instead of being escalated to a fallback optimizer.
    """
    pass


class TSOptimizer:
    """Wrapper for the Sella transition state optimizer."""

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

    def find_ts(
        self,
        atoms: Atoms,
        calculator,
        fmax: Optional[float] = None,
        max_steps: Optional[int] = None,
        trajectory_dir: Optional[str] = None,
        timeout_seconds: float = 7200,
        stall_window: int = 20,
        delta0: float = 0.01,
        v0: Optional[np.ndarray] = None,
        internal: bool = False,
    ) -> Atoms:
        """
        Perform a transition state search.

        Parameters
        ----------
        v0 : array, optional
            Initial eigenvector (TS mode) for Sella to follow.  Shape
            ``(3*N,)`` where *N* is the number of atoms.  When supplied
            Sella uses this direction as its initial guess for the
            imaginary mode instead of relying on Hessian updates alone.
        internal : bool
            If True, Sella works in automatically-detected internal
            (delocalized) coordinates, which improves mode identification
            for bond-breaking / bond-forming reactions.

        Raises TSOptimizationStall on wall-time timeout or force-norm stall.
        """
        if Sella is None:
            available, reason = self.probe_availability()
            if available:
                reason = "Sella is unavailable in the active environment."
            raise ImportError(reason)

        fmax = fmax or self.fmax
        max_steps = max_steps or self.max_steps
        atoms.calc = calculator

        # Trajectory and logging
        logfile = None
        traj_path = None
        if trajectory_dir:
            td = Path(trajectory_dir)
            td.mkdir(parents=True, exist_ok=True)
            logfile = str(td / "sella.log")
            traj_path = str(td / "sella.traj")

        sella_kwargs: dict = dict(
            logfile=logfile,
            trajectory=traj_path,
            delta0=delta0,
            order=1,
            internal=internal,
        )
        if v0 is not None:
            sella_kwargs["v0"] = np.asarray(v0, dtype=float)
            logger.info("    [SELLA] Using supplied v0 (initial TS mode) for eigenvector-following.")
        dyn = Sella(atoms, **sella_kwargs)  # type: ignore

        # Watchdog: wall-time + stall detection via cached forces
        start = time.monotonic()
        force_history: list[float] = []

        def _watchdog():
            elapsed = time.monotonic() - start
            if elapsed > timeout_seconds:
                raise TSOptimizationStall(
                    f"TS search exceeded {timeout_seconds:.0f}s wall time ({elapsed:.0f}s elapsed)"
                )
            cached = getattr(atoms.calc, "results", {}).get("forces")
            if cached is not None:
                fmax_now = float(np.sqrt((cached**2).sum(axis=1).max()))
                force_history.append(fmax_now)
                logger.info(f"    [SELLA] step {len(force_history):3d}  fmax={fmax_now:.5f}  t={elapsed:.0f}s")
                if len(force_history) >= stall_window:
                    recent = force_history[-stall_window:]
                    if min(recent) > 0.95 * recent[0]:
                        raise TSOptimizationStall(
                            f"TS search stalled: fmax={fmax_now:.5f} not improving "
                            f"over {stall_window} steps"
                        )

        dyn.attach(_watchdog)
        dyn.run(fmax=fmax, steps=max_steps)
        return atoms

    def is_converged(self, atoms: Atoms) -> bool:
        forces = atoms.get_forces()
        f_max = float((forces**2).sum(axis=1).max()**0.5)
        return f_max <= self.fmax
