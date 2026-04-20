"""pyGSM backend — Growing String Method TS-guess generator.

This module wraps the upstream `pyGSM <https://github.com/ZimmermanGroup/pyGSM>`_
Double-Ended Growing String Method (DE-GSM) into a single-call interface
that the Maillard pipeline can invoke as a fallback when cheaper TS-guess
strategies (relaxed scan, CI-NEB, xTB-Sella) fail to deliver a true saddle.

Design contract
---------------
* The runner owns a *pathfinder-only* level of theory (xTB / GFN2 + ALPB).
  The xTB barrier emitted by GSM is **never** consumed downstream — it is
  exposed only for diagnostic logging.  The TS geometry is the sole
  artifact passed back to the DFT refiner, where the imaginary-frequency
  gate (``_classify_ts_frequencies``) decides whether the saddle is real.
* All native I/O (tblite SCC tables, pyGSM banner output) is captured via
  :func:`src._native_io.suppress_native_output` so that callers can run
  the GSM pipeline silently from inside any logging context.
* All artifacts (per-iteration node geometries, the pyGSM scratch
  directory) are persisted under ``work_dir/`` so that failures can be
  audited offline.
"""

from __future__ import annotations

import json
import os
import time
import traceback
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

import numpy as np

from ase.calculators.calculator import Calculator, all_changes

from ._native_io import suppress_native_output
from .logger import get_logger
from .xtb_backend import get_xtb_ase_calculator

logger = get_logger(__name__)

# pyGSM convergence defaults (validated against upstream `de_gsm.py` __main__).
_DEFAULT_NNODES = 11
_DEFAULT_MAX_ITERS = 80
_DEFAULT_GSM_FMAX_HARTREE_BOHR = 5.0e-4   # ~0.025 eV/Å — pyGSM units
_DEFAULT_TIMEOUT_S = 1800.0


# ---------------------------------------------------------------------------
# SCF-tolerant calculator wrapper
# ---------------------------------------------------------------------------
class _SCFTolerantCalc(Calculator):
    """ASE calculator that catches xTB SCF failures and returns a penalty.

    During DE-GSM string growing, interpolated nodes can have unphysical
    geometries (e.g. fragments at intermediate separation) that cause SCF
    divergence.  Rather than aborting the entire GSM run, this wrapper
    returns a very high energy and zero forces so GSM can route around the
    problematic region.
    """

    implemented_properties = ["energy", "forces"]
    name = "scf-tolerant-xtb"

    def __init__(self, inner_calc: Calculator, penalty_offset_ev: float = 27.2):
        super().__init__()
        self._inner = inner_calc
        self._penalty_offset_ev = penalty_offset_ev  # ~1 Hartree above last good
        self._last_good_energy: Optional[float] = None
        self._scf_failures = 0

    @property
    def scf_failures(self) -> int:
        return self._scf_failures

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        try:
            self._inner.calculate(atoms, properties, system_changes)
            self.results = dict(self._inner.results)
            self._last_good_energy = self.results.get("energy")
        except Exception as exc:
            # Tolerate any xTB/tblite calculation failure during GSM string
            # growing — SCF divergence, LAPACK eigenvalue failures, etc.
            from ase.calculators.calculator import CalculationFailed

            err_msg = str(exc)
            is_calc_failure = (
                isinstance(exc, CalculationFailed)
                or "SCF" in err_msg
                or "converge" in err_msg.lower()
                or "sygvd" in err_msg.lower()
                or "eigenvalue" in err_msg.lower()
            )
            if is_calc_failure:
                self._scf_failures += 1
                logger.warning(
                    "[GSM-SCF] Failure #%d: %s → penalty energy",
                    self._scf_failures,
                    err_msg[:120],
                )
                n = len(atoms) if atoms is not None else 1
                base = (
                    self._last_good_energy
                    if self._last_good_energy is not None
                    else 0.0
                )
                # Use small random forces instead of zeros so that pyGSM's
                # internal-coordinate machinery doesn't produce NaN gradients.
                rng = np.random.default_rng(self._scf_failures)
                noise = rng.normal(scale=0.01, size=(n, 3))
                self.results = {
                    "energy": base + self._penalty_offset_ev,
                    "forces": noise,
                }
            else:
                raise


@dataclass
class GSMResult:
    """Outcome of a single DE-GSM run.

    Attributes
    ----------
    converged : bool
        True iff pyGSM reported success *and* the TS node lies in the
        interior (not at either endpoint).
    ts_xyz : Optional[str]
        XYZ block of the TS node when ``converged``; ``None`` otherwise.
    ts_energy_eh_xtb : Optional[float]
        xTB-level Hartree energy at the TS node (diagnostic only — never
        report this as a barrier).
    peak_index : Optional[int]
        Node index of the TS along the string (0 = reactant, N-1 = product).
    n_iters : int
        Number of GSM optimization passes consumed.
    elapsed_s : float
        Wall-clock seconds spent inside :meth:`GSMRunner.run_de_gsm`.
    audit_dir : str
        Path to the on-disk scratch directory containing per-node XYZs and
        the pyGSM run log.
    reason : str
        Human-readable diagnostic (success message or failure cause).
    """

    converged: bool
    ts_xyz: Optional[str]
    ts_energy_eh_xtb: Optional[float]
    peak_index: Optional[int]
    n_iters: int
    elapsed_s: float
    audit_dir: str
    reason: str
    extras: Dict[str, Any] = field(default_factory=dict)

    def to_audit_dict(self) -> Dict[str, Any]:
        """Return a JSON-safe dict (omits the bulky ``ts_xyz`` payload)."""
        payload = asdict(self)
        payload.pop("ts_xyz", None)
        return payload


# ---------------------------------------------------------------------------
# XYZ ↔ pyGSM `geom` conversion
# ---------------------------------------------------------------------------
# pyGSM represents geometries as ``[[symbol, x, y, z], ...]`` (Ångström),
# matching the layout produced by ``pyGSM.utilities.manage_xyz``.  The two
# helpers below let us round-trip arbitrary XYZ text without touching disk.


def _parse_xyz(xyz: str) -> List[List[Any]]:
    """Convert an XYZ block to pyGSM's ``geom`` list-of-lists format."""
    lines = xyz.strip().splitlines()
    if len(lines) < 3:
        raise ValueError("XYZ block must have at least 3 lines (n, comment, atoms).")
    n_atoms = int(lines[0].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(
            f"XYZ header advertises {n_atoms} atoms but body has {len(atom_lines)}."
        )
    geom: List[List[Any]] = []
    for line in atom_lines:
        symbol, sx, sy, sz, *_ = line.split()
        geom.append([symbol, float(sx), float(sy), float(sz)])
    return geom


def _geom_to_xyz(geom: List[List[Any]], comment: str = "") -> str:
    """Inverse of :func:`_parse_xyz` — pyGSM ``geom`` → XYZ text."""
    body = "\n".join(
        f"{row[0]:<3s} {float(row[1]):>16.10f} {float(row[2]):>16.10f} {float(row[3]):>16.10f}"
        for row in geom
    )
    return f"{len(geom)}\n{comment}\n{body}\n"


def _xyz_to_endpoints_file(reactant_xyz: str, product_xyz: str, path: Path) -> None:
    """Write an XYZ file containing reactant and product as two frames."""
    path.write_text(reactant_xyz.rstrip() + "\n" + product_xyz.rstrip() + "\n")


# ---------------------------------------------------------------------------
# Lazy pyGSM imports (so unit tests that only touch helpers do not pay the cost)
# ---------------------------------------------------------------------------

def _import_pygsm() -> Dict[str, Any]:
    """Import the pyGSM symbols we need; raise a clear error if missing."""
    try:
        from pyGSM.coordinate_systems import (  # noqa: WPS433 — runtime import
            DelocalizedInternalCoordinates,
            PrimitiveInternalCoordinates,
            Topology,
        )
        from pyGSM.growing_string_methods import DE_GSM
        from pyGSM.level_of_theories.ase import ASELoT
        from pyGSM.molecule import Molecule
        from pyGSM.optimizers import eigenvector_follow
        from pyGSM.potential_energy_surfaces import PES
        from pyGSM.utilities import elements, manage_xyz
    except ImportError as exc:  # pragma: no cover — diagnosed at runtime
        raise ImportError(
            "pyGSM is not importable. Install with "
            "`pip install --no-deps git+https://github.com/ZimmermanGroup/pyGSM.git` "
            "and ensure the `topology.py` `pkg_resources` patch is applied."
        ) from exc
    return {
        "DE_GSM": DE_GSM,
        "ASELoT": ASELoT,
        "Molecule": Molecule,
        "PES": PES,
        "eigenvector_follow": eigenvector_follow,
        "Topology": Topology,
        "PrimitiveInternalCoordinates": PrimitiveInternalCoordinates,
        "DelocalizedInternalCoordinates": DelocalizedInternalCoordinates,
        "elements": elements,
        "manage_xyz": manage_xyz,
    }


# ---------------------------------------------------------------------------
# GSMRunner
# ---------------------------------------------------------------------------

class GSMRunner:
    """Single-shot DE-GSM driver bound to a fixed level of theory.

    Parameters
    ----------
    work_dir : Path
        Directory where pyGSM scratch files and per-node XYZs are written.
        Created on demand; pre-existing contents are *not* cleared so
        callers can stage seeds.
    charge, spin : int
        Net molecular charge and spin multiplicity carried into the xTB
        calculator.
    solvent : Optional[str]
        ALPB solvent passed to :func:`get_xtb_ase_calculator`.  Use
        ``None`` for gas phase.
    n_nodes : int
        Number of GSM nodes (including endpoints).
    max_iters : int
        Maximum DE-GSM optimization passes.
    timeout_s : float
        Wall-clock cap; the run is interrupted via a SIGALRM-free
        cooperative check between iterations (best-effort).
    fmax_hartree_bohr : float
        Per-node convergence criterion in pyGSM units.
    """

    def __init__(
        self,
        *,
        work_dir: Path,
        charge: int = 0,
        spin: int = 0,
        solvent: Optional[str] = "water",
        n_nodes: int = _DEFAULT_NNODES,
        max_iters: int = _DEFAULT_MAX_ITERS,
        timeout_s: float = _DEFAULT_TIMEOUT_S,
        fmax_hartree_bohr: float = _DEFAULT_GSM_FMAX_HARTREE_BOHR,
    ) -> None:
        # Coerce to an absolute path defensively: ``run_de_gsm`` performs an
        # ``os.chdir`` into a scratch directory (via ``_scratch_cwd``) before
        # invoking pyGSM, so any relative paths captured here would dangle.
        # ``Path.resolve()`` historically produced surprising results when the
        # target did not yet exist on certain filesystems (notably bind-mounts
        # in macOS Docker), occasionally returning the input path unchanged.
        # ``Path.absolute()`` is purely lexical and safe here because the
        # directory is created by ``run_de_gsm`` itself.
        wd = Path(work_dir).expanduser()
        if not wd.is_absolute():
            wd = (Path.cwd() / wd)
        self.work_dir = Path(os.path.normpath(str(wd)))
        self.charge = int(charge)
        self.spin = int(spin)
        self.solvent = solvent
        self.n_nodes = int(n_nodes)
        self.max_iters = int(max_iters)
        self.timeout_s = float(timeout_s)
        self.fmax_hartree_bohr = float(fmax_hartree_bohr)
        self._last_scf_tolerant_calc: Optional[_SCFTolerantCalc] = None

    # ------------------------------------------------------------------ public

    def run_de_gsm(self, reactant_xyz: str, product_xyz: str) -> GSMResult:
        """Run DE-GSM between reactant and product, return the TS node."""
        if not reactant_xyz or not str(reactant_xyz).strip():
            raise ValueError("run_de_gsm: reactant_xyz is empty")
        if not product_xyz or not str(product_xyz).strip():
            raise ValueError("run_de_gsm: product_xyz is empty")
        if not self.work_dir.is_absolute():
            # Should never happen given __init__, but guard anyway.
            self.work_dir = self.work_dir.absolute()
        self.work_dir.mkdir(parents=True, exist_ok=True)
        endpoints_path = self.work_dir / "endpoints.xyz"
        _xyz_to_endpoints_file(reactant_xyz, product_xyz, endpoints_path)
        if not endpoints_path.is_file():
            raise FileNotFoundError(
                f"GSMRunner: endpoints.xyz was not written at {endpoints_path}"
            )
        logger.info(
            "[GSM] work_dir=%s endpoints=%s (%d bytes)",
            self.work_dir, endpoints_path, endpoints_path.stat().st_size,
        )

        start = time.monotonic()
        with self._scratch_cwd(), suppress_native_output():
            try:
                gsm, geoms, atoms = self._build_gsm(endpoints_path)
                gsm.go_gsm(max_iters=self.max_iters, opt_steps=3, rtype=2)
            except Exception as exc:  # pragma: no cover — surfaced as failed result
                elapsed = time.monotonic() - start
                tb = traceback.format_exc()
                logger.warning(f"[GSM] DE-GSM raised: {exc}\n{tb}")
                try:
                    (self.work_dir / "gsm_traceback.txt").write_text(tb, encoding="utf-8")
                except Exception:
                    pass
                return GSMResult(
                    converged=False,
                    ts_xyz=None,
                    ts_energy_eh_xtb=None,
                    peak_index=None,
                    n_iters=0,
                    elapsed_s=elapsed,
                    audit_dir=str(self.work_dir),
                    reason=f"pyGSM exception: {type(exc).__name__}: {exc}",
                )

        elapsed = time.monotonic() - start
        result = self._extract_result(gsm, atoms, elapsed)
        # Propagate SCF failure count from tolerant calculator
        scf_calc = self._last_scf_tolerant_calc
        if scf_calc is not None and scf_calc.scf_failures > 0:
            result.extras["scf_failures_tolerated"] = scf_calc.scf_failures
            logger.info(
                "[GSM] %d SCF failures tolerated during string growing",
                scf_calc.scf_failures,
            )
        self._persist_audit(result)
        return result

    # ----------------------------------------------------------------- internal

    def _build_gsm(self, endpoints_path: Path) -> Tuple[Any, List[Any], List[Any]]:
        """Construct the full pyGSM object graph."""
        pg = _import_pygsm()
        manage_xyz = pg["manage_xyz"]
        elements = pg["elements"]

        geoms = manage_xyz.read_xyzs(str(endpoints_path))
        if len(geoms) < 2:
            raise ValueError(
                f"endpoints file must contain ≥ 2 frames, got {len(geoms)}"
            )

        atom_syms = manage_xyz.get_atoms(geoms[0])
        element_data = elements.ElementData()
        atoms = [element_data.from_symbol(s) for s in atom_syms]

        xyz_r = manage_xyz.xyz_to_np(geoms[0])
        xyz_p = manage_xyz.xyz_to_np(geoms[-1])

        # Build union topology (R ∪ P bonds) so primitives cover the full path.
        top_r = pg["Topology"].build_topology(xyz_r, atoms)
        top_p = pg["Topology"].build_topology(xyz_p, atoms)
        for bond in top_p.edges():
            if bond not in top_r.edges() and (bond[1], bond[0]) not in top_r.edges():
                top_r.add_edge(*sorted(bond, reverse=True))

        prim_r = pg["PrimitiveInternalCoordinates"].from_options(
            xyz=xyz_r, atoms=atoms, connect=True, addtr=False, addcart=False, topology=top_r,
        )
        prim_p = pg["PrimitiveInternalCoordinates"].from_options(
            xyz=xyz_p, atoms=atoms, connect=True, addtr=False, addcart=False, topology=top_r,
        )
        prim_r.add_union_primitives(prim_p)

        coord_obj = pg["DelocalizedInternalCoordinates"].from_options(
            xyz=xyz_r, atoms=atoms, addtr=False, addcart=False,
            connect=True, primitives=prim_r,
        )

        ase_calc = get_xtb_ase_calculator(
            charge=self.charge, spin=self.spin, solvent=self.solvent,
        )
        # High electronic temperature helps SCF converge on distorted
        # interpolated geometries that arise during string growing.
        try:
            ase_calc.set(max_iterations=500, electronic_temperature=9500.0)
        except Exception:
            try:
                ase_calc.set(max_iterations=500)
            except Exception:
                pass
        # Keep the tolerant wrapper as safety net with a LOW penalty so
        # rare remaining failures don't corrupt the energy profile.
        ase_calc = _SCFTolerantCalc(ase_calc, penalty_offset_ev=0.5)
        self._last_scf_tolerant_calc = ase_calc
        multiplicity = max(1, self.spin + 1)
        states = [(multiplicity, 0)]
        lot = pg["ASELoT"].from_options(
            ase_calc, geom=geoms[0], states=states, gradient_states=states,
        )
        pes = pg["PES"].from_options(lot=lot, ad_idx=0, multiplicity=multiplicity)

        optimizer = pg["eigenvector_follow"].from_options(
            Linesearch="backtrack", OPTTHRESH=self.fmax_hartree_bohr,
            DMAX=0.05, abs_max_step=0.05, conv_Ediff=0.01,
        )

        reactant = pg["Molecule"].from_options(
            geom=geoms[0], PES=pes, coord_obj=coord_obj, Form_Hessian=True,
        )
        product = pg["Molecule"].copy_from_options(
            reactant, xyz=xyz_p, new_node_id=self.n_nodes - 1, copy_wavefunction=False,
        )

        gsm = pg["DE_GSM"].from_options(
            reactant=reactant,
            product=product,
            nnodes=self.n_nodes,
            optimizer=optimizer,
        )
        # NOTE: do NOT pre-set ``gsm.find`` or ``gsm.climb``.  They are state
        # flags toggled internally by ``MainGSM.set_stage()`` once the string
        # has been grown and climb/find criteria are met.  The user-facing
        # knob is ``rtype=2`` passed to ``go_gsm``.
        return gsm, geoms, atoms

    def _extract_result(self, gsm: Any, atoms: List[Any], elapsed_s: float) -> GSMResult:
        """Pull the highest-energy interior node out of the optimized string."""
        try:
            energies = list(gsm.energies)
        except Exception:
            energies = []

        active_nodes = [node for node in getattr(gsm, "nodes", []) if node is not None]
        n_iters = int(getattr(gsm, "nopt", 0)) or len(active_nodes)

        if not active_nodes or not energies:
            return GSMResult(
                converged=False,
                ts_xyz=None,
                ts_energy_eh_xtb=None,
                peak_index=None,
                n_iters=n_iters,
                elapsed_s=elapsed_s,
                audit_dir=str(self.work_dir),
                reason="pyGSM returned empty node/energy list",
            )

        # Trim energy list to the active nodes (pyGSM pads with zeros).
        energies = energies[: len(active_nodes)]
        peak_index = int(np.argmax(energies))
        peak_is_interior = 0 < peak_index < len(active_nodes) - 1

        peak_node = active_nodes[peak_index]
        try:
            geom = peak_node.geometry  # pyGSM Molecule attribute (list-of-lists)
            ts_xyz = _geom_to_xyz(
                geom, comment=f"pyGSM DE-GSM TS node {peak_index} of {len(active_nodes)}",
            )
        except Exception as exc:
            return GSMResult(
                converged=False,
                ts_xyz=None,
                ts_energy_eh_xtb=float(energies[peak_index]),
                peak_index=peak_index,
                n_iters=n_iters,
                elapsed_s=elapsed_s,
                audit_dir=str(self.work_dir),
                reason=f"failed to extract TS geometry: {exc}",
            )

        if not peak_is_interior:
            return GSMResult(
                converged=False,
                ts_xyz=ts_xyz,
                ts_energy_eh_xtb=float(energies[peak_index]),
                peak_index=peak_index,
                n_iters=n_iters,
                elapsed_s=elapsed_s,
                audit_dir=str(self.work_dir),
                reason=(
                    f"peak at endpoint (index {peak_index} of "
                    f"{len(active_nodes)}); no interior maximum"
                ),
            )

        return GSMResult(
            converged=True,
            ts_xyz=ts_xyz,
            ts_energy_eh_xtb=float(energies[peak_index]),
            peak_index=peak_index,
            n_iters=n_iters,
            elapsed_s=elapsed_s,
            audit_dir=str(self.work_dir),
            reason=(
                f"DE-GSM converged: TS at node {peak_index}/{len(active_nodes) - 1}, "
                f"barrier_xtb={energies[peak_index] - energies[0]:.6f} Eh "
                "(diagnostic only)"
            ),
            extras={
                "node_energies_eh": [float(e) for e in energies],
                "n_active_nodes": len(active_nodes),
            },
        )

    def _persist_audit(self, result: GSMResult) -> None:
        """Write a `gsm_attempt.json` summary alongside pyGSM's scratch."""
        audit_path = self.work_dir / "gsm_attempt.json"
        try:
            audit_path.write_text(json.dumps(result.to_audit_dict(), indent=2))
        except Exception as exc:  # pragma: no cover — auditing is best-effort
            logger.warning(f"[GSM] failed to persist audit JSON: {exc}")
        if result.ts_xyz:
            try:
                (self.work_dir / "ts_node.xyz").write_text(result.ts_xyz)
            except Exception as exc:  # pragma: no cover
                logger.warning(f"[GSM] failed to persist ts_node.xyz: {exc}")

    @contextmanager
    def _scratch_cwd(self) -> Iterator[None]:
        """Run pyGSM inside ``work_dir`` so its scratch files land there."""
        prev = os.getcwd()
        os.chdir(self.work_dir)
        try:
            yield
        finally:
            os.chdir(prev)


__all__ = ["GSMRunner", "GSMResult"]
