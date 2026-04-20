"""P0 — TS-guess quality gate & refinement orchestrator.

Evaluates an incoming TS guess at xTB level and, if it fails the quality
gate (fmax > threshold or n_imag ≠ 1), attempts to improve it via a
cascade of lightweight methods before handing it to the expensive DFT/Sella
pipeline.

Cascade order:
  1. Relaxed scan along the forming bond (xTB, most robust for bimolecular)
  2. CI-NEB seeded with existing xTB path frames (xTB, fallback)
  3. Sella eigenvector-following TS search (xTB) to polish the best candidate
  4. pyGSM Double-Ended Growing String Method (opt-in, only when 1–3
     fail to deliver a sane candidate)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from .logger import get_logger
from .xtb_backend import (
    detect_forming_bond,
    probe_ts_guess_xtb,
    refine_ts_sella_xtb,
    run_xtb_cineb,
    run_xtb_relaxed_scan,
)

logger = get_logger(__name__)

# Sanity bounds (kcal/mol relative to reactant)
_ENERGY_FLOOR_KCAL = -5.0
_ENERGY_CEILING_KCAL = 200.0
_EV_TO_KCAL = 23.0609


@dataclass
class RefinementAttempt:
    method: str
    fmax_ev_ang: Optional[float] = None
    n_imag: Optional[int] = None
    energy_eh: Optional[float] = None
    accepted: bool = False
    reason: str = ""


@dataclass
class RefinementOutcome:
    """Result of the TS-guess refinement pipeline."""
    xyz: str
    accepted: bool
    final_method: str  # "original", "relaxed_scan", "cineb"
    initial_probe: Dict[str, Any] = field(default_factory=dict)
    attempts: List[RefinementAttempt] = field(default_factory=list)
    reason: str = ""



def _energy_above_ceiling(
    probe: Dict[str, Any],
    reactant_energy_eh: Optional[float],
) -> bool:
    """Return True if the candidate energy is absurdly high (> ceiling above reactant)."""
    if reactant_energy_eh is None or probe.get("energy_eh") is None:
        return False
    delta_kcal = (probe["energy_eh"] - reactant_energy_eh) * 627.509
    return delta_kcal > _ENERGY_CEILING_KCAL


def _is_sane_candidate(
    probe: Dict[str, Any],
    reactant_energy_eh: Optional[float],
    *,
    require_n_imag: bool = True,
) -> bool:
    """Check whether probe results are physically reasonable.

    Parameters
    ----------
    require_n_imag : bool
        If True (default), reject candidates with n_imag != 1.
        Set to False for relaxed-scan candidates where the scan maximum
        locates the TS region by construction and the xTB Hessian may
        not resolve the negative curvature (common for bimolecular).
    """
    n_imag = probe.get("n_imag")
    if n_imag is None:
        return False  # couldn't even compute Hessian
    if require_n_imag and n_imag != 1:
        return False
    if reactant_energy_eh is not None and probe.get("energy_eh") is not None:
        delta_kcal = (probe["energy_eh"] - reactant_energy_eh) * 627.509
        if delta_kcal < _ENERGY_FLOOR_KCAL or delta_kcal > _ENERGY_CEILING_KCAL:
            return False
    return True


def refine_if_needed(
    ts_xyz: str,
    reactant_xyz: str,
    product_xyz: Optional[str] = None,
    charge: int = 0,
    spin: int = 0,
    *,
    threshold_fmax: float = 1.0,
    path_frames_dir: Optional[str] = None,
    policy: str = "auto",  # "auto", "always", "off"
    enable_gsm_fallback: bool = False,
    gsm_work_dir: Optional[str] = None,
    gsm_nodes: int = 11,
    gsm_max_iters: int = 80,
    gsm_timeout_s: float = 1800.0,
) -> RefinementOutcome:
    """Evaluate and optionally refine a TS guess.

    Parameters
    ----------
    ts_xyz : str
        XYZ content of the TS guess.
    reactant_xyz : str
        XYZ content of the optimized reactant (for energy reference).
    product_xyz : str, optional
        XYZ content of the product (needed for forming-bond detection and NEB).
    charge, spin : int
        Molecular charge and spin multiplicity.
    threshold_fmax : float
        fmax (eV/Å) above which refinement is triggered.
    path_frames_dir : str, optional
        Directory containing ``xtbpath_*.xyz`` frames for NEB seeding.
    policy : str
        ``"auto"`` (refine only if fmax > threshold), ``"always"``
        (always refine), ``"off"`` (skip refinement, probe only).
    enable_gsm_fallback : bool
        When True, run pyGSM DE-GSM as a final cascade step if every
        cheaper method failed to produce a sane candidate.  Requires
        ``product_xyz``.
    gsm_work_dir : str, optional
        Scratch directory for pyGSM artifacts (per-node XYZs, audit JSON).
        Defaults to ``./gsm_scratch`` next to the current working dir.
    gsm_nodes, gsm_max_iters, gsm_timeout_s : int / float
        GSM tuning knobs forwarded to :class:`src.gsm_backend.GSMRunner`.

    Returns
    -------
    RefinementOutcome
        Contains the (possibly improved) XYZ, acceptance flag, and audit trail.
    """
    # ── 1. Probe the incoming guess ──────────────────────────────────────
    logger.info(">>> [P0] Probing TS guess quality at xTB level...")
    initial_probe = probe_ts_guess_xtb(ts_xyz, charge, spin)

    # Probe the reactant for energy reference
    reactant_probe = probe_ts_guess_xtb(reactant_xyz, charge, spin)
    reactant_energy_eh = reactant_probe.get("energy_eh")

    outcome = RefinementOutcome(
        xyz=ts_xyz,
        accepted=True,
        final_method="original",
        initial_probe=initial_probe,
    )

    if policy == "off":
        outcome.reason = "policy=off; skipping refinement"
        logger.info(f"    [P0] {outcome.reason}")
        return outcome

    fmax = initial_probe.get("fmax_ev_ang", float("inf"))
    needs_refinement = (
        policy == "always"
        or fmax > threshold_fmax
        or not _is_sane_candidate(initial_probe, reactant_energy_eh)
    )

    if not needs_refinement:
        outcome.reason = (
            f"TS guess passes gate: fmax={fmax:.3f} ≤ {threshold_fmax}, "
            f"n_imag={initial_probe.get('n_imag')}"
        )
        logger.info(f"    [P0] {outcome.reason}")
        return outcome

    logger.info(
        f"    [P0] TS guess needs refinement: fmax={fmax:.3f} eV/Å, "
        f"n_imag={initial_probe.get('n_imag')}, policy={policy}"
    )

    # ── 2. Cascade: relaxed scan → CI-NEB ────────────────────────────────
    best_candidate_xyz = ts_xyz
    best_candidate_fmax = fmax
    best_candidate_method = "original"
    best_candidate_probe = initial_probe

    # 2a. Relaxed scan (most robust for bimolecular)
    if product_xyz:
        logger.info(">>> [P0] Attempting relaxed-scan refinement...")
        try:
            scan_result = run_xtb_relaxed_scan(
                reactant_xyz, product_xyz, charge, spin,
            )
            if scan_result is not None:
                scan_probe = probe_ts_guess_xtb(scan_result["xyz"], charge, spin)
                attempt = RefinementAttempt(
                    method="relaxed_scan",
                    fmax_ev_ang=scan_probe.get("fmax_ev_ang"),
                    n_imag=scan_probe.get("n_imag"),
                    energy_eh=scan_probe.get("energy_eh"),
                )
                scan_fmax = scan_probe.get("fmax_ev_ang", float("inf"))

                if _is_sane_candidate(scan_probe, reactant_energy_eh):
                    # Best case: n_imag=1 and sane energy
                    if scan_fmax < best_candidate_fmax:
                        best_candidate_xyz = scan_result["xyz"]
                        best_candidate_fmax = scan_fmax
                        best_candidate_method = "relaxed_scan"
                        best_candidate_probe = scan_probe
                    attempt.accepted = True
                    attempt.reason = "sane candidate (n_imag=1)"
                elif scan_fmax < best_candidate_fmax and not _energy_above_ceiling(scan_probe, reactant_energy_eh):
                    # Relaxed-scan maximum locates the TS region by
                    # construction.  The xTB Hessian often fails to
                    # resolve n_imag for bimolecular reactions, and the
                    # energy delta vs an unoptimized reactant is
                    # meaningless.  Accept if fmax improved — Sella will
                    # do the real saddle-point search at DFT level.
                    best_candidate_xyz = scan_result["xyz"]
                    best_candidate_fmax = scan_fmax
                    best_candidate_method = "relaxed_scan"
                    best_candidate_probe = scan_probe
                    attempt.accepted = True
                    attempt.reason = (
                        f"scan-maximum accepted (n_imag="
                        f"{scan_probe.get('n_imag')}, fmax improved "
                        f"{fmax:.3f}→{scan_fmax:.3f})"
                    )
                    logger.info(
                        f"    [P0] Relaxed-scan candidate accepted "
                        f"(scan max, fmax {fmax:.3f}→{scan_fmax:.3f})"
                    )
                else:
                    attempt.reason = (
                        f"fmax not improved ({scan_fmax:.3f} ≥ "
                        f"{best_candidate_fmax:.3f})"
                    )
                outcome.attempts.append(attempt)
        except Exception as exc:
            logger.warning(f"    [P0] Relaxed scan failed: {exc}")
            outcome.attempts.append(
                RefinementAttempt(method="relaxed_scan", reason=f"exception: {exc}")
            )

    # 2b. CI-NEB (always tried — may find better TS region than scan)
    if product_xyz:
        logger.info(">>> [P0] Attempting CI-NEB refinement...")
        try:
            neb_result = run_xtb_cineb(
                reactant_xyz, product_xyz, charge, spin,
                path_frames_dir=path_frames_dir,
                max_steps=200,
                spring_constant=0.5,
            )
            if neb_result is not None:
                neb_probe = probe_ts_guess_xtb(neb_result["xyz"], charge, spin)
                attempt = RefinementAttempt(
                    method="cineb",
                    fmax_ev_ang=neb_probe.get("fmax_ev_ang"),
                    n_imag=neb_probe.get("n_imag"),
                    energy_eh=neb_probe.get("energy_eh"),
                )
                neb_fmax = neb_probe.get("fmax_ev_ang", float("inf"))

                if _is_sane_candidate(neb_probe, reactant_energy_eh):
                    if neb_fmax < best_candidate_fmax:
                        best_candidate_xyz = neb_result["xyz"]
                        best_candidate_fmax = neb_fmax
                        best_candidate_method = "cineb"
                        best_candidate_probe = neb_probe
                    attempt.accepted = True
                    attempt.reason = "sane candidate (n_imag=1)"
                elif neb_fmax < best_candidate_fmax and not _energy_above_ceiling(neb_probe, reactant_energy_eh):
                    # NEB climbing image — same relaxed gate as scan
                    best_candidate_xyz = neb_result["xyz"]
                    best_candidate_fmax = neb_fmax
                    best_candidate_method = "cineb"
                    best_candidate_probe = neb_probe
                    attempt.accepted = True
                    attempt.reason = (
                        f"NEB climbing-image accepted (n_imag="
                        f"{neb_probe.get('n_imag')}, fmax improved "
                        f"{fmax:.3f}→{neb_fmax:.3f})"
                    )
                    logger.info(
                        f"    [P0] CI-NEB candidate accepted "
                        f"(fmax {fmax:.3f}→{neb_fmax:.3f})"
                    )
                else:
                    attempt.reason = (
                        f"fmax not improved ({neb_fmax:.3f} ≥ "
                        f"{best_candidate_fmax:.3f}) or energy above ceiling"
                    )
                outcome.attempts.append(attempt)
        except Exception as exc:
            logger.warning(f"    [P0] CI-NEB failed: {exc}")
            outcome.attempts.append(
                RefinementAttempt(method="cineb", reason=f"exception: {exc}")
            )

    # 2c. xTB-level Sella refinement of best candidate so far
    if best_candidate_method != "original":
        logger.info(">>> [P0] Polishing best candidate with xTB Sella TS search...")
        try:
            sella_result = refine_ts_sella_xtb(
                best_candidate_xyz, charge, spin,
                fmax=0.05,
                max_steps=150,
                timeout_seconds=300.0,
            )
            if sella_result is not None:
                sella_probe = probe_ts_guess_xtb(sella_result["xyz"], charge, spin)
                sella_fmax = sella_probe.get("fmax_ev_ang", float("inf"))
                attempt = RefinementAttempt(
                    method="xtb_sella",
                    fmax_ev_ang=sella_fmax,
                    n_imag=sella_probe.get("n_imag"),
                    energy_eh=sella_probe.get("energy_eh"),
                )
                if sella_fmax < best_candidate_fmax and not _energy_above_ceiling(sella_probe, reactant_energy_eh):
                    best_candidate_xyz = sella_result["xyz"]
                    best_candidate_fmax = sella_fmax
                    best_candidate_method = f"{best_candidate_method}+xtb_sella"
                    best_candidate_probe = sella_probe
                    attempt.accepted = True
                    attempt.reason = (
                        f"xTB-Sella polished: fmax→{sella_fmax:.3f}, "
                        f"n_imag={sella_probe.get('n_imag')}"
                    )
                    logger.info(f"    [P0] {attempt.reason}")
                else:
                    attempt.reason = (
                        f"xTB-Sella did not improve fmax ({sella_fmax:.3f} "
                        f"≥ {best_candidate_fmax:.3f})"
                    )
                outcome.attempts.append(attempt)
        except Exception as exc:
            logger.warning(f"    [P0] xTB-Sella polishing failed: {exc}")
            outcome.attempts.append(
                RefinementAttempt(method="xtb_sella", reason=f"exception: {exc}")
            )

    # 2d. pyGSM DE-GSM — last-resort generator when 2a–2c failed to deliver
    #     a sane candidate.  Off by default (opt-in via enable_gsm_fallback)
    #     because it is significantly more expensive than the rest of the
    #     cascade and only pays off for reactions where the topology of the
    #     PES needs a multi-coordinate path search.
    needs_gsm = (
        enable_gsm_fallback
        and product_xyz
        and not _is_sane_candidate(best_candidate_probe, reactant_energy_eh)
    )
    if needs_gsm:
        logger.info(">>> [P0] Attempting pyGSM DE-GSM fallback...")
        try:
            from .gsm_backend import GSMRunner  # local import: heavy module
            gsm_dir = Path(gsm_work_dir) if gsm_work_dir else Path("gsm_scratch")
            runner = GSMRunner(
                work_dir=gsm_dir,
                charge=charge,
                spin=spin,
                n_nodes=gsm_nodes,
                max_iters=gsm_max_iters,
                timeout_s=gsm_timeout_s,
            )
            gsm_result = runner.run_de_gsm(reactant_xyz, product_xyz)
            attempt = RefinementAttempt(method="pygsm")
            if gsm_result.converged and gsm_result.ts_xyz:
                gsm_probe = probe_ts_guess_xtb(gsm_result.ts_xyz, charge, spin)
                attempt.fmax_ev_ang = gsm_probe.get("fmax_ev_ang")
                attempt.n_imag = gsm_probe.get("n_imag")
                attempt.energy_eh = gsm_probe.get("energy_eh")
                gsm_fmax = gsm_probe.get("fmax_ev_ang", float("inf"))
                if (
                    _is_sane_candidate(gsm_probe, reactant_energy_eh)
                    or (
                        gsm_fmax < best_candidate_fmax
                        and not _energy_above_ceiling(gsm_probe, reactant_energy_eh)
                    )
                ):
                    best_candidate_xyz = gsm_result.ts_xyz
                    best_candidate_fmax = gsm_fmax
                    best_candidate_method = "pygsm"
                    best_candidate_probe = gsm_probe
                    attempt.accepted = True
                    attempt.reason = (
                        f"pyGSM peak adopted (node {gsm_result.peak_index}, "
                        f"n_imag={gsm_probe.get('n_imag')}, fmax={gsm_fmax:.3f}); "
                        f"audit={gsm_result.audit_dir}"
                    )
                    logger.info(f"    [P0] {attempt.reason}")
                else:
                    attempt.reason = (
                        f"pyGSM ran but candidate failed sanity check "
                        f"(n_imag={gsm_probe.get('n_imag')}, fmax={gsm_fmax:.3f})"
                    )
            else:
                attempt.reason = f"pyGSM did not converge: {gsm_result.reason}"
            outcome.attempts.append(attempt)
        except Exception as exc:
            logger.warning(f"    [P0] pyGSM fallback failed: {exc}")
            outcome.attempts.append(
                RefinementAttempt(method="pygsm", reason=f"exception: {exc}")
            )

    # ── 3. Decision ──────────────────────────────────────────────────────
    if best_candidate_method != "original":
        outcome.xyz = best_candidate_xyz
        outcome.final_method = best_candidate_method
        outcome.accepted = True
        outcome.reason = (
            f"Replaced via {best_candidate_method}: "
            f"fmax {fmax:.3f}→{best_candidate_fmax:.3f}, "
            f"n_imag={best_candidate_probe.get('n_imag')}"
        )
        logger.info(f"    [P0] {outcome.reason}")
    else:
        # No refinement improved the guess.  Apply the "better bad guess
        # than no guess" policy: accept the original UNLESS it is clearly
        # a minimum (0 imaginary) or has absurd energy.
        is_clearly_invalid = (
            initial_probe.get("n_imag") == 0
            or (
                reactant_energy_eh is not None
                and initial_probe.get("energy_eh") is not None
                and (initial_probe["energy_eh"] - reactant_energy_eh) * 627.509 > _ENERGY_CEILING_KCAL
            )
        )
        if is_clearly_invalid:
            outcome.accepted = False
            outcome.reason = (
                f"All refinements failed and original is invalid "
                f"(n_imag={initial_probe.get('n_imag')}, "
                f"fmax={fmax:.3f})"
            )
            logger.warning(f"    [P0] REJECTED: {outcome.reason}")
        else:
            outcome.accepted = True
            outcome.reason = (
                f"No refinement improved guess; accepting original "
                f"(fmax={fmax:.3f}, n_imag={initial_probe.get('n_imag')})"
            )
            logger.info(f"    [P0] {outcome.reason}")

    return outcome
