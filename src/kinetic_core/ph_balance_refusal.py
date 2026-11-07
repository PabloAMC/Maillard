"""
The outputs that depend on the falsified pH balance are REFUSED (S3 application of M1b part (i), 2026-11-07).

WHAT DECIDED IT. M1b part (i) (``results/validation/kinetic_core_m1b_ph_balance_pentose_steps_prereg.md`` secs. 20 -
21, PR #94) found the shipped charge balance -- B0 + the shipped centre ledger (``ph_state``, ``sulfur.CENTRE_LEDGER``:
the amine fates and the carboxyl routing into CBX, CF-134's CO2 included) -- falsified under STATUS_AND_PLAN sec. 7's
rule S3 on CF-084's controls (Meynier 1995's buffered cysteine pots; Holm p <= 3e-3 across the defensible noise range,
independently confirmed, sec. 21.9). No rival is credited: B1 (amine survival) and B13 (B1 + acid supply) are carried,
not S3-tested at a converged fit. S3's refusal clause: the outputs that depend on a falsified incumbent "(those whose
sensitivity to the falsified structure's switch exceeds the lane's within-lab term) are refused, never shown at the
incumbent's point".

THE RULE (prereg sec. 21.10, rulings R1 - R3), applied per pot, per output. On every pot whose pH the balance solves
(the sulfur lane, any buffer but a clamped one), :func:`assess` runs the carried rivals beside B0 (the shipped point) at
their screen-passing draws' extremes: for each rival, the passing draw that holds the minimum and the one that holds
the maximum of each of its free coordinates, plus its best draw (``results/validation/ph_balance_s3_rival_points.json``,
written by ``build_ph_balance_s3_sensitivity.py --points`` from the screen's stored draws and A2's noise choices in D).
Each output is judged by the uncertainty it is shipped with (R1):

* an ABSOLUTE concentration, shipped as a level: its spread over B0 and every rival point, log10(max / min), against
  the lane's between-laboratory term tau (H0's 0.964 dex; prereg sec. 21.8.6). Above it the LEVEL is refused by name,
  as H0 refuses a thiol level: the number is still computed, kept in the arm, and serves ratios;
* a RELATIVE output (a ratio between two pots: ``compare``, the contrast diagnostic, the directional panel's steps):
  refused when the spread of its log ratio over the points exceeds the within-laboratory term omega (0.419 dex) or the
  points disagree in sign (S3), see :func:`ratio_refusal`;
* the cooled and in-situ pH the balance solves (R2): withheld when its spread exceeds the pH-level term, A1's run offset
  and per-point shape in quadrature, sqrt(0.45^2 + 0.16^2) = 0.478 pH (sec. 20.2), in any pot.

A value below ``ZERO_FLOOR_UG_PER_L`` is a zero: all zero is no spread, a zero beside a number an infinite one. A rival
whose integration exceeds the deterministic budget ``ph_balance_rivals.RIVAL_RHS_BUDGET`` cannot show independence, and
every output of that pot is refused with that reason. The rivals run only in the NumPy engine on the inline coupling
(not in JAX; W3). The verdict is decided once per (pot, targets), at the shipped centre, and cached; a draw reuses it.

WHAT IS NOT REFUSED. A pot whose pH is stated and held (``BufferSpec(kind="clamped")``: a pH-stat, or a measured pH
supplied as the held value; the measured-trajectory path) does not solve the balance and is unchanged. Nothing on the
trunk, acrylamide or lipid lanes reads the sulfur lane's pH. Numbers already sent outside the repository are not
re-scored (the WUR blind-test predictions of 2026-10-20 ran with a clamped pH).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from statistics import NormalDist
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from src import data_access, data_paths

#: THE SWITCH (config key ``ph.balance_refusal``): ``engine.predict``'s default. True ships; False only where the
#: model's own numbers are computed rather than shown (a frozen fit's objective, the sensitivity run).
ENABLED = True

FALSIFICATION_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_falsification.json"
SCREEN_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_screen.json"
A2_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_a2.json"
#: The rivals' points (R3): each rival's screen-passing draws' extremes plus its best draw.
POINTS_RECORD = data_paths.VALIDATION_DIR / "ph_balance_s3_rival_points.json"
#: The sensitivity record (``build_ph_balance_s3_sensitivity.py``): the rule's audit across the recipe space.
SENSITIVITY_RECORD = data_paths.VALIDATION_DIR / "ph_balance_s3_sensitivity.json"

#: The prefix every refusal of a balance-dependent output starts with (the scorers match on it).
BALANCE_REFUSED_TAG = "REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134)"
#: B3's latent acetic release barrier: ``ea_dpo_aa_kj_mol`` (Martins 2005's hexose analogue r_odg_aa), as the screen.
EA_ISO_KJ_MOL = 76.0
#: A log ratio this small has no sign: it is inside the numerical difference between the shipped fixed-point coupling
#: (B0) and the inline coupling the rivals run on, at most 0.0061 dex over the sensitivity record's concentration
#: outputs (``b0_vs_inline``; tests/unit/test_ph_balance_refusal.py keeps the floor above that maximum).
SIGN_FLOOR_DEX = 0.01
#: Below this a concentration is zero, not a number (``directional.ZERO_FLOOR_UG_PER_L``).
ZERO_FLOOR_UG_PER_L = 1e-6
#: The balance's pH outputs in a run's segment metadata (the in-situ START is the label's temperature shift before any
#: product exists, the shift the measured-trajectory path also uses, and is kept).
WITHHELD_PH_KEYS: Tuple[str, ...] = ("ph_in_situ_end", "ph_cooled_end")
#: What a failed or over-budget arm reads as, in place of its point values (JSON-safe: the scorecard stores it).
FAILED = "failed"

LIFT_TEXT = (
    "WHAT WOULD LIFT IT (STATUS_AND_PLAN sec. 7, S3's lift clause, pending independent confirmation): an S3 change of "
    "its own (items 6 and 9 replaced as S3 / S4 changes replace them) that ships the unweighted union of the remaining "
    "carried rivals' envelopes. A rival enters the union only if its converged best fit has passed S3's test on the "
    "data that falsified the incumbent (CF-084's controls); the union credits no rival. For B1 and B13 that also waits "
    "on the CO2 ledger audit (CF-134), which both inherit. Until then the refusal stands."
)


# ----------------------------------------------------------------------------------------------------------------
# Terms and points, from the records
# ----------------------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def thresholds() -> Dict[str, float]:
    """R1: tau for absolute levels, omega for relative outputs; R2: the pH-level term. With the two sensitivities the
    record reports (omega on absolute levels; the five-pot bound sigma_max,5 = 0.265 pH on the pH, prereg sec. 21.8.1)."""
    from . import thiol_relative_only

    r = thiol_relative_only.ruling()
    fals = data_access.load_json(FALSIFICATION_RECORD)
    five_pot = 0.2 / NormalDist().inv_cdf((1.0 + 0.05 ** (1.0 / 5.0)) / 2.0)
    return {"tau_dex": float(r.tau_abs_dex), "omega_dex": float(r.omega_dex),
            "ph_level": math.hypot(float(fals["s_off"]), float(fals["s_shape"])),
            "ph_five_pot_sensitivity": float(five_pot)}


@lru_cache(maxsize=1)
def rival_points() -> Dict[str, Tuple[Dict[str, Any], Any]]:
    """{point name: (Config changes, PhDrift)}: every rival point of the points record (R3), each built as
    ``fit_m1b_ph_balance.rival_changes`` builds a screen draw (checked equal in tests/unit/test_ph_balance_refusal.py)."""
    from . import ph_balance_rivals as pbr
    from .ph_state import PhDrift

    record = data_access.load_json(POINTS_RECORD)
    breakpoints = tuple(float(v) for v in record["davidek_breakpoints_in_situ"])
    out: Dict[str, Tuple[Dict[str, Any], Any]] = {}
    for rival, points in record["rivals"].items():
        for p in points:
            v = p["coordinates"]
            changes: Dict[str, Any] = {"ph.coupling": "inline",
                                       "ph.amine_survival": pbr.AmineSurvival(True, float(v["phi_cys"]),
                                                                              float(v["phi_arp"]))}
            if rival == "B13":
                changes["ph.acid_supply"] = pbr.AcidSupply(True, float(v["y5"]), float(v["y6"]), float(v["y_hi"]),
                                                           breakpoints, float(v["k_iso"]), EA_ISO_KJ_MOL, "hold_total")
                drift = PhDrift(acid_yield=0.0, arp_amine_pka=float(v["arp_amine_pka"]))
            elif rival == "B1":
                drift = PhDrift(acid_yield=float(v["acid_yield"]), arp_amine_pka=float(v["arp_amine_pka"]))
            else:
                raise ValueError(f"{POINTS_RECORD.name}: unknown rival {rival!r}")
            out[str(p["name"])] = (changes, drift)
    return out


def rival_point_names() -> Tuple[str, ...]:
    return tuple(rival_points())


def all_points() -> Tuple[str, ...]:
    return ("B0",) + rival_point_names()


# ----------------------------------------------------------------------------------------------------------------
# The verdict
# ----------------------------------------------------------------------------------------------------------------
def balance_solved(spec: Any, lanes: Sequence[str]) -> bool:
    """Whether this pot's pH is solved from the charge balance: the sulfur lane, with any buffer but a clamped one."""
    from .engine import SULFUR

    if SULFUR not in tuple(lanes or ()):
        return False
    buffer = getattr(spec.process, "buffer", None)
    return not (buffer is not None and buffer.is_clamped)


def log_spread(values: Sequence[float]) -> float:
    """log10(max / min) over the values; 0 when all are zero, inf when a zero sits beside a number."""
    positive = [float(v) for v in values if float(v) > ZERO_FLOOR_UG_PER_L]
    if not positive:
        return 0.0
    if len(positive) < len(values):
        return math.inf
    return math.log10(max(positive) / min(positive))


@dataclass
class Verdict:
    """The rule's reading of one (pot, targets) at the shipped centre."""

    #: {point: {compound: ug/L}} for B0 and every rival point (in memory only: it holds absolute levels)
    values: Dict[str, Dict[str, float]] = field(default_factory=dict)
    #: {point: {"ph_cooled_end": ..., "ph_in_situ_end": ...}}
    ph: Dict[str, Dict[str, Optional[float]]] = field(default_factory=dict)
    spread_dex: Dict[str, float] = field(default_factory=dict)
    #: the compounds whose absolute LEVEL is refused (spread above tau, or a failure)
    refused: Tuple[str, ...] = ()
    ph_spread: Dict[str, float] = field(default_factory=dict)
    ph_refused: bool = False
    failure: Optional[str] = None

    def summary(self) -> Dict[str, Any]:
        """What a run's metadata may carry: spreads and verdicts, no level."""
        th = thresholds()
        return {"rule": "S3 (STATUS_AND_PLAN sec. 7) on the M1b part (i) falsification; prereg sec. 21.10",
                "n_points": 1 + len(rival_point_names()), "spread_dex": dict(self.spread_dex),
                "tau_dex": th["tau_dex"], "level_refused": list(self.refused), "ph_spread": dict(self.ph_spread),
                "ph_level_term": th["ph_level"], "ph_refused": self.ph_refused, "failure": self.failure}


_CACHE: Dict[Tuple[str, Tuple[str, ...]], Verdict] = {}
_CACHE_MAX = 4096


def _key(spec: Any, targets: Sequence[str]) -> Tuple[str, Tuple[str, ...]]:
    return (repr(spec), tuple(str(t) for t in targets))


def cached_verdict(spec: Any, targets: Sequence[str]) -> Optional[Verdict]:
    return _CACHE.get(_key(spec, targets))


def _ph_of(run) -> Dict[str, Optional[float]]:
    seg = (run.run_metadata.get("segments") or [{}])[-1]
    return {k: seg.get(k) for k in WITHHELD_PH_KEYS}


def assess(spec: Any, targets: Sequence[str], b0_run: Any = None) -> Verdict:
    """The verdict for ``spec`` asked for ``targets``: B0 (``b0_run`` when the caller already integrated the shipped
    centre, else integrated here) against every rival point. Cached by (spec, targets)."""
    key = _key(spec, targets)
    if key in _CACHE:
        return _CACHE[key]
    from . import config, engine

    verdict = Verdict()
    if b0_run is None:
        b0_run = engine.predict(spec, list(targets), balance_refusal=False)
    names = [c for c in b0_run.concentrations_ug_per_l]
    verdict.values["B0"] = {c: float(v) for c, v in b0_run.concentrations_ug_per_l.items()}
    verdict.ph["B0"] = _ph_of(b0_run)
    try:
        base = config.Config.from_globals()
        for point, (changes, drift) in rival_points().items():
            with config.applied(base.replace(changes)):
                run = engine.predict(spec, list(targets), draw=engine.CoreDraw(ph_drift=drift), balance_refusal=False)
            verdict.values[point] = {c: float(v) for c, v in run.concentrations_ug_per_l.items()}
            verdict.ph[point] = _ph_of(run)
    except Exception as exc:  # noqa: BLE001 - a rival that cannot be run cannot show independence: refuse
        verdict.failure = f"{type(exc).__name__}: {exc}"[:300]
        verdict.refused = tuple(names)
        verdict.ph_refused = True
        _store(key, verdict)
        return verdict
    th = thresholds()
    for c in names:
        vals = [verdict.values[p].get(c) for p in all_points()]
        verdict.spread_dex[c] = math.inf if any(v is None for v in vals) else log_spread(vals)
    verdict.refused = tuple(c for c in names if verdict.spread_dex[c] > th["tau_dex"])
    for k in WITHHELD_PH_KEYS:
        vals = [verdict.ph[p].get(k) for p in all_points()]
        if all(v is not None for v in vals):
            verdict.ph_spread[k] = max(vals) - min(vals)
    verdict.ph_refused = any(s > th["ph_level"] for s in verdict.ph_spread.values())
    _store(key, verdict)
    return verdict


def _store(key, verdict: Verdict) -> None:
    if len(_CACHE) >= _CACHE_MAX:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = verdict


def point_values(spec: Any, targets: Sequence[str], run: Any, compound: str) -> Any:
    """{point: ug/L} of ``compound`` in one arm: the cached verdict's values where the balance is solved, else the
    run's own value at every point (a held pH does not move with the balance). ``FAILED`` for an arm whose rivals could
    not be run; ``None`` when nothing is available."""
    v = cached_verdict(spec, targets)
    if v is not None:
        if v.failure:
            return FAILED
        out = {p: v.values.get(p, {}).get(compound) for p in all_points()}
        return None if any(x is None for x in out.values()) else {p: float(x) for p, x in out.items()}
    if run is None or compound not in getattr(run, "concentrations_ug_per_l", {}):
        return None
    lanes = run.declaration.lanes or (run.declaration.lane,)
    if balance_solved(spec, lanes):
        return None
    value = float(run.concentrations_ug_per_l[compound])
    return {p: value for p in all_points()}


def offsets_dex(spec: Any, targets: Sequence[str], compound: str) -> Any:
    """{point: log10(value at point / value at B0)} for ``compound`` on a balance-solved pot whose verdict is cached:
    the rivals' values RELATIVE to B0, which a within-study ratio needs and which carry no absolute level. ``FAILED``
    when the rivals could not be run; ``None`` when there is no verdict or a value is zero beside a number (the ratio is
    then decided on zeros, below); all zero reads 0 at every point."""
    v = cached_verdict(spec, targets)
    if v is None:
        return None
    if v.failure:
        return FAILED
    vals = {p: v.values.get(p, {}).get(compound) for p in all_points()}
    if any(x is None for x in vals.values()):
        return None
    zero = [float(x) <= ZERO_FLOOR_UG_PER_L for x in vals.values()]
    if all(zero):
        return {p: 0.0 for p in vals}
    if any(zero):
        return FAILED
    return {p: math.log10(float(x) / float(vals["B0"])) for p, x in vals.items()}


def offsets_ratio_refusal(log_ratio_b0: float, off_a: Any, off_b: Any) -> Optional[str]:
    """:func:`ratio_refusal` from a B0 log ratio and each arm's :func:`offsets_dex` (``None`` = an arm the balance does
    not move: zero at every point)."""
    if off_a == FAILED or off_b == FAILED:
        return f"{BALANCE_REFUSED_TAG}: an arm's carried rivals could not be run or hit a zero. " + short_reason()
    names = tuple((off_a or off_b or {p: 0.0 for p in all_points()}).keys())
    zero = {p: 0.0 for p in names}
    a, b = off_a or zero, off_b or zero
    lin = {p: 10.0 ** (float(log_ratio_b0) + float(a[p]) - float(b[p])) for p in names}
    return ratio_refusal(lin, {p: 1.0 for p in names})


def ratio_refusal(a: Any, b: Any) -> Optional[str]:
    """Why the a / b ratio is refused, given each arm's {point: value} (:func:`point_values`): an arm's rivals failed,
    the points disagree in sign (S3: the direction is refused), or the log ratio's spread exceeds omega. ``None`` when
    it stands, or when an arm carries no point values (nothing to compare)."""
    if a == FAILED or b == FAILED:
        return (f"{BALANCE_REFUSED_TAG}: an arm's carried rivals could not be run, so the ratio's independence from the "
                "balance is not shown. " + short_reason())
    if a is None or b is None:
        return None
    logs = []
    for point in a:
        x, y = float(a[point]), float(b[point])
        zx, zy = x <= ZERO_FLOOR_UG_PER_L, y <= ZERO_FLOOR_UG_PER_L
        if zx and zy:
            logs.append(0.0)
        elif zx or zy:
            logs.append(math.inf if zy else -math.inf)
        else:
            logs.append(math.log10(x / y))
    nonzero = [v for v in logs if abs(v) > SIGN_FLOOR_DEX]
    if nonzero and min(nonzero) < 0.0 < max(nonzero):
        return (f"{BALANCE_REFUSED_TAG}: the DIRECTION of this ratio is refused: B0 and the carried rivals disagree in "
                "its sign (S3). " + short_reason())
    finite = all(math.isfinite(v) for v in logs)
    spread = (max(logs) - min(logs)) if finite else (0.0 if len(set(logs)) == 1 else math.inf)
    if spread > thresholds()["omega_dex"]:
        return (f"{BALANCE_REFUSED_TAG}: the spread of this log ratio over B0 and the carried rivals "
                f"({'inf' if spread == math.inf else f'{spread:.2f}'} dex) exceeds omega. " + short_reason())
    return None


def level_refused(declaration: Any, compound: str) -> bool:
    """Whether ``compound``'s absolute level is refused by this rule on ``declaration`` (a dataclass or its dict)."""
    names = (declaration.get("balance_refused_targets") if isinstance(declaration, Mapping)
             else getattr(declaration, "balance_refused_targets", ()))
    return str(compound) in tuple(names or ())


# ----------------------------------------------------------------------------------------------------------------
# The words
# ----------------------------------------------------------------------------------------------------------------
def short_reason() -> str:
    """The one-line form the scorecards carry per refused row."""
    th = thresholds()
    return (
        f"{BALANCE_REFUSED_TAG}: computed on the sulfur lane's pH solved from the shipped charge balance, which S3 "
        "falsified on CF-084's controls with no rival credited (M1b prereg sec. 21); the output moves with the carried "
        f"rivals B1 / B13 beyond the term it is shipped with (a level: tau = {th['tau_dex']:.3f} dex; a ratio or "
        f"direction: omega = {th['omega_dex']:.3f} dex or a sign; the pH: {th['ph_level']:.3f} pH; prereg sec. 21.10). "
        "Lifted by an S3 change shipping the union of the carried rivals' envelopes once each has passed S3's test on "
        "CF-084's controls (sec. 7, S3's lift clause)."
    )


def refusal_text(compounds: Sequence[str] = (), verdict: Optional[Verdict] = None) -> str:
    """The full reason for refused levels: what decided it, the spread, what still answers, and what would lift it."""
    th = thresholds()
    named = ""
    if compounds:
        parts = []
        for c in compounds:
            s = verdict.spread_dex.get(c) if verdict is not None else None
            parts.append(f"{c} (spread {'inf' if s == math.inf else f'{s:.2f}'} dex)" if s is not None else str(c))
        named = ", ".join(parts) + ": "
    failure = (f" The carried rivals could not be run on this pot ({verdict.failure}), so its independence from the "
               "balance is not shown." if verdict is not None and verdict.failure else "")
    return (
        f"{BALANCE_REFUSED_TAG}: {named}the absolute level is computed with the sulfur lane's pH solved from the shipped "
        "charge balance, and that structure, B0 + the shipped centre ledger (its amine fates and its carboxyl routing "
        "into CBX, CO2 included: CF-134), is FALSIFIED under STATUS_AND_PLAN sec. 7's rule S3 on CF-084's controls: "
        "Meynier 1995's buffered cysteine pots, Holm p <= 3e-3 across the defensible noise range, independently "
        "confirmed (M1b part (i) prereg sec. 21; results/validation/kinetic_core_m1b_part_i_a2_s3.md). No rival is "
        "credited; B1 (amine survival) and B13 (B1 + acid supply) are carried. S3 never shows an output that depends "
        "on the falsified structure at the incumbent's point: this level's spread over B0 and the carried rivals' "
        f"screen-passing extremes exceeds tau = {th['tau_dex']:.3f} dex, the term a level is shipped with (prereg sec. "
        "21.10, R1)." + failure + " Ratios of it between pots are answered by `compare` where they stand under omega "
        "and agree in sign. A pot whose pH is measured or held is answered: state the held pH as a clamped buffer "
        "(buffer kind 'clamped'), and the balance is not solved. " + LIFT_TEXT
    )


def ph_refusal_text() -> str:
    """Why the balance's pH is withheld from a run."""
    th = thresholds()
    return ("the pH solved from the charge balance is withheld: its spread over B0 and the carried rivals exceeds the "
            f"pH-level term {th['ph_level']:.3f} pH (prereg sec. 20.2; sec. 21.10, R2). " + short_reason())


__all__ = ["A2_RECORD", "BALANCE_REFUSED_TAG", "EA_ISO_KJ_MOL", "ENABLED", "FAILED", "FALSIFICATION_RECORD", "LIFT_TEXT",
           "POINTS_RECORD", "SCREEN_RECORD", "SENSITIVITY_RECORD", "SIGN_FLOOR_DEX", "Verdict", "WITHHELD_PH_KEYS",
           "ZERO_FLOOR_UG_PER_L", "all_points", "assess", "balance_solved", "cached_verdict", "level_refused",
           "log_spread", "offsets_dex", "offsets_ratio_refusal", "ph_refusal_text", "point_values", "ratio_refusal",
           "refusal_text", "rival_point_names", "rival_points", "short_reason", "thresholds"]
