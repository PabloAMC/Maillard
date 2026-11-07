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

THE RULE, APPLIED PER OUTPUT. Whether an output depends on the balance is a property of the pot: Meynier's cysteine
pots move by 0.2 - 2 dex between B0 and the rivals, a glycine + glucose pot or a pot buffered far from any titratable
product by less than 1e-3 (``results/validation/ph_balance_s3_sensitivity.md``, which ran every sulfur-lane panel pot,
every sulfur-lane directional arm and a stratified grid of CLI recipes). So the rule is not a list of species: on every
pot whose pH the balance solves (the sulfur lane, any buffer but a clamped one), :func:`assess` runs the two carried
rivals at their representative points -- B1 at phi_cys = 1 (its other coordinates at its best draw) and B13 at its best
draw, both from ``kinetic_core_m1b_part_i_falsification.json``, B3's in-situ breakpoints from the screen record, the
configuration the screen used -- and compares each output with B0's (the shipped point):

* a concentration is refused when its spread over {B0, B1, B13}, log10(max / min), exceeds the lane's within-lab term
  omega (H0's 0.419 dex; sec. 7, the S1 - S4 preamble). A value below ``ZERO_FLOOR_UG_PER_L`` is a zero: all zero is no
  spread, a zero beside a number is an infinite one;
* a ratio between two pots (``engine.compare``; the directional panel's steps) is refused when the spread of its log
  ratio exceeds omega or the points disagree in sign (S3: "rivals that disagree in sign make the lever's direction
  refused"), see :func:`ratio_refusal`;
* the cooled and in-situ pH the balance solves are withheld when their spread exceeds the prereg's pH term:
  sqrt(s_off^2 + s_shape^2) = 0.478 pH for an unbuffered pot (sec. 20.2), sigma_max = 0.135 pH for a declared buffer
  (sec. 21.8.1).

The rivals run only in the NumPy engine on the inline coupling (not in JAX; W3), at about 0.1 - 15 s per pot. The
verdict is decided once per (pot, targets), at the shipped centre, and cached; a Monte-Carlo draw reuses it.

HOW IT IS REFUSED, through the machinery the engine already has. A request whose every answerable target is refused
gets the reason below in ``declaration.reasons`` (``out_of_envelope``: no number, ``CorePrediction.require`` raises;
the scorecard, the envelope and the directional panel record the row as refused with it; ``compare`` reports "not
comparable"). A request with some targets refused answers the rest and names the refused ones, as the engine names an
unreachable target ("NOT ANSWERED, BY NAME", ``run_metadata["refused_targets"]``).

WHAT IS NOT REFUSED. A pot whose pH is stated and held (``BufferSpec(kind="clamped")``: a pH-stat, or a measured pH
supplied as the held value; the measured-trajectory path) does not solve the balance and is unchanged. Nothing on the
trunk, acrylamide or lipid lanes reads the sulfur lane's pH. Numbers already sent outside the repository are not
re-scored (the WUR blind-test predictions of 2026-10-20 ran with a clamped pH).
"""
from __future__ import annotations

import contextlib
import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from src import data_access, data_paths

#: THE SWITCH (config key ``ph.balance_refusal``): ``engine.predict``'s default. True ships; False only where the
#: model's own numbers are computed rather than shown (a frozen fit's objective, the sensitivity run).
ENABLED = True

FALSIFICATION_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_falsification.json"
SCREEN_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_screen.json"
A2_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_m1b_part_i_a2.json"
#: The sensitivity record (``build_ph_balance_s3_sensitivity.py``): the rule's audit across the recipe space.
SENSITIVITY_RECORD = data_paths.VALIDATION_DIR / "ph_balance_s3_sensitivity.json"

#: The prefix every refusal of a balance-dependent output starts with (the scorers match on it).
BALANCE_REFUSED_TAG = "REFUSED, FALSIFIED pH BALANCE (S3; M1b part (i), CF-084, CF-134)"
#: The carried rivals' representative points, in the order they are run.
RIVAL_POINTS: Tuple[str, ...] = ("B1_phi_cys_1", "B13_best")
#: B3's latent acetic release barrier: ``ea_dpo_aa_kj_mol`` (Martins 2005's hexose analogue r_odg_aa), as the screen.
EA_ISO_KJ_MOL = 76.0
#: A log ratio this small has no sign: it is inside the numerical difference between the shipped fixed-point coupling
#: (B0) and the inline coupling the rivals run on, at most 0.0061 dex over the sensitivity record's 3523 concentration
#: outputs (``b0_vs_inline``; tests/unit/test_ph_balance_refusal.py keeps the floor above that maximum).
SIGN_FLOOR_DEX = 0.01
#: Below this a concentration is zero, not a number (``directional.ZERO_FLOOR_UG_PER_L``).
ZERO_FLOOR_UG_PER_L = 1e-6
#: The balance's pH outputs in a run's segment metadata (the in-situ START is the label's temperature shift before any
#: product exists, the shift the measured-trajectory path also uses, and is kept).
WITHHELD_PH_KEYS: Tuple[str, ...] = ("ph_in_situ_end", "ph_cooled_end")

LIFT_TEXT = (
    "WHAT WOULD LIFT IT (STATUS_AND_PLAN sec. 7, S3's lift clause, pending independent confirmation): an S3 change of "
    "its own (items 6 and 9 replaced as S3 / S4 changes replace them) that ships the unweighted union of the remaining "
    "carried rivals' envelopes. A rival enters the union only if its converged best fit has passed S3's test on the "
    "data that falsified the incumbent (CF-084's controls); the union credits no rival. For B1 and B13 that also waits "
    "on the CO2 ledger audit (CF-134), which both inherit. Until then the refusal stands."
)


# ----------------------------------------------------------------------------------------------------------------
# Thresholds and points, from the records
# ----------------------------------------------------------------------------------------------------------------
@lru_cache(maxsize=1)
def thresholds() -> Dict[str, float]:
    """omega (H0) for concentrations and ratios; the prereg's pH terms (A1 sec. 20.2, A2 sec. 21.8.1)."""
    from . import thiol_relative_only

    fals = data_access.load_json(FALSIFICATION_RECORD)
    a2 = data_access.load_json(A2_RECORD)
    return {"omega_dex": float(thiol_relative_only.ruling().omega_dex),
            "ph_unbuffered": math.hypot(float(fals["s_off"]), float(fals["s_shape"])),
            "ph_buffered": float(a2["sigma_max"]["0.2"]["sigma_max"])}


@lru_cache(maxsize=1)
def rival_points() -> Dict[str, Tuple[Dict[str, Any], Any]]:
    """{point: (Config changes, PhDrift)} for the two carried rivals, as ``fit_m1b_ph_balance.rival_changes`` builds
    them (checked equal in tests/unit/test_ph_balance_refusal.py)."""
    from . import ph_balance_rivals as pbr
    from .ph_state import PhDrift

    best = data_access.load_json(FALSIFICATION_RECORD)["best_draws"]
    breakpoints = tuple(float(v) for v in data_access.load_json(SCREEN_RECORD)["davidek_breakpoints_in_situ"])
    b1 = dict(best["B1"]["coordinates"])
    b13 = dict(best["B13"]["coordinates"])
    return {
        "B1_phi_cys_1": (
            {"ph.coupling": "inline", "ph.amine_survival": pbr.AmineSurvival(True, 1.0, float(b1["phi_arp"]))},
            PhDrift(acid_yield=float(b1["acid_yield"]), arp_amine_pka=float(b1["arp_amine_pka"]))),
        "B13_best": (
            {"ph.coupling": "inline",
             "ph.amine_survival": pbr.AmineSurvival(True, float(b13["phi_cys"]), float(b13["phi_arp"])),
             "ph.acid_supply": pbr.AcidSupply(True, float(b13["y5"]), float(b13["y6"]), float(b13["y_hi"]),
                                              breakpoints, float(b13["k_iso"]), EA_ISO_KJ_MOL, "hold_total")},
            PhDrift(acid_yield=0.0, arp_amine_pka=float(b13["arp_amine_pka"]))),
    }


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


def ph_tolerance(spec: Any) -> float:
    buffer = getattr(spec.process, "buffer", None)
    th = thresholds()
    return th["ph_buffered"] if (buffer is not None and buffer.kind != "none") else th["ph_unbuffered"]


@dataclass
class Verdict:
    """The rule's reading of one (pot, targets) at the shipped centre."""

    #: {point: {compound: ug/L}} for B0 and the rivals (in memory only: never serialised, it holds absolute levels)
    values: Dict[str, Dict[str, float]] = field(default_factory=dict)
    #: {point: {"ph_cooled_end": ..., "ph_in_situ_end": ...}}
    ph: Dict[str, Dict[str, Optional[float]]] = field(default_factory=dict)
    spread_dex: Dict[str, float] = field(default_factory=dict)
    refused: Tuple[str, ...] = ()
    ph_spread: Dict[str, float] = field(default_factory=dict)
    ph_refused: bool = False
    failure: Optional[str] = None

    def summary(self) -> Dict[str, Any]:
        """What a run's metadata may carry: spreads and verdicts, no level."""
        return {"rule": "S3 (STATUS_AND_PLAN sec. 7) on the M1b part (i) falsification; ph_balance_refusal",
                "points": ["B0"] + list(RIVAL_POINTS), "spread_dex": dict(self.spread_dex),
                "omega_dex": thresholds()["omega_dex"], "refused": list(self.refused),
                "ph_spread": dict(self.ph_spread), "ph_refused": self.ph_refused, "failure": self.failure}


#: A rival run longer than this is stopped and the pot's outputs refused (its independence is not shown). In the
#: sensitivity run every rival point that finished took at most 49.5 s (median 0.43 s; Xu 2007's pH 8 arm), and eight
#: ground past the generator's 300 s limit (Chen 2000 / 2002's pH 5 ribose + cysteine pots, 180 C): 120 s separates the
#: two with room and gives every recorded pot the record's verdict. Held only in the main thread (SIGALRM).
RIVAL_TIMEOUT_S = 120


@contextlib.contextmanager
def _time_limit(seconds: int):
    import signal
    import threading

    if threading.current_thread() is not threading.main_thread() or not hasattr(signal, "SIGALRM"):
        yield
        return

    def _raise(signum, frame):
        raise TimeoutError("a carried rival's integration ran past RIVAL_TIMEOUT_S")

    previous = signal.signal(signal.SIGALRM, _raise)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


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
    centre, else integrated here) against the two carried rivals. Cached by (spec, targets)."""
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
            with config.applied(base.replace(changes)), _time_limit(RIVAL_TIMEOUT_S):
                run = engine.predict(spec, list(targets), draw=engine.CoreDraw(ph_drift=drift), balance_refusal=False)
            verdict.values[point] = {c: float(v) for c, v in run.concentrations_ug_per_l.items()}
            verdict.ph[point] = _ph_of(run)
    except Exception as exc:  # noqa: BLE001 - a rival that cannot be run cannot show independence: refuse
        verdict.failure = f"{type(exc).__name__}: {exc}"[:300]
        verdict.refused = tuple(names)
        verdict.ph_refused = True
        _store(key, verdict)
        return verdict
    omega = thresholds()["omega_dex"]
    for c in names:
        vals = [verdict.values[p].get(c) for p in ("B0",) + RIVAL_POINTS]
        spread = math.inf if any(v is None for v in vals) else log_spread(vals)
        verdict.spread_dex[c] = spread
    verdict.refused = tuple(c for c in names if verdict.spread_dex[c] > omega)
    tol = ph_tolerance(spec)
    for k in WITHHELD_PH_KEYS:
        vals = [verdict.ph[p].get(k) for p in ("B0",) + RIVAL_POINTS]
        if all(v is not None for v in vals):
            verdict.ph_spread[k] = max(vals) - min(vals)
    verdict.ph_refused = any(s > tol for s in verdict.ph_spread.values())
    _store(key, verdict)
    return verdict


def _store(key, verdict: Verdict) -> None:
    if len(_CACHE) >= _CACHE_MAX:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = verdict


def point_values(spec: Any, targets: Sequence[str], run: Any, compound: str) -> Optional[Dict[str, float]]:
    """{point: ug/L} of ``compound`` in one arm: the cached verdict's values where the balance is solved, else the
    run's own value at every point (a held pH does not move with the balance). ``None`` when not available."""
    v = cached_verdict(spec, targets)
    if v is not None:
        out = {p: v.values.get(p, {}).get(compound) for p in ("B0",) + RIVAL_POINTS}
        return None if any(x is None for x in out.values()) else {p: float(x) for p, x in out.items()}
    if run is None or compound not in getattr(run, "concentrations_ug_per_l", {}):
        return None
    lanes = run.declaration.lanes or (run.declaration.lane,)
    if balance_solved(spec, lanes):
        return None
    value = float(run.concentrations_ug_per_l[compound])
    return {p: value for p in ("B0",) + RIVAL_POINTS}


def offsets_dex(spec: Any, targets: Sequence[str], compound: str) -> Optional[Dict[str, float]]:
    """{point: log10(value at point / value at B0)} for ``compound`` on a balance-solved pot whose verdict is cached:
    the rivals' values RELATIVE to B0, which a within-study ratio needs and which carry no absolute level. ``None``
    when there is no verdict or a value is zero beside a number; all zero reads 0 at every point."""
    v = cached_verdict(spec, targets)
    if v is None:
        return None
    vals = {p: v.values.get(p, {}).get(compound) for p in ("B0",) + RIVAL_POINTS}
    if any(x is None for x in vals.values()):
        return None
    zero = [float(x) <= ZERO_FLOOR_UG_PER_L for x in vals.values()]
    if all(zero):
        return {p: 0.0 for p in vals}
    if any(zero):
        return None
    return {p: math.log10(float(x) / float(vals["B0"])) for p, x in vals.items()}


def offsets_ratio_refusal(log_ratio_b0: float, off_a: Optional[Mapping[str, float]],
                          off_b: Optional[Mapping[str, float]]) -> Optional[str]:
    """:func:`ratio_refusal` from a B0 log ratio and each arm's :func:`offsets_dex` (``None`` = an arm the balance
    does not move: zero at every point)."""
    zero = {p: 0.0 for p in ("B0",) + RIVAL_POINTS}
    a, b = off_a or zero, off_b or zero
    lin = {p: 10.0 ** (float(log_ratio_b0) + float(a[p]) - float(b[p])) for p in zero}
    return ratio_refusal(lin, {p: 1.0 for p in zero})


def ratio_refusal(a: Optional[Mapping[str, float]], b: Optional[Mapping[str, float]]) -> Optional[str]:
    """Why the a / b ratio is refused, given each arm's {point: value} (:func:`point_values`): the points disagree in
    sign (S3: the direction is refused) or the log ratio's spread exceeds omega. ``None`` when it stands, or when an arm
    carries no point values (nothing to compare)."""
    if a is None or b is None:
        return None
    logs = []
    for point in ("B0",) + RIVAL_POINTS:
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


# ----------------------------------------------------------------------------------------------------------------
# The words
# ----------------------------------------------------------------------------------------------------------------
def short_reason() -> str:
    """The one-line form the scorecards carry per refused row."""
    th = thresholds()
    return (
        f"{BALANCE_REFUSED_TAG}: computed on the sulfur lane's pH solved from the shipped charge balance, which S3 "
        "falsified on CF-084's controls with no rival credited (M1b prereg sec. 21); this output's spread over B0 and "
        f"the carried rivals B1 / B13 exceeds omega = {th['omega_dex']:.3f} dex "
        "(results/validation/ph_balance_s3_sensitivity.md). Lifted by an S3 change shipping the union of the carried "
        "rivals' envelopes once each has passed S3's test on CF-084's controls (sec. 7, S3's lift clause)."
    )


def refusal_text(compounds: Sequence[str] = (), verdict: Optional[Verdict] = None) -> str:
    """The full reason: what decided it, the spread that refused it, what still answers, and what would lift it."""
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
        f"{BALANCE_REFUSED_TAG}: {named}computed with the sulfur lane's pH solved from the shipped charge balance, and "
        "that structure, B0 + the shipped centre ledger (its amine fates and its carboxyl routing into CBX, CO2 "
        "included: CF-134), is FALSIFIED under STATUS_AND_PLAN sec. 7's rule S3 on CF-084's controls: Meynier 1995's "
        "buffered cysteine pots, Holm p <= 3e-3 across the defensible noise range, independently confirmed (M1b part "
        "(i) prereg sec. 21; results/validation/kinetic_core_m1b_part_i_a2_s3.md). No rival is credited; B1 (amine "
        "survival) and B13 (B1 + acid supply) are carried. S3 refuses an output whose sensitivity to that switch "
        "exceeds the lane's within-laboratory term and never shows it at the incumbent's point: here the spread of the "
        f"output over B0, B1 (phi_cys = 1) and B13 exceeds omega = {th['omega_dex']:.3f} dex "
        "(results/validation/ph_balance_s3_sensitivity.md runs the rule across the recipe space)." + failure
        + " A pot whose pH is measured or held is answered: state the held pH as a clamped buffer (buffer kind "
        "'clamped'), and the balance is not solved. " + LIFT_TEXT
    )


def ph_refusal_text() -> str:
    """Why the balance's pH is withheld from a run."""
    th = thresholds()
    return ("the pH solved from the charge balance is withheld: its spread over B0 and the carried rivals exceeds the "
            f"prereg's pH term ({th['ph_unbuffered']:.3f} pH unbuffered, {th['ph_buffered']:.3f} pH buffered). "
            + short_reason())


__all__ = ["A2_RECORD", "BALANCE_REFUSED_TAG", "EA_ISO_KJ_MOL", "ENABLED", "FALSIFICATION_RECORD", "LIFT_TEXT", "RIVAL_POINTS",
           "SCREEN_RECORD", "SENSITIVITY_RECORD", "SIGN_FLOOR_DEX", "Verdict", "WITHHELD_PH_KEYS", "ZERO_FLOOR_UG_PER_L", "assess",
           "balance_solved", "cached_verdict", "log_spread", "offsets_dex", "offsets_ratio_refusal", "point_values", "ph_refusal_text", "ph_tolerance", "ratio_refusal",
           "refusal_text", "rival_points", "short_reason", "thresholds"]
