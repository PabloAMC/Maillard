"""
Outputs that depend on thiamine's MFT step in a pot that also holds a thiol donor are REFUSED (STEP5_RECORD.md, second
amendment of 2026-11-26, owner; B66, ``results/validation/kinetic_core_b66_cysteine_thiamine_refusal_prereg.md``).

WHY. The thiamine route (``r_thi_hmp``, ``r_hmp_mft``) reads neither cysteine nor glutathione, so adding either to a
thiamine pot moves 2-methyl-3-furanthiol (MFT) by x1.0 - x1.3 in the engine. Measured (Grosch & Zeiler-Hilgart 1992,
100 C, pH 5.7, thiamine 60 uM): x7.0 with cysteine and no sugar (Table II), x3.7 with cysteine and ribose, x4.4 with
glutathione and ribose, x2.7 at the lowest cysteine dose and x0.38 at a 91-fold excess (Table III). Zhang et al. 2024
(Fig. 1: thiamine 44.5 mM + xylose, 115 C, pH 4.9) report x19 with cysteine and x15 with glutathione, in a pot where
the sugar route also makes MFT. A with / without ratio is judged by omega (0.42 dex): Grosch's x2.7 lies inside it
(0.33 dex from the engine), x3.7 and x4.4 just outside, x7.0 (0.85) and Zhang (1.2) far outside, and the fall at a
large excess has the other sign. No step inside cited bands reproduces these (B65's plug-in screen). So the incumbent
structure is falsified with no rival credited, and S3's clause for that case applies (STATUS_AND_PLAN, S3): the
outputs sensitive to it are refused, never shown at the incumbent's point.

THE RULE: a sensitivity rival, decided once per (pot, targets) at the shipped centre, in the form of
``thiamine_h2s_refusal``. On a sulfur-lane pot whose charge holds thiamine AND a thiol donor (cysteine, glutathione,
or TTCA, which releases cysteine in the pot) the pot is integrated a second time with the rate of ``r_hmp_mft``
multiplied by ``parameters_sulfur.CYSTEINE_THIAMINE_MFT_RIVAL.factor`` (7.0: the largest ratio measured where thiamine
is the pot's only MFT source). The multiplier is a sensitivity probe, not a model of the effect. Each requested output
whose log10 spread between the shipped point and the rival exceeds ``threshold_dex()`` (omega / 2, the thiamine rule's
per-arm threshold, for its reason) is REFUSED in that pot: dropped from the answer with a reason naming the sources and
what would lift it; a request whose every output is refused is ``out_of_envelope``.

DECLARED CHOICES (prereg sec. 3 and its amendment). Any amount of a donor triggers the rule: nothing is measured below
0.17 mM cysteine, so a trace over-refuses. There is no ratio or concentration domain: the measured effect follows the
cysteine concentration, not its ratio to thiamine (CF-089), and changes sign at a large excess. With the factor at 7.0
an output is refused where thiamine's route carries more than about 10 % of the MFT; below that share the published
effect cannot move it past the threshold. The rule states no direction for the error: in Bolton 1994's pot the engine
is already above the measured level. The three per-arm rules (S3's, the thiamine H2S rule, this one) are judged
separately; the spreads seen are bimodal (0.01 dex or less, or 0.3 and above), so their sum is not an issue in
practice.

WHAT IS NOT REFUSED. A pot without thiamine, a pot without a thiol donor, an output the rival does not move
(2-furfurylthiol, hydrogen sulfide, furfural: 0.00 dex in the prereg's probe), and every lane but the sulfur lane: the
shipped numbers there are bit-identical. The rival runs only in the NumPy engine, never in JAX (its config key is held
in the port).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Sequence, Tuple

#: THE SWITCH (config key ``sulfur.cysteine_thiamine_mft_refusal``): ``engine.predict``'s default. True ships; False
#: only where the model's own numbers are computed rather than shown (a rule's own integrations, a fit's objective).
ENABLED = True

#: The prefix every refusal by this rule starts with (the scorers match on it).
CYSTEINE_THIAMINE_REFUSED_TAG = (
    "REFUSED, A THIOL DONOR'S EFFECT ON THIAMINE'S MFT ROUTE UNDECIDED (Grosch & Zeiler-Hilgart 1992 and Zhang 2024 "
    "vs the engine)")
#: The species keys the rule needs in the charge: thiamine, and one of the thiol donors (cysteine; glutathione; TTCA,
#: which the engine opens to cysteine in the pot).
THIAMINE_KEY = "THI"
THIOL_DONOR_KEYS = ("Cys", "GSH", "TTCA")

LIFT_TEXT = (
    "WHAT WOULD LIFT IT: a cysteine series at fixed thiamine in one medium, with thiamine and cysteine followed over "
    "time, with and without a pentose (docs/reference/decisive_measurements.yml: "
    "cysteine_boost_thiol_yield_concentrated), or a step that reproduces Grosch's and Zhang's dose curves inside "
    "cited bands (STEP5_RECORD.md, second amendment of 2026-11-26)."
)


def threshold_dex() -> float:
    """The per-arm spread above which an output is refused: omega / 2, as ``thiamine_h2s_refusal`` and for its reason
    (an arm judged alone must be safe in any ratio it enters)."""
    from . import thiamine_h2s_refusal

    return thiamine_h2s_refusal.threshold_dex()


def applies(spec: Any, declaration: Any) -> bool:
    """Whether the rule must integrate the rival: the sulfur lane, thiamine and a thiol donor both in the charge."""
    from .engine import SULFUR

    lanes = tuple(getattr(declaration, "lanes", None) or (declaration.lane,))
    if SULFUR not in lanes:
        return False
    charge = dict(declaration.mapped_precursors or {})
    return (float(charge.get(THIAMINE_KEY, 0.0)) > 0.0
            and any(float(charge.get(k, 0.0)) > 0.0 for k in THIOL_DONOR_KEYS))


@dataclass
class Verdict:
    """The rule's reading of one (pot, targets) at the shipped centre."""

    #: {compound: log10 spread between the shipped point and the rival} (inf: a zero beside a number)
    spread_dex: Dict[str, float] = field(default_factory=dict)
    refused: Tuple[str, ...] = ()
    failure: Optional[str] = None

    def summary(self) -> Dict[str, Any]:
        """What a run's metadata may carry: spreads and verdicts, no level."""
        from .parameters_sulfur import CYSTEINE_THIAMINE_MFT_RIVAL as rival

        return {"rule": "a thiol donor on thiamine's MFT route, sensitivity rival (STEP5_RECORD second amendment "
                        "2026-11-26; B66)",
                "rival_factor": float(rival.factor),
                "spread_dex": {c: ("inf" if s == math.inf else s) for c, s in self.spread_dex.items()},
                "threshold_dex": threshold_dex(), "refused": list(self.refused), "failure": self.failure}


_CACHE: Dict[Tuple[str, Tuple[str, ...], str], Verdict] = {}
_CACHE_MAX = 4096


def _key(spec: Any, targets: Sequence[str]) -> Tuple[str, Tuple[str, ...], str]:
    """(spec, targets, the installed engine Config's digest): a verdict decided under one configuration (an S3 rival
    point, a fit's switches) is never read under another."""
    from . import config

    return (repr(spec), tuple(str(t) for t in targets), config.Config.from_globals().digest())


def cached_verdict(spec: Any, targets: Sequence[str]) -> Optional[Verdict]:
    return _CACHE.get(_key(spec, targets))


def _rival_config():
    from dataclasses import replace

    from . import config
    from . import parameters_sulfur as ps

    rival = replace(ps.CYSTEINE_THIAMINE_MFT_RIVAL, enabled=True)
    return config.Config.from_globals().replace({"sulfur.cysteine_thiamine_mft_rival": rival})


def _own_run(spec: Any, targets: Sequence[str]) -> Any:
    """The model's own number for this rule: no refusal of any kind applied."""
    from . import engine

    return engine.predict(spec, list(targets), balance_refusal=False, thiamine_refusal=False,
                          cysteine_thiamine_refusal=False)


def assess(spec: Any, targets: Sequence[str], shipped_run: Any = None) -> Verdict:
    """The verdict for ``spec`` asked for ``targets``: the shipped centre (``shipped_run`` when the caller already
    integrated it, else integrated here) against the rival. Cached by (spec, targets, Config)."""
    from . import config
    from .ph_balance_refusal import log_spread

    key = _key(spec, targets)
    if key in _CACHE:
        return _CACHE[key]
    if shipped_run is None:
        shipped_run = _own_run(spec, targets)
    shipped = {c: float(v) for c, v in shipped_run.concentrations_ug_per_l.items()}
    verdict = Verdict()
    try:
        with config.applied(_rival_config()):
            rival_run = _own_run(spec, targets)
        rival = {c: float(v) for c, v in rival_run.concentrations_ug_per_l.items()}
    except Exception as exc:  # noqa: BLE001 - a rival that cannot be run cannot show independence: refuse
        verdict.failure = f"{type(exc).__name__}: {exc}"[:300]
        verdict.spread_dex = {c: math.inf for c in shipped}
        verdict.refused = tuple(shipped)
        _store(key, verdict)
        return verdict
    limit = threshold_dex()
    for c, v in shipped.items():
        verdict.spread_dex[c] = math.inf if c not in rival else log_spread([v, rival[c]])
    verdict.refused = tuple(c for c in shipped if verdict.spread_dex[c] > limit)
    _store(key, verdict)
    return verdict


def _store(key, verdict: Verdict) -> None:
    if len(_CACHE) >= _CACHE_MAX:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = verdict


def refusal_text(compounds: Sequence[str] = (), verdict: Optional[Verdict] = None) -> str:
    """The reason for refused outputs: the conflict, the sources, the spread, and what would lift it."""
    from .parameters_sulfur import CYSTEINE_THIAMINE_MFT_RIVAL as rival

    named = ""
    if compounds:
        parts = []
        for c in compounds:
            s = verdict.spread_dex.get(c) if verdict is not None else None
            parts.append(f"{c} (spread {'inf' if s == math.inf else f'{s:.2f}'} dex)" if s is not None else str(c))
        named = ", ".join(parts) + ": "
    failure = (f" The rival could not be run on this pot ({verdict.failure}), so the output's independence from the "
               "step is not shown." if verdict is not None and verdict.failure else "")
    return (
        f"{CYSTEINE_THIAMINE_REFUSED_TAG}: {named}this pot is charged with thiamine and a thiol donor (cysteine, "
        "glutathione or TTCA), and this output depends on thiamine's step to 2-methyl-3-furanthiol (MFT), which the "
        "engine runs as if no donor were present (a donor moves its MFT by x1.0 - x1.3). Measured, from thiamine at "
        "100 C and pH 5.7 (Grosch & Zeiler-Hilgart 1992): x7.0 with cysteine and no sugar, x3.7 with cysteine and "
        "ribose, x4.4 with glutathione and ribose, x2.7 at the lowest cysteine dose (0.17 mM) and x0.38 at a 91-fold "
        "cysteine excess; in a concentrated thiamine + xylose pot at 115 C and pH 4.9 (Zhang et al. 2024): x19 with "
        "cysteine and x15 with glutathione. No step inside cited bands reproduces these ratios, so the mechanism is "
        "undecided and neither the size nor, at a large excess, the sign of the engine's error is known. With the "
        f"step's rate multiplied by {rival.factor:g} (the largest ratio measured where thiamine is the only MFT "
        f"source) the output moves by more than {threshold_dex():.2f} dex (omega / 2, so that no ratio it enters can "
        "move beyond omega), so it is not shown at the shipped value." + failure + " " + LIFT_TEXT
    )


__all__ = ["CYSTEINE_THIAMINE_REFUSED_TAG", "ENABLED", "LIFT_TEXT", "THIAMINE_KEY", "THIOL_DONOR_KEYS", "Verdict",
           "applies", "assess", "cached_verdict", "refusal_text", "threshold_dex"]
