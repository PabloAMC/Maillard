"""
Outputs that depend on thiamine's H2S share where its zero is UNMEASURED are REFUSED (STEP5_RECORD.md, amendment of
2026-11-17, owner).

WHY. ``parameters_sulfur.thiamine_h2s_share`` (B50) is 0 below pH 6.25 at every temperature. That zero is measured only
near 100 C (Dwivedi & Arnold 1971, a ~96 - 100 C bath). At 120 C the two assays of one laboratory disagree by more than
25x at pH 6.0: Dwivedi & Arnold 1971 (methylene blue) read ~0.36 of the degraded thiamine as H2S, Dwivedi & Arnold
1972 (35S) bound it at <= 0.014 (``dwivedi1972_extraction.md`` sec. 4 item 4). Nothing is measured above 121 C. A
fitted share(pH, T) is not identifiable (the premortem of 2026-11-17), so the cell waits for data, and until then an
output that depends on the zero there is not shown at the zero.

THE RULE: a sensitivity rival, decided once per (pot, targets) at the shipped centre, as S3's (``ph_balance_refusal``).
On a sulfur-lane pot whose charge holds thiamine and whose program goes above
``ThiamineH2sUnmeasuredRival.above_c`` (105 C), the pot is integrated a second time with the share set to the carried
rival (``parameters_sulfur.THIAMINE_H2S_UNMEASURED_RIVAL``: 1971's ~0.36) wherever the share is evaluated below pH
6.25 above 105 C -- read from the program's T(t) and the pH trajectory the pot actually runs (solved, held or
prescribed; the share's own label-scale table), inside ``sulfur.sulfur_rate_constants_at``. Where the pot never enters
that region the second integration is the first, bit for bit, and nothing is refused. Each requested output whose
log10 spread between the shipped point and the rival exceeds ``threshold_dex()`` is REFUSED in that pot: dropped
from the answer, level and ratio alike, with a reason naming the two sources and what would lift it; a request
whose every output is refused is ``out_of_envelope``.

ONE THRESHOLD, PER ARM, AND WHY (the simplest honest rule). S3 judges a level by tau (0.964 dex) and a ratio by
omega (0.419 dex) over both arms' point values, which needs a ratio hook in every surface that forms a ratio
(``compare``, the panel's within-study rows, the contrast diagnostic, the null baselines, the directional claims).
This rule refuses in the arm instead, so every surface inherits it from ``engine.predict`` and none is touched. An arm
judged alone has to be safe for any partner: a ratio's shift is the difference of its arms' shifts, so it stays
within omega only when each arm moves by at most omega / 2. ``threshold_dex()`` is therefore omega / 2 (0.21 dex),
applied to levels too. Its cost is over-refusal, stated and counted rather than hidden: a level that moves between
omega / 2 and tau is refused though S3's R1 would show it, and a ratio of two arms that move together is refused
though it cancels. What it meets (2026-11-17): on the P3 design's corners and the panel's thiamine pot the spreads
fall in two groups, <= 0.03 dex and >= 0.52 dex, so there omega / 2 and omega refuse the same outputs; in the
directional panel Cerny 2007's pH-5 thiamine arm moves furfural by 0.26 dex, refused here and answered under omega
(the claim is not evaluable either way: S3 refuses its steps). Levels between omega and tau are refused though S3's R1
would show them (furfural 0.52 - 0.81 dex in thiamine pots at pH 6; Bolton 1994's FFT and H2S, 0.61 / 0.64 dex).
That over-refusal is the price of the one-hook rule.

WHAT IS NOT REFUSED. A thiamine-free pot, a pot that stays at or below 105 C, a pot whose share-reading pH never
falls below 6.25 above 105 C, and every lane but the sulfur lane: the shipped numbers there are bit-identical. The
pH the balance solves is not judged by this rule (S3 judges it). The rival runs only in the NumPy engine, never in JAX (its config key is held in the port).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Sequence, Tuple

#: THE SWITCH (config key ``sulfur.thiamine_h2s_refusal``): ``engine.predict``'s default. True ships; False only where
#: the model's own numbers are computed rather than shown (the rule's own integrations, a fit's objective).
ENABLED = True

#: The prefix every refusal by this rule starts with (the scorers match on it).
THIAMINE_REFUSED_TAG = "REFUSED, THIAMINE H2S SHARE UNMEASURED (pH < 6.25 above 105 C; Dwivedi & Arnold 1971 vs 1972)"
#: The species key the rule needs in the charge.
THIAMINE_KEY = "THI"

LIFT_TEXT = (
    "WHAT WOULD LIFT IT: a measurement of thiamine's H2S (or 2-furfurylthiol) above 100 C below pH 6.25 -- the P3 "
    "factorial study's thiamine arms, or a partner laboratory's thiamine series -- that decides between Dwivedi & "
    "Arnold 1971's ~0.36 and Dwivedi & Arnold 1972's <= 0.014, or replaces both (STEP5_RECORD.md, amendment "
    "2026-11-17)."
)


def thresholds() -> Dict[str, float]:
    """S3's terms (``ph_balance_refusal.thresholds``: tau, omega) and this rule's per-arm threshold, omega / 2."""
    from . import ph_balance_refusal as s3

    th = s3.thresholds()
    return {"tau_dex": th["tau_dex"], "omega_dex": th["omega_dex"], "arm_dex": 0.5 * th["omega_dex"]}


def threshold_dex() -> float:
    """The per-arm spread above which an output is refused: omega / 2 (module docstring)."""
    return thresholds()["arm_dex"]


def _max_temperature_c(spec: Any) -> float:
    thermal = spec.process.thermal
    temps = [float(c) for _, c in thermal.segments]
    if getattr(thermal, "profile", None):
        temps += [float(c) for _, c in thermal.profile]
    return max(temps)


def applies(spec: Any, declaration: Any) -> bool:
    """Whether the rule must integrate the rival: the sulfur lane, thiamine in the charge, the program above 105 C."""
    from .engine import SULFUR
    from .parameters_sulfur import THIAMINE_H2S_UNMEASURED_RIVAL as rival

    lanes = tuple(getattr(declaration, "lanes", None) or (declaration.lane,))
    if SULFUR not in lanes:
        return False
    if float(dict(declaration.mapped_precursors or {}).get(THIAMINE_KEY, 0.0)) <= 0.0:
        return False
    return _max_temperature_c(spec) > float(rival.above_c)


@dataclass
class Verdict:
    """The rule's reading of one (pot, targets) at the shipped centre."""

    #: {compound: log10 spread between the shipped point and the rival} (inf: a zero beside a number)
    spread_dex: Dict[str, float] = field(default_factory=dict)
    refused: Tuple[str, ...] = ()
    failure: Optional[str] = None

    def summary(self) -> Dict[str, Any]:
        """What a run's metadata may carry: spreads and verdicts, no level."""
        from .parameters_sulfur import THIAMINE_H2S_UNMEASURED_RIVAL as rival

        return {"rule": "thiamine H2S share, carried rival where the zero is unmeasured (STEP5_RECORD amendment "
                        "2026-11-17)",
                "rival_share": float(rival.share), "above_c": float(rival.above_c),
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

    rival = replace(ps.THIAMINE_H2S_UNMEASURED_RIVAL, enabled=True)
    return config.Config.from_globals().replace({"sulfur.thiamine_h2s_unmeasured_rival": rival})


def assess(spec: Any, targets: Sequence[str], shipped_run: Any = None) -> Verdict:
    """The verdict for ``spec`` asked for ``targets``: the shipped centre (``shipped_run`` when the caller already
    integrated it, else integrated here) against the rival. Cached by (spec, targets)."""
    from . import config, engine
    from .ph_balance_refusal import log_spread

    key = _key(spec, targets)
    if key in _CACHE:
        return _CACHE[key]
    if shipped_run is None:
        shipped_run = engine.predict(spec, list(targets), balance_refusal=False, thiamine_refusal=False,
                                     cysteine_thiamine_refusal=False)
    shipped = {c: float(v) for c, v in shipped_run.concentrations_ug_per_l.items()}
    verdict = Verdict()
    try:
        with config.applied(_rival_config()):
            rival_run = engine.predict(spec, list(targets), balance_refusal=False, thiamine_refusal=False,
                                       cysteine_thiamine_refusal=False)
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
    """The reason for refused outputs: the region, the two sources, the spread, and what would lift it."""
    from .parameters_sulfur import THIAMINE_H2S_UNMEASURED_RIVAL as rival

    named = ""
    if compounds:
        parts = []
        for c in compounds:
            s = verdict.spread_dex.get(c) if verdict is not None else None
            parts.append(f"{c} (spread {'inf' if s == math.inf else f'{s:.2f}'} dex)" if s is not None else str(c))
        named = ", ".join(parts) + ": "
    failure = (f" The rival could not be run on this pot ({verdict.failure}), so the output's independence from the "
               "share is not shown." if verdict is not None and verdict.failure else "")
    return (
        f"{THIAMINE_REFUSED_TAG}: {named}this output depends on thiamine's H2S share in a region where the shipped "
        f"share is zero (below pH 6.25, above {rival.above_c:g} C) and that zero is not measured: it is measured "
        "only near 100 C (Dwivedi & Arnold 1971); at 120 C, pH 6.0, Dwivedi & Arnold 1971 (methylene blue) read "
        f"~{rival.share:g} of the degraded thiamine as H2S and Dwivedi & Arnold 1972 (35S, 121 C) at most 0.014; "
        "nothing is measured above 121 C. With the share set to 1971's reading there, the output moves by more than "
        f"{threshold_dex():.2f} dex (omega / 2, so that no ratio it enters can move beyond omega), so it is not shown "
        "at the zero." + failure + " " + LIFT_TEXT
    )


__all__ = ["ENABLED", "LIFT_TEXT", "THIAMINE_KEY", "THIAMINE_REFUSED_TAG", "Verdict", "applies", "assess",
           "cached_verdict", "refusal_text", "threshold_dex", "thresholds"]
