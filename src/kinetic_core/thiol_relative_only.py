"""
The thiol lane ships RELATIVE within-laboratory predictions only (H0 / EX-059, D2 branch (i); applied 2026-11-05).

H0 (``results/validation/kinetic_core_thiol_lolo_nested_tau.{json,md}``; pre-registration
``kinetic_core_between_lab_tau_h0_prereg.md`` v3, amendments A1 and A2, review ``..._review.md``) scored the thiol
lane's admitted structure (R3En) by strict leave-two-laboratories-out nesting. Its ship rule D2 read
tau_ABS(all) = 0.964 dex (the larger of its two sets), with a one-sided 95 % upper bound U >= 1.60 (the grid top),
so branch (i) fired: the lane ships within-study ratios, ranks and time-course shapes with omega in their
intervals, and absolute concentrations for a new laboratory are refused. The owner approved applying it to the
shipped model on 2026-11-05.

WHICH OBSERVABLES (read from H0's own rows, not from memory). The ship-rule quantities are computed on EX-056's
PAIRED rows of the five scorable laboratories (score_thiol_lolo_nested_tau.py choice 19; the ``paired`` flag of
each ``kinetic_core_thiol_lolo_<lab>.json`` row). Those rows hold these species of the lane (laboratory: rows):

  =========  =========================================  ======================  ==========================
  species    ABS (sida, calibrated; the D2 quantity)    FIG (semi_figure)       FREE (zhou_free, shape)
  =========  =========================================  ======================  ==========================
  MFT        DFA 7 (sida)                               --                      Dalian 4, Jiangnan 33
  FFT        DFA 7 (sida)                               --                      Dalian 4, Jiangnan 34
  MFTD       --                                         --                      Dalian 4, Jiangnan 4
  FUR        DFA 1 (calibrated)                         Dalian 4                Jiangnan 33
  NF         DFA 1, Ochanomizu 2 (calibrated)           Reading 4               Jiangnan 1
  MGO        Jiangnan 81, Ochanomizu 3 (calibrated)     --                      Jiangnan 32
  TDP        --                                         Jiangnan 123, Reading 4 --
  =========  =========================================  ======================  ==========================

The ruling is applied to every one of those observables that the SHIPPED sulfur lane emits as a target:
MFT, FFT, MFTD, FUR, NF. The rest need nothing: MGO is already refused on the sulfur lane (B13: the dicarbonyl
steps are trunk-only, so a sulfur-lane MGO is never emitted) and the trunk's MGO is another lane's network, which
H0 did not score; TDP (3-deoxypentosone) is not a target. NOT in scope, because H0 did not score them: ACTZ
(2-acetylthiazole: one Jiangnan row, unpaired -- no null answers it -- so it entered only the D-all diagnostic),
FFTD and DPO (unpaired; neither is a target), cysteine, pentose and acetic acid (Nestle's rows: unpaired, D-all
only), H2S and methanethiol (not in the fit table). Norfuraneol and furfural are in: H0 scored both in the absolute
class (they are R3En's own products, not reactants). MFTD is in although H0 scored it on shape only: D2 rules the
lane's observables, and an observable H0 tested only relatively has no absolute evidence at all.

A NUMBER IS STILL COMPUTED. The engine integrates these species as before (the between-condition ratios need
them), declares them in ``EnvelopeDeclaration.absolute_refused_targets`` with a warning carrying
:data:`ABSOLUTE_REFUSED_TAG`, and every output that would emit an absolute level -- ``predict``, the per-arm
absolutes and OAV tables of ``compare``, the ``score`` verb, the panel scorecard and the envelope -- refuses it by
name with :func:`refusal_text`. The one exception is a per-laboratory calibration that carries a response factor
for the compound: its level is then set by that laboratory's own measurements, which is the within-laboratory case
the ruling allows, and omega is added to its interval (calibration.Calibration.apply_factors).

OMEGA ENTERS ONCE. No interval in the engine carried a within-laboratory shape term before this ruling (no ratio
had an interval at all; the dispersion band ``reliability_band_x`` is the measurement's same-sample scatter and
decides only "resolved"). A between-condition ratio of one species is two rows of one level group, so it carries
omega in each arm: half-width z90 * sqrt(2) * omega. omega is the LARGER of H0's two sets' omega(all) (0.419 dex,
warm; reset 0.390): the review notes the declared leaks bias omega low, and the ship rule reads the larger tau for
the same reason (choices register CR-19).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Mapping, Optional, Tuple

from src import data_paths

#: The frozen H0 scoring record (EX-059). Read, never rewritten (scripts/ci/artifact_freshness_gate.py pins it).
H0_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_thiol_lolo_nested_tau.json"

#: The sulfur-lane target species the ruling covers (module docstring: the table, and why each is in or out).
RELATIVE_ONLY_SPECIES = frozenset({"MFT", "FFT", "MFTD", "FUR", "NF"})

#: The warning prefix every consumer matches on (as the scorers match ``HEXOSE ENTRY UNIDENTIFIED``).
ABSOLUTE_REFUSED_TAG = "ABSOLUTE LEVEL REFUSED (H0 / EX-059)"

#: Two-sided 90 % normal quantile (statistics.NormalDist().inv_cdf(0.95)); the envelope's and the
#: calibration's 90 % convention.
Z90 = 1.6448536269514722


@dataclass(frozen=True)
class H0Ruling:
    """The numbers the ruling quotes, read from the frozen record."""

    tau_abs_dex: float           # tau_ABS(all) the ship rule read (the larger set)
    tau_abs_set: str
    upper_bound_text: str        # "U >= 1.60" (the grid top)
    upper_bound_dex: float
    omega_dex: float             # the larger of the two sets' omega(all)
    omega_set: str
    branch: str
    verdict: str
    relative_branch: Tuple[Tuple[str, float, str], ...]   # (set, engine - null nats/row, label)


@lru_cache(maxsize=1)
def ruling() -> H0Ruling:
    """H0's D2 numbers from its record. Raises if the record no longer says branch (i): the ruling would then
    have no basis, and silently keeping it would be worse than failing."""
    record = json.loads(H0_RECORD.read_text(encoding="utf-8"))
    d2 = record["D2"]
    if d2.get("branch") != "(i)":
        raise RuntimeError(f"{data_paths.rel(H0_RECORD)}: D2 branch is {d2.get('branch')!r}, not '(i)'")
    per_set = d2["per_set"]
    tau = float(d2["tau_ABS_all_used"])
    tau_set = next(s for s, v in per_set.items() if float(v["tau_ABS_all"]) == tau)
    omega_set = max(per_set, key=lambda s: float(per_set[s]["omega_all"]))
    relative = tuple(
        (s, float(v["relative_branch"]["typical_yield"]["mean"]), str(v["relative_branch"]["typical_yield"]["label"]))
        for s, v in record["sets"].items()
    )
    return H0Ruling(
        tau_abs_dex=tau, tau_abs_set=tau_set, upper_bound_text=str(d2["U_text"]),
        upper_bound_dex=float(d2["U_used"]), omega_dex=float(per_set[omega_set]["omega_all"]),
        omega_set=omega_set, branch=str(d2["branch"]), verdict=str(record["verdict"]["reported"]),
        relative_branch=relative,
    )


def in_scope(species_key: Optional[str]) -> bool:
    return species_key in RELATIVE_ONLY_SPECIES


def refusal_text() -> str:
    """The named refusal of an absolute level: what decided it, and what would lift it."""
    r = ruling()
    return (
        f"{ABSOLUTE_REFUSED_TAG}: the thiol lane ships relative within-laboratory predictions only. H0 scored the lane "
        "by leave-two-laboratories-out and its pre-registered ship rule D2 took branch (i): the level of a lane "
        f"observable in a laboratory the model has never seen is off by tau_ABS = {r.tau_abs_dex:.3f} dex (RMS), "
        f"with a one-sided 95 % upper bound {r.upper_bound_text} dex (the test grid's top); "
        "the absolute evidence for the thiols comes from one laboratory (DFA Garching), and nominal 90 % intervals "
        "covered about 75 - 80 % out of sample. The concentration is computed and refused; between-condition "
        "ratios, ranks and time-course shapes from `compare` are answered (with the within-laboratory shape error "
        f"omega = {r.omega_dex:.2f} dex per arm in their interval). WHAT WOULD LIFT IT: absolute measurements of "
        "these compounds from a second independent laboratory (for example the pending WUR blind-test samples), "
        "re-scored under H0's rule (H0', results/validation/kinetic_core_between_lab_tau_h0_prereg.md). A "
        "laboratory's own calibration with a response factor for the compound (`calibrate`) is answered: its level "
        "is then that laboratory's."
    )


def short_reason() -> str:
    """The one-line form the scorecards carry per refused row (the full text is :func:`refusal_text`)."""
    r = ruling()
    return (
        f"{ABSOLUTE_REFUSED_TAG}: the thiol lane ships relative predictions only (D2 branch (i): tau_ABS "
        f"{r.tau_abs_dex:.3f} dex, {r.upper_bound_text}; absolute thiol evidence from one laboratory). Lifted by a "
        "second independent laboratory's absolute measurements re-scored under H0's rule."
    )


def relative_caveat() -> str:
    """The owner's note, stated wherever a relative thiol-lane prediction is shown (H0 owner note)."""
    r = ruling()
    scores = "; ".join(f"{s} {m:+.2f} nats/row ({label})" for s, m, label in r.relative_branch)
    return (
        "RELATIVE ONLY, AND NOT SHOWN BETTER THAN A NULL (H0 / EX-059): within-study ratios, ranks and shapes are "
        "what the thiol lane ships, but H0's relative-branch score, engine minus a chemistry-free null (the two nulls "
        f"are identical on this branch), was {scores}: these relative predictions have no evidence of beating a "
        f"null either. The interval carries omega = {r.omega_dex:.2f} dex in each arm (a two-arm ratio: a factor of "
        f"{ratio_interval_factor():.1f} either way at 90 %); "
        "it protects only a ratio of ONE compound between two conditions in one study, not a ratio across compounds."
    )


def ratio_half_width_dex() -> float:
    """90 % half-width of a two-arm, one-species ratio: z90 * sqrt(2) * omega (omega once per arm)."""
    return Z90 * math.sqrt(2.0) * ruling().omega_dex


def ratio_interval_factor() -> float:
    return 10.0 ** ratio_half_width_dex()


def ratio_interval(ratio: float) -> Optional[Tuple[float, float]]:
    """The 90 % interval of a positive finite ratio under omega; ``None`` for an undefined ratio."""
    if not isinstance(ratio, (int, float)) or not math.isfinite(float(ratio)) or float(ratio) <= 0.0:
        return None
    f = ratio_interval_factor()
    return (float(ratio) / f, float(ratio) * f)


def calibrated_extra_decades(rf_sigma_log10: float) -> float:
    """The 90 % half-width a calibrated level of an in-scope compound carries: the response factor's own sigma and
    one row's omega, in quadrature. The level comes from the laboratory; the new condition's shape error does not."""
    return Z90 * math.hypot(float(rf_sigma_log10), ruling().omega_dex)


def declared(declaration: Any, compound: str) -> bool:
    """Whether ``compound`` (a requested name) is declared absolute-refused on this declaration."""
    return str(compound) in tuple(getattr(declaration, "absolute_refused_targets", ()) or ())


def refused_for(declaration: Any, compound: str, calibration: Any = None) -> bool:
    """Declared, and not rescued by a calibration carrying this compound's response factor."""
    if not declared(declaration, compound):
        return False
    factors: Mapping[str, Any] = getattr(calibration, "response_factors", None) or {}
    return str(compound) not in factors


__all__ = [
    "ABSOLUTE_REFUSED_TAG", "H0_RECORD", "H0Ruling", "RELATIVE_ONLY_SPECIES", "Z90", "calibrated_extra_decades",
    "declared", "in_scope", "ratio_half_width_dex", "ratio_interval", "ratio_interval_factor", "refusal_text",
    "refused_for", "relative_caveat", "ruling",
]
