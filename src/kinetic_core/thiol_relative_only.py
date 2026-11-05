"""
The thiol lane ships RELATIVE within-laboratory predictions only (H0 / EX-059, D2 branch (i); applied 2026-11-05).

H0 (``results/validation/kinetic_core_thiol_lolo_nested_tau.{json,md}``; pre-registration
``kinetic_core_between_lab_tau_h0_prereg.md`` v3, amendments A1 and A2, review ``..._review.md``) scored the thiol
lane's admitted structure (R3En) by strict leave-two-laboratories-out nesting. Its ship rule D2 read
tau_ABS(all) = 0.964 dex (the larger of its two sets), with a one-sided 95 % upper bound U >= 1.60 (the grid top),
so branch (i) fired: the lane ships within-study ratios, ranks and time-course shapes with omega in their
intervals, and absolute concentrations for a new laboratory are refused. The owner approved applying it to the
shipped model on 2026-11-05.

R3En IS NOT THE SHIPPED VECTOR. H0 scored R3En (the joint deoxyosone + thiol fit, s042_r3en_fixed); the shipped
sulfur constants are B60's. Refusing the shipped lane's absolute levels is the conservative transfer (nothing says
B60 transfers between laboratories better than R3En did; the review put P(tau_ABS > 0.8 for the shipped model)
at ~0.75), and the omega in the ratio intervals is BORROWED from R3En's scoring.

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

The ruling is applied to every one of those observables that the SHIPPED sulfur lane answers for a pentose +
cysteine pot (re-derived through declare_envelope by tests/unit/test_thiol_relative_only.py): MFT, FFT, MFTD,
FUR, NF and MGO. MGO counts only when the sulfur lane runs (it is most of tau_ABS's rows: Jiangnan 81,
Ochanomizu 3); a trunk MGO, from a hexose and an amine with no sulfur charge, is another lane's network, which H0
did not score, and stays answered. TDP (3-deoxypentosone) is not a target. NOT in scope, because H0 did not score
them: ACTZ (2-acetylthiazole: one Jiangnan row, unpaired -- no null answers it -- so it entered only the D-all
diagnostic), FFTD and DPO (unpaired; neither is a target), Nestle's acetic acid and pentose rows (unpaired, D-all
only), H2S and methanethiol (not in the fit table). Norfuraneol and furfural are in: H0 scored both in the absolute
class (they are R3En's own products, not reactants). MFTD is in although H0 scored it on shape only: D2 rules the
lane's observables, and an observable H0 tested only relatively has no absolute evidence at all.

A NUMBER IS STILL COMPUTED. The engine integrates these species as before (the between-condition ratios need
them) and names them in the declaration's FIELD ``EnvelopeDeclaration.absolute_refused_targets`` (not a warning:
the envelope state is unchanged). Every output that would emit an absolute level -- ``predict``, the per-arm
absolutes and OAV tables of ``compare``, ``score``, ``calibrate``'s shipped-model hold-out rows, the panel
scorecard and the envelope -- refuses it by name with :func:`refusal_text`. The one exception is a per-laboratory
calibration made by ``calibrate`` (:func:`calibration_rescues`): the level is then that laboratory's own, and the
interval carries the shape error of a two-row ratio (:func:`calibrated_extra_decades`).

THE RATIO INTERVAL. No interval in the engine carried a within-laboratory shape term before this ruling (no ratio
had an interval at all; the dispersion band ``reliability_band_x`` is the measurement's same-sample scatter). H0's
relative predictive gives a row centred on its group the variance s_i^2 - 1/W + c_i * omega^2 (prereg sec. 2 and
sec. 5, score 6): a two-row, one-species ratio carries 2 * omega^2 + s_a^2 + s_b^2, so the 90 % half-width is
z90 * sqrt(2 * omega^2 + 2 * s^2). s is the smallest class noise of EX-056's largest training set (the Reading fold,
513 rows: calibrated 0.097 dex), because a ratio from ``compare`` names no measurement method. omega is the LARGER of
H0's two sets' omega(all) (0.419 dex, warm; reset 0.390): the ship rule reads the larger set (prereg sec. 3.6 and
D2; amendment A1 item 3), and the review notes the declared leaks bias omega low (choices register CR-19).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Mapping, Optional, Tuple

from src import data_access, data_paths

#: The frozen H0 scoring record (EX-059). Read, never rewritten (scripts/ci/artifact_freshness_gate.py pins it).
H0_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_thiol_lolo_nested_tau.json"
#: The EX-056 fold record whose class noise the ratio interval reads (the largest training set).
NOISE_RECORD = data_paths.VALIDATION_DIR / "kinetic_core_thiol_lolo_reading_mottram.json"

#: The target species the ruling covers on the SULFUR lane (module docstring: the table, and why each is in or out).
RELATIVE_ONLY_SPECIES = frozenset({"MFT", "FFT", "MFTD", "FUR", "NF", "MGO"})

#: The prefix every refusal of an absolute level starts with (the scorers match it as they match
#: ``HEXOSE ENTRY UNIDENTIFIED``).
ABSOLUTE_REFUSED_TAG = "ABSOLUTE LEVEL REFUSED (H0 / EX-059)"

#: Two-sided 90 % normal quantile (statistics.NormalDist().inv_cdf(0.95)); the envelope's and the
#: calibration's 90 % convention.
Z90 = 1.6448536269514722

#: What ``user_fit.calibrate`` writes as its calibration's provenance; a calibration without it was not made by it.
CALIBRATE_GENERATOR = "src/kinetic_core/user_fit.py"


@dataclass(frozen=True)
class H0Ruling:
    """The numbers the ruling quotes, read from the frozen records."""

    tau_abs_dex: float           # tau_ABS(all) the ship rule read (the larger set)
    tau_abs_set: str
    upper_bound_text: str        # "U >= 1.60" (the grid top)
    upper_bound_dex: float
    omega_dex: float             # the larger of the two sets' omega(all)
    omega_set: str
    class_noise_dex: float       # the smallest class noise of the largest EX-056 training set
    branch: str
    verdict: str
    relative_branch: Tuple[Tuple[str, float, str], ...]   # (set, engine - null nats/row, label)


@lru_cache(maxsize=1)
def ruling() -> H0Ruling:
    """H0's D2 numbers from its record. Raises if the record no longer says branch (i): the ruling would then
    have no basis, and silently keeping it would be worse than failing."""
    record = data_access.load_json(H0_RECORD)
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
    noise = data_access.load_json(NOISE_RECORD)["models"]["engine_warm"]["terms"]["class_noise"]
    return H0Ruling(
        tau_abs_dex=tau, tau_abs_set=tau_set, upper_bound_text=str(d2["U_text"]),
        upper_bound_dex=float(d2["U_used"]), omega_dex=float(per_set[omega_set]["omega_all"]),
        omega_set=omega_set, class_noise_dex=float(min(noise.values())), branch=str(d2["branch"]),
        verdict=str(record["verdict"]["reported"]), relative_branch=relative,
    )


def in_scope(species_key: Optional[str]) -> bool:
    """On the sulfur lane (the caller gates on the lane)."""
    return species_key in RELATIVE_ONLY_SPECIES


def _relative_labels() -> str:
    return "; ".join(f"{s} {m:+.2f} nats/row ({label})" for s, m, label in ruling().relative_branch)


def refusal_text() -> str:
    """The named refusal of an absolute level: what decided it, and what would lift it."""
    r = ruling()
    return (
        f"{ABSOLUTE_REFUSED_TAG}: the thiol lane ships relative within-laboratory predictions only. H0 scored the lane "
        "by leave-two-laboratories-out and its pre-registered ship rule D2 took branch (i): the level of a lane "
        f"observable in a laboratory the model has never seen is off by tau_ABS = {r.tau_abs_dex:.3f} dex (RMS), "
        f"with a one-sided 95 % upper bound {r.upper_bound_text} dex (the test grid's top); the absolute evidence for "
        "the thiols comes from one laboratory (DFA Garching), and nominal 90 % intervals covered about 75 - 80 % out "
        "of sample. H0 scored the admitted structure R3En; the shipped sulfur constants are B60's, so refusing their "
        "absolute levels is the conservative transfer. The concentration is computed and refused; between-condition "
        "ratios, ranks and time-course shapes from `compare` are answered, with the within-laboratory shape error "
        f"omega = {r.omega_dex:.2f} dex per arm (borrowed from R3En's scoring) in their interval. WHAT WOULD LIFT "
        "IT: absolute measurements of these compounds from a second independent laboratory (for example the "
        "pending WUR blind-test samples), re-scored under H0's rule (H0', "
        "results/validation/kinetic_core_between_lab_tau_h0_prereg.md). A calibration made by `calibrate` from "
        "your own measured levels of the compound is answered: its level is then your laboratory's."
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
    return (
        "RELATIVE ONLY; NOT SHOWN BETTER THAN A NULL, AND IN H0's RESET SET THE NULL IS BETTER (H0 / EX-059): "
        "within-study ratios, ranks and shapes are what the thiol lane ships, but H0's relative-branch score, engine "
        "minus a chemistry-free null (the two nulls are identical on this branch), was "
        f"{_relative_labels()} ('the null is better' is D3's label: below -2 SE with sign-flip p <= 2/32). These "
        "relative predictions have no evidence of beating a null either. The 90 % interval carries "
        f"omega = {r.omega_dex:.2f} dex in each arm and the measurement noise s = {r.class_noise_dex:.2f} dex of each "
        f"(a two-arm ratio: a factor of {ratio_interval_factor():.1f} either way); it protects only a ratio of ONE "
        "compound between two conditions in one study, not a ratio across compounds."
    )


def ratio_half_width_dex() -> float:
    """90 % half-width of a two-arm, one-species ratio: z90 * sqrt(2 omega^2 + s_a^2 + s_b^2)."""
    r = ruling()
    return Z90 * math.sqrt(2.0 * r.omega_dex ** 2 + 2.0 * r.class_noise_dex ** 2)


def ratio_interval_factor() -> float:
    return 10.0 ** ratio_half_width_dex()


def ratio_interval(ratio: Any) -> Optional[Tuple[float, float]]:
    """The 90 % interval of a positive finite ratio; ``None`` for an undefined ratio."""
    if not isinstance(ratio, (int, float)) or not math.isfinite(float(ratio)) or float(ratio) <= 0.0:
        return None
    f = ratio_interval_factor()
    return (float(ratio) / f, float(ratio) * f)


def calibrated_extra_decades(rf_sigma_log10: float, n_rows: int) -> float:
    """The 90 % half-width of a calibrated level of an in-scope compound. The level is itself a two-row ratio -- the
    new condition against the mean of the n fit rows the response factor came from -- so it carries the factor's
    own sigma, the new row's omega and the anchor mean's omega / sqrt(n), in quadrature."""
    omega = ruling().omega_dex
    return Z90 * math.sqrt(float(rf_sigma_log10) ** 2 + omega ** 2 + omega ** 2 / max(1, int(n_rows)))


def declared(declaration: Any, compound: str) -> bool:
    """Whether ``compound`` (a requested name) is declared absolute-refused on this declaration (a dataclass or
    its ``as_dict``)."""
    names = (declaration.get("absolute_refused_targets") if isinstance(declaration, Mapping)
             else getattr(declaration, "absolute_refused_targets", ()))
    return str(compound) in tuple(names or ())


def calibration_rescues(calibration: Any, compound: str) -> bool:
    """A calibration rescues a compound's level only when ``calibrate`` made it from data: a response factor
    with at least one row and a positive sigma, fit records named, and calibrate's provenance. A hand-written or
    empty file does not (review M2, 2026-11-05)."""
    if calibration is None:
        return False
    rf = (getattr(calibration, "response_factors", None) or {}).get(str(compound))
    if rf is None or int(getattr(rf, "n_rows", 0)) < 1 or not float(getattr(rf, "sigma", 0.0)) > 0.0:
        return False
    if not tuple(getattr(calibration, "fit_records", ()) or ()):
        return False
    provenance = getattr(calibration, "provenance", None) or {}
    return str(provenance.get("generated_by", "")) == CALIBRATE_GENERATOR


def refused_for(declaration: Any, compound: str, calibration: Any = None) -> bool:
    """Declared, and not rescued by a calibration :func:`calibration_rescues` accepts."""
    return declared(declaration, compound) and not calibration_rescues(calibration, compound)


__all__ = [
    "ABSOLUTE_REFUSED_TAG", "CALIBRATE_GENERATOR", "H0_RECORD", "H0Ruling", "NOISE_RECORD", "RELATIVE_ONLY_SPECIES",
    "Z90", "calibrated_extra_decades", "calibration_rescues", "declared", "in_scope", "ratio_half_width_dex",
    "ratio_interval", "ratio_interval_factor", "refusal_text", "refused_for", "relative_caveat", "ruling",
    "short_reason",
]
