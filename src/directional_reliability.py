"""Per-axis reliability tags, read from the core's directional scorecard.

WHY THIS MODULE EXISTS
----------------------
The measured skill of this model is **ordinal**, not absolute: its absolute concentrations
are mostly wrong out of sample, while its direction sense on a few axes is measured to beat
a coin. A user-facing interface that leads with a ppb number is therefore reporting the part
of the model that was measured to be wrong. The comparative CLI leads with ratios instead --
and a ratio is only as good as the model's direction sense **on the axis the two arms differ
along**, which the directional panel measures per category.

So this module does exactly one thing: it turns ``results/validation/core_directional_scores.json``
(the artifact ``kinetic_core.directional.score_panel`` writes; step B7, 2026-09-03) into
``(agree, evaluable, verdict)`` triples per axis, and works out which axes a given pair of
formulation specs actually exercises. No count is typed here: when the panel is re-scored the
artifact moves and every tag in the CLI and in the README model card moves with it.

**Nothing here is a measurement.** It reads one and applies a threshold.

THE THRESHOLDS, AND WHY THEY SIT WHERE THEY DO
----------------------------------------------
``trust``       >= 0.80 agreement on >= 3 evaluable claims AND the 95 % Wilson lower bound of
                the rate above 0.50 (demonstrably better than a coin; added 2026-09-03).
``caution``     >= 0.60 agreement, or a rate at/above 0.80 on fewer than 3 claims (an axis
                measured twice and right twice is encouraging, not established).
``do-not-use``  < 0.60 agreement.

The 0.60 floor is not arbitrary and it is not tuned: most panel rows are binary orderings, so
a coin scores ~0.50. An axis that cannot clear 0.60 has not been shown to beat guessing on
this evidence. The thresholds live here rather than in the artifact so that the artifact
stays a record of measurement only; they have not been moved since 2026-08-28. A verdict that
changes because a refit moved the counts (B9 took pH on the sulfur lane to 4/5, ``trust`` at
the boundary) is recorded in the tests as a re-pin, not hidden by moving a threshold.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

from src import data_paths

#: The core's directional scorecard (step B7, 2026-09-03). The retired lane's markdown
#: report at ``docs/validation/directional_accuracy_report.md`` is history and is not read.
DIRECTIONAL_SCORES_PATH = data_paths.CORE_DIRECTIONAL_SCORES
DIRECTIONAL_PANEL_PATH = data_paths.DIRECTIONAL_CLAIMS_PANEL

TRUST_MIN_RATE = 0.80
CAUTION_MIN_RATE = 0.60
MIN_EVALUABLE_FOR_TRUST = 3

VERDICT_TRUST = "trust"
VERDICT_CAUTION = "caution"
VERDICT_DO_NOT_USE = "do-not-use"


@dataclass(frozen=True)
class AxisReliability:
    """What the directional panel measured about one comparison axis."""

    axis: str
    agree: int
    evaluable: int
    verdict: str
    note: str = ""

    @property
    def rate(self) -> float:
        return 0.0 if self.evaluable <= 0 else self.agree / self.evaluable

    @property
    def score(self) -> str:
        return f"{self.agree}/{self.evaluable}"

    def render(self) -> str:
        return f"{self.verdict} ({self.score})"


#: Axes the panel measures but which no pair of formulation specs can "differ along" in the
#: mechanical sense -- they are properties of the claim, not of the inputs. Kept out of the
#: spec-diffing logic and reported separately where relevant.
_NON_SPEC_AXES = frozenset({"scope", "ranking", "process_heating"})

#: Compound-class axes. A comparison of two systems always carries these when the compound
#: being compared belongs to the class, regardless of what the two specs differ in.
LIPID_LANE_COMPOUNDS = frozenset({"hexanal", "nonanal", "2-pentylfuran", "1-hexanol", "octanal"})

#: The sulfur branch. Wave S2c established it has ZERO absolute literature anchors -- both
#: values on its only benchmark are marked ``no_verifiable_source``. That is a statement about
#: absolute scale, not about direction, and the two must not be conflated: the same branch
#: predicted an unseen isotope-dilution MFT measurement to 1.52x (Wave U).
SULFUR_COMPOUNDS = frozenset(
    {
        "furanthiol",
        "furfurylthiol",
        "furyl) disulfide",
        "hydrogen sulfide",
        "methanethiol",
        "thiophene",
        "thiazol",
        "thiol",
    }
)


def load_directional_artifact(path: Path | str = DIRECTIONAL_SCORES_PATH) -> Dict[str, Any]:
    """The core's directional scorecard (``results/validation/core_directional_scores.json``).

    Raises rather than falling back: a silent default here would republish a number after
    the evidence behind it stopped being read. (Until 2026-09-03 this module parsed the
    CURRENT STANDING table of the retired lane's markdown report; that report is history.)
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"core directional scorecard not found at {path}. Regenerate it with "
            "scripts/generators/generate_core_directional_scores.py; there is no fallback copy."
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path}: not a JSON artifact ({exc})") from exc
    if not isinstance(payload, dict) or payload.get("artifact") != "core_directional_scores":
        raise ValueError(f"{path}: not a core_directional_scores artifact")
    return payload


def load_panel_counts(path: Path | str = DIRECTIONAL_SCORES_PATH) -> Dict[str, tuple]:
    """``{category: (agree, evaluable)}`` over the STRICTLY INDEPENDENT claims."""
    payload = load_directional_artifact(path)
    table = payload["summary"]["independent"]["by_category"]
    if not table:
        raise ValueError(f"{path}: the independent split has no categories")
    return {str(k): (int(v["agree"]), int(v["evaluable"])) for k, v in table.items()}


def load_axis_notes(path: Path | str = DIRECTIONAL_SCORES_PATH) -> Dict[str, str]:
    """One sentence per category, derived from the artifact rather than typed."""
    payload = load_directional_artifact(path)
    notes: Dict[str, str] = {}
    for axis, bucket in payload["summary"]["independent"]["by_category"].items():
        parts = []
        if bucket["misses"]:
            parts.append(f"misses: {', '.join(bucket['misses'])}")
        if bucket.get("mechanism_absent"):
            parts.append(
                f"{bucket['mechanism_absent']} of them are identical predictions across the moved "
                "axis (the lane carries no term for it)"
            )
        if bucket["not_evaluable"]:
            parts.append(f"{bucket['not_evaluable']} claim(s) not evaluable on the core")
        notes[str(axis)] = "; ".join(parts) if parts else "every evaluable independent claim agrees"
    return notes


def verdict_for(agree: int, evaluable: int) -> str:
    """Apply the documented thresholds. Pure arithmetic, no tuning."""
    if evaluable <= 0:
        return VERDICT_DO_NOT_USE
    rate = agree / evaluable
    if rate < CAUTION_MIN_RATE:
        return VERDICT_DO_NOT_USE
    if rate >= TRUST_MIN_RATE and evaluable >= MIN_EVALUABLE_FOR_TRUST and wilson_lower(agree, evaluable) > COIN:
        return VERDICT_TRUST
    return VERDICT_CAUTION


#: 2026-09-03 (B9 took pH on the sulfur lane to 4/5 = 0.80 exactly). A point rate on five claims
#: cannot distinguish the model from a coin: the 95 % Wilson lower bound of 4/5 is 0.38. ``trust``
#: now also requires that lower bound to clear 0.50, i.e. the axis is demonstrably better than a
#: coin at 95 %. 8/8 clears it (0.68), 7/7 does (0.65), 4/5 and 5/6 do not. This ADDS a condition
#: to ``trust``; the 0.80 / 0.60 thresholds themselves were not moved.
COIN = 0.50


def wilson_lower(agree: int, evaluable: int, z: float = 1.959964) -> float:
    """95 % Wilson score interval, lower bound; 0.0 when nothing was evaluable."""
    if evaluable <= 0:
        return 0.0
    p = agree / evaluable
    denom = 1.0 + z * z / evaluable
    centre = p + z * z / (2.0 * evaluable)
    half = z * math.sqrt(p * (1.0 - p) / evaluable + z * z / (4.0 * evaluable * evaluable))
    return (centre - half) / denom




def reliability_for_axis(
    axis: str,
    counts: Optional[Mapping[str, tuple]] = None,
    notes: Optional[Mapping[str, str]] = None,
) -> AxisReliability:
    table = dict(counts if counts is not None else load_panel_counts())
    if notes is None:
        try:
            notes = load_axis_notes()
        except (FileNotFoundError, ValueError):
            notes = {}
    if axis not in table:
        return AxisReliability(
            axis=axis,
            agree=0,
            evaluable=0,
            verdict=VERDICT_DO_NOT_USE,
            note=(
                "This axis is not measured by the directional panel at all. Unmeasured is "
                "reported as do-not-use deliberately: absence of evidence is not evidence."
            ),
        )
    agree, evaluable = table[axis]
    return AxisReliability(
        axis=axis,
        agree=agree,
        evaluable=evaluable,
        verdict=verdict_for(agree, evaluable),
        note=dict(notes).get(axis, ""),
    )


# ---------------------------------------------------------------------------------------
# Which axes does a given comparison actually exercise?
# ---------------------------------------------------------------------------------------

_SUGAR_TOKENS = ("ribose", "glucose", "fructose", "xylose", "arabinose", "maltose", "lactose", "sucrose")
_CYSTEINE_TOKENS = ("cysteine", "cystine", "glutathione")


def _precursor_set(spec: Mapping[str, Any]) -> set:
    return {str(name).strip().lower() for name in (spec.get("precursors") or {})}


def _sugars_in(spec: Mapping[str, Any]) -> set:
    return {name for name in _precursor_set(spec) if any(tok in name for tok in _SUGAR_TOKENS)}


def _num(spec: Mapping[str, Any], key: str) -> Optional[float]:
    value = spec.get(key)
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _differs(a: Optional[float], b: Optional[float], tolerance: float = 1e-9) -> bool:
    if a is None or b is None:
        return False
    return abs(a - b) > tolerance


def axes_exercised(spec_a: Mapping[str, Any], spec_b: Mapping[str, Any]) -> List[str]:
    """The panel categories along which two comparison arms differ.

    Deliberately over-inclusive: if a comparison moves three knobs at once, all three axes are
    reported, and the weakest verdict among them governs. A user who changes pH *and* sugar
    gets the pH tag, because the pH miss does not go away by being bundled.
    """
    axes: List[str] = []

    if _sugars_in(spec_a) != _sugars_in(spec_b):
        axes.append("sugar_identity")

    cys_a = any(tok in name for name in _precursor_set(spec_a) for tok in _CYSTEINE_TOKENS)
    cys_b = any(tok in name for name in _precursor_set(spec_b) for tok in _CYSTEINE_TOKENS)
    if cys_a != cys_b:
        axes.append("additive_cysteine")

    if _differs(_num(spec_a, "temp_C"), _num(spec_b, "temp_C")):
        axes.append("temperature")
    if _differs(_num(spec_a, "time_min"), _num(spec_b, "time_min")):
        axes.append("time")
    if _differs(_num(spec_a, "ph"), _num(spec_b, "ph")):
        axes.append("ph")
    if _differs(_num(spec_a, "aw"), _num(spec_b, "aw")):
        axes.append("moisture_aw")

    matrix_a = str(spec_a.get("protein_type", "free") or "free").lower()
    matrix_b = str(spec_b.get("protein_type", "free") or "free").lower()
    if matrix_a != matrix_b:
        axes.append("matrix_identity")

    return axes


def _matches(compound: str, vocabulary: frozenset) -> bool:
    """Substring match against a compound vocabulary.

    Substring rather than equality because the pipeline carries three aliases for every
    compound -- canonical SMILES, species name and a display label like
    "2-Methyl-3-furanthiol (MFT)" -- and a caveat that only fires on one spelling is a caveat
    that silently does not fire.
    """
    name = str(compound).strip().lower()
    return any(token in name for token in vocabulary)


def compound_axes(compound: str) -> List[str]:
    """Compound-class axes that apply however the two arms differ."""
    axes: List[str] = []
    if _matches(compound, LIPID_LANE_COMPOUNDS):
        axes.append("lipid_lane")
    return axes


def is_sulfur_compound(compound: str) -> bool:
    return _matches(compound, SULFUR_COMPOUNDS)


def weakest(reliabilities: Sequence[AxisReliability]) -> Optional[AxisReliability]:
    """The governing tag for a multi-axis comparison: the worst one."""
    if not reliabilities:
        return None
    order = {VERDICT_DO_NOT_USE: 0, VERDICT_CAUTION: 1, VERDICT_TRUST: 2}
    return min(reliabilities, key=lambda item: (order.get(item.verdict, 0), item.rate))


def describe_comparison(
    spec_a: Mapping[str, Any],
    spec_b: Mapping[str, Any],
    *,
    counts: Optional[Mapping[str, tuple]] = None,
) -> Dict[str, Any]:
    """Full reliability picture for one A-vs-B comparison."""
    table = dict(counts if counts is not None else load_panel_counts())
    try:
        notes = load_axis_notes() if counts is None else {}
    except (FileNotFoundError, ValueError):
        notes = {}
    axes = axes_exercised(spec_a, spec_b)
    reliabilities = [reliability_for_axis(axis, table, notes) for axis in axes]
    governing = weakest(reliabilities)
    return {
        "axes": axes,
        "per_axis": reliabilities,
        "governing": governing,
        "no_axis_differs": not axes,
    }
