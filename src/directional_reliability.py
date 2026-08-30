"""Per-axis reliability tags, read from the directional-accuracy artifact.

WHY THIS MODULE EXISTS
----------------------
2026-08-28 (Wave S5). The measured skill of this model is **ordinal**, not absolute: the
hold-out medians are 6x on the free-precursor lane and 67-94x on the matrix lane, while the
directional panel scores 24/35 on strictly independent claims (Wave W, 2026-08-28; was
21/29). A user-facing interface that
leads with a ppb number is therefore reporting the part of the model that was measured to be
wrong. The comparative CLI leads with ratios instead -- and a ratio is only as good as the
model's direction sense **on the axis the two arms differ along**, which the panel measures
per category.

So this module does exactly one thing: it turns
``docs/validation/directional_accuracy_report.md`` -- specifically its CURRENT STANDING
section, the one place in the repo that states the live per-category counts -- into
``(agree, evaluable, verdict)`` triples, and works out which categories a given pair of
formulation specs actually exercises.

**Nothing here is a measurement.** It parses one and applies a threshold. If the panel is
re-scored, edit the table in that report and every tag in the CLI and in the README model
card moves with it; there is no second copy of the numbers to forget.

THE THRESHOLDS, AND WHY THEY SIT WHERE THEY DO
----------------------------------------------
``trust``       >= 0.80 agreement on >= 3 evaluable claims.
``caution``     >= 0.60 agreement, or a rate at/above 0.80 on fewer than 3 claims (an axis
                measured twice and right twice is encouraging, not established).
``do-not-use``  < 0.60 agreement.

The 0.60 floor is not arbitrary and it is not tuned: most panel rows are binary orderings, so
a coin scores ~0.50, and the report's own §2 says the model's edge over a coin is "real but
modest". An axis that cannot clear 0.60 has not been shown to beat guessing on this evidence,
and the report's §A6 says in terms that the model "licenses no pH or moisture recommendation".
Under these thresholds ``moisture_aw`` (0/3) lands on ``do-not-use``.

UPDATED 2026-08-28 (Wave W). ``ph`` was 4/7 = 0.571 and therefore ``do-not-use``. Two new
independent pH rows from Mottram & Nobrega 2002 both AGREE, taking it to 6/9 = 0.667, which
clears the 0.60 floor and moves the tag to ``caution``. THE THRESHOLD WAS NOT MOVED TO GET
THIS -- the evidence moved. Read the change conservatively: 6/9 is two rows above the floor,
both from one source, and the report's SS8/SSA6 licence still says the model "licenses no pH
recommendation". ``caution`` here means "no longer demonstrably worse than a coin", not
"trustworthy". The thresholds live here rather than in the report so that the report stays a
record of measurement only.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]

#: The artifact this module reads. Its CURRENT STANDING section is the live panel score.
DIRECTIONAL_REPORT_PATH = ROOT / "docs" / "validation" / "directional_accuracy_report.md"

#: The panel data itself, used only to check that the categories we tag are categories the
#: panel actually contains -- so a typo in an axis name is caught rather than silently
#: reported as "unmeasured".
DIRECTIONAL_PANEL_PATH = ROOT / "docs" / "validation" / "directional_claims_panel.yml"

_CURRENT_STANDING_HEADING = "# CURRENT STANDING"

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


def _parse_count_tables(text: str) -> Dict[str, tuple]:
    """Pull every ``| label | int | int |`` row out of a markdown blob."""
    counts: Dict[str, tuple] = {}
    row = re.compile(r"^\|\s*([^|]+?)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*$")
    for line in text.splitlines():
        match = row.match(line.strip())
        if not match:
            continue
        label = match.group(1).strip().strip("*").strip()
        if not label or label.lower() in {"category", "bucket"}:
            continue
        counts[label] = (int(match.group(2)), int(match.group(3)))
    return counts


def load_panel_counts(report_path: Path | str = DIRECTIONAL_REPORT_PATH) -> Dict[str, tuple]:
    """``{label: (agree, evaluable)}`` from the report's CURRENT STANDING section.

    Raises rather than falling back to a baked-in default. A silent default here would be the
    exact failure this repository has been removing all campaign: a number that keeps being
    published after the evidence behind it stops being read.
    """
    path = Path(report_path)
    if not path.exists():
        raise FileNotFoundError(
            f"directional-accuracy report not found at {path}. The reliability tags are read "
            "from its CURRENT STANDING section and there is no fallback copy of the numbers."
        )
    text = path.read_text(encoding="utf-8")
    index = text.find(_CURRENT_STANDING_HEADING)
    if index < 0:
        raise ValueError(
            f"{path} has no '{_CURRENT_STANDING_HEADING}' section. That section is the only "
            "place the live per-category panel counts are recorded; see its own "
            "'How to update this table' subsection."
        )
    counts = _parse_count_tables(text[index:])
    if not counts:
        raise ValueError(f"{path}: CURRENT STANDING section contains no parseable count table.")
    return counts


def verdict_for(agree: int, evaluable: int) -> str:
    """Apply the documented thresholds. Pure arithmetic, no tuning."""
    if evaluable <= 0:
        return VERDICT_DO_NOT_USE
    rate = agree / evaluable
    if rate < CAUTION_MIN_RATE:
        return VERDICT_DO_NOT_USE
    if rate >= TRUST_MIN_RATE and evaluable >= MIN_EVALUABLE_FOR_TRUST:
        return VERDICT_TRUST
    return VERDICT_CAUTION


_AXIS_NOTES = {
    "ph": (
        "pH is 6/10 on the panel and 6/13 combined with water activity -- BELOW chance when the "
        "two axes are pooled, and EXACTLY ON the caution floor on its own: one more miss makes it "
        "do-not-use (Wave X, 2026-08-28; was 6/9 and 6/12 at Wave W, and 4/7 and 4/10 before "
        "that). The Wave X row is the sharpest of them: Hofmann's hydroxyacetaldehyde + "
        "mercapto-2-propanone system, which has no amino acid and no sugar to confound it, "
        "measures MFT rising 20x from pH 3 to pH 7 and the model predicts the IDENTICAL value at "
        "pH 3, 5 and 7, because that lane carries no pH term at all. The model licenses no pH "
        "recommendation."
    ),
    "moisture_aw": (
        "Water activity is 0/3 and the miss is STRUCTURAL, not a wiring bug: HMF and furfural "
        "share one reaction family, so they always carry the same aw factor, and the "
        "projection budget has no moisture dependence at all."
    ),
    "temperature": (
        "6/8 on the panel, but the free-precursor hold-out found acrylamide ~40x "
        "under-responsive in time and the sulfur branch's temperature dependence BACKWARDS "
        "between 100 and 130 C. Treat a temperature direction as a hypothesis."
    ),
    "lipid_lane": (
        "2/2 on the panel, but on only two claims, and the matrix lane's absolute lipid "
        "aldehydes carry the campaign's largest errors (up to ~94x). Direction only."
    ),
    "matrix_identity": (
        "One evaluable claim. The pea-vs-soy question is genuinely open: binding physics and "
        "one literature dataset disagree on the SIGN of the pea-vs-soy hexanal difference."
    ),
}


def reliability_for_axis(
    axis: str, counts: Optional[Mapping[str, tuple]] = None
) -> AxisReliability:
    table = dict(counts if counts is not None else load_panel_counts())
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
        note=_AXIS_NOTES.get(axis, ""),
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
    axes = axes_exercised(spec_a, spec_b)
    reliabilities = [reliability_for_axis(axis, table) for axis in axes]
    governing = weakest(reliabilities)
    return {
        "axes": axes,
        "per_axis": reliabilities,
        "governing": governing,
        "no_axis_differs": not axes,
    }
