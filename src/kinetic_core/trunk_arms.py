"""
src/kinetic_core/trunk_arms.py

THE TRUNK'S OPTIONAL ARMS, AS ONE TABLE (2026-09-10, on review).
================================================================

Waves B13, B18, B20, B22 and B24 each added a set of trunk-only target species and each wrote its
own checks by hand into ``engine.declare_envelope``: a not-shipped refusal, a missing-precursor
refusal, a lane-conflict clause, a caveat on every answer and a note about the amine being charged
as glycine. By the fifth wave that function had grown from 370 lines and 40 branches to 461 and 59,
and a sixth wave would have added three more of each. This table is what those blocks had in
common; ``declare_envelope`` iterates it.

WHAT CHANGED AND WHAT DID NOT. Nothing a user sees. Every reason and warning string is the same
byte for byte, and they are emitted in the same order as before, which is why each arm carries
three explicit order keys: the refusal blocks, the lane-conflict clause and the warning blocks were
each written in a different order over five waves, and preserving that order is cheaper and
safer than re-pinning every scorecard string that quotes it. A later wave may canonicalise the
order and re-pin; this one is a refactor and refuses to change output.

WHAT STAYS OUTSIDE THE TABLE, and why. Two checks are genuinely irregular: the glycation arm
refuses on the MATRIX LAYER's charged amine sites rather than on a precursor key, and the
dimethyl-trisulfide refusal matches the RAW target string because the compound has no species key.
Both are still in ``declare_envelope`` as named special cases, next to the loop, rather than forced
into fields they do not fit.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, FrozenSet, Optional, Sequence, Tuple

DICARBONYL_TARGET_KEYS: FrozenSet[str] = frozenset({"G", "GO", "DA"})
PYRAZINE_TARGET_KEYS: FrozenSet[str] = frozenset({"PZ", "DMP", "MPZ"})
GLYCATION_TARGET_KEYS: FrozenSet[str] = frozenset({"CML", "CEL", "FLP", "LYSP"})
METHIONINE_TARGET_KEYS: FrozenSet[str] = frozenset({"MTAL", "MSH", "DMDS"})
PROLINE_TARGET_KEYS: FrozenSet[str] = frozenset({"AP", "PYRL"})


@dataclass(frozen=True)
class TrunkArm:
    """One optional arm of the trunk and everything the envelope says about it."""

    #: The capitalised label every message starts with, e.g. "METHIONINE CHAIN TARGETS".
    label: str
    #: The species keys the arm answers for.
    target_keys: FrozenSet[str]
    #: The wave tag the lane-conflict clause quotes, e.g. "B22".
    wave: str
    #: Returns (shipped, not_shipped_reason). ``None`` for an arm that always ships.
    shipped: Optional[Callable[[], Tuple[bool, str]]]
    #: Precursor keys of which at least one must be charged, and the message appended after the
    #: label and target list when none is. ``None`` for an arm with no precursor requirement.
    required_precursors: Optional[Tuple[str, ...]]
    missing_precursor_message: Optional[str]
    #: Returns the caveat strings appended to every answer that names one of the targets.
    target_caveats: Optional[Callable[[], Sequence[str]]]
    #: (precursor key, message on the trunk, message on another lane or None). The messages are
    #: format strings taking ``amount`` and ``lane``.
    charged_as_glycine: Optional[Tuple[str, str, Optional[str]]]
    #: Emission orders, kept explicit so the output is byte-identical to the hand-written blocks.
    refusal_order: int
    conflict_order: int
    warning_order: int


def _proline_shipped() -> Tuple[bool, str]:
    from .parameters_proline import PROLINE_NOT_SHIPPED_REASON, PROLINE_SHIPPED

    return bool(PROLINE_SHIPPED), PROLINE_NOT_SHIPPED_REASON


def _methionine_shipped() -> Tuple[bool, str]:
    from .parameters_methionine import METHIONINE_NOT_SHIPPED_REASON, METHIONINE_SHIPPED

    return bool(METHIONINE_SHIPPED), METHIONINE_NOT_SHIPPED_REASON


def _proline_caveats() -> Sequence[str]:
    from .parameters_proline import PROLINE_CAVEAT

    return (PROLINE_CAVEAT,)


def _methionine_caveats() -> Sequence[str]:
    from .parameters_methionine import METHIONINE_CAVEAT

    return (METHIONINE_CAVEAT,)


def _glycation_caveats() -> Sequence[str]:
    from .parameters_glycation import GLYCATION_AVAILABILITY_CAVEAT

    return (GLYCATION_AVAILABILITY_CAVEAT,)


def _pyrazine_caveats() -> Sequence[str]:
    from .parameters_pyrazine import PYRAZINE_SINK_CAVEAT, PYRAZINE_SUPPLY_CAVEAT

    return (PYRAZINE_SUPPLY_CAVEAT, PYRAZINE_SINK_CAVEAT)


TRUNK_ARMS: Tuple[TrunkArm, ...] = (
    TrunkArm(
        label="DICARBONYL TARGETS", target_keys=DICARBONYL_TARGET_KEYS, wave="B13",
        shipped=None, required_precursors=None, missing_precursor_message=None,
        target_caveats=None, charged_as_glycine=None,
        refusal_order=99, conflict_order=0, warning_order=99,
    ),
    TrunkArm(
        label="PYRAZINE TARGETS", target_keys=PYRAZINE_TARGET_KEYS, wave="B18",
        shipped=None, required_precursors=None, missing_precursor_message=None,
        target_caveats=_pyrazine_caveats, charged_as_glycine=None,
        refusal_order=99, conflict_order=1, warning_order=0,
    ),
    TrunkArm(
        label="GLYCATION TARGETS", target_keys=GLYCATION_TARGET_KEYS, wave="B20",
        shipped=None, required_precursors=None, missing_precursor_message=None,
        target_caveats=_glycation_caveats, charged_as_glycine=None,
        refusal_order=99, conflict_order=2, warning_order=3,
    ),
    TrunkArm(
        label="METHIONINE CHAIN TARGETS", target_keys=METHIONINE_TARGET_KEYS, wave="B22",
        shipped=_methionine_shipped,
        required_precursors=("MET",),
        missing_precursor_message=(
            " need methionine in the charge (wave B22): methional is methionine's Strecker aldehyde and "
            "methanethiol and the disulfide are made from it. Refused rather than answered with a structural zero."),
        target_caveats=_methionine_caveats,
        charged_as_glycine=(
            "MET",
            "methionine ({amount:g} mmol/L) is charged as GLYCINE at the same molarity for the sugar "
            "path's Amadori chemistry (declared, wave B22); its own Strecker chain did not ship and carries no flux.",
            "methionine ({amount:g} mmol/L) is carried by the sugar path only (wave B22); "
            "the {lane} lane's network has no methionine step, so it is recorded and not charged.",
        ),
        refusal_order=1, conflict_order=3, warning_order=2,
    ),
    TrunkArm(
        label="2-ACETYL-1-PYRROLINE TARGETS", target_keys=PROLINE_TARGET_KEYS, wave="B24",
        shipped=_proline_shipped,
        required_precursors=("PRO", "PYRL"),
        missing_precursor_message=(
            " need proline (or fed 1-pyrroline) in the charge (wave B24). Refused rather than answered with a structural zero."),
        target_caveats=_proline_caveats,
        charged_as_glycine=(
            "PRO",
            "proline ({amount:g} mmol/L) is charged as GLYCINE at the same molarity for the sugar "
            "path's Amadori chemistry (declared, wave B24; a secondary amine, so an upper bound on the supply it makes).",
            None,
        ),
        refusal_order=0, conflict_order=4, warning_order=1,
    ),
)
TRUNK_ONLY_TARGET_KEYS: FrozenSet[str] = frozenset().union(*(arm.target_keys for arm in TRUNK_ARMS))


def named_targets(arm: TrunkArm, mapped_targets) -> list:
    """The user's compound names that map onto this arm, sorted as the hand-written blocks sorted them."""
    return sorted(c for c, key in mapped_targets.items() if key in arm.target_keys)
