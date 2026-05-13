"""S24 — Sensory readout panel.

Translate the kernel's per-compound predicted ppb (with 90 % CI) into Odour
Activity Values (OAV = predicted / odour_threshold) and roll those into three
axes (`meaty`, `off_note`, `safety`) using the same keyword vocabulary the VoI
ranker already uses. Pure presentation layer: no kernel, projection, or
prior changes.

Inputs
------
- ``FormulationResult`` with ``predicted_ppb`` and ``uncertainty_envelopes``.
- Curated ``CompoundSpec.odour_threshold_ug_per_kg`` from
  ``data/species/{desirable_targets,off_flavour_targets}.yml`` (loaded via
  ``src.experiment_value.load_compound_specs``).

Conventions
-----------
- ``ppb`` here is μg/kg (mass), which matches the curated ODT units, so OAV is
  unitless: ``OAV = predicted_ppb / odour_threshold_ug_per_kg``.
- A compound is "above threshold" when ``OAV >= 1``.
- ``None`` OAV means the compound has no curated ODT; we surface it but do not
  let it inflate or deflate any axis.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional

from src.experiment_value import (
    _MEATY_KEYWORDS,
    _OFFNOTE_KEYWORDS,
    _SAFETY_KEYWORDS,
    CompoundSpec,
    load_compound_specs,
    lookup_spec,
)

AXIS_MEATY = "meaty"
AXIS_OFF_NOTE = "off_note"
AXIS_SAFETY = "safety"
AXIS_UNCLASSIFIED = "unclassified"

_AXIS_ORDER = (AXIS_MEATY, AXIS_OFF_NOTE, AXIS_SAFETY)


@dataclass(frozen=True)
class CompoundOAV:
    """Per-compound OAV with the same 90 % CI envelope as the kernel."""

    compound: str
    axis: str
    odour_threshold_ug_per_kg: Optional[float]
    predicted_ppb: float
    predicted_p5: Optional[float]
    predicted_p50: Optional[float]
    predicted_p95: Optional[float]
    oav: Optional[float]
    oav_p5: Optional[float]
    oav_p95: Optional[float]
    above_threshold: bool


@dataclass(frozen=True)
class AxisRollup:
    """Axis-level summary built from per-compound OAVs."""

    axis: str
    compounds_with_odt: int
    compounds_without_odt: int
    max_oav: Optional[float]
    sum_oav: Optional[float]
    above_threshold_count: int
    top_contributor: Optional[str]


@dataclass(frozen=True)
class SensoryReadout:
    """Top-level readout consumed by the report renderer."""

    per_compound: List[CompoundOAV] = field(default_factory=list)
    axes: Dict[str, AxisRollup] = field(default_factory=dict)
    headline: str = ""
    notes: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# OAV math
# ---------------------------------------------------------------------------

def compute_oav(
    predicted_ppb: Optional[float],
    odour_threshold_ug_per_kg: Optional[float],
) -> Optional[float]:
    """Return ``predicted / odour_threshold`` or ``None`` when undefined.

    Returns ``None`` when ODT is missing/non-positive or when ``predicted_ppb``
    is ``None``/non-finite. Negative predictions clamp to ``0.0`` rather than
    propagating sign noise.
    """
    if odour_threshold_ug_per_kg is None:
        return None
    try:
        odt = float(odour_threshold_ug_per_kg)
    except (TypeError, ValueError):
        return None
    if odt <= 0.0 or odt != odt:  # non-positive or NaN
        return None
    if predicted_ppb is None:
        return None
    try:
        ppb = float(predicted_ppb)
    except (TypeError, ValueError):
        return None
    if ppb != ppb:  # NaN
        return None
    if ppb < 0.0:
        ppb = 0.0
    return ppb / odt


def classify_axis(compound_name: str) -> str:
    """Return one of ``meaty | off_note | safety | unclassified``.

    Mirrors ``_suggest_template`` in ``src/experiment_value.py`` — safety wins
    over meaty wins over off-note when keywords overlap.
    """
    name = (compound_name or "").lower()
    if any(k in name for k in _SAFETY_KEYWORDS):
        return AXIS_SAFETY
    if any(k in name for k in _MEATY_KEYWORDS):
        return AXIS_MEATY
    if any(k in name for k in _OFFNOTE_KEYWORDS):
        return AXIS_OFF_NOTE
    return AXIS_UNCLASSIFIED


# ---------------------------------------------------------------------------
# Roll-ups
# ---------------------------------------------------------------------------

def _roll_axis(axis: str, rows: List[CompoundOAV]) -> AxisRollup:
    with_odt = [r for r in rows if r.oav is not None]
    without_odt = [r for r in rows if r.oav is None]
    if with_odt:
        max_row = max(with_odt, key=lambda r: (r.oav or 0.0))
        max_oav = max_row.oav
        top_contributor = max_row.compound
        sum_oav = sum(r.oav or 0.0 for r in with_odt)
        above_count = sum(1 for r in with_odt if r.above_threshold)
    else:
        max_oav = None
        top_contributor = None
        sum_oav = None
        above_count = 0
    return AxisRollup(
        axis=axis,
        compounds_with_odt=len(with_odt),
        compounds_without_odt=len(without_odt),
        max_oav=max_oav,
        sum_oav=sum_oav,
        above_threshold_count=above_count,
        top_contributor=top_contributor,
    )


def roll_up_axes(rows: List[CompoundOAV]) -> Dict[str, AxisRollup]:
    """Group per-compound OAVs by axis and return a stable, ordered dict."""
    grouped: Dict[str, List[CompoundOAV]] = {axis: [] for axis in _AXIS_ORDER}
    for row in rows:
        if row.axis in grouped:
            grouped[row.axis].append(row)
    return {axis: _roll_axis(axis, rs) for axis, rs in grouped.items()}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def _format_headline(axes: Mapping[str, AxisRollup]) -> str:
    parts: List[str] = []
    for axis_key, label_above, label_below, label_unknown in (
        (AXIS_MEATY, "meaty above threshold", "meaty below threshold", "meaty no ODT"),
        (AXIS_OFF_NOTE, "off-notes above threshold", "off-notes below threshold", "off-notes no ODT"),
        (AXIS_SAFETY, "safety above threshold", "safety clear", "safety no ODT"),
    ):
        rollup = axes.get(axis_key)
        if rollup is None or rollup.compounds_with_odt == 0:
            parts.append(label_unknown)
            continue
        parts.append(label_above if rollup.above_threshold_count > 0 else label_below)
    return "; ".join(parts)


def build_sensory_readout(
    result,
    specs: Optional[Mapping[str, CompoundSpec]] = None,
) -> SensoryReadout:
    """Build a ``SensoryReadout`` from a ``FormulationResult``.

    ``specs`` defaults to ``load_compound_specs()``. Pass an explicit mapping
    in tests to avoid disk I/O.
    """
    if specs is None:
        specs = load_compound_specs()

    predicted: Mapping[str, float] = getattr(result, "predicted_ppb", {}) or {}
    envelopes: Mapping[str, object] = getattr(result, "uncertainty_envelopes", {}) or {}

    rows: List[CompoundOAV] = []
    for compound, ppb in predicted.items():
        spec = lookup_spec(compound, specs)
        odt = spec.odour_threshold_ug_per_kg if spec else None
        env = envelopes.get(compound)
        p5 = getattr(env, "predicted_p5", None) if env is not None else None
        p50 = getattr(env, "predicted_p50", None) if env is not None else None
        p95 = getattr(env, "predicted_p95", None) if env is not None else None
        oav = compute_oav(p50 if p50 is not None else ppb, odt)
        oav_p5 = compute_oav(p5, odt)
        oav_p95 = compute_oav(p95, odt)
        rows.append(
            CompoundOAV(
                compound=compound,
                axis=classify_axis(compound),
                odour_threshold_ug_per_kg=float(odt) if odt is not None else None,
                predicted_ppb=float(ppb) if ppb is not None else 0.0,
                predicted_p5=float(p5) if p5 is not None else None,
                predicted_p50=float(p50) if p50 is not None else None,
                predicted_p95=float(p95) if p95 is not None else None,
                oav=oav,
                oav_p5=oav_p5,
                oav_p95=oav_p95,
                above_threshold=bool(oav is not None and oav >= 1.0),
            )
        )

    axes = roll_up_axes(rows)

    notes: List[str] = []
    no_odt_count = sum(1 for r in rows if r.oav is None)
    if no_odt_count:
        notes.append(
            f"{no_odt_count}/{len(rows)} predicted compounds have no curated odour threshold; "
            "they appear in the per-compound table but do not contribute to axis roll-ups."
        )
    if not rows:
        notes.append(
            "FormulationResult.predicted_ppb is empty — no compounds to score."
        )

    headline = _format_headline(axes) if rows else "no predicted compounds"
    return SensoryReadout(per_compound=rows, axes=axes, headline=headline, notes=notes)
