"""
src/kinetic_core/vessel.py -- THE PHYSICAL STATE OF THE POT (programme step R1, 2026-09-06).

WHY THIS EXISTS
===============
The 2026-09-06 review (tasks/data_restructure_plan.md, "Reaction-modelling programme") found
that the largest between-laboratory term in the sulfur lane's residuals is not noise but an
INPUT nobody had recorded: how much oxygen the closed vessel held above the sample. Yiltirak
2026 heated 3 mL under a 17 mL air headspace (0.148 mmol O2 over 0.075 mmol cysteine, ratio
2.0); Bolton 1994 heated 33 g in a 125 mL sealed vial (ratio ~2); Hofmann & Schieberle 1998
heated 100 mL in a 200 mL autoclave (0.87 mmol O2 over 3.3 mmol cysteine, ratio 0.26). The core
is near-unbiased on Hofmann and over-predicts the thiols 9-100x on the other two. Thiols are the
oxidation-labile products, and the network's oxidant pool ``OX`` is charged from cystine only.

WHAT THIS MODULE DOES, AND DOES NOT DO
======================================
It reads the ``conditions.vessel`` block a bundle carries (completed from the source paper by
``scripts/generators/complete_benchmark_vessel_fields.py``, with per-field provenance, under the
same licence as the buffer completion: FIT_HOLDOUT_DECLARATION.md Amendment 19) and turns it
into numbers the scorecard can print next to every row: headspace volume, oxygen in it at
closure, thiol charged, and their ratio. It changes NO prediction. Charging the network's
oxidant pool from this block is programme step R2(a), a pre-registered wave, not this module.

THE ARITHMETIC, STATED
======================
Headspace = vessel volume - fill volume. Air at closure is taken at 1 atm and 20 C, 20.946 %
O2 by volume: n(O2) = 0.20946 * P * V / (R * T) = 8.706e-3 mmol per mL of headspace. Dissolved
O2 in the liquid (~0.27 mM air-saturated at 20 C) is added for completeness; it is <1 % of the
headspace term in every bundle here. "Thiol charged" sums the precursors whose molecule carries
a free -SH: cysteine, hydrogen sulfide, mercapto-2-propanone. Thiamine is sulfur but not a
thiol and is not counted. A fill volume given in grams is taken as millilitres at 1 g/mL and the
block's note says so.

STATUSES
========
``computed``      -- fill and vessel volumes stated; the numbers are printed.
``ambiguous``     -- a closed vessel, but the source leaves the volumes ambiguous or unstated.
``open``          -- an open vessel (oven, dry bath): oxygen is not limited by a headspace.
``continuous``    -- an extruder or UHT line: no fixed headspace; oxygen access is process-specific.
``not_applicable``-- unheated ingredient, commercial product, or synthetic snapshot.
``missing``       -- the bundle carries no vessel block (should not happen on the panel; tested).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

#: mmol of O2 per mL of air-filled headspace at 1 atm, 20 C (0.20946 * 101325 / (8.314 * 293.15) / 1000).
O2_MMOL_PER_ML_HEADSPACE = 0.20946 * 101325.0 / (8.314 * 293.15) / 1000.0
#: mM of dissolved O2 in air-saturated water at 20 C (Henry's law; ~8.7 mg/L).
DISSOLVED_O2_MM_AT_20C = 0.27

ATMOSPHERES = (
    "air_by_default",     # closed under laboratory air; no inert-gas step in the methods
    "air_stated",         # the source says air explicitly
    "inert_gas",          # flushed with N2 / argon BEFORE or DURING heating
    "open",               # open vessel: oven, dry bath, uncovered
    "continuous_process", # extruder, UHT line, retort belt
    "not_applicable",     # unheated / commercial product / synthetic
    "unstated",           # closed vessel whose atmosphere the source does not describe
)
WATER_SOURCES = ("deionised", "distilled", "tap", "unstated", "not_applicable")
PROVENANCE_CLASSES = (
    "primary_source_pdf", "repo_verbatim_methods_quote", "unknown", "not_applicable",
)
REQUIRED_FIELDS = (
    "fill_mL", "vessel_mL", "closure", "atmosphere", "water_source", "stirred",
    "provenance_class", "provenance_note",
)

#: Precursor names (lower-case substrings) whose molecule carries a free thiol.
THIOL_PRECURSOR_MARKERS = ("cystein", "hydrogen sulfide", "mercapto")


@dataclass(frozen=True)
class VesselSpec:
    fill_mL: Optional[float]
    vessel_mL: Optional[float]
    closure: str
    atmosphere: str
    water_source: str
    stirred: Optional[bool]
    provenance_class: str
    provenance_note: str

    @property
    def headspace_mL(self) -> Optional[float]:
        if self.fill_mL is None or self.vessel_mL is None:
            return None
        return max(float(self.vessel_mL) - float(self.fill_mL), 0.0)

    def o2_mmol(self) -> Optional[float]:
        """O2 available to the cook at closure: headspace air plus dissolved O2."""
        hs = self.headspace_mL
        if hs is None or self.atmosphere not in ("air_by_default", "air_stated"):
            return None
        dissolved = DISSOLVED_O2_MM_AT_20C * float(self.fill_mL) / 1000.0
        return hs * O2_MMOL_PER_ML_HEADSPACE + dissolved


def vessel_from_bundle(bench: Mapping[str, Any]) -> Optional[VesselSpec]:
    block = (bench.get("conditions") or {}).get("vessel")
    if not isinstance(block, dict):
        return None
    return VesselSpec(
        fill_mL=None if block.get("fill_mL") is None else float(block["fill_mL"]),
        vessel_mL=None if block.get("vessel_mL") is None else float(block["vessel_mL"]),
        closure=str(block.get("closure", "")),
        atmosphere=str(block.get("atmosphere", "unstated")),
        water_source=str(block.get("water_source", "unstated")),
        stirred=block.get("stirred"),
        provenance_class=str(block.get("provenance_class", "unknown")),
        provenance_note=str(block.get("provenance_note", "")),
    )


def thiol_mmol(bench: Mapping[str, Any], fill_mL: Optional[float]) -> Optional[float]:
    """mmol of free-thiol precursor charged, from the bundle's mM and the fill volume."""
    if fill_mL is None:
        return None
    total = 0.0
    for name, spec in (bench.get("precursors") or {}).items():
        if any(marker in str(name).lower() for marker in THIOL_PRECURSOR_MARKERS):
            conc = (spec or {}).get("concentration_mM")
            if conc is not None:
                total += float(conc) * float(fill_mL) / 1000.0
    return total


def oxygen_record(bench: Mapping[str, Any]) -> Dict[str, Any]:
    """
    The scorecard's per-benchmark oxygen line. Never raises; every branch returns a
    ``status`` and a one-sentence ``basis`` a reader can act on.
    """
    spec = vessel_from_bundle(bench)
    if spec is None:
        return {"status": "missing", "basis": "the bundle carries no conditions.vessel block"}
    common = {
        "atmosphere": spec.atmosphere,
        "water_source": spec.water_source,
        "provenance_class": spec.provenance_class,
        "fill_mL": spec.fill_mL,
        "vessel_mL": spec.vessel_mL,
        "headspace_mL": spec.headspace_mL,
    }
    if spec.atmosphere == "not_applicable":
        return {**common, "status": "not_applicable", "basis": spec.closure or "no cook"}
    if spec.atmosphere == "open":
        return {**common, "status": "open",
                "basis": "open vessel: oxygen is not limited by a headspace"}
    if spec.atmosphere == "continuous_process":
        return {**common, "status": "continuous",
                "basis": "continuous process: no fixed headspace; oxygen access is process-specific"}
    if spec.atmosphere == "inert_gas":
        return {**common, "status": "computed", "o2_mmol": 0.0,
                "thiol_mmol": thiol_mmol(bench, spec.fill_mL), "o2_to_thiol": 0.0,
                "basis": "inert-gas atmosphere stated: O2 taken as zero"}
    o2 = spec.o2_mmol()
    if o2 is None:
        return {**common, "status": "ambiguous",
                "basis": "closed vessel with volumes or atmosphere the source leaves unstated"}
    thiol = thiol_mmol(bench, spec.fill_mL)
    ratio = None if not thiol else o2 / thiol
    return {
        **common,
        "status": "computed",
        "o2_mmol": o2,
        "thiol_mmol": thiol,
        "o2_to_thiol": ratio,
        "basis": (
            f"{spec.headspace_mL:.1f} mL headspace of air at closure = {o2:.3f} mmol O2"
            + (f" over {thiol:.3f} mmol thiol charged = {ratio:.2f} mol/mol" if ratio is not None
               else "; no free-thiol precursor charged")
            + f" ({spec.atmosphere.replace('_', ' ')}; water: {spec.water_source})"
        ),
    }


def format_o2_to_thiol(record: Mapping[str, Any]) -> str:
    """One cell for the scorecard table."""
    status = record.get("status")
    if status == "computed":
        ratio = record.get("o2_to_thiol")
        if ratio is None:
            return f"{record['o2_mmol']:.3f} mmol O2, no thiol"
        return f"{ratio:.2f}"
    return str(status or "-")


__all__ = [
    "ATMOSPHERES", "DISSOLVED_O2_MM_AT_20C", "O2_MMOL_PER_ML_HEADSPACE", "PROVENANCE_CLASSES",
    "REQUIRED_FIELDS", "THIOL_PRECURSOR_MARKERS", "VesselSpec", "WATER_SOURCES",
    "format_o2_to_thiol", "oxygen_record", "thiol_mmol", "vessel_from_bundle",
]
