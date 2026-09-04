"""How a benchmark's number was quantified, and which envelope bands that implies.

2026-09-03: lifted out of ``src/kinetic_core/panel.py`` so that ``scripts/ci/schema_gate.py``
(which runs in CI with only pyyaml and jsonschema installed) can classify a bundle without
importing the kinetic core and its numpy dependency. ``panel`` re-exports every name.
"""
from __future__ import annotations

from typing import Any, Mapping, Optional, Tuple

#: The bundle quantified the analyte in the HEADSPACE (HS-SPME, static
#: headspace): the core's declared K_aw band and the HS-SPME same-sample
#: dispersion are measurement facts of THIS number.
QUANTIFICATION_HEADSPACE = "headspace"
#: The bundle quantified the analyte in the LIQUID/EXTRACT (SIDA after solvent
#: extraction, HPLC-UV, LC-MS/MS, internal-standard GC/MS of an extract): no
#: air/water partition step and no SPME fibre enters the number, so neither
#: headspace band applies to it.
QUANTIFICATION_EXTRACTION = "extraction"
#: The bundle declares no ``quantification_class``. The envelope then keeps
#: the core's own convention (``matrix_oav.absolute_concentration``: every
#: absolute ppb carries the HS-SPME dispersion) and says so on the row.
QUANTIFICATION_UNDECLARED = "undeclared"

_HEADSPACE_MARKERS = ("spme", "headspace", "hs-gc", "hs_gc", "dhs", "dynamic_headspace")
_EXTRACTION_MARKERS = (
    "isotope_dilution", "sida", "hplc", "lcms", "lc-ms", "internal_standard",
    "external_standard", "response_factor", "calibration_curve",
)


def _first_value(node: Any, key: str) -> Optional[str]:
    """Depth-first first occurrence of ``key`` in a nested bundle."""
    if isinstance(node, Mapping):
        if key in node and isinstance(node[key], str):
            return node[key]
        for value in node.values():
            found = _first_value(value, key)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = _first_value(value, key)
            if found is not None:
                return found
    return None


def quantification_family(bench: Mapping[str, Any]) -> Tuple[str, str]:
    """
    ``(family, why)`` for a bundle: one of :data:`QUANTIFICATION_HEADSPACE`,
    :data:`QUANTIFICATION_EXTRACTION`, :data:`QUANTIFICATION_UNDECLARED`.

    Read from the bundle's own ``quantification_class`` (written during the
    primary-source verification waves under ``content_verification``); the
    headspace markers are checked first because a headspace class can also
    name its calibration ("... internal standard ..., SPME-GC-MS").
    """
    declared = _first_value(bench, "quantification_class")
    if not declared:
        return QUANTIFICATION_UNDECLARED, "no quantification_class in the bundle"
    lowered = declared.lower()
    if any(marker in lowered for marker in _HEADSPACE_MARKERS):
        return QUANTIFICATION_HEADSPACE, f"quantification_class={declared[:80]!r}"
    if any(marker in lowered for marker in _EXTRACTION_MARKERS):
        return QUANTIFICATION_EXTRACTION, f"quantification_class={declared[:80]!r}"
    return QUANTIFICATION_UNDECLARED, f"quantification_class={declared[:80]!r} matches no known family"


__all__ = [
    "QUANTIFICATION_EXTRACTION",
    "QUANTIFICATION_HEADSPACE",
    "QUANTIFICATION_UNDECLARED",
    "quantification_family",
]
