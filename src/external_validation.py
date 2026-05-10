"""S20.5 — External literature validation candidate inventory.

Lane A.1 of the trust-loop-actionability sprint.

Purpose
-------
The 19 benchmarks under :file:`data/benchmarks/` are the *calibration* panel:
the framework is fitted (or at least diagnosed) against them via
:mod:`src.uncertainty_propagation` and :mod:`src.cross_validation`. The 39
flavor-reference anchors stored in
:file:`data/lit/flavor_reference_payloads.json` are quantitative ppb / OAV
measurements in plant-protein matrices that are **never** used in calibration
— they enter the runtime only as `optimization_constraint` /
`reference_anchor` priors. They are therefore the natural pool of *external*
validation candidates: scoring the panel against them gives an honest
out-of-fit accuracy signal.

This module produces an eligibility classification artifact, not a scoring
artifact. Lane A.2 turns the eligible rows into synthetic benchmark JSON
payloads (kept under :file:`data/benchmarks/external_validation/`) and Lane
A.3 runs :func:`src.uncertainty_propagation.propagate_benchmarks` against
that hold-out set with the *frozen* current calibration.

The module is read-only relative to the calibration surface — it only
inspects the existing payload + benchmark JSON files. No prior is upgraded,
no offset is fitted, no benchmark file is mutated.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
FLAVOR_REFERENCE_PATH = ROOT / "data" / "lit" / "flavor_reference_payloads.json"
BENCHMARK_DIR = ROOT / "data" / "benchmarks"

# Compound aliases used to detect overlap between flavor-anchor compound names
# and the calibration panel compound names. Keep this list short and only
# add entries that have been manually verified.
_COMPOUND_ALIASES: Mapping[str, Tuple[str, ...]] = {
    "2-methyl-3-furanthiol": ("mft", "2-methyl-3-furanthiol (mft)"),
    "2-furfurylthiol": ("fft", "2-furfurylthiol (fft)", "furfuryl mercaptan"),
    "bis(2-methyl-3-furyl) disulfide": ("mft disulfide", "bis(2-methyl-3-furyl) disulphide"),
    "nε-(carboxymethyl)lysine (cml)": ("cml", "carboxymethyl lysine"),
    "nε-(carboxyethyl)lysine (cel)": ("cel", "carboxyethyl lysine"),
}

# Matrix tags that we currently know how to express as an executable
# benchmark precursor / process-state pair. Anchors whose matrix_context is
# in this set are flagged ``executable_candidate``; the rest are
# ``narrative_only`` for now and become eligible once matrix mapping is
# extended in Lane A.2.
_EXECUTABLE_MATRIX_TAGS: frozenset[str] = frozenset(
    {
        "pea_protein_isolate",
        "soy_protein_isolate",
        "pbma_vs_meat_panel",
        "pbma_vs_beef_comparator",
        "spi_wheat_gluten_hme",
        "raw_pea_flour",
        "boiled_beef_reference",
        "cooked_pbma_oil_comparison",
        "high_oleic_oil_model_system",
        "yeast_extract_reaction_flavor",
    }
)


@dataclass(frozen=True)
class ExternalCandidate:
    """One quantitative external-validation candidate point.

    Attributes
    ----------
    anchor_id, citation, doi:
        Provenance of the literature anchor.
    section:
        ``flavor_reference_payloads.json`` top-level section the anchor
        was loaded from (e.g. ``sulfur_reference_anchors``).
    compound:
        Compound string as it appears in the anchor.
    canonical_compound:
        Lowercased / aliased compound string used for panel-overlap checks.
    matrix_context:
        Matrix tag string from the anchor.
    units:
        Unit string declared on the anchor (e.g. ``ng_per_g``,
        ``ppb``, ``oav``).
    point_value, band_low, band_high:
        Numeric extraction. For point measurements ``point_value`` is set;
        for bands ``band_low`` / ``band_high`` are set and
        ``point_value`` carries the geometric midpoint.
    eligibility:
        One of ``executable_candidate`` (synthesizable as benchmark JSON
        in Lane A.2), ``narrative_only`` (matrix not yet mapped), or
        ``redundant_with_panel`` (compound + matrix already in calibration
        panel — would not be an external check).
    eligibility_reason:
        Human-readable note explaining the classification.
    """

    anchor_id: str
    citation: str
    doi: str
    section: str
    compound: str
    canonical_compound: str
    matrix_context: str
    units: str
    point_value: Optional[float]
    band_low: Optional[float]
    band_high: Optional[float]
    eligibility: str
    eligibility_reason: str

    def to_jsonable(self) -> Dict[str, Any]:
        return {
            "anchor_id": self.anchor_id,
            "citation": self.citation,
            "doi": self.doi,
            "section": self.section,
            "compound": self.compound,
            "canonical_compound": self.canonical_compound,
            "matrix_context": self.matrix_context,
            "units": self.units,
            "point_value": self.point_value,
            "band_low": self.band_low,
            "band_high": self.band_high,
            "eligibility": self.eligibility,
            "eligibility_reason": self.eligibility_reason,
        }


@dataclass(frozen=True)
class CandidateInventory:
    panel_compounds: Tuple[str, ...]
    candidates: Tuple[ExternalCandidate, ...]

    def by_eligibility(self) -> Dict[str, List[ExternalCandidate]]:
        groups: Dict[str, List[ExternalCandidate]] = {}
        for c in self.candidates:
            groups.setdefault(c.eligibility, []).append(c)
        return groups

    def to_jsonable(self) -> Dict[str, Any]:
        groups = self.by_eligibility()
        return {
            "schema_version": "1",
            "panel_compounds": list(self.panel_compounds),
            "summary": {
                "total_candidates": len(self.candidates),
                "executable_candidate": len(groups.get("executable_candidate", [])),
                "narrative_only": len(groups.get("narrative_only", [])),
                "redundant_with_panel": len(groups.get("redundant_with_panel", [])),
            },
            "candidates": [c.to_jsonable() for c in self.candidates],
        }


# --- numeric extraction --------------------------------------------------

def _geo_mid(low: float, high: float) -> Optional[float]:
    if low <= 0 or high <= 0:
        return None
    import math
    return float(math.sqrt(low * high))


def _extract_numeric(numeric: Mapping[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Return (point_value, band_low, band_high) from a numeric_band_or_point block.

    Supported numeric block shapes mirror what
    :file:`data/lit/flavor_reference_payloads.json` actually uses today:

    * ``{"type": "point", "value": <float>}``
    * ``{"type": "band", "low": <float>, "high": <float>}``
    * ``{"type": "band_comparison", "meat_reference": <float>,
        "pbma_band": [<low>, <high>]}`` — band is the pbma band; the
        meat_reference is recorded as ``point_value`` only when no band
        is available.
    * ``{"type": "oav_band", "low": <float>, "high": <float>}``

    Anything we don't recognise yields ``(None, None, None)``.
    """

    if not isinstance(numeric, Mapping):
        return None, None, None
    kind = str(numeric.get("type", "")).lower()
    if kind == "point":
        v = numeric.get("value")
        if isinstance(v, (int, float)):
            return float(v), None, None
        return None, None, None
    if kind in {"band", "oav_band"}:
        lo = numeric.get("low")
        hi = numeric.get("high")
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)):
            mid = _geo_mid(float(lo), float(hi))
            return mid, float(lo), float(hi)
        return None, None, None
    if kind == "band_comparison":
        band = numeric.get("pbma_band") or numeric.get("band")
        if isinstance(band, Sequence) and len(band) == 2 and all(isinstance(x, (int, float)) for x in band):
            lo = float(band[0])
            hi = float(band[1])
            return _geo_mid(lo, hi), lo, hi
        ref = numeric.get("meat_reference") or numeric.get("reference")
        if isinstance(ref, (int, float)):
            return float(ref), None, None
        return None, None, None
    # Unknown numeric shape — return raw mid if a `value` happens to exist.
    v = numeric.get("value")
    if isinstance(v, (int, float)):
        return float(v), None, None
    return None, None, None


# --- panel overlap -------------------------------------------------------

def _load_panel_compounds() -> Tuple[str, ...]:
    """Return the lowercased set of measured volatile names in the
    calibration panel."""

    out: set[str] = set()
    if not BENCHMARK_DIR.exists():
        return tuple()
    for path in sorted(BENCHMARK_DIR.glob("*.json")):
        # Skip aggregate / non-benchmark JSON files; the per-benchmark
        # files we care about expose ``measured_volatiles``.
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        measured = data.get("measured_volatiles") if isinstance(data, Mapping) else None
        if not isinstance(measured, Mapping):
            continue
        for compound in measured.keys():
            out.add(str(compound).strip().lower())
    return tuple(sorted(out))


def _canonicalise_compound(name: str) -> str:
    """Lowercase and apply known alias mapping so anchor compound
    strings line up with panel compound strings."""

    needle = name.strip().lower()
    for canonical, aliases in _COMPOUND_ALIASES.items():
        if needle == canonical:
            return canonical
        for alias in aliases:
            if needle == alias.lower():
                return canonical
    return needle


def _compound_in_panel(canonical: str, panel: Iterable[str]) -> bool:
    panel_norm = {p.strip().lower() for p in panel}
    if canonical in panel_norm:
        return True
    # Allow contains-match for paren-suffixed panel names like
    # ``2-methyl-3-furanthiol (mft)``.
    return any(canonical in p for p in panel_norm)


# --- inventory builder ---------------------------------------------------

def build_inventory(
    *,
    flavor_reference_path: Path = FLAVOR_REFERENCE_PATH,
    benchmark_dir: Path = BENCHMARK_DIR,
) -> CandidateInventory:
    """Walk the flavor-reference payload and emit a per-anchor
    eligibility classification."""

    panel = _load_panel_compounds() if benchmark_dir == BENCHMARK_DIR else _load_panel_compounds()
    raw = json.loads(flavor_reference_path.read_text(encoding="utf-8"))
    candidates: List[ExternalCandidate] = []

    for section, body in raw.items():
        if not isinstance(body, list):
            continue
        if not section.endswith("_anchors"):
            continue
        for anchor in body:
            if not isinstance(anchor, Mapping):
                continue
            anchor_id = str(anchor.get("id", ""))
            citation = str(anchor.get("source_citation", ""))
            doi = str(anchor.get("doi", ""))
            compound = str(anchor.get("compound", ""))
            matrix_context = str(anchor.get("matrix_context", ""))
            units = str(anchor.get("units", ""))
            numeric = anchor.get("numeric_band_or_point", {}) or {}
            point, lo, hi = _extract_numeric(numeric)
            canonical = _canonicalise_compound(compound)

            in_panel = _compound_in_panel(canonical, panel)
            executable_matrix = matrix_context in _EXECUTABLE_MATRIX_TAGS

            # Decide eligibility. The strictest classification wins.
            if point is None and lo is None and hi is None:
                eligibility = "narrative_only"
                reason = "numeric block has no recognised point or band shape"
            elif not in_panel:
                eligibility = "narrative_only"
                reason = (
                    f"compound '{compound}' (canonical '{canonical}') is not in the "
                    f"current calibration panel; cannot be scored against it"
                )
            elif not executable_matrix:
                eligibility = "narrative_only"
                reason = (
                    f"matrix_context '{matrix_context}' is not yet mapped to an "
                    f"executable benchmark precursor / process-state set"
                )
            else:
                # Compound is in panel AND matrix is executable: it is an
                # external validation candidate. We do not flag it as
                # ``redundant_with_panel`` purely on compound match because
                # the existing panel benchmarks at this matrix are tied to
                # specific (T, time, pH) conditions that the anchor may
                # not duplicate; full redundancy detection is Lane A.2's
                # responsibility once a synthesized payload exists.
                eligibility = "executable_candidate"
                reason = (
                    f"compound '{canonical}' is in panel and matrix "
                    f"'{matrix_context}' is mappable to an executable "
                    f"precursor / process-state set"
                )

            candidates.append(
                ExternalCandidate(
                    anchor_id=anchor_id,
                    citation=citation,
                    doi=doi,
                    section=section,
                    compound=compound,
                    canonical_compound=canonical,
                    matrix_context=matrix_context,
                    units=units,
                    point_value=point,
                    band_low=lo,
                    band_high=hi,
                    eligibility=eligibility,
                    eligibility_reason=reason,
                )
            )

    return CandidateInventory(
        panel_compounds=panel,
        candidates=tuple(candidates),
    )


# --- markdown rendering --------------------------------------------------

_BANNER = (
    '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python '
    'scripts/generators/generate_external_validation_inventory.py". '
    "Manual edits will be overwritten. -->"
)


def render_markdown(inventory: CandidateInventory) -> str:
    """Render the inventory as a scientist-readable markdown report.

    The first line is always the auto-regeneration banner so the
    artifact's provenance survives drift.
    """

    groups = inventory.by_eligibility()
    lines: List[str] = [
        _BANNER,
        "",
        "# External Validation Candidate Inventory (Lane A.1)",
        "",
        f"Total anchors scanned: **{len(inventory.candidates)}**.",
        "",
        "| Eligibility | Count |",
        "| --- | ---: |",
        f"| executable_candidate | {len(groups.get('executable_candidate', []))} |",
        f"| narrative_only | {len(groups.get('narrative_only', []))} |",
        f"| redundant_with_panel | {len(groups.get('redundant_with_panel', []))} |",
        "",
        "**executable_candidate** rows are eligible for Lane A.2 synthesis "
        "into hold-out benchmark JSONs. **narrative_only** rows need either "
        "a matrix mapping or a panel compound expansion before they can be "
        "scored. **redundant_with_panel** rows would not be an external check.",
        "",
        f"Calibration panel compounds ({len(inventory.panel_compounds)}): "
        + ", ".join(f"`{c}`" for c in inventory.panel_compounds),
        "",
    ]

    def _row(c: ExternalCandidate) -> str:
        if c.point_value is not None and c.band_low is not None:
            num = f"{c.band_low:g}–{c.band_high:g} ({c.units})"
        elif c.point_value is not None:
            num = f"{c.point_value:g} ({c.units})"
        else:
            num = "—"
        return (
            f"| `{c.anchor_id}` | {c.compound} | `{c.matrix_context}` | "
            f"{num} | {c.citation} |"
        )

    for label, header in (
        ("executable_candidate", "## Executable candidates"),
        ("narrative_only", "## Narrative-only (need mapping)"),
        ("redundant_with_panel", "## Redundant with panel"),
    ):
        rows = groups.get(label, [])
        lines.append(header)
        lines.append("")
        if not rows:
            lines.append("_None._")
            lines.append("")
            continue
        lines.append("| Anchor | Compound | Matrix | Numeric | Citation |")
        lines.append("| --- | --- | --- | --- | --- |")
        for c in rows:
            lines.append(_row(c))
        lines.append("")
        if label == "narrative_only":
            # Show distinct unmapped matrix tags so a follow-up can
            # extend ``_EXECUTABLE_MATRIX_TAGS`` if appropriate.
            unmapped = sorted({c.matrix_context for c in rows if c.matrix_context not in _EXECUTABLE_MATRIX_TAGS})
            if unmapped:
                lines.append("Unmapped matrix tags (extend `_EXECUTABLE_MATRIX_TAGS` if any of these become executable):")
                lines.append("")
                for tag in unmapped:
                    lines.append(f"- `{tag}`")
                lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_artifact(
    inventory: CandidateInventory,
    *,
    output_dir: Path,
    basename: str = "external_validation_candidates",
) -> Dict[str, Path]:
    """Write the inventory as paired ``.md`` and ``.json`` artifacts."""

    output_dir.mkdir(parents=True, exist_ok=True)
    md_path = output_dir / f"{basename}.md"
    json_path = output_dir / f"{basename}.json"
    md_path.write_text(render_markdown(inventory), encoding="utf-8")
    json_path.write_text(
        json.dumps(inventory.to_jsonable(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {"markdown": md_path, "json": json_path}


__all__ = [
    "ExternalCandidate",
    "CandidateInventory",
    "build_inventory",
    "render_markdown",
    "write_artifact",
    "FLAVOR_REFERENCE_PATH",
    "BENCHMARK_DIR",
]
