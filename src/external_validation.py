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
import math
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml

ROOT = Path(__file__).resolve().parents[1]
FLAVOR_REFERENCE_PATH = ROOT / "data" / "lit" / "flavor_reference_payloads.json"
BENCHMARK_INTAKE_REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"
BENCHMARK_DIR = ROOT / "data" / "benchmarks"
EXTERNAL_VALIDATION_PROTOCOL_DIR = ROOT / "data" / "protocols" / "external_validation"
EXTERNAL_VALIDATION_BENCHMARK_DIR = ROOT / "data" / "benchmarks" / "external_validation"
EXTERNAL_VALIDATION_EVIDENCE_CLASS = "external_validation_only"

# Read the live uncalibrated matrix sigma so the report's methodology disclosure
# can never drift from the tier the envelopes are actually drawn with (it did:
# the text asserted 2.0 for a day after the value was raised to 2.86).
from src.uncertainty_propagation import (  # noqa: E402
    DEFAULT_UNCALIBRATED_OBSERVABLE_PRIORS as _UNCALIBRATED_PRIORS,
)

_UNCALIBRATED_MATRIX_SIGMA = float(_UNCALIBRATED_PRIORS["matrix_headspace"])
_UNCALIBRATED_MATRIX_FOLD = math.exp(1.645 * _UNCALIBRATED_MATRIX_SIGMA)

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

# Registry-backed Lane A.2 synthesis specs. Each bundle becomes one intake YAML
# plus one isolated benchmark JSON under data/benchmarks/external_validation/.
_HOLDOUT_BUNDLE_SPECS: Tuple[Dict[str, Any], ...] = (
    {
        "bundle_id": "external_validation_bi_2020_raw_pea_hexanal",
        "anchor_ids": ("bi_2020_raw_pea_hexanal_point",),
        "matrix_context": "raw_pea_flour",
        "protein_type": "pea_iso",
        "process_state": "ambient_slurry",
        "matrix_format": "raw pea flour external hold-out, proxied onto the pea-isolate matrix-only lane",
        "conditions": {
            "temp_C": 40.0,
            "ph": 6.0,
            "water_activity": 0.95,
            "time_min": 10.0,
        },
        "precursors": {
            "Pea Protein Isolate": {"concentration_mM": 1000.0},
        },
        "benchmark_alignment": {
            "target_benchmark_id": "pea_isolate_40C_PratapSingh2021",
            "notes": "Proxy the raw-pea literature baseline onto the closest executable pea matrix-only lane. evidence_class=external_validation_only keeps it out of calibration.",
        },
        "analytical_context": {
            "headspace_method": "HS-SPME-GC-SIM",
            "quantification_mode": "internal_standard_calibrated",
            "replicates": 3,
            "notes": "Raw-pea flour baseline from Bi et al. (2020). Benchmark conditions are copied from the closest executable pea matrix-only lane because the literature point is a structural hold-out, not a closure benchmark.",
        },
        "denaturation_state": 0.0,
    },
    {
        "bundle_id": "external_validation_bi_2020_roasted_pea_hexanal",
        "anchor_ids": ("bi_2020_roasted_pea_hexanal_point",),
        "matrix_context": "roasted_pea_flour_160c_30min",
        "protein_type": "pea_iso",
        "process_state": "heated_matrix",
        "matrix_format": "roasted pea flour external hold-out, proxied onto the pea-isolate heated matrix lane",
        "conditions": {
            "temp_C": 160.0,
            "ph": 6.0,
            "water_activity": 0.35,
            "time_min": 30.0,
        },
        "precursors": {
            "Pea Protein Isolate": {"concentration_mM": 1000.0},
        },
        "benchmark_alignment": {
            "target_benchmark_id": "pea_isolate_40C_PratapSingh2021",
            "notes": "Roasted-pea carryover point from Bi et al. (2020). The thermal condition comes from the paper; the matrix proxy remains the pea-isolate matrix-only lane and stays external_validation_only.",
        },
        "analytical_context": {
            "headspace_method": "HS-SPME-GC-SIM",
            "quantification_mode": "internal_standard_calibrated",
            "replicates": 3,
            "notes": "Roasted-pea flour anchor after 160 C / 30 min roasting. Water activity is carried as an extrusion-style low-moisture proxy because the paper does not report a benchmark-ready aw.",
        },
        "denaturation_state": 0.2,
    },
    {
        "bundle_id": "external_validation_liu_2023_ppi_offnote_baseline",
        "anchor_ids": (
            "liu_2023_ppi_hexanal_band",
            "liu_2023_ppi_nonanal_band",
        ),
        "matrix_context": "native_pea_protein_isolate_aqueous_slurry",
        "protein_type": "pea_iso",
        "process_state": "ambient_slurry",
        "matrix_format": "native commercial pea protein isolate aqueous slurry external hold-out",
        "conditions": {
            "temp_C": 40.0,
            "ph": 6.0,
            "water_activity": 0.95,
            "time_min": 10.0,
        },
        "precursors": {
            "Pea Protein Isolate": {"concentration_mM": 1000.0},
        },
        "benchmark_alignment": {
            "target_benchmark_id": "pea_isolate_40C_PratapSingh2021",
            "notes": "Lot-to-lot PPI off-note baseline from Liu (2023). Conditions are proxied onto the existing ambient pea-isolate matrix-only lane so the bundle stays executable but external-only.",
        },
        "analytical_context": {
            "headspace_method": "HS-SPME-GC-MS plus GC-O",
            "quantification_mode": "internal_standard_calibrated",
            "replicates": 6,
            "notes": "Band midpoints are used as hold-out reference values; uncertainty is derived from the reported band width.",
        },
        "denaturation_state": 0.0,
    },
    {
        "bundle_id": "external_validation_li_2026_spi_wg_hme_control",
        "anchor_ids": (
            "li_2026_spi_wg_hme_hexanal_control_point",
            "li_2026_spi_wg_hme_nonanal_control_point",
            "li_2026_spi_wg_hme_1_hexanol_control_point",
            "li_2026_spi_wg_hme_2_pentylfuran_control_point",
        ),
        "matrix_context": "spi_wheat_gluten_hme_control_57pct_moisture",
        "protein_type": "soy_iso",
        "process_state": "extrusion_structured",
        "matrix_format": "SPI:wheat gluten 6:4 dry-basis HME control, proxied onto the soy-isolate extrusion lane",
        "conditions": {
            "temp_C": 160.0,
            "ph": 7.0,
            "water_activity": 0.35,
            "time_min": 0.4167,
        },
        "precursors": {
            "Soy Protein Isolate": {"concentration_mM": 1000.0},
        },
        "benchmark_alignment": {
            "notes": "Li et al. (2026) HME control anchor. The paper does not publish a final-blend pH or water-activity closure, so this bundle stays external_validation_only and never enters calibration.",
        },
        "analytical_context": {
            "headspace_method": "HS-SPME-GC-MS with isotopic internal standards",
            "quantification_mode": "isotope_dilution_calibrated",
            "replicates": 3,
            "notes": "The control matrix includes 40% wheat gluten hydrolysate on a dry basis. The executable proxy uses soy_iso only so the hold-out remains a narrow external check, not a benchmark-closure claim.",
        },
        "denaturation_state": 0.8,
    },
)
_HOLDOUT_ANCHOR_IDS: frozenset[str] = frozenset(
    anchor_id
    for spec in _HOLDOUT_BUNDLE_SPECS
    for anchor_id in spec["anchor_ids"]
)
_MAPPED_MATRIX_TAGS: frozenset[str] = frozenset(spec["matrix_context"] for spec in _HOLDOUT_BUNDLE_SPECS)
_TOP_LEVEL_BENCHMARK_EQUIVALENTS: Mapping[str, str] = {
    # The Hernandez PBMA furfural panel point is already represented by the
    # executable Resconi benchmark subset at the top-level calibration panel.
    "hernandez_2023_furfural_ratio_anchor": "resconi_2023_pbma_beef_identity_benchmark",
}


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


@dataclass(frozen=True)
class HoldoutBundle:
    bundle_id: str
    anchor_ids: Tuple[str, ...]
    intake_payload: Dict[str, Any]
    benchmark_payload: Dict[str, Any]

    def matched_compound_count(self) -> int:
        return len(self.benchmark_payload.get("measured_volatiles", {}))


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
    if kind in {"band", "oav_band", "range"}:
        lo = numeric.get("low")
        hi = numeric.get("high")
        if lo is None and hi is None:
            lo = numeric.get("min")
            hi = numeric.get("max")
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)):
            mid = _geo_mid(float(lo), float(hi))
            return mid, float(lo), float(hi)
        return None, None, None
    if kind == "band_comparison":
        band = numeric.get("pbma_band") or numeric.get("band")
        if isinstance(band, (list, tuple)) and len(band) == 2 and all(isinstance(x, (int, float)) for x in band):
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


def _load_flavor_anchor_lookup(
    flavor_reference_path: Path = FLAVOR_REFERENCE_PATH,
) -> Dict[str, Dict[str, Any]]:
    raw = json.loads(flavor_reference_path.read_text(encoding="utf-8"))
    lookup: Dict[str, Dict[str, Any]] = {}
    for section, body in raw.items():
        if not isinstance(body, list) or not section.endswith("_anchors"):
            continue
        for anchor in body:
            if not isinstance(anchor, Mapping):
                continue
            anchor_id = str(anchor.get("id", "")).strip()
            if not anchor_id:
                continue
            lookup[anchor_id] = {"section": section, **anchor}
    return lookup


def _load_registry_artifact_lookup(
    registry_path: Path = BENCHMARK_INTAKE_REGISTRY_PATH,
) -> Dict[str, Dict[str, Any]]:
    raw = json.loads(registry_path.read_text(encoding="utf-8"))
    lookup: Dict[str, Dict[str, Any]] = {}
    for row in raw.get("eligible_references", []):
        if not isinstance(row, Mapping):
            continue
        for artifact in row.get("runtime_artifacts", []) or []:
            if not isinstance(artifact, Mapping):
                continue
            artifact_id = str(artifact.get("artifact_id", "")).strip()
            if artifact_id:
                lookup[artifact_id] = dict(row)
    return lookup


def _anchor_has_top_level_benchmark_artifact(row: Mapping[str, Any]) -> bool:
    for artifact in row.get("runtime_artifacts", []) or []:
        if not isinstance(artifact, Mapping):
            continue
        if str(artifact.get("artifact_type", "")).strip() != "benchmark":
            continue
        path = str(artifact.get("path", "")).strip()
        if path.startswith("data/benchmarks/") and "/external_validation/" not in path:
            return True
    return False


def _bundle_spec_for_anchor(anchor_id: str) -> Optional[Dict[str, Any]]:
    for spec in _HOLDOUT_BUNDLE_SPECS:
        if anchor_id in spec["anchor_ids"]:
            return dict(spec)
    return None


def _publication_year_proxy(citation: str) -> str:
    digits = [part for part in str(citation).split() if part.strip("(),")[:4].isdigit()]
    for part in digits:
        cleaned = "".join(ch for ch in part if ch.isdigit())
        if len(cleaned) >= 4:
            return f"{cleaned[:4]}-01-01"
    return "unspecified"


def _measurement_row(candidate: ExternalCandidate) -> Dict[str, Any]:
    if candidate.point_value is None:
        raise ValueError(f"Candidate {candidate.anchor_id} has no numeric point value")
    row: Dict[str, Any] = {"conc_ppb": float(candidate.point_value)}
    if (
        candidate.band_low is not None
        and candidate.band_high is not None
        and candidate.point_value > 0.0
    ):
        half_span = max(
            abs(candidate.point_value - candidate.band_low),
            abs(candidate.band_high - candidate.point_value),
        )
        row["uncertainty_pct"] = round((half_span / candidate.point_value) * 100.0, 1)
    return row


def _unit_is_ppb_equivalent(units: str) -> bool:
    return str(units).strip().lower() in {"ppb", "ug_per_kg", "ng_per_g"}


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
    registry_lookup = _load_registry_artifact_lookup()
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
            registry_row = registry_lookup.get(anchor_id, {})
            has_top_level_benchmark = _anchor_has_top_level_benchmark_artifact(registry_row) or (
                anchor_id in _TOP_LEVEL_BENCHMARK_EQUIVALENTS
            )
            executable_bundle = _bundle_spec_for_anchor(anchor_id)

            # Decide eligibility. The strictest classification wins.
            if point is None and lo is None and hi is None:
                eligibility = "narrative_only"
                reason = "numeric block has no recognised point or band shape"
            elif has_top_level_benchmark:
                eligibility = "redundant_with_panel"
                reason = (
                    f"anchor '{anchor_id}' already has a top-level benchmark runtime artifact in "
                    f"data/benchmarks ({_TOP_LEVEL_BENCHMARK_EQUIVALENTS.get(anchor_id, 'registry-backed benchmark')}); "
                    "it is not an external hold-out"
                )
            elif not _unit_is_ppb_equivalent(units):
                eligibility = "narrative_only"
                reason = (
                    f"units '{units}' are not directly ppb-equivalent; an explicit conversion contract "
                    "is required before this anchor can be scored"
                )
            elif executable_bundle is not None:
                eligibility = "executable_candidate"
                reason = (
                    f"anchor '{anchor_id}' is wired to hold-out bundle '{executable_bundle['bundle_id']}' "
                    "and stays isolated from the default calibration panel"
                )
            elif not in_panel:
                eligibility = "narrative_only"
                reason = (
                    f"compound '{compound}' (canonical '{canonical}') is not in the "
                    f"current calibration panel and does not yet have a curated hold-out bundle"
                )
            else:
                eligibility = "narrative_only"
                reason = (
                    f"anchor '{anchor_id}' has no Lane A.2 synthesis bundle yet; matrix_context "
                    f"'{matrix_context}' may still need a dedicated executable mapping"
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
            # extend the Lane A.2 bundle table if appropriate.
            unmapped = sorted({c.matrix_context for c in rows if c.matrix_context not in _MAPPED_MATRIX_TAGS})
            if unmapped:
                lines.append("Unmapped matrix tags (extend the hold-out bundle specs if any of these become executable):")
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


def build_holdout_bundles(
    *,
    flavor_reference_path: Path = FLAVOR_REFERENCE_PATH,
) -> List[HoldoutBundle]:
    """Build the Lane A.2 external-validation-only intake + benchmark bundles.

    The resulting benchmark JSONs are meant to live under
    data/benchmarks/external_validation/ so they remain outside the default
    calibration panel loaded by get_benchmark_files().
    """

    from src.matrix_experiment_intake import materialize_matrix_experiment_benchmark

    inventory = build_inventory(flavor_reference_path=flavor_reference_path)
    by_anchor = {row.anchor_id: row for row in inventory.candidates}
    anchor_lookup = _load_flavor_anchor_lookup(flavor_reference_path)
    bundles: List[HoldoutBundle] = []

    for spec in _HOLDOUT_BUNDLE_SPECS:
        measured_volatiles: Dict[str, Dict[str, Any]] = {}
        citations: List[str] = []
        dois: List[str] = []
        notes: List[str] = []

        for anchor_id in spec["anchor_ids"]:
            candidate = by_anchor.get(anchor_id)
            if candidate is None:
                raise KeyError(f"Unknown external-validation anchor: {anchor_id}")
            if candidate.eligibility != "executable_candidate":
                raise ValueError(
                    f"Anchor {anchor_id} must be executable_candidate to materialize; got {candidate.eligibility}"
                )
            anchor_row = anchor_lookup.get(anchor_id, {})
            measured_volatiles[candidate.compound] = _measurement_row(candidate)
            citations.append(str(anchor_row.get("source_citation", candidate.citation)))
            doi = str(anchor_row.get("doi", candidate.doi) or "")
            if doi:
                dois.append(doi)
            target_direction = str(anchor_row.get("target_direction", "")).strip()
            if target_direction:
                notes.append(f"{anchor_id}: {target_direction}")

        primary_citation = citations[0] if citations else spec["bundle_id"]
        intake_payload: Dict[str, Any] = {
            "experiment_id": spec["bundle_id"],
            "source_kind": "external_literature",
            "evidence_class": EXTERNAL_VALIDATION_EVIDENCE_CLASS,
            "protein_type": spec["protein_type"],
            "process_state": spec["process_state"],
            "matrix_format": spec["matrix_format"],
            "conditions": dict(spec["conditions"]),
            "formulation": {"precursors": dict(spec["precursors"]),},
            "measured_volatiles": measured_volatiles,
            "provenance": {
                "origin": "external_literature",
                "source_reference": primary_citation,
                "source_doi": dois[0] if dois else "",
                "measurement_date": _publication_year_proxy(primary_citation),
                "notes": (
                    "External validation hold-out synthesized from flavor_reference_payloads.json; "
                    "never use for calibration or promotion. "
                    + " ".join(notes)
                ).strip(),
            },
            "benchmark_alignment": dict(spec.get("benchmark_alignment", {})),
            "analytical_context": {
                **dict(spec.get("analytical_context", {})),
                "citation_provenance": citations,
            },
            "denaturation_state": float(spec.get("denaturation_state", 0.5)),
        }
        benchmark_payload = materialize_matrix_experiment_benchmark(intake_payload)
        bundles.append(
            HoldoutBundle(
                bundle_id=spec["bundle_id"],
                anchor_ids=tuple(spec["anchor_ids"]),
                intake_payload=intake_payload,
                benchmark_payload=benchmark_payload,
            )
        )
    return bundles


def write_holdout_bundles(
    bundles: Sequence[HoldoutBundle],
    *,
    protocol_dir: Path = EXTERNAL_VALIDATION_PROTOCOL_DIR,
    benchmark_dir: Path = EXTERNAL_VALIDATION_BENCHMARK_DIR,
) -> Dict[str, List[Path]]:
    protocol_dir.mkdir(parents=True, exist_ok=True)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    protocol_paths: List[Path] = []
    benchmark_paths: List[Path] = []
    for bundle in bundles:
        protocol_path = protocol_dir / f"{bundle.bundle_id}.yaml"
        benchmark_path = benchmark_dir / f"{bundle.bundle_id}.json"
        protocol_path.write_text(
            yaml.safe_dump(bundle.intake_payload, sort_keys=False, allow_unicode=False),
            encoding="utf-8",
        )
        benchmark_path.write_text(
            json.dumps(bundle.benchmark_payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        protocol_paths.append(protocol_path)
        benchmark_paths.append(benchmark_path)
    return {"protocols": protocol_paths, "benchmarks": benchmark_paths}


def get_holdout_benchmark_files(
    benchmark_dir: Path = EXTERNAL_VALIDATION_BENCHMARK_DIR,
) -> List[Path]:
    """Return the isolated external-validation benchmark JSON files.

    These files live under data/benchmarks/external_validation so they stay
    outside the default top-level calibration panel.
    """

    if not benchmark_dir.exists():
        return []
    return sorted(benchmark_dir.glob("*.json"))


def build_external_validation_report(
    *,
    benchmark_files: Optional[Iterable[Path | str]] = None,
    n_samples: int = 200,
    seed: int = 0,
    target_tag: str = "meaty",
) -> Dict[str, Any]:
    """Run the uncertainty propagation over the isolated hold-out panel.

    Returns a JSON-serialisable report with the same per-compound envelope
    statistics as S20.1 plus external-hold-out accuracy summaries computed from
    the P50 prediction against the measured value.
    """

    from src.uncertainty_propagation import default_priors, propagate_benchmarks

    holdout_files = [Path(item) for item in (benchmark_files or get_holdout_benchmark_files())]
    # S27 Workstream B: the hold-out bundles are evidence_class=external_validation_only,
    # i.e. UNCALIBRATED by construction — their matrix process-states are not pinned by
    # the calibration registry. Use the wide structural-ignorance observable priors so
    # the reported CIs reflect the model's genuine (large) uncertainty on these out-of-
    # calibration predictions rather than the tight in-registry priors. The sigmas are a
    # stated physical prior, not fitted to the measured hold-out values.
    envelope_payload = propagate_benchmarks(
        benchmark_files=holdout_files,
        n_samples=n_samples,
        seed=seed,
        priors=default_priors(matrix_tier="uncalibrated"),
        target_tag=target_tag,
        execution_paths=("free_precursor", "matrix_precursor_augmented", "matrix_only"),
    )

    benchmark_rows: List[Dict[str, Any]] = []
    abs_log10_errors: List[float] = []
    fold_errors: List[float] = []

    for benchmark in envelope_payload.get("benchmarks", []):
        compound_rows: List[Dict[str, Any]] = []
        inside_ci_count = 0
        local_errors: List[float] = []
        local_folds: List[float] = []

        for compound in benchmark.get("compounds", []):
            measured = float(compound.get("measured_ppb", 0.0) or 0.0)
            predicted_p50 = float(compound.get("predicted_p50", 0.0) or 0.0)
            abs_log10_error = None
            fold_error = None
            if measured > 0.0 and predicted_p50 > 0.0 and math.isfinite(measured) and math.isfinite(predicted_p50):
                ratio = predicted_p50 / measured
                abs_log10_error = abs(math.log10(ratio))
                fold_error = max(ratio, 1.0 / ratio)
                abs_log10_errors.append(abs_log10_error)
                fold_errors.append(fold_error)
                local_errors.append(abs_log10_error)
                local_folds.append(fold_error)
            inside_ci = bool(compound.get("inside_ci", False))
            inside_ci_count += int(inside_ci)
            compound_rows.append(
                {
                    **dict(compound),
                    "abs_log10_error": abs_log10_error,
                    "fold_error": fold_error,
                }
            )

        matched_compounds = int(benchmark.get("matched_compounds", 0) or 0)
        median_abs_log10_error = statistics.median(local_errors) if local_errors else None
        median_accuracy_fold = (10 ** median_abs_log10_error) if median_abs_log10_error is not None else None
        benchmark_rows.append(
            {
                "benchmark_id": benchmark.get("benchmark_id"),
                "bench_file": benchmark.get("bench_file"),
                "execution_path": benchmark.get("execution_path"),
                "matched_compounds": matched_compounds,
                "inside_ci_count": inside_ci_count,
                "inside_ci_rate": (inside_ci_count / matched_compounds) if matched_compounds else None,
                "median_abs_log10_error": median_abs_log10_error,
                "median_accuracy_fold": median_accuracy_fold,
                "compounds": compound_rows,
            }
        )

    summary = dict(envelope_payload.get("summary", {}))
    median_abs_log10_error = statistics.median(abs_log10_errors) if abs_log10_errors else None
    median_accuracy_fold = (10 ** median_abs_log10_error) if median_abs_log10_error is not None else None
    summary.update(
        {
            "holdout_benchmark_count": len(holdout_files),
            "median_abs_log10_error": median_abs_log10_error,
            "median_accuracy_fold": median_accuracy_fold,
            "max_fold_error": max(fold_errors) if fold_errors else None,
            "evidence_class": EXTERNAL_VALIDATION_EVIDENCE_CLASS,
            "benchmark_dir": str(EXTERNAL_VALIDATION_BENCHMARK_DIR),
        }
    )

    # Lane F (sprint 2026-05-10b): per-compound failing aggregate. A compound is
    # "external_failing" when its mean |log10 error| across the hold-out exceeds
    # 1.0 dex (~10x error). The runtime report layer reads this list to colour
    # the affected compounds amber regardless of inside-CI status.
    compound_errors: Dict[str, List[float]] = {}
    for benchmark in benchmark_rows:
        for compound in benchmark.get("compounds", []) or []:
            err = compound.get("abs_log10_error")
            if err is None:
                continue
            name = str(compound.get("compound", "")).strip().lower()
            if not name:
                continue
            compound_errors.setdefault(name, []).append(float(err))

    failing_compounds: List[Dict[str, Any]] = []
    failure_threshold_dex = 1.0
    for name, errors in sorted(compound_errors.items()):
        if not errors:
            continue
        mean_err = sum(errors) / len(errors)
        if mean_err > failure_threshold_dex:
            failing_compounds.append(
                {
                    "compound": name,
                    "mean_abs_log10_error": mean_err,
                    "n_holdout_observations": len(errors),
                }
            )

    return {
        "source": {
            "benchmark_dir": str(EXTERNAL_VALIDATION_BENCHMARK_DIR),
            "evidence_class": EXTERNAL_VALIDATION_EVIDENCE_CLASS,
            "description": "Hold-out matrix benchmarks synthesized from flavor_reference_payloads.json and excluded from the default calibration panel.",
        },
        "summary": summary,
        "external_failing_compounds": failing_compounds,
        "external_failing_threshold_dex": failure_threshold_dex,
        "priors": list(envelope_payload.get("priors", [])),
        "benchmarks": benchmark_rows,
    }


def render_external_validation_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {})
    coverage_rate = summary.get("ci_coverage_rate")
    coverage_str = f"{coverage_rate * 100:.1f}%" if coverage_rate is not None else "n/a"
    median_abs_log10_error = summary.get("median_abs_log10_error")
    median_abs_str = f"{float(median_abs_log10_error):.3f}" if median_abs_log10_error is not None else "n/a"
    median_accuracy_fold = summary.get("median_accuracy_fold")
    median_accuracy_str = f"{float(median_accuracy_fold):.2f}x" if median_accuracy_fold is not None else "n/a"
    lines = [
        '<!-- Auto-regenerated by ./scripts/docker_maillard.sh run "python scripts/generators/generate_external_validation_report.py". Manual edits will be overwritten. -->',
        "",
        "# External Validation Report",
        "",
        "_Monte Carlo envelope evaluation on isolated hold-out matrix bundles that are explicitly excluded from calibration via evidence_class = external_validation_only._",
        "",
        f"**Headline external trust metric**: measured value lies inside 90% CI for **{summary.get('ci_coverage_hits', 0)} / {summary.get('matched_compound_count', 0)}** matched hold-out compounds (**{coverage_str}**).",
        "",
        f"**Median accuracy on hold-outs**: **{median_accuracy_str}** median fold error (median |log10 error| = **{median_abs_str}** dex).",
        "",
        f"Samples per hold-out bundle: {summary.get('n_samples', 0)}; seed {summary.get('seed', 0)}; bundles evaluated: {summary.get('holdout_benchmark_count', 0)}.",
        "",
        "> **Methodology disclosure (audit 2026-08-26, sigma refreshed 2026-08-27).**",
        "> Hold-out envelopes are computed with the `uncalibrated` matrix prior tier",
        f"> (matrix_headspace ln-sigma {_UNCALIBRATED_MATRIX_SIGMA:.2f}, ~±{_UNCALIBRATED_MATRIX_FOLD:.0f}x at 90% CI),",
        "> which is substantially wider than the calibrated tier used",
        "> for the in-panel headline — coverage here is therefore not comparable to the",
        "> in-panel coverage number. Because that sigma was raised from 2.0 to",
        f"> {_UNCALIBRATED_MATRIX_SIGMA:.2f} on 2026-08-26 (residual-derived, see",
        "> `results/validation/matrix_sigma_residual_derivation.md`), any coverage gain",
        "> against earlier runs of this report reflects a WIDER interval, not a more",
        "> accurate prediction — read the median fold error, which is unaffected by the",
        "> prior, alongside it. Additionally, bundles whose executable conditions",
        "> are copied from an in-panel calibration benchmark re-score that anchor's",
        "> prediction at its own conditions rather than testing extrapolation; only",
        "> bundles at genuinely new process states (e.g. HME extrusion, roasting) test",
        "> out-of-envelope transfer.",
        "",
        "## Hold-out bundles",
        "",
        "| Benchmark | Matched compounds | Inside 90% CI | Median accuracy |",
        "| --- | ---: | ---: | ---: |",
    ]
    for benchmark in payload.get("benchmarks", []):
        inside_rate = benchmark.get("inside_ci_rate")
        inside_str = (
            f"{benchmark.get('inside_ci_count', 0)}/{benchmark.get('matched_compounds', 0)} ({inside_rate * 100:.1f}%)"
            if inside_rate is not None
            else f"{benchmark.get('inside_ci_count', 0)}/{benchmark.get('matched_compounds', 0)}"
        )
        benchmark_median = benchmark.get("median_accuracy_fold")
        benchmark_median_str = f"{float(benchmark_median):.2f}x" if benchmark_median is not None else "n/a"
        lines.append(
            f"| {benchmark.get('benchmark_id')} | {benchmark.get('matched_compounds', 0)} | {inside_str} | {benchmark_median_str} |"
        )
    lines.extend([
        "",
        "## Per-compound envelopes",
        "",
    ])
    for benchmark in payload.get("benchmarks", []):
        lines.append(f"### {benchmark.get('benchmark_id')}")
        lines.append("")
        lines.append(f"- Execution path: {benchmark.get('execution_path')}")
        lines.append(f"- Matched compounds: {benchmark.get('matched_compounds', 0)}")
        lines.append("")
        lines.append("| Compound | Measured (ppb) | P5 | P50 | P95 | Fold error | Inside 90% CI |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for compound in benchmark.get("compounds", []):
            inside = "yes" if compound.get("inside_ci") else "no"
            fold_error = compound.get("fold_error")
            fold_error_str = f"{float(fold_error):.2f}x" if fold_error is not None else "n/a"
            lines.append(
                f"| {compound.get('compound')} | {float(compound.get('measured_ppb', 0.0)):.3g} | "
                f"{float(compound.get('predicted_p5', 0.0)):.3g} | {float(compound.get('predicted_p50', 0.0)):.3g} | "
                f"{float(compound.get('predicted_p95', 0.0)):.3g} | {fold_error_str} | {inside} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_external_validation_artifact(
    payload: Mapping[str, Any],
    *,
    output_dir: Path | str = ROOT / "results" / "validation",
    basename: str = "external_validation_report",
) -> Dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{basename}.json"
    md_path = output_dir / f"{basename}.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    md_path.write_text(render_external_validation_markdown(payload), encoding="utf-8")
    # Lane F sidecar: small file the runtime report layer reads to flag
    # external-failing compounds without parsing the full report.
    failing_path = output_dir / "external_failing_compounds.json"
    failing_payload = {
        "external_failing_compounds": list(payload.get("external_failing_compounds", []) or []),
        "external_failing_threshold_dex": payload.get("external_failing_threshold_dex", 1.0),
    }
    failing_path.write_text(json.dumps(failing_payload, indent=2, sort_keys=True), encoding="utf-8")
    return {"json": json_path, "md": md_path, "failing_sidecar": failing_path}


__all__ = [
    "BENCHMARK_INTAKE_REGISTRY_PATH",
    "ExternalCandidate",
    "CandidateInventory",
    "EXTERNAL_VALIDATION_BENCHMARK_DIR",
    "EXTERNAL_VALIDATION_EVIDENCE_CLASS",
    "EXTERNAL_VALIDATION_PROTOCOL_DIR",
    "build_inventory",
    "build_external_validation_report",
    "build_holdout_bundles",
    "get_holdout_benchmark_files",
    "render_markdown",
    "render_external_validation_markdown",
    "write_artifact",
    "write_external_validation_artifact",
    "write_holdout_bundles",
    "FLAVOR_REFERENCE_PATH",
    "BENCHMARK_DIR",
]
