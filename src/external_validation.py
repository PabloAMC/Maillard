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

import copy
import json
import math
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml

from src import data_paths

FLAVOR_REFERENCE_PATH = data_paths.FLAVOR_REFERENCE_PAYLOADS
BENCHMARK_INTAKE_REGISTRY_PATH = data_paths.BENCHMARK_INTAKE_REGISTRY
BENCHMARK_DIR = data_paths.BENCHMARKS_DIR
# The intake YAML that each frozen hold-out benchmark was materialized from lives NEXT TO
# the benchmark JSON (moved from data/protocols/external_validation/ on 2026-09-01). Both
# are frozen evidence; see write_holdout_bundles.
EXTERNAL_VALIDATION_PROTOCOL_DIR = data_paths.EXTERNAL_VALIDATION_INTAKE_DIR
EXTERNAL_VALIDATION_BENCHMARK_DIR = data_paths.EXTERNAL_VALIDATION_DIR
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
        "provenance_overrides": {
            "shared_anchor": True,
            "citation_audit_note": {
                "date": "2026-08-26",
                "status": "verified_correct_shared_anchor",
                "basis": (
                    "10.1021/acs.jafc.9b07711 (Bi et al. 2020, key aroma compounds in raw and "
                    "roasted peas) is the correct anchor here and for "
                    "external_validation_bi_2020_raw_pea_hexanal / "
                    "external_validation_bi_2020_roasted_pea_hexanal / intake id "
                    "acs_2020_raw_pea_hexanal_baseline. It was ALSO reused in "
                    "data/lit/benchmark_intake_registry.json for the unrelated EGCG/deoxyosone "
                    "claim (jafc_2020_egcg_deoxyosone_trapping); that anchor has been repointed "
                    "to 10.1021/acs.jafc.0c05098, so this DOI is no longer double-counted across "
                    "unrelated claims. 2026-08-27 (Wave I): this note used to live only in the "
                    "generated benchmark file, where a regeneration deleted it -- it is now "
                    "carried in the bundle spec so it survives. See doi_repair_reversal on the "
                    "bi_2020 anchors in data/lit/flavor_reference_payloads.json."
                ),
            },
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
        # 2026-08-27 (Wave R): was "native commercial pea protein ISOLATE aqueous slurry".
        # The source's quantified panel is 5 isolates + 4 concentrates (Table 2.3), so
        # calling it isolate-only overstated the match to the `pea_iso` lane.
        "matrix_format": (
            "commercial pea protein rehydrated to 10% solids (w/w) in deionized water; external "
            "hold-out. The 9 quantified samples of Table 2.7 are 5 ISOLATES (4, 12, 13, 14, 16) "
            "and 4 CONCENTRATES (1, 2, 19, 22). protein_type stays `pea_iso` because that is the "
            "executable lane this bundle is proxied onto, not because the panel is isolate-only."
        ),
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
            "notes": (
                "Lot-to-lot pea-protein off-note baseline from Liu (2021 thesis / 2023 paper). "
                "Conditions are proxied onto the existing ambient pea-isolate matrix-only lane so "
                "the bundle stays executable but external-only. 2026-08-27 (Wave R): the proxied "
                "pH 6.0 is NOT the source pH -- Table 2.3 reports the rehydrated 10%-solids "
                "slurries at pH 6.3-7.3 (mean ~6.8) across the nine quantified samples. The proxy "
                "is left unchanged because editing an executable condition would change the "
                "prediction, and Wave R corrects reference data only. [P] for a later wave."
            ),
        },
        "analytical_context": {
            # 2026-08-27 (Wave R): retyped against the primary source, which was read in full.
            "headspace_method": (
                "HS-SPME-GC-MS/MS (triple quadrupole, MRM) for both scored rows; the same thesis "
                "also ran DSE-SAFE-GC-MS, GC-O and AEDA, which are NOT the source of these numbers"
            ),
            # Was `internal_standard_calibrated`. There ARE internal standards
            # (2-methyl-3-heptanone, 2-ethoxy-3-isopropyl pyrazine, heptanal-d14), but
            # quantification is by a five-point EXTERNAL standard curve whose response factor was
            # measured in DI WATER rather than in the protein matrix, so matrix binding is
            # uncorrected and the reported concentrations are LOWER BOUNDS on total analyte.
            "quantification_mode": "external_standard_curve_water_calibrated",
            # Was 6. The Table 2.7 footnote reads "Average concentration +/- standard deviation
            # (n = 2)" -- duplicate lots. The triplicate extractions are within-lot and do not
            # make six independent observations.
            "replicates": 2,
            "panel_composition": (
                "Table 2.7 quantifies 9 of the 24 commercial pea proteins: samples 4, 12, 13, 14, "
                "16 (isolates) and 1, 2, 19, 22 (concentrates), per Table 2.3."
            ),
            "notes": (
                "2026-08-27 (Wave R): BOTH scored rows were content-corrected against the primary "
                "source after the previous values were found to match nothing in it. See the "
                "per-row `content_correction_note`. Band midpoints are used as hold-out reference "
                "values; the uncertainty column is derived from the band width, and a geometric "
                "midpoint of an 18-21x band is a construction of ours, not a measurement."
            ),
        },
        # 2026-08-27 (Wave I): typed identifier carried in the SPEC so a regeneration
        # cannot drop it. 2026-08-27 (Wave R): the thesis still has no DOI, but the
        # peer-reviewed version of the same dataset does, and it now flows into
        # `source_doi` from the anchors' `doi` field.
        "provenance_overrides": {
            "identifier": (
                "Liu, Y. (2021), \"Flavor Chemistry of Pea Proteins\", MS thesis, North Carolina "
                "State University Institutional Repository (item "
                "db647868-5ffe-4621-9f11-bbc4db357406)"
            ),
            "identifier_scheme": "citation",
            "identifier_note": (
                "2026-08-27 (Wave I): retyped from source_doi 'LiuThesis2023', which was never a "
                "DOI. 2026-08-27 (Wave R): the thesis carries no DOI, so the typed `identifier` "
                "stays a citation -- and the YEAR is corrected from 2023 to 2021 (the 2023 date "
                "belongs to the derived Food Chemistry paper, not the thesis). `source_doi` now "
                "carries 10.1016/j.foodchem.2022.134998, the DOI of the PEER-REVIEWED VERSION OF "
                "THE SAME DATASET (Liu, Cadwallader & Drake 2023, Food Chemistry 406:134998, "
                "CrossRef-verified). This REVERSES the Wave I note that said the DOI was "
                "'deliberately NOT used here': that reasoning confused CITING the paper with "
                "PROMOTING the bundle. The bundle stays evidence_class external_validation_only, "
                "stays out of every fit, and its classification is unchanged by carrying the "
                "correct DOI."
            ),
            "content_correction_note": (
                "2026-08-27 (Wave R): BOTH scored rows of this bundle were content-corrected. The "
                "prior citation string 'Liu, Y. (2023 thesis)' was doubly wrong -- the thesis is "
                "2021, and the numbers attached to it (hexanal 15-180 ppb, nonanal 5-50 ppb) "
                "matched nothing in either the thesis or the paper. The full thesis text was read "
                "for this wave; the corrected values come from Table 2.7. Two further anchors on "
                "the same source, liu_2023_ppi_heptadienal_band and liu_2023_ppi_ibmp_band, were "
                "marked no_verifiable_source because their compounds are ABSENT from the source; "
                "neither was ever part of this scored bundle."
            ),
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
    # 2026-08-27 (Wave I). What the scored number IS. See VALUE_PROVENANCE_NOTES.
    value_provenance: str = "reported_point_value"
    value_provenance_note: str = ""

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
            "value_provenance": self.value_provenance,
            "value_provenance_note": self.value_provenance_note,
        }


# 2026-08-27 (Wave I) — hold-out value provenance.
#
# The cold-start red team found that hold-out "measurements" were being reported with an
# air of measurement they had not earned: 17-significant-figure values next to an invented
# `measurement_date`, and one value that the repo had computed from its own constants.
# Nothing here is fabricated -- the underlying anchors are real -- but a scored hold-out
# point must say which of these it is, because they support very different claims.
VALUE_PROVENANCE_NOTES = {
    "reported_point_value": (
        "A single concentration reported by the source. Scoring this against a prediction "
        "is a genuine test."
    ),
    "band_geometric_midpoint": (
        "The source reports a RANGE, not a value. The scored number is the geometric "
        "midpoint sqrt(min*max) of that range -- a construction of ours, not a measurement. "
        "The honest uncertainty is the band itself, which is typically a factor of several; "
        "a prediction landing inside such a band is weak evidence."
    ),
    "derived_from_oav_and_repo_threshold": (
        "The scored number is NOT independently measured: it is an odour-activity value "
        "from the source multiplied by an odour-detection threshold taken from THIS REPO. "
        "It therefore partly encodes one of our own constants, and it moves if we change "
        "that constant. Treat it as a consistency check, never as an external measurement."
    ),
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
        if path.startswith(f"{data_paths.rel(BENCHMARK_DIR)}/") and "/external_validation/" not in path:
            return True
    return False


def _bundle_spec_for_anchor(anchor_id: str) -> Optional[Dict[str, Any]]:
    for spec in _HOLDOUT_BUNDLE_SPECS:
        if anchor_id in spec["anchor_ids"]:
            return dict(spec)
    return None


def _publication_year(citation: str) -> str:
    """The source's publication year as a bare year string, or "unspecified"."""
    for part in str(citation).split():
        cleaned = "".join(ch for ch in part if ch.isdigit())
        if len(cleaned) >= 4 and cleaned[:4].isdigit():
            return cleaned[:4]
    return "unspecified"


def _measurement_row(candidate: ExternalCandidate) -> Dict[str, Any]:
    """Build one `measured_volatiles` row, labelled with what the number actually is.

    2026-08-27 (Wave I). Three changes, all forced by the cold-start red team:

    1.  The value is rounded to 4 significant figures. It used to be written at full
        float precision -- 51.96152422706632 for a number that is sqrt(15 * 180) -- which
        advertised a measurement precision that does not exist anywhere in the chain.
    2.  Every row carries `value_provenance` (and its note) so a reader cannot mistake a
        band midpoint or a repo-derived number for a reported measurement.
    3.  For band-derived rows the band itself is carried through as
        `band_min_ppb`/`band_max_ppb`. `uncertainty_pct` is kept for schema compatibility
        but is explicitly labelled: it is the UPWARD half-width of a geometric band, not a
        symmetric analytical uncertainty, and the band is asymmetric in linear space by
        construction.

    2026-08-27 (Wave R). The Wave I 4-significant-figure rounding was applied to EVERY
    number, including ones the source actually reports. That is the opposite of the honesty
    it was written for: rounding sqrt(15 * 180) removes precision we invented, but rounding
    a `reported_point_value` removes precision the SOURCE has. Two measured points were
    being silently altered on every regeneration -- Li 2026 2-pentylfuran 5625.8 -> 5626.0
    and (once Wave R corrected it) the Liu hexanal band max 52454 -> 52450 -- which is also
    why the committed li_2026 artifacts had drifted away from their own generator. Rounding
    is now applied ONLY to values this repo CONSTRUCTED (band midpoints and OAV-derived
    points), never to reported values and never to the band endpoints, which are per-lot
    measurements quoted verbatim from the source table.
    """
    if candidate.point_value is None:
        raise ValueError(f"Candidate {candidate.anchor_id} has no numeric point value")

    is_reported = str(candidate.value_provenance) == "reported_point_value"
    row: Dict[str, Any] = {
        "conc_ppb": float(candidate.point_value)
        if is_reported
        else _round_sig(float(candidate.point_value), 4)
    }
    row["value_provenance"] = candidate.value_provenance
    note = candidate.value_provenance_note or VALUE_PROVENANCE_NOTES.get(
        candidate.value_provenance, ""
    )
    if note:
        row["value_provenance_note"] = note

    if (
        candidate.band_low is not None
        and candidate.band_high is not None
        and candidate.point_value > 0.0
    ):
        # 2026-08-27 (Wave R): NOT rounded. These are per-lot values quoted from the source
        # table (Liu 2021 Table 2.7 hexanal max is 52454, not 52450); only the midpoint
        # between them is a construction of ours.
        row["band_min_ppb"] = float(candidate.band_low)
        row["band_max_ppb"] = float(candidate.band_high)
        row["band_span_fold"] = _round_sig(
            float(candidate.band_high) / float(candidate.band_low), 3
        ) if candidate.band_low > 0 else None
        half_span = max(
            abs(candidate.point_value - candidate.band_low),
            abs(candidate.band_high - candidate.point_value),
        )
        row["uncertainty_pct"] = round((half_span / candidate.point_value) * 100.0, 1)
        row["uncertainty_pct_basis"] = (
            "upward half-width of the reported band relative to its geometric midpoint; "
            "NOT a symmetric analytical uncertainty. Read band_min_ppb/band_max_ppb instead."
        )
    return row


def _round_sig(value: float, digits: int) -> float:
    """Round to `digits` significant figures. 2026-08-27 (Wave I)."""
    if value == 0.0 or not math.isfinite(value):
        return float(value)
    exponent = math.floor(math.log10(abs(value)))
    return float(round(value, -(exponent - digits + 1)))


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

            # 2026-08-27 (Wave I). An anchor may DECLARE how its number was obtained
            # (`value_provenance`); otherwise it is inferred from the numeric block shape.
            # A band's scored value is a geometric midpoint we constructed, and saying so
            # is the whole point.
            declared_provenance = str(anchor.get("value_provenance", "") or "").strip()
            if declared_provenance:
                value_provenance = declared_provenance
            elif lo is not None and hi is not None:
                value_provenance = "band_geometric_midpoint"
            else:
                value_provenance = "reported_point_value"
            value_provenance_note = str(anchor.get("value_provenance_note", "") or "").strip()

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
                    f"{data_paths.rel(BENCHMARK_DIR)} ({_TOP_LEVEL_BENCHMARK_EQUIVALENTS.get(anchor_id, 'registry-backed benchmark')}); "
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
                    value_provenance=value_provenance,
                    value_provenance_note=value_provenance_note,
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
            # 2026-08-27 (Wave R). A row whose value was CORRECTED against the primary source
            # must carry the correction with it into the generated bundle, or the only place
            # the old-vs-new record survives is the ledger -- and a regeneration would then
            # silently present the corrected number as if it had always been right.
            correction_note = str(anchor_row.get("content_correction_note", "") or "").strip()
            if correction_note:
                measured_volatiles[candidate.compound]["content_correction_note"] = correction_note
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
                # 2026-08-27 (Wave I): was `_publication_year_proxy(...)`, which wrote
                # "<year>-01-01" -- an invented date that made a synthesized bundle read
                # like a dated lab record. There is no measurement date for these bundles.
                # 2026-08-27 (Wave R): that function is now DELETED, not just unused. Wave I
                # left it in place "so importers do not break"; a repo-wide search found no
                # importer, and a retired date-fabricator sitting in the module is a loaded
                # gun. Its retirement rationale is preserved in this comment.
                "measurement_date": "not_applicable",
                "measurement_date_note": (
                    "This bundle is SYNTHESIZED from published anchors; no measurement date "
                    "exists. Until 2026-08-27 this field carried a fabricated "
                    "'<publication year>-01-01'."
                ),
                "source_publication_year": _publication_year(primary_citation),
                "notes": (
                    "External validation hold-out synthesized from flavor_reference_payloads.json; "
                    "never use for calibration or promotion. "
                    + " ".join(notes)
                ).strip(),
                # 2026-08-27 (Wave I). Typed identifiers and hand-written audit annotations
                # live in the bundle SPEC, not only in the generated file, so a regeneration
                # cannot silently delete them. This closes the "materialize_external_validation
                # typed-identifier drift" item recorded as FOUND-NOT-FIXED at the end of
                # Wave G3.
                **dict(spec.get("provenance_overrides", {}) or {}),
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


# The committed hold-out bundles are FROZEN EVIDENCE, not a build product. They started
# life as `_HOLDOUT_BUNDLE_SPECS` renderings, but have since been corrected by hand from
# primary sources (Wave W: Bi 2020 Table 4 values and uncertainties; the
# `conditions.buffer` blocks read from the PDFs by
# scripts/generators/complete_benchmark_buffer_fields.py). None of that lives in the spec,
# so a naive regeneration silently reverted it -- which is how `citation_audit_note` was
# lost once (see the spec comment above). Rules, enforced by
# tests/unit/test_external_validation_inventory.py:
#   * by default an existing bundle file is never rewritten;
#   * with ``overwrite=True`` the regenerated payload is written, but every key the
#     committed file has and the regenerated payload lacks is carried forward.
def _fill_missing_from_existing(payload: Any, existing: Any) -> Any:
    """Deep-merge: keys present in ``existing`` but absent from ``payload`` are copied over."""
    if not (isinstance(payload, dict) and isinstance(existing, dict)):
        return payload
    merged = dict(payload)
    for key, value in existing.items():
        if key not in merged:
            merged[key] = copy.deepcopy(value)
        else:
            merged[key] = _fill_missing_from_existing(merged[key], value)
    return merged


def _load_existing_json(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def write_holdout_bundles(
    bundles: Sequence[HoldoutBundle],
    *,
    protocol_dir: Path = EXTERNAL_VALIDATION_PROTOCOL_DIR,
    benchmark_dir: Path = EXTERNAL_VALIDATION_BENCHMARK_DIR,
    overwrite: bool = False,
) -> Dict[str, List[Path]]:
    """Materialize bundles to disk. Existing files are frozen unless ``overwrite``.

    Returns the paths written under ``protocols`` / ``benchmarks`` and the paths left
    untouched under ``skipped``.
    """
    protocol_dir.mkdir(parents=True, exist_ok=True)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    protocol_paths: List[Path] = []
    benchmark_paths: List[Path] = []
    skipped: List[Path] = []
    for bundle in bundles:
        protocol_path = protocol_dir / f"{bundle.bundle_id}.yaml"
        benchmark_path = benchmark_dir / f"{bundle.bundle_id}.json"

        if protocol_path.exists() and not overwrite:
            skipped.append(protocol_path)
        else:
            protocol_path.write_text(
                yaml.safe_dump(bundle.intake_payload, sort_keys=False, allow_unicode=False),
                encoding="utf-8",
            )
            protocol_paths.append(protocol_path)

        existing = _load_existing_json(benchmark_path)
        if existing is not None and not overwrite:
            skipped.append(benchmark_path)
            continue
        benchmark_payload = bundle.benchmark_payload
        if existing is not None:
            benchmark_payload = _fill_missing_from_existing(benchmark_payload, existing)
        benchmark_path.write_text(
            json.dumps(benchmark_payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        benchmark_paths.append(benchmark_path)
    return {"protocols": protocol_paths, "benchmarks": benchmark_paths, "skipped": skipped}


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
            # 2026-08-27 (Wave I). Carry the hold-out point's value provenance through to
            # the report. Half of these eight "measurements" are not measurements: two are
            # geometric midpoints of 10-12x bands, and two are the source's OAV multiplied
            # by this repo's own hexanal odour threshold. A hold-out score means different
            # things for each, and the report has to say which it is scoring.
            provenance = _holdout_value_provenance(
                Path(str(benchmark.get("bench_file", ""))), str(compound.get("compound", ""))
            )
            compound_rows.append(
                {
                    **dict(compound),
                    "abs_log10_error": abs_log10_error,
                    "fold_error": fold_error,
                    "value_provenance": provenance,
                    "is_direct_measurement": provenance == "reported_point_value",
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

    # 2026-08-27 (Wave I). Two honesty blocks the red team asked for.
    #
    # (a2) GENUINE EXTRAPOLATION vs RE-SCORING. A bundle whose executable conditions are
    #      COPIED from an in-panel benchmark is not testing extrapolation at all: it
    #      re-scores that anchor's prediction at that anchor's own conditions. Only bundles
    #      at a genuinely new process state (roasting, HME extrusion) test out-of-envelope
    #      transfer. This was stated in prose before; it is computed here so the headline
    #      cannot drift from it.
    rescoring_rows: List[Dict[str, Any]] = []
    extrapolation_rows: List[Dict[str, Any]] = []
    extrapolation_ids: List[str] = []
    for benchmark in benchmark_rows:
        is_rescoring = _holdout_is_condition_copy(Path(str(benchmark.get("bench_file", ""))))
        benchmark["holdout_kind"] = "in_panel_rescoring" if is_rescoring else "genuine_extrapolation"
        if is_rescoring:
            rescoring_rows.extend(benchmark.get("compounds", []) or [])
        else:
            extrapolation_rows.extend(benchmark.get("compounds", []) or [])
            extrapolation_ids.append(str(benchmark.get("benchmark_id")))

    # (a) MEASUREMENT vs DERIVED. Four of the eight hold-out points are not measurements:
    #     two Liu bands whose scored value is a geometric midpoint we constructed, and two
    #     Bi points computed as the source's OAV times this repo's own hexanal odour
    #     threshold. The latter partly encode one of our own constants and MOVE if that
    #     constant is corrected, so they can only ever be a consistency check.
    provenance_counts: Dict[str, int] = {}
    direct_rows: List[Dict[str, Any]] = []
    derived_rows: List[Dict[str, Any]] = []
    for benchmark in benchmark_rows:
        for compound in benchmark.get("compounds", []) or []:
            key = str(compound.get("value_provenance", "reported_point_value"))
            provenance_counts[key] = provenance_counts.get(key, 0) + 1
            (direct_rows if compound.get("is_direct_measurement") else derived_rows).append(compound)

    def _coverage(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
        total = len(rows)
        hits = sum(1 for row in rows if row.get("inside_ci"))
        return {"hits": hits, "total": total, "rate": (hits / total) if total else None}

    summary["holdout_kind_split"] = {
        "genuine_extrapolation": {
            **_coverage(extrapolation_rows),
            "benchmarks": extrapolation_ids,
            "meaning": (
                "Bundles at a process state the calibration panel does not contain "
                "(roasting, HME extrusion). These are the only rows that test transfer."
            ),
        },
        "in_panel_rescoring": {
            **_coverage(rescoring_rows),
            "meaning": (
                "Bundles whose executable conditions are COPIED from an in-panel benchmark. "
                "Scoring them re-runs that anchor's prediction at that anchor's own "
                "conditions; it is a reproducibility comparison, not an extrapolation test."
            ),
        },
    }

    summary["value_provenance_split"] = {
        "counts": provenance_counts,
        "direct_measurement": _coverage(direct_rows),
        "derived_or_constructed": _coverage(derived_rows),
        "note": (
            "`band_geometric_midpoint` rows score against a number WE constructed as "
            "sqrt(min*max) of a reported range; `derived_from_oav_and_repo_threshold` rows "
            "score against the source's OAV multiplied by this repo's own (compilation-"
            "level, primary-table-unverified) hexanal odour threshold of 4.5 ppb. Only "
            "`reported_point_value` rows are external measurements in the ordinary sense."
        ),
    }

    # (b) PRE-WIDENING COVERAGE. The uncalibrated matrix sigma was raised 2.0 -> 2.86 on
    #     2026-08-26. Any coverage improvement since then is a WIDER INTERVAL, not a better
    #     prediction. Re-run the same hold-out at the pre-widening sigma so the report can
    #     lead with the number that is not an artifact of the prior. Same seed, same
    #     samples, only the matrix prior differs.
    narrow_summary: Dict[str, Any] = {}
    try:
        narrow_payload = propagate_benchmarks(
            benchmark_files=holdout_files,
            n_samples=n_samples,
            seed=seed,
            priors=default_priors(matrix_tier="uncalibrated", matrix_sigma_override=2.0),
            target_tag=target_tag,
            execution_paths=("free_precursor", "matrix_precursor_augmented", "matrix_only"),
        )
    except TypeError:
        narrow_payload = None
    if narrow_payload is not None:
        narrow_hits = 0
        narrow_total = 0
        narrow_extrap_hits = 0
        narrow_extrap_total = 0
        for benchmark in narrow_payload.get("benchmarks", []):
            is_extrapolation = str(benchmark.get("benchmark_id")) in set(extrapolation_ids)
            for compound in benchmark.get("compounds", []):
                narrow_total += 1
                narrow_hits += int(bool(compound.get("inside_ci")))
                if is_extrapolation:
                    narrow_extrap_total += 1
                    narrow_extrap_hits += int(bool(compound.get("inside_ci")))
        narrow_summary = {
            "matrix_sigma": 2.0,
            "ci_coverage_hits": narrow_hits,
            "matched_compound_count": narrow_total,
            "ci_coverage_rate": (narrow_hits / narrow_total) if narrow_total else None,
            "genuine_extrapolation_hits": narrow_extrap_hits,
            "genuine_extrapolation_total": narrow_extrap_total,
            "why": (
                "The uncalibrated matrix ln-sigma was raised 2.0 -> "
                f"{_UNCALIBRATED_MATRIX_SIGMA:.2f} on 2026-08-26. Coverage at the SHIPPED "
                "sigma is therefore not comparable with any earlier run of this report. "
                "This row re-scores the identical hold-out at the pre-widening sigma so "
                "the difference between 'the model got better' and 'the interval got "
                "wider' is visible rather than inferred. Read the median fold error, which "
                "no prior can change, alongside both."
            ),
        }

    summary.update(
        {
            "holdout_benchmark_count": len(holdout_files),
            "median_abs_log10_error": median_abs_log10_error,
            "median_accuracy_fold": median_accuracy_fold,
            "max_fold_error": max(fold_errors) if fold_errors else None,
            "evidence_class": EXTERNAL_VALIDATION_EVIDENCE_CLASS,
            "benchmark_dir": str(EXTERNAL_VALIDATION_BENCHMARK_DIR),
            "pre_widening_coverage": narrow_summary,
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


def _render_holdout_honest_lead(summary: Mapping[str, Any]) -> str:
    """Lead the hold-out report with the pre-widening number. 2026-08-27 (Wave I)."""
    narrow = summary.get("pre_widening_coverage") or {}
    if not narrow:
        return (
            "**Headline external trust metric**: unavailable — this artifact predates the "
            "Wave I pre-widening comparison. Regenerate it."
        )
    hits = narrow.get("ci_coverage_hits", 0)
    total = narrow.get("matched_compound_count", 0)
    shipped_hits = summary.get("ci_coverage_hits", 0)
    shipped_total = summary.get("matched_compound_count", 0)

    extrap_hits = narrow.get("genuine_extrapolation_hits")
    extrap_total = narrow.get("genuine_extrapolation_total")
    kinds = summary.get("holdout_kind_split", {}) or {}
    shipped_extrap = kinds.get("genuine_extrapolation", {}) or {}
    rescoring = kinds.get("in_panel_rescoring", {}) or {}

    lead = (
        "**Headline external trust metric — GENUINE EXTRAPOLATIONS ONLY, at the "
        "PRE-WIDENING prior (matrix ln-sigma 2.0)**: "
        f"**{extrap_hits} / {extrap_total}**.\n\n"
        if extrap_total
        else "**Headline external trust metric**: no genuine-extrapolation rows.\n\n"
    )
    return (
        lead
        + "Everything else on this page is a weaker claim than that number, in two separate "
        "ways, and both were being read as if they were not:\n\n"
        f"1. **Prior width.** At the currently shipped sigma "
        f"({_UNCALIBRATED_MATRIX_SIGMA:.2f}, ~±{_UNCALIBRATED_MATRIX_FOLD:.0f}× at 90%) the "
        f"*same predictions* score {shipped_extrap.get('hits', 0)} / "
        f"{shipped_extrap.get('total', 0)} on those same extrapolation rows, and "
        f"{shipped_hits} / {shipped_total} over all hold-out rows (vs {hits} / {total} at "
        "ln-sigma 2.0). **Nothing about the model changed between those numbers — only the "
        "width of the interval drawn around it.**\n"
        f"2. **Which rows test anything.** {rescoring.get('total', 0)} of the "
        f"{shipped_total} rows come from bundles whose executable conditions are *copied "
        "from an in-panel benchmark*: scoring them re-runs an existing anchor at its own "
        f"conditions. They score {rescoring.get('hits', 0)} / {rescoring.get('total', 0)}, "
        "and that is a reproducibility comparison, not evidence of transfer.\n"
    )


def _render_holdout_provenance_split(split: Mapping[str, Any]) -> str:
    """What the eight hold-out 'measurements' actually are. 2026-08-27 (Wave I)."""
    if not split:
        return ""
    counts = split.get("counts", {}) or {}
    direct = split.get("direct_measurement", {}) or {}
    derived = split.get("derived_or_constructed", {}) or {}

    def _cell(bucket: Mapping[str, Any]) -> str:
        total = bucket.get("total", 0)
        hits = bucket.get("hits", 0)
        rate = bucket.get("rate")
        return f"{hits}/{total}" + (f" ({rate * 100:.0f}%)" if rate is not None else "")

    lines = [
        "> **What is actually being scored here.** These are not eight external",
        "> measurements. The hold-out points divide as follows:",
        ">",
        "> | Provenance | Rows | Inside 90% CI | What a score against it means |",
        "> | --- | ---: | ---: | --- |",
    ]
    descriptions = {
        "reported_point_value": (
            "A concentration the source reports. A genuine external test."
        ),
        "band_geometric_midpoint": (
            # 2026-08-27 (Wave R): the parenthetical read "(10-12x end to end)", which was
            # the span of the two Liu bands BEFORE they were corrected against Table 2.7 of
            # their own source. The real inter-lot spans are 21.45x (hexanal) and 18.19x
            # (nonanal). Each row's own `band_span_fold` is authoritative; this prose no
            # longer hard-codes a number it cannot keep in step with.
            "The source reports a RANGE; the scored value is sqrt(min*max), a number we "
            "constructed. The honest uncertainty is the band itself -- read each row's "
            "band_min_ppb/band_max_ppb and band_span_fold -- so landing inside it is weak "
            "evidence."
        ),
        "derived_from_oav_and_repo_threshold": (
            "The source's odour-activity value multiplied by THIS REPO'S OWN hexanal "
            "odour threshold (4.5 ppb, compilation-level and never verified against a "
            "primary table). Partly encodes one of our own constants and moves if that "
            "constant is corrected. A consistency check, not an external measurement."
        ),
    }
    for key, description in descriptions.items():
        n = counts.get(key, 0)
        if not n:
            continue
        bucket = direct if key == "reported_point_value" else None
        cell = _cell(bucket) if bucket is not None else "see combined row"
        lines.append(f"> | `{key}` | {n} | {cell} | {description} |")
    lines.extend(
        [
            ">",
            f"> Combined: direct measurements **{_cell(direct)}**, "
            f"derived or constructed **{_cell(derived)}**.",
        ]
    )
    return "\n".join(lines)


def _holdout_is_condition_copy(bench_file: Path) -> bool:
    """True when this bundle's conditions are copied from its in-panel target benchmark.

    2026-08-27 (Wave I). Such a bundle re-scores an existing anchor at that anchor's own
    conditions -- a reproducibility comparison, not an extrapolation test. Computed by
    comparing the four executable condition fields against the benchmark named in the
    bundle spec's `benchmark_alignment.target_benchmark_id`.
    """
    try:
        bench = json.loads(Path(bench_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    bundle_id = str(bench.get("benchmark_id", ""))
    spec = next(
        (item for item in _HOLDOUT_BUNDLE_SPECS if item["bundle_id"] == bundle_id),
        None,
    )
    if spec is None:
        return False
    target_id = str((spec.get("benchmark_alignment") or {}).get("target_benchmark_id", "")).strip()
    if not target_id:
        return False
    target_path = BENCHMARK_DIR / f"{target_id}.json"
    if not target_path.exists():
        return False
    try:
        target = json.loads(target_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    fields = ("temp_C", "ph", "water_activity", "time_min")
    here = bench.get("conditions", {}) or {}
    there = target.get("conditions", {}) or {}
    return all(
        abs(float(here.get(field, float("nan"))) - float(there.get(field, float("inf")))) < 1e-9
        for field in fields
        if field in here and field in there
    ) and all(field in here and field in there for field in fields)


def _holdout_value_provenance(bench_file: Path, compound: str) -> str:
    """`value_provenance` recorded on this hold-out row. 2026-08-27 (Wave I)."""
    try:
        bench = json.loads(Path(bench_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "unknown"
    measured = bench.get("measured_volatiles", {}) or {}
    row = measured.get(compound)
    if row is None:
        wanted = str(compound).strip().lower()
        for key, value in measured.items():
            if str(key).strip().lower() == wanted:
                row = value
                break
    if not isinstance(row, Mapping):
        return "unknown"
    return str(row.get("value_provenance", "reported_point_value"))


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
        # 2026-08-27 (Wave I): the honest lead. Two things were being read out of this
        # report that it does not support: that eight external MEASUREMENTS were being
        # scored (four are constructed or repo-derived), and that coverage had improved
        # (the interval was widened 2.0 -> 2.86 on 2026-08-26). Both are stated first now.
        _render_holdout_honest_lead(summary),
        "",
        f"**Median accuracy on hold-outs**: **{median_accuracy_str}** median fold error (median |log10 error| = **{median_abs_str}** dex). *This number is unaffected by the prior width and is the one to track across runs.*",
        "",
        _render_holdout_provenance_split(summary.get("value_provenance_split") or {}),
        "",
        f"_Secondary, prior-dependent figure: at the SHIPPED uncalibrated sigma the measured value lies inside the 90% CI for {summary.get('ci_coverage_hits', 0)} / {summary.get('matched_compound_count', 0)} matched hold-out compounds ({coverage_str})._",
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
        # 2026-08-27 (Wave I): `Reference value` replaces `Measured (ppb)` and carries a
        # provenance column -- half of these are not measurements (see the split above).
        lines.append(
            "| Compound | Reference value (ppb) | Provenance | P5 | P50 | P95 | Fold error | Inside 90% CI |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for compound in benchmark.get("compounds", []):
            inside = "yes" if compound.get("inside_ci") else "no"
            fold_error = compound.get("fold_error")
            fold_error_str = f"{float(fold_error):.2f}x" if fold_error is not None else "n/a"
            provenance = str(compound.get("value_provenance", "unknown"))
            lines.append(
                f"| {compound.get('compound')} | {float(compound.get('measured_ppb', 0.0)):.3g} | "
                f"`{provenance}` | "
                f"{float(compound.get('predicted_p5', 0.0)):.3g} | {float(compound.get('predicted_p50', 0.0)):.3g} | "
                f"{float(compound.get('predicted_p95', 0.0)):.3g} | {fold_error_str} | {inside} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_external_validation_artifact(
    payload: Mapping[str, Any],
    *,
    output_dir: Path | str = data_paths.VALIDATION_DIR,
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
