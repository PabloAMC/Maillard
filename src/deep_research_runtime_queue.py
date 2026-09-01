from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Set

from src import data_paths

RUNTIME_LANE_METADATA = {
    "process_state_calibration": {
        "landing_path": data_paths.rel(data_paths.PROCESS_STATE_CALIBRATIONS),
        "title": "Process-State Calibration",
    },
    "safety_reference": {
        "landing_path": data_paths.rel(data_paths.SAFETY_REFERENCE_PAYLOADS),
        "title": "Safety Reference",
    },
    "computational_prior": {
        "landing_path": data_paths.rel(data_paths.COMPUTATIONAL_PRIORS),
        "title": "Computational Prior",
    },
}

RUNTIME_LANE_ORDER = {
    "process_state_calibration": 0,
    "safety_reference": 1,
    "computational_prior": 2,
}

CURATED_RUNTIME_BATCH_01 = [
    {
        "citation": "ACS JAFC 3c05991 / PMC10739987",
        "runtime_lane": "process_state_calibration",
        "theme": "protein_volatile_partitioning",
        "target_modules": ["matrix_correction", "headspace", "literature_runtime"],
        "why": "Commercial PPI/SPI partition constants plus temperature dependence can land as a volatile-binding calibration layer without pretending benchmark closure.",
        "repo_next_action": "Extract PPI/SPI partition coefficients and temperature slope into a volatile-binding calibration entry; keep it as calibration, not benchmark.",
    },
    {
        "citation": "Malia et al. (2025)",
        "runtime_lane": "process_state_calibration",
        "theme": "pea_accessibility_state",
        "target_modules": ["matrix_correction"],
        "why": "Free-SH data from pea protein is a direct process-state anchor for cysteine accessibility and belongs in the runtime calibration layer.",
        "repo_next_action": "Extract free-SH versus process state into a pea cysteine-accessibility calibration entry.",
    },
    {
        "citation": "Raman/SDS",
        "runtime_lane": "process_state_calibration",
        "theme": "extrusion_disulfide_severity",
        "target_modules": ["matrix_correction", "headspace"],
        "why": "SME-dependent disulfide growth and thiol loss can become an extrusion-severity calibration before any mechanistic retention model is expanded.",
        "repo_next_action": "Encode SME-dependent disulfide growth and Ellman thiol loss as extrusion disulfide-severity calibration metadata.",
    },
    {
        "citation": "PMC PMCID:PMC12648097 (Ref. 5)",
        "runtime_lane": "safety_reference",
        "theme": "acrylamide_mitigation_reference",
        "target_modules": ["safety"],
        "why": "Strong quantitative suppression of acrylamide by amino-acid capping should land as a mitigation-side safety reference, not as a kinetic benchmark.",
        "repo_next_action": "Encode as a mitigation-side safety reference for acrylamide suppression by amino-acid capping and keep it separate from kinetic benchmarking.",
    },
    {
        "citation": "Kocadağlı & Gökmen (2016)",
        "runtime_lane": "safety_reference",
        "theme": "acrylamide_process_modifier",
        "target_modules": ["safety"],
        "why": "The saponin-linked acrylamide and pyrazine shift can land as a bounded process modifier for safety interpretation without pretending full kinetic closure.",
        "repo_next_action": "Encode the saponin-linked acrylamide modifier as an extended safety/process reference with explicit DSC and pyrazine-side caveats.",
    },
    {
        "citation": "ACS JAFC 3c02618",
        "runtime_lane": "computational_prior",
        "theme": "disulfide_retention_prior",
        "target_modules": ["literature_runtime", "headspace", "matrix_correction"],
        "why": "MFT disulfide trapping and pyrazine stacking belong in a bounded retention prior before any broader bidirectional crosstalk model is opened.",
        "repo_next_action": "Encode MFT trapping and pyrazine stacking as a bounded retention prior with explicit uncertainty and matrix applicability.",
    },
    {
        "citation": "ACS JAFC 0c01925",
        "runtime_lane": "computational_prior",
        "theme": "protein_volatile_binding_prior",
        "target_modules": ["literature_runtime", "headspace"],
        "why": "The Lys-versus-Cys adduct preference and headspace-recovery hierarchy improve runtime realism but still belong to a bounded prior rather than a benchmark lane.",
        "repo_next_action": "Encode Lys-versus-Cys adduct preference and headspace-recovery hierarchy as a bounded protein-binding prior.",
    },
    {
        "citation": "Nakagawa et al. (2004)",
        "runtime_lane": "computational_prior",
        "theme": "isoflavone_dicarbonyl_sink",
        "target_modules": ["safety", "literature_runtime", "matrix_correction"],
        "why": "Soy isoflavone inhibition of AGE chemistry is a bounded dicarbonyl-pool modifier that should land as a prior, not as endpoint closure.",
        "repo_next_action": "Encode soy-isoflavone inhibition as a bounded dicarbonyl-pool suppression prior with explicit uncertainty posture.",
    },
]


CURATED_RUNTIME_BATCH_02 = [
    {
        "citation": "Rizzello et al. (2024)",
        "runtime_lane": "process_state_calibration",
        "theme": "fermentation_hexanal_cleanup",
        "target_modules": ["literature_runtime", "matrix_correction", "headspace"],
        "why": "LAB-driven hexanal depletion with simultaneous sulfur uplift is a high-value process-state anchor that can calibrate fermentation cleanup without claiming benchmark closure.",
        "repo_next_action": "Extract fermentation-driven hexanal depletion and sulfur uplift multipliers into a process-state calibration entry for fermented plant matrices.",
    },
    {
        "citation": "Zhao et al. (2022)",
        "runtime_lane": "process_state_calibration",
        "theme": "fermentation_precursor_release",
        "target_modules": ["literature_runtime", "matrix_correction"],
        "why": "Koji and moromi amino-nitrogen plus nucleotide release can tighten the upstream donor-accessibility lane before any new benchmark promotion.",
        "repo_next_action": "Encode fermentation-driven precursor-release multipliers for amino-type nitrogen and nucleotide-rich systems as process-state calibration metadata.",
    },
    {
        "citation": "J. Agric. Food Chem. 2019 (Ref. 24)",
        "runtime_lane": "computational_prior",
        "theme": "polyphenol_thiol_capping_prior",
        "target_modules": ["literature_runtime", "safety"],
        "why": "Quantified cysteine adduct burden and thiol loss in the polyphenol lane belong in a bounded prior because they sharpen precursor-sink behavior without claiming endpoint closure.",
        "repo_next_action": "Encode the quantified cysteine adduct and thiol-loss burden as a bounded polyphenol-thiol capping prior for Family 13 runtime lanes.",
    },
    {
        "citation": "Liardon, de Weck-Gaudard & Philippossian (1991)",
        "runtime_lane": "computational_prior",
        "theme": "pentose_potency_prior",
        "target_modules": ["literature_runtime", "family_upstream_contract"],
        "why": "R5P and ribose donor-strength multipliers are directly encodable and immediately improve upstream sugar-routing realism.",
        "repo_next_action": "Encode R5P-versus-glucose donor-strength multipliers plus the ribose reaction-order prior into the reducing-sugar runtime lane.",
    },
    {
        "citation": "Huang et al. (2021)",
        "runtime_lane": "computational_prior",
        "theme": "sulfur_oav_support_prior",
        "target_modules": ["literature_runtime", "reporting"],
        "why": "High-OAV sulfur markers in cysteine-xylose-glutamate systems strengthen bounded sulfur-support expectations without needing a benchmark payload.",
        "repo_next_action": "Encode sulfur-positive OAV support priors for pentose plus sulfur systems, focused on MFT and FFT support rather than absolute closure.",
    },
    {
        "citation": "Blank et al. (2001)",
        "runtime_lane": "computational_prior",
        "theme": "lipid_epoxide_offnote_prior",
        "target_modules": ["literature_runtime", "lipid_oxidation"],
        "why": "C20:4-derived epoxydecenal with an ultra-low odor threshold is a high-leverage oxidative guardrail that belongs in the bounded-prior lane.",
        "repo_next_action": "Encode the C20:4 epoxydecenal off-note prior as a bounded oxidative guardrail for lipid-rich plant matrices.",
    },
]


CURATED_RUNTIME_BATCH_03 = [
    {
        "citation": "Ordoudi et al. (2014 / PMC12484514)",
        "runtime_lane": "process_state_calibration",
        "theme": "caramelization_hmf_peak_window",
        "target_modules": ["literature_runtime", "matrix_correction", "reporting"],
        "why": "The pH- and time-resolved HMF peak plus secondary decline is a direct process-state calibration surface for caramelization severity rather than a benchmark payload.",
        "repo_next_action": "Encode the HMF peak-at-pH-5.0 and secondary-decline window as caramelization process-state calibration metadata with amino-acid co-catalysis bounds.",
    },
    {
        "citation": "Hrncirik & Zeelenberg (2014)",
        "runtime_lane": "process_state_calibration",
        "theme": "coconut_oil_thermal_profile",
        "target_modules": ["lipid_oxidation", "matrix_correction", "reporting"],
        "why": "The coconut-oil thermal profile gives a usable process-state anchor for lipid-rich co-matrices without pretending broad benchmark closure in that family.",
        "repo_next_action": "Extract coconut-oil thermal marker behavior and lactone OAV support into a lipid-rich co-matrix process-state calibration entry.",
    },
    {
        "citation": "Aliani & Farmer (2005)",
        "runtime_lane": "computational_prior",
        "theme": "donor_potency_and_nucleotide_synergy_prior",
        "target_modules": ["literature_runtime", "family_upstream_contract", "reporting"],
        "why": "The ribose, G6P, and glucose donor multipliers plus nucleotide context are directly encodable as bounded upstream priors and fit existing donor-hierarchy lanes.",
        "repo_next_action": "Encode ribose/G6P/glucose potency multipliers and nucleotide-support context as bounded donor-priority and sulfur-support priors.",
    },
    {
        "citation": "Blank, Devaud & Grosch (2003)",
        "runtime_lane": "computational_prior",
        "theme": "phosphorylated_sugar_hdmf_prior",
        "target_modules": ["literature_runtime", "family_upstream_contract", "reporting"],
        "why": "The G6P-to-HDMF uplift is a clean phosphorylated-donor prior that sharpens reducing-sugar routing without opening a new benchmark lane.",
        "repo_next_action": "Encode G6P-versus-glucose HDMF uplift as a bounded phosphorylated-sugar prior in the donor hierarchy runtime lane.",
    },
    {
        "citation": "Glomb & Monnier (1995)",
        "runtime_lane": "computational_prior",
        "theme": "dicarbonyl_fragmentation_stoichiometry_prior",
        "target_modules": ["literature_runtime", "safety", "lipid_oxidation"],
        "why": "The 3-DG fragmentation stoichiometry gives deterministic branching fractions for MGO, GO, and diketones and belongs in the low-cost prior layer.",
        "repo_next_action": "Encode 3-DG retro-aldol branching fractions into a bounded dicarbonyl fragmentation prior for AGE and off-note pathways.",
    },
    {
        "citation": "Hidalgo & Zamora (2004)",
        "runtime_lane": "computational_prior",
        "theme": "4hne_pyrrole_offnote_prior",
        "target_modules": ["literature_runtime", "lipid_oxidation", "reporting"],
        "why": "The quantified 4-HNE plus phenylalanine route to 2-pentylpyrrole is a high-value off-note prior that can land without inventing a new benchmark surface.",
        "repo_next_action": "Encode the 4-HNE-to-2-pentylpyrrole route as a bounded lipid-derived off-note prior with absolute concentration support.",
    },
]


CURATED_RUNTIME_BATCHES = [
    {
        "batch_id": "2026-04-05-runtime-first-batch-01",
        "batch_label": "Batch 01",
        "candidates": CURATED_RUNTIME_BATCH_01,
    },
    {
        "batch_id": "2026-04-05-runtime-first-batch-02",
        "batch_label": "Batch 02",
        "candidates": CURATED_RUNTIME_BATCH_02,
    },
    {
        "batch_id": "2026-04-05-runtime-first-batch-03",
        "batch_label": "Batch 03",
        "candidates": CURATED_RUNTIME_BATCH_03,
    },
]


def _load_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _derive_slr_families(files: List[str]) -> List[str]:
    family_ids = set()
    for filename in files:
        match = re.match(r"^(\d{2})_", str(filename))
        if match:
            family_ids.add(match.group(1))
    return sorted(family_ids)


def _runtime_lane_row(runtime_lane: str) -> Mapping[str, str]:
    return RUNTIME_LANE_METADATA[runtime_lane]


def _normalize_citation(value: str) -> str:
    normalized = str(value).strip().casefold()
    normalized = normalized.replace("`", "")
    normalized = re.sub(r"\s+", " ", normalized)
    return normalized


def _under_root(root: Path, path: Path) -> Path:
    """``path`` (a ``data_paths`` constant) re-rooted under a scratch ``root``."""
    if Path(root).resolve() == data_paths.REPO_ROOT:
        return path
    return root / data_paths.rel(path)


def _iter_registry_citations(root: Path) -> Iterable[str]:
    process_payload = _load_json(_under_root(root, data_paths.PROCESS_STATE_CALIBRATIONS))
    for entry in process_payload.get("entries", []):
        if isinstance(entry, Mapping) and str(entry.get("source_citation", "")).strip():
            yield str(entry.get("source_citation", ""))

    safety_payload = _load_json(_under_root(root, data_paths.SAFETY_REFERENCE_PAYLOADS))
    for entry in safety_payload.get("entries", []):
        if isinstance(entry, Mapping) and str(entry.get("source_citation", "")).strip():
            yield str(entry.get("source_citation", ""))

    prior_payload = _load_json(_under_root(root, data_paths.COMPUTATIONAL_PRIORS))
    for section_name, entries in prior_payload.items():
        if section_name == "section_family_metadata" or not isinstance(entries, list):
            continue
        for entry in entries:
            if isinstance(entry, Mapping) and str(entry.get("source", entry.get("source_citation", ""))).strip():
                yield str(entry.get("source", entry.get("source_citation", "")))


def _load_landed_runtime_citations(root: Path) -> Set[str]:
    return {
        _normalize_citation(citation)
        for citation in _iter_registry_citations(root)
        if _normalize_citation(citation)
    }


def _citation_is_landed(citation: str, landed_citations: Set[str]) -> bool:
    normalized = _normalize_citation(citation)
    if not normalized:
        return False
    return any(normalized in landed or landed in normalized for landed in landed_citations)


def _select_curated_batch(landed_citations: Set[str]) -> Mapping[str, Any]:
    for batch in CURATED_RUNTIME_BATCHES:
        candidates = batch.get("candidates", []) or []
        if any(not _citation_is_landed(str(spec.get("citation", "")), landed_citations) for spec in candidates):
            return batch
    return CURATED_RUNTIME_BATCHES[-1]


def _next_curated_batch(batch_id: str) -> Mapping[str, Any]:
    for index, batch in enumerate(CURATED_RUNTIME_BATCHES):
        if str(batch.get("batch_id", "")) == str(batch_id):
            if index + 1 < len(CURATED_RUNTIME_BATCHES):
                return CURATED_RUNTIME_BATCHES[index + 1]
            return {}
    return {}


def _selected_candidate(spec: Mapping[str, Any], item: Mapping[str, Any]) -> Dict[str, Any]:
    runtime_lane = str(spec["runtime_lane"])
    lane_metadata = _runtime_lane_row(runtime_lane)
    return {
        "citation": str(item.get("citation", spec["citation"])),
        "runtime_lane": runtime_lane,
        "runtime_lane_title": lane_metadata["title"],
        "landing_path": lane_metadata["landing_path"],
        "score": str(item.get("score", "unknown")),
        "score_value": int(item.get("score_value", 0)),
        "status": str(item.get("status", "unknown")),
        "source_files": [str(value) for value in item.get("files", []) or []],
        "slr_families": _derive_slr_families([str(value) for value in item.get("files", []) or []]),
        "description": str((item.get("descriptions", []) or [""])[0]),
        "theme": str(spec["theme"]),
        "target_modules": [str(value) for value in spec.get("target_modules", []) or []],
        "why": str(spec["why"]),
        "repo_next_action": str(spec["repo_next_action"]),
    }


def _prepared_candidate(spec: Mapping[str, Any], backlog_items: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    citation = str(spec.get("citation", ""))
    item = backlog_items.get(citation, {}) if isinstance(backlog_items.get(citation, {}), Mapping) else {}
    runtime_lane = str(spec.get("runtime_lane", "unknown"))
    return {
        "citation": citation,
        "runtime_lane": runtime_lane,
        "runtime_lane_title": _runtime_lane_row(runtime_lane)["title"] if runtime_lane in RUNTIME_LANE_METADATA else "unknown",
        "status": str(item.get("status", "unknown")),
        "score_value": int(item.get("score_value", 0) or 0),
        "slr_families": _derive_slr_families([str(value) for value in item.get("files", []) or []]),
        "theme": str(spec.get("theme", "unknown")),
        "target_modules": [str(value) for value in spec.get("target_modules", []) or []],
        "repo_next_action": str(spec.get("repo_next_action", "")),
    }


def build_deep_research_runtime_queue(root: Path = data_paths.REPO_ROOT) -> Dict[str, Any]:
    payload = _load_json(_under_root(root, data_paths.DEEP_RESEARCH_BACKLOG))
    backlog_items = {
        str(item.get("citation", "")): item
        for item in payload.get("items", [])
        if isinstance(item, Mapping)
    }
    landed_citations = _load_landed_runtime_citations(root)
    batch = _select_curated_batch(landed_citations)
    next_batch = _next_curated_batch(str(batch.get("batch_id", "")))
    curated_specs: Sequence[Mapping[str, Any]] = batch.get("candidates", []) or []

    selected_candidates: List[Dict[str, Any]] = []
    excluded_candidates: List[Dict[str, Any]] = []
    for spec in curated_specs:
        citation = str(spec["citation"])
        if _citation_is_landed(citation, landed_citations):
            excluded_candidates.append(
                {
                    "citation": citation,
                    "reason": "already_landed_in_runtime_registry",
                }
            )
            continue

        item = backlog_items.get(citation)
        if item is None:
            excluded_candidates.append(
                {
                    "citation": citation,
                    "reason": "not_found_in_backlog",
                }
            )
            continue

        if str(item.get("status", "")) != "BACKLOG":
            excluded_candidates.append(
                {
                    "citation": citation,
                    "reason": f"already_{str(item.get('status', 'unknown')).casefold()}",
                    "registry_id": item.get("registry_id"),
                }
            )
            continue

        selected_candidates.append(_selected_candidate(spec, item))

    selected_candidates.sort(
        key=lambda row: (
            RUNTIME_LANE_ORDER.get(str(row.get("runtime_lane", "")), 9),
            -int(row.get("score_value", 0)),
            str(row.get("citation", "")).casefold(),
        )
    )

    lane_counts: Dict[str, int] = {lane: 0 for lane in RUNTIME_LANE_METADATA}
    themes: Dict[str, int] = {}
    for row in selected_candidates:
        lane = str(row.get("runtime_lane", "unknown"))
        lane_counts[lane] = lane_counts.get(lane, 0) + 1
        theme = str(row.get("theme", "unknown"))
        themes[theme] = themes.get(theme, 0) + 1

    prepared_next_batch_candidates = [
        _prepared_candidate(spec, backlog_items)
        for spec in (next_batch.get("candidates", []) or [])
    ]
    prepared_next_batch_lane_counts: Dict[str, int] = {lane: 0 for lane in RUNTIME_LANE_METADATA}
    for row in prepared_next_batch_candidates:
        lane = str(row.get("runtime_lane", "unknown"))
        prepared_next_batch_lane_counts[lane] = prepared_next_batch_lane_counts.get(lane, 0) + 1

    return {
        "schema_version": "1.0",
        "batch_id": str(batch.get("batch_id", "unknown")),
        "batch_label": str(batch.get("batch_label", "unknown")),
        "source": data_paths.rel(data_paths.DEEP_RESEARCH_BACKLOG),
        "selection_policy": "Select only Deep Research citations still marked BACKLOG and stage them into non-benchmark runtime lanes: process_state_calibration, safety_reference, and computational_prior. Exclude anything already landed in runtime registries, already runtime-bound, or benchmark-first.",
        "selected_candidates": selected_candidates,
        "excluded_candidates": excluded_candidates,
        "prepared_next_batch": {
            "batch_id": str(next_batch.get("batch_id", "")),
            "batch_label": str(next_batch.get("batch_label", "")),
            "candidates": prepared_next_batch_candidates,
            "summary": {
                "candidate_count": len(prepared_next_batch_candidates),
                "process_state_calibration_count": prepared_next_batch_lane_counts.get("process_state_calibration", 0),
                "safety_reference_count": prepared_next_batch_lane_counts.get("safety_reference", 0),
                "computational_prior_count": prepared_next_batch_lane_counts.get("computational_prior", 0),
            },
        },
        "summary": {
            "selected_candidate_count": len(selected_candidates),
            "process_state_calibration_count": lane_counts.get("process_state_calibration", 0),
            "safety_reference_count": lane_counts.get("safety_reference", 0),
            "computational_prior_count": lane_counts.get("computational_prior", 0),
            "landed_runtime_citation_count": len([spec for spec in curated_specs if _citation_is_landed(str(spec.get("citation", "")), landed_citations)]),
            "excluded_candidate_count": len(excluded_candidates),
            "themes": dict(sorted(themes.items())),
        },
    }


def render_deep_research_runtime_queue_markdown(payload: Mapping[str, Any]) -> str:
    summary = payload.get("summary", {}) or {}
    prepared_next_batch = payload.get("prepared_next_batch", {}) or {}
    prepared_next_summary = prepared_next_batch.get("summary", {}) or {}
    lines = [
        "# Deep Research Runtime Queue",
        "",
        f"Batch: {payload.get('batch_id', 'unknown')}",
        f"Source: {payload.get('source', 'unknown')}",
        f"Selection policy: {payload.get('selection_policy', 'unknown')}",
        "",
        f"Selected runtime candidates: {summary.get('selected_candidate_count', 0)}",
        f"Process-state calibrations: {summary.get('process_state_calibration_count', 0)}",
        f"Safety references: {summary.get('safety_reference_count', 0)}",
        f"Computational priors: {summary.get('computational_prior_count', 0)}",
        f"Already landed in runtime registries: {summary.get('landed_runtime_citation_count', 0)}",
        f"Excluded because already encoded or absent: {summary.get('excluded_candidate_count', 0)}",
        "",
        "Out of scope for this batch: benchmark payload promotion.",
    ]

    if prepared_next_batch.get("batch_id"):
        lines.extend([
            "",
            f"Prepared next batch: {prepared_next_batch.get('batch_id', 'unknown')}",
            f"Prepared candidate count: {prepared_next_summary.get('candidate_count', 0)}",
            f"Prepared process-state calibrations: {prepared_next_summary.get('process_state_calibration_count', 0)}",
            f"Prepared safety references: {prepared_next_summary.get('safety_reference_count', 0)}",
            f"Prepared computational priors: {prepared_next_summary.get('computational_prior_count', 0)}",
        ])

    for runtime_lane in ["process_state_calibration", "safety_reference", "computational_prior"]:
        lane_rows = [
            row for row in payload.get("selected_candidates", [])
            if str(row.get("runtime_lane", "")) == runtime_lane
        ]
        lane_metadata = _runtime_lane_row(runtime_lane)
        lines.extend([
            "",
            f"## {lane_metadata['title']}",
            "",
            f"Landing registry: {lane_metadata['landing_path']}",
            "",
            "| Citation | Families | Score | Theme | Target Modules | Next Action |",
            "| --- | --- | ---: | --- | --- | --- |",
        ])
        for row in lane_rows:
            lines.append(
                f"| {row['citation']} | {', '.join(row.get('slr_families', [])) or 'unknown'} | {row['score_value']} | {row['theme']} | {', '.join(row.get('target_modules', [])) or 'none'} | {row['repo_next_action']} |"
            )
        if not lane_rows:
            lines.append("| none | n/a | 0 | n/a | n/a | no curated candidates selected for this lane |")

    if payload.get("excluded_candidates"):
        lines.extend([
            "",
            "## Excluded Candidates",
            "",
            "| Citation | Reason |",
            "| --- | --- |",
        ])
        for row in payload.get("excluded_candidates", []):
            lines.append(f"| {row.get('citation', 'unknown')} | {row.get('reason', 'unknown')} |")

    if prepared_next_batch.get("candidates"):
        lines.extend([
            "",
            "## Prepared Next Batch",
            "",
            "| Citation | Families | Status | Lane | Theme | Next Action |",
            "| --- | --- | --- | --- | --- | --- |",
        ])
        for row in prepared_next_batch.get("candidates", []):
            lines.append(
                f"| {row.get('citation', 'unknown')} | {', '.join(row.get('slr_families', [])) or 'unknown'} | {row.get('status', 'unknown')} | {row.get('runtime_lane', 'unknown')} | {row.get('theme', 'unknown')} | {row.get('repo_next_action', 'none')} |"
            )

    return "\n".join(lines) + "\n"