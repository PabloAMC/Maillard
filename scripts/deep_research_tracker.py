#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


ROOT = Path(__file__).resolve().parents[1]
DEEP_RESEARCH_DIR = ROOT / "data" / "Gemini_Deep_Research"
REGISTRY_FILE = ROOT / "data" / "lit" / "benchmark_intake_registry.json"
OUTPUT_JSON = ROOT / "data" / "lit" / "deep_research_backlog.json"
OUTPUT_MD = ROOT / "reports" / "deep_research_gap_analysis.md"
NUMBERED_REPORT_PATTERN = re.compile(r"^\d{2}_.+\.md$")
ENTRY_PATTERN = re.compile(
    r"^\d+\.\s+(?:`(?P<citation_bt>[^`]+)`|(?P<citation_plain>[^—-]+?))\s*[-—]\s*(?P<description>.*?)\s*(?P<score>[1-8]/8)(?:\.\s*(?P<trailing>.*))?$",
    re.IGNORECASE,
)

PRIORITY_53_TARGETS = {
    "spi_hvp_xylose": ["pmc9905368", "spi-hvp + xylose", "mft oav 450", "fft oav 84"],
    "wheat_gluten_hvp_xylose": ["wheat gluten hvp", "mft oav 850"],
    "acrylamide_fast_kinetics": ["22.36", "62.62", "130°c", "130c"],
    "cml_cel_commercial_pba_ranges": ["foods 12(10):1967", "cml", "cel", "commercial pbma"],
    "furosine_crossover": ["8.7 mg/100 g", "furosine", "fall above 150", "degradation dominates"],
}


def _normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", str(value).casefold())).strip()


def _author_year_key(citation: str) -> str:
    lowered = str(citation).casefold().replace(" et al.", "").replace(" & ", " ")
    match = re.match(r"(?P<author>.+?)\s*\((?P<year>\d{4})\)", lowered)
    if not match:
        return _normalize_text(citation)
    author = _normalize_text(match.group("author")).split(" ")[0]
    return f"{author} {match.group('year')}".strip()


def parse_score(score: str) -> int:
    return int(str(score).split("/", 1)[0])


def iter_markdown_paths(directory: Path = DEEP_RESEARCH_DIR) -> Iterable[Path]:
    for path in sorted(directory.iterdir()):
        if path.is_file() and NUMBERED_REPORT_PATTERN.match(path.name):
            yield path


def load_registry_entries(
    registry_file: Path = REGISTRY_FILE,
    *,
    root: Path = ROOT,
) -> List[Dict[str, Any]]:
    if not registry_file.exists():
        return []
    payload = json.loads(registry_file.read_text(encoding="utf-8"))
    entries: List[Dict[str, Any]] = []
    for row in payload.get("eligible_references", []):
        citation = str(row.get("citation", "")).strip()
        runtime_artifacts = list(row.get("runtime_artifacts", []) or [])
        runtime_bound = bool(runtime_artifacts) and all(
            (root / str(artifact.get("path", ""))).exists()
            for artifact in runtime_artifacts
            if str(artifact.get("path", ""))
        )
        entries.append(
            {
                "id": str(row.get("id", "unknown")),
                "citation": citation,
                "normalized_citation": _normalize_text(citation),
                "author_year_key": _author_year_key(citation),
                "runtime_artifacts": runtime_artifacts,
                "runtime_bound": runtime_bound,
            }
        )
    return entries


def parse_markdowns(
    directory: Path = DEEP_RESEARCH_DIR,
    *,
    min_score: int = 6,
) -> List[Dict[str, Any]]:
    occurrences: List[Dict[str, Any]] = []
    for filepath in iter_markdown_paths(directory):
        with open(filepath, "r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                match = ENTRY_PATTERN.search(line)
                if not match:
                    continue
                citation = str(match.group("citation_bt") or match.group("citation_plain") or "").strip()
                description = str(match.group("description") or "").strip()
                score = str(match.group("score") or "").strip()
                trailing = str(match.group("trailing") or "").strip()
                if trailing:
                    description = f"{description} {trailing}".strip()
                score_value = parse_score(score)
                if score_value < int(min_score):
                    continue
                occurrences.append(
                    {
                        "file": filepath.name,
                        "line": line_number,
                        "citation": citation,
                        "description": description,
                        "score": score,
                        "score_value": score_value,
                        "normalized_citation": _normalize_text(citation),
                        "author_year_key": _author_year_key(citation),
                        "raw_line": line,
                    }
                )
    return occurrences


def resolve_registry_match(item: Dict[str, Any], registry_entries: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    normalized = str(item.get("normalized_citation", ""))
    author_year = str(item.get("author_year_key", ""))
    for entry in registry_entries:
        if normalized and (
            normalized == entry["normalized_citation"]
            or normalized in entry["normalized_citation"]
            or entry["normalized_citation"] in normalized
        ):
            return entry
        if author_year and author_year == entry["author_year_key"]:
            return entry
    return None


def collapse_occurrences(
    occurrences: List[Dict[str, Any]],
    registry_entries: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    grouped: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"occurrences": [], "files": set(), "descriptions": []})
    for item in occurrences:
        key = str(item.get("normalized_citation", "")) or str(item.get("citation", ""))
        bucket = grouped[key]
        bucket["citation"] = item["citation"]
        bucket["normalized_citation"] = item["normalized_citation"]
        bucket["author_year_key"] = item["author_year_key"]
        bucket["score_value"] = max(int(bucket.get("score_value", 0)), int(item["score_value"]))
        bucket["score"] = f"{bucket['score_value']}/8"
        bucket["occurrences"].append(
            {
                "file": item["file"],
                "line": item["line"],
                "description": item["description"],
                "raw_line": item["raw_line"],
            }
        )
        bucket["files"].add(item["file"])
        if item["description"] not in bucket["descriptions"]:
            bucket["descriptions"].append(item["description"])

    results: List[Dict[str, Any]] = []
    for bucket in grouped.values():
        match = resolve_registry_match(bucket, registry_entries)
        status = "BACKLOG"
        if match is not None:
            status = "RUNTIME_BOUND" if bool(match.get("runtime_bound", False)) else "REGISTRY_ONLY"
        results.append(
            {
                "citation": bucket["citation"],
                "score": bucket["score"],
                "score_value": int(bucket["score_value"]),
                "status": status,
                "occurrence_count": len(bucket["occurrences"]),
                "files": sorted(bucket["files"]),
                "descriptions": bucket["descriptions"],
                "occurrences": bucket["occurrences"],
                "registry_id": None if match is None else match.get("id"),
                "runtime_artifact_count": 0 if match is None else len(match.get("runtime_artifacts", [])),
            }
        )
    return sorted(results, key=lambda item: (item["status"] != "BACKLOG", -int(item["score_value"]), item["citation"].casefold()))


def summarize_priority_targets(items: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    summary: Dict[str, Dict[str, Any]] = {}
    for target, needles in PRIORITY_53_TARGETS.items():
        matched_items: List[Dict[str, Any]] = []
        for item in items:
            haystack = " ".join(
                [
                    str(item.get("citation", "")),
                    " ".join(str(desc) for desc in item.get("descriptions", [])),
                    " ".join(str(occ.get("raw_line", "")) for occ in item.get("occurrences", [])),
                ]
            ).casefold()
            if any(needle.casefold() in haystack for needle in needles):
                matched_items.append(item)
        statuses = {str(item.get("status", "BACKLOG")) for item in matched_items}
        if "RUNTIME_BOUND" in statuses:
            target_status = "RUNTIME_BOUND"
        elif "REGISTRY_ONLY" in statuses:
            target_status = "REGISTRY_ONLY"
        elif matched_items:
            target_status = "BACKLOG"
        else:
            target_status = "NOT_FOUND"
        summary[target] = {
            "status": target_status,
            "matches": [item["citation"] for item in matched_items],
        }
    return summary


def build_audit_payload(*, min_score: int = 6) -> Dict[str, Any]:
    registry_entries = load_registry_entries()
    occurrences = parse_markdowns(min_score=min_score)
    items = collapse_occurrences(occurrences, registry_entries)
    summary = {
        "total_occurrences": len(occurrences),
        "unique_citations": len(items),
        "runtime_bound": sum(1 for item in items if item["status"] == "RUNTIME_BOUND"),
        "registry_only": sum(1 for item in items if item["status"] == "REGISTRY_ONLY"),
        "backlog": sum(1 for item in items if item["status"] == "BACKLOG"),
        "min_score": int(min_score),
    }
    return {
        "summary": summary,
        "priority_53_targets": summarize_priority_targets(items),
        "items": items,
    }


def write_outputs(payload: Dict[str, Any]) -> None:
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    summary = payload["summary"]
    target_summary = payload["priority_53_targets"]
    items = payload["items"]
    OUTPUT_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_MD, "w", encoding="utf-8") as handle:
        handle.write("# Deep Research Gap Analysis\n\n")
        handle.write(f"**Eligible occurrence count (>= {summary['min_score']}/8):** {summary['total_occurrences']}\n")
        handle.write(f"**Unique citations:** {summary['unique_citations']}\n")
        handle.write(f"**Runtime-bound in code database:** {summary['runtime_bound']}\n")
        handle.write(f"**Tracked in registry only:** {summary['registry_only']}\n")
        handle.write(f"**Backlog:** {summary['backlog']}\n\n")

        handle.write("## Priority 5.3 Target Status\n")
        handle.write("| Target | Status | Matches |\n")
        handle.write("| :--- | :--- | :--- |\n")
        for target, row in target_summary.items():
            matches = "; ".join(row.get("matches", [])) or "-"
            handle.write(f"| {target} | {row.get('status', 'UNKNOWN')} | {matches} |\n")
        handle.write("\n")

        handle.write("## Backlog Items\n")
        for item in items:
            if item["status"] != "BACKLOG":
                continue
            handle.write(f"- **{item['citation']}** ({item['score']})\n")
            handle.write(f"  - files: {', '.join(item['files'])}\n")
            for description in item.get("descriptions", [])[:3]:
                handle.write(f"  - {description}\n")

        handle.write("\n## Runtime-bound Items\n")
        for item in items:
            if item["status"] != "RUNTIME_BOUND":
                continue
            handle.write(f"- [x] **{item['citation']}** ({item['score']})\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit benchmark-eligible Deep Research citations against the current code database.")
    parser.add_argument("--min-score", type=int, default=6, help="Minimum eligibility score to include from the markdown corpus (default: 6).")
    args = parser.parse_args()

    payload = build_audit_payload(min_score=args.min_score)
    write_outputs(payload)
    summary = payload["summary"]
    print("Starting Deep Research Tracker...")
    print(
        "Extraction complete. "
        f"Found {summary['total_occurrences']} eligible occurrences across {summary['unique_citations']} unique citations; "
        f"runtime-bound={summary['runtime_bound']}, registry-only={summary['registry_only']}, backlog={summary['backlog']}."
    )
    print(f"Reports saved to {OUTPUT_JSON.relative_to(ROOT)} and {OUTPUT_MD.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
