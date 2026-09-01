#!/usr/bin/env python3
"""
ingest_deep_research_markdown.py

Parses Deep Research Markdown reports (e.g. data/Gemini_Deep_Research/*_*.md),
extracts the 'Consolidated entries' sections which have been scored via the 8-point SLR,
and appends them to a dedicated backlog file instead of the operational
benchmark intake registry.

Includes --human-review interactive gate specifically to support S12.2.
"""

import os
import glob
import json
import re
import argparse
import sys
from datetime import datetime

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from src import data_paths  # noqa: E402

DEEP_RESEARCH_DIR = str(data_paths.RESEARCH_CORPUS_DIR)
REGISTRY_PATH = str(data_paths.BENCHMARK_INTAKE_REGISTRY)
# Orphan output: this file does not exist in the repo and nothing reads it.
BACKLOG_PATH = str(data_paths.DEEP_RESEARCH_CANDIDATE_REGISTRY)


def _load_json(path: str) -> dict:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _default_backlog_payload() -> dict:
    return {
        "generated_at": "",
        "source": f"{data_paths.rel(data_paths.RESEARCH_CORPUS_DIR)}/*.md",
        "notes": "Markdown-derived Deep Research candidate backlog. This file is not consumed by the runtime or benchmark registry.",
        "candidate_entries": [],
    }

def sanitize_id(citation: str) -> str:
    """Creates a basic snake case id from a citation string."""
    cleaned = re.sub(r'[^a-zA-Z0-9_\s]', '', citation)
    cleaned = cleaned.replace(" ", "_").lower()
    return cleaned[:30].strip("_")

def parse_markdown_files():
    parsed_entries = []
    
    md_files = glob.glob(os.path.join(DEEP_RESEARCH_DIR, "*.md"))
    for md_file in md_files:
        if "cross_family" in md_file or "maillard_chemical_families" in md_file:
            continue
            
        with open(md_file, "r", encoding="utf-8") as f:
            content = f.read()
            
        # Try to find the SLR family number from filename e.g. 06_alternative_proteins.md
        family_match = re.search(r'(\d+)_([a-zA-Z_]+)\.md', os.path.basename(md_file))
        if not family_match:
            continue
            
        family_num = family_match.group(1)
        family_name = family_match.group(2)
        
        # Look for the consolidated section
        consolidated_match = re.search(r'## Consolidated entries for `benchmark_schema\.json`.*?([\s\S]+)', content)
        if not consolidated_match:
            continue
            
        section_content = consolidated_match.group(1)
        
        # Split by the bold headers
        categories = {
            "REFERENCE_ANCHOR": "reference_anchor",
            "PRIMARY": "benchmark_eligible",
            "CALIBRATION": "conditional_calibration"
        }
        
        current_cat = None
        for line in section_content.split('\n'):
            line = line.strip()
            
            # Check for category shifts
            cat_found = False
            for k, v in categories.items():
                if k in line and line.startswith("**"):
                    current_cat = v
                    cat_found = True
                    break
                    
            if cat_found:
                continue
                
            # Match list items like `1. `Citation` — details. Score/8.`
            match = re.match(r'^\d+\.\s+`([^`]+)`\s+[—\-]\s+(.*)', line)
            if match and current_cat:
                citation = match.group(1)
                desc = match.group(2)
                
                # Check for score
                score_match = re.search(r'(\d+)/8', desc)
                score = score_match.group(1) if score_match else "??"
                
                entry = {
                    "id": sanitize_id(citation),
                    "citation": citation,
                    "doi": "",
                    "kind": current_cat,
                    "chemistry_family": family_name,
                    "slr_family_source": family_num,
                    "payload_role": "benchmark_intake" if current_cat in ["reference_anchor", "benchmark_eligible"] else "calibration",
                    "observable_panel_tags": ["extracted_from_markdown"],
                    "process_state_scope": [],
                    "supporting_families": [],
                    "target_modules": ["matrix_correction", "headspace"] if family_num == "06" else [],
                    "matrix_family": "unknown_see_details",
                    "eligibility": current_cat,
                    "status": "pending_json_payload",
                    "what_it_supports": [desc],
                    "what_it_does_not_support": [],
                    "key_values": {"slr_score": f"{score}/8"},
                    "repo_next_action": "Generate canonical benchmark_payload for exact ppb values",
                    "runtime_artifacts": []
                }
                parsed_entries.append(entry)
                
    return parsed_entries


def main():
    parser = argparse.ArgumentParser(description="Ingest Deep Research Markdown")
    parser.add_argument("--dry-run", action="store_true", help="Do not save, just print")
    parser.add_argument("--human-review", action="store_true", help="Prompt for each entry")
    args = parser.parse_args()
    
    entries = parse_markdown_files()
    print(f"Parsed {len(entries)} references from Gemini_Deep_Research markdown.")
    
    if args.dry_run:
        for e in entries:
            print(json.dumps(e, indent=2))
        return

    registry = _load_json(REGISTRY_PATH)
    backlog = _load_json(BACKLOG_PATH) or _default_backlog_payload()

    existing_ids = {
        str(ref.get("id", "")).strip(): True
        for ref in registry.get("eligible_references", [])
        if str(ref.get("id", "")).strip()
    }
    existing_ids.update(
        {
            str(ref.get("id", "")).strip(): True
            for ref in backlog.get("candidate_entries", [])
            if str(ref.get("id", "")).strip()
        }
    )
    
    added_count = 0
    for e in entries:
        if e["id"] in existing_ids:
            print(f"Skipping {e['id']} - already exists.")
            continue
            
        if args.human_review:
            print("\n----------")
            print(json.dumps(e, indent=2))
            ans = input("Approve ingestion? (y/n): ")
            if ans.lower() != 'y':
                print("Skipped.")
                continue
                
        backlog.setdefault("candidate_entries", []).append(e)
        existing_ids[e["id"]] = True
        added_count += 1

    backlog["generated_at"] = datetime.now().strftime("%Y-%m-%d")
    
    if added_count > 0:
        with open(BACKLOG_PATH, "w", encoding="utf-8") as f:
            json.dump(backlog, f, indent=2)
        print(f"\nAdded {added_count} new entries to {BACKLOG_PATH}.")
    else:
        print("\nNo new entries added.")


if __name__ == "__main__":
    main()
