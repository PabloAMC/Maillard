import json
import re
from pathlib import Path

ROOT = Path(".").resolve()
DATA_LIT_DIR = ROOT / "data" / "lit"
GDR_DIR = ROOT / "data" / "Gemini_Deep_Research"
OUTPUT_REPORT_PATH = ROOT / "results" / "validation" / "key_value_trace_report.md"

def load_json(path):
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def extract_numeric_values(data, found_vals=None):
    if found_vals is None:
        found_vals = set()
    if isinstance(data, dict):
        for k, v in data.items():
            # Skip keys that are just metadata like years, line numbers, ids
            if k in ["year", "id", "paper_id", "registry_id", "line", "replicates", "slr_family_source", "implementation_wave", "order_in_wave"]:
                continue
            extract_numeric_values(v, found_vals)
    elif isinstance(data, list):
        for item in data:
            extract_numeric_values(item, found_vals)
    elif isinstance(data, (int, float)):
        # Ignore common boilerplate numbers like 0, 1, 2, 3
        if data not in [0, 1, 2, 3, 4, 6]:
            found_vals.add(data)
    return found_vals

def clean_citation(cit):
    # E.g. "Trikusuma et al. (2019)" -> "Trikusuma"
    # E.g. "Mottram & Wedzicha (1990)" -> ["Mottram", "Wedzicha"]
    # E.g. "Brands et al. (2002)" -> "Brands"
    if not cit:
        return []
    cit_clean = re.sub(r'\(.*?\)', '', cit) # remove years
    cit_clean = re.sub(r'et al\.?', '', cit_clean) # remove et al.
    parts = re.split(r'[&,;]', cit_clean)
    names = []
    for p in parts:
        name = p.strip()
        if len(name) > 2:
            names.append(name)
    return names

def search_value_in_text(text_window, val):
    # Check for float and int versions
    val_str = str(val)
    if val_str in text_window:
        return True
    
    # If float ends with .0, check for the int version
    if isinstance(val, float) and val.is_integer():
        if str(int(val)) in text_window:
            return True
            
    # Check with comma as decimal separator (European style sometimes in raw files)
    val_comma = val_str.replace('.', ',')
    if val_comma in text_window:
        return True
        
    return False

def trace():
    # 1. Load markdown files
    md_files = list(GDR_DIR.glob("**/*.md"))
    md_contents = {}
    for fp in md_files:
        try:
            with open(fp, "r", encoding="utf-8") as f:
                md_contents[fp.relative_to(ROOT)] = f.read().splitlines()
        except Exception as e:
            print(f"Error reading {fp}: {e}")

    # 2. Load database files
    registries = {
        "benchmark_intake_registry.json": load_json(DATA_LIT_DIR / "benchmark_intake_registry.json").get("eligible_references", []),
        "flavor_reference_payloads.json": [item for sec in load_json(DATA_LIT_DIR / "flavor_reference_payloads.json").values() if isinstance(sec, list) for item in sec],
        "computational_priors.json": [item for sec in load_json(DATA_LIT_DIR / "computational_priors.json").values() if isinstance(sec, list) for item in sec],
        "retention_reference_payloads.json": [item for sec in load_json(DATA_LIT_DIR / "retention_reference_payloads.json").values() if isinstance(sec, list) for item in sec],
        "safety_reference_payloads.json": load_json(DATA_LIT_DIR / "safety_reference_payloads.json").get("entries", []),
    }

    # Gather all unique entries to verify
    entries_to_trace = {}
    for db_name, items in registries.items():
        for item in items:
            if not isinstance(item, dict):
                continue
            entry_id = item.get("id") or item.get("paper_id")
            if not entry_id:
                continue
            
            citation = item.get("citation") or item.get("source") or item.get("source_citation")
            if not citation:
                # Try to construct citation from ID
                citation = entry_id.split("_")[0].capitalize()

            numeric_vals = extract_numeric_values(item)
            if not numeric_vals:
                continue

            if entry_id not in entries_to_trace:
                entries_to_trace[entry_id] = {
                    "citation": citation,
                    "values": set(),
                    "sources": set()
                }
            entries_to_trace[entry_id]["values"].update(numeric_vals)
            entries_to_trace[entry_id]["sources"].add(db_name)

    print(f"Loaded {len(entries_to_trace)} entries with numeric values to trace.")

    # 3. Perform the tracing
    results = {}
    for entry_id, info in entries_to_trace.items():
        cit_names = clean_citation(info["citation"])
        vals = info["values"]
        
        matches = []
        unmatched_vals = set(vals)

        for filepath, lines in md_contents.items():
            # Find citation occurrences
            for idx, line in enumerate(lines):
                # Check if any of the citation names appear in this line
                if any(name.lower() in line.lower() for name in cit_names):
                    # Check window of +/- 30 lines
                    start_win = max(0, idx - 30)
                    end_win = min(len(lines), idx + 30)
                    window_text = "\n".join(lines[start_win:end_win])

                    # Check each unmatched value in this window
                    found_in_win = []
                    for val in list(unmatched_vals):
                        if search_value_in_text(window_text, val):
                            found_in_win.append(val)
                            unmatched_vals.remove(val)
                    
                    if found_in_win:
                        matches.append({
                            "file": str(filepath),
                            "line_num": idx + 1,
                            "matched_vals": found_in_win
                        })

        results[entry_id] = {
            "citation": info["citation"],
            "all_vals": list(vals),
            "unmatched_vals": list(unmatched_vals),
            "matches": matches,
            "sources": list(info["sources"])
        }

    # 4. Generate the MD Report
    total_entries = len(results)
    fully_verified = sum(1 for r in results.values() if len(r["unmatched_vals"]) == 0)
    partially_verified = sum(1 for r in results.values() if 0 < len(r["unmatched_vals"]) < len(r["all_vals"]))
    unverified = sum(1 for r in results.values() if len(r["unmatched_vals"]) == len(r["all_vals"]))

    report_lines = [
        "# Numeric Value Traceability and Verification Report",
        f"This report traces the numeric values (`key_values`, `Ea`, `yields`, etc.) in the Maillard database registries back to the raw and structured Deep Research markdown corpus.",
        "",
        "## Summary Metrics",
        f"- **Total Entries Evaluated**: {total_entries}",
        f"- **Fully Verified (All values matched)**: {fully_verified} ({fully_verified/total_entries*100:.1f}%)",
        f"- **Partially Verified (Some values matched)**: {partially_verified} ({partially_verified/total_entries*100:.1f}%)",
        f"- **Unverified (No values matched)**: {unverified} ({unverified/total_entries*100:.1f}%)",
        "",
        "---",
        "",
        "## Unverified & Partially Verified Entries (Action Required)",
        "The following registry entries have numerical values that could not be verified in the surrounding context of their citation in the markdown files.",
        ""
    ]

    report_lines.append("| Entry ID | Citation | Sources | Mismatched / Unmatched Values |")
    report_lines.append("|---|---|---|---|")
    
    for entry_id, r in sorted(results.items()):
        if r["unmatched_vals"]:
            sources_str = ", ".join(r["sources"])
            vals_str = ", ".join(map(str, sorted(r["unmatched_vals"])))
            report_lines.append(f"| `{entry_id}` | {r['citation']} | {sources_str} | {vals_str} |")

    report_lines.extend([
        "",
        "---",
        "",
        "## Fully Verified Entries",
        "The following registry entries are 100% matched and verified in the Deep Research markdown files.",
        ""
    ])

    report_lines.append("| Entry ID | Citation | Sources | Verified Values |")
    report_lines.append("|---|---|---|---|")
    
    for entry_id, r in sorted(results.items()):
        if not r["unmatched_vals"]:
            sources_str = ", ".join(r["sources"])
            vals_str = ", ".join(map(str, sorted(r["all_vals"])))
            report_lines.append(f"| `{entry_id}` | {r['citation']} | {sources_str} | {vals_str} |")

    # Write report
    OUTPUT_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    print(f"Trace report successfully written to {OUTPUT_REPORT_PATH}")
    print(f"Verified: {fully_verified}/{total_entries}")

if __name__ == "__main__":
    trace()
