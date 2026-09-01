"""LLM-DIGEST ECHO CENSUS. **THIS SCRIPT DOES NOT VERIFY ANYTHING.**

2026-08-27 (Wave T3) -- DISARMED. Read this before believing its output.

WHAT THIS SCRIPT ACTUALLY MEASURES
----------------------------------
For each registry entry in ``data/lit/*.json`` it collects every numeric leaf, splits
the entry's citation string into surnames, then searches the **LLM-generated markdown
under ``data/Gemini_Deep_Research/``** for lines containing one of those surnames. If a
number appears as a **substring** anywhere in a +/-30-line window around such a line,
that number is recorded as "echoed".

That is the entire test. In full:

* the corpus searched is **machine-generated text**, not literature;
* the match is textual proximity, not attribution -- nothing checks that the digest
  sentence is *about* the quantity, the compound, the units, or the conditions;
* the match is a bare substring, so ``4.5`` matches ``14.52``, ``0.045`` and ``2004.5``;
* a surname match is case-insensitive and unanchored, so ``Bi`` matches ``binding``,
  ``combination`` and every other word containing "bi";
* a +/-30-line window in these digests routinely spans several unrelated references.

WHAT IT DOES **NOT** DO
-----------------------
It never opens a paper. It never contacts CrossRef, PubMed, Unpaywall or a publisher.
It cannot distinguish a number the digest copied correctly from one the digest invented.
It has no notion of a primary source at all.

WHY THE VOCABULARY CHANGED
--------------------------
This script previously printed its top-line result as **"Fully Verified (All values
matched): 153 (57.5%)"**, under the report title *"Numeric Value Traceability and
**Verification** Report"*. Those 153 rows were a census of values whose only known
upstream is an LLM digest -- i.e. the **laundering census** -- published under a heading
asserting the opposite. It directly contradicted this repository's own rule
(``data/Gemini_Deep_Research/README.md``):

    "The deep-research report says so" is not provenance.

The old and new vocabulary, one-for-one:

    OLD                                      NEW
    "Fully Verified (All values matched)" -> "DIGEST-ECHO (NOT VERIFICATION)"
    "Partially Verified (Some values ...)"-> "PARTIAL DIGEST-ECHO (NOT VERIFICATION)"
    "Unverified (No values matched)"      -> "NO DIGEST ECHO (origin unaccounted for)"
    "Fully Verified Entries"              -> "Entries whose every number echoes an LLM digest"
    "...verified in the Deep Research
       markdown files"                    -> "...found only inside LLM-generated text;
                                              no paper was opened"
    report title: "... Verification Report"->"LLM-Digest Echo Census -- NOT a verification
                                              result"

A HIGH ECHO COUNT IS A BAD SIGN, NOT A GOOD ONE
-----------------------------------------------
Read the output backwards from how the old wording invited. An entry in the DIGEST-ECHO
class is one for which this repository can point at nothing but its own machine-generated
research dump. That is the condition the audit is trying to eliminate. An entry in the
NO DIGEST ECHO class is *not* thereby better -- it merely has no identified upstream at
all, which may be worse.

Nothing in this output may be cited as evidence that a value is verified, sourced,
anchored, or safe to ship. Verification means a human opened the primary source; see
``results/validation/citation_verification_ledger.md`` for the (identity-only, CrossRef)
pass that does exist, and note that even that one checks the DOI, not the number.
"""

import json
import re
import sys
from pathlib import Path

# 2026-08-27 (Wave T3): was `Path(".").resolve()`, i.e. the process CWD, so the script
# silently produced an empty census unless invoked from the repo root. Anchored to the
# file instead.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

DATA_LIT_DIR = data_paths.LIT_DIR
GDR_DIR = data_paths.RESEARCH_CORPUS_DIR
OUTPUT_REPORT_PATH = data_paths.VALIDATION_DIR / "key_value_trace_report.md"

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
    #
    # 2026-08-27 (Wave T3): the three counters below are UNCHANGED in what they
    # compute. Only their names and their printed labels changed, because the old
    # names asserted a fact the computation cannot establish. `fully_verified` was
    # never a count of verified entries; it is the count of entries all of whose
    # numbers can be found inside LLM-generated text near a surname.
    total_entries = len(results)
    full_digest_echo = sum(1 for r in results.values() if len(r["unmatched_vals"]) == 0)
    partial_digest_echo = sum(1 for r in results.values() if 0 < len(r["unmatched_vals"]) < len(r["all_vals"]))
    no_digest_echo = sum(1 for r in results.values() if len(r["unmatched_vals"]) == len(r["all_vals"]))

    report_lines = [
        "# LLM-Digest Echo Census — **NOT a verification result**",
        "",
        "> ## ⚠ READ THIS BEFORE READING ANY NUMBER BELOW",
        "> ",
        "> **This report verifies nothing. No paper was opened to produce it.**",
        "> ",
        "> It measures ONE thing: whether a registry entry's numbers appear as text within",
        "> ±30 lines of one of its citation surnames inside the **machine-generated markdown**",
        "> under `data/Gemini_Deep_Research/`. Matching is bare substring matching (`4.5`",
        "> matches `14.52`) and surname matching is unanchored (`Bi` matches `binding`).",
        "> No primary source, DOI resolver, publisher or index is contacted at any point.",
        "> ",
        "> Per this repository's own rule (`data/Gemini_Deep_Research/README.md`):",
        "> ",
        "> > \"The deep-research report says so\" **is not provenance.**",
        "> ",
        "> **Therefore the DIGEST-ECHO class below is the LAUNDERING CENSUS, not a clean bill",
        "> of health.** A high echo count is a bad sign: it means that for those entries this",
        "> repository can point at nothing but its own LLM research dump. Entries in the",
        "> NO DIGEST ECHO class are not thereby better — they simply have no identified",
        "> upstream at all, which may be worse.",
        "> ",
        "> Nothing in this file may be cited as evidence that a value is verified, sourced,",
        "> anchored or safe to ship. The only verification pass that exists in this repo is",
        "> `results/validation/citation_verification_ledger.md`, and even that checks the DOI's",
        "> identity, not the number.",
        "> ",
        "> *History: until 2026-08-27 this file was titled \"Numeric Value Traceability and",
        "> **Verification** Report\" and printed its top line as \"**Fully Verified** (All values",
        "> matched): 153 (57.5%)\". Same computation, opposite claim. Retitled and re-worded by",
        "> Wave T3 of the audit remediation; see `scripts/trace_key_values.py` for the full",
        "> old→new vocabulary map and `tasks/audit_remediation.md` § Wave T3.*",
        "",
        "---",
        "",
        "## Summary Metrics",
        f"- **Total Entries Evaluated**: {total_entries}",
        f"- **DIGEST-ECHO (NOT VERIFICATION) — every number echoes the LLM digest corpus**: "
        f"{full_digest_echo} ({full_digest_echo/total_entries*100:.1f}%)",
        f"- **PARTIAL DIGEST-ECHO (NOT VERIFICATION) — some numbers echo the corpus**: "
        f"{partial_digest_echo} ({partial_digest_echo/total_entries*100:.1f}%)",
        f"- **NO DIGEST ECHO — origin unaccounted for even by the digests**: "
        f"{no_digest_echo} ({no_digest_echo/total_entries*100:.1f}%)",
        "",
        "---",
        "",
        "## Entries with NO / PARTIAL digest echo",
        "The numbers listed for each entry were **not** found near their citation surname in the",
        "LLM digest corpus. This says nothing about whether they are right: it says only that the",
        "digests do not account for them, so their origin is unidentified by this script.",
        ""
    ]

    report_lines.append("| Entry ID | Citation | Sources | Values NOT echoed anywhere in the digest corpus |")
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
        "## DIGEST-ECHO entries — every number found only inside LLM-generated text",
        "**NOT VERIFICATION.** For each entry below, every numeric value was located inside",
        "`data/Gemini_Deep_Research/**` within ±30 lines of a citation surname. No paper was",
        "opened. This is the list of entries whose only demonstrated upstream is a machine-",
        "generated digest — i.e. the remediation worklist, not the safe list.",
        ""
    ])

    report_lines.append("| Entry ID | Citation | Sources | Values echoed in the digest corpus (unverified) |")
    report_lines.append("|---|---|---|---|")

    for entry_id, r in sorted(results.items()):
        if not r["unmatched_vals"]:
            sources_str = ", ".join(r["sources"])
            vals_str = ", ".join(map(str, sorted(r["all_vals"])))
            report_lines.append(f"| `{entry_id}` | {r['citation']} | {sources_str} | {vals_str} |")

    report_lines.append("")

    # Write report
    OUTPUT_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    print(f"Digest-echo census written to {OUTPUT_REPORT_PATH}")
    print("NOTE: this script verifies nothing. It measures echo into LLM-generated text.")
    print(f"DIGEST-ECHO (not verified): {full_digest_echo}/{total_entries}")

if __name__ == "__main__":
    trace()
