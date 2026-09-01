import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

DATA_LIT_DIR = data_paths.LIT_DIR

json_files = [
    "benchmark_intake_registry.json",
    "slr_incorporation_matrix.json",
    "flavor_reference_payloads.json",
    "safety_reference_payloads.json",
    "retention_reference_payloads.json",
    "process_state_calibrations.json"
]

def find_references_recursive(data):
    """
    Recursively finds all dictionaries representing reference entries.
    An entry is identified by having an 'id' or 'paper_id' and 'doi' or 'citation'.
    """
    entries = []
    if isinstance(data, dict):
        # Check if this dict itself is a reference entry
        is_ref = ("id" in data or "paper_id" in data) and ("doi" in data or "citation" in data or "source_citation" in data)
        if is_ref:
            entries.append(data)
        # Recurse into all values
        for val in data.values():
            entries.extend(find_references_recursive(val))
    elif isinstance(data, list):
        for val in data:
            entries.extend(find_references_recursive(val))
    return entries

def scan_file(filename):
    filepath = DATA_LIT_DIR / filename
    if not filepath.exists():
        print(f"File {filename} does not exist.")
        return []
        
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    items = find_references_recursive(data)
                
    gaps = []
    for idx, item in enumerate(items):
        ref_id = item.get("id") or item.get("paper_id") or f"index_{idx}"
        citation = item.get("citation") or item.get("source_citation") or "Unknown Citation"
        
        # Check if 'doi' key is present
        if "doi" not in item:
            gaps.append((ref_id, citation, "MISSING_KEY"))
        else:
            doi = item["doi"]
            if doi is None:
                gaps.append((ref_id, citation, "NULL_VALUE"))
            elif not isinstance(doi, str):
                gaps.append((ref_id, citation, f"NON_STR_VALUE ({type(doi)})"))
            elif doi.strip() == "":
                gaps.append((ref_id, citation, "EMPTY_STRING"))
            elif "xxxx" in doi or (
                "10." not in doi
                and not doi.startswith("US")
                and not doi.startswith("PAS_")  # Philippine Agricultural Scientist handle
                and not doi.startswith("LiuThesis")
                and not doi.startswith("ClemsonThesis")
                and not doi.startswith("synthesis_calibration")
            ):
                # Identify placeholders, but allow patents starting with US
                gaps.append((ref_id, citation, f"PLACEHOLDER/OTHER ({doi})"))
                
    return gaps

def main():
    print(f"{'Registry File':<35} | {'Reference ID':<45} | {'Citation':<30} | {'Status/Gap'}")
    print("-" * 140)
    
    total_gaps = 0
    for filename in json_files:
        gaps = scan_file(filename)
        if gaps:
            for ref_id, citation, gap_type in gaps:
                print(f"{filename:<35} | {ref_id:<45} | {citation:<30} | {gap_type}")
                total_gaps += 1
                
    print("-" * 140)
    print(f"Scan complete. Found {total_gaps} items with missing, empty, or non-standard DOIs.")

if __name__ == "__main__":
    main()
