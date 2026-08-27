import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_LIT_DIR = ROOT / "data" / "lit"
GDR_DIR = ROOT / "data" / "Gemini_Deep_Research"

INTAKE_REGISTRY_PATH = DATA_LIT_DIR / "benchmark_intake_registry.json"
SLR_MATRIX_PATH = DATA_LIT_DIR / "slr_incorporation_matrix.json"
BACKLOG_PATH = DATA_LIT_DIR / "deep_research_backlog.json"


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def extract_references_from_file(file_path):
    refs = []
    print(f"Parsing {file_path.name}...")
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Find bibliography or works cited sections
    bib_started = False
    for line in lines:
        if "structured bibliography" in line.lower() or "structured kinetic bibliography" in line.lower():
            bib_started = True
            continue
        if bib_started and ("works cited" in line.lower() or "references" in line.lower() and not "structured" in line.lower()):
            # We can still parse, but this marks the end of the structured bib
            pass

        # Match lines with DOI
        doi_match = re.search(r'\b10\.\d{4,9}/[^\s—,\]\)\#]+', line)
        if doi_match:
            doi = doi_match.group(0).strip().rstrip(".").rstrip(",")
            
            # Get citation
            cit_match = re.search(r'^\d+[\.\\\s]+([^\(]+(?:\(\d+\))?)', line)
            cit = cit_match.group(1).strip() if cit_match else line[:50].strip()

            refs.append({
                "line": line.strip(),
                "citation": cit,
                "doi": doi
            })
    return refs


def analyze():
    # Load registries
    intake = load_json(INTAKE_REGISTRY_PATH)
    matrix = load_json(SLR_MATRIX_PATH)
    backlog = load_json(BACKLOG_PATH)

    registered_dois = {ref.get("doi").lower().strip() for ref in intake["eligible_references"] if ref.get("doi")}
    matrix_dois = {ent.get("doi").lower().strip() for ent in matrix["entries"] if ent.get("doi")}
    backlog_dois = {item.get("doi", "").lower().strip() for item in backlog["items"] if item.get("doi")}
    backlog_citations = {item.get("citation", "").lower().strip() for item in backlog["items"]}

    files = ["Gemini_deep_research_A.md", "Gemini_deep_research_B.md", "Gemini_deep_research_C.md", "Gemini_deep_research_D.md", "Gemini_deep_research_E.md"]
    
    all_refs = {}
    for fn in files:
        fp = GDR_DIR / fn
        if fp.exists():
            all_refs[fn] = extract_references_from_file(fp)

    print("\n=== ANALYSIS OF DEEP RESEARCH REFERENCES ===")
    for fn, refs in all_refs.items():
        print(f"\nFile: {fn} (found {len(refs)} references)")
        for idx, ref in enumerate(refs, 1):
            doi = ref["doi"]
            cit = ref["citation"]
            doi_lower = doi.lower().strip() if doi else ""
            
            in_registry = doi_lower in registered_dois
            in_matrix = doi_lower in matrix_dois
            in_backlog_doi = doi_lower in backlog_dois if doi_lower else False
            in_backlog_cit = cit.lower().strip() in backlog_citations
            
            status = []
            if in_registry: status.append("REGISTRY")
            if in_matrix: status.append("MATRIX")
            if in_backlog_doi or in_backlog_cit: status.append("BACKLOG")
            
            status_str = "+".join(status) if status else "NEW"
            print(f"  {idx}. Citation: {cit[:60]}")
            print(f"     DOI: {doi}")
            print(f"     Status: {status_str}")


if __name__ == "__main__":
    analyze()
