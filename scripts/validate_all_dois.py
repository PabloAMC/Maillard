import json
import urllib.request
import urllib.error
import time
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402

REGISTRY_PATH = data_paths.BENCHMARK_INTAKE_REGISTRY

def clean_doi(doi):
    if not doi:
        return ""
    # Strip url prefix if present
    doi = doi.strip()
    for prefix in ["https://doi.org/", "http://doi.org/", "doi.org/"]:
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix):]
    return doi.strip()

def parse_citation(citation):
    """
    Parses a citation string like 'Ma et al. (2024)' or 'Morel et al. (2002)'
    Returns (author_name, year) or (None, None)
    """
    if not citation:
        return None, None
    # Match Author (Year) or Author et al. (Year)
    match = re.search(r"^([A-Za-z\-'\u00C0-\u017F]+)(?:\s+et\s+al\.)?\s*\((\d{4})\)", citation)
    if match:
        return match.group(1), int(match.group(2))
    return None, None

def query_crossref(doi):
    """
    Queries Crossref for DOI metadata.
    """
    url = f"https://api.crossref.org/works/{doi}"
    headers = {
        "User-Agent": "MaillardDOIVerifier/1.0 (mailto:pablo@maillard.org; pairing as Antigravity)"
    }
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode("utf-8"))
            return data.get("message", {})
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {"error": "Not Found (404)"}
        return {"error": f"HTTP Error {e.code}"}
    except Exception as e:
        return {"error": f"Connection Error: {str(e)}"}

def validate():
    if not REGISTRY_PATH.exists():
        print(f"Error: Registry file not found at {REGISTRY_PATH}")
        return

    with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
        registry = json.load(f)

    references = registry.get("eligible_references", [])
    print(f"Loaded {len(references)} references from registry.")

    # Track unique DOIs to query once
    doi_to_refs = {}
    for ref in references:
        doi = clean_doi(ref.get("doi", ""))
        # Ignore empty, thesis, or patent IDs
        if not doi or not ("/" in doi):
            continue
        doi_to_refs.setdefault(doi, []).append(ref)

    print(f"Found {len(doi_to_refs)} unique valid DOIs to check.\n")
    print(f"{'DOI':<30} | {'Citation':<30} | {'Status':<10} | {'Details'}")
    print("-" * 100)

    success_count = 0
    warning_count = 0
    failure_count = 0

    for doi, refs in sorted(doi_to_refs.items()):
        # Query API (polite rate limiting)
        time.sleep(0.2)
        
        metadata = query_crossref(doi)
        
        if "error" in metadata:
            status = "FAILED"
            details = metadata["error"]
            failure_count += 1
            print(f"{doi:<30} | {refs[0]['citation']:<30} | {status:<10} | {details}")
            continue

        # Extract title, authors, and year from Crossref metadata
        title = metadata.get("title", [""])[0]
        # Clean title HTML tags
        title = re.sub(r"<[^>]*>", "", title)
        
        authors = metadata.get("author", [])
        author_surnames = [a.get("family", "") for a in authors if "family" in a]
        
        # Get year
        pub_parts = metadata.get("published-print", {}).get("date-parts", [[None]])[0]
        if not pub_parts or pub_parts[0] is None:
            pub_parts = metadata.get("published-online", {}).get("date-parts", [[None]])[0]
        year = pub_parts[0] if pub_parts else None

        # Compare with each reference referencing this DOI
        for ref in refs:
            citation = ref.get("citation", "")
            exp_author, exp_year = parse_citation(citation)
            
            mismatches = []
            
            # Author matching (case-insensitive check of last names)
            if exp_author:
                author_matched = False
                for surname in author_surnames:
                    if exp_author.lower() in surname.lower() or surname.lower() in exp_author.lower():
                        author_matched = True
                        break
                if not author_matched:
                    mismatches.append(f"Author mismatch (Expected: {exp_author}, Got: {', '.join(author_surnames[:3])})")
            
            # Year matching
            if exp_year and year:
                if abs(exp_year - year) > 1: # Allow +/- 1 year for publication vs online posting
                    mismatches.append(f"Year mismatch (Expected: {exp_year}, Got: {year})")

            if mismatches:
                status = "WARNING"
                details = "; ".join(mismatches)
                warning_count += 1
            else:
                status = "OK"
                details = f"Verified: '{title[:45]}...' by {author_surnames[0] if author_surnames else 'Unknown'} ({year})"
                success_count += 1
                
            print(f"{doi:<30} | {citation:<30} | {status:<10} | {details}")

    print("\n" + "=" * 100)
    print(f"Validation finished: {success_count} OK, {warning_count} Warnings, {failure_count} Failures.")
    print("=" * 100)

if __name__ == "__main__":
    validate()
