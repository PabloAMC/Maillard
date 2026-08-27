import json
import re
import urllib.request
import urllib.parse
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
    data = json.load(f)

# Matches standard format: LeadAuthor et al. (year) or LeadAuthor & CoAuthor (year) or LeadAuthor (year)
# E.g. "Sun et al. (2026)" or "Hofmann et al. (1996)" or "Martins & van Boekel (2003)"
pattern = re.compile(r"^[A-Z\u00C0-\u024F][a-zA-Z\s\u00C0-\u024F\-\']+(?:(\set\sal\.|& [A-Z\u00C0-\u024F][a-zA-Z\s\u00C0-\u024F\-\']+))?\s\(\d{4}\)$")

def fetch_openalex_metadata(identifier):
    # identifier can be doi (e.g. "10.1021/jf960062o"), pmid (e.g. "pmid:36878579"), pmcid (e.g. "pmcid:PMC9905368")
    url = f"https://api.openalex.org/works/{identifier}"
    headers = {
        "User-Agent": "MaillardCitationStandardizer/1.0 (mailto:pablo@maillard.org)"
    }
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            if response.status == 200:
                return json.loads(response.read().decode("utf-8"))
    except Exception as e:
        print(f"  Error fetching {identifier}: {e}")
    return None

def format_authors_and_year(work_data):
    if not work_data:
        return None
    
    # Extract year
    year = work_data.get("publication_year")
    if not year:
        return None
    
    # Extract authors
    authorships = work_data.get("authorships", [])
    if not authorships:
        return None
    
    author_names = []
    for auth in authorships:
        author = auth.get("author", {})
        display_name = author.get("display_name", "")
        if display_name:
            # Try to get family name or last word
            parts = display_name.strip().split()
            if parts:
                family_name = parts[-1]
                # Clean up punctuation
                family_name = re.sub(r"[^\w\s\-']", "", family_name)
                author_names.append(family_name)
    
    if not author_names:
        return None
    
    num_authors = len(author_names)
    if num_authors == 1:
        citation = f"{author_names[0]} ({year})"
    elif num_authors == 2:
        citation = f"{author_names[0]} & {author_names[1]} ({year})"
    else:
        citation = f"{author_names[0]} et al. ({year})"
    
    return citation

mismatches = []
for ref in data.get("eligible_references", []):
    citation = ref.get("citation", "").strip()
    ref_id = ref.get("id", "")
    doi = ref.get("doi", "").strip()
    if not pattern.match(citation):
        mismatches.append(ref)

print(f"Found {len(mismatches)} references to check/standardize.")

fixed_count = 0
for ref in mismatches:
    ref_id = ref["id"]
    current_citation = ref.get("citation", "")
    doi = ref.get("doi", "").strip()
    new_citation = None
    
    print(f"\nProcessing {ref_id} (current: {current_citation!r})")
    
    # 1. Try to query using DOI
    if doi:
        # Normalize DOI to openalex format
        doi_clean = doi
        if not doi_clean.startswith("http"):
            doi_clean = "https://doi.org/" + doi_clean
        print(f"  Querying OpenAlex by DOI: {doi_clean}")
        work_data = fetch_openalex_metadata(doi_clean)
        new_citation = format_authors_and_year(work_data)
        if new_citation:
            print(f"  Resolved from DOI -> {new_citation!r}")
            time.sleep(0.2)
            
    # 2. Try to query by PMC / PMID in ID if DOI wasn't present or failed
    if not new_citation:
        pmc_match = re.search(r"pmc(\d+)", ref_id, re.IGNORECASE)
        pmid_match = re.search(r"pmid(\d+)", ref_id, re.IGNORECASE)
        if pmc_match:
            pmcid = f"pmcid:PMC{pmc_match.group(1)}"
            print(f"  Querying OpenAlex by PMCID: {pmcid}")
            work_data = fetch_openalex_metadata(pmcid)
            new_citation = format_authors_and_year(work_data)
            if new_citation:
                print(f"  Resolved from PMCID -> {new_citation!r}")
                time.sleep(0.2)
        elif pmid_match:
            pmid = f"pmid:{pmid_match.group(1)}"
            print(f"  Querying OpenAlex by PMID: {pmid}")
            work_data = fetch_openalex_metadata(pmid)
            new_citation = format_authors_and_year(work_data)
            if new_citation:
                print(f"  Resolved from PMID -> {new_citation!r}")
                time.sleep(0.2)

    # 3. Heuristic clean up based on current citation or ID
    if not new_citation:
        # Examples: "Wang, Z. et al. (2012)" -> "Wang et al. (2012)"
        # "Hofmann, Schieberle & Grosch (1996), JAFC 44:2873" -> "Hofmann et al. (1996)"
        # "Voelker, Taylor & Mauer (2021)" -> "Voelker et al. (2021)"
        # "Resconi et al. (2023 / PMC10096055)" -> "Resconi et al. (2023)"
        
        # Match standard et al with year
        et_al_year_match = re.search(r"^([A-Z\u00C0-\u024F][a-zA-Z\s\u00C0-\u024F\-\',]+)\set\sal\.\s*\((\d{4})\)", current_citation)
        if et_al_year_match:
            author = et_al_year_match.group(1).split(",")[0].strip()
            year = et_al_year_match.group(2)
            new_citation = f"{author} et al. ({year})"
            print(f"  Cleaned up et al pattern -> {new_citation!r}")
        else:
            # Check for multiple authors listed: "Author1, Author2 & Author3 (year)"
            amp_year_match = re.search(r"^([A-Z\u00C0-\u024F][a-zA-Z\s\u00C0-\u024F\-\',&]+)\s*\((\d{4})\)", current_citation)
            if amp_year_match:
                authors_str = amp_year_match.group(1)
                year = amp_year_match.group(2)
                # Split by comma or &
                authors_parts = re.split(r",|&", authors_str)
                authors_parts = [a.strip() for a in authors_parts if a.strip()]
                if len(authors_parts) == 1:
                    new_citation = f"{authors_parts[0]} ({year})"
                elif len(authors_parts) == 2:
                    new_citation = f"{authors_parts[0]} & {authors_parts[1]} ({year})"
                else:
                    new_citation = f"{authors_parts[0]} et al. ({year})"
                print(f"  Cleaned up multiple authors pattern -> {new_citation!r}")
                
    # 4. Special manual corrections for specific IDs that are known
    if not new_citation:
        known_mappings = {
            "acs_apts_ref24_3dg_arrhenius_anchor": "ACS APTS (2022)",
            "mottram_nobrega_2002_furanone_bridge": "Mottram & Nobrega (2002)",
            "acs_jafc_3c02618_binding_prior": "ACS JAFC (2023)",
            "acs_jafc_0c01925_protein_binding_hierarchy": "ACS JAFC (2020)",
            "acs_ref3_spi_acrylamide_fast_kinetics": "ACS ResearchGate (2024)",
            "ref41_ppi_sulfur_binding": "JAFC (2019)",
            "scielo_brasil_aa_crosslink_hierarchy_anchor": "SciELO Brasil (2018)",
            "uspto_ptacts_2023_yeast_extract_anchor": "USPTO PTACTS (2023)",
            "vtechworks_2022_fava_hydrolysis": "VTechWorks (2022)",
        }
        if ref_id in known_mappings:
            new_citation = known_mappings[ref_id]
            print(f"  Assigned known mapping -> {new_citation!r}")
            
    # Apply change if resolved
    if new_citation:
        ref["citation"] = new_citation
        fixed_count += 1
    else:
        print(f"  Could not resolve {ref_id} citation!")

print(f"\nSuccessfully updated {fixed_count} out of {len(mismatches)} citation fields.")

# Save data back to registry
with open(REGISTRY_PATH, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)
print("Saved data back to benchmark_intake_registry.json.")
