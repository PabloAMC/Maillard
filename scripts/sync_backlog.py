import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKLOG_PATH = ROOT / "data" / "lit" / "deep_research_backlog.json"
REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

def clean_string(s):
    if not s:
        return ""
    return str(s).strip().lower().replace(" ", "").replace("&", "and").replace(",", "").replace("-", "").replace("_", "").replace(".", "").replace("/", "").replace(":", "").replace("(", "").replace(")", "")

def main():
    print(f"Loading {BACKLOG_PATH.name}...")
    with open(BACKLOG_PATH, "r", encoding="utf-8") as f:
        backlog = json.load(f)
        
    print(f"Loading {REGISTRY_PATH.name}...")
    with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
        registry = json.load(f)

    # Map different identifiers to registry entry
    reg_map = {}
    
    # Custom manually-defined mapping overrides for special cases
    custom_map = {
        "jafcdoi101021jf034037p2003": "blank_devaud_grosch_2003_g6p_hdmf_prior",
        "foods121019672023": "foods_2023_cml_cel_proxy_benchmark",
        "foods202312101967": "foods_2023_cml_cel_proxy_benchmark",
        "foods2023121019672023": "foods_2023_cml_cel_proxy_benchmark",
    }
    
    for ref in registry.get("eligible_references", []):
        ref_id = ref["id"]
        cit = ref.get("citation", "")
        doi = ref.get("doi", "")
        aliases = ref.get("citation_aliases", [])
        
        reg_map[clean_string(ref_id)] = ref
        reg_map[clean_string(cit)] = ref
        if doi:
            reg_map[clean_string(doi)] = ref
        for alias in aliases:
            reg_map[clean_string(alias)] = ref

    print("Syncing backlog items...")
    updated_count = 0
    warnings_count = 0
    
    for item in backlog.get("items", []):
        cit = item.get("citation", "")
        doi = item.get("doi", "")
        status = item.get("status", "BACKLOG")
        reg_id = item.get("registry_id", "")
        
        # Try matching
        matched_ref = None
        cleaned_cit = clean_string(cit)
        cleaned_doi = clean_string(doi)
        
        # 1. Try custom map first
        if cleaned_cit in custom_map:
            target_id = custom_map[cleaned_cit]
            for ref in registry.get("eligible_references", []):
                if ref["id"] == target_id:
                    matched_ref = ref
                    break
        
        # 2. Try clean citation
        if not matched_ref and cleaned_cit in reg_map:
            matched_ref = reg_map[cleaned_cit]
            
        # 3. Try clean DOI
        if not matched_ref and cleaned_doi and cleaned_doi in reg_map:
            matched_ref = reg_map[cleaned_doi]
            
        # 4. Substring match try
        if not matched_ref:
            for key, val in reg_map.items():
                if len(key) > 8 and len(cleaned_cit) > 8:
                    if cleaned_cit in key or key in cleaned_cit:
                        matched_ref = val
                        break
                        
        if matched_ref:
            target_id = matched_ref["id"]
            artifact_count = len(matched_ref.get("runtime_artifacts", []))
            
            # Update backlog entry if needed
            needs_update = False
            if status != "RUNTIME_BOUND":
                item["status"] = "RUNTIME_BOUND"
                needs_update = True
            if item.get("registry_id") != target_id:
                item["registry_id"] = target_id
                needs_update = True
            if item.get("runtime_artifact_count") != artifact_count:
                item["runtime_artifact_count"] = artifact_count
                needs_update = True
                
            if needs_update:
                print(f"  [UPDATED] '{cit}' -> status: RUNTIME_BOUND, registry_id: '{target_id}', artifacts: {artifact_count}")
                updated_count += 1
        else:
            if status == "RUNTIME_BOUND":
                print(f"  [WARNING] '{cit}' is marked RUNTIME_BOUND but not found in the registry!")
                warnings_count += 1

    if updated_count > 0:
        with open(BACKLOG_PATH, "w", encoding="utf-8") as f:
            json.dump(backlog, f, indent=2, ensure_ascii=False)
        print(f"Successfully synced backlog: updated {updated_count} items.")
    else:
        print("Backlog is already in sync with the registry.")
        
    if warnings_count > 0:
        print(f"Sync finished with {warnings_count} warnings.")

if __name__ == "__main__":
    main()
