import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

files_to_fix = [
    ROOT / "data" / "lit" / "benchmark_intake_registry.json",
    ROOT / "data" / "lit" / "slr_incorporation_matrix.json",
    ROOT / "data" / "lit" / "safety_reference_payloads.json",
    ROOT / "data" / "lit" / "deep_research_backlog.json",
    ROOT / "scripts" / "ingest_batch_protein_damage_markers.py",
    ROOT / "scripts" / "ingest_deep_research_campaigns.py"
]

def main():
    for filepath in files_to_fix:
        if not filepath.exists():
            print(f"Skipping non-existent: {filepath.name}")
            continue
            
        print(f"Processing: {filepath.name}")
        content = filepath.read_text(encoding="utf-8")
        
        # Replace citations
        new_content = content.replace("Urugo et al. (2024)", "Mandelli et al. (2025)")
        
        # Also clean up any accidental 13030409 DOI references for Urugo if they exist
        if "foods13030409" in new_content and "urugo" in filepath.name:
            new_content = new_content.replace("10.3390/foods13030409", "10.3390/foods14193295")
            
        if new_content != content:
            filepath.write_text(new_content, encoding="utf-8")
            print(f"  -> Fixed citation in {filepath.name}")
        else:
            print("  -> No changes needed")

if __name__ == "__main__":
    main()
