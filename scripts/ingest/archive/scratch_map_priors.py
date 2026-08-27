import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PRIORS_PATH = ROOT / "data" / "lit" / "computational_priors.json"

def main():
    with open(PRIORS_PATH) as f:
        data = json.load(f)
    
    # Let's print out all ids / keys in the sections to find kinetic ones
    for sec, payload in data.items():
        if sec in ["generated_at", "description", "section_family_metadata"]:
            continue
        print(f"=== {sec} ===")
        lst = payload if isinstance(payload, list) else payload.get("entries", [])
        for entry in lst:
            id_val = entry.get("id") or entry.get("reaction_key")
            # Look for energy / rate keys
            keys = [k for k in entry.keys() if any(term in k.lower() for term in ["ea", "barrier", "rate", "k_"])]
            if keys:
                print(f"  {id_val}:")
                for k in keys:
                    print(f"    {k}: {entry[k]}")

if __name__ == "__main__":
    main()
