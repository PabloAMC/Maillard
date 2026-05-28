import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
    data = json.load(f)

pattern = re.compile(r"^[A-Z][a-zA-Z\s\u00C0-\u024F\-\']+(?:(\set\sal\.|& [A-Z][a-zA-Z\s\u00C0-\u024F\-\']+))?\s\(\d{4}\)")

mismatches = []
for ref in data.get("eligible_references", []):
    citation = ref.get("citation", "")
    ref_id = ref.get("id", "")
    doi = ref.get("doi", "")
    if not pattern.match(citation):
        mismatches.append({"id": ref_id, "citation": citation, "doi": doi})

print(json.dumps(mismatches, indent=2))
