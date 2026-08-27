import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REGISTRY_PATH = ROOT / "data" / "lit" / "benchmark_intake_registry.json"

with open(REGISTRY_PATH, "r", encoding="utf-8") as f:
    data = json.load(f)

# Regex matching:
# - Starts with an uppercase letter (including unicode accented letters like Š)
# - Followed by name characters (letters, spaces, dashes, apostrophes)
# - Optionally followed by " et al." or " & SecondAuthor" (allowing lowercase names/particles)
# - Ends with " (year)"
pattern = re.compile(r"^[A-Z\u00C0-\u024F][a-zA-Z\s\u00C0-\u024F\-\']*(?:\set\sal\.|&\s+[a-zA-Z\s\u00C0-\u024F\-\']+)?\s+\(\d{4}\)$")

mismatches = []
for ref in data.get("eligible_references", []):
    citation = ref.get("citation", "").strip()
    ref_id = ref.get("id", "")
    if not pattern.match(citation):
        mismatches.append({"id": ref_id, "citation": citation})

print(f"Total eligible references: {len(data.get('eligible_references', []))}")
print(f"Total mismatches: {len(mismatches)}")
if mismatches:
    print("Mismatched citations:")
    for m in mismatches:
        print(f"  - {m['id']}: {m['citation']!r}")
else:
    print("All citations are perfectly standardized!")
