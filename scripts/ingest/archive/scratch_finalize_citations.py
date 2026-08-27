import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# File paths
registry_path = ROOT / "data" / "lit" / "benchmark_intake_registry.json"
matrix_path = ROOT / "data" / "lit" / "slr_incorporation_matrix.json"
priors_path = ROOT / "data" / "lit" / "computational_priors.json"
calibrations_path = ROOT / "data" / "lit" / "process_state_calibrations.json"

# Helper function to recursively find and update entries in json
def update_recursive(data):
    count = 0
    if isinstance(data, dict):
        if "id" in data:
            if data["id"] == "jafc_2019_ref21_pea_gum_arabic_architecture_v1" and "source" in data:
                data["source"] = "Zha et al. (2019)"
                print("Updated jafc_2019_ref21_pea_gum_arabic_architecture_v1 in priors")
                count += 1
            elif data["id"] == "jafc_2022_melanoidin_thiol_binding" and "source" in data:
                data["source"] = "Gigl et al. (2021)"
                print("Updated jafc_2022_melanoidin_thiol_binding in priors")
                count += 1
        for k, v in data.items():
            count += update_recursive(v)
    elif isinstance(data, list):
        for item in data:
            count += update_recursive(item)
    return count

# 1. Update benchmark_intake_registry.json
with open(registry_path, "r", encoding="utf-8") as f:
    registry = json.load(f)

for ref in registry.get("eligible_references", []):
    if ref.get("id") == "jafc_2019_ref21_pea_gum_arabic_architecture_anchor":
        ref["citation"] = "Zha et al. (2019)"
        print("Updated jafc_2019_ref21_pea_gum_arabic_architecture_anchor in registry")
    elif ref.get("id") == "jafc_2022_melanoidin_thiol_binding":
        ref["citation"] = "Gigl et al. (2021)"
        print("Updated jafc_2022_melanoidin_thiol_binding in registry")

with open(registry_path, "w", encoding="utf-8") as f:
    json.dump(registry, f, indent=2, ensure_ascii=False)

# 2. Update slr_incorporation_matrix.json
with open(matrix_path, "r", encoding="utf-8") as f:
    matrix = json.load(f)

for entry in matrix.get("entries", []):
    if entry.get("paper_id") == "jafc_2022_melanoidin_thiol_binding":
        entry["citation"] = "Gigl et al. (2021)"
        print("Updated jafc_2022_melanoidin_thiol_binding in matrix")

with open(matrix_path, "w", encoding="utf-8") as f:
    json.dump(matrix, f, indent=2, ensure_ascii=False)

# 3. Update computational_priors.json
with open(priors_path, "r", encoding="utf-8") as f:
    priors = json.load(f)

updated_priors = update_recursive(priors)
print(f"Traversed and updated {updated_priors} priors.")

with open(priors_path, "w", encoding="utf-8") as f:
    json.dump(priors, f, indent=2, ensure_ascii=False)

# 4. Update process_state_calibrations.json
with open(calibrations_path, "r", encoding="utf-8") as f:
    calibrations = json.load(f)

for entry in calibrations.get("entries", []):
    if entry.get("id") == "jafc_2019_ref21_pea_gum_arabic_architecture_state":
        entry["source_citation"] = "Zha et al. (2019)"
        print("Updated jafc_2019_ref21_pea_gum_arabic_architecture_state in calibrations")

with open(calibrations_path, "w", encoding="utf-8") as f:
    json.dump(calibrations, f, indent=2, ensure_ascii=False)

print("All updates completed successfully!")
