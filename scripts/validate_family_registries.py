#!/usr/bin/env python3
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import data_paths  # noqa: E402


def validate():
    scope_path = data_paths.CHEMISTRY_FAMILY_SCOPE_REGISTRY
    ingestion_path = data_paths.FAMILY_INGESTION_PLAN

    if not scope_path.exists():
        print(f"Error: {scope_path} not found")
        sys.exit(1)
    if not ingestion_path.exists():
        print(f"Error: {ingestion_path} not found")
        sys.exit(1)

    with open(scope_path, "r") as f:
        scope = json.load(f)
    with open(ingestion_path, "r") as f:
        ingestion = json.load(f)

    scope_families = {f["family_id"]: f for f in scope["families"]}
    ingestion_families = {f["family_id"]: f for f in ingestion["families"]}

    errors = []

    # Check for missing families in scope
    for fid in ingestion_families:
        if fid not in scope_families:
            errors.append(f"Family '{fid}' found in ingestion plan but missing in scope registry.")

    # Check for missing families in ingestion
    for fid in scope_families:
        if fid not in ingestion_families:
            errors.append(f"Family '{fid}' found in scope registry but missing in ingestion plan.")

    # Validate strategic roles in scope registry
    allowed_roles = {
        "first_class_runtime_lane",
        "bounded_support_lane",
        "matrix_scope_lane",
        "guardrail_lane",
        "structural_gap_only",
        "computational_closure_candidate"
    }

    for fid, fdoc in scope_families.items():
        role = fdoc.get("strategic_role")
        if role not in allowed_roles:
            errors.append(f"Family '{fid}' has invalid strategic_role: '{role}'. Must be one of {allowed_roles}")

    if errors:
        print("Registry Validation Failed:")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)
    else:
        print(f"Registry Validation Passed: {len(scope_families)} families synchronized.")
        sys.exit(0)

if __name__ == "__main__":
    validate()
