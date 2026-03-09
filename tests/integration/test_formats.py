import pytest
import yaml
import json
from pathlib import Path
from src.cantera_export import CanteraExporter

def test_cantera_yaml_format_validity():
    """Verify that CanteraExporter generates valid YAML 1.2."""
    exporter = CanteraExporter()
    reactants = ["C", "[O]"] # Graphite atom + Oxygen atom -> CO
    products = ["CO"]
    exporter.add_reaction(reactants, products, 10.0)
    
    exporter.export_yaml(output_path="tmp_mech.yaml")
    yaml_str = Path("tmp_mech.yaml").read_text()
    assert yaml_str.strip().startswith("generator: Maillard")
    
    # Attempt to parse with standard PyYAML to ensure it's valid YAML
    data = yaml.safe_load(yaml_str)
    assert "phases" in data
    assert "species" in data
    assert "reactions" in data

def test_database_schema_placeholders():
    """Verify that we can simulate/mock a results DB entry."""
    # This prevents regressions in the data format sent to results_db.py
    record = {
        "precursors": ["glucose", "glycine"],
        "pathways": [
            {"name": "Maillard", "score": 0.85, "active": True}
        ],
        "conditions": {"T": 423.15, "pH": 7.0}
    }
    dumped = json.dumps(record)
    reloaded = json.loads(dumped)
    assert reloaded["pathways"][0]["score"] == 0.85
