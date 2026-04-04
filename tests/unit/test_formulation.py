import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.formulation import Formulation  # noqa: E402


def test_formulation_from_dict_normalizes_legacy_keys():
    formulation = Formulation.from_dict(
        {
            "name": "legacy_formulation",
            "PH": "6.2",
            "temp": "145.0",
            "AW": "0.72",
            "protein": "soy_iso",
            "matrix": "extruded",
        }
    )

    assert formulation.ph == 6.2
    assert formulation.temperature == 145.0
    assert formulation.water_activity == 0.72
    assert formulation.protein_type == "soy_iso"
    assert formulation.matrix_type == "extruded"


def test_formulation_from_dict_prefers_canonical_keys_over_legacy_aliases():
    formulation = Formulation.from_dict(
        {
            "name": "canonical_formulation",
            "temperature": 130.0,
            "temp": 155.0,
            "water_activity": 0.61,
            "AW": 0.93,
            "protein_type": "pea_iso",
            "protein": "soy_iso",
        }
    )

    assert formulation.temperature == 130.0
    assert formulation.water_activity == 0.61
    assert formulation.protein_type == "pea_iso"
