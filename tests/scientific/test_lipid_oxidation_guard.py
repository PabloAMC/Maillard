import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_free_amino_acid_benchmark_does_not_inject_lipid_oxidation_products(hofmann_meaty_result):
    result = hofmann_meaty_result

    assert result.predicted_ppb.get("Hexanal") is None
    assert result.predicted_ppb.get("Nonanal") is None
    assert result.predicted_ppb.get("2-Pentylfuran") is None


def test_free_amino_acid_benchmark_predicted_ppb_excludes_input_precursors(mottram_meaty_result):
    result = mottram_meaty_result

    assert result.predicted_ppb.get("D-Ribose") is None
    assert result.predicted_ppb.get("L-Cysteine") is None