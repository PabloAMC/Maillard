import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.benchmark_validation import evaluate_benchmark
from src.headspace import HeadspaceModel


def test_pouvreau_pea_ph_family_is_reproduced_as_an_acidic_release_trend():
    """
    Pouvreau 2021 reports roughly 55-65% higher headspace release at pH 4.5
    vs pH 6.5 for pea-isolate aldehydes/furans.

    We validate the trend directly at the headspace layer rather than forcing an
    absolute benchmark file with incomplete experimental metadata.
    """
    model = HeadspaceModel()
    matrix = {
        "Hexanal": 205.0,
        "Nonanal": 75.0,
        "2-Pentylfuran": 310.0,
        "2,5-Dimethylpyrazine": 10.0,
    }

    air_acid = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=4.5)
    air_less_acid = model.predict_headspace(matrix, 40.0, protein_type="pea_iso", pH=6.5)

    for compound in ["Hexanal", "Nonanal", "2-Pentylfuran"]:
        ratio = air_acid[compound] / air_less_acid[compound]
        assert 1.5 <= ratio <= 1.7, f"{compound} acidic/basic ratio drifted to {ratio:.3f}"

    pyrazine_ratio = air_acid["2,5-Dimethylpyrazine"] / air_less_acid["2,5-Dimethylpyrazine"]
    assert pyrazine_ratio == pytest.approx(1.0, rel=0.02)


def test_matrix_only_pratap_singh_baselines_remain_stable_at_reference_ph():
    pea = evaluate_benchmark(ROOT / "data" / "benchmarks" / "pea_isolate_40C_PratapSingh2021.json")
    soy = evaluate_benchmark(ROOT / "data" / "benchmarks" / "soy_isolate_40C_PratapSingh2021.json")

    for evaluation in [pea, soy]:
        assert evaluation.supported is True
        for comparison in evaluation.comparisons:
            assert comparison.ratio <= 1.01