from src.benchmark_validation import benchmark_to_conditions, benchmark_to_formulation


def test_benchmark_to_formulation_routes_support_family_precursors_into_additives():
    bench = {
        "benchmark_id": "family03_support_route",
        "precursors": {
            "Thiamine": {"concentration_mM": 0.1},
            "D-Ribose": {"concentration_mM": 1.0},
            "L-Cysteine": {"concentration_mM": 1.0},
        },
        "conditions": {
            "temp_C": 100.0,
            "ph": 5.5,
            "water_activity": 0.95,
            "time_min": 30.0,
        },
    }

    formulation = benchmark_to_formulation(bench)

    assert formulation["sugars"] == ["D-Ribose"]
    assert formulation["amino_acids"] == ["L-Cysteine"]
    assert formulation["additives"] == ["Thiamine"]


def test_benchmark_condition_alias_accepts_aw_for_conditions_and_formulation():
    bench = {
        "benchmark_id": "family03_aw_alias",
        "precursors": {
            "Thiamine": {"concentration_mM": 0.1},
            "D-Ribose": {"concentration_mM": 1.0},
            "L-Cysteine": {"concentration_mM": 1.0},
        },
        "conditions": {
            "temp_C": 100.0,
            "ph": 5.5,
            "aw": 0.98,
            "time_min": 30.0,
        },
    }

    conditions = benchmark_to_conditions(bench)
    formulation = benchmark_to_formulation(bench)

    assert conditions.water_activity == 0.98
    assert formulation["aw"] == 0.98