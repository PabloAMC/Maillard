from src.lipid_oxidation import build_lipid_input_proxy_loads, build_lipid_oxidation_context


def test_build_lipid_input_proxy_loads_falls_back_to_generic_oils_when_ratios_are_missing():
    loads = build_lipid_input_proxy_loads(["sunflower oil", "flax oil"], {})

    assert loads["sunflower oil"] == 1.0
    assert loads["flax oil"] == 1.0


def test_lipid_oxidation_context_surfaces_named_benchmark_ready_markers():
    context = build_lipid_oxidation_context(
        {"sunflower oil": 1.0},
        temp_C=120.0,
        time_min=15.0,
        water_activity=0.75,
    )

    assert context["dominant_marker"] == "Hexanal"
    assert "Hexanal" in context["generated_markers"]
    assert "2-Pentylfuran" in context["benchmark_ready_targets"]