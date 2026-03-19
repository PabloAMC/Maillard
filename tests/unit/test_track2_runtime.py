from src.literature_runtime import build_flavor_axis_summary, describe_retention_runtime


def test_build_flavor_axis_summary_penalizes_sulfur_positive_but_strecker_poor_profile():
    projection_metadata = {
        "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 12.0},
        "furfural": {"compound": "Furfural", "observable_ppb": 200.0},
        "pyraz": {"compound": "2,5-Dimethylpyrazine", "observable_ppb": 40.0},
    }

    summary = build_flavor_axis_summary(
        projection_metadata=projection_metadata,
        sugars=["fructose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=8.5,
    )

    assert summary["strecker_balance_score"] < 0.4
    assert summary["strecker_gap_penalty"] > 0.0
    assert summary["pyrazine_burden"] > summary["pyrazine_signal_ppb"]
    assert summary["pyrazine_penalty"] > 0.0


def test_describe_retention_runtime_for_soy_hexanal_is_temperature_and_time_sensitive():
    cool = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=25.0,
        time_minutes=2.0,
        process_state="ambient_slurry",
    )
    hot_short = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=95.0,
        time_minutes=2.0,
        process_state="heated_matrix",
    )
    hot_long = describe_retention_runtime(
        "Hexanal",
        protein_type="soy_iso",
        temperature_celsius=95.0,
        time_minutes=30.0,
        process_state="heated_matrix",
    )

    assert hot_short["dynamic_retention_factor"] > cool["dynamic_retention_factor"]
    assert hot_long["dynamic_retention_factor"] < hot_short["dynamic_retention_factor"]
    assert "Ince 2024" in " ".join(hot_short["retention_reference_sources"])