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


def test_build_flavor_axis_summary_penalizes_expected_but_missing_furanones():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["alanine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
    )

    assert summary["furanone_expected"] == ["HEMF", "DMHF"]
    assert summary["furanone_observed"] == []
    assert summary["furanone_support_score"] == 0.0
    assert summary["furanone_penalty"] > 0.0
    assert summary["furanone_missing"] == ["HEMF", "DMHF"]


def test_build_flavor_axis_summary_keeps_native_soy_thiamine_inactive_by_default():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
    )

    assert summary["thiamine_pathway_active"] is False
    assert summary["thiamine_provenance_mode"] == "inactive"
    assert summary["thiamine_availability_source"] == "native_matrix_default_inactive"
    assert summary["thiamine_availability_explicit"] is False


def test_build_flavor_axis_summary_uses_explicit_thiamine_metadata_for_pbma_like_context():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "mft": {"compound": "2-Methyl-3-furanthiol (MFT)", "observable_ppb": 8.0},
        },
        sugars=["ribose"],
        amino_acids=["cysteine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=5.6,
        thiamine_availability={"available": True, "source": "pbma_fortified"},
    )

    assert summary["thiamine_pathway_active"] is True
    assert summary["thiamine_provenance_mode"] == "mixed_thiamine_plus_pentose"
    assert summary["thiamine_availability_source"] == "pbma_fortified"
    assert summary["thiamine_availability_explicit"] is True
    assert summary["thiamine_mft_fraction_estimate"] > 0.0


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
    assert "Ince et al. (2024)" in " ".join(hot_short["retention_reference_sources"])


def test_describe_retention_runtime_for_pea_markers_carries_structured_ph_release_provenance():
    hexanal = describe_retention_runtime(
        "Hexanal",
        protein_type="pea_iso",
        temperature_celsius=25.0,
        time_minutes=10.0,
        process_state="ambient_slurry",
    )
    pentylfuran = describe_retention_runtime(
        "2-Pentylfuran",
        protein_type="pea_iso",
        temperature_celsius=25.0,
        time_minutes=10.0,
        process_state="ambient_slurry",
    )

    assert hexanal["retention_runtime_mode"] == "direct_binding_plus_ph_release_reference"
    assert "Karolkowski et al. (2021)" in " ".join(hexanal["retention_reference_sources"])
    assert pentylfuran["retention_runtime_mode"] == "ph_release_reference"
    assert "Karolkowski et al. (2021)" in " ".join(pentylfuran["retention_reference_sources"])


def test_build_flavor_axis_summary_surfaces_secondary_strecker_reference_markers():
    summary = build_flavor_axis_summary(
        projection_metadata={
            "three_mb": {"compound": "3-Methylbutanal", "observable_ppb": 11.0},
            "benz": {"compound": "Benzaldehyde", "observable_ppb": 30.0},
        },
        sugars=["glucose"],
        amino_acids=["leucine", "phenylalanine"],
        additives=[],
        lipids=[],
        protein_type="soy_iso",
        pH=6.0,
    )

    assert summary["strecker_signals_ppb"]["3_methylbutanal"] == 11.0
    assert summary["strecker_signals_ppb"]["benzaldehyde"] == 30.0
    assert summary["strecker_reference_targets_ppb"]["3_methylbutanal"] > 0.0
    assert summary["strecker_reference_targets_ppb"]["benzaldehyde"] > 0.0