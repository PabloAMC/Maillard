import pytest
from src.matrix_correction import (
    ACCESSIBILITY_LITERATURE_WINDOWS,
    DENATURATION_HEURISTICS,
    MATRIX_CORRECTIONS,
    ProteinType,
    apply_matrix_correction,
    build_matrix_explainability,
    classify_volatile_matrix_family,
    estimate_denaturation_state,
    get_volatile_class_retention_factor,
    get_accessibility_literature_window,
    resolve_compound_matrix_retention,
    resolve_effective_denaturation_state,
    resolve_matrix_correction,
)

def test_apply_matrix_correction_free_aa():
    """Verify that FREE_AMINO_ACID type applies no correction."""
    concs = {"hexanal": 100.0}
    amino_acids = {"lysine": 1.0, "cysteine": 1.0}
    
    vol, aa = apply_matrix_correction(concs, amino_acids, ProteinType.FREE_AMINO_ACID)
    
    assert vol["hexanal"] == 100.0
    assert aa["lysine"] == 1.0
    assert aa["cysteine"] == 1.0

def test_apply_matrix_correction_pea_native():
    """Verify native pea protein scaling (high burial)."""
    concs = {"hexanal": 100.0}
    amino_acids = {"lysine": 1.0, "cysteine": 1.0}
    
    # 0.0 denaturation = native state
    vol, aa = apply_matrix_correction(concs, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    
    corr = MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE]
    
    # Volatiles should be multiplied by the native retention factor.
    assert vol["hexanal"] == pytest.approx(corr.volatile_retention_native * 100.0)
    # Amino acids should be multiplied by accessibility factors
    assert aa["lysine"] == pytest.approx(corr.lysine_accessibility_native)
    assert aa["cysteine"] == pytest.approx(corr.cysteine_accessibility_native)

def test_apply_matrix_correction_pea_denatured():
    """Verify that denaturation increases accessibility."""
    concs = {"hexanal": 100.0}
    amino_acids = {"lysine": 1.0, "cysteine": 1.0}
    
    # 1.0 denaturation = fully denatured (approaches 1.0 accessibility)
    vol, aa = apply_matrix_correction(concs, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=1.0)
    
    corr = MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE]

    assert aa["lysine"] == pytest.approx(corr.lysine_accessibility_denatured)
    assert aa["cysteine"] == pytest.approx(corr.cysteine_accessibility_denatured)
    
    assert vol["hexanal"] == pytest.approx(
        MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE].volatile_retention_denatured * 100.0
    )

def test_apply_matrix_correction_unmapped_aa():
    """Verify that unmapped amino acids use the average factor."""
    amino_acids = {"leucine": 1.0}
    vol, aa = apply_matrix_correction({}, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    
    corr = MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE]
    expected = (corr.lysine_accessibility_native + corr.cysteine_accessibility_native) / 2.0
    assert aa["leucine"] == pytest.approx(expected)


def test_apply_matrix_correction_recognizes_canonical_smiles_for_lysine_and_cysteine():
    amino_acids = {
        "NCCCCC(N)C(=O)O": 1.0,
        "NC(CS)C(=O)O": 1.0,
        "NCC(=O)O": 1.0,
    }

    _vol, aa = apply_matrix_correction({}, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    corr = MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE]

    assert aa["NCCCCC(N)C(=O)O"] == pytest.approx(corr.lysine_accessibility_native)
    assert aa["NC(CS)C(=O)O"] == pytest.approx(corr.cysteine_accessibility_native)
    assert aa["NCC(=O)O"] == pytest.approx(
        (corr.lysine_accessibility_native + corr.cysteine_accessibility_native) / 2.0
    )


def test_resolve_matrix_correction_preserves_midpoint_baseline_and_interpolates_retention():
    corr = MATRIX_CORRECTIONS[ProteinType.PEA_ISOLATE]

    native = resolve_matrix_correction(ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    midpoint = resolve_matrix_correction(ProteinType.PEA_ISOLATE, denaturation_state=0.5)
    denatured = resolve_matrix_correction(ProteinType.PEA_ISOLATE, denaturation_state=1.0)

    assert native.lysine_accessibility == pytest.approx(corr.lysine_accessibility_native)
    assert midpoint.lysine_accessibility == pytest.approx(corr.lysine_accessibility)
    assert denatured.lysine_accessibility == pytest.approx(corr.lysine_accessibility_denatured)
    assert native.cysteine_accessibility == pytest.approx(corr.cysteine_accessibility_native)
    assert midpoint.cysteine_accessibility == pytest.approx(corr.cysteine_accessibility)
    assert denatured.cysteine_accessibility == pytest.approx(corr.cysteine_accessibility_denatured)
    assert native.volatile_retention == pytest.approx(corr.volatile_retention_native)
    assert midpoint.volatile_retention == pytest.approx(corr.volatile_retention)
    assert denatured.volatile_retention == pytest.approx(corr.volatile_retention_denatured)
    assert native.lysine_accessibility < midpoint.lysine_accessibility < denatured.lysine_accessibility
    assert native.cysteine_accessibility < midpoint.cysteine_accessibility < denatured.cysteine_accessibility
    assert denatured.cysteine_accessibility < 0.1


@pytest.mark.parametrize(
    "protein_type",
    [ProteinType.PEA_ISOLATE, ProteinType.SOY_ISOLATE],
)
def test_isolate_accessibility_midpoints_stay_within_literature_windows(protein_type):
    corr = MATRIX_CORRECTIONS[protein_type]
    window = ACCESSIBILITY_LITERATURE_WINDOWS[protein_type]
    resolved = resolve_matrix_correction(protein_type, denaturation_state=0.5)

    assert window.lysine_min <= corr.lysine_accessibility <= window.lysine_max
    assert window.cysteine_min <= corr.cysteine_accessibility <= window.cysteine_max
    assert window.lysine_min <= resolved.lysine_accessibility <= window.lysine_max
    assert window.cysteine_min <= resolved.cysteine_accessibility <= window.cysteine_max


@pytest.mark.parametrize(
    ("protein_type", "native_limit", "denatured_limit"),
    [
        (ProteinType.PEA_ISOLATE, 0.35, 0.45),
        (ProteinType.SOY_ISOLATE, 0.40, 0.50),
    ],
)
def test_native_and_denatured_endpoints_remain_in_calibrated_accessibility_envelope(
    protein_type,
    native_limit,
    denatured_limit,
):
    native = resolve_matrix_correction(protein_type, denaturation_state=0.0)
    denatured = resolve_matrix_correction(protein_type, denaturation_state=1.0)
    window = get_accessibility_literature_window(protein_type)

    assert window is not None
    assert native.lysine_accessibility < native_limit
    assert denatured.lysine_accessibility <= denatured_limit
    assert window.cysteine_min <= native.cysteine_accessibility <= window.cysteine_max
    assert window.cysteine_min <= denatured.cysteine_accessibility <= window.cysteine_max


@pytest.mark.parametrize("denaturation_state", [0.0, 0.5, 1.0])
def test_legume_isolate_accessibility_hierarchy_stays_explicit(denaturation_state):
    free = resolve_matrix_correction(ProteinType.FREE_AMINO_ACID, denaturation_state=denaturation_state)
    soy = resolve_matrix_correction(ProteinType.SOY_ISOLATE, denaturation_state=denaturation_state)
    pea = resolve_matrix_correction(ProteinType.PEA_ISOLATE, denaturation_state=denaturation_state)

    assert free.lysine_accessibility > soy.lysine_accessibility > pea.lysine_accessibility
    assert free.cysteine_accessibility > soy.cysteine_accessibility > pea.cysteine_accessibility
    assert free.volatile_retention > soy.volatile_retention > pea.volatile_retention


@pytest.mark.parametrize("denaturation_state", [0.0, 0.5, 1.0])
def test_concentrates_remain_more_buried_than_isolates_within_each_legume_family(denaturation_state):
    soy_conc = resolve_matrix_correction(ProteinType.SOY_CONCENTRATE, denaturation_state=denaturation_state)
    soy_iso = resolve_matrix_correction(ProteinType.SOY_ISOLATE, denaturation_state=denaturation_state)
    pea_conc = resolve_matrix_correction(ProteinType.PEA_CONCENTRATE, denaturation_state=denaturation_state)
    pea_iso = resolve_matrix_correction(ProteinType.PEA_ISOLATE, denaturation_state=denaturation_state)

    assert soy_conc.lysine_accessibility < soy_iso.lysine_accessibility
    assert soy_conc.cysteine_accessibility < soy_iso.cysteine_accessibility
    assert soy_conc.volatile_retention < soy_iso.volatile_retention

    assert pea_conc.lysine_accessibility < pea_iso.lysine_accessibility
    assert pea_conc.cysteine_accessibility < pea_iso.cysteine_accessibility
    assert pea_conc.volatile_retention < pea_iso.volatile_retention


def test_literature_windows_are_exposed_only_for_quantitatively_anchored_isolates():
    assert get_accessibility_literature_window(ProteinType.PEA_ISOLATE) is not None
    assert get_accessibility_literature_window(ProteinType.SOY_ISOLATE) is not None
    assert get_accessibility_literature_window(ProteinType.PEA_CONCENTRATE) is None
    assert get_accessibility_literature_window(ProteinType.FREE_AMINO_ACID) is None


def test_denaturation_heuristics_are_exposed_for_matrix_types_only():
    assert DENATURATION_HEURISTICS[ProteinType.PEA_ISOLATE].midpoint_celsius > 70.0
    assert DENATURATION_HEURISTICS[ProteinType.SOY_ISOLATE].midpoint_celsius > 70.0
    assert ProteinType.FREE_AMINO_ACID not in DENATURATION_HEURISTICS


def test_estimate_denaturation_state_tracks_temperature_time_and_ph_for_pea_isolate():
    cool = estimate_denaturation_state(ProteinType.PEA_ISOLATE, 40.0, time_minutes=10.0, pH=6.0)
    warm = estimate_denaturation_state(ProteinType.PEA_ISOLATE, 90.0, time_minutes=30.0, pH=6.0)
    hot = estimate_denaturation_state(ProteinType.PEA_ISOLATE, 140.0, time_minutes=30.0, pH=5.5)
    acidic = estimate_denaturation_state(ProteinType.PEA_ISOLATE, 90.0, time_minutes=30.0, pH=4.5)

    assert 0.02 <= cool <= 0.15
    assert 0.55 <= warm <= 0.80
    assert hot >= 0.95
    assert cool < warm < hot
    assert acidic > warm


def test_estimate_denaturation_state_keeps_soy_slightly_more_structured_than_pea_at_same_conditions():
    pea = estimate_denaturation_state(ProteinType.PEA_ISOLATE, 90.0, time_minutes=30.0, pH=6.0)
    soy = estimate_denaturation_state(ProteinType.SOY_ISOLATE, 90.0, time_minutes=30.0, pH=6.0)

    assert soy < pea
    assert 0.45 <= soy <= 0.75


def test_resolve_effective_denaturation_state_respects_explicit_override():
    inferred = resolve_effective_denaturation_state(
        ProteinType.PEA_ISOLATE,
        temperature_celsius=140.0,
        time_minutes=30.0,
        pH=5.5,
    )
    forced = resolve_effective_denaturation_state(
        ProteinType.PEA_ISOLATE,
        temperature_celsius=140.0,
        time_minutes=30.0,
        pH=5.5,
        explicit_denaturation_state=0.2,
    )

    assert inferred > 0.9
    assert forced == pytest.approx(0.2)


def test_compound_class_retention_penalizes_aldehydes_more_than_furans_and_sulfur_volatiles():
    aldehyde = get_volatile_class_retention_factor("Hexanal", ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    furan = get_volatile_class_retention_factor("Furfural", ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    sulfur = get_volatile_class_retention_factor("2-Furfurylthiol (FFT)", ProteinType.PEA_ISOLATE, denaturation_state=0.0)

    assert classify_volatile_matrix_family("Hexanal") == "aldehyde"
    assert classify_volatile_matrix_family("Furfural") == "furan"
    assert classify_volatile_matrix_family("2-Furfurylthiol (FFT)") == "sulfur"
    assert aldehyde < furan < sulfur


def test_resolve_compound_matrix_retention_stays_bounded_and_class_aware():
    aldehyde = resolve_compound_matrix_retention("Hexanal", ProteinType.SOY_ISOLATE, denaturation_state=0.0)
    alcohol = resolve_compound_matrix_retention("1-Hexanol", ProteinType.SOY_ISOLATE, denaturation_state=0.0)

    assert 0.01 <= aldehyde <= 1.0
    assert 0.01 <= alcohol <= 1.0
    assert aldehyde < alcohol


def test_build_matrix_explainability_surfaces_effective_accessibility_context():
    payload = build_matrix_explainability(
        protein_type=ProteinType.PEA_ISOLATE,
        effective_denaturation_state=0.65,
        temperature_celsius=95.0,
        time_minutes=30.0,
        pH=5.5,
    )

    assert payload["protein_type"] == "pea_iso"
    assert payload["effective_denaturation_state"] == pytest.approx(0.65)
    assert payload["lysine_accessibility"] > payload["cysteine_accessibility"]
    assert payload["literature_window"] is not None
