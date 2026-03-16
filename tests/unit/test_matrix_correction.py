import pytest
from src.matrix_correction import apply_matrix_correction, ProteinType, MATRIX_CORRECTIONS, resolve_matrix_correction

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
