import pytest
from src.matrix_correction import apply_matrix_correction, ProteinType, MATRIX_CORRECTIONS

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
    
    # Volatiles should be multiplied by retention factor (0.50)
    assert vol["hexanal"] == pytest.approx(50.0)
    # Amino acids should be multiplied by accessibility factors
    assert aa["lysine"] == pytest.approx(0.40)
    assert aa["cysteine"] == pytest.approx(0.25)

def test_apply_matrix_correction_pea_denatured():
    """Verify that denaturation increases accessibility."""
    concs = {"hexanal": 100.0}
    amino_acids = {"lysine": 1.0, "cysteine": 1.0}
    
    # 1.0 denaturation = fully denatured (approaches 1.0 accessibility)
    vol, aa = apply_matrix_correction(concs, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=1.0)
    
    assert aa["lysine"] == pytest.approx(1.0)
    assert aa["cysteine"] == pytest.approx(1.0)
    
    # Volatile retention stays the same as it's a matrix property, not state-dependent in this version
    assert vol["hexanal"] == pytest.approx(50.0)

def test_apply_matrix_correction_unmapped_aa():
    """Verify that unmapped amino acids use the average factor."""
    amino_acids = {"leucine": 1.0}
    vol, aa = apply_matrix_correction({}, amino_acids, ProteinType.PEA_ISOLATE, denaturation_state=0.0)
    
    expected = (0.40 + 0.25) / 2.0
    assert aa["leucine"] == pytest.approx(expected)
