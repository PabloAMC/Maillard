"""
tests/unit/test_sensory_full.py

Verifies the unified SensoryDatabase and psychophysical mixing.
"""

import pytest
from src.sensory import SensoryDatabase, SensoryPredictor  # noqa: E402
from src.headspace import HeadspaceModel  # noqa: E402

def test_sensory_database_loading():
    """Verify that we load compounds from multiple YAML files."""
    db = SensoryDatabase()
    
    # Check FFT (from desirable_targets.yml)
    fft = db.find_entry("2-Furfurylthiol (FFT)")
    assert fft is not None
    assert fft["threshold_ppm"] == 0.00001 # 0.01 ppb / 1000
    
    # Check Hexanal (from off_flavour_targets.yml)
    hexanal = db.find_entry("Hexanal")
    assert hexanal is not None
    assert hexanal["threshold_ppm"] == 0.0045 # 4.5 ppb / 1000
    
    # Check CML (from toxic_markers.yml - might not have threshold)
    cml = db.find_entry("Nε-(Carboxymethyl)lysine (CML)")
    assert cml is not None

def test_psychophysical_scaling():
    """Verify Stevens' Law scaling: Intensity = (C/Threshold)^0.5"""
    predictor = SensoryPredictor()
    
    # threshold = 1.0 ppb, conc = 100.0 ppb -> OAV = 100 -> Intensity = 10
    mock_entry = {
        "name": "TestCompound",
        "threshold_ppb": 1.0,
        "descriptors": ["test"],
        "smiles": "C_TEST"
    }
    predictor.db.compounds["TestCompound"] = mock_entry
    predictor.db.smiles_map["C_TEST"] = mock_entry
    predictor.db.chemical_to_descriptor["C_TEST"] = {"odt": 1.0, "descriptor": "test"}
    
    # predict_profile expects SMILES 
    profile = predictor.predict_profile({"C_TEST": 100.0})
    assert profile["TestCompound"][0] == pytest.approx(10.0)

def test_radar_aggregation():
    """Verify that intensities are correctly grouped into categories."""
    predictor = SensoryPredictor()
    
    # Mock concentrations in ppb
    # 2-Furfurylthiol (FFT) threshold ~0.01 ppb
    # Hexanal threshold ~4.5 ppb
    mock_conc = {
        "2-Furfurylthiol (FFT)": 0.1,  # OAV = 10, Intensity = sqrt(10) ~ 3.16
    }
    
    radar = predictor.get_radar_data(mock_conc)
    
    # radar[tag] is now (score, uncertainty)
    assert radar["roasted"][0] == pytest.approx(3.162, rel=1e-3)

def test_sensory_headspace_integration():
    """Verify end-to-end headspace aware sensory scoring."""
    hs = HeadspaceModel()
    predictor = SensoryPredictor(headspace=hs)
    
    # 1000 ppb (1 ppm) Hexanal in 10% fat matrix
    # Expected intensity ~0.27 (see calculations in original test)
    profile = predictor.predict_profile({"Hexanal": 1000.0}, temp_c=25.0, fat_fraction=0.1)
    # Using .get because sub-threshold compounds might be omitted from profile
    assert profile.get("Hexanal", (0.0, 0.0))[0] < 0.4
    
    # Same concentration, but 0% fat -> much higher intensity.
    # Was ~1.82 when hexanal Kaw_25c was 0.015; the audit-remediation correction of
    # data/lit/henry_constants.yml to the Sander/Brockbank value 0.009 moves it to ~1.41.
    # The assertion pins the *relationship* (fat suppresses release), not the constant.
    profile_pure = predictor.predict_profile({"Hexanal": 1000.0}, temp_c=25.0, fat_fraction=0.0)
    assert profile_pure.get("Hexanal", (0.0, 0.0))[0] > 1.2
    assert profile_pure.get("Hexanal", (0.0, 0.0))[0] > 3.0 * profile.get("Hexanal", (0.0, 0.0))[0]


# ---------------------------------------------------------------------------
# Audit remediation (wave G5, module-hunt TIER 4) regression tests.
# ---------------------------------------------------------------------------


def _mock_predictor_with_test_compound(threshold_ppb: float = 1.0) -> SensoryPredictor:
    predictor = SensoryPredictor()
    mock_entry = {
        "name": "TestCompound",
        "threshold_ppb": threshold_ppb,
        "descriptors": ["test"],
        "smiles": "C_TEST",
    }
    predictor.db.compounds["TestCompound"] = mock_entry
    predictor.db.smiles_map["C_TEST"] = mock_entry
    return predictor


def test_masking_uses_consistent_token_aware_family_membership():
    """
    TIER 4.1: 'Intense meaty' must satisfy the 'meaty' family.

    The old code tested exact list membership on descriptors for the mask application and
    substring on the name for the totals, so most meaty compounds -- MFT among them --
    were never masked at all.
    """
    predictor = SensoryPredictor()
    meaty = predictor._family_members("meaty")

    # Every compound the radar tags as meaty is now a masking family member...
    for display_name in predictor.db.tags["meaty"]:
        entry = predictor.db.find_entry(display_name)
        assert entry is not None, display_name
        assert entry["name"] in meaty, display_name

    # ...including the ones whose free-text descriptor only *contains* the token.
    assert "2-Methyl-3-furanthiol (MFT)" in meaty      # "Intense meaty, boiled meat, ..."
    assert "Dimethyl disulfide" in meaty               # "Onion-like, meaty base note"
    assert "Hydrogen Sulfide" in meaty                 # "Sulfurous, boiled egg, meaty base"

    # Token matching, not substring: 'meat' must not drag in unrelated words.
    assert "Hexanal" not in meaty
    assert predictor._family_members("beany") >= {"Hexanal", "2-Pentylfuran", "Nonanal"}


def test_beany_masking_actually_reduces_meaty_intensity_end_to_end():
    """TIER 4.1: hexanal (beany) must measurably suppress MFT (meaty)."""
    predictor = SensoryPredictor()

    mft_alone = predictor.predict_profile({"2-Methyl-3-furanthiol (MFT)": 1.0})
    with_beany = predictor.predict_profile(
        {"2-Methyl-3-furanthiol (MFT)": 1.0, "Hexanal": 500.0}
    )

    masked = with_beany["2-Methyl-3-furanthiol (MFT)"][0]
    unmasked = mft_alone["2-Methyl-3-furanthiol (MFT)"][0]
    assert masked < unmasked
    # Uncertainty is scaled by the same factor.
    assert with_beany["2-Methyl-3-furanthiol (MFT)"][1] == pytest.approx(
        mft_alone["2-Methyl-3-furanthiol (MFT)"][1] * masked / unmasked
    )
    # The beany compound itself is not masked.
    assert with_beany["Hexanal"][0] == pytest.approx(
        predictor.predict_profile({"Hexanal": 500.0})["Hexanal"][0]
    )
    # Mask is bounded by k_mask = 0.35.
    assert masked >= 0.65 * unmasked


def test_synergy_requires_two_distinct_compounds():
    """
    TIER 4.2: 2,5-dimethylpyrazine used to satisfy both halves of the
    ('2,5-dimethylpyrazine', 'pyrazine') pair on its own, inflating itself by 30%.
    """
    predictor = SensoryPredictor()

    # 180000 ppb / 1800 ppb ODT = OAV 100 -> intensity 10 for a single compound.
    solo = predictor.get_radar_data({"2,5-Dimethylpyrazine": 180000.0})["roasted"][0]
    assert solo == pytest.approx(10.0), "self-synergy must not boost a lone compound"

    # A genuine second pyrazine partner does earn the 1.3x boost.
    pair = predictor.get_radar_data(
        {"2,5-Dimethylpyrazine": 180000.0, "2,3-Dimethylpyrazine": 250000.0}
    )["roasted"][0]
    unboosted = pow(pow(10.0, 1.2) + pow(10.0, 1.2), 1.0 / 1.2)
    assert pair == pytest.approx(unboosted * 1.3)


def test_compound_without_threshold_is_not_scored_with_a_fabricated_odt():
    """
    TIER 4.3: toxic markers carry no odour threshold. They must be skipped with a note,
    never assigned the fabricated 0.1 ppb default (which made them wildly odour-active).
    """
    predictor = SensoryPredictor()
    toxic = predictor.db.find_entry("Acrylamide")
    assert toxic is not None
    assert toxic["threshold_ppm"] is None

    profile = predictor.predict_profile({"Acrylamide": 1_000_000.0})
    assert "Acrylamide" not in profile

    reasons = {row["reason"] for row in predictor.unscored_compounds}
    assert "no_odour_threshold" in reasons
    assert any(row.get("name") == "Acrylamide" for row in predictor.unscored_compounds)

    # And every toxic marker behaves the same way (all 8 were affected).
    toxic_names = [
        name for name, entry in predictor.db.compounds.items()
        if entry["source"] == "toxic_markers.yml"
    ]
    assert len(toxic_names) >= 8
    profile = predictor.predict_profile({name: 1_000_000.0 for name in toxic_names})
    assert profile == {}


def test_zero_threshold_is_honoured_rather_than_treated_as_missing():
    """TIER 4.3 (falsy-0.0 bug): a stored 0.0 threshold must not fall through to a default."""
    predictor = SensoryPredictor()
    entry = {"name": "ZeroThreshold", "threshold_ppm": 0.0, "descriptors": ["test"], "smiles": None}
    assert SensoryPredictor._resolve_odt_ppb(entry) == 0.0
    predictor.db.compounds["ZeroThreshold"] = entry
    # 0.0 is not a usable divisor -> skipped with a note, not silently defaulted.
    profile = predictor.predict_profile({"ZeroThreshold": 5.0})
    assert "ZeroThreshold" not in profile
    assert predictor.unscored_compounds[0]["reason"] == "no_odour_threshold"


def test_find_entry_requires_exact_or_alias_match():
    """
    TIER 4.4: the fuzzy fallback resolved any partial name to an arbitrary first hit and
    filed the caller's concentration under that compound's name.
    """
    db = SensoryDatabase()

    # Declared aliases still resolve.
    assert db.find_entry("MFT")["name"] == "2-Methyl-3-furanthiol (MFT)"
    assert db.find_entry("2-Furfurylthiol")["name"] == "2-Furfurylthiol (FFT)"
    assert db.find_entry("methional")["name"] == "Methional (3-(methylthio)propanal)"
    assert db.find_entry("3-(methylthio)propanal")["name"] == "Methional (3-(methylthio)propanal)"
    assert db.find_entry("isobutyraldehyde")["name"] == "2-Methylpropanal (isobutyraldehyde)"
    assert db.find_entry("hexanal")["name"] == "Hexanal"

    # Ambiguous partial names no longer silently become one specific compound.
    assert db.find_entry("dimethylpyrazine") is None
    assert db.find_entry("pyrazine") is None
    assert db.find_entry("thiazole") is None
    assert db.find_entry("Totally Made Up Compound") is None
    assert "dimethylpyrazine" in db.unresolved_identifiers


def test_unresolved_identifier_does_not_contaminate_another_compound():
    """TIER 4.4: an unknown name must not have its concentration filed under a real one."""
    predictor = SensoryPredictor()
    profile = predictor.predict_profile({"dimethylpyrazine": 1.0e7})
    assert profile == {}
    assert predictor.unscored_compounds[0]["reason"] == "unresolved_identifier"


def test_intensity_is_continuous_at_the_odour_threshold():
    """
    TIER 4.5: conc <= ODT scored 0 while ODT + eps scored ~1.0 -- a unit-sized jump in the
    optimizer objective at every threshold.
    """
    predictor = _mock_predictor_with_test_compound(threshold_ppb=1.0)

    just_below = predictor.predict_profile({"C_TEST": 0.999999})
    assert just_below == {}

    for oav, bound in [(1.000001, 1e-5), (1.001, 2e-2), (1.01, 2e-1)]:
        value = predictor.predict_profile({"C_TEST": oav})["TestCompound"][0]
        assert 0.0 < value < bound, (oav, value)

    # Monotone across the onset band and into the suprathreshold region.
    grid = [1.0001, 1.01, 1.1, 1.5, 1.9, 2.0, 4.0, 100.0]
    values = [predictor.predict_profile({"C_TEST": c})["TestCompound"][0] for c in grid]
    assert values == sorted(values)


def test_suprathreshold_behaviour_is_unchanged_by_the_onset_ramp():
    """TIER 4.5: the ramp saturates at the band edge; nothing above it moves."""
    predictor = _mock_predictor_with_test_compound(threshold_ppb=1.0)
    for oav in [2.0, 4.0, 25.0, 100.0]:
        value = predictor.predict_profile({"C_TEST": oav})["TestCompound"][0]
        assert value == pytest.approx(pow(oav, 0.5))


def test_export_qda_profile_is_absolute_not_self_normalized():
    """TIER 4.5/4.6: the top note must not always be 10.0 on an 'absolute' 0-10 scale."""
    predictor = SensoryPredictor()

    weak = predictor.export_qda_profile(predictor.predict_profile({"2-Furfurylthiol (FFT)": 0.1}))
    assert weak["roasted"] < 10.0
    assert weak["roasted"] == pytest.approx(3.162, rel=1e-3)

    strong = predictor.export_qda_profile(predictor.predict_profile({"2-Furfurylthiol (FFT)": 10.0}))
    assert strong["roasted"] > weak["roasted"]

    # Still clipped at the top of the scale.
    saturating = predictor.export_qda_profile({"2-Furfurylthiol (FFT)": (1000.0, 0.0)})
    assert max(saturating.values()) == pytest.approx(10.0)


if __name__ == "__main__":
    pytest.main([__file__])
