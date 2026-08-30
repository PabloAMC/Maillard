import pytest

from src.extrusion import (
    autoclave_damage_increment,
    build_extrusion_process_profile,
    build_sequential_isothermal_zones,
    classify_volatile_compound,
    compute_extrusion_headspace_adjustment,
    compute_sme_temperature_offset,
    estimate_die_exit_temperature_celsius,
    estimate_mean_residence_time_seconds,
    estimate_relative_rtd_spread,
    normalize_moisture_regime,
    pre_extrusion_damage_load,
)


def test_compute_sme_temperature_offset_stays_inside_report_window():
    # Expectations updated 2026-08-27: the SME response was recalibrated from
    # min(40, max(5, slope*SME)) — which saturated at 250 (LME) / 400 (HME)
    # kJ/kg and jumped to 5 C for any SME > 0 — to a melt energy balance with a
    # tanh soft ceiling. There is no longer a 5 C floor or a 40 C cap.
    assert compute_sme_temperature_offset(0.0, "hme", 0.7) == 0.0
    assert 0.0 < compute_sme_temperature_offset(60.0, "hme", 0.7) < 5.0
    assert 15.0 < compute_sme_temperature_offset(250.0, "lme", 0.35) < 25.0


def test_sme_temperature_offset_discriminates_across_the_real_window():
    """300-800 kJ/kg is the real twin-screw window; it must not be saturated."""
    for regime, aw in (("hme", 0.7), ("lme", 0.35)):
        offsets = [compute_sme_temperature_offset(sme, regime, aw) for sme in (300.0, 450.0, 600.0, 800.0)]
        assert offsets == sorted(offsets)
        # Strictly increasing, with a real gradient (not a saturated plateau).
        for lower, higher in zip(offsets, offsets[1:]):
            assert higher - lower > 1.0
        assert offsets[-1] > 1.5 * offsets[0]


def test_normalize_moisture_regime_infers_from_aw():
    assert normalize_moisture_regime(None, 0.30) == "lme"
    assert normalize_moisture_regime(None, 0.85) == "hme"


def test_build_extrusion_process_profile_combines_damage_and_zones():
    profile = build_extrusion_process_profile(
        base_temperature_celsius=150.0,
        water_activity=0.55,
        protein_type="pea_iso",
        sme_kj_per_kg=150.0,
        moisture_regime="hme",
        sterilization_temperature_celsius=124.0,
        sterilization_time_minutes=25.0,
    )

    base_load = pre_extrusion_damage_load("pea_iso")
    assert profile["effective_temperature_celsius"] > 150.0
    assert len(profile["zone_profile"]) == 3
    assert profile["total_damage_load"]["furosine_mg_per_kg"] > base_load["furosine_mg_per_kg"]
    assert profile["total_damage_load"]["lal_mg_per_kg"] > base_load["lal_mg_per_kg"]
    assert profile["die_exit_temperature_celsius"] < profile["effective_temperature_celsius"]
    assert profile["rtd_assessment"]["decision"] == "sequential_zone_sufficient_for_current_use_case"
    assert "Hexanal" in profile["volatile_transport"]["panel"]


def test_extrusion_headspace_adjustment_penalizes_hexanal_more_than_pyrazine_at_the_die():
    profile = build_extrusion_process_profile(
        base_temperature_celsius=145.0,
        water_activity=0.75,
        protein_type="soy_iso",
        sme_kj_per_kg=180.0,
        moisture_regime="hme",
        screw_speed_rpm=280.0,
        feed_rate_kg_per_h=4.6,
    )

    hexanal = compute_extrusion_headspace_adjustment("Hexanal", profile)
    pyrazine = compute_extrusion_headspace_adjustment("2,5-Dimethylpyrazine", profile)

    assert hexanal["die_stripping_fraction"] > pyrazine["die_stripping_fraction"]
    assert hexanal["combined_headspace_factor"] < pyrazine["combined_headspace_factor"]


def test_explicit_die_exit_temperature_override_is_respected():
    assert estimate_die_exit_temperature_celsius(165.0, "hme", die_exit_temperature_celsius=60.0) == 60.0

# ---------------------------------------------------------------------------
# Audit remediation 2026-08-27 (item 3.1: SME desaturation and related fixes)
# ---------------------------------------------------------------------------

REAL_SME_WINDOW = (300.0, 450.0, 600.0, 800.0)


def _strictly_monotone(values, *, increasing: bool) -> bool:
    pairs = list(zip(values, values[1:]))
    if increasing:
        return all(higher > lower for lower, higher in pairs)
    return all(higher < lower for lower, higher in pairs)


def test_residence_time_discriminates_across_the_real_sme_window():
    values = [
        estimate_mean_residence_time_seconds(sme_kj_per_kg=sme, moisture_regime="hme")
        for sme in REAL_SME_WINDOW
    ]
    assert _strictly_monotone(values, increasing=False)
    assert (values[0] - values[-1]) / values[0] > 0.10


def test_rtd_spread_discriminates_across_the_real_sme_window():
    for regime in ("hme", "lme"):
        values = [
            estimate_relative_rtd_spread(sme_kj_per_kg=sme, moisture_regime=regime)
            for sme in REAL_SME_WINDOW
        ]
        assert _strictly_monotone(values, increasing=True)


def test_headspace_shear_release_discriminates_across_the_real_sme_window():
    factors = []
    for sme in REAL_SME_WINDOW:
        profile = build_extrusion_process_profile(
            base_temperature_celsius=145.0,
            water_activity=0.75,
            protein_type="soy_iso",
            sme_kj_per_kg=sme,
            moisture_regime="hme",
        )
        factors.append(compute_extrusion_headspace_adjustment("Hexanal", profile)["shear_release_factor"])
    assert _strictly_monotone(factors, increasing=True)


def test_die_stripping_responds_to_water_activity_above_140c():
    """The old min() capped the fraction, so aw had no effect at a hot die."""
    fractions = []
    for aw in (0.35, 0.55, 0.75, 0.85):
        profile = build_extrusion_process_profile(
            base_temperature_celsius=160.0,
            water_activity=aw,
            protein_type="soy_iso",
            sme_kj_per_kg=500.0,
            moisture_regime="hme",
        )
        assert profile["die_exit_temperature_celsius"] > 140.0
        fractions.append(compute_extrusion_headspace_adjustment("Hexanal", profile)["die_stripping_fraction"])
    assert _strictly_monotone(fractions, increasing=True)


def test_autoclave_damage_has_no_step_discontinuity_at_121c():
    def furosine(temp):
        return autoclave_damage_increment(
            temp, 25.0, protein_type="pea_iso", moisture_regime="hme"
        )["furosine_mg_per_kg"]

    below, at, above = furosine(120.99), furosine(121.0), furosine(121.01)
    assert below < at < above
    # A 0.02 C excursion must not move the marker by more than a few percent.
    assert (above - below) / at < 0.05
    # Deep-ambient conditions still contribute nothing.
    assert furosine(25.0) == 0.0


def test_autoclave_moisture_factor_applies_to_lal_as_well_as_furosine():
    wet = autoclave_damage_increment(130.0, 20.0, protein_type="soy_iso", moisture_regime="hme")
    dry = autoclave_damage_increment(130.0, 20.0, protein_type="soy_iso", moisture_regime="lme")
    assert wet["furosine_mg_per_kg"] > dry["furosine_mg_per_kg"]
    assert wet["lal_mg_per_kg"] > dry["lal_mg_per_kg"]


def test_mismatched_zone_lists_raise_instead_of_truncating():
    with pytest.raises(ValueError):
        build_sequential_isothermal_zones(
            120.0, 160.0,
            zone_temperatures=[120.0, 140.0, 150.0, 160.0, 170.0],
            zone_time_fractions=[0.3, 0.3, 0.4],
        )


def test_classify_volatile_handles_aldehyde_spellings():
    assert classify_volatile_compound("Benzaldehyde") == "aldehyde"
    assert classify_volatile_compound("Phenylacetaldehyde") == "aldehyde"
    assert classify_volatile_compound("Acetaldehyde") == "aldehyde"
    assert classify_volatile_compound("Hexanal") == "aldehyde"
    # Furanic carbonyls stay with the furan class.
    assert classify_volatile_compound("Furfural") == "furan"
    assert classify_volatile_compound("2-Pentylfuran") == "furan"
    assert classify_volatile_compound("2-furfurylthiol") == "sulfur"


def test_profile_is_inactive_without_thermal_exposure_and_adds_no_damage():
    ambient = build_extrusion_process_profile(
        base_temperature_celsius=25.0,
        water_activity=0.8,
        protein_type="pea_iso",
    )
    baseline = pre_extrusion_damage_load("pea_iso")

    assert ambient["active"] is False
    assert ambient["thermal_exposure"] is False
    assert ambient["thermal_exposure_basis"] == []
    assert ambient["process_damage_increment"] == {"furosine_mg_per_kg": 0.0, "lal_mg_per_kg": 0.0}
    # Ingredient provenance is preserved and clearly separated from processing.
    assert ambient["total_damage_load"] == baseline
    assert ambient["damage_load_attribution"]["ingredient_baseline"] == baseline


def test_profile_is_active_once_there_is_real_thermal_exposure():
    hot = build_extrusion_process_profile(
        base_temperature_celsius=145.0,
        water_activity=0.6,
        protein_type="soy_iso",
        sme_kj_per_kg=400.0,
    )
    assert hot["active"] is True
    assert "barrel_temperature_above_ambient_threshold" in hot["thermal_exposure_basis"]
    assert "mechanical_energy_input" in hot["thermal_exposure_basis"]

    sterilized_cold = build_extrusion_process_profile(
        base_temperature_celsius=25.0,
        water_activity=0.8,
        protein_type="pea_iso",
        sterilization_temperature_celsius=124.0,
        sterilization_time_minutes=25.0,
    )
    assert sterilized_cold["active"] is True
    assert sterilized_cold["process_damage_increment"]["furosine_mg_per_kg"] > 0.0


def test_zero_process_time_gates_the_profile_off():
    profile = build_extrusion_process_profile(
        base_temperature_celsius=145.0,
        water_activity=0.6,
        protein_type="soy_iso",
        sme_kj_per_kg=400.0,
        process_time_minutes=0.0,
    )
    assert profile["active"] is False
    assert profile["total_damage_load"] == pre_extrusion_damage_load("soy_iso")
