from src.extrusion import (
    build_extrusion_process_profile,
    compute_extrusion_headspace_adjustment,
    compute_sme_temperature_offset,
    estimate_die_exit_temperature_celsius,
    normalize_moisture_regime,
    pre_extrusion_damage_load,
)


def test_compute_sme_temperature_offset_stays_inside_report_window():
    assert compute_sme_temperature_offset(0.0, "hme", 0.7) == 0.0
    assert 5.0 <= compute_sme_temperature_offset(60.0, "hme", 0.7) <= 40.0
    assert compute_sme_temperature_offset(250.0, "lme", 0.35) == 40.0


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