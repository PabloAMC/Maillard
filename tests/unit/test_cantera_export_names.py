from src.cantera_export import CanteraExporter, sanitize_cantera_name


def test_sanitize_cantera_name_replaces_parser_unsafe_characters():
    assert sanitize_cantera_name("bis(2-methyl-3-furyl) disulfide") == "bis_2_methyl_3_furyl_disulfide"
    assert sanitize_cantera_name("A/B + C, D") == "A_B_C_D"


def test_add_species_disambiguates_name_collisions_after_sanitization():
    exporter = CanteraExporter()

    first = exporter.add_species("CC=O", name="alpha-beta")
    second = exporter.add_species("C=CO", name="alpha beta")

    assert first == "alpha_beta"
    assert second == "alpha_beta_2"


def test_add_species_preserves_existing_name_for_same_species():
    exporter = CanteraExporter()

    first = exporter.add_species("SCc1ccco1", name="2-furfurylthiol")
    second = exporter.add_species("SCc1ccco1", name="2 furfurylthiol")

    assert first == "2_furfurylthiol"
    assert second == first