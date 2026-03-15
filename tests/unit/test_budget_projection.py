from pathlib import Path

import pytest
import yaml

from src.pathway_extractor import ElementaryStep, Species
from src.recommend import (
    Recommender,
    _apply_output_projection,
    _canon,
    _headspace_observability_factor,
    _project_weighted_flux_to_ppb,
    _is_budget_relevant_species,
    _is_observable_target_species,
    _select_accumulating_projection_species,
)


ROOT = Path(__file__).resolve().parents[2]


def test_budget_projection_filters_bookkeeping_coproducts():
    target_lookup = {}

    assert not _is_budget_relevant_species(Species("water", "O"), target_lookup)
    assert not _is_budget_relevant_species(Species("CO2", "O=C=O"), target_lookup)
    assert not _is_budget_relevant_species(Species("elemental-sulfur", "[S]"), target_lookup)
    assert not _is_budget_relevant_species(Species("formaldehyde", "C=O"), target_lookup)
    assert _is_budget_relevant_species(Species("furfural", "O=Cc1ccco1"), target_lookup)


def test_budget_projection_filters_reactive_intermediates():
    target_lookup = {}

    assert not _is_budget_relevant_species(Species("dehydroalanine", "C=C(N)C(=O)O"), target_lookup)
    assert not _is_budget_relevant_species(Species("imine-adduct", "CC(=O)CN=CCS"), target_lookup)
    assert _is_budget_relevant_species(Species("2-furfurylthiol", "SCc1ccco1"), target_lookup)


def test_budget_projection_keeps_curated_targets_even_if_small():
    dmds_canon = _canon("CSSC")
    target_lookup = {dmds_canon: {"name": "Dimethyl disulfide", "type": "desirable", "data": {}}}

    assert _is_budget_relevant_species(Species("DMDS", "CSSC"), target_lookup)


def test_observability_marks_nonvolatile_toxic_targets_as_non_headspace_observable():
    acrylamide_canon = _canon("C=CC(=O)N")
    target_lookup = {
        acrylamide_canon: {"name": "Acrylamide", "type": "toxic", "data": {}},
    }

    assert not _is_observable_target_species(Species("Acrylamide", "C=CC(=O)N"), target_lookup)


def test_budget_projection_keeps_curated_observable_furan_targets():
    furfural_canon = _canon("O=Cc1ccco1")
    target_lookup = {
        furfural_canon: {"name": "Furfural", "type": "desirable", "data": {}},
    }

    assert _is_budget_relevant_species(Species("furfural", "O=Cc1ccco1"), target_lookup)


def test_observability_marks_hmf_as_low_headspace_even_if_budget_keeps_it_for_now():
    hmf_canon = _canon("O=Cc1ccc(CO)o1")
    target_lookup = {
        hmf_canon: {"name": "5-Hydroxymethylfurfural (HMF)", "type": "desirable", "data": {}},
    }

    assert not _is_observable_target_species(Species("HMF", "O=Cc1ccc(CO)o1"), target_lookup)
    assert _is_budget_relevant_species(Species("HMF", "O=Cc1ccc(CO)o1"), target_lookup)


def test_budget_projection_excludes_inorganic_targets_from_budget():
    h2s_canon = _canon("S")
    target_lookup = {h2s_canon: {"name": "Hydrogen Sulfide", "type": "desirable", "data": {}}}

    assert not _is_budget_relevant_species(Species("H2S", "S"), target_lookup)


def test_disulfide_target_smiles_matches_generated_species():
    generated = _canon("Cc1cc(SSc2ccoc2C)co1")
    with open(ROOT / "data" / "species" / "desirable_targets.yml", "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    disulfide = next(item for item in data["compounds"] if item["name"] == "Bis(2-methyl-3-furyl) disulfide")
    curated = _canon(disulfide["smiles"])

    assert curated == generated


def test_budget_projection_uses_curated_targets_only():
    recommender = Recommender()
    desirable = recommender._load_desirable()
    competing = recommender._load_off_flavours()
    curated = {
        _canon(item["smiles"])
        for item in [*desirable.values(), *competing.values()]
        if item.get("smiles")
    }

    assert _canon("O=Cc1ccco1") in curated
    assert _canon("Cc1occc1S") in curated
    assert _canon("CC(=O)C=O") not in curated


def test_projection_prefers_terminal_budget_relevant_endpoints():
    ribose = Species("D-Ribose", "O=CC(O)C(O)C(O)CO")
    furfural = Species("furfural", "O=Cc1ccco1")
    mft = Species("2-methyl-3-furanthiol", "Cc1occc1S")
    disulfide = Species("bis(2-methyl-3-furyl) disulfide", "Cc1cc(SSc2ccoc2C)co1")
    thiohemiacetal = Species("furfural-thiohemiacetal", "OC(S)c1ccco1")

    steps = [
        ElementaryStep(reactants=[ribose], products=[furfural], reaction_family="Enolisation_1_2"),
        ElementaryStep(reactants=[furfural], products=[thiohemiacetal], reaction_family="Thiohemiacetal_Formation"),
        ElementaryStep(reactants=[furfural], products=[mft], reaction_family="Thiol_Addition"),
        ElementaryStep(reactants=[mft], products=[disulfide], reaction_family="Thiol_Oxidation"),
    ]

    tracked_species = {
        _canon(ribose.smiles): (0.0, 10.0, 0, 10.0, 0.0),
        _canon(furfural.smiles): (24.0, 10.0, 1, 0.8, 0.0),
        _canon(mft.smiles): (24.8, 10.0, 2, 0.4, 0.0),
        _canon(disulfide.smiles): (27.4, 10.0, 3, 0.2, 0.0),
        _canon(thiohemiacetal.smiles): (25.1, 10.0, 2, 0.1, 0.0),
    }
    species_catalog = {
        _canon(ribose.smiles): ribose,
        _canon(furfural.smiles): furfural,
        _canon(mft.smiles): mft,
        _canon(disulfide.smiles): disulfide,
        _canon(thiohemiacetal.smiles): thiohemiacetal,
    }
    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
        _canon(mft.smiles): {"name": "2-Methyl-3-furanthiol", "type": "desirable", "data": {}},
        _canon(disulfide.smiles): {"name": "Bis(2-methyl-3-furyl) disulfide", "type": "desirable", "data": {}},
    }

    selected = _select_accumulating_projection_species(
        steps,
        tracked_species,
        species_catalog,
        target_lookup,
        exogenous_reactants={_canon(ribose.smiles)},
    )

    assert _canon(ribose.smiles) not in selected
    assert _canon(furfural.smiles) in selected
    assert _canon(mft.smiles) in selected
    assert _canon(disulfide.smiles) in selected
    assert _canon(thiohemiacetal.smiles) not in selected


def test_projection_downweights_low_headspace_targets_in_budget_allocation():
    precursor = Species("glucose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    furfural = Species("furfural", "O=Cc1ccco1")
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")

    tracked_species = {
        _canon(precursor.smiles): (0.0, 10.0, 0, 10.0, 0.0),
        _canon(furfural.smiles): (25.0, 10.0, 2, 0.5, 0.0),
        _canon(hmf.smiles): (25.0, 10.0, 2, 0.5, 0.0),
    }
    species_catalog = {
        _canon(precursor.smiles): precursor,
        _canon(furfural.smiles): furfural,
        _canon(hmf.smiles): hmf,
    }
    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "toxic", "data": {}},
    }

    projected = _project_weighted_flux_to_ppb(
        steps=[],
        tracked_species=tracked_species,
        best_paths={},
        species_catalog=species_catalog,
        corrected_initial={_canon(precursor.smiles): 10.0},
        target_lookup=target_lookup,
        exogenous_reactants={_canon(precursor.smiles)},
        temperature_kelvin=423.15,
        time_minutes=30.0,
    )

    assert projected[_canon(furfural.smiles)] > projected[_canon(hmf.smiles)] * 50.0


def test_low_headspace_factor_is_temperature_aware_but_stays_conservative():
    furfural = Species("furfural", "O=Cc1ccco1")
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")

    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "toxic", "data": {}},
    }

    hmf_25c = _headspace_observability_factor(hmf, target_lookup, 298.15)
    hmf_150c = _headspace_observability_factor(hmf, target_lookup, 423.15)
    furfural_150c = _headspace_observability_factor(furfural, target_lookup, 423.15)

    assert 1.0e-6 <= hmf_25c < hmf_150c < 0.05
    assert furfural_150c == 1.0


def test_output_projection_exposes_proxy_and_observable_channels():
    furfural = Species("furfural", "O=Cc1ccco1")
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")
    raw = {
        _canon(furfural.smiles): 100.0,
        _canon(hmf.smiles): 100.0,
    }
    species_catalog = {
        _canon(furfural.smiles): furfural,
        _canon(hmf.smiles): hmf,
    }
    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "toxic", "data": {}},
    }

    observable, metadata = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=423.15,
        protein_type="free",
        fat_fraction=0.0,
        protein_fraction=1.0,
    )

    assert metadata[_canon(furfural.smiles)]["proxy_ppb"] == pytest.approx(100.0)
    assert metadata[_canon(furfural.smiles)]["observable_ppb"] == pytest.approx(observable[_canon(furfural.smiles)])
    assert observable[_canon(furfural.smiles)] == pytest.approx(100.0)
    assert observable[_canon(hmf.smiles)] < 10.0


def test_output_projection_uses_matrix_retention_fallback_when_fractions_are_unspecified():
    furfural = Species("furfural", "O=Cc1ccco1")
    raw = {
        _canon(furfural.smiles): 100.0,
    }
    species_catalog = {
        _canon(furfural.smiles): furfural,
    }
    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
    }

    observable, metadata = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=313.15,
        protein_type="pea_iso",
        fat_fraction=0.0,
        protein_fraction=1.0,
    )

    assert metadata[_canon(furfural.smiles)]["matrix_factor"] == pytest.approx(0.5)
    assert observable[_canon(furfural.smiles)] == pytest.approx(50.0)