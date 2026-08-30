from pathlib import Path

import pytest
import yaml

from src import recommend as recommend_module
from src.literature_runtime import build_family_upstream_contract
from src.pathway_extractor import ElementaryStep, Species
from src.recommend import (
    DEFAULT_PROJECTION_STRATEGY,
    Recommender,
    _apply_output_projection,
    _canon,
    _estimate_projection_budget,
    _headspace_observability_factor,
    _mw_from_smiles,
    _project_weighted_flux_to_ppb,
    _is_budget_relevant_species,
    _is_observable_target_species,
    _resolve_upstream_observability_factor,
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
    # RETARGETED 2026-08-27 (cause: Henry's-law table rebuilt against Sander 5.0.0,
    # owner-approved 2026-08-26). Acrylamide's Kaw was raised 1.0e-9 -> 6.0e-8 on a
    # type-V Sander cluster (Mackay 2006d / Duchowicz 2020 / HSDB 2015, all landing
    # 5.9-7.3e-8), which puts it ABOVE the 1.0e-8 `_NON_OBSERVABLE_KAW_THRESHOLD`.
    # It is therefore no longer a valid example of a gated species -- that is the
    # intended effect of the correction, not a regression, and it is numerically
    # inert on the panel because acrylamide is produced by src/safety.py rather than
    # through the Henry-gated projection loop.
    # The CONTRACT under test is unchanged: genuinely non-volatile toxic targets must
    # not be treated as headspace-observable. CML is now the honest example (Kaw
    # 1.0e-12, zwitterionic at food pH), so the assertion moves to it and acrylamide
    # keeps an explicit assertion recording its new classification.
    cml_canon = _canon("NCCCCC(N)C(=O)O")
    acrylamide_canon = _canon("C=CC(=O)N")
    target_lookup = {
        cml_canon: {"name": "Nε-(Carboxymethyl)lysine (CML)", "type": "toxic", "data": {}},
        acrylamide_canon: {"name": "Acrylamide", "type": "toxic", "data": {}},
    }

    assert not _is_observable_target_species(
        Species("Nε-(Carboxymethyl)lysine (CML)", "NCCCCC(N)C(=O)O"), target_lookup
    )
    assert _is_observable_target_species(Species("Acrylamide", "C=CC(=O)N"), target_lookup), (
        "Acrylamide should now sit above the 1e-8 observability gate (Kaw 6.0e-8, Sander "
        "5.0.0). If this fails, the Henry table has been reverted -- check "
        "data/lit/henry_constants.yml before touching this test."
    )


def test_budget_projection_keeps_curated_observable_furan_targets():
    furfural_canon = _canon("O=Cc1ccco1")
    target_lookup = {
        furfural_canon: {"name": "Furfural", "type": "desirable", "data": {}},
    }

    assert _is_budget_relevant_species(Species("furfural", "O=Cc1ccco1"), target_lookup)


def test_hmf_moved_above_the_headspace_observability_gate():
    # RENAMED + INVERTED 2026-08-27 (cause: Henry's-law table rebuilt against Sander
    # 5.0.0, owner-approved 2026-08-26). HMF's Kaw was raised 1.0e-10 -> 5.0e-8. No
    # direct measurement of HMF's Henry constant exists anywhere, but every available
    # route lands above the gate (Sander/HSDB 2015 -> 2.2e-8; OPERA 2.6 -> 9.4e-8;
    # vapour-pressure/solubility calculation 1.7e-8 to 1.4e-6), so the shipped 1e-10
    # was 200-1000x low. Crossing 1.0e-8 flips HMF from `low_headspace` to
    # `observable`, which is the whole point of the correction.
    # The old assertion (`not observable`) is therefore now FALSE BY DESIGN, and it is
    # inverted rather than deleted so the classification stays pinned: if HMF ever
    # reads as low_headspace again, the Henry table has been reverted.
    # Consequence recorded for the panel: there are now ZERO low_headspace rows in the
    # benchmark panel.
    hmf_canon = _canon("O=Cc1ccc(CO)o1")
    target_lookup = {
        hmf_canon: {"name": "5-Hydroxymethylfurfural (HMF)", "type": "desirable", "data": {}},
    }

    assert _is_observable_target_species(Species("HMF", "O=Cc1ccc(CO)o1"), target_lookup)
    assert _is_budget_relevant_species(Species("HMF", "O=Cc1ccc(CO)o1"), target_lookup)


def test_budget_projection_excludes_inorganic_targets_from_budget():
    h2s_canon = _canon("S")
    target_lookup = {h2s_canon: {"name": "Hydrogen Sulfide", "type": "desirable", "data": {}}}

    assert not _is_budget_relevant_species(Species("H2S", "S"), target_lookup)


def test_disulfide_target_smiles_matches_generated_species():
    # 2026-08-27: pinned to the engine constant rather than a literal, after
    # the audit found the literal here, the template layer and the species
    # registry all carrying the WRONG regioisomer while the curated layer had
    # the right one. A hardcoded string cannot catch that class of drift.
    from src.smirks_engine import _FURYL_DISULFIDE_CANONICAL

    generated = _canon(_FURYL_DISULFIDE_CANONICAL)
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
    disulfide = Species("bis(2-methyl-3-furyl) disulfide", "Cc1c(SSc2ccoc2C)cco1")
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

    # Wave T4 2026-08-27: `steps` and `downstream_margin_kcal` dropped from the
    # signature — AST-verified unused in the body since Wave S1 replaced the
    # selection heuristic that read them.
    selected = _select_accumulating_projection_species(
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

    expected_ratio = _mw_from_smiles(furfural.smiles) / _mw_from_smiles(hmf.smiles)
    observed_ratio = projected[_canon(furfural.smiles)] / projected[_canon(hmf.smiles)]
    assert observed_ratio == pytest.approx(expected_ratio, rel=1.0e-3)


def test_low_headspace_factor_is_temperature_aware_but_stays_conservative(monkeypatch):
    # REWRITTEN 2026-08-27 (cause: Henry's-law table rebuilt against Sander 5.0.0,
    # owner-approved 2026-08-26). This test used HMF as its low-headspace exemplar.
    # HMF's Kaw moved 1.0e-10 -> 5.0e-8, i.e. above the 1.0e-8 gate, so the function
    # now short-circuits to 1.0 for it and the old `hmf_25c < hmf_150c < 0.05` band is
    # unreachable through that species.
    #
    # HONEST CONSEQUENCE, recorded here because it is easy to miss: after the refresh
    # the ONLY registry entries still below the gate are CML and CEL at Kaw 1.0e-12,
    # and both clamp to the function's 1.0e-6 intrinsic floor at every temperature.
    # So no SHIPPED compound currently exercises the temperature ramp at all. Rather
    # than delete the contract (the ramp is still live code and still governs any
    # future sub-gate compound), the exemplar is now a synthetic registry entry placed
    # in the sensitive band. The test is not relaxed -- it asserts the same three
    # properties as before: floor respected, monotone in temperature, stays conservative.
    furfural = Species("furfural", "O=Cc1ccco1")
    hmf = Species("HMF", "O=Cc1ccc(CO)o1")
    cml = Species("Nε-(Carboxymethyl)lysine (CML)", "NCCCCC(N)C(=O)O")

    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
        _canon(hmf.smiles): {"name": "5-Hydroxymethylfurfural (HMF)", "type": "toxic", "data": {}},
        _canon(cml.smiles): {"name": "Nε-(Carboxymethyl)lysine (CML)", "type": "toxic", "data": {}},
    }

    # Above-gate species are passed through untouched.
    assert _headspace_observability_factor(furfural, target_lookup, 423.15) == 1.0
    assert _headspace_observability_factor(hmf, target_lookup, 423.15) == 1.0

    # Deep sub-gate species clamp to the conservative floor at every temperature.
    cml_25c = _headspace_observability_factor(cml, target_lookup, 298.15)
    cml_150c = _headspace_observability_factor(cml, target_lookup, 423.15)
    assert cml_25c == pytest.approx(1.0e-6)
    assert cml_150c == pytest.approx(1.0e-6)

    # The temperature ramp itself: a synthetic entry inside the sensitive band
    # (below the 1e-8 gate, above the floor once scaled against furfural).
    probe = Species("probe-low-headspace", "O=Cc1ccc(CO)o1")
    probe_entry = {"name": "5-Hydroxymethylfurfural (HMF)", "Kaw_25c": 1.0e-10}
    monkeypatch.setattr(
        recommend_module,
        "_henry_entry_for_species",
        lambda species, lookup: probe_entry if species is probe else None,
    )
    probe_25c = _headspace_observability_factor(probe, target_lookup, 298.15)
    probe_150c = _headspace_observability_factor(probe, target_lookup, 423.15)
    assert 1.0e-6 <= probe_25c < probe_150c < 0.05


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

    # UPDATED 2026-08-27 (cause: Henry's-law table rebuilt against Sander 5.0.0,
    # owner-approved 2026-08-26). HMF's Kaw moved 1.0e-10 -> 5.0e-8, above the 1.0e-8
    # observability gate, so it is no longer suppressed relative to its proxy value and
    # the old `observable < 10.0` / `observable < proxy_ppb` assertions are false by
    # design. What this test actually exists to prove -- that the projection exposes a
    # proxy channel and an observable channel and that the metadata row agrees with the
    # returned mapping -- is unchanged and still asserted for both species.
    assert metadata[_canon(furfural.smiles)]["proxy_ppb"] == pytest.approx(100.0)
    assert metadata[_canon(furfural.smiles)]["observable_ppb"] == pytest.approx(observable[_canon(furfural.smiles)])
    assert observable[_canon(furfural.smiles)] == pytest.approx(100.0)
    assert metadata[_canon(hmf.smiles)]["proxy_ppb"] == pytest.approx(100.0)
    assert metadata[_canon(hmf.smiles)]["observable_ppb"] == pytest.approx(observable[_canon(hmf.smiles)])
    assert observable[_canon(hmf.smiles)] == pytest.approx(100.0)


def test_projection_budget_increases_with_temperature_time_and_load():
    low_temp = _estimate_projection_budget({"ribose": 100.0, "cysteine": 100.0}, 373.15, 5.0)
    high_temp = _estimate_projection_budget({"ribose": 100.0, "cysteine": 100.0}, 423.15, 5.0)
    long_time = _estimate_projection_budget({"ribose": 100.0, "cysteine": 100.0}, 423.15, 60.0)
    high_load = _estimate_projection_budget({"ribose": 400.0, "cysteine": 400.0}, 423.15, 60.0)

    assert high_temp.temperature_factor > low_temp.temperature_factor
    assert high_temp.total_volatile_budget_molar > low_temp.total_volatile_budget_molar
    assert long_time.time_factor > high_temp.time_factor
    assert long_time.total_volatile_budget_molar > high_temp.total_volatile_budget_molar
    assert high_load.limiting_precursor_molar > long_time.limiting_precursor_molar
    assert high_load.total_volatile_budget_molar > long_time.total_volatile_budget_molar


def test_output_projection_metadata_exposes_budget_context():
    furfural = Species("furfural", "O=Cc1ccco1")
    raw = {
        _canon(furfural.smiles): 25.0,
    }
    species_catalog = {
        _canon(furfural.smiles): furfural,
    }
    target_lookup = {
        _canon(furfural.smiles): {"name": "Furfural", "type": "desirable", "data": {}},
    }
    budget = _estimate_projection_budget({"ribose": 100.0, "cysteine": 100.0}, 423.15, 30.0)

    _observable, metadata = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=423.15,
        protein_type="free",
        projection_budget=budget,
    )

    projection = metadata[_canon(furfural.smiles)]
    assert projection["total_volatile_budget_molar"] == pytest.approx(budget.total_volatile_budget_molar)
    assert projection["projection_temperature_factor"] == pytest.approx(budget.temperature_factor)
    assert projection["projection_time_factor"] == pytest.approx(budget.time_factor)
    assert projection["projection_strategy_name"] == DEFAULT_PROJECTION_STRATEGY.name
    assert projection["projection_ppb_basis"] == DEFAULT_PROJECTION_STRATEGY.ppb_basis


def test_projection_budget_uses_explicit_strategy_constants():
    # REWRITTEN 2026-08-27 (cause: projection budget retune, audit remediation Part 1).
    # v1 computed `volatile_yield_fraction = baseline + severity_slope * severity`,
    # where severity was a sigmoid centred on 110 C with width 18 C. That sigmoid
    # saturates by construction (0.966 at 170 C, 0.988 at 190 C), so the budget grew
    # only ~1.11x from 150 C to 190 C against a real Arrhenius drive of ~20x -- the
    # cause of the spurious furfural maximum near 145-150 C.
    # v2 replaces the slope term with a first-order conversion extent under an apparent
    # Arrhenius rate; `severity_volatile_yield_slope` no longer exists on the strategy.
    # `severity` itself survives, unchanged, as the bounded process-state index consumed
    # by melanoidin trapping / depth bias / the direct-sulfur bonus, and is still
    # asserted below so nothing silently re-wires it into the budget.
    import math

    budget = _estimate_projection_budget({"ribose": 100.0, "cysteine": 100.0}, 423.15, 30.0)
    strategy = DEFAULT_PROJECTION_STRATEGY

    # severity is still the plain product of the two bounded factors.
    assert budget.severity == pytest.approx(budget.temperature_factor * budget.time_factor)

    # The budget scale now follows k(T) * t, referenced to T_ref with the literature
    # apparent activation energy carried on the strategy.
    gas_constant = 8.314462618
    expected_drive = (30.0 / strategy.reference_conversion_time_min) * math.exp(
        -(strategy.apparent_activation_energy_kj_mol * 1.0e3 / gas_constant)
        * (1.0 / 423.15 - 1.0 / strategy.reference_temperature_kelvin)
    )
    expected_extent = 1.0 - math.exp(-expected_drive)
    expected_yield = strategy.baseline_volatile_yield_fraction + (
        (strategy.conversion_ceiling_fraction - strategy.baseline_volatile_yield_fraction)
        * expected_extent
    )

    assert budget.kinetic_drive == pytest.approx(expected_drive)
    assert budget.conversion_extent == pytest.approx(expected_extent)
    assert budget.volatile_yield_fraction == pytest.approx(expected_yield)


def test_projection_budget_temperature_dependence_tracks_arrhenius_not_a_sigmoid():
    """Regression guard for the 2026-08-27 retune (audit remediation Part 1).

    The v1 budget saturated: 150 C -> 190 C grew it by only 1.108x. The whole point of
    v2 is that the budget inherits the network's exponential temperature dependence in
    the low-conversion regime the panel occupies. Pin that, so a future edit cannot
    quietly reintroduce a saturating scale.
    """
    import math

    strategy = DEFAULT_PROJECTION_STRATEGY
    precursors = {"ribose": 100.0, "cysteine": 100.0}
    budget_150 = _estimate_projection_budget(precursors, 423.15, 60.0)
    budget_170 = _estimate_projection_budget(precursors, 443.15, 60.0)
    budget_190 = _estimate_projection_budget(precursors, 463.15, 60.0)

    gas_constant = 8.314462618
    expected_150_to_170 = math.exp(
        -(strategy.apparent_activation_energy_kj_mol * 1.0e3 / gas_constant)
        * (1.0 / 443.15 - 1.0 / 423.15)
    )

    observed = budget_170.total_volatile_budget_molar / budget_150.total_volatile_budget_molar
    # ~4.6x, matching the Arrhenius factor to within the first-order rollover.
    assert observed == pytest.approx(expected_150_to_170, rel=0.02)
    assert observed > 4.0, (
        f"150 C -> 170 C budget growth is {observed:.3f}x. The v1 sigmoid gave 1.07x here; "
        "anything back in that range means the saturating parameterisation has returned."
    )
    assert budget_190.total_volatile_budget_molar > budget_170.total_volatile_budget_molar


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

    assert metadata[_canon(furfural.smiles)]["base_matrix_factor"] == pytest.approx(0.5)
    assert metadata[_canon(furfural.smiles)]["class_matrix_factor"] == pytest.approx(0.945)
    assert metadata[_canon(furfural.smiles)]["matrix_factor"] == pytest.approx(0.4725)
    assert observable[_canon(furfural.smiles)] == pytest.approx(47.25)


def test_output_projection_respects_denaturation_state_in_matrix_retention_fallback():
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

    observable_native, metadata_native = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=313.15,
        protein_type="pea_iso",
        denaturation_state=0.0,
        fat_fraction=0.0,
        protein_fraction=1.0,
    )
    observable_denatured, metadata_denatured = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=313.15,
        protein_type="pea_iso",
        denaturation_state=1.0,
        fat_fraction=0.0,
        protein_fraction=1.0,
    )

    assert metadata_native[_canon(furfural.smiles)]["base_matrix_factor"] == pytest.approx(0.42)
    assert metadata_denatured[_canon(furfural.smiles)]["base_matrix_factor"] == pytest.approx(0.58)
    assert metadata_native[_canon(furfural.smiles)]["class_matrix_factor"] == pytest.approx(0.92)
    assert metadata_denatured[_canon(furfural.smiles)]["class_matrix_factor"] == pytest.approx(0.97)
    assert metadata_native[_canon(furfural.smiles)]["matrix_factor"] == pytest.approx(0.3864)
    assert metadata_denatured[_canon(furfural.smiles)]["matrix_factor"] == pytest.approx(0.5626)
    assert observable_native[_canon(furfural.smiles)] < observable_denatured[_canon(furfural.smiles)]


def test_melanoidin_trapping_penalizes_fft_more_than_mft_in_heated_soy_matrix():
    fft = Species("2-furfurylthiol", "SCc1ccco1")
    mft = Species("2-methyl-3-furanthiol", "Cc1occc1S")
    raw = {
        _canon(fft.smiles): 100.0,
        _canon(mft.smiles): 100.0,
    }
    species_catalog = {
        _canon(fft.smiles): fft,
        _canon(mft.smiles): mft,
    }
    target_lookup = {
        _canon(fft.smiles): {"name": "2-Furfurylthiol (FFT)", "type": "desirable", "data": {}},
        _canon(mft.smiles): {"name": "2-Methyl-3-furanthiol (MFT)", "type": "desirable", "data": {}},
    }

    observable, metadata = _apply_output_projection(
        raw,
        species_catalog,
        target_lookup,
        temperature_kelvin=393.15,
        protein_type="soy_iso",
        time_minutes=20.0,
        denaturation_state=0.6,
        fat_fraction=0.0,
        protein_fraction=1.0,
    )

    fft_row = metadata[_canon(fft.smiles)]
    mft_row = metadata[_canon(mft.smiles)]

    assert fft_row["process_state"] == "heated_matrix"
    assert fft_row["melanoidin_trapping_factor"] < mft_row["melanoidin_trapping_factor"] < 1.0
    assert observable[_canon(fft.smiles)] < observable[_canon(mft.smiles)]


def test_melanoidin_trapping_can_penalize_free_thiamine_systems_when_family_16_is_active():
    mft = Species("2-methyl-3-furanthiol", "Cc1occc1S")
    canon = _canon(mft.smiles)
    contract = build_family_upstream_contract(
        sugars=["xylose"],
        amino_acids=["cysteine", "thiamine"],
        additives=[],
        protein_type="free",
        pH=6.0,
        process_state="heated_matrix",
        temperature_celsius=145.0,
        time_minutes=20.0,
        water_activity=0.98,
        molar_ratios={"xylose": 10.0, "cysteine": 10.0, "thiamine": 10.0},
        thiamine_availability={"available": True, "source": "benchmark_native_default"},
    )

    observable, metadata = _apply_output_projection(
        {canon: 100.0},
        {canon: mft},
        {canon: {"name": "2-Methyl-3-furanthiol (MFT)", "type": "desirable", "data": {}}},
        temperature_kelvin=418.15,
        protein_type="free",
        time_minutes=20.0,
        water_activity=0.98,
        family_upstream_contract=contract,
    )

    assert contract["family_lanes"]["16"]["active"] is True
    assert metadata[canon]["melanoidin_trapping_factor"] < 1.0
    assert observable[canon] < 100.0


def test_hydrolysate_supported_sulfur_observability_uses_peptide_release_uplift_but_preserves_source_ranking():
    contract = build_family_upstream_contract(
        sugars=["xylose"],
        amino_acids=["cysteine", "methionine", "glycine"],
        additives=[],
        support_cues=["soy protein isolate hydrolysate"],
        protein_type="free",
        pH=6.0,
        process_state="heated_matrix",
        temperature_celsius=120.0,
        time_minutes=30.0,
        water_activity=0.92,
        molar_ratios={"xylose": 50.0, "cysteine": 3.0, "methionine": 4.0, "glycine": 8.0},
    )

    expected_hydrolysate_uplift = 1.0 + (1.6 - 1.0) * 0.92

    soy_fft = _resolve_upstream_observability_factor(
        "2-Furfurylthiol",
        protein_source="soy_isolate",
        family_upstream_contract=contract,
    )
    wheat_fft = _resolve_upstream_observability_factor(
        "2-Furfurylthiol",
        protein_source="wheat_gluten",
        family_upstream_contract=contract,
    )
    soy_methional = _resolve_upstream_observability_factor(
        "Methional",
        protein_source="soy_isolate",
        family_upstream_contract=contract,
    )

    assert soy_fft == pytest.approx(0.13 * expected_hydrolysate_uplift * 1.25)
    assert wheat_fft == pytest.approx(soy_fft * (0.58 / 1.25))
    # RE-PINNED 2026-08-27 (Wave I) — REVERTED to the pre-Wave-H value, 0.05623 -> 0.0045.
    # Wave H had re-derived this factor 0.0045 -> 0.05623 against the two literature xylose
    # HVP benchmarks. The cold-start red team then established that BOTH of those
    # benchmarks are fabricated: the cited paper 10.1007/s10068-022-01194-w reacts protein
    # hydrolysates with glucose/fructose at pH 7.5 for 90 min and reports only RELATIVE
    # PEAK AREAS; it never mentions FFT or MFT and gives no absolute concentration for any
    # analyte, so the conc_ppb values that fit was solved against have no possible source.
    # Both files are now quarantined, the re-derivation record is RETRACTED, and the
    # constant is back to its pre-Wave-H value.
    # See: data/benchmarks/quarantined/README.md (Wave I section);
    #      results/validation/hydrolysate_observability_rederivation.md (RETRACTED banner);
    #      src/recommend.py::_HYDROLYSATE_SULFUR_OBSERVABILITY_PROFILES.
    # The property this test is named for -- that Methional is NOT source-sensitive while
    # the thiols are -- is unchanged, and is what the assertion below still checks.
    assert soy_methional == pytest.approx(0.0045)
    assert wheat_fft < soy_fft


def test_output_projection_surfaces_extrusion_process_state_and_surrogates():
    species = Species("Hexanal", "CCCCCC=O")
    canon = _canon(species.smiles)

    hydrated_observable, hydrated_meta = _apply_output_projection(
        {canon: 100.0},
        {canon: species},
        {},
        150.0 + 273.15,
        protein_type="soy_iso",
        time_minutes=3.0,
        water_activity=0.60,
        denaturation_state=0.8,
    )
    dry_observable, dry_meta = _apply_output_projection(
        {canon: 100.0},
        {canon: species},
        {},
        165.0 + 273.15,
        protein_type="soy_iso",
        time_minutes=3.0,
        water_activity=0.35,
        denaturation_state=0.95,
    )

    assert hydrated_meta[canon]["process_state"] == "aqueous_pre_extrusion_model"
    assert dry_meta[canon]["process_state"] == "extrusion_structured"
    assert dry_meta[canon]["extrusion_moisture_factor"] < hydrated_meta[canon]["extrusion_moisture_factor"]
    assert dry_meta[canon]["extrusion_structure_factor"] < hydrated_meta[canon]["extrusion_structure_factor"]
    assert dry_observable[canon] < hydrated_observable[canon]