"""Wave X (2026-08-28) — the step-level Hofmann tables, and the norfuraneol resolution.

These tests pin three things that Wave X decided and that a later wave could undo by
accident:

1. **The isotope constraint survives the step's return.** Wave N removed
   `norfuraneol + H2S -> MFT` because Cerny & Davidek 2003 showed by spiking that
   norfuraneol is unimportant as an IN-SITU intermediate. Wave X put the step back as a
   SLOW PARALLEL channel, because Hofmann & Schieberle 1998 Table 4 measures it directly
   when norfuraneol is FED. The two are compatible only while the channel stays a
   MINORITY of in-situ MFT flux, so that condition is asserted here rather than left to
   the reader's memory. If a future barrier change lets the retired route take the lane
   back, `test_norfuraneol_route_stays_a_minor_share_of_ribose_mft_flux` fails.

2. **The retired constant is not quietly reused.** 26.85 was fitted through a
   contradicted route and against values that are not in the paper. The new key must
   never be assigned it.

3. **Executability is measured, not assumed.** Six step-level rows are on the panel
   because the network can run them; nine are in `data/benchmarks/step_level_unreachable/`
   because it cannot. Both halves are asserted, so a row cannot drift from one to the
   other without a test noticing.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from src import smirks_engine
from src.barrier_constants import FAST_BARRIERS
from src.benchmark_validation import evaluate_benchmark, get_benchmark_files
from src.conditions import ReactionConditions
from src.precursor_resolver import resolve, resolve_many
from src.smirks_engine import (
    SmirksEngine,
    _MERCAPTO_2_PROPANONE,
    _MFT_CANONICAL,
    _NORFURANEOL_CANONICAL,
)

ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = ROOT / "data" / "benchmarks"
UNREACHABLE_DIR = BENCH_DIR / "step_level_unreachable"

MFT_COMPOUND = "2-Methyl-3-furanthiol (MFT)"
RIBOSE_BENCHMARK = BENCH_DIR / "hofmann1998_ribose_cysteine_145C_20min_pH5.json"

#: Cerny & Davidek 2003 (10.1021/jf026123f), verbatim (quoted in full in
#: docs/validation/isotope_topology_evidence.md): "The resulting 2-methyl-3-furanthiol was
#: mainly 13C5-labeled, suggesting that it stems from ribose and that
#: 4-hydroxy-5-methyl-3(2H)-furanone is unimportant as an intermediate."
#:
#: "Mainly" supports exactly one quantitative reading -- the non-norfuraneol fraction is
#: the MAJORITY -- so the norfuraneol share is below one half. The paper prints no
#: percentage and a tighter ceiling would be invented precision. The same constant is
#: decision rule 5 of
#: scripts/generators/fit_furanone_reductive_sulfhydrylation_hofmann.py, and the two must
#: agree; `test_the_ceiling_here_matches_the_one_the_fit_script_enforced` checks that.
ISOTOPE_SHARE_CEILING = 0.50

#: The six step-level rows Wave X put ON the panel, with the table each comes from.
STEP_LEVEL_PANEL_ROWS = {
    "hofmann1998_furan2aldehyde_h2s_145C_20min_pH5": "Table 3, furan-2-aldehyde row",
    "hofmann1998_norfuraneol_h2s_145C_20min_pH5": "Table 4, NF row (FIT TARGET)",
    "hofmann1998_norfuraneol_cysteine_145C_20min_pH5": "Table 10, NF/cysteine row",
    "hofmann1998_c2c3_recombination_145C_20min_pH3": "Table 8, experiment 1",
    "hofmann1998_c2c3_recombination_145C_20min_pH5": "Table 8, experiment 2 = Table 10 row 1",
    "hofmann1998_c2c3_recombination_145C_20min_pH7": "Table 8, experiment 3",
}

FIT_TARGET_ID = "hofmann1998_norfuraneol_h2s_145C_20min_pH5"

#: name -> (expected molecular formula, the engine constant this must equal after
#: canonicalisation). A precursor that canonicalises to a different string from the
#: engine's own copy is a SILENT DUPLICATE: the pool would carry two molecules where
#: chemistry has one, every template matches on canonical SMILES, and the fed copy would
#: be inert while looking perfectly correct in the YAML.
NEW_SPECIES = {
    "Hydroxyacetaldehyde": ("C2H4O2", "O=CCO"),
    "2-Oxopropanal": ("C3H4O2", "CC(=O)C=O"),
    "Mercapto-2-propanone": ("C3H6OS", _MERCAPTO_2_PROPANONE),
    "Norfuraneol": ("C5H6O3", _NORFURANEOL_CANONICAL),
    "Furan-2-aldehyde": ("C5H4O2", "O=Cc1ccco1"),
    "3-Deoxyribosulose": ("C5H8O4", "O=CC(=O)CC(O)CO"),
    "Hydrogen sulfide": ("H2S", "S"),
}


def _canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def _atom_counts(species_list):
    counts: dict[str, int] = {}
    charge = 0
    for species in species_list:
        mol = Chem.AddHs(Chem.MolFromSmiles(species.smiles))
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
        charge += Chem.GetFormalCharge(mol)
    return counts, charge


def _enumerate(names, ph=5.0):
    conditions = ReactionConditions(temperature_celsius=145.0, pH=ph, water_activity=0.98)
    return SmirksEngine(conditions).enumerate(resolve_many(names), max_generations=4)


def _ribose_mft_ppb(*, with_nf_channel: bool) -> float:
    """Predicted MFT on the ribose/cysteine row with the norfuraneol channel on or off.

    The channel is disabled by replacing the template function, NOT by raising its
    barrier: a barrier can be made large but never absent, and "large" is the very
    quantity this measurement is trying to bound.
    """
    saved = smirks_engine._norfuraneol_mft_parallel_route
    try:
        if not with_nf_channel:
            smirks_engine._norfuraneol_mft_parallel_route = lambda pool: []
        evaluation = evaluate_benchmark(RIBOSE_BENCHMARK)
        for comparison in evaluation.comparisons:
            if comparison.compound == MFT_COMPOUND:
                return float(comparison.predicted_ppb)
    finally:
        smirks_engine._norfuraneol_mft_parallel_route = saved
    raise AssertionError(f"{MFT_COMPOUND!r} is not among the ribose row's comparisons")


# ──────────────────────────────────────────────────────────────────────────────
# 1. The isotope constraint, as a test rather than as an absence
# ──────────────────────────────────────────────────────────────────────────────
def test_norfuraneol_route_stays_a_minor_share_of_ribose_mft_flux():
    """The re-added norfuraneol channel must stay a MINORITY of in-situ MFT flux.

    This is Wave N's isotope evidence, kept as a structural property instead of as a
    deleted step. The share is measured as

        share = 1 - MFT(channel off) / MFT(channel on)

    on the RIBOSE/CYSTEINE row, where norfuraneol is made in situ. That subtraction is a
    share ONLY because the Wave S1 flux propagator is ADDITIVE -- under a max-channel or
    multiplicative rule it would mean nothing, and re-adding the step would not have been
    defensible in the first place.
    """
    without = _ribose_mft_ppb(with_nf_channel=False)
    with_channel = _ribose_mft_ppb(with_nf_channel=True)

    assert with_channel > without > 0.0, (
        "the norfuraneol channel must ADD flux (additive propagator) and the incumbent "
        f"routes must still carry some: off={without}, on={with_channel}"
    )
    share = 1.0 - (without / with_channel)
    assert share < ISOTOPE_SHARE_CEILING, (
        f"the norfuraneol channel now supplies {share:.1%} of predicted MFT from "
        f"ribose/cysteine, a MAJORITY. Cerny & Davidek 2003 spiked authentic norfuraneol "
        f"into a [13C5]ribose/cysteine system and found the MFT 'mainly 13C5-labeled', i.e. "
        f"NOT from norfuraneol. Wave X re-added this step only as a SLOW PARALLEL channel; "
        f"if it has become the dominant one, either a barrier moved or the propagator "
        f"changed, and the Wave X resolution is falsified rather than merely regressed. "
        f"Read tasks/audit_remediation.md '## Wave X' (a) before changing this threshold."
    )


def test_the_ceiling_here_matches_the_one_the_fit_script_enforced():
    """The gate that rejected the fit and the gate that guards the shipped tree agree.

    A ceiling enforced in two places with two values is a ceiling enforced nowhere.
    """
    record = json.loads(
        (ROOT / "results" / "validation"
         / "furanone_reductive_sulfhydrylation_refit_hofmann.json").read_text(encoding="utf-8")
    )
    assert record["isotope_gate"]["ceiling"] == ISOTOPE_SHARE_CEILING


def test_the_table4_fit_was_rejected_and_the_constant_kept_its_unfitted_seed():
    """The shipped barrier is the class seed, because the fit failed the isotope gate.

    This pins the WAVE'S RESULT, not a number: Hofmann Table 4 and Cerny & Davidek's
    labelling experiment admit disjoint barrier ranges, so the model cannot reproduce
    both. If a later wave makes them compatible, this test should fail and be rewritten
    with the new reason -- it must not be deleted quietly.
    """
    record = json.loads(
        (ROOT / "results" / "validation"
         / "furanone_reductive_sulfhydrylation_refit_hofmann.json").read_text(encoding="utf-8")
    )
    assert record["isotope_gate"]["verdict"].startswith("VIOLATED"), record["isotope_gate"]["verdict"]
    assert record["adopted"] == record["incumbent"], "a rejected fit must not move the constant"
    assert float(FAST_BARRIERS["furanone_reductive_sulfhydrylation"][0]) == pytest.approx(
        record["adopted"]
    ), "the shipped constant disagrees with the fit record that decided it"
    # Disjointness, restated as an assertion so the finding cannot rot into prose.
    limiting = record["isotope_gate"]["isotope_limiting_barrier"]
    assert limiting is not None and limiting > record["selected_by_rules_1_to_4"], (
        "the two constraints are supposed to be INCOMPATIBLE here; if they now overlap, "
        "the fit should be re-run and this test rewritten"
    )


# ──────────────────────────────────────────────────────────────────────────────
# 2. The retired constant stays retired
# ──────────────────────────────────────────────────────────────────────────────
def test_retired_norfuraneol_barrier_is_neither_reused_nor_resurrected():
    retired = FAST_BARRIERS["thiol_addition_norfuraneol"]
    assert float(retired[0]) == pytest.approx(26.85), (
        "the retired constant is a PROVENANCE RECORD; changing its value destroys the record"
    )
    assert "RETIRED" in retired[1]

    new = FAST_BARRIERS["furanone_reductive_sulfhydrylation"]
    assert float(new[0]) != pytest.approx(26.85), (
        "the Wave X norfuraneol step must NOT inherit the retired 26.85. That value was "
        "fitted THROUGH a route the isotope literature contradicts AND against MFT 342 / "
        "FFT 200 ppb, which Wave S2b traced to this repository's own arithmetic and Wave W "
        "confirmed appear nowhere in 10.1021/jf9705983."
    )
    for token in ("26.85", "thiol_addition_norfuraneol"):
        assert token in new[1], (
            "the new constant's rationale must cross-reference the retired one, so the two "
            "can never be confused by a future reader"
        )


def test_the_retired_family_is_still_emitted_by_nothing():
    for names in (["D-Ribose", "L-Cysteine"], ["Norfuraneol", "Hydrogen sulfide"]):
        families = {step.reaction_family for step in _enumerate(names)}
        assert "Thiol_Addition_Norfuraneol" not in families, (
            f"[{'+'.join(names)}] the RETIRED family name is being emitted again; Wave X's "
            "step is `Furanone_Reductive_Sulfhydrylation` and carries its own barrier"
        )


# ──────────────────────────────────────────────────────────────────────────────
# 3. The new species and the new steps
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("name", sorted(NEW_SPECIES))
def test_new_precursor_species_are_verified_and_are_not_silent_duplicates(name):
    expected_formula, engine_copy = NEW_SPECIES[name]
    species = resolve(name)
    mol = Chem.MolFromSmiles(species.smiles)
    assert mol is not None, f"{name} does not parse"
    assert rdMolDescriptors.CalcMolFormula(mol).rstrip("+-") == expected_formula
    assert _canonical(species.smiles) == _canonical(engine_copy), (
        f"{name} canonicalises to {_canonical(species.smiles)!r} but the engine's own copy "
        f"is {_canonical(engine_copy)!r}. Every template matches on canonical SMILES, so a "
        "mismatch makes the FED molecule inert while the internally-made one reacts -- two "
        "copies of one molecule with different reactivity."
    )


def test_norfuraneol_mft_step_exists_and_is_exactly_balanced():
    steps = [
        step for step in _enumerate(["Norfuraneol", "Hydrogen sulfide"])
        if step.reaction_family == "Furanone_Reductive_Sulfhydrylation"
    ]
    assert len(steps) == 1, "expected exactly one norfuraneol -> MFT step"
    step = steps[0]
    reactant_counts, reactant_charge = _atom_counts(step.reactants)
    product_counts, product_charge = _atom_counts(step.products)
    assert reactant_counts == product_counts, (
        f"norfuraneol + H2S + 2[H] -> MFT + 2 H2O is unbalanced: "
        f"{reactant_counts} vs {product_counts}"
    )
    assert reactant_charge == product_charge
    assert reactant_counts == {"C": 5, "H": 10, "O": 3, "S": 1}
    assert _canonical(_MFT_CANONICAL) in {_canonical(p.smiles) for p in step.products}


def test_h2s_reducing_equivalent_couple_is_balanced_and_h2s_gated():
    """`2 H2S -> HSSH + 2[H]`, the second donor of the model's reducing-equivalent token.

    Sourced to Hofmann & Schieberle 1998's own sentence: the intermediate "might be
    reduced by a reductone or BY A FURTHER MOLECULE OF HYDROGEN SULFIDE". Without it every
    reductive sulfur step is a hidden downstream dependent of cysteine being present, and
    the paper's cysteine-free Tables 3, 4 and 7 are structurally zero.
    """
    steps = [
        step for step in _enumerate(["Norfuraneol", "Hydrogen sulfide"])
        if step.reaction_family == "Thiol_Oxidation"
        and all(r.smiles == "S" for r in step.reactants)
    ]
    assert len(steps) == 1
    reactant_counts, _ = _atom_counts(steps[0].reactants)
    product_counts, _ = _atom_counts(steps[0].products)
    assert reactant_counts == product_counts == {"H": 4, "S": 2}
    assert not any(species.smiles == "[S]" for species in steps[0].products), (
        "the retired elemental-sulfur balance token must not come back"
    )

    # It must not fire where there is no H2S at all.
    families = [
        step for step in _enumerate(["Hydroxyacetaldehyde", "Mercapto-2-propanone"])
        if step.reaction_family == "Thiol_Oxidation"
        and all(r.smiles == "S" for r in step.reactants)
    ]
    assert not families


def test_c2_c3_lane_no_longer_needs_free_h2s_to_close_the_ring():
    """Steps 2 and 3 of the C2+C3 lane contain no sulfur reagent; only step 1 needs H2S.

    Hofmann Table 8 / Table 10 row 1 is hydroxyacetaldehyde + mercapto-2-propanone with no
    H2S and no cysteine, and it is the single most effective MFT system in the paper. The
    lane Wave P built FROM that measurement used to abort before it on an H2S guard.
    """
    families = {step.reaction_family for step in _enumerate(["Hydroxyacetaldehyde", "Mercapto-2-propanone"])}
    assert "Mercaptoketone_Aldol_Addition" in families
    assert "Mercaptoketone_Cyclodehydration" in families
    assert "Mercaptoketone_Formation" not in families, (
        "step 1 needs H2S and there is none in this system; it must NOT fire"
    )


def test_thiol_addition_selects_its_substrate_by_structure_not_by_spelling():
    """A fed furan-2-aldehyde must react exactly like an internally-made 'furfural'."""
    families = {step.reaction_family for step in _enumerate(["Furan-2-aldehyde", "Hydrogen sulfide"])}
    assert {"Thiol_Addition_H2", "Thiohemiacetal_Formation", "Thiol_Dehydration"} <= families, (
        "the FFT lane did not fire on a precursor resolved under the paper's own name "
        "('furan-2-aldehyde'). If this fails, the substrate match has regressed to a LABEL "
        "lookup and the network's best-measured FFT route is spelling-dependent."
    )


# ──────────────────────────────────────────────────────────────────────────────
# 4. The ingestion: what is scored, what is not, and why
# ──────────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("benchmark_id", sorted(STEP_LEVEL_PANEL_ROWS))
def test_step_level_panel_rows_are_executable_and_produce_a_matched_prediction(benchmark_id):
    path = BENCH_DIR / f"{benchmark_id}.json"
    assert path.exists(), f"{benchmark_id} is missing from the panel"
    evaluation = evaluate_benchmark(path)
    assert evaluation.supported, evaluation.reason
    matched = [
        comparison for comparison in evaluation.comparisons
        if comparison.matched_name is not None and comparison.predicted_ppb > 0.0
    ]
    assert matched, (
        f"{benchmark_id} is on the SCORED panel but the model predicts nothing for it. A row "
        "the network cannot execute belongs in data/benchmarks/step_level_unreachable/, "
        "where it is named rather than silently skipped by the log-error aggregate."
    )


def test_the_fit_target_is_declared_in_the_file_and_in_the_fit_record():
    payload = json.loads((BENCH_DIR / f"{FIT_TARGET_ID}.json").read_text(encoding="utf-8"))
    assert "fit_target_declaration" in payload
    record = json.loads(
        (ROOT / "results" / "validation"
         / "furanone_reductive_sulfhydrylation_refit_hofmann.json").read_text(encoding="utf-8")
    )
    assert FIT_TARGET_ID in record["fit_target"], (
        "the fit record must NAME its target, because that string is how "
        "src/fit_target_index.py and scripts/ci/fit_target_gate.py remove the row from "
        "scored evidence"
    )
    assert record["parameters_fitted"] == 1


def test_no_step_level_row_is_ingested_twice():
    """Table 4's NF row and Table 10's NF/H2S row are one measurement; so are Table 8
    experiment 2 and Table 10's top row. Each must appear on the panel once."""
    seen: dict[float, list[str]] = {}
    for benchmark_id in STEP_LEVEL_PANEL_ROWS:
        payload = json.loads((BENCH_DIR / f"{benchmark_id}.json").read_text(encoding="utf-8"))
        for compound, data in payload["measured_volatiles"].items():
            key = (compound, float(data["conc_ppb"]), payload["conditions"]["ph"])
            seen.setdefault(key, []).append(benchmark_id)
    duplicates = {key: ids for key, ids in seen.items() if len(ids) > 1}
    assert not duplicates, f"the same measurement is on the panel twice: {duplicates}"


def test_unreachable_rows_are_off_the_panel_and_each_names_its_blocker():
    panel_ids = {path.stem for path in get_benchmark_files()}
    unreachable = sorted(UNREACHABLE_DIR.glob("*.json"))
    assert len(unreachable) == 9, (
        f"expected the nine Wave X unexecutable rows, found {len(unreachable)}"
    )
    for path in unreachable:
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["benchmark_id"] not in panel_ids, (
            f"{path.name} is being scored; the log-error aggregate SKIPS non-positive "
            "predictions, so a row the network cannot execute would cost nothing and the "
            "total miss would be invisible"
        )
        blocker = payload["not_executable"]
        assert blocker["blocker_class"] in {
            "MISSING_STEP", "MISSING_OBSERVABLE", "MISSING_SPECIES_IN_THE_FLUX_TRACKER"
        }
        assert blocker["what_is_missing"] and blocker["step_needed"], (
            f"{path.name} must name what is missing AND what would have to be written; a "
            "blocker without a specification is just a complaint"
        )
        assert payload["measured_volatiles"], "the measurement itself must be preserved"


def test_every_step_level_row_carries_the_measurements_own_precision_contract():
    """Hofmann Table 1 footnote b (+/-10%) is referenced by every table in the paper."""
    for benchmark_id in STEP_LEVEL_PANEL_ROWS:
        payload = json.loads((BENCH_DIR / f"{benchmark_id}.json").read_text(encoding="utf-8"))
        thresholds = payload["validation_contract"]["scale_thresholds"]
        assert thresholds["max_ratio"] == pytest.approx(1.10)
        assert thresholds["mean_abs_log10_error"] == pytest.approx(math.log10(1.10), abs=5e-5)
        for data in payload["measured_volatiles"].values():
            assert data["uncertainty_pct"] == 10
