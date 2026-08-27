"""Wave P (2026-08-27) — the six evidence-grounded chemistry changes, pinned.

Every assertion here pins a number the wave MEASURED, not a number it wanted. Where a
number got worse it is pinned worse, with the cause in the comment. The point of this
file is that a future change to any of the six lanes has to move a named pin and say why.

The six items:
  1. `thiol_addition_pentodiulose` refitted against Hofmann1998 (28.60 -> 26.35)
  2. the Hofmann & Schieberle C2 + C3 recombination lane to MFT
  3. norfuraneol -> 2,3-pentanedione -> 2-mercapto-3-pentanone
  4. nonanal cleaved from the OLEATE pool, not the linoleate pool
  5. fructose reaches HMF by its own ring-retained dehydration
  6. the DMHF `[HH]` pool gate removed in favour of the amino-acid reduction
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest
from rdkit import Chem

from src.barrier_constants import (
    DEFAULT_BARRIER,
    FAST_BARRIERS,
    _canonical_fast_family,
    get_barrier,
)
from src.benchmark_validation import evaluate_benchmark
from src.conditions import ReactionConditions
from src.lipid_oxidation import (
    MARKER_HYDROPEROXIDE_POOL,
    PEA_LIPID_PROFILE,
    SOY_LIPID_PROFILE,
    predict_hexanal_generation,
)
from src.smirks_engine import SmirksEngine, Species

_COND = ReactionConditions(pH=5.5, temperature_celsius=150.0)

_RIBOSE = Species("D-Ribose", "O=CC(O)C(O)C(O)CO")
_GLUCOSE = Species("D-Glucose", "O=CC(O)C(O)C(O)C(O)CO")
_FRUCTOSE = Species("D-Fructose", "OCC(O)C(O)C(O)C(=O)CO")
_CYSTEINE = Species("L-Cysteine", "NC(CS)C(=O)O")
_GLYCINE = Species("Glycine", "NCC(=O)O")

_HOFMANN = Path("data/benchmarks/cys_ribose_140C_Hofmann1998.json")

#: The six families Wave P introduced.
_NEW_FAMILIES = (
    "Mercaptoketone_Formation",
    "Mercaptoketone_Aldol_Addition",
    "Mercaptoketone_Cyclodehydration",
    "Furanone_Reductive_Opening",
    "Furanone_Amino_Acid_Reduction",
    "Fructofuranosyl_Dehydration",
)


def _enumerate(precursors, generations: int = 3):
    return SmirksEngine(conditions=_COND).enumerate(list(precursors), generations)


def _atoms(species_list):
    counts: collections.Counter = collections.Counter()
    for sp in species_list:
        mol = Chem.MolFromSmiles(sp.smiles)
        assert mol is not None, f"unparseable SMILES {sp.smiles!r}"
        for atom in Chem.AddHs(mol).GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts


# ── family plumbing ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("family", _NEW_FAMILIES)
def test_new_families_have_explicit_barrier_keys_and_do_not_fall_through(family):
    """The canonicalisation fallthrough is a KNOWN defect class in this table.

    Wave G1 fix 8 found eight emitted families silently inheriting
    ``DEFAULT_BARRIER = 45.0`` (a ~39,000-year half-life at 150 C), i.e. switched off
    without anyone noticing. `_canonical_fast_family` also has substring rules that
    would have caught two of these names in the wrong bucket. Both halves are pinned:
    the key resolves to ITSELF (not to a class analogue), and the value is not the
    default.
    """
    key = _canonical_fast_family(family)
    assert key == family.lower(), (
        f"{family} canonicalises to {key!r}, i.e. it is inheriting another family's "
        "barrier through a substring rule instead of carrying its own explicit key"
    )
    assert key in FAST_BARRIERS, f"{family} has no explicit FAST_BARRIERS entry"
    assert get_barrier(family)[0] != DEFAULT_BARRIER, (
        f"{family} falls through to DEFAULT_BARRIER"
    )
    rationale = FAST_BARRIERS[key][1]
    assert "ESTIMATED" in rationale, (
        f"{family}'s barrier must be labelled ESTIMATED — none of the Wave P values "
        "has a kinetic measurement behind it"
    )


def test_new_families_collect_no_accidental_ph_or_water_activity_correction():
    """`src/conditions.py` classifies reactions for the pH/aw corrections by SUBSTRING.

    A family name containing "strecker", "amadori", "schiff", "pyrazine" or "furfural"
    silently picks up corrections that can be worth hundreds of x. The Wave P families
    are named so that they receive the SAME treatment as the sibling steps they extend
    or replace (`Furanone_Cyclisation`, `Furanone_Formation`), which receive none.
    This test is the guard on that: it fails the day someone renames one of them into a
    trigger word and changes a prediction by accident.
    """
    for family in _NEW_FAMILIES:
        fam = family.lower()
        assert _COND._ionization_correction(fam) == 1.0, (
            f"{family} picked up an amine-ionisation correction from its NAME"
        )
        assert _COND._water_activity_correction(fam) == 1.0, (
            f"{family} picked up a Labuza water-activity correction from its NAME"
        )


# ── item 1: the refit ──────────────────────────────────────────────────────────

def test_pentodiulose_barrier_is_the_wave_p_fit_and_carries_the_conversion_caveat():
    """26.35 kcal/mol, FITTED, against an anchor whose unit conversion is UNVERIFIED.

    Wave N shipped 28.60 (the un-fitted `thiol_addition` class value) and declined to
    refit, because refitting re-couples the only sulfur anchor to a fitted constant.
    Wave P did it with owner approval, AFTER the wave's chemistry additions so the fit
    sees the network that ships. The load-bearing part is not the number: it is that
    the rationale must carry the fit target's own `content_verification_note` verbatim,
    so nobody can read 26.35 as a measured barrier.
    """
    value, rationale = FAST_BARRIERS["thiol_addition_pentodiulose"]
    assert value == pytest.approx(26.35, abs=1e-9)
    assert "FITTED 2026-08-27 (Wave P item 1)" in rationale
    assert "cys_ribose_140C_Hofmann1998" in rationale
    # The verbatim caveat, and the sentence that says what it means for THIS constant.
    assert "mol%->ppb conversion" in rationale
    assert "NOT documented anywhere in this repo" in rationale
    assert "UNVERIFIED" in rationale
    assert "LOCALISED HERE" in rationale
    # The boundary honesty: the profile minimum is at the range floor, so the adopted
    # value is the conservative edge and the residual is NOT removable by this barrier.
    assert "PROFILE MINIMUM SITS AT THE RANGE FLOOR" in rationale
    assert "sulfur_barrier_refit_pentodiulose" in rationale


def test_hofmann1998_after_the_refit():
    """Both rows improved; both are still misses; the benchmark is still a fit target.

    Movement across Wave P, measured 2026-08-27:
        MFT  151.87 -> 242.38 ppb vs 342   (2.2519x under -> 1.4110x under)
        FFT  243.72 -> 217.98 ppb vs 200   (1.2186x over  -> 1.0900x over)
    FFT was NOT fitted. It co-moves because the two lanes draw on the same upstream
    sugar flux — which is exactly why ONE knob was fitted against these two rows and
    not two.
    """
    evaluation = evaluate_benchmark(_HOFMANN)
    predicted = {c.compound: c.predicted_ppb for c in evaluation.comparisons}
    assert predicted["2-methyl-3-furanthiol"] == pytest.approx(242.38, rel=1e-3)
    assert predicted["2-furfurylthiol"] == pytest.approx(217.99, rel=1e-3)
    # Still under / still over — the fit did not flip either direction.
    assert predicted["2-methyl-3-furanthiol"] < 342.0
    assert predicted["2-furfurylthiol"] > 200.0


# ── item 2: the C2 + C3 lane, and the finding it produced ──────────────────────

def test_c2_c3_recombination_lane_exists_and_is_exactly_balanced():
    """Hofmann & Schieberle 1998's highest-yielding MFT system, finally expressible."""
    steps = _enumerate([_RIBOSE, _CYSTEINE])
    families = [s.reaction_family for s in steps]
    for family in (
        "Mercaptoketone_Formation",
        "Mercaptoketone_Aldol_Addition",
        "Mercaptoketone_Cyclodehydration",
    ):
        assert family in families, f"{family} does not fire in ribose + cysteine"
    for step in steps:
        if step.reaction_family.startswith("Mercaptoketone"):
            assert _atoms(step.reactants) == _atoms(step.products), (
                f"{step.reaction_family} unbalanced"
            )


def test_hexose_system_also_reaches_the_c2_c3_lane():
    """Measured, and it corrected a guess.

    The lane needs glycolaldehyde, which `_retro_aldol_fragmentation` emits only from
    the PENTOSE deoxyosone — so the lane was expected to be pentose-only. It is not:
    in glucose + cysteine the identical species arrives by
    glyoxal -> Strecker -> 2-aminoethanal -> hydrolytic deamination. Pool identity is
    by canonical SMILES, not by label, which is what connects the longer route.
    """
    families = {s.reaction_family for s in _enumerate([_GLUCOSE, _CYSTEINE])}
    assert "Mercaptoketone_Cyclodehydration" in families


def test_the_second_mft_channel_contributes_exactly_zero_to_the_prediction():
    """THE FINDING OF THIS WAVE, and it is about the propagator, not the chemistry.

    `src/recommend.py` relaxes to the LOWEST-SPAN path per product; it does not sum
    fluxes over parallel channels. Measured on cys_ribose_140C_Hofmann1998 at the
    shipped constants, 2026-08-27:

        pentodiulose lane alone   242.38 ppb
        C2+C3 lane alone           71.02 ppb
        both lanes                242.38 ppb   <- the max, not the sum

    Had the two added, MFT would read 313.39 ppb, i.e. 1.09x under instead of 1.41x
    under. So adding a real, literature-evidenced second channel to the flagship
    compound moved the flagship number by exactly nothing. That is reported, not
    fixed: changing the propagator to an additive one is a model-wide recalibration
    event and belongs to its own wave.

    This test pins the ZERO. If someone makes the propagator additive, it fails, and
    that failure is the notification.
    """
    import src.smirks_engine as engine_module

    baseline = {
        c.compound: c.predicted_ppb for c in evaluate_benchmark(_HOFMANN).comparisons
    }
    original = engine_module._c2_c3_mft_recombination
    try:
        engine_module._c2_c3_mft_recombination = lambda pool: []
        without = {
            c.compound: c.predicted_ppb
            for c in evaluate_benchmark(_HOFMANN).comparisons
        }
    finally:
        engine_module._c2_c3_mft_recombination = original

    assert without["2-methyl-3-furanthiol"] == pytest.approx(
        baseline["2-methyl-3-furanthiol"], rel=1e-12
    ), (
        "the C2+C3 lane now contributes to predicted MFT. If that is because the flux "
        "propagator became additive, this is the expected failure — re-pin it and say "
        "so in the ledger."
    )


# ── item 3: norfuraneol's real sulfur fate ─────────────────────────────────────

def test_norfuraneol_has_consumers_again_and_they_are_the_evidenced_ones():
    """Wave N retired norfuraneol from the MFT lane and left it a dead end.

    Cerny & Davidek 2003 (10.1021/jf026123f) assign it 2-mercapto-3-pentanone and NOT
    MFT; Whitfield & Mottram 1999 (10.1021/jf980980v) supply the 2,3-pentanedione
    intermediate. The isomer matters: 3-mercapto-2-pentanone has a second,
    ribose-derived origin, so this lane must NOT emit it.
    """
    steps = _enumerate([_RIBOSE, _CYSTEINE])
    consumers = {
        s.reaction_family
        for s in steps
        if any(r.label == "norfuraneol" for r in s.reactants)
    }
    assert consumers == {"Furanone_Reductive_Opening"}

    products = {p.label for s in steps for p in s.products}
    assert "2,3-pentanedione" in products
    assert "2-mercapto-3-pentanone" in products
    assert "3-mercapto-2-pentanone" not in products, (
        "this lane must not emit the isomer Cerny & Davidek assign partly to ribose"
    )

    # And nothing on this lane may reach MFT: that was the Wave N correction.
    for step in steps:
        if step.reaction_family in {"Furanone_Reductive_Opening"}:
            assert not any(p.label == "2-methyl-3-furanthiol" for p in step.products)


# ── item 4: nonanal comes off the oleate pool ──────────────────────────────────

def test_nonanal_scales_with_oleate_and_hexanal_with_linoleate():
    """`oleic_acid_pct` was declared, populated, registered — and read by nothing."""
    assert MARKER_HYDROPEROXIDE_POOL["Nonanal"] == "oleate"
    assert MARKER_HYDROPEROXIDE_POOL["Hexanal"] == "linoleate"
    assert MARKER_HYDROPEROXIDE_POOL["2-Pentylfuran"] == "linoleate"

    out = predict_hexanal_generation(PEA_LIPID_PROFILE, temp_C=100.0, time_min=45.0)
    ratio = PEA_LIPID_PROFILE.oleic_acid_pct / PEA_LIPID_PROFILE.linoleic_acid_pct
    assert out["total_hydroperoxide_oleate"] == pytest.approx(
        out["total_hydroperoxide"] * ratio, rel=1e-12
    )
    # The branching ratio itself was NOT touched; only the pool it multiplies.
    assert out["nonanal"] == pytest.approx(out["total_hydroperoxide_oleate"] * 0.15, rel=1e-12)


def test_trikusuma_nonanal_recovery_broke_by_exactly_the_oleic_linoleic_ratio():
    """An honest regression, pinned as one, and NOT refitted away.

    `pea_iso/heated_matrix/nonanal` observable_factor was back-solved so this benchmark
    reproduced its own measured 24 ppb — against the pre-Wave-P linoleate-pool nonanal.
    After the substrate correction the row reads 10.56 ppb, i.e. 2.2727x under, which is
    EXACTLY 1 / (22.0 / 50.0). That exact arithmetic identity is the evidence that the
    factor was absorbing a substrate-assignment error and nothing else. Refitting it
    would put the error straight back into the same constant, so it was left alone; see
    the dated note on the record in src/matrix_calibration_registry.py.
    """
    evaluation = evaluate_benchmark(
        Path("data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json")
    )
    predicted = {c.compound: c.predicted_ppb for c in evaluation.comparisons}
    expected_ratio = PEA_LIPID_PROFILE.oleic_acid_pct / PEA_LIPID_PROFILE.linoleic_acid_pct
    assert predicted["nonanal"] == pytest.approx(24.0 * expected_ratio, rel=1e-4)
    # The two rows whose factors were NOT substrate-affected still recover exactly.
    assert predicted["hexanal"] == pytest.approx(782.0, rel=1e-4)
    assert predicted["2-pentylfuran"] == pytest.approx(163.0, rel=1e-4)


def test_soy_nonanal_uses_the_soy_oleic_fraction():
    out = predict_hexanal_generation(SOY_LIPID_PROFILE, temp_C=100.0, time_min=45.0)
    ratio = SOY_LIPID_PROFILE.oleic_acid_pct / SOY_LIPID_PROFILE.linoleic_acid_pct
    assert ratio == pytest.approx(23.0 / 53.0)
    assert out["total_hydroperoxide_oleate"] == pytest.approx(
        out["total_hydroperoxide"] * ratio, rel=1e-12
    )


# ── item 5: fructose reaches HMF by its own route ──────────────────────────────

def test_fructose_reaches_hmf_without_the_glucose_3_deoxyosone():
    """Perez Locas & Yaylayan 2008 exclude 3-deoxyglucosone for fructose.

    The distinction is TOPOLOGICAL: fructose keeps its own ring (Antal 1990). The Heyns
    product and the 3-deoxyosone are still formed — they feed the retro-aldol and
    Strecker lanes — but they no longer carry fructose's HMF.
    """
    steps = _enumerate([_FRUCTOSE, _GLYCINE])
    hmf_producers = {
        s.reaction_family for s in steps if any(p.label == "HMF" for p in s.products)
    }
    assert hmf_producers == {"Fructofuranosyl_Dehydration"}

    direct = next(s for s in steps if s.reaction_family == "Fructofuranosyl_Dehydration")
    assert [r.label for r in direct.reactants] == ["D-Fructose"]
    assert _atoms(direct.reactants) == _atoms(direct.products)

    # The glucose limb is untouched: it still reaches HMF through the 3-deoxyosone.
    glucose_hmf = {
        s.reaction_family
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "HMF" for p in s.products)
    }
    assert glucose_hmf == {"Enolisation_1_2"}


# ── item 6: the DMHF [HH] gate is gone ─────────────────────────────────────────

def test_furaneol_from_glucose_glycine_survives_disabling_the_pyrazine_lane():
    """Red-team H4, second half — the exact measurement Wave L1 published as a defect.

    Wave L1, monkey-patching `_aminoketone_condensation` to return no steps:
        [baseline]                    16 steps | DMHF steps 1 | H2 producers: [pyrazine]
        [aminoketone cond. DISABLED]  14 steps | DMHF steps 0 | H2 producers: []
    i.e. predicted furaneol from glucose was a downstream dependent of pyrazine
    chemistry, while Wang & Ho 2008's CAMOLA experiment establishes the route as real.

    After Wave P the token is not re-sourced, it is GONE: the amino acid is the
    reductant (Blank & Fay 1996; Kerler et al. 2010), the step balances exactly with
    no `[HH]`, and DMHF survives.
    """
    import src.smirks_engine as engine_module

    def dmhf_steps():
        return [
            s
            for s in _enumerate([_GLUCOSE, _GLYCINE])
            if any(p.label == "DMHF" for p in s.products)
        ]

    baseline = dmhf_steps()
    assert len(baseline) == 1

    original = engine_module._aminoketone_condensation
    try:
        engine_module._aminoketone_condensation = lambda pool: []
        without_pyrazine = dmhf_steps()
    finally:
        engine_module._aminoketone_condensation = original

    assert len(without_pyrazine) == 1, (
        "furaneol from glucose/glycine is contingent on pyrazine chemistry again"
    )


def test_the_dmhf_step_consumes_no_reducing_equivalent_token():
    steps = [
        s
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "DMHF" for p in s.products)
    ]
    assert steps, "no DMHF step at all"
    for step in steps:
        assert step.reaction_family == "Furanone_Amino_Acid_Reduction"
        assert not any(r.smiles == "[HH]" for r in step.reactants), (
            "the DMHF step is pool-gated on the reducing-equivalent token again"
        )
        assert _atoms(step.reactants) == _atoms(step.products)
        # The amino acid is the reductant: it leaves as its Strecker aldehyde + CO2 + NH3.
        labels = {p.label for p in step.products}
        assert {"DMHF", "formaldehyde", "CO2", "ammonia", "water"} <= labels


def test_hexose_dmhf_carbons_all_come_from_the_sugar():
    """Wang & Ho 2008's CAMOLA constraint: the C6 skeleton stays intact.

    A 1:1 [13C6]/[12C6]glucose mixture gave a 1:1 [13C6]/[12C6]DMHF mixture with no
    intermediate isotopologues. The amino-acid-coupled step preserves that: the amino
    acid's carbon leaves as its own aldehyde and as CO2, so all six DMHF carbons are
    still the hexose's.
    """
    steps = [
        s
        for s in _enumerate([_GLUCOSE, _GLYCINE])
        if any(p.label == "DMHF" for p in s.products)
    ]
    step = steps[0]
    sugar_carbons = _atoms(
        [r for r in step.reactants if r.label == "hexose-1-deoxyosone"]
    )["C"]
    assert sugar_carbons == 6
    dmhf_carbons = _atoms([p for p in step.products if p.label == "DMHF"])["C"]
    assert dmhf_carbons == 6
