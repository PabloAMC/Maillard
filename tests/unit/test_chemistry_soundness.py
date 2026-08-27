"""Chemistry-soundness regression tests (audit remediation, Wave G1, 2026-08-27).

These guard the findings of the mechanism-level SMIRKS review recorded in
tasks/audit_remediation.md under "SMIRKS chemistry soundness review":

  * structures — several molecules carried the right LABEL and the wrong
    MOLECULE (bis(2-methyl-3-furyl) disulfide was a mixed 3-furyl/4-furyl
    regioisomer in the template layer while the curated layer had it right;
    DMHF/HEMF had the enol OH on the ring-oxygen carbon). Pinned by InChIKey,
    which is the only thing a wrong regioisomer cannot pass.
  * the accepted MFT route exists and is atom-balanced.
  * free H2 appears only where it is a documented oxidation lumping.
  * elemental sulfur is never ejected as a balance token.
  * atoms and unpaired electrons balance across a whole enumeration.
"""

import collections

import pytest
from rdkit import Chem

from src.conditions import ReactionConditions
from src.curated_pathways import PATHWAYS
from src.pathway_extractor import Species
from src.smirks_engine import (
    SmirksEngine,
    _DMHF_CANONICAL,
    _FURYL_DISULFIDE_CANONICAL,
    _HEMF_CANONICAL,
    _MFT_CANONICAL,
    _NORFURANEOL_CANONICAL,
    _DEOXYOSONE_1_PENTOSE,
    _DEOXYOSONE_1_HEXOSE,
)


# ── Reference structures ────────────────────────────────────────────────────
# InChIKey skeletons (first block) of the verified reference molecules.
_REFERENCE_STRUCTURES = {
    "2-methyl-3-furanthiol": (_MFT_CANONICAL, "RUYNUXHHUVUINQ", "C5H6OS"),
    "bis(2-methyl-3-furyl) disulfide": (
        _FURYL_DISULFIDE_CANONICAL,
        "OHDFENKFSKIFBJ",
        "C10H10O2S2",
    ),
    # 4-hydroxy-2,5-dimethyl-3(2H)-furanone (furaneol)
    "DMHF": (_DMHF_CANONICAL, "INAXVXBDKKUCGI", "C6H8O3"),
    # 2-ethyl-4-hydroxy-5-methyl-3(2H)-furanone
    "HEMF": (_HEMF_CANONICAL, "GWCRPYGYVRXVLI", "C7H10O3"),
    # 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol)
    "norfuraneol": (_NORFURANEOL_CANONICAL, "DLVYTANECMRFGX", "C5H6O3"),
    "1-deoxy-2,3-pentodiulose": (_DEOXYOSONE_1_PENTOSE, "UYTRITJAZOPLCZ", "C5H8O4"),
    "1-deoxy-2,3-hexodiulose": (_DEOXYOSONE_1_HEXOSE, "JPNOFZPGVHGIAU", "C6H10O5"),
}

# Structures that were WRONG before 2026-08-27. Nothing in the codebase may
# carry these skeletons again under the corresponding label.
_RETIRED_WRONG_STRUCTURES = {
    "JQDSBCSEZSRHAY": "mixed 3-furyl/4-furyl disulfide regioisomer",
    "CFMSLBFKBKWCIX": "old DMHF (enol OH on the ring-oxygen carbon)",
    "GMYXRSOIZDSJIT": "old HEMF (enol OH on the ring-oxygen carbon)",
}


def _formula(smiles: str) -> str:
    from rdkit.Chem import rdMolDescriptors

    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles))


def _skeleton(smiles: str) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)).split("-")[0]


@pytest.mark.parametrize(
    "name,smiles,inchikey,formula",
    [(n, s, k, f) for n, (s, k, f) in _REFERENCE_STRUCTURES.items()],
)
def test_reference_structure_identity(name, smiles, inchikey, formula):
    assert _skeleton(smiles) == inchikey, f"{name} is the wrong molecule"
    assert _formula(smiles) == formula, f"{name} has the wrong formula"


def test_template_and_curated_layers_agree_on_the_disulfide():
    """One label, one molecule. The two layers used to disagree."""
    curated = next(
        p.smiles
        for step in PATHWAYS["C_S_Maillard_FFT"]
        for p in step.products
        if p.label == "bis(2-methyl-3-furyl) disulfide"
    )
    assert _skeleton(curated) == _skeleton(_FURYL_DISULFIDE_CANONICAL)
    assert _skeleton(curated) not in _RETIRED_WRONG_STRUCTURES


# ── Network invariants ──────────────────────────────────────────────────────

_COND = ReactionConditions(pH=5.5, temperature_celsius=150.0)

_RIBOSE = Species("D-ribose", "O=CC(O)C(O)C(O)CO")
_GLUCOSE = Species("D-glucose", "O=CC(O)C(O)C(O)C(O)CO")
_CYSTEINE = Species("L-cysteine", "NC(CS)C(=O)O")
_GLYCINE = Species("glycine", "NCC(=O)O")
_LYSINE = Species("L-lysine", "NCCCCC(N)C(=O)O")
_THIAMINE = Species("thiamine", "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1")
_HPODE_13 = Species("13-HPODE", r"CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)O")
_O2 = Species("O2", "[O]=[O]")

_SYSTEMS = {
    "ribose_cysteine": [_RIBOSE, _CYSTEINE],
    "glucose_cysteine_glycine_lysine": [_GLUCOSE, _CYSTEINE, _GLYCINE, _LYSINE],
    "thiamine_cysteine_ribose": [_THIAMINE, _CYSTEINE, _RIBOSE],
    "lipid_oxidation": [_HPODE_13, _O2],
}


def _enumerate(key):
    return SmirksEngine(conditions=_COND).enumerate(_SYSTEMS[key])


def _atoms(species_list):
    counts = collections.Counter()
    charge = 0
    for sp in species_list:
        mol = Chem.MolFromSmiles(sp.smiles)
        assert mol is not None, f"unparseable SMILES {sp.smiles!r}"
        charge += Chem.GetFormalCharge(mol)
        mol = Chem.AddHs(mol)
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] += 1
    return counts, charge


def _radicals(species_list):
    total = 0
    for sp in species_list:
        mol = Chem.MolFromSmiles(sp.smiles)
        total += sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    return total


# Only these families may create or destroy an unpaired-electron PAIR: a
# homolysis makes two radicals from a closed shell, a termination consumes two.
_RADICAL_PAIR_FAMILIES = {"Lipid_Homolysis", "Radical_Termination"}

# Families allowed to emit free H2. Every one of them is an OXIDATION whose
# electron acceptor (O2, a quinone, a dehydro-reductone) the model carries no
# species for; the 2[H] are a documented reducing-equivalent token. See the
# LUMPING NOTEs in src/reaction_templates.py. Anything outside this list means
# hydrogen is being manufactured to patch a stoichiometry.
_H2_PRODUCING_WHITELIST = {
    "Aminoketone_Condensation",   # dihydropyrazine -> pyrazine aromatisation
    "Thiol_Oxidation",            # 2 RSH -> RSSR
    "Lipid_Strecker_Synergy",     # alkylthiazole aromatisation
    "Furan_Ring_Aromatisation",   # thiamine -> MFT furan ring closure
}


@pytest.mark.parametrize("system", sorted(_SYSTEMS))
def test_every_step_conserves_atoms_and_charge(system):
    for step in _enumerate(system):
        r_atoms, r_charge = _atoms(step.reactants)
        p_atoms, p_charge = _atoms(step.products)
        assert r_atoms == p_atoms, (
            f"[{system}] {step.reaction_family} unbalanced: "
            f"{[r.smiles for r in step.reactants]} -> "
            f"{[p.smiles for p in step.products]} ({r_atoms} vs {p_atoms})"
        )
        assert r_charge == p_charge, f"[{system}] {step.reaction_family} charge leak"


@pytest.mark.parametrize("system", sorted(_SYSTEMS))
def test_unpaired_electrons_balance(system):
    for step in _enumerate(system):
        delta = _radicals(step.products) - _radicals(step.reactants)
        if step.reaction_family in _RADICAL_PAIR_FAMILIES:
            assert delta % 2 == 0, (
                f"[{system}] {step.reaction_family} changes the radical count by "
                f"{delta}, which is not a pair"
            )
        else:
            assert delta == 0, (
                f"[{system}] {step.reaction_family} creates/destroys "
                f"{delta} unpaired electron(s): "
                f"{[r.smiles for r in step.reactants]} -> "
                f"{[p.smiles for p in step.products]}"
            )


@pytest.mark.parametrize("system", sorted(_SYSTEMS))
def test_free_hydrogen_only_from_whitelisted_oxidations(system):
    offenders = collections.Counter(
        step.reaction_family
        for step in _enumerate(system)
        if any(p.smiles == "[HH]" for p in step.products)
        and step.reaction_family not in _H2_PRODUCING_WHITELIST
    )
    assert not offenders, (
        f"[{system}] free H2 produced outside the documented oxidation "
        f"lumping whitelist: {dict(offenders)}"
    )


@pytest.mark.parametrize("system", sorted(_SYSTEMS))
def test_no_elemental_sulfur_balance_token(system):
    """`[S]` was ejected purely to balance the old thiol/MFT steps."""
    for step in _enumerate(system):
        for sp in list(step.reactants) + list(step.products):
            assert sp.smiles != "[S]", (
                f"[{system}] elemental sulfur reappeared in "
                f"{step.reaction_family}"
            )


def test_radical_propagation_only_consumes_real_radicals():
    for step in _enumerate("lipid_oxidation"):
        if step.reaction_family != "Radical_Propagation_O2":
            continue
        assert _radicals(step.reactants) > 0, (
            "Radical_Propagation_O2 fired on a closed-shell substrate: "
            f"{[r.smiles for r in step.reactants]}"
        )


# ── The accepted MFT route ──────────────────────────────────────────────────

# 2026-08-27 (Wave N) — RE-PINNED from the norfuraneol route to the
# 1,4-dideoxyosone route. CAUSE: isotope evidence contradicts norfuraneol as
# the in-situ MFT intermediate (Cerny & Davidek 2003, 10.1021/jf026123f:
# norfuraneol spiked into a [13C5]ribose/cysteine system leaves MFT mainly
# 13C5-labelled; Cerny & Davidek 2004, 10.1021/jf035265m: 1,4-dideoxypento-
# 2,3-diulose positionally confirmed). Norfuraneol is still produced
# (Furanone_Cyclisation) but no longer feeds MFT.
def test_pentose_system_builds_mft_through_the_dideoxyosone_route():
    steps = _enumerate("ribose_cysteine")
    families = collections.Counter(s.reaction_family for s in steps)

    assert families["Enolisation_2_3_Amadori"] >= 1, "1-deoxyosone step missing"
    assert families["Furanone_Cyclisation"] >= 1, "norfuraneol step missing"
    assert families["Deoxyosone_Reduction"] >= 1, "1,4-dideoxyosone step missing"
    assert families["Thiol_Addition_Pentodiulose"] >= 1, "dideoxyosone + H2S missing"
    assert families["Thiol_Addition_Norfuraneol"] == 0, (
        "retired norfuraneol->MFT step re-appeared (contradicted by "
        "Cerny & Davidek 2003)"
    )

    labels = {p.label for s in steps for p in s.products}
    assert "norfuraneol" in labels
    assert "1,4-dideoxypentodiulose" in labels
    assert "2-methyl-3-furanthiol" in labels


def test_dideoxyosone_route_steps_are_exactly_balanced():
    # 2026-08-27 (Wave N): family set updated with the route correction; the
    # balance requirement itself is unchanged.
    wanted = {
        "Enolisation_2_3_Amadori",
        "Furanone_Cyclisation",
        "Deoxyosone_Reduction",
        "Thiol_Addition_Pentodiulose",
    }
    route = [s for s in _enumerate("ribose_cysteine") if s.reaction_family in wanted]
    assert len(route) == len(wanted)
    for step in route:
        assert _atoms(step.reactants)[0] == _atoms(step.products)[0], step.reaction_family


def test_pentose_reaches_norfuraneol_without_a_reduction_hexose_does_not():
    """The structural reason pentoses beat hexoses into the MFT branch.

    1-deoxy-2,3-pentodiulose (C5H8O4) -> norfuraneol (C5H6O3) + H2O is a pure
    cyclodehydration. The hexose analogue needs an extra reduction to reach
    furaneol, which is why the classic DMHF precursor is the 6-deoxy sugar
    rhamnose rather than glucose.
    """
    # 2026-08-27 (Wave N): function renamed `_norfuraneol_mft_route` ->
    # `_furanone_and_mft_route` with the MFT route correction (Cerny & Davidek
    # 2003/2004); the structural pentose-vs-hexose claim tested here is
    # unchanged — the reduction-free path to norfuraneol is still pentose-only.
    from src.reaction_templates import _furanone_and_mft_route

    pentose_only = _furanone_and_mft_route(
        [Species("pentose-1-deoxyosone", _DEOXYOSONE_1_PENTOSE)]
    )
    assert [s.reaction_family for s in pentose_only] == ["Furanone_Cyclisation"]

    hexose_no_reductant = _furanone_and_mft_route(
        [Species("hexose-1-deoxyosone", _DEOXYOSONE_1_HEXOSE)]
    )
    assert hexose_no_reductant == []

    hexose_with_reductant = _furanone_and_mft_route(
        [
            Species("hexose-1-deoxyosone", _DEOXYOSONE_1_HEXOSE),
            Species("H2", "[HH]"),
        ]
    )
    assert [s.reaction_family for s in hexose_with_reductant] == ["Furanone_Cyclisation"]


def test_legacy_mft_shortcut_is_labelled_as_such():
    """The fabricated one-step lump must stay distinguishable in any breakdown."""
    families = {s.reaction_family for s in _enumerate("ribose_cysteine")}
    families |= {s.reaction_family for s in _enumerate("glucose_cysteine_glycine_lysine")}
    shortcut = {f for f in families if "Legacy_Shortcut" in f}
    assert shortcut, "the demoted 3-deoxyosone shortcut lost its distinct family name"


# ── Families that used to be silently disabled ──────────────────────────────

def test_previously_dead_families_now_have_real_barriers():
    """All eight fell through to DEFAULT_BARRIER = 45 kcal/mol (~39,000 y half-life)."""
    from src.barrier_constants import DEFAULT_BARRIER, get_barrier

    for family in (
        "Safety_Risk_Acrylamide",
        "Safety_Risk_AGE",
        "Furanone_Formation",
        "Furanone_Cyclisation",
        "Additive_Thermal_Degradation",
        "Generalized_Deamination",
        "Radical_Propagation_O2",
        "Peroxy_H_Abstraction",
        "Radical_Termination",
        "Thiol_Addition_Norfuraneol",
        "Furan_Ring_Aromatisation",
    ):
        barrier, _ = get_barrier(family)
        assert barrier != DEFAULT_BARRIER, f"{family} still falls through to the default"


def test_no_emitted_family_falls_through_to_the_default_barrier():
    from src.barrier_constants import DEFAULT_BARRIER, get_barrier

    emitted = {
        step.reaction_family
        for system in _SYSTEMS
        for step in _enumerate(system)
    }
    emitted |= {
        step.reaction_family for steps in PATHWAYS.values() for step in steps
    }
    fallthrough = sorted(f for f in emitted if get_barrier(f)[0] == DEFAULT_BARRIER)
    assert not fallthrough, (
        "families emitted with no FAST barrier entry (silently kinetically "
        f"disabled): {fallthrough}"
    )


def test_glyoxal_is_reachable_so_cml_and_thiazoles_are_not_dead():
    steps = _enumerate("glucose_cysteine_glycine_lysine")
    labels = {p.label for s in steps for p in s.products}
    assert "glyoxal" in labels, "glyoxal is unreachable; CML/thiazole lanes are dead"

    families = {s.reaction_family for s in steps}
    assert "Safety_Risk_AGE" in families
