"""Lipid radical chain regression tests.

REWRITTEN 2026-08-27 (audit remediation, Wave G1 — SMIRKS chemistry soundness).
The chemistry review found the radical chain non-functional in three ways at
once, and two of the tests in this file were passing *because* of the bugs:

  * `radical_propagation_o2` matched `[C;X3]`, i.e. ANY sp2 carbon, so every
    aldehyde, ketone, imine and alkene in the pool was "peroxidised". The old
    `test_full_autoxidation_cycle` asserted that propagation fires — and it did,
    on fabricated substrates. It now uses a model system small enough that the
    genuine peroxy radical survives the MAX_MW = 300 Da prune (see the note in
    `test_propagation_requires_a_real_carbon_radical`).
  * radical flags were lost in `_apply_smirks_rule`'s ordering, so ROO. was
    serialised as ROOH with an invented hydrogen.
  * `Radical_Crosstalk` was retired: it existed to consume the fictitious
    elemental sulfur that the old thiol steps ejected, and its second branch
    quenched radicals by CONSUMING MFT. `test_radical_crosstalk` and
    `test_mft_quenching` are deleted with it.
"""

from rdkit import Chem

from src.smirks_engine import SmirksEngine, Species, ReactionConditions


# 13-hydroperoxy-9,11-octadecadienoic acid (13-HPODE), the linoleate
# hydroperoxide whose alkoxy radical is the accepted hexanal source.
HPODE_13 = r"CCCCCC(OO)/C=C/C=C\CCCCCCCC(=O)O"
# 9-HPODE, whose scission gives 2,4-decadienal and 9-oxononanoic acid instead.
HPODE_9 = r"CCCCC/C=C\C=C/C(OO)CCCCCCCC(=O)O"


def _canonical_products(steps):
    out = set()
    for step in steps:
        for product in step.products:
            mol = Chem.MolFromSmiles(product.smiles)
            if mol is not None:
                out.add(Chem.MolToSmiles(mol))
    return out


def test_lipid_hydroperoxide_scission():
    lipid_h = Species(label="13-HPODE", smiles=HPODE_13)

    engine = SmirksEngine()
    steps = engine.enumerate([lipid_h])

    homolysis_steps = [s for s in steps if s.reaction_family == "Lipid_Homolysis"]
    assert len(homolysis_steps) > 0

    alkoxy = homolysis_steps[0].products[0]
    assert "alkoxy-radical" in alkoxy.label
    assert "[O]" in alkoxy.smiles

    scission_steps = [s for s in steps if s.reaction_family == "Beta_Scission"]
    assert len(scission_steps) > 0


def test_hexanal_is_reachable_from_the_13_hydroperoxide():
    """THE regression this file exists for.

    Before 2026-08-27 the beta-scission SMARTS pinned the beta carbon to
    `[CX4:3]` (sp3). In the 13-alkoxy radical of linoleate the bond whose
    scission releases hexanal is C12-C13, and C12 is sp2, so HEXANAL — the
    model's headline off-note — was STRUCTURALLY UNREACHABLE in the reaction
    network. In production this was masked by a fixed branching ratio applied
    downstream, so nothing failed loudly.
    """
    engine = SmirksEngine()
    steps = engine.enumerate([Species(label="13-HPODE", smiles=HPODE_13)])

    products = _canonical_products(steps)
    assert "CCCCCC=O" in products, (
        "hexanal is not reachable from the 13-hydroperoxide; the allylic "
        "beta-scission channel has regressed"
    )


def test_nine_hydroperoxide_gives_its_own_scission_products_not_hexanal():
    """Positional selectivity: the 9-isomer must NOT give hexanal."""
    engine = SmirksEngine()
    steps = engine.enumerate([Species(label="9-HPODE", smiles=HPODE_9)])

    products = _canonical_products(steps)
    assert "CCCCCC=O" not in products
    # 2,4-decadienal and 9-oxononanoic acid are the accepted 9-HPODE products.
    assert Chem.MolToSmiles(Chem.MolFromSmiles(r"CCCCC/C=C\C=C/C=O")) in products
    assert Chem.MolToSmiles(Chem.MolFromSmiles("O=CCCCCCCCC(=O)O")) in products


def test_beta_scission_fragment_stays_a_radical():
    """The alkyl fragment must be open-shell, not a closed-shell alkane.

    The 1-reactant branch of `_apply_smirks_rule` used to serialise the product
    SMILES BEFORE repairing the radical flags, so the pentyl radical was
    recorded as pentane and the chain terminated silently.
    """
    engine = SmirksEngine()
    steps = engine.enumerate([Species(label="13-HPODE", smiles=HPODE_13)])

    fragments = [
        p.smiles
        for s in steps
        if s.reaction_family == "Beta_Scission"
        for p in s.products
    ]
    radical_fragments = [
        smi
        for smi in fragments
        if any(
            a.GetNumRadicalElectrons() > 0
            for a in Chem.MolFromSmiles(smi).GetAtoms()
        )
    ]
    assert radical_fragments, f"no open-shell beta-scission fragment in {fragments}"


def test_propagation_requires_a_real_carbon_radical():
    """R. + O2 -> ROO. must fire, and ONLY on genuine carbon radicals.

    NOTE for the network: the peroxy radical of a C18 substrate is ~311 Da and
    is pruned by MAX_MW = 300, so the propagation/termination cycle only runs on
    the smaller fragments the chain itself produces. The model system below is
    sized accordingly.
    """
    rad = Species(label="hexyl-radical", smiles="CCCCC[CH2]")
    alkene = Species(label="1-octene", smiles="C=CCCCCCC")
    aldehyde = Species(label="hexanal", smiles="CCCCCC=O")
    o2 = Species(label="O2", smiles="[O]=[O]")

    engine = SmirksEngine()
    steps = engine.enumerate([rad, alkene, aldehyde, o2], max_generations=2)

    propagation = [s for s in steps if s.reaction_family == "Radical_Propagation_O2"]
    assert propagation, "R. + O2 did not fire"

    # The peroxy radical must survive as a radical, not as a hydroperoxide.
    peroxy_smiles = {p.smiles for s in propagation for p in s.products}
    assert any("O[O]" in smi for smi in peroxy_smiles), peroxy_smiles

    # And nothing closed-shell may be peroxidised: every propagation substrate
    # has to carry an unpaired electron.
    for step in propagation:
        carbon_reactant = next(
            r for r in step.reactants if r.smiles not in ("O=O", "[O][O]")
        )
        mol = Chem.MolFromSmiles(carbon_reactant.smiles)
        assert any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms()), (
            f"{carbon_reactant.smiles} is closed-shell but was peroxidised"
        )


def test_abstraction_and_russell_termination_fire_with_o2_present():
    rad = Species(label="hexyl-radical", smiles="CCCCC[CH2]")
    substrate = Species(label="1-octene", smiles="C=CCCCCCC")
    o2 = Species(label="O2", smiles="[O]=[O]")

    engine = SmirksEngine()
    steps = engine.enumerate([rad, substrate, o2], max_generations=3)

    abstraction = [s for s in steps if s.reaction_family == "Peroxy_H_Abstraction"]
    assert abstraction, "ROO. + allylic C-H did not fire"

    termination = [s for s in steps if s.reaction_family == "Radical_Termination"]
    assert termination, "Russell termination did not fire"

    # Russell products: alcohol + carbonyl + O2 (NOT the old `O2 + 2 ROH`,
    # which invented two hydrogens).
    for step in termination:
        products = {p.smiles for p in step.products}
        assert "O=O" in products, products
        assert any("=O" in p and p != "O=O" for p in products), products
        assert any(
            p.endswith("O") and "=" not in p for p in products
        ), f"no alcohol among {products}"
