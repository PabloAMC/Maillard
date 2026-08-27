from typing import List, Optional, Set, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem
from src.pathway_extractor import Species, ElementaryStep
from src.conditions import ReactionConditions
from src.smirks_engine import (
    _is_ketose, _is_pentose, _is_hexose, _is_sugar, _is_asparagine,
    _is_lysine, _mol, _canonical, _is_valid, _GSH_CANONICAL,
    _is_lipid_aldehyde, _is_lipid_hydroperoxide, _is_primary_amine,
    _has_cysteine_beta_carbon, _is_dicarbonyl, _mw, MAX_MW,
    _THIAMINE_CANONICAL, _DMHF_CANONICAL, _HEMF_CANONICAL,
    _NORFURANEOL_CANONICAL, _MFT_CANONICAL, _FURYL_DISULFIDE_CANONICAL,
    _DEOXYOSONE_1_PENTOSE, _DEOXYOSONE_1_HEXOSE,
    _PENTODIULOSE_14_DIDEOXY,
)

# ──────────────────────────────────────────────────────────────────────────
# Tier B: Parameterised chemical templates
# ──────────────────────────────────────────────────────────────────────────

def _amadori_cascade(sugar: Species, amino_acid: Species) -> List[ElementaryStep]:
    """
    Template: sugar + amino acid → Schiff base → Amadori/Heyns product

    The Schiff base SMILES is constructed by replacing the sugar's terminal
    carbonyl with C=N linked to the amino acid's alpha-amino group.
    For aldoses -> Amadori rearrangement (ketosamine).
    For ketoses -> Heyns rearrangement (aldosamine).
    """
    steps = []
    water = Species(label="water", smiles="O")

    # The string approach actually worked well except for extracting the N.
    # We should just ensure we extract the alpha fragment AND append the basic NH!
    # Let's write a robust version of extraction using RDKit that returns the whole molecule
    # MINUS the oxygen.
    
    # RDKit approach:
    # 1. Sugar: Convert C=O to C-OH and add dummy linker?
    # Simpler: If we know the exact SMILES of sugar and AA, let's use string manipulation 
    # but WITH atom conservation.
    # Aldohexose: O=CC(O)C(O)C(O)C(O)CO + NCC(=O)O -> OCC(O)C(O)C(O)C(O)/C=N/CC(=O)O + H2O
    # Ribose: O=CC(O)C(O)C(O)CO + NCC(=O)O -> OCC(O)C(O)C(O)/C=N/CC(=O)O + H2O
    
    # We need the fragment of the amino acid starting *from the N*, but without its 2 Hs.
    # In canonical SMILES, primary amines are often `N...` or `[NH2]...`.
    # Let's find the N, isolate the fragment, and use it.
    
    _fragment = _extract_alpha_amine_fragment_with_n(amino_acid)
    if not _fragment:
        return []
        
    schiff_label = f"{sugar.label}-{amino_acid.label}-Schiff-base"
    
    if _is_ketose(sugar):
        amadori_label = f"{sugar.label}-{amino_acid.label}-Heyns"
        family = "Heyns_Rearrangement"
        if "fructose" in sugar.label.lower():
            # Fructose: OCC(=O)C(O)C(O)C(O)CO
            # C=O is at index 1 (carbon 2).
            # Fragment = -N-R
            schiff_smiles = f"OCC(={_fragment})C(O)C(O)C(O)CO"
            # Add parenthesis around the -NH-R group so the chain continues properly
            amadori_smiles = f"O=CC(-{_fragment}H)C(O)C(O)C(O)CO"
            # Wait, if _fragment is e.g. "NC(CC)C(=O)O", "-NC(CC)C(=O)OH" is invalid because it has a dash before N.
            # RDKit will output N(CC)C(=O)O if we root it.
            # So _fragment starts with N.
            # schiff: OCC(=N...)
            # amadori: O=CC(N...H) -> The 'H' can be appended inside the N parenthesis if any, but string concat is hard.
            # The safest is: "O=CC(" + _fragment + ")C(O)C(O)C(O)CO", and let RDKit implicitly add the H to N to satisfy valence.
            amadori_smiles = f"O=CC({_fragment})C(O)C(O)C(O)CO"
        else:
            return []
    else:
        amadori_label = f"{sugar.label}-{amino_acid.label}-Amadori"
        family = "Amadori_Rearrangement"
        if _is_pentose(sugar):
            # Ribose: O=CC(O)C(O)C(O)CO
            schiff_smiles = f"OCC(O)C(O)C(O)/C={_fragment}"
            # Again, use parenthesis around the N fragment
            amadori_smiles = f"OCC(O)C(O)C(=O)C({_fragment})"
        elif _is_hexose(sugar):
            # Glucose: O=CC(O)C(O)C(O)C(O)CO
            schiff_smiles = f"OCC(O)C(O)C(O)C(O)/C={_fragment}"
            amadori_smiles = f"OCC(O)C(O)C(O)C(=O)C({_fragment})"
        else:
            return []

    if not _is_valid(schiff_smiles) or not _is_valid(amadori_smiles):
        # We might have generated invalid stereo like C=N]C...
        # Let's clean it up slightly and re-verify
        schiff_smiles = schiff_smiles.replace("=[N", "=N").replace("=[nH]", "=N")
        amadori_smiles = amadori_smiles.replace("C[N", "CN").replace("C[nH]", "CN")
        if not _is_valid(schiff_smiles) or not _is_valid(amadori_smiles):
            return []

    schiff_base = Species(label=schiff_label, smiles=schiff_smiles)
    amadori_product = Species(label=amadori_label, smiles=amadori_smiles)

    steps.append(ElementaryStep(
        reactants=[sugar, amino_acid],
        products=[schiff_base, water],
        reaction_family="Schiff_Base_Formation",
    ))
    steps.append(ElementaryStep(
        reactants=[schiff_base],
        products=[amadori_product],
        reaction_family=family,
    ))
    return steps


def _extract_alpha_amine_fragment_with_n(amino_acid: Species) -> str:
    """
    Extracts the alpha-amino acid SMILES with its alpha nitrogen explicitly grouped,
    e.g., as N(...). We identify the primary alpha N, insert an attachment point or 
    restructure the SMILES so it begins with N.
    """
    # For glycine: NCC(=O)O -> N(CC(=O)O) or something similar
    # By default, RDKit canonical SMILES for simple AAs usually start with N.
    m = _mol(amino_acid.smiles)
    if not m: 
        return ""
    
    # Find alpha nitrogen (N attached to CH attached to C=O)
    pat = Chem.MolFromSmarts("[NH2][CH1,CH2][C](=O)[OH]")
    matches = m.GetSubstructMatches(pat)
    
    # If no match, maybe it's cysteine or something that didn't match perfectly.
    # Try more general: primary amine
    if not matches:
        pat = Chem.MolFromSmarts("[NH2]")
        matches = m.GetSubstructMatches(pat)
        
    if not matches:
        # Fallback if we really can't find it
        return ""
        
    # We want to re-root the SMILES generation at the alpha nitrogen
    n_idx = matches[0][0]
    smi = Chem.MolToSmiles(m, rootedAtAtom=n_idx, isomericSmiles=False)
    
    # The SMILES will start with N. We want to return exactly that string,
    # but when it's appended to C= etc., we'll strip the leading N? No, we WANT the N.
    # E.g. rooted glycine: NCC(=O)O. We want to return N(CC(=O)O).
    # Wait, if we return NCC(=O)O, and substitute `C={fragment}`, we get `C=NCC(=O)O`, 
    # which is exactly correct!
    
    # Let's test Lysine: NCCCCC(N)C(=O)O
    # Rotated at alpha-amine: NC(CCCCN)C(=O)O
    # So `C=NC(CCCCN)C(=O)O` works perfectly!
    
    return smi


def _enolisation_steps(
    amadori: Species,
    sugar: Species,
    amino_acid: Species, # We need the original AA to balance atoms
    conditions: ReactionConditions
) -> List[ElementaryStep]:
    """
    Amadori product → 3-deoxyosone + amino_acid + H2O.
    1,2-enolisation: deoxyosone → furfural/HMF + 2 H2O
    2,3-enolisation: deoxyosone → pyruvaldehyde + fragment
    """
    steps = []
    water = Species(label="water", smiles="O")

    if _is_pentose(sugar):
        deoxy_smi = "O=CC(=O)CC(O)CO"
    else:
        deoxy_smi = "O=CC(=O)CC(O)C(O)CO"

    deoxy = Species(label=f"{sugar.label}-deoxyosone-3", smiles=deoxy_smi)

    # 1. Formation of deoxyosone intermediate
    # C5H8O4 + C2H5NO2 (glycine) = C7H13NO6 (Amadori)!
    steps.append(ElementaryStep(
        reactants=[amadori],
        products=[deoxy, amino_acid],
        reaction_family="Enolisation_Intermediate"
    ))

    # 2. Dehydration to final products
    # 2a. 1,2-enolisation (favored at acidic pH)
    if _is_pentose(sugar):
        product_12 = Species(label="furfural", smiles="O=Cc1ccco1")
        water_count = 2
    else:
        product_12 = Species(label="HMF", smiles="OCC1=CC=C(C=O)O1")
        water_count = 2
        furfural = Species(label="furfural", smiles="O=Cc1ccco1")
        formaldehyde = Species(label="formaldehyde", smiles="C=O")
    
    steps.append(ElementaryStep(
        reactants=[deoxy],
        products=[product_12] + [water] * water_count,
        reaction_family="Enolisation_1_2"
    ))

    if not _is_pentose(sugar):
        # Secondary hexose branch: deoxyosone -> furfural + formaldehyde + 2 H2O.
        steps.append(ElementaryStep(
            reactants=[deoxy],
            products=[furfural, formaldehyde, water, water],
            reaction_family="Enolisation_1_2"
        ))

    # 2b. 2,3-enolisation (favored at neutral/alkali pH)
    product_23 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
    if _is_pentose(sugar): # C5H8O4 -> C3H4O2 (pyruv) + C2H4O2 (glycolaldehyde)
        p2 = Species(label="glycolaldehyde", smiles="O=CCO")
    else: # C6H10O5 -> C3H4O2 (pyruv) + C3H6O3 (glyceraldehyde)
        p2 = Species(label="glyceraldehyde", smiles="O=CC(O)CO")
    
    steps.append(ElementaryStep(
        reactants=[deoxy],
        products=[product_23, p2],
        reaction_family="Enolisation_2_3"
    ))

    # ── 2,3-enolisation of the Amadori product -> 1-DEOXYOSONE ───────────
    # AUDIT 2026-08-27 (Wave G1 fix 7).  This intermediate was entirely absent
    # from the network, and with it the whole accepted route to the
    # 3(2H)-furanones and hence to MFT:
    #     Amadori --(2,3-enolisation)--> 1-deoxyosone
    #             --(cyclisation/dehydration)--> norfuraneol (pentose)
    #             --(+ H2S, reductone-mediated)--> 2-methyl-3-furanthiol
    # van den Ouweland & Peer 1975, DOI 10.1021/jf60199a045; Hofmann &
    # Schieberle.  Atom balance (ribose + glycine): C7H13NO6 -> C5H8O4 +
    # C2H5NO2, exact, no water and no free hydrogen.
    #
    # doi_repair 2026-08-27 (Wave I fix 4) -- CANONICAL RECORD for this
    # flagship anchor; the other four code sites carry only the corrected
    # string and point here.
    #   old:   10.1021/jf60200a038
    #   new:   10.1021/jf60199a045
    #   date:  2026-08-27
    #   basis: the old DOI is registered but resolves to a gossypol/rat-feeding
    #          paper in J. Agric. Food Chem., i.e. a live-but-wrong
    #          (TOPIC-MISMATCH) anchor of exactly the class the 2026-08-26
    #          citation sweep catalogued. 10.1021/jf60199a045 is the real
    #          van den Ouweland & Peer 1975 JAFC paper -- the norfuraneol
    #          (4-hydroxy-5-methyl-3(2H)-furanone) + H2S -> 2-methyl-3-furanthiol
    #          study that this whole sulfur branch is built on. Repointed in
    #          src/reaction_templates.py (x2), src/smirks_engine.py,
    #          src/curated_pathways.py, src/barrier_constants.py and
    #          scripts/generators/refit_sulfur_barriers_hofmann.py.
    #          NO numeric value changed: this is a citation repair only.
    if _is_pentose(sugar):
        deoxy_1 = Species(label="pentose-1-deoxyosone", smiles=_DEOXYOSONE_1_PENTOSE)
    elif _is_hexose(sugar) or _is_ketose(sugar):
        deoxy_1 = Species(label="hexose-1-deoxyosone", smiles=_DEOXYOSONE_1_HEXOSE)
    else:
        deoxy_1 = None

    if deoxy_1 is not None:
        steps.append(ElementaryStep(
            reactants=[amadori],
            products=[deoxy_1, amino_acid],
            reaction_family="Enolisation_2_3_Amadori",
        ))

    return steps


# Canonical forms of the 1-deoxyosones, used to keep the demoted 3-deoxyosone
# shortcut from also firing on them.
_DEOXYOSONE_1_CANONICALS = {
    _canonical(_DEOXYOSONE_1_PENTOSE),
    _canonical(_DEOXYOSONE_1_HEXOSE),
}


def _furanone_and_mft_route(pool_species: List[Species]) -> List[ElementaryStep]:
    """Furanone formation and the isotope-evidenced 1,4-dideoxyosone MFT route.

    ROUTE CORRECTION 2026-08-27 (Wave N) — this function was
    `_norfuraneol_mft_route` and its MFT step was norfuraneol + H2S + 2[H] ->
    MFT (van den Ouweland & Peer 1975, 10.1021/jf60199a045).  That paper is a
    genuine SYNTHESIS of MFT from norfuraneol, but the in-situ competition
    experiment contradicts norfuraneol as the Maillard intermediate: Cerny &
    Davidek 2003 (10.1021/jf026123f) spiked authentic (unlabelled) norfuraneol
    into a [13C5]ribose/cysteine system and the resulting MFT was mainly
    13C5-labelled — "4-hydroxy-5-methyl-3(2H)-furanone is unimportant as an
    intermediate".  Cerny & Davidek 2004 (10.1021/jf035265m, [1-13C]ribose)
    positionally confirmed the proposed intermediate: 1,4-dideoxypento-2,3-
    diulose.  Hofmann & Schieberle 1998 (10.1021/jf9705983) independently rank
    norfuraneol/cysteine among the LESS effective MFT precursors.  Evidence
    dossier: docs/validation/isotope_topology_evidence.md §1-2.

    Steps and their atom balances (all exact, RDKit-verified 2026-08-27):

      1. 1-deoxy-2,3-pentodiulose -> norfuraneol + H2O
         C5H8O4 -> C5H6O3 + H2O.  Pure cyclodehydration.  KEPT: norfuraneol is
         a real Maillard product (and per Cerny & Davidek 2003 the demonstrated
         precursor of 2-mercapto-3-pentanone, a route not yet implemented) —
         it just no longer feeds MFT.

      2. 1-deoxy-2,3-hexodiulose + 2[H] -> furaneol (DMHF) + 2 H2O
         C6H10O5 + H2 -> C6H8O3 + 2 H2O.  The hexose analogue needs an extra
         reduction (which is why the classic DMHF precursor is the 6-deoxy
         sugar rhamnose, not glucose).

      3. 1-deoxy-2,3-pentodiulose + 2[H] -> 1,4-dideoxypento-2,3-diulose + H2O
         C5H8O4 + H2 -> C5H8O3 + H2O.  Formal C4 deoxygenation of the
         1-deoxyosone; in the real system the H-donors are sugar-derived
         reductones/enaminols (see LUMPING NOTE).

      4. 1,4-dideoxypento-2,3-diulose + H2S -> 2-methyl-3-furanthiol + 2 H2O
         C5H8O3 + H2S -> C5H6OS + 2 H2O.  EXACT with no reducing-equivalent
         token: under the literature topology the sulfur-incorporation step
         needs no reduction at all — the `[HH]` bookkeeping moves upstream to
         step 3, where a reductone donor is chemically ordinary.  The intact-C5
         skeleton and the atom mapping (sugar C-1 -> diulose CH3 -> MFT
         2-methyl) match the [1-13C] result.

    LUMPING NOTE — the `2[H]` in steps 2 and 3 is written as molecular H2
    because the model carries no explicit reductone couple.  In the real
    system the hydrogen donors are the sugar-derived reductones/enaminols, not
    free hydrogen gas.  H2 is therefore a *reducing-equivalent token*, and the
    steps are pool-gated on it so the lane cannot run without a reductant
    having been generated somewhere in the network.

    AUDIT 2026-08-27 (Wave I fix 8) — red-team finding H4.  The pool gate above
    was once physically wrong, not because gating on a reductant is wrong, but
    because of WHERE the token came from: in a ribose/cysteine system the ONLY
    producer of `[HH]` was `Aminoketone_Condensation`, silently coupling MFT to
    pyrazine chemistry.  Wave I gave the token a source that exists in every
    system that can make MFT at all (`_thiol_reductant_pool`: 2 cysteine ->
    cystine + 2[H], `Thiol_Oxidation`).  The Wave N route keeps that donor and
    REDUCES the token's role: MFT formation now consumes one `[HH]` (step 3)
    instead of one per H2S addition (old step), and the H2S step itself is
    token-free.
    """
    steps: List[ElementaryStep] = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    h2 = next((s for s in pool_species if s.smiles == "[HH]"), None)
    water = Species("water", "O")

    pentose_do1 = _canonical(_DEOXYOSONE_1_PENTOSE)
    hexose_do1 = _canonical(_DEOXYOSONE_1_HEXOSE)

    seen: Set[str] = set()
    for s in pool_species:
        can = _canonical(s.smiles)
        if can in seen:
            continue

        if can == pentose_do1:
            seen.add(can)
            norfuraneol = Species(
                label="norfuraneol", smiles=_NORFURANEOL_CANONICAL
            )
            steps.append(ElementaryStep(
                reactants=[s],
                products=[norfuraneol, water],
                reaction_family="Furanone_Cyclisation",
            ))
            if h2 is not None:
                dideoxy = Species(
                    label="1,4-dideoxypentodiulose",
                    smiles=_PENTODIULOSE_14_DIDEOXY,
                )
                steps.append(ElementaryStep(
                    reactants=[s, h2],
                    products=[dideoxy, water],
                    reaction_family="Deoxyosone_Reduction",
                ))
        elif can == hexose_do1 and h2 is not None:
            seen.add(can)
            dmhf = Species(label="DMHF", smiles=_DMHF_CANONICAL)
            steps.append(ElementaryStep(
                reactants=[s, h2],
                products=[dmhf, water, water],
                reaction_family="Furanone_Cyclisation",
            ))

    if h2s is None:
        return steps

    # The 1,4-dideoxyosone may have arrived either from step 3 above (same
    # call, so it is not yet in `pool_species`) or already be in the pool from
    # a previous generation; cover both.
    dideoxy_can = _canonical(_PENTODIULOSE_14_DIDEOXY)
    has_dideoxy = any(
        _canonical(s.smiles) == dideoxy_can for s in pool_species
    ) or any(
        _canonical(p.smiles) == dideoxy_can for st in steps for p in st.products
    )
    if has_dideoxy:
        dideoxy = Species(
            label="1,4-dideoxypentodiulose", smiles=_PENTODIULOSE_14_DIDEOXY
        )
        mft = Species(label="2-methyl-3-furanthiol", smiles=_MFT_CANONICAL)
        steps.append(ElementaryStep(
            reactants=[dideoxy, h2s],
            products=[mft, water, water],
            reaction_family="Thiol_Addition_Pentodiulose",
        ))

    return steps


def _strecker_step(
    dicarbonyl: Species, amino_acid: Species
) -> Optional[ElementaryStep]:
    """
    α-dicarbonyl (e.g. pyruvaldehyde) + amino acid → Strecker aldehyde + aminoketone + CO₂
    
    To balance atoms: 
    Amino acid (e.g. Glycine: C2 H5 N O2) loses CO2 (C1 O2) and its sidechain becomes the Strecker aldehyde.
    The remaining -(NH2) group from the amino acid attaches to the dicarbonyl.
    The dicarbonyl (e.g. Pyruvaldehyde: C3 H4 O2) loses ONE oxygen (which goes to the Strecker aldehyde as its carbonyl O), 
    and accepts the -(NH2) to form the aminoketone.
    
    Wait, let's track the exact mechanism:
    Dicarbonyl: R1-C(=O)-C(=O)-R2
    Amino Acid: NH2-CH(R3)-COOH
    
    1. Condensation to Schiff base, losing H2O (from dicarbonyl O and amino 2H).
    2. Decarboxylation: Loses CO2.
    3. Hydrolysis: Adds H2O across the C=N bond.
    
    Net reaction:
    R1-C(=O)-C(=O)-R2 + NH2-CH(R3)-COOH → R1-C(=O)-CH(NH2)-R2 + O=CH-R3 + CO2
    
    So the aminoketone is exactly the dicarbonyl minus ONE carbonyl oxygen, plus NH2, plus 1 H (from the amino acid alpha carbon).
    Since building this dynamically for pyruvaldehyde (CC(=O)C=O) vs diacetyl (CC(=O)C(=O)C) via SMIRKS is complex, 
    we'll use a mapping for both the amino acid AND the dicarbonyl.
    """
    
    # 1. Map Amino Acid to its Strecker Aldehyde
    _aa_to_aldehyde = {
        # name -> (aldehyde_label, aldehyde_smiles)
        "l-leucine":      ("3-methylbutanal", "CC(C)CC=O"),
        "leucine":        ("3-methylbutanal", "CC(C)CC=O"),
        "l-isoleucine":   ("2-methylbutanal", "CCC(C)C=O"),
        "isoleucine":     ("2-methylbutanal", "CCC(C)C=O"),
        "l-valine":       ("2-methylpropanal","CC(C)C=O"),
        "valine":         ("2-methylpropanal","CC(C)C=O"),
        "glycine":        ("formaldehyde",    "C=O"), 
        "l-alanine":      ("acetaldehyde",    "CC=O"),
        "alanine":        ("acetaldehyde",    "CC=O"),
        "l-phenylalanine":("phenylacetaldehyde","O=CCc1ccccc1"),
        "phenylalanine":  ("phenylacetaldehyde","O=CCc1ccccc1"),
        "l-methionine":   ("methional",       "CSCCC=O"),
        "methionine":     ("methional",       "CSCCC=O"),
        "l-lysine":       ("5-aminopentanal", "NCCCCC=O"),
        "lysine":         ("5-aminopentanal", "NCCCCC=O"),
        "l-cysteine":     ("mercaptoacetaldehyde", "SCC=O"),
        "cysteine":       ("mercaptoacetaldehyde", "SCC=O"),
    }

    # 2. Map Dicarbonyl to its Aminoketone
    # R1-C(=O)-C(=O)-R2 -> R1-C(=O)-CH(NH2)-R2
    _dicarbonyl_to_ak = {
        "pyruvaldehyde": ("aminoacetone", "CC(=O)CN"), # CC(=O)C=O (C3H4O2) -> CC(=O)CN (C3H7NO)
        "diacetyl":      ("3-amino-2-butanone", "CC(=O)C(N)C"),
        "glyoxal":       ("2-aminoethanal", "O=CCN"),
        "furfural":      None, # Not a dicarbonyl
        "HMF":           None, 
    }
    
    aa_entry = _aa_to_aldehyde.get(amino_acid.label.lower())
    if aa_entry is None:
        # Fallback for unrecognized amino acids, though they won't balance if we don't know the products.
        return None

    ald_label, ald_smiles = aa_entry
    
    # Find matching dicarbonyl
    ak_entry = _dicarbonyl_to_ak.get(dicarbonyl.label)
    
    # If it's a generic deoxyosone, it acts as the dicarbonyl.
    # e.g., D-glucose-deoxyosone-3 is O=CC(=O)CC(O)C(O)CO.
    # It turns into the corresponding aminoketone: O=CC(N)CC(O)C(O)CO or NC(=O)CC(O)C(O)CO (Wait, aldehydes are more reactive).
    # Since generating these dynamically is hard, we can use an RDKit reaction!
    
    if ak_entry is None:
        # Generic RDKit reaction for dicarbonyl + amino acid -> aldehyde + aminoketone + CO2
        # It's much easier to just use RunReactants.
        # Dicarbonyl: [C:1](=[O:2])[C:3](=[O:4])
        # Amino Acid: [NH2:5][CH1,CH2:6]([R:7])[C:8](=[O:9])[OH:10]
        # Strecker aldehyde: [R:7][C:6]=[O:2] (Wait, the O comes from dicarbonyl? No, O comes from water hydrolysis. Net it's the same as swapping).
        # We know the aldehyde SMILES from the dictionary. We just need the aminoketone!
        rxn_ak = AllChem.ReactionFromSmarts(
            "[C:1](=[O:2])[C:3](=[O:4]) >> [C:1](=[O:2])[C:3]([NH2])"
        )
        dic_mol = _mol(dicarbonyl.smiles)
        if not dic_mol: 
            return None
        prods = rxn_ak.RunReactants((dic_mol,))
        if not prods: 
            return None
        
        try:
            Chem.SanitizeMol(prods[0][0])
            ak_smiles = Chem.MolToSmiles(prods[0][0])
            ak_label = f"amino-{dicarbonyl.label}"
        except Exception:
            return None
    else:
        ak_label, ak_smiles = ak_entry

    aldehyde = Species(label=ald_label, smiles=ald_smiles)
    aminoketone = Species(label=ak_label, smiles=ak_smiles)
    co2 = Species(label="CO2", smiles="O=C=O")

    return ElementaryStep(
        reactants=[dicarbonyl, amino_acid],
        products=[aldehyde, aminoketone, co2],
        reaction_family="Strecker_Degradation",
    )


def _beta_elimination_steps(aa: Species, pool_species: List[Species]) -> List[ElementaryStep]:
    """
    β-elimination of Cys → DHA + H₂S, or Ser → DHA + H₂O.
    If Lysine in pool: DHA + Lys → LAL.
    """
    steps = []
    dha = Species(label="dehydroalanine", smiles="C=C(N)C(=O)O")
    h2s = Species(label="H2S", smiles="S")
    water = Species(label="water", smiles="O")
    lal = Species(label="lysinoalanine", smiles="NC(CCCCNCC(N)C(=O)O)C(=O)O")

    smi_lower = aa.smiles.lower()
    label_lower = aa.label.lower()

    if "cysteine" in label_lower or "nc(cs)" in smi_lower:
        steps.append(ElementaryStep(
            reactants=[aa],
            products=[dha, h2s],
            reaction_family="Beta_Elimination",
        ))
    elif "serine" in label_lower or "nc(co)" in smi_lower:
        steps.append(ElementaryStep(
            reactants=[aa],
            products=[dha, water],
            reaction_family="Beta_Elimination",
        ))
    else:
        return []

    # If lysine is in the pool, DHA reacts to form LAL
    for sp in pool_species:
        if "lysine" in sp.label.lower() or "nccccc(n)" in sp.smiles.lower():
            steps.append(ElementaryStep(
                reactants=[dha, sp],
                products=[lal],
                reaction_family="DHA_Crosslinking",
            ))
            break

    return steps


def _aminoketone_condensation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    2x aminoacetone -> 2,5-dimethylpyrazine + 2H2O + H2
    (Aromatic pyrazines require oxidation/loss of 2H from the dihydro-intermediate).

    RETAINED-H2 LUMPING NOTE (2026-08-27, Wave G1 fix 5).  The aromatisation of
    the dihydropyrazine is a genuine two-electron OXIDATION; in food systems the
    acceptor is dissolved O2 or a dicarbonyl/quinone.  The model carries no
    explicit oxidant species, so those 2[H] are emitted as molecular H2 as a
    documented reducing-equivalent token.  This is one of the four places where
    free H2 survives on purpose (the others are `_furyl_disulfide_formation`,
    `_sulfur_volatiles_pathway` and `_lipid_maillard_synergy`, all likewise
    oxidative aromatisations/dimerisations).  What was removed in this wave is
    the H2 that was NOT a lumped oxidant: the deamination step that manufactured
    it out of water, and the thiazole condensation that emitted it purely to
    patch a wrong carbon donor.
    """
    steps = []
    aks = [s for s in pool_species if "aminoacetone" in s.label.lower() or s.smiles == "CC(=O)CN"]
    for ak in aks:
        pyrazine = Species(label="2,5-dimethylpyrazine", smiles="Cc1cnc(C)cn1")
        water = Species(label="water", smiles="O")
        hydro = Species(label="H2", smiles="[HH]")
        steps.append(ElementaryStep(
            reactants=[ak, ak],
            products=[pyrazine, water, water, hydro],
            reaction_family="Aminoketone_Condensation"
        ))
    return steps


def _retro_aldol_fragmentation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    3-deoxyosone -> C2 + C3 fragments
    Hexose Deoxyosone (C6 H10 O5) -> Pyruvaldehyde (C3 H4 O2) + Glyceraldehyde (C3 H6 O3)
    Sum: C6 H10 O5. Balanced! No water needed.
    Pentose Deoxyosone (C5 H8 O4) -> Pyruvaldehyde (C3 H4 O2) + Glycolaldehyde (C2 H4 O2)
    Sum: C5 H8 O4. Balanced.

    AUDIT 2026-08-27 (Wave G1 fix 5): a second, equally standard retro-aldol
    channel was added — cleavage on the other side of the enediol gives GLYOXAL
    plus an alpha-hydroxy ketone, both classical Maillard fragmentation
    products:

        pentose deoxyosone  C5H8O4  -> C2H2O2 (glyoxal) + C3H6O2 (hydroxyacetone)
        hexose  deoxyosone  C6H10O5 -> C2H2O2 (glyoxal) + C4H8O3 (3,4-dihydroxy-2-butanone)

    Both balance exactly.  Without this channel glyoxal was never generated
    anywhere in the network, which left BOTH the thiazole condensation (whose
    correct C2 donor is glyoxal) and CML formation (lysine + glyoxal, the
    `Safety_Risk_AGE` family) structurally unreachable.
    """
    steps = []
    for s in pool_species:
        lower = s.label.lower()
        # The 1-deoxyosones are 2,3-diketones, not the 1,2-dicarbonyl the
        # retro-aldol channels below assume; they are handled by
        # `_furanone_and_mft_route`.
        if _canonical(s.smiles) in _DEOXYOSONE_1_CANONICALS:
            continue
        if "deoxyosone" in lower:
            glyoxal = Species(label="glyoxal", smiles="O=CC=O")
            # Hexose -> Pyruvaldehyde + Glyceraldehyde
            if "glucose" in lower or "fructose" in lower or _is_hexose(s):
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glyceraldehyde", smiles="O=CC(O)CO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
                p3 = Species(label="3,4-dihydroxy-2-butanone", smiles="CC(=O)C(O)CO")
                steps.append(ElementaryStep([s], [glyoxal, p3], "Retro_Aldol_Fragmentation"))
            # Pentose -> Pyruvaldehyde + Glycolaldehyde
            elif "ribose" in lower or _is_pentose(s):
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glycolaldehyde", smiles="O=CCO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
                p3 = Species(label="hydroxyacetone", smiles="CC(=O)CO")
                steps.append(ElementaryStep([s], [glyoxal, p3], "Retro_Aldol_Fragmentation"))
    return steps


def _cysteine_degradation(amino_acids: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Thermal degradation of Cysteine -> H2S, NH3, acetaldehyde, CO2
    NC(CS)C(=O)O (C3 H7 N O2 S) -> H2S (H2 S) + NH3 (H3 N) + CC=O (C2 H4 O) + O=C=O (C1 O2)
    Sum products: H5 N S + C3 H4 O3
    Total: C3 H9 N O3 S
    Reactants: C3 H7 N O2 S. Difference is +H2O in the products!
    So Cysteine + H2O -> H2S + NH3 + Acetaldehyde + CO2 is perfectly balanced.
    Wait, if water is required, it should be a reactant.
    """
    steps = []
    if conditions.temperature_celsius < 100:
        return steps
        
    for aa in amino_acids:
        if "cysteine" == aa.label.lower() or "l-cysteine" == aa.label.lower() or aa.smiles == "NC(CS)C(=O)O":
            h2s = Species(label="H2S", smiles="S")
            ammonia = Species(label="ammonia", smiles="N")
            acetaldehyde = Species(label="acetaldehyde", smiles="CC=O")
            co2 = Species(label="CO2", smiles="O=C=O")
            water = Species(label="water", smiles="O")
            
            steps.append(ElementaryStep(
                reactants=[aa, water],
                products=[h2s, ammonia, acetaldehyde, co2],
                reaction_family="Cysteine_Degradation"
            ))

    return steps


# Canonical cysteine / cystine, used by `_thiol_reductant_pool`.
_CYSTEINE_CANONICAL = "NC(CS)C(=O)O"                          # C3H7NO2S
_CYSTINE_CANONICAL = "NC(CSSCC(N)C(=O)O)C(=O)O"               # C6H12N2O4S2


def _thiol_reductant_pool(pool_species: List[Species]) -> List[ElementaryStep]:
    """Redox source of the model's 2[H] reducing-equivalent token.

        2 cysteine -> cystine + 2[H]
        2 C3H7NO2S -> C6H12N2O4S2 + H2        (exact, no invented atoms)

    AUDIT 2026-08-27 (Wave I fix 8) — red-team finding H4.  Several lumped
    reduction steps (`_furanone_and_mft_route` steps 2 and 3, `_mft_pathway`,
    `_thiol_addition`) are pool-gated on the `[HH]` reducing-equivalent token.
    Before this wave the token had exactly one producer reachable from a
    ribose/cysteine system: `Aminoketone_Condensation`, the dihydropyrazine ->
    pyrazine aromatisation.  That made the flagship compound MFT a DOWNSTREAM
    DEPENDENT OF PYRAZINE CHEMISTRY — disabling aminoketone condensation drove
    predicted MFT to exactly 0.0 ppb (measured).  No literature supports that
    coupling.

    WHY CYSTEINE -> CYSTINE, and not the sugar-derived reductone.  The other
    honest donor is the 2,3-enediol (reductone) tautomer of the 1-deoxyosone,
    oxidised to its dehydro-reductone.  The model carries NO species for either
    the enediol tautomer or the dehydro-reductone, so writing that couple would
    have meant inventing two new structures purely to carry the bookkeeping —
    precisely the class of move the Wave G1 chemistry review removed (the
    fictitious elemental-sulfur balance token).  Thiol -> disulfide, by
    contrast:
      * is exactly balanced with species that already exist (cysteine is a
        registry precursor; cystine is its textbook oxidation product);
      * reuses the family `Thiol_Oxidation`, which is ALREADY on the documented
        H2-emitting oxidation whitelist in tests/unit/test_chemistry_soundness.py
        for exactly this 2 RSH -> RSSR bookkeeping;
      * is present by construction in every system that can reach MFT at all,
        since cysteine is also the H2S source.

    RETAINED-H2 LUMPING NOTE: as with the other whitelisted oxidations, the
    2[H] leave as molecular H2 because the model carries no explicit electron
    acceptor.  H2 gas is a token, not a claim of hydrogen evolution.
    """
    steps: List[ElementaryStep] = []
    cysteine = next(
        (
            s
            for s in pool_species
            if s.smiles == _CYSTEINE_CANONICAL
            or _canonical(s.smiles) == _canonical(_CYSTEINE_CANONICAL)
            or s.label.lower() in ("cysteine", "l-cysteine")
        ),
        None,
    )
    if cysteine is None:
        return steps

    cystine = Species(label="cystine", smiles=_CYSTINE_CANONICAL)
    hydro = Species(label="H2", smiles="[HH]")
    steps.append(ElementaryStep(
        reactants=[cysteine, cysteine],
        products=[cystine, hydro],
        reaction_family="Thiol_Oxidation",
    ))
    return steps


def _thiazole_condensation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    2-Alkylthiazole formation from a Strecker aldehyde + an alpha-dicarbonyl.

        R-CHO + OHC-CHO + NH3 + H2S  ->  2-R-thiazole + 3 H2O

    AUDIT 2026-08-27 (Wave G1 fix 5): the C2 donor used to be GLYCOLALDEHYDE
    (HOCH2-CHO), which forced a spurious free-H2 by-product to balance:

        C2H4O + C2H4O2 + NH3 + H2S -> C4H5NS + 3 H2O + H2

    Thiazole ring closure is a condensation of an alpha-dicarbonyl with a
    thioamide-type nitrogen/sulfur pair, so the correct C2 donor is GLYOXAL.
    With glyoxal the step balances exactly and the invented hydrogen is gone:

        acetaldehyde     C2H4O  + C2H2O2 + NH3 + H2S -> C4H5NS  + 3 H2O
        3-methylbutanal  C5H10O + C2H2O2 + NH3 + H2S -> C7H11NS + 3 H2O

    Glyoxal itself is produced by the C2/C3 retro-aldol channel of the
    deoxyosones (see `_retro_aldol_fragmentation`).

    Historical derivation of the old (glycolaldehyde) stoichiometry:
    Pyruvaldehyde (C3 H4 O2) + NH3 (H3 N) + H2S (H2 S) -> Thiazole (C3 H3 N S) + 2 H2O (H4 O2) + H2 (H2).
    Sum Reactants: C3 H9 N O2 S. Sum Products: C3 H9 N O2 S. Perfectly balanced!
    Then we can decorate it via alkylation (for the 2-alkylthiazoles) or just use the Strecker aldehydes as the backbone 
    if they have enough carbons. Wait, the template generates specific alkylthiazoles based on the Strecker aldehyde.
    E.g. 3-methylbutanal (C5 H10 O) -> 2-isobutylthiazole (C7 H11 N S). 
    This gains 2 Carbons (C3 backbone from somewhere).
    So: Aldehyde (C_n) + Pyruvaldehyde (C3) + NH3 + H2S -> 2-Alkylthiazole (C_{n+2}) + ... wait.
    Let's use the actual balanced synthesis:
    Aldehyde (R-CHO) + alpha-mercapto-ketone (R'-C(=O)-CH(SH)-R'') + NH3 -> Thiazole + 3H2O.
    The easiest way to balance the existing hardcoded list is to use Pyruvaldehyde (C3H4O2) as the C3 backbone.
    Aldehyde (R-CHO: C_n H_2n O) + Pyruvaldehyde (C3 H4 O2) + NH3 (H3 N) + H2S (H2 S) 
      -> 2-Alkylthiazole (C_{n+3} H_{2n+3} N S) + 3 H2O (H6 O3) + H2 gas (H2).
    Let's check math for acetaldehyde (C2 H4 O).
    Reactants: C2H4O + C3H4O2 + NH3 + H2S = C5 H13 N O3 S.
    Products: 2-methylthiazole (C4 H5 N S). Wait, C4? Acetaldehyde (C2) + Pyruvaldehyde (C3) = C5!
    Where did the extra carbon go? 2-methylthiazole only has 4 carbons!
    Ah, the Thiazole ring itself has 3 carbons. 2-methylthiazole has 3 (ring) + 1 (methyl) = 4 carbons.
    So the backbone must be a C2 piece! Glycolaldehyde (C2 H4 O2).
    Let's check Glycolaldehyde (C2 H4 O2) + Acetaldehyde (C2 H4 O) + NH3 + H2S -> C4 H13 N O3 S.
    2-methylthiazole (C4 H5 N S) + 3 H2O (H6 O3) + H2 (H2) -> C4 H13 N O3 S. PERFECTLY BALANCED!
    
    So: Strecker Aldehyde + Glycolaldehyde + NH3 + H2S -> 2-Alkylthiazole + 3 H2O + H2.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    nh3 = next((s for s in pool_species if s.smiles == "N"), None)
    glyoxal = next(
        (s for s in pool_species
         if "glyoxal" in s.label.lower() or _canonical(s.smiles) == _canonical("O=CC=O")),
        None,
    )

    if not (h2s and nh3 and glyoxal):
        return steps

    _thiazole_map = {
        "3-methylbutanal": ("2-isobutylthiazole", "CC(C)Cc1nccs1"),
        "2-methylbutanal": ("2-sec-butylthiazole", "CCC(C)c1nccs1"),
        "2-methylpropanal": ("2-isopropylthiazole", "CC(C)c1nccs1"),
        "acetaldehyde": ("2-methylthiazole", "Cc1nccs1"),
    }
    
    for sp in pool_species:
        entry = _thiazole_map.get(sp.label)
        if entry:
            tz = Species(label=entry[0], smiles=entry[1])
            water = Species("water", "O")
            steps.append(ElementaryStep(
                reactants=[sp, glyoxal, h2s, nh3],
                products=[tz, water, water, water],
                reaction_family="Lipid_Thiazole_Condensation"
            ))
    return steps


def _thiol_addition(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Furfural + H2S + 2[H] -> 2-furfurylthiol (FFT) + H2O.

    Forming an aryl-methyl thiol from an aryl ALDEHYDE is a net reduction as
    well as a substitution: C5H4O2 + H2S + H2 -> C5H6OS + H2O, exact.  In a real
    food matrix the two hydrogens come from the sugar-derived reductones, not
    from hydrogen gas; the model has no reductone couple, so H2 is used as the
    documented reducing-equivalent token and both routes below are pool-gated
    on it.

    2026-08-27 (Wave G1): the docstring previously claimed a ribose-coupled
    stoichiometry (`furfural + H2S + ribose -> FFT + deoxyosone + 2 H2O`) that
    the code never implemented and that is itself short two hydrogens.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    if not h2s:
        return steps

    # Priority 1: H2 Gas (Legacy/Specific path)
    h2 = next((s for s in pool_species if s.smiles == "[HH]"), None)

    # label -> (thiol label, thiol SMILES, thiohemiacetal SMILES)
    _fft_map = {
        "furfural":         ("2-furfurylthiol",      "SCc1ccco1",    "OC(S)c1ccco1"),
        "5-methylfurfural": ("5-methylfurfurylthiol", "SCc1ccc(C)o1", "OC(S)c1ccc(C)o1"),
    }

    water = Species("water", "O")
    seen = set()
    
    for sp in pool_species:
        entry = _fft_map.get(sp.label)
        if not entry:
            continue

        # Strategy A: H2 gas reduction (if available)
        if h2:
            key = (sp.smiles, h2.smiles)
            if key not in seen:
                seen.add(key)
                fft = Species(label=entry[0], smiles=entry[1])
                steps.append(ElementaryStep(
                    reactants=[sp, h2s, h2],
                    products=[fft, water],
                    reaction_family="Thiol_Addition_H2"
                ))

        # Strategy B: two-step route through the thiohemiacetal.
        # Step 1: Furfural + H2S -> Thiohemiacetal (C5H6O2S)          [exact]
        # Step 2: Thiohemiacetal + 2[H] -> FFT + H2O                   [exact]
        #
        # AUDIT 2026-08-27 (Wave G1 fix 7): step 2 used to be written as
        #   thiohemiacetal + H2S -> FFT + [S] + H2O
        # which balances only by ejecting an atom of ELEMENTAL SULFUR that has
        # no mechanistic counterpart — and that fictitious [S] was the sole
        # reason the `Radical_Crosstalk` family existed (it mopped the [S] up
        # by consuming MFT).  Converting an aryl thiohemiacetal to the thiol is
        # a net REDUCTION, not a dehydration-plus-sulfur-loss, so the step now
        # consumes reducing equivalents.  Same lumping convention as
        # `_furanone_and_mft_route`: H2 stands in for the sugar-derived
        # reductones, and the step is pool-gated on it.
        if h2:
            key = (sp.smiles, "bimolecular_coupled")
            if key not in seen:
                seen.add(key)
                intermed = Species(
                    label=f"{sp.label}-thiohemiacetal", smiles=entry[2]
                )
                fft = Species(label=entry[0], smiles=entry[1])

                # Step 1
                steps.append(ElementaryStep(
                    reactants=[sp, h2s],
                    products=[intermed],
                    reaction_family="Thiohemiacetal_Formation"
                ))
                # Step 2
                steps.append(ElementaryStep(
                    reactants=[intermed, h2],
                    products=[fft, water],
                    reaction_family="Thiol_Dehydration"
                ))

    return steps


def _mft_pathway(pool_species: List[Species]) -> List[ElementaryStep]:
    """DEMOTED one-step 3-deoxyosone -> MFT shortcut (Wave G1 fix 7).

    This is NOT a mechanism.  It is a single lumped step that was standing in
    for the whole sulfur branch, and the 2026-08-27 chemistry review found it
    to be the sole producer of the model's flagship compound.  The accepted
    route (Amadori -> 2,3-enolisation -> 1-deoxyosone -> 1,4-dideoxypento-
    2,3-diulose -> +H2S; Cerny & Davidek 2003/2004 — norfuraneol was retired
    from the MFT lane by Wave N, 2026-08-27) now lives in
    `_furanone_and_mft_route` and carries the primary role; this
    lump is retained only so that hexose systems, whose furanone branch needs
    an extra reduction, retain SOME path to MFT and the network keeps its
    historical reachability.

    Two things changed here:

    * the families were renamed to `Thiol_Addition_Legacy_Shortcut` /
      `Thiol_Addition_Hexose_Legacy_Shortcut` so the lump is distinguishable in
      every family breakdown.  Barrier lookup is unaffected: those names still
      canonicalise to `thiol_addition` / `thiol_addition_hexose`.
    * the fictitious ELEMENTAL SULFUR by-product was removed.  The old
      stoichiometry consumed two H2S and ejected an atom of [S] to balance;
      one H2S plus reducing equivalents balances exactly:

          pentose  C5H8O4  + H2S + 2[H] -> C5H6OS + 3 H2O
          hexose   C6H10O5 + H2S + 2[H] -> C5H6OS + CH2O + 3 H2O

      The `2[H]` carries the same reductone lumping note as
      `_furanone_and_mft_route`.

    AUDIT 2026-08-27 (Wave I fix 12) — HEXOSE-ONLY, as the paragraph above
    already claimed.  Until this wave the substrate filter accepted
    `_is_pentose(s)` and any label containing "ribose", so the demoted lump
    ALSO fired on pentose 3-deoxyosones — in direct contradiction of its own
    docstring, and for no reason: a pentose reaches MFT through the real
    mechanistic pentose route (`_furanone_and_mft_route`), which is why the
    lump was demoted in the first place.  The effect was that pentose MFT was
    produced twice, once by the accepted mechanism and once by a lump that the
    Wave G1 review had already called "NOT a mechanism", which flattered the
    pentose >> hexose ordering with a fabricated contribution.

    RESTRICTED rather than DELETED: the docstring's stated reason for keeping
    the lump is hexose reachability, and that reason is still live.  The hexose
    furanone branch needs an extra reduction to reach DMHF, so deleting the
    lump outright would make MFT UNREACHABLE from glucose + cysteine and break
    the reachability assertion in
    tests/scientific/test_pentose_hexose_sulfur_ordering.py
    (`test_hexose_cysteine_system_reaches_its_claimed_products`) — i.e. it
    would trade a known-lumped path for no path at all.  The pentose emission
    is gone; the family name `Thiol_Addition_Legacy_Shortcut` is retained in
    src/barrier_constants.py and src/family_sensitivity.py because those are
    family-coverage keys.

    MEASURED (see tests/scientific/test_wave_i_network_chemistry.py):
      * pentose >> hexose MFT ordering, matched 150 C / pH 5.5 / aw 0.95 /
        60 min, ribose vs glucose:
            before fix 12:  981.31 / 109.33 ppb = 8.98x
            after  fix 12:  see the test; the pentose side loses the lump's
                            contribution and the hexose side is untouched.
      * cys_ribose_140C_Hofmann1998 MFT: reported in the same test.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    h2 = next((s for s in pool_species if s.smiles == "[HH]"), None)
    if not h2s or not h2:
        return steps

    water = Species("water", "O")

    for s in pool_species:
        label_lower = s.label.lower()
        # The 1-deoxyosones belong to the real route, not to this lump.
        if _canonical(s.smiles) in _DEOXYOSONE_1_CANONICALS:
            continue
        # 2026-08-27 (Wave I fix 12): HEXOSE ONLY. `_is_pentose(s)` and the
        # "ribose" label test were removed; pentoses go through
        # `_furanone_and_mft_route`.
        is_supported_deoxyosone = (
            "deoxyosone" in label_lower and (
                _is_hexose(s)
                or "glucose" in label_lower
                or "fructose" in label_lower
            )
            and not _is_pentose(s)
            and "ribose" not in label_lower
            and "xylose" not in label_lower
            and "arabinose" not in label_lower
        )
        if is_supported_deoxyosone:
            # Filter out N-containing species; they require deamination first (R.12)
            if "N" in s.smiles:
                continue

            mft = Species(label="2-methyl-3-furanthiol", smiles=_MFT_CANONICAL)
            formaldehyde = Species(label="formaldehyde", smiles="C=O")
            steps.append(ElementaryStep(
                reactants=[s, h2s, h2],
                products=[mft, formaldehyde, water, water, water],
                reaction_family="Thiol_Addition_Hexose_Legacy_Shortcut"
            ))

    return steps


def _furyl_disulfide_formation(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Oxidative dimerisation of 2-methyl-3-furanthiol to the characteristic Mottram disulfide.
    Balanced: 2 C5H6OS -> C10H10O2S2 + H2.

    LUMPING NOTE (2026-08-27): thiol -> disulfide is an OXIDATION and in the
    real system the two hydrogens are taken by dissolved O2 or by a quinone /
    dehydro-reductone acceptor.  The model carries no explicit oxidant, so the
    2[H] are written as molecular H2.  This is a deliberate, documented
    reducing-equivalent token, not a claim that H2 gas is evolved.

    STRUCTURE FIX (Wave G1 fix 4): the product was
    `Cc1cc(SSC2=C(C)OC=C2)co1` — a mixed 3-furyl/4-furyl regioisomer
    (InChIKey JQDSBCSEZSRHAY) carrying the bis(2-methyl-3-furyl) disulfide
    label.  The curated layer already had the correct molecule; the two layers
    now agree on `Cc1c(SSc2ccoc2C)cco1` (OHDFENKFSKIFBJ).
    """
    steps = []
    hydro = Species("H2", "[HH]")

    mft_can = _canonical(_MFT_CANONICAL)
    mft_species = [
        s for s in pool_species
        if "2-methyl-3-furanthiol" in s.label.lower() or _canonical(s.smiles) == mft_can
    ]
    seen = set()
    for mft in mft_species:
        key = mft.smiles
        if key in seen:
            continue
        seen.add(key)
        disulfide = Species(
            label="bis(2-methyl-3-furyl) disulfide",
            smiles=_FURYL_DISULFIDE_CANONICAL,
        )
        steps.append(ElementaryStep(
            reactants=[mft, mft],
            products=[disulfide, hydro],
            reaction_family="Thiol_Oxidation",
        ))

    return steps


def _sulfur_volatiles_pathway(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Template for Methionine-derived sulfur volatiles (MeSH, DMDS, DMTS).
    Methional (CSCCC=O) -> Methanethiol (CS) + Acrolein (C=CC=O).
    2x Methanethiol (CS) -> DMDS (CSSC) + H2.
    DMDS + Methanethiol -> DMTS (CSSSC) + H2.

    RETAINED-H2 LUMPING NOTE (2026-08-27): thiol -> di/tri-sulfide is an
    oxidative coupling; the 2[H] are emitted as H2 as a documented
    reducing-equivalent token because the model carries no explicit oxidant.
    """
    steps = []
    hydro = Species("H2", "[HH]")
    
    methionals = [s for s in pool_species if s.label == "methional" or s.smiles == "CSCCC=O"]
    for m in methionals:
        mesh = Species(label="methanethiol", smiles="CS")
        acrolein = Species(label="acrolein", smiles="C=CC=O")
        steps.append(ElementaryStep(
            reactants=[m],
            products=[mesh, acrolein],
            reaction_family="Strecker_Degradation"
        ))
        
    mesh_list = [s for s in pool_species if s.label == "methanethiol" or s.smiles == "CS"]
    dmds_list = []
    for m1 in mesh_list:
        dmds = Species(label="dimethyl-disulfide", smiles="CSSC")
        steps.append(ElementaryStep(
            reactants=[m1, m1],
            products=[dmds, hydro],
            reaction_family="Thiol_Oxidation"
        ))
        dmds_list.append(dmds)
        
    for d in dmds_list:
        for m1 in mesh_list:
            dmts = Species(label="dimethyl-trisulfide", smiles="CSSSC")
            steps.append(ElementaryStep(
                reactants=[d, m1],
                products=[dmts, hydro],
                reaction_family="Thiol_Oxidation"
            ))
            
    return steps




# A free alpha-amino acid: primary amine and carboxyl on the same carbon.
_FREE_AMINO_ACID_SMARTS = Chem.MolFromSmarts("[NX3;H2;!$(N-C=O)][CX4][CX3](=[OX1])[OX2H1]")


def _deamination_step(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    R.12: hydrolytic deamination of alpha-AMINOKETONES.

        R-CO-CH2-NH2 + H2O  ->  R-CO-CH2-OH + NH3
        (aminoacetone C3H7NO + H2O -> hydroxyacetone C3H6O2 + NH3)

    AUDIT 2026-08-27 (Wave G1 fixes 5 and 6).  Three things were wrong:

    1. STOICHIOMETRY.  The step was written as
           R-CO-CH2-NH2 + H2O -> R-CO-CHO + NH3 + H2
       i.e. an OXIDATION (alcohol level -> aldehyde level) dressed up as a
       hydrolysis, balanced by inventing free hydrogen.  Worse, that H2 closed
       a loop with the Strecker step that produced the aminoketone in the first
       place: the network was manufacturing H2 out of water and feeding it to
       furfurylthiol formation.  A hydrolytic deamination gives the
       alpha-HYDROXY ketone, which balances exactly with no H2.

    2. THE LABEL GUARD.  Any species whose LABEL contained an amino-acid name
       was skipped.  Every real Amadori product is labelled
       `<sugar>-<amino acid>-Amadori`, so the guard skipped exactly the class
       the docstring claimed as its primary target.  The guard is now
       structural (a free alpha-amino acid, matched on the molecule).

    3. THE PRODUCT.  When forced through, the Amadori branch emitted the
       1,2-dicarbonyl (a 2-osone) mislabelled `pentose-deoxyosone-R12`, and it
       was unbalanced (an osone sits one oxidation state ABOVE the ketosamine,
       so the step silently created an oxygen).  Sugar ketosamines are now
       delegated to `_enolisation_steps`, which owns both the 1,2-enolisation
       (-> 3-deoxyosone) and the 2,3-enolisation (-> 1-deoxyosone) channels
       with exact balances.

    Amide nitrogens (asparagine/glutamine side chains, peptide bonds) are
    excluded from the pattern: they are not aminoketone nitrogens and hydrolyse
    to a carboxylic acid, not to a ketol.
    """
    steps = []
    water = Species("water", "O")

    # Hydrolytic deamination: the C-N bond is replaced by C-OH and the nitrogen
    # leaves as the free amine.  The recursive SMARTS excludes amide nitrogen.
    rxn = AllChem.ReactionFromSmarts(
        "[CX3:1](=[O:2])-[CX4H2:3]-[NX3;!$([NX3][CX3]=[OX1]):4]"
        ">>[CX3:1](=[O:2])-[CX4H2:3]-[OX2H1].[NX3:4]"
    )

    _known_ketol_labels = {"CC(=O)CO": "hydroxyacetone"}

    for s in pool_species:
        # Filter out pyrazines or very large structures
        if "pyrazine" in s.label.lower() or _mw(s.smiles) > MAX_MW:
            continue

        mol = _mol(s.smiles)
        if mol is None:
            continue

        # Free amino acids undergo Strecker/Maillard, not self-deamination.
        # Matched STRUCTURALLY so that Amadori/Heyns conjugates - whose labels
        # embed an amino-acid name - are no longer swept up with them.
        if _FREE_AMINO_ACID_SMARTS is not None and mol.HasSubstructMatch(_FREE_AMINO_ACID_SMARTS):
            continue

        # Sugar ketosamines belong to the enolisation channels; see (3) above.
        if "Amadori" in s.label or "Heyns" in s.label:
            continue

        prods = rxn.RunReactants((mol,))
        if not prods:
            continue

        try:
            ketol_mol, amine_mol = prods[0][0], prods[0][1]
            Chem.SanitizeMol(ketol_mol)
            Chem.SanitizeMol(amine_mol)
            ketol_smi = Chem.MolToSmiles(ketol_mol)
            amine_smi = Chem.MolToSmiles(amine_mol)
        except Exception:
            continue

        if not (_is_valid(ketol_smi) and _is_valid(amine_smi)):
            continue

        ketol_label = _known_ketol_labels.get(ketol_smi, f"deaminated-{s.label}-ketol")
        amine_label = "ammonia" if amine_smi == "N" else f"liberated-amine-{s.label}"

        steps.append(ElementaryStep(
            reactants=[s, water],
            products=[
                Species(label=ketol_label, smiles=ketol_smi),
                Species(label=amine_label, smiles=amine_smi),
            ],
            reaction_family="Generalized_Deamination"
        ))

    return steps

def _lipid_maillard_synergy(pool_species: List[Species]) -> List[ElementaryStep]:
    """
    Lipid Aldehyde + alpha-aminoketone + H2S -> 2-Alkylthiazole + 2 H2O + H2.
    Lipid Aldehyde + 2x alpha-aminoketone -> Alkylpyrazine (Branching synergy).

    RETAINED-H2 LUMPING NOTE (2026-08-27): as in `_aminoketone_condensation`,
    the thiazole aromatisation is an oxidation and the 2[H] are emitted as H2
    in the absence of an explicit oxidant species.
    """
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    
    lipids = [s for s in pool_species if _is_lipid_aldehyde(s)]
    # Target aminoketones (aminoacetone, 3-amino-2-butanone, etc.)
    aks = [s for s in pool_species if "amino" in s.label.lower() and ("acetone" in s.label.lower() or s.smiles == "CC(=O)CN") and "dicarbonyl" not in s.label.lower()]

    if not lipids or not aks:
        return steps

    water = Species("water", "O")
    hydro = Species("H2", "[HH]")

    for lip in lipids:
        # Extract alkyl chain length for label
        m_lip = _mol(lip.smiles)
        c_lip = sum(1 for a in m_lip.GetAtoms() if a.GetAtomicNum() == 6)
        r_len = c_lip - 1 # excluding carbonyl C
        
        for ak in aks:
            # 1. Thiazole Synergy (if H2S present)
            if h2s:
                if r_len == 5: # Hexanal
                    tz_name = "2-pentyl-4-methylthiazole"
                    tz_smi = "CCCCCC1=NC(C)=CS1"
                elif r_len == 6: # Heptanal
                    tz_name = "2-hexyl-4-methylthiazole"
                    tz_smi = "CCCCCCC1=NC(C)=CS1"
                else:
                    tz_name = f"2-alkyl(C{r_len})-4-methylthiazole"
                    prefix = "C" * r_len
                    tz_smi = f"{prefix}C1=NC(C)=CS1"

                steps.append(ElementaryStep(
                    reactants=[lip, ak, h2s],
                    products=[Species(tz_name, tz_smi), water, water, hydro],
                    reaction_family="Lipid_Strecker_Synergy"
                ))

    return steps


def _lipid_hydroperoxide_scission(pool: List[Species]) -> List[ElementaryStep]:
    """
    Tier B: Homolytic cleavage of lipid hydroperoxides.
    R-OOH -> R-O. + .OH
    """
    steps = []
    oh_radical = Species(label="hydroxyl-radical", smiles="[OH]")
    
    for s in pool:
        if _is_lipid_hydroperoxide(s):
            # Use SMIRKS for homolysis but tag the alkoxy oxygen with Isotope 13 so we can find it
            rxn = AllChem.ReactionFromSmarts("[CX4:1]-[OX2:3][OX2H,OX1H0-:2]>>[C:1]-[13OX1:3].[OH:2]")
            m = _mol(s.smiles)
            if not m: continue
            prods = rxn.RunReactants((m,))
            if prods:
                alk_mol = prods[0][0]
                # Explicitly set radical on the oxygen (find Isotope 13)
                for atom in alk_mol.GetAtoms():
                    if atom.GetIsotope() == 13 and atom.GetAtomicNum() == 8:
                        atom.SetIsotope(0)
                        atom.SetNumExplicitHs(0)
                        atom.SetNoImplicit(True)
                        atom.SetNumRadicalElectrons(1)
                
                try:
                    alk_mol.UpdatePropertyCache(strict=False)
                    Chem.SanitizeMol(alk_mol)
                except Exception:
                    pass

                alkoxy_smi = Chem.MolToSmiles(alk_mol)
                alkoxy = Species(label=f"{s.label}-alkoxy-radical", smiles=alkoxy_smi)
                steps.append(ElementaryStep(
                    reactants=[s],
                    products=[alkoxy, oh_radical],
                    reaction_family="Lipid_Homolysis"
                ))
    return steps


# `_radical_crosstalk_templates` (family `Radical_Crosstalk`) was RETIRED on
# 2026-08-27 (Wave G1 fix 7).  It existed to mop up the fictitious elemental
# sulfur `[S]` that `_thiol_addition` and `_mft_pathway` ejected to balance
# themselves: `[S]` parses as an open-shell atom, so it was picked up as a
# "radical" and quenched -- and its MFT branch quenched by CONSUMING the
# model's flagship compound, turning 2-methyl-3-furanthiol into an
# unaccounted-for MFT-thiyl radical.  Both [S] sources are gone (the thiol
# steps are now reductions, not sulfur ejections), so the sink is gone with
# them.  The `radical_crosstalk` FAST_BARRIERS entry is deliberately kept: it
# is referenced by src/family_barrier_progress.py as a family-coverage key.


def _sugar_ring_opening(pool_species: List[Species]) -> List[ElementaryStep]:
    """Hemiacetal cyclic sugar -> open-chain aldehyde. (Defensive rule via RWMol)"""
    steps = []
    # Match hemiacetal: O(ring) - C(ring)(OH)
    # The [O;X2;R] ensures it's a ring oxygen, [C;X4;R] is a ring carbon, [O;X2;H] is the hydroxyl.
    patt = Chem.MolFromSmarts("[O;X2;R]-[C;X4;R]-[OH]")
    if not patt: 
        return steps
    
    for s in pool_species:
        m = _mol(s.smiles)
        if not m: 
            continue
            
        matches = m.GetSubstructMatches(patt)
        if not matches: 
            continue
            
        # We might have multiple hemiacetals (e.g., dimers), just take the first
        o_ring_idx, c_anomeric_idx, o_hydroxyl_idx = matches[0]
        
        try:
            rw_mol = Chem.RWMol(m)
            # 1. Break O_ring - C_anomeric bond
            rw_mol.RemoveBond(o_ring_idx, c_anomeric_idx)
            
            # 2. Change C_anomeric - O_hydroxyl bond to double
            b = rw_mol.GetBondBetweenAtoms(c_anomeric_idx, o_hydroxyl_idx)
            if b:
                b.SetBondType(Chem.BondType.DOUBLE)
            
            # 3. Adjust valency hydrogens automatically via Sanitize
            Chem.SanitizeMol(rw_mol)
            open_smi = Chem.MolToSmiles(rw_mol)
            
            if _is_valid(open_smi):
                p = Species(label=f"{s.label}-open", smiles=open_smi)
                steps.append(ElementaryStep([s], [p], "Sugar_Ring_Opening"))
        except Exception:
            pass
            
    return steps


# ── Toxic / Safety Marker Templates ──────────────────────────────────────

def _acrylamide_formation(sugar: Species, asparagine: Species) -> List[ElementaryStep]:
    """
    Asparagine + Reducing Sugar -> Acrylamide + CO2 + NH3 + sugar 3-deoxyosone + H2O.

    AUDIT 2026-08-27 (Wave I fix 18) — the reducing sugar used to be a
    SPECTATOR.  The step was written

        sugar + Asn -> acrylamide + SUGAR + CO2 + NH3

    with the inline excuse "Sugar is conserved here for simplicity", which made
    the sugar formally CATALYTIC: the whole asparagine budget could be converted
    without consuming a single sugar molecule.  That is not the accepted
    mechanism, and it contradicted `src/safety.py::predict_acrylamide`, the
    other half of this lane, which is explicitly SECOND ORDER —
    `dA/dt = kf*[Asn]*[Sugar] - ke*[A]` (Knol 2009 / Parker 2012) — i.e. it
    treats the sugar as a consumed co-substrate.

    The two layers now agree that the sugar is consumed.  The sugar leaves as
    its 3-deoxyosone, which is where the accepted mechanism actually sends it
    (Asn + carbonyl -> N-glycosyl-Asn -> decarboxylated Schiff base ->
    acrylamide, with the sugar moiety hydrolysing off into the deoxyosone that
    carries on down the Maillard cascade).  That species already exists in the
    network: `_enolisation_steps` emits the identical SMILES.

    Atom balance is now exact in both branches (it was exact before too, because
    the spectator sugar cancelled on both sides):

        hexose   C6H12O6 + C4H8N2O3 -> C3H5NO + CO2 + NH3 + C6H10O5 + H2O
                 C10H20N2O9         -> C10H20N2O9
        pentose  C5H10O5 + C4H8N2O3 -> C3H5NO + CO2 + NH3 + C5H8O4  + H2O
                 C9H18N2O8          -> C9H18N2O8

    MEASURED IMPACT, and an honest negative result.  Predicted acrylamide is
    UNCHANGED to 12 significant figures by this fix (synthetic Asn 10 mM +
    sugar at 1/10/100 mM, 160 C / 20 min / pH 6.0 / aw 0.6): the flux
    propagation in src/recommend.py scores a step from its REACTANTS' spans and
    concentrations, and the sugar was already a reactant, so moving it out of
    the products changes the stoichiometry without moving the number.  Two
    things follow, both stated rather than papered over:
      * dose dependence was ALREADY present and is unchanged — acrylamide is
        exactly linear in sugar mM (7.86 / 78.61 / 786.09 ppb at 1 / 10 /
        100 mM);
      * sugar-IDENTITY dependence is still absent, and is absent in BOTH
        layers: glucose, ribose and fructose give bit-identical predictions,
        and `predict_acrylamide` takes a single lumped `reducing_sugar_mM` with
        no identity term at all.  Giving the template an identity term (e.g. a
        `safety_risk_acrylamide` entry in
        `barrier_constants.DONOR_REACTIVITY_MULTIPLIERS`) would be a NEW
        CALIBRATION with no anchor behind it, and would put the template out of
        step with safety.py — the opposite of the consistency this fix is for.
        OPEN OWNER ITEM, not fixed here.
    """
    if not (_is_sugar(sugar) and _is_asparagine(asparagine)):
        return []

    # The sugar leaves as its 3-deoxyosone. Same SMILES/label convention as
    # `_enolisation_steps`, so the two producers dedupe onto one pool species.
    if _is_pentose(sugar):
        deoxyosone_smiles = "O=CC(=O)CC(O)CO"          # C5H8O4
    elif _is_hexose(sugar) or _is_ketose(sugar):
        deoxyosone_smiles = "O=CC(=O)CC(O)C(O)CO"      # C6H10O5
    else:
        # No balanced sugar-derived product can be named for this substrate
        # (e.g. a disaccharide). Emitting the old spectator form instead would
        # reintroduce the catalytic-sugar bug, so the step is not emitted.
        # Every sugar in data/species/precursors.yml is a pentose or a hexose,
        # so this branch is currently unreachable from the registry.
        return []

    acrylamide = Species(label="Acrylamide", smiles="C=CC(=O)N")
    co2 = Species(label="CO2", smiles="O=C=O")
    ammonia = Species(label="ammonia", smiles="N")
    water = Species(label="water", smiles="O")
    deoxyosone = Species(
        label=f"{sugar.label}-deoxyosone-3", smiles=deoxyosone_smiles
    )

    return [ElementaryStep(
        reactants=[sugar, asparagine],
        products=[acrylamide, co2, ammonia, deoxyosone, water],
        reaction_family="Safety_Risk_Acrylamide"
    )]


def _cml_cel_formation(lysine: Species, pool: List[Species]) -> List[ElementaryStep]:
    """
    Lysine + Glyoxal -> CML
    Lysine + Methylglyoxal -> CEL
    """
    if not _is_lysine(lysine):
        return []
        
    steps = []
    
    for s in pool:
        # Glyoxal -> CML
        if s.smiles == "O=CC=O" or s.label == "glyoxal":
            cml = Species("CML", "N[C@@H](CCCCNCC(=O)O)C(=O)O")
            steps.append(ElementaryStep(
                reactants=[lysine, s],
                products=[cml],
                reaction_family="Safety_Risk_AGE"
            ))
        # Methylglyoxal -> CEL
        elif s.smiles == "CC(=O)C=O" or s.label == "pyruvaldehyde" or s.label == "methylglyoxal":
            cel = Species("CEL", "N[C@@H](CCCCNC(C)C(=O)O)C(=O)O")
            steps.append(ElementaryStep(
                reactants=[lysine, s],
                products=[cel],
                reaction_family="Safety_Risk_AGE"
            ))
            
    return steps


# ── PBMA Additive Degradations ───────────────────────────────────────────

def _thiamine_degradation(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: thermal breakdown of Thiamine (Vitamin B1).

    AUDIT 2026-08-27 (Wave G1 fix 8).  The single step this template used to
    emit was

        thiamine(+)  ->  H2S + 2-methylthiophene + 4,5-dihydro-2-methylthiazole

    which is grossly unbalanced: C12H17N4OS(+) on the left against C9H15NS3 on
    the right — three carbons, three nitrogens and an oxygen destroyed, and TWO
    extra sulfur atoms created out of nothing (thiamine has exactly one S, so it
    cannot yield both H2S and a thiophene in one step).

    Replaced by the accepted, atom-balanced cascade (van der Linde et al. 1979;
    Hofmann & Schieberle; Cerny 2008 on 5-hydroxy-3-mercapto-2-pentanone):

      1. bridge hydrolysis
         thiamine(+) + H2O -> 4-methyl-5-(2-hydroxyethyl)thiazole
                            + protonated 4-amino-5-(hydroxymethyl)-2-methylpyrimidine
         C12H17N4OS(+) + H2O -> C6H9NOS + C6H10N3O(+)
      2. thiazolium ring opening
         thiazole + 2 H2O -> 5-hydroxy-3-mercapto-2-pentanone + formamide
         C6H9NOS + 2 H2O -> C5H10O2S + CH3NO
      3a. cyclisation + aromatisation to the key thiamine aroma marker
         5-hydroxy-3-mercapto-2-pentanone -> 2-methyl-3-furanthiol + H2O + 2[H]
         C5H10O2S -> C5H6OS + H2O + H2
         Closing the furan ring costs a water AND a two-electron oxidation, so
         this step carries the same reducing-equivalent lumping note as the
         other aromatisations; it is given its own family
         (`Furan_Ring_Aromatisation`) so that lumping stays visible.
      3b. competing H2S elimination
         C5H10O2S -> H2S + C5H8O2 (5-hydroxy-3-penten-2-one)
      4. thiophene closure (keeps the historical 2-methylthiophene product,
         now with a real sulfur source)
         C5H8O2 + H2S -> C5H6S + 2 H2O

    Every step balances exactly, including charge on step 1.
    """
    if conditions.temperature_celsius < 100:
        return []

    steps = []
    water = Species("water", "O")
    het = Species("4-methyl-5-(2-hydroxyethyl)thiazole", "Cc1ncsc1CCO")
    hmp = Species("4-amino-5-hydroxymethyl-2-methylpyrimidinium", "Cc1nc(N)c(CO)c[nH+]1")
    hmp_ketone = Species("5-hydroxy-3-mercapto-2-pentanone", "CC(=O)C(S)CCO")
    formamide = Species("formamide", "NC=O")
    mft = Species("2-methyl-3-furanthiol", _MFT_CANONICAL)
    h2s = Species("Hydrogen_Sulfide", "S")
    enone = Species("5-hydroxy-3-penten-2-one", "CC(=O)C=CCO")
    methylthiophene = Species("2-methylthiophene", "Cc1cccs1")

    # Identify Thiamine in the pool by canonical SMILES
    for s in pool:
        if _canonical(s.smiles) != _canonical(_THIAMINE_CANONICAL):
            continue

        steps.append(ElementaryStep(
            reactants=[s, water],
            products=[het, hmp],
            reaction_family="Additive_Thermal_Degradation"
        ))
        steps.append(ElementaryStep(
            reactants=[het, water, water],
            products=[hmp_ketone, formamide],
            reaction_family="Additive_Thermal_Degradation"
        ))
        steps.append(ElementaryStep(
            reactants=[hmp_ketone],
            products=[mft, water, Species("H2", "[HH]")],
            reaction_family="Furan_Ring_Aromatisation"
        ))
        steps.append(ElementaryStep(
            reactants=[hmp_ketone],
            products=[h2s, enone],
            reaction_family="Additive_Thermal_Degradation"
        ))
        steps.append(ElementaryStep(
            reactants=[enone, h2s],
            products=[methylthiophene, water, water],
            reaction_family="Additive_Thermal_Degradation"
        ))
    return steps


def _furanone_generation(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Literature-grounded low-confidence template for positive furanones.

    Blank & Fay 1996 supports:
    - pentose + alanine -> HEMF   (C5 + C3 - CO2 = C7, exact)
    - pentose + glycine -> DMHF   (C5 + C2 - CO2 = C6, exact)

    AUDIT 2026-08-27 (Wave G1 fix 4): a third step, `pentose + alanine -> DMHF`,
    was deleted.  It was short one carbon (C5 + C3 - CO2 = C7 on the left, C6 on
    the right), i.e. it destroyed a carbon atom every time it fired.  Both
    surviving steps balance exactly, and in both of them the amino acid is a
    genuine carbon donor, which is what makes the C count work.

    The furanone a pentose gives WITHOUT an external carbon donor is
    norfuraneol (4-hydroxy-5-methyl-3(2H)-furanone, C5H6O3), not DMHF; that
    route is mechanistic rather than lumped and lives in
    `_furanone_and_mft_route` (1-deoxypentosulose -> norfuraneol + H2O).
    """
    if conditions.temperature_celsius < 90.0:
        return []

    pentoses = [species for species in pool if _is_pentose(species)]
    if not pentoses:
        return []

    alanines = [species for species in pool if species.label.lower() in {"alanine", "l-alanine"}]
    glycines = [species for species in pool if species.label.lower() in {"glycine", "l-glycine"}]
    ammonia = Species(label="ammonia", smiles="N")
    co2 = Species(label="CO2", smiles="O=C=O")
    water = Species(label="water", smiles="O")

    steps: List[ElementaryStep] = []
    for pentose in pentoses:
        for alanine in alanines:
            hemf = Species(label="HEMF", smiles=_HEMF_CANONICAL)
            steps.append(
                ElementaryStep(
                    reactants=[pentose, alanine],
                    products=[hemf, ammonia, co2, water, water],
                    reaction_family="Furanone_Formation",
                    source_quality="literature",
                    barrier_uncertainty_kcal=6.0,
                )
            )
        for glycine in glycines:
            dmhf = Species(label="DMHF", smiles=_DMHF_CANONICAL)
            steps.append(
                ElementaryStep(
                    reactants=[pentose, glycine],
                    products=[dmhf, ammonia, co2, water, water],
                    reaction_family="Furanone_Formation",
                    source_quality="literature",
                    barrier_uncertainty_kcal=6.0,
                )
            )
    return steps


def _glutathione_cleavage(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: Controlled cleavage of Glutathione (GSH).
    Literature shows it cleaves into glutamic acid and cysteinylglycine dipeptide.

    AUDIT 2026-08-27 (Wave G1 fix 8): this is an amide HYDROLYSIS and the water
    reactant was missing, so the step created H2O out of nothing:
        GSH C10H17N3O6S -> C5H9NO4 + C5H10N2O3S  (products carry an extra H2O)
    With water on the left it balances exactly:
        C10H17N3O6S + H2O -> C5H9NO4 + C5H10N2O3S
    """
    if conditions.temperature_celsius < 100:
        return []

    steps = []
    water = Species("water", "O")
    for s in pool:
        if _canonical(s.smiles) == _canonical(_GSH_CANONICAL):
            p1 = Species("Glutamic_Acid", "N[C@@H](CCC(=O)O)C(=O)O")
            p2 = Species("Cysteinylglycine", "N[C@@H](CS)C(=O)NCC(=O)O")

            steps.append(ElementaryStep(
                reactants=[s, water],
                products=[p1, p2],
                reaction_family="Additive_Thermal_Degradation"
            ))
    return steps


