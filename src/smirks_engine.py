"""
src/smirks_engine.py — Hybrid SMIRKS + Template Maillard Rule Engine

Phase 6.1: Replaces the static curated_pathways.py dictionary with a
rule-based engine that automatically enumerates Maillard reaction pathways
for any arbitrary sugar + amino acid precursor pool.

Architecture:
  Tier B (Templates): Complex multi-step rearrangements encoded as Python
    functions — Amadori cascade, Strecker degradation, enolisation
    branching, beta-elimination (DHA). Exact and chemically grounded.

  Tier A (SMIRKS): Simple, tight functional-group transforms applied
    iteratively on the growing intermediate pool — Schiff base formation
    from lipid aldehydes, thiol addition to generate FFT-type thiols.
    Guarded by abundant MW capping and canonical SMILES deduplication.

Output: List[ElementaryStep] — fully compatible with xtb_screener.py
        and the existing Tier 1 / recommend.py pipeline.
"""

import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple
from dataclasses import dataclass, field

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Suppress RDKit atom-mapping warnings from SMIRKS rules that use unmapped atoms
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.warning")

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

from src.pathway_extractor import Species, ElementaryStep
from src.conditions import ReactionConditions

# ──────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────

MAX_MW = 300.0  # Daltons — prune products above this (volatiles are small)

# Additive Canonical SMILES (for exact matching)
_THIAMINE_CANONICAL = "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"
_GSH_CANONICAL = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"

# Tier A SMIRKS rules: (name, reaction_family, smirks, ph_gate)
# ph_gate: "any" | "acid" (pH<6) | "neutral" (pH>=6)
#
# These rules are deliberately narrow:
# - Lipid_Schiff_Base: only aliphatic C3+ aldehydes with a plain alkyl carbon (no OH)
#   adjacent to the carbonyl. This excludes sugars, which have C(O) alpha-carbons
#   and are handled exclusively by the Tier B Amadori template.
# - Thiol_Addition: only aromatic (furanyl) aldehydes + H2S → thiol.
_SMIRKS_RULES: List[Tuple[str, str, str, str]] = [
    (
        "schiff_base_lipid",
        "Lipid_Schiff_Base",
        # C3+ aliphatic aldehyde whose alpha-carbon has NO hydroxyl (excludes sugars).
        # The amine donor can be anything with a primary amine on an sp3 carbon (like amino acids).
        "[CX4H2,CX4H3:2][CH1:1]=O.[NH2:3][CX4:4]>>[C:2][C:1]=[N:3]-[C:4].O",
        "any",
    ),
    (
        "thiol_addition_furfural",
        "Thiol_Addition",
        # Aromatic aldehyde (furfural-type) + H2S -> thiol + water
        "[c:1][CH:2]=O.[SH2]>>[c:1][CH2:2]S.O",
        "any",
    ),
]


# ──────────────────────────────────────────────────────────────────────────
# Structural classification helpers
# ──────────────────────────────────────────────────────────────────────────

# SMARTS for classifying precursors
_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[CH1]=O")          # aliphatic aldehyde
_POLYOL_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[CH1]=O")   # reused on sugars
_AMINO_SMARTS = Chem.MolFromSmarts("[NH2][CX4]")          # primary amine, not amide
_THIOL_SMARTS = Chem.MolFromSmarts("[SH]")                # thiol or H2S
_DICARBONYL_SMARTS = Chem.MolFromSmarts("[C](=O)[C](=O)")  # adjacent carbonyls
_AROMATIC_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[c][CH]=O") # furfural-type


def _mol(smi: str) -> Optional[Chem.Mol]:
    return Chem.MolFromSmiles(smi) if smi else None


def _canonical(smi: str) -> Optional[str]:
    m = _mol(smi)
    return Chem.MolToSmiles(m) if m else None


def _mw(smi: str) -> float:
    m = _mol(smi)
    return Descriptors.MolWt(m) if m else 9999.0


def _is_valid(smi: str) -> bool:
    """Return True if SMILES is parseable and MW is below cap."""
    return _mol(smi) is not None and _mw(smi) <= MAX_MW


def _is_sugar(s: Species) -> bool:
    """Heuristic: has an aldehyde OR ketone AND at least 2 hydroxyl groups."""
    m = _mol(s.smiles)
    if m is None:
        return False
    has_aldehyde = m.HasSubstructMatch(_ALDEHYDE_SMARTS)
    has_ketone = m.HasSubstructMatch(Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]"))
    # Count OH groups (O with at least one H)
    oh_count = sum(
        1 for atom in m.GetAtoms()
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1
        and atom.GetDegree() == 1  # terminal OH
    )
    return (has_aldehyde or has_ketone) and oh_count >= 2


def _is_ketose(s: Species) -> bool:
    """Heuristic: has a ketone C=O and multiple OH."""
    m = _mol(s.smiles)
    if not m: return False
    pat = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]")
    has_ketone = m.HasSubstructMatch(pat)
    oh_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1)
    return has_ketone and oh_count >= 2


def _is_hexose(s: Species) -> bool:
    """Heuristic: 6 carbons + is a sugar."""
    m = _mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return _is_sugar(s) and c_count == 6


def _is_pentose(s: Species) -> bool:
    """Heuristic: 5 carbons + is a sugar."""
    m = _mol(s.smiles)
    if m is None:
        return False
    c_count = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6)
    return _is_sugar(s) and c_count == 5


def _is_primary_amine(s: Species) -> bool:
    """Has a primary amine not in amide context (amino acids, amines)."""
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_AMINO_SMARTS)


def _has_thiol(s: Species) -> bool:
    """Has thiol group (–SH) or is H₂S."""
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_THIOL_SMARTS)


def _has_cysteine_beta_carbon(s: Species) -> bool:
    """
    Detects cysteine-like β-carbon with thiol (Cys, Ser).
    Pattern: N[C][C][SH or OH] (2-aminol/thiol-3-carbon skeleton).
    """
    pat = Chem.MolFromSmarts("[NH2][CH1][CH2][SX2H,OX2H]")
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(pat)


def _is_dicarbonyl(s: Species) -> bool:
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_DICARBONYL_SMARTS)


def _is_aromatic_aldehyde(s: Species) -> bool:
    m = _mol(s.smiles)
    return m is not None and m.HasSubstructMatch(_AROMATIC_ALDEHYDE_SMARTS)


def _species_from_pool(pool: Set[str], label: str, smiles: str) -> Species:
    """Create a Species, canonicalise its SMILES, and add to pool."""
    can = _canonical(smiles)
    if can:
        pool.add(can)
    return Species(label=label, smiles=smiles)


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

    # Build Schiff base label
    schiff_label = f"{sugar.label}-{amino_acid.label}-Schiff-base"
    
    if _is_ketose(sugar):
        # Heyns route
        amadori_label = f"{sugar.label}-{amino_acid.label}-Heyns"
        family = "Heyns_Rearrangement"
        # Fructose pattern: OCC(=O)C(O)C(O)C(O)CO
        if "fructose" in sugar.label.lower():
            schiff_smiles = f"OCC(=N{_extract_alpha_amine_fragment(amino_acid)})C(O)C(O)C(O)CO"
            amadori_smiles = f"O=CC(N{_extract_alpha_amine_fragment(amino_acid)})C(O)C(O)C(O)CO"
            deoxyosone_smiles = "O=CC(=O)CC(O)C(O)CO"
        else:
            return []
    else:
        # Amadori route
        amadori_label = f"{sugar.label}-{amino_acid.label}-Amadori"
        family = "Amadori_Rearrangement"
        if _is_pentose(sugar):
            schiff_smiles = f"OCC(O)C(O)C(O)/C=N/{_extract_alpha_amine_fragment(amino_acid)}"
            amadori_smiles = f"OCC(O)C(O)C(=O)CN{_extract_alpha_amine_fragment(amino_acid)}"
            deoxyosone_smiles = "O=CC(=O)CC(O)CO"
        elif _is_hexose(sugar):
            schiff_smiles = f"OCC(O)C(O)C(O)C(O)/C=N/{_extract_alpha_amine_fragment(amino_acid)}"
            amadori_smiles = f"OCC(O)C(O)C(O)C(=O)CN{_extract_alpha_amine_fragment(amino_acid)}"
            deoxyosone_smiles = "O=CC(=O)CC(O)C(O)CO"
        else:
            return []

    schiff_base = Species(label=schiff_label, smiles=schiff_smiles)
    amadori_product = Species(label=amadori_label, smiles=amadori_smiles)
    deoxyosone = Species(label=f"{sugar.label}-deoxyosone-3", smiles=deoxyosone_smiles)

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


def _extract_alpha_amine_fragment(amino_acid: Species) -> str:
    """
    Extract the fragment attached to the amine of an alpha-amino acid
    for SMILES construction. Returns the SMILES of R–CH(NH–)–COOH minus NH2.
    E.g. glycine NCC(=O)O → CC(=O)O (alpha side chain + carboxyl).
    This is used to append to the imine carbon in Schiff base SMILES.
    """
    # Simple heuristic: remove the terminal NH2 from the SMILES
    smi = amino_acid.smiles
    # Replace first occurrence of NH2-connected C with just the carbon
    # This is a string approximation — chemically valid for linear amino acids
    for prefix in ["N[C@@H]", "N[C@H]", "NC"]:
        if smi.startswith(prefix):
            return smi[len("N"):]  # drop leading N
    return smi  # fallback


def _enolisation_steps(
    amadori: Species,
    sugar: Species,
    conditions: ReactionConditions
) -> List[ElementaryStep]:
    """
    Amadori product → 3-deoxyosone → dehydrated product.
    1,2-enolisation (pH<6): pentose → furfural, hexose → HMF
    2,3-enolisation (pH>=6): both → pyruvaldehyde (dicarbonyl)
    """
    steps = []
    water = Species(label="water", smiles="O")

    if _is_pentose(sugar):
        deoxy_smi = "O=CC(=O)CC(O)CO"
    else:
        deoxy_smi = "O=CC(=O)CC(O)C(O)CO"

    deoxy = Species(label=f"{sugar.label}-deoxyosone-3", smiles=deoxy_smi)

    # 1. Formation of deoxyosone intermediate
    steps.append(ElementaryStep(
        reactants=[amadori],
        products=[deoxy, water],
        reaction_family="Enolisation_Intermediate"
    ))

    # 2. Dehydration to final product
    if conditions.pH < 6:
        if _is_pentose(sugar):
            product = Species(label="furfural", smiles="O=Cc1ccco1")
        else:
            product = Species(label="HMF", smiles="OCC1=CC=C(C=O)O1")
        family = "Enolisation_1_2"
    else:
        product = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
        family = "Enolisation_2_3"

    steps.append(ElementaryStep(
        reactants=[deoxy],
        products=[product, water, water], # 2 more waters lost
        reaction_family=family,
    ))

    return steps


def _strecker_step(
    dicarbonyl: Species, amino_acid: Species
) -> Optional[ElementaryStep]:
    """
    α-dicarbonyl (e.g. pyruvaldehyde) + amino acid → Strecker aldehyde + aminoketone + CO₂

    The Strecker aldehyde is one carbon shorter than the amino acid's side chain.
    We use a lookup table for the standard amino acids.
    """
    _strecker_map = {
        # amino acid smiles pattern → (aldehyde_label, aldehyde_smiles, aminoketone_smiles)
        "l-leucine":      ("3-methylbutanal", "CC(C)CC=O",    "CC(=O)CN"),
        "l-isoleucine":   ("2-methylbutanal", "CCC(C)C=O",    "CC(=O)CN"),
        "l-valine":       ("2-methylpropanal","CC(C)C=O",     "CC(=O)CN"),
        "glycine":        ("acetaldehyde",    "CC=O",          "CC(=O)CN"),
        "l-alanine":      ("acetaldehyde",    "CC=O",          "CC(=O)CN"),
        "l-phenylalanine":("phenylacetaldehyde","O=CCc1ccccc1","CC(=O)CN"),
        "l-methionine":   ("methional",       "O=CCCS",        "CC(=O)CN"),
    }

    entry = _strecker_map.get(amino_acid.label.lower())
    if entry is None:
        return None  # amino acid not in Strecker map, skip

    ald_label, ald_smiles, ak_smiles = entry
    aldehyde = Species(label=ald_label, smiles=ald_smiles)
    aminoketone = Species(label="aminoacetone", smiles=ak_smiles)
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
    """2x aminoacetone -> 2,5-dimethylpyrazine + 3H2O"""
    steps = []
    aks = [s for s in pool_species if "aminoacetone" in s.label.lower() or s.smiles == "CC(=O)CN"]
    for ak in aks:
        pyrazine = Species(label="2,5-dimethylpyrazine", smiles="Cc1cnc(C)cn1")
        water = Species(label="water", smiles="O")
        steps.append(ElementaryStep(
            reactants=[ak, ak],
            products=[pyrazine, water, water, water],
            reaction_family="Aminoketone_Condensation"
        ))
    return steps


def _retro_aldol_fragmentation(pool_species: List[Species]) -> List[ElementaryStep]:
    """3-deoxyosone -> C2 + C3 fragments"""
    steps = []
    for s in pool_species:
        lower = s.label.lower()
        if "deoxyosone" in lower:
            # C6 -> C3 + C3
            if "glucose" in lower or s.smiles == "O=CC(=O)CC(O)C(O)CO":
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glyceraldehyde", smiles="O=CC(O)CO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
            # C5 -> C2 + C3
            elif "ribose" in lower or s.smiles == "O=CC(=O)CC(O)CO":
                p1 = Species(label="pyruvaldehyde", smiles="CC(=O)C=O")
                p2 = Species(label="glycolaldehyde", smiles="O=CCO")
                steps.append(ElementaryStep([s], [p1, p2], "Retro_Aldol_Fragmentation"))
    return steps


def _cysteine_degradation(pool_species: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """Thermal degradation of Cysteine to H2S, NH3, acetaldehyde, CO2."""
    steps = []
    if conditions.temperature_celsius < 100:
        return steps
    for s in pool_species:
        if "cysteine" in s.label.lower() or s.smiles == "NC(CS)C(=O)O":
            h2s = Species("H2S", "S")
            nh3 = Species("ammonia", "N")
            aca = Species("acetaldehyde", "CC=O")
            co2 = Species("CO2", "O=C=O")
            steps.append(ElementaryStep([s], [h2s, nh3, aca, co2], "Cysteine_Degradation"))
    return steps


def _thiazole_condensation(pool_species: List[Species]) -> List[ElementaryStep]:
    """Strecker aldehyde + H2S + NH3 -> Thiazole"""
    steps = []
    h2s = next((s for s in pool_species if s.smiles == "S"), None)
    nh3 = next((s for s in pool_species if s.smiles == "N"), None)
    if not (h2s and nh3):
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
                reactants=[sp, h2s, nh3], 
                products=[tz, water, water], 
                reaction_family="Lipid_Thiazole_Condensation"
            ))
    return steps


def _sugar_ring_opening(pool_species: List[Species]) -> List[ElementaryStep]:
    """Hemiacetal cyclic sugar -> open-chain aldehyde. (Defensive rule via RWMol)"""
    steps = []
    # Match hemiacetal: O(ring) - C(ring)(OH)
    # The [O;X2;R] ensures it's a ring oxygen, [C;X4;R] is a ring carbon, [O;X2;H] is the hydroxyl.
    patt = Chem.MolFromSmarts("[O;X2;R]-[C;X4;R]-[OH]")
    if not patt: return steps
    
    for s in pool_species:
        m = _mol(s.smiles)
        if not m: continue
            
        matches = m.GetSubstructMatches(patt)
        if not matches: continue
            
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


# ── PBMA Additive Degradations ───────────────────────────────────────────

def _thiamine_degradation(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: Thermal breakdown of Thiamine (Vitamin B1).
    Literature shows this yields H2S, 2-methylthiophene, and thiazoles.
    """
    if conditions.temperature_celsius < 100:
        return []
        
    steps = []
    # Identify Thiamine in the pool by canonical SMILES
    for s in pool:
        if _canonical(s.smiles) == _canonical(_THIAMINE_CANONICAL):
            p1 = Species("Hydrogen_Sulfide", "S")
            p2 = Species("2-methylthiophene", "Cc1cccs1")
            p3 = Species("4,5-dihydro-2-methylthiazole", "CC1=NCCS1")
            
            steps.append(ElementaryStep(
                reactants=[s],
                products=[p1, p2, p3],
                reaction_family="Additive_Thermal_Degradation"
            ))
    return steps


def _glutathione_cleavage(pool: List[Species], conditions: ReactionConditions) -> List[ElementaryStep]:
    """
    Tier B Template: Controlled cleavage of Glutathione (GSH).
    Literature shows it cleaves into glutamic acid and cysteinylglycine dipeptide.
    """
    if conditions.temperature_celsius < 100:
        return []
        
    steps = []
    for s in pool:
        if _canonical(s.smiles) == _canonical(_GSH_CANONICAL):
            p1 = Species("Glutamic_Acid", "N[C@@H](CCC(=O)O)C(=O)O")
            p2 = Species("Cysteinylglycine", "N[C@@H](CS)C(=O)NCC(=O)O")
            
            steps.append(ElementaryStep(
                reactants=[s],
                products=[p1, p2],
                reaction_family="Additive_Thermal_Degradation"
            ))
    return steps


# ──────────────────────────────────────────────────────────────────────────
# Tier A: SMIRKS application
# ──────────────────────────────────────────────────────────────────────────

def _apply_smirks_rule(
    name: str, family: str, smirks: str, ph_gate: str,
    pool: List[Species], conditions: ReactionConditions
) -> List[ElementaryStep]:
    """Apply a single SMIRKS rule to all relevant species pairs in the pool."""
    # pH gate check
    if ph_gate == "acid" and conditions.pH >= 6:
        return []
    if ph_gate == "neutral" and conditions.pH < 6:
        return []

    rxn = AllChem.ReactionFromSmarts(smirks)
    if rxn is None:
        return []

    steps = []
    pool_smiles = list({s.smiles for s in pool})

    n_reactants = rxn.GetNumReactantTemplates()

    if n_reactants == 1:
        for smi in pool_smiles:
            m = _mol(smi)
            if m is None:
                continue
            try:
                prods = rxn.RunReactants((m,))
            except Exception:
                continue
            for prod_tuple in prods:
                prod_smiles = []
                for p in prod_tuple:
                    try:
                        ps = Chem.MolToSmiles(p)
                        if _is_valid(ps):
                            prod_smiles.append(ps)
                    except Exception:
                        pass
                if prod_smiles:
                    reactant_sp = next((s for s in pool if s.smiles == smi), Species(smi, smi))
                    steps.append(ElementaryStep(
                        reactants=[reactant_sp],
                        products=[Species(ps, ps) for ps in prod_smiles],
                        reaction_family=family,
                    ))

    elif n_reactants == 2:
        for i, smi1 in enumerate(pool_smiles):
            for smi2 in pool_smiles:
                m1, m2 = _mol(smi1), _mol(smi2)
                if m1 is None or m2 is None:
                    continue
                try:
                    prods = rxn.RunReactants((m1, m2))
                except Exception:
                    continue
                for prod_tuple in prods:
                    prod_smiles = []
                    for p in prod_tuple:
                        try:
                            ps = Chem.MolToSmiles(p)
                            if _is_valid(ps):
                                prod_smiles.append(ps)
                        except Exception:
                            pass
                    if prod_smiles:
                        r1 = next((s for s in pool if s.smiles == smi1), Species(smi1, smi1))
                        r2 = next((s for s in pool if s.smiles == smi2), Species(smi2, smi2))
                        steps.append(ElementaryStep(
                            reactants=[r1, r2],
                            products=[Species(ps, ps) for ps in prod_smiles],
                            reaction_family=family,
                        ))
    return steps


# ──────────────────────────────────────────────────────────────────────────
# Main Engine
# ──────────────────────────────────────────────────────────────────────────

class SmirksEngine:
    """
    Hybrid rule-based Maillard pathway generator.

    Usage:
        engine = SmirksEngine(conditions=ReactionConditions(pH=5.5, temperature_celsius=150))
        steps = engine.enumerate([ribose, cysteine, glycine])
        # steps: List[ElementaryStep] — feed into xtb_screener or recommend.py
    """

    def __init__(self, conditions: ReactionConditions = None):
        self.conditions = conditions or ReactionConditions()

    def enumerate(
        self,
        precursors: List[Species],
        max_generations: int = 3
    ) -> List[ElementaryStep]:
        """
        Enumerate Maillard pathways from the given precursor pool.

        Phase 1: Apply Tier B templates (Amadori cascade, enolisation,
                 Strecker, beta-elimination).
        Phase 2: Apply Tier A SMIRKS rules iteratively on the growing pool.
        """
        all_steps: List[ElementaryStep] = []

        # Working pool: canonical SMILES → Species
        pool_dict: dict = {_canonical(s.smiles): s for s in precursors if _canonical(s.smiles)}

        def pool_list() -> List[Species]:
            return list(pool_dict.values())

        def add_to_pool(sp: Species):
            can = _canonical(sp.smiles)
            if can and can not in pool_dict:
                pool_dict[can] = sp

        def add_step_products(step: ElementaryStep):
            for p in step.products:
                if _is_valid(p.smiles):
                    add_to_pool(p)

        # ── Pre-Phase: Sugar Ring Opening ────────────────────────────────
        ring_steps = _sugar_ring_opening(pool_list())
        for step in ring_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Pre-Phase: PBMA Additive Degradations ─────────────────────────
        thiamine_steps = _thiamine_degradation(pool_list(), self.conditions)
        for step in thiamine_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)
                
        gsh_steps = _glutathione_cleavage(pool_list(), self.conditions)
        for step in gsh_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Tier B Phase 1: Amadori / Heyns cascade ──────────────────────
        sugars = [s for s in pool_list() if _is_sugar(s)]
        amines = [s for s in pool_list() if _is_primary_amine(s)]

        for sugar in sugars:
            for amine in amines:
                cascade = _amadori_cascade(sugar, amine)
                for step in cascade:
                    # Dedup by reactant+product labels
                    if not _step_exists(step, all_steps):
                        all_steps.append(step)
                        add_step_products(step)

                amadori_label = next((p.label for s in cascade for p in s.products if "Amadori" in p.label or "Heyns" in p.label), None)
                amadori_sp = next((s for s in pool_list() if s.label == amadori_label), None)
                if amadori_sp:
                    enols = _enolisation_steps(amadori_sp, sugar, self.conditions)
                    for enol in enols:
                        if not _step_exists(enol, all_steps):
                            all_steps.append(enol)
                            add_step_products(enol)

        # ── Tier B Phase 2: Retro-aldol and Strecker degradation ─────────
        ra_steps = _retro_aldol_fragmentation(pool_list())
        for step in ra_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        dicarbonyls = [s for s in pool_list() if _is_dicarbonyl(s)]
        amines_now = [s for s in pool_list() if _is_primary_amine(s)]

        for dc in dicarbonyls:
            for amine in amines_now:
                s_step = _strecker_step(dc, amine)
                if s_step and not _step_exists(s_step, all_steps):
                    all_steps.append(s_step)
                    add_step_products(s_step)

        # ── Tier B Phase 3: Secondary Condensations & Eliminations ────────
        # 3a. Beta-elimination (DHA pathway)
        beta_candidates = [s for s in pool_list() if _has_cysteine_beta_carbon(s)]
        for aa in beta_candidates:
            be_steps = _beta_elimination_steps(aa, pool_list())
            for step in be_steps:
                if not _step_exists(step, all_steps):
                    all_steps.append(step)
                    add_step_products(step)

        # 3b. Cysteine thermal degradation
        cys_steps = _cysteine_degradation(pool_list(), self.conditions)
        for step in cys_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3c. Aminoketone Condensation (Pyrazines)
        ak_steps = _aminoketone_condensation(pool_list())
        for step in ak_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # 3d. Lipid Thiazole Condensation
        tz_steps = _thiazole_condensation(pool_list())
        for step in tz_steps:
            if not _step_exists(step, all_steps):
                all_steps.append(step)
                add_step_products(step)

        # ── Tier A: SMIRKS rules, iterative ──────────────────────────────
        seen_step_keys: Set[str] = {_step_key(s) for s in all_steps}

        for _gen in range(max_generations):
            new_steps_this_gen = []
            current_pool = pool_list()

            for name, family, smirks, gate in _SMIRKS_RULES:
                candidates = _apply_smirks_rule(
                    name, family, smirks, gate, current_pool, self.conditions
                )
                for step in candidates:
                    k = _step_key(step)
                    if k not in seen_step_keys:
                        new_steps_this_gen.append(step)
                        seen_step_keys.add(k)

            if not new_steps_this_gen:
                break  # No new reactions found — converged

            for step in new_steps_this_gen:
                add_step_products(step)
            all_steps.extend(new_steps_this_gen)

        return all_steps


def _step_key(step: ElementaryStep) -> str:
    """Stable hash key for deduplication: sorted reactants + sorted products."""
    reacts = tuple(sorted(_canonical(r.smiles) or r.label for r in step.reactants))
    prods = tuple(sorted(_canonical(p.smiles) or p.label for p in step.products))
    return str((reacts, prods))


def _step_exists(step: ElementaryStep, existing: List[ElementaryStep]) -> bool:
    key = _step_key(step)
    return any(_step_key(s) == key for s in existing)


# ──────────────────────────────────────────────────────────────────────────
# CLI demo
# ──────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from data.reactions.curated_pathways import (
        RIBOSE, GLUCOSE, GLYCINE, CYSTEINE, LEUCINE, LYSINE, HEXANAL, H2S
    )

    SYSTEMS = [
        ("Ribose + Glycine @ pH 5",  [RIBOSE, GLYCINE],            ReactionConditions(pH=5.0)),
        ("Glucose + Glycine @ pH 7", [GLUCOSE, GLYCINE],           ReactionConditions(pH=7.0)),
        ("Ribose + Cysteine @ pH 6", [RIBOSE, CYSTEINE, H2S],      ReactionConditions(pH=6.0)),
        ("Ribose + Cys + Leu",       [RIBOSE, CYSTEINE, LEUCINE, H2S], ReactionConditions(pH=6.0)),
    ]

    for label, precursors, conds in SYSTEMS:
        print(f"\n{'='*60}")
        print(f"System: {label}")
        print(f"Input:  {[p.label for p in precursors]}")
        engine = SmirksEngine(conditions=conds)
        steps = engine.enumerate(precursors)
        print(f"Generated {len(steps)} elementary steps:")
        for step in steps:
            print(f"  {step}")
