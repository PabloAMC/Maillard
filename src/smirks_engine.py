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

from typing import List, Optional, Set, Tuple
from functools import lru_cache

# Suppress RDKit atom-mapping warnings from SMIRKS rules that use unmapped atoms
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402

from src.pathway_extractor import Species, ElementaryStep  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from src.chem_utils import canonicalize_smiles as _canonical  # noqa: E402

# Suppress RDKit atom-mapping warnings
RDLogger.DisableLog("rdApp.warning")

# ──────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────

MAX_MW = 300.0  # Daltons — prune products above this (volatiles are small)
# AUDIT 2026-08-27 (Wave H) — NOTE ONLY, deliberately NOT changed in this wave.
# This cap is right for the Maillard trunk (every volatile of interest is under
# 200 Da) and WRONG for the lipid lane, where it silently truncates the radical
# chain Wave G1 had just repaired. Measured, at this exact threshold:
#     linoleic acid            280.5   kept
#     linoleyl alkoxy radical  295.4   kept
#     linoleyl peroxyl radical 311.4   PRUNED
#     13-HPODE                 312.5   PRUNED
#     9-HPODE                  326.5   PRUNED
#     arachidonic acid         304.5   PRUNED
# So `Radical_Propagation_O2` (R. + O2 -> ROO.) and `Peroxy_H_Abstraction`
# (ROO. + RH -> ROOH) cannot deposit their C18 products in the pool: the chain
# propagates only as far as the alkoxy radical, and a hydroperoxide seed supplied
# as a PRECURSOR survives (precursors are not pruned) while the same molecule
# made by the network does not. It also means arachidonic acid cannot be used as
# a precursor product at all. Hexanal is unaffected (100 Da) because it arrives by
# beta-scission of the alkoxy radical, which stays under the cap.
# Raising the cap is a real change to network size and cost and belongs to a
# dedicated lipid workstream with its own regeneration — see
# tasks/audit_remediation.md. It is recorded here so the limit stops being invisible.

# Additive Canonical SMILES (for exact matching)
_THIAMINE_CANONICAL = "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"
_GSH_CANONICAL = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"

# ── Verified reference structures (AUDIT 2026-08-27, Wave G1 fix 4) ───────
# The previous _DMHF_/_HEMF_ strings placed the enol OH on the ring-oxygen
# carbon, which is a different heterocycle from the furaneol family they were
# labelled as.  Both are now the correct 4-hydroxy-3(2H)-furanones, checked by
# InChIKey (see tests/unit/test_chemistry_soundness.py).
_DMHF_CANONICAL = "CC1OC(C)=C(O)C1=O"       # furaneol, C6H8O3,  INAXVXBDKKUCGI
_HEMF_CANONICAL = "CCC1OC(C)=C(O)C1=O"      # HEMF,     C7H10O3, GWCRPYGYVRXVLI
_NORFURANEOL_CANONICAL = "CC1=C(O)C(=O)CO1"  # 4-hydroxy-5-methyl-3(2H)-furanone
#                                              (norfuraneol), C5H6O3, DLVYTANECMRFGX
_MFT_CANONICAL = "Cc1c(S)cco1"              # 2-methyl-3-furanthiol, RUYNUXHHUVUINQ
# bis(2-methyl-3-furyl) disulfide.  The template layer previously carried the
# mixed 3-furyl/4-furyl regioisomer `Cc1cc(SSC2=C(C)OC=C2)co1`
# (JQDSBCSEZSRHAY); the curated layer already had the correct molecule.
_FURYL_DISULFIDE_CANONICAL = "Cc1c(SSc2ccoc2C)cco1"  # C10H10O2S2, OHDFENKFSKIFBJ

# 1-deoxyosones (2,3-enolisation products) — the real MFT/furanone precursors.
_DEOXYOSONE_1_PENTOSE = "CC(=O)C(=O)C(O)CO"      # 1-deoxy-2,3-pentodiulose, C5H8O4
_DEOXYOSONE_1_HEXOSE = "CC(=O)C(=O)C(O)C(O)CO"   # 1-deoxy-2,3-hexodiulose,  C6H10O5
# 1,4-dideoxypento-2,3-diulose (CH3-CO-CO-CH2-CH2OH, C5H8O3) — the isotope-
# evidenced in-situ MFT precursor (Cerny & Davidek 2003, 10.1021/jf026123f;
# positionally confirmed with [1-13C]ribose in Cerny & Davidek 2004,
# 10.1021/jf035265m). Added 2026-08-27 (Wave N) when the norfuraneol->MFT step
# was retired; see `_furanone_and_mft_route` in reaction_templates.py.
_PENTODIULOSE_14_DIDEOXY = "CC(=O)C(=O)CCO"      # 1,4-dideoxypento-2,3-diulose, C5H8O3

# ── Wave P (2026-08-27): the C2 + C3 recombination lane to MFT ────────────
# Hofmann & Schieberle 1998 (10.1021/jf9705983) got their HIGHEST MFT yield of
# all — 1.4 mol %, 6 min, 180 C, in the absence of water — from hydroxyacetaldehyde
# (= glycolaldehyde, already a species) plus MERCAPTO-2-PROPANONE, a topology the
# network could not express because the C3 partner did not exist.  Cerny 2015
# (10.1016/b978-1-78242-103-0.00009-6, full text) renders Hofmann & Schieberle's
# own scheme: "The formation of 2-methyl-3-furanthiol (16) starts with aldol
# reaction of hydroxyacetaldehyde (18) and mercaptopropanone (19) to give
# 4,5-dihydroxy-3-mercapto-2-pentanone (20). The intermediate 20 cyclises and
# dehydrates to yield 2-methyl-3-furanthiol (16)."  Both intermediates below are
# NAMED in that scheme; neither is invented here.
_MERCAPTO_2_PROPANONE = "CC(=O)CS"               # 1-mercapto-2-propanone, C3H6OS
_DIHYDROXY_MERCAPTO_PENTANONE = "CC(=O)C(S)C(O)CO"
#                                                # 4,5-dihydroxy-3-mercapto-2-pentanone,
#                                                # C5H10O3S — Cerny 2015 intermediate 20

# ── Wave P (2026-08-27): the norfuraneol -> 2-mercapto-3-pentanone lane ───
# Cerny & Davidek 2003 (10.1021/jf026123f) showed norfuraneol is NOT the MFT
# intermediate but IS the demonstrated precursor of 2-mercapto-3-pentanone:
# "Whereas 2-mercapto-3-pentanone was found unlabeled and hence originated from
# 4-hydroxy-5-methyl-3(2H)-furanone, its isomer 3-mercapto-2-pentanone was formed
# from both 4-hydroxy-5-methyl-3(2H)-furanone and ribose."  That paper proposes no
# mechanism; Whitfield & Mottram 1999 (10.1021/jf980980v) supply it, reporting
# 2,3-PENTANEDIONE as one of the "main non-sulfur compounds" of exactly this
# system (norfuraneol + cysteine or H2S, pH 4.5, 140 C, 60 min).
_PENTANE_2_3_DIONE = "CCC(=O)C(C)=O"             # 2,3-pentanedione, C5H8O2
_MERCAPTO_3_PENTANONE_2 = "CCC(=O)C(C)S"         # 2-mercapto-3-pentanone, C5H10OS

# Isotope tag used inside SMIRKS PRODUCT templates to mark the atom that must
# survive sanitisation as an open-shell radical centre.  Without it RDKit's
# sanitiser fills the unsatisfied valence with an INVENTED implicit hydrogen
# (ROO. silently becomes ROOH).  `_apply_smirks_rule` strips the tag and pins
# NoImplicit + one radical electron BEFORE sanitising.  Mirrors the isotope-13
# trick already used in `_lipid_hydroperoxide_scission`.
RADICAL_TAG_ISOTOPE = 99

# Tier A SMIRKS rules: (name, reaction_family, smirks, ph_gate)
# ph_gate: "any" | "acid" (pH<6) | "neutral" (pH>=6)
_SMIRKS_RULES: List[Tuple[str, str, str, str]] = [
    (
        "schiff_base_lipid",
        "Lipid_Schiff_Base",
        # C3+ aliphatic aldehyde whose alpha-carbon has NO hydroxyl (excludes sugars).
        # The amine donor can be anything with a primary amine on an sp3 carbon (like amino acids).
        #
        # 2026-08-27: `!$(C[#7])` excludes alpha-AMINO aldehydes.  Those are
        # bifunctional (aldehyde AND primary amine), so every condensation
        # product is itself a substrate for the next generation and the rule
        # ran away into amino-aldehyde oligomers that are not volatiles.  This
        # became load-bearing once glyoxal entered the pool (its Strecker
        # partner is 2-aminoethanal).
        #
        # AUDIT 2026-08-27 (Wave I fix 9) — red-team finding H3.  The pattern
        # was far broader than the comment two lines above it claimed, and had
        # been for as long as the comment existed:
        #   * "alpha-carbon has NO hydroxyl (excludes sugars)" was NOT
        #     implemented.  `[CX4H2]` matches the alpha carbon of an
        #     alpha-hydroxy aldehyde perfectly well — glycolaldehyde (OCC=O),
        #     glyceraldehyde, and every sugar-derived alpha-hydroxy carbonyl in
        #     the pool were legitimate substrates for a rule named
        #     `Lipid_Schiff_Base`.
        #   * "C3+" was NOT implemented either: `CX4H3` is a methyl, so
        #     acetaldehyde (C2, a Strecker product, not a lipid aldehyde) matched.
        # The two exclusions the comment already claimed are now actually in the
        # SMARTS:
        #   `!$(C[OX2])`      alpha-hydroxy/alkoxy carbonyls are out (no sugars);
        #   `$(C([#6])[#6])`  the alpha carbon carries a second carbon besides
        #                     the carbonyl, i.e. the aldehyde is C3 or longer;
        #   `CX4H3` dropped   which is the C2 (acetaldehyde) case.
        # MEASURED, core no-lipid network (D-Glucose + D-Ribose + Glycine +
        # L-Cysteine + L-Lysine + L-Asparagine, 150 C / pH 5.5 / aw 0.95,
        # max_generations=4; reproduced by
        # tests/scientific/test_wave_i_network_chemistry.py):
        #   before: 228 Lipid_Schiff_Base steps of 303 total = 75.2%
        #   after:   28 Lipid_Schiff_Base steps of 103 total = 27.2%
        # (Wave I fix 12 later removed one unrelated step from the same network,
        # so the end-of-wave figure the regression test pins is 28 of 102 = 27.5%.
        # The 303-step "before" is itself post-fix-8; the pre-Wave-I baseline was
        # 224 of 298 = 75.2%, i.e. the same share.)
        # Genuine lipid aldehydes are unaffected: hexanal (CCCCCC=O), nonanal
        # and propanal still match; acetaldehyde (C2) and glycolaldehyde
        # (alpha-hydroxy) no longer do.  The same test pins both directions.
        #
        # RESIDUAL DEFECT, NOT FIXED HERE (stated so it is not mistaken for a
        # clean result).  All 28 surviving steps come from ONE substrate,
        # 5-aminopentanal (NCCCCC=O), the Strecker aldehyde of lysine.  It is
        # bifunctional in exactly the way the 2026-08-27 note above describes,
        # but its amine sits on the OMEGA carbon, not the alpha carbon, so
        # `!$(C[#7])` — which is a constraint on the alpha carbon only — does
        # not see it.  Each condensation product still carries a free aldehyde
        # and a free amine, so the rule oligomerises: 12 first-generation steps
        # become 28 over four generations.  Killing that needs a whole-molecule
        # "no primary amine anywhere on the aldehyde" condition, which is not
        # expressible in this per-atom SMARTS and would need a rule gate; it is
        # an OPEN OWNER ITEM, deliberately out of scope for Wave I fix 9, whose
        # brief was to implement the two exclusions the comment already claimed.
        "[CX4H2;!$(C[#7]);!$(C[OX2]);$(C([#6])[#6]):2][CH1:1]=[O:5].[NH2:3][CX4:4]"
        ">>[*:2][CH1:1]=[N:3][*:4].[O:5]",
        "any",
    ),
    (
        "beta_scission_alkoxy",
        "Beta_Scission",
        # Alkoxy radical (R-O.) beta-scission: R-C(O.)(H)-R' -> R-CHO + .R'
        #
        # AUDIT 2026-08-27 (Wave G1 fix 3): the beta carbon was pinned to
        # [CX4:3], i.e. sp3.  The 13-alkoxy radical of linoleate is
        # C5H11-CH(O.)-CH=CH-CH=CH-(CH2)7-COOH, so the bond whose scission
        # releases HEXANAL (C13-C18 fragment) is C12-C13 with an sp2 beta
        # carbon — structurally unmatchable, which made hexanal UNREACHABLE
        # in the lipid network.  The beta carbon is now [C;X3,X4], covering
        # both the sp3 (pentyl-radical) and the allylic/vinylic (hexanal)
        # channels.  The leaving fragment is tagged so it stays a radical.
        "[C;X3,X4:3][CX4:1][OX1H0;v1:2]>>[C:1]=[O:2].[99C:3]",
        "any",
    ),
    (
        "radical_propagation_o2",
        "Radical_Propagation_O2",
        # Alkyl radical + O2 -> peroxy radical.
        #
        # AUDIT 2026-08-27 (Wave G1 fix 1): the carbon was [C;X3:1], which
        # matches ANY sp2 carbon — every aldehyde, ketone, imine and alkene in
        # the pool was being "peroxidised", fabricating the majority of the
        # lipid network.  [#6;v3;+0] is total-valence 3 neutral carbon, i.e.
        # an actual carbon radical (a C=O or C=C carbon has valence 4).
        # O2 is accepted in either of its two RDKit spellings (O=O / [O][O]).
        "[#6;v3;+0:1].[O;X1;+0:2]-,=[O;X1;+0:3]>>[#6:1]-[O:2]-[99O:3]",
        "any",
    ),
    (
        "peroxy_h_abstraction",
        "Peroxy_H_Abstraction",
        # Peroxy radical + allylic C-H -> hydroperoxide + allylic carbon radical.
        # [O;X1;v1]-[O;X2] is a genuine ROO. (a hydroperoxide's terminal O is X2).
        # The abstracting carbon is tagged so the product stays a radical.
        # `X4` keeps the abstraction on genuinely allylic sp3 C-H; without it the
        # rule also stripped VINYLIC hydrogens, seeding vinyl radicals that then
        # fed the propagation rule with fabricated peroxides.
        "[O;X1;v1;+0:1]-[O;X2:2].[C;X4;H1,H2,H3;$(C-C=C);!$(C=O):3]>>[O;X2:1]([H])[O;X2:2].[99C:3]",
        "any",
    ),
    (
        "radical_termination_russell",
        "Radical_Termination",
        # AUDIT 2026-08-27 (Wave G1 fix 2): the old rule emitted
        # `O2 + 2 ROH`, i.e. it invented two hydrogens out of nothing.
        # Replaced by the Russell mechanism, the accepted bimolecular
        # termination of secondary/primary peroxy radicals:
        #     R2CH-OO. + R'2CH-OO.  ->  R2CH-OH + R'2C=O + O2
        # The H that closes the alcohol is the one lost by the carbon that
        # becomes the carbonyl, so the step is balanced with no free H.
        "[C;!H0:5]-[O;X2:2]-[O;X1;v1;+0:1].[C:6]-[O;X2:4]-[O;X1;v1;+0:3]"
        ">>[C:5]=[O:2].[C:6]-[O:4][H].[O:1]=[O:3]",
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
_DICARBONYL_SMARTS = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)")  # adjacent carbonyls
_AROMATIC_ALDEHYDE_SMARTS = Chem.MolFromSmarts("[c][CH]=O") # furfural-type


@lru_cache(maxsize=4096)
def _mol_cached(smi: str) -> Optional[Chem.Mol]:
    """Internal cached Mol parsing."""
    return Chem.MolFromSmiles(smi) if smi else None


def _mol(smi: str) -> Optional[Chem.Mol]:
    """Returns a CLONE of the cached Mol to prevent in-place mutation issues."""
    m = _mol_cached(smi)
    return Chem.Mol(m) if m else None





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
    if not m: 
        return False
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


def _is_asparagine(s: Species) -> bool:
    """Detects strictly free asparagine."""
    if s.label.lower() in ["l-asparagine", "asparagine"]: 
        return True
    return s.smiles == "NC(CC(N)=O)C(=O)O"


def _is_lysine(s: Species) -> bool:
    """Detects strictly free lysine."""
    if s.label.lower() in ["l-lysine", "lysine"]: 
        return True
    return s.smiles == "NCCCCC(N)C(=O)O"

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


def _is_lipid_aldehyde(s: Species) -> bool:
    """C5+ aliphatic monocarbonyl aldehyde without multiple OH (excludes sugars/dicarbonyls)."""
    m = _mol(s.smiles)
    if not m: 
        return False
    # Has aldehyde
    if not m.HasSubstructMatch(_ALDEHYDE_SMARTS): 
        return False
    # NOT aromatic
    if m.HasSubstructMatch(Chem.MolFromSmarts("a")): 
        return False
    # NOT dicarbonyl
    if m.HasSubstructMatch(_ALDEHYDE_SMARTS) and len(m.GetSubstructMatches(_ALDEHYDE_SMARTS)) > 1:
        return False
    # NOT nitrogenous (excludes amino-aldehydes like 5-aminopentanal)
    if any(atom.GetAtomicNum() == 7 for atom in m.GetAtoms()): 
        return False
    # Not a sugar (oh_count < 2)
    oh_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1)
    # C5+ (typically lipid-derived volatiles like pentanal, hexanal)
    c_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    return oh_count < 2 and c_count >= 5


def _is_lipid_hydroperoxide(s: Species) -> bool:
    """Heuristic: has a hydroperoxide group [OX2H,OX2-] and a lipid chain."""
    m = _mol(s.smiles)
    if not m:
        return False
    # Match R-O-OH or anion
    pat = Chem.MolFromSmarts("[OX2,OX1-][OX2H,OX1H0-]")
    if not m.HasSubstructMatch(pat):
        return False
    # C8+ (typical PUFA derivatives)
    c_count = sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() == 6)
    return c_count >= 8



from src.reaction_templates import (
    _amadori_cascade, _enolisation_steps, _strecker_step, _beta_elimination_steps,
    _aminoketone_condensation, _retro_aldol_fragmentation, _cysteine_degradation,
    _thiazole_condensation, _thiol_addition, _mft_pathway, _furyl_disulfide_formation,
    _sulfur_volatiles_pathway, _deamination_step, _lipid_maillard_synergy,
    _lipid_hydroperoxide_scission, _sugar_ring_opening,
    _acrylamide_formation, _cml_cel_formation, _thiamine_degradation, _furanone_generation,
    _glutathione_cleavage, _furanone_and_mft_route, _thiol_reductant_pool,
    _c2_c3_mft_recombination, _norfuraneol_mercaptopentanone_route,
)
# ──────────────────────────────────────────────────────────────────────────
# Tier A: SMIRKS application
# ──────────────────────────────────────────────────────────────────────────

def _fix_radicals(mol: Chem.Mol, family: str, clear_all: bool = False):
    """Ensure atoms with unsatisfied valence are marked as radicals, and clear satisfied ones."""
    if clear_all:
        for atom in mol.GetAtoms():
            atom.SetNumRadicalElectrons(0)
        return mol

    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        # Use explicit valence (bond order sum) + implicit H
        # This is more accurate for double/triple bonds.
        try:
            val = atom.GetExplicitValence() + atom.GetNumImplicitHs()
        except Exception:
            # Fallback if property cache is not updated
            val = atom.GetTotalDegree() + atom.GetTotalNumHs()
        
        # Determine if we should set or clear radicals
        if an == 6: # Carbon
            if val == 3:
                atom.SetNumRadicalElectrons(1)
            elif val >= 4:
                atom.SetNumRadicalElectrons(0)
        elif an == 8: # Oxygen
            if val == 1:
                if any(f in family for f in ["Radical_Propagation", "Lipid_Homolysis", "Beta_Scission"]):
                    atom.SetNumRadicalElectrons(1)
            elif val >= 2:
                atom.SetNumRadicalElectrons(0)
    return mol


def _pin_tagged_radicals(mol: Chem.Mol) -> None:
    """Convert isotope-tagged product atoms into pinned radical centres.

    MUST run BEFORE ``Chem.SanitizeMol``.  Sanitisation is what fills an
    unsatisfied valence with implicit hydrogens, so a peroxy oxygen that is
    only marked afterwards has already become a hydroperoxide with an invented
    H (AUDIT 2026-08-27, Wave G1 fix 2).  ``SetNoImplicit(True)`` freezes the
    hydrogen count one below the closed-shell fill, so the deficit is carried
    as a radical electron instead of an invented hydrogen.
    """
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != RADICAL_TAG_ISOTOPE:
            continue
        atom.SetIsotope(0)
        # The tagged atom either lost a bond (beta-scission, homolysis) or lost
        # an H (abstraction).  In both cases RDKit's default implicit-H fill
        # would close the shell; the radical is exactly one H short of that.
        radical_hs = max(atom.GetTotalNumHs() - 1, 0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(radical_hs)
        atom.SetNumRadicalElectrons(1)


def _finalize_smirks_product(mol: Chem.Mol, family: str) -> Optional[str]:
    """Sanitise one SMIRKS product and return its SMILES, or None on failure.

    Order is load-bearing: tagged radical centres are pinned first, then the
    molecule is sanitised, then ``_fix_radicals`` is given a last pass to catch
    open valences that the rule did not tag explicitly (e.g. fragments
    inherited from the reactant).  The SMILES is emitted LAST so that the
    string stored on the ElementaryStep is the radical one — previously the
    1-reactant branch serialised the product BEFORE ``_fix_radicals`` ran, so
    every radical fragment was recorded as its closed-shell analogue.
    """
    try:
        _fix_radicals(mol, family, clear_all=True)
        _pin_tagged_radicals(mol)
        Chem.SanitizeMol(mol)
        # Product templates that spell an added hydrogen as `[H]` leave it as a
        # separate heavy-atom node; fold it back into the implicit count so the
        # emitted SMILES canonicalises identically to the same species arriving
        # from a template.
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass
        _fix_radicals(mol, family)
        mol.UpdatePropertyCache(strict=False)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

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
                # AUDIT 2026-08-27 (Wave G1 fix 2): the SMILES used to be
                # serialised here, BEFORE the radical flags were repaired, so
                # every radical fragment this branch produced (notably the
                # beta-scission alkyl fragment) was recorded closed-shell.
                prod_smiles = []
                valid_step = True
                for p in prod_tuple:
                    ps = _finalize_smirks_product(p, family)
                    if ps is not None and _is_valid(ps):
                        prod_smiles.append(ps)
                    else:
                        valid_step = False
                        break
                if valid_step and len(prod_smiles) == len(prod_tuple):
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
                    valid_step = True
                    for p in prod_tuple:
                        # AUDIT 2026-08-27 (Wave G1 fix 2): sanitisation used to
                        # run before the radical flags were re-checked, which
                        # handed the terminal peroxy oxygen an invented implicit
                        # H (ROO. -> ROOH).  `_finalize_smirks_product` pins the
                        # tagged radical centres first.
                        ps = _finalize_smirks_product(p, family)
                        if ps is not None and _is_valid(ps):
                            prod_smiles.append(ps)
                        else:
                            valid_step = False
                            break

                    # Only append if ALL products were successfully generated and are valid
                    # This guarantees mass conservation. The RDKit reaction MUST output everything.
                    if valid_step and len(prod_smiles) == len(prod_tuple):
                        r1 = next((s for s in pool if s.smiles == smi1), Species(smi1, smi1))
                        r2 = next((s for s in pool if s.smiles == smi2), Species(smi2, smi2))
                        steps.append(ElementaryStep(
                            reactants=[r1, r2],
                            products=[Species(ps, ps) for ps in prod_smiles],
                            reaction_family=family,
                        ))

    elif n_reactants == 3:
        for smi1 in pool_smiles:
            for smi2 in pool_smiles:
                for smi3 in pool_smiles:
                    m1, m2, m3 = _mol(smi1), _mol(smi2), _mol(smi3)
                    if m1 is None or m2 is None or m3 is None:
                        continue
                    try:
                        prods = rxn.RunReactants((m1, m2, m3))
                    except Exception:
                        continue
                    for prod_tuple in prods:
                        prod_smiles = []
                        valid_step = True
                        for p in prod_tuple:
                            ps = _finalize_smirks_product(p, family)
                            if ps is not None and _is_valid(ps):
                                prod_smiles.append(ps)
                            else:
                                valid_step = False
                                break

                        if valid_step and len(prod_smiles) == len(prod_tuple):
                            r1 = next((s for s in pool if s.smiles == smi1), Species(smi1, smi1))
                            r2 = next((s for s in pool if s.smiles == smi2), Species(smi2, smi2))
                            r3 = next((s for s in pool if s.smiles == smi3), Species(smi3, smi3))
                            steps.append(ElementaryStep(
                                reactants=[r1, r2, r3],
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

    def __init__(self, conditions: Optional[ReactionConditions] = None):
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

        def _add_steps(steps_to_add: List[ElementaryStep]):
            for s in steps_to_add:
                if not _step_exists(s, all_steps):
                    all_steps.append(s)
                    add_step_products(s)

        # ── Pre-Phase: Sugar Ring Opening ────────────────────────────────
        ring_steps = _sugar_ring_opening(pool_list())
        _add_steps(ring_steps)

        # ── Pre-Phase: PBMA Additive Degradations ─────────────────────────
        thiamine_steps = _thiamine_degradation(pool_list(), self.conditions)
        _add_steps(thiamine_steps)
                
        gsh_steps = _glutathione_cleavage(pool_list(), self.conditions)
        _add_steps(gsh_steps)

        furanone_steps = _furanone_generation(pool_list(), self.conditions)
        _add_steps(furanone_steps)

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

                amadori_smiles = next((p.smiles for s in cascade for p in s.products if "Amadori" in p.label or "Heyns" in p.label), None)
                amadori_can = _canonical(amadori_smiles) if amadori_smiles else None
                amadori_sp = pool_dict.get(amadori_can) if amadori_can else None
                
                if amadori_sp:
                    enols = _enolisation_steps(amadori_sp, sugar, amine, self.conditions)
                    for enol in enols:
                        if not _step_exists(enol, all_steps):
                            all_steps.append(enol)
                            add_step_products(enol)

        # ── Tier B Phase 2: Retro-aldol and Strecker degradation ─────────
        ra_steps = _retro_aldol_fragmentation(pool_list())
        _add_steps(ra_steps)

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
            _add_steps(be_steps)

        # 3b. Cysteine thermal degradation
        cys_steps = _cysteine_degradation(pool_list(), self.conditions)
        _add_steps(cys_steps)

        # 3b-1. Reducing-equivalent pool (Wave I fix 8, 2026-08-27; red-team H4).
        # `2 cysteine -> cystine + 2[H]`. Must run BEFORE every lane that is
        # pool-gated on the `[HH]` token (3e, 3e-1, 3e-2), otherwise the token's
        # only producer reachable from a cysteine/sugar system is the pyrazine
        # aromatisation at 3c and MFT becomes a downstream dependent of pyrazine
        # chemistry. See `_thiol_reductant_pool`.
        _add_steps(_thiol_reductant_pool(pool_list()))

        # 3c. Aminoketone Condensation (Pyrazines)
        ak_steps = _aminoketone_condensation(pool_list())
        _add_steps(ak_steps)

        # 3d. Lipid Thiazole Condensation
        tz_steps = _thiazole_condensation(pool_list())
        _add_steps(tz_steps)

        # 3e-0. R.12: Deamination (Must happen before volatile templates)
        deam_steps = _deamination_step(pool_list())
        _add_steps(deam_steps)

        # 3e. Thiol Addition (Furfural + H2S + H2 -> FFT)
        ta_steps = _thiol_addition(pool_list())
        _add_steps(ta_steps)

        # 3e-1. PRIMARY MFT route (Wave G1 fix 7; route corrected Wave N,
        # 2026-08-27 on isotope evidence — Cerny & Davidek 2003
        # 10.1021/jf026123f, 2004 10.1021/jf035265m):
        #   Amadori -(2,3-enolisation)-> 1-deoxyosone      [in _enolisation_steps]
        #          -(+2[H], C4 deoxygenation)-> 1,4-dideoxypento-2,3-diulose
        #          -(+ H2S)-> 2-methyl-3-furanthiol
        # Norfuraneol is still produced (cyclodehydration branch;
        # van den Ouweland & Peer 1975, 10.1021/jf60199a045, is its genuine
        # SYNTHESIS route to MFT) but no longer feeds MFT in situ; see
        # `_furanone_and_mft_route`'s docstring.
        furanone_steps = _furanone_and_mft_route(pool_list())
        _add_steps(furanone_steps)
        # Second pass so the +H2S step can see a 1,4-dideoxyosone that only
        # entered the pool on the first pass.
        _add_steps(_furanone_and_mft_route(pool_list()))

        # 3e-1b. SECOND MFT channel (Wave P item 2, 2026-08-27): the C2 + C3
        # recombination that Hofmann & Schieberle 1998 (10.1021/jf9705983) measured
        # as their HIGHEST-yielding MFT system (1.4 mol %, 6 min, 180 C, dry) and
        # that the network could not express at all until this wave, because
        # mercapto-2-propanone was not a species:
        #   pyruvaldehyde + H2S + 2[H] -> 1-mercapto-2-propanone + H2O
        #   glycolaldehyde + 1-mercapto-2-propanone
        #                            -> 4,5-dihydroxy-3-mercapto-2-pentanone
        #                            -> MFT + 2 H2O
        # Runs AFTER 3e-1 so that MFT from both channels lands in one pool species,
        # and after 3e-0/Phase 2 so the C2/C3 fragments exist. Deliberately NOT
        # moisture-gated; see the function docstring for why, and for the measured
        # reachability limit (glycolaldehyde is emitted only by the pentose
        # retro-aldol channel).
        _add_steps(_c2_c3_mft_recombination(pool_list()))
        # Second pass so the aldol step can see a mercaptopropanone that only
        # entered the pool on the first pass.
        _add_steps(_c2_c3_mft_recombination(pool_list()))

        # 3e-1c. NORFURANEOL's real sulfur fate (Wave P item 3, 2026-08-27). Wave N
        # retired norfuraneol from the MFT lane and left it with ZERO consumers;
        # Cerny & Davidek 2003 (10.1021/jf026123f) demonstrated what it does
        # instead, and Whitfield & Mottram 1999 (10.1021/jf980980v) supply the
        # 2,3-pentanedione intermediate:
        #   norfuraneol + 2 x 2[H] -> 2,3-pentanedione + H2O
        #   2,3-pentanedione + H2S + 2[H] -> 2-mercapto-3-pentanone + H2O
        _add_steps(_norfuraneol_mercaptopentanone_route(pool_list()))
        _add_steps(_norfuraneol_mercaptopentanone_route(pool_list()))

        # 3e-2. DEMOTED one-step 3-deoxyosone shortcut, kept only for hexose
        #       reachability; see `_mft_pathway`'s docstring.
        mft_steps = _mft_pathway(pool_list())
        _add_steps(mft_steps)

        # 3e-2b. Furyl disulfide formation from MFT
        furyl_disulfide_steps = _furyl_disulfide_formation(pool_list())
        _add_steps(furyl_disulfide_steps)

        # 3e-3. Methionine Sulfur Volatiles (Phase R.2 Fix)
        sulf_steps = _sulfur_volatiles_pathway(pool_list())
        _add_steps(sulf_steps)

        # 3f. Safety / Toxic Markers (Acrylamide, CML, CEL)
        for sug in sugars:
            for amine in amines_now:
                acry_steps = _acrylamide_formation(sug, amine)
                _add_steps(acry_steps)
        
        for amine in amines_now:
            age_steps = _cml_cel_formation(amine, pool_list())
            _add_steps(age_steps)

        # 3g. Lipid-Maillard Synergy (Lipid Aldehyde + Strecker AK)
        syn_steps = _lipid_maillard_synergy(pool_list())
        _add_steps(syn_steps)

        # 3h. Lipid Oxidation Radicals (Phase 19)
        hooh_steps = _lipid_hydroperoxide_scission(pool_list())
        _add_steps(hooh_steps)

        # `_radical_crosstalk_templates` was retired here on 2026-08-27
        # (Wave G1 fix 7): it existed to consume the fictitious elemental
        # sulfur ejected by the old thiol steps, and its second branch
        # quenched radicals by CONSUMING MFT.  See the retirement note in
        # src/reaction_templates.py.

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

