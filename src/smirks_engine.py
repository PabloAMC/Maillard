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
from typing import Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Suppress RDKit atom-mapping warnings from SMIRKS rules that use unmapped atoms
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import AllChem, Descriptors  # noqa: E402
from rdkit.Chem.rdChemReactions import ChemicalReaction  # noqa: E402

from src.pathway_extractor import Species, ElementaryStep  # noqa: E402
from src.conditions import ReactionConditions  # noqa: E402
from src.chem_utils import canonicalize_smiles as _canonical, parse_mol as _mol, calculate_mw as _mw  # noqa: E402
from src.sugar_classifier import is_ketose, is_pentose, is_hexose, is_sugar

# Suppress RDKit atom-mapping warnings
RDLogger.DisableLog("rdApp.warning")

# ──────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────

MAX_MW = 300.0  # Daltons — prune products above this (volatiles are small)

# Additive Canonical SMILES (for exact matching)
_THIAMINE_CANONICAL = "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"
_GSH_CANONICAL = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"
_DMHF_CANONICAL = "CC1=C(O)OC(C)C1=O"
_HEMF_CANONICAL = "CCC1=C(O)OC(C)C1=O"

# Tier A SMIRKS rules: (name, reaction_family, smirks, ph_gate)
# ph_gate: "any" | "acid" (pH<6) | "neutral" (pH>=6)
_SMIRKS_RULES: List[Tuple[str, str, str, str]] = [
    (
        "schiff_base_lipid",
        "Lipid_Schiff_Base",
        # C3+ aliphatic aldehyde whose alpha-carbon has NO hydroxyl (excludes sugars).
        # The amine donor can be anything with a primary amine on an sp3 carbon (like amino acids).
        "[CX4H2,CX4H3:2][CH1:1]=[O:5].[NH2:3][CX4:4]>>[*:2][CH1:1]=[N:3][*:4].[O:5]",
        "any",
    ),
    (
        "beta_scission_alkoxy",
        "Beta_Scission",
        # Alkoxy radical (R-O.) fragmentation. 
        # Pattern: [C:1]([O:2])-[C:3]. >> [C:1]=[O:2] + [C:3].
        "[CX4:3][CX4:1][OX1H0:2]>>[C:1]=[O:2].[C:3]",
        "any",
    ),
    (
        "radical_propagation_o2",
        "Radical_Propagation_O2",
        # Alkyl radical + O2 -> Peroxy radical
        "[C;X3:1].[O;X1:2]=[O;X1:3]>>[C;X4:1]-[O;X2:2]-[O;X1:3]",
        "any",
    ),
    (
        "peroxy_h_abstraction",
        "Peroxy_H_Abstraction",
        # Peroxy radical + Reactive H (Allylic or specific lipid H) -> Hydroperoxide + Alkyl radical
        # Pattern: [O.] + [H-C-C=C] or [H-C-O]
        # We'll use a broad but restricted pattern to avoid matching EVERY sugar carbon.
        "[O;X1:1]-[O;X2:2].[C;H1,H2,H3;$(C-C=C);!$(C=O):3]>>[O;X2:1]([H])[O;X2:2].[C;X3:3]",
        "any",
    ),
    (
        "radical_termination_disproportionation",
        "Radical_Termination",
        # Two peroxy radicals -> O2 + ... (simplified)
        "[O;X1:1]-[O;X2:2].[O;X1:3]-[O;X2:4]>>[O;X2:1]=[O;X2:3].[O;X2:2].[O;X2:4]",
        "any",
    ),
]

# Pre-compile SMIRKS rules at module load to avoid repeated ReactionFromSmarts() calls
# in the hot enumerate() loop (eliminates ~5 redundant compilations per rule per call).
_COMPILED_SMIRKS_RULES: List[Tuple[str, str, ChemicalReaction, str]] = [
    (name, family, AllChem.ReactionFromSmarts(smirks), gate)
    for name, family, smirks, gate in _SMIRKS_RULES
]

# Module-level cache for SmirksEngine.enumerate() results.
# Key includes canonical precursor pool, max generations, and the full
# ReactionConditions state so cache hits stay semantically correct as the
# rule engine grows new condition-sensitive branches.
_ENUMERATE_CACHE: Dict[tuple, Tuple[ElementaryStep, ...]] = {}


def _normalize_cache_value(value):
    if isinstance(value, float):
        return round(value, 6)
    return value


def _conditions_cache_key(conditions: ReactionConditions) -> Tuple[Tuple[str, object], ...]:
    return tuple(
        sorted((name, _normalize_cache_value(value)) for name, value in vars(conditions).items())
    )


def _clone_species(species: Species) -> Species:
    return Species(label=species.label, smiles=species.smiles)


def _clone_step(step: ElementaryStep) -> ElementaryStep:
    return ElementaryStep(
        reactants=[_clone_species(species) for species in step.reactants],
        products=[_clone_species(species) for species in step.products],
        reaction_family=step.reaction_family,
        rate_constant_k=step.rate_constant_k,
        source_quality=step.source_quality,
        barrier_uncertainty_kcal=step.barrier_uncertainty_kcal,
    )


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


def _is_valid(smi: str) -> bool:
    """Return True if SMILES is parseable and MW is below cap."""
    return _mol(smi) is not None and _mw(smi) <= MAX_MW





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


def _species_from_pool(pool: Set[str], label: str, smiles: str) -> Species:
    """Create a Species, canonicalise its SMILES, and add to pool."""
    can = _canonical(smiles)
    if can:
        pool.add(can)
    return Species(label=label, smiles=smiles)



from src.reaction_templates import (
    _amadori_cascade, _enolisation_steps, _strecker_step, _beta_elimination_steps,
    _aminoketone_condensation, _retro_aldol_fragmentation, _cysteine_degradation,
    _thiazole_condensation, _thiol_addition, _mft_pathway, _furyl_disulfide_formation,
    _sulfur_volatiles_pathway, _deamination_step, _lipid_maillard_synergy,
    _lipid_hydroperoxide_scission, _radical_crosstalk_templates, _sugar_ring_opening,
    _acrylamide_formation, _cml_cel_formation, _thiamine_degradation, _furanone_generation,
    _glutathione_cleavage
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

def _apply_smirks_rule(
    name: str, family: str, rxn: ChemicalReaction, ph_gate: str,
    pool: List[Species], conditions: ReactionConditions
) -> List[ElementaryStep]:
    """Apply a single SMIRKS rule to all relevant species pairs in the pool."""
    # pH gate check
    if ph_gate == "acid" and conditions.pH >= 6:
        return []
    if ph_gate == "neutral" and conditions.pH < 6:
        return []

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
                    for p in prod_tuple:
                        _fix_radicals(p, family, clear_all=True)
                        try:
                            Chem.SanitizeMol(p)
                        except Exception:
                            pass
                        _fix_radicals(p, family)
                    
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
                        try:
                            _fix_radicals(p, family, clear_all=True)
                            # Sanitize to catch valence issues
                            Chem.SanitizeMol(p)
                            _fix_radicals(p, family)
                            ps = Chem.MolToSmiles(p)
                            if _is_valid(ps):
                                prod_smiles.append(ps)
                            else:
                                valid_step = False
                                break
                        except Exception:
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
                            try:
                                Chem.SanitizeMol(p)
                                ps = Chem.MolToSmiles(p)
                                if _is_valid(ps):
                                    prod_smiles.append(ps)
                                else:
                                    valid_step = False
                                    break
                            except Exception:
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

        # Module-level cache: skip re-enumeration for identical precursor pools and conditions.
        # Return clones so callers cannot mutate shared cached objects.
        _cache_key = (
            tuple(sorted(_canonical(s.smiles) for s in precursors if _canonical(s.smiles))),
            max_generations,
            _conditions_cache_key(self.conditions),
        )
        if _cache_key in _ENUMERATE_CACHE:
            return [_clone_step(step) for step in _ENUMERATE_CACHE[_cache_key]]

        # O(1) deduplication set — maintained throughout all Tier B and Tier A phases
        seen_step_keys: Set[str] = set()

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
                k = _step_key(s)
                if k not in seen_step_keys:
                    seen_step_keys.add(k)
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
        sugars = [s for s in pool_list() if is_sugar(s)]
        amines = [s for s in pool_list() if _is_primary_amine(s)]

        for sugar in sugars:
            for amine in amines:
                cascade = _amadori_cascade(sugar, amine)
                for step in cascade:
                    k = _step_key(step)
                    if k not in seen_step_keys:
                        seen_step_keys.add(k)
                        all_steps.append(step)
                        add_step_products(step)

                amadori_smiles = next((p.smiles for s in cascade for p in s.products if "Amadori" in p.label or "Heyns" in p.label), None)
                amadori_can = _canonical(amadori_smiles) if amadori_smiles else None
                amadori_sp = pool_dict.get(amadori_can) if amadori_can else None
                
                if amadori_sp:
                    enols = _enolisation_steps(amadori_sp, sugar, amine, self.conditions)
                    for enol in enols:
                        k = _step_key(enol)
                        if k not in seen_step_keys:
                            seen_step_keys.add(k)
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
                if s_step:
                    k = _step_key(s_step)
                    if k not in seen_step_keys:
                        seen_step_keys.add(k)
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

        # 3e-2. MFT Formation (Phase R.2 Fix)
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
        
        cross_steps = _radical_crosstalk_templates(pool_list())
        _add_steps(cross_steps)

        # ── Tier A: SMIRKS rules, iterative ──────────────────────────────
        for _gen in range(max_generations):
            new_steps_this_gen = []
            current_pool = pool_list()

            for name, family, rxn, gate in _COMPILED_SMIRKS_RULES:
                candidates = _apply_smirks_rule(
                    name, family, rxn, gate, current_pool, self.conditions
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

        _ENUMERATE_CACHE[_cache_key] = tuple(_clone_step(step) for step in all_steps)
        return all_steps


def _step_key(step: ElementaryStep) -> str:
    """Stable hash key for deduplication: sorted reactants + sorted products."""
    reacts = tuple(sorted(_canonical(r.smiles) or r.label for r in step.reactants))
    prods = tuple(sorted(_canonical(p.smiles) or p.label for p in step.products))
    return str((reacts, prods))


# ──────────────────────────────────────────────────────────────────────────
# CLI demo
# ──────────────────────────────────────────────────────────────────────────

