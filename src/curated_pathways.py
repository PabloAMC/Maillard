"""Hand-curated Maillard reaction pathways used by the screening pipeline.

═══════════════════════════════════════════════════════════════════════════════
PARITY NOTE — 2026-08-27 (Wave T4).  THIS LAYER IS ONE WAVE BEHIND THE ENGINE.
═══════════════════════════════════════════════════════════════════════════════

`PATHWAYS` is a HAND-MAINTAINED MIRROR of the chemistry the SMIRKS engine
enumerates (`src/reaction_templates.py` + `src/smirks_engine.py`).  Nothing keeps
the two in step: no `src/` module imports this file, and no test compares the two
family vocabularies.  Drift here is therefore silent, and it has happened.

MEASURED STATUS (families emitted by `PATHWAYS`, counted 2026-08-27):

  * WAVE N — MIRRORED.  The corrected MFT route is present
    (`Deoxyosone_Reduction` + `Thiol_Addition_Pentodiulose`), the retired
    `Thiol_Addition_Norfuraneol` step is gone, norfuraneol survives as a terminal
    furanone product, and the block cites both Cerny & Davidek DOIs.

  * WAVE P — NOT MIRRORED.  NONE of the six Wave P families appear here:
        Mercaptoketone_Formation        Mercaptoketone_Aldol_Addition
        Mercaptoketone_Cyclodehydration  (the C2+C3 recombination lane to MFT,
                                          Hofmann & Schieberle 1998)
        Fructofuranosyl_Dehydration      (fructose's own ring-retained HMF route)
        Furanone_Reductive_Opening       Furanone_Amino_Acid_Reduction
    The Wave P ledger section contains ZERO occurrences of the string "curated",
    so this was not a decision to skip the mirror — it appears simply not to have
    been considered.

  * The ketose lane is likewise absent: the engine emits `Heyns_Rearrangement`
    for fructose + an amino acid (Wave T4); this layer has only the aldose
    `Amadori_Rearrangement`.

CONSEQUENCE THE READER MUST KNOW ABOUT.  `scripts/generate_reaction_network.py`
draws `docs/assets/reaction_network.pdf` FROM THIS FILE.  THE PUBLISHED
ARCHITECTURE FIGURE THEREFORE DEPICTS PRE-WAVE-P CHEMISTRY: it shows neither the
C2+C3 recombination lane to MFT nor fructose's independent HMF route.  It is a
schematic of the core lanes, not a picture of the shipped network.  Read it that
way until the figure is regenerated from a re-synced layer.

THIS LAYER IS NOT DEAD — it is the input to that figure and to the chemistry
invariants in `tests/unit/test_chemistry_soundness.py` and
`tests/unit/test_data_integrity.py`.

OWNER DECISION, filed as [P] in `tasks/audit_remediation.md` under Wave T4:
either (i) mirror Wave P here and regenerate the figure, or (ii) declare this
layer DELIBERATELY FROZEN at the Wave N topology and say so in the figure's own
caption.  Wave T4 deliberately did NOT add the chemistry: hand-adding six
families to a mirror nothing tests is how the drift got here, and the right fix
is probably a parity TEST rather than another hand-sync.
"""

from src.barrier_constants import get_barrier
from src.pathway_extractor import ElementaryStep, Species

# Computational priors that are LOADED but deliberately NOT WIRED to any step,
# with the reason. Recorded explicitly so a parked prior cannot be mistaken for
# an oversight and quietly re-attached to the nearest-looking family.
_PARKED_COMPUTATIONAL_PRIORS = {
    "quinone_cys_michael": (
        "SLR family 13 (cysteine thiol Michael addition to a polyphenol "
        "quinone). Parked 2026-08-27: it was being averaged into the barrier "
        "for furfural + H2S -> furfurylthiol, an unrelated reaction, where its "
        "6.93 kcal/mol dominated the mean. No family-13 quinone/cysteine step "
        "exists in the curated layer, so there is nothing to route it to."
    ),
}


def _species(label: str, smiles: str) -> Species:
    return Species(label=label, smiles=smiles)


RIBOSE = _species("D-ribose", "O=CC(O)C(O)C(O)CO")
GLUCOSE = _species("D-glucose", "O=CC(O)C(O)C(O)C(O)CO")
GLYCINE = _species("glycine", "NCC(=O)O")
CYSTEINE = _species("L-cysteine", "NC(CS)C(=O)O")
LEUCINE = _species("L-leucine", "CC(C)CC(N)C(=O)O")
LYSINE = _species("L-lysine", "NCCCCC(N)C(=O)O")
ASPARAGINE = _species("L-asparagine", "NC(=O)CC(N)C(=O)O")
ACRYLAMIDE = _species("acrylamide", "C=CC(=O)N")
MFT = _species("2-methyl-3-furanthiol", "Cc1c(S)cco1")
MFT_DIMER = _species("bis(2-methyl-3-furyl) disulfide", "Cc1c(SSc2ccoc2C)cco1")
DIMETHYLPYRAZINE_25 = _species("2,5-dimethylpyrazine", "Cc1cnc(C)cn1")
GLUCOSE_ASN_SCHIFF = _species("glucose-asparagine-Schiff-base", "OCC(O)C(O)C(O)C(O)/C=N/C(CC(N)=O)C(=O)O")
GLUCOSE_REMAINDER = _species("sugar-derived-amine-fragment", "NCC(O)C(O)C(O)C(O)C=O")
WATER = _species("water", "O")
CO2 = _species("CO2", "O=C=O")
H2S = _species("H2S", "S")
NH3 = _species("ammonia", "N")
PYRUVALDEHYDE = _species("pyruvaldehyde", "CC(=O)C=O")
FURFURAL = _species("furfural", "O=Cc1ccco1")
FFT = _species("2-furfurylthiol", "SCc1ccco1")
METHYLBUTANAL_3 = _species("3-methylbutanal", "CC(C)CC=O")
ACETALDEHYDE = _species("acetaldehyde", "CC=O")
DHA = _species("dehydroalanine", "C=C(N)C(=O)O")
LAL = _species("lysinoalanine", "NC(CCCCNCC(N)C(=O)O)C(=O)O")
RIBOSE_GLY_SCHIFF = _species("ribose-glycine-Schiff-base", "OCC(O)C(O)C(O)/C=N/CC(=O)O")
RIBOSE_GLY_AMADORI = _species("ribose-glycine-Amadori", "OCC(O)C(O)C(=O)CNCC(=O)O")
RIBOSE_CYS_SCHIFF = _species("ribose-cysteine-Schiff-base", "OCC(O)C(O)C(O)/C=N/C(CS)C(=O)O")
RIBOSE_CYS_AMADORI = _species("ribose-cysteine-Amadori", "OCC(O)C(O)C(=O)CNC(CS)C(=O)O")
GLUCOSE_GLY_SCHIFF = _species("glucose-glycine-Schiff-base", "OCC(O)C(O)C(O)C(O)/C=N/CC(=O)O")
GLUCOSE_GLY_AMADORI = _species("glucose-glycine-Amadori", "OCC(O)C(O)C(O)C(=O)CNCC(=O)O")
DEOXYOSONE_3 = _species("3-deoxyosone", "O=CC(=O)CC(O)CO")
# 1-deoxy-2,3-pentodiulose and the 1,4-dideoxyosone: the isotope-evidenced MFT
# precursors (Cerny & Davidek 2003, 10.1021/jf026123f; 2004, 10.1021/jf035265m).
# Norfuraneol (van den Ouweland & Peer 1975, 10.1021/jf60199a045 — a synthesis
# route, not the in-situ one) was retired from the MFT lane by Wave N 2026-08-27;
# the species is kept as a genuine furanone product.
DEOXYOSONE_1 = _species("pentose-1-deoxyosone", "CC(=O)C(=O)C(O)CO")
NORFURANEOL = _species("norfuraneol", "CC1=C(O)C(=O)CO1")
DIDEOXYOSONE_14 = _species("1,4-dideoxypentodiulose", "CC(=O)C(=O)CCO")
GLUCOSE_DEOXYOSONE_3 = _species("glucose-3-deoxyosone", "O=CC(=O)CC(O)C(O)CO")
HMF = _species("HMF", "OCC1=CC=C(C=O)O1")
HEXANAL = _species("hexanal", "CCCCCC=O")
HEXANAL_GLY_SCHIFF = _species("hexanal-glycine-Schiff-base", "CCCCC/C=N/CC(=O)O")
HEXANAL_LYS_SCHIFF = _species("hexanal-lysine-Schiff-base", "CCCCC/C=N/CCCCC(N)C(=O)O")
AMINOACETONE = _species("aminoacetone", "CC(=O)CN")
HYDROGEN = _species("hydrogen", "[HH]")


def _step(reactants: list[Species], products: list[Species], reaction_family: str) -> ElementaryStep:
    return ElementaryStep(reactants=reactants, products=products, reaction_family=reaction_family)


# AUDIT 2026-08-27 (Wave G1 fix 9). Two things were repaired here.
#
# (a) FAMILY VOCABULARY. The curated layer used family strings the engine never
#     emits ("Enolisation", "Sugar_Dehydration", and a pyrazine step mislabelled
#     "Strecker_Degradation"). That split the two layers into separate
#     vocabularies: `_wire_computational_priors` below keyed on the curated
#     names, while the FAST barrier table keys on the engine names -- and
#     "Sugar_Dehydration" in particular canonicalises to nothing, so it fell
#     through to DEFAULT_BARRIER whenever the FAST path was used. Both layers
#     now speak the engine's vocabulary.
#
# (b) PATHWAY C was not self-contained: it consumed DEOXYOSONE_3, which nothing
#     in it produced, and reached furfural through a one-step ribose -> furfural
#     lump. It now carries its own Schiff/Amadori/enolisation trunk (built on
#     cysteine, so its declared `consumes` list is unchanged) and its MFT step
#     follows the accepted norfuraneol route instead of the retired one-step
#     3-deoxyosone + H2S shortcut.
#
# HYDROGEN in pathways B and C is a documented REDUCING-EQUIVALENT TOKEN, not a
# pathway intermediate: the real donors/acceptors are the sugar-derived
# reductones and dissolved O2, for which this model carries no species. Same
# lumping convention as src/reaction_templates.py.
#
# 2026-08-27 (Wave I fix 8, red-team H4) — PARITY NOTE, deliberately NOT acted on
# here. The engine lane had a real bug from this token: in src/smirks_engine.py the
# MFT step is POOL-GATED on `[HH]`, whose only producer reachable from a
# ribose/cysteine system was the pyrazine aromatisation, so disabling the pyrazine
# lane drove predicted MFT to exactly 0.0 ppb. That is fixed engine-side by
# `src/reaction_templates._thiol_reductant_pool` (2 cysteine -> cystine + 2[H]).
# THIS layer cannot have that failure mode -- a curated pathway is a fixed list of
# steps, not a reachability search, so nothing here can be gated off -- but note
# that pathway C's only HYDROGEN producer (`MFT + MFT -> dimer + H2`) sits
# DOWNSTREAM of the MFT step that consumes it. That is a bookkeeping circularity,
# harmless today, and it would become load-bearing the moment this layer is ever
# fed through a flux propagator. Adding the cysteine -> cystine step here for
# parity is an OPEN OWNER ITEM; it was left out of Wave I because it changes the
# curated step inventory, which is pinned elsewhere.
PATHWAYS = {
    "A_Core_Maillard_Ribose_Gly": [
        _step([RIBOSE, GLYCINE], [RIBOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([RIBOSE_GLY_SCHIFF], [RIBOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([RIBOSE_GLY_AMADORI], [DEOXYOSONE_3, GLYCINE], "Enolisation_Intermediate"),
        _step([DEOXYOSONE_3], [FURFURAL, WATER, WATER], "Enolisation_1_2"),
    ],
    "A_Core_Maillard_Glucose_Gly": [
        _step([GLUCOSE, GLYCINE], [GLUCOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([GLUCOSE_GLY_SCHIFF], [GLUCOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([GLUCOSE_GLY_AMADORI], [GLUCOSE_DEOXYOSONE_3, GLYCINE], "Enolisation_Intermediate"),
        _step([GLUCOSE_DEOXYOSONE_3], [HMF, WATER, WATER], "Enolisation_1_2"),
    ],
    "B_Strecker_Leu": [
        _step([PYRUVALDEHYDE, LEUCINE], [METHYLBUTANAL_3, AMINOACETONE, CO2], "Strecker_Degradation"),
        # This is the aminoketone self-condensation to the pyrazine, not a
        # Strecker degradation; it was carrying the Strecker family (and hence
        # the Strecker barrier) purely by mislabelling.
        _step([AMINOACETONE, AMINOACETONE], [DIMETHYLPYRAZINE_25, WATER, WATER, HYDROGEN], "Aminoketone_Condensation"),
    ],
    "C_S_Maillard_FFT": [
        _step([CYSTEINE, WATER], [ACETALDEHYDE, H2S, NH3, CO2], "Cysteine_Degradation"),
        _step([RIBOSE, CYSTEINE], [RIBOSE_CYS_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([RIBOSE_CYS_SCHIFF], [RIBOSE_CYS_AMADORI], "Amadori_Rearrangement"),
        # 1,2-enolisation limb -> 3-deoxyosone -> furfural -> FFT
        _step([RIBOSE_CYS_AMADORI], [DEOXYOSONE_3, CYSTEINE], "Enolisation_Intermediate"),
        _step([DEOXYOSONE_3], [FURFURAL, WATER, WATER], "Enolisation_1_2"),
        _step([FURFURAL, H2S, HYDROGEN], [FFT, WATER], "Thiol_Addition"),
        # 2,3-enolisation limb -> 1-deoxyosone -> 1,4-dideoxyosone -> MFT.
        # ROUTE CORRECTION 2026-08-27 (Wave N): the norfuraneol -> MFT step was
        # retired on isotope evidence (Cerny & Davidek 2003, 10.1021/jf026123f:
        # spiked norfuraneol does NOT label the MFT; 2004, 10.1021/jf035265m:
        # 1,4-dideoxypento-2,3-diulose positionally confirmed). Norfuraneol
        # stays as a terminal furanone product of this pathway. Balances
        # (RDKit-verified): C5H8O4 + H2 -> C5H8O3 + H2O; C5H8O3 + H2S ->
        # C5H6OS + 2 H2O (no reducing token on the H2S step). Mirrors
        # `_furanone_and_mft_route` in src/reaction_templates.py.
        _step([RIBOSE_CYS_AMADORI], [DEOXYOSONE_1, CYSTEINE], "Enolisation_2_3_Amadori"),
        _step([DEOXYOSONE_1], [NORFURANEOL, WATER], "Furanone_Cyclisation"),
        _step([DEOXYOSONE_1, HYDROGEN], [DIDEOXYOSONE_14, WATER], "Deoxyosone_Reduction"),
        _step([DIDEOXYOSONE_14, H2S], [MFT, WATER, WATER], "Thiol_Addition_Pentodiulose"),
        _step([MFT, MFT], [MFT_DIMER, HYDROGEN], "Thiol_Oxidation"),
    ],
    "D_Offflavour_Trapping_Gly": [
        _step([HEXANAL, GLYCINE], [HEXANAL_GLY_SCHIFF, WATER], "Lipid_Schiff_Base"),
    ],
    "D_Offflavour_Trapping_Lys": [
        _step([HEXANAL, LYSINE], [HEXANAL_LYS_SCHIFF, WATER], "Lipid_Schiff_Base"),
    ],
    "E_DHA_Competition": [
        _step([CYSTEINE], [DHA, H2S], "Beta_Elimination"),
        _step([DHA, LYSINE], [LAL], "DHA_Crosslinking"),
    ],
    "F_Safety_Acrylamide": [
        _step([GLUCOSE, ASPARAGINE], [GLUCOSE_ASN_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([GLUCOSE_ASN_SCHIFF], [ACRYLAMIDE, CO2, GLUCOSE_REMAINDER], "Beta_Elimination"),
    ],
}


PATHWAY_METADATA = {
    "A_Core_Maillard_Ribose_Gly": {
        "target": FURFURAL,
        "type": "desirable",
        "consumes": ["D-ribose", "glycine"],
        "toxicity_flag": None,
    },
    "A_Core_Maillard_Glucose_Gly": {
        "target": HMF,
        "type": "desirable",
        "consumes": ["D-glucose", "glycine"],
        "toxicity_flag": "5-Hydroxymethylfurfural (HMF)",
    },
    "B_Strecker_Leu": {
        "target": METHYLBUTANAL_3,
        "type": "desirable",
        "consumes": ["pyruvaldehyde", "L-leucine"],
        "toxicity_flag": None,
    },
    "C_S_Maillard_FFT": {
        "target": FFT,
        "type": "desirable",
        "consumes": ["L-cysteine", "D-ribose"],
        "toxicity_flag": None,
    },
    "D_Offflavour_Trapping_Gly": {
        "target": HEXANAL_GLY_SCHIFF,
        "type": "masking",
        "consumes": ["hexanal", "glycine"],
        "toxicity_flag": None,
    },
    "D_Offflavour_Trapping_Lys": {
        "target": HEXANAL_LYS_SCHIFF,
        "type": "masking",
        "consumes": ["hexanal", "L-lysine"],
        "toxicity_flag": None,
    },
    "E_DHA_Competition": {
        "target": LAL,
        "type": "competing",
        "consumes": ["L-cysteine", "L-lysine"],
        "toxicity_flag": "Lysinoalanine (LAL)",
    },
    "F_Safety_Acrylamide": {
        "target": ACRYLAMIDE,
        "type": "competing",
        "consumes": ["D-glucose", "L-asparagine"],
        "toxicity_flag": "Acrylamide",
    },
}


def _wire_computational_priors():
    import json
    import math
    from pathlib import Path

    # 1. Load computational priors
    root = Path(__file__).resolve().parents[1]
    priors_path = root / "data" / "lit" / "computational_priors.json"
    if not priors_path.exists():
        return

    try:
        with open(priors_path, "r", encoding="utf-8") as f:
            priors_data = json.load(f)
    except Exception:
        return

    # Helper to calculate rate constant k at 150C (423.15 K) using Eyring equation
    # k = (kB * T / h) * exp(-Delta G‡ / RT)
    def calculate_k(barrier_kcal):
        T = 423.15  # 150 C
        kb = 1.380649e-23
        h = 6.62607015e-34
        R = 8.314462618
        pre_exponential = (kb * T) / h
        j_per_mol = barrier_kcal * 4184.0
        exponent = -j_per_mol / (R * T)
        return pre_exponential * math.exp(exponent)

    # Extract prior kinetic values dynamically
    pe_schiff_base_barrier = 22.21
    pe_schiff_base_unc = 5.0
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "pe_schiff_base":
            if e.get("barrier_kcal_mol") is not None:
                pe_schiff_base_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                pe_schiff_base_unc = float(e["uncertainty_kj"]) / 4.184

    glucose_ea_kcal = 25.33  # 106 kJ/mol — van Boekel (2001) fallback
    ribose_ea_kcal = 27.25   # 114 kJ/mol — van Boekel (2001) fallback
    for e in priors_data.get("carbonyl_donor_priors", []):
        if e.get("id") == "van_boekel_2001_maillard_kinetics_review_v1":
            ea_by_sugar = e.get("ea_kj_per_mol_by_sugar", {})
            if "glucose" in ea_by_sugar:
                glucose_ea_kcal = float(ea_by_sugar["glucose"]) / 4.184
            if "ribose" in ea_by_sugar:
                ribose_ea_kcal = float(ea_by_sugar["ribose"]) / 4.184
            break

    pe_amadori_barrier = 19.81
    pe_amadori_unc = 5.0
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "pe_amadori":
            if e.get("barrier_kcal_mol") is not None:
                pe_amadori_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                pe_amadori_unc = float(e["uncertainty_kj"]) / 4.184

    egcg_ea_kcal = 18.6
    for e in priors_data.get("dicarbonyl_sink_priors", []):
        if e.get("id") == "jafc_2020_egcg_deoxyosone_trapping_v1":
            if e.get("baseline_arp_formation_ea_kj_per_mol") is not None:
                egcg_ea_kcal = float(e["baseline_arp_formation_ea_kj_per_mol"]) / 4.184

    huang_ea_kcal = 27.5
    for e in priors_data.get("carbonyl_donor_priors", []):
        if e.get("id") == "huang_2024_dixyl_arp_degradation":
            if e.get("bifurcated_degradation_kinetics", {}).get("three_deoxyglucosone_mediated_furosine_ea_kj_mol") is not None:
                huang_ea_kcal = float(e["bifurcated_degradation_kinetics"]["three_deoxyglucosone_mediated_furosine_ea_kj_mol"]) / 4.184

    furosine_from_3dg_barrier = 19.53
    furosine_from_3dg_unc = 3.35
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "furosine_from_3dg":
            if e.get("barrier_kcal_mol") is not None:
                furosine_from_3dg_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                furosine_from_3dg_unc = float(e["uncertainty_kj"]) / 4.184

    pyrraline_from_3dg_barrier = 12.77
    pyrraline_from_3dg_unc = 0.96
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "pyrraline_from_3dg":
            if e.get("barrier_kcal_mol") is not None:
                pyrraline_from_3dg_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                pyrraline_from_3dg_unc = float(e["uncertainty_kj"]) / 4.184

    martins_furosine_ea = 28.5
    martins_pyrraline_ea = 31.2
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("id") == "martins_2003_lys_glucose_kinetics_v1":
            params = e.get("kinetic_parameters", {})
            if params.get("furosine_formation_ea_kcal_mol") is not None:
                martins_furosine_ea = float(params["furosine_formation_ea_kcal_mol"])
            if params.get("pyrraline_formation_ea_kcal_mol") is not None:
                martins_pyrraline_ea = float(params["pyrraline_formation_ea_kcal_mol"])

    frontiers_ea_kcal = 5.66
    for e in priors_data.get("ascorbic_pathway_priors", []):
        if e.get("id") == "frontiers_2022_hcw_aa_arrhenius_v1":
            ea_by_ph = e.get("ea_kj_per_mol_by_ph", {})
            if "5.0" in ea_by_ph and "7.0" in ea_by_ph:
                frontiers_ea_kcal = (float(ea_by_ph["5.0"]) + float(ea_by_ph["7.0"])) / 2.0 / 4.184

    aa_ring_open_barrier = 7.58
    aa_ring_open_unc = 4.78
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "aa_ring_open_dicarbonyl":
            if e.get("barrier_kcal_mol") is not None:
                aa_ring_open_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                aa_ring_open_unc = float(e["uncertainty_kj"]) / 4.184

    lagrain_elim_ea_kcal = 21.08
    lagrain_lanth_ea_kcal = 1.56
    for e in priors_data.get("crosslink_kinetics_priors", []):
        if e.get("id") == "lagrain_2010_cystine_elimination_lanthionine":
            elim = e.get("elimination_kinetics", {})
            if "ea_ph6_kj_mol" in elim:
                lagrain_elim_ea_kcal = float(elim["ea_ph6_kj_mol"]) / 4.184
            lanth = e.get("lanthionine_formation_kinetics", {})
            if "ea_ph6_kj_mol" in lanth:
                lagrain_lanth_ea_kcal = float(lanth["ea_ph6_kj_mol"]) / 4.184

    # PARKED 2026-08-27 -- still read so the value stays visible in the load
    # path and in any future re-routing, but deliberately not wired to a step.
    # See `_PARKED_COMPUTATIONAL_PRIORS`.
    quinone_cys_barrier = 6.93
    quinone_cys_unc = 3.58
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "quinone_cys_michael":
            if e.get("barrier_kcal_mol") is not None:
                quinone_cys_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                quinone_cys_unc = float(e["uncertainty_kj"]) / 4.184

    lal_crosslink_barrier = 16.0
    lal_crosslink_unc = 10.0
    for e in priors_data.get("dft_kinetic_priors", {}).get("entries", []):
        if e.get("reaction_key") == "lysinoalanine_crosslink":
            if e.get("barrier_kcal_mol") is not None:
                lal_crosslink_barrier = float(e["barrier_kcal_mol"])
            if e.get("uncertainty_kj") is not None:
                lal_crosslink_unc = float(e["uncertainty_kj"]) / 4.184

    # Map reactions
    for path_name, steps in PATHWAYS.items():
        for step in steps:
            fam = step.reaction_family
            reactants = [r.label for r in step.reactants]
            
            barrier = None
            uncertainty = 5.0
            
            if fam == "Schiff_Base_Formation":
                if "D-glucose" in reactants:
                    barrier = (pe_schiff_base_barrier * 2 + glucose_ea_kcal * 2) / 4.0
                    uncertainty = (pe_schiff_base_unc * 2 + 3.0 * 2) / 4.0
                elif "D-ribose" in reactants:
                    barrier = (pe_schiff_base_barrier * 2 + ribose_ea_kcal * 2) / 4.0
                    uncertainty = (pe_schiff_base_unc * 2 + 3.0 * 2) / 4.0
                else:
                    barrier = pe_schiff_base_barrier
                    uncertainty = pe_schiff_base_unc
            
            elif fam == "Amadori_Rearrangement":
                barrier = (pe_amadori_barrier * 2 + egcg_ea_kcal * 1 + huang_ea_kcal * 2) / 5.0
                uncertainty = (pe_amadori_unc * 2 + 3.5 * 1 + 3.0 * 2) / 5.0
                
            elif fam == "Enolisation_Intermediate":
                # 2026-08-27: key renamed from "Enolisation" (a family string the
                # engine never emits) to the engine's own family name. Numbers
                # unchanged.
                if "glucose-glycine-Amadori" in reactants:
                    barrier = (furosine_from_3dg_barrier * 2 + martins_furosine_ea * 3) / 5.0
                    uncertainty = (furosine_from_3dg_unc * 2 + 1.5 * 3) / 5.0
                else:
                    barrier = (pyrraline_from_3dg_barrier * 2 + martins_pyrraline_ea * 3) / 5.0
                    uncertainty = (pyrraline_from_3dg_unc * 2 + 1.5 * 3) / 5.0

            elif fam == "Enolisation_1_2":
                # 2026-08-27: key renamed from "Sugar_Dehydration". Numbers
                # unchanged, but flagged: [P] CARRIED FORWARD -- both priors
                # averaged here (`frontiers_2022_hcw_aa_arrhenius_v1` and
                # `aa_ring_open_dicarbonyl`) are ASCORBIC ACID priors, and the
                # ~6.6 kcal/mol they produce for deoxyosone -> furfural
                # dehydration is far below the FAST table's literature-informed
                # 28.0 for the same step. This is a provenance mismatch of the
                # same kind as the quinone_cys_michael one fixed below; it was
                # NOT in this wave's mandate and is left for an owner call
                # rather than silently re-anchored.
                barrier = (frontiers_ea_kcal * 2 + aa_ring_open_barrier * 2) / 4.0
                uncertainty = (3.0 * 2 + aa_ring_open_unc * 2) / 4.0

            elif fam == "Enolisation_2_3_Amadori":
                # New in the curated layer 2026-08-27 (the accepted MFT route).
                # Inherits the enolisation family constant, as instructed.
                barrier, uncertainty = get_barrier("Enolisation_2_3_Amadori")

            elif fam == "Furanone_Cyclisation":
                barrier, uncertainty = get_barrier("Furanone_Cyclisation")

            elif fam == "Strecker_Degradation":
                barrier = 22.67
                uncertainty = 1.5

            elif fam == "Aminoketone_Condensation":
                # Previously mislabelled "Strecker_Degradation" and therefore
                # silently given the Strecker barrier. Now takes the FAST
                # table's `aminoketone_condensation` value (itself flagged there
                # as an unconstrained legacy fit).
                barrier, uncertainty = get_barrier("Aminoketone_Condensation")

            elif fam == "Cysteine_Degradation":
                barrier = lagrain_elim_ea_kcal
                uncertainty = 1.5

            elif fam == "Thiol_Addition":
                # 2026-08-27: the `quinone_cys_michael` DFT prior was PARKED
                # here. It is an SLR family-13 observable -- Michael addition of
                # a cysteine thiol to a POLYPHENOL QUINONE -- and was being
                # averaged into the barrier for furfural + H2S, a completely
                # different reaction; at 6.93 kcal/mol it dominated the average
                # and pulled the step to 13.6 kcal/mol. There is no family-13
                # quinone/cysteine step in the curated layer to route it to, so
                # it is parked rather than re-pointed: see
                # `_PARKED_COMPUTATIONAL_PRIORS` at the top of this module. The
                # step now takes the
                # FAST table's `thiol_addition` value, which is the only
                # sulfur-branch constant with a surviving literature constraint
                # (Hofmann 1998, 10.1021/jf9705983).
                # CORRECTED 2026-08-27 (Wave S2c): the two lines above are FALSE and are
                # kept because they record what the repo believed when this routing was
                # written. `thiol_addition` has NO literature constraint. Its "Hofmann
                # window" [28.10, 28.85] was derived from cys_ribose_140C_Hofmann1998,
                # whose MFT 342 / FFT 200 ppb Wave S2b traced to
                # data/benchmarks/maillard_validation_benchmarks.md section 1.3 -- an
                # abstract-reconstructed table committed in the same commit as the
                # benchmark, with both values interior points of two invented mol % bands.
                # THE SULFUR BRANCH HAS ZERO ABSOLUTE LITERATURE ANCHORS.
                # (Superseded 2026-08-28, Wave W: it now has three, from the real
                # Table 1 of 10.1021/jf9705983 -- and the model fails all three by
                # 12-30x. Nothing here was refitted to them. The sentence is kept
                # because the provenance point it was making is unchanged: the
                # constant this comment guards has never been anchored to any
                # literature value, before or after Wave W.) The ROUTING
                # DECISION here is unchanged and is still the right one: `thiol_addition`
                # remains the closest sulfur-addition class analogue for this parked step,
                # and parking it beats re-pointing it at a family that does not exist. What
                # changes is only what may be claimed for the number it returns -- it is a
                # class estimate, not a measured barrier. See src/barrier_constants.py.
                barrier, uncertainty = get_barrier("Thiol_Addition")

            # Wave N 2026-08-27: the two families of the corrected MFT route
            # (Deoxyosone_Reduction, Thiol_Addition_Pentodiulose) both have
            # explicit FAST_BARRIERS entries, so the plain get_barrier lookup
            # below resolves them exactly. Thiol_Addition_Norfuraneol no longer
            # appears in any pathway.
            elif fam in ("Deoxyosone_Reduction", "Thiol_Addition_Pentodiulose"):
                barrier, uncertainty = get_barrier(fam)

            elif fam == "Thiol_Oxidation":
                barrier, uncertainty = get_barrier("Thiol_Oxidation")

            elif fam == "Lipid_Schiff_Base":
                barrier = (pe_schiff_base_barrier * 2 + 23.66 * 3) / 5.0
                uncertainty = (pe_schiff_base_unc * 2 + 1.5 * 3) / 5.0
                
            elif fam == "Beta_Elimination":
                barrier = lagrain_elim_ea_kcal
                uncertainty = 1.5
                
            elif fam == "DHA_Crosslinking":
                barrier = (lal_crosslink_barrier * 2 + lagrain_lanth_ea_kcal * 3) / 5.0
                uncertainty = (lal_crosslink_unc * 2 + 1.5 * 3) / 5.0

            if barrier is not None:
                step.barrier_kcal_mol = round(barrier, 3)
                step.barrier_uncertainty_kcal = round(uncertainty, 3)
                step.rate_constant_k = round(calculate_k(barrier), 7)
                step.source_quality = "literature"


# Run the wiring
_wire_computational_priors()

