"""Hand-curated Maillard reaction pathways used by the screening pipeline."""

from src.pathway_extractor import ElementaryStep, Species


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
GLUCOSE_GLY_SCHIFF = _species("glucose-glycine-Schiff-base", "OCC(O)C(O)C(O)C(O)/C=N/CC(=O)O")
GLUCOSE_GLY_AMADORI = _species("glucose-glycine-Amadori", "OCC(O)C(O)C(O)C(=O)CNCC(=O)O")
DEOXYOSONE_3 = _species("3-deoxyosone", "O=CC(=O)CC(O)CO")
GLUCOSE_DEOXYOSONE_3 = _species("glucose-3-deoxyosone", "O=CC(=O)CC(O)C(O)CO")
HMF = _species("HMF", "OCC1=CC=C(C=O)O1")
HEXANAL = _species("hexanal", "CCCCCC=O")
HEXANAL_GLY_SCHIFF = _species("hexanal-glycine-Schiff-base", "CCCCC/C=N/CC(=O)O")
HEXANAL_LYS_SCHIFF = _species("hexanal-lysine-Schiff-base", "CCCCC/C=N/CCCCC(N)C(=O)O")
AMINOACETONE = _species("aminoacetone", "CC(=O)CN")
HYDROGEN = _species("hydrogen", "[HH]")


def _step(reactants: list[Species], products: list[Species], reaction_family: str) -> ElementaryStep:
    return ElementaryStep(reactants=reactants, products=products, reaction_family=reaction_family)


PATHWAYS = {
    "A_Core_Maillard_Ribose_Gly": [
        _step([RIBOSE, GLYCINE], [RIBOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([RIBOSE_GLY_SCHIFF], [RIBOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([RIBOSE_GLY_AMADORI], [DEOXYOSONE_3, GLYCINE], "Enolisation"),
        _step([DEOXYOSONE_3], [FURFURAL, WATER, WATER], "Sugar_Dehydration"),
    ],
    "A_Core_Maillard_Glucose_Gly": [
        _step([GLUCOSE, GLYCINE], [GLUCOSE_GLY_SCHIFF, WATER], "Schiff_Base_Formation"),
        _step([GLUCOSE_GLY_SCHIFF], [GLUCOSE_GLY_AMADORI], "Amadori_Rearrangement"),
        _step([GLUCOSE_GLY_AMADORI], [GLUCOSE_DEOXYOSONE_3, GLYCINE], "Enolisation"),
        _step([GLUCOSE_DEOXYOSONE_3], [HMF, WATER, WATER], "Sugar_Dehydration"),
    ],
    "B_Strecker_Leu": [
        _step([PYRUVALDEHYDE, LEUCINE], [METHYLBUTANAL_3, AMINOACETONE, CO2], "Strecker_Degradation"),
        _step([AMINOACETONE, AMINOACETONE], [DIMETHYLPYRAZINE_25, WATER, WATER, HYDROGEN], "Strecker_Degradation"),
    ],
    "C_S_Maillard_FFT": [
        _step([CYSTEINE, WATER], [ACETALDEHYDE, H2S, NH3, CO2], "Cysteine_Degradation"),
        _step([RIBOSE], [FURFURAL, WATER, WATER, WATER], "Sugar_Dehydration"),
        _step([FURFURAL, H2S, HYDROGEN], [FFT, WATER], "Thiol_Addition"),
        _step([DEOXYOSONE_3, H2S, HYDROGEN], [MFT, WATER, WATER, WATER], "Thiol_Addition"),
        _step([MFT, MFT], [MFT_DIMER, HYDROGEN], "Thiol_Addition"),
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
                
            elif fam == "Enolisation":
                if "glucose-glycine-Amadori" in reactants:
                    barrier = (furosine_from_3dg_barrier * 2 + martins_furosine_ea * 3) / 5.0
                    uncertainty = (furosine_from_3dg_unc * 2 + 1.5 * 3) / 5.0
                else:
                    barrier = (pyrraline_from_3dg_barrier * 2 + martins_pyrraline_ea * 3) / 5.0
                    uncertainty = (pyrraline_from_3dg_unc * 2 + 1.5 * 3) / 5.0
                    
            elif fam == "Sugar_Dehydration":
                barrier = (frontiers_ea_kcal * 2 + aa_ring_open_barrier * 2) / 4.0
                uncertainty = (3.0 * 2 + aa_ring_open_unc * 2) / 4.0
                
            elif fam == "Strecker_Degradation":
                barrier = 22.67
                uncertainty = 1.5
                
            elif fam == "Cysteine_Degradation":
                barrier = lagrain_elim_ea_kcal
                uncertainty = 1.5
                
            elif fam == "Thiol_Addition":
                barrier = (quinone_cys_barrier * 2 + 18.0 * 3) / 5.0
                uncertainty = (quinone_cys_unc * 2 + 1.5 * 3) / 5.0
                
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

