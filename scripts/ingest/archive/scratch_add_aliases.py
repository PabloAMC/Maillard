import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# File paths
registry_path = ROOT / "data" / "lit" / "benchmark_intake_registry.json"
matrix_path = ROOT / "data" / "lit" / "slr_incorporation_matrix.json"
priors_path = ROOT / "data" / "lit" / "computational_priors.json"
calibrations_path = ROOT / "data" / "lit" / "process_state_calibrations.json"

aliases_to_add = {
    "pmc_4419266_pe_interfacial_maillard_kinetics": ["PMC PMCID:PMC4419266 (Ref. 18)"],
    "pmc5992167_amadori_pe_burden_anchor": ["PMC PMCID:PMC5992167 (Refs. 16/17)"],
    "pmc12155365_sunflower_roasted_anchor": ["PMC12155365 (2025)"],
    "pmc11049305_spirulina_offnote_anchor": ["PMC11049305 (2024)"],
    "pmc9905368_spi_hvp_xylose_benchmark": ["PMC9905368 (2023)"],
    "pmc9905368_wheat_gluten_hvp_xylose_benchmark": ["PMC9905368 (2023)"],
    "jafc_2019_ref24_polyphenol_thiol_capping": ["J. Agric. Food Chem. 2019 (Ref. 24)"],
    "pmid_1904866_pentosidine_equivalence_anchor": ["PubMed PMID:1904866 (Ref. 5)"],
    "pmc9351765_crosspy_trapping_anchor": ["PMC PMCID:PMC9351765 (Ref. 12)"],
    "pmc_2024_pba_cml_cel_ranges_anchor": ["PMC 2024 (PMC12451096)"],
    "pmc_12648097_acrylamide_mitigation_anchor": ["PMC PMCID:PMC12648097 (Ref. 5)"],
    "wageningen_ref9_hme_rework_hydration_anchor": ["Wageningen Ref. 9"],
    "jafc_2019_ref21_pea_gum_arabic_architecture_anchor": ["J. Agric. Food Chem. 2019 (Ref. 21)"],
    "jafc_2022_melanoidin_thiol_binding": ["J. Agric. Food Chem. 2022 (Ref. 41) + Ref. 38"],
    "acrylamide_heat_trapping_2024": ["ResearchGate ref. 10 (ACRYLAMIDE/HEAT)"],
    "cga_cys_adduct_sida_2024": ["ResearchGate (Ref. 14)"],
    "pmc3199460_ascorbic_dicarbonyl": ["PMC PMCID:PMC3199460 (Ref. 2)"],
    "pmc11889959_spi_tvp_volatiles": ["PMC11889959 (2025)"]
}

# 1. Update benchmark_intake_registry.json
with open(registry_path, "r", encoding="utf-8") as f:
    registry = json.load(f)

for ref in registry.get("eligible_references", []):
    ref_id = ref.get("id")
    if ref_id in aliases_to_add:
        # Add to citation_aliases (ensure unique and non-empty list exists)
        if "citation_aliases" not in ref or ref["citation_aliases"] is None:
            ref["citation_aliases"] = []
        for alias in aliases_to_add[ref_id]:
            if alias not in ref["citation_aliases"]:
                ref["citation_aliases"].append(alias)
                print(f"Added alias '{alias}' to {ref_id} in registry")

    # Fix Matoba
    if ref_id == "matoba_1988_nucleotide_hydrolysis":
        ref["citation"] = "Matoba et al. (1988)"
        ref["doi"] = ""
        print("Standardized matoba_1988_nucleotide_hydrolysis in registry")
    # Fix Cerny
    elif ref_id == "cerny_guntz_dubini_2008":
        ref["citation"] = "Cerny & Guntz-Dubini (2008)"
        ref["doi"] = "10.1021/jf801762c"
        print("Standardized cerny_guntz_dubini_2008 in registry")

with open(registry_path, "w", encoding="utf-8") as f:
    json.dump(registry, f, indent=2, ensure_ascii=False)

# 2. Update slr_incorporation_matrix.json
with open(matrix_path, "r", encoding="utf-8") as f:
    matrix = json.load(f)

for entry in matrix.get("entries", []):
    paper_id = entry.get("paper_id")
    if paper_id == "cerny_guntz_dubini_2008":
        entry["citation"] = "Cerny & Guntz-Dubini (2008)"
        entry["doi"] = "10.1021/jf801762c"
        print("Standardized cerny_guntz_dubini_2008 in matrix")
    elif paper_id == "matoba_1988_nucleotide_hydrolysis":
        entry["citation"] = "Matoba et al. (1988)"
        entry["doi"] = ""
        print("Standardized matoba_1988_nucleotide_hydrolysis in matrix")

with open(matrix_path, "w", encoding="utf-8") as f:
    json.dump(matrix, f, indent=2, ensure_ascii=False)

# 3. Update computational_priors.json
with open(priors_path, "r", encoding="utf-8") as f:
    priors = json.load(f)

def update_recursive(data):
    count = 0
    if isinstance(data, dict):
        if "id" in data:
            if data["id"] == "matoba_1988_nucleotide_hydrolysis" and "source" in data:
                data["source"] = "Matoba et al. (1988)"
                print("Updated matoba_1988_nucleotide_hydrolysis in priors")
                count += 1
            elif data["id"] == "cerny_guntz_dubini_2008" and "source" in data:
                data["source"] = "Cerny & Guntz-Dubini (2008)"
                print("Updated cerny_guntz_dubini_2008 in priors")
                count += 1
        for k, v in data.items():
            count += update_recursive(v)
    elif isinstance(data, list):
        for item in data:
            count += update_recursive(item)
    return count

update_recursive(priors)

with open(priors_path, "w", encoding="utf-8") as f:
    json.dump(priors, f, indent=2, ensure_ascii=False)

# 4. Update process_state_calibrations.json
with open(calibrations_path, "r", encoding="utf-8") as f:
    calibrations = json.load(f)

for entry in calibrations.get("entries", []):
    entry_id = entry.get("id")
    if entry_id == "matoba_1988_nucleotide_hydrolysis_state":
        entry["source_citation"] = "Matoba et al. (1988)"
        entry["doi"] = ""
        print("Standardized matoba_1988_nucleotide_hydrolysis_state in calibrations")
    elif entry_id == "cerny_guntz_dubini_2008_state":
        entry["source_citation"] = "Cerny & Guntz-Dubini (2008)"
        entry["doi"] = "10.1021/jf801762c"
        print("Standardized cerny_guntz_dubini_2008_state in calibrations")

with open(calibrations_path, "w", encoding="utf-8") as f:
    json.dump(calibrations, f, indent=2, ensure_ascii=False)

print("Finished applying alias and standardization mappings.")
