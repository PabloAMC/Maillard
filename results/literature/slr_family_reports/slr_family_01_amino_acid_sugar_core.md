# Systematic Literature Review — Family 1: Amino acid-sugar Maillard core
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering amino acid-sugar maillard core.  
**Strategic Posture:** `first_class_core`  
**Runtime Concept:** `core_reaction_network`  
**Scope & Targets:** Covers target compounds/variables: `2-methyl-3-furanthiol`, `2-furfurylthiol`, `methional`, `2,5-dimethylpyrazine`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`. Target runtime modules: `src/precursor_resolver.py`, `src/literature_runtime.py`, `src/benchmark_validation.py`.  

---

## Acceptance Checklist (8 criteria per paper)

| # | Criterion |
|---|---|
| C1 | Reactant/substrate identities with supplier, purity, or CAS |
| C2 | Reactant concentrations / precursor loads specified |
| C3 | Temperature, time, pH, and water activity (Aw) reported |
| C4 | Analytical method with identified internal standard (IS) |
| C5 | Absolute yields/concentrations (not relative peak areas only) |
| C6 | ≥ 3 replicates or per-compound uncertainty reported |
| C7 | LOD/LOQ or non-detect policy stated |
| C8 | Odor thresholds, TAV, or OAV reported |

**Primary benchmark threshold:** ≥ 6/8 → enters `benchmark_schema.json`  
**Calibration threshold:** 3–5/8 → useful for individual parameter adjustment  
**Rejection:** < 3/8  

---

## SECTION 1 — Active Repository Anchors (from SLR Incorporation Matrix)

### Mottram & Nobrega (2002), JAFC 50:1383
- **DOI:** [10.1021/jf0200826](https://doi.org/10.1021/jf0200826)
- **Matrix Family:** free_model_system
- **Confidence Tier:** medium
- **Incorporation Status:** `modeled_not_payloaded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ❌ | LOD/LOQ or incorporation status confirmed (not verified in matrix entry) |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** 2-methyl-3-furanthiol, 2-furfurylthiol
- **Parameters Supported:** carbon_skeleton_assignment, ribose_fft_suppression
- **Exact Numeric Anchors:** pH 5, 95 C, 4 h
- **Notes on Limits:** *Mechanistic anchor only; not an absolute concentration benchmark.*
- **Next Action:** Keep as mechanistic anchor in the matrix and only add a dedicated payload if reporting needs explicit isotope provenance.

---

### Hofmann & Schieberle (1998), JAFC 46:235
- **DOI:** [10.1021/jf9705983](https://doi.org/10.1021/jf9705983)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_not_yet_scored`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** 2-methyl-3-furanthiol, 2-furfurylthiol
- **Parameters Supported:** relative_pathway_weighting, sugar_family_selectivity
- **Exact Numeric Anchors:** max MFT 1.4 mol% at 180 C
- **Notes on Limits:** *High-confidence free-system calibration, not a protein-matrix benchmark.*
- **Next Action:** Expose as a free-system calibration anchor in benchmark-oriented reporting only.

---

### Nishimura & Abe (2024)
- **DOI:** [10.1016/j.foodchem.2024.141599](https://doi.org/10.1016/j.foodchem.2024.141599)
- **Matrix Family:** soy_hydrolysate
- **Confidence Tier:** medium
- **Incorporation Status:** `encoded_qualitative_only`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ❌ | LOD/LOQ or incorporation status confirmed (not verified in matrix entry) |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** sulfur compounds, pyrazines, methional
- **Parameters Supported:** peptide_bound_cysteine_reactivity, protease_dependence
- **Exact Numeric Anchors:** Cys 16.5 mM, Ribose 16.5 mM, 95 C, 90 min
- **Notes on Limits:** *No absolute ppb values for MFT or FFT.*
- **Next Action:** Retain as qualitative soy-hydrolysate anchor; do not promote to absolute benchmark payloads.

---

### Hofmann & Schieberle (1997/1998), JAFC
- **DOI:** [10.1021/jf970892v](https://doi.org/10.1021/jf970892v)
- **Matrix Family:** cooked_meat
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** 2-methyl-3-furanthiol, 2-furfurylthiol
- **Parameters Supported:** meat_reference_target_band
- **Exact Numeric Anchors:** Beef MFT 7-28 ug/kg, Beef FFT 13-42 ug/kg
- **Notes on Limits:** *Reference anchor only; not a plant-matrix formation calibration.*
- **Next Action:** Use as the main meat-side sulfur anchor for reference reporting.

---

### Hernandez et al. (2023), Molecules 28:3151
- **DOI:** [10.3390/molecules28073151](https://doi.org/10.3390/molecules28073151)
- **Matrix Family:** pbma_vs_meat
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_partially_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** 2-methyl-3-furanthiol, furfural, methylpyrazine, methional, 2-methylbutanal, phenylacetaldehyde, benzaldehyde
- **Parameters Supported:** mft_furfural_ratio, pbma_vs_meat_strecker_gap, pbma_pyrazine_overproduction
- **Exact Numeric Anchors:** GEN MFT 41.09 ng/g, PBMA furfural 987-1093 ng/g, RGB 2-methylbutanal 23.51 ng/g, Beyond methional 8.38 ng/g, Impossible phenylacetaldehyde 18.21 ng/g
- **Notes on Limits:** *Commercial endpoint anchor; not a controlled precursor kinetics dataset.*
- **Next Action:** Promote Strecker and pyrazine sub-anchors into user-visible comparison reporting after Track 2 runtime wiring.

---

### Davidek et al. (2006)
- **DOI:** [10.1021/jf053246m](https://doi.org/10.1021/jf053246m)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** furan, methylfuran
- **Parameters Supported:** safety_furan_yield_ceiling
- **Exact Numeric Anchors:** furan yield 0.045 mol% at 120 C
- **Notes on Limits:** *Free model system, not an extrusion-matrix benchmark.*
- **Next Action:** Use as safety ceiling for furan formation from threonine.

---

### Martins & van Boekel (2003)
- **DOI:** [10.1021/jf021111m](https://doi.org/10.1021/jf021111m)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** furosine, pyrraline
- **Parameters Supported:** furosine_formation_ea, pyrraline_formation_ea
- **Exact Numeric Anchors:** furosine Ea 28.5 kcal/mol, pyrraline Ea 31.2 kcal/mol
- **Notes on Limits:** *Free model system kinetics.*
- **Next Action:** Use as prior bounds for furosine and pyrraline activation energies.

---

### Liu et al. (2023/2024)
- **DOI:** [10.1021/acs.jafc.9b05898](https://doi.org/10.1021/acs.jafc.9b05898)
- **Matrix Family:** free_model_system
- **Confidence Tier:** medium
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** Amadori compound, deoxyosone
- **Parameters Supported:** baseline_arp_formation_ea_kj_per_mol, egcg_arp_formation_ea_kj_per_mol
- **Exact Numeric Anchors:** Ea reduction 15.0 kJ/mol
- **Notes on Limits:** *Free model system, not a plant-protein matrix benchmark.*
- **Next Action:** Use to model EGCG dicarbonyl trapping and ARP formation catalysis.

---

### Mottram & Wedzicha (1990)
- **DOI:** [10.1016/0308-8146(90)90009-Q](https://doi.org/10.1016/0308-8146(90)90009-Q)
- **Matrix Family:** free_model_system
- **Confidence Tier:** medium
- **Incorporation Status:** `modeled_not_payloaded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ❌ | LOD/LOQ or incorporation status confirmed (not verified in matrix entry) |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** methional
- **Parameters Supported:** met_mgo_to_methional_conversion_pct
- **Exact Numeric Anchors:** conversion yield 12%
- **Notes on Limits:** *Free model system kinetics.*
- **Next Action:** Use to model methionine-dicarbonyl Strecker kinetics.

---

### Yu et al. (2021)
- **DOI:** [10.1016/j.lwt.2021.111802](https://doi.org/10.1016/j.lwt.2021.111802)
- **Matrix Family:** free_model_system
- **Confidence Tier:** medium
- **Incorporation Status:** `modeled_not_payloaded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ❌ | LOD/LOQ or incorporation status confirmed (not verified in matrix entry) |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** pyrazines, strecker aldehydes
- **Parameters Supported:** corn_hydrolysate_volatile_yield
- **Exact Numeric Anchors:** absolute volatiles at 100 C / 60 min
- **Notes on Limits:** *Hydrolysate model system.*
- **Next Action:** Use to model corn protein hydrolysate kinetics.

---

### Perdiguero et al. (2004)
- **DOI:** [10.1021/jf0494452](https://doi.org/10.1021/jf0494452)
- **Matrix Family:** yeast_ye
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** GMP
- **Parameters Supported:** yeast_autolysis_kinetics
- **Exact Numeric Anchors:** GMP 2.8 mg/g DW at 55 C, 24h
- **Notes on Limits:** *Yeast cell autolysis matrix.*
- **Next Action:** Expose as a fermentation process state calibration payload.

---

### Bi et al. (2020)
- **DOI:** [10.1021/acs.jafc.9b07711](https://doi.org/10.1021/acs.jafc.9b07711)
- **Matrix Family:** raw_pea_flour
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** hexanal
- **Parameters Supported:** hexanal_baseline_concentration
- **Exact Numeric Anchors:** hexanal 1260 ug/kg
- **Notes on Limits:** *Raw pea flour matrix.*
- **Next Action:** Track raw pea off-note baseline in lipid-Maillard crosstalk model.

---

### Fu et al. (2023)
- **DOI:** [10.3390/foods12101967](https://doi.org/10.3390/foods12101967)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** acrylamide
- **Parameters Supported:** high_temperature_acrylamide_range
- **Exact Numeric Anchors:** acrylamide 31.81–186.70 ug/kg
- **Notes on Limits:** *High-temperature thermal degradation.*
- **Next Action:** Expose as a high-temperature acrylamide safety reference payload.

---

### Poojary et al. (2023)
- **DOI:** [10.1016/j.foodchem.2022.134406](https://doi.org/10.1016/j.foodchem.2022.134406)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** cysteine
- **Parameters Supported:** cga_cys_adduct_stability
- **Exact Numeric Anchors:** stability at 90 C, MRM m/z 474.1
- **Notes on Limits:** *Polyphenol adduct kinetics.*
- **Next Action:** Encode CGA-Cys adduct SIDA stability priors.

---

### Nemet & Monnier (2011)
- **DOI:** [10.1074/jbc.M111.245100](https://doi.org/10.1074/jbc.M111.245100)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** 3-deoxythreosone, xylosone
- **Parameters Supported:** ascorbic_dicarbonyl_yields
- **Exact Numeric Anchors:** 3-deoxythreosone 9.1 pmol/mg, xylosone 0.5 pmol/mg
- **Notes on Limits:** *Ascorbic acid breakdown.*
- **Next Action:** Encode ascorbic dicarbonyl release priors.

---

### Fujimoto et al. (2023)
- **DOI:** [10.1021/acs.jafc.2c08283](https://doi.org/10.1021/acs.jafc.2c08283)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** phosphatidylethanolamine
- **Parameters Supported:** pe_glycation_stoichiometry
- **Exact Numeric Anchors:** 1:2 sugar:PE stoichiometry, absorbance 350 nm
- **Notes on Limits:** *Lipid glycation stoichiometry.*
- **Next Action:** Encode PE glycation stoichiometry priors.

---

### Gigl et al. (2021)
- **DOI:** [10.1021/acs.jafc.1c06163](https://doi.org/10.1021/acs.jafc.1c06163)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** 2-furfurylthiol
- **Parameters Supported:** thiol_melanoidin_depletion
- **Exact Numeric Anchors:** FFT depletion 16-fold, IC50 183 mg/L
- **Notes on Limits:** *Melanoidin radical trapping.*
- **Next Action:** Encode melanoidin thiol trapping priors.

---

### De Vleeschouwer et al. (2006)
- **DOI:** [10.1021/jf051197n](https://doi.org/10.1021/jf051197n)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** acrylamide
- **Parameters Supported:** acrylamide_formation_elimination_kinetics
- **Exact Numeric Anchors:** Ea formation 140.81-168.25 kJ/mol, Ea elimination 148.23-167.21 kJ/mol
- **Notes on Limits:** *Aqueous model system.*
- **Next Action:** Expose as acrylamide aqueous kinetics safety payload.

---

### Zhu et al. (2022)
- **DOI:** [10.3389/fnut.2022.940202](https://doi.org/10.3389/fnut.2022.940202)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** acrylamide, 5-HMF
- **Parameters Supported:** lipid_crosstalk_activation_energy
- **Exact Numeric Anchors:** Ea acrylamide 12.87 kJ/mol, Ea 5-HMF 14.85 kJ/mol
- **Notes on Limits:** *Glucose-asparagine-linoleic acid system.*
- **Next Action:** Expose as lipid-crosstalk acrylamide safety payload.

---

### Sun et al. (2015)
- **DOI:** [10.1016/j.foodchem.2014.09.129](https://doi.org/10.1016/j.foodchem.2014.09.129)
- **Matrix Family:** beef
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, CEL
- **Parameters Supported:** beef_ages_formation_kinetics
- **Exact Numeric Anchors:** Ea CML 61.01 kJ/mol, Ea CEL 29.21 kJ/mol
- **Notes on Limits:** *Pasteurized ground beef.*
- **Next Action:** Expose as beef AGEs safety reference payload.

---

### Zhu et al. (2021)
- **DOI:** [10.1002/jsfa.10528](https://doi.org/10.1002/jsfa.10528)
- **Matrix Family:** chicken
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, CEL
- **Parameters Supported:** free_vs_bound_ages_kinetics
- **Exact Numeric Anchors:** Ea free CML 74.87 kJ/mol, Ea bound CML 68.21 kJ/mol
- **Notes on Limits:** *Braised chicken matrix.*
- **Next Action:** Expose as chicken AGEs safety reference payload.

---

### Hamzalioglu et al. (2026)
- **DOI:** [10.1021/acs.jafc.5c14296](https://doi.org/10.1021/acs.jafc.5c14296)
- **Matrix Family:** milk
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, CEL
- **Parameters Supported:** milk_uht_lactulosyllysine_conversion_kinetics
- **Exact Numeric Anchors:** Ea LacLys formation 52.1 kJ/mol
- **Notes on Limits:** *Ultra-High-Temperature treated milk.*
- **Next Action:** Expose as milk UHT AGEs safety reference payload.

---

### Krause et al. (2003)
- **DOI:** [10.1007/s00217-003-0685-6](https://doi.org/10.1007/s00217-003-0685-6)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** furosine
- **Parameters Supported:** acid_hydrolysis_conversion_yields
- **Exact Numeric Anchors:** fructosyllysine 32% yield (6M HCl), lactulosyllysine 34% yield (6M HCl)
- **Notes on Limits:** *Acid hydrolysis reaction stoichiometry.*
- **Next Action:** Encode furosine hydrolysis yields computational priors.

---

### Hidalgo & Pompei (2000)
- **DOI:** [10.1021/jf990120u](https://doi.org/10.1021/jf990120u)
- **Matrix Family:** tomato
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** furosine
- **Parameters Supported:** tomato_furosine_kinetics
- **Exact Numeric Anchors:** Ea formation 93.9 kJ/mol
- **Notes on Limits:** *Heated tomato products.*
- **Next Action:** Expose as tomato furosine safety reference payload.

---

### Cantre et al. (2007)
- **DOI:** [10.1111/j.1745-4549.2007.00109.x](https://doi.org/10.1111/j.1745-4549.2007.00109.x)
- **Matrix Family:** beef
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** furosine
- **Parameters Supported:** beef_sterilization_furosine_kinetics
- **Exact Numeric Anchors:** Ea formation 95.7 kJ/mol
- **Notes on Limits:** *Corned beef thermal processing.*
- **Next Action:** Expose as beef furosine safety reference payload.

---

### Yu et al. (2017)
- **DOI:** [10.1016/j.tifs.2020.01.021](https://doi.org/10.1016/j.tifs.2020.01.021)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, CEL
- **Parameters Supported:** cml_cel_activation_energy
- **Exact Numeric Anchors:** CML formation Ea 61.01 kJ/mol, CEL formation Ea 29.21 kJ/mol
- **Notes on Limits:** *Complex isolates and meat products review.*
- **Next Action:** Expose as CML/CEL activation energy safety reference payload.

---

### Charissou et al. (2007)
- **DOI:** [10.1021/jf063024j](https://doi.org/10.1021/jf063024j)
- **Matrix Family:** cookie_model
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, furosine
- **Parameters Supported:** cookie_glycation_kinetics
- **Exact Numeric Anchors:** zero-order CML accumulation, transient furosine decay 180 to 220 C
- **Notes on Limits:** *Low moisture cookie model system.*
- **Next Action:** Expose as cookie CML/furosine safety reference payload.

---

### Fratianni et al. (2016)
- **DOI:** [10.1016/j.foodres.2016.12.009](https://doi.org/10.1016/j.foodres.2016.12.009)
- **Matrix Family:** apricot
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** furosine
- **Parameters Supported:** apricot_furosine_kinetics
- **Exact Numeric Anchors:** zero-order furosine formation, Ea 83.3 kJ/mol
- **Notes on Limits:** *Apricot tissue matrix.*
- **Next Action:** Expose as apricot furosine safety reference payload.

---

### Ma et al. (2024)
- **DOI:** [10.3390/ijms25168668](https://doi.org/10.3390/ijms25168668)
- **Matrix Family:** plant_based_meat_analogue
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, CEL, acrylamide
- **Parameters Supported:** extrusion_chemical_damage_kinetics
- **Exact Numeric Anchors:** CML acceleration via glyoxal, CEL 39-101% spike at barrel temp 160 C
- **Notes on Limits:** *High-moisture PBMA twin-screw extrusion.*
- **Next Action:** Expose as PBMA extrusion safety reference payload.

---

### Knol et al. (2005)
- **DOI:** [10.1021/jf050504m](https://doi.org/10.1021/jf050504m)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** acrylamide
- **Parameters Supported:** acrylamide_formation_kinetics, acrylamide_degradation_kinetics
- **Exact Numeric Anchors:** formation Ea 52.1 kJ/mol, degradation Ea 72.9 kJ/mol
- **Notes on Limits:** *Glucose-asparagine model system.*
- **Next Action:** Encode as acrylamide safety reference payload.

---

### De Vleeschouwer et al. (2008)
- **DOI:** [10.1021/bp060389f](https://doi.org/10.1021/bp060389f)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** acrylamide
- **Parameters Supported:** acrylamide_moisture_dependent_kinetics
- **Exact Numeric Anchors:** minimum elimination rate at Aw 0.82
- **Notes on Limits:** *Low-moisture starch systems.*
- **Next Action:** Encode as water activity dependent acrylamide safety reference payload.

---

### Ishak et al. (2022)
- **DOI:** [10.1016/j.foodchem.2022.132372](https://doi.org/10.1016/j.foodchem.2022.132372)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** PhIP, HCAs
- **Parameters Supported:** phip_formation_kinetics, activation_energy
- **Exact Numeric Anchors:** Ea Phe model: 95.36 kJ/mol, Ea Pro model: 114.12 kJ/mol
- **Notes on Limits:** *Creatinine-phenylalanine/proline models.*
- **Next Action:** Encode as PhIP HCA safety reference payload.

---

### Shu et al. (2019)
- **DOI:** [10.1016/j.freeradbiomed.2019.04.026](https://doi.org/10.1016/j.freeradbiomed.2019.04.026)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** cysteine-quinone Michael adducts
- **Parameters Supported:** quinone_thiol_conjugation_kinetics, rate_constant
- **Exact Numeric Anchors:** second-order rate constant k: 10^5 to 10^6 M-1 s-1
- **Notes on Limits:** *Recombinant protein models.*
- **Next Action:** Encode cysteine-quinone Michael adduct kinetics.

---

### Hidalgo et al. (2007)
- **DOI:** [10.1021/jf070527w](https://doi.org/10.1021/jf070527w)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** styrene, 2-pentylpyridine, acrylamide, met-leghemoglobin
- **Parameters Supported:** decadienal_amino_crosstalk_kinetics, leghemoglobin_oxygen_dissociation
- **Exact Numeric Anchors:** Ea styrene: 150.4 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode lipid-amino crosstalk prior.

---

### Zamora et al. (2010)
- **DOI:** [10.1021/jf102026c](https://doi.org/10.1021/jf102026c)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** styrene, 2-pentylpyridine, acrylamide, met-leghemoglobin
- **Parameters Supported:** decadienal_amino_crosstalk_kinetics, leghemoglobin_oxygen_dissociation
- **Exact Numeric Anchors:** Ea asparagine: 81.0 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode lipid-catalyzed decarboxylation prior.

---

### Ding et al. (2020)
- **DOI:** [10.1021/acs.jafc.0c04738](https://doi.org/10.1021/acs.jafc.0c04738)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** styrene, 2-pentylpyridine, acrylamide, met-leghemoglobin
- **Parameters Supported:** decadienal_amino_crosstalk_kinetics, leghemoglobin_oxygen_dissociation
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode micellar Schiff base prior.

---

### Richards et al. (2009)
- **DOI:** [10.1021/jf9013394](https://doi.org/10.1021/jf9013394)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** styrene, 2-pentylpyridine, acrylamide, met-leghemoglobin
- **Parameters Supported:** decadienal_amino_crosstalk_kinetics, leghemoglobin_oxygen_dissociation
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode Hb-mediated oxidation prior.

---

### Smagghe et al. (2006)
- **DOI:** [10.1021/bi051902l](https://doi.org/10.1021/bi051902l)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** styrene, 2-pentylpyridine, acrylamide, met-leghemoglobin
- **Parameters Supported:** decadienal_amino_crosstalk_kinetics, leghemoglobin_oxygen_dissociation
- **Exact Numeric Anchors:** dissociation rate k: 5.6 s-1
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode leghemoglobin dissociation prior.

---

### Urugo et al. (2024)
- **DOI:** [10.3390/foods14193295](https://doi.org/10.3390/foods14193295)
- **Matrix Family:** plant_protein_matrix
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** HCAs, CML, CEL, acrylamide
- **Parameters Supported:** HCA_formation_kinetics, CML_CEL_accumulation_kinetics, browning_activation_energy
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode HCA safety reference prior.

---

### Chen et al. (2016)
- **DOI:** [10.1007/s10068-016-0185-5](https://doi.org/10.1007/s10068-016-0185-5)
- **Matrix Family:** high_moisture_melt
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** HCAs, CML, CEL, acrylamide
- **Parameters Supported:** HCA_formation_kinetics, CML_CEL_accumulation_kinetics, browning_activation_energy
- **Exact Numeric Anchors:** glucose doping acceleration: 4.0x
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode sterilization AGE safety prior.

---

### Pruteanu et al. (2023)
- **DOI:** [10.1016/j.lwt.2024.117316](https://doi.org/10.1016/j.lwt.2024.117316)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** HCAs, CML, CEL, acrylamide
- **Parameters Supported:** HCA_formation_kinetics, CML_CEL_accumulation_kinetics, browning_activation_energy
- **Exact Numeric Anchors:** Ea Phe model: 145 kJ/mol, Ea Gly model: 109 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode browning activation energy safety prior.

---

### Zhu et al. (2020)
- **DOI:** [10.1021/acs.jafc.0c01761](https://doi.org/10.1021/acs.jafc.0c01761)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Exact Numeric Anchors:** MGO k: 1.6 M-1 s-1, GO k: 0.059 M-1 s-1
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode epicatechin dicarbonyl trapping prior.

---

### Zhu et al. (2020)
- **DOI:** [10.1016/j.foodchem.2020.126500](https://doi.org/10.1016/j.foodchem.2020.126500)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode structural polyphenol trapping prior.

---

### Liu et al. (2022)
- **DOI:** [10.1016/j.foodres.2022.112187](https://doi.org/10.1016/j.foodres.2022.112187)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Exact Numeric Anchors:** Ea 4MBQ-Lys: 19.00 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode benzoquinone-lysine conjugation prior.

---

### Cömert & Gökmen (2019)
- **DOI:** [10.1016/j.foodres.2019.03.046](https://doi.org/10.1016/j.foodres.2019.03.046)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode digestive dicarbonyl scavenging prior.

---

### Munoz et al. (2007)
- **DOI:** [10.1021/jf062081+](https://doi.org/10.1021/jf062081+)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Exact Numeric Anchors:** k CQA_Q: 2.73 M-1 s-1
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode chlorogenic o-quinone prior.

---

### Cömert & Gökmen (2021)
- **DOI:** [10.1016/j.foodchem.2020.128670](https://doi.org/10.1016/j.foodchem.2020.128670)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode EC-Cys synergistic scavenging prior.

---

### Song et al. (2009)
- **DOI:** [10.1073/pnas.0810352106](https://doi.org/10.1073/pnas.0810352106)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Exact Numeric Anchors:** k forward: 1547 M-1 s-1
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode quinone-GSH conjugation prior.

---

### Li et al. (2017)
- **DOI:** [10.1021/acs.jafc.6b05811](https://doi.org/10.1021/acs.jafc.6b05811)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode quercetin MGO adduction prior.

---

### Monforte et al. (2018)
- **DOI:** [10.1021/acs.jafc.7b00264](https://doi.org/10.1021/acs.jafc.7b00264)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** epicatechin-adducts, quinone-conjugates, phenylacetaldehyde
- **Parameters Supported:** quinone_Michael_addition_kinetics, dicarbonyl_trapping_rates, Strecker_catalysis
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode polyphenol-catalyzed Strecker prior.

---

### Solís-Calero et al. (2015)
- **DOI:** [10.1039/C4CP05360E](https://doi.org/10.1039/C4CP05360E)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Exact Numeric Anchors:** Ea dehydration: 17.50 kcal/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode phospholipid surface dehydration catalyst prior.

---

### Solís-Calero et al. (2013)
- **DOI:** [10.1021/jp401488j](https://doi.org/10.1021/jp401488j)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Exact Numeric Anchors:** Ea condensation: 8.76 kcal/mol, Ea enaminol: 16.78 kcal/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode PE surface sugar condensation prior.

---

### Lertsiri et al. (1998)
- **DOI:** [10.1271/bbb.62.893](https://doi.org/10.1271/bbb.62.893)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode Amadori-PE accumulation prior.

---

### Hidalgo et al. (2005)
- **DOI:** [10.1007/s00217-004-1065-x](https://doi.org/10.1007/s00217-004-1065-x)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Exact Numeric Anchors:** Ea browning: 66.5 kJ/mol, Ea fluorescence: 50.0 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode PE-ribose-lysine global browning prior.

---

### Zamora et al. (2020)
- **DOI:** [10.3390/molecules25020373](https://doi.org/10.3390/molecules25020373)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode dihydropyridine adduct radical scavenging prior.

---

### Hidalgo et al. (2006)
- **DOI:** [10.1021/jf060848s](https://doi.org/10.1021/jf060848s)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode PE-lysine Rancimat oxidative synergy prior.

---

### Vilanova et al. (2012)
- **DOI:** [10.1021/jp2116033](https://doi.org/10.1021/jp2116033)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Exact Numeric Anchors:** Ea dehydration: 13.08 kcal/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode Schiff base dehydration rate prior.

---

### Biondi et al. (2010)
- **DOI:** [10.1016/j.lwt.2010.02.016](https://doi.org/10.1016/j.lwt.2010.02.016)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** phosphatidylethanolamine, glyoxal, Amadori-PE, dihydropyridines
- **Parameters Supported:** glycation_kinetics, dehydration_barrier, radical_scavenging
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode microwave thermal stress lipid oxidation prior.

---

### Brands et al. (2002)
- **DOI:** [10.1021/jf010789c](https://doi.org/10.1021/jf010789c)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Exact Numeric Anchors:** Ea polymerization: 128 kJ/mol, rate k: 1.14 h-1
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode casein-sugar melanoidin polymerization prior.

---

### Gigl et al. (2021)
- **DOI:** [10.1021/acs.jafc.1c06163](https://doi.org/10.1021/acs.jafc.1c06163)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode melanoidin covalent thiol staling prior.

---

### Suzuki & Philp (1990)
- **DOI:** [10.1016/0146-6380(90)90162-s](https://doi.org/10.1016/0146-6380(90)90162-s)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode sulfur-incorporated melanoidin prior.

---

### Mundt & Wedzicha (2007)
- **DOI:** [10.1016/j.lwt.2006.07.014](https://doi.org/10.1016/j.lwt.2006.07.014)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Exact Numeric Anchors:** Ea browning: 105 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode biscuit surface browning rate prior.

---

### Cao et al. (2024)
- **DOI:** [10.1111/1750-3841.17378](https://doi.org/10.1111/1750-3841.17378)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode myoglobin-MRP lipid oxidation protection prior.

---

### Hofmann et al. (2001)
- **DOI:** [10.1021/jf001302l](https://doi.org/10.1021/jf001302l)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode irreversible melanoidin-thioether staling prior.

---

### Adams et al. (2005)
- **DOI:** [10.1021/jf047903m](https://doi.org/10.1021/jf047903m)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** melanoidins, thiols, pyrazines
- **Parameters Supported:** polymerization_kinetics, covalent_binding, aroma_staling
- **Exact Numeric Anchors:** Ea histidine: 35.31 kJ/mol, Ea lysine: 54.94 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode ascorbic-basic-amino browning prior.

---

### Smuda & Glomb (2013)
- **DOI:** [10.1002/anie.201300399](https://doi.org/10.1002/anie.201300399)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** oxidative alpha-fragmentation: 31%, beta-cleavage: 32%
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode vitamin C degradation mass balance prior.

---

### Serpen & Gökmen (2007)
- **DOI:** [10.1016/j.foodchem.2006.11.073](https://doi.org/10.1016/j.foodchem.2006.11.073)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** copper acceleration: 88-fold, iron acceleration: 14-fold
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode metal-catalyzed ascorbic acid oxidation prior.

---

### Yang et al. (2021)
- **DOI:** [10.47836/ifrj.28.3.16](https://doi.org/10.47836/ifrj.28.3.16)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** Ea (4:1 AA:Gly): 60.76 kJ/mol, Ea (1:4 AA:Gly): 70.16 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode AA-glycine browning activation energy prior.

---

### Yu et al. (2018)
- **DOI:** [10.1590/1678-457x.08717](https://doi.org/10.1590/1678-457x.08717)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** Ea lysine: 54.94 kJ/mol, Ea arginine: 50.08 kJ/mol, Ea histidine: 35.31 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode basic amino acid AA-browning prior.

---

### Manso et al. (2001)
- **DOI:** [10.1046/j.1365-2621.2001.t01-1-00460.x](https://doi.org/10.1046/j.1365-2621.2001.t01-1-00460.x)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** Ea degradation: 55.6 kJ/mol
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode orange juice ascorbic degradation Weibull prior.

---

### Takase et al. (2025)
- **DOI:** [10.1021/acsfoodscitech.4c00956](https://doi.org/10.1021/acsfoodscitech.4c00956)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode anaerobic DHAA degradation prior.

---

### Hendrickx et al. (1998)
- **DOI:** [10.1021/jf9708251](https://doi.org/10.1021/jf9708251)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Exact Numeric Anchors:** Ea degradation: 54.8 kJ/mol, D-value 121C: 246 min
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode squeezed tomato isobaric degradation prior.

---

### Jian et al. (2012)
- **DOI:** [10.1021/jf3032342](https://doi.org/10.1021/jf3032342)
- **Matrix Family:** free_model_system
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations / loads specified (not verified in matrix entry) |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 4/8 → Calibration reference**

- **Compounds Supported:** ascorbic_acid, dehydroascorbic_acid, CML, xylosone, threose
- **Parameters Supported:** degradation_kinetics, mass_balance, browning_activation_energy
- **Notes on Limits:** *Surrogate literature transfer limits.*
- **Next Action:** Encode ethanolic ascorbic degradation prior.

---

### Wehrmaker et al. (2022)
- **DOI:** [10.1021/acsfoodscitech.2c00242](https://doi.org/10.1021/acsfoodscitech.2c00242)
- **Matrix Family:** pbma
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** lysine
- **Parameters Supported:** lysine_loss_extrusion
- **Exact Numeric Anchors:** Lys loss 22%
- **Notes on Limits:** *See safety_reference_payloads.json.*
- **Next Action:** None — fully encoded.

---

### Ma et al. (2024)
- **DOI:** [10.3390/ijms25168668](https://doi.org/10.3390/ijms25168668)
- **Matrix Family:** spi_heated
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** acrylamide
- **Parameters Supported:** acrylamide_spi
- **Exact Numeric Anchors:** acrylamide 112 ug/kg at 150°C
- **Notes on Limits:** *See safety_reference_payloads.json.*
- **Next Action:** None — fully encoded.

---

### Fu et al. (2023)
- **DOI:** [10.3390/foods12101967](https://doi.org/10.3390/foods12101967)
- **Matrix Family:** pbma
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** CML
- **Parameters Supported:** cml_pbma_proxy
- **Exact Numeric Anchors:** CML 38 ug/g protein ELISA
- **Notes on Limits:** *See safety_reference_payloads.json.*
- **Next Action:** None — fully encoded.

---

### Ramírez-Jiménez et al. (2000)
- **DOI:** [10.1021/jf9907687](https://doi.org/10.1021/jf9907687)
- **Matrix Family:** toasted_bread
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** furosine
- **Parameters Supported:** furosine_bread_benchmark
- **Exact Numeric Anchors:** furosine 1250 mg/100g protein toast
- **Notes on Limits:** *See safety_reference_payloads.json.*
- **Next Action:** None — fully encoded.

---

### Feng et al. (2023)
- **DOI:** [10.3389/fnut.2022.1022254](https://doi.org/10.3389/fnut.2022.1022254)
- **Matrix Family:** hcw_model
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** ascorbic acid
- **Parameters Supported:** hcw_aa_degradation
- **Exact Numeric Anchors:** Ea 82 kJ/mol HCW
- **Notes on Limits:** *See target artifact.*
- **Next Action:** None — fully encoded.

---

### Grandhee & Monnier (1991)
- **DOI:** [10.1074/jbc.266.18.11644](https://doi.org/10.1074/jbc.266.18.11644)
- **Matrix Family:** protein_model
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 5/8 → Calibration reference**

- **Compounds Supported:** pentosidine
- **Parameters Supported:** pentosidine_aa_equivalence
- **Exact Numeric Anchors:** pentosidine 13.2-17.0 pmol/mg
- **Notes on Limits:** *See target artifact.*
- **Next Action:** None — fully encoded.

---

### Yu et al. (2018)
- **DOI:** [10.1590/1678-457X.08717](https://doi.org/10.1590/1678-457X.08717)
- **Matrix Family:** protein_model
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ❌ | Absolute yields reported (not verified in matrix entry) |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 6/8 → Benchmark-eligible**

- **Compounds Supported:** CML, crosslinks
- **Parameters Supported:** aa_crosslink_amino_acid_hierarchy
- **Exact Numeric Anchors:** Lys Ea 54.94, Arg 50.08, His 35.31 kJ/mol
- **Notes on Limits:** *See target artifact.*
- **Next Action:** None — fully encoded.

---

### Marquez-Ruiz et al. (2014)
- **DOI:** [10.1021/jf502636m](https://doi.org/10.1021/jf502636m)
- **Matrix Family:** oleic_model
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Conditions T, pH, time specified (not verified in matrix entry) |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ✅ | Odor thresholds / sensory parameters reported |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** nonanal
- **Parameters Supported:** nonanal_oav_anchor
- **Exact Numeric Anchors:** nonanal ODT 1.0 ppb, OAV 1-28
- **Notes on Limits:** *See target artifact.*
- **Next Action:** None — fully encoded.

---

### Li et al. (2026)
- **DOI:** [10.3390/foods15050912](https://doi.org/10.3390/foods15050912)
- **Matrix Family:** hme_matrix
- **Confidence Tier:** high
- **Incorporation Status:** `encoded_shown`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output compounds specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates / confidence tier verified |
| ✅ | LOD/LOQ or incorporation status confirmed |
| ❌ | Odor thresholds / sensory parameters reported (not verified in matrix entry) |

**Score: 7/8 → Benchmark-eligible**

- **Compounds Supported:** hexanal
- **Parameters Supported:** hme_hexanal_baseline
- **Exact Numeric Anchors:** hexanal 1.8 ug/g HME
- **Notes on Limits:** *See target artifact.*
- **Next Action:** None — fully encoded.

---

## SECTION 2 — Backlog Research Candidates (from Deep Research Backlog)

*No backlog research candidates found for Family 01.*
## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Amino acid-sugar Maillard core` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `2-methyl-3-furanthiol`, `2-furfurylthiol`, `methional`, `2,5-dimethylpyrazine` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 1

Total papers analyzed: **82** (Benchmark-eligible: **46**, Calibration references: **36**, Rejected: **0**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Mottram & Nobrega (2002), JAFC 50:1383 | 5/8 | ⚠️ Calibration | free_model_system |
| Hofmann & Schieberle (1998), JAFC 46:235 | 7/8 | ✅ Eligible | free_model_system |
| Nishimura & Abe (2024) | 5/8 | ⚠️ Calibration | soy_hydrolysate |
| Hofmann & Schieberle (1997/1998), JAFC | 6/8 | ✅ Eligible | cooked_meat |
| Hernandez et al. (2023), Molecules 28:3151 | 7/8 | ✅ Eligible | pbma_vs_meat |
| Davidek et al. (2006) | 7/8 | ✅ Eligible | free_model_system |
| Martins & van Boekel (2003) | 6/8 | ✅ Eligible | free_model_system |
| Liu et al. (2023/2024) | 6/8 | ✅ Eligible | free_model_system |
| Mottram & Wedzicha (1990) | 5/8 | ⚠️ Calibration | free_model_system |
| Yu et al. (2021) | 5/8 | ⚠️ Calibration | free_model_system |
| Perdiguero et al. (2004) | 7/8 | ✅ Eligible | yeast_ye |
| Bi et al. (2020) | 7/8 | ✅ Eligible | raw_pea_flour |
| Fu et al. (2023) | 7/8 | ✅ Eligible | free_model_system |
| Poojary et al. (2023) | 6/8 | ✅ Eligible | free_model_system |
| Nemet & Monnier (2011) | 6/8 | ✅ Eligible | free_model_system |
| Fujimoto et al. (2023) | 6/8 | ✅ Eligible | free_model_system |
| Gigl et al. (2021) | 7/8 | ✅ Eligible | free_model_system |
| De Vleeschouwer et al. (2006) | 6/8 | ✅ Eligible | free_model_system |
| Zhu et al. (2022) | 6/8 | ✅ Eligible | free_model_system |
| Sun et al. (2015) | 6/8 | ✅ Eligible | beef |
| Zhu et al. (2021) | 6/8 | ✅ Eligible | chicken |
| Hamzalioglu et al. (2026) | 6/8 | ✅ Eligible | milk |
| Krause et al. (2003) | 6/8 | ✅ Eligible | free_model_system |
| Hidalgo & Pompei (2000) | 5/8 | ⚠️ Calibration | tomato |
| Cantre et al. (2007) | 5/8 | ⚠️ Calibration | beef |
| Yu et al. (2017) | 6/8 | ✅ Eligible | free_model_system |
| Charissou et al. (2007) | 6/8 | ✅ Eligible | cookie_model |
| Fratianni et al. (2016) | 5/8 | ⚠️ Calibration | apricot |
| Ma et al. (2024) | 6/8 | ✅ Eligible | plant_based_meat_analogue |
| Knol et al. (2005) | 5/8 | ⚠️ Calibration | free_model_system |
| De Vleeschouwer et al. (2008) | 6/8 | ✅ Eligible | free_model_system |
| Ishak et al. (2022) | 6/8 | ✅ Eligible | free_model_system |
| Shu et al. (2019) | 6/8 | ✅ Eligible | free_model_system |
| Hidalgo et al. (2007) | 5/8 | ⚠️ Calibration | free_model_system |
| Zamora et al. (2010) | 5/8 | ⚠️ Calibration | free_model_system |
| Ding et al. (2020) | 4/8 | ⚠️ Calibration | free_model_system |
| Richards et al. (2009) | 4/8 | ⚠️ Calibration | free_model_system |
| Smagghe et al. (2006) | 6/8 | ✅ Eligible | free_model_system |
| Urugo et al. (2024) | 4/8 | ⚠️ Calibration | plant_protein_matrix |
| Chen et al. (2016) | 6/8 | ✅ Eligible | high_moisture_melt |
| Pruteanu et al. (2023) | 6/8 | ✅ Eligible | free_model_system |
| Zhu et al. (2020) | 5/8 | ⚠️ Calibration | free_model_system |
| Zhu et al. (2020) | 4/8 | ⚠️ Calibration | free_model_system |
| Liu et al. (2022) | 5/8 | ⚠️ Calibration | free_model_system |
| Cömert & Gökmen (2019) | 4/8 | ⚠️ Calibration | free_model_system |
| Munoz et al. (2007) | 6/8 | ✅ Eligible | free_model_system |
| Cömert & Gökmen (2021) | 4/8 | ⚠️ Calibration | free_model_system |
| Song et al. (2009) | 5/8 | ⚠️ Calibration | free_model_system |
| Li et al. (2017) | 4/8 | ⚠️ Calibration | free_model_system |
| Monforte et al. (2018) | 4/8 | ⚠️ Calibration | free_model_system |
| Solís-Calero et al. (2015) | 6/8 | ✅ Eligible | free_model_system |
| Solís-Calero et al. (2013) | 6/8 | ✅ Eligible | free_model_system |
| Lertsiri et al. (1998) | 4/8 | ⚠️ Calibration | free_model_system |
| Hidalgo et al. (2005) | 6/8 | ✅ Eligible | free_model_system |
| Zamora et al. (2020) | 4/8 | ⚠️ Calibration | free_model_system |
| Hidalgo et al. (2006) | 4/8 | ⚠️ Calibration | free_model_system |
| Vilanova et al. (2012) | 6/8 | ✅ Eligible | free_model_system |
| Biondi et al. (2010) | 4/8 | ⚠️ Calibration | free_model_system |
| Brands et al. (2002) | 6/8 | ✅ Eligible | free_model_system |
| Gigl et al. (2021) | 4/8 | ⚠️ Calibration | free_model_system |
| Suzuki & Philp (1990) | 4/8 | ⚠️ Calibration | free_model_system |
| Mundt & Wedzicha (2007) | 5/8 | ⚠️ Calibration | free_model_system |
| Cao et al. (2024) | 4/8 | ⚠️ Calibration | free_model_system |
| Hofmann et al. (2001) | 4/8 | ⚠️ Calibration | free_model_system |
| Adams et al. (2005) | 6/8 | ✅ Eligible | free_model_system |
| Smuda & Glomb (2013) | 6/8 | ✅ Eligible | free_model_system |
| Serpen & Gökmen (2007) | 6/8 | ✅ Eligible | free_model_system |
| Yang et al. (2021) | 5/8 | ⚠️ Calibration | free_model_system |
| Yu et al. (2018) | 6/8 | ✅ Eligible | free_model_system |
| Manso et al. (2001) | 5/8 | ⚠️ Calibration | free_model_system |
| Takase et al. (2025) | 4/8 | ⚠️ Calibration | free_model_system |
| Hendrickx et al. (1998) | 6/8 | ✅ Eligible | free_model_system |
| Jian et al. (2012) | 4/8 | ⚠️ Calibration | free_model_system |
| Wehrmaker et al. (2022) | 5/8 | ⚠️ Calibration | pbma |
| Ma et al. (2024) | 7/8 | ✅ Eligible | spi_heated |
| Fu et al. (2023) | 7/8 | ✅ Eligible | pbma |
| Ramírez-Jiménez et al. (2000) | 6/8 | ✅ Eligible | toasted_bread |
| Feng et al. (2023) | 6/8 | ✅ Eligible | hcw_model |
| Grandhee & Monnier (1991) | 5/8 | ⚠️ Calibration | protein_model |
| Yu et al. (2018) | 6/8 | ✅ Eligible | protein_model |
| Marquez-Ruiz et al. (2014) | 7/8 | ✅ Eligible | oleic_model |
| Li et al. (2026) | 7/8 | ✅ Eligible | hme_matrix |

---

## Consolidated entries for benchmark_schema.json — Family 1

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Hofmann & Schieberle (1998), JAFC 46:235` (Score: 7/8)
- `Hofmann & Schieberle (1997/1998), JAFC` (Score: 6/8)
- `Hernandez et al. (2023), Molecules 28:3151` (Score: 7/8)
- `Davidek et al. (2006)` (Score: 7/8)
- `Martins & van Boekel (2003)` (Score: 6/8)
- `Liu et al. (2023/2024)` (Score: 6/8)
- `Perdiguero et al. (2004)` (Score: 7/8)
- `Bi et al. (2020)` (Score: 7/8)
- `Fu et al. (2023)` (Score: 7/8)
- `Poojary et al. (2023)` (Score: 6/8)
- `Nemet & Monnier (2011)` (Score: 6/8)
- `Fujimoto et al. (2023)` (Score: 6/8)
- `Gigl et al. (2021)` (Score: 7/8)
- `De Vleeschouwer et al. (2006)` (Score: 6/8)
- `Zhu et al. (2022)` (Score: 6/8)
- `Sun et al. (2015)` (Score: 6/8)
- `Zhu et al. (2021)` (Score: 6/8)
- `Hamzalioglu et al. (2026)` (Score: 6/8)
- `Krause et al. (2003)` (Score: 6/8)
- `Yu et al. (2017)` (Score: 6/8)
- `Charissou et al. (2007)` (Score: 6/8)
- `Ma et al. (2024)` (Score: 6/8)
- `De Vleeschouwer et al. (2008)` (Score: 6/8)
- `Ishak et al. (2022)` (Score: 6/8)
- `Shu et al. (2019)` (Score: 6/8)
- `Smagghe et al. (2006)` (Score: 6/8)
- `Chen et al. (2016)` (Score: 6/8)
- `Pruteanu et al. (2023)` (Score: 6/8)
- `Munoz et al. (2007)` (Score: 6/8)
- `Solís-Calero et al. (2015)` (Score: 6/8)
- `Solís-Calero et al. (2013)` (Score: 6/8)
- `Hidalgo et al. (2005)` (Score: 6/8)
- `Vilanova et al. (2012)` (Score: 6/8)
- `Brands et al. (2002)` (Score: 6/8)
- `Adams et al. (2005)` (Score: 6/8)
- `Smuda & Glomb (2013)` (Score: 6/8)
- `Serpen & Gökmen (2007)` (Score: 6/8)
- `Yu et al. (2018)` (Score: 6/8)
- `Hendrickx et al. (1998)` (Score: 6/8)
- `Ma et al. (2024)` (Score: 7/8)
- `Fu et al. (2023)` (Score: 7/8)
- `Ramírez-Jiménez et al. (2000)` (Score: 6/8)
- `Feng et al. (2023)` (Score: 6/8)
- `Yu et al. (2018)` (Score: 6/8)
- `Marquez-Ruiz et al. (2014)` (Score: 7/8)
- `Li et al. (2026)` (Score: 7/8)

### CALIBRATION (Score 3-5)
- `Mottram & Nobrega (2002), JAFC 50:1383` (Score: 5/8)
- `Nishimura & Abe (2024)` (Score: 5/8)
- `Mottram & Wedzicha (1990)` (Score: 5/8)
- `Yu et al. (2021)` (Score: 5/8)
- `Hidalgo & Pompei (2000)` (Score: 5/8)
- `Cantre et al. (2007)` (Score: 5/8)
- `Fratianni et al. (2016)` (Score: 5/8)
- `Knol et al. (2005)` (Score: 5/8)
- `Hidalgo et al. (2007)` (Score: 5/8)
- `Zamora et al. (2010)` (Score: 5/8)
- `Ding et al. (2020)` (Score: 4/8)
- `Richards et al. (2009)` (Score: 4/8)
- `Urugo et al. (2024)` (Score: 4/8)
- `Zhu et al. (2020)` (Score: 5/8)
- `Zhu et al. (2020)` (Score: 4/8)
- `Liu et al. (2022)` (Score: 5/8)
- `Cömert & Gökmen (2019)` (Score: 4/8)
- `Cömert & Gökmen (2021)` (Score: 4/8)
- `Song et al. (2009)` (Score: 5/8)
- `Li et al. (2017)` (Score: 4/8)
- `Monforte et al. (2018)` (Score: 4/8)
- `Lertsiri et al. (1998)` (Score: 4/8)
- `Zamora et al. (2020)` (Score: 4/8)
- `Hidalgo et al. (2006)` (Score: 4/8)
- `Biondi et al. (2010)` (Score: 4/8)
- `Gigl et al. (2021)` (Score: 4/8)
- `Suzuki & Philp (1990)` (Score: 4/8)
- `Mundt & Wedzicha (2007)` (Score: 5/8)
- `Cao et al. (2024)` (Score: 4/8)
- `Hofmann et al. (2001)` (Score: 4/8)
- `Yang et al. (2021)` (Score: 5/8)
- `Manso et al. (2001)` (Score: 5/8)
- `Takase et al. (2025)` (Score: 4/8)
- `Jian et al. (2012)` (Score: 4/8)
- `Wehrmaker et al. (2022)` (Score: 5/8)
- `Grandhee & Monnier (1991)` (Score: 5/8)

### REJECTED (Score < 3)
- *No rejected references.*

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/precursor_resolver.py`, `src/literature_runtime.py`, `src/benchmark_validation.py` must explicitly account for `core_reaction_network` as a modifier when predicting `2-methyl-3-furanthiol`, `2-furfurylthiol`, `methional`, `2,5-dimethylpyrazine` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `first_class_core` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
