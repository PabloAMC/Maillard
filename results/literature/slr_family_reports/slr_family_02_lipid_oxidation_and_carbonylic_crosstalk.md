# Systematic Literature Review — Family 2: Lipid oxidation and carbonylic crosstalk
**Last updated:** 2026-10-27  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering lipid oxidation and carbonylic crosstalk.  
**Strategic Posture:** `immediate_expansion_lane`  
**Runtime Concept:** `lipid_crosstalk_lane`  
**Scope & Targets:** Covers target compounds/variables: `hexanal`, `2-pentylfuran`, `nonanal`, `carbonyl_competition_factor`. Preferred payload types: `benchmark_payload`, `retention_payload`, `computational_prior`. Target runtime modules: `src/lipid_oxidation.py`.  

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

## SECTION 1 — Curated References

### Trikusuma et al. (2020)
- **DOI:** [10.1016/j.foodchem.2019.126082](https://doi.org/10.1016/j.foodchem.2019.126082)
- **Matrix Family:** pea_uht_beverage
- **Kind:** `conditional_calibration`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ✅ | Analytical method / output mode specified |
| ✅ | Absolute yields reported |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 6/8 → Benchmark-eligible**

**Key Values:**
| Parameter | Value |
|---|---|
| `pea_protein_wt_pct` | `3.0` |
| `pea_protein_fat_wt_pct` | `8.0` |
| `pea_protein_sugar_wt_pct` | `2.0` |
| `homogenization_mpa` | `17.2` |
| `preheat_temp_C` | `80.0` |
| `uht_temp_C` | `140.0` |
| `uht_hold_seconds` | `6.0` |
| `storage_temp_C` | `5.0` |
| `storage_duration_weeks` | `7` |
| `volatile_output_mode` | `absolute_ug_per_l` |
| `tracked_uht_markers_ug_per_l` | `{'hexanal': 782.0, '2-pentylfuran': 163.0, 'nonanal': 24.0, 'methional': 3.1, '2,5-dimethylpyrazine': 2.29}` |
| `replicates` | `3` |

**What it supports:**
- Quantitative UHT pea-protein aroma panel with 21 DHS-GC-MS-MS quantified compounds
- Mixed Maillard and lipid-oxidation context in a real heated pea-protein beverage
- Executable matrix-only validation subset for heated pea off-flavour markers

**What it does not support:**
- Strict-gate free-precursor validation
- Direct MFT or FFT calibration in pea matrices
- A full meaty-positive matrix benchmark

**Next Action:** Use the extracted UHT quantitative table as a heated pea matrix-only benchmark subset while keeping the full 21-compound panel available as intake evidence.

---

### Lin, Fay, Welti & Blank (2001), Lipids 36:749
- **DOI:** [10.1007/s11745-001-0781-x](https://doi.org/10.1007/s11745-001-0781-x)
- **Matrix Family:** lipid_rich_model_system
- **Kind:** `directional_prior`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_directional_prior_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `target_offnote` | `trans-4,5-epoxy-(E)-2-decenal` |
| `odor_threshold_ng_per_l` | `0.07` |
| `precursor_fatty_acid` | `arachidonic_acid` |

**What it supports:**
- An ultra-low-threshold epoxydecenal guardrail for arachidonic-acid branches
- Bounded off-note reporting support for highly unsaturated lipid matrices

**What it does not support:**
- Matrix-specific concentration closure
- A full oxidation mechanism

**Next Action:** Use as a lipid off-note guardrail in bounded runtime and reporting lanes.

---

### Hrncirik & Zeelenberg (2013), J. Am. Oil Chem. Soc. 91:591
- **DOI:** [10.1007/s11746-013-2401-2](https://doi.org/10.1007/s11746-013-2401-2)
- **Matrix Family:** coconut_oil_co_matrix
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_calibration_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ✅ | LOD/LOQ stated |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `journal_volume_page` | `62:9666` |
| `delta_decalactone_oav_reported` | `True` |
| `comparative_to_pufa_oils` | `True` |

**What it supports:**
- Thermal-stability context for coconut-oil-rich co-matrices against PUFA-rich oils
- Bounded delta-decalactone support for lipid-rich aroma interpretation

**What it does not support:**
- Exact plant-matrix volatile concentration closure
- A full lipid degradation mechanism outside coconut-oil-like systems

**Next Action:** Use as a coconut-oil co-matrix thermal-profile calibration so lipid-rich formulations are not interpreted with PUFA-only oxidation context.

---

### Belitz, Grosch & Schieberle (2009), Food Chemistry, 4th rev. ed., Springer (aroma-threshold tables, Ch. 5)
- **DOI:** [10.1007/978-3-540-69934-7](https://doi.org/10.1007/978-3-540-69934-7)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `hexanal_odt_ppb` | `4.5` |
| `aw_dependent_scission_rate_factor` | `1.2` |

**What it supports:**
- Hexanal Odor Detection Threshold (ODT) anchor at 4.5 ppb with water activity-dependent rate.

**What it does not support:**
- Commercial patty texture effects
- Trapping inside soy or pea matrices

**Next Action:** Encode ODT hexanal value for off-note scoring.

---

### Frankel (1983), Prog. Lipid Res. 22:1
- **DOI:** [10.1016/0163-7827(83)90002-4](https://doi.org/10.1016/0163-7827(83)90002-4)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `c182_to_hexanal_beta_scission_yield_pct` | `14.5` |
| `isotopic_recovery_efficiency` | `0.88` |

**What it supports:**
- C18:2 to hexanal beta-scission coefficients and isotopic validation.

**What it does not support:**
- Enzymatic lipoxygenase inhibition
- Aqueous phase pH effects

**Next Action:** Encode lipid hexanal beta-scission yield prior.

---

### Farmer & Patterson (1991), Food Chemistry 40:201
- **DOI:** [10.1016/0308-8146(91)90103-u](https://doi.org/10.1016/0308-8146(91)90103-u)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `alkyl_thiazole_oav_ceiling` | `50.0` |
| `aldehyde_precursor_reactivity_ratio` | `1.2` |

**What it supports:**
- Aliphatic aldehyde + H2S + NH3 thiazole cyclisation kinetics and relative yields.
- Confirms alkyl thiazole OAV > 1 in lipid-crosstalk modeling.

**What it does not support:**
- Protein-matrix retention constants
- pH-dependent oxidation rate limits

**Next Action:** Expose as a lipid-Maillard thiazole cross-coupling mechanistic reference.

---

### Esterbauer, Schaur & Zollner (1991), Free Radic. Biol. Med. 11:81
- **DOI:** [10.1016/0891-5849(91)90192-6](https://doi.org/10.1016/0891-5849(91)90192-6)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `cys_Michael_addition_rate_M_1_s_1` | `1.33` |
| `his_Michael_addition_rate_M_1_s_1` | `0.001` |
| `lys_Schiff_addition_rate_M_1_s_1` | `0.0001` |

**What it supports:**
- Unsaturated aldehyde (4-HNE) scavenging rate hierarchy: Cys >> His > Lys.
- Michael addition rate constants for reactive amino acid sidechains.

**What it does not support:**
- Saturated aldehyde trapping kinetics
- Denaturation-induced accessibility values

**Next Action:** Encode 4-HNE reactive amino acid scavenging rate priors.

---

### Kamal-Eldin, Velasco & Dobarganes (2003), Eur. J. Lipid Sci. Technol. 105:329
- **DOI:** [10.1002/ejlt.200390032](https://doi.org/10.1002/ejlt.200390032)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `c181_to_octanal_beta_scission_yield_pct` | `3.5` |
| `c181_to_nonanal_beta_scission_yield_pct` | `4.2` |

**What it supports:**
- Triolein oxidation and oleic acid (C18:1) beta-scission yield parameters.

**What it does not support:**
- Phospholipid fraction contributions
- LOX inactivation temperatures

**Next Action:** Encode oleic (C18:1) scission volatile split priors.

---

### Tan, Che Man, Selamat & Yusoff (2001), J. Am. Oil Chem. Soc. 78:1133
- **DOI:** [10.1007/s11746-001-0401-1](https://doi.org/10.1007/s11746-001-0401-1)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `computational_prior`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |

**What it supports:**
- Isothermal peroxide formation activation energy parameters derived via DSC.

**What it does not support:**
- Extrusion high-shear physical matrix transformations.

**Next Action:** Encode lipid oxidation Arrhenius prior.

---

### Kavran, Yucel, Bakkalbasi & Gulec (2024), Grasas y Aceites 75:e546
- **DOI:** [10.3989/gya.0320231](https://doi.org/10.3989/gya.0320231)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `computational_prior`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |

**What it supports:**
- Arrhenius hexanal kinetics modulated by antioxidants like ascorbyl palmitate (Ea 1.62 to 89.40 kJ/mol).

**What it does not support:**
- High-moisture extrusion structures without oil phase isolation.

**Next Action:** Encode hexanal accumulation Arrhenius prior.

---

### West, Adams & Zabarnick (2011), Energy & Fuels 25:897
- **DOI:** [10.1021/ef101678s](https://doi.org/10.1021/ef101678s)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `computational_prior`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |
| `Fe_Ea_kcal_mol` | `23.0` |
| `Zn_Ea_kcal_mol` | `26.0` |
| `Cu_Ea_kcal_mol` | `14.0` |
| `Mn_Ea_kcal_mol` | `11.0` |

**What it supports:**
- Apparent Ea for transition-metal catalyzed peroxide decomposition (11.0 to 26.0 kcal/mol).

**What it does not support:**
- Heme-catalyzed lipid breakdown kinetics.

**Next Action:** Encode catalytic peroxide decomposition prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Lipid oxidation and carbonylic crosstalk` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `hexanal`, `2-pentylfuran`, `nonanal`, `carbonyl_competition_factor` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 2

Total papers analyzed: **11** (Benchmark-eligible: **1**, Calibration references: **6**, Rejected: **4**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Trikusuma et al. (2020) | 6/8 | ✅ Eligible | pea_uht_beverage |
| Lin, Fay, Welti & Blank (2001), Lipids 36:749 | 3/8 | ⚠️ Calibration | lipid_rich_model_system |
| Hrncirik & Zeelenberg (2013), J. Am. Oil Chem. Soc. 91:591 | 3/8 | ⚠️ Calibration | coconut_oil_co_matrix |
| Belitz, Grosch & Schieberle (2009), Food Chemistry, 4th rev. ed., Springer (aroma-threshold tables, Ch. 5) | 3/8 | ⚠️ Calibration | free_model_system |
| Frankel (1983), Prog. Lipid Res. 22:1 | 2/8 | ❌ Rejected | free_model_system |
| Farmer & Patterson (1991), Food Chemistry 40:201 | 2/8 | ❌ Rejected | free_model_system |
| Esterbauer, Schaur & Zollner (1991), Free Radic. Biol. Med. 11:81 | 2/8 | ❌ Rejected | free_model_system |
| Kamal-Eldin, Velasco & Dobarganes (2003), Eur. J. Lipid Sci. Technol. 105:329 | 2/8 | ❌ Rejected | free_model_system |
| Tan, Che Man, Selamat & Yusoff (2001), J. Am. Oil Chem. Soc. 78:1133 | 3/8 | ⚠️ Calibration | free_model_system |
| Kavran, Yucel, Bakkalbasi & Gulec (2024), Grasas y Aceites 75:e546 | 3/8 | ⚠️ Calibration | free_model_system |
| West, Adams & Zabarnick (2011), Energy & Fuels 25:897 | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 2

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Trikusuma et al. (2020)` (Score: 6/8)

### CALIBRATION (Score 3-5)
- `Lin, Fay, Welti & Blank (2001), Lipids 36:749` (Score: 3/8)
- `Hrncirik & Zeelenberg (2013), J. Am. Oil Chem. Soc. 91:591` (Score: 3/8)
- `Belitz, Grosch & Schieberle (2009), Food Chemistry, 4th rev. ed., Springer (aroma-threshold tables, Ch. 5)` (Score: 3/8)
- `Tan, Che Man, Selamat & Yusoff (2001), J. Am. Oil Chem. Soc. 78:1133` (Score: 3/8)
- `Kavran, Yucel, Bakkalbasi & Gulec (2024), Grasas y Aceites 75:e546` (Score: 3/8)
- `West, Adams & Zabarnick (2011), Energy & Fuels 25:897` (Score: 3/8)

### REJECTED (Score < 3)
- `Frankel (1983), Prog. Lipid Res. 22:1` (Score: 2/8)
- `Farmer & Patterson (1991), Food Chemistry 40:201` (Score: 2/8)
- `Esterbauer, Schaur & Zollner (1991), Free Radic. Biol. Med. 11:81` (Score: 2/8)
- `Kamal-Eldin, Velasco & Dobarganes (2003), Eur. J. Lipid Sci. Technol. 105:329` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/lipid_oxidation.py` must explicitly account for `lipid_crosstalk_lane` as a modifier when predicting `hexanal`, `2-pentylfuran`, `nonanal`, `carbonyl_competition_factor` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `immediate_expansion_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-10-27 using Maillard Ingestion Pipeline.*
