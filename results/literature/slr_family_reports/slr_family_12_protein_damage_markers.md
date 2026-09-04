# Systematic Literature Review — Family 12: Protein Damage Markers
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering protein damage markers.  
**Strategic Posture:** `first_class_runtime_lane`  
**Runtime Concept:** `safety_damage_lane`  
**Scope & Targets:** Covers target compounds/variables: `CML`, `CEL`, `furosine`, `lysinoalanine`. Preferred payload types: `safety_payload`, `guardrail_entry`. Target runtime modules: `src/safety.py`, `src/extrusion_model.py`.  

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

### Wehrmaker et al. (2022)
- **DOI:** [10.1021/acsfoodscitech.2c00242](https://doi.org/10.1021/acsfoodscitech.2c00242)
- **Matrix Family:** spi_ppi_petfood
- **Kind:** `quantitative_benchmark`
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
| `spi_lysine_loss_pct` | `31.4` |
| `ppi_lysine_loss_pct` | `15.4` |
| `sterilization_temp_C` | `125.5` |
| `sterilization_time_min` | `26` |

**What it supports:**
- Definitive lysine loss differential: SPI (31.4%) vs PPI (15.4%)
- Total AA loss quantification after severe shearing+sterilization
- Calibration of lysine_accessibility_factor under extreme conditions

**Next Action:** Use the already-encoded safety reference payload as the primary reactive lysine baseline for SPI vs PPI resilience modeling, and keep the intake row linked to that artifact.

---

### Ma et al. (2024)
- **DOI:** [10.3390/ijms25168668](https://doi.org/10.3390/ijms25168668)
- **Matrix Family:** spi_extrusion_short_residence
- **Kind:** `quantitative_benchmark`
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
| `temp_C` | `130.0` |
| `time_seconds` | `25.0` |
| `acrylamide_start_ug_per_kg` | `22.36` |
| `acrylamide_end_ug_per_kg` | `62.62` |
| `moisture_regime` | `lme` |

**What it supports:**
- Short-residence acrylamide benchmark for SPI extrusion with SME-aware effective temperature
- 22.36 -> 62.62 ug/kg increase in 20-30 s at 130 C as a fast-kinetics closure point

**What it does not support:**
- Raw-flour baseline reconstruction inside the same mechanistic solver
- Generalization to microwave or infrared processing without separate process terms

**Next Action:** Keep encoded as the fastest extrusion acrylamide closure point and route through effective extrusion temperature.

---

### Fu et al. (2023)
- **DOI:** [10.3390/foods12101967](https://doi.org/10.3390/foods12101967)
- **Matrix Family:** commercial_pbma_products
- **Kind:** `quantitative_benchmark`
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
| `cml_range_mg_per_kg` | `[16.46, 47.61]` |
| `cel_range_mg_per_kg` | `[25.21, 86.23]` |
| `representative_cml_mg_per_kg` | `32.0` |
| `representative_cel_mg_per_kg` | `55.0` |

**What it supports:**
- Representative executable CML and CEL proxy benchmark tied to the Foods 2023 commercial product survey
- Keeps the exact Foods 2023 and PMC 2024 AGE range anchors attached to the runtime artifact

**What it does not support:**
- Exact reconstruction of any one commercial product process history
- Substitute for the broader reference bands in safety reporting

**Next Action:** Use as an executable AGE proxy benchmark while retaining Foods 2023 and PMC 2024 as explicit reporting bands.

---

### Ramírez-Jiménez et al. (2000)
- **DOI:** [10.1021/jf9907687](https://doi.org/10.1021/jf9907687)
- **Matrix Family:** mild_legume_extrudate
- **Kind:** `quantitative_benchmark`
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
| `peak_furosine_mg_per_100g_protein` | `8.7` |
| `peak_temp_C` | `140.0` |
| `crossover_temp_C` | `150.0` |

**What it supports:**
- Executable furosine peak anchor near 140 C with explicit crossover semantics above 150 C
- Non-monotonic damage benchmark so furosine decline is not misread as improvement

**What it does not support:**
- Single-marker safety closure for severe extrusion above 160 C
- Exact SPI/PPI source attribution without matrix-specific calibration

**Next Action:** Keep encoded as the reference crossover benchmark and pair it with CEL/CML reporting above the degradation threshold.

---

### Arinzechukwu et al. (2025)
- **DOI:** [10.3390/foods14010096](https://doi.org/10.3390/foods14010096)
- **Matrix Family:** commercial_pba_products
- **Kind:** `safety_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_reference_encoding`

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
| `cml_range_mg_per_kg` | `[16.0, 110.0]` |
| `cel_range_mg_per_kg` | `[25.0, 110.0]` |

**What it supports:**
- A broad commercial AGE validation band for CML and CEL in plant-based analog products
- Family 12 guardrail visibility beyond the narrower Foods 2023 benchmark cluster

**What it does not support:**
- A formulation-specific kinetic reconstruction
- Standalone process attribution without an executable severity proxy

**Next Action:** Keep encoded as the broad Family 12 commercial AGE band and pair it with executable CML/CEL proxy benchmarks for severity interpretation.

---

### Henle (2005)
- **DOI:** [10.1007/s00726-005-0187-z](https://doi.org/10.1007/s00726-005-0187-z)
- **Matrix Family:** soy_isolate
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
| `glycinin_cml_umol_g` | `2.3` |

**What it supports:**
- Quantified CML formation (2.3 umol/g protein) from heated soy glycinin.
- Dicarbonyl consumption modeling.

**What it does not support:**
- Acrylamide kinetics on wheat gluten matrices
- Free volatile headspace partition factors

**Next Action:** Encode soy glycinin dicarbonyl consumption safety reference.

---

### Fu et al. (2023)
- **DOI:** [10.3390/foods12101967](https://doi.org/10.3390/foods12101967)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `acrylamide_range_ug_kg` | `[31.81, 186.7]` |

**What it supports:**
- Acrylamide concentration range of 31.81-186.70 ug/kg under heating above 150 C
- Pictet-Spengler mechanism pathway validation above 150 C
- In-matrix trapping and mitigation dynamics

**What it does not support:**
- Low temperature ambient slurry kinetics

**Next Action:** Expose as a high-temperature acrylamide safety reference payload.

---

### De Vleeschouwer et al. (2006)
- **DOI:** [10.1021/jf051197n](https://doi.org/10.1021/jf051197n)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `kinetic_order_formation` | `first_order` |
| `kinetic_order_elimination` | `first_order` |
| `replicates` | `3` |

**What it supports:**
- Kinetics of acrylamide formation and elimination in aqueous model system under varying pH and temperature.

**What it does not support:**
- Low temperature ambient slurry kinetics

**Next Action:** Expose as acrylamide aqueous kinetics safety payload.

---

### Zhu et al. (2022)
- **DOI:** [10.3389/fnut.2022.940202](https://doi.org/10.3389/fnut.2022.940202)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `acrylamide_formation_ea_kj_mol` | `12.87` |
| `hmf_formation_ea_kj_mol` | `14.85` |
| `replicates` | `3` |

**What it supports:**
- Acrylamide and HMF kinetics in glucose-asparagine-linoleic acid model system showing lipid-Maillard crosstalk acceleration.

**What it does not support:**
- Aqueous buffer without lipids

**Next Action:** Expose as lipid-crosstalk acrylamide safety payload.

---

### Sun et al. (2015)
- **DOI:** [10.1016/j.foodchem.2014.09.129](https://doi.org/10.1016/j.foodchem.2014.09.129)
- **Matrix Family:** beef
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `cml_formation_ea_kj_mol` | `61.01` |
| `cel_formation_ea_kj_mol` | `29.21` |
| `replicates` | `3` |

**What it supports:**
- CML and CEL zero-order formation kinetics in ground beef under pasteurisation conditions.

**What it does not support:**
- High temperature dry extrusion kinetics

**Next Action:** Expose as beef AGEs safety reference payload.

---

### Zhu et al. (2021)
- **DOI:** [10.1002/jsfa.10528](https://doi.org/10.1002/jsfa.10528)
- **Matrix Family:** chicken
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `free_cml_ea_kj_mol` | `74.87` |
| `bound_cml_ea_kj_mol` | `68.21` |
| `free_cel_ea_kj_mol` | `42.45` |
| `bound_cel_ea_kj_mol` | `41.52` |
| `replicates` | `3` |

**What it supports:**
- Free and bound CML/CEL kinetics in braised chicken matrix.

**What it does not support:**
- Dairy model systems

**Next Action:** Expose as chicken AGEs safety reference payload.

---

### Hamzalioglu et al. (2026)
- **DOI:** [10.1021/acs.jafc.5c14296](https://doi.org/10.1021/acs.jafc.5c14296)
- **Matrix Family:** milk
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `lactulosyllysine_formation_ea_md_mol` | `52.1` |
| `replicates` | `3` |

**What it supports:**
- Multi-response kinetic modeling of lactulosyllysine conversion to CML/CEL in UHT milk.

**What it does not support:**
- Solid meat analog matrices

**Next Action:** Expose as milk UHT AGEs safety reference payload.

---

### Krause et al. (2003)
- **DOI:** [10.1007/s00217-003-0685-6](https://doi.org/10.1007/s00217-003-0685-6)
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
| `fructosyllysine_yield_6M_pct` | `32.0` |
| `fructosyllysine_yield_8M_pct` | `46.0` |
| `lactulosyllysine_yield_6M_pct` | `34.0` |
| `lactulosyllysine_yield_8M_pct` | `50.0` |
| `tagatosyllysine_yield_6M_pct` | `42.0` |
| `tagatosyllysine_yield_8M_pct` | `42.0` |
| `replicates` | `3` |

**What it supports:**
- Stoichiometric conversion yields of fructosyllysine, lactulosyllysine, tagatosyllysine, and maltulosyllysine to furosine under 6M vs 8M HCl hydrolysis.

**What it does not support:**
- In vivo biological degradation kinetics

**Next Action:** Encode furosine hydrolysis yields computational priors.

---

### Hidalgo & Pompei (2000)
- **DOI:** [10.1021/jf990120u](https://doi.org/10.1021/jf990120u)
- **Matrix Family:** tomato
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `furosine_formation_ea_kj_mol` | `93.9` |
| `replicates` | `3` |

**What it supports:**
- Furosine zero-order kinetics in tomato products.

**What it does not support:**
- Meat matrices

**Next Action:** Expose as tomato furosine safety reference payload.

---

### Cantre et al. (2007)
- **DOI:** [10.1111/j.1745-4549.2007.00109.x](https://doi.org/10.1111/j.1745-4549.2007.00109.x)
- **Matrix Family:** beef
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `furosine_formation_ea_kj_mol` | `95.7` |
| `replicates` | `3` |

**What it supports:**
- Furosine first-order kinetics during commercial sterilization of corned beef.

**What it does not support:**
- Fluid milk model systems

**Next Action:** Expose as beef furosine safety reference payload.

---

### Yu et al. (2017)
- **DOI:** [10.1016/j.tifs.2020.01.021](https://doi.org/10.1016/j.tifs.2020.01.021)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `cml_formation_ea_kj_mol` | `61.01` |
| `cel_formation_ea_kj_mol` | `29.21` |
| `replicates` | `3` |

**What it supports:**
- Differential thermodynamics for advanced glycation endproducts: CML activation energy of 61.01 kJ/mol versus CEL activation energy of 29.21 kJ/mol.

**What it does not support:**
- Extrusion mechanical shear effects

**Next Action:** Expose as CML/CEL activation energy safety reference payload.

---

### Charissou et al. (2007)
- **DOI:** [10.1021/jf063024j](https://doi.org/10.1021/jf063024j)
- **Matrix Family:** cookie_model
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `temp_range_C` | `[180.0, 220.0]` |
| `replicates` | `3` |

**What it supports:**
- CML acts as a stable, linear (zero-order) accumulator marker under low moisture, while furosine shows a transient, bell-shaped decay curve dependent on temperature (180 to 220 C).

**What it does not support:**
- High moisture extrusion mechanical shear effects

**Next Action:** Expose as cookie CML/furosine safety reference payload.

---

### Fratianni et al. (2016)
- **DOI:** [10.1016/j.foodres.2016.12.009](https://doi.org/10.1016/j.foodres.2016.12.009)
- **Matrix Family:** apricot
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `furosine_formation_ea_kj_mol` | `83.3` |
| `replicates` | `3` |

**What it supports:**
- Zero-order reaction kinetics for furosine formation in low-moisture apricot matrices governed by an activation energy of 83.3 kJ/mol.

**What it does not support:**
- High temperature dry extrusion kinetics

**Next Action:** Expose as apricot furosine safety reference payload.

---

### Ma et al. (2024)
- **DOI:** [10.3390/ijms25168668](https://doi.org/10.3390/ijms25168668)
- **Matrix Family:** plant_based_meat_analogue
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_payload`
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
| `temp_threshold_C` | `160.0` |
| `cel_increase_pct_range` | `[39, 101]` |
| `replicates` | `3` |

**What it supports:**
- High-moisture extrusion selective acceleration of CML via glyoxal transport.
- Elevated barrel temperatures up to 160 C spikes CEL concentration by 39-101% due to rapid thermochemical cycling of methylglyoxal.

**What it does not support:**
- Low temperature ambient slurry kinetics

**Next Action:** Expose as PBMA extrusion safety reference payload.

---

### Knol et al. (2005)
- **DOI:** [10.1021/jf050504m](https://doi.org/10.1021/jf050504m)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
| `formation_Ea_kj_mol` | `52.1` |
| `degradation_Ea_kj_mol` | `72.9` |
| `replicates` | `3` |

**What it supports:**
- Bifurcated acrylamide formation (Ea 52.1 kJ/mol) and thermal degradation (Ea 72.9 kJ/mol) kinetics.

**What it does not support:**
- Low moisture solid-state matrices without water activity correction

**Next Action:** Encode as acrylamide safety reference payload.

---

### De Vleeschouwer et al. (2008)
- **DOI:** [10.1021/bp060389f](https://doi.org/10.1021/bp060389f)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
| `min_elimination_Aw` | `0.82` |
| `replicates` | `3` |

**What it supports:**
- Acrylamide moisture-dependent kinetics showing minimum elimination rate at Aw of 0.82.

**What it does not support:**
- Aqueous solution systems without starch-like structural water restrictions

**Next Action:** Encode as water activity dependent acrylamide safety reference payload.

---

### Ishak et al. (2022)
- **DOI:** [10.1016/j.foodchem.2022.132372](https://doi.org/10.1016/j.foodchem.2022.132372)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
| `phenylalanine_Ea_kj_mol` | `95.36` |
| `proline_Ea_kj_mol` | `114.12` |
| `replicates` | `3` |

**What it supports:**
- PhIP heterocyclic amine formation kinetics from phenylalanine (Ea 95.36 kJ/mol) vs. proline (Ea 114.12 kJ/mol) precursor pools.

**What it does not support:**
- Low-temperature extrusion texturization below 100 C

**Next Action:** Encode as PhIP HCA safety reference payload.

---

### Urugo et al. (2024)
- **DOI:** [10.3390/foods14193295](https://doi.org/10.3390/foods14193295)
- **Matrix Family:** plant_protein_matrix
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
- Quantitation and extraction profile of 10 HCA variants in soy and almond milk beverages.

**What it does not support:**
- Dry extrusion HCA profiles.

**Next Action:** Encode HCA safety reference prior.

---

### Chen et al. (2016)
- **DOI:** [10.1007/s10068-016-0185-5](https://doi.org/10.1007/s10068-016-0185-5)
- **Matrix Family:** high_moisture_melt
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
| `glucose_doping_multiplier` | `4.0` |

**What it supports:**
- Zero-order CML and CEL formation rates during high-moisture sterilization accelerated 4x by glucose doping.

**What it does not support:**
- Low-moisture starch-dominated systems.

**Next Action:** Encode sterilization AGE safety prior.

---

### Pruteanu et al. (2023)
- **DOI:** [10.1016/j.lwt.2024.117316](https://doi.org/10.1016/j.lwt.2024.117316)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `safety_reference_payload`
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
| `glucose_glycine_Ea_kj_mol` | `109.0` |
| `glucose_phenylalanine_Ea_kj_mol` | `145.0` |

**What it supports:**
- Apparent Arrhenius activation energy comparison for glucose/glycine (Ea = 109 kJ/mol) and glucose/phenylalanine (Ea = 145 kJ/mol).

**What it does not support:**
- Complex solid-state protein matrices with native starch.

**Next Action:** Encode browning activation energy safety prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Protein Damage Markers` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `CML`, `CEL`, `furosine`, `lysinoalanine` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 12

Total papers analyzed: **25** (Benchmark-eligible: **0**, Calibration references: **18**, Rejected: **7**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Wehrmaker et al. (2022) | 2/8 | ❌ Rejected | spi_ppi_petfood |
| Ma et al. (2024) | 2/8 | ❌ Rejected | spi_extrusion_short_residence |
| Fu et al. (2023) | 2/8 | ❌ Rejected | commercial_pbma_products |
| Ramírez-Jiménez et al. (2000) | 2/8 | ❌ Rejected | mild_legume_extrudate |
| Arinzechukwu et al. (2025) | 2/8 | ❌ Rejected | commercial_pba_products |
| Henle (2005) | 2/8 | ❌ Rejected | soy_isolate |
| Fu et al. (2023) | 2/8 | ❌ Rejected | free_model_system |
| De Vleeschouwer et al. (2006) | 3/8 | ⚠️ Calibration | free_model_system |
| Zhu et al. (2022) | 3/8 | ⚠️ Calibration | free_model_system |
| Sun et al. (2015) | 3/8 | ⚠️ Calibration | beef |
| Zhu et al. (2021) | 3/8 | ⚠️ Calibration | chicken |
| Hamzalioglu et al. (2026) | 3/8 | ⚠️ Calibration | milk |
| Krause et al. (2003) | 3/8 | ⚠️ Calibration | free_model_system |
| Hidalgo & Pompei (2000) | 3/8 | ⚠️ Calibration | tomato |
| Cantre et al. (2007) | 3/8 | ⚠️ Calibration | beef |
| Yu et al. (2017) | 3/8 | ⚠️ Calibration | free_model_system |
| Charissou et al. (2007) | 3/8 | ⚠️ Calibration | cookie_model |
| Fratianni et al. (2016) | 3/8 | ⚠️ Calibration | apricot |
| Ma et al. (2024) | 3/8 | ⚠️ Calibration | plant_based_meat_analogue |
| Knol et al. (2005) | 3/8 | ⚠️ Calibration | free_model_system |
| De Vleeschouwer et al. (2008) | 3/8 | ⚠️ Calibration | free_model_system |
| Ishak et al. (2022) | 3/8 | ⚠️ Calibration | free_model_system |
| Urugo et al. (2024) | 3/8 | ⚠️ Calibration | plant_protein_matrix |
| Chen et al. (2016) | 3/8 | ⚠️ Calibration | high_moisture_melt |
| Pruteanu et al. (2023) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 12

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `De Vleeschouwer et al. (2006)` (Score: 3/8)
- `Zhu et al. (2022)` (Score: 3/8)
- `Sun et al. (2015)` (Score: 3/8)
- `Zhu et al. (2021)` (Score: 3/8)
- `Hamzalioglu et al. (2026)` (Score: 3/8)
- `Krause et al. (2003)` (Score: 3/8)
- `Hidalgo & Pompei (2000)` (Score: 3/8)
- `Cantre et al. (2007)` (Score: 3/8)
- `Yu et al. (2017)` (Score: 3/8)
- `Charissou et al. (2007)` (Score: 3/8)
- `Fratianni et al. (2016)` (Score: 3/8)
- `Ma et al. (2024)` (Score: 3/8)
- `Knol et al. (2005)` (Score: 3/8)
- `De Vleeschouwer et al. (2008)` (Score: 3/8)
- `Ishak et al. (2022)` (Score: 3/8)
- `Urugo et al. (2024)` (Score: 3/8)
- `Chen et al. (2016)` (Score: 3/8)
- `Pruteanu et al. (2023)` (Score: 3/8)

### REJECTED (Score < 3)
- `Wehrmaker et al. (2022)` (Score: 2/8)
- `Ma et al. (2024)` (Score: 2/8)
- `Fu et al. (2023)` (Score: 2/8)
- `Ramírez-Jiménez et al. (2000)` (Score: 2/8)
- `Arinzechukwu et al. (2025)` (Score: 2/8)
- `Henle (2005)` (Score: 2/8)
- `Fu et al. (2023)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/safety.py`, `src/extrusion_model.py` must explicitly account for `safety_damage_lane` as a modifier when predicting `CML`, `CEL`, `furosine`, `lysinoalanine` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `first_class_runtime_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
