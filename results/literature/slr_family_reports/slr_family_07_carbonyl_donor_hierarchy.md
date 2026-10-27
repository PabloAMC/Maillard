# Systematic Literature Review — Family 7: Reducing sugar and carbonyl donor hierarchy
**Last updated:** 2026-10-27  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering reducing sugar and carbonyl donor hierarchy.  
**Strategic Posture:** `immediate_expansion_lane`  
**Runtime Concept:** `carbonyl_donor_hierarchy`  
**Scope & Targets:** Covers target compounds/variables: `ribose`, `xylose`, `glucose`, `fructose`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`. Target runtime modules: .  

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

### Zhang et al. (2024), J. Agric. Food Chem. 72:22461
- **DOI:** [10.1021/acs.jafc.4c05736](https://doi.org/10.1021/acs.jafc.4c05736)
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
| `three_dg_mediated_ea_kj_mol` | `115.0` |
| `mgo_mediated_ea_kj_mol` | `85.0` |
| `replicates` | `3` |

**What it supports:**
- Bifurcated, first-order Amadori degradation kinetics: 3-deoxyglucosone-mediated pathway (Ea 115.0 kJ/mol) versus methylglyoxal-mediated pathway (Ea 85.0 kJ/mol).

**What it does not support:**
- High moisture extrusion mechanical effects

**Next Action:** Encode bifurcated Amadori degradation pathway priors.

---

### Martins & van Boekel (2003), Carbohydr. Res. 338:1665
- **DOI:** [10.1016/s0008-6215(03)00174-5](https://doi.org/10.1016/s0008-6215(03)00174-5)
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
- Multiresponse kinetic modeling of DFG Amadori enolization and degradation.

**What it does not support:**
- Solid-state extrusion with low water activity.

**Next Action:** Encode Amadori degradation prior.

---

### Nashalian & Yaylayan (2014), J. Agric. Food Chem. 62:8518
- **DOI:** [10.1021/jf502751n](https://doi.org/10.1021/jf502751n)
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `decarboxylation_k_s1` | `0.41` |

**What it supports:**
- Copper-catalyzed enolization (rate k = 0.41 s-1) and decarboxylation kinetics in model systems.

**What it does not support:**
- Non-transition-metal catalyzed Strecker pathways.

**Next Action:** Encode catalyzed Strecker enolization prior.

---

### Cerny & Davidek (2004), J. Agric. Food Chem. 52:958
- **DOI:** [10.1021/jf035265m](https://doi.org/10.1021/jf035265m)
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `labeling_efficiency_pct` | `90.0` |

**What it supports:**
- Branching ratios from [1-13C]ribose to 3-mercapto-2-pentanone showing 90% labeling

**What it does not support:**
- Enzymatic deamidation shifts

**Next Action:** Encode ribose-to-mercaptoketone branching ratio priors.

---

### Buera, Chirife, Resnik & Wetzler (1987), J. Food Sci. 52:1063
- **DOI:** [10.1111/j.1365-2621.1987.tb14276.x](https://doi.org/10.1111/j.1365-2621.1987.tb14276.x)
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `maillard_ea_kcal_mol_range` | `[25.0, 35.0]` |
| `caramelisation_glucose_ea_kcal_mol_range` | `[33.0, 36.0]` |

**What it supports:**
- Arrhenius activation energies (Ea) for fructose, xylose, and glucose Maillard browning with glycine (25-35 kcal/mol)
- Ea for caramelisation (fructose/xylose 25-30 kcal/mol, glucose 33-36 kcal/mol, maltose/lactose 35-48 kcal/mol)

**What it does not support:**
- Salivary mucin kinetics

**Next Action:** Encode sugar browning and caramelisation activation energy priors.

---

### Fu, Ma, Wang, Sun, Chen, Cheng & Liu (2023), Foods 12:1967
- **DOI:** [10.3390/foods12101967](https://doi.org/10.3390/foods12101967)
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `cml_range_mg_kg` | `[16.0, 48.0]` |
| `cel_range_mg_kg` | `[25.0, 86.0]` |

**What it supports:**
- Quantified safety endpoints in commercial plant-based meat analogues: CML 16-48 mg/kg, CEL 25-86 mg/kg

**What it does not support:**
- NADES extraction parameters

**Next Action:** Expose as a commercial PBMA safety reference payload.

---

### van Boekel (2001)
- **DOI:** [10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9](https://doi.org/10.1002/1521-3803(20010601)45:3<150::AID-FOOD150>3.0.CO;2-9)
- **Matrix Family:** lysine_sugar_model_system
- **Kind:** `directional_prior`
- **Payload Role:** `benchmark_intake`
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
| `ribose_k_obs_at_60c_milli_per_hour` | `98.4` |
| `glucose_k_obs_at_60c_milli_per_hour` | `12.0` |
| `ribose_vs_glucose_relative_reactivity` | `8.2` |
| `glucose_ea_kj_per_mol` | `106.0` |
| `ribose_ea_kj_per_mol` | `114.0` |
| `replicates` | `3` |

**What it supports:**
- Quantitative k_obs and Ea for 9 sugars at pH 6, 50-70°C in lysine model system
- Pentose >> hexose >> disaccharide reactivity hierarchy with temperature-dependence

**What it does not support:**
- Plant-matrix endpoint closure
- A universal Maillard rate under extruder shear or low-moisture conditions

**Next Action:** Wired as the authoritative sugar-reactivity prior for Family 07. Activation energies (glucose 106 kJ/mol, ribose 114 kJ/mol) are now hardcoded in curated_pathways.py.

---

### Martins, Jongen & van Boekel (2000), Trends Food Sci. Technol. 11:364
- **DOI:** [10.1016/s0924-2244(01)00022-x](https://doi.org/10.1016/s0924-2244(01)00022-x)
- **Matrix Family:** lysine_sugar_model_system
- **Kind:** `directional_prior`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_directional_prior_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ❌ | Precursor concentrations not specified |
| ❌ | Temperature, time or pH missing from key_values |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 1/8 → Rejected**

**What it supports:**
- Mechanistic kinetic modelling framework for the Maillard reaction hierarchy

**What it does not support:**
- Standalone plant-matrix benchmark closure

**Next Action:** Keep as co-reference to van_boekel_2001 for mechanistic context in Family 07.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Reducing sugar and carbonyl donor hierarchy` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `ribose`, `xylose`, `glucose`, `fructose` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 7

Total papers analyzed: **8** (Benchmark-eligible: **0**, Calibration references: **3**, Rejected: **5**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Zhang et al. (2024), J. Agric. Food Chem. 72:22461 | 3/8 | ⚠️ Calibration | free_model_system |
| Martins & van Boekel (2003), Carbohydr. Res. 338:1665 | 3/8 | ⚠️ Calibration | free_model_system |
| Nashalian & Yaylayan (2014), J. Agric. Food Chem. 62:8518 | 2/8 | ❌ Rejected | free_model_system |
| Cerny & Davidek (2004), J. Agric. Food Chem. 52:958 | 2/8 | ❌ Rejected | free_model_system |
| Buera, Chirife, Resnik & Wetzler (1987), J. Food Sci. 52:1063 | 2/8 | ❌ Rejected | free_model_system |
| Fu, Ma, Wang, Sun, Chen, Cheng & Liu (2023), Foods 12:1967 | 2/8 | ❌ Rejected | plant_based_meat_analogue |
| van Boekel (2001) | 3/8 | ⚠️ Calibration | lysine_sugar_model_system |
| Martins, Jongen & van Boekel (2000), Trends Food Sci. Technol. 11:364 | 1/8 | ❌ Rejected | lysine_sugar_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 7

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Zhang et al. (2024), J. Agric. Food Chem. 72:22461` (Score: 3/8)
- `Martins & van Boekel (2003), Carbohydr. Res. 338:1665` (Score: 3/8)
- `van Boekel (2001)` (Score: 3/8)

### REJECTED (Score < 3)
- `Nashalian & Yaylayan (2014), J. Agric. Food Chem. 62:8518` (Score: 2/8)
- `Cerny & Davidek (2004), J. Agric. Food Chem. 52:958` (Score: 2/8)
- `Buera, Chirife, Resnik & Wetzler (1987), J. Food Sci. 52:1063` (Score: 2/8)
- `Fu, Ma, Wang, Sun, Chen, Cheng & Liu (2023), Foods 12:1967` (Score: 2/8)
- `Martins, Jongen & van Boekel (2000), Trends Food Sci. Technol. 11:364` (Score: 1/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in  must explicitly account for `carbonyl_donor_hierarchy` as a modifier when predicting `ribose`, `xylose`, `glucose`, `fructose` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `immediate_expansion_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-10-27 using Maillard Ingestion Pipeline.*
