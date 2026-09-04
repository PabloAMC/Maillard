# Systematic Literature Review — Family 7: Reducing sugar and carbonyl donor hierarchy
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering reducing sugar and carbonyl donor hierarchy.  
**Strategic Posture:** `immediate_expansion_lane`  
**Runtime Concept:** `carbonyl_donor_hierarchy`  
**Scope & Targets:** Covers target compounds/variables: `ribose`, `xylose`, `glucose`, `fructose`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`. Target runtime modules: `src/precursor_resolver.py`, `src/literature_runtime.py`.  

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

### Cerny & Davídek (2003)
- **DOI:** [10.1021/jf026123f](https://doi.org/10.1021/jf026123f)
- **Matrix Family:** aqueous_model_system
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `g6p_hdmf_multiplier_vs_glucose` | `4.4` |
| `analytical_method` | `SIDA_[2H3]-HDMF` |
| `replicates` | `3` |

**What it supports:**
- Bounded G6P-versus-glucose uplift for HDMF-rich phosphorylated-sugar routing
- Directional support for phosphorylated donor prioritization when furanone output matters

**What it does not support:**
- Direct plant-matrix closure
- Universal HDMF endpoint concentrations

**Next Action:** Encode as a bounded phosphorylated-sugar prior so G6P can outrank glucose when HDMF-support signals matter.

---

### Maillard & van Boekel (1992)
- **DOI:** [10.1016/0308-8146(92)90216-Q](https://doi.org/10.1016/0308-8146(92)90216-Q)
- **Matrix Family:** aqueous_sugar_lysine_model_system
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
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `ribose_vs_glucose_multiplier` | `8.2` |
| `xylose_vs_glucose_multiplier` | `5.9` |
| `fructose_vs_glucose_multiplier` | `0.65` |
| `maltose_vs_glucose_multiplier` | `0.1` |
| `glucose_ea_kj_per_mol` | `106.0` |
| `ribose_ea_kj_per_mol` | `114.0` |

**What it supports:**
- Quantitative sugar-reactivity hierarchy across ribose, pentoses, hexoses, and disaccharides
- A temperature-aware donor-priority anchor so ribose and other pentoses are not treated as fixed scalar uplifts relative to glucose

**What it does not support:**
- Direct plant-matrix closure for PPI or SPI systems
- A full downstream sulfur-yield benchmark without matching amino-acid and matrix context

**Next Action:** Keep as the quantitative donor-hierarchy prior so Family 07 can reason over pentose-versus-hexose reactivity with explicit k_obs and Ea support rather than only R5P-specific evidence.

---

### Blank et al. (1997)
- **DOI:** [10.1021/jf960777q](https://doi.org/10.1021/jf960777q)
- **Matrix Family:** aqueous_model_system
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `hdmf_molar_yield_fraction_lower_bound` | `0.4` |
| `hdmf_odt_ug_per_l` | `0.6` |
| `hdmf_oav_lower_bound` | `80.0` |
| `temperature_c` | `145.0` |
| `replicates` | `3` |

**What it supports:**
- A bounded rhamnose-plus-proline HDMF routing anchor with unusually high molar yield
- Scientist-facing context for why minor rhamnose supplementation can matter when proline is available

**What it does not support:**
- A universal rhamnose donor hierarchy across all amino-acid partners
- Direct plant-matrix HDMF benchmark closure

**Next Action:** Encode as a bounded rhamnose-plus-proline HDMF prior so Family 07 can surface this donor-specific shortcut without rewriting the broader sugar hierarchy.

---

### Liu et al. (2022)
- **DOI:** [10.1016/j.carbpol.2022.120468](https://doi.org/10.1016/j.carbpol.2022.120468)
- **Matrix Family:** pea_protein_isolate
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
| `pea_maltose_liberation_umol_g` | `18.4` |
| `waxy_maize_maltose_liberation_umol_g` | `31.8` |
| `ti_starch_reduction_pct` | `66.0` |

**What it supports:**
- Starch reducing-sugar liberation during HME (maltose generation)
- Compares pea starch (18.4 umol/g) versus waxy maize starch (31.8 umol/g) yield parameters

**What it does not support:**
- Free ribose or glucose addition rates
- Volatile recovery yields under vacuum packaging

**Next Action:** Encode starch reducing-sugar liberation priors during HME.

---

### Huang et al. (2024)
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

### Martins & van Boekel (2003)
- **DOI:** [10.1021/jf025830v](https://doi.org/10.1021/jf025830v)
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

### Nashalian & Yaylayan (2014)
- **DOI:** [10.1021/acsomega.7b00321](https://doi.org/10.1021/acsomega.7b00321)
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

### Wang et al. (2024)
- **DOI:** [10.3390/foods13213453](https://doi.org/10.3390/foods13213453)
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
| `amino_depletion_umol_ml` | `0.69` |
| `tvc_yield_ug_g` | `28.3` |

**What it supports:**
- Glucosamine synergy and amino group depletion kinetics
- Glucosamine MRPs exhibit 0.69 umol/mL amino group depletion and 28.3 ug/g TVC yield with low bitterness

**What it does not support:**
- GSH sulfur release rates

**Next Action:** Encode glucosamine synergy priors.

---

### Blank & Mottram (Mercaptoketone labelling 2002)
- **DOI:** [10.1021/jf020582d](https://doi.org/10.1021/jf020582d)
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

### Buera et al. (1987)
- **DOI:** [10.1111/j.1365-2621.1987.tb14251.x](https://doi.org/10.1111/j.1365-2621.1987.tb14251.x)
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

### Poulsen et al. (2023)
- **DOI:** [10.1016/j.foodchem.2023.136200](https://doi.org/10.1016/j.foodchem.2023.136200)
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

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Reducing sugar and carbonyl donor hierarchy` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `ribose`, `xylose`, `glucose`, `fructose` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 7

Total papers analyzed: **11** (Benchmark-eligible: **0**, Calibration references: **4**, Rejected: **7**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Cerny & Davídek (2003) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Maillard & van Boekel (1992) | 2/8 | ❌ Rejected | aqueous_sugar_lysine_model_system |
| Blank et al. (1997) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Liu et al. (2022) | 2/8 | ❌ Rejected | pea_protein_isolate |
| Huang et al. (2024) | 3/8 | ⚠️ Calibration | free_model_system |
| Martins & van Boekel (2003) | 3/8 | ⚠️ Calibration | free_model_system |
| Nashalian & Yaylayan (2014) | 2/8 | ❌ Rejected | free_model_system |
| Wang et al. (2024) | 2/8 | ❌ Rejected | free_model_system |
| Blank & Mottram (Mercaptoketone labelling 2002) | 2/8 | ❌ Rejected | free_model_system |
| Buera et al. (1987) | 2/8 | ❌ Rejected | free_model_system |
| Poulsen et al. (2023) | 2/8 | ❌ Rejected | plant_based_meat_analogue |

---

## Consolidated entries for benchmark_schema.json — Family 7

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Cerny & Davídek (2003)` (Score: 3/8)
- `Blank et al. (1997)` (Score: 3/8)
- `Huang et al. (2024)` (Score: 3/8)
- `Martins & van Boekel (2003)` (Score: 3/8)

### REJECTED (Score < 3)
- `Maillard & van Boekel (1992)` (Score: 2/8)
- `Liu et al. (2022)` (Score: 2/8)
- `Nashalian & Yaylayan (2014)` (Score: 2/8)
- `Wang et al. (2024)` (Score: 2/8)
- `Blank & Mottram (Mercaptoketone labelling 2002)` (Score: 2/8)
- `Buera et al. (1987)` (Score: 2/8)
- `Poulsen et al. (2023)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/precursor_resolver.py`, `src/literature_runtime.py` must explicitly account for `carbonyl_donor_hierarchy` as a modifier when predicting `ribose`, `xylose`, `glucose`, `fructose` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `immediate_expansion_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
