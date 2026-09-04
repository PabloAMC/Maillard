# Systematic Literature Review — Family 14: Ascorbic Acid Maillard
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering ascorbic acid maillard.  
**Strategic Posture:** `bounded_upstream_source`  
**Runtime Concept:** `stealth_dicarbonyl_source`  
**Scope & Targets:** Covers target compounds/variables: `ascorbic_acid_loss`, `3-deoxythreosone`, `pentosidine_load`. Preferred payload types: `bounded_source_term_payload`, `computational_prior`. Target runtime modules: `src/precursor_resolver.py`, `src/flavor.py`, `src/safety.py`.  

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

### Yu et al. (2020)
- **DOI:** [10.1016/j.foodchem.2020.126458](https://doi.org/10.1016/j.foodchem.2020.126458)
- **Matrix Family:** aqueous_model_system
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
| `pyrraline_ea_kj_mol` | `53.45` |
| `pyrraline_ea_uncertainty_kj` | `4.02` |
| `furosine_ea_kj_mol` | `81.7` |
| `furosine_ea_uncertainty_kj` | `14.01` |

**What it supports:**
- Arrhenius anchors for 3-DG routing into pyrraline and furosine
- Bounded AGE-lane calibration where pyrraline outruns furosine under the same 3-DG burden

**What it does not support:**
- A complete AGE kinetic mechanism
- Matrix-specific endpoint closure for one commercial PBMA

**Next Action:** Keep the two 3-DG damage routes encoded as bounded Arrhenius anchors rather than opening a new full AGE mechanism.

---

### Feng et al. (2023)
- **DOI:** [10.3389/fnut.2022.1022254](https://doi.org/10.3389/fnut.2022.1022254)
- **Matrix Family:** ascorbic_hcw_model
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
| `ea_kj_per_mol_by_ph` | `{'5.0': 15.77, '7.0': 31.7, '9.5': 47.53}` |
| `kinetic_order` | `pseudo_first_order` |
| `temperature_window_c` | `[150.0, 190.0]` |

**What it supports:**
- Primary Family 14 Arrhenius anchor for pH-sensitive ascorbic-acid degradation under extrusion-like thermal conditions
- Explicit HCW calibration for AA-driven dicarbonyl release instead of treating AA as a static precursor burden

**What it does not support:**
- Direct AGE endpoint closure in a real SPI or PPI extrudate
- A full oxygen-transition model without additional barrel-state evidence

**Next Action:** Keep encoded as a bounded Arrhenius prior until a direct SPI/PPI ascorbate extrusion benchmark exists.

---

### Grandhee & Monnier (1991)
- **DOI:** [10.1074/jbc.266.18.11644](https://doi.org/10.1074/jbc.266.18.11644)
- **Matrix Family:** bsa_ascorbic_crosslink_model
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
| `pentosidine_glucose_pmol_per_mg_range` | `[13.2, 17.0]` |
| `pentosidine_ascorbic_acid_pmol_per_mg_range` | `[13.2, 17.0]` |
| `molar_equivalence_to_glucose` | `True` |
| `temperature_c` | `37.0` |

**What it supports:**
- A bounded Family 14 pentosidine ceiling showing that ascorbic acid can match glucose on a molar basis as a cross-link precursor
- An explicit absolute-yield anchor for AA-driven damage reporting rather than a purely narrative AGE warning

**What it does not support:**
- A direct SPI or PPI extrusion endpoint benchmark for pentosidine
- A full oxygen-dependent AA transition network in one payload

**Next Action:** Keep encoded as the Family 14 pentosidine-equivalence anchor so antioxidant AA additions remain visible as a bounded AGE burden source in reporting and safety lanes.

---

### Yu et al. (2018)
- **DOI:** [10.1590/1678-457X.08717](https://doi.org/10.1590/1678-457X.08717)
- **Matrix Family:** ascorbic_amino_crosslink_model
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
| `ea_kj_per_mol_by_amino_acid` | `{'lysine': 54.94, 'arginine': 50.08, 'histidine': 35.31}` |
| `observed_rate_order` | `['lysine', 'arginine', 'histidine']` |

**What it supports:**
- Amino-acid hierarchy for AA-derived cross-link routing with explicit Lys, Arg, and His Arrhenius anchors
- Safety-facing Family 14 logic that keeps lysine as the dominant practical sink even when histidine has the lowest Ea

**What it does not support:**
- Direct pentosidine or pyrraline closure in plant-protein extrudates
- A universal cross-link law outside the reported AA model systems

**Next Action:** Keep encoded as a bounded cross-link hierarchy prior until same-run AA damage-marker panels are benchmarked in plant matrices.

---

### Bi et al. (2020)
- **DOI:** [10.1021/acs.jafc.9b07711](https://doi.org/10.1021/acs.jafc.9b07711)
- **Matrix Family:** xylose_alanine_model_system
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
| `baseline_arp_formation_ea_kj_per_mol` | `77.8` |
| `egcg_arp_formation_ea_kj_per_mol` | `62.8` |
| `ea_reduction_kj_per_mol` | `15.0` |
| `max_arp_yield_improvement_percent` | `95.0` |

**What it supports:**
- Activation energy reduction for ARP formation from 77.8 to 62.8 kJ/mol in the presence of EGCG
- EGCG mediated trapping of reactive deoxyosones, improving yield of N-(1-deoxy-D-xylulos-1-yl)alanine up to 95%

**What it does not support:**
- Direct plant-protein extrusion kinetics
- Generalized polyphenol-sugar reactions without amino acid partners

**Next Action:** Keep as the Family 14 EGCG deoxyosone trapping and ARP formation catalysis prior.

---

### Nemet & Monnier (2011)
- **DOI:** [10.1074/jbc.M111.245100](https://doi.org/10.1074/jbc.M111.245100)
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
| `three_deoxythreosone_pmol_mg` | `9.1` |
| `xylosone_pmol_mg` | `0.5` |
| `dha_molar_yield_pct` | `28.0` |
| `dkg_23_molar_yield_pct` | `4.7` |
| `xylosone_molar_yield_pct` | `5.8` |

**What it supports:**
- Absolute dicarbonyl yields from ascorbic acid: 3-deoxythreosone 9.1 pmol/mg, xylosone 0.5 pmol/mg
- Molar yields: DHA 28%, 2,3-DKG 4.7%, xylosone 5.8%

**What it does not support:**
- Acrylamide accumulation rates at 160 C

**Next Action:** Encode ascorbic dicarbonyl release priors.

---

### Smuda & Glomb (2013)
- **DOI:** [10.1002/anie.201300399](https://doi.org/10.1002/anie.201300399)
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
| `oxidative_alpha_fragmentation_pct` | `31.0` |
| `beta_cleavage_pct` | `32.0` |
| `decarboxylation_pct` | `12.0` |

**What it supports:**
- Degradation mass balance of vitamin C (75% accounted for) producing CML (carboxymethyl-lysine) via threose and xylosone intermediates.

**What it does not support:**
- Non-enzymatic browning from stable aldohexose sugars.

**Next Action:** Encode vitamin C degradation mass balance prior.

---

### Serpen & Gökmen (2007)
- **DOI:** [10.1016/j.foodchem.2006.11.073](https://doi.org/10.1016/j.foodchem.2006.11.073)
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
| `oxidation_acceleration_copper_fold` | `88.0` |
| `oxidation_acceleration_iron_fold` | `14.0` |

**What it supports:**
- Reversible first-order consecutive degradation kinetics of ascorbic acid, showing 88-fold copper and 14-fold iron catalytic acceleration.

**What it does not support:**
- Anaerobic degradation kinetics without catalytic metals.

**Next Action:** Encode metal-catalyzed ascorbic acid oxidation prior.

---

### Yang et al. (2021)
- **DOI:** [10.47836/ifrj.28.3.16](https://doi.org/10.47836/ifrj.28.3.16)
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
| `Ea_excess_ascorbate_kj_mol` | `60.76` |
| `Ea_excess_glycine_kj_mol` | `70.16` |

**What it supports:**
- Pseudo-first-order non-enzymatic browning kinetics between AA and glycine (Ea = 60.76 kJ/mol at 4:1 AA:Gly, 70.16 kJ/mol at 1:4 AA:Gly).

**What it does not support:**
- Polypeptide matrix denaturation shifts.

**Next Action:** Encode AA-glycine browning activation energy prior.

---

### Yu et al. (2018)
- **DOI:** [10.1590/1678-457x.08717](https://doi.org/10.1590/1678-457x.08717)
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
| `Ea_lysine_kj_mol` | `54.94` |
| `Ea_arginine_kj_mol` | `50.08` |
| `Ea_histidine_kj_mol` | `35.31` |

**What it supports:**
- Zero-order macroscopic browning kinetics of AA with basic amino acids, with lysine producing the fastest rate despite a higher Ea of 54.94 kJ/mol.

**What it does not support:**
- Non-basic aliphatic amino acid browning rates.

**Next Action:** Encode basic amino acid AA-browning prior.

---

### Manso et al. (2001)
- **DOI:** [10.1046/j.1365-2621.2001.t01-1-00460.x](https://doi.org/10.1046/j.1365-2621.2001.t01-1-00460.x)
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
| `Ea_ascorbic_degradation_kj_mol` | `55.6` |

**What it supports:**
- Weibullian degradation kinetics of ascorbic acid in orange juice (Ea = 55.6 kJ/mol) under aerobic storage conditions.

**What it does not support:**
- High-temperature short-time sterilization limits.

**Next Action:** Encode orange juice ascorbic degradation Weibull prior.

---

### Takase et al. (2025)
- **DOI:** [10.1021/acsfoodscitech.4c00956](https://doi.org/10.1021/acsfoodscitech.4c00956)
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
- Anaerobic degradation pathway of dehydroascorbic acid to xylosone and threosone in lemon juice under oxygen-excluded conditions.

**What it does not support:**
- Transition metal-catalyzed aerobic oxidation rates.

**Next Action:** Encode anaerobic DHAA degradation prior.

---

### Hendrickx et al. (1998)
- **DOI:** [10.1021/jf9708251](https://doi.org/10.1021/jf9708251)
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
| `Ea_ascorbic_degradation_kj_mol` | `54.8` |
| `D_value_121C_min` | `246.0` |
| `z_value_C` | `27.15` |

**What it supports:**
- Isobaric-isothermal degradation kinetics of L-ascorbic acid, yielding a D-value at 121C of 246 min, z-value of 27.15C, and Ea of 54.8 kJ/mol.

**What it does not support:**
- Alkaline solution phase isomerization pathways.

**Next Action:** Encode squeezed tomato isobaric degradation prior.

---

### Jian et al. (2012)
- **DOI:** [10.1021/jf3032342](https://doi.org/10.1021/jf3032342)
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
| `Ea_ethanolic_min_kj_mol` | `43.3` |
| `Ea_ethanolic_max_kj_mol` | `96.6` |

**What it supports:**
- Ascorbic acid degradation kinetics in ethanolic solutions, where lower water activity accelerates L-xylosone dehydration (Ea = 43.3 to 96.6 kJ/mol).

**What it does not support:**
- High-moisture extrusion structures.

**Next Action:** Encode ethanolic ascorbic degradation prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Ascorbic Acid Maillard` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `ascorbic_acid_loss`, `3-deoxythreosone`, `pentosidine_load` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 14

Total papers analyzed: **14** (Benchmark-eligible: **0**, Calibration references: **10**, Rejected: **4**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Yu et al. (2020) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Feng et al. (2023) | 2/8 | ❌ Rejected | ascorbic_hcw_model |
| Grandhee & Monnier (1991) | 2/8 | ❌ Rejected | bsa_ascorbic_crosslink_model |
| Yu et al. (2018) | 2/8 | ❌ Rejected | ascorbic_amino_crosslink_model |
| Bi et al. (2020) | 3/8 | ⚠️ Calibration | xylose_alanine_model_system |
| Nemet & Monnier (2011) | 2/8 | ❌ Rejected | free_model_system |
| Smuda & Glomb (2013) | 3/8 | ⚠️ Calibration | free_model_system |
| Serpen & Gökmen (2007) | 3/8 | ⚠️ Calibration | free_model_system |
| Yang et al. (2021) | 3/8 | ⚠️ Calibration | free_model_system |
| Yu et al. (2018) | 3/8 | ⚠️ Calibration | free_model_system |
| Manso et al. (2001) | 3/8 | ⚠️ Calibration | free_model_system |
| Takase et al. (2025) | 3/8 | ⚠️ Calibration | free_model_system |
| Hendrickx et al. (1998) | 3/8 | ⚠️ Calibration | free_model_system |
| Jian et al. (2012) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 14

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Yu et al. (2020)` (Score: 3/8)
- `Bi et al. (2020)` (Score: 3/8)
- `Smuda & Glomb (2013)` (Score: 3/8)
- `Serpen & Gökmen (2007)` (Score: 3/8)
- `Yang et al. (2021)` (Score: 3/8)
- `Yu et al. (2018)` (Score: 3/8)
- `Manso et al. (2001)` (Score: 3/8)
- `Takase et al. (2025)` (Score: 3/8)
- `Hendrickx et al. (1998)` (Score: 3/8)
- `Jian et al. (2012)` (Score: 3/8)

### REJECTED (Score < 3)
- `Feng et al. (2023)` (Score: 2/8)
- `Grandhee & Monnier (1991)` (Score: 2/8)
- `Yu et al. (2018)` (Score: 2/8)
- `Nemet & Monnier (2011)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/precursor_resolver.py`, `src/flavor.py`, `src/safety.py` must explicitly account for `stealth_dicarbonyl_source` as a modifier when predicting `ascorbic_acid_loss`, `3-deoxythreosone`, `pentosidine_load` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `bounded_upstream_source` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
