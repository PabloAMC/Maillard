# Systematic Literature Review — Family 3: Thiamine degradation and sulfur support
**Last updated:** 2026-10-27  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering thiamine degradation and sulfur support.  
**Strategic Posture:** `high_value_support_lane`  
**Runtime Concept:** `thiamine_fragmentation_support`  
**Scope & Targets:** Covers target compounds/variables: `thiamine_availability`, `2-methyl-3-furanthiol`. Preferred payload types: `computational_prior`, `flavor_reference_payload`. Target runtime modules: .  

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

### Voelker et al. (2021)
- **DOI:** [10.1186/s13065-021-00773-y](https://doi.org/10.1186/s13065-021-00773-y)
- **Matrix Family:** aqueous_thiamine_model
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
| ✅ | Replicates reported |
| ✅ | LOD/LOQ stated |
| ❌ | No sensory or odor threshold tags present |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `temp_c_points` | `[40.0, 55.0, 70.0]` |
| `ph_points` | `[3.0, 4.0, 5.0, 6.0, 7.0]` |
| `tmn_low_conc_k_obs_per_h_ph3_70c` | `0.0011` |
| `tmn_low_conc_k_obs_per_h_ph6_70c` | `0.0435` |
| `tmn_low_conc_ea_kcal_per_mol_ph3` | `28.4` |
| `tmn_low_conc_ea_kcal_per_mol_ph6_midpoint` | `18.5` |
| `neutral_ph_high_concentration_multiplier` | `9.3` |
| `storage_retention_pct_25c_ph3_day392` | `91.0` |
| `replicates` | `3` |

**What it supports:**
- Primary Family 03 Arrhenius anchor for pH-sensitive thiamine disappearance
- Bounded concentration dependence at neutral pH instead of a single global disappearance rate

**What it does not support:**
- Direct MFT-forming step kinetics separate from overall thiamine disappearance
- Real plant-matrix metal effects without an explicit cofactor correction

**Next Action:** Use as the primary Family 03 Arrhenius calibration anchor while keeping Cerny as the default flavor-yield prior for the runtime lane.

---

### Huang (2022)
- **Matrix Family:** metal_catalyzed_thiamine_model
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_calibration_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ✅ | LOD/LOQ stated |
| ❌ | No sensory or odor threshold tags present |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `ph` | `6.0` |
| `temp_c_points` | `[25.0, 55.0]` |
| `fe3_acceleration_vs_control` | `19.0` |
| `cu_plus_acceleration_vs_control` | `12.0` |
| `fe2_acceleration_vs_control` | `8.0` |
| `zn2_acceleration_vs_control` | `1.5` |
| `fe3_loss_pct_55c_day7` | `95.98` |
| `replicates` | `3` |

**What it supports:**
- Bounded Fe3+ and Cu+ acceleration factors for thiamine loss in metal-bearing matrices
- A direct matrix-correction bridge from aqueous thiamine kinetics to iron-rich plant systems

**What it does not support:**
- Compound-specific sulfur volatile closure under the same run
- Extrusion shear effects beyond metal-assisted degradation

**Next Action:** Use as the compact Family 03 matrix-correction anchor for iron- and copper-rich systems without changing the default thiamine flavor prior.

---

### Arabshahi & Lund (1988), J. Food Sci. 53:199
- **DOI:** [10.1111/j.1365-2621.1988.tb10208.x](https://doi.org/10.1111/j.1365-2621.1988.tb10208.x)
- **Matrix Family:** starch_based_low_moisture_model
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
| ✅ | Replicates reported |
| ✅ | LOD/LOQ stated |
| ❌ | No sensory or odor threshold tags present |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `water_activity_points` | `[0.3, 0.65, 0.9]` |
| `temperature_c_points` | `[110.0, 130.0, 150.0]` |
| `starch_ea_kcal_per_mol_by_aw` | `{'0.3': 26.8, '0.65': 23.4, '0.9': 20.1}` |
| `cellulose_ea_kcal_per_mol_aw_0_65` | `24.1` |
| `replicates` | `3` |

**What it supports:**
- A bounded Family 03 aw-dependent Arrhenius anchor for thiamine loss in starch-like matrices
- Runtime-visible attenuation of thiamine support under wetter extrusion-like conditions instead of treating water activity as a generic scalar

**What it does not support:**
- Direct MFT-forming step kinetics separate from overall thiamine disappearance
- A universal plant-matrix correction outside starch-like or starch-containing systems

**Next Action:** Use as the bounded Family 03 aw-dependent support anchor so wetter extrusion-like thiamine systems no longer inherit the same effective Arrhenius posture as drier starch-rich systems.

---

### Cerny & Guntz-Dubini (2008)
- **DOI:** [10.1021/jf801762c](https://doi.org/10.1021/jf801762c)
- **Matrix Family:** free_model_system
- **Kind:** `quantitative_benchmark`
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
| `temp_C` | `145.0` |
| `time_min` | `20.0` |
| `ph_points` | `[4.0, 5.0, 6.0, 7.0, 8.0]` |
| `thiamine_alone_mft_ug_per_g` | `{'4.0': 0.02, '5.0': 0.14, '6.0': 0.58, '7.0': 0.42, '8.0': 0.19}` |
| `mixed_system_mft_ug_per_g` | `{'5.5': 3.11, '6.0': 2.47}` |
| `mixed_system_synergy_factor_vs_thiamine_alone` | `4.3` |
| `replicates` | `3` |
| `volatile_output_mode` | `absolute_ug_per_g` |

**What it supports:**
- pH-resolved quantitative MFT surface for thiamine alone at 145 C / 20 min
- Mixed thiamine plus cysteine plus xylose uplift near the pH 5.5 to 6.0 optimum window
- Benchmark-facing anchor for Family 03 sulfur support instead of generic availability-only heuristics

**What it does not support:**
- SPI or PPI matrix retention effects
- Extrusion survival before cooking
- A complete consumer-cook kinetic reconstruction outside the reported pH and temperature window

**Next Action:** Expose as the pH-calibrated Family 03 benchmark anchor for thiamine-driven MFT support in the runtime lane.

---

### Bolton, Reineccius, Liardon & Huynh-Ba (1994), ACS Symp. Ser. 543:270-278 -- REPLACES the nonexistent "Hofmann et al. (1996)"
- **DOI:** [10.1021/bk-1994-0543.ch022](https://doi.org/10.1021/bk-1994-0543.ch022)
- **Matrix Family:** free_model_system
- **Kind:** `quantitative_benchmark`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output mode specified |
| ✅ | Absolute yields reported |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 6/8 → Benchmark-eligible**

**Key Values:**
| Parameter | Value |
|---|---|
| `temp_C` | `120.0` |
| `time_min` | `60.0` |
| `ph` | `None` |
| `ph_note` | `NOT REPORTED in the source. The rebuilt benchmark JSON pins 5.0 only because the loader requires a value; it is not an anchor.` |
| `thiamine_mM` | `None` |
| `cysteine_mM` | `None` |
| `sugar` | `D-glucose` |
| `molarity_note` | `NOT REPORTED in the indexed abstract. The rebuilt benchmark JSON uses assumed process-flavor loadings (thiamine 1.0 mM, cysteine 10.0 mM, glucose 10.0 mM), marked assumed_not_from_source.` |
| `mft_ppb_range` | `[11.7, 14.7]` |
| `mft_ppb_target` | `13.0` |
| `mft_sulfur_from_cysteine_pct` | `[7.5, 8.0]` |
| `model_variants` | `4` |
| `volatile_output_mode` | `absolute_ppb` |

**What it supports:**
- A quantitative MFT anchor for the thiamine + cysteine + hexose system at 120 C / 60 min (GC-MS-SIM, 389-489 ng per 33.3 g = 11.7-14.7 ppb)
- The qualitative claim that thiamine is the dominant sulfur source for MFT in this system (only 7.5-8.0 percent of MFT sulfur came from 34S-cysteine; no thiamine gave no MFT)
- An order-of-magnitude Family 03 scale check independent of the high-concentration Cerny 2008 pentose optimum

**What it does not support:**
- Any pH response -- the source does not report pH
- A thiamine-alone vs cysteine+sugar synergy factor -- the previously recorded 2.64 had no source and was withdrawn 2026-08-26
- Extrusion survival or encapsulation protection
- Plant-protein matrix trapping or release

**Next Action:** Use as the Family 03 hexose-system MFT scale check. NOTE (2026-08-26): the model currently under-predicts this anchor by ~163x (0.080 ppb vs 13 ppb) and the gap is robust (35-745x) across the assumed-value sweep -- this is now the panel's honest thiamine failure signal, not a passing row.

---

### Tang, Jiang, Yuan & Ho (2013), J. Sulfur Chem. 34:326 (CrossRef issue year 2012)
- **DOI:** [10.1080/17415993.2012.708933](https://doi.org/10.1080/17415993.2012.708933)
- **Matrix Family:** free_model_system
- **Kind:** `calibration_reference`
- **Payload Role:** `flavor_reference_payload`
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
| `mft_odt_ppb` | `0.05` |

**What it supports:**
- Compiled MFT table across meat systems and odor detection threshold (ODT) parameters.

**What it does not support:**
- Radical scavenging rate constants
- Extruded fiber size distribution

**Next Action:** Expose as a thiamine-MFT flavor reference payload.

---

### Mauri, Alzamora, Chirife & Tomio (1989), Int. J. Food Sci. Technol. 24:1
- **DOI:** [10.1111/j.1365-2621.1989.tb00613.x](https://doi.org/10.1111/j.1365-2621.1989.tb00613.x)
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
| `dilute_solution_ea_kj_mol` | `78.5` |
| `low_aw_ea_kj_mol` | `105.0` |

**What it supports:**
- Compiled activation energies (Ea) for thiamine degradation over wide T and pH range
- Low-Aw Ea values (105.0 kJ/mol) compared to dilute solutions (78.5 kJ/mol)

**What it does not support:**
- Cysteine-accessible nucleophilicity indices

**Next Action:** Encode thiamine activation energy priors.

---

### Voelker et al. (2018)
- **DOI:** [10.1016/j.foodres.2018.06.056](https://doi.org/10.1016/j.foodres.2018.06.056)
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
| `TMN_Ea_range_kj_mol` | `[88.0, 105.0]` |
| `TClHCl_Ea_range_kj_mol` | `[90.0, 135.0]` |
| `replicates` | `3` |

**What it supports:**
- First-order thiamine degradation kinetics where TMN degrades faster (Ea 88-105 kJ/mol) than TClHCl (Ea 90-135 kJ/mol) due to pyrimidine protonation.

**What it does not support:**
- Solid-state high-moisture extrusion without pH validation

**Next Action:** Encode thiamine salt form Arrhenius kinetics.

---

### Mottram (1998)
- **DOI:** [10.1016/s0308-8146(98)00076-4](https://doi.org/10.1016/s0308-8146(98)00076-4)
- **Matrix Family:** methionine_model_system
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
| `methional_strecker_yield_mol_pct` | `38.0` |
| `temperature_c` | `150.0` |
| `amino_acid` | `methionine` |

**What it supports:**
- Methional Strecker yield anchor: ~38 mol% from methionine at ~150°C
- General sulfur volatile formation narrative from amino acid Strecker pathways

**What it does not support:**
- Plant-matrix methional quantification
- A universal Strecker yield under variable Aw or pH conditions

**Next Action:** Replaces the ghost mottram_wedzicha_1990_methional entry. Use as the methional Strecker yield anchor for Family 03 sulfur reporting.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Thiamine degradation and sulfur support` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `thiamine_availability`, `2-methyl-3-furanthiol` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 3

Total papers analyzed: **9** (Benchmark-eligible: **2**, Calibration references: **4**, Rejected: **3**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Voelker et al. (2021) | 4/8 | ⚠️ Calibration | aqueous_thiamine_model |
| Huang (2022) | 5/8 | ⚠️ Calibration | metal_catalyzed_thiamine_model |
| Arabshahi & Lund (1988), J. Food Sci. 53:199 | 4/8 | ⚠️ Calibration | starch_based_low_moisture_model |
| Cerny & Guntz-Dubini (2008) | 6/8 | ✅ Eligible | free_model_system |
| Bolton, Reineccius, Liardon & Huynh-Ba (1994), ACS Symp. Ser. 543:270-278 -- REPLACES the nonexistent "Hofmann et al. (1996)" | 6/8 | ✅ Eligible | free_model_system |
| Tang, Jiang, Yuan & Ho (2013), J. Sulfur Chem. 34:326 (CrossRef issue year 2012) | 2/8 | ❌ Rejected | free_model_system |
| Mauri, Alzamora, Chirife & Tomio (1989), Int. J. Food Sci. Technol. 24:1 | 2/8 | ❌ Rejected | free_model_system |
| Voelker et al. (2018) | 3/8 | ⚠️ Calibration | free_model_system |
| Mottram (1998) | 2/8 | ❌ Rejected | methionine_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 3

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Cerny & Guntz-Dubini (2008)` (Score: 6/8)
- `Bolton, Reineccius, Liardon & Huynh-Ba (1994), ACS Symp. Ser. 543:270-278 -- REPLACES the nonexistent "Hofmann et al. (1996)"` (Score: 6/8)

### CALIBRATION (Score 3-5)
- `Voelker et al. (2021)` (Score: 4/8)
- `Huang (2022)` (Score: 5/8)
- `Arabshahi & Lund (1988), J. Food Sci. 53:199` (Score: 4/8)
- `Voelker et al. (2018)` (Score: 3/8)

### REJECTED (Score < 3)
- `Tang, Jiang, Yuan & Ho (2013), J. Sulfur Chem. 34:326 (CrossRef issue year 2012)` (Score: 2/8)
- `Mauri, Alzamora, Chirife & Tomio (1989), Int. J. Food Sci. Technol. 24:1` (Score: 2/8)
- `Mottram (1998)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in  must explicitly account for `thiamine_fragmentation_support` as a modifier when predicting `thiamine_availability`, `2-methyl-3-furanthiol` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `high_value_support_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-10-27 using Maillard Ingestion Pipeline.*
