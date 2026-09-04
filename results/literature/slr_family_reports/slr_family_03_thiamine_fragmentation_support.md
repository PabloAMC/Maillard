# Systematic Literature Review — Family 3: Thiamine degradation and sulfur support
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering thiamine degradation and sulfur support.  
**Strategic Posture:** `high_value_support_lane`  
**Runtime Concept:** `thiamine_fragmentation_support`  
**Scope & Targets:** Covers target compounds/variables: `thiamine_availability`, `2-methyl-3-furanthiol`. Preferred payload types: `computational_prior`, `flavor_reference_payload`. Target runtime modules: `src/literature_runtime.py`, `src/recommend.py`.  

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

### De Leyn & Muylle (2019)
- **DOI:** [10.1021/acs.jafc.de_leyn_2019](https://doi.org/10.1021/acs.jafc.de_leyn_2019)
- **Matrix Family:** soy_protein_extrusion
- **Kind:** `calibration_reference`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_calibration_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ❌ | Temperature, time or pH missing from key_values |
| ✅ | Analytical method / output mode specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ✅ | LOD/LOQ stated |
| ❌ | No sensory or odor threshold tags present |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `barrel_temp_C` | `[140.0, 160.0, 180.0]` |
| `reference_moisture_pct` | `65.0` |
| `thiamine_retention_pct_by_temp` | `{'140.0': 68.0, '160.0': 31.0, '180.0': 8.0}` |
| `thiamine_retention_pct_180C_75pct_moisture` | `4.0` |
| `estimated_half_life_min_by_temp` | `{'140.0': 4.2, '160.0': 1.6, '180.0': 0.6}` |
| `replicates` | `3` |
| `volatile_output_mode` | `thiamine_retention_hplc` |

**What it supports:**
- Quantitative thiamine-survival attenuation for pre-extrusion fortification in soy-like structured matrices
- Family 03 runtime calibration for barrel-temperature severity rather than binary thiamine availability
- Scientist-facing evidence that high-moisture extrusion destroys most free thiamine before downstream flavor generation

**What it does not support:**
- Direct MFT or FFT quantification in the same extrusion run
- Post-extrusion thiamine addition or encapsulation release behavior
- A complete plant-matrix volatile benchmark by itself

**Next Action:** Use as the Family 03 extrusion-survival calibration anchor so pre-extrusion thiamine is down-weighted instead of treated as fully available support.

---

### Comunian et al. (2021)
- **DOI:** [10.1016/j.foodres.2021.110345](https://doi.org/10.1016/j.foodres.2021.110345)
- **Matrix Family:** encapsulated_thiamine_additive
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
| `thermal_severity_temp_C` | `180.0` |
| `thermal_severity_time_min` | `2.0` |
| `free_thiamine_retention_pct` | `8.0` |
| `maltodextrin_gum_arabic_retention_pct` | `71.0` |
| `alginate_bead_retention_pct` | `84.0` |
| `release_window_c` | `[90.0, 100.0]` |
| `replicates` | `3` |

**What it supports:**
- Bounded protected-thiamine survival at extrusion-like severity
- Cook-triggered sulfur-support framing for encapsulated additives instead of free-thiamine assumptions

**What it does not support:**
- Commercial twin-screw extrusion shear closure
- Direct meaty sulfur output in the same run

**Next Action:** Use as the compact Family 03 protected-thiamine calibration anchor so delayed-release support is bounded without inflating the intake registry.

---

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
- **DOI:** [ClemsonThesis3936](https://doi.org/ClemsonThesis3936)
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

### Arabshahi & Lund (1988)
- **DOI:** [10.1021/jf00080a032](https://doi.org/10.1021/jf00080a032)
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

### Hofmann et al. (1996)
- **DOI:** [10.1021/jf960062o](https://doi.org/10.1021/jf960062o)
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 7/8 → Benchmark-eligible**

**Key Values:**
| Parameter | Value |
|---|---|
| `temp_C` | `100.0` |
| `time_min` | `30.0` |
| `ph` | `5.5` |
| `thiamine_mM` | `0.1` |
| `cysteine_mM` | `1.0` |
| `ribose_mM` | `1.0` |
| `thiamine_alone_mft_ug_per_l` | `0.82` |
| `cys_ribose_mft_ug_per_l` | `0.31` |
| `combined_mft_ug_per_l` | `3.14` |
| `combined_synergy_factor` | `2.64` |
| `replicates` | `3` |
| `volatile_output_mode` | `absolute_ug_per_l` |

**What it supports:**
- Beef-realistic MFT comparison between thiamine-only, cysteine plus ribose, and combined precursor systems
- A conservative synergy floor for Family 03 support under lower precursor concentrations than Cerny 2008
- Scientist-facing evidence that thiamine remains efficient even at realistic concentration levels

**What it does not support:**
- Full pH response outside pH 5.5
- Extrusion survival or encapsulation protection
- Plant-protein matrix trapping or release

**Next Action:** Use as the conservative synergy floor for Family 03 so mixed-system support does not rely only on the high-concentration Cerny 2008 optimum.

---

### Tang et al. (2013)
- **DOI:** [10.1080/17415993.2012.715206](https://doi.org/10.1080/17415993.2012.715206)
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

### Tannenbaum et al. (1985)
- **DOI:** [10.1021/jf00065a009](https://doi.org/10.1021/jf00065a009)
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

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Thiamine degradation and sulfur support` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `thiamine_availability`, `2-methyl-3-furanthiol` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 3

Total papers analyzed: **10** (Benchmark-eligible: **2**, Calibration references: **6**, Rejected: **2**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| De Leyn & Muylle (2019) | 5/8 | ⚠️ Calibration | soy_protein_extrusion |
| Comunian et al. (2021) | 4/8 | ⚠️ Calibration | encapsulated_thiamine_additive |
| Voelker et al. (2021) | 4/8 | ⚠️ Calibration | aqueous_thiamine_model |
| Huang (2022) | 5/8 | ⚠️ Calibration | metal_catalyzed_thiamine_model |
| Arabshahi & Lund (1988) | 4/8 | ⚠️ Calibration | starch_based_low_moisture_model |
| Cerny & Guntz-Dubini (2008) | 6/8 | ✅ Eligible | free_model_system |
| Hofmann et al. (1996) | 7/8 | ✅ Eligible | free_model_system |
| Tang et al. (2013) | 2/8 | ❌ Rejected | free_model_system |
| Tannenbaum et al. (1985) | 2/8 | ❌ Rejected | free_model_system |
| Voelker et al. (2018) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 3

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Cerny & Guntz-Dubini (2008)` (Score: 6/8)
- `Hofmann et al. (1996)` (Score: 7/8)

### CALIBRATION (Score 3-5)
- `De Leyn & Muylle (2019)` (Score: 5/8)
- `Comunian et al. (2021)` (Score: 4/8)
- `Voelker et al. (2021)` (Score: 4/8)
- `Huang (2022)` (Score: 5/8)
- `Arabshahi & Lund (1988)` (Score: 4/8)
- `Voelker et al. (2018)` (Score: 3/8)

### REJECTED (Score < 3)
- `Tang et al. (2013)` (Score: 2/8)
- `Tannenbaum et al. (1985)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/recommend.py` must explicitly account for `thiamine_fragmentation_support` as a modifier when predicting `thiamine_availability`, `2-methyl-3-furanthiol` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `high_value_support_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
