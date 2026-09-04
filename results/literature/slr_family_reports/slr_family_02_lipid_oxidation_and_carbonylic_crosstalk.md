# Systematic Literature Review — Family 2: Lipid oxidation and carbonylic crosstalk
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering lipid oxidation and carbonylic crosstalk.  
**Strategic Posture:** `immediate_expansion_lane`  
**Runtime Concept:** `lipid_crosstalk_lane`  
**Scope & Targets:** Covers target compounds/variables: `hexanal`, `2-pentylfuran`, `nonanal`, `carbonyl_competition_factor`. Preferred payload types: `benchmark_payload`, `retention_payload`, `computational_prior`. Target runtime modules: `src/lipid_oxidation.py`, `src/literature_runtime.py`.  

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

### Trikusuma et al. (2019)
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

### Lincoln & Girard (2025)
- **DOI:** [10.1002/sfp2.1044](https://doi.org/10.1002/sfp2.1044)
- **Matrix Family:** ppi_spi_glucose_lipid_polyphenol
- **Kind:** `directional_prior`
- **Payload Role:** `benchmark_intake`
- **Status:** `ready_for_directional_prior_encoding`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output mode specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `ph` | `7.0` |
| `glucose_to_protein_w_w` | `1.0` |
| `corn_oil_wt_pct` | `0.55` |
| `spi_temp_C` | `90.0` |
| `ppi_temp_C` | `95.0` |
| `heating_time_h` | `1.0` |
| `storage_temp_C` | `37.0` |
| `storage_duration_days` | `14` |
| `polyphenol_examples` | `['catechin', 'tannic acid', 'green tea extract', 'grape seed extract']` |
| `volatile_output_mode` | `directional_volatiles` |

**What it supports:**
- Directional Strecker suppression prior for glucose plus lipids plus polyphenols in SPI and PPI systems
- Joint Maillard-lipid crosstalk framing for formulation reporting
- Scientist-facing provenance for polyphenol-driven Strecker moderation

**What it does not support:**
- Absolute concentration targets
- Strict quantitative optimizer calibration
- Benchmark promotion into the authoritative free-precursor surface

**Next Action:** Encode as a structured crosstalk prior and keep its directional status explicit until full quantitative extraction is available.

---

### Blank et al. (2001)
- **DOI:** [10.1021/jf000900i](https://doi.org/10.1021/jf000900i)
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

### Hrncirik & Zeelenberg (2014)
- **DOI:** [10.1016/j.foodchem.2013.11.082](https://doi.org/10.1016/j.foodchem.2013.11.082)
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

### Marquez-Ruiz et al. (2014)
- **DOI:** [10.1021/jf502636m](https://doi.org/10.1021/jf502636m)
- **Matrix Family:** high_oleic_oil_model_system
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `nonanal_oav_range` | `[8.0, 22.0]` |
| `hexanal_oav_ceiling` | `1.5` |
| `temperature_c` | `180.0` |
| `frying_cycle_max` | `20` |
| `replicates` | `3` |

**What it supports:**
- Compact high-oleic oil OAV anchor showing nonanal dominance with hexanal kept near or below threshold
- Scientist-facing routing context for why oleic-rich oils are materially safer than linoleic-rich oils for PBMA off-note control

**What it does not support:**
- Direct plant-matrix benchmark closure for one burger or extrudate
- A full fatty-acid oxidation kinetic mechanism across all frying conditions

**Next Action:** Keep as the oleic-rich reference anchor so oil-selection logic does not treat nonanal-dominant systems as equivalent to hexanal-heavy linoleic oxidation.

---

### Messina et al. (2022)
- **DOI:** [10.1016/j.foodchem.2021.131102](https://doi.org/10.1016/j.foodchem.2021.131102)
- **Matrix Family:** cooked_pbma_oil_comparison
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `canola_hexanal_oav` | `28.0` |
| `sunflower_high_oleic_nonanal_oav` | `12.0` |
| `coconut_delta_decalactone_oav` | `14.0` |
| `coconut_canola_3_1_hexanal_oav` | `4.0` |
| `replicates` | `6` |

**What it supports:**
- Comparative cooked-PBMA oil OAV panel for canola, high-oleic sunflower, coconut, and coconut:canola blend
- A compact reference surface for oil-choice ranking that keeps both green-rancid and sweet-lactone failure modes visible

**What it does not support:**
- Executable benchmark closure for one final commercial burger formulation
- A universal oil-ranking rule outside the reported SPI-based burger context

**Next Action:** Keep as the comparative oil-choice anchor so recommendation logic can distinguish high-oleic improvement from coconut-lactone drift without opening a new benchmark family.

---

### Hidalgo & Zamora (2004)
- **DOI:** [10.1021/jf035223w](https://doi.org/10.1021/jf035223w)
- **Matrix Family:** lipid_oxidation_model_system
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
| `target_offnote` | `2-pentylpyrrole` |
| `reactant_pairing` | `4-HNE + phenylalanine` |
| `absolute_concentrations_reported` | `True` |
| `loq_reported` | `True` |

**What it supports:**
- Bounded 4-HNE-to-2-pentylpyrrole off-note routing in oxidized lipid systems
- Absolute-concentration and LOQ-backed guardrail support for pyrrolic off-notes

**What it does not support:**
- Matrix-specific concentration closure
- A complete secondary lipid-oxidation mechanism

**Next Action:** Use as a bounded 4-HNE pyrrole guardrail so oxidized-lipid systems can surface pyrrolic off-note risk explicitly.

---

### Mottram et al. (2001)
- **DOI:** [10.1021/jf010391+](https://doi.org/10.1021/jf010391+)
- **Matrix Family:** sulfur_model_system
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
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `hexanal_to_cysteine_ratio_threshold` | `0.5` |
| `mft_quench_fraction_range` | `[0.15, 0.4]` |
| `optimal_carnosine_concentration_M` | `0.1` |
| `carnosine_mft_uplift_fraction_vs_unbuffered` | `0.61` |

**What it supports:**
- Competitive quenching of MFT by lipid aldehydes once hexanal pressure rises into the sulfur-limited regime
- A bounded carnosine-buffered process-state uplift for MFT-positive sulfur systems

**What it does not support:**
- Exact final-product closure in intact plant-protein matrices
- A complete lipid-Maillard kinetic mechanism across all aldehydes and buffers

**Next Action:** Keep the Mottram landing compact by splitting sulfur-quench and buffering support into one bounded prior plus one process-state anchor instead of opening a new benchmark family.

---

### Yeo & Mottram (2023)
- **DOI:** [10.1016/j.foodchem.2023.136009](https://doi.org/10.1016/j.foodchem.2023.136009)
- **Matrix Family:** plant_protein_plus_lecithin_model
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
| `optimal_soy_lecithin_pct_w_w_range` | `[0.5, 1.0]` |
| `alkyl_thiophene_uplift_fold` | `2.4` |
| `sensory_confirmation` | `True` |

**What it supports:**
- A bounded soy-lecithin dose window that lifts alkyl-thiophene meaty markers in plant-protein sulfur systems
- Compact phospholipid-crosstalk support for formulation ranking and reporting

**What it does not support:**
- Direct closure for any one commercial PBMA product
- A universal phospholipid-addition rule outside the reported model systems

**Next Action:** Treat soy lecithin as a bounded phospholipid-crosstalk calibration rather than as a generic oil addition or a new executable benchmark.

---

### Zhang et al. (2022)
- **DOI:** [10.1016/j.foodchem.2022.132009](https://doi.org/10.1016/j.foodchem.2022.132009)
- **Matrix Family:** aldehyde_model_system
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
| `unsaturation_effect_direction` | `unsaturated_greater_than_saturated` |
| `oav_reported` | `True` |
| `lod_loq_reported` | `True` |
| `triplicates` | `True` |

**What it supports:**
- A directional potency hierarchy where unsaturated aldehydes exert stronger crosstalk and off-note pressure than saturated analogues
- LOD/LOQ-backed off-note guardrail support for unsaturated aldehyde branches

**What it does not support:**
- Matrix-specific concentration closure
- A full aldehyde reaction-network mechanism

**Next Action:** Use as a directional off-note prior so unsaturated aldehyde branches stay visibly higher risk than saturated analogues without opening a new benchmark lane.

---

### Grosch (1982)
- **DOI:** [10.1021/jf00111a008](https://doi.org/10.1021/jf00111a008)
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

### Frankel et al. (1982)
- **DOI:** [10.1007/bf02534608](https://doi.org/10.1007/bf02534608)
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

### Grosch & Wieser (1999)
- **DOI:** [10.1021/jf990111u](https://doi.org/10.1021/jf990111u)
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
| `c183_to_propanal_beta_scission_yield_pct` | `8.2` |
| `c183_to_heptadienal_beta_scission_yield_pct` | `5.4` |

**What it supports:**
- C18:3 to propanal/2,4-heptadienal beta-scission yield splits.

**What it does not support:**
- C18:1/C18:2 scission crosstalk
- pH sensitivity

**Next Action:** Encode lipid c183 scission split yield prior.

---

### Farmer & Patterson (1991)
- **DOI:** [10.1016/0308-8146(91)90013-I](https://doi.org/10.1016/0308-8146(91)90013-I)
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

### Esterbauer et al. (1991)
- **DOI:** [10.1016/0891-5849(91)90081-A](https://doi.org/10.1016/0891-5849(91)90081-A)
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

### Kamal-Eldin et al. (2003)
- **DOI:** [10.1007/s11745-003-1082-x](https://doi.org/10.1007/s11745-003-1082-x)
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

### Xu et al. (2024)
- **DOI:** [10.3390/foods13091393](https://doi.org/10.3390/foods13091393)
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
| `soy_pbma_storage_hexanal_ppb` | `185.0` |
| `storage_temp_C` | `4.0` |
| `storage_time_days` | `14.0` |

**What it supports:**
- Soybean PBMA storage hexanal accumulation trends and baseline concentration values.

**What it does not support:**
- Free solution kinetic activation energies
- Furanone reaction branching ratios

**Next Action:** Expose as a soybean PBMA storage flavor reference payload.

---

### Yeo & Shibamoto (1991)
- **DOI:** [10.1021/jf00002a024](https://doi.org/10.1021/jf00002a024)
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
| `phospholipid_hexanal_fraction_pct_range` | `[82.0, 88.0]` |
| `lox_d_value_70C_min` | `28.0` |
| `lox_z_value_C` | `10.0` |

**What it supports:**
- Phospholipid fraction contribution (82-88%) to warmed-over flavor (WOF) hexanal.
- LOX denaturation D-value table (70 C D=28 min, z~10 C).

**What it does not support:**
- Direct peptide-cysteine nucleophilicity indices
- Absolute volatile release in non-meat systems

**Next Action:** Encode phospholipid hexanal contribution and LOX denaturation priors.

---

### Tan et al. (2001)
- **DOI:** [10.1007/s11746-001-0402-y](https://doi.org/10.1007/s11746-001-0402-y)
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

### Bayram et al. (2023)
- **DOI:** [10.3989/gya.1051211](https://doi.org/10.3989/gya.1051211)
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

### Spier et al. (2010)
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

Total papers analyzed: **21** (Benchmark-eligible: **1**, Calibration references: **13**, Rejected: **7**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Trikusuma et al. (2019) | 6/8 | ✅ Eligible | pea_uht_beverage |
| Lincoln & Girard (2025) | 5/8 | ⚠️ Calibration | ppi_spi_glucose_lipid_polyphenol |
| Blank et al. (2001) | 3/8 | ⚠️ Calibration | lipid_rich_model_system |
| Hrncirik & Zeelenberg (2014) | 3/8 | ⚠️ Calibration | coconut_oil_co_matrix |
| Marquez-Ruiz et al. (2014) | 4/8 | ⚠️ Calibration | high_oleic_oil_model_system |
| Messina et al. (2022) | 4/8 | ⚠️ Calibration | cooked_pbma_oil_comparison |
| Hidalgo & Zamora (2004) | 3/8 | ⚠️ Calibration | lipid_oxidation_model_system |
| Mottram et al. (2001) | 4/8 | ⚠️ Calibration | sulfur_model_system |
| Yeo & Mottram (2023) | 3/8 | ⚠️ Calibration | plant_protein_plus_lecithin_model |
| Zhang et al. (2022) | 3/8 | ⚠️ Calibration | aldehyde_model_system |
| Grosch (1982) | 3/8 | ⚠️ Calibration | free_model_system |
| Frankel et al. (1982) | 2/8 | ❌ Rejected | free_model_system |
| Grosch & Wieser (1999) | 2/8 | ❌ Rejected | free_model_system |
| Farmer & Patterson (1991) | 2/8 | ❌ Rejected | free_model_system |
| Esterbauer et al. (1991) | 2/8 | ❌ Rejected | free_model_system |
| Kamal-Eldin et al. (2003) | 2/8 | ❌ Rejected | free_model_system |
| Xu et al. (2024) | 2/8 | ❌ Rejected | soy_isolate |
| Yeo & Shibamoto (1991) | 2/8 | ❌ Rejected | free_model_system |
| Tan et al. (2001) | 3/8 | ⚠️ Calibration | free_model_system |
| Bayram et al. (2023) | 3/8 | ⚠️ Calibration | free_model_system |
| Spier et al. (2010) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 2

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Trikusuma et al. (2019)` (Score: 6/8)

### CALIBRATION (Score 3-5)
- `Lincoln & Girard (2025)` (Score: 5/8)
- `Blank et al. (2001)` (Score: 3/8)
- `Hrncirik & Zeelenberg (2014)` (Score: 3/8)
- `Marquez-Ruiz et al. (2014)` (Score: 4/8)
- `Messina et al. (2022)` (Score: 4/8)
- `Hidalgo & Zamora (2004)` (Score: 3/8)
- `Mottram et al. (2001)` (Score: 4/8)
- `Yeo & Mottram (2023)` (Score: 3/8)
- `Zhang et al. (2022)` (Score: 3/8)
- `Grosch (1982)` (Score: 3/8)
- `Tan et al. (2001)` (Score: 3/8)
- `Bayram et al. (2023)` (Score: 3/8)
- `Spier et al. (2010)` (Score: 3/8)

### REJECTED (Score < 3)
- `Frankel et al. (1982)` (Score: 2/8)
- `Grosch & Wieser (1999)` (Score: 2/8)
- `Farmer & Patterson (1991)` (Score: 2/8)
- `Esterbauer et al. (1991)` (Score: 2/8)
- `Kamal-Eldin et al. (2003)` (Score: 2/8)
- `Xu et al. (2024)` (Score: 2/8)
- `Yeo & Shibamoto (1991)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/lipid_oxidation.py`, `src/literature_runtime.py` must explicitly account for `lipid_crosstalk_lane` as a modifier when predicting `hexanal`, `2-pentylfuran`, `nonanal`, `carbonyl_competition_factor` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `immediate_expansion_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
