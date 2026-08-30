# Systematic Literature Review — Family 4: Nucleotide degradation and ribose support
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering nucleotide degradation and ribose support.  
**Strategic Posture:** `high_value_support_lane`  
**Runtime Concept:** `nucleotide_and_ribose_support`  
**Scope & Targets:** Covers target compounds/variables: `IMP`, `GMP`, `ribose-5-phosphate`. Preferred payload types: `flavor_reference_payload`, `computational_prior`. Target runtime modules: `src/literature_runtime.py`, `src/recommend.py`.  

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

### Tanaka et al. (2025)
- **DOI:** [10.1093/chemse/bjaf043](https://doi.org/10.1093/chemse/bjaf043)
- **Matrix Family:** aqueous_model_system
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
| ✅ | Odor thresholds / sensory reported |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `nucleotide_concentration_mM` | `3.0` |
| `imp_fold_reduction_vs_msg` | `45.2` |
| `gmp_fold_reduction_vs_msg` | `29.8` |
| `amp_fold_reduction_vs_msg` | `8.9` |
| `cmp_effect_direction` | `negative_modulator` |
| `panelists` | `30` |
| `assay_mode` | `2_afc_threshold` |

**What it supports:**
- Quantitative fold-reduction anchor for IMP and GMP as upstream umami support markers
- A family-04 runtime distinction between real nucleotide support and non-potentiating nucleotide-like cues
- Scientist-facing support that small nucleotide additions can matter even when they are not treated as direct volatile precursors

**What it does not support:**
- Thermal degradation or ribose release kinetics
- PBMA-specific matrix trapping or retention
- A complete EUC reconstruction for every formulation without concentration data

**Next Action:** Use as the Family 04 umami-support anchor so IMP and GMP are scored as bounded upstream support rather than generic additive cues.

---

### Yamaguchi & Ninomiya (2000)
- **DOI:** [10.1093/jn/130.4.921S](https://doi.org/10.1093/jn/130.4.921S)
- **Matrix Family:** aqueous_and_food_matrix_review
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
| `umami_synergy_constant` | `1218.0` |
| `imp_relative_umami_unit` | `1.0` |
| `gmp_relative_umami_unit` | `2.3` |
| `amp_relative_umami_unit` | `0.18` |
| `validated_food_matrices` | `15` |

**What it supports:**
- The EUC synergy constant used to interpret IMP or GMP as multiplicative umami support rather than additive flavor mass
- A conservative GMP relative-potency anchor for Family 04 scoring
- Scientist-facing explanation for why low-dose nucleotide enrichment still matters

**What it does not support:**
- Thermal hydrolysis or ribose release kinetics
- Matrix-specific nucleotide retention
- A direct volatile benchmark

**Next Action:** Carry the EUC constant into the bounded Family 04 prior so the lane reports preserved umami support explicitly.

---

### Soladoye et al. (2020)
- **DOI:** [10.3390/foods9030251](https://doi.org/10.3390/foods9030251)
- **Matrix Family:** animal_meat_low_temperature_reference
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
| `raw_euc_percent_msg` | `0.18` |
| `imp_mg_per_100g_by_condition` | `{'60C_2h': 38.0, '60C_12h': 6.0, '70C_12h': 51.0}` |
| `glutamate_mg_per_100g_by_condition` | `{'60C_2h': 14.0, '60C_12h': 18.0, '70C_12h': 15.0}` |
| `euc_percent_msg_by_condition` | `{'60C_2h': 0.13, '60C_12h': 0.04, '70C_12h': 0.15}` |

**What it supports:**
- A low-temperature Family 04 reference window for preserved IMP-driven EUC under long mild heating
- Scientist-facing context that 60 C and 70 C holding do not collapse into the same nucleotide-support outcome
- A bounded rationale for PBMA post-extrusion nucleotide addition or protection without claiming direct transfer from beef

**What it does not support:**
- Direct PBMA endpoint closure under extrusion or slurry conditions
- Universal low-temperature kinetics for IMP and GMP across plant matrices
- A standalone ribose-release model

**Next Action:** Use as a Family 04 mild-heating reference so long low-temperature windows report preserved-versus-collapsed EUC context explicitly while keeping PBMA transfer bounded.

---

### Ahlberg & Mohammadi (2021)
- **DOI:** [10.4014/jmb.2207.07057](https://doi.org/10.4014/jmb.2207.07057)
- **Matrix Family:** commercial_yeast_extract
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
| `standard_autolysate_imp_mg_per_100g_dw_range` | `[120.0, 250.0]` |
| `standard_autolysate_gmp_mg_per_100g_dw_range` | `[45.0, 90.0]` |
| `high_nucleotide_imp_mg_per_100g_dw_range` | `[1200.0, 2400.0]` |
| `high_nucleotide_gmp_mg_per_100g_dw_range` | `[300.0, 650.0]` |
| `amp_deaminase_euc_uplift_factor_vs_standard_range` | `[5.0, 6.0]` |

**What it supports:**
- A bounded Family 04 grade window for yeast-extract nucleotide support instead of treating all yeast extracts as compositionally equivalent
- Scientist-facing context that AMP-deaminase treatment can convert weak-umami autolysates into high-IMP support inputs

**What it does not support:**
- A universal grade assignment for any yeast extract label that lacks compositional detail
- Direct PBMA endpoint closure for combined nucleotide plus glutamate systems

**Next Action:** Use as a Family 04 source-profile anchor so yeast-extract additions expose bounded grade-driven nucleotide support rather than behaving like a single default composition.

---

### Cui et al. (2022)
- **DOI:** [10.31665/JFB.2022.18309](https://doi.org/10.31665/JFB.2022.18309)
- **Matrix Family:** commercial_dried_mushroom
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
| `shiitake_gmp_mg_per_100g_dw` | `108.0` |
| `shiitake_total_euc_g_msg_per_100g` | `2.1` |
| `porcini_gmp_mg_per_100g_dw` | `82.0` |
| `porcini_total_euc_g_msg_per_100g` | `1.8` |
| `enoki_gmp_mg_per_100g_dw` | `45.0` |
| `oyster_gmp_mg_per_100g_dw` | `28.0` |

**What it supports:**
- A bounded Family 04 clean-label nucleotide support anchor for shiitake, porcini, enoki, and oyster powders
- Species-specific GMP and EUC context so mushroom additions can be interpreted as more than a generic savoury cue

**What it does not support:**
- Combined PBMA endpoint closure for mushroom plus yeast-extract blends
- Thermal retention kinetics for mushroom nucleotides under extrusion or consumer cooking

**Next Action:** Use as a Family 04 source-profile anchor so mushroom additions expose bounded species-specific nucleotide support, especially when shiitake or porcini are present.

---

### Matoba et al. (1988)
- **DOI:** [10.1111/j.1365-2621.1988.tb13551.x](https://doi.org/10.1111/j.1365-2621.1988.tb13551.x)
- **Matrix Family:** aqueous_model_system
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `imp_half_life_h_at_100c_by_ph` | `{'4.0': 8.7, '5.0': 12.0, '7.0': 13.1, '9.0': 46.2}` |
| `gmp_half_life_h_at_100c_by_ph` | `{'4.0': 6.4, '5.0': 9.0, '7.0': 8.2, '9.0': 38.5}` |
| `imp_half_life_h_at_ph7_by_temp` | `{'100.0': 13.1, '110.0': 4.4, '130.0': 0.5, '150.0': 0.0583}` |
| `gmp_half_life_h_at_ph7_by_temp` | `{'100.0': 8.2, '110.0': 2.7, '130.0': 0.3, '150.0': 0.0367}` |
| `replicates` | `3` |
| `assay_mode` | `hplc_uv` |

**What it supports:**
- Quantitative IMP and GMP survival under heating and extrusion-like severity
- A benchmark-facing tradeoff between preserved umami support and hydrolysis-driven ribose delivery
- Separate treatment of IMP and GMP instead of a pooled nucleotide bucket

**What it does not support:**
- Plant-matrix retention or encapsulation
- Direct MFT or FFT yields in the same experiment
- Complete ribose liberation kinetics after inosine cleavage

**Next Action:** Use as the Family 04 thermal-severity anchor so nucleotide support survives mildly but shifts toward ribose delivery under extrusion-like conditions.

---

### Epstein et al. (1988)
- **DOI:** [10.1021/jf00083a026](https://doi.org/10.1021/jf00083a026)
- **Matrix Family:** aqueous_and_food_matrix_comparison
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
| `imp_to_inosine_k_h_inv_at_100c_ph6_5` | `0.053` |
| `inosine_to_ribose_k_h_inv_at_100c_ph6_5` | `0.018` |
| `effective_ribose_release_fraction_at_100c_30min` | `0.06` |
| `dominant_constraint` | `n_glycosidic_cleavage_rate_limiting` |

**What it supports:**
- A bounded ribose-release ceiling so Family 04 does not over-convert hydrolyzed IMP into immediate pentose delivery
- Scientist-facing context that inosine cleavage remains the slow step under consumer-cook conditions
- A distinction between preserved nucleotide pool and actual ribose-delivery shift

**What it does not support:**
- Extrusion residence-time reconstruction
- IMP or GMP sensory synergy on its own
- Matrix-specific nucleotide retention

**Next Action:** Use as the bounded ribose-delivery anchor so Family 04 adds only modest ribose under mild cooking and larger support only under extrusion-like severity.

---

### Blank & Grosch (1991)
- **DOI:** [10.1021/jf00009a012](https://doi.org/10.1021/jf00009a012)
- **Matrix Family:** boiled_beef_reference
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `hdmf_ug_per_kg_range` | `[28.0, 67.0]` |
| `hdmf_odt_ug_per_kg` | `6.0` |
| `hdmf_oav_range` | `[4.7, 11.2]` |
| `fd_factor` | `256.0` |
| `replicates` | `3` |

**What it supports:**
- Beef-relevant HDMF target band for brothy-caramel support
- OAV-bounded donor-routing support showing HDMF is sensorially active well below overt caramelization overload

**What it does not support:**
- A full PBMA benchmark by itself
- Direct proof that endogenous IMP is required for HDMF closure in plant matrices

**Next Action:** Use as the compact Family 04 donor-routing and reporting anchor for HDMF without adding a large inline benchmark block.

---

### Liardon et al. (1991)
- **DOI:** [10.1021/acs.jafc.liardon_1991](https://doi.org/10.1021/acs.jafc.liardon_1991)
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `r5p_vs_glucose_multiplier_floor` | `100.0` |
| `ribose_vs_glucose_multiplier` | `7.2` |
| `r5p_reaction_order` | `1.5` |

**What it supports:**
- Large donor-strength separation between R5P, ribose, and glucose
- A bounded reaction-order prior for high-potency reducing sugar routing

**What it does not support:**
- Matrix-specific retention or trapping
- A direct plant benchmark

**Next Action:** Encode as a donor-potency prior that keeps R5P and ribose above glucose in sulfur-positive routing decisions.

---

### Aliani & Farmer (2005)
- **DOI:** [10.1021/jf0401831](https://doi.org/10.1021/jf0401831)
- **Matrix Family:** meat_model_system
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
| `ribose_mft_multiplier` | `3.8` |
| `glucose_6_phosphate_mft_multiplier` | `3.2` |
| `glucose_mft_multiplier` | `1.4` |
| `nucleotide_pool_measured` | `True` |

**What it supports:**
- Bounded donor-priority multipliers for ribose, G6P, and glucose in sulfur-positive chemistry
- Nucleotide-context support for interpreting ribose-rich cooked systems

**What it does not support:**
- Plant-matrix closure
- A universal sulfur yield multiplier

**Next Action:** Encode as a bounded donor-potency prior and keep the nucleotide context explicit in reporting.

---

### Finnigan et al. (2019)
- **DOI:** [10.1039/c9fo00878e](https://doi.org/10.1039/c9fo00878e)
- **Matrix Family:** mycoprotein
- **Kind:** `calibration_reference`
- **Payload Role:** `process_state_calibration`
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
| `rna_reduction_temp_C` | `65.0` |
| `residual_nucleotides_dry_wt_pct` | `1.5` |

**What it supports:**
- Mycoprotein RNA reduction process parameters
- Residual IMP and GMP concentrations at 60-70 C and EUC impact mapping

**What it does not support:**
- Volatile release in soy matrices

**Next Action:** Expose as mycoprotein RNA reduction process state calibration.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Nucleotide degradation and ribose support` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `IMP`, `GMP`, `ribose-5-phosphate` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 4

Total papers analyzed: **11** (Benchmark-eligible: **0**, Calibration references: **8**, Rejected: **3**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Tanaka et al. (2025) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Yamaguchi & Ninomiya (2000) | 3/8 | ⚠️ Calibration | aqueous_and_food_matrix_review |
| Soladoye et al. (2020) | 3/8 | ⚠️ Calibration | animal_meat_low_temperature_reference |
| Ahlberg & Mohammadi (2021) | 3/8 | ⚠️ Calibration | commercial_yeast_extract |
| Cui et al. (2022) | 3/8 | ⚠️ Calibration | commercial_dried_mushroom |
| Matoba et al. (1988) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Epstein et al. (1988) | 3/8 | ⚠️ Calibration | aqueous_and_food_matrix_comparison |
| Blank & Grosch (1991) | 3/8 | ⚠️ Calibration | boiled_beef_reference |
| Liardon et al. (1991) | 2/8 | ❌ Rejected | aqueous_model_system |
| Aliani & Farmer (2005) | 2/8 | ❌ Rejected | meat_model_system |
| Finnigan et al. (2019) | 2/8 | ❌ Rejected | mycoprotein |

---

## Consolidated entries for benchmark_schema.json — Family 4

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Tanaka et al. (2025)` (Score: 3/8)
- `Yamaguchi & Ninomiya (2000)` (Score: 3/8)
- `Soladoye et al. (2020)` (Score: 3/8)
- `Ahlberg & Mohammadi (2021)` (Score: 3/8)
- `Cui et al. (2022)` (Score: 3/8)
- `Matoba et al. (1988)` (Score: 3/8)
- `Epstein et al. (1988)` (Score: 3/8)
- `Blank & Grosch (1991)` (Score: 3/8)

### REJECTED (Score < 3)
- `Liardon et al. (1991)` (Score: 2/8)
- `Aliani & Farmer (2005)` (Score: 2/8)
- `Finnigan et al. (2019)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/recommend.py` must explicitly account for `nucleotide_and_ribose_support` as a modifier when predicting `IMP`, `GMP`, `ribose-5-phosphate` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `high_value_support_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
