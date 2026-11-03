# Systematic Literature Review — Family 4: Nucleotide degradation and ribose support
**Last updated:** 2026-11-03  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering nucleotide degradation and ribose support.  
**Strategic Posture:** `high_value_support_lane`  
**Runtime Concept:** `nucleotide_and_ribose_support`  
**Scope & Targets:** Covers target compounds/variables: `IMP`, `GMP`, `ribose-5-phosphate`. Preferred payload types: `flavor_reference_payload`, `computational_prior`. Target runtime modules: .  

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

### Amado, Hanselman, Harmon, Deng et al. (2024), Chemical Senses 49:bjad049
- **DOI:** [10.1093/chemse/bjad049](https://doi.org/10.1093/chemse/bjad049)
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

### Aliani & Farmer (2005), J. Agric. Food Chem. 53:6067
- **DOI:** [10.1021/jf050085t](https://doi.org/10.1021/jf050085t)
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

### Finnigan, Wall, Wilde, Stephens et al. (2019), Curr. Dev. Nutr. 3:nzz021
- **DOI:** [10.1093/cdn/nzz021](https://doi.org/10.1093/cdn/nzz021)
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

Total papers analyzed: **3** (Benchmark-eligible: **0**, Calibration references: **1**, Rejected: **2**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Amado, Hanselman, Harmon, Deng et al. (2024), Chemical Senses 49:bjad049 | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Aliani & Farmer (2005), J. Agric. Food Chem. 53:6067 | 2/8 | ❌ Rejected | meat_model_system |
| Finnigan, Wall, Wilde, Stephens et al. (2019), Curr. Dev. Nutr. 3:nzz021 | 2/8 | ❌ Rejected | mycoprotein |

---

## Consolidated entries for benchmark_schema.json — Family 4

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Amado, Hanselman, Harmon, Deng et al. (2024), Chemical Senses 49:bjad049` (Score: 3/8)

### REJECTED (Score < 3)
- `Aliani & Farmer (2005), J. Agric. Food Chem. 53:6067` (Score: 2/8)
- `Finnigan, Wall, Wilde, Stephens et al. (2019), Curr. Dev. Nutr. 3:nzz021` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in  must explicitly account for `nucleotide_and_ribose_support` as a modifier when predicting `IMP`, `GMP`, `ribose-5-phosphate` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `high_value_support_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-11-03 using Maillard Ingestion Pipeline.*
