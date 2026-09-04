# Systematic Literature Review — Family 8: Plant off-notes and Maillard suppression
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering plant off-notes and maillard suppression.  
**Strategic Posture:** `guardrail_lane`  
**Runtime Concept:** `off_note_and_maillard_suppression`  
**Scope & Targets:** Covers target compounds/variables: `dicarbonyl_trapping_factor`, `amino_group_blocking_factor`. Preferred payload types: `benchmark_payload`, `computational_prior`, `safety_payload`. Target runtime modules: `src/literature_runtime.py`, `src/safety.py`.  

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

### Squeo et al. (2023)
- **DOI:** [10.3390/foods12061331](https://doi.org/10.3390/foods12061331)
- **Matrix Family:** commercial_pbpi
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `soy_wet_extraction_acrylamide_ug_per_kg` | `[185, 748]` |
| `wet_extraction_mean_ug_per_kg` | `451` |
| `replicates` | `3` |
| `lod_ng_per_ml` | `7` |
| `loq_ng_per_ml` | `24` |

**What it supports:**
- Primary quantitative endpoint anchor for acrylamide in commercial soy and pea-derived protein ingredients
- Upper-bound safety reference for processed SPI/pea systems

**What it does not support:**
- Kinetic reconstruction across temperature, time, and free asparagine
- A free-precursor or controlled-process acrylamide benchmark

**Next Action:** Expose as a safety-reference payload without misclassifying it as a kinetic benchmark.

---

### Wang et al. (2022)
- **DOI:** [10.1016/j.lwt.2022.113009](https://doi.org/10.1016/j.lwt.2022.113009)
- **Matrix Family:** lab_fermented_plant_protein
- **Kind:** `calibration_reference`
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
| ✅ | Odor thresholds / sensory reported |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `hexanal_depletion_fraction` | `0.957` |
| `residual_hexanal_oav` | `0.47` |
| `adh_aldh_mechanism_confirmed` | `True` |
| `sensory_panel_reported` | `True` |

**What it supports:**
- A compact residual-hexanal OAV target after LAB cleanup
- Scientist-facing reporting support for pushing fermentation-treated systems below the hexanal threshold

**What it does not support:**
- A universal fermentation kinetics model
- Exact final-product closure across all plant proteins and starter cultures

**Next Action:** Keep Wang 2022 as a compact hexanal-cleanup reference anchor that complements, rather than duplicates, the already-landed fermentation cleanup state package.

---

### Bhandari et al. (1998)
- **DOI:** [10.1016/S0260-8774(98)00085-4](https://doi.org/10.1016/S0260-8774(98)00085-4)
- **Matrix Family:** aqueous_encapsulation_model_system
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
| `beta_cyclodextrin_loading_wt_pct_reference` | `1.0` |
| `hexanal_headspace_reduction_fraction` | `0.62` |
| `nonanal_headspace_reduction_fraction` | `0.68` |
| `e_2_nonenal_headspace_reduction_fraction` | `0.74` |
| `hexanal_k_assoc_m_inverse` | `1840.0` |

**What it supports:**
- A bounded beta-cyclodextrin sequestration surface for hexanal, nonanal, and (E)-2-nonenal
- Multi-hurdle off-note reasoning that treats cyclodextrin as a partial suppressor rather than a complete fix

**What it does not support:**
- A standalone encapsulation solver
- Finished-product sensory closure without complementary off-note mitigation

**Next Action:** Encode as a bounded cyclodextrin-binding prior so Family 08 can expose partial aldehyde suppression without conflating it with protein-binding retention.

---

### Liu (2022)
- **DOI:** [10.1016/j.foodchem.2022.134998](https://doi.org/10.1016/j.foodchem.2022.134998)
- **Matrix Family:** pea_isolate
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
| ✅ | Odor thresholds / sensory reported |

**Score: 3/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `methoxypyrazine_odt_ppb_range` | `[0.002, 0.006]` |
| `volatiles_count` | `12` |

**What it supports:**
- Baseline OAV table for 12 key pea protein isolate (PPI) volatiles
- Methoxypyrazines odor detection threshold (ODT) of 0.002-0.006 ppb

**What it does not support:**
- Free solution kinetic activation energies
- Browning index evolution rates

**Next Action:** Expose as a pea isolate off-note flavor reference payload.

---

### Rawel et al. (2002)
- **DOI:** [10.1021/jf020082z](https://doi.org/10.1021/jf020082z)
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
| `cys_blocking_rate_vs_lys` | `2.9` |
| `lys_blocked_pct_at_10_1` | `43.0` |

**What it supports:**
- Chlorogenic acid (CGA) Cys blocking k_obs is 2.9x faster than Lys
- Lysine blocking efficiency of 43% at a 10:1 polyphenol-to-protein ratio

**What it does not support:**
- Matrix-assisted cysteine release rates
- Radical scavenging mechanism validation

**Next Action:** Encode CGA Cys and Lys blocking rate priors.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Plant off-notes and Maillard suppression` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `dicarbonyl_trapping_factor`, `amino_group_blocking_factor` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 8

Total papers analyzed: **5** (Benchmark-eligible: **0**, Calibration references: **4**, Rejected: **1**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Squeo et al. (2023) | 3/8 | ⚠️ Calibration | commercial_pbpi |
| Wang et al. (2022) | 3/8 | ⚠️ Calibration | lab_fermented_plant_protein |
| Bhandari et al. (1998) | 3/8 | ⚠️ Calibration | aqueous_encapsulation_model_system |
| Liu (2022) | 3/8 | ⚠️ Calibration | pea_isolate |
| Rawel et al. (2002) | 2/8 | ❌ Rejected | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 8

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Squeo et al. (2023)` (Score: 3/8)
- `Wang et al. (2022)` (Score: 3/8)
- `Bhandari et al. (1998)` (Score: 3/8)
- `Liu (2022)` (Score: 3/8)

### REJECTED (Score < 3)
- `Rawel et al. (2002)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/safety.py` must explicitly account for `off_note_and_maillard_suppression` as a modifier when predicting `dicarbonyl_trapping_factor`, `amino_group_blocking_factor` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `guardrail_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
