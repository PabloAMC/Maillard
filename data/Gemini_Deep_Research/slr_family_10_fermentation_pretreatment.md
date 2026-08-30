# Systematic Literature Review — Family 10: Microbial fermentation pretreatment
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering microbial fermentation pretreatment.  
**Strategic Posture:** `upstream_pretreatment_lane`  
**Runtime Concept:** `fermentation_pretreatment_node`  
**Scope & Targets:** Covers target compounds/variables: `free_amino_acid_enrichment`, `nucleotide_enrichment`, `pretreatment_pH_shift`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`, `computational_prior`. Target runtime modules: `src/literature_runtime.py`, `src/recommend.py`, `src/matrix_prior_registry.py`.  

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

### Rizzello et al. (2024)
- **DOI:** [10.1021/acs.jafc.3c08432](https://doi.org/10.1021/acs.jafc.3c08432)
- **Matrix Family:** fermented_plant_based_matrix
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
| `hexanal_native_ug_per_g` | `28.4` |
| `hexanal_fermented_ug_per_g` | `1.42` |
| `mft_fermented_ug_per_g` | `3.82` |
| `mft_uplift_fold` | `382.0` |
| `sensory_panel_n` | `15` |

**What it supports:**
- Fermentation-driven hexanal cleanup prior for plant-based meat analogues
- Bounded sulfur uplift context after LAB fermentation

**What it does not support:**
- Executable benchmark closure for a single matrix benchmark
- General fermentation kinetics outside the reported cleanup window

**Next Action:** Keep encoded as a process-state calibration for fermentation cleanup and sulfur-positive uplift rather than as a benchmark payload.

---

### Zhao et al. (2022)
- **DOI:** [10.3390/molecules27196182](https://doi.org/10.3390/molecules27196182)
- **Matrix Family:** moromi_like_fermentation
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
| `glutamate_multiplier_range` | `[1.99, 2.49]` |
| `glycine_multiplier_range` | `[2.38, 3.73]` |
| `gmp_multiplier_range` | `[3.02, 6.61]` |
| `amino_type_nitrogen_initial_mg_pct` | `1175.0` |
| `amino_type_nitrogen_final_mg_pct` | `2221.0` |

**What it supports:**
- Fermentation-driven amino-acid and nucleotide release multipliers
- Process-state context for moromi-like precursor enrichment

**What it does not support:**
- Absolute volatile endpoint closure
- Universal pretreatment transfer outside koji or moromi systems

**Next Action:** Keep encoded as fermentation-state metadata for donor release rather than as a standalone flavor benchmark.

---

### Perdiguero et al. (2004)
- **DOI:** [10.1021/jf0494452](https://doi.org/10.1021/jf0494452)
- **Matrix Family:** yeast_ye
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
| `gmp_mg_per_g_dw` | `2.8` |
| `temp_C` | `55.0` |
| `time_h` | `24.0` |

**What it supports:**
- Yeast autolysis nucleotide release kinetics yielding GMP 2.8 mg/g DW at 55 C after 24h

**What it does not support:**
- Free amino acid extraction scaling parameters
- Browning intensity at high temperature

**Next Action:** Expose as a fermentation process state calibration payload.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Microbial fermentation pretreatment` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `free_amino_acid_enrichment`, `nucleotide_enrichment`, `pretreatment_pH_shift` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 10

Total papers analyzed: **3** (Benchmark-eligible: **0**, Calibration references: **2**, Rejected: **1**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Rizzello et al. (2024) | 4/8 | ⚠️ Calibration | fermented_plant_based_matrix |
| Zhao et al. (2022) | 3/8 | ⚠️ Calibration | moromi_like_fermentation |
| Perdiguero et al. (2004) | 2/8 | ❌ Rejected | yeast_ye |

---

## Consolidated entries for benchmark_schema.json — Family 10

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Rizzello et al. (2024)` (Score: 4/8)
- `Zhao et al. (2022)` (Score: 3/8)

### REJECTED (Score < 3)
- `Perdiguero et al. (2004)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/recommend.py`, `src/matrix_prior_registry.py` must explicitly account for `fermentation_pretreatment_node` as a modifier when predicting `free_amino_acid_enrichment`, `nucleotide_enrichment`, `pretreatment_pH_shift` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `upstream_pretreatment_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
