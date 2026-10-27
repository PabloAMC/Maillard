# Systematic Literature Review — Family 9: Carbohydrate pyrolysis and caramelization
**Last updated:** 2026-10-27  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering carbohydrate pyrolysis and caramelization.  
**Strategic Posture:** `failure_mode_lane`  
**Runtime Concept:** `carbohydrate_pyrolysis_lane`  
**Scope & Targets:** Covers target compounds/variables: `HMF`, `furfural`, `caramelization_severity`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`, `computational_prior`. Target runtime modules: .  

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

### Hernandez et al. (2023), Molecules 28:3151
- **DOI:** [10.3390/molecules28073151](https://doi.org/10.3390/molecules28073151)
- **Matrix Family:** pbma_vs_beef_comparator
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
| ❌ | No sensory or odor threshold tags present |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `diacetyl_pbma_midpoint_ng_per_g` | `12.9` |
| `diacetyl_beef_midpoint_ng_per_g` | `38.1` |
| `acetoin_pbma_midpoint_ng_per_g` | `7.65` |
| `acetoin_beef_midpoint_ng_per_g` | `30.45` |
| `furfural_pbma_midpoint_ng_per_g` | `1040.0` |
| `furfural_beef_midpoint_ng_per_g` | `19.83` |
| `dimethylpyrazine_pbma_midpoint_ng_per_g` | `31.65` |
| `dimethylpyrazine_beef_midpoint_ng_per_g` | `15.72` |
| `acetylpyrrole_pbma_range_ng_per_g` | `[560.0, 726.0]` |
| `acetylpyrrole_beef_range_ng_per_g` | `[141.0, 170.0]` |
| `beef_flavour_identity_pbma_range` | `[11.5, 17.5]` |
| `beef_flavour_identity_beef_range` | `[46.8, 51.0]` |
| `diacetyl_deficit_fold_vs_beef` | `3.0` |
| `acetoin_deficit_fold_vs_beef` | `4.0` |
| `furfural_excess_fold_vs_beef` | `50.0` |
| `pyrrole_excess_fold_vs_beef_range` | `[4.0, 5.0]` |
| `replicates` | `4` |
| `volatile_output_mode` | `absolute_ng_per_g_cooked_product` |

**What it supports:**
- Executable external-literature benchmark subset for the PBMA-versus-beef browning gap using absolute furfural and 2,5-dimethylpyrazine values
- Scientist-facing benchmark framing for under-meaty and over-browned plant-based volatile signatures while the full comparator panel remains attached in intake metadata

**What it does not support:**
- Full executable closure of the diacetyl, acetoin, and 2-acetylpyrrole comparator panel in the current matrix runtime surface
- One-family mechanistic attribution of the full PBMA-versus-beef gap

**Next Action:** Keep the executable furfural plus dimethylpyrazine subset active while broader mixed-matrix comparator coverage for diacetyl, acetoin, and pyrroles is still unresolved.

---

### Glomb & Monnier (1995), J. Biol. Chem. 270:10017
- **DOI:** [10.1074/jbc.270.17.10017](https://doi.org/10.1074/jbc.270.17.10017)
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
| `methylglyoxal_fraction` | `0.41` |
| `glyoxal_fraction` | `0.28` |
| `diacetyl_fraction` | `0.18` |
| `2_3_pentanedione_fraction` | `0.06` |

**What it supports:**
- Deterministic 3-DG fragmentation fractions for MGO, GO, and diketones
- Bounded routing context for AGE and off-note reasoning

**What it does not support:**
- A complete kinetic mechanism
- Matrix-specific concentration closure

**Next Action:** Use as a bounded dicarbonyl fragmentation prior for AGE and diketone routing instead of opening a new benchmark lane.

---

### Quintas, Guimaraes, Baylina, Brandao & Silva (2007), Innov. Food Sci. Emerg. Technol. 8:306
- **DOI:** [10.1016/j.ifset.2007.02.002](https://doi.org/10.1016/j.ifset.2007.02.002)
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
| `hmf_mg_per_l_at_ph5_120c_60min` | `48.0` |

**What it supports:**
- Sucrose caramelisation k_obs and Ea versus pH and Temperature
- HMF generation of 48 mg/L at pH 5.0, 120 C after 60 min

**What it does not support:**
- Free amino acid condensation mechanisms
- Extruded texture impact parameters

**Next Action:** Encode sucrose caramelisation kinetics priors.

---

### de Bruijn, Kieboom & van Bekkum (1987), Recl. Trav. Chim. Pays-Bas 106:35
- **DOI:** [10.1002/recl.19871060201](https://doi.org/10.1002/recl.19871060201)
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
- Alkaline-driven retro-aldolization isomerization kinetics of hexose monosaccharides.

**What it does not support:**
- Acid-catalyzed sugar dehydration pathways.

**Next Action:** Encode alkaline retro-aldol prior.

---

### Luna & Aguilera (2013), Food Biophysics 9:61
- **DOI:** [10.1007/s11483-013-9317-0](https://doi.org/10.1007/s11483-013-9317-0)
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
- Arrhenius activation energies for color development in caramelizing molten hexose systems.

**What it does not support:**
- Low-temperature aqueous sugar degradation.

**Next Action:** Encode caramelization color kinetics prior.

---

### Cui et al. (2022), J. Agric. Food Chem. 70:11556
- **DOI:** [10.1021/acs.jafc.2c03427](https://doi.org/10.1021/acs.jafc.2c03427)
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
- Glycylglycine catalytic cleavage of Amadori products and 3-deoxyosone accumulation at 120 C.

**What it does not support:**
- Non-dipeptide free amino acid reactions.

**Next Action:** Encode dipeptide Amadori cleavage prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Carbohydrate pyrolysis and caramelization` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `HMF`, `furfural`, `caramelization_severity` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 9

Total papers analyzed: **6** (Benchmark-eligible: **0**, Calibration references: **4**, Rejected: **2**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Hernandez et al. (2023), Molecules 28:3151 | 5/8 | ⚠️ Calibration | pbma_vs_beef_comparator |
| Glomb & Monnier (1995), J. Biol. Chem. 270:10017 | 2/8 | ❌ Rejected | aqueous_model_system |
| Quintas, Guimaraes, Baylina, Brandao & Silva (2007), Innov. Food Sci. Emerg. Technol. 8:306 | 2/8 | ❌ Rejected | free_model_system |
| de Bruijn, Kieboom & van Bekkum (1987), Recl. Trav. Chim. Pays-Bas 106:35 | 3/8 | ⚠️ Calibration | free_model_system |
| Luna & Aguilera (2013), Food Biophysics 9:61 | 3/8 | ⚠️ Calibration | free_model_system |
| Cui et al. (2022), J. Agric. Food Chem. 70:11556 | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 9

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Hernandez et al. (2023), Molecules 28:3151` (Score: 5/8)
- `de Bruijn, Kieboom & van Bekkum (1987), Recl. Trav. Chim. Pays-Bas 106:35` (Score: 3/8)
- `Luna & Aguilera (2013), Food Biophysics 9:61` (Score: 3/8)
- `Cui et al. (2022), J. Agric. Food Chem. 70:11556` (Score: 3/8)

### REJECTED (Score < 3)
- `Glomb & Monnier (1995), J. Biol. Chem. 270:10017` (Score: 2/8)
- `Quintas, Guimaraes, Baylina, Brandao & Silva (2007), Innov. Food Sci. Emerg. Technol. 8:306` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in  must explicitly account for `carbohydrate_pyrolysis_lane` as a modifier when predicting `HMF`, `furfural`, `caramelization_severity` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `failure_mode_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-10-27 using Maillard Ingestion Pipeline.*
