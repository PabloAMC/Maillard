# Systematic Literature Review — Family 9: Carbohydrate pyrolysis and caramelization
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering carbohydrate pyrolysis and caramelization.  
**Strategic Posture:** `failure_mode_lane`  
**Runtime Concept:** `carbohydrate_pyrolysis_lane`  
**Scope & Targets:** Covers target compounds/variables: `HMF`, `furfural`, `caramelization_severity`. Preferred payload types: `benchmark_payload`, `flavor_reference_payload`, `computational_prior`. Target runtime modules: `src/literature_runtime.py`, `src/recommend.py`.  

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

### Resconi et al. (2023)
- **DOI:** [10.3390/foods12071511](https://doi.org/10.3390/foods12071511)
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

### Mottram & Nobrega (2002)
- **DOI:** [10.1021/jf0200826](https://doi.org/10.1021/jf0200826)
- **Matrix Family:** caramelization_sulfur_model_system
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
| `norfuraneol_to_mft_mol_fraction` | `0.014` |
| `hdmf_cysteine_thiophene_ph_optimum` | `5.1` |
| `diacetyl_mercaptoketone_bridge_reported` | `True` |

**What it supports:**
- Bounded bridge from norfuraneol and HDMF-like intermediates into sulfur-positive notes
- Scientist-facing explanation for why caramel-heavy systems can still retain sulfur-side opportunity when sulfur donors are present

**What it does not support:**
- Exact MFT or thiophene endpoint closure in one PBMA matrix
- A standalone sulfur benchmark

**Next Action:** Keep encoded as a bounded sulfur-bridge prior so furanone-rich systems do not get treated as sulfur-disconnected by default.

---

### Ordoudi et al. (2014)
- **DOI:** [10.1016/j.foodchem.2014.05.097](https://doi.org/10.1016/j.foodchem.2014.05.097)
- **Matrix Family:** caramelization_model_system
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
| `hmf_peak_mg_per_l` | `140.0` |
| `hmf_residual_mg_per_l_at_360_min` | `20.0` |
| `secondary_decline_horizon_min` | `360.0` |
| `amino_acid_cocatalysis_multiplier` | `7.8` |

**What it supports:**
- An explicit HMF peak-and-decline process-state window at pH 5.0 and 125 C
- Bounded amino-acid co-catalysis context for caramelization severity

**What it does not support:**
- Direct HDMF closure
- A universal time-course outside the reported window

**Next Action:** Use as a caramelization process-state calibration rather than as an executable benchmark.

---

### Brands & van Boekel (2002)
- **DOI:** [10.1021/jf015902m](https://doi.org/10.1021/jf015902m)
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
| `mgo_mM` | `50.0` |
| `temperature_c` | `120.0` |
| `pH` | `6.0` |
| `time_minutes` | `120.0` |
| `hdmf_reference_yield_ug_per_g` | `2.9` |
| `replicates` | `3` |

**What it supports:**
- A bounded Family 09 C3-route support anchor showing that HDMF can form directly from methylglyoxal-rich fragmentation chemistry
- Scientist-facing context for alkaline or fragmentation-heavy systems where HDMF support is not purely tied to intact sugar or amino-acid pairing

**What it does not support:**
- A full standalone caramelization solver driven by explicit methylglyoxal mass balance
- Direct PBMA endpoint closure without additional matrix-specific benchmark anchors

**Next Action:** Use as a bounded Family 09 methylglyoxal-to-HDMF support prior so the runtime can surface an explicit C3 route without replacing the default Blank & Fay furanone expectation prior.

---

### Glomb & Monnier (1995)
- **DOI:** [10.1021/jf00058a011](https://doi.org/10.1021/jf00058a011)
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

### Ramaswamy & Zareifard (2000)
- **DOI:** [10.1016/S0260-8774(00)00047-9](https://doi.org/10.1016/S0260-8774(00)00047-9)
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

### Hauck & Tressl (1999)
- **DOI:** [10.1021/bk-1999-0740.ch012](https://doi.org/10.1021/bk-1999-0740.ch012)
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
| `hdmf_from_iso_oligosaccharides_ug_per_g` | `3.8` |
| `hdmf_from_maltose_ug_per_g` | `1.4` |
| `hdmf_from_mgo_self_condensation_ug_per_g` | `2.9` |
| `reaction_temp_C` | `140.0` |

**What it supports:**
- HDMF generation from iso-oligosaccharides (3.8 ug/g) and maltose (1.4 ug/g) at 140 C without amino acids
- Methylglyoxal (MGO) self-condensation HDMF route yield (2.9 ug/g)

**What it does not support:**
- Amadori rearrangement activation barriers
- Peptide-bound lysine trapping kinetics

**Next Action:** Expose as a non-amino acid dependent furanone flavor reference payload.

---

### CHITTENDEN & Regeling (1987)
- **DOI:** [10.1002/recl.19871060202](https://doi.org/10.1002/recl.19871060202)
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

### Luna & Aguilera (2014)
- **DOI:** [10.1016/j.jfoodeng.2014.06.032](https://doi.org/10.1016/j.jfoodeng.2014.06.032)
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

### Bao et al. (2022)
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

Total papers analyzed: **10** (Benchmark-eligible: **0**, Calibration references: **6**, Rejected: **4**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Resconi et al. (2023) | 5/8 | ⚠️ Calibration | pbma_vs_beef_comparator |
| Mottram & Nobrega (2002) | 2/8 | ❌ Rejected | caramelization_sulfur_model_system |
| Ordoudi et al. (2014) | 3/8 | ⚠️ Calibration | caramelization_model_system |
| Brands & van Boekel (2002) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Glomb & Monnier (1995) | 2/8 | ❌ Rejected | aqueous_model_system |
| Ramaswamy & Zareifard (2000) | 2/8 | ❌ Rejected | free_model_system |
| Hauck & Tressl (1999) | 2/8 | ❌ Rejected | free_model_system |
| CHITTENDEN & Regeling (1987) | 3/8 | ⚠️ Calibration | free_model_system |
| Luna & Aguilera (2014) | 3/8 | ⚠️ Calibration | free_model_system |
| Bao et al. (2022) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 9

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Resconi et al. (2023)` (Score: 5/8)
- `Ordoudi et al. (2014)` (Score: 3/8)
- `Brands & van Boekel (2002)` (Score: 3/8)
- `CHITTENDEN & Regeling (1987)` (Score: 3/8)
- `Luna & Aguilera (2014)` (Score: 3/8)
- `Bao et al. (2022)` (Score: 3/8)

### REJECTED (Score < 3)
- `Mottram & Nobrega (2002)` (Score: 2/8)
- `Glomb & Monnier (1995)` (Score: 2/8)
- `Ramaswamy & Zareifard (2000)` (Score: 2/8)
- `Hauck & Tressl (1999)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/recommend.py` must explicitly account for `carbohydrate_pyrolysis_lane` as a modifier when predicting `HMF`, `furfural`, `caramelization_severity` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `failure_mode_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
