# Systematic Literature Review — Family 16: Melanoidin Polymerization
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering melanoidin polymerization.  
**Strategic Posture:** `trapping_burden_modifier`  
**Runtime Concept:** `matrix_trapping_node`  
**Scope & Targets:** Covers target compounds/variables: `melanoidin_mass`, `TH_scavenging_factor`. Preferred payload types: `trapping_burden_modifier_payload`, `computational_prior`. Target runtime modules: `src/extrusion_model.py`, `src/flavor.py`.  

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

### Weerawatanakorn et al. (2015)
- **DOI:** [10.1021/acs.jafc.5b02009](https://doi.org/10.1021/acs.jafc.5b02009)
- **Matrix Family:** melanoidin_trapping_model
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
| `mft_rfd_baseline` | `4.0` |
| `mft_rfd_with_melanoidins` | `2.0` |
| `melanoidin_loading_reference_mg` | `125.0` |
| `mechanism` | `crosspy_radical_thioether_adduction` |

**What it supports:**
- Irreversible CROSSPY-driven MFT trapping support for high-melanoidin process states
- Family 16 justification for non-recoverable thiol loss instead of purely reversible sulfur attenuation

**What it does not support:**
- A full texture-to-flavor co-model for WHC, WSI, and G-prime collapse
- Species-resolved closure for every thiol in every extrusion condition

**Next Action:** Keep encoded as a bounded Family 16 trapping prior until the extrusion severity surface and same-run sulfur retention benchmarks are wider.

---

### Snel et al. (2023)
- **DOI:** [10.3390/foods12132543](https://doi.org/10.3390/foods12132543)
- **Matrix Family:** spi_ppi_hme_rework
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
| `spi_whc_change_pct_by_cycle` | `[-17.0, -63.0]` |
| `ppi_whc_change_pct_by_cycle` | `[-65.0, -90.0]` |
| `g_prime_recovery_after_rework` | `False` |

**What it supports:**
- A direct HME rework anchor for supra-additive WHC collapse in SPI and PPI
- Family 16 reporting of melanoidin-rich extrusion damage beyond simple browning proxies

**What it does not support:**
- A coupled rheology-and-volatile closure model
- Same-run sulfur retention quantification

**Next Action:** Keep encoded as the primary HME rework hydration-collapse anchor until same-run flavor-plus-texture datasets exist.

---

### Cruz et al. (2025)
- **DOI:** [10.1021/acsfoodscitech.4c00677](https://doi.org/10.1021/acsfoodscitech.4c00677)
- **Matrix Family:** spi_rice_hme_burger_analogues
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
| `breaking_force_n_range` | `[10.0, 104.0]` |
| `optimal_moisture_pct_range` | `[20.0, 25.0]` |

**What it supports:**
- A bounded HME firmness envelope for SPI-rich analogues across the practical moisture window
- Family 16 process-state visibility for texture optimization under melanoidin-heavy extrusion

**What it does not support:**
- A universal firmness model across all protein sources
- An Ogden-based constitutive solver in the current runtime

**Next Action:** Keep encoded as the first direct HME firmness calibration anchor while the extrusion-mechanics solver remains bounded and report-first.

---

### J. Agric. Food Chem. 2019 (Ref. 21)
- **DOI:** [10.1021/acs.jafc.9b04099](https://doi.org/10.1021/acs.jafc.9b04099)
- **Matrix Family:** pea_protein_hydrolysate_gum_arabic_conjugate
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
| `native_mw_kda` | `290.0` |
| `conjugated_mw_kda_range` | `[340.0, 359.0]` |
| `native_radius_of_gyration_nm` | `29.5` |
| `conjugated_radius_of_gyration_nm_range` | `[30.5, 33.5]` |
| `native_mark_houwink_exponent` | `0.53` |
| `conjugated_mark_houwink_exponent` | `0.4` |
| `native_pdi` | `1.35` |
| `conjugated_pdi` | `1.74` |

**What it supports:**
- A compact Family 16 architecture anchor for Maillard-conjugated pea hydrolysate systems
- Runtime visibility that gum arabic plus hydrolysate conditions can shift polymer size, gyration radius, and compaction toward higher trapping pressure

**What it does not support:**
- A universal hydrocolloid rule across all PBMA formulations
- A full texture-and-headspace closure model in the current runtime

**Next Action:** Use as a compact Family 16 architecture anchor so gum-arabic-mediated polymer growth becomes runtime-visible without opening a broad rheology solver.

---

### Gigl et al. (2021)
- **DOI:** [10.1021/acs.jafc.1c06163](https://doi.org/10.1021/acs.jafc.1c06163)
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
| `fft_depletion_fold` | `16.0` |
| `binding_ic50_mg_thiol_per_L_melanoidin` | `183.0` |

**What it supports:**
- Thiol (FFT) depletion by melanoidins (16-fold reduction)
- IC50 values (183 mg thiol / L melanoidin) and binding saturation parameters

**What it does not support:**
- Lysine accessibility changes

**Next Action:** Encode melanoidin thiol trapping priors.

---

### Brands et al. (2002)
- **DOI:** [10.1021/jf010789c](https://doi.org/10.1021/jf010789c)
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
| `rate_k_h1` | `1.14` |
| `Ea_kj_mol` | `128.0` |
| `extinction_coefficient_glucose_casein` | `1200.0` |

**What it supports:**
- Zero-order terminal polymerization kinetics (k = 1.14 h-1, Ea = 128 kJ/mol) and molar extinction coefficients for casein-sugar systems.

**What it does not support:**
- Low-temperature non-polymeric browning pathways.

**Next Action:** Encode casein-sugar melanoidin polymerization prior.

---

### Gigl et al. (2021)
- **DOI:** [10.1021/acs.jafc.1c06163](https://doi.org/10.1021/acs.jafc.1c06163)
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
| `thiol_loss_20min_pct` | `50.0` |
| `line_broadening_hz` | `2.14` |

**What it supports:**
- Covalent trapping of volatile thiols (FWHM line broadening of 2.14 Hz, ~50% loss in 20 min) by coffee melanoidins.

**What it does not support:**
- Reversible non-covalent aroma-matrix stacking equilibria.

**Next Action:** Encode melanoidin covalent thiol staling prior.

---

### Suzuki & Philp (1990)
- **DOI:** [10.1016/0146-6380(90)90162-s](https://doi.org/10.1016/0146-6380(90)90162-s)
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
- Zero-order accumulation of sulfur-incorporated high-molecular-weight melanoidin precipitates in the presence of hydrogen sulfide.

**What it does not support:**
- Short-chain flavor heterocyclic generation.

**Next Action:** Encode sulfur-incorporated melanoidin prior.

---

### Mundt & Wedzicha (2007)
- **DOI:** [10.1016/j.lwt.2006.07.014](https://doi.org/10.1016/j.lwt.2006.07.014)
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
| `Ea_browning_kj_mol` | `105.0` |

**What it supports:**
- First-order surface browning rate (Ea = 105 kJ/mol) independent of water activity (Aw 0.04-0.4) during dough baking.

**What it does not support:**
- Liquid solution phase color development.

**Next Action:** Encode biscuit surface browning rate prior.

---

### Cao et al. (2024)
- **DOI:** [10.1111/1750-3841.17378](https://doi.org/10.1111/1750-3841.17378)
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
| `secondary_metabolite_inhibition_pct` | `68.19` |
| `aggregation_reduction_pct` | `36.95` |

**What it supports:**
- Melanoidins (1% MRPs) inhibit lipid oxidation secondary metabolites by 68.19% and protect myoglobin stability by reducing aggregation by 36.95%.

**What it does not support:**
- Protein-free vegetable oil systems.

**Next Action:** Encode myoglobin-MRP lipid oxidation protection prior.

---

### Hofmann et al. (2001)
- **DOI:** [10.1021/jf001302l](https://doi.org/10.1021/jf001302l)
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
| `thiol_loss_30min_pct` | `90.0` |

**What it supports:**
- Irreversible covalent thioether formation of volatile thiols (FWHM >90% loss in 30 min) directly into coffee melanoidins.

**What it does not support:**
- Reversible hydrogen-bonding flavor entrapment.

**Next Action:** Encode irreversible melanoidin-thioether staling prior.

---

### Adams et al. (2005)
- **DOI:** [10.1021/jf047903m](https://doi.org/10.1021/jf047903m)
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
| `Ea_histidine_kj_mol` | `35.31` |
| `Ea_lysine_kj_mol` | `54.94` |

**What it supports:**
- Browning product (A420) formation rate activation energies for ascorbic acid systems with basic amino acids (Ea = 35.31 kJ/mol for histidine, 54.94 kJ/mol for lysine).

**What it does not support:**
- Sugar-amino acid Maillard core pathways.

**Next Action:** Encode ascorbic-basic-amino browning prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Melanoidin Polymerization` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `melanoidin_mass`, `TH_scavenging_factor` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 16

Total papers analyzed: **12** (Benchmark-eligible: **0**, Calibration references: **7**, Rejected: **5**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Weerawatanakorn et al. (2015) | 2/8 | ❌ Rejected | melanoidin_trapping_model |
| Snel et al. (2023) | 2/8 | ❌ Rejected | spi_ppi_hme_rework |
| Cruz et al. (2025) | 2/8 | ❌ Rejected | spi_rice_hme_burger_analogues |
| J. Agric. Food Chem. 2019 (Ref. 21) | 2/8 | ❌ Rejected | pea_protein_hydrolysate_gum_arabic_conjugate |
| Gigl et al. (2021) | 2/8 | ❌ Rejected | free_model_system |
| Brands et al. (2002) | 3/8 | ⚠️ Calibration | free_model_system |
| Gigl et al. (2021) | 3/8 | ⚠️ Calibration | free_model_system |
| Suzuki & Philp (1990) | 3/8 | ⚠️ Calibration | free_model_system |
| Mundt & Wedzicha (2007) | 3/8 | ⚠️ Calibration | free_model_system |
| Cao et al. (2024) | 3/8 | ⚠️ Calibration | free_model_system |
| Hofmann et al. (2001) | 3/8 | ⚠️ Calibration | free_model_system |
| Adams et al. (2005) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 16

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Brands et al. (2002)` (Score: 3/8)
- `Gigl et al. (2021)` (Score: 3/8)
- `Suzuki & Philp (1990)` (Score: 3/8)
- `Mundt & Wedzicha (2007)` (Score: 3/8)
- `Cao et al. (2024)` (Score: 3/8)
- `Hofmann et al. (2001)` (Score: 3/8)
- `Adams et al. (2005)` (Score: 3/8)

### REJECTED (Score < 3)
- `Weerawatanakorn et al. (2015)` (Score: 2/8)
- `Snel et al. (2023)` (Score: 2/8)
- `Cruz et al. (2025)` (Score: 2/8)
- `J. Agric. Food Chem. 2019 (Ref. 21)` (Score: 2/8)
- `Gigl et al. (2021)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/extrusion_model.py`, `src/flavor.py` must explicitly account for `matrix_trapping_node` as a modifier when predicting `melanoidin_mass`, `TH_scavenging_factor` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `trapping_burden_modifier` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
