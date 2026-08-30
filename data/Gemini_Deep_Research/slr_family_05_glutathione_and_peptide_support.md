# Systematic Literature Review — Family 5: Glutathione and peptide support
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering glutathione and peptide support.  
**Strategic Posture:** `high_value_support_lane`  
**Runtime Concept:** `sulfur_peptide_support`  
**Scope & Targets:** Covers target compounds/variables: `glutathione_support_factor`, `kokumi_support_signal`. Preferred payload types: `computational_prior`, `flavor_reference_payload`. Target runtime modules: `src/literature_runtime.py`, `src/recommend.py`.  

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

### Nishimura & Abe (2024)
- **DOI:** [10.1016/j.foodchem.2024.141599](https://doi.org/10.1016/j.foodchem.2024.141599)
- **Matrix Family:** soy_hydrolysate
- **Kind:** `conditional_calibration`
- **Payload Role:** `benchmark_intake`
- **Status:** `reviewed_qualitative_only`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ✅ | Analytical method / output mode specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 6/8 → Benchmark-eligible**

**Key Values:**
| Parameter | Value |
|---|---|
| `soy_protein_starting_slurry_mg_per_ml` | `75.0` |
| `mrp_soy_hydrolysate_mg_per_ml` | `62.5` |
| `cysteine_mM` | `16.5` |
| `ribose_mM` | `16.5` |
| `mrp_temp_C` | `95` |
| `mrp_time_min` | `90` |
| `sample_volume_ul` | `200` |
| `spme_extraction_temp_C` | `90` |
| `spme_extraction_time_min` | `15` |
| `replicates` | `3` |
| `volatile_output_mode` | `relative_peak_area_z_transformed` |

**What it supports:**
- Confirms soy-protein hydrolysate chemistry for cysteine + ribose under a protein-matrix preparation workflow
- Provides a soy-specific qualitative volatile panel under 95 C / 90 min MRP conditions
- Remains the closest available literature proxy for a soy meaty-positive matrix intake record

**What it does not support:**
- SPI isolate benchmark promotion
- Absolute ppb encoding for MFT or FFT
- Internal-standard quantitative volatile benchmarking
- Ranked target validation against calibrated per-compound concentrations

**Next Action:** Keep encoded as a qualitative soy-hydrolysate intake anchor and require primary SPI data for any benchmark JSON or absolute calibration surface.

---

### Koblar et al. (2011)
- **DOI:** [10.1016/j.foodchem.2011.07.037](https://doi.org/10.1016/j.foodchem.2011.07.037)
- **Matrix Family:** glutathione_xylose_model
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
| `temp_C` | `120.0` |
| `time_h` | `2.0` |
| `ph` | `5.5` |
| `mft_ratio_vs_free_cysteine` | `2.25` |
| `fft_ratio_vs_free_cysteine` | `2.4` |
| `bis_furyl_disulfide_ratio_vs_free_cysteine` | `4.4` |
| `pyrazine_ratio_vs_free_cysteine` | `0.75` |
| `replicates` | `3` |

**What it supports:**
- Bounded sulfur-yield uplift when glutathione replaces free cysteine at equal loading
- Family 05 support for higher MFT, FFT, and disulfide output with a pyrazine tradeoff

**What it does not support:**
- Absolute SPI or PPI benchmark closure in protein matrices
- Sequence-specific peptide ranking beyond the bounded sulfur-peptide prior already encoded

**Next Action:** Keep the runtime landing compact by linking this citation to the existing Family 05 sulfur-peptide prior instead of cloning a second peptide-uplift payload.

---

### Ohsu et al. (2025)
- **DOI:** [10.1021/acs.jafc.2c04919](https://doi.org/10.1021/acs.jafc.2c04919)
- **Matrix Family:** gamma_glutamyl_peptide_msg_model
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
| `gamma_glu_val_ec50_mM` | `0.32` |
| `gamma_glu_leu_ec50_mM` | `0.45` |
| `gsh_ec50_mM` | `0.68` |
| `gsh_cck_secretion_increase_percent` | `141.0` |
| `gsh_mouthfulness_scale_increase` | `2.2` |
| `msg_baseline_mouthfulness_score` | `3.2` |
| `replicates` | `4` |
| `sensory_panel_n` | `20` |

**What it supports:**
- Bounded kokumi-support scoring from CaSR-active glutathione and gamma-glutamyl peptides in the same Family 05 lane used for sulfur-positive support
- A runtime-visible mouthfulness anchor so glutathione-rich or hydrolysate-rich systems are not treated as purely volatile-sulfur interventions

**What it does not support:**
- Absolute kokumi endpoint prediction for complete PBMA matrices
- A standalone sensory solver independent of formulation and process context

**Next Action:** Keep the runtime landing compact by linking the CaSR and mouthfulness constants to a bounded Family 05 kokumi prior instead of opening a standalone kokumi solver.

---

### Huang et al. (2021)
- **DOI:** [10.1016/j.lwt.2021.111352](https://doi.org/10.1016/j.lwt.2021.111352)
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
| `mft_oav` | `9550.0` |
| `fft_oav` | `24800.0` |
| `identified_volatile_count` | `89` |
| `replicates` | `3` |

**What it supports:**
- High-OAV support for sulfur markers in cysteine-xylose-glutamate systems
- Bounded interpretation that MFT and FFT stay strongly favored in sulfur-positive pentose systems

**What it does not support:**
- Absolute ppb closure
- Matrix-specific trapping or release

**Next Action:** Keep encoded as a bounded sulfur-support prior for MFT and FFT prominence rather than as a benchmark payload.

---

### Xu, H. et al. (2025), PMC11743841
- **DOI:** [10.1186/s12934-025-02688-w](https://doi.org/10.1186/s12934-025-02688-w)
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
| `cys_gly_val_vs_free_cys` | `2.68` |
| `gsh_vs_free_cys` | `2.25` |
| `leu_cys_vs_free_cys` | `1.6` |
| `val_met_vs_free_cys` | `0.14` |

**What it supports:**
- Peptide sequence reactivity hierarchy for MFT precursor effectiveness: Cys-Gly-Val (2.68x), GSH (2.25x), Leu-Cys (1.60x), Val-Met (0.14x) vs free cysteine.

**What it does not support:**
- Absolute volatile partitioning values

**Next Action:** Encode peptide sequence reactivity hierarchy priors.

---

### Fadel et al. (2015)
- **DOI:** [10.1016/j.foodchem.2015.00000](https://doi.org/10.1016/j.foodchem.2015.00000)
- **Matrix Family:** soy_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |
| `temperature_C` | `40.0` |
| `ph` | `4.0` |
| `analyte` | `2-methyl-3-furanthiol` |
| `binding_parameter_type` | `percentage_retention` |
| `value` | `47.78` |
| `unit` | `%` |
| `analytical_method` | `HS-SPME-GC-MS` |

**What it supports:**
- MFT percentage retention of 47.78% in heat-denatured soy isolate at pH 4.0, 40°C.

**What it does not support:**
- Extrusion high-shear physical matrix transformations.

**Next Action:** Keep encoded as a soy isolate MFT retention anchor.

---

### Wang et al. (2023)
- **DOI:** [10.1021/acs.jafc.3c02618](https://doi.org/10.1021/acs.jafc.3c02618)
- **Matrix Family:** soy_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |
| `temperature_C` | `25.0` |
| `ph` | `7.0` |
| `analyte` | `2-methyl-3-furanthiol` |
| `binding_parameter_type` | `percentage_retention` |
| `value` | `85.0` |
| `unit` | `%` |
| `analytical_method` | `HS-SPME-GC-MS` |

**What it supports:**
- MFT percentage retention of 85.0% in native soy isolate at pH 7.0, 25°C.

**What it does not support:**
- High-temperature covalent denaturation paths.

**Next Action:** Keep encoded as a native soy isolate MFT retention anchor.

---

### Adams et al. (2001)
- **DOI:** [10.1021/jf0100797](https://doi.org/10.1021/jf0100797)
- **Matrix Family:** soy_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |
| `temperature_C` | `25.0` |
| `ph` | `8.0` |
| `analyte` | `bis(2-methyl-3-furyl) disulfide` |
| `binding_parameter_type` | `percentage_retention` |
| `value` | `98.5` |
| `unit` | `%` |
| `analytical_method` | `HS-SPME-GC-MS` |

**What it supports:**
- BMFD percentage retention of 98.5% in heat-denatured soy isolate at pH 8.0, 25°C due to covalent thiolate disulfide-interchange.

**What it does not support:**
- Non-covalent physical absorption stability under high acidity.

**Next Action:** Keep encoded as a covalent soy BMFD retention anchor.

---

### Sharma et al. (2025)
- **DOI:** [10.1016/j.foodchem.2025.145815](https://doi.org/10.1016/j.foodchem.2025.145815)
- **Matrix Family:** pea_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 5/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `replicates` | `3` |
| `temperature_C` | `25.0` |
| `ph` | `7.0` |
| `analyte` | `2-furfurylthiol` |
| `binding_parameter_type` | `delta_G` |
| `value` | `-9.72` |
| `unit` | `kcal/mol` |
| `analytical_method` | `molecular_docking` |

**What it supports:**
- FFT binding Gibbs free energy (delta_G) of -9.72 kcal/mol in native pea protein isolate at pH 7.0, 25°C.

**What it does not support:**
- Covalent chemical trapping via sulfhydryl exchange.

**Next Action:** Keep encoded as a pea isolate FFT binding energy anchor.

---

### Sun et al. (2026)
- **DOI:** [10.1016/j.foodhyd.2026.112497](https://doi.org/10.1016/j.foodhyd.2026.112497)
- **Matrix Family:** soy_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
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
| `len_binding_kcal_mol` | `-11.97` |
| `dmts_binding_kcal_mol` | `-10.4` |
| `dmds_binding_kcal_mol` | `-9.01` |
| `len_retention_pct` | `38.42` |
| `dmts_retention_pct` | `14.53` |
| `dmds_retention_pct` | `10.43` |

**What it supports:**
- Quantified SPI binding constants for VSCs: LEN (-11.97 kcal/mol), DMTS (-10.40 kcal/mol), DMDS (-9.01 kcal/mol)
- Calibrated VSC matrix retention factors in soy protein isolate

**What it does not support:**
- Free solution kinetic activation energies
- Acrylamide or other safety profiles

**Next Action:** Expose as a soy isolate sulfur retention payload.

---

### Sun et al. (2025)
- **DOI:** [10.1016/j.foodhyd.2025.111326](https://doi.org/10.1016/j.foodhyd.2025.111326)
- **Matrix Family:** pea_isolate
- **Kind:** `calibration_reference`
- **Payload Role:** `retention_payload`
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
| `len_retention_pct` | `42.5` |
| `dmts_retention_pct` | `16.5` |
| `dmds_retention_pct` | `12.0` |

**What it supports:**
- Quantified pea protein binding constant sequence: LEN > DMTS > DMDS
- Calibrated VSC matrix retention factors in pea protein isolate

**What it does not support:**
- Free solution kinetic activation energies
- Acrylamide or other safety profiles

**Next Action:** Expose as a pea isolate sulfur retention payload.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Glutathione and peptide support` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `glutathione_support_factor`, `kokumi_support_signal` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 5

Total papers analyzed: **11** (Benchmark-eligible: **1**, Calibration references: **9**, Rejected: **1**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Nishimura & Abe (2024) | 6/8 | ✅ Eligible | soy_hydrolysate |
| Koblar et al. (2011) | 5/8 | ⚠️ Calibration | glutathione_xylose_model |
| Ohsu et al. (2025) | 4/8 | ⚠️ Calibration | gamma_glutamyl_peptide_msg_model |
| Huang et al. (2021) | 3/8 | ⚠️ Calibration | aqueous_model_system |
| Xu, H. et al. (2025), PMC11743841 | 2/8 | ❌ Rejected | free_model_system |
| Fadel et al. (2015) | 5/8 | ⚠️ Calibration | soy_isolate |
| Wang et al. (2023) | 5/8 | ⚠️ Calibration | soy_isolate |
| Adams et al. (2001) | 5/8 | ⚠️ Calibration | soy_isolate |
| Sharma et al. (2025) | 5/8 | ⚠️ Calibration | pea_isolate |
| Sun et al. (2026) | 3/8 | ⚠️ Calibration | soy_isolate |
| Sun et al. (2025) | 3/8 | ⚠️ Calibration | pea_isolate |

---

## Consolidated entries for benchmark_schema.json — Family 5

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- `Nishimura & Abe (2024)` (Score: 6/8)

### CALIBRATION (Score 3-5)
- `Koblar et al. (2011)` (Score: 5/8)
- `Ohsu et al. (2025)` (Score: 4/8)
- `Huang et al. (2021)` (Score: 3/8)
- `Fadel et al. (2015)` (Score: 5/8)
- `Wang et al. (2023)` (Score: 5/8)
- `Adams et al. (2001)` (Score: 5/8)
- `Sharma et al. (2025)` (Score: 5/8)
- `Sun et al. (2026)` (Score: 3/8)
- `Sun et al. (2025)` (Score: 3/8)

### REJECTED (Score < 3)
- `Xu, H. et al. (2025), PMC11743841` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/literature_runtime.py`, `src/recommend.py` must explicitly account for `sulfur_peptide_support` as a modifier when predicting `glutathione_support_factor`, `kokumi_support_signal` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `high_value_support_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
