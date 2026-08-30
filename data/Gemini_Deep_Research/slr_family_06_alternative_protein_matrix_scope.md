# Systematic Literature Review — Family 6: Alternative protein matrix scope
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering alternative protein matrix scope.  
**Strategic Posture:** `matrix_scope_lane`  
**Runtime Concept:** `matrix_family_scope_extension`  
**Scope & Targets:** Covers target compounds/variables: `matrix_family_support_posture`, `sulfur_deficiency_risk`. Preferred payload types: `computational_prior`, `process_state_calibration`. Target runtime modules: `src/matrix_family_coverage.py`, `src/matrix_prior_registry.py`.  

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

### Asen & Aluko (2022)
- **DOI:** [10.3389/fnut.2022.852225](https://doi.org/10.3389/fnut.2022.852225)
- **Matrix Family:** pea_protein
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
| `base_td_celsius` | `74.45` |
| `heated_fraction_td_celsius_range` | `[124, 206]` |
| `replicates` | `3` |

**What it supports:**
- Pea denaturation-state calibration as a function of pH and heating
- DSC/Td anchor for thermal state heuristics

**What it does not support:**
- Full Ellman + OPA + DSC trifecta in one experiment
- Exact PPI benchmark condition at pH 5.5 and 95C

**Next Action:** Encode as explicit denaturation-state calibration payload for pea process-state modeling.

---

### Li et al. (2025)
- **DOI:** [10.1016/j.crfs.2025.101173](https://doi.org/10.1016/j.crfs.2025.101173)
- **Matrix Family:** pea_protein
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
| `assay` | `DTNB / Ellman` |
| `signal_units` | `nmol per mg protein` |
| `replicates` | `3` |

**What it supports:**
- Free-SH / Ellman anchor for pea cysteine accessibility under heating
- Triplicate process-state calibration input

**What it does not support:**
- OPA-based lysine accessibility in the same experiment
- Exact pH 5.5 PPI condition used by the candidate benchmark

**Next Action:** Encode as explicit cysteine-accessibility calibration payload for pea process-state modeling.

---

### Liu (2023)
- **DOI:** [10.1016/j.foodchem.2022.134998](https://doi.org/10.1016/j.foodchem.2022.134998)
- **Matrix Family:** pea_protein_isolate
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
| `hexanal_ppb_range` | `[15.0, 180.0]` |
| `nonanal_ppb_range` | `[5.0, 50.0]` |
| `ee_2_4_heptadienal_ppb_range` | `[0.5, 8.0]` |
| `methoxypyrazine_ppb_max` | `0.08` |
| `mft_detected` | `False` |
| `meaty_potential_multiplier` | `0.12` |
| `replicate_lot_count` | `6` |

**What it supports:**
- Compact PPI off-note baseline for hexanal, nonanal, unsaturated aldehydes, and methoxypyrazine
- A matrix-scope penalty anchor showing that native PPI starts strongly off-note before any meaty supplementation

**What it does not support:**
- An executable cooked PPI benchmark with matched heat history
- A positive meaty target without exogenous sulfur and pentose support

**Next Action:** Use as the lean PPI matrix-baseline anchor so off-note penalties are explicit before any meaty-positive optimization.

---

### Huseynli et al. (2025)
- **DOI:** [10.3390/foods14111940](https://doi.org/10.3390/foods14111940)
- **Matrix Family:** sunflower_protein_roasted
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
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `roasting_temp_c` | `140.0` |
| `roasting_time_min` | `20.0` |
| `2_methylbutanal_fd_factor` | `2048.0` |
| `3_methylbutanal_fd_factor` | `2048.0` |
| `2_acetylpyrrole_fd_factor` | `512.0` |
| `2_5_dimethylpyrazine_fd_factor` | `256.0` |
| `hexanal_fd_factor` | `128.0` |
| `phenylacetaldehyde_fd_factor` | `256.0` |
| `4_vinylguaiacol_fd_factor` | `512.0` |
| `chlorogenic_acid_dw_pct_range` | `[1.0, 3.0]` |
| `meaty_potential_multiplier` | `0.25` |

**What it supports:**
- A compact roasted-sunflower matrix anchor with strong Strecker aldehydes and explicit chlorogenic-acid-derived phenolic interference
- A Family 06 reference surface for nutty roasted analogues where sunflower is viable but not a direct beef-like proxy

**What it does not support:**
- An executable sunflower meaty benchmark with matched pH and water activity
- A direct lysine-blocking kinetic prior without a dedicated chlorogenic-acid sink payload

**Next Action:** Keep encoded as the compact sunflower matrix anchor so roasted Strecker strength and chlorogenic-acid phenolic interference stay explicit in Family 06 reporting.

---

### Paraskevopoulou et al. (2024)
- **DOI:** [10.3390/foods13081257](https://doi.org/10.3390/foods13081257)
- **Matrix Family:** spirulina_dried_supplement
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
| `alkane_volatiles_pct` | `47.2` |
| `ketone_volatiles_pct` | `25.7` |
| `aldehyde_volatiles_pct` | `10.9` |
| `alcohol_volatiles_pct` | `8.4` |
| `pyrazine_volatiles_pct` | `0.8` |
| `beta_ionone_oav_floor` | `100.0` |
| `one_octen_3_ol_oav` | `10.0` |
| `hexanal_oav_range` | `[2.0, 8.0]` |
| `meaty_potential_multiplier` | `0.04` |

**What it supports:**
- A compact Spirulina matrix anchor where marine and fishy off-notes dominate despite a nominally viable sulfur amino-acid profile
- A Family 06 lower-bound reference for alternative proteins that are chemically misaligned with meat-like flavor targets without lipid extraction

**What it does not support:**
- An executable heated Spirulina meaty benchmark
- A direct sulfur-positive runtime surface without an explicit defatting or extraction correction

**Next Action:** Keep encoded as the lower-bound Spirulina matrix anchor so marine and fishy penalties remain explicit before any alternative-protein recommendation promotes it into meaty search space.

---

### PMC9905368 (2023)
- **DOI:** [10.1007/s10068-022-01194-w](https://doi.org/10.1007/s10068-022-01194-w)
- **Matrix Family:** spi_hvp_xylose
- **Kind:** `quantitative_benchmark`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `temp_C` | `120.0` |
| `time_min` | `30.0` |
| `ph` | `6.0` |
| `mft_ug_per_g` | `0.18` |
| `fft_ug_per_g` | `0.42` |
| `methional_ug_per_g` | `1.88` |
| `mft_oav` | `450` |
| `fft_oav` | `84` |

**What it supports:**
- Executable sulfur-focused SPI-HVP plus xylose benchmark surface at 120 C / 30 min / pH 6.0
- Quantitative MFT OAV 450 and FFT OAV 84 anchor for soy hydrolysate supplementation logic

**What it does not support:**
- Exact hydrolysate peptide-sequence reconstruction
- Whole-panel closure for all 28 reported volatiles in a single mechanistic run

**Next Action:** Keep encoded as a sulfur-focused executable benchmark subset tied to the protein-source registry.

---

### Cho et al. (2023)
- **DOI:** [10.1007/s10068-022-01194-w](https://doi.org/10.1007/s10068-022-01194-w)
- **Matrix Family:** wheat_gluten_hvp_xylose
- **Kind:** `quantitative_benchmark`
- **Payload Role:** `benchmark_intake`
- **Status:** `encoded`

| Criterion | Assessment |
|---|---|
| ✅ | Matrix family / reactant identity specified |
| ✅ | Precursor concentrations / loads specified |
| ✅ | Conditions T, pH, time specified |
| ❌ | Analytical output mode not specified |
| ❌ | Yields reported in relative/qualitative mode only |
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `temp_C` | `120.0` |
| `time_min` | `30.0` |
| `ph` | `6.0` |
| `mft_ug_per_g` | `0.34` |
| `fft_ug_per_g` | `0.61` |
| `methional_ug_per_g` | `3.44` |
| `mft_oav` | `850` |
| `fft_oav` | `122` |

**What it supports:**
- Executable wheat-gluten HVP plus xylose benchmark with the strongest reported non-yeast plant MFT anchor
- Highest-MFT plant-source reference: MFT OAV 850 and FFT OAV 122 at 120 C / 30 min

**What it does not support:**
- Exact gluten peptide-sequence reconstruction
- Direct transfer to non-hydrolyzed wheat gluten formulations without process correction

**Next Action:** Expose as the highest-MFT plant-source sulfur benchmark for protein-source-aware formulation ranking.

---

### Fraser et al. (2018)
- **DOI:** [US9943096B2](https://doi.org/US9943096B2)
- **Matrix Family:** yeast_extract_reaction_flavor
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
| `mft_oav_band` | `[2000.0, 8000.0]` |
| `thiamine_context` | `high_internal_thiamine` |
| `analytical_method` | `sida` |

**What it supports:**
- A high-ceiling sulfur-positive anchor for thiamine-rich yeast-extract systems outside the soy or beef reference bands
- Family 06 matrix-scope support showing that precursor-rich yeast extracts can legitimately target very high MFT OAVs

**What it does not support:**
- A complete yeast-extract benchmark payload in the current runtime surface
- Automatic transfer of the reported OAV band to dilute SPI or PPI systems

**Next Action:** Keep encoded as a sulfur-positive reference anchor until a full executable yeast-extract benchmark package is built.

---

### Liu et al. (2022)
- **DOI:** [10.1016/j.foodchem.2022.134998](https://doi.org/10.1016/j.foodchem.2022.134998)
- **Matrix Family:** pea_protein_isolate
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
| `ee_2_4_heptadienal_fd_factor` | `4096.0` |
| `methoxypyrazine_fd_factor` | `256.0` |

**What it supports:**
- Pea protein heated slurry AEDA and GC-O volatile ranking.
- Identifies (E,E)-2,4-heptadienal (FD 4096) and methoxypyrazine (FD 256) as active off-notes.

**What it does not support:**
- Absolute concentration values for meaty sulfur markers
- Extrusion process-shear texturization parameters

**Next Action:** Expose as a pea protein heating volatile AEDA reference payload.

---

### Malia et al. (2025)
- **DOI:** [10.1016/j.foodhyd.2024.110509](https://doi.org/10.1016/j.foodhyd.2024.110509)
- **Matrix Family:** pea_protein_isolate
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
| `assay` | `DTNB / Ellman` |
| `free_sh_reported` | `True` |

**What it supports:**
- Pea protein free-SH levels and DSC thermal denaturation crosscheck.
- Cysteine accessibility baseline in heated matrices.

**What it does not support:**
- Absolute volatile yield curves
- Lysine-loss kinetics

**Next Action:** Use as a supportive pea free-SH behavior crosscheck.

---

### Lagrain et al. (2010)
- **DOI:** [10.1021/jf102575r](https://doi.org/10.1021/jf102575r)
- **Matrix Family:** wheat_gliadin
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
| `cystine_elimination_ea_ph8_kj_mol` | `119.0` |
| `cystine_elimination_ea_ph6_kj_mol` | `88.2` |
| `lanthionine_formation_ea_ph8_kj_mol` | `12.3` |
| `lanthionine_formation_ea_ph6_kj_mol` | `6.51` |
| `cystine_elimination_k_ph8_min` | `0.054` |
| `replicates` | `3` |

**What it supports:**
- Derived strict first-order kinetics for cystine beta-elimination (Ea 119 kJ/mol at pH 8.0, 88.2 kJ/mol at pH 6.0).
- Second-order kinetics for lanthionine formation (Ea 12.3 kJ/mol at pH 8.0, 6.51 kJ/mol at pH 6.0).

**What it does not support:**
- High-moisture extrusion shear cells

**Next Action:** Encode lanthionine and cystine elimination kinetics priors.

---

### Morel et al. (2002)
- **DOI:** [10.1021/bm015639p](https://doi.org/10.1021/bm015639p)
- **Matrix Family:** wheat_gluten
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
| `shear_reduced_ea_kj_mol` | `33.7` |
| `replicates` | `3` |

**What it supports:**
- Mechanochemical shift reducing apparent Arrhenius activation energy for protein solubility loss under intense mechanical shear from thermally dominated down to 33.7 kJ/mol.

**What it does not support:**
- Free amino acid kinetics in solution

**Next Action:** Encode shear-induced aggregation priors.

---

### Rombouts et al. (2012)
- **DOI:** [10.1021/jf3024672](https://doi.org/10.1021/jf3024672)
- **Matrix Family:** wheat_gluten
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
| `gluten_polymerization_ea_kj_mol` | `142.0` |
| `thermal_threshold_temp_C` | `90.0` |
| `replicates` | `3` |

**What it supports:**
- Overall reaction rate constants for wheat gluten heat-induced polymerization and non-disulfide cross-linking (lanthionine/lysinoalanine) with activation energy of 142.0 kJ/mol above thermal threshold of 90.0 C.

**What it does not support:**
- Low moisture extrusion Apparent Viscosity

**Next Action:** Encode gluten cross-linking priors.

---

### Ilo et al. (1996)
- **DOI:** [10.1006/fstl.1996.0090](https://doi.org/10.1006/fstl.1996.0090)
- **Matrix Family:** maize_grits
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
| `kinetic_order` | `first_order` |
| `replicates` | `3` |

**What it supports:**
- Specific mechanical energy (SME) inputs correlate to chemical breakdown.
- Lysine damage adheres to first-order kinetics and is suppressed by higher moisture acting as a lubricant.

**What it does not support:**
- Free amino acid kinetics in solution

**Next Action:** Encode maize extrusion shear damage priors.

---

### Liu et al. (2022)
- **DOI:** [10.1016/j.foodchem.2022.134998](https://doi.org/10.1016/j.foodchem.2022.134998)
- **Matrix Family:** pea_protein_isolate
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
| `ee_2_4_heptadienal_fd_factor` | `4096.0` |
| `methoxypyrazine_fd_factor` | `256.0` |

**What it supports:**
- Pea protein heated slurry AEDA and GC-O volatile ranking.
- Identifies (E,E)-2,4-heptadienal (FD 4096) and methoxypyrazine (FD 256) as active off-notes.

**What it does not support:**
- Absolute concentration values for meaty sulfur markers
- Extrusion process-shear texturization parameters

**Next Action:** Expose as a pea protein heating volatile AEDA reference payload.

---

### Zhang et al. (1993)
- **DOI:** [10.1111/j.1745-4530.1993.tb00179.x](https://doi.org/10.1111/j.1745-4530.1993.tb00179.x)
- **Matrix Family:** soy_protein
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
| `asn_deamidation_Ea_range_kj_mol` | `[25.1, 113.0]` |
| `gln_deamidation_Ea_range_kj_mol` | `[120.0, 268.0]` |
| `replicates` | `3` |

**What it supports:**
- Deamidation of asparagine (Ea 25.1-113 kJ/mol) vs. glutamine (Ea 120-268 kJ/mol) in protein melts yielding free ammonia.

**What it does not support:**
- Free amino acid solutions without protein structural constraints

**Next Action:** Encode thermal deamidation kinetics.

---

### Li et al. (2010)
- **DOI:** [10.1016/j.mineng.2010.05.003](https://doi.org/10.1016/j.mineng.2010.05.003)
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
| `chelation_Ea_kj_mol` | `84.54` |
| `replicates` | `3` |

**What it supports:**
- Adsorption and coordinate covalent chelation activation energy of phytic acid binding metal ions (Ea 84.54 kJ/mol).

**What it does not support:**
- Divalent cation complexes without phosphate ester coordination

**Next Action:** Encode phytate chelation kinetics prior.

---

### Zha et al. (2020)
- **DOI:** [10.1021/acs.jafc.0c04281](https://doi.org/10.1021/acs.jafc.0c04281)
- **Matrix Family:** plant_protein_matrix
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
- Pea protein isolate glycation aggregation kinetics and lysine blockage parameters.

**What it does not support:**
- Cereal matrix gliadin beta-elimination kinetics.

**Next Action:** Encode PPI glycation prior.

---

### Kutzli et al. (2020)
- **DOI:** [10.1039/D0FO00292E](https://doi.org/10.1039/D0FO00292E)
- **Matrix Family:** plant_protein_matrix
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
- Pseudo-first order glycation within electrospun pea protein-maltodextrin fibers under relative humidity control.

**What it does not support:**
- Free amino acid solution systems without structural fiber constraints.

**Next Action:** Encode fiber glycation prior.

---

### Nguyen et al. (2025)
- **DOI:** [10.1016/j.foodchem.2025.146396](https://doi.org/10.1016/j.foodchem.2025.146396)
- **Matrix Family:** plant_protein_matrix
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
- Pea protein isolate physical barrier coefficients limiting microencapsulated oil self-oxidation.

**What it does not support:**
- Bulk liquid phase systems without microencapsulation boundaries.

**Next Action:** Encode oil stabilization barrier prior.

---

### Pereira et al. (2020)
- **DOI:** [10.3390/antiox9080756](https://doi.org/10.3390/antiox9080756)
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
| `Cu_PM_formation_free_energy_kcal_mol` | `-35.8` |
| `Fe_PM_formation_free_energy_kcal_mol` | `-58.9` |

**What it supports:**
- Transition metal coordination free energy calculations for Cu(II) and Fe(III) complexes (Cu-PM: -35.8 kcal/mol, Fe-PM: -58.9 kcal/mol).

**What it does not support:**
- Biological matrix-level enzyme-mediated chelation paths.

**Next Action:** Encode metal chelation thermodynamic prior.

---

### Sun et al. (2020)
- **DOI:** [10.1016/j.meatsci.2020.108151](https://doi.org/10.1016/j.meatsci.2020.108151)
- **Matrix Family:** high_moisture_melt
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
| `Ea_free_CML_kj_mol` | `44.158` |
| `Ea_free_CEL_kj_mol` | `40.971` |

**What it supports:**
- Zero-order free CML (Ea = 44.158 kJ/mol) and free CEL (Ea = 40.971 kJ/mol) accumulation in a solid matrix during sterilization.

**What it does not support:**
- Extrusion processes with low water activities.

**Next Action:** Encode solid matrix CML/CEL accumulation prior.

---

### Shirai et al. (2015)
- **DOI:** [10.1021/acs.est.5b02902](https://doi.org/10.1021/acs.est.5b02902)
- **Matrix Family:** plant_protein_matrix
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
- Diffusion-limited dityrosine formation modeling using kinetic multilayer surface/bulk approaches.

**What it does not support:**
- Disulfide-based crosslinking mechanisms.

**Next Action:** Encode dityrosine formation prior.

---

### Chen et al. (2024)
- **DOI:** [10.3390/plants13020274](https://doi.org/10.3390/plants13020274)
- **Matrix Family:** hemp_protein
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
| ✅ | Replicates reported |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ✅ | Odor thresholds / sensory reported |

**Score: 4/8 → Calibration reference**

**Key Values:**
| Parameter | Value |
|---|---|
| `e_2_octenal_oav` | `35.0` |
| `replicates` | `3` |

**What it supports:**
- Hemp volatile profile comparing NADES versus alkaline extraction
- Identifies (E)-2-octenal OAV of 35 and beta-caryophyllene as active notes

**What it does not support:**
- Cysteine depletion kinetics

**Next Action:** Expose as a hemp protein off-note flavor reference payload.

---

### Chambers et al. (2018)
- **DOI:** [10.3390/foods7080126](https://doi.org/10.3390/foods7080126)
- **Matrix Family:** soy_isolate
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
| `glutamic_acid_fold_increase` | `2.49` |
| `glycine_fold_increase` | `3.73` |

**What it supports:**
- Aspergillus oryzae soybean fermentation amino acid liberation
- Liberates Glutamic acid (+2.49x) and Glycine (+3.73x) relative to unfermented control

**What it does not support:**
- Hydrolysis shear damage

**Next Action:** Expose as a soybean fermentation process state calibration payload.

---

### Williams (2025)
- **DOI:** [10.1016/j.foodchem.2025.146396](https://doi.org/10.1016/j.foodchem.2025.146396)
- **Matrix Family:** fava_bean
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
| `native_2_pentylfuran_fd` | `8192` |
| `hydrolysed_2_pentylfuran_fd` | `16384` |

**What it supports:**
- Fava bean protein hydrolysate off-note evolution
- 2-pentylfuran FD factor increases from 8192 native up to 16384 post-hydrolysis

**What it does not support:**
- Acrylamide safety limits

**Next Action:** Expose as a fava bean off-note flavor reference payload.

---

### Šimková et al. (2023)
- **DOI:** [10.3390/foods12061349](https://doi.org/10.3390/foods12061349)
- **Matrix Family:** rubisco
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
| `met_content_pct` | `2.2` |
| `cys_content_pct` | `1.9` |
| `gly_content_pct` | `9.5` |

**What it supports:**
- Rubisco storage protein amino acid composition: Met 2.2%, Cys 1.9%, Gly 9.5%
- Exhibits rapid Amadori product formation at low temperatures (80 C)

**What it does not support:**
- Extrusion structuring viscoelastic profiles

**Next Action:** Encode Rubisco amino acid composition and Amadori kinetics priors.

---

### Streule et al. (2024)
- **DOI:** [10.3390/foods13162590](https://doi.org/10.3390/foods13162590)
- **Matrix Family:** lentil_protein
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
| `deflavoring_efficiency_pct` | `78.0` |

**What it supports:**
- Lentil protein isolate volatile profiles and thermal de-flavoring efficiency
- Quantifies residual hexanal OAV post-treatment

**What it does not support:**
- GSH kokumi mouthfulness

**Next Action:** Expose as a lentil isolate off-note flavor reference payload.

---

### PMC11889959 (2025)
- **DOI:** [10.1016/j.foodchem.2025.142222](https://doi.org/10.1016/j.foodchem.2025.142222)
- **Matrix Family:** soy_isolate
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
| `pyrazines_detected` | `True` |

**What it supports:**
- SPI-based textured vegetable protein (TVP) volatile profile
- Quantified pyrazines and 2-pentylfuran baselines

**What it does not support:**
- Thiamine thermal scission split

**Next Action:** Expose as an SPI TVP off-note flavor reference payload.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Alternative protein matrix scope` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `matrix_family_support_posture`, `sulfur_deficiency_risk` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 6

Total papers analyzed: **29** (Benchmark-eligible: **0**, Calibration references: **21**, Rejected: **8**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Asen & Aluko (2022) | 4/8 | ⚠️ Calibration | pea_protein |
| Li et al. (2025) | 4/8 | ⚠️ Calibration | pea_protein |
| Liu (2023) | 3/8 | ⚠️ Calibration | pea_protein_isolate |
| Huseynli et al. (2025) | 2/8 | ❌ Rejected | sunflower_protein_roasted |
| Paraskevopoulou et al. (2024) | 3/8 | ⚠️ Calibration | spirulina_dried_supplement |
| PMC9905368 (2023) | 4/8 | ⚠️ Calibration | spi_hvp_xylose |
| Cho et al. (2023) | 4/8 | ⚠️ Calibration | wheat_gluten_hvp_xylose |
| Fraser et al. (2018) | 3/8 | ⚠️ Calibration | yeast_extract_reaction_flavor |
| Liu et al. (2022) | 2/8 | ❌ Rejected | pea_protein_isolate |
| Malia et al. (2025) | 2/8 | ❌ Rejected | pea_protein_isolate |
| Lagrain et al. (2010) | 3/8 | ⚠️ Calibration | wheat_gliadin |
| Morel et al. (2002) | 3/8 | ⚠️ Calibration | wheat_gluten |
| Rombouts et al. (2012) | 3/8 | ⚠️ Calibration | wheat_gluten |
| Ilo et al. (1996) | 3/8 | ⚠️ Calibration | maize_grits |
| Liu et al. (2022) | 2/8 | ❌ Rejected | pea_protein_isolate |
| Zhang et al. (1993) | 3/8 | ⚠️ Calibration | soy_protein |
| Li et al. (2010) | 3/8 | ⚠️ Calibration | free_model_system |
| Zha et al. (2020) | 3/8 | ⚠️ Calibration | plant_protein_matrix |
| Kutzli et al. (2020) | 3/8 | ⚠️ Calibration | plant_protein_matrix |
| Nguyen et al. (2025) | 3/8 | ⚠️ Calibration | plant_protein_matrix |
| Pereira et al. (2020) | 2/8 | ❌ Rejected | free_model_system |
| Sun et al. (2020) | 3/8 | ⚠️ Calibration | high_moisture_melt |
| Shirai et al. (2015) | 3/8 | ⚠️ Calibration | plant_protein_matrix |
| Chen et al. (2024) | 4/8 | ⚠️ Calibration | hemp_protein |
| Chambers et al. (2018) | 2/8 | ❌ Rejected | soy_isolate |
| Williams (2025) | 3/8 | ⚠️ Calibration | fava_bean |
| Šimková et al. (2023) | 2/8 | ❌ Rejected | rubisco |
| Streule et al. (2024) | 3/8 | ⚠️ Calibration | lentil_protein |
| PMC11889959 (2025) | 2/8 | ❌ Rejected | soy_isolate |

---

## Consolidated entries for benchmark_schema.json — Family 6

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Asen & Aluko (2022)` (Score: 4/8)
- `Li et al. (2025)` (Score: 4/8)
- `Liu (2023)` (Score: 3/8)
- `Paraskevopoulou et al. (2024)` (Score: 3/8)
- `PMC9905368 (2023)` (Score: 4/8)
- `Cho et al. (2023)` (Score: 4/8)
- `Fraser et al. (2018)` (Score: 3/8)
- `Lagrain et al. (2010)` (Score: 3/8)
- `Morel et al. (2002)` (Score: 3/8)
- `Rombouts et al. (2012)` (Score: 3/8)
- `Ilo et al. (1996)` (Score: 3/8)
- `Zhang et al. (1993)` (Score: 3/8)
- `Li et al. (2010)` (Score: 3/8)
- `Zha et al. (2020)` (Score: 3/8)
- `Kutzli et al. (2020)` (Score: 3/8)
- `Nguyen et al. (2025)` (Score: 3/8)
- `Sun et al. (2020)` (Score: 3/8)
- `Shirai et al. (2015)` (Score: 3/8)
- `Chen et al. (2024)` (Score: 4/8)
- `Williams (2025)` (Score: 3/8)
- `Streule et al. (2024)` (Score: 3/8)

### REJECTED (Score < 3)
- `Huseynli et al. (2025)` (Score: 2/8)
- `Liu et al. (2022)` (Score: 2/8)
- `Malia et al. (2025)` (Score: 2/8)
- `Liu et al. (2022)` (Score: 2/8)
- `Pereira et al. (2020)` (Score: 2/8)
- `Chambers et al. (2018)` (Score: 2/8)
- `Šimková et al. (2023)` (Score: 2/8)
- `PMC11889959 (2025)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/matrix_family_coverage.py`, `src/matrix_prior_registry.py` must explicitly account for `matrix_family_scope_extension` as a modifier when predicting `matrix_family_support_posture`, `sulfur_deficiency_risk` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `matrix_scope_lane` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
