# Systematic Literature Review — Family 15: PE Stealth Sugar Sink
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering pe stealth sugar sink.  
**Strategic Posture:** `upstream_precursor_sink`  
**Runtime Concept:** `sugar_depletion_sink`  
**Scope & Targets:** Covers target compounds/variables: `PE_glycation_fraction`, `available_sugar_subtraction`. Preferred payload types: `upstream_modifier_payload`, `computational_prior`. Target runtime modules: `src/precursor_resolver.py`, `src/flavor.py`.  

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

### Solís-Calero et al. (2015)
- **DOI:** [10.1371/journal.pone.0124658](https://doi.org/10.1371/journal.pone.0124658)
- **Matrix Family:** phospholipid_interface
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
| `pe_schiff_base_ea_kj_mol` | `92.9` |
| `pe_amadori_ea_kj_mol` | `82.9` |
| `protein_lys_maillard_ea_kj_mol` | `118.0` |
| `interfacial_pseudocatalysis_reported` | `True` |

**What it supports:**
- Interfacial PE Schiff-base and Amadori Arrhenius anchors already landed in the runtime DFT-prior surface
- Scientist-facing provenance that phospholipid-rich systems can accelerate amine Maillard initiation below the protein-only barrier

**What it does not support:**
- One-condition phospholipid endpoint benchmark closure
- Full mixed-lipid membrane kinetics

**Next Action:** Expose the already-landed Family 15 Arrhenius anchors through intake so phospholipid acceleration stops looking unused in the literature tracker.

---

### Kodate et al. (2018)
- **DOI:** [10.1038/s41598-018-27010-2](https://doi.org/10.1038/s41598-018-27010-2)
- **Matrix Family:** aminophospholipid_food_matrix
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
| `amadori_pe_nmol_per_ml_range` | `[10.0, 100.0]` |
| `pe_glycation_fraction_ceiling_pct` | `30.0` |
| `lac_pe_recovery_pct_range` | `[73.0, 87.8]` |
| `glc_pe_positive_ion_mz` | `326.1` |
| `lac_pe_positive_ion_mz` | `488.1` |

**What it supports:**
- A food-matrix burden anchor for Amadori-PE accumulation that complements the interfacial PE Arrhenius priors
- A bounded Family 15 sugar-sequestration ceiling before dedicated soy-lecithin HMEC data exists

**What it does not support:**
- A direct HMEC lecithin endpoint benchmark in SPI or PPI extrudates
- A full oxidation or ALE progression model in the same run

**Next Action:** Keep encoded as the first Family 15 food-matrix burden anchor so the PE sugar sink has magnitude support in addition to mechanistic DFT-derived priors.

---

### Fujimoto et al. (2023)
- **DOI:** [10.1021/acs.jafc.2c08283](https://doi.org/10.1021/acs.jafc.2c08283)
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
| `sugar_pe_molar_stoichiometry` | `0.5` |
| `uv_absorbance_nm` | `350.0` |

**What it supports:**
- 1:2 sugar to phospholipid (PE) stoichiometry and pyridinium derivative UV absorption at 350 nm
- Calcium stearate inhibition mechanism parameters in 160 C oil models

**What it does not support:**
- Alternative headspace retention coefficients

**Next Action:** Encode PE glycation stoichiometry priors.

---

### Solís-Calero et al. (2015)
- **DOI:** [10.1039/C4CP05360E](https://doi.org/10.1039/C4CP05360E)
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
| `Ea_dehydration_kcal_mol` | `17.5` |
| `Ea_dehydration_kj_mol` | `73.2` |
| `lysine_dehydration_Ea_diff_kcal_mol` | `-5.3` |

**What it supports:**
- Kinetic dehydration barrier reduction (by 5.3 kcal/mol / 22.2 kJ/mol vs free L-lysine) on organized PE surface during CM-PE formation.

**What it does not support:**
- High-moisture extrusion structures without lipid bilayer phase organization.

**Next Action:** Encode phospholipid surface dehydration catalyst prior.

---

### Solís-Calero et al. (2013)
- **DOI:** [10.1021/jp401488j](https://doi.org/10.1021/jp401488j)
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
| `Ea_condensation_kcal_mol` | `8.76` |
| `Ea_condensation_kj_mol` | `36.7` |
| `Ea_enaminol_kcal_mol` | `16.78` |
| `Ea_enaminol_kj_mol` | `70.2` |

**What it supports:**
- Catalytic reduction of the initial sugar-PE condensation barrier to 8.76 kcal/mol (36.7 kJ/mol), with subsequent 1,2-enaminol formation at 16.78 kcal/mol (70.2 kJ/mol).

**What it does not support:**
- Bulk aqueous phase glycation without interfacial alignment.

**Next Action:** Encode PE surface sugar condensation prior.

---

### Lertsiri et al. (1998)
- **DOI:** [10.1271/bbb.62.893](https://doi.org/10.1271/bbb.62.893)
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
| `glycation_rate_mol_pct_min` | `0.05` |
| `glycation_rate_mol_pct_max` | `0.12` |

**What it supports:**
- Pseudo-first-order accumulation of Amadori-PE under excess glucose, with physiological accumulation rates of 0.05 to 0.12 mol%.

**What it does not support:**
- Dicarbonyl-mediated direct advanced glycation pathways.

**Next Action:** Encode Amadori-PE accumulation prior.

---

### Hidalgo et al. (2005)
- **DOI:** [10.1007/s00217-004-1065-x](https://doi.org/10.1007/s00217-004-1065-x)
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
| `Ea_browning_kj_mol` | `66.5` |
| `Ea_fluorescence_kj_mol` | `50.0` |

**What it supports:**
- Global activation energy for non-enzymatic browning (66.5 kJ/mol) and advanced fluorescence development (50.0 kJ/mol) in lipid/amine systems.

**What it does not support:**
- Metal-catalyzed lipid oxidation without amino acid presence.

**Next Action:** Encode PE-ribose-lysine global browning prior.

---

### Zamora et al. (2020)
- **DOI:** [10.3390/molecules25020373](https://doi.org/10.3390/molecules25020373)
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
| `dpph_scavenging_increase_pct` | `85.1` |

**What it supports:**
- Secondary radical scavenging capacity increase of 85.1% via 2-pentyl-3,5-dibutyl-dihydropyridine adduct formation in soybean oil emulsions.

**What it does not support:**
- Protein-bound amine crosstalk without phospholipid interface partition.

**Next Action:** Encode dihydropyridine adduct radical scavenging prior.

---

### Hidalgo et al. (2006)
- **DOI:** [10.1021/jf060848s](https://doi.org/10.1021/jf060848s)
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
| `induction_period_synergy_pct` | `185.0` |

**What it supports:**
- Synergistic 185% increase in oxidative induction period utilizing 300 ppm PE and 100 ppm lysine, illustrating polar paradox kinetics.

**What it does not support:**
- Purely aqueous phase radical scavenging without lipid interfaces.

**Next Action:** Encode PE-lysine Rancimat oxidative synergy prior.

---

### Vilanova et al. (2012)
- **DOI:** [10.1021/jp2116033](https://doi.org/10.1021/jp2116033)
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
| `Ea_dehydration_kcal_mol` | `13.08` |
| `Ea_dehydration_kj_mol` | `54.7` |

**What it supports:**
- Second-order rates and dehydration barrier evaluation of 13.08 kcal/mol (54.7 kJ/mol) for Schiff base formation.

**What it does not support:**
- Long-chain polysaccharide steric hindrance effects.

**Next Action:** Encode Schiff base dehydration rate prior.

---

### Biondi et al. (2010)
- **DOI:** [10.1016/j.lwt.2010.02.016](https://doi.org/10.1016/j.lwt.2010.02.016)
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
- Microwave thermal-stress degradation rates of vegetable oils correlated with volatile hexanal accumulation and diene formation.

**What it does not support:**
- Non-thermal oxidation induction phase modeling.

**Next Action:** Encode microwave thermal stress lipid oxidation prior.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `PE Stealth Sugar Sink` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `PE_glycation_fraction`, `available_sugar_subtraction` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 15

Total papers analyzed: **11** (Benchmark-eligible: **0**, Calibration references: **10**, Rejected: **1**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| Solís-Calero et al. (2015) | 3/8 | ⚠️ Calibration | phospholipid_interface |
| Kodate et al. (2018) | 3/8 | ⚠️ Calibration | aminophospholipid_food_matrix |
| Fujimoto et al. (2023) | 2/8 | ❌ Rejected | free_model_system |
| Solís-Calero et al. (2015) | 3/8 | ⚠️ Calibration | free_model_system |
| Solís-Calero et al. (2013) | 3/8 | ⚠️ Calibration | free_model_system |
| Lertsiri et al. (1998) | 3/8 | ⚠️ Calibration | free_model_system |
| Hidalgo et al. (2005) | 3/8 | ⚠️ Calibration | free_model_system |
| Zamora et al. (2020) | 3/8 | ⚠️ Calibration | free_model_system |
| Hidalgo et al. (2006) | 3/8 | ⚠️ Calibration | free_model_system |
| Vilanova et al. (2012) | 3/8 | ⚠️ Calibration | free_model_system |
| Biondi et al. (2010) | 3/8 | ⚠️ Calibration | free_model_system |

---

## Consolidated entries for benchmark_schema.json — Family 15

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Solís-Calero et al. (2015)` (Score: 3/8)
- `Kodate et al. (2018)` (Score: 3/8)
- `Solís-Calero et al. (2015)` (Score: 3/8)
- `Solís-Calero et al. (2013)` (Score: 3/8)
- `Lertsiri et al. (1998)` (Score: 3/8)
- `Hidalgo et al. (2005)` (Score: 3/8)
- `Zamora et al. (2020)` (Score: 3/8)
- `Hidalgo et al. (2006)` (Score: 3/8)
- `Vilanova et al. (2012)` (Score: 3/8)
- `Biondi et al. (2010)` (Score: 3/8)

### REJECTED (Score < 3)
- `Fujimoto et al. (2023)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/precursor_resolver.py`, `src/flavor.py` must explicitly account for `sugar_depletion_sink` as a modifier when predicting `PE_glycation_fraction`, `available_sugar_subtraction` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `upstream_precursor_sink` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
