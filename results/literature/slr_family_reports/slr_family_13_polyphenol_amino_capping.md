# Systematic Literature Review — Family 13: Polyphenol-Amino Capping
**Last updated:** 2026-05-29  
**Objective:** Identify papers with quantitative data usable as benchmarks or calibration references covering polyphenol-amino capping.  
**Strategic Posture:** `upstream_precursor_sink`  
**Runtime Concept:** `precursor_depletion_sink`  
**Scope & Targets:** Covers target compounds/variables: `quinone_budget`, `cysteine_depletion_factor`. Preferred payload types: `upstream_modifier_payload`, `computational_prior`. Target runtime modules: `src/precursor_resolver.py`, `src/literature_runtime.py`.  

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

### J. Agric. Food Chem. 2019 (Ref. 24)
- **DOI:** [10.1021/acs.jafc.9b07752](https://doi.org/10.1021/acs.jafc.9b07752)
- **Matrix Family:** polyphenol_reactive_model_system
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
| `cys_4mc_adduct_nmol_per_mg_range` | `[2.2, 8.1]` |
| `thiol_loss_fraction` | `0.62` |
| `digestive_stability_umol_per_g` | `1.7` |

**What it supports:**
- Quantified cysteine-adduct burden in the polyphenol capping lane
- Bounded free-thiol depletion prior for sulfur-support interpretation

**What it does not support:**
- Endpoint volatile closure
- Universal transfer across all polyphenol classes

**Next Action:** Use as a bounded polyphenol-thiol capping prior in runtime lanes that need explicit cysteine sink behavior.

---

### Mesías et al. (2020)
- **DOI:** [10.3390/foods9081026](https://doi.org/10.3390/foods9081026)
- **Matrix Family:** amino_acid_capped_model_system
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
| ❌ | Replicates < 3 or not specified |
| ⚠️ | LOD/LOQ not explicitly confirmed (qualitative/provisional status) |
| ❌ | No sensory or odor threshold tags present |

**Score: 2/8 → Rejected**

**Key Values:**
| Parameter | Value |
|---|---|
| `acrylamide_suppression_fraction` | `0.978` |
| `mitigation_pair` | `['cysteine', 'glycine']` |

**What it supports:**
- An upper-bound Family 13 mitigation anchor for cysteine-plus-glycine acrylamide suppression
- Safety-side interpretation of aggressive amino-acid capping strategies

**What it does not support:**
- Generalized kinetic closure across arbitrary baking or extrusion profiles
- Direct finished-product prediction without benchmark context

**Next Action:** Keep encoded as the strongest Family 13 safety mitigation anchor without misclassifying it as a universal kinetic rule.

---

### Hosry et al. (2025)
- **DOI:** [10.3390/foods14111881](https://doi.org/10.3390/foods14111881)
- **Matrix Family:** coffee_melanoidin_model
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
| `fft_reduction_factor` | `16.0` |
| `volatile_adduct_range_ug_per_l` | `[5.0, 270.0]` |
| `crosspy_mediated` | `True` |

**What it supports:**
- FFT 16x reduction via CROSSPY-radical mediated covalent binding to melanoidins under warm storage conditions
- Quantified QA/CQA/QAL-volatile adduct bounds ranging from 5 to 270 ug/L

**What it does not support:**
- Direct pea-protein or soy-protein matrix volatile benchmarking
- Reversible non-covalent binding kinetics

**Next Action:** Keep as the Family 13 covalent thiol-depletion anchor for CROSSPY-mediated capping modeling.

---

### Poojary et al. (2023)
- **DOI:** [10.1016/j.foodchem.2022.134406](https://doi.org/10.1016/j.foodchem.2022.134406)
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
| `cga_cys_mrm_mz` | `474.1` |
| `stability_temp_C` | `90.0` |

**What it supports:**
- CGA-Cys and CA-Cys adduct SIDA quantification (m/z 474.1 MRM) showing stability at 90 C

**What it does not support:**
- Alternative sensory descriptor mapping

**Next Action:** Encode CGA-Cys adduct SIDA stability priors.

---

### Shu et al. (2019)
- **DOI:** [10.1016/j.freeradbiomed.2019.04.026](https://doi.org/10.1016/j.freeradbiomed.2019.04.026)
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
| `quinone_thiol_k_range` | `[100000.0, 1000000.0]` |
| `replicates` | `3` |

**What it supports:**
- Second-order rate constants of 10^5 to 10^6 M-1 s-1 for quinone-thiol Michael additions representing rapid cysteine depletion by oxidized polyphenols.

**What it does not support:**
- Non-oxidizing environments without electron-acceptor or metal catalysts

**Next Action:** Encode cysteine-quinone Michael adduct kinetics.

---

### Zhu et al. (2020)
- **DOI:** [10.1021/acs.jafc.0c01761](https://doi.org/10.1021/acs.jafc.0c01761)
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
| `replicates` | `4` |
| `MGO_trapping_k_M1_s1` | `1.6` |
| `GO_trapping_k_M1_s1` | `0.059` |

**What it supports:**
- Second-order epicatechin methylglyoxal (k = 1.6 M-1 s-1) and glyoxal (k = 0.059 M-1 s-1) trapping in stored milk models.

**What it does not support:**
- Enzymatic polyphenol oxidase-catalyzed capping.

**Next Action:** Encode epicatechin dicarbonyl trapping prior.

---

### Zhu et al. (2020)
- **DOI:** [10.1016/j.foodchem.2020.126500](https://doi.org/10.1016/j.foodchem.2020.126500)
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
| `kaempferol_MGO_k_M1_s1` | `0.063` |
| `resveratrol_MGO_k_M1_s1` | `0.027` |

**What it supports:**
- Structure-dependent MGO trapping kinetics for kaempferol (k = 0.063 M-1 s-1) and resveratrol (k = 0.027 M-1 s-1).

**What it does not support:**
- Catechin-class flavanol trapping mechanisms.

**Next Action:** Encode structural polyphenol trapping prior.

---

### Liu et al. (2022)
- **DOI:** [10.1016/j.foodres.2022.112187](https://doi.org/10.1016/j.foodres.2022.112187)
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
| `Ea_4MBQ_Lysine_kj_mol` | `19.0` |

**What it supports:**
- Second-order 4-methylbenzoquinone-lysine conjugation kinetics (Ea = 19.00 kJ/mol).

**What it does not support:**
- Michael additions to soft thiol sulfur nucleophiles.

**Next Action:** Encode benzoquinone-lysine conjugation prior.

---

### Cömert & Gökmen (2019)
- **DOI:** [10.1016/j.foodres.2019.03.046](https://doi.org/10.1016/j.foodres.2019.03.046)
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
| `egg_MGO_rate_k_L_mol_min` | `26.6` |

**What it supports:**
- Methylglyoxal scavenging kinetics under simulated gastrointestinal digestion fluid (egg-MGO rate k = 26.6 L/mol/min).

**What it does not support:**
- High-temperature dry extrusion caramelization.

**Next Action:** Encode digestive dicarbonyl scavenging prior.

---

### Munoz et al. (2007)
- **DOI:** [10.1021/jf062081+](https://doi.org/10.1021/jf062081+)
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
| `CQA_Q_reaction_rate_k_M1_s1` | `2.73` |

**What it supports:**
- Chlorogenic acid enzymatically oxidized to o-quinone reacting with substrate (k = 2.73 M-1 s-1).

**What it does not support:**
- Non-catechol structural polyphenol trapping paths.

**Next Action:** Encode chlorogenic o-quinone prior.

---

### Cömert & Gökmen (2021)
- **DOI:** [10.1016/j.foodchem.2020.128670](https://doi.org/10.1016/j.foodchem.2020.128670)
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
- Synergistic rate augmentation between epicatechin and cysteine in trapping methylglyoxal.

**What it does not support:**
- Antagonistic pathways driven by gallic acid variants.

**Next Action:** Encode EC-Cys synergistic scavenging prior.

---

### Song et al. (2009)
- **DOI:** [10.1073/pnas.0810352106](https://doi.org/10.1073/pnas.0810352106)
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
| `forward_rate_k_M1_s1` | `1547.0` |
| `reverse_Ea_kcal_mol` | `11.2` |

**What it supports:**
- Benzoquinone conjugation with GSH rate parameters (forward k = 1547 M-1 s-1, reverse Ea = 11.2 kcal/mol).

**What it does not support:**
- Reactions catalyzed by polyphenol oxidases.

**Next Action:** Encode quinone-GSH conjugation prior.

---

### Li et al. (2017)
- **DOI:** [10.1021/acs.jafc.6b05811](https://doi.org/10.1021/acs.jafc.6b05811)
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
- Quercetin-MGO dicarbonyl trapping kinetics mitigating browning in lysine/glucose Maillard models.

**What it does not support:**
- Purely non-enzymatic lipid oxidation paths.

**Next Action:** Encode quercetin MGO adduction prior.

---

### Monforte et al. (2018)
- **DOI:** [10.1021/acs.jafc.7b00264](https://doi.org/10.1021/acs.jafc.7b00264)
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
| `Ea_hydroquinone_catalyzed_kj_mol` | `32.9` |
| `Ea_trihydroxybenzene_catalyzed_kj_mol` | `31.5` |

**What it supports:**
- Phenylacetaldehyde formation activation energy catalyzed by hydroquinone (Ea = 32.9 kJ/mol) or trihydroxybenzene (Ea = 31.5 kJ/mol).

**What it does not support:**
- Strecker degradation in non-polyphenol systems.

**Next Action:** Encode polyphenol-catalyzed Strecker prior.

---

### Rawel et al. (2005)
- **DOI:** [10.1021/jf0480290](https://doi.org/10.1021/jf0480290)
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
| `phenolic_compounds` | `['chlorogenic acid', 'ferulic acid', 'gallic acid', 'quercetin', 'rutin']` |
| `protein_targets` | `['soy glycinin', 'human serum albumin', 'bovine serum albumin', 'lysozyme']` |

**What it supports:**
- Binding constants for polyphenols (chlorogenic, ferulic, gallic acids) to soy glycinin
- Quantifies non-covalent protein-phenolic interactions as a function of pH and temperature

**What it does not support:**
- Volatile sulfur compounds binding
- Reaction pathways kinetics

**Next Action:** Expose as a soy protein phenolic binding reference.

---

## Confirmed Gaps and Structural Limitations

- **Matrix Transferability Gap:** Most kinetic and volatile data for `Polyphenol-Amino Capping` are derived from simple aqueous buffers or free model systems. Transferability to complex plant-based meat matrices (e.g., Pea Protein Isolate or Soy Protein Isolate) under high-moisture extrusion remains a major gap.
- **Quantitative Yield Gap:** There is a scarcity of absolute concentrations (e.g., in ng/g or ppb) for the target analytes `quinone_budget`, `cysteine_depletion_factor` under realistic cooking profiles.
- **Analytical Standardisation Gap:** Many literature sources rely on relative peak area percentages rather than stable isotope dilution assays (SIDA) or absolute external calibrations, limiting their direct ingestion as hard ODE benchmarks.

---

## EXECUTIVE SUMMARY: COVERAGE MAP FOR FAMILY 13

Total papers analyzed: **15** (Benchmark-eligible: **0**, Calibration references: **12**, Rejected: **3**)

| Reference | Score | Status | Focus |
|---|---|---|---|
| J. Agric. Food Chem. 2019 (Ref. 24) | 2/8 | ❌ Rejected | polyphenol_reactive_model_system |
| Mesías et al. (2020) | 2/8 | ❌ Rejected | amino_acid_capped_model_system |
| Hosry et al. (2025) | 3/8 | ⚠️ Calibration | coffee_melanoidin_model |
| Poojary et al. (2023) | 2/8 | ❌ Rejected | free_model_system |
| Shu et al. (2019) | 3/8 | ⚠️ Calibration | free_model_system |
| Zhu et al. (2020) | 3/8 | ⚠️ Calibration | free_model_system |
| Zhu et al. (2020) | 3/8 | ⚠️ Calibration | free_model_system |
| Liu et al. (2022) | 3/8 | ⚠️ Calibration | free_model_system |
| Cömert & Gökmen (2019) | 3/8 | ⚠️ Calibration | free_model_system |
| Munoz et al. (2007) | 3/8 | ⚠️ Calibration | free_model_system |
| Cömert & Gökmen (2021) | 3/8 | ⚠️ Calibration | free_model_system |
| Song et al. (2009) | 3/8 | ⚠️ Calibration | free_model_system |
| Li et al. (2017) | 3/8 | ⚠️ Calibration | free_model_system |
| Monforte et al. (2018) | 3/8 | ⚠️ Calibration | free_model_system |
| Rawel et al. (2005) | 3/8 | ⚠️ Calibration | soy_isolate |

---

## Consolidated entries for benchmark_schema.json — Family 13

### REFERENCE_ANCHOR / PRIMARY (Score >= 6)
- *No primary benchmark-eligible references found.*

### CALIBRATION (Score 3-5)
- `Hosry et al. (2025)` (Score: 3/8)
- `Shu et al. (2019)` (Score: 3/8)
- `Zhu et al. (2020)` (Score: 3/8)
- `Zhu et al. (2020)` (Score: 3/8)
- `Liu et al. (2022)` (Score: 3/8)
- `Cömert & Gökmen (2019)` (Score: 3/8)
- `Munoz et al. (2007)` (Score: 3/8)
- `Cömert & Gökmen (2021)` (Score: 3/8)
- `Song et al. (2009)` (Score: 3/8)
- `Li et al. (2017)` (Score: 3/8)
- `Monforte et al. (2018)` (Score: 3/8)
- `Rawel et al. (2005)` (Score: 3/8)

### REJECTED (Score < 3)
- `J. Agric. Food Chem. 2019 (Ref. 24)` (Score: 2/8)
- `Mesías et al. (2020)` (Score: 2/8)
- `Poojary et al. (2023)` (Score: 2/8)

---

## Model corrections identified during review

1. **Precursor Sourcing & Translation:** The precursor resolution logic in `src/precursor_resolver.py`, `src/literature_runtime.py` must explicitly account for `precursor_depletion_sink` as a modifier when predicting `quinone_budget`, `cysteine_depletion_factor` yields.
2. **Dynamic Calibration Offsets:** Apply the identified `upstream_precursor_sink` parameters to set the relative boundaries and calibration offsets for the chemical species.
3. **Environmental Sensitivities:** Implement temperature and pH scaling factors to simulate matrices where buffering or thermal degradation deviates from optimal model system pathways.

*Programmatically generated on 2026-05-29 using Maillard Ingestion Pipeline.*
