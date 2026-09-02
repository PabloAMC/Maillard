# Maillard Simulation Report - Ribose + cysteine + leucine

**Date:** 2026-09-01 22:11:12

## 1. Input Formulation & Conditions
| Parameter | Value |
| :--- | :--- |
| name | Ribose + cysteine + leucine |
| protein_type | pea_iso |
| process_state | aqueous_pre_extrusion_model |
| target | meaty |
| minimize | beany |
| temp | 105.0 |
| ph | 5.5 |
| time_minutes | 45.0 |
| aw | 0.85 |
| generated_by | scripts/generators/generate_report_visual_examples.py |
| note | REAL pipeline output, not a rendering fixture. Regenerate with scripts/generators/generate_report_visual_examples.py. |

## 2. Decision Summary
```text

════════════════════════════════════════════════════════════════════════════════
                         📊 MAILLARD DECISION SUMMARY
════════════════════════════════════════════════════════════════════════════════

  🟡 DIRECTIONAL ONLY
      Confidence is low. Use for hypothesis generation, not decision-grade prioritization.
  ────────────────────────────────────────────────────────────────────────────

  [0] DECISION CONFIDENCE: LOW (47/100)
      Benchmark Basis  : matrix_intake_only
      Decision Mode    : directional_hypothesis
      Prediction Mode : directional_only
      Accessibility    : free_like ✅
      Recommended Use  : Use directionally only; absolute concentrations should be treated as provisional.
      Why              : Plant-matrix support is still intake/headspace validated rather than target-ranking validated.
      Why              : Average pathway uncertainty is moderate (3.5 kcal/mol).
      Calibration      : Recommendation extrapolates beyond the strongest support on: benchmark_neighborhood, matrix, precursors.
      Extrapolation    : benchmark_neighborhood, matrix, precursors

  [1] SCIENTIFIC ENVELOPE: LIMITED ⚠️
    [!!] MATRIX         : Matrix 'pea_iso' uses speculative accessibility scaling; PRIMARY benchmarks are free-precursor only.
    [!] PRECURSORS     : Sparse benchmark analogies: leucine lack PRIMARY quantitative validation.

  [2] TOP-LINE METRICS:
      Target Score     : 0.27
      Off-Flavour Risk : 0.01
      Safety Score     : 1.04  (2x band-relative risk; higher is worse; >1.0 = above the action band)

════════════════════════════════════════════════════════════════════════════════

```

## 3. Detailed Results
- **Target Score:** 0.27
- **Off-Flavour Risk:** 0.01
- **Safety Score:** 1.04 *(2× band-relative risk, range [0, 2]; **higher is worse**; >1.0 = above the action band)*

- **MFT/Furfural Ratio:** 0.0023
- **Meaty Quality Penalty:** 0.79

- **Strecker Balance Score:** 0.00
- **Strecker Gap Penalty:** 0.00
- **Pyrazine Propensity:** 0.65
- **Pyrazine Burden:** 0.00
- **Pyrazine Penalty:** 0.00
- **Furanone Penalty:** 0.35

### Confidence & Support
- **tier:** low
- **score:** 47.0
- **benchmark_neighborhood:** matrix_intake_only
- **decision_mode:** directional_hypothesis
- **prediction_mode:** directional_only
- **recommended_posture:** Use directionally only; absolute concentrations should be treated as provisional.
- **factor:** Plant-matrix support is still intake/headspace validated rather than target-ranking validated.
- **factor:** Average pathway uncertainty is moderate (3.5 kcal/mol).
- **factor:** Domain-of-validity warnings indicate extrapolation beyond the most trusted scientific envelope.

### Calibration Diagnostics
- **supported_envelope:** False
- **summary:** Recommendation extrapolates beyond the strongest support on: benchmark_neighborhood, matrix, precursors.
- **extrapolation_axes:** benchmark_neighborhood, matrix, precursors

### Benchmark Neighborhood
- **benchmark_neighborhood:** matrix_intake_only
- **category:** matrix_transfer
- **prediction_mode:** directional_only
- **summary:** Run uses matrix intake/headspace support and transferred accessibility priors rather than direct ranking benchmarks.

### Calibration Summary
| Source | Support Origin | Evidence | Fallback | Compounds | Observable ppb |
| :--- | :--- | :--- | :--- | :--- | ---: |
| Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) | standard_matrix_support | class_anchored | class_level | furfural | 0.11 |
| Interpolated base sulfur yield matching internal benchmark limits | standard_matrix_support | directional_transferred | class_level | 2-furfurylthiol, 2-methyl-3-furanthiol, bis(2-methyl-3-furyl) disulfide | 0.00 |
| Pratap-Singh 2021 pea isolate ambient slurry baseline (generic aldehyde transfer) | standard_matrix_support | class_anchored | class_level | 3-methylbutanal | 0.00 |
| Interpolated base pyrazine yield matching internal benchmark limits | standard_matrix_support | directional_transferred | class_level | 2,5-dimethylpyrazine | 0.00 |

### Extrusion Process Model
- **model:** sequential_isothermal_zones
- **moisture_regime:** hme
- **jacket_temperature_celsius:** 105.0
- **effective_temperature_celsius:** 105.0
- **die_exit_temperature_celsius:** 70.0
- **sme_kj_per_kg:** 0.0
- **mean_residence_time_seconds:** 44.1
- **rtd_decision:** sequential_zone_sufficient_for_current_use_case
- **furosine_mg_per_kg:** 18.0
- **lal_mg_per_kg:** 68.0

| Extrusion volatile transport | Class | Shear-release | Die stripping | Combined factor |
| :--- | :--- | ---: | ---: | ---: |
| 2-methyl-3-furanthiol | sulfur | 1.00 | 0.08 | 0.92 |
| 2-furfurylthiol | sulfur | 1.00 | 0.08 | 0.92 |
| Hexanal | aldehyde | 1.00 | 0.18 | 0.82 |
| 2-Pentylfuran | furan | 1.00 | 0.09 | 0.91 |
| 2,5-Dimethylpyrazine | pyrazine | 1.00 | 0.01 | 0.99 |

### Compound Confidence
| Compound | Predicted | Tier | Score | Mode | Reachability | Calibration Source | Observable Assumption |
| :--- | :--- | :---: | ---: | :--- | :--- | :--- | :--- |
| furfural | 0.108 ppb [0.000973-11.9, 90% CI] | exploratory | 43.0 | hypothesis_only | chemically_reachable | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) | static_class_profile \| class_level \| standard_matrix_support |
| 2-furfurylthiol | 0.0021 ppb [4.6e-06-0.232, 90% CI] | exploratory | 33.0 | hypothesis_only | chemically_reachable | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior \| class_level \| standard_matrix_support |
| 3-methylbutanal | 0.000903 ppb | exploratory | 43.0 | hypothesis_only | chemically_reachable | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic aldehyde transfer) | static_class_profile \| class_level \| standard_matrix_support |
| 2-methyl-3-furanthiol | 0.000252 ppb [2.28e-06-0.0278, 90% CI] | exploratory | 33.0 | hypothesis_only | chemically_reachable | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior \| class_level \| standard_matrix_support |
| bis(2-methyl-3-furyl) disulfide | 5.51e-05 ppb [5.26e-08-0.00609, 90% CI] | exploratory | 33.0 | hypothesis_only | conditionally_reachable | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior \| class_level \| standard_matrix_support |

### Compound Confidence Overlay
- **figure:** compound_confidence_overlay.png
- **evidence classes:** calibration = green, external = amber, surrogate = slate.

### Intervention Waterfall
- **baseline:** Ribose + leucine (no cysteine)
- **current:** Ribose + cysteine + leucine
- **Furans:** lowered 0.07 ppb (-41%)
- **Sulfur / Thiols:** raised 0.00 ppb
- **Aldehydes:** lowered 0.00 ppb (-41%)
- **figure:** intervention_waterfall.png

### Compound Evidence Ladder
| Compound | Class | Evidence State | Direct Anchor | Transferred Prior | Mechanistic Surrogate | Computational Refinement | Support Origin | Source |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- |
| furfural | severity_markers | conditional_calibration | - | yes | - | - | standard_matrix_support | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) |
| 2-furfurylthiol | sulfur_positives | externally_benchmarked | yes | - | yes | - | standard_matrix_support | Interpolated base sulfur yield matching internal benchmark limits |
| 3-methylbutanal | strecker_aldehydes | externally_benchmarked | yes | yes | - | - | standard_matrix_support | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic aldehyde transfer) |
| 2-methyl-3-furanthiol | sulfur_positives | externally_benchmarked | yes | - | yes | - | standard_matrix_support | Interpolated base sulfur yield matching internal benchmark limits |
| bis(2-methyl-3-furyl) disulfide | sulfur_positives | transferred_prior | - | yes | yes | - | standard_matrix_support | Interpolated base sulfur yield matching internal benchmark limits |
| 2,5-dimethylpyrazine | pyrazines | externally_benchmarked | yes | - | - | - | standard_matrix_support | Interpolated base pyrazine yield matching internal benchmark limits |

### Missing Data
- HEMF: mechanistically expected but still unobserved in the current prediction surface.
- DMHF: mechanistically expected but still unobserved in the current prediction surface.
- Benchmark neighborhood remains extrapolative relative to the primary free-precursor validation envelope.

### Safety Reference Context
| Visibility | Kind | Source | Summary |
| :--- | :--- | :--- | :--- |
| default | industrial_endpoint_reference | Squeo et al. (2023), Foods 12(6):1331 | Upper-bound reference context for acrylamide in processed plant-protein ingredients |
| default | finished_product_reference | Fu et al. (2023) | Finished-product PBMA acrylamide reporting context |
| default | precursor_correlation_reference | Jung et al. (2024), Food Science and Biotechnology 33:2333 | Precursor-correlation explanation for pea safety outputs |
| default | industrial_endpoint_reference | Fu et al. (2023) | Acrylamide concentration range of 31.81-186.70 ug/kg under heating above 150 C |
| default | kinetic_model_reference | De Vleeschouwer, Van der Plancken & Van Loey (2006), JAFC 54:7847 | Pseudo-first-order formation/elimination kinetics model |
| default | kinetic_model_reference | Ma et al. (2022), Frontiers in Nutrition 9:940202 | Demonstration of lipid oxidation carbonyl crosstalk in accelerating acrylamide formation |
| default | industrial_endpoint_reference | Ma et al. (2024) | CML and CEL selective generation under extrusion barrel temperatures of 160 C. |
| default | industrial_endpoint_reference | Knol et al. (2005) | Bifurcated acrylamide formation (Ea 94.4 +/- 11 kJ/mol) and degradation (Ea 85.1 +/- 14 kJ/mol) kinetic parameters. |
| default | industrial_endpoint_reference | De Vleeschouwer et al. (2008) | Moisture-dependent acrylamide degradation with minimum elimination at Aw 0.82. |

- Extended safety provenance entries available in JSON: 7

- 6 safety reference entries were EXCLUDED from this acrylamide context: they carry no `analyte` field, and their legacy `category` value does not resolve to the requested analyte. Their evidence is not represented above; see excluded_entries.

### Flavor Reference Policy
| Compound | Pipeline Role | Benchmark Role | Source |
| :--- | :--- | :--- | :--- |
| 2,3-butanedione | diagnostic_marker | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| 2-methylbutanal | diagnostic_marker | reference_anchor | Huseynli et al. (2025) |
| 3-methylbutanal | diagnostic_marker | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| 3-methylbutanal | diagnostic_marker | reference_anchor | Huseynli et al. (2025) |
| acetoin | diagnostic_marker | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| methional | diagnostic_marker | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| phenylacetaldehyde | diagnostic_marker | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| (E,E)-2,4-heptadienal | no_verifiable_source | retired_unverifiable | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 |
| 3-isobutyl-2-methoxypyrazine | no_verifiable_source | retired_unverifiable | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 |
| 2-ethyl-3,5-dimethylpyrazine | optimization_constraint | directional_comparison_anchor | Arsa & Puechkamutr (2022), Journal of Food Science and Technology 59:890 |
| DMHF | optimization_constraint | low_confidence_mechanistic_anchor | Watanabe et al. (2015), Meat Science |
| HDMF | optimization_constraint | reference_anchor | Blank & Grosch (1991) |
| HEMF | optimization_constraint | low_confidence_mechanistic_anchor | Blank & Fay (1996), JAFC 44:531 |
| furfural | optimization_constraint | reference_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| hexanal | optimization_constraint | reference_anchor | Tao et al. (2022), Front. Microbiol. 13:1070773 |
| methylpyrazine | optimization_constraint | reference_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| nonanal | optimization_constraint | reference_anchor | Marquez-Ruiz et al. (2014), JAFC 62:10295 |
| pyrazine_family | optimization_constraint | directional_comparison_anchor | Laemont & Barringer (2023), Foods 12(22):4155 |
| pyrazine_family | optimization_constraint | low_confidence_mechanistic_anchor | Wang et al. (2021), Foods 10:273 |
| 2-furfurylthiol | primary_target | reference_anchor | Kerscher & Grosch (1998), Journal of Agricultural and Food Chemistry 46:1954 |
| 2-methyl-3-furanthiol | primary_target | reference_anchor | Kerscher & Grosch (1998), Journal of Agricultural and Food Chemistry 46:1954 |
| 2-methyl-3-furanthiol | primary_target | reference_anchor | US 9,943,096 B2 (2018), Impossible Foods Inc. |
| (E)-2-octenal | reference_only | reference_anchor | Chen, Oliveira, Dias & Ismail (2025), Plants 14(2):274 |
| (E,E)-2,4-heptadienal | reference_only | reference_anchor | Liu et al. (2023), Food Chemistry 406:134998 |
| (E,E)-2,4-heptadienal | reference_only | reference_anchor | Liu et al. (2023), Food Chemistry 406:134998 |
| 1-hexanol | reference_only | reference_anchor | Li et al. (2026), Foods 15(5):912 |
| 1-octen-3-ol | reference_only | reference_anchor | Paraskevopoulou et al. (2024) |
| 2-furfurylthiol | reference_only | reference_anchor | Siripitakpong, Wongprasert, Rungrotmongkol & Suppavorasatit (2026), Food Chem. X 34:103712 |
| 2-methyl-3-furanthiol | reference_only | pbma_counterexample | Hernandez et al. (2023), Molecules 28:3151 |
| 2-methyl-3-furanthiol | reference_only | reference_anchor | Tang, Jiang, Yuan & Ho (2013), J. Sulfur Chem. 34:38-47 |
| 2-methyl-3-furanthiol | reference_only | reference_anchor | Wang et al. (2023) |
| 2-methyl-3-furanthiol | reference_only | reference_anchor | Fadel et al. (2015) |
| 2-pentylfuran | reference_only | reference_anchor | Li et al. (2026), Foods 15(5):912 |
| 2-pentylfuran | reference_only | reference_anchor | Williams, A. D. (2025), MS thesis, Virginia Tech (VTechWorks hdl:10919/137837) |
| 2-pentylfuran | reference_only | reference_anchor | Park et al. (2025), Current Research in Food Science 10:100999 |
| 4-vinylguaiacol | reference_only | reference_anchor | Huseynli et al. (2025) |
| HDMF | reference_only | reference_anchor | Hauck & Tressl (1999) |
| benzaldehyde | reference_only | directional_comparison_anchor | Hernandez et al. (2023), Molecules 28:3151 |
| beta-ionone | reference_only | reference_anchor | Paraskevopoulou et al. (2024) |
| bis(2-methyl-3-furyl)disulfide | reference_only | reference_anchor | Adams et al. (2001) |
| chlorogenic acid | reference_only | reference_anchor | Rawel, Czajka, Rohn & Kroll (2002), Int. J. Biol. Macromol. 30:137 |
| heptanal | reference_only | reference_anchor | Li et al. (2026), Foods 15(5):912 |
| hexanal | reference_only | reference_anchor | Bi et al. (2020), J. Agric. Food Chem. 68:2718-2727 |
| hexanal | reference_only | reference_anchor | Bi et al. (2020), J. Agric. Food Chem. 68:2718-2727 |
| hexanal | reference_only | reference_anchor | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 |
| hexanal | reference_only | reference_anchor | Paraskevopoulou et al. (2024) |
| hexanal | reference_only | reference_anchor | Li et al. (2026), Foods 15(5):912 |
| hexanal | reference_only | reference_anchor | Grosch, W. (1982), in Food Flavours Part A (Morton & MacLeod, eds.), Elsevier, pp. 325-398 |
| hexanal | reference_only | reference_anchor | Xu et al. (2024) |
| hexanal | reference_only | reference_anchor | Vurro, De Angelis, Squeo, Caponio, Summo & Pasqualone (2024), Foods 13:2608 |
| hexanal | reference_only | reference_anchor | Liu et al. (2023), Food Chemistry 406:134998 |
| hexanal | reference_only | reference_anchor | Liu et al. (2023), Food Chemistry 406:134998 |
| hexanal | reference_only | reference_anchor | Yeo & Shibamoto (1991) |
| methional | reference_only | reference_anchor | Mottram (1998), Food Chemistry |
| methoxypyrazines | reference_only | reference_anchor | Liu et al. (2023), Food Chemistry 406:134998 |
| nonanal | reference_only | reference_anchor | Liu, Y. (2021), "Flavor Chemistry of Pea Proteins", MS thesis, North Carolina State University; published as Liu, Cadwallader & Drake (2023), Food Chemistry 406:134998 |
| nonanal | reference_only | reference_anchor | Li et al. (2026), Foods 15(5):912 |
| nonanal | reference_only | reference_anchor | Marquez-Ruiz et al. (2014) |
| 2-methylbutanal | secondary_marker | reference_anchor | Hernandez et al. (2023), Molecules 28:3151 |

### Literature Evidence Summary
- **source:** data/lit/benchmark_intake_registry.json
- **eligible_reference_count:** 205
- **ready_reference_count:** 52
- **closable_without_primary_data_count:** 205
- **structural_gap_count:** 5
- **ready_reference_ids:** wang_2012_gsh_xylose_sulfur_uplift, ohsu_2025_kokumi_casr_anchor, squeo_2023, asen_2022, li_2025, de_leyn_2019, comunian_2021_thiamine_encapsulation, voelker_2021_thiamine_kinetics
- **structural_gap_ids:** ppi_meaty_positive_matrix_benchmark, spi_meaty_positive_matrix_benchmark, mft_fft_matrix_retention, ppi_spi_time_series, meaty_off_flavour_safety_tradeoff_panel

### Literature Learning Loop Summary
- **ready_reference_count:** 197
- **encoded_runtime_reference_count:** 197
- **template_queue_count:** 0
- **matrix_family_count:** 83
- **intake_structural_gap_count:** 5
- **process_gap_count:** 5
- **matrix_prior_families:** pea_iso, soy_iso, myco

### Family Runtime Support Summary
| SLR | Family | Active | Posture | Evidence Posture | Primary Payloads | Supporting Payloads | Priors |
| :--- | :--- | :---: | :--- | :--- | ---: | ---: | ---: |
| 01 | Amino acid-sugar Maillard core | yes | first_class_core | core_benchmarked_chemistry | 40 | 77 | 1 |
| 02 | Lipid oxidation and carbonylic crosstalk | - | immediate_expansion_lane | structural_gap_extrapolation | 66 | 59 | 34 |
| 03 | Thiamine degradation and sulfur support | - | high_value_support_lane | structural_gap_extrapolation | 20 | 51 | 8 |
| 04 | Nucleotide degradation and ribose support | yes | high_value_support_lane | structural_gap_extrapolation | 25 | 20 | 10 |
| 05 | Glutathione and peptide support | - | high_value_support_lane | structural_gap_extrapolation | 26 | 48 | 3 |
| 06 | Alternative protein matrix scope | - | matrix_scope_lane | structural_gap_extrapolation | 46 | 33 | 11 |
| 07 | Reducing sugar and carbonyl donor hierarchy | yes | immediate_expansion_lane | structural_gap_extrapolation | 29 | 58 | 17 |
| 08 | Plant off-notes and Maillard suppression | - | guardrail_lane | structural_gap_extrapolation | 8 | 107 | 1 |
| 09 | Carbohydrate pyrolysis and caramelization | yes | failure_mode_lane | directional_priors | 31 | 31 | 7 |
| 10 | Microbial fermentation pretreatment | - | upstream_pretreatment_lane | structural_gap_extrapolation | 5 | 85 | 0 |
| 11 | Maillard/Lipid Crosstalk | - | first_class_runtime_lane | structural_gap_extrapolation | 53 | 47 | 8 |
| 12 | Protein Damage Markers | yes | first_class_runtime_lane | structural_gap_extrapolation | 28 | 55 | 1 |
| 13 | Polyphenol-Amino Capping | - | upstream_precursor_sink | structural_gap_extrapolation | 34 | 30 | 17 |
| 14 | Ascorbic Acid Maillard | - | bounded_upstream_source | structural_gap_extrapolation | 31 | 26 | 14 |
| 15 | PE Stealth Sugar Sink | - | upstream_precursor_sink | structural_gap_extrapolation | 21 | 17 | 8 |
| 16 | Melanoidin Polymerization | yes | trapping_burden_modifier | structural_gap_extrapolation | 22 | 21 | 10 |

### Family Evidence Ladder
| SLR | Family | Evidence Posture | Direct Anchors | Transferred Priors | Surrogates | Compounds |
| :--- | :--- | :--- | ---: | ---: | ---: | :--- |
| 01 | Amino acid-sugar Maillard core | core_benchmarked_chemistry | 4 | 2 | 3 | 2,5-dimethylpyrazine, 2-furfurylthiol, 2-methyl-3-furanthiol, 3-methylbutanal, bis(2-methyl-3-furyl) disulfide |
| 09 | Carbohydrate pyrolysis and caramelization | directional_priors | 0 | 1 | 0 | furfural |

### Family Specific Open Gaps
- **Nucleotide degradation and ribose support:** active_runtime_lane_without_direct_benchmark_or_calibration_grade_payload_closure
- **Reducing sugar and carbonyl donor hierarchy:** active_runtime_lane_without_direct_benchmark_or_calibration_grade_payload_closure
- **Protein Damage Markers:** active_runtime_lane_without_direct_benchmark_or_calibration_grade_payload_closure
- **Melanoidin Polymerization:** active_runtime_lane_without_direct_benchmark_or_calibration_grade_payload_closure
- **Amino acid plus sugar core Maillard chemistry:** expected_marker_missing:HEMF
- **Amino acid plus sugar core Maillard chemistry:** expected_marker_missing:DMHF

### Family Lane Sensitivity
| SLR | Family Lane | Active | Target Δ | Closure Δ | Off-flavour Δ | Toggle Magnitude |
| :--- | :--- | :---: | ---: | ---: | ---: | ---: |
| 07 | Reducing sugar and carbonyl donor hierarchy | yes | +0.11 | +0.00 | -0.01 | 0.13 |
| 01 | Amino acid-sugar Maillard core | yes | +0.06 | +0.00 | +0.00 | 0.06 |
| 04 | Nucleotide degradation and ribose support | yes | +0.06 | +0.00 | +0.00 | 0.06 |
| 09 | Carbohydrate pyrolysis and caramelization | yes | +0.06 | +0.00 | +0.00 | 0.06 |
| 16 | Melanoidin Polymerization | yes | -0.02 | -0.00 | +0.02 | 0.04 |
| 12 | Protein Damage Markers | yes | -0.01 | -0.00 | +0.01 | 0.02 |
| 02 | Lipid oxidation and carbonylic crosstalk | - | +0.00 | +0.00 | +0.00 | 0.00 |
| 03 | Thiamine degradation and sulfur support | - | +0.00 | +0.00 | +0.00 | 0.00 |
| 05 | Glutathione and peptide support | - | +0.00 | +0.00 | +0.00 | 0.00 |
| 06 | Alternative protein matrix scope | - | +0.00 | +0.00 | +0.00 | 0.00 |

### Sensitivity Summary
- **mode:** local_oat
- **evaluated_perturbations:** 9
- **ranking_driver:** ribose | remove_sugars | Δdecision -1.28 | Δsafety +0.00
- **ranking_driver:** temp | decrease_temp | Δdecision +1.08 | Δsafety -1.04
- **ranking_driver:** temp | increase_temp | Δdecision -0.09 | Δsafety +0.00
- **safety_driver:** temp | decrease_temp | Δsafety -1.04
- **safety_driver:** ribose | remove_sugars | Δsafety +0.00
- **safety_driver:** cysteine | remove_amino_acids | Δsafety +0.00

### Predicted Desirable Targets
| Compound | Barrier |
| :--- | :--- |
| Hydrogen Sulfide | - |
| 2-Furfurylthiol (FFT) | - |
| 2-Methyl-3-furanthiol (MFT) | - |

### Projection Calibration

| Compound | Panel Class | Role | Kind | Evidence State | Reachability | Process State | Retention | Calibration Source | Observable Assumption | Evidence | Browning | Observable ppb |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | ---: |
| furfural | severity_markers | constrained | volatile | conditional_calibration | chemically_reachable | intermediate_matrix | static_class_profile | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer) | static_class_profile | class_level | standard_matrix_support | class_anchored | 0.36 | 0.11 |
| 2-furfurylthiol | sulfur_positives | scored | volatile | externally_benchmarked | chemically_reachable | intermediate_matrix | sulfur_binding_prior | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior | class_level | standard_matrix_support | directional_transferred | 0.36 | 0.00 |
| 3-methylbutanal | strecker_aldehydes | scored | volatile | externally_benchmarked | chemically_reachable | intermediate_matrix | static_class_profile | Pratap-Singh 2021 pea isolate ambient slurry baseline (generic aldehyde transfer) | static_class_profile | class_level | standard_matrix_support | class_anchored | 0.36 | 0.00 |
| 2-methyl-3-furanthiol | sulfur_positives | scored | volatile | externally_benchmarked | chemically_reachable | intermediate_matrix | sulfur_binding_prior | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior | class_level | standard_matrix_support | directional_transferred | 0.36 | 0.00 |
| bis(2-methyl-3-furyl) disulfide | sulfur_positives | scored | volatile | transferred_prior | conditionally_reachable | intermediate_matrix | sulfur_binding_prior | Interpolated base sulfur yield matching internal benchmark limits | sulfur_binding_prior | class_level | standard_matrix_support | directional_transferred | 0.36 | 0.00 |
| 2,5-dimethylpyrazine | pyrazines | constrained | volatile | externally_benchmarked | chemically_reachable | intermediate_matrix | pyrazine_binding_prior | Interpolated base pyrazine yield matching internal benchmark limits | pyrazine_binding_prior | class_level | standard_matrix_support | directional_transferred | 0.36 | 0.00 |
### Trust Surface
- **decision_mode:** directional_hypothesis
- **benchmark_neighborhood:** matrix_intake_only
- **extrapolation_axes:** benchmark_neighborhood, matrix, precursors
- **top_calibration_source:** Pratap-Singh 2021 pea isolate ambient slurry baseline (generic furan transfer)
- **top_observable_assumption:** static_class_profile | class_level | standard_matrix_support
- **top_reachability_status:** chemically_reachable

### Flavor Axis Diagnostics
- **strecker_balance_score:** 0.00
- **strecker_gap_penalty:** 0.00
- **pyrazine_signal_ppb:** 0.00
- **pyrazine_propensity:** 0.65
- **pyrazine_burden:** 0.00
- **pyrazine_penalty:** 0.00
- **furanone_support_score:** 0.00
- **furanone_penalty:** 0.35
- **thiamine_pathway_active:** False
- **thiamine_availability_source:** native_matrix_default_inactive
- **thiamine_availability_explicit:** False
- **thiamine_provenance_mode:** inactive
- **effective_runtime_ph:** 5.50
- **dominant_donor_class:** pentose
- **donor_limited:** False
- **donor_pool_factors:** ribose=1.20
- **furanone_expected:** HEMF, DMHF
- **furanone_missing:** HEMF, DMHF
- **lincoln_crosstalk_prior:** Lincoln 2025 qualitative prior inactive for this formulation context.
- **active_family_lanes:** Amino acid-sugar Maillard core, Reducing sugar and carbonyl donor hierarchy, Protein Damage Markers, Nucleotide degradation and ribose support, Melanoidin Polymerization, Carbohydrate pyrolysis and caramelization
- **family_lane_01:** Amino acid-sugar Maillard core | posture=first_class_core | The amino acid-sugar core lane is active and already supports sulfur-positive routing signals.
- **family_lane_07:** Reducing sugar and carbonyl donor hierarchy | posture=immediate_expansion_lane | Dominant donor class is pentose, and donor identity is now constrained by pyrazine-control references plus a bounded prior instead of a fixed sugar heuristic.
- **family_lane_12:** Protein Damage Markers | posture=first_class_runtime_lane | Protein damage markers now travel as a first-class guardrail lane, so reactive lysine loss, AGE proxies, and furosine crossover remain visible alongside aroma-positive lanes instead of staying buried inside aggregate safety output.
- **family_lane_04:** Nucleotide degradation and ribose support | posture=high_value_support_lane | Explicit ribose delivery is active, so Family 04 is supporting downstream chemistry without preserved nucleotide-driven umami support.
- **family_lane_16:** Melanoidin Polymerization | posture=trapping_burden_modifier | Browning burden is high enough that melanoidin build-up should be treated as a bounded trapping lane, reducing sulfur-positive observability rather than being interpreted as uniformly helpful browning.
- **family_lane_09:** Carbohydrate pyrolysis and caramelization | posture=failure_mode_lane | Carbohydrate pyrolysis or caramelization markers are active, and the lane now carries explicit furanone and carbonyl anchors so browning support is separated from benchmark-visible over-furan drift.
- **family_target_score_delta:** 0.28
- **family_maillard_closure_delta:** -0.01
- **family_off_flavour_risk_delta:** 0.01
- **family_prior_bundle:** alternative_protein_matrix_scope=11, amino_acid_sugar_core=1, ascorbic_acid_maillard=14, carbohydrate_pyrolysis_and_caramelization=7, carbonyl_donor_hierarchy=17, glutathione_and_peptide_support=3, lipid_maillard_crosstalk=8, lipid_oxidation_and_carbonylic_crosstalk=34, melanoidin_polymerization=10, nucleotide_and_ribose_support=10, off_note_and_maillard_suppression=1, phospholipid_amine_sink=8, polyphenol_amino_capping=17, protein_damage_markers=1, thiamine_fragmentation_support=8
- **lipid_benchmark_ready_targets:** 1-Hexanol, 2-Pentylfuran, Hexanal, Nonanal
- **lipid_benchmark_marker_targets_ug_per_l:** hexanal=782.0, 2-pentylfuran=163.0, nonanal=24.0, methional=3.1, 2,5-dimethylpyrazine=2.3
- **lipid_crosstalk_priors:** lincoln_2025_polyphenol_crosstalk_v1, ding_2020_schiff_base_amadori_emulsion_rates, hidalgo_2007_decadienal_phenylalanine_styrene, richards_2009_hemoglobin_liposome_oxidation, smagghe_2006_leghemoglobin_oxygen_dissociation, zamora_2010_decadienal_asparagine_decarboxylation
- **lipid_maillard_closure_pressure:** 0.00
- **lipid_maillard_kinetic_priors:** hexanal_radical_quench
- **lipid_hexanal_suppression_fraction:** 0.00
- **protein_damage_benchmark_anchors:** acs_2022_pba_lysine_loss_benchmark, acs_ref3_spi_acrylamide_fast_kinetics, foods_2023_cml_cel_proxy_benchmark, ramirez_jimenez_2000_furosine_crossover_benchmark
- **protein_damage_burden_score:** 0.05
- **state_marker_thiamineavailability:** Thiamine availability | role=diagnostic | kind=state_variable | influence=upstream_state_only | summary=available=False, source=native_matrix_default_inactive, mode=inactive, mft_fraction_estimate=0.00
- **state_marker_nucleotideenrichment:** Nucleotide enrichment | role=diagnostic | kind=state_variable | influence=upstream_state_only | summary=support_score=0.40, nucleotide_active=False, ribose_delivery_active=True
- **state_marker_freeaminoacidenrichment:** Free amino acid enrichment | role=diagnostic | kind=state_variable | influence=upstream_state_only | summary=precursor_release_active=False, pretreatment_support_score=0.18
- **state_marker_pretreatmentphshift:** Pretreatment pH shift | role=diagnostic | kind=state_variable | influence=upstream_state_only | summary=pretreatment_pH_shift_active=False, fermentation_cues_active=False
- **state_marker_protein damage burden:** protein_damage_burden | role=report_only | kind=state_variable | influence=guardrail_lane | summary=damage_burden_score=0.05, reactive_lysine_loss_pct=4.2, selected_anchor=acs_2022_pba_lysine_loss_benchmark
- **state_marker_polyphenol precursor sink:** polyphenol_precursor_sink | role=report_only | kind=state_variable | influence=upstream_precursor_sink | summary=quinone_budget=0.00, cysteine_depletion_factor=0.00, lysine_depletion_factor=0.00
- **state_marker_ascorbic dicarbonyl source:** ascorbic_dicarbonyl_source | role=report_only | kind=state_variable | influence=bounded_upstream_source | summary=dicarbonyl_source_pressure=0.00, ascorbic_acid_loss=0.00, pentosidine_load=0.00
- **state_marker_pe glycation fraction:** pe_glycation_fraction | role=report_only | kind=state_variable | influence=upstream_precursor_sink | summary=pe_glycation_fraction=0.00, sugar_sink_fraction=0.00, sugar_retention_factor=1.00
- **state_marker_melanoidin trapping burden:** melanoidin_trapping_burden | role=report_only | kind=state_variable | influence=bounded_trapping_lane | summary=melanoidin_mass=0.15, thiol_scavenging_factor=0.04, fft_fold_reduction_anchor=0.0
- **state_marker_caramelizationseverity:** Caramelization severity | role=diagnostic | kind=state_variable | influence=upstream_state_plus_marker_panel | summary=severity_signal_ppb=0.11, severity_penalty_factor=0.00
## 4. Analytical Metadata
### Matrix Explainability
- **protein_type:** pea_iso
- **effective_denaturation_state:** 0.9214887890914122
- **temperature_celsius:** 105.0
- **time_minutes:** 45.0
- **pH:** 5.5
- **lysine_accessibility:** 0.43057865469096945
- **cysteine_accessibility:** 0.06528932734548473
- **bulk_volatile_retention:** 0.5674382062546259
- **literature_window:** {'lysine_min': 0.3, 'lysine_max': 0.45, 'cysteine_min': 0.0, 'cysteine_max': 0.08, 'source': 'Asen 2022 DSC/DTNB + Li 2025 Ellman SH envelope, retained as a conservative pea-isolate interpolation while exact benchmark-condition values remain unmeasured'}
- **denaturation_source:** Asen 2022 DSC/DTNB thermal window + Li 2025 free-SH response; calibrated so pea isolate stays mostly native at 40C but opens progressively by 90-140C
- **prior_summary:** {'accessibility_window': {'parameter': 'accessibility_window', 'source': 'Asen 2022 DSC/DTNB + Li 2025 Ellman SH envelope, retained as a conservative pea-isolate interpolation while exact benchmark-condition values remain unmeasured', 'provenance_tier': 'literature_derived_transfer', 'confidence_tier': 'medium', 'process_state_applicability': ['ambient_slurry', 'heated_matrix', 'aqueous_pre_extrusion_model', 'extrusion_structured'], 'notes': 'Directly informed by validated pea process-state anchors, but still transferred to the benchmark condition rather than measured there.'}, 'denaturation_heuristic': {'parameter': 'denaturation_heuristic', 'source': 'Asen 2022 DSC/DTNB thermal window + Li 2025 free-SH response; calibrated so pea isolate stays mostly native at 40C but opens progressively by 90-140C', 'provenance_tier': 'literature_derived_transfer', 'confidence_tier': 'medium', 'process_state_applicability': ['ambient_slurry', 'heated_matrix', 'aqueous_pre_extrusion_model', 'extrusion_structured'], 'notes': 'This is a low-cost interpolation heuristic derived from pea literature anchors, not a direct benchmark fit.'}, 'volatile_class_profile': {'parameter': 'volatile_class_profile', 'source': 'Literature-calibrated class-aware trapping: VSCs bind strongly to pea proteins as verified by Sun 2025/2026 and Lozano 2009.', 'provenance_tier': 'literature_derived_transfer', 'confidence_tier': 'medium', 'process_state_applicability': ['ambient_slurry', 'heated_matrix', 'aqueous_pre_extrusion_model', 'extrusion_structured'], 'notes': 'Used as a class-level surrogate until compound-specific matrix retention data exists for meaty sulfur targets.'}, 'matrix_correction': {'parameter': 'matrix_correction', 'source': 'Asen 2022 DSC/DTNB + Li 2025 Ellman SH envelope; values remain a conservative PPI interpolation because the exact 95C pH 5.5 benchmark condition still lacks direct wet-lab measurement', 'provenance_tier': 'literature_derived_transfer', 'confidence_tier': 'medium', 'process_state_applicability': ['ambient_slurry', 'heated_matrix', 'aqueous_pre_extrusion_model', 'extrusion_structured'], 'notes': 'This is a reliability improvement over a generic heuristic because the numbers now explicitly track the validated pea process-state literature.'}}
- **matrix_prior_process_state_applicability:** ['ambient_slurry', 'aqueous_pre_extrusion_model', 'extrusion_structured', 'heated_matrix']
- **accessibility_profile:** free_like
- **accessibility_warning:** False
- **accessibility_dominant_source:** estimated_from_conditions
- **jacket_temperature_celsius:** 105.0
- **effective_temperature_celsius:** 105.0
- **moisture_regime:** hme
- **sme_kj_per_kg:** 0.0
- **sterilization_temperature_celsius:** None
- **sterilization_time_minutes:** 0.0
- **extrusion_process:** {'active': True, 'thermal_exposure': True, 'thermal_exposure_basis': ['barrel_temperature_above_ambient_threshold'], 'ambient_thermal_threshold_celsius': 50.0, 'model': 'sequential_isothermal_zones', 'moisture_regime': 'hme', 'water_activity': 0.85, 'sme_kj_per_kg': 0.0, 'screw_speed_rpm': None, 'feed_rate_kg_per_h': None, 'jacket_temperature_celsius': 105.0, 'sme_temperature_offset_celsius': 0.0, 'effective_temperature_celsius': 105.0, 'die_exit_temperature_celsius': 70.0, 'mean_residence_time_seconds': 44.1, 'relative_rtd_spread': 0.24188920068170988, 'pre_extrusion_damage': {'furosine_mg_per_kg': 18.0, 'lal_mg_per_kg': 68.0}, 'sterilization': {'enabled': False, 'temperature_celsius': None, 'time_minutes': 0.0, 'damage_increment': {'furosine_mg_per_kg': 0.0, 'lal_mg_per_kg': 0.0}}, 'process_damage_increment': {'furosine_mg_per_kg': 0.0, 'lal_mg_per_kg': 0.0}, 'damage_load_attribution': {'ingredient_baseline': {'furosine_mg_per_kg': 18.0, 'lal_mg_per_kg': 68.0}, 'process_increment': {'furosine_mg_per_kg': 0.0, 'lal_mg_per_kg': 0.0}, 'note': 'total_damage_load = ingredient_baseline + process_increment. With no thermal exposure the process increment is exactly zero, so any non-zero total is ingredient provenance, not processing.'}, 'total_damage_load': {'furosine_mg_per_kg': 18.0, 'lal_mg_per_kg': 68.0}, 'zone_profile': [{'zone_index': 1.0, 'temperature_celsius': 85.0, 'time_fraction': 0.2}, {'zone_index': 2.0, 'temperature_celsius': 95.0, 'time_fraction': 0.35}, {'zone_index': 3.0, 'temperature_celsius': 105.0, 'time_fraction': 0.45}], 'rtd_assessment': {'decision': 'sequential_zone_sufficient_for_current_use_case', 'reason': 'For the current HME scientist workflow, the sequential-zone model already captures the dominant thermal history while screw-geometry-free RTD would add parameters the repo cannot yet benchmark.', 'mean_residence_time_seconds': 44.1, 'relative_spread': 0.24188920068170988}, 'volatile_transport': {'panel': {'2-methyl-3-furanthiol': {'compound': '2-methyl-3-furanthiol', 'compound_class': 'sulfur', 'relative_vapor_pressure_proxy': 0.92, 'shear_release_factor': 0.9971599560994293, 'die_stripping_fraction': 0.0785743159037375, 'combined_headspace_factor': 0.9188087947023157}, '2-furfurylthiol': {'compound': '2-furfurylthiol', 'compound_class': 'sulfur', 'relative_vapor_pressure_proxy': 0.92, 'shear_release_factor': 0.9971599560994293, 'die_stripping_fraction': 0.0785743159037375, 'combined_headspace_factor': 0.9188087947023157}, 'Hexanal': {'compound': 'Hexanal', 'compound_class': 'aldehyde', 'relative_vapor_pressure_proxy': 1.0, 'shear_release_factor': 0.9973403006423134, 'die_stripping_fraction': 0.182048443050568, 'combined_headspace_factor': 0.8157760517187949}, '2-Pentylfuran': {'compound': '2-Pentylfuran', 'compound_class': 'furan', 'relative_vapor_pressure_proxy': 0.78, 'shear_release_factor': 0.9983654250148758, 'die_stripping_fraction': 0.08512217556238229, 'combined_headspace_factor': 0.9133823880313471}, '2,5-Dimethylpyrazine': {'compound': '2,5-Dimethylpyrazine', 'compound_class': 'pyrazine', 'relative_vapor_pressure_proxy': 0.36, 'shear_release_factor': 0.9990750025127441, 'die_stripping_fraction': 0.012131883789173688, 'combined_headspace_factor': 0.986954340685591}}, 'dominant_tradeoff': 'die_loss_vs_precursor_release'}}

## 5. Provenance
- **artifact_kind:** single_run_report
- **generated_at:** 2026-09-01T22:10:58.462156
- **generator:** scripts/generators/generate_report_visual_examples.py --output-dir results/validation/report_visual_examples --docs-asset-dir docs/assets
- **repository:** workspace | branch cleaning | commit 2208fd0 | dirty True
- **input_fingerprint_sha256:** b63663346583a45dfd5cb20eeecab6c9a8d6d27baa6bb1a9b6b2324eeeee6e70
- **scientific_surface:**
  - scientific_reference: docs/reference/SCIENTIFIC_REFERENCE.md
  - benchmark_summary: results/validation/benchmark_summary.md
  - validated_envelope: results/validation/validated_envelope.md
  - validation_overview: results/validation/validation_overview.md
  - matrix_decision_panel: data/lit/matrix_decision_panel.json
  - chemistry_family_scope_registry: data/lit/chemistry_family_scope_registry.json
  - family_ingestion_plan_registry: data/lit/family_ingestion_plan.json
  - family_identifier_contract: results/validation/family_identifier_contract.md
  - family_identifier_contract_json: results/validation/family_identifier_contract.json
  - family_strategy_policy: results/validation/family_strategy_policy.md
  - family_strategy_policy_json: results/validation/family_strategy_policy.json
  - family_payload_coverage: results/validation/family_payload_coverage.md
  - family_payload_coverage_json: results/validation/family_payload_coverage.json
  - matrix_family_coverage_registry: data/lit/matrix_family_coverage_registry.json
  - benchmark_intake_registry: data/lit/benchmark_intake_registry.json
  - computational_priors: data/lit/computational_priors.json
  - slr_incorporation_matrix: results/literature/slr_incorporation_matrix.json
  - flavor_reference_payloads: data/lit/flavor_reference_payloads.json
  - process_state_calibrations: data/lit/process_state_calibrations.json
  - retention_reference_payloads: data/lit/retention_reference_payloads.json
  - process_gap_registry: data/lit/process_gap_registry.json
  - safety_reference_payloads: data/lit/safety_reference_payloads.json
  - primary_benchmark_protocol: docs/protocols/PPI_SPI_PRIMARY_BENCHMARK_PROTOCOL.md
  - primary_benchmark_contract: data/protocols/ppi_spi_primary_benchmark_contract.json
  - literature_learning_loop: results/validation/literature_learning_loop.md
  - literature_learning_loop_json: results/validation/literature_learning_loop.json
  - literature_backlog: results/validation/literature_backlog.md
  - literature_backlog_json: results/validation/literature_backlog.json
  - deep_research_runtime_queue: results/validation/deep_research_runtime_queue.md
  - deep_research_runtime_queue_json: results/validation/deep_research_runtime_queue.json
  - literature_runtime_templates: results/validation/literature_runtime_templates.json
  - family_ingestion_plan: results/validation/family_ingestion_plan.md
  - family_ingestion_plan_json: results/validation/family_ingestion_plan.json
  - matrix_target_status: results/validation/matrix_target_status.md
  - matrix_target_status_json: results/validation/matrix_target_status.json
  - chemistry_family_scope: results/validation/chemistry_family_scope.md
  - chemistry_family_scope_json: results/validation/chemistry_family_scope.json
  - matrix_family_coverage: results/validation/matrix_family_coverage.md
  - matrix_family_coverage_json: results/validation/matrix_family_coverage.json
  - family_sensitivity: results/validation/family_sensitivity.md
  - family_sensitivity_json: results/validation/family_sensitivity.json
  - family_lane_validation: results/validation/family_lane_validation.md
  - family_lane_validation_json: results/validation/family_lane_validation.json
  - refinement_surrogate_patches: data/lit/refinement_surrogate_patches.json

## 6. Glossary
Plain-language meaning of the labels used above. The model is honest about *how* it knows what it claims; this section names that vocabulary.

**Three different tier vocabularies appear in this report. They are not interchangeable.** Each answers a different question:

**1. `tier` — how well benchmark-supported is *this run*?** Emitted on the run-level *Confidence & Support* block, on every *Compound Confidence* row, and on every *Aggregate Sensory Confidence* row. It is a band on a 0-100 confidence score, and it always travels with a `prediction_mode`:

| `tier` | score band | paired `prediction_mode` | what it licenses |
| :--- | :--- | :--- | :--- |
| `high` | >= 85 | `benchmark_supported_quantitative` | quantitative prioritisation before wet-lab confirmation |
| `medium` | 65-84 | `ranking_supported` | ranking and triage; verify absolute levels experimentally |
| `low` | 45-64 | `directional_only` | direction only; absolute ppb is provisional |
| `exploratory` | < 45 | `hypothesis_only` | hypothesis generation, not decision-grade |

**2. `calibration_evidence_strength` (shown as *Evidence* in the Calibration Summary) — what kind of anchor stands behind a compound's projection?** Strongest to weakest, and only these values are emitted at compound level:
- `literature_anchored` — a published measurement on a directly comparable system backs this compound's retention/partition treatment.
- `conditional_literature_anchored` — literature-anchored, but only under stated conditions (pH / process-state caveats attached).
- `class_anchored` — anchored at compound-class level (e.g. "sulfur volatiles"), not for this molecule specifically.
- `directional_transferred` — a prior transferred from an adjacent matrix or process state; direction is meaningful, magnitude is not.
- `process_state_mismatch` — an anchor exists, but for a different process state than the one you asked about; the nearest state was substituted.
- `heuristic` — no anchor at all; a built-in class default. Ranking use only.

  When the run is out of calibration scope every one of these is demoted one notch (see the scope banner below), and the pre-demotion value is preserved as `scope_demoted_from` in `report.json`.

**3. `confidence_tier` — how strong is the *literature prior* behind a chemistry lane?** This one comes from the curated literature registries (`data/lit/`), not from your run, and uses a five-point scale: `high`, `medium_high`, `medium`, `medium_low`, `low`. It grades the source, not the prediction: a `high` `confidence_tier` prior can still feed an `exploratory` `tier` prediction, because your formulation may sit far from where that prior was measured.

  *Name collision, stated plainly:* the campaign/comparison JSON also carries a key spelled `confidence_tier` that holds the run-level `tier` value (vocabulary 1), kept as a legacy alias. Prefer the `tier` key alongside it; only `confidence_tier` inside `data/lit/` payloads means vocabulary 3.

**Scope banner.** A `⚠️ Out of calibration scope` banner at the top of the report means the matrix or process you asked about lies outside the convex hull of formulations the model has been calibrated against. The predictions are still emitted, but every compound's evidence strength is demoted one notch.

**Reachability** (`reachability_status`). `chemically_reachable` — the compound is produced by an enumerated, barrier-scored pathway from your precursors. `conditionally_reachable` — reachable only under an assumption the run had to make. `merely_plausible` — no enumerated route; the number is a class-level projection.

**Observable assumption** (`observable_assumption_summary`). A pipe-joined triple: retention runtime mode, calibration fallback mode, support origin — e.g. `static_class_profile | class_level | standard_matrix_support` means the volatile's matrix retention came from a static class profile, its calibration fell back to class level, and no special matrix-support route was used.

**Confidence envelope.** `0.038 ppb [0.012-0.089, 90% CI]` means the p50 (median) Monte-Carlo prediction is `0.038 ppb`, with the central 90 % of samples between `0.012` and `0.089 ppb`. A compound printed without an interval had no envelope sampled. Wide envelopes make coverage cheap — read coverage and width together.

**Intervention waterfall.** When two formulations are compared, the per-compound delta is decomposed into class-level (e.g. "thiols", "aldehydes") and per-precursor (e.g. "cysteine", "glutathione") contributions. Per-precursor attribution sums to the compound delta and is explicit about attribution mode.

Full machine-readable trust artifacts: `results/validation/`. Per-compound 90 % envelope: `results/validation/prediction_uncertainty.md`. External hold-out: `results/validation/external_validation_report.md`.

## 7. Recommended next experiment
Ranked by value-of-information (envelope miss × CI width × ODT-anchored decision relevance). Source: `results/validation/experiment_value_ranking.json`.

Scoped to matrix `pea_iso` (filtered from the global ranking).

| Rank | VoI | Benchmark | Matrix | Compound | DoE template | Why this one |
| ---: | ---: | --- | --- | --- | --- | --- |
| 21 | 0.97 | `pea_isolate_ribose_cysteine_100C_45min_Internal2026` | `pea_iso` | 2,5-dimethylpyrazine | `missing_kinetic_dataset` | CI width 6.47 dex; ≈2e-08× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |
| 22 | 0.97 | `pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026` | `pea_iso` | 2,5-dimethylpyrazine | `missing_kinetic_dataset` | CI width 6.47 dex; ≈2e-08× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |
| 25 | 0.68 | `pea_isolate_ribose_cysteine_100C_45min_Internal2026` | `pea_iso` | bis(2-methyl-3-furyl) disulfide | `missing_kinetic_dataset` | CI width 4.56 dex; ≈0.06× ODT (decision_relevance=0.50); wide envelope — time-course narrows the rate-limiting step |

_How to use this: run `./scripts/docker_maillard.sh next-experiment --top 3` to materialise pre-filled intake YAMLs and protocol Markdown for each row. Ingest the resulting measurement via `./scripts/docker_maillard.sh ingest --file results.csv ...`._

## 8. Sensory readout
Per-compound odour activity value (OAV = predicted ppb ÷ curated odour threshold). An axis is *above threshold* when at least one of its compounds reaches OAV ≥ 1. Compounds without a curated odour threshold are listed but excluded from axis roll-ups. Source thresholds: `data/species/desirable_targets.yml`, `data/species/off_flavour_targets.yml`.

**Headline:** meaty below threshold; off-notes no ODT; safety no ODT.

### Axis roll-up
| Axis | Compounds (with ODT) | Above threshold | Max OAV | Top contributor |
| --- | ---: | ---: | ---: | --- |
| meaty | 4 | 0 | 0.210 | 2-furfurylthiol |
| off-note | 0 | 0 | n/a | _no compound with curated ODT in this run_ |
| safety | 0 | 0 | n/a | _no compound with curated ODT in this run_ |

### Per-compound OAV (90 % CI)
| Compound | Axis | ODT (μg/kg) | Predicted ppb (p50) | OAV (p50) | OAV p5 | OAV p95 | ≥1? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| 2-furfurylthiol | meaty | 0.01 | 0.0021 | 0.210 | 4.60e-04 | 23.2 | · |
| 2-Furfurylthiol (FFT) | meaty | 0.01 | 0.0021 | 0.210 | 4.60e-04 | 23.2 | · |
| 2-methyl-3-furanthiol | meaty | 0.007 | 0.000252 | 0.036 | 3.25e-04 | 3.97 | · |
| 2-Methyl-3-furanthiol (MFT) | meaty | 0.007 | 0.000252 | 0.036 | 3.25e-04 | 3.97 | · |
| bis(2-methyl-3-furyl) disulfide | unclassified | 0.02 | 5.51e-05 | 2.76e-03 | 2.63e-06 | 0.305 | · |
| Bis(2-methyl-3-furyl) disulfide | unclassified | 0.02 | 5.51e-05 | 2.76e-03 | 2.63e-06 | 0.305 | · |
| 3-methylbutanal | unclassified | 1.5 | 0.000903 | 6.02e-04 | n/a | n/a | · |
| 3-Methylbutanal | unclassified | 1.5 | 0.000903 | 6.02e-04 | n/a | n/a | · |
| furfural | unclassified | 3000 | 0.108 | 3.58e-05 | 3.24e-07 | 3.96e-03 | · |
| Furfural | unclassified | 3000 | 0.108 | 3.58e-05 | 3.24e-07 | 3.96e-03 | · |
| 2,5-dimethylpyrazine | unclassified | 1800 | 4.3e-06 | 2.39e-09 | 6.61e-13 | 1.93e-06 | · |
| 2,5-Dimethylpyrazine | unclassified | 1800 | 4.3e-06 | 2.39e-09 | 6.61e-13 | 1.93e-06 | · |
| O=Cc1ccco1 | unclassified | n/a | 0.108 | n/a | n/a | n/a | — |
| CC(C)CC=O | unclassified | n/a | 0.000903 | n/a | n/a | n/a | — |
| Cc1cnc(C)cn1 | unclassified | n/a | 4.3e-06 | n/a | n/a | n/a | — |
| SCc1ccco1 | unclassified | n/a | 0.0021 | n/a | n/a | n/a | — |
| Cc1occc1S | unclassified | n/a | 0.000252 | n/a | n/a | n/a | — |
| Cc1occc1SSc1ccoc1C | unclassified | n/a | 5.51e-05 | n/a | n/a | n/a | — |

_6/18 predicted compounds have no curated odour threshold; they appear in the per-compound table but do not contribute to axis roll-ups._

