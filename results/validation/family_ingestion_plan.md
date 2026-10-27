# Family Ingestion Plan

| SLR | Family | Strategic Posture | Runtime Concept | Payload Types | Wave | Next Build Action |
| --- | --- | --- | --- | --- | ---: | --- |
| 01 | amino_acid_sugar_core | first_class_core | core_reaction_network | benchmark_payload, flavor_reference_payload | 0 | Keep as the quantitative trunk and dependency map. |
| 02 | lipid_oxidation_and_carbonylic_crosstalk | immediate_expansion_lane | lipid_crosstalk_lane | benchmark_payload, retention_payload, computational_prior | 1 | Split into adverse-marker and carbonyl-competition sub-lanes. |
| 07 | carbonyl_donor_hierarchy | immediate_expansion_lane | carbonyl_donor_hierarchy | benchmark_payload, flavor_reference_payload | 1 | Encode donor identity as a first-class formulation variable. |
| 10 | fermentation_pretreatment | upstream_pretreatment_lane | fermentation_pretreatment_node | benchmark_payload, flavor_reference_payload, computational_prior | 1 | Implement fermentation as an upstream state modifier. |
| 08 | off_note_and_maillard_suppression | guardrail_lane | off_note_and_maillard_suppression | benchmark_payload, computational_prior, safety_payload | 1 | Represent inhibitor loads as family-aware priors. |
| 11 | lipid_maillard_crosstalk | first_class_runtime_lane | lipid_maillard_competition | benchmark_payload, computational_prior | 1 | Implement the hexanal/MFT competition surface. |
| 12 | protein_damage_markers | first_class_runtime_lane | safety_damage_lane | safety_payload, guardrail_entry | 1 | Integrate AGE/ALE kinetics into the safety module. |
| 13 | polyphenol_amino_capping | upstream_precursor_sink | precursor_depletion_sink | upstream_modifier_payload, computational_prior | 1 | Implement the pH-dependent quinone-thiol kinetics. |
| 03 | thiamine_fragmentation_support | high_value_support_lane | thiamine_fragmentation_support | computational_prior, flavor_reference_payload | 2 | Model thiamine first as an availability modifier. |
| 04 | nucleotide_and_ribose_support | high_value_support_lane | nucleotide_and_ribose_support | flavor_reference_payload, computational_prior | 2 | Encode nucleotide enrichment state variables explicitly. |
| 05 | glutathione_and_peptide_support | high_value_support_lane | sulfur_peptide_support | computational_prior, flavor_reference_payload | 2 | Model peptide-bound sulfur as a support lane. |
| 14 | ascorbic_acid_maillard | bounded_upstream_source | stealth_dicarbonyl_source | bounded_source_term_payload, computational_prior | 2 | Implement the biphasic AA functionality. |
| 15 | phospholipid_amine_sink | upstream_precursor_sink | sugar_depletion_sink | upstream_modifier_payload, computational_prior | 2 | Implement the 1:2 sugar-PE stoichiometry multiplier. |
| 16 | melanoidin_polymerization | trapping_burden_modifier | matrix_trapping_node | trapping_burden_modifier_payload, computational_prior | 2 | Implement pseudo-zero-order melanoidin accumulation. |
| 09 | carbohydrate_pyrolysis_and_caramelization | failure_mode_lane | carbohydrate_pyrolysis_lane | benchmark_payload, flavor_reference_payload, computational_prior | 3 | Use first as a severity and failure-mode lane. |
| 06 | alternative_protein_matrix_scope | matrix_scope_lane | matrix_family_scope_extension | computational_prior, process_state_calibration | 3 | Use to extend matrix scope artifacts and bounded priors. |

## Deep Research Priority Surface

| Rank | SLR | Family | Backlog Citations | 8/8 Citations | Impact | Ease | Priority | Next Slice |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 12 | Protein Damage Markers | 2 | 2 | 1.00 | 0.95 | 85.4 | expand safety damage proxies with benchmark-backed guardrail payloads |
| 2 | 11 | Maillard/Lipid Crosstalk | 2 | 5 | 1.00 | 0.93 | 84.8 | wire lipid-Maillard competition constants into the runtime competition surface |
| 3 | 10 | Microbial fermentation pretreatment | 5 | 2 | 0.88 | 0.98 | 84.0 | bind fermentation pretreatment state to precursor release and pH-shift effects |
| 4 | 08 | Plant off-notes and Maillard suppression | 3 | 1 | 0.78 | 0.98 | 77.9 | Represent inhibitor loads as family-aware priors. |
| 5 | 05 | Glutathione and peptide support | 5 | 2 | 0.82 | 0.83 | 76.8 | bound peptide and glutathione support with reusable sulfur-lane priors |
| 6 | 03 | Thiamine degradation and sulfur support | 1 | 4 | 0.82 | 0.78 | 70.6 | promote thiamine from availability cue to calibrated sulfur-support modifier |
| 7 | 04 | Nucleotide degradation and ribose support | 1 | 3 | 0.82 | 0.78 | 70.6 | encode nucleotide and ribose enrichment as explicit upstream state variables |
| 8 | 13 | Polyphenol-Amino Capping | 4 | 2 | 0.60 | 0.93 | 69.9 | Implement the pH-dependent quinone-thiol kinetics. |
| 9 | 09 | Carbohydrate pyrolysis and caramelization | 3 | 4 | 0.68 | 0.80 | 68.0 | Use first as a severity and failure-mode lane. |
| 10 | 15 | PE Stealth Sugar Sink | 4 | 2 | 0.60 | 0.78 | 65.4 | Implement the 1:2 sugar-PE stoichiometry multiplier. |
| 11 | 06 | Alternative protein matrix scope | 1 | 5 | 0.68 | 0.78 | 64.3 | extend matrix-family scope without over-claiming transferred closure |
| 12 | 16 | Melanoidin Polymerization | 4 | 3 | 0.54 | 0.78 | 62.7 | Implement pseudo-zero-order melanoidin accumulation. |
| 13 | 14 | Ascorbic Acid Maillard | 3 | 4 | 0.58 | 0.73 | 61.4 | bound ascorbic-acid dicarbonyl source terms before they leak into safety claims |

## Target Modules

| SLR | Target Modules | Target Observables or State Variables |
| --- | --- | --- |
| 01 | none | 2-methyl-3-furanthiol; 2-furfurylthiol; methional; 2,5-dimethylpyrazine |
| 02 | src/lipid_oxidation.py | hexanal; 2-pentylfuran; nonanal; carbonyl_competition_factor |
| 07 | none | ribose; xylose; glucose; fructose |
| 10 | src/matrix_prior_registry.py | free_amino_acid_enrichment; nucleotide_enrichment; pretreatment_pH_shift |
| 08 | none | dicarbonyl_trapping_factor; amino_group_blocking_factor |
| 11 | src/lipid_oxidation.py | hexanal; MDA; 4-HNE; radical_quenching_rate |
| 12 | src/extrusion.py | CML; CEL; furosine; lysinoalanine |
| 13 | none | quinone_budget; cysteine_depletion_factor |
| 03 | none | thiamine_availability; 2-methyl-3-furanthiol |
| 04 | none | IMP; GMP; ribose-5-phosphate |
| 05 | none | glutathione_support_factor; kokumi_support_signal |
| 14 | none | ascorbic_acid_loss; 3-deoxythreosone; pentosidine_load |
| 15 | none | PE_glycation_fraction; available_sugar_subtraction |
| 16 | src/extrusion.py | melanoidin_mass; TH_scavenging_factor |
| 09 | none | HMF; furfural; caramelization_severity |
| 06 | src/matrix_family_coverage.py; src/matrix_prior_registry.py | matrix_family_support_posture; sulfur_deficiency_risk |

Families tracked: 16
Recommended first wave: 02, 07, 10, 08, 11, 12, 13
Active build sequence: 01, 02, 07, 10, 08, 11, 12, 13, 03, 04, 05, 14, 15, 16, 09, 06
Backlog-bearing families: 13
Backlog citations mapped to family plan: 38
Recommended runtime queue: 11, 12, 03, 04, 05
Recommended next slice: 11 Maillard/Lipid Crosstalk -> wire lipid-Maillard competition constants into the runtime competition surface
Mapped scope families: 16
Unmapped scope families: none
Policy: extend_the_quantitative_core_by_explicit_family_lanes_with_machine_readable_payloads_not_narrative_only_docs
Identifier policy: scope_family_id_uses_the_same_canonical_chemistry_family_id_as_payloads_and_validation_outputs
Axis policy: chemistry_family_and_matrix_family_are_separate_axes_and_should_not_share_identifier_spaces
