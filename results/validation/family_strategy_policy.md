# Family Strategy Policy

Quantitative trunk: amino_acid_sugar_core
Benchmarked foundation: keep_amino_acid_plus_sugar_core_with_strecker_and_sulfur_scoring_surfaces_as_the_benchmarked_foundation
Default next expansion family: lipid_oxidation_and_carbonylic_crosstalk
Default next expansion reason: In plant-based matrices this family often dominates early off-notes and interacts directly with Maillard intermediates.

## Family Classification

First-class core: amino_acid_sugar_core
High-priority partial lanes: lipid_oxidation_and_carbonylic_crosstalk
Bounded lanes: thiamine_fragmentation_support, nucleotide_and_ribose_support, glutathione_and_peptide_support, alternative_protein_matrix_scope, carbonyl_donor_hierarchy, off_note_and_maillard_suppression, carbohydrate_pyrolysis_and_caramelization, fermentation_pretreatment
Open gaps: none

## Ingestion Contract

Machine-readable only: True
Narrative-only workflow allowed: False
Required surfaces: benchmark_intake_registry, runtime_payloads, process_gap_registry
Policy: new_families_use_the_same_machine_readable_ingestion_contract_as_the_core_and_do_not_create_a_parallel_markdown_only_workflow

## Lipid Dual Lane

Family: lipid_oxidation_and_carbonylic_crosstalk
Observable lane payloads: benchmark_payload, retention_payload
Competition lane payloads: computational_prior, structural_gap_entry
Policy: treat_lipid_oxidation_as_a_dual_lane_of_observable_adverse_markers_plus_carbonyl_competition_and_crosstalk_priors

## Compute Policy

P4 status: infrastructure_not_active_product_objective
DFT policy: selective_dft_is_reserved_for_benchmark_visible_sulfur_carbonyl_transfer_and_ts_sensitive_gaps_after_cheap_first_screening
MLP policy: mlps_remain_bounded_offline_accelerators_until_local_geometry_or_ts_benchmarks_accept_them_on_maillard_relevant_systems

## Lane Positions

| SLR | Family | Classification | Priority | Payload Types | Wave |
| --- | --- | --- | --- | --- | ---: |
| 01 | amino_acid_sugar_core | first_class_core | high | benchmark_payload, flavor_reference_payload | 0 |
| 02 | lipid_oxidation_and_carbonylic_crosstalk | partially_encoded_high_priority | high | benchmark_payload, retention_payload, computational_prior | 1 |
| 03 | thiamine_fragmentation_support | bounded_lane | medium | computational_prior, flavor_reference_payload | 2 |
| 04 | nucleotide_and_ribose_support | bounded_lane | medium | flavor_reference_payload, computational_prior | 2 |
| 05 | glutathione_and_peptide_support | bounded_lane | medium | computational_prior, flavor_reference_payload | 2 |
| 06 | alternative_protein_matrix_scope | bounded_lane | medium | computational_prior, process_state_calibration | 3 |
| 07 | carbonyl_donor_hierarchy | bounded_lane | high | benchmark_payload, flavor_reference_payload | 1 |
| 08 | off_note_and_maillard_suppression | bounded_lane | high | benchmark_payload, computational_prior, safety_payload | 1 |
| 09 | carbohydrate_pyrolysis_and_caramelization | bounded_lane | medium | benchmark_payload, flavor_reference_payload, computational_prior | 3 |
| 10 | fermentation_pretreatment | bounded_lane | high | benchmark_payload, flavor_reference_payload, computational_prior | 1 |
| 11 | lipid_maillard_crosstalk | open_gap | high | benchmark_payload, computational_prior | 1 |
| 12 | protein_damage_markers | open_gap | high | safety_payload, guardrail_entry | 1 |
| 13 | polyphenol_amino_capping | open_gap | medium | upstream_modifier_payload, computational_prior | 1 |
| 14 | ascorbic_acid_maillard | open_gap | medium | bounded_source_term_payload, computational_prior | 2 |
| 15 | phospholipid_amine_sink | open_gap | high | upstream_modifier_payload, computational_prior | 2 |
| 16 | melanoidin_polymerization | open_gap | high | trapping_burden_modifier_payload, computational_prior | 2 |
