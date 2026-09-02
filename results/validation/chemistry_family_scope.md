# Chemistry Family Scope

| Family | Status | Strategic Role | Preferred Ingestion Mode | Priority | Next Best Action |
| --- | --- | --- | --- | --- | --- |
| amino_acid_sugar_core | first_class_core | first_class_runtime_lane | benchmark_or_reference_payload | high | Keep as the strict core and use it as the dependency map for every non-core family extension. |
| lipid_oxidation_and_carbonylic_crosstalk | partially_encoded_high_priority | computational_closure_candidate | dual_lane_benchmark_plus_retention_payload | high | Promote this to the first non-core family sprint by adding explicit benchmark and calibration lanes. |
| thiamine_fragmentation_support | bounded_lane | bounded_support_lane | directional_prior_then_benchmark | medium | Keep thiamine as an upstream availability and routing support lane. |
| nucleotide_and_ribose_support | bounded_lane | bounded_support_lane | reference_payload_plus_priors | medium | Keep nucleotide support explicit as an upstream state modifier and decision-panel lane. |
| glutathione_and_peptide_support | bounded_lane | bounded_support_lane | directional_prior_then_benchmark | medium | Keep peptide-bound sulfur explicit as a support lane. |
| alternative_protein_matrix_scope | bounded_lane | matrix_scope_lane | process_state_calibration_plus_priors | medium | Keep this lane explicit as matrix-scope support and trust-surface logic. |
| carbonyl_donor_hierarchy | bounded_lane | first_class_runtime_lane | benchmark_payload_plus_priors | high | Keep donor identity first-class in runtime and reporting. |
| off_note_and_maillard_suppression | bounded_lane | guardrail_lane | guardrail_payload_plus_priors | high | Keep this family co-equal with positive flavor lanes. |
| carbohydrate_pyrolysis_and_caramelization | bounded_lane | guardrail_lane | safety_or_structural_gap_registry | medium | Keep this family as a severity and failure-mode lane. |
| fermentation_pretreatment | bounded_lane | bounded_support_lane | process_state_calibration_plus_priors | high | Keep pretreatment explicit as an upstream state-modifier lane. |
| lipid_maillard_crosstalk | active_expansion | first_class_runtime_lane | benchmark_linked_suppression_surface | high | Implement the competition surface between sulfur volatiles and lipid off-notes. |
| protein_damage_markers | active_expansion | first_class_runtime_lane | guardrail_severity_surface | high | Promote CML, CEL, and Furosine to first-class safety metrics. |
| polyphenol_amino_capping | active_expansion | computational_closure_candidate | upstream_modifier_payload | medium | Implement the quinone-budget sink with Cys-first priority mapping. |
| ascorbic_acid_maillard | active_expansion | computational_closure_candidate | bounded_source_term_payload | medium | Add the AA-driven dicarbonyl source term with pH and aw-dependent routing. |
| phospholipid_amine_sink | active_expansion | computational_closure_candidate | mechanistic_to_bounded_closure | high | Subtract PE sink capacity using 1:2 sugar:PE stoichiometry multiplier. |
| melanoidin_polymerization | active_expansion | computational_closure_candidate | trapping_burden_modifier_payload | high | Implement pseudo-zero-order melanoidin accumulation with IC50-based flavor scavenging. |

## Why These Families Matter

| Family | Why It Matters | Current Runtime Assets | Ingestion Surfaces |
| --- | --- | --- | --- |
| amino_acid_sugar_core | This is the generative trunk that creates the dicarbonyl and Amadori or Heyns intermediates feeding most downstream meaty chemistry. | data/lit/flavor_reference_payloads.json; src/benchmark_validation.py | benchmark_intake_registry; flavor_reference_payloads; benchmark_validation |
| lipid_oxidation_and_carbonylic_crosstalk | In plant-based matrices this family often dominates early off-notes and interacts directly with Maillard intermediates. | src/lipid_oxidation.py; data/lit/retention_reference_payloads.json; data/lit/computational_priors.json | benchmark_intake_registry; retention_reference_payloads; computational_priors |
| thiamine_fragmentation_support | Fortified precursor systems can unlock sulfur chemistry that the native plant matrix lacks. | data/lit/computational_priors.json; data/lit/flavor_reference_payloads.json | flavor_reference_payloads; computational_priors |
| nucleotide_and_ribose_support | Nucleotides and ribose-delivery systems can improve savory support and alter precursor availability. | data/lit/flavor_reference_payloads.json; data/lit/matrix_decision_panel.json | benchmark_intake_registry; flavor_reference_payloads; matrix_decision_panel |
| glutathione_and_peptide_support | Low-molecular-weight peptides and glutathione can supply bounded sulfur support and kokumi-like effects. | data/lit/flavor_reference_payloads.json; data/lit/computational_priors.json | flavor_reference_payloads; computational_priors |
| alternative_protein_matrix_scope | States clearly where chemistry knowledge is being transferred into alternative matrices under bounded confidence. | data/lit/matrix_family_coverage_registry.json; data/lit/computational_priors.json; src/matrix_prior_registry.py | benchmark_intake_registry; computational_priors; process_gap_registry |
| carbonyl_donor_hierarchy | Flavor recommendations are wrong if ribose, xylose, glucose, and fructose are treated as interchangeable. | data/lit/flavor_reference_payloads.json; data/lit/matrix_decision_panel.json; src/recommend.py | benchmark_intake_registry; flavor_reference_payloads; computational_priors; matrix_decision_panel |
| off_note_and_maillard_suppression | Flavor predictions are unsafe if plant off-note loads and dicarbonyl trapping are not modeled. | data/lit/safety_reference_payloads.json; data/lit/computational_priors.json; data/lit/matrix_decision_panel.json | benchmark_intake_registry; computational_priors; data/lit/safety_reference_payloads.json |
| carbohydrate_pyrolysis_and_caramelization | Severity markers like HMF and furfural help separate useful browning from over-processed drift. | data/lit/flavor_reference_payloads.json; data/lit/matrix_decision_panel.json | benchmark_intake_registry; flavor_reference_payloads; safety_reference_payloads |
| fermentation_pretreatment | Fermentation changes precursor pools and off-note burden before thermal chemistry starts. | data/lit/computational_priors.json; data/lit/process_state_calibrations.json; data/lit/flavor_reference_payloads.json | benchmark_intake_registry; flavor_reference_payloads; computational_priors |
| lipid_maillard_crosstalk | Lipid-derived dicarbonyls (MDA, 4-HNE) compete for protein amines, while Maillard intermediates quench lipid radicals. | data/research_corpus/raw/Gemini_deep_research_11.md; src/lipid_oxidation.py | computational_priors; benchmark_intake_registry |
| protein_damage_markers | Extrusion intensity drives the formation of CML, CEL, Furosine, and LAL, critical for safety labeling. | data/research_corpus/12_protein_damage_markers.md; src/safety.py | safety_reference_payloads; process_gap_registry |
| polyphenol_amino_capping | Oxidized polyphenols (quinones) rapidly cap cysteine and lysine residues, depleting flavor precursors. | data/research_corpus/raw/Gemini_deep_research_13.md | computational_priors; matrix_decision_panel |
| ascorbic_acid_maillard | Ascorbic acid acts as a dicarbonyl and AGE initiator at extrusion temperatures (>120°C). | data/research_corpus/raw/Gemini_deep_research_14.md | computational_priors; safety_reference_payloads |
| phospholipid_amine_sink | PE in lecithin competes with protein amines for sugar with lower Ea and interfacial acceleration. | data/research_corpus/raw/Gemini_deep_research_15.md | computational_priors; matrix_decision_panel |
| melanoidin_polymerization | HMW melanoidins scavenge sulfur volatiles via irreversible CROSSPY-radical thioether formation. | data/research_corpus/raw/Gemini_deep_research_16.md; src/extrusion.py | process_gap_registry; computational_priors |

Chemistry families tracked: 16
First-class families: amino_acid_sugar_core
Expansion candidates: lipid_oxidation_and_carbonylic_crosstalk, thiamine_fragmentation_support, nucleotide_and_ribose_support, glutathione_and_peptide_support, alternative_protein_matrix_scope, carbonyl_donor_hierarchy, off_note_and_maillard_suppression, carbohydrate_pyrolysis_and_caramelization, fermentation_pretreatment
Open-gap families: none
Recommended next family: lipid_oxidation_and_carbonylic_crosstalk
Policy: expand_product_scope_by_adding_benchmark_visible_family_lanes_not_by_diluting_the_core
Ingestion policy: beyond_amino_acid_sugar_use_the_same_runtime_contract_but_choose_family_specific_payload_types_benchmark_payloads_for_observable_panels_priors_for_transfer_and_gap_registries_for_non_closable_scope
