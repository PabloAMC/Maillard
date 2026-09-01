# Matrix Family Coverage

| Matrix Family | Category | Runtime Posture | Evidence Surface | Importance | Next Best Action |
| --- | --- | --- | --- | --- | --- |
| free_precursor | chemistry_core | quantitative_core | strict_ready_benchmarks | critical | preserve_as_reference_core |
| pea_isolate | plant_protein_matrix | directional_matrix | matrix_benchmark_plus_calibration | critical | promote_primary_meaty_positive_pea_benchmark |
| soy_isolate | plant_protein_matrix | directional_matrix | matrix_benchmark_plus_calibration | critical | promote_primary_meaty_positive_soy_benchmark |
| soy_hydrolysate | plant_protein_hydrolysate | qualitative_intake_only | qualitative_intake | high | preserve_as_qualitative_design_anchor |
| mycoprotein | alternative_protein_matrix | directional_matrix | bounded_calibration_prior | high | land_first_family_specific_mycoprotein_measurement_package |
| extrusion_heavy_systems | process_regime | directional_matrix | exploratory_mode | critical | add_family_specific_process_intake_when_literature_allows |
| coconut_oil_co_matrix | lipid_rich_co_matrix | indirect_generic_support | generic_fat_fraction_plus_generic_lipid_chemistry | high | encode_coconut_oil_as_explicit_open_family_gap |
| other_plant_proteins | alternative_protein_matrix | open_gap | narrative_only_or_unencoded_literature | high | rank_next_family_by_product_impact_and_literature_closability |

## Expansion Gates

| Matrix Family | Support Class | Expansion Status | Primary Blocker | Artifacts |
| --- | --- | --- | --- | ---: |
| free_precursor | explicit_supported | reference_core | matrix accessibility effects | 2 |
| pea_isolate | explicit_supported | promote_primary_benchmark | external quantitative meaty-positive matrix benchmark | 3 |
| soy_isolate | explicit_supported | promote_primary_benchmark | external quantitative SPI meaty-positive benchmark | 3 |
| soy_hydrolysate | explicit_supported | hold_intake_only | absolute ppb validation | 1 |
| mycoprotein | explicit_supported | bounded_expansion_candidate | executable benchmark family | 2 |
| extrusion_heavy_systems | explicit_supported | hold_process_regime_only | benchmark-grade quantitative closure | 1 |
| coconut_oil_co_matrix | indirect_only | blocked_on_family_specific_evidence | coconut-specific lipid profile | 3 |
| other_plant_proteins | open_gap | blocked_on_runtime_prior_and_benchmark | runtime priors | 1 |

## Scope Notes

| Matrix Family | Supported Today | Not Supported Today |
| --- | --- | --- |
| free_precursor | quantitative free-precursor screening; reaction-family sensitivity and selective offline refinement; reaction-centered chemistry benchmarking | matrix accessibility effects; family-specific fat-phase partitioning; protein-bound precursor release |
| pea_isolate | matrix-only off-flavour intake and headspace benchmarks; pea accessibility and process-state priors; internal meaty-positive ranking harnesses | external quantitative meaty-positive matrix benchmark; direct retention benchmark for MFT or FFT; broad closure beyond pea-focused literature families |
| soy_isolate | matrix-only off-flavour intake benchmarks; soy accessibility and process-state priors; qualitative and internal meaty-positive harnesses | external quantitative SPI meaty-positive benchmark; benchmark-closed sulfur time series; broad closure beyond soy-focused literature families |
| soy_hydrolysate | qualitative soy sulfur-chemistry intake; design support for future benchmark conditions | absolute ppb validation; strict benchmark promotion; direct transfer to soy isolate quantitative claims |
| mycoprotein | first-class matrix family in computational priors; directional accessibility and process-state support; runtime-executable bounded reference windows | executable benchmark family; quantitative headspace closure; meaty-positive matrix benchmark |
| extrusion_heavy_systems | exploratory recommendation mode; shared accessibility and process-structured surrogates; explicit trust warnings | benchmark-grade quantitative closure; family-specific fat-phase calibration; strict decision-gate use |
| coconut_oil_co_matrix | generic fat-fraction partitioning in headspace and recommendation layers; generic lipid aldehyde formation and trapping chemistry; generic lipid-Maillard crosstalk at the chemistry-family level | coconut-specific lipid profile; coconut-specific oxidation calibration; coconut-specific retention or release benchmark; family-specific tradeoff benchmark combining sulfur positives and lipid off-notes |
| other_plant_proteins | high-level scientific awareness that these families matter | runtime priors; executable benchmark families; family-specific accessibility or observability support |

Matrix families tracked: 8
Explicitly supported families: free_precursor, pea_isolate, soy_isolate, soy_hydrolysate, mycoprotein, extrusion_heavy_systems
Indirect-only families: coconut_oil_co_matrix
Open-gap families: other_plant_proteins
Bounded expansion candidates: mycoprotein
Scope-hold families: soy_hydrolysate, extrusion_heavy_systems
Evidence-blocked families: coconut_oil_co_matrix, other_plant_proteins
Policy: matrix_family_scope_must_distinguish_explicit_support_from_generic_indirect_support
Expansion policy: do_not_broaden_matrix_scope_beyond_bounded_candidates_until_the_next_family_has_runtime_evidence_and_a_named_benchmark_or_calibration_landing
