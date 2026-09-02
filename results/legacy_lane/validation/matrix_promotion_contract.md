# Matrix Promotion Contract

Contract id: matrix_external_decision_ready_v1
Minimum quantitative closed targets: 2
Requires measured_volatiles: yes
Requires external quantitative origin: yes
Requires mixed or meaty-positive target profile: yes
Requires passing ranking contract: yes
Disallow internal-candidate support: yes
Disallow directional support: yes

## Rule Notes

- External decision readiness is a benchmark-level promotion state, not a generic matrix-family claim.
- Internal reproducibility candidates and transferred priors can strengthen triage, but they do not by themselves unlock promotion.
- Mechanistic refinement stays secondary until the observable audit says the remaining blocker is no longer external evidence or transfer dependence.

## Selected Promotion Target

- benchmark: pea_isolate_ribose_cysteine_100C_45min_Internal2026
- protein: pea_iso
- process_state: aqueous_pre_extrusion_model
- target_profile: mixed
- selection_policy: prefer_mixed_matrix_lanes_with_the_broadest_same_protein_external_anchor_span_before_mechanistic_escalation
- rationale: same_protein_external_anchor_count=2
- rationale: same_protein_external_process_states=2
- rationale: quantitative_closed=0
- rationale: internal_candidate=6
- rationale: internal_measured_candidate=0
- rationale: internal_reference_candidate=6

## Benchmark Assessment

| Benchmark | Protein | Process State | Target Profile | Promotion Ready | Blocker | Passed Requirements |
| --- | --- | --- | --- | --- | --- | ---: |
| pea_isolate_40C_PratapSingh2021 | pea_iso | ambient_slurry | adverse_only | no | benchmark lacks meaty-positive targets | 5/6 |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | pea_iso | aqueous_pre_extrusion_model | mixed | no | insufficient externally measured target closure; current comparator is internal reference-only | 2/6 |
| pea_isolate_uht_140C_Trikusuma2019 | pea_iso | heated_matrix | adverse_only | no | benchmark lacks meaty-positive targets | 5/6 |
| soy_isolate_40C_PratapSingh2021 | soy_iso | ambient_slurry | adverse_only | no | benchmark lacks meaty-positive targets | 5/6 |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | soy_iso | aqueous_pre_extrusion_model | mixed | no | insufficient externally measured target closure; current comparator is internal reference-only | 2/6 |

## Requirement Details

| Benchmark | Requirement | Passed | Detail |
| --- | --- | --- | --- |
| pea_isolate_40C_PratapSingh2021 | Target profile includes meaty-positive compounds | no | adverse_only |
| pea_isolate_40C_PratapSingh2021 | Ranking contract passes | yes | pass |
| pea_isolate_40C_PratapSingh2021 | Comparator signal is wet-lab measured_volatiles | yes | measured_volatiles |
| pea_isolate_40C_PratapSingh2021 | Source is externally quantitative | yes | external_quantitative |
| pea_isolate_40C_PratapSingh2021 | At least two compounds are quantitatively closed | yes | 2 |
| pea_isolate_40C_PratapSingh2021 | No internal-candidate or directional dependencies remain | yes | internal=0; directional=0 |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | Target profile includes meaty-positive compounds | yes | mixed |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | Ranking contract passes | yes | pass |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | Comparator signal is wet-lab measured_volatiles | no | reference_volatiles |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | Source is externally quantitative | no | internal_reference_only |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | At least two compounds are quantitatively closed | no | 0 |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | No internal-candidate or directional dependencies remain | no | internal=6; directional=0 |
| pea_isolate_uht_140C_Trikusuma2019 | Target profile includes meaty-positive compounds | no | adverse_only |
| pea_isolate_uht_140C_Trikusuma2019 | Ranking contract passes | yes | pass |
| pea_isolate_uht_140C_Trikusuma2019 | Comparator signal is wet-lab measured_volatiles | yes | measured_volatiles |
| pea_isolate_uht_140C_Trikusuma2019 | Source is externally quantitative | yes | external_quantitative |
| pea_isolate_uht_140C_Trikusuma2019 | At least two compounds are quantitatively closed | yes | 3 |
| pea_isolate_uht_140C_Trikusuma2019 | No internal-candidate or directional dependencies remain | yes | internal=0; directional=0 |
| soy_isolate_40C_PratapSingh2021 | Target profile includes meaty-positive compounds | no | adverse_only |
| soy_isolate_40C_PratapSingh2021 | Ranking contract passes | yes | pass |
| soy_isolate_40C_PratapSingh2021 | Comparator signal is wet-lab measured_volatiles | yes | measured_volatiles |
| soy_isolate_40C_PratapSingh2021 | Source is externally quantitative | yes | external_quantitative |
| soy_isolate_40C_PratapSingh2021 | At least two compounds are quantitatively closed | yes | 2 |
| soy_isolate_40C_PratapSingh2021 | No internal-candidate or directional dependencies remain | yes | internal=0; directional=0 |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | Target profile includes meaty-positive compounds | yes | mixed |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | Ranking contract passes | yes | pass |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | Comparator signal is wet-lab measured_volatiles | no | reference_volatiles |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | Source is externally quantitative | no | internal_reference_only |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | At least two compounds are quantitatively closed | no | 0 |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | No internal-candidate or directional dependencies remain | no | internal=6; directional=0 |
