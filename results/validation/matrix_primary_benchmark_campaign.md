# Matrix Primary Benchmark Campaign

Protocol id: ppi_spi_ribose_cysteine_primary_benchmark_2026
Selected matrix: pea_iso
Selected benchmark id: pea_isolate_ribose_cysteine_100C_45min_Internal2026
Fallback matrix: soy_iso
Primary data blocker: external_quantitative_measured_volatiles_for_mixed_matrix_lane
Next best action: run_primary_benchmark_protocol_and_land_results_as_benchmark_json_plus_process_state_calibration

## Campaign Arms

| Matrix | Benchmark | Temp C | pH | Time Points | Calibration Route | Promotion Blocker | Would Close | Remaining After Protocol |
| --- | --- | ---: | ---: | --- | --- | --- | --- | --- |
| pea_iso | pea_isolate_ribose_cysteine_100C_45min_Internal2026 | 95.0 | 5.5 | 0, 30, 60, 120, 240 | no_internal_comparator | insufficient externally measured target closure; current comparator is internal reference-only | comparator_is_measured_volatiles, external_quantitative_origin, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |
| soy_iso | soy_isolate_ribose_cysteine_100C_45min_Internal2026 | 120.0 | 5.8 | 0, 30, 60, 120, 240 | no_internal_comparator | insufficient externally measured target closure; current comparator is internal reference-only | comparator_is_measured_volatiles, external_quantitative_origin, minimum_quantitative_closed_targets | no_internal_or_directional_dependencies |

## Required Panel

| Matrix | Transfer-Ready Targets | Evidence/Calibration Blockers | Mechanistic Blockers | Hexanal Ratio | Nonanal Ratio | Companion Assays | Replicates |
| --- | --- | --- | --- | ---: | ---: | --- | ---: |
| pea_iso | 2-furfurylthiol, 2-methyl-3-furanthiol, bis(2-methyl-3-furyl) disulfide, 2,5-dimethylpyrazine | Hexanal, Nonanal | none | n/a | n/a | Ellman free SH, OPA free amino groups, DSC or equivalent denaturation proxy, post-heating pH measurement | 3 |
| soy_iso | 2-furfurylthiol, 2-methyl-3-furanthiol, bis(2-methyl-3-furyl) disulfide, 2,5-dimethylpyrazine | Hexanal, Nonanal | none | n/a | n/a | Ellman free SH, OPA free amino groups, DSC or equivalent denaturation proxy, post-heating pH measurement | 3 |

## Promotion Delta

| Matrix | Requirement | Passed Today | Detail |
| --- | --- | --- | --- |
| pea_iso | Comparator signal is wet-lab measured_volatiles | False | reference_volatiles |
| pea_iso | Source is externally quantitative | False | internal_reference_only |
| pea_iso | At least two compounds are quantitatively closed | False | 0 |
| pea_iso | No internal-candidate or directional dependencies remain | False | internal=6; directional=0 |
| soy_iso | Comparator signal is wet-lab measured_volatiles | False | reference_volatiles |
| soy_iso | Source is externally quantitative | False | internal_reference_only |
| soy_iso | At least two compounds are quantitatively closed | False | 0 |
| soy_iso | No internal-candidate or directional dependencies remain | False | internal=6; directional=0 |
