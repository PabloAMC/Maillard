# Family Lane Validation

| SLR | Family | Posture | Benchmarks | Strict Ready | Supported | Payload Roles | Execution Paths |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| 01 | amino_acid_sugar_core | first_class_core | 7 | 0 | 7 | benchmark_intake, benchmark_payload, calibration_payload, directional_prior | free_precursor=3, matrix_precursor_augmented=4 |
| 02 | lipid_oxidation_and_carbonylic_crosstalk | immediate_expansion_lane | 7 | 0 | 7 | benchmark_intake, benchmark_payload, calibration_payload, directional_prior | matrix_only=3, matrix_precursor_augmented=4 |
| 03 | thiamine_fragmentation_support | high_value_support_lane | 2 | 0 | 2 | benchmark_intake, benchmark_payload | free_precursor=2 |
| 08 | off_note_and_maillard_suppression | guardrail_lane | 1 | 0 | 1 | benchmark_intake, directional_prior | free_precursor=1 |
| 09 | carbohydrate_pyrolysis_and_caramelization | failure_mode_lane | 5 | 0 | 5 | benchmark_intake, benchmark_payload, calibration_payload, directional_prior | free_precursor=1, matrix_precursor_augmented=4 |
| 12 | protein_damage_markers | first_class_runtime_lane | 3 | 0 | 3 | benchmark_intake, directional_prior | free_precursor=3 |

## Lane Summary

| Execution Path | Benchmarks | Strict Ready | Supported | Chemistry Families | Payload Roles |
| --- | ---: | ---: | ---: | --- | --- |
| free_precursor | 7 | 0 | 7 | amino_acid_sugar_core, carbohydrate_pyrolysis_and_caramelization, off_note_and_maillard_suppression, protein_damage_markers, thiamine_fragmentation_support | benchmark_intake, benchmark_payload, calibration_payload, directional_prior |
| matrix_only | 3 | 0 | 3 | lipid_oxidation_and_carbonylic_crosstalk | benchmark_intake, benchmark_payload |
| matrix_precursor_augmented | 4 | 0 | 4 | amino_acid_sugar_core, carbohydrate_pyrolysis_and_caramelization, lipid_oxidation_and_carbonylic_crosstalk | benchmark_payload, calibration_payload, directional_prior |

Benchmarks summarized: 14
Chemistry families summarized: 6
Execution lanes summarized: 3
