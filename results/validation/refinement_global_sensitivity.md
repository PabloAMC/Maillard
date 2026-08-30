# P3 Global Sensitivity

## Barrier Families

| Reaction Family | Dominant Direction | Max Impact | Affected Benchmarks |
| --- | --- | ---: | --- |
| amadori_rearrangement | down | 17.76 | cys_ribose_140C_Hofmann1998, hofmann1998_glucose_cysteine_145C_20min_pH5, hofmann1998_ribose_cysteine_145C_20min_pH5, pea_isolate_ribose_cysteine_100C_45min_Internal2026 |
| 1,2-enolisation | up | 8.02 | cys_ribose_140C_Hofmann1998, hofmann1998_fructose_cysteine_145C_20min_pH5, hofmann1998_glucose_cysteine_145C_20min_pH5, hofmann1998_ribose_cysteine_145C_20min_pH5 |
| thiol_oxidation | up | 6.37 | cys_ribose_140C_Hofmann1998, hofmann1998_c2c3_recombination_145C_20min_pH3, hofmann1998_c2c3_recombination_145C_20min_pH5, hofmann1998_c2c3_recombination_145C_20min_pH7 |
| aminoketone_condensation | down | 4.08 | pea_isolate_ribose_cysteine_100C_45min_Internal2026, pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026, resconi_2023_pbma_beef_identity_benchmark, soy_isolate_ribose_cysteine_100C_45min_Internal2026 |
| thiol_addition | down | 2.02 | cys_ribose_140C_Hofmann1998, hofmann1998_fructose_cysteine_145C_20min_pH5, hofmann1998_glucose_cysteine_145C_20min_pH5, hofmann1998_ribose_cysteine_145C_20min_pH5 |
| retro_aldol | down | 0.04 | none |
| schiff_condensation | up | 0.00 | none |
| strecker_degradation | up | 0.00 | none |
| cysteine_thermolysis | none | 0.00 | none |

## Process Axes

| Axis | Direction | Magnitude | Weighted Drift | Touched Benchmarks |
| --- | --- | ---: | ---: | ---: |
| time_minutes | up | 15.00 | 34.05 | 18 |
| temperature_celsius | up | 10.00 | 24.22 | 20 |
| temperature_celsius | down | -10.00 | 8.83 | 20 |
| time_minutes | down | -15.00 | 8.76 | 18 |
| pH | up | 0.40 | 6.17 | 12 |
| pH | down | -0.40 | 5.14 | 12 |
| water_activity | down | -0.05 | 1.00 | 9 |
| water_activity | up | 0.05 | 0.47 | 4 |

## Formulation Axes

| Bucket | Direction | Factor | Weighted Drift | Touched Benchmarks |
| --- | --- | ---: | ---: | ---: |
| amino_acids | down | 0.75 | 1.98 | 15 |
| amino_acids | up | 1.25 | 1.89 | 15 |
| sugars | down | 0.75 | 1.05 | 9 |
| sugars | up | 1.25 | 0.91 | 9 |

Benchmarks evaluated: 20
Top barrier family: amadori_rearrangement
Top process axis: time_minutes
Top formulation axis: amino_acids
