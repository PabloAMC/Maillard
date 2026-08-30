# Matrix Benchmark Assertions

| Benchmark | Protein | Path | Process State | Target Profile | Ranking Contract | Coverage | Min Coverage | Top-k | Top-k Hits | Top-k Status | Adverse Order | Max Ratio | Ratio Tol. | Ratio Status | Overall | Strict Gate Blocked | Blocker |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |
| pea_isolate_40C_PratapSingh2021 | pea_iso | matrix_only | ambient_slurry | adverse_only | pass | 1.000 | 1.000 | 2 | 2 | pass | pass | 1.011 | 2.000 | pass | pass | yes | benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present |
| pea_isolate_ribose_cysteine_100C_45min_Internal2026 | pea_iso | matrix_precursor_augmented | aqueous_pre_extrusion_model | mixed | pass | 1.000 | 1.000 | 6 | 6 | pass | fail | 1.000 | 2.500 | pass | fail | yes | missing external quantitative matrix evidence for meaty-positive targets; current comparator is internal reference-only |
| pea_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 | pea_iso | matrix_precursor_augmented | aqueous_pre_extrusion_model | mixed | pass | 1.000 | 1.000 | 6 | 6 | pass | fail | 1.000 | 2.000 | pass | fail | yes | missing external quantitative matrix evidence for meaty-positive targets; current comparator is synthetic model output (diagnostic only) |
| pea_isolate_uht_140C_Trikusuma2019 | pea_iso | matrix_only | heated_matrix | adverse_only | pass | 1.000 | 1.000 | 3 | 3 | pass | pass | 2.273 | 2.000 | fail | fail | yes | benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present |
| soy_isolate_40C_PratapSingh2021 | soy_iso | matrix_only | ambient_slurry | adverse_only | pass | 1.000 | 1.000 | 2 | 2 | pass | pass | 1.011 | 2.000 | pass | pass | yes | benchmark only anchors adverse/off-flavour markers; no external meaty-positive targets are present |
| soy_isolate_ribose_cysteine_100C_45min_Internal2026 | soy_iso | matrix_precursor_augmented | aqueous_pre_extrusion_model | mixed | pass | 1.000 | 1.000 | 6 | 6 | pass | pass | 1.000 | 2.500 | pass | pass | yes | missing external quantitative matrix evidence for meaty-positive targets; current comparator is internal reference-only |
| soy_isolate_ribose_cysteine_100C_45min_ProtocolPilot2026 | soy_iso | matrix_precursor_augmented | aqueous_pre_extrusion_model | mixed | pass | 1.000 | 1.000 | 6 | 6 | pass | pass | 1.000 | 2.000 | pass | pass | yes | missing external quantitative matrix evidence for meaty-positive targets; current comparator is synthetic model output (diagnostic only) |

Benchmarks asserted: 7
Assertion passes: 4
Strict gate remains blocked for all matrix benchmarks by contract until external evidence exists.
