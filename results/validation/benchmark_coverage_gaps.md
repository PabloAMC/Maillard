# Benchmark Coverage Gaps

| Dimension | Category | Benchmarks | Status | Note |
| --- | --- | ---: | --- | --- |
| protein_system | free | 7 | covered | present in current benchmark set |
| protein_system | pea_iso | 4 | covered | present in current benchmark set |
| protein_system | soy_iso | 3 | covered | present in current benchmark set |
| protein_system | pea_conc | 0 | gap | no benchmark currently covers this protein system |
| protein_system | soy_conc | 0 | gap | no benchmark currently covers this protein system |
| protein_system | myco | 0 | gap | no benchmark currently covers this protein system |
| process_state | ambient_slurry | 2 | covered | present in current benchmark set |
| process_state | aqueous_pre_extrusion_model | 4 | covered | present in current benchmark set |
| process_state | heated_matrix | 2 | covered | present in current benchmark set |
| execution_path | free_precursor | 7 | covered | free-AA sulfur and safety envelope |
| execution_path | matrix_only | 3 | covered | matrix off-flavour intake/headspace envelope |
| execution_path | matrix_precursor_augmented | 4 | covered | matrix meaty-positive reproducibility harness envelope |
| scientific_gap | external_matrix_meaty_positive | 0 | gap | wet-lab quantitative meaty-positive matrix benchmarks are the main blocker for broad alt-protein validation |
| structural_gap | ppi_meaty_positive_matrix_benchmark | 0 | gap | No benchmark-eligible paper exists; the SLR confirms this is a structural literature gap rather than an encoding backlog. Near misses: nishimura_abe_2024, pmc9905368_spi_hvp_xylose_benchmark. |
| structural_gap | spi_meaty_positive_matrix_benchmark | 0 | gap | The SPI-side literature bifurcates into soy hydrolysate sulfur chemistry and separate off-flavour papers; no single benchmark-eligible package closes the mixed meaty-positive contract. Near misses: nishimura_abe_2024, trikusuma_2019, pmc9905368_spi_hvp_xylose_benchmark. |
| structural_gap | meaty_off_flavour_safety_tradeoff_panel | 0 | gap | The literature still splits meaty chemistry, off-flavour suppression, and safety endpoints across different studies, so there is no external benchmark package to promote honestly. Near misses: pmc9905368_spi_hvp_xylose_benchmark, trikusuma_2019, squeo_2023. |

Rows: 16
Gaps: 7
