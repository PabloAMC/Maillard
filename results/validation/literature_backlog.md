# Literature Backlog

Queue policy: ready queues are exclusive to non-encoded intake rows; wet_lab_blocked comes only from structural gaps with closure_outcome=wet_lab_only

Encoded references: 205
Ready runtime: 0
Ready benchmark: 0
Wet-lab blocked: 4
Queue conflicts: 0

## Ready Runtime

| ID | Template | Family | Matrix | Triage | Encoding | Next Action |
| --- | --- | --- | --- | --- | --- | --- |
| none | n/a | n/a | n/a | n/a | n/a | no runtime-only backlog remains in the curated intake registry |

## Ready Benchmark

| ID | Template | Family | Matrix | Triage | Encoding | Next Action |
| --- | --- | --- | --- | --- | --- | --- |
| none | n/a | n/a | n/a | n/a | n/a | no benchmark-ready backlog remains in the curated intake registry |

## Wet-Lab Blocked

| Gap | Priority | Decision | Near Misses | Missing Contract |
| --- | --- | --- | --- | --- |
| ppi_meaty_positive_matrix_benchmark | critical | no_external_candidate_package | nishimura_abe_2024, pmc9905368_spi_hvp_xylose_benchmark | intact aqueous PPI matrix, exogenous ribose plus cysteine in the same run, absolute simultaneous MFT or FFT plus hexanal quantification |
| spi_meaty_positive_matrix_benchmark | critical | no_external_candidate_package | nishimura_abe_2024, trikusuma_2019, pmc9905368_spi_hvp_xylose_benchmark | intact aqueous SPI matrix, absolute MFT or FFT quantification with internal standards, same-run off-flavour markers such as hexanal or 2-pentylfuran |
| ppi_spi_time_series | medium |  | none | none |
| meaty_off_flavour_safety_tradeoff_panel | high | split_across_multiple_studies | pmc9905368_spi_hvp_xylose_benchmark, trikusuma_2019, squeo_2023 | same-run meaty sulfur markers, same-run adverse lipid markers, same-run safety endpoint coverage |

## Minimum Primary Experiment

Matrices: PPI 5% buffered slurry, SPI 5% buffered slurry
Exogenous precursors: D-ribose_mM=1.0, L-cysteine_mM=1.0
Analytical panel: 2-methyl-3-furanthiol, 2-furfurylthiol, hexanal, 2-pentylfuran, furfural, 2,5-dimethylpyrazine
Companion assays: Ellman / free SH, OPA / free amino, DSC or equivalent denaturation proxy
Instrumentation: HS-SPME-GC-MS with isotopically labelled internal standard; add acrylamide LC-MS/MS for the high-temperature SPI arm.

| Matrix | pH | Temp C | Time Points min |
| --- | ---: | ---: | --- |
| PPI | 5.5 | 95 | 0, 30, 60, 120, 240 |
| SPI | 5.8 | 120 | 0, 30, 60, 120, 240 |
