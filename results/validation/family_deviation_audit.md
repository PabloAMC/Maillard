# Family Deviation Audit

This artifact identifies which benchmark points dominate family-level error tails and where targeted closure work should start.

- Quantitative points analyzed: 39
- Families tracked: 16
- Families with quantitative points: 4
- High-ratio threshold: 5.00
- High-ratio points: 3
- High-log10-error threshold: 0.50
- High-log10-error points: 4
- Experimental points (measured-volatiles only): 24
- Experimental high-ratio points: 3

## Root-Cause Distribution

- internal_synthetic_reference_match: 14
- model_or_mapping_mismatch: 4
- non_experimental_reference_comparator: 1
- within_expected_error: 20

## Family Deviation Table

| SLR | Family | Benchmarks | Quant Points | Median Ratio | P90 Ratio | Max Ratio | Mean |log10 error| | Trimmed Mean |log10 error| | High Ratio Pts |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 | Amino acid-sugar Maillard core | 7 | 20 | 1.000 | 4.380 | 6730.854 | 0.302 | 0.051 | 2 |
| 08 | Plant off-notes and Maillard suppression | 1 | 1 | 15.387 | 15.387 | 15.387 | 1.187 | 1.187 | 1 |
| 09 | Carbohydrate pyrolysis and caramelization | 5 | 3 | 1.000 | 5.457 | 5.457 | 0.246 | 0.246 | 1 |
| 02 | Lipid oxidation and carbonylic crosstalk | 7 | 15 | 1.000 | 1.011 | 2.273 | 0.024 | 0.001 | 0 |
| 03 | Thiamine degradation and sulfur support | 2 | 0 | - | - | - | - | - | 0 |
| 04 | Nucleotide degradation and ribose support | 0 | 0 | - | - | - | - | - | 0 |
| 05 | Glutathione and peptide support | 0 | 0 | - | - | - | - | - | 0 |
| 06 | Alternative protein matrix scope | 0 | 0 | - | - | - | - | - | 0 |
| 07 | Reducing sugar and carbonyl donor hierarchy | 0 | 0 | - | - | - | - | - | 0 |
| 10 | Microbial fermentation pretreatment | 0 | 0 | - | - | - | - | - | 0 |
| 11 | Maillard/Lipid Crosstalk | 0 | 0 | - | - | - | - | - | 0 |
| 12 | Protein Damage Markers | 3 | 0 | - | - | - | - | - | 0 |
| 13 | Polyphenol-Amino Capping | 0 | 0 | - | - | - | - | - | 0 |
| 14 | Ascorbic Acid Maillard | 0 | 0 | - | - | - | - | - | 0 |
| 15 | PE Stealth Sugar Sink | 0 | 0 | - | - | - | - | - | 0 |
| 16 | Melanoidin Polymerization | 0 | 0 | - | - | - | - | - | 0 |

## Fix Queue

| Ticket | Benchmark | Root Cause | Action | Max Ratio | Compounds |
| --- | --- | --- | --- | ---: | --- |
| FD-001 | Thiamine + cysteine + glucose, 120 C (Bolton et al., 1994) | model_or_mapping_mismatch | triage_mapping_projection_calibration | 6730.854 | 2-Methyl-3-furanthiol (MFT) |
| FD-002 | SPI extrusion, acrylamide benchmark, 130 C (ACS Ref. 3) | model_or_mapping_mismatch | triage_mapping_projection_calibration | 15.387 | acrylamide |
| FD-003 | PBMA vs beef comparator, 150 C (Resconi et al., 2023) | model_or_mapping_mismatch | triage_mapping_projection_calibration | 5.457 | furfural |

## Top Outliers by Ratio

| Benchmark | Family | Compound | Lane | Ratio | |log10 error| |
| --- | --- | --- | --- | ---: | ---: |
| Thiamine + cysteine + glucose, 120 C (Bolton et al., 1994) | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 6730.854 | 3.828 |
| Thiamine + cysteine + xylose, 145 C (Cerny, 2008) [reference anchor] | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 25.741 | 1.411 |
| SPI extrusion, acrylamide benchmark, 130 C (ACS Ref. 3) | Plant off-notes and Maillard suppression | acrylamide | free_precursor | 15.387 | 1.187 |
| PBMA vs beef comparator, 150 C (Resconi et al., 2023) | Carbohydrate pyrolysis and caramelization | furfural | free_precursor | 5.457 | 0.737 |
| Cysteine + ribose, 140 C (Hofmann, 1998) | Amino acid-sugar Maillard core | 2-methyl-3-furanthiol | free_precursor | 4.380 | 0.641 |
| Pea isolate UHT, 140 C (Trikusuma, 2020) | Lipid oxidation and carbonylic crosstalk | nonanal | matrix_only | 2.273 | 0.357 |
| Cysteine + ribose, 140 C (Hofmann, 1998) | Amino acid-sugar Maillard core | 2-furfurylthiol | free_precursor | 1.468 | 0.167 |
| Soy isolate, 40 C (Pratap Singh, 2021) | Lipid oxidation and carbonylic crosstalk | hexanal | matrix_only | 1.011 | 0.005 |
| Pea isolate, 40 C (Pratap Singh, 2021) | Lipid oxidation and carbonylic crosstalk | hexanal | matrix_only | 1.011 | 0.005 |
| Pea isolate, 40 C (Pratap Singh, 2021) | Lipid oxidation and carbonylic crosstalk | 2-pentylfuran | matrix_only | 1.000 | 0.000 |
| Soy isolate, 40 C (Pratap Singh, 2021) | Lipid oxidation and carbonylic crosstalk | 2-pentylfuran | matrix_only | 1.000 | 0.000 |
| Pea isolate UHT, 140 C (Trikusuma, 2020) | Lipid oxidation and carbonylic crosstalk | hexanal | matrix_only | 1.000 | 0.000 |
| Pea isolate UHT, 140 C (Trikusuma, 2020) | Lipid oxidation and carbonylic crosstalk | 2-pentylfuran | matrix_only | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | 2,5-dimethylpyrazine | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | 2-furfurylthiol | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | 2-methyl-3-furanthiol | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Lipid oxidation and carbonylic crosstalk | Hexanal | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Lipid oxidation and carbonylic crosstalk | Nonanal | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | bis(2-methyl-3-furyl) disulfide | matrix_precursor_augmented | 1.000 | 0.000 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Carbohydrate pyrolysis and caramelization | furfural | matrix_precursor_augmented | 1.000 | 0.000 |
