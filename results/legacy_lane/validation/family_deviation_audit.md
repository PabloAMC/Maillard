# Family Deviation Audit

This artifact identifies which benchmark points dominate family-level error tails and where targeted closure work should start.

- Quantitative points analyzed: 39
- Families tracked: 16
- Families with quantitative points: 4
- High-ratio threshold: 5.00
- High-ratio points: 8
- High-log10-error threshold: 0.50
- High-log10-error points: 11
- Experimental points (measured-volatiles only): 24
- Experimental high-ratio points: 8

## Root-Cause Distribution

- internal_synthetic_reference_match: 14
- model_or_mapping_mismatch: 11
- non_experimental_reference_comparator: 1
- within_expected_error: 13

## Family Deviation Table

| SLR | Family | Benchmarks | Quant Points | Median Ratio | P90 Ratio | Max Ratio | Mean |log10 error| | Trimmed Mean |log10 error| | High Ratio Pts |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 01 | Amino acid-sugar Maillard core | 14 | 24 | 2.806 | 24.037 | 6730.854 | 0.678 | 0.539 | 7 |
| 08 | Plant off-notes and Maillard suppression | 1 | 1 | 15.387 | 15.387 | 15.387 | 1.187 | 1.187 | 1 |
| 09 | Carbohydrate pyrolysis and caramelization | 3 | 3 | 1.268 | 5.457 | 5.457 | 0.314 | 0.314 | 1 |
| 02 | Lipid oxidation and carbonylic crosstalk | 5 | 11 | 1.000 | 1.011 | 2.273 | 0.033 | 0.001 | 0 |
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
| FD-002 | hofmann1998 glucose cysteine 145C 20min pH5 | model_or_mapping_mismatch | triage_mapping_projection_calibration | 29.584 | 2-Furfurylthiol (FFT) |
| FD-003 | SPI extrusion, acrylamide benchmark, 130 C (ACS Ref. 3) | model_or_mapping_mismatch | triage_mapping_projection_calibration | 15.387 | acrylamide |
| FD-004 | hofmann1998 fructose cysteine 145C 20min pH5 | model_or_mapping_mismatch | triage_mapping_projection_calibration | 14.458 | 2-Furfurylthiol (FFT), 2-Methyl-3-furanthiol (MFT) |
| FD-005 | hofmann1998 ribose cysteine 145C 20min pH5 | model_or_mapping_mismatch | triage_mapping_projection_calibration | 11.808 | 2-Furfurylthiol (FFT) |
| FD-006 | hofmann1998 c2c3 recombination 145C 20min pH3 | model_or_mapping_mismatch | triage_mapping_projection_calibration | 7.033 | 2-Methyl-3-furanthiol (MFT) |
| FD-007 | PBMA vs beef comparator, 150 C (Resconi et al., 2023) | model_or_mapping_mismatch | triage_mapping_projection_calibration | 5.457 | furfural |

## Top Outliers by Ratio

| Benchmark | Family | Compound | Lane | Ratio | |log10 error| |
| --- | --- | --- | --- | ---: | ---: |
| Thiamine + cysteine + glucose, 120 C (Bolton et al., 1994) | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 6730.854 | 3.828 |
| hofmann1998 glucose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Furfurylthiol (FFT) | free_precursor | 29.584 | 1.471 |
| Thiamine + cysteine + xylose, 145 C (Cerny, 2008) [reference anchor] | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 24.037 | 1.381 |
| SPI extrusion, acrylamide benchmark, 130 C (ACS Ref. 3) | Plant off-notes and Maillard suppression | acrylamide | free_precursor | 15.387 | 1.187 |
| hofmann1998 fructose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 14.458 | 1.160 |
| hofmann1998 ribose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Furfurylthiol (FFT) | free_precursor | 11.808 | 1.072 |
| hofmann1998 c2c3 recombination 145C 20min pH3 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 7.033 | 0.847 |
| hofmann1998 fructose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Furfurylthiol (FFT) | free_precursor | 5.939 | 0.774 |
| PBMA vs beef comparator, 150 C (Resconi et al., 2023) | Carbohydrate pyrolysis and caramelization | furfural | free_precursor | 5.457 | 0.737 |
| hofmann1998 glucose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 4.749 | 0.677 |
| hofmann1998 furan2aldehyde h2s 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Furfurylthiol (FFT) | free_precursor | 4.564 | 0.659 |
| hofmann1998 ribose cysteine 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 3.605 | 0.557 |
| Cysteine + ribose, 140 C (Hofmann, 1998) | Amino acid-sugar Maillard core | 2-methyl-3-furanthiol | free_precursor | 2.872 | 0.458 |
| hofmann1998 c2c3 recombination 145C 20min pH7 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 2.857 | 0.456 |
| hofmann1998 norfuraneol h2s 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 2.755 | 0.440 |
| Soy isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | 2-methyl-3-furanthiol | matrix_precursor_augmented | 2.470 | 0.393 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | 2-methyl-3-furanthiol | matrix_precursor_augmented | 2.470 | 0.393 |
| hofmann1998 c2c3 recombination 145C 20min pH5 | Amino acid-sugar Maillard core | 2-Methyl-3-furanthiol (MFT) | free_precursor | 2.459 | 0.391 |
| Pea isolate UHT, 140 C (Trikusuma, 2020) | Lipid oxidation and carbonylic crosstalk | nonanal | matrix_only | 2.273 | 0.357 |
| Pea isolate + ribose + cysteine, 100 C (Internal, 2026) | Amino acid-sugar Maillard core | bis(2-methyl-3-furyl) disulfide | matrix_precursor_augmented | 1.823 | 0.261 |
