# PREREG A6 audit: the keyed pH control, final pH and base mixture (STATUS_AND_PLAN sec. 0, step 8)

*audit on sources already read; no reading of the lever table moves; PREREG A6 rules are scored only on sources keyed after the A6 commit. Spec: [STEP8_FIELDS_SPEC.md](../../docs/reference/lever_table/STEP8_FIELDS_SPEC.md).*

## Coverage

- pH contrasts: 623; `ph_control`: adjusted_unbuffered 261, buffered 348, unadjusted 9, unstated 5 (check 69, derived 237, keyed 317).
- final pH on every arm: 146 contrasts; final-pH spread of at least 0.5 pH unit over the arms the slope uses (A6-R2 slope reported): 135.
- categorical pairs: 2322; base mixture differs beyond the lever (A6-R3 drops them): 37; arms' base from check 498, derived 3302.
- reader-checked sources: Zhang2024, ames2001b, feng2022, grosch1992, hincelin1992, hofmann2000, hwang1995, martin2001, mottram2002, park2022mealworm, shibamoto1977, wei2019flax, xu2007, zhang2018soypep, zhou2023, zhou2023c.

## pH cells with two or more laboratories

| odorant | band | watched | labs | control (keyed) | main reading | R1 buffered | buffered set = free text's? | R1 unbuffered | R2 final pH (contrasts) | changes |
|---|---|---|---|---|---|---|---|---|---|---|
| (E)-2-Nonenal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 6 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (1) | none |
| (E)-2-Octenal | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 8, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| (E,E)-2,4-Decadienal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 6 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (1) | none |
| 1-Octen-3-ol | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 6, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 3 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (4) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2,3-Dimethylpyrazine | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 3 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| 2,5-Dimethylpyrazine | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 4 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2,5-Dimethylpyrazine | pH 6-7.5, initial or label |  | 2 | buffered 2 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | too few laboratories (1) | none |
| 2,5-Dimethylpyrazine | pH >7.5, initial or label | answered, fragile | 4 | adjusted_unbuffered 2, buffered 5 | answered by the literature | too few laboratories (2 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2,6-Dimethylpyrazine | pH 6-7.5, initial or label |  | 4 | adjusted_unbuffered 1, buffered 3 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | direction unclear (2) | R1 unbuffered band: too few laboratories |
| 2,6-Dimethylpyrazine | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| 2-Acetylthiazole | pH 4.5-6, initial or label | fragile | 4 | adjusted_unbuffered 2, buffered 1, unadjusted 1 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (3 labs) | direction unclear (2) | R1 buffered band: too few laboratories |
| 2-Acetylthiazole | pH 6-7.5, initial or label | answered | 5 | adjusted_unbuffered 2, buffered 6 | answered by the literature | answered by the literature (4 labs) | yes | too few laboratories (1 labs) | answered by the literature (4) | R1 unbuffered band: too few laboratories |
| 2-Ethyl-3,5-dimethylpyrazine | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (1) | none |
| 2-Ethyl-3,6-dimethylpyrazine | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 2 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 3 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (4) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Ethylpyrazine | pH >7.5, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 2 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Furanmethanol | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 3 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | direction unclear (3) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Furanmethanol | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (3) | none |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | fragile | 7 | adjusted_unbuffered 2, buffered 14 | direction unclear | answered by the literature (5 labs) | yes | too few laboratories (2 labs) | direction unclear (7) | R1 buffered band: answered by the literature; R1 unbuffered band: too few laboratories |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | answered | 5 | adjusted_unbuffered 3, buffered 6 | answered by the literature | answered by the literature (4 labs) | yes | too few laboratories (1 labs) | answered by the literature (5) | R1 unbuffered band: too few laboratories |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | sign changes inside one laboratory | 5 | adjusted_unbuffered 1, buffered 10 | direction unclear | direction unclear (4 labs) | yes | too few laboratories (1 labs) | direction unclear (4) | R1 unbuffered band: too few laboratories |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label |  | 4 | adjusted_unbuffered 2, buffered 4 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | direction unclear (3) | R1 unbuffered band: too few laboratories |
| 2-Methylbutanal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 8 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 2-Methylpropanal (isobutyraldehyde) | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 4 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| 2-Methylpyrazine | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 4 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Methylpyrazine | pH 6-7.5, initial or label |  | 9 | adjusted_unbuffered 3, buffered 7 | direction unclear | direction unclear (6 labs) | yes | direction unclear (3 labs) | direction unclear (3) | none |
| 2-Methylpyrazine | pH >7.5, initial or label |  | 5 | adjusted_unbuffered 3, buffered 4 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (2 labs) | - (0) | R1 unbuffered band: too few laboratories |
| 2-Pentylfuran | pH 4.5-6, initial or label | fragile | 4 | adjusted_unbuffered 17, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (3 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| 3-Methylbutanal | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 3-Methylbutanal | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Benzaldehyde | pH 4.5-6, initial or label | fragile | 4 | adjusted_unbuffered 17, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (3 labs) | - (0) | R1 buffered band: too few laboratories |
| Benzaldehyde | pH 6-7.5, initial or label | fragile | 4 | adjusted_unbuffered 1, buffered 2, unadjusted 1 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Bis(2-furfuryl) disulfide | pH 4.5-6, initial or label |  | 2 | buffered 2 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | too few laboratories (1) | none |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | conflict between laboratories | 4 | adjusted_unbuffered 1, buffered 4 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | direction unclear (3) | R1 unbuffered band: too few laboratories |
| DMHF | pH 6-7.5, initial or label |  | 3 | buffered 5 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (0 labs) | - (0) | R1 unbuffered band: too few laboratories |
| Furfural | pH 4.5-6, initial or label | fragile | 6 | adjusted_unbuffered 17, buffered 12 | direction unclear | direction unclear (4 labs) | yes | direction unclear (3 labs) | direction unclear (8) | none |
| Furfural | pH 6-7.5, initial or label | answered, fragile | 7 | adjusted_unbuffered 2, buffered 6 | answered by the literature | direction unclear (5 labs) | yes | too few laboratories (2 labs) | direction unclear (3) | R1 buffered band: direction unclear; R1 unbuffered band: too few laboratories; R2 final pH: direction unclear |
| Heptanal | pH 4.5-6, initial or label | fragile | 5 | adjusted_unbuffered 7, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (4 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Hexanal | pH 4.5-6, initial or label | fragile | 4 | adjusted_unbuffered 17, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (3 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Hexanal | pH 6-7.5, initial or label |  | 2 | buffered 3 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Methional (3-(methylthio)propanal) | pH 4.5-6, initial or label | conflict between laboratories | 2 | adjusted_unbuffered 2 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| Nonanal | pH 4.5-6, initial or label | fragile | 5 | adjusted_unbuffered 18, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | direction unclear (4 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Nonanal | pH 6-7.5, initial or label | conflict between laboratories | 4 | buffered 3, unadjusted 1 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | direction unclear (1) | R1 unbuffered band: too few laboratories |
| Octanal | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | direction unclear (1) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Octanal | pH 6-7.5, initial or label |  | 2 | buffered 1, unadjusted 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Pentanal | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 3, buffered 6 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Pentanal | pH 6-7.5, initial or label | fragile | 3 | buffered 2, unadjusted 1 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH 4.5-6, initial or label | fragile | 4 | adjusted_unbuffered 16, buffered 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH 6-7.5, initial or label |  | 4 | buffered 6, unstated 5 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (0 labs) | - (0) | R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH <4.5, initial or label |  | 2 | buffered 9 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Phenylacetaldehyde | pH >7.5, initial or label | sign changes inside one laboratory | 2 | adjusted_unbuffered 1, buffered 12 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Pyrazines | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 4 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Pyrazines | pH 6-7.5, initial or label | answered, fragile | 4 | adjusted_unbuffered 1, buffered 3 | answered by the literature | direction unclear (3 labs) | yes | too few laboratories (1 labs) | - (0) | R1 buffered band: direction unclear; R1 unbuffered band: too few laboratories |
| Pyrazines | pH >7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Tetramethylpyrazine | pH 6-7.5, final | sign changes inside one laboratory | 2 | adjusted_unbuffered 4 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | too few laboratories (3) | none |
| Tetramethylpyrazine | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Trimethylpyrazine | pH 6-7.5, final |  | 2 | adjusted_unbuffered 5 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | too few laboratories (5) | none |
| Trimethylpyrazine | pH 6-7.5, initial or label | conflict between laboratories | 5 | adjusted_unbuffered 2, buffered 4 | direction unclear | direction unclear (4 labs) | yes | too few laboratories (2 labs) | direction unclear (3) | R1 unbuffered band: too few laboratories |
| Unsubstituted pyrazine | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (1) | none |
| Unsubstituted pyrazine | pH >7.5, initial or label | fragile | 4 | adjusted_unbuffered 3, buffered 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |

### Contrasts that carry a final pH (initial range against the pH the pots ended at)

| odorant | band | source | laboratory | initial or label pH | final pH by arm | slope, keyed basis | slope on final pH |
|---|---|---|---|---|---|---|---|
| (E)-2-Nonenal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.67 | -0.80 |
| (E,E)-2,4-Decadienal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.33 | -0.38 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.06 | +0.23 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.39 | -0.45 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6.2, 6.8 | +0.09 | +0.17 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.55 | +1.65 |
| 2,5-Dimethylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.83 | -0.98 |
| 2,6-Dimethylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.10 | -0.12 |
| 2,6-Dimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.89 | +1.61 |
| 2-Acetylthiazole | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | +0.23 | +0.25 |
| 2-Acetylthiazole | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | +0.35 | +0.28 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.12 | +0.46 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 4.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | +0.79 | +0.82 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | +0.49 | - |
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.17 | +0.32 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.16 | +0.50 |
| 2-Ethyl-3,5-dimethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.06 | +0.24 |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.03 | +0.12 |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.32 | -0.38 |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.14 | +0.27 |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.57 | +1.72 |
| 2-Furanmethanol | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | +0.33 | +0.26 |
| 2-Furanmethanol | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.82, 6.38 | -0.10 | -0.08 |
| 2-Furanmethanol | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.72, 6.15 | +0.02 | +0.01 |
| 2-Furanmethanol | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 3.86, 6.05 | +0.00 | +0.00 |
| 2-Furanmethanol | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.60 | -1.14 |
| 2-Furanmethanol | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.19 | +0.58 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | +0.40 | +0.50 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 4.09, 6.42 | -0.71 | -0.61 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.7, 6.38 | -0.23 | -0.17 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | -0.37 | -0.29 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.82, 6.38 | +0.08 | +0.06 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.51, 6.11 | -0.50 | -0.38 |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.72, 6.15 | -0.73 | -0.60 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 3.86, 6.05 | -0.24 | -0.39 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 5.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.32 | -0.35 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.11 | - |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -1.09 | -1.97 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.46 | -1.34 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 8 | 3.2, 3.4, 5 | -0.11 | -0.22 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5 | 4.4, 4.9, 5.4, 5.8 | -0.65 | -0.65 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 4.09, 6.42 | -1.02 | -0.87 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | -0.16 | -0.13 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.72, 6.15 | +0.57 | +0.47 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 4.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.20 | -0.21 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.09 | - |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.96 | -1.81 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.27 | -0.80 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 7 | 3.2, 3.4, 5 | +0.36 | - |
| 2-Methylpyrazine | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.35 | - |
| 2-Methylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.02 | +0.08 |
| 2-Methylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | +0.07 | +0.08 |
| 2-Methylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.35 | +0.64 |
| 2-Pentylfuran | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.28 | -0.29 |
| Bis(2-furfuryl) disulfide | pH 4.5-6, initial or label | hou2017ph | btbu_beijing | 4.5 - 5.5 | 4.13, 5.08, 6.05, 6.9 | +0.28 | +0.30 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 5.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -1.00 | -1.07 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.71 | -1.33 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6, 6.3 | +0.09 | +0.24 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 7 | 3.2, 3.4, 5 | +0.24 | - |
| Furfural | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.64 | -0.67 |
| Furfural | pH 4.5-6, initial or label | mottram2002 | reading_mottram_parker | 4.2 - 5.6 | 5.5, 3.9 | -2.93 | -2.56 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 4.09, 6.42 | -0.30 | -0.25 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.7, 6.38 | -0.02 | -0.01 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | -0.03 | -0.02 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.82, 6.38 | -0.71 | -0.55 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.51, 6.11 | -0.28 | -0.22 |
| Furfural | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.72, 6.15 | -0.54 | -0.44 |
| Furfural | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 4.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.44 | -0.46 |
| Furfural | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.18 | - |
| Furfural | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.66 | -1.25 |
| Furfural | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.24 | -0.71 |
| Heptanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.40 | -0.45 |
| Hexanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.18 | -0.19 |
| Nonanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | +0.03 | +0.05 |
| Nonanal | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 4.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | +0.25 | +0.27 |
| Octanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.05 | -0.05 |
| Tetramethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 5.6 - 6.8 | 5.76, 6.62, 7.36 | +0.14 | +0.20 |
| Tetramethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 6.8 - 7.4 | 5.6, 6.85, 7.64 | -1.48 | -1.13 |
| Tetramethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 6.8 - 7.4 | 5.72, 6.98, 7.28 | +0.80 | - |
| Tetramethylpyrazine | pH 6-7.5, final | huang2010 | ohio_state_barringer | 6.43 - 7.87 | 5.9, 6.43, 6.98, 7.39, 7.87 | +0.43 | +0.43 |
| Trimethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 5.6 - 6.8 | 5.76, 6.62, 7.36 | -0.01 | -0.02 |
| Trimethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 5.6 - 6.8 | 5.6, 6.85, 7.64 | -0.10 | -0.09 |
| Trimethylpyrazine | pH 6-7.5, final | ames2001b | reading_mottram_parker | 5.6 - 7.4 | 5.72, 6.98, 7.28 | -0.23 | -0.23 |
| Trimethylpyrazine | pH 6-7.5, final | huang2010 | ohio_state_barringer | 6.52 - 8.01 | 5.33, 6.02, 6.52, 6.98, 7.6, 8.01 | +0.36 | +0.36 |
| Trimethylpyrazine | pH 6-7.5, final | huang2010 | ohio_state_barringer | 6.52 - 8.01 | 5.33, 6.02, 6.52, 6.98, 7.6, 8.01 | +0.28 | +0.28 |
| Trimethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.09 | +0.33 |
| Trimethylpyrazine | pH 6-7.5, initial or label | hou2017ph | btbu_beijing | 6.5 - 7.5 | 4.13, 5.08, 6.05, 6.9 | -0.29 | -0.34 |
| Trimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6.2, 6.8 | +1.02 | +2.04 |
| Unsubstituted pyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.03 | +0.13 |

## Categorical cells: base mixture (A6-R3)

| odorant | lever | watched | pairs | base differs | reference pot by laboratory | changes |
|---|---|---|---|---|---|---|
| 2,3-Dimethylpyrazine | lipid |  | 6 | 3 | reading_mottram_parker: none: amine+carbonyl, none: amine+carbonyl+sulfur; shou_lai_keqiang: none: amine+carbonyl | none |
| 2,5-Dimethylpyrazine | lipid |  | 13 | 4 | chonbuk_hwang: none: amine+carbonyl; reading_mottram_parker: none: amine+carbonyl, none: amine+carbonyl+sulfur; shou_lai_keqiang: none: amine+carbonyl | none |
| 2-Acetylfuran | sugar identity | sign changes inside one laboratory | 14 | 0 | astar_sifbi_chiang: xylose: amine; hfut_wei_zhaojun: fructose: amine+other+sulfur; kyunghee_lee: glucose: amine, glucose: nothing else | none |
| 2-Ethyl-3,6-dimethylpyrazine | lipid |  | 5 | 2 | reading_mottram_parker: none: amine+carbonyl, none: amine+carbonyl+sulfur; shou_lai_keqiang: none: amine+carbonyl | none |
| 2-Furanmethanol | amino acid identity | sign changes inside one laboratory | 45 | 0 | kyunghee_lee: glycine: carbonyl, lysine: carbonyl, none: carbonyl; reading_mottram_parker: glycine: carbonyl; wonkwang_cho_inhee: none: amine+carbonyl | none |
| 2-Furanmethanol | sugar identity | conflict between laboratories | 24 | 0 | astar_sifbi_chiang: glucose: amine; kyunghee_lee: glucose: amine, glucose: nothing else; nus_liu_shaoquan: none: amine | none |
| 2-Furfurylthiol (FFT) | lipid |  | 12 | 2 | chonbuk_hwang: none: amine+carbonyl; reading_mottram_parker: none: amine+carbonyl+sulfur | none |
| 2-Methyl-3-furanthiol (MFT) | lipid |  | 11 | 2 | chonbuk_hwang: none: amine+carbonyl; reading_mottram_parker: none: amine+carbonyl+sulfur | none |
| 2-Methyl-3-furanthiol (MFT) | sugar identity |  | 27 | 1 | dfa_garching: ribose: amine+sulfur, ribose: sulfur; hfut_wei_zhaojun: xylose: amine+other+sulfur; reading_mottram_parker: none: amine+sulfur, ribose 5-phosphate: amine+sulfur | none |
| 2-Methyl-3-furanthiol (MFT) | sulfur source | answered (level) | 18 | 0 | btbu_beijing: gsh + cysteine + thiamine: carbonyl, gsh + cysteine: carbonyl, none: amine+carbonyl+sulfur; dfa_garching: hydrogen sulfide: other, none: amine+carbonyl+sulfur, none: amine+sulfur, thiamine: amine+sulfur, thiamine: nothing else; ho_network: glutathione: amine+carbonyl, none: amine+carbonyl+sulfur; reading_mottram_parker: cysteine: other | none |
| 2-Methylpyrazine | amino acid identity | answered (level) | 143 | 0 | caas_bi_jinfeng: arginine: carbonyl+other; ghent_de_kimpe: free amino acids (lysine): carbonyl, gly-gly: carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl, hydrolysed protein: amine+carbonyl; ho_network: cys: amine+carbonyl+sulfur, glutathione: amine+carbonyl, none: amine+carbonyl; njau_zhao_liyan: methionine: amine+carbonyl; reading_mottram_parker: asparagine: nothing else, glycine: carbonyl, glycine: carbonyl+lipid, none: amine+carbonyl+lipid; wonkwang_cho_inhee: arginine: amine+carbonyl | none |
| 2-Methylpyrazine | lipid |  | 15 | 6 | chonbuk_hwang: none: amine+carbonyl; reading_mottram_parker: none: amine+carbonyl, none: amine+carbonyl+sulfur; shou_lai_keqiang: none: amine+carbonyl | none |
| Benzaldehyde | sugar identity | conflict between laboratories | 19 | 0 | astar_sifbi_chiang: glucose: amine; hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine; ohio_state_barringer: glucose: nothing else | none |
| Furfural | sugar identity | sign changes inside one laboratory | 25 | 0 | astar_sifbi_chiang: glucose: amine; kyunghee_lee: glucose: amine; nus_liu_shaoquan: none: amine, xylose: amine; ohio_state_barringer: glucose: nothing else; reading_mottram_parker: ribose 5-phosphate: amine+sulfur, ribose 5-phosphate: nothing else | none |
| Hexanal | sugar identity | conflict between laboratories | 20 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: arabinose: amine, none: amine, xylose: amine; ohio_state_barringer: glucose: nothing else | none |
| Methional (3-(methylthio)propanal) | sugar identity | answered (level) | 14 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine, xylose: amine; ohio_state_barringer: glucose: nothing else | none |
| Nonanal | lipid |  | 18 | 1 | caas_beijing: oleic acid: carbonyl+other+sulfur; chonbuk_hwang: c18:1n-9: amine+carbonyl, c18:1n-9: nothing else; huazhong_pan: oxidized lard: amine+carbonyl+sulfur; reading_mottram_parker: beef triglyceride: amine+carbonyl+other+sulfur, beef triglyceride: other; shou_lai_keqiang: none: amine+carbonyl | none |
| Phenylacetaldehyde | sugar identity | answered (level) | 9 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; nus_liu_shaoquan: none: amine; ohio_state_barringer: glucose: nothing else; quest_international: glucose: amine, xylose: amine | none |
| Pyrazines | amino acid identity | answered (level) | 138 | 0 | caas_bi_jinfeng: arginine: carbonyl+other; ghent_de_kimpe: free amino acids (lysine): carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl; ho_network: glutathione: amine+carbonyl, none: amine+carbonyl; minnesota_reineccius: lysine: carbonyl; reading_mottram_parker: asparagine: nothing else, none: amine+carbonyl+lipid; wonkwang_cho_inhee: none: amine+carbonyl | none |
| Pyrazines | lipid |  | 12 | 6 | chonbuk_hwang: none: amine+carbonyl; shou_lai_keqiang: none: amine+carbonyl | none |
| Pyrazines | sugar identity | sign changes inside one laboratory | 5 | 1 | hfut_wei_zhaojun: arabinose: amine+other+sulfur; minnesota_reineccius: glucose: amine; reading_mottram_parker: glucose: amine+carbonyl+lipid | none |
| Trimethylpyrazine | amino acid identity | answered (level) | 97 | 0 | caas_bi_jinfeng: arginine: carbonyl+other; ghent_de_kimpe: free amino acids (lysine): carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl; ho_network: none: amine+carbonyl; reading_mottram_parker: glycine: carbonyl, glycine: carbonyl+lipid; wonkwang_cho_inhee: none: amine+carbonyl | none |
| Trimethylpyrazine | lipid |  | 5 | 2 | reading_mottram_parker: none: amine+carbonyl, none: amine+carbonyl+sulfur; shou_lai_keqiang: none: amine+carbonyl | none |
| Unsubstituted pyrazine | lipid |  | 12 | 6 | chonbuk_hwang: none: amine+carbonyl; shou_lai_keqiang: none: amine+carbonyl | none |

R1 - R3 are PREREG amendment A6's frozen candidate rules, read here on contrasts already seen: a change is a finding about this table's sensitivity, not evidence for the rule. A reading needs three laboratories (PREREG sec. 5); 'changes' lists only readings that differ from the main table's.
