# PREREG A6 audit: the keyed pH control, final pH and base mixture (STATUS_AND_PLAN sec. 0, step 8)

*audit on sources already read; no reading of the lever table moves; PREREG A6 rules are scored only on sources keyed after the A6 commit. Spec: [STEP8_FIELDS_SPEC.md](../../docs/reference/lever_table/STEP8_FIELDS_SPEC.md).*

## Coverage

- pH contrasts: 388; `ph_control`: adjusted_unbuffered 205, buffered 178, unstated 5 (check 69, derived 237, keyed 82).
- final pH on every arm: 109 contrasts; final-pH spread of at least 0.5 pH unit over the arms the slope uses (A6-R2 slope reported): 100.
- categorical pairs: 1998; base mixture differs beyond the lever (A6-R3 drops them): 7; arms' base from check 498, derived 2735.
- reader-checked sources: Zhang2024, ames2001b, feng2022, grosch1992, hincelin1992, hofmann2000, hwang1995, martin2001, mottram2002, park2022mealworm, shibamoto1977, wei2019flax, xu2007, zhang2018soypep, zhou2023, zhou2023c.

## pH cells with two or more laboratories

| odorant | band | watched | labs | control (keyed) | main reading | R1 buffered | buffered set = free text's? | R1 unbuffered | R2 final pH (contrasts) | changes |
|---|---|---|---|---|---|---|---|---|---|---|
| (E)-2-Octenal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 8 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 1-Octen-3-ol | pH 4.5-6, initial or label | sign changes inside one laboratory | 2 | adjusted_unbuffered 6 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (3) | none |
| 2-Acetylthiazole | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 1 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | answered by the literature (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories; R2 final pH: answered by the literature |
| 2-Acetylthiazole | pH 6-7.5, initial or label | answered | 4 | adjusted_unbuffered 2, buffered 5 | answered by the literature | answered by the literature (3 labs) | yes | too few laboratories (1 labs) | answered by the literature (3) | R1 unbuffered band: too few laboratories |
| 2-Ethylpyrazine | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (3) | none |
| 2-Furanmethanol | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 3 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (3) | none |
| 2-Furanmethanol | pH 6-7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (3) | none |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label | conflict between laboratories | 6 | adjusted_unbuffered 2, buffered 10 | direction unclear | direction unclear (4 labs) | yes | too few laboratories (2 labs) | direction unclear (7) | R1 unbuffered band: too few laboratories |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | answered | 4 | adjusted_unbuffered 3, buffered 5 | answered by the literature | answered by the literature (3 labs) | yes | too few laboratories (1 labs) | answered by the literature (4) | R1 unbuffered band: too few laboratories |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | sign changes inside one laboratory | 4 | adjusted_unbuffered 1, buffered 6 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | direction unclear (4) | R1 unbuffered band: too few laboratories |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 3 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Methylbutanal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 8 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 2-Methylpyrazine | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 2 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| 2-Methylpyrazine | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 2 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Methylpyrazine | pH >7.5, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 2 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Pentylfuran | pH 4.5-6, initial or label | conflict between laboratories | 3 | adjusted_unbuffered 17 | direction unclear | too few laboratories (0 labs) | yes | direction unclear (3 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| 3-Methylbutanal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Benzaldehyde | pH 4.5-6, initial or label | sign changes inside one laboratory | 3 | adjusted_unbuffered 17 | direction unclear | too few laboratories (0 labs) | yes | direction unclear (3 labs) | - (0) | R1 buffered band: too few laboratories |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 3 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| DMHF | pH 6-7.5, initial or label |  | 2 | buffered 4 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Furfural | pH 4.5-6, initial or label | fragile | 5 | adjusted_unbuffered 17, buffered 8 | direction unclear | answered by the literature (3 labs) | yes | direction unclear (3 labs) | direction unclear (8) | R1 buffered band: answered by the literature |
| Furfural | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 1, buffered 3 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | direction unclear (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Heptanal | pH 4.5-6, initial or label |  | 4 | adjusted_unbuffered 7 | direction unclear | too few laboratories (0 labs) | yes | direction unclear (4 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Hexanal | pH 4.5-6, initial or label |  | 3 | adjusted_unbuffered 17 | direction unclear | too few laboratories (0 labs) | yes | direction unclear (3 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Nonanal | pH 4.5-6, initial or label | sign changes inside one laboratory | 4 | adjusted_unbuffered 18 | direction unclear | too few laboratories (0 labs) | yes | direction unclear (4 labs) | direction unclear (1) | R1 buffered band: too few laboratories |
| Octanal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 2 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | too few laboratories (1) | none |
| Pentanal | pH 4.5-6, initial or label |  | 2 | adjusted_unbuffered 3 | too few laboratories | too few laboratories (0 labs) | yes | too few laboratories (2 labs) | - (0) | none |
| Phenylacetaldehyde | pH 4.5-6, initial or label | fragile | 3 | adjusted_unbuffered 16, buffered 1 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH 6-7.5, initial or label | fragile | 3 | buffered 5, unstated 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH <4.5, initial or label |  | 2 | buffered 9 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Pyrazines | pH 6-7.5, initial or label |  | 2 | buffered 2 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Pyrazines | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Trimethylpyrazine | pH 6-7.5, initial or label | fragile | 3 | adjusted_unbuffered 2, buffered 1 | direction unclear | too few laboratories (1 labs) | yes | too few laboratories (2 labs) | direction unclear (2) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Unsubstituted pyrazine | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 2, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |

### Contrasts that carry a final pH (initial range against the pH the pots ended at)

| odorant | band | source | laboratory | initial or label pH | final pH by arm | slope, keyed basis | slope on final pH |
|---|---|---|---|---|---|---|---|
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.06 | +0.23 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6.2, 6.8 | +0.09 | +0.17 |
| 2,3-Dimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.55 | +1.65 |
| 2-Acetylthiazole | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | +0.23 | +0.25 |
| 2-Acetylthiazole | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | +0.35 | +0.28 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.12 | +0.46 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | +0.49 | - |
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.17 | +0.32 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.16 | +0.50 |
| 2-Ethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.03 | +0.12 |
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
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.11 | - |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -1.09 | -1.97 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.46 | -1.34 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 8 | 3.2, 3.4, 5 | -0.11 | -0.22 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5 | 4.4, 4.9, 5.4, 5.8 | -0.65 | -0.65 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 4.09, 6.42 | -1.02 | -0.87 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.85, 6.37 | -0.16 | -0.13 |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | yang2011 | btbu_beijing | 4.7 - 6.7 | 3.72, 6.15 | +0.57 | +0.47 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.09 | - |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.96 | -1.81 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.27 | -0.80 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 7 | 3.2, 3.4, 5 | +0.36 | - |
| 2-Methylpyrazine | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.35 | - |
| 2-Methylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.02 | +0.08 |
| 2-Methylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.35 | +0.64 |
| 2-Pentylfuran | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.28 | -0.29 |
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
| Furfural | pH 6-7.5, initial or label | kang2026 | ho_network | 5.5 - 8 | 5.1, 4.9, 4.7 | -0.18 | - |
| Furfural | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.66 | -1.25 |
| Furfural | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.24 | -0.71 |
| Heptanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.40 | -0.45 |
| Hexanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.18 | -0.19 |
| Nonanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | +0.03 | +0.05 |
| Octanal | pH 4.5-6, initial or label | madruga1995 | reading_mottram_parker | 4.5 - 5.6 | 4.4, 4.9, 5.4, 5.8 | -0.05 | -0.05 |
| Trimethylpyrazine | pH 6-7.5, initial or label | chen2000 | ho_network | 5 - 8.5 | 7.43, 8.35 | +0.09 | +0.33 |
| Trimethylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6.2, 6.8 | +1.02 | +2.04 |

## Categorical cells: base mixture (A6-R3)

| odorant | lever | watched | pairs | base differs | reference pot by laboratory | changes |
|---|---|---|---|---|---|---|
| 2-Methyl-3-furanthiol (MFT) | sugar identity |  | 27 | 1 | dfa_garching: ribose: amine+sulfur, ribose: sulfur; hfut_wei_zhaojun: xylose: amine+other+sulfur; reading_mottram_parker: none: amine+sulfur, ribose 5-phosphate: amine+sulfur | none |
| 2-Methyl-3-furanthiol (MFT) | sulfur source | answered (level) | 16 | 0 | btbu_beijing: none: amine+carbonyl+sulfur; dfa_garching: hydrogen sulfide: other, none: amine+carbonyl+sulfur, none: amine+sulfur, thiamine: amine+sulfur, thiamine: nothing else; ho_network: glutathione: amine+carbonyl, none: amine+carbonyl+sulfur; reading_mottram_parker: cysteine: other | none |
| 2-Methylpyrazine | amino acid identity | fragile | 141 | 0 | ghent_de_kimpe: free amino acids (lysine): carbonyl, gly-gly: carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl, hydrolysed protein: amine+carbonyl; ho_network: cys: amine+carbonyl+sulfur, glutathione: amine+carbonyl, none: amine+carbonyl; njau_zhao_liyan: methionine: amine+carbonyl; reading_mottram_parker: asparagine: nothing else, glycine: carbonyl, glycine: carbonyl+lipid, none: amine+carbonyl+lipid; wonkwang_cho_inhee: arginine: amine+carbonyl | none |
| Benzaldehyde | sugar identity | conflict between laboratories | 18 | 0 | astar_sifbi_chiang: glucose: amine; hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine | none |
| Hexanal | sugar identity | conflict between laboratories | 19 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: arabinose: amine, none: amine, xylose: amine | none |
| Methional (3-(methylthio)propanal) | sugar identity | answered (level) | 13 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine, xylose: amine | none |
| Pyrazines | sugar identity | sign changes inside one laboratory | 5 | 1 | hfut_wei_zhaojun: arabinose: amine+other+sulfur; minnesota_reineccius: glucose: amine; reading_mottram_parker: glucose: amine+carbonyl+lipid | none |
| Trimethylpyrazine | amino acid identity | conflict between laboratories | 95 | 0 | ghent_de_kimpe: free amino acids (lysine): carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl; ho_network: none: amine+carbonyl; reading_mottram_parker: glycine: carbonyl, glycine: carbonyl+lipid; wonkwang_cho_inhee: none: amine+carbonyl | none |

R1 - R3 are PREREG amendment A6's frozen candidate rules, read here on contrasts already seen: a change is a finding about this table's sensitivity, not evidence for the rule. A reading needs three laboratories (PREREG sec. 5); 'changes' lists only readings that differ from the main table's.
