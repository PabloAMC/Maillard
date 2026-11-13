# PREREG A6 audit: the keyed pH control, final pH and base mixture (STATUS_AND_PLAN sec. 0, step 8)

*audit on sources already read; no reading of the lever table moves; PREREG A6 rules are scored only on sources keyed after the A6 commit. Spec: [STEP8_FIELDS_SPEC.md](../../docs/reference/lever_table/STEP8_FIELDS_SPEC.md).*

## Coverage

- pH contrasts: 306; `ph_control`: adjusted_unbuffered 153, buffered 148, unstated 5 (check 69, derived 237).
- final pH on every arm: 63 contrasts; final-pH spread of at least 0.5 pH unit over the arms the slope uses (A6-R2 slope reported): 59.
- categorical pairs: 1949; base mixture differs beyond the lever (A6-R3 drops them): 7; arms' base from check 498, derived 2633.
- reader-checked sources: Zhang2024, ames2001b, feng2022, grosch1992, hincelin1992, hofmann2000, hwang1995, martin2001, mottram2002, park2022mealworm, shibamoto1977, wei2019flax, xu2007, zhang2018soypep, zhou2023, zhou2023c.

## pH cells with two or more laboratories

| odorant | band | watched | labs | control (keyed) | main reading | R1 buffered | buffered set = free text's? | R1 unbuffered | R2 final pH (contrasts) | changes |
|---|---|---|---|---|---|---|---|---|---|---|
| 2-Acetylthiazole | pH 6-7.5, initial or label |  | 2 | buffered 4 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | too few laboratories (2) | none |
| 2-Furfurylthiol (FFT) | pH 4.5-6, initial or label |  | 4 | adjusted_unbuffered 1, buffered 4 | direction unclear | direction unclear (3 labs) | yes | too few laboratories (1 labs) | - (0) | R1 unbuffered band: too few laboratories |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | answered, fragile | 3 | adjusted_unbuffered 1, buffered 4 | answered by the literature | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | answered by the literature (3) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| 2-Methyl-3-furanthiol (MFT) | pH 4.5-6, initial or label | sign changes inside one laboratory | 2 | buffered 3 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | conflict between laboratories | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (2) | none |
| 2-Methylpyrazine | pH 6-7.5, initial or label |  | 2 | buffered 2 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | too few laboratories (1) | none |
| 2-Methylpyrazine | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 2, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | conflict between laboratories | 2 | adjusted_unbuffered 1, buffered 2 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | too few laboratories (2) | none |
| DMHF | pH 6-7.5, initial or label |  | 2 | buffered 4 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Furfural | pH 4.5-6, initial or label | answered, fragile | 3 | adjusted_unbuffered 15, buffered 2 | answered by the literature | too few laboratories (2 labs) | yes | too few laboratories (1 labs) | answered by the literature (1) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH 4.5-6, initial or label | sign changes inside one laboratory | 2 | adjusted_unbuffered 15, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |
| Phenylacetaldehyde | pH 6-7.5, initial or label | fragile | 3 | buffered 5, unstated 5 | direction unclear | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | R1 buffered band: too few laboratories; R1 unbuffered band: too few laboratories |
| Phenylacetaldehyde | pH <4.5, initial or label |  | 2 | buffered 9 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Pyrazines | pH 6-7.5, initial or label |  | 2 | buffered 2 | too few laboratories | too few laboratories (2 labs) | yes | too few laboratories (0 labs) | - (0) | none |
| Pyrazines | pH >7.5, initial or label |  | 2 | adjusted_unbuffered 1, buffered 1 | too few laboratories | too few laboratories (1 labs) | yes | too few laboratories (1 labs) | - (0) | none |

### Contrasts that carry a final pH (initial range against the pH the pots ended at)

| odorant | band | source | laboratory | initial or label pH | final pH by arm | slope, keyed basis | slope on final pH |
|---|---|---|---|---|---|---|---|
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.17 | +0.32 |
| 2-Acetylthiazole | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | +0.16 | +0.50 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -1.09 | -1.97 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.46 | -1.34 |
| 2-Furfurylthiol (FFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 8 | 3.2, 3.4, 5 | -0.11 | -0.22 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.96 | -1.81 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6, 6.3 | -0.27 | -0.80 |
| 2-Methyl-3-furanthiol (MFT) | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 7 | 3.2, 3.4, 5 | +0.36 | - |
| 2-Methylpyrazine | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | +0.35 | +0.64 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 8 | 5.5, 6.2, 6.8 | -0.71 | -1.33 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | xu2007 | cau_gao_yanxiang | 5.6 - 7 | 5.5, 6, 6.3 | +0.09 | +0.24 |
| Bis(2-methyl-3-furyl) disulfide | pH 6-7.5, initial or label | zhou2023 | ho_network | 6 - 7 | 3.2, 3.4, 5 | +0.24 | - |
| Furfural | pH 4.5-6, initial or label | mottram2002 | reading_mottram_parker | 4.2 - 5.6 | 5.5, 3.9 | -2.93 | -2.56 |

## Categorical cells: base mixture (A6-R3)

| odorant | lever | watched | pairs | base differs | reference pot by laboratory | changes |
|---|---|---|---|---|---|---|
| 2-Methyl-3-furanthiol (MFT) | sugar identity |  | 27 | 1 | dfa_garching: ribose: amine+sulfur, ribose: sulfur; hfut_wei_zhaojun: xylose: amine+other+sulfur; reading_mottram_parker: none: amine+sulfur, ribose 5-phosphate: amine+sulfur | none |
| 2-Methyl-3-furanthiol (MFT) | sulfur source | answered (level) | 15 | 0 | btbu_beijing: none: amine+carbonyl+sulfur; dfa_garching: hydrogen sulfide: other, none: amine+carbonyl+sulfur, none: amine+sulfur, thiamine: amine+sulfur, thiamine: nothing else; ho_network: glutathione: amine+carbonyl, none: amine+carbonyl+sulfur | none |
| 2-Methylpyrazine | amino acid identity | fragile | 136 | 0 | ghent_de_kimpe: free amino acids (lysine): carbonyl, gly-gly: carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl, hydrolysed protein: amine+carbonyl; ho_network: glutathione: amine+carbonyl, none: amine+carbonyl; njau_zhao_liyan: methionine: amine+carbonyl; reading_mottram_parker: asparagine: nothing else, glycine: carbonyl, glycine: carbonyl+lipid, none: amine+carbonyl+lipid; wonkwang_cho_inhee: arginine: amine+carbonyl | none |
| Benzaldehyde | sugar identity | conflict between laboratories | 18 | 0 | astar_sifbi_chiang: glucose: amine; hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine | none |
| Hexanal | sugar identity | conflict between laboratories | 19 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: arabinose: amine, none: amine, xylose: amine | none |
| Methional (3-(methylthio)propanal) | sugar identity | answered (level) | 13 | 0 | hfut_wei_zhaojun: xylose: amine+other+sulfur; njau_zhao_liyan: ribose: amine+sulfur; nus_liu_shaoquan: none: amine, xylose: amine | none |
| Pyrazines | sugar identity | sign changes inside one laboratory | 5 | 1 | hfut_wei_zhaojun: arabinose: amine+other+sulfur; minnesota_reineccius: glucose: amine; reading_mottram_parker: glucose: amine+carbonyl+lipid | none |
| Trimethylpyrazine | amino acid identity | conflict between laboratories | 95 | 0 | ghent_de_kimpe: free amino acids (lysine): carbonyl, glycine (2 mmol): carbonyl, glycine (3 mmol): carbonyl, glycine + lysine: carbonyl; ho_network: none: amine+carbonyl; reading_mottram_parker: glycine: carbonyl, glycine: carbonyl+lipid; wonkwang_cho_inhee: none: amine+carbonyl | none |

R1 - R3 are PREREG amendment A6's frozen candidate rules, read here on contrasts already seen: a change is a finding about this table's sensitivity, not evidence for the rule. A reading needs three laboratories (PREREG sec. 5); 'changes' lists only readings that differ from the main table's.
