# Cerny 1994 — EXTRACTION (beef water-solubles, alanine + hexoses, and alanine + methylglyoxal on kieselguhr in oil, pH 5.6, 180 C / 7 min; ethyldimethylpyrazines by isotope dilution)
### Alanine is the precursor of the ethyldimethylpyrazines in roasted beef; 2-oxopropanal + alanine is the most productive pair; the 3,6-isomer outnumbers the odour-active 3,5-isomer ten to one.

**Source on disk:** `data/articles/cerny1994.pdf` (owner's download, 2026-09-08). Read-only extraction from
the pypdf text layer in the scratchpad (OCR'd scan; "gg" = µg, "retool" = mmol, "rain" = min throughout);
Tables 1-5 are clean and re-typed in full below, with the OCR corrected where the meaning is unambiguous.
Figure 5 (mechanism) is a scheme, Figures 3-4 chromatograms; nothing numeric is on them. Repo status
before this dossier: roadmap §5b lists 2-ethyl-3,5-dimethylpyrazine with "the aldehyde-addition step"
missing; `2_ethyl_3_5_dimethylpyrazine` is a declared desirable target.

## 0. Identity

| field | value |
|---|---|
| Title | "Precursors of ethyldimethylpyrazine isomers and 2,3-diethyl-5-methylpyrazine formed in roasted beef" |
| Authors | Christoph Cerny, Werner Grosch (Deutsche Forschungsanstalt für Lebensmittelchemie, Garching) |
| Venue | Z. Lebensm. Unters. Forsch. 198 (1994) 210-214 |
| DOI | 10.1007/BF01192597 |
| Numerical key | **I** = 2-ethyl-3,6-dimethylpyrazine; **II** = 2-ethyl-3,5-dimethylpyrazine (repo `2_ethyl_3_5_dimethylpyrazine`, see flag 8); **III** = 2-ethyl-5,6-dimethylpyrazine; **IV** = 2,3-diethyl-5-methylpyrazine. "2-oxopropanal" = methylglyoxal (MGO). |
| Companions | Cerny & Grosch 1992, 1993 (roasted beef odorants and their SIDA); Arnoldi 1988 (fructose + eight amino acids, 120 C, 3 h: alanine, aspartate and valine as precursors of II) |

## 1. Why it matters

B18 makes pyrazine, methylpyrazine and 2,5-dimethylpyrazine from glyoxal, methylglyoxal and glycine
and declares the condensation fast. The roadmap's next step is the aldehyde addition that makes the
trialkylpyrazines, and this paper is the cleanest single-pair evidence for it: fed **MGO + alanine**
gives I, II, III, IV in that order, I : II about 10 : 1, in a matrix and at a temperature that are a
roast (180 C, 7 min, pH 5.6). The mechanism the authors give is exactly the rule the repository needs:
the Strecker of alanine on MGO yields aminoacetone, 2-aminopropanal and **acetaldehyde**; two
aminoketones condense to a dimethyl-dihydropyrazine; acetaldehyde adds to it and dehydrates. Which
dihydropyrazine (2,5- from the homo pair, 2,6- from the cross pair) sets which ethyl isomer forms, and
the 10 : 1 split is a within-study measurement of the aminoketone pool's composition. The paper also
shows the sugar side (glucose, fructose, their 6-phosphates all give the same pattern at 5x lower
level than MGO) and the enhancers (carnosine, lactic acid ~3.5x) a meat-like matrix carries.

## 2. Methods as they matter to a model

- **Reaction system:** reactants dissolved in **5 mL Na/K-phosphate 0.07 mol/L, pH 5.6** (checked and
  corrected); the solution **absorbed on 10 g kieselguhr** (HCl-washed), poured into **100 g hydrogenated
  peanut oil at 180 C**, stirred **7 min**, cooled 10 min at 10-15 C. This is a wet-solid roast in hot
  oil, not an aqueous pot: water leaves during the 7 min. Nominal starting concentrations for the
  2 mmol charges are 2 mmol / 5 mL = **400 mmol/L** (4 mmol: 800 mmol/L), before drying.
- **Charges:** Table 2 footnote a (beef LMF composition, referring to 300 g meat for expts 1-2 and a
  smaller basis for 3-7): glucose 214 mg (1.19 mmol), fructose 57 mg (0.32), ribose 20 mg (0.13),
  glucose 6-phosphate Na 429 mg (1.52), fructose 6-phosphate Na 118 mg (0.42); amino acids alanine
  63 mg (**0.71 mmol**), arginine 15 (0.09), aspartic acid 1.5 (0.01), cysteine 7.4 (0.06), glutamic acid
  20.6 (0.15), glycine 19.2 (0.26), histidine 9.7 (0.06), isoleucine 8.9 (0.07), leucine 16 (0.12),
  lysine·HCl 13.9 (0.08), methionine 4.6 (0.03), proline 7.8 (0.07), phenylalanine 10.3 (0.06), serine
  7.2 (0.07), threonine 6.3 (0.05), tyrosine 9.5 (0.05), valine 16 (0.14); further: glutamine 140 mg
  (0.96), carnosine 600 mg (2.65), creatine·H2O 570 mg (3.82), lactic acid 90 % 1 g (10 mmol).
  Tables 3-4: **2 mmol of each reactant** (lactic acid 10 mmol). Table 5: **2 mmol or 4 mmol each** of
  MGO and alanine.
- **Work-up:** labelled internal standards **d-II and d-IV** in 200 mL diethyl ether added after cooling;
  filtered; high-vacuum distillation (5 mPa) of the volatiles; acid extraction (0.1 M HCl 3 x 50 mL),
  alkalised to pH 12, back-extracted into ether (3 x 100 mL), dried, concentrated to 200 µL.
- **Quantification: stable isotope dilution by mass chromatography** (ion-trap GC-MS, DB-Wax; ions
  m/z 137, 140, 151, 154 in Figure 4). **I, II and III are calculated against d-II; IV against d-IV.**
  So only II and IV have their own labelled standard; I and III assume II's response (flag 1).
- **Replicates:** "mean values of duplicates" (Tables 2-4); Table 5 does not say.
- **Reporting unit:** **µg of pyrazine per pot** (the whole 5 mL charge). Conversion to mol % of the
  2 mmol alanine (MW I, II, III 136.19; IV 150.22): mol % = µg / MW / 2000 x 100 = µg / (20 x MW).
  Example: 256 µg I = 1.880 µmol = 0.094 mol % of alanine.
- **Odour thresholds (Table 1):** GC-olfactometry, three assessors, reference (E)-2-decenal 2.7 ng/L air.

## 3. Tables re-typed

### Table 1. Odour thresholds of ethyldimethylpyrazines (ng/L air; mean of three assessors)

| pyrazine | threshold |
|---|---:|
| 2-ethyl-3,6-dimethyl (I) | 2.45 |
| 2-ethyl-3,5-dimethyl (II) | 0.009 |
| 2-ethyl-5,6-dimethyl (III) | > 95 |

II is 270x below I and > 10 000x below III.

### Table 2. "Water-soluble components of beef as precursors of pyrazines I to IV" (µg formed after 7 min at 180 C; means of duplicates)

| expt | reaction system | I | II | III | IV |
|---:|---|---:|---:|---:|---:|
| 1 | water-soluble fraction (WF) of beef (300 g basis) | 43 | 15 | 5.5 | 7.4 |
| 2 | low-molecular-mass fraction (LMF, Mr <= 1000) | 57 | 23 | 10 | 4.4 |
| 3 | monosaccharides + amino acids | 11 | 3.1 | 0.6 | 2.6 |
| 4 | system 3 without alanine | < 0.1 | < 0.1 | < 0.1 | < 0.1 |
| 5 | monosaccharides + amino acids + carnosine + glutamine + creatine + lactic acid | 78 | 14 | 2.4 | 5.2 |
| 6 | monosaccharides + carnosine + glutamine + creatine + lactic acid (no free amino acids) | 4.9 | 1.5 | 0.1 | < 0.1 |
| 7 | system 6 + alanine | 80 | 6.9 | 1.1 | 5.3 |

### Table 3. "Alanine and different monosaccharides as precursors" (2 mmol each; µg after 7 min at 180 C; means of duplicates)

| expt | alanine plus | I | II | III | IV |
|---:|---|---:|---:|---:|---:|
| 8 | fructose | 49 | 3.9 | 0.5 | 8.8 |
| 9 | glucose | 26 | 2.1 | 0.4 | 5.3 |
| 10 | fructose 6-phosphate | 28 | 6.3 | 2.1 | 7.9 |
| 11 | glucose 6-phosphate | 25 | 3.0 | 1.6 | 4.4 |

### Table 4. "Fructose and different nitrogen sources as precursors" (2 mmol each, lactic acid 10 mmol; µg after 7 min at 180 C; means of duplicates)

| expt | fructose plus | I | II | III | IV |
|---:|---|---:|---:|---:|---:|
| 12 | alanine | 49 | 3.9 | 0.5 | 8.8 |
| 13 | carnosine | > 0.15 | < 0.15 | < 0.15 | < 0.15 |
| 14 | creatine | < 0.1 | < 0.1 | < 0.1 | < 0.1 |
| 15 | glutamine | < 0.1 | < 0.1 | < 0.1 | < 0.1 |
| 16 | alanine + carnosine | 174 | 8.3 | 0.5 | 32 |
| 17 | alanine + glutamine | 77 | 4.9 | 0.6 | 12 |
| 18 | alanine + lactic acid | 172 | 17 | 1.5 | 30 |
| 19 | alanine + ammonium formate | 40 | 3.6 | 0.4 | 15 |

### Table 5. "Reaction system 2-oxopropanal and alanine" (µg after 7 min at 180 C)

| expt | amount of each reactant | I | II | III | IV | sum | I share |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 2 mmol | 256 | 27 | 2.6 | 18 | 303.6 | 84.3 % |
| 21 | 4 mmol | 837 | 45 | 1.8 | 54 | 937.8 | 89.3 % |

Text: I accounted for 84 % of the pyrazines; doubling the charge raised I 3.3-fold and IV 3.0-fold,
II by 70 %.

### Tables 3 and 5 in mol % of alanine (mine; µg / (20 x MW))

| system | I | II | III | IV |
|---|---:|---:|---:|---:|
| MGO + Ala, 2 mmol each | 0.094 | 0.0099 | 0.00095 | 0.0060 |
| MGO + Ala, 4 mmol each (÷ 4000 µmol) | 0.154 | 0.0083 | 0.00033 | 0.0090 |
| fructose + Ala | 0.018 | 0.0014 | 0.00018 | 0.0029 |
| glucose + Ala | 0.0095 | 0.00077 | 0.00015 | 0.0018 |
| fructose 6-P + Ala | 0.0103 | 0.0023 | 0.00077 | 0.0026 |
| glucose 6-P + Ala | 0.0092 | 0.0011 | 0.00059 | 0.0015 |
| fructose + Ala + carnosine | 0.064 | 0.0030 | 0.00018 | 0.0107 |
| fructose + Ala + lactic acid | 0.063 | 0.0062 | 0.00055 | 0.0100 |

## 4. Numbers and steps the repository can use

Registry keys: II = `2_ethyl_3_5_dimethylpyrazine` (exists; SMILES caveat in flag 8); I
(2-ethyl-3,6-dimethylpyrazine), III (2-ethyl-5,6-dimethylpyrazine), IV (2,3-diethyl-5-methylpyrazine),
methylglyoxal, alanine, carnosine: **not in registry**. `acetaldehyde` exists.

| quantity or step | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| MGO + Ala -> I (2-ethyl-3,6-dimethylpyrazine) | 256 µg = 0.094 mol % of Ala | µg per pot | 2 + 2 mmol on kieselguhr in oil, pH 5.6 buffer, 180 C, 7 min, SIDA (against d-II) | Table 5 expt 20 | fed_intermediate_yield |
| MGO + Ala -> II (2-ethyl-3,5-dimethylpyrazine) | 27 µg = 0.0099 mol % | µg per pot | same | Table 5 expt 20 | fed_intermediate_yield (own labelled standard) |
| MGO + Ala -> III | 2.6 µg = 0.00095 mol % | µg per pot | same | Table 5 | fed_intermediate_yield |
| MGO + Ala -> IV (2,3-diethyl-5-methylpyrazine) | 18 µg = 0.0060 mol % | µg per pot | same | Table 5 | fed_intermediate_yield (own labelled standard) |
| MGO + Ala at double charge | I 837, II 45, III 1.8, IV 54 (0.154 / 0.0083 / 0.00033 / 0.0090 mol %) | µg per pot | 4 + 4 mmol, same | Table 5 expt 21 | fed_intermediate_yield |
| charge dependence 4 mmol / 2 mmol | I 3.3x; IV 3.0x; II 1.7x; III 0.7x | — | | Table 5 | within_study_ratio (I and IV rise faster than linearly in charge, II slower: the two isomers do not share one rate law in this pot) |
| I : II from MGO + Ala | 9.5 (2 mmol); 18.6 (4 mmol) | — | | Table 5 | within_study_ratio |
| I : II from sugar + Ala | 12.6 (Fru); 12.4 (Glc); 4.4 (Fru-6-P); 8.3 (Glc-6-P) | — | 2 + 2 mmol | Table 3 | within_study_ratio |
| implied aminoketone split (mine) | ~95-97 % aminoacetone if acetaldehyde adds equally fast to the 2,5- and 2,6-dimethyldihydropyrazines: homo : cross = (p² + (1-p)²) : 2p(1-p) = 9.5 -> p = 0.95; 18.6 -> 0.97; 12.5 -> 0.96 | — | | derived from Table 5 / Table 3 | derived; supports one aminoketone per dicarbonyl (B18) |
| MGO vs fructose as the carbonyl, same Ala | I 5.2x; II 6.9x; III 5.2x; IV 2.0x | — | 2 + 2 mmol | Table 5 vs Table 4 expt 12 | within_study_ratio (the fed dicarbonyl vs the sugar's own supply of it in 7 min at 180 C) |
| hexose identity | I 49 / 26 / 28 / 25 µg (Fru / Glc / F6P / G6P) | µg per pot | 2 + 2 mmol | Table 3 | level_only (within 2x; "without preference for one of the hexoses") |
| fructose + Ala -> I, II, III, IV | 49, 3.9, 0.5, 8.8 (0.018 / 0.0014 / 0.00018 / 0.0029 mol %) | µg per pot | 2 + 2 mmol, 180 C, 7 min | Table 3 expt 8 | level_only (end of a 7 min roast; validation) |
| glucose + Ala -> I, II, III, IV | 26, 2.1, 0.4, 5.3 (0.0095 / 0.00077 / 0.00015 / 0.0018 mol %) | µg per pot | same | Table 3 expt 9 | level_only |
| alanine requirement | all four < 0.1 µg without alanine vs 11 / 3.1 / 0.6 / 2.6 with it, in the full amino-acid + sugar mix | µg per pot | Table 2 expts 3 vs 4 | | within_study_ratio (>= 110x for I; >= 31x for II) |
| enhancers on fructose + Ala | carnosine: I 3.6x, II 2.1x, IV 3.6x; lactic acid (10 mmol): I 3.5x, II 4.4x, IV 3.4x; glutamine: I 1.6x, II 1.3x, IV 1.4x; ammonium formate: 0.8x / 0.9x / 1.7x | — | Table 4 expts 16-19 vs 12 | | within_study_ratio (not a nitrogen effect — NH4 formate does nothing; carnosine's imidazole and lactate act as catalysts) |
| carnosine, creatine, glutamine alone with fructose | I, II, III, IV all <= 0.15 µg | µg per pot | Table 4 expts 13-15 | | level_only (non-detects: a peptide/imidazole without alanine gives none of these) |
| beef water-solubles | I 43 / 57, II 15 / 23, III 5.5 / 10, IV 7.4 / 4.4 (WF / LMF, 300 g meat) | µg per pot | Table 2 expts 1-2 | | level_only (the food-relevant absolute level) |
| odour thresholds | I 2.45; II 0.009; III > 95 | ng/L air | GC-O, 3 assessors | Table 1 | threshold |

**The step in the authors' words (Figure 5).** Strecker degradation of alanine by MGO yields
aminoacetone, 2-aminopropanal and acetaldehyde. Condensation of two aminoacetones (or two
2-aminopropanals) gives 2,5-dimethyl-dihydropyrazine; its reaction with acetaldehyde and dehydration
gives I. Condensation of aminoacetone with 2-aminopropanal gives the 2,6-dimethyl-dihydropyrazine,
which with acetaldehyde gives II. "In the model reactions, pyrazine I was always formed as the major
and II as the minor product ... more aminoacetone is produced from 2-oxopropanal, because the aldehyde
group of the latter is more reactive than the keto group." IV cannot be explained from MGO + alanine
by the Bemis-Young route (2-aminopropanal + 2-hydroxy-3-amino-4-hexanone) and is left open.

**What this gives the aldehyde-addition step.** One pair (MGO + alanine), one temperature, one time,
absolute molar yields of all four trialkylpyrazines against the same alanine charge, the isomer split,
and the sign of the charge dependence. Against the B18 output (2,5-dimethylpyrazine from 2 aminoacetone)
this is the branch ratio "dihydropyrazine + acetaldehyde -> I" versus "dihydropyrazine -> oxidation ->
2,5-dimethylpyrazine", but 2,5-dimethylpyrazine itself is NOT quantified here, so the branch ratio
needs another paper (Adams 2008 gives the 2,5-DMP side from the same pair in water).

## 5. Flags

1. **I and III are quantified against d-II** (no labelled I or III); their levels assume II's response.
   II and IV have their own deuterated standards and are the reliable numbers.
2. **The medium is a kieselguhr roast in 180 C oil**, not an aqueous pot; the water in the 5 mL buffer
   leaves during the 7 min, so "400 mmol/L" is nominal and the effective temperature history of the
   solid is unknown. A rate cannot be read from this; a yield per charge can.
3. **One time point (7 min).** No time series, no temperature ladder.
4. **Table 5 does not state duplicates**; Tables 2-4 are means of duplicates without spread.
5. **The doubled charge (4 mmol) is not the same pot at twice the concentration** — twice the
   material on the same 10 g kieselguhr in the same oil; the > 3x rise of I and IV is real but its
   reading as reaction order is not clean.
6. **Acetaldehyde's origin is inferred** (Strecker aldehyde of alanine); there is no labelling and no
   run with added acetaldehyde or with another amino acid + MGO. Adams 2008 shows I-type products from
   amino acids that make no acetaldehyde (arginine, lysine), so the aldehyde pool is wider than the
   amino acid's own.
7. **pH 5.6 is the buffer's starting pH**; 0.07 M phosphate is weak against 2 mmol alanine and 10 mmol
   lactic acid (expt 18).
8. **Registry SMILES check (mine, RDKit on the host):** `compounds.yml` gives `2_ethyl_3_5_dimethylpyrazine`
   (CAS 13360-64-0, correct for II) the SMILES `CCc1nc(C)cnc1C`, which is isomer **I**
   (2-ethyl-3,6-dimethylpyrazine, InChIKey WHMWOHBXYIZFPF); II is `CCc1ncc(C)nc1C` (JZBCTZLGKSYRSF); the
   entry's InChI string is a C7 species. Since this paper's whole point is that I is 270x less potent
   than II, the key must be corrected before any I / II number is attached to it. Not edited here.
9. **Beef basis of Table 2 shifts** (300 g for expts 1-2, "a smaller amount (200 g)" for expt 3), so
   expts 1-2 and 3-7 are not on one basis; the authors say the difference alone cannot explain the drop.
