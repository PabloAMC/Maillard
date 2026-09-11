# Amrani-Hemaimi, Cerny & Fay 1995 — EXTRACTION (2 mmol sugar + 2 mmol amino acid in 5 mL pH 5.6 phosphate on 5 g kieselguhr, 7 min at 180 C in coconut oil; alanine or glycine x glucose or fructose, unlabelled and 13C-labelled)

### Eleven alkylpyrazines with their within-pot percentage shares and forty 13C-labelling fractions, and no absolute quantity of anything: the paper is a carbon-origin bookkeeping study, not a yield study.

**Source on disk:** `data/articles/amrani-hemaimi1995.pdf` (5 pp., J. Agric. Food Chem. 43:2818-2822,
1995). The text layer is an Acrobat Capture OCR layer with heavy glyph damage in prose
("[3-l3C1alanine", "9Oe", "3100", "x 1") and Table 3's ion columns were lost entirely. **Tables 1
and 2 were verified cell by cell against a 200-dpi raster of printed page 2819**
(`scratchpad/img/ah95-2.png`) and **Table 3 against a raster of page 2820**
(`scratchpad/img/ah95t3-3.png`); every number below matches the raster. Figure 1 (a CI mass spectrum
at three labelling contents) and Figures 2-5 (proposed mechanism schemes) are FIGURE-ONLY and carry
no numbers; the schemes were read for structure, not for values.

**Repository status before this dossier:** an `amrani-hemaimi1995_extraction.md` existed at wave Z3
and is quoted verbatim in `k3_final_parameter_inventory.md` C.15 ("there is NO absolute
quantification anywhere in this paper"), but the file is not in `data/lit/extraction_dossiers/`
today. `k3` carries the paper only through second-hand rows B10.13, B10.14, the FIT/HOLD-OUT ledger
lines and flag 28. `kinetic_core_b19_prereg_draft.md` section 3 names Table 2 as "on disk, never
used" and section 5 names it as one of three ratio sources; `kinetic_core_b22_prereg.md` section 6
asks for exactly such ratios after its own fit was refuted.

## 0. Identity

| field | value |
|---|---|
| Title | "Mechanisms of Formation of Alkylpyrazines in the Maillard Reaction" |
| Authors | Miriam Amrani-Hemaimi (Institute of Organic Chemistry, University of Lausanne) and Christoph Cerny (corresponding) and Laurent B. Fay (Nestle Research Centre, Nestec Ltd., Lausanne) |
| Venue | J. Agric. Food Chem. 1995, 43 (11), 2818-2822. Received December 23 1994; revised July 31 1995; accepted September 11 1995; abstract in Advance ACS Abstracts October 15 1995 |
| Article id | JF9407228 |
| Naming | Experiments are numbered 1-8 and the tables use those numbers as column heads: **1** = alanine + glucose, **2** = alanine + fructose, **3** = glycine + glucose, **4** = glycine + fructose (unlabelled, Table 1); **5** = [3-13C]alanine + glucose, **6** = [3-13C]alanine + fructose, **7** = [2-13C]glycine + glucose, **8** = [2-13C]glycine + fructose (Table 2). Compounds are numbered **1-11** in Table 1's order and referred to by number throughout the prose. "2-oxopropanal" is the paper's name for methylglyoxal |
| Same laboratory | Cerny & Grosch 1994 (Z. Lebensm. Unters. Forsch. 198:210) supplies the apparatus and the beefsteak-frying conditions; `cerny1994_extraction.md`, `cerny2004_extraction.md`, `cerny2007_extraction.md` are the same group |
| Cited in repo | `k3_final_parameter_inventory.md` B10.13, B10.14, C.15, flag 28, and the FIT / HOLD-OUT ledger lines at 1349-1351 |

## 1. Why it matters

This is the paper `kinetic_core_b19_prereg_draft.md` section 5 names as one of the three
within-study identity sources for an amino-acid layer on the fitted Strecker step, and the one it
calls "stranded since B2". It matters in two different ways, and it is important not to confuse
them.

**What it settles (structure).** The engine's pyrazine step
(`src/kinetic_core/parameters_pyrazine.py`) is built on the claim, stated in that module's own
docstring, that "the ring carbons are the dicarbonyl's, the amino acid gives the nitrogen". This
paper is the direct isotope test of that claim, and it both confirms and bounds it:
methylpyrazine and 2,5-dimethylpyrazine take **0 %** of their carbon from either labelled amino acid
in all four labelled pots, so for the two compounds the trunk actually makes (`AKG + AKG` and
`AKM + AKM`) the declared topology is measured, not assumed. But 2,6-dimethylpyrazine takes
**25-30 %** of one methyl from glycine's C-2, and the trimethyl- and ethyl-substituted pyrazines
take **20-100 %** of an alkyl group from the amino acid — so the aldehyde-addition step that
`kinetic_core_b19_prereg_draft.md` section 2 wants to build is real, is measured here as a partition
fraction, and applies to the compounds beyond the trunk's three.

**What it does not settle (rate or yield).** There is no internal standard for absolute
quantification, no response factor, no concentration, no µg, no ppb, no mol %, one experiment per
condition and no replicates. Table 1 is percent of the sum of the eleven peaks **within one
experiment**, from NPD peak areas. Section 4 states the consequence for the identity ratio the two
pre-registrations ask for; the short version is that a per-amino-acid yield ratio **cannot** be
taken from this paper, and that what can be taken is a set of within-pot shares and pathway
partitions.

Its other limits are worth stating with the numbers: one temperature (180 C), one time (7 min), one
pH (5.6 initial), and a pot whose water evaporates during the heat, so it is neither an aqueous nor
a controlled-a_w system. Nothing here can be transferred as a rate to the trunk's 100-145 C window.

## 2. Methods as they matter to a model

- **Reactants and loading.** Eight pots, each "sugar and amino acid ... dissolved in 5 mL of
  sodium/potassium phosphate buffer (pH 5.6; 0.07 mol/L)", at **2 mmol of amino acid + 2 mmol of
  sugar** per pot. That is **400 mmol/L of each before the water leaves** (mine: 2 mmol / 5 mL),
  1:1 molar, in **70 mmol/L phosphate at pH 5.6**. The paper prints the millimoles, not the
  molarity.
- **Support and heating.** "The solution was adsorbed on kieselguhr (5 g) and heated for 7 min in
  coconut oil at 180 C, using the conditions and the apparatus recently described (Cerny and
  Grosch, 1994). **During heating, the water evaporated from the mixture.** After heating, the
  mixture was cooled to 15 C for 10 min." So: a thin film on 5 g of diatomaceous earth, oil bath at
  180 C, 7 minutes, drying as it goes. There is no stirring, no headspace control, no a_w
  measurement, and no statement of the sample temperature (only the bath's).
- **Isolation.** 100 mL diethyl ether added to each heated mixture; filtered to remove kieselguhr
  and insolubles; high-vacuum distillation at 0.1 Pa, dropwise, into a liquid-nitrogen trap (Jung
  1992); a **basic fraction** taken by extracting the distillate with 0.1 mol/L HCl (3 x 50 mL),
  raising the aqueous phase to pH 12 and re-extracting with ether (3 x 100 mL); dried over Na2SO4;
  concentrated to 0.2 mL by Vigreux then Bemelmans microdistillation. **Every step of this is a
  recovery loss of unknown and compound-dependent size, and none of it is calibrated.**
- **Gas chromatography.** Carlo Erba Mega 2, cold on-column injector, DB-Wax 30 m x 0.32 mm i.d.,
  0.25 µm film; **NPD and FID running simultaneously**; He 2.5 mL/min; 1 µL injected; detector
  240 C; 35 C for 2 min, 40 C/min to 50 C, 6 C/min to 180 C, 10 C/min to 240 C, 240 C for 20 min.
  Retention indices by van den Dool & Kratz.
- **Quantification: NONE, and the table says so.** Table 1's footnote c reads "Based on the
  comparison of **NPD peak areas** in one experiment". So each column is normalised inside its own
  experiment; the four columns are four separate normalisations. There is no absolute area, no
  internal standard, no response factor, and **no total** printed for any pot. This is the single
  most important methodological fact in the paper for a kinetic repository, and section 5 flag 1
  restates it.
- **Identification.** Retention index plus EI (70 eV) and positive-CI (ammonia) mass spectra against
  reference substances, on a HP-5890 / Finnigan TSQ 700 (DB-Wax, cold on-column; 50 C 1 min, 6 C/min
  to 150 C, 30 C/min to 240 C, 240 C 1 min). Compound 11 (3,5-diethyl-2-methylpyrazine) was
  "identified only by comparison with data from the library of mass spectra" — no authentic
  standard.
- **Labelling measurement.** GC-MS after chemical ionisation, "to minimize the fragmentation".
  The distribution was computed "from the intensities of the ions [M + 1]+ (protonated molecular
  ion), [M + 2]+ (protonated molecular ion singly labeled), and [M + 3]+ (protonated molecular ion
  doubly labeled)", corrected twice: for "the natural 13C content of the corresponding unlabeled
  reference compounds" and for "the small extent of M+ found in all spectra". **No uncertainty is
  printed on any labelling fraction**, and the figures are quoted to two significant digits at most
  (0, 20, 25, 30, 40, 50, 70, 80, 90, 100).
- **Labelling position.** GC-MS/MS on the CI-generated protonated molecular ion, collision energy
  29 eV in the laboratory frame, argon at 0.1 Pa. The assignment logic is printed for
  3-ethyl-2,5-dimethylpyrazine: the labelled compound loses 13CH3, "This fragment comes preferentially
  from the breakdown of the bond C1-C2 of the ethyl group rather than from the breakdown between the
  methyl groups and the pyrazine cycle ... corroborated by the presence of the unlabeled ions at
  m/z 107, 80, and 42, which do not carry the ethyl group."
- **Reaction model the authors propose (Figures 2-5, structures only).** Two a-aminocarbonyl
  compounds (from the Strecker reaction of an amino acid with an a-dicarbonyl such as
  2-oxopropanal) condense to a **dihydropyrazine**; a Strecker **aldehyde** then adds to the
  dihydropyrazine and the adduct dehydrates and oxidises to the alkylpyrazine. Formaldehyde (from
  glycine) + 2,5-dimethyldihydropyrazine -> trimethylpyrazine; acetaldehyde (from alanine) +
  methyldihydropyrazine -> the ethylmethylpyrazines; acetaldehyde + 2,5- or
  2,6-dimethyldihydropyrazine -> 3-ethyl-2,5- and 2-ethyl-3,5-dimethylpyrazine; acetaldehyde +
  2-ethyl-5-methyldihydropyrazine -> the diethylmethylpyrazines. **This is the step
  `kinetic_core_b19_prereg_draft.md` section 2 calls "the aldehyde-addition step to the aminoketone
  pair", and this paper is its structural source and its only measured branch fraction.**

## 3. Tables re-typed

### Table 1. "Pyrazines Formed in the Model Reactions with Alanine, Glycine, Glucose, and Fructose"

Columns as printed. RI = retention index on DB-Wax. The four numeric columns are the four
experiments: **1** alanine + glucose, **2** alanine + fructose, **3** glycine + glucose, **4**
glycine + fructose. Values are "% of total pyrazines", footnote c: "Based on the comparison of NPD
peak areas in one experiment."

| no. | compound | RI | 1 (Ala+Glc) | 2 (Ala+Fru) | 3 (Gly+Glc) | 4 (Gly+Fru) |
|---:|---|---:|---:|---:|---:|---:|
| 1 | methylpyrazine | 1160 | 8 | 11 | 24 | 12 |
| 2 | 2,5-dimethylpyrazine | 1319 | 25 | 35 | 30 | 33 |
| 3 | 2,6-dimethylpyrazine | 1325 | 4 | 6 | 9 | 11 |
| 4 | 2-ethyl-5-methylpyrazine | 1380 | 10 | 6 | <1 | <1 |
| 5 | 2-ethyl-6-methylpyrazine | 1386 | 9 | 7 | <1 | 1 |
| 6 + 7 | trimethylpyrazine + 2-ethyl-3-methylpyrazine (co-eluted, footnote d) | 1400 | 12 | 7 | 35 | 42 |
| 8 | 3-ethyl-2,5-dimethylpyrazine | 1441 | 20 | 19 | 0 | 0 |
| 9 | 2-ethyl-3,5-dimethylpyrazine | 1455 | 3 | 2 | 0 | 0 |
| 10 | 2,3-diethyl-5-methylpyrazine | 1487 | 6 | 4 | 0 | 0 |
| 11 | 3,5-diethyl-2-methylpyrazine (footnote e) | 1505 | 3 | 1 | 0 | 0 |

Footnotes as printed: *a* identification by comparison with reference substance on retention index
and MS (EI) and MS (CI); *b* retention indices on DB-Wax per van den Dool and Kratz (1963); *c* as
quoted above; *d* "Peaks were not separated; given values represent the sum of both pyrazines";
*e* "The compound was identified only by comparison with data from the library of mass spectra."

**Column sums (mine):** experiment 1 = 100; experiment 2 = **98**; experiment 3 = 98 plus two "<1"
cells; experiment 4 = 99 plus one "<1" cell. Experiment 2's two-point shortfall is either rounding
across ten cells or a dropped unit in one cell; it is printed as re-typed and the raster confirms
it. **There is no total row: the four columns cannot be compared to one another as amounts.**

**Derived shares (mine, all within one column, so response-factor-immune only to the extent that
the NPD responds equally to all eleven pyrazines — see flag 2):**

| quantity | Ala+Glc | Ala+Fru | Gly+Glc | Gly+Fru |
|---|---:|---:|---:|---:|
| ethyl-bearing share (4+5+7 is unresolvable; 4+5+8+9+10+11) | 51 % | 39 % | <2 % | <2 % |
| methyl-only share (1+2+3+6, with 6 inseparable from 7) | 49 % | 59 % | 98 % | 98 % |
| 2,5-DMP / methylpyrazine | 3.1 | 3.2 | 1.3 | 2.8 |
| 2,5-DMP / 2,6-DMP | 6.3 | 5.8 | 3.3 | 3.0 |

### Table 2. "13C-Labeling Content and Labeling Position of Pyrazines Formed from 13C-Labeled Amino Acids"

Columns 5-8 are the four labelled experiments: **5** [3-13C]alanine + glucose, **6** [3-13C]alanine
+ fructose, **7** [2-13C]glycine + glucose, **8** [2-13C]glycine + fructose. Values are "% of
13C-labeled pyrazine". ND = "Not detected by MS". Footnote c: for row 6+7, "Peaks were not
separated; values represent the average of both pyrazines." Footnote e: "Including 20% of doubly
labeled compound."

| no. | compound | 5 (Ala+Glc) | 6 (Ala+Fru) | 7 (Gly+Glc) | 8 (Gly+Fru) | labeling position |
|---:|---|---:|---:|---:|---:|---|
| 1 | methylpyrazine | 0 | 0 | 0 | 0 | — |
| 2 | 2,5-dimethylpyrazine | 0 | 0 | 0 | 0 | — |
| 3 | 2,6-dimethylpyrazine | 0 | 0 | 30 | 25 | methyl group (7, 8) |
| 4 | 2-ethyl-5-methylpyrazine | 70 | 70 | 30 | 0 | C-2 of ethyl group (5, 6) |
| 5 | 2-ethyl-6-methylpyrazine | 30 | 20 | 0 | 0 | C-2 of ethyl group (5, 6) |
| 6 + 7 | trimethylpyrazine + 2-ethyl-3-methylpyrazine (co-eluted, footnote c) | 50 | 40 | 80 | 100 | — |
| 8 | 3-ethyl-2,5-dimethylpyrazine | 100 | 100 | ND | ND | C-2 of ethyl group (5, 6) |
| 9 | 2-ethyl-3,5-dimethylpyrazine | 70 | 70 | ND | ND | C-2 of ethyl group (5, 6) |
| 10 | 2,3-diethyl-5-methylpyrazine | 90 (footnote e) | 90 (footnote e) | ND | ND | C-2 of one or both ethyl groups (5, 6) |
| 11 | 3,5-diethyl-2-methylpyrazine | 90 (footnote e) | 90 (footnote e) | ND | ND | C-2 of one or both ethyl groups (5, 6) |

The "labeling position" column's parenthesised numbers are the footnote-b device: "Corresponding
experiments are given in parentheses", i.e. the position statement was established in experiments
5 and 6 (alanine) or 7 and 8 (glycine).

**Prose that fixes the reading of Table 2** (page 2820, verbatim fragments): "The ethylmethylpyrazines
4 and 5 formed during the reaction with [3-13C]alanine were 20-70% 13C-labeled; pyrazine 8 was even
100% 13C-labeled. Both 10 and 11 were 70% mono-13C-labeled and 20% bi-13C-labeled." So the "90"
entries for 10 and 11 are 70 % singly labelled plus 20 % doubly labelled, which is what footnote e
says. "Obviously for the formation of 8 one single reaction route exists, since it was 100%
13C-labeled. In contrast, 4, 5, and 9 are formed by at least one additional pathway without the
amino acid carbon." And for glycine: "Compound 6+7 from the reaction with [2-13C]glycine was 80-100%
13C-labeled, and 3 was 25-30% 13C-labeled ... for the formation of 3 there are at least two routes,
one with and another without incorporation of the amino acid carbon."

**A discrepancy in row 4 that the raster confirms is printed:** 2-ethyl-5-methylpyrazine shows
**30 %** label in experiment 7 ([2-13C]glycine + glucose), while Table 1 records that same compound
at **<1 %** of the glycine + glucose slate. A 30 % labelling fraction on a peak at under one percent
of the total is at the edge of what the CI ion-intensity method can support, and the paper offers no
comment. Flag 5.

### Table 3. "Characteristic Ions in the Mass Spectra of the Labeled Pyrazines and Their Corresponding Unlabeled Analogues"

The text layer lost these columns entirely; re-typed from the page-2820 raster. "[M + H]+ diagnostic
ions, m/z (intensity, %)". Footnote a: "The spectra were obtained by GC-MS/MS after
collision-induced dissociation of the protonated molecular ion generated by positive chemical
ionization." *b* "Not analyzed. Diagnostic ions of unlabeled 2-ethyl-6-methylpyrazine were used."
*c* "Not analyzed. Diagnostic ions of unlabeled 2,3-diethyl-5-methylpyrazine were used."

| compound | no. | [M + H]+ diagnostic ions, m/z (intensity, %) |
|---|---:|---|
| unlabeled 2-ethyl-5-methylpyrazine | 4 | na (footnote b) |
| labeled 2-ethyl-5-methylpyrazine | | 124 (70), 109 (30), 108 (100), 107 (40), 83 (10), 80 (15) |
| unlabeled 2-ethyl-6-methylpyrazine | 5 | 123 (100), 108 (39), 107 (3), 82 (18), 80 (5), 67 (7), 42 (26) |
| labeled 2-ethyl-6-methylpyrazine | | 124 (71), 109 (100), 108 (29), 83 (11), 80 (4), 68 (8), 42 (18) |
| unlabeled 3-ethyl-2,5-dimethylpyrazine | 8 | 137 (39), 122 (100), 121 (70), 107 (8), 94 (9), 80 (20), 42 (11) |
| labeled 3-ethyl-2,5-dimethylpyrazine | | 138 (31), 122 (100), 121 (70), 107 (8), 95 (11), 80 (22), 42 (13) |
| unlabeled 2-ethyl-3,5-dimethylpyrazine | 9 | 137 (37), 122 (100), 121 (67), 107 (7), 94 (6), 80 (20), 42 (12) |
| labeled 2-ethyl-3,5-dimethylpyrazine | | 138 (30), 122 (100), 121 (71), 107 (5), 95 (7), 80 (22), 42 (13) |
| unlabeled 2,3-diethyl-5-methylpyrazine | 10 | 151 (27), 136 (34), 135 (100), 134 (0), 121 (7), 108 (26), 107 (16), 82 (3) |
| labeled 2,3-diethyl-5-methylpyrazine | | 152 (41), 137 (39), 136 (100), 135 (36), 121 (38), 108 (17), 107 (23), 82 (9) |
| unlabeled 3,5-diethyl-2-methylpyrazine | 11 | na (footnote c) |
| labeled 3,5-diethyl-2-methylpyrazine | | 152 (24), 137 (20), 136 (100), 135 (63), 121 (61) |

These are identification evidence, not kinetic numbers, and none of them enters section 4.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`):** methylpyrazine -> `methylpyrazine`;
2,5-dimethylpyrazine -> `2_5_dimethylpyrazine`; 2,6-dimethylpyrazine -> `2_6_dimethylpyrazine`;
trimethylpyrazine -> `trimethylpyrazine`; 2-ethyl-3,5-dimethylpyrazine ->
`2_ethyl_3_5_dimethylpyrazine`. **Not in the registry:** 2-ethyl-5-methylpyrazine,
2-ethyl-6-methylpyrazine, 2-ethyl-3-methylpyrazine, 3-ethyl-2,5-dimethylpyrazine,
2,3-diethyl-5-methylpyrazine, 3,5-diethyl-2-methylpyrazine — six of the eleven, and five of the six
are the ethyl-bearing compounds this paper exists to explain. Alanine, glycine, glucose, fructose,
formaldehyde and acetaldehyde are not in `compounds.yml` either.

**All rows share these conditions:** 400 mmol/L amino acid + 400 mmol/L sugar (before drying),
70 mmol/L Na/K phosphate initially pH 5.6, adsorbed on kieselguhr, 7 min at 180 C in coconut oil,
water evaporating during the heat; basic fraction after high-vacuum distillation; NPD peak areas
normalised within one experiment; one experiment per condition, no replicates, no error bars.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| pyrazine distribution, alanine + glucose | 8 / 25 / 4 / 10 / 9 / 12 / 20 / 3 / 6 / 3 (compounds 1-11, 6+7 merged) | % of the eleven-peak sum | as above | Table 1 col. 1, p. 2819 | within_study_ratio |
| pyrazine distribution, alanine + fructose | 11 / 35 / 6 / 6 / 7 / 7 / 19 / 2 / 4 / 1 | % of sum | as above | Table 1 col. 2 | within_study_ratio (column sums to 98) |
| pyrazine distribution, glycine + glucose | 24 / 30 / 9 / <1 / <1 / 35 / 0 / 0 / 0 / 0 | % of sum | as above | Table 1 col. 3 | within_study_ratio |
| pyrazine distribution, glycine + fructose | 12 / 33 / 11 / <1 / 1 / 42 / 0 / 0 / 0 / 0 | % of sum | as above | Table 1 col. 4 | within_study_ratio |
| ethyl-bearing share of the slate: alanine vs glycine | 51 % and 39 % (alanine, glucose and fructose) against <2 % (glycine, both sugars) | % of sum (mine, from Table 1) | as above | Table 1 | within_study_ratio |
| ON/OFF switch: compounds 8-11 with glycine | 0, 0, 0, 0 in both glycine pots against 20/3/6/3 and 19/2/4/1 with alanine | % of sum | as above | Table 1 + text p. 2819 ("The pyrazines 8-11 were not formed when glycine was used as the nitrogen source") | within_study_ratio (a structural zero, not a small number) |
| carbon origin of methylpyrazine and 2,5-dimethylpyrazine | 0 | % of the molecule 13C-labelled, from the amino acid | all four labelled pots | Table 2 rows 1-2 | within_study_ratio (a measured zero) |
| carbon origin of 2,6-dimethylpyrazine, from glycine C-2 | 30 (glucose), 25 (fructose); 0 from alanine C-3 | % 13C-labelled, one methyl group | as above | Table 2 row 3 | within_study_ratio |
| carbon origin of 2-ethyl-5-methylpyrazine, from alanine C-3 | 70, 70 | % 13C-labelled at C-2 of the ethyl group | glucose, fructose | Table 2 row 4 | within_study_ratio |
| carbon origin of 2-ethyl-6-methylpyrazine, from alanine C-3 | 30, 20 | % 13C-labelled at C-2 of the ethyl group | glucose, fructose | Table 2 row 5 | within_study_ratio |
| carbon origin of trimethylpyrazine + 2-ethyl-3-methylpyrazine (co-eluted) | 50, 40 (alanine); 80, 100 (glycine) | % 13C-labelled, average of two compounds | glucose, fructose | Table 2 row 6+7 | within_study_ratio — **DO NOT USE as a single species** (flag 4) |
| carbon origin of 3-ethyl-2,5-dimethylpyrazine, from alanine C-3 | 100, 100 | % 13C-labelled at C-2 of the ethyl group | glucose, fructose | Table 2 row 8 | within_study_ratio (a measured unity: "one single reaction route exists") |
| carbon origin of 2-ethyl-3,5-dimethylpyrazine, from alanine C-3 | 70, 70 | % 13C-labelled at C-2 of the ethyl group | glucose, fructose | Table 2 row 9 | within_study_ratio |
| carbon origin of 2,3-diethyl-5-methylpyrazine and 3,5-diethyl-2-methylpyrazine | 90, 90 each (70 singly + 20 doubly labelled) | % 13C-labelled | glucose, fructose | Table 2 rows 10-11 + footnote e + text p. 2820 | within_study_ratio |
| sugar effect (glucose vs fructose) on the labelling fractions | "The nature of the sugar had no pronounced effect" — the pairs differ by 0-20 points | — | four labelled pots | Table 2 + text p. 2820 | within_study_ratio |
| any absolute yield, concentration, mass or rate | **NOT PRINTED ANYWHERE** | — | — | — | — |
| Figure 1 (CI spectra of 2-ethyl-5-methylpyrazine at 0, 30 and 70 % label); Figures 2-5 (mechanism schemes) | — | — | — | Figures 1-5 | figure_only (Figures 2-5 carry no numbers at all) |

### Does this paper support an amino-acid IDENTITY RATIO? — the answer, stated plainly

**No, not in the form `kinetic_core_b19_prereg_draft.md` section 5 and `kinetic_core_b22_prereg.md`
section 6 need.** An identity ratio in that sense is one amino acid's yield of a product over
another amino acid's yield of the same product in a comparable pot. Taking it here would need the
alanine pot's total and the glycine pot's total, and **neither is printed**: Table 1's columns are
each normalised to their own hundred from NPD peak areas, and there is no total area, no internal
standard, no response factor and no concentration anywhere in the five pages. The paper's own words
close the door — footnote c, "Based on the comparison of NPD peak areas **in one experiment**". The
repository's own earlier reading, quoted in `k3_final_parameter_inventory.md` C.15, says the same
thing.

**What the ratio would be if one assumed the two pots make the same total pyrazine, and why that
assumption must be refused.** Under that assumption the Table 1 percentages become yield ratios
directly:

| compound | glycine / alanine, glucose | glycine / alanine, fructose |
|---|---:|---:|
| methylpyrazine | 24/8 = 3.0 | 12/11 = 1.1 |
| 2,5-dimethylpyrazine | 30/25 = 1.2 | 33/35 = 0.94 |
| 2,6-dimethylpyrazine | 9/4 = 2.3 | 11/6 = 1.8 |
| trimethylpyrazine (+ co-eluate) | 35/12 = 2.9 | 42/7 = 6.0 |
| 3-ethyl-2,5-dimethylpyrazine | 0/20 = 0 | 0/19 = 0 |

Three things refuse it. (i) The same ratio moves by up to a factor of 2.7 (methylpyrazine 3.0 vs
1.1; trimethylpyrazine 2.9 vs 6.0) when only the sugar is changed, which a genuine amino-acid
identity ratio should not do; the normalisation is doing the moving. (ii) The corpus already
measures how badly the assumption fails: Leahy & Reineccius 1989 (`leahy1989_extraction.md`
Table I, 95 C, pH 9, glucose) has lysine over asparagine at **34.9x** for pyrazine, **6.4x** for
methylpyrazine and **2.1x** for 2,5-dimethylpyrazine, and Table III's totals span **27x** across five
sugar-amino-acid pairs (0.74 to 19.9 ppm). Amino acid identity changes the total by more than an
order of magnitude, so equating two pots' totals cannot be assumed to a factor of two, let alone to
the precision an identity layer would need. (iii) Amrani's own columns disagree about even the
direction: glycine looks 3x "better" than alanine for methylpyrazine and worse by an infinite factor
for 3-ethyl-2,5-dimethylpyrazine.

**What CAN be taken, and it is not nothing.**

1. **Two structural constraints on the trunk's own pyrazine step.** The zero label on
   methylpyrazine and 2,5-dimethylpyrazine in all four pots is a direct confirmation of the topology
   `parameters_pyrazine.py` declares — the ring and its methyls come from the dicarbonyl, the amino
   acid gives the nitrogen. The 25-30 % glycine label on 2,6-dimethylpyrazine is the matching
   refutation of that topology for the 2,6 isomer, which the trunk does not carry and which
   Leahy 1989 shows is the majority dimethylpyrazine in asparagine + fructose (32.3 % of total).
   A future wave that adds 2,6-dimethylpyrazine cannot make it by AKG + AKM condensation alone.
2. **The first measured branch fractions for the aldehyde-addition step** that
   `kinetic_core_b19_prereg_draft.md` section 2 wants: of the 3-ethyl-2,5-dimethylpyrazine formed,
   **100 %** came through acetaldehyde addition; of 2-ethyl-3,5-dimethylpyrazine, **70 %**; of
   2-ethyl-5-methylpyrazine, **70 %**; of 2-ethyl-6-methylpyrazine, **30 / 20 %**; of
   2,6-dimethylpyrazine from glycine's formaldehyde, **25-30 %**. Each is a partition between the
   aldehyde route and everything else, inside one pot, immune to response factors and to the missing
   total. **These are the numbers to fit an aldehyde-addition step against**, and they are the only
   ones of their kind on disk.
3. **An unfakeable qualitative identity statement:** the ethyl-bearing pyrazines are 39-51 % of the
   slate with alanine and below 2 % with glycine, and compounds 8-11 are exactly zero with glycine.
   `k3_final_parameter_inventory.md` line 1350 already declares this the star HOLD-OUT for the
   pyrazine lane, and this dossier does not change that: a fitted continuous identity ratio cannot
   reproduce an on/off switch, so it remains a test rather than a fit row.

### Against the repository's fitted pyrazine step

`src/kinetic_core/parameters_pyrazine.py` `FROZEN_B18` holds glycine's two Strecker constants,
`log10 k_go_ak = -6.5415` and `log10 k_mgo_ak = -7.5297` (L/(mmol*min) at 100 C, reference pH 6.8)
with barriers 103.1 and 114.9 kJ/mol, fitted on Zhou 2024's **fed-dicarbonyl** rates with alanine and
transferred to glycine under a declared +/- 0.5 dex band. **Amrani-Hemaimi cannot constrain any of
those four numbers**: one temperature, one time, no concentration, no rate. It touches the module in
three other places. (a) It confirms the structural claim in the module's docstring for the two
compounds the module makes (item 1 above). (b) It bears on `K_COND`, the shared "declared fast"
condensation constant, only negatively — the fact that all three pairings' products appear in one pot
is consistent with a shared constant but does not measure it. (c) It is silent on the statistical
mixed-pyrazine rule `rate_MPZ = 2 sqrt(rate_PZ rate_DMP)`, because **parent pyrazine is not among
the eleven compounds** at all, so the rule's two anchors are not both present. Flag 3 says why that
absence is itself informative.

### Against Leahy 1989's amino-acid ranking

`leahy1989_extraction.md` and this paper agree on the one thing they both measure — that amino-acid
identity restructures the pyrazine slate, not merely its size — and they cannot be combined:

| | Leahy & Reineccius 1989 | Amrani-Hemaimi 1995 |
|---|---|---|
| amino acids | lysine, asparagine (and a cysteine null) | alanine, glycine |
| sugars | glucose, fructose, ribose | glucose, fructose |
| pot | 0.1 M + 0.1 M, 0.1 M borate pH 9.0, aqueous, capped tube | 0.4 M + 0.4 M, 0.07 M phosphate pH 5.6, on kieselguhr, drying |
| temperature x time | 75 / 85 / 95 C, up to 24 h | 180 C, 7 min |
| quantification | absolute, internal standard, empirical response factors, ppm | none; within-pot NPD area % |
| compounds resolved | pyrazine, methyl-, three dimethyl- | eleven, including all the ethyl- and trimethyl- species |
| what it gives | rates (ppm/h), 3-point Ea, a total (0.74-19.9 ppm), amino-acid yield ratios | slate shares and carbon-origin fractions |

Two substantive comparisons survive the mismatch. **Parent pyrazine:** 36-56 % of the total in
Leahy's lysine pots, and not reported at all by Amrani. **The dimethylpyrazine isomer split:** Leahy's
2,5 / 2,6 ratio is effectively infinite in the lysine systems (2,6 not detected except with ribose at
0.3 %) and 3.6 in asparagine + glucose and 1.16 in asparagine + fructose; Amrani's is 6.3 and 5.8
(alanine) and 3.3 and 3.0 (glycine). So on the one axis both papers measure, the amino acid moves
the 2,5 / 2,6 split by a factor of a few in the same direction in both — the small, primary-amine
amino acids favour 2,6 relative to lysine. That is a directional cross-check, not a number to carry.

## 5. Flags

1. **There is no absolute quantification of anything, and no total.** No internal standard for
   quantification, no response factor, no concentration, no mass, no mol %, one experiment per
   condition and no replicates or error bars. Any downstream row that carries an Amrani-Hemaimi
   number in ppm, µg or mmol/L has been fabricated in transit — the same warning
   `k3_final_parameter_inventory.md` C.15 already carries. **The per-amino-acid identity ratio the
   two open pre-registrations ask for is not in this paper**, and section 4 gives the arithmetic of
   what it would take and why it fails.
2. **Table 1's percentages assume the NPD responds equally to eleven different alkylpyrazines.**
   The paper does not say so and does not correct for it. NPD response scales roughly with nitrogen
   count, which is 2 for all eleven, so the assumption is more defensible here than an FID
   normalisation would be — but it is still an assumption of the reader's, not a measurement, and
   the isolation train (ether extraction, 0.1 Pa distillation, acid/base partition, two
   concentration steps) has compound-dependent recovery that nothing corrects. Ratios between
   compounds of very different volatility (methylpyrazine RI 1160 against diethylmethylpyrazine
   RI 1505) inherit that.
3. **Parent pyrazine is absent from the compound list entirely** — the eleven begin at
   methylpyrazine. Leahy 1989 makes 36-56 % parent pyrazine from lysine + glucose at 95 C, and
   Zhou 2024's fed-glyoxal pot is the trunk's parent-pyrazine anchor. Whether Amrani's pots made
   none, or made it below detection, or lost it in the isolation, or the authors chose not to report
   it, is not stated. Do not read the absence as a measured zero.
4. **Row 6+7 is an unresolved co-elution the authors declined to resolve.** "The pyrazines 6 and 7
   were not separated on DB-Wax and showed therefore the same retention index (1400) ... Since our
   attention was more directed toward the ethyl-substituted pyrazines, we did not try to verify this
   point." The 35 % and 42 % glycine-column values are therefore **not** trimethylpyrazine's share,
   and Table 2's 50 / 40 / 80 / 100 labelling figures for that row are, by footnote c, "the average
   of both pyrazines". The authors' own guess — "the percentage of 6 plus 7 in experiments 3 and 4
   is probably mainly due to 6" — is a guess and is labelled as one. `k3` flag 28 already records
   this; it is repeated here because trimethylpyrazine IS in the registry and this is the most
   quotable-looking number in the paper.
5. **A single internally awkward cell:** 2-ethyl-5-methylpyrazine is 30 % 13C-labelled in the
   glycine + glucose pot (Table 2 row 4, column 7) while Table 1 puts that compound at <1 % of the
   glycine + glucose slate. Both are printed and raster-confirmed. A labelling fraction measured on
   a sub-percent peak by CI ion intensities carries an uncertainty the paper does not state; and
   [2-13C]glycine's carbon entering the C-2 of an *ethyl* group has no route in any of the paper's
   own schemes (Figure 3 makes that ethyl from alanine's acetaldehyde). Treat the cell as
   unexplained and do not fit it.
6. **No uncertainty anywhere on the labelling fractions.** They are quoted at 0, 20, 25, 30, 40, 50,
   70, 80, 90, 100 — a grid of about 5-10 points, which is probably the method's real resolution.
   A fit against them should carry a sigma of at least +/- 10 percentage points, declared, and
   should not treat 100 % and 90 % as distinguishable from 95 %.
7. **The pot is not aqueous and its temperature is not the sample's.** 5 mL of buffer on 5 g of
   kieselguhr, dried during a 7-minute exposure to a 180 C oil bath: a_w falls from 1 to
   something unmeasured, the concentration rises without bound as the water leaves, and the 180 C is
   the bath. `parameters_pyrazine.py` stores `AW_OF_MEASUREMENT` and a pH of measurement on every
   constant; nothing from this paper can be given either. **No rate transfer is licensed in any
   direction.**
8. **Initial pH 5.6 is printed; final pH is not.** Maillard pots acidify, and this one has 70 mmol/L
   of phosphate against 400 mmol/L of reactants, so the buffer is outmatched by nearly sixfold.
   `parameters_pyrazine.py`'s pH term (knot at 7, slopes 0.197 and 0.580 decades per unit) would put
   these pots about 0.7 decades below the trunk reference before the drift; that arithmetic is
   pointless here because there is no rate to scale, and it is recorded only so no one is tempted.
9. **Registry gaps against `data/keys/compounds.yml`:** six of the eleven compounds have no id —
   2-ethyl-5-methylpyrazine, 2-ethyl-6-methylpyrazine, 2-ethyl-3-methylpyrazine,
   3-ethyl-2,5-dimethylpyrazine, 2,3-diethyl-5-methylpyrazine, 3,5-diethyl-2-methylpyrazine. Five of
   those six are exactly the ethyl-bearing compounds that carry the alanine on/off switch and the
   aldehyde-addition branch fractions, so the repository currently cannot even name the evidence.
   Note also that the registry HAS `2_ethylpyrazine` and `tetramethylpyrazine`, neither of which
   this paper reports. Alanine, glycine, glucose, fructose, acetaldehyde and formaldehyde are not in
   `compounds.yml`.
10. **What the paper does NOT contain, listed so a future wave stops looking:** no rate constant of
    any kind; no activation energy; no second temperature; no time course (one 7-minute point); no
    dicarbonyl measurement (glyoxal and 2-oxopropanal appear only in the schemes); no aminoketone or
    dihydropyrazine measurement (both are "postulated as intermediates" in the abstract's own
    words); no Strecker aldehyde measurement (acetaldehyde and formaldehyde are inferred from the
    label positions, never quantified); no amino-acid or sugar loss; no browning; no pH after
    heating; no water activity; no yield of the eleven pyrazines individually or together.
11. **What to request.** (i) The same eight pots with an internal standard and response factors —
    that single addition would convert Table 1 into the amino-acid identity ratios both open
    pre-registrations need, and it is the one experiment that closes the gap. (ii) A fed-dicarbonyl
    version (glyoxal or methylglyoxal + alanine or glycine at a matched loading, in water) so the
    aldehyde-addition branch fractions can be tied to `FROZEN_B18`'s measured Strecker rates rather
    than to an evaporating film. (iii) Cerny & Grosch 1994 (Z. Lebensm. Unters. Forsch. 198:210),
    the companion this paper takes its apparatus from, which quantified 2-ethyl-3,5-dimethylpyrazine
    and 2,3-diethyl-5-methylpyrazine in roasted beef and may print amounts. CORRECTION (2026-09-09, the same
    day, by a reading audit): it IS on disk and read -- `data/articles/cerny1994.pdf`, dossier
    `cerny1994_extraction.md`, Z. Lebensm. Unters. Forsch. 198 (1994) 210-214, Cerny & Grosch, which quantifies
    the ethyldimethylpyrazines by isotope dilution and names alanine as their precursor. Read the two together.
