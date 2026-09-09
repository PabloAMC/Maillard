# Hwang, Hartman & Ho 1995 — EXTRACTION (glucose + 15N-glycine + one tested amino acid, equimolar, on 20 g wheat starch at 12-14 % moisture, pH 7, 180 C / 1 h; 56 pyrazines by purge-and-trap GC-MS)

### The corpus's only two-amino-acids-in-one-pot competition experiment: it prints a complete 56-compound yield table in µg per g of glucose for nine pots, but the glycine-versus-tested-amino-acid split — the identity ratio itself — exists only as percentage labels drawn on Figures 1-4.

**Source on disk:** `data/articles/hwang1995.pdf` (6 pp., owner's download, 2026-09-08). J. Agric. Food
Chem. 1995, 43, 179-184. The text layer is an OCR layer with glyph damage in prose ("Wazines",
"ly~ine-a-amine-'~N") and it **dropped the entire `ctrl` column of Table 1**. Table 1 was therefore
re-typed cell by cell from a 200-dpi raster of printed page 181 (`scratchpad/img/h95-3.png`) and
**verified by arithmetic: all nine column sums of the 56 re-typed rows reproduce the printed totals row
exactly** (598.63 / 958.30 / 2595.68 / 1347.72 / 1793.10 / 1185.42 / 923.98 / 1170.01 / 500.43), and the
row count is 56, as the abstract states. Figures 1-4 (bar charts with percentage labels) are
FIGURE-ONLY; their labels are reported in section 3 under an explicit warning and are NOT typed as
numbers in section 4.

## 0. Identity

| field | value |
|---|---|
| Title | "Relative Reactivities of Amino Acids in Pyrazine Formation" |
| Authors | Hui-Ing Hwang, Thomas G. Hartman, Chi-Tang Ho (Dept. of Food Science and Center for Advanced Food Technology, Rutgers, New Brunswick NJ) |
| Venue | J. Agric. Food Chem. 1995, 43, 179-184; received 21 June 1994, accepted 21 October 1994 |
| Article id | JF940332L (no DOI printed) |
| Naming | ctrl = glycine only; Gln / Lys / Asn / Phe / Glu / Asp / Ile / Arg each mean **labelled glycine PLUS that amino acid** in the same pot. "Ref" in Figures 1-4 is the same pot as "ctrl" in Table 1 |
| Companion | `hwang1995b_extraction.md` — same pots, same nine systems, the pyridines / pyrroles / oxazoles from the identical experiment (JAFC 1995, 43, 2917-2921) |
| Cited by | `martin2001_extraction.md` cites this paper as its ref. 33 for mixed-amino-acid pyrazine yields |

## 1. Why it matters

`results/validation/kinetic_core_b19_prereg_draft.md` section 5 and
`results/validation/kinetic_core_b22_prereg.md` section 6 both refuse a per-amino-acid Strecker rate:
no such rate exists in the corpus, and the only honest structure left is an **identity-ratio layer on
glycine's fitted step** — the two constants of `FROZEN_B18` in `src/kinetic_core/parameters_pyrazine.py`
(`log10_k_go_ak_100C` -6.542 with 103.1 kJ/mol, `log10_k_mgo_ak_100C` -7.530 with 114.9 kJ/mol) anchor
the rate, and every other amino acid enters as a partition ratio of the same dicarbonyl pool, fitted on
**within-study ratios of one amino acid against another in the same pot**.

This paper is the ideal shape for that layer and the wrong medium for it. Its shape is exactly right:
one pot, equimolar glycine and a competitor, a 15N label that assigns each pyrazine nitrogen to one of
the two, nine pots sharing a single reference. Its medium is wrong in two ways that section 4 states
plainly: the reference amino acid is glycine, which is what the repository wants, but the split is
figure-borne, and the pot is a 180 C low-moisture starch bed, not the aqueous pH 5-9 regime in which
`FROZEN_B18` was fitted. What IS printed is the pot-level yield of every pyrazine, which gives a
different and weaker within-study ratio: how much more (or less) total pyrazine a pot makes when a
second amino acid joins glycine.

## 2. Methods as they matter to a model

- **Charge (verbatim):** "Twenty grams of wheat starch, as well as an equal mole (**2.66 µmol of each**)
  of glucose, L-glycine-α-amine-15N, and the tested amino acid ... were mixed with **150 mL of deionized
  water and adjusted to pH 7** by using hydrochloric acid or sodium hydroxide." The control pot
  ("ctrl contained glycine only") therefore carries glucose + glycine and **half the total amine** of
  every other pot. See flag 1: 2.66 µmol of glucose is 0.479 mg in 20 g of starch, which does not
  square with the reported yields; 2.66 mmol does.
- **Concentrations, both readings (mine).** As printed (2.66 µmol): 17.7 µmol/L in the 150 mL steeping
  water, and **0.116 mmol/kg** of each in the ~23 g rehydrated bed. On the 2.66 mmol reading:
  17.7 mmol/L in water and **~116 mmol/kg** in the bed (glucose ~2.1 % w/w). Neither is printed as a
  concentration by the authors.
- **Drying and rehydration:** freeze-dried, then held over 20 mL of water in a desiccator to bring the
  moisture back to **12-14 %** (AOAC air-oven method). This is a low-moisture solid, not a solution: no
  buffer, and the pH 7 is the pH of the pre-drying solution, not of the reacting bed.
- **Heating:** the solid was "transferred into a reaction vessel and **heated at 180 °C for 1 h**". The
  vessel is not described (open or closed, headspace, geometry). One time point, one temperature; no
  ladder of either.
- **Isolation:** 2 g of the heated sample packed between silanized glass wool in a glass tube; **1 µL of
  1.001 mg/mL deuterated toluene** (= 1.001 µg) spiked in as internal standard; sealed into an SIS solid
  sample purge-and-trap; purged with **nitrogen at 40 mL/min at 80 °C for 1 h** onto Tenax-TA +
  Carbotrap desorption tubes.
- **Quantification:** GC-MS after Hwang et al. (1993). Linear retention indices against a C5-C25
  n-paraffin standard; identification by the NIST library or published literature. **No response factors
  are mentioned anywhere**, no calibration curve, no LOD, no recovery, no replicate count and no error
  bar on any number in Table 1. The reported unit is **µg per g of glucose**, i.e. the peak areas were
  converted against the single toluene-d8 internal standard and then divided by the glucose charge.
  Treat the absolute values as internal-standard-normalised amounts with an assumed unit response
  (flag 3), and the ratios within a column or between columns as the trustworthy part.
- **The 15N bookkeeping (verbatim):** each pyrazine can be W1 (two 14N in the ring, both from the tested
  amino acid), W2 (one 14N and one 15N) or W3 (two 15N, both from labelled glycine). The three W's are
  solved from a three-equation system in the M-1, M, M+1, M+2 abundances of the labelled and unlabelled
  runs, then

      % contribution of tested amino acid = [(W1 + ½W2)/(W1 + W2 + W3)] x 100 %
      % contribution of labelled glycine  = [(W3 + ½W2)/(W1 + W2 + W3)] x 100 %

  The abundances behind this are the paper's **supplementary material (13 pages), which is not on disk**
  (flag 5). Only the solved percentages survive, and only as figure labels.
- **Naming caution:** the label is on the **α-amine of glycine only**. Any nitrogen in the tested amino
  acid — α-amino AND side chain — counts as "tested amino acid". Lysine, arginine, asparagine and
  glutamine bring two nitrogens to the pot against glycine's one (flag 4).

## 3. Tables re-typed

### Table 1. "Pyrazines Identified in the Reaction of Glucose, Glycine-α-amine-15N, and Tested Amino Acids" — yield (µg/g of glucose)

Column meanings from the printed footnote: ctrl contained glycine only; Gln contained labeled glycine
and glutamine; Lys labeled glycine and lysine; Asn labeled glycine and asparagine; Phe labeled glycine
and phenylalanine; Glu labeled glycine and glutamic acid; Asp labeled glycine and aspartic acid; Ile
labeled glycine and isoleucine; Arg labeled glycine and arginine. "–" = not observed.

| compound | ctrl | Gln | Lys | Asn | Phe | Glu | Asp | Ile | Arg |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pyrazine | 6.71 | 41.99 | 99.94 | 78.85 | 31.65 | 75.09 | 25.91 | 8.43 | 11.49 |
| methylpyrazine | 146.77 | 279.85 | 1253.43 | 226.53 | 366.44 | 253.66 | 147.92 | 148.95 | 61.69 |
| 2,5(6)-dimethylpyrazine | 132.31 | 115.47 | 162.71 | 260.33 | 571.07 | 249.77 | 112.30 | 140.69 | 99.47 |
| 2,3-dimethylpyrazine | 45.68 | 60.66 | 274.89 | 79.64 | 162.51 | 81.29 | 61.86 | 43.28 | 33.68 |
| ethylpyrazine | – | – | – | 14.25 | – | – | – | – | 16.59 |
| vinylpyrazine | 1.23 | 12.19 | 10.04 | 9.48 | 9.51 | 2.79 | 6.04 | 2.95 | 0.89 |
| 2-ethyl-6-methylpyrazine | 39.01 | 49.23 | 66.58 | 84.83 | 93.15 | 80.54 | 49.00 | 22.08 | 8.29 |
| 2-ethyl-5-methylpyrazine | – | 17.92 | 66.58 | 14.25 | 103.96 | 125.49 | 30.10 | 22.08 | 8.29 |
| trimethylpyrazine | 133.57 | – | 281.27 | 148.48 | 207.91 | 41.49 | 41.56 | 105.54 | 73.43 |
| propylpyrazine | – | 14.81 | 12.69 | 9.18 | 7.52 | – | 8.73 | – | 2.64 |
| isopropylpyrazine | – | – | – | 0.96 | – | – | 1.39 | – | – |
| 2-vinyl-5-methylpyrazine | 5.74 | 23.29 | 20.63 | 25.44 | 10.04 | 12.60 | 14.79 | 6.81 | 2.03 |
| 2-vinyl-6-methylpyrazine | 2.79 | 9.59 | 12.69 | 9.73 | 7.52 | 7.27 | 11.05 | 5.50 | 1.16 |
| 2,6-diethylpyrazine | – | – | – | – | – | – | 40.81 | – | – |
| 3-ethyl-2,5-dimethylpyrazine | 10.40 | 34.54 | 28.79 | 53.75 | 39.29 | 36.76 | 20.40 | 20.91 | 13.25 |
| 2-ethyl-3,5-dimethylpyrazine | – | 23.33 | – | – | – | 79.41 | 63.49 | 33.38 | 4.11 |
| 2-methyl-6-propylpyrazine | – | 8.53 | – | – | – | 8.67 | 1.66 | – | – |
| tetramethylpyrazine | 61.21 | 17.70 | 80.36 | 71.28 | 85.62 | – | 28.33 | 33.38 | 38.81 |
| dimethylvinylpyrazine | – | 8.53 | – | 11.66 | – | – | 19.12 | – | 2.14 |
| 2-methyl-6-(1-propenyl)pyrazine | 3.89 | 4.33 | – | 0.96 | – | 2.23 | 1.66 | 6.37 | 2.14 |
| 2-methyl-5-(1-propenyl)pyrazine | 1.64 | 11.89 | – | 5.06 | – | 3.42 | 2.73 | – | 2.29 |
| 3,5-diethyl-2-methylpyrazine | 7.68 | 17.72 | 24.46 | 83.60 | 10.12 | 40.48 | 33.41 | 59.37 | 23.13 |
| 2,5-diethyl-3-methylpyrazine | – | 10.91 | 33.22 | 8.63 | 3.98 | 9.70 | 19.43 | 6.82 | 5.19 |
| 2,3-diethyl-5-methylpyrazine | – | 10.27 | 7.89 | 7.29 | 11.29 | 4.85 | 41.55 | 4.63 | 2.29 |
| 2,5-dimethyl-3-propylpyrazine | – | 7.70 | – | – | – | 2.07 | – | – | – |
| 2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 27.76 | – |
| 2,6-diethyl-3,5-dimethylpyrazine | – | 5.06 | – | 8.54 | – | 3.48 | 18.71 | – | 2.68 |
| 2,5-diethyl-3,6-dimethylpyrazine | – | – | – | 2.39 | – | – | 20.71 | – | – |
| 2,3-dimethyl-5-(1-methylpropyl)pyrazine | – | 1.71 | – | 2.32 | – | – | 10.72 | – | – |
| 5-methyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 51.99 | – |
| 6-methyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 39.87 | – |
| 3-methyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 36.27 | – |
| ethyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 53.08 | – |
| dimethyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 29.90 | – |
| dimethyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 29.90 | – |
| ethyl-2-(2-methylbutyl)pyrazine | – | – | – | – | – | – | – | 37.93 | – |
| C8-alkylpyrazine | – | – | – | – | – | – | – | 50.01 | – |
| C8-alkylpyrazine | – | – | – | – | – | – | – | 53.41 | – |
| C9-alkylpyrazine | – | – | – | – | – | – | – | 19.75 | – |
| 5-methyl-5H-cyclopentapyrazine | – | 10.91 | 15.74 | – | – | – | – | – | – |
| 5-methyl-6,7-dihydro-5H-cyclopentapyrazine | – | 7.51 | 28.48 | 17.29 | 18.45 | 9.38 | 14.19 | 13.60 | 11.62 |
| 2-ethyl-6,7-dihydro-5H-cyclopentapyrazine | – | – | – | 3.23 | – | – | 6.89 | 6.37 | 1.82 |
| 5,7-dimethyl-6,7-dihydro-5H-cyclopentapyrazine | – | 17.86 | 22.74 | 21.83 | – | 10.26 | 12.09 | 21.39 | 19.02 |
| 2,6-dimethyl-6,7-dihydro-5H-cyclopentapyrazine | – | 4.22 | – | 2.39 | – | – | 6.16 | – | 1.32 |
| 2-ethyl-5-methyl-6,7-dihydro-5H-cyclopentapyrazine | – | 4.93 | – | 5.19 | – | 1.43 | – | – | – |
| 2-(2'-furanyl)pyrazine | – | 21.18 | 30.5 | 14.02 | 4.26 | 1.81 | 4.15 | – | 10.94 |
| 2-(2'-furanyl)-6-methylpyrazine | – | 24.79 | 18.34 | 17.35 | 6.01 | 1.87 | 3.89 | – | 12.96 |
| 2-(2'-furanyl)-5-methylpyrazine | – | 24.79 | 18.34 | 8.68 | 6.01 | 1.87 | 2.82 | – | 8.64 |
| 2-(2'-furanyl)-5,6-dimethylpyrazine | – | 11.16 | 8.43 | 5.58 | 1.75 | 2.09 | 2.89 | – | 5.73 |
| 2-(2'-furanyl)-methylethylpyrazine | – | 2.01 | – | – | – | – | 4.82 | – | 2.34 |
| 2-(2'-phenylethyl)-5(6)-methylpyrazine | – | – | – | – | 4.85 | – | – | – | – |
| 2-(2'-phenylethyl)-3,5(6)-dimethylpyrazine | – | – | – | – | 20.47 | – | – | – | – |
| 6-methyl-2-acetylpyrazine | – | 18.99 | 8.47 | 14.15 | 4.86 | 15.27 | 15.14 | – | 5.87 |
| 5-methyl-2-acetylpyrazine | – | 14.19 | 8.47 | 13.29 | 4.86 | 7.14 | 10.10 | – | 1.38 |
| 3,5(6)-dimethyl-2-acetylpyrazine | – | 5.04 | – | – | – | 5.22 | 1.78 | – | – |
| 3,5,6-trimethyl-2-acetylpyrazine | – | 3.50 | – | 7.29 | – | 8.02 | 5.68 | 27.61 | 3.11 |
| **totals** | **598.63** | **958.30** | **2595.68** | **1347.72** | **1793.10** | **1185.42** | **923.98** | **1170.01** | **500.43** |

56 compound rows, as the abstract states. Every column sum reproduces the printed total exactly (mine).

### Table 2. "Mass Spectra of Furanylpyrazines and Acetylpyrazines Tentatively Identified..."

Only the MW and retention-index columns survive the text layer; the m/z (relative intensity) column is
empty in the extracted text and was not needed here, so it was not rastered. Values printed:

| compound | MW | RI (DB-1, C5-C25 n-paraffins) |
|---|---:|---:|
| 2-(2'-furanyl)pyrazine | 146 | 1476 |
| 2-(2'-furanyl)-6-methylpyrazine | 160 | 1728 |
| 2-(2'-furanyl)-5-methylpyrazine | 160 | 1732 |
| 2-(2'-furanyl)-5,6-dimethylpyrazine | 174 | 2024 |
| 2-(2'-furanyl)-methylethylpyrazine | 188 | 2214 |
| 6-methyl-2-acetylpyrazine | 136 | 1141 |
| 5-methyl-2-acetylpyrazine | 136 | 1149 |
| 3,5(6)-dimethyl-2-acetylpyrazine | 150 | 1328 |
| 3,5,6-trimethyl-2-acetylpyrazine | 164 | 1587 |

The mass-spectral intensities themselves are garbled in the text layer and are not transcribed
(identification data, of no kinetic use).

### Percentage labels drawn on Figures 1 and 2 — FIGURE_ONLY, NOT to be fitted

Read from 200-dpi rasters of printed page 182 (`scratchpad/img/h95f-4.png`). These are printed digits
placed above the bars, not values interpolated off an axis; they are recorded here so the information is
not lost, but under the house rule they are figure-borne and section 4 does not carry them as numbers.
The caption defines them: "The numbers on the tops of columns show the percent contributions of tested
amino acids", i.e. the 15N-solved share of pyrazine ring nitrogen coming from the tested amino acid
rather than from glycine.

| pot | Figure 1 (all 56 pyrazines) | Figure 2 (the 28 alkylpyrazines) |
|---|---:|---:|
| Gln | 55 % | 55 % |
| Glu | 50 % | 56 % |
| Asn | 76 % | 77 % |
| Asp | 65 % | 67 % |
| Lys | 64 % | 66 % |
| Arg | 68 % | 70 % |
| Phe | 71 % | 72 % |
| Ile | 70 % | 72 % |

The bar heights in Figures 1 and 2 are the Table 1 totals (Ref ~600, Lys ~2600, Phe ~1800, Arg ~500),
which is how the mapping above was checked. Figures 3 and 4 carry the same kind of labels for the
bicyclic and acetylpyrazine subsets and were not rastered; nothing in this dossier depends on them.

### Statements printed in prose

- "In the presence of glycine, glutamine and glutamic acid are the smallest contributors, while
  asparagine is the highest contributor to pyrazine formation among these tested amino acids."
- "the yield of pyrazines from the reaction mixture containing lysine is the highest, whereas the amount
  of pyrazines from the reaction mixture containing arginine is the lowest" — and glycine's own yield is
  highest in the lysine pot and lowest in the arginine pot: "lysine acts as a synergist to increase the
  reactivity of other amino acids (glycine in this case) ... arginine acts like an inhibitor".
- "the yield of pyrazines from the reaction mixture containing phenylalanine is the second highest."
- "Twenty-eight alkylpyrazines were detected in this study."
- Amino-acid-specific products: two phenylethylpyrazines from phenylalanine (via phenylacetaldehyde);
  ten 2-methylbutyl-substituted pyrazines from isoleucine (via 2-methylbutanal). The isoleucine pot was
  the only one with "a roasted cocoa-like flavor".
- Mechanistic assertion: "The most direct route for their formation results from the interaction of
  α-dicarbonyls and amines through Strecker degradation" — the same step as `FROZEN_B18`.

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`): `pyrazines`, `methylpyrazine`, `2_3_dimethylpyrazine`,
`2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`, `2_ethylpyrazine`, `trimethylpyrazine`,
`tetramethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine` are present. **Vinylpyrazine, the ethyl-methyl and
diethyl pyrazines, the propyl/propenyl pyrazines, all six cyclopentapyrazines, all five
furanylpyrazines, both phenylethylpyrazines, all four acetylpyrazines and every 2-methylbutylpyrazine
have no id.** No amino acid has an id (the rules file uses short names).

**The pot for every row below:** glucose + 15N-glycine + one tested amino acid, equimolar
(2.66 µmol each as printed — see flag 1), on 20 g of wheat starch rehydrated to 12-14 % moisture, the
pre-drying solution adjusted to pH 7, **180 C for 1 h**, one time point, no replicates stated. The
control pot carries glycine alone and therefore **half the total amine**.

### The identity ratios this paper actually supports

| ratio | value | unit | what it is a ratio OF | conditions | source location | evidence class |
|---|---:|---|---|---|---|---|
| total pyrazine, Gly+Gln pot ÷ Gly-only pot | 1.601 | – | a **yield** (µg/g glucose), pot against pot | as above | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Lys ÷ Gly-only | 4.336 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Asn ÷ Gly-only | 2.251 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Phe ÷ Gly-only | 2.995 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Glu ÷ Gly-only | 1.980 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Asp ÷ Gly-only | 1.543 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Ile ÷ Gly-only | 1.954 | – | yield | " | Table 1 totals (mine) | within_study_ratio |
| total pyrazine, Gly+Arg ÷ Gly-only | 0.836 | – | yield (the only pot BELOW the control) | " | Table 1 totals (mine) | within_study_ratio |
| the same eight ratios, halved to correct for the doubled amine charge | 0.800 / 2.168 / 1.126 / 1.498 / 0.990 / 0.772 / 0.977 / 0.418 (Gln / Lys / Asn / Phe / Glu / Asp / Ile / Arg) | – | yield per unit total amine | " | Table 1 totals (mine) | derived_assumption (assumes yield is first order in total amine, which the lysine synergy contradicts — flag 6) |
| **per-amino-acid share of pyrazine ring nitrogen (tested AA vs glycine)** | **not typed — see §3, Figures 1 and 2** | % | the quantity the identity layer actually wants | " | Figures 1-4 labels | **figure_only** |
| individual compound yields, all 56 x 9 cells | see Table 1 | µg per g of glucose | absolute amount, internal-standard-normalised | " | Table 1 | level_only (no response factors; flag 3) |
| molar yield of the Gly-only pot, total pyrazines | 0.10 | mol % of glucose (nominal, at MW 108.14) | – | " | Table 1 totals (mine) | derived_assumption |
| molar yield, methylpyrazine, Gly+Lys pot | 0.24 | mol % of glucose | – | " | Table 1 (mine) | derived_assumption |
| molar yield, 2,5(6)-dimethylpyrazine, Gly+Phe pot | 0.095 | mol % of glucose | – | " | Table 1 (mine) | derived_assumption |
| lysine synergy on glycine | glycine's own yield is highest in the lysine pot; arginine's the lowest | – | a **direction**, not a number | " | text, p. 181-182 | level_only |
| ordering of pot totals | Lys > Phe > Asn > Glu > Ile > Asp > Gln > ctrl > Arg | – | ordinal | " | Table 1 totals | within_study_ratio (ordinal) |

**Molar-yield arithmetic (mine).** 1 g of glucose is 5.551 mmol; a yield of Y µg per g of glucose at
molar mass M is 100·(Y·10⁻⁶/M)/5.551·10⁻³ mol % of glucose. Worked for methylpyrazine (M = 94.11) in the
lysine pot: 1253.43 µg/g → 1.332·10⁻⁵ mol per g glucose → **0.240 mol %**. The mol % column is
independent of whether the charge was 2.66 µmol or 2.66 mmol, because the denominator is the glucose
charge itself; that is the one quantity flag 1 does not damage.

### What this does and does not give the identity layer

**The reference amino acid is glycine**, in every one of the nine pots, which is exactly the anchor
`FROZEN_B18` provides. **The span of the printed ratios is 0.84 to 4.34** on the pot totals (a factor of
5.2 from arginine to lysine), or 0.42 to 2.17 once the doubled amine charge is divided out. **But the
printed ratio is the wrong ratio.** A pot total answers "how much pyrazine does glycine plus X make
against glycine alone" — it mixes the competitor's own reactivity, its effect on glycine's reactivity
(the paper's own headline: lysine catalyses, arginine inhibits), and its effect on sugar fragmentation,
into one number. The quantity a partition layer needs — X's share of the ring nitrogen against
glycine's, in the same pot, from the 15N solve — is printed **only as labels on Figures 1-4**, and under
the house rule a figure-borne number is not typed as a value. Section 3 records those labels so the size
of the prize is visible (Asn 76 %, Glu 50 %, i.e. partition ratios against glycine of roughly 3.2 down
to 1.0 if they were usable); flag 5 says how to make them usable.

**Two further reasons this paper cannot carry the layer alone.** (i) The pot is a **180 C, 12-14 %
moisture starch bed**, and `FROZEN_B18` was fitted on fed dicarbonyls in water at 100-120 C, initial
pH 8, with a pH term knotted at 7 (`PYRAZINE_PH_SLOPES`). Nothing here is at the trunk's temperature,
water activity or pH, and there is no ladder in any of the three to bridge. A ratio fitted here is a
dry-bed ratio, and `PYRAZINE_SUPPLY_CAVEAT` already records that the dicarbonyl supply, not the Strecker
step, is what moves between media. (ii) The nitrogen count is asymmetric: glycine has one nitrogen,
lysine, arginine, asparagine and glutamine two, and the paper's own equations assign **all** the tested
amino acid's nitrogen to it (flag 4), so even the figure percentages are not a clean per-α-amino-group
partition.

## 5. Flags

1. **The charge is internally inconsistent by three orders of magnitude.** "2.66 µmol of each" of
   glucose is 0.479 mg spread through 20 g of starch. On that reading the whole batch would contain
   0.29 µg of total pyrazines (598.63 µg/g × 4.79·10⁻⁴ g), of which the 2 g analysed would hold 29 ng,
   against a 1.001 µg internal standard — 56 compounds could not be identified at that level. On the
   2.66 **mmol** reading the numbers are ordinary (glucose ~2 % w/w of the bed; ~29 µg of pyrazines in
   the analysed portion). The companion `hwang1995b_extraction.md` prints the identical "2.66 µmol".
   The repository must not choose: the paper prints µmol, and every absolute concentration derived from
   it is therefore unsafe. The **µg per g of glucose** basis and every ratio within it are unaffected.
2. **The control pot is not the same pot.** It carries glycine only, so it has half the amine and none
   of the competitor. A ratio of pot totals to it is not a per-amino-acid reactivity.
3. **No response factors, no calibration, no LOD, no replicates, no error bars.** One deuterated toluene
   internal standard for 56 analytes spanning MW 80 to ~190, quantified as "µg/g of glucose". Comparing
   a furanylpyrazine with pyrazine across a row is comparing peak responses, not moles.
4. **The 15N label is on glycine's α-amine only**, and the paper's equations count every nitrogen that
   is not that label as "tested amino acid" — including the side-chain nitrogens of lysine, arginine,
   asparagine and glutamine. The paper says so itself ("the participation of side-chain nitrogen of
   glutamine and lysine in pyrazine formation has been proved"). Any partition ratio taken from the
   percentages is per amino acid **molecule**, not per α-amino group.
5. **The supplementary material is not on disk.** "Experimental data from the reaction of glycine, the
   tested amino acid, and glucose (13 pages)" — the M-1 / M / M+1 / M+2 abundances from which W1, W2 and
   W3 were solved. **To request:** that supplement, or the solved W1/W2/W3 table, from ACS or the
   authors. With it the per-amino-acid partition ratios become printed numbers rather than figure
   labels, and the identity layer's rows exist.
6. **The pot totals are not additive and the paper knows it.** Glycine's own yield goes up in the lysine
   pot and down in the arginine pot, so a second amino acid changes the first's rate. A partition model
   in which each amine competes for a fixed dicarbonyl pool cannot produce that; a term for
   amine-catalysed sugar fragmentation can. Halving the pot ratios (row 9 of §4) assumes the additivity
   the paper refutes and is marked derived_assumption for that reason.
7. **One temperature, one time, no kinetics.** 180 C for 1 h, single point. Nothing here identifies a
   rate or a barrier, and no time course exists to tell formation from subsequent destruction — the
   authors note pyrazines can themselves be consumed.
8. **The reacting pH is unknown.** pH 7 was set in the 150 mL of water before freeze-drying; the bed
   then reacted at 12-14 % moisture. `PYRAZINE_PH_STEPS` scales the two Strecker constants by pH, and
   this pot cannot be placed on that axis.
9. **Two rows of Table 1 are printed twice with different values** ("dimethyl-2-(2-methylbutyl)pyrazine"
   29.90 twice, and "ethyl-2-(2-methylbutyl)pyrazine" 53.08 and 37.93): isomers the authors could not
   distinguish, listed under one name. They are transcribed as printed.
10. **Registry gaps:** of the 56 compounds, 9 are registered. The unregistered families that carry real
    yield here are the cyclopentapyrazines (six), the furanylpyrazines (five), the acetylpyrazines
    (four, and the paper notes 5- and 6-methyl-2-acetylpyrazine are popcorn-like) and the ten
    isoleucine-specific 2-methylbutylpyrazines.
