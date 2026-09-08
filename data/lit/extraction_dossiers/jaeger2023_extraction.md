# Jaeger et al. 2023 — EXTRACTION (commercial pea and soy protein isolates from one supplier: proximate composition, total and FREE amino acids, total free sugars and FODMAPs, pH, solubility; the brewer's-spent-grain isolate EverPro is the paper's subject and gets one line here)
### The only paper on disk that prints, for a commercial pea and a commercial soy isolate, both the total amino-acid table and a free-amino-acid table on the same powders — and finds essentially no free amino acids in either.

**Source on disk:** `data/articles/Jaeger2023.pdf` (owner's download, 2026-09-08). Read from the
scratchpad text layer (`Jaeger2023.txt`, clean; the MDPI PDF carries every page twice — a
"FOR PEER REVIEW" layer under the typeset one — so figures and paragraphs repeat in the text
file; values were taken from the typeset copy); pypdf layout mode on pp. 7-8 confirmed Table 1
and both halves of Table 2 row by row. Three tables (Table 1 p. 7, Table 2 pp. 7-8, Table 3
p. 9), seven figures. Repo status before this dossier: `data/species/protein_matrices.yml`
carries the amine pool of `pea_isolate` and `soy_isolate` as a lysine content taken from FLOUR
compositions (pea 0.524, soy 0.379 mmol per g protein, USDA via Xiao 2024 §6) and states nothing
about the free amino acids or free sugars an isolate brings; Programme 7 of
`tasks/roadmap_for_scientists.md` §5d asks for exactly those. No dossier cites this paper.

## 0. Identity

| field | value as printed |
|---|---|
| Title | "Functional Properties of Brewer's Spent Grain Protein Isolate: The Missing Piece in the Plant Protein Portfolio" |
| Authors | Alice Jaeger, Aylin W. Sahin, Laura Nyhan, Emanuele Zannini, Elke K. Arendt (corresponding) — School of Food and Nutritional Science, University College Cork; Zannini also Sapienza Rome; Arendt also APC Microbiome Institute |
| Venue | Foods 2023, 12, 798; received 13 Jan 2023, revised 9 Feb, accepted 10 Feb, published 13 Feb 2023; open access CC BY |
| DOI | 10.3390/foods12040798 |
| Products | **PPI** = pea protein isolate and **SPI** = soy protein isolate, both "obtained from Naturz Organics (Helmond, Netherlands)" — no product name, grade or lot number printed. **EverPro** = barley-rice protein isolate from brewer's spent grain, EverGrain Ingredients (St. Louis, MO, USA); rice was a brewing adjunct, so EverPro is a barley + rice protein. |
| Naming | "Pea" / "Soy" in Table 3 and the figures = PPI / SPI. "g/100 g DM" = per 100 g dry matter of the powder. "n.d." = not detected or below LoQ (< 0.025 g/100 g for the HPAEC-PAD carbohydrates). "% Requirement" = each essential amino acid per g protein over the WHO (2007) adult requirement. |
| Nitrogen factor | **6.25** for all three (Kjeldahl, AACC 46-12) |
| Compound registry | none of the paper's analytes is a volatile; lysine / arginine / free amino acids / sucrose / glucose / fructose / maltose map to the engine's reactant pools, not to the volatile registry |

## 1. Why it matters

Programme 7 says the engine charges free amino acids at tens of millimoles per litre while a
pea or soy isolate holds few free amino acids and much protein-bound lysine. This paper is the
direct measurement of that sentence for two commercial isolates: **free amino acids total
0.079 +/- 0.014 g/100 g DM in the pea isolate (arginine, glutamine, leucine only) and 0 in the
soy isolate**, against total (protein-bound) lysine of 6.399 and 5.343 g/100 g DM. In molar terms
the pea isolate brings about 0.005 mmol free amino acids per gram of powder — at 100 g isolate
per litre that is 0.5 mM, two orders below what the engine charges — and 0.44 mmol bound lysine
per gram of powder (0.54 mmol per g protein). The same table gives arginine (the second amine
pool) and the sulfur amino acids, the latter only as a methionine + cysteine sum. Table 1 gives
the free sugars as one total (0.19 / 0.06 g/100 g DM, i.e. 0.006-0.011 / 0.002-0.003 mmol per g
powder depending on whether it is sucrose or hexose) and total FODMAPs (0.53 / 0.04 g/100 g DM,
mostly the galacto-oligosaccharides of pulses). So the isolate is, in this paper's numbers, a
protein-bound-lysine reservoir with almost no free reducing sugar and almost no free amine: the
Maillard chemistry of an isolate-only recipe is the glycation of protein lysine by whatever
sugar the recipe adds, not the free-precursor chemistry the engine runs today. The lysine per
gram of protein measured here (pea 0.539, soy 0.410 mmol/g) is what would replace the flour-derived
0.524 / 0.379 in `protein_matrices.yml`. The paper also prints the pH of the powders in water
(6.36 / 6.77), their solubility at pH 7 (22 / 52 %), colour, particle size, gelation and foaming
(Table 3, figures) — matrix data, not reactant data.

## 2. Methods as they matter to a model

- **Materials.** PPI and SPI from Naturz Organics (Helmond, NL); EverPro from EverGrain. Nothing
  on the isolates' process; the discussion says "Alkaline extraction followed by iso-electric
  point precipitation ... is likely the method used to produce PPI and SPI, although the exact
  parameters are unknown", and from the Bioanalyzer profiles "it appears that the isolates are in
  their native state and have not been subjected to protein degradation."
- **Proximate.** Protein Kjeldahl (AACC 46-12) x **6.25**; fat AACC 30-25.01; moisture oven
  (AACC 44-15.02); total starch Megazyme K-RAPRS. Ash NOT measured (flag 4).
- **Sugars (verbatim core):** "sugar content was determined by HPLC, using the extraction method
  of Hoehnel et al. (2020) ... Sugars were quantified by HPLC on an Infinity 1260 system with a
  refractive index detector ... using a Sugar-Pak I column (300 mm x 6.5 mm ...), and an eluent of
  0.0001 M CaEDTA at a flow rate of 0.5 mL/min and a column temperature of 80 C. Maltotriose,
  sucrose, lactose, glucose, fructose and mannitol were used as external standards." Table 1
  prints ONE row, "Total Sugars (Sucrose, Glucose, Fructose, Maltose)"; the individual sugars are
  not printed (flag 2).
- **FODMAPs.** HPAEC-PAD (Dionex ICS-5000+), Ispiryan et al. 2019: mono-, di-saccharides,
  galacto-oligosaccharides, fructans (by enzymatic difference, mixtures A and B), polyols; LoQ
  0.025 g/100 g; "All extractions were carried out in triplicate"; results on a dry-weight basis.
  Only the total is printed.
- **Amino acids (verbatim core):** "The amino acid composition was determined externally by
  Chelab S.r.l. (Resana, Italy), using ion chromatography with post-column derivatisation with
  ninhydrin, or HPLC-UV analysis in the case of tryptophan." Total amino acids after hydrolysis
  (conditions not stated; a contract-laboratory method) and free amino acids on the same
  powders, both "on a dry matter basis" (Table 2 caption). Methionine and cysteine are reported
  only as a sum; phenylalanine and tyrosine likewise. The "% Requirement" column is each essential
  amino acid per g protein (using Table 1's protein content) over WHO 2007 — the arithmetic
  reproduces exactly (see §3), which pins the units of Table 2 as g per 100 g dry powder.
- **Replicates and uncertainty.** "All experiments were performed in triplicate"; one-way ANOVA +
  Tukey, p < 0.05, SPSS 26. The +/- on every total-amino-acid row of PPI and SPI is a fixed
  13.9 % of the value (EverPro: 10 %; Trp 10.8 %; the two summed pairs 10-12 %), so these are a
  laboratory-assigned relative uncertainty, not replicate scatter (flag 1).
- **pH and TTA.** 10 g powder + 95 mL water + 5 mL acetone, stirred; TTA to pH 8.5 with 0.1 M
  NaOH (mL).
- **Protein solubility.** 1 % (w/v) protein at pH 7, overnight 4 C, 22 C, 4893 g 30 min, Kjeldahl
  on the supernatant, % of total protein.
- **Other functional assays** (foaming 2 % w/v; fat absorption; Bioanalyzer profiles in 2 % SDS /
  2 M thiourea / 6 M urea, +/- DTT; particle size dry and 1 % w/v wet; minimum gelation 6-22 %,
  90 C 30 min; rheology 20 -> 90 C at 2 C/min, 30 min hold, cool; ANS surface hydrophobicity;
  zeta at 0.1 % w/v pH 7; emulsions 1.2 % w/v with 10 % sunflower oil, LUMiSizer; Minolta colour;
  SEM) — matrix descriptors only.

## 3. Tables re-typed

### Table 1. "Compositional analysis of PPI, SPI and EverPro." (g/100 g DM; row letters = Tukey groups within the row)

| row | PPI | SPI | EverPro |
|---|---|---|---|
| Moisture | 7.29 +/- 0.11 a | 6.62 +/- 0.08 b | 6.10 +/- 0.05 c |
| Protein | 81.22 +/- 0.43 a | 89.22 +/- 0.33 b | 83.22 +/- 1.11 c |
| Fat | 8.51 +/- 0.07 a | 1.72 +/- 0.17 b | 0.37 +/- 0.06 c |
| Total Sugars (Sucrose, Glucose, Fructose, Maltose) | 0.19 +/- 0.00 a | 0.06 +/- 0.00 b | 0.30 +/- 0.00 c |
| Total Starch (digestible) | 2.83 +/- 0.06 a | 1.67 +/- 0.03 b | 1.42 +/- 0.04 c |
| FODMAPs (Total) | 0.53 +/- 0 a | 0.04 +/- 0 b | n.d.* c |

\* n.d.: not detected or below LoQ (< 0.025 %). Protein by Kjeldahl N x 6.25. The column head
says g/100 g DM for every row including moisture (flag 3). Ash is not a row (flag 4).

### Table 2 (first half). "Complete amino acid analysis ... total amino acids quantification on a dry matter basis" (g/100 g DM of powder; "% Requirement" = per g protein vs WHO 2007)

Every row re-typed for PPI and SPI; EverPro summarised below the table.

| amino acid | PPI (g/100 g DM) | PPI % req. | SPI (g/100 g DM) | SPI % req. |
|---|---|---|---|---|
| Histidine | 2.227 +/- 0.311 | 182.8 | 2.255 +/- 0.314 | 168.5 |
| Isoleucine | 3.955 +/- 0.550 | 162.3 | 3.116 +/- 0.434 | 116.4 |
| Leucine | 7.422 +/- 1.032 | 154.9 | 6.689 +/- 0.931 | 127.1 |
| **Lysine** | **6.399 +/- 0.891** | 175.1 | **5.343 +/- 0.743** | 133.1 |
| **Methionine + cysteine** | **0.605 +/- 0.075** | 33.9 | **0.942 +/- 0.118** | 48.0 |
| Phenylalanine + tyrosine | 8.456 +/- 0.843 | 273.9 | 8.16 +/- 0.814 | 240.8 |
| Threonine | 3.137 +/- 0.437 | 167.9 | 3.281 +/- 0.456 | 159.9 |
| Tryptophan | 0.398 +/- 0.043 | 81.7 | 0.690 +/- 0.075 | 128.8 |
| Valine | 4.433 +/- 0.617 | 139.9 | 3.358 +/- 0.468 | 96.5 |
| Alanine | 3.66 +/- 0.509 | — | 3.571 +/- 0.497 | — |
| **Arginine** | **8.422 +/- 1.172** | — | **7.327 +/- 1.019** | — |
| Aspartic acid | 9.286 +/- 1.292 | — | 8.963 +/- 1.248 | — |
| Glutamic acid | 14.549 +/- 2.025 | — | 15.555 +/- 2.164 | — |
| Proline | 4.035 +/- 0.562 | — | 4.511 +/- 0.628 | — |
| Serine | 4.694 +/- 0.654 | — | 4.907 +/- 0.683 | — |
| Glycine | 3.842 +/- 0.535 | — | 3.688 +/- 0.514 | — |

Checks. Column sums: PPI 85.52 g amino acids / 100 g DM against 81.22 g protein (ratio 1.05);
SPI 82.36 against 89.22 (0.92) — both within the usual hydrolysate-recovery band, consistent with
g per 100 g dry POWDER. WHO 2007 lysine requirement 45 mg/g protein: PPI 63.99 mg/g powder /
0.8122 = 78.8 mg/g protein, 78.8 / 45 = 175.1 % (printed 175.1); SPI 53.43 / 0.8922 = 59.9,
59.9 / 45 = 133.1 % (printed 133.1). Met + Cys requirement 22 mg/g: PPI 7.45 / 22 = 33.9 %, SPI
10.56 / 22 = 48.0 % (both as printed). So the paper's own per-g-protein arithmetic uses Table 1's
protein content, and the conversions in §4 follow it. Relative +/-: 0.891/6.399 = 0.139,
1.172/8.422 = 0.139, 0.311/2.227 = 0.140, 2.025/14.549 = 0.139 ... a fixed 13.9 % (flag 1).

EverPro (one line): lysine 3.035 +/- 0.277 (81.1 % of requirement), arginine 4.707 +/- 0.469,
Met + Cys 2.769 +/- 0.296 (151.3 %), glutamic acid 20.234, proline 8.019, Trp 1.171 g/100 g DM.

### Table 2 (second half). "Free Amino Acids" (g/100 g DM of powder)

| free amino acid | PPI | SPI | EverPro |
|---|---|---|---|
| Glutamic acid | 0 | 0 | 0.168 +/- 0.028 |
| Alanine | 0 | 0 | 0.085 +/- 0.016 |
| **Arginine** | **0.053 +/- 0.009** | 0 | 0.018 +/- 0.004 |
| Asparagine | n.d. | n.d. | 0.012 +/- 0.003 |
| Citrulline | 0 | 0 | 0.044 +/- 0.010 |
| Phenylalanine | 0 | 0 | 0.132 +/- 0.022 |
| Glycine | n.d. | n.d. | 0.007 +/- 0.002 |
| **Glutamine** | **0.013 +/- 0.008** | 0 | 0.038 +/- 0.008 |
| Isoleucine | 0 | 0 | 0.086 +/- 0.016 |
| Histidine | 0 | 0 | 0.014 +/- 0.003 |
| **Leucine** | **0.013 +/- 0.008** | 0 | 0.217 +/- 0.036 |
| **Lysine** | **0** | **0** | 0.024 +/- 0.005 |
| Methionine | 0 | 0 | 0.054 +/- 0.012 |
| Ornithine | n.d. | n.d. | 0.01 +/- 0.002 |
| Serine | n.d. | n.d. | 0.039 +/- 0.009 |
| Tyrosine | n.d. | n.d. | 0.108 +/- 0.019 |
| Threonine | n.d. | n.d. | 0.026 +/- 0.006 |
| Aspartic acid | n.d. | n.d. | 0.081 +/- 0.015 |
| Valine | 0 | 0 | 0.085 +/- 0.016 |
| **Total free amino acid** | **0.079 +/- 0.014** | **0 +/- 0** | 1.248 +/- 0.066 |

The table prints both "0" and "n.d." for absent analytes without saying what distinguishes them
(flag 5). Free cysteine, free proline and free tryptophan are not rows. Text: "No free amino
acids were detected in SPI, while only low levels of free arginine, glutamine and leucine were
detected in PPI." EverPro's free amino acids (1.25 g/100 g, with free lysine 0.024) are the
brewing-hydrolysis signature the paper is about.

### Table 3. "Functional properties of pea, soy and EverPro." (PPI and SPI columns)

| row | Pea | Soy |
|---|---|---|
| pH (10 g in 95 mL water + 5 mL acetone) | 6.360 +/- 0.030 a | 6.770 +/- 0.010 b |
| TTA (mL 0.1 M NaOH to pH 8.5) | 11.470 +/- 0.020 a | 12.060 +/- 0.260 a |
| Protein solubility at pH 7 (%) | 22.267 +/- 1.457 a | 51.960 +/- 3.354 b |
| Surface hydrophobicity (a.u.) | 4292.467 +/- 500 a | 7471.367 +/- 324 b |
| Zeta potential at pH 7 (mV) | -22.600 +/- 1.633 a | -33.778 +/- 1.524 b |
| Separation rate (%/min) | 1.327 +/- 0.110 a | 0.990 +/- 0.008 a |
| Fat absorption (%) | 157.723 +/- 3.202 a | 120.050 +/- 17.847 b |
| Colour L* / a* / b* | 84.688 +/- 1.164 / 2.180 +/- 0.036 / 22.338 +/- 0.215 | 83.821 +/- 0.603 / 0.734 +/- 0.049 / 18.498 +/- 0.129 |

EverPro: pH 7.90, TTA 17.85, solubility 101.7 %, zeta -30.0 mV, L* 57.4 (dark: the paper
attributes it to Maillard browning during kilning and extraction).

### Numbers printed only in the running text (figures otherwise FIGURE-ONLY)

| item | as printed | where |
|---|---|---|
| Particle size, dry powder D[4,3] | PPI 58.63 +/- 0.058 um; SPI 56.57 +/- 1.069 um | §3.3.4 (Fig. 3) |
| Particle size, 1 % dispersion D[4,3] / D[3,2] | PPI 70.48 +/- 2.59 / 32.22 +/- 1.15 um; SPI 144.33 +/- 8.65 / 88.17 +/- 6.09 um | §3.3.4 |
| Foaming capacity / stability | PPI 38.19 +/- 1.20 % / 80.12 +/- 5.57 %; SPI 70.14 +/- 3.18 % / 74.31 +/- 5.97 % | §3.3.6 (Fig. 4) |
| Minimum gelation concentration | PPI 14 %, SPI 8 % (w/w), 90 C 30 min | §3.3.7 |
| Initial tan delta at 20 C | PPI 2.99, SPI 0.36 (SPI is a cold-set paste before heating; its G' falls during the 90 C hold) | §3.3.7 (Fig. 5) |
| Emulsion droplet D[3,2] | PPI 8.77, SPI 10.16 um; D[4,3] 22.78-27.74 um for all | §3.3.8 |
| Bioanalyzer bands | pea: ~65 kDa (convicilin), ~48 kDa (vicilin), ~40 / ~20 kDa (legumin acidic / basic); soy: ~70 / ~80 / ~50 kDa (beta-conglycinin alpha / alpha' / beta), ~40 kDa (glycinin acidic) | §3.2.2, §4 (Fig. 1) |
| Another commercial pea isolate's GOS | 1.16 g/100 g DM (Ispiryan 2020, cited) | §4 |

## 4. Numbers the repository can use

Molar masses used: Lys 146.19, Arg 174.20, Met 149.21, Cys 121.16, Gln 146.15, Leu 131.17,
glucose / fructose 180.16, sucrose / maltose 342.30 g/mol. "Powder" = dry matter of the isolate
(Table 1 basis). Per g protein = per g powder / (protein g/100 g DM / 100), the paper's own
convention (§3 check).

| product | quantity | value +/- sd | unit as printed | mmol per g protein (arithmetic) | method | source | evidence class |
|---|---|---|---|---|---|---|---|
| PPI (Naturz Organics) | protein content | 81.22 +/- 0.43 | g/100 g DM | — (N x 6.25) | Kjeldahl | Table 1 | measured |
| SPI (Naturz Organics) | protein content | 89.22 +/- 0.33 | g/100 g DM | — (N x 6.25) | Kjeldahl | Table 1 | measured |
| PPI | moisture | 7.29 +/- 0.11 | g/100 g DM (sic) | — | oven AACC 44-15.02 | Table 1 | measured (flag 3) |
| SPI | moisture | 6.62 +/- 0.08 | g/100 g DM (sic) | — | oven | Table 1 | measured (flag 3) |
| PPI, SPI | ash | NOT MEASURED | — | — | — | — | — |
| PPI | fat | 8.51 +/- 0.07 | g/100 g DM | — | AACC 30-25.01 | Table 1 | measured; the lipid the isolate carries into a cook (Programme 7 ii) |
| SPI | fat | 1.72 +/- 0.17 | g/100 g DM | — | AACC 30-25.01 | Table 1 | measured |
| **PPI** | **total lysine (protein-bound; free = 0)** | **6.399 +/- 0.891** | g/100 g DM | 63.99 mg/g powder / 0.8122 = 78.79 mg/g protein; / 146.19 = **0.539 +/- 0.075 mmol/g protein** (0.438 mmol/g powder) | ion chromatography, ninhydrin (Chelab) | Table 2 | measured; candidate `amine` for `pea_isolate` (now 0.524 from flour) |
| **SPI** | **total lysine (protein-bound; free = 0)** | **5.343 +/- 0.743** | g/100 g DM | 53.43 / 0.8922 = 59.89 mg/g protein; / 146.19 = **0.410 +/- 0.057 mmol/g protein** (0.366 mmol/g powder) | as above | Table 2 | measured; candidate `amine` for `soy_isolate` (now 0.379 from flour) |
| PPI | total arginine | 8.422 +/- 1.172 | g/100 g DM | 84.22 / 0.8122 = 103.69 mg/g protein; / 174.20 = **0.595 +/- 0.083 mmol/g protein** (0.484 mmol/g powder) | as above | Table 2 | measured; guanidino pool (methylglyoxal / hydroimidazolone chemistry), not the epsilon-amine pool |
| SPI | total arginine | 7.327 +/- 1.019 | g/100 g DM | 73.27 / 0.8922 = 82.12 mg/g protein; / 174.20 = **0.471 +/- 0.066 mmol/g protein** (0.421 mmol/g powder) | as above | Table 2 | measured |
| PPI | methionine + cysteine (sum only) | 0.605 +/- 0.075 | g/100 g DM | 6.05 / 0.8122 = 7.45 mg/g protein; **0.050 mmol/g if all Met (/149.21) to 0.062 if all Cys (/121.16)**; the split is not printed | as above | Table 2 | measured as a sum; neither Met nor Cys alone (flag 2) |
| SPI | methionine + cysteine (sum only) | 0.942 +/- 0.118 | g/100 g DM | 9.42 / 0.8922 = 10.56 mg/g protein; **0.071 (all Met) to 0.087 (all Cys) mmol/g protein** | as above | Table 2 | measured as a sum (flag 2) |
| **PPI** | **free arginine** | **0.053 +/- 0.009** | g/100 g DM | 0.53 mg/g powder / 174.20 = 3.0 umol/g powder = **0.0030 mmol/g powder; 0.0037 mmol/g protein** | ion chromatography (free fraction) | Table 2 | measured; the only free amino acid above 0.02 g/100 g in either isolate |
| PPI | free glutamine | 0.013 +/- 0.008 | g/100 g DM | 0.13 / 146.15 = 0.9 umol/g powder (0.0011 mmol/g protein) | as above | Table 2 | measured (sd 60 % of value) |
| PPI | free leucine | 0.013 +/- 0.008 | g/100 g DM | 0.13 / 131.17 = 1.0 umol/g powder (0.0012 mmol/g protein) | as above | Table 2 | measured (sd 60 % of value) |
| PPI | free lysine, free methionine, free histidine, free Glu, Ala, Phe, Ile, Val, citrulline | 0 | g/100 g DM | 0 | as above | Table 2 | measured absent ("0"); Asn, Gly, Orn, Ser, Tyr, Thr, Asp "n.d." (flag 5) |
| **PPI** | **total free amino acids** | **0.079 +/- 0.014** | g/100 g DM | molar sum of the three = 4.9 umol/g powder = **0.0049 mmol/g powder; 0.0061 mmol/g protein**; = 1.1 % of the bound lysine (0.438) on a molar basis | as above | Table 2 | measured; at 100 g isolate/L this is 0.49 mM total free amino acid, 0.30 mM free Arg |
| **SPI** | **total free amino acids (every row)** | **0 +/- 0** (each row 0 or n.d.) | g/100 g DM | **0** | as above | Table 2 | measured absent; at any dose the soy isolate charges no free amino acid (upper bound not stated; flag 5) |
| **PPI** | **total free sugars (sucrose + glucose + fructose + maltose)** | **0.19 +/- 0.00** | g/100 g DM | 1.9 mg/g powder: **0.0055 mmol/g if all sucrose / maltose (/342.30) to 0.0105 if all glucose / fructose (/180.16)**; per g protein 0.0068-0.0130 | HPLC-RI, Sugar-Pak I, external standards | Table 1 | measured as a total; split not printed (flag 2) |
| **SPI** | **total free sugars** | **0.06 +/- 0.00** | g/100 g DM | 0.6 mg/g powder: **0.0018 (disaccharide) to 0.0033 (hexose) mmol/g powder**; per g protein 0.0020-0.0037 | as above | Table 1 | measured as a total (flag 2) |
| PPI | total FODMAPs (mono-, di-, GOS, fructans, polyols) | 0.53 +/- 0 | g/100 g DM | not convertible (mixed oligomers; raffinose 504, stachyose 667 g/mol would give 0.008-0.011 mmol/g if all GOS) | HPAEC-PAD, triplicate | Table 1 | measured as a total; exceeds "total sugars", so ~0.34 g/100 g is oligosaccharide (GOS / fructan), non-reducing until hydrolysed |
| SPI | total FODMAPs | 0.04 +/- 0 | g/100 g DM | — | HPAEC-PAD | Table 1 | measured as a total |
| PPI / SPI | digestible starch | 2.83 +/- 0.06 / 1.67 +/- 0.03 | g/100 g DM | 0.157 / 0.093 mmol anhydroglucose per g powder if fully hydrolysed (162.14 g/mol) | Megazyme K-RAPRS | Table 1 | measured; a latent glucose source only under amylase or acid hydrolysis, not a free sugar |
| PPI / SPI | pH of the powder in water | 6.360 +/- 0.030 / 6.770 +/- 0.010 | — (10 g in 100 mL, 5 % acetone) | — | pH meter | Table 3 | measured; the isolate's own pH before any buffer |
| PPI / SPI | protein solubility at pH 7 | 22.267 +/- 1.457 / 51.960 +/- 3.354 | % of protein | — | 1 % w/v, 4893 g, Kjeldahl | Table 3 | measured; how much of the amine pool is in solution vs particles at 22 C |
| PPI / SPI | zeta potential at pH 7 | -22.6 +/- 1.6 / -33.8 +/- 1.5 | mV | — | 0.1 % w/v | Table 3 | measured |
| PPI / SPI | colour L* a* b* | 84.7 / 2.18 / 22.3 ; 83.8 / 0.73 / 18.5 | — | — | Minolta, powder | Table 3 | measured; baseline colour of the unheated powders for any browning observable |
| EverPro (BSG barley-rice isolate) | one line | protein 83.22, fat 0.37, sugars 0.30, FODMAPs n.d., lysine 3.035 (0.250 mmol/g protein), free amino acids 1.248 g/100 g DM (free Lys 0.024), pH 7.90, solubility 101.7 %, L* 57.4 | g/100 g DM | 30.35 / 0.8322 / 146.19 = 0.250 | as above | Tables 1-3 | measured; a hydrolysed cereal protein, not a candidate matrix here |

Comparison the matrix layer can print (arithmetic, not a fit): the repo's flour-derived amine
pools are pea 0.524 and soy 0.379 mmol/g protein; this paper's isolates give 0.539 +/- 0.075 and
0.410 +/- 0.057 — 3 % and 8 % higher, inside the stated uncertainty, so the flour values were
not far wrong but the isolate values are the right basis. Per gram of POWDER (what a recipe
weighs) the pools are 0.438 (pea) and 0.366 (soy) mmol/g, i.e. 43.8 and 36.6 mM per 100 g
isolate/L, against a free-amino-acid charge of 0.49 mM (pea) and 0 (soy) and a free-sugar charge
of 0.55-1.05 mM (pea) and 0.18-0.33 mM (soy) at the same dose.

## 5. Flags

1. **The +/- on the amino-acid table is a fixed relative uncertainty, not replicate scatter.**
   Every total-amino-acid row of PPI and SPI carries +/- 13.9 % of its value (EverPro 10 %,
   tryptophan 10.8 %, the two summed pairs 10-12 %), which is what a contract laboratory
   (Chelab) attaches as measurement uncertainty. The paper's "All experiments were performed in
   triplicate" cannot be assumed to cover the external amino-acid analysis; the number of
   hydrolysates is not stated. Enter 0.539 / 0.410 mmol/g with the 14 % band as a stated
   uncertainty, not as a standard deviation of n = 3.
2. **Sums where the model wants parts.** Methionine + cysteine is one number (so no cysteine per
   g protein for the sulfur lane from this paper — only a 0.050-0.062 / 0.071-0.087 mmol/g
   protein band for the pair), phenylalanine + tyrosine is one number, and "Total Sugars" lumps
   sucrose, glucose, fructose and maltose although the HPLC-RI method resolves them. The
   reducing fraction of the 0.19 / 0.06 g/100 g is unknown; the band in §4 spans all-sucrose to
   all-hexose.
3. **Moisture printed under "g/100 g DM".** Moisture on a dry-matter basis is a contradiction;
   the 7.29 / 6.62 are almost certainly g per 100 g powder as received, and whether the other
   rows were actually corrected to dry matter is not verifiable (the protein values 81.22 / 89.22
   would be 75.3 / 83.3 on an as-is basis if the correction was applied, and the reverse if not).
   The per-g-protein ratios in §4 are unaffected because numerator and denominator share the
   basis; the per-g-powder numbers carry up to a 7 % basis ambiguity.
4. **No ash.** Ash (and hence the salt load of the isolate, which sets ionic strength and the
   Na/K/Ca/P the isolate brings) is not measured; Sägesser 2024 has it for other products but
   only as a figure.
5. **"0" versus "n.d." in the free-amino-acid table** are not defined; no LoD or LoQ for the free
   amino acids is printed (the 0.025 g/100 g LoQ belongs to the FODMAP method). "0" for SPI is
   therefore an upper bound of unknown size, plausibly of order 0.01 g/100 g DM (the smallest
   printed PPI value is 0.013). Free cysteine, free proline and free tryptophan are not rows at
   all.
6. **Products are anonymous.** "Naturz Organics (Helmond, Netherlands)" is a distributor; no
   product name, lot, cultivar, origin or process is given, and the isolates' extraction route
   is inferred ("likely" alkaline / isoelectric). The 8.51 % fat in the PPI is high for a pea
   isolate (Gao 2020's lab isolate: 1.5 %) and marks it as an unfat-extracted product; the
   hexanal / 2-pentylfuran carry-over of Programme 7 (ii) would scale with it.
7. **Total-amino-acid recovery differs between the two isolates** (sum 105 % of protein for PPI,
   92 % for SPI). With N x 6.25 for both, and pea's true factor nearer 5.4-5.7, the PPI protein
   content is overstated and its per-g-protein lysine correspondingly understated relative to
   soy; nothing to correct, but the two products are not on an identical footing.
8. **Table 2's "% Requirement" is per g protein, the amino-acid values are per g powder.** Easy to
   misread as one basis; the §3 check (175.1 %, 133.1 %, 33.9 %, 48.0 % all reproduce from Table
   1's protein content) settles it.
9. **The moisture-corrected free-sugar totals are at the method's floor** (0.06 +/- 0.00 for SPI
   with a printed sd of 0.00) — three significant figures should not be read into them.
10. **The MDPI PDF duplicates every page** (peer-review layer under the typeset one); a text search
    finds each value twice. Values here were taken from the typeset layer and cross-checked with
    pypdf layout extraction of pp. 7-8.
