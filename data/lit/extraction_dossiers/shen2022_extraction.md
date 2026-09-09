# Shen, Hong, Singh, Koppel & Li 2022 — EXTRACTION (lab pea protein isolate: free SH 13.5 µmol/g, free NH2 "8.44 mmol/g", H0 202 000; enzyme / gum modifications)
### One baseline row for unmodified pea protein isolate; method, protein basis and S-S all delegated to a paper not on disk.

**Source on disk:** `data/articles/shen2022.pdf` (owner's download, 2026-09-08). Read from the pdftotext
layer in the scratchpad (`articles/shen2022.txt`); Tables 1-3 have a clean text layer (cell values and
significance letters interleaved line by line) and Table 1 is re-typed in full below. Supplementary
document (sensory references) not on disk. Repo status before this dossier: no pea site density on
file.

## 0. Identity

| field | value |
|---|---|
| Title | "Improving functional properties of pea protein through 'green' modifications using enzymes and polysaccharides" |
| Authors | Yanting Shen, Shan Hong, Gaganpreet Singh, Kadri Koppel, Yonghui Li (Kansas State University, Grain Science and Industry / Food, Nutrition, Dietetics and Health) |
| Venue | Food Chemistry 385 (2022) 132687; received 6 Dec 2021, accepted 10 Mar 2022, online 12 Mar 2022 |
| DOI | 10.1016/j.foodchem.2022.132687 |
| Naming | "Control" = unmodified lab pea protein isolate; PG = protein-glutaminase deamidated; TG = transglutaminase cross-linked; Guar / Arabic = conjugated with guar gum / gum arabic (60 C, 24 h); PG-Guar etc. = sequential. "Free SH" and "Free NH2" as in Table 1. |
| Method parent | Shen & Li 2021, Food Hydrocolloids 117, 106686 (all SH, NH2, H0, FTIR, digestibility and functionality methods "without any modification"); NOT on disk |

## 1. Why it matters

The matrix layer needs a pea free-thiol density. This paper prints one for an unmodified laboratory pea
protein isolate, 13.52 ± 0.09 µmol/g, roughly 6x the pea-globulin value in Chihi 2016 (2.1 µmol/g) and
close to soy (7.5-8 µmol/g). It prints no S-S. It also prints a free amino-group content of 8.44 mmol/g
which, if taken at face value, is ~17x the lysine content of pea protein and is therefore flagged as
unusable without the parent paper's units. The modifications are of secondary interest: they show
that any heat step (100 C enzyme kill; 60 C / 24 h conjugation) lowers the titratable free SH of pea
protein by 12-64 %.

## 2. Methods as they matter to a model

- **Protein (§2.2):** laboratory isolate from ADM yellow pea flour: hexane-defatted; 10 % solids in
  water, pH 8.5 (1 M NaOH), 1 h room temperature, 8000 g 20 min 4 C; supernatant to pH 4.5 (1 M HCl),
  2 h 4 C; 8000 g; washed twice; re-adjusted to pH 7.0; lyophilised; 4 C. **Protein content of the
  isolate is NOT stated** in this paper (no Kjeldahl / Dumas figure, no nitrogen factor). Alkaline
  extraction + isoelectric precipitation (same route as Ruan's soy and Shen's own commercial-PPI
  comparisons), no albumin removal, so this is a whole pea protein isolate (globulins + albumins), not
  a globulin fraction.
- **Modifications (§2.3):** 10 % protein in water; PG 1 % (protein basis) pH 6.5, 55 C, 3 h; TG 1 %,
  40 C, 3 h; both then "heated to 100 C to inactivate the enzyme" (§3.4 says "boiling the protein
  slurries at 100 C for 10 min"); conjugation 5 % gum (protein basis), 60 C, 24 h; sequential = enzyme,
  kill, then gum. All lyophilised.
- **SH method, verbatim (the whole of it):** "Protein physicochemical properties, including free
  sulfhydryl group content, free amino group content, protein secondary structures, surface
  hydrophobicity, and in vitro gastrointestinal digestibility were determined following previous methods
  without any modification (Shen & Li, 2021)." Reagent, buffer, denaturant, whether "exposed" or
  "total", protein basis of "µmol/g": ALL UNSTATED here. No total-SH / no reducing agent / no S-S
  measurement anywhere in the paper.
- **Surface hydrophobicity:** method also delegated; H0 values of order 10^5 (arbitrary fluorescence
  slope units), not comparable in scale with Chihi's H0 x 10^6 of order 1-4.
- **Replicates:** "All the tests were conducted in at least duplicates, and the results were presented
  as mean ± standard deviation (SD)."
- **Units in the repo:** µmol/g -> mmol/g by / 1000. Basis (per g isolate powder vs per g protein) is
  ambiguous; see §4 for the two readings.

## 3. Tables re-typed

### Table 1. "Physicochemical properties including free sulfhydryl group content, free amino group content, secondary structures of pea and modified pea proteins."
Footnotes: "*Means with different letters in each column indicate significant differences (p < 0.05).
** ND: not detected."

| sample | free SH (µmol/g) | free NH2 (mmol/g) | α-helix (%) | β-sheet (%) | β-turn (%) | random coil (%) | hydrophobicity H0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Control** | **13.52 ± 0.09 a** | **8.44 ± 0.06 a** | 18.64 ± 0.09 cd | 27.52 ± 4.37 bc | 11.48 ± 2.27 bc | 42.37 ± 6.55 a | **202,096 ± 12,306 b** |
| PG | 9.91 ± 0.06 c | 7.53 ± 0.13 b | 53.72 ± 0.48 a | 26.67 ± 2.96 bc | 19.61 ± 2.48 a | ND | 161,826 ± 1,274 d |
| TG | 11.86 ± 0.08 b | 5.30 ± 0.06 c | 21.97 ± 1.60 c | 60.69 ± 3.30 a | 17.33 ± 1.70 ab | ND | 73,910 ± 1,500 f |
| Guar | 7.68 ± 0.02 d | 7.31 ± 0.16 b | 37.62 ± 1.56 b | 52.54 ± 0.78 a | 9.84 ± 0.77 cd | ND | 93,342 ± 1,099 e |
| Arabic | 6.46 ± 0.03 ef | 7.56 ± 0.22 b | 41.20 ± 10.39 ab | 48.08 ± 9.75 ab | 10.71 ± 0.64 cd | ND | 105,724 ± 1,995 e |
| PG-Guar | 5.61 ± 0.00 g | 7.34 ± 0.29 b | 9.66 ± 0.38 cd | 47.96 ± 1.33 ab | 6.88 ± 0.10 cd | 35.51 ± 1.05 a | 186,742 ± 3,243 c |
| PG-Arabic | 4.89 ± 0.04 h | 7.41 ± 0.22 b | 7.39 ± 1.02 d | 58.92 ± 1.72 a | 4.85 ± 0.28 d | 28.84 ± 2.46 a | 230,281 ± 1,223 a |
| TG-Guar | 6.68 ± 0.04 e | 5.39 ± 0.06 c | 10.29 ± 2.14 cd | 56.58 ± 12.13 a | 6.05 ± 0.80 cd | 27.08 ± 15.06 a | 28,158 ± 1,846 g |
| TG-Arabic | 6.28 ± 0.15 f | 5.19 ± 0.10 c | 20.90 ± 0.54 cd | 22.33 ± 4.19 c | 8.86 ± 2.34 cd | 47.91 ± 5.99 a | 22,563 ± 1,098 g |

Text restatement: control free SH "13.5 µmol/g"; control free amino "8.44 mmol/g". Authors attribute
every SH decrease to air oxidation to disulfide during mixing, more so at the higher conjugation
temperature and with two heat steps (sequential samples).

### Table 2. Functional properties (condensed; control row full, modifications one line each)

| sample | WHC (g/g) | OHC (g/g) | EC (%) | ES (%) | LGC (%) |
|---|---:|---:|---:|---:|---:|
| **Control** | 2.66 ± 0.06 f | 2.76 ± 0.05 c | 58.58 ± 2.21 c | 48.14 ± 1.77 d | 11 d |
| PG | 3.62 ± 0.04 d | 2.68 ± 0.08 c | 63.46 ± 4.95 bc | 51.91 ± 0.95 cd | 15 a |
| TG | 5.31 ± 0.08 b | 3.08 ± 0.03 b | 94.51 ± 0.33 a | 57.69 ± 1.39 b | 11 d |
| Guar | 3.62 ± 0.04 d | 2.62 ± 0.04 cd | 97.94 ± 0.34 a | 96.31 ± 0.95 a | 9 e |
| Arabic | 2.66 ± 0.01 f | 2.50 ± 0.06 d | 57.79 ± 4.05 c | 52.11 ± 2.81 c | 13 b |
| PG-Guar | 5.06 ± 0.02 c | 3.36 ± 0.05 a | 100.00 ± 0.00 a | 97.74 ± 0.08 a | 12 c |
| PG-Arabic | 3.27 ± 0.03 e | 2.75 ± 0.04 c | 67.57 ± 1.48 b | 56.71 ± 2.15 b | 15 a |
| TG-Guar | 5.62 ± 0.04 a | 2.98 ± 0.07 b | 100.00 ± 0.00 a | 100.00 ± 0.00 a | 9 e |
| TG-Arabic | 5.21 ± 0.06 b | 2.70 ± 0.02 c | 66.51 ± 4.65 b | 54.62 ± 1.97 bc | 9 e |

(Text: control WHC 2.8 g/g in the abstract vs 2.66 in Table 2; commercial PPI comparators from Shen &
Li 2021: OHC 1.03 g/g, LGC 18 %.)

### Table 3. Descriptive sensory (0-15 scale), control row only
Beany 6, starchy 6, grain 5, green 3, pulpy 0, powdery mouthfeel 5.5, umami 2, astringent 2.5, bitter
2.5, metallic 1.5. Modified samples within ±1 except pulpy 3-5 for TG samples and umami 0 for PG,
PG-Guar, TG-Guar.

### Fig. 3 (digestibility, DH %): FIGURE-ONLY; control highest (letter a), all modified lower.

## 4. Site densities the repository can use

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN | conditions | source | evidence |
|---|---|---:|---|---:|---|---|---|
| pea protein isolate, lab-made (alkaline extraction / pI precipitation, whole isolate) | free SH | 13.52 ± 0.09 | µmol/g | **0.0135 if "g" = g protein; 0.0135 / f_p if "g" = g powder** (f_p = protein fraction, unstated; for a typical lab pea isolate of 80-90 % protein this is 0.015-0.017) | native, lyophilised, method conditions unstated | Table 1 | measured, basis ambiguous |
| same | S-S | not measured | — | — | — | — | absent |
| same | free amino (NH2) | 8.44 ± 0.06 | mmol/g | **do not use** (see flag 2) | | Table 1 | measured, unit implausible |
| same | surface hydrophobicity H0 | 202,096 ± 12,306 | arbitrary (ANS slope) | n/a | | Table 1 | measured, scale not transferable |
| same, after PG 55 C 3 h + 100 C 10 min | free SH | 9.91 ± 0.06 | µmol/g | 0.0099 (same basis caveat) | enzyme + heat confounded | Table 1 | measured |
| same, after TG 40 C 3 h + 100 C 10 min | free SH | 11.86 ± 0.08 | µmol/g | 0.0119 | enzyme + heat confounded | Table 1 | measured |
| same, after 60 C / 24 h with 5 % guar / gum arabic | free SH | 7.68 ± 0.02 / 6.46 ± 0.03 | µmol/g | 0.0077 / 0.0065 | Maillard conjugation + heat confounded | Table 1 | measured |

Nothing here is a clean "after heating only" value: every heated sample also carries an enzyme or a
reducing sugar polymer. The closest to a pure heat effect is TG (mild enzyme, 40 C, then 100 C / 10
min): -12 % free SH.

## 5. Flags

1. **The SH method is not in this paper.** Reagent (presumably Ellman's DTNB, by the authors' other
   work, but unverified), buffer, denaturant, and whether the 13.5 µmol/g is surface or unfolded SH are
   all in Shen & Li 2021 (not on disk). Until that paper is on disk the value cannot be placed on the
   same footing as Chihi (1.5 M GdnCl) or Shimada (6 M urea + SDS).
2. **Free amino group 8.44 mmol/g is physically implausible as a per-gram-protein lysine count.** Pea
   protein carries ~7 g lysine / 100 g protein = ~0.48 mmol ε-NH2 per g protein (plus ~0.03 mmol
   N-termini); total residues per gram are only ~9 mmol. 8.44 mmol/g would mean nearly every residue
   carries a free amine. Most likely a unit or calibration artefact (e.g. µmol/g mislabeled, or a
   leucine-equivalent OPA/TNBS scale). The repo's amine site for pea must NOT be taken from this table;
   the BLG comparator in `protein_matrices.yml` is 0.817 mmol/g.
3. **Protein basis of "µmol/g" is unstated** and the isolate's protein content is unstated. §4 gives
   both readings; the per-powder reading raises the density by 10-25 %.
4. **No S-S, no total SH, no half-cystine.** The paper's disulfide statements are inferences from the
   SH decrease and SEC peak shifts, not measurements.
5. **Baseline is not a globulin fraction.** Whole isolate including albumins (pea albumin PA1/PA2 are
   cysteine-rich; lipoxygenase and Bowman-Birk-type inhibitors also contribute), which is one plausible
   reason the free SH (13.5) is 6x Chihi's albumin-depleted globulin value (2.1). The two numbers are
   for different matrices.
6. **Abstract vs table**: control WHC 2.8 (abstract) vs 2.66 (Table 2); EC 58 / ES 48 consistent.
7. **Heating effects are all confounded** with enzyme or polysaccharide; the "after any heat step" rows
   in §4 are bounded, not clean.
8. **H0 units** are instrument-specific slopes (10^5 scale) and cannot be compared with Chihi's
   (10^0 scale after x 10^6) or used quantitatively.
