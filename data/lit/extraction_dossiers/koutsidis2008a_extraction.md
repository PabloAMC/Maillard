# Koutsidis et al. 2008a — EXTRACTION (beef, 10 days conditioned: effect of diet and breed on sugars, nucleotides, free amino acids)

**Source on disk:** `data/articles/koutsidis2008a.pdf` (the publisher's PDF). Tables 1 to 3 read from the
PDF and checked by eye on 2026-09-14; this dossier was written from the PDF directly, there was no
earlier automated version. Written for the cultivated-tissue composition box
(`src/cultivated_tissue_invariance.py`), not for the core fit.

| field | value |
|---|---|
| Title | "Water-soluble precursors of beef flavour: I. Effect of diet and breed" |
| Authors | G. Koutsidis, J.S. Elmore, M.J. Oruna-Concha, M.M. Campo, J.D. Wood, D.S. Mottram |
| Venue | Meat Science 2008, 79(1), 124-130 |
| DOI | 10.1016/j.meatsci.2007.08.008 |
| Companion | Part II, effect of post-mortem conditioning: `koutsidis2008b_extraction.md` |

## 1. Methods (section 2 of the paper)

32 steers, 16 Aberdeen Angus × Holstein-Friesian (AA, a beef cross) and 16 Holstein-Friesian (HF, a dairy
breed), allocated at ~9 months to grass silage (with ~15 % sugar-beet pulp) or a restricted concentrate
(barley, molassed sugar-beet pulp, full-fat soya; 70:30 with barley straw), 8 per breed per diet; two HF
steers died, so n = 8 / 8 / 7 / 7. Slaughtered at 24 months; M. longissimus lumborum vacuum packed at 4 °C
for **10 days**, then steaks frozen at −20 °C. Cold-water extraction, 3000 Da ultrafiltration, rhamnose /
norvaline / purine internal standards. Sugars as trimethylsilyl ethers by GC-MS (Leblanc & Ball 1978);
nucleotides, creatine, creatinine, carnosine by capillary electrophoresis with diode-array detection;
free amino acids (except arginine) by EZ-Faast derivatisation and GC-MS. Tables 1 to 3 print group means
± standard errors in **mmol per kg** (wet basis) with a "concentration range" column, which is the
spread across individual animals over all four groups.

## 2. Findings that matter (mmol/kg meat; group mean ± SE, AA-concentrate / AA-silage / HF-concentrate / HF-silage; then the range over animals)

Table 1, sugars and sugar phosphates:

| compound | AA conc. | AA silage | HF conc. | HF silage | range over animals |
|---|---|---|---|---|---|
| glucose | 8.48 ± 0.18 | 8.24 ± 0.27 | 9.51 ± 0.38 | 8.12 ± 0.30 | 6.94–10.6 |
| fructose | 1.68 ± 0.13 | 1.64 ± 0.15 | 1.84 ± 0.24 | 1.57 ± 0.21 | 0.69–2.53 |
| mannose | 1.38 ± 0.06 | 1.31 ± 0.09 | 1.56 ± 0.10 | 1.27 ± 0.16 | 0.64–1.89 |
| **ribose** | 0.71 ± 0.03 | 0.67 ± 0.02 | 0.82 ± 0.06 | 0.78 ± 0.06 | **0.57–1.08** |
| glucose 6-phosphate | 6.21 ± 0.25 | 6.31 ± 0.29 | 7.22 ± 0.40 | 5.90 ± 0.48 | 3.71–8.56 |
| fructose 6-phosphate | 1.60 ± 0.18 | 1.47 ± 0.16 | 2.02 ± 0.26 | 1.34 ± 0.23 | 0.46–3.13 |
| mannose 6-phosphate | 1.82 ± 0.15 | 1.66 ± 0.17 | 2.39 ± 0.26 | 1.58 ± 0.28 | 0.46–3.49 |
| maltose | 0.08 ± 0.01 | 0.07 ± 0.01 | 0.11 ± 0.03 | 0.06 ± 0.01 | 0.02–0.23 |
| total reducing sugars | 22.0 ± 0.75 | 21.4 ± 0.97 | 25.5 ± 1.54 | 20.6 ± 1.62 | 13.5–31.1 |

Ribose 5-phosphate is not in this table (it is in Part II's). Ribose was higher in HF than AA (p < 0.05);
glucose higher on concentrate (p < 0.05); the diet effect on sugars was otherwise small.

Table 2, nucleotides and amino compounds:

| compound | AA conc. | AA silage | HF conc. | HF silage | range |
|---|---|---|---|---|---|
| IMP | 3.54 ± 0.14 | 3.46 ± 0.23 | 3.45 ± 0.22 | 3.51 ± 0.12 | 2.51–4.44 |
| inosine | 1.37 ± 0.08 | 1.77 ± 0.09 | 1.73 ± 0.10 | 1.73 ± 0.10 | 1.01–2.23 |
| hypoxanthine | 1.58 ± 0.03 | 1.73 ± 0.10 | 1.76 ± 0.08 | 1.89 ± 0.04 | 1.41–2.33 |
| creatine | 47.2 ± 1.10 | 45.1 ± 1.90 | 44.9 ± 0.76 | 46.2 ± 1.11 | 40.3–55.8 |
| creatinine | 0.70 ± 0.02 | 0.72 ± 0.04 | 0.63 ± 0.04 | 0.65 ± 0.03 | 0.52–0.92 |
| carnosine | 30.1 ± 1.15 | 26.8 ± 1.64 | 29.3 ± 1.36 | 29.0 ± 0.70 | 22.8–34.3 |

Table 3, free amino acids (the rows the box or the sulfur lane cares about; 23 rows in the paper):

| compound | AA conc. | AA silage | HF conc. | HF silage | range |
|---|---|---|---|---|---|
| **cysteine** | 0.09 ± 0.01 | 0.12 ± 0.01 | 0.10 ± 0.01 | 0.10 ± 0.01 | **0.05–0.17** |
| methionine | 0.23 ± 0.01 | 0.30 ± 0.01 | 0.23 ± 0.02 | 0.29 ± 0.02 | 0.15–0.40 |
| leucine | 1.22 ± 0.08 | 1.63 ± 0.09 | 1.33 ± 0.08 | 1.55 ± 0.17 | 0.78–2.40 |
| glycine | 0.57 ± 0.02 | 0.71 ± 0.02 | 0.63 ± 0.03 | 0.77 ± 0.03 | 0.46–0.91 |
| alanine | 3.35 ± 0.21 | 3.44 ± 0.15 | 3.41 ± 0.11 | 3.19 ± 0.09 | 2.46–4.36 |
| glutamic acid | 1.04 ± 0.06 | 1.15 ± 0.09 | 0.86 ± 0.10 | 1.08 ± 0.09 | 0.58–1.58 |
| lysine | 0.85 ± 0.05 | 1.08 ± 0.07 | 0.91 ± 0.07 | 1.04 ± 0.09 | 0.60–1.38 |
| total free amino acids | 18.5 ± 0.90 | 20.3 ± 0.68 | 19.1 ± 1.17 | 18.6 ± 0.63 | 14.3–23.2 |

Silage-fed animals had higher free amino acids (most rows p < 0.01); cysteine showed no diet or breed
effect. Nothing on thiamine. Nothing on cultured cells.

## 3. What the repo takes

Beef-side values for the composition box, in mM in tissue water at 75 % moisture (mmol/kg ÷ 0.75),
taken over the paper's range across individual animals (30 steers, two breeds, two diets, one
conditioning time of 10 days):

- **ribose (free)**: 0.57 to 1.08 mmol/kg → **0.76 to 1.44 mM**. Inside Part II's 21-day span (0.33 to
  2.2 mM), and consistent with it: Part II's day-7 and day-14 means (0.76 and 1.26 mmol/kg) bracket
  these 10-day values.
- **glucose (free)**: 6.94 to 10.6 mmol/kg → **9.3 to 14.1 mM**. Inside the box's glucose range.
- **cysteine (free)**: 0.05 to 0.17 mmol/kg → **0.067 to 0.23 mM**. This range across animals sets the
  box's beef cysteine upper corner (Muroya 2019's D14 mean gives 0.14 mM, Part II's day-21 mean 0.21).
- **leucine (free)**: 0.78 to 2.40 mmol/kg → 1.04 to 3.2 mM: sets the box's beef leucine upper corner.
  **IMP**: 2.51 to 4.44 mmol/kg → 3.3 to 5.9 mM, inside the box.

The value of this paper for the box is the animal-to-animal spread at one conditioning time, which the
other beef sources (three steers over time; sixteen steers over time; fourteen bulls over time) do not
print as a range.
