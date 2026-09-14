# Koutsidis et al. 2008b — EXTRACTION (beef conditioning, 21 days: sugars, nucleotides, free amino acids)

**Source on disk:** `data/articles/koutsidis2008b.pdf` (the publisher's PDF). Tables 1 to 4 read from the
PDF and checked by eye on 2026-09-14. This dossier replaces the abstract-only version written earlier the
same day from the AGRIS record: that version carried only the "sixfold" ribose statement and no absolute
value; every number below is from the paper's own tables. Written for the cultivated-tissue composition
box (`src/cultivated_tissue_invariance.py`), not for the core fit.

| field | value |
|---|---|
| Title | "Water-soluble precursors of beef flavour. Part II: Effect of post-mortem conditioning" |
| Authors | G. Koutsidis, J.S. Elmore, M.J. Oruna-Concha, M.M. Campo, J.D. Wood, D.S. Mottram |
| Venue | Meat Science 2008, 79(2), 270-277 |
| DOI | 10.1016/j.meatsci.2007.09.010 |
| Companion | Part I, effect of diet and breed: `koutsidis2008a_extraction.md` |

## 1. Methods (section 2 of the paper)

Sixteen Charolais steers on a concentrate diet (IGER Aberystwyth), slaughtered at Bristol, electrically
stimulated (90 V DC, 1 min), chilled below 7 °C. M. longissimus lumborum vacuum packed at 4 °C; steaks
cut after 1, 3, 7, 14 and 21 days (glycogen after 2, 4, 7, 14, 21 days), blast frozen, stored at −18 °C.
Cold-water extraction, 3000 Da ultrafiltration. Sugars and sugar phosphates as trimethylsilyl ethers by
GC-MS (Leblanc & Ball 1978, modified); nucleotides, carnosine, creatine, creatinine and arginine by
capillary electrophoresis with diode-array detection; free amino acids (except arginine) by EZ-Faast
derivatisation and GC-MS; free phosphate by the ascorbic acid method; glycogen after perchloric acid
extraction. Tables 2 to 4 print means ± standard errors of 16 replicates, in **mmol per kg of meat**
(wet basis); different letters in a row mark p < 0.05.

## 2. Findings that matter (mmol/kg meat, mean ± SE, n = 16; conditioning day 1 / 3 / 7 / 14 / 21)

Table 2, sugars and sugar phosphates:

| compound | d1 | d3 | d7 | d14 | d21 |
|---|---|---|---|---|---|
| glucose | 7.33 ± 0.20 | 7.73 ± 0.25 | 7.80 ± 0.26 | 9.18 ± 0.36 | 10.3 ± 0.36 |
| fructose | 1.81 ± 0.17 | 1.98 ± 0.13 | 2.66 ± 0.18 | 3.18 ± 0.21 | 3.81 ± 0.18 |
| mannose | 1.21 ± 0.06 | 1.46 ± 0.06 | 2.00 ± 0.12 | 2.37 ± 0.12 | 2.93 ± 0.14 |
| **ribose** | **0.25 ± 0.02** | 0.43 ± 0.02 | 0.76 ± 0.04 | 1.26 ± 0.06 | **1.67 ± 0.09** |
| ribose 5-phosphate | 0.04 ± 0.003 | 0.04 ± 0.003 | 0.04 ± 0.004 | 0.04 ± 0.003 | 0.04 ± 0.003 |
| ribulose 5-phosphate | 0.07 ± 0.004 | 0.08 ± 0.006 | 0.07 ± 0.005 | 0.07 ± 0.005 | 0.06 ± 0.006 |
| xylulose 5-phosphate | 0.17 ± 0.01 | 0.17 ± 0.01 | 0.15 ± 0.01 | 0.15 ± 0.01 | 0.13 ± 0.01 |
| glucose 6-phosphate | 8.79 ± 0.51 | 8.83 ± 0.62 | 6.74 ± 0.44 | 7.05 ± 0.53 | 6.06 ± 0.41 |
| fructose 6-phosphate | 2.06 ± 0.12 | 2.06 ± 0.15 | 1.57 ± 0.10 | 1.73 ± 0.14 | 1.46 ± 0.09 |
| mannose 6-phosphate | 3.02 ± 0.16 | 3.08 ± 0.20 | 2.34 ± 0.14 | 2.44 ± 0.17 | 2.05 ± 0.12 |
| free phosphate | 26.0 ± 0.90 | 29.1 ± 0.89 | 32.1 ± 0.85 | 33.1 ± 1.09 | 35.3 ± 0.69 |
| total reducing sugars | 24.8 ± 0.97 | 25.9 ± 1.19 | 24.1 ± 0.98 | 27.5 ± 1.25 | 28.5 ± 1.12 |

Ribose rose 6.7-fold over the 21 days (the abstract's "sixfold"); the text says the rise was linear
(0.25 mmol/kg at 24 h to 1.67 at 21 d) and pairs it with the IMP fall (Fig. 2). Ribose 5-phosphate did
not change. Glycogen (Table 1, g/kg) was 1.85 to 1.98 at every time from 2 to 21 days: no further
glycogen breakdown after 48 h.

Table 3, nucleotide degradation products and amino compounds:

| compound | d1 | d3 | d7 | d14 | d21 |
|---|---|---|---|---|---|
| IMP | 6.27 ± 0.17 | 5.79 ± 0.14 | 4.61 ± 0.10 | 3.70 ± 0.11 | 2.63 ± 0.12 |
| inosine | 1.14 ± 0.04 | 1.36 ± 0.06 | 1.47 ± 0.07 | 1.87 ± 0.09 | 1.91 ± 0.09 |
| hypoxanthine | 0.84 ± 0.02 | 1.13 ± 0.03 | 1.52 ± 0.03 | 2.09 ± 0.04 | 2.66 ± 0.05 |
| GMP | 0.11 ± 0.005 | 0.10 ± 0.006 | 0.09 ± 0.005 | 0.08 ± 0.004 | 0.06 ± 0.004 |
| creatine | 40.5 ± 0.70 | 42.3 ± 1.78 | 43.1 ± 1.14 | 36.9 ± 1.32 | 38.8 ± 0.99 |
| creatinine | 0.65 ± 0.02 | 0.72 ± 0.03 | 0.87 ± 0.02 | 0.86 ± 0.03 | 0.99 ± 0.03 |
| carnosine | 33.8 ± 0.93 | 32.7 ± 0.74 | 31.1 ± 0.89 | 30.5 ± 0.55 | 28.5 ± 0.44 |

Table 4, free amino acids (the rows the box or the sulfur lane cares about; the table has 23 rows):

| compound | d1 | d3 | d7 | d14 | d21 |
|---|---|---|---|---|---|
| **cysteine** | **0.05 ± 0.005** | 0.06 ± 0.008 | 0.07 ± 0.009 | 0.12 ± 0.02 | **0.16 ± 0.02** |
| methionine | 0.05 ± 0.005 | 0.07 ± 0.005 | 0.11 ± 0.007 | 0.24 ± 0.01 | 0.35 ± 0.02 |
| leucine | 0.43 ± 0.03 | 0.46 ± 0.03 | 0.59 ± 0.03 | 1.29 ± 0.07 | 1.75 ± 0.11 |
| glycine | 0.53 ± 0.02 | 0.53 ± 0.02 | 0.58 ± 0.02 | 0.66 ± 0.01 | 0.73 ± 0.02 |
| alanine | 1.82 ± 0.14 | 1.76 ± 0.12 | 1.88 ± 0.07 | 2.40 ± 0.08 | 2.52 ± 0.10 |
| glutamic acid | 0.43 ± 0.04 | 0.39 ± 0.02 | 0.49 ± 0.02 | 0.79 ± 0.04 | 0.97 ± 0.05 |
| lysine | 0.31 ± 0.02 | 0.30 ± 0.02 | 0.36 ± 0.02 | 0.62 ± 0.03 | 0.79 ± 0.05 |
| total free amino acids | 10.8 ± 0.59 | 11.0 ± 0.52 | 11.9 ± 0.36 | 16.2 ± 0.52 | 18.9 ± 0.58 |

Cysteine rose threefold, methionine sevenfold; most of the amino-acid rise fell between days 7 and 14.
Nothing on thiamine. Nothing on cultured cells: the paper is about slaughtered beef only.

## 3. What the repo takes

Beef-side ranges for the composition box, converted to mM in tissue water at 75 % moisture
(mmol/kg ÷ 0.75), spanning the day-1 and day-21 means:

- **ribose (free)**: 0.25 to 1.67 mmol/kg → **0.33 to 2.2 mM**. The first PRIMARY beef ribose range on
  file; it replaces the earlier construction (a cited 0.26 mg/g point from Aliani 2013 via Hwang 2026 as
  the upper corner and one sixth of it as the lower). That cited point, 2.3 mM, sits at the top edge of
  the primary range. The companion paper's 10-day values (`koutsidis2008a`, 0.76 to 1.44 mM across 30
  animals) sit inside it.
- **glucose (free)**: 7.33 to 10.3 mmol/kg → **9.8 to 13.7 mM**. Above Bischof 2023's NMR range
  (`bischof2023`, corrected 2026-09-14: 5.9 to 13.3 mM over means) at the low end and consistent with it
  at the top.
- **cysteine (free)**: 0.05 to 0.16 mmol/kg → **0.067 to 0.21 mM**. Muroya 2019's CE-TOFMS span
  (0.002 to 0.14 mM) starts lower because it includes a pre-rigor D0 point; at matched times the two
  methods agree within twofold (Koutsidis d14 0.12 vs Muroya D14 0.107 mmol/kg).
- **leucine (free)**: 0.43 to 1.75 mmol/kg → 0.57 to 2.3 mM. **IMP**: 2.63 to 6.27 mmol/kg → 3.5 to
  8.4 mM. **ribose 5-phosphate**: 0.04 mmol/kg → 0.053 mM. All three inside the ranges already in the box
  from Muroya 2019; the box cites this paper on them as a second laboratory.

Sixteen animals, one breed, one diet, one laboratory. Every range above is a span over conditioning
time; it brackets what a beef reference could be, it does not estimate a population.
