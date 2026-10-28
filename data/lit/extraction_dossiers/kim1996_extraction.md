# Kim, Taub, Choi & Prakash 1996 — EXTRACTION (review-plus-data chapter on intrinsic chemical markers of sterility: M-1 = DDMP, M-2 = norfuraneol, M-3 = HMF; spiking of beef extract with 1 % glucose / fructose / ribose / glucose-6-P / ribose-5-P at 121 °C, 20-100 min; M-1 / M-3 vs pH 2-12 in beef extract + glucose, 121 °C / 15 min; ohmic-heated meatballs and microwave-heated ham as applications)

### NO KINETICS AND NO Ea FOR M-2 IN THIS CHAPTER; ITS VALUE IS MECHANISTIC AND SEMI-QUANTITATIVE: (i) ribose AND ribose-5-phosphate both give M-2 in beef extract at 121 °C (ribose ~10× the control); (ii) the ribose → furfural vs ribose → M-2 switch lies at **pH ≈ 4.5** (text only, no numbers); (iii) the hexose analogue M-1 rises monotonically from pH 2 to ~9 and **collapses above pH 10** (121 °C, 15 min) — the only stability-vs-pH signal; (iv) the **only apparent M-2 loss in the four marker papers**: with ribose-5-P, M-2 is flat at 20-60 min then **~30 % lower at 80 and 100 min** (121 °C) — a step, not a first-order decay, while the free-ribose arm keeps rising. No amine dependence is measured; no temperature dependence of M-2.

**Source on disk:** `data/articles/kim1996.pdf` (16 pp., ACS Symp. Ser. publisher scan with OCR text layer). Prose from
`pdftotext -layout`; Figs. 2-4 read from 120/250-dpi renders (`scratch/lit8/kim-07.png`, `kim_f3.png`, `kim-09.png`) and
tagged [D]. Figs. 2-3 are 3-D bar charts: read error **±0.1 absorbance units (AU) on Fig. 3, ±0.01 on Fig. 2**.

## 0. Identity

| field | value |
|---|---|
| Title | "Principles and Applications of Chemical Markers of Sterility in High-Temperature—Short-Time Processing of Particulate Foods" |
| Authors | Hie-Joon Kim, Irwin A. Taub (US Army Natick RD&E Center), Yang-Mun Choi (Korea University), Anuradha Prakash (Ohio State) |
| Venue | ACS Symposium Series **631** (Chemical Markers for Processed and Stored Foods; Lee & Kim eds.), ch. 6, pp. 54-69, 1996 (publication date 5 May 1996) |
| DOI as printed | **10.1021/bk-1996-0631.ch006** ✔ matches the brief. **Identity verified.** |

## 1. Why it matters

Brief item 4: the M-2 (NF) mechanism from ribose / ribose-5-P via 2,3-enolisation, its pH window, the amine requirement,
any kinetics, and any M-2 stability vs pH or T (decay data bear on B61's declared hydroxide-decay barrier).

## 2. Methods as they matter to a model

- **Matrix.** Beef extract (preparation, dilution and pH **NOT described** in this chapter; heated meats are "pH 5.4",
  p. 61). Spikes of 1 % sugar (Figs. 2-3); 1 % (Fig. 4 caption) or 10 % (text, p. 61) glucose for the pH series, pH set with
  6 N HCl / NaOH.
- **Heating.** 121 °C (pressurised steam or sealed; vessel not described), 20-100 min (Figs. 2-3), 15 min (Fig. 4).
- **Analysis.** Anion-exclusion HPLC with photodiode array; M-2 at 5.8 min, λmax 285 nm; M-1 at 4.2 min, 298 nm (p. 57).
  **Results as absorbance at 285 nm (peak height/area unstated); no calibration to concentration.** Heated samples or
  extracts "could be frozen and stored for several days without affecting the analysis".
- **Replicates / error bars:** none shown.

## 3. Numbers

### Text facts

- p. 55: thiamin destruction D = 244 min at 122 °C → k = 0.0094 min⁻¹ (quoted, ref. 8) — the motivating "slow
  destruction" example.
- p. 57 (Fig. 1, beef, steam 121 °C, 1/6/15 min): M-2 "is formed quite rapidly and approaches a limiting value after
  15 min heating"; M-1 is slower and still rising at 15 min.
- p. 57: M-1 forms in meats and vegetables; **M-2 in meats only** (ribose / ribose-5-P are meat constituents).
- p. 59: "both ribose and ribose-5-phosphate participate in the formation of M-2 in meats"; the authors propose the
  phosphate of ribose-5-P is hydrolysed and free ribose reacts, whereas glucose-6-P is not hydrolysed and gives no M-1.
- p. 61: 1,2-enolisation (protonated Amadori N) in strong acid → furfural / HMF; "under weakly acidic or alkaline
  conditions 2,3-enolization is favored leading to furanones and pyranones"; Hicks & Feather 1975: Amadori → furfural in
  2 N H2SO4, M-2 at pH 7. **"When 1 % D-ribose was added to the beef extract, formation of 2-furaldehyde dominated at low pH
  and M-2 formation became more important at pH above 4.5"** — no data shown for the ribose series.
- Amine requirement: asserted through the Amadori pathway (Fig. 5); Anderegg & Neukom 1974 (cited) report M-2 from
  pentoses and nucleotides **without amines**. Not tested here.

### Fig. 3 (p. 60), M-2 (A285) in beef extract + 1 % sugar, 121 °C [D, ±0.1 AU]

| arm | 20 min | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|
| control | 0.13 | 0.17 | 0.18 | 0.18 | 0.22 |
| + glucose | 0.10 | 0.17 | 0.17 | 0.17 | 0.17 |
| + fructose | 0.13 | 0.20 | 0.17 | 0.17 | 0.17 |
| + ribose | 1.2 | 2.05 | 2.25 | 2.45 | 3.05 |
| + ribose-5-P | 1.25 | 1.25 | 1.25 | **0.88** | **0.88** |

### Fig. 2 (p. 60), M-1 (A285) in beef extract + 1 % sugar, 121 °C [D, ±0.01 AU]

| arm | 20 | 40 | 60 | 80 | 100 min |
|---|---|---|---|---|---|
| control | 0.045 | 0.080 | 0.085 | 0.090 | 0.110 |
| + ribose | 0.045 | 0.065 | 0.075 | 0.090 | 0.110 |
| + fructose | 0.065 | 0.105 | 0.135 | 0.145 | 0.185 |
| + glucose | 0.095 | 0.185 | 0.235 | 0.29 | 0.33 |
| + glucose-6-P | 0.025 | 0.055 | 0.060 | 0.080 | 0.105 |

### Fig. 4 (p. 62), M-1 and M-3 (A285) vs pH, beef extract + glucose, 121 °C / 15 min [D, ±0.003 AU, pH ±0.2]

| pH | ~2.1 | ~3 | ~4 | ~5 | ~6 | ~7 | ~8 | ~9 | ~10.5 | ~12 |
|---|---|---|---|---|---|---|---|---|---|---|
| M-1 (DDMP) | 0.006 | 0.009 | 0.019 | 0.034 | 0.047 | 0.052 | 0.080 | **0.095** | **0.007** | — |
| M-3 (HMF) | 0.075 | 0.075 | 0.018 | 0.005 | ~0 | ~0 | ~0 | ~0 | ~0 | ~0 |

### Arithmetic (mine)

1. **Ribose-5-P arm, apparent loss:** 1.25 → 0.88 AU between 60 and 80 min (−30 %), then flat to 100 min. A first-order
   reading over 60→100 min gives k ≈ ln(1.25/0.88)/40 = **8.8e-3 min⁻¹**, but the shape (flat, step, flat) is not
   first-order; batch or sampling artefact cannot be excluded (no replicates shown). In the same figure the free-ribose
   arm rises 2.25 → 3.05 over the same interval, so whatever removed M-2 in the ribose-5-P arm was not a matrix-wide M-2
   decay of that size.
2. **B61 comparison:** B61's hydroxide law at 121 °C and pH 5.4-5.8 (meat) gives 1.0-1.7e-3 min⁻¹ (with the neutral term),
   i.e. 2-3 % loss in 20 min — **an order of magnitude too small to produce the ribose-5-P step, and the ribose arm shows
   no loss at all**, so Kim neither supports nor contradicts B61.
3. **M-1 collapse above pH 10** (0.095 → 0.007, −93 % in 15 min at 121 °C): M-1 is the hexose pyranone enol analogue of NF.
   A base-driven loss of the enol at pH ≥ 10 is in the direction B61 assumes for NF, but (a) the collapse may be a
   formation switch (the authors cite Shaw 1968: M-1 is not formed at pH 11.5 but forms as pH falls), and (b) M-1 rises
   1.8× from pH 6 to pH 9, which a hydroxide decay of 10^(pH−5) with B61's NF magnitude would hardly allow for an enol of
   similar reactivity (B61 NF at 121 °C, pH 9: ~1.9 min⁻¹, t½ ≈ 22 s). **DDMP ≠ NF; no transfer is licensed.**

## 4. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| ribose → furfural vs → M-2 crossover | pH | ~4.5 | — | beef extract + 1 % ribose, 121 °C | text p. 61 | qualitative (no data shown) |
| glucose → HMF (M-3) vs → DDMP (M-1) crossover | pH | ~4 | — | beef extract + glucose, 121 °C, 15 min | Fig. 4 [D] | semi-quantitative |
| DDMP yield collapse | pH | > 10 | — | same | Fig. 4 [D] | semi-quantitative |
| M-2 from ribose vs control | ratio at 60 min | ~12 | — | beef extract + 1 % ribose, 121 °C | Fig. 3 [D] | within-study ratio |
| M-2 apparent loss (ribose-5-P arm) | step | −30 % (60→80 min) | — | 121 °C | Fig. 3 [D] | not first-order; unreplicated |
| M-2 kinetics / Ea | — | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** the qualitative pH switch between 1,2-enolisation (furfural/HMF, pH ≲ 4-4.5) and 2,3-enolisation (NF/DDMP,
pH ≳ 4.5-5) in a real meat matrix at 121 °C; ribose and ribose-5-P as NF precursors; hexoses giving no extra M-2.
**Cannot:** no rate constants, no Ea, no concentrations, no matrix pH, no amine levels, no replicates; no NF decay series
that can be fitted; nothing on T dependence of M-2 (Lau 2003's citation of "Kim & Taub 1993" for M-1 Ea 23.7 kcal/mol is
a different paper).

## 6. Flags

1. Fig. 4 caption says **1 %** glucose; text (p. 61) says **10 %** glucose. Unresolved.
2. Absorbance units throughout; whether peak height or area, and the injection volume, are not stated.
3. The 2,3- vs 1,2-enolisation crossover for ribose (pH 4.5) is asserted in the text without a figure.
4. Figures 2-3 are 3-D bar charts; ±0.1 AU reading error on Fig. 3 is material for the control-level arms.

## 7. Use in the model

- **(a) Furfural barrier swap:** nothing on temperature. Indirectly, its pH ~4.5 crossover (ribose in beef, 121 °C) is a
  qualitative target for the engine's 1,2- vs 2,3-enolisation pH law (SULFUR_FLUX_MAP sec. 5, pH 3 arm: the engine sends
  100 % of pentose to TDP at pH 3 and 0.4 % to DPO; a crossover near 4.5 is consistent with that direction, and does not
  by itself show the law is too steep).
- **(b) B61:** the only apparent M-2 loss (ribose-5-P arm) is step-shaped, unreplicated and absent from the ribose arm in
  the same figure; it **does not bear on the 103.7 kJ/mol barrier** and is ~10× larger than B61 predicts at meat pH. The
  DDMP collapse above pH 10 is directionally consistent with base-promoted enol destruction but concerns a different
  compound and may be a formation effect.
- **Mechanism record:** supports routing ribose-5-P → ribose → NF (the engine's PENT entry) for meat systems, and NF
  formation from ribose being a pH ≥ ~5 phenomenon. No amine-dependence data.
