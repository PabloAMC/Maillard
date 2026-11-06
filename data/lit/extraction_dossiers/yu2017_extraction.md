# Yu, Xu & Yu 2017 (J. Sci. Food Agric. 97:3159-3165) — EXTRACTION (file stem `yu2017`)
**DOI:** 10.1002/jsfa.8158 (printed in the Accepted-Article banner and the PDF keywords; Crossref 2026-11-05: authors Pei Yu, Xian-Bing Xu, Shu-Juan Yu, J. Sci. Food Agric. 97(10) 3159-3165, online 2017-01-27, issue 2017-08) — matches the brief.
### *"The effect of pH and amino acids on the formation of methylglyoxal in glucose-amino acid model system"*
### Glucose + Lys, Arg, Gly or Pro, **40 mM each (1:1)**, 4 mL in screw-capped Pyrex tubes, **oil bath 120 °C**; (i) 0.2 M phosphate pH 7.4, 0.5-4 h; (ii) 0.2 M Na2HPO4-citric acid **pH 5.0 / 6.0 / 7.0 / 7.4 / 8.0**, 0.5 h ("nearly no significant pH change" after heating). MGO by OPD / HPLC-DAD (LOQ range 0.002-3.2 mM, RSD ±9 %, recovery 96 %); glucose and amino acids by HPAEC-ED (peak-height loss); acetic acid by HPAEC-CD. **Table 2: MGO at 0.5 h rises 0.50-0.76 dex per pH unit from pH 5 to 7 and 0.18-0.36 from 7 to 8; at fixed pH, Lys gives 2-7× the MGO of Gly and Pro 0.2-0.8×.** No thiol, no MFT, one temperature.

**Year and key.** The file on disk is the **Accepted Article** (21 pp., PDF created 2016-11-25; no volume, issue or year
printed in it). The queue listed it as 2016 from that date. Crossref: online 2017-01-27, issue 97(10), August 2017. Both
publication dates are 2017, so the key is `yu2017` (the file stem). `yu2016` in the repo is a different paper (Yu & Tang
2016, ascorbic acid + cysteine).

**Read 2026-11-05 (first reader) from `data/articles/yu2017.pdf`.** Prose and Table 2 from `pdftotext -raw`. **Figs. 1,
3a, 3b and 4 read by eye** from 90 dpi renders of pp. 18-21 (scratchpad `m4img/yu2017_90-18.png` … `-21.png`): **read
error ±0.0007 mol/mol (Fig. 1), ±2 % (Fig. 3), ±0.02 mmol/mmol (Fig. 4).** Fig. 2 (A294, A420) not digitised. Arithmetic:
scratchpad `m4tools/m4_arith.py`. Page numbers below are those of the Accepted Article (p. 1-21).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Pei Yu, Xian-Bing Xu, **Shu-Juan Yu** (corresponding, shujuanyu8@gmail.com) |
| Affiliation | College of Light Industry and Food Sciences, South China University of Technology, Guangzhou; National Engineering Research Center of Seafood, Dalian Polytechnic University (Xu); Guangdong Province Key Laboratory for Green Processing of Natural Products and Product Safety |
| Venue | *J. Sci. Food Agric.* **97**(10) (2017) 3159-3165 (Crossref; the Accepted Article prints none of these) |
| DOI | 10.1002/jsfa.8158 |

**Laboratory.** Shu-Juan Yu's group, SCUT Guangzhou. Not in `data/keys/laboratories.yml`; proposed id `scut_yu_sj`.

## 1. What this paper answers for the model

| model need | what Yu 2017 gives | tag |
|---|---|---|
| MGO from glucose + amino acid vs pH | Table 2: four amino acids × five pH, 0.5 h, 120 °C, mean ± SD | [M] |
| MGO vs time | Fig. 1: four amino acids, pH 7.4 phosphate, 0.5-4 h, mol / mol glucose | [F] |
| glucose loss, amino-acid loss | Fig. 3a / 3b (% of t = 0, by HPAEC-ED peak height) | [F] |
| acetic acid | Fig. 4 (mmol / mmol glucose) | [F] |
| amine identity at fixed sugar | Lys > Gly > Arg ≈ Pro for MGO; Arg consumed fastest; Arg gives most acetic acid | [M] / [F] |
| thiol, MFT, NF, temperature | **none** (120 °C only) | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Pots | 0.8 mmol amino acid + 0.8 mmol glucose in 20 mL buffer [Z: 40 mM each]; 4 mL per screw-capped Pyrex tube |
| Time series | 0.2 M phosphate pH 7.4; oil bath 120 °C; 0.5, 1, 2, 3, 4 h; ice-water cooling |
| pH series | 0.2 M Na2HPO4-citric acid buffer, pH 5.0, 6.0, 7.0, 7.4, 8.0; 120 °C, 0.5 h |
| pH scale | buffer as prepared (label); after heating "Nearly no significant pH change was observed" (p > 0.05) [V] — so **held** within an unstated tolerance |
| Headspace | 4 mL in a tube of unstated volume; air presumed; no degassing |
| Replicates | "All experiments were carried out in triplicates" [V]; mean ± SD |
| MGO derivatisation | 1 mL sample + 1 mL OPD (0.5 mg/mL in 0.2 M PBS pH 7.4), 12 h at room temperature |

## 3. Assays [M]

- **MGO:** 2-methylquinoxaline by HPLC-DAD 315 nm (Atlantis T3), method of Glomb & Tschirnich 2001 [C]; calibration
  a* = 10⁷ × C1 (mM), R² 0.9999, **"experimental range was 0.002 to 3.2 mmol/L"** — read as the reporting floor
  0.002 mM; instrument RSD ±9.0 %; recovery 96 %.
- **Glucose, amino acids:** HPAEC-ED (CarboPac PA1), loss % = (1 − H_t / H_0) × 100 on peak heights (nC).
- **Acetic acid:** HPAEC-CD (IonPac AS-23), range 0.3-10 mM, RSD ±4.1 %, recovery 97 %.
- **Replicate noise:** Table 2 SDs are 0.0002-0.0011 mM (0.1-17 % RSD; at the high values 0.1-0.4 %) [!] — far below
  the method's own ±9 % instrument RSD; treat the printed SDs as injection repeats at best, and ±9 % as the floor of the
  real noise. Fig. 1 error bars ±0.0002-0.0010 mol/mol.

## 4. Data

### 4.1 Table 2 (verbatim) — MGO (mmol/L) after 0.5 h at 120 °C, 40 + 40 mM, 0.2 M citrate-phosphate

| pH | Pro-Glc | Gly-Glc | Lys-Glc | Arg-Glc |
|---|---|---|---|---|
| 5 | 0.0041 ± 0.0007 | 0.0050 ± 0.0007 | 0.0369 ± 0.0002 | 0.0086 ± 0.0003 |
| 6 | 0.0218 ± 0.0003 | 0.0346 ± 0.0007 | 0.1275 ± 0.0005 | 0.0528 ± 0.0004 |
| 7 | 0.0572 ± 0.0011 | 0.1649 ± 0.0003 | 0.3660 ± 0.0005 | 0.1484 ± 0.0004 |
| 7.4 | 0.0679 ± 0.0008 | 0.2717 ± 0.0010 | 0.6243 ± 0.0006 | 0.2250 ± 0.0011 |
| 8 | 0.0855 ± 0.0002 | 0.3785 ± 0.0004 | 0.7967 ± 0.0008 | 0.2867 ± 0.0006 |

### 4.2 Fig. 1 — MGO (mol / mol glucose), 0.2 M phosphate pH 7.4, 120 °C [F ±0.0007]

| t (h) | 0.5 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Lys | 0.0161 | 0.0220 | 0.0226 | 0.0193 | 0.0190 |
| Gly | 0.0098 | 0.0131 | 0.0146 | 0.0155 | 0.0142 |
| Arg | 0.0053 | 0.0068 | 0.0087 | 0.0095 | 0.0101 |
| Pro | 0.0031 | 0.0047 | 0.0075 | 0.0084 | 0.0090 |

### 4.3 Figs. 3 and 4 — pH 7.4 phosphate, 120 °C [F]

| t (h) | 0.5 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| glucose lost (%): Pro / Gly / Lys / Arg | 66 / 67 / 83 / 85 | 79 / 89 / 94 / 95 | 89 / 96 / 98 / 98 | 95 / 98 / 99 / 99 | 97 / 99 / 100 / 100 |
| amino acid lost (%): Pro / Gly / Lys / Arg | 3 / 11 / 25 / 30 | 12 / 16 / 35 / 42 | 15 / 22 / 40 / 48 | 14 / 27 / 44 / 54 | 16 / 45 / 44 / 57 |
| acetic acid (mmol/mmol glucose): Pro / Gly / Lys / Arg | 0.18 / 0.28 / 0.31 / 0.63 | 0.32 / 0.42 / 0.43 / 0.72 | 0.45 / 0.52 / 0.56 / 0.74 | 0.55 / 0.70 / 0.70 / 0.80 | 0.59 / 0.80 / 0.79 / 0.92 |

### 4.4 Arithmetic [Z] (`m4_arith.py`)

1. **pH slope of MGO at 0.5 h** (dex per pH unit): pH 5 → 7: Pro 0.57, **Gly 0.76**, Lys 0.50, Arg 0.62; pH 7 → 8: Pro
   0.18, **Gly 0.36**, Lys 0.34, Arg 0.29; pH 5 → 8: 0.44-0.63. Concave: steeper in acid-to-neutral.
2. **Amine identity at fixed glucose and pH (ratio to Gly):** Lys 7.4 / 3.7 / 2.2 / 2.3 / 2.1; Arg 1.7 / 1.5 / 0.90 / 0.83 /
   0.76; Pro 0.82 / 0.63 / 0.35 / 0.25 / 0.23 at pH 5 / 6 / 7 / 7.4 / 8.
3. MGO yield at pH 7, 0.5 h: 0.0014 (Pro) to 0.0092 (Lys) mol/mol glucose.
4. **Catalysis vs consumption at 0.5 h (Fig. 3):** glucose lost 66-85 % while amino acid lost 3-30 % (Gly 11 %): most
   glucose is lost without amine consumption.

## 5. What it can and cannot test in the model

**Can:**
- The pH law of MGO standing level from glucose + α-amino acid at 120 °C in a held buffer (five pH points), per amino
  acid. Gly-Glc is the trunk's own pair (`r_schiff` Glc + Gly), at 40 mM.
- Amine identity at fixed sugar and pH (four free α-amino acids); amine-catalysed glucose loss (Fig. 3).
- Glucose loss vs time at pH 7.4 / 120 °C (trunk observable), acetic acid yield (a 1-DG-exit marker).

**Cannot:**
- MGO formation and loss separately (standing level only); no fed MGO arm.
- Temperature dependence (120 °C only); thiol chemistry, MFT, NF; amine classes beyond free α-amino acids.
- The pH-series buffer differs from the time-series buffer (citrate-phosphate vs phosphate) [!]: at pH 7.4, 0.5 h,
  Table 2 vs Fig. 1 agree for Lys (0.0156 vs 0.0161 mol/mol) and Arg (0.0056 vs 0.0053) but not Gly (0.0068 vs 0.0098,
  ×1.4) or Pro (0.0017 vs 0.0031, ×1.8).

## 6. Conflicts and agreements with the repo's records

1. [!] Table 2 SDs are smaller than the stated instrument RSD (sec. 3); not real replicate noise.
2. [!] Table 2 vs Fig. 1 at pH 7.4 differ ×1.4-1.8 for Gly and Pro (different buffers; sec. 5).
3. Agrees in sign with `hayashi1986` (C3 fragment rising with pH) and with `kanzler2026` / `kanzler2021` that amines act
   largely catalytically on dicarbonyl chemistry (here, glucose loss far exceeds amine loss).
4. [!] Acetic acid up to 0.6-0.9 mol/mol glucose in 4 h at pH 7.4 is large; it is of the size `davidek2006a` reports for
   fed 1-DG at pH ≥ 6.5 (0.80 mol/mol), but here per glucose charged, with 97-100 % glucose lost. Plausible only if
   most lost glucose goes through 1-DG.

## 7. Provenance tiers (first reader only; second reader pending)

- Table 2: `direct_measurement` (printed), confidence **medium** (calibrated, recovery stated; SDs implausibly small).
- Figs. 1, 3, 4: `direct_measurement`, figure-read by eye at 90 dpi, confidence low-medium.
- Slopes and ratios (sec. 4.4): [Z], `derived`.

## M4 role (case 8660c005)

**Role fixed before reading:** part A evidence (Lys / Gly / Arg / Pro at fixed glucose and pH); band source for MGO
supply vs pH; measurement that would count: "MGO, glucose and acid vs pH for >= 2 amines"; P ~0.4.

- **The named measurement is partly met:** MGO vs pH for four amines (Table 2); glucose and acetic acid only at pH 7.4
  (Figs. 3-4), not vs pH.
- **Band (R3's C3 supply vs pH; sec. 5.2 "0 - 1 dex per pH unit, rising").** MGO standing level after 0.5 h at 120 °C,
  0.2 M citrate-phosphate (held), 40 mM: **0.50-0.76 dex/pH (pH 5 → 7), 0.18-0.36 (7 → 8)**, Gly 0.76 / 0.36. Inside the
  band, at its upper half on Hofmann's D2 limb (5 → 7). Transfer: α-amino acids (Gly is the trunk's amine; no Cys), 120
  vs 145 °C, a standing level (formation minus loss) rather than a supply constant, 0.2 M citrate-phosphate vs 0.5 M
  phosphate. A better-matched source than `hayashi1986` for the 5 → 7 limb.
- **Part A (D7-type evidence).** Within one class (free α-amino acids) at fixed glucose and pH, MGO varies ×4-30 (Lys
  2.1-7.4× Gly; Pro 0.23-0.82× Gly), and the spread shrinks with rising pH. All four are one class in the case's scheme
  (free α-amino acid, cysteine, peptide, protein lysine), so this is not a class contrast; it shows that **within-class
  spread is as large as the between-class contrasts D7 relies on** (Xia's ×5-7). It also shows amine catalysis of sugar
  loss without consumption (glucose −66 to −85 %, amino acid −3 to −30 % at 0.5 h), consistent with the case's "amine
  catalysis of dicarbonyl loss declared separately from amine consumption".
- **Decisive criterion:** none met (no thiol, no MFT).
