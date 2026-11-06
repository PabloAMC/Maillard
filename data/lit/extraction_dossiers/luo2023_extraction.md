# Luo, Zhu, Peng, Cui, Huang, Xu & Ho 2023 (J. Agric. Food Chem. 72:657-669, issue 2024) — EXTRACTION (file stem `luo2023`)
**DOI:** 10.1021/acs.jafc.3c06250 (printed on every page footer; Crossref 2026-11-05: same title, authors Yue Luo, Siyue Zhu, Jie Peng, Heping Cui, Qingrong Huang, Baojun Xu, Chi-Tang Ho, JAFC 72(1) 657-669, online 2023-12-18, issue 2024-01-10) — matches the brief.
### *"Feasibility Study of Amadori Rearrangement Products of Glycine, Diglycine, Triglycine, and Glucose as Potential Food Additives for Production, Stability, and Flavor Formation"*
### Purified Gly-ARP, Digly-ARP, Trigly-ARP (N-(1-deoxy-D-fructos-1-yl)-glycine / -diglycine / -triglycine, ≥ 96 %), **1 mg/mL in deionized water, initial pH 3 / 5 / 7 / 9, heated at 100 / 130 / 160 °C for 30-180 min** (Fig. 5: % ARP left and % amine released, LC-MS/MS); dicarbonyls by OPD-HPLC (Table S1, not on disk); volatiles from ARPs vs glucose + amine (20 mM, pH 7, 100-160 °C, 1-2 h; Tables 1-2). **Gly-ARP left after 3 h: ~75 % (100 °C, pH 3-7), 43 % (pH 9); 34-54 % (130 °C, pH 3-7), 4 % (pH 9); ~0 at 160 °C. Peptide ARPs are less stable: 53-73 % at 100 °C, < 50 % within 30 min at 130 °C.** Decay is NOT first order (a fast first 30 min, then a slow tail).

**Key and year.** Online 18 Dec 2023, issue JAFC 72(1), 2024 (the PDF's citation line reads 2024). The stem `luo2023`
follows the file name, the brief and case M4's "Luo et al. 2023", as `aylonparra2026` keeps its online year. `luo2024` in
the repo is a different paper (Luo, Tian et al. 2024, a sulfur-volatile review).

**Read 2026-11-05 (first reader) from `data/articles/luo2023.pdf`** (13 pp., ACS publisher PDF, clean text layer).
Prose from `pdftotext -raw`. **Fig. 5a and 5b digitised by eye** from a 300 dpi render of p. 6 (scratchpad
`m4img/luo2023_300-06.png`, crop `m4img/luo_f5ab.png`; 2.52 px per %): **read error ±2 %** (marker centres; the circle
and square series overlap at 100 °C). Panels 5c-5i read at 130 dpi (±4 %). Table 2's furfural row read from a 130 dpi
render of p. 9 (rotated table; the text layer's column order is scrambled). Figs. 1, 2, 4, 6, 7 not digitised.
Arithmetic: scratchpad `m4tools/m4_arith.py`. SI (NMR, Table S1 dicarbonyls, Table S2 volatiles, LODs) not on disk.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure by me ·
**[Q]** qualitative statement by the authors · **[C]** cited from another paper · **[P]** proposed, not shown ·
**[Z]** derived by me, never printed · **[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Yue Luo, Siyue Zhu, Jie Peng, Heping Cui, Qingrong Huang, **Baojun Xu** (corresponding), **Chi-Tang Ho** (corresponding) |
| Affiliation | (Author Information, p. 667) Dept. of Food Science, Rutgers University, New Brunswick NJ (Luo, Zhu, Peng, Huang, Ho); BNU-HKBU United International College, Zhuhai (Xu); Jiangnan University, Wuxi (Cui) |
| Venue | *J. Agric. Food Chem.* **72**(1) (2024) 657-669 |
| DOI | 10.1021/acs.jafc.3c06250 |
| Dates | received 3 Sep 2023; revised 29 Nov 2023; accepted 5 Dec 2023; published 18 Dec 2023 |

**Laboratory.** Same group as `luo2022` (Xu / Ho; Heping Cui of Jiangnan now a co-author). See `luo2022` sec. 0 on the
relation to `xia2022` (Jiangnan); owner call for leave-one-laboratory-out.

## 1. What this paper answers for the model

| model need | what Luo 2023 gives | tag |
|---|---|---|
| amino-acid ARP decay vs T and pH | Fig. 5a-c: Gly-ARP left (%) at 30-180 min, 100 / 130 / 160 °C, pH 3 / 5 / 7 / 9 | [F] / [M] |
| peptide ARP decay | Fig. 5d-i: Digly-ARP, Trigly-ARP; text: 52.61-72.43 % and 60.42-73.42 % left after 3 h at 100 °C | [F] / [M] |
| amine returned per ARP lost | Fig. 5 light series: glycine (or peptide) released, % of initial ARP | [F] |
| bond-cleavage site (Trigly-ARP, 100 °C) | triglycine 10.02-12.44 % (30 min) → 12.71-24.26 % (3 h); Digly-ARP 4.07-4.88 %; Gly-ARP 0 → 1.89-4.37 %; Trigly-ARP ~85 % (30 min), 48.90-74.63 % (3 h) | [M] |
| amine stability alone | glycine and diglycine keep "over 90, 80, and 75%" at 100 / 130 / 160 °C; triglycine 75 / 40 / 10 % "in the worst case" | [M] |
| ARP preparation yields | freeze-dried, 70 °C: Gly-ARP 54.24 %, Digly-ARP 57.64 %, Trigly-ARP 98.00 %; two-step vacuum 18.30 / 34.44 / 78.70 % | [M] |
| dicarbonyls from ARPs | Table S1 (not on disk); text: 3-DG main in acid, GS and 1-DG at neutral / alkaline; total dicarbonyls at pH 3 vs pH 9 (100 °C) ×1.37 / 7.89 / 8.45 (Gly- / Digly- / Trigly-ARP) | [M] / [Q] |
| furfural from glucose + amine and from ARPs | Table 2 (sec. 4.3), µg/L | [M] |
| thiol, MFT, NF | **none** | — |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Stability pots | **1 mg ARP** (or an equimolar amount of the amino compound) in **1 mL deionized water** [Z: Gly-ARP 237 g/mol → 4.2 mM; Digly-ARP ~3.4 mM; Trigly-ARP ~2.8 mM], pH adjusted to 3.0 / 5.0 / 7.0 / 9.0; heated at 100 / 130 / 160 °C for 30, 60, 90, 120, 180 min; ice-water quench |
| Stability vessel | **not stated** for this experiment [!] (the volatile series uses 20 mL pressure vessels) |
| Volatile pots | pure ARP 0.02 M, or glucose + gly / digly / trigly 0.02 M each, Milli-Q water, pH 7.0 (6 N NaOH), 20 mL pressure vessels with stirring bars, oil bath 100 / 130 / 160 °C, 60 and 120 min; triplicate |
| pH | label pH; no post-heating pH reported |
| Headspace | not stated |
| Replicates | volatiles: triplicate; stability series: not stated, no error bars in Fig. 5 |

## 3. Assays [M]

- **ARPs and amino compounds:** HILIC-LC-MS/MS, MRM as in `luo2022`; purified standards (≥ 96 %; NMR in SI). No LOD,
  recovery or replicate statistic printed in the main text.
- **Dicarbonyls:** OPD 0.5 %, overnight dark; HPLC-DAD 315 nm (Innoval C18); standards GS, 1-DG, 3-DG, GO, MGO; "standard
  curve equations, limits of detection, and quantification" in the SI (not on disk).
- **Volatiles:** HS-SPME (DVB/CAR/PDMS), GC-MS DB-Wax 60 m; pyrazines, furfural, benzaldehyde by external standard
  curves; the rest semi-quantified vs 1,2-dichlorobenzene with RF = 1.
- [!] 1-DG as a commercial / available OPD standard is unusual; its source is not given.

## 4. Data

### 4.1 Fig. 5a / 5b — Gly-ARP left (% of initial) [F ±2]

| | 30 | 60 | 90 | 120 | 180 min |
|---|---|---|---|---|---|
| 100 °C pH 3 | 81.7 | 77 | 76.2 | 76 | 76 |
| 100 °C pH 5 | 77.4 | 77 | 76 | 76 | 76 |
| 100 °C pH 7 | 79.8 | 76 | 75.5 | 74.5 | 72.6 |
| 100 °C pH 9 | 63.9 | 59.5 | 55.6 | 52.4 | 42.9 |
| 130 °C pH 3 | 75.0 | 69.4 | 64.3 | 58.7 | 54.0 |
| 130 °C pH 5 | 61.5 | 61.1 | 58.3 | 58.3 | 49.2 |
| 130 °C pH 7 | 67.5 | 60.3 | 54.0 | 40.5 | 34.1 |
| 130 °C pH 9 | 26.2 | 12.3 | 9.1 | 5.2 | 4.0 |
| 160 °C (5c), pH 3 / 5 / 7 / 9 | ~29 / ~22 / ~30 / ~0 | ~0 / ~8 / ~17 / 0 | 0 | 0 | 0 |

Glycine released (% of initial ARP) at 100 °C: pH 3 and 5 ≈ 5 / 7.5 / 10 / 11.5 / 13.5; pH 7 11.1 / 17 / 20.6 / 21.4 / 24.6;
pH 9 27.8 / 36.5 / 41.3 / 43.7 / 46.0 (30 … 180 min). At 130 °C 23-53 % at 30 min, 44-50 % at 180 min; at 160 °C a
plateau at 33-45 % from 30 min.

Text: Gly-ARP "dropped to around 75%, 45%, and 0%" after 3 h at 100 / 130 / 160 °C [M] (the pH-pooled statement; Fig. 5
shows pH 9 well below at each temperature).

### 4.2 Fig. 5d-i — peptide ARPs [F ±4 unless printed]

- Digly-ARP: 100 °C, 180 min: 52.61-72.43 % [M]; 130 °C: 20-45 % at 30 min, 4-25 % at 180 min; 160 °C: 5-15 % at 30 min,
  ~0 from 90 min. Diglycine released at 130-160 °C: 30-60 %.
- Trigly-ARP: 100 °C, 180 min: 60.42-73.42 % [M]; 130 °C: 17-53 % at 30 min, 0-22 % at 180 min; 160 °C: 0-3 % at 30 min.
  Triglycine released: 10-25 %.

### 4.3 Table 2 — 2-furfural (µg/L, mean ± SD, n = 3), pH 7 initial, 20 mM [M; read from the page image]

| system | 130 °C 1 h | 130 °C 2 h | 160 °C 1 h | 160 °C 2 h |
|---|---|---|---|---|
| glucose + glycine | 0.12 ± 0.02 | 0.52 ± 0.48 | 1.18 ± 1.00 | 15.99 ± 0.77 |
| glucose + diglycine | 0.15 ± 0.04 | 1.82 ± 1.36 | 16.94 ± 1.86 | 19.35 ± 2.79 |
| glucose + triglycine | 0.16 ± 0.09 | 0.49 ± 0.06 | 16.40 ± 1.28 | 19.60 ± 2.52 |
| Gly-ARP | 0.93 ± 0.56 | 3.58 ± 0.58 | 24.85 ± 2.75 | 35.19 ± 3.11 |
| Digly-ARP | 1.22 ± 0.34 | 0.47 ± 0.29 | 19.34 ± 1.39 | 22.34 ± 2.40 |
| Trigly-ARP | 0.43 ± 0.51 | 1.54 ± 0.26 | 10.31 ± 1.79 | 12.23 ± 0.61 |

At 100 °C furfural is blank or ≤ 0.14 µg/L in every system. Other furans and the pyrazines (Tables 1-2) are not
transcribed (no engine species; semi-quantitative RF = 1 except the pyrazines).

### 4.4 Arithmetic [Z] (`m4_arith.py`)

1. **Gly-ARP apparent first-order k (min⁻¹):** 100 °C, pH 3 / 5 / 7 / 9: k(0-30) 0.0067 / 0.0085 / 0.0075 / 0.0149;
   k(30-180) 0.00048 / 0.00012 / 0.00063 / 0.0027; k(0-180) 0.0015 / 0.0015 / 0.0018 / 0.0047. 130 °C: k(0-30) 0.0096 /
   0.0162 / 0.0131 / 0.0446; k(0-180) 0.0034 / 0.0039 / 0.0060 / 0.0179. **Not first order**: at 100 °C the 0-30 min
   rate is 10-70× the 30-180 min rate.
2. **Apparent Ea (100 → 130 °C), Gly-ARP:** 34 / 40 / 51 / 56 kJ/mol (pH 3 / 5 / 7 / 9, on k(0-180)); 15-46 on k(0-30).
3. **Amine returned per ARP lost** (glycine released / ARP lost, 180 min): 100 °C pH 3 0.54, pH 5 0.58, pH 7 0.90, pH 9
   0.81; 130 °C ~1.0 (pH 3), ~0.86 (pH 5), ~0.70 (pH 7), ~0.48 (pH 9).
4. Peptide ARPs: 100 °C k(0-180) 1.7-3.6 × 10⁻³; 130 °C k(0-30) ≥ 0.023; 160 °C k(0-30) ≥ 0.063 min⁻¹.
5. Furfural yield, glucose + glycine, 160 °C, 2 h: 16 µg/L = 0.17 µM, 8 × 10⁻⁶ mol/mol glucose.

## 5. What it can and cannot test in the model

**Can:**
- **The decay of a hexose amino-acid ARP (Gly-ARP, the trunk's DFG class) vs temperature and pH**, fed, with the amine
  return measured in the same pots. Gly-ARP is the Amadori compound of the trunk's own pair (Glc + Gly).
- A peptide-ARP vs amino-acid-ARP stability contrast at fixed conditions (part A).
- Furfural from glucose + Gly and from fed Gly-ARP at pH 7, 130-160 °C (small; a check on `r_glc_fur` / `r_arp_tdp`
  size at pH 7, triplicate, external standard).

**Cannot:**
- A clean rate constant: the fast-then-slow shape (sec. 4.4) is not first order, and its cause (reversible ARP
  formation, come-up in an unstated vessel, pH drift in unbuffered water, or assay) is not resolved.
- Dicarbonyl yields (Table S1 not on disk). Thiol chemistry, Cys-ARPs, MFT, NF.
- Replicate noise for Fig. 5 (no error bars).

## 6. Conflicts and agreements with the repo's records

1. **[!] Temperature dependence far below the repo's hexose-ARP barriers.** Apparent Ea 34-56 kJ/mol for Gly-ARP decay
   (100 → 130 °C) vs Martins 2005 (97.1 / 107.3 kJ/mol, case M4 sec. 2.7), Zhang 2026 (85.7) and the case's band (85.7-107).
   On the fast 0-30 min window it is 15-46 kJ/mol. With a non-first-order shape and no vessel stated, a heat-transfer
   artefact cannot be excluded.
2. **Size at 100 °C, pH 5:** Gly-ARP k(0-180) 1.5 × 10⁻³ min⁻¹, ×0.31 `deng2021`'s Met-glucose ARP 0.0048 (100 °C, label
   pH 5.5); k(0-30) 8.5 × 10⁻³, ×1.8. Brackets Deng's value.
3. **Amine return:** 0.54-1.0 (pH 3-9, 100-130 °C) vs the case's band 0.71-1.0 (Martins 2003b, Davidek 2002 / 2005,
   Deng 2021). Extends the band downward at pH 3-5 / 100 °C (0.54-0.58) and pH 9 / 130 °C (0.48).
4. Agrees with `luo2022` and `xia2022` that peptide ARPs degrade faster than Gly-ARP above 100 °C; agrees with
   Martins / Davidek that 3-DG dominates ARP degradation in acid and 1-DG / glucosone at neutral-alkaline pH [Q].
5. [!] Text "75%, 45%, and 0%" is a pooled statement; pH 9 is 43 / 4 % at 100 / 130 °C (Fig. 5).

## 7. Provenance tiers (first reader only; second reader pending)

- Fig. 5 values: `direct_measurement`, figure-digitised, LC-MS/MS with purified standards, confidence **low-medium** (no
  error bars, vessel unstated, non-first-order shape).
- Text percentages and Table 2: `direct_measurement` (printed); Table 2 with SDs (n = 3), confidence medium.
- k, Ea and return ratios in sec. 4.4: [Z], `derived`.

## M4 role (case 8660c005)

**Role fixed before reading:** band source for peptide-ARP decay vs T (100 / 130 / 160 °C); measurement that would count:
"ARP decay at >= 2 T"; P ~0.4.

- **The named measurement is met:** Gly-, Digly- and Trigly-ARP decay at three temperatures and four pH.
- **What it supplies.**
  - *Peptide-ARP decay vs T (part A, the role as fixed):* Digly / Trigly-ARP k(0-180) 1.7-3.6 × 10⁻³ min⁻¹ at 100 °C;
    ≥ 0.023 min⁻¹ (0-30 min) at 130 °C; ≥ 0.063 at 160 °C; water, label pH 3-9, ~3 mM. Peptide ARPs decay faster than
    Gly-ARP at ≥ 130 °C.
  - *Bearing on the part G band "Cys-hexose Amadori → 1-DG + Cys, decay" (sec. 5.2: Deng 2021's 0.0048 min⁻¹ at 100 °C,
    ×0.3-3, Ea 85.7-107):* the Gly-ARP decay at 100 °C (1.5 × 10⁻³ to 8.5 × 10⁻³ min⁻¹ by window) sits inside the ×0.3-3
    width; its apparent Ea (34-56 kJ/mol) is **below** the band's 85.7-107. Transfer: same class (hexose α-amino-acid
    ARP) but glycine, not cysteine; water, unbuffered; vessel unstated.
  - *Amine returned per decay (sec. 5.2 band 0.71-1.0):* 0.54-1.0, lower than the band at pH 3-5 / 100 °C and pH 9 / 130 °C.
- **Unexpected (not decisive):** the decay is fast-then-slow at every temperature and pH, so a single first-order
  constant with an Arrhenius law does not describe a fed hexose ARP here; the band's Ea range is not supported by this
  paper's apparent values. A required-agreement sensitivity at a lower Ea (~40-55 kJ/mol) would have a measured basis.
- **Decisive criterion:** none met (no thiol, no hexose MFT).
