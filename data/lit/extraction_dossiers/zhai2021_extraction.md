# Zhai et al. 2021 — EXTRACTION (TTCA degradation kinetics, with and without extra xylose)
### 2-Threityl-thiazolidine-4-carboxylic acid 10 mM at pH 5.5 / 6 / 7 / 8, 100 / 120 / 140 C, 20-120 min.

**Source on disk:** `data/articles/zhai2021.pdf` (owner's download, 2026-09-07). Read-only extraction,
plan item W6 (a cysteine-xylose intermediate the core can charge by name).

## 0. Identity

| field | value |
|---|---|
| Title | "Degradation of 2-Threityl-Thiazolidine-4-Carboxylic Acid and Corresponding Browning Accelerated by Trapping Reaction between Extra-Added Xylose and Released Cysteine during Maillard Reaction" |
| Authors | Yun Zhai, Heping Cui, Khizar Hayat, Xiaoming Zhang, Chi-Tang Ho et al. (Jiangnan) |
| Venue | J. Agric. Food Chem. 2021, 69, 10648-10656 |
| DOI | 10.1021/acs.jafc.1c03727 |

## 1. Methods

- TTCA synthesised from Xyl (1.0 g) + Cys (0.8 g), pH 7.4, 90 C, then purified (98 %); mass/NMR in the SI.
- Aqueous TTCA 10 mmol/L with or without Xyl 10 mmol/L, pH 5.5 / 6 / 7 / 8 (NaOH/HCl), oil bath at
  100 / 120 / 140 C for 20 / 40 / 60 / 80 / 100 / 120 min, ice quench; TTCA, ARP, Cys, Xyl by
  HPLC-ELSD; 13C5-xylose tracer for the trapping reaction; dicarbonyls (3-DX, 1-DX, MGO, GO) at 100 C
  vs pH (Fig. 5); browning A420.

## 2. The kinetics (sec. "Reaction Kinetics" and Fig. 3a; TTCA alone, pH 7)

Zero-order fits c = c0 - k t (c in mmol/L, t in min):

| T (C) | fitted line | k (mmol/L/min) | R2 | first-order equivalent at 10 mM (/min) |
|---:|---|---:|---:|---:|
| 100 | y = -0.0271 x + 10.331 | 0.0271 | 0.9516 | 0.0027 |
| 120 | y = -0.0651 x + 9.9718 | 0.0651 | 0.9949 | 0.0065 |
| 140 | y = -0.0813 x + 9.3375 | 0.0813 | 0.9802 | 0.0081 |

Authors: Ea of TTCA degradation **80.99 kJ/mol** (Arrhenius on the three k). ⚠ Re-derived from the
three printed constants the slope gives ~35 kJ/mol (100 -> 120 C: 53; 120 -> 140 C: 15), so the printed
Ea is not reproducible from the printed k; the dossier carries both and flags the discrepancy.
Degradation is faster at higher pH (Fig. 4b: pH 7 and 8 notably faster than 5.5 and 6 at 100 C) and
is accelerated by extra xylose through re-formation of TTCA from released cysteine.

## 3. Against the core

The sulfur lane's `k_ttca_deg` is a FITTED coordinate (B9: log10 k at 145 C = -1.698, i.e. 0.020 /min,
route barrier 64.1 kJ/mol), giving 0.0022 /min at 100 C and 0.016 /min at 140 C. Zhai's first-order
equivalents: 0.0027 and 0.0081. The core matches within 1.25x at 100 C and overshoots 2x at 140 C; the
core's 100 -> 140 C ratio (7.3x) is steeper than the measured 3.0x. **These three constants are
primary evidence (a measured rate at three temperatures on a species the lane carries) and qualify as
FIT rows under the owner's rule** — the first direct temperature series on the TTCA step in the corpus.
Recorded for the next sulfur wave; not installed here.

## 4. What else the paper holds

Fig. 5: 3-deoxypentosone, 1-deoxypentosone, MGO and GO vs time at 100 C and four pH values (figure
only); Table 1: LC-MS/MS fragments of 13C5-TTCA. The pH dependence of TTCA degradation is the missing
input for Wang 2026's series (pH unstated there).
