# Dworkin & Miller 2000 — EXTRACTION (relative free-aldehyde (open-chain) content of 12 aldoses from their second-order rate of reaction with urazole, 0.2 M + 0.2 M in D2O, pD ≈ 4.5, sealed NMR tubes, 25-120 °C, up to 2 years; calibrated to literature 13C-NMR values; 2-deoxyribose free aldehyde measured directly by 1H NMR at 40-100 °C)

### THE CORPUS'S BEST SET OF OPEN-CHAIN FRACTIONS AT 25 °C (mole % × 10³): **ribose 39.6, arabinose 19.2, xylose 7.97, lyxose 6.79, glucose 1.38** (i.e. 0.040 %, 0.019 %, 0.0080 %, 0.0068 %, 0.0014 %), ±~10 % relative; ribose : arabinose : xylose : lyxose : glucose = 29 : 14 : 5.8 : 4.9 : 1. **Temperature dependence is NOT measured for the aldoses**: the paper adjusts with glucose's ΔH = 9.98 kcal/mol (cited NMR) and states that the 100 °C / 25 °C rate ratios give ΔΔH vs ribose "within experimental error" for seven aldoses — i.e. a common van 't Hoff slope is assumed, not fitted. Only 2-deoxyribose has its own van 't Hoff line (log X = −1360/T + 1.82; ΔH 6.2 kcal/mol; 0.18 % at 25 °C).

**Source on disk:** `data/articles/dworkin2000.pdf` (7 pp., Elsevier publisher PDF, clean text layer). Prose from
`pdftotext -layout`; Table 1 re-read from a 200-dpi render (`scratch/lit8/dw_t1.png`) because the text layer turns "±"
into "9". Figs. 1-2 are log-log/van 't Hoff plots whose values are given in Table 1 or the text; not digitised.

## 0. Identity

| field | value |
|---|---|
| Title | "A kinetic estimate of the free aldehyde content of aldoses" |
| Authors | Jason P. Dworkin, **Stanley L. Miller** (corresponding) — UC San Diego |
| Venue | Carbohydrate Research **329** (2000) 359-365. Received 19 Jul 1999, accepted 20 Jun 2000 |
| DOI | PII S0008-6215(00)00204-4 as printed → **10.1016/S0008-6215(00)00204-4** ✔ matches the brief. **Identity verified.** |

## 1. Why it matters

Open-chain (free aldehyde) fractions set the reactive pool for amine condensation and enolisation; they bear on the
engine's sugar-reactivity ordering (pentose ≫ hexose; ribose > xylose) and on whether that ordering needs a fitted factor.
The brief asks for ribose, xylose, arabinose, lyxose, glucose vs temperature, with any van 't Hoff parameters.

## 2. Methods (pp. 360-361)

- **Reaction.** Urazole (1,2,4-triazole-3,5-dione) + aldose, **0.2 M each in D2O**, sealed NMR tubes, heated at 25, 40, 60,
  80, 100, 120 °C for up to 2 years; **pD ≈ 4.5** ("where the rate constants are independent of pH", refs. 13-14).
  Followed by 1H NMR of anomeric protons; fitted to a second-order reversible equation (A + B ⇌ C, Keq 0.5-3.1).
- **25 °C rates** measured directly for all aldoses except **glucose and mannose, which are Arrhenius extrapolations from
  60-120 °C** (p. 361). Some aldoses also measured at 100 °C (used only for ΔΔH, footnote 3).
- **Calibration.** log X_aldehyde = a + log k with a = 3.686, fitted so the urazole rates reproduce 13C-NMR free-aldehyde
  values (adjusted to 25 °C) for ribose (0.050 %), glucose (1.25 × 10⁻³ %), idose (0.0757 %), talose (0.0254 %).
  Assumes rate ∝ mole fraction of free aldehyde.
- **2-Deoxyribose.** 1.0 M in D2O; free-aldehyde 1H signal (10.4 ppm) integrated vs anomerics at 40, 60, 80, 100 °C.
- Ketoses (fructose, sorbose, tagatose) do not react detectably (< 5 %); Keq < 0.1 for ketoses and ketones.

## 3. Numbers

### Table 1 (p. 362), verbatim (% aldehyde × 10³)

| aldose | k (M⁻¹ s⁻¹) × 10⁹ at 25 °C | CD at 20 °C (Hayward & Angyal) | CD adjusted to 25 °C + NMR | **urazole (this work)** |
|---|---|---|---|---|
| idose | 160 ± 20 | 182 | 96.9 | 78.6 |
| **ribose** | 81.6 ± 0.7 | 46.9 | 25.0 | **39.6** |
| talose | 58 ± 2 | 33.6 | 19.2 | 28.0 |
| **arabinose** | 39.6 ± 0.5 | 33.3 | 17.9 | **19.2** |
| altrose | 36 ± 1 | 36.0 | 17.7 | 17.3 |
| **xylose** | 16.5 ± 0.7 | 15.6 | 13.6 | **7.97** |
| gulose | 14.4 ± 0.2 | 19.5 | 10.4 | 7.00 |
| **lyxose** | 14.0 ± 0.3 | 25.5 | 9.05 | **6.79** |
| galactose | 12.6 ± 0.6 | 17.0 | 8.30 | 6.10 |
| allose | 11.2 ± 0.4 | 10.9 | 5.80 | 5.41 |
| mannose | 6.2 ± 1 | 5.35 | 2.85 | 3.04 |
| **glucose** | 2.8 ± 0.4 | 2.22 | 1.18 | **1.38** |

Abstract: relative values "within about 10 %". Rate data: typically 24 points over 2 years; t½ 85 days (idose) to ~6.7
years (glucose) at 25 °C.

### Temperature dependence (all that the paper gives)

- **Glucose:** ΔH = **9.98 kcal/mol** (41.8 kJ/mol), from Maple & Allerhand's 13C NMR (ref. 10), used to move 20 °C CD
  values to 25 °C (×1.33).
- **Aldoses vs ribose:** footnote 3 — the 100 °C / 25 °C rate ratios for arabinose, lyxose, allose, galactose, glucose,
  mannose and talose relative to ribose give ΔΔH "within experimental error" (no numbers printed).
- **2-Deoxyribose (measured):** **log X = −1360/T + 1.82** → X(25 °C) = 0.18 %, ΔH = **6.2 kcal/mol**; consistent with
  0.3 % at 43.3 °C (Bauer 1986).
- Model compounds: acetaldehyde + urazole Keq 18.3 (N1), 1.3 (N1+N2); glycolaldehyde Keq 5.1 ± 0.2; glyceraldehyde ~791
  (total); formaldehyde ~600 (single adduct).

### Arithmetic (mine): extrapolation with a COMMON ΔH = 9.98 kcal/mol (glucose's), as the paper's footnote 3 licenses only loosely

Factor vs 25 °C: ×29.5 at 100 °C, ×126 at 145 °C (van 't Hoff, ΔH 41.8 kJ/mol).

| aldose | 25 °C (%) | 100 °C (%) | 145 °C (%) |
|---|---|---|---|
| ribose | 0.0396 | 1.17 | 5.0 |
| arabinose | 0.0192 | 0.57 | 2.4 |
| xylose | 0.00797 | 0.235 | 1.0 |
| lyxose | 0.00679 | 0.20 | 0.85 |
| glucose | 0.00138 | 0.041 | 0.17 |

These are **declared extrapolations** (one borrowed ΔH, 120 K beyond 25 °C; 2-deoxyribose's own ΔH is 6.2, not 10, so
the slope is sugar-dependent at least for deoxy sugars). Ratios between aldoses are T-independent under the common-ΔH
assumption, which is the only part the paper supports.

## 4. Numbers the repository can use

| quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|
| open-chain fraction, ribose | 0.0396 | mol % | 25 °C, D2O, pD 4.5 | Table 1 | measured (kinetic, NMR-calibrated) |
| arabinose / xylose / lyxose / glucose | 0.0192 / 0.00797 / 0.00679 / 0.00138 | mol % | same | Table 1 | same (glucose extrapolated from 60-120 °C) |
| relative open-chain (Rib = 1) | Ara 0.48, Xyl 0.20, Lyx 0.17, Glc 0.035 | — | same | Table 1 | ratio, ±~10 % |
| glucose ΔH (open ⇌ ring) | 9.98 | kcal/mol | — | cited (Maple & Allerhand 1987) | secondary |
| 2-deoxyribose van 't Hoff | log X = −1360/T + 1.82 | — | D2O, 40-100 °C | Fig. 2, text | measured |
| aldose-specific van 't Hoff | **NOT REPORTED** | — | — | — | — |

## 5. What this can and cannot anchor

**Can:** a measured ordering and ratio of reactive open-chain fractions among the pentoses and glucose at 25 °C, with a
justified (not proven) assumption that the ratios persist at 100 °C.
**Cannot:** aldose-specific ΔH; any value under Maillard conditions (pD 4.5 D2O, no phosphate, no amine); open-chain
fractions of ketoses (urazole does not react).

## 6. Flags

1. Glucose and mannose 25 °C rates are Arrhenius extrapolations, not measurements.
2. Calibration rests on four NMR values; ribose's urazole value (0.0396 %) is below the NMR anchor used (0.050 %).
3. The CD column and the urazole column disagree by up to 2× (xylose 13.6 vs 7.97; lyxose 9.05 vs 6.79).
4. Rizzi 2010 quotes ribose 0.05 %, arabinose 0.03 %, glucose 0.002 % (older Angyal values); use these instead.
5. D2O, not H2O; solvent isotope effects on ring-chain equilibria are not discussed.

## 7. Use in the model

- **Sugar-reactivity ordering without a fitted factor.** Rib : Ara : Xyl : Glc open-chain = 1 : 0.48 : 0.20 : 0.035. This
  matches the ordering of Rizzi 2010's amine-free reductant production (1 : 0.56 : 0.38 : 0.08) and Mikami 2017's NF yields
  (ribose > xylose ≈ arabinose) in direction, and gives a hexose : pentose open-chain ratio of ~1 : 6-29, consistent with
  the engine's structural "pentose above hexose" claim. It does NOT explain Garrett 1969's acid-furfural ranking (Rib > Xyl
  ≈ Lyx > Ara): arabinose is 2.4× more open than xylose but dehydrates slower, so furfural formation is not open-chain-limited.
- **Temperature:** a common ΔH ≈ 42 kJ/mol for ring opening means the open fraction rises ~4× from 100 to 145 °C
  (×29.5 → ×126 vs 25 °C). Any engine step written on the open-chain sugar would carry ~42 kJ/mol of its apparent barrier
  from this equilibrium — relevant when splitting apparent Ea's (e.g. Lau/Pandit's 81-93 kJ/mol for ribose → NF) into
  ring-opening and chemical steps. Declared extrapolation, not a measurement.
- No direct bearing on the furfural barrier swap or B61.
