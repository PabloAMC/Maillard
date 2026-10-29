# Yu & Tang 2016 — EXTRACTION (L-ascorbic acid 0.1 M + L-cysteine 0.1 M in 0.2 M phosphate, pH 8.0 (NaOH-adjusted), sealed 15 mL pressure vials, oil bath 110 / 125 / 135 / 140 / 145 / 155 °C, 10-120 min; A420, A294, residual ASA (UV 243 nm), residual Cys (DNFB-HPLC))

### TIME-RESOLVED CYSTEINE AND ASCORBIC ACID LOSS AT SIX TEMPERATURES IN ONE BUFFERED POT — but only in figures. The paper prints one rate law and one Ea: **browning (A420) Ea = 114.33 kJ/mol (125-155 °C)**, with k(125-155 °C) = 0.0466-0.539 (units not printed; not reproducible from Fig. 1A, see Flag 1). From the digitised Fig. 2 (mine): **Cys loss k ≈ 0.010 / 0.039 / 0.12 / 0.16 / 0.25 min⁻¹ at 110 / 125 / 135 / 140 / 145 °C, Ea ≈ 109-124 kJ/mol**; **ASA loss k ≈ 0.0012 / 0.0031 / 0.0042 / 0.0070 / 0.0117 min⁻¹ at 125-155 °C, Ea ≈ 98-109 kJ/mol**. No H2S, no volatiles, no pH after heating.

**Source on disk:** `data/articles/Yu2016.pdf` (8 pp., publisher PDF, clean text layer). Prose and Tables
1-2 from `pdftotext -layout`. **Figs. 1A and 2A/2B digitised** from 300 dpi crops (scratchpad
`lit5/yu_f1a-5.png`, `yu_f2a-6.png`, `yu_f2b-6.png`) and tagged **[D]**. Fits in scratchpad
`lit5/yu_fit.py` (maillard env). No SI.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of Non-Enzymatic Browning Reaction from the L-Ascorbic Acid/L-Cysteine Model System" |
| Authors | **Ai-Nong Yu** (corresponding, anyufoodchem@hotmail.com) and Le-Pan Tang — Hubei University for Nationalities, Enshi, China |
| Venue | Czech J. Food Sci. **34** (6), 503-510, 2016. Received 12 Mar 2016, accepted 16 Nov 2016 |
| DOI as printed | **10.17221/86/2016-CJFS** ✔ matches the brief. **Identity verified.** |
| Funding | NSFC 31360408; Hubei discipline group 2015-5 |
| Same lab, in corpus | `yu2010_extraction.md` (Yu & Zhang 2010: ASA + Cys 0.1 M each, phosphate pH 5-9, 141 °C, 2 h, volatiles as TIC areas) — **same pot recipe at pH 8**, so this paper is that pot's kinetic companion. Yu, Tan & Wang 2012a (Food Chem. 132:1316, labelled sulfur-volatile mechanism) and Tang, Yang & Yu 2015 (Food Science 36:63, volatile formation kinetics, in Chinese) are cited, **not on disk**. (`yu2012_extraction.md` in the repo is a different paper: Yu, Meng, Ramaswamy & Boye, extrusion.) |

## 1. Why it matters

Ascorbic acid is a non-sugar C6 reductone/dicarbonyl source (via dehydroascorbic acid, 2,3-diketogulonic
acid, and decarboxylation to pentose-like C5 fragments — the "PT" of this paper). With cysteine it is a
Maillard-type sulfur pot. The engine's cysteine sinks are thermolysis (`k_cys_h2s`, Zheng & Ho 1994,
2.8 × 10⁻³ min⁻¹ at 145 °C, Ea 133) and `k_cys_thermal` (fitted, 7.5 × 10⁻³ min⁻¹ at 145 °C, Ea 55.1)
plus the carbonyl reactions. This paper gives the total Cys loss rate in a carbonyl-rich pot at six
temperatures.

## 2. Methods as they matter to a model

- **Charge:** ASA 1.0 mmol + Cys 1.0 mmol in 10 mL 0.2 M phosphate (pH 8) → **0.10 M each** (ASA
  17.6 mg/mL, Cys 12.1 mg/mL nominal). **pH adjusted to 8.0 with NaOH after dissolution** (ASA is a
  0.1 M acid; NaOH demand not reported). Set "according to Yu & Zhang 2010".
- **Vessel:** 15 mL Synthware pressure glass vial (headspace ~5 mL, air — "sealed ... for blocking off
  oxygen" but not degassed). Oil bath with stirring; **110, 125, 135, 140, 145, 155 °C**; 10, 20, 30, 40,
  60, 90, 120 min; cold-water quench. Heat-up time not reported.
- **pH after heating: not measured.**
- **ASA:** 50 µL → 50 mL water, A243, external standard (0-40 µg/mL), RSD ≤ 5.4 %. **Not selective**:
  A243 includes any UV-absorbing product (Flag 2).
- **Cys:** pre-column DNFB derivatisation, C18 HPLC, 360 nm, external standard, RSD ≤ 11.9 %. Measures
  free thiol-amine Cys (cystine and Cys adducts not measured).
- **A420 / A294:** room temperature, diluted 0-25×, back-calculated; RSD ≤ 14.5 %.
- **Kinetics:** zero/first/second-order fits "by linear regression of [A], ln[A], 1/[A] versus t";
  Arrhenius on ln k vs 1/T. Triplicates.

## 3. Printed numbers

### Table 1 (verbatim) — "rate constants of browning development at 420 nm"

| T (°C) | 110 | 125 | 135 | 140 | 145 | 155 |
|---|---|---|---|---|---|---|
| k | — (rejected, not first order) | 0.0466 | 0.126 | 0.170 | 0.260 | 0.539 |
| R² | — | 0.9299 | 0.9889 | 0.9695 | 0.9319 | 0.9185 |

**Arrhenius (p. 506): ln k = −13 751.79/T + 31.54, R² 0.9948 → Ea = 114.33 kJ/mol.** I reproduce Ea
114.5 and intercept 31.57 from Table 1. Units of k **not printed** (min⁻¹ implied by the time axis).

### Table 2 (verbatim) — Pearson r (10-120 min)

| T (°C) | 110 | 125 | 135 | 140 | 145 | 155 |
|---|---|---|---|---|---|---|
| A294 vs A420 | 0.67 | 0.76 | 0.87 | 0.95 | 0.96 | 0.92 |
| A420 vs residual ASA | −0.72 | −0.90 | −0.97 | −0.99 | −0.99 | −0.95 |

### Text facts

- Cys alone (no ASA), same system, 140 °C, 2 h: A420 = 0.08 ("light yellow") — Cys self-browning
  negligible.
- "Almost complete degradation [of Cys] ... after 30 min at a temperature above 135 °C."
- Rogacheva 1999 (cited): ASA + Glu / Gly / Lys melanoidin Ea 136.96 / 106.07 / 99.79 kJ/mol.
- Below 110 °C "the non-enzymatic browning reaction did not occur at all or the reaction was slow";
  Tang 2015 (cited) found no flavour volatiles at 110 °C.
- Proposed rate law: d[browning]/dt = k1[ASA]{1 + K3K4[Cys]/k2} → k1[ASA] once Cys is depleted.

## 4. Figures digitised

### Fig. 2B — residual Cys (mg/mL) [D], read error ± 0.1 mg/mL (≈ 0.8 mM)

| T \ t (min) | 10 | 20 | 30 | 40 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|
| 110 °C | 9.93 | 8.73 | 7.76 | 6.88 | 5.35 | 3.78 | 3.41 |
| 125 °C | 8.12 | 5.27 | 2.94 | 1.75 | 0.78 | 0.27 | 0.12 |
| 135 °C | 4.88 | 1.49 | 0.41 | 0.16 | ~0.08 | ~0.05 | ~0.05 |
| 140 °C | 3.94 | 0.78 | ~0.1 | ~0.05 | ~0.05 | ~0.05 | ~0.05 |
| 145 °C | 1.53 | ~0.12 | ~0.05 | … | | | |
| 155 °C | ~0.12 | ~0.05 | … | | | | |

Nominal start 12.1 mg/mL; **no t = 0 point is plotted**; the 110 °C pot has already lost 18 % by 10 min.

### Fig. 2A — residual ASA (mg/mL, by A243) [D], read error ± 0.15 mg/mL

| T \ t (min) | 10 | 20 | 30 | 40 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|
| 110 °C | 18.4 | 18.0 | 18.0 | 18.4 | 17.8 | 17.6 | 17.3 |
| 125 °C | 18.2 | 17.5 | 17.5 | 17.2 | 16.8 | 16.2 | 15.8 |
| 135 °C | 17.9 | 17.4 | 16.6 | 16.5 | 15.4 | 13.3 | 13.0 |
| 140 °C | 17.6 | 17.3 | 16.4 | 15.7 | 15.0 | 12.8 | 11.1 |
| 145 °C | 17.3 | 16.5 | 15.1 | 13.5 | 12.0 | 9.8 | 8.2 |
| 155 °C | 16.7 | 14.0 | 12.0 | 10.4 | 7.6 | 5.6 | 4.7 |

Nominal start 17.6 mg/mL; the 110 °C series reads 18.0-18.4, **above** nominal (Flag 2). Marker
assignment at 10-30 min where series overlap is ± 1 series.

### Fig. 1A — A420 [D], read error ± 0.2

| T \ t (min) | 10 | 20 | 30 | 40 | 60 | 90 | 120 |
|---|---|---|---|---|---|---|---|
| 110 °C | ~0.05 | ~0.05 | ~0.05 | ~0.05 | ~0.1 | ~0.1 | 0.2 |
| 125 °C | 0.1 | 0.2 | 0.2 | 0.3 | 0.85 | 2.2 | 3.6 |
| 135 °C | 0.1 | 0.3 | 0.85 | 2.0 | 4.8 | 8.7 | 12.0 |
| 140 °C | 0.15 | 0.5 | 2.5 | 4.3 | 7.7 | 12.5 | 15.5 |
| 145 °C | 0.2 | 2.0 | 5.8 | 8.2 | 12.6 | 16.7 | 18.3 |
| 155 °C | 1.9 | 7.8 | 12.6 | 15.1 | 18.2 | 18.2 | 18.4 |

## 5. Arithmetic (mine)

1. **Cys loss, apparent first order** (ln-linear over the points above the ~0.1 mg/mL floor):
   **k = 0.010 (110 °C, curved, decelerating), 0.039 (125), 0.115 (135), 0.162 (140), ~0.25 (145, two
   points, floor-limited) min⁻¹; 155 °C ≥ 0.46 min⁻¹** (nominal 12.1 → 0.12 mg/mL in 10 min).
   **Ea = 124 kJ/mol** (110-145 °C ln-fit) or **109 kJ/mol** (125-145 °C, from the nominal start to the
   10 min point). Take **Ea ≈ 110-125 kJ/mol, k(145 °C) ≈ 0.2-0.25 min⁻¹** at 0.1 M ASA, pH 8 (initial).
   The 110 °C curve flattens towards ~3 mg/mL (k falls from 0.012 to 0.0075 min⁻¹), so the 110 °C point
   is not a clean first-order constant.
2. **Against the engine:** at 145 °C the engine's two cysteine-alone sinks sum to ~0.010 min⁻¹
   (`k_cys_h2s` 0.0028 + `k_cys_thermal` 0.0075). Cys disappears **~25× faster** here, so **this is not
   thermolysis**: Cys is consumed by ASA/dehydroascorbate-derived carbonyls (thiazolidine/Strecker-type
   chemistry). **Do not use these numbers for `k_cys_h2s`.** They are a total-loss observable for an
   ASA + Cys pot.
3. **ASA loss, first order, 10-120 min:** k = 4.9e-4 (110, noise-level), 1.18e-3 (125), 3.1e-3 (135),
   4.2e-3 (140), 7.0e-3 (145), 1.17e-2 (155) min⁻¹; **Ea = 109 kJ/mol (125-155 °C)** or 98 (135-155).
   Upper bound on the true decay rate's error: A243 is not specific (Flag 2), so true ASA loss is likely
   faster.
4. **Browning k is not reproducible from Fig. 1A** (Flag 1): exponential-growth slopes of ln A420 over
   the early window are ≤ 0.14 min⁻¹ at 155 °C, and a first-order approach to the 18.2 plateau gives
   ~0.055 min⁻¹ at 155 °C; Table 1 prints 0.539. The **Ea (114 kJ/mol) is internally consistent with
   Table 1** but the k's definition is unknown, so only the Ea (and only as "A420 apparent") is usable.
5. Molar picture at 145 °C: Cys gone (> 99 %) by 20 min while only ~6 % of ASA is lost; ASA keeps
   browning for 2 h after. So **≥ 0.94 mol of ASA per mol Cys survives the Cys-consuming phase**; the
   Cys:ASA consumption ratio in the first 20 min is ≥ 10 : 1 (if A243 reads true). That points to Cys
   reacting with a minor ASA-derived species (DHA or its products) faster than ASA is degraded, not 1:1
   with ASA.

## 6. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| Cys total loss, ASA pot | first-order k | 0.039 / 0.115 / 0.162 / ~0.25 | min⁻¹ | 0.1 M Cys + 0.1 M ASA, 0.2 M phosphate pH 8 (initial), 125 / 135 / 140 / 145 °C | Fig. 2B [D]; arithmetic 1 | derived_from_digitised |
| Cys total loss, ASA pot | Ea | ≈ 110-125 | kJ/mol | 110-145 °C | arithmetic 1 | derived |
| ASA loss (A243) | first-order k | 1.2e-3 … 1.17e-2 | min⁻¹ | 125-155 °C, same | Fig. 2A [D]; arithmetic 3 | derived_from_digitised, non-specific assay |
| ASA loss | Ea | ≈ 98-109 | kJ/mol | 125-155 °C | arithmetic 3 | derived |
| browning A420 | apparent Ea | **114.33** | kJ/mol | 125-155 °C | Table 1 / text p. 506 | printed (k definition unclear) |
| browning A420 | k | 0.0466-0.539 | (min⁻¹?) | 125-155 °C | Table 1 | printed, not reproducible |
| Cys self-browning | A420 | 0.08 | — | Cys alone, 140 °C, 2 h | text p. 506 | printed |
| H2S, volatiles, pH drift, DHA | — | **NOT REPORTED** | — | — | — | — |

## 7. What is NOT reported

- No H2S, no sulfur volatiles (those are in Yu & Zhang 2010 and Tang 2015), no thiazolidine, no
  dehydroascorbic acid, no 2,3-diketogulonate, no C5 "PT" quantity.
- No pH after heating (a 0.1 M ASA pot adjusted to pH 8 with NaOH will drift; 0.2 M phosphate holds
  some of it).
- No t = 0 values, no heat-up time; no Cys at 110 °C below 3 mg/mL.
- Units of the browning k; the rate-law constants k1, K3, K4, k2 are not fitted.

## 8. Flags

1. **Table 1's k values cannot be reproduced from Fig. 1A** under either reading of "first order"
   (arithmetic 4). Possibly a fit on transformed or normalised data not described.
2. **ASA by A243 reads above nominal** at 110 °C (18.0-18.4 vs 17.6 mg/mL) — the assay picks up products
   absorbing at 243 nm; residual ASA is therefore overstated and its loss rate understated. Treat
   ASA k's as lower bounds.
3. The proposed mechanism ("Cys would firstly undertake thermal self-degradation ... to form DC") is
   contradicted by the paper's own Cys-alone control (barely any browning) and by the engine's measured
   thermolysis being 25× slower than this Cys loss; Cys is reacting with ASA-derived carbonyls.
4. Fig. 1B (A294) y-axis 0-120 after dilution; not digitised (qualitative use only).

## 9. Use in the model

- **No pentose, no furanone, no thiol volatile** — it does not bear on B61 or the MFT/FFT branch
  directly. The engine has no ascorbic-acid lane; using this pot needs ASA → dicarbonyl chemistry the
  engine lacks.
- **One transferable, qualitative point:** in a carbonyl-rich pot at pH 8, **free cysteine disappears with
  an apparent Ea of ~110-125 kJ/mol and k(145 °C) ~0.2 min⁻¹ — 25× the engine's thermolysis** —
  consistent with the flux map's picture that free cysteine is consumed by carbonyl chemistry, not
  thermolysis, whenever a reactive carbonyl is present. The Ea is similar to the engine's
  `k_cys_h2s` (133) and well above `k_cys_thermal`'s fitted 55.1, which is a weak argument that the
  fitted 55 kJ/mol cysteine-sink barrier is too low (different pot, different carbonyl; weak).
- **Browning Ea 114 kJ/mol** (A420, ASA + Cys, pH 8) sits in the corpus's carbohydrate-class range
  (B55/B61 median 103.7) — usable as one more point in `ea_class_prior` only if browning-Ea rows are
  admitted there (they are an apparent, multi-step observable).
