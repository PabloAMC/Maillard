# Larralde, Robertson & Miller 1995 (PNAS) — EXTRACTION (free ribose, no amine, 0.05 M deuterated phosphate or acetate, pD 4.4-8.6, 60-118 °C; first-order loss of the furanose 1-H NMR signal; a pD-rate profile at four temperatures, an Arrhenius line at pD 7.4, a phosphate general-base term at pD 7.4 / 100 °C, and twelve other aldoses at one condition)

### THE ONLY MEASURED pH LAW FOR AMINE-FREE PENTOSE LOSS IN THE CORPUS, AND THE ONLY BUFFER-CATALYSIS TERM. Between pD 4 and 6 free ribose disappears ~1 decade faster per pD unit (base-, not acid-catalysed), levelling off near pD 7; at pD 7.4 / 100 °C the phosphate term is already 1.6x the buffer-free rate at 0.05 M and would be ~16x it at Hofmann's 0.5 M. The model's two amine-free pentose routes are both pH-independent.

**Source on disk:** `data/articles/larralde1995.pdf` (3 pp., scanned journal pages with an OCR text
layer). Prose read from `pdftotext -layout`; **Figs. 1, 2 and 3 digitised from 300 dpi renders of
p. 8159** (scratchpad `pent/lar_fig1.png`, `lar_fig2.png`, `lar_fig3.png`). Digitised values are
tagged **[D]** with their read error. No supplementary material exists.

## 0. Identity

| field | value |
|---|---|
| Title | "Rates of decomposition of ribose and other sugars: Implications for chemical evolution" ✔ |
| Authors | Rosa Larralde, Michael P. Robertson, **Stanley L. Miller** (contributed by S. L. Miller, 30 May 1995) — Dept. of Chemistry and Biochemistry, University of California at San Diego |
| Venue | Proc. Natl. Acad. Sci. USA **92**, 8158-8160, August 1995 (Chemistry) ✔ |
| DOI as printed | **not printed** (1995 PNAS page). Volume, issue month and pages match 10.1073/pnas.92.18.8158 |
| Paper type | Origin-of-life argument; the kinetics are a means. **No product was identified** ("No attempt was made to identify the decomposition products"). |
| Verdict | **Identity verified** (title, authors, journal, volume, pages). |

## 1. Why it matters

The model's amine-free pentose lane (`sulfur.py`) has two first-order steps: `r_pent_caramel`
(PENT → TDP, the 1,2-enolisation; operative pH factor `""` per `REACTION_PH_FACTOR`) and
`r_pent_thermal_dpo` (PENT → DPO; `""`). B60 values at 145 °C: **k_pent_caramel = 10^-2.249 =
5.6 × 10⁻³ min⁻¹, k_pent_thermal = 3.7 × 10⁻⁴ min⁻¹** (`kinetic_core_b60_fit_report.json`,
`frozen_parameters.log10_k_ref_at_145C`; units mmol/L and min per `parameters_sulfur.py` l. 52), both
on the lumped formation Ea of 64.1 kJ/mol, and **no buffer term anywhere**. This paper measures the
same quantity (amine-free aldopentose loss) and finds:

1. **A strong pH law that rises with pD** (slope ≈ +1 per pD unit from pD 4 to 6, plateau above
   pD ~7.5), fitted by the authors as k = k₀/(1 + [H⁺]/K_a). The rate is "approximately
   proportional to the OD⁻ ion concentration in the pD range 4-6". Acid catalysis appears only
   "at pD < 4 … (data not shown)".
2. **General-base catalysis by phosphate** at pD 7.4 / 100 °C (Section 3, rate law).
3. **Ea = 107 kJ/mol** at pD 7.4 (25.7 kcal/mol), against the model's lumped 64 kJ/mol.

It does not say which deoxyosone is made: it measures sugar disappearance only.

## 2. Methods as they matter to a model

- **Medium.** Sugar "dissolved in aqueous 0.05 M deuterated phosphate [pD 6.4-8.4] or acetate
  (pD 4.4-5.4) buffer and sealed in NMR tubes". So the **low-pD points are in acetate, the high-pD
  points in phosphate** — buffer identity changes across the profile (Flag 3).
- **Sugar concentration: NOT reported.**
- **pD.** Measured at 25 °C, pD = pH + 0.4; "The pD was corrected for the temperature by ΔpD =
  ΔpK" (Fig. 2 caption: points "adjusted for the change in the buffer pD with temperature"). Whether
  Fig. 1's abscissa is the 25 °C or the corrected pD is not stated explicitly (Flag 4).
- **Observable.** 300 MHz ¹H NMR, decrease of the 1-H signal of the α- and β-furanose forms against
  an external CHCl₃/CCl₄ capillary. Other sugars: the 1-H signals furthest from water.
- **Isotope effect.** H₂O pH 7.0 vs D₂O pD 7.4 at 100 °C: k_H2O/k_D2O = **1.0 ± 0.2**, so the D₂O
  numbers transfer to water at pH ≈ pD − 0.4.
- **Buffer choice.** "the lowest buffer concentration that maintained the pH at a constant value was
  used (0.05 M phosphate or acetate)" — decomposition products are acidic.
- **Temperatures.** Abstract says 40-120 °C; figures show **60, 80, 100, 118 °C** only. No 40 °C
  datum is printed.
- Replicates, error bars, standard deviations: **not reported**.

## 3. Printed numbers (verbatim)

| quantity | value | conditions | locator |
|---|---|---|---|
| rate law | **rate = 2.86 × 10⁻³ + 9.36 × 10⁻² ([DPO₄²⁻] + [D₂PO₄⁻])**, printed unit "(s⁻¹)" | pD 7.4, 100 °C, ribose | p. 8158 col. 2 |
| Arrhenius | **log k (s⁻¹) = 11.23 − 5608/T**; "heat of activation of 25.7 kcal·mol⁻¹" | pD 7.4, 0.05 M phosphate, 60-118 °C | p. 8159 col. 2, Fig. 2 |
| ribose t½ | **73 min** | pD 7.4 (pH 7.0), 100 °C | abstract; p. 8159 |
| ribose t½ (extrapolated) | 300 days at 25 °C; 44 years at 0 °C | pD 7.4 | p. 8159 |
| bicarbonate vs phosphate | rate in 0.05 M bicarbonate = **0.52 ×** that in 0.05 M phosphate | pD 7.4, (100 °C implied) | p. 8159 col. 1 |
| 2-deoxyribose | t½ **225 min** (2.6× slower than ribose) | pD 7.4, 100 °C | p. 8159 |
| ribose 5-phosphate | t½ **7 min** (12× faster than ribose) | pD 7.4, 100 °C | p. 8159 |
| ribose 2,4-bisphosphate | t½ **31 min** | pD 7.4, 100 °C | p. 8159 |
| pH 8.2 vs 7.0 | "would increase the rate by about 50%" | — | p. 8159-8160 |
| Fig. 1 fit form | k = k₀/(1 + [H⁺]/K_a), "K_a and k₀ are adjusted constants" — **values not printed** | — | Fig. 1 caption |

**Unit check on the rate law (mine).** At 0.05 M phosphate, read in s⁻¹ the law gives
k = 7.54 × 10⁻³ s⁻¹, t½ = 1.5 min — 48x faster than the paper's own 73 min. Read in **min⁻¹** it gives
7.54 × 10⁻³ min⁻¹, **t½ = 92 min**, within 26 % of 73 min and of Fig. 2's 100 °C point (84 min). The
printed "s⁻¹" is almost certainly a typo for **min⁻¹** (same conclusion as
`docs/guides/LITERATURE_SEARCH_2026-10_PATHS.md`). Carry: k₀ = 2.86 × 10⁻³ min⁻¹ (buffer-free
intercept), k_phos = 9.36 × 10⁻² L mol⁻¹ min⁻¹ (total phosphate), pD 7.4, 100 °C. Evidence class
**printed_with_unit_correction**.

**Arrhenius check (mine).** 5608 × ln10 × R = **107.4 kJ/mol** (= 25.7 kcal ✔); A = 10^11.23 =
1.7 × 10¹¹ s⁻¹. The equation gives t½ = 72.7 min at 100 °C ✔, 516 min at 80 °C, 12.5 min at 120 °C,
**1.75 min at 145 °C** (extrapolated). Note the Arrhenius line is **in 0.05 M phosphate**, so it is
buffer-term plus intercept, not the intercept alone.

## 4. Digitised figures [D]

### Fig. 1 — log k (s⁻¹) vs pD, ribose, four temperatures (read error ± 0.1 pD, ± 0.05 log k)

| T | (pD, log k) points |
|---|---|
| 118 °C | (4.71, −5.20), (5.71, −4.52), (6.63, −3.50), (7.64, −2.97), (8.64, −2.72) |
| 100 °C | (4.61, −6.17), (5.62, −5.09), (6.56, −4.41), (7.57, −3.94), (8.58, −3.74) |
| 80 °C | (4.52, −7.01), (5.53, −6.17), (6.46, −5.13), (7.46, −4.72), (8.47, −4.64) |
| 60 °C | (5.46, −6.98), (6.21, −6.28), (7.40, −5.54), (8.41, −5.38) |

Fitted plateaus (curve, [D]): **118 °C −2.95, 100 °C −3.80, 80 °C −4.60, 60 °C −5.33**. The slope-1
asymptotes meet the plateaus at an apparent **pK_a ≈ 7.0 (pD)** at both 100 and 118 °C (mine, from the
low-pD points and the plateaus; ± 0.2). Low-pD asymptotes: log k ≈ pD − 10.78 (100 °C), pD − 9.91
(118 °C), pD − 11.53 (80 °C).

### Fig. 2 — Arrhenius, pD 7.4, 0.05 M phosphate (± 0.05 log k)

118 °C −3.03; 100 °C −3.87; 80 °C −4.71; 60 °C −5.54 (log k, s⁻¹).

### Fig. 3 — other aldoses, pD 7.4, 100 °C, 0.05 M phosphate (± 0.04 log k)

| sugar | log k (s⁻¹) | k (min⁻¹) | t½ |
|---|---|---|---|
| ribose | −3.94 | 6.9 × 10⁻³ | 101 min |
| lyxose | −4.17 | 4.0 × 10⁻³ | 2.9 h |
| **xylose** | **−4.60** | **1.5 × 10⁻³** | **7.7 h** |
| arabinose | −4.60 | 1.5 × 10⁻³ | 7.7 h |
| glucose (extrapolated from 80 °C, 1 M phosphate) | −5.28 | 3.1 × 10⁻⁴ | — |

Note ribose's Fig. 3 point (−3.94) sits 0.14 log below the text's 73 min (−3.80) (Flag 5).
Pentose order at pD 7.4 (relative to ribose, mine): **ribose 1 : lyxose 0.59 : xylose 0.22 :
arabinose 0.22**. Compare Garrett 1969's specific-acid order 1 : 0.41 : 0.41 : 0.19 (lyxose, xylose,
arabinose): xylose is relatively slower under base than under acid. The authors attribute the order to the fraction
of open-chain aldehyde (line drawn through ribose); xylose and arabinose fall **below** the line.

## 5. Arithmetic (all mine)

**A. Phosphate acceleration at pD 7.4 / 100 °C** (min⁻¹ reading). Buffer-free 2.86 × 10⁻³ min⁻¹;
0.05 M: ×2.6; **0.5 M (Hofmann's phosphate): ×17.4**. The law sums D₂PO₄⁻ and DPO₄²⁻ at one pD,
so it **cannot apportion activity between H₂PO₄⁻ (the species at pH 3) and HPO₄²⁻**.

**B. Extrapolated amine-free ribose loss at 145 °C** (Fig. 1 100 °C curve, k₀/(1+[H⁺]/K_a) with
plateau −3.80 and pK_a 6.98; pD = pH + 0.4; scaled 100 → 145 °C with the pD-7.4 Ea, ×41.4; in 0.05 M
buffer):

| pH | k (100 °C, s⁻¹) | k (145 °C, min⁻¹) | model k_pent_caramel + k_pent_thermal (145 °C, min⁻¹) | model / Larralde |
|---|---|---|---|---|
| 3 | 4.2 × 10⁻⁸ | 1.0 × 10⁻⁴ | 6.0 × 10⁻³ | **~60x too fast** |
| 5 | 4.1 × 10⁻⁶ | 1.0 × 10⁻² | 6.0 × 10⁻³ | 0.6x (agrees) |
| 7 | 1.1 × 10⁻⁴ | 0.29 | 6.0 × 10⁻³ | **~50x too slow** |

Assumptions: the pD-7.4 Ea holds at every pD; acid catalysis below pD 4 (not shown) is ignored
(Garrett 1969 sizes it — `garrett1969_extraction.md` puts the specific-acid term at ~45 % of the
total at pH 3, 145 °C, ribose); ribose ≠ xylose (xylose ~0.22× at pD 7.4). Order-of-magnitude.
Evidence class **derived_extrapolated**.

## 6. Kinetic numbers the repository can use

| step / observable | quantity | value | unit | conditions | source | evidence class |
|---|---|---|---|---|---|---|
| ribose loss, buffer-free intercept | k₀ | 2.86 × 10⁻³ | min⁻¹ (printed s⁻¹) | pD 7.4, 100 °C | rate law p. 8158 | printed_with_unit_correction |
| ribose loss, phosphate general base | k_phos | 9.36 × 10⁻² | L mol⁻¹ min⁻¹ (total P_i) | pD 7.4, 100 °C | same | printed_with_unit_correction |
| ribose loss, Arrhenius | Ea, log A | **107.4 kJ/mol**, 11.23 (s⁻¹) | — | pD 7.4, 0.05 M phosphate, 60-118 °C | p. 8159 | printed (single pD) |
| ribose loss pH law | slope d log k / d pD | **≈ +1** (pD 4-6), plateau above pD ~7.5, pK_a,app ≈ 7.0 | — | 60-118 °C, 0.05 M acetate/phosphate | Fig. 1 [D] | digitised |
| ribose loss vs pD | 19 (pD, log k) points | Section 4 | s⁻¹ | four T | Fig. 1 [D] | digitised time-free rate constants |
| xylose loss | k | 1.5 × 10⁻³ (t½ 7.7 h) | min⁻¹ | pD 7.4, 100 °C, 0.05 M phosphate | Fig. 3 [D] | digitised |
| ribose 5-phosphate loss | t½ | 7 min | — | pD 7.4, 100 °C | text | printed |
| bicarbonate / phosphate | ratio | 0.52 | — | pD 7.4, 0.05 M | text | within_study_ratio |

## 7. What this can and cannot anchor

**Can:** a declared band (or hold-out) on the **total amine-free pentose loss** vs pH and T; a
**phosphate general-base term** at pD 7.4; an **Ea ≈ 107 kJ/mol** for free-sugar loss (one pD).

**Cannot:** the 1,2- vs 2,3-enolisation split (no products); the acid branch below pD 4 (data not
shown); phosphate activity at pH 3-5 (term measured only at pD 7.4); anything amine-catalysed.

## 8. Flags

1. **Rate-law unit almost certainly misprinted** (s⁻¹ for min⁻¹); see Section 3.
2. **Sugar concentration not reported**; first-order assumed, never tested in the paper.
3. **Buffer changes across the pD profile** (acetate ≤ pD 5.4, phosphate ≥ pD 6.4); part of the
   slope between pD 5.4 and 6.4 could be a buffer effect. Acetate catalysis is not measured.
4. **pD abscissa**: Fig. 2 points are temperature-corrected; Fig. 1's are not stated to be.
5. **Internal inconsistency of ~0.1-0.15 log**: text t½ 73 min (log k −3.80), Fig. 2 point −3.87,
   Fig. 3 point −3.94, rate law at 0.05 M (min⁻¹) −3.90.
6. **The plateau may be HPO₄²⁻ catalysis, not an ionisation** — the authors offer both readings.
7. Abstract's 40 °C data are not shown.
8. No product identification: "decomposition" may include products other than either deoxyosone.

## Use in the model

- **Gap (c), pentose enolisation pH law at low pH; phosphate catalysis — directly.** The amine-free
  pentose loss rises ~10x per pH unit from pH ~3.6 to ~5.6 (pD 4-6) and is NOT acid-catalysed above
  pD 4. The model's amine-free routes (`r_pent_caramel`, `r_pent_thermal_dpo`) are both flat in pH,
  which on the extrapolation of Section 5B makes them **~60x too fast at pH 3 and ~50x too slow at
  pH 7** at 145 °C while agreeing at pH 5 (where the fit rows sit). This **contradicts** the model's
  assumption that the amine-free route is a pH-independent "caramelisation" whose dominant product
  is the 3-deoxyosone (`k_pent_caramel` is 15x `k_pent_thermal`): Garrett 1969's specific-acid
  furfural route is ~10⁴x too slow at pH 5 to be that flux, so whatever Larralde measures at pH 5-7
  is a base-catalysed route of unidentified product (general Maillard chemistry would put it on
  the 2,3-enolisation/retro-aldol side, but this paper does not show that).
- **Phosphate.** A **measured general-base term**: at 0.5 M phosphate (every Hofmann arm) and pD 7.4
  it would accelerate free-sugar loss ~17x. It supports the flux map's suspicion that phosphate
  moves the Hofmann pots, but it **cannot support the specific hypothesis that phosphate rescues
  the MFT branch at pH 3**: there the phosphate is H₂PO₄⁻, whose activity this paper does not
  separate, and a general-base term would favour the base-type route whose products are unknown.
- **Temperature.** Ea 107 kJ/mol for free-sugar loss (one pD) is well above the lumped 64 kJ/mol
  and inside the flux map's "100-150 kJ/mol typical"; usable as a declared band on an un-lumped
  trunk Ea in B59.
- **Proposed use:** a pre-registered hold-out row set — Fig. 1's 19 points (ribose, no amine) against
  the engine's PENT decay in an amine-free, buffered pot — only after the amine-free routes acquire
  a pH law; scored as **total loss**, not products. Caveats: ribose not xylose (xylose ~0.22× at
  pD 7.4); D₂O (isotope effect 1.0 ± 0.2); 0.05 M buffer; sugar concentration unknown; ≤ 118 °C.
