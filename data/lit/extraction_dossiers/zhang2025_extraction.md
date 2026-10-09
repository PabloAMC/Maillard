# Zhang H, Cui, Xia, Hussain, Hayat, Zhang X & Ho 2025 — EXTRACTION (Nα,Nε-di(1-deoxy-D-xylulos-1-yl)lysine, 20 mmol/L, pH 7.5 initial, heated at 100/120/140 °C for 30-120 min with and without 20 mmol/L xylose: deoxypentosones, glyoxal, methylglyoxal, released lysine, pH at 120 °C; furfural and other volatiles by HS-SPME)

**Source on disk:** `data/articles/Zhang2025.pdf` (the publisher's PDF, 12 pp.). Pages 1-10 read and
checked by eye on 2026-10-09 from the page images, cross-checked against `pdftotext -layout`; Figure 3
(p. 7) re-read from a 220 dpi crop. Tables 1-3 matched the text layer cell for cell (the layer drops the
methyl-pyrazine row of Table 2; the image has it). The Supplementary Material (Tables S1-S4, Figs. S1-S7)
is **not on disk**. Written for the `k_arp_dpo` / `k_arp_tdp` comparison. Not to be confused with
Zhang et al. 2026 (`zhang2026_extraction.md`), the source of the Amadori barrier the engine uses.

| field | value |
|---|---|
| Title | "Dual role of exogenous xylose in regulating pyrazines and furans formation during the thermal degradation of Nα,Nε-di(1-deoxy-D-xylulos-1-yl)lysine through temperature, reaction time, and xylose concentration control" |
| Authors | Han Zhang, Heping Cui, Xue Xia, Shahzad Hussain, Khizar Hayat, Xiaoming Zhang, Chi-Tang Ho (Jiangnan / King Saud / Alabama A&M / Rutgers) |
| Venue | Food Chemistry 2025, 479, 143828 |
| DOI | 10.1016/j.foodchem.2025.143828 |

## 1. Methods

- **The Amadori compound (p. 2-3, Fig. 1).** **Nα,Nε-di-Xul-Lys ARP**: lysine carrying a
  1-deoxy-D-xylulos-1-yl residue on BOTH the α- and the ε-nitrogen (two xylose-derived residues per
  lysine). Made from xylose : lysine 5 : 1, pH 7.5, 80 °C, 40 min under vacuum, purified on Dowex 50WX8
  H⁺, purity ≥ 95 %. It degrades stepwise (Fig. 1) to the mono-glycated Nα-Xul-Lys and Nε-Xul-Lys ARPs
  plus deoxypentosone, then to lysine plus 3-deoxypentosone (3-DX) and 1-deoxypentosone (1-DX); 3-DX
  retro-aldolises to methylglyoxal (MGO) and glycolaldehyde, which oxidises to glyoxal (GO). Exogenous
  xylose can re-glycate the mono-ARPs back to the di-ARP.
- **Heating (p. 3, §2.3-2.4).** 20 mmol/L di-ARP, alone or + 20 mmol/L xylose, **pH set to 7.5 at the
  start, unbuffered**; 100, 120, 140 °C for 30, 60, 90, 120 min; sealed vials, stirred oil bath
  (120 ± 1 °C stated for the xylose-dose series). Xylose dose series 20-100 mmol/L at 120 °C, 60 and
  120 min; xylose added at 0-120 min into a 120 min run.
- **Analysis.** ARPs, xylose, lysine by HPLC (amide column), method in Zhang 2023. α-Dicarbonyls by
  o-phenylenediamine derivatisation, HPLC-DAD, external standards. Volatiles (furfural included) by
  HS-SPME-GC/MS at 60 °C, DB-WAX, external calibration curves (Table S1, SI). n = 3.

## 2. Findings that matter

**Figure 3 (p. 7), 120 °C, 20 mmol/L di-ARP, pH 7.5 initial. All values read from graph, approx., mmol/L
unless stated.** Black = di-ARP alone; red = di-ARP + 20 mmol/L xylose.

| t (min) | 3-DX alone / +Xyl | 1-DX alone / +Xyl | GO alone / +Xyl | MGO alone / +Xyl | released Lys alone / +Xyl | pH alone / +Xyl | A420 alone / +Xyl |
|---|---|---|---|---|---|---|---|
| 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 7.5 / 7.5 | — |
| 30 | 0.12 / 0.44 | 0.24 / 0.33 | 0.84 / 2.31 | 0.032 / 0.130 | 3.35 / 3.79 | 6.62 / 7.00 | 0.514 / 0.509 |
| 60 | 0.20 / 0.63 | 0.31 / 0.27 | 1.19 / 1.89 | 0.066 / 0.097 | 4.13 / 4.96 | 5.86 / 6.16 | 0.603 / 0.611 |
| 90 | 0.35 / 0.41 | 0.27 / 0.14 | 1.14 / 1.22 | 0.051 / 0.102 | 3.27 / 2.92 | 5.31 / 5.11 | 0.618 / 0.645 |
| 120 | 0.39 / 0.28 | 0.11 / 0.09 | 1.26 / 1.09 | 0.070 / 0.080 | 2.90 / 1.56 | 4.82 / 4.43 | 0.654 / 0.701 |

Text anchors (p. 5-6) agree with the graph: GO 0.84-1.26 and MGO 0.03-0.07 mmol/L (alone); 3-DX 0.63 and
1-DX 0.33 mmol/L maxima with xylose, "3.15 and 1.38 times" the di-ARP-alone values (these ratios are
time-matched, 0.63/0.20 at 60 min and 0.33/0.24 at 30 min, derived here). **Discrepancy:** the text gives
the +Xyl GO peak as 2.49 mmol/L; Figure 3C plots ~2.31.

Derived here (alone, 120 °C): at 60 min the free lysine (4.13 mM) means at least 2 × 4.13 = 8.3 mM of the
40 mM xylulosyl residues have left (≥ 21 %), while 3-DX + 1-DX stand at only 0.51 mM and GO at 1.19 mM,
so the deoxyosones turn over quickly to C2/C3 fragments at this temperature. pH falls 2.7 units in
120 min, so no rate here is at a fixed pH.

**2-Furfural, µg/L, di-ARP + xylose (Tables 1-3, p. 4, 6, 8), mean ± SD:**

| T | 30 min | 60 min | 90 min | 120 min |
|---|---|---|---|---|
| 100 °C | 11.66 ± 1.35 | 26.87 ± 3.02 | 35.61 ± 3.98 | 47.47 ± 2.43 |
| 120 °C | 14.86 ± 1.11 | 36.15 ± 3.24 | 70.27 ± 5.21 | 164.39 ± 12.95 |
| 140 °C | 494.07 ± 13.84 | 327.03 ± 17.69 | 1129.67 ± 34.13 | 1182.30 ± 21.99 |

Di-ARP alone: only the 100 °C furfural peak, **7.51 µg/L**, survives in the main text (p. 4; Table S2 in
the SI). Scale (derived here): 47.47 µg/L = 0.49 µmol/L, about 2.5 × 10⁻⁵ mol per mol di-ARP.

**Totals (Fig. 2, p. 5; text values printed, others read from graph, approx.), µg/L:** furans alone
24.91 (100 °C, 30 min), 36.23 (100 °C peak, 90 min), 23.56 (100 °C, 120 min), 185.42 (120 °C, 120 min),
2415.56 (140 °C, 60 min); furans + Xyl 87.96 -> 179.59 (100 °C), 333.37 (120 °C, 120 min), 1245.63
(140 °C peak, 90 min). Pyrazines alone 4.15 -> 8.27 -> 5.86 (100 °C), 19.21 (120 °C, 120 min), 209.14
(140 °C, 60 min); + Xyl 54.60 (120 °C, 60 min), 400.95 (140 °C, 60 min) -> 40.37 (120 min).

**Internal inconsistencies to know about.** The same nominal condition (1 : 1 xylose, 120 °C, 120 min)
reads 333.37 µg/L furans in Fig. 2C, 323.37 in Fig. 4A (xylose at 0 min) and 303.37 in Fig. 5C (text
p. 10): three runs, not one. Fig. 5's y-axes are labelled "mmol/L" while the text uses µg/L for the same
numbers.

## 3. What it means for the model

**Different compound from the engine's ARP.** The sulfur lane's `ARP` is
"N-(1-deoxy-D-xylulos-1-yl)-alanine" (`src/kinetic_core/species_sulfur.py` line 91), a mono-glycated
alanine Amadori fed by Zhou 2023. This paper's precursor is a **di-glycated lysine**, whose first step is
loss of one residue to a mono-ARP and whose free amine is regenerated only after both leave. The engine's
`DPO` = 1-deoxypentosone (1-DX here), `TDP` = 3-deoxypentosone (3-DX here).

Live values from `results/validation/core_prediction_uncertainty.json` (first order, 1/min, at 145 °C;
`k_arp_dpo` carries the base pH factor, `k_arp_tdp` the acid one):

| key | centre | distribution / band | reason |
|---|---|---|---|
| `b8.k_arp_dpo.log10_k_ref_145C` | -2.310 (4.9e-3 min⁻¹, derived) | normal_log10, σ 1.805, band [-10, 0.5] | laplace_covariance_at_b8_optimum; **bound_limited** in data_wishlist §1 |
| `b8.k_arp_tdp.log10_k_ref_145C` | -1.727 (1.9e-2 min⁻¹, derived) | normal_log10, σ 0.607, band [-10, 0.5] | laplace_covariance_at_b8_optimum |
| `b8.k_arp_dpo_th.log10_k_ref_145C` | -1.290 | fixed | frozen in the sulfur fit |
| `b8.k_arp_tdp_th.log10_k_ref_145C` | -4.354 | fixed | frozen in the sulfur fit |
| `b8.k_tdp_fur.log10_k_ref_145C` (3-DX -> furfural) | -3.037 | normal_log10, σ 0.340 | laplace_covariance_at_b8_optimum |
| `b8.k_osone_decay.log10_k_ref_145C` | -1.187 | normal_log10, σ 0.449 | laplace_covariance_at_b8_optimum |

The barrier on `k_arp_dpo` and `k_arp_tdp` is **85.7 kJ/mol, a measured override the fit cannot move**
(`ZHANG_EA_CYS_AMADORI_TO_ALPHA_DC_KJ_MOL`, Zhang 2026 k16 refit; `kinetic_core_b9_fit_report.json`
`t_structure.measured_barriers_the_fit_cannot_move`).

**Can it pin `k_arp_dpo` / `k_arp_tdp`?** No. The main text prints **no ARP loss time course at any
temperature, no rate constant and no activation energy**. Deoxyosone time courses exist only at 120 °C,
four points, with pH drifting from 7.5 to 4.8; at 100 and 140 °C they are in SI figures S5/S6 (not on
disk). Furfural is reported as µg/L by HS-SPME, five orders below the precursor, with the di-ARP-alone
series in SI Tables S2-S4. The compound is a di-glycated lysine, not the mono-glycated alanine ARP the
engine carries. What it does give is qualitative: at 120 °C the deoxyosone pools stay below 0.4 mM while
GO reaches ~1.2 mM within 60 min, so 3-DX/1-DX consumption (including retro-aldol to C2/C3, the engine's
`k_dpo_c2c3` / `k_osone_decay` territory) is fast against their formation from this ARP. And furfural
rises steeply between 120 and 140 °C (164 -> ~1180 µg/L at 120 min with xylose), which is a shape check on
the TDP -> FUR leg, not a number for it.

## What it does not give

- The loss of the di-ARP (or the mono-ARPs) against time: not in the main text; the authors' earlier
  paper (Zhang et al. 2024b, cited as showing xylose accelerating di-ARP degradation) and the SI hold it.
- Any rate constant, Ea, or Arrhenius treatment.
- 3-DX, 1-DX, GO, MGO, lysine or pH at 100 or 140 °C (SI only), or furfural for the di-ARP alone beyond
  the 7.51 µg/L peak at 100 °C.
- Molar furfural yields with a stated calibration (Table S1 is in the SI); buffered pH.
