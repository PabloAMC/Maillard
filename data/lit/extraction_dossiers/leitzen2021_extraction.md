# Leitzen et al. 2021 — EXTRACTION (glucose degradation products in autoclaved 10 % glucose)

**Source on disk:** `data/articles/Leitzen2021.pdf` (owner's download, 2026-09-07; open access). Read-only.

| field | value |
|---|---|
| Title | "Quantification of Degradation Products Formed during Heat Sterilization of Glucose Solutions by LC-MS/MS: Impact of Autoclaving Temperature and Duration on Degradation" |
| Venue | Pharmaceuticals 2021, 14, 1121 |
| DOI | 10.3390/ph14111121 |

## 1. Methods

10 % (w/v) glucose in water (555 mmol/L) in PP bottles (volume and headspace not stated), autoclaved at
111 / 116 / 121 C under scheme A (autoclaving times 180 / 57 / 18 min, F0-based "overkill") and scheme B
(233 / 85 / 30 min), plus 121 C for 350 min; LC-MS/MS of the o-PDA quinoxalines (Mittelmaier method,
modified). pH before 4.98, after 4.4-5.2 (scheme A) / 4.1 (scheme B). Non-autoclaved reference in Table 6.

## 2. Table 4 (scheme A) and Table 5 — ug/mL (n = 27; n = 9 for Table 5)

| T (C) / time | GO | MGO | glucosone | 3-DG (+3-DGal) | 3,4-DGE | 5-HMF |
|---|---:|---:|---:|---:|---:|---:|
| 111 / 180 min | 4.4 +/- 2.7 | 3.0 +/- 0.3 | 5.9 +/- 2.2 | 56.0 +/- 10.2 | 59.6 +/- 14.5 | 81.9 +/- 29.5 |
| 116 / 57 min | 4.2 +/- 0.7 | 2.5 +/- 0.2 | 7.1 +/- 0.6 | 55.0 +/- 1.6 | 50.9 +/- 1.7 | 31.6 +/- 0.5 |
| 121 / 18 min | 5.6 +/- 1.3 | 2.6 +/- 0.2 | 7.5 +/- 1.4 | 52.2 +/- 4.0 | 55.5 +/- 1.7 | 17.4 +/- 3.9 |
| 121 / 350 min | 8.3 +/- 0.0 | 1.2 +/- 0.1 | 1.1 +/- 0.1 | 13.7 +/- 0.1 | 12.5 +/- 0.3 | 41.1 +/- 0.1 |
| reference (unheated) | 1.0 +/- 0.5 | 0.9 +/- 0.1 | 0.1 | n.d. | n.d. | 0.1 |

Scheme B (longer): GO 18-23, MGO 12-13, glucosone 5-7, 3-DG 60-73, DGE 60-74, HMF 37-94.

## 3. Reading

Levels are credible (3-DG 56 ppm from 555 mM glucose = 0.06 % conversion; the unit is ug/mL of
solution). The ORDERING at every scheme-A point: 3-DG ~ 3,4-DGE > HMF > glucosone > GO > MGO, with
glucosone / 3-DG = 0.11-0.14. The long 121 C hold (350 min) consumes 3-DG, 3,4-DGE, glucosone and MGO
(down 75-85 %) while GO and HMF rise: the deoxyosones are intermediates, not end products, on that
timescale. Diacetyl not measured.

## 4. What the repo takes

Directional claim DIC-03: glucose alone, 121 C, 18 min, aqueous: 3-DG > glucosone > glyoxal > methylglyoxal
(the four species the trunk carries), independent of the B7/B13 sources. The same check DIC-01 made from
Zhang 2020 — whose absolute levels are unusable — now stands on credible numbers. Recorded also: the
repo's `mp_holdout_glucose_only_autoclave_121C_Steinhagen2021` bundle cites this DOI for HMF.
