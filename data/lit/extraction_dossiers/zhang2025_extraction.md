# Zhang H, Cui, Xia, Hussain, Hayat, Zhang X & Ho 2025 — EXTRACTION (Nα,Nε-di(1-deoxy-D-xylulos-1-yl)lysine, 20 mmol/L, pH 7.5 initial, heated at 100/120/140 °C for 30-120 min with and without 20 mmol/L xylose: deoxypentosones, glyoxal, methylglyoxal, released lysine, pH at 120 °C; furfural and other volatiles by HS-SPME)

**Source on disk:** `data/articles/Zhang2025.pdf` (the publisher's PDF, 12 pp.). Pages 1-10 read and
checked by eye on 2026-10-09 from the page images, cross-checked against `pdftotext -layout`; Figure 3
(p. 7) re-read from a 220 dpi crop. Tables 1-3 matched the text layer cell for cell (the layer drops the
methyl-pyrazine row of Table 2; the image has it). The Supplementary Material (Tables S1-S4, Figs. S1-S8)
is now on disk as `data/articles/Zhang2025_supplementary.docx` (saved first as `Wang2025_supplementary.docx` by mistake and renamed the same day; its title page carries this title and these authors). Read in
full on 2026-10-09; see section 2x. Written for the `k_arp_dpo` / `k_arp_tdp` comparison. Not to be confused with
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

Di-ARP alone: only the 100 °C furfural peak, **7.51 µg/L**, survives in the main text (p. 4). The full
di-ARP-alone series at 100/120/140 °C is in SI Tables S2-S4 (section 2x.3). Scale (derived here): 47.47 µg/L = 0.49 µmol/L, about 2.5 × 10⁻⁵ mol per mol di-ARP.

**Totals (Fig. 2, p. 5; text values printed, others read from graph, approx.), µg/L:** furans alone
24.91 (100 °C, 30 min), 36.23 (100 °C peak, 90 min), 23.56 (100 °C, 120 min), 185.42 (120 °C, 120 min),
2415.56 (140 °C, 60 min); furans + Xyl 87.96 -> 179.59 (100 °C), 333.37 (120 °C, 120 min), 1245.63
(140 °C peak, 90 min). Pyrazines alone 4.15 -> 8.27 -> 5.86 (100 °C), 19.21 (120 °C, 120 min), 209.14
(140 °C, 60 min); + Xyl 54.60 (120 °C, 60 min), 400.95 (140 °C, 60 min) -> 40.37 (120 min).

**Internal inconsistencies to know about.** The same nominal condition (1 : 1 xylose, 120 °C, 120 min)
reads 333.37 µg/L furans in Fig. 2C, 323.37 in Fig. 4A (xylose at 0 min) and 303.37 in Fig. 5C (text
p. 10): three runs, not one. Fig. 5's y-axes are labelled "mmol/L" while the text uses µg/L for the same
numbers. The SI adds two more: Figure S7's 1:1 bars disagree with Figure 3's "+Xyl" series at 60 min
(3-DX 0.24 vs ~0.63), and several Figure S5/S6 error bars are copies (section 2x.1-2x.2). The furan totals
include 2,3-butanedione and 2,3-pentanedione (section 2x.3).

## 2x. Supplementary information (read 2026-10-09)

**Source:** `data/articles/Zhang2025_supplementary.docx` (renamed from `Wang2025_supplementary.docx`; its title page is
this paper: Zhang H., Cui, Xia, Hussain, Hayat, Zhang X. & Ho, "Dual Role of Exogenous Xylose ...
Nα,Nε-di(1-deoxy-d-xylulos-1-yl)lysine ..."). Read in full on 2026-10-09. Contents: Tables S1-S4 (text),
Figures S1-S4 (structures and 1H/13C NMR of the three ARPs, no kinetics), **Figures S5-S7 (α-dicarbonyl and
lysine data)** and Figure S8 (TIC and MS/MS of the purified ARPs). Figures S5-S7 are embedded as vector
EMF drawings (Origin). Their values were read from the drawings' own marker coordinates, calibrated
against the printed axis ticks, and checked by eye against a rendering; they are still **read from graph,
approx.** The plotted values land on round numbers (two decimals, three for MGO; e.g. 0.07997 → 0.08), which
suggests the authors plotted rounded means.

**What the SI does not contain:** no time course of the di-ARP (or of the mono-ARPs Nα-/Nε-Xul-Lys) at any
temperature, no rate constant, no activation energy, no Arrhenius treatment, no pH at 100 or 140 °C, and
no furfural or other volatile for the di-ARP + xylose system beyond what the main text prints.

### 2x.1 Figures S5 and S6: α-dicarbonyls and released lysine at 100 and 140 °C

Caption (S5; S6 identical at 140 °C): "Concentration of 3-DX (a), 1-DX (b), GO (c), MGO (d), and released
Lys (e) in Nα,Nε-di-Xul-Lys ARP and Nα,Nε-di-Xul-Lys ARP/Xyl model system at 100 °C. (concentrations of
Nα,Nε-di-Xul-Lys ARP and exogenous xylose were 20 mmol/L, pH 7.5)". Axes "Concentration (mmol/L)" vs
"Reaction time (min)". All values mmol/L, **read from graph, approx.**; all series start at 0 at t = 0.

**Figure S5, 100 °C** (alone / + 20 mmol/L xylose):

| t (min) | 3-DX alone / +Xyl | 1-DX alone / +Xyl | GO alone / +Xyl | MGO alone / +Xyl | released Lys alone / +Xyl |
|---|---|---|---|---|---|
| 30 | 0.08 / 0.51 | 0.20 / 0.29 | 1.01 / 1.21 | 0.045 / 0.066 | 1.52 / 1.09 |
| 60 | 0.18 / 0.39 | 0.26 / 0.30 | 1.08 / 1.57 | 0.040 / 0.050 | 3.84 / 3.06 |
| 90 | 0.29 / 0.46 | 0.29 / 0.25 | 1.52 / 1.63 | 0.110 / 0.070 | 4.27 / 4.92 |
| 120 | 0.30 / 0.38 | 0.20 / 0.21 | 1.28 / 1.84 | 0.061 / 0.086 | 3.71 / 2.86 |

**Figure S6, 140 °C** (alone / + 20 mmol/L xylose):

| t (min) | 3-DX alone / +Xyl | 1-DX alone / +Xyl | GO alone / +Xyl | MGO alone / +Xyl | released Lys alone / +Xyl |
|---|---|---|---|---|---|
| 30 | 0.45 / 0.99 | 0.09 / 0.27 | 1.31 / 1.60 | 0.080 / 0.058 | 3.96 / 4.29 |
| 60 | 0.77 / 0.23 | 0.03 / 0.11 | 0.93 / 0.48 | 0.086 / 0.065 | 3.44 / 2.96 |
| 90 | 0.34 / 0.21 | 0.11 / 0.08 | 0.88 / 0.29 | 0.084 / 0.088 | 2.57 / 1.42 |
| 120 | 0.20 / 0.13 | 0.07 / 0.02 | 0.58 / 0.12 | 0.073 / 0.052 | 1.91 / 1.06 |

**Error bars are not independent (derived here from the vector data).** The error-bar half-lengths, in
drawing units, are identical to within 1 unit (i) between the two series in the 1-DX and MGO panels, at
both temperatures, and (ii) between Figure S5 and Figure S6 for the 1-DX, MGO and released-Lys panels. Two
series at two temperatures cannot have the same four SDs by chance, so the bars were evidently copied
between panels. Only the means are transcribed here.

Reading of the shapes (with Figure 3 at 120 °C from the main text, section 2):
- 100 °C: 3-DX (alone) climbs slowly to 0.30 at 120 min, and 1-DX peaks at 0.29 (90 min). Free lysine is
  still rising to 4.27 at 90 min.
- 140 °C: everything peaks by 30-60 min and then falls. Free lysine falls from 3.96 (30 min) to 1.91
  (120 min), so lysine is consumed after release. 1-DX stays ≤ 0.11 alone.
- With xylose, 3-DX and 1-DX are highest at the first point (30 min) and lowest late at 140 °C, matching
  the main text's "dual role".

### 2x.2 Figure S7: xylose dose series at 120 °C

Caption: "Concentration of 3-DX, 1-DX, GO, MGO in Nα,Nε-di-Xul-Lys ARP model system after adding different
molar concentrations of exogenous Xyl and reacting at 120 °C for 60 (a ~ d) or 120 min (e ~ h)". Bar charts;
x axis "Molar ratio of Nα,Nε-di-Xul-Lys ARP to exogenous xylose" 1:1 to 1:5 (20 mmol/L di-ARP, so xylose 20
to 100 mmol/L, derived here). mmol/L, **read from graph, approx.**; letters as printed above the bars.

| di-ARP : Xyl | 3-DX 60 min | 1-DX 60 min | GO 60 min | MGO 60 min | 3-DX 120 min | 1-DX 120 min | GO 120 min | MGO 120 min |
|---|---|---|---|---|---|---|---|---|
| 1:1 | 0.24 a | 0.10 a | 1.71 a | 0.07 b | 0.28 a | 0.11 a | 1.64 a | 0.07 a |
| 1:2 | 0.31 b | 0.10 a | 2.18 b | 0.09 c | 0.35 b | 0.15 b | 2.65 e | 0.09 b |
| 1:3 | 0.42 c | 0.19 b | 2.98 d | 0.11 a | 0.51 c | 0.23 c | 2.36 d | 0.19 c |
| 1:4 | 0.48 d | 0.22 b | 2.83 e | 0.07 a | 0.54 c | 0.24 c | 2.17 c | 0.21 d |
| 1:5 | 0.56 e | 0.31 c | 2.38 c | 0.08 d | 0.76 d | 0.22 c | 1.39 b | 0.11 c |

Some letters contradict the bar heights (MGO 60 min: 0.11 is "a", 0.07 is "b"; GO 60 min: 2.83 is "e",
2.98 is "d"; MGO 120 min: 0.11 is "c" like 0.19), as printed.

**Inconsistency with main-text Figure 3.** The 1:1 bars are the same nominal condition as Figure 3's
"+Xyl" series at 120 °C, but do not match at 60 min: 3-DX 0.24 (S7) vs ~0.63 (Fig. 3, and the text's
maximum "0.63"), 1-DX 0.10 vs ~0.27, GO 1.71 vs ~1.89, MGO 0.07 vs ~0.097. At 120 min 3-DX agrees (0.28 vs
~0.28); GO does not (1.64 vs ~1.09). Like the three different furan totals in section 2, this points to
separate runs, or to figures that do not come from one dataset.

### 2x.3 Tables S2-S4: volatiles from the di-ARP alone (no xylose)

Captions: "Main volatile compounds produced from Nα,Nε-di-Xul-Lys ARP degradation under [100 / 120 / 140]
°C for 30, 60, 90 and 120 min (initial pH was 7.5)". Units "Concentration (µg/L)". Footnotes as printed:
a, mean ± SD (n = 3), letters = Duncan's test p < 0.05 within a row; b, linear RI against C7-C30 n-alkanes
on a 30 m × 0.25 mm × 0.25 µm DB-WAX column; c, MS = NIST library, RI = retention index, **S = authentic
standard**; d, ND = not detected. Every row is coded "MS, RI, S". Transcribed in full.

**Table S2, 100 °C:**

| Compounds | RI^b | Identification methods^c | 30 min | 60 min | 90 min | 120 min |
|---|---|---|---|---|---|---|
| furan | 797 | MS, RI, S | 0.09±0.00a | 0.80±0.11b | 1.45±0.76d | 1.25±0.34c |
| 2-methyl-furan | 851 | MS, RI, S | 5.25±0.52a | 8.07±1.01c | 11.47±1.28d | 6.48±0.67b |
| 2-furfural | 1457 | MS, RI, S | 3.19±0.88b | 7.51±1.12c | 3.58±0.95b | 1.28±0.29a |
| 2,3-butanedione | 979 | MS, RI, S | 10.55±1.46a | 14.53±2.01b | 14.66±2.24b | 10.28±1.92a |
| 2,3-pentanedione | 1065 | MS, RI, S | 0.42±0.01a | 0.74±0.02c | 0.95±0.12d | 0.61±0.23b |
| 2-vinyl-furan | 1096 | MS, RI, S | 5.41±0.61c | 2.20±0.35b | 1.65±0.26a | 2.08±0.31b |
| 3-acetyl-2,5-dimethyl-furan | 1093 | MS, RI, S | ND^d | 1.17±0.15a | 2.48±0.31c | 1.59±0.18b |
| 2-furanmethanol | 1603 | MS, RI, S | ND | ND | ND | ND |
| 2(5H)-furanone | 1748 | MS, RI, S | ND | 0.46±0.13a | ND | ND |
| pyrazine | 1215 | MS, RI, S | 1.88±0.57a | 3.72±0.62c | 3.99±0.41c | 2.43±0.22b |
| methyl-pyrazine | 1265 | MS, RI, S | 2.27±0.61a | 2.65±0.35b | 3.83±0.88d | 3.08±0.57c |
| 2,5-dimethyl-pyrazine | 1305 | MS, RI, S | ND | ND | 0.45±0.16b | 0.35±0.28a |

**Table S3, 120 °C:**

| Compounds | RI^b | Identification methods^c | 30 min | 60 min | 90 min | 120 min |
|---|---|---|---|---|---|---|
| furan | 811 | MS, RI, S | 0.03±0.00a | 0.52±0.02b | 0.86±0.05c | 3.40±0.47d |
| 2-methyl-furan | 871 | MS, RI, S | 6.03±1.35a | 7.35±0.77b | 10.65±2.38c | 25.73±3.72d |
| 2-furfural | 1457 | MS, RI, S | 7.59±1.22a | 18.30±1.95b | 25.49±3.03c | 61.49±3.88d |
| 2,3-butanedione | 989 | MS, RI, S | 37.80±1.62a | 56.02±3.01b | 67.30±4.44c | 71.68±3.71c |
| 2,3-pentanedione | 1073 | MS, RI, S | 1.41±0.03a | 1.42±0.05a | 1.49±0.07b | 2.58±0.11c |
| 2-vinyl-furan | 1096 | MS, RI, S | 0.72±0.03a | 1.42±0.07b | 1.59±0.03c | 5.72±0.47d |
| 3-acetyl-2,5-dimethyl-furan | 1093 | MS, RI, S | ND^d | 2.79±0.09a | 3.80±0.17b | 6.26±0.24c |
| 2-furanmethanol | 1685 | MS, RI, S | 2.47±0.12a | 4.61±0.68b | 6.09±0.39c | 7.69±1.06d |
| 2(5H)-furanone | 1769 | MS, RI, S | ND | 0.44±0.05a | 0.70±0.08b | 0.88±0.06c |
| pyrazine | 1211 | MS, RI, S | 2.03±0.28a | 3.45±0.35b | 4.64±0.57c | 6.94±1.02d |
| methyl-pyrazine | 1295 | MS, RI, S | 4.15±0.16a | 6.96±0.83b | 8.80±1.26c | 11.77±1.33d |
| 2,5-dimethyl-pyrazine | 1345 | MS, RI, S | 0.35±0.02a | 0.34±0.09a | 0.39±0.03b | 0.50±0.06c |

**Table S4, 140 °C:**

| Compounds | RI^b | Identification methods^c | 30 min | 60 min | 90 min | 120 min |
|---|---|---|---|---|---|---|
| furan | 813 | MS, RI, S | 10.26±2.14c | 13.80±3.42d | 5.15±0.87b | 3.36±0.28a |
| 2-methyl-furan | 853 | MS, RI, S | 34.45±4.93c | 90.89±8.72d | 21.24±1.76a | 25.73±3.30b |
| 2-furfural | 1457 | MS, RI, S | 163.52±9.34a | 682.61±13.49b | 1165.37±31.92c | 1081.45±21.63c |
| 2,3-butanedione | 989 | MS, RI, S | 643.04±13.90c | 1495.44±27.93d | 487.04±11.02b | 273.20±9.01a |
| 2,3-pentanedione | 1051 | MS, RI, S | 12.97±1.20b | 11.77±1.39b | 1.91±0.03a | 2.14±0.09a |
| 2-vinyl-furan | 1094 | MS, RI, S | 14.83±1.72b | 22.51±2.76d | 16.56±1.03c | 6.88±0.97a |
| 3-acetyl-2,5-dimethyl-furan | 1093 | MS, RI, S | 52.68±3.19c | 52.02±4.24c | 26.52±1.18b | 22.31±2.52a |
| 2-furanmethanol | 1695 | MS, RI, S | 39.67±3.29b | 46.54±2.81c | 26.64±3.10a | 24.06±1.87a |
| 2(5H)-furanone | 1773 | MS, RI, S | ND^d | ND | ND | ND |
| pyrazine | 1213 | MS, RI, S | 30.10±2.91c | 49.14±4.13d | 20.49±2.44b | 5.20±0.49a |
| methyl-pyrazine | 1275 | MS, RI, S | 56.62±2.19c | 72.40±5.73d | 28.87±3.23b | 17.00±1.08a |
| 2,3-dimethyl-pyrazine | 1343 | MS, RI, S | 3.11±0.83c | 2.92±0.32b | 2.34±0.32a | ND |
| 2,5-dimethyl-pyrazine | 1355 | MS, RI, S | 11.08±1.38a | 14.15±1.93b | 20.35±2.49c | 11.76±1.19a |
| 2,6-dimethyl-pyrazine | 1339 | MS, RI, S | 3.59±0.94c | 5.65±0.39d | 1.82±0.27b | 1.03±0.12a |
| 3-ethyl-2,5-dimethyl-pyrazine | 1137 | MS, RI, S | 6.16±0.75a | 6.71±0.53a | 51.66±4.29b | ND |
| 2-ethyl-3,5-dimethyl-pyrazine | 1459 | MS, RI, S | ND | 44.05±3.21a | ND | ND |
| trimethyl-pyrazine | 1393 | MS, RI, S | 10.87±2.18b | 14.13±1.22c | 10.48±1.05b | 5.85±0.38a |

As printed: 3-ethyl-2,5-dimethyl-pyrazine has RI 1137 on DB-WAX, below pyrazine (1213), which cannot be
right for that compound; 2-furanmethanol is ND at all times at 100 °C.

**Table S1, calibration equations** ("The calibration equations for different volatile compounds in the
Xyl-Lys-ARPs model system"; x and y are not defined in the SI):

| Compounds | Calibration equations |
|---|---|
| furan | y = 0.0211 x - 0.0016, R2 = 0.9941 |
| 2-methyl-furan | y = 0.0256 x + 0.0033, R2 = 0.9927 |
| 2-vinyl-furan | y = 0.0324 x - 0.0132, R2 = 0.9904 |
| 2-furfural | y = 0.0188 x - 0.3366, R2 = 0.9901 |
| 2-acetylfuran | y = 0.0112 x + 0.0843, R2 = 0.9936 |
| 3-acetyl-2,5-dimethyl-furan | y = 0.0629 x + 0.0024, R2 = 0.9948 |
| 2(5H)-furanone | y = 0.0598 x - 0.0011, R2 = 0.9993 |
| 2,3-butanedione | y = 0.0132 x + 0.0461, R2 = 0.9913 |
| 2,3-pentanedione | y = 0.0276 x - 0.0021, R2 = 0.9913 |
| pyrazine | y = 0.0023 x + 0.0031, R2 = 0.9988 |
| methyl-pyrazine | y = 0.0041 x + 0.0027, R2 = 0.9976 |
| 2,3-dimethyl-pyrazine | y = 0.0026 x + 0.0048, R2 = 0.9955 |
| 2,5-dimethyl-pyrazine | y = 0.0031 x + 0.0043, R2 = 0.9981 |
| 2,6-dimethyl-pyrazine | y = 0.0036 x + 0.0048, R2 = 0.9970 |
| trimethyl-pyrazine | y = 0.0083 x + 0.0037, R2 = 0.9983 |
| 2-ethyl-3-methyl-pyrazine | y = 0.013 x + 0.0142, R2 = 0.9994 |

Table S1 has no equation for 2-furanmethanol, 3-ethyl-2,5-dimethyl-pyrazine or 2-ethyl-3,5-dimethyl-pyrazine,
although Tables S2-S4 quantify them. It lists 2-acetylfuran, which appears in no SI table.

Derived here from Tables S2-S4:
- **Furfural from the di-ARP alone**: 100 °C peaks at 7.51 µg/L (60 min) and falls to 1.28; 120 °C rises
  monotonically to 61.49 (120 min); 140 °C reaches 1165.37 (90 min). 1165.37 µg/L / 96.08 g/mol = 12.1
  µmol/L = 6.1 × 10⁻⁴ mol per mol di-ARP. With xylose at 140 °C (main-text Table 3) the 30 min value is
  higher (494.07 vs 163.52) but the 60 min value lower (327.03 vs 682.61).
- **The main text's "furans" totals include 2,3-butanedione and 2,3-pentanedione.** 140 °C, 60 min: 13.80
  + 90.89 + 682.61 + 1495.44 + 11.77 + 22.51 + 52.02 + 46.54 = 2415.58 µg/L, against the printed 2415.56.
  100 °C, 30 min: 0.09 + 5.25 + 3.19 + 10.55 + 0.42 + 5.41 = 24.91, as printed. 120 °C, 120 min: 185.43
  against 185.42. Without the two diketones, the 140 °C 60 min furans would be 908.37 µg/L. The pyrazine
  totals match the pyrazine rows (e.g. 120 °C 120 min: 6.94 + 11.77 + 0.50 = 19.21).

### 2x.4 What the lysine data can and cannot say about ARP turnover (derived here)

Free lysine requires both xylulosyl residues to have left, so in the di-ARP-alone system the fraction of
the 40 mmol/L residues lost is **at least** 2·[Lys]/40. This is a lower bound: mono-ARPs that have lost one
residue are not counted, and released lysine is itself consumed (it falls after 30 min at 140 °C). If loss
is first order per residue, then k ≥ −ln(1 − 2[Lys]/40)/t. A stronger but model-dependent bound assumes the
two residues leave independently with the same k and lysine is not consumed: [Lys]/20 ≤ (1 − e^(−kt))², so
k ≥ −ln(1 − √([Lys]/20))/t. Largest bound per temperature, min⁻¹:

| T | from | loose bound (residue fraction) | independent-residue bound |
|---|---|---|---|
| 100 °C | Lys 3.84 at 60 min (loose); 1.52 at 30 min (independent) | 3.6e-3 (−ln(1 − 0.192)/60) | 1.1e-2 (−ln(1 − 0.276)/30) |
| 120 °C | Lys 3.35 at 30 min (Fig. 3, main text) | 6.1e-3 (−ln(1 − 0.1675)/30) | 1.8e-2 (−ln(1 − 0.409)/30) |
| 140 °C | Lys 3.96 at 30 min | 7.4e-3 (−ln(1 − 0.198)/30) | 2.0e-2 (−ln(1 − 0.445)/30) |

The bounds barely rise with temperature because lysine is consumed faster at higher T, so they are
weakest where they matter most. They do not separate the 1-DX branch from the 3-DX branch: lysine release
counts both enolisations, plus any other route that frees the amine.

Engine comparison (live `results/validation/core_prediction_uncertainty.json`, `priors`): `k_arp_dpo` +
`k_arp_tdp` at 145 °C = 10^−2.310 + 10^−1.727 = 4.9e-3 + 1.88e-2 = 2.37e-2 min⁻¹ (derived). The engine's ARP
also leaves by the uncatalysed `k_arp_dpo_th` (centre −1.290, fixed; 5.1e-2 min⁻¹) and `k_arp_tdp_th`
(centre −4.354, fixed), whose barrier is the lumped trunk one and is not quoted here. Scaling the two
catalysed channels to 100 °C with the 85.7 kJ/mol override gives a factor exp[(85700/8.314)(1/373.15 −
1/418.15)] = 19.5 (derived), i.e. 1.2e-3 min⁻¹. That is below this paper's loose 100 °C bound (3.6e-3).
At 140 °C (factor 1.35, i.e. 1.76e-2 min⁻¹) it is above the bound. The tension would have to be at the low
temperature, and it does not survive the caveats: different compound (di-glycated lysine vs the engine's
mono-glycated alanine ARP), unbuffered and unmeasured pH at 100 °C, the engine's pH factors not applied,
and the uncatalysed channels omitted. **Not usable as a constraint; recorded as a shape check only.**

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

**Can it pin `k_arp_dpo` / `k_arp_tdp`?** Still **no**, after reading the SI (verdict re-checked
2026-10-09). Neither the main text nor the SI prints an ARP loss time course at any temperature, a rate
constant or an activation energy. The SI adds the 3-DX, 1-DX, GO, MGO and released-lysine time courses at
100 and 140 °C (Figures S5/S6, section 2x.1) and the di-ARP-alone volatiles at all three temperatures
(Tables S2-S4). These give two things, neither of them a number for either constant:
- **A lower bound on residue loss from released lysine** (section 2x.4): at least ~3.6e-3 min⁻¹ at 100 °C,
  or ~1.1e-2 min⁻¹ under an independent-residue model (derived here). The bound covers both enolisation
  branches together, applies to a di-glycated lysine at drifting, unmeasured pH, and becomes weaker with
  temperature because lysine is consumed. Against the engine's catalysed pair scaled to 100 °C with the
  85.7 kJ/mol override (1.2e-3 min⁻¹, derived from `b8.k_arp_dpo.log10_k_ref_145C` −2.310 and
  `b8.k_arp_tdp.log10_k_ref_145C` −1.727), the loose bound is ~3× higher. Given the compound, pH and
  missing-channel caveats, this is a shape check, not a constraint.
- **Fast deoxyosone turnover.** At 100 °C alone, 60 min, at least 2 × 3.84 = 7.7 mmol/L of residues have
  left, while 3-DX + 1-DX stand at 0.18 + 0.26 = 0.44 mmol/L (derived here). At 140 °C the deoxyosones peak
  by 30-60 min and decay. With the 120 °C main-text data this supports fast 3-DX/1-DX consumption
  (retro-aldol to C2/C3, the engine's `k_dpo_c2c3` / `k_osone_decay` territory), but it still gives no
  separable rate.

Furfural from the di-ARP alone (Tables S2-S4) peaks at 7.51 / rises to 61.49 / reaches 1165.37 µg/L at 100 /
120 / 140 °C, at most 6.1 × 10⁻⁴ mol/mol (derived here). That is a shape check on the TDP -> FUR leg's
temperature response, not a number for it. The ARP is still a di-glycated lysine, not the mono-glycated
alanine ARP the engine carries. One more reason for caution with the SI: some error bars are copied between
panels and temperatures, and Figure S7 contradicts Figure 3 at the shared condition.

## What it does not give

- The loss of the di-ARP (or the mono-ARPs) against time: **not in the main text and not in the SI**
  (checked 2026-10-09). The authors' earlier paper (Zhang et al. 2024b, cited as showing xylose
  accelerating di-ARP degradation) may hold it.
- Any rate constant, Ea, or Arrhenius treatment (main text and SI).
- pH at 100 or 140 °C (the SI's Figures S5/S6 have no pH or A420 panel); buffered pH.
- Branch-resolved formation of 1-DX vs 3-DX. Only standing pools are given, and these are net of fast
  consumption.
- Volatiles for the di-ARP + xylose system beyond main-text Tables 1-3. Molar furfural yields can be
  derived from the µg/L values, but the Table S1 calibration does not define x or y.
- Usable SDs for several Figure S5/S6 panels (copied error bars, section 2x.1).
