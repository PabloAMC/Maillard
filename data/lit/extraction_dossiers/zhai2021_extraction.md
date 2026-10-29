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

## 2. The kinetics (sec. "Reaction Kinetics" and Fig. 3a; TTCA alone, natural pH 5.5)

> **Correction 2026-10-29.** This heading said "pH 7" until 2026-10-29. The paper runs the kinetic series on TTCA
> "(pH 5.5, 10 mmol/L) without the adjustment of pH" and calls it "the natural pH of 5.5" (text layer, lines
> 220 and 238; the Fig. 3a points coincide with the Fig. 4b pH 5.5 series, sec. 5). B16 and the route
> calibration (stage 1) used these three rows as pH 7 pots; at pH 7 the measured TTCA falls about ten times faster
> (sec. 5). The 100 °C slope is probably 0.0274, not 0.0271 (sec. 5).

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

## 5. Digitised Figs. 4 and 5 (2026-10-29)

**Method.** Both figures are vector graphics in the PDF (`pdfimages -list` shows no raster image on pages
4-5). Page 4 (journal p. 10651, Figs. 3-4) and page 5 (p. 10652, Fig. 5) were converted with `pdftocairo -svg`;
every marker, error-bar half and line segment is a filled path. Marker centres = bounding-box centres of the
marker paths (diamond, square, triangle, circle; Excel triangles are box-centred, checked against the shared
t = 0 point); values by linear interpolation between the outermost tick-mark centres. Coordinates below are
SVG page points (y down). Series identity = Excel draw order (markers drawn in legend order) **and** checked
geometrically: each marker group lies on its own line (max marker-to-line distance 0.21 pt; dashed lines =
TTCA alone, solid = TTCA-Xyl, as in the legends). Error bars: each bar is two filled halves joined at the mean;
SD = stem join to cap centre, read from each half's vertices (cap thickness 0.45 pt in Figs. 3-4, 0.61 pt in
Fig. 5); a zero-length bar is drawn as a bare cap. All numbers below are **[D]** (digitised).

**Conditions as stated in Methods** (verbatim fragments in quotes):
- Systems: "Aqueous TTCA (10 mmol/L) model systems with the presence or absence of Xyl (10 mmol/L)".
- pH: "adjusted to pH 5.5, 6, 7, and 8 using NaOH (2 mol/L) or HCl (2 mol/L)". **No buffer** (unbuffered
  water; pH 5.5 is TTCA's "natural pH" with no adjustment). ⚠ pH not re-measured after heating (not stated);
  for an unbuffered 10 mM system the pH will drift as acids form, so the pH labels are initial pH only.
- Heating: "heated in a thermostatic oil bath at 100, 120, and 140 °C for various durations (20, 40, 60, 80,
  100, and 120 min)"; "rapidly transferred to an ice bath". ⚠ Vessel (type, volume, sealing, headspace) not
  stated; heat-up time not stated (t counts from immersion, presumably). Figs. 4-5 are 100 °C only.
- TTCA, ARP, Cys, Xyl: HPLC-ELSD, XBridge BEH Amide column; "purified compounds prepared in our laboratory as
  external standards" (TTCA, ARP), analytical-grade Cys and Xyl as external standards. TTCA LOD/LOQ =
  2.21/6.63 mg/L.
- Dicarbonyls: OPD derivatisation, HPLC-DAD (SunFire C18). "3-DX quinoxaline (3-DXqx) and 1-DX quinoxaline
  (1-DXqx), were quantitated using diacetyl as an internal standard, and MGO quinoxaline (MGOqx) and GO
  quinoxaline (GOqx) were detected using external standards." ⚠ No authentic 3-DX/1-DX standard: their
  absolute amounts rest on an assumed response factor relative to diacetyl-quinoxaline (MGO and GO are
  calibrated against the commercial 40 % solutions).
- A420: UV-1800; "diluted to reasonable multiples ... (0.2−0.8)". ⚠ Plotted A420 values (0.08-0.98) are
  presumably back-multiplied by the dilution factor; not stated.
- Statistics: "All experiments ... in triplicate ... mean value ± standard deviation". Error bars = SD (n = 3).

### 5.1 Fig. 4 (page 4, lower right)

Caption: "Color development (a) and thermal degradation of TTCA (b) in the absence or presence of
extra-added Xyl under different pH conditions (5.5, 6, 7, and 8), reaction temperature of 100 °C."
Only 100 °C is shown (no 120/140 °C panel); no Cys, Xyl or ARP series in Fig. 4 (Cys and Xyl vs time are in
Fig. 6, ARP in no figure).

**Axis calibration (tick centres, pt).** 4a: x = 0 / 120 min at 351.422 / 487.000 (22.596 pt per 20 min);
y = 0 / 1.0 A420 at 574.324 / 489.311 (17.0025 pt per 0.2). 4b: x = 0 / 120 min at 351.326 / 488.440
(22.852 pt per 20 min); y = 0 / 12 mmol/L at 697.071 / 612.063 (14.168 pt per 2 mmol/L). All intermediate
ticks are equally spaced to < 0.005 pt.

**Read error.** Recovered values land within ±0.0003 mmol/L of two-decimal numbers and ±0.0005 of
three-decimal A420 numbers (the plotted data were evidently entered to 2 / 3 decimals), so the read error is
**≤ ±0.005 mmol/L (4b) and ≤ ±0.001 A420 (4a)**, SDs to the same precision. Times land on 0-120 min within
±0.002 min.

**Fig. 4b — TTCA (mmol/L), 100 °C, mean ± SD [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 10.00 | 10.00 | 10.00 | 10.00 | 10.00 | 10.00 | 10.00 | 10.00 |
| 20 | 9.81 ± 0.095 | 9.78 ± 0.145 | 9.14 ± 0.385 | 7.79 | 9.79 | 9.68 | 8.57 | 8.12 |
| 40 | 9.34 ± 0.395 | 9.18 ± 0.240 | 8.21 ± 0.210 | 5.18 | 9.19 | 8.97 | 7.93 | 6.07 |
| 60 | 9.10 ± 0.290 | 8.79 ± 0.360 | 7.02 ± 0.290 | 3.79 | 8.79 | 8.36 | 6.97 | 2.16 |
| 80 | 8.19 ± 0.245 | 7.97 ± 0.480 | 5.52 ± 0.260 | 1.47 | 7.99 | 7.29 | 4.29 | 0.00 |
| 100 | 7.68 ± 0.165 | 5.98 ± 0.140 | 2.08 ± 0 ⚠ | 0.00 | 7.34 | 5.12 | 0.00 | 0.00 |
| 120 | 6.69 ± 0.370 | 4.37 ± 0.165 | 0.89 ± 0 ⚠ | 0.00 | 6.59 | 3.89 | 0.00 | 0.00 |

Flags (4b):
- ⚠ **Error bars exist only for TTCA-alone pH 5.5, 6 and 7.** The vector file draws exactly three error-bar
  sets, repeated once per series (8 × 3 identical copies); TTCA pH 8 and all four +Xyl series carry none.
  pH 7 at 100 and 120 min is drawn with zero-length bars (SD 0.00 printed, i.e. probably not measured as SD).
- ⚠ "0.00" values are plotted at exactly 0 (LOD 2.21 mg/L ≈ 0.01 mmol/L): read as "< LOD", not as measured
  zero. The pH 8 and +Xyl pH 7 lines stop at 100 min (no 100 -> 120 segment drawn); the 120-min markers sit
  on the axis.
- ⚠ +Xyl pH 5.5 is nearly identical to TTCA pH 6 (9.79/9.78, 9.19/9.18, 8.79/8.79, 7.99/7.97); recorded, not
  interpreted.
- The text says that at pH 8 "the TTCA concentration of the TTCA−Xyl system was even higher than that in the
  TTCA system during the first 40 min" (digitised: 8.12 vs 7.79, 6.07 vs 5.18 ✓) and "reached zero at 80 min"
  (✓); at pH 7 +Xyl "consumed fully at 100 min" (✓). The text's "minor increase first" at pH 7 +Xyl is **not**
  in the plotted data (monotone 10 -> 8.57 -> 7.93 -> 6.97) ⚠.

**Fig. 4a — A420, 100 °C, mean ± SD [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.079 | 0.079 | 0.079 | 0.079 | 0.079 | 0.079 | 0.079 | 0.079 |
| 20 | 0.099 ± 0.013 | 0.103 ± 0.019 | 0.113 ± ⚠ | 0.197 ± 0.038 | 0.100 ± 0.009 | 0.113 ± ⚠ | 0.153 ± 0.018 | 0.229 ± 0.018 |
| 40 | 0.102 ± 0.017 | 0.105 ± 0.011 | 0.154 ± 0.020 | 0.287 ± 0.017 | 0.107 ± 0.019 | 0.121 ± 0.011 | 0.203 ± 0.012 | 0.315 ± 0.013 |
| 60 | 0.132 ± 0.020 | 0.128 ± 0.019 | 0.208 ± 0.019 | 0.318 ± 0.029 | 0.146 ± 0.019 | 0.196 ± 0.019 | 0.248 ± 0.013 | 0.478 ± 0.013 |
| 80 | 0.177 ± 0.019 | 0.197 ± 0.016 | 0.279 ± 0.014 | 0.346 ± 0.026 | 0.218 ± 0.014 | 0.246 ± 0.016 | 0.403 ± 0.012 | 0.590 ± 0.012 |
| 100 | 0.219 ± 0.012 | 0.223 ± 0.006 | 0.309 ± 0.016 | 0.407 ± 0.008 | 0.266 ± 0.017 | 0.298 ± 0.006 | 0.579 ± 0.020 | 0.718 ± 0.014 |
| 120 | 0.248 ± 0.002 | 0.298 ± 0.018 | 0.348 ± 0.014 | 0.498 ± 0.013 | 0.300 ± 0.014 | 0.379 ± 0.018 | 0.745 ± 0.012 | 0.978 ± 0.013 |

Flags (4a): unlike 4b, every series has its own error bars (48 distinct bars). ⚠ At 20 min TTCA pH 7 and +Xyl
pH 6 share the mean 0.113; the two bars there have SD 0.009 and 0.018 and cannot be assigned. All series
start at the same 0.079 at t = 0 (a single shared t = 0 value, zero SD).

### 5.1a Cross-check against the §2 zero-order fit (y = −0.0271 x + 10.331)

- ⚠ **The §2 fits are at pH 5.5, not pH 7.** The text says Fig. 3 is "under the natural pH of 5.5" ("The
  aqueous solution of TTCA (pH 5.5, 10 mmol/L) without the adjustment of pH was heated under different
  temperatures"); §2's heading "TTCA alone, pH 7" is a mislabel (left unedited above).
- **pH 7 / 100 °C (as asked):** the fit predicts 9.79 / 9.25 / 8.70 / 8.16 / 7.62 / 7.08 mmol/L at 20-120 min;
  Fig. 4b pH 7 gives 9.14 / 8.21 / 7.02 / 5.52 / 2.08 / 0.89. Incompatible (residual −6.2 mmol/L at 120 min);
  pH 7 decay is strongly accelerating, not zero-order.
- **pH 5.5 / 100 °C:** Fig. 4b TTCA pH 5.5 means are identical, point for point, to the 100 °C series of
  Fig. 3a (digitised the same way: 10.00, 9.81, 9.34, 9.10, 8.19, 7.68, 6.69 [D]). A least-squares line through
  these 7 points (t = 0 included) gives **y = −0.02739 x + 10.3310, R² = 0.9517**: intercept and R² reproduce
  the printed fit, the slope does not (−0.0274 vs −0.0271). Fig. 3a refits at 120 and 140 °C reproduce the
  printed fits to 4 figures (−0.06513 x + 9.9723, R² 0.9949; −0.08134 x + 9.3380, R² 0.9802). ⚠ The printed
  100 °C slope "0.0271" is most likely a transcription of 0.0274 (the printed intercept is only consistent
  with 0.0274 given the data mean 8.687 at t̄ = 60 min). Effect on Ea: negligible (1 %).
- ⚠ Same means, different SDs: Fig. 3a (100 °C) SDs are 0.048 / 0.279 / 0.179 / 0.214 / 0.098 / 0.119 at 20-120 min,
  Fig. 4b pH 5.5 SDs 0.095 / 0.395 / 0.290 / 0.245 / 0.165 / 0.370 [D]. One data set, two error-bar sets.

Fig. 3a (TTCA, pH 5.5 natural, digitised for this cross-check; tick centres x 85.502 / 167.389 pt for
0 / 120 min, y 451.526 / 378.971 pt for 0 / 12 mmol/L; read error ≤ ±0.005 mmol/L) [D]:

| t (min) | 100 °C | 120 °C | 140 °C |
|---:|---:|---:|---:|
| 0 | 10.00 | 10.00 | 10.00 |
| 20 | 9.81 ± 0.048 | 8.77 ± 0.115 | 7.59 ± 0.285 |
| 40 | 9.34 ± 0.279 | 7.27 ± 0.220 | 5.87 ± 0.210 |
| 60 | 9.10 ± 0.179 | 6.18 ± 0.190 | 3.98 ± 0.180 |
| 80 | 8.19 ± 0.214 | 4.59 ± 0.290 | 2.29 ± 0.190 |
| 100 | 7.68 ± 0.098 | 3.17 ± 0.299 | 1.20 ± 0.220 |
| 120 | 6.69 ± 0.119 | 2.47 ± 0.165 | 0.27 ± 0.180 |

### 5.2 Fig. 5 (page 5, top; four panels)

Caption: "Effect of extra-added Xyl on the generation of characteristic intermediates derived at different
steps of TTCA degradation under different pH conditions (5.5, 6, 7, and 8): (a) 3-DX, (b) 1-DX, (c) MGO,
(d) GO, reaction temperature of 100 °C." Y axes are labelled "Concentration (mmol/L)"; **tables below are
in µmol/L (= mmol/L × 1000)** to avoid long decimals. Legend (each panel): dashed = TTCA-5.5/6/7/8
(diamond/square/triangle/circle), solid = TTCA-Xyl-5.5/6/7/8, same symbols. Series identity confirmed by
draw order and by line geometry (each marker group lies on its own line, max distance 0.29 pt; dashed lines
1-4, solid 5-8). Every series has its own error bars (48 distinct bars per panel, each drawn 8×).

**Axis calibration (tick centres, pt; x = 0 / 120 min, y = 0 / top tick).**

| panel | x0 | x120 | y0 | y top (value) | pt per x tick | pt per y tick |
|---|---:|---:|---:|---|---:|---:|
| a (3-DX) | 173.791 | 298.824 | 199.006 | 96.439 (1.2 mmol/L) | 20.839 | 25.642 (0.3) |
| b (1-DX) | 346.418 | 471.439 | 199.205 | 98.066 (0.2 mmol/L) | 20.837 | 25.285 (0.05) |
| c (MGO) | 174.209 | 303.619 | 365.191 | 263.008 (0.08 mmol/L) | 21.568 | 25.546 (0.02) |
| d (GO) | 348.201 | 472.748 | 365.186 | 264.865 (0.005 mmol/L) | 20.758 | 20.064 (0.001) |

**Read error.** Tick spacing is uniform to ≤ 0.03 pt; recovered times fall on 0-120 min within ±0.003 min and
most means land within ±0.03 µmol/L of round 3-4-decimal mmol/L values. Taking ±0.02 pt as the
calibration-limited error: **a ±0.25, b ±0.04, c ±0.02, d ±0.001 µmol/L** (SDs to the same precision). SDs are
displayed rounded; the ⚠{x/y} cells are points where two series share the same mean and the two bars there
cannot be assigned.

**Fig. 5a — 3-DX (µmol/L), 100 °C, mean ± SD [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 280 ± 17 | 330 ± 14 | 150 ± 23 | 62 ± 22 | 290 ± 19 | 380 ± 8 | 170 ± 19 | 48 ± 74 ⚠ (upper half only) |
| 40 | 370 ± 18 | 480 ± 12 | 210 ± 11 | 72 ± 24 | 400 ± 18 | 510 ± 9 | 230 ± 30 | 76 ± 13 |
| 60 | 520 ± 12 | 590 ± 10 | 380 ± 10 | 86 ± 63 | 580 ± 18 | 630 ± 20 | 360 ± 30 | 107 ± 17 |
| 80 | 650 ± 17 | 750 ± 22 | 470 ± 18 | 98 ± 16 | 700 ± 22 | 810 ± 28 | 570 ± 19 | 147 ± 19 |
| 100 | 640 ± 20 | 890 ± 22 | 510 ± 18 | 107 ± 19 | 750 ± 7 | 1030 ± 33 | 770 ± 19 | 258 ± 32 |
| 120 | 620 ± 21 | 930 ± ⚠{18/30} | 420 ± 19 | 75 ± 23 | 770 ± 18 | 1090 ± 27 | 930 ± ⚠{18/30} | 359 ± 9 |

**Fig. 5b — 1-DX (µmol/L) [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 10.0 ± 2.1 | 17.2 ± 1.8 | 21.0 ± 2.4 | 32.0 ± 2.8 | 12.0 ± 2.4 | 17.9 ± 1.0 | 26.0 ± 2.4 | 31.0 ± 9.3 |
| 40 | 13.0 ± 2.2 | 21.5 ± 1.5 | 29.0 ± 3.7 | 52.0 ± 3.0 | 13.0 ± 2.2 | 23.0 ± 1.1 | 31.0 ± 3.7 | 53.0 ± 1.6 |
| 60 | 25.0 ± 1.5 | 33.7 ± 1.2 | 48.0 ± 3.7 | 74.0 ± 7.8 | 28.0 ± 2.2 | 34.0 ± 2.5 | 53.0 ± 3.7 | 82.0 ± 2.2 |
| 80 | 37.0 ± 2.1 | 49.8 ± 2.7 | 61.0 ± 2.4 | 99.0 ± 2.0 | 40.0 ± 2.7 | 51.0 ± 3.5 | 75.0 ± 2.4 | 119.0 ± 2.4 |
| 100 | 36.0 ± 2.5 | 62.8 ± 2.7 | 79.0 ± 2.4 | 117.0 ± 2.3 | 43.0 ± 0.9 | 68.8 ± 4.1 | 94.0 ± 2.4 | 172.0 ± 4.0 |
| 120 | 32.0 ± 2.6 | 61.9 ± 2.2 | 69.0 ± 3.7 | 106.0 ± 2.9 | 47.0 ± 2.2 | 79.0 ± 3.3 | 125.0 ± 3.7 | 187.0 ± 1.2 |

**Fig. 5c — MGO (µmol/L) [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 8.1 ± 0.5 | 12.9 ± 0.5 | 21.9 ± 0.7 | 30.2 ± 0.7 | 10.2 ± 0.6 | 16.9 ± 0.2 | 19.8 ± 0.6 | 37.4 ± 2.3 |
| 40 | 10.7 ± 0.6 | 17.8 ± 0.4 | 27.9 ± 0.3 | 45.2 ± 0.7 | 12.7 ± 0.6 | 19.9 ± 0.3 | 23.1 ± 0.9 | 55.3 ± 0.4 |
| 60 | 13.2 ± 0.4 | 24.7 ± 0.3 | 34.8 ± 0.3 | 57.4 ± 2.0 | 14.4 ± 0.6 | 21.4 ± 0.6 | 27.9 ± 0.9 | 63.8 ± 0.5 |
| 80 | 13.9 ± 0.5 | 25.2 ± 0.7 | 31.1 ± 0.6 | 49.9 ± 0.5 | 19.9 ± 0.7 | 24.9 ± 0.9 | 31.8 ± 0.6 | 67.2 ± 0.6 |
| 100 | 13.7 ± 0.6 | 24.8 ± 0.7 | 29.9 ± 0.6 | 34.7 ± 0.6 | 21.3 ± 0.15 | 27.8 ± 1.0 | 37.4 ± 0.6 | 59.7 ± 1.0 |
| 120 | 11.4 ± 0.7 | 20.1 ± 0.6 | 26.9 ± 0.6 | 27.9 ± 0.7 | 24.7 ± 0.6 | 30.0 ± 0.8 | 42.5 ± 0.9 | 56.9 ± 0.3 |

**Fig. 5d — GO (µmol/L) [D]**

| t (min) | TTCA pH 5.5 | TTCA pH 6 | TTCA pH 7 | TTCA pH 8 | +Xyl pH 5.5 | +Xyl pH 6 | +Xyl pH 7 | +Xyl pH 8 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 2.49 ± ⚠{0.04/0.08} | 2.39 ± 0.07 | 1.10 ± 0.11 | 1.90 ± 0.11 | 2.63 ± 0.10 | 2.49 ± ⚠{0.04/0.08} | 1.94 ± 0.09 | 2.24 ± 0.37 |
| 40 | 3.06 ± 0.09 | 2.78 ± 0.06 | 1.90 ± 0.05 | 2.15 ± 0.12 | 3.26 ± 0.09 | 2.84 ± 0.04 | 2.33 ± 0.15 | 2.70 ± 0.06 |
| 60 | 3.48 ± 0.06 | 2.94 ± 0.05 | 2.18 ± 0.05 | 2.57 ± 0.11 | 3.89 ± 0.09 | 3.27 ± 0.10 | 2.81 ± 0.15 | 3.15 ± 0.09 |
| 80 | 3.78 ± 0.08 | 3.19 ± 0.11 | 2.34 ± 0.09 | 2.89 ± 0.08 | 4.18 ± 0.11 | 3.68 ± 0.14 | 2.98 ± 0.10 | 3.12 ± 0.09 |
| 100 | 4.08 ± 0.10 | 3.58 ± 0.11 | 2.89 ± 0.09 | 3.47 ± 0.09 | 4.37 ± 0.03 | 3.89 ± 0.16 | 2.87 ± 0.09 | 2.98 ± 0.16 |
| 120 | 4.13 ± ⚠{0.11/0.13} | 3.82 ± 0.09 | 3.09 ± 0.10 | 3.63 ± 0.11 | 4.79 ± 0.09 | 4.13 ± ⚠{0.11/0.13} | 2.73 ± 0.15 | 2.53 ± 0.05 |

Flags (Fig. 5):
- ⚠ 5a, +Xyl pH 8 at 20 min: only the upper half of the bar is drawn (the lower half would cross zero and is
  clipped); the SD 74 µmol/L exceeds the mean (48). 5a TTCA pH 8 at 60 min also has an outsized SD (86 ± 63).
- ⚠ 3-DX and 1-DX are absolute only up to the diacetyl-internal-standard assumption (Methods); MGO and GO are
  externally calibrated. Treat 3-DX/1-DX levels as ±(unknown response-factor) and their **ratios across pH
  and ±Xyl** as the robust content.
- ⚠ 5b SDs are suspiciously repetitive (2.4 µmol/L at eight points, 3.7 at six, 2.2 at six, across different
  series and times); 5b TTCA pH 5.5 and +Xyl pH 5.5 share 13.0 ± 2.2 at 40 min. Recorded, not corrected.
- Text cross-checks: "3-DX was the dominant product" (✓, peak 1.09 mmol/L vs 1-DX ≤ 0.19, MGO ≤ 0.067,
  GO ≤ 0.0048); 3-DX "relatively low under alkaline (pH 8)" (✓, ≤ 0.11 mmol/L without Xyl); pH 5.5 3-DX "lower
  than that under 6 and even tended to decrease" (✓, 650 -> 620); 1-DX "higher under pH 8" (✓); "MGO ...
  much higher than ... GO" (✓, 4-22× depending on pH); +Xyl raises GO "during the initial 0−60 min under pH 7 and 8" — true
  for pH 8 (2.24/2.70/3.15 vs 1.90/2.15/2.57) and pH 7 (1.94/2.33/2.81 vs 1.10/1.90/2.18), then +Xyl pH 7/8 GO
  falls below TTCA-alone at 100-120 min (✓ "decrease in GO after ... 60 min").
- Yield scale (for the model): at pH 6, 120 min, 3-DX 0.93 mmol/L against 5.63 mmol/L TTCA lost (Fig. 4b)
  = 17 % molar; at pH 5.5, 0.62 vs 3.31 = 19 %; at pH 7, 0.42 vs 9.11 = 5 %; at pH 8, 0.075 vs 10.0 = 0.8 %
  [derived from D].

### 5.3 Derived from the digitised TTCA curves (Fig. 4b) [derived from D]

| series (100 °C) | t½ (min, linear interp.) | −ln(c/c0)/t, 0-20 min (/min) | −ln(c/c0)/t, 0-40 min (/min) | mean rate 0-40 min (mmol/L/min) |
|---|---:|---:|---:|---:|
| TTCA pH 5.5 | > 120 | 0.0010 | 0.0017 | 0.0165 |
| TTCA pH 6 | 112 | 0.0011 | 0.0021 | 0.0205 |
| TTCA pH 7 | 83 | 0.0045 | 0.0049 | 0.0447 |
| TTCA pH 8 | 43 | 0.0125 | 0.0164 | 0.121 |
| +Xyl pH 5.5 | > 120 | 0.0011 | 0.0021 | 0.0203 |
| +Xyl pH 6 | 102 | 0.0016 | 0.0027 | 0.0257 |
| +Xyl pH 7 | 75 | 0.0077 | 0.0058 | 0.0518 |
| +Xyl pH 8 | 45 | 0.0104 | 0.0125 | 0.098 |

Every TTCA-alone curve **accelerates** (apparent rate constant rises with time; pH 7: 0.0045 /min over 0-20
min, 0.049 /min over 80-100 min), i.e. neither zero- nor first-order: an induction phase (TTCA ⇌ open-chain
Schiff base/ARP before an irreversible step), autocatalysis, or pH drift in the unbuffered solution can each
produce this; the paper does not distinguish them. ⚠ The "zero-order" fits of §2 hold only at pH 5.5, where
curvature is mild. Over pH 5.5 -> 8 the 0-40 min first-order-equivalent rises ~10× (0.0017 -> 0.0164 /min),
roughly ∝ [OH⁻]^0.4 (2.5 pH units). Extra xylose speeds loss modestly at pH 5.5-7 and slightly *slows* it at
pH 8 during the first 40 min (8.12/6.07 vs 7.79/5.18), then speeds it.

### 5.4 What this can constrain

Fig. 4b gives the pH dependence (5.5-8, 100 °C) of TTCA loss from a clean 10 mM start, with and without
equimolar xylose — directly usable to fit a pH factor on the sulfur lane's `k_ttca_deg` (≈10× from pH 5.5 to 8
in the first 40 min) and to test whether the lane's TTCA ⇌ Cys + Xyl reversibility reproduces the induction
and the small +Xyl effect; Fig. 3a (pH 5.5) supplies the temperature series, now recovered as raw points, with
the printed 100 °C slope corrected to 0.0274. Fig. 5 constrains the branching of TTCA carbon into 3-DX (5-19 %
molar of TTCA lost, falling steeply with pH), 1-DX, MGO and GO at each pH, i.e. the 1,2- vs 2,3-enolisation
split and C2/C3 fragmentation vs pH. Limits: single temperature (100 °C) for all pH data, unbuffered and
unre-measured pH, unstated vessel/heat-up, 3-DX/1-DX absolute scale tied to a diacetyl surrogate standard,
and error bars missing for five of eight TTCA series.


