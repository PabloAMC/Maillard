# Xia 2022 — EXTRACTION (glucose + glycine / diglycine / triglycine, 0.2 + 0.2 mol/L in water, pH 7.5 initial, 100-130 C, 5-80 min; Amadori formation rates at six temperatures with Ea; 1-DG and 3-DG formation rates at 120 and 130 C; glyoxal and methylglyoxal at 130 C figure-only)
### The trunk's own pot (glucose + glycine, 0.2 M each) from the Zhou/Zhang laboratory: a six-temperature Amadori ladder with a barrier, deoxyglucosone rates at two temperatures, and the statement that glyoxal far exceeds methylglyoxal in water at 130 C — with the dicarbonyl curves themselves not printed.

**Source on disk:** `data/articles/Xia2022.pdf` (12 pp., owner's download, 2026-09-08). Read from
the text layer (`scratchpad/articles/Xia2022.txt`, clean); Tables 1 and 2 were re-extracted with
`pypdf` from page 8 and matched the text layer digit for digit. Page 9 (journal p. 14915) was
rasterised at 75 dpi to read the axis labels of Figure 5 only (all concentration panels
"Concentration (mmol/L)", x "Time (min)"; panel E "A420"). No value was read off any figure. The
Supporting Information (Figures S1-S5: 1-DG/3-DG MS identification, NMR spectra, the 120-130 C
deoxyglucosone scatter plots, the Arrhenius plot) is NOT on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Glycine, Diglycine, and Triglycine Exhibit Different Reactivities in the Formation and Degradation of Amadori Compounds" |
| Authors | Xue Xia, Yun Zhai, Heping Cui, Han Zhang, Khizar Hayat, Xiaoming Zhang*, Chi-Tang Ho* (Jiangnan University, Wuxi; Miami University; Rutgers) |
| Venue | J. Agric. Food Chem. 2022, 70, 14907-14918. Received September 24 2022, revised October 31 2022, accepted November 1 2022, published November 15 2022 |
| DOI | 10.1021/acs.jafc.2c06639 |
| Naming | G-ARP / DiG-ARP / TriG-ARP = the glucose Amadori compounds of glycine, diglycine, triglycine; 1-DG, 3-DG = 1- and 3-deoxyglucosone; GO = glyoxal; MGO = methylglyoxal; Gly / DiGly / TriGly = the free amines |
| Companions | same laboratory as Zhou 2023 / 2024 / 2025 (`zhou2023_extraction.md`, `zhou2024_extraction.md`, `zhou2025b_extraction.md`); ref 21 = Xia 2022 Food Chem. 390, 133144 (the OPD method's origin, alanyl-glutamine ARP); ref 18 = Zhan 2020 (glycine-ribose ARP kinetics, Ea 30 kJ/mol quoted); ref 24 = Yu 2018 (G-ARP Ea 64.8 at pH 10 quoted) |

## 1. Why it matters

Wave B18 (`results/validation/kinetic_core_b18_prereg.md` section 6) says the sugar path makes far
too little glyoxal and methylglyoxal in water at 70 to 120 C and makes its glyoxal only through a
dry-glass glucosone entry. This paper heats the trunk's exact pot — **glucose + glycine, 0.2 + 0.2
mol/L, water, pH 7.5 initial** — at 100 to 130 C and reports, for the glycine arm:

- the Amadori formation rate at 100 / 105 / 110 / 115 / 120 / 130 C with an Arrhenius barrier of
  **84.76 kJ/mol** (the trunk's `k_schiff` carries 96.8 +/- 2.8 from Martins 2005 at pH 6.8);
- 1-DG and 3-DG net formation rates at 120 and 130 C (a two-point barrier each, mine);
- the statement that at 130 C **glyoxal "was generated far more than" methylglyoxal**, with the GO
  panel of Figure 5 drawn to a 15 mmol/L full scale and the MGO panel to 0.5 mmol/L.

The last point is the direct check the brief asks for, and it comes back as a direction, not a
number: the GO and MGO time courses are figure-only, the paper prints no GO:MGO ratio, and the
absolute GO levels implied by the axis scale (mmol/L from 200 mmol/L glucose) are large enough to
need the derivatisation caveat in Flag 3 before any of it is believed.

## 2. Methods as they matter to a model

- **Kinetic pots.** "Glucose and Gly/DiGly/TriGly were added to 5 mL of Milli-Q water with the same
  molarity of 0.2 mol/L. The reaction solution was controlled at pH 7.5 +/- 0.1 at 25 C with NaOH
  solution (6 M)." So **[Glc] = [amine] = 200 mmol/L, water, no buffer, pH 7.5 initial at 25 C**.
  10 mL pressure-resistant glass bottle with PTFE cork, magnetic stirring, oil bath; **100, 105,
  110, 115, 120, 130 C**; **5, 10, 15, 20, 25, 30, 35, 45, 60, 80 min**; ice quench. Heat-up time
  not stated; pH during or after heating not reported.
- **Amadori compounds (HPLC-ELSD, no derivatisation).** Waters 1525 + Alltech 3300 ELSD, Xbridge
  Amide 4.6 x 150 mm 3.5 um, water:acetonitrile 13:7 with 0.1 % formic acid, 0.042 L/h = 0.7
  mL/min, 10 uL, drift tube 55 C. Log-log calibrations (A = peak area, C in mmol/L): G-ARP
  lg A = 2.9013 lg C + 12.2480 (R2 0.9954, LOD 0.25, LOQ 0.40 mmol/L); DiG-ARP 1.5041 lg C +
  10.2210 (LOD 0.23, LOQ 0.30); TriG-ARP 2.6911 lg C + 12.2380 (LOD 0.10, LOQ 0.18); Gly 2.7068
  lg C + 11.2910 (LOD 0.15, LOQ 0.26); DiGly 1.3908 lg C + 9.6736 (LOD 0.07, LOQ 0.13); TriGly
  2.6557 lg C + 11.9030 (LOD 0.05, LOQ 0.13); glucose 2.6842 lg C + 11.6520 (LOD 0.09, LOQ 0.15).
  ELSD response is power-law, hence the log fits. Standards were the laboratory's own purified ARPs
  (purity by NMR: G-ARP 96.40 %, DiG-ARP 98.03 %, TriG-ARP 96.07 %).
- **ARP preparation** (not the kinetic run): glucose:amine 2:1, amine 0.1 mol/L, 100 mL water,
  pH 7.5, 80 C water bath 60-100 min, rotary evaporation 80 C, Dowex column, 0.2 mol/L ammonia
  elution, lyophilised. Yields 77.67 / 80.97 / 56.33 %.
- **Dicarbonyls (OPD-quinoxalines, HPLC-PDA).** Reagent: OPD 4.6 mmol + DTPA 1.1 mmol in 100 mL
  of 1 mol/L HEPES (46 mmol/L OPD, 11 mmol/L DTPA). **0.2 mL sample + 0.2 mL reagent, dark, 25 C,
  12 h.** Waters e2695 + 2998 PDA at 315 nm, SunFire C18 4.6 x 150 mm 5 um, 35 C, 10 uL, water
  (0.1 % formic acid) : methanol 7:3 -> 2:3 over 3-10 min, 0.8 mL/min. Linear calibrations (x in
  mmol/L, y = peak area): **GO y = 1.4377e6 x + 1.0039e5 (R2 0.9987, LOD 0.033, LOQ 0.041 mmol/L);
  MGO y = 2.7247e6 x + 2.1453e4 (R2 0.9998, LOD 0.013, LOQ 0.025 mmol/L); 3-DG y = 8.8101e5 x +
  7.9959e3 (R2 0.9965, LOD 0.0012, LOQ 0.0038 mmol/L).** GO and MGO from 40 % commercial
  solutions; 3-DG 95 %. **1-DG has no standard and was semi-quantified on the 3-DG curve**
  ("similar proton-accepting groups"); its identity rests on UPLC-MS/MS of the quinoxaline (Figure
  S1). Note the LOQ of GO (0.041 mmol/L) is a thousand times the LOQ of 3-DG.
- **Browning**: A420 (dilution not stated).
- **Kinetic treatment.** "Pseudo-first-order fit model" over an initial window, but the printed unit
  of every k is **mmol/L per min** — a slope of concentration against time, i.e. a zero-order
  (initial-rate) quantity at fixed 200 + 200 mmol/L. The window over which the line was fitted:
  G-ARP 45 min at 100-115 C, 35 min at 120 and 130 C; DiG-ARP 45 min (100-115), 35 min (120),
  25 min (130); TriG-ARP 45 min (100-110), 25 min (115, 120), 20 min (130). ARP maxima: G-ARP at
  60 min (120 C) and 45 min (130 C); DiG-ARP 45 min (120 C), 25 min (130 C); TriG-ARP 80 min
  (100, 105 C), 45 min (110-120 C), 25 min (130 C). Arrhenius: ln k vs 1/T on the six Table 1
  points (Figure S5G, not on disk).
- **Unit conversion.** 1 mmol/L/min = 1e-3 mol/L/min = 1.667e-5 mol L-1 s-1.
- **Replication.** Mean +/- SD; the number of replicates is not stated.

## 3. Tables re-typed

### Table 1. "Kinetic Parameters (k and Ea) for the Formation of G-ARP, DiG-ARP, and TriG-ARP"

The printed column header reads "100 105 110 105 120 130" — the fourth column is **115 C** (the
text lists 100, 105, 110, 115, 120, 130 C). k in mmol/L·min-1.

| system | quantity | 100 C | 105 C | 110 C | 115 C | 120 C | 130 C | Ea (kJ/mol) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| G-ARP (glycine) | k | 0.0604 | 0.0845 | 0.1692 | 0.2148 | 0.3026 | 0.4300 | 84.76 |
| | R2 | 0.9695 | 0.9573 | 0.9876 | 0.9918 | 0.9916 | 0.9765 | |
| DiG-ARP | k | 0.0954 | 0.1606 | 0.2079 | 0.2660 | 0.4032 | 0.5610 | 72.84 |
| | R2 | 0.9279 | 0.9954 | 0.9551 | 0.9901 | 0.9636 | 0.9369 | |
| TriG-ARP | k | 0.4153 | 0.6159 | 0.6354 | 0.9989 | 1.1130 | 2.1561 | 63.48 |
| | R2 | 0.9663 | 0.9602 | 0.9929 | 0.9623 | 0.9368 | 0.9694 | |

Printed Arrhenius lines (text): G-ARP ln k = -10194/T + 24.6219 (R2 0.9541); DiG-ARP ln k =
-8760/T + 21.2537 (R2 0.9764); TriG-ARP ln k = -7635/T + 19.6254 (R2 0.9682).

**Arithmetic.** 10194 x 8.314 = 84.75 kJ/mol; 8760 x 8.314 = 72.83; 7635 x 8.314 = 63.48 — the
printed Ea are the printed slopes. **Re-fit (mine, unweighted, six points):** G-ARP slope -10194,
intercept 24.6219, R2 0.9541 — identical to the printed line; DiG-ARP slope -8773 (Ea 72.94),
R2 0.9756 (printed -8760, R2 0.9764: rounding); TriG-ARP slope -7928 (Ea 65.9), intercept 20.37,
R2 0.9723 (printed -7635, Ea 63.48, R2 0.9682) — the TriG line was not fitted to exactly these six
values, but the difference is 2.4 kJ/mol. The glycine line reproduces exactly, which is the one
the trunk cares about.

**G-ARP unit table:**

| T (C) | k (mmol/L/min) | mol L-1 s-1 | k / ([Glc][Gly]) = k / (200 x 200), L mmol-1 min-1 (my second-order re-expression) |
|---:|---:|---:|---:|
| 100 | 0.0604 | 1.01e-6 | 1.51e-6 |
| 105 | 0.0845 | 1.41e-6 | 2.11e-6 |
| 110 | 0.1692 | 2.82e-6 | 4.23e-6 |
| 115 | 0.2148 | 3.58e-6 | 5.37e-6 |
| 120 | 0.3026 | 5.04e-6 | 7.56e-6 |
| 130 | 0.4300 | 7.17e-6 | 1.08e-5 |

### Table 2. "Formation Rate (k) of 1-DG and 3-DG in Glucose and Gly, DiGly, or TriGly Systems"

k in mmol/L·min-1; the 120 C scatter plots are in Figure S5 (not on disk), the 130 C plots in
Figure 5A, B.

| system | quantity | 1-DG 120 C | 1-DG 130 C | 3-DG 120 C | 3-DG 130 C |
|---|---|---:|---:|---:|---:|
| Gly | k | 0.0055 | 0.0094 | 0.0131 | 0.0357 |
| | R2 | 0.9991 | 0.9962 | 0.9790 | 0.9959 |
| DiGly | k | 0.0081 | 0.0437 | 0.1002 | 0.2944 |
| | R2 | 0.9661 | 0.9674 | 0.9920 | 0.9972 |
| TriGly | k | 0.0057 | 0.0199 | 0.0894 | 0.1717 |
| | R2 | 0.9882 | 0.9876 | 0.9955 | 0.9982 |

**Arithmetic.** Two-point barriers (mine, 120 -> 130 C): Gly 3-DG **132 kJ/mol**, Gly 1-DG **71
kJ/mol**; DiGly 3-DG 142, 1-DG 222; TriGly 3-DG 86, 1-DG 165. Two points ten degrees apart with a
semi-quantified 1-DG: treat as order-of-magnitude only. Within-study ratios as printed: 3-DG rate
at 130 C DiGly/Gly = 8.25, TriGly/Gly = 4.81; 1-DG DiGly/Gly = 4.65, TriGly/Gly = 2.12; 1-DG
130/120 C DiGly 5.40 (text "4.40 times"; the table gives 0.0437/0.0081 = 5.40), TriGly 3.49 (text
"2.49"); Gly 3-DG/1-DG = 2.38 (120 C), 3.80 (130 C). The text's "increased by 4.40 times and 2.49
times" reads as (ratio - 1), i.e. an increase *by* 4.40x = 5.40 times; kept as printed with this
note.

### Numbers and directions in the text (curves FIGURE-ONLY)

- GO and MGO were measured at 130 C only (Figure 5C, D). Glycine arm: GO "at a low level at 0-35
  min and gradually increased at 45-80 min"; MGO "increased slowly with time". Orders: GO TriGly >
  DiGly > Gly; MGO DiGly > TriGly > Gly; "in Gly, DiGly, and TriGly systems, the formation of MGO
  was lower than that of GO"; conclusion: "GO was generated far more than MGO due to the multiple
  formation pathways, while the generation of MGO was more related to the degradation of 3-DG."
- Figure 5 axis extents (raster, labels only): panel C (GO) 0-15 mmol/L; panel D (MGO) 0-0.5
  mmol/L; panels A (1-DG) 0-3 and B (3-DG) 0-15 mmol/L; x 0-80 min. No point was read.
- Browning A420 order DiGly > TriGly > Gly at 130 C.
- pK2 quoted: Gly 9.77, DiGly 8.25, TriGly 7.91 (the authors' explanation of the order).
- Comparators quoted: G-ARP Ea 64.8 kJ/mol at pH 10 (Yu 2018); glycine-ribose ARP Ea 30 kJ/mol at
  pH 7.5 (Zhan 2020); glutathione-xylose ARP 79.35 kJ/mol at 50-80 C (Tang 2019).

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): **glyoxal, methylglyoxal, 3-deoxyglucosone,
1-deoxyglucosone, the Amadori compounds -> not in registry**; browning -> no molecule row.

| step | quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|---|
| Glc + Gly -> G-ARP (net accumulation) | k at 100/105/110/115/120/130 C | 0.0604 / 0.0845 / 0.1692 / 0.2148 / 0.3026 / 0.4300 | mmol L-1 min-1 | 200 + 200 mmol/L, water, pH 7.5 initial, initial window 35-45 min | Table 1 | measured_rate (initial slope; see Flag 1) |
| Glc + Gly -> G-ARP | Ea | 84.76 | kJ/mol | same, 100-130 C, six points | Table 1 / text line | measured_barrier (R2 0.954; no interval printed); trunk `k_schiff` 96.8 +/- 2.8 (Martins, pH 6.8, phosphate) |
| Glc + DiGly -> DiG-ARP | k, Ea | 0.0954 ... 0.5610; 72.84 | mmol L-1 min-1; kJ/mol | same | Table 1 | measured_rate / measured_barrier (peptide comparator) |
| Glc + TriGly -> TriG-ARP | k, Ea | 0.4153 ... 2.1561; 63.48 (refit 65.9) | same | same | Table 1 | measured_rate / measured_barrier (peptide comparator) |
| Glc + Gly -> G-ARP | second-order re-expression (mine) | 1.51e-6 (100 C) ... 1.08e-5 (130 C) | L mmol-1 min-1 | assumes rate = k2[Glc][Gly], both 200 mmol/L unconsumed | derived | derived_assumption; Martins `k_schiff` X = 1.6e-5 at 100 C (Flag 2) |
| Amadori -> 3-DG (net, Gly arm) | k | k(120 C) = 0.0131; k(130 C) = 0.0357 | mmol L-1 min-1 | same pot | Table 2 | measured_rate (net slope) |
| Amadori -> 1-DG (net, Gly arm) | k | k(120 C) = 0.0055; k(130 C) = 0.0094 | mmol L-1 min-1 (on the 3-DG calibration) | same pot | Table 2 | measured_rate, semi-quantitative |
| Amadori -> 3-DG vs -> 1-DG | rate ratio, Gly arm | 2.38 (120 C), 3.80 (130 C) | ratio | same | derived from Table 2 | within_study_ratio (1-DG semi-quantitative) |
| Amadori -> 3-DG (net, Gly) | two-point barrier (mine) | 132 | kJ/mol | 120-130 C | derived | derived (two points; order of magnitude) |
| Amadori -> 1-DG (net, Gly) | two-point barrier (mine) | 71 | kJ/mol | 120-130 C | derived | derived (two points; semi-quantitative species) |
| glyoxal vs methylglyoxal supply | GO vs MGO in the Gly pot at 130 C | GO >> MGO ("far more"); GO panel drawn to 15 mmol/L, MGO panel to 0.5 mmol/L | direction | 130 C, 0-80 min | text + Figure 5 axis extents | level_only (direction); the ratio itself figure_only |
| glyoxal, methylglyoxal | concentrations vs time, all three amines | — | mmol/L (axis) | 130 C | Figure 5C, D | figure_only |
| 1-DG, 3-DG | concentrations vs time | — | mmol/L (axis) | 130 C (Fig 5A, B); 120 C (Fig S5, not on disk) | figure | figure_only |
| Amadori | ARP concentrations vs time, six temperatures | — | mmol/L | 100-130 C | Figure 4 | figure_only |
| method | GO / MGO / 3-DG calibration LOQ | 0.041 / 0.025 / 0.0038 | mmol/L | OPD-HEPES method | Methods | method fact |

Cross-reference inside the repo: `src/kinetic_core/parameters.py` `MARTINS_M4` step 1 (`k_schiff`,
X = 1.6e-5 L mmol-1 min-1 at 100 C, Ea 96.8 +/- 2.8) is the trunk's Amadori entry. This paper's
glycine barrier (84.8) sits 12 kJ/mol below it and 20 above Yu 2018's pH 10 value (64.8): the
three barriers fall monotonically with pH (6.8 / 7.5 / 10 -> 96.8 / 84.8 / 64.8), but the three
pots also differ in buffer (0.1 M phosphate / none / none) and Yu 2018 used ultrasound, so the pH
reading is suggestive only (Flag 2 on the rate difference). For B18's question the usable content is the direction GO >> MGO in
the trunk pot at 130 C, and the 3-DG and 1-DG net rates at 120 / 130 C against which the trunk's
`k_ama_tdg` x [DFG] and `k_ama_odg` x [DFG] fluxes can be checked once the trunk's DFG at 120 C is
integrated.

## 5. Flags

1. **"Pseudo-first-order" in words, zero-order in units.** Every k is printed in mmol L-1 min-1,
   the slope of concentration against time over an initial window, at fixed 200 + 200 mmol/L. The
   Table 1 numbers are initial *net accumulation* rates of ARP (formation minus degradation over
   35-45 min); at 130 C the ARP peaks at 45 min, so even the 35 min window is not free of
   degradation. The Ea of 84.76 kJ/mol is the barrier of that net rate. Store the k values with the
   concentrations attached; the second-order re-expression in section 4 is my assumption.
2. **Ten times slower than Martins' pot.** Martins 2005 step 1 gives 1.6e-5 x 200 x 200 = 0.64
   mmol/L/min into the Amadori pool at 100 C, pH 6.8, 0.1 M phosphate; with Martins' own DFG loss
   (0.034 /min, half-life 20 min) the net DFG slope over the first 45 min is still ~ 0.33
   mmol/L/min. Xia's 0.0604 is 5 to 10 times lower in unbuffered water at pH 7.5. Phosphate is a
   known catalyst of the Amadori rearrangement; the difference is consistent with buffer catalysis
   and is a reason not to move the trunk's constant on this paper alone.
3. **Glyoxal at mmol/L from a glucose + glycine pot is a large number.** The GO panel's 15 mmol/L
   scale (with the peptide arms high on it) would mean several per cent of the glucose as free
   glyoxal within 80 min at 130 C, far above the umol/L-to-low-mmol/L glyoxal that other aqueous
   glucose-amine studies on disk report. The OPD derivatisation here runs 12 h at 25 C on an
   unquenched pot (sugar, ARP, deoxyglucosones all present, 46 mmol/L OPD, no ice-cold quench or
   pH lowering described), a condition known to generate GO and MGO from precursors during
   derivatisation; DTPA is present, which suppresses metal-catalysed autoxidation but not the
   retro-aldol routes. The GO >> MGO direction may partly be a derivatisation-time artefact
   (glyoxal is the dicarbonyl most readily produced from glycolaldehyde/ARP oxidation during a long
   incubation). Take the direction as a hypothesis, not a measurement, until a short-incubation
   study in the same pot exists.
4. **1-DG is semi-quantitative** (3-DG calibration), so the 1-DG rates and the 3-DG/1-DG ratio carry
   an unknown response-factor error.
5. **Unbuffered; pH drift not reported.** Initial pH 7.5 at 25 C set with NaOH; acids formed during
   80 min at 130 C will have lowered it substantially, and the ARP degradation split (1,2- vs
   2,3-enolisation) follows that drift. The authors themselves attribute the 3-DG > 1-DG order to
   "the decreasing pH as the reaction progressed".
6. **Table 1 header typo** (the fourth temperature is 115 C, printed "105"). **Text-vs-table**: the
   1-DG 130/120 C increase is quoted as "4.40 times" and "2.49 times" where the table ratios are
   5.40 and 3.49. Table 1 and 2 numbers are the ones to keep.
7. **TriG-ARP Arrhenius line does not reproduce from the six printed k** (refit Ea 65.9 vs printed
   63.48); the glycine and diglycine lines do.
8. **Number of replicates not stated**; "mean +/- SD" throughout; no confidence intervals on k or
   Ea. Heat-up time in a 10 mL bottle in an oil bath not stated.
9. **Not measured**: glucose consumption in the kinetic runs, any acid, any Strecker product, GO/MGO
   below 130 C, glucosone.
10. **SI not on disk** (Figure S5 has the 120 C deoxyglucosone points and the Arrhenius plot).
11. Registry gap: glyoxal, methylglyoxal, 3-deoxyglucosone and 1-deoxyglucosone have no molecule
    row in `compounds.yml`.
