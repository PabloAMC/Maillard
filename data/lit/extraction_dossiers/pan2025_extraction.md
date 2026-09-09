# Pan 2025 — EXTRACTION (methionine 40 mg/L + fructose / glucose / sucrose 20 / 15 / 15 g/L in 50 mmol/L citrate pH 6.2, 100 / 120 / 140 C, 30-600 s; zero-order rate constants for methanethiol, methional, dimethyl disulfide and dimethyl trisulfide; muskmelon pectin fractions 0.05-0.2 % at 120 C / 4 min)
### The first printed rate constants for the whole chain methional -> methanethiol -> DMDS / DMTS at three cooking temperatures from a methionine + sugar pot — with the unit of the constants left unprinted.

**Source on disk:** `data/articles/Pan2025.pdf` (10 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Pan2025.txt`); Tables 1 (calibration) and 2 (kinetics) came through
clean and are re-typed below. Figures 1-6 (time courses, pectin dose responses, sugar and methionine
consumption, dicarbonyls, A294/A420) were not read; the few numbers the text quotes from them are
recorded as text values. Supplementary Figs. S1-S3 (pectin composition, molecular weight, A235) are
NOT on disk. "Data will be made available on request."

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of three pectin fractions from muskmelon on the formation of volatile sulfur compounds produced from methionine-sugars Maillard reaction" |
| Authors | Xin Pan, Shuang Bi, Yingying Xu, Yanpei Cai, Fei Lao, Jingya Ai, Wenjiang Dong*, Jihong Wu** (China Agricultural University, Beijing; BTBU; Ningxia University; CATAS Wanning) |
| Venue | Food Chemistry 493 (2025) 145821 (PII S0308814625030729). Received 6 April 2025, revised 20 July 2025, accepted 3 August 2025, online 6 August 2025 |
| DOI | 10.1016/j.foodchem.2025.145821 |
| Naming | "methanthiol" (sic) = methanethiol (MeSH); 3-DG = 3-deoxyglucosone; MGO = methylglyoxal; WSP / CSP / NSP = water-, chelator- (CDTA) and sodium-carbonate-soluble pectin; A = concentration in the zero-order equation A = Kt + A0 |
| Companions | Pan 2021 (Food Chem 343, 128459; the VSC method and the juice precursor study, not on disk); Zhang, Wang & Cao 2023 (`zhang2023_extraction.md`); Chin & Lindsay 1994 (`chin1994_extraction.md`, cited for the metal-mediated MeSH oxidation); Ajandouz & Puigserver 1999 (cited for pseudo-zero-order browning of glucose + methionine) |

## 1. Why it matters

The B19 methionine-chain log (`tasks/data_restructure_plan.md` section 7) closes with "Still to
fetch: Pan 2025 (zero-order rates at 120 C, on the list)". This is that paper. Table 2 prints a
zero-order rate constant K at 100, 120 and 140 C for each of methanethiol, methional, DMDS and
DMTS in a buffered methionine + sugar pot at juice-like concentrations, over 30-600 s. Nothing else
on disk gives the four links of the chain in one run at cooking temperature. The rates are at
[Met] = 0.268 mmol/L with a 240 mmol/L hexose pool, pH 6.2 citrate, so they transport to the
sugar-path lane as a Strecker rate per unit methionine plus three downstream ratios; the three
temperatures give an apparent barrier for each (not printed; derived below). Two things the paper
does not settle: the unit of K (never printed; section 3 reconstructs it from the two printed
end-point concentrations) and the oxidant for MeSH -> DMDS -> DMTS (Chin & Lindsay's metal
mediation is cited; the pot has none added and the reactor atmosphere is not described). The pectin
part (120 C, 4 min) supplies the no-pectin control levels of the four volatiles, methionine
consumption (8 %) and the 3-DG / MGO levels (0.8 / 0.45 µg/mL) — all but the last two figure-only.

## 2. Methods as they matter to a model

- **Pot.** "The mixture that consisted of 40 mg/L methionine and sugars, including 20 mg/mL
  fructose, 15 mg/mL glucose, and 15 mg/mL sucrose, was dissolved in a 50 mM citrate buffer at pH
  6.2." Conversions: **L-methionine 40 mg/L = 0.268 mmol/L**; fructose 20 g/L = 111.0 mmol/L;
  glucose 15 g/L = 83.3 mmol/L; sucrose 15 g/L = 43.8 mmol/L (reducing hexoses 194 mmol/L, plus
  what sucrose hydrolysis adds: "sucrose exhibited a higher consumption rate of approximately 10 %"
  at 120 C / 4 min); citrate 50 mmol/L, pH 6.2. Note the methionine is in mg/L and the sugars in
  mg/mL (Flags 2). Designed "by referencing the authentic muskmelon juice composition".
- **Heating.** "heated at 100 C, 120 C, and 140 C using a four-position parallel reactor
  (YZMR-4*100, Yanzheng Instrument Co., Ltd., Shanghai, China), with sampling at intervals of 30,
  60, 120, 180, 240, 360, 480, and 600 s. All samples were immediately cooled in an ice bath and
  stored at -20 C before analysis." Vessel volume, fill, headspace, pressure and heat-up time are
  not stated (the YZMR-4*100 is a 4 x 100 mL sealed parallel synthesiser). Triplicate.
- **Pectin runs.** WSP / CSP / NSP from muskmelon (Pan 2022 protocol) at 0.05, 0.1, 0.2 % w/v in the
  same pot; commercial low-methoxyl sunflower pectin (GalA > 80 %) as comparison; **120 C for 4 min**
  in the parallel synthesiser, ice bath. The **no-pectin control** is the 0 % member of each Fig. 2
  panel (figure-only); the same pot as the kinetic run at 120 C / 240 s.
- **Volatile sulfur compounds.** 5 mL sample + 2.5 g NaCl in a 20-mL vial, PTFE-silicone septum, 45 C
  10 min equilibration with agitation, DVB/CAR/PDMS 50/30 µm 40 min at 45 C, desorption 250 C 5
  min; Agilent 7890/5975C, DB-Wax 30 m x 0.25 mm x 0.25 µm, He 1 mL/min, 35 C (2 min) -> 45 C at 2
  C/min -> 130 C at 5 C/min -> 225 C at 10 C/min (5 min). Identification: NIST10, RI vs authentic
  standards and NIST WebBook. **Quantification: external standard, SIM**, stock 1500 µg/L in
  methanol, working standards diluted in the **citrate buffer (50 mM, pH 6.2)** — buffer-matched,
  not sugar-matched (Flags 6). Table 1 gives the curves and ranges.
- **Kinetics.** "Concentrations versus reaction time data were plotted to fit zero, first, and
  second-order reaction models, respectively"; zero order retained; K and R2 in Table 2 as
  A = Kt + A0. A0 not printed. **The unit of K is not printed anywhere** (see section 3).
- **Sugars** HPLC-RI (XBridge Amide), external calibration. **Methionine** UPLC-MS/MS, MRM 150 > 133
  and 150 > 104, SPE clean-up. **3-DG and MGO** OPD (1.5 % in methanol), 37 C 24 h dark, HPLC-PDA 315
  nm, external calibration (standards in mobile phase; 3-DG standard >= 75 % purity). **A294, A420**
  after four-fold dilution. All in triplicate, Duncan p < 0.05.
- **Unit conversions used below.** M (g/mol): methanethiol 48.11, methional 104.17, DMDS 94.20, DMTS
  126.26, 3-DG 162.14, MGO 72.06. 1 µmol L-1 s-1 = 60 µmol L-1 min-1 = 6e-2 mmol L-1 min-1.

## 3. Tables re-typed

### Table 1. "Standard curve of volatiles sulfur compounds used for quantification"

| compound | threshold (µg/kg, water, van Gemert) | odor note | quantitative ions (m/z) | standard curve | R2 | range |
|---|---:|---|---|---|---:|---|
| Methanthiol | 0.02 | cooked cabbage | 48 | Y = 201 x + 30,318 | 0.9957 | 0.005-40 µg/L |
| Dimethyl disulfide | 0.16 | cooked potato | 79, 94 | Y = 37,601 x + 169,804 | 0.9990 | 0.005-5 µg/L |
| Dimethyl trisulfide | 0.01 | cooked onion | 79, 126 | Y = 90,949 x - 363,833 | 0.9970 | 0.002-0.05 µg/L |
| Methional | 0.2 | cooked potato | 76, 104 | Y = 256 x - 205,076 | 0.9945 | 1-800 µg/L |

(The odour notes for DMDS and methional are as printed; the DMDS one is presumably swapped with
"cooked cabbage/onion", irrelevant to the numbers.) No LOD/LOQ printed; the lower range limits are
the effective LOQs.

### Table 2. "Kinetic of volatiles sulfur compounds in the chemical model under different heating conditions"

Equation A = Kt + A0 for all four; K printed as "K x 10^-4" (DMTS: "K x 10^-6"), i.e. the tabulated
number times 1e-4 (1e-6) is K; **unit of K and of A not printed**. Time points 30-600 s.

| compound | parameter | 100 C | 120 C | 140 C |
|---|---|---:|---:|---:|
| Methanthiol | K x 10^-4 | 1.579 | 2.342 | 6.335 |
| | R2 | 0.911 | 0.851 | 0.980 |
| Methional | K x 10^-4 | 1.823 | 16.8 | 89.9 |
| | R2 | 0.908 | 0.993 | 0.987 |
| Dimethyl disulfide | K x 10^-4 | 0.015 | 0.070 | 0.187 |
| | R2 | 0.781 | 0.909 | 0.941 |
| Dimethyl trisulfide | K x 10^-6 | 0.026 | 0.121 | 0.054 |
| | R2 | 0.965 | 0.758 | 0.801 |

**Unit reconstruction (mine; the paper is silent).** The text prints two end points: at 140 C /
600 s "methanthiol and methional concentrations peaked at 26.43 µg/L and 626.31 µg/L" = 0.549 and
6.012 µmol/L. Test each candidate unit for K against K x 600 s (an over-estimate of the rise if A0 >
0, an under-estimate if the curve flattens after 4 min as the text says it does):

| candidate unit of K | methional K x 600 s | methanethiol K x 600 s | verdict |
|---|---|---|---|
| µg L-1 s-1 | 5.4 µg/L vs 626 measured | 0.38 µg/L vs 26.4 | no |
| µg L-1 min-1 | 0.09 µg/L | 0.006 µg/L | no |
| mg L-1 min-1 (= µg mL-1 min-1) | 90 µg/L (x 0.14) | 6.3 µg/L (x 0.24) | poor |
| **µmol L-1 s-1** (= nmol mL-1 s-1) | 5.39 µmol/L = **562 µg/L (x 0.90)** | 0.380 µmol/L = **18.3 µg/L (x 0.69)** | **consistent** |
| µmol L-1 min-1 | 9.4 µg/L | 0.6 µg/L | no |

With µmol L-1 s-1 the DMDS 140 C rise is 0.187e-4 x 600 = 0.0112 µmol/L = 1.06 µg/L (inside the
0.005-5 µg/L calibration range) and the DMTS 120 C rise 0.121e-6 x 600 = 7.3e-5 µmol/L = 0.0092 µg/L
(inside 0.002-0.05 µg/L); with any mass unit DMTS falls below its calibration range. **Working
reading: K in µmol L-1 s-1, A in µmol/L, t in s.** This is an inference from two printed numbers,
to be confirmed with the authors or the raw data before the constants enter a fit; the table below
carries the printed numbers and the inferred unit side by side.

**K in the inferred unit, re-expressed per minute (K x 60):**

| compound | 100 C | 120 C | 140 C | unit |
|---|---:|---:|---:|---|
| methional | 1.09e-2 | 0.101 | 0.539 | µmol L-1 min-1 |
| methanethiol | 9.47e-3 | 1.41e-2 | 3.80e-2 | µmol L-1 min-1 |
| dimethyl disulfide | 9.0e-5 | 4.2e-4 | 1.12e-3 | µmol L-1 min-1 |
| dimethyl trisulfide | 1.56e-6 | 7.26e-6 | 3.24e-6 | µmol L-1 min-1 |

(1 µmol L-1 min-1 = 1e-3 mmol L-1 min-1.)

**Apparent Arrhenius barriers (mine, from the three printed K; the authors print none):** methional
125 kJ/mol (R2 0.997); methanethiol 44 kJ/mol (R2 0.93); DMDS 81 kJ/mol (R2 0.99); DMTS
non-monotonic (K falls from 120 to 140 C: "thermal instability ... disproportionation" per the
authors). These are unit-independent (ratios of K) and hold whatever the unit turns out to be.

**Within-study ratios at 120 C (unit-independent if the four rows share a unit):** K(MeSH) /
K(methional) = 2.342 / 16.8 = 0.139; K(DMDS) / K(MeSH) = 0.070 / 2.342 = 0.030; K(DMTS) / K(DMDS) =
0.121e-2 / 0.070 = 0.017. At 100 C: 0.866, 0.0095, 0.0017; at 140 C: 0.070, 0.030, 0.0003. If the
unit is molar, these are molar flux ratios; if it were a mass unit they would need the M ratios
(0.46, 1.96, 1.34 respectively) applied.

### Numbers in the running text (Figs. 1-6 are FIGURE-ONLY except these)

- 140 C, 600 s: methanethiol 26.43 µg/L (0.549 µmol/L), methional 626.31 µg/L (6.01 µmol/L) — the
  maxima of the whole study. DMDS and DMTS "always detected at a low level in all models".
- "the most significant changes were observed within the first 4 min of reaction time in all
  models" (Fig. 1).
- Pectin runs, 120 C, 4 min: fructose and glucose consumption "around 4 %", sucrose "approximately
  10 %", unchanged by pectin (glucose "only a slight change"); **methionine consumption without
  pectin "around 8 %"** (= 21 µmol/L of the 268 µmol/L), up to 16 % with pectin; **3-DG 0.8 µg/mL
  (4.9 µmol/L) and MGO 0.45 µg/mL (6.2 µmol/L) without pectin**; CSP and NSP at 0.2 % cut 3-DG by
  ~86 % and ~83 %; WSP "directly reduce methylglyoxal concentrations by 50 % (data not shown)";
  A294 down with CSP/NSP, A420 up with pectin dose; A235 of WSP / CSP up 13 % / 35 % after heating;
  pectin linearity 0.79 / 1.17 / 0.25 (WSP / CSP / NSP) rising after heating; Rha/GalA of NSP 0.44 ->
  0.31. The volatile levels of the no-pectin control at 120 C / 4 min are not printed (Fig. 2).

## 4. Kinetic numbers the repository can use

Registry mapping: methanethiol -> `methanethiol`; methional -> `methional`; dimethyl disulfide ->
`dimethyl_disulfide`; dimethyl trisulfide -> `dimethyl_trisulfide`; methionine, fructose, glucose,
sucrose, 3-deoxyglucosone, methylglyoxal, pectin -> not in registry.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| methional formation K, 100 / 120 / 140 C | 1.823e-4 / 16.8e-4 / 89.9e-4 (printed); = 1.09e-2 / 0.101 / 0.539 µmol L-1 min-1 if the unit is µmol L-1 s-1 | unit not printed; inferred µmol L-1 s-1 | [Met] 0.268 mmol/L, fructose 111 + glucose 83 + sucrose 44 mmol/L, 50 mmol/L citrate pH 6.2, sealed parallel reactor, 30-600 s | zero order (A = Kt + A0), R2 0.908 / 0.993 / 0.987 | Table 2, p. 4 | measured_rate (unit unresolved) |
| methanethiol formation K, 100 / 120 / 140 C | 1.579e-4 / 2.342e-4 / 6.335e-4 (printed); = 9.5e-3 / 1.41e-2 / 3.80e-2 µmol L-1 min-1 (inferred unit) | as above | same | zero order, R2 0.911 / 0.851 / 0.980 | Table 2 | measured_rate (unit unresolved) |
| DMDS formation K, 100 / 120 / 140 C | 0.015e-4 / 0.070e-4 / 0.187e-4 (printed); = 9.0e-5 / 4.2e-4 / 1.12e-3 µmol L-1 min-1 (inferred) | as above | same | zero order, R2 0.781 / 0.909 / 0.941 | Table 2 | measured_rate (unit unresolved) |
| DMTS formation K, 100 / 120 / 140 C | 0.026e-6 / 0.121e-6 / 0.054e-6 (printed); = 1.6e-6 / 7.3e-6 / 3.2e-6 µmol L-1 min-1 (inferred) | as above | same | zero order, R2 0.965 / 0.758 / 0.801 | Table 2 | measured_rate (unit unresolved; non-monotonic in T) |
| apparent Ea, methional / methanethiol / DMDS | 125 / 44 / 81 | kJ/mol | 100-140 C, three points, from the printed K | Arrhenius on zero-order K | derived from Table 2 (mine) | derived (unit-independent) |
| K(MeSH)/K(methional), K(DMDS)/K(MeSH), K(DMTS)/K(DMDS) at 120 C | 0.139, 0.030, 0.017 | — | same pot | — | derived from Table 2 | within_study_ratio |
| methanethiol, methional at 140 C / 600 s | 26.43, 626.31 (0.549, 6.01 µmol/L) | µg/L | same pot | — | text 3.1 | level_only (the two anchors for the unit) |
| methional yield per methionine at 140 C / 600 s (mine) | 6.01 / 268 = 2.2 % of the initial Met | — | same | — | derived | derived_assumption |
| methionine consumed, no pectin | ~8 % (~21 µmol/L) | — | 120 C, 240 s | — | text 3.3 (Fig. 4) | level_only (rounded text value) |
| 3-DG, MGO, no pectin | 0.8, 0.45 (4.9, 6.2 µmol/L) | µg/mL | 120 C, 240 s | — | text 3.4 (Fig. 5) | level_only |
| fructose, glucose, sucrose consumed | ~4, ~4, ~10 % | — | 120 C, 240 s | — | text 3.3 | level_only |
| the four volatiles vs time at three T (points behind Table 2); no-pectin control levels at 120 C / 4 min; pectin dose responses | — | µg/L | — | — | Figs. 1, 2 | figure_only |
| pectin fraction effects on Met, sugars, 3-DG, MGO, A294, A420 | — | — | 120 C, 4 min | — | Figs. 3-6 | figure_only |

**Reading for programme 6.** At 120 C the methional rate (0.10 µmol L-1 min-1, inferred unit) sits
over [Met] = 0.268 mmol/L and ~200 mmol/L reducing sugar: a pseudo-first-order 3.8e-4 min-1 in
methionine if the sugar-derived dicarbonyl supply is treated as steady (my re-expression, an
assumption). Methanethiol accrues at 14 % of the methional rate, DMDS at 3 % of methanethiol's, DMTS
at 2 % of DMDS's — the disulfides are a small sink at 120 C in a pot with no added metal or oxidant,
which is what B17's oxidant-limited picture predicts for a sealed pot, and Chin & Lindsay 1994
(Cu(II), air) shows what an oxidant does to the same step at 30 C.

## 5. Flags

1. **The unit of K is not printed.** Section 3 reconstructs µmol L-1 s-1 from the two printed end
   points; the alternatives fail by factors of 7 to 100. Treat the unit as inferred until confirmed
   (data "available on request"). The three apparent barriers and all within-study ratios are
   unit-independent and usable now.
2. **Methionine "40 mg/L" against sugars in "mg/mL".** Free methionine in melon juice is tens of
   mg/L, so 40 mg/L (0.27 mmol/L) is the plausible reading and the 2.2 % methional yield is
   consistent with it; but the mixed units in one sentence are a transcription risk. Confirm with
   Pan 2021 or the authors before charging the pot.
3. **Zero order is a 10-min initial-rate statement.** Methionine falls ~8 % and the sugars ~4 % over
   the run, so the data cannot distinguish orders in the reactants; the K are rates at the stated
   concentrations. The authors' own "most significant changes within the first 4 min" and the R2 of
   0.76-0.99 say the curves bend; K over 30-600 s under-reads the early rate.
4. **DMTS falls from 120 to 140 C** (K 0.121e-6 -> 0.054e-6) — the authors invoke DMTS
   disproportionation to DMDS. DMTS needs a sulfur source beyond two MeSH (H2S or polysulfide; Chin
   & Lindsay 1994 saw no DMTS without H2S); this pot has no cysteine, so the DMTS here is the
   smallest quantity in the paper (0.002-0.05 µg/L calibration range).
5. **No oxidant is named or controlled.** The paper cites Chin & Lindsay for transition-metal
   mediation of MeSH oxidation and says "This work focused on the model system and ignored the
   presence of ions". Citrate at 50 mmol/L chelates trace metals; the reactor atmosphere (air in the
   headspace, or not) is not stated. The DMDS / MeSH ratio of 3 % is therefore a no-added-oxidant
   figure.
6. **Calibration matrix.** Standards in citrate buffer without the 50 g/L sugar and 8 % NaCl (the
   NaCl is added to samples and, presumably, standards alike); sugars change headspace partition.
   External standard, no internal standard; systematic level error plausible at tens of %.
7. **Methanethiol is a gas (b.p. 6 C)**: samples were cooled, stored at -20 C, then 5 mL transferred
   to a fresh vial. Losses of MeSH between reactor and vial are not addressed; the MeSH K may be a
   lower bound, and the MeSH/methional ratio with it.
8. **Reactor** volume, headspace, pressure and heat-up time not stated; the 30-s point at 140 C is
   inside the heat-up of a 100-mL vessel.
9. **R2 claim** "exceeding 0.80" is not met by DMDS at 100 C (0.781) and DMTS at 120 C (0.758).
10. **No SD on K**, no A0, no time-zero blank printed; triplicate stated.
11. **Sucrose hydrolysis** (~10 %) adds ~4 mmol/L glucose + fructose during the run: the reducing
    pool is not constant, though the change is small against 194 mmol/L.
12. **No-pectin control levels** of the four volatiles at 120 C / 4 min are figure-only (Fig. 2);
    the kinetic run's 240-s point at 120 C is the same condition and is figure-only too (Fig. 1).
13. **3-DG standard >= 75 % pure** (external calibration) — the 0.8 µg/mL is uncertain by up to a
    third in the same direction.
14. **Registry**: the four volatiles are keyed; methionine, the sugars and the dicarbonyls are not.
