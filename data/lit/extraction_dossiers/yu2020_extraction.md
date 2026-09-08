# Yu 2020 — EXTRACTION (D-glucose + L-lysine, 0.15 mol/L each, 0.1 mol/L phosphate pH 7.0, 100/110/120 C, 0-60 min; nine-step multiresponse model with 3-deoxyglucosone and methylglyoxal as the only dicarbonyls; furosine, pyrraline, melanoidins as products)
### A glucose + lysine multiresponse fit in water at the trunk's temperatures with a printed barrier for every step, but every concentration is figure-only and the time unit of the rate constants is not printed.

**Source on disk:** `data/articles/yu2020.pdf` (8 pp., owner's download, 2026-09-08). The text layer
(`scratchpad/articles/yu2020.txt`) has its words run together; every number below was checked
against a fresh `pypdf` extraction of the same pages, and Table 1 came through identically both
ways. Page 5 (Fig. 3, Table 2) was rasterised at 75 dpi to read the axis labels of Fig. 3 only
(y "Concentration (mmol/L)", x "Reaction time (s)"). No value was read off any figure. No
Supplementary data are on disk (the paper points to the DOI for them).

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetic study on the generation of furosine and pyrraline in a Maillard reaction model system of D-glucose and L-lysine" |
| Authors | Hang Yu, Qili Zhong, Yunfei Xie, Yahui Guo, Yuliang Cheng, Weirong Yao* (Jiangnan University, Wuxi) |
| Venue | Food Chemistry 317 (2020) 126458. Received 2 August 2019, revised 3 February 2020, accepted 19 February 2020, online 20 February 2020 |
| DOI | 10.1016/j.foodchem.2020.126458 (PII S0308814620303204 per the task brief) |
| Naming | Amd = "Amadori compounds" (the paper's single pool for mono- and difructosyl-lysine, MFL and DFL); 3-DG = 3-deoxyglucosone; MG = methylglyoxal; GO = glyoxal (named in the introduction, never measured); "CEL" is mis-expanded as N-epsilon-(hydroxyethyl)lysine in the abbreviations |
| Companions | Yu, Seow, Ong & Zhou 2018, Food Chem. 269, 628 (glucose + glycine, pH 10, ultrasound; the same reparameterised Arrhenius and the same OPD method; in the repo as `yu2018_extraction.md`, whose Flag 1, "time unit of k not printed", recurs here); Martins & van Boekel 2005 (the trunk, `martins2005_extraction.md`); Kocadagli & Gokmen 2016 (cited for the method's name only) |

## 1. Why it matters

Wave B18 (`results/validation/kinetic_core_b18_prereg.md` section 6) found that from a sugar + amine
pot in water at 70 to 120 C the trunk makes far too little glyoxal and methylglyoxal. This paper is
the nearest thing on disk to the trunk's own conditions with lysine instead of glycine: glucose +
lysine, 0.15 + 0.15 mol/L, phosphate pH 7.0, 100 / 110 / 120 C, with methylglyoxal and
3-deoxyglucosone quantified by OPD derivatisation and fitted in a multiresponse scheme. It gives:

- a barrier for Amadori -> 3-deoxyglucosone (83.4 kJ/mol) and Amadori -> methylglyoxal (112.8 kJ/mol),
  the same two steps the trunk carries as `k_ama_tdg` (Ea 97.1) and `k_ama_mgo` (Ea 124.5) from
  Martins 2005, so a second laboratory's barrier for each;
- a barrier for the methylglyoxal + lysine sink (110.2 kJ/mol) and the 3-DG + lysine sink (103.9);
- the qualitative facts that 1-deoxyglucosone was **not detected** at pH 7.0 and that 3-DG ran
  higher than methylglyoxal at all three temperatures (direction only; the curves are figure-only).

What it does not give: any printed concentration of either dicarbonyl, any glyoxal at all, and a
time unit for the rate constants. The rate constants are therefore barrier-and-ratio evidence, not
transportable rates, unless the unit is resolved (Flag 1).

## 2. Methods as they matter to a model

- **Pot.** "Equal molar (0.15 mol) of both reactants, i.e. D-glucose (>= 99.0 %) and L-lysine
  (>= 99.0 %), were completely dissolved into 1 L phosphate buffer (0.1 mol/L, pH = 7.0)". So
  **[Glc] = [Lys] = 150 mmol/L, 100 mmol/L phosphate, pH 7.0 initial**. Whether lysine was the
  free base or the hydrochloride is not stated; pH after heating is not reported.
- **Heating.** 15 mL in a screw-capped glass tube (26 x 125 mm), oil bath at 100, 110 or 120 C
  (+/- 1 C); a 3 min pre-heat brings the tube to temperature and **reaction time is counted after
  the 3 min pre-heat**; "predetermined heating times ranging from 0 to 60 min"; ice-water quench;
  stored at -20 C. Headspace volume in the tube is not stated (26 x 125 mm tube holds ~ 55 mL, so
  ~ 40 mL headspace over 15 mL of liquid, my estimate).
- **Glucose**: dinitrosalicylic acid colorimetry at 540 nm (Miller 1959), standards 0.1 to 10 mg/mL
  (0.56 to 56 mmol/L). Not specific to glucose: fructose and other reducing species read as glucose.
- **Lysine**: TCA precipitation (10 %, 1:1), amino-acid analyser (Agilent 1260 Infinity II),
  ninhydrin, 570 nm.
- **3-DG and methylglyoxal** (the dicarbonyl method): 1 mL sample + 1 mL water + 2 mL of
  **1 mol/L OPD in methanol**, dark, 25 C, **24 h**; 0.45 um filter; 50 uL onto a Sunfire C18
  (4.6 x 250 mm, 5 um), 40 C, mobile phase A 20 mmol/L ammonium acetate pH 3.50, B acetonitrile,
  90 -> 70 % A over 50 min at 0.5 mL/min, PDA at 312 nm (the methods paragraph also says 210 nm;
  312 nm is the quinoxaline wavelength). Products named: 2-methylquinoxaline (from MG) and
  2-(2,3,4-trihydroxybutyl)quinoxaline (from 3-DG). **No calibration curve, LOD, LOQ, recovery or
  standard source for the quinoxalines is printed**; "quantified accordingly" refers to Yu 2018.
  Glyoxal (quinoxaline) and 1-DG were not reported; "only MG and 3-DG were detected". A 24 h OPD
  incubation of an unquenched sugar/Amadori mixture can generate dicarbonyls during derivatisation
  (a known artefact of long OPD incubations); no blank is described.
- **Furosine**: HPLC-PDA 280 nm, Sunfire C18, 30 C, A 0.4 % acetic acid, B 0.27 % KCl, 1.2 mL/min,
  external standard (Neosystem). **No acid-hydrolysis step is described**: the "furosine" is the
  280 nm peak at the standard's retention time in the unhydrolysed reaction mixture (Flag 4).
- **Pyrraline**: HPLC-PDA 297 nm, same column, 32 C, same solvents, external standard (Neosystem).
- **Melanoidins**: A420, converted to mmol/L "based on Lambert-Beer's with an extinction coefficient
  of 0.60 L mmol-1 cm-1 as reported in Kim (2010)" — a declared convention, not a measurement.
- **Model.** Scheme Fig. 1 (re-typed from Eqs 1-8):

  | step | reaction | rate law (concentrations in mmol/L) | order |
  |---|---|---|---|
  | 1 | Glc + Lys -> Amd | k1 [Glc][Lys] | 2 |
  | 2 | Amd -> 3-DG (+ Lys regenerated) | k2 [Amd] | 1 |
  | 3 | Amd -> MG (+ Lys regenerated) | k3 [Amd] | 1 |
  | 4 | 3-DG + Lys -> melanoidins | k4 [3-DG][Lys] | 2 |
  | 5 | MG + Lys -> melanoidins | k5 [MG][Lys] | 2 |
  | 6 | 3-DG + Lys -> furosine | k6 [3-DG][Lys] | 2 |
  | 7 | 3-DG + Lys -> pyrraline | k7 [3-DG][Lys] | 2 |
  | 8 | 2 MG + Lys -> furosine | k8 [MG]^2 [Lys] | 3 |
  | 9 | 2 MG + Lys -> pyrraline | k9 [MG]^2 [Lys] | 3 |

  Full balances: r_Glc = -k1[Glc][Lys]; r_Lys = -k1[Glc][Lys] + (k2 + k3)[Amd] - (k4 + k6 + k7)[3-DG][Lys]
  - k5[MG][Lys] - (k8 + k9)[MG]^2[Lys]; r_Amd = k1[Glc][Lys] - (k2 + k3)[Amd]; r_3-DG = k2[Amd] -
  (k4 + k6 + k7)[3-DG][Lys] (the printed Eq. 4 omits the "[Lys]" factor on the loss term, an
  inconsistency with Eq. 2; I read it as a typo); r_MG = k3[Amd] - k5[MG][Lys] - (k8 + k9)[MG]^2[Lys];
  r_furosine = k6[3-DG][Lys] + k8[MG]^2[Lys]; r_pyrraline = k7[3-DG][Lys] + k9[MG]^2[Lys];
  r_melanoidins = k4[3-DG][Lys] + k5[MG][Lys]. **There is no glucose isomerisation, no sugar-only
  degradation, no 1-DG, no glyoxal, no dicarbonyl sink other than lysine adducts, and no acid
  formation.** Solved with MATLAB `dsolve`, fitted with `nlinfit`, one temperature at a time.
- **Temperature dependence.** Reparameterised Arrhenius k = A exp(-Y Ea), Y = (1/R)(1/T - 1/T_av),
  T_av = sum(T)/n = **110 C = 383.15 K**, so the printed "A" is **k at 110 C**, not a
  pre-exponential factor (same convention as Yu 2018).
- **Units.** Concentrations are in mmol/L (Eq. text). So k1, k4-k7 are L mmol-1 time-1; k2, k3 are
  time-1; k8, k9 are L^2 mmol-2 time-1. **The time unit is printed nowhere**; the Methods count
  reaction time in minutes, the figure x-axes are in seconds (raster of Fig. 3: "Reaction time
  (s)", 0-4000 s). See Flag 1 for what each choice implies.
- **Replication.** Triplicate experiments; Duncan's test for the Ea/A "significance letters" (not
  reproduced in the text layer).

## 3. Tables re-typed

### Table 1. "Kinetic model parameters of Ea (kJ mol-1) and its corresponding A for steps 1-8 (in Fig. 1) in the MR model system" (nine rows are printed despite the "1-8" in the caption)

Footnote: "The unit of A is the same as that of the corresponding k." +/- is as printed (the
paper does not say whether it is SD, SE or a confidence interval).

| step | Ea (kJ/mol) | A = k(110 C) | unit of A (by rate law; time unit unresolved) |
|---|---:|---:|---|
| 1 Glc + Lys -> Amd | 89.57 +/- 4.33 | 7.68e-5 +/- 7.54e-6 | L mmol-1 time-1 |
| 2 Amd -> 3-DG | 83.37 +/- 5.86 | 5.37e-6 +/- 4.93e-6 | time-1 |
| 3 Amd -> MG | 112.82 +/- 5.29 | 3.39e-9 +/- 1.68e-9 | time-1 |
| 4 3-DG + Lys -> Mel | 103.85 +/- 1.22 | 7.68e-10 +/- 5.02e-10 | L mmol-1 time-1 |
| 5 MG + Lys -> Mel | 110.23 +/- 2.32 | 4.79e-10 +/- 1.03e-10 | L mmol-1 time-1 |
| 6 3-DG + Lys -> furosine | 81.70 +/- 14.01 | 3.54e-6 +/- 3.40e-6 | L mmol-1 time-1 |
| 7 3-DG + Lys -> pyrraline | 53.45 +/- 4.02 | 3.34e-3 +/- 3.22e-3 | L mmol-1 time-1 |
| 8 2 MG + Lys -> furosine | 52.08 +/- 4.48 | 2.27e-3 +/- 1.98e-3 | L^2 mmol-2 time-1 |
| 9 2 MG + Lys -> pyrraline | 110.22 +/- 18.77 | 4.40e-8 +/- 4.01e-8 | L^2 mmol-2 time-1 |

Text-vs-table discrepancies: the text prints step 2 as "83.38" and step 3 as "112.83" (table
83.37 / 112.82, rounding); the melanoidin paragraph quotes "4.40e-8 +/- 4.01e-8" as the exponential
factor of step 5, which is step 9's value (table: step 5 A = 4.79e-10). Table 1 is the number to keep.

**k at 100 and 120 C (mine, from A and Ea through the paper's own form, same unresolved time
unit):**

| step | k(100 C) | k(110 C) = A | k(120 C) | k120/k100 |
|---|---:|---:|---:|---:|
| 1 | 3.61e-5 | 7.68e-5 | 1.57e-4 | 4.34 |
| 2 | 2.66e-6 | 5.37e-6 | 1.04e-5 | 3.92 |
| 3 | 1.31e-9 | 3.39e-9 | 8.35e-9 | 6.36 |
| 4 | 3.21e-10 | 7.68e-10 | 1.76e-9 | 5.49 |
| 5 | 1.89e-10 | 4.79e-10 | 1.16e-9 | 6.10 |
| 6 | 1.78e-6 | 3.54e-6 | 6.80e-6 | 3.82 |
| 7 | 2.13e-3 | 3.34e-3 | 5.12e-3 | 2.40 |
| 8 | 1.46e-3 | 2.27e-3 | 3.44e-3 | 2.35 |
| 9 | 1.74e-8 | 4.40e-8 | 1.06e-7 | 6.09 |

**Unit-free ratios at 110 C (the transportable part):** k2/k3 (Amd -> 3-DG over Amd -> MG) =
5.37e-6 / 3.39e-9 = **1584**; k7/k6 (pyrraline over furosine from 3-DG + Lys) = 944; k8/k9
(furosine over pyrraline from 2 MG + Lys) = 5.2e4; k4/k5 (melanoidin from 3-DG over from MG, per
mmol/L of dicarbonyl) = 1.60.

### Table 2. "Relative RMSE and R2 for evaluating kinetic model fitting performance"

| compound | 100 C RMSE | R2 | 110 C RMSE | R2 | 120 C RMSE | R2 |
|---|---:|---:|---:|---:|---:|---:|
| D-glucose | 1.59 % | 0.98 | 2.86 % | 0.99 | 2.91 % | 0.99 |
| L-lysine | 0.62 % | 0.99 | 3.06 % | 0.98 | 3.47 % | 0.97 |
| Amadori compounds | 1.14 % | 0.99 | 4.01 % | 0.99 | 4.08 % | 0.99 |
| 3-DG | 0.89 % | 0.99 | 1.82 % | 0.97 | 0.74 % | 0.99 |
| MG | 0.29 % | 0.99 | 0.56 % | 0.99 | 1.28 % | 0.99 |
| Furosine | 0.33 % | 0.99 | 0.88 % | 0.99 | 0.69 % | 0.99 |
| Pyrraline | 0.96 % | 0.96 | 0.97 % | 0.97 | 0.68 % | 0.98 |
| Melanoidins | 0.24 % | 0.99 | 2.70 % | 0.99 | 3.12 % | 0.99 |

### Numbers in the text (the underlying curves are FIGURE-ONLY: Figs 2, 3, 4)

- "the depleted concentration of D-glucose was about 1.67 times higher than that of L-lysine at
  the same reaction time and temperature" (all T; cf. 1.7 in the xylose-lysine paper).
- "the concentration of 3-DG was normally higher than that of MG"; "a decreasing trend of 3-DG
  after a 40 min-reaction at 120 C was observed"; 1-DG absent ("may be attributed to a relatively
  low pH condition that suppressed ... 2,3-enolization").
- Furosine at 120 C "approximately two times higher than that at 100 C" at the same time; the
  pyrraline paragraph says "the concentration of furosine generated at 120 C was approximately four
  times higher than that at 100 C" — from context this second sentence is about pyrraline.
- Fig. 3 axes (raster, labels only): y "Concentration (mmol/L)" with a full scale of order 1
  mmol/L for both 3-DG and MG; x "Reaction time (s)" to 4000 s. No point was read.
- Comparators quoted by the authors: Martins 2005 step 1 Ea 96.80 +/- 2.80, X 1.6e-5 +/- 3.3e-7;
  Yu 2018 step 1 (pH 10, ultrasound) 72.9 +/- 4.1, 8.5e-5 +/- 4.8e-6; furosine Ea in foods 83.3
  (apricot), 93.9 (tomato), 88.9 (infant formula); Liang 2016 peptide-bound pyrraline from 3-DG Ea
  72.45 +/- 17.32, 63.42 +/- 20.21, 80.22 +/- 17.69 kJ/mol at peptide:glucose 1:1, 4:1, 1:4.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): furosine -> `furosine`; **3-deoxyglucosone,
methylglyoxal, glyoxal, 1-deoxyglucosone, pyrraline, the Amadori pool -> not in registry**;
melanoidins -> no molecule row (browning readout).

| step | quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|---|
| Glc + Lys -> Amadori | Ea | 89.57 +/- 4.33 | kJ/mol | 150 + 150 mmol/L, 0.1 M phosphate pH 7.0, 100-120 C | Table 1 step 1 | measured_barrier (multiresponse fit) |
| Glc + Lys -> Amadori | k(110 C) | 7.68e-5 +/- 7.54e-6 | L mmol-1 per (s or min, unresolved) | same | Table 1 ("A") | measured_rate, time unit unresolved |
| Amadori -> 3-DG | Ea | 83.37 +/- 5.86 | kJ/mol | same | Table 1 step 2 | measured_barrier (multiresponse fit); trunk `k_ama_tdg` carries 97.1 +/- 1.7 (Martins) |
| Amadori -> 3-DG | k(110 C) | 5.37e-6 +/- 4.93e-6 | per (s or min) | same | Table 1 | measured_rate, time unit unresolved; SE 92 % of value |
| Amadori -> methylglyoxal | Ea | 112.82 +/- 5.29 | kJ/mol | same | Table 1 step 3 | measured_barrier (multiresponse fit); trunk `k_ama_mgo` carries 124.5 +/- 4.7 (Martins) |
| Amadori -> methylglyoxal | k(110 C) | 3.39e-9 +/- 1.68e-9 | per (s or min) | same | Table 1 | measured_rate, time unit unresolved |
| Amadori -> 3-DG vs Amadori -> MG | k2 / k3 | 1584 | ratio (unit-free) | 110 C, same pot | derived from Table 1 | within_study_ratio (Martins 2005 at 100 C: 1.1e-2 / 7.1e-3 = 1.55) |
| 3-DG + Lys -> melanoidin | Ea | 103.85 +/- 1.22 | kJ/mol | same | Table 1 step 4 | measured_barrier (multiresponse fit); trunk `k_tdg_mel` (3-DG + Gly) carries 95.2 +/- 2.3 |
| MG + Lys -> melanoidin | Ea | 110.23 +/- 2.32 | kJ/mol | same | Table 1 step 5 | measured_barrier (multiresponse fit) — the only methylglyoxal-sink barrier in a lysine pot on disk |
| MG + Lys -> melanoidin | k(110 C) | 4.79e-10 +/- 1.03e-10 | L mmol-1 per (s or min) | same | Table 1 | measured_rate, time unit unresolved |
| 3-DG + Lys -> furosine | Ea | 81.70 +/- 14.01 | kJ/mol | same | Table 1 step 6 | measured_barrier (model-conditional; see Flag 4) |
| 3-DG + Lys -> pyrraline | Ea | 53.45 +/- 4.02 | kJ/mol | same | Table 1 step 7 | measured_barrier (multiresponse fit) |
| 2 MG + Lys -> furosine | Ea | 52.08 +/- 4.48 | kJ/mol | same | Table 1 step 8 | measured_barrier (model-conditional, third order; Flag 4) |
| 2 MG + Lys -> pyrraline | Ea | 110.22 +/- 18.77 | kJ/mol | same | Table 1 step 9 | measured_barrier (model-conditional, third order) |
| Amadori -> 1-DG | 1-deoxyglucosone level | not detected | — | pH 7.0, 100-120 C, 0-60 min | text 3.2 | level_only (verified negative in this method) |
| dicarbonyl levels | 3-DG vs MG | 3-DG > MG at every time and temperature; 3-DG falls after 40 min at 120 C | direction | same | text 3.2 | level_only |
| Glc + Lys -> Amadori | glucose loss / lysine loss | 1.67 | ratio | all runs | text 3.1 | within_study_ratio |
| all | 3-DG, MG, furosine, pyrraline, melanoidin, Glc, Lys, Amd concentrations vs time | — | mmol/L (axis) | 100/110/120 C | Figs 2-4 | figure_only |
| glyoxal | any | never measured | — | — | — | not in this paper |

Cross-reference inside the repo: the trunk's Amadori -> 3-DG and Amadori -> methylglyoxal barriers
(`src/kinetic_core/parameters.py`, `MARTINS_M4` steps 4 and 6: 97.1 and 124.5 kJ/mol, glucose +
glycine, pH 6.8) are each 12 to 14 kJ/mol above this paper's lysine-pot values (83.4 and 112.8);
the gap between the two barriers (Martins 27.4, Yu 29.5 kJ/mol) is the same in both laboratories.
That gap, not the absolute constants, is what a dicarbonyl wave can take from here.

## 5. Flags

1. **The time unit of every k is not printed** (as in Yu 2018, Flag 1 there). Concentrations are
   mmol/L. The two candidate readings give very different pots: with **time in minutes**, k1 x 150
   x 150 = 1.7 mmol/L/min of glucose at 110 C, i.e. ~ 70 % of the glucose in 60 min; with **time
   in seconds** (the Fig. 3 axis), the same product is 1.7 mmol/L/s and the glucose would be gone in
   about 2 min, which the 4000 s axis of the figure contradicts. Either the time unit is minutes
   and the figure axis is in seconds for display only, or the concentrations were fitted in mol/L.
   Until resolved, store the k values as unit-unresolved and use only the barriers and the ratios.
2. **The Amadori pool is nearly inert in the fit.** k2 + k3 = 5.4e-6 (per min or s) at 110 C is an
   Amadori half-life of 90 days (min) or 36 h (s); Martins 2005 has the DFG half-life at 100 C at
   20 min (k4 + k6 + k7 = 0.034 /min). So in this fit almost all of the reacted glucose stays as
   "Amadori compounds", and 3-DG and MG are tiny fluxes off a large pool. This is a consequence of
   the scheme (no sugar-only route, no acids, no other sinks, Amadori measured only by difference?
   — how [Amd] was measured is not stated anywhere in the Methods; it is fitted, apparently as a
   balance). The dicarbonyl constants are therefore heavily model-conditional.
3. **Standard errors of the order of the values**: step 2 (92 %), step 6 (96 %), step 7 (96 %),
   step 8 (87 %), step 9 (91 %). Only the barriers of steps 1-5 are well determined.
4. **Furosine is the wrong species for the scheme.** Furosine is formed by acid hydrolysis of
   fructosyl-lysine; here it is quantified by direct HPLC of the unhydrolysed pot with no
   hydrolysis step, and modelled as a product of 3-DG + Lys and of 2 MG + Lys (a third-order step
   the authors justify by "MG is more likely to react with another MG ... and form a furan ring").
   The steps 6, 8 barriers describe whatever the 280 nm peak is; do not carry them as furosine
   kinetics. Pyrraline from 3-DG + lysine (step 7) is the chemically expected route (Paal-Knorr).
5. **Dicarbonyl method has no printed calibration**: no standard source, curve, LOD/LOQ or recovery
   for the quinoxalines; 1 mol/L OPD in methanol for 24 h at 25 C on an unquenched mixture (a
   condition under which sugars and Amadori compounds keep releasing dicarbonyls); detection
   wavelength stated twice (210 and 312 nm). Glyoxal-quinoxaline elutes in this system but is not
   reported, so "not detected" for glyoxal cannot be claimed either way.
6. **Melanoidin concentration is a convention** (A420 / 0.60 L mmol-1 cm-1 from Kim 2010).
7. **pH not followed**; 0.1 M phosphate is a modest buffer against 150 mM lysine and the acids a
   60 min run makes. Lysine form (free base vs HCl) not stated.
8. **No glucose isomerisation, no 1-DG, no glyoxal, no acids in the scheme**; the DNS glucose assay
   counts fructose as glucose, so "glucose" is really total reducing sugar.
9. Eq. 4 as printed lacks the [Lys] factor that Eq. 2 carries on the same term; read as a typo.
10. **Supplementary data not on disk**; the concentration tables, if any exist, would be there.
11. Registry gap: 3-deoxyglucosone, methylglyoxal, glyoxal, 1-deoxyglucosone and pyrraline have no
    molecule row in `compounds.yml`.
