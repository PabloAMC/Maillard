# Gursul Aktag 2020 — EXTRACTION (laboratory-made apple juice, orange juice and peach nectar at pH 3.4, stored 4 / 27 / 37 C for 24 weeks with sampling every 2 weeks; sugars, free amino acids, five alpha-dicarbonyls and HMF followed in time; a 15-step multiresponse network fitted per juice per temperature with 95 % HPD intervals and Arrhenius barriers)

### This paper is in the plant-protein cluster by accident of the file name: it is a fruit-juice sugar-degradation kinetics paper with no protein, no lipid and no volatile in it — but it is the only member of the cluster that resolves anything chemical over a real storage time series, and the repository already holds a fuller dossier on it under a different name.

**Source on disk:** `data/articles/Gursul2020.pdf` (585,691 bytes, 11 pp., Food Chemistry 320
(2020) 126620). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Gursul2020.txt`), whole file. **Table 1 came through with its column
header damaged**: the six k/HPD column pairs are all present and legible, but the header row prints
"37 C" for the Apple Juice block and then "37 C / 27 C" for Orange and Peach, losing the Apple
27 C label; and the last two rows (k14 and k15) have their HPD values truncated to a bare "±"
(see Flags 4). The column order is recovered below from the running text and is verified against
four independently quoted values. Table 2 is clean. **Figure 1** (the three production flowcharts)
came through as readable text and is transcribed. **Figures 2 and 3** (the comprehensive and the
proposed reaction networks) and **Figure 4a-c** (the model fit over the concentration-time data for
all three juices) are images: **every concentration-time datum in this paper is figure-only**.
**Supplementary Material — Tables S1-S4 and Figures S1-S4** — is online-only and is **not on
disk**; it holds the 4 C data (S1), the free amino acids (S2), the pH and Brix (S3), the
comprehensive model's rate constants (S4) and the mass balances.

**Repo status before this dossier: this paper already has a dossier**, written on 2026-08-29 for
wave K5a, as `data/lit/extraction_dossiers/gursulaktag2020_extraction.md`. That file is longer, is
organised around a different question (whether the abstract's claim about the fructose-versus-3-DG
partition survives contact with Table 1 — it does not), and reaches the same Table 1 numbers this
dossier re-types. **The two files are about the same PDF and neither supersedes the other**; this
one exists because the file name `Gursul2020.pdf` was put into the carried-volatile reading list,
and its job is to record what the paper is and, for that list, what it is not.

## 0. Identity

| field | value |
|---|---|
| Title | "Multiresponse kinetic modelling of alpha-dicarbonyl compounds formation in fruit juices during storage" |
| Authors | Isil Gursul Aktag, Vural Gokmen (corresponding, vgokmen@hacettepe.edu.tr) — Food Quality and Safety (FoQuS) Research Group, Food Engineering Department, Hacettepe University, 06800 Beytepe, Ankara, Turkey |
| Venue | **Food Chemistry 320 (2020) 126620**. Received 26 Oct 2019, revised 17 Feb 2020, accepted 15 Mar 2020, online 17 Mar 2020 |
| DOI | **10.1016/j.foodchem.2020.126620** |
| Materials | **Apple juice** (Golden Delicious, clear), **orange juice** (Washington, cloudy), **peach nectar** (Bursa, added-sugar), all made in the laboratory from fruit bought at local markets |
| Naming | SUC sucrose; GLU glucose; FRU fructose; 1,2-ED the 1,2-enediol intermediate (**never quantified**); FFC fructofuranosyl cation (**never quantified**); 3-DG 3-deoxyglucosone; G glucosone; MGO methylglyoxal; GO glyoxal; T threosone; HMF 5-hydroxymethyl-2-furfural; P1 and P2 unmeasured product sinks |
| Existing dossier | `gursulaktag2020_extraction.md` (2026-08-29, wave K5a) — same paper, same tables, different framing |
| Funding | "This research did not receive any specific grant from funding agencies in the public, commercial ..." |

**The name check.** The task's reading list gives this as "Food Chemistry 320:126620". That is
exactly what is printed on page 1. There is no ambiguity and no companion paper: `Gursul2020.pdf`
is this article and nothing else.

## 1. Why it matters — and, for the carried-volatile programme, it does not

`tasks/roadmap_for_scientists.md` section 5d, Programme 7, is about a pea or soy isolate carrying
hexanal, 2-pentylfuran, 1-octen-3-ol and methoxypyrazines into a cook. **This paper contains none
of that.** There is no plant protein, no legume, no lipid, no lipoxygenase and no volatile analysis
of any kind. Its analytes are sucrose, glucose, fructose, twenty free amino acids, five
alpha-dicarbonyls and HMF, all by HPLC in fruit juice. Answering the cluster's three questions
directly:

- **(i) Does it print a concentration of a carried volatile in a named plant material?** No. It
  prints no volatile at all. HMF is the closest thing to an aroma compound in it, and HMF is a
  browning marker measured by diode-array absorbance at 285 nm, not an odorant here — it carries no
  odour threshold and no sensory measurement in this paper.
- **(ii) Does it print anything with TIME?** **Yes, and this is the paper's whole point.** Three
  juices x three temperatures x 24 weeks with sampling every two weeks in triplicate, fitted as a
  15-step multiresponse network with 95 % HPD intervals and Arrhenius activation energies. It is by
  a wide margin the most time-resolved paper in the cluster. But the concentration-time data
  themselves are drawn in Figure 4 and are **figure-only**; what is printed is the fitted constants.
- **(iii) Is the material a flour, an isolate, a concentrate or a whole seed?** **None of these.**
  It is a **liquid fruit juice or nectar**: clear apple juice (enzyme-treated and clarified), cloudy
  orange juice, and a peach nectar reconstituted from pulp with a 66 % sugar syrup and citric acid.
  For the roadmap's levels table, which needs the material named exactly, the honest entry is "not
  a plant-protein material".

What it **does** give the repository is real and already recorded elsewhere: a second laboratory's
multiresponse network for sugar degradation at **food temperatures and food pH** (3.4), on the
**storage** timescale rather than the cooking timescale, with a full set of rate constants and
barriers. `k5a_hmf_synthesis.md` and the K5a wave already use it, and
`gursulaktag2020_extraction.md` carries that analysis. The one thing this dossier adds to the
carried-volatile programme is a negative entry with a reason attached.

## 2. Methods as they matter to a model

- **The materials, exactly as described.** "Golden delicious variety of apples, Washington variety
  of oranges and Bursa variety of peaches were obtained from local markets. Apple juices, orange
  juices and peach nectars were **produced in the laboratory** using the flowchart shown in Fig. 1."
  The three were chosen "as typical examples of **clear, cloudy and added-sugar** products".
- **Production, from Figure 1 (transcribed from the text layer).**
  - *Apple*: washing -> grinding -> **heating 85 C x 5 min** -> pressing -> **enzymatic treatment at
    50 C (1 mL/L pectinase, 0.2 mL/L amylase)** -> clarification at 50 C (gelatine-bentonite
    flocculation) -> filtration -> filling into falcon tubes -> **pasteurisation in a water bath,
    90 C x 10 min**.
  - *Orange*: washing -> grinding -> pulper -> filling into falcon tubes -> **pasteurisation 90 C x
    10 min**. No heating step before the pulper, no enzyme, no clarification.
  - *Peach*: washing -> grinding -> **heating 85 C x 5 min** -> pulper -> **nectar production
    (66 % sugar syrup, citric acid, water, pulp)** -> filling into falcon tubes -> **pasteurisation
    90 C x 10 min**.
- **The storage design.** "The juice samples were stored at **4, 27 and 37 C for 24 weeks**.
  Sub-samples were taken from the stored samples **3 parallel in every 2 weeks**, and kept frozen at
  **-18 C** prior to analysis." So thirteen time points (0 to 24 weeks) x three temperatures x three
  juices x three replicates. **The 4 C arm showed no significant change in any reactant or product
  (Table S1, not on disk) and was excluded from the kinetic analysis**, so every fitted constant
  rests on 27 and 37 C only — two temperatures for an Arrhenius fit (Flags 3).
- **Sample clean-up.** Carrez clarification (1 mL of juice + 50 uL Carrez I + 50 uL Carrez II) for
  HMF and sugars; acetonitrile clarification (500 uL juice + 500 uL acetonitrile) for the
  alpha-dicarbonyls and free amino acids. Both centrifuged 10,000 g for 5 min.
- **Sugars.** Supernatant through a preconditioned OASIS HLB cartridge (first 8 drops discarded);
  Agilent 1200 HPLC with refractive index detection; **Shodex Sugar SH-1011 (300 x 8 mm, 6 um) at
  50 C**; 5 mM H2SO4 at 1 mL/min; 10 uL injected. **External calibration curves for each of sucrose,
  glucose and fructose over 0.25-2.5 g/L (five points).** LOD 1.0 mg/L, LOQ 3.0 mg/L.
- **HMF.** Filtered through 0.45 um; Shimadzu UFLC with diode-array detection; **Atlantis dC18**;
  isocratic 10 mM aqueous formic acid / acetonitrile 90:10 at 1.0 mL/min, 25 C; **285 nm**.
  **External calibration 1-10 mg/L (four points).** "The LOD and LOQ values for HMF were **10 mg/L**
  and **30 ug/L**" — as printed, which is internally impossible (Flags 5).
- **Alpha-dicarbonyls: derivatised, and only two have their own standard.** 500 uL supernatant +
  150 uL of 0.2 % **o-phenylenediamine** containing 11 mM DETAPAC + 150 uL of 0.5 M sodium phosphate
  buffer pH 7; filtered; held **2 h at room temperature in the dark**. Agilent 1200 HPLC + Agilent
  6130 single quadrupole MS; **Merck Purospher Star RP-18e (150 x 4.6 mm, 5 um)**; gradient of 1 %
  formic acid in water and 1 % formic acid in methanol, 0.7 mL/min, 30 C, 30 % B rising to 60 % B in
  12 min then back to 30 % B in 3 min, 15 min run; 10 uL injected. Electrospray positive, drying gas
  11.0 L/min at 320 C, nebuliser 400 psig, capillary 4000 V, fragmentor 100 V. **SIM [M+H]+**:
  glucosone quinoxaline 251; 1- or 3-DG 235; MGO 145; GO 131; **threosone 191**; dwell 97 ms.
  **Authentic standards and response factors:** quinoxaline, 2-methylquinoxaline and
  2,3-dimethylquinoxaline (the derivatives of GO, MGO and 2,3-butanedione respectively) each have an
  **external calibration curve over 0.1-1.0 mg/L (five points)**; 3-DG (75 % purity) and glucosone
  (>= 98 %) were derivatised and calibrated the same way. **Threosone has no standard**: "the
  calibration curve of glucosone was used for **semi-quantitation** of threosone derivatives, since
  both have same proton-accepting groups." LOD 2.5-15 ug/L, LOQ 8.3-50 ug/L.
- **Free amino acids.** Waters Acquity UPLC + triple quadrupole, positive electrospray; **Atlantis
  HILIC (150 x 2.1 mm, 3 um) at 30 C**; gradient of 0.1 % formic acid in water and in acetonitrile
  at 0.4 mL/min; autosampler 10 C; capillary 3.5 kV, cone 20 V, extractor 3 V, source 120 C,
  desolvation 350 C at 900 L/h. **External calibration for all amino acids over 0.1-5.0 mg/L (six
  points).** LOD 0.2-5 ug/L, LOQ 0.7-16.7 ug/L.
- **pH and Brix.** PHM210 pH meter; Atago Pocket Pal-3 refractometer. Values in Table S3, **not on
  disk**; the text gives initial Brix as **15, 12 and 16** for apple, orange and peach, and states
  the juices are at **pH 3.4** and that pH did not change during storage.
- **Basis of every concentration.** "**Data used for modelling was expressed in mmol/L.**" Sugars
  and dicarbonyls in mmol/L of juice; HMF quoted in the text in mg/L and modelled in mmol/L; amino
  acids in mg/L. All rate constants in **week^-1**, printed in Table 1 as **week^-1 x 10^3**.
- **Fitting.** "All individual analytical measurements were used to estimate model parameters"
  (not the means). Differential equations from the network of Fig. 3, numerically integrated;
  **non-linear regression with the determinant criterion (van Boekel 1996)** in **Athena Visual
  Studio v14.2**. **Parameter estimation was performed separately for each storage temperature** —
  so the constants are not a simultaneous multi-temperature fit, and the Arrhenius barriers of
  Table 2 are computed afterwards from two points. Model discrimination by goodness of fit and by
  the HPD intervals. **All steps first order** ("Since it is known that sugar degradation reactions
  obey first degree reaction kinetics, these elementary reaction steps in the proposed model were
  defined by differential equations in accordance with first degree kinetics").
- **Statistics.** Verbatim and complete: "**Free amino acid data** were subjected to analysis of
  variance (one-way ANOVA). The SPSS 17.0 statistical package was used for the evaluation of
  statistical significance of the differences between mean values by Tukey test. P < 0.05 was
  considered to be statistically significant for the results." **No significance test on any rate
  constant is described anywhere** (Flags 6).

## 3. Tables re-typed

### Table 1. "Estimated reaction rate constants (k, week^-1 x 10^3) with 95 % highest posterior density (HPD) intervals according to the proposed kinetic model shown in Fig. 3 for sugar degradation reactions in apple juice, orange juice and peach nectar during storage at different temperatures"

Footnote as printed: "*ind, indeterminate, which means a large uncertainty in the estimated
parameter within 95 % confidence interval."

**Column order.** The printed header reads "Apple Juice | Orange Juice | Peach Nectar" over six
k/HPD pairs, but only prints "37 C" once for Apple and then "37 C / 27 C" for the other two,
so the Apple 27 C label is missing from the text layer. The order below is **Apple 37, Apple 27,
Orange 37, Orange 27, Peach 37, Peach 27**, and it is verified against four values the running text
quotes independently: k1 at 27 C is quoted as "0.04, 0.05 and 0.05 week^-1 ... for apple juice,
orange juice and peach nectar" against columns 2, 4 and 6 (36.3, 50.7, 44.8 x 10^-3 — 0.036, 0.051,
0.045); k6 at 37 C is quoted as "6.7 x 10^-3, 0.7 x 10^-3, and 0.3 x 10^-3 for orange juice, apple
juice, and peach nectar" against columns 3, 1 and 5 (6.7, 0.7, 0.3); k12 at 37 C is quoted as
"0.010, 0.009 and 0.023 week^-1 for apple juice, orange juice, and peach nectar" against columns 1,
3 and 5 (10.0, 9.2, 23.2); and k13 in peach is quoted as "0.022 week^-1 and 0.029 week^-1 at 37 and
27 C" against columns 5 and 6 (22.1, 29.4). **Four independent confirmations; the ordering is
secure.**

| step | reaction | Apple 37 C, k | HPD | Apple 27 C, k | HPD | Orange 37 C, k | HPD | Orange 27 C, k | HPD | Peach 37 C, k | HPD | Peach 27 C, k | HPD |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SUC -> FUR + GLU | 123.7 | ± 10.03 | 36.3 | ± 3.06 | 147.4 | ± 9.68 | 50.7 | ± 3.29 | 147.2 | ± 13.1 | 44.8 | ± 3.79 |
| 2 | GLU -> 1,2-ED | 598.4 | **ind\*** | 33.1 | ± 22.02 | 58.3 | ± 29.84 | 37.5 | ± 8.48 | 3690.7 | **ind\*** | 595.5 | ± 556.3 |
| 3 | 1,2-ED -> GLU | 781.8 | ± 111.70 | 66.6 | ± 104.20 | 495.4 | ± 277.1 | 99.3 | ± 29.90 | 2047.3 | ± 161.5 | 441.6 | ± 531.4 |
| 4 | 1,2-ED -> FRU | 544.3 | ± 390.00 | 90.5 | ± 83.33 | 28.5 | ± 30.53 | 39.2 | ± 34.08 | 407.8 | ± 282 | 65.0 | ± 68.81 |
| 5 | FRU -> 1,2-ED | 143.1 | ± 98.33 | 0.0 | ± 0.0 | 0.9 | **fixed** | 17.8 | ± 9.34 | 684.2 | ± 463.6 | 105.0 | ± 66.6 |
| 6 | 1,2-ED -> 3-DG | 0.7 | ± 0.08 | 2.1 | ± 1.22 | 6.7 | ± 2.14 | 0.4 | ± 0.05 | 0.3 | ± 0.037 | 36.4 | ± 33.03 |
| 7 | 1,2-ED -> G | 2.0 | ± 0.27 | 9.6 | ± 5.73 | 20.4 | ± 6.59 | 1.2 | ± 0.14 | 0.0001 | ± 0.002 | 0.00002 | ± 0.009 |
| 8 | 3-DG -> HMF | 0.7 | ± 1.98 | 0.0 | **fixed** | 2.0 | ± 0.68 | 0.0 | **fixed** | 11.8 | ± 2.79 | 0.0 | **fixed** |
| 9 | 3-DG -> MGO | 4.6 | ± 0.41 | 3.5 | ± 0.43 | 6.0 | ± 0.60 | 11.7 | ± 2.51 | 0.0 | **fixed** | 0.0 | **fixed** |
| 10 | FRU -> FFC | 0.2 | ± 0.04 | 0.0 | ± 0 | 7.1 | ± 2.26 | 0.0 | ± 0.0 | 22.6 | ± 12.69 | 0.0 | **fixed** |
| 11 | FFC -> HMF | 5.3 | **ind\*** | 0.0 | **fixed** | 0.0 | ± 0.02 | 0.0 | **fixed** | 0.4 | ± 0.22 | 0.0 | **fixed** |
| 12 | G -> GO | 10.0 | ± 1.02 | 5.1 | ± 1.36 | 9.2 | ± 0.74 | 8.7 | ± 1.49 | 23.2 | ± 1.57 | 4.8 | ± 1.25 |
| 13 | G -> T | 0.0 | **fixed** | 0.0 | **fixed** | 0.0 | **fixed** | 0.0 | **fixed** | 22.1 | ± 1.19 | 29.4 | ± 1 |
| 14 | G -> P1 | 56.1 | **fixed** | 562.0 | **fixed** | 453.2 | **fixed** | 21.1 | **fixed** | 0.0 | **fixed** | 29.7 | ± *(truncated)* |
| 15 | 3-DG -> P2 | 107.6 | **fixed** | 123.3 | **fixed** | 228.2 | **fixed** | 71.7 | **fixed** | 427.2 | **fixed** | 94365.5 | ± *(truncated)* |

The HPD values for the Peach 27 C entries of k14 and k15 are printed as a bare "±" with nothing
after it in the text layer. **They are not transcribed and are not guessed** (Flags 4).

### Table 2. "Estimated activation energies (Ea, kJ/mol) according to the proposed kinetic model shown in Fig. 3 for sugar degradation reactions in apple juice, orange juice and peach nectar during storage at different temperatures"

| step | reaction | Apple Ea | Apple R2 | Orange Ea | Orange R2 | Peach Ea | Peach R2 |
|---|---|---|---|---|---|---|---|
| 1 | SUC -> FUR + GLU | **54.16** | 0.96 | **55.51** | 0.98 | **48.23** | 0.97 |
| 2 | GLU -> 1,2-ED | 11.10 | 0.03 | 48.19 | 0.99 | 34.95 | 0.43 |
| 3 | 1,2-ED -> GLU | -5.94 | 0.01 | 56.29 | 0.90 | 79.16 | 0.86 |
| 4 | 1,2-ED -> FRU | 2.41 | 0.00 | -81.09 | 0.96 | -18.61 | 0.51 |
| 5 | FRU -> 1,2-ED | 65.09 | 0.53 | -145.77 | 0.97 | 19.99 | 0.18 |
| 6 | 1,2-ED -> 3-DG | -20.02 | 0.52 | **112.67** | 0.94 | -141.35 | 0.90 |
| 7 | 1,2-ED -> G | 10.12 | 0.07 | 94.62 | 0.88 | 23.43 | 0.72 |
| 8 | 3-DG -> HMF | -127.02 | 0.51 | -107.71 | 0.51 | -77.04 | 0.51 |
| 9 | 3-DG -> MGO | 18.43 | 1.00 | -24.21 | 0.90 | – |  |
| 10 | FRU -> FFC | -30.71 | 0.03 | 65.16 | 0.21 | -65.77 | 0.51 |
| 11 | FFC -> HMF | -90.94 | 0.51 | **-230.31** | 0.51 | -188.99 | 0.51 |
| 12 | G -> GO | -9.82 | 0.72 | -0.15 | 0.02 | -97.29 | 0.71 |
| 13 | G -> T | – |  | – |  | 26.54 | 0.79 |
| 14 | G -> P1 | -36.53 | 0.46 | 15.13 | 0.06 | 54.27 | 0.34 |
| 15 | 3-DG -> P2 | -26.52 | 0.97 | -41.38 | 0.56 | **121.55** | 0.41 |

**Twenty-one of the forty-three printed activation energies are negative** (mine, counted). The
authors address this: "The calculated negative activation energy values might indicate that no
energy barriers were presented in these reaction steps due to the accumulation of intermediate
compounds (Vyazovkin, 2016)." The text states the range as "-230 and 122 kJ/mol", which matches
Table 2's extremes (-230.31 for orange FFC -> HMF, 121.55 for peach 3-DG -> P2).

### Appendix A — the differential equations, as printed

```
d[SUC]/dt      = -k1[SUC]
d[GLU]/dt      =  k1[SUC] + k3[1,2-enediol] - k2[GLU]
d[FRU]/dt      =  k1[SUC] + k4[1,2-enediol] - (k5 + k10)[FRU]
d[HMF]/dt      =  k11[FFC] + k8[3-DG]
d[3-DG]/dt     =  k6[1,2-enediol] - (k8 + k9 + k15)[3-DG]
d[G]/dt        =  k7[1,2-enediol] - (k12 + k13 + k14)[G]
d[GO]/dt       =  k12[G]
d[T]/dt        =  k13[G]
d[1,2-enediol]/dt = k2[GLU] + k5[FRU] - (k3 + k4 + k6 + k7)[1,2-enediol]
d[FFC]/dt      =  k10[FRU] - k11[FFC]
d[P1]/dt       =  k14[G]
d[P2]/dt       =  k15[3-DG]
```

Note that `d[MGO]/dt` is **not printed in Appendix A**, although k9 (3-DG -> MGO) is in the network
and in Table 1.

### Numbers printed only in the running text

| quantity | value | where |
|---|---|---|
| pH of the juices | **3.4**; "no change was observed in the pH values of samples during storage" | section 3.1 |
| initial Brix | apple **15**, orange **12**, peach **16** | section 3.1 (values in Table S3, not on disk) |
| 4 C arm | "there were no significant changes in concentrations of reactants and reaction products in the samples stored at 4 C ... Therefore, kinetic analysis was limited with the data observed for 27 and 37 C" | section 3.1 |
| free amino acids | "remained relatively stable in all samples during storage at all temperatures" | section 3.1 (values in Table S2, not on disk) |
| sucrose loss after 24 weeks | **50 %, 62 % and 54 % at 27 C**; **93 %, 89 % and 96 % at 37 C** for apple, orange, peach | section 3.1 |
| maximum HMF | **16.2 +/- 0.7, 3.8 +/- 0.2 and 12.2 +/- 0.5 mg/L** in apple juice, orange juice and peach nectar | section 3.1 |
| the regulatory comparator | AIJN maximum HMF for fruit juices, **10 mg/L** (AIJN 1996) — apple and peach exceed it | section 3.1 |
| HMF shape | "formation of HMF followed a typical kinetic pattern in juices stored at 37 C, while there was **no accumulation of HMF at 27 C**" | section 3.1 |
| dominant dicarbonyl | **glucosone** in apple and orange; **3-DG** in peach nectar | section 3.1 |
| maximum glucosone | **3.5 mmol/L in apple juice after 16 weeks**; **1.7 mmol/L in orange juice after 14 weeks**, both at 37 C | section 3.1 |
| maximum 3-DG | **0.4 mmol/L in peach nectar after 14 weeks at 37 C**; 3-DG "started to increase dramatically with time after 4 weeks of storage at 37 C and 8 weeks at 27 C" | section 3.1 |
| maximum GO at 37 C, 24 weeks | **0.4 mmol/L apple, 0.2 mmol/L orange, 0.08 mmol/L peach** | section 3.1 |
| which short-chain dicarbonyls appear where | **GO in all three**; **MGO only in apple and orange**; **threosone only in peach nectar** | section 3.1 |
| mannose | looked for, "could not be detected", so glucose-mannose epimerisation was left out of the model | section 3.2(i) |
| 1-deoxyglucosone | "not detected in the acidic samples", so that pathway was left out | section 3.2(iii) |
| k1, sucrose hydrolysis | "0.04, 0.05 and 0.05 week^-1 at 27 C; **0.1, 0.2 and 0.2 week^-1 at 37 C** for apple, orange and peach" — the 37 C figures do not round from Table 1 (Flags 5) | section 3.2(i) |
| k6, 3-DG formation at 37 C | 6.7e-3 (orange), 0.7e-3 (apple), 0.3e-3 (peach) week^-1 | section 3.2(iii) |
| k7, glucosone formation | "0.020 week^-1, 0.002 week^-1, and **0.03 x 10^-3 week^-1** for orange juice, apple juice, and peach nectar" — the peach figure does not match Table 1 (Flags 5) | section 3.2(iii) |
| k12, GO formation at 37 C | 0.010, 0.009 and 0.023 week^-1 for apple, orange, peach | section 3.2(iii) |
| k13, threosone formation, peach only | **0.022 week^-1 at 37 C and 0.029 week^-1 at 27 C** | section 3.2(iii) |
| the HMF-partition claim | "HMF formation from fructose were found almost **7, 4 and 2 times higher** in apple juice, orange juice and peach nectar" | section 3.2(ii) |
| glucose content contrast | orange juice held "approximately **2.4 times and 3.2 times more glucose** than apple juice and peach nectar" | section 3.2(iii) |
| Arrhenius range | "the activation energies for each reaction step was found to be in the range of **-230 and 122 kJ/mol**" | section 3.2 |
| sucrose hydrolysis barrier | "fairly temperature dependent in juices (Ea; **48-56 kJ/mol**)" | section 3.2 |
| mass balance | "the recovery values were approximately 100 % or slightly higher through the storage period for all samples" | section 3.2 (Fig. S1, not on disk) |
| the comprehensive model | tried first, but "the confidence intervals of rate constants were not well estimated for comprehensive model (Table S4 and Fig. S2-S4)" — hence the simplification to Fig. 3 | section 3.2 |
| the enolisation test | removing the unquantified 1,2-enediol made "the reaction rate of each step, especially those of sucrose hydrolysis" unestimable — so the enediol stays in as an unmeasured lump | section 3.2 |
| the MGO test | removing MGO also degraded the fit, so MGO stays in "imprecisely" | section 3.2(iii) |

**Concentration-time data: FIGURE-ONLY.** Figure 4a (apple), 4b (orange) and 4c (peach) carry
sucrose, glucose, fructose, 3-DG, MGO, HMF, glucosone and GO against week, observed and predicted,
at 27 and 37 C. Per house rule they are not typed as numbers; only the maxima the text quotes are
recorded above.

### Arithmetic on the printed numbers (all mine)

1. **The Arrhenius fits have two points each.** The 4 C arm was excluded, so every Ea in Table 2 is
   computed from k(27 C) and k(37 C) alone. **A two-point Arrhenius over a 10 K interval has no
   residual degrees of freedom**, which is why so many R^2 values in Table 2 are exactly 0.51 or
   0.00 — those cannot be regression R^2 over two points and must come from something else the
   paper does not describe (Flags 7). Twenty-one of forty-three barriers are negative.
2. **The 3-DG-versus-fructose partition at the HMF-forming step, from Table 1 at 37 C.** k8
   (3-DG -> HMF) against k11 (FFC -> HMF): apple **0.7 vs 5.3**, but k8's HPD is ±1.98 (2.8x its own
   estimate) and k11 is author-declared indeterminate; orange **2.0 ± 0.68 vs 0.0 ± 0.02**; peach
   **11.8 ± 2.79 vs 0.4 ± 0.22**. **In two of the three juices the 3-DG limb is the larger one, and
   in orange the fructose limb is exactly zero.** The abstract's claim that the fructose route is
   "significantly higher (p < 0.05)" than the 3-DG route is not supported by this table and no
   significance test on a rate constant is described in the Methods. This is the central finding of
   `gursulaktag2020_extraction.md` and it reproduces here from the primary table.
3. **The "7, 4 and 2 times higher" claim mixes two different steps.** k11/k8 at 37 C gives apple
   5.3/0.7 = **7.6**, orange 0.0/2.0 = **0**, peach 0.4/11.8 = **0.034**. Only the apple figure
   matches the claimed "7"; the orange "4" and peach "2" do not come from k11/k8, and I could not
   reproduce them from any ratio in Table 1.
4. **Sucrose hydrolysis is the one clean result.** k1 is the only step whose HPD is small in all six
   columns (2.6 % to 8.9 % of the estimate, mine), whose barrier is consistent across all three
   juices (54.16, 55.51, 48.23 kJ/mol) and whose R^2 is high in all three (0.96-0.98). **If one
   number is taken from this paper, it should be this one.** Its 37/27 C rate ratio is 3.41
   (apple), 2.91 (orange) and 3.29 (peach), i.e. a Q10 of about 3.2 (mine).
5. **Half-life of sucrose from k1 (mine, first order).** At 37 C: ln2/0.1237 = **5.6 weeks**
   (apple), 4.7 weeks (orange), 4.7 weeks (peach). At 27 C: 19.1, 13.7 and 15.5 weeks. Against the
   measured 24-week losses of 93 / 89 / 96 % at 37 C, first order predicts 1 - exp(-0.1237 x 24) =
   **95 %**, 97 % and 97 % — apple and orange agree well, peach's measured 96 % agrees, so the
   first-order sucrose step is self-consistent.
6. **Some constants are enormous and unphysical as printed.** Peach 27 C k15 (3-DG -> P2) is
   **94,365.5 x 10^-3 week^-1 = 94.4 week^-1**, a half-life of 7 minutes for 3-DG in a refrigerated
   nectar, against a peach 37 C value of 0.427 week^-1 — a 221-fold *decrease* with a 10 K
   *increase* in temperature. Peach 27 C k2 (595.5) exceeds peach 37 C k2 only by being smaller
   (3690.7), but peach 27 C k6 (36.4) is **121 times larger** than peach 37 C k6 (0.3). **Several
   steps run backwards in temperature**, which is what the negative barriers in Table 2 are
   recording. These constants should not be transported anywhere.

## 4. Numbers the repository can use

All rows: laboratory-made juice at **pH 3.4**, stored in falcon tubes after pasteurisation at
90 C x 10 min, sampled every 2 weeks in triplicate over 24 weeks, frozen at -18 C before analysis;
constants fitted per juice per temperature in Athena Visual Studio by the determinant criterion,
all steps first order, unit **week^-1 x 10^3** as printed.

| quantity | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| **any carried volatile (hexanal, 2-pentylfuran, 1-octen-3-ol, a methoxypyrazine)** | **not measured** | — | — | — | **absent — no volatile analysis of any kind in this paper** |
| the material | apple juice (clear, enzyme-treated, clarified), orange juice (cloudy), peach nectar (66 % sugar syrup, citric acid, water, pulp) | — | Golden Delicious / Washington / Bursa | Fig. 1, section 2.2 | **not a flour, an isolate, a concentrate or a whole seed** |
| the fifteen rate constants, six columns each | see Table 1 above | week^-1 x 10^3 | pH 3.4, 27 and 37 C, 24 weeks | Table 1 p. 8 | measured_rate (with the caveats in Flags 2, 3 and 7) |
| **k1, sucrose hydrolysis** | 123.7 ± 10.03 / 36.3 ± 3.06 (apple 37 / 27 C); 147.4 ± 9.68 / 50.7 ± 3.29 (orange); 147.2 ± 13.1 / 44.8 ± 3.79 (peach) | week^-1 x 10^3, first order in sucrose | as above | Table 1 | measured_rate — **the best-determined constant in the paper** |
| Ea(k1) | 54.16 (apple, R2 0.96); 55.51 (orange, 0.98); 48.23 (peach, 0.97) | kJ/mol | two-point Arrhenius over 27-37 C | Table 2 p. 8 | measured_barrier (**two temperatures only**) |
| the fifteen activation energies, three juices | see Table 2 above | kJ/mol | as above | Table 2 | measured_barrier — **21 of 43 are negative** (Flags 7) |
| sucrose loss at 24 weeks | 50 / 62 / 54 % at 27 C; 93 / 89 / 96 % at 37 C | % of initial | apple / orange / peach | text, section 3.1 | measured_level |
| maximum HMF | 16.2 ± 0.7 / 3.8 ± 0.2 / 12.2 ± 0.5 | mg/L of juice | apple / orange / peach, 37 C over 24 weeks | text, section 3.1 | measured_level |
| maximum glucosone | 3.5 (apple, 16 weeks) and 1.7 (orange, 14 weeks) | mmol/L of juice | 37 C | text, section 3.1 | measured_level |
| maximum 3-DG | 0.4 | mmol/L of juice | peach nectar, 37 C, 14 weeks | text, section 3.1 | measured_level |
| maximum GO | 0.4 / 0.2 / 0.08 | mmol/L of juice | apple / orange / peach, 37 C, 24 weeks | text, section 3.1 | measured_level |
| MGO, threosone | present only in apple and orange (MGO) and only in peach (threosone); "of only minor importance" quantitatively | — | 27 and 37 C | text, section 3.1 | level_only |
| initial Brix | 15 / 12 / 16 | degrees Brix | apple / orange / peach | text, section 3.1 (Table S3 **not on disk**) | measured_level |
| pH | 3.4, unchanged over 24 weeks | — | all three juices | text, section 3.1 (Table S3 **not on disk**) | measured_level |
| free amino acids | "no statistically significant change during storage" | mg/L | all three, all temperatures | Table S2 (**not on disk**) | level_only |
| the 4 C arm | no significant change in any reactant or product over 24 weeks | — | all three juices | Table S1 (**not on disk**) | level_only |
| Q10 of sucrose hydrolysis | 3.41 / 2.91 / 3.29 | dimensionless, 27 -> 37 C | apple / orange / peach | derived from Table 1 (mine) | within_study_ratio |
| sucrose half-life | 5.6 / 4.7 / 4.7 weeks at 37 C; 19.1 / 13.7 / 15.5 weeks at 27 C | weeks, first order | apple / orange / peach | derived (mine) | derived_assumption |
| k11/k8 (fructose vs 3-DG route to HMF) at 37 C | 7.6 (apple, both terms unsafe) / 0 (orange) / 0.034 (peach) | dimensionless | as above | derived (mine) | within_study_ratio — **contradicts the abstract** (Flags 6) |
| every concentration-time course | — | mmol/L vs week | three juices, two temperatures, thirteen time points | Figs. 4a, 4b, 4c | **figure_only** |
| the comprehensive (rejected) model's constants | — | — | — | Table S4 (**not on disk**) | not available |
| the reaction networks | comprehensive (Fig. 2) and proposed (Fig. 3) | — | — | Figs. 2, 3 | **figure_only** (structures; the differential equations of Appendix A carry the same information and are transcribed above) |

### What this can and cannot be used for

**Can:** be a second laboratory's storage-timescale sugar-degradation network at food pH, which is
what wave K5a already uses it for; supply a well-determined sucrose hydrolysis constant and barrier
in three real acidic matrices; supply measured HMF, glucosone, 3-DG and GO maxima with the week at
which they occur; and supply a documented case of a multiresponse fit whose barriers come out
negative, which is a useful cautionary example for any wave doing the same thing on two
temperatures.

**Cannot:** contribute anything to Programme 7's carried-volatile levels table. Cannot supply an
Arrhenius barrier that should be trusted outside k1. Cannot supply a benchmark row with
concentrations, because every concentration-time datum is in a figure.

## 5. Flags

1. **This paper is in the wrong cluster.** It has no plant protein, no legume, no lipid, no
   lipoxygenase and no volatile. Its presence in a carried-volatile reading list is a file-name
   coincidence. Recorded so the PDF is not re-opened for that purpose. Its real home in the
   repository is the HMF and alpha-dicarbonyl work, where `gursulaktag2020_extraction.md` already
   covers it in more depth than this dossier does.
2. **A duplicate dossier already exists.** `data/lit/extraction_dossiers/gursulaktag2020_extraction.md`
   (2026-08-29, wave K5a) is about the same PDF, re-types the same Table 1, and reaches the same
   conclusion about the abstract's fructose-versus-3-DG claim. **Two dossiers on one paper is a
   maintenance hazard**: if either is corrected, the other goes stale silently. The owner should
   decide whether to merge them or to make one a pointer to the other. This file does not delete or
   modify the existing one.
3. **The kinetics rest on two temperatures.** The 4 C arm produced no measurable change and was
   dropped, leaving 27 and 37 C. Every activation energy in Table 2 is therefore a two-point slope
   over a 10 K interval, with **no residual degrees of freedom and no possibility of detecting
   curvature or a bad point**. Combined with the fact that each temperature was fitted separately
   (unlike Knol 2005 or the Leuven lane, which fit all temperatures simultaneously), this is the
   weakest possible Arrhenius design.
4. **Table 1's text layer is damaged in two places.** The header loses the Apple Juice 27 C label,
   and the HPD entries for the peach 27 C values of k14 (29.7) and k15 (94,365.5) are printed as a
   bare "±" with nothing following. The column order was recovered from four values the running text
   quotes independently and is secure; **the two truncated HPD values are not transcribed and must
   be read from the page image before anyone uses those rows.**
5. **Three printed inconsistencies.** (i) "The LOD and LOQ values for HMF were **10 mg/L** and
   **30 ug/L**, respectively" — an LOD a thousand times above its own LOQ, and in a different unit;
   the LOD is almost certainly 10 ug/L. (ii) The text gives k1 at 37 C as "0.1, 0.2 and 0.2
   week^-1"; Table 1 gives 0.124, 0.147 and 0.147, all of which round to 0.1. (iii) The text gives
   k7 for peach nectar as **0.03 x 10^-3 week^-1**; Table 1 gives **0.0001 x 10^-3** at 37 C and
   **0.00002 x 10^-3** at 27 C — a mismatch of 300x or 1500x. The other quoted values (k6, k12, k13
   and k1 at 27 C) all reconcile.
6. **The abstract's headline claim is not supported by the paper's own table, and no test exists to
   support it.** The abstract says "The contribution of fructose dehydration through fructofuranosyl
   cation on the formation of 5-hydroxymethylfurfural was significantly higher (**p < 0.05**) than
   3-deoxyglucosone dehydration." The Methods describe an ANOVA on **the free amino acid data only**;
   no significance test on any rate constant is described anywhere. And at 37 C, k8 (3-DG -> HMF) is
   **larger** than k11 (FFC -> HMF) in orange (2.0 vs 0.0) and in peach (11.8 vs 0.4); only apple
   favours the fructose route, and there k11 is author-declared indeterminate while k8's HPD is 2.8
   times its own estimate. **Do not carry the abstract's claim forward.** This is the finding of
   `gursulaktag2020_extraction.md` section 1 and it is confirmed here from the primary table.
7. **Twenty-one of forty-three activation energies are negative and several constants fall with
   rising temperature.** The authors offer an explanation (accumulation of intermediates, no energy
   barrier), but a negative Arrhenius slope from a two-point fit on an unquantified intermediate is
   more simply read as the fit being under-determined. The R^2 column reinforces this: many entries
   are exactly **0.51** or **0.00**, values that cannot arise as a regression R^2 over two points,
   and the paper never says what the R^2 refers to. **No barrier from Table 2 except k1's should
   enter a registry.**
8. **Three parameters have no independent measurement behind them.** The **1,2-enediol** and the
   **fructofuranosyl cation** are never quantified — they are lumps whose concentrations exist only
   inside the model — and **P1 and P2 are unmeasured sinks**, so k14 and k15 are fitted against
   nothing and are marked "fixed" in every column but one. **Threosone has no authentic standard**
   and is semi-quantified against the glucosone curve on the argument that the two derivatives share
   proton-accepting groups; every threosone number is therefore an estimate, and k13 inherits that.
9. **What this paper does NOT contain, and what to request.** No volatile, no protein, no lipid, no
   water activity (the juices are aqueous and it is never stated), no oxygen or headspace
   description in the falcon tubes — which matters because the authors attribute glucosone to
   **oxidation** and then note there is "no evidence of the relationship between oxygen levels in
   fruit juice and the formation of reactive carbonyl species"; no tabulated concentration; and no
   simultaneous multi-temperature fit. **To request:** (a) the **Supplementary Material** — Tables
   S1 (the 4 C arm), S2 (free amino acids), S3 (pH and Brix) and S4 (the comprehensive model), and
   Figures S1-S4 — none of which is on disk and which together hold most of the paper's raw
   numbers; (b) the numeric data behind Figure 4 (3 juices x 2 temperatures x 13 weeks x 8 responses
   x 3 replicates); (c) the dissolved-oxygen and headspace conditions of the stored tubes;
   (d) confirmation of the HMF LOD.
10. **Registry gaps against `data/keys/compounds.yml`.** Present and keyable: `hmf`. **Absent:**
    glucosone, 3-deoxyglucosone, threosone, methylglyoxal, glyoxal, the 1,2-enediol, the
    fructofuranosyl cation, sucrose, glucose and fructose. The registry carries `2_3_butanedione`
    (whose o-phenylenediamine derivative, 2,3-dimethylquinoxaline, is one of this paper's three
    calibration standards) but **not the two alpha-dicarbonyls this paper actually measures most**,
    glyoxal and methylglyoxal — which is a real gap, since those two are the reactive carbonyls the
    acrylamide and AGE lanes both depend on. `browning_intermediates` exists as a class id and is
    the nearest thing the registry has to a home for them.
11. **Registry gaps against `data/species/off_flavour_targets.yml`.** That file holds six compounds:
    hexanal, nonanal, 1-octen-3-ol, 2-pentylfuran, 1-hexanol and furfural. **None of the six is
    measured in this paper**, and the file has no entry for HMF or for any alpha-dicarbonyl. The
    file is about off-flavour in plant-based meat alternatives; this paper is about browning
    chemistry in acidic beverages. There is no overlap and no action.
