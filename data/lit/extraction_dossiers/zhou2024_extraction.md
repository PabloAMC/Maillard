# Zhou 2024 — EXTRACTION (fed alanine + glyoxal or + methylglyoxal, 1:1 at 20 mmol/L, pH 8, 100/110/120 C, 0-120 min; plus Ala-xylose Amadori compound ± 15N-alanine at 120 C)
DOI: 10.1021/acs.jafc.4c03706 (printed on every page; Crossref 2026-11-01 resolves to this title, J. Agric. Food Chem. 72, 18630-18637)
### The only paper in the corpus with a printed rate, three temperatures and an Ea for "amino acid + alpha-dicarbonyl -> pyrazine".

**Source on disk:** `data/articles/Zhou2024.pdf` (8 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Zhou2024.txt`); Tables 1 and 2 came through clean and are re-typed
below. Page 5 (journal p. 18634) was rasterised at 130 dpi to confirm Table 2 and to read the axis
labels of Figure 2 (y in µmol/L, x in min). No value was read off any figure. The Supporting
Information (Table S1, Figures S1-S2, including the Arrhenius plot) is NOT on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Disclosing the Nitrogen Sources via Isotope Labeling Technique and the Formation Mechanism of Pyrazine and Alkylpyrazines during the Heat Treatment of N-(1-Deoxy-D-xylulos-1-yl)-alanine and Exogenous Alanine" |
| Authors | Tong Zhou, Meigui Huang, Heping Cui, Shahzad Hussain, Khizar Hayat, Xiaoming Zhang*, Chi-Tang Ho* (Jiangnan University, Wuxi; Rutgers) |
| Venue | J. Agric. Food Chem. 2024, 72, 18630-18637. Received April 28 2024, accepted August 5 2024, published August 8 2024 |
| DOI | 10.1021/acs.jafc.4c03706 |
| Naming | Ala-ARP = N-(1-deoxy-D-xylulos-1-yl)-alanine (the alanine-xylose Amadori compound); GO = glyoxal; MGO = methylglyoxal; DXs = deoxyxylosones (1-DX, 3-DX); "pyrazine" in this paper means the unsubstituted parent C4H4N2, not the family |
| Companions | ref 19 = Zhou 2023 (JAFC 71, 2472; already in the repo as `zhou2023_extraction.md`, the source of rule R28); ref 18 = Zhou 2022 (JAFC 70, 15202; ARP preparation, GC-MS conditions); ref 20 = Zhou 2024a (JAFC 72(11), 5878; the alpha-dicarbonyl and 14N/15N-Ala UPLC methods) — none of the method details delegated to refs 18/20 are re-printed here |

## 1. Why it matters

The engine has no pyrazine lane. Rule R28 (two alpha-aminoketones -> 2,5-disubstituted pyrazine)
exists in `data/lit/reaction_rules.yml` as a net rule sourced from Zhou 2023, but Zhou 2023 has no
rate. This paper feeds alanine and one alpha-dicarbonyl at known concentrations and prints, in a
table, a formation rate for the pyrazine at 100, 110 and 120 C plus an Arrhenius Ea:

- alanine + glyoxal -> pyrazine (via R07 Strecker to aminoacetaldehyde, then R28), and
- alanine + methylglyoxal -> 2,5-dimethylpyrazine (via aminoacetone, then R28; the R28 positive
  control `CC(=O)CN` x2 -> `Cc1cnc(C)cn1` is exactly this product).

That is the lumped R07 + R28 sequence with a measured rate, units, conditions and a barrier, i.e.
the minimum the repository's rule requires before a pyrazine lane can be pre-registered. The
isotope part (Ala-ARP + 15N-Ala) gives nitrogen-source fractions that are response-factor-immune
ratios; the ARP time courses (Figures 1, 3, 4) are figure-only.

## 2. Methods as they matter to a model

- **Fed-dicarbonyl kinetic runs (the rate source).** "The reaction solutions for the Ala-GO and
  Ala-MGO models were prepared with the reactant concentrations set at 20 mmol/L, and their pH
  values were adjusted to 8.0 ± 0.02 with a 1 mol/L NaOH solution." Molar ratio 1:1 ("With a
  molar ratio of 1:1, Ala was reacted with GO, and MGO"). So **[Ala] = [GO] = 20 mmol/L** and
  **[Ala] = [MGO] = 20 mmol/L**, in ultrapure water, **no buffer** (NaOH only). Vessels:
  "high-temperature/pressure-resistant reaction vessels, undergoing continuous stirring and heating
  at different temperatures (100, 110, and 120 C) for various durations (0, 30, 60, 90, and 120
  min)". Volume, headspace and heat-up time are not stated. Quench in ice-water.
- **Reagents.** Glyoxal "40% in H2O" and methylglyoxal "40% in H2O" (Sinopharm). The 20 mmol/L is
  nominal on the commercial solution; no assay of the MGO solution is reported.
- **ARP runs.** Ala-ARP 20 mmol/L, with or without 20 mmol/L 14N-Ala or 15N-Ala; pH 8.0 ± 0.02
  with 5 mmol/L NaOH; oil bath 120 C; 30, 60, 90, 120 min; ice-water quench. The Gly variant
  (Ala-ARP + 15N-Gly, 120 C, 60 min) is in the SI only.
- **ARP preparation.** Xyl 0.4 mol/L + Ala 0.2 mol/L in 100 mL water, pH 7.5, reflux 80 C 60 min,
  vacuum evaporation 80 C 15 min, Dowex 50WX8 (H+) column, water then 0.1 mol/L ammonia elution,
  freeze-dried; identity by UPLC-Q-TOF/MS and NMR (as in ref 19). Purity not quantified here.
- **Pyrazine quantification.** HS-SPME-GC/MS: 3 g of reaction product + 1.2 g NaCl in a 20 mL
  vial + 5 µL internal standard (1,2-dichlorobenzene, 0.0018 µg/µL in methanol, i.e. **0.009 µg
  per vial**); DVB/CAR/PDMS fibre, 60 C water bath, 30 min adsorption; desorption 250 C, 10 min;
  DB-Wax 30 m x 0.25 mm x 0.25 µm. Identification: MS library + retention index against C7-C30
  alkanes (RI 1205 / 1265 / 1321 for pyrazine / methylpyrazine / 2,5-dimethylpyrazine vs NIST KI
  1210 / 1263 / 1315) + authentic standards. **Calibration: absolute, external standards at several
  concentrations with a fixed internal-standard amount**; "x (the concentration of flavor compounds
  in µg/L) against y (the ratio of the peak area of the standard compound to the peak area of the
  internal standard)":
  - pyrazine: y = 0.0088x + 0.0064, R2 = 0.9968
  - methylpyrazine: y = 0.037x + 0.0425, R2 = 0.992
  - 2,5-dimethylpyrazine: y = 0.0377x + 0.0174, R2 = 0.9992
  The matrix of the standards (water? buffer? with 20 mmol/L reactants?) is not stated. No LOD/LOQ
  printed. Concentrations are then plotted in µmol/L (Figure 2 axis), so µg/L were divided by the
  molar masses (pyrazine 80.09, methylpyrazine 94.12, 2,5-dimethylpyrazine 108.14 g/mol).
- **Other analytes (ARP runs only).** Ala-ARP by HPLC-ELSD with a purified-ARP calibration
  (XBridge BEH amide, RT 20.1 min). GO, MGO, 3-DX, 1-DX by o-phenylenediamine derivatisation (1 %
  OPD + 11 mmol/L DTPA in HEPES; 0.5 mL reagent + 0.5 mL sample, 12 h dark, room temperature) and
  HPLC-DAD on Sunfire C18 (method of ref 20). 14N-Ala and 15N-Ala by UPLC-TQD, ESI+, MRM, external
  standards. pH by meter.
- **Replicates / statistics.** "All experimental results were obtained through three repetitions,
  and the results were presented as mean ± standard deviation." ANOVA, p < 0.05.
- **Kinetic treatment (verbatim).** "Figure 2a-c indicates that the formation of pyrazine followed
  the characteristics of the zero-order kinetic model (expressed as c = kt + b). The reaction rate
  constant (k) was 0.0279, 0.0791, and 0.1507 µmol/L·min-1 at 100, 110, and 120 C, respectively
  (Table 2), and each R2 value was greater than 0.96." Also: "The time-concentration curve of the
  Ala-GO model at three temperatures was further fitted simultaneously" — but Table 2 gives one k
  and one R2 per temperature, so the fit is per temperature in effect. Arrhenius: "By plotting ln k
  versus 1/T ... For pyrazine: ln k = -12100/T + 28.87 (R2 = 0.9972), and for 2,5-dimethylpyrazine:
  ln k = -13430/T + 30.33 (R2 = 0.9997)." Ea from the slope: 100.59 and 111.66 kJ/mol. The
  intercept b of c = kt + b is never printed.

## 3. Tables re-typed

### Table 1. "Ratio of Labeled 15N to 14N in Pyrazines Formed during the Thermal Treatment of Ala-ARP and 15N-Labeled Ala"

Columns: proportions (%) of the isotopomers carrying 0, 1 or 2 15N atoms (footnote a: "Number of 15N
atoms in the isotopomers"); last column the atom ratio 15N:14N. Conditions: Ala-ARP 20 mmol/L +
15N-Ala 20 mmol/L, pH 8.0, 120 C.

| compound | time (min) | 0 x 15N | 1 x 15N | 2 x 15N | 15N:14N |
|---|---:|---:|---:|---:|---:|
| pyrazine | 30 | 20.26 | 50.54 | 29.20 | 1.20 |
| | 60 | 20.97 | 50.06 | 28.97 | 1.17 |
| | 90 | 21.26 | 50.03 | 28.71 | 1.16 |
| | 120 | 21.18 | 48.97 | 29.85 | 1.19 |
| methylpyrazine | 30 | 39.73 | 46.29 | 13.98 | 0.59 |
| | 60 | 37.75 | 47.53 | 14.72 | 0.63 |
| | 90 | 35.62 | 47.95 | 16.43 | 0.68 |
| | 120 | 33.94 | 48.00 | 18.06 | 0.73 |
| 2,5-dimethylpyrazine | 30 | 35.97 | 41.56 | 22.47 | 0.76 |
| | 60 | 33.34 | 45.51 | 21.15 | 0.78 |
| | 90 | 32.11 | 44.48 | 23.41 | 0.84 |
| | 120 | 31.60 | 44.37 | 24.03 | 0.86 |

Arithmetic check (30 min rows): 15N atoms per 100 molecules of pyrazine = 50.54 + 2 x 29.20 =
108.94; 14N = 2 x 20.26 + 50.54 = 91.06; ratio 1.196 -> printed 1.20. 15N share of N =
108.94 / 200 = 54.5 % -> the abstract's "55% of the formed pyrazine originated from exogenous Ala".
Methylpyrazine 14N share = (2 x 39.73 + 46.29)/200 = 62.9 % -> "63%"; 2,5-dimethylpyrazine 14N
share = (2 x 35.97 + 41.56)/200 = 56.8 % -> "57%". The table is internally consistent. The
Gly variant: "The ratio of 15N atoms to 14N atoms in pyrazine even reached 1.53 (Table S1)" — SI
not on disk.

### Table 2. "Kinetic Model Parameters (k and Ea) for the Formation of Pyrazine and 2,5-Dimethylpyrazine"

Conditions (Figure 2 caption): "The concentrations of Ala, GO, and MGO were all 20 mmol/L, initial
pH was 8.0, and reaction time was 0-120 min", 100/110/120 C. k is the slope of c = kt + b.

| compound | quantity | 100 C | 110 C | 120 C | Ea (kJ/mol) |
|---|---|---:|---:|---:|---:|
| pyrazine (Ala + GO) | k (µmol/L·min-1) | 0.0279 | 0.0791 | 0.1507 | 100.59 |
| | R2 | 0.9789 | 0.9918 | 0.9634 | |
| 2,5-dimethylpyrazine (Ala + MGO) | k (µmol/L·min-1) | 0.0035 | 0.0100 | 0.0230 | 111.66 |
| | R2 | 0.9891 | 0.9946 | 0.9853 | |

**Unit reconciliation.** 1 µmol/L·min-1 = 1e-3 mmol/L/min = 1e-6 mol/L/min = 1.667e-8 mol/L/s.

| compound | T (C) | k printed (µmol/L/min) | mmol/L/min | mol L-1 s-1 |
|---|---:|---:|---:|---:|
| pyrazine | 100 | 0.0279 | 2.79e-5 | 4.65e-10 |
| pyrazine | 110 | 0.0791 | 7.91e-5 | 1.32e-9 |
| pyrazine | 120 | 0.1507 | 1.507e-4 | 2.51e-9 |
| 2,5-dimethylpyrazine | 100 | 0.0035 | 3.5e-6 | 5.83e-11 |
| 2,5-dimethylpyrazine | 110 | 0.0100 | 1.00e-5 | 1.67e-10 |
| 2,5-dimethylpyrazine | 120 | 0.0230 | 2.30e-5 | 3.83e-10 |

**Arrhenius arithmetic.** Slope x R: 12100 K x 8.314 J/mol/K = 100.6 kJ/mol; 13430 x 8.314 =
111.7 kJ/mol — the printed Ea are the slopes of the printed lines. Prefactors from the intercepts:
exp(28.87) = 3.45e12 µmol/L/min; exp(30.33) = 1.49e13 µmol/L/min (lumped zero-order prefactors,
not molecular quantities). Ratio of the two k at 120 C: 0.1507 / 0.0230 = **6.55** (the text says
"5.55 times higher"; see Flags). At 110 C: 7.91; at 100 C: 7.97.

**Re-fit check (mine).** An unweighted least-squares line through the three printed (1/T, ln k)
pairs gives Ea = 103.1 kJ/mol, intercept 29.69, R2 = 0.986 for pyrazine and Ea = 114.9 kJ/mol,
intercept 31.41, R2 = 0.997 for 2,5-dimethylpyrazine. The printed pyrazine line (R2 = 0.9972)
predicts k(110 C) = 0.0665, not the printed 0.0791; so the printed lines were not fitted to the
three Table 2 k values alone (perhaps to all replicate points, or weighted). The difference is 3
kJ/mol, within any reasonable uncertainty, but the repo should carry the printed values and note
this.

### Numbers in the text (ARP runs; the underlying curves are FIGURE-ONLY)

- "the amounts of pyrazine, methylpyrazine, and 2,5-dimethylpyrazine formed in the Ala-ARP/Ala
  model were 1.40, 1.48, and 1.21 times that in the Ala-ARP model at the 120th min" (120 C, pH 8,
  20 + 20 mmol/L). Within-study ratios — usable as such.
- 2,6-dimethylpyrazine "was not detected" in either ARP model.
- Figure 3 (Ala-ARP, pH, 3-DX, 1-DX, GO, MGO vs time, both models) and Figure 4 (14N-Ala and
  15N-Ala vs time): FIGURE-ONLY; no numbers printed. Qualitative statements: ARP and DXs lower with
  added Ala; pH falls more slowly with added Ala; GO and MGO higher than in the ARP-alone model
  during 0-30 min and lower after 30 min; regenerated Ala released faster with added Ala.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): 2,5-dimethylpyrazine -> `2_5_dimethylpyrazine`;
methylpyrazine -> `methylpyrazine`; **pyrazine (the parent) -> not in registry as a molecule** (the
class row `pyrazines` carries the alias "pyrazine" with the explicit note that the bare word means
the family; a new molecule row is needed before this rate can be keyed); glyoxal, methylglyoxal,
alanine, glycine, Ala-ARP, 1-DX, 3-DX -> not in `compounds.yml` (`reaction_rules.yml` uses the short
species names GO, MGO, Ala, Gly).

| quantity | value | unit | conditions | reaction order (authors) | source location | evidence class |
|---|---|---|---|---|---|---|
| pyrazine formation rate, Ala + GO | 0.0279 / 0.0791 / 0.1507 | µmol L-1 min-1 | [Ala] = [GO] = 20 mmol/L, water, initial pH 8.0 (NaOH, unbuffered), 100 / 110 / 120 C, 0-120 min, stirred sealed vessel | zero order in product vs time (c = kt + b) | Table 2, p. 18634 | measured_rate |
| Ea, pyrazine formation, Ala + GO | 100.59 | kJ/mol | same, 100-120 C, 3 points | Arrhenius on the zero-order k | Table 2; text "ln k = -12100/T + 28.87 (R2 = 0.9972)" | measured_barrier |
| 2,5-dimethylpyrazine formation rate, Ala + MGO | 0.0035 / 0.0100 / 0.0230 | µmol L-1 min-1 | [Ala] = [MGO] = 20 mmol/L, same conditions | zero order in product vs time | Table 2 | measured_rate |
| Ea, 2,5-dimethylpyrazine, Ala + MGO | 111.66 | kJ/mol | same | Arrhenius on the zero-order k | Table 2; text "ln k = -13430/T + 30.33 (R2 = 0.9997)" | measured_barrier |
| k(Ala+GO->pyrazine) / k(Ala+MGO->2,5-DMP) at 120 / 110 / 100 C | 6.55 / 7.91 / 7.97 (text says 5.55 at 120 C) | — | same | — | derived from Table 2 | within-study ratio |
| second-order re-expression (mine, NOT the authors') k2 = k / ([Ala][GO]) | 3.77e-4 (120 C), 1.98e-4 (110 C), 6.98e-5 (100 C) | L mol-1 min-1 | assumes rate = k2[Ala][GO] with both at 0.020 mol/L and unconsumed | assumed second order | derived | derived_assumption |
| second-order re-expression, Ala + MGO -> 2,5-DMP | 5.75e-5 (120 C), 2.50e-5 (110 C), 8.75e-6 (100 C) | L mol-1 min-1 | same assumption | assumed second order | derived | derived_assumption |
| 15N share of pyrazine N (Ala-ARP + 15N-Ala) | 54.5 / 53.9 / 53.7 / 54.3 % at 30/60/90/120 min | atom % | 20 + 20 mmol/L, pH 8, 120 C | — | Table 1 | measured ratio (response-factor-immune) |
| 14N (regenerated-Ala) share of methylpyrazine N | 62.9 / 61.5 / 59.6 / 57.9 % | atom % | same | — | Table 1 | measured ratio |
| 14N share of 2,5-dimethylpyrazine N | 56.8 / 56.1 / 54.4 / 53.8 % | atom % | same | — | Table 1 | measured ratio |
| Ala-ARP/Ala vs Ala-ARP, 120 min: pyrazine, methylpyrazine, 2,5-DMP | 1.40, 1.48, 1.21 | ratio | 120 C, pH 8 | — | text p. 18632 | within-study ratio (levels themselves FIGURE-ONLY) |
| pyrazine, methylpyrazine, 2,5-DMP concentrations in ARP models vs time | — | µmol/L (axis) | 120 C | — | Figure 1 | figure_only |
| Ala-ARP, pH, 3-DX, 1-DX, GO, MGO vs time; 14N/15N-Ala vs time | — | — | 120 C | — | Figures 3, 4 | figure_only |
| pyrazine, 2,5-DMP concentration vs time in the fed runs (the raw points behind Table 2) | — | µmol/L | 100/110/120 C | — | Figure 2 | figure_only |

Cross-reference inside the repo: `data/lit/arrhenius_params.yml` carries `pyrazine_condensation`
with Ea 138.072 kJ/mol tagged `estimated`. The two measured lumped barriers here (100.6 and 111.7
kJ/mol, R07 + R28 together, pH 8, 100-120 C) are the first measured numbers for that entry's
substrate class in the corpus; they are not the same quantity as a single-step R28 barrier.

## 5. Flags

1. **"Zero order" is an initial-rate statement, not a determined order.** Over 120 min at 120 C the
   pyrazine line rises by k x t = 0.1507 x 120 = 18.1 µmol/L; two Ala and two GO are consumed per
   pyrazine, so 36 µmol/L of each reactant = **0.18 % of 20 mmol/L**. Reactant concentrations are
   constant to within analytical error, so the data cannot distinguish zero, first or second order
   in the reactants. The k values are rates at [Ala] = [dicarbonyl] = 20 mmol/L, and should be
   stored with those concentrations attached. The second-order re-expression in §4 is my
   assumption for transport into a mass-action lane, not a measurement of order. The authors'
   inference that "the formation rates ... were not influenced by the concentration of their
   respective reactants" is not supported by a single-concentration experiment.
2. **6.55 vs "5.55".** 0.1507 / 0.0230 = 6.55; the text prints 5.55. Table 2 is the number to keep.
3. **Printed Arrhenius lines do not reproduce Table 2's 110 C pyrazine k** (0.0665 predicted vs
   0.0791 printed; my 3-point refit gives Ea 103.1 kJ/mol, R2 0.986, not 0.9972). The SI Arrhenius
   plot (Figure S2) is not on disk to resolve which points were fitted.
4. **Unbuffered, initial pH only.** pH 8.0 set with NaOH; no buffer; the pH drift in the GO/MGO runs
   is not reported (only the ARP runs' pH is followed, Figure 3b, figure-only). Strecker degradation
   releases CO2 and the dicarbonyl solutions are acidic, so the working pH is unknown after t = 0.
5. **Commercial 40 % dicarbonyl solutions, nominal concentration.** MGO solutions contain hydrates,
   oligomers and formaldehyde/acetol impurities; GO solutions are largely oligomeric hydrates. The
   effective monomeric dicarbonyl concentration was not assayed.
6. **Vessel volume and headspace unstated**, and the SPME sample is 3 g of product with 1.2 g NaCl.
   Calibration standards' matrix not stated; single internal standard (1,2-dichlorobenzene) for
   three analytes of different polarity; the calibration is external-with-IS, so the numbers are
   absolute concentrations subject to matrix effects. Intercept b of each fit and any t = 0 blank
   are not printed.
7. **Three temperatures, five time points, triplicates** — R2 of the zero-order fits 0.963-0.995. No
   confidence intervals on k or Ea.
8. **Not measured in the fed runs:** methylpyrazine in the Ala+GO or Ala+MGO runs (only reported for
   the ARP models); any GO/MGO consumption; Strecker aldehyde (acetaldehyde); any mixed GO+MGO run.
   So there is no rate for the methylpyrazine channel (needs both GO and MGO) from this paper.
9. **SI not on disk**: Table S1 (Gly isotopomers), Figure S1 (Gly-model pyrazine levels), Figure S2
   (Arrhenius plots).
10. Registry gap: the parent pyrazine has no molecule row in `compounds.yml`; the class alias
    "pyrazine" must not be used for it.

## Supplementary Information (read 2026-11-01)

**Source on disk:** `data/articles/ zhou2024_supplementary.pdf` (**the filename begins with a space** [!]; 7 pp.,
Aspose-generated, ACS download stamp of 01 Nov 2026, the SI of DOI 10.1021/acs.jafc.4c03706, file jf4c03706_si_001).
Read from `pdftotext -layout`; pp. 4, 6, 7 rendered at 150 dpi and checked cell by cell (Table S1) and against the axes
(Figs. S1, S2). Contents: title and authors (pp. 1-2), a one-paragraph contents note (p. 3), **Table S1** (p. 4), figure
captions (p. 5), **Fig. S1** (p. 6), **Fig. S2** (p. 7). Legend as the other dossiers ([V] verbatim, [M] printed,
[F] figure, [Z] mine, [!] flag). Computations in the `maillard` env
(scratchpad `pyr_scripts/zhou_labels.py`, `leak_shape.py`).

### S.1 The key question (EX-052): were the 15N isotopologue data corrected for [M-H]+?

**Not stated anywhere. The SI contains no methods text at all.** Its contents note lists only Table S1, Fig. S1 and
Fig. S2 [M]. The main paper's only description of the isotopomer analysis is that isotopomers "devoid of any 15N atom
were deemed to be pyrazines generated from regenerated Ala" [V, p. 18633]. It gives no ions, no correction for the
[M-H]+ fragment, and none for 13C / 15N natural abundance. GC and MS conditions are delegated to ref 18 (Zhou 2022;
`zhou2022.pdf` on disk: EI 70 eV, full scan m/z 35-450), which has no isotope work. So **the correction is
undocumented in all three Ho-lab documents on disk**, and EX-052's S3 primary band (uncorrected leakage
`op.mh_leak_mpz` in [0, 0.30]) stays the applicable one. The SI does not decide EX-052. P(the authors corrected for
[M-H]+) ~0.4 [Z]: the lab's own Hwang 1993 method did solve an M-1 / M / M+1 / M+2 system (prereg note, hwang1995
dossier); here nothing is said.

**What the printed numbers say about leakage anyway [Z]** (Table 1, Ala pot; isotopomer σ ~1 point per fraction
assumed, since no SD is printed):
- **Uniform leakage cannot close the label deficit.** Inverting the same uncorrected [M-H]+ share r for PZ, MPZ and DMP,
  the 30-min deficit is 11.7 / 12.5 / 13.5 / 14.4 points at r = 0 / 0.1 / 0.2 / 0.3, and 12.4-15.1 with a 13C / 15N
  natural-abundance correction added. At 120 min: 8.2 -> 9.9. (The red team's version, leakage on MPZ and DMP only,
  gave 9-10 points: same conclusion.)
- **Closing it needs an MPZ-specific leak:** r_MPZ = 0.39 / 0.42 / 0.44 / 0.48 when PZ and DMP leak r = 0 / 0.05 / 0.1
  / 0.2. That is an [M-H]+ / M+ ratio of ~0.4 for methylpyrazine above 2,5-DMP's, and EI spectra of methylpyrazines do
  not obviously give that (2,5-DMP has two methyls to lose H from). Check against reference EI spectra (NIST m/z 93 / 94
  vs 107 / 108) before relying on this.
- **Isotopomer shape is weakly informative:** a binomial with leakage r shows a homozygote excess F of only ~0.01-0.06
  for r = 0.05-0.3, comparable to noise. Observed F: PZ -0.02 / 0.01, MPZ 0.008 / 0.015, DMP 0.15 / 0.11 (30 / 120 min).

### S.2 The label-average test, recomputed (every pairing rival)

Pair-specific constants (P1x) and pair-specific dihydropyrazine fates (P1y) do not move labels, so under P1, P1x and
P1y, with the A and M pools each homogeneous in time: label(MPZ) = (label(PZ) + label(DMP)) / 2, label = 15N atom
share.

| pot | t (min) | PZ | DMP | MPZ observed | predicted | deficit (points) |
|---|---:|---:|---:|---:|---:|---:|
| Ala-ARP + 15N-Ala (Table 1) | 30 | 0.545 | 0.432 | 0.371 | 0.489 | **11.7** |
| | 60 | 0.540 | 0.439 | 0.385 | 0.490 | 10.5 |
| | 90 | 0.537 | 0.457 | 0.404 | 0.497 | 9.3 |
| | 120 | 0.543 | 0.462 | 0.421 | 0.503 | 8.2 |
| Ala-ARP + 15N-Gly (Table S1), DMP from the printed ratio 0.59 | 60 | 0.606 | 0.371 | 0.461 | 0.488 | 2.7 |
| same, DMP with n2 = 100 - n0 - n1 | 60 | 0.606 | 0.421 | 0.461 | 0.513 | 5.2 |
| same, DMP row as printed, renormalised | 60 | 0.606 | 0.528 | 0.461 | 0.567 | 10.6 |

**Interval labels** (Table 1 labels x Fig. 1 amounts, µg/L read off the main paper's Fig. 1, converted to µmol/L; read
error ±0.3-1 µg/L, so ±0.1-0.2 on the later intervals) [F][Z]:

| interval (min) | PZ | DMP | MPZ | predicted MPZ | deficit |
|---|---:|---:|---:|---:|---:|
| 0-30 | 0.545 | 0.432 | 0.371 | 0.489 | 11.7 |
| 30-60 | 0.52 | 0.45 | 0.42 | 0.49 | ~7 |
| 60-120 | 0.62 | 0.50 | 0.60 | 0.56 | ~-4 (within error) |

**The deficit is front-loaded:** it sits in MPZ made in the first 30-60 min; MPZ made after 60 min carries the label the
identity predicts, within the figure-read error.

**Does it separate pairing (P1x / P1y) from an extra MPZ route (P2b)?** Partly, and only through shape and timing; the
mean-label test alone does not.

1. **Mean label:** the identity fails at 30-120 min in the Ala pot, robust to any uniform leakage. Two escapes remain:
   (a) P2b, an extra MPZ route whose N label is lower than the pairing average; (b) pairing with pool heterogeneity in
   time, where M halves made while A is abundant (early) carry a low label (~0.20 needed at 30 min vs 0.43 for DMP's
   M halves). DMP's own homozygote excess (F 0.15 at 30 min) shows the M pool is heterogeneous, so (b) is not excluded.
2. **Isotopomer shape** (σ 1 point per fraction assumed):
   - MPZ is binomial at its own label (39.7 / 46.3 / 14.0 vs 39.5 / 46.7 / 13.8 for p = 0.371).
   - Escape (b) keeps the A half at PZ's label (0.545, which is itself binomial and stable over time), so MPZ = A
     (0.545) x M (free label). That gives a heterozygote excess: best fit χ² 60 / 32 / 21 / 21 at 30 / 60 / 90 / 120 min
     for 1 dof.
   - Uncorrected leakage only moves F by +0.01-0.06, so it cannot rescue this.
   - P2b with an extra route of label q (both N from one mixed pool) fits exactly: w 0.6-1.0, q 0.35-0.39. But that
     is 2 parameters for 2 dof, so it is admitted by flexibility, not tested.
   - P2b with an unlabelled extra route (q = 0, N from ARP or regenerated Ala only) fits poorly: χ² 30-51.
3. **Glycine pot (Table S1)** gives the S3 statistic D = P2 / f_A + P0 / (1 - f_A) - 1 = **-0.085**. Pairing gives
   D = 0, uncorrected leakage pushes D positive, and the Ala pot gives +0.13. So MPZ has a heterozygote excess (F -0.17)
   that no pairing rival with PZ's A label produces and leakage cannot cause. With the DMP row unusable (S.4) and one
   time point, this row carries little weight.

**Consequence.**
- **P1x / P1y as the sole explanation: disfavoured** (red team's "leaning rejected" stands). The shape argument adds
  a second, partly independent strike. It assumes ~1 point isotopomer noise, which is unprinted.
- **P2b: still provisional.** It fits, but with free parameters, and the fitted extra route is a single mixed-label N
  pool (q ~0.37), not an ARP-N-only route.
- **Time structure:** the extra MPZ, or the label heterogeneity, is confined to the first 30-60 min, where ARP decay
  and GO / MGO release peak (main paper Fig. 3, figure-only).
- **The decisive datum is still the authors' raw ion intensities** (m/z 79-82, 93-96, 107-110, per time point) or
  their correction method. Ask Zhang / Ho before EX-052 is re-scored.

### S.3 Table S1 (Ala-ARP 20 mM + 15N-Gly 20 mM, pH 8.0 initial, 120 °C, 60 min) [M]

| compound | 0 x 15N (%) | 1 x 15N (%) | 2 x 15N (%) | 15N:14N printed | row sum | 15N:14N recomputed [Z] |
|---|---:|---:|---:|---:|---:|---:|
| pyrazine | 14.65 | 49.59 | 35.75 | 1.53 | 99.99 | 1.535 |
| methylpyrazine | 24.92 | 57.94 | 17.14 | 0.86 | 100.00 | 0.856 |
| 2,5-dimethylpyrazine | 37.52 | 40.84 | **44.48** | 0.59 | **122.84** | 1.12 as printed; 0.73 with n2 = 21.64 |

### S.4 Flags

1. **Table S1's 2,5-DMP row is internally inconsistent** [!]: it sums to 122.84 %, and no single-cell correction
   reproduces the printed ratio 0.59. n2 = 100 - n0 - n1 = 21.64 gives 0.73. Keeping n0 and n1, a ratio of 0.59 needs
   n2 = 13.8 (sum 92). The render confirms the printed digits. The Gly-pot DMP label is therefore unknown (0.37-0.53).
2. **Fig. S1's y axis reads "Concentration (mmol/L)"** [!]. That is impossible: ~70 mmol/L pyrazine from 20 mmol/L ARP.
   The ARP bars (~25.5 / ~5.5 / ~2.0 [F]) match the main paper's Fig. 1 ARP values at 60 min in µg/L (25.5 / 5.7 / 2.1),
   so the unit is µg/L.
   - ARP + Gly at 60 min: PZ ~69 (±4), MPZ ~13 (±1), DMP ~7.0 (±0.2) µg/L [F].
   - Molar R: ARP 0.58, ARP + Gly 0.34 [Z].
   - Fig. 1 gives molar R for the Ala pot (ARP + 15N-Ala) of 0.61 / 0.64 / 0.62 / 0.53 at 30-120 min, and for ARP
     alone 1.10 / 0.59 / 0.45 / 0.41 [Z].
   - Every one of these is below 4 and several are at or below the rivals' floor 4h_any = 0.4, the same low-R pattern
     as CF-087.
   - A caveat specific to labelled pots: the SPME calibration used unlabelled standards. If the quantifier ion was the
     molecular ion of the unlabelled compound, 15N isotopologues would be under-counted, in proportion to label. The
     quantifier ions are not stated.
3. **Fig. S2 (Arrhenius) shows exactly three points per compound**, at ln k = ln of Table 2's k [F]: pyrazine -1.89 /
   -2.54 / -3.58; 2,5-DMP -3.77 / -4.61 / -5.65. They are joined point to point, not by a fitted line. This answers
   flag 3 of this dossier: no extra replicate points lie behind the printed equations. So the printed line
   ln k = -12100/T + 28.87 (R² 0.9972) is still not reproduced by the only data shown; the 3-point refit is
   Ea 103.1 kJ/mol. Carry Table 2's k values; treat the printed Ea as ±3 kJ/mol.
4. Flag 9 above ("SI not on disk") is closed by this section.

### S.5 Disagreements found (new)

1. **Table S1 internal** (S.4.1): 2,5-DMP isotopomers sum to 122.84 % and do not reproduce the printed 15N:14N.
2. **Fig. S1 unit** (S.4.2): mmol/L printed, µg/L by cross-check with the main paper's Fig. 1.
3. **Printed Arrhenius line vs the plotted points** (S.4.3): Ea 100.59 printed vs 103.1 refit from the only three
   points; R² 0.9972 printed vs 0.986.
4. **Glycine pots, low R, two labs:** here 0.34 (Ho, SPME); bemisyoung1993 0.02-0.54 at pH 6-9 (UC Davis, calibrated
   capillary FID). Against CF-087's "low-level data" explanation.

### S.6 Ledger-ready (additions)

1. Table S1 PZ and MPZ isotopomers (Gly pot, 60 min). Key the DMP row only with the sum flag.
2. Fig. S1 bars (µg/L, unit corrected), two readers.
3. Table 1 + Fig. 1 interval labels: derived, not keyable; record in the EX-052 re-score instead.
