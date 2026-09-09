# Zhou 2025b — EXTRACTION (Ala-xylose Amadori compound ± exogenous xylose, 20 + 20 mmol/L, initial pH 8, 120 C, 30-120 min; 13C5-xylose CAMOLA; zero-order k for Ala + Xyl -> AX-ARP and Ala + GO -> pyrazine at 70/80/90 C)
### The exogenous-XYLOSE twin of Zhou 2024 (exogenous alanine): same lab, same ARP, same 20 + 20 mmol/L / pH 8 / 120 C design, same SPME calibration. Adds the first printed rate and Ea for "Ala + Xyl -> Amadori compound" and a second, lower-temperature ladder for "Ala + GO -> pyrazine" that does NOT sit on the 2024 ladder.

**Source on disk:** `data/articles/Zhou2025b.pdf` (10 pp., owner's download, 2026-09-08). Read from the
text layer (`scratchpad/articles/Zhou2025b.txt`); Tables 1 and 2 came through clean and are re-typed
below. Page 3 was re-extracted with `pdftotext -layout` to recover the calibration-line exponents (the
plain layer drops the minus signs), and pages 4-7 were rasterised at 110 dpi to read the axis labels
of Figures 1-4 and the layout of Table 1. No value was read off any figure. The paper has no
Supporting Information. NOTE: `data/articles/Zhou2025.pdf` / `Zhou2025_extraction.md` is a different
paper (coffee); this dossier is the Food Chemistry ARP paper and is keyed `zhou2025b`.

## 0. Identity

| field | value |
|---|---|
| Title | "Role of exogenous xylose during its co-heating with alanine-xylose Amadori compound: Competitively promoting 2-furfural formation and limiting pyrazines generation" |
| Authors | Tong Zhou, Qingqing Chai, Man Zhang, Heping Cui, Ahmad Mohammad, Khizar Hayat, Xiaoming Zhang*, Chi-Tang Ho* (Jiangnan University, Wuxi; Qilu Univ. of Technology; Yangzhou Univ.; Anhui Qiangwang; King Saud Univ.; Parkland College; Rutgers) |
| Venue | Food Chemistry 487 (2025) 144764. Received 8 March 2025, revised 10 May 2025, accepted 12 May 2025, online 14 May 2025 |
| DOI | 10.1016/j.foodchem.2025.144764 (PII S0308-8146(25)02015-1) |
| Naming | AX-ARP = N-(1-deoxy-D-xylulos-1-yl)-alanine (= Ala-ARP of Zhou 2024; repo key `ARP` in `data/species/structures.yml`); Xyl = D-xylose; GO = glyoxal; MGO = methylglyoxal; 1-DX / 3-DX = 1-/3-deoxyxylosone (repo `DPO` / `TDP`); "pyrazine" = the unsubstituted parent (repo `PZ`), "pyrazines" = the family; CAMOLA = carbon module labelling |
| Companions | Zhou 2024 (JAFC 72, 18630; `zhou2024_extraction.md`; the exogenous-Ala twin and the source of the B18 pyrazine rates); Zhou, Huang, Cui, Chen et al. 2024 (JAFC 72(11), 5878; exogenous Ala inhibits 2-furfural; the dicarbonyl method); Zhou 2023 (JAFC 71, 2472; `zhou2023_extraction.md`; ARP preparation, GC-MS conditions); Zhou 2022 (JAFC 70, 15202; cysteine); Zhan 2020 (Food Chem 311, 125877; glycine-ribose ARP zero-order formation, cited as precedent) |

## 1. Why it matters

Wave B18 (`results/validation/kinetic_core_b18_prereg.md` §6) fitted the pyrazine step on Zhou 2024's
six fed-dicarbonyl rates (Ala + GO -> pyrazine, Ala + MGO -> 2,5-dimethylpyrazine; 20 + 20 mmol/L, pH 8,
100/110/120 C) and named this group's later work as the corpus for a follow-up wave. This paper is
that follow-up from the group's side. It contains:

1. **A second zero-order ladder for Ala + GO -> pyrazine** at 70, 80, 90 C with an Arrhenius Ea of
   104.97 kJ/mol (2024: 100.59 kJ/mol at 100-120 C). The Ea repeats 2024 within 4 kJ/mol, but the
   ABSOLUTE rates are ~75x higher than the 2024 line extrapolated to the same temperatures, and the
   reactant concentration of this run is not printed (§5, flag 1-2). It cannot be merged into the
   B18 fit rows as they stand.
2. **The first printed rate and Ea for Ala + Xyl -> AX-ARP** (zero order, 70-90 C, Ea 83.12 kJ/mol),
   i.e. the lumped Schiff + Amadori step for a pentose + alanine — the quantity the repo's
   `arrhenius_params.yml` `amadori` / `schiff_condensation` entries and the Martins glucose-glycine
   refit argue about, now for xylose. Concentration again unstated.
3. **The mechanism claim** the model would need to represent: free xylose captures the alanine
   regenerated from the ARP (Ala + Xyl -> ARP outcompetes Ala + GO -> pyrazine; the authors' ground is
   the lower Ea and the faster rate), so Ala cycles catalytically through the ARP, DXs accumulate,
   the pH falls faster, 3-DX dehydrates to 2-furfural rather than retro-aldolising to GO/MGO, and
   the less-nucleophilic Ala at lower pH makes fewer pyrazines even though GO is HIGHER. Note that
   the competition is for the AMINE (upstream of the Strecker step), not for the dicarbonyl.
4. A 13C5-xylose isotopomer table (Table 2) that fixes, as response-factor-immune ratios, how much
   of the 2-furfural and pyrazine carbon comes from the free sugar (2-furfural: 3 -> 14 % of
   molecules over 30-120 min; pyrazine and methylpyrazine: ~6-8 % of carbon, flat in time).

Everything else (pyrazine, methylpyrazine, 2,5-DMP, 2,6-DMP, 2-furfural, 3-DX, 1-DX, GO, MGO, AX-ARP,
Ala, A420 and pH against time in the two 120 C models; the Xyl-alone dicarbonyl run) is FIGURE-ONLY.

## 2. Methods as they matter to a model

Conventions shared with Zhou 2024 are marked (= 2024); differences are marked (NEW / differs).

- **AX-ARP preparation (= 2024, with small differences).** "Ala (1.78 g) and Xyl (6 g) were mixed
  in 100 mL deionized water at a 1:1 M ratio" — the masses are 20.0 mmol Ala (MW 89.09) and 40.0
  mmol Xyl (MW 150.13), i.e. **0.2 mol/L Ala + 0.4 mol/L Xyl, 1:2**, exactly the charge the 2024
  dossier records; the printed "1:1" is wrong. pH 7.5 with 3 mol/L NaOH; 80 C water bath 60 min;
  rotary-evaporator dehydration 80 C **20 min** (2024: 15 min); ice bath; Dowex 50WX8 H+, water at
  3 mL/min then **0.18 mol/L ammonia** at 1 mL/min (2024: 0.1 mol/L); fractions checked by HPLC-ELSD
  (XBridge BEH Amide 4.6 x 150 mm, 3.5 µm; 0.1 % formic acid / acetonitrile gradient; ELSD 55 C,
  N2 1.5 L/min); pooled and lyophilised. Purity not quantified; identity by reference to Zhou 2023.
- **ARP runs (= 2024 design).** "The initial concentrations of AX-ARP and Xyl were both set at 20
  mmol/L. The pH of each system was adjusted to 8 before transferring the solution into the
  pressure-resistant flask." Oil bath **120 C**, **30, 60, 90, 120 min**, ice-water quench. Arms:
  AX-ARP alone; AX-ARP + Xyl; AX-ARP + 13C5-Xyl (98 atom % 13C, replacing 12C-Xyl "under identical
  conditions"); and Xyl alone (Fig. 2a; "identical reaction conditions", concentration therefore
  presumably 20 mmol/L, not restated). The base used to set pH 8 is not stated for these runs (2024:
  5 mmol/L NaOH); no buffer; solvent water. Volume, headspace, heat-up time not stated.
- **Kinetic runs (the rate source; differs from 2024).** "The reaction solutions for the Ala-GO and
  Ala-Xyl models were prepared with **equimolar amounts** of Ala and GO or Ala and Xyl. The initial
  pH of both models was adjusted to 8 using a 1 mol/L NaOH solution." Stirred high-temperature /
  pressure-resistant vessels at **70, 80 or 90 C**, sampled at **0, 15, 30, 45, 60, 75, 90 and 120
  min**, ice-water quench. Figure 4 caption: "The molar ratios of Xyl to Ala and GO to Ala were both
  1:1, and reaction time was 0-120 min". **The absolute concentration is never printed** (2024
  printed 20 mmol/L for the Ala-GO / Ala-MGO runs at 100-120 C). Glyoxal is the commercial "40 % in
  H2O" solution (Sigma), nominal. No MGO kinetic run in this paper.
- **pH.** Meter (Mettler-Toledo). Followed in the ARP runs only (Fig. 3d, axis 4-8); not reported
  for the kinetic runs.
- **Alpha-dicarbonyls (= 2024 method, IS now specified).** OPD derivatisation: reagent = HEPES 2.38 g
  + OPD 100 mg + DTPA 43.3 mg in 10 mL water; 0.5 mL sample + 0.5 mL reagent + **10 µL 3,4-hexanedione
  86 mmol/L as internal standard**; 25 C dark 12 h; HPLC-UV 315 nm, SunFire C18 5 µm 4.6 x 150 mm,
  water/0.1 % formic acid vs methanol gradient, 1 mL/min; RT 15.6 min 1-DX, 16.1 min 3-DX, 18.2 min
  GO, 24.1 min MGO. Calibration not described here (delegated to Zhou, Huang, Cui, Chen et al. 2024).
- **Pyrazines and 2-furfural (= 2024 SPME method; furfural and 2,6-DMP calibrations NEW).**
  HS-SPME-GC/MS: **3 mL** sample (2024: 3 g) + 1.2 g NaCl in a 20 mL vial + 5 µL 1,2-dichlorobenzene
  0.0018 mg/mL in methanol (= **0.009 µg IS per vial**, = 2024); DVB/CAR/PDMS 50/30 µm, fibre 1 cm
  above the liquid, 60 C, 30 min, stirred; desorption 250 C 10 min splitless; Agilent 7890B/5977B;
  column and oven per Zhou 2023 (DB-Wax in 2024). Identification: LRI vs C7-C30 alkanes, NIST 17,
  authentic standards (2-furfural >= 99.5 %, pyrazine >= 99 %, methylpyrazine, 2,5-DMP, 2,6-DMP >= 98 %).
  **Calibration: "combined internal and external standard method" (Fan et al. 2018)**, y = peak-area
  ratio analyte/IS, x = standard concentration (unit not stated in this paper; 2024 states µg/L, and
  Fig. 1 is plotted in µg/L):
  - 2-furfural: y = 8.9e-3 x + 2.9e-1, R2 = 0.9988 (NEW)
  - pyrazine: y = 8.8e-3 x + 6.4e-3, R2 = 0.9968 (identical to 2024's 0.0088x + 0.0064)
  - methylpyrazine: y = 3.7e-2 x + 4.25e-2, R2 = 0.992 (identical to 2024)
  - 2,5-dimethylpyrazine: y = 3.77e-2 x + 1.74e-2, R2 = 0.9992 (identical to 2024)
  - 2,6-dimethylpyrazine: y = 3.41e-2 x + 3.52e-2, R2 = 0.999 (NEW)
  The three shared lines are the SAME calibration as 2024, not a re-determination. No LOD/LOQ; the
  standards' matrix is not stated. The 2-furfural intercept (0.29) equals the signal of ~33 µg/L on
  its own slope, i.e. a non-trivial blank or offset (my arithmetic).
- **AX-ARP in the Ala-Xyl kinetic run (NEW).** UPLC-TQD, ESI+, MRM m/z 222 -> 204, cone 20 V,
  collision 8 V; BEH Amide 1.7 µm 2.1 x 100 mm; acetonitrile / 0.1 % formic acid gradient, 0.2
  mL/min; **external standard = the purified AX-ARP**. How AX-ARP, regenerated Ala and A420 were
  measured in the 120 C ARP runs (Fig. 3a-c) is NOT described in this paper (2024: HPLC-ELSD for
  ARP, UPLC-TQD MRM for 14N/15N-Ala).
- **Replicates / statistics.** "All experiments were performed in triplicate", mean ± SD; ANOVA with
  Duncan, p < 0.05. Error bars are drawn in every figure.
- **Kinetic treatment (verbatim where it matters).** "The fitted solid line (expressed as: c = kt + b)
  ... the formation of AX-ARP during the initial stage of the reaction was observed to follow a
  zero-order kinetic model ... The rate constant k at 70, 80 and 90 C was 0.0034, 0.0060, and 0.017
  mmol/(L·min)"; "ln k = -9997/T + 23.37, R2 = 0.9664 ... Ea value for AX-ARP formation was 83.12
  kJ/mol". Pyrazine: "The rate constant (k) at 70, 80 and 90 C was 0.1165, 0.2806, and 0.8867
  µmol/(L·min) ... As the reaction progressed, the formation rate of pyrazine showed a decreasing
  trend under all three temperatures, and the deviation from zero-order kinetics was observed ...
  Focusing on the initial stage ... ln k = -12,626/T + 34.59, with an R2 value of 0.9913. The Ea ...
  was ... 104.97 kJ/mol". At 90 C AX-ARP "reached a peak at the 60th min, after which it declined".
  By inspection of Fig. 4 the red fitted segments cover only the initial linear part and get shorter
  as temperature rises; the fitted time window is not stated numerically. The intercept b is never
  printed.

## 3. Tables re-typed

### Table 1. "Kinetic model parameters (k and Ea) for the formation of AX-ARP and pyrazine."

Conditions: Ala : Xyl = 1 : 1 and Ala : GO = 1 : 1 (absolute concentration unstated), initial pH 8
(1 mol/L NaOH, unbuffered), stirred sealed vessels, 0-120 min, zero-order fit c = kt + b to the
initial stage. Note the two k rows are printed in DIFFERENT units.

| system | quantity | 70 C | 80 C | 90 C | Ea (kJ/mol) |
|---|---|---:|---:|---:|---:|
| AX-ARP (Ala + Xyl) | k (mmol/L·min-1) | 0.0034 | 0.0060 | 0.0170 | 83.12 |
| | R2 | 0.9950 | 0.9817 | 0.9786 | |
| pyrazine (Ala + GO) | k (µmol/L·min-1) | 0.1165 | 0.2806 | 0.8867 | 104.97 |
| | R2 | 0.9803 | 0.9653 | 0.9649 | |

**Unit reconciliation** (1 µmol/L·min-1 = 1e-6 mol/L/min = 1.667e-8 mol L-1 s-1):

| system | T (C) | k printed | µmol/L/min | mol L-1 s-1 |
|---|---:|---:|---:|---:|
| AX-ARP | 70 | 0.0034 mmol/L/min | 3.4 | 5.67e-8 |
| AX-ARP | 80 | 0.0060 mmol/L/min | 6.0 | 1.00e-7 |
| AX-ARP | 90 | 0.0170 mmol/L/min | 17.0 | 2.83e-7 |
| pyrazine | 70 | 0.1165 µmol/L/min | 0.1165 | 1.94e-9 |
| pyrazine | 80 | 0.2806 µmol/L/min | 0.2806 | 4.68e-9 |
| pyrazine | 90 | 0.8867 µmol/L/min | 0.8867 | 1.48e-8 |

k(AX-ARP) / k(pyrazine) in the same units: 29.2 (70 C), 21.4 (80 C), 19.2 (90 C).

**Arrhenius arithmetic (mine).** 9997 K x 8.314 = 83.12 kJ/mol; 12626 x 8.314 = 104.97 kJ/mol: the
printed Ea are the printed slopes. Unlike 2024, an unweighted least-squares line through the three
printed (1/T, ln k) pairs reproduces the printed lines exactly: AX-ARP slope -9997, intercept 23.37,
R2 0.9663; pyrazine slope -12627, intercept 34.59, R2 0.9913. So the Table 1 k values ARE the fitted
points. Prefactors: exp(23.37) = 1.41e10 mmol/L/min; exp(34.59) = 1.05e15 µmol/L/min (lumped
zero-order prefactors at the unstated concentration; not molecular quantities).

**Cross-paper check against Zhou 2024 (mine).** 2024's printed line for the same reaction, ln k =
-12100/T + 28.87 (k in µmol/L/min at 20 + 20 mmol/L), predicts k = 0.00168 / 0.00455 / 0.0117
µmol/L/min at 70 / 80 / 90 C; this paper prints 0.1165 / 0.2806 / 0.8867, i.e. **69x / 62x / 76x
higher**. Conversely this paper's line predicts 2.12 / 5.14 / 11.9 µmol/L/min at 100 / 110 / 120 C
against 2024's printed 0.0279 / 0.0791 / 0.1507 (76x / 65x / 79x). The two ladders have the same
slope within 4 kJ/mol and an intercept offset of exp(34.59 - 28.87) = 305 on the intercept, ~70x at
the measured temperatures. See §5 flag 2.

### Table 2. "Isotope distribution patterns of 2-furfural and pyrazines formed in the thermal reaction of AX-ARP-[13C5-Xyl] model at 120 C and pH 8."

Footnote a: reaction time. Footnote b: "Proportions of isotopomers of the given compound were
calculated by normalization of peak areas of detected ions from [M] to [M + n] ... with correction
for the naturally occurring abundance of 13C (1.10 %) subtracted from the measured intensities. The
values were also corrected based on the isotopic purity (98 atom %) of the 13C5-labeled xylose."
M+ = molecular ion. "–" = not observed. Every row re-sums to 100.0 ± 0.01 %.

| compound (M+, m/z) | time (min) | M | M+1 | M+2 | M+3 | M+4 | M+5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2-furfural (96) | 30 | 96.87 | – | – | – | – | 3.13 |
| | 60 | 92.74 | – | – | – | – | 7.26 |
| | 90 | 88.79 | – | – | – | – | 11.21 |
| | 120 | 86.34 | – | – | – | – | 13.66 |
| pyrazine (80) | 30 | 85.10 | – | 12.94 | – | 1.96 | – |
| | 60 | 86.31 | – | 12.74 | – | 0.95 | – |
| | 90 | 87.78 | – | 11.37 | – | 0.86 | – |
| | 120 | 88.10 | – | 11.24 | – | 0.65 | – |
| methylpyrazine (94) | 30 | 84.97 | – | 15.03 | – | – | – |
| | 60 | 87.26 | – | 12.74 | – | – | – |
| | 90 | 84.98 | – | 10.46 | 4.56 | – | – |
| | 120 | 87.00 | – | 8.64 | 4.36 | – | – |

Derived (mine): share of the compound's carbon that came from free xylose = pyrazine (2 x M+2 + 4 x
M+4)/4 = **8.4 / 7.3 / 6.5 / 6.3 %** at 30/60/90/120 min; methylpyrazine (2 x M+2 + 3 x M+3)/5 =
**6.0 / 5.1 / 6.9 / 6.1 %**; 2-furfural = the M+5 column (all-or-none skeleton) = **3.1 -> 13.7 %**.
For pyrazine, a single well-mixed GO pool with labelled fraction f gives M : M+2 : M+4 = (1-f)^2 :
2f(1-f) : f^2; from the M column f = 0.078 / 0.071 / 0.063 / 0.061, predicting M+2 = 14.3 / 13.2 /
11.8 / 11.5 and M+4 = 0.60 / 0.50 / 0.40 / 0.38 — M+2 matches, M+4 is 3x high at 30 min and ~1.7x
high later. Roughly consistent with one GO pool that is 6-8 % xylose-derived and not changing in
time, which is the authors' reading ("did not show obvious fluctuations"). 2,5-DMP and 2,6-DMP are
absent from Table 2 because they were "entirely absent in the AX-ARP-Xyl model".

### Numbers printed in the text (the underlying curves are FIGURE-ONLY)

- Xyl alone, 120 C, 120 min: "Only 1.1 % of the initial Xyl was converted into α-dicarbonyl
  compounds" (3-DX + GO + MGO; GO the most abundant). At 20 mmol/L that is ~0.22 mmol/L total (my
  arithmetic; the charge is inferred from "identical conditions").
- "at 120 min, the concentration of GO in the AX-ARP-Xyl model was **0.36 mmol/L higher** than that
  in the AX-ARP model, whereas the degradation of Xyl itself contributed only **0.12 mmol/L**."
- Pyrazine, methylpyrazine, 2,5-DMP, 2,6-DMP all lower with Xyl; **2,5-DMP and 2,6-DMP "entirely
  absent"** in AX-ARP-Xyl; in AX-ARP-Xyl the pyrazine curves "stabilized with a slight decrease".
  2-furfural higher with Xyl, "rapid increase over time".
- 3-DX and 1-DX higher with Xyl; GO "noticeable rise"; MGO "largely unchanged". AX-ARP "slightly
  higher" with Xyl; regenerated Ala "no significant difference"; A420 "no significant difference";
  pH falls in both, "obviously faster" with Xyl.
- Xyl-Ala vs GO-Ala: Ea 83.12 vs 104.97 kJ/mol; k(ARP) > k(pyrazine) at every temperature.

### Figure inventory (axis labels read from the rasterised pages; no values taken)

| figure | content | axes | status |
|---|---|---|---|
| Fig. 1a-d | pyrazine, methylpyrazine, 2,5-DMP, 2,6-DMP vs time, AX-ARP and AX-ARP-Xyl, 120 C | µg/L (ordinate spans 0-15 for a-c, 0-10 for d) vs min (0-120) | figure_only |
| Fig. 1e | 2-furfural vs time, both models | µg/L (ordinate spans 0-2500) vs min | figure_only; 2-furfural is two orders of magnitude above any pyrazine on these axes |
| Fig. 2a | GO, 3-DX, MGO from Xyl alone, 120 C | mmol/L (0-0.8) vs min | figure_only |
| Fig. 2b-e | 3-DX, 1-DX, GO, MGO vs time, both models | mmol/L (0-0.8; MGO 0-0.6) vs min | figure_only |
| Fig. 3a | AX-ARP vs time, both models (starts at 20 mmol/L) | mmol/L (0-20) vs min | figure_only |
| Fig. 3b | regenerated Ala vs time | mmol/L (0-20) vs min | figure_only; method not described here |
| Fig. 3c | A420 vs time | absorbance (0-1.0) vs min | figure_only |
| Fig. 3d | pH vs time (starts at 8) | pH (4-8) vs min | figure_only |
| Fig. 4a | AX-ARP vs time, Ala-Xyl, 70/80/90 C, with zero-order fit segments | mmol/L (0-2.0) vs min (0-120) | figure_only (the fits are Table 1) |
| Fig. 4b | pyrazine vs time, Ala-GO, 70/80/90 C, with fit segments | µmol/L (0-40) vs min | figure_only |
| Fig. 4c-d | ln k vs 1/T | ln k (-6 to -3.5; -2.5 to 0) vs 1/T (0.0027-0.0030 K-1) | the lines are printed in the text |

## 4. Numbers / routes the repository can use

Registry mapping: pyrazine -> `PZ` (structures.yml, B18); methylpyrazine -> `MPZ` / `methylpyrazine`;
2,5-DMP -> `DMP` / `2_5_dimethylpyrazine`; 2,6-DMP -> `2_6_dimethylpyrazine` (registry; no structures
key); 2-furfural -> `furfural`; GO, MGO, Ala, ARP, DPO (1-DX), TDP (3-DX), PENT (xylose as open-chain
aldopentose) -> `data/species/structures.yml`.

| quantity or route | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| AX-ARP formation rate, Ala + Xyl | 0.0034 / 0.0060 / 0.0170 | mmol L-1 min-1 | Ala : Xyl 1 : 1, **absolute concentration unstated**, water, initial pH 8 (NaOH, unbuffered), 70 / 80 / 90 C, zero-order fit to the initial stage of 0-120 min | Table 1, p. 7; text §3.3 | measured_rate (concentration-anchored, order not determined) |
| Ea, AX-ARP formation, Ala + Xyl | 83.12 | kJ/mol | same, 70-90 C, three points; ln k = -9997/T + 23.37, R2 0.9664 | Table 1; text | measured_barrier |
| pyrazine formation rate, Ala + GO | 0.1165 / 0.2806 / 0.8867 | µmol L-1 min-1 | Ala : GO 1 : 1, **absolute concentration unstated**, initial pH 8, 70 / 80 / 90 C, zero-order initial stage | Table 1 | measured_rate — **~70x above the Zhou 2024 ladder for the same reaction; do not merge (flag 2)** |
| Ea, pyrazine formation, Ala + GO | 104.97 | kJ/mol | same; ln k = -12626/T + 34.59, R2 0.9913 | Table 1 | measured_barrier (second determination; 2024: 100.59 at 100-120 C) |
| k(Ala + Xyl -> ARP) / k(Ala + GO -> pyrazine) | 29.2 / 21.4 / 19.2 at 70 / 80 / 90 C | — | same runs, same lab, same day-design; both zero-order initial rates | derived from Table 1 | within_study_ratio (holds only if the two runs used the same reactant concentration, which is not stated) |
| 2-furfural molecules carrying the intact free-xylose skeleton (M+5) | 3.13 / 7.26 / 11.21 / 13.66 % at 30/60/90/120 min | % of molecules | AX-ARP 20 + 13C5-Xyl 20 mmol/L, pH 8, 120 C | Table 2 | within_study_ratio (response-factor-immune) |
| pyrazine carbon from free xylose | 8.4 / 7.3 / 6.5 / 6.3 % | atom % C | same | derived from Table 2 | within_study_ratio |
| methylpyrazine carbon from free xylose | 6.0 / 5.1 / 6.9 / 6.1 % | atom % C | same | derived from Table 2 | within_study_ratio |
| GO excess in AX-ARP-Xyl over AX-ARP at 120 min | 0.36 | mmol/L | 20 + 20 mmol/L, pH 8, 120 C | text §3.2 | level_only (difference of two figure-only curves, printed) |
| GO from Xyl alone at 120 min | 0.12 | mmol/L | Xyl (presumably 20 mmol/L), pH 8, 120 C | text §3.2 | level_only |
| Xyl alone -> total α-dicarbonyls (3-DX + GO + MGO) at 120 min | 1.1 % of initial Xyl | mol % | same | text §3.2 | level_only |
| 2,5-DMP and 2,6-DMP in AX-ARP-Xyl | not detected at any time point | — | 20 + 20 mmol/L, pH 8, 120 C | text §3.1, Fig. 1c-d | level_only (non-detect; LOD not printed) |
| ordering claims with Xyl added: pyrazines lower, 2-furfural higher, 3-DX / 1-DX / GO higher, MGO same, ARP slightly higher, Ala same, A420 same, pH lower | direction only | — | same | text §3.1-3.4, Figs 1-3 | figure_only (directions printed, magnitudes not) |
| mechanism: free Xyl + regenerated Ala -> ARP outcompetes Ala + GO -> pyrazine; Ala cycles; DXs accumulate; faster pH drop routes 3-DX to 2-furfural and lowers Ala nucleophilicity | — | — | authors' scheme, grounded on Table 1 (Ea, k) and Table 2 | §3.3-3.4, Conclusion | mechanism_drawn |
| all time courses in Figs 1-4 | — | µg/L, mmol/L, µmol/L, pH, A420 | 120 C ARP runs; 70-90 C kinetic runs | Figures 1-4 | figure_only |

### 4b. New versus repeated relative to Zhou 2024

| item | status |
|---|---|
| k and Ea for Ala + Xyl -> AX-ARP (70/80/90 C) | **NEW** quantity; no counterpart in 2024 |
| k for Ala + GO -> pyrazine at 70/80/90 C | **NEW temperatures**, but discrepant with the 2024 ladder (flag 2); the 2024 100/110/120 C values are the B18 fit rows |
| Ea for Ala + GO -> pyrazine | second determination (104.97 vs 100.59 kJ/mol); slopes agree, intercepts do not |
| Ala + MGO -> 2,5-DMP rate | NOT repeated here (2024 only) |
| 13C5-Xyl isotopomer table | **NEW** (2024 had the 15N-Ala table) |
| 2,6-DMP and 2-furfural calibration lines; 2,6-DMP measured | **NEW** |
| pyrazine / methylpyrazine / 2,5-DMP calibration lines | identical coefficients to 2024 (same calibration reused) |
| AX-ARP-alone control at 20 mmol/L, pH 8, 120 C, 30-120 min (pyrazines, DXs, GO, MGO, ARP, Ala, pH) | same design as 2024's Ala-ARP arm; whether re-measured or re-plotted is not stated; figure-only in both papers (2024 Fig. 1 in µmol/L, here Fig. 1 in µg/L) |
| GO excess 0.36 mmol/L, Xyl-alone GO 0.12 mmol/L, 1.1 % Xyl conversion | **NEW** printed numbers |
| ARP preparation, OPD method, SPME method, 20 + 20 mmol/L / pH 8 / 120 C design, triplicates | repeat (IS for the OPD method — 3,4-hexanedione — is printed here and not in the 2024 dossier) |

## 5. Flags

1. **Reactant concentration of the kinetic runs is unstated.** Section 2.7 says "equimolar amounts";
   Figure 4's caption says "molar ratios ... 1:1". Nothing gives mmol/L. Zero-order k values are
   rates at a concentration, so these four numbers cannot be transported into a mass-action lane
   without it; 20 mmol/L by analogy with §2.3 and with 2024 is an inference, and flag 2 argues
   against it.
2. **The Ala + GO -> pyrazine rates are ~70x above Zhou 2024's for the nominally same reaction.**
   Same lab, same product, same SPME calibration line; Ea agrees within 4 kJ/mol; the intercepts
   differ by exp(5.72). If the true rate law were second order in [Ala][GO], a 70x offset needs
   ~8x higher concentrations here (~0.17 mol/L each); if first order, ~70x. Alternatively a
   µmol/mmol slip in one of the two papers, a different vessel/headspace, or a different fitted
   window (this paper fits the initial segment only and says the rate falls off; 2024 fitted 0-120
   min). The paper does not mention the 2024 rates. **B18's fit rows are 2024's; these must not be
   added to them, nor used to check them, until the concentration is known.** Recorded here as a
   discrepancy, not resolved.
3. **"Zero order" is an initial-rate statement** (as in 2024 flag 1). Here the authors themselves
   report curvature: pyrazine "deviation from zero-order kinetics ... under all three temperatures";
   AX-ARP peaks at 60 min at 90 C. The fitted window is not printed; from the rasterised Fig. 4 the
   fitted segments are shorter at higher temperature. The inference "the reaction rate ... was
   independent of the concentrations of substrates" does not follow from a single-concentration
   experiment.
4. **The two k rows are in different units** (mmol vs µmol per L·min). The text's "lower than that for
   the formation of AX-ARP" is correct only after conversion (factor 19-29, §3).
5. **Unbuffered, initial pH only.** ARP runs: pH 8 set with an unstated base; kinetic runs: 1 mol/L
   NaOH. pH followed only in the 120 C ARP runs and only as a figure (axis 4-8; both arms fall from
   8). The kinetic runs' working pH is unknown after t = 0; the Xyl-Ala run generates acids.
6. **ARP preparation ratio misprinted**: "1:1 M ratio" for 1.78 g Ala + 6 g Xyl, which is 20 : 40 mmol
   (1 : 2), the charge the 2024 dossier records.
7. **2,6-dimethylpyrazine is detected in the AX-ARP-alone model here (Fig. 1d) but "was not
   detected" in either ARP model in Zhou 2024** under the same design. Either a sensitivity difference
   (a 2,6-DMP calibration exists only here) or a different batch; not discussed.
8. **Calibration.** x unit not stated (µg/L in 2024; Fig. 1 in µg/L); single IS for five analytes;
   standards' matrix not stated; no LOD/LOQ, so "entirely absent" for 2,5-/2,6-DMP has no numeric
   bound; the 2-furfural line's intercept is ~33 µg/L-equivalent.
9. **Methods missing for Fig. 3a-c**: how AX-ARP, regenerated Ala and A420 were measured in the 120 C
   runs is not in this paper (§2.8 covers only the 70-90 C Ala-Xyl run). Presumably 2024's HPLC-ELSD
   and UPLC-TQD.
10. **Xyl-alone charge inferred** (20 mmol/L from "identical reaction conditions"); the 1.1 %
    conversion and the 0.12 mmol/L GO rest on it.
11. **Isotope table**: pyrazine M+4 at 30 min (1.96 %) is ~3x what a single well-mixed GO pool
    predicts from the M column (my arithmetic, §3); later points fit. Methylpyrazine M+3 (MGO from
    labelled xylose) appears only at 90 and 120 min. Natural-abundance and 98 atom % corrections
    applied by the authors; the raw intensities are not printed.
12. **Sampling starts at 30 min**; with ~80 % of the ARP gone by 30 min (Fig. 3a, figure-only), the
    early kinetics of the 120 C ARP runs are not resolved, as in 2024.
13. **Commercial 40 % glyoxal**, nominal, oligomeric hydrates; no assay (as 2024 flag 5).
14. **No SI; data "available on request".** Triplicates, mean ± SD, error bars drawn; no confidence
    intervals on k or Ea.
15. Registry: 2,6-dimethylpyrazine has a registry id but no `structures.yml` key; xylose is carried
    only as the generic open-chain aldopentose `PENT`; 13C isotopomers have no representation (none
    needed for a rule).
