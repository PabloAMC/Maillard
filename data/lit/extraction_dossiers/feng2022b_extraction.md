# Feng, Yang, Xie, Liu, Peng, Hu & Yu 2023 — EXTRACTION (file stem `feng2022b`)
DOI: 10.3389/fnut.2022.1022254 (printed on p. 1 of every page header; Crossref 2026-11-01: "Kinetics of L-ascorbic
acid degradation and non-enzymatic browning development in hot-compressed water", Front. Nutr. **9**:1022254, authors
Feng, Yang, Xie, Liu, Peng, Hu, Yu; published online / issued **2023-01-12**. **Matches the brief.**)
[!] **Year: the file stem says 2022b, the paper is 2023** (received 20 Aug 2022, accepted 20 Dec 2022, published
12 Jan 2023; volume 9 is the journal's 2022 volume, hence the "2022" in the DOI). The paper's own citation line reads
"(2023)". The stem is kept so the coverage gate matches.
**Not to be confused with** `feng2022_extraction.md` (an unrelated Jiangnan paper). **Same laboratory as**
`feng2023_extraction.md` (Feng et al. 2024 JSFA, ASA + Cys volatiles in the same autoclave).

### L-ascorbic acid (ASA) alone, 0.03 / 0.05 / 0.07 M in 10 mL 0.20 M phosphate, initial pH 5.0 / 7.0 / 9.5 (NaOH), N2-purged sealed 316L autoclave with quartz liner heated in an oven, 150 / 170 / 180 / 190 °C, 0-150 min; residual ASA by HPLC-DAD (243 nm, external standard), A294 and A420. Pseudo-first-order k with ± for 3 pH × 3 [ASA]0 × 4 T (Table 1) and an "average apparent Ea" per pH (15.77 / 31.70 / 47.53 kJ/mol). **No dicarbonyl, no DHA, no furfural, no pH after heating, no 160 °C, no time courses at 170 / 180 °C.**

**Read 2026-11-01 from `data/articles/feng2022b.pdf`** (10 pp., Frontiers, born-digital LaTeX). Prose from
`pdftotext -layout`; **Table 1 (rotated, p. 5) read from a 170 dpi render** (scratchpad
`ingest/feng2022b_work/p5-05.png`; the text layer scrambles the rotated cells). **Fig. 2 (p. 4) and Fig. 3 (p. 6)
read from 400 dpi crops** (`fig2hi-04.png`, `fig3top-06.png`, `fig3bot-06.png`), read error about ±0.002 M on
[ASA] and ±1 absorbance unit on the A294/A420 axis (0-60). Fig. 1 A420 values are printed on the figure [M].
Arrhenius re-fits of Table 1 by `awk` on `feng2022b_work/table1_k.tsv` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** School of Chemistry and Environmental Engineering, **Hubei Minzu University**, Enshi (and the Key
Laboratory of Biologic Resources Protection and Utilization of Hubei Province); corresponding author **Yan Yang**,
senior author **Ai-nong Yu**; Peng X at the Hubei Provincial Key Laboratory of Rheumatic Diseases (same university).
**Not in `data/keys/laboratories.yml`: NEW**, proposed id **`hubei_minzu_yu_an`**, lab string
"Hubei Minzu Univ. (Yang Y / Yu AN)" — the same proposal as `feng2023_extraction.md`. Evidence: identical
institution, corresponding and senior authors, first author, autoclave (BE100-HP, 316L, quartz liner, N2 purge) and
ASA chemistry as Feng 2024 (`feng2023`); same senior author and ASA pots as `yu2016_extraction.md` ("Hubei Univ. for
Nationalities (Yu AN)", the university's former name) and `yu2010_extraction.md`. The "previous work" whose Ea
the paper prints as Ea* (ref. 9, Li, Yang, Yu & Wang 2016 Food Sci Biotechnol 25:97, 10.1007/s10068-016-0014-x) and
the method reference (ref. 4, Feng et al. 2022 Food Chem 131155, ASA self-degradation volatiles in HCW) are by the
same group; neither is on disk.

## 1. What this paper answers for the model

| model need | what Feng 2023 gives | tag |
|---|---|---|
| ASA loss rate vs T, pH, [ASA]0 (anaerobic, aqueous) | Table 1: k (min⁻¹) ± and R² for 36 cells, 150-190 °C | [M] |
| Ea of ASA loss by pH | 15.77 ± 1.08 / 31.7 ± 2.41 / 47.53 ± 8.14 kJ/mol (pH 5.0 / 7.0 / 9.5); ± not defined | [M] |
| pre-exponential factors | **not printed**; per-cell ln A from my re-fit (sec. 4) | [Z] |
| ASA time courses | Fig. 2 / Fig. 3: 150 and 190 °C only, 8 time points | [F] |
| dicarbonyls (DHA, DKG, 3-DP, C4/C5 osones) from ASA | **none measured**; mechanism drawn via aldopentose → furfural | [P] |
| browning (A420) and UV intermediates (A294) | Fig. 1 (190 °C, 90 min, printed A420), Fig. 3 time courses, text endpoints; **no fitted browning k** | [M] / [F] |
| browning rate law | Eq. 11, d[BP]/dt ∝ [H3O+][ASA]T, derived by steady state; no constant fitted | [P] |
| pH during / after heating | **not measured** | [!] |
| O2 | N2-purged, "oxygen-free" [Q]; anaerobic route only | [M] |

**Answers.**
1. **Are the repo's three Ea values in the paper? Yes, as printed** (confidence ~99 %): 15.77, 31.70, 47.53 kJ/mol at
   pH 5.0, 7.0, 9.5 (abstract, Table 1, p. 5 text). **But only the pH 5.0 value reproduces from the paper's own
   Table 1** (my mean of three per-concentration Arrhenius slopes = 15.77; pH 7.0 gives 29.50, pH 9.5 gives 42.17;
   pH 9.5's printed 47.53 equals the 0.07 M slope alone, 47.51) [Z]. Possibly the authors used data not in Table 1
   (the abstract's 150-190 °C may include a 160 °C run the table omits); not resolvable from the paper.
2. **Are these Ea usable as an ASA-degradation activation energy? Poorly** (confidence ~75 % that they understate
   the true chemical Ea). Reasons: (a) measured [ASA] **rises 6-27 % above the initial value at 30 min** in all
   150 °C arms (e.g. 0.050 → 0.062 M, pH 7.0) and in several 190 °C arms, outside the error bars — impossible for a
   reactant with no source, so either the HPLC peak picks up a 243 nm co-eluting product or the t = 0 basis is off;
   (b) at 190 °C the ASA hardly moves in the first 30 min though k(190) would predict ~40 % loss, i.e. the
   oven-heated autoclave was probably not at temperature for much of 0-30 min (heat-up not reported); (c) the
   printed k reproduce from the 30-150 min window, not from t = 0 [Z]; (d) within that window the decay
   decelerates 2.5-3 × (not first-order; the authors say ASA "reached equilibrium" [Q]); (e) k rises only 1.4-1.6 ×
   across 150 → 190 °C at pH 5, an Ea of 12-18 kJ/mol, below anything expected for covalent chemistry and far below
   the same lab's ASA loss Ea at 125-155 °C (Yu & Tang 2016, ASA + Cys, pH 8: 98-109 kJ/mol, `yu2016_extraction.md`)
   and the group's own earlier ambient-pressure values printed here as Ea* (77.79 / 46.88 / 93.33 kJ/mol).
3. **Direction of pH: acid faster — yes, at 150 °C; nearly gone at 190 °C** (confidence ~80 % for the 150 °C
   ordering). k(pH 5)/k(pH 9.5) = 1.8-2.3 at 150 °C but 1.11-1.16 at 190 °C [Z from Table 1]. Because the low Ea is
   at low pH, the acid advantage *shrinks* as T rises. "Acidic high-temperature conditions accelerate" is not what
   Table 1 shows at the high-temperature end.
4. **Dicarbonyl release: not shown** (confidence ~99 %). No α-dicarbonyl, DHA, DKG, furfural or any product was
   quantified; the browning precursor pool is A294 "uncoloured intermediate products" (UIPs). The mechanism (Fig. 4:
   lactone hydrolysis → (Z)-pentahydroxyhexenoic acid → decarboxylation → aldopentose → furfural → furan polymers)
   is [P], resting on volatiles from ref. 4.
5. **Browning kinetics: no constants** (confidence ~99 %). A420 at 190 °C / 90 min is printed (Fig. 1); the rate law
   (Eq. 11) predicts browning ∝ [H3O+], i.e. a 10^4.5 span from pH 9.5 to 5.0, whereas Fig. 1 shows a ×1.22-1.33
   span [Z]. The model is contradicted by the paper's own data.

## 2. Conditions

- **Reactant** [M]: L-ascorbic acid (≥ 99.7 %, Sigma-Aldrich), **0.30 / 0.50 / 0.70 mmol in 10 mL** → 0.03 / 0.05 /
  0.07 M. ASA is the only reactant ("only ASA participated in the reaction as the material" [V]).
- **Buffer / pH** [M]: 0.20 M sodium phosphate (Na2HPO4 / NaH2PO4); pH adjusted to 5.0, 7.0, 9.5 with NaOH
  (Sartorius meter). **Initial room-temperature pH only; no pH during or after heating.** [!] Phosphate buffers
  between ~pH 6.2 and 8.2 (pKa2 7.2): the **pH 5.0 and pH 9.5 arms are weakly buffered** [Z], and 0.07 M ASA plus its
  acid products is a third of the phosphate concentration. The authors treat [H3O+] as 10^-5.0 to 10^-9.5 at reaction
  temperature (p. 7), ignoring the temperature shifts of both pKw and phosphate pKa2.
- **Vessel** [M]: BE100-HP high-temperature high-pressure autoclave (Shanghai LABE), 316L stainless steel, quartz
  lining; 10 mL charge (vessel volume not printed; "BE100" suggests ~100 mL, i.e. ~90 mL N2 headspace [Z,
  unverified]).
- **Headspace / O2** [M]: N2 flowed through the inlet to expel air from vessel and solution; "oxygen-free" [Q].
  Anaerobic ASA degradation only.
- **Heating** [M]: **in an oven** ("reaction temperatures were controlled using the oven" [V]); pressure = saturated
  vapour pressure, 0.791-1.254 MPa. **Heat-up time not reported** [!]. Quench in ice-bath water. One run per time
  point presumably (terminated "when each reaction time was reached").
- **T / t** [M]: 150, 170, 180, 190 °C in Table 1 (**no 160 °C**, though the text says 150-190 °C); time 30-150 min
  (abstract: 0.50, 0.75, 1.00, 1.25, 1.50, 2.00, 2.50 h; figures add t = 0). Time courses shown for 150 and 190 °C
  only.
- **Replicates** [M]: "at least in triplicate"; RSD "lower than 5 %" — repeated until the RSD target was met [V/Q]
  [!].

## 3. Analytical method

- **ASA** [M]: RP-HPLC (Agilent 1260, DAD), Agilent HC-C18 3.5 µm 4.6 × 100 mm, 30 °C, 1.0 mL/min, 10 µL, **243 nm**;
  mobile phase water / methanol / 0.1 wt % metaphosphoric acid, gradient 8/8/84 → 50/50/0 (20 min) → 8/8/84.
  External standard, Y = 30.10757 X + 6.64471, R² = 0.9998 (X units not printed). No LOD, no peak-purity check, no
  resolution from degradation products reported [!]. DHA (no 243 nm absorbance) would be invisible.
- **A294 (UIPs) and A420 (BPs)** [M]: Shimadzu UV-2550 after dilution "30x-100x / 0x-25x" (so A < 1.0); values
  reported dilution-corrected (up to ~56). BPs were "measured in a neutral environment because the solution has been
  greatly diluted" [Q].
- **Kinetics** [M]: Formula A, −(1/n) d[ASA]/dt = k[ASA]^n with n = 1 chosen; Arrhenius from ln k vs 1/T; Origin 8.0.
  The meaning of the ± on k and on Ea is **not stated** (fit SE or replicate SD).

## 4. Data

### 4a. Table 1 (p. 5): pseudo-first-order k (min⁻¹), ± (undefined), (R²) [M]

| pH | [ASA]0 | 150 °C | 170 °C | 180 °C | 190 °C |
|---|---|---|---|---|---|
| 5.0 | 0.03 M | 0.01279 ± 0.00049 (0.9922) | 0.01466 ± 0.00093 (0.987) | 0.01523 ± 0.00077 (0.9782) | 0.01768 ± 0.00020 (0.9734) |
| 5.0 | 0.05 M | 0.01049 ± 0.00014 (0.9777) | 0.01218 ± 0.00094 (0.8688) | 0.01329 ± 0.00087 (0.9089) | 0.01667 ± 0.00049 (0.9187) |
| 5.0 | 0.07 M | 0.0101 ± 0.00039 (0.9281) | 0.01175 ± 0.00049 (0.9427) | 0.01303 ± 0.00098 (0.9611) | 0.01594 ± 0.00078 (0.925) |
| 7.0 | 0.03 M | 0.0085 ± 0.00054 (0.9693) | 0.01378 ± 0.00068 (0.9679) | 0.01436 ± 0.00014 (0.9505) | 0.01739 ± 0.00061 (0.8943) |
| 7.0 | 0.05 M | 0.00736 ± 0.00039 (0.996) | 0.01181 ± 0.00099 (0.9925) | 0.01235 ± 0.00027 (0.9886) | 0.01571 ± 0.00076 (0.853) |
| 7.0 | 0.07 M | 0.0069 ± 0.00040 (0.9462) | 0.01089 ± 0.00092 (0.9574) | 0.0118 ± 0.00016 (0.9717) | 0.0149 ± 0.00070 (0.9542) |
| 9.5 | 0.03 M | 0.00697 ± 0.00060 (0.9992) | 0.01146 ± 0.00041 (0.9894) | 0.01361 ± 0.00023 (0.9817) | 0.01587 ± 0.00099 (0.9961) |
| 9.5 | 0.05 M | 0.0049 ± 0.00094 (0.9945) | 0.01031 ± 0.00017 (0.9867) | 0.0121 ± 0.00091 (0.9824) | 0.01478 ± 0.00015 (0.9901) |
| 9.5 | 0.07 M | 0.00439 ± 0.00077 (0.9925) | 0.00919 ± 0.00026 (0.9888) | 0.01162 ± 0.00041 (0.9904) | 0.0138 ± 0.00038 (0.9773) |

Ea (kJ/mol, this work) [M]: **pH 5.0 15.77 ± 1.08; pH 7.0 31.7 ± 2.41; pH 9.5 47.53 ± 8.14.**
Ea* ("previous work (9)", ambient-pressure ASA self-degradation, Li et al. 2016) [M/C]: 77.79 / 46.88 / 93.33.
Pre-exponential factors: **not printed**.

### 4b. My Arrhenius re-fit of Table 1 (4 temperatures per row, unweighted ln k vs 1/T) [Z]

| pH | [ASA]0 | Ea (kJ/mol) | ln A (A in min⁻¹) | A (min⁻¹) |
|---|---|---|---|---|
| 5.0 | 0.03 / 0.05 / 0.07 | 12.25 / 17.51 / 17.54 | −0.890 / 0.385 / 0.360 | 0.41 / 1.47 / 1.43 |
| 7.0 | 0.03 / 0.05 / 0.07 | 28.41 / 29.69 / 30.39 | 3.339 / 3.548 / 3.676 | 28 / 35 / 40 |
| 9.5 | 0.03 / 0.05 / 0.07 | 33.82 / 45.19 / 47.51 | 4.670 / 7.579 / 8.123 | 107 / 1.96e3 / 3.37e3 |

Mean ± SD of the three per-row Ea: pH 5.0 **15.77 ± 3.04** (matches the printed mean; the printed ± 1.08 matches
neither the SD nor the SE 1.75); pH 7.0 **29.50 ± 1.01** (printed 31.7); pH 9.5 **42.17 ± 7.33** (printed 47.53 =
the 0.07 M row). Pooled 12-point fits give the same 15.77 / 29.50 / 42.17.

Derived contrasts [Z]: k(pH 5.0)/k(pH 9.5) at 150 °C = 1.84 / 2.14 / 2.30 (0.03 / 0.05 / 0.07 M); at 190 °C = 1.11 /
1.13 / 1.16. k(190)/k(150) = 1.38-1.59 at pH 5.0, 2.05-2.16 at pH 7.0, 2.28-3.14 at pH 9.5. k falls 10-37 % from
0.03 to 0.07 M at fixed T, pH (the authors' "excessive concentration ... no obvious effect" [V] vs "k decreased
slightly" [Q], p. 5).

### 4c. Fig. 2 = Fig. 3 [ASA] time courses (M), 150 and 190 °C, read ±0.002 M [F]

pH 5.0 (Fig. 2A):

| arm | 0 | 30 | 45 | 60 | 75 | 90 | 120 | 150 min |
|---|---|---|---|---|---|---|---|---|
| 150 °C, 0.07 | 0.070 | 0.079 | 0.061 | 0.048 | 0.039 | 0.032 | 0.030 | 0.022 |
| 190 °C, 0.07 | 0.070 | 0.071 | 0.059 | 0.035 | 0.025 | 0.018 | 0.015 | 0.011 |
| 150 °C, 0.05 | 0.050 | 0.056 | 0.053 | 0.044 | 0.034 | 0.028 | 0.023 | 0.017 |
| 190 °C, 0.05 | 0.050 | 0.052 | 0.041 | 0.022 | 0.017 | 0.012 | 0.009 | 0.007 |
| 150 °C, 0.03 | 0.030 | 0.038 | 0.036 | 0.029 | 0.022 | 0.018 | 0.014 | 0.010 |
| 190 °C, 0.03 | 0.030 | 0.036 | 0.031 | 0.021 | 0.014 | 0.009 | 0.008 | 0.005 |

pH 7.0 (Fig. 2B), 30 min / 150 min: 150 °C 0.074 / 0.032 (0.07), 0.062 / 0.026 (0.05), 0.038 / 0.014 (0.03);
190 °C 0.071 / 0.012, 0.050 / 0.009, 0.033 / 0.005.
pH 9.5 (Fig. 2C), 30 min / 150 min: 150 °C 0.078 / 0.046, 0.055 / 0.030, 0.032 / 0.014; 190 °C 0.071 / 0.014,
0.051 / 0.009, 0.031 / 0.005.
Text [M]: at pH 5.0 / 190 °C, "degradation rate" (fraction lost by the end) 84.4 / 86.2 / 87.3 % for 0.07 / 0.05 /
0.03 M (Fig. 2A endpoints give 85 / 86 / 84 % [Z]).

Check of the k basis [Z]: ln(C30/C150)/120 min reproduces Table 1 (pH 5.0 / 190 °C / 0.07: 0.0156 vs 0.01594;
pH 9.5 / 150 °C / 0.05: 0.0050 vs 0.0049; pH 7.0 / 150 °C / 0.05: 0.0074 vs 0.00736); ln(C0/C150)/150 does not
(0.0125, 0.0034, 0.0045). So the fits apparently start at 30 min (~80 % confidence).

### 4d. Browning and UV intermediates

**Fig. 1 (p. 4), A420 after 190 °C, 90 min** [M] (dilution-corrected):

| pH | 0.03 M | 0.05 M | 0.07 M |
|---|---|---|---|
| 5.0 | 2.70 | 5.28 | 7.83 |
| 7.0 | 2.40 | 5.03 | 6.90 |
| 9.5 | 2.03 | 4.32 | 6.15 |

**Text p. 4, pH 5.0, A420 at 150 min** [M]: 150 °C 1.58 / 3.15 / 4.11; 190 °C 4.10 / 6.40 / 10.33 (order of the
three concentrations not stated; Fig. 3 legend implies 0.03 / 0.05 / 0.07 M).
**Fig. 3, A294 at 150 min** [F] (±1): pH 5.0 150 °C 15 / 21 / 31; 190 °C 30 / 48 / 56 (0.03 / 0.05 / 0.07 M);
pH 7.0 190 °C 15 / 35 / 42; pH 9.5 190 °C 12 / 27 / 36. A420 time courses rise roughly linearly after a ~30 min lag;
A294 rises sigmoidally at 190 °C and keeps rising after ASA has levelled off. The 150 °C pH 7.0 / 9.5 A420 and A294
series are figure-only (not digitised).
**No browning rate constant, order or Ea is printed**; Eq. 11 is a derived proportionality.

## 5. Mechanism evidence and its kind

Within-study factor contrasts (T × pH × [ASA]0), time-resolved reactant loss, anaerobic, buffered at the start.
Products: none identified in this paper; the Fig. 4 pathways and Scheme 1 (ASA → aldopentose → UIPs → BPs, plus a
"caramelisation" step from the pentose in acid) are [P], drawn from ref. 4 and Kurata & Sakurai 1967 [C]. The
kinetic model (Eqs. 1-11) treats water as a reactant and H3O+ as a catalyst; it fits no constant.

## 6. Check of the repo's existing citations

Records checked: `data/lit/computational_priors.json` → `ascorbic_pathway_priors` → `frontiers_2022_hcw_aa_arrhenius_v1`;
`data/lit/benchmark_intake_registry.json` → `frontiers_2022_hcw_aa_arrhenius_anchor`; and, found on the way,
the **twin record of the same id in `computational_priors.json`** (`frontiers_2022_hcw_aa_arrhenius_anchor`,
"Ea_hcw_ascorbic_kj_mol": 82.0 at T_C 160, confidence "high"), which `results/validation/key_value_trace_report.md`
attributes to this paper. Also `scripts/generators/generate_gap_heatmap.py` l. 234 (Family 14 card text).

### 6a. `computational_priors.json` → `frontiers_2022_hcw_aa_arrhenius_v1`

| field | repo value | paper | verdict |
|---|---|---|---|
| source | "Frontiers in Nutrition 2022, DOI:10.3389/fnut.2022.1022254" | Front. Nutr. 9:1022254, published 12 Jan **2023**; Feng et al. | **wrong year** (2023; vol. 9 is the 2022 volume); journal and DOI match; no author named |
| doi | 10.3389/fnut.2022.1022254 | same | **matches** |
| provenance_tier | direct_measurement | Ea is the authors' Arrhenius fit of their own pseudo-first-order k | **overstated**: a fitted secondary quantity, and pH 7.0 / 9.5 do not reproduce from Table 1 (sec. 4b) |
| confidence_tier | medium_high | heat-up / supra-initial [ASA] / non-first-order artefacts (sec. 1, ans. 2) | **wrong** — `low` is defensible |
| matrix_context | hot_compressed_water | 0.20 M phosphate in a sealed N2-purged autoclave at saturated vapour pressure | **matches** (record omits: phosphate, anaerobic) |
| temperature_window_c | [150, 190] | 150-190 °C; Table 1 has 150 / 170 / 180 / 190 (no 160) | **matches** |
| kinetic_order | pseudo_first_order | authors' choice, n = 1; R² 0.853-0.999; fits start at 30 min; decay decelerates | **matches the authors' claim** (data only loosely first-order) |
| ea 5.0 | 15.77 | 15.77 ± 1.08 | **matches** (reproduces from Table 1) |
| ea 7.0 | 31.7 | 31.7 ± 2.41 | **matches as printed** (Table 1 gives 29.5 [Z]) |
| ea 9.5 | 47.53 | 47.53 ± 8.14 | **matches as printed** (Table 1 gives 42.2; 47.5 = 0.07 M row only [Z]) |
| recommended_interpolation | linear_between_pH_5_and_9_5 | three discrete pH values; no interpolation proposed | **not in paper** — the repo's choice (Ea is roughly linear in pH, 8.0 then 6.3 kJ/mol per pH unit [Z]; but that is linearity of an artefact-laden Ea) |
| effect_direction (acid faster) | "acidic … accelerate ascorbic breakdown" | k(pH 5) > k(pH 9.5) in every cell; authors say acid favours degradation | **matches** |
| effect_direction ("high-temperature" × acid) | "acidic high-temperature conditions accelerate" | acid advantage ×1.8-2.3 at 150 °C, ×1.1-1.2 at 190 °C | **misleading**: the pH effect fades as T rises |
| effect_direction ("dicarbonyl release") | "…and dicarbonyl release" | no dicarbonyl measured; mechanism via aldopentose / furfural [P] | **not in paper** |
| notes ("Primary Family 14 Arrhenius anchor for extrusion-like AA degradation") | — | dilute aqueous ASA (99 % water), 30-150 min | **not in paper**; the transfer to extrusion is ours and crosses a water-activity gap |
| notes ("keep ascorbic-acid dicarbonyl release explicitly pH-sensitive") | — | no dicarbonyl | **not in paper** |
| doi_propagation.crossref_verified (title, authors, year 2023, container) | — | Crossref 2026-11-01 agrees | **matches** (and contradicts the record's own "2022" source string) |

### 6b. `benchmark_intake_registry.json` → `frontiers_2022_hcw_aa_arrhenius_anchor`

| field | repo value | paper | verdict |
|---|---|---|---|
| citation | Feng et al. (2023) | Feng et al. 2023 | **matches** |
| doi | 10.3389/fnut.2022.1022254 | same | **matches** |
| chemistry_family / kind | ascorbic_acid_maillard / calibration_reference | ASA alone, no amine: ASA self-degradation and browning, not a Maillard pot | **imprecise** (no amino compound) |
| observable_panel_tags | arrhenius_anchor, ascorbic_support, **dicarbonyl** | no dicarbonyl | **dicarbonyl: not in paper**; others fine |
| process_state_scope | heated_matrix, aqueous_pre_extrusion_model, extrusion_structured | dilute aqueous buffer, sealed autoclave | **extrusion_structured: not in paper** |
| what_it_supports [0] | pH-sensitive AA degradation "under extrusion-like thermal conditions" | pH-sensitive at 150 °C, weakly at 190 °C; temperatures overlap extrusion, moisture does not | **partly** |
| what_it_supports [1] | "Explicit HCW calibration for AA-driven dicarbonyl release" | no dicarbonyl | **not in paper** |
| key_values.ea_kj_per_mol_by_ph | 15.77 / 31.7 / 47.53 | same | **matches as printed** (see 6a) |
| key_values.kinetic_order | pseudo_first_order | same | **matches** |
| key_values.temperature_window_c | [150, 190] | same | **matches** |
| eligibility / status | runtime_eligible / encoded | — | a repo decision; given sec. 1 ans. 2, not justified |
| requires_primary_data | false | — | repo decision |

### 6c. `computational_priors.json` → `frontiers_2022_hcw_aa_arrhenius_anchor` (twin, not named in the brief)

| field | repo value | paper | verdict |
|---|---|---|---|
| source | "Frontiers (2022) HCW AA Arrhenius anchor" | — | **wrong year**, no author, no DOI |
| provenance_tier / confidence_tier | literature_derived_transfer / **high** | — | **not supported** |
| reference_conditions.T_C | 160 | no 160 °C run in Table 1 or figures | **not in paper** |
| Ea_hcw_ascorbic_kj_mol | **82.0** | printed Ea: 15.77 / 31.7 / 47.53 (this work), 77.79 / 46.88 / 93.33 (Ea*, previous work, ambient pressure) | **not in paper** — no 82 anywhere; nearest is Ea*(pH 5.0) 77.79, which is not an HCW value |
| notes ("HCW AA degradation: Ea=82 kJ/mol") | — | — | **not in paper** |

### 6d. `scripts/generators/generate_gap_heatmap.py` l. 234 (Family 14 card)

"Dehydroascorbic Acid (DHA) thermal degradation: pseudo-first-order DHA ring-opening hydration to
2,3-diketogulonic acid (Ea 15.77 / 31.7 / 47.53 …; Frontiers Nutr. 2022)" — **wrong compound and wrong mechanism**:
the paper measured reduced L-ascorbic acid by HPLC at 243 nm (DHA does not absorb there), under N2; DHA and
2,3-diketogulonic acid are never mentioned. Year wrong as above.

### 6e. Who reads these records (grep 2026-11-01)

- `src/kinetic_core/`: **no reference** to ascorbic acid, `ascorbic_pathway_priors`, `ea_kj_per_mol_by_ph` or
  `Ea_hcw_ascorbic_kj_mol`. No fit generator reads them.
- `src/matrix_prior_registry.py` loads the whole priors file; `summarize_family_prior_bundle` copies only
  source / tier / tag metadata and has no caller outside the module.
- `src/deep_research_runtime_queue.py` reads every priors record's `source` string into the set of "landed
  runtime citations" (so editing the string changes that bookkeeping); `tests/scientific/test_literature_backlog.py`
  asserts `frontiers_2022_hcw_aa_arrhenius_anchor` is "encoded_runtime_artifact".
- `scripts/generators/generate_gap_heatmap.py` hard-codes the three Ea in report text (6d).
- `results/validation/key_value_trace_report.md` and `citation_verification_ledger.md` (the latter calls the DOI a
  "MATCH … as cited", which held for identity only).
- `data/keys/papers.yml` lists only `_v1` under `computational_priors.json` record_ids, not the 82 kJ/mol twin.

## 7. Bearing on open questions

- **No conflicts-register entry concerns ascorbic acid** (`grep -i "ascorb|vitamin c|dehydroascorb"` of
  `docs/reference/conflicts_register.yml`: no hits). The engine has no ascorbic-acid species, so nothing here is
  simulable today.
- **CF-054 / CF-084-style pH drift** — another NaOH-adjusted, initial-pH-only pot, here with 0.2 M phosphate that
  barely buffers at pH 5.0 and 9.5. Example of the problem, not a datum. P(helps) ~0 %.
- **Any future ASA carbonyl lane** — the Table 1 k at 150 °C (0.0044-0.0128 min⁻¹, anaerobic, phosphate) are the
  most usable numbers, as a weak upper / lower bracket on ASA loss; the Ea should not be carried. P(changes a
  shipped parameter) ~0 % (nothing reads them); P(the citation fix matters for reporting) ~90 %.
- **Agcam 2022 (`agcam2022_extraction.md`) "Q10 flattening of the most reactive precursor"** — same signature
  (ASA, flattest temperature response); a measurement-limited or heat-transfer-limited regime is a shared candidate
  cause.

## 8. Disagreements found

1. **vs Yu & Tang 2016 (same lab; ASA + Cys 0.10 M, 0.2 M phosphate pH 8, 125-155 °C, oil-bath vials)**: ASA loss
   k(155 °C) ≈ 0.0117 min⁻¹ and Ea ≈ 98-109 kJ/mol (digitised, `yu2016_extraction.md`). Feng 2023 at 150 °C, pH 7.0 /
   9.5, ASA alone: k = 0.0044-0.0085 (similar magnitude) but Ea 30-48. Extrapolating Yu's Ea to 190 °C gives
   k ≈ 0.1 min⁻¹ [Z], ~6 × Feng's 0.014-0.016. Candidate causes: oven-heated autoclave heat-up (Feng) vs oil-bath
   vials (Yu); Cys present (Yu); ASA by A243 total UV (Yu) vs HPLC (Feng); fit window.
2. **vs the group's own ambient-pressure Ea* (Li et al. 2016, printed in Table 1)**: 77.79 / 46.88 / 93.33 vs 15.77 /
   31.7 / 47.53. The authors attribute the drop to water participating in HCW [Q]; no test of that is offered.
3. **Ea vs pH ordering**: this paper has Ea rising monotonically 5.0 → 9.5; the earlier work has pH 7.0 lowest.

## 9. Flags

1. [!] **Supra-initial [ASA] at 30 min** in all nine 150 °C arms (+6-27 %) and some 190 °C arms; error bars
   ~±0.002 M do not cover it. Either co-elution at 243 nm or a t = 0 basis problem; either way the HPLC series is
   not a clean reactant trajectory.
2. [!] **Heat-up not reported; oven-heated steel autoclave**; ASA nearly constant for 0-30 min at 190 °C. The k are
   fitted from 30 min [Z]; the reaction time axis is not isothermal time.
3. [!] **Ea for pH 7.0 and 9.5 do not reproduce from Table 1** (29.5 and 42.2 vs 31.7 and 47.53); the Ea "±" matches
   neither SD nor SE; the ± on k is undefined.
4. [!] **Table 1 has no 160 °C column** although the methods say 150-190 °C; 170 / 180 °C time courses not shown.
5. [!] **pH initial only**; phosphate weak at pH 5.0 and 9.5; [H3O+] taken as 10^-pH(25 °C) at 150-190 °C.
6. [!] **Browning model contradicts the data**: Eq. 11 makes the browning rate ∝ [H3O+] (10^4.5 span) yet Fig. 1
   spans ×1.3; the text says browning is "independent of the concentration of ASA" [V] while Eq. 11 is ∝ [ASA]T
   and Fig. 1 A420 rises ~×2.9 from 0.03 to 0.07 M.
7. [!] Text p. 5 sentence on A294 / A420 vs pH is garbled (can be read either way); Fig. 1 has pH 5.0 highest.
8. [!] Abstract and conclusion: ASA degradation "more favorable" at the lower concentration — true of k (−10-37 %),
   while the authors also say concentration had "no obvious effect" on k (p. 5).
9. [!] Replicates run until RSD < 5 % (selection to a target); dilution factors for A294/A420 given as ranges.
10. [!] File stem 2022b, paper 2023.

## 10. Ledger-ready measurements

- Table 1 k (min⁻¹) ± and R², 3 pH × 3 [ASA]0 × 4 T, plus the printed Ea — trunk ledger, rate_constant,
  calibrated-ish (external-standard HPLC) but artefact-flagged; not simulable (no ASA species).
- Fig. 2 / Fig. 3 [ASA] time courses (150, 190 °C; 3 pH × 3 [ASA]0; 8 times), digitised for pH 5.0 and as
  30 / 150 min pairs for pH 7.0 / 9.5 — trunk ledger, time_series.
- Fig. 1 A420 (190 °C, 90 min, printed) + text A420 endpoints + Fig. 3 A294 / A420 curves — trunk ledger, browning,
  quant_class other.
- Not keyable: any dicarbonyl, DHA, furfural, pH trajectory, browning rate constant (none measured or fitted).
