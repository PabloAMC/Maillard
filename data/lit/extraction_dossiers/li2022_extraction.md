# Li, Paulsen, Akıllıoğlu, Nielsen, Engholm-Keller & Lund 2022 — EXTRACTION
### *"Cysteine residues are responsible for the sulfurous off-flavor formed in heated whey protein solutions"*
### H2S by static-headspace GC-FPD from heated whey proteins and from FREE L-CYSTEINE (5.7 mM, 10 mM phosphate pH 7.0, 1 mL in a sealed crimp vial) at 90 °C / 10 min, 90 °C / 120 min and 160 °C / 160 s. The free-cysteine arm is the only sealed-vial, buffered, low-concentration cysteine → H2S measurement in the corpus. At 10 min it matches the engine's Zheng & Ho constant (1.1-1.5×); at 120 min it is 3-4× BELOW it. Read 2026-10-29.

**Full read of `data/articles/li2022.pdf`** (9 pp., born-digital, text layer clean; Table 1 checked
row by row against the layout text). Supporting Information (Figs. S1-S2, Table S1) is **not on disk**.

**Provenance:** **[M]** printed value · **[Z]** derived by me from printed values, never printed ·
**[!]** flag.

---

## 0. IDENTITY

| field | value as printed |
|---|---|
| Authors | Chengkang Li, Peter A. Paulsen, Halise Gül Akıllıoğlu, Søren B. Nielsen, Kasper Engholm-Keller, Marianne N. Lund |
| Affiliations | Dept. of Food Science, Univ. Copenhagen; Arla Foods Ingredients, Videbæk; Dept. of Biomedical Sciences, Univ. Copenhagen |
| Venue | *Food Chemistry: Molecular Sciences* **5** (2022) 100120 (open access, CC BY) |
| DOI (printed p. 1) | **10.1016/j.fochms.2022.100120** |
| Dates | received 3 April 2022, revised 30 May 2022, accepted 9 July 2022, online 12 July 2022 |
| Competing interests | two authors employed by Arla Foods Ingredients (declared) |

---

## 1. THE ONE-PARAGRAPH ANSWER

H2S is the only sulfur volatile GC-FPD detects in any of the heated whey-protein or amino-acid
solutions. For **free L-cysteine, 5.7 mM** in 10 mM phosphate pH 7.0, 1 mL sealed: **3.8 ± 0.3 µM
(90 °C, 10 min; n ≥ 2), 10.7 µM (90 °C, 120 min; n = 1), 44.8 µM (160 °C, 160 s; n = 1)**. The Zheng &
Ho first-order constant (pH 7, 90 °C) predicts 3.4 µM (Table I k) or 2.6 µM (the engine's Arrhenius
refit) at 10 min, and **41 or 31 µM at 120 min**: the 10-min point agrees, the 120-min point is
**0.26-0.34× the no-sink prediction**. Because the vial is sealed and sampled cold against an external
Na2S standard curve, headspace partitioning at 90 °C cannot remove sulfide from the count, and a
first-order storage loss would scale both time points alike. The 120-min deficit therefore points to
a **time-dependent loss of H2S or of its precursor during heating** (air headspace, so cysteine
oxidation and sulfide oxidation are both open), with the caveat that the 120-min value is a single
determination.

---

## 2. SYSTEM AND CONDITIONS [M], quoted

- **Free-cysteine arm (sec. 2.2):** "α-LA (2.45 mg/mL), β-LG (15.88 mg/mL), a whey model (2.45 mg/mL
  α-LA and 15.88 mg/mL β-LG), Cys (0.69 mg/mL), Met (0.54 mg/mL), sulfur amino acid mixture
  (0.69 mg/mL Cys and 0.54 mg/mL Met) and 3 % WPI (w/w) solutions were prepared in 10 mM phosphate
  buffer (pH 7.0) for GC-FPD analysis." Table 1 footnote converts: "Cys: 0.69 mg/mL (5.7 mM); Met:
  0.54 mg/mL (3.6 mM)". L-Cys purity "(97 %)"; supplier Merck.
- **Concentrations chosen** "to reflect their concentrations in a 3 % WPI (w/w) solution".
- **Vessel, fill, seal (sec. 2.3/2.4):** "each 1 mL of sample solution in glass vials (item no. ML
  33003 V, Mikrolab Aarhus A/S, Aarhus, Denmark) were sealed with silicone/polytetrafluoroethylene
  (PTFE) septa and crimp caps (item no. ML 33034L ...)". GC-FPD: "One milliliter of each untreated
  protein or amino acid sample in phosphate buffer was transferred to a sealed glass vial as
  described in the '2.3 Heat Treatments for Proteomic Analysis' section above."
- **Heating:** "Samples were then heated following the abovementioned temperature/time combinations
  in the preheated aluminum block" (VWR cat. 460-0008). Combinations for GC: 80 °C / 10 min (not for
  Cys), 90 °C / 10 min, 90 °C / 120 min, and "160 °C for 160 s ... to mimic an ultra-high-temperature
  (UHT) treatment".
- **Quench:** stated only for the UHT-like treatment: "The samples were first cooled at room
  temperature for 1 min and then immersed in ice-water."
- **Storage before analysis:** "stored at 10 °C before analysis". Duration **not stated**.
- **Headspace gas:** **not stated** (no flushing described → presumably air).
- **Vial volume / headspace volume:** **not stated**; only the part number ML 33003 V.

---

## 3. ANALYTICAL METHOD [M]

- **Sampling:** "After five min incubation at 32 °C in the autosampler (G6501B GC Sampler 80 ...),
  250 μL of headspace was injected through a splitless injection liner at 250 °C". So the headspace
  is sampled **after cooling and storage, at 32 °C**, not at the heating temperature.
- **GC:** Agilent GC7890B; DB-WAX Ultra Inert 30 m × 0.25 mm × 0.25 µm; H2 1.2 mL/min; 32 °C 6 min,
  50 °C/min to 230 °C (5 min), 5 min hold. FPD Plus in sulfur mode at 255 °C.
- **Calibration:** "A H2S stock solution was prepared by dissolving 385 mg of sodium sulfide in 50 g
  of phosphate buffer (pH 8.5) ... The stock solution was diluted into a series of different
  concentrations in phosphate buffer at pH 7.0 for the preparation of the standard curve."
  **Not stated:** number of levels, range, curve form (FPD sulfur response is classically
  non-linear), whether standards were 1 mL in the same vials with the same 32 °C equilibration,
  and the Na2S hydrate. **[Z]** 385 mg / 50 g = 98.7 mM if anhydrous Na2S, 32.1 mM if Na2S·9H2O.
- **LOD:** "LOD = 3.3 × (standard deviation of the standard curve/slope of the standard curve)";
  **LOD = 1.3 µM** (Table 1 footnote; text "LOD = 1.3 µM"). **LOQ, recovery: not stated.**
- **Replicates:** "All heating (or control) experiments were performed as independent triplicates,
  unless otherwise mentioned, and results were shown as mean values ± standard deviation (SD)."
  Table 1 footnote overrides: "tests were conducted with at least two independent replicates, while
  only one measurement was performed for the values without standard deviation."
- **Units:** µM, i.e. per volume of sample liquid (the paper does not say "liquid-equivalent", but
  the calibration is against Na2S solutions in pH 7.0 buffer).

---

## 4. TABLE 1 — GC-FPD H2S (µM) in heat-treated samples [M], verbatim

| sample | unheated | 80 °C (10 min) | 90 °C (10 min) | 90 °C (120 min) | UHT-like (160 °C, 160 s) |
|---|---|---|---|---|---|
| α-LA | <LOD | NT | <LOD | <LOD | <LOD |
| β-LG | <LOD | 2.0 ± 0.1 *; a | 7.3 ± 0.2 **; c | 15.4 ± 0.7 ***; c | 22.2 ± 1.8 ****; a |
| Whey model | <LOD | 2.9 ± 0.5 *; b | 9.5 ± 0.6 **; d | 18.7 ± 1.0 ***; d | 30.1 ± 3.0 ****; a,b |
| WPI | <LOD | <LOD | 4.0 ± 0.1 *; b | 23.3 ± 0.4 **; e | 38.1 ± 8.6 ***; a,b |
| Concentrated α-LA | <LOD | NT | 1.5 ± 0.2 *; a | 4.0 ± 0.3 *; a | 51.6 ± 17.9 **; c |
| **Cys** | <LOD | NT | **3.8 ± 0.3 b** | **10.7** (n = 1) | **44.8** (n = 1) |
| Met | <LOD | NT | <LOD | <LOD | <LOD |
| Mixed amino acid (Cys + Met) | <LOD | NT | 4.3 ± 0.6 *; b | 11.6 ± 0.8 **; b | NT |

Footnote: α-LA 2.45 mg/mL (0.17 mM); β-LG 15.88 mg/mL (0.87 mM); whey model both; WPI 30 mg/mL;
conc. α-LA 16 mg/mL (1.11 mM); Cys 5.7 mM; Met 3.6 mM; mixed 5.7 mM Cys + 3.6 mM Met. Asterisks:
row-wise differences between treatments; letters: column-wise between samples (one-way ANOVA,
Tukey-Kramer, p < 0.05). NT not tested. LOD 1.3 µM.

**[!]** The Cys 120-min and UHT values carry no SD → single measurements. The mixed-amino-acid arm
(same 5.7 mM Cys plus Met, which itself gives no H2S) is the replicated near-twin: 4.3 ± 0.6 and
11.6 ± 0.8 µM, consistent with the Cys arm.

---

## 5. TABLES 3, 4, 5 — protein thiols, lanthionine, derived extents [M]

### Table 3. Free thiol (mM, DTNB, Cys standards 0.1-1.0 mM in 55 mM HEPES pH 7.0; n = 3; LOD 0.003 mM)

| sample | unheated | 70 °C 10 min | 80 °C 10 min | 90 °C 10 min | 90 °C 120 min | UHT-like |
|---|---|---|---|---|---|---|
| α-LA | <LOD | NT | <LOD | 0.00 ± 0.01 * | 0.01 ± 0.01 * | 0.02 ± 0.01 ** |
| Conc. α-LA | 0.01 ± 0.01 * | NT | 0.01 ± 0.01 * | 0.01 ± 0.01 * | 0.05 ± 0.01 ** | 0.12 ± 0.01 *** |
| β-LG | 0.78 ± 0.06 ** | 0.77 ± 0.06 ** | 0.75 ± 0.06 ** | 0.73 ± 0.06 **; a,b | 0.48 ± 0.06 *; a | 1.20 ± 0.14 ***; b |
| Whey model | 0.77 ± 0.05 ** | 0.77 ± 0.05 ** | 0.74 ± 0.05 ** | 0.70 ± 0.05 **; a | 0.45 ± 0.05 *; a | 1.07 ± 0.12 ***; a |
| WPI | 0.77 ± 0.06 ** | 0.76 ± 0.06 ** | 0.77 ± 0.06 ** | 0.78 ± 0.06 **; b | 0.60 ± 0.06 *; b | 1.09 ± 0.07 ***; a,b |

No free-cysteine row: **residual cysteine in the Cys arm was not measured.** Thiol samples were the
200 µL HEPES / PCR-tube preparations (sec. 2.3), not the GC vials.

### Table 4. Lanthionine (µM) after acid hydrolysis (LC-MS, lysine-d4 IS, 10-4,000 ng/mL; LOD 0.4 µM)

| sample | unheated | 80 °C 10 min | 90 °C 10 min | 90 °C 120 min | UHT-like |
|---|---|---|---|---|---|
| α-LA | <LOD | <LOD | <LOD | 1.0 | 1.5 |
| β-LG | <LOD | <LOD | <LOD | 3.9 ± 0.5 **; a | 4.2 ± 0.9 **; a |
| Whey model | <LOD | 1.6 ± 0.6 *; a | 2.0 ± 0.4 *; a | 11.8 ± 1.7 **; b | 14.4 ± 2.1 **; b |
| WPI | <LOD | 1.3 ± 0.1 *; a | 2.6 ± 0.6 *; a | 11.5 ± 1.3 **; b | 21.4 ± 5.0 ***; b |

α-LA values "obtained from the mixture of five independent heat treatments".

### Table 5. (a) disulfide cleavage, Eq. 3; (b) Cys converted to H2S, Eqs. 1-2 (%)

| (a) cleavage | 80 °C 10 min | 90 °C 10 min | 90 °C 120 min | UHT-like |
|---|---|---|---|---|
| α-LA | N/A | 0.10 % | 0.42 % | 1.44 % |
| Conc. α-LA | 0.03 % | 0.09 % | 0.45 % | 1.34 % |
| β-LG | N/A | N/A | N/A | 12.15 % |
| Whey model | N/A | N/A | N/A | 6.13 % |
| WPI | N/A | N/A | N/A | 6.57 % |

| (b) Cys → H2S | 80 °C 10 min | 90 °C 10 min | 90 °C 120 min | UHT-like |
|---|---|---|---|---|
| β-LG | 0.25 % | 0.93 % | 1.93 % | 1.81 % |
| Whey model | 0.37 % | 1.21 % | 2.37 % | 2.74 % |
| WPI | N/A | 0.52 % | 2.94 % | 3.38 % |
| Conc. α-LA | N/A | 10.74 % | 8.16 % | 29.31 % |

Eq. 1: extent = [H2S]heated / ([free thiol]unheated + [H2S]heated) × 100 %; Eq. 2 uses
[free thiol]heated when it exceeds the unheated value. Text: H2S "only accounted for a small (<5%)
fraction of all Cys reactions".

### Table 2 (not transcribed in full)

Peptide-bound DHA and persulfide signals, **relative** (each residue's normalised peak area divided by
its highest value across treatments; n = 3). Not absolute; not transcribed beyond this: β-LG (Glu-C)
Cys66 DHA ND / 10.0 ± 1.8 / 14.9 ± 3.2 / 46.3 ± 7.9 / 100.0 ± 23.0 % (unheated / 80 / 90 10 min /
90 120 min / UHT); Cys160 DHA 3.4 ± 0.9 / 3.9 ± 1.0 / 6.1 ± 1.6 / 25.2 ± 6.4 / 100.0 ± 35.7 %.

---

## 6. THE FREE-CYSTEINE ARM AGAINST THE ENGINE [Z]

Engine constant: `k_cys_h2s` = `cysteine_thermolysis_k(pH, T)` in
`src/kinetic_core/parameters_sulfur.py` (Zheng & Ho 1994; pH 7 pair A = 2.36e13 s⁻¹ = 1.416e15 min⁻¹,
Ea = 135.5 kJ/mol; Zheng & Ho Table I measured k(90 °C, pH 7) = 5.99e-5 min⁻¹; zheng1994 dossier
sec. 3, 5b). Under B51 the table reads the 25 °C label pH, here 7.0.

No-sink first-order release: [H2S](t) = C0 (1 − e^(−kt)), C0 = 5.7 mM = 5700 µM.

| condition | k (min⁻¹) | k·t | predicted µM | measured µM | measured / predicted |
|---|---|---|---|---|---|
| 90 °C, 10 min, Table I k | 5.99e-5 | 5.99e-4 | 5700 × 5.99e-4 = **3.41** | 3.8 ± 0.3 | 1.11 |
| 90 °C, 10 min, engine Arrhenius | 4.586e-5 | 4.586e-4 | **2.61** | 3.8 ± 0.3 | 1.46 |
| 90 °C, 120 min, Table I k | 5.99e-5 | 7.19e-3 | 5700 × 0.00716 = **40.8** | 10.7 (n = 1) | **0.26** |
| 90 °C, 120 min, engine Arrhenius | 4.586e-5 | 5.50e-3 | 5700 × 0.00549 = **31.3** | 10.7 (n = 1) | **0.34** |
| 160 °C, 160 s (2.667 min), engine Arrhenius | 6.473e-2 | 0.1726 | 5700 × 0.1585 = **904** | 44.8 (n = 1) | 0.050 |

(engine k computed in the `maillard` env with `cysteine_thermolysis_k(7.0, T)`; the 160 °C value is
an extrapolation 50 °C beyond Zheng & Ho's 80-110 °C and assumes the liquid is at 160 °C for the whole
160 s, which a 1 mL vial in a block cannot be — it is an upper bound, not a prediction.)

**Shape.** Measured 120-min / 10-min ratio = 10.7 / 3.8 = **2.8**; first-order release with no sink
gives 12.0. The increment 10 → 120 min is 6.9 µM against a predicted 28.7-37.4 µM.

**The engine's own sink.** The shipped lumped `r_h2s_loss` (k = 0.0446 min⁻¹ at 145 °C, Ea
64.1 kJ/mol → 2.74e-3 min⁻¹ at 90 °C) brings the 120-min prediction to 26.7 µM (still 2.5× high).
A first-order H2S sink that reproduces 10.7 µM at 120 min needs k ≈ **0.023 min⁻¹** (engine source)
or **0.031 min⁻¹** (Table I source), t½ 22-30 min at 90 °C, ~8-11× the shipped lump; it then
predicts 2.3-2.9 µM at 10 min (measured 3.8 ± 0.3), so a pure sink fits the end point but undershoots
the start. Precursor loss (cysteine → cystine by headspace O2) is the other reading; the paper cannot
separate the two.

**Why partitioning and storage do not explain it [Z].** (i) The vial stays sealed from heating to
injection and is sampled at 32 °C; whatever went to the gas at 90 °C returns on cooling, and the
calibration is an external Na2S curve in the same buffer, so partitioning is (if standards were
handled identically, not stated) calibrated out. It is also time-independent and cannot bend the
10 → 120 min ratio. (ii) Storage at 10 °C is the same for every vial; a first-order storage loss
scales both points by the same factor. It can bend the ratio only if the sink acting during storage
was itself made during heating (e.g. cystine or polysulfide accumulated in the 120-min vial).
(iii) Septum (silicone/PTFE) sorption or permeation of H2S at 90 °C for 2 h is not addressed and is
time-dependent.

---

## 7. WHAT THE PAPER DOES NOT STATE THAT A KINETIC MODEL NEEDS

1. **Vial total volume and headspace volume** (only "ML 33003 V"); hence no headspace O2 budget.
2. **Headspace gas** (air assumed; no N2 flush described). With any air headspace of ≥ 1 mL, O2
   (≥ ~8.7 µmol per mL of air at 20 °C) exceeds what is needed to oxidise all 5.7 µmol cysteine to
   cystine (1.4 µmol O2) **[Z]**.
3. **Thermal history of 1 mL in a block**: no heat-up or cool-down profile; quench stated only for
   UHT-like.
4. **Storage time at 10 °C before injection.**
5. **pH after heating** (not measured); cysteine stock pH not re-adjusted (not stated).
6. **Residual cysteine / cystine** in the Cys arm (Table 3 has no Cys row).
7. **Calibration details**: levels, range, curve form, Na2S hydrate, standards' vial/volume/incubation.
8. **LOQ and recovery.** Single determinations for Cys at 120 min and UHT-like.
9. **Trace metals / chelator**: none added (10 mM phosphate only), so metal-catalysed thiol and
   sulfide autoxidation is uncontrolled.

---

## 8. FLAGS

1. **[!]** Cys 90 °C / 120 min = 10.7 µM and 160 °C / 160 s = 44.8 µM are n = 1.
2. **[!]** Headspace sampled cold (32 °C after 10 °C storage); the paper's µM is total sulfide per
   liquid volume only insofar as the standards were treated like the samples.
3. **[!]** UHT-like 160 s in a 1 mL sealed vial: the liquid will lag the block; any kinetic use needs
   a declared thermal history.
4. Protein arms are matrix (β-LG, WPI): H2S from protein-bound Cys via β-elimination per the
   authors; not engine species (`PROT_SS` is a disulfide site input only).
5. Registry: paper absent from `data/keys/papers.yml`.

## Use in the model

- **Check on `k_cys_h2s` at 90 °C, pH 7, 5.7 mM:** the 10-min point passes (1.1-1.5×). It is the
  only sealed, buffered, mM-scale point in the corpus, 18× more dilute than Zheng & Ho's 0.1 M.
- **Evidence for a time-dependent H2S (or precursor) loss in a sealed air-headspace vial at 90 °C**:
  0.26-0.34× of no-sink release at 120 min; needs a sink t½ ~22-30 min, ~10× the shipped lumped
  `k_h2s_loss` at 90 °C. Single determination; the replicated mixed-amino-acid twin (11.6 ± 0.8 µM)
  corroborates the level.
- **Not usable for headspace partitioning**: the measurement is cold and externally calibrated.
- **160 °C / 160 s**: not simulable without a thermal history.
</content>
</invoke>
