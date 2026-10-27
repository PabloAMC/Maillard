# Zhai et al. 2023 (JAFC 71:16260) — EXTRACTION (purified TTCA vs purified Cys-ARP, 50 mmol/L each, pH 7.4, 120 °C, 0-140 min: total volatiles and A420 time courses; Cys + Xyl 82.7 mmol/L at 90 °C under vacuum dehydration 0-30 min and spray-drying at 170-190 °C inlet: 3-DX, MGO, GO, water states — a PROCESS / WATER-STATE paper)

### No TTCA or ARP concentration is plotted against time anywhere in this paper. The one kinetic comparison it carries — TTCA vs its Amadori isomer heated side by side at 120 °C — survives only as total-volatile and A420 curves, and the headline claim that "TTCA and ARP were completely consumed at 130 and 100 min" is attributed to a reference (Zhai 2019) that contains no such experiment.

**Source on disk:** `data/articles/zhai2023c.pdf` (10 pp., typeset JAFC; download stamp "DOKUZ EYLUL
UNIV user on 30 August 2026"). Prose from `pdftotext` (raw and `-layout`); **Figs. 1a, 1c, 2b, 5, 6a,
6b digitised from 110-400 dpi renders** (scratchpad `img/z23d_f1a-04.png`, `z23d_f1c-04.png`,
`z23d_f5-06.png`, `z23d_p7-07.png`, `zhai23d-05.png`), tagged **[D]**. **There is no table in the main
text.** Supporting Information (Fig. S1, **Table S1**, **Table S2** calibration curves, **Fig. S2**
TTCA/ARP yields during vacuum dehydration, **Fig. S3** TTCA/ARP distribution after spray drying) is
**not on disk** — and Figs. S2/S3 are exactly the TTCA/ARP data the brief is after.

## 0. Identity

| field | value |
|---|---|
| Title | "Intrinsic Molecular Mechanisms of Transformation between Isomeric Intermediates Formed at Different Stages of Cysteine−Xylose Maillard Reaction Model through Dehydration" |
| Authors | Yun Zhai, Khizar Hayat, Teng Li, **Yuying Fu** (corresponding, webfu@126.com), **Chi-Tang Ho** (corresponding, ctho@sebs.rutgers.edu) |
| Affiliations | Zhejiang Gongshang University (Hangzhou); Miami University (Ohio); Rutgers University |
| Venue | J. Agric. Food Chem. **2023, 71 (43), 16260-16269**. Received 11 Sep 2023, revised 2 Oct 2023, accepted 9 Oct 2023, published 19 Oct 2023 |
| DOI as printed | **10.1021/acs.jafc.3c06506** ✔ matches the brief |
| Funding | NSFC 32001821, 32172330; National Key R&D 2017YFD0400105; Jiangsu KYCX20_1880; King Saud RSP-2023R437 |
| Group references used | ref 7 = Zhai 2019 (`zhai2019_extraction.md`); ref 13 = Zhai 2020 (`zhai2020_extraction.md`, spray-drying TTCA/ARP ratio); ref 19 = Zhai 2021 (`zhai2021_extraction.md`); ref 10 = Cui et al. (Xyl-Phe thermal reaction-vacuum dehydration) |
| File-stem note | This is the **fourth** Zhai 2023 item in the corpus: `zhai2023b` / `zhai2023jafc` (JAFC 71:14300, TTCA + Cys), `zhai2023foodchem`, `zhai2023c`. Different paper, different DOI. |

## 1. Why it matters (and the limits)

**(a) TTCA vs Cys-ARP reactivity at 120 °C, same loading, same pH.** Fig. 1a/1c heat purified TTCA
and purified Cys-ARP side by side at **50 mmol/L, pH 7.4, 120 °C**. The ARP pot's total volatiles peak
at **60 min (40.63)** and the TTCA pot's at **100 min (38.87)** [printed]; A420 rises faster from ARP
at every time [D]. The authors explain the TTCA lag by "the equilibrium conversion from TTCA to ARP
dominat[ing] the reaction process during the initial stage" — i.e. **TTCA reaches its products
through the Amadori isomer**. The engine's TTCA has two one-way fates (`r_ttca_cys`, `r_ttca_deg`
→ DPO) and **no Cys-Amadori**; this is the third paper (with `cao2017_extraction.md` and
`zhai2019_extraction.md`) to put TTCA ⇌ Cys-ARP on the main path.

**(b) The lag is a shape the engine can be tested against.** TTCA's volatile curve is sigmoidal —
2.4 / 5.0 / 11.5 at 20 / 40 / 60 min, then 29.6 at 80 min [D] — a ~40 min induction that a one-step
TTCA → DPO → … chain produces only weakly. A two-step TTCA ⇌ Cys-ARP → DX chain gives it naturally.

**(c) Dehydration shifts the condensation equilibrium toward TTCA.** The vacuum-dehydration and
spray-drying arms say TTCA forms by an equilibrium dehydration ("the decrease of moisture content
could make the equilibrium migrate from the dehydration of N-substituted D-xylosamine to TTCA") —
consistent with Cys + Xyl ⇌ TTCA + H2O being reversible. **This is the water-activity axis the sulfur
lane does not model**; nothing here transfers to an aqueous pot.

What this paper does NOT give: any TTCA or ARP concentration against time in the main text, any rate
constant, any activation energy, any free cysteine or xylose measurement, any individual volatile,
any 1-DX value, or its SI.

### How it relates to the companion dossiers

| | `zhai2019_extraction.md` | `zhai2020_extraction.md` | **this paper** | `zhai2023b_extraction.md` |
|---|---|---|---|---|
| question | TTCA formation vs T; Cys effect | spray-drying shifts TTCA/ARP ratio | TTCA vs ARP reactivity; water states during dehydration | added Cys vs TTCA |
| TTCA vs t | **yes**, 70-130 °C | no | **no** (claimed, not shown) | no |
| ARP fed | no | no | **yes, 50 mM** | no |
| a_w / dehydration | no | spray drying | **vacuum + spray drying, LF-NMR** | no |
| ARP share after drying | — | 20.83 % (190 °C inlet), 47.23 % (pH 9.5 stock) | **23.75 %** (190 °C inlet) | — |

## 2. Methods as they matter to a model

- **TTCA/ARP synthesis.** "An **equimolar** solution of Cys and Xyl (**0.0827 mol/L**) ... maintained at
  pH 7.4 ± 0.01 using NaOH (6 mol/L), followed by heat treatment for **40 min at 90 °C**"; Dowex 50WX4
  (H⁺), semi-prep RP-HPLC XBridge Amide; identity by UPLC-ESI-MS and NMR (Fig. S1, Table S1, not on
  disk). **No purity figure printed.** (Zhai 2019, cited as the method, uses Xyl : Cys = 2 : 1 — see
  Flag 4 of `zhai2019_extraction.md`.)
- **TTCA/ARP thermal models.** "Model reaction systems of **TTCA (50 mmol/L) and ARP (50 mmol/L)** were
  prepared and adjusted to a pH value of **7.4 ± 0.01** with NaOH solution. Then ... treated under
  **120 °C in the pressure-resistant bottles** for different times (**20, 40, 60, 80, 100, 120, and
  140 min**)"; ice bath. Volume and fill not stated; unbuffered.
- **Quantitation (HPLC-ELSD, external standards).** Calibrations printed (x, y as in the paper; ELSD
  in log-log form):

  | analyte | calibration | R² |
  |---|---|---|
  | TTCA | lg y = 0.942 × lg x + 4.078 | 0.998 |
  | ARP | lg y = 0.801 × lg x + 4.241 | 0.991 |
  | GO | y = 2.221 × 10⁶ x + 4.631 × 10⁴ | 0.999 |
  | MGO | y = 1.618 × 10⁶ x + 3.470 × 10⁴ | 0.999 |
  | DA (diacetyl) | y = 3.173 × 10⁶ x + 1.931 × 10⁴ | 0.998 |
  | 3-DX | y = 2.730 × 10⁶ x + 3.292 × 10⁴ | 0.994 |
  | 1-DX | y = 2.785 × 10⁶ x + 6.828 × 10⁴ | 0.994 |

  Units of x and y are **not stated** for these curves; ranges not stated.
- **Thermal reaction-vacuum dehydration.** The Cys-Xyl solution, 90 °C water bath 30 min at
  atmospheric pressure, then **10 mbar within 2 min** (Büchi R-215), dehydration **at 90 °C** for
  0, 5, 10, 15, 20, 25, 30 min, ice bath.
- **Spray drying.** Step 1: 90 °C, 30 min; step 2: vacuum dehydration **at 40 °C**, 10 mbar, 40 min →
  concentrate **193.33 mg/mL**, pH re-set to **7.4 ± 0.1**; Labplant SD-06 at inlet / flow / outlet
  **170 °C-450 mL/h-65 °C, 180 °C-500 mL/h-70 °C, 190 °C-550 mL/h-75 °C**.
- **Water.** Moisture content by mass loss; a_w by Novasina; LF-NMR (0.5 T, 23 MHz nominal, CPMG:
  SF1 21 MHz, offset 401.952 kHz, P1 8.00 µs, P2 15.04 µs, TD 48 036, TW 1000 ms, NS 6, NECH 800);
  T₂₁ 1-10 ms bound, T₂₂ 10-100 ms immobilised, T₂₃ 100-1000 ms free.
- **A420** after dilution "to a certain multiple" (factor not stated).
- **Volatiles.** 3 mL sample + 2 mL saturated NaCl + 3 µL 1,2-dichlorobenzene (**25 mg/mL**) in a
  15 mL vial; DVB/CAR/PDMS 75 µm, 60 °C 20 min; Agilent 7890B/5977B, DB-Wax; quantified by
  calibration curves (Table S2, not on disk) in µg/L. Oven ramps printed **in °C/s** (Flag 6).
- **Statistics.** Triplicate experiments, each analysed in triplicate; mean ± SD; p < 0.05.

## 3. Values re-typed (no main-text table)

### 3.1 Values printed in the running text

| quantity | value | conditions | where |
|---|---|---|---|
| time to complete consumption, **TTCA** | **130 min** | (attributed to ref 7) 120 °C, pH 7.4 implied | Results, first section |
| time to complete consumption, **ARP** | **100 min** | same | same |
| total volatile maximum, TTCA | **38.87 "µg/mL" at 100 min** | 50 mM, 120 °C, pH 7.4 | same (Fig. 1a axis says µg/L — Flag 2) |
| total volatile maximum, ARP | **40.63 "µg/mL" at 60 min** | same | same |
| PCA of volatiles | PC1 98.845 %, PC2 1.111 % | TTCA vs ARP | Fig. 1b |
| bound water, vacuum dehydration | **57.21 → 75.39 %** (20 → 30 min) | 90 °C | Fig. 3b text |
| free water | **7.26 → 5.21 %** (20 → 30 min) | 90 °C | same |
| immobilised water | "decreased rapidly in the early stage and reached a balance after 26 min" | | same |
| A420 during vacuum dehydration | **0.018 → 0.179** (0 → 20 min) | 90 °C | Fig. 5 text |
| 3-DX, vacuum dehydration | **199.21 → 208.49 µmol/L** (0 → 10 min) | 90 °C | Fig. 5 text |
| **TTCA yield at 20 min of dehydration** | **79.78 %** | 90 °C, 10 mbar | Results |
| ARP / 3-DX / MGO / GO as % of TTCA, 20 min | **3.03 / 0.45 / 0.21 / 0.18 %** | same | Results |
| TTCA yield at 30 min | **73.06 %** | same | Results |
| ARP / 3-DX / MGO / GO at 30 min relative to 20 min | **1.05 / 1.32 / 1.07 / 1.12 ×** | same | Results |
| ARP share of MRI after spray drying, 190 °C inlet | **23.75 %** | stock pH 7.4 | Results (Fig. S3, not on disk) |
| 3-DX after spray drying | **0.199 → 0.208 / 0.219 / 0.249 mmol/L** (stock → 170 / 180 / 190 °C) | | Results / Fig. 6a |
| ARP yield 190 vs 180 °C inlet | **1.90 ×**; 3-DX **1.14 ×** | | Results |
| a_w of spray-dried products | **"10.172"**, 0.167, 0.139 (170 / 180 / 190 °C) | | Results (Flag 3) |
| bound water after spray drying | **62.29 → 49.07 %** (170 → 190 °C) | | Results |

### 3.2 Digitised Fig. 1a — total volatile flavour compounds, 120 °C, 50 mmol/L, pH 7.4 [D]

Axis label "Total contents of flavor compounds (µg/L)", 0-50; read error **± 0.5**.

| time (min) | 0 | 20 | 40 | 60 | 80 | 100 | 120 | 140 |
|---|---|---|---|---|---|---|---|---|
| TTCA | 0 | 2.4 | 5.0 | 11.5 | 29.6 | **39.0** | 34.4 | 30.2 |
| ARP | 0 | 7.9 | 20.3 | **40.8** | 38.5 | 31.2 | 28.3 | 21.3 |

### 3.3 Digitised Fig. 1c — A420, same pots [D], read error ± 0.01

| time (min) | 0 | 20 | 40 | 60 | 80 | 100 | 120 | 140 |
|---|---|---|---|---|---|---|---|---|
| TTCA | 0.081 | 0.146 | 0.279 | 0.318 | 0.395 | 0.439 | 0.502 | 0.599 |
| ARP | 0.081 | 0.263 | 0.423 | 0.530 | 0.640 | 0.680 | 0.733 | — (turbid) |

### 3.4 Digitised Fig. 5 — vacuum dehydration at 90 °C [D]

Left axis A420 (± 0.003), right axis µmol/L (± 4).

| dehydration time (min) | 0 | 5 | 10 | 15 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|---|
| A420 | 0.018 | 0.023 | 0.045 | 0.072 | 0.179 | 0.198 | 0.234 |
| 3-DX | 199 | 203 | 208 | 217 | 296 | 355 | 389 |
| MGO | 121 | 117 | 119 | 119 | 135 | 140 | 145 |
| GO | 113 | 118 | 118 | 117 | 120 | 132 | 136 |

### 3.5 Digitised Figs. 2b, 6a, 6b [D] (coarse)

Fig. 2b a_w during vacuum dehydration (20-30 min): 0.208, 0.204, 0.192, 0.190, 0.188, 0.187 (± 0.002).
Fig. 6a (µmol/L, ± 3): control (stock) 3-DX 199 / MGO 122 / GO 112; 170 °C 208 / 122 / 118; 180 °C
219 / 128 / 121; 190 °C 249 / 131 / 128. Fig. 6b: a_w ≈ 0.172 / 0.167 / 0.139; combined / immobilised /
free water ≈ 62 / 35 / 3, 54 / 43 / 2, 49 / 50 / 1 % at 170 / 180 / 190 °C.

### Arithmetic (mine)

**1. Cross-check of the printed percentages against Fig. 5.** If 3-DX at 20 min (296 µmol/L [D]) is
0.45 % of TTCA, TTCA = 66 mmol/L; and 79.78 % × 82.7 mmol/L Cys = **66.0 mmol/L** ✔. MGO 0.21 % →
139 µmol/L vs 135 [D] ✔; GO 0.18 % → 119 vs 120 [D] ✔; 3-DX 30/20 min = 389/296 = 1.31 vs printed
1.32 ✔. **The figure and the text agree, on a basis of the original (pre-evaporation) volume.**
ARP at 20 min = 3.03 % of 66 = **2.0 mmol/L**; ARP / TTCA ≈ 0.03 — the same small isomer share as
Zhai 2019 at pH 7.4 (≤ 0.06), against Cao 2017's 0.56 at pH ~4.2.

**2. What "complete consumption by 130 min" would imply.** For TTCA 50 mM at 120 °C, pH 7.4: ≥ 99 %
loss in 130 min needs a first-order k ≥ ln(100)/130 = **≥ 0.035 min⁻¹**, or a zero-order k ≥
0.38 mmol L⁻¹ min⁻¹. Zhai 2021's purified TTCA 10 mM at 120 °C, pH 7 fits a zero-order **0.0651
mmol L⁻¹ min⁻¹** (first-order equivalent 0.0065 min⁻¹), leaving 6.07 mM at 60 min. The two are
**5-6× apart in first-order terms** at the same temperature and pH. Either the claim is wrong, or
TTCA loss is strongly concentration-dependent (higher order, or autocatalytic through pH drift), or
the "complete consumption" refers to a different experiment. **Treat as an unverified bound.**

**3. Volatile yield, if the axis unit (µg/L) is right.** 39 µg/L total from 50 mmol/L (12.7 g/L)
TTCA is ~3 × 10⁻⁶ by mass — the same trace-branch order as Zhai 2023b's 87 µg/L from 10 mM. If the
text unit (µg/mL) were right it would be 1000× larger and inconsistent with every other Zhai paper.

## 4. Kinetic numbers the repository can use

**No rate constant, no activation energy, no TTCA/ARP time series in the main text.**

| step / observable | quantity | value | unit | conditions | order | source | evidence class |
|---|---|---|---|---|---|---|---|
| TTCA vs Cys-ARP reactivity | time of volatile maximum | **100 vs 60** | min | 50 mM, pH 7.4, 120 °C | — | text + Fig. 1a | **within_study_ratio** (ARP ~1.7× faster to peak) |
| total volatiles vs t | Fig. 1a table | µg/L (axis) | 50 mM, 120 °C | — | Fig. 1a [D] | level_only (unit conflict, Flag 2; calibration SI not on disk) |
| browning vs t | A420, Fig. 1c table | — | same | — | Fig. 1c [D] | level_only (dilution not stated) |
| TTCA consumed | "completely" by **130 min** | — | 50 mM, 120 °C, pH 7.4 | — | text, cites ref 7 | **threshold, unverified** (Flag 1) |
| Cys-ARP consumed | "completely" by **100 min** | — | same | — | same | threshold, unverified |
| TTCA lag | volatiles ≤ 11.5 until 60 min, then rise | — | same | — | Fig. 1a [D] | shape (sigmoid) |
| TTCA yield under vacuum dehydration | 79.78 % (20 min), 73.06 % (30 min) | % of Cys | 90 °C, 10 mbar | — | text | level_only — **dehydrated system, outside the aqueous lane** |
| ARP / TTCA under dehydration | 0.030 → 0.032 | — | same | — | text, arithmetic 1 | within_study_ratio |
| 3-DX, MGO, GO vs dehydration time | Fig. 5 table | µmol/L | same | — | Fig. 5 [D] | level_only (concentration basis = original volume) |

## 5. What this can and cannot anchor in the model

**Can:**
1. **A structural argument for a Cys-Amadori node between TTCA and the deoxyosones.** Same loading,
   same pH, same temperature: ARP reaches its volatile peak in 60 min, TTCA in 100 min, with a
   ~40 min induction. Together with Cao 2017 (TTCA → Cys-Amadori with conserved cysteine pool) and
   Zhai 2019 (3-DX > 1-DX from fed TTCA), it argues that `r_ttca_deg` (TTCA → DPO in one step) is a
   lumped surrogate for TTCA ⇌ CysARP → DPO + TDP.
2. **A shape-only hold-out** on fed TTCA at 120 °C / pH 7.4 (sigmoid onset of downstream products),
   once an engine run can report a total-volatile or A420 proxy — weak.

**Cannot:**
- Any rate: there is no TTCA or ARP concentration against time in the main text; the "130 / 100 min"
  claim is unsupported by its citation and inconsistent with Zhai 2021 (arithmetic 2).
- Anything from the dehydration / spray-drying arms: they are low-a_w systems (a_w 0.14-0.21), which
  the sulfur lane explicitly does not model.
- Absolute volatile levels (unit conflict; SI calibrations not on disk).

## 6. Flags

1. **The decay claim is cited to a paper that does not contain it.** "According to our previous
   experiments,⁷ the concentration of TTCA/ARP showed a gradual downward trend during the thermal
   reaction, and TTCA and ARP were completely consumed at 130 and 100 min". Ref 7 is Zhai 2019
   (JAFC 67:8632), which never heats purified ARP and never follows TTCA at 50 mmol/L. The nearest
   real data (Zhai 2021, 10 mM, 120 °C) show 61 % TTCA left at 60 min. Do not use these two times as
   measurements.
2. **Unit conflict, ×1000.** Text: "100 (38.87 µg/mL) and 60 min (40.63 µg/mL)"; Fig. 1a axis: "Total
   contents of flavor compounds (µg/L)". The figure unit is consistent with the rest of the series
   (Zhai 2023b: 13.7-87.3 µg/L); the text unit is very probably a typo.
3. **a_w = "10.172"** for the 170 °C spray-dried product is impossible (a_w ≤ 1); Fig. 6b reads
   ≈ 0.172.
4. **The 3-DX baseline is the same number in two different experiments.** Vacuum dehydration starts
   at **199.21 µmol/L** (Fig. 5, after 30 min at 90 °C, before evacuation); the spray-drying stock
   ("control") is **0.199 mmol/L** after 90 °C 30 min + 40 °C vacuum 40 min concentration to
   193.33 mg/mL. Identical to three figures despite a ~2-4× concentration step between them;
   either the spray-drying "control" is the same pre-vacuum solution, or one value was copied.
5. **The dehydration time axis is ambiguous.** Text: vacuum dehydration sampled at 0-30 min; Figs. 2b
   and 3b start at 20 min ("the MRI solution with a moisture content of 0% was treated under further
   vacuum dehydration"). The 20 min point is where free water is gone.
6. **GC oven ramps printed in °C/s** ("raised to 80 °C at a rate of 3 °C/s ... 4 °C/s ... 10 °C/s"),
   physically implausible for a capillary GC programme; read °C/min. The UGC-E-nose programme is
   legitimately in °C/s.
7. **"Equimolar 0.0827 mol/L" vs Zhai 2019's 2 : 1 Xyl : Cys** for the cited preparation (see
   `zhai2019_extraction.md` Flag 3).
8. **Calibration units missing:** none of the seven calibration lines states x or y units or range;
   the TTCA/ARP lines are log-log (ELSD), so an unknown range matters for extrapolation.
9. **Cys-ARP ≠ engine `ARP`** (alanine Amadori). Never map.
10. **Not on disk:** Figs. S1-S3, Tables S1-S2 — including **Fig. S2 (TTCA and ARP yields against
    vacuum-dehydration time)** and **Fig. S3 (TTCA/ARP distribution after spray drying)**.
11. **Registry gaps:** TTCA, Cys-ARP, 1-DX, 3-DX, cysteine, xylose absent; `glyoxal`,
    `methylglyoxal` and `2_3_butanedione` (diacetyl) present. Paper not in `papers.yml`.
12. **What to request:** TTCA and ARP concentrations vs time in the two 50 mM pots at 120 °C (the
    data behind the "130 / 100 min" claim); the SI; per-compound volatile data behind Fig. 1a; the
    A420 dilution factor; final pH of the 120 °C pots.
