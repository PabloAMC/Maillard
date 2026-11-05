# Chen, Hu, Dou, Meng & Deng 2026 (Curr. Res. Food Sci. 12:101433) — EXTRACTION (file stem `chen2026`)
### *"Exogenous xylose trapped threonine to regenerate Amadori compounds and accelerated pH decline impeding pyrazines formation upon thermal treatment of Amadori compounds revealed by sotope labeling"* (sic, title)
### Purified threonine-xylose Amadori compound (**Thr-ARP 20 mM**) alone or + D-xylose 10-40 mM, water, **initial pH 7.5 (NaOH), no buffer**, pressure-resistant glass bottles, oil bath **120 °C, 0-150 min**: Thr-ARP and released Thr (HPLC-ELSD), **3-DX, 1-DX, GO, MGO (OPD quinoxalines, HPLC-DAD)**, furfural (HPLC-UV), **pH measured on cooled samples**, five pyrazines at 150 min (HS-SPME-GC-MS); [¹³C5]-xylose shows labelled Thr-ARP re-forms. **Printed LODs: ELSD 0.20 µM (Thr-ARP) / 0.42 µM (Thr); OPD 0.47-0.59 µM; furfural 0.53 µM** — the ELSD values are implausibly low (sec. 5).

**Read 2026-11-05 (first reader) from `data/articles/chen2026.pdf`** (8 pp., Elsevier open access, clean text layer).
Prose and calibrations from `pdftotext -raw`; **Figs. 1, 2, 3, 5 digitised from 110 dpi page renders**, Fig. 1 re-read at
250 dpi (scratchpad `img/chen26-5.png`, `-6.png`, `-7.png`, `img/chen26_fig1-5.png`). Read errors: 3-DX ±0.003 mM,
1-DX ±0.002 mM, GO ±0.01 mM, MGO ±0.0007 mM, Thr / ARP ±0.2 mM, pH ±0.03. Supplementary (Fig. S1 furfural) **NOT on
disk**. Only the 0 / 20 / 40 mM xylose arms are plotted (10 and 30 mM appear in Table 1 only).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative statement by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me ·
**[!]** flag.

## 0. Identity

| field | value as printed |
|---|---|
| Authors | Pusen Chen (Chuzhou Polytechnic), Baishun Hu, Zuman Dou (Jimei Univ.), **Fei Meng*** (Chuzhou Polytechnic), **Qiong Deng*** (Enshi Academy of Agricultural Sciences) |
| Venue | *Current Research in Food Science* **12** (2026) 101433 |
| DOI | **10.1016/j.crfs.2026.101433** (printed) — matches the brief |
| Dates | received 4 Mar 2026; revised 27 Apr 2026; accepted 5 May 2026; online 6 May 2026 |

**Laboratory.** Corresponding authors at Chuzhou Polytechnic / Enshi Academy. The first author's methods are the Jiangnan
(Zhang X / Cui H / Ho) ones: Thr-ARP synthesis and purification "adapted from" Chen et al. 2022 (JAFC, Jiangnan + Ho
[C]), ELSD calibrations cited to Chen et al. 2024 (Jiangnan [C]), the OPD / HEPES / DTPA derivatisation, and the
Thr-ARP + Thr precursor paper Chen 2023 (JAFC, Jiangnan [C]). **Not in `data/keys/laboratories.yml`.** Proposal: a new
group `chuzhou_meng_deng` with the evidence "methods and first author from `jiangnan_zhang_x`"; for leave-one-lab-out
the conservative choice is to **merge with `jiangnan_zhang_x`** (P ≈ 0.6 that the analytical chain is Jiangnan's).
Owner call.

## 1. What this paper answers for the model

| model need | what Chen 2026 gives | tag |
|---|---|---|
| LOD floors (ELSD) | Thr-ARP **LOD 0.199 / LOQ 0.601 µM**; Thr **0.421 / 1.223 µM**; linear calibrations | [M] / [!] |
| LOD floors (OPD quinoxalines, DAD 315 nm) | **3-DX 0.59 / 2.03; 1-DX 0.58 / 2.29; GO 0.57 / 2.32; MGO 0.47 / 2.08 µM** (LOD / LOQ) | [M] |
| LOD floor (furfural, UV 277 nm) | **0.531 / 1.438 µM** | [M] |
| 1-DX / 3-DX vs time | Fig. 1, 0-150 min, Thr-ARP ± 20 / 40 mM Xyl; **3-DX peaks at 60 min, 1-DX at 30 min** | [F] |
| 1-DX / 3-DX vs pH | **No pH series** (one initial pH, 7.5); the pH drifts inside each pot (Fig. 5), so ratio-vs-measured-pH within pots only | [F] / [Z] |
| in-situ pH drift (hub H7) | **Measured (cooled to 25 °C ± 0.2)**: Thr-ARP alone 7.50 → 5.75 in 150 min; +20 mM Xyl → 4.88; +40 mM → 4.28 | [F] |
| ARP lifetime | Thr-ARP alone: apparent first-order **k ≈ 0.022 min⁻¹ over 0-90 min** at 120 °C, slowing after (pH falling) | [F] / [Z] |
| amine release | Thr up to 16.5 mM at 60 min — **exceeds the ARP consumed** [!] | [F] / [Z] |
| GO / MGO | GO 0.31-0.39 mM, MGO 0.021-0.028 mM (ARP alone, 30-150 min) | [F] |
| pyrazines | Table 1, 150 min, µg/L (IS-normalised external calibration) | [M] |

## 2. Conditions [M]

| item | as printed |
|---|---|
| Thr-ARP | xylose 0.5 M + Thr 0.25 M, pH 7.50, 90 °C 80 min, pH re-set 7.50, rotary evaporation 95 °C 12 min; Dowex 50WX4 (H⁺) then semi-preparative RP-HPLC; ¹H NMR and UPLC-ESI-MS. **Purity not stated** [!] |
| Model pots (sec. 2.3) | Thr-ARP **20 mM** with D-xylose **0 / 10 / 20 / 30 / 40 mM**; pH **7.5 ± 0.02** with "5 mmol/L sodium hydroxide" [V]; **no buffer**; pressure-resistant glass bottles (volume and fill not stated); oil bath **120 °C**; 0 / 30 / 60 / 90 / 120 / 150 min; ice-bath quench |
| Headspace / O2 / metals | not stated; DTPA only in the OPD reagent |
| ¹³C arm (sec. 2.8) | [¹³C5]-xylose 20 mM + Thr-ARP 20 mM, pH 7.5, 120 °C, 20 min; UPLC-ESI(−)-MS/MS |
| pH (sec. 2.9) | cooled in ice to ~25 °C (± 0.2), then measured: **label scale (cooled), not in situ** |
| Replicates | three independent replicates; mean ± SD (error bars mostly inside the markers) |

## 3. Assays and printed floors [M]

| analyte | method | calibration as printed | LOD (µM) | LOQ (µM) |
|---|---|---|---|---|
| Thr-ARP | HPLC-ELSD, XBridge Amide 4.6 × 150 mm, 0.1 % HCOOH water / MeCN, 8 µL | y = 3.429 × 10⁶ x − 2.548 × 10⁴ (R² 0.998) | 0.199 | 0.601 |
| Thr | same | y = 1.109 × 10⁶ x − 5.23 × 10⁵ (R² 0.998) | 0.421 | 1.223 |
| 3-DX | OPD (0.5 % OPD, 1 M HEPES, 12 mM DTPA; 0.4 + 0.4 mL; 30 °C 4 h dark); SunFire C18; DAD 315 nm; 5 µL | y = 4.79 × 10⁶ x − 2.19 × 10⁴ | 0.59 | 2.03 |
| 1-DX | same | y = 4.78 × 10⁶ x − 3.22 × 10⁴ | 0.58 | 2.29 |
| GO | same | y = 2.47 × 10³ x − 7.83 × 10⁴ | 0.57 | 2.32 |
| MGO | same | y = 3.46 × 10³ x − 1.44 × 10⁴ | 0.47 | 2.08 |
| furfural | HPLC-UV 277 nm, SunFire C18 | y = 9.01 × 10⁻⁸ x − 5.81 × 10⁻² (R² 0.996), x = area, y = mM | 0.531 | 1.438 |

The paper says, for the dicarbonyls, that "x denoted the peak area, and y ... concentration in mmol/L" [V, abridged].
[!] For 3-DX / 1-DX that reading is impossible (4.79 × 10⁶ mM per area unit); the slopes only make sense as area per mM.
GO / MGO slopes are **~1500-2000× smaller** than 3-DX / 1-DX on the same detector and derivatisation; quinoxaline molar
absorptivities at 315 nm do not differ by three orders [Z]. One of the two pairs has its units or axes mislabelled; the
GO / MGO concentrations in Fig. 2 therefore carry an **unresolved scale doubt**. 1-DX and 3-DX slopes agree to 0.2 %,
which suggests 1-DX was quantified on the 3-DX-quinoxaline standard (only 3-DX-quinoxaline is listed as synthesised)
[Z]: **1-DX is semi-quantitative** (as in Hou 2026 / Wei 2022). The LOD / LOQ method (S/N, or SD of blank) and whether
they refer to the derivatised solution or the original sample (1:1 dilution) are **not stated**.

## 4. Data [F]

### 4.1 Fig. 1 — 3-DX and 1-DX (mM), 120 °C, initial pH 7.5

| t (min) | 30 | 60 | 90 | 120 | 150 |
|---|---|---|---|---|---|
| 3-DX, ARP | 0.021 | **0.178** | 0.113 | 0.065 | 0.053 |
| 3-DX, +20 Xyl | 0.044 | 0.187 | 0.135 | 0.101 | 0.078 |
| 3-DX, +40 Xyl | 0.069 | 0.219 | 0.159 | 0.130 | 0.118 |
| 1-DX, ARP | **0.033** | 0.031 | 0.020 | 0.018 | 0.010 |
| 1-DX, +20 Xyl | 0.054 | 0.038 | 0.027 | 0.024 | 0.019 |
| 1-DX, +40 Xyl | 0.135 | 0.058 | 0.033 | 0.028 | 0.025 |

### 4.2 Fig. 5 — pH (cooled, 25 °C)

| t (min) | 0 | 30 | 60 | 90 | 120 | 150 |
|---|---|---|---|---|---|---|
| ARP | 7.50 | 7.12 | 6.77 | 6.35 | 6.00 | 5.75 |
| +20 Xyl | 7.50 | 6.59 | 6.02 | 5.75 | 5.26 | 4.88 |
| +40 Xyl | 7.50 | 5.78 | 5.42 | 5.07 | 4.66 | 4.28 |

### 4.3 Fig. 3 — Thr-ARP left and Thr released (mM)

| t (min) | 0 | 30 | 60 | 90 | 120 | 150 |
|---|---|---|---|---|---|---|
| Thr-ARP, ARP | 20 | 10.5 | 5.3 | 2.9 | 2.4 | 1.8 |
| Thr-ARP, +20 Xyl | 20 | 13.0 | 6.2 | 4.4 | 3.8 | 2.9 |
| Thr-ARP, +40 Xyl | 20 | 15.8 | 7.7 | 5.5 | 4.8 | 4.5 |
| Thr, ARP | 0 | 12.6 | **16.5** | 15.8 | 15.1 | 14.6 |
| Thr, +20 Xyl | 0 | 10.1 | 14.6 | 14.7 | 14.5 | 14.2 |
| Thr, +40 Xyl | 0 | 8.9 | 13.5 | 13.9 | 13.0 | 12.9 |

[!] Fig. 3b's legend swaps the colours used in every other figure (blue = +20, red = +40); the series were assigned by the
legend, which agrees with the text (more xylose → more ARP retained).

### 4.4 Fig. 2 — GO and MGO (mM; GO scale doubtful, sec. 3)

| t (min) | 30 | 60 | 90 | 120 | 150 |
|---|---|---|---|---|---|
| GO, ARP | 0.31 | 0.32 | 0.33 | 0.37 | 0.39 |
| GO, +20 Xyl | 0.42 | 0.46 | 0.48 | 0.49 | 0.50 |
| GO, +40 Xyl | 0.49 | 0.66 | 0.73 | 0.78 | 0.84 |
| MGO, ARP | 0.0206 | 0.0230 | 0.0237 | 0.0242 | 0.0280 |
| MGO, +20 Xyl | 0.0299 | 0.0304 | 0.0318 | 0.0356 | 0.0378 |
| MGO, +40 Xyl | 0.0311 | 0.0330 | 0.0344 | 0.0383 | 0.0411 |

### 4.5 Table 1 [M] — pyrazines at 150 min (µg/L), Xyl 0 / 10 / 20 / 30 / 40 mM

methylpyrazine 481.97 / 288.54 / 37.68 / 9.73 / ND; 2,5-dimethyl 651.88 / 457.17 / 314.49 / 264.11 / 185.83;
2,6-dimethyl 323.48 / 240.26 / 207.31 / 178.32 / 167.66; 2,3-dimethyl 226.06 / 136.85 / 28.96 / 4.62 / ND;
trimethyl 286.89 / 206.07 / 97.95 / 34.56 / ND. [!] The text says tetramethylpyrazine was also undetectable at 40 mM,
but the table has no tetramethylpyrazine row. No LOD for the pyrazines.

### 4.6 Arithmetic [Z]

1. **Thr-ARP loss, alone:** ln(20/C)/t = 0.0215 (30), 0.0221 (60), 0.0215 (90), 0.0177 (120), 0.0161 min⁻¹ (150) →
   **k ≈ 0.022 min⁻¹ over 0-90 min**, then slower as the pH falls below ~6.3. Hou 2026's GX-ARP (5 mM, pH 7.0, 120 °C):
   0.022 min⁻¹ over 0-20 min — the same number for a different ARP and pot.
2. **[!] Amine closure > 100 %:** Thr released / ARP consumed = 12.6 / 9.5 = **1.33 (30 min)**, 16.5 / 14.7 = **1.12
   (60 min)**, 15.8 / 17.1 = 0.92 (90 min). More free Thr than ARP lost early cannot be (Thr-ARP carries one Thr); the
   Thr and ARP calibrations (both ELSD, linear) are mutually inconsistent by ≳ 30 % at 30 min, or the ARP stock contained
   free Thr. Either way, **neither ELSD series is a reliable absolute scale** in the early window.
3. **3-DX / 1-DX within the ARP-alone pot** (cooled pH in brackets): 0.64 (7.12), 5.7 (6.77), 5.7 (6.35), 3.6 (6.00), 5.3
   (5.75). The jump between 30 and 60 min coincides with only a 0.35 pH-unit fall; no monotone ratio-vs-pH law here,
   unlike Hou 2026 (ratio 2.6 → 10.9 as pH 5.05 → 3.55). At 30 min 1-DX > 3-DX — the only point on disk from this
   design where the 1-DX side dominates.
4. **Peak deoxyosone yields (ARP alone):** 3-DX 0.178 mM = 0.9 % of the ARP (60 min); 1-DX 0.033 mM = 0.17 % (30 min).
   ~15× below Hou 2026's GX-ARP 3-DX (15.5 %) and ~30× below its 1-DX (4.8 %). Different ARP (Thr vs Gly), 4× the
   concentration, pH 7.5 vs 7.0, and a much slower pH fall here (7.50 → 7.12 in 30 min vs 7.0 → 5.05 in 20 min in Hou):
   **the same laboratory family reports deoxyosone levels an order of magnitude apart** for similar pots.
5. **Xylose adds deoxyosones**: at 30 min +40 mM Xyl gives 1-DX ×4.1 and 3-DX ×3.3 over ARP alone, while ARP loss is
   slower; the authors attribute the extra DX to xylose's own degradation [Q] — at pH 7.5 → 5.8, 120 °C, free xylose
   (40 mM) supplying 0.1 mM 1-DX in 30 min is a 0.25 % yield [Z], not tested against the engine.

## 5. What it can and cannot test

**LOD floors as transferable numbers — mostly NO.**
- **ELSD:** 0.199 µM × 8 µL = 1.6 pmol = **~0.4 ng Thr-ARP on column**. ELSD detection limits are typically tens of ng on
  column; `hu2012_extraction.md` (xylose, ELSD, S/N = 3) prints **19 mg/L = 127 µM, ~95 ng on column**, and
  `deng2022_extraction.md` (methionine, ELSD, same XBridge Amide phase) prints **LOD 89.7 µM**. Chen's ELSD LODs are
  **~200-600× lower** than either, with LINEAR calibrations where ELSD is a power law (Hu 2012: exponent 1.76-2.17;
  Deng 2022: 1.47). **Do not transfer Chen's ELSD LODs.** A defensible ELSD floor for sugars / ARPs / amino acids on this
  class of method is **~0.05-0.15 mM** (Hu, Deng) — owner call before any of it becomes a censoring rule.
- **OPD quinoxalines (DAD 315 nm):** 0.47-0.59 µM is ~5-20× below Deng 2022's same-family values (GO 8.8, MGO 3.9, DA
  2.9 µM) and plausible for a UV-active quinoxaline at 5 µL injection [Z]. **Usable as an optimistic floor (~0.5 µM in
  the derivatised solution, ~1 µM in the sample after the 1:1 dilution) for 3-DX / GO / MGO from this method family; a
  conservative floor is Deng's 3-9 µM.** The LOQs (~2 µM) are the better censoring threshold.
- **Furfural UV:** 0.53 / 1.44 µM; plausible.

**1-DX / 3-DX vs time:** yes (Fig. 1, five times, three arms); **vs pH:** only as within-pot co-variation with the
measured (cooled) pH trajectory; no designed pH series.

**Can also:** give a second measured pH trajectory for an unbuffered Jiangnan-style ARP pot (with Hou 2026) — here the
drift is ~5× slower in the first 30 min; give Thr-ARP lifetime at 120 °C; show xylose re-forming ARP (¹³C5 label;
qualitative MS only — no labelled-ARP quantity).

**Cannot:** separate pH from xylose effects (both change together); give absolute GO / MGO (scale doubt); give amine
release reliably (closure > 100 %); give any pyrazine LOD; give O2 status.

## 6. Conflicts and agreements with the repo's records

1. **[!] Deoxyosone levels ×15-30 below `hou2026` / `wei2022`** for a similar ARP pot of the same method family (4.6
   item 4). Same family as CF-060 / CF-062 (ARP-pot scale spread; semi-quantitative DX scale). Not registered here.
2. **[!] pH drift ~5× slower than `hou2026`** in the first 30 min (Δ 0.38 vs 1.95 units), opposite to what a 4× more
   concentrated ARP would suggest; the measured Thr release (a zwitterion, little buffering at pH 6-7 [Z]) does not
   obviously explain it. Supports hub H7 that every unbuffered "pH 7" ARP pot drifts, but shows the drift rate is
   pot-specific: **an in-situ-pH operator cannot be transferred from Hou's record to other ARPs**.
3. **[!] Calibration-unit inconsistency** (sec. 3) and **amine closure > 100 %** (4.6 item 2).
4. `hou2026`'s apparent ARP k (0.022 min⁻¹ at 120 °C) is matched here — an agreement on ARP lifetime, not on products.

## 7. Proposed ledger entries (trunk.yml / aroma_other.yml; NOT written — first reader only, second reader pending)

```yaml
  - id: chen2026_fig1_fig5_deoxyxylosones_ph_thr_arp_xylose_120c
    dossier: chen2026_extraction.md
    paper: "Chen, Hu, Dou, Meng & Deng 2026 Curr Res Food Sci 12:101433"
    lab: "Chuzhou / Enshi (Meng F / Deng Q; Jiangnan methods)"
    system: "threonine-xylose Amadori compound 20 mM, + D-xylose 0 / 20 / 40 mM"
    medium: "water, no buffer, initial pH 7.5 (NaOH); pressure-resistant glass bottles, oil bath"
    temperature_c: [120]
    ph: {initial: [7.5], measured_during: true, final: [5.75, 4.88, 4.28]}
    time_min: [0, 30, 60, 90, 120, 150]
    species: ["3-deoxyxylosone", "1-deoxyxylosone", glyoxal, methylglyoxal, "Thr-Amadori", threonine]
    kind: time_series
    location: "Figs. 1, 2, 3, 5 pp. 5-7 (digitised in dossier sec. 4)"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    read_by: "first reader 2026-11-05; second reader pending"
    condition_quote: "All previously mentioned solutions were brought to a pH of 7.5"
    flags: [unbuffered, ph_measured_cooled, o2_not_stated, dx_semi_quantitative, go_mgo_calibration_units_inconsistent, amine_closure_above_100pct, elsd_lod_implausible]
    note: "pH measured on cooled samples at every time. 3-DX peak 0.9 % of ARP (60 min), 1-DX 0.17 % (30 min), x15-30 below hou2026. ARP k ~0.022 1/min 0-90 min."
  - id: chen2026_methods_lod_loq_elsd_opd_furfural
    dossier: chen2026_extraction.md
    paper: "Chen, Hu, Dou, Meng & Deng 2026 Curr Res Food Sci 12:101433"
    lab: "Chuzhou / Enshi (Meng F / Deng Q; Jiangnan methods)"
    system: "method validation: Thr-ARP, Thr (HPLC-ELSD); 3-DX, 1-DX, GO, MGO (OPD, HPLC-DAD 315 nm); furfural (HPLC-UV 277 nm)"
    medium: "standards (matrix not stated)"
    temperature_c: []
    ph: {initial: [], measured_during: false, final: null}
    time_min: []
    species: ["Thr-Amadori", threonine, "3-deoxyxylosone", "1-deoxyxylosone", glyoxal, methylglyoxal, furfural]
    kind: other
    location: "Methods secs. 2.4-2.6 pp. 2-3 (dossier sec. 3)"
    status: text_only
    quant_class: other
    simulable: false
    read_by: "first reader 2026-11-05; second reader 2026-11-05 AGREES on the OPD LODs and LOQs (see Second reader section)"
    condition_quote: "detection limit (LOD) and quantitation limit (LOQ) values were determined to be 0.199"
    flags: [lod_method_not_stated, elsd_linear_calibration, elsd_lod_implausible_vs_hu2012_deng2022]
    note: "LOD/LOQ uM: Thr-ARP 0.199/0.601, Thr 0.421/1.223, 3-DX 0.59/2.03, 1-DX 0.58/2.29, GO 0.57/2.32, MGO 0.47/2.08, furfural 0.531/1.438. ELSD values ~200-600x below hu2012/deng2022: do not use as floors."
```

## 8. Supplementary material

Not on disk: Fig. S1 (furfural vs time, all five xylose levels). The 10 and 30 mM xylose arms' time courses are not in the
paper at all.

## Second reader (2026-11-05)

**Independent second reading** by `agent-second-reader-2026-11-05`, PDF read before this dossier. Source:
`data/articles/chen2026.pdf` (Curr. Res. Food Sci. 12:101433), sec. 2.5 (OPD quinoxalines, HPLC-DAD 315 nm). Printed numbers.

| analyte | first reader LOD / LOQ (µM) | second reader LOD / LOQ (µM) | calibration as printed (both readers) | verdict |
|---|---|---|---|---|
| 3-DX | 0.59 / 2.03 | **0.59 / 2.03** | y = 4.79 × 10⁶ x − 2.19 × 10⁴ | AGREE |
| 1-DX | 0.58 / 2.29 | **0.58 / 2.29** | y = 4.78 × 10⁶ x − 3.22 × 10⁴ | AGREE |
| GO | 0.57 / 2.32 | **0.57 / 2.32** | y = 2.47 × 10³ x − 7.83 × 10⁴ | AGREE |
| MGO | 0.47 / 2.08 | **0.47 / 2.08** | y = 3.46 × 10³ x − 1.44 × 10⁴ | AGREE |

Range 0.47-0.59 µM: AGREE. I see the same defects the first reader flags: x is printed as peak area and y as mmol/L, which
cannot be right for 4.79 × 10⁶; the 3-DX / 1-DX and GO / MGO slopes differ by 10³; the LOD criterion (S/N or blank SD) and
the basis (pot or derivatised solution, diluted 1:1 with reagent) are not printed. I also note that sec. 2.5 gives the
column temperature as 30 °C and, a few lines later, the oven as 35 °C.
