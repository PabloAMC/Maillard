# Wang et al. 2026 — Cys-Amadori / Glu-Amadori from cysteine–xylose–glutamic acid: a FIVE-RUNG thiol ladder

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★★ THE CORPUS'S FIRST INDEPENDENT MFT/FFT TEMPERATURE LADDER — five rungs, a different lab, a different pot, and it shows the OPPOSITE of the Kang switch-on: both thiols peak in mid-ladder and FALL at high temperature.

---

## §0. IDENTITY `[M]` — **the file is what the brief expected**

| item | value | anchor |
|---|---|---|
| file | `data/articles/Wang2026.pdf`, 6 358 283 B, 11 pages, Elsevier/Distiller, born-digital | `pdfinfo` |
| title | **"Impact of Amadori rearrangement products of cysteine-xylose-glutamic acid on meat flavor development"** | p. 1 |
| authors | **Xiao-Lan Wang**, **Chao-Kun Wei\***, Dai-Wei Wei, Chun-Mei Yang, **Yuan Liu\*** | p. 1 |
| affiliations | School of Food Science and Engineering, **Ningxia University**, Yinchuan; School of Agriculture & Biology, **Shanghai Jiao Tong University**; Ningxia Ningyang Food Co. Ltd | p. 1 |
| **DOI** | ★ **`10.1016/j.foodres.2026.118535`** ✔ exactly as the brief expected | p. 1 |
| journal | ***Food Research International* 230 (2026) 118535** | running head |
| dates | Received **8 Oct 2025** · Revised **24 Jan 2026** · Accepted **1 Feb 2026** · Online **5 Feb 2026** | p. 1 |
| SI | **NOT ON DISK.** Table S1 (MS fragments), **Table S2 — the 31-compound volatile data**, Figs S1–S2 | referenced pp. 3–5 |

★ **This is a DIFFERENT LABORATORY from the Zhai/Kang chain** (Ningxia + SJTU, no author overlap with Jiangnan/Zhejiang Gongshang), which is why its ladder matters: it is the corpus's **only** independent MFT/FFT temperature series.

---

## §1. SYSTEM AND CONDITIONS `[M]`

### 1.1 The two intermediates (prepared, purified, then used as substrates)
| | Glu-Amadori | Cys-Amadori |
|---|---|---|
| recipe | glutamic acid **1.47 g** + xylose **3 g** in **170 mL** water, pH 7 (NaOH 6 M) | xylose **3.3 g** + L-cysteine·HCl·H₂O **3.5 g** + Na₂HPO₄ **7.2 g** in **5 mL** water |
| two-stage heating | **90 °C / 150 min → 130 °C / 90 min** | **50 °C / 60 min → 70 °C / 180 min** |
| work-up | ice quench, −80 °C 24 h, freeze-dry, 4 °C | same |
| purification | rotary evaporation, ethanol extraction, **Dowex 50WX4 (H⁺)**, 1 M ammonia elution, semi-prep RP-HPLC, **C18 10 μm 22 × 200 mm**, lyophilised | same |
| identity | **C₁₀H₁₇NO₈, m/z 280** (MS² 262/244/226/180/148/130/102/84) | **C₈H₁₅NO₆S, m/z 254** |
| ⚠ purity | **not reported for either** `[NEG]` | |

### 1.2 The model reaction — **the ladder**
| parameter | value |
|---|---|
| charge | ★ **1.0 g Glu-Amadori + 1.0 g Cys-Amadori** (both, together) in **20.0 mL PBS** |
| ★ buffer | ★ **PBS 0.2 mol/L** — **BUFFERED**, unlike every other sulfur paper in the corpus except Zhang 2026 and Chan 1994 |
| pH | **5, 6, 7, 8, 9** |
| **temperature** | ★ **85, 95, 105, 115, 125 °C** — **FIVE RUNGS** |
| time | **80, 90, 100, 110, 120 min** (Methods) ⚠ **but Fig. 7 plots 60, 80, 100, 120, 140 min** — see §5.1 |
| vessel | thick-walled pressure-resistant bottle, **precision oil bath, continuous magnetic stirring** |
| quench | ice bath |
| replication | **n = 3**, mean ± SD, one-way ANOVA, SPSS 20.0, p < 0.05 |
| ⚠ **the fixed levels of the other two factors are NEVER STATED** | the design is one-factor-at-a-time, but the paper never says at what pH and time the temperature series was run, nor at what temperature the pH and time series were run. **`[!]` This is the single largest defect in the paper for modelling purposes** — the ladder's other coordinates are unknown |

**Inference (marked as inference, not fact):** the Conclusions say *"At a reaction temperature of **105 °C** it is beneficial to promote the formation of 2-acetylthiazole"* and *"when the reaction time is between **100 and 140 min**…"*, and Fig. 2's pH panel and Fig. 7's time panel show mutually consistent magnitudes with Fig. 5's 105 °C column. **Best guess: the pH and time series were run at ~105 °C, and the temperature series at pH 7 for 120 min. NOT VERIFIED.**

---

## §2. ANALYTICAL METHOD AND QUANTIFICATION BASIS

| parameter | value |
|---|---|
| method | HS-SPME−GC−MS |
| sample | ★ **0.5 g of LYOPHILISED sample** (not solution) |
| **internal standard** | **1,2-dichlorobenzene, 2 μL of 50 μg/mL in methanol → 0.1 μg absolute** |
| fibre | **CAR/PDMS 50/35 μm, 2 cm** |
| conditions | equilibrate **50 °C / 10 min**, adsorb **20 min** |
| column | ★ **HP-5 MS UI 30 m × 250 μm × 0.25 μm — a NON-POLAR column**, unlike the DB-Wax used by the whole Zhai/Kang/Feng family |
| oven | 40 °C (2 min) → 2 °C/min → 100 °C → 4 °C/min → 150 °C (2 min) → 20 °C/min → 280 °C |
| MS | EI 70 eV, source 230 °C, quad 150 °C, **full scan 40–450 m/z** |
| ID | **NIST 08 library + retention indices (C₇–C₃₀)** |
| **quantification** | ⚠ **NEVER STATED.** The Methods describe identification only. No calibration curve, no response factor, no external standard, no equation. Values are printed as **μg/L** |

### 2.1 ★ THE VERDICT ON THE QUANTIFICATION BASIS

> **⚠ THE QUANTIFICATION BASIS IS NOT DESCRIBED. Given an internal standard of known mass and no
> calibration curves anywhere, this is almost certainly SINGLE-IS SEMI-QUANT with the response
> factor assumed = 1 — the same basis as Zhai 2022 — but the paper does not say so, and that
> inference must be labelled as an inference.**

⚠️ **And a unit problem the paper never addresses:** the internal standard is spiked into **0.5 g of
freeze-dried solid**, yet the results are reported in **μg/L**. There is no stated dilution or
reconstitution volume connecting a mass basis to a volume basis. **`[!]` The absolute μg/L scale is
therefore not reconstructible.**

**★ THE WAVE RULE, APPLIED.** A constant unknown response factor and a constant unknown mass→volume
conversion both cancel in a ratio, and therefore in an Arrhenius slope and in any within-figure
fold-change. **⇒ Every temperature, pH and time RATIO in this paper is usable; no absolute value is.**
This paper is `RATIO-ONLY` throughout on the volatile side.

★ **The ARP data are on a still weaker footing: Fig. 5B and 7B report LC-MS PEAK AREA, not
concentration** — explicitly labelled "Peak area" on the y-axis. Those are ratio-only by construction.

---

## §3. ★★ FIGURE 5 — THE FIVE-RUNG TEMPERATURE LADDER `[D]` (digitised)

**Digitisation method:** page 7 rendered at **900 dpi**; the panel-C y-axis was calibrated from its
tick grid (**12.443 px per μg/L**, 220 → 0 over 2 737.5 px) and bar tops were located by matching
each series' exact fill RGB (85 °C `#FDC798`, 95 °C `#9DD79D`, 105 °C `#C2B3D4`, 115 °C `#FEFEA4`,
125 °C `#9BBBF4`) and walking up each bar's column from the baseline.
**Quoted precision ± 5 % of reading, with a hard floor: bars below ≈ 2.5 μg/L are at the limit of
this method and several small ones were missed by the automatic pass and recovered by direct
inspection of a 900 dpi crop. Those are marked ~.**
⚠️ **The authoritative numbers are in Table S2, which is not on disk. Everything in §3.1 is a
digitisation and should be replaced when the SI is obtained.**

### 3.1 Panel C — content of the eleven main flavour compounds, μg/L

| compound | **85 °C** | **95 °C** | **105 °C** | **115 °C** | **125 °C** | shape |
|---|---:|---:|---:|---:|---:|---|
| Pyrazine | n.d. | n.d. | ~1 | **2.0** | **3.9** | ↑ monotone, on above ~105 °C |
| Methylpyrazine | n.d. | n.d. | ~1 | **1.8** | **3.0** | ↑ monotone |
| 2-Methylthiophene | ~1 | **4.9** | **10.1** | **13.1** | **17.0** | ★ ↑ monotone, ×17 over the ladder |
| **Furfural** | **27.6** | **67.3** | **136.3** | **150.2** | **155.5** | ★ ↑ **saturating** — the classic plateau |
| **2-Methyl-3-furanthiol (MFT)** | n.d. | ~1 | **2.5** | **6.6** | **~2** | ★★ **PEAKS at 115 °C, then falls** |
| **2-Furfurylthiol (FFT)** | **4.5** | **12.0** | **6.7** | **~1** | n.d. | ★★ **PEAKS at 95 °C, then collapses to zero** |
| 2-Thiophenethiol | **7.8** | **6.2** | **6.7** | **5.7** | **4.4** | ★ **monotonically DOWN** |
| 2-Thiophenecarboxaldehyde | n.d. | n.d. | ~1 | **7.8** | **14.1** | ↑ |
| **2-Acetylthiazole** | **92.1** | **140.2** | **183.4** | **119.8** | **155.4** | ⚠ **non-monotone; peaks at 105 °C** |
| Nonanal | ~3 | **2.6** | ~4 | **6.9** | **9.4** | ↑ (a lipid-oxidation marker, see §6) |
| Difurfuryl disulfide | **2.7** | ~1 | **1.5** | **2.1** | **5.3** | ↑ noisy |

### 3.2 Panels A and B `[D]`
| T (°C) | **A₄₂₀** | **Cys-Amadori peak area** | **Glu-Amadori peak area** |
|---:|---:|---:|---:|
| 85 | **0.44** | 1.3 × 10⁵ | 1.6 × 10⁵ |
| 95 | **0.67** | 3.25 × 10⁵ | 1.9 × 10⁵ |
| 105 | **0.805** | ★ **6.2 × 10⁵ ← max** | 2.4 × 10⁵ |
| 115 | **0.94** | 5.5 × 10⁵ | ★ **4.9 × 10⁵ ← max** |
| 125 | ★ **1.35** | 3.1 × 10⁵ | 4.5 × 10⁵ |

★ **The two Amadori compounds have DIFFERENT optimum temperatures — Cys-Amadori at 105 °C,
Glu-Amadori at 115 °C — in the same pot.** That is a clean, useful structural fact: the sulfur-bearing
Amadori is the less thermally stable of the two by ~10 K.

`[D]` **Browning Eₐ over the five rungs**, ln A₄₂₀ vs 1/T: **Eₐ = 30.6 kJ mol⁻¹, R² = 0.971** —
★ **within 8 % of the 28.3 kJ mol⁻¹ (R² 0.969) that `zhai2023jafc_extraction.md` §5.1 derives from a
completely different lab's seven-rung A₄₂₀ ladder in an unbuffered TTCA pot.** Two independent
browning ladders, two labs, two systems, two buffer states: **28.3 and 30.6 kJ mol⁻¹.**
**That is the best-replicated activation energy in the K6a corpus** — and note that it replicates
*across* the buffered/unbuffered divide that §4 blames for the thiol disagreement, which is itself
informative: browning is far less pH-trajectory-sensitive than the free-thiol channel.

---

## §4. ★★ THE KANG SWITCH-ON CROSS-CHECK — **REFUTED IN AN INDEPENDENT LABORATORY**

Kang's claim (`kang2026_SI_extraction.md` §7a): MFT and FFT are flat from 100→120 °C (×1.12, ×1.10)
and then **jump** to 140 °C (×4.26, ×2.78), with apparent Eₐ climbing from ~6 to 70–98 kJ mol⁻¹,
**against** a sulfur class that decelerates.

Wang's ladder covers **85–125 °C**, i.e. it brackets Kang's low leg completely and reaches to within
15 K of Kang's top rung. `[D]` leg-by-leg fold changes:

| species | ×(85→95) | ×(95→105) | ★ ×(105→115) | ★ ×(115→125) |
|---|---:|---:|---:|---:|
| **MFT** | — | ≈ ×2.5 | ★ **×2.6** | ★ **×0.30 ↓** |
| **FFT** | **×2.7** | ★ **×0.56 ↓** | ★ **×0.15 ↓** | ★ **→ 0** |
| 2-Thiophenethiol | ×0.79 ↓ | ×1.08 | ×0.85 ↓ | ×0.77 ↓ |
| Furfural | ×2.44 | ×2.03 | ×1.10 | ×1.03 |
| 2-Methylthiophene | ×4.9 | ×2.05 | ×1.30 | ×1.30 |
| 2-Acetylthiazole | ×1.52 | ×1.31 | ×0.65 ↓ | ×1.30 |
| Pyrazine | — | — | ×2.0 | ×1.9 |

### ★★ THE VERDICT
1. **No switch-on.** Over 105→125 °C — the interval that contains Kang's claimed threshold from
   below — **MFT falls by 3.3× and FFT falls to below detection.** Whatever happens to free thiols
   between 105 and 125 °C in a buffered Cys-Amadori pot, it is not a step increase.
2. **The free thiols behave OPPOSITELY to their own precursor.** Furfural — FFT's direct precursor
   on the paper's own Fig. 3B scheme — rises monotonically and **saturates** (×2.44, ×2.03, ×1.10,
   ×1.03), while FFT collapses. **⇒ The loss is on the thiol, not on the precursor supply. This is
   a SINK observation, and it is the second independent one this wave produced** (the first being
   Zhai's 140 °C FFT turnover, `zhai2023jafc_extraction.md` §8.5).
3. **The one thing that does replicate is the CLASS-level shape** — the ring-closure products
   (2-methylthiophene, 2-thiophenecarboxaldehyde, pyrazines) rise monotonically while the free
   thiols do not. Kang's contrast between "free thiols" and "the sulfur class" is real; **its sign
   at the top of the ladder is not.**
4. ⚠️ **The honest caveat.** Wang's system is a *purified-Amadori* pot in **0.2 M phosphate buffer**
   with **no free xylose and no free cysteine**, and its thiol precursors are limited by what the
   two ARPs carry. Kang/Zhai's is an unbuffered TTCA pot whose pH crashes to 4.9. **These are not
   the same chemistry, and the disagreement could be a pH-trajectory effect rather than a
   temperature effect.** That is itself the most useful reading: **the buffered pot does not show
   the switch-on and the unbuffered one does, which points at pH, not at a barrier.**

---

## §5. THE TIME AND pH AXES

### 5.1 ★ Figure 7 — the TIME series `[D]` (digitised, ±10 %), μg/L
⚠️ **Methods say 80/90/100/110/120 min; Figure 7 plots 60/80/100/120/140 min. The figure is used.**
⚠️ **The temperature and pH of this series are not stated (§1.2).**

| compound | 60 min | 80 min | 100 min | 120 min | 140 min |
|---|---:|---:|---:|---:|---:|
| Pyrazine | ~1 | 2 | 5.5 | 7 | **15** |
| Methylpyrazine | n.d. | n.d. | 3.5 | 4.5 | **7** |
| 2-Methylthiophene | n.d. | 4 | **138** | 148 | **173** |
| **Furfural** | 75 | **57** | 108 | **178** | **103** | ⚠ non-monotone, saw-toothed |
| **MFT** | 6 | 14 | ★ **20** | 17 | 18 |
| **FFT** | ~1 | **5** | **5** | ~1 | ~1 | ★ **peaks at 80–100 min, then collapses** |
| 2-Thiophenethiol | 9 | 14 | ★ **18** | 16 | 15 |
| 2-Thiophenecarboxaldehyde | ~1 | ~1 | ~1 | ~1 | 2.5 |
| **2-Acetylthiazole** | 69 | 62 | ★ **76** | 70 | **43** | ★ falls 1.8× after 100 min |
| Nonanal | 9 | 9 | 11 | 12 | 14 |
| Difurfuryl disulfide | ~0.5 | ~0.5 | 2 | ~1 | ~1 |
| **A₄₂₀** (panel A) | **0.175** | **0.19** | **0.285** | **0.48** | **0.495** |
| **Cys-Amadori peak area** (panel B) | 2.5e5 | 2.55e5 | 3.35e5 | ★ **7.5e5** | 6.0e5 |
| **Glu-Amadori peak area** | 1.55e5 | 4.3e5 | ★ **5.95e5** | 4.35e5 | 3.25e5 |

★★ **THE TIME AXIS SHOWS THE SAME TURNOVER AS THE TEMPERATURE AXIS, AND THIS IS DIRECTLY THE
YILTIRAK FAILURE MODE.** FFT peaks at 80–100 min and is back to its 60-min value by 120 min;
2-acetylthiazole peaks at 100 min and falls 1.8× by 140 min; MFT peaks at 100 min.
**A model with formation-only kinetics predicts monotone growth here and will over-predict every
long-hold point** — which is exactly what `cutover_final_exam.md` reports for the Yiltirak family
(0/8, over-prediction at 4/4 rungs, worst at long holds).

### 5.2 The pH axis (Fig. 2, **not digitised** — reported from the authors' own text) `[M]`
> *"with increasing pH, the content of the main sulfur-containing compounds, such as
> **2-furfurylthiol, 2-methyl-3-furanthiol and 2-thiophenethiol, gradually decreased**. … The
> content of ARPs showed the opposite trend."* (p. 6)
> *"the content of methylpyrazine and pyrazine in the reaction system was relatively low, but the
> overall trend **increased with the pH**. … **Pyrazine was not detected at pH 5 or 6.**"* (p. 6)
> *"the content of furfural **decreased significantly with the increase of the pH** … especially
> under alkaline conditions with pH > 7, where its content dropped sharply."* (p. 6)
> *"at pH < 7, the reaction between xylose and cysteine mainly produced **2-furanethiol**, while at
> pH > 7, it inclined to form **2-methyl-3-furanthiol**."* (p. 7)

★ **The last sentence is a pH-dependent BRANCHING claim between FFT and MFT — the exact quantity the
repo carries as a fixed FFT:MFT ratio.** It is stated qualitatively and its data are in Fig. 2C.
`[NEG]` **no numbers are printed in the text; digitising Fig. 2C is a declared gap (§9).**

⚠️ **This contradicts the Kang finding** that FFT:MFT is *invariant to pH* (2.73–3.02 over pH 5.5–8,
`kang2026_SI_extraction.md` §7b). **Two papers, opposite claims, different systems.** Recorded as an
open contradiction, not adjudicated.

---

## §6. STRUCTURAL AND MECHANISTIC CLAIMS `STRUCTURAL`

| # | claim | evidence in this paper |
|---|---|---|
| 1 | ★ **Cys-Amadori and Glu-Amadori have different thermal optima (105 vs 115 °C) in the same pot** | Fig. 5B |
| 2 | ★ **Free thiols peak and decline on BOTH the temperature and the time axis while their furan precursor saturates** | Figs. 5C, 7C |
| 3 | **Pyrazine formation switches on above pH 6 and above ~105 °C**; not detected at pH 5–6 at any temperature | Fig. 2C, 5C |
| 4 | **2-Acetylthiazole formation is "suppressed by high temperature"** — the authors annotate their own Fig. 6A scheme with this in red, and Fig. 5C shows the 115 °C dip | Fig. 6A, 5C |
| 5 | FFT from furfural + H₂S; MFT from 1-deoxypentosone → 4-hydroxy-5-methyl-3(2H)-furanone + H₂S; thiophene from mercaptoacetaldehyde | Fig. 3B, 6B — **the same scheme Zhai's Fig. 6 model reactions verified experimentally** |
| 6 | ⚠ **Nonanal rises monotonically with temperature (2.6 → 9.4 μg/L)** in a system with **no lipid**. Nonanal is the corpus's lipid-oxidation marker | Fig. 5C — ★ **either an artefact (septum/fibre bleed, plasticiser) or evidence that nonanal has a non-lipid route. Worth flagging to the lipid lane** |

---

## §7. VERIFIED NEGATIVES `[NEG]`

| # | sought | verdict |
|---|---|---|
| 1 | Any rate constant, Eₐ, half-life, or kinetic model | **ABSENT.** The paper reports no kinetics. (A sentence on p. 8 about *"rate constants for the final 17 reactions … exponential growth with increasing temperature"* is a **citation to other work**, not this paper's result) |
| 2 | The quantification basis for the volatiles | ★ **ABSENT** (§2.1) |
| 3 | The fixed levels of the other two factors in each one-factor series | ★ **ABSENT** (§1.2) — the most damaging omission |
| 4 | ARP concentrations | **ABSENT — peak area only** |
| 5 | ARP purity | **ABSENT** |
| 6 | Measured pH after heating | **ABSENT** (the pot is buffered at 0.2 M, so drift should be small — but it is not measured) |
| 7 | H₂S measurement | **ABSENT** |
| 8 | α-Dicarbonyl measurement | **ABSENT** — the paper's mechanism is entirely inferential on the dicarbonyl side |
| 9 | Sensory or odour-activity values | **ABSENT** |
| 10 | Any temperature outside 85–125 °C | **ABSENT** |
| 11 | Numeric volatile table in the main text | **ABSENT — Table S2, not on disk** |

---

## §8. CONSOLIDATED PARAMETER TABLE

**Common conditions:** 1.0 g Cys-Amadori + 1.0 g Glu-Amadori in 20 mL **0.2 M PBS**, pressure bottle,
oil bath with stirring, ice quench, HS-SPME (CAR/PDMS 50/35 μm, 2 cm; 50 °C/10 min + 20 min) GC-MS on
**HP-5 MS UI**, IS = 1,2-dichlorobenzene 0.1 μg, n = 3, **μg/L (basis undescribed — §2.1)**.

| # | parameter | value | units | condition | prov | anchor |
|---:|---|---|---|---|---|---|
| 1 | ★ **MFT temperature ladder** | **n.d. / ~1 / 2.5 / 6.6 / ~2** | μg/L | 85/95/105/115/125 °C | **D** (digitised) | Fig. 5C |
| 2 | ★ **FFT temperature ladder** | **4.5 / 12.0 / 6.7 / ~1 / n.d.** | μg/L | " | **D** | Fig. 5C |
| 3 | 2-Thiophenethiol ladder | 7.8 / 6.2 / 6.7 / 5.7 / 4.4 | μg/L | " | **D** | Fig. 5C |
| 4 | Furfural ladder | 27.6 / 67.3 / 136.3 / 150.2 / 155.5 | μg/L | " | **D** | Fig. 5C |
| 5 | 2-Methylthiophene ladder | ~1 / 4.9 / 10.1 / 13.1 / 17.0 | μg/L | " | **D** | Fig. 5C |
| 6 | 2-Acetylthiazole ladder | 92.1 / 140.2 / 183.4 / 119.8 / 155.4 | μg/L | " | **D** | Fig. 5C |
| 7 | **A₄₂₀ ladder** | **0.44 / 0.67 / 0.805 / 0.94 / 1.35** | — | " | **D** | Fig. 5A |
| 8 | ★★ **Browning Eₐ** | ★ **30.6** (R² 0.971) | kJ mol⁻¹ | 85–125 °C, buffered | **D** | §3.2 — ★ **replicates Zhai's 28.3 from another lab** |
| 9 | **Cys-Amadori thermal optimum** | ★ **105 °C** | — | peak area basis | **M** | Fig. 5B |
| 10 | **Glu-Amadori thermal optimum** | ★ **115 °C** | — | " | **M** | Fig. 5B |
| 11 | **FFT time optimum** | ★ **80–100 min**, then collapse | — | T unstated | **D** | Fig. 7C |
| 12 | **MFT time optimum** | ★ **100 min** | — | " | **D** | Fig. 7C |
| 13 | 2-Acetylthiazole time optimum | **100 min**, falls 1.8× by 140 min | — | " | **D** | Fig. 7C |
| 14 | Cys-Amadori time optimum | **120 min** | — | peak area | **D** | Fig. 7B |
| 15 | Glu-Amadori time optimum | **100 min** | — | " | **D** | Fig. 7B |
| 16 | A₄₂₀ vs time | 0.175 / 0.19 / 0.285 / 0.48 / 0.495 | — | 60–140 min | **D** | Fig. 7A |
| 17 | **pH response, thiols** | FFT, MFT, 2-thiophenethiol **all decrease** with pH 5→9 | — | qualitative | **M** | p. 6 |
| 18 | **pH response, pyrazines** | **increase**; **not detected at pH 5 or 6** | — | " | **M** | p. 6 |
| 19 | **pH response, furfural** | **decreases**, sharply above pH 7 | — | " | **M** | p. 6 |
| 20 | ⚠ **FFT:MFT branching vs pH** | *"at pH < 7 … mainly 2-furanethiol; at pH > 7 … 2-methyl-3-furanthiol"* | — | qualitative | **M** | p. 7 — ★ **contradicts Kang's pH-invariance** |
| 21 | Compound count | **31 volatiles**: 14 S, 3 pyrazines, 4 thiazoles+thiophenes, 6 O | — | across 5 pH levels | **M** | p. 5 |

---

## §9. USABILITY VERDICTS

### `USE-Q`
- **Browning Eₐ 30.6 kJ mol⁻¹** — five rungs, buffered, R² 0.971; qualified only by the unstated
  fixed pH/time and the 5× dilution used for the A₄₂₀ read.

### `RATIO-ONLY` — **the default verdict for this paper**
Every volatile number (§8 rows 1–6, 11–13). The quantification basis is undescribed and the
mass→volume conversion is unreconstructible; **only fold-changes within a figure transfer.**

### `STRUCTURAL` — the durable content
1. ★ **MFT and FFT peak and decline across 85–125 °C, while furfural saturates.**
2. ★ **The same peak-and-decline appears on the time axis at 80–140 min.**
3. **Cys-Amadori (105 °C) is less thermally stable than Glu-Amadori (115 °C).**
4. **Pyrazines require pH > 6 and T ≳ 105 °C.**
5. **2-Acetylthiazole is suppressed at high temperature** (the authors' own annotation).

### `PRIOR-ONLY`
- The pH-dependent FFT↔MFT branching claim (§8 row 20) — qualitative, contradicts Kang, and its
  data are in a figure this dossier did not digitise.

### ★ `REFUSE`
| item | reason |
|---|---|
| **Any absolute μg/L value** | basis undescribed; 0.5 g solid reported as μg/L with no volume (§2.1) |
| **Any ARP magnitude** | peak area, not concentration |
| **Cross-comparison of these μg/L to Zhai/Kang/Feng μg/L** | different column (HP-5 vs DB-Wax), different fibre, different matrix (lyophilisate vs solution), undescribed basis |
| **Any Eₐ from the MFT/FFT ladders** | ★ **non-monotone in T — no Arrhenius form exists.** This is a refusal with content: the refusal *is* the finding |

---

## §10. DRAFT FIT / HOLD-OUT ROLES — **FOR THE ORCHESTRATOR, NOT A DECLARATION EDIT**

### 10a. Recommended **HOLD-OUT** (strong)
| candidate | what it tests | why it is a real test |
|---|---|---|
| ★ **"FFT declines from 105 to 125 °C while furfural saturates"** — as a SIGN test, not a magnitude test | the thiol sink's temperature dependence | the current core has no mechanism that can turn FFT over while its precursor rises; a pass would be evidence, a fail is diagnostic |
| ★ **"FFT peaks at 80–100 min and returns to its 60-min level by 120 min"** — sign + peak-time test | the time axis | ★ **this is the Yiltirak failure in miniature, in a clean two-substrate pot** |
| **"Cys-Amadori peaks at 105 °C, Glu-Amadori at 115 °C"** | differential ARP stability | the repo carries one Amadori stability parameter; this says there are at least two |
| **"Pyrazine = 0 at pH 5 and 6 at every temperature"** | the amine/dicarbonyl gate | a hard on/off |

### 10b. **NOT suitable as FIT rows**
Nothing here should be fitted: no absolute scale, no stated fixed conditions, no numeric table.

### 10c. **DO NOT USE**
- Absolute μg/L for calibration.
- Any Arrhenius fit to the thiol ladders.

---

## §11. DECLARED GAPS

| # | gap | closure |
|---|---|---|
| G1 | ★ **Table S2 — the 31-compound numeric volatile data across five pH levels** — is the paper's only real table and it is not on disk. **Every number in §3.1 and §5.1 is a digitisation standing in for it.** | order the SI of `10.1016/j.foodres.2026.118535` — ★ **highest-value SI request of this wave** |
| G2 | ★ **The fixed levels of the untested factors in each one-factor series are unstated**, so the ladder's other coordinates are unknown | ask the authors; not recoverable from the paper |
| G3 | **The quantification basis is undescribed** | SI may state it; otherwise the paper is permanently ratio-only |
| G4 | **Fig. 2C (the pH panel) was not digitised** in this pass; it carries the FFT↔MFT branching claim that contradicts Kang | digitise Fig. 2C, or use Table S2 when obtained |
| G5 | **No α-dicarbonyl, no H₂S, no pH-after-heating** | not obtainable |
| G6 | ⚠ **Nonanal in a lipid-free system** rising 3.6× with temperature — unexplained | flag to the lipid lane; check whether the repo's nonanal anchors assume a lipid-only route |
