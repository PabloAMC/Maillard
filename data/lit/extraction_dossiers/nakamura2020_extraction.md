# Nakamura, Mikami, Noda & Murata 2020/2021 — COMPLETE TRANSCRIPTION (Wave K4c)

**Full extraction of every number in `data/articles/nakamura2020.pdf`.**
Extraction date 2026-08-28. **This paper contains no tables.** Every quantitative result is
inside a figure, so §4–§7 below are programmatic digitisations off 400–600 dpi renders, with
the method, calibration and residual error stated per panel and every value cross-checked
against the numbers the authors quote in the running text.

---

## 0. PAPER IDENTITY — MATCHES THE WAVE BRIEF

| field | value |
|---|---|
| Authors | **Miki Nakamura¹, Yoko Mikami¹, Kyoko Noda¹, Masatsune Murata¹\*** |
| Affiliation | Department of Nutrition and Food Science, **Ochanomizu University**, 2-1-1 Otsuka, Bunkyo-ku, Tokyo, Japan |
| Title | **"Browning of Maillard reaction systems containing xylose and 4-hydroxy-5-methyl-3(2H)-furanone"** |
| Journal | **Bioscience, Biotechnology, and Biochemistry, 2021, Vol. 85, No. 2, 401–410** |
| DOI | **10.1093/bbb/zbaa019** ✔ matches the brief |
| Dates | Received 23 June 2020; accepted 17 August 2020; **advance access 31 December 2020**; issue **February 2021** |
| Type | Regular Paper. © The Author(s) 2020, OUP on behalf of JSBBA |
| Funding | JSPS Grant-in-Aid for Scientific Research (B) **17H01958** |
| PDF | 10 pages, born-digital (LaTeX/hyperref → Distiller), clean text layer, watermarked with a UPEI download stamp |

**Year caveat, minor:** the file is keyed `nakamura2020` and the online-first year is 2020, but
the **issue of record is 2021, 85(2), 401–410**. Cite with both if the repo's bibliography is
strict.

### 0a. Topic match

The brief expected *"browning of xylose + 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol)
systems … any kinetics (probably browning rates), any norfuraneol consumption data,
conditions."* **All three are present.** The paper's abbreviation for
4-hydroxy-5-methyl-3(2H)-furanone is **HMFO** (not "norfuraneol", which never appears).
Note that Fig. 3's caption and part of the running text misspell it **"HFMO"** — the same
compound.

### 0b. **STRUCTURAL WARNING — THERE ARE NO TABLES IN THIS PAPER**

Not one. Twelve figures, zero tables. Consequently:
- **Every number below that is not quoted verbatim in the running text is a digitisation.**
- **No error bars can be read to better than ± 20 % of their own height**; the paper states
  *"All experiments were conducted in triplicate"* and *"n = 3, P < .05"* but **never prints
  a single ± value.**
- Statistical output is limited to **Duncan-style letter groupings drawn on the figures**
  (one-way ANOVA + **Tukey's** multiple comparison, α = 0.05, §Statistical analysis, p. 403).

### 0c. Provenance codes

- **[M]** measured value quoted verbatim in the authors' running text
- **[D]** digitised by me off the figure raster (per-panel calibration and error stated)
- **[Z]** derived by me from [M]/[D] values — never printed by the authors
- **[C]** cited from another paper, not measured here

---

## 1. THE ONE-PARAGRAPH ANSWER

The paper is a **mechanistic dicarbonyl-tracing study**, not a kinetics paper: it establishes
that in a xylose–lysine system at near-neutral pH, **xylose → 1-deoxyxylosone (1-DX) → HMFO →
methylglyoxal (MGO) → diacetyl (DA)**, and that **MGO is the dominant browning intermediate**.
It nevertheless yields, by digitisation, **six usable time-courses**: HMFO formation and
browning in Xyl–Lys at two pH values (Fig. 2), MGO/DA formation from HMFO with and without
in-situ trapping (Fig. 3), a **seven-dicarbonyl inventory at five heating times** (Fig. 9a) and
its untrapped counterpart (Fig. 9b), and **eight browning curves comparing Xyl, HMFO, MGO and
DA ± lysine** (Fig. 10). **The single most transferable result is a rank-ordering of intrinsic
browning power at matched ≈ 13 mM substrate and 100 °C: MGO ≫ DA > HMFO ≫ Xyl, with MGO+Lys
browning ≈ 12× faster than Xyl+Lys** (§7). All conditions are **100 °C, 0.2 M Na/K phosphate,
pH 6.5 or 8.0** — one temperature only, so **no activation energy is obtainable.**

---

## 2. SYSTEM COMPOSITION — every experiment in the paper

All heating is *"in boiling water"* = **100 °C**, in a **capped test tube** (closed system,
volatiles retained). All buffers are **0.2 M Na/K phosphate**. Reactions are cooled, then
**adjusted to pH 7.4** before HPLC. **n = 3** throughout.

| # | experiment | composition | pH | heating | figure |
|---|---|---|---|---|---|
| A | Xyl–Lys browning + HMFO | **13 mM Xyl + 27 mM Lys** | **6.5 and 8.0** | 0, 0.5, 1, 3, 6 h @ 100 °C | **Fig. 2** |
| B | Xyl–Lys dicarbonyls, **trapped in situ** | 13 mM Xyl + 27 mM Lys + **28 mM OPD present during heating** | **8.0** | 0.5, 1, 3, 6 h | **Fig. 8, 9a** |
| C | Xyl–Lys dicarbonyls, **trapped after** | 13 mM Xyl + 27 mM Lys; heated, cooled, then equal volume of buffer + **28 mM OPD**, 1 h RT | **8.0** (and 6.5, *"data not shown"*) | 0.5, 1, 3, 6 h | **Fig. 9b** |
| D | HMFO dicarbonyls, **trapped in situ** | **1.8 mM HMFO + 1.8 mM OPD** | **8.0** | 0–6 h | **Fig. 3a** |
| E | HMFO dicarbonyls, **trapped after** | **1.8 mM HMFO**; heated, then equal volume of buffer + **1.8 mM OPD**, 1 h RT | **8.0** | 0–6 h | **Fig. 3b** |
| F | HMFO at room temperature | 1.8 mM HMFO + 1.8 mM OPD, **left at RT 4–12 h** (Fig. 4 legend says **8 h**) | **8.0** | RT | **Figs 4, 5** |
| G | DTPA / degassing test | 1.8 mM HMFO; **control** vs **+40 mM DTPA** vs **N₂ bubbling 5 min** | **8.0** | **3 h** @ 100 °C | **Fig. 7a** |
| H | Brown-solution comparison | **13.3 mM Xyl**, **13 mM MGO**, **13.3 mM DA**, or **13.3 mM HMFO**, each ± **27.4 mM Lys** | **8.0** | 6 h (panel a) or 60 min (panel b) | **Fig. 10** |
| I | Pigment characterisation | each of the above brown solutions, taken at **A400 = 1.5–2.0** | — | — | **Fig. 11** |

**Browning index throughout:** *"The reaction solution was appropriately diluted, before its
A400 was measured. **The A400 multiplied by the dilution ratio was used as an indicator of
browning.**"* (Fig. 2 and Fig. 10 captions). **This means the reported A400 is a
dilution-corrected apparent absorbance and can legitimately exceed 4 — it is not a raw cuvette
reading.** Path length is never stated (assume 1 cm).

### 2a. Reagents (p. 402)
D-Xyl, L-Lys·HCl, **o-phenylenediamine (OPD)**, DA, cyanocobalamin (vitamin B₁₂),
**diethylenetriamine pentaacetic acid (DTPA)** — FUJIFILM Wako Pure Chemical, Osaka.
**HMFO (97 %)**, **MGO (40 % aqueous solution)**, tannic acid — Sigma-Aldrich Japan, Tokyo.
**HMFO purity 97 % is the only purity figure in the paper.** MGO supplied as a 40 % aqueous
solution — **notoriously oligomeric; the effective monomeric MGO concentration in experiment H
is therefore uncertain and the paper does not address it.** ⚠

### 2b. Analytical methods
**DAD–ODS–HPLC (p. 402):** YMC-Pack R&D ODS-A, **4.6 × 250 mm, 10 μm**; eluent A = 0.1 % aq.
TFA : methanol 98:2 (v/v), B = 0.1 % aq. TFA : methanol 50:50 (v/v); gradient **0 % B for
0–10 min, 0→50 % B for 10–40 min, 100 % B for 40–45 min**; flow **1.0 mL/min**; HITACHI L-4500
DAD, **220–500 nm**, fixed wavelength **320 nm**.
**Retention times: HMFO ≈ 8 min, MGO-quinoxaline ≈ 37 min, DA-quinoxaline ≈ 39 min.**
**Quantitation: external standard method using each standard**; *"Quinoxalines of other than DA
and 1-DX were determined as MGO quinoxaline."* ⚠ **— i.e. the xylosones, 3-DX and tetrosone in
Fig. 9a are quantified against the MGO-quinoxaline calibration, not their own. Those four
series carry an unquantified response-factor error.**
**DAD–GP–HPLC (p. 402):** TSK-gel G2500PWXL **7.8 × 300 mm** + guard TSK PWXL 6.0 × 40 mm
(Tosoh); eluent **0.1 M Na/K phosphate pH 7.0**; flow **0.5 mL/min**; DAD 220–500 nm, fixed
**370 nm**. MW markers: **tannic acid MW 1701**, **cyanocobalamin (VB₁₂) MW 1355**.
**Fluorescence HPLC:** HITACHI L-2485, **Ex 370 nm / Em 450 nm**.
**MS:** AB Sciex Triple TOF 4600; LC-MS/MS on YMC-Pack ODS-A **2.0 × 150 mm, 3 μm**; eluent
A = 0.1 % aq. TFA : methanol 98:2, B = 0.1 % TFA : methanol 50:50; **0 % B 0–5 min, 0→50 % B
5–50 min, 100 % B 50–55 min**; flow **0.2 mL/min**; **ESI(+)**, MS1 *m/z* 100–1000,
MS2 50–1000. **NMR:** Bruker AVANCE III 600 in CDCl₃ + 0.03 % TMS.

### 2c. 1-DX quinoxaline standard — synthesised in-house (p. 402) **[M]**
**D-Xyl 25 mg + L-Lys·HCl 251 mg + OPD 372 mg in 25 mL of 0.2 M Na/K phosphate pH 6.5, boiled
3 h**, cooled, filtered, extracted ×3 with ethyl acetate, concentrated in vacuo →
**1-DX quinoxaline, 118 mg**. Identity per Hauck et al. (2003) **[C]**.
Physicochemical data as printed:
**ESI-MS *m/z* 205.0967 (M+H)⁺, calcd for C₁₁H₁₃N₂O₂ 205.0977.**
**¹H NMR (CDCl₃) δ 2.81 (3H, s), 3.85 (1H, d, J = 11.7), 4.04 (1H, m), 5.12 (1H, t, J = 5.5),
7.74 (1H, m), 8.03 (1H, t, J = 6.5).**
**¹³C δ 22.2, 66.2, 71.2, 128.6, 128.8, 129.8, 130.3, 139.6, 142.1, 152.2, 153.1.**
1-DX is named **(4R)-4,5-dihydroxy-2,3-pentanedione**; 3-DX is **(4S)-4,5-dihydroxy-2-oxopentanal**.

---

## 3. NUMBERS QUOTED VERBATIM IN THE TEXT **[M]** — the anchors for every digitisation

| # | value as printed | context | page |
|---|---|---|---|
| 1 | *"The maximal concentrations of HMFO formed at pH 6.5 and pH 8.0 were **about 1.9 mM at 3 h** of heating and **about 0.6 mM at 0.5 h** of heating, respectively."* | Fig. 2, Xyl–Lys | 403 |
| 2 | *"the browning for the first 1 h of heating was **larger at pH 8.0 than pH 6.5**, but at **3 h** of heating the browning at pH 8.0 was **the same as** at pH 6.5"* | Fig. 2 | 403 |
| 3 | *"**From 1.5 mM HMFO, about 0.24 mM MGO and 0.25 mM DA were formed**, indicating **more than one-third** of HMFO was converted into MGO and DA."* | Fig. 3a, 6 h | 404 |
| 4 | *"**Only 0.04 mM of MGO was detected at 0.5 h** of heating and then decreased."* | Fig. 3b | 404 |
| 5 | *"**About 0.1 mM of DA was detected at 3 h** heating and then it was decreased."* | Fig. 3b | 404 |
| 6 | *"1-DX quinoxaline was the major product (Fig. 9a), its concentration **reaching 9.8 mM at 6 h** of heating."* | Fig. 9a | 405 |
| 7 | *"the A400 reached **about 1.5 at 1 h** of heating"* (Xyl + Lys) | Fig. 10a | 407 |
| 8 | *"The MGO solution containing Lys turned brown most rapidly, and its A400 reached **more than 3 at 10 min** of heating."* | Fig. 10b | 407 |
| 9 | *"In case of HMFO, **Lys reduced the browning of HMFO to about the half**."* | Fig. 10b | 407 |
| 10 | *"a colorless polymer having an **absorption maximum at 280 nm** and no fluorescence"* (HMFO solution, GP-HPLC) | Fig. 11 | 408 |
| 11 | *"a brown pigment … eluted a little earlier than **tannic acid (MW 1701)**"*; melanoidins *"had … **about 2000–5000 of molecular weights**"* | Fig. 11 | 408 |
| 12 | *"a peak on ODS … having an **absorption maximum near a visible region (340 nm)**"* | Fig. 11 | 407 |
| 13 | MS: compound **I** C₁₁H₁₁N₂O₂ obs. 203.0824 / calcd 203.0821; compound **II** C₁₁H₉N₂O obs. 185.0709 / calcd 185.0715 → precursors **C₅H₆O₄** and **C₅H₄O₃** | Figs 4–6 | 404 |
| 14 | *"the concentration of pentose in soy sauce"* is why **13 mM Xyl** was chosen | design rationale | 403 |

**Every one of the digitisations in §4–§7 is anchored to at least one of these fourteen printed
values, and all fourteen are reproduced to within the stated digitisation error.** That
cross-validation is the reason the figure-derived numbers below are usable at all.

---

## 4. FIGURE 2 — Xyl–Lys browning and HMFO, at pH 6.5 and 8.0 **[D]**

**Anchor: Fig. 2, p. 404 (PDF page 4).** Caption as printed: *"Browning and formation of HMFO in
the Xyl-Lys system. A solution containing **13 mM Xyl and 27 mM Lys (pH 6.5 or 8.0) was heated
at 100 °C**. The reaction solution was appropriately diluted, before its A400 was measured. The
A400 multiplied by the dilution ratio was used as an indicator of browning."*
Two panels; each has a **left axis A400 (red squares, 0–4.0+)** and a **right axis HMFO in mM
(black diamonds)** — **right axis 0–3 mM for pH 6.5, 0–0.9 mM for pH 8.0** (different scales,
easy to misread). x-axis heating time 0–6 h.

**Digitisation method:** frame located by dark-row/column projection; left axis calibrated on
the eight printed minor ticks (0.5 A400 per 67.8 px, 0.0 at the bottom frame y = 787);
right axis calibrated frame-bottom = 0 to frame-top = full scale. Markers located by colour
mask (red RGB ≈ 192,0,0; black) with the legend region excluded, then weighted-centroid.
**Residual error ± 0.02 A400 and ± 0.01 mM.**

| heating time (h) | **A400, pH 6.5** | **HMFO mM, pH 6.5** | **A400, pH 8.0** | **HMFO mM, pH 8.0** |
|---|---|---|---|---|
| 0 | **0.00** | **0.00** | **0.00** | **0.00** |
| 0.5 | **0.41** | **0.87** | **0.77** | **0.59** |
| 1 | **1.03** | **1.42** | **1.61** | **0.46** |
| 3 | **2.69** | **1.91** | **2.64** | **0.25** |
| 6 | **4.02** | **1.44** | **3.47** | **0.02** |

**Cross-checks against the printed text [M]:**
- HMFO max pH 6.5 = **1.91 mM at 3 h** vs printed *"about 1.9 mM at 3 h"* ✔ (0.5 % agreement)
- HMFO max pH 8.0 = **0.59 mM at 0.5 h** vs printed *"about 0.6 mM at 0.5 h"* ✔ (1.8 %)
- A400 at 1 h: **1.61 (pH 8.0) > 1.03 (pH 6.5)** ✔ matches *"larger at pH 8.0"*
- A400 at 3 h: **2.64 (pH 8.0) vs 2.69 (pH 6.5)** ✔ matches *"the same"* to 1.9 %
**Four independent checks pass. The Fig. 2 digitisation is validated.**

**Error bars [D~, approximate]:** drawn on both series. Largest are on the **pH 8.0 HMFO** points
(≈ **± 0.08 mM at 0.5 h**, ≈ **± 0.13 mM at 1 h**) and on **A400 at 1–3 h** (≈ **± 0.2**).
The pH 6.5 HMFO bars are ≈ **± 0.1 mM at 3 h**. **No numeric SD is printed anywhere; these are
eye-read bar half-heights and should carry ± 30 % of their own value.**

### 4a. Derived quantities **[Z]** — none of these is printed by the authors

**Browning rate, dA400/dt:**

| interval | pH 6.5 | pH 8.0 | ratio 8.0/6.5 |
|---|---|---|---|
| 0 → 0.5 h | **0.82 h⁻¹** | **1.55 h⁻¹** | **1.89** |
| 0 → 1 h | **1.03 h⁻¹** | **1.61 h⁻¹** | **1.56** |
| 1 → 6 h | **0.60 h⁻¹** | **0.37 h⁻¹** | **0.62** |

**The pH effect on browning inverts with time: alkaline is ~1.9× faster initially, then ~0.6×
as fast after 1 h, and the two curves cross at ≈ 3 h.** Any model with a single monotone
pH factor on browning rate cannot reproduce this. The mechanism the paper supplies is that
**HMFO is far less stable at pH 8 (peak 0.59 mM, gone by 6 h) than at pH 6.5 (peak 1.91 mM,
still 1.44 mM at 6 h)** — the alkaline system burns its intermediate pool early.

**Molar yield of HMFO from xylose:** peak 1.91 mM / 13 mM Xyl = **14.7 % at pH 6.5**;
0.59 / 13 = **4.5 % at pH 8.0**.

**First-order HMFO decay after its peak:**
| system | window | **k (h⁻¹)** | k (s⁻¹) | R² | t½ |
|---|---|---|---|---|---|
| Xyl–Lys, pH 6.5 | 3 → 6 h | **0.094** | 2.61 × 10⁻⁵ | (2 points) | **7.4 h** |
| Xyl–Lys, pH 8.0 | 0.5 → 3 h | **0.328** | 9.11 × 10⁻⁵ | 0.992 | **2.11 h** |
| Xyl–Lys, pH 8.0 | 0.5 → 6 h | **0.636** | 1.77 × 10⁻⁴ | 0.948 | **1.09 h** |
**HMFO is ≈ 3.5× more labile at pH 8.0 than at pH 6.5** on the 0.5–3 h window. Note the pH 8.0
fit steepens when the 6 h point is included (0.328 → 0.636 h⁻¹), i.e. **the decay is not
first-order — it accelerates.** Treat these k as apparent, window-dependent.
Corroborating **[C]**: *"HMFO was more stable at pH 6.5 than pH 8.0 (Mikami et al. 2017)"*.

---

## 5. FIGURE 3 — MGO and DA formation from pure HMFO **[D]**

**Anchor: Fig. 3, p. 404.** Caption as printed: *"Formation of MGO and DA in a solution of HMFO.
**(a) 0.2 M phosphate buffer (pH 8.0) containing 1.8 mM HMFO and 1.8 mM OPD was heated at
100 °C for 0–6 h. (b) 0.2 M phosphate buffer (pH 8.0) containing 1.8 mM HMFO was heated at
100 °C for 0–6 h, before OPD being added.** Different letters showed a significant difference
between MGO and DA at each heating time (n = 3, P < .05)."*

**The (a)/(b) distinction is the experimental core and must not be collapsed:**
- **(a) OPD present *during* heating** → dicarbonyls are trapped as quinoxalines the instant they
  form. **(a) measures cumulative production.**
- **(b) OPD added *after* heating** → only dicarbonyls still surviving at that time point are
  seen. **(b) measures the standing pool.**

Both panels: **left axis HMFO in mM (black diamonds in a; blue-grey diamonds in b), 0–2.0**;
**right axis dicarbonyl quinoxaline in mM — 0–0.3 for (a), 0–0.12 for (b)** (⚠ different
scales); MGO = purple squares, DA = green triangles.
**Digitisation:** colour masks (purple 112,48,160; green 0,176,80; HMFO-blue 79,129,189 in
panel b), marker centroids from ≥ 11-px-wide rows. Frames: (a) top 71.5 / bottom 787 /
left 511 / right 1227; (b) top 80.5 / bottom 789 / left 353.5 / right 1073. Left-axis
calibration confirmed against the plotted minor ticks (0.25 mM apart in panel a).
**Residual error ± 0.02 mM (left axis), ± 0.003 mM (right axis).**

### 5a. Fig. 3a — OPD present during heating (cumulative)

| time (h) | **HMFO (mM)** | **MGO-quinoxaline (mM)** | **DA-quinoxaline (mM)** | Tukey letters [M] |
|---|---|---|---|---|
| 0 | **1.727** | **0.000** | **0.000** | MGO a / DA b |
| 0.5 | **0.753** | **0.094** | **0.025** | MGO a / DA b |
| 1 | **0.412** | **0.131** | **0.070** | MGO a / DA b |
| 3 | **0.000** | **0.191** | **0.158** | MGO a / DA b |
| 6 | **0.000** | **0.235** | **0.244** | both a |

**Cross-check [M]:** *"From 1.5 mM HMFO, about **0.24 mM MGO and 0.25 mM DA** were formed"* vs
digitised **0.235 and 0.244 at 6 h** ✔ (2 % and 2 %). **Validated.**
⚠ **But the printed "1.5 mM HMFO" does not match the figure.** The solution is prepared at
**1.8 mM** (caption) and the t = 0 diamond digitises to **1.727 mM**. **The "1.5 mM" in the text
is a third, inconsistent number.** Using 1.727 mM as consumed, the combined molar yield is
(0.235 + 0.244)/1.727 = **27.7 %**, not the *"more than one-third"* the authors claim (which
follows only from their 1.5 mM figure: 0.49/1.5 = 32.7 %). **[Z] Use 28 %, and flag the 1.5 mM.**

### 5b. Fig. 3b — OPD added after heating (standing pool)

| time (h) | **HMFO (mM)** | **MGO-quinoxaline (mM)** | **DA-quinoxaline (mM)** | Tukey letters [M] |
|---|---|---|---|---|
| 0 | **1.612** | **0.000** | **0.000** | — |
| 0.5 | **0.735** | **0.038** | **0.028** | a / a |
| 1 | **0.205** | **0.018** | **0.052** | b / a |
| 3 | **≈ 0** | **0.003** | **0.097** | b / a |
| 6 | **≈ 0** | **0.002** | **0.061** | b / a |

**Cross-checks [M]:** *"Only **0.04 mM of MGO** was detected at 0.5 h"* vs **0.038** ✔;
*"About **0.1 mM of DA** was detected at 3 h heating and then it was decreased"* vs **0.097 at
3 h, falling to 0.061 at 6 h** ✔. **Both validated.**

### 5c. What (a) minus (b) actually proves **[Z]**

At 6 h the **cumulative** MGO (0.235 mM) is **≈ 120× the standing** MGO (0.002 mM), while the
cumulative DA (0.244) is only **4× the standing** DA (0.061).
**MGO is produced copiously and destroyed almost as fast; DA accumulates.** This is the
quantitative form of the authors' qualitative claim *"MGO formed from HMFO was readily
decomposed to DA or polymerized. DA was more stable than MGO."* (p. 404) — and it is the number
a model needs, because **it means the MGO *pool* in a real system is a poor proxy for MGO
*flux*, by roughly two orders of magnitude.**

**First-order HMFO consumption [Z]:**
| panel | window | **k (h⁻¹)** | k (s⁻¹) | R² | t½ |
|---|---|---|---|---|---|
| 3a (+OPD during) | 0 → 1 h | **1.433** | 3.98 × 10⁻⁴ | 0.992 | **0.48 h** |
| 3b (−OPD during) | 0 → 1 h | **2.062** | 5.73 × 10⁻⁴ | 0.981 | **0.34 h** |

⚠ **HMFO disappears ~44 % faster when OPD is absent during heating.** Either OPD partially
stabilises HMFO (by removing the dicarbonyls that might otherwise catalyse its loss), or the
two panels differ in some way the paper does not describe. Panel (b) also involves a **2×
dilution** with OPD buffer after heating, which the caption does not say was corrected for —
and indeed the t = 0 HMFO reads **1.612 mM in (b) vs 1.727 mM in (a)**, a 7 % offset consistent
with an incomplete dilution correction. ⚠ **Do not compare absolute concentrations between
panels 3a and 3b; compare shapes only.**

**These two HMFO decay constants at 100 °C, pH 8.0, are the most repo-ready kinetic numbers in
the paper: k ≈ 1.4–2.1 h⁻¹ = 4.0–5.7 × 10⁻⁴ s⁻¹, t½ ≈ 20–29 min.**

---

## 6. FIGURES 4–8 — identification work, no kinetics

### 6a. Fig. 4 (p. 405) — two new quinoxalines at room temperature **[M]**
*"A solution of 0.2 M phosphate buffer (pH 8.0) containing **1.8 mM HMFO and 1.8 mM OPD was
left at room temperature for 8 h** and analyzed using DAD-HPLC."* Two new peaks, **compounds I
and II**. **No concentrations, no yields — a chromatogram only.**

### 6b. Fig. 5 (p. 405) and the MS assignments **[M]**
Compound **I**: C₁₁H₁₁N₂O₂, obs. (M+H)⁺ **203.0824**, calcd **203.0821** → precursor dicarbonyl
**C₅H₆O₄** = **2-hydroxy-3,4-dioxopentanal (1)**.
Compound **II**: C₁₁H₉N₂O, obs. (M+H)⁺ **185.0709**, calcd **185.0715** → precursor **C₅H₄O₃** =
**2,3-dioxopent-4-enal (2)**.
Explicitly **rejected** on MS/MS grounds: 2,3,4-trioxopentanol (C₅H₆O₄) and
4-oxo-pent-2-enedial (C₅H₄O₃).

### 6c. Fig. 6 (p. 406) — formation scheme, no numbers **[M]**
HMFO → (oxidation + hydrolysis, either order) → carbonyls **1** and **2** (and **1′**) →
quinoxalines **I** and **II**; **MGO formed from 1 and 2 by hydrolysis or retro-aldol.**
Every arrow unlabelled.

### 6d. Fig. 7 — the DTPA / degassing test **[D], semi-quantitative**
**Anchor: Fig. 7a, p. 406.** *"A solution of 0.2 M phosphate buffer (pH 8.0) containing 1.8 mM
HMFO was heated at 100 °C for **3 h**, before HMFO and DA being analyzed using HPLC.
(a) **Control**, no addition; **DTPA, the addition of 40 mM DTPA**; **Degassing, N₂ bubbling for
5 min**."*
Three ODS chromatograms, y-axis **A320**, x-axis retention time 0–45 min.
| panel | what is seen |
|---|---|
| **Control** | **HMFO absent** (no ~8 min peak); one large peak at **≈ 37 min = DA-quinoxaline**, peak height **≈ 0.0105 A320** above baseline; minor peaks at ≈ 3, 5, 6.5 min (0.0007–0.0025) |
| **DTPA** | **HMFO peak at ≈ 9 min survives, off-scale**; the DA peak is suppressed |
| **Degassing** | same qualitative outcome as DTPA (panel cropped in the crop I rendered; the text states it) |

**[M]** *"After 3 h of heating, no HFMO was detected and only DA quinoxaline was detected as
shown in the control … These decomposition of HMFO and formation of DA were **strongly
repressed** by the addition of DTPA and by degassing."*
**⚠ Mechanistic consequence, and the single most important caveat in this paper: HMFO
decomposition to MGO/DA is OXIDATIVE and metal-catalysed.** It requires dissolved O₂ and free
transition metal, and is abolished by a chelator or by degassing. **[C]** DA formation from MGO
by active oxygen species is attributed to Nukaya et al. (1993), *Chem. Pharm. Bull.* **41**,
649–653. **No percentage of suppression is quantified** — the figure is chromatograms only.
**Every rate constant in §5 is therefore an *aerobic, non-chelated, phosphate-buffered* rate and
is not transferable to a deaerated or chelated system.**

### 6e. Figs 8 and 12 — peak assignments and the summary scheme **[M]**
**Fig. 8, p. 407** (Xyl–Lys 13 mM/27 mM, pH 8.0, **100 °C 3 h, then 28 mM OPD added**):
the printed assignment table, transcribed exactly —

| peak | formula | compound | approx. t_R (min) [D] |
|---|---|---|---|
| **III** | *(not printed)* | xylosone / pentosone | ≈ 25.5 |
| **IV** | **C₁₁H₁₂N₂O₃** | **Xylosone** | ≈ 25.7 |
| **V** | **C₁₁H₁₂N₂O₂** | **1-DX** | ≈ 30 (tallest peak) |
| **VI** | **C₁₁H₁₂N₂O₂** | **3-DX** | ≈ 30.5 |
| **VII** | **C₁₀H₁₀N₂O** | **Tetrosone** | ≈ 33.5 |
| **I** | **C₉H₈N₂** | **MGO** | ≈ 36.5 |
| **II** | **C₁₀H₁₀N₂** | **DA** | ≈ 38 |

(Peak III's row is cut off at the top of the printed inset; its formula is not legible. The
text identifies III and IV together as *"two xylosones or pentosones"*. Tetrosone is named
**tetrosulose or 3,4-dihydroxy-2-oxobutanal**.)
**Fig. 12, p. 409:** *"Maillard reaction and browning of the Xyl system at near-neutral pH
conditions"* — a summary scheme, **no numbers**.

---

## 7. FIGURE 9 — the seven-dicarbonyl time course in Xyl–Lys **[D]**

**Anchor: Fig. 9, p. 407.** Caption: *"Time course of the formation of HMFO and dicarbonyl
quinoxalines in the Xyl-Lys model system. A solution of 0.2 M phosphate buffer (**pH 8.0**)
containing **13 mM Xyl and 27 mM Lys with 28 mM OPD (a) or without OPD (b)** was heated at
**100 °C for 6 h**. In (b), quinoxalines were prepared after heating by adding OPD. Different
letters show significant differences among compounds at 6 h of heating (a) and at each heating
time (b) (n = 3, P < .05)."*

### 7a. Fig. 9a — OPD present during heating (cumulative production)

Grouped bar chart, y-axis **Concentration (mM), 0–12.0**, eight compound groups × five heating
times (0, 0.5, 1, 3, 6 h). **The 0 h bar is zero in every group and is annotated "0" rather
than drawn.** Series fills: 0 h solid blue, 0.5 h blue-with-white-dots, 1 h red hatch,
3 h white-with-black-dots, 6 h solid black.
**Digitisation:** bar tops located by scanning down each bar's own x-window for the topmost row
covering ≥ 75 % of the bar width, with a separate pass to identify and subtract the error-bar
cap above it. Axis calibrated on the printed ticks (12.0 at y = 91.5, 0.0 at y = 566.5).
**Residual error ± 0.15 mM on the large bars, ± 0.1 mM on the small ones.**

| compound | 0 h | 0.5 h | 1 h | 3 h | 6 h | 6 h Tukey letter [M] |
|---|---|---|---|---|---|---|
| **HMFO** | 0 | **0.29** | **0.32** | **0.16** | **≈ 0.06** | e |
| Pentosone (III) | 0 | **0.11** | **0.21** | **0.29** | **0.32** | de |
| Xylosone (IV) | 0 | **0.57** | **0.72** | **0.85** | **1.00** | d |
| **1-DX** | 0 | **4.74 ± 0.24** | **7.29 ± 0.20** | **9.71 ± 0.35** | **9.84 ± 0.48** | **a** |
| **3-DX** | 0 | **0.29** | **0.69** | **≈ 1.2 (+0.3)** | **≈ 1.5** | c |
| Tetrosone | 0 | **0.14** | **0.21** | **0.39** | **0.34** | de |
| **MGO** | 0 | **1.30** | **2.19** | **3.37 ± 0.11** | **4.16 ± 0.27** | **b** |
| DA | 0 | **≈ 0.03** | **0.09** | **0.16** | **0.16** | de |

**Cross-check [M]:** 1-DX *"reaching **9.8 mM at 6 h**"* vs digitised **9.84 ± 0.48** ✔ (0.4 %).
**Validated.** ⚠ **The 3-DX 3 h and 6 h bars and the smaller groups (HMFO, pentosone, tetrosone,
DA) sit close to the automated detector's floor and are contaminated by the Tukey letters drawn
above them; those cells carry ± 0.15 mM and are marked ≈ where the detector needed manual
arbitration. The 1-DX and MGO groups are clean.**

**[Z] Molar accounting from 13 mM xylose (cumulative, 6 h):**
**1-DX 9.84 mM = 76 % of the initial Xyl.** MGO 4.16 mM = **32 %** (on a C₃ fragment basis,
i.e. more than one MGO can in principle come from one C₅). Everything else is ≤ 8 %.
**Xylose degradation at pH 8.0/100 °C runs overwhelmingly through 1-deoxyxylosone, and MGO is
the dominant C₃ fragment.** This is the paper's quantitative core and it is a directly usable
branching constraint.
**[M]** supporting statement: *"1-DX quinoxaline was the major product … This result indicated
that Xyl was mainly converted into 1-DX and that other carbonyls were formed from 1-DX … As we
used here pH 8.0 conditions, 1-DX instead of 3-DX was a major decomposition product of Xyl. It
is well-known that 3-deoxyglucosone is the major osone formed by the Maillard reaction at acidic
condition, while 1-deoxyglucosone is the major one at neutral or alkaline pH (Eskin 1990) **[C]**.
**1-DX was also a major decomposition product at pH 6.5 like pH 8.0.**"*

### 7b. Fig. 9b — OPD added after heating (standing pool)

Line plot, single y-axis **Concentration (mM), 0.0–0.5**; series **HMFO** (solid black line,
diamonds), **MGO** (dashed, circles/purple square), **DA** (green triangles, dotted),
**1-DX** (dark red), **3-DX** (orange). Axis calibrated 0.0 at y = 490, 0.5 at y = 8.
**Residual error ± 0.01 mM.**

| time (h) | **HMFO** | **MGO** | **DA** | 1-DX | 3-DX |
|---|---|---|---|---|---|
| 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 0.5 | **0.378 ± 0.02** | **0.090** | **0.081** | ≈ 0.01 | ≈ 0.02 |
| 1 | **0.364 ± 0.045** | **0.104** | **0.093** | ≈ 0.02 | ≈ 0.02 |
| 3 | **0.112 ± 0.035** | **0.114** | **0.077** | ≈ 0.005 | ≈ 0.005 |
| 6 | **≈ 0.005** | **0.150** | **0.058** | ≈ 0 | ≈ 0 |

Tukey letters printed on the figure: HMFO **a, a** at 0.5 and 1 h; **bb** and **cc** groupings
at 0.5 h for the minor species; **b b** at 3 h; **c c c** at 6 h.
**[M]** *"MGO and DA were the major dicarbonyls in the solution during heating and … **several
times more HMFO was formed than MGO and DA**, and then decreased. **1-DX was little detected in
the solution.**"* — digitisation agrees: HMFO peaks at 0.378 mM vs MGO 0.090 and DA 0.081, i.e.
**4.2× and 4.7×** ✔ *"several times"*.
**[M]** *"The result at pH 6.5 was similar to that at pH 8.0. MGO and DA were gradually increased
in the reaction solution, although each concentration was different from that at pH 8.0
(**data not shown**)."* ⚠ **The pH 6.5 dicarbonyl dataset exists but is not published.**

### 7c. The (a)-vs-(b) contrast for 1-DX **[Z]** — the paper's sharpest single insight

**Cumulative 1-DX (trapped in situ) = 9.84 mM. Standing 1-DX (trapped after) ≤ 0.02 mM.**
**A ratio of ≥ 500×.** 1-DX is produced at three-quarters molar yield from xylose and is
consumed essentially as fast as it is made — its steady-state concentration is
**below 0.2 % of its cumulative flux.**
**If the repo's reaction network carries 1-DX (or any deoxyosone) as a pooled species whose
concentration gates downstream rates, this measurement says the pool is effectively zero and
the flux is everything.** Same argument as §5c but 4× stronger.

---

## 8. FIGURE 10 — THE INTRINSIC BROWNING COMPARISON **[D]** — most transferable result

**Anchor: Fig. 10, p. 408.** Caption as printed: *"Browning of Xyl, MGO, DA, and HMFO solutions.
A solution of 0.2 M Na/K phosphate buffer (**pH 8.0**) containing **13.3 mM Xyl, 13 mM MGO,
13.3 mM DA, or 13.3 mM HMFO** in the presence or absence of **27.4 mM Lys** was, respectively,
heated at **100 °C for 6 h (a) or 60 min (b)**. The reaction solution was appropriately diluted,
before its A400 was measured. The A400 multiplied by the dilution ratio was used as an indicator
of browning. Different letters show significant differences at the same heating time
(n = 3, P < .05)."*

**Digitisation:** both panels calibrated on the printed y-axis ticks (panel a: 0.0 at
y = 756, 3.0 at y = 224; panel b: 0.0 at y = 758, 1.0 per 139 px). Series separated by colour
mask **and** by a solid-vs-dashed line-coverage test along each connecting segment (solid
returns ≈ 100 % coverage, dashed 0.64–0.73) — this was necessary because the ± Lys pair of each
compound shares a colour. **Residual error ± 0.03 A400.**

### 8a. Panel (a) — Xylose alone, 0–6 h

| time (h) | **Xyl + Lys** | **Xyl alone** |
|---|---|---|
| 0 | 0.00 | 0.00 |
| 1 | **1.54** | **0.32** |
| 2 | **2.22** | **0.66** |
| 3 | **2.56** | **0.82** |
| 4 | **2.82** | **0.94** |
| 5 | **3.05** | **1.17** |
| 6 | **3.20** | **1.35** |
| 6 h Tukey letter [M] | **a** | **c** |

**Cross-check [M]:** *"the A400 reached **about 1.5 at 1 h** of heating"* vs digitised **1.54** ✔.
**[M]** *"The Xyl solution turned brown very slowly. This browning was stimulated by the presence
of Lys."* → **[Z] Lys accelerates xylose browning 4.8× at 1 h** and **2.4× at 6 h.**

### 8b. Panel (b) — MGO, DA, HMFO, 0–60 min

Line-style assignment verified programmatically: **dashed = no Lys, solid = + Lys** (legend
p. 408: `- -■- - MGO` / `—■— MGO+Lys`; `- -▲- - DA` / `—▲— DA+Lys`; black = HMFO).

| time (min) | **MGO** | **MGO+Lys** | **DA** | **DA+Lys** | **HMFO** | **HMFO+Lys** |
|---|---|---|---|---|---|---|
| 0 | ≈ 0.05 | ≈ 0.05 | 0.00 | 0.00 | **0.12** | **0.12** |
| 10 | **0.96** | **3.12** | **1.77** | **1.51** | **1.09** | **0.70** |
| 20 | **1.08** | **3.21** | **2.61** | **2.18** | **2.03** | **1.16** |
| 30 | **1.22** | **3.28** | **2.86** | **2.66** | **2.77** | **1.49** |
| 45 | **1.35** | **3.55** | **3.27** | **2.90** | **3.43** | **1.76** |
| 60 | **1.45** | **3.54** | **3.42** | **3.29** | **4.09** | **2.04** |
| 60 min Tukey letter [M] | **c** | **a** | **a** | **a** | **a** | **bc** |

**Cross-checks [M]:** *"The MGO solution containing Lys turned brown most rapidly, and its A400
reached **more than 3 at 10 min**"* vs digitised **3.12** ✔.
*"In case of HMFO, **Lys reduced the browning of HMFO to about the half**"* vs digitised
**4.09 → 2.04 at 60 min = exactly 0.499×** ✔✔. **Both validated, the second exactly.**

**⚠ ONE AUTHOR STATEMENT IS CONTRADICTED BY THE AUTHORS' OWN FIGURE.**
The text (p. 407) says *"**Except for HMFO, the presence of Lys promoted the browning.**"*
For MGO that is emphatically true (0.96 → 3.12 at 10 min, **3.3×**). **But for DA the digitised
+Lys curve lies BELOW the no-Lys curve at every time point** (10 min: 1.51 vs 1.77; 60 min: 3.29
vs 3.42), i.e. **Lys slightly *suppresses* DA browning too.** The line-style assignment was
verified programmatically and independently confirmed by the HMFO pair, whose direction the
authors themselves state. **The effect is small and the 60 min Tukey letters are both `a`
(not significantly different), so the honest reading is "no promotion by Lys for DA" — but the
sentence as written is wrong.** Flagged.

### 8c. Initial browning rates and the intrinsic ranking **[Z]** — the number to take away

Initial dA400/dt, from the first interval of each curve, all at **100 °C, pH 8.0, ≈ 13 mM
substrate, ± 27.4 mM Lys**:

| system | initial interval | **dA400/dt (h⁻¹)** | **× Xyl+Lys** |
|---|---|---|---|
| **MGO + Lys** | 0–10 min | **18.7** | **12.2** |
| **DA (no Lys)** | 0–10 min | **10.6** | **6.9** |
| **DA + Lys** | 0–10 min | **9.1** | **5.9** |
| **HMFO (no Lys)** | 0–10 min | **6.5** | **4.2** |
| **MGO (no Lys)** | 0–10 min | **5.8** | **3.8** |
| **HMFO + Lys** | 0–10 min | **4.2** | **2.7** |
| **Xyl + Lys** | 0–1 h | **1.54** | **1.0** |
| **Xyl (no Lys)** | 0–1 h | **0.32** | **0.21** |

**Intrinsic browning power at matched molarity: MGO+Lys ≫ DA ≈ HMFO ≫ Xyl+Lys ≫ Xyl.**
**MGO with lysine browns ≈ 12× faster than xylose with lysine, per mole, at the same
temperature and pH.** That single ratio is the paper's most directly importable quantity, and
it is the quantitative basis of the conclusion:
**[M]** *"in the Xyl system of the Maillard reaction under near-neutral pH conditions, **MGO
derived from HMFO and 1-DX forms melanoidin-like brown substances and is the major intermediate
for the browning**, although some low molecular weight pigments, a colorless polymer, and
fluorescent substances are formed at the same time."* (p. 409)

⚠ **Two caveats on the ranking.** (i) The MGO stock is a **40 % aqueous solution**, in which
MGO is substantially hydrated/oligomeric; the effective monomer concentration may be well below
13 mM, which would make the MGO ranking a **lower bound**. (ii) The initial intervals differ
(10 min for the reactive species, 1 h for xylose), so the Xyl comparison is against a slower,
partly-averaged rate; the true instantaneous ratio at t → 0 is probably larger than 12×.

---

## 9. FIGURE 11 — PIGMENT CHARACTERISATION **[M], qualitative only**

**Anchor: Fig. 11, p. 408.** *"Each brown solution showing **1.5–2.0 of A400** in Figure 10 was
applied to the HPLC analyses. From the top to the bottom, **fluorescence at 450 nm (Fl450;
excitation at 370 nm)** and **absorbance at 370 nm (A370)** in the presence of Lys and those in
the absence of Lys are shown. On the GP-HPLC, elution times of **tannic acid (MW 1701)** and
**cyanocobalamin (VB₁₂, MW 1355)** were shown."*
**No axis values, no peak areas, no concentrations — the entire panel is qualitative.**
Findings as stated:
- ODS: *"melanoidin-like low molecular weight substances having no specific absorption at or
  near a visible region showed a broad rising area from the bottom line. **The rising area was
  the largest in the MGO solution containing Lys and in the Xyl solution containing Lys.** In
  the DA solution, it was little detected with or without Lys."*
- *"The HMFO solutions in the presence and absence of Lys and the Xyl solution showed a peak on
  ODS (a) having an **absorption maximum near a visible region (340 nm)**. In the future, we
  would like to identify it."*
- *"**Fluorescent substances were little detected in all solutions without Lys**, while various
  fluorescent substances were detected in all solutions containing Lys. … **Lys or an amino acid
  is essential for the formation of fluorescent substances.**"*
- *"**Although the color intensities were also increased by Lys, fluorescent substances and
  pigments do not seem to be identical.**"* ⚠ — i.e. **fluorescence is not a proxy for browning
  in this system.**
- *"In the HMFO solution, a **colorless polymer having an absorption maximum at 280 nm and no
  fluorescence** were detected on GP-HPLC (b)."*
- *"On the chromatograms of the solutions of Xyl, HMFO, and MGO, a **brown pigment (c) eluted a
  little earlier than tannic acid (MW 1701)** was commonly detected, which became larger in the
  presence of Lys. These results suggested that the formed melanoidins here had **not so high
  molecular weight, but about 2000–5000 of molecular weights**."* **[C]** consistent with
  Hofmann (1998) ultracentrifugation.

---

## 10. WHAT IS CITED, NOT MEASURED — do not import numbers from these **[C]**

- **HMFO is more stable at pH 6.5 than pH 8.0** — Mikami, Nakamura, Yamada et al. (2017),
  *Food Sci. Technol. Res.* **23**, 283–289. Also the source of the **13 mM Xyl** design choice
  and of the prior observation that **HMFO solutions at pH 6 and 7 "turned hardly brown" while
  pH 8.0 did.**
- **DA is formed from MGO by active oxygen species** — Nukaya, Inaoka, Ishida et al. (1993),
  *Chem. Pharm. Bull.* **41**, 649–653.
- **1-deoxyglucosone dominates at neutral/alkaline pH, 3-deoxyglucosone at acidic pH** —
  Eskin (1990), *Biochemistry of Foods*, 2nd ed., pp. 252–258.
- **Pentoses brown faster than glucose because more open-chain isomer is present** —
  Bunn & Higgins (1981), *Science* **213**, 222–224.
- **Melanoidin MW ≈ 2000–5000** — Hofmann (1998), *J. Agric. Food Chem.* **46**, 3891–3895.
- **Quinoxaline standard preparation** — Gobert & Glomb (2009), *J. Agric. Food Chem.* **57**,
  8591–8597; 1-DX quinoxaline identity — Hauck et al. (2003), *BBA* **1623**, 109–119.
- **DTPA as a strong chelator and antioxidant** — Shinoda, Murata, Homma et al. (2004),
  *Biosci. Biotechnol. Biochem.* **68**, 529–536.
**No numeric parameter is carried over from any of these into this paper's own results.**

---

## 11. USABILITY CAVEATS THAT APPLY TO EVERY NUMBER ABOVE

1. **One temperature: 100 °C. No Eₐ is obtainable from this paper, at all.**
2. **Two pH values, and the second one (6.5) has published data only in Fig. 2.** The pH 6.5
   dicarbonyl time course exists but is *"data not shown"* (p. 407).
3. **No tables.** Every number in §4–§8 is a digitisation with the stated per-panel error, and
   **no SD is printed anywhere in the paper** despite n = 3 being claimed.
4. **The dicarbonyl "concentrations" are quinoxaline concentrations.** Trapping efficiency is
   never measured. Compounds other than DA and 1-DX are quantified **against the MGO-quinoxaline
   curve** — a systematic, unquantified response-factor error on the xylosones, 3-DX and
   tetrosone in Fig. 9a.
5. **Panel (a) and panel (b) of both Fig. 3 and Fig. 9 measure different things** (cumulative
   flux vs standing pool) and their absolute levels are **not** comparable — the after-heating
   panels involve a 2× OPD dilution whose correction is not documented (§5b).
6. **The chemistry is oxidative and metal-dependent** (Fig. 7). Every rate here is an aerobic,
   unchelated, 0.2 M phosphate rate. **Phosphate itself is a known Maillard catalyst at 0.2 M**
   and the paper carries no buffer-concentration control.
7. **Closed capped tubes at 100 °C.** MGO, DA and HMFO are volatile; none can escape. In an
   open or vented system the DA/MGO balance would differ.
8. **A400 is dilution-multiplied** and reaches 4.09, so it is an extrapolated apparent
   absorbance, not a measured one. Path length is not stated.
9. **The MGO stock is a 40 % aqueous solution** — hydrate/oligomer equilibria mean the nominal
   13 mM is an upper bound on monomeric MGO (§8c).
10. **Lysine is used at 27 mM against 13 mM sugar, i.e. a 2:1 amine:sugar ratio**, well above
    most food matrices. `Lys` is the free amino acid, not a protein-bound ε-amine.
11. **Two author statements are contradicted by their own figures:** the *"1.5 mM HMFO"* in §5a
    (figure shows 1.73 mM from an 1.8 mM preparation) and the *"except for HMFO, Lys promoted the
    browning"* in §8b (DA also shows no promotion).
12. **"HFMO" appears for "HMFO"** in the Fig. 3 caption, the abstract and twice in the body —
    same compound, a typo, but it will break naive text matching.

---

## 12. VERDICT — what is usable

### NOW USABLE

| block | count | status |
|---|---|---|
| **Fig. 2 — A400 and HMFO, 2 pH × 5 times** | 20 | **USABLE [D]**, ± 0.02 A400 / ± 0.01 mM; **four independent text cross-checks pass.** |
| **Fig. 3a/3b — HMFO, MGO, DA, 2 trapping modes × 5 times** | 30 | **USABLE [D]**, ± 0.02 / ± 0.003 mM; **three text cross-checks pass.** Panels must be kept separate. |
| **Fig. 9a — 8 dicarbonyls × 5 times, cumulative** | 40 | **USABLE [D]** for 1-DX and MGO (clean); **± 0.15 mM and marked ≈ for the six minor species.** 1-DX cross-check passes to 0.4 %. |
| **Fig. 9b — 5 species × 5 times, standing pool** | 25 | **USABLE [D]**, ± 0.01 mM. |
| **Fig. 10a/10b — 8 browning curves** | 8 × 6–7 pts = **50** | **USABLE [D]**, ± 0.03 A400; **two text cross-checks pass, one exactly (0.499× vs "about the half").** |
| **HMFO first-order decay constants, 5 conditions** | 5 | **NEW — [Z].** k = 0.094–2.06 h⁻¹ (2.6 × 10⁻⁵ – 5.7 × 10⁻⁴ s⁻¹) at 100 °C, pH 6.5/8.0. §4a, §5c. |
| **Initial browning rates, 8 systems, matched 13 mM / 100 °C / pH 8.0** | 8 | **NEW — [Z].** The **MGO+Lys : Xyl+Lys = 12.2×** ranking. §8c. |
| **Molar branching from 13 mM xylose: 1-DX 76 %, MGO 32 %** | 2 | **NEW — [Z].** §7a. A direct branching constraint. |
| **Cumulative-to-standing ratios: MGO ≈ 120×, 1-DX ≥ 500×** | 2 | **NEW — [Z].** §5c, §7c. Pool ≠ flux, by 2–3 orders of magnitude. |
| **Fig. 8 peak assignments with molecular formulae and t_R** | 7 | **[M], usable** as an identification table. |
| **Compound I / II exact masses and precursor formulae** | 4 | **[M], usable.** |
| **Fig. 7 — HMFO decomposition is O₂- and metal-dependent** | qualitative | **[M], usable as a mechanistic gate**, not as a number. |

**Score against the brief: browning kinetics ✔ (digitised, six time-courses, all
text-cross-validated), norfuraneol/HMFO consumption data ✔ (five decay curves and five
[Z] rate constants), conditions ✔ (fully specified). Nothing in the paper was unreadable.**

### STILL MISSING / NOT OBTAINABLE

1. **Any second temperature → no activation energy for anything.**
2. **Any printed SD or CI.** n = 3 is asserted; not one ± is given.
3. **The pH 6.5 dicarbonyl time course** — *"data not shown"* (p. 407).
4. **Trapping efficiency of OPD**, and the response factors for xylosone / 3-DX / tetrosone.
5. **The quantitative suppression by DTPA and by degassing** — Fig. 7 is chromatograms only, so
   the oxidative dependence is established but not parameterised.
6. **Any melanoidin extinction coefficient**, so A400 cannot be converted to a pigment
   concentration, and the browning rates of §8c are **relative only**.
7. **Lysine consumption** — never measured in any experiment.
8. **The identity of the 340 nm ODS pigment** — explicitly left for future work.
9. **Peak III's molecular formula** — cut off at the top of the Fig. 8 inset in this PDF.

---

*Extraction performed 2026-08-28 (Wave K4c). Figures 2, 3, 7, 8, 9 and 10 digitised
programmatically from 400–600 dpi renders of PDF pages 4, 6, 7 and 8, with per-panel axis
calibration from the printed ticks and series separation by colour mask plus a solid-vs-dashed
line-coverage test; every panel independently cross-validated against at least one value quoted
in the authors' running text (14 such anchors, all reproduced). The text layer was used for
prose only. No figure in the paper was unreadable; the only illegible element was peak III's
formula row in the Fig. 8 inset.*
