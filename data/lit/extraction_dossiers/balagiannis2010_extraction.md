# Balagiannis, Howard, Parker, Desforges & Mottram 2010 (10.1021/bk-2010-1042.ch002) — per-paper extraction 2026-09-07

**Source PDF:** `data/articles/balagiannis2010.pdf` (13 pp., 333,897 bytes, SHA-256 `7e95e2a091aca9d0bc42faf8d63d21215643a3f460cf8aa3011aaceb002b0a45`). ACS book-chapter PDF with a clean, born-digital text layer (producer SPDF / iText 4.2.0; not OCR).
Read method: **both** — full text layer read directly, **plus** 200 dpi rasters of PDF pages 7, 8, 10, 11 (`bal_p7-07.png` … `bal_p11-11.png`, scratchpad) for Figure 1, Figure 2 + equations 1–9, Figures 3–4, and Table I. Figures 1–4 and the nine ODEs are **raster images with no text layer**; Table I is in the text layer and was raster-verified cell by cell.

**Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or any other dossier was touched.**

---

## 0. IDENTITY

| field | value | how verified |
|---|---|---|
| Authors | Dimitrios P. Balagiannis¹, Jack Howard¹, Jane K. Parker¹, Neil Desforges², Donald S. Mottram\*¹ | p.13 raster + text |
| Affiliations | ¹ Dept. of Food and Nutritional Sciences, University of Reading, Whiteknights, Reading RG6 6AP, UK; ² Waltham Centre for Pet Nutrition, Waltham-on-the-Wold, Melton Mowbray, Leicestershire LE14 4RS, UK; \*d.s.mottram@reading.ac.uk | p.13 |
| Title | **"Kinetic Modeling of the Formation of Volatile Compounds in Heated Beef Muscle Extracts Containing Added Ribose"** | p.13 |
| Venue | *Controlling Maillard Pathways To Generate Flavors*; Mottram, D. S., Taylor, A. J., Eds.; **ACS Symposium Series 1042**, **Chapter 2, pp. 13–25** (13 pp.); American Chemical Society, Washington DC, 2010 | running footer every page |
| DOI | **10.1021/bk-2010-1042.ch002** | page-footer stamp every page |
| Dates | Publication Date (Web): **August 10, 2010**. **[NEG] no received/accepted dates printed** | footer stamp; whole-document sweep |
| Provenance stamp | "Downloaded by UNIV OF MICHIGAN ANN ARBOR on September 15, 2015 \| http://pubs.acs.org" | every page |
| Funding | BBSRC research studentship for DB "in conjunction with Waltham Centre for Pet Nutrition (Mars Petcare UK Ltd)" | p.24 Acknowledgments |
| Conflicts | **[NEG] none stated** | — |
| PDF character | born-digital; text layer complete and reliable for all prose and Table I; **Figures 1–4 and equations (1)–(9) are raster only**; **no SI, no tabulated time-course data anywhere in the chapter** | — |
| Companion paper | ref. (11): Balagiannis et al., *J. Agric. Food Chem.* **2009**, 57, 9916–9922 — the ox-liver-extract model (Figure 3) whose parameters are said to be "of the same order of magnitude" (p.22) but are **not reprinted here** [NEG] | p.17, p.22 |

**Correct file for its expected identity.**

---

## 1. ONE-PARAGRAPH VERDICT — READ THIS BEFORE USING ANY NUMBER HERE

This is a **multiresponse (Athena Visual Studio) kinetic fit of a nine-step lumped scheme to one temperature (130 °C), one matrix (1:1 beef-sirloin/water extract), two ribose levels (native 1.16 mmol/kg and ≈10×), and 8 nominal time points from 0 to 90 min**. It fits **only five measured responses**: glucose, ribose, furfural, 3-methylbutanal, 2-methylbutanal — all in **mmol/kg extract**, the three volatiles by **standard-addition calibration (absolute)**. It yields **nine rate constants with 95 % HPD intervals (Table I)**, of which the two cleanest are **k5 (ribose disappearance, 5.20 × 10⁻² min⁻¹, ±5 %)** and **k1 (glucose disappearance, 7.52 × 10⁻³ min⁻¹, ±22 %)**, both **pseudo-first-order in sugar with amino compounds assumed in excess**. **No sulfur volatile of any kind is quantified** — no thiol, no cysteine, no thiophene, no methional; "thiazoles" and "dimethyl sulfides" are only named as classes present (p.16). **No pyrazine is quantified** (one is discussed qualitatively). **No pH is reported anywhere** [NEG]. **No temperature ladder, hence no Ea.** INT1/INT2 are **unquantified fitted latent pools**; the authors say so in "Limitations" (p.23). Two internal inconsistencies to carry forward (§7): the printed **units of k3/k4 are dimensionally inverted for the second-order rate law they multiply**, and **Rleu/Rile are announced as Table I entries but are absent from Table I, with the prose values disagreeing with the slopes printed in Figure 2d/e**. Net: a **validation-shaped dataset for the trunk's sugar-consumption and Strecker-aldehyde lanes at 130 °C in a real-food matrix**, contributing exactly one first-class number the repository lacks — a **ribose disappearance rate at 130 °C in an aqueous meat matrix** — and nothing to the sulfur, acrylamide, or lipid lanes.

---

## 2. SYSTEM DEFINITION `[M]` — verbatim

### 2.1 Meat extract (p.14, "Preparation and Cooking of Meat Extract")

> "Beef muscle taken from sirloin steak was sliced, mixed with **an equal quantity of deionized water** and homogenized by blending it for 1 min using a domestic food processor (Megamix Cuisine Système 5100). The mixture was then centrifuged for 20 min at 14,000 rpm at 4°C and the supernatant was filtered (Whatman filter number 3) under vacuum. **The ribose concentration was determined (1.16 mmol/kg)** and the extract was split into two portions. **D-Ribose was added to one portion to increase its concentration by 10-fold**; the other portion was unchanged."

| variable | value as printed | note |
|---|---|---|
| Matrix | aqueous supernatant of 1:1 (w/w) beef sirloin : deionised water, centrifuged and filtered | "an equal quantity" — mass basis implied, not stated |
| Native ribose | **1.16 mmol/kg** (extract basis) | the only precursor concentration printed as a number |
| Added-ribose portion | "increase its concentration by 10-fold" → nominal **≈ 11.6 mmol/kg** `[Z]` | Figure 1a t=0 triangle reads **≈ 9.5–10 mmol/kg** `[fig]` — the two do not quite agree; see §7 |
| Native glucose | **not printed as a number**; Figure 1a t=0 reads ≈ **5.3 mmol/kg** `[fig]` | "glucose had a higher concentration than ribose" (p.15) |
| Other sugars | mannose, fructose "present in low concentration which did not change significantly during the heating period" (p.15) | **no values** [NEG] |
| Free amino acids | total of 20 individual FAAs by GC-MS; Figure 2a t=0 ≈ **12–13 mmol/kg** `[fig]` | "at higher molar concentrations than the sugars" (p.16) |
| **pH** | **[NEG] never measured, adjusted, or mentioned** | text-layer sweep: zero occurrences of "pH" |
| Water activity / moisture | not stated (dilute aqueous extract, a_w ≈ 1 implied) | — |

### 2.2 Heating protocol (p.14)

> "Aliquots from both meat extracts (**20 mL**) were sealed in **30 mL glass ampoules**, immersed in an **oil bath at 130°C** and heated for different time intervals from **5 to 90 min**. **At least two replicates** were prepared for each heating time. After heating, each ampoule was immersed in a **coolant at –50°C** to stop the reaction."

| variable | value |
|---|---|
| Temperature | **130 °C**, single level; **no ladder → no Ea in this paper** [NEG] |
| Vessel | sealed 30 mL glass ampoule, 20 mL fill (≈ 33 % headspace), oil bath |
| Time points | prose: "5 to 90 min". Figure 1/2 marker positions: **0, 5, 10, 15, 20, 40, 60, 90 min** `[fig]` (8 nominal points; not listed in the text) |
| Come-up time / thermal lag | **not reported** [NEG] |
| Replication | "at least two replicates" per time (duplicate markers visible in Fig 1) |
| Quench | −50 °C coolant |

---

## 3. ANALYTICAL METHOD AND QUANTIFICATION BASIS `[M]`

### 3.1 Volatiles (p.14–15)

> "Dynamic headspace analysis … Homogenized heated extract (**5 g**) and HPLC grade water (**10 mL**) were placed in a 250 mL conical flask with a Dreschel head and the volatiles were swept onto **Tenax** absorbent by nitrogen gas. … Perkin-Elmer Clarus 500 GC-MS … automated thermal desorption unit (Turbomatrix ATD, using a **DB5** non-polar column (**60 m × 0.32 mm i.d., 1 µm film**) …"

> "**Quantification of the volatiles, which were used for modeling, was performed by generating calibration curves using the standard addition method.** For all the other compounds, approximate quantification was obtained by comparison with the internal standard, as described by Methven et al. (13)." (p.15)

| item | status |
|---|---|
| Modelled volatiles (furfural, 3-methylbutanal, 2-methylbutanal) | **absolute**, standard-addition calibration, reported in **mmol/kg** (Figure 1 axes) |
| All other volatiles (>100 identified) | **semi-quantitative** vs. an internal standard; **the internal standard is not named in this chapter** [NEG] (delegated to ref. 13); **no values are printed for any of them** [NEG] |
| Identification | MS + LRI (C6–C25 alkanes) vs authentic compounds / published data / NIST 2.0a |

### 3.2 Sugars and free amino acids (p.15)

| analyte | method | standard |
|---|---|---|
| Sugars | Dionex 8220i IC, CarboPac PA10, isocratic 96 % H₂O / 4 % 400 mM NaOH for 30 min, PAD (420 ms 0.05 V / 180 ms 0.75 V / 420 ms −0.15 V) | **trehalose internal standard**; standards of glucose, fructose, sucrose, ribose, maltose, mannose |
| Free amino acids | 1 g sample + 10 mL 0.01 M HCl, 15 min stir, settle, centrifuge 7200 g; **EZ-Faast** derivatisation, GC-MS Agilent 5975 EI (ref. 14) | 20 individual amino acids summed for "total" |

### 3.3 Modelling (p.15, p.19–21)

> "Multiresponse modeling was performed using **Athena Visual Studio** software package" — parameters "estimated along with their **95% highest posterior density (HPD) intervals**" (p.21). **[NEG] No objective function, no residual plots, no goodness-of-fit statistic is printed**; the only fit-quality statement is "the fit of the model to the experimental data is satisfactory. Only glucose presented a less acceptable fit and this is the reason that the HPD intervals for k2 and k7 were high" (p.21).

---

## 4. THE REACTION SCHEME `[M]` (Figures 3 and 4, p.22, raster-read)

### 4.1 Species groups, verbatim definitions (p.17–18)

- **INT1**: "a group of kinetically important intermediates (INT1) was formed. Probably, these intermediates were **Amadori type products and/or their breakdown products such as deoxyhexosuloses** (26)." … "compounds that contain the **whole carbon skeleton of their parent sugar**". Split into **INT1Glu** (from glucose, k1) and **INT1Rib** (from ribose, k5).
- **INT2**: "a second group of intermediates (INT2), probably **short chain dicarbonyl compounds**" … "**very short-lived**, because they react very fast, through **diffusion controlled reactions**". A single shared pool fed by both INT1Glu (k2) and INT1Rib (k6).
- **M**: "amino acid-specific Maillard products" from INT2 + all other free amino acids — "Strecker compounds from other amino acids … as well as melanoidins, which incorporate both free amino acids and INT2 into the molecule in a ratio of approx 1:1".
- **M′**: "other Maillard type products (M′), such as protein and/or peptide bound products and condensation products with sugars, sugar fragments or more advanced stage Maillard reaction products" — INT2 reacting **without** amino acids.
- **Furfural**: from **both** INT1Glu (k7) and INT1Rib (k8); degrades by pseudo-first-order k9.
- **3-methylbutanal / 2-methylbutanal**: from INT2 + leucine / isoleucine ("FAST"); degrade by **second-order** reaction with INT1 (k3, k4) — "3-methylbutanal and 2-methylbutanal were found to give better models when their reaction to other products changed from first order to second order, with the participation of both ribose- and glucose-derived INT1 type intermediates" (p.18); a version with INT2 as the aldehyde sink "responded with undetermined estimates" and was dropped (p.18–19).

### 4.2 Step table (Figure 4, p.22)

| step | reaction as drawn | rate constant | order | note |
|---|---|---|---|---|
| 1 | Glucose (+ R-NH₂) → INT1Glu | k1 | pseudo-1st in [Glu] | amino compounds "assumed to be in excess during the whole cooking period" (p.23) |
| 2 | INT1Glu → INT2 | k2 | 1st | INT2 then partitions FAST |
| 3 | 3-methylbutanal + INT1 → other products | k3 | **2nd** ([3MeBut]·([INT1Glu]+[INT1Rib])) | |
| 4 | 2-methylbutanal + INT1 → other products | k4 | **2nd** | |
| 5 | Ribose → INT1Rib | k5 | pseudo-1st in [Rib] | |
| 6 | INT1Rib → INT2 | k6 | 1st | |
| 7 | INT1Glu → furfural | k7 | 1st | |
| 8 | INT1Rib → furfural | k8 | 1st | |
| 9 | Furfural → other products | k9 | pseudo-1st | "contribute significantly to the color and flavor development" (p.18) |
| FAST | INT2 + Leu → 3-methylbutanal (fraction Rleu·Fleu); INT2 + Ile → 2-methylbutanal (Rile·File); INT2 + other AA → M (1 − RleuFleu − RileFile); INT2 → M′ | none | — | not rate-limiting; INT2 never appears as a state variable |

### 4.3 The ODE system as printed (equations 1–9, p.19–20, raster-read `[F]`)

```
(1) d[Glu]/dt      = −k1[Glu]
(2) d[Rib]/dt      = −k5[Rib]
(3) d[Int1Glu]/dt  = k1[Glu] − k2[Int1Glu]RleuFleu − k2[Int1Glu]RileFile − k2[Int1Glu](1 − RleuFleu − RileFile)
                     − k2[Int1Glu] − k3[3MeBut][Int1Glu] − k4[2MeBut][Int1Glu] − k7[Int1Glu]
(4) d[Int1Rib]/dt  = k5[Rib] − k6[Int1Rib]RleuFleu − k6[Int1Rib]RileFile − k6[Int1Rib](1 − RleuFleu − RileFile)
                     − k6[Int1Rib] − k3[3MeBut][Int1Rib] − k4[2MeBut][Int1Rib] − k8[Int1Rib]
(5) d[3MeBut]/dt   = k2[Int1Glu]RleuFleu + k6[Int1Rib]RleuFleu − k3[3MeBut]([Int1Glu] + [Int1Rib])
(6) d[2MeBut]/dt   = k2[Int1Glu]RileFile + k6[Int1Rib]RileFile − k4[2MeBut]([Int1Glu] + [Int1Rib])
(7) d[M]/dt        = k2[Int1Glu](1 − RleuFleu − RileFile) + k6[Int1Rib](1 − RleuFleu − RileFile)
(8) d[M′]/dt       = k2[Int1Glu] + k6[Int1Rib]
(9) d[Furf]/dt     = k7[Int1Glu] + k8[Int1Rib] − k9[Furf]
```

`[D]` **Read (3) and (4) literally**: the three branch terms sum to k2[Int1Glu], and the separate M′ term is a *second* k2[Int1Glu]. So the **total first-order drain of INT1Glu is 2·k2 + k7 (≈ 0.337 min⁻¹)** and of **INT1Rib is 2·k6 + k8 (≈ 7.04 × 10⁻³ min⁻¹)** `[Z]`. Any re-implementation must reproduce this doubled sink or the fitted k2/k6 will not reproduce the Figure 1 curves. The authors do not comment on it.

`[M]` Rleu, Rile (p.20–21): "3-methylbutanal and 2-methylbutanal are formed in proportion to the ratio of the parent amino acids to the total amino acid concentration … the constants Rleu and Rile express these ratios … As shown in Figure 2(d,e), these were determined to be **0.0778 and 0.0425 for the extract with no added ribose and 0.0774 and 0.0432 for the extract with added ribose**." Fleu, File: "the ratio of leucine and isoleucine that were converted to the corresponding methylbutanals." The factor 1 − RleuFleu − RileFile: "the ratio of the concentration of all the amino acids present in the studied system, except the leucine and isoleucine, to the total amino acid concentration."

---

## 5. FITTED PARAMETERS — Table I, p.23, re-typed `[F]` (text layer, raster-verified all 11 rows)

"Table I. Optimal estimates and 95% Higher Posterior Density intervals for the parameters which comprise the model of Figure 4"

| parameter | unit as printed | optimal estimate ± 95 % HPD (rel.) | derived `[Z]` |
|---|---|---|---|
| **k1** | min⁻¹ | **7.52E-03 ± 1.68E-03 (22 %)** | glucose t½ = **92 min** (HPD band 75–119 min); 50.8 % glucose remaining at 90 min |
| **k2** | min⁻¹ | **1.68E-01 ± 1.29E-01 (77 %)** | INT1Glu drained at 2k2+k7 = 0.337 min⁻¹ → t½ ≈ 2.1 min (quasi-steady) |
| **k3** | mmol kg⁻¹ min⁻¹ *(sic; see §7)* | **4.28E-03 ± 1.17E-03 (27 %)** | 3-MeBut sink |
| **k4** | mmol kg⁻¹ min⁻¹ *(sic)* | **1.16E-03 ± 6.09E-04 (53 %)** | 2-MeBut sink; k3/k4 = 3.69 |
| **k5** | min⁻¹ | **5.20E-02 ± 2.73E-03 (5 %)** | **ribose t½ = 13.3 min** (13.0–13.9); 0.93 % remaining at 90 min; **k5/k1 = 6.9** |
| **k6** | min⁻¹ | **3.50E-03 ± 1.06E-03 (30 %)** | INT1Rib drained at 2k6+k8 = 7.04E-3 min⁻¹ → t½ ≈ 98 min (**accumulates**); k2/k6 = 48 |
| **k7** | min⁻¹ | **9.39E-04 ± 7.65E-04 (81 %)** | furfural from INT1Glu |
| **k8** | min⁻¹ | **4.16E-05 ± 9.68E-06 (23 %)** | furfural from INT1Rib; k7/k8 = 22.6 |
| **k9** | min⁻¹ | **2.73E-02 ± 8.58E-03 (31 %)** | furfural t½ = 25.4 min |
| **Fleu** | — | **7.89E-02 ± 1.37E-02 (17 %)** | Rleu·Fleu = 0.0774 × 0.0789 = **6.1 × 10⁻³** |
| **File** | — | **5.56E-02 ± 9.87E-03 (18 %)** | Rile·File = 0.0432 × 0.0556 = **2.4 × 10⁻³** |
| Rleu | — | **NOT IN TABLE I** — prose only: 0.0778 (no added Rib) / 0.0774 (added) | Fig 2d prints y = 0.0775x; Fig 2e prints y = 0.0774x `[fig]` |
| Rile | — | **NOT IN TABLE I** — prose only: 0.0425 (no added Rib) / 0.0432 (added) | Fig 2d prints y = **0.044x**; Fig 2e prints y = **0.0428x** `[fig]` |

`[Z]` Only **0.85 %** of INT2 flux goes to the two methylbutanals (1 − RleuFleu − RileFile = 0.9915 goes to M). Which R pair (with/without ribose) was used in the single fitted model is **not stated** [NEG].

The authors' own reading of the ordering (p.21–22): k5 > k1 ("ribose is more reactive than glucose"); k6 < k2 ("INT1rib are less reactive than the corresponding glucose related intermediates and accumulate", citing Hofmann 1999 on xylose vs glucose fragments); k3 > k4 ("3-methylbutanal degrades faster than 2-methylbutanal"); k7 > k8 ("this model predicts that more furfural is derived from glucose than ribose. This is an interesting estimation but needs further investigation, since it was expected that furfural would be formed faster from ribose").

---

## 6. MEASURED TIME COURSES — **all figure-only, none tabulated** `[fig]`

**No concentration–time table exists in the chapter.** Values below are digitised from the 200 dpi raster of Figure 1 (p.19) and Figure 2 (p.20) by eye against the printed axis ticks; treat as **±10 % of full scale**. Duplicate markers at some times are both listed. Nominal time grid: 0, 5, 10, 15, 20, 40, 60, 90 min.

### 6.1 Figure 1a — sugars, mmol/kg extract

| t (min) | ribose, +Rib (▲) | glucose, +Rib (◆) | ribose, native (∗) | glucose, native (■) |
|---|---|---|---|---|
| 0 | ≈ 9.7 (and ≈ 8.0) | ≈ 5.3 | ≈ 1.2 | ≈ 5.3 (and 4.9) |
| 5 | ≈ 8 | ≈ 4.8 | ≈ 1.0 | ≈ 4.5 |
| 10 | ≈ 5.3 | ≈ 4.6 | ≈ 0.8 | ≈ 4.4 |
| 15 | ≈ 4.5 | ≈ 4.4 | ≈ 0.6 | ≈ 3.8 |
| 20 | ≈ 3.5 | ≈ 4.3 | ≈ 0.5 | ≈ 3.5 |
| 40 | ≈ 1.5 | ≈ 4.4 | ≈ 0.1 | ≈ 2.5 |
| 60 | ≈ 0.4 | ≈ 4.3 | ≈ 0 | ≈ 2.3 |
| 90 | ≈ 0.1 | ≈ 3.5 (and 2.7) | ≈ 0 | ≈ 2.8 |

Model curves: ribose (both extracts) exponential to ≈ 0 by 60 min; glucose model falls to ≈ 2.7 at 90 min in both. The +Rib glucose points at 40–60 min sit **visibly above** the model — this is the "less acceptable fit" the authors concede.

### 6.2 Figure 1b–d — modelled volatiles, mmol/kg extract

| t (min) | furfural +Rib (◆) | furfural native (■) | 3-MeBut +Rib (■) | 3-MeBut native (◆) | 2-MeBut +Rib (◆) | 2-MeBut native (■) |
|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5 | ≈ 0.0003 | ≈ 0.0001 | ≈ 0.0005 | ≈ 0.0003 | ≈ 0.0002 | ≈ 0.0001 |
| 10 | ≈ 0.0008 | ≈ 0.0003 | ≈ 0.0010 | ≈ 0.0008 | ≈ 0.0005 | ≈ 0.0003 |
| 15 | ≈ 0.0015 | ≈ 0.0005 | ≈ 0.0020 | ≈ 0.0012 | ≈ 0.0008 | ≈ 0.0005 |
| 20 | ≈ 0.0045 | ≈ 0.0015 | ≈ 0.0030 | ≈ 0.0020 | ≈ 0.0015 | ≈ 0.0008 |
| 40 | ≈ 0.0090 | ≈ 0.0035 | ≈ 0.0055 | ≈ 0.0045 | ≈ 0.0028 | ≈ 0.0015 |
| 60 | ≈ 0.0085–0.0095 | ≈ 0.0040 | ≈ 0.0062 | ≈ 0.0058 | ≈ 0.0038 | ≈ 0.0020 |
| 90 | ≈ 0.0105–0.0120 | ≈ 0.0035 | ≈ 0.0070 | ≈ 0.0078 | ≈ 0.0055 | ≈ 0.0035–0.0040 |

Caption note `[D]`: the Figure 1c caption text is scrambled ("with (■ values, ─ model) and without (◆ values, ···· model) with added ribose"); the prose (p.16–17: added ribose ⇒ 3-MeBut levels off and ends **lower**) fixes ■ = added ribose, ◆ = native, as tabulated. Unit conversion `[Z]` (MW 86.13 for both methylbutanals, 96.08 for furfural): 0.0075 mmol/kg 3-MeBut ≈ 0.65 mg/kg; 0.0055 mmol/kg 2-MeBut ≈ 0.47 mg/kg; 0.012 mmol/kg furfural ≈ 1.15 mg/kg. Molar furfural yield on consumed ribose (+Rib, 90 min) ≈ 0.012 / ≈ 9.5 ≈ **0.13 %** `[Z, fig-based]`.

### 6.3 Figure 2a–c — free amino acids, mmol/kg extract

| t (min) | total FAA +Rib (◆) | total FAA native (■) | Leu +Rib | Leu native | Ile +Rib | Ile native |
|---|---|---|---|---|---|---|
| 0 | ≈ 13 | ≈ 11.5 | ≈ 1.0 | ≈ 0.9 | ≈ 0.50 | ≈ 0.42 |
| 5 | ≈ 12.5 | ≈ 11.5 | ≈ 0.9 | ≈ 0.8 | ≈ 0.48 | ≈ 0.45 |
| 10 | ≈ 11.5 | ≈ 11 | ≈ 0.9 | ≈ 0.8 | ≈ 0.50 | ≈ 0.48 |
| 15 | ≈ 8.5 | ≈ 9 | ≈ 0.7 | ≈ 0.75 | ≈ 0.40 | ≈ 0.45 |
| 20 | ≈ 10 | ≈ 10 | ≈ 0.85 | ≈ 0.85 | ≈ 0.50 | ≈ 0.45 |
| 40 | ≈ 10 | ≈ 7 | ≈ 0.8 | ≈ 0.6 | ≈ 0.45 | ≈ 0.33 |
| 60 | ≈ 9.5 | ≈ 7.5 | ≈ 0.75 | ≈ 0.65 | ≈ 0.42 | ≈ 0.38 |
| 90 | ≈ 10 | ≈ 10 | ≈ 0.8 | ≈ 0.8 | ≈ 0.47 | ≈ 0.45 |

Authors' reading (p.16): "initial slow decrease in total amino acid concentration was followed by a leveling off and/or a small increase" — attributed to sugar reaction, regeneration from Amadori products, reaction with advanced products, and "it is probable that **proteolysis** affected the amino acid concentration as well." **Free amino acids are not state variables in the model** — they enter only through the fixed ratios Rleu/Rile.

### 6.4 Pyrazines and everything else — **[NEG] no numbers**

"More than 100 compounds were identified … pyrazines, thiazoles, pyrroles, dimethyl sulfides, furans and furanones." (p.16). Only two qualitative statements: lipid-derived compounds "not affected by the addition of ribose"; Maillard compounds increased, "An exception was **2-ethyl-3,6-dimethylpyrazine** where the addition of ribose caused a reduction in its concentration." **No pyrazine, thiazole, sulfide, furanone or lipid-oxidation value is printed anywhere.**

---

## 7. INTERNAL INCONSISTENCIES AND DEFECTS `[D]`

| # | issue | evidence | severity |
|---|---|---|---|
| C1 | **k3, k4 units printed as "mmol kg⁻¹ min⁻¹"**, but eqs (5)/(6) use them as second-order constants multiplying [aldehyde]·[INT1] (both mmol/kg). Dimensional consistency requires **kg mmol⁻¹ min⁻¹** — the printed unit is the reciprocal. | Table I raster vs eq (5)/(6) raster | ★★ — correct the unit before reuse; the numbers themselves are presumably right in kg mmol⁻¹ min⁻¹ |
| C2 | Prose (p.21) says Rleu, Rile "were estimated along with their 95% HPD intervals (Table I)". **Table I has no Rleu/Rile rows.** Prose gives them as fixed regression slopes from Fig 2d/e instead. | Table I raster (11 rows) | ★★ — R values are fixed inputs, not fitted; which pair the single model used is unstated |
| C3 | Prose Rile for the no-added-ribose extract = **0.0425**; Figure 2d prints **y = 0.044x** (3.5 % off). Prose Rile (added) = 0.0432; Fig 2e prints 0.0428. Rleu: 0.0778 vs printed 0.0775. | Fig 2d/e raster | ★ — small, but the paper contradicts its own figure |
| C4 | "Increase … by 10-fold" of 1.16 mmol/kg → 11.6 (or 12.8 if "add 10×"); Figure 1a t=0 ribose reads ≈ 9.5–10 `[fig]`. | Fig 1a raster | ★ — ~15 % gap between nominal and plotted initial ribose |
| C5 | Equations (3)/(4) contain the M′ sink as a **second full k2 (k6) term**, so INT1 → INT2 flux is effectively 2k2 (2k6). Not discussed. | eq (3)/(4) raster | ★★ — matters for any re-implementation (see §4.3) |
| C6 | Figure 1c caption legend is scrambled (see §6.2). | Fig 1 caption | ★ |
| C7 | Reported HPD for k2 (77 %) and k7 (81 %) — the authors themselves flag these as poorly determined; "Only glucose presented a less acceptable fit". | p.21, p.23 | informational — do not use k2 or k7 as fixed values |

---

## 8. VERIFIED NEGATIVES `[NEG]` — do not re-open this paper for any of these

| item | status |
|---|---|
| Any thiol (2-methyl-3-furanthiol, 2-furfurylthiol, mercaptoketones), methional, any thiophene, cysteine, H₂S | **absent** — zero text-layer hits for "thiol", "cystein", "thiophen"; sulfur classes named only as "thiazoles" and "dimethyl sulfides" without values |
| pH (initial, final, buffer) | **absent** |
| Temperature other than 130 °C; Ea; Q₁₀ | **absent** |
| HMF, 5-hydroxymethylfurfural quantitation | absent (HMF mentioned once, p.17, as a literature furfural precursor) |
| Acrylamide | absent (appears only in reference 7's title) |
| Amadori compounds, deoxyosones, any INT1/INT2 measurement | absent — "the compounds INT1 and INT2 have not been quantified" (p.23) |
| Pyrazine, furanone, thiazole, lipid-oxidation concentrations | absent (qualitative only) |
| Any tabulated concentration–time data | absent — Figures 1 and 2 only |
| Parameter values of the 2009 liver model for comparison | absent ("same order of magnitude", p.22) |
| Goodness-of-fit statistics, residuals, covariance / correlation matrix | absent |
| Come-up time; headspace gas; moisture / a_w; meat pH; meat provenance beyond "sirloin steak" | absent |

---

## 9. WHAT THE REPOSITORY CAN USE

**Lane by lane:**

- **Trunk (glucose/fructose/glycine → Amadori, deoxyosones, HMF)** — *partial, validation-grade.* Gives a **glucose pseudo-first-order disappearance k1 = 7.52 × 10⁻³ min⁻¹ (±22 %) at 130 °C in an amino-rich aqueous meat matrix** (FAA ≈ 12 mmol/kg, i.e. ≈ 2.3× glucose), plus one **furfural** formation/degradation pair (k7 from a hexose pool, k9 = 2.73 × 10⁻² min⁻¹ furfural decay). No Amadori, deoxyosone or HMF measurement. Fructose is present but flat and unmodelled. The k1 value is a **whole-matrix lumped rate** (sugar consumption by all amines including peptides/protein), so it is an **upper-bound-ish aggregate** for the trunk's glucose + amine step at 130 °C, not an elementary constant. Usable as a **hold-out sanity check** on the mass-action model's total glucose consumption at 130 °C, not as a fit target.
- **Sulfur (ribose/cysteine → MFT, FFT)** — *precursor side only.* **Confirmed: NO thiol is quantified** (§8). What the paper *does* pin is the **ribose side of the MFT precursor pair**: **k5 = 5.20 × 10⁻² min⁻¹ (±5 %) at 130 °C, t½ ≈ 13 min, ribose ≈ 99 % gone by 90 min**, and k5/k1 ≈ 7 (pentose : hexose reactivity ratio in the same matrix). A mass-action sulfur lane that consumes ribose via an amine-catalysed step at 130 °C must reproduce a ribose disappearance rate of this order **in the presence of ≈ 12 mmol/kg FAA and no added cysteine**; the paper also gives the qualitative fact that **ribose-derived INT1 (pentose deoxyosone-type pool) accumulates** (k6 ≪ k2), which is the pool MFT chemistry draws on. This is the single most transferable number in the chapter.
- **Acrylamide** — nothing (asparagine not reported; no acrylamide).
- **Lipid** — nothing quantitative; one qualitative statement that lipid-derived volatiles were unaffected by ribose.

**Temperature–time-resolved series a mass-action model could be validated against** (all at 130 °C only, all `[fig]`, §6): (i) ribose decay at two initial levels (≈ 1.2 and ≈ 9.7 mmol/kg) — 8 points each; (ii) glucose decay, two extracts; (iii) furfural rise-and-plateau, two extracts; (iv) 3-methylbutanal and 2-methylbutanal accumulation, two extracts — the **only Strecker-aldehyde time courses in a real meat matrix in the corpus**, with the striking feature that 10× ribose *increases* 2-MeBut but *not* the 90-min 3-MeBut; (v) total FAA, Leu, Ile — flat-ish, not for fitting. **No temperature ladder exists, so nothing here constrains an Ea.**

**What the k values pin:** k5 (and to a lesser degree k1) pin the **sugar-consumption timescale at 130 °C in a high-amine aqueous matrix**; k9 pins furfural's own lifetime at 130 °C (t½ 25 min) — relevant wherever the trunk treats furfural/HMF as terminal; k3/k4 give a rare **Strecker-aldehyde loss** ordering (3-MeBut degrades ≈ 3.7× faster than 2-MeBut, second-order in an INT1-like pool) — but with the unit defect (C1) and the latent-pool dependence, treat only the *ratio* as portable. k2, k6, k7, k8 are latent-pool constants with no independent meaning outside this scheme.

**Proposed role:** hold-out / plausibility check for the trunk and sulfur-precursor lanes at 130 °C; **not** a fit source. If any single number is taken into `k1_kinetic_parameters.md`, it should be k5 with its 5 % HPD, tagged 130 °C / meat extract / pH unreported.

---

## 10. CAVEATS

1. **Undefined matrix.** A 1:1 beef-sirloin/water supernatant: sugars, 20 FAAs, peptides, proteins, nucleotides, phosphate, minerals — none but ribose (1.16 mmol/kg) is printed as a number. "Amino compounds … assumed to be in excess" (p.23) is the authors' own disclaimer; proteolysis is invoked to explain FAA drift.
2. **pH unreported** — and meat-extract pH (typically ≈ 5.5–5.8, *not from this paper*) is neither buffered nor tracked over 90 min at 130 °C. Every k is pH-conditional on an unknown, drifting value.
3. **Single temperature.** No Arrhenius information; nothing here transfers to another T without an external Ea.
4. **Units.** Concentrations are **mmol/kg of extract**, not of meat (extract ≈ 50 % meat by mass, so meat-basis values would be roughly 2× if extraction were complete — an assumption `[Z]`, not the paper's). k3/k4 printed units are dimensionally inverted (C1).
5. **Lumped/latent species.** INT1, INT2, M, M′ are unmeasured; k2, k6, k7, k8, k3, k4 are identifiable only within this scheme and its doubled-sink equations (C5). The authors: "the compounds INT1 and INT2 have not been quantified", HPD for k2/k7 "revealed the need for manipulation of glucose concentration" (p.23).
6. **Data are figure-only.** All §6 values are digitised by eye at 200 dpi; no replicate scatter, no error bars, no table. Do not put them in a benchmark file without re-digitisation at higher resolution and an explicit tolerance.
7. **Time grid is inferred** from marker positions (0, 5, 10, 15, 20, 40, 60, 90 min), not listed by the authors; no come-up time is given for a 20 mL ampoule in a 130 °C oil bath, so early points (5–15 min) carry unknown thermal lag.
8. **Quantitation basis mixed.** The three modelled volatiles are absolute (standard addition); everything else semi-quantitative against an internal standard named only in ref. 13 — and none of it is printed anyway.
9. **The authors' own hedge** (p.23): "all models are wrong but some are useful" … "the present model does have some limitations" … "it would be important to conduct additional experiments in order to achieve a more robust and accurate model".
10. **Rleu/Rile ambiguity** (C2–C3): fixed inputs, two candidate pairs, figure and prose disagree at the 1–4 % level, and it is unstated which pair the fit used.

### Provenance summary for this dossier
`[F]` Table I (11 rows), equations (1)–(9), Figure 3/4 topology, all quoted prose — raster-verified. `[M]` §2–3 methods verbatim. `[fig]` every value in §6 and the initial-concentration estimates in §2.1. `[Z]` half-lives, ratios, products RF, unit conversions, yields — trivial arithmetic on `[F]` or `[fig]` inputs. `[D]` §4.3 doubled-sink reading, §7 C1–C7, lane mapping in §9. `[NEG]` §8 — whole-text-layer sweeps.
