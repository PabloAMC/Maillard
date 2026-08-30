# Kang et al. 2026 — SUPPLEMENTARY INFORMATION — COMPLETE TRANSCRIPTION
### Tables S1–S6 + Figs S1–S4 of `data/articles/Kang2026_supplementary.pdf` (DOI 10.1039/D5FB00932D)
### **THE SULFUR BRANCH'S FIRST TEMPERATURE LADDER — Table S4 recovered in full, Tier A, μg L⁻¹.**

Wave **K4d**, 2026-08-28. Read-only extraction: no repo file outside
`data/lit/extraction_dossiers/` was written or modified. Companion to
`kang2026_extraction.md` (main paper, Wave K4c) — **read that first**; this dossier does not
re-transcribe the main paper's Table 1, Fig. 1c–f or §3 text, and instead **closes the four
gaps that dossier left open** (its §7 "not in this PDF" list).

**Provenance codes:** **[M]** measured and printed · **[D]** digitised by this wave from a
figure · **[Z]** derived by this wave · **[Q]** qualitative · **[!]** integrity flag.

---

## 0. IDENTITY — **THIS IS THE EXPECTED DOCUMENT.**

| field | value |
|---|---|
| Header | *"Supplementary Information (SI) for **Sustainable Food Technology**. This journal is © The Royal Society of Chemistry 2026"* |
| Title | **"Unlocking the Potential of Cysteine-Xylose Maillard Reaction Intermediates as Natural Flavor Precursor: A Comprehensive Study on Flavor Regulation, Storage Stability, and Antioxidant Activity"** |
| Authors | **Di Kang§, Lin Jiang⊥, Songjin Zheng§, Yuan Hu⊥, Haifeng Wang⊥, Teng Li⊥, Yuying Fu⊥\*, Yun Zhai⊥\*** |
| Affiliations | ⊥ School of Food Science and Biotechnology, **Zhejiang Gongshang University**, Hangzhou 310018, P. R. China · § **China Tobacco Hebei Industrial Co., LTD**, Shijiazhuang 050051 |
| Pages | **16**, A4 |
| Producer | **Aspose.Pdf for .NET 9.3.0** (author metadata: `dengshibin`); CreationDate 2022-09-19, ModDate 2026-03-10 |
| Contents | **Tables S1–S6, Figs S1–S4** |

**⇒ Title, journal and author list match `Kang2026.pdf` exactly. Correct SI, correct paper.**
Table S4 and Table S5 — the two tables `kang2026_extraction.md` §7 named as the corpus's #1
and #2 missing artefacts — **are both here, complete.**

---

## 0b. ⚠️ THE μ→m HAZARD (Amendment 4) — RESOLVED, AND THE RESOLUTION IS ASYMMETRIC

**Binding requirement discharged. Every unit in this dossier was verified against a raster
render, not the text layer.** Method: `pdftoppm -r 600 -png` (Tables S2–S5, pages 5–10) and
`-r 400` (Figs S1–S4, pages 13–16), read as images.

| artefact | text-layer says | **raster (600 dpi) says** | verdict |
|---|---|---|---|
| **Table S4 caption** | `μg/L` | **`μg/L` — true Greek mu, verified glyph-by-glyph** | ✅ **NOT corrupted** |
| **Table S4 column header** | `Concentration (μg/L)` | **`Concentration (μg/L)`** | ✅ **NOT corrupted** |
| **Table S5 caption** | `μg/L` | **`μg/L`** | ✅ **NOT corrupted** |
| **Table S5 column header** | `Concentration (μg/L)` | **`Concentration (μg/L)`** | ✅ **NOT corrupted** |
| **S4/S5 footnote a** (column film thickness) | `0.25 μm` | **`0.25 μm`** | ✅ **NOT corrupted** |
| **Fig. S3 y-axis** | *(image, no text layer)* | **`Cys Concentration (mmol/L)`** | ✅ raster-only |
| **Fig. S4 y-axis** | *(image, no text layer)* | **`Cys concentration (mmol/L)`** | ✅ raster-only |

**Per-table basis of verification:**
- **Table S2** — raster, page 5 top. **Table S3** — raster, page 6 top (header row + first two
  data rows); remaining 29 rows from text layer, which was confirmed byte-identical to raster
  on the verified rows.
- **Table S4** — **raster, pages 7 and 8: caption, both header rows, all 16 rows of the first
  sulfur block, all 11 rows of the second sulfur block, the sulfur Subtotal, the full
  N-heterocycle block, and the O-heterocycle Subtotal + Total + all four footnotes.**
  Every raster reading matched the text layer to the last decimal. **Table S4 is verified on a
  raster basis end-to-end.**
- **Table S5** — **raster, pages 9 and 10: caption, both header rows, the first two thiol rows,
  and the complete O-heterocycle block + "Total content of sulfur-containing compounds" + Total
  + footnotes.** Matched text layer exactly. Interior rows from text layer, on a text basis
  validated by the arithmetic check in §4b (every one of the nine subtotals reconciles to
  ±0.001 from the transcribed rows — a self-check that would fail on any digit error).
- **Table S6, Figs S1–S4** — raster.

### **⚠️ THE HAZARD IS REAL, BUT IT LIVES IN THE MAIN PDF, NOT THIS ONE.**
This SI was produced by **Aspose.Pdf**, which embeds a real `μ`. The main article
`Kang2026.pdf` is RSC-typeset and **is** corrupted — confirmed this wave by re-running
`pdftotext -layout` on it:

- main text p. 3242 renders as **`reaching 11.039 mg L−1 at 140 °C`** — the true value is
  **11.039 μg L⁻¹** (Table S4, furfural at 140 °C, raster-verified here). **A 1000× error.**
- main text §2.7 renders the SPME fibre as **`CAR/PDMS/DVB, 75 mm`** — the true value is
  **75 μm**.

**⇒ `kang2026_extraction.md`'s μ-correction of the main paper is independently vindicated by
this SI. Any future re-read of `Kang2026.pdf` must apply the same correction; this SI is the
authoritative unit source.**

---

## 1. THE ONE-PARAGRAPH ANSWER

**Table S4 delivers exactly what the corpus's top declared gap asked for: 2-methyl-3-furanthiol
(MFT) and 2-furfurylthiol (FFT), absolutely quantified in μg L⁻¹, at three temperatures
(100 / 120 / 140 °C) at fixed pH, time, and substrate loading, with n = 3 and SDs, and with
their external calibration curves printed in Table S3 (MFT R² = 0.9989, FFT R² = 0.9992).**
Both are Tier A by the main paper's asterisk footnote **and** by the independent test of
appearing in Table S3. **The gap is closed for MFT and FFT** — with three material caveats
(§7): (i) the ladder is a **yield-at-120-min** ladder, not a rate ladder, and the precursor is
demonstrably not exhausted at 100 °C (Fig. S4, §5b), so the apparent Eₐ mixes rate and extent;
(ii) MFT and FFT are **strongly non-Arrhenius** — nearly flat from 100→120 °C then a 4.3× /
2.8× jump to 140 °C (§6a), so a single Eₐ cannot reproduce the ladder; (iii) **the reported SDs
are not trustworthy** — Table S5's pH-7 column reproduces Table S4's 100 °C column to the last
decimal in every mean but with *different* standard deviations (§4c **[!]**). **Use the means;
do not use the SDs as measurement uncertainty.** Separately, **Table S4 confirms the Fig. 1a
corruption diagnosed in `kang2026_extraction.md` §5a and supplies the true values** (§3), and
**Figs S3/S4 are not H₂S evidence — no H₂S is measured anywhere in this SI** (§5) — but they
are something better: two digitisable Cys time-courses at three temperatures, from which this
wave derives **Eₐ(Cys depletion) = 55.1 kJ mol⁻¹, R² = 0.994** (§5b).

---

## 2. TABLES S1, S2, S3, S6 — the supporting tables

### 2a. Table S1 — *"Derivative fragmentation in the LC-MS/MS spectrum and ¹³C NMR data of TTCA"*
**[!] TABLE S1 IS AN EMPTY TITLE.** The caption is printed on page 4 and **no table body
follows it** — the next printed object is Table S2's caption. Verified on raster (page 4) and
by text layer. **The ¹³C NMR assignment and the MS/MS fragmentation of TTCA — the only
structural proof of the intermediate's identity — are not in this SI.** The main paper cites
"Fig. S1 and Table S1" for TTCA identity (p. 3240); Fig. S1 is present (an HPLC chromatogram,
§5a), Table S1 is not. **TTCA's structure is therefore unverified from the published record.**

### 2b. Table S2 — water activity of the saturated salt solutions **[M]**
*Caption as printed (raster): "Water activity corresponding to each saturated salt solution at 25 °C."*

| saline | **a_w at 25 °C** |
|---|---|
| lithium chloride | **0.113 ± 0.003** |
| potassium carbonate | **0.432 ± 0.004** |
| dicalcium phosphate | **0.843 ± 0.003** |

**Only three salts, only 25 °C.** These are the a_w levels of the storage-stability block
(Table S6). Values agree with the Greenspan reference data for LiCl (0.113) and K₂CO₃ (0.432).
**⚠️ "dicalcium phosphate" at a_w 0.843 is not a standard saturated-salt a_w standard**; the
0.843 level is normally KCl (0.843) — **almost certainly a mis-naming of KCl.** [!] Use the
a_w value, not the salt name.

### 2c. Table S3 — the calibration curves **[M]** — *this is the Tier A/B adjudicator*
*Caption as printed (raster): "The calibration curves of volatile flavor compounds"*.
31 compounds. `y` = response ratio vs the 1,2-dichlorobenzene internal standard; `x` =
concentration (units not stated, but must be μg L⁻¹ for consistency with S4/S5).

| # | compound | calibration curve | R² | x-intercept **[Z]** |
|---|---|---|---|---|
| 1 | Thiophene | y = 0.0223 x + 0.0011 | 0.9991 | −0.049 |
| 2 | 2-Methylthiophene | y = 0.0472 x + 0.0009 | 0.9994 | −0.019 |
| 3 | Thiazole | y = 0.0601 x − 0.0012 | 0.9982 | +0.020 |
| 4 | **2-Methylthiazole** | **y = 0.0084 x − 0.0521** | 0.9923 | **+6.20 ⚠ [!]** |
| 5 | 3-Mercapto-2-butanone | y = 0.1216 x − 0.0007 | 0.9999 | +0.006 |
| 6 | 2,5-Dimethyl-thiazole | y = 0.0637 x + 0.0012 | 0.9979 | −0.019 |
| 7 | **2-Methyl-3-furanthiol (MFT)** | **y = 0.0039 x + 0.0011** | **0.9989** | −0.282 |
| 8 | 2-Ethyl-thiazole | y = 0.0629 x + 0.0041 | 0.9994 | −0.065 |
| 9 | 3-Mercapto-2-pentanone | y = 0.0379 x + 0.0131 | 0.9997 | −0.346 |
| 10 | 2,4,5-Trimethylthiazole | y = 0.0571 x − 0.0007 | 0.9994 | +0.012 |
| 11 | 4,5-Dimethyl-thiazole | y = 0.0297 x − 0.0042 | 0.9969 | +0.141 |
| 12 | **2-Furfurylthiol (FFT)** | **y = 0.0287 x − 0.0003** | **0.9992** | +0.010 |
| 13 | 5-Ethyl-2,4-dimethyl-thiazole | y = 0.0184 x + 0.0008 | 0.9981 | −0.043 |
| 14 | 2-Ethyl-4-methyl-thiazole | y = 0.0430 x + 0.0024 | 0.9923 | −0.056 |
| 15 | 3-Thiophenethiol | y = 0.0285 x + 0.0033 | 0.9954 | −0.116 |
| 16 | 2-Acetylthiazole | y = 0.0364 x − 0.0014 | 0.9967 | +0.038 |
| 17 | 2-Thiophenecarboxaldehyde | y = 0.0446 x + 0.0009 | 0.9998 | −0.020 |
| 18 | 5-Methyl-2-thiophenecarboxaldehyde | y = 0.0296 x − 0.0007 | 0.9999 | +0.024 |
| 19 | 2-Thiophenemethanethiol | y = 0.0517 x + 0.0012 | 0.9954 | −0.023 |
| 20 | 1-(3-Thienyl)-ethanone | y = 0.0728 x − 0.00012 | 0.9981 | +0.002 |
| 21 | 2,5-Thiophenedicarboxaldehyde | y = 0.0423 x − 0.0031 | 0.9997 | +0.073 |
| 22 | Methylpyrazine | y = 0.0398 x + 0.0012 | 0.9949 | −0.030 |
| 23 | Pyrazine | y = 0.0294 x + 0.0012 | 0.9989 | −0.041 |
| 24 | 3-Methyl-pyridine | y = 0.0563 x + 0.0032 | 0.9998 | −0.057 |
| 25 | 2,5-Dimethyl-pyrazine | y = 0.0319 x + 0.0013 | 0.9978 | −0.041 |
| 26 | Furan | y = 0.0253 x − 0.0012 | 0.9968 | +0.047 |
| 27 | 2-Methyl-furan | y = 0.0336 x + 0.0021 | 0.9997 | −0.063 |
| 28 | **Furfural** | y = 0.0347 x − 0.0012 | 0.9923 | +0.035 |
| 29 | 1-(2-Furanyl)-ethanone | y = 0.0629 x + 0.0024 | 0.9948 | −0.038 |
| 30 | 2(5H)-Furanone | y = 0.0598 x − 0.0011 | 0.9993 | +0.018 |
| 31 | 4-Hydroxy-5-methyl-3(2H)-furanone | y = 0.0098 x + 0.0002 | 0.9971 | −0.020 |

**What Table S3 does deliver:** the definitive Tier A roster (§4d), and confirmation that
**MFT and FFT are genuine external-standard quantifications**, not IS-only semi-quant.

**What Table S3 does NOT deliver [!]:** **no calibration range, no LOD, no LOQ, no recovery, no
replicate count, no units on x, no concentration levels used.** Without a range it is
impossible to tell whether the reported 0.02–11.4 μg L⁻¹ values are inside or extrapolated
beyond the calibrated span. **This is the single largest unquantified uncertainty in S4/S5.**

**⚠️ [!] Row 4 (2-methylthiazole) is internally impossible.** Its curve has an x-intercept of
**+6.20 μg L⁻¹**, but every 2-methylthiazole concentration reported in S4 (0.142 / 0.194 /
0.254) and S5 (0.023 / 0.142 / 1.779) is **below** that intercept — back-substituting them
gives a *negative* response ratio. Its slope (0.0084) is also 4–7× below every other thiazole's
(0.030–0.064). Either the intercept is a typo (−0.00521?) or the slope is. **Exclude
2-methylthiazole from Tier A; treat it as Tier B or unusable.** It contributes ≤ 0.6 % of any
S4 sulfur subtotal, so nothing else is affected.

### 2d. Table S6 — retention rates of TTCA / ARP under storage **[M]**
*Caption: "Retention rates of TTCA/ARP under different conditions"*. Units: **%**. No time
axis is stated in the SI; the main paper's storage block is the reference for duration.
ARP = the Amadori rearrangement product comparator.

| block | condition | **TTCA retention (%)** | **ARP retention (%)** |
|---|---|---|---|
| Temperature | **4 °C** | **99.12** | **95.52** |
| | **25 °C** | **95.53** | **94.18** |
| | **40 °C** | **92.94** | **87.83** |
| pH | **5.5** | **95.51** | **94.23** |
| | **7.0** | **96.70** | **90.52** |
| | **9.0** | **88.81** | **78.75** |
| a_w | **0.113** | **86.49** | **78.48** |
| | **0.432** | **72.47** | **57.27** |
| | **0.843** | **64.23** | **39.39** |

**[Z] Structural readings:**
1. **TTCA is more storage-stable than ARP under every one of the nine conditions**, by
   1.5 pp (pH 5.5) to **24.8 pp (a_w 0.843)**. Mean advantage **+8.4 pp**.
2. **a_w is by far the dominant stressor**, not temperature and not pH: over the a_w ladder
   TTCA falls 86.5 → 64.2 % (−22.3 pp) and ARP 78.5 → 39.4 % (−39.1 pp), against −6.2 pp and
   −7.7 pp over the 4→40 °C ladder. **The intermediate's storage risk is a water-activity
   problem.** The a_w 0.113 row already sits ~9 pp below the 25 °C row, so the a_w block is run
   on a drier, harsher basis than the temperature block and the two are not directly comparable.
3. **Alkaline pH is the second stressor** (pH 9: −7.9 pp TTCA, −11.8 pp ARP vs pH 7).
4. **[!] The blocks share no common control point** — no cell is replicated across blocks and
   no reference condition is defined, so the three ladders cannot be composed into a single
   response surface. **Ratios within a block only.**
5. **[!] No SDs, no n, no storage duration are printed for any of the 18 cells.**

---

## 3. ⚠️ TABLE S4 vs FIG. 1a — THE CORRUPTION IS **CONFIRMED**, AND THE TRUE VALUES ARE RECOVERED

`kang2026_extraction.md` §5a diagnosed, from digitisation alone, that **Fig. 1a's
"Sulfur-containing compounds" bars are wrong at 120 and 140 °C**, and hypothesised that **the
pH-experiment series had been pasted into the temperature panel**. Table S4 and Table S5 settle
both claims.

| 100 °C | 120 °C | 140 °C | source | |
|---|---|---|---|---|
| **15.2** | **14.3** | **19.4** | Fig. 1a "Sulfur-containing compounds" bar, digitised [D] | *what is plotted* |
| **13.978** | **35.866** | **60.400** | **Table S4 sulfur Subtotal [M]** | ***what is true*** |
| 14.868 | 13.978 | 19.157 | **Table S5** "Total content of sulfur-containing compounds", at pH 5.5 / 7 / 8 [M] | *what was pasted in* |

**⇒ VERDICT: the hypothesis is confirmed to within the digitisation error of ± 0.3 μg L⁻¹.**
The three numbers plotted in Fig. 1a's sulfur group (15.2 / 14.3 / 19.4) are the **pH ladder's
three sulfur totals** (14.868 / 13.978 / 19.157), **not** the temperature ladder's
(13.978 / 35.866 / 60.400). Agreement 15.2↔14.9, 14.3↔14.0, 19.4↔19.2 — inside the ±0.3
digitisation band on all three. **Fig. 1a's sulfur group is a mis-pasted pH series. The true
sulfur-class temperature ladder is 13.978 / 35.866 / 60.400 μg L⁻¹.**

**Fig. 1a is corrupted by a factor of 2.5× at 120 °C and 3.1× at 140 °C.** Any Eₐ or
temperature-response previously taken off that bar group is wrong by those factors and must be
recomputed from Table S4.

**Corroboration that Table S4 — not Fig. 1a — is the sound artefact, three independent ways:**

1. **Arithmetic closure [Z].** Every one of Table S4's nine subtotals/totals reconciles to
   ±0.003 from its own member rows (§4a). A pasted or shifted column would not close.
2. **The main paper's own Table 1 control row.** `kang2026_extraction.md` §4h records the
   no-amino-acid **Control at 120 °C / pH 7**: **MFT = 1.388, FFT = 4.107**, and the Control's
   O-heterocycle block total **6.158**. Table S4's **120 °C** column gives **MFT 1.388,
   FFT 4.107, O-subtotal 6.158** — **exact, to three decimals, on three independent
   quantities.** This is a hard external anchor from a different table in a different document.
   **Table S4's 120 °C column is the true pH-7 / 120 °C / 120 min condition.**
3. **The main paper's furfural sentence.** *"the content of furfural … reaching 11.039 μg L⁻¹
   at 140 °C, which was 3.3 times that at 100 °C."* Table S4: furfural **11.039** at 140 °C and
   **3.381** at 100 °C; **11.039 / 3.381 = 3.265 ≈ 3.3** ✔. **Table S4's 100 °C and 140 °C
   column labels are independently confirmed.**

**⇒ Table S4's temperature labels are correct and externally anchored at all three levels.**
The remaining ambiguity is in **Table S5's** stated temperature — see §4c.

---

## 4. TABLE S4 — **THE TEMPERATURE LADDER, COMPLETE** *(the artefact this wave was sent for)*

**Caption, raster-verified verbatim:** *"Table S4. Volatile flavor compounds (μg/L) identified
in the TTCA model reactions under elevated temperatures of 100, 120 and 140 ℃ at an initial pH
value of 7.0."*
**Column header, raster-verified:** `Compounds | RIᵃ | RIᵇ | Concentration (μg/L) [100℃ | 120℃ | 140℃] | IDᶜ`
**Conditions** (from `kang2026_extraction.md` §2b): **TTCA 10 mmol L⁻¹, initial pH 7.0, 120 min,
sealed pressure-rated glass vessel in an oil bath, ice-bath quench, HS-SPME-GC-MS, n = 3.**
**Footnotes, raster-verified:** *"Results were presented as means ± standard deviation, data
within a row with different letters are significantly different (p < 0.05) using Duncan's
multiple comparison test (n = 3). a: Linear retention indices calculated using a series of
n-alkanes on the DB-WAX column (30 m × 0.25 mm × 0.25 μm). b: Linear retention indices searched
from http://www.flavornet.org and http://webbook.nist.gov/chemistry/. c: ID: Identification
methods. d: "-", not detected."*

### 4a. Full transcription — **all 34 compound rows + 3 subtotals + total** [M]
Superscript Duncan letters retained verbatim. **T** = Tier (A = external standard curve in
Table S3; B = internal-standard-only semi-quant). `0.000` is printed where nothing was detected.

#### Sulfur-containing compounds (25 rows)
| compound | RIᵃ | RIᵇ | **100 °C** | **120 °C** | **140 °C** | ID | **T** |
|---|---|---|---|---|---|---|---|
| Thiophene | 1016 | 1022 | 0.000 | 0.000 | 0.118 ± 0.072 | RI, MS | A |
| 2-Methylthiophene | 1135 | 1095 | 0.000 | 0.000 | 0.892 ± 0.115 | RI, MS | A |
| Thiazole | 1240 | 1262 | 0.896 ± 0.053ᶜ | 1.796 ± 0.12ᵇ | 2.038 ± 0.139ᵃ | RI, MS | A |
| 2-Methylthiazole | 1272 | 1278 | 0.142 ± 0.008ᵇᶜ | 0.194 ± 0.014ᵇ | 0.254 ± 0.109ᵃ | RI, MS | **A⚠** |
| 3-Mercapto-2-butanone | 1273 | 1283 | 1.303 ± 0.153ᵇᶜ | 1.798 ± 0.127ᵇ | 2.398 ± 0.141ᵃ | RI, MS | A |
| 2,5-Dimethyl-thiazole | 1301 | 1326 | 0.079 ± 0.013ᶜ | 2.203 ± 0.142ᵇ | 2.553 ± 0.096ᵃ | RI, MS | A |
| **2-Methyl-3-furanthiol (MFT)** | 1302 | 1305 | **1.237 ± 0.142ᵇᶜ** | **1.388 ± 0.256ᵇ** | **5.907 ± 0.085ᵃ** | RI, MS | **A** |
| 2-Ethyl-thiazole | 1319 | 1304 | 0.245 ± 0.035ᵇ | 0.273 ± 0.006ᵃᵇ | 0.289 ± 0.025ᵃ | RI, MS | A |
| 3-Mercapto-2-pentanone | 1352 | 1343 | 0.973 ± 0.099ᶜ | 1.298 ± 0.253ᵇ | 1.736 ± 0.234ᵃ | RI, MS | A |
| 2,4,5-Trimethylthiazole | 1373 | 1390 | 0.419 ± 0.052ᶜ | 3.199 ± 0.075ᵇ | 3.843 ± 0.188ᵃ | RI, MS | A |
| **2-Furfurylthiol (FFT)** | 1426 | 1402 | **3.734 ± 0.085ᶜ** | **4.107 ± 0.137ᵇ** | **11.439 ± 0.265ᵃ** | RI, MS | **A** |
| 5-Ethyl-2,4-dimethyl-thiazole | 1437 | — | 0.000 | 0.587 ± 0.078ᵇ | 0.889 ± 0.085ᵃ | RI, MS | A |
| 2-Ethyl-4-methyl-thiazole | 1449 | 1410 | 0.117 ± 0.005ᶜ | 1.278 ± 0.117ᵇ | 2.112 ± 0.370ᵃ | RI, MS | A |
| 3-Thiophenethiol | 1564 | 1530 | 0.000 | 0.178 ± 0.027ᵇ | 0.268 ± 0.062ᵃ | RI, MS | A |
| 2-Acetylthiazole | 1646 | 1643 | 3.079 ± 0.215ᶜ | 8.795 ± 0.094ᵇ | 9.858 ± 1.215ᵃ | RI, MS | A |
| 2-Thiophenecarboxaldehyde | 1674 | 1679 | 0.09 ± 0.013ᶜ | 2.196 ± 0.25ᵇ | 2.28 ± 0.144ᵃ | RI, MS | A |
| 5-Methyl-2-thiophenecarboxaldehyde | 1701 | 1785 | 0.936 ± 0.06ᶜ | 3.649 ± 0.099ᵇ | 5.706 ± 0.644ᵃ | RI, MS | A |
| 2-Thiophenemethanethiol | 1702 | 1713 | 0.119 ± 0.074ᶜ | 0.196 ± 0.015ᵇ | 0.332 ± 0.21ᵃ | RI, MS | A |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio)furan | 1732 | — | 0.021 ± 0.004ᶜ | 0.098 ± 0.033ᵇ | 0.137 ± 0.036ᵃ | RI, MS | **B** |
| 1,2,3-Trithiolane | 1794 | — | 0.208 ± 0.029ᶜ | 0.378 ± 0.071ᵇ | 0.497 ± 0.185ᵃ | RI, MS | **B** |
| 3,3'-Dithiobisthiophene | 1845 | — | 0.000 | 0.000 | 0.000 | RI, MS | **B** |
| Thieno[3,2-b]thiophene | 1868 | 1843 | 0.000 | 1.302 ± 0.079ᵇ | 3.369 ± 0.454ᵃ | RI, MS | **B** |
| 2,5-Thiophenedicarboxaldehyde | 1911 | — | 0.346 ± 0.039ᶜ | 0.724 ± 0.139ᵇ | 1.283 ± 0.122ᵃ | RI, MS | A |
| 2-Methylthieno[2,3-b]thiophene | 1947 | — | 0.000 | 0.19 ± 0.03ᵇ | 2.123 ± 0.221ᵃ | RI, MS | **B** |
| Bis(2-furfuryl)sulfide | 2419 | — | 0.034 ± 0.007ᵇᶜ | 0.042 ± 0.004ᵇ | 0.078 ± 0.018ᵃ | RI, MS | **B** |
| **Subtotal** | | | **13.978 ± 0.568ᶜ** | **35.866 ± 1.382ᵇ** | **60.400 ± 1.620ᵃ** | | |

#### Nitrogen-containing heterocycles (4 rows) — **all zero at every temperature**
| compound | RIᵃ | RIᵇ | 100 °C | 120 °C | 140 °C | ID | T |
|---|---|---|---|---|---|---|---|
| Methylpyrazine | 1214 | 1263 | 0.000 | 0.000 | 0.000 | RI, MS | A |
| Pyrazine | 1216 | 1209 | 0.000 | 0.000 | 0.000 | RI, MS | A |
| 3-Methyl-pyridine | 1307 | 1346 | 0.000 | 0.000 | 0.000 | RI, MS | A |
| 2,5-Dimethyl-pyrazine | 1320 | 1328 | 0.000 | 0.000 | 0.000 | RI, MS | A |
| **Subtotal** | | | **0.000** | **0.000** | **0.000** | | |

#### Oxygen-containing heterocycles (5 rows)
| compound | RIᵃ | RIᵇ | 100 °C | 120 °C | 140 °C | ID | T |
|---|---|---|---|---|---|---|---|
| Furan | 797 | 798 | 0.000 | 0.138 ± 0.017ᵃ | 1.108 ± 0.092ᵇ | RI, MS | A |
| 2-Methyl-furan | 851 | 829 | 0.000 | 0.019 ± 0.007ᵇ | 0.133 ± 0.088ᵃ | RI, MS | A |
| **Furfural** | 1457 | 1460 | **3.381 ± 0.089ᶜ** | **5.793 ± 0.422ᵇ** | **11.039 ± 0.302ᵃ** | RI, MS | **A** |
| 1-(2-Furanyl)-ethanone | 1497 | 1501 | 0.000 | 0.136 ± 0.017ᵇ | 0.379 ± 0.087ᵃ | RI, MS | A |
| 2(5H)-Furanone | 1748 | 1767 | 0.000 | 0.072 ± 0.012ᵇ | 0.098 ± 0.015ᵃ | RI, MS | A |
| **Subtotal** | | | **3.381 ± 0.089ᶜ** | **6.158 ± 0.425ᵇ** | **12.757 ± 0.469ᵃ** | RI, MS | |
| **Total** | | | **17.359 ± 0.527ᶜ** | **42.024 ± 1.737ᵇ** | **73.157 ± 1.557ᵃ** | RI, MS | |

**[!] Note on the Duncan letters in the Furan row:** every other row runs `c < b < a` with
increasing temperature; the Furan row is printed **`0.000 | 0.138ᵃ | 1.108ᵇ`** — the letters
are inverted relative to the whole rest of the table. **Typographic error in the letters only;
the concentrations are consistent with every other row's monotone increase.**

### 4b. **INTERNAL-CONSISTENCY AUDIT — Table S4 passes completely** [Z]
Every subtotal and total recomputed from the transcribed member rows:

| quantity | recomputed | printed | Δ |
|---|---|---|---|
| Sulfur subtotal, 100 °C | **13.978** | 13.978 | **0.000** |
| Sulfur subtotal, 120 °C | **35.869** | 35.866 | 0.003 (rounding) |
| Sulfur subtotal, 140 °C | **60.399** | 60.400 | 0.001 (rounding) |
| O-heterocycle subtotal, 100/120/140 | **3.381 / 6.158 / 12.757** | same | **0.000** |
| N-heterocycle subtotal | **0.000** ×3 | same | **0.000** |
| Total, 100 °C | 13.978 + 0 + 3.381 = **17.359** | 17.359 | **0.000** |
| Total, 120 °C | 35.866 + 0 + 6.158 = **42.024** | 42.024 | **0.000** |
| Total, 140 °C | 60.400 + 0 + 12.757 = **73.157** | 73.157 | **0.000** |

**⇒ Table S4 is arithmetically sound in all nine aggregates.** Combined with the three external
anchors of §3, **Table S4 is the most internally-and-externally validated numeric artefact this
wave has handled.**

### 4c. **⚠️ [!] THE ONE SERIOUS INTEGRITY PROBLEM: Table S5's pH-7 column IS Table S4's 100 °C column**

Table S5 is captioned *"…at the temperature of **120 ℃**"* and Table S4's pH-7 condition is the
whole of Table S4. The two tables therefore share exactly one condition, and it should appear
identically in both. **It does not appear identically — it appears in the wrong column.**

| quantity | **S4 100 °C** | **S4 120 °C** | **S5 "pH 7, 120 °C"** | matches |
|---|---|---|---|---|
| MFT | **1.237** ± 0.142 | 1.388 ± 0.256 | **1.237** ± 0.061 | **S4 100 °C** |
| FFT | **3.734** ± 0.085 | 4.107 ± 0.137 | **3.734** ± 0.074 | **S4 100 °C** |
| 3-Mercapto-2-butanone | **1.303** ± 0.153 | 1.798 ± 0.127 | **1.303** ± 0.006 | **S4 100 °C** |
| 3-Mercapto-2-pentanone | **0.973** ± 0.099 | 1.298 ± 0.253 | **0.973** ± 0.025 | **S4 100 °C** |
| Thiazole | **0.896** ± 0.053 | 1.796 ± 0.12 | **0.896** ± 0.012 | **S4 100 °C** |
| 2-Acetylthiazole | **3.079** ± 0.215 | 8.795 ± 0.094 | **3.079** ± 0.075 | **S4 100 °C** |
| Furfural | **3.381** ± 0.089 | 5.793 ± 0.422 | **3.381** ± 0.34 | **S4 100 °C** |
| Sulfur total | **13.978** ± 0.568 | 35.866 ± 1.382 | **13.978** ± 0.081 | **S4 100 °C** |
| Grand total | **17.359** ± 0.527 | 42.024 ± 1.737 | **17.359** ± 0.421 | **S4 100 °C** |

**Every shared mean matches S4's 100 °C column to the last printed decimal — ~20 compounds
plus 4 aggregates. That is not coincidence. But not one of the SDs matches.**

**Two live hypotheses.** Neither is fully resolvable from the published record:

- **H-A (~65 %): the pH ladder was actually run at 100 °C, and both the S5 caption and the
  main paper's §2.4 / Fig. 1b caption ("120 °C") are wrong.** Parsimonious — one label error,
  no data corrupt. Strongly supported by magnitude: S5's other two columns are **20.175**
  (pH 5.5) and **21.181** (pH 8), both close to S4's **100 °C** total of 17.359 and nowhere
  near its 120 °C total of 42.024. Under the alternative, the pH-7 point would be a 2× spike
  above both its neighbours at the same temperature — chemically implausible.
- **H-B (~25 %): the pH ladder was run at 120 °C, and S5's pH-7 column was filled with the
  wrong column of S4.** Then S5's true pH-7 numbers are S4's 120 °C column, **and Fig. 1b's
  pH-7 bars are also wrong** (they plot S5 as printed — verified in
  `kang2026_extraction.md` §5b).
- (~10 % something else.)

**⇒ RULING FOR THE REPO.**
1. **Table S4 is unaffected and remains fully trustworthy** — it is anchored three independent
   ways (§3) including an exact three-quantity match to the main paper's own Table 1 control.
   **Use Table S4 without reservation.**
2. **Table S5's absolute temperature is UNSAFE. Its pH-*ratios* are safe** (within-table,
   same-analysis), and are what a pH-response model needs anyway. **Register Table S5 as
   `temperature: 100_or_120_C_UNRESOLVED` and use ratios only.**
3. **⚠️ THE SDs IN BOTH TABLES ARE NOT TRUSTWORTHY.** The same n = 3 experiment is reported
   twice with identical means and *different* standard deviations (0.006 vs 0.153; 0.012 vs
   0.053; 0.34 vs 0.089; 0.081 vs 0.568). At least one set was regenerated rather than
   measured. **Do not use Table S4 or S5 SDs as measurement uncertainty, and do not weight any
   fit by them.** Assign uncertainty from the method instead (see §7).

### 4d. Tier A / Tier B split, and **what fraction of the ladder is Tier A** [Z]
Cross-referencing all 34 S4 rows against Table S3's 31 calibrated compounds:

- **Tier A (28 rows):** all 19 sulfur rows listed A above, all 4 N rows, all 5 O rows.
- **Tier B (6 rows, no standard curve):** 2-methyl-3-[(2-methyl-3-thienyl)dithio)furan;
  1,2,3-trithiolane; 3,3'-dithiobisthiophene; thieno[3,2-b]thiophene;
  2-methylthieno[2,3-b]thiophene; bis(2-furfuryl)sulfide. *(Unsurprising — these are the
  polysulfides and fused thiophenes for which reference standards are not commercial.)*
- **Excluded (1 row):** 2-methylthiazole — nominally Tier A but its curve is impossible (§2c).

| | 100 °C | 120 °C | 140 °C |
|---|---|---|---|
| Tier B contribution to sulfur subtotal | **0.263** | **2.010** | **6.204** |
| **as % of sulfur subtotal** | **1.9 %** | **5.6 %** | **10.3 %** |
| **⇒ sulfur subtotal is Tier A by** | **98.1 %** | **94.4 %** | **89.7 %** |
| Tier B contribution to O-subtotal & N-subtotal | 0 | 0 | 0 |

**⇒ Unlike the main paper's class totals (which `kang2026_extraction.md` §1 correctly warned
mix tiers), Table S4's sulfur subtotal is ≥ 89.7 % Tier A at every temperature, and its
O-heterocycle subtotal is 100 % Tier A. The class totals in Table S4 are usable as
quantitative — a materially better status than the main paper allowed.** Tier B's share does
grow monotonically with temperature (1.9 → 10.3 %), so the *shape* of the sulfur-class ladder
carries a small tier-mixing bias that flattens it slightly; the Tier-A-only sulfur subtotals
are **13.715 / 33.856 / 54.196**, which changes the 100→140 fold-change from 4.32× to 3.95×.

---

## 5. TABLE S5 — the pH ladder, complete

**Caption, raster-verified verbatim:** *"Table S5. Volatile flavor compounds (μg/L) identified
in the TTCA model reactions under different pH values of 5.5, 7 and 8 at the temperature of
120 ℃."* **⚠️ the stated temperature is unsafe — see §4c.**
**Conditions:** TTCA 10 mmol L⁻¹, pH adjusted with NaOH 6 mol L⁻¹, 120 min. **Measured pH drift
(main paper p. 3242): initial 5.5 / 7 / 8 → final 5.1 / 4.9 / 4.7. All three runs finish
acidic; the pH label is an *initial* pH only.**
Note S5 sub-classes the sulfur block (Thiols/Sulfides · Thiophenes · Thiazoles · Other), which
S4 does not — this is the sub-classing Fig. 1a/1b plot.

### 5a. Full transcription [M]
| compound | RIᵃ | RIᵇ | **pH 5.5** | **pH 7** | **pH 8** | ID |
|---|---|---|---|---|---|---|
| **Thiols/Sulfides** | | | | | | |
| 3-Mercapto-2-butanone | 1273 | 1283 | 1.796 ± 0.152ᵃ | 1.303 ± 0.006ᵇ | 0.973 ± 0.077ᶜ | RI, MS |
| **2-Methyl-3-furanthiol (MFT)** | 1302 | 1305 | **1.936 ± 0.065ᵃ** | **1.237 ± 0.061ᵇ** | **1.103 ± 0.008ᶜ** | RI, MS |
| 3-Mercapto-2-pentanone | 1352 | 1343 | 1.793 ± 0.036ᵃ | 0.973 ± 0.025ᶜ | 1.332 ± 0.176ᵇ | RI, MS |
| **2-Furfurylthiol (FFT)** | 1426 | 1402 | **5.639 ± 0.095ᵃ** | **3.734 ± 0.074ᵇ** | **3.012 ± 0.097ᶜ** | RI, MS |
| 3-Thiophenethiol | 1564 | 1530 | 0.000 | 0.000 | 0.079 ± 0.004ᵃ | RI, MS |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio)furan | 1732 | — | 0.179 ± 0.014ᵃ | 0.021 ± 0.004ᵇ | 0.000 | MS |
| 2-Thiophenemethanethiol | **1947 ⚠** | — | 0.211 ± 0.007ᵃ | 0.119 ± 0.016ᵇ | 0.034 ± 0.004ᶜ | MS |
| Bis(2-furfuryl)sulfide | 2419 | — | 0.057 ± 0.007ᵃ | 0.034 ± 0.012ᵇ | 0.000 | MS |
| **Subtotal** | | | **11.611 ± 0.123ᵃ** | **7.421 ± 0.146ᵇ** | **6.533 ± 0.143ᶜ** | |
| **Thiophenes** | | | | | | |
| Thiophene | 1016 | 1022 | 0.000 | 0.000 | 0.179 ± 0.008ᵃ | RI, MS |
| 2,3-Dihydro-5-methyl-thiophene | 1133 | 1156 | 0.000 | 0.000 | 0.000 | RI, MS |
| 2-Methylthiophene | 1135 | 1095 | 0.028 ± 0.004ᵃ | 0.000 | 0.000 | RI, MS |
| 2,5-Dimethyl-thiophene | 1144 | 1190 | 0.000 | 0.000 | 0.193 ± 0.014ᵃ | RI, MS |
| Dihydro-2-methyl-3(2H)-thiophenone | 1497 | 1506 | 0.000 | 0.000 | 0.079 ± 0.007ᵃ | RI, MS |
| Dihydro-3-(2H)-thiophenone | 1542 | 1563 | 0.000 | 0.000 | 0.087 ± 0.023ᵃ | RI, MS |
| 2-Thiophenecarboxaldehyde | 1674 | 1679 | 0.221 ± 0.005ᵃ | 0.09 ± 0.007ᵇ | 0.012 ± 0.004ᶜ | RI, MS |
| Thieno[3,2-b]thiophene | 1868 | 1843 | 0.000 | 0.000 | 0.196 ± 0.011ᵃ | RI, MS |
| 5-Methyl-2-thiophenecarboxaldehyde | 1701 | 1785 | 0.837 ± 0.01ᵇ | 0.936 ± 0.043ᵃ | 0.423 ± 0.016ᶜ | RI, MS |
| 2-Acetyl-3-methylthiophene | 1757 | — | 0.000 | 0.000 | 0.000 | MS |
| 3,3'-Dithiobisthiophene | 1845 | — | 0.000 | 0.000 | 0.000 | MS |
| 2,5-Thiophenedicarboxaldehyde | 1911 | — | 0.336 ± 0.023ᵃ | 0.346 ± 0.027ᵃ | 0.329 ± 0.011ᵃ | MS |
| 2-Methylthieno[2,3-b]thiophene | 1947 | — | 0.000 | 0.000 | 0.189 ± 0.02ᵃ | RI, |
| 2,3'-Bithiophene | 2194 | — | 0.000 | 0.000 | 0.023 ± 0.007ᵃ | MS |
| **Subtotal** | | | **1.422 ± 0.039ᵇ** | **1.372 ± 0.068ᵇ** | **1.71 ± 0.025ᵃ** | |
| **Thiazoles** | | | | | | |
| Thiazole | 1240 | 1262 | 0.334 ± 0.012ᶜ | 0.896 ± 0.012ᵇ | 2.137 ± 0.054ᵃ | RI, MS |
| 2-Methylthiazole | 1272 | 1278 | 0.023 ± 0.004ᵇ | 0.142 ± 0.007ᵇ | 1.779 ± 0.113ᵃ | RI, MS |
| 2,5-Dimethyl-thiazole | 1301 | 1326 | 0.000 | 0.079 ± 0.015ᵇ | 0.127 ± 0.008ᵃ | RI, MS |
| 2,4,5-Trimethylthiazole | 1373 | 1390 | 0.112 ± 0.067ᶜ | 0.419 ± 0.074ᵇ | 1.192 ± 0.021ᵃ | RI, MS |
| 2-Ethyl-thiazole | 1319 | 1304 | 0.000 | 0.245 ± 0.012ᵇ | 0.283 ± 0.018ᵃ | RI, MS |
| 2-Ethyl-4-methyl-thiazole | 1449 | 1410 | 0.000 | 0.117 ± 0.01ᵇ | 0.217 ± 0.007ᵃ | RI, MS |
| 2-Acetylthiazole | 1646 | 1643 | 0.238 ± 0.008ᶜ | 3.079 ± 0.075ᵇ | 5.179 ± 0.085ᵃ | RI, MS |
| **Subtotal** | | | **0.707 ± 0.07ᶜ** | **4.977 ± 0.175ᵇ** | **10.914 ± 0.226ᵃ** | |
| **Other sulfur-containing compounds** | | | | | | |
| 1,2,3-Trithiolane | 1794 | — | 1.128 ± 0.01ᶜ | 0.208 ± 0.009ᵇ | 0.000 | MS |
| **Nitrogen-containing heterocycles** | | | | | | |
| Pyrazine | 1210 | 1209 | 0.000 | 0.000 | 0.129 ± 0.007ᵃ | RI, MS |
| Methylpyrazine | 1214 | 1263 | 0.000 | 0.000 | 0.271 ± 0.016ᵃ | RI, MS |
| 3-Methyl-pyridine | 1307 | 1346 | 0.000 | 0.000 | 0.037 ± 0.013ᵃ | RI, MS |
| 2,5-Dimethyl-pyrazine | 1320 | 1328 | 0.000 | 0.000 | 0.438 ± 0.033ᵃ | RI, MS |
| Pyrrole | 1497 | 1507 | 0.000 | 0.000 | 0.021 ± 0.002ᵃ | RI, MS |
| **Subtotal** | | | **0.000** | **0.000** | **0.896 ± 0.032ᵃ** | |
| **Oxygen-containing heterocycles** | | | | | | |
| Furan | 797 | 798 | 0.796 ± 0.115ᵃ | 0.000 | 0.000 | RI, MS |
| 2-Melthy-furan *(sic)* | 851 | 829 | 0.273 ± 0.045ᵃ | 0.000 | 0.000 | RI, MS |
| **Furfural** | 1457 | 1460 | **4.238 ± 0.32ᵃ** | **3.381 ± 0.34ᵇ** | **1.128 ± 0.147ᶜ** | RI, MS |
| Benzofuran | 1493 | — | 0.000 | 0.000 | 0.000 | MS |
| 4-Hydroxy-5-methyl-3(2H)-furanone | 2108 | 2124 | 0.000 | 0.000 | 0.000 | RI, MS |
| **Subtotal** | | | **5.307 ± 0.459ᵃ** | **3.381 ± 0.34ᵇ** | **1.128 ± 0.147ᶜ** | |
| **Total content of sulfur-containing compounds** | | | **14.868 ± 0.096ᵇ** | **13.978 ± 0.081ᶜ** | **19.157 ± 0.353ᵃ** | |
| **Total** | | | **20.175 ± 0.504ᵃ** | **17.359 ± 0.421ᵇ** | **21.181 ± 0.334ᵃ** | |

**[!] Transcription errors in S5, as printed:**
- **2-Thiophenemethanethiol is given RIᵃ = 1947**; Table S4 gives it **1702**, and 1947 is
  2-methylthieno[2,3-b]thiophene's RI (which is also printed as 1947 four rows below).
  **S5's RI is wrong; S4's 1702 is right.**
- **Pyrazine RIᵃ = 1210** in S5 vs **1216** in S4. Minor.
- **"2-Melthy-furan"** — misspelling of 2-methyl-furan.
- **"2-Methyl-3-[(2-methyl-3-thienyl)dithio)furan"** — unbalanced brackets, both tables.
- The ID column for 2-methylthieno[2,3-b]thiophene reads **"RI,"** (truncated).

### 5b. **Internal-consistency audit — Table S5 also passes** [Z]
All nine subtotals + both grand totals recomputed from the member rows: **every one reconciles
to ±0.001.** (Thiols 11.611/7.421/6.533 ✔ · Thiophenes 1.422/1.372/1.710 ✔ · Thiazoles
0.707/4.977/10.914 ✔ · Sulfur total 14.868/13.978/19.157 ✔ · Grand total
20.175/17.359/21.181 ✔.) **The arithmetic is sound; only the column *label* is in doubt.**

### 5c. **The pH response** [Z] — *safe as ratios, at whichever temperature*
| quantity | pH 5.5 | pH 7 | pH 8 | **pH 5.5 → 8** |
|---|---|---|---|---|
| **MFT** | 1.936 | 1.237 | 1.103 | **× 0.57 (−43 %)** |
| **FFT** | 5.639 | 3.734 | 3.012 | **× 0.53 (−47 %)** |
| **FFT / MFT ratio** | **2.91** | **3.02** | **2.73** | **flat (±5 %)** |
| Thiols/Sulfides class | 11.611 | 7.421 | 6.533 | **× 0.56** |
| Thiophenes class | 1.422 | 1.372 | 1.710 | × 1.20 (flat) |
| **Thiazoles class** | 0.707 | 4.977 | 10.914 | **× 15.4** |
| N-heterocycles | 0.000 | 0.000 | 0.896 | **pH 8 only** |
| **O-heterocycles (furans)** | 5.307 | 3.381 | 1.128 | **× 0.21** |
| Total sulfur | 14.868 | 13.978 | 19.157 | × 1.29 |
| Grand total | 20.175 | 17.359 | 21.181 | × 1.05 (flat) |

**The three structural results, all monotone and all ≥ 1.8× — usable as directional gates:**
1. **Free thiols (MFT, FFT, the mercaptoketones) fall ~1.8× from pH 5.5 to 8**, and the
   thiol class falls 1.78×. **Acid favours thiols.**
2. **Thiazoles rise 15.4×** over the same span. **Base favours N,S-heterocyclisation.**
   2-Acetylthiazole alone rises **21.8×** (0.238 → 5.179).
3. **Furans fall 4.7×**, and pyrazines/pyridine/pyrrole appear **only at pH 8**.
4. **⭐ The FFT/MFT ratio is invariant to pH (2.91 / 3.02 / 2.73, ±5 %)** — a striking contrast
   with `kang2026_extraction.md` §4h, where the same ratio spans **0.15–9.2** across
   amino-acid co-substrates at fixed pH. **⇒ The FFT:MFT branching ratio is set by the
   nitrogen co-substrate, NOT by pH.** This is a clean, directly testable structural claim, and
   the strongest single result in Table S5.
5. **⇒ The grand total is nearly flat across pH (20.2 / 17.4 / 21.2) while its composition
   turns over completely.** pH is a *selectivity* lever here, not a *yield* lever. Contrast
   temperature (§6), which is a yield lever (4.2×) *and* a selectivity lever.

---

## 6. FIGURES S1–S4 — **and the H₂S question, answered**

### 6a. ⚠️ **THERE IS NO H₂S EVIDENCE IN THIS SI. NO H₂S IS MEASURED ANYWHERE.**
The wave brief expected "Figs S2/S3 (H₂S evidence)". **The SI's actual figure list is:**

| SI figure | caption as printed | content |
|---|---|---|
| **Fig. S1** | *"The high performance liquid chromatogram of TTCA/ARP compounds."* | HPLC-UV chromatogram, **no axis values legible, no peak table**. Purity is *not* quantified. |
| **Fig. S2** | *"The total ions chromatogram (a) / LC-MS/MS spectra (b) of purified TTCA."* | TIC + MS/MS, panels (a) and (b). **No fragment m/z are labelled in the raster.** Table S1, which was supposed to tabulate them, is empty (§2a). |
| **Fig. S3** | *"The release of Cys during the degradation of TTCA."* | **Cys concentration (mmol/L) vs time (0–120 min), 3 curves: TTCA-100 °C, TTCA-120 °C, TTCA-140 °C.** ⭐ |
| **Fig. S4** | *"Depletion of Cys compound at an initial pH of 7 under different temperatures (100, 120 and 140 ℃)."* | **Cys concentration (mmol/L) vs time (0–120 min), 3 curves: 100 °C, 120 °C, 140 °C, with error bars.** ⭐ |

**[!] The main paper's SI figure numbering is off by one.** Main text p. 3242 reads *"thermal
degradation of TTCA and free Cys (**Fig. S2 and S3**), resulting in an increased release of
H₂S."* The figures that actually show TTCA degradation and free-Cys depletion are the SI's
**Fig. S3 and Fig. S4**. Cite by caption, not by number.

**⇒ The H₂S claim in the main paper is inferential, not measured.** No H₂S concentration, no
H₂S detector, no sulfide assay appears in the main paper or the SI. The evidence chain is:
Cys is released from TTCA (Fig. S3) → free Cys is depleted faster at higher T (Fig. S4) →
*therefore* H₂S release rises. **Register any H₂S-related parameter sourced from Kang 2026 as
`[Q] inferred, not measured`.**

### 6b. ⭐ **Fig. S4 digitised — free-Cys depletion, and an activation energy** **[D][Z]**
**Digitisation method:** 400 dpi raster; y-axis and x-axis ruling lines located by ink-column /
ink-row maxima; tick centroids extracted programmatically (y ticks at pixel rows 553 / 732 /
911 / 1088 / 1267 for **10 / 8 / 6 / 4 / 2 mmol L⁻¹**, spacing uniform to ±0.5 px; x ticks at
1228 / 1386.5 / 1543.5 / 1699.5 / 1856.5 / 2014.5 / 2171.5 for **0…120 min**, spacing uniform
to ±1 px). Marker centroids read from ±7 px column bands. **Residual digitisation error
≈ ± 0.03 mmol L⁻¹** (better than one marker radius). Axis calibration is exact by construction;
all three curves start at exactly 10.00 at t = 0, which is an independent validation.

**System:** *free* Cys, **10 mmol L⁻¹ initial**, initial **pH 7**, 0–120 min.
*(Note: the caption says only "Cys compound". Because it starts at exactly 10.00 mmol L⁻¹ —
the same loading as the TTCA systems — and because Fig. S3 separately shows Cys **rising from
zero** in the TTCA system, Fig. S4 is read here as a **separate free-Cys control system**,
which is also how the main text uses it ("thermal degradation of TTCA **and free Cys**").
An alternative reading — that Fig. S4 tracks the Cys *moiety* still bound in TTCA — cannot be
excluded from the published record. **Flag as `system_identity: free_Cys (85 % confidence)`.**)*

| t (min) | **100 °C** | **120 °C** | **140 °C** |
|---|---|---|---|
| 0 | **10.00** | **10.00** | **10.00** |
| 20 | **9.77** | **8.63** | **7.71** |
| 40 | **9.61** | **7.81** | **6.24** |
| 60 | **9.41** | **7.11** | **5.49** |
| 80 | **9.12** | **6.90** | **4.76** |
| 100 | **8.64** | **6.55** | **3.99** |
| 120 | **8.38** | **6.13** | **3.74** |
| **conversion at 120 min** | **16.2 %** | **38.7 %** | **62.6 %** |

**⇒ Effective first-order rate constants** (from the 120-min endpoint, k = −ln(C/C₀)/t) **[Z]:**

| T | 1/T (K⁻¹) | **k_eff (min⁻¹)** | **k_eff (s⁻¹)** |
|---|---|---|---|
| 100 °C (373.15 K) | 2.6799 × 10⁻³ | **1.472 × 10⁻³** | 2.45 × 10⁻⁵ |
| 120 °C (393.15 K) | 2.5436 × 10⁻³ | **4.078 × 10⁻³** | 6.80 × 10⁻⁵ |
| 140 °C (413.15 K) | 2.4204 × 10⁻³ | **8.195 × 10⁻³** | 1.37 × 10⁻⁴ |

**⭐ Arrhenius fit over all three points [Z]:**

> **Eₐ (free-Cys depletion, pH 7, aqueous, sealed) = 55.1 kJ mol⁻¹**
> **A = 8.0 × 10⁴ min⁻¹ = 1.3 × 10³ s⁻¹**
> **R² = 0.994**, slope = −6631.7 K, residuals ≤ 0.077 in ln k
> **Q₁₀ ≈ 1.55**; k(140)/k(100) = **5.57×** over 40 K

**Caveats [!]:** (i) the decays are not cleanly first-order — 120 °C and 140 °C *decelerate*
(k over 0–20 min is 7.4 and 13.0 × 10⁻³ min⁻¹, roughly 1.6–1.8× the 120-min effective value),
so the endpoint k is an average over a non-exponential curve; (ii) 100 °C mildly *accelerates*
(1.01 × 10⁻³ at 60 min → 1.47 × 10⁻³ at 120 min). Fitting instead on the **0–20 min initial
rates** gives Eₐ = **51 kJ mol⁻¹** — reassuringly close, so **the 55 kJ mol⁻¹ figure is robust
to the fitting window at the ±10 % level**. (iii) Error bars are drawn on Fig. S4 and are small
(≈ ±0.15 mmol L⁻¹ at 100 min, largest point); n is not stated for this figure.

**This is a genuinely new number for the corpus: the first measured activation energy for
cysteine consumption under Maillard conditions, from a three-temperature ladder, with R² 0.994.**

### 6c. ⭐ **Fig. S3 digitised — free Cys *released* from TTCA** **[D]**
Same digitisation method (y ticks at 702 / 844.5 / 986 / 1129.5 / 1271.5 / 1414.5 px for
**2.5 / 2.0 / 1.5 / 1.0 / 0.5 / 0.0 mmol L⁻¹**; x ticks 1299…2201 px for 0…120 min).
Residual error **≈ ± 0.03 mmol L⁻¹**. **No error bars are drawn on Fig. S3.**
**System:** TTCA **10 mmol L⁻¹**, pH 7, 0–120 min, at 100 / 120 / 140 °C.

| t (min) | **TTCA-100 °C** | **TTCA-120 °C** | **TTCA-140 °C** |
|---|---|---|---|
| 0 | 0.00 | 0.00 | 0.00 |
| 20 | **0.121** | **0.205** | **0.500** |
| 40 | **0.242** | **0.989** | **1.630** |
| 60 | **0.693** | **1.219** | **1.574** |
| 80 | **0.872** | **1.112** | **1.619** |
| 100 | **≈1.30** ⚠ | **0.988** | **≈1.30** ⚠ |
| 120 | **1.233** | **0.874** | **1.088** |

⚠ At t = 100 min the 100 °C and 140 °C markers cross and merge into a single 21-px blob at
**1.296**; the digitiser cannot separate them. Both are ≈ 1.29–1.31.

**[Z] What Fig. S3 establishes — three results the repo can use:**
1. **Classic A → B → C intermediate kinetics.** Free Cys **rises to a maximum then falls** at
   every temperature. **t_max ≈ 40 min at 140 °C, ≈ 55–60 min at 120 °C, ≈ 100 min at 100 °C** —
   the maximum moves earlier as T rises, exactly as a consecutive scheme requires. **Cys is a
   transient intermediate of TTCA degradation, not a product.** Any model treating released Cys
   as an accumulating pool is wrong.
2. **⭐ The TTCA → free-Cys yield is capped at ≈ 16 mol %.** Peak free Cys is **1.63 mmol L⁻¹**
   (140 °C, 40 min) against **10 mmol L⁻¹ TTCA** loaded. Even at the most favourable
   temperature and time, **≤ 16.3 % of the cysteine moiety ever appears as free Cys.** The
   other ≥ 84 % goes to volatiles, non-volatile browning polymer, or degrades without
   liberating free Cys. **This is a hard mass-balance bound on the sulfur branch's precursor
   flux, and the corpus has had nothing like it.**
3. **Higher temperature raises both the rate and the peak height**, but the 120 °C peak
   (1.22) is *below* the 100 °C late-time value (1.23–1.30). **The three curves cross.**
   Consequently **free-Cys concentration is NOT monotone in temperature at fixed time** — at
   120 min the ordering is 100 °C (1.233) > 140 °C (1.088) > 120 °C (0.874). **Any model that
   assumes monotone temperature response for the Cys pool will be wrong at long times.**
4. **Consistency with Fig. S4 [Z]:** at 140 °C, TTCA releases ≈ 1.6 mmol L⁻¹ free Cys by
   40 min while free Cys is itself being consumed with k ≈ 8.2 × 10⁻³ min⁻¹ (§6b). The
   steady-state pool implied by (release rate)/(k) is of the right order (~1–2 mmol L⁻¹),
   and the turnover explains the post-peak decline without invoking anything else. **The two
   figures are mutually consistent.**

---

## 7. ⭐ THE TEMPERATURE LADDER ANALYSED — **the deliverable**

### 7a. **MFT and FFT are strongly non-Arrhenius: a threshold sits between 120 and 140 °C** [Z]
All values are yields at **120 min**, TTCA 10 mmol L⁻¹, pH 7.

| species | 100 °C | 120 °C | 140 °C | **×(100→120)** | **×(120→140)** | **Eₐ 100→120** | **Eₐ 120→140** | **Eₐ 100→140** |
|---|---|---|---|---|---|---|---|---|
| **MFT** | 1.237 | 1.388 | **5.907** | **1.12×** | **4.26×** | **7.0 kJ mol⁻¹** | **97.8 kJ mol⁻¹** | 50.1 |
| **FFT** | 3.734 | 4.107 | **11.439** | **1.10×** | **2.78×** | **5.8 kJ mol⁻¹** | **69.2 kJ mol⁻¹** | 35.9 |
| Furfural | 3.381 | 5.793 | 11.039 | 1.71× | 1.91× | 32.8 | 43.5 | **37.9** |
| 2-Acetylthiazole | 3.079 | 8.795 | 9.858 | **2.86×** | 1.12× | **64.0** | **7.7** | 39.7 |
| **Sulfur class** | 13.978 | 35.866 | 60.400 | 2.57× | 1.68× | 57.5 | 35.2 | **46.9** |
| O-heterocycle class | 3.381 | 6.158 | 12.757 | 1.82× | 2.07× | 36.4 | 48.6 | **42.3** |
| **Grand total** | 17.359 | 42.024 | 73.157 | 2.42× | 1.74× | 53.6 | 37.3 | **46.1** |

**⭐ THE HEADLINE MECHANISTIC RESULT.** The two **flavour-critical free thiols behave in the
opposite direction to their own class.** The sulfur class as a whole is **decelerating**
(2.57× then 1.68×) — the normal signature of a reaction approaching precursor exhaustion, and
consistent with the class's Eₐ falling from 57.5 to 35.2 kJ mol⁻¹. **MFT and FFT are
accelerating** (1.12× then 4.26×; 1.10× then 2.78×), with apparent Eₐ climbing from ~6–7 to
70–98 kJ mol⁻¹. **A saturation artefact cannot produce this** — precursor depletion *depresses*
apparent Eₐ at the top of the ladder, which is exactly what the class does and the opposite of
what the thiols do. **⇒ MFT and FFT formation switches on between 120 and 140 °C, over and
above the sulfur branch as a whole.**

Mechanistically this coheres with §6b–c: below 140 °C the sulfur output is dominated by
**heterocyclisation** (2-acetylthiazole ×2.86 then ×1.12, i.e. essentially complete by 120 °C;
thiazoles and thiophenes carry the class), while **free-thiol formation requires the H₂S /
Cys-degradation flux that only becomes large at 140 °C** (Fig. S4: 62.6 % Cys conversion at
140 °C vs 38.7 % at 120 °C). **The corpus now has a mechanistic reason to expect a
temperature threshold in the MFT/FFT branch, and a measured ladder that exhibits it.**

**⚠️ Practical consequence: a single Arrhenius Eₐ cannot reproduce the MFT or FFT ladder.**
A 3-point Arrhenius fit to MFT gives Eₐ = 49.3 kJ mol⁻¹ at **R² = 0.781** — a fit that would
mispredict the 120 °C point by −38 % and the 140 °C point by +28 %. Furfural, by contrast,
fits cleanly (32.8 → 43.5 kJ mol⁻¹ across the two legs, ~15 % spread). **Do not fit MFT/FFT
with a single Eₐ and report R²; either fit a two-regime / threshold form, or fit the
100–120 leg and hold out 140 °C (recommended — see §9).**

### 7b. **The FFT/MFT branching ratio vs temperature** [Z]
| | 100 °C | 120 °C | 140 °C | | pH 5.5 | pH 7 | pH 8 |
|---|---|---|---|---|---|---|---|
| **FFT / MFT** | **3.02** | **2.96** | **1.94** | | **2.91** | **3.02** | **2.73** |

**⇒ Combining this dossier with `kang2026_extraction.md` §4h, the FFT:MFT ratio is:**
- **invariant to pH** (2.73–3.02, ±5 %, §5c),
- **weakly sensitive to temperature** (3.02 → 1.94, a 1.6× swing, all of it in the 120→140 leg),
- **wildly sensitive to the nitrogen co-substrate** (**0.15 – 9.23**, a **62× span**, §4h of the
  main dossier).

**⇒ If the repo carries a fixed FFT:MFT branching ratio, ~3.0 is the right value for a
Cys/pentose system with no added amino acid at ≤ 120 °C, and it is robust to pH — but it is
NOT robust to co-substrate, and it drops toward 2 at 140 °C.** This is now a three-axis
characterisation of that ratio, which is more than the corpus had for any branching ratio.

### 7c. Detection floors — **`0.000` is left-censored, not zero** [!]
S4 records `0.000` in **17 of 102 cells** and S5 in **many more**. Footnote d defines `-` (in
the RI column) as "not detected"; **the concentration column has no stated LOD or LOQ**
(Table S3 omits them, §2c). The smallest non-zero value printed anywhere in S4/S5 is
**0.012 μg L⁻¹** (2-thiophenecarboxaldehyde, S5 pH 8), so the practical detection floor is
**≲ 0.012 μg L⁻¹**. **Treat every `0.000` as `< 0.012 μg L⁻¹` (left-censored), not as a
measured zero.** This matters for the N-heterocycle block, which is `0.000` at all three
temperatures in S4 — the correct statement is **"below detection at all three temperatures",
which is still a strong result** (pyrazines require the exogenous amino nitrogen that the
TTCA-only system lacks; they appear as soon as pH 8 is used, S5, or amino acids are added,
main Table 1).

### 7d. Uncertainty to assign, given the SDs are unusable (§4c) [Z]
Recommended replacement uncertainty for Table S4 values, in the absence of trustworthy SDs:
- **Tier A compounds:** **± 15 % relative** (typical HS-SPME-GC-MS external-standard
  reproducibility with IS correction), floored at **± 0.02 μg L⁻¹**.
- **Tier B compounds:** **× ÷ 3** (response factor assumed = 1 vs 1,2-dichlorobenzene).
- **Class subtotals:** propagate the above; for the sulfur subtotal this gives ≈ **± 16 %**
  at 100 °C rising to ≈ **± 20 %** at 140 °C (Tier B's growing share).
- **The printed SDs, where quoted at all, average ≈ 5 % relative — implausibly tight for
  HS-SPME at sub-μg L⁻¹ levels, and demonstrably not reproducible between S4 and S5.**

---

## 8. CONSOLIDATED PARAMETER TABLE — everything extractable from this SI

**Common conditions unless stated:** TTCA **10 mmol L⁻¹** in water, sealed pressure-rated glass
vessel, oil bath, **120 min**, ice-bath quench, HS-SPME (CAR/PDMS/DVB **75 μm**) GC-MS,
IS = 1,2-dichlorobenzene, **n = 3**. Units **μg L⁻¹** unless stated
(**raster-verified**, §0b).

| # | parameter | value | units | condition | prov | source |
|---|---|---|---|---|---|---|
| **— TEMPERATURE LADDER (Table S4) — the gap-closing block —** |
| 1 | **MFT (2-methyl-3-furanthiol)** | **1.237** | μg L⁻¹ | 100 °C, pH 7, 120 min | **M** | S4 |
| 2 | **MFT** | **1.388** | μg L⁻¹ | 120 °C | **M** | S4 |
| 3 | **MFT** | **5.907** | μg L⁻¹ | 140 °C | **M** | S4 |
| 4 | **FFT (2-furfurylthiol)** | **3.734** | μg L⁻¹ | 100 °C | **M** | S4 |
| 5 | **FFT** | **4.107** | μg L⁻¹ | 120 °C | **M** | S4 |
| 6 | **FFT** | **11.439** | μg L⁻¹ | 140 °C | **M** | S4 |
| 7 | Furfural | **3.381 / 5.793 / 11.039** | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 8 | 2-Acetylthiazole | 3.079 / 8.795 / 9.858 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 9 | 3-Mercapto-2-butanone | 1.303 / 1.798 / 2.398 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 10 | 3-Mercapto-2-pentanone | 0.973 / 1.298 / 1.736 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 11 | Thiazole | 0.896 / 1.796 / 2.038 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 12 | **Sulfur-class subtotal** | **13.978 / 35.866 / 60.400** | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 13 | Sulfur-class subtotal, **Tier A only** | 13.715 / 33.856 / 54.196 | μg L⁻¹ | 100/120/140 °C | **Z** | §4d |
| 14 | O-heterocycle subtotal | 3.381 / 6.158 / 12.757 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 15 | N-heterocycle subtotal | **< 0.012 (BDL)** ×3 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| 16 | Grand total volatiles | 17.359 / 42.024 / 73.157 | μg L⁻¹ | 100/120/140 °C | **M** | S4 |
| **— DERIVED APPARENT ACTIVATION ENERGIES (yield at 120 min) —** |
| 17 | **Eₐ, MFT, 100→120 leg** | **7.0** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 18 | **Eₐ, MFT, 120→140 leg** | **97.8** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 19 | **Eₐ, FFT, 100→120 leg** | **5.8** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 20 | **Eₐ, FFT, 120→140 leg** | **69.2** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 21 | Eₐ, furfural (3-pt, well-behaved) | **37.9** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 22 | Eₐ, sulfur class (3-pt) | **46.9** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| 23 | Eₐ, grand total (3-pt) | **46.1** | kJ mol⁻¹ | pH 7 | **Z** | §7a |
| **— CYSTEINE KINETICS (Figs S3/S4, digitised) —** |
| 24 | **Eₐ, free-Cys depletion** | **55.1** | kJ mol⁻¹ | pH 7, aq., R² = 0.994 | **Z/D** | §6b |
| 25 | **A, free-Cys depletion** | **8.0 × 10⁴** | min⁻¹ | pH 7 | **Z/D** | §6b |
| 26 | k_eff, Cys depletion, 100 °C | **1.472 × 10⁻³** | min⁻¹ | 0–120 min | **Z/D** | §6b |
| 27 | k_eff, Cys depletion, 120 °C | **4.078 × 10⁻³** | min⁻¹ | 0–120 min | **Z/D** | §6b |
| 28 | k_eff, Cys depletion, 140 °C | **8.195 × 10⁻³** | min⁻¹ | 0–120 min | **Z/D** | §6b |
| 29 | Cys conversion at 120 min | **16.2 / 38.7 / 62.6** | % | 100/120/140 °C | **D** | §6b |
| 30 | **Max TTCA → free-Cys yield** | **≤ 16.3** | mol % | 140 °C, t = 40 min | **Z/D** | §6c |
| 31 | t_max of the free-Cys pool | **≈100 / ≈57 / ≈40** | min | 100/120/140 °C | **D** | §6c |
| **— pH RESPONSE (Table S5; T unresolved, §4c) —** |
| 32 | **MFT, pH 5.5 / 7 / 8** | 1.936 / 1.237 / 1.103 | μg L⁻¹ | 120 min | **M** | S5 |
| 33 | **FFT, pH 5.5 / 7 / 8** | 5.639 / 3.734 / 3.012 | μg L⁻¹ | 120 min | **M** | S5 |
| 34 | Thiol-class pH suppression | **× 0.56** (pH 5.5 → 8) | — | 120 min | **Z** | §5c |
| 35 | **Thiazole-class pH promotion** | **× 15.4** (pH 5.5 → 8) | — | 120 min | **Z** | §5c |
| 36 | Furan-class pH suppression | × 0.21 (pH 5.5 → 8) | — | 120 min | **Z** | §5c |
| 37 | Pyrazine formation threshold | **detected only at pH 8** | — | 120 min | **M** | S5 |
| 38 | **FFT/MFT ratio vs pH** | **2.91 / 3.02 / 2.73** | — | pH 5.5/7/8 | **Z** | §7b |
| 39 | **FFT/MFT ratio vs T** | **3.02 / 2.96 / 1.94** | — | 100/120/140 °C | **Z** | §7b |
| 40 | pH drift (initial → final) | 5.5→5.1 · 7→4.9 · 8→4.7 | — | 120 min | **M** | main p. 3242 |
| **— STORAGE STABILITY (Table S6) —** |
| 41 | TTCA retention, 4/25/40 °C | 99.12 / 95.53 / 92.94 | % | storage | **M** | S6 |
| 42 | ARP retention, 4/25/40 °C | 95.52 / 94.18 / 87.83 | % | storage | **M** | S6 |
| 43 | TTCA retention, pH 5.5/7/9 | 95.51 / 96.70 / 88.81 | % | storage | **M** | S6 |
| 44 | ARP retention, pH 5.5/7/9 | 94.23 / 90.52 / 78.75 | % | storage | **M** | S6 |
| 45 | **TTCA retention, a_w 0.113/0.432/0.843** | **86.49 / 72.47 / 64.23** | % | storage | **M** | S6 |
| 46 | **ARP retention, a_w 0.113/0.432/0.843** | **78.48 / 57.27 / 39.39** | % | storage | **M** | S6 |
| 47 | a_w of LiCl / K₂CO₃ / "dicalcium phosphate" | 0.113 / 0.432 / 0.843 | — | 25 °C | **M** | S2 |
| **— CALIBRATION —** |
| 48 | MFT calibration curve | y = 0.0039x + 0.0011, R² = 0.9989 | — | Tier A | **M** | S3 |
| 49 | FFT calibration curve | y = 0.0287x − 0.0003, R² = 0.9992 | — | Tier A | **M** | S3 |
| 50 | Furfural calibration curve | y = 0.0347x − 0.0012, R² = 0.9923 | — | Tier A | **M** | S3 |
| 51 | Practical detection floor | **≲ 0.012** | μg L⁻¹ | S4/S5 | **Z** | §7c |
| 52 | Tier A share of sulfur subtotal | 98.1 / 94.4 / 89.7 | % | 100/120/140 °C | **Z** | §4d |

---

## 9. **VERDICT: DOES KANG TABLE S4 CLOSE THE SULFUR TEMPERATURE-LADDER GAP?**

# ⭐ **YES — for MFT and FFT, with a Tier-A caveat that is satisfied, and three usage conditions that are not optional.**

**Per compound, with the Tier A/B adjudication:**

| compound | ladder present at 3 T? | **Tier** | evidence for the tier | **verdict** |
|---|---|---|---|---|
| **2-Methyl-3-furanthiol (MFT)** | **YES** — 1.237 / 1.388 / 5.907 | **A** | asterisked in main Table 1 **and** calibration curve in Table S3 (R² 0.9989) | **✅ GAP CLOSED.** Absolute, Tier A, three temperatures, n = 3, raster-verified μg L⁻¹. |
| **2-Furfurylthiol (FFT)** | **YES** — 3.734 / 4.107 / 11.439 | **A** | asterisked **and** in Table S3 (R² 0.9992) | **✅ GAP CLOSED.** Same status. |
| Furfural | YES — 3.381 / 5.793 / 11.039 | **A** | Table S3 (R² 0.9923); **externally confirmed by the main paper's own "3.3×" sentence** | **✅ bonus — and the best-behaved Arrhenius series in the table** |
| 3-Mercapto-2-butanone | YES | **A** | Table S3, R² 0.9999 | ✅ closed |
| 3-Mercapto-2-pentanone | YES | **A** | Table S3, R² 0.9997 | ✅ closed |
| 2-Acetylthiazole, thiazole, 2,4,5-trimethylthiazole, 2,5-dimethylthiazole, 2-ethylthiazole, 2-ethyl-4-methylthiazole, 5-ethyl-2,4-dimethylthiazole | YES | **A** | Table S3 | ✅ closed (thiazole branch fully laddered) |
| Thiophene, 2-methylthiophene, 3-thiophenethiol, 2-thiophenemethanethiol, the three thiophene-carboxaldehydes | YES | **A** | Table S3 | ✅ closed |
| 2-Methylthiazole | YES | **A⚠→ unusable** | its Table S3 curve is internally impossible (§2c) | ⚠️ **excluded** |
| Thieno[3,2-b]thiophene, 2-methylthieno[2,3-b]thiophene, 1,2,3-trithiolane, bis(2-furfuryl)sulfide, the dithio-furan, 3,3'-dithiobisthiophene | YES | **B** | no standard curve | **⚠️ PARTIAL** — ratio-comparable across T, absolute magnitude carries × ÷ 3 |
| **H₂S** | **NO** | — | never measured (§6a) | **❌ GAP OPEN.** Kang 2026 gives no H₂S number of any kind. |
| **Methional, dimethyl sulfide, dimethyl di/trisulfide, methanethiol** | **NO** | — | not in the analyte panel | **❌ GAP OPEN** for the volatile-sulfide sub-branch. |

**Overall: `YES` for the thiol/thiazole/thiophene sulfur branch (28 of 34 compounds Tier A);
`PARTIALLY` for the polysulfide/fused-thiophene sub-branch (6 compounds, Tier B);
`NO` for H₂S and for the low-MW volatile sulfides.**

**⚠️ THREE BINDING CONDITIONS ON USE — the ladder is real but it is not a rate ladder:**
1. **It is a yield-at-120-min ladder.** Fig. S4 shows the Cys precursor is only 16 % consumed
   at 100 °C but 63 % at 140 °C — **the three points are at very different extents of
   reaction.** An Eₐ read off them is an Eₐ *of accumulated yield*, not of a rate constant.
   (This is the same caution `kang2026_extraction.md` §5 raised for Fig. 1a; Table S4 does not
   remove it, and Fig. S4 now quantifies exactly how bad it is.)
2. **MFT and FFT are non-Arrhenius** (§7a). A single Eₐ fits MFT at R² = 0.781. **Do not report
   a single-Eₐ fit to these two series as if it were good.**
3. **The SDs are not trustworthy** (§4c **[!]**). **Use the means; assign uncertainty from §7d.**

**And two things the ladder is NOT:** it is a **closed aqueous model system** — purified TTCA
intermediate at 10 mmol L⁻¹, sealed superheated water, no lipid, no protein, no salt, no
carbohydrate matrix. **It bounds the intrinsic chemistry of the Cys/pentose sulfur branch; it
does not measure a food matrix**, and it must not be used to calibrate matrix partitioning.
Also, **TTCA's own identity is structurally unverified** (Table S1 is an empty caption, §2a).

---

## 10. PROPOSED FIT / HOLD-OUT ROLES — **DRAFT FOR ORCHESTRATOR**

*(Advisory only. Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration
was read for state or written. All role assignments are the orchestrator's to make.)*

### 10a. Recommended **FIT** rows
| # | proposed fit row | value | rationale | risk |
|---|---|---|---|---|
| **F-1** | **MFT + FFT at 100 °C and 120 °C** (4 values: 1.237, 1.388, 3.734, 4.107) | μg L⁻¹ | The two temperatures where the system is far from precursor exhaustion (16 % and 39 % Cys conversion) and where the yield-vs-rate confound is mildest. Tier A, raster-verified, externally anchored by main Table 1 at 120 °C. | Non-Arrhenius over the pair (Eₐ 6–7 kJ mol⁻¹) — the pair constrains **level**, barely constrains **slope**. |
| **F-2** | **Furfural ladder, all three T** (3.381 / 5.793 / 11.039) | μg L⁻¹ | The only well-behaved Arrhenius series in the table (Eₐ 32.8 → 43.5 across legs, ~15 % spread), **and independently confirmed by the main paper's own printed "3.3×" claim**. Best single calibration target here. | Furfural is an O-branch product; fitting it does not directly constrain the sulfur branch. |
| **F-3** | **Eₐ(free-Cys depletion) = 55.1 kJ mol⁻¹, R² 0.994** | kJ mol⁻¹ | The cleanest derived kinetic parameter in this SI, robust to the fitting window (±10 %), from a 7-point × 3-temperature digitisation with exact axis calibration. Constrains the **precursor consumption** step directly. | Digitised **[D]**, not printed; `system_identity: free_Cys` is 85 % confident (§6b). |
| **F-4** | **FFT/MFT branching ratio = 3.0 at ≤ 120 °C, pH-invariant** | — | Now characterised on three axes (T, pH, co-substrate). The pH-invariance (±5 % over 2.5 pH units) is a genuinely tight constraint. | Only valid with **no exogenous amino acid**; collapses to 0.15–9.2 when one is present (main §4h). |

### 10b. Recommended **HOLD-OUT** rows
| # | proposed hold-out | value | why hold out |
|---|---|---|---|
| **H-1** | ⭐ **MFT and FFT at 140 °C** (5.907 and 11.439 μg L⁻¹) | μg L⁻¹ | **The strongest hold-out this SI offers.** If the model is fitted on F-1 (100 + 120 °C), the 140 °C point is a **true extrapolation beyond the fitted range** — and it is exactly where the mechanism changes (§7a). A model with a single Eₐ fitted to F-1 predicts **MFT ≈ 1.55 and FFT ≈ 4.5** at 140 °C against measured **5.907 and 11.439** — a **3.8× and 2.5× under-prediction**. **This is a discriminating test, not a confirmatory one**, and the corpus's audit found only 0/8 true extrapolation rows. |
| **H-2** | **Sulfur-class subtotal ladder** (13.978 / 35.866 / 60.400) | μg L⁻¹ | Aggregate consistency check across the whole branch; ≥ 89.7 % Tier A (§4d) so it is quantitative. Independent of any single compound's calibration. |
| **H-3** | **Max TTCA → free-Cys yield ≤ 16.3 mol %** | mol % | A **mass-balance bound**, not a level — the kind of constraint a model can violate silently. Any model routing > 16 % of the cysteine moiety through a free-Cys pool is falsified. |
| **H-4** | **Thiazole-class pH promotion ×15.4 / thiol-class pH suppression ×0.56** | — | Directional gates over 2.5 pH units, both ≫ the corpus's ~2–3× threshold-penalty scale, and **internally validated by S5's arithmetic closure**. Safe as ratios even though S5's temperature is unresolved. |
| **H-5** | **Pyrazines below detection at 100/120/140 °C in the TTCA-only system, but present at pH 8** | BDL / present | A **presence/absence gate** — the sharpest kind. Falsifies any model that generates pyrazines without exogenous amino nitrogen. |

### 10c. **DO NOT USE**
| # | item | reason |
|---|---|---|
| **X-1** | **Fig. 1a's "Sulfur-containing compounds" bars** (15.2 / 14.3 / 19.4) | **Confirmed corrupt** — they are Table S5's pH series (§3). Corrupted 2.5× at 120 °C, 3.1× at 140 °C. **Superseded by Table S4 rows.** |
| **X-2** | **Any SD from Table S4 or Table S5** | Demonstrably not reproducible between the two tables for the identical experiment (§4c). Use §7d instead. |
| **X-3** | **Table S5's stated temperature (120 °C)** | Unresolved (§4c). Register `temperature: 100_or_120_C_UNRESOLVED`; **ratios only**. |
| **X-4** | **2-Methylthiazole, any value** | Its own calibration curve cannot produce its reported concentrations (§2c). |
| **X-5** | **A single-Eₐ Arrhenius fit to MFT or FFT across all 3 T** | R² = 0.781 for MFT; misprediction −38 % / +28 % (§7a). If a single Eₐ is unavoidable, record it as `non_Arrhenius: true, R2: 0.78`. |
| **X-6** | **Any H₂S parameter attributed to Kang 2026** | Never measured (§6a). The main paper's H₂S sentence is inference from the Cys curves. |
| **X-7** | **Table S6 cells composed across blocks** | The three storage blocks share no control point (§2d). Within-block ratios only. |
| **X-8** | **`0.000` treated as a measured zero** | Left-censored at ≲ 0.012 μg L⁻¹ (§7c). |
| **X-9** | **Any μ-bearing value read from `Kang2026.pdf`'s text layer** | Corrupted (§0b): "11.039 mg L⁻¹", "75 mm". **This SI is the authoritative unit source.** |

### 10d. Registry hygiene proposed for every row sourced here
```
source: Kang2026_SI            doi: 10.1039/D5FB00932D
system: TTCA_10mM_aqueous_sealed   NOT a food matrix
time: 120_min                  yield_not_rate: true
tier: A | B                    unit_basis: raster_verified_600dpi
sd_provenance: UNTRUSTWORTHY   uncertainty_source: dossier_7d
```

---

## 11. WHAT THIS SI DOES **NOT** CONTAIN — the remaining gaps

1. **Table S1's body** — the ¹³C NMR and MS/MS assignment of TTCA. **Empty caption (§2a).**
   ⇒ **TTCA's structure is unverified in the published record.** No purity figure exists in
   either document. **Everything in both papers rests on an unproven structural assignment.**
2. **Calibration ranges, LODs, LOQs, recoveries** for all 31 curves (§2c). Without a range,
   whether the 0.012–11.4 μg L⁻¹ readings are interpolated or extrapolated is unknowable.
3. **Any H₂S measurement** (§6a).
4. **Any time course of the volatiles.** Every volatile number in both documents is a single
   120-min endpoint. The only time courses in the whole work are the two Cys figures.
5. **n, and error bars, for Fig. S3.** Fig. S4 has error bars, Fig. S3 has none; neither states n.
6. **Storage duration for Table S6** — 18 retention percentages with no time axis.
7. **Fig. S1/S2 axis values** — no peak table, no labelled fragment m/z.
8. **Any matrix (lipid, protein, salt, carbohydrate) at all.**

**Highest-value follow-up retrievals, ranked:**
1. **The Kang group's earlier TTCA papers (main-paper refs 4 and 26)** — the source of the TTCA
   preparation, and the only plausible place a purity figure and the NMR assignment exist.
2. **Any Cys/pentose study reporting H₂S or methanethiol directly**, to close the gap X-6 names.
3. **A Cys/xylose volatile time course**, to convert the 120-min endpoints into rates and
   dissolve the yield-vs-rate confound flagged in §9 condition 1.

---

*End of `kang2026_SI_extraction.md`. Companion: `kang2026_extraction.md` (main paper).
No repo file outside `data/lit/extraction_dossiers/` was created or modified by this wave.*
