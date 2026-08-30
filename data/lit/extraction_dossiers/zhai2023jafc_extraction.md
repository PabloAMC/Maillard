# Zhai et al. 2023 (JAFC) — TTCA + additional cysteine: colour inhibition and flavour regulation

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★★ THIS PAPER IS THE PARENT OF KANG 2026's TABLE S4. THE REPO'S ONLY SULFUR TEMPERATURE LADDER IS NOT AN INDEPENDENT MEASUREMENT — IT IS THIS PAPER'S FIGURE 3, COLUMN E, RE-PUBLISHED. SEE §7.

---

# §0. ★ WRONG-FILE IDENTITY — REPORT THIS FIRST

**THE TWO ZHAI FILES ON DISK ARE SWAPPED RELATIVE TO THE K6a BRIEF.**

| brief said | file on disk | what the file actually is |
|---|---|---|
| `Zhai2023.pdf` = `10.1021/acs.jafc.3c04166` (TTCA + Cys, seven temperatures) | **`Zhai2023b.pdf`** (3 935 568 B) | ★ **THIS paper** — JAFC `10.1021/acs.jafc.3c04166` |
| `Zhai2023b.pdf` = `10.1016/j.foodchem.2022.134420` (TTCA + Xyl) | **`Zhai2023.pdf`** (1 567 301 B) | Food Chem. `10.1016/j.foodchem.2022.134420` — dossier `zhai2023foodchem_extraction.md` |

Any ingestion keyed on the filename stem has them the wrong way round. **This dossier is named
by paper identity, not by file stem** — the same convention wave K5a adopted for the two
Kocadağlı files.

### 0.1 ⚠️ AND THE BRIEF'S SCIENTIFIC DESCRIPTION IS ALSO WRONG IN ONE DECISIVE RESPECT

The brief expected *"SEVEN temperatures 90–150 °C … every MFT/FFT/volatile cell at every
temperature."* **The seven-temperature ladder exists but it is a BROWNING (A₄₂₀) ladder only.**

| measurement | temperatures actually run |
|---|---|
| Browning A₄₂₀, 120 min, pH 7, Cys 1:1 | ★ **90, 100, 110, 120, 130, 140, 150 °C** (Fig. 1b) — 7 rungs |
| **Volatiles (all 39 compounds, HS-SPME-GC-MS)** | **100, 120, 140 °C ONLY** (Fig. 3a–c) — 3 rungs |
| ARP / 3-DX / 1-DX / MGO / GO | **100, 120, 140 °C only** (Fig. 5) |
| A₄₂₀ time course | **100, 120, 140 °C only** (Fig. 2) |
| pH ladder (5.5–8.5, 7 values) | **120 °C only**, browning only (Fig. 1c) |

**⇒ There is NO seven-rung MFT or FFT ladder in this paper, and none exists in the corpus.
What this paper does contain, and what nobody in the repo has yet used, is a
3 temperature × 4 hold-time × 2 systems × 39 compounds numeric grid — 936 printed values —
which is worth far more than a seven-rung single-time ladder would have been (§6, §7, §8).**

---

## §1. IDENTITY `[M]`

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/Zhai2023b.pdf`, 3 935 568 B, 12 pages | `ls`, `pdfinfo` |
| title | **"Regulated Formation of Inhibited Color and Enhanced Flavor Derived from Heated 2-Threityl-Thiazolidine-4-Carboxylic Acid with Additional Cysteine Targeting at Different Degradation Stages"** | p. 14300 |
| authors | **Yun Zhai**, Heping Cui, Khizar Hayat, Teng Li, Xian Wu, **Yuying Fu\***, **Xiaoming Zhang\***, **Chi-Tang Ho\*** | p. 14300 |
| affiliations | Zhejiang Gongshang Univ. (Zhai, Li, Fu); Jiangnan Univ. State Key Lab of Food Science and Resources (Cui, Zhang); Miami Univ. Ohio (Hayat, Wu); Rutgers (Ho) | p. 14309 |
| **DOI** | ★ **`10.1021/acs.jafc.3c04166`** ✔ matches the brief's expected DOI | p. 14300 footer |
| journal | ***J. Agric. Food Chem.* 2023, 71 (39), 14300−14311** | running head, all 12 pages |
| dates | Received **20 Jun 2023** · Revised **31 Aug 2023** · Accepted **13 Sep 2023** · Published **25 Sep 2023** | p. 14300 |
| version | version of record, © 2023 American Chemical Society | p. 14300 |
| funding | NSFC **32172330**; National Key R&D Program of China **2017YFD0400105**; Jiangsu Postgraduate Research & Practice Innovation Program **KYCX20_1880**; National First-Class Discipline Program of Food Science and Technology **JUFSTR20180204** | p. 14309 |
| competing interests | *"The authors declare no competing financial interest."* | p. 14309 |
| SI | **NOT ON DISK.** Declared contents: UPLC-MS/MS + NMR data for TTCA/ARP (Fig. S1, Table S1); **calibration curves of volatile compounds (Table S2)**; main flavour compounds (μg/L) from the thermal-treated **MGO−GO−Cys** model; proposed generation pathway (Fig. S2) | p. 14309 |

### 1.1 ⚠️ An internal cross-reference defect in the paper
`Table S2` is introduced on p. 14302 as *"the corresponding calibration curves (Table S2)"* and
then re-used on p. 14308 as *"the Cys−MGO−GO thermal reaction model, the results shown in
**Table S2** also revealed its excellent ability of typical flavor formation."* The ASSOCIATED
CONTENT list distinguishes the two objects (calibration curves; MGO−GO−Cys flavour table), so
**one of the two in-text pointers is mislabelled.** Anyone ordering the SI should expect the
MGO−GO−Cys table under a different number.

---

## §2. ⚠️ THE μ → m RASTER CHECK — **PASSES, but the check was NOT optional**

This is an **Arbortext Publishing Engine / PDFlib** ACS PDF — the exact producer stack whose text
layer silently converted `μ` to `m` in `Kang2026.pdf` (Amendment 4 hazard). **Every unit below was
read off a 500–900 dpi raster (`pdftoppm -r 600/-r 900 -png`), not the text layer.**

| artefact | text layer says | **raster says** | verdict |
|---|---|---|---|
| IS loading, p. 14302 | `3 μL of internal standard (1,2-dichlorobenzene, 0.018 μg/μL in methanol)` | **`3 μL … 0.018 μg/μL`** — true Greek mu, glyph-verified | ✅ not corrupted |
| SPME fibre, p. 14302 | `CAR/PDMS/DVB, 75 μm` | **`75 μm`** | ✅ |
| column film, p. 14302 | `0.25 μm film thickness` | **`0.25 μm`** | ✅ |
| calibration formula, p. 14302 | `concentration of flavor compounds (μg/L)` | **`(μg/L)`** | ✅ |
| **all volatile values, Fig. 3a–c** | *(figure — no text layer at all)* | **raster-only; the whole grid in §6 is raster-read** | ✅ raster-only |
| totals quoted in text, p. 14305 | `29.565, 63.249, and 87.262 μg/L` | **`μg/L`** | ✅ |
| Fig. 5 y-axes | *(figure)* | **`Cncentration (mmol/L)`** — note the paper's own typo, "Cncentration", on all five panels | ✅ raster-only |
| Fig. 6 y-axes | *(figure)* | ★ **`Concentration (g/L)`** — **g/L, not μg/L: three orders of magnitude above the Fig. 3 scale, and it is correct, because Fig. 6 is a neat-substrate model reaction** | ✅ raster-only |

**⇒ This PDF's text layer is clean. But the Fig. 3 grid — the artefact this whole dossier is
built on — has NO text layer, so 100 % of §6 is a raster transcription and is subject to the
audit in §6.4 rather than to a text-layer diff.**

---

## §3. SYSTEM AND CONDITIONS `[M]`

### 3.1 Substrate preparation
- **TTCA synthesis:** Cys + Xyl, **0.0827 mol/L each** (equimolar), deionised water, adjusted to
  **pH 7.4 ± 0.01** with NaOH (2 and 6 mol/L), heated **90 °C / 40 min**, ice-bath quenched.
- **Purification:** Dowex 50WX4 ion-exchange resin, H⁺ form, 200−400 mesh; then semi-preparative
  RP-HPLC, XBridge amide column 4.6 × 150 mm, 3.5 μm.
- **Identity confirmation:** UPLC-ESI-MS (Waters Synapt MALDI Q-TOF) + NMR (Bruker DRX 400 MHz);
  spectra in Fig. S1 / Table S1 (**not on disk**).
- ⚠️ **No purity figure is printed anywhere in the main text.** The phrase used is *"which was used
  to produce TTCA purity"* — a sentence that does not parse and does not carry a number.
  **`[NEG]` TTCA purity: not reported.**

### 3.2 The reaction models
| parameter | value |
|---|---|
| TTCA | **10 mmol/L** |
| added Cys | **0.1, 0.5, 1, 2, 10 mol per mol TTCA** (dose series, browning only); **1:1 for every other experiment** |
| pH | **5.5, 6, 6.5, 7, 7.5, 8, 8.5** (browning only); **7.0 for every other experiment** |
| temperature | **90, 100, 110, 120, 130, 140, 150 °C** (browning only); **100, 120, 140 °C** for volatiles, α-dicarbonyls, and the A₄₂₀ time course |
| time | **20, 40, 60, 80, 100, 120, 140, 160, 180 min** (α-dicarbonyls, A₄₂₀); **40, 80, 120, 180 min** for volatiles |
| vessel | *"appropriate pressure-resistant bottles"* with a **collector-type magnetic stirrer (DF-101S)** — ★ **stirred**, unlike most sealed-vial ladders in the corpus |
| quench | ice bath |
| replication | **n = 3**, *"All measurements were conducted in triplicate"*; ANOVA + p < 0.05, SPSS 21 |
| ⚠️ buffer | ★ **NONE. The pH is set with NaOH and is unbuffered.** No buffer is named anywhere in the Materials section. **`[NEG]`** |
| ⚠️ measured pH drift | ★ **NOT REPORTED.** Unlike Kang 2026, which prints a measured post-reaction pH of 4.9, this paper reports **no pH at all after heating**. **`[NEG]`** |

### 3.3 The validation models (Fig. 6)
Neat two-component systems, all **pH 7.0 ± 0.1, 120 °C, pressure-resistant bottles, ice-quenched**:
- (NH₄)₂S **10 mmol/L** + furfural / furan / 2-methylfuran / 4-hydroxy-5-methyl-3(2H)-furanone,
  each at the **same molar ratio** (i.e. 10 mmol/L), sampled at **0, 10, 20 min** (Fig. 6a–d);
- (NH₄)₂S + **MGO + GO**, sampled at **0, 30, 60 min** (Fig. 6e).
- ⚠️ The Methods say the models were *"treated at 120 °C for a defined reaction time (10, 20, 30,
  and 40 min)"* but **the figure plots 0/10/20 for panels a–d and 0/30/60 for panel e.**
  **Methods and figure disagree on the time grid; the figure is used here.** `[!]`

---

## §4. ANALYTICAL METHOD AND QUANTIFICATION BASIS — **REQUIRED SECTION**

### 4.1 Volatiles
| parameter | value |
|---|---|
| method | HS-SPME−GC−MS |
| sample | **3 g** in a 20 mL headspace vial, PTFE/BYTL septum |
| **internal standard** | ★ **1,2-dichlorobenzene, 3 μL of a 0.018 μg/μL methanolic solution → 0.054 μg absolute per 3 g sample** |
| fibre | CAR/PDMS/DVB **75 μm**, exposed 2.5 cm into headspace |
| extraction | **60 °C, 20 min**, thermostatic water bath |
| desorption | 250 °C, 10 min |
| column | **DB-Wax 30 m × 0.25 mm × 0.25 μm** |
| instrument | Agilent 7890B GC + 5977B MSD |
| identification | NIST 17 + Kovats RI (C₇−C₃₀ n-alkanes) + WILEY 07 + literature |
| **quantification** | ★ **"The flavor compounds were quantified by the corresponding calibration curves (Table S2) constructed by the standards. x and y in the formula represented the concentration of flavor compounds (μg/L) and the ratio of peak area between the flavor standards and internal standard"** (p. 14302) |

### 4.2 ★ THE VERDICT ON THE QUANTIFICATION BASIS

> **This is ABSOLUTE quantification by external calibration with authentic standards, expressed
> as an IS-normalised peak-area ratio — NOT single-IS semi-quant — for every compound that has a
> calibration curve in Table S2.**

The formula is `y = f(x)` with `y` = *area(analyte)/area(IS)* and `x` = μg/L, built *"by the
standards"*. That is a genuine response-factor-carrying calibration; the IS is used only to
normalise injection/extraction variance. It is the same protocol Kang 2026 used, where the
Kang SI's Table S3 lists 29 calibration curves and the K5-wave dossier used them to split the
table into **Tier A** (has a curve) and **Tier B** (no curve, response factor assumed = 1).

⚠️ **BUT: Table S2 is not on disk, so the Tier A / Tier B split CANNOT be reproduced for this
paper.** Kang's Table S3 covered **29 of its 34 rows**; this paper reports **39** compounds, of
which at least the four extra ones (2,3-dihydro-5-methylthiophene, 4,5-dimethylthiazole,
1-(3-thienyl)-ethanone, 1-(2-thienyl)-1-propanone) have no counterpart in Kang's curve list.
**Working assumption, stated as an assumption: the compounds shared with Kang's Table S3 are
Tier A; the rest must be treated as Tier B (× ÷ 3) until Table S2 is obtained.**

**★ THE WAVE RULE, APPLIED EXPLICITLY.** Even for the Tier-B rows, **a constant unknown response
factor cancels in a ratio and therefore in an Arrhenius slope.** Every Ea in §8 is computed from
a ratio of two cells in the *same* row of the *same* table, so:
- **Ea derivations are legitimate for every row, Tier A and Tier B alike** — they are within-study
  shape, and the response factor divides out exactly.
- **What the Tier-B rows do NOT license:** absolute μg/L yields, cross-compound magnitude
  comparison, class subtotals treated as physical totals, and any comparison against another
  paper's absolute concentration.
- **What NOTHING here licenses:** comparing the μg/L headspace numbers to a liquid-phase
  concentration. These are **headspace** values from a 60 °C/20 min SPME equilibration of a
  3 g aliquot; the partition coefficient is not measured and is compound-specific.

### 4.3 Non-volatiles
| analyte | method | quantification |
|---|---|---|
| **A₄₂₀ browning** | Shimadzu 2100 UV−vis, 420 nm, samples *"diluted to an appropriate concentration"* | ⚠️ **dilution factor NOT reported → the A₄₂₀ values are NOT comparable between panels unless the dilution was constant, which is not stated.** `[!]` |
| **ARP** | HPLC-ELSD, XBridge BEH amide 3.5 μm 4.6 × 150 mm, 10 mM ammonium formate (pH 6) / MeCN gradient, 1 mL/min, 10 μL, 25 °C; ELSD drift tube **45 °C**, N₂ **1.5 L/min** | ★ **external standards** (standard ARP prepared in-lab) → absolute mmol/L |
| **α-dicarbonyls (3-DX, 1-DX, MGO, GO)** | OPD derivatisation (OPD 2 g/100 mL + DTPA 11 mmol/L in HEPES pH 7; 1:1 v/v; **dark, 25 °C, 4 h**), HPLC-DAD, SunFire C18 5 μm 4.6 × 150 mm | ★ **quinoxaline external standards** → absolute mmol/L |
| ⚠️ gradients | *"Linear gradient elution process is not presented here"* — **twice.** Neither HPLC method is reproducible from the paper. `[NEG]` |

---

## §5. FIGURE 1b — ★ THE SEVEN-TEMPERATURE BROWNING LADDER `[D]` (digitised)

**Conditions:** TTCA 10 mmol/L ± Cys (1 mol/mol), **pH 7, 120 min**, A₄₂₀.
**Digitisation method:** page 3 rendered at **600 dpi**; the y-axis tick rows were located
programmatically (ticks at 654.5 / 798.5 / 942.5 / 1086.5 / 1230.5 / 1374.5 px = 1.0 / 0.8 / 0.6 /
0.4 / 0.2 / 0.0, exactly 144 px per 0.2 A₄₂₀ unit), the x ticks likewise (2803.5 → 3729.5 px for
90 → 150 °C, 154.3 px per 10 °C), and marker centroids were extracted by dark-pixel clustering in
a ±12 px column band. **Quoted uncertainty ±0.02 A₄₂₀ units, which is the marker half-height.**

| T (°C) | **TTCA** A₄₂₀ | **TTCA−Cys** A₄₂₀ | inhibition (1 − ratio) |
|---:|---:|---:|---:|
| 90 | 0.24 ⚠ | 0.22 ⚠ | ~8 % |
| 100 | 0.25 ⚠ | 0.23 ⚠ | ~8 % |
| 110 | **0.346** | **0.267** | 23 % |
| 120 | **0.497** | **0.274** | 45 % |
| 130 | **0.572** | **0.361** | 37 % |
| 140 | **0.629** | **0.418** | 34 % |
| 150 | **0.867** | **0.463** | 47 % |

⚠️ **At 90 and 100 °C the two series' markers physically overlap** (merged dark blob 42–49 px tall
against a 26–30 px single marker), so those four values are the least certain; they are reported
as ±0.02 but should be read as "the two systems are indistinguishable at 90–100 °C", which is
also exactly what the authors claim in the text.

**`[M]` The only error bar printed in the whole panel** is on the 100 °C point (whiskers at
A₄₂₀ ≈ 0.294 and 0.163, i.e. **±0.065**). Every other point is drawn without one. **⇒ Treat
±0.065 (≈ ±25 % relative) as the panel's replicate scatter, not ±0.02, which is only the
digitisation error.**

### 5.1 `[D]` Apparent Ea of browning over seven rungs — **the cleanest Arrhenius object in the paper**
Least squares on ln(A₄₂₀) vs 1/T, all 7 rungs:

| system | **Ea (kJ mol⁻¹)** | R² | leg 90→120 | leg 120→150 |
|---|---:|---:|---:|---:|
| **TTCA** | **28.3** | **0.969** | 28.8 | 25.7 |
| **TTCA−Cys** | **16.9** | **0.948** | 8.7 | 24.2 |

**Caveat chain, in force:** (i) A₄₂₀ at a **fixed 120 min** is a yield, not a rate constant, and
the TTCA pot is visibly past its exponential phase by 120 min at 140–150 °C (Fig. 2c), which
**depresses** the top-leg slope; (ii) the dilution factor is unreported (§4.3), so if it varied
with condition the ladder is not on one scale; (iii) unbuffered — the pot's own acidification is
temperature-dependent and unmeasured; (iv) n = 3 with a ±25 % replicate scatter on the one point
that carries a bar. **Verdict `USE-Q`: usable as a browning-Ea prior with ±8 kJ mol⁻¹, and as a
`STRUCTURAL` statement that added Cys LOWERS the browning barrier's apparent value while
suppressing its magnitude — i.e. Cys inhibits browning more at low temperature than at high.**

**★ The Cys inhibition is NOT constant in temperature** — 8 % at 90–100 °C rising to 45 % at
120 °C. The authors' own summary ("temperature was the key parameter influencing the color
inhibition behavior of Cys") is supported by their own digitised curve.

---

## §6. ★★ FIGURE 3 — THE COMPLETE 936-VALUE GRID, TRANSCRIBED

**This is the artefact.** Figure 3a/b/c are **numeric heat maps**: every cell carries a printed
number in **μg/L**. 39 compounds × 8 columns × 3 temperatures = **936 values**, none of which has
ever been in the repo.

**Column key** (printed under every panel):
`A` = T-40 min · `B` = T−C-40 min · `C` = T-80 min · `D` = T−C-80 min ·
`E` = T-120 min · `F` = T−C-120 min · `G` = T-180 min · `H` = T−C-180 min
where **T** = TTCA alone (10 mmol/L, pH 7) and **T−C** = TTCA + Cys 1:1.

**Row key:** the identical 39-compound list is printed beside all three panels; I raster-verified
the list separately for panel a, panel b and panel c and **all three are the same list in the
same order.** Rows 1−29 = sulfur compounds, 30−33 = nitrogen heterocycles, 34−39 = oxygen
heterocycles (the grouping is the authors', per Fig. 4's class totals).

### 6a. Panel (a) — **100 °C**, μg/L `[M]`

| # | compound | A T-40 | B T−C-40 | C T-80 | D T−C-80 | **E T-120** | F T−C-120 | G T-180 | H T−C-180 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Thiophene | 0 | 0 | 0 | 0 | 0 | 0.004 | 0 | 0.015 |
| 2 | 2-Methylthiophene | 0 | 0 | 0 | 0 | 0 | 0.011 | 0 | 0.118 |
| 3 | 2,3-Dihydro-5-methylthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 4 | Thiazole | 0.097 | 0 | 0.289 | 0.149 | **0.896** | 0.279 | 0.423 | 0.578 |
| 5 | 2-Methylthiazole | 0.001 | 0 | 0.022 | 0 | **0.142** | 0.072 | 0.096 | 1.172 |
| 6 | 3-Mercapto-2-butanone | 0 | 0 | 0.198 | 0.142 | **1.303** | 0.742 | 1.298 | 1.768 |
| 7 | 2,5-Dimethyl-thiazole | 0 | 0 | 0.028 | 0 | **0.079** | 0.041 | 0.058 | 0.135 |
| **8** | **2-Methyl-3-furanthiol (MFT)** | 0.078 | 0.021 | 0.509 | 0.218 | **1.237** | 1.178 | 0.923 | 3.078 |
| 9 | 2-Ethyl-thiazole | 0 | 0 | 0.163 | 0 | **0.245** | 0.076 | 0.203 | 0.468 |
| 10 | 3-Mercapto-2-pentanone | 0.127 | 0 | 0.554 | 0.119 | **0.973** | 0.217 | 1.145 | 1.783 |
| 11 | 2,4,5-Trimethylthiazole | 0.004 | 0 | 0.392 | 0.028 | **0.419** | 0.179 | 0.457 | 0.936 |
| 12 | 4,5-Dimethyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.112 |
| **13** | **2-Furfurylthiol (FFT)** | 0.708 | 0.117 | 1.206 | 0.526 | **3.734** | 2.034 | 3.119 | 4.596 |
| 14 | 5-Ethyl-2,4-dimethyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0.078 | 0 | 0.362 |
| 15 | 2-Ethyl-4-methyl-thiazole | 0 | 0 | 0.012 | 0 | **0.117** | 0.098 | 0.213 | 0.308 |
| 16 | 3-Thiophenethiol | 0 | 0 | 0 | 0 | 0 | 0.038 | 0 | 0.179 |
| 17 | 2-Acetylthiazole | 0.634 | 0.079 | 0.814 | 0.126 | **3.079** | 2.589 | 1.519 | 4.309 |
| 18 | 2-Thiophenecarboxaldehyde | 0 | 0 | 0 | 0 | **0.090** | 0.129 | 0.140 | 0.232 |
| 19 | 5-Methyl-2-thiophenecarboxaldehyde | 0 | 0 | 0.235 | 0.072 | **0.936** | 0.738 | 0.427 | 0.936 |
| 20 | 2-Thiophenemethanethiol | 0 | 0 | 0 | 0 | **0.119** | 0.108 | 0.078 | 0.372 |
| 21 | 1-(3-Thienyl)-ethanone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.078 |
| 22 | 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0 | 0 | 0 | 0 | **0.021** | 0 | 0.059 | 0.063 |
| 23 | 1,2,3-Trithiolane | 0 | 0 | 0 | 0.018 | **0.208** | 0.117 | 0.154 | 0.213 |
| 24 | 3,3′-Dithiobisthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0.078 | 0.119 |
| 25 | Thieno[3,2-b]thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.179 |
| 26 | 1-(2-Thienyl)-1-propanone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 27 | 2,5-Thiophenedicarboxaldehyde | 0 | 0 | 0.095 | 0.024 | **0.346** | 0.179 | 0.079 | 0.472 |
| 28 | 2-Methylthieno[2,3-b]thiophene | 0 | 0 | 0 | 0 | 0 | 0.338 | 0 | 2.428 |
| 29 | Bis(2-furfuryl)sulfide | 0 | 0 | 0 | 0.019 | **0.034** | 0.026 | 0.248 | 0.217 |
| 30 | Methylpyrazine | 0 | 0 | 0 | 0 | 0 | 0 | 0.119 | 0.218 |
| 31 | Pyrazine | 0 | 0 | 0 | 0 | 0 | 0 | 0.024 | 0.019 |
| 32 | 3-Methyl-pyridine | 0 | 0 | 0 | 0 | 0 | 0 | 0.238 | 0.029 |
| 33 | 2,5-Dimethyl-pyrazine | 0 | 0 | 0 | 0 | 0 | 0.009 | 0.358 | 0.937 |
| 34 | Furan | 0 | 0 | 0 | 0.039 | 0 | 0.048 | 0 | 0 |
| 35 | 2-Methyl-furan | 0 | 0 | 0 | 0.140 | 0 | 0.198 | 0 | 0 |
| **36** | **Furfural** | 0.307 | 0 | 2.891 | 1.129 | **3.381** | 2.897 | 2.293 | 2.136 |
| 37 | 1-(2-Furanyl)-ethanone | 0 | 0 | 0 | 0.029 | 0 | 0 | 0 | 0 |
| 38 | 2(5H)-Furanone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 39 | 4-Hydroxy-5-methyl-3(2H)-furanone | 0.103 | 0.049 | 0 | 0 | 0 | 0 | 0 | 0 |
| — | **SULFUR subtotal `[D]`** | 1.649 | 0.217 | 4.517 | 1.441 | **13.978** | 9.271 | 10.717 | 26.226 |
| — | **NITROGEN subtotal `[D]`** | 0 | 0 | 0 | 0 | **0** | 0.009 | **0.739** | **1.203** |
| — | **OXYGEN subtotal `[D]`** | 0.410 | 0.049 | 2.891 | 1.337 | **3.381** | 3.143 | 2.293 | 2.136 |
| — | **GRAND TOTAL `[D]`** | 2.059 | 0.266 | 7.408 | 2.778 | **17.359** | 12.423 | **13.749** | **29.565** |

### 6b. Panel (b) — **120 °C**, μg/L `[M]`

| # | compound | A T-40 | B T−C-40 | C T-80 | D T−C-80 | **E T-120** | F T−C-120 | G T-180 | H T−C-180 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Thiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0.028 | 0.013 |
| 2 | 2-Methylthiophene | 0.736 | 0.075 | 1.109 | 0.243 | **1.128** | 0.237 | 1.134 | 0.317 |
| 3 | 2,3-Dihydro-5-methylthiophene | 0 | 0 | 0 | 0 | 0 | 0.058 | 0 | 0.117 |
| 4 | Thiazole | 0.632 | 0.037 | 1.147 | 0.343 | **1.796** | 0.461 | 1.213 | 0.928 |
| 5 | 2-Methylthiazole | 0.024 | 0 | 0.079 | 0.487 | **0.194** | 0.579 | 0.172 | 1.106 |
| 6 | 3-Mercapto-2-butanone | 0.034 | 0 | 1.217 | 0 | **1.798** | 0 | 1.824 | 0.228 |
| 7 | 2,5-Dimethyl-thiazole | 0.272 | 0 | 1.148 | 0.102 | **2.203** | 0.127 | 1.918 | 0.179 |
| **8** | **MFT** | 0.432 | 1.104 | 0.937 | 1.741 | **1.388** | 4.083 | 2.479 | 5.032 |
| 9 | 2-Ethyl-thiazole | 0.079 | 0 | 0.178 | 0.142 | **0.273** | 0.152 | 0.324 | 0.197 |
| 10 | 3-Mercapto-2-pentanone | 0.927 | 0 | 1.173 | 0 | **1.298** | 0.438 | 1.136 | 0.527 |
| 11 | 2,4,5-Trimethylthiazole | 1.178 | 0 | 2.928 | 0.821 | **3.199** | 0.928 | 3.012 | 1.376 |
| 12 | 4,5-Dimethyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0.012 | 0 | 0.027 |
| **13** | **FFT** | 1.129 | 2.064 | 3.019 | 5.073 | **4.107** | 11.591 | 4.279 | 12.937 |
| 14 | 5-Ethyl-2,4-dimethyl-thiazole | 0.127 | 0 | 0.217 | 0.061 | **0.587** | 0.174 | 0.419 | 0.398 |
| 15 | 2-Ethyl-4-methyl-thiazole | 0.239 | 0 | 0.739 | 0.225 | **1.278** | 0.284 | 1.304 | 0.727 |
| 16 | 3-Thiophenethiol | 0 | 0 | 0.028 | 0 | **0.178** | **12.892** | 0.192 | **13.728** |
| 17 | 2-Acetylthiazole | 2.392 | 0 | 4.298 | 0 | **8.795** | 0.538 | 7.298 | 1.178 |
| 18 | 2-Thiophenecarboxaldehyde | 0.732 | 0 | 1.178 | 0 | **2.196** | 1.060 | 1.689 | 1.893 |
| 19 | 5-Methyl-2-thiophenecarboxaldehyde | 0.172 | 0 | 2.217 | 1.584 | **3.649** | 2.535 | 4.027 | 3.783 |
| 20 | 2-Thiophenemethanethiol | 0 | 0 | 0.078 | 0 | **0.196** | 0.483 | 0.278 | 0.794 |
| 21 | 1-(3-Thienyl)-ethanone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.089 |
| 22 | 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0 | 0 | 0 | 0 | **0.098** | 0 | 0.112 | 0.172 |
| 23 | 1,2,3-Trithiolane | 0.021 | 0 | 0.178 | 0 | **0.378** | 0 | 0.217 | 0.042 |
| 24 | 3,3′-Dithiobisthiophene | 0 | 0 | 0 | 0 | 0 | 0 | 0.117 | 1.783 |
| 25 | Thieno[3,2-b]thiophene | 0 | 0 | 0.273 | 1.727 | **1.302** | 7.967 | 1.402 | 8.737 |
| 26 | 1-(2-Thienyl)-1-propanone | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.028 |
| 27 | 2,5-Thiophenedicarboxaldehyde | 0.127 | 0 | 0.289 | 0 | **0.724** | 1.261 | 0.637 | 1.378 |
| 28 | 2-Methylthieno[2,3-b]thiophene | 0 | 0 | 0.089 | 0.334 | **0.190** | 1.886 | 0.171 | 2.279 |
| 29 | Bis(2-furfuryl)sulfide | 0 | 0 | 0.017 | 0.146 | **0.042** | 0.972 | 0.058 | 1.172 |
| 30 | Methylpyrazine | 0 | 0 | 0 | 0.214 | 0 | 0.343 | 0 | 0.793 |
| 31 | Pyrazine | 0 | 0 | 0 | 0 | 0 | 0.063 | 0 | 0.119 |
| 32 | 3-Methyl-pyridine | 0 | 0 | 0 | 0 | 0 | 0.076 | 0 | 0.124 |
| 33 | 2,5-Dimethyl-pyrazine | 0 | 0 | 0 | 0 | 0 | 0.679 | 0.028 | 1.048 |
| 34 | Furan | 0 | 0 | 0 | 0 | **0.138** | 0 | 0.214 | 0 |
| 35 | 2-Methyl-furan | 0.213 | 0 | 0.192 | 1.280 | **0.019** | 0 | 0.017 | 0 |
| **36** | **Furfural** | 3.179 | 0 | 4.378 | 0.348 | **5.793** | 0 | 3.289 | 0 |
| 37 | 1-(2-Furanyl)-ethanone | 0.021 | 0 | 0.078 | 0 | **0.136** | 0 | 0.099 | 0 |
| 38 | 2(5H)-Furanone | 0.012 | 0 | 0.042 | 0 | **0.072** | 0 | 0.081 | 0 |
| 39 | 4-Hydroxy-5-methyl-3(2H)-furanone | 0.067 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| — | **SULFUR subtotal `[D]`** | 9.253 | 3.280 | 22.536 | 13.029 | **36.997** | 48.718 | 35.440 | 61.165 |
| — | **NITROGEN subtotal `[D]`** | 0 | 0 | 0 | 0.214 | 0 | 1.161 | 0.028 | 2.084 |
| — | **OXYGEN subtotal `[D]`** | 3.492 | 0 | 4.690 | 1.628 | **6.158** | 0 | 3.700 | 0 |
| — | **GRAND TOTAL `[D]`** | 12.745 | 3.280 | 27.226 | 14.871 | **43.155** | 49.879 | **39.168** | **63.249** |

### 6c. Panel (c) — **140 °C**, μg/L `[M]`

| # | compound | A T-40 | B T−C-40 | C T-80 | D T−C-80 | **E T-120** | F T−C-120 | G T-180 | H T−C-180 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | Thiophene | 0 | 0 | 0.238 | 0 | **0.118** | 0.072 | 0.092 | 2.054 |
| 2 | 2-Methylthiophene | 1.112 | 0.172 | 2.129 | 0.472 | **0.892** | 0.791 | 0.542 | 2.079 |
| 3 | 2,3-Dihydro-5-methylthiophene | 0 | 0 | 0 | 0.012 | 0 | 0.127 | 0 | 0.308 |
| 4 | Thiazole | 0.832 | 0.116 | 1.992 | 0.478 | **2.038** | 2.429 | 2.012 | 3.038 |
| 5 | 2-Methylthiazole | 0.078 | 0.019 | 0.372 | 0.112 | **0.254** | 1.108 | 0.241 | 2.482 |
| 6 | 3-Mercapto-2-butanone | 1.124 | 0.024 | 2.799 | 0.708 | **2.398** | 1.713 | 2.038 | 2.014 |
| 7 | 2,5-Dimethyl-thiazole | 0.912 | 0.213 | 2.439 | 0.458 | **2.553** | 2.136 | 2.515 | 3.148 |
| **8** | **MFT** | 1.129 | 1.576 | 3.043 | 3.126 | **5.907** | 7.189 | 4.012 | 5.037 |
| 9 | 2-Ethyl-thiazole | 0.304 | 0 | 0.598 | 0.078 | **0.289** | 0.279 | 0.178 | 0.408 |
| 10 | 3-Mercapto-2-pentanone | 1.274 | 0.502 | 2.986 | 1.102 | **1.736** | 1.378 | 1.539 | 1.889 |
| 11 | 2,4,5-Trimethylthiazole | 1.179 | 0.081 | 3.783 | 1.386 | **3.843** | 4.029 | 3.628 | 4.478 |
| 12 | 4,5-Dimethyl-thiazole | 0 | 0 | 0 | 0 | 0 | 0.028 | 0 | 0.432 |
| **13** | **FFT** | 5.079 | 2.728 | **12.439** | 6.963 | **11.439** | 17.763 | **4.906** | 16.628 |
| 14 | 5-Ethyl-2,4-dimethyl-thiazole | 0.578 | 0.027 | 0.772 | 0.118 | **0.889** | 1.137 | 0.804 | 1.189 |
| 15 | 2-Ethyl-4-methyl-thiazole | 1.179 | 0.162 | 2.326 | 0.773 | **2.112** | 1.127 | 1.837 | 2.436 |
| 16 | 3-Thiophenethiol | 0.035 | 0 | 0.178 | 0.027 | **0.268** | 0.179 | 0.241 | 1.038 |
| 17 | 2-Acetylthiazole | 5.869 | 0.098 | 7.376 | 1.936 | **9.858** | 4.473 | 10.035 | 6.069 |
| 18 | 2-Thiophenecarboxaldehyde | 2.389 | 0.112 | 3.079 | 0.578 | **2.280** | 3.378 | 2.185 | 4.809 |
| 19 | 5-Methyl-2-thiophenecarboxaldehyde | 1.178 | 0.716 | 3.239 | 1.345 | **5.706** | 5.628 | 5.281 | 6.439 |
| 20 | 2-Thiophenemethanethiol | 0 | 0 | 0.069 | 0.018 | **0.332** | 0.978 | 0.428 | 1.172 |
| 21 | 1-(3-Thienyl)-ethanone | 0 | 0 | 0 | 0 | 0 | 0.021 | 0 | 0.179 |
| 22 | 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0 | 0 | 0.049 | 0 | **0.137** | 0.048 | 0.213 | 0.378 |
| 23 | 1,2,3-Trithiolane | 0.234 | 0 | 0.519 | 0.103 | **0.497** | 0.528 | 0.517 | 0.727 |
| 24 | 3,3′-Dithiobisthiophene | 0 | 0 | 0 | 0 | 0 | 0.028 | 0.149 | 0.112 |
| 25 | Thieno[3,2-b]thiophene | 0.089 | 0 | 1.479 | 1.489 | **3.369** | 2.035 | 3.439 | 4.138 |
| 26 | 1-(2-Thienyl)-1-propanone | 0 | 0 | 0 | 0 | 0 | 0.102 | 0 | 0.209 |
| 27 | 2,5-Thiophenedicarboxaldehyde | 0.179 | 0.052 | 1.178 | 0.771 | **1.283** | 1.045 | 1.246 | 2.256 |
| 28 | 2-Methylthieno[2,3-b]thiophene | 0.489 | 0 | 1.197 | 0.428 | **2.123** | 1.039 | 2.039 | 2.802 |
| 29 | Bis(2-furfuryl)sulfide | 0 | 0 | 0.048 | 0.372 | **0.078** | 2.198 | 1.126 | 3.028 |
| 30 | Methylpyrazine | 0 | 0.102 | 0 | 0.121 | 0 | 2.892 | 0.372 | 3.047 |
| 31 | Pyrazine | 0 | 0 | 0 | 0 | 0 | 0.386 | 0 | 0.582 |
| 32 | 3-Methyl-pyridine | 0 | 0 | 0 | 0 | 0 | 0.218 | 0 | 0.309 |
| 33 | 2,5-Dimethyl-pyrazine | 0 | 0 | 0 | 0.189 | 0 | 1.178 | 0.568 | 2.348 |
| 34 | Furan | 0.432 | 0 | 1.328 | 0 | **1.108** | 0 | 0.119 | 0 |
| 35 | 2-Methyl-furan | 0.719 | 0 | 1.127 | 0 | **0.133** | 0 | 0 | 0 |
| **36** | **Furfural** | 4.497 | 0.073 | 10.089 | 0 | **11.039** | 0 | 12.386 | 0 |
| 37 | 1-(2-Furanyl)-ethanone | 0.126 | 0 | 0.425 | 0 | **0.379** | 0 | 0.221 | 0 |
| 38 | 2(5H)-Furanone | 0 | 0 | 0 | 0 | **0.098** | 0 | 0.126 | 0 |
| 39 | 4-Hydroxy-5-methyl-3(2H)-furanone | 0 | 0.124 | 0 | 0 | 0 | 0 | 0 | 0 |
| — | **SULFUR subtotal `[D]`** | 25.242 | 6.598 | 54.327 | 22.853 | **60.399** | 62.986 | 51.243 | 80.976 |
| — | **NITROGEN subtotal `[D]`** | 0 | 0.102 | 0 | 0.310 | 0 | 4.674 | 0.940 | 6.286 |
| — | **OXYGEN subtotal `[D]`** | 5.774 | 0.197 | 12.969 | 0 | **12.757** | 0 | 12.852 | 0 |
| — | **GRAND TOTAL `[D]`** | 31.016 | 6.897 | 67.296 | 23.163 | **73.156** | 67.660 | **65.035** | **87.262** |

### 6.4 ★ INTERNAL-CONSISTENCY AUDIT OF MY OWN TRANSCRIPTION — **five of six checks pass exactly**

The paper quotes six column totals in its running text (p. 14305). I recomputed all six from my
raster transcription **before** reading the sentence back:

| quantity | **my recomputed sum** | **paper's printed text** | Δ |
|---|---:|---:|---|
| Total volatiles, TTCA, 100 °C, 180 min | **13.749** | 13.749 | **0.000** ✔ |
| Total volatiles, TTCA, 120 °C, 180 min | **39.168** | 39.168 | **0.000** ✔ |
| Total volatiles, TTCA, 140 °C, 180 min | **65.035** | **63.035** | ⚠️ **2.000** |
| Total volatiles, TTCA−Cys, 100 °C, 180 min | **29.565** | 29.565 | **0.000** ✔ |
| Total volatiles, TTCA−Cys, 120 °C, 180 min | **63.249** | 63.249 | **0.000** ✔ |
| Total volatiles, TTCA−Cys, 140 °C, 180 min | **87.262** | 87.262 | **0.000** ✔ |
| **Pyrazine total (4 rows), TTCA, 100 °C, 180 min** | **0.739** | 0.739 | **0.000** ✔ |
| **Pyrazine total (4 rows), TTCA−Cys, 100 °C, 180 min** | **1.203** | 1.203 | **0.000** ✔ |

**Eight independent checks, seven exact to the last decimal.** The one failure is a **2.000
exactly** discrepancy at 140 °C/T/180 min. Column G of panel (c) was re-rastered and re-read at
900 dpi cell by cell; the transcription is confirmed. **Fig. 4c's stacked bar for that cell also
sits at ≈ 65, above the 140 °C/T-120 bar's neighbour and consistent with 65.035.**
**⇒ `[!]` The paper's text figure 63.035 μg/L is a typographical error for 65.035 μg/L. The
figure is right; the sentence is wrong. Do not ingest 63.035.**

⚠️ **`[NEG]` NO STANDARD DEVIATIONS ANYWHERE IN FIGURE 3.** The Methods promise mean ± SD, n = 3,
but the heat map prints bare means. **The 936-value grid carries no uncertainty at all.**
Recommended replacement uncertainty, by the same logic wave K5 applied to Kang Table S4:
**±15 % relative for compounds shared with Kang's calibration list (Tier A), × ÷ 3 for the rest
(Tier B), floored at ±0.02 μg/L.** `0` should be read as **left-censored below ≈ 0.001 μg/L**
(the smallest non-zero value printed anywhere in the grid is 0.001, panel a row 5 col A).

---

# §7. ★★★ THE HEADLINE FINDING — KANG 2026's TABLE S4 **IS** THIS FIGURE'S COLUMN E

The repo's only sulfur temperature ladder — the one wave K4c/K5 recovered from
`Kang2026_supplementary.pdf` Table S4, the one Amendment 8/9 built the `kang_140C_*` and
`kang_switch_on_*` hold-out rows on — **is a re-publication of the TTCA-alone, 120-minute column
of this 2023 figure.**

### 7.1 The evidence, cell by cell

Kang Table S4 has **34 compound rows**; this figure has **39**. Kang's 34 are a strict subset:
the five it omits are rows **3** (2,3-dihydro-5-methylthiophene), **12** (4,5-dimethylthiazole),
**21** (1-(3-thienyl)-ethanone), **26** (1-(2-thienyl)-1-propanone) and **39**
(4-hydroxy-5-methyl-3(2H)-furanone) — **all five of which are 0.000 in column E at all three
temperatures**, so the omission changes no subtotal.

Over the **102 shared cells** (34 rows × 3 temperatures):

| | agreement |
|---|---|
| cells identical to the last printed decimal | ★ **101 of 102** |
| the single exception | **2-Methylthiophene, 120 °C: Kang prints `0.000`, Zhai prints `1.128`** |
| Sulfur subtotal, 100 °C | Zhai col E **13.978** = Kang **13.978** ✔ |
| Sulfur subtotal, 140 °C | Zhai col E **60.399** = Kang's recomputed **60.399** (printed 60.400) ✔ |
| Sulfur subtotal, 120 °C | Zhai **36.997** vs Kang **35.869** — **the difference is 1.128, i.e. exactly the one discrepant cell** ✔ |
| O-heterocycle subtotal, all three | **3.381 / 6.158 / 12.757** = Kang, exactly ✔ |
| N-heterocycle subtotal, all three | **0 / 0 / 0** = Kang's BDL ✔ |
| Grand total, 100 °C | **17.359** = Kang **17.359** ✔ |
| Grand total, 140 °C | **73.156** = Kang **73.157** (0.001 rounding) ✔ |
| MFT | **1.237 / 1.388 / 5.907** = Kang, exactly ✔ |
| FFT | **3.734 / 4.107 / 11.439** = Kang, exactly ✔ |
| Furfural | **3.381 / 5.793 / 11.039** = Kang, exactly ✔ |
| FFT/MFT ratio | **3.02 / 2.96 / 1.94** = the K5-wave's derived ratio for Kang, exactly ✔ |

**The probability that two independent HS-SPME runs, in different years, agree to three decimal
places in 101 of 102 cells is nil.** These are the same measurements.

### 7.2 And the authorship confirms it
| | Zhai 2023 JAFC | Kang 2026 *Sust. Food Technol.* |
|---|---|---|
| corresponding authors | **Yuying Fu**, Xiaoming Zhang, Chi-Tang Ho | **Yuying Fu**, **Yun Zhai** |
| first author of the other paper | **Yun Zhai** | (Zhai is corresponding author) |
| institution of the ladder work | Zhejiang Gongshang University | Zhejiang Gongshang University |
| system | TTCA 10 mmol/L, pH 7, sealed pressure bottle, HS-SPME/DB-Wax, IS = 1,2-dichlorobenzene 0.018 μg/μL | **identical in every particular** |

**Kang 2026 is a 2026 re-publication of a 2023 data set by the same group. It is not a
replication.**

### 7.3 ★★ WHAT THIS COSTS THE REPO — three consequences, all binding

1. **THE SULFUR LADDER HAS n = 1, NOT n = 2.** `research_round4_nulls.md` §C.2 recorded
   "NULL SURVIVES for MFT and FFT — no second sulfur temperature ladder." **That null is now
   worse than it looked: the corpus's *first* ladder is also its only one, and one of the two
   papers that appeared to carry it is a duplicate.** Any uncertainty budget that treated
   Kang-140 °C as corroborated by a second source must be re-opened.
2. **A FIT/HOLD-OUT COLLISION IS NOW POSSIBLE.** If any wave ingests this paper's Fig. 3 column E
   as a new anchor while `kang_140C_MFT` / `kang_140C_FFT` / `kang_switch_on_*` remain
   hold-outs, **the model would be fitted on its own hold-out.** ★ **The 120-min column of this
   figure must be marked as ALREADY-SEEN and must never be declared as an independent row.**
   Columns A/B/C/D/F/G/H (the other three hold times and the whole +Cys arm) are genuinely new
   and unpublished elsewhere — **they are the fresh material, and they are 819 of the 936 values.**
3. **The one discrepant cell is a data-integrity flag, not a rounding difference.**
   2-Methylthiophene at 120 °C is `1.128` in the 2023 figure and `0.000` in the 2026 SI table.
   One of the two has been altered. It propagates into Kang's printed sulfur subtotal
   (35.866 vs the 36.997 the parent figure implies) and therefore into every ratio computed from
   that subtotal. **The K5 dossier's audit of Kang S4 could not have caught this, because Kang's
   table is internally consistent with its own altered cell.**

### 7.4 ⚠️ AND IT COMPOUNDS A DUPLICATION THE CORPUS ALREADY KNEW ABOUT
`kang2026_SI_extraction.md` §4c already found that **Kang's Table S5 pH-7 column IS Kang's
Table S4 100 °C column.** With this dossier, the same 34-number vector is now known to appear in
**three** places: Zhai 2023 Fig. 3a col E → Kang Table S4 100 °C → Kang Table S5 pH 7.
**Any weighting scheme that counts those as three observations is wrong by a factor of three.**

---

# §8. ★★ THE Ea ANALYSIS — AND THE KANG "SWITCH-ON" DOES **NOT** SURVIVE THE TIME AXIS

## 8.1 What the switch-on claim was
`kang2026_SI_extraction.md` §7a, the finding this whole build lane was designed around:

> "MFT and FFT are strongly non-Arrhenius: a threshold sits between 120 and 140 °C … MFT 1.12×
> then 4.26×; FFT 1.10× then 2.78× … apparent Eₐ climbing from ~6–7 to 70–98 kJ mol⁻¹ … **A
> saturation artefact cannot produce this** — precursor depletion *depresses* apparent Eₐ at the
> top of a ladder, which is exactly what the class does and the opposite of what the thiols do."

That claim rests on **one column of one table at one hold time (120 min)**. This paper contains
that column **and the three others**.

## 8.2 ★ The same experiment at four hold times `[D]`

**MFT (row 8), TTCA alone, μg/L, and the two-leg fold changes and apparent Eₐ:**

| hold time | 100 °C | 120 °C | 140 °C | ×(100→120) | ×(120→140) | Eₐ low leg | **Eₐ high leg** | Eₐ overall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 min | 0.078 | 0.432 | 1.129 | **5.54×** | 2.61× | 104.4 | 64.9 | 85.6 |
| 80 min | 0.509 | 0.937 | 3.043 | 1.84× | 3.25× | 37.2 | 79.5 | 57.3 |
| **120 min** ← *Kang's column* | **1.237** | **1.388** | **5.907** | **1.12×** | **4.26×** | **7.0** | **97.8** | **50.1** |
| 180 min | 0.923 | 2.479 | 4.012 | 2.69× | 1.62× | 60.3 | **32.5** | 47.1 |

**FFT (row 13), TTCA alone:**

| hold time | 100 °C | 120 °C | 140 °C | ×(100→120) | ×(120→140) | Eₐ low leg | **Eₐ high leg** | Eₐ overall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 min | 0.708 | 1.129 | 5.079 | 1.59× | **4.50×** | 28.5 | 101.5 | 63.1 |
| 80 min | 1.206 | 3.019 | 12.439 | 2.50× | 4.12× | 56.0 | 95.6 | 74.8 |
| **120 min** ← *Kang's column* | **3.734** | **4.107** | **11.439** | **1.10×** | **2.79×** | **5.8** | **69.2** | **35.9** |
| 180 min | 3.119 | 4.279 | **4.906** | 1.37× | **1.15×** | 19.3 | **9.2** | **14.5** |

## 8.3 ★★ THE VERDICT — **the switch-on is a hold-time artefact, and this is the wave's most consequential result**

- **It is not reproducible across hold time in its own experiment.** The "flat low leg then jump"
  signature appears at **120 min and nowhere else**. At 40 min MFT does the *opposite*
  (5.54× then 2.61× — decelerating). At 180 min **both** thiols decelerate, and FFT's high-leg
  apparent Eₐ collapses from 69.2 to **9.2 kJ mol⁻¹**.
- **The mechanism is visible in the raw time courses**, because at 140 °C the thiols have already
  **peaked and turned over** inside the measurement window:

| series (TTCA alone) | 40 min | 80 min | 120 min | 180 min | peak |
|---|---:|---:|---:|---:|---|
| **FFT @ 140 °C** | 5.079 | **12.439** | 11.439 | **4.906** | ★ **peaks at 80 min, then falls 2.5×** |
| FFT @ 120 °C | 1.129 | 3.019 | 4.107 | 4.279 | still rising at 180 min |
| FFT @ 100 °C | 0.708 | 1.206 | **3.734** | 3.119 | peaks ~120 min |
| **MFT @ 140 °C** | 1.129 | 3.043 | **5.907** | 4.012 | peaks ~120 min |
| MFT @ 120 °C | 0.432 | 0.937 | 1.388 | **2.479** | still rising at 180 min |
| MFT @ 100 °C | 0.078 | 0.509 | **1.237** | 0.923 | peaks ~120 min |

  **A fixed-time slice across a family of curves whose maxima move left as temperature rises does
  not measure an activation energy — it measures where each curve happens to be relative to its
  own peak.** At 120 min the 140 °C FFT curve is just past its peak while the 120 °C curve is
  still climbing, which *inflates* the 120→140 leg; at 180 min the 140 °C curve has fallen well
  past its peak while the 120 °C one is still climbing, which *deflates* it to 1.15×.
- **The K5 dossier's own falsification argument is what fails.** It argued the switch-on could not
  be a saturation artefact because saturation *depresses* the top-leg Eₐ. That reasoning is
  correct for a monotone-rising system. **It does not hold when the top-temperature curve has a
  maximum inside the window, because then the fixed-time slice can sit either side of the peak
  and the sign of the bias flips with hold time — which is precisely what this table shows.**

## 8.4 ★ WHAT THIS SAYS ABOUT THE B2.2 / B2.3 KANG BLOCK AND THE YILTIRAK FAILURE

`kinetic_core_b2_2_diagnosis.md` §3 already refused to count the `kang_140C_MFT` pass, on the
independent ground that the model's pH trajectory was six units wrong at the two fitted rungs.
**This dossier adds a second, prior reason: the target itself is not a barrier.**

- ★ **`kang_switch_on_*` should be retired as a scored row, not merely failed.** It asks the model
  to reproduce a feature that the source experiment does not reproduce at any other hold time.
  Fitting a two-regime or threshold Ea form to it — which
  `kang2026_SI_extraction.md` §7a explicitly recommended — **would be fitting a sampling
  artefact into the physics.**
- ★ **This is the same failure mode as the Yiltirak family, seen from the other side.**
  `cutover_final_exam.md` diagnoses Yiltirak as a TIME-axis failure: *"Yiltirak's protocol
  compensates lower temperature with longer holds, and over a 4 h hold at 100 °C the core
  accumulates thiol far beyond the measurement."* **Zhai's grid is the direct measurement of what
  the model is missing: thiols in a Cys/pentose pot PEAK and then DECLINE, and the peak time
  moves from ≈ 80 min at 140 °C to ≥ 180 min at 120 °C.** A model that reproduces the 120-min
  slice but has no turnover will over-predict at long hold exactly as the exam reports.
- ★ **And it supplies the missing decay constraint** (§8.5), which is the thing the B2.2 fit had
  to invent (`Ea_decay_thiol_sink` = 248.0 kJ mol⁻¹ pressed against its 250 ceiling, F-5
  FALSIFIED).

## 8.5 ★★ NET THIOL DECAY — the first measured turnover in the corpus `[D]`

Treating the post-peak limb as apparent first order (**a lower bound on the sink, because
formation is still running**):

| species | T | window | from → to (μg/L) | **k_net (min⁻¹)** | apparent t₁⁄₂ |
|---|---:|---|---|---:|---:|
| **FFT** | **140 °C** | 120 → 180 min | 11.439 → 4.906 | ★ **0.0141** | **49 min** |
| FFT | 140 °C | 80 → 180 min | 12.439 → 4.906 | 0.0093 | 75 min |
| FFT | 120 °C | 120 → 180 min | 4.107 → 4.279 | **no net decay** | — |
| FFT | 100 °C | 120 → 180 min | 3.734 → 3.119 | 0.0030 | 231 min |
| **MFT** | **140 °C** | 120 → 180 min | 5.907 → 4.012 | **0.0065** | 108 min |
| MFT | 120 °C | 120 → 180 min | 1.388 → 2.479 | **no net decay** | — |
| MFT | 100 °C | 120 → 180 min | 1.237 → 0.923 | 0.0049 | 142 min |

**★ The single most useful number here: FFT loses more than half its headspace concentration in
60 minutes at 140 °C, in a pot where its precursor (furfural) is still increasing
(11.039 → 12.386 μg/L over the same window).** That is a real sink running fast at 140 °C, and it
is incompatible with the B2.2 fit's picture of *"a sink whose absolute rate is tiny across the
whole 100–145 °C window"* (`kinetic_core_b2_2_diagnosis.md` §3).

**Caveat chain, mandatory when quoting these:** (i) these are **net** rates — formation continues,
so the true sink is faster, and the 120 °C "no net decay" cells do **not** mean the sink is off,
only that formation still exceeds it; (ii) headspace concentration is not liquid concentration,
and thiol partitioning is temperature-sensitive; (iii) the interval carries no replicate SD;
(iv) two points cannot establish first order — **the shape is an assumption, and at 140 °C FFT the
three post-peak points (12.439 / 11.439 / 4.906) are visibly NOT single-exponential**, which is
itself informative (an accelerating loss, consistent with a bimolecular sink against a rising
melanoidin/dicarbonyl pool rather than a first-order one).
**⇒ Verdict `USE-Q` as an ORDER-OF-MAGNITUDE FLOOR on the thiol sink at 140 °C, and `STRUCTURAL`
for the claim that the sink is fast at 140 °C and comparable to formation at 120 °C.**

## 8.6 The full derived-Ea inventory `[D]` — column E (Kang's column) vs the other three

**Complete two-leg Eₐ table for the seven best-populated compounds, all four hold times, TTCA
alone**, kJ mol⁻¹, from `Eₐ = R·ln(y₂/y₁)/(1/T₁ − 1/T₂)`:

| compound | 40 min low / high | 80 min low / high | **120 min low / high** | 180 min low / high |
|---|---|---|---|---|
| **MFT** | 104.4 / 64.9 | 37.2 / 79.5 | **7.0 / 97.8** | 60.3 / 32.5 |
| **FFT** | 28.5 / 101.5 | 56.0 / 95.6 | **5.8 / 69.2** | 19.3 / 9.2 |
| 2-Acetylthiazole | 81.0 / 60.6 | 101.5 / 36.5 | **64.0 / 7.7** | 95.7 / 21.5 |
| Furfural | 142.6 / 23.4 | 25.3 / 56.4 | **32.8 / 43.5** | 22.0 / 89.5 |
| 3-Mercapto-2-butanone | — / 236.2 | 110.7 / 56.2 | **19.6 / 19.4** | 20.7 / 7.5 |
| 3-Mercapto-2-pentanone | 121.2 / 21.5 | 45.8 / 63.1 | **17.6 / 19.6** | −0.5 / 20.5 |
| Thiazole | 114.3 / 18.6 | 84.1 / 37.3 | **42.4 / 8.5** | 64.2 / 34.2 |
| **SULFUR CLASS** | 105.2 / 67.8 | 98.0 / 59.4 | **59.4 / 33.1** | 72.9 / 24.9 |
| **O-HETEROCYCLE CLASS** | 130.6 / 34.0 | 29.5 / 68.7 | **36.6 / 49.2** | 29.2 / 84.1 |
| **GRAND TOTAL** | 111.2 / 60.1 | 79.4 / 61.1 | **55.5 / 35.6** | 63.8 / 34.2 |

**Reading:** the **class-level** numbers are far more stable across hold time than any individual
compound — the sulfur class's overall Eₐ (100→140) is **87.4 / 79.7 / 46.9 / 50.1** kJ mol⁻¹ at
40 / 80 / 120 / 180 min, a 1.9× spread, while MFT's high-leg value spans 32.5–97.8, a 3.0× spread
with no monotone trend. **⇒ If anything from this experiment is to be fitted, fit the CLASS, not
the individual thiol; and quote the range across hold times, not one hold time's value.**

**★ The most defensible single number this paper supports:** the **sulfur-class apparent Eₐ over
100–140 °C, averaged over the four hold times = 66 ± 20 kJ mol⁻¹** (values 87.4, 79.7, 46.9,
50.1). Every caveat above applies; it is a *yield-response coefficient*, not a barrier.

---

# §9. THE +CYS ARM — what a 1:1 cysteine addition does `[M]`, `[D]`

The +Cys columns are **entirely new to the corpus** (Kang published only the TTCA-alone column).

| finding | evidence | verdict |
|---|---|---|
| ★ **Cys SUPPRESSES total volatiles early and REVERSES late.** | 100 °C: T−C/T = **0.13 / 0.37 / 0.72 / 2.15** at 40/80/120/180 min. Same crossover at 120 °C (0.26 → 0.55 → 1.16 → 1.61) and 140 °C (0.22 → 0.34 → 0.93 → 1.34). | ★ `STRUCTURAL`, replicated at 3 temperatures |
| ★ **The crossover time moves EARLIER with temperature** — between 120 and 180 min at 100 °C, between 80 and 120 min at 140 °C. | the ratios above | `STRUCTURAL` — this is the authors' central claim and their own grid supports it |
| ★★ **Cys ABOLISHES furfural at 120 and 140 °C.** | Furfural in T−C: 120 °C **0 / 0.348 / 0 / 0**; 140 °C **0.073 / 0 / 0 / 0**, against 3.179–12.386 in the TTCA arm. **Every O-heterocycle is 0 in the +Cys arm at 120 and 140 °C from 120 min onward.** | ★★ `STRUCTURAL`, very strong: the H₂S released by Cys consumes the furan pool essentially completely |
| ★ **and the sulfur it goes into is identifiable.** | 3-Thiophenethiol at 120 °C jumps from **0.178** (T) to **12.892** (T−C) at 120 min — a **72×** increase, the largest single-cell effect in the grid. Thieno[3,2-b]thiophene 1.302 → 7.967 (6.1×). FFT 4.107 → 11.591 (2.8×). | ★ `USE-Q` as a branching constraint: **furan-ring oxygen → sulfur substitution is the dominant fate of the furan pool under excess H₂S** |
| **Cys switches the pyrazine channel on.** | N-heterocycle subtotal, 140 °C, 180 min: **0.940** (T) → **6.286** (T−C), 6.7×. At 100 °C, 180 min: 0.739 → 1.203 (the two numbers the paper quotes in text). | `USE-Q` |
| ⚠ **Cys SUPPRESSES 2-acetylthiazole**, against the trend of every other sulfur compound. | 120 °C, 120 min: 8.795 (T) → 0.538 (T−C), a **16× decrease**; 140 °C: 9.858 → 4.473. | ★ `STRUCTURAL` and **counter-intuitive** — added cysteine *reduces* the thiazole that cysteine chemistry is supposed to make. Worth an explicit note: the model should not assume "more Cys ⇒ more of every S-heterocycle." |

⚠️ **A caution on the +Cys arm's Ea values:** several +Cys cells are `0` at 120 °C but non-zero at
100 and 140 °C (2-acetylthiazole, 3-mercapto-2-butanone, furfural), which produces **negative
apparent Eₐ on the low leg and enormous positive ones on the high leg** (e.g. 2-acetylthiazole
T−C-120: −95.8 then +143.0). These are **non-monotone in temperature**, i.e. unfittable by any
single Arrhenius form. **`REFUSE` every +Cys Eₐ that spans a zero cell.**

---

# §10. FIGURE 5 — the α-dicarbonyl and ARP time courses `[D]` (digitised, coarse)

**Digitisation basis:** page 7 at 600 dpi; panels (a) and (b) read marker-by-marker; panels
(c)–(e) read at 150–300 dpi for peak value and peak time only. **Series were assigned to
temperatures by marker shape where resolvable (□ = 100 °C, ◇ = 120 °C, △ = 140 °C; dashed = TTCA,
solid = TTCA+Cys) and by peak-time ordering where not. Quoted precision ±0.02 mmol/L on panel (b),
±10–15 % elsewhere. Units mmol/L, verified on raster (the y-axis label is misspelled
"Cncentration" in all five panels — the paper's typo, not mine).**

### 10a. Panel (b), 3-DX — the best-resolved panel

| t (min) | T-100 | T-120 | T-140 | T+C-100 | T+C-120 | T+C-140 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 20 | 0.15 | 0.21 | **0.325** | 0.09 | **0 (BDL)** | **0 (BDL)** |
| 40 | 0.21 | 0.32 | 0.49 | 0.105 | **0 (BDL)** | **0 (BDL)** |
| 60 | 0.38 | 0.47 | 0.70 | 0.18 | 0.13 | **0 (BDL)** |
| 80 | 0.47 | 0.69 | **0.81 ← peak** | 0.215 | 0.155 | 0.15 |
| 100 | **0.51 ← peak** | **0.69 ← peak** | 0.79 | 0.23 | 0.18 | **0.30 ← peak** |
| 120 | 0.42 | 0.59 | 0.42 | 0.32 | 0.17 | 0.22 |
| 140 | 0.39 | 0.43 | 0.21 | 0.33 | 0.17 | 0.075 |
| 160 | 0.255 | 0.345 | ~0.01 | 0.25 | 0.11 | 0 |
| 180 | 0.19 | 0.29 | 0 | 0.08 | 0.02 | 0 |

★ **The structural facts, which are robust to the digitisation error:**
1. **3-DX peaks and is consumed at every temperature.** Peak time **≈ 80 min at 140 °C,
   ≈ 100 min at 120 °C, ≈ 100 min at 100 °C**; by 180 min the 140 °C pot has **none left**.
2. ★ **In the +Cys arm 3-DX is BELOW DETECTION for the first 40–60 min at 120 and 140 °C** —
   the authors' own headline (*"3-DX was undetectable before the reaction period of 40−60 min in
   the TTCA−Cys model"*, p. 14307), and it is the cleanest available evidence that **the free
   thiol scavenges the deoxyosone before it can go downstream.**
3. Peak 3-DX is **≈ 0.81 mmol/L against 10 mmol/L TTCA loaded, i.e. ≈ 8 % of the carbon appears
   as free 3-DX at the maximum.**

### 10b. Panels (a), (c), (d), (e) — peak values and times `[D]`, coarse

| panel | species | y-scale | **peak, TTCA arm** | peak time | note |
|---|---|---|---|---|---|
| (a) | **ARP** | 0–2.5 mmol/L | **1.93** (120 °C) · **0.93** (140 °C) · **1.11 and still rising at 180 min** (100 °C) | 100 min (120 °C); **40 min (140 °C)**; > 180 min (100 °C) | ★ classic intermediate: higher T ⇒ earlier, lower peak |
| (c) | **1-DX** | 0–0.2 mmol/L | **≈ 0.11** | ≈ 100 min | an order of magnitude below 3-DX |
| (d) | **MGO** | 0–0.1 mmol/L | **≈ 0.072** | ≈ 60–80 min | |
| (e) | **GO** | 0–0.01 mmol/L | **≈ 0.0057** | ≈ 80–100 min | ★ **GO is ~13× below MGO and ~140× below 3-DX** |

★ **The α-dicarbonyl hierarchy in a Cys/pentose pot, measured on one panel by one lab:**
**3-DX (0.81) ≫ 1-DX (0.11) > MGO (0.072) ≫ GO (0.0057) mmol/L at peak.** That is a
**140-fold span from 3-DX to GO** and it is a directly usable `USE-Q` constraint on the trunk's
dicarbonyl partition — the corpus has very few same-panel dicarbonyl ratios.

⚠️ **The ARP panel's series assignment is the least secure item in this dossier.** The two dashed
curves peaking at 40 min (0.93) and 100 min (1.93) were assigned to 140 °C and 120 °C
respectively **on the physical argument that a higher temperature gives an earlier, lower
intermediate peak**, corroborated by panel (b) where the marker shapes ARE resolvable and give
exactly that ordering. **If a later wave needs the ARP numbers quantitatively, re-verify the
marker glyphs at ≥ 900 dpi.**

---

# §11. FIGURE 6 — ★ THE MECHANISM VALIDATION MODELS, and the only molar yields in the paper

**Conditions:** (NH₄)₂S 10 mmol/L + each substrate at equimolar loading, **pH 7.0, 120 °C**.
**Units: g/L** (raster-verified — three orders above Fig. 3's μg/L scale, because these are neat
model reactions, not a Maillard pot). 3-D perspective bars → **quoted precision ±0.03 g/L.**

| panel | substrate (blue) | 0 min | 10 min | 20 min | product (red) | 0 min | 10 min | 20 min |
|---|---|---:|---:|---:|---|---:|---:|---:|
| (a) | **Furfural** | 1.00 | 0.22 | ≈0.02 | **2-Furfurylthiol (FFT)** | 0.02 | 0.17 | **0.28** |
| (b) | **4-Hydroxy-5-methyl-3(2H)-furanone** | 1.20 | 0.03 | 0.03 | **2-Methyl-3-furanthiol (MFT)** | 0.01 | 0.39 | **0.42** |
| (c) | **2-Methylfuran** | 0.83 | 0.06 | ≈0.01 | **2-Methylthiophene** | 0.01 | 0.15 | **0.16** |
| (d) | **Furan** | 0.69 | ≈0.02 | ≈0.02 | **Thiophene** | 0.01 | 0.09 | **0.25** |

Panel (e), time grid **0 / 30 / 60 min**: **MGO+GO 0.67 → 0.27 → 0.02**; **thiazoles
0.01 → 0.24 → 0.27**; **thiophenes 0.01 → 0.15 → 0.20**; **N-heterocycles 0.01 → 0.05 → 0.08**.

### 11.1 `[D]` **Molar conversion yields — the numbers the repo's furan→thiol edge has never had**

| reaction | substrate consumed | product formed | **molar yield** |
|---|---|---|---|
| **Furfural + H₂S → FFT** | 0.98 g/L ÷ 96.08 = **10.2 mmol/L** | 0.26 g/L ÷ 114.16 = **2.3 mmol/L** | ★ **≈ 22 %** |
| **4-OH-5-Me-furanone + H₂S → MFT** | 1.17 g/L ÷ 114.10 = **10.3 mmol/L** | 0.41 g/L ÷ 114.16 = **3.6 mmol/L** | ★ **≈ 35 %** |
| **2-Methylfuran + H₂S → 2-methylthiophene** | 0.82 g/L ÷ 82.10 = **10.0 mmol/L** | 0.15 g/L ÷ 98.17 = **1.5 mmol/L** | ≈ **15 %** |
| **Furan + H₂S → thiophene** | 0.67 g/L ÷ 68.07 = **9.8 mmol/L** | 0.24 g/L ÷ 84.14 = **2.9 mmol/L** | ≈ **29 %** |

★★ **Four independent measurements of the same reaction class, all at 120 °C, all landing in
15–35 %.** This is the corpus's **first direct constraint on the O→S substitution branch
fraction**, and it says the edge is **efficient but far from quantitative**: three quarters of
the furan carbon goes somewhere else (polymerisation/melanoidin, per the authors).

★ **And the rate is fast:** every furan substrate is **> 78 % consumed in 10 min** and
essentially **exhausted by 20 min at 120 °C** — pseudo-first-order k ≳ 0.15 min⁻¹ against
10 mmol/L sulfide. **The corresponding thiol keeps rising from 10 to 20 min after its substrate
is gone** (FFT 0.17 → 0.28; thiophene 0.09 → 0.25), which means **the thiol's immediate precursor
is not the parent furan itself but an intermediate**, exactly as the authors' Fig. 3d scheme
draws it (furfural → [H₂S addition] → ... → FFT).

**Caveats:** (i) these are neat two-component systems at ~10 mmol/L, three orders of magnitude
above the μg/L levels in the Maillard pot — **the branch fraction may not transfer**; (ii) the
substrate is (NH₄)₂S, a sulfide *source*, not the H₂S released by cysteine degradation, and its
free-H₂S activity at pH 7 is not stated; (iii) 3-D bars, no error bars, n unstated for this
experiment; (iv) the Methods/figure time-grid disagreement (§3.3). **Verdict `USE-Q` for the
15–35 % branch fraction as a PRIOR with a factor-2 interval; `STRUCTURAL` for "the furan pool is
consumed within 20 min at 120 °C under excess sulfide."**

---

# §12. VERIFIED NEGATIVES `[NEG]` — do not re-open these

| # | what a reader might hope for | verdict |
|---|---|---|
| 1 | A seven-temperature MFT/FFT ladder | ★ **ABSENT.** Seven temperatures exist for **A₄₂₀ only** (§0.1) |
| 2 | Any rate constant, half-life, Eₐ, or reaction order printed by the authors | **ABSENT.** The paper contains **no kinetic model and no fitted parameter of any kind**, despite calling pH/T/Cys-dose "kinetic parameters" |
| 3 | Standard deviations on the volatile grid | **ABSENT** (§6.4). Methods promise mean ± SD; Fig. 3 prints bare means |
| 4 | Measured pH after heating, or any buffer | **ABSENT** (§3.2). Unbuffered, drift unmeasured |
| 5 | H₂S measured directly | **ABSENT.** H₂S is inferred throughout and **surrogated** by (NH₄)₂S in Fig. 6. No H₂S concentration anywhere |
| 6 | Free cysteine or TTCA concentration vs time | **ABSENT.** Neither the substrate nor the released Cys is quantified — only ARP, the two DXs, MGO, GO. (Kang 2026 Figs. S3/S4 do carry these, for the TTCA-alone system) |
| 7 | TTCA purity | **ABSENT** (§3.1) |
| 8 | HPLC gradient programmes | **ABSENT**, twice, explicitly (§4.3) |
| 9 | LOD/LOQ for any analyte | **ABSENT** |
| 10 | The A₄₂₀ dilution factor | **ABSENT** — so A₄₂₀ magnitudes are not transferable between panels (§4.3) |
| 11 | Sensory data of any kind | **ABSENT** |
| 12 | Melanoidin or high-MW characterisation | **ABSENT** — browning is A₄₂₀ only |
| 13 | Any measurement above 150 °C or below 90 °C | **ABSENT** |

---

# §13. CONSOLIDATED PARAMETER TABLE

**Common conditions unless stated:** TTCA 10 mmol/L in water, **unbuffered**, pH set to 7.0,
sealed pressure-resistant bottle with magnetic stirring, ice-quench, HS-SPME (CAR/PDMS/DVB 75 μm,
60 °C/20 min) GC-MS on DB-Wax, IS = 1,2-dichlorobenzene, external calibration, n = 3, **μg/L**.

| # | parameter | value | units | condition | prov | anchor |
|---:|---|---|---|---|---|---|
| **— the grid —** |
| 1 | **Full 39 × 8 volatile grid, 100 °C** | see §6a | μg/L | 4 hold times × ±Cys | **M** | Fig. 3a |
| 2 | **Full 39 × 8 volatile grid, 120 °C** | see §6b | μg/L | " | **M** | Fig. 3b |
| 3 | **Full 39 × 8 volatile grid, 140 °C** | see §6c | μg/L | " | **M** | Fig. 3c |
| 4 | ⚠ of which column E (T-120 min) | **already published as Kang Table S4** | — | — | **M** | §7 |
| **— thiols, TTCA arm —** |
| 5 | MFT, 100/120/140 °C, 40 min | 0.078 / 0.432 / 1.129 | μg/L | | **M** | Fig. 3 |
| 6 | MFT, 80 min | 0.509 / 0.937 / 3.043 | μg/L | | **M** | " |
| 7 | MFT, 120 min | 1.237 / 1.388 / 5.907 | μg/L | ⚠ = Kang S4 | **M** | " |
| 8 | MFT, 180 min | 0.923 / 2.479 / 4.012 | μg/L | | **M** | " |
| 9 | FFT, 40 min | 0.708 / 1.129 / 5.079 | μg/L | | **M** | " |
| 10 | FFT, 80 min | 1.206 / 3.019 / **12.439** | μg/L | ★ FFT maximum of the whole grid, TTCA arm | **M** | " |
| 11 | FFT, 120 min | 3.734 / 4.107 / 11.439 | μg/L | ⚠ = Kang S4 | **M** | " |
| 12 | FFT, 180 min | 3.119 / 4.279 / 4.906 | μg/L | | **M** | " |
| **— derived kinetics —** |
| 13 | **MFT apparent Eₐ, 100→140** | **85.6 / 57.3 / 50.1 / 47.1** | kJ mol⁻¹ | at 40 / 80 / 120 / 180 min | **D** | §8.2 |
| 14 | **FFT apparent Eₐ, 100→140** | **63.1 / 74.8 / 35.9 / 14.5** | kJ mol⁻¹ | " | **D** | §8.2 |
| 15 | **Sulfur-class apparent Eₐ, 100→140** | **87.4 / 79.7 / 46.9 / 50.1**; mean **66 ± 20** | kJ mol⁻¹ | " | **D** | §8.6 |
| 16 | **O-heterocycle-class Eₐ, 100→140** | 84.8 / 48.1 / 42.6 / 55.2 | kJ mol⁻¹ | " | **D** | §8.6 |
| 17 | ★ **k_net, FFT decay, 140 °C** | **0.0141** | min⁻¹ | 120→180 min, net of formation | **D** | §8.5 |
| 18 | ★ **k_net, MFT decay, 140 °C** | **0.0065** | min⁻¹ | " | **D** | §8.5 |
| 19 | k_net, FFT decay, 100 °C | 0.0030 | min⁻¹ | " | **D** | §8.5 |
| 20 | FFT peak time, 140 °C | **≈ 80** | min | TTCA arm | **D** | §8.3 |
| 21 | MFT peak time, 140 °C | ≈ 120 | min | " | **D** | §8.3 |
| **— browning —** |
| 22 | A₄₂₀, TTCA, 90–150 °C | 0.24 / 0.25 / 0.346 / 0.497 / 0.572 / 0.629 / 0.867 | — | 120 min, pH 7 | **D** (digitised) | Fig. 1b |
| 23 | A₄₂₀, TTCA−Cys 1:1, 90–150 °C | 0.22 / 0.23 / 0.267 / 0.274 / 0.361 / 0.418 / 0.463 | — | " | **D** | Fig. 1b |
| 24 | **Browning apparent Eₐ, TTCA** | **28.3** (R² 0.969; legs 28.8 / 25.7) | kJ mol⁻¹ | 7 rungs, 120 min | **D** | §5.1 |
| 25 | **Browning apparent Eₐ, TTCA−Cys** | **16.9** (R² 0.948; legs 8.7 / 24.2) | kJ mol⁻¹ | " | **D** | §5.1 |
| 26 | Cys colour-inhibition, dose optimum | **1 mol Cys per mol TTCA** (2:1 no better) | — | 120 °C, 120 min, pH 7 | **M** | Fig. 1a |
| 27 | A₄₂₀ vs pH, **TTCA** | **0.50⚠ / 0.482 / 0.41⚠ / 0.503 / 0.550 / 0.827 / 0.940** | — | pH 5.5 / 6 / 6.5 / 7 / 7.5 / 8 / 8.5; 120 °C, 120 min | **D** (digitised, ±0.02; axis calibrated at 735 px per A₄₂₀ unit) | Fig. 1c |
| 27b | A₄₂₀ vs pH, **TTCA−Cys 1:1** | **0.455⚠ / 0.408 / 0.36⚠ / 0.397 / 0.415 / 0.496 / 0.603** | — | " | **D** | Fig. 1c |
| 27c | ★ shape of the pH response | **U-shaped with a minimum at pH 6.5, flat 5.5–7.5, then a sharp rise above 7.5** — TTCA rises **2.0×** from pH 7.5 to 8.5, TTCA−Cys only **1.45×** | — | " | **D** | §14 `PRIOR-ONLY` |
| **— dicarbonyls —** |
| 28 | **3-DX peak** | **0.81 / 0.69 / 0.51** | mmol/L | 140 / 120 / 100 °C, TTCA arm | **D** | §10a |
| 29 | 3-DX peak time | **80 / 100 / 100** | min | 140 / 120 / 100 °C | **D** | §10a |
| 30 | ★ **dicarbonyl hierarchy at peak** | 3-DX 0.81 ≫ 1-DX 0.11 > MGO 0.072 ≫ GO 0.0057 | mmol/L | 140 °C, TTCA arm | **D** | §10b |
| 31 | ARP peak | 1.93 (120 °C, 100 min); 0.93 (140 °C, 40 min) | mmol/L | | **D** | §10b |
| 32 | ★ 3-DX below detection in +Cys arm | **first 40–60 min** at 120 and 140 °C | — | | **M** | p. 14307 |
| **— the H₂S branch (Fig. 6) —** |
| 33 | ★ **furfural → FFT molar yield** | **≈ 22 %** | — | (NH₄)₂S 10 mM, 120 °C, 20 min | **D** | §11.1 |
| 34 | ★ **4-OH-5-Me-furanone → MFT molar yield** | **≈ 35 %** | — | " | **D** | §11.1 |
| 35 | 2-methylfuran → 2-methylthiophene yield | ≈ 15 % | — | " | **D** | §11.1 |
| 36 | furan → thiophene yield | ≈ 29 % | — | " | **D** | §11.1 |
| 37 | furan-pool consumption rate | > 78 % in 10 min; complete by 20 min | — | " | **M** | Fig. 6 |
| 38 | MGO+GO consumption under sulfide | 0.67 → 0.27 → 0.02 g/L at 0/30/60 min | g/L | 120 °C | **M** | Fig. 6e |
| **— branching ratios —** |
| 39 | **FFT/MFT, TTCA arm, 120 min** | **3.02 / 2.96 / 1.94** | — | 100/120/140 °C ⚠ = Kang's | **D** | §8 |
| 40 | ★ **FFT/MFT at other hold times** | 40 min: **9.08 / 2.61 / 4.50**; 80 min: 2.37 / 3.22 / 4.09; 180 min: **3.38 / 1.73 / 1.22** | — | 100/120/140 °C | **D** | ★ **new** |
| 41 | ★ **FFT/MFT, +Cys arm** | 1.49–5.57 across the grid | — | | **D** | ★ **new** |

★ **Row 40 is important on its own:** the corpus has been carrying **FFT:MFT ≈ 3.0** as a
Cys/pentose branching constant (K5's three-axis characterisation). **Across this paper's full
grid the ratio runs 1.22 – 9.08, a 7.4× span, driven by HOLD TIME as much as by temperature.**
The "≈ 3.0, robust to pH, weakly sensitive to T" statement was derived from the 120-min slice and
does not survive the time axis either.

---

# §14. USABILITY VERDICTS

### `USE` — ingestible as measured
- The **936-value grid** (§6a–c) as measured headspace concentrations, **excluding column E**,
  with the ±15 %/×÷3 uncertainty of §6.4 and the headspace caveat of §4.2.
- The **column totals** (7 of 8 confirmed against the authors' own text).

### `USE-Q` — usable with the stated qualification
| item | qualification |
|---|---|
| Sulfur-class apparent Eₐ **66 ± 20 kJ mol⁻¹** | yield-response coefficient over 100–140 °C, not a barrier; quote the across-hold-time range |
| **k_net(FFT, 140 °C) = 0.0141 min⁻¹** | **lower bound** on the thiol sink; net of ongoing formation; headspace basis; shape assumed |
| Browning Eₐ **28.3 (TTCA) / 16.9 (TTCA−Cys)** | ±8 kJ mol⁻¹; dilution factor unreported. ★ **Independently replicated: Wang 2026's five-rung buffered A₄₂₀ ladder in a different lab and a different pot gives 30.6 kJ mol⁻¹ (R² 0.971) — see `wang2026_extraction.md` §3.2** |
| **Furan→thiol molar yields 15–35 %** | neat 10 mM model, three decades above pot concentrations; factor-2 interval |
| Dicarbonyl hierarchy 3-DX ≫ 1-DX > MGO ≫ GO | digitised, ±15 %; same-panel so the RATIOS are the trustworthy part |

### `RATIO-ONLY`
- Every **+Cys / TTCA ratio** (§9). The absolute +Cys values inherit an unquantified matrix effect
  (added Cys changes the headspace matrix as well as the chemistry), but the ratios are within-run.
- **FFT/MFT** at any fixed (T, time) cell.

### `STRUCTURAL` — shape constraints with no number to fit
1. ★ **Thiols in a TTCA/Cys pot PEAK and then DECLINE, and the peak moves earlier with
   temperature** (≈ 80 min at 140 °C for FFT).
2. ★ **Added Cys suppresses total volatiles early and enhances them late**, with the crossover
   moving earlier as temperature rises.
3. ★ **Excess H₂S abolishes the furan pool** (furfural = 0 in the +Cys arm at 120/140 °C).
4. ★ **Added Cys DECREASES 2-acetylthiazole** while increasing every other S-heterocycle.
5. **3-DX is scavenged to below detection for the first 40–60 min when free thiol is present.**
6. **Pyrazines require exogenous amino nitrogen** — 0 in the TTCA-only arm at 40–120 min at every
   temperature, non-zero as soon as Cys is added (corroborates the same finding in Kang S4/S5).

### `PRIOR-ONLY`
- The pH-browning curve (Fig. 1c): U-shaped with a minimum at pH 6.5 and a sharp rise above 7.5;
  digitised, browning only, one temperature.

### ★ `REFUSE` — with the reason recorded so a later wave does not re-ingest
| item | reason |
|---|---|
| **`63.035 μg/L`** (the paper's 140 °C/T/180 min total) | **arithmetic error in the text**; the figure gives 65.035 (§6.4) |
| **Any Kang-derived "switch-on" Eₐ as a physical barrier** (7.0 → 97.8 MFT; 5.8 → 69.2 FFT) | **hold-time artefact**; not reproduced at 40, 80 or 180 min in the same experiment (§8.3) |
| **Every +Cys Eₐ that spans a zero cell** | non-monotone in T; unfittable by any Arrhenius form (§9) |
| **Absolute A₄₂₀ magnitudes across panels** | dilution factor unreported |
| **Fig. 3 column E as a NEW anchor** | ★ **already in the repo as Kang Table S4 — declaring it again would fit the model on its own hold-out** (§7.3) |

---

# §15. DRAFT FIT / HOLD-OUT ROLES — **FOR THE ORCHESTRATOR, NOT A DECLARATION EDIT**

*(No declaration file was opened or edited by this dossier.)*

### 15a. ★ FIRST, A CORRECTION TO EXISTING ROWS
- **`kang_switch_on_*` (both rows): recommend RETIREMENT, not re-scoring.** The feature is not a
  property of the chemistry (§8.3). Retiring a hold-out on evidence about the *source* rather
  than about model performance is exactly the case the `k3` §D rules permit.
- **`kang_140C_MFT` / `kang_140C_FFT`: keep, but re-label the provenance to
  `zhai2023jafc_fig3_colE` and record that Kang is not an independent source.**

### 15b. Recommended **FIT** rows (all genuinely unpublished)
| candidate | why |
|---|---|
| **FFT and MFT full time courses at 100 and 120 °C, TTCA arm** (8 points) | constrains formation on the two rungs the pH state is least wrong about; contains no 140 °C information |
| **Furan→thiol molar yields (4 reactions, 120 °C)** | the O→S branch fraction, currently unconstrained |
| **3-DX peak height and peak time at 120 °C** | trunk dicarbonyl |

### 15c. Recommended **HOLD-OUT** rows
| candidate | why it is a real test |
|---|---|
| ★ **FFT at 140 °C, four hold times (5.079 / 12.439 / 11.439 / 4.906)** | **the turnover the model has never been asked to reproduce.** A model that only makes thiol will fail it; the current core would predict a monotone rise |
| ★ **k_net(FFT, 140 °C) ≥ 0.01 min⁻¹** as a one-sided test | directly contradicts the B2.2 fit's 248 kJ mol⁻¹ near-inert sink |
| **The +Cys crossover time at each temperature** | tests whether the model's thiol-scavenging channel has the right *timing*, not just the right sign |
| **Furfural = 0 in the +Cys arm at 120 and 140 °C** | a hard structural test of the H₂S–furan edge |
| **Dicarbonyl hierarchy 3-DX ≫ MGO ≫ GO (140×)** | trunk partition |

### 15d. **DO NOT USE**
- Column E at any temperature (already seen, §7.3).
- The 63.035 total (§6.4).
- Any +Cys Eₐ spanning a zero.

---

# §16. DECLARED GAPS FROM THIS PAPER

| # | gap | what would close it |
|---|---|---|
| G1 | **Table S2 — the calibration curves** — is not on disk, so the Tier A/Tier B split cannot be reproduced and the 39-compound grid has no per-compound uncertainty basis | order the SI of `10.1021/acs.jafc.3c04166` |
| G2 | The **MGO−GO−Cys model's flavour table** (SI) is the only quantitative link between the dicarbonyl pool and the pyrazine channel in this paper | same SI |
| G3 | **No free-Cys or TTCA depletion data** here; Kang's Figs. S3/S4 have them for the TTCA-alone system, but **not for the +Cys system**, which is the arm this paper is about | not available anywhere in the corpus |
| G4 | **No SDs on the volatile grid**, so all 936 values carry an assumed uncertainty | SI may or may not have them |
| G5 | **No measured pH after heating** in an unbuffered pot that the B2.2 diagnosis shows is the model's weakest state variable | not obtainable from this paper |
| G6 | **The 90/100/110/130/150 °C rungs exist only for A₄₂₀** — the volatile ladder cannot be extended | would require new experiments |
| G7 | ★ **The independence of the corpus's sulfur ladder is now a live question.** Before any further weighting, someone should check whether Kang 2026 cites this paper as the source of Table S4 | read Kang's reference list against §7 |
