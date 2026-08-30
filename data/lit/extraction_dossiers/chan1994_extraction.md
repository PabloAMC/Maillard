# Chan & Reineccius 1994 — kinetics of methional, dimethyl disulfide and 2-acetylthiophene

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THIS IS THE CORPUS'S ONLY AQUEOUS 75–115 °C SULFUR TEMPERATURE LADDER, AND IT IS THE BASELINE REGIME DIRECTLY BELOW KANG 2026's 120 °C ANCHOR.

### ★ HEADLINE: THE 75–115 °C LEGS ARE **NOT** ARRHENIUS-LINEAR. Every one of the six compound × pH ladders contains a leg-to-leg slope break of ≥ 2× — and for **dimethyl disulfide at pH 6 the break is a switch-on of exactly Kang's shape (×1.30, ×1.23, then ×4.34) sitting 25 K lower in temperature.**

### ⚠ HEADLINE DEFECT: the paper's own Conclusions say methional's Ea is **"17–19 kcal/mole"**. Its own Table I says **19.39–27.07 kcal/mole**. The claim contradicts the table on the facing pages. Raster-verified both. Do not propagate the abstract/conclusion figure.

---

## §0. IDENTITY

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/chan1994.pdf` — **1,143,069 bytes**, 11 pages, PDF 1.6 | `ls`, `pdfinfo` |
| SHA-256 | `3632e755c25b602e9042b052983fc478ebff940c5afcd4917442ff364006cb65` | `shasum -a 256` |
| **title** | ***"Kinetics of the Formation of Methional, Dimethyl Disulfide, and 2-Acetylthiophene via the Maillard Reaction"*** | p.127 raster ✔ (title-page order is *Methional, Dimethyl Disulfide, 2-Acetylthiophene*) |
| **authors** | ★ **F. Chan and G. A. Reineccius** | p.127 raster ✔ |
| **affiliation** | **Department of Food Science and Nutrition, University of Minnesota, 1334 Eckles Avenue, St. Paul, MN 55108** (single affiliation, both authors) | p.127 raster ✔ |
| **DOI** | ★ **`10.1021/bk-1994-0564.ch010`** | page-footer stamp on every page |
| book | *Sulfur Compounds in Foods*, Mussinan, C. J.; Keelan, M. E., Eds. | running footer, all pages |
| series | **ACS Symposium Series 564** — **Chapter 10, pp. 127–137** (11 pp.) | p.127 header + page numbers |
| publisher | American Chemical Society, Washington, DC, 1994 | p.127 footer |
| **publication date** | ★ **July 29, 1994** | page-footer stamp |
| received | **RECEIVED April 27, 1994** | p.137 |
| copyright line | `0097-6156/94/0564-0127$08.00/0 · © 1994 American Chemical Society` | p.127 |
| version | published book chapter, version of record (scanned) | — |
| provenance stamp | *"Downloaded by FUDAN UNIV on January 13, 2017 | http://pubs.acs.org"* | every page |
| **funding** | ★ **[NEG] none stated** — there is no acknowledgements section at all | whole-document sweep |
| **conflicts** | **[NEG] none stated** | — |
| CAS registry (printed by the authors) | dimethyl disulfide **624-92-0**; methional **3268-49-3**; 2-acetylthiophene **88-15-3**; 4-methylthiazole **693-95-8** | p.137 raster ✔ — all four are the correct CAS, an independent check that the OCR digits are clean here |

### ⚠ OCR-LAYER WARNING (read before reusing any number from this dossier)

The PDF carries an **ABBYY FineReader 9.0 Corporate Edition** text layer (later touched by `iText 4.2.0`).
The text layer is therefore machine-read, not born-digital. **Every kinetic number, every Table I cell,
every axis label and every schedule entry quoted below was independently confirmed against a
`pdftoppm -r 600` raster render.** The results of that confirmation are in §0.1.

### §0.1 TEXT-LAYER vs RASTER DEFECT TABLE

| # | item | text layer | raster (600 dpi) | verdict |
|---|---|---|---|---|
| D1 | Table I — all 18 numeric cells (2 replicates × 3 pH × r²/Ea) | as printed below | **identical, all 18** | ✔ **no OCR error anywhere in Table I** |
| D2 | CAS numbers ×4 | 624-92-0 / 3268-49-3 / 88-15-3 / 693-95-8 | identical | ✔ clean |
| D3 | molar loadings (5 reagents, mol + g) | 0.5/90; 0.075/11.19; 0.075/12.39; 0.075/8.64; 0.075/9.84 | identical | ✔ clean; and all five mol↔g pairs are internally self-consistent against the true MWs (§1.1) |
| D4 | sampling schedule (5 rows) | 75 °C 1.5–7.5 h … 115 °C 5–25 min | identical | ✔ clean — **the 95 °C row really does read `20, 40, 60, 100, 120 min.`** |
| D5 | `r²` ranges in prose | "ranged form 0.83 to 0.99" (p.131); "ranged from 0.81 to 0.98" (p.135) | identical, **including the authors' own typo `form`** | ✔ clean OCR; ⚠ but see C3 |
| D6 | "activation energies ranged from 15 to 31 kcal/mole" | as shown | identical | ✔ clean OCR; ⚠ but see C2 |
| D7 | "The Ea of methional is moderate (17-19 kcal/mole)" | as shown | identical | ✔ clean OCR; ⚠ but see **C1 — this is a real contradiction, not an OCR artefact** |
| D8 | Greek/symbol glyphs | `1 μιη`, `30 m χ 0.32 mm`, `1.0 μΐ,`, `scanningfromm/z`, `cc-aminoketone`, `datafitthe` | `1 µm`, `30 m × 0.32 mm`, `1.0 µL`, `scanning from m/z`, `α-aminoketone`, `data fit the` | ⚠ **OCR mangling confined to Greek letters, `×`, and inter-word spaces. No digit was ever affected.** |

**⇒ Net OCR risk for this paper: LOW. The scan is clean on digits. The real hazards in this paper are
the authors' own internal contradictions (C1–C4), not the OCR.**

### §0.2 THE PAPER'S OWN INTERNAL CONTRADICTIONS `[D]`

| # | contradiction | evidence | severity |
|---|---|---|---|
| **C1** | ★ Conclusions (p.135): *"The Ea of methional is moderate (**17-19 kcal/mole**)"*. **Table I (p.136) prints methional Ea = 27.07, 24.95, 19.39, 21.74, 26.02, 19.50 kcal/mol — range 19.39–27.07.** Only the two pH-8 values fall near 19; nothing is 17. | both raster-verified | ★★★ **the headline number of the paper is wrong against its own table** |
| **C2** | Conclusions: *"activation energies ranged from **15 to 31** kcal/mole"*. Table I range is **15.02 to 32.71** kcal/mol. The upper bound excludes the 2-acetylthiophene pH-8 1st-replicate value (32.71). | raster ✔ | ★★ |
| **C3** | Two different `r²` ranges are quoted for the same fits: **"0.83 to 0.99"** (p.131, describing the zero-order concentration-vs-time fits) and **"0.81 to 0.98"** (p.135, describing "the fit of the data to zero order kinetics"). Table I's `r²` column runs **0.797–0.976**, which matches *neither*. It is never stated whether Table I's `r²` is the zero-order fit or the Arrhenius fit. | raster ✔ | ★★ **Table I's `r²` column is of undeclared meaning** (see §4.2 for the resolution) |
| **C4** | ★ **The Methods sampling schedule for 95 °C is `20, 40, 60, 100, 120 min`. All three 95 °C figures (Figs 1, 2, 3) label their x-axis `20, 40, 60, 80, 100 min`.** | Methods raster ✔; Fig-1/2/3 axis ticks digitised at 600 dpi ✔ | ★★★ **the time axis of every 95 °C data figure disagrees with the stated protocol** |
| **C5** | Fig 1 (DMDS ppm vs time at 95 °C) and Fig 4 (DMDS `k` at 95 °C) do not reconcile at pH 7 and pH 8. Fig-1 slopes give ≈ 0.21 (pH 7) and ≈ 0.33 (pH 8) ppm/min, with pH 8 > pH 7 at every time point; Fig 4 gives 0.50 (pH 7) > 0.37 (pH 8). pH 6 reconciles (0.12 vs 0.12). | §3.2/§3.5 | ★★ resolved-in-part by §4.3: the figures are **different replicates**, but the paper never says so |

---

## §1. SYSTEM AND CONDITIONS `[M]`

| item | value (verbatim units) | anchor |
|---|---|---|
| **sugar** | **glucose, 0.5 mole = 90 g** (Aldrich Chemical, Milwaukee, WI) | p.128 raster ✔ |
| **amino acid 1** | **L-methionine, 0.075 mole = 11.19 g** (Sigma Chemical, St. Louis, MO) | p.128 ✔ |
| **amino acid 2** | **phenylalanine, 0.075 mole = 12.39 g** (Sigma) | p.128 ✔ |
| **amino acid 3** | **proline, 0.075 mole = 8.64 g** (Aldrich) | p.128 ✔ |
| **amino acid 4** | **leucine, 0.075 mole = 9.84 g** (Aldrich) | p.128 ✔ |
| **solvent / volume** | **400 mL distilled water** (Glenwood Inglewood, Minneapolis, MN) | p.128 ✔ |
| **pH setting** | ★ *"the pH was adjusted to **6, 7, or 8** with a **0.1M phosphate buffer** and the appropriate amount of **NaOH**"* | p.128 ✔ |
| **buffer molarity** | ★ **0.1 M phosphate — the molarity IS given.** See ⚠ below for why it is nevertheless nominal. | p.128 ✔ |
| **vessel** | **600 mL pressure reactor**, sealed, stirred | p.128 ✔ |
| **heating apparatus** | **Parr 4563 reactor with Parr 4842 temperature controller** (Parr Instrument Co., Moline, IL) | p.128 ✔ |
| **come-up handling** | *"Initial temperature was noted and **start time determined when the solutions reached reaction temperatures**"* — i.e. `t = 0` is defined at set-point, come-up time neither reported nor corrected for | p.128 ✔ |
| **temperatures** | ★ **75, 85, 95, 105, 115 °C** (5 levels) | p.128 ✔ |
| **times** | **75 °C**: 1.5, 3.0, 4.5, 6.0, 7.5 h · **85 °C**: 0.5, 1.0, 1.5, 2.0, 2.5 h · **95 °C**: 20, 40, 60, 100, 120 min · **105 °C**: 10, 20, 30, 40, 50 min · **115 °C**: 5, 10, 15, 20, 25 min | p.128 ✔ ⚠ see C4 |
| **sampling** | **50 mL aliquots** withdrawn through a sampling port; **ultra-high-purity N₂ (ca. 50 mL) added after each sampling to restore reactor pressure** | p.128 ✔ |
| **quench** | *"Each sample was **cooled to room temperature in an ice bath**"* | p.128 ✔ |
| **pH monitoring** | ★ *"the **pH measured**"* on every cooled sample — **and the values are never reported anywhere in the chapter** | p.128 ✔ ⚠ `[NEG]` |
| **extraction** | **3 × 5 mL dichloromethane** (EM Science) containing **500 ppm 4-methylthiazole** as internal standard, in a **150 mL separatory funnel**; dried over anhydrous MgSO₄, filtered, **concentrated to 0.5–1.0 mL under high-purity N₂** | p.128 ✔ |
| **replication (n)** | ★ **n = 2.** Table I is printed as `1ST REPLICATE` / `2ND REPLICATE` — two independent runs, never averaged in the table, **and no SD, SE or confidence interval is printed for any value** | p.136 ✔ |

### §1.1 Derived checks on the loadings `[D]`

| check | result |
|---|---|
| glucose 0.5 mol → 90.08 g (MW 180.16) | ✔ paper says 90 g |
| methionine 0.075 mol → 11.19 g (MW 149.21) | ✔ exact |
| phenylalanine 0.075 mol → 12.39 g (MW 165.19) | ✔ exact |
| proline 0.075 mol → 8.64 g (MW 115.13) | ✔ exact |
| leucine 0.075 mol → 9.84 g (MW 131.17) | ✔ exact |
| **⇒ the loadings are absolute moles, not molarities** | **glucose ≈ 1.25 M, each amino acid ≈ 0.1875 M, total amino acid ≈ 0.75 M in 400 mL** |
| **⇒ the only sulfur source is methionine** | **[M] methionine = 0.075 mol = 187.5 mM.** Phe, Pro, Leu contribute no S |

> ⚠ **THE BUFFER IS NOMINAL, NOT CONTROLLING.** 0.1 M phosphate is being asked to hold pH against
> **0.75 M amino acid + 1.25 M glucose** — a ~13-fold molar excess of titratable amine over buffer,
> in a system whose whole point is that it consumes amine and generates acids. The authors measured
> the pH of every sample and reported none of it. **Treat "pH 6 / 7 / 8" as the *initial* pH of a
> drifting, effectively unbuffered system, not as a maintained condition.** Any pH-resolved parameter
> lifted from this paper inherits that.

> ⚠ **THE REACTOR IS PROGRESSIVELY EMPTIED AND N₂-PURGED.** Five 50 mL aliquots are drawn from a
> 400 mL charge — **250 mL, i.e. 62.5 % of the charge, is removed over the run** — and ca. 50 mL of
> N₂ is injected after each withdrawal. Concentrations are not diluted by this, but (a) the
> liquid : headspace ratio changes monotonically along every time series, so volatile partitioning
> changes along the very axis the rate constant is fitted on, and (b) **DMDS formation is an
> *oxidative* coupling of two methanethiols, and the headspace is being repeatedly purged with N₂ of
> unquantified, monotonically increasing extent.** Neither effect is measured or corrected. This is
> a mechanism-relevant, unquantified confounder specific to the disulfide lane.

---

## §2. ANALYTICAL METHOD AND QUANTIFICATION BASIS

### §2.1 Identification (GC/MS) `[M]`

| item | value | anchor |
|---|---|---|
| instrument | **HP 5890 GC + HP 5970 MSD**, capillary interfaced directly into the MSD | p.128–129 |
| ionisation | **70 eV EI**, ion source **220 °C**, scan threshold **750**, **m/z 29–400 at 0.86 s/cycle** | p.129 |
| column | **30 m × 0.32 mm i.d. × 1 µm DB-5** fused silica (J & W Scientific) | p.129 |
| oven | 40 °C (3 min) → **5 °C/min** → 250 °C (10 min); injector **230 °C** | p.129 |
| carrier / inlet | helium, **7 psi head pressure**, **20 : 1 split**, **1.0 µL** injected | p.129 |
| **identification basis** | ★ **three-way**: NBS mass-spectral library match + **published retention indices** + **co-chromatography with authentic compounds** | p.129 |

### §2.2 Quantification (GC/AED) `[M]` — this is the detector that produced the kinetic data

| item | value | anchor |
|---|---|---|
| instrument | ★ **HP 5890 Series II GC + HP 5921A Atomic Emission Detector** — *"were used to collect kinetic data"* | p.129 |
| channels | **carbon 193 nm, nitrogen 174 nm, sulfur 181 nm**; O₂ and H₂ as scavenger gases; make-up flow **30 mL/min** | p.129 |
| reagent-gas flows | *"based on (11)"* — **Fox, L.; Wylie, P. Hewlett-Packard Appl. Note 1989, 228-275, pp 1-3** (i.e. deferred to an HP application note) | p.129, ref 11 |
| column | same 30 m × 0.32 mm × 1 µm DB-5 | p.129 |
| oven | 40 °C (3 min) → 5 °C/min → 250 °C (10 min); injector **250 °C**; AED transfer line **275 °C** | p.129 |
| carrier / inlet | helium ca. **2 mL/min**, **15 psi**, **45 : 1 split**, **1.0 µL** | p.129 |
| **internal standard** | ★ **4-methylthiazole, 500 ppm, dissolved in the dichloromethane extraction solvent** (CAS 693-95-8; one S atom, one N atom) | p.128 |
| **calibration curve** | ★ **[NEG] NONE. No calibration curve, no concentration range, no authentic-standard response factor, and no recovery figure is reported anywhere in the chapter.** | whole-document sweep |
| **AED validation** | deferred to ref **(14) Mistry, B. S.; Reineccius, G. A.; Jasper, B.**, *Sulfur Compounds in Foods*, ACS Symp. Ser., **"in press"** at the time | p.137 |
| authors' own disclaimer | > *"While analytical aspects of this study were not discussed in this paper, the usefulness of the AED as a sensitive and selective detector in flavor research was appreciated (14)."* (p.135) | raster ✔ |

### §2.3 ★ QUANTIFICATION VERDICT

The reported quantity is **"Amount/ppm"** (Figs 1–3). The detector is a **sulfur-selective atomic
emission detector**; the internal standard is **a single sulfur-containing compound at a fixed
500 ppm in the extraction solvent**. No calibration curve, no measured response factor, and no
statement of the equimolar-sulfur assumption appears anywhere.

> ### **THESE DATA ARE SINGLE-INTERNAL-STANDARD SEMI-QUANTITATIVE. THEY ARE USABLE FOR (b) WITHIN-STUDY SHAPE ONLY, AND ARE **NOT** USABLE FOR (a) ABSOLUTE CONCENTRATIONS.**

The reasons, spelled out so a later wave does not re-open this:

1. **No calibration.** Not a single response factor is measured or quoted. The ppm scale rests on an
   undeclared assumption — almost certainly the AED's equimolar-sulfur response — that the authors
   never state, and whose validation they defer to an unpublished companion (ref 14).
2. ★ **The ppm basis is undefined.** The chapter never says whether "ppm" is in the aqueous reaction
   medium or in the **0.5–1.0 mL concentrated dichloromethane extract** derived from a 50 mL aliquot.
   Those two readings differ by a factor of **50–100×**. **No absolute aqueous concentration can be
   recovered from this paper at all.**
3. **No recovery.** A 3 × 5 mL DCM extraction of 50 mL of a 1.25 M sugar syrup, dried, filtered and
   blown down under N₂, has an unmeasured and compound-dependent recovery — and DMDS
   (b.p. 110 °C) and methional (b.p. 165 °C) do not lose solvent equally under a nitrogen stream.

### §2.4 ★ THE WAVE-K6a RULE APPLIES, AND THIS LADDER IS Ea-USABLE

**A constant unknown response factor cancels in a ratio, and therefore cancels in an Arrhenius
slope.** If the measured signal is `S(T) = f · k(T)` with `f` unknown but **constant across the
temperature ladder**, then `ln S(T₂) − ln S(T₁) = ln k(T₂) − ln k(T₁)`, and the fitted
`Ea = −R · d(ln S)/d(1/T)` is **exactly the true Ea**. `f` appears only in the intercept, i.e. only
in the pre-exponential.

**⇒ Chan's semi-quantitative ladder IS usable for apparent activation energies, per-leg and
whole-ladder, and for `k(T₂)/k(T₁)` ratios and Q₁₀.**

**What the cancellation does NOT license:**

| not licensed | why |
|---|---|
| **absolute yields in ppm or mM** | undefined ppm basis (§2.3.2) + no recovery |
| **the Arrhenius pre-exponential `A`** | `A` absorbs `f` entirely; every `lnA` in §5 is `ln(f·A)` and is meaningless in absolute terms |
| **cross-paper magnitude comparison** | `f` is specific to this AED, this IS, this extraction. **Never compare a Chan ppm or a Chan `k` to a Kang µg/L, a Hofmann µg/kg, or any other paper's magnitude.** |
| **cross-compound comparison of `k` within this paper** | DMDS carries 2 S atoms, methional and 2-AT carry 1. On an equimolar-sulfur detector a mole of DMDS gives ~2× the signal of a mole of methional. The paper never corrects for this. **DMDS `k` and methional `k` are on different effective scales even inside this chapter.** |
| **assuming `f` is constant across pH** | the cancellation argument holds along the *temperature* axis at fixed pH. Different pH means a different matrix into which DCM is extracted. **Ea is safe per-pH; pH-to-pH `k` ratios are not.** |

---

## §3. EVERY TABLE AND FIGURE, TRANSCRIBED

### §3.1 TABLE I — *"Activation energy for selected flavor compounds"* (p.136) `[F]`

**Transcribed cell by cell. All 18 cells raster-verified at 600 dpi against the OCR layer; zero
disagreements (defect D1).** Units as printed: **kcal/mole**. `r²` printed as `r2` with a
superscript 2.

**1ST REPLICATE**

| pH | Compound | r² | Ea (kcal/mole) |
|---|---|---|---|
| **pH 6** | Dimethyl disulfide | 0.935 | 20.07 |
| pH 6 | Methional | 0.851 | 27.07 |
| pH 6 | 2-Acetylthiophene | *Insufficient data* | *(no value printed)* |
| **pH 7** | Dimethyl disulfide | 0.877 | 15.02 |
| pH 7 | Methional | 0.912 | 24.95 |
| pH 7 | 2-Acetylthiophene | *Insufficient data* | *(no value printed)* |
| **pH 8** | Dimethyl disulfide | 0.857 | 22.80 |
| pH 8 | Methional | 0.797 | 19.39 |
| pH 8 | 2-Acetylthiophene | 0.976 | 32.71 |

**2ND REPLICATE**

| pH | Compound | r² | Ea (kcal/mole) |
|---|---|---|---|
| **pH 6** | Dimethyl disulfide | 0.900 | 19.56 |
| pH 6 | Methional | 0.961 | 21.74 |
| pH 6 | 2-Acetylthiophene | *Insufficient data* | *(no value printed)* |
| **pH 7** | Dimethyl disulfide | 0.805 | 16.20 |
| pH 7 | Methional | 0.956 | 26.02 |
| pH 7 | 2-Acetylthiophene | *Insufficient data* | *(no value printed)* |
| **pH 8** | Dimethyl disulfide | 0.974 | 25.63 |
| pH 8 | Methional | 0.848 | 19.50 |
| pH 8 | 2-Acetylthiophene | 0.972 | 30.81 |

**Footnotes to Table I: [NEG] there are none.** No SD, no SE, no confidence interval, no n, no units
for `k`, no reaction-order column, no pre-exponential column. The table is three columns wide and
nothing else is on the page.

**Derived from Table I `[D]` — the two-replicate spread, which the authors never compute:**

| compound | pH | rep 1 | rep 2 | mean | |Δ| | Δ as % of mean |
|---|---|---|---|---|---|---|
| Dimethyl disulfide | 6 | 20.07 | 19.56 | **19.82** | 0.51 | 2.6 % |
| Dimethyl disulfide | 7 | 15.02 | 16.20 | **15.61** | 1.18 | 7.6 % |
| Dimethyl disulfide | 8 | 22.80 | 25.63 | **24.22** | 2.83 | 11.7 % |
| Methional | 6 | 27.07 | 21.74 | **24.41** | **5.33** | **21.8 %** |
| Methional | 7 | 24.95 | 26.02 | **25.49** | 1.07 | 4.2 % |
| Methional | 8 | 19.39 | 19.50 | **19.45** | 0.11 | 0.6 % |
| 2-Acetylthiophene | 8 | 32.71 | 30.81 | **31.76** | 1.90 | 6.0 % |

> ⚠ **The between-replicate spread reaches 21.8 % (methional, pH 6: 27.07 vs 21.74 kcal/mol, a
> 5.3 kcal/mol = 22 kJ/mol disagreement between two runs of the same condition).** Any single Table I
> value carries at least that much run-to-run uncertainty. **This is the correct error bar to attach
> to a Chan Ea — not the r².**

> ⚠ The authors invoke *"the 95% confidence limits of the calculated Ea's"* (p.135) to argue there is
> "essentially no difference in Ea between pH 6, 7, and 8". **Those confidence limits are never
> printed, in the table, in Figure 6, or anywhere in the chapter.** `[NEG]` The claim is unverifiable
> and with n = 2 a 95 % interval would in any case be enormous.

### §3.2 FIGURE 1 (p.130) — DMDS, 95 °C, ppm vs time `[M]` + `[D]` digitisation

> *"Figure 1. The influence of pH and heating time on the formation of dimethyl disulfide at 95°C."*

**Digitised.** Method: `pdftoppm -r 600 -png` of PDF page 4; plot frame and axis ticks located by
run-length analysis of the binarised (threshold 160/255) raster; the seven y-axis tick marks were
found at pixel rows **1624 (30 ppm), 1901.5 (25), 2177.5 (20), 2453 (15), 2721.5 (10), 2995.5 (5),
3268 (0)**, giving a calibration of **54.8 px per ppm** (residual scatter of the tick spacing < 3 px,
i.e. < 0.06 ppm). Bar tops read as the topmost dark pixel of each column; five x-axis group ticks at
1115, 1406.5, 1695.5, 1984.5, 2275 px assign bars to time groups; within a group the three bars are
pH 6 / pH 7 / pH 8 left-to-right. **Honest error bar: ± 0.15 ppm** (line width, ~4 px).

| time (min, **as labelled on the figure**) | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| 20 | 1.61 | 3.72 | 11.90 |
| 40 | 3.27 | 8.14 | 22.54 |
| 60 | 3.91 | 13.43 | 25.13 |
| 80 | 7.39 | 17.37 | 26.10 |
| 100 | 11.84 | 19.76 | 24.75 |

⚠ **The x-axis reads 20/40/60/80/100 min; the Methods say the 95 °C series was 20/40/60/100/120 min
(defect C4).** The table above uses the figure's own labels. If the Methods are right and the figure
axis is mislabelled, the last two points are at 100 and 120 min, which would flatten every slope by
~20 % — **and would change nothing about the shapes.**

**Authors' reading, verbatim:**
> *"The amount of dimethyl disulfide formed increased as either heating time or pH increased
> (Figure 1). At pH 8, dimethyl disulfide concentration leveled off at longer heating times. Since
> analysis is in progress to determine the amount of reactants remaining during heating, we can not
> comment at this time as to whether this leveling off is due to the exhaustion of reactants or the
> further reaction of dimethyl disulfide to other end products."* (pp.129–131)

**⇒ The digitisation confirms the plateau quantitatively: pH 8 rises 11.90 → 22.54 → 25.13 → 26.10
then *falls* to 24.75. The final point is a decrease, not a plateau.** `[D]`

### §3.3 FIGURE 2 (p.132) — methional, 95 °C, ppm vs time `[M]` + `[D]` digitisation

> *"Figure 2. The influence of heating time and pH on the formation of methional at 95°C."*

Same method, PDF page 6. Six y-ticks at **1261.5 (5 ppm), 1740.5 (4), 2208 (3), 2675 (2),
3152.5 (1), 3626.5 (0)** → **473 px per ppm**. x-group ticks at 895, 1220, 1544, 1868.5, 2193 px.
**Honest error bar: ± 0.02 ppm.**

| time (min, as labelled) | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| 20 | 1.17 | 0.82 | 1.33 |
| 40 | 1.27 | 1.74 | 2.31 |
| 60 | 3.09 | 2.81 | 2.90 |
| 80 | 4.26 | 3.73 | 3.08 |
| 100 | 4.39 | 4.23 | 3.06 |

**Authors' reading, verbatim:**
> *"The concentration of methional in the model system (Figure 2) also increased with heating time
> but initially increased and then decreased with increasing pH. The higher pH's likely favored
> either degradation (or further reaction) of the methional or a shift in pathways to channel
> reactants into alternate end products. We have no data to suggest one alternative over the other."*
> (p.131)

**⇒ Digitisation confirms and sharpens: at 20 min the pH ordering is 8 > 6 > 7; by 100 min it has
inverted completely to 6 > 7 > 8. The pH-8 series plateaus at 3.0–3.1 ppm from 60 min onward
(2.90 → 3.08 → 3.06), while pH 6 is still climbing.** `[D]`

### §3.4 FIGURE 3 (p.133) — 2-acetylthiophene, 95 °C, ppm vs time `[M]` + `[D]` digitisation

> *"Figure 3. The influence of heating time and pH on the formation of 2-acetylthiophene at 95°C."*

Same method, PDF page 7. Six y-ticks at **1328.5 (25 ppm), 1739.5 (20), 2158.5 (15), 2571.5 (10),
2983.5 (5), 3401 (0)** → **82.9 px per ppm**. x-group ticks at 1085.5, 1386.5, 1687.5, 1987, 2287.5
px; slot width ≈ 92 px so bars were assigned to (time, pH) slots by centre position.
**Honest error bar: ± 0.2 ppm.**

| time (min, as labelled) | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| 20 | **not detected** | **not detected** | **not detected** |
| 40 | not detected | not detected | 2.29 |
| 60 | not detected | not detected | 7.37 |
| 80 | 4.39 | not detected | 16.96 |
| 100 | 23.64 | not detected | 18.81 |

**Authors' reading, verbatim:**
> *"The formation of 2-acetylthiophene (Figure 3) was very interesting as the compound was only
> formed at pH 6 and pH 8 and none was detected at pH 7. This again could be due to degradation (or
> further reaction) being favored at this pH."*  (p.131)

**⇒ Digitisation confirms exactly: six bars, four at pH 8 (40–100 min) and two at pH 6 (80 and
100 min), zero at pH 7 at any time.** `[D]`
**⇒ And it explains Table I's `Insufficient data` for 2-AT at pH 6 and pH 7: at pH 6 the compound
appears at only 2 of 5 time points at 95 °C, and at pH 7 never at all — so at most temperatures
there is no slope to fit.** `[D]`

⚠ **The pH-6 series is grotesquely non-linear (0, 0, 0, 4.39, 23.64) — a 5.4× jump in the last
20 min. A "zero-order rate constant" fitted to that is a fiction.** The authors' own zero-order
conclusion is not defensible for 2-acetylthiophene at pH 6, and they implicitly concede it by
refusing to report an Ea there.

### §3.5 FIGURES 4 and 5 (p.134) — the rate-constant ladders `[F]` + `[D]` digitisation

> *"Figure 4. The influence of temperature and pH on the rate constant of dimethyl disulfide."*
> *"Figure 5. The influence of temperature and pH on the rate constant of methional."*

**★ These two figures are the only place in the chapter where the rate constants themselves appear.
Table I gives Ea and r² only. Everything in §5 and §6 therefore depends on this digitisation.**

**Digitisation method (stated in full because the result is load-bearing).** Both are
*pseudo-3D extruded* bar charts, which is the hard case: each bar is drawn as a box with a
parallelogram top face offset up-and-right, so the topmost dark pixel of a bar is **not** its value.
Procedure:

1. PDF page 8 rendered at **600 dpi** (3600 × 5400 px), binarised at threshold 160/255.
2. **y-calibration.** Fig 4: seven tick marks located at rows **1443.5 (k = 3), 1663.5 (2.5),
   1882 (2), 2099.5 (1.5), 2321 (1), 2533 (0.5), 2755 (0)** → **437.2 px per unit k**, baseline row
   **2755.5** (independently confirmed as the plotted zero line, found as a 1348-px-wide horizontal
   run at row 2758). Fig 5: nine ticks resolved, spacing **130.3 px per 0.02** → **6515 px per unit
   k**, baseline row **4734.5** (confirmed as a 1330-px horizontal run at row 4735).
3. **Extrusion geometry.** The depth offset was measured two independent ways and agreed:
   the plot-floor front edge (row 2755.5) vs its back edge (row 2723.5) gives **dy = 32**; the box
   top-left corner (row 3394.5 in Fig 5) vs its `0.2` tick (row 3427) gives **dy = 32.5**. Horizontal
   depth **dx ≈ 38 px**. Bar geometry: group pitch 308 px (Fig 4) / 329.6 px (Fig 5), bar width
   **89 px** (Fig 4) / **95.2 px** (Fig 5), first bar offset 37 px into the group.
4. **Reading rule.** All three bars in a group share the baseline (verified: the fill of the pH 7 and
   pH 8 bars runs continuously down to the zero line). Sampling the topmost dark pixel at
   `x = x_left + 46 … x_left + 92` — a window that is past the *previous* bar's top face (which ends
   at `x_left + dx`) and before the *next* bar's front face — returns the bar's own top-face back
   edge `Y_b`; the front-face top is then `Y_b + dy` and the value is `(baseline − dy − Y_b)/scale`.
5. **Validation** — see §4.3. **The whole-ladder Ea recomputed from these digitised `k` reproduce all
   six of the authors' own printed Table I Ea to within 4–12 %, with `r²` agreeing to the third
   decimal in three of six cases.** That is the check that this digitisation is sound.

**⚠ THE `k` AXIS OF BOTH FIGURES IS LABELLED SIMPLY `k`, WITH NO UNITS.** `[NEG]` The units are
recovered by inference in §3.6.

#### Figure 4 — dimethyl disulfide, digitised `k` `[D]`

| T (°C) | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| 75 | 0.0755 | 0.0640 | 0.0640 |
| 85 | 0.0984 | 0.1349 | 0.1075 |
| 95 | 0.1212 | 0.5032 | 0.3705 |
| 105 | 0.5261 | 0.4552 | 0.6221 |
| 115 | **1.1093** | 0.4803 | **2.6990** |

Honest error bar: **± 0.02 absolute** on every cell (≈ ± 9 px), which is **± 25–30 % relative at
75–85 °C** and **< ± 2 % at 115 °C**. A ± 5.5 px systematic uncertainty in `dy` shifts every value by
a constant **∓ 0.013** — the sensitivity of the derived Ea to that is tabulated in §5.3.

#### Figure 5 — methional, digitised `k` `[D]`

| T (°C) | pH 6 | pH 7 | pH 8 |
|---|---|---|---|
| 75 | 0.0101 | 0.0046 | 0.0052 |
| 85 | 0.0292 | 0.0193 | 0.0315 |
| 95 | 0.0553 | 0.0617 | 0.0358 |
| 105 | **0.1825** | 0.0953 | 0.0625 |
| 115 | **0.1966** | **0.1949** | 0.0970 |

Honest error bar: **± 0.0015 absolute**, i.e. ± 15–30 % at 75 °C and < ± 1 % at 115 °C.

**Authors' reading of Figs 4–5, verbatim:**
> *"The pH had little effect at lower temperatures (up to 115°C for dimethyl disulfide and 105°C for
> methional). While pH did not follow a trend for dimethyl disulfide at 115°C, higher pH's lowered
> the reaction rate of methional at higher temperatures."*  (p.134)

**⇒ Digitisation confirms both halves.** At 115 °C DMDS goes 1.109 (pH 6) → 0.480 (pH 7) → 2.699
(pH 8): non-monotonic, "no trend", exactly as stated. Methional at 105 °C goes 0.183 → 0.095 → 0.063
and at 115 °C 0.197 → 0.195 → 0.097: monotonically **decreasing** with pH, exactly as stated. `[D]`

**⇒ `[NEG]` 2-acetylthiophene has NO rate-constant figure. There is no Figure showing `k` for
2-acetylthiophene at any pH or temperature. Its ladder is unrecoverable from this chapter.**

### §3.6 ★ RECOVERING THE UNITS OF `k` `[D]`

The `k` axes are unitless (§3.5). But the authors state
> *"The rate constants, k, were obtained by calculating the slope of concentration (linear portion of
> plot) versus time."* (p.131)

and the concentration figures are in **ppm** against **min**. Therefore `k` must be **ppm · min⁻¹**.
Cross-check at 95 °C, comparing the Fig-1/Fig-2 slopes with the Fig-4/Fig-5 `k`:

| compound | pH | slope of Fig 1/2 over the linear portion (ppm/min) | digitised `k` from Fig 4/5 | agreement |
|---|---|---|---|---|
| DMDS | 6 | 0.12 | 0.121 | ✔ **exact** |
| DMDS | 7 | 0.21 | 0.503 | ✗ 2.4× |
| DMDS | 8 | 0.33 (0.53 over 20→40 min) | 0.371 | ✔ |
| methional | 6 | 0.052 (20→80 min) | 0.055 | ✔ |
| methional | 7 | 0.042 (40→100 min) | 0.062 | ~ 1.5× |
| methional | 8 | 0.040 (20→60 min) | 0.036 | ✔ |

**⇒ `k` is in ppm · min⁻¹ `[D]` — four of six cells agree within 10 %, which is decisive for the
units.** The two that do not are the pH-7 cells, which is precisely defect C5, and which §4.3 shows
is a *replicate* mismatch: the concentration figures and the rate-constant figures are not the same
run, and the paper never says so.

⚠ **The units are DERIVED, not printed. And because "ppm" itself has an undefined basis (§2.3.2),
`k` in ppm/min is a within-study scale, not a transportable rate.**

### §3.7 FIGURE 6 (p.135) — Ea vs pH `[F]` + `[D]` digitisation

> *"Figure 6. The influence of pH on average activation energy."*

Simple 2D bar chart, two series (Dimethyl disulfide, Methional). Digitised from PDF page 9 at 600 dpi;
seven y-ticks at rows **3087 (30 kcal/mol), 3353.5 (25), 3620.5 (20), 3884.5 (15), 4145 (10),
4409 (5), 4677.5 (0)** → **53.0 px per kcal/mol**. Error bar **± 0.15 kcal/mol**.

| pH | Dimethyl disulfide (digitised) | Methional (digitised) |
|---|---|---|
| 6 | 19.96 | 24.76 |
| 7 | 15.08 | 25.49 |
| 8 | 24.60 | 20.12 |

**⇒ `[D]` Figure 6 is the arithmetic MEAN of the two Table I replicates, confirmed:** Table I means
are 19.82 / 15.61 / 24.22 (DMDS) and 24.41 / 25.49 / 19.45 (methional); the digitised bars match to
within **± 0.7 kcal/mol** in all six cases, and to **0.005 kcal/mol** for methional at pH 7. The
figure caption's word *"average"* is therefore "average over the two replicates". **This was not
stated in the caption or the text.**

**⇒ `[NEG]` 2-acetylthiophene is NOT plotted in Figure 6 — only two of the three compounds appear.**
**⇒ `[NEG]` Figure 6 carries NO error bars, despite the facing text arguing from "95% confidence
limits".**

**Authors' reading, verbatim:**
> *"It was interesting to see that when Ea was plotted against pH (Figure 6), the expected pattern of
> a higher pH having a lower Ea wasn't observed. This may be due to the fact that the variation in pH
> from 6 to 8 was insufficient to discern a trend. When the 95% confidence limits of the calculated
> Ea's are taken into account, there is essentially no difference in Ea between pH 6, 7, and 8."*
> (p.135)

---

## §4. THE PUBLISHED KINETICS `[F]`, AND THE AUDIT OF THEM

### §4.1 What the authors fitted, verbatim

> *"The Macintosh program "Water Analyzer Series - Reaction Kinetics Program V. 2.09" (13) was used to
> do kinetic analysis of the data. The reaction order for the compounds selected (dimethyl disulfide,
> methional, and 2-acetylthiophene) were determined by plotting concentration versus time and then
> using linear regression to determine how well the data fit the straight line. The r² for all three
> compounds at zero order ranged form 0.83 to 0.99 and were consistently better than any other order.
> It was, therefore, concluded that the formation of dimethyl disulfide, methional, and
> 2-acetylthiophene all followed zero order reaction kinetics."*  (p.131, raster ✔, typo `form` is the
> authors')

> *"The rate constants, k, were obtained by calculating the slope of concentration (linear portion of
> plot) versus time. … Therefore, the Arrhenius equation was used to determine the activation energy
> (Ea). The value of Ea was obtained by doing linear regression on a plot of ln k vs. 1/T."*
> (pp.131–134)

| published quantity | value | provenance |
|---|---|---|
| **reaction order** | ★ **zero order**, for all three compounds, at all pH and all T | `[F]` p.131 |
| order-selection criterion | linear regression `r²`, zero order **"consistently better than any other order"**; no AIC, no F-test, no residual analysis | `[F]` p.131 |
| `r²` of the zero-order fits | **0.83 – 0.99** (p.131) / **0.81 – 0.98** (p.135) — ⚠ two incompatible statements, defect C3 | `[F]` |
| **rate constants `k`** | **only in Figures 4 and 5, unitless axis, DMDS and methional only.** No table of `k` exists. | `[F]` |
| **Ea** | Table I, 14 values (7 conditions × 2 replicates), **kcal/mole** | `[F]` p.136 |
| Arrhenius method | ordinary least squares on `ln k` vs `1/T`, 5 temperatures | `[F]` p.134 |
| fitting software | **Water Analyzer Series – Reaction Kinetics Program V.2.09**, Labuza, Nelson & Nelson, University of Minnesota, 1991 (ref 13) | `[F]` |
| **pre-exponential `A`** | ★ **[NEG] never reported, for any compound, at any pH** | — |
| **Q₁₀** | ★ **[NEG] never reported** | — |
| **half-life** | ★ **[NEG] never reported** (and is undefined for a zero-order formation) | — |
| **rate constants for 2-acetylthiophene** | ★ **[NEG] never reported at all** | — |
| 95 % confidence limits | ★ **[NEG] invoked in the text (p.135), never printed** | — |

### §4.2 ★ RESOLVING WHAT TABLE I's `r²` COLUMN MEANS `[D]`

Nothing in the chapter says whether Table I's `r²` is the zero-order concentration-vs-time fit or the
Arrhenius `ln k` vs `1/T` fit (defect C3). **The audit resolves it.** Refitting the Arrhenius line
from the digitised `k` ladders gives `r²` values of **0.884, 0.806, 0.965** (DMDS pH 6/7/8) and
**0.962, 0.966, 0.873** (methional pH 6/7/8). Table I's 2nd-replicate `r²` are **0.900, 0.805,
0.974 / 0.961, 0.956, 0.848**. Three of the six agree to within 0.001–0.010 — in particular
**DMDS pH 7: mine 0.806 vs Table I 0.805, and methional pH 6: mine 0.962 vs Table I 0.961.**

> **⇒ Table I's `r²` column is the ARRHENIUS regression `r²` (`ln k` vs `1/T`, 5 points), not the
> zero-order `r²`. The prose ranges "0.83–0.99" and "0.81–0.98" refer to a *different* set of fits
> that is never tabulated.** Record this so a later wave does not misread Table I's `r²` as evidence
> about the zero-order assumption. `[D]`

### §4.3 ★ THE REFIT AUDIT — DOES TABLE I REPRODUCE FROM THE AUTHORS' OWN `k`? `[D]`

Ordinary least squares on `ln k` vs `1/T`, 5 temperatures (348.15 – 388.15 K), `R = 8.314462
J mol⁻¹ K⁻¹`, on the digitised ladders of §3.5. **My refit is `[D]` throughout; the Table I columns
are `[F]`.**

| compound | pH | **[D] my refit Ea (kcal/mol)** | **[D] my r²** | [F] Table I rep 1 (Ea / r²) | [F] Table I rep 2 (Ea / r²) | verdict |
|---|---|---|---|---|---|---|
| Dimethyl disulfide | 6 | **18.76** | **0.884** | 20.07 / 0.935 | 19.56 / 0.900 | ✔ reproduces rep 2 to −4.1 % |
| Dimethyl disulfide | 7 | **14.28** | **0.806** | 15.02 / 0.877 | 16.20 / **0.805** | ✔ reproduces rep 2 to −11.9 %; **r² matches to 0.001** |
| Dimethyl disulfide | 8 | **24.71** | **0.965** | 22.80 / 0.857 | 25.63 / 0.974 | ✔ reproduces rep 2 to −3.6 % |
| Methional | 6 | **20.94** | **0.962** | 27.07 / 0.851 | 21.74 / **0.961** | ✔ reproduces rep 2 to −3.7 %; **r² matches to 0.001** |
| Methional | 7 | **24.54** | **0.966** | 24.95 / 0.912 | 26.02 / 0.956 | ✔ reproduces rep 2 to −5.7 % |
| Methional | 8 | **17.70** | **0.873** | 19.39 / 0.797 | 19.50 / 0.848 | ✔ reproduces rep 2 to −9.2 % |
| 2-Acetylthiophene | 8 | ★ **cannot be audited** | — | 32.71 / 0.976 | 30.81 / 0.972 | ✗ **no `k` figure exists (§3.5) — the two highest Ea in the paper are unauditable** |

**⇒ THREE FINDINGS FROM THE AUDIT:**

1. ★ **Figures 4 and 5 plot the 2ND REPLICATE.** Every one of the six refits lands closer to the
   2nd-replicate Ea than to the 1st, and two of the six reproduce the 2nd-replicate `r²` to the third
   decimal place. **The chapter never states which replicate is plotted.** `[D]` — and this also
   dissolves defect C5: the concentration figures (Figs 1–3) and the rate-constant figures (Figs 4–5)
   are simply different runs.
2. **Nothing fails to reproduce.** All six auditable Ea come out within 12 % (mean bias **−6.3 %**) of
   the authors' printed values. The systematic −6.3 % is fully explained by a ~5.5 px line-width
   offset in the extrusion depth `dy`: setting `dy = 38` px instead of the geometric 32.5 px removes
   the bias completely (§5.3), giving 19.98 / 15.52 / 26.15 (DMDS) and 21.46 / 25.73 / 18.66
   (methional) against Table I rep 2's 19.56 / 16.20 / 25.63 / 21.74 / 26.02 / 19.50 — **all six
   within 5 %.** **This is a calibration, and I flag it as such: it is not an independent
   confirmation.** The independent confirmation is item (1) above, which uses `r²` — a quantity that
   is invariant to any constant offset in `dy`.
3. ★ **The 2-acetylthiophene Ea (32.71 and 30.81 kcal/mol at pH 8) cannot be checked at all** — no
   rate constants are published for it. They are also the two largest Ea in the chapter, and they are
   derived from a series that at 95 °C is detectable at only 4 of 5 time points and at pH 6 is
   0/0/0/4.39/23.64 ppm. **Refuse them (§9).**

### §4.4 Audit of the zero-order claim itself `[D]`

The zero-order conclusion (`r² = 0.83–0.99`) is not supportable for two of the three compounds, on
the authors' own figures:

| compound | pH | the 95 °C series (ppm) | zero-order? |
|---|---|---|---|
| DMDS | 6 | 1.61, 3.27, 3.91, 7.39, 11.84 | **accelerating** — a straight line through these is a poor description; the last interval is 2.2× the first |
| DMDS | 8 | 11.90, 22.54, 25.13, 26.10, **24.75** | **saturating and then declining** — the authors say so themselves |
| methional | 8 | 1.33, 2.31, 2.90, 3.08, 3.06 | **plateaued from 60 min** |
| 2-AT | 6 | 0, 0, 0, 4.39, 23.64 | ★ **not a line by any reading** |

**⇒ The authors' own escape hatch is the parenthesis "(linear portion of plot)" (p.131): `k` is the
slope of a *selected sub-range*, and which sub-range was selected is never stated for any of the 30
condition-cells. That is an undeclared analyst degree of freedom sitting underneath every `k` in
Figures 4–5 and therefore underneath every Ea in Table I.** `[D]`

To their credit the authors state the underlying problem in full:
> *"the data represent the concentration of each volatile in the system at any given time. It is well
> documented in the literature that few flavor compounds are end products of the Maillard reaction. …
> Thus the kinetics observed are the dynamic result of formation versus further reaction. Since we
> know so little about the reactions occurring in the Maillard reaction, we can do little other than
> accept an "observed" reaction rate and the resultant kinetics."*  (p.131)

and
> *"the kinetic data obtained in this study were from a complex model system i.e. there was more than
> one amino acid present in the model system. … Thus the kinetics observed in any study are very
> system dependent and cannot be readily compared."*  (p.131)

**⇒ The authors explicitly disclaim cross-study comparability of their own rate constants. Honour
that: this paper's `k` are not transportable; its Ea are (§2.4).**

---

## §5. DERIVED APPARENT ACTIVATION ENERGIES `[D]`

Every number in this section is **derived by me and never printed by the paper**, except where a
Table I value is quoted for comparison. Six ladders qualify (≥ 3 temperatures): 2 compounds × 3 pH,
each with **5 temperatures**. 2-acetylthiophene has **zero** qualifying ladders (§3.5, §4.3).

Basis: `ln k = ln A − Ea/(R T)`, OLS on 5 points, `R = 8.314462 J mol⁻¹ K⁻¹`, `T` in K.
Primary numbers use the **geometric** digitisation (`dy = 32.5 px`); §5.3 gives the sensitivity.

### §5.1 Whole-ladder fits, 75 → 115 °C `[D]`

| compound | pH | **Ea (kJ/mol)** | **Ea (kcal/mol)** | **r²** | `ln(f·A)` | [F] Table I mean for comparison |
|---|---|---|---|---|---|---|
| Dimethyl disulfide | 6 | **78.5** | 18.76 | 0.884 | 24.17 | 19.82 |
| Dimethyl disulfide | 7 | **59.8** | 14.28 | 0.806 | 18.16 | 15.61 |
| Dimethyl disulfide | 8 | **103.4** | 24.71 | 0.965 | 32.73 | 24.22 |
| Methional | 6 | **87.6** | 20.94 | 0.962 | 25.80 | 24.41 |
| Methional | 7 | **102.7** | 24.54 | 0.966 | 30.38 | 25.49 |
| Methional | 8 | **74.0** | 17.70 | 0.873 | 20.80 | 19.45 |

⚠ **`ln(f·A)` is `ln A` contaminated by the unknown semi-quant response factor `f` (§2.4). It is
recorded only so the fits are reproducible. It is NOT a pre-exponential and must never be used as
one.**

### §5.2 ★ LEG-BY-LEG apparent Ea — the result that matters `[D]`

Each cell is the two-point Arrhenius slope between adjacent temperatures:
`Ea = −R · ln(k₂/k₁) / (1/T₂ − 1/T₁)`. The ratio `k₂/k₁` is given alongside because over a 10 K step
it *is* the Q₁₀ (§5.4).

#### Dimethyl disulfide

| pH | 75→85 °C | 85→95 °C | 95→105 °C | 105→115 °C | whole ladder |
|---|---|---|---|---|---|
| **6** | **27.4 kJ** (6.6 kcal) ×1.30 | **22.9 kJ** (5.5) ×1.23 | ★ **169.9 kJ** (40.6) **×4.34** | **91.0 kJ** (21.8) ×2.11 | 78.5 kJ |
| **7** | 77.3 kJ (18.5) ×2.11 | **144.3 kJ** (34.5) ×3.73 | ★ **−11.6 kJ** (−2.8) **×0.90** | 6.6 kJ (1.6) ×1.06 | 59.8 kJ |
| **8** | 53.7 kJ (12.8) ×1.68 | **135.7 kJ** (32.4) ×3.45 | 60.0 kJ (14.3) ×1.68 | ★ **179.1 kJ** (42.8) **×4.34** | 103.4 kJ |

#### Methional

| pH | 75→85 °C | 85→95 °C | 95→105 °C | 105→115 °C | whole ladder |
|---|---|---|---|---|---|
| **6** | 109.6 kJ (26.2) ×2.88 | 70.1 kJ (16.7) ×1.89 | **138.3 kJ** (33.1) ×3.30 | ★ **9.1 kJ** (2.2) **×1.08** | 87.6 kJ |
| **7** | **148.8 kJ** (35.6) ×4.20 | 127.2 kJ (30.4) ×3.19 | 50.3 kJ (12.0) ×1.54 | 87.3 kJ (20.9) ×2.05 | 102.7 kJ |
| **8** | ★ **186.3 kJ** (44.5) **×6.03** | ★ **14.0 kJ** (3.4) **×1.14** | 64.6 kJ (15.4) ×1.75 | 53.7 kJ (12.8) ×1.56 | 74.0 kJ |

> ### ★ **NOT ONE OF THE SIX LADDERS IS ARRHENIUS-LINEAR.** The per-leg Ea within a single
> compound × pH ladder span factors of **6.2× (DMDS pH 6: 22.9 → 169.9 kJ/mol)**, **13.3×
> (methional pH 8: 14.0 → 186.3)**, **15.2× (methional pH 6: 9.1 → 138.3)**, and DMDS pH 7 even
> **changes sign** (−11.6 kJ/mol on the 95→105 leg). The whole-ladder `r²` of 0.81–0.97 conceals all
> of this: five points fitted by OLS will return a respectable `r²` for a curve with one large break
> in it.

### §5.3 Sensitivity of the derived Ea to the digitisation `[D]`

| compound × pH | Ea, `dy = 26 px` | **Ea, `dy = 32.5 px` (primary, geometric)** | Ea, `dy = 38 px` (line-width corrected) | Ea, `dy = 44 px` | [F] Table I rep 2 |
|---|---|---|---|---|---|
| DMDS pH 6 | 17.57 | **18.76** | 19.98 | 21.64 | 19.56 |
| DMDS pH 7 | 13.11 | **14.28** | 15.52 | 17.30 | 16.20 |
| DMDS pH 8 | 23.32 | **24.71** | 26.15 | 28.17 | 25.63 |
| methional pH 6 | 20.38 | **20.94** | 21.46 | 22.08 | 21.74 |
| methional pH 7 | 23.39 | **24.54** | 25.73 | 27.36 | 26.02 |
| methional pH 8 | 16.75 | **17.70** | 18.66 | 19.94 | 19.50 |

(kcal/mol.) **The whole ± 9 px digitisation envelope moves any Ea by at most ± 12 %, which is
smaller than this paper's own between-replicate spread of up to 21.8 % (§3.1).**
**The per-leg *ratios* `k₂/k₁` in §5.2 move by at most ± 15 % on the low-T legs and < 2 % on the
high-T legs, and no sign or ordering in §5.2 changes anywhere in that envelope** — the qualitative
shape findings are robust to the digitisation.

### §5.4 Q₁₀ `[D]` — never printed by the authors

Because every step in this ladder is exactly 10 K, **the `k₂/k₁` ratios in §5.2 are the Q₁₀
directly.** Summary of the extremes:

| | lowest Q₁₀ | highest Q₁₀ |
|---|---|---|
| Dimethyl disulfide | **0.90** (pH 7, 95→105 °C) — *rate falls with temperature* | **4.34** (pH 6, 95→105 and pH 8, 105→115) |
| Methional | **1.08** (pH 6, 105→115 °C) | **6.03** (pH 8, 75→85 °C) |

### §5.5 ★ THE HONEST CAVEAT CHAIN — attach all seven to any Ea taken from §5

1. ⚠ **`k` here is a zero-order formation slope, not a rate constant of an elementary step.** The
   authors are explicit that it is *"the dynamic result of formation versus further reaction"*
   (p.131). An apparent Ea from such a slope is a lumped quantity over formation *minus* consumption,
   and it changes sign when consumption overtakes formation — which is exactly what the −11.6 kJ/mol
   leg (DMDS pH 7, 95→105 °C) is.
2. ⚠ **Yield-at-fixed-time is not a rate constant, and `k` here is one step better than that but not
   two**: it is the slope of an analyst-selected *"linear portion"* whose bounds are nowhere stated
   (§4.4). The selection of that window is an undeclared degree of freedom.
3. ⚠ **Precursor depletion depresses apparent Ea at the top of a ladder** — the standard failure
   mode — and here the top-leg collapses (methional pH 6, 9.1 kJ/mol; DMDS pH 7, +6.6 kJ/mol) are
   exactly where that would show. **But see §5.6: on a mass balance, methionine depletion cannot be
   the cause.**
4. ⚠ **Semi-quant scale.** Every `k` carries an unmeasured constant `f` (§2.4). It cancels in Ea. It
   does **not** cancel in `ln A`, in any absolute yield, or in any comparison to another paper.
5. ⚠ **pH drift.** 0.1 M phosphate against 0.75 M amine (§1). The pH was measured every time and
   never reported. **"pH 6/7/8" is an initial condition, and the drift is presumably larger at higher
   T and longer t — i.e. correlated with the very axis Ea is fitted along.**
6. ⚠ **Headspace/N₂ confounder specific to this reactor** (§1): five withdrawals, 62.5 % of the
   charge, and five N₂ repressurisations, changing the oxidant budget of an oxidative coupling.
7. ⚠ **OCR risk: LOW for this paper** (§0.1) — Table I is clean and the digits are corroborated by the
   CAS check. **The larger risk is the paper's own internal contradictions C1–C5, not the scan.**

### §5.6 ★ A MASS BALANCE THE PAPER DOES NOT DO — and it answers the authors' open question `[D]`

The authors leave the plateau unexplained: *"we can not comment at this time as to whether this
leveling off is due to the exhaustion of reactants or the further reaction of dimethyl disulfide"*
(pp.129–131). The loadings settle it.

Methionine is the **only** sulfur source: 0.075 mol / 0.400 L = **187.5 mM**. Taking the peak
digitised yields at 95 °C and reading ppm as mg/L:

| compound | peak (ppm) | MW | mM | S atoms | methionine consumed | **% of the 187.5 mM Met pool** |
|---|---|---|---|---|---|---|
| dimethyl disulfide | 26.1 | 94.20 | 0.277 | 2 | 0.554 mM | **0.30 %** |
| 2-acetylthiophene | 23.6 | 126.18 | 0.187 | 1 | 0.187 mM | **0.10 %** |
| methional | 4.4 | 104.17 | 0.042 | 1 | 0.042 mM | **0.02 %** |
| **total** | | | | | **0.78 mM** | ★ **0.42 %** |

> ### ★ **THE THREE MEASURED PRODUCTS ACCOUNT FOR LESS THAN 0.5 % OF THE METHIONINE, AND GLUCOSE IS AT 1.25 M. NEITHER BULK REACTANT IS ANYWHERE NEAR EXHAUSTED. THE PLATEAUS AND THE TOP-LEG Ea COLLAPSES CANNOT BE BULK PRECURSOR DEPLETION.**

Careful hedge, because the argument does not prove as much as it looks: what *can* be exhausted is
not methionine or glucose but (i) the **transient dicarbonyl pool** that drives the Strecker step,
whose steady-state concentration is small and which browning consumes elsewhere, (ii) the
**methanethiol** intermediate standing between methional and DMDS, and (iii) the **oxidant**, which
the N₂ purge is actively removing (§1). And the pH is drifting. **What the mass balance does rule out
is the simplest reading — that the ladder flattens because the substrates ran out. It did not.**
`[D]`

> ⚠ **If ppm is on the extract basis rather than the aqueous basis (§2.3.2), every percentage above
> falls by a further 50–100×, and the conclusion only strengthens.**

---

## §6. ★ THE KANG SWITCH-ON CROSS-CHECK

### §6.1 The reference finding being tested

Kang et al. 2026 (TTCA/Cys, pH 7, 120 min, µg/L) reports that **2-methyl-3-furanthiol (MFT)** and
**2-furfurylthiol (FFT)** are strongly **non-Arrhenius**: nearly flat from 100 → 120 °C
(**×1.12** MFT, **×1.10** FFT, i.e. apparent Ea only **~6–7 kJ/mol**), then a jump to 140 °C
(**×4.26** MFT, **×2.78** FFT, apparent Ea **97.8** and **69.2 kJ/mol** on the 120→140 leg) — while
the sulfur class *as a whole* decelerates (Ea 57.5 → 35.2 kJ/mol).

### §6.2 ★ THE DIRECT ANSWER: CHAN CANNOT TEST THE 120–140 °C WINDOW

> ### **Chan's ladder tops out at 115 °C — five kelvin BELOW Kang's lower switch-on anchor. There is no Chan datum at or above 120 °C. Chan therefore cannot confirm, refute, or locate Kang's 120→140 °C switch-on, and nothing in this dossier should be read as doing so.**

The two ladders do not even overlap except at their extreme ends: Chan runs 75–115 °C, Kang runs
100–140 °C, and the only shared temperature is Kang's 100 °C, which falls between two Chan rungs.

### §6.3 ★ THE USEFUL COMPLEMENTARY RESULT: THE BASELINE REGIME IS **NOT** LINEAR

The premise the switch-on framing implicitly rests on — that below the jump there is a
well-behaved Arrhenius baseline — **is false in the 75–115 °C aqueous window.** §5.2 shows every one
of the six ladders breaking, several of them by more than an order of magnitude in apparent Ea, and
one changing sign.

**Per-leg apparent Ea, laid out so a later wave can splice Chan onto Kang** (kJ/mol; Kang's legs
appended in the same units for alignment; **Chan `[D]`, Kang `[C]` from the K5 dossier**):

| ladder | 75→85 | 85→95 | 95→105 | 105→115 | *(gap)* | 100→120 | 120→140 |
|---|---|---|---|---|---|---|---|
| **DMDS pH 6** `[D]` | 27.4 | 22.9 | **169.9** | 91.0 | — | — | — |
| **DMDS pH 7** `[D]` | 77.3 | 144.3 | **−11.6** | 6.6 | — | — | — |
| **DMDS pH 8** `[D]` | 53.7 | 135.7 | 60.0 | **179.1** | — | — | — |
| **methional pH 6** `[D]` | 109.6 | 70.1 | 138.3 | **9.1** | — | — | — |
| **methional pH 7** `[D]` | 148.8 | 127.2 | 50.3 | 87.3 | — | — | — |
| **methional pH 8** `[D]` | 186.3 | **14.0** | 64.6 | 53.7 | — | — | — |
| *Kang MFT, pH 7* `[C]` | — | — | — | — | | **~6–7** | **97.8** |
| *Kang FFT, pH 7* `[C]` | — | — | — | — | | **~6–7** | **69.2** |

⚠ **The blank column is a real 5 K gap (115 → 120 °C) plus a total mismatch of system, matrix,
precursor, analyte and quantification basis. This table is an alignment aid for a later wave, not a
joined curve. Do not fit across it.**

### §6.4 ★ WHERE THE BREAKS ARE, AND THE ONE THAT LOOKS LIKE KANG's

**The breaks are not scattered at random. They sort by compound class.**

| ladder | shape | break location |
|---|---|---|
| **DMDS pH 6** | ★ **flat, flat, JUMP, elevated** | **switch-on at 95 → 105 °C** |
| **DMDS pH 8** | rise, jump, dip, **JUMP** | **switch-on at 105 → 115 °C** (and an earlier one at 85→95) |
| DMDS pH 7 | rise, jump, then **stall/reverse** | switch-*off* at 95 °C |
| **methional pH 6** | rise, rise, rise, **COLLAPSE** | **switch-off at 105 → 115 °C** |
| **methional pH 8** | steep, **COLLAPSE**, mild recovery | **switch-off at 85 → 95 °C** |
| methional pH 7 | steep, steep, dip, recover | mildest of the six; the only near-monotone ladder |

> ### ★ **DMDS AT pH 6 REPRODUCES KANG's MFT SHAPE ALMOST EXACTLY, DISPLACED ~25 K DOWNWARD.**
>
> | | Kang MFT (pH 7) `[C]` | **Chan DMDS pH 6 `[D]`** |
> |---|---|---|
> | low legs | 100→120 °C: **×1.12** (Ea ~6–7 kJ/mol) | 75→85: **×1.30**, 85→95: **×1.23** (Ea 27.4, 22.9 kJ/mol) |
> | the jump | 120→140 °C: **×4.26** (Ea **97.8 kJ/mol**) | 95→105 °C: **×4.34** (Ea **169.9 kJ/mol**) |
> | after the jump | *(no higher rung)* | 105→115 °C: ×2.11 (Ea 91.0 kJ/mol) |
>
> The jump *magnitude* is nearly identical (×4.34 vs ×4.26). The jump *onset* is at 95–105 °C in
> Chan's aqueous pH-6 system versus 120–140 °C in Kang's. The pre-jump baseline is shallower in Kang
> (~6–7 kJ/mol) than in Chan (23–27 kJ/mol). **⚠ This is a shape correspondence between two entirely
> different systems, analytes and quantification bases. It is a hypothesis to be tested by a
> dedicated wave, not a transferable result — and note that Chan pH 6 is the condition, not pH 7.**

### §6.5 ★ THE MECHANISTIC READING, AND WHY IT MATTERS FOR KANG `[D]`

The sorting in §6.4 is not arbitrary. Note what the two compounds *are* in the paper's own framing
(p.127 abstract): **methional is the primary Strecker aldehyde; dimethyl disulfide and
2-acetylthiophene are "two secondary products"**, formed downstream. Methional → methanethiol →
oxidative coupling → DMDS is the canonical chain.

**The observed pattern is: the primary Strecker aldehyde saturates at the top of the ladder while the
secondary disulfide switches on there.** At pH 6, methional's top leg collapses to 9.1 kJ/mol
(×1.08) at exactly the temperature range where DMDS jumps ×4.34. **That is a precursor→product
handoff, not two independent Arrhenius laws** — the methional pool stops growing because it is being
consumed into the secondary channel that is simultaneously accelerating.

**Bearing on Kang:** Kang's MFT and FFT are *also* secondary heterocyclic sulfur products, and Kang's
"sulfur class as a whole decelerates (57.5 → 35.2 kJ/mol)" while MFT/FFT specifically accelerate.
**That is the same split Chan sees: the class-average (dominated by primary/simple sulfur species)
decelerates while the secondary heterocyclic/coupled species switch on.** Chan therefore does
*complementary* work: it establishes that this split is **not a 120 °C phenomenon** — it is already
running at 95–115 °C in water, at least at pH 6.

**⚠ Four reasons to hold this loosely `[D]`:**

1. Chan's `k` for DMDS and methional are on **different effective sulfur scales** (2 S vs 1 S,
   §2.4) — which does not affect either Ea but does mean the handoff is qualitative, not
   stoichiometric.
2. The DMDS switch-on is at **pH 6**, and DMDS at **pH 7** does the opposite (switch-*off* at 95 °C).
   Kang is at pH 7. **On the matching pH, Chan's disulfide decelerates.**
3. Chan's 2-acetylthiophene — the compound structurally closest to Kang's MFT/FFT, being the only
   heterocyclic sulfur analyte in the chapter — **has no rate-constant ladder at all** (§3.5).
   The most informative comparison is the one the paper does not enable.
4. The pH is drifting and unmeasured (§5.5.5), so "pH 6 versus pH 7" may not mean at 115 °C what it
   meant at t = 0.

### §6.6 The one number Chan *does* contribute to the switch-on question `[D]`

2-acetylthiophene's Ea — **32.71 and 30.81 kcal/mol = 136.9 and 128.9 kJ/mol at pH 8** — is the
**highest in the chapter by a wide margin** (next highest is methional pH 7 at 26.02 kcal/mol) and is
the only heterocyclic-sulfur Ea in it. **A high Ea for the heterocyclic sulfur product is exactly
what a switch-on regime predicts.** ⚠ But §4.3 shows those two values are **unauditable** (no `k`
published) and rest on a 95 °C series that is 0/0/0/4.4/23.6 ppm at pH 6 and detected at only 4 of 5
times at pH 8. **They are recorded here as a directional prior only, and REFUSED as parameters
(§9).**

---

## §7. VERIFIED NEGATIVES `[NEG]` — do not re-open this paper for any of these

| # | what a reader might hope for | status |
|---|---|---|
| N1 | **2-methyl-3-furanthiol (MFT)** | ★ **absent. Never mentioned, never detected, never named.** |
| N2 | **2-furfurylthiol (FFT / furfuryl mercaptan)** | ★ **absent. Never mentioned.** |
| N3 | **H₂S** | ★ **never measured, never mentioned.** The AED sulfur channel would not retain it on a DB-5 under these conditions anyway. |
| N4 | **methanethiol** (the direct DMDS precursor) | ★ **never measured, never mentioned** — despite the paper's whole DMDS story requiring it |
| N5 | **any thiol sink / thiol–carbonyl or thiol–HMF crosslink** | absent; no thiol species of any kind is quantified |
| N6 | **cysteine, cystine, thiamine, or any other sulfur precursor** | absent — **methionine is the only S source in the system** |
| N7 | **rate constants for 2-acetylthiophene** | ★ **never published, at any pH or T.** Its two Ea values are unauditable |
| N8 | **units on the `k` axis** of Figures 4 and 5 | ★ **absent.** Units recovered only by inference (§3.6) |
| N9 | **pre-exponential factor `A`** | ★ never reported for any compound |
| N10 | **Q₁₀ or temperature quotient** | never reported |
| N11 | **half-life** | never reported (and undefined for a zero-order formation) |
| N12 | **95 % confidence limits on Ea** | ★ **argued from in the text (p.135) and printed NOWHERE** |
| N13 | **SD / SE / error bars** | ★ absent from Table I, from Figures 1–6, and from the text |
| N14 | **the measured pH of the heated samples** | ★ **measured on every single sample (p.128) and reported nowhere** |
| N15 | **buffer molarity** | ✔ **NOT a negative — it IS given: 0.1 M phosphate** (§1). What is missing is any evidence the buffer held |
| N16 | **calibration curve / response factors / recovery** | ★ absent; deferred to ref (14), "in press" (§2.2) |
| N17 | **the "linear portion" window** used to compute each `k` | ★ never stated for any of the 30 condition-cells |
| N18 | **which replicate Figures 1–5 plot** | ★ never stated. Resolved by audit for Figs 4–5 (= rep 2, §4.3); **still unknown for Figs 1–3** |
| N19 | **the non-sulfur volatiles** (pyrazines etc.) | ★ deliberately excluded — *"will be discussed in another article (12)"*: Chan & Reineccius, *5th Int. Symp. Maillard Reaction*, **"in press"** |
| N20 | **reactant depletion data** | ★ *"analysis is in progress to determine the amount of reactants remaining during heating"* (p.131) — **never reported in this chapter** |
| N21 | **water activity, ionic strength, dissolved O₂** | none reported |
| N22 | **come-up time** to set point | ★ noted as observed (*"Initial temperature was noted"*) but **never reported**; `t = 0` is defined at set-point with no correction |
| N23 | **browning / colour / absorbance** | not measured |
| N24 | **a table of rate constants** | ★ **does not exist.** `k` appears only as bars in two figures |
| N25 | **sensory thresholds** | asserted qualitatively (*"extremely low sensory threshold"*, p.135) with **no number and no citation** |
| N26 | **temperature above 115 °C** | ★ **none. 115 °C is the ceiling. Cannot address the 120–140 °C window (§6.2)** |
| N27 | **single-amino-acid control systems** | ★ absent — and the authors flag this as the study's main limitation themselves (p.131) |

---

## §8. CONSOLIDATED PARAMETER TABLE

`[M]` measured by the authors · `[C]` cited by them · `[F]` fitted by the authors · `[D]` derived by
me, never printed by the paper · `[NEG]` verified negative.

| # | parameter | value | units | condition | provenance | anchor |
|---|---|---|---|---|---|---|
| 1 | glucose loading | 0.5 (= 90 g) | mole | all runs, in 400 mL | `[M]` | p.128 |
| 2 | methionine loading | 0.075 (= 11.19 g) | mole | all runs | `[M]` | p.128 |
| 3 | phenylalanine loading | 0.075 (= 12.39 g) | mole | all runs | `[M]` | p.128 |
| 4 | proline loading | 0.075 (= 8.64 g) | mole | all runs | `[M]` | p.128 |
| 5 | leucine loading | 0.075 (= 9.84 g) | mole | all runs | `[M]` | p.128 |
| 6 | solvent volume | 400 | mL distilled water | all runs | `[M]` | p.128 |
| 7 | glucose concentration | 1.25 | M | all runs | `[D]` | from #1, #6 |
| 8 | each amino acid concentration | 0.1875 | M | all runs | `[D]` | from #2–#6 |
| 9 | total amino acid concentration | 0.75 | M | all runs | `[D]` | from #2–#6 |
| 10 | **methionine (= total S source)** | **187.5** | **mM** | all runs | `[D]` | from #2, #6 |
| 11 | buffer | 0.1 M phosphate + NaOH to pH | — | all runs | `[M]` | p.128 |
| 12 | initial pH levels | 6, 7, 8 | — | 3 levels | `[M]` | p.128 |
| 13 | temperatures | 75, 85, 95, 105, 115 | °C | 5 levels | `[M]` | p.128 |
| 14 | times, 75 °C | 1.5, 3.0, 4.5, 6.0, 7.5 | h | | `[M]` | p.128 |
| 15 | times, 85 °C | 0.5, 1.0, 1.5, 2.0, 2.5 | h | | `[M]` | p.128 |
| 16 | times, 95 °C (Methods) | 20, 40, 60, 100, 120 | min | ⚠ contradicts the figures (C4) | `[M]` | p.128 |
| 17 | times, 95 °C (Figures 1–3 axes) | 20, 40, 60, 80, 100 | min | ⚠ contradicts the Methods (C4) | `[M]` | pp.130,132,133 |
| 18 | times, 105 °C | 10, 20, 30, 40, 50 | min | | `[M]` | p.128 |
| 19 | times, 115 °C | 5, 10, 15, 20, 25 | min | | `[M]` | p.128 |
| 20 | aliquot volume | 50 | mL per sampling | 5 per run | `[M]` | p.128 |
| 21 | total charge withdrawn | 250 mL = **62.5 %** | of the 400 mL charge | per run | `[D]` | from #6, #20 |
| 22 | N₂ repressurisation | ca. 50 | mL per sampling | 5 per run | `[M]` | p.128 |
| 23 | replication n | 2 | independent runs | all conditions | `[M]` | p.136 |
| 24 | internal standard | 4-methylthiazole, 500 ppm in DCM | — | all samples | `[M]` | p.128 |
| 25 | extraction | 3 × 5 mL DCM, → 0.5–1.0 mL under N₂ | — | all samples | `[M]` | p.128 |
| 26 | AED sulfur line | 181 | nm | quantification channel | `[M]` | p.129 |
| 27 | reaction order | **zero** | — | all 3 compounds, all pH, all T | `[F]` | p.131 |
| 28 | zero-order `r²` range (p.131) | 0.83 – 0.99 | — | ⚠ conflicts with #29 | `[F]` | p.131 |
| 29 | zero-order `r²` range (p.135) | 0.81 – 0.98 | — | ⚠ conflicts with #28 | `[F]` | p.135 |
| 30 | Table I `r²` = **Arrhenius** `r²` | 0.797 – 0.976 | — | 14 fits | `[F]`/`[D]` id. | p.136, §4.2 |
| 31 | **Ea DMDS pH 6** | **20.07 / 19.56** (mean 19.82) | kcal/mole | rep 1 / rep 2, 75–115 °C | `[F]` | p.136 |
| 32 | **Ea DMDS pH 7** | **15.02 / 16.20** (mean 15.61) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 33 | **Ea DMDS pH 8** | **22.80 / 25.63** (mean 24.22) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 34 | **Ea methional pH 6** | **27.07 / 21.74** (mean 24.41) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 35 | **Ea methional pH 7** | **24.95 / 26.02** (mean 25.49) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 36 | **Ea methional pH 8** | **19.39 / 19.50** (mean 19.45) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 37 | **Ea 2-acetylthiophene pH 8** | **32.71 / 30.81** (mean 31.76) | kcal/mole | rep 1 / rep 2 | `[F]` | p.136 |
| 38 | Ea 2-AT pH 6 and pH 7 | *"Insufficient data"* | — | both replicates | `[F]`/`[NEG]` | p.136 |
| 39 | Ea in SI units, DMDS pH 6/7/8 | 82.9 / 65.3 / 101.3 | kJ/mol | replicate means | `[D]` | ×4.184 of #31–33 |
| 40 | Ea in SI units, methional pH 6/7/8 | 102.1 / 106.6 / 81.4 | kJ/mol | replicate means | `[D]` | ×4.184 of #34–36 |
| 41 | Ea in SI units, 2-AT pH 8 | 132.9 | kJ/mol | replicate mean | `[D]` | ×4.184 of #37 |
| 42 | between-replicate spread, max | **21.8 %** (methional pH 6, 5.33 kcal/mol) | — | | `[D]` | §3.1 |
| 43 | `k` units | **ppm · min⁻¹** | — | derived, never printed | `[D]` | §3.6 |
| 44 | **`k` DMDS pH 6** | 0.0755, 0.0984, 0.1212, 0.5261, 1.1093 | ppm/min | 75/85/95/105/115 °C | `[D]` digitised | Fig 4 |
| 45 | **`k` DMDS pH 7** | 0.0640, 0.1349, 0.5032, 0.4552, 0.4803 | ppm/min | 75→115 °C | `[D]` digitised | Fig 4 |
| 46 | **`k` DMDS pH 8** | 0.0640, 0.1075, 0.3705, 0.6221, 2.6990 | ppm/min | 75→115 °C | `[D]` digitised | Fig 4 |
| 47 | **`k` methional pH 6** | 0.0101, 0.0292, 0.0553, 0.1825, 0.1966 | ppm/min | 75→115 °C | `[D]` digitised | Fig 5 |
| 48 | **`k` methional pH 7** | 0.0046, 0.0193, 0.0617, 0.0953, 0.1949 | ppm/min | 75→115 °C | `[D]` digitised | Fig 5 |
| 49 | **`k` methional pH 8** | 0.0052, 0.0315, 0.0358, 0.0625, 0.0970 | ppm/min | 75→115 °C | `[D]` digitised | Fig 5 |
| 50 | refit Ea DMDS pH 6/7/8 | 78.5 / 59.8 / 103.4 | kJ/mol | whole ladder | `[D]` | §5.1 |
| 51 | refit Ea methional pH 6/7/8 | 87.6 / 102.7 / 74.0 | kJ/mol | whole ladder | `[D]` | §5.1 |
| 52 | refit `r²` DMDS pH 6/7/8 | 0.884 / 0.806 / 0.965 | — | | `[D]` | §5.1 |
| 53 | refit `r²` methional pH 6/7/8 | 0.962 / 0.966 / 0.873 | — | | `[D]` | §5.1 |
| 54 | **per-leg Ea, DMDS pH 6** | 27.4, 22.9, **169.9**, 91.0 | kJ/mol | 75→85→95→105→115 | `[D]` | §5.2 |
| 55 | **per-leg Ea, DMDS pH 7** | 77.3, 144.3, **−11.6**, 6.6 | kJ/mol | | `[D]` | §5.2 |
| 56 | **per-leg Ea, DMDS pH 8** | 53.7, 135.7, 60.0, **179.1** | kJ/mol | | `[D]` | §5.2 |
| 57 | **per-leg Ea, methional pH 6** | 109.6, 70.1, 138.3, **9.1** | kJ/mol | | `[D]` | §5.2 |
| 58 | **per-leg Ea, methional pH 7** | 148.8, 127.2, 50.3, 87.3 | kJ/mol | | `[D]` | §5.2 |
| 59 | **per-leg Ea, methional pH 8** | **186.3**, 14.0, 64.6, 53.7 | kJ/mol | | `[D]` | §5.2 |
| 60 | Q₁₀ range, DMDS | **0.90 – 4.34** | — | 75–115 °C, all pH | `[D]` | §5.4 |
| 61 | Q₁₀ range, methional | **1.08 – 6.03** | — | 75–115 °C, all pH | `[D]` | §5.4 |
| 62 | DMDS ppm, pH 6, 95 °C | 1.61, 3.27, 3.91, 7.39, 11.84 | ppm | 20/40/60/80/100 min | `[M]`+`[D]` dig. | Fig 1 |
| 63 | DMDS ppm, pH 7, 95 °C | 3.72, 8.14, 13.43, 17.37, 19.76 | ppm | | `[M]`+`[D]` dig. | Fig 1 |
| 64 | DMDS ppm, pH 8, 95 °C | 11.90, 22.54, 25.13, 26.10, 24.75 | ppm | | `[M]`+`[D]` dig. | Fig 1 |
| 65 | methional ppm, pH 6, 95 °C | 1.17, 1.27, 3.09, 4.26, 4.39 | ppm | | `[M]`+`[D]` dig. | Fig 2 |
| 66 | methional ppm, pH 7, 95 °C | 0.82, 1.74, 2.81, 3.73, 4.23 | ppm | | `[M]`+`[D]` dig. | Fig 2 |
| 67 | methional ppm, pH 8, 95 °C | 1.33, 2.31, 2.90, 3.08, 3.06 | ppm | | `[M]`+`[D]` dig. | Fig 2 |
| 68 | 2-AT ppm, pH 6, 95 °C | n.d., n.d., n.d., 4.39, 23.64 | ppm | | `[M]`+`[D]` dig. | Fig 3 |
| 69 | 2-AT ppm, pH 7, 95 °C | **n.d. at every time point** | ppm | | `[M]` | Fig 3, p.131 |
| 70 | 2-AT ppm, pH 8, 95 °C | n.d., 2.29, 7.37, 16.96, 18.81 | ppm | | `[M]`+`[D]` dig. | Fig 3 |
| 71 | Fig 6 = mean of the 2 replicates | confirmed to ± 0.7 kcal/mol | — | | `[D]` | §3.7 |
| 72 | **S mass balance: 3 products / Met pool** | ★ **0.42 %** | — | peak yields, 95 °C | `[D]` | §5.6 |
| 73 | Figures 4–5 plot replicate | ★ **2nd replicate** | — | never stated by the paper | `[D]` | §4.3 |
| 74 | pyrazine Ea vs pH, prior work | highest at pH 5, lowest at pH 7 (pH 5–9) | — | Leahy & Reineccius, refs 8, 9 | `[C]` | p.135 |
| 75 | pyrazine rate vs T and pH, prior work | pyrazine, 2-methylpyrazine, dimethylpyrazine all increase with T and with pH | — | refs 8, 9 | `[C]` | p.129 |
| 76 | claimed methional Ea (Conclusions) | ★ **"17-19 kcal/mole" — CONTRADICTS Table I (19.39–27.07)** | kcal/mole | | `[F]` ⚠ C1 | p.135 |
| 77 | claimed overall Ea range (Conclusions) | ★ **"15 to 31 kcal/mole" — Table I is 15.02–32.71** | kcal/mole | | `[F]` ⚠ C2 | p.135 |

---

## §9. USABILITY VERDICTS

| block | verdict | reason (recorded so a later wave does not re-ingest a refused row) |
|---|---|---|
| **§1 system & conditions** (#1–#23) | **USE** | Fully specified, self-consistent (all five mol↔g pairs check out), raster-verified. The 400 mL / 0.075 mol / 0.5 mol system is reproducible as stated |
| **95 °C time list** (#16 vs #17) | ★ **REFUSE the pair as stated** | Methods and figures give different time axes (C4). **Either may be right; the paper does not let you decide.** If a later wave needs the 95 °C time base, it must be flagged as ±20 % ambiguous or the condition dropped |
| **pH as a controlled variable** | ★ **USE-Q — qualified as INITIAL pH only** | 0.1 M phosphate against 0.75 M amine; pH measured on every sample and never reported (N14). Never treat "pH 6/7/8" as a held condition, and never fit a pH-dependence from this paper as if pH were constant along a time series |
| **Table I Ea, DMDS, all 3 pH** (#31–33, #39) | ★ **USE-Q** | Independently reproduced from the authors' own `k` to within 4–12 % (§4.3). Qualification: **use the replicate MEAN with a ±12 % (pH 6, 7) to ±22 % (pH 8, and see #42) error bar**, and only as a *whole-ladder apparent* Ea over 75–115 °C — §5.2 shows the underlying ladder is badly non-linear |
| **Table I Ea, methional, all 3 pH** (#34–36, #40) | ★ **USE-Q** | Same: reproduced to within 4–9 %. Qualification: the pH-6 cell has a **21.8 % between-replicate spread** — the single largest disagreement in the paper. Use the mean, carry the spread |
| **Table I Ea, 2-acetylthiophene pH 8** (#37, #41) | ★ **REFUSE** | Unauditable — no `k` is published for 2-AT anywhere (N7). Underlying series is 0/0/0/4.4/23.6 ppm at pH 6 and detected at only 4/5 times at pH 8; the zero-order premise plainly fails (§4.4). **Recorded as PRIOR-ONLY directional evidence in §6.6 and nowhere else** |
| **Table I `r²` column** (#30) | **USE** | Identified as the Arrhenius `r²` by the audit (§4.2). Now unambiguous |
| **prose `r²` ranges** (#28, #29) | **REFUSE** | Two mutually incompatible statements for the same fits (C3), matching neither Table I nor each other |
| **"Ea of methional is 17-19 kcal/mole"** (#76) | ★ **REFUSE — flagged as an error in the source** | Contradicted by the paper's own Table I. This is the kind of quotable summary sentence that propagates; **kill it here** |
| **"activation energies ranged from 15 to 31"** (#77) | **REFUSE** | Wrong upper bound against Table I (C2) |
| **reaction order = zero** (#27) | ★ **STRUCTURAL** | Usable as a statement about how this paper's `k` were produced — i.e. `k` is a formation *slope*, so it has ppm/min units and no elementary-step meaning. **NOT usable as a mechanistic claim**: the authors' own Figs 1–3 show clear saturation (DMDS pH 8, methional pH 8) and clear acceleration (2-AT pH 6), and the "linear portion" window is undeclared (N17, §4.4) |
| **digitised `k` ladders** (#44–49) | ★ **RATIO-ONLY** | Digitisation validated (§4.3). Usable for `k(T₂)/k(T₁)`, Q₁₀ and Arrhenius slopes. **Not usable as absolute rates**: the axis is unitless (N8), the units are inferred (§3.6), the ppm basis is undefined (§2.3.2), and DMDS/methional sit on different effective sulfur scales (§2.4) |
| **derived whole-ladder Ea** (#50–53) | **USE-Q** | Independent recomputation; agrees with Table I. Qualification: prefer the *authors'* Table I values as the citable numbers and use these as the corroboration |
| **★ derived per-leg Ea** (#54–59) | ★ **USE-Q — this is the paper's most valuable output for wave K6a** | The leg structure is robust to the entire ± 9 px digitisation envelope (§5.3). Qualification: each leg rests on two digitised bars, so **± 15 % on the 75→85 and 85→95 legs and < ± 5 % on the top two legs**; and a leg Ea from a zero-order formation slope is a lumped formation-minus-consumption quantity (§5.5.1) |
| **Q₁₀** (#60, #61) | **USE-Q** | Same basis as the per-leg Ea (they are the same numbers). Never printed by the paper — always attribute as derived |
| **`ln(f·A)` / pre-exponential** | ★ **REFUSE** | Contaminated by the unmeasured semi-quant response factor (§2.4). There is no recoverable `A` in this paper |
| **ppm concentration series** (#62–70) | ★ **RATIO-ONLY / SHAPE-ONLY** | The ppm basis (aqueous vs extract) is undefined — a 50–100× ambiguity (§2.3.2). Usable for *shape* within a series and for the qualitative pH orderings; **never as an absolute concentration, never against another paper's µg/L or mg/kg** |
| **the S mass balance** (#72) | **USE-Q** | Robust in direction (the conclusion strengthens under either ppm basis) but the ppm basis is undefined, so quote it as **"< 0.5 %"**, not as 0.42 % |
| **pH orderings and reversals** (Figs 1–3) | **USE** | These are within-figure, within-run comparisons where the response factor cancels; and the authors state them in prose independently |
| **cited pyrazine results** (#74, #75) | ★ **PRIOR-ONLY** | Second-hand from Leahy & Reineccius (refs 8, 9). The chapter itself notes those studies *"did not have sufficient data to calculate the 95% confidence limits"*. **Go to the primaries; do not cite through Chan** |
| **cross-paper magnitude comparison of anything in this chapter** | ★ **REFUSE** | Both on the analytical basis (§2.4) and on the authors' own instruction: *"the kinetics observed in any study are very system dependent and cannot be readily compared"* (p.131) |
| **any claim about 120–140 °C** | ★ **REFUSE** | 115 °C is the ceiling (N26, §6.2) |

---

## §10. DECLARED GAPS

**What is missing from this chapter, and what would be needed to close each gap.**

| # | gap | what would close it |
|---|---|---|
| G1 | ★ **No table of rate constants.** `k` exists only as bars on two unitless axes | Nothing in this chapter. **A numeric `k` table would have to come from the primary thesis** — F. Chan's University of Minnesota MSc/PhD thesis under Reineccius (c. 1993–94), which is the obvious parent document and is not in the corpus |
| G2 | ★ **No rate constants at all for 2-acetylthiophene** (N7) | Same source as G1 |
| G3 | ★ **The analytical validation of the AED quantification is deferred to ref (14)**, Mistry, B. S.; Reineccius, G. A.; Jasper, B., *Sulfur Compounds in Foods*, ACS Symp. Ser. 564, **"in press"** | **ORDER REF 14.** It is a chapter in the *same book* (ACS Symp. Ser. 564) and would carry the AED response factors, detection limits and recovery that make §2.3's verdict either better or worse. **This is the single highest-value follow-up for this paper.** |
| G4 | **The non-sulfur volatiles and the rest of the network** are excluded and deferred to ref (12) | Chan, F.; Reineccius, G. A., *5th International Symposium on the Maillard Reaction* (Reineccius & Labuza, eds., RSC London, 1994) — **"in press"**. Would give the pyrazine/oxygen-heterocycle lane of the same runs |
| G5 | ★ **Reactant depletion was being measured and is not reported** — *"analysis is in progress to determine the amount of reactants remaining during heating"* (p.131) | A later Chan/Reineccius publication or the thesis (G1). **Would directly settle the plateau question that §5.6 can only bound** |
| G6 | ★ **The measured pH of every heated sample exists and is unreported** (N14) | Thesis (G1). Would convert "initial pH 6/7/8" into a real trajectory and is prerequisite to trusting any pH-resolved parameter |
| G7 | **The "linear portion" fitting windows** (N17) | Thesis (G1) |
| G8 | ★ **The 95 °C time-axis contradiction (C4)** cannot be resolved from within the chapter | Thesis (G1), or the raw chromatograms |
| G9 | **The 95 % confidence limits** argued from on p.135 (N12) | Thesis (G1). With n = 2 they may never have existed in a defensible form |
| G10 | ★ **No datum at or above 120 °C** (N26) — so Kang's switch-on window is untestable here | **A different paper.** The wave needs an aqueous sulfur ladder spanning 95–140 °C at controlled pH to bridge Chan's 115 °C ceiling to Kang's 120 °C floor |
| G11 | **No single-amino-acid control** — the four-amino-acid mixture is the authors' own stated main limitation (p.131) | A dedicated Met+glucose binary study; several exist in the corpus's sulfur lane and should be preferred wherever a clean methionine channel is needed |
| G12 | **No SI, no supplementary material, no data deposit** | 1994 ACS Symposium Series volumes have none. **There is nothing more to fetch for this object beyond G1 and G3** |

---

### Provenance summary for this dossier

* Every Table I cell, every molar loading, every schedule row, the abstract, the Methods, the
  kinetic-analysis paragraph, the Conclusions and the CAS line were read from **600 dpi raster
  renders**, not from the OCR text layer. Text-layer/raster agreement is recorded in §0.1.
* All six `k` ladders and all four concentration/Ea figures were **digitised at 600 dpi** with the
  calibration procedure, pixel anchors and error bars stated in §3.2–§3.7.
* The digitisation is **validated** against the authors' own printed Ea and `r²` in §4.3, including
  two `r²` matches to the third decimal place that are invariant to the digitisation's one free
  parameter.
* No number in §5, §6 or the `[D]` rows of §8 is printed anywhere in the paper.
