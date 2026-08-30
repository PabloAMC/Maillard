# Meng et al. — 2-methyl-3-furanthiol in fermented soy sauce: a three-temperature MFT ladder in a real food matrix

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THE CORPUS'S ONLY **THREE**-TEMPERATURE MFT LADDER — AND THEREFORE THE ONLY PAPER IN THIS WAVE THAT CAN ACTUALLY TEST FOR A SLOPE BREAK.

---

## ⚠️ HEADLINE — TWO THINGS THE BRIEF DID NOT EXPECT

> **① THE BRIEF SAYS "two-temperature MFT data: 95 and 120 °C". THAT IS WRONG — AND WRONG IN THE
> USEFUL DIRECTION. Table 3 carries FOUR levels: unheated, 80 °C, 95 °C and 120 °C, each at two
> hold times (5 and 20 min), for BOTH MFT and FFT, with n = 4 and SDs. This is a genuine
> three-point Arrhenius ladder with one degree of freedom — the only one in wave K6a — and it
> CAN detect curvature, which `feng2022_extraction.md` explicitly cannot.**
>
> **② WHAT THE CURVATURE SHOWS IS THE OPPOSITE SIGN FROM KANG'S SWITCH-ON. Apparent Eₐ FALLS
> as temperature rises, in all four series, with an identical residual pattern.**

| series (raw soy sauce A) | **80 → 95 °C** | **95 → 120 °C** | **3-point Eₐ** | **R²** | residual signs |
|---|---|---|---|---|---|
| **MFT, 5 min** | **× 5.500 · Eₐ 122.8** | **× 1.682 · Eₐ 25.0** | 61.1 | **0.839** | **− + −** |
| **MFT, 20 min** | × 2.833 · Eₐ 75.0 | × 2.691 · Eₐ 47.7 | 57.8 | 0.983 | **− + −** |
| **FFT (2FM), 5 min** | × 4.167 · Eₐ 102.8 | × 1.940 · Eₐ 31.9 | 58.1 | **0.900** | **− + −** |
| **FFT (2FM), 20 min** | × 3.714 · Eₐ 94.6 | × 3.500 · Eₐ 60.3 | 72.9 | 0.984 | **− + −** |
| *(reference: Kang 2026, TTCA/Cys)* | *100→120: Eₐ **6.9** (MFT), **5.8** (FFT)* | *120→140: Eₐ **97.9**, **69.0*** | *n/a* | *n/a* | *+ − +* |

**⇒ Kang's ladder is CONVEX in Arrhenius coordinates (flat then steep, Eₐ 6.9 → 97.9).
Meng's is CONCAVE (steep then flat, Eₐ 122.8 → 25.0). All four Meng series show the identical
`− + −` residual pattern about a straight line — a four-fold-replicated, same-sign curvature
signal.** ⚠️ **But see §7d: the 95 °C point is time-saturated (5 → 20 min gain of only × 1.03),
and that alone can manufacture exactly this curvature. The curvature is real in the data and
confounded in its interpretation.**

**On the specific 95 → 120 °C leg — the closest thing Meng has to Kang's low leg:**

| | **Meng 95 → 120 °C, MFT** | **Kang 100 → 120 °C, MFT** | replicates? |
|---|---|---|---|
| 5 min | **× 1.682** (Eₐ 25.0) | × 1.12 (Eₐ 6.9) | **partly — same order, 1.5 × steeper** |
| 20 min | × 2.691 (Eₐ 47.7) | × 1.12 | ✗ |
| | **Meng 95 → 120 °C, FFT** | **Kang 100 → 120 °C, FFT** | |
| 5 min | **× 1.940** (Eₐ 31.9) | × 1.10 (Eₐ 5.8) | **partly** |
| 20 min | × 3.500 (Eₐ 60.3) | × 1.10 | ✗ |

**⇒ VERDICT: Kang's flat low leg is CORROBORATED IN DIRECTION AND ORDER OF MAGNITUDE at short
hold, CONTRADICTED at long hold, and the flatness is time-dependent in a way Kang never tested.
Kang's *switch-on* — Eₐ rising with T — is CONTRADICTED by the only three-point ladder available,
with the caveat of §7d.** Full argument in §7.

**Provenance codes:** **[M]** measured and printed · **[C]** cited from another work ·
**[F]** fitted by the authors · **[D]** derived by this wave, never printed · **[NEG]** verified
negative · **[!]** integrity flag.

---

## §0. IDENTITY — **CONFIRMED, with ONE flag: THE CITE YEAR**

| field | value | how verified |
|---|---|---|
| file on disk | `data/articles/meng2017.pdf`, **359,316 bytes**, 6 pages, PDF 1.4 | `ls`, `pdfinfo` |
| SHA-256 | `b88900d73b621919c56fe0aaa89505894ef88f4a6c1279dfe7dff864a9a17129` | `shasum -a 256` |
| **title** | ***"Contribution of 2-methyl-3-furanthiol to the cooked meat-like aroma of fermented soy sauce"*** | p. 1 title block, raster-verified; matches PDF `Title` metadata exactly |
| **authors** | **Qi Meng¹, Riho Kitagawa², Miho Imamura³, Hiroshi Katayama³, Akio Obata³, Etsuko Sugawara¹\*** | p. 1, raster-verified |
| affiliations | ¹ **Faculty of Education, Iwate University, Morioka, Japan** · ² **Graduate School of Agriculture, Iwate University, Morioka, Japan** · ³ **Research and Development Division, Kikkoman Corporation, Noda, Japan** | p. 1 |
| corresponding author | **Etsuko Sugawara**, `etsukos@iwate-u.ac.jp` | p. 1 footnote |
| **DOI** | **`10.1080/09168451.2016.1238295`** | p. 1 and the T&F cover page |
| **journal** | ***Bioscience, Biotechnology, and Biochemistry***; ISSN 0916-8451 (print) / 1347-6947 (online) | cover page, running heads |
| **dates** | **Received April 28, 2016 · Accepted September 7, 2016 · Published online 3 October 2016** | p. 1, raster-verified; cover page |
| **© / society** | **© 2016 Japan Society for Bioscience, Biotechnology, and Agrochemistry** | p. 1 |
| producer | **Arbortext Advanced Print Publisher 9.1.531/W Unicode** → iText 4.2.0. CreationDate 2016-10-03; ModDate 2016-10-21 | `pdfinfo` |
| provenance stamp | *"Download by: [Ryerson University Library] Date: 20 October 2016"*; *"Article views: 26"* | cover page |

### ⚠️ **THE CITE-YEAR FLAG — REPORTED AS ASKED, AND NOT RESOLVED FROM THIS PDF**

**This document self-identifies as 2016, everywhere, and carries NO volume, NO issue and NO
page numbers.**

- PDF `Subject` metadata: *"Bioscience, Biotechnology, and Biochemistry, **2016**. doi:10.1080/09168451.2016.1238295"*
- Running head on p. 1: *"Bioscience, Biotechnology, and Biochemistry, **2016**"*
- T&F cover-page citation block: *"Qi Meng, … & Etsuko Sugawara (**2016**): Contribution of 2-methyl-3-furanthiol…, Bioscience, Biotechnology, and Biochemistry, **DOI: 10.1080/09168451.2016.1238295**"* — **the DOI stands in place of a volume/page range**, which is the T&F "Latest Articles" signature.
- Page footers read `2`, `3`, `4`, `5` — **article-relative page numbers, not journal pagination.**
- The DOI slug itself contains **`2016`**.

**⇒ This is the online-first / "Latest Articles" version of record, published 3 Oct 2016, before
issue assignment.** The repo filename `meng2017.pdf` implies the journal later assigned it to a
2017 issue — which is the normal T&F pattern for an October online-first — **but that final
volume/issue/page assignment is NOT present anywhere in this file and I have not verified it.**

**RECOMMENDATION for the repo:** cite as
> **Meng, Q.; Kitagawa, R.; Imamura, M.; Katayama, H.; Obata, A.; Sugawara, E. "Contribution of
> 2-methyl-3-furanthiol to the cooked meat-like aroma of fermented soy sauce." *Biosci.
> Biotechnol. Biochem.*, published online 3 October 2016. DOI 10.1080/09168451.2016.1238295.**

and **verify the final volume/issue/pages against Crossref before any bibliography is frozen.**
⚠️ **Do not assert "2017, 81(1), 168–172" or any other volume from this dossier — that number is
not in the document.** The filename stem `meng2017` should be kept as the repo key (it is
already referenced in `data/lit/slr_incorporation_matrix.json` and
`data/lit/benchmark_intake_registry.json`), but the *printed* citation must not silently gain a
volume this file cannot support. **Everything else — title, all six authors, all three
affiliations, the DOI, the received/accepted dates — matches the expected identity exactly.**

**⚠️ No funding statement appears anywhere in the paper.** `[NEG]`
**Competing interests, verbatim:** *"No potential conflict of interest was reported by the
authors."* (p. 5)
**Author contributions, verbatim:** *"Qi Meng (QM) and Etsuko Sugawara designed the experiment.
QM and Riho Kitagawa performed the experiments. Miho Imamura designed the sensory analysis and
performed it by herself. QM wrote the manuscript. All authors contributed to the development of
the manuscript. All authors have read and approved the final manuscript."* (p. 4)
**⚠️ Note the industry involvement:** three of six authors are at **Kikkoman Corporation**, which
also **prepared the soy sauce samples A, B and BH** used for every heating experiment (§1). Not a
declared conflict, but a material one for a paper whose conclusion is that *fermented* soy sauce
generates a desirable aroma that *acid-hydrolysed* soy sauce does not. `[!]`

**⚠️ There is no Supporting Information.** No SI is referenced anywhere; every table and figure
is in the 6-page PDF. Nothing is missing to a supplement — the gaps in §8 are permanent.

---

## §0b. ⚠️ THE μ → m RASTER CHECK (Amendment 4) — **DISCHARGED. MENG IS CLEAN.**

**Binding requirement satisfied. Every unit that carries a number used anywhere in this dossier
was verified against a raster render, not the text layer.** Method:
`pdftoppm -r 500/600/800 -png -f N -l N meng2017.pdf out`, cropped to the artefact and read as an
image, glyph by glyph. This file is **Arbortext Advanced Print Publisher**-produced, the same
family as the ACS PDFs where the hazard has fired before — so the check was not a formality.

| artefact | text-layer says | **raster says** | verdict |
|---|---|---|---|
| **Table 1 column super-header** (p. 3) | `Concentration (μg/L)` | **`Concentration (μg/L)` — true Greek mu, 600 dpi** | ✅ **NOT corrupted** |
| **Table 1 threshold column** (p. 3) | `Threshold¹¹˒¹⁷⁾ (ng/L)` | **`Threshold (ng/L)`** — **`n`, not `μ`; the thresholds really are nanograms** | ✅ **NOT corrupted** |
| **Table 2 row label** (p. 4) | `…for addition to soy sauce (μg/L)` | **`(μg/L)`** | ✅ **NOT corrupted** |
| ★ **Table 3 column headers** (p. 4) | `2M3F (μg/L)` and `2FM (μg/L)` | ★ **`2M3F (μg/L)` and `2FM (μg/L)`** | ✅ **NOT corrupted** |
| ★ **Table 4 column headers** (p. 4) | `Raw soy sauce (B, μg/L)` / `Heat-treated soy sauce (BH, μg/L)` | ★ **`(B, μg/L)` / `(BH, μg/L)`** | ✅ **NOT corrupted** |
| **Abstract** (p. 1) | `with the addition of only 0.2 μg/L of 2M3F` | ★ **`0.2 μg/L`, 800 dpi, bold face** | ✅ **NOT corrupted** |
| **Methods, internal standard** (p. 2) | `10.72 μg/L of 4-methoxy-2-methyl-2-mercaptobutane` | **`10.72 μg/L`** | ✅ **NOT corrupted** |
| **Methods, concentration step** (p. 2) | `concentrated to 25 μL under a nitrogen flow` | **`25 μL`** | ✅ **NOT corrupted** |
| **Methods, column film** (p. 2, twice) | `60 m × 0.25 mm (i.d.), 0.25 μm film thickness` | **`0.25 μm`** | ✅ **NOT corrupted** |
| **Results, standard curve** (p. 3) | `At concentrations ranging from 0 to 20 μg/L` | **`0 to 20 μg/L`** | ✅ **NOT corrupted** |
| **Results, sensory** (p. 3) | `0.5 μg/L`, `0.2, 0.4, or 0.8 μg/L` | **all four mus true** | ✅ **NOT corrupted** |
| **Discussion** (p. 4) | `greater than 10 μg/L of 2M3F was formed in FSS` | **`10 μg/L`** | ✅ **NOT corrupted** |
| **threshold citation** (p. 3) | `the threshold for 2M3F … is 4 ng/L` | **`4 ng/L` — `n`, confirmed** | ✅ **NOT corrupted** |

**Per-table basis of verification.**
- **Table 1** — raster, p. 3 at 600 dpi: caption, the `Concentration (μg/L)` super-header, both
  sub-headers (`Raw soy sauce` / `Heat-treated soy sauce`), the `Threshold (ng/L)` column head,
  **all four compound rows** and the footnote. Every raster reading matched the text layer to the
  last decimal. **Verified end-to-end.**
- **Table 2** — raster, p. 4 at 600 dpi: caption, the `(μg/L)` row label, **all three data rows
  and the p-value row.** **Verified end-to-end.**
- **Table 3** — raster, p. 4, in two overlapping crops at 600 dpi: caption, both `(μg/L)` column
  headers, the `Temperature / 5 min / 20 min` sub-header row, **all four data rows (unheated,
  80, 95, 120 °C) and all five footnotes.** **Verified end-to-end — this is the paper's key
  artefact and every one of its 14 numbers was read off a raster.**
- **Table 4** — raster, p. 4 at 600 dpi: caption, both `(μg/L)` column groups, all four
  sub-headers, **all four data rows and all three footnotes.** **Verified end-to-end.**
- **Abstract, Methods and Results prose** — raster at 500 and 800 dpi for every clause carrying a
  unit, as tabulated above.

**⇒ Despite the Arbortext producer, this file embeds a real Greek `μ` throughout. No μ→m
correction is required anywhere in Meng. The hazard was checked exhaustively and did not fire.**
⚠️ The asymmetry noted in `kang2026_SI_extraction.md` §0b still stands and still binds: **only
`Kang2026.pdf` is corrupted**, so any cross-comparison of Meng's numbers with Kang's must use the
SI dossier's corrected units.

**⚠️ ONE GENUINE TYPOGRAPHIC ERROR, unrelated to μ, found on the raster** `[!]`: p. 2, Methods:
*"7.7 mM of ribose or 47 mM of glucose and **1.8 m M** of cysteine were dissolved in a phosphate
buffer"* — a spurious space inside `mM`. The **same sentence three lines earlier** correctly reads
`1.8 mM`, and **Table 4's footnote a** independently confirms **`cysteine (1.8 mM)`**. **The value
is 1.8 mmol L⁻¹; the space is a typesetting slip, not a unit error.**

---

## §1. SYSTEM AND CONDITIONS `[M]`

### 1a. The samples — **five distinct materials, and they are NOT interchangeable**

| label | what it is | provenance | used in |
|---|---|---|---|
| **five raw FSS** | five **raw (unpasteurised) fermented** soy sauces | purchased, local market, Japan | Table 1 survey |
| **five heat-treated FSS** | five **pasteurised fermented** soy sauces | purchased, local market, Japan | Table 1 survey |
| ★ **A** | **raw** soy sauce, prepared for this study | **Kikkoman Corporation** | ★ **Table 3 — the temperature ladder** |
| **B** | **raw** soy sauce, prepared for this study, **different lot from A** | **Kikkoman Corporation** | Table 4 |
| **BH** | **B, heated by Kikkoman for pasteurisation** | **Kikkoman Corporation** | Table 4 |
| **three ASS** | **acid-hydrolysed-vegetable-protein-mixed** soy sauce | purchased, **USA and Brazil** | the negative control |

> *"For heating experiments, a raw soy sauce sample (A), raw sample (B) and heat-treated soy
> sauce sample (BH) were prepared by Kikkoman Corporation for this study. **Prepared A and B had
> different lot numbers.** BH was prepared by heating B for pasteurization by Kikkoman
> Corporation."* — p. 2

**⚠️ [!] Tables 3 and 4 use DIFFERENT raw soy sauces (A vs B) and the paper says so explicitly.**
At the shared condition — **120 °C, 5 min, no addition** — Table 3 gives **11.1 μg L⁻¹** of MFT
in A and Table 4 gives **13.5 μg L⁻¹** in B: **a 21.6 % discrepancy between two nominally
identical treatments of two lots of raw soy sauce.** ⇒ **That 21.6 % is the honest lot-to-lot
reproducibility floor for this study**, and it is *larger* than most of the printed SDs
(11.1 ± 0.4 is 3.6 % relative). **Do not pool Tables 3 and 4, and do not treat the printed SDs as
capturing between-lot variability.** `[!]`

### 1b. Heating conditions — **the temperature ladder**

> *"In the first series of experiments, 200 mL of A was loaded into a glass cylinder and heated
> for 5 or 20 min at 80 °C or 95 °C in a dry stirring bath and at 120 °C in a laboratory
> autoclave. Next, the samples were cooled in water to 25 °C after attaining the final
> temperature."* — p. 2

| item | value | anchor |
|---|---|---|
| **substrate** | **raw fermented soy sauce A, neat** (undiluted for heating) | p. 2 |
| **volume** | **200 mL** in a **glass cylinder** | p. 2 |
| ★ **temperatures** | ★ **unheated control, 80 °C, 95 °C, 120 °C — FOUR levels, THREE of them heated** | p. 2, Table 3 |
| ★ **times** | ★ **5 and 20 min at every heated temperature** | Table 3 |
| **apparatus, 80 and 95 °C** | **dry stirring bath** — ⚠️ **no statement that the cylinder is sealed; a glass cylinder in a dry bath is presumptively OPEN and at atmospheric pressure** `[!]` | p. 2 |
| **apparatus, 120 °C** | ★ **laboratory autoclave — SEALED and PRESSURISED** | p. 2 |
| ⚠️ **vessel-type discontinuity** | ★ **THE 120 °C POINT IS TAKEN IN A DIFFERENT VESSEL AND A DIFFERENT PRESSURE REGIME FROM THE 80 AND 95 °C POINTS.** See §6b — this is a first-order confound on the ladder. | p. 2 |
| **quench** | *"cooled in water to 25 °C after attaining the final temperature"* — ⚠️ **water cooling to 25 °C, not an ice bath; no cooling time given, and for 200 mL of liquid this is slow** `[!]` | p. 2 |
| ⚠️ **come-up time** | ★ **NOT STATED for any temperature.** *"heated for 5 or 20 min … after attaining the final temperature"* implies the clock starts at temperature — **but for 200 mL in a dry bath and in an autoclave the come-up times differ by minutes and are unequal across the three levels.** At a nominal 5 min hold this is a large fractional error. `[NEG]` | — |
| **pH** | ★ **"The pH of the FSS sample was approximately 4.7."** — a single approximate value; **no buffer, no adjustment, no drift measurement, no per-sample pH** `[M]`/`[NEG]` | p. 4 |
| **replication** | **n = 4** — *"Values correspond to the mean ± SD (standard deviation) of four analyses"* | Table 3 footnote c, Table 4 footnote d |
| statistics | **Student's t test** for the heating and buffer-model concentrations; **χ² test** for the sensory panel | p. 2 |

**⚠️ The matrix is a real food and its composition is essentially unmeasured.** Soy sauce is
≈ 16–18 % NaCl, 2–3 % ethanol (the paper states *"Alcoholic fermentation generates 2–3 % ethanol
during the aging of soy sauce"*, p. 3), pH ≈ 4.7, and contains an unmeasured pool of free amino
acids and reducing sugars. **None of the precursor concentrations is measured. `[NEG]`** This is
the defining limitation of the paper for the repo (§10).

### 1c. The precursor-addition arms

> *"In the second series of experiments, 7.7 mM of ribose or 47 mM of glucose and 1.8 mM of
> cysteine were added to B and BH, respectively, followed by immediately heating at 120 °C for
> 5 min. In addition, 7.7 mM of ribose or 47 mM of glucose and 1.8 mM of cysteine were dissolved
> in a phosphate buffer (0.5 mol/L, 180 mL) at pH 5.0 and thermally treated for 5 min at 120 °C
> in a laboratory autoclave. The concentrations of ribose, glucose, and cysteine added in the
> second series of experiments represented the final concentrations in the sample solutions."*
> — p. 2

| item | value | anchor |
|---|---|---|
| **added ribose** | **7.7 mmol L⁻¹**, final concentration | p. 2, Table 4 footnote a |
| **added glucose** | **47 mmol L⁻¹**, final concentration | p. 2, Table 4 footnote a |
| **added cysteine** | **1.8 mmol L⁻¹**, final concentration | p. 2, Table 4 footnote a |
| ⚠️ **sugar : Cys molar ratio** | **ribose arm 4.28 : 1** · **glucose arm 26.1 : 1** — ★ **the two sugar arms are NOT matched on molar loading; glucose is at 6.1 × the ribose molarity.** Any ribose-vs-glucose comparison in this paper is confounded with a 6.1-fold concentration difference. `[!]` `[D]` | `[D]` from p. 2 |
| **volume** | **180 mL** (soy sauce arms and buffer arm alike) | p. 2, Table 4 footnote a |
| **temperature / time** | **120 °C, 5 min, laboratory autoclave** — all arms | p. 2 |
| ★ **THE BUFFER CONTROL** | ★ **phosphate buffer, 0.5 mol L⁻¹, pH 5.0, 180 mL** — **the only buffered system in the paper, and the only one with a stated buffer molarity** | p. 2 |
| ⚠️ **the soy-sauce arms are UNBUFFERED** | pH ≈ 4.7, native, no adjustment | p. 4 |
| **pH mismatch** | ★ **buffer arm at pH 5.0 vs soy sauce at ≈ 4.7 — a 0.3-unit offset that is never reconciled** `[!]` | pp. 2, 4 |

### 1d. Sensory analysis conditions

| item | value | anchor |
|---|---|---|
| design | **paired comparison test**, two samples per pair | p. 2 |
| base sample | **heat-treated soy sauce containing 0.5 μg L⁻¹ of 2M3F** (chosen to match the lowest of the five heat-treated FSS) | pp. 2–3 |
| spiked levels | **+0.2, +0.4, +0.8 μg L⁻¹** of 2M3F | p. 2 |
| panel | ★ **34 female panelists, aged 27–55**, selected per **ISO 8586** and by soy-sauce sensitivity (ref. 14) | p. 2 |
| presentation | **15 mL in a 60 mL black plastic cup**, 3-digit code | p. 2 |
| task | *"to smell both soy sauce samples and to choose the one that exhibited stronger cooked meat-like aroma"* | p. 2 |
| statistic | **χ² test** | p. 2 |
| ⚠️ | **All-female panel, no rationale given, no gender-balance statement.** `[!]` **No replicate presentations, no order-balancing statement, no forced-choice/no-difference option stated.** `[NEG]` | |

---

## §2. ANALYTICAL METHOD AND QUANTIFICATION BASIS — **REQUIRED SECTION**

### 2a. Extraction and instrumentation

| item | value | anchor |
|---|---|---|
| ★ **internal standard** | ★ **4-methoxy-2-methyl-2-mercaptobutane** (Oxford Chemicals, Bound Brook NJ), added at ★ **10.72 μg L⁻¹** | p. 2, raster-verified |
| sample preparation | **sodium sulfite (1 g)** + **0.5 mL ethyl acetate** added to **200 mL of a 40 % soy sauce solution**, then IS added | p. 2 |
| **⚠️ the sample is DILUTED to 40 %** | ★ **All reported concentrations refer to the 40 % solution basis unless the authors back-corrected — and they never say whether they did.** See §2c item 2. `[!]` | p. 2 |
| Na₂SO₃ role | not stated; standard practice is **antioxidant, to prevent thiol oxidation during workup** `[D]` inference | p. 2 |
| extraction | **stirred twice with 5 mL dichloromethane at 1000 rpm for 5 min**; organic phase dried over **anhydrous sodium sulfate**; **concentrated to 25 μL under nitrogen** in a glass tube | p. 2 |
| ⚠️ **NOT SPME** | ★ **This is a liquid–liquid solvent extraction with a concentration step, not headspace SPME.** Unlike Feng 2022 and Kang 2026, **there is no headspace-partition step and no partition-coefficient concern** — but there is a solvent-recovery concern instead, and **no recovery is reported.** `[NEG]` | p. 2 |
| **GC-O** | Shimadzu **GC-17A**; DB-XLB, **60 m × 0.25 mm i.d., 0.25 μm**; oven **40 °C hold 10 min → 220 °C at 3 °C/min → hold 1 min**; He **1.5 mL/min**, linear velocity **26.5 cm/s**; injector **210 °C**, FID **230 °C**; effluent split **1 : 1** FID / sniffing port | p. 2 |
| **GC-MS** | Shimadzu **QP-2010**; **same DB-XLB column and same oven programme**; **EI 70 eV**, ion source **200 °C**; **splitless, splitless time 1 min** | p. 2 |
| **acquisition mode** | ★ **SIM.** 2M3F confirmed on m/z **45, 85, 86, 113, 114**; ★ **quantifier ions: m/z 114 for 2M3F AND 2FM, m/z 134 for ET2MP, m/z 124 for BM, m/z 134 for the internal standard** | p. 2 |
| identification | **odour quality + retention time vs the authentic reference compound**, plus the SIM spectrum | pp. 2–3 |
| **detection anchor** | **cooked meat-like aroma sniffed at 20.6 min** on DB-XLB, matched to the 2M3F reference | p. 2 |
| ⚠️ **no retention indices** | ★ **No RI, no Kovats index, no n-alkane series anywhere in the paper.** Identification rests on absolute retention time against a standard on a single column. `[NEG]` | — |

**⚠️ [!] The internal standard and ET2MP share the same quantifier ion, m/z 134.** 4-methoxy-2-
methyl-2-mercaptobutane (C₆H₁₄OS, M = 134) and ethyl 2-mercaptopropionate (C₅H₁₀O₂S, M = 134)
are **nominal-mass isobars**, both quantified on m/z 134 in the same SIM run on the same column.
This is safe **only if they are fully chromatographically resolved**, which the paper does not
demonstrate — no chromatogram is shown and no resolution figure is given. **⇒ ET2MP's numbers
carry a co-elution risk that MFT's and FFT's do not.** ET2MP is peripheral to this wave, so the
consequence is limited, but flag it: **do not use ET2MP from this paper.**

### 2b. **THE QUANTIFICATION BASIS — one bolded sentence, as required**

> *"Hence, a standard curve for a 2M3F assay is prepared by GC-MS in the SIM mode. At
> concentrations ranging from 0 to 20 μg/L, the regression equation was linear in the form
> **y = 0.5002x, r = 0.999** (y = 2M3F, μg/L; x represents the height of the thiol peak/the
> height of the internal standard peak)."* — p. 3, raster-verified

**⇒ FOR 2M3F, MENG'S QUANTIFICATION IS ABSOLUTE: an authentic 2M3F standard, an external
calibration curve spanning 0–20 μg L⁻¹ with r = 0.999, and a measured response factor of 0.5002
μg L⁻¹ per unit peak-height ratio against the 10.72 μg L⁻¹ internal standard — so no response
factor is assumed equal to 1 for the paper's headline analyte. FOR 2FM, ET2MP AND BM, HOWEVER,
NO CALIBRATION CURVE IS PRINTED IN THIS PAPER; their concentrations rest on the method of
reference 8 (Meng, Kakuta & Sugawara 2012) and are, from the evidence available here, effectively
SINGLE-IS SEMI-QUANT WITH AN UNVERIFIED RESPONSE FACTOR.**

**Consequences, spelled out:**
- **2M3F is Tier A** and is the only analyte in this paper that is.
- **2FM (= FFT) is Tier B *within this document*.** ⚠️ **This matters for §7, because the Kang
  cross-check needs FFT.** The paper does assert *"the 2FM, ET2MP, and BM concentrations in this
  study were in agreement with those previously reported"* (p. 3), which is a consistency claim
  rather than a calibration, and it applies only to Table 1's survey — **not to the heated
  samples of Table 3, which are 5–15 × above the survey range.**
- **⚠️ The curve is forced through the origin: `y = 0.5002x` has NO intercept.** A zero intercept
  is an assumption, not a fitted result, and it is the assumption most likely to fail at the
  bottom of the range — exactly where the 0.3 ± 0.1 μg L⁻¹ unheated control sits. **Treat the
  unheated row of Table 3 as an upper-bounded estimate, not a measurement.** `[!]`
- **⚠️ The reported r = 0.999 is `r`, not `R²`** (r² = 0.998). Minor, but the paper's own
  notation.
- **⚠️ The range is 0–20 μg L⁻¹ and Table 3's largest value is 18.3.** **In range — just.** The
  headroom is 9 %. **Nothing in this paper is extrapolated beyond its calibration**, which is
  better than either Feng 2022 (curves unavailable) or Kang 2026 (range unstated).

### 2c. **THIS WAVE'S RULE, APPLIED EXPLICITLY**

The wave rule: *a single-IS semi-quant ladder is still usable for activation energies as
within-study SHAPE, because a constant unknown response factor cancels in a ratio and therefore
in an Arrhenius slope.*

**Applied here.** For 2FM the reported concentration is `C = f · (h_2FM / h_IS)` with `f`
unverified. The 120 °C / 95 °C ratio is
`C₁₂₀/C₉₅ = (h/h_IS)₁₂₀ / (h/h_IS)₉₅` — **`f` cancels exactly.** So **the 2FM fold-changes and
apparent Eₐ in §6 are valid even though its absolute scale is not established here.** For 2M3F
the same holds, and its absolute scale *is* established. **⇒ Every Arrhenius quantity in this
dossier is licensed for both analytes.**

**What that does NOT license — spelled out:**
1. **Not absolute 2FM yields.** Use 2FM's μg L⁻¹ values as a *shape*, not as a level. If an
   absolute 2FM number is needed, it must come from ref. 8, which is not on disk. `[NEG]`
2. **Not cross-paper magnitude comparison.** Meng's MFT at 120 °C (11.1 μg L⁻¹ in soy sauce) and
   Kang's (1.388 μg L⁻¹ in a 10 mmol L⁻¹ TTCA model) and Feng's (7.22 in a 20 mmol L⁻¹ ARP model)
   are three different extractions from three different matrices at three different precursor
   loadings. **Those ratios are meaningless and must never enter a benchmark.** Only the
   *within-paper* temperature ratios are comparable across the three papers.
3. **Not any compound whose response factor could plausibly change with matrix — and here that is
   a live, not a theoretical, concern.** The extraction is from a **40 % soy sauce solution**
   carrying ≈ 7 % NaCl and 1 % ethanol into dichloromethane. **Salting-out raises DCM extraction
   efficiency for volatile thiols, and the effect is analyte-specific.** The standard curve was
   presumably built in a clean solvent — the paper does not say — so **the response factor
   measured on standards may not be the response factor operating in the sample.** ⚠️ **This
   applies to 2M3F too, and is the reason its Tier A status covers the *ladder* but not
   necessarily the *absolute level in soy sauce*.** ⇒ Marked **USE-Q** rather than **USE** in §10.
   **Critically, it cancels in the temperature ratio** — all three heated samples and the
   unheated control come from the same 40 % soy sauce, so the matrix is constant down each column
   of Table 3.
4. **Not a rate.** See §6a.

### 2d. Everything about the calibration that is **not** stated `[NEG]`

**No LOD. No LOQ. No number of calibration levels. No replicate structure on the curve. No
recovery. No spike experiment. No matrix-matched calibration statement. No blank. No
carry-over check. No chromatogram.** For a paper whose central claim rests on a **negative**
(2M3F "not detected" in ASS and in the phosphate buffer), **the absence of an LOD is the single
most consequential omission in the document** — see §5b.

---

## §3. EVERY TABLE, CELL BY CELL

### 3a. **TABLE 1 — the ten-sample survey of volatile thiols in soy sauce** `[M]`

**Caption, raster-verified verbatim:** *"Table 1. Concentrations of volatile thiols in FSS
samples."* — p. 3.
**Column structure, raster-verified:** `Compound | Threshold¹¹˒¹⁷⁾ (ng/L) | Concentration (μg/L)
[ Raw soy sauce: Rangeᵃ | Average ‖ Heat-treated soy sauce: Range | Average ]`
**Footnote, verbatim:** *"ᵃRange of volatile thiol concentrations in five soy sauce samples."*

| compound | **threshold, ng L⁻¹** `[C]` | **raw FSS range** | **raw avg** | **heat-treated FSS range** | **HT avg** |
|---|---|---|---|---|---|
| ★ **2M3F** (2-methyl-3-furanthiol = **MFT**) | **4** | **0.4–1.0** | **0.7** | **0.5–1.9** | **1.3** |
| ★ **2FM** (2-furanmethanethiol = 2-furfurylthiol = **FFT**) | **0.4** | **0.8–1.9** | **1.3** | **1.2–2.8** | **2.0** |
| ET2MP (ethyl 2-mercaptopropionate) | **500** | 5.4–27.4 | 13.6 | 13.7–40.9 | 22.8 |
| BM (benzenemethanethiol) | **0.3** | 0.1–0.8 | 0.3 | 0.3–0.7 | 0.6 |

**⚠️ The thresholds are `[C]`, not `[M]`, and their basis matters:** they are cited to refs 11
and 17 (Tominaga & Dubourdieu 2006; Tominaga, Guimbertau & Dubourdieu 2003) and are, as the paper
itself states, *"the threshold for 2M3F in a **model dilute alcohol solution** is 4 ng/L"*
(p. 3). ★ **These are wine-model thresholds transplanted into a soy-sauce paper, and the paper
says so.** It also concedes: *"the perception threshold of 2M3F in soy sauce could not be
directly determined"* (p. 3). **⇒ Threshold values from this paper are PRIOR-ONLY. `[C]`**

**Derived odour activity values `[D]` — the paper prints none:**

| compound | **OAV, raw FSS** | **OAV, heat-treated FSS** | **OAV at 120 °C / 20 min (Table 3)** |
|---|---|---|---|
| **2M3F** | **100 – 250** | **125 – 475** | ★ **4 575** |
| **2FM** | **2 000 – 4 750** | **3 000 – 7 000** | ★ **45 500** |
| ET2MP | 10.8 – 54.8 | 27.4 – 81.8 | — |
| BM | 333 – 2 667 | 1 000 – 2 333 | — |

**⇒ Every thiol in every soy sauce is far above its cited threshold. 2FM dominates on OAV by
≈ 20 × over 2M3F throughout — a fact the paper never states and which materially qualifies its
title claim.** `[D]` **⚠️ These OAVs inherit the wine-model thresholds and are therefore
PRIOR-ONLY, not evidence of sensory dominance in soy sauce.**

**⚠️ [!] Table 1 gives no n, no SD and no per-sample values** — only a range and an average over
five samples. **Whether the "range" is over samples or over analytical replicates is not stated**
(footnote a implies over samples). **No sample identities, no brands, no fermentation details.**
`[NEG]`

**⚠️ [!] The two "sources" are confounded.** The five raw and five heat-treated soy sauces are
**different commercial products**, not paired before/after samples of the same product. The
raw → heat-treated increase (2M3F 0.7 → 1.3, × 1.86) is therefore **a between-product
comparison, not a heating effect.** The paper's controlled heating evidence is Tables 3 and 4,
not Table 1. **Do not read Table 1 as a heating experiment.**

### 3b. **TABLE 2 — the sensory paired-comparison test** `[M]`

**Caption, raster-verified verbatim:** *"Table 2. Number of panelists who recorded stronger
cooked meat-like aroma soy sauce (N = 34)."* — p. 4.

| | **+0.2 μg L⁻¹** | **+0.4 μg L⁻¹** | **+0.8 μg L⁻¹** |
|---|---|---|---|
| chose **soy sauce** (unspiked, base = 0.5 μg L⁻¹ 2M3F) | **9** | **2** | **4** |
| chose **soy sauce + 2M3F** | **25** | **32** | **30** |
| **p-value (χ² test)** | **< 0.01** | **< 0.001** | **< 0.001** |

**Audit `[D]`:** all three columns sum to 34 ✔. χ² recomputed against a 17 : 17 null:
**+0.2 → χ² = 7.53, p = 0.0061** ✔ consistent with "< 0.01"; **+0.4 → χ² = 26.47, p = 2.7 × 10⁻⁷**
✔; **+0.8 → χ² = 19.88, p = 8.2 × 10⁻⁶** ✔. **All three printed p-values reproduce.**

**⚠️ [!] The dose–response is NON-MONOTONE:** 25 / 32 / 30 correct out of 34 at +0.2 / +0.4 /
+0.8 μg L⁻¹. **The highest dose performs worse than the middle dose.** With n = 34 the 32-vs-30
difference is well inside binomial noise (± 2.7 at p ≈ 0.9), **so the honest reading is that the
test saturates by +0.4 μg L⁻¹ and carries no dose information above it.** `[D]`

**⇒ What Table 2 does establish, and it is a clean result: a 0.2 μg L⁻¹ increment of 2M3F on top
of a 0.5 μg L⁻¹ background is detectable by a trained panel at p < 0.01.** That is a
**just-noticeable-difference in a real matrix**, and it is the only sensory number in this wave
that is measured *in soy sauce* rather than transplanted from a wine model. **⇒ The paper's own
data therefore imply an in-matrix JND of ≈ 0.2 μg L⁻¹ = 200 ng L⁻¹, which is 50 × the cited
4 ng L⁻¹ wine-model absolute threshold.** `[D]` ⚠️ A JND on a 0.5 μg L⁻¹ background is not an
absolute threshold and the two must not be conflated — but the 50-fold gap is a strong caution
against importing the 4 ng L⁻¹ figure into a food matrix.

### 3c. ★ **TABLE 3 — THE TEMPERATURE LADDER. THE PAPER'S KEY ARTEFACT.** `[M]`

**Caption, raster-verified verbatim:** *"Table 3. Concentrations of 2M3F and 2FM in raw soy sauce
A after heating."* — p. 4.
**Column headers, raster-verified: `Compound | 2M3F (μg/L) | 2FM (μg/L)`, each split into
`5 min` and `20 min`.** Row axis: `Temperature`.
**Footnotes, raster-verified verbatim:** *"ᵃHeating time. ᵇBefore heat treatment. ᶜValues
correspond to the mean ± SD (standard deviation) of four analyses. ᵈ˒ᵉValues in the same column
are significantly different from each other (p < 0.05). ᶠ˒ᵍ˒ʰ˒ⁱValues with different superscripts
indicate significance (p < 0.01)."*
**Conditions:** raw soy sauce A, 200 mL, glass cylinder, dry stirring bath (80, 95 °C) or
laboratory autoclave (120 °C), pH ≈ 4.7 native, **unbuffered**, cooled in water to 25 °C, **n = 4**.
**Units μg L⁻¹ throughout.**

| **temperature** | **2M3F (MFT), 5 min** | **2M3F, 20 min** | **2FM (FFT), 5 min** | **2FM, 20 min** |
|---|---|---|---|---|
| **–** *(before heat treatment)* | **0.3 ± 0.1ᶜ** *(single value spanning both time columns)* | ← | **0.6 ± 0.1** *(single value spanning both)* | ← |
| **80 °C** | **1.2 ± 0.1** | **2.4 ± 0.3** | **1.2 ± 0.1** | **1.4 ± 0.1** |
| **95 °C** | **6.6 ± 0.4** | **6.8 ± 0.3** | **5.0 ± 0.3** | **5.2 ± 0.1** |
| **120 °C** | **11.1ᶠ ± 0.4** | **18.3ᵍ ± 0.8** | **9.7ʰ ± 0.7** | **18.2ⁱ ± 0.6** |

**⚠️ Read the significance letters carefully — they are not what they look like.**
- The `d,e` pair marks the **5 min columns** of 2M3F and 2FM respectively as headers
  (`5 minᵃ˒ᵈ` and `5 minᵉ`), and footnote d,e says *"Values in the same column are significantly
  different from each other (p < 0.05)."* ⇒ **within each 5-min column, the four temperature
  levels differ at p < 0.05.**
- The `f,g,h,i` superscripts appear **only on the four 120 °C cells**, all four distinct, with
  footnote *"Values with different superscripts indicate significance (p < 0.01)."*
- **⚠️ [!] There is therefore NO printed pairwise test between 95 °C and 120 °C, and none between
  80 °C and 95 °C.** The only stated comparisons are "within a column, some pair differs" and
  "the four 120 °C cells differ from each other". **⇒ Every leg-wise fold-change in §6 is
  statistically untested by the authors.** The paper's prose claim — *"The 2M3F and 2FM
  concentrations significantly increased with temperature"* (p. 3) — is asserted, not tabulated.

**Author's own reading, verbatim (p. 3):**
> *"The 2M3F and 2FM concentrations significantly increased with temperature. Although their
> concentrations significantly increased with heating time at 120 °C, **the concentrations of
> these thiols did not increase with heating time at 80–95 °C.** Moreover, the concentration of
> 2FM in raw soy sauce was approximately two times greater than that of 2M3F, but as compared to
> 2FM, 2M3F was formed in higher amounts by heating. Notably, at 95 °C, almost 7 μg/L of 2M3F and
> 5 μg/L of 2FM were detected in the sample. These results imply that as compared to 2FM, 2M3F is
> formed more easily in FSS under the same heating condition."*

**⚠️ The clause I have bolded is the single most important sentence in the paper for this wave —
the authors themselves state the 80–95 °C points are time-independent, i.e. AT PLATEAU. See
§6b and §7d; it is the confound that governs every Eₐ here.**

**Audit of the author's claims `[D]`:** *"almost 7 μg/L of 2M3F and 5 μg/L of 2FM"* at 95 °C ✔
(6.6–6.8 and 5.0–5.2). *"2FM … approximately two times greater than that of 2M3F"* in raw soy
sauce ✔ (Table 1: 1.3 / 0.7 = **1.86**; Table 3 unheated: 0.6 / 0.3 = **2.00**). *"did not
increase with heating time at 80–95 °C"* — ⚠️ **true at 95 °C (× 1.03, × 1.04) but FALSE at
80 °C for 2M3F, which doubles (1.2 → 2.4, × 2.00).** The claim over-generalises; **the plateau
is at 95 °C only.** `[!]`

**★ THE TIME-RATIO TABLE — the diagnostic that governs everything below** `[D]`

| temperature | **2M3F, 20 min ÷ 5 min** | **2FM, 20 min ÷ 5 min** | reading |
|---|---|---|---|
| **80 °C** | **× 2.000** | × 1.167 | 2M3F still climbing steeply; 2FM nearly flat |
| ★ **95 °C** | ★ **× 1.030** | ★ **× 1.040** | ★ **BOTH AT PLATEAU — the reaction is essentially complete by 5 min** |
| **120 °C** | **× 1.649** | **× 1.876** | **both climbing hard, nowhere near plateau** |

**⚠️ [!] THIS IS CHEMICALLY ANOMALOUS AND MUST BE FLAGGED.** A system that is at plateau at 95 °C
by 5 min but is *still rising at 20 min* at 120 °C is not behaving like a single reaction with a
single precursor pool. Three readings, none excluded by the data:
- **(H1, most likely)** **A first, low-barrier precursor pool is exhausted by 5 min at 95 °C, and
  a second, higher-barrier pool only opens above 95 °C.** This is chemically natural in a food
  matrix: free ribose/pentose + free Cys first, then slower routes from bound sugar, protein-bound
  cysteine, thiamine, or from HMF/furfural formed *in situ*.
- **(H2)** **A come-up-time artefact.** If the 95 °C samples spend much longer approaching
  temperature than the 120 °C autoclave samples, the "5 min" 95 °C point already carries a long
  effective thermal history and the "20 min" point adds proportionally little. **The paper states
  no come-up times.** `[NEG]`
- **(H3)** **A vessel artefact.** 80 and 95 °C are in an open glass cylinder, 120 °C is in a
  sealed autoclave. **Volatile loss from the open cylinder over 20 min would cap the 20-min
  values at 80–95 °C and not at 120 °C — producing exactly this pattern.** ⚠️ **For a compound
  with the volatility of MFT, over 20 min, from 200 mL of stirred hot liquid in an open cylinder,
  this is not a small effect.**

**⇒ H3 in particular would inflate the 80→95 leg's apparent Eₐ and depress the 95→120 leg's,
which is precisely the curvature observed. THE CURVATURE CANNOT BE SEPARATED FROM THE VESSEL
CHANGE FROM THE PUBLISHED RECORD.** This is the governing caveat of §7d.

### 3d. **TABLE 4 — the precursor-addition experiment at 120 °C** `[M]`

**Caption, raster-verified verbatim:** *"Table 4. Changes of 2M3F and 2FM in FSS samples at
120 °C for 5 minᵃ."* — p. 4.
**Footnotes, raster-verified verbatim:** *"ᵃSoy sauce samples (180 mL) were heated in a laboratory
autoclave. For these heating experiments, cysteine (1.8 mM) and the corresponding carbohydrate
(47 mM of glucose or 7.7 mM of ribose) were added. ᵇ˒ᶜ˒ᵉ˒ᶠValues with different superscripts
indicate significance in raw or heat-treated soy sauce (p < 0.01). ᵈValues correspond to the mean
± SD (standard deviation) of four analyses."*

| | **B (raw), 2M3F** | **B (raw), 2FM** | **BH (heat-treated), 2M3F** | **BH, 2FM** |
|---|---|---|---|---|
| **before heating** | **0.7ᵇ ± 0.1ᵈ** | **1.2 ± 0.1** | **1.3ᵉ ± 0.2** | **2.4 ± 0.2** |
| **after heating, not added** | **13.5ᶜ ± 0.3** | **8.9 ± 0.3** | **5.6ᶠ ± 0.4** | **4.4 ± 0.1** |
| **after heating, + ribose + Cys** | **15.7 ± 1.1** | **12.1 ± 0.1** | **6.6 ± 0.2** | **6.8 ± 0.3** |
| **after heating, + glucose + Cys** | **12.4 ± 1.4** | **11.4 ± 0.1** | **4.3 ± 0.4** | **6.1 ± 0.5** |

**Author's own reading, verbatim (p. 3):**
> *"As shown in Table 4, the 2M3F concentration significantly increased by approximately four
> times in BH and by 19 times in B at 120 °C for 5 min. Furthermore, the 2M3F concentration
> **tended to increase** in B and BH with the addition of ribose and cysteine
> **(0.05 < p < 0.1)**, but the 2M3F concentration did not increase by the addition of glucose and
> cysteine. In contrast, 2M3F was not detected in a pH 5 phosphate buffer by the addition of
> ribose or glucose and cysteine at 120 °C for 5 min (data are not shown)."*

**Audit `[D]`:** *"approximately four times in BH"* → 5.6/1.3 = **× 4.31** ✔. *"19 times in B"*
→ 13.5/0.7 = **× 19.29** ✔.

**★ Full derived fold-change table `[D]` — the paper prints none of these:**

| effect | **B (raw)** | **BH (pasteurised)** |
|---|---|---|
| **heating alone, 2M3F** | **× 19.29** | **× 4.31** |
| **heating alone, 2FM** | **× 7.42** | **× 1.83** |
| **+ ribose + Cys, on 2M3F** | × 1.163 | × 1.179 |
| ★ **+ glucose + Cys, on 2M3F** | ★ **× 0.919 (a DECREASE)** | ★ **× 0.768 (a DECREASE)** |
| **+ ribose + Cys, on 2FM** | × 1.360 | × 1.545 |
| ★ **+ glucose + Cys, on 2FM** | ★ **× 1.281 (an INCREASE)** | ★ **× 1.386 (an INCREASE)** |

**⭐ THE CLEANEST STRUCTURAL RESULT IN THE PAPER, AND THE AUTHORS DO NOT STATE IT NUMERICALLY:**

> **Adding glucose + cysteine at 120 °C DECREASES 2M3F (× 0.92, × 0.77 — in both samples) while
> simultaneously INCREASING 2FM (× 1.28, × 1.39 — in both samples). Adding ribose + cysteine
> raises both (2M3F × 1.16–1.18, 2FM × 1.36–1.55).**

**⇒ The pentose and the hexose route to the two thiols are not merely different in efficiency —
they have opposite signs on the MFT branch.** Ribose feeds MFT; glucose competes it away while
still feeding FFT. **This is a directly testable branching statement for the sulfur lane, it
replicates across two independent soy-sauce samples with the same sign and similar magnitude, and
it does not depend on 2FM's absolute calibration (it is a within-sample ratio).**

**⚠️ FOUR CAVEATS on that result, all material:**
1. ★ **The two sugar arms are not molar-matched: glucose 47 mM vs ribose 7.7 mM, a 6.1-fold
   difference (§1c).** The comparison is confounded with concentration. Whether ribose would
   still beat glucose at equimolar loading is untested here. `[!]`
2. **The authors themselves report only `0.05 < p < 0.1` for the ribose effect** — *"tended to
   increase"*. **The ribose enhancement is NOT significant at p < 0.05.** The glucose *decrease*
   is not tested at all.
3. **Every "after heating" number already contains the endogenous soy-sauce reaction** (13.5 and
   5.6 μg L⁻¹ with nothing added). The additions perturb a system that is already producing MFT
   at high yield, so these are **increments on a large baseline**, not de novo yields.
4. **Sample B and BH give the same signs but different magnitudes**, which is reassuring, but
   they are the same starting soy sauce (BH = pasteurised B), **not independent replicates.**

**⭐ A second derived result the paper does not draw `[D]`: pasteurisation depletes the
precursor pool.** B (raw) makes **13.5** μg L⁻¹ on heating; BH (already pasteurised once) makes
only **5.6** — **× 0.41.** The paper's own explanation (p. 3–4): *"when raw soy sauce was heated
for pasteurization, the Maillard reaction caused the reduction in the amino acids and sugar
concentrations. Thus, a higher amount of 2M3F is formed in B as compared with that in BH."*
**⇒ A 59 % loss of MFT-forming capacity from a single prior pasteurisation. That is a quantified
precursor-depletion effect in a real matrix, and it is the paper's best evidence that the
formation is precursor-limited rather than rate-limited.** It also **directly parallels
`feng2022_extraction.md` §6c**, where precursor depletion was shown to dominate the apparent
temperature response.

### 3e. **Derived branching ratio, FFT : MFT (2FM : 2M3F), across every condition** `[D]`

| condition | **2FM / 2M3F** |
|---|---|
| raw FSS survey average (Table 1) | **1.857** |
| heat-treated FSS survey average (Table 1) | **1.538** |
| soy sauce A, unheated | **2.000** |
| A, 80 °C / 5 min | 1.000 |
| A, 80 °C / 20 min | **0.583** |
| A, 95 °C / 5 min | 0.758 |
| A, 95 °C / 20 min | 0.765 |
| A, 120 °C / 5 min | 0.874 |
| A, 120 °C / 20 min | **0.995** |
| B, before / after heating | **1.714 → 0.659** |
| BH, before / after heating | **1.846 → 0.786** |

**⭐ A clean, five-fold-replicated result: FFT : MFT ≈ 1.7–2.0 in unheated soy sauce and collapses
to ≈ 0.6–1.0 on heating, in every sample tested.** Heating inverts the branching in favour of
MFT. **⇒ The FFT : MFT ratio is NOT a constant of the chemistry; it is a marker of thermal
history.** `[D]`

**⚠️ This is in tension with — but not contradicted by — the corpus's other two measurements.**
`kang2026_SI_extraction.md` §5c found FFT : MFT **invariant to pH** at ≈ 2.9–3.0 in the
xylose–Cys/TTCA system; `feng2022_extraction.md` §6d item 12 found it **≈ 1.15 and nearly
invariant to added nucleophile** in the ribose–GSH system. **Meng now adds a third system
(fermented soy sauce) with a third value (≈ 0.6–2.0) that is NOT invariant — it moves with
temperature and with thermal history.** Consolidating all three:

| system | FFT : MFT | invariant to | varies with |
|---|---|---|---|
| xylose–Cys (TTCA), Kang | **≈ 2.9–3.0** | **pH** (±5 %) | **nitrogen co-substrate** (0.15–9.2) |
| ribose–GSH (ARP), Feng | **≈ 1.15–1.45** | **added nucleophile** (0.2 %) | temperature (+26 % over 20 K) |
| **fermented soy sauce, Meng** | ★ **0.58–2.00** | — | ★ **temperature and thermal history (× 0.3 swing)** |

**⇒ The branching ratio is set by the SUGAR/PRECURSOR IDENTITY (a ≈ 2.6-fold spread between the
two model systems) and, in a real matrix where multiple precursor pools coexist, by which pool is
being consumed. It is not a universal constant. STRUCTURAL, three-paper triangulated, and a
genuine addition to the corpus.**

---

## §4. FIGURES

### ⚠️ **THERE ARE NO FIGURES IN THIS PAPER.** `[NEG]`

**Verified exhaustively:** the 6-page PDF (1 T&F cover page + 5 article pages) contains **four
tables and no figures.** There is **no "Figure 1"**, no scheme, no chromatogram, no mass spectrum,
no aroma-extract-dilution plot, no time course, no calibration plot, and **no figure caption of
any kind** in the text layer or on any raster render.

**Consequently: nothing was digitised, at any dpi, and no pixel-to-value calibration exists or is
needed.** There are no `[D]` figure values anywhere in this dossier — every `[D]` here is a
derivation from the printed tables, not a digitisation.

**⚠️ What the absence of figures costs, specifically:**
- **No chromatogram** — so the ET2MP / internal-standard m/z 134 co-elution risk (§2a) cannot be
  assessed, and the *"cooked meat-like aroma at 20.6 min"* identification cannot be seen.
- **No mass spectrum of the 2M3F peak in the sample** — identification rests on retention time +
  SIM ion set against a standard, described in prose only.
- **No time course beyond the two points (5 and 20 min)** — the plateau diagnosis of §3c rests on
  exactly two points per temperature.
- **No calibration plot** — the forced-through-origin `y = 0.5002x` cannot be inspected for
  curvature or leverage.
- **No aroma extract dilution analysis (AEDA) and no flavour dilution factors anywhere.** The
  wave brief anticipated *"every aroma extract dilution or OAV value"* — **there are none; the
  paper's only aroma-strength data are the cited thresholds of Table 1 and the paired-comparison
  counts of Table 2.** `[NEG]` The OAVs in §3a are `[D]`, computed by this wave.

---

## §5. PUBLISHED KINETICS `[F]` — **THERE ARE NONE. `[NEG]`**

**⚠️ Meng contains no rate constant, no reaction order, no activation energy, no half-life, no
rate law and no kinetic model.** Nothing to audit, nothing to refit. Every kinetic quantity in
this dossier is `[D]`.

### 5a. The paper's own quantitative claims, audited

| author claim | printed | **my audit `[D]`** | verdict |
|---|---|---|---|
| *"increased by approximately four times in BH"* | × 4 | 5.6 / 1.3 = **× 4.31** | ✔ |
| *"and by 19 times in B"* | × 19 | 13.5 / 0.7 = **× 19.29** | ✔ |
| *"at 95 °C, almost 7 μg/L of 2M3F and 5 μg/L of 2FM"* | ≈ 7, ≈ 5 | 6.6–6.8, 5.0–5.2 | ✔ |
| *"2FM … approximately two times greater than that of 2M3F"* (raw) | × 2 | 1.3/0.7 = **1.86**; 0.6/0.3 = **2.00** | ✔ |
| *"greater than 10 μg/L of 2M3F was formed in FSS"* | > 10 | B: **13.5** ✔; ⚠️ **BH: 5.6 — the claim holds for B only** | ✔ **partly** `[!]` |
| *"tended to increase … (0.05 < p < 0.1)"* (ribose effect) | 0.05 < p < 0.1 | not recomputable (per-replicate data absent) | untestable |
| *"the concentrations of these thiols did not increase with heating time at 80–95 °C"* | — | ⚠️ **95 °C: × 1.03, × 1.04 ✔. 80 °C: 2M3F × 2.00 ✗** | ✗ **over-generalised** `[!]` |

**⇒ Six of seven claims audit correctly; one over-generalises the plateau from 95 °C to
"80–95 °C". The paper's internal arithmetic is otherwise sound.**

### 5b. ★ **The two verified negatives the paper DOES establish — and what they are worth** `[M]`

**Negative 1 — the phosphate-buffer null.**
> *"2M3F was not detected in a pH 5 phosphate buffer by the addition of ribose or glucose and
> cysteine at 120 °C for 5 min (data are not shown)."* (p. 3)
> *"we believe that 2M3F is not formed in the phosphate buffer because of the low temperature and
> a short heating time."* (p. 4)

**Conditions, fully specified:** ribose **7.7 mM** *or* glucose **47 mM**, **plus cysteine
1.8 mM**, in **0.5 mol L⁻¹ phosphate buffer, pH 5.0, 180 mL, 120 °C, 5 min, laboratory
autoclave.**

**⇒ This is a genuinely useful bounded negative for the repo — the only one in wave K6a.** It
says the MFT formation rate from a clean, buffered ribose + Cys system at 120 °C is small enough
that 5 minutes produces less than the method's detection limit. **⚠️ AND THAT IS EXACTLY WHERE IT
BREAKS DOWN: the paper reports NO LOD (§2d).** Without an LOD the negative cannot be converted
into a number. The best that can be said, from the calibration range starting at 0 and the
smallest reported value being **0.3 ± 0.1 μg L⁻¹**, is:

> **`[D]` MFT yield from ribose (7.7 mM) + Cys (1.8 mM), 0.5 M phosphate pH 5.0, 120 °C, 5 min
> ≲ 0.3 μg L⁻¹ ≈ 2.6 nmol L⁻¹ — an ORDER-OF-MAGNITUDE upper bound, not a measured limit.**

Against a 1.8 mmol L⁻¹ cysteine loading that is a **conversion of ≲ 1.5 × 10⁻⁶ of the sulfur**.
**⇒ USE as a one-sided bound; REFUSE as a measurement.**

**Negative 2 — the ASS null.**
> *"2M3F was not detected in the aroma concentrates of three ASS samples analyzed by GC-O and
> GC-MS. Furthermore, 2M3F was not detected in any of the three ASS samples heated at 120 °C for
> 5 min."* (p. 4)

Acid-hydrolysed-vegetable-protein-mixed soy sauce, three commercial samples from the USA and
Brazil, unheated and heated. **⇒ The paper's central conclusion rests on this: fermentation is
required. But ASS differs from FSS in composition in *many* ways at once** (no microbial
metabolism, different free-amino-acid profile, different sugar profile, different pH, different
processing history), **and none of them is measured.** The paper's inference — *"the largest
difference between FSS and ASS or the phosphate buffer model is fermentation by
micro-organisms. Hence, fermentation by micro-organisms in FSS possibly plays a key role"* — is
**a difference-of-one-uncontrolled-variable argument across three uncontrolled variables.**
**PRIOR-ONLY.**

### 5c. Cited external kinetics `[C]` — not this paper's data, but relevant to the wave

| claim | source as cited | relevance |
|---|---|---|
| *"heating a phosphate buffer system by the addition of ribose and cysteine in the pH range of 3–7 at **145 °C for 20 min** results in the **decrease of 2M3F with increasing pH**. The greatest yields for glucose and cysteine were obtained at **pH 5**."* | **Hofmann & Schieberle 1998**, *J. Agric. Food Chem.* **46**, 235–241 (ref 12) | ★ **Independent confirmation of the acid-favours-thiols direction that `kang2026_SI_extraction.md` §5c measured (MFT × 0.57 from pH 5.5 → 8).** Two independent sources now agree on the sign. `[C]` |
| *"2M3F is detected from glucose and cysteine at **pH 5.5 and 150 °C** in a phosphate buffer, but 2M3F is **not detected at 120 °C**"* | **Ames, Guy & Kipping 2001**, *J. Agric. Food Chem.* **49**, 1885–1894 (ref 16, cited in the text as *"Jennifer et al."*) | ★ **An independent 120 °C null in buffer, corroborating Negative 1 above.** `[C]` ⚠️ **[!] The in-text attribution "Jennifer et al." uses the FIRST name of the first author (Jennifer M. Ames) — a citation-style error. The reference-list entry reads "Jennifer MA, Robin CEG, Gary JK", i.e. all three authors are given first-name-first. Do not propagate this attribution; the correct short form is Ames et al. 2001.** |
| *"2M3F … exhibits a cooked meat-like aroma at an extremely low odor threshold level of **4 ng/L** in a model dilute alcohol solution"* | Hofmann & Schieberle 1995 (ref 10); Tominaga & Dubourdieu 2006 (ref 11) | the Table 1 threshold; **wine-model matrix** `[C]` |
| *"the 2FM concentration has been reported to be increased by greater than two times by the Maillard reaction of cysteine and furfural"* | Meng, Kakuta & Sugawara 2012 (ref 8); Meng, Hatakeyama & Sugawara 2014 (ref 13) | mechanism for 2FM; **also the source of the 2FM quantification method (§2b)** `[C]` |
| *"2FM is generated by yeast during the fermentation of soy sauce"* | ref 13 | the fermentation hypothesis `[C]` |
| *"Sugar and cysteine have been reported as effective precursors of 2M3F"* | Mottram & Whitfield 1994 (ref 15) | `[C]` |

**⚠️ The Ames 2001 citation is load-bearing and pulls in the opposite direction from Meng's own
data.** Ames reports 2M3F **not detected** from glucose + Cys at 120 °C in buffer; Meng reports
**11.1–18.3 μg L⁻¹** from soy sauce at 120 °C. **The matrix, not the temperature, is doing the
work — and that is the paper's real finding.** It also means **⚠️ Meng's absolute MFT levels
cannot be transferred to a buffered model system at all** (§10).

---

## §6. DERIVED `[D]` — the three-point ladder

### 6a. **⚠️ THE CAVEAT CHAIN — read this before any number in §6b**

Two-point Eₐ from `Eₐ = R · ln(C₂/C₁) / (1/T₁ − 1/T₂)`. Three-point from an OLS fit of
`ln C` on `1/T`. *(Arithmetic validated: the same code reproduces Kang's published
120 → 140 °C legs, 97.8 kJ mol⁻¹ for MFT and 69.2 for FFT, as **97.86** and **69.04**.)*

1. **A two-point Eₐ has zero degrees of freedom** — no residual, **R² undefined and meaningless**,
   **no curvature detectable by construction**. This applies to Meng's individual legs.
2. ★ **A two-point Eₐ can CORROBORATE a slope magnitude. It can NEVER test for a slope BREAK.**
   **Meng, uniquely in this wave, has THREE heated temperatures and therefore ONE degree of
   freedom, so its 3-point fit CAN test for a break — weakly (n = 3, one residual pattern).**
   That single degree of freedom is the entire reason this paper matters to the wave.
3. **Yield-at-fixed-time is not a rate constant.** Both time columns are fixed holds; Table 3 is
   a yield table.
4. ★ **Precursor depletion depresses apparent Eₐ at the top of the ladder — and here it does the
   opposite of the usual thing.** The 95 °C system is demonstrably **saturated in time**
   (× 1.03 from 5 to 20 min), while the 120 °C system is **not** (× 1.65). **Saturation at the
   LOWER temperature inflates that temperature's yield relative to its rate, which DEPRESSES the
   apparent Eₐ of the leg ABOVE it and INFLATES the leg BELOW it.** ⇒ **The observed
   `122.8 → 25.0` collapse is exactly the artefact this mechanism produces.** See §7d.
5. **Semi-quant scale** — applies to 2FM only; cancels in every ratio (§2c).
6. **Unbuffered, and the pH is a food's pH.** Soy sauce at pH ≈ 4.7, unbuffered, no drift
   measurement, **and no per-sample pH at all.** Hofmann & Schieberle (ref 12, cited by Meng
   itself) report MFT yield **falling with increasing pH over 3–7**, so pH sensitivity is real
   and unmonitored here.
7. ★ **The vessel changes between the 95 and 120 °C rungs** (open glass cylinder → sealed
   autoclave), so the top rung of the ladder differs from the lower two in **pressure, headspace
   and volatile-retention regime**. **This is a first-order, unquantified confound on exactly the
   leg the Kang cross-check needs.** ⚠️ **No paper in this wave has a worse structural confound
   than this one.**
8. **Matrix.** Soy sauce is ≈ 16–18 % NaCl, 2–3 % ethanol, and its precursor pool is entirely
   unmeasured. Any Eₐ here is an **effective, whole-matrix** quantity, not an elementary barrier.

### 6b. **Two-point apparent Eₐ, leg by leg**

**Constants:** `1/353.15 − 1/368.15 = 1.15369 × 10⁻⁴ K⁻¹` (80→95 leg) ·
`1/368.15 − 1/393.15 = 1.72734 × 10⁻⁴ K⁻¹` (95→120 leg) ·
`1/353.15 − 1/393.15 = 2.88103 × 10⁻⁴ K⁻¹` (80→120, the full span).

| series | **80 → 95 °C** | | **95 → 120 °C** | | **80 → 120 °C (full span)** | |
|---|---|---|---|---|---|---|
| | fold | **Eₐ, kJ mol⁻¹** | fold | **Eₐ, kJ mol⁻¹** | fold | **Eₐ, kJ mol⁻¹** |
| ★ **MFT (2M3F), 5 min** | **× 5.500** | ★ **122.8** | **× 1.682** | ★ **25.0** | × 9.250 | **64.2** |
| **MFT (2M3F), 20 min** | × 2.833 | **75.0** | × 2.691 | **47.7** | × 7.625 | **58.6** |
| ★ **FFT (2FM), 5 min** | **× 4.167** | ★ **102.8** | **× 1.940** | ★ **31.9** | × 8.083 | **60.3** |
| **FFT (2FM), 20 min** | × 3.714 | **94.6** | × 3.500 | **60.3** | × 13.000 | **74.0** |

**Also computable, from the unheated control** (⚠️ **the unheated value is 0.3 ± 0.1, a 33 %
relative SD on a forced-through-origin curve — the least trustworthy cell in the table**):
unheated → 80 °C gives MFT × 4.0 (5 min) and × 8.0 (20 min); FFT × 2.0 and × 2.33. **No
temperature can be assigned to "unheated", so no Eₐ exists for this step. It bounds the
pre-existing pool, nothing more.**

### 6c. ★ **Three-point Arrhenius fits — the one thing only this paper can do**

OLS of `ln C` on `1/T` over the three heated temperatures (80, 95, 120 °C), **n = 3, 1 degree of
freedom**:

| series | **3-point Eₐ, kJ mol⁻¹** | **R²** | **residuals (80, 95, 120 °C)** | sign pattern |
|---|---|---|---|---|
| **MFT, 5 min** | **61.1** | ★ **0.839** | −0.321, **+0.536**, −0.214 | **− + −** |
| **MFT, 20 min** | **57.8** | 0.983 | −0.090, **+0.150**, −0.060 | **− + −** |
| **FFT, 5 min** | **58.1** | ★ **0.900** | −0.233, **+0.388**, −0.156 | **− + −** |
| **FFT, 20 min** | **72.9** | 0.984 | −0.112, **+0.188**, −0.075 | **− + −** |

**⭐ ALL FOUR SERIES SHOW THE IDENTICAL `− + −` RESIDUAL PATTERN.** A middle point lying *above* a
straight Arrhenius line with both ends *below* it means `ln C` vs `1/T` is **CONCAVE** — the local
slope magnitude falls as 1/T falls, i.e. **the apparent Eₐ DECREASES with increasing temperature.**
This is a **four-fold replication of the same-signed curvature**, across two analytes and two hold
times, and it is **the opposite sign** of the convex `+ − +` pattern a Kang-style switch-on would
produce.

**⚠️ How much weight this bears — stated honestly.**
- With **n = 3 and 1 residual degree of freedom, no formal curvature test is available.** R² of
  0.839 and 0.900 for the 5-min series are poor for a 3-point fit, which is itself informative;
  R² of 0.983 for the 20-min series is good, and there the curvature is small.
- **The four series are NOT independent** — they come from the same four heating runs on the same
  soy sauce A, so the replication is of the *pattern*, not of the *experiment*. **This is one
  experiment read four ways, not four experiments.**
- ★ **The 5-min series carry the strongest curvature AND the strongest saturation confound
  (95 °C at plateau by 5 min). The 20-min series carry weaker curvature and a weaker confound.
  The curvature scales with the confound.** That co-variation is exactly what an artefact would
  look like, and it is why §7d refuses to call this a clean contradiction of Kang.

**⇒ The honest statement: Meng shows a real, four-fold-consistent concave curvature whose sign
opposes Kang's, and whose magnitude tracks a known saturation/vessel artefact so closely that the
two cannot be separated from the published record.**

### 6d. Other derived quantities `[D]`

| # | quantity | value | basis | note |
|---|---|---|---|---|
| 1 | ★ **MFT 95 → 120 °C fold, 5 min** | ★ **× 1.682** (Eₐ 25.0) | Table 3 | **the closest analogue of Kang's low leg** |
| 2 | ★ **FFT 95 → 120 °C fold, 5 min** | ★ **× 1.940** (Eₐ 31.9) | Table 3 | ditto |
| 3 | **95 °C time-saturation index, MFT / FFT** | ★ **× 1.030 / × 1.040** | Table 3 | ★ **the governing confound** |
| 4 | **120 °C time index, MFT / FFT** | × 1.649 / × 1.876 | Table 3 | not saturated |
| 5 | **80 °C time index, MFT / FFT** | × 2.000 / × 1.167 | Table 3 | MFT not saturated |
| 6 | **pasteurisation precursor depletion** | ★ **× 0.41** (13.5 → 5.6 μg L⁻¹) | Table 4 | ★ **59 % capacity lost to one prior heat** |
| 7 | **glucose + Cys effect on MFT** | ★ **× 0.919 (B), × 0.768 (BH) — a DECREASE** | Table 4 | ★ **replicated sign** |
| 8 | **glucose + Cys effect on FFT** | ★ **× 1.281 (B), × 1.386 (BH) — an INCREASE** | Table 4 | ★ **opposite sign to MFT** |
| 9 | **ribose + Cys effect on MFT** | × 1.163 (B), × 1.179 (BH) | Table 4 | ⚠️ 0.05 < p < 0.1 |
| 10 | **ribose + Cys effect on FFT** | × 1.360 (B), × 1.545 (BH) | Table 4 | |
| 11 | **sugar : Cys molar ratio mismatch** | ★ **ribose 4.28 : 1, glucose 26.1 : 1** | §1c | ★ **the arms are not molar-matched** |
| 12 | **FFT : MFT branching, unheated → heated** | ★ **1.7–2.0 → 0.6–1.0** | §3e | ★ **heating inverts it, 5 samples** |
| 13 | **in-matrix JND for MFT** | ★ **0.2 μg L⁻¹ on a 0.5 μg L⁻¹ base, p < 0.01** | Table 2 | ★ **50 × the cited 4 ng L⁻¹ threshold** |
| 14 | **OAV, 2M3F in FSS** | 100–475 | §3a | ⚠️ wine-model threshold |
| 15 | **OAV, 2FM in FSS** | 2 000–7 000 | §3a | ⚠️ **2FM dominates 2M3F by ≈ 20 ×** |
| 16 | **lot-to-lot reproducibility floor** | ★ **21.6 %** (11.1 in A vs 13.5 in B at the same treatment) | §1a | ★ **exceeds every printed SD** |
| 17 | **MFT bound, buffer, 120 °C 5 min** | ★ **≲ 0.3 μg L⁻¹** | §5b | one-sided, no LOD |
| 18 | **sulfur conversion in the buffer null** | ★ **≲ 1.5 × 10⁻⁶ of added Cys** | §5b `[D]` | order of magnitude only |

---

## §7. ⭐ **THE KANG SWITCH-ON CROSS-CHECK — REQUIRED SECTION**

### 7a. The reference finding, restated

Kang et al. 2026 (TTCA 10 mmol L⁻¹ / Cys, pH 7.0, 120 min, sealed pressure vessel, oil bath,
HS-SPME-GC-MS, external calibration, n = 3 — `kang2026_SI_extraction.md` Table S4):

| | 100 °C | 120 °C | 140 °C | **low leg 100→120** | **high leg 120→140** |
|---|---|---|---|---|---|
| MFT, μg L⁻¹ | 1.237 | 1.388 | 5.907 | **× 1.12, Eₐ 6.9** | **× 4.26, Eₐ 97.9** |
| FFT, μg L⁻¹ | 3.734 | 4.107 | 11.439 | **× 1.10, Eₐ 5.8** | **× 2.78, Eₐ 69.0** |
| sulfur class | 13.978 | 35.866 | 60.400 | × 2.566, Eₐ 57.5 | × 1.684, Eₐ 35.2 |
| 2-acetylthiazole | 3.079 | 8.795 | 9.858 | × 2.856, Eₐ 64.0 | × 1.121, Eₐ 7.7 |

**Two separable claims: (K1) the free thiols SWITCH ON (Eₐ climbs 6.9 → 97.9); (K2) the sulfur
class as a whole DECELERATES (Eₐ falls 57.5 → 35.2).**

### 7b. **How well matched is Meng to Kang? — worse than Feng, and it must be said**

| design element | **Kang 2026** | **Meng** | comparable? |
|---|---|---|---|
| temperatures | 100 / 120 / 140 °C | ★ **80 / 95 / 120 °C — three levels** | ⚠️ **overlapping only at 120 °C; the spans abut rather than coincide** |
| **degrees of freedom for curvature** | **1** | ★ **1** | ✅ **the only two 3-point ladders in the corpus** |
| hold time | 120 min | ★ **5 and 20 min** | ⛔ **6–24 × shorter** |
| pH | **7.0**, NaOH-set, unbuffered | ★ **≈ 4.7, native, unbuffered** | ⛔ **2.3 pH units apart, and MFT is pH-sensitive** |
| matrix | clean aqueous model, 10 mmol L⁻¹ TTCA | ★ **fermented soy sauce: ≈ 17 % NaCl, 2–3 % ethanol, unmeasured precursor pool** | ⛔ **not comparable** |
| sulfur source | Cys bound in TTCA, 10 mmol L⁻¹ | **endogenous, unmeasured** (or +1.8 mmol L⁻¹ added) | ⛔ |
| vessel | sealed pressure-rated glass, all temperatures | ⚠️ **open glass cylinder at 80/95 °C, sealed autoclave at 120 °C** | ⛔ **changes mid-ladder** |
| extraction | HS-SPME | ★ **liquid–liquid DCM, from a 40 % dilution, concentrated to 25 μL** | ⚠️ different but internally constant |
| internal standard | 1,2-dichlorobenzene | 4-methoxy-2-methyl-2-mercaptobutane | different |
| calibration | external, per compound | **external for 2M3F only**; 2FM inherited from ref 8 | ⚠️ |
| replication | n = 3 | **n = 4** | ✅ |

**⇒ Meng is a WEAK replication of Kang's design and a STRONG complement to it.** It shares the
one thing Feng lacks — a third temperature — and shares almost nothing else. **Read it as an
independent probe of whether the *phenomenon* (a temperature-dependent apparent Eₐ for MFT/FFT)
exists in a completely different chemical setting, not as a repeat of Kang's measurement.**

### 7c. **THE FOLD-CHANGE TABLES**

**① The 95 → 120 °C leg — Meng's nearest analogue of Kang's 100 → 120 °C low leg**

| compound | **Meng, 5 min** | **Meng, 20 min** | **Kang, 100→120** | **Feng ARP-alone, 100→120** |
|---|---|---|---|---|
| **MFT** | **× 1.682** (Eₐ 25.0) | × 2.691 (Eₐ 47.7) | **× 1.12** (Eₐ 6.9) | **× 1.645** (Eₐ 30.3) |
| **FFT** | **× 1.940** (Eₐ 31.9) | × 3.500 (Eₐ 60.3) | **× 1.10** (Eₐ 5.8) | **× 2.071** (Eₐ 44.4) |
| **FFT ÷ MFT fold** | 1.153 | 1.301 | 0.982 | 1.259 |

**⭐ Meng's 5-min leg and Feng's ARP-alone leg agree to within 4 % on MFT (× 1.682 vs × 1.645) and
6 % on FFT (× 1.940 vs × 2.071) — across a real food matrix at pH 4.7 held 5 min, and a clean
aqueous model at pH 7.0 held 120 min.** That two such different systems land on the same
low-leg fold-change is a **strong, non-trivial corroboration** and the best cross-paper agreement
in the wave. **Kang's × 1.10–1.12 sits below both, by a factor of ≈ 1.5–1.8.**

**② The 80 → 95 °C leg — a lower leg no other paper in the corpus reaches**

| compound | **5 min** | **20 min** |
|---|---|---|
| **MFT** | **× 5.500** (Eₐ **122.8**) | × 2.833 (Eₐ 75.0) |
| **FFT** | **× 4.167** (Eₐ **102.8**) | × 3.714 (Eₐ 94.6) |

★ **This is the only sub-100 °C thiol-formation data in the entire corpus.** It says that below
95 °C the MFT/FFT channel carries an apparent barrier of **75–123 kJ mol⁻¹** — **much steeper**
than anything measured above 95 °C by anyone.

**③ The full three-point picture, all four series**

| series | 80→95 | 95→120 | **direction of Eₐ change** |
|---|---|---|---|
| MFT, 5 min | **122.8** | **25.0** | ★ **FALLS by 98 kJ mol⁻¹** |
| MFT, 20 min | 75.0 | 47.7 | **FALLS by 27** |
| FFT, 5 min | 102.8 | 31.9 | ★ **FALLS by 71** |
| FFT, 20 min | 94.6 | 60.3 | **FALLS by 34** |
| *(Kang, MFT)* | *6.9 (100→120)* | *97.9 (120→140)* | *★ **RISES by 91*** |
| *(Kang, FFT)* | *5.8* | *69.0* | *★ **RISES by 63*** |
| *(Kang, sulfur class)* | *57.5* | *35.2* | *FALLS by 22* |
| *(Kang, 2-acetylthiazole)* | *64.0* | *7.7* | *FALLS by 56* |

### 7d. **⚠️ THE VERDICT, IN THE REQUIRED VOCABULARY**

**On the 95 → 120 °C leg magnitude: ✅ CORROBORATED at short hold, ❌ CONTRADICTED at long hold.**
At 5 min, Meng gives MFT × 1.682 and FFT × 1.940 — the same order as Kang's × 1.12 / × 1.10, and
almost exactly Feng's × 1.645 / × 2.071. **Three independent systems place the MFT low-leg
fold-change in the band × 1.1 – × 2.1 (Eₐ 7 – 44 kJ mol⁻¹). That band is now triangulated and
should be treated as the corpus's low-leg prior.** ⚠️ At 20 min the same experiment gives
× 2.691 and × 3.500 (Eₐ 48 and 60). **⇒ THE "FLATNESS" IS HOLD-TIME-DEPENDENT, and Kang tested
only one hold time. That is a limitation of Kang's design that Meng exposes and that neither
Kang nor Feng could have revealed.**

**On the SLOPE BREAK — Meng CAN test it, unlike Feng: ⚠️ TESTED, AND THE SIGN OPPOSES KANG, BUT
THE TEST IS CONFOUNDED.**
**Meng's apparent Eₐ FALLS going up the ladder in all four series (122.8 → 25.0; 75.0 → 47.7;
102.8 → 31.9; 94.6 → 60.3), with a four-fold-consistent `− + −` residual pattern about a straight
Arrhenius line. Kang's RISES (6.9 → 97.9; 5.8 → 69.0). These are opposite-signed curvatures.**

**⚠️ But three confounds, each individually capable of producing the observed sign, and none
excludable from the published record:**

1. ★ **TIME SATURATION AT THE MIDDLE RUNG.** The 95 °C point gains only **× 1.03 (MFT) / × 1.04
   (FFT)** from 5 to 20 min — the authors say so themselves. The 120 °C point gains **× 1.65 /
   × 1.88**. **A middle rung sitting at its ceiling while the top rung is still climbing lifts the
   middle point above the Arrhenius line — producing exactly the `− + −` residual pattern
   observed, with no curvature in the underlying chemistry at all.** And the effect scales
   correctly: the 5-min series (worst saturation, R² 0.839/0.900) show the strongest curvature;
   the 20-min series (R² 0.983/0.984) show the weakest. ⚠️ **The curvature co-varies with the
   confound. That is the signature of an artefact, not of chemistry.**
2. ★ **THE VESSEL CHANGES BETWEEN THE 95 AND 120 °C RUNGS.** Open glass cylinder (dry bath) → sealed
   autoclave. **Volatile loss from an open 200 mL cylinder over 5–20 min at 80–95 °C would
   suppress the two lower rungs relative to the sealed top rung — which flattens the upper leg
   and steepens the lower one, again producing `− + −`.** ⚠️ **This confound is perfectly aligned
   with the leg boundary. It cannot be separated.**
3. **COME-UP TIME.** Never stated, and certainly unequal between a dry bath and an autoclave. At a
   nominal 5-minute hold this is a large fractional perturbation, and it acts differently on the
   two legs.

**⇒ THE HONEST VERDICT:**

> **Meng's three-point ladder is the only test of curvature available in wave K6a. It shows a
> clear, four-fold-consistent CONCAVE curvature — apparent Eₐ falling from ~100–123 to ~25–32
> kJ mol⁻¹ — which is the OPPOSITE SIGN to Kang's switch-on. Taken at face value it
> CONTRADICTS (K1).**
>
> **But the curvature is confounded with (i) a demonstrated time-plateau at the middle rung and
> (ii) a vessel/pressure change exactly at the leg boundary, either of which reproduces the
> observed pattern without any real change in barrier. ⚠️ THE CONTRADICTION MUST THEREFORE BE
> RECORDED AS UNRESOLVED, NOT AS A REFUTATION.**
>
> **⇒ Register as: STRUCTURAL — "the sign of the curvature disagrees between Kang and Meng, and
> the disagreement cannot be adjudicated from the published record of either." NEITHER PAPER
> SHOULD BE USED TO SET A SLOPE BREAK.**

**On (K2), class-level deceleration: ⛔ UNTESTABLE.** **Meng measures four thiols and no
heterocycles at all** — no thiazoles, no thiophenes, no pyrazines, no furans. **There is no
"sulfur class" here to accelerate or decelerate.** Meng bears on (K1) only.

**Two further things Meng contributes that neither Kang nor Feng can:**
- ★ **A sub-100 °C leg.** MFT's apparent Eₐ over 80 → 95 °C is **75–123 kJ mol⁻¹**, far steeper
  than anything above 95 °C in any paper. **If real, the corpus's picture becomes: steep below
  95 °C, flat 95–120 °C, and then (per Kang alone) steep again 120–140 °C — a non-monotone,
  two-break shape that NO SINGLE PAPER OBSERVES and that the corpus should treat as a hypothesis,
  not a finding.** ⚠️ **The 95–120 °C flat window is now supported by three papers (Kang, Feng,
  Meng-5-min); the steep legs on either side rest on ONE paper each. That asymmetry is the
  corpus's single most important open structural question.**
- ★ **A hold-time dimension.** Kang and Feng both fix the hold at 120 min. **Meng shows the
  apparent Eₐ of the same leg moving from 25.0 to 47.7 kJ mol⁻¹ (MFT) purely by changing the hold
  from 5 to 20 min.** ⚠️ **⇒ "The apparent Eₐ of MFT formation" is not a well-defined quantity
  without specifying the hold time. Any repo parameter of that name MUST carry its hold time, or
  it is meaningless.** This is the most directly actionable finding in the dossier.

---

## §8. VERIFIED NEGATIVES `[NEG]` — what a reader might hope is here, and is not

| # | what is missing | why someone would look for it |
|---|---|---|
| 1 | ★ **Any figure at all.** Four tables, zero figures. | No chromatogram, no spectrum, no time course, no calibration plot. §4. |
| 2 | ★ **Any aroma extract dilution analysis, flavour dilution factor, or printed OAV.** | The wave brief anticipated these. **They do not exist.** The OAVs in §3a are `[D]`. |
| 3 | ★ **Any LOD or LOQ.** | **The paper's two central results are NEGATIVES (buffer, ASS) and neither can be converted into a bound without an LOD.** §5b. The single most consequential omission. |
| 4 | ★ **Any rate constant, order, Eₐ, half-life or kinetic model.** | §5. |
| 5 | ★ **Any measurement of the precursor pool.** No sugar, no free amino acid, no cysteine, no thiamine, no furfural quantified in any soy sauce. | The matrix that does all the work is chemically unmeasured. |
| 6 | ★ **Any per-sample pH.** One approximate figure, "approximately 4.7", for FSS as a class. No pH for A, B, BH, or any ASS; no drift; no buffer in the soy-sauce arms. | |
| 7 | ★ **Any come-up time or cooling time.** | At a 5-min nominal hold this is first-order. §3c H2. |
| 8 | ★ **Any statement that the 80/95 °C glass cylinder was sealed.** | If open, volatile loss confounds the whole lower ladder. §3c H3. |
| 9 | ★ **Any temperature between 95 and 120 °C, and any above 120 °C.** | The gap where Kang's break is claimed to be (120–140 °C) is entirely unprobed. |
| 10 | **Any time point other than 5 and 20 min.** No t = 0 kinetics, no intermediate points, no plateau characterisation. | The saturation diagnosis rests on two points per temperature. |
| 11 | **Any calibration curve for 2FM, ET2MP or BM.** | §2b. 2FM — the FFT the Kang cross-check needs — is Tier B here. |
| 12 | **Any recovery, spike, matrix-matched calibration, blank or carry-over check.** | §2d. The 40 % soy-sauce/DCM extraction efficiency is unvalidated. |
| 13 | **Any retention index.** No RI, no Kovats, no n-alkane ladder. | Identification is retention-time-only on one column. |
| 14 | **Any chromatographic resolution evidence for the ET2MP / internal-standard m/z 134 isobar.** | §2a `[!]`. |
| 15 | **Any numeric data for the ASS samples.** Only "not detected". | The paper's fermentation conclusion rests on an unbounded null. |
| 16 | **Any numeric data for the phosphate-buffer model.** "(data are not shown)". | Same. §5b. |
| 17 | **Any threshold determined in soy sauce.** *"the perception threshold of 2M3F in soy sauce could not be directly determined"* (p. 3) — the paper says so explicitly. | All thresholds are wine-model `[C]` imports. |
| 18 | **Any measurement of H₂S, ammonia, furfural, HMF or 4-hydroxy-5-methyl-3(2H)-furanone.** | The mechanistic precursors of MFT and FFT are discussed and never measured. **No sulfur paper in the corpus measures H₂S.** |
| 19 | **Any browning, colour or absorbance measurement.** | Nothing at any wavelength. |
| 20 | **Any microbiological or enzymatic data**, despite the fermentation conclusion. The paper concedes: *"no studies have been reported on the microbial or enzymatic aspects associated with the production of 2M3F"* (p. 4). | The headline mechanism is entirely inferential. |
| 21 | **Any molar-matched sugar comparison.** Ribose 7.7 mM vs glucose 47 mM. §1c. | The ribose-vs-glucose conclusion is concentration-confounded. |
| 22 | **Any funding statement.** | §0. |
| 23 | **Any pairwise significance test between 80 and 95 °C, or between 95 and 120 °C.** | ⚠️ **Every leg-wise fold-change in §6 is statistically untested by the authors.** §3c. |
| 24 | **Any replicate presentation, order-balancing, or no-difference option in the sensory test; any rationale for the all-female panel.** | §1d. |
| 25 | **Any Supporting Information.** None exists; nothing is recoverable from a supplement. | Every gap above is permanent. |

---

## §9. CONSOLIDATED PARAMETER TABLE

| # | parameter | value | units | condition | provenance | source anchor |
|---|---|---|---|---|---|---|
| 1 | matrix | fermented soy sauce A, neat | — | 200 mL, glass cylinder / autoclave | `[M]` | p. 2 |
| 2 | **pH** | **≈ 4.7** | — | FSS, native, **unbuffered**, no drift measured | `[M]` | p. 4 |
| 3 | **temperatures** | **unheated, 80, 95, 120** | °C | ★ **three heated levels** | `[M]` | Table 3 |
| 4 | **times** | **5 and 20** | min | at every heated temperature | `[M]` | Table 3 |
| 5 | apparatus | dry stirring bath (80, 95 °C); **laboratory autoclave (120 °C)** | — | ⚠️ **vessel changes mid-ladder** | `[M]` | p. 2 |
| 6 | quench | water-cooled to 25 °C | — | no time given | `[M]` | p. 2 |
| 7 | come-up time | **not stated** | — | — | `[NEG]` | — |
| 8 | replication | **4** | — | Tables 3, 4 | `[M]` | footnotes |
| 9 | internal standard | 4-methoxy-2-methyl-2-mercaptobutane, **10.72** | μg L⁻¹ | added to 200 mL of 40 % soy sauce | `[M]` | p. 2 |
| 10 | extraction | 2 × 5 mL CH₂Cl₂, 1000 rpm, 5 min; + Na₂SO₃ 1 g + EtOAc 0.5 mL; conc. to 25 μL under N₂ | — | ⚠️ **liquid–liquid, not SPME** | `[M]` | p. 2 |
| 11 | column | DB-XLB, 60 m × 0.25 mm × 0.25 μm | — | GC-O and GC-MS, same programme | `[M]` | p. 2 |
| 12 | **2M3F calibration** | **y = 0.5002 x, r = 0.999, range 0–20** | μg L⁻¹ | forced through origin; x = peak-height ratio | `[F]` | p. 3 |
| 13 | 2FM / ET2MP / BM calibration | **not printed** | — | method inherited from ref 8 | `[NEG]` | — |
| 14 | quantifier ions | 114 (2M3F **and** 2FM), 134 (ET2MP), 124 (BM), 134 (IS) | m/z | ⚠️ IS/ET2MP isobar | `[M]` | p. 2 |
| 15 | ★ **MFT ladder, 5 min** | ★ **0.3 / 1.2 / 6.6 / 11.1** | μg L⁻¹ | unheated / 80 / 95 / 120 °C | `[M]` | Table 3 |
| 16 | ★ **MFT ladder, 20 min** | ★ **0.3 / 2.4 / 6.8 / 18.3** | μg L⁻¹ | ditto | `[M]` | Table 3 |
| 17 | ★ **FFT ladder, 5 min** | ★ **0.6 / 1.2 / 5.0 / 9.7** | μg L⁻¹ | ditto | `[M]` | Table 3 |
| 18 | ★ **FFT ladder, 20 min** | ★ **0.6 / 1.4 / 5.2 / 18.2** | μg L⁻¹ | ditto | `[M]` | Table 3 |
| 19 | ★ **MFT fold, 95 → 120 °C** | ★ **× 1.682** (5 min) / × 2.691 (20 min) | — | | `[D]` | §6b |
| 20 | ★ **FFT fold, 95 → 120 °C** | ★ **× 1.940** (5 min) / × 3.500 (20 min) | — | | `[D]` | §6b |
| 21 | **MFT fold, 80 → 95 °C** | × 5.500 (5 min) / × 2.833 (20 min) | — | | `[D]` | §6b |
| 22 | **FFT fold, 80 → 95 °C** | × 4.167 (5 min) / × 3.714 (20 min) | — | | `[D]` | §6b |
| 23 | ★ **MFT apparent Eₐ, 95 → 120 °C** | ★ **25.0** (5 min) / **47.7** (20 min) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 24 | ★ **FFT apparent Eₐ, 95 → 120 °C** | ★ **31.9** (5 min) / **60.3** (20 min) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 25 | ★ **MFT apparent Eₐ, 80 → 95 °C** | ★ **122.8** (5 min) / **75.0** (20 min) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 26 | ★ **FFT apparent Eₐ, 80 → 95 °C** | ★ **102.8** (5 min) / **94.6** (20 min) | kJ mol⁻¹ | 2-point | `[D]` | §6b |
| 27 | **MFT 3-point Eₐ** | **61.1** (5 min, R² 0.839) / **57.8** (20 min, R² 0.983) | kJ mol⁻¹ | OLS over 80/95/120 °C | `[D]` | §6c |
| 28 | **FFT 3-point Eₐ** | **58.1** (5 min, R² 0.900) / **72.9** (20 min, R² 0.984) | kJ mol⁻¹ | ditto | `[D]` | §6c |
| 29 | ★ **curvature sign** | ★ **CONCAVE (`− + −` residuals), all 4 series** | — | Eₐ falls with T | `[D]` | §6c |
| 30 | ★ **95 °C time-saturation index** | ★ **× 1.030 (MFT), × 1.040 (FFT)** | — | 20 min ÷ 5 min | `[D]` | §3c |
| 31 | 120 °C time index | × 1.649 (MFT), × 1.876 (FFT) | — | 20 min ÷ 5 min | `[D]` | §3c |
| 32 | 80 °C time index | × 2.000 (MFT), × 1.167 (FFT) | — | 20 min ÷ 5 min | `[D]` | §3c |
| 33 | MFT, raw FSS survey | 0.4–1.0, avg **0.7** | μg L⁻¹ | 5 commercial samples | `[M]` | Table 1 |
| 34 | MFT, heat-treated FSS survey | 0.5–1.9, avg **1.3** | μg L⁻¹ | 5 commercial samples | `[M]` | Table 1 |
| 35 | FFT, raw / HT FSS survey | 0.8–1.9 (1.3) / 1.2–2.8 (2.0) | μg L⁻¹ | | `[M]` | Table 1 |
| 36 | ET2MP, raw / HT | 5.4–27.4 (13.6) / 13.7–40.9 (22.8) | μg L⁻¹ | ⚠️ isobar risk | `[M]` | Table 1 |
| 37 | BM, raw / HT | 0.1–0.8 (0.3) / 0.3–0.7 (0.6) | μg L⁻¹ | | `[M]` | Table 1 |
| 38 | **MFT odour threshold** | **4** | ng L⁻¹ | ⚠️ **model dilute alcohol solution** | `[C]` refs 11, 17 | Table 1 |
| 39 | **FFT odour threshold** | **0.4** | ng L⁻¹ | ⚠️ same caveat | `[C]` | Table 1 |
| 40 | ET2MP / BM thresholds | 500 / 0.3 | ng L⁻¹ | ⚠️ same caveat | `[C]` | Table 1 |
| 41 | OAV, MFT in FSS | **100–475** | — | from ×38 | `[D]` | §3a |
| 42 | OAV, FFT in FSS | **2 000–7 000** | — | ★ **≈ 20 × MFT** | `[D]` | §3a |
| 43 | ★ **in-matrix JND, MFT** | ★ **0.2** on a 0.5 base, p < 0.01, N = 34 | μg L⁻¹ | trained panel, ISO 8586 | `[M]`/`[D]` | Table 2 |
| 44 | sensory panel | 34 female, 27–55 y, ISO 8586 | — | paired comparison, χ² | `[M]` | p. 2 |
| 45 | added ribose / glucose / Cys | **7.7 / 47 / 1.8** | mmol L⁻¹ | final conc., 180 mL, 120 °C 5 min | `[M]` | Table 4 fn a |
| 46 | ★ **sugar : Cys ratio** | ★ **4.28 : 1 (ribose) vs 26.1 : 1 (glucose)** | — | ⚠️ **arms not molar-matched** | `[D]` | §1c |
| 47 | ★ **glucose + Cys effect on MFT** | ★ **× 0.919 (B), × 0.768 (BH)** | — | **a DECREASE, both samples** | `[D]` | §3d |
| 48 | ★ **glucose + Cys effect on FFT** | ★ **× 1.281 (B), × 1.386 (BH)** | — | **an INCREASE, both samples** | `[D]` | §3d |
| 49 | ribose + Cys effect on MFT | × 1.163 (B), × 1.179 (BH) | — | ⚠️ 0.05 < p < 0.1 | `[D]`/`[M]` | §3d |
| 50 | ribose + Cys effect on FFT | × 1.360 (B), × 1.545 (BH) | — | | `[D]` | §3d |
| 51 | ★ **pasteurisation precursor depletion** | ★ **× 0.41** (13.5 → 5.6) | — | B vs BH, 120 °C 5 min | `[D]` | §3d |
| 52 | ★ **buffer null, MFT** | ★ **not detected; ≲ 0.3** | μg L⁻¹ | **ribose 7.7 mM or glucose 47 mM + Cys 1.8 mM, 0.5 M phosphate pH 5.0, 180 mL, 120 °C, 5 min** | `[M]`/`[D]` | p. 3 |
| 53 | ASS null, MFT | **not detected**, unheated and after 120 °C / 5 min | — | 3 commercial ASS, USA + Brazil | `[M]` | p. 4 |
| 54 | ★ **FFT : MFT branching** | ★ **1.7–2.0 unheated → 0.6–1.0 heated** | — | 5 sample sets | `[D]` | §3e |
| 55 | ★ **lot-to-lot reproducibility floor** | ★ **21.6 %** | — | A (11.1) vs B (13.5), same treatment | `[D]` | §1a |
| 56 | Hofmann & Schieberle: MFT falls with rising pH | pH 3–7, 145 °C, 20 min, ribose + Cys | — | corroborates Kang's pH direction | `[C]` ref 12 | p. 4 |
| 57 | Ames et al. 2001: MFT not detected at 120 °C in buffer, detected at 150 °C | glucose + Cys, pH 5.5 | — | corroborates the buffer null | `[C]` ref 16 | p. 4 |

---

## §10. USABILITY VERDICTS

| artefact | **verdict** | reason |
|---|---|---|
| ★ **Table 3 MFT ladder (80/95/120 °C × 5/20 min), fold-changes and within-study shape** | ★ **USE-Q** | Absolute external calibration for 2M3F (r = 0.999, in range), raster-verified end-to-end, n = 4, and the **only three-temperature MFT ladder in the corpus**. Qualified by: vessel change mid-ladder; 95 °C time plateau; no come-up times; no author significance test between rungs; real-food matrix. |
| **Table 3 FFT (2FM) fold-changes** | **USE-Q** | Fold-changes are valid under the wave rule even though 2FM's absolute scale is not calibrated in this paper (§2c). |
| **Table 3 FFT absolute μg L⁻¹ values** | **RATIO-ONLY** | No calibration curve printed here; the scale is inherited from ref. 8, which is not on disk. |
| **Table 3 MFT absolute μg L⁻¹ values** | **USE-Q** | Genuine external standard, in range. But extracted from a 40 % soy-sauce/DCM system with **no recovery and no matrix-matched calibration** — salting-out could shift the response factor (§2c item 3). Sound as a level *in this matrix by this method*; not portable. |
| **Absolute magnitude compared with Kang 2026 or Feng 2022** | ⛔ **REFUSE** | Three different matrices, extractions, internal standards and loadings. Only the within-paper temperature ratios are comparable. |
| **Two-point apparent Eₐ per leg (§6b)** | **PRIOR-ONLY** | Zero degrees of freedom per leg; yield-not-rate; unbuffered food matrix; hold-time dependent (25.0 vs 47.7 for the same leg). |
| ★ **The 95 → 120 °C leg at 5 min (MFT × 1.682, FFT × 1.940)** | ★ **USE-Q** | Agrees with Feng's ARP-alone leg to within 4 % / 6 %. **The strongest cross-paper agreement in the wave.** Contributes to the triangulated low-leg band × 1.1 – × 2.1. |
| ★ **Three-point Arrhenius fits and the concave curvature (§6c)** | ★ **STRUCTURAL** | Real, four-fold-consistent in sign, and the only curvature test available. **But confounded with time saturation at the middle rung and a vessel change at the leg boundary — the two cannot be separated. Record the disagreement with Kang as UNRESOLVED, never as a refutation.** |
| **Meng as evidence for or against Kang's 120 → 140 °C switch-on specifically** | ⛔ **REFUSE** | Meng's top rung is 120 °C. It has nothing above it. The break Kang claims is entirely outside Meng's range. |
| ★ **The hold-time dependence of apparent Eₐ (25.0 → 47.7 kJ mol⁻¹, same leg, 5 vs 20 min)** | ★ **USE** | Directly measured, internally controlled, no cross-paper transfer required. **⇒ Any repo parameter named "Eₐ of MFT formation" MUST carry its hold time.** The most actionable finding here. |
| ★ **The 80 → 95 °C leg (Eₐ 75–123 kJ mol⁻¹)** | **STRUCTURAL** | The corpus's only sub-100 °C thiol data. Carries the strongest vessel confound (open cylinder) and is unreplicated. **Direction only: steeper below 95 °C than above.** |
| ★ **Table 4 glucose-decreases-MFT / glucose-increases-FFT branching** | ★ **STRUCTURAL** | Same sign in two samples, within-sample ratio so calibration-independent. **⚠️ Confounded by the 6.1-fold molar mismatch between the sugar arms; direction only, never magnitude.** |
| **Table 4 ribose + Cys enhancement (× 1.16–1.18 on MFT)** | ⛔ **REFUSE** as a measurement; **PRIOR-ONLY** as a direction | The authors themselves report **0.05 < p < 0.1**. Not significant. |
| ★ **Pasteurisation precursor depletion, × 0.41** | ★ **STRUCTURAL** | Large, unambiguous, and mechanistically explained by the authors. Parallels `feng2022_extraction.md` §6c's extent effect in a real matrix. |
| ★ **The phosphate-buffer null at 120 °C / 5 min** | ★ **USE as a ONE-SIDED BOUND** (≲ 0.3 μg L⁻¹); ⛔ **REFUSE as a measurement** | Fully specified conditions, corroborated independently by Ames et al. 2001 `[C]`. **But no LOD exists, so the bound is order-of-magnitude only.** |
| **The ASS null and the fermentation conclusion** | **PRIOR-ONLY** | A difference-of-one-variable argument across at least four uncontrolled variables, with no LOD and no compositional data. |
| ★ **The in-matrix JND: 0.2 μg L⁻¹ detectable on a 0.5 μg L⁻¹ base, p < 0.01, N = 34** | ★ **USE** | Measured **in soy sauce**, trained ISO 8586 panel, χ² reproduces. The only in-matrix sensory number in the wave. ⚠️ **A JND, not an absolute threshold — the two must never be conflated.** |
| **Table 2 dose–response above +0.4 μg L⁻¹** | ⛔ **REFUSE** | Non-monotone (32 then 30 of 34); the test saturates. |
| **The 4 ng L⁻¹ MFT and 0.4 ng L⁻¹ FFT thresholds** | **PRIOR-ONLY** `[C]` | Wine-model values, cited not measured; the paper states the soy-sauce threshold *"could not be directly determined"*. **⚠️ The paper's own JND implies 200 ng L⁻¹ in this matrix — 50 × larger.** |
| **The `[D]` OAVs of §3a** | **PRIOR-ONLY** | Inherit the wine-model thresholds entirely. |
| **Table 1 raw-vs-heat-treated increase** | ⛔ **REFUSE** as a heating effect | Different commercial products, not paired samples. §3a `[!]`. |
| **ET2MP anything** | ⛔ **REFUSE** | Shares quantifier ion m/z 134 with the internal standard, with no resolution evidence. §2a `[!]`. |
| **Any LOD, LOQ, recovery, RI, browning, H₂S, precursor-pool or microbiological parameter** | ⛔ **REFUSE** | Not in this document, and there is no SI. §8. |
| **Pooling Tables 3 and 4** | ⛔ **REFUSE** | Different soy-sauce lots; 21.6 % discrepancy at the shared condition, exceeding every printed SD. §1a. |

---

## §11. DECLARED GAPS

### 11a. ★ **There is NO Supporting Information — every gap below is permanent**

Unlike `feng2022_extraction.md` (whose §11a lists five recoverable SI artefacts) and
`kang2026_SI_extraction.md`, **this paper has no supplement of any kind. Nothing in §8 can be
closed by ordering a file.** The only recoverable external artefact is:

- **Meng, Kakuta & Sugawara 2012**, *Food Sci. Technol. Res.* **18**, 429–436 (ref 8) — **the
  source of the 2FM / ET2MP / BM quantification method and, presumably, their calibration
  curves.** ★ **Ordering this would upgrade 2FM from Tier B to Tier A in this paper and give the
  FFT ladder an absolute scale.** It is the single highest-value follow-up for wave K6a, and it
  is the *only* one.

### 11b. What would actually be needed to resolve the wave's central question

| # | what is needed | why | obtainable? |
|---|---|---|---|
| 1 | ★ **An LOD for the 2M3F method** | Converts the buffer null and the ASS null from unbounded to bounded. | ⛔ not in the record |
| 2 | ★ **A statement of whether the 80/95 °C glass cylinder was sealed** | Determines whether the concave curvature of §6c is chemistry or volatile loss. **This single sentence would decide the paper's contribution to the Kang question.** | ⛔ |
| 3 | ★ **Come-up times for the dry bath and the autoclave** | At a 5-min hold, first-order. | ⛔ |
| 4 | ★ **A third hold time, or any intermediate time point** | Would separate the 95 °C plateau from a genuine barrier change. | ⛔ requires re-running |
| 5 | ★ **A temperature between 120 and 140 °C, in ANY system, with 3+ rungs and constant vessel** | ★ **The corpus's #1 gap. Kang's 140 °C column is the ONLY high-leg measurement anywhere, and it is a single unreplicated point.** | ⛔ **requires new experiment** |
| 6 | Precursor-pool quantification in soy sauce (sugars, free Cys, thiamine) | Would let the matrix effect be modelled rather than accepted. | ⛔ |
| 7 | Molar-matched ribose/glucose arms | Would de-confound the branching result of §3d. | ⛔ |

### 11c. What this dossier hands the next wave

- ★ **The corpus's only three-temperature MFT ladder, with one degree of freedom, raster-verified
  cell by cell.**
- ★ **A triangulated low-leg band: MFT × 1.1 – × 2.1 over ~100 → 120 °C, from three independent
  systems (Kang × 1.12, Feng × 1.645, Meng × 1.682 at 5 min).** Meng and Feng agree to 4 %.
- ★ **A curvature test whose sign OPPOSES Kang's switch-on, and an honest account of why it
  cannot be trusted to adjudicate.** ⇒ **Register as UNRESOLVED.**
- ★ **The finding that apparent Eₐ for the SAME leg moves 25.0 → 47.7 kJ mol⁻¹ purely with hold
  time.** ⇒ **No repo parameter called "Eₐ of MFT formation" is well-defined without a hold time.
  This should be checked against every sulfur-lane parameter currently in the registry.**
- ★ **The corpus's only sub-100 °C thiol data (80 → 95 °C, Eₐ 75–123 kJ mol⁻¹), suggesting a
  steep-flat-steep shape that NO SINGLE PAPER OBSERVES and whose two steep legs each rest on one
  paper.**
- ★ **A bounded negative (MFT ≲ 0.3 μg L⁻¹ from ribose + Cys in 0.5 M phosphate, pH 5.0, 120 °C,
  5 min), independently corroborated by Ames et al. 2001.**
- ★ **An in-matrix sensory JND of 0.2 μg L⁻¹ — 50 × the wine-model threshold the corpus currently
  carries for MFT.** ⚠️ **This should be checked against `k4b_paired_thresholds_and_browning.md`
  before any MFT threshold is used in a food-matrix prediction.**
- ★ **A three-system view of the FFT : MFT branching ratio (≈ 3.0 xylose–Cys, ≈ 1.15 ribose–GSH,
  0.6–2.0 soy sauce and moving with thermal history) showing it is set by precursor identity, not
  a universal constant.**
- ⚠️ **A citation flag: this PDF carries no volume, issue or page numbers and self-identifies as
  2016. Verify the final assignment against Crossref before freezing any bibliography.**
