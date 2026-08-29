# Ames, Guy & Kipping 2001 — pH and temperature in extruded cysteine/reducing sugar/starch

### Wave K6a per-paper extraction. 2026-08-29. **Nothing in `src/`, `tests/`, `results/`, `data/benchmarks/` or the declaration was touched.**

### ★ THE EXTRUSION LANE'S FIRST TEMPERATURE LADDER — AND THE FIRST MATRIX-PHASE TEST OF THE KANG SWITCH-ON.

---

## §0. IDENTITY — CONFIRMED, DOI PRESENT ON THE PAGE

| item | value | how verified |
|---|---|---|
| file on disk | `data/articles/ames2001.pdf` — **225,179 bytes** | `ls -la` |
| SHA-256 | `c4f9f16c762fafe5772bcf4c6bad46813c1ee555e5411109173d8653233316fc` | `shasum -a 256` |
| title | *"Effect of pH and Temperature on the Formation of Volatile Compounds in Cysteine/Reducing Sugar/Starch Mixtures during Extrusion Cooking"* | p.1885, raster |
| **authors** | ★ **Jennifer M. Ames**,\*,† **Robin C. E. Guy**,‡ **Gary J. Kipping**† | p.1885, raster |
| affiliation † | School of Food Biosciences, The University of Reading, Whiteknights, Reading **RG6 6AP**, United Kingdom | p.1885, raster |
| affiliation ‡ | Campden and Chorleywood Food Research Association, Chipping Campden, Gloucestershire **GL55 6LD**, United Kingdom | p.1885, raster |
| corresponding | J. M. Ames, tel +44 118 931 8730, fax +44 118 931 0080, `j.m.ames@afnovell.reading.ac.uk` | p.1885 footnote |
| **DOI** | ★ **`10.1021/jf0012547`** — ✔ **PRESENT**, printed in the p.1885 footer as `10.1021/jf0012547 CCC: $20.00`; article ID `JF0012547` repeated at the very end of p.1894 | raster-verified, both places |
| journal | ***J. Agric. Food Chem.* 2001, 49 (4), 1885–1894** (10 pages) | running heads on all 10 pages |
| dates | **Received for review October 18, 2000 · Revised manuscript received January 16, 2001 · Accepted January 16, 2001** · **Published on Web 03/07/2001** | p.1894 / p.1885 footer |
| **funding** | ★ **Ministry of Agriculture, Fisheries and Food (UK) and a consortium of food companies, Food Processing Sciences LINK Scheme, Project CSA 2172** | p.1894 |
| acknowledgment | Hazel Ling (extruder operation); **Donald Mottram and Harry Nursten** (advice) | p.1893 |
| version | published version of record, © 2001 American Chemical Society | p.1885 |
| producer | Acrobat Distiller Command 3.01 for Solaris (SPARC); creator Parlance Publisher 5.0/Xyvision PostScript Formatter; PDF 1.2, created 2001-04-05 | `pdfinfo` |

**⇒ Identity fully matches the brief. Revision was accepted the same day it was received (16 Jan 2001) — a same-day accept, worth noting only as a light signal on review depth; it is not a defect.**

---

## §0b. RASTER / UNIT CHECK — **THE ACS µ→m HAZARD DOES NOT BITE HERE**

Method: `pdftoppm -r 600 -png -f N -l N`, then read the PNG; plus `pdftotext -bbox` word coordinates for
column assignment. Pages rastered: **1, 2, 4, 5, 6** at 600 dpi; **7, 8** at 300 dpi.

| # | item | text layer says | raster says | verdict |
|---|---|---|---|---|
| 1 | internal standard mass/volume | `65 ng in 1 µL of diethyl ether` | **`65 ng in 1 µL of diethyl ether`** — the glyph is a true **µ** | ✔ **no µ→m corruption** |
| 2 | GC column film thickness | `0.5 µm film of BPX5` | **`0.5 µm film`** — true µ | ✔ agree |
| 3 | GC column dimensions | `50 m long × 0.32 mm i.d.` | `50 m long × 0.32 mm i.d.` | ✔ agree |
| 4 | ionization current | `50 µA` | **`50 µA`** — true µ | ✔ agree |
| 5 | Tenax mass | `6 mg of Tenax TA` | `6 mg` | ✔ agree |
| 6 | purge gas flow | `40 mL/min` | `40 mL/min` | ✔ agree |
| 7 | carrier flow | `1.6 mL/min` | `1.6 mL/min` | ✔ agree |
| 8 | feed moisture | `moisture content, 18%` | `moisture content, 18%` | ✔ agree |
| 9 | Table 1 header | `product temp (°C)`, `die pressure (MPa)`, `SME (kJ/kg)` | identical | ✔ agree |
| 10 | Table 4/5 unit footnote | `Amounts of compounds are quoted in ng/10 g of extrudate` | identical | ✔ agree — **ng per 10 g, not per g** |
| 11 | Table 4/5 replication footnote | `Figures quoted are the means of triplicate analyses` | identical | ✔ agree |
| 12 | column order inside each pH block | `180 °C, 150 °C, 120 °C` | identical — **DESCENDING** | ✔ agree ⚠ easy to invert |
| 13 | odour threshold units | `0.0032–0.0128 ng/L`, `0.0025 ng/L` | identical | ✔ agree |
| 14 | Table 1, all 18 rows × 4 numeric columns | see §1 | **all 72 cells identical** | ✔ agree |

**Numeric cell sample check.** Rather than sample, I verified **every** numeric cell of Tables 4 and 5 three
independent ways: (a) `pdftotext -layout`, (b) 600 dpi raster read, (c) `pdftotext -bbox` word coordinates with
each number assigned to a column by its right edge against the header `180/150/120 °C` centroids. **All three
agree on every cell and on every blank.** The bbox pass was necessary and it earned its keep — it resolved the
one genuinely ambiguous placement in the paper (3-acetylthiophene, §3.3).

⚠ **One raster-only finding, invisible in the text layer:** blanks in Tables 4 and 5 are **true absences**, not
zeros — and this is provable, because the printed class subtotals reproduce exactly when blanks are treated as
absent (see §3.4). Three cells are printed as an explicit `0` in total rows; those are real measured zeros.

---

## §1. SYSTEM AND CONDITIONS `[M]`

### 1.1 Feed formulation

| item | value | anchor |
|---|---|---|
| starch | **soft wheat starch (Abrastarch), ABR Foods, Corby, U.K.**; moisture **11.1 % w/w**; nitrogen **<0.15 % w/v** | Materials |
| batch basis | cysteine + sugar premixed with **1 kg of starch**, then blended into the remainder of the starch in a **ribbon blender** | Prep. of Extrudates |
| **cysteine** | **L-cysteine (98 %)**, Sigma — final **0.044 mol/kg** | Prep. of Extrudates |
| **sugar** | **D-glucose (99 %)** *or* **D-xylose (99.5 %)**, Sigma — final **0.044 mol/kg** | Prep. of Extrudates |
| molar ratio | ★ **1 : 1 cysteine : sugar**, both at 0.044 mol/kg of feed | derived from the above `[D]` |
| **pH adjustment** | ★ **sodium hydroxide added to the extruder WATER feed** — **34.2 / 23.4 / 19.7 g NaOH per litre of extruder water feed** for target feed pH **7.5 / 6.5 / 5.5** | Prep. of Extrudates |
| pH-5.5 point is not un-dosed | ⚠ note that **even the "pH 5.5" feed received 19.7 g/L NaOH** — there is no zero-alkali control; the ladder is three alkali doses, not alkali vs none | `[D]` from the above |
| moisture in | **18 %** (in-barrel target) | Prep. of Extrudates |
| **moisture out** | **`[NEG]` not reported** | — |

### 1.2 Extruder and processing

| item | value | anchor |
|---|---|---|
| extruder | **APV Baker (Peterborough, U.K.) MPF 50D co-rotating twin-screw** | Prep. of Extrudates |
| feed rate | **600 g/min** | Prep. of Extrudates |
| screw speed | **350 rpm** | Prep. of Extrudates |
| **modal residence time** | ★ **35 s** — quoted **once, for the whole study**, among "significant operating conditions" | Prep. of Extrudates |
| target die temperature | **120, 150 or 180 °C** | Prep. of Extrudates |
| **screw configuration** | **`[NEG]` not given** — deferred to refs *(18, 19)*: *"Full experimental details have been reported previously (18, 19)."* | Prep. of Extrudates |
| **barrel zone temperatures** | **`[NEG]` not given** — only the die target | — |
| **die dimensions** | **`[NEG]` not given** | — |
| **expansion ratio** | **`[NEG]` not quantified** — only *"to give expanded products"* | Prep. of Extrudates |
| **number of replicate EXTRUSIONS** | ★ **`[NEG]` — one extrusion per condition.** Replication is analytical only: *"**Triplicate isolates** were prepared for each set of extrusion conditions"* and Tables 4/5 report *"means of **triplicate analyses**"* | Isolation of Volatiles; Table 4/5 footnote *a* |
| storage before analysis | ground to powder, packed in polyamide–polyethylene laminate (Optivac), hermetically sealed, **−20 °C for ≤ 2 months** | Prep. of Extrudates |

★ **This is the single most important methodological fact in the paper for our purposes: `n = 3` is three
headspace isolates from ONE extrudate, not three extrusion runs. Every point on every ladder is one process
realisation. There is therefore no estimate of run-to-run process variance anywhere in the study, and §5 shows
that this is not academic — three of the six ladders contain a cell that is almost certainly a bad run.**

### 1.3 ★ TABLE 1 — MEASURED PROCESSING PARAMETERS `[M]` (p.1886, Table 1)

**The die temperatures are MEASURED, not target-only.** Table 1 gives a measured **product temperature** for
every one of the 18 conditions, alongside the measured **extrudate pH**, die pressure and SME.

| sugar, target die temp, and pH | extrudate pH | product temp (°C) | die pressure (MPa) | SME (kJ/kg) | **Δ vs target (°C)** `[D]` |
|---|--:|--:|--:|--:|--:|
| **glucose** | | | | | |
| 120 °C, pH 7.5 | 6.10 | **120** | 2.73 | 673 | 0 |
| 150 °C, pH 7.5 | 5.91 | **148** | 2.69 | 603 | −2 |
| 180 °C, pH 7.5 | 5.83 | **182** | 2.31 | 534 | +2 |
| 120 °C, pH 6.5 | 4.06 | **125** | 2.17 | 700 | +5 |
| 150 °C, pH 6.5 | 4.36 | **154** | 1.53 | 479 | +4 |
| 180 °C, pH 6.5 | 4.67 | **194** | 1.33 | 415 | **+14** |
| 120 °C, pH 5.5 | 3.55 | **117** | 2.29 | 578 | −3 |
| 150 °C, pH 5.5 | 3.26 | **146** | 1.09 | 425 | −4 |
| 180 °C, pH 5.5 | 3.28 | **173** | 0.97 | 371 | −7 |
| **xylose** | | | | | |
| 120 °C, pH 7.5 | 5.60 | **118** | 2.87 | 748 | −2 |
| 150 °C, pH 7.5 | 5.50 | **153** | 2.05 | 598 | +3 |
| 180 °C, pH 7.5 | 5.40 | **174** | 1.24 | 501 | −6 |
| 120 °C, pH 6.5 | 3.98 | **128** | 2.20 | 723 | +8 |
| 150 °C, pH 6.5 | 3.96 | **164** | 1.28 | 484 | **+14** |
| 180 °C, pH 6.5 | 3.90 | **182** | 1.20 | 462 | +2 |
| 120 °C, pH 5.5 | 3.10 | **115** | 2.30 | 696 | −5 |
| 150 °C, pH 5.5 | 3.09 | **150** | 0.50 | 417 | 0 |
| 180 °C, pH 5.5 | 3.13 | **177** | 0.54 | 393 | −3 |

Authors' own gloss, verbatim: *"Product temperature and pH differed from target values and varied between the
two sugars. **Product temperature ranged from −5 to +14 °C of the target value.** The pH of the product always
decreased compared to that of the feed; the pH drop was more pronounced at higher temperatures and was usually
greater for xylose."*

### 1.4 ⚠ THE LADDER'S X-AXIS IS **NOT** 120/150/180 — BUT THE DAMAGE IS MODERATE, NOT FATAL

**Flagging this prominently as instructed.** The nominal ladder is 120/150/180 °C. The achieved ladder differs
in every one of the 18 runs but two, and the deviations are **not symmetric across the ladder**:

| ladder | achieved (°C), 120→150→180 | achieved span | nominal span | span error |
|---|---|--:|--:|--:|
| glucose pH 7.5 | 120 → 148 → 182 | 62 | 60 | +3 % |
| glucose pH 6.5 | 125 → 154 → **194** | **69** | 60 | **+15 %** |
| glucose pH 5.5 | 117 → 146 → 173 | 56 | 60 | −7 % |
| xylose pH 7.5 | 118 → 153 → 174 | **56** | 60 | −7 % |
| xylose pH 6.5 | 128 → **164** → 182 | 54 | 60 | −10 % |
| xylose pH 5.5 | 115 → 150 → 177 | 62 | 60 | +3 % |

**Consequences, stated plainly:**

1. **Every Ea in §5 is reported on BOTH axes** — a target-temperature axis and a measured-product-temperature
   axis. Where they disagree materially I say so.
2. The distortion is **worst on the upper leg**. For xylose pH 7.5 — the ladder that carries the entire
   thiol result — the achieved upper leg is **153 → 174 °C = 21 °C**, not 30 °C. Fitting the same yield jump
   across a 30 % narrower temperature interval **inflates the apparent Ea on that leg by ~40 %**
   (146 → 208 kJ/mol for FFT). **Any Ea taken off the nominal axis for this ladder is too small; any taken off
   the measured axis is the larger and the more nearly correct of the two.**
3. **The three ladders are not on a common x-axis**, so the three pH ladders are not directly stackable.
4. ★ **SME and die pressure fall monotonically as die temperature rises, in all six ladders** (e.g. xylose
   pH 7.5: 748 → 598 → 501 kJ/kg; 2.87 → 2.05 → 1.24 MPa). **Mechanical energy input and shear are therefore
   ANTI-correlated with temperature across every ladder.** This is a confound with the same footprint as the
   temperature axis itself and it cannot be separated: the 180 °C runs are also the lowest-shear runs.

---

## §2. ANALYTICAL METHOD AND QUANTIFICATION BASIS — **SEMI-QUANTITATIVE, CONFIRMED**

### 2.1 Headspace isolation `[M]`

| item | value |
|---|---|
| sample | **10 g powdered extrudate + 20 mL distilled water**, 250 mL conical flask with Drechsel head |
| trap | glass-lined stainless steel tube, **6 mg Tenax TA** (SGE) |
| purge | **oxygen-free nitrogen, 40 mL/min, 1 h**, flask held at **37 °C** |
| dry-down | trap connected to N₂ alone for **5 min** to remove residual water |
| **internal standard** | ★ **1,2-dichlorobenzene — 65 ng in 1 µL diethyl ether**, injected onto the **front end of the Tenax tube under slight vacuum, immediately prior to analysis** |
| replication | **triplicate isolates** per condition; **blank isolations** with an empty flask |

⚠ **The IS is added to the trap after trapping, not to the sample before purging.** It therefore corrects for
desorption, transfer, GC injection and MS response drift — but **NOT** for headspace-partitioning or trapping
efficiency, which is exactly where analyte-to-analyte discrimination lives. This is a real limit on the
cross-compound comparability of the numbers and it is not acknowledged in the paper.

⚠ **The extrudate is re-wetted (10 g + 20 mL water) and purged at 37 °C for an hour.** The reported numbers are
therefore **not the extrudate's composition** but *what a re-hydrated extrudate releases into a nitrogen stream
at 37 °C over 1 h*. Any compound that hydrolyses, oxidises or reacts during that hour is measured after it.

### 2.2 GC-MS `[M]`

| item | value |
|---|---|
| instrument | **Hewlett-Packard 5972A MS + HP 5890 GC**, HP G1034C (C.01.05) workstation |
| desorption | modified Unijector (SGE), **250 °C for 5 min**, GC oven held at **0 °C** (subambient accessory), He **1.6 mL/min** |
| column | **50 m × 0.32 mm i.d., 0.5 µm BPX5** (SGE), ultralow-bleed |
| program | after desorption ramp to **60 °C over 1 min**; then **60 °C (5 min) → 250 °C at 4 °C/min → 250 °C (10 min)** |
| MS | **EI, 70 eV**; ion source **175 °C**; ionization current **50 µA**; **1.44 scans/s**; **m/z 29–300** |
| retention indices | experimental **LRI** vs **C6–C22 n-alkanes** run under the same conditions |

### 2.3 Identification basis `[M]`

Verbatim: *"Compounds were identified by comparing their mass spectra to those held in the mass spectrometer
data system library and other libraries (e.g., refs 21 and 22). … When both the MS and LRI data were consistent
with those in the literature or obtained for authentic compounds, identities were considered to be **positive**.
When only MS data were available, identities were considered to be **tentative**."*

⇒ **The `lit.` LRI column of Tables 4 and 5 is the identification-confidence flag.** A compound with a literature
LRI is *positive*; a compound with a blank `lit.` column is **tentative (MS-only)**. `[D]` **Counted: 28 of the
82 glucose rows and 8 of the 41 xylose rows are MS-only, i.e. tentative.**

★ **Both priority thiols are POSITIVE identifications** — MFT (873 / 870) and FFT (918 / 913) both carry a
literature LRI. **This matters: the compounds the wave cares about most are the ones the paper is most confident
about.**

⚠ **But two heavily-weighted glucose contributors are tentative:** `2(or 3)-thiophenethiol` — which is the
single largest thiophene in the glucose table (244 ng/10 g at 180 °C / pH 7.5, 28 % of the class) — and
`2-methylthiazole`, plus `4,5-dimethyloxazole`, which alone *is* the glucose oxazole class. The glucose
thiophene, thiazole and oxazole ladders therefore rest partly on MS-only assignments; note also that the
compound is written `2(or 3)-thiophenethiol`, i.e. **the authors could not place the substituent.**

### 2.4 ★ QUANTIFICATION BASIS — THE BOLDED SENTENCE

The authors' complete statement on quantification is one sentence:

> *"**Semiquantitative data were obtained from the mass spectral integration report, with reference to the
> internal standard.**"*

**⇒ QUANTIFICATION IS SEMI-QUANTITATIVE, NOT ABSOLUTE. There are no authentic-standard calibrations, no
measured response factors, and the word "response factor" does not appear anywhere in the paper. Every value in
Tables 4 and 5 is a total-ion-current peak area scaled to the 65 ng 1,2-dichlorobenzene internal standard and
then reported in the *units* of ng/10 g — i.e. an equal-response assumption is silently imposed on all 118
compounds.**

Corroborating verified negatives: `grep -ic "response factor"` → **0**; `"standard deviation"` → **0**;
`"calibration"` → **0**. No compound is stated to have been calibrated against an authentic standard for
*quantity* (authentic standards were used only for *LRI identity*).

### 2.5 ★ THE WAVE RULE, APPLIED EXPLICITLY

**A semi-quant ladder is STILL usable for activation energies, as within-study SHAPE.** For a fixed compound
measured on one instrument by one method, the unknown response factor `f` is a constant multiplier:
`Y_measured = f · Y_true`. An Arrhenius slope is taken on `ln Y` differences,

  `Ea = −R · Δ ln Y / Δ(1/T)`,  and  `ln(f·Y₂) − ln(f·Y₁) = ln(Y₂/Y₁)`,

so **`f` cancels identically**. The same holds for any ratio of the same compound across conditions. **Ea and
fold-change are therefore licensed by this dataset. This is the basis on which §5, §6 and §7 proceed.**

**What it does NOT license — enumerated, and each of these is a refusal:**

1. **Absolute yields in ppb / ng·g⁻¹.** The ng/10 g figures are area-ratio units wearing mass clothing. They
   must never enter the repo as concentrations, and never be compared to a measured ppb from any other paper.
2. **Cross-paper magnitude comparison.** `f` is method-specific; nothing here is commensurable with Kang's
   aqueous SIDA/absolute numbers or with anyone else's.
3. **Cross-compound comparison within this paper.** MFT vs FFT vs thiophene at the same condition differ by
   their unknown `f`s. ⚠ **This directly voids the paper's own headline "levels of FFT were about twice as
   great as those of MFT"** — that is an area ratio between two different compounds with two unknown response
   factors, not a molar statement. It reproduces from the tables (§3.5) but it is not a quantitative claim.
4. ★ **Comparing the glucose and xylose runs against each other.** The two sugars are reported in two separate
   tables with two separate grand totals and there is **no statement anywhere that they were analysed in one
   batch or sequence**. `[NEG]` — no run order, no batch statement, no QC sample, no bracketing standard. The
   paper nevertheless makes cross-sugar magnitude claims (*"The grand total amount of all compounds was often
   higher for xylose extrudates"*). **Those claims are outside what the method supports and this dossier does
   not carry them forward.** Within-sugar shape is what survives.
5. **Class subtotals as chemistry.** A "total thiophenes" of 876 is a sum of 22 unknown-`f` areas; its
   *changes* are usable, its *value* is not.

---

## §3. EVERY TABLE, CELL BY CELL `[M]`

### 3.1 Table inventory

| table | page | content | transcribed below |
|---|---|---|---|
| Table 1 | 1886 | Measured processing parameters, 18 conditions × 4 measures | §1.3 — **complete** |
| Table 2 | 1887 | Free-choice aroma descriptors, glucose extrudates | §3.6 — **complete** |
| Table 3 | 1887 | Free-choice aroma descriptors, xylose extrudates | §3.7 — **complete** |
| Table 4 | 1888–1889 | **80 compounds + 2 misc, glucose**, 9 conditions | §3.2 — **complete, all 82 rows** |
| Table 5 | 1890 | **38 compounds + 3 misc, xylose**, 9 conditions | §3.3 — **complete, all 41 rows** |

**Column key for Tables 4 and 5**, in the paper's own left-to-right order — note it runs **downward** in
temperature inside each pH block:

`7.5·180 | 7.5·150 | 7.5·120 || 6.5·180 | 6.5·150 | 6.5·120 || 5.5·180 | 5.5·150 | 5.5·120`

**Units, verbatim from footnote *a* of both tables:** *"Amounts of compounds are quoted in **ng/10 g of
extrudate**. Figures quoted are the means of triplicate analyses."*
Footnote *b*: *"Calculated linear retention indices for identified compounds."*
Footnote *c*: *"Linear retention indices obtained for authentic compounds analyzed on the same GC column or from
the literature (23, 24)."*

**`—` = blank in the original (compound not detected under that condition). `0` = an explicit printed zero,
which occurs only in total rows.** ⚠ **One fidelity note: the paper prints no `total oxazoles` row for xylose,
because that class has a single member; the row shown in §3.3 is marked `[D]` and is my own sum, identical to
the single compound.** The blank/zero distinction is not decorative: §3.4 proves the blanks are true
absences.

### 3.2 ★ TABLE 4 — CYSTEINE/GLUCOSE/STARCH EXTRUDATES `[M]` (p.1888–1889)

**All 82 rows. ng /10 g extrudate.**

| class / compound | LRI exptl / lit. | 7.5·180 | 7.5·150 | 7.5·120 | 6.5·180 | 6.5·150 | 6.5·120 | 5.5·180 | 5.5·150 | 5.5·120 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| **carbonyls** | | | | | | | | | | |
| &nbsp;&nbsp;2-butanone | 604 / 604 | 4 | 2 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-methylbutanal | 625 / 652 | 1 | 2 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-pentanone | 680 / 683 | 15 | 11 | 6 | 4 | 1 | 2 | 7 | — | — |
| &nbsp;&nbsp;pentanal | 690 / 689 | — | — | — | — | 1 | — | — | 2 | — |
| &nbsp;&nbsp;2,3-pentanedione | 709 / 702 | 6 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;3-hydroxy-2-butanone | 720 / 711 | — | 4 | 4 | 3 | — | — | — | — | — |
| &nbsp;&nbsp;*total carbonyls* | | **26** | **19** | **10** | **7** | **2** | **2** | **7** | **2** | **0** |
| **aliphatic S compds** | | | | | | | | | | |
| &nbsp;&nbsp;3-mercapto-2-butanone | 823 / 821 | 438 | 560 | 128 | 46 | 11 | 20 | 90 | 190 | 67 |
| &nbsp;&nbsp;3-mercapto-2-pentanone | 913 / 898 | 39 | 27 | 5 | 5 | — | 4 | — | — | — |
| &nbsp;&nbsp;2-mercapto-3-pentanone | 921 / 904 | 18 | 25 | 6 | — | — | — | — | — | — |
| &nbsp;&nbsp;bis(1-methyl-2-oxopropyl) disulfide | 1478 / 1476 | 1 | 2 | 2 | 2 | — | — | — | — | — |
| &nbsp;&nbsp;*total aliphatic S compds* | | **496** | **614** | **141** | **53** | **11** | **24** | **90** | **190** | **67** |
| **alicyclic S compds** | | | | | | | | | | |
| &nbsp;&nbsp;3-methyl-1,2-dithiolan-4-one | 1098 / 1102 | 39 | 23 | 9 | 21 | 1 | 9 | 15 | 1 | — |
| &nbsp;&nbsp;1,2,4-trithiolane | 1157 / — | 8 | 6 | 3 | — | — | — | — | — | — |
| &nbsp;&nbsp;3,5-dimethyl-1,2,4-trithiolane | 1172 / — | 13 | 3 | — | — | — | — | — | — | 1 |
| &nbsp;&nbsp;1,3,5-trithiane | 1271 / — | 4 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;3-methyl-1,2,4-trithiane | 1287 / — | 8 | 1 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total alicyclic S compds* | | **72** | **33** | **12** | **21** | **1** | **9** | **15** | **1** | **1** |
| **furans (non-S-containing)** | | | | | | | | | | |
| &nbsp;&nbsp;2-methylfuran | 608 / 606 | 1 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2,5-dimethylfuran | 701 / 707 | 2 | 2 | — | — | 1 | — | — | 2 | 1 |
| &nbsp;&nbsp;2-ethylfuran | 707 / 705 | 10 | 2 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-vinylfuran | 710 / 720 | 10 | 3 | 1 | — | — | — | — | — | — |
| &nbsp;&nbsp;2-furfural | 838 / 848 | — | — | — | 18 | 5 | 11 | 24 | 8 | 2 |
| &nbsp;&nbsp;2-furanmethanol | 863 / 867 | 10 | 2 | — | 3 | — | — | — | — | — |
| &nbsp;&nbsp;2-pentylfuran | 990 / 995 | 38 | 68 | 70 | 17 | 8 | 53 | 25 | 34 | — |
| &nbsp;&nbsp;*total furans (non-S-containing)* | | **71** | **77** | **70** | **39** | **14** | **64** | **49** | **44** | **3** |
| **pyrroles** | | | | | | | | | | |
| &nbsp;&nbsp;pyrrole | 744 / 748 | 212 | 2 | — | 7 | 2 | 4 | 34 | 4 | 13 |
| &nbsp;&nbsp;1-ethyl-1H-pyrrole | 818 / 817 | 18 | 13 | — | 5 | — | 6 | — | — | — |
| &nbsp;&nbsp;2-methyl-1H-pyrrole | 850 / 853 | 76 | 136 | 74 | 16 | 4 | 6 | — | — | — |
| &nbsp;&nbsp;2-ethyl-1H-pyrrole | 934 / 929 | 8 | 3 | 5 | — | — | — | — | — | — |
| &nbsp;&nbsp;3-ethyl-1H-pyrrole | 950 / — | 5 | 4 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;a trimethylpyrrole | 1026 / — | 10 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;1-(2-furanylmethyl)-1H-pyrrole | 1195 / 1196 | 4 | 1 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total pyrroles* | | **333** | **159** | **79** | **28** | **6** | **16** | **34** | **4** | **13** |
| **thiophenes** | | | | | | | | | | |
| &nbsp;&nbsp;thiophene | 665 / 672 | 72 | 26 | 5 | 21 | 8 | 9 | — | 2 | 1 |
| &nbsp;&nbsp;2,3-dihydrothiophene | 768 / — | 46 | 8 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-methylthiophene | 771 / 780 | 112 | 56 | 16 | 13 | 4 | 8 | 12 | 6 | 3 |
| &nbsp;&nbsp;3-methylthiophene | 783 / 784 | 124 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-ethylthiophene | 872 / 874 | 16 | 10 | 7 | 15 | 3 | 6 | 43 | — | — |
| &nbsp;&nbsp;2,5-dimethylthiophene | 877 / 878 | 20 | 25 | 13 | 13 | 3 | 13 | 10 | — | — |
| &nbsp;&nbsp;2,3-dimethylthiophene | 898 / 896 | 28 | 25 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-vinylthiophene | 915 / — | 9 | 4 | — | 4 | — | 3 | — | — | — |
| &nbsp;&nbsp;2,4-dimethylthiophene | 982 / — | 5 | 3 | — | 3 | — | — | — | — | — |
| &nbsp;&nbsp;2(or 3)-thiophenethiol | 988 / — | 244 | 196 | 22 | 77 | 25 | 69 | 34 | 5 | 12 |
| &nbsp;&nbsp;2-thiophenecarboxaldehyde | 1001 / 1000 | 32 | 21 | 5 | — | — | — | — | — | — |
| &nbsp;&nbsp;2-methyl-3-thiophenethiol | 1066 / — | 47 | 48 | 3 | 6 | — | 6 | — | — | — |
| &nbsp;&nbsp;2(or 3)-(methylthio)thiophene | 1104 / — | 15 | 5 | 4 | — | — | — | — | — | — |
| &nbsp;&nbsp;3-acetylthiophene | 1107 / — | 7 | 4 | 4 | — | — | — | — | — | 6 |
| &nbsp;&nbsp;3-methyl-2-thiophenecarboxaldehyde | 1109 / 1119 | 7 | 4 | 3 | — | — | — | — | — | — |
| &nbsp;&nbsp;2-acetylthiophene | 1112 / 1113 | 13 | — | — | — | — | — | 2 | — | — |
| &nbsp;&nbsp;5-methyl-2-thiophenecarboxaldehyde | 1141 / 1133 | 50 | 14 | 9 | 16 | — | 12 | 23 | 3 | — |
| &nbsp;&nbsp;2-pentylthiophene | 1165 / 1165 | — | 3 | 4 | — | — | — | — | — | — |
| &nbsp;&nbsp;2,3-dihydro-6-methylthieno[2,3c]furan | 1213 / — | 3 | — | — | 10 | — | — | — | — | — |
| &nbsp;&nbsp;thieno[2,3b]thiophene | 1239 / — | 20 | 4 | — | 46 | — | — | — | — | — |
| &nbsp;&nbsp;2-(2-thienyl)furan | 1244 / — | 3 | — | — | 8 | — | — | 3 | — | — |
| &nbsp;&nbsp;thieno[3,2b]thiophene | 1282 / 1248 | 3 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total thiophenes* | | **876** | **456** | **95** | **232** | **43** | **126** | **133** | **16** | **16** |
| **S-substituted furans** | | | | | | | | | | |
| &nbsp;&nbsp;2-methyl-3-furanthiol | 873 / 870 | 31 | 37 | 15 | 7 | 3 | 6 | — | 6 | — |
| &nbsp;&nbsp;2-furanmethanethiol | 918 / 913 | 106 | 77 | 18 | 18 | 4 | 18 | 32 | 16 | 12 |
| &nbsp;&nbsp;*total S-substituted furans* | | **137** | **114** | **33** | **25** | **7** | **24** | **32** | **22** | **12** |
| **pyridines** | | | | | | | | | | |
| &nbsp;&nbsp;pyridine | 741 / 750 | 7 | 6 | 4 | 1 | — | — | — | — | — |
| &nbsp;&nbsp;2-methylpyridine | 824 / 818 | 22 | 16 | 6 | 4 | — | — | — | 12 | — |
| &nbsp;&nbsp;*total pyridines* | | **29** | **22** | **10** | **5** | **0** | **0** | **0** | **12** | **0** |
| **pyrazines** | | | | | | | | | | |
| &nbsp;&nbsp;pyrazine | 736 / 743 | 82 | 17 | 4 | 5 | — | — | — | — | — |
| &nbsp;&nbsp;methylpyrazine | 830 / 826 | 274 | 138 | 25 | 18 | 4 | 8 | 12 | 12 | 5 |
| &nbsp;&nbsp;ethylpyrazine | 923 / 924 | 266 | 140 | 19 | 47 | 2 | 12 | 108 | 32 | 8 |
| &nbsp;&nbsp;2,3-dimethylpyrazine | 931 / 920 | 54 | 25 | 5 | — | — | — | — | — | — |
| &nbsp;&nbsp;vinylpyrazine | 946 / 948 | 31 | 5 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2-ethyl-5-methylpyrazine | 1006 / 1014 | 50 | 14 | 4 | 3 | — | — | — | — | — |
| &nbsp;&nbsp;2-ethyl-6-methylpyrazine | 1008 / 1016 | 54 | 21 | — | 2 | — | — | — | — | — |
| &nbsp;&nbsp;2,6-diethylpyrazine | 1089 / — | 8 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2,5-diethylpyrazine | 1091 / — | 4 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;2,3-dimethyl-5-ethylpyrazine | 1097 / — | 16 | 4 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;5-methyl-(5H)-6,7-dihydrocyclopentapyrazine | 1158 / — | 12 | 5 | 4 | 4 | — | 2 | — | — | — |
| &nbsp;&nbsp;3,5-diethyl-2-methylpyrazine | 1165 / — | 4 | 3 | — | — | — | — | — | 21 | — |
| &nbsp;&nbsp;*total pyrazines* | | **855** | **372** | **61** | **79** | **6** | **22** | **124** | **65** | **13** |
| **oxazoles** | | | | | | | | | | |
| &nbsp;&nbsp;4,5-dimethyloxazole | 750 / — | 154 | 223 | 74 | 26 | 3 | 6 | 23 | 1 | 1 |
| &nbsp;&nbsp;trimethyloxazole | 852 / 858 | 10 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total oxazoles* | | **164** | **223** | **74** | **26** | **3** | **6** | **23** | **1** | **1** |
| **thiazoles** | | | | | | | | | | |
| &nbsp;&nbsp;thiazole | 730 / 715 | 82 | 17 | 4 | 5 | 2 | 4 | 5 | — | — |
| &nbsp;&nbsp;2-methylthiazole | 815 / — | 42 | 51 | 15 | 16 | 8 | 18 | 23 | — | — |
| &nbsp;&nbsp;4-methylthiazole | 861 / — | 31 | 16 | 7 | 4 | — | — | — | — | — |
| &nbsp;&nbsp;2-methylthiazolidine | 893 / — | 8 | 5 | 4 | — | — | — | — | — | — |
| &nbsp;&nbsp;2-ethylthiazole | 906 / — | 14 | 9 | 4 | 4 | — | — | — | — | — |
| &nbsp;&nbsp;4,5-dimethylthiazole | 947 / 943 | 25 | 6 | 5 | 3 | 4 | 24 | 5 | 2 | — |
| &nbsp;&nbsp;5-ethylthiazole | 955 / 959 | 31 | 6 | 6 | 8 | 2 | 4 | 6 | — | — |
| &nbsp;&nbsp;trimethylthiazole | 1007 / 1005 | 27 | 11 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;5-ethyl-2-methylthiazole | 1022 / — | 84 | 15 | 5 | 15 | — | — | — | — | — |
| &nbsp;&nbsp;2-acetylthiazole | 1036 / 1027 | 40 | 9 | 5 | 20 | — | — | — | — | — |
| &nbsp;&nbsp;2,5-dimethyl-4-ethylthiazole | 1086 / — | 13 | 12 | 13 | — | — | — | 7 | — | — |
| &nbsp;&nbsp;*total thiazoles* | | **397** | **157** | **68** | **75** | **16** | **50** | **46** | **2** | **0** |
| **miscellaneous** | | | | | | | | | | |
| &nbsp;&nbsp;hexanenitrile | 884 / 884 | 17 | 62 | 22 | 9 | 7 | 23 | 4 | 5 | 8 |
| &nbsp;&nbsp;limonene | 1032 / 1033 | 6 | 11 | 7 | 7 | — | 6 | 7 | — | — |
| &nbsp;&nbsp;*total miscellaneous* | | **23** | **73** | **29** | **16** | **7** | **29** | **11** | **5** | **8** |
| **grand total** | | **3479** | **2319** | **682** | **606** | **116** | **372** | **564** | **364** | **134** |

### 3.3 ★ TABLE 5 — CYSTEINE/XYLOSE/STARCH EXTRUDATES `[M]` (p.1890)

**All 41 rows. ng /10 g extrudate.**

| class / compound | LRI exptl / lit. | 7.5·180 | 7.5·150 | 7.5·120 | 6.5·180 | 6.5·150 | 6.5·120 | 5.5·180 | 5.5·150 | 5.5·120 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| **carbonyls** | | | | | | | | | | |
| &nbsp;&nbsp;2-pentanone | 680 / 683 | 98 | 12 | 6 | 40 | 12 | 26 | 11 | 7 | 26 |
| &nbsp;&nbsp;2,3-pentanedione | 709 / 702 | 151 | 14 | 10 | 48 | — | 13 | — | — | 26 |
| &nbsp;&nbsp;*total carbonyls* | | **249** | **26** | **16** | **88** | **12** | **39** | **11** | **7** | **52** |
| **aliphatic S compds** | | | | | | | | | | |
| &nbsp;&nbsp;dimethyl disulfide | 744 / 756 | 13 | — | — | — | — | — | — | — | 10 |
| &nbsp;&nbsp;1-pentanethiol | 816 / 822 | — | — | — | — | — | 16 | — | — | — |
| &nbsp;&nbsp;3-mercapto-2-butanone | 823 / 821 | 721 | 43 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;3-mercapto-2-pentanone | 913 / 898 | 736 | 36 | 4 | 106 | — | 122 | — | — | — |
| &nbsp;&nbsp;*total aliphatic S compds* | | **1470** | **79** | **4** | **106** | **0** | **138** | **0** | **0** | **10** |
| **furans (non-S-containing)** | | | | | | | | | | |
| &nbsp;&nbsp;2,5-dimethylfuran | 701 / 707 | — | — | — | — | — | — | 17 | — | 23 |
| &nbsp;&nbsp;2-ethylfuran | 707 / 705 | 98 | 7 | 13 | 25 | 28 | 39 | 11 | 6 | 28 |
| &nbsp;&nbsp;2-furfural | 838 / 848 | 1670 | 66 | 21 | 1804 | 93 | 580 | 844 | 76 | 522 |
| &nbsp;&nbsp;2-pentylfuran | 990 / 995 | 517 | 28 | 116 | 377 | 49 | 438 | 246 | 54 | 305 |
| &nbsp;&nbsp;*total furans (non-S-containing)* | | **2285** | **101** | **150** | **2206** | **170** | **1057** | **1118** | **136** | **878** |
| **pyrroles** | | | | | | | | | | |
| &nbsp;&nbsp;2-methyl-1H-pyrrole | 850 / 853 | 380 | 26 | 13 | — | — | — | — | — | — |
| &nbsp;&nbsp;3-ethyl-1H-pyrrole | 950 / — | 35 | 2 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total pyrroles* | | **415** | **28** | **13** | **0** | **0** | **0** | **0** | **0** | **0** |
| **thiophenes** | | | | | | | | | | |
| &nbsp;&nbsp;thiophene | 665 / 672 | 178 | 10 | 4 | 67 | 38 | 19 | — | 2 | 20 |
| &nbsp;&nbsp;2-methylthiophene | 771 / 780 | 1487 | 106 | 62 | 280 | 65 | 264 | 17 | — | 135 |
| &nbsp;&nbsp;3-methylthiophene | 783 / 784 | 57 | — | — | 768 | 5 | — | — | 31 | — |
| &nbsp;&nbsp;2-ethylthiophene | 872 / 874 | 32 | 1 | — | 20 | — | — | — | — | — |
| &nbsp;&nbsp;2(or 3)-thiophenethiol | 988 / — | 103 | — | — | — | — | — | — | — | — |
| &nbsp;&nbsp;4,5-dihydro-2-methyl-3(2H)-thiophenone | 996 / 1000 | — | 1 | — | 71 | 10 | 71 | 53 | 9 | 54 |
| &nbsp;&nbsp;2-thiophenecarboxaldehyde | 1001 / 1000 | — | — | — | 69 | — | — | — | — | — |
| &nbsp;&nbsp;2-methyl-3-thiophenethiol | 1066 / — | 49 | 1 | 2 | 39 | — | 27 | — | — | 45 |
| &nbsp;&nbsp;2-pentylthiophene | 1165 / 1165 | 143 | — | 5 | 59 | 8 | 37 | 38 | 10 | 46 |
| &nbsp;&nbsp;*total thiophenes* | | **2049** | **119** | **73** | **1373** | **126** | **418** | **108** | **52** | **300** |
| **S-substituted furans** | | | | | | | | | | |
| &nbsp;&nbsp;2-methyl-3-furanthiol | 873 / 870 | 3038 | 171 | 44 | 604 | — | 764 | — | — | 400 |
| &nbsp;&nbsp;2-furanmethanethiol | 918 / 913 | 5346 | 341 | 66 | 1169 | 124 | — | — | — | 362 |
| &nbsp;&nbsp;2-methyl-3-[(2-furanylmethyl)thio]furan | 1503 / 1501 | — | — | 5 | — | 7 | 29 | — | — | 21 |
| &nbsp;&nbsp;bis(2-methyl-3-furanyl)disulfide | 1541 / 1547 | 20 | — | — | — | 4 | — | — | — | — |
| &nbsp;&nbsp;2-methyl-3-[(2-furanylmethyl)dithio]furan | 1650 / 1649 | — | — | 8 | — | 7 | 13 | — | — | — |
| &nbsp;&nbsp;*total S-substituted furans* | | **8404** | **512** | **123** | **1773** | **142** | **806** | **0** | **0** | **783** |
| **pyrazines** | | | | | | | | | | |
| &nbsp;&nbsp;pyrazine | 736 / 743 | 86 | 5 | 6 | 26 | — | — | — | — | — |
| &nbsp;&nbsp;methylpyrazine | 830 / 826 | 872 | 56 | 9 | 135 | — | 7 | — | — | — |
| &nbsp;&nbsp;ethylpyrazine | 923 / 924 | 708 | 4 | 7 | 80 | — | — | — | — | — |
| &nbsp;&nbsp;2-ethyl-6-methylpyrazine | 1008 / 1016 | 16 | 2 | 2 | 70 | 5 | — | — | — | — |
| &nbsp;&nbsp;propylpyrazine | 1010 / — | 316 | 9 | 21 | 66 | 4 | 133 | — | — | — |
| &nbsp;&nbsp;isopropylpyrazine | 1060 / — | 42 | — | 6 | — | — | — | — | — | — |
| &nbsp;&nbsp;isopropenylpyrazine | 1075 / 1092 | — | — | 9 | 100 | 9 | 50 | — | — | 47 |
| &nbsp;&nbsp;3-ethyl-2,5-dimethylpyrazine | 1081 / 1088 | 47 | — | — | 57 | 7 | 29 | 23 | 6 | 27 |
| &nbsp;&nbsp;*total pyrazines* | | **2087** | **76** | **60** | **534** | **25** | **219** | **23** | **6** | **74** |
| **oxazoles** | | | | | | | | | | |
| &nbsp;&nbsp;4,5-dimethyloxazole | 750 / — | 664 | 46 | 24 | 176 | — | — | — | — | — |
| &nbsp;&nbsp;*total oxazoles* `[D]` | | *664* | *46* | *24* | *176* | *0* | *0* | *0* | *0* | *0* |
| **thiazoles** | | | | | | | | | | |
| &nbsp;&nbsp;thiazole | 730 / 715 | 86 | 5 | 6 | 29 | — | — | — | — | — |
| &nbsp;&nbsp;2-methylthiazole | 815 / — | 195 | 14 | 29 | 135 | 23 | 161 | — | — | 82 |
| &nbsp;&nbsp;4-methylthiazole | 861 / — | 106 | 7 | — | — | — | — | — | — | — |
| &nbsp;&nbsp;*total thiazoles* | | **387** | **26** | **35** | **164** | **23** | **161** | **0** | **0** | **82** |
| **miscellaneous** | | | | | | | | | | |
| &nbsp;&nbsp;hexanal | 796 / 787 | 1028 | 90 | 75 | 826 | 49 | 406 | 604 | 74 | 534 |
| &nbsp;&nbsp;hexanenitrile | 884 / 884 | 180 | 12 | 8 | 87 | 10 | 112 | 47 | 15 | 79 |
| &nbsp;&nbsp;limonene | 1032 / 1033 | 42 | 2 | 14 | 235 | 29 | 246 | 54 | 17 | 256 |
| &nbsp;&nbsp;*total miscellaneous* | | **1250** | **104** | **97** | **1148** | **88** | **764** | **705** | **106** | **869** |
| **grand total** | | **19260** | **1117** | **595** | **7568** | **586** | **3602** | **1965** | **307** | **3048** |

### 3.4 ★ ARITHMETIC AUDIT OF THE PUBLISHED SUBTOTALS `[D]` — TABLE 5 IS CLEAN, TABLE 4 IS NOT

The paper prints class subtotals and a grand total for all 9 conditions in each table. I re-summed every one
from the transcribed compound rows, treating blanks as absent. **216 independent checks.**

| table | class-subtotal checks | grand-total checks | mismatches |
|---|--:|--:|--:|
| **Table 5 (xylose)** | 90 | 9 | ★ **0 — arithmetically flawless** |
| Table 4 (glucose) | 108 | 9 | **5** |

**The five glucose mismatches, in full:**

| # | class | condition | my sum | printed | Δ | reading |
|---|---|---|--:|--:|--:|---|
| 1 | furans (non-S) | 120 °C / pH 7.5 | 71 | **70** | +1 | rounding |
| 2 | furans (non-S) | 180 °C / pH 6.5 | 38 | **39** | −1 | rounding |
| 3 | **thiophenes** | **180 °C / pH 5.5** | 127 | **133** | **−6** | ★ see below |
| 4 | **thiophenes** | **120 °C / pH 5.5** | 22 | **16** | **+6** | ★ see below |
| 5 | pyrazines | 180 °C / pH 5.5 | 120 | **124** | −4 | unattributed |

**Mismatches 3 and 4 are one error, not two, and I can name it.** They are exactly complementary (−6 / +6)
within the same class and the same pH block. The only 6 in the glucose thiophene block at pH 5.5 is
**3-acetylthiophene**. I checked its placement by word bounding box: its `6` has right edge x ≈ 546 pt, against
column centroids of 473 pt (180 °C) and 539 pt (120 °C) — **it is typeset in the 120 °C column**, and the
600 dpi raster shows the same. But the authors' own subtotals only balance if that 6 is summed into the
**180 °C** column. **⇒ The published Table 4 contains a one-cell column-shift typesetting error:
3-acetylthiophene's 6 ng belongs at 180 °C / pH 5.5 and was set one block to the right.** Correcting it makes
both subtotals balance exactly (133 and 16).

Mismatches 1 and 2 are ±1 ng and are consistent with the compound rows being rounded to whole ng while the
totals were computed from unrounded areas. **This is a useful calibration: single-digit cells in these tables
carry ~±0.5 ng of rounding, i.e. up to ±50 % relative error at a printed value of 1.** Mismatch 5 (−4 ng in
pyrazines) I cannot attribute to any row and I record it as unresolved.

**Two consequences that matter more than the errors themselves:**

1. ★ **The blanks are true absences.** The subtotals reproduce to 211/216 *only* when blanks are excluded
   rather than read as zero. Where a compound is blank at one rung of a ladder, it was genuinely not detected —
   not merely unreported. This is what makes the "<3 points" refusals in §5 honest refusals.
2. ★ **The xylose table — which carries MFT, FFT and the entire thiol result — is arithmetically perfect
   across all 99 checks, while the glucose table is not.** The dataset the wave most wants to use is the more
   trustworthy of the two. That is a fortunate alignment and it is recorded as evidence, not as an assumption.

### 3.5 ★ INDEPENDENT VALIDATION AGAINST THE AUTHORS' OWN PROSE `[D]`

The paper states fourteen relative-abundance and total figures in the Results text. I recomputed each from my
transcription. **Thirteen of fourteen reproduce EXACTLY; the fourteenth reveals an error in the paper's text.**

| # | authors' statement | my value from the transcription | verdict |
|---|---|---|---|
| 1 | glucose S-compounds total range **78 → 1978** ng/10 g | 78 (150/6.5), 1978 (180/7.5) | ✔ exact |
| 2 | glucose S-compound RA **51 % → 72 %** | 51.2 % (120/7.5), 71.6 % (120/5.5) | ✔ exact |
| 3 | **44** S-containing compounds, glucose | 44 | ✔ exact |
| 4 | xylose S-compounds total range **52 → 12310** ng/10 g | 52 (150/5.5), 12310 (180/7.5) | ✔ exact |
| 5 | **21** S-containing compounds, xylose | 21 | ✔ exact |
| 6 | xylose S-compound RA min **5 %** (180 °C, pH 5.5) | 5.5 % at 180/5.5 | ✔ exact |
| 7 | xylose S-compound RA max **66 % (150 °C, pH 6.5)** | 65.9 % — but at **pH 7.5**, not 6.5 (pH 6.5 gives 49.7 %) | ⚠ **paper's pH label is wrong** |
| 8 | S-substituted furans up to **45.8 %**, xylose | 45.8 % (150/7.5) | ✔ exact |
| 9 | glucose thiophene RA **4.4–38.3 %** | 4.4 % (150/5.5), 38.3 % (180/6.5) | ✔ exact |
| 10 | xylose thiophene RA **5.5–21.5 %** | 5.5 % (180/5.5), 21.5 % (150/6.5) | ✔ exact |
| 11 | glucose aliphatic-S RA **6.5–52 %** | 6.5 % (120/6.5), 52.2 % (150/5.5) | ✔ exact |
| 12 | xylose aliphatic-S RA **0–7.6 %**, absent under **three** conditions | 7.6 % max; zero at 150/6.5, 180/5.5, 150/5.5 = **three** | ✔ exact |
| 13 | pyrazine RA **5.2–24.6 %** glucose, **1.2–10.8 %** xylose | 5.2/24.6 and 1.2/10.8 | ✔ exact |
| 14 | thiazole RA up to **13.8 %** glucose, **5.9 %** xylose; 2-furfural **43 %** of total at 180/5.5 xylose | 13.8 %, 5.9 %, 42.9 % | ✔ exact |

**⇒ My transcription is validated to the level of the authors' own derived statistics, and the compound counts
(80 in 11 classes + 3 miscellaneous, glucose; 38 in 9 classes, xylose; 44 and 21 sulfur-containing) reproduce
exactly. Finding #7 is a genuine erratum in the published text: the 66 % maximum sulfur RA occurs at pH 7.5,
not pH 6.5 as printed.**

The paper's FFT/MFT claim also reproduces: *"higher levels of FFT compared to MFT were obtained for all samples,
with the exception of the 120 °C, pH 6.5 and 5.5, xylose extrudates"* — confirmed (at xylose 120/6.5 FFT is
absent while MFT is 764; at xylose 120/5.5 FFT 362 < MFT 400; everywhere else FFT/MFT = 1.2–3.4). ⚠ But see
§2.5 item 3 — as a cross-compound area ratio this is **not** a quantitative claim.

### 3.6 TABLE 2 — DESCRIPTORS, CYSTEINE/GLUCOSE/STARCH EXTRUDATES `[M]` (p.1887)

Free-choice profiling, **five untrained but flavour-experienced assessors**, no intensity scale, no statistics.

| target product temp | pH 7.5 dry | pH 7.5 wet | pH 6.5 dry | pH 6.5 wet | pH 5.5 dry | pH 5.5 wet |
|---|---|---|---|---|---|---|
| **120 °C** | sulfur, egg, biscuity, perished rubber, cake-like, dough, warm, onion | puffed wheat, onion, cereal, toasted popcorn, sweet, strong, warm, musty, roasted, grilled corn | sulfur, egg, acidic, apples, sweet, burnt | yeasty/fruity, acrid, plasticine/varnish, sweet, strong, meaty, acid, toast, cereal | egg, sulfur, sweet, candy, cheesey | roasted meat, stewed fruit, sweet, plasticine, cheese, savory, acrid, bread, cereal, popcorn |
| **150 °C** | sulfur, egg, biscuity, sweet, toast, burnt, bread, puffed wheat, cornflakes, popcorn roasted | puffed wheat, onion, cereal, biscuity, bread, acrid, weak meaty | sulfur, egg, acid, apples | yeast/fruity, acrid, sweet, plasticine/varnish, sulfury, onion, burnt, acid, meaty, rubber, dentist drilling | egg, sulfur, sweet, roasted meat, rotten, savory, mildly acidic, apples, dentist drilling, acrid | roasted meat, stewed fruit, perished rubber, yeast, plasticine, sulfur, mildly burnt |
| **180 °C** | sulfur, egg, biscuity, meaty, fruity, puffed wheat, gravy browning, sweet, popcorn | puffed wheat, onion, meaty, old cooking oil, roasted/nutty, acrid, cabbagey, stale, sulfur, garlic, burnt, popcorn, bread, grilled corn | sulfur, egg, garlic, meat-like | yeast/fruity, acrid, sulfur, burnt, cereal, bread, stale, raw onion, roasted/nutty | egg, sulfur, sweet, acidic, apples, dentist drilling, acrid, savory | roasted meat, stewed fruit, cereal, sulfur, wet washing, stale, acidic, popcorn, bread, acrid, dentist drilling, warm, fresh yeast, paint/varnish |
| **common terms** | sulfur, egg, biscuity | puffed wheat, onion | sulfur, egg | yeast/fruity, acrid | egg, sulfur, sweet | roasted meat, stewed fruit |

### 3.7 TABLE 3 — DESCRIPTORS, CYSTEINE/XYLOSE/STARCH EXTRUDATES `[M]` (p.1887)

| target product temp | pH 7.5 dry | pH 7.5 wet | pH 6.5 dry | pH 6.5 wet | pH 5.5 dry | pH 5.5 wet |
|---|---|---|---|---|---|---|
| **120 °C** | puffed wheat, cereal, sulfur, egg, rubber, acidic | puffed wheat, cooked meat, slightly acidic, varnish, sweet, cereal, bread, apples | burnt, meaty, acid, sulfur, sweet, egg | puffed wheat, burnt varnish, sharp, sulfur, burnt, corn, bread, roasted | meat, burnt, weak, green apples, wet paint, acidic, roast, sweet, egg | puffed wheat, meat burnt, toast, weak acid, sulfur, old cooking oil, fruity, apples, corn, sweet, putty, roast, acrid |
| **150 °C** | puffed wheat, cereal, sulfur, biscuit, burnt, acrid, sharp, onion-like, meaty, sweet | puffed wheat, cooked meat, grilled corn, raw onion, burnt, bread, sharp, roast | burnt, meaty, acid, sulfur, garlic, cereal, acrid, paint | puffed wheat, burnt varnish, sharp, varnish, rubber, fruity/yeasty, roast, bread, sulfur, corn | meat, burnt, weak, sulfur, fatty | puffed wheat, meat, burnt toast, *wak* [sic — "weak"], sulfur, roasted, onion, bread |
| **180 °C** | puffed wheat, cereal, sulfur, egg, meaty, biscuit, rubber, slightly apples | puffed wheat, cooked meat, sharp, onion, burnt, acrid, stale, acid, metallic | burnt, meaty, acid, toast, acrid, garlic, wet paint, metallic, rubbery, egg | puffed wheat, burnt varnish, sharp, acrid, fruity/yeasty, cereal, onion, toast, stale | meat, burnt, weak, sulfur, biscuit, rubber, onion | puffed wheat, meat, burnt, toast, weak, onion, sharp, rubber |
| **common terms** | puffed wheat, cereal, sulfur | puffed wheat, cooked meat | burnt, meaty, acid | puffed wheat, burnt varnish, sharp | meat, burnt, weak | puffed wheat, meat, burnt, toast, weak |

Authors' summary of the sensory axis, verbatim: *"Meaty terms were used more frequently for samples processed at
higher temperature or lower pH or prepared using xylose rather than glucose… The intensity of the odor increased
with pH… Moistening samples reduced the use of the terms 'egg' and 'sulfur' and increased the frequency of the
terms 'onion', 'roast', 'toast', and 'burnt'."*

⚠ **This is free-choice profiling by an untrained panel with no replication, no intensity scaling and no
statistics. It is `STRUCTURAL` evidence only — direction of travel, never a magnitude.** Note in particular that
the sensory "meatiness" rises toward **lower** pH while MFT and FFT rise sharply toward **higher** pH — the
descriptors and the chemistry point in opposite directions on the pH axis, which is itself a caution against
reading the panel as an assay for the thiols.

---

## §4. EVERY FIGURE `[D]`

| figure | page | content |
|---|---|---|
| Figure 1 | 1891 | RAs of **sulfur-containing** classes, cysteine/glucose/starch |
| Figure 2 | 1891 | RAs of **non-sulfur-containing** classes, cysteine/glucose/starch |
| Figure 3 | 1891 | RAs of **sulfur-containing** classes, cysteine/xylose/starch |
| Figure 4 | 1891 | RAs of **non-sulfur-containing** classes, cysteine/xylose/starch |
| Figure 5 | 1892 | RAs **normalised to thiophenes** — aliphatic S, thiophenes, pyrazines — glucose |
| Figure 6 | 1892 | RAs **normalised to thiophenes** — non-S furans, thiophenes, S-containing furans — xylose |

### 4.1 ★ I DID NOT DIGITISE THESE FIGURES — AND DIGITISING THEM WOULD HAVE BEEN THE WRONG CHOICE

**Stated plainly, with the reasoning, as the brief requires.**

Rendered at 300 dpi (`pdftoppm -r 300`) and inspected. All six are **3-D perspective bar charts** with receding
depth axes, drawn with the bar tops floating off the plotted grid plane. Digitising a 3-D perspective bar to its
value axis requires reconstructing the projection, and the residual error on such charts is routinely 10–25 % —
far larger than any signal I would be extracting.

**But the decisive point is that digitisation is unnecessary, because these figures are pure restatements of
Tables 4 and 5.** Every quantity plotted is `class total ÷ grand total × 100`, and both inputs are printed. I
therefore **reconstructed the figures exactly by computation rather than approximately by digitisation**, and
§3.5 confirms the reconstruction against fourteen RA values the authors quote in the text — thirteen match to
the last decimal.

**⇒ The table below IS Figures 1–4, at full arithmetic precision, with zero digitisation error.** Figures 5 and
6 are the same numbers divided by the thiophene row and are not reproduced separately.

### 4.2 FIGURES 1–4 RECONSTRUCTED — relative abundance (% of grand total) `[D]`

**Glucose (Figures 1 and 2):**

| class | 7.5·180 | 7.5·150 | 7.5·120 | 6.5·180 | 6.5·150 | 6.5·120 | 5.5·180 | 5.5·150 | 5.5·120 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| carbonyls | 0.7 | 0.8 | 1.5 | 1.2 | 1.7 | 0.5 | 1.2 | 0.5 | 0.0 |
| aliphatic S | 14.3 | 26.5 | 20.7 | 8.7 | 9.5 | 6.5 | 16.0 | 52.2 | 50.0 |
| alicyclic S | 2.1 | 1.4 | 1.8 | 3.5 | 0.9 | 2.4 | 2.7 | 0.3 | 0.7 |
| furans (non-S) | 2.0 | 3.3 | 10.3 | 6.4 | 12.1 | 17.2 | 8.7 | 12.1 | 2.2 |
| pyrroles | 9.6 | 6.9 | 11.6 | 4.6 | 5.2 | 4.3 | 6.0 | 1.1 | 9.7 |
| thiophenes | 25.2 | 19.7 | 13.9 | 38.3 | 37.1 | 33.9 | 23.6 | 4.4 | 11.9 |
| S-substituted furans | 3.9 | 4.9 | 4.8 | 4.1 | 6.0 | 6.5 | 5.7 | 6.0 | 9.0 |
| pyridines | 0.8 | 0.9 | 1.5 | 0.8 | 0.0 | 0.0 | 0.0 | 3.3 | 0.0 |
| pyrazines | 24.6 | 16.0 | 8.9 | 13.0 | 5.2 | 5.9 | 22.0 | 17.9 | 9.7 |
| oxazoles | 4.7 | 9.6 | 10.9 | 4.3 | 2.6 | 1.6 | 4.1 | 0.3 | 0.7 |
| thiazoles | 11.4 | 6.8 | 10.0 | 12.4 | 13.8 | 13.4 | 8.2 | 0.5 | 0.0 |
| miscellaneous | 0.7 | 3.1 | 4.3 | 2.6 | 6.0 | 7.8 | 2.0 | 1.4 | 6.0 |
| **Σ sulfur classes** | **56.9** | **59.2** | **51.2** | **67.0** | **67.2** | **62.6** | **56.0** | **63.5** | **71.6** |
| *grand total, ng/10 g* | *3479* | *2319* | *682* | *606* | *116* | *372* | *564* | *364* | *134* |

**Xylose (Figures 3 and 4):**

| class | 7.5·180 | 7.5·150 | 7.5·120 | 6.5·180 | 6.5·150 | 6.5·120 | 5.5·180 | 5.5·150 | 5.5·120 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| carbonyls | 1.3 | 2.3 | 2.7 | 1.2 | 2.0 | 1.1 | 0.6 | 2.3 | 1.7 |
| aliphatic S | 7.6 | 7.1 | 0.7 | 1.4 | 0.0 | 3.8 | 0.0 | 0.0 | 0.3 |
| furans (non-S) | 11.9 | 9.0 | 25.2 | 29.1 | 29.0 | 29.3 | 56.9 | 44.3 | 28.8 |
| pyrroles | 2.2 | 2.5 | 2.2 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| thiophenes | 10.6 | 10.7 | 12.3 | 18.1 | 21.5 | 11.6 | 5.5 | 16.9 | 9.8 |
| **S-substituted furans** | **43.6** | **45.8** | **20.7** | **23.4** | **24.2** | **22.4** | **0.0** | **0.0** | **25.7** |
| pyrazines | 10.8 | 6.8 | 10.1 | 7.1 | 4.3 | 6.1 | 1.2 | 2.0 | 2.4 |
| oxazoles | 3.4 | 4.1 | 4.0 | 2.3 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| thiazoles | 2.0 | 2.3 | 5.9 | 2.2 | 3.9 | 4.5 | 0.0 | 0.0 | 2.7 |
| miscellaneous | 6.5 | 9.3 | 16.3 | 15.2 | 15.0 | 21.2 | 35.9 | 34.5 | 28.5 |
| **Σ sulfur classes** | **63.9** | **65.9** | **39.5** | **45.1** | **49.7** | **42.3** | **5.5** | **16.9** | **38.5** |
| *Σ sulfur classes, ng* | *12310* | *736* | *235* | *3416* | *291* | *1523* | *108* | *52* | *1175* |
| *grand total, ng/10 g* | *19260* | *1117* | *595* | *7568* | *586* | *3602* | *1965* | *307* | *3048* |

**Error bar on these values `[D]`:** none from digitisation — they are exact ratios of printed integers. The
only propagated uncertainty is the ±0.5 ng rounding of the compound rows (§3.4), which is negligible for class
totals above ~100 ng and reaches ±3 % relative for the smallest class totals.
---

## §5. DERIVED APPARENT Ea `[D]` — THE DELIVERABLE, AND ITS LIMITS

**Method.** For each compound and class with non-zero values at ≥ 3 die temperatures at fixed pH and fixed
sugar: leg-by-leg apparent Ea from `Ea = −R·ln(Y₂/Y₁)/(1/T₂ − 1/T₁)`, plus an ordinary least-squares
`ln Y` vs `1/T` fit over all three rungs with R². **Every value is computed twice** — once on the nominal
120/150/180 °C axis and once on the measured product-temperature axis of Table 1 (§1.3).

### 5.0 ★ READ THIS FIRST: THREE OF THE SIX LADDERS ARE BROKEN AT 150 °C

Before any Ea can be quoted, the shape of the raw data has to be stated. Grand totals, ordered 120 → 150 → 180 °C:

| ladder | 120 °C | 150 °C | 180 °C | shape |
|---|--:|--:|--:|---|
| glucose pH 7.5 | 682 | 2319 | 3479 | ✔ monotonic |
| glucose pH 6.5 | 372 | **116** | 606 | ✘ **150 °C collapses to 31 % of the 120 °C value** |
| glucose pH 5.5 | 134 | 364 | 564 | ✔ monotonic |
| **xylose pH 7.5** | 595 | 1117 | 19260 | ✔ **monotonic — the one clean ladder** |
| xylose pH 6.5 | 3602 | **586** | 7568 | ✘ **150 °C collapses to 16 %** |
| xylose pH 5.5 | 3048 | **307** | 1965 | ✘ **150 °C collapses to 10 %, and 180 < 120** |

★ **In three of the six ladders the entire 150 °C column collapses — not one compound, but every class at once,
by roughly an order of magnitude.** A whole-column, all-class collapse is the signature of **one bad extrusion
run or one bad isolate batch**, not of chemistry: no Maillard mechanism removes carbonyls, pyrazines, thiophenes,
thiazoles and furans simultaneously and then restores them all 30 °C later. Recall from §1.2 that **there is
exactly one extrusion per condition and no process replication**, so there is no way to test this and no error
bar that would have caught it.

**⇒ Consequence, applied throughout: every ladder containing a collapsed 150 °C cell yields a meaningless
"negative Ea then large positive Ea" pair. Those are NOT reported as activation energies below. They are
reported as REFUSED, with the collapse named as the reason. This removes glucose pH 6.5, xylose pH 6.5 and
xylose pH 5.5 from the Ea deliverable entirely — half the design.**

### 5.1 ★ PRIORITY COMPOUNDS — XYLOSE pH 7.5, THE ONLY CLEAN LADDER

Measured product temperatures **118 → 153 → 174 °C** (targets 120/150/180).

| # | compound / class | y₁₂₀ | y₁₅₀ | y₁₈₀ | Ea lo-leg (kJ/mol) | Ea hi-leg | whole fit | R² | Ea lo (meas. axis) | Ea hi (meas. axis) | whole (meas.) | R² |
|--:|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | ★ **2-furanmethanethiol (FFT)** | 66 | 341 | 5346 | **75.7** | **146.3** | 107.7 | 0.966 | **65.0** | **207.6** | 107.6 | 0.899 |
| 2 | ★ **2-methyl-3-furanthiol (MFT)** | 44 | 171 | 3038 | **62.6** | **152.9** | 103.5 | 0.941 | **53.8** | **217.1** | 102.6 | 0.861 |
| 3 | **total S-substituted furans** | 123 | 512 | 8404 | 65.8 | 148.7 | 103.3 | 0.950 | 56.5 | 211.1 | 102.7 | 0.874 |
| 4 | **total aliphatic S compds** | 4 | 79 | 1470 | **137.5** | **155.4** | **145.6** | **0.999** | 118.1 | 220.6 | 148.7 | 0.971 |
| 5 | **total thiophenes** | 73 | 119 | 2049 | 22.5 | 151.2 | 80.8 | 0.827 | 19.4 | 214.7 | 77.7 | 0.713 |
| 6 | **total thiazoles** | 35 | 26 | 387 | −13.7 | 143.5 | 57.5 | 0.619 | −11.8 | 203.7 | 52.6 | 0.483 |
| 7 | **total pyrazines** | 60 | 76 | 2087 | 10.9 | 176.1 | 85.7 | 0.766 | 9.4 | 249.9 | 81.2 | 0.641 |
| 8 | 2-furfural | 21 | 66 | 1670 | 52.8 | 171.7 | 106.6 | 0.907 | 45.3 | 243.8 | 104.6 | 0.814 |
| 9 | thiophene | 4 | 10 | 178 | 42.2 | 153.0 | 92.4 | 0.894 | 36.3 | 217.2 | 90.3 | 0.796 |
| 10 | 2-methylthiophene | 62 | 106 | 1487 | 24.7 | 140.4 | 77.1 | 0.844 | 21.2 | 199.3 | 74.4 | 0.733 |
| 11 | 2-methyl-3-thiophenethiol | 2 | 1 | 49 | −32.0 | 206.8 | 76.2 | 0.553 | −27.4 | 293.6 | 68.5 | 0.416 |
| 12 | 2-methylthiazole | 29 | 14 | 195 | −33.6 | 140.0 | 45.0 | 0.450 | −28.8 | 198.7 | 39.2 | 0.317 |
| 13 | total non-S furans | 150 | 101 | 2285 | −18.2 | 165.8 | 65.1 | 0.603 | −15.7 | 235.3 | 59.3 | 0.467 |
| 14 | 2-pentylfuran | 116 | 28 | 517 | −65.5 | 155.0 | 34.3 | 0.227 | −56.3 | 220.0 | 26.3 | 0.124 |
| 15 | hexanal | 75 | 90 | 1028 | 8.4 | 129.4 | 63.2 | 0.768 | 7.2 | 183.8 | 60.0 | 0.644 |
| 16 | total carbonyls | 16 | 26 | 249 | 22.4 | 120.1 | 66.6 | 0.850 | 19.2 | 170.5 | 64.4 | 0.740 |
| — | *GRAND TOTAL (reference)* | *595* | *1117* | *19260* | *29.0* | *151.3* | *84.4* | *0.853* | *24.9* | *214.8* | *81.7* | *0.744* |

**Refused on this ladder for <3 non-zero rungs:** FFT and MFT at pH 6.5 and 5.5, aliphatic S at pH 6.5 and 5.5,
total S-substituted furans at pH 5.5, thiazoles at pH 5.5, 2-methylthiophene and 2-methyl-3-thiophenethiol at
pH 5.5. `[NEG]` These are true absences (§3.4), not missing data.

### 5.2 GLUCOSE pH 7.5 AND pH 5.5 — the two surviving glucose ladders

| compound / class | pH | y₁₂₀ | y₁₅₀ | y₁₈₀ | Ea lo | Ea hi | whole | R² | whole (meas.) |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| **2-furanmethanethiol (FFT)** | 7.5 | 18 | 77 | 106 | **67.0** | **17.0** | 44.4 | 0.905 | 42.3 |
| **2-methyl-3-furanthiol (MFT)** | 7.5 | 15 | 37 | 31 | **41.6** | **−9.4** | 18.5 | 0.616 | 17.2 |
| total S-substituted furans | 7.5 | 33 | 114 | 137 | 57.2 | 9.8 | 35.7 | 0.873 | 34.0 |
| total aliphatic S compds | 7.5 | 141 | 614 | 496 | 67.8 | −11.3 | 32.0 | 0.665 | 29.9 |
| total alicyclic S compds | 7.5 | 12 | 33 | 72 | 46.6 | 41.5 | **44.3** | **0.999** | 42.9 |
| total thiophenes | 7.5 | 95 | 456 | 876 | 72.3 | 34.7 | 55.3 | 0.963 | 53.1 |
| total thiazoles | 7.5 | 68 | 157 | 397 | 38.6 | 49.3 | **43.4** | **0.995** | 42.4 |
| total pyrazines | 7.5 | 61 | 372 | 855 | 83.4 | 44.2 | 65.6 | 0.972 | 63.2 |
| total pyrroles | 7.5 | 79 | 159 | 333 | 32.2 | 39.3 | **35.4** | **0.997** | 34.5 |
| total oxazoles | 7.5 | 74 | 223 | 164 | 50.9 | −16.3 | 20.4 | 0.529 | 18.8 |
| 3-mercapto-2-butanone | 7.5 | 128 | 560 | 438 | 68.0 | −13.1 | 31.3 | 0.645 | 29.2 |
| 2(or 3)-thiophenethiol | 7.5 | 22 | 196 | 244 | 100.8 | 11.6 | 60.4 | 0.848 | 57.4 |
| thiazole | 7.5 | 4 | 17 | 82 | 66.7 | 83.6 | **74.4** | **0.996** | 72.5 |
| thiophene | 7.5 | 5 | 26 | 72 | 76.0 | 54.1 | 66.1 | 0.991 | 63.9 |
| methylpyrazine | 7.5 | 25 | 138 | 274 | 78.8 | 36.4 | 59.6 | 0.960 | 57.2 |
| ethylpyrazine | 7.5 | 19 | 140 | 266 | 92.1 | 34.1 | 65.8 | 0.940 | 63.1 |
| *GRAND TOTAL* | *7.5* | *682* | *2319* | *3479* | *56.4* | *21.6* | *40.6* | *0.943* | *38.9* |
| **2-furanmethanethiol (FFT)** | 5.5 | 12 | 16 | 32 | 13.3 | 36.8 | **23.9** | **0.926** | 24.9 |
| total S-substituted furans | 5.5 | 12 | 22 | 32 | 27.9 | 19.9 | **24.3** | **0.991** | 25.5 |
| total thiophenes | 5.5 | 16 | 16 | 133 | 0.0 | 112.5 | 51.0 | 0.714 | 52.7 |
| total pyrazines | 5.5 | 13 | 65 | 124 | 74.2 | 34.3 | 56.1 | 0.960 | 58.9 |
| total aliphatic S compds | 5.5 | 67 | 190 | 90 | 48.1 | −39.7 | 8.3 | 0.098 | 9.2 |
| ethylpyrazine | 5.5 | 8 | 32 | 108 | 63.9 | 64.6 | **64.2** | **1.000** | 67.2 |
| methylpyrazine | 5.5 | 5 | 12 | 12 | 40.4 | 0.0 | 22.1 | 0.785 | 23.3 |
| 2-methylthiophene | 5.5 | 3 | 6 | 12 | 32.0 | 36.8 | **34.2** | **0.998** | 35.7 |
| *GRAND TOTAL* | *5.5* | *134* | *364* | *564* | *46.1* | *23.3* | *35.7* | *0.968* | *37.5* |

**Refused (collapsed 150 °C column, §5.0): the whole of glucose pH 6.5 and xylose pH 6.5 and pH 5.5.**

### 5.3 ★ THE CAVEAT CHAIN — LONGER FOR EXTRUSION THAN FOR A SEALED VESSEL

**(i) Yield at the die exit is not a rate constant, and the time axis is only partly controlled.**
An Arrhenius slope needs `k`, i.e. yield **per unit time** with time held fixed. What is measured here is the
accumulated yield at the die of a flowing, reacting, expanding melt.

Partial credit where due, correcting the brief's expectation: **residence time IS reported — "modal residence
time, 35 s"** — and because **feed rate (600 g/min) and screw speed (350 rpm) were held constant across all 18
runs**, the mean residence time is plausibly near-constant. That is materially better than an uncontrolled time
axis. **But:** the value is quoted **once for the whole study, not per condition**; it is a *modal* value from
an RTD, with no spread given; and melt viscosity — the main determinant of fill and hence of RTD — changes
strongly with temperature and moisture, which is precisely why die pressure falls from 2.87 to 1.24 MPa across
a single ladder. **⇒ Time is approximately, not demonstrably, constant. These remain apparent yield-response
coefficients, not barriers.**

**(ii) ★ Die volatilisation is temperature-dependent and biases apparent Ea DOWNWARD.**
The authors raise this themselves, in their own opening sentence — verbatim:

> *"Unfortunately, the short cooking times limit flavor formation and, in addition, **a proportion of flavor
> compounds that are formed are volatilized at the die** (1, 2)."*

At the die the melt flashes from ~2–3 MPa to atmospheric pressure and superheated water boils off explosively;
this is what expands the product. **The hotter the die, the more violent the flash and the greater the fraction
of volatiles lost.** Loss therefore rises with temperature, suppressing `Y₁₈₀` relative to `Y₁₂₀` and
**biasing every apparent Ea in this paper DOWNWARD. Every Ea in §5.1 and §5.2 is a LOWER BOUND.**

★ **The direction of this bias is the single most useful fact in the section, because it is favourable:
the large positive thiol Ea values survive a bias that pushes the other way, and §7 turns this into the
central argument.**

**(iii) Moisture and shear change with temperature — and shear runs backwards.**
Feed moisture is fixed at 18 %, but in-barrel and post-die moisture are not (`[NEG]` product moisture never
reported), and more water flashes off at higher die temperature, so the 180 °C extrudates are the driest.
Worse, **SME falls monotonically as die temperature rises in all six ladders** (§1.4): the hottest runs are also
the lowest-shear, lowest-pressure runs. Temperature and mechanical energy are anti-correlated by design, and
nothing in the study separates them.

**(iv) Semi-quantitative scale.** Per §2.5, the constant `f` cancels in every ratio and therefore in every Ea
here. This caveat costs the Ea values nothing — but it forbids ever attaching a ng or ppb magnitude to them.

**(v) Precursor depletion.** Cysteine and sugar are both at 0.044 mol/kg and both are consumed. At the top of a
ladder the reaction is further along, so the incremental yield per degree understates the true rate ratio — a
second downward bias on apparent Ea, in the same direction as (ii).

**(vi) The extrudate pH is not the feed pH, and it moves along the ladder.** The x-axis of a "fixed-pH" ladder
is not fixed: glucose pH 6.5 runs at extrudate pH 4.06 → 4.36 → 4.67, a **0.6-unit drift**, and glucose pH 7.5
drifts 6.10 → 5.83 the other way. Every Ea here is a *joint* response to temperature **and** a correlated pH
drift. (Xylose pH 7.5 is the best behaved: 5.60 → 5.50 → 5.40.)

**(vii) One run per condition (§1.2).** No process replication anywhere; §5.0 shows three of six ladders carry
a probable bad run.

### 5.4 ★ VERDICT ON THE Ea BLOCK

**These are the weakest class of activation energies in the wave, and I recommend `PRIOR-ONLY` for the absolute
values and `RATIO-ONLY` for the contrasts.** Concretely:

- **No Ea from this paper should be fitted to, or used to set, a rate constant in the repo.** They are not
  barriers. They are yield-response coefficients contaminated by die loss, shear anti-correlation, pH drift and
  a partly-uncontrolled time axis, from a design with no process replication.
- **What is specifically rescued, and why:** the *ratio* results of §5.1 and §7 — the excess Ea of one class
  over the bulk on the same leg of the same ladder. In a ratio taken at one condition pair, the response
  factor, the die-loss factor common to all volatiles, the residence time, the shear and the pH drift **all
  cancel to first order**, because they act on numerator and denominator alike. That is a genuinely stronger
  statement than any single Ea in the table, and it is the only thing here I would carry forward.
- **The one absolute number I would even offer as a weak prior:** xylose pH 7.5 total aliphatic sulfur
  compounds, whole-ladder **145.6 kJ/mol at R² = 0.999** — the only near-perfectly log-linear sulfur ladder in
  the paper, and one whose two legs (137.5, 155.4) agree within 13 %. Even this is a lower bound per (ii) and
  it is `PRIOR-ONLY`.
---

## §6. THE pH AXIS `[D]` — REPORTED AS RATIOS

**The repo has almost no matrix-phase pH evidence, so this is the block with the broadest reach.** Reported as
fold-changes at fixed target temperature, relative to the pH 5.5 feed. Ratios of the same compound at the same
temperature cancel the response factor exactly (§2.5).

⚠ **Two framing warnings before the numbers.**
**(a) The axis is FEED pH, and the extrudate pH is very different** (§1.3). The three glucose "pH" levels
correspond to measured extrudate pH of roughly 5.8–6.1 / 4.1–4.7 / 3.3–3.6, and the xylose ones to 5.4–5.6 /
3.9–4.0 / 3.1–3.1. **The real span is about 2.5 pH units of extrudate acidity, all of it on the acid side of
neutral. Nothing here is a neutral or alkaline matrix measurement.**
**(b) The pH 7.5 / pH 5.5 column at 150 °C and 120 °C is contaminated by the collapsed-cell problem** (§5.0),
because the collapse falls in the pH 6.5 and 5.5 ladders. **Only the 180 °C row is clean for both sugars**, and
it is the row to trust.

### 6.1 Xylose — fold change vs pH 5.5, same target temperature

| class | 180 °C: 7.5/5.5 | 180 °C: 6.5/5.5 | 150 °C: 7.5/5.5 | 150 °C: 6.5/5.5 | 120 °C: 7.5/5.5 | 120 °C: 6.5/5.5 |
|---|--:|--:|--:|--:|--:|--:|
| ★ **S-substituted furans (MFT+FFT)** | **∞** (8404 vs 0) | **∞** (1773 vs 0) | **∞** (512 vs 0) | **∞** (142 vs 0) | 0.16 | 1.03 |
| **thiophenes** | **18.97** | **12.71** | 2.29 | 2.42 | 0.24 | 1.39 |
| **pyrazines** | **90.74** | **23.22** | 12.67 | 4.17 | 0.81 | 2.96 |
| aliphatic S compds | ∞ (1470 vs 0) | ∞ (106 vs 0) | ∞ (79 vs 0) | — (0 vs 0) | 0.40 | 13.80 |
| thiazoles | ∞ (387 vs 0) | ∞ (164 vs 0) | ∞ (26 vs 0) | ∞ (23 vs 0) | 0.43 | 1.96 |
| furans (non-S) | 2.04 | 1.97 | 0.74 | 1.25 | 0.17 | 1.20 |
| carbonyls | 22.64 | 8.00 | 3.71 | 1.71 | 0.31 | 0.75 |
| miscellaneous | 1.77 | 1.63 | 0.98 | 0.83 | 0.11 | 0.88 |
| **GRAND TOTAL** | **9.80** | **3.85** | 3.64 | 1.91 | 0.20 | 1.18 |

★ **The headline, from the clean 180 °C row: raising the feed pH from 5.5 to 7.5 in the xylose matrix turns MFT
and FFT on from literally zero to 3038 and 5346 ng/10 g.** They are **not detected at all** at pH 5.5 at either
150 or 180 °C. Even against the 9.8-fold rise in the grand total, the sulfur-substituted furans rise
**without bound** — an infinite ratio, because the denominator is a true measured absence (§3.4).

⚠ **The 120 °C row inverts** (most ratios < 1) — at 120 °C the pH 5.5 xylose run out-yields the pH 7.5 run
across the board. That row is one of the suspect ladders and I would not build on it.

### 6.2 Glucose — fold change vs pH 5.5, same target temperature

| class | 180 °C: 7.5/5.5 | 180 °C: 6.5/5.5 | 150 °C: 7.5/5.5 | 150 °C: 6.5/5.5 | 120 °C: 7.5/5.5 | 120 °C: 6.5/5.5 |
|---|--:|--:|--:|--:|--:|--:|
| ★ **S-substituted furans (MFT+FFT)** | **4.28** | 0.78 | 5.18 | 0.32 | 2.75 | 2.00 |
| **thiophenes** | **6.59** | 1.74 | 28.50 | 2.69 | 5.94 | 7.88 |
| **thiazoles** | **8.63** | 1.63 | 78.50 | 8.00 | ∞ (68 vs 0) | ∞ (50 vs 0) |
| **pyrazines** | **6.90** | 0.64 | 5.72 | 0.09 | 4.69 | 1.69 |
| aliphatic S compds | 5.51 | 0.59 | 3.23 | 0.06 | 2.10 | 0.36 |
| alicyclic S compds | 4.80 | 1.40 | 33.00 | 1.00 | 12.00 | 9.00 |
| furans (non-S) | 1.45 | 0.80 | 1.75 | 0.32 | 23.33 | 21.33 |
| pyrroles | 9.79 | 0.82 | 39.75 | 1.50 | 6.08 | 1.23 |
| oxazoles | 7.13 | 1.13 | 223.00 | 3.00 | 74.00 | 6.00 |
| carbonyls | 3.71 | 1.00 | 9.50 | 1.00 | — (10 vs 0) | — (2 vs 0) |
| miscellaneous | 2.09 | 1.45 | 14.60 | 1.40 | 3.62 | 3.62 |
| **GRAND TOTAL** | **6.17** | 1.07 | 6.37 | 0.32 | 5.09 | 2.78 |

### 6.3 ★ WHAT THE pH AXIS ACTUALLY SUPPORTS

1. ★ **pH 7.5 vs 5.5 at 180 °C is a robust, large, same-signed effect in both sugars**: grand total ×9.8
   (xylose) and ×6.2 (glucose). **This is the most transferable single number in the paper** — it is a ratio,
   at the cleanest temperature row, reproduced independently in two chemistries.
2. ★ **The thiols are the most pH-sensitive class of all, and in the matrix the effect runs OPPOSITE to
   aqueous systems.** The authors flag this themselves: *"In contrast to the extrudates, levels of these
   compounds **decrease** with increasing pH in aqueous systems (6, 29)."* And on FFT/MFT specifically:
   *"On heating aqueous systems in the pH range 3–7 at 145 °C for 20 min, levels of both compounds **decreased
   with increasing pH**, when ribose was the sugar."* **⇒ Any aqueous-fitted pH dependence for MFT/FFT in the
   repo has the WRONG SIGN for a low-moisture matrix. This is the single most consequential transfer finding
   in the dossier.**
3. **pH 6.5 is not intermediate — it is frequently the worst cell**, with several sub-1.0 ratios against
   pH 5.5. This is the collapsed-column artefact (§5.0) leaking into the pH axis, not a real inversion.
4. **Mechanistic gloss the authors offer, `[C]` from their citations, for the classes that rise with pH:**
   *"The formation of aldehydes via fragmentation of any sugar is favored by neutral/alkaline conditions"*
   (alicyclic S), and for thiazoles *"both the Strecker degradation of amino acids and sugar fragmentation lead
   to these precursors and are favored by elevated pH."* Conversely for 2-furfural, *"levels of 2-furfural
   increase with decreasing pH in aqueous systems (6, 36)"* — and indeed 2-furfural is the one major xylose
   species whose share is largest at pH 5.5 (43 % of total volatiles at 180 °C/pH 5.5).

---

## §7. ★ THE KANG SWITCH-ON CROSS-CHECK — REQUIRED SECTION

### 7.1 The reference finding

Kang et al. 2026 (aqueous TTCA/cysteine, pH 7, 120 min, sealed) found MFT and FFT strongly **non-Arrhenius**:
nearly flat 100→120 °C (×1.12 MFT, ×1.10 FFT), then a jump to 140 °C (×4.26 MFT, ×2.78 FFT); apparent Ea
climbing from ~6–7 kJ/mol on the low leg to **97.8 (MFT) / 69.2 (FFT)** kJ/mol on the 120→140 leg — **while the
sulfur class as a whole DECELERATED** (class Ea 57.5 → 35.2 kJ/mol).

Ames' ladder **brackets that window from above**: its lowest rung (118–128 °C measured) sits at Kang's upper
switch-on boundary, and it extends to 174–194 °C. Only the xylose pH 7.5 ladder is clean enough to ask the
question (§5.0), so the whole of §7 rests on that one ladder — **a real limitation, stated up front.**

### 7.2 Do MFT and FFT accelerate or decelerate across 120→150→180 °C?

**In raw yield: they accelerate, hard.** FFT 66 → 341 → 5346; MFT 44 → 171 → 3038. Leg Ea rises from 75.7 to
146.3 (FFT) and 62.6 to 152.9 (MFT) on the target axis; 65.0 → 207.6 and 53.8 → 217.1 on the measured axis.
**Both thiols roughly DOUBLE their apparent Ea between the two legs. There is a slope break, and it is upward.**

**But that comparison alone is worthless, because so does everything else.** The grand total does the same
thing (29.0 → 151.3), as do pyrazines (10.9 → 176.1), non-sulfur furans (−18.2 → 165.8) and hexanal
(8.4 → 129.4) — species with no cysteine chemistry at all. **A break shared by every class in the sample is a
property of the process, not of thiol chemistry**: the 180 °C run differs from the 150 °C run in shear
(598 → 501 kJ/kg), die pressure (2.05 → 1.24 MPa) and expansion as well as in temperature.

**The question therefore has to be asked as a ratio**, which is also all the semi-quant scale licenses (§2.5).

### 7.3 ★ THE RIGHT TEST: EXCESS Ea OVER THE BULK, xylose pH 7.5

Defining **excess Ea = (class leg Ea) − (grand-total leg Ea)** on the same leg of the same ladder. In this
difference the response factor, the residence time, the shear change, the pH drift **and the die-loss factor
common to all volatiles cancel to first order.** Positive means the class genuinely outruns the bulk.

| class / compound | excess Ea, **120→150 leg** | excess Ea, **150→180 leg** | reading |
|---|--:|--:|---|
| ★ **total aliphatic S compds** | ★ **+108.5** | **+4.0** | switches on hard on the low leg, then tracks the bulk |
| ★ **2-furanmethanethiol (FFT)** | ★ **+46.7** | **−5.1** | same |
| ★ **total S-substituted furans** | ★ **+36.7** | **−2.6** | same |
| ★ **2-methyl-3-furanthiol (MFT)** | ★ **+33.5** | **+1.6** | same |
| total thiophenes | −6.5 | −0.1 | tracks the bulk on both legs |
| total pyrazines | −18.1 | **+24.7** | the mirror image — accelerates only on the HIGH leg |
| total carbonyls | −6.7 | −31.2 | decelerates |
| total thiazoles | **−42.7** | −7.8 | *decelerates* on the low leg |
| total non-S furans | **−47.3** | +14.4 | decelerates on the low leg |
| miscellaneous (incl. hexanal) | −25.8 | −19.2 | decelerates |

*(measured-temperature axis gives the same signs and ordering throughout: FFT +40.1/−7.2, MFT +28.8/+2.3,
aliphatic S +93.2/+5.7, thiazoles −36.7/−11.1.)*

### 7.4 ★ THE ANSWERS

**(a) Is there a slope break, and where?**
**Yes — and once the global process trend is divided out, the break is on the 120→150 °C leg, not the
150→180 °C leg.** The thiols and mercaptoketones carry a large positive excess Ea across 120→150 (+33 to +109
kJ/mol) and then **fall to within ±5 kJ/mol of the bulk across 150→180**. In ratio terms the switch is thrown
between 118 and 153 °C measured, and it is finished by 153 °C. **This is the same window Kang located
(120→140 °C), reached by a completely different route: different phase, different moisture, different sugar,
different vessel, 35 seconds instead of 120 minutes.**

**(b) Does the thiol/class contrast Kang found have an analogue here?**
★ **Yes, and with the same sign.** Kang's contrast was *thiols accelerating while the sulfur class as a whole
decelerated*. Here, on the switch-on leg, **MFT (+33.5), FFT (+46.7) and the aliphatic mercaptoketones (+108.5)
accelerate, while the other sulfur classes do not** — thiazoles **−42.7** and thiophenes **−6.5**. The sulfur
branch splits in the matrix exactly as it split in water: the **direct H₂S + carbonyl products** (mercaptoketones
and the furanthiols) switch on, while the **ring-closure products** (thiazoles, thiophenes) do not. That the
split falls along a mechanistic seam — and not simply along "sulfur vs non-sulfur" — is what makes it credible.

**(c) ⚠ Could the die-volatilisation bias have manufactured or destroyed this break, and in which direction?**
**This is the decisive check, and it comes out favourably — the bias works AGAINST the finding, so the finding
survives it as a lower bound.**

- Die loss rises with die temperature (§5.3 ii), so it **suppresses high-temperature yields and biases all
  apparent Ea DOWNWARD**.
- Loss at the flashing die is governed by volatility. **MFT (LRI 873) and FFT (LRI 918) are among the *more*
  volatile species in the sample** — more volatile than the median of the thiophene, thiazole and pyrazine
  classes they are being compared against (class members run to LRI 1100–1650). So the thiols should be lost
  **preferentially**, and their *excess* Ea should be biased **more negative** than the bulk's.
- **⇒ The observed large POSITIVE thiol excess on the 120→150 leg is therefore a lower bound. The bias could
  not have manufactured it; it can only have shrunk it. The true switch-on is at least as sharp as measured.**
- **⇒ Conversely, the near-zero thiol excess on the 150→180 leg is NOT safe to interpret.** That is exactly
  where the bias is strongest and where it acts most selectively against the volatile thiols. The honest
  statement is *"the switch is thrown by ~150 °C"*, **not** *"the thiols stop outrunning the bulk above 150 °C"*
  — the latter may well be the die eating the signal.
- One weak counter-consideration, recorded for fairness: **hexanal (LRI 796) is the most volatile species in
  the xylose table and it tracks the bulk closely** (excess −25.8 / −19.2), which argues that volatility-ordered
  loss is not overwhelming. But hexanal is a lipid-oxidation product with an independent source term, so it is
  not a clean control.

### 7.5 What this cross-check is, and is not

**It IS:** the first matrix-phase, low-moisture, non-aqueous corroboration in the corpus that the MFT/FFT
switch-on window sits near 120–150 °C, and the first evidence that the thiol-vs-ring-closure split within the
sulfur branch is **not** an artefact of aqueous sealed-vessel conditions.

**It is NOT:** a measurement of the switch-on temperature (the ladder has three rungs, so the break can only be
localised to a 35 °C-wide leg), nor a confirmation of Kang's Ea magnitudes (Kang's 120→140 leg gives 97.8/69.2
kJ/mol; Ames' 120→150 leg gives 62.6/75.7 — the same order, but the two are not commensurable, per §2.5), nor
independent of the one clean ladder it rests on. **`RATIO-ONLY`.**

⚠ **And one genuine disagreement with Kang, recorded rather than smoothed over: the pH direction is opposite.**
Kang's aqueous system and the cited aqueous literature give MFT/FFT *falling* with rising pH; this matrix gives
them rising from zero to maximum across pH 5.5 → 7.5 (§6.3). **The switch-on transfers; the pH dependence does
not.**
---

## §8. VERIFIED NEGATIVES `[NEG]` AND FOUND POSITIVES

Each item from the brief, answered explicitly.

| # | question | answer | evidence |
|---|---|---|---|
| 1 | **Any absolute concentration anywhere?** | ★ **NO — verified negative.** All 118 compounds are semi-quantitative peak areas referenced to one internal standard and *expressed* in ng/10 g. No authentic-standard quantification, no response factors, no calibration curves. `grep -ic "response factor"` = 0, `"calibration"` = 0 | §2.4 |
| 2 | **Any residence time?** | ⚠ **FOUND POSITIVE, but only as a single global value** — *"modal residence time, 35 s"*, quoted once for the whole study among the operating conditions. **Not per condition**, no RTD spread, no measurement per run. Feed rate and screw speed were constant, so it is plausibly near-constant | §1.2, §5.3(i) |
| 3 | **Any measured (not target) die temperature?** | ★ **FOUND POSITIVE — this is the paper's most valuable methodological feature.** Table 1 gives a **measured product temperature for all 18 conditions**. ⚠ Note the header reads **"product temp"**, not "die temp" — it is the temperature of the product, measured at the die. Deviations run **−7 to +14 °C** from target | §1.3 |
| 4 | **Any replicate SD?** | ★ **NO — verified negative.** *"means of triplicate analyses"* with **no SD, no RSD, no error bars, no confidence intervals anywhere in the paper.** `grep -ic "standard deviation"` = 0. Worse, the triplicates are **isolates from one extrudate**, so even a hidden SD would be analytical, not process, variance. **No statistical test of any kind is applied to the volatile data** | §1.2, Table 4/5 footnote *a* |
| 5 | **Any H₂S?** | ★ **NO — verified negative.** Hydrogen sulfide appears 6 times, **exclusively as cited mechanism** from refs 26, 27, 30 (*"3-mercapto-2-butanone was first identified from butanedione and hydrogen sulfide by Takken et al."*). **Never measured.** The MS scan range begins at *m/z* 29, which excludes H₂S (*m/z* 34) from routine detection in any case | §2.2 |
| 6 | **Any α-dicarbonyl measurement?** | ⚠ **PARTIAL FOUND POSITIVE.** **2,3-pentanedione is measured as a headspace volatile in both tables** (glucose 6 ng at 180/7.5; xylose 151/14/10/48/13/26 across six cells) — a genuine α-dicarbonyl. **3-hydroxy-2-butanone** (acetoin, the reduction product) is also measured, glucose only. **But: no glyoxal, no methylglyoxal, no diacetyl/2,3-butanedione, no glucosone, and no derivatisation-based dicarbonyl assay.** The word "dicarbonyl" occurs once, mechanistically. **⇒ There is an α-dicarbonyl *signal* but not an α-dicarbonyl *measurement* in the repo's sense** | Tables 4, 5 |
| 7 | **Any pH measured in the extrudate rather than the feed?** | ★ **FOUND POSITIVE — and it is a large effect.** Table 1 gives measured **extrudate pH** for all 18 conditions, by the AOAC flour method (5 g extrudate + 20 mL tap water at pH 7.3, stirred 30 min). **The extrudate pH is 1.4–2.6 units BELOW the feed pH in every run** (e.g. feed 6.5 → extrudate 4.06). *"The pH of the product always decreased compared to that of the feed; the pH drop was more pronounced at higher temperatures and was usually greater for xylose"* | §1.3, §6 |

**Additional verified negatives found in the course of extraction:**

| item | status |
|---|---|
| screw configuration / profile | `[NEG]` not given — deferred to refs 18, 19 |
| barrel zone temperatures | `[NEG]` not given; `grep -ic "barrel"` = 0 |
| die dimensions / number of die holes | `[NEG]` not given |
| product moisture (out) | `[NEG]` not given; only feed moisture 18 % |
| expansion ratio / density | `[NEG]` not quantified; `grep -ic "expansion"` = 0 |
| residual cysteine / residual sugar | `[NEG]` not measured — no precursor mass balance |
| browning / colour / melanoidin | `[NEG]` not measured in this paper (colour is in ref 17, not here) |
| run order, batch structure, QC or bracketing samples | `[NEG]` none stated — the basis of the cross-sugar refusal, §2.5 item 4 |
| trapping/recovery efficiency for any analyte | `[NEG]` not determined; blanks were run but no recovery |
| odour thresholds | `[C]` cited only, from refs 7 and 8 (FFT/MFT 0.0025 ng/L in air; 2-methyl-3-thiophenethiol 0.0032–0.0128 ng/L; 3-mercapto-2-butanone 0.2–0.8; 3-mercapto-2-pentanone 0.05–0.2) |
| any kinetic model, rate constant or Ea fitted by the authors | `[NEG]` **none — the paper contains no kinetics at all.** Every Ea in this dossier is `[D]`, derived here and never printed |

---

## §9. CONSOLIDATED PARAMETER TABLE

| # | parameter | value | units | condition | provenance | source anchor |
|--:|---|---|---|---|---|---|
| 1 | cysteine loading | 0.044 | mol/kg feed | all runs | `[M]` | p.1885, Prep. of Extrudates |
| 2 | sugar loading (glucose or xylose) | 0.044 | mol/kg feed | all runs | `[M]` | p.1885 |
| 3 | cysteine : sugar molar ratio | 1 : 1 | — | all runs | `[D]` | from #1, #2 |
| 4 | starch moisture | 11.1 | % w/w | feed material | `[M]` | p.1885, Materials |
| 5 | starch nitrogen | <0.15 | % w/v | feed material | `[M]` | p.1885 |
| 6 | in-barrel moisture | 18 | % | all runs | `[M]` | p.1885 |
| 7 | NaOH dose for feed pH 7.5 / 6.5 / 5.5 | 34.2 / 23.4 / 19.7 | g per L water feed | all runs | `[M]` | p.1885 |
| 8 | feed rate | 600 | g/min | all runs | `[M]` | p.1885 |
| 9 | screw speed | 350 | rpm | all runs | `[M]` | p.1885 |
| 10 | modal residence time | 35 | s | all runs (single global value) | `[M]` | p.1885 |
| 11 | measured product temperature | 115–194 (18 values) | °C | per condition | `[M]` | Table 1, p.1886 |
| 12 | product-temp deviation from target | −7 to +14 | °C | per condition | `[M]`/`[D]` | Table 1; authors quote −5 to +14 |
| 13 | measured extrudate pH | 3.09–6.10 (18 values) | — | per condition | `[M]` | Table 1 |
| 14 | feed→extrudate pH drop | 1.4–2.6 | pH units | per condition | `[D]` | from Table 1 |
| 15 | die pressure | 0.50–2.87 | MPa | per condition | `[M]` | Table 1 |
| 16 | specific mechanical energy (SME) | 371–748 | kJ/kg | per condition | `[M]` | Table 1 |
| 17 | volatile yields, glucose | 82 rows × 9 conditions | ng/10 g extrudate (semi-quant) | per condition | `[M]` | Table 4, p.1888–89 |
| 18 | volatile yields, xylose | 41 rows × 9 conditions | ng/10 g extrudate (semi-quant) | per condition | `[M]` | Table 5, p.1890 |
| 19 | compound counts | 80 in 11 classes (glucose); 38 in 9 (xylose); +3 misc | — | — | `[M]`/`[D]` | p.1887; recount confirms |
| 20 | sulfur-containing compound counts | 44 (glucose), 21 (xylose) | — | — | `[M]`/`[D]` | p.1887; recount confirms |
| 21 | sulfur-class RA range, glucose | 51.2 – 71.6 | % of total | across 9 conditions | `[D]` | §4.2; authors quote 51–72 |
| 22 | sulfur-class RA range, xylose | 5.5 – 65.9 | % of total | across 9 conditions | `[D]` | §4.2; authors quote 5–66 |
| 23 | ★ apparent Ea, **FFT**, xylose pH 7.5 | **75.7** (120→150), **146.3** (150→180), 107.7 whole (R²=0.966) | kJ/mol | target axis | `[D]` | §5.1 |
| 24 | ★ apparent Ea, **FFT**, xylose pH 7.5 | **65.0** / **207.6** / 107.6 (R²=0.899) | kJ/mol | measured axis | `[D]` | §5.1 |
| 25 | ★ apparent Ea, **MFT**, xylose pH 7.5 | **62.6** / **152.9** / 103.5 (R²=0.941) | kJ/mol | target axis | `[D]` | §5.1 |
| 26 | ★ apparent Ea, **MFT**, xylose pH 7.5 | **53.8** / **217.1** / 102.6 (R²=0.861) | kJ/mol | measured axis | `[D]` | §5.1 |
| 27 | apparent Ea, total aliphatic S, xylose pH 7.5 | 137.5 / 155.4 / **145.6 (R²=0.999)** | kJ/mol | target axis | `[D]` | §5.1 — best-conditioned fit in the paper |
| 28 | apparent Ea, total thiophenes, xylose pH 7.5 | 22.5 / 151.2 / 80.8 (R²=0.827) | kJ/mol | target axis | `[D]` | §5.1 |
| 29 | apparent Ea, total thiazoles, xylose pH 7.5 | −13.7 / 143.5 / 57.5 (R²=0.619) | kJ/mol | target axis | `[D]` | §5.1 |
| 30 | apparent Ea, total pyrazines, xylose pH 7.5 | 10.9 / 176.1 / 85.7 (R²=0.766) | kJ/mol | target axis | `[D]` | §5.1 |
| 31 | apparent Ea, grand total, xylose pH 7.5 | 29.0 / 151.3 / 84.4 (R²=0.853) | kJ/mol | target axis | `[D]` | §5.1 |
| 32 | ★ **excess Ea (class − bulk), 120→150 leg, xylose pH 7.5** | aliphatic S **+108.5**, FFT **+46.7**, S-subst furans **+36.7**, MFT **+33.5**; thiophenes −6.5, thiazoles **−42.7** | kJ/mol | ratio, target axis | `[D]` | §7.3 |
| 33 | ★ **excess Ea, 150→180 leg, xylose pH 7.5** | all thiol/aliphatic-S classes within **±5**; pyrazines +24.7 | kJ/mol | ratio, target axis | `[D]` | §7.3 |
| 34 | ★ pH effect, grand total, 180 °C | **×9.80** (xylose), **×6.17** (glucose) | fold, pH 7.5 / pH 5.5 | 180 °C | `[D]` | §6.1, §6.2 |
| 35 | ★ pH effect, MFT + FFT, xylose, 150 and 180 °C | **0 → 8404** (180 °C), **0 → 512** (150 °C) — unbounded ratio | ng/10 g | pH 5.5 → 7.5 | `[D]` | §6.1 |
| 36 | FFT : MFT area ratio | 1.2 – 3.4 (both sugars, all detected cells) | — | — | `[D]` | §3.5 ⚠ cross-compound, not molar |
| 37 | 2-furfural share of total volatiles | 42.9 | % | xylose, 180 °C, pH 5.5 | `[D]` | §3.5; authors quote 43 % |
| 38 | number of extrusion replicates | **1** | runs per condition | all | `[M]`/`[NEG]` | §1.2 |
| 39 | number of analytical replicates | 3 | isolates per condition | all | `[M]` | §1.2 |

---

## §10. USABILITY VERDICTS

Recorded per block **with the reason**, so a later wave does not re-ingest a refused row.

| block | verdict | reason |
|---|---|---|
| **§1.3 Table 1 — measured product temperature, extrudate pH, die pressure, SME (18×4)** | ★ **USE** | Directly measured process instrumentation, fully tabulated, internally consistent, raster-verified cell by cell. The only unconditional USE in the dossier. **The measured-temperature column is what lets every Ea here be recomputed on a real axis.** |
| **§1.3 extrudate-pH column specifically** | ★ **USE** | The corpus's only measured low-moisture-matrix pH offset: feed→extrudate drops 1.4–2.6 units, larger at higher T. Use it to correct any "feed pH" assumption elsewhere in the extrusion lane. |
| **§3.2 / §3.3 Tables 4 and 5 — 123 compound rows × 9 conditions** | **USE-Q** (qualified) | Complete and triple-verified, and Table 5 passes all 99 subtotal checks. **Qualified because the values are semi-quant area ratios wearing ng clothing (§2.4).** Usable as *shape*: ratios, fold-changes, Ea, RA. **Never as concentrations.** Table 4 additionally carries 5 known internal arithmetic errors (§3.4). |
| **§4.2 relative abundances (Figures 1–4 reconstructed)** | ★ **USE** | Exact ratios of printed integers, zero digitisation error, and validated against 14 RA figures the authors quote independently (13 exact, 1 revealing a pH mislabel in their text). |
| **§5.1 Ea, xylose pH 7.5 ladder** | **PRIOR-ONLY** for absolute values | The one clean ladder, but still: die-loss bias (downward), shear anti-correlated with T, one run per condition, ~21 °C real upper leg vs 30 °C nominal. **Order-of-magnitude prior only — never fit to.** |
| **§5.2 Ea, glucose pH 7.5 and pH 5.5 ladders** | **PRIOR-ONLY**, weaker still | Same caveats, plus Table 4's arithmetic errors and 28 MS-only identifications including the largest thiophene. |
| **§5 Ea, glucose pH 6.5 / xylose pH 6.5 / xylose pH 5.5** | ★ **REFUSE** | **Whole-column collapse of the 150 °C cell across every chemical class simultaneously (§5.0).** Almost certainly one bad extrusion run or isolate batch; unfalsifiable because there is no process replication. Any Ea from these three ladders is an artefact. **Do not re-ingest.** |
| **§6 pH axis, 180 °C row, both sugars** | ★ **RATIO-ONLY — the highest-value transferable block** | Clean row, large same-signed effect in two independent chemistries (×9.8 and ×6.2). Response factor cancels. The repo's best matrix-phase pH evidence. |
| **§6 pH axis, 150 °C and 120 °C rows** | **REFUSE** | Contaminated by the same collapsed 150 °C cells; the 120 °C row inverts implausibly. |
| **§6.3(2) — pH SIGN REVERSAL for MFT/FFT vs aqueous** | ★ **USE (structural)** | Authors state it explicitly and cite aqueous work showing the opposite sign. **Actionable: an aqueous-fitted MFT/FFT pH dependence has the wrong sign in a low-moisture matrix.** |
| **§7 Kang switch-on cross-check** | ★ **RATIO-ONLY** | Genuine, and it survives the one bias that could have faked it (die loss pushes the other way, §7.4c). But it rests on a single ladder with three rungs, so it localises the break to a 35 °C-wide leg and no better. **Corroborative, never quantitative.** |
| **§7.4(b) thiol-vs-ring-closure split** | ★ **STRUCTURAL** | The mechanistic seam (direct H₂S+carbonyl products switch on; ring-closure products do not) reproduces in a matrix. Use as a structural constraint on the sulfur branch's topology, not as a number. |
| **§2.4 quantification basis** | **STRUCTURAL** | Records *why* everything above is ratio-only. |
| **§3.6 / §3.7 Tables 2 and 3 — sensory descriptors** | **STRUCTURAL** | Free-choice profiling, 5 untrained assessors, no intensity scale, no replication, no statistics. Direction only — and note it contradicts the chemistry on the pH axis (§3.7). |
| **cross-sugar comparisons (xylose vs glucose magnitudes)** | ★ **REFUSE** | No batch, run-order or QC statement anywhere (§2.5 item 4); two tables, two grand totals, two unknown response-factor sets. The paper makes such claims; this dossier does not carry them. |
| **FFT/MFT "about twice" claim** | **REFUSE as quantitative** | Cross-compound area ratio with two unknown response factors (§2.5 item 3). Reproduces arithmetically; means nothing molar. |
| **any absolute ng or ppb from this paper** | ★ **REFUSE** | Semi-quant (§2.4). **Must never enter the repo as a concentration or a benchmark row.** |
| **any H₂S, α-dicarbonyl, browning or precursor-depletion parameter** | ★ **REFUSE — not measured** | §8 items 5, 6 and the negatives table. |

---

## §11. DECLARED GAPS

1. ★ **No process replication anywhere.** One extrusion per condition; `n = 3` is analytical. This is the root
   cause of §5.0 and it caps everything the paper can support. **Half the temperature design (3 of 6 ladders)
   is lost to a single suspected bad run that cannot be identified or excluded on evidence.**
2. ★ **No absolute quantification.** No response factors, no authentic-standard calibration, no recovery
   determination. The ng/10 g unit is an area ratio in disguise.
3. ★ **Only one usable ladder for the wave's central question.** The entire §7 switch-on answer rests on
   xylose pH 7.5. A single-ladder result cannot be internally replicated.
4. **Three rungs is the minimum for a slope break and too few to locate one.** The break can be placed on the
   120→150 leg but not within it. Kang's window (120–140) is *inside* that leg; this paper cannot resolve
   further.
5. **Temperature is confounded with shear and pressure by construction.** SME and die pressure fall
   monotonically as die temperature rises in all six ladders. Nothing separates them. A "temperature effect"
   here is a temperature-plus-lower-shear effect.
6. **The pH axis is feed pH; extrudate pH is 1.4–2.6 units lower and drifts along each ladder** (up to 0.6
   units within a single "fixed-pH" ladder). No ladder is genuinely isothermal in pH.
7. **Die volatilisation is acknowledged but never quantified.** The authors cite it in their first paragraph
   and never return to it. No die-loss fraction is measured, at any temperature, for any compound — so the
   downward Ea bias can be signed but not corrected.
8. **The headspace measurement is of a re-hydrated extrudate purged at 37 °C for 1 h**, not of the extrudate.
   Anything that changes during that hour is measured after it.
9. **The internal standard is spiked onto the trap after trapping**, so it corrects desorption and MS drift but
   not headspace partitioning or trapping efficiency — the very step where analytes discriminate (§2.1).
10. **28 of 82 glucose rows and 8 of 41 xylose rows are MS-only (tentative)**, including the largest single
    glucose thiophene, `2(or 3)-thiophenethiol`, whose substitution position the authors could not even assign.
11. **Table 4 contains 5 internal arithmetic errors**, one of them a diagnosable column-shift typesetting error
    (§3.4). Table 5 is clean.
12. **The published text misattributes its maximum xylose sulfur RA to pH 6.5; it is pH 7.5** (§3.5 #7).
13. **No H₂S, no α-dicarbonyl assay, no browning, no precursor mass balance** — so nothing here constrains the
    sulfur branch's *inputs*, only its volatile *outputs*.
14. **Screw configuration, barrel profile, die geometry, product moisture and expansion are all absent**,
    deferred to refs 18 and 19 (Bredie et al. 1997, 1998). ⇒ **Follow-up lead: those two papers are where the
    extrusion lane's missing process parameters live**, and ref 19 (*J. Agric. Food Chem.* **1998**, 46,
    1479–1487, extruded maize flour) is a plausible K6 candidate in its own right.
15. **No kinetics of any kind are reported by the authors.** Every Ea in this dossier is `[D]` — derived here,
    never printed — and must be labelled as such wherever it travels.

---

### Provenance key

`[M]` measured · `[C]` cited · `[F]` fitted by the authors *(nothing in this paper carries `[F]` — the authors
fitted nothing)* · `[D]` derived by me and never printed · `[NEG]` verified negative.

### Reproducibility of this extraction

Text layer `pdftotext -layout` and `pdftotext -bbox`; rasters `pdftoppm -r 600` (pp. 1, 2, 4, 5, 6) and
`-r 300` (pp. 7, 8). Every numeric cell of Tables 1, 4 and 5 was read three independent ways (text layer,
raster, bounding-box column assignment) and all three agree. All derived quantities were computed from the
transcribed integers and cross-validated against fourteen statistics the authors state independently in their
Results text (§3.5).
