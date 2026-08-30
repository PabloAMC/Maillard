# Hofmann & Schieberle 1998 (10.1021/jf9705983) — Wave K4a IDENTITY RECONCILIATION + complete table record, 2026-08-28

**Source PDF:** `data/articles/hofmann1998.pdf` (7 pp.). Born-digital ACS PDF, clean text layer,
two-column; the pdftotext `-layout` extraction interleaves the two columns but every table
survives intact and was read verbatim.

> **This file answers Wave K4a cross-check #1 and does NOT supersede the existing extraction.**
> Hofmann 1998's ten tables are already transcribed in
> `data/lit/extraction_dossiers/k3_final_parameter_inventory.md` §A.3.1–A.3.2 (lines ~210–260).
> I re-read the paper end to end and **verified every value in that existing record against the
> PDF — all of them are correct.** What follows is (a) the identity answer, (b) the
> temperature/water-content answer, (c) a complete standalone transcription of all ten tables so
> the record is verifiable without the inventory, and (d) **three discrepancies** the
> reconciliation turned up.

---

## 0. ANSWER TO CROSS-CHECK #1 — RESOLVED, NO DUPLICATE ON DISK

**Question asked:** is DOI 10.1021/jf9705983, Hofmann & Schieberle 1998, *"Quantitative Model
Studies on the Effectiveness of Different Precursor Systems…"*, already on disk under another
name, and is it the same paper as the corpus's `hofmann1998` extraction?

**Answer: YES to the second, NO to the first. `data/articles/hofmann1998.pdf` IS that paper.
There is no second copy under another filename, and there is no ambiguity to resolve.**

| field | value |
|---|---|
| Authors | **T. Hofmann and P. Schieberle\*** |
| Title | **"Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol"** |
| Venue | *J. Agric. Food Chem.* **1998**, **46** (1), **235−241** |
| DOI | **10.1021/jf9705983** |
| Affiliation | Deutsche Forschungsanstalt für Lebensmittelchemie, Lichtenbergstrasse 4, D-85748 Garching, Germany |
| Tables | **10** (Tables 1–10) |
| Figures | 3 (all mechanism schemes — **no numeric figures, nothing to digitise**) |

**Filename audit of the whole corpus.** The five other `hofmann*` PDFs are distinct papers and
none of them is jf9705983:

| file | what it actually is | relationship to jf9705983 |
|---|---|---|
| `hofmann1998.pdf` | **jf9705983 itself** | — |
| `hofmann1995.pdf` | Hofmann & Schieberle 1995, ribose/cysteine AEDA, *JAFC* 43, 2187−2194 | cited as a reference by jf9705983 |
| `hofmann1996.pdf` | the 6 °C organic-solvent artefact-control study | declared "neither" in the FIT declaration |
| `hofmann1996b.pdf` | the sotolon / oxidant-series paper | **its Tables 1 and 2 are verbatim re-publications of jf9705983 Tables 2, 8 and 10** — already flagged as a double-counting hazard at `FIT_HOLDOUT_DECLARATION.md:142`. **Confirmed correct: that flag is real and must stay** |
| `hofmann2000.pdf`, `hofmann2000b.pdf`, `hofmann2001.pdf` | later papers | unrelated |

**No action needed on identity.** The corpus's `hofmann1998` extraction is correctly attributed.

---

## 1. ★ ANSWER TO CROSS-CHECK #1b — THE TEMPERATURE / WATER-CONTENT SERIES **DOES NOT EXIST**

**This is the important finding, and it is negative.** The brief hoped that extracting this
paper's "temperature and water-content tables" would close a declared top gap. **It cannot,
because the paper has no such series.** The abstract's phrase *"model systems varying in
temperature, pH value, or water content"* (p. 235) promises more than the Experimental delivers.

### 1.1 There are exactly TWO process conditions in the entire paper, and they are confounded

**Anchor: Experimental Procedures, "Model Experiments", p. 236.** Verbatim, complete:

> "**In a first series of experiments**, the reactants (cf. amounts detailed in the tables) were
> dissolved in **phosphate buffer (0.5 mol/L) at pH 3.0, 5.0, or 7.0** and thermally treated for
> **20 min at 145 °C** in a laboratory autoclave (200 mL; Type II; Roth, Germany). **In a second
> series of experiments**, the reactants were intimately mixed with **silica gel (3.0 g
> containing 300 µL of the respective phosphate buffer; 2 mol/L)** and reacted for **6 min at
> 180 °C**."

That is the whole of it. Tabulated:

| | **Arm A — "aqueous"** | **Arm B — "dry-heated"** |
|---|---|---|
| Medium | phosphate buffer solution, **0.5 mol/L** | **3.0 g silica gel** wetted with **300 µL** of **2 mol/L** phosphate buffer |
| Volume | **100 mL** (cysteine/carbohydrate systems) or **50 mL** (fed-precursor systems) | 300 µL liquid on 3.0 g solid |
| **Temperature** | **145 °C** | **180 °C** |
| **Time** | **20 min** | **6 min** |
| pH | **3.0, 5.0 or 7.0** | **5.0 only** |
| Vessel | laboratory autoclave, 200 mL, Roth Type II | not stated for the dry arm |
| Which tables use it | T1, T2 (pH cols), T3, T4, T5, T6, T7, T8 (expts 1–3), T9, T10 | **T2 (the "5.0ᶜ" column) and T8 (experiment 4) — and nothing else** |

### 1.2 Why this cannot serve as a water-content axis

Between Arm A and Arm B, **four variables change simultaneously**:

| variable | Arm A | Arm B | change |
|---|---|---|---|
| water content | ~100 mL bulk water | 300 µL on silica | **~330× less liquid, and a solid support introduced** |
| temperature | 145 °C | 180 °C | **+35 K** |
| time | 20 min | 6 min | **3.3× shorter** |
| buffer molarity | 0.5 mol/L | 2 mol/L | **4× more concentrated** |

**There is no third condition and no intermediate point.** A two-level contrast in which
a_w, T, t and buffer strength all move together **cannot identify a water-content effect, a
temperature effect, or a time effect separately.** Any model term fitted to the A→B jump is
fitting the *combination*, and would be mis-attributed if labelled "a_w" or "temperature".

**Consequence for the declaration.** `FIT_HOLDOUT_DECLARATION.md:26` reads:

> | Hofmann 1998 T2 **dry 180 °C / 6 min** rows | **HOLD-OUT (new)** | a genuinely different regime (silica, 6 min, 180 °C); the only a_w extrapolation test the sulfur branch has |

The **role** (HOLD-OUT) is right and the parenthetical is honest about the regime change. But
the phrase **"the only a_w extrapolation test"** overstates what the rows can test: it is a
**confounded four-variable regime jump**, not an a_w axis. **Recommend the orchestrator amend
the wording** (I did not edit the file). The rows remain valuable as a *regime* hold-out — a
model that gets the A→B direction right is doing something real — but a **pass or fail there
localises nothing**, because four candidate causes are aliased.

### 1.3 The authors themselves describe it as a combined variable

> "To examine the influence of **water and the reaction temperature**, the three carbohydrates
> were heated under dry-heating conditions in the presence of cysteine. A small amount of
> concentrated buffer was added to maintain pH 5.0." (p. 237)

Note "water **and** the reaction temperature" — the authors never claim to separate them. Their
conclusion is a direction, not a coefficient:

> "In summary, the **highest amounts of FFT were obtained by dry-heating the ribose/cysteine
> mixture** (Table 2), whereas **MFT was preferentially formed from this model mixture under
> aqueous conditions at pH 3.0**. Substitution of ribose by hexoses generally gave different
> results. Besides lower yields, the formation of **FFT from both hexose/cysteine systems was
> highest under aqueous conditions at pH 5.0**, whereas the **MFT was formed from these mixtures
> preferentially under dry-heating conditions**." (pp. 237–238)

**★ That is a genuine and useful structural constraint — a SIGN REVERSAL between pentose and
hexose in response to the same regime change** — and it is worth ingesting *as a sign
constraint*. It is not a quantitative a_w term.

---

## 2. COMPLETE TABLE RECORD — all ten tables, verbatim

Shared footnote **b** (Table 1, and referenced by every other table): *"Data are mean values of
**triplicates and differed by not more than 10 %**."* — **this ±10 % is the paper's only error
statement and it applies to every number below.** Quantitation is **SIDA** (stable isotope
dilution assay) with [²H₂]-FFT and [²H₃]-MFT internal standards, 10 µg each in 0.5 mL pentane,
10 min equilibration, diethyl ether extraction (3 × , total 60 mL), Na₂SO₄ drying, mass
chromatography with response factors per Guth et al. 1995. **This is the highest-quality
quantitation tier in the corpus.** Provenance for every cell below: **[M]**.

### 2.1 Table 1 — FFT and MFT from cysteine + various carbohydrates
**Anchor: Table 1, p. 237.** Title as printed: *"Amounts of FFT and MFT Generated from Cysteine
and Various Carbohydratesᵃ"*. Footnote a: *"A solution (**100 mL**) of **cysteine (3.3 mmol)**
and the corresponding **carbohydrate (10.0 mmol)** was reacted at **pH 5.0** in a laboratory
autoclave."* Column header as printed: `amountᵇ (µg)`.
(⇒ cysteine 33 mM, sugar 100 mM, **1:3 cys:sugar**, 145 °C / 20 min.)

| carbohydrate | **FFT (µg)** | **MFT (µg)** |
|---|---|---|
| ribose | **12.1** | **19.8** |
| xylose | **9.6** | **14.3** |
| fructose | **3.2** | **2.5** |
| glucose | **2.8** | **1.9** |
| glucose-6-phosphate | **0.9** | **0.6** |
| rhamnose | **0.8** | **0.8** |
| maltose | **0.6** | **0.3** |

**[Z]** ribose/glucose ratio: FFT **4.32×**, MFT **10.42×**. MFT/FFT at pH 5, ribose: **1.636**.

### 2.2 Table 2 — influence of pH, plus the single dry-heated column
**Anchor: Table 2, p. 237.** Title as printed: *"Influence of pH on the Amounts of FFT and MFT
Generated from Cysteine and Carbohydratesᵃ"*. Footnote a: *"A solution (**100 mL**) of cysteine
(3.3 mmol) and the corresponding carbohydrate (10.0 mmol) was reacted in a laboratory
autoclave."* Footnote c: *"A mixture of **cysteine (3.3 mmol)** and the corresponding
**carbohydrate (10.0 mmol)** was **dry-heated for 6 min at 180 °C**."*
Column headers as printed: `FFTᵇ (µg) at pH  3.0 | 5.0 | 7.0 | 5.0ᶜ` and the same four for MFT.

| carbohydrate | **FFT pH 3.0** | **FFT pH 5.0** | **FFT pH 7.0** | **FFT 5.0ᶜ (DRY)** | **MFT pH 3.0** | **MFT pH 5.0** | **MFT pH 7.0** | **MFT 5.0ᶜ (DRY)** |
|---|---|---|---|---|---|---|---|---|
| **ribose** | **22.9** | **12.1** | **1.2** | **97.2** | **55.3** | **19.8** | **2.5** | **25.1** |
| **glucose** | **0.7** | **2.8** | **0.6** | **1.4** | **0.3** | **1.9** | **0.4** | **4.2** |
| **rhamnose** | **0.2** | **0.8** | **0.1** | **0.4** | **0.1** | **0.8** | **0.1** | **3.1** |

**The four "5.0ᶜ" cells per analyte are the ENTIRE water-content/temperature dataset of this
paper**, together with Table 8 experiment 4. Six numbers in all (3 carbohydrates × 2 analytes).

**[Z] the sign reversal, quantified** — dry/aqueous ratio at pH 5.0:

| carbohydrate | **FFT dry/aqueous** | **MFT dry/aqueous** |
|---|---|---|
| ribose (pentose) | **8.03×** ↑ | **1.27×** ↑ |
| glucose (hexose) | **0.50×** ↓ | **2.21×** ↑ |
| rhamnose (deoxyhexose) | **0.50×** ↓ | **3.88×** ↑ |

**Pentose FFT goes UP 8× while hexose FFT goes DOWN 2×, under the same regime change. MFT goes
up for all three but by very different factors.** This is the constraint worth keeping.

### 2.3 Table 3 — FFT from 5-carbon precursors + H₂S
**Anchor: Table 3, p. 237.** Title: *"Generation of FFT from 5-Carbon Precursors in the Presence
of Hydrogen Sulfideᵃ"*. Footnote a: *"An **aqueous solution (50 mL)** of **hydrogen sulfide
(1 mmol)** and the **precursor (1 mmol)** was reacted at **pH 5.0** in a laboratory
autoclave."* Column headers: `FFT | µgᵇ | %`.

| precursor | **FFT µg** | **FFT %** |
|---|---|---|
| ribose | **9.2** | **0.008** |
| 3-deoxyribosulose | **78.6** | **0.08** |
| **furan-2-aldehyde** | **550.8** | **0.48** |

**[Z]** furan-2-aldehyde is **59.9× ribose** and **7.0× the 3-deoxyosone**.

### 2.4 Table 4 — MFT from ribose or NF + H₂S
**Anchor: Table 4, p. 237.** Title: *"Generation of MFT from Ribose or NF in the Presence of
Hydrogen Sulfideᵃ"*. Footnote a: *"An aqueous solution (**50 mL**) of hydrogen sulfide
(**1 mmol**) and the precursor was heated at **pH 5.0** in a laboratory autoclave."*
Column header: `precursor (1 mmol) | MFT µgᵇ | %`.

| precursor (1 mmol) | **MFT µg** | **MFT %** |
|---|---|---|
| ribose | **15.1** | **0.01** |
| **NF** (norfuraneol) | **211.2** | **0.19** |

**[Z]** NF is **13.99× ribose** (the text says "by a factor of ≈14").

### 2.5 Table 5 — the FA/FFT and NF/MFT pairs
**Anchor: Table 5, p. 238.** Title: *"Amounts (Micrograms)ᵃ of the Pairs Furan-2-aldehyde
(FA)/FFT and NF/MFT Generated from Carbohydrates in the Presence of Cysteineᵇ"*. Footnote b:
*"An aqueous solution (**100 mL**) of cysteine (3.3 mmol) and the corresponding carbohydrate
(10.0 mmol) was reacted at **pH 5.0** in a laboratory autoclave."*

| carbohydrate | **FA (µg)** | **FFT (µg)** | **NF (µg)** | **MFT (µg)** |
|---|---|---|---|---|
| ribose | **67.5** | **12.1** | **54530.0** | **19.8** |
| glucose | **6.9** | **2.8** | **13.5** | **1.9** |
| rhamnose | **6.2** | **0.8** | **19.1** | **0.8** |

**★ The authors' own falsification of the NF→MFT-only hypothesis**, verbatim (p. 239):

> "the amounts of NF were **not well correlated** with those of MFT. From ribose, by a factor of
> **3000** more NF was formed compared to MFT, whereas in the glucose system only by a factor of
> **10** more NF was present compared to MFT. This discrepancy implied that the **MFT formation
> might not run exclusively via NF as the key intermediate.**"

**[Z]** exact ratios: ribose NF/MFT = **2754×**; glucose = **7.1×**; rhamnose = **23.9×**.
(The text's "3000" and "10" are the authors' rounding of 2754 and 7.1.)

### 2.6 Table 6 — FA and NF from hydroxyacetaldehyde + 2-oxopropanal (the C₂ + C₃ route)
**Anchor: Table 6, p. 239.** Title: *"Amounts of Furan-2-aldehyde (FA) and NF Generated from
Hydroxyacetaldehyde and 2-Oxopropanalᵃ"*. Footnote a: *"An **aqueous solution (50 mL)** of
**hydroxyacetaldehyde (1 mmol)** and **2-oxopropanal (1 mmol)** was reacted in a laboratory
autoclave."*

| pH | **FA µg** | **FA %** | **NF µg** | **NF %** |
|---|---|---|---|---|
| **3.0** | **101.2** | **0.11** | **153.1** | **0.13** |
| **5.0** | **364.5** | **0.40** | **885.1** | **0.78** |
| **7.0** | **1443.1** | **1.60** | **3610.5** | **3.18** |

**★ Monotonically increasing with pH, and the ONLY monotone-in-pH panel in the paper** — every
cysteine/carbohydrate panel peaks at pH 5 or rises toward pH 3. Authors' reading:
> "in all reaction mixtures **more NF than FA** had been formed … these data indicate
> **hydroxyacetaldehyde is the better nucleophile**."

### 2.7 Table 7 — mercapto-2-propanone from 2-oxopropanal + H₂S (the H₂S dose–response)
**Anchor: Table 7, p. 239.** Title: *"Generation of Mercapto-2-propanone from 2-Oxopropanal and
Hydrogen Sulfideᵃ"*. Footnote a: *"An **aqueous solution (50 mL)** of hydrogen sulfide and
2-oxopropanal was reacted at **pH 5.0** in a laboratory autoclave."*
Column headers as printed: `2-oxopropanal (mmol) | H2S (mmol) | µgᵇ | %`.

| 2-oxopropanal (mmol) | H₂S (mmol) | **mercapto-2-propanone µg** | **%** |
|---|---|---|---|
| **1** | **1** | **1650** | **1.8** |
| **1** | **2** | **3600** | **4.0** |

**★ The only dose–response in the paper. [Z]: doubling H₂S gives 2.18× product on a mass basis
and 2.22× on a mol% basis — SUPER-LINEAR in H₂S.** Authors: > "it might be assumed that
**H₂S is the limiting factor** in this reaction."

### 2.8 Table 8 — FFT and MFT from hydroxyacetaldehyde + mercapto-2-propanone (pH series + the second dry point)
**Anchor: Table 8, p. 240.** Title: *"Amounts of FFT and MFT Generated from Hydroxyacetaldehyde
and Mercapto-2-propanoneᵃ"*. Footnote a: *"An **aqueous solution (50 mL)** of hydroxyacetaldehyde
(**1 mmol**) and mercapto-2-propanone (**1 mmol**) was reacted in a laboratory autoclave."*
Footnote c: *"A mixture of hydroxyacetaldehyde (1 mmol) and mercapto-2-propanone (1 mmol) was
**intimately mixed with silica gel (3.0 g containing 300 µL of phosphate buffer, pH 5.0)** and
was **dry-heated for 6 min at 180 °C**."*

| expt | pH | **FFT µg** | **FFT %** | **MFT µg** | **MFT %** |
|---|---|---|---|---|---|
| **1** | **3.0** | **26.1** | **0.02** | **15.5** | **0.01** |
| **2** | **5.0** | **40.5** | **0.04** | **268.1** | **0.23** |
| **3** | **7.0** | **58.2** | **0.05** | **311.5** | **0.27** |
| **4** | **5.0ᶜ (DRY, 180 °C / 6 min)** | **51.3** | **0.05** | **1553.9** | **1.39** |

**This is the second and last dry-heating measurement in the paper.** **[Z]** expt 4 vs expt 2
(the matched aqueous pH 5.0 control): **FFT 1.27× ↑, MFT 5.80× ↑.** Authors:
> "compared to an aqueous system (cf. experiments 2 and 4; Table 8), the **FFT formation
> increased slightly** during dry-heating. However, the **MFT generation was enhanced by a
> factor of nearly 6.**"

**★ Note this is a fed-intermediate system with no sugar and no cysteine, and it reproduces the
same qualitative dry-heating signature as the hexose rows of Table 2 (FFT flat-to-down, MFT
strongly up). Two independent systems, same direction.** That is the strongest thing the
water-content data can support.

### 2.9 Table 9 — mercapto-2-propanone from glucose + cysteine ⚠️ **NOT IN THE EXISTING INVENTORY**
**Anchor: Table 9, p. 240.** Title: *"Formation of Mercapto-2-propanone from Glucose and
Cysteineᵃ"*. Footnote a: *"A solution (**100 mL**) of **glucose (10.0 mmol)** and **cysteine
(3.3 mmol)** was reacted in a laboratory autoclave."*

| pH | **mercapto-2-propanone (µg)** |
|---|---|
| **3.0** | **11.4** |
| **5.0** | **59.6** |
| **7.0** | **26.5** |

> "The data (Table 9) confirmed that the mercapto-2-propanone is indeed formed when
> carbohydrates are reacted in the presence of cysteine, with its **formation going through a
> maximum at pH 5.0**."

**★ THIS TABLE IS A GAP.** I could not find these three values anywhere in
`k3_final_parameter_inventory.md`, and **Table 9 is absent from
`FIT_HOLDOUT_DECLARATION.md` lines 23–26**, which name only T1/T2 and T3/T4/T6/T7/T8/T10.
**Table 9 is the only in-situ measurement of a C₃ sulfur intermediate generated from a real
sugar/amino-acid pair, at three pH values** — it is the missing link between the fed-precursor
Table 8 chemistry and the actual glucose/cysteine system, and it shares the pH-5 maximum of the
hexose rows in Table 2. See the proposal below.

### 2.10 Table 10 — MFT from various precursor systems (the route comparison)
**Anchor: Table 10, p. 240.** Title: *"Amounts of MFT Liberated from Various Precursors or
Precursor Systemsᵃ"*. Footnote a: *"An **aqueous solution (50 mL)** of the precursors was reacted
at **pH 5.0** in a laboratory autoclave."* Column header: `precursor (1.0 mmol) each | MFT µgᵇ | %`.

| precursor (1.0 mmol each) | **MFT µg** | **MFT %** |
|---|---|---|
| **hydroxyacetaldehyde / mercapto-2-propanone** (C₂ + C₃) | **268.1** | **0.24** |
| NF / H₂S | **211.2** | **0.19** |
| NF / cysteine | **50.8** | **0.05** |
| thiamin | **8.2** | **0.01** |

**The C₂ + C₃ route is the strongest MFT route at 145 °C.** (Already ingested; the Cerny 2004
conflict at 95 °C is recorded at `k3_final_parameter_inventory.md:720`.)

⚠️ **INTERNAL INCONSISTENCY [Z]:** the hydroxyacetaldehyde/mercapto-2-propanone row appears in
**both Table 8 (expt 2) and Table 10** with the same µg (**268.1**) but **different mol %**:
**0.23 % in Table 8, 0.24 % in Table 10.** Same experiment, same mass, two percentages. A
rounding artefact of the mol% basis, not a second measurement — **do not ingest it twice.**

---

## 3. THREE DISCREPANCIES FOUND BY THE RECONCILIATION

### 3.1 ⚠️ Benchmark filename says 140 °C; the paper says 145 °C
`docs/reference/VALIDATION_CONTRACT.md:93` and `:158`, and
`k3_final_parameter_inventory.md:1209`, all reference a benchmark named
**`cys_ribose_140C_Hofmann1998`**. **The paper's aqueous condition is 145 °C / 20 min**
(Experimental, p. 236, quoted in §1.1), and the inventory's own condition strings correctly say
**145 °C** in five places (lines 215, 235, 304, 740, 755).

**So the benchmark's *name* says 140 and its *documented condition* says 145.** Either the name
is a legacy typo (most likely) or the benchmark JSON is actually parameterised at 140 °C, in
which case **every simulation against it has been run 5 K cold.** I have **not** opened or
touched `data/benchmarks/` (another agent is building there, per the brief). **Orchestrator
action required: check the `temperature` field inside
`data/benchmarks/cys_ribose_140C_Hofmann1998.json`.** If it reads 140, that is a live
correctness bug, not a cosmetic one; if it reads 145, the filename should be renamed or a note
added.

### 3.2 ⚠️ Table 9 is uningested (see §2.9)
Three measured values, SIDA tier, absent from both the inventory and the declaration.

### 3.3 ⚠️ The declaration's "only a_w extrapolation test" wording overstates the data (see §1.2)
The six dry-heated cells vary a_w, temperature, time and buffer molarity together.

**Everything else checks out.** All of Tables 1–8 and 10 as recorded in
`k3_final_parameter_inventory.md` §A.3.1–A.3.2 match the PDF exactly, including the ±10 %
triplicate footnote, the 100 mL / 50 mL volume split, the 1:3 cys:sugar stoichiometry, and the
0.5 mol/L buffer. The `342 / 200 ppb` values the inventory marks **REFUSE** as fabricated do
indeed appear **nowhere** in this paper — **confirmed by full-text search.**

---

## NEW-PARAMETER TABLE (consolidated)

Only genuinely **new** rows are listed — the ~60 values already correctly held in
`k3_final_parameter_inventory.md` §A.3 are re-transcribed in §2 above for verification but are
**not** new parameters and must not be re-ingested.

| parameter | value | units (as printed) | conditions | anchor | provenance |
|---|---|---|---|---|---|
| **mercapto-2-propanone from glucose + cysteine, pH 3.0** | **11.4** | µg | 100 mL, glucose 10.0 mmol + cysteine 3.3 mmol, 0.5 M phosphate, **145 °C / 20 min** | **Table 9, p. 240** | **[M] — NEW** |
| **mercapto-2-propanone from glucose + cysteine, pH 5.0** | **59.6** | µg | ditto | **Table 9, p. 240** | **[M] — NEW** |
| **mercapto-2-propanone from glucose + cysteine, pH 7.0** | **26.5** | µg | ditto | **Table 9, p. 240** | **[M] — NEW** |
| mercapto-2-propanone pH profile shape | **maximum at pH 5.0** (5.23× pH 3.0, 2.25× pH 7.0) | — | ditto | Table 9 + p. 240 text | **[M]/[Z] — NEW** |
| dry/aqueous ratio at pH 5.0, **ribose** | **FFT 8.03× ↑, MFT 1.27× ↑** | — | T2 col 5.0ᶜ ÷ col 5.0 | Table 2, p. 237 | **[Z]** |
| dry/aqueous ratio at pH 5.0, **glucose** | **FFT 0.50× ↓, MFT 2.21× ↑** | — | ditto | Table 2, p. 237 | **[Z]** |
| dry/aqueous ratio at pH 5.0, **rhamnose** | **FFT 0.50× ↓, MFT 3.88× ↑** | — | ditto | Table 2, p. 237 | **[Z]** |
| dry/aqueous ratio, **fed C₂+C₃ system** | **FFT 1.27× ↑, MFT 5.80× ↑** | — | T8 expt 4 ÷ expt 2 | Table 8, p. 240 | **[Z]** |
| H₂S dose–response exponent | **2.18× product for 2× H₂S** (mass); 2.22× (mol %) | — | 2-oxopropanal 1 mmol, pH 5.0, 145 °C | Table 7, p. 239 | **[Z]** |
| NF/MFT ratio, ribose / glucose / rhamnose | **2754× / 7.1× / 23.9×** | — | pH 5.0, 145 °C, aqueous | Table 5, p. 238 | **[Z]** |
| FA/FFT ratio, ribose / glucose / rhamnose | **5.58× / 2.46× / 7.75×** | — | ditto | Table 5, p. 238 | **[Z]** |
| dry-arm process definition | **silica gel 3.0 g + 300 µL of 2 mol/L phosphate; 180 °C; 6 min** | — | — | Experimental, p. 236 | **[M]** — the exact a_w-arm spec, needed to parameterise any dry-regime run |
| aqueous-arm process definition | **0.5 mol/L phosphate; 145 °C; 20 min; 200 mL Roth Type II autoclave** | — | — | Experimental, p. 236 | **[M]** |
| global error convention | **triplicates, "differed by not more than 10 %"** | % | every table | Table 1 footnote b, p. 237 | **[M]** |
| Table 8 expt 2 vs Table 10 row 1 mol % | **0.23 % vs 0.24 %** for the same 268.1 µg | mol % | — | Tables 8 and 10 | **[M]** — ⚠️ same experiment printed twice; **do not double-count** |

---

## PROPOSED FIT / HOLD-OUT ROLE — DRAFT FOR ORCHESTRATOR

> **Hofmann 1998 is ALREADY in `docs/reference/FIT_HOLDOUT_DECLARATION.md` (lines 23–26, 168),
> and that existing partition must be preserved** — `k3_final_parameter_inventory.md:1200`
> states the rule explicitly. **The proposals below concern ONLY the newly surfaced Table 9 and
> the two wording corrections.** No existing Hofmann row's role is proposed for change, and the
> declaration was not edited.

| dataset (specific rows) | proposed role | cut axis | rationale |
|---|---|---|---|
| **Table 9, pH 5.0 row** (mercapto-2-propanone 59.6 µg from glucose + cysteine) | **FIT** | pH | **Respects the existing Hofmann axis cut exactly** (`k3_final_parameter_inventory.md:1200–1202`: pH 5.0 rows FIT, pH 3/7 held out). This is the corpus's only *in-situ* C₃-sulfur-intermediate measurement from a real sugar/amino-acid pair; it pins the mercapto-2-propanone node that Table 8 and Table 10 currently reach only by feeding the intermediate |
| **Table 9, pH 3.0 and pH 7.0 rows** (11.4 and 26.5 µg) | **★ HOLD-OUT** | **pH** | Same cut axis as every other Hofmann row. These two are a **strong** hold-out: the pH-5 maximum (5.2× over pH 3, 2.3× over pH 7) is a *non-monotone* shape, which is the hardest thing for a pH term to reproduce, and the model has never seen a mercapto-2-propanone pH profile |
| **Table 2 dry-heated column, all six cells** | **HOLD-OUT** *(already declared, role unchanged)* | **regime (confounded)** | Role is right. **Recommend the rationale text be amended** from "the only a_w extrapolation test" to something like *"a confounded regime jump (a_w × 180 °C × 6 min × 2 M buffer); informative on DIRECTION only — a fail does not localise to a_w"*. See §1.2 |
| **Table 8 experiment 4** (dry, FFT 51.3 / MFT 1553.9 µg) | **HOLD-OUT** *(already declared under T8)* | regime | Same caveat. Its value is that it **replicates the Table 2 hexose direction in a sugar-free, cysteine-free fed system** — two independent systems agreeing on the sign is worth more than either alone |
| **The pentose/hexose SIGN REVERSAL under dry heating** (§2.2 [Z] table) | **qualitative ingest, HIGH VALUE** | — | Pentose FFT **+8×** while hexose FFT **−2×** under the identical regime change. A single scalar a_w or temperature multiplier **cannot** produce opposite signs for two feedstocks; reproducing this requires the two routes to be genuinely distinct in the network. **This is a structural constraint, and it is free** |
| **Table 8 expt 2 / Table 10 row 1 duplicate** | **★ neither — DEDUPLICATE** | — | Same experiment, 268.1 µg, printed as 0.23 % and 0.24 %. Ingest once |
| **Hofmann 1996b Tables 1 and 2** | **neither — DUPLICATES** *(already declared at line 142)* | — | **Confirmed independently by this reconciliation**: they are verbatim re-publications of jf9705983 Tables 2, 8 and 10. The existing flag is correct and must stay |
| **All of Tables 1–8, 10** | **roles UNCHANGED** | pH (existing cut) | Verified correct against the PDF. No amendment proposed |

### Actions the orchestrator must take that are outside this wave's write scope

1. **★ Check the temperature field inside `data/benchmarks/cys_ribose_140C_Hofmann1998.json`.**
   The paper is **145 °C**. If the JSON says 140, it is a live 5 K error affecting a benchmark
   the validation contract calls *"the last strict-ready free-amino-acid benchmark"*
   (`VALIDATION_CONTRACT.md:158`). I did not open `data/benchmarks/` per the brief's scope
   restriction.
2. **Amend the declaration's line-26 rationale** to drop "the only a_w extrapolation test"
   (§1.2). The role stays HOLD-OUT; only the claim about what it tests is wrong.
3. **Add Table 9** to the declaration with the pH cut proposed above.
4. **Record the negative result**: the "Hofmann 1998 temperature/water-content series" named as
   a top gap **does not exist in the paper**. The gap cannot be closed from this source and
   should be reopened against a different one (a paper with ≥3 water-activity levels at fixed
   T and t). Leaving it marked "closed" on the strength of six confounded cells would be the
   kind of circular certification the validation contract exists to prevent.
