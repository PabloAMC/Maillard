# Wave K — Benchmark CONTENT verification (not citation identity)

**Date:** 2026-08-27 · **Scope:** `data/benchmarks/*.json` (live panel), `data/benchmarks/external_validation/*.json` (hold-out).
**Method:** for each benchmark, the exact numeric claim was extracted from the JSON, then the *cited paper's own text/tables* were
retrieved (Europe PMC full-text XML, PMC-rendered PDFs, ACS figshare supporting information, publisher-indexed abstracts,
CrossRef/OpenAlex/Unpaywall) and compared. **Repo was not modified.** Every CONFIRMED / MISMATCH below carries a retrieved quote.

**Prior state.** `results/validation/citation_verification_ledger.md` states its own source of truth explicitly:
*"**Source of truth:** CrossRef REST API (`api.crossref.org`) plus CrossRef bibliographic title search"* and
*"**No numeric value anywhere in the repo was changed by this pass.**"* It verified DOI→paper identity only.
Wave K is the first pass to open the papers.

---

## Headline

**The Pratap-Singh suspicion is confirmed, and it is worse than a 4.4× scaling error.**
Both `PratapSingh2021` files misreport hexanal by ~4.3×, **and** both invent a hexanol value for a compound the paper
reports as **not detected**. Separately, **two of the four Li 2026 hold-out points are wrong table rows** — the
external-validation set, the one lane that is supposed to be untouched hold-out evidence, contains misattributed numbers.

| # | Finding | Severity |
|---|---|---|
| 1 | Pratap-Singh hexanal 260/380 ppb vs paper's **1138.00 / 1621.71 ppb** | FATAL (4.38× / 4.27× low) |
| 2 | Pratap-Singh hexanol 80/120 ppb vs paper's **n.d.** for both | FATAL (fabricated) |
| 3 | Li 2026 hold-out "nonanal 29.42" is the paper's **Decanal** row (nonanal = 74.37) | FATAL (wrong compound, hold-out) |
| 4 | Li 2026 hold-out "2-pentylfuran 221.5" is the paper's **Maltol** row (2-PF = 5625.80) | FATAL (wrong compound, 25× low, hold-out) |
| 5 | Acrylamide 62.62 ppb not locatable in the cited paper (130 °C bar ≈ 150 µg/kg) | FATAL-probable |
| 6 | Resconi/Hernandez furfural 1040 ppb silently excludes 1 of the 3 PBMAs | RECONCILABLE but undeclared |
| 7 | Cerny2008 file's stored citation string names a *different* paper than its own DOI | METADATA |

---

## (a) Verdict table

| File | Claim | Source (as cited) | Verdict | Evidence |
|---|---|---|---|---|
| `pea_isolate_40C_PratapSingh2021.json` | hexanal **260** ppb, pea isolate | 10.3390/molecules26134104 | **MISMATCH — fatal** | Table 1: `Hexanal (in ppb): Pea 1138.00 ± 297.30` |
| `pea_isolate_40C_PratapSingh2021.json` | hexanol **80** ppb | same | **MISMATCH — fatal (fabricated)** | Table 1: `Hexanol (in ppb): Pea n.d.`; text: *"pea proteins contained no alcohol compounds"* |
| `pea_isolate_40C_PratapSingh2021.json` | 2-pentylfuran **638** ppb (±8%) | same | **CONFIRMED** | *"2-pentyl furan was found in large quantities in soy and pea proteins (2492 ± 199 and 638 ± 49 ppb equivalents of hexanal, respectively…)"* (49/638 = 7.7% ≈ the file's 8%) |
| `soy_isolate_40C_PratapSingh2021.json` | hexanal **380** ppb, soy isolate | same | **MISMATCH — fatal** | Table 1: `Hexanal (in ppb): Soy 1621.71 ± 159.69` |
| `soy_isolate_40C_PratapSingh2021.json` | hexanol **120** ppb | same | **MISMATCH — fatal (fabricated)** | Table 1: `Hexanol (in ppb): Soy n.d.`; text: *"Total alcohol concentration of soy protein was only 40 ± 9 equivalents of hexanol in form of 1-octen-3-ol"* |
| `soy_isolate_40C_PratapSingh2021.json` | 2-pentylfuran **2492** ppb (±8%) | same | **CONFIRMED** | same sentence; 199/2492 = 8.0% |
| `external_validation_li_2026_spi_wg_hme_control.json` | hexanal **605.6** ppb | 10.3390/foods15050912 | **CONFIRMED** | Table 2, HMPE-0 min: `Hexanal 780 5 605.64 ± 6.50 a` |
| `…li_2026…` | 1-hexanol **20.04** ppb | same | **CONFIRMED** | Table 2: `1-Hexanol 868 5.6 20.04 ± 0.66 b` |
| `…li_2026…` | nonanal **29.42** ppb | same | **MISMATCH — fatal (wrong row)** | Table 2: `Nonanal 1102 1.1 74.37 ± 0.11 a`; **`Decanal 1205 3 29.42 ± 1.08 a`** |
| `…li_2026…` | 2-pentylfuran **221.5** ppb | same | **MISMATCH — fatal (wrong row, 25×)** | Table 2: `2-Pentylfuran 994 5.8 5625.80 ± 63.75 c`; **`Maltol 1118 1.24 221.51 ± 5.74 d`** |
| `…li_2026…` | matrix: SPI:WG 6:4 dry basis, 160 °C | same | **CONFIRMED (conditions)** | *"soy protein isolate and wheat gluten hydrolysates at a dry basis ratio of 6:4… temperature profile set at 30, 90, 120, 140, 150, and 160 °C"* — but see caveat on `water_activity` below |
| `acrylamide_spi_extrusion_130C_ACSRef3.json` | acrylamide **62.62** ppb @130 °C | 10.3390/ijms25168668 (Ma 2024) | **MISMATCH — probably fatal** | Fig. 2D, 130 °C bar ≈ **150 µg/kg** (highest, letter `a`); text: *"The lowest and highest acrylamide content were found in PBMAs extruded by 20% moisture and 130 °C, respectively"*. No 62.62 anywhere; no table of acrylamide values exists in the paper |
| `cml_cel_commercial_pbma_Foods2023.json` | CML range **16.46–47.61** mg/kg; CEL **25.21–86.23** mg/kg | 10.3390/foods12101967 (Fu 2023) | **CONFIRMED (range)** | Table 1: CML min `O 16.46 ± 0.67`, max `F 47.61 ± 3.82`; CEL min `L 25.21 ± 3.12`, max `G 86.23 ± 2.11` |
| `cml_cel_commercial_pbma_Foods2023.json` | scored values CML **32 000** / CEL **55 000** ppb | same | **PARTIAL (repo-derived)** | Range midpoints: (16.46+47.61)/2 = 32.04 ✓; (25.21+86.23)/2 = 55.72 ≈ 55. Not reported values; sample means are CML 32.6, **CEL 44.3** mg/kg, so "55" is a range midpoint, not a central tendency |
| `resconi_2023_pbma_beef_identity_benchmark.json` | furfural **1040** ppb, "midpoint PBMA" | 10.3390/molecules28073151 | **MISMATCH — reconcilable, undeclared** | Table 3 (ng/g): `Furfural … Beyond Meat 987.41 a, Impossible Burger 64.71 b, Third Retail Brand 1093.54 a`. (987.41+1093.54)/2 = **1040.5**. The third PBMA (64.71) is silently dropped; the true 3-product mean is 715.2 |
| `external_validation_bi_2020_raw_pea_hexanal.json` | hexanal **1260** µg/kg raw pea | 10.1021/acs.jafc.9b07711 | **PARTIAL** | Abstract, verbatim: *"3-methylbutanoic acid (OAV = 382) and hexanal (**OAV = 280**) significantly contributed to the aroma of peas."* OAV 280 CONFIRMED; the 1260 concentration itself not retrieved (see below) |
| `external_validation_bi_2020_roasted_pea_hexanal.json` | hexanal **324** µg/kg roasted pea | same | **UNVERIFIABLE** | Roasted-pea hexanal appears in neither the abstract nor the SI; ACS full text 403 |
| `external_validation_bi_2020_*` | DOI points at real Bi 2020 pea paper (EGCG mis-repair reversed) | — | **CONFIRMED** | CrossRef `10.1021/acs.jafc.9b07711` = Bi, Xu, Luo, Lao, Pang, Shen, *"Characterization of Key Aroma Compounds in Raw and Roasted Peas (Pisum sativum L.)"*, JAFC 68:2718-2727. Both files and the intake row `acs_2020_raw_pea_hexanal_baseline` carry this DOI. **Reversal is correctly applied.** |
| `external_validation_liu_2023_ppi_offnote_baseline.json` | hexanal band **15–180**, nonanal **5–50** µg/kg | Liu 2023 NCSU thesis / Liu, Cadwallader & Drake 2023 Food Chem 406:134998 | **UNVERIFIABLE** | Thesis repository behind bot-detection; journal version closed (Unpaywall `oa_status: closed`, no repository copy) |
| `pea_isolate_uht_140C_Trikusuma2019.json` | hexanal **782**, 2-PF **163**, nonanal **24** µg/L | 10.1016/j.foodchem.2019.126082 | **UNVERIFIABLE** | Closed, no OA copy anywhere (Unpaywall/OpenAlex both `closed`, `has_repository_copy: false`); abstract confirms *"twenty-one aroma compounds were identified and further quantified using dynamic headspace-GC/MS/MS"* — i.e. a quantitative table of the right shape exists, contents unseen |
| `cys_ribose_140C_Hofmann1998.json` | MFT **342**, FFT **200** ppb (ribose+Cys, 140 °C, pH 5, 30 min) | 10.1021/jf9705983 | **UNVERIFIABLE — with a warning** | Closed. Abstract reports yields in **mol %** (*"the highest yields for MFT (1.4 mol %)"*), not ppb. A ppb value requires an undocumented conversion; the repo file records no table/figure locus and no derivation |
| `thiamine_cys_glucose_120C_Bolton1994.json` | MFT **13** ppb (389–489 ng/33.3 g) | 10.1021/bk-1994-0543.ch022 | **UNVERIFIABLE (already self-declared)** | CrossRef confirms chapter identity + pages 270-278; no abstract is served by CrossRef or Semantic Scholar (*"fields have been elided by the publisher"*). The file already says composition/values come from an indexed abstract and that pH and molarities are assumed |
| `thiamine_cys_xylose_145C_Cerny2008.json` | MFT **2.47** ppb | 10.1021/jf801762c | **UNVERIFIABLE + METADATA-MISMATCH** | CrossRef for that DOI = *"Identification of 5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine, Cysteine, and Xylose"*, JAFC 56:**10679-10682**. The file's own `source_metadata.citation` says *"Formation of Cysteine-S-Conjugates…"*, 56(23), **10237-10241** — a different title and page range from its own DOI. The file already flags the value as unverified |
| `furosine_extrusion_crossover_140C_RamirezJimenez2000.json` | furosine **8.7 mg/100 g protein** → 17 400 ppb | 10.1021/jf9907687 | **UNVERIFIABLE (value) / MISMATCH (matrix, self-declared)** | CrossRef confirms *"Browning Indicators in **Bread**"*, JAFC 48:4176-4181. Closed. The file's own audit note concedes the matrix is toasted bread, not a legume extrudate |
| `pea/soy_isolate_ribose_cysteine_…_Internal2026.json`, `…ProtocolPilot2026.json` | 6–7 volatiles at 0.003–0.28 ppb | none (self-declared synthetic) | **N/A — correctly labelled** | `source_metadata.origin: internal_reproducibility_candidate` / `synthetic_diagnostic`; *"frozen output of the current model"*. No literature claim to check |

---

## (b) Ranked fatal findings

### 1. Pratap-Singh hexanal is ~4.3× low — the matrix reference lane is anchored low (FATAL)

Molecules 2021, 26, 4104 (PMC8271896), **Table 1**, caption *"Fully quantified concentration of the standard compound
(hexanal, 2-nonanone, and hexanol) in plant protein extracts"*:

```
Concentrations        Pea                Brown Rice              Soy
Hexanal (in ppb)      1138.00 ± 297.30   22,590.24 ± 1643.70     1621.71 ± 159.69
2-Nonanone (in ppb)      6.382 ± 0.62        94.02 ± 12.38       n.d.
Hexanol (in ppb)      n.d.                  102.04 ± 9.30        n.d.
```

Repo: pea 260, soy 380. Ratios **1138.00/260 = 4.38×**, **1621.71/380 = 4.27×**.
This is *not* a unit error — the table header says "in ppb", the same unit the repo uses, and the 2-pentylfuran values in the
same paper were transcribed exactly right (638, 2492) with the paper's own ±% preserved. The hexanal numbers are the only
ones in these two files that required opening Table 1 rather than reading the running text.

**Probable mechanism (hypothesis, clearly labelled as such).** The paper's *running text* gives total aldehydes as
*"1359 ± 321 [pea] … 1998 ± 201 [soy] ppb equivalents of hexanal"*. The repo values are ≈19% of those totals in both cases
(260/1359 = 0.1913; 380/1998 = 0.1902 — the same fraction to three digits). A single share-of-total factor applied to the
text numbers reproduces both repo values; reading two different proteins off Figure 2's "% aldehyde composition" would not
give the same percentage twice. I did not retrieve Figure 2's numbers, so this remains a hypothesis — but the two verified
mismatched values are established regardless.

**Implication.** `results/validation/projection_constant_refit.md` shows the projection constant was fitted to these exact
numbers (predicted 260.6 vs 260, 379.9 vs 380 — residuals ~0.001 dex). The matrix-only lane is therefore **calibrated to
values ~4.3× below the source**, and `matrix_sigma_residual_derivation.json` derives the matrix σ from residuals against
them. Correcting the anchors moves the lipid-aldehyde calibration by ~0.63 dex in the *opposite* direction to the S27
over-prediction story: the model's matrix hexanal is not 36× too high against truth, it is 36× too high against a reference
that is itself 4.3× too low. The "matrix path over-predicts" framing in `memory/matrix-path-accuracy-plan.md` needs to be
re-derived from scratch.

### 2. Pratap-Singh hexanol 80/120 ppb is fabricated content behind a valid citation (FATAL)

Table 1 gives `Hexanol (in ppb): Pea n.d. … Soy n.d.` The results text is explicit:

> *"Total alcohol concentration of brown rice protein was 4852 ± 458 equivalents of hexanol. Total alcohol concentration of
> soy protein was only **40 ± 9** equivalents of hexanol **in form of 1-octen-3-ol**, while **pea proteins contained no
> alcohol compounds**."*

The repo's soy hexanol (120 ppb) is **3× the paper's entire soy alcohol fraction**, and its pea hexanol (80 ppb) is a
compound class the paper says is absent from pea. Both files also encode hexanol as `expected_rank 3` in their
`matrix_ranking_contract`, so a non-existent measurement is load-bearing in the ranking contract.

*Caveat, stated for honesty:* the paper is internally inconsistent — its Discussion says *"the alcohol profile derived from
soy protein… was higher in 1-hexanol, along with a small proportion of 1-heptanol"*, contradicting its own Results and
Table 1. Even on the most generous reading, soy 1-hexanol cannot exceed the 40 ± 9 ppb total, so 120 ppb is unsupported;
pea 80 ppb is unsupported on every reading. The paper's Table S1 (semi-quantified 76-compound profile) is hosted only on
mdpi.com and could not be retrieved — a human should check whether S1 lists a soy 1-hexanol figure, but it cannot exceed 40.

### 3–4. Li 2026 hold-out: two values are the wrong table rows (FATAL — this is hold-out evidence)

Foods 2026, 15, 912 (PMC12984281), **Table 2**, *"Relative contents (µg/kg) of volatile odour compounds in HMPEs"*,
column HMPE-0 min (the control the repo claims):

```
Nonanal        1102  1.1    74.37 ± 0.11 a
Decanal        1205  3      29.42 ± 1.08 a      <-- repo calls this "nonanal"
2-Pentylfuran   994  5.8  5625.80 ± 63.75 c
Maltol         1118  1.24   221.51 ± 5.74 d     <-- repo calls this "2-pentylfuran"
```

Repo `nonanal: 29.42` and `2-pentylfuran: 221.5`. Both are the **row immediately below** the correct one in the table's own
ordering (Nonanal→Decanal; 2-Pentylfuran→Maltol) — an off-by-one transcription, twice, in the same table. The 2-pentylfuran
error is **25.4× low**. Hexanal (605.64) and 1-hexanol (20.04) in the same column were transcribed correctly.

These four points carry `"value_provenance": "reported_point_value"` and
`"value_provenance_note": "A single concentration reported by the source. Scoring this against a prediction is a genuine test."`
Two of the four are not what the source reports. Any hold-out score computed on this bundle is invalid for those two markers.

Secondary condition mismatch in the same file: the paper runs high-moisture extrusion at *"moisture content maintained at
approximately 57%"*; the benchmark encodes `water_activity: 0.35`, which is a low-moisture-extrusion value. (Whether this is
a modelling knob or a claim about the source is a repo-internal question, but as a description of the source it is wrong.)

### 5. Acrylamide 62.62 ppb cannot be located in the cited paper (FATAL-probable)

Ma et al. 2024, IJMS 25:8668 reports acrylamide **only in Figure 2** (bar charts; no acrylamide table exists — the paper's
single table is water/protein/TAA). Reading the PMC-rendered Figure 2D at 260 dpi: the 130 °C bar sits at ≈**150 µg/kg**
(the highest bar, letter `a`), 150 °C ≈ 120, 170 °C ≈ 82, unextruded control ≈ 38. The lowest extruded bar anywhere in the
figure is 20% moisture at ≈68 µg/kg. **62.62 appears nowhere**, and the paper's own 130 °C point is ~2.4× the benchmark value.

I also tested whether 62.62 is a derivation from the companion paper (Fu 2023, Table 1, acrylamide 31.81–186.70 µg/kg):
mean 68.55, median 63.07, geometric mean 61.39, range midpoint 109.26 — **none is 62.62**.

Note the file's own id: `acrylamide_spi_extrusion_130C_**ACSRef3**` against an **MDPI IJMS** DOI. The value plausibly
predates the current DOI and was never re-checked against it — the exact failure mode this wave was chartered to find.
The file already declares itself `expected_to_fail` for an unrelated reason (SME desaturation), which would mask this.

### 6. Resconi/Hernandez furfural: an undeclared cherry-pick (reconcilable, but fix the label)

Table 3 of Molecules 28:3151 gives furfural for **three** PBMAs: 987.41, 64.71, 1093.54 ng/g. The benchmark's 1040.0 is
exactly `(987.41 + 1093.54)/2`, i.e. the mean of the two high brands with Impossible Burger dropped. The file describes it as
*"Midpoint PBMA furfural value extracted from the absolute ng/g cooked-product table"* — which reads as a midpoint over the
PBMA set. Over the actual PBMA set the mean is 715.2 (1.45× lower) and the range is 64.71–1093.54 (17×). Units are fine
(ng/g = ppb) and these are least-squares means over cooked commercial product. Not fabrication; an undeclared subset choice
presented as a midpoint, inside a benchmark whose tolerance is ±5%.

### 7. Cerny2008 file's citation string names a different paper than its own DOI (metadata)

`source_metadata.citation` = *"Formation of Cysteine-S-Conjugates in the Maillard Reaction of Cysteine, Xylose and Thiamine.
JAFC 2008, 56 (23), 10237-10241."* CrossRef for the file's DOI `10.1021/jf801762c` = *"Identification of
5-Hydroxy-3-mercapto-2-pentanone in the Maillard Reaction of Thiamine, Cysteine, and Xylose"*, JAFC 56, **10679-10682**.
Different title, different pages. The 2026-08-26 DOI repair fixed the identifier but left the human-readable citation
describing the old (or a third) paper. The file's `VALUES_NEED_RE_EXAMINATION` warning already covers the numeric risk.

---

## (c) The Pratap-Singh answer — and the route that worked

**Answer: the higher numbers are right.** Pea hexanal **1138.00 ± 297.30 ppb**, soy hexanal **1621.71 ± 159.69 ppb**,
Table 1 of Molecules 2021, 26, 4104. The repo is low by 4.38× (pea) and 4.27× (soy). It is **not** a units artefact:
the table header reads "Hexanal (in ppb)" — the same unit the repo declares — and the paper's method section confirms the
matrix basis (protein isolate powder dispersed for HS-SPME, standards spiked into the pea/soy protein sample itself:
*"The pea protein sample was spiked with 1, 2, 3, 4, 5 µL of stock standard"*). No wet/dry or µg/kg-vs-ng/g reconciliation
is available or needed.

**Routes tried, in order:**

| Route | Result |
|---|---|
| `WebFetch https://www.mdpi.com/1420-3049/26/13/4104` | HTTP 403 |
| `curl` MDPI figure/PDF paths with a browser UA | HTTP 403 (MDPI blocks this environment entirely) |
| Europe PMC search API on the DOI | **Hit** — `PMC8271896`, `isOpenAccess: Y` |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8271896/fullTextXML` | **80 KB of JATS XML including full Table 1 — this is the route that worked** |
| WebSearch (independent corroboration) | Returned the 2-pentylfuran sentence verbatim, confirming the XML matches the published text |

MDPI deposits into PMC; the Europe PMC `fullTextXML` endpoint returns the tables as markup and bypasses the MDPI block
entirely. **Same route worked for every other MDPI benchmark** (Li 2026 → PMC12984281, Fu 2023 → PMC10217484,
Ma 2024 → PMC11354377, Hernandez 2023 → PMC10096055). Recommend it as the repo's standard MDPI verification path.
For figures, `https://europepmc.org/articles/<PMCID>?pdf=render` → `pdftoppm` → visual read works (used for Ma 2024 Fig. 2).

---

## (d) Still UNVERIFIABLE — what a human with library access should pull

Ranked by consequence:

1. **Trikusuma, Paravisini & Peterson (2020), Food Chem 312:126082** — the DHS-GC-MS-MS quantitative table.
   Check hexanal 782, 2-pentylfuran 163, nonanal 24, methional 3.1, 2,5-DMP 2.29 µg/L. **This is the only remaining
   heated-matrix quantitative anchor and it is now the last unchecked pillar of the matrix lane.**
   Tried: Unpaywall (`closed`, no repository copy), OpenAlex (2 closed locations), Europe PMC (no PMC ID),
   ScienceDirect (paywall), OhioLINK ETD for Trikusuma's OSU thesis carrying the same table (search endpoints returned
   only the shell page; no accession recoverable), fatcat/scholar.archive.org (API timed out).
2. **Bi et al. (2020), JAFC 68:2718-2727** — the raw *and roasted* pea hexanal concentration rows.
   Confirmed from the abstract: *"hexanal (OAV = 280)"*. Not confirmed: that the paper prints 1260 µg/kg, and nothing at all
   on the roasted point (OAV 72 / 324 µg/kg). **New evidence on the repo's own open question:** the intake registry's
   companion value, 3-methylbutanoic acid 4580 µg/kg, divides by the published OAV 382 to give exactly 12.0 µg/kg — a
   standard literature ODT for that acid, just as 4.5 is for hexanal. So *both* compounds divide exactly by their
   conventional thresholds. That is equally consistent with (i) the paper printing concentrations that Bi computed OAVs
   from using those thresholds, and (ii) the repo back-computing concentrations from OAVs. The ambiguity the Wave I note
   raised is real and cannot be resolved from outside the paper. Tried: ACS 403, figshare SI (retrieved — 3 pages,
   figure captions only, no data tables), Unpaywall closed.
3. **Hofmann & Schieberle (1998), JAFC 46:235-241, Tables** — the ribose/cysteine 140 °C point.
   The abstract reports **mol %**, so 342/200 ppb must come from a conversion the repo does not document. Ask for the
   ribose+cysteine row and the units of the tabulated yields. Note the repo already knows this citation was previously
   used to back a claim it could not support (`slr_incorporation_matrix.json`: *"I could NOT confirm 0.021% or a 2.3x
   cysteine factor in it… re-extract the yields from the primary tables"*). This benchmark carries the tightest contract
   in the panel (1.45× / 0.09 dex) on a number nobody has read.
   **RESOLVED 2026-08-27 (Waves S2b/S2c), and the answer is worse than this item assumed.**
   The ILL request below is still worth making, but not to validate 342/200: those two numbers
   are not from the paper at all. Wave S2b traced them to
   `data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
   table committed in `c7efbbc`, the same commit that created the benchmark JSON — whose row
   gives MFT `~0.02–0.05` mol % and FFT `~0.01–0.03` mol %. On the file's declared, unattested
   10 mM basis with MW 114.17: `0.0300 mol % → 342.5 → 342 ppb`, and the FFT band's geometric
   mean `0.017321 mol % → 197.8 → 200 ppb`. Both are interior points of two invented,
   overlapping bands (~90% confidence, arithmetic exact). Wave S2c retired the 1.45× / 0.09 dex
   contract, demoted the tier `PRIMARY → REFERENCE`, marked both values `no_verifiable_source`
   and reverted the barrier that had been fitted to them. **The sulfur branch has zero absolute
   literature anchors.** The ILL ask is now a *rebuild* ask: get the paper's own yields table in
   native mol %, with the Experimental section's volume and molar amounts. Full pack:
   `tasks/audit_remediation.md` "## Wave S2b" §(f).
4. **Liu, Cadwallader & Drake (2023), Food Chem 406:134998** (or the NCSU thesis) — the commercial-PPI hexanal 15–180 and
   nonanal 5–50 µg/kg bands. Tried: NCSU repository (Anubis bot-detection challenge), Unpaywall closed.
5. **Bolton et al. (1994), ACS Symp. Ser. 543, 270-278** — the 389–489 ng/33.3 g figure and the actual molarities.
   No abstract is served by CrossRef or Semantic Scholar (publisher-elided). The file's own caveats are accurate.
6. **Ramírez-Jiménez et al. (2000), JAFC 48:4176-4181** — the 8.7 mg/100 g protein furosine peak, and whether it is a
   140 °C extrusion point at all (the paper is titled *"Browning Indicators in Bread"*).
7. **Pratap-Singh Table S1** (mdpi.com/article/10.3390/molecules26134104/s1, not in the PMC deposit) — only to close out
   whether any soy 1-hexanol number exists there. Cannot exceed 40 ± 9 ppb.

---

## (e) Confidence per verdict

| Verdict | Confidence | Why |
|---|---|---|
| Pratap-Singh hexanal 4.3× low | **~99%** | Publisher-deposited JATS table read directly; independently corroborated by the 2-pentylfuran sentence in the same file matching the repo to the digit |
| Pratap-Singh hexanol fabricated | **~97%** | Table says n.d.; results text says pea has no alcohols and soy total is 40 ppb. Residual 3%: the paper's self-contradictory discussion sentence + unread Table S1 |
| The ~19%-of-total-aldehydes mechanism | **~65%** | Two data points agreeing to three digits is suggestive, not proof; Figure 2 not read |
| Li 2026 nonanal = Decanal row | **~97%** | Exact value match on the adjacent row; correct nonanal value present in the same column |
| Li 2026 2-pentylfuran = Maltol row | **~97%** | Same |
| Li 2026 hexanal / 1-hexanol correct | **~99%** | Exact match including the paper's ± |
| Acrylamide 62.62 unsupported | **~85%** | Confident the paper has no acrylamide table and that 130 °C ≈ 150 µg/kg; the 15% is figure read-off error and the unread supplementary (`…/ijms25168668/s1`, MDPI-hosted, unreachable) |
| Fu 2023 CML/CEL ranges | **~99%** | Table 1 read directly |
| Hernandez furfural = 2-of-3 mean | **~95%** | Arithmetic is exact (1040.475); the only judgement is whether the omission was intentional |
| Bi 2020 OAV 280 | **~99%** | Publisher abstract quoted verbatim |
| Bi 2020 concentration 1260 | **~50/50** | Genuinely undecidable from outside the paper — see (d)2 |
| Cerny2008 citation-string mismatch | **~99%** | CrossRef title and page range both differ from the stored string |
| All UNVERIFIABLE rows | n/a | No full text retrieved; nothing asserted about their contents |

**Not checked by this wave:** `data/lit/*.json` payloads (flavor_reference_payloads, computational_priors,
slr_incorporation_matrix, benchmark_intake_registry) beyond the provenance chains feeding the benchmarks above.
Given a 4-fatal-finding hit rate on ~24 checkable benchmark claims, the ~200 remaining `data/lit` anchors should be
assumed to carry a similar content-error rate until someone opens them.
