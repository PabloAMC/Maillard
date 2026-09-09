# Frankel, Neff & Selke 1981 — EXTRACTION (pure hydroperoxides of methyl OLEATE, LINOLEATE and LINOLENATE, from both autoxidation and photosensitized oxidation, thermolysed neat in a GC injector port at 210 C; four tables, ~60 quantified product shares, plus an isomer-composition table before and after pyrolysis)

### BOTH REFUSED COMPOUNDS ARE PRINTED IN THIS PAPER: **nonanal = 15 % (autoxidised oleate) and 10 % (photosensitized oleate)** of total volatile peak area, and **2-pentylfuran = 2.4 % (autoxidised linoleate) and 0.6 % (photosensitized linoleate)** — the alkylfuran the engine refuses, measured, in the same laboratory and by the same injector-port method as Frankel & Gardner 1989, and with hexanal in the same column so the ratio is denominator-free: **2-pentylfuran / hexanal = 0.16 (mine)**.

**Source on disk:** `data/articles/frankel1981.pdf` (7 pp., *Lipids* **16** (5), 279-285, 1981). The
`pdftotext -layout` text layer is an OCR layer whose prose is glyph-spaced and which substitutes
"S" for "5" in several table cells ("1.S", "0.S", "S.0") and "tel %" for "rel %". **All four tables
were re-read off 300-dpi rasters of printed pages 280, 281, 282 and 283**
(`scratchpad/img/fr81t1.png`, `fr81t2.png`, `fr81t3.png`, `fr81t4.png`) and every cell below matches
the raster. Scheme I (the alkoxy-radical scission scheme), Scheme II (the hydroperoxy cyclic
peroxide scheme) and the two inline mechanism sketches are **FIGURE-ONLY**; no number is read off
them. **There is no figure with data in this paper** — the entire quantitative content is Tables
I-IV.

## 0. Identity

| field | value |
|---|---|
| Title | "Analysis of Autoxidized Fats by Gas Chromatography-Mass Spectrometry: VII. Volatile Thermal Decomposition Products of Pure Hydroperoxides from Autoxidized and Photosensitized Oxidized Methyl Oleate, Linoleate and Linolenate" |
| Authors | **Edwin N. Frankel, William E. Neff and Edward Selke** — Northern Regional Research Center, Agricultural Research, Science and Education Administration, U.S. Department of Agriculture, Peoria, IL 61604. **Three authors, the same three as Selke 1978 in a different order** |
| Venue | *Lipids* **16** (5), 279-285 (1981). Received November 24, 1980 |
| **DOI** | **NO DOI IS PRINTED IN THE PDF.** No DOI, no CrossRef stamp, no copyright line anywhere in the seven pages. Cite by volume/page only |
| Naming | "rel %" = relative percent of the **total volatile peak area of that injection**, including the "Unidentified peaks" row — every one of the six product columns sums to 99.5-100.9 (mine, checked). "Auto" = free-radical autoxidation in air; "Photo" = photosensitized oxidation in O2, methanol solution, methylene blue. "Origin" = the hydroperoxide isomer the product is assigned to **by Scheme I**, not by measurement; "?" means the authors could not assign one |
| Its own lineage | This is **the full paper Selke 1978 promised** ("the implication of these results will be discussed in the full paper that will deal also with the decomposition of linoleate and linolenate hydroperoxides"). Its ref (21) is Selke, Frankel & Neff, *Lipids* 13:511 (1978), and **Table II's entire autoxidation column is footnoted "Data from ref. 21"** — i.e. it is Selke 1978's Table I reprinted, not a new measurement. See flag 2 |
| Repo status before this dossier | **Not cited anywhere in `src/kinetic_core/`.** `parameters_lipid.PROHIBITED_DERIVATIONS` refuses an alkylfuran branch fraction and refuses importing a propanal (linolenate) number; both refusals were correct about the corpus as constituted |

## 1. Why it matters

**This paper contains, printed, both of the two numbers the scoring panel refuses.**

*Nonanal.* `src/kinetic_core/lipid.py` refuses `NONANAL` on the ground that "the oleate -> nonanal
branch fraction is measured NOWHERE in the fit corpus". Table II gives it twice: **15 % of total
volatile peak area from autoxidised methyl oleate hydroperoxides (8-+9-+10-+11-OOH), and 10 % from
photosensitized-oxidised methyl oleate hydroperoxides (9-+10-OOH)**, both at a 210 C injector port,
with the isomer composition of each feed printed in Table I. The first of those two is the same
single measurement already carried in `selke1978_extraction.md` — it is republished here, not
replicated. **The photosensitized column is genuinely independent** and is the only second oleate
determination in the corpus.

*2-Pentylfuran.* The panel refuses it because "2-pentylfuran is NOT in Frankel 1989's six-product
slate and no branch fraction for the linoleate -> alkylfuran route is measured anywhere in the fit
corpus". Both clauses were true. Table III gives **2-pentylfuran at 2.4 % of total volatile peak
area from autoxidised methyl linoleate hydroperoxides and 0.6 % from photosensitized ones** — and,
crucially, **hexanal is in the same column at 15 % and 17 %**, so the repository can take the
**denominator-free ratio 2-pentylfuran / hexanal = 2.4/15 = 0.16 (autoxidation) and 0.6/17 = 0.035
(photosensitized)** (both mine) and hang the alkylfuran off a product it already models. That is
the single most directly usable number in this paper.

**But the alkylfuran's parent is unknown, and the paper says so in three places.** Table III's
Origin column gives 2-pentylfuran "**?**". The Results say "**2-Pentylfuran is a unique product of
autoxidation hydroperoxide, but its origin is not well established.**" The Discussion says
"**Although their origin is not clear** … 2-pentylfuran [is speculated to come] **from a
10-hydroperoxide intermediate (27)**", citing Chang et al. 1966. So a *branch fraction* exists in
the sense of "share of the volatile slate from an autoxidised linoleate hydroperoxide pool"; a
*mechanistic branch* from a named isomer does not. Section 4 keeps those apart, and flag 6 records
that **the paper's own Table I argues against the 10-hydroperoxide origin it speculates**.

**Relation to `frankel1989_extraction.md`, stated explicitly as required.**

| axis | Frankel, Neff & Selke 1981 (this paper) | Frankel & Gardner 1989 |
|---|---|---|
| laboratory | Northern Regional Research Center, USDA, Peoria IL | **the same** |
| shared author | **E. N. Frankel** (first author) | **E. N. Frankel** (first author) |
| method | injector-port thermolysis, GC-MS | **the same** injector-port thermolysis, GC |
| **injector temperature** | **210 C** | **180 C** |
| sample | **neat, 4-10 µL** | hexane solution, 1 µL of 200 µL |
| trapping / column start | **25 C**, programmed to 275 C at 2 C/min | **-65 C** cryo-trap, then 5 C/min to 260 C |
| GC column | **packed**: glass 14 ft x 4 mm i.d., 10 % OV-101 on Chromosorb G | **capillary**: DB-5, 60 m x 0.315 mm, 1 µm film |
| substrates | oleate **and** linoleate **and** linolenate; autoxidation **and** photosensitized | linoleate only, autoxidation and lipoxygenase |
| denominator | % of **total** peak area incl. an explicit "Unidentified peaks" row | % of a **selected six-peak** sum |
| internal standard | **none** | methyl hexanoate |
| replication | **none stated** | duplicate GC, RSD ±3.9-4.8 % |
| additive arm | none | alpha-tocopherol and 1,4-cyclohexadiene |
| rate, Ea, absolute yield | **none** | **none** |

**Do the product distributions agree? — a real cross-check is possible, and this is its result.**
Five of Frankel 1989's six products appear in this paper's linoleate table: pentane, hexanal, methyl
octanoate, 2,4-decadienal and methyl 9-oxononanoate. **The sixth, methyl 13-oxo-9,11-tridecadienoate,
is absent from the 1981 slate on purpose** — the Discussion says "All of the products expected from
the autoxidation hydroperoxides by Scheme I were detected **except the 12- and 13-carbon unsaturated
aldehyde esters, for which no authentic references were available**." So 1989 is 1981 plus the C13
ester reference compound. Renormalising both papers' autoxidation columns over the **same five
products** (mine, since neither paper prints this):

| product | 1981, 210 C, packed, 25 C start (share of the five) | 1989, 180 C, capillary, -65 C trap (share of the five) | 1981 / 1989 |
|---|---:|---:|---:|
| pentane | 13.6 | 20.0 | **0.68** |
| hexanal | 20.6 | 13.8 | **1.50** |
| methyl octanoate | 20.6 | 21.3 | **0.97** |
| 2,4-decadienal | 19.2 | 28.8 | **0.67** |
| methyl 9-oxononanoate | 26.1 | 16.3 | **1.60** |

**They agree to within a factor of 1.6 in either direction, and the disagreement has a shape.**
Methyl octanoate is identical (0.97). The two products Frankel 1989 identifies as reachable by
**both** the Hock and the homolytic-B route — hexanal and methyl 9-oxononanoate — are **higher** in
1981, and the two reachable only by homolytic pathway A — pentane and 2,4-decadienal — are
**lower**. Frankel 1989's own summary statistic makes this legible: **(hexanal + Me
9-oxononanoate)/(pentane + Me octanoate) is 1.37 here (mine, autoxidation) against 0.73 there** — a
factor of **1.87**. That ratio is **denominator-free** (both terms come from the same column), which
is why it is the only honest quantitative comparison between the two papers.

**Two readings, and the paper cannot settle between them.** (i) A real temperature effect: 210 C
favours the heterolytic/pathway-B products over pathway A, the same direction Frankel 1989's
hydrogen-donor arm pushed at 180 C. (ii) A method artefact: pentane boils at 36 C and 1981 started
the column at **25 C with no cryo-trap**, so the light end is the part most likely under-recovered —
and pentane is indeed the largest single discrepancy going down. **The 1981 paper itself warns of
exactly this**: "the volatile compounds analyzed are only those best amenable to our GC-MS detection
technique, and **their relative concentration is probably also affected by their thermal
stability**." Two studies, 30 C apart, with confounded methods, is not an activation energy and must
not be turned into one.

**And it changes what `LOOH_OL` means.** Table I shows that a photosensitized 9-/10-hydroperoxide
pool **isomerises inside the injector port** into 18 % 8- / 26 % 9- / 31 % 10- / 25 % 11-OOH — very
nearly the autoxidation composition (27/23/23/27). `species_lipid.LOOH_OL` lumps the four oleate
isomers; **this paper is the direct experimental justification for that lump**, at least under
injector-port conditions. The linoleate pools do **not** behave that way: a 50:50 9-/13- pool comes
out 47/2/4/47, i.e. essentially unchanged. So lumping is licensed for oleate and is **not** licensed
for linoleate — which is exactly the asymmetry the lane already has.

## 2. Methods as they matter to a model

- **The pot.** Again there is no pot: **reaction chromatography**. No solvent, no reaction time, no
  reactor volume, no conversion, no concentration. "Direct injection of different purified
  hydroperoxides onto the GC system affords **immediate anaerobic pyrolysis** necessary for
  efficient separation of volatile products."
- **Substrates.** "The same **methyl oleate (98 %), linoleate (99 %) and linolenate (99.9 %)** were
  used as described previously (12)."
- **Two oxidations.** Footnote *b* to Table I: "**Auto: free radical autoxidation in air (12); PV =
  peroxide value. Photo: photosensitized oxidation in O2, methanol solution in presence of methylene
  blue (15).**" The oxidation conditions, printed per row in Table I: oleate auto **40 C, PV 1051**;
  oleate photo **0 C, PV 1727**; linoleate auto **27 C, PV 2970**; linoleate photo **0 C, PV 1124**;
  linolenate auto **27 C, PV 790**; linolenate photo **0 C, PV 1566**.
- **Exactly which hydroperoxide isomers were isolated — the answer is "pools, characterised".**
  "**Two methods were used to purify hydroperoxides. Silicic acid column chromatography with
  methanolic benzene (22) was used to isolate the hydroperoxides of autoxidized oleate.
  Reverse-phase high pressure liquid chromatography was used with a 1:1 (w/w) H2O/acetonitrile
  solvent system (23) to isolate all other hydroperoxides.**" Purity checked by TLC (silica gel,
  ether/petroleum ether 60:40). **No single positional isomer was isolated in this paper.** Six
  mixed pools were prepared and each pool's isomer distribution was measured (Table I) by GC-MS of
  the hydroxystearate TMS ethers after hydrogenation and silylation. Contrast Frankel 1989, which
  did isolate a pure cis,trans-13 hydroperoxide by lipoxygenase.
- **The thermolysis.** "Hydroperoxides were decomposed in the injector port of a gas chromatograph
  (on column injection Model 7400, Packard Instrument Co.) **at 210 C**. The He carrier flow was
  decreased to ca. 5 mL/min just prior to and 1 min after injection of a **neat sample (4 to 10
  µL)**. The flow was then increased (ca. 36 mL/min) and **temperature programming started from 25
  to 275 C at 2 C/min**. The GC column (**glass 14 ft x 4 mm i.d.**) was packed with **10 % OV-101 on
  Chromosorb G.**"
- **The separate isomerisation experiment, at a different temperature.** A glass liner packed with
  the same OV-101 and plugged with siliconized glass wool was inserted into the injector of a
  Hewlett-Packard 5711A; after injection the carrier flow was cut for 1 min, the port disassembled
  and the liner **cooled immediately in Dry Ice**, then rewarmed and purged with absolute ethanol;
  the recovered material was hydrogenated (PtO2, atmospheric H2), silylated and analysed by GC-MS.
  **Table I's column head says "Pyrolysis (200 C)" while the decomposition experiments were run at
  210 C** — see flag 8.
- **Identification.** "…to identify volatile compounds by **matching mass spectra with those of
  reference compounds and confirming by GC-retention data**." Tentative identifications are
  footnoted per table.
- **Quantification.** "rel %", with no further statement of method, no internal standard, no
  response factors, no replicates and no error bar. Every product column carries an explicit
  **"Unidentified peaks"** row and sums to ~100 (mine: oleate auto **100.0**, oleate photo **99.5**,
  linoleate auto **99.9**, linoleate photo **99.9**, linolenate auto **100.9**, linolenate photo
  **99.6**), which fixes the denominator as total volatile peak area.
- **The mechanism (Scheme I), as stated.** Carbon-carbon scission on either side (A and B) of the
  alkoxy radical intermediate, plus the same 1-enol/hydroxyl-radical tautomerisation Selke 1978
  invoked for the saturated aldehydes. **Scheme II is new in this paper**: hydroperoxy cyclic
  peroxides, formed by cyclisation of the internal 10- and 12-hydroperoxides of linoleate and the
  12-/13-/15-hydroperoxides of linolenate, are proposed as a **second precursor class** that also
  yields pentane, hexanal, methyl octanoate and methyl 9-oxononanoate. **This matters for the lane's
  branch model**: if it is right, hexanal and methyl 9-oxononanoate have parents outside the
  9-/13-hydroperoxide pair that `species_lipid` models. The paper offers it as a proposal supported
  by an absence ("the absence of these isomers (Table I) supports cyclization"), not a measurement.
- **What the paper says about rates.** Nothing of its own. It repeats an *unpublished* claim
  attributed to Chan and Levett (34) that "**the 9- and 13-linolenate hydroperoxides (prepared by
  lipoxygenase action) decompose at the same rate**", and its Discussion says the two competing
  processes (decomposition and rearrangement) "**may be competitive and controlled kinetically**".
  **No rate constant, no half-life, no activation energy, no time axis exists in this paper.**

## 3. Tables re-typed

Marks: `[M]` measured in this paper, `[C]` cited/reprinted from another paper, `[F]` fitted/assigned
by the authors.

### TABLE I (p. 280). "GC-MS Analysis (12) of Isomeric Hydroxystearates from Hydroperoxides before and after Partial GC Pyrolysis*a*"

Footnotes: *a* "Conditions given in Experimental Methods." *b* "Auto: free radical autoxidation in
air (12); PV = peroxide value. Photo: photosensitized oxidation in O2, methanol solution in presence
of methylene blue (15)." Column head over the numbers: "**Pyrolysis (200 C)**" and "**Relative
percent**".

| Hydroperoxide esters | Oxidation (conditions)*b* | Pyrolysis (200 C) | 8-OH | 9-OH | 10-OH | 11-OH |
|---|---|---|---:|---:|---:|---:|
| Oleate | Auto (40 C, PV 1051) | Before | 27 `[C]` | 23 `[C]` | 23 `[C]` | 27 `[C]` |
| Oleate | Photo (0 C, PV 1727) | Before | | 50 `[M]` | 50 `[M]` | |
| Oleate | Photo (0 C, PV 1727) | **After** | **18** `[M]` | **26** `[M]` | **31** `[M]` | **25** `[M]` |

| Hydroperoxide esters | Oxidation (conditions)*b* | Pyrolysis (200 C) | 9-OH | 10-OH | 12-OH | 13-OH |
|---|---|---|---:|---:|---:|---:|
| Linoleate | Auto (27 C, PV 2970) | Before | 50 `[M]` | | | 50 `[M]` |
| Linoleate | Auto (27 C, PV 2970) | After | 47 `[M]` | 2 `[M]` | 4 `[M]` | 47 `[M]` |
| Linoleate | Photo (0 C, PV 1124) | Before | 32 `[M]` | 17 `[M]` | 17 `[M]` | 34 `[M]` |
| Linoleate | Photo (0 C, PV 1124) | After | 28 `[M]` | 19 `[M]` | 21 `[M]` | 32 `[M]` |

| Hydroperoxide esters | Oxidation (conditions)*b* | Pyrolysis (200 C) | 9-OH | 10-OH | 12-OH | 13-OH | 15-OH | 16-OH |
|---|---|---|---:|---:|---:|---:|---:|---:|
| Linolenate | Auto (27 C, PV 790) | Before | 32 `[M]` | | 11 `[M]` | 11 `[M]` | | 46 `[M]` |
| Linolenate | Auto (27 C, PV 790) | After | 31 `[M]` | | 10 `[M]` | 14 `[M]` | | 45 `[M]` |
| Linolenate | Photo (0 C, PV 1566) | Before | 23 `[M]` | 13 `[M]` | 12 `[M]` | 14 `[M]` | 13 `[M]` | 25 `[M]` |
| Linolenate | Photo (0 C, PV 1566) | After | 22 `[M]` | 14 `[M]` | 8 `[M]` | 13 `[M]` | 15 `[M]` | 28 `[M]` |

**Arithmetic check (mine): every one of the eleven populated rows sums to exactly 100.**

### TABLE II (p. 281). "GC-MS Analysis of Volatiles from Thermally Decomposed Methyl Oleate Hydroperoxides"

Footnotes: *a* "**Data from ref. 21.**" *b* "Based on Scheme I (5,21)." *c* "Tentative
identification."

| Compound | Elution temp (C) | Autoxidation*a* (rel %) | Photosensitized oxidation (rel %) | Origin*b* |
|---|---:|---:|---:|---|
| Heptane | 106 | 4.4 `[C]` | 4.6 `[M]` | 11-OOH `[F]` |
| Octane | 121 | 2.7 `[C]` | 10 `[M]` | 10-OOH `[F]` |
| Heptanal | 151 | 0.5 `[C]` | 0.5 `[M]` | ? |
| 1-Heptanol | 161 | 0.4 `[C]` | 0.4 `[M]` | 11-OOH `[F]` |
| Octanal | 169 | 11 `[C]` | 3.8 `[M]` | 11-OOH `[F]` |
| Me heptanoate | 170 | 1.5 `[C]` | 4.9 `[M]` | 8-OOH `[F]` |
| 1-Octanol | 181 | 0.4 `[C]` | 1.0 `[M]` | 10-OOH `[F]` |
| **Nonanal** | **186** | **15** `[C]` | **10** `[M]` | **9-/10-OOH** `[F]` |
| Me octanoate | 189 | 5.0 `[C]` | 9.7 `[M]` | 9-OOH `[F]` |
| 2-Nonenal | 197 | 0.5 `[C]` | 0.7 `[M]` | ? |
| Decanal | 201 | 3.9 `[C]` | 2.0 `[M]` | 8-OOH `[F]` |
| Me nonanoate | 203 | 1.5 `[C]` | 0.8 `[M]` | ? |
| 2-Decenal | 211 | 5.4 `[C]` | 12 `[M]` | 9-OOH `[F]` |
| 2-Undecenal | 225 | 1.7 `[C]` | 7.1 `[M]` | 8-OOH `[F]` |
| Me 8-Oxooctanoate | 230 | 3.5 `[C]` | 3.0 `[M]` | 8-OOH `[F]` |
| Me 9-Oxononanoate | 245 | 15 `[C]` | 11 `[M]` | 9-/10-OOH `[F]` |
| Me 10-Oxodecanoate | 256 | 12 `[C]` | 1.7 `[M]` | 11-OOH `[F]` |
| Me 10-Oxo-8-decenoate*c* | 265 | 3.4 `[C]` | 5.0 `[M]` | 10-OOH `[F]` |
| Me 11-Oxo-9-undecenoate*c* | 275 | 5.8 `[C]` | 4.6 `[M]` | 11-OOH `[F]` |
| Unidentified peaks | | 6.4 `[C]` | 6.7 `[M]` | |

Column sums (mine): **100.0** and **99.5**. **The autoxidation column is Selke 1978's Table I,
cell for cell** — including Me nonanoate 1.5, Me 9-oxononanoate 15 and Me 10-oxodecanoate 12, which
independently confirms the raster reading in `selke1978_extraction.md` against the OCR corruption
there.

### TABLE III (p. 282). "GC-MS Analysis of Volatiles from Thermally Decomposed Methyl Linoleate Hydroperoxides"

Footnotes: *a* "Based on Scheme I (5,21); \*isomerized." *b* "**Not separated by GC, estimated by
MS.**" *c* "Different peaks due to geometric isomers." *d* "Tentative identification."

| Compound | Elution temp (C) | Autoxidation (rel %) | Photosensitized oxidation (rel %) | Origin*a* |
|---|---:|---:|---:|---|
| Acetaldehyde | 70 | 0.3 `[M]` | 0.4 `[M]` | ? |
| Pentane | 88 | 9.9 `[M]` | 4.3 `[M]` | 13-OOH `[F]` |
| Pentanal | 117 | 0.8 `[M]` | 0.3 `[M]` | 13-OOH `[F]` |
| 1-Pentanol | 129 | 1.3 `[M]` | 0.3 `[M]` | 13-OOH `[F]` |
| **Hexanal** | **136** | **15** `[M]` | **17** `[M]` | **12-/13-OOH** `[F]` |
| 2-Heptenal*b* | 165 | Tr `[M]` | 9.9 `[M]` | 12-OOH `[F]` |
| 1-Octen-3-ol*b* | 165 | Tr `[M]` | 1.9 `[M]` | 10-OOH `[F]` |
| **2-Pentylfuran*b*** | **165** | **2.4** `[M]` | **0.6** `[M]` | **?** |
| Me heptanoate | 170 | 1.0 `[M]` | 0.3 `[M]` | ? |
| 2-Octenal | 182 | 2.7 `[M]` | 1.5 `[M]` | ? |
| Me octanoate | 189 | 15 `[M]` | 7.6 `[M]` | 9-OOH `[F]` |
| 2-Nonenal | 195, 197*c* | 1.4 `[M]` | 1.6 `[M]` | 9-/10-OOH\* `[F]` |
| 2,4-Nonadienal | 208 | 0.3 `[M]` | 0.3 `[M]` | ? |
| 2,4-Decadienal | 219, 223*c* | 14 `[M]` | 4.3 `[M]` | 9-OOH `[F]` |
| Me 8-Oxooctanoate | 230 | 1.3 `[M]` | 0.9 `[M]` | ? |
| Me 9-Oxononanoate | 245 | 19 `[M]` | 22 `[M]` | 9-/10-OOH `[F]` |
| Me 10-Oxodecanoate | 256 | 0.7 `[M]` | 0.7 `[M]` | ? |
| Me 10-Oxo-8-decenoate*d* | 265 | 4.9 `[M]` | 14 `[M]` | 10-OOH `[F]` |
| Unidentified peaks | | 9.9 `[M]` | 12 `[M]` | |

Column sums (mine, counting "Tr" as 0): **99.9** and **99.9**. **"Tr" = trace, unquantified.**
**Note that 2-heptenal, 1-octen-3-ol and 2-pentylfuran all elute at 165 C, were NOT separated by
GC, and were estimated by MS** — see flag 5.

### TABLE IV (p. 283). "GC-MS Analysis of Volatiles from Thermally Decomposed Methyl Linolenate Hydroperoxides"

Footnotes: *a* "Based on Scheme I (5,21)." *b* "Different peaks due to geometric isomers." *c*
"Identified as 2,6-nonadienal." *d* "Tentative identification."

| Compound | Elution temp (C) | Autoxidation (rel %) | Photosensitized oxidation (rel %) | Origin*a* |
|---|---:|---:|---:|---|
| Ethane/ethene | 65 | 10 `[M]` | 3.2 `[M]` | 16-OOH `[F]` |
| Acetaldehyde | 70 | 0.8 `[M]` | 0.6 `[M]` | ? |
| **Propanal/acrolein** | **80** | **7.7** `[M]` | **9.0** `[M]` | **15-/16-OOH** `[F]` |
| Butanal | 97 | 0.1 `[M]` | 0.8 `[M]` | ? |
| 2-Butenal | 109 | 0.5 `[M]` | 11 `[M]` | 15-OOH `[F]` |
| 2-Pentenal | 131 | 1.6 `[M]` | 1.2 `[M]` | 13-OOH `[F]` |
| 2-/3-Hexenal | 137 | 1.4 `[M]` | 3.4 `[M]` | 12-/13-OOH `[F]` |
| **2-Butylfuran** | **158** | **0.5** `[M]` | **0.3** `[M]` | **?** |
| Me heptanoate | 170 | 1.8 `[M]` | 1.0 `[M]` | ? |
| 2,4-Heptadienal | 174, 178*b* | 9.3 `[M]` | 8.8 `[M]` | 12-OOH `[F]` |
| Me octanoate | 189 | 22 `[M]` | 15 `[M]` | 9-OOH `[F]` |
| 4,5-Epoxyhepta-2-enal | 194 | 0.2 `[M]` | 0.2 `[M]` | ? |
| 3,6-Nonadienal*c* | 196, 198*b* | 0.5 `[M]` | 1.1 `[M]` | 9-/10-OOH `[F]` |
| Me Nonanoate | 203 | 0.7 `[M]` | 0.3 `[M]` | ? |
| Decatrienal | 219, 226*b* | 14 `[M]` | 4.8 `[M]` | 9-OOH `[F]` |
| Me 8-Oxooctanoate | 230 | 0.6 `[M]` | 0.4 `[M]` | ? |
| Me 9-Oxononanoate | 245 | 13 `[M]` | 12 `[M]` | 9-/10-OOH `[F]` |
| Me 10-Oxodecanoate | 256 | 1.0 `[M]` | 1.5 `[M]` | ? |
| Me 10-Oxo-8-decenoate*d* | 267 | 4.2 `[M]` | 13 `[M]` | 10-OOH `[F]` |
| Unidentified peaks | | 11 `[M]` | 12 `[M]` | |

Column sums (mine): **100.9** and **99.6**.

### Statements printed in the running text, transcribed

- **On oleate:** "Although the two starting hydroperoxide mixtures are isomerically different, **both
  samples formed the same volatile products**." "**Photooxidized oleate hydroperoxides produced not
  only all the volatiles expected from the 9- and 10-isomers but also those expected from the 8- and
  11-isomers.** The photosensitized oxidation-derived hydroperoxides produced much more octane,
  2-decenal, 2-undecenal, 1-octanol, methyl heptanoate and octanoate, and much less octanal and
  methyl 10-oxodecanoate than the autoxidation-derived hydroperoxides."
- **The isomerisation result:** "the mixture of 9- and 10-hydroperoxides from photosensitized
  oxidation of oleate **isomerized into a mixture of 8-, 9-, 10- and 11-hydroperoxides** (Table I)."
- **On linoleate:** "The autoxidized linoleate hydroperoxides produced **much more pentane,
  2-pentylfuran, 2,4-decadienal, and methyl octanoate** and much less methyl 10-oxo-8-decenoate and
  2-heptenal than the photooxidized linoleate hydroperoxides."
- "Under our thermal decomposition conditions, **very little interconversion of linoleate
  hydroperoxide mixtures occurred.** The 50:50 mixture of 9- and 13-hydroperoxides in the
  autoxidation sample, after heating in the injector port, produced a mixture containing also **2 %
  10- and 4 % 12-hydroperoxides**."
- "**2-Pentylfuran is a unique product of autoxidation hydroperoxide, but its origin is not well
  established.** 2-Heptenal is a unique product of photosensitized oxidation hydroperoxides and
  would originate from the 12-hydroperoxide by cleavage A according to Scheme I."
- "In contrast to Chan et al. (9), **we identified 2-enals and 2-pentylfuran among the volatiles from
  the autoxidation hydroperoxides.** Although their origin is not clear, 2-octenal is speculated to
  come from the nonconjugated 11-hydroperoxide (16-18) and **2-pentylfuran from a 10-hydroperoxide
  intermediate (27)**."
- "All of the products expected from the autoxidation hydroperoxides by Scheme I were detected
  **except the 12- and 13-carbon unsaturated aldehyde esters, for which no authentic references were
  available**."
- **On linolenate:** "The most distinguishing volatiles include **ethane for autoxidation
  hydroperoxides and 2-butenal for the photooxidation hydroperoxides**." Products expected but not
  identified include "the methyl C-12, C-13, C-14, C-15 and C-16 unsaturated aldehyde esters … as
  well as the unsaturated hydrocarbons (**2-pentene and 2,5-octadiene**)".
- **On malondialdehyde:** "**Malonaldehyde is another expected product** of cleavage on either side
  of the cyclic peroxide … **However, this dialdehyde was not identified under our GC-MS
  conditions.**"
- **The paper's own limitation statement:** "the volatile compounds analyzed are **only those best
  amenable to our GC-MS detection technique**, and their relative concentration is **probably also
  affected by their thermal stability**."
- **Its closing sentence:** "Therefore, contrary to the view of Chan et al. (9), these **secondary
  products may make an important contribution even under anaerobic pyrolysis of pure
  hydroperoxides. Further study is necessary to establish a mechanism of decomposition involving the
  secondary products of hydroperoxides.**"

### Derived numbers (mine, arithmetic on the printed tables — NOT the paper's)

- **2-pentylfuran / hexanal, same column:** autoxidised linoleate **2.4/15 = 0.160**; photosensitized
  linoleate **0.6/17 = 0.0353**. Ratio between the two oxidations: **4.5x** more alkylfuran per
  hexanal from autoxidation.
- **2-pentylfuran as a share of the five Frankel-1989-slate products present** (pentane 9.9 +
  hexanal 15 + Me octanoate 15 + 2,4-decadienal 14 + Me 9-oxononanoate 19 = 72.9): **2.4/72.9 =
  3.29 %** (autoxidation).
- **The five-product renormalisation and the 1981/1989 comparison**: as tabulated in section 1.
- **Frankel 1989's Figure-4 ratio, (hexanal + Me 9-oxononanoate)/(pentane + Me octanoate)**:
  1981 autoxidised linoleate **(15+19)/(9.9+15) = 1.37**; 1981 photosensitized linoleate
  **(17+22)/(4.3+7.6) = 3.28**; against 1989's **0.73 / 0.93 / 2.69** at zero additive. The 1981
  photosensitized value sits between 1989's tt-isomer value and beyond it.
- **pentane / hexanal**: 1981 auto **9.9/15 = 0.66**; 1989 mixed auto **16/11 = 1.45**. The single
  largest inter-paper discrepancy, factor **2.2**, and the one most plausibly a light-end recovery
  artefact of the 25 C start.
- **nonanal / methyl 9-oxononanoate (the strict 1:1 scission pairing)**: autoxidised oleate
  **15/15 = 1.00**; photosensitized oleate **10/11 = 0.91**. Both essentially exact — the tightest
  mechanistic closure anywhere in the Frankel-lab corpus.
- **Identified fraction of each chromatogram**: oleate auto **93.6 %**, oleate photo **93.3 %**,
  linoleate auto **90.1 %**, linoleate photo **88 %**, linolenate auto **~89.9 %**, linolenate photo
  **88 %** (100 minus each printed "Unidentified peaks" row).
- **Sum of the isomer distributions**: all eleven populated Table I rows close to exactly 100.

## 4. Numbers the repository can use

**All rows share: neat hydroperoxide (4-10 µL) thermolysed in a GC injector port at 210 C, packed
OV-101 column programmed from 25 C at 2 C/min, GC-MS identification against reference compounds,
quantified as relative peak area of the total chromatogram with NO internal standard, NO response
factors, NO replicates and NO stated error.**

### The two refused compounds

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **linoleate LOOH → 2-pentylfuran, share of the whole chromatogram** | **2.4** | % of total volatile peak area | **autoxidised** methyl linoleate hydroperoxides (27 C, PV 2970; 50:50 9-/13-OOH before pyrolysis) | Table III, p. 282 | **branch_fraction** — the number the panel says does not exist. **Co-eluting, MS-estimated; see flag 5** |
| **linoleate LOOH → 2-pentylfuran, photosensitized** | **0.6** | % of total volatile peak area | **photosensitized** methyl linoleate hydroperoxides (0 C, PV 1124; 32/17/17/34) | Table III, p. 282 | **branch_fraction** — a second, independent pool |
| **2-pentylfuran / hexanal, same column** | **0.160** (autoxidation), **0.0353** (photosensitized) — both mine | mol-free peak-area ratio | as above | derived from Table III | **within_study_ratio — THE RECOMMENDED FORM.** Denominator-free, and hexanal is already a modelled species, so this hangs the alkylfuran on an existing node without importing 1981's denominator |
| **oleate LOOH → nonanal, autoxidation** | **15** | % of total volatile peak area | autoxidised methyl oleate LOOH, 27/23/23/27 8-/9-/10-/11-OOH (40 C, PV 1051) | Table II, p. 281 | **branch_fraction — but `[C]`, footnoted "Data from ref. 21".** This is Selke 1978's single measurement republished, NOT a replicate |
| **oleate LOOH → nonanal, photosensitized** | **10** | % of total volatile peak area | photosensitized methyl oleate LOOH, 50:50 9-/10-OOH before pyrolysis, isomerising in the port to 18/26/31/25 | Table II, p. 281 | **branch_fraction — genuinely new here.** The only independent second oleate nonanal determination in the corpus |
| nonanal / methyl 9-oxononanoate | 1.00 (auto), 0.91 (photo) — mine | — | as above | derived from Table II | within_study_ratio |
| nonanal, assigned origin | **9-/10-OOH** | — | — | Table II Origin column | **`[F]` assigned by Scheme I**, not measured |
| 2-pentylfuran, assigned origin | **"?"** — unassigned; Discussion speculates a **10-hydroperoxide intermediate**, citing ref (27) = Chang et al., *Chem. Ind.* 1926 (1966) | — | — | Table III Origin column; Discussion p. 284 | **NOT a measurement, and contradicted by the paper's own Table I — see flag 6** |

### The rest of the linoleate slate (the lane's own substrate)

| product | autoxidation (rel %) | photosensitized (rel %) | assigned origin | evidence class |
|---|---:|---:|---|---|
| pentane | 9.9 | 4.3 | 13-OOH | branch_fraction |
| hexanal | 15 | 17 | 12-/13-OOH | branch_fraction |
| methyl octanoate | 15 | 7.6 | 9-OOH | branch_fraction |
| 2,4-decadienal | 14 | 4.3 | 9-OOH | branch_fraction |
| methyl 9-oxononanoate | 19 | 22 | 9-/10-OOH | branch_fraction |
| **methyl 13-oxo-9,11-tridecadienoate** | **ABSENT** | **ABSENT** | — | not identified: **"no authentic references were available"** for the C12/C13 unsaturated aldehyde esters |
| pentanal | 0.8 | 0.3 | 13-OOH | branch_fraction |
| 1-pentanol | 1.3 | 0.3 | 13-OOH | branch_fraction |
| 2-heptenal | Tr | 9.9 | 12-OOH | branch_fraction; **Tr is unquantified**; co-eluted at 165 C |
| 1-octen-3-ol | Tr | 1.9 | 10-OOH | branch_fraction; co-eluted at 165 C |
| 2-octenal | 2.7 | 1.5 | ? | branch_fraction |
| 2-nonenal | 1.4 | 1.6 | 9-/10-OOH, isomerized | branch_fraction — **the Hock partner Frankel 1989 names and never measures, measured here from the right substrate** |
| 2,4-nonadienal | 0.3 | 0.3 | ? | branch_fraction |
| acetaldehyde | 0.3 | 0.4 | ? | branch_fraction |
| methyl heptanoate | 1.0 | 0.3 | ? | branch_fraction |
| methyl 8-oxooctanoate | 1.3 | 0.9 | ? | branch_fraction |
| methyl 10-oxodecanoate | 0.7 | 0.7 | ? | branch_fraction |
| methyl 10-oxo-8-decenoate | 4.9 | 14 | 10-OOH | branch_fraction; tentative identification |
| **unidentified peaks** | **9.9** | **12** | — | **measured_bound — the direct analogue of `LIPID_FRAG_C`** |

### The oleate slate (the `LOOH_OL` pool)

Full column as re-typed in section 3, Table II. The autoxidation column is `[C]` from Selke 1978;
the photosensitized column is `[M]`. Headline entries: octanal 11 / 3.8, decanal 3.9 / 2.0,
2-decenal 5.4 / 12, 2-undecenal 1.7 / 7.1, methyl octanoate 5.0 / 9.7, methyl 9-oxononanoate 15 / 11,
methyl 10-oxodecanoate 12 / 1.7, unidentified 6.4 / 6.7.

### The linolenate slate — new species the corpus has never had a number for

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **linolenate LOOH → propanal/acrolein** | **7.7** (auto), **9.0** (photo) | % of total volatile peak area | 210 C, autoxidised (27 C, PV 790) / photosensitized (0 C, PV 1566) linolenate LOOH | Table IV, p. 283 | **branch_fraction — but the peak is "Propanal/acrolein", TWO compounds in one row.** `parameters_lipid.PROHIBITED_DERIVATIONS` refuses importing a propanal; this measures the pair, not propanal alone. The refusal should stand in its current form |
| linolenate LOOH → 2-butylfuran | 0.5 (auto), 0.3 (photo) | % of total volatile peak area | as above | Table IV | branch_fraction — the **C4 alkylfuran homologue**, unassigned origin ("?"), 4.8x smaller than 2-pentylfuran from linoleate |
| linolenate LOOH → ethane/ethene | 10 (auto), 3.2 (photo) | % | as above | Table IV | branch_fraction — again **two compounds in one row** |
| linolenate LOOH → 2,4-heptadienal | 9.3, 8.8 | % | as above | Table IV | branch_fraction |
| linolenate LOOH → methyl octanoate | 22, 15 | % | as above | Table IV | branch_fraction |
| linolenate LOOH → decatrienal | 14, 4.8 | % | as above | Table IV | branch_fraction |
| linolenate LOOH → methyl 9-oxononanoate | 13, 12 | % | as above | Table IV | branch_fraction |
| linolenate LOOH → 2-butenal | 0.5, **11** | % | as above | Table IV | branch_fraction — the paper's "distinguishing volatile" for photosensitized linolenate |
| linolenate LOOH → 2-/3-hexenal | 1.4, 3.4 | % | as above | Table IV | branch_fraction — two compounds in one row |
| linolenate LOOH → 4,5-epoxyhepta-2-enal | 0.2, 0.2 | % | as above | Table IV | branch_fraction |
| **malondialdehyde** | **NOT IDENTIFIED** — "expected … However, this dialdehyde was not identified under our GC-MS conditions" | — | — | Discussion p. 284 | measured_bound (a stated non-detection, with no detection limit) |

### Hydroperoxide pool compositions — direct parameters for `LOOH_OL` and its linoleate sibling

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| oleate autoxidation isomer distribution | 27 / 23 / 23 / 27 | % 8- / 9- / 10- / 11-OOH | air, 40 C, PV 1051 | Table I | measured_ratio `[C]` (same as Selke 1978) |
| oleate photosensitized isomer distribution, **before** pyrolysis | 50 / 50 | % 9- / 10-OOH | O2/methanol/methylene blue, 0 C, PV 1727 | Table I | measured_ratio |
| **oleate photosensitized isomer distribution, AFTER pyrolysis** | **18 / 26 / 31 / 25** | % 8- / 9- / 10- / 11-OOH | 200 C liner pyrolysis | Table I | **measured_ratio — the experimental licence for lumping `LOOH_OL`** |
| linoleate autoxidation, before / after | 50/-/-/50 → 47 / 2 / 4 / 47 | % 9- / 10- / 12- / 13-OOH | air, 27 C, PV 2970; 200 C liner | Table I | measured_ratio — **essentially no isomerisation; lumping is NOT licensed for linoleate** |
| linoleate photosensitized, before / after | 32/17/17/34 → 28/19/21/32 | % 9- / 10- / 12- / 13-OOH | 0 C, PV 1124 | Table I | measured_ratio |
| linolenate autoxidation, before / after | 32/-/11/11/-/46 → 31/-/10/14/-/45 | % 9-/10-/12-/13-/15-/16-OOH | air, 27 C, PV 790 | Table I | measured_ratio |
| linolenate photosensitized, before / after | 23/13/12/14/13/25 → 22/14/8/13/15/28 | % as above | 0 C, PV 1566 | Table I | measured_ratio |
| peroxide values of the six feeds | 1051, 1727, 2970, 1124, 790, 1566 | (meq/kg — **unit not printed**) | as above | Table I | level_only |

### Absent quantities

| quantity | status |
|---|---|
| **rate constant at any temperature** | **NOT PRESENT** |
| **activation energy** | **NOT PRESENT** — one thermolysis temperature (210 C), one pyrolysis temperature (200 C), no time axis |
| **absolute yield (mol/mol, mass, mmol/L)** | **NOT PRESENT** |
| **total volatile output, comparable between substrates** | **NOT PRESENT** — no internal standard; every column is normalised to 100 independently, so **oleate output cannot be compared to linoleate output** |
| **any error bar, replicate or n** | **NOT PRESENT** |
| the "same rate" claim for 9- vs 13-linolenate hydroperoxide | `[C]` and explicitly **unpublished**, attributed to Chan & Levett (34). Not usable |

### **WHICH KIND OF QUANTITY IS "2.4 %" AND "15 %"? — the question the brief demands be answered**

**Both are relative GC peak areas, expressed as a percent of the total peak area of that
injection's whole chromatogram.** Specifically:

1. **NOT a fraction of the hydroperoxide consumed.** No conversion is measured. Nothing in this
   paper says how much hydroperoxide decomposed or how much went to non-volatile material. (Selke
   1978's introduction cites a 90 % polymeric / 10 % volatile split for heated fats; **that citation
   does not reappear in this paper and must not be imported through it**.)
2. **It IS a share of total volatiles in the peak-area sense**, and it is a well-posed one: each
   column carries an explicit "Unidentified peaks" row and closes to 100 within rounding (mine:
   100.0, 99.5, 99.9, 99.9, 100.9, 99.6). This is a *better-posed* denominator than Frankel 1989's
   six-peak sum, which excludes everything outside the six.
3. **It is an AREA share, not molar or mass.** No response factors, no internal standard, no
   calibration. Over a slate running from ethane to a C14 oxo-ester on a packed OV-101 column, area
   is not moles. **Do not convert 2.4 % or 15 % to mmol/L.**
4. **The only quantity in this paper that survives transfer to another study is a RATIO taken within
   one column** — e.g. 2-pentylfuran/hexanal = 0.160, or (hexanal + Me 9-oxononanoate)/(pentane + Me
   octanoate) = 1.37. Those are denominator-free and response-factor-free only to the extent that
   the two compounds' FID/MS responses are similar, which for hexanal against 2-pentylfuran (C6H12O
   vs C9H14O) they are not exactly. Treat the ratio as good to a factor of ~1.5, not better.

**Recommended form for the two refused branches.**
`linoleate_looh_to_2_pentylfuran_per_hexanal = 0.16` (autoxidation pool; `= 0.035` photosensitized),
method `injector_port_thermolysis_210C`, `co_eluting = True`, `ms_estimated = True`,
`parent_isomer = unassigned`, `n = 1`, `error = none stated`.
`oleate_looh_to_nonanal = 0.15` (autoxidation, `[C]` = Selke 1978) and `= 0.10` (photosensitized,
`[M]`), denominator = total volatile peak area including unidentified, same method flags.
**Whether that is enough to lift the refusals is a decision for the lane owner, not for this
dossier**; sections 5's flags are the case against, and they are not small.

## 5. Flags

1. **The paper contains no rate, no Ea and no absolute yield — the same hole as Frankel 1989.** One
   injector temperature. The lane's Q10 assumption is untouched. Two papers 30 C apart with
   confounded columns, trapping and substrates do **not** make a two-point Arrhenius fit, and any
   attempt to build one from 180 C and 210 C should be refused.
2. **The oleate autoxidation column is not an independent measurement.** Table II footnote *a* reads
   "Data from ref. 21" — Selke 1978. **Nonanal = 15 % appears in two papers and was measured once.**
   Counting it as two determinations would be double-counting. The photosensitized column (nonanal
   = 10 %) is the only genuinely second oleate determination, and it is from a *different* isomer
   pool that isomerised in the port.
3. **n = 1 throughout, no error bars, no replicates.** Single neat injections, 1981 packed-column
   GC-MS. No RSD is quoted anywhere — this paper is *weaker* on that axis than Frankel 1989, which
   at least reports duplicate-analysis RSDs of ±3.9-4.8 %.
4. **No internal standard, so nothing can be compared between tables.** Each column is normalised to
   its own 100. **The oleate, linoleate and linolenate slates are three separate normalisations of
   three separate injections at unstated loadings.** Do not compute "how much more hexanal than
   nonanal a mixed fat gives" from Tables II and III.
5. **2-Pentylfuran was NOT resolved by gas chromatography.** Footnote *b* to Table III:
   "**Not separated by GC, estimated by MS.**" It shares an elution temperature of **165 C** with
   **2-heptenal** and **1-octen-3-ol**, and all three carry that footnote. So the 2.4 % is a
   mass-spectral deconvolution of a co-eluting triplet, on 1981 instrumentation, with no stated
   deconvolution method and no uncertainty. **In the autoxidation column the two co-eluting partners
   are both "Tr" (trace)**, which makes the 2.4 % the cleanest of the three — but in the
   photosensitized column 2-heptenal is 9.9 and 1-octen-3-ol 1.9, so the 0.6 % there is a small
   signal extracted from a 16x larger neighbour. **The photosensitized 2-pentylfuran value is the
   less trustworthy of the two and the 4.5x auto/photo contrast may be partly a deconvolution
   artefact.**
6. **The paper's speculated origin for 2-pentylfuran is contradicted by the paper's own Table I.**
   The Discussion attributes 2-pentylfuran to "a 10-hydroperoxide intermediate (27)". But Table I
   shows the **autoxidised** linoleate pool contains **0 % 10-OOH before pyrolysis and 2 % after**,
   while the **photosensitized** pool contains **17 % before and 19 % after** — and the autoxidised
   pool yields **4x MORE** 2-pentylfuran (2.4 vs 0.6), not less. **The 10-hydroperoxide hypothesis
   predicts the opposite of what the tables show** (this inference is mine, from Tables I and III;
   the paper does not make it and does not notice the tension). Consequence for the model: **there
   is no defensible way to attach the alkylfuran to a named hydroperoxide isomer.** It can only be
   carried as a share of the whole autoxidised-linoleate pool, or as a ratio to hexanal.
7. **Scheme II undermines the single-parent branch model the lane assumes.** The Discussion proposes
   that hydroperoxy cyclic peroxides — formed by cyclisation of the internal 10- and
   12-hydroperoxides — also produce **pentane, hexanal, methyl octanoate and methyl
   9-oxononanoate**. If true, four of the lane's six modelled products have a second precursor class
   that `species_lipid` does not represent, and their shares are not a clean function of the 9-/13-
   split. The evidence offered is an absence (the 8- and 14-hydroperoxides that isomerisation would
   have produced were not seen). **Recorded as a declared mechanistic gap, not as a number.**
8. **A printed temperature inconsistency.** Table I's column head says "**Pyrolysis (200 C)**"; the
   Experimental Methods say the hydroperoxides were decomposed in the injector port "**at 210 C**".
   The liner-isomerisation experiment used a *different instrument* (HP 5711A vs Packard 7400), so
   200 C may be that instrument's setting, but the paper never reconciles them. **The isomer
   redistribution that licenses lumping `LOOH_OL` was measured at 200 C; the product slates were
   measured at 210 C.**
9. **Three table rows are two compounds each.** "Propanal/acrolein" (7.7 / 9.0), "Ethane/ethene"
   (10 / 3.2), "2-/3-Hexenal" (1.4 / 3.4). **None of these can be assigned to a single species.**
   In particular, `parameters_lipid.PROHIBITED_DERIVATIONS`'s refusal to import a propanal number
   should be **retained**: this paper measures propanal+acrolein together, and acrolein is a
   separate hazard-relevant species the trunk would want kept apart.
10. **Six tentative or unconfirmed identifications.** Methyl 10-oxo-8-decenoate and methyl
    11-oxo-9-undecenoate (Table II, footnote *c*), methyl 10-oxo-8-decenoate (Tables III and IV,
    footnote *d*), 3,6-nonadienal ("**Identified as 2,6-nonadienal**", Table IV footnote *c* — the
    row label and the footnote disagree about the double-bond positions). And "2-Nonenal" in Table
    III carries the origin "9-/10-OOH\*" where \* means "isomerized", i.e. the assignment itself is
    contingent.
11. **"Tr" is not zero and is not a number.** 2-Heptenal and 1-octen-3-ol in the autoxidation
    linoleate column are "Tr". Treat as `> 0, unquantified`; my column sums count them as 0, which
    is why the column reads 99.9 rather than 100.
12. **The 25 C column start is the largest known method difference from Frankel 1989 and it acts on
    exactly the light end.** Pentane (b.p. 36 C), heptane (98 C) and ethane/ethene are the products
    a 25 C start recovers worst; pentane is 0.68x of Frankel 1989's renormalised share and is the
    biggest downward discrepancy in the five-product comparison. **Any inter-paper ratio involving
    pentane or a hydrocarbon should be treated as unusable; ratios among the C6+ oxygenates are much
    safer.**
13. **No aqueous phase, no amine, no pH, no water activity, no matrix, no antioxidant.** Neat
    hydroperoxide in a hot metal port under helium. Contributes nothing to the aldehyde-lysine
    channel, matrix retention, or any pH-dependent term. **And no additive arm at all** — this paper
    has no hydrogen-donor experiment, so it cannot be used as a hold-out companion to Frankel 1989's
    tocopherol columns.
14. **No DOI exists to cite.** See section 0.
15. **A hold-out hygiene note.** `frankel1989_extraction.md` records that its tocopherol columns are
    a declared hold-out. **This paper is a different paper and its numbers are not covered by that
    declaration** — but it shares the FIT-side substrate (autoxidised methyl linoleate
    hydroperoxides) and gives an independent slate for it, so **using Frankel 1981's linoleate
    column to fit anything makes the 1989 zero-additive column partly redundant rather than
    independent.** If the lane refits, that overlap should be declared.
16. **What to request.** (i) Chang, Smouse, Krishnamurthy, Mookherjee & Reddy, *Chem. Ind.* 1926
    (1966), ref (27) — the only source in the chain that claims a mechanism for 2-pentylfuran, and
    the one whose claim this paper's Table I appears to contradict. (ii) Neff, Frankel & Weisleder,
    AOCS New Orleans abstract, May 1981, ref (31) — the hydroperoxy cyclic peroxide identification
    behind Scheme II. (iii) Still, and as in every Frankel-lab dossier: **a time-resolved,
    response-factor-corrected decomposition at two or more temperatures below 180 C**, which is the
    one experiment that would replace the lane's Q10 assumption with a measurement and which no
    paper in this corpus supplies.
