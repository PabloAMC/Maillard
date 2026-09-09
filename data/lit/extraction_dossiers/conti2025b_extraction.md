# Conti 2025b — EXTRACTION (Part II of the same soy-protein-concentrate + 1.5 % thiamine extrusion study: sodium and glutamic acid released into saliva across three mastication times by ion chromatography, plus a full SAFE-GC-MS volatile inventory and GC-O of the two extreme extrudates, 30 % M / 180 C and 38 % M / 140 C)

### The volatile half of the Conti pair: 98 and 71 compounds inventoried in a thiamine-flavoured textured soy protein, with 2-pentylpyridine and eight alkylthiazoles among them — but every "ug/g" in its big table is a peak area scaled by one deuterated hexanal, so nothing in it is a concentration.

**Source on disk:** `data/articles/Conti2025b.pdf` (11 pp., Food Research International 218 (2025)
116938). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Conti2025b.txt`), whole file. Table 1 came through clean. **Table 2 is very
large (98 rows across seven printed blocks with two page-break continuations) and came through in
pieces with the column alignment shifting between blocks**; every row is legible and is re-typed
below, but the block boundaries are reconstructed from the printed sub-headings and one row
(2,6-dimethoxyphenol) is printed twice by the page-break repeat (Flags 6). Figures 1 (sensory
intensities), 2 (sodium and glutamic acid in saliva), 3 (multiple factor analysis) and 4 (principal
component analysis) are images and are **figure-only**. **Supplementary Tables 1, 2 and 3** (panel
characteristics; the sensory intensity values; the sodium and glutamic acid values) are online-only
and are **not on disk** — which matters, because the paper's own numeric release data live there.
Repo status before this dossier: neither Conti paper had a dossier; see the companion
`conti2025_extraction.md` for Part I.

## 0. Identity

| field | value |
|---|---|
| Title | "Oral processing of meat-flavour textured soy proteins — **Part II: Influence on taste compounds release and sensory perception**" |
| Authors | Ana Carolina Conti (corresponding, ac.conti@unesp.br), Chantal Septier, **Karine Gourrat**, Helene Laboure, Christian Salles |
| Affiliations | (a) Sao Paulo State University (Unesp), Ibilce, Sao Jose do Rio Preto, Brazil; (b) Centre des Sciences du Gout et de l'Alimentation, CNRS / INRAE / Institut Agro / Universite Bourgogne Europe, Dijon; (c) **CNRS, INRAE PROBE Research Infrastructure, Chemosens Facility, Dijon** — the analytical platform, and the reason Part II has a GC and Part I does not |
| Venue | **Food Research International 218 (2025) 116938**. Received 18 Sep 2024, revised 13 Jun 2025, accepted 20 Jun 2025, online 21 Jun 2025 |
| DOI | **10.1016/j.foodres.2025.116938** |
| Companion | Part I is `Conti2025.pdf`, FRI 208:116169, cited throughout as "Conti et al. (2025)" — see `conti2025_extraction.md` |
| Same design | Same soy protein concentrate, same 1.5 % thiamine, same three extrusion conditions, same hydration, same salt and MSG, same oil, same eight panellists, same 1/3-2/3-3/3 mastication design. Part I measured the bolus; Part II measures what leaves it |
| Ethics | Inserm Ethics Evaluation Committee No. 20-754bis, March 2021 |
| Funding | FAPESP 2019/20911-4; CNPq 303602/2022-8; Institut Agro Dijon |

**Both Conti papers are now identified.** `Conti2025.pdf` = Part I, FRI 208:116169 (bolus and
texture). `Conti2025b.pdf` = Part II, FRI 218:116938 (tastant release, sensory, and the volatile
inventory). Only Part II contains chemistry.

## 1. Why it matters

Two things in this paper touch `tasks/roadmap_for_scientists.md` section 5d directly, and one of
them touches a part of Programme 7 that no other paper in this cluster reaches.

**(a) Programme 7 part (iii) — the lipid-Maillard cross products.** The roadmap names them:
"The lipid-Maillard cross products (**2-pentylpyridine, the alkylthiazoles**) as rules first, waves
when a rate exists." This paper prints **2-pentylpyridine** in the 30 % M / 180 C extrudate at
0.13 +/- 0.01 (their ug/g units), described at the sniff port as grass/floral by 38 % of assessors,
and absent from the 38 % M / 140 C extrudate. It prints **eight thiazoles**: 2-methylthiazole,
4,5-dimethylthiazole, 2,4,5-trimethylthiazole, 5-ethyl-4-methylthiazole, 5-ethenyl-4-methylthiazole,
2-acetylthiazole, benzothiazole and 4-methyl-5-ethanolthiazole. **This is the only paper on disk in
this cluster that carries both of the roadmap's named cross-product families in one table, from a
real extruded plant-protein product with a stated process.** The catch is severe and is the theme
of this dossier: the numbers are not concentrations (Flags 1).

**(b) The thiamine lane.** The whole product is soy protein concentrate carrying **1.5 % w/w
thiamine** through an extruder, and the authors read the volatile inventory as thiamine degradation
plus lipid oxidation. They name, from Dreher 2003, four compounds "formed from thiamine
degradation" that they find here: 4,5-dimethylthiazole, dihydro-2-methyl-3(2H)-thiophenone,
2-carboxaldehyde-5-methylthiophene and pentan-1-ol; and they trace the mechanism through
**4-methyl-5-(2-hydroxyethyl)thiazole** — printed in Table 2 as 4-methyl-5-ethanolthiazole, and by
some margin the largest single peak in both samples (20.25 and 13.58 in their units). The
repository carries `thiamine_availability` as a modifier and works the sulfur ladder in
`k6a_sulfur_ladders_synthesis.md`; this is a real product-level inventory of what a 1.5 % thiamine
charge produces in a soy concentrate at two extrusion severities.

**(c) A carried-volatile signal, in the wrong units.** Hexanal, 1-octen-3-ol, 2-pentylfuran,
nonanal, hexan-1-ol, pentan-1-ol, (E,E)-nona-2,4-dienal and (E,E)-deca-2,4-dienal are all in the
table, and the authors read them as the soy's own lipid-oxidation load. **They also print the
lipid content of the concentrate: "lower than 3 % (data provided by the supplier)"** — which is the
only lipid figure in the whole Conti pair and puts this concentrate inside the roadmap's assumed
"1 to 3 % lipid" band, unlike Trikusuma's 8 %-fat isolate.

What it does not give: no true concentration of anything volatile; no odour threshold and therefore
no odour-activity value; no time course of any volatile; and, for the taste compounds it does
resolve in time, the actual numbers are in a supplementary table that is not on disk.

## 2. Methods as they matter to a model

- **Material and product: identical to Part I.** Soy protein concentrate **Arcon SM (ADM Foods &
  Wellness), minimum 70 g protein/100 g dry basis**; thiamine hydrochloride (Sigma-Aldrich,
  > 99 %) at **1.5 % w/w** on the concentrate; RXPQ Labor 24 single-screw extruder, compression
  3.3:1, L/D 15.5:1, pre-die 5.8 mm, die 3.6 mm, feed 170 g/min, 216 rpm, zones 1-3 at ~40 / 60 /
  80 C; three conditions **30 % M / 180 C, 34 % M / 160 C, 38 % M / 140 C** (moisture dry basis,
  temperature at zone 5). Hydration: 2 cm pieces in water taken off the boil at **1:4 for 15 min**;
  then **1.0 g salt and 0.4 g monosodium glutamate per 100 g** hydrated; half with **7.0 g
  vegetable oil per 100 g**; six products; frozen at -18 C; served reheated to **50 C in the
  middle**.
- **The one composition number Part I did not print.** Section 3.2: "**the lipid content in the soy
  protein concentrate used in this study was lower than 3 % (data provided by the supplier)**". A
  supplier figure, an upper bound, not a measurement made here.
- **The in-vivo design.** Eight panellists (6 F, 2 M, 23-60), the same panel as Part I, with
  salivary flow 0.72-3.15 mL/min and carrot d50 3.58-7.43 mm. **6 FTSPs x 3 mastication times x 3
  replicates = 54 samples.** Times are 1/3, 2/3 and 3/3 of each panellist's own time-to-swallow for
  that product, which Part I reports as group means of 14.5 to 19.8 s. A saliva sample with no food
  serves as the per-panellist blank and **is subtracted**. After chewing, the bolus is spat into a
  180 mL cup, an aliquot of saliva taken immediately, the bolus moved to a 2 mL tube on ice; all
  saliva centrifuged **6000 g for 20 min**; supernatant held at **-80 C**.
- **Sodium and glutamic acid: ion chromatography with authentic standard curves.** Saliva defrosted,
  **diluted 100-fold** in Milli-Q water, centrifuged 5 min at 10,000 g. Thermo Fisher HPIC.
  *Sodium*: CG16-4UM (4 x 50 mm) guard, CS16-4UM (4 x 250 mm) column, **sulfuric acid 13 N as
  eluent at 0.64 mL/min**, conductivity detection. *Glutamic acid*: CG AminoPac PA10 (2 x 50 mm)
  guard, AminoPac PA10 (2 x 250 mm) column, Milli-Q water / 250 mM NaOH / sodium acetate at
  0.25 mL/min, **pulsed amperometric detection**. Chromeleon 6.8. **Standard curves 0.00625 to
  0.1 g sodium/L and 0.003 to 0.035 g glutamic acid/L.** These two are genuine external-standard
  quantifications with a stated calibration range; the volatiles are not.
- **Volatile extraction: SAFE, and the conditions were changed from the published method.**
  Only two products were extracted — **30 % M / 180 C and 38 % M / 140 C, both WITHOUT vegetable
  oil**. "First, 100 mL of ultrapure water was stirred with **200 uL [175 uL] of an internal
  standard (hexanal D12 marked with deuterium at 109.8 ng/uL)** for 30 min, followed by its
  addition to **20 g (+ 0.4 g) of sample (previously ground while frozen)** and stirring for 10 min
  (always agitated under ice). The mixture was introduced into the SAFE apparatus, and vacuum
  distillation (**2 Pa**) was performed for **90 min at 30 C** [10^-2 Pa for 2 h]." The bracketed
  values are Thomsen et al. 2014's originals: **the vacuum used here is 200 times poorer and the
  distillation is half as long**, though at a slightly lower body temperature. The aqueous
  distillate was extracted three times with 10 mL [15 mL] of distilled dichloromethane (> 99.9 %),
  filtered through glass wool, dried over anhydrous sodium sulphate, concentrated on a
  Kuderna-Danish in a **70 C** water bath to ~250 uL, and made up to **300 uL [500 uL]**; frozen at
  -20 C. **Three extractions per condition**; the internal standard's coefficient of variation over
  the whole procedure was **11 % (30 % M / 180 C)** and **12.3 % (38 % M / 140 C)**.
  **Internal-standard loading (mine): 200 uL x 109.8 ng/uL = 21.96 ug into 20 g of sample =
  1.098 ug/g**, which is the scale against which every printed number in Table 2 is set.
- **GC-MS.** 1 uL injected; Agilent 7890A; **DB-WAX 30 m x 0.25 mm x 0.5 um**; helium 1.2 mL/min;
  40 C then 4 C/min to 240 C, 10 min isotherm; Agilent 5973 MSD, 70 eV, scan 29-350 amu.
  Identification by mass spectrum against Wiley, NIST and an **INRAE internal database built from
  standard compounds**, plus retention indices from a C10-C30 alkane series compared with NIST 2021
  published values. **No authentic standard is injected for the individual analytes as a
  quantification step**; the internal database is an identification aid.
- **Quantification: it is not one, and the paper says so.** Verbatim, section 2.6.2:
  "**Semiquantitative data for each compound were obtained by automatic or manual integration of the
  total ion count peak area (arbitrary units) and expressed as a concentration (ug/g of sample)
  based on the concentration of a hexanal D12 standard of 109.8 ng/uL.**" So every number in
  Table 2 is a total-ion-count peak area divided by the deuterated hexanal's peak area and
  multiplied by 1.098 ug/g — **a single response factor of 1 applied to 98 compounds spanning
  pyrazines, thiazoles, thiophenes, alcohols, aldehydes, ketones, furans and phenols**. Under the
  house rule this is `peak_area_only` and is never a concentration (Flags 1).
- **GC-O: detection frequency, no dilution series.** Agilent 6890A, FID plus a Gerstel olfactometric
  detection port, **Y-splitter 1:1**, same DB-WAX column and same oven program. **Eight assessors
  (18-65, average 45; 7 F, 1 M), a different panel from the sensory one**, sensitive to odour and
  experienced in GC-O. A 40-minute run each; free description in their own words; AcquiSniff
  software with a button and a microphone; Openlab for the chromatography. Odourant zones
  characterised by (i) the number of simultaneous detections, (ii) the panel's description and
  (iii) the polar retention index. **An odourant zone counts as identified only if at least 3 of 8
  assessors describe the same or a similar odour.** The percentages in Table 2's OD column are the
  percentage of assessors who detected that odour.
- **Sensory evaluation.** Salty taste, umami taste and meat aroma scored at each mastication time on
  a **5-point scale** (1 "weak", 3 "moderate", 5 "strong"). **No training, no references given** —
  "All panelists were experienced in sensory evaluation and performant in flavour rating and
  recognition; thus, no particular training was required. Moreover, no references for the attributes
  were provided." Eighteen samples per session.
- **Basis of the numbers.** Sodium and glutamic acid in **g/L of saliva**, blank-subtracted, from a
  100-fold dilution. Volatiles in **ug/g of sample** as defined above (semi-quantitative). Sensory
  in scale points 1-5. Retention indices dimensionless on a wax column.
- **Statistics.** General linear model on the sensory and release data with panellist, extrusion
  condition, oil and mastication time as factors, then Tukey; Statistica 7.0; alpha = 0.05. The two
  volatile samples compared by **Student's t-test on n = 3**, with normality first checked by
  Anderson-Darling, Lilliefors and Jarque-Bera. MFA on sensory plus release; PCA at time 3/3 only,
  adding Part I's texture and mastication time; **KMO = 0.677, 89.9 % cumulative variance on the
  first two axes**; only loadings >= 0.70 or <= -0.70 interpreted.

## 3. Tables re-typed

### Table 1. "Values of p (F-ratio) of factors from general linear model to sensory intensities and compound release from the flavoured textured soy proteins"

Footnote as printed: "p-values lower than the significance level used (0.05) indicate significant
effect of the factor."

| factor | salty taste | umami taste | meat aroma | sodium | glutamic acid |
|---|---|---|---|---|---|
| panelist | 0.000 (29.01) | 0.000 (37.27) | 0.000 (31.14) | 0.000 (24.33) | 0.000 (13.93) |
| extrusion condition (EC) | 0.418 (0.87) | **0.000 (17.45)** | **0.039 (3.29)** | **0.000 (11.63)** | **0.000 (15.58)** |
| vegetable oil (VO) | 0.352 (0.87) | 0.586 (0.30) | 0.090 (2.89) | 0.175 (1.85) | 0.839 (0.04) |
| time (T) | 0.890 (0.12) | 0.200 (1.62) | 0.085 (2.48) | 0.451 (0.80) | **0.002 (6.34)** |
| EC x VO | 0.862 (0.15) | 0.877 (0.13) | 0.902 (0.10) | 0.527 (0.64) | **0.017 (4.15)** |
| EC x T | **0.000 (5.33)** | 0.134 (1.78) | 0.573 (0.73) | **0.011 (3.32)** | **0.023 (2.88)** |
| VO x T | 0.979 (0.02) | 0.324 (1.13) | 0.699 (0.36) | **0.031 (3.53)** | 0.444 (0.81) |

### Table 2. "Volatile compounds (VC) identified on the flavoured textured soy proteins, their quantities and odour description"

Footnotes exactly as printed: 1 = linear retention index (LRI) from NIST 2021; 2 = M (moisture of
the soy protein concentrate, dry basis) before extrusion, temperature at zone 5 of the extruder
barrel; 3 = experimental LRI calculated by GC-MS; 4 = **mean +/- standard deviation (n = 3)**;
5 = odour description by the panellists through GC-O (% of panellists that detected odour).
"Different letters in the same line indicate different statistical means by the Student's t-test
(p < 0.05). The dash (-) means volatile compound not identified, or relative area non-existent, or
odoriferous description non-existent. Trace means that the volatile compound is present in trace
quantities in the sample." `*` = CAS number not found at NIST 2021; `**` = LRI not found on wax
columns at NIST 2021; `***` = the compound has no CAS number in the mass spectra libraries.

**The "ug/g sample" columns are semi-quantitative peak-area ratios against deuterated hexanal, not
concentrations (section 2 and Flags 1).**

#### Pyrazines

| compound | LRI lit | LRI exp, 30 % M / 180 C | ug/g, 30 % M / 180 C | odour (% detection), 30 % M / 180 C | LRI exp, 38 % M / 140 C | ug/g, 38 % M / 140 C | odour, 38 % M / 140 C |
|---|---|---|---|---|---|---|---|
| pyrazine | 1219 | 1220 | 1.73 +/- 0.04 a | – | 1220 | 0.16 +/- 0.02 b | – |
| methylpyrazine | 1276 | 1276 | 7.47 +/- 0.36 a | – | 1276 | 0.62 +/- 0.00 b | – |
| 2,5-dimethylpyrazine | 1333 | 1334 | 5.99 +/- 0.54 a | grilled (50) | 1334 | 0.94 +/- 0.07 b | – |
| 2,6-dimethylpyrazine | 1339 | 1340 | 3.85 +/- 0.33 a | popcorn/biscuit (38) | 1340 | 0.30 +/- 0.02 b | – |
| ethylpyrazine | 1344 | 1346 | 2.49 +/- 0.24 a | peanut (63) | 1346 | 0.28 +/- 0.04 b | grilled (50) |
| 2,3-dimethylpyrazine | 1357 | 1358 | 1.84 +/- 0.18 a | nutty/croissant (38) | 1358 | 0.16 +/- 0.04 b | – |
| 2-ethyl-6-methylpyrazine | 1395 | 1397 | 2.52 +/- 0.26 a | – | 1397 | 0.30 +/- 0.01 b | – |
| 2-ethyl-5-methylpyrazine | 1402 | 1404 | 2.43 +/- 0.25 a | green/cardboard/vegetable (75) | 1404 | 0.31 +/- 0.01 b | green/vegetable (63) |
| trimethylpyrazine | 1421 | 1418 | 9.60 +/- 0.94 | mushroom (38) | – | – | – |
| 2-(n-propyl)pyrazine | 1428 | 1431 | 0.29 +/- 0.13 | – | – | – | – |
| 2,6-diethylpyrazine | 1444 | 1447 | 0.34 +/- 0.06 a | grilled (38) | 1445 | 0.21 +/- 0.02 b | – |
| 3-ethyl-2,5-dimethylpyrazine | 1458 | 1458 | 4.18 +/- 0.41 | – | – | – | – |
| 5-ethyl-2,3-dimethylpyrazine | 1460 | 1475 | 2.32 +/- 0.29 | – | – | – | – |
| tetramethylpyrazine | 1478 | 1489 | 1.11 +/- 0.10 | – | – | – | – |
| 3,5-diethyl-2-methylpyrazine | 1505 | 1507 | 1.01 +/- 0.06 | – | – | – | – |
| 2-acetylpyrazine | 1638 | 1639 | 1.66 +/- 0.15 | – | – | – | – |
| 2-acetyl-6-methylpyrazine | 1688 | 1696 | 0.71 +/- 0.05 | – | – | – | – |
| 2-isobutyl-3-methylpyrazine (printed "pytazine") | 1490 | **1701** | 0.21 +/- 0.01 a | – | **1701** | 0.10 +/- 0.01 b | – |

The abstract and section 3.2 also mention **2,5-dimethyl-3-(3-methylbutyl)pyrazine** as one of the
compounds that differed between the samples and was odour-active in both; **that compound does not
appear anywhere in the transcribed Table 2** (Flags 5).

#### Thiazoles

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| 2-methylthiazole | 1250 | 1245 | 0.06 +/- 0.00 | – | – | – | – |
| **4,5-dimethylthiazole** | 1372 | 1387 | **3.52 +/- 0.31 a** | **sulfuric (50)** | 1387 | 0.39 +/- 0.04 b | – |
| 2,4,5-trimethylthiazole | 1390 | 1392 | 0.11 +/- 0.01 | – | – | – | – |
| 5-ethyl-4-methylthiazole | 1467 | 1451 | 0.46 +/- 0.05 | – | – | – | – |
| 5-ethenyl-4-methylthiazole | 1512 | 1539 | 5.98 +/- 0.50 | – | – | – | – |
| 2-acetylthiazole | 1660 | 1660 | 0.29 +/- 0.03 | – | – | – | – |
| benzothiazole | 1973 | 1975 | 0.24 +/- 0.09 | cooked/solvent (38) | – | – | – |
| **4-methyl-5-ethanolthiazole** (= 4-methyl-5-(2-hydroxyethyl)thiazole) | 2311 | 2329 | **20.25 +/- 2.83 a** | – | 2329 | **13.58 +/- 1.21 b** | – |

#### Pyridines

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | 38 % M / 140 C |
|---|---|---|---|---|---|
| **2-pentylpyridine** | 1572 | 1589 | **0.13 +/- 0.01** | **grass/floral (38)** | not detected |
| 3-methoxypyridine | 1581 | 1595 | 0.18 +/- 0.03 | – | not detected |
| 2-acetylpyridine | 1603 | 1615 | 0.45 +/- 0.04 a | – | 1615; 0.06 +/- 0.01 b; no odour |

#### Pyrrole, oxazole, indole

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | 38 % M / 140 C |
|---|---|---|---|---|---|
| 2-acetylpyrrole | 1983 | 1985 | 1.12 +/- 0.08 a | sweet/grilled/chemical (88) | 1985; 0.25 +/- 0.03 b; no odour |
| 2,4,5-trimethyl-oxazole | 1206 | 1207 | 0.11 +/- 0.01 | – | not detected |
| 1H-indole | 2460 | 2464 | 0.93 +/- 0.07 | – | not detected |

#### Thiophenes

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| 2-methylthiophene | 1100 | 1099 | 0.56 +/- 0.37 | – | – | – | – |
| 2-pentylthiophene | 1460 | 1467 | 0.11 +/- 0.01 a | – | 1467 | 0.09 +/- 0.01 b | – |
| dihydro-2-methyl-3(2H)-thiophenone | 1542 | 1541 | 2.76 +/- 0.34 | – | – | – | – |
| dihydro-3-(2H)-thiophenone | 1570 | 1575 | 0.10 +/- 0.01 | – | – | – | – |
| 2- or 3-carboxaldehydethiophene | 1693 | 1706 | 2.53 +/- 0.22 a | – | 1705 | 0.21 +/- 0.01 b | grilled/grass/floral (50) |
| 2-carboxaldehyde-5-methylthiophene | 1759 | 1741 | 0.43 +/- 0.03 | – | – | – | – |
| 2-propionylthiophene | 1842 | 1852 | 0.27 +/- 0.01 | – | – | – | – |
| 2-ethyl-5-isopentylthiophene | * | 2177 | 0.36 +/- 0.04 a | – | 2177 | 0.14 +/- 0.01 b | – |

#### Other sulfur-containing, and one nitrogen-and-sulfur compound

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | LRI exp 38 % | ug/g 38 % |
|---|---|---|---|---|---|
| 2-methylthiocyclohexanedi-1,3-one | *** | 1900 | 1.36 +/- 0.05 a | 1900 | 0.43 +/- 0.06 b |
| 1,4-dithian-2-one | * | 1944 | 1.74 +/- 0.11 a | 1944 | 0.31 +/- 0.02 b |
| 3,6-dithione-hexahydro-1,2,4,5-tetrazine | * | 1909 | 3.68 +/- 0.22 **b** | 1910 | **4.63 +/- 0.22 a** |

#### Aldehydes

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| pentanal | 984 | < 1000 | 0.49 +/- 0.02 a | – | < 1000 | 0.46 +/- 0.02 a | – |
| **hexanal** | 1088 | 1087 | **1.95 +/- 0.12 b** | fruity (38) | 1087 | **4.70 +/- 0.22 a** | green/floral (50) |
| heptanal | 1194 | 1192 | 0.13 +/- 0.03 b | – | 1191 | 0.22 +/- 0.03 a | – |
| octanal | 1295 | 1296 | 0.24 +/- 0.06 a | – | 1296 | 0.23 +/- 0.06 a | – |
| **nonanal** | 1400 | 1401 | 0.44 +/- 0.04 a | cosmetic/solvent (50) | 1401 | 0.42 +/- 0.04 a | – |
| (E)-oct-2-enal | 1437 | 1439 | 0.14 +/- 0.02 a | – | 1439 | 0.14 +/- 0.01 a | grilled (50) |
| **methional** | 1463 | 1463 | **trace** | potato (50) | 1463 | **trace** | potato (75) |
| decanal | 1507 | – | – | – | 1507 | 0.27 +/- 0.03 | – |
| 2-butyloct-2-enal | 1659 | – | – | – | 1677 | 0.16 +/- 0.01 | **meat/peanut (88)** |
| (E,E)-nona-2,4-dienal | 1712 | – (dash printed) | 0.03 +/- 0.00 b | – | 1712 | 0.21 +/- 0.01 a | – |
| (E,Z)-deca-2,4-dienal | 1770 | – | – | – | 1775 | 0.06 +/- 0.00 | – |
| (1-methylethyl)-4-benzaldehyde | 1794 | – | – | – | 1794 | 0.08 +/- 0.00 | – |
| (E,E)-deca-2,4-dienal | 1824 | 1822 | 0.38 +/- 0.03 b | – | 1822 | 0.70 +/- 0.03 a | fruit/floral (88) |

#### Ketones

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| butan-2,3-dione (diacetyl) | 980 | < 1000 | 0.28 +/- 0.02 a | butter (63) | < 1000 | 0.04 +/- 0.00 b | butter (38) |
| pentan-2,3-dione | 1060 | 1059 | 0.18 +/- 0.01 | – | – | – | – |
| pent-3-en-2-one | 1132 | 1131 | 0.04 +/- 0.00 | – | – | – | – |
| heptan-2-one | 1184 | 1189 | 0.97 +/- 0.09 a | – | 1189 | 0.32 +/- 0.03 b | – |
| cyclopentanone | 1187 | 1193 | 0.70 +/- 0.69 | – | – | – | – |
| 3-hydroxybutan-2-one (acetoin) | 1291 | 1293 | 0.97 +/- 0.10 a | – | 1293 | 0.27 +/- 0.01 b | – |
| oct-1-en-3-one | 1306 | 1307 | 0.46 +/- 0.04 a | mushroom/vegetable (88) | 1307 | 0.03 +/- 0.01 b | mushroom (88) |
| octa-2,3-dione | 1325 | 1329 | 0.09 +/- 0.03 a | – | 1329 | 0.21 +/- 0.08 a | – |
| 3-hydroxypentan-2-one | 1344 | 1352 | 0.63 +/- 0.04 | – | – | – | – |
| 2-cyclopenten-1-one | 1383 | 1367 | 0.71 +/- 0.06 | – | – | – | – |
| 4-hydroxy-4-methylpentan-2-one | 1376 | 1372 | 0.11 +/- 0.01 | – | – | – | – |
| 2-methyl-2-cyclopenten-1-one | 1395 | 1381 | 0.19 +/- 0.02 | – | – | – | – |
| 1-hydroxybutan-2-one | 1381 | 1382 | 0.10 +/- 0.01 a | – | 1382 | 0.02 +/- 0.01 b | floral (38) |
| 1-(acetyloxy)propan-2-one | 1470 | 1470 | 3.43 +/- 0.30 | – | – | – | – |
| (2-octenyl)-cyclopentan-2-one | *** | – | – | – | 1477 | 0.11 +/- 0.01 | – |
| decan-2-one | 1503 | 1502 | 0.28 +/- 0.03 a | – | 1502 | 0.18 +/- 0.01 b | – |
| non-3-en-2-one | 1518 | – | – | – | 1522 | 0.08 +/- 0.00 | – |
| (E,Z)-octa-3,5-dien-2-one | 1534 | 1534 | 3.18 +/- 0.29 a | musty/green/vegetable (63) | 1534 | 1.78 +/- 0.05 b | vegetable (75) |
| (E,E)-octa-3,5-dien-2-one | 1585 | – | – | – | 1581 | 0.21 +/- 0.02 | vegetable/plastic (38) |
| butyrolactone | 1640 | – | – | – | 1640 | 0.12 +/- 0.02 | – |
| 2-hydroxy-3-methylcyclopent-2-en-1-one | 1837 | 1843 | 0.31 +/- 0.00 | grilled/grass/spice (88) | – | – | – |
| (E)-6,10-dimethylundecadi-5,9-en-2-one | 1865 | 1864 | 1.17 +/- 0.04 a | fruity/floral/chemical (88) | 1864 | 0.17 +/- 0.01 b | fruity/floral/chemical (63) |

#### Alcohols and phenols

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| butan-1-ol | 1150 | – | – | – | 1149 | 0.06 +/- 0.01 | – |
| pent-1-en-3-ol | 1164 | 1164 | 0.09 +/- 0.01 b | – | 1164 | 0.13 +/- 0.01 a | – |
| pentan-1-ol | 1256 | 1258 | 0.76 +/- 0.03 b | – | 1258 | 1.63 +/- 0.03 a | – |
| **hexan-1-ol** | 1363 | 1362 | 0.45 +/- 0.07 a | green (38) | 1362 | 0.46 +/- 0.02 a | grilled (50) |
| **oct-1-en-3-ol** | 1458 | – | – | – | 1458 | **2.62 +/- 0.06** | cardboard/green/nutty (50) |
| heptan-1-ol | 1465 | 1465 | 0.16 +/- 0.03 b | vegetable (38) | 1465 | 0.35 +/- 0.03 a | potato (75) |
| 2-ethylhexan-1-ol | 1499 | 1498 | 0.23 +/- 0.07 a | musty/solvent (38) | 1498 | 0.27 +/- 0.06 a | musty/green (75) |
| octan-1-ol | 1568 | – | – | – | 1568 | 0.55 +/- 0.02 | grass/floral/fruity (63) |
| 2-methyldecan-2-ol | ** | – | – | – | 1593 | 0.13 +/- 0.02 | – |
| (E)-oct-2-en-1-ol | 1626 | 1626 | 0.38 +/- 0.04 a | – | 1626 | 0.41 +/- 0.01 a | – |
| undecan-6-ol or non-1-en-4-ol | 1640 | 1658 | 0.30 +/- 0.03 b | – | 1658 | 0.41 +/- 0.02 a | – |
| 2-methoxyphenol (guaiacol) | 1872 | 1872 | 4.60 +/- 0.37 | – | – | – | – |
| phenylmethanol (benzyl alcohol) | 1893 | 1890 | 0.54 +/- 0.03 a | – | 1890 | 0.44 +/- 0.02 b | – |
| 2-phenylethanol | 1927 | 1926 | 0.17 +/- 0.13 a | – | 1926 | 0.17 +/- 0.13 a | – |
| phenol | 2018 | 2017 | 0.54 +/- 0.10 a | sweet/floral/solvent (38) | 2017 | 0.08 +/- 0.00 b | sweet/floral/solvent (50) |
| **2-methoxy-4-vinylphenol** | 2210 | 2209 | 5.87 +/- 0.45 a | floral (63) | 2209 | 0.26 +/- 0.08 b | floral/soap (75) |
| 2,6-dimethoxyphenol | 2273 | 2280 | 0.78 +/- 0.07 | – | 2280 | trace | – |
| hexadecan-1-ol | 2385 | 2388 | 0.29 +/- 0.02 | – | – | – | – |

#### Esters, furans, pyrone, unidentified

| compound | LRI lit | LRI exp 30 % | ug/g 30 % | odour 30 % | LRI exp 38 % | ug/g 38 % | odour 38 % |
|---|---|---|---|---|---|---|---|
| pentylhexanoate | 1525 | 1519 | 0.05 +/- 0.01 | – | – | – | – |
| methyl-4-oxopentanoate | 1560 | 1579 | 0.13 +/- 0.01 | – | – | – | – |
| **2-pentylfuran** | 1235 | 1237 | 0.26 +/- 0.05 a | – | 1237 | 0.24 +/- 0.11 a | – |
| 2-methyl-dihydrofuran-3(2H)-one | 1270 | 1272 | 2.63 +/- 1.71 a | – | 1272 | 0.49 +/- 0.01 a | – |
| **2-acetylfuran** | 1514 | 1514 | 1.31 +/- 0.13 a | – | 1514 | 0.12 +/- 0.01 b | – |
| dihydro-5-methylfuran-2(3H)-one | 1619 | 1623 | 0.27 +/- 0.02 | – | – | – | – |
| 2-methanolfuran (furfuryl alcohol) | 1670 | 1670 | 6.34 +/- 0.29 a | **meat/sulfuric/solvent (50)** | 1672 | 2.04 +/- 0.11 b | **meat (88)** |
| gamma-hexalactone | 1715 | 1715 | 0.32 +/- 0.04 a | – | 1715 | 0.36 +/- 0.02 a | – |
| 5-methyl-2-furfuryl alcohol | 1729 | 1732 | 0.75 +/- 0.08 | – | – | trace | – |
| furan-2(5H)-one | 1767 | 1766 | 0.25 +/- 0.05 a | – | 1766 | 0.16 +/- 0.01 b | – |
| 5-pentyl-5(H)-furan-2-one | 2076 | 2091 | 0.55 +/- 0.03 a | – | 2091 | 0.58 +/- 0.05 a | – |
| **maltol** | 1991 | 1991 | 3.68 +/- 0.37 a | – | 1998 | 1.81 +/- 0.23 b | – |
| unidentified compound | – | 1783 | 0.38 +/- 0.03 | – | – | – | – |
| unidentified compound | – | 1788 | 1.15 +/- 0.08 | mushroom (38) | – | – | – |

### Numbers printed only in the running text

| quantity | value | where |
|---|---|---|
| compounds identified | **98** in the 30 % M / 180 C extrudate, **71** in the 38 % M / 140 C extrudate | section 3.2 |
| odourant zones detected by GC-O | **50** (30 % M / 180 C) and **38** (38 % M / 140 C), at detection frequency > 3/8 | section 3.2 |
| odourant zones matched to a GC-MS compound | **26** and **22** respectively | section 3.2 |
| lipid content of the soy protein concentrate | **< 3 %** ("data provided by the supplier") | section 3.2 |
| internal-standard reproducibility over the whole SAFE procedure | CV **11 %** (30 % M / 180 C) and **12.3 %** (38 % M / 140 C) | section 2.6.1 |
| sodium in saliva after mastication | **4.10 +/- 1.76 to 6.91 +/- 5.57 g/L** (range over the six products x three times) | section 3.1, full values in Supplementary Table 3 (**not on disk**) |
| glutamic acid in saliva | **0.46 +/- 0.33 to 1.49 +/- 1.19 g/L** | section 3.1, Supplementary Table 3 (**not on disk**) |
| salty taste intensity | 2.3 +/- 1.0 to 3.0 +/- 1.2 | section 3.1, Supplementary Table 2 (**not on disk**) |
| umami taste intensity | 2.6 +/- 1.1 to 3.6 +/- 1.1 | as above |
| meat aroma intensity | 2.5 +/- 1.0 to 3.2 +/- 1.2 | as above |
| sodium release, significant contrasts | at 1/3, higher for 38 % M / 140 C with oil than 34 % M / 160 C without oil; at 2/3, higher for 38 % M / 140 C without oil than 30 % M / 180 C with oil; **at 3/3 no significant difference between any samples** | section 3.1, Fig. 2A |
| glutamic acid release, significant contrasts | **only at 1/3**: 38 % M / 140 C with oil higher than 30 % M / 180 C without oil, 30 % M / 180 C with oil, and 34 % M / 160 C without oil | section 3.1, Fig. 2B |
| PCA quality | KMO 0.677; **89.9 %** cumulative variance on the first two axes | section 2.7 |
| salivary flow of the panel | 0.72 to 3.15 mL/min | section 2.3 (from Part I) |
| carrot d50 of the panel | 3.58 to 7.43 mm | section 2.3 (from Part I) |
| mastication times used | group means 14.5 to 19.8 s (from Part I) | section 3.1 |
| compounds the authors attribute to thiamine degradation and find here | 4,5-dimethylthiazole (skunky/earthy), 2-methyl-4,5-dihydro-3(2H)-thiophenone (sour-fruity/musty/green; printed in Table 2 as dihydro-2-methyl-3(2H)-thiophenone), 2-formyl-5-methylthiophene (meaty; printed as 2-carboxaldehyde-5-methylthiophene), pentan-1-ol (fruity/green) — after Dreher et al. 2003 | section 3.2 |
| the first thiamine degradation product | 4-methyl-5-(2-hydroxyethyl)thiazole, "one of the first compounds generated from thiamine degradation, subsequently forming thiazoles and other sulfur compounds, such as **5-hydroxy-3-mercaptopentan-2-one**, and this one yields other sulfur-containing compounds (such as thiophenes) and furans"; **5-hydroxy-3-mercaptopentan-2-one was not found** | section 3.2 |
| compounds the authors read as off-flavour markers of soy | (E,E)-nona-2,4-dienal and oct-1-en-3-ol "found only in the FTSPs at 38 % M/140 C"; greater hexanal, heptanal, (E,E)-deca-2,4-dienal, pentan-1-ol and pent-1-en-3-ol in the same sample | section 3.2 |
| compounds higher at the severe condition | octa-3,5-dien-2-one, heptan-2-one, oct-1-en-3-one | section 3.2 |
| the literature contrast the authors note | Maga & Kim 1989 found the **least** severe conditions retained the most volatiles; here the **most** severe condition gives more, which the authors attribute to thiamine degradation | section 3.2 |

**Figure-only:** Figures 1A-C (salty, umami and meat aroma against extrusion condition x time),
2A-B (sodium and glutamic acid in saliva against condition, oil and time), 3A-B (MFA) and 4A-B
(PCA). The paper's actual release and intensity values live in Supplementary Tables 2 and 3, which
are **not on disk**, so the only numbers available are the ranges quoted above.

### Arithmetic on the printed numbers (all mine)

1. **The internal-standard level sets the scale.** 200 uL x 109.8 ng/uL = 21.96 ug added to 20 g of
   sample = **1.098 ug/g**. So a table entry of 1.098 means "the same total-ion-count peak area as
   the deuterated hexanal", 20.25 means "18.4 times that area", and 0.03 means "1/37 of it".
   **Reading the table this way is the only defensible reading**, because no compound-specific
   response factor was determined.
2. **The severe condition dominates the nitrogen and sulfur heterocycles.** Comparing the two
   samples where both are present: methylpyrazine **12.0x**, 2,5-dimethylpyrazine **6.4x**,
   2,6-dimethylpyrazine **12.8x**, pyrazine **10.8x**, ethylpyrazine **8.9x**,
   2,3-dimethylpyrazine **11.5x**, 2-ethyl-6-methylpyrazine **8.4x**, 2-ethyl-5-methylpyrazine
   **7.8x**, 4,5-dimethylthiazole **9.0x**, 2-acetylpyrrole **4.5x**, 2-acetylpyridine **7.5x**,
   2-methanolfuran **3.1x**, 2-acetylfuran **10.9x**, 2-methoxy-4-vinylphenol **22.6x**, phenol
   **6.8x**, maltol **2.0x**, 4-methyl-5-ethanolthiazole **1.5x**. Ten pyrazines, five thiazoles,
   two pyridines, six thiophenes, guaiacol, indole and several ketones are **present at the severe
   condition and absent at the mild one**.
3. **The lipid-oxidation markers go the other way.** hexanal **0.41x** (i.e. 2.4x higher at the mild
   condition), heptanal 0.59x, (E,E)-nona-2,4-dienal 0.14x, (E,E)-deca-2,4-dienal 0.54x, pentan-1-ol
   0.47x, pent-1-en-3-ol 0.69x, heptan-1-ol 0.46x; **oct-1-en-3-ol, octan-1-ol, decanal,
   (E,Z)-deca-2,4-dienal and (E,E)-octa-3,5-dien-2-one appear only at the mild condition.** Against
   that, oct-1-en-3-one is **15x higher** at the severe condition and (E,Z)-octa-3,5-dien-2-one
   1.8x higher, so the lipid set does not move as one block. **2-Pentylfuran is statistically
   unchanged** (0.26 vs 0.24, letters both "a"), and so are nonanal, octanal, (E)-oct-2-enal,
   hexan-1-ol, 2-ethylhexan-1-ol, (E)-oct-2-en-1-ol, pentanal, octa-2,3-dione, gamma-hexalactone,
   5-pentyl-5(H)-furan-2-one, 2-methyl-dihydrofuran-3(2H)-one and 2-phenylethanol.
4. **One compound is higher at the mild condition among the heterocycles**:
   3,6-dithione-hexahydro-1,2,4,5-tetrazine, 4.63 vs 3.68 (letters a and b respectively), a factor
   0.79. It is also a compound with no CAS number found at NIST.
5. **A cross-check against the flour and isolate numbers, in the units the paper gives.** If the
   table's hexanal figures were taken at face value as concentrations, 1.95 and 4.70 ug/g would be
   **1950 and 4700 ug/kg of product**, against Bi 2020's raw pea flour at 1260 ug/kg and Trikusuma's
   pea isolate beverage at 331 ug/L. The order of magnitude is right, which is mildly reassuring
   about the scaling but is **not** evidence that the response factor of 1 is adequate for the other
   97 compounds — hexanal is the one compound for which the deuterated internal standard is an
   exact analogue, so hexanal is precisely the row where the semi-quantitative method is at its
   best and every other row at its worst.
6. **The GC-O yield is low.** 50 odourant zones detected but only 26 matched to a GC-MS compound at
   the severe condition (52 %); 38 and 22 at the mild condition (58 %). Nearly half the odour the
   panel smelled is unassigned, and the authors say so: "Some OZs remained unidentified because of
   the presence of compounds that were not detected through GC-MS or because there was no
   information in the literature to link with the GC-MS data."
7. **A retention-index outlier.** 2-isobutyl-3-methylpyrazine is listed with a literature index of
   **1490** and an experimental index of **1701** on the same wax column, a gap of **211 units**.
   Every other row in the table agrees within about 30 units. **That assignment should be treated as
   unsafe** (Flags 5).

## 4. Numbers the repository can use

All volatile rows share: soy protein concentrate Arcon SM (>= 70 % protein dry basis, **< 3 %
lipid**) with 1.5 % w/w thiamine hydrochloride, single-screw extruded, hydrated 15 min at 1:4 in
off-boil water, salted 1.0 g/100 g and MSG 0.4 g/100 g, **no vegetable oil**, frozen, ground frozen,
20 g extracted by SAFE (2 Pa, 90 min, 30 C) into dichloromethane, n = 3 extractions, GC-MS on
DB-WAX. Only the two extreme extrusion conditions were extracted.

| compound | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| **every compound in Table 2, all 98 rows** | see the tables above | printed as "ug/g sample"; **actually a total-ion-count peak area divided by that of deuterated hexanal and multiplied by 1.098 ug/g** | FTSP without oil, 30 % M / 180 C and 38 % M / 140 C | Table 2 pp. 9-11 | **peak_area_only — none of these is a concentration** |
| hexanal | 1.95 +/- 0.12 (severe) and 4.70 +/- 0.22 (mild) | as above | as above | Table 2 | peak_area_only (**but the internal standard IS deuterated hexanal, so this one row's response factor is essentially exact** — see Flags 2) |
| 2-pentylfuran | 0.26 +/- 0.05 and 0.24 +/- 0.11, not significantly different | as above | as above | Table 2 | peak_area_only |
| 1-octen-3-ol (oct-1-en-3-ol) | not detected (severe) and 2.62 +/- 0.06 (mild) | as above | as above | Table 2 | peak_area_only |
| nonanal | 0.44 +/- 0.04 and 0.42 +/- 0.04, not significantly different | as above | as above | Table 2 | peak_area_only |
| 1-hexanol (hexan-1-ol) | 0.45 +/- 0.07 and 0.46 +/- 0.02, not significantly different | as above | as above | Table 2 | peak_area_only |
| **2-pentylpyridine** | 0.13 +/- 0.01 (severe); not detected (mild) | as above | as above | Table 2, pyridines block | **peak_area_only — the roadmap's named lipid-Maillard cross product, present and odour-active** |
| **the eight thiazoles** | 2-methylthiazole 0.06; 4,5-dimethylthiazole 3.52 / 0.39; 2,4,5-trimethylthiazole 0.11; 5-ethyl-4-methylthiazole 0.46; 5-ethenyl-4-methylthiazole 5.98; 2-acetylthiazole 0.29; benzothiazole 0.24; 4-methyl-5-ethanolthiazole 20.25 / 13.58 | as above | as above | Table 2, thiazoles block | **peak_area_only — the roadmap's named alkylthiazole family** |
| methional | **trace** in both, odour-active (potato) at 50 % and 75 % detection | — | as above | Table 2 | level_only (an explicit "trace", no value) |
| the pyrazine set (18 compounds) | see the pyrazines block | as above | as above | Table 2 | peak_area_only |
| the thiophene set (8 compounds) | see the thiophenes block | as above | as above | Table 2 | peak_area_only |
| maltol | 3.68 +/- 0.37 and 1.81 +/- 0.23 | as above | as above | Table 2 | peak_area_only |
| 2-acetylfuran, furfuryl alcohol, acetoin, diacetyl, 4-vinylguaiacol, guaiacol, benzaldehyde derivatives | see the tables | as above | as above | Table 2 | peak_area_only |
| severe-vs-mild ratios, all shared compounds | see section 3 arithmetic 2, 3 and 4 | dimensionless | one extrusion severity against another | derived from Table 2 (mine) | within_study_ratio (**valid even though the levels are not** — a ratio of the same compound's peak area between two runs cancels the unknown response factor, provided the extraction recovery is equal, which the 11 / 12.3 % internal-standard CVs partly support) |
| GC-O odour descriptions and detection frequencies | see the odour columns above | % of 8 assessors | free description, threshold 3/8 | Table 2 | level_only (**no dilution series, no flavour dilution factor, no threshold, so no odour activity can be computed**) |
| odourant zones | 50 and 38 detected; 26 and 22 assigned | counts | as above | section 3.2 | level_only |
| lipid content of the concentrate | **< 3 %** | % w/w | Arcon SM, supplier data | section 3.2 | level_only (a supplier bound, not measured here) |
| internal standard loading | 21.96 ug into 20 g = 1.098 ug/g | ug/g | hexanal-d12 | derived from section 2.6.1 (mine) | derived (arithmetic on printed quantities) |
| SAFE recovery reproducibility | CV 11 % and 12.3 % | % on the internal standard | whole procedure, n = 3 | section 2.6.1 | measured_level |
| sodium in saliva | 4.10 +/- 1.76 to 6.91 +/- 5.57 | g/L of saliva, blank-subtracted | six products x three mastication times, n = 24 | section 3.1; per-sample values in Supplementary Table 3 (**not on disk**) | measured_level (range only) |
| glutamic acid in saliva | 0.46 +/- 0.33 to 1.49 +/- 1.19 | g/L of saliva, blank-subtracted | as above | as above | measured_level (range only) |
| salty / umami / meat aroma intensities | 2.3-3.0 / 2.6-3.6 / 2.5-3.2 | 5-point scale | as above | section 3.1; Supplementary Table 2 (**not on disk**) | level_only (range only) |
| significance of every factor on the five responses | see Table 1 above | p (F) | GLM, alpha = 0.05 | Table 1 p. 4 | level_only |
| release curves and multivariate maps | — | — | — | Figs. 1, 2, 3, 4 | **figure_only** |

### What is resolved in time, and what is not

**The volatiles are not.** Each extrudate was extracted once per replicate, three replicates, one
state. There is no time inside the extruder (no residence time is printed in either Conti paper),
no storage series and no time course of any volatile.

**The tastants are, coarsely.** Sodium and glutamic acid in saliva are measured at 1/3, 2/3 and 3/3
of each panellist's own chewing time, which is a real in-mouth release axis spanning roughly 5 to
20 s. Glutamic acid shows a significant time effect (p = 0.002) and a condition-by-time interaction
(p = 0.023); sodium shows condition-by-time (p = 0.011) and oil-by-time (p = 0.031) interactions but
no main time effect. **The values behind those tests are in Supplementary Table 3, which is not on
disk**, so this dossier can record that the effect exists and its direction, not its magnitude.

## 5. Flags

1. **Table 2's "ug/g" is not a concentration and must never be entered as one.** The paper's own
   Methods call the data "semiquantitative" and define them as total-ion-count peak areas scaled by
   one deuterated hexanal at 1.098 ug/g of sample. A single response factor of 1 is applied across
   pyrazines, thiazoles, thiophenes, phenols, furans, alcohols, aldehydes and ketones, whose GC-MS
   total-ion responses differ by an order of magnitude or more. Under the house rule every one of
   the 98 rows is `peak_area_only`. **The within-sample ratios are the usable quantity**, because
   the unknown response factor cancels between the two runs of the same compound.
2. **Hexanal is the one exception worth naming, and only partly.** The internal standard is
   **hexanal-d12**, so hexanal's response factor against it is close to exact (isotope effects
   aside) and its printed values, 1.95 and 4.70 ug/g, are much better grounded than the rest of the
   table. They are still not a validated quantification — there is no calibration curve, no
   recovery check on hexanal itself, and the SAFE was run at 200x poorer vacuum than the method it
   cites — but of all 98 rows this is the one where a cautious reader might treat the number as an
   estimate rather than an area. **Even so it is recorded as `peak_area_only` above**, because the
   paper never claims otherwise.
3. **The SAFE conditions were substantially weakened from the cited method and the paper flags this
   itself with brackets.** 2 Pa instead of 10^-2 Pa (**200x poorer vacuum**) and 90 min instead of
   2 h, at 30 C instead of the original temperature. A poorer vacuum recovers less of the
   less-volatile fraction; the high-boiling rows (4-methyl-5-ethanolthiazole at LRI 2329, indole at
   2464, hexadecan-1-ol at 2388) are the most affected, and those are among the largest numbers in
   the table. The 11 / 12.3 % internal-standard CVs measure **reproducibility, not recovery**.
4. **Only two of the six products were analysed for volatiles, and neither had oil.** The 34 % M /
   160 C condition — the one Milani et al. 2022 called optimal, and the one on which the whole
   product design rests — was **not extracted**. Neither was any oiled sample. So the table says
   nothing about the intermediate condition and nothing about what 7 % vegetable oil does to the
   volatile profile, which for a lipid-oxidation question is the obvious missing arm.
5. **Three internal inconsistencies, one of them substantive.** (i) **The meat-aroma result is
   stated three different ways.** Section 2.6.1 says meat aroma was "more intense for FTSPs
   obtained under 38 % M/140 C compared to 30 % M/180 C (data not shown)"; section 3.1 says
   "Fig. 1C shows that neither the duration of mastication nor the extrusion conditions had an
   effect on the aroma of the meat"; the abstract and section 3.2 say the 30 % M / 180 C products
   "stood out for their salty taste, umami taste, and **meat aroma**". Table 1 gives meat aroma a
   significant extrusion-condition effect (p = 0.039). **These cannot all be true**, and the
   direction of the meat-aroma effect is the paper's headline claim. (ii)
   **2,5-dimethyl-3-(3-methylbutyl)pyrazine** is named twice in section 3.2 as a compound that
   differed significantly between samples and was odour-active in both, but **it does not appear in
   Table 2 at all**. (iii) Section 3.2 says "(E,E)-nona-2,4-dienal and oct-1-en-3-ol were found
   **only** in the FTSPs at 38 % M/140 C", but Table 2 gives (E,E)-nona-2,4-dienal a value of
   **0.03 +/- 0.00 b at 30 % M / 180 C**. Separately, section 3.2 lists "almost all the aldehydes,
   namely, heptan-2-one, decan-2-one, (E,E)-octa-3,5-dien-2-one, hexan-1-ol and oct-1-en-3-ol",
   none of which is an aldehyde.
6. **Table 2 has transcription hazards.** It spans three printed pages with two "(continued)"
   repeats; **2,6-dimethoxyphenol is printed twice with identical values** (once in the alcohols
   block, once at the head of the pyrone continuation), and the "other nitrogen and sulfur-
   containing compound" sub-heading is repeated over the aldehyde and ketone blocks. The
   (E,E)-nona-2,4-dienal row has a dash in the 30 % experimental-LRI cell but a value in the
   quantity cell. The retention index of 2-isobutyl-3-methylpyrazine differs from its literature
   value by **211 units**, which is far outside every other row's agreement and marks that
   identification as unsafe. Compound names in the table are non-standard throughout
   (4-methyl-5-ethanolthiazole for 4-methyl-5-(2-hydroxyethyl)thiazole; 2-methanolfuran for
   furfuryl alcohol; 2-carboxaldehydethiophene for thiophene-2-carbaldehyde), and one is a
   typographical error ("2-isobutyl-3-methylpytazine").
7. **No odour thresholds and no dilution series, so no odour activity.** The GC-O uses detection
   frequency only. There is no aroma extract dilution analysis, no flavour dilution factor and no
   threshold anywhere in the paper. The odour descriptions are qualitative and the percentages are
   panel agreement, not potency. Nothing in this paper can be ranked by odour activity.
8. **What this paper does NOT contain, and what to request.** No true concentration of any volatile;
   no measurement of thiamine remaining after extrusion (so the conversion of the 1.5 % charge is
   unknown); no fatty acid profile; no lipoxygenase status of the concentrate; no water activity;
   no extruder residence time, die pressure or melt temperature; no volatile time course; no
   analysis of the intermediate condition or of any oiled sample; no reported limit of detection.
   **To request:** (a) **Supplementary Tables 1, 2 and 3**, which hold the panel characteristics and
   the actual sodium, glutamic acid and sensory intensity values that section 3.1 only summarises as
   ranges — these are the paper's real quantitative results and they are not on disk; (b) the raw
   peak areas and the identity of the GC-MS response used, so a proper response-factor correction
   could be attempted; (c) **Milani, Menis-Henrique & Conti (2022)** and **Milani & Conti (2024)**,
   the two papers that developed this product and are cited by both Conti papers, neither on disk;
   (d) an extraction of the 34 % M / 160 C condition and of an oiled sample.
9. **Registry gaps against `data/keys/compounds.yml`.** Present and keyable: `hexanal`, `nonanal`,
   `heptanal`, `1_hexanol`, `1_octen_3_ol`, `2_pentylfuran`, `methional`, `2_acetylfuran`,
   `2_methylthiophene`, `acetoin` (3-hydroxybutan-2-one), `2_3_butanedione` (diacetyl),
   `4_vinylguaiacol` (2-methoxy-4-vinylphenol), `methylpyrazine`, `2_5_dimethylpyrazine`,
   `2_6_dimethylpyrazine`, `2_3_dimethylpyrazine`, `2_ethylpyrazine`, `trimethylpyrazine`,
   `tetramethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine` (the table's 3-ethyl-2,5-dimethylpyrazine
   is a different isomer — check before mapping), the class id `pyrazines`, and
   `2_methyltetrahydrofuran_3_one` (probably the table's 2-methyl-dihydrofuran-3(2H)-one).
   **Absent and worth adding, because the roadmap names them:** **2-pentylpyridine** — Programme 7
   part (iii) names it explicitly and the registry has no id for it; and **the alkylthiazoles found
   here** — the registry carries `2_hexyl_4_methylthiazole`, `2_pentyl_4_methylthiazole` and
   `4_5_dihydro_2_methylthiazole`, none of which is any of this paper's eight
   (2-methylthiazole, 4,5-dimethylthiazole, 2,4,5-trimethylthiazole, 5-ethyl-4-methylthiazole,
   5-ethenyl-4-methylthiazole, 2-acetylthiazole, benzothiazole,
   4-methyl-5-(2-hydroxyethyl)thiazole). **Also absent:** thiamine as a molecule (the registry has
   only the `thiamine_availability` modifier), maltol, guaiacol, pentan-1-ol, 1-octen-3-one,
   2-acetylpyrrole, 2-acetylpyrazine, 2-acetylpyridine, 2-acetylthiazole, furfuryl alcohol,
   (E,E)-nona-2,4-dienal, (E,E)-deca-2,4-dienal, (E)-oct-2-enal, the octadienones and
   3,6-dithione-hexahydro-1,2,4,5-tetrazine.
10. **Registry gaps against `data/species/off_flavour_targets.yml`.** That file's six compounds are
    hexanal, nonanal, 1-octen-3-ol, 2-pentylfuran, 1-hexanol and furfural. **Five of the six appear
    in this paper's Table 2** (all but furfural, which the paper does not list, though it lists
    furfuryl alcohol and several other furans). The file names the mitigation for 1-octen-3-ol as
    "lipoxygenase inactivation (blanching, ultrasound)"; this paper offers a different lever —
    **extrusion severity** — with 1-octen-3-ol present at the mild condition and absent at the
    severe one, and hexanal 2.4x lower at the severe condition. That is directionally useful
    evidence for the file's `mitigation` field, but it is peak-area evidence and cannot set a
    number. The file has no entry for 2-pentylpyridine or for any thiazole, and no matrix entry for
    textured vegetable protein.
