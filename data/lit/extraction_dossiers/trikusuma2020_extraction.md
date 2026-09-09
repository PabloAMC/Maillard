# Trikusuma 2020 — EXTRACTION (commercial pea protein ISOLATE, 3 % w/w beverage with carrageenan, pH 7.1; unprocessed control vs indirect-tubular UHT at 140 C for 6 s vs the same aged 7 weeks at 5 C; 21 odorants quantified by dynamic headspace GC/MS-QQQ with standard addition, plus a stable-isotope dilution assay for 2-acetyl-1-pyrroline)

### The one paper in this cluster that measures a protein ISOLATE, and it prints all four volatiles the roadmap names — hexanal, 2-pentylfuran, 1-octen-3-ol and a methoxypyrazine — as real concentrations in an unheated pea-isolate beverage, and again after a stated heat treatment and again after seven weeks of cold storage.

**Source on disk:** `data/articles/trikusuma2020.pdf` (8 pp., Food Chemistry 312 (2020) 126082).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/trikusuma2020.txt`), whole file. Tables 1, 2, 3 and 4 came through clean and
are re-typed in full below. Figure 1 (proposed IBMP formation scheme) and Figure 2a/2b (sensory
spider plots for the three beverages and for the two recombination models) are images and are
**figure-only**. **Supplementary Table S1** (the five-point standard-addition curves) exists online
and is **not on disk**. Repo status before this dossier: the benchmark
`data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json` is built from this paper's Table 2 and
its `content_verification` block quotes it, but this paper had **no extraction dossier**.

## 0. Identity — and yes, this IS the benchmark's paper

| field | value |
|---|---|
| Title | "Identification of aroma compounds in pea protein UHT beverages" |
| Authors | Mariana Trikusuma, Laurianne Paravisini, Devin G. Peterson (corresponding, peterson.892@osu.edu) — Department of Food Science and Technology, The Ohio State University, Columbus, OH |
| Venue | **Food Chemistry 312 (2020) 126082**. Received 8 Sep 2019, revised 13 Dec 2019, accepted 17 Dec 2019, online 26 Dec 2019 |
| DOI | **10.1016/j.foodchem.2019.126082** |
| Ethics | OSU IRB #2017H0072 |
| Funding | Flavor Research and Education Center, The Ohio State University |

**Identity check against the benchmark, item by item.** The bundle
`data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json` declares `source_doi`
`10.1016/j.foodchem.2019.126082`. That is the DOI printed on page 1 of this PDF. **This is that
work, not a companion.** Nothing else by this group is on disk.

| benchmark field | benchmark value | what this paper prints | verdict |
|---|---|---|---|
| `source_doi` | 10.1016/j.foodchem.2019.126082 | same, p. 1 | **matches** |
| `conditions.temp_C` | 140 | "processed at a final temperature and pressure of 140 C and 0.55 MPa for 6 s" (section 2.2) | **matches** |
| `conditions.time_min` | 0.1 | 6 s = 0.1 min | **matches** |
| `conditions.ph` | 7.1 | "The pH of the control beverage was measured at 7.1, and there were no significant pH changes observed on the UHT processed and aged samples" | **matches** |
| `precursors."Pea Protein Isolate".concentration_mM` | 1000.0 | "Pea protein isolate (3 % w/w)" — the paper states a mass fraction, **not a molarity**; 1000 "mM" is a placeholder, not a printed number (Flags 6) | **not printed here** |
| `measured_volatiles.hexanal` | 782.0 ppb, 8 % | Table 2, UHT Processed: **782 +/- 58.6 ug/L**; 58.6/782 = 7.5 % | **matches** |
| `measured_volatiles."2-pentylfuran"` | 163.0 ppb, 9 % | Table 2, UHT Processed: **163 +/- 15.1 ug/L**; 15.1/163 = 9.3 % | **matches** |
| `measured_volatiles.nonanal` | 24.0 ppb, 4 % | Table 2, UHT Processed: **24.0 +/- 0.80 ug/L**; 0.80/24.0 = 3.3 % | **matches** |
| `quantification_class` | dynamic_headspace_gcmsms | "Aroma compounds were quantified using a standard addition technique and DHS coupled with GC/MS-QQQ ... multiple reaction monitoring (MRM) ... 5-point standard addition curves (r2 > 0.960)" | **matches** |
| `process_metadata.preheating` | 80C | "the samples were first homogenized at 17.2 MPa, **preheated to 80 C**" | **matches** |
| `protein_type` | pea_iso | "Pea protein isolate (PPI, 80 % protein, 8 % fat, 2 % sugar w/w) was obtained from a commercial source" | **matches** |
| `conditions.water_activity` | 0.98 | **not printed anywhere in this paper** | assumed by the bundle |
| `benchmark_id` year "2019" | — | the article is 2020; the DOI's `2019` is Elsevier's accepted-manuscript year | the bundle already records this as an opaque key; **confirmed wrong as a citation year, correctly retained as an identifier** |

The bundle's own note says the values were verified against the 2018 OSU MS thesis (OhioLINK
`osu1531495328317918`, Table 6 p. 35), whose Processed column reads 781.72 +/- 58.59, 163.16 +/-
15.06, 23.98 +/- 0.80. **The journal article rounds those to 782 +/- 58.6, 163 +/- 15.1 and 24.0
+/- 0.80.** So the thesis is the more precise source of the same three numbers, and the bundle's
loaded warning — that the Processed column must not be confused with the Control or Aged column —
is confirmed here from the article itself: Table 2's three columns for hexanal read **331 / 782 /
683 ug/L** and for 2-pentylfuran **59.4 / 163 / 197 ug/L**, which are materially different.

## 1. Why it matters

`tasks/roadmap_for_scientists.md` section 5d, Programme 7, opens by saying a pea or soy isolate
"carries 1 to 3 % lipid and its own volatiles (hexanal, 2-pentylfuran, 1-octen-3-ol,
methoxypyrazines) into every cook", and part (ii) asks that those volatiles be charged as declared
inputs "with their measured levels and bands". **This paper is the closest thing in the corpus to
that measurement done properly:**

1. **The material is an isolate, named as such**, with its composition printed: 80 % protein, 8 %
   fat, 2 % sugar w/w. Note that the fat is **8 %**, well above the roadmap's "1 to 3 % lipid" —
   worth carrying, because the carried-volatile budget scales with the lipid that comes with it
   (Flags 2).
2. **All four of the roadmap's named volatiles are quantified in the UNPROCESSED isolate
   beverage**, which no other paper in this cluster manages: hexanal **331 +/- 81.3 ug/L**,
   2-pentylfuran **59.4 +/- 1.93 ug/L**, 1-octen-3-ol **7.13 +/- 0.19 ug/L**, and the
   methoxypyrazine 2-isobutyl-3-methoxypyrazine **0.031 +/- 0.003 ug/L**. Zhang 2020b
   (`zhang2020b_extraction.md`) measures a seed-based pea *milk*, not an isolate; Fischer 2021
   (`fischer2021_extraction.md`) is the other roadmap-named source; `bi2020_extraction.md` is a
   *flour* and prints no 2-pentylfuran and no methoxypyrazine at all.
3. **It prints what heat does to each of them at a stated cook** (140 C, 6 s), and **what seven
   weeks at 5 C then do**, all against the same control, in the same units, from the same
   laboratory. That is three states of the same material — the axis a "carried input plus a loss
   or gain term" needs.
4. **It prints an odour threshold in water for all 21 compounds**, sourced (Leffingwell for most,
   Buttery 1988 for 2-acetyl-1-pyrroline, Sterckx 2011 retronasal for 4-hydroxybenzaldehyde), and
   marks with an asterisk which compounds are above threshold in each state.
5. **The Maillard side is measured in the same run**: methional, 2-acetyl-1-pyrroline,
   2,5-dimethylpyrazine, maltol, p-vinylguaiacol, sulfurol. Methional and 2-acetyl-1-pyrroline in
   a **6-second** heat treatment at 140 C are a hard, short-time test for any Maillard lane.
6. Its **2-acetyl-1-pyrroline in the unheated control (0.29 ug/L)** is explicitly interpreted by the
   authors as having been formed "during protein extraction from the pea or during storage" — i.e.
   a *Maillard* product carried in by the isolate, not only lipid-oxidation products. Programme 7
   part (ii) as written charges only the lipoxygenase volatiles; this paper says the carried load
   includes Maillard markers too.

What it does not give: no rate constant, no second temperature, no intermediate time inside the
6-second hold, no storage point between 0 and 7 weeks, no lipid or free-fatty-acid measurement, no
enzyme activity, no water activity, and no measurement of the dry isolate itself — everything is
measured in the 3 % beverage.

## 2. Methods as they matter to a model

- **The material, exactly as described.** "**Pea protein isolate (PPI, 80 % protein, 8 % fat, 2 %
  sugar w/w) was obtained from a commercial source and stored at -80 C.**" No supplier, no
  cultivar, no extraction route, no batch, no lot. It is a **protein isolate** — not a flour, not a
  concentrate, not a whole seed. The remaining 10 % w/w of the isolate is unaccounted for in the
  printed composition (ash, moisture, fibre); the paper does not say (Flags 2).
- **The beverage.** "0.03 % w/w carrageenan (stabilizer) was added to nano-filtered water and
  mixed using a bench-top high shear mixer on medium speed for 5 min ... **Pea protein isolate
  (3 % w/w)** was then added slowly while mixing and covered for foam reduction at **room
  temperature for approximately 40 min** until no clumps were observed and mixture was
  homogeneous." So: 3 g isolate + 0.03 g carrageenan per 100 g of water-based beverage,
  **2.4 % w/w protein and 0.24 % w/w fat by my arithmetic**, held stirred and open (covered, but
  not sealed) for about 40 min at room temperature before anything else happens. **That 40-minute
  hydration window is itself a reactive window** for whatever lipoxygenase survived isolate
  manufacture, and it is inside the "non-UHT control" (Flags 3).
- **The cook, exactly as described.** MicroThermics HTST/UHT **indirect tubular** processing system
  with homogenizer and aseptic filling station. "During the UHT treatment, the samples were first
  **homogenized at 17.2 MPa, preheated to 80 C**, then processed at a final temperature and
  pressure of **140 C and 0.55 MPa for 6 s** to ensure sterilization, **rapidly cooled to 10 C**,
  and filled in a sterilized laminar flow hood into pre-sterilized 500-mL bottles." No temperature
  -time profile is logged or printed; the come-up from 80 C to 140 C and the cool-down to 10 C are
  not quantified, so the effective thermal load exceeds 6 s at 140 C by an unstated amount
  (Flags 4).
- **The three states.** (i) **non-UHT processed (control)** — the beverage as mixed, no heat, taken
  straight to -80 C; (ii) **UHT processed** — the 140 C / 6 s product, "immediately stored at
  -80 C"; (iii) **UHT aged** — the same product held **5 C for 7 weeks**, then moved to -80 C
  before analysis. Note the control is **not** a time-zero of the UHT line: it never went through
  the homogenizer or the 80 C preheat (Flags 4).
- **Microbiology.** 3M Petrifilm aerobic count and coliform count plates on all three samples; "no
  detectable microbial growth". So the storage change is chemical, not microbial — a real point in
  the paper's favour for using the aged column as a storage series.
- **pH.** 7.1 in the control, "no significant pH changes ... on the UHT processed and aged samples
  relative to the control beverage". A rare case where the pH is stated to be stable across the
  whole design.
- **Identification: SAFE + GC/MS/O.** 200 g of beverage spiked with 2-methyl-3-heptanone
  (5.5 ug/mL in DCM) as internal standard, 70 g NaCl, extracted three times with DCM (2 x 50 mL,
  1 x 20 mL), 1 h on an orbital shaker at room temperature, centrifuged 30 min at 2700 g / 4 C;
  pooled organics dried over Na2SO4, filtered, then **solvent-assisted flavour evaporation** (Engel
  1999) at 40 C body/legs, liquid-nitrogen receiver, ~1.3 mPa vacuum, ~90 min; concentrated to
  500 uL on a 60 cm Vigreux. GC/MS/O on an Agilent 6890N + HP-5 (30 m x 0.25 mm x 0.25 um) + 5973
  MSD with a Gerstel ODP2 sniff port, 1:1 split; CIS 20 C (0.2 min) ramped 5 C/min to 250 C; oven
  40 C, 5 C/min to 250 C, hold 4 min; EI, 30-300 amu. **Three GC/O panelists, duplicate; a compound
  is selected only if detected by at least 50 % of panelists and above dilution 4 (flavour dilution
  FD > 8).** Identification by mass spectrum, linear retention index (n-alkanes C6-C26) **and
  injection of a pure standard**.
- **Quantification: DHS-GC/MS-QQQ with standard addition — and authentic standards throughout.**
  "Aroma compounds were quantified using a **standard addition** technique and DHS coupled with
  GC/MS-QQQ. The twenty-one selected volatile compounds from GC/MS/O had a wide range of
  concentrations and trap affinities; therefore, **different DHS trapping and drying methods were
  applied**, with methods ranging from 50 to 500 mL/min trapping and 100-1000 mL/min drying." CIS
  -50 C, 1 min equilibration, 0.1 min hold, 12 C/min to 240 C, 3 min hold; TDU 40 C (0.5 min),
  300 C/min to 250 C, 1 min. Columns HP-5MS and VF-WAXMS (60 m), "depending on their polarity and
  volatility". MRM; the transitions, collision energies and detector gains are printed in Table 1
  and were "optimized by injection of pure standards". "Compounds were quantified using **5-point
  standard addition curves (r2 > 0.960)** ... Concentrations are given as an average of
  triplicates."
  **Standard addition is the right method here**: it builds the calibration inside the sample's own
  protein matrix, so the matrix suppression a 2.4 % protein dispersion exerts on headspace release
  is corrected for compound by compound. This is a stronger basis than the external- or
  matrix-blank calibrations used elsewhere in this cluster.
- **2-Acetyl-1-pyrroline is different: a stable isotope dilution assay.** "One hundred uL of
  0.05 ug/L of [2H3]-2AP were added to 5 mL of sample in a 20-mL glass vial for analysis.
  **Concentration of 2AP was calculated assuming an equal response factor with the deuterated
  standard.**" So one compound out of 21 has an isotopic internal standard, with an *assumed*
  response factor of 1 rather than a measured one.
- **Response factors generally.** For the other 20, the standard-addition slope *is* the response
  factor in that matrix; no separate response factor is reported or needed. All 21 named standards
  were bought (Sigma-Aldrich, except nonanal from Alfa Aesar), and the deuterated 2AP was
  synthesised internally.
- **Basis of every concentration.** **ug L^-1 of beverage** (i.e. ppb w/v of the 3 % w/w pea protein
  isolate drink), n = 3, mean +/- standard deviation. There is **no per-isolate and no dry-matter
  basis anywhere in the paper**; converting to a per-isolate basis is arithmetic I do in section 3
  and mark as mine.
- **Sensory.** Ten panelists (4 M, 6 F) from the OSU Flavor Research and Education Center; term
  generation and training over ten 45-min sessions; seven orthonasal attributes with physical
  references (beany = crushed boiled chickpeas; cooked green beans = green beans boiled 4 min; saw
  dust; pasta; potato; cardboard; oxidized = oxidized walnut); 10-point line scale anchored 0 "not
  present" to 10 "strong"; duplicate on two days; Compusense Cloud.
- **Recombination.** Built on the **non-UHT beverage as the base**, not on water, "in an effort to
  account for physicochemical interactions between the protein and aroma compounds". Two 5 %
  ethanolic stock solutions (one for the UHT state, one for the aged state) containing every
  compound that was both above its water threshold and significantly changed; 10 mL of control
  beverage + 100 uL of stock in a 60-mL amber bottle; equilibrated 3 h at room temperature.
- **Statistics.** One-way ANOVA per compound with Tukey HSD (alpha = 0.05) on the quantitative data;
  three-way ANOVA (sample, panelist, replicate and interactions) on the sensory data; SPSS.

## 3. Tables re-typed

### Table 1. "Optimized Multiple Reaction Monitoring (MRM) parameters for the selected aroma compounds"

Footnote as printed: "* Internal standard."

| identified compound | MW (g/mol) | precursor ion | product ion, quantifier | product ion, qualifier | CE (eV) | gain |
|---|---|---|---|---|---|---|
| 1-Pentanol | 88.2 | 70 | 55 | 42 | 5 | 1 |
| Hexanal | 100.2 | 82 | 41 | 67 | 20 | 1 |
| Isovaleric acid | 102.1 | 60 | 45 | 42 | 8 | 75 |
| 4-Heptanone* | 114.2 | 114 | 99 | 71 | 8 | 2 |
| 2-Heptanone | 114.2 | 114 | 99 | 71 | 5 | 0.5 |
| Heptanal | 114.2 | 96 | 81 | 79 | 10 | 0.5 |
| Methional | 104.2 | 104 | 48 | 61 | 8 | 50 |
| 2,5-Dimethylpyrazine | 108.1 | 108 | 42 | 81 | 8 | 35 |
| 2-Acetyl-1-pyrroline | 111.1 | 111 | 83 | 41 | 5 | 100 |
| [2H3]-2-Acetyl-1-pyrroline | 114.1 | 114 | 86 | 41 | 5 | 100 |
| 1-Octen-3-ol | 128.2 | 72 | 43 | 57 | 5 | 10 |
| 2-Pentylfuran | 138.2 | 138 | 81 | 94 | 5 | 10 |
| (E)-2-Octenal | 126.2 | 83 | 49 | 55 | 25 | 35 |
| Nonanal | 142.2 | 98 | 41 | 56 | 5 | 25 |
| Maltol | 126.1 | 126 | 97 | 71 | 15 | 75 |
| Octanoic acid | 144.2 | 73 | 55 | 45 | 10 | 75 |
| 2-Isobutyl-3-methoxypyrazine | 166.2 | 124 | 81 | 95 | 15 | 35 |
| (E,E)-2,4-Nonadienal | 138.2 | 81 | 53 | 51 | 25 | 5 |
| Sulfurol | 143.2 | 112 | 85 | 45 | 5 | 100 |
| (E,E)-2,4-Decadienal | 152.2 | 81 | 53 | 51 | 20 | 5 |
| p-Vinylguaiacol | 150.2 | 150 | 77 | 107 | 20 | 10 |
| 4-Hydroxybenzaldehyde | 122.1 | 121 | 65 | 93 | 5 | 10 |
| gamma-Nonalactone | 156.2 | 85 | 57 | 41 | 5 | 10 |

**4-Heptanone is the DHS internal standard**, distinct from the 2-methyl-3-heptanone used for the
SAFE identification extract.

### Table 2. "Average concentrations (n = 3) of aroma-active compounds in non-UHT processed, UHT processed and UHT aged pea protein beverages"

Footnotes exactly as printed: 1 "In each row, values with different letter are significantly
different according to Tukey's HSD (alpha = 0.05)"; 2 "Detected at Flavor Dilution >= 8"; 3 "Odor
threshold value in water; retrieved from Leffingwell & Associates except where noted"; 4 "Buttery
et al. (1988)"; 5 "Retronasal odor threshold, Sterckx, Missiaen, Saison, and Delvaux (2011)";
"* Compounds reported above odor threshold in water."

**This is the load-bearing table of the paper. Rows are in the printed order (HP-5 LRI ascending).**

| HP-5 LRI | descriptor | compound | odour threshold (ug/L, water) | non-UHT processed (ug/L) | UHT processed (ug/L) | UHT aged (ug/L) |
|---|---|---|---|---|---|---|
| 768 | Floral | 1-Pentanol | 4000 | 38.6 +/- 0.73 a | 115 +/- 3.75 b | 113 +/- 4.65 b |
| 802 | Cut grass | **Hexanal** | **4.5** | **331 +/- 81.3 a\*** | **782 +/- 58.6 b\*** | **683 +/- 58.1 b\*** |
| 827 | Animal, yeasty | Isovaleric acid | 120 | 1010 +/- 215 * | 1220 +/- 241 * | 1380 +/- 293 * |
| 890 | Brothy | 2-Heptanone | 140 | 26.8 +/- 1.40 a | 32.5 +/- 1.36 b | 37.7 +/- 1.31 c |
| 902 | Brothy, roasted | Heptanal | 3 | 6.95 +/- 0.56 a* | 21.3 +/- 2.39 c* | 12.8 +/- 2.55 b* |
| 910 | Meaty, brothy | 2,5-Dimethylpyrazine | 800 | 2.46 +/- 0.12 b | 2.29 +/- 0.13 b | 1.74 +/- 0.11 a |
| 911 | Boiled potato | **Methional** | 0.2 | **0.55 +/- 0.04 a\*** | **3.10 +/- 0.25 b\*** | **3.21 +/- 0.30 b\*** |
| 935 | Corn chip, pasta | **2-Acetyl-1-pyrroline** | 0.14 | **0.29 +/- 0.01 a\*** | **0.41 +/- 0.05 a\*** | **0.78 +/- 0.12 b\*** |
| 977 | Mushroom | **1-Octen-3-ol** | 1 | **7.13 +/- 0.19 a\*** | **14.0 +/- 1.32 b\*** | **13.0 +/- 0.77 b\*** |
| 1003 | Solvent-like, floral | **2-Pentylfuran** | 6 | **59.4 +/- 1.93 a\*** | **163 +/- 15.1 b\*** | **197 +/- 6.01 c\*** |
| 1060 | Dusty, musty, moldy | (E)-2-Octenal | 3 | 0.31 +/- 0.09 a | 13.5 +/- 0.49 b* | 13.1 +/- 0.89 b* |
| 1103 | Plastic, waxy | **Nonanal** | 1 | **8.24 +/- 0.44 a\*** | **24.0 +/- 0.80 c\*** | **22.4 +/- 0.28 b\*** |
| 1141 | Burnt sugar | Maltol | 35,000 | 1440 +/- 192 | 1820 +/- 181 | 1980 +/- 429 |
| 1160 | Plastic, green | Octanoic acid | 3000 | 1010 +/- 109 a | 1000 +/- 249 a | 2840 +/- 401 b |
| 1181 | Green pepper | **2-Isobutyl-3-methoxypyrazine** | **0.002** | **0.031 +/- 0.003 a\*** | **0.072 +/- 0.001 b\*** | **0.115 +/- 0.018 c\*** |
| 1215 | Brothy, meaty, plastic | (E,E)-2,4-Nonadienal | 0.09 | 0.71 +/- 0.05 a* | 4.84 +/- 0.35 c* | 1.30 +/- 0.14 b* |
| 1286 | Brothy, nutty | Sulfurol | 10,800 | 113 +/- 17.9 a | 226 +/- 70.3 ab | 278 +/- 61.3 b |
| 1327 | Nutty | (E,E)-2,4-Decadienal | 0.07 | 0.06 +/- 0.005 a | 46.9 +/- 2.53 c* | 16.6 +/- 1.70 b* |
| 1336 | Roasted peanut | p-Vinylguaiacol | 3 | 33.8 +/- 1.40 a* | 200 +/- 12.13 c* | 133 +/- 3.99 b* |
| 1369 | Sweet milk, plastic | 4-Hydroxybenzaldehyde | 10,000 (footnote 5, retronasal) | 2.30 +/- 0.23 c | 1.05 +/- 0.07 b | 0.56 +/- 0.03 a |
| 1381 | Sweet, peachy | gamma-Nonalactone | 65 | 4.12 +/- 0.36 a | 8.56 +/- 0.06 b | 9.07 +/- 0.49 b |

Footnote 4 ("Buttery et al. 1988") carries no visible superscript marker in the text layer; from the
running text — "2AP ... exhibiting a very low odor threshold of 0.1 ug/L in water (Buttery,
Turnbaugh, & Ling, 1988)" — it attaches to the 2-acetyl-1-pyrroline row, whose printed threshold is
**0.14**, not 0.1 (Flags 5).

### Table 3. "Summary of Variance Analysis (ANOVA) results on sensory evaluation data for the pea protein beverage samples"

Footnote a: "= Non-UHT processed, UHT processed, and UHT aged."

| attribute | p, sample | p, panelist | p, rep | p, panelist*sample |
|---|---|---|---|---|
| Beany | 0.024 | 0.009 | 0.900 | 0.283 |
| Cooked green bean | 0.830 | 0.426 | 0.769 | 0.069 |
| Saw dust | < 0.001 | 0.012 | 0.579 | 0.114 |
| Pasta | 0.800 | 0.320 | 0.792 | 0.434 |
| Potato | 0.064 | 0.616 | 0.484 | 0.970 |
| Cardboard | 0.002 | < 0.001 | 0.623 | 0.210 |
| Oxidized/painty | 0.001 | 0.002 | 0.893 | 0.487 |

### Table 4. "Summary of Variance Analysis (ANOVA) results on sensory evaluation data for the pea protein beverage recombination model samples"

Footnote a: "= Non-UHT processed, UHT processed recombination model (non-UHT processed + UHT
processed compounds), and UHT aged recombination model (non-UHT processed + UHT aged compounds)."

| attribute | p, sample | p, panelist | p, rep | p, panelist*sample |
|---|---|---|---|---|
| Beany | 0.005 | 0.002 | 0.876 | 0.040 |
| Cooked green bean | 0.782 | 0.281 | 0.740 | 0.157 |
| Saw dust | 0.001 | 0.041 | 0.602 | 0.690 |
| Pasta | 0.783 | 0.312 | 0.086 | 0.405 |
| Potato | 0.063 | 0.568 | 0.470 | 0.363 |
| Cardboard | 0.005 | 0.003 | 0.641 | 0.657 |
| Oxidized/painty | 0.006 | 0.016 | 0.904 | 0.832 |

### Numbers printed only in the running text

| quantity | value | where |
|---|---|---|
| compounds identified | 21, across alcohols, aldehydes, ketones, pyrroles, carboxylic acids, pyrazines, furans, lactone and phenols | section 3.1 |
| first reports in pea protein | 2-pentylfuran, 2-heptanone, (E,E)-2,4-nonadienal, maltol, octanoic acid; also 4-hydroxybenzaldehyde; **2-acetyl-1-pyrroline "reported here for the first time in a pea protein isolate"** | section 3.1 |
| compounds significantly changed by processing and/or storage | 19 of 21; only isovaleric acid and maltol did not change | section 3.1 |
| compounds significantly increased by UHT | 15: 1-pentanol, hexanal, heptanal, methional, 1-octen-3-ol, (E)-2-octenal, nonanal, IBMP, sulfurol, (E,E)-2,4-decadienal, p-vinylguaiacol, gamma-nonalactone, 2-pentylfuran, 2-heptanone, (E,E)-2,4-nonadienal — **11 of the 15 are fatty-acid oxidation products** | section 3.1 |
| largest increase on UHT | "**up to 800-fold** for (E,E)-2,4-decadienal" | section 3.1 |
| heptanal and nonanal on UHT | "increased about **3-fold**" | section 3.1 |
| methional on UHT | text says "**~3-fold**" in section 3.1 and "**6-fold**" in section 3.2; the table gives 3.10/0.55 = 5.6 (mine) | sections 3.1, 3.2 (Flags 5) |
| compounds significantly increased in storage | 2-heptanone, 2-pentylfuran, octanoic acid, "approximately 1.2 to 3-fold"; plus IBMP and 2-acetyl-1-pyrroline | section 3.1 |
| compounds that fell in storage | heptanal, nonanal, (E,E)-2,4-nonadienal, (E,E)-2,4-decadienal, p-vinylguaiacol, 4-hydroxybenzaldehyde | section 3.1 |
| main pea fatty acid | linoleic acid, named as the precursor of 1-pentanol, hexanal, 2-heptanone, 1-octen-3-ol, 2-pentylfuran, (E)-2-octenal, octanoic acid, (E,E)-2,4-nonadienal, (E,E)-2,4-decadienal and gamma-nonalactone | section 3.1 |
| second pea fatty acid | oleic acid, named as the precursor of heptanal and nonanal | section 3.1 |
| beany percept | 1-octen-3-one with hexanal at a ratio **1:100** gives an intense beany note (from Bott & Chambers 2006, not measured here) | Introduction |
| beany attribute intensity | "significantly decreased from **3.3** to below **2** in the UHT processed sample" | section 3.2 |
| potato attribute intensity | "significantly decreased from **2.7** to below **2**" | section 3.2 |
| cardboard / oxidized / saw dust intensities | "about **1** in the non-UHT sample and between **2.5 and 3.2** in the UHT processed and UHT aged samples" | section 3.2 |
| non-UHT profile | beany, potato, pasta, cooked green bean, intensities **2.1 to 3.3** | section 3.2 |
| recombination fit | no significant difference (t-test, alpha = 0.05) for **6 of 7** attributes; only **pasta** was rated higher in the models ("data not shown") | section 3.2 |

**Figure-only:** Figure 1 is a proposed reaction scheme for IBMP via 2-isobutyl-3-hydroxypyrazine
(a mechanism, no numbers). Figure 2a (three beverages) and 2b (two recombination models) carry the
seven attribute intensities; only the few values quoted in the text above are printed, the rest are
**figure_only**.

### Arithmetic on the printed numbers (all mine)

1. **Fold changes on UHT (140 C, 6 s), UHT / non-UHT.** (E,E)-2,4-decadienal **782x** (the paper's
   "up to 800-fold"); (E)-2-octenal **43.5x**; p-vinylguaiacol **5.92x**; (E,E)-2,4-nonadienal
   **6.82x**; methional **5.64x**; heptanal **3.06x**; nonanal **2.91x**; 2-pentylfuran **2.74x**;
   **hexanal 2.36x**; IBMP 2.32x; sulfurol 2.00x; 1-octen-3-ol 1.96x; gamma-nonalactone 2.08x;
   1-pentanol 2.98x; 2-heptanone 1.21x; maltol 1.26x; isovaleric acid 1.21x; 2-acetyl-1-pyrroline
   1.41x (not significant); octanoic acid 0.99x; 2,5-dimethylpyrazine 0.93x; 4-hydroxybenzaldehyde
   0.46x.
2. **Fold changes over 7 weeks at 5 C, aged / UHT.** octanoic acid **2.84x**; 2-acetyl-1-pyrroline
   **1.90x**; IBMP **1.60x**; 2-pentylfuran 1.21x; 2-heptanone 1.16x; sulfurol 1.23x; maltol 1.09x;
   isovaleric acid 1.13x; gamma-nonalactone 1.06x; methional 1.04x; hexanal 0.87x; nonanal 0.93x;
   1-octen-3-ol 0.93x; (E)-2-octenal 0.97x; 1-pentanol 0.98x; heptanal **0.60x**; p-vinylguaiacol
   **0.67x**; 4-hydroxybenzaldehyde 0.53x; (E,E)-2,4-nonadienal **0.27x**; (E,E)-2,4-decadienal
   **0.35x**; 2,5-dimethylpyrazine 0.76x.
3. **Average storage slopes over the 7 weeks (two-point, zero-order reading — NOT rates).**
   2-acetyl-1-pyrroline +0.053 ug/L per week; IBMP +0.0061 ug/L per week; 2-pentylfuran
   +4.86 ug/L per week; octanoic acid +263 ug/L per week; hexanal **-14.1 ug/L per week**;
   (E,E)-2,4-decadienal **-4.33 ug/L per week**. Two points is not a shape; these are averages over
   the interval and cannot distinguish a fast early change from a linear one.
4. **Per-isolate basis (mine, and the assumption is load-bearing).** The beverage is 3 % w/w
   isolate. Taking the beverage density as 1.00 kg/L and assuming **every molecule measured came in
   with the isolate and none was lost**, a level of C ug/L in the beverage is C / 0.03 ug per kg of
   isolate. On the **non-UHT control** column that gives, per kilogram of pea protein isolate:
   hexanal **11,000 ug/kg**, 2-pentylfuran **1980 ug/kg**, nonanal **275 ug/kg**, 1-octen-3-ol
   **238 ug/kg**, IBMP **1.03 ug/kg**, 1-pentanol 1290 ug/kg, heptanal 232 ug/kg, methional
   18.3 ug/kg, 2-acetyl-1-pyrroline 9.7 ug/kg, isovaleric acid 33,700 ug/kg, maltol 48,000 ug/kg,
   octanoic acid 33,700 ug/kg, sulfurol 3770 ug/kg, p-vinylguaiacol 1130 ug/kg,
   (E,E)-2,4-nonadienal 23.7 ug/kg, (E)-2-octenal 10.3 ug/kg, (E,E)-2,4-decadienal 2.0 ug/kg,
   gamma-nonalactone 137 ug/kg, 2-heptanone 893 ug/kg, 2,5-dimethylpyrazine 82 ug/kg,
   4-hydroxybenzaldehyde 76.7 ug/kg. **These are `derived_assumption`, not measurements** — the
   40-minute hydration at room temperature sits between the dry isolate and the control sample and
   may itself generate some of them (Flags 3).
5. **The cross-study check the levels table wants.** On a per-protein basis: this control beverage
   is 2.4 % w/w protein (3 % isolate x 80 % protein), so hexanal is 331/0.024 = **13,800 ug per kg
   of protein**. Zhang 2020b's raw pea milk is 164.18 ug/L at 2.00 % protein =
   **8210 ug per kg of protein**. **The two laboratories agree within a factor 1.7** on the hexanal
   a pea protein preparation carries per unit protein, from completely different materials
   (commercial isolate vs laboratory seed milk) and different methods (DHS standard addition vs
   matrix-matched external curves). That is the strongest single argument in this cluster that a
   carried-hexanal input of order 10^4 ug per kg of protein is a real, transferable number.
   Bi 2020's raw pea flour gives 1260 ug/kg of flour, which at a nominal 22 % protein is
   ~5700 ug/kg of protein (mine, the flour's protein content is **not** printed in Bi 2020, so 22 %
   is an assumption from typical pea composition and this third comparison is the weakest of the
   three).
6. **Odour-activity values at the control state (mine, concentration / printed water threshold).**
   IBMP **15.5**; hexanal **73.6**; 1-octen-3-ol **7.1**; 2-pentylfuran **9.9**; nonanal **8.2**;
   (E,E)-2,4-nonadienal **7.9**; methional 2.75; 2-acetyl-1-pyrroline 2.07; heptanal 2.32;
   isovaleric acid 8.4; p-vinylguaiacol 11.3; (E)-2-octenal 0.10; (E,E)-2,4-decadienal 0.86;
   2-heptanone 0.19; 1-pentanol 0.010; maltol 0.041; octanoic acid 0.34; sulfurol 0.010;
   gamma-nonalactone 0.063; 2,5-dimethylpyrazine 0.003; 4-hydroxybenzaldehyde 0.0002. **Ten of the
   21 are already above threshold before any heat is applied**, which is the quantitative form of
   the roadmap's claim that the isolate carries its own aroma into every cook. The paper's own
   asterisks mark the same ten.
7. **At the UHT state the OAV ranking changes.** hexanal 174, (E,E)-2,4-decadienal **670**,
   (E,E)-2,4-nonadienal 53.8, IBMP 36, 2-pentylfuran 27.2, nonanal 24.0, methional 15.5,
   1-octen-3-ol 14.0, p-vinylguaiacol 66.7, heptanal 7.1, (E)-2-octenal 4.5, 2-acetyl-1-pyrroline
   2.9, isovaleric acid 10.2 (mine). The compound with the highest odour activity after the cook is
   **(E,E)-2,4-decadienal**, which is not in the repository's registry at all (Flags 9).

## 4. Numbers the repository can use

All rows share: commercial pea protein isolate (80 % protein, 8 % fat, 2 % sugar w/w) at 3 % w/w
with 0.03 % w/w carrageenan in nano-filtered water, **pH 7.1**, n = 3, DHS-GC/MS-QQQ in MRM with
five-point standard addition (r2 > 0.960) except 2-acetyl-1-pyrroline by SIDA; basis **ug L^-1 of
beverage**. Three states: **control** = unheated; **UHT** = homogenized 17.2 MPa, preheated 80 C,
140 C / 0.55 MPa for 6 s, cooled to 10 C; **aged** = the UHT product held 7 weeks at 5 C.

| compound | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| hexanal | 331 +/- 81.3 | ug/L of beverage | pea protein isolate 3 % w/w, **unheated**, pH 7.1 | Table 2 p. 4 | measured_level |
| hexanal | 782 +/- 58.6 | ug/L of beverage | same, UHT 140 C / 6 s | Table 2 | measured_level (**the benchmark's value**) |
| hexanal | 683 +/- 58.1 | ug/L of beverage | same, +7 weeks at 5 C | Table 2 | measured_level |
| 2-pentylfuran | 59.4 +/- 1.93 / 163 +/- 15.1 / 197 +/- 6.01 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level (the UHT value is the benchmark's) |
| 1-octen-3-ol | 7.13 +/- 0.19 / 14.0 +/- 1.32 / 13.0 +/- 0.77 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| nonanal | 8.24 +/- 0.44 / 24.0 +/- 0.80 / 22.4 +/- 0.28 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level (the UHT value is the benchmark's) |
| 2-isobutyl-3-methoxypyrazine (IBMP) | 0.031 +/- 0.003 / 0.072 +/- 0.001 / 0.115 +/- 0.018 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level — **the only methoxypyrazine concentration in this cluster** |
| 1-pentanol | 38.6 +/- 0.73 / 115 +/- 3.75 / 113 +/- 4.65 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| isovaleric acid (3-methylbutanoic acid) | 1010 +/- 215 / 1220 +/- 241 / 1380 +/- 293 | ug/L of beverage | control / UHT / aged; **not significantly changed** | Table 2 | measured_level |
| 2-heptanone | 26.8 +/- 1.40 / 32.5 +/- 1.36 / 37.7 +/- 1.31 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| heptanal | 6.95 +/- 0.56 / 21.3 +/- 2.39 / 12.8 +/- 2.55 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| 2,5-dimethylpyrazine | 2.46 +/- 0.12 / 2.29 +/- 0.13 / 1.74 +/- 0.11 | ug/L of beverage | control / UHT / aged — **falls**, never above threshold | Table 2 | measured_level |
| methional | 0.55 +/- 0.04 / 3.10 +/- 0.25 / 3.21 +/- 0.30 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level — a Strecker marker made in **6 s at 140 C** |
| 2-acetyl-1-pyrroline | 0.29 +/- 0.01 / 0.41 +/- 0.05 / 0.78 +/- 0.12 | ug/L of beverage | control / UHT / aged; UHT not significantly different from control | Table 2, by SIDA | measured_level (**response factor assumed equal to the d3 standard**) |
| (E)-2-octenal | 0.31 +/- 0.09 / 13.5 +/- 0.49 / 13.1 +/- 0.89 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| maltol | 1440 +/- 192 / 1820 +/- 181 / 1980 +/- 429 | ug/L of beverage | control / UHT / aged; **not significantly changed** | Table 2 | measured_level |
| octanoic acid | 1010 +/- 109 / 1000 +/- 249 / 2840 +/- 401 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| (E,E)-2,4-nonadienal | 0.71 +/- 0.05 / 4.84 +/- 0.35 / 1.30 +/- 0.14 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| sulfurol | 113 +/- 17.9 / 226 +/- 70.3 / 278 +/- 61.3 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level (a thiamine degradation marker) |
| (E,E)-2,4-decadienal | 0.06 +/- 0.005 / 46.9 +/- 2.53 / 16.6 +/- 1.70 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level — **the largest change in the paper, 782x** |
| p-vinylguaiacol (4-vinylguaiacol) | 33.8 +/- 1.40 / 200 +/- 12.13 / 133 +/- 3.99 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| 4-hydroxybenzaldehyde | 2.30 +/- 0.23 / 1.05 +/- 0.07 / 0.56 +/- 0.03 | ug/L of beverage | control / UHT / aged — **falls monotonically** | Table 2 | measured_level |
| gamma-nonalactone | 4.12 +/- 0.36 / 8.56 +/- 0.06 / 9.07 +/- 0.49 | ug/L of beverage | control / UHT / aged | Table 2 | measured_level |
| odour thresholds in water, all 21 | see the threshold column above (hexanal 4.5; 1-octen-3-ol 1; 2-pentylfuran 6; nonanal 1; IBMP 0.002; methional 0.2; 2-acetyl-1-pyrroline 0.14; heptanal 3; (E)-2-octenal 3; (E,E)-2,4-nonadienal 0.09; (E,E)-2,4-decadienal 0.07; p-vinylguaiacol 3; 1-pentanol 4000; isovaleric acid 120; 2-heptanone 140; 2,5-dimethylpyrazine 800; maltol 35,000; octanoic acid 3000; sulfurol 10,800; gamma-nonalactone 65; 4-hydroxybenzaldehyde 10,000 retronasal) | ug/L in water | — | Table 2, from Leffingwell & Associates except footnotes 4 and 5 | threshold |
| per-isolate levels, control column | hexanal 11,000; 2-pentylfuran 1980; nonanal 275; 1-octen-3-ol 238; IBMP 1.03 (full list in section 3) | ug per kg of isolate | 3 % w/w beverage, density assumed 1.00 kg/L, all volatile assumed carried in | derived from Table 2 (mine) | **derived_assumption** |
| per-protein hexanal, control | 13,800 | ug per kg of protein | 2.4 % w/w protein | derived (mine) | derived_assumption |
| UHT fold changes, 21 compounds | see section 3 arithmetic 1 | dimensionless | 140 C / 6 s against unheated | derived from Table 2 (mine) | within_study_ratio |
| storage fold changes, 21 compounds | see section 3 arithmetic 2 | dimensionless | 7 weeks at 5 C against the UHT product | derived from Table 2 (mine) | within_study_ratio |
| average storage slopes | e.g. 2-acetyl-1-pyrroline +0.053, IBMP +0.0061, 2-pentylfuran +4.86, hexanal -14.1 | ug/L per week | 5 C, two points 7 weeks apart | derived (mine) | **derived_assumption** — two points, no shape |
| odour-activity values, control and UHT | see section 3 arithmetic 6 and 7 | dimensionless | beverage concentration over water threshold | derived (mine) | within_study_ratio |
| MRM transitions, collision energies, gains | 22 compounds + 1 deuterated standard | m/z, eV, dimensionless | — | Table 1 p. 3 | level_only (instrument settings) |
| sensory ANOVA p-values, 7 attributes x 3 samples and x 3 models | see Tables 3 and 4 above | p | 10 panelists, duplicate | Tables 3, 4 p. 4 | level_only |
| attribute intensities | beany 3.3 -> below 2; potato 2.7 -> below 2; cardboard/oxidized/saw dust ~1 -> 2.5-3.2 | 10-point scale | control vs UHT and aged | text, section 3.2 | level_only (the rest of Fig. 2 is **figure_only**) |
| IBMP formation scheme | via 2-isobutyl-3-hydroxypyrazine, non-enzymatic methylation, trigonelline proposed as the methyl donor | — | proposed, not demonstrated | Fig. 1 | **figure_only** (a mechanism, no numbers) |

### What this can and cannot be used for

**Can:** be the anchor of Programme 7's levels table. It is the only source in this cluster with
the material the roadmap names (a pea protein isolate), all four volatiles the roadmap names, real
concentrations with standard deviations from triplicates, a stated pH, and a stated process. Its
control column is the "carried in" charge; its UHT column is the same charge after a stated cook;
its aged column is the same after a stated storage. It can also supply the odour thresholds that
`data/species/off_flavour_targets.yml` currently carries uncited for four of its six compounds.

**Cannot:** supply a rate. Six seconds is one point, seven weeks is one point, and the come-up is
unlogged. Cannot be put on a dry-isolate basis without the assumption in section 3 item 4. Cannot
be transferred to a different isolate: the supplier, extraction route and batch are all withheld,
and the fat content (8 %) is at the top of, or above, the roadmap's assumed 1-3 % band, so the
carried lipid-oxidation load of this isolate is probably at the high end of the population.

## 5. Flags

1. **The benchmark rests on the middle column of a three-column table, and the columns differ by up
   to a factor 800.** Confirmed here from the article: the Processed column is 782 / 163 / 24.0
   ug/L, the Control column 331 / 59.4 / 8.24, the Aged column 683 / 197 / 22.4. The bundle's
   `process_metadata.state = "heated_matrix"` is the only thing selecting Processed. That guard is
   correct and load-bearing and this dossier confirms it independently.
2. **The isolate is 8 % fat, not the 1-3 % the roadmap assumes.** "PPI, 80 % protein, 8 % fat, 2 %
   sugar w/w". Ten per cent w/w of the isolate is unaccounted for (moisture, ash, fibre). If the
   carried-volatile load scales with the lipid that comes with the isolate, this isolate is a
   high-lipid example and its levels should not be treated as the population median. The paper
   names **no supplier, no cultivar, no extraction route (isoelectric vs salt vs dry
   fractionation), no batch and no manufacture date**, so it cannot be located on any such
   population.
3. **The "non-UHT control" is not a zero-time of the isolate.** Between the dry isolate at -80 C
   and the control sample there is a 5-min carrageenan pre-mix, a slow addition under high shear,
   and **about 40 minutes covered at room temperature**. If any lipoxygenase activity survived
   isolate manufacture, part of the control's hexanal and 1-octen-3-ol was made in that window
   rather than carried in by the powder. The paper does not assay any enzyme and does not measure
   the dry isolate. **The control column is therefore an upper bound on the truly carried level,
   and my per-isolate conversion in section 3 inherits that bound.**
4. **The thermal history is not logged and the control did not experience it.** The UHT sample was
   homogenized at 17.2 MPa, preheated to 80 C, taken to 140 C for 6 s and cooled to 10 C; no
   temperature-time trace is printed, so the come-up and cool-down contribute an unknown extra
   thermal load. Separately, the control **never went through the homogenizer or the 80 C
   preheat**, so the control-to-UHT difference conflates the 140 C hold with homogenization and
   with the preheat. A model fitting this pair to a 6 s isothermal step will over-attribute.
5. **Three internal inconsistencies in the printed text.** (i) The 2-acetyl-1-pyrroline threshold
   is **0.14 ug/L** in Table 2 but **0.1 ug/L** in the running text of section 3.1, both credited to
   Buttery 1988. (ii) Methional's rise on UHT is called "~3-fold" in section 3.1 and "6-fold" in
   section 3.2; the table gives **5.6x** (mine), so the section 3.2 statement is the closer one and
   the section 3.1 statement is wrong. (iii) The abstract says "twenty-one aroma compounds were
   identified"; Table 1 lists 22 compounds plus the deuterated standard, because it includes
   **4-heptanone, which is the internal standard, not an identified odorant**. Table 2 correctly
   lists 21.
6. **The benchmark's `concentration_mM: 1000.0` for "Pea Protein Isolate" is not a printed
   number.** The paper states a mass fraction (3 % w/w), never a molarity, and a protein isolate
   has no single molar mass. The 1000 mM is a placeholder in the bundle. Recorded so that no one
   later reads it back as a measurement.
7. **The benchmark's `water_activity: 0.98` is not printed here either.** No water activity, no
   moisture content and no dry-matter figure appears anywhere in this paper. 0.98 is an assumption
   about a 3 % w/w aqueous dispersion, and it is a reasonable one, but it is not this paper's
   number.
8. **What this paper does NOT contain, and what to request.** No lipid or fatty-acid profile of the
   isolate or the beverage (linoleic and oleic acid are named as precursors from the literature,
   never measured); no lipoxygenase or any enzyme activity; no free amino acid or sugar profile
   beyond "2 % sugar"; no measurement of the dry isolate; no second temperature; no time point
   inside the 6 s hold or between 0 and 7 weeks; no water activity; no replicate isolate or
   supplier; no reported detection or quantification limit for any of the 21 compounds; no
   recovery figure for the DHS traps. **To request:** (a) **Supplementary Table S1** (the 21
   standard-addition curves), which is not on disk and which would give the slopes and hence the
   effective matrix suppression for each compound in a 2.4 % pea protein dispersion — that is
   exactly the "matrix binding" quantity Programme 7 wants to apply, and it is one table away;
   (b) the isolate's supplier and specification sheet; (c) intermediate storage points; (d) the
   underlying intensity values behind Figure 2a and 2b; (e) Murat et al. 2013, *Food Res. Int.*
   53:31-41, this paper's own key reference for the flour-to-protein-extract volatile series,
   which is **not on disk** and which `bi2020_extraction.md` flags for the same reason.
9. **Registry gaps against `data/keys/compounds.yml`.** Present and directly keyable: `hexanal`,
   `nonanal`, `heptanal`, `1_octen_3_ol`, `2_pentylfuran`, `e_2_octenal`, `methional`,
   `2_acetyl_1_pyrroline`, `2_5_dimethylpyrazine`, `4_vinylguaiacol` (the paper's
   p-vinylguaiacol), and `3_isobutyl_2_methoxypyrazine` — the last one **is** this paper's
   2-isobutyl-3-methoxypyrazine, the same molecule under the other locant ordering (registry
   InChIKey UXFSPRAGHGMRSQ-UHFFFAOYSA-N), so the registry can take the IBMP numbers directly. The
   class id `methoxypyrazines` also applies. **Absent from the registry:** 1-pentanol, isovaleric
   acid (3-methylbutanoic acid), 2-heptanone, maltol, octanoic acid, **(E,E)-2,4-nonadienal**,
   sulfurol, **(E,E)-2,4-decadienal**, 4-hydroxybenzaldehyde and gamma-nonalactone. The two
   dienals matter most: **(E,E)-2,4-decadienal is the highest-odour-activity compound in the cooked
   beverage** (OAV 670, mine) and the repository cannot name it; the registry's nearest entry is
   `e_e_2_4_heptadienal`, a different molecule.
10. **Registry gaps against `data/species/off_flavour_targets.yml`.** That file holds six
    compounds. This paper prints a level for four of them — hexanal, nonanal, 1-octen-3-ol,
    2-pentylfuran — in **all three states**, and is silent on 1-hexanol and furfural. Its
    thresholds agree with the file exactly for hexanal (4.5), 1-octen-3-ol (1.0), 2-pentylfuran
    (6.0) and nonanal (1.0), from a different compilation (Leffingwell) than the file's own
    Belitz/Czerny attribution — an independent corroboration of all four. **The file has no entry
    for any methoxypyrazine**, even though the roadmap names methoxypyrazines as a carried volatile
    and this paper measures one at an odour activity of 15.5 before any heat and 36 after
    (mine). Adding IBMP to `off_flavour_targets.yml`, with threshold 0.002 ug/L, is the single
    clearest registry action this dossier supports. The same file also has no entry for
    (E,E)-2,4-decadienal, (E,E)-2,4-nonadienal or (E)-2-octenal, three lipid-oxidation aldehydes
    this paper shows rising by 43x to 782x across a six-second cook.
