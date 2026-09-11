# Guo 2020 — EXTRACTION (soy protein isolate + 0-40 % wheat gluten at 50/60/70/80 % material moisture, high-moisture twin-screw extrusion to 150 C, then headspace GC-MS retention of a 16-compound spice flavouring against the unextruded raw material; plus LF-NMR water distribution and FT-IR secondary structure)

### NOT A BINDING STUDY: this paper contains no binding constant, no partition coefficient and no equilibrium of any kind — it is a PROCESS-RETENTION study, and what it delivers to the repository is the corpus's only measurement of how much of a dosed flavour survives high-moisture extrusion of a soy matrix, with the aldehyde class retaining only 7.8-21.2 % and the whole flavour load falling 2.3x from 50 % to 80 % material moisture.

**Source on disk:** `data/articles/guo2020.pdf` (11 pp., Food Hydrocolloids 105 (2020) 105752).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/guo2020.txt`); **Tables 1-6 all
came through clean** and are re-typed in full below. The extraction renders the plus-or-minus sign
as a stray character throughout ("59.56 � 1.15") and the degree sign the same way; both are restored
here. Superscript minus signs on wavenumbers are lost ("4 cm 1" for 4 cm^-1) and are restored.
Figures 1 and 2 (total ion chromatograms), 3-6 (scanning electron micrographs), 7 and 8 (T2
relaxation distributions) and 9 and 10 (amide-I deconvolutions) are images and carry no tabulated
number that is not already in Tables 1-6. There is no supplementary material.

**Identity warning, resolved first.** This is **not** the paper the flavour-binding literature
usually means by "Guo 2020". Bi 2022 cites *Guo, J.; He, Z.; Wu, S.; Zeng, M.; Chen, J. (2020),
"Effects of concentration of flavor compounds on interaction between soy protein isolate and flavor
compounds", Food Hydrocolloids 100, Article 105388* — a Klotz-model binding study on soy protein
isolate. **The PDF on disk is a different paper by different authors**: Guo, Teng, Huang, Lv, Lv,
Babich, Yu, Li, Wang, Jiang, Food Hydrocolloids **105** (2020) **105752**. The volume and article
number in the task match the PDF, so this dossier describes the extrusion paper; **the Klotz soy
binding paper is a separate document and is not on disk** (Flags 1).

Repo status before this dossier: no `guo2020` source id appears in
`src/kinetic_core/parameters_matrix.py`, `data/lit/binding_constants.yml` or
`data/species/protein_matrices.yml`, and there is no extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of material characteristics on the structural characteristics and flavor substances retention of meat analogs" |
| Authors | Zengwang Guo (a, equal first), Fei Teng (a, equal first), Zhaoxian Huang, Bo Lv, Xiqiao Lv, Olga Babich (b), Wenhua Yu (c), Yang Li (a, corresponding), Zhongjiang Wang (a,d, corresponding), Lianzhou Jiang (a,c,d, corresponding) |
| Affiliations | a College of Food Science, Northeast Agricultural University, Harbin, Heilongjiang 150030, China; b Institute of Life Systems Research, Immanuel Kant Baltic Federal University, Kaliningrad, Russia; c Shandong Wan De Fu Industrial Group Co. Ltd., Dongying, Shandong 257500; d Linyi Yuwang Plant Protein Industry Co. Ltd., Dezhou, Shandong 253000 |
| Venue | Food Hydrocolloids 105 (2020) 105752. Received 8 September 2019, revised 4 February 2020, accepted 7 February 2020, online 9 February 2020 |
| DOI | 10.1016/j.foodhyd.2020.105752 |
| Funding | Heilongjiang Provincial Science Foundation C2018024; NSFC 31671807; Shandong Taishan Industry Leader Talent LJNY201607; Shandong key R&D 2018YYSP026 and 2018YYSP021; Heilongjiang GA17B002; China Postdoctoral 2018M641798; Heilongjiang Postdoctoral LBH-Z18030. **Two of the ten authors are employed by the two industrial partners that supplied the soy isolate (Yuwang Group / Linyi Yuwang Plant Protein) — declared affiliations, and the competing-interest statement says "No conflict of interest."** |
| The matrix | **soybean protein isolate, 94.2 % protein, 2.8 % moisture, 1.7 % carbohydrate, 0.6 % fat** (Yuwang Group Ltd., Shandong), blended with **wheat gluten, 85.9 % protein, 2.6 % moisture, 4.8 % carbohydrate, 2.6 % fat** (Zhengzhou Mindtek) |
| The flavour | **"Natural flavor powder (spice flavoring extract, β-cyclodextrin, 11 % fat)"** from Guangzhou Chenyi Trade Co. Ltd., dosed at **1 % of the dry mass of the soy isolate**. The sixteen compounds resolved are an anise/fennel/citrus spice profile, not Maillard products |
| Naming | "meat analogs" = the high-moisture extrudate; "retention rate (%)" = the ratio of a compound in the extrudate to the same compound in the unextruded raw material (the formula itself is never printed — Flags 3); T2b / T21 / T22 = LF-NMR spin-spin relaxation times of strongly / moderately / weakly bound protons, M2b / M21 / M22 their signal shares |
| Companions on disk | `bi2022_extraction.md` (pea, binding constants), `barallatperez2024_extraction.md` (lupin, in-mouth), `bornhorst2017_extraction.md` / `bornhorst2017b_extraction.md` (whey), and on the soy side `damodaran1981_extraction.md`, `ruan2014_extraction.md`, `shimada1988_extraction.md`, `xiao2024_extraction.md`, `xiao2025_extraction.md` (the last of which is the other extrusion source in the corpus) |

## 1. Why it matters

**Every number here is a NON-COVALENT retention, and it is not even a binding measurement.** The
paper is explicit about the mechanism it believes in — "In most cases, the interactions between
protein and flavor compounds are reversible, involving hydrophobic and hydrogen bonding" — and it
attributes its results to microstructure, mass-transfer resistance and steam stripping rather than
to chemistry: "The resistance to mass transfer is affected by the texture and microstructure of the
product, not the partitioning." **There is no covalent measurement anywhere in it**, so nothing here
touches `src/kinetic_core/matrix_sites.py`. Nor is there any equilibrium constant, so nothing here
can enter `REVERSIBLE_BINDING` in `src/kinetic_core/parameters_matrix.py` either.

**So where does it go?** It is a **process-retention** source. Its natural home is a benchmark or a
process-loss term, not a parameter registry: `data/lit/binding_constants.yml`'s
`percent_bound_at_conditions` record type is the closest existing shape, but it would be a misuse —
that field means "bound at equilibrium at a stated protein loading", and this paper's percentage
means "survived a 150 C extruder and a steam flash". The two are different physical quantities and
must not share a column.

**Is soy a protein the matrix table carries? Yes — `soy_isolate` is one of the three entries in
`data/species/protein_matrices.yml`** (free thiol 0.0067, disulfide 0.0457, amine 0.36 mmol per gram
of protein). And the isolate here is characterised the way that table wants: **94.2 % protein**. But
this paper measures no site density and no binding constant on it, so the pairing gains the matrix
table nothing directly.

**What it does give the repository, and it is worth having.** Three things the matrix layer does not
currently contain at all:

1. **A measured aldehyde survival fraction through a real thermal process.** The named aldehyde class
   retains **9.9-13.3 %** across the gluten series and **7.8-21.2 %** across the moisture series
   (Tables 1 and 2). The matrix layer today can only say what fraction of an aroma compound is bound
   at equilibrium in a vial; it has no term at all for what fraction survives the process. These are
   aromatic aldehydes, not the aliphatic hexanal-class the engine emits, so the number transfers by
   analogy only (Flags 6).
2. **A moisture axis.** Retention of the whole flavour load falls monotonically from
   **52.46 % at 50 % material moisture to 23.04 % at 80 %** — a **2.28x span (mine)** driven by one
   variable. `src/kinetic_core/acrylamide_conditions.py` carries water-activity windows on the
   Maillard side; the matrix layer has no moisture term, and this is a direct measurement of one on
   the retention side.
3. **A structural covariate series on soy.** Tables 5 and 6 print the FT-IR secondary-structure
   composition of the same extrudates whose retention is in Tables 1 and 2, and Tables 3 and 4 print
   the LF-NMR water distribution. So the paper offers a matched (retention, structure, water) triple
   on nine soy-based samples — the kind of covariate set a later wave would need to test whether a
   structural descriptor predicts retention at all. **The paper itself does not fit any such
   relation and reports no correlation coefficient.**

What this paper does NOT give the repository: any binding constant; any partition coefficient; any
protein loading in g/L; any water activity; any residence time; any product temperature (only barrel
set-points); any absolute concentration of any compound; any aliphatic aldehyde; any Maillard
product; any sensory measurement; any replicate-level data; and — critically — **any statement of
the formula by which "retention rate" was computed**.

## 2. Methods as they matter to a model

- **The proteins and their exact description.** Soybean protein isolate, **94.2 % protein, 2.8 %
  moisture, 1.7 % carbohydrate, 0.6 % fat**, from Yuwang Group Ltd., Shandong. Wheat gluten,
  **85.9 % protein, 2.6 % moisture, 4.8 % carbohydrate, 2.6 % fat**, from Zhengzhou Mindtek
  Biotechnology. **The abstract calls the soy material a "soy protein concentrate"; the Materials
  section calls it an isolate and gives 94.2 % protein, which is an isolate specification.** No
  further characterisation — no thiol assay, no lysine content, no molar mass, no solubility.
- **The flavour and how it was delivered.** "Natural flavor powder (spice flavoring extract,
  **β-cyclodextrin**, 11 % fat)", dosed at **1 % of the dry mass of the soy isolate**. **The flavour
  is a β-cyclodextrin inclusion powder**, i.e. the compounds enter the extruder already encapsulated
  in a host molecule whose entire purpose is to hold volatiles. This is never mentioned again and is
  never controlled for (Flags 4). No individual compound concentration is stated anywhere.
- **The blend.** Soy isolate + wheat gluten + flavour powder mixed with deionised water for 20 min.
  **Wheat gluten at 0, 10, 20, 30 and 40 % (by dry mass, added to the soy isolate); material moisture
  fixed at 50, 60, 70 and 80 %.** The design is one-factor-at-a-time around a shared centre point
  (Flags 5).
- **The extruder and the thermal programme.** FT-36 twin-screw extruder (Ji'nan Delun Machinery),
  **eight independently heated zones at 20, 50, 80, 150, 140, 100, 80 and 60 C** from zone 1 to
  zone 8. Barrel **L/D 26:1**, screw **compression ratio 4.6:1**, **die diameter 1.00 cm**. **The
  feed rate is stated twice and inconsistently: "fed to the extruder at a rate of 6 kg/h" and "The
  feed rate was kept constant at 30 g/min" (= 1.8 kg/h), a 3.3x discrepancy** (Flags 2). **No screw
  speed, no residence time, no measured melt temperature, no die pressure, no specific mechanical
  energy is reported** — so the thermal history of the flavour cannot be reconstructed. The peak
  set-point is **150 C**.
- **Method 1 — retention by headspace GC-MS (what it measures).** After Menis et al. 2013. The
  extraction head was activated for **1 h at 270 C** at the gas-phase inlet. **2 g of sample cut to
  2 mm**, and separately **the raw material for extrusion**, were put into vials; **1 uL of
  2,4,6-trimethylpyridine (TMP) internal standard** added; sealed; loaded into a **5 mL automatic
  headspace injector** (Hamilton-CTC 203182). Analysed on an **Agilent 7890A GC / 5975C MSD** with a
  fused-silica **Elite 5MS, 30 m x 0.25 mm x 1.4 um** column, helium at **1 mL/min**. Compounds
  accepted at a **library match above 80 %**, and quantified as **"the ratio of volatile flavor
  substances to the area of the TMP peak"**. **The retention rate itself is never defined by a
  formula.** The only reasonable reading is (TMP-normalised area in the extrudate) / (TMP-normalised
  area in the raw material) x 100, since the raw material is the only other thing put in a vial; that
  reading is an inference, not a statement (Flags 3). **The headspace temperature and equilibration
  time are not given** — only the SPME/head conditioning at 270 C. This is a **headspace abundance
  ratio, not a partition coefficient and not a mass balance**: nothing here weighs the compound.
- **Method 2 — scanning electron microscopy.** TM-5570 (Hitachi). Slices 1 x 4 x 4 mm, fixed in
  2.5 % glutaraldehyde at pH 7.4 (4 C, 1.5 h), washed twice in 0.1 M pH 7.4 phosphate (10 min),
  dehydrated through 50 / 70 / 90 % alcohol (10 min each) and 100 % ethanol twice (10 min each),
  then 1:1 ethanol/tert-butanol and tert-butanol (15 min), freeze-dried 4 h, gold-coated (1500 nm),
  imaged at 500x. Longitudinal and cross sections.
- **Method 3 — low-field NMR water distribution.** Bruker Mq20-NMR, **Carr-Purcell-Meiboom-Gill**
  sequence for the spin-spin relaxation time T2; three components T2b (strongly bound), T21
  (moderately bound), T22 (weakly bound) and their signal shares M2b, M21, M22. **Duplicate
  measurements.** This measures where the water is, not where the flavour is.
- **Method 4 — FT-IR secondary structure.** Bruker Vertex 70, **64 scans at 4 cm^-1 resolution,
  4000-400 cm^-1, 25 C**. Amide I deconvoluted with **PeakFit 4.12**, Gaussian fitting. Band
  assignments taken from the literature: alpha-helix **1646-1662 cm^-1**, antiparallel beta-sheet
  **1608-1622 and 1682-1700 cm^-1**, intramolecular beta-sheet **1622-1639 cm^-1**, beta-turn
  **1662-1681 cm^-1**, random coil **1637-1645 cm^-1**. Note that the random-coil and intramolecular
  beta-sheet windows **overlap** (1637-1639 cm^-1) as printed.
- **Statistics.** "All statistical analyses were performed in triplicate with duplicate samples."
  One-way ANOVA, SPSS 20.0, significance at **p <= 0.05**. Letters in every table compare **within a
  row** for Tables 1-4 and **within a column** for Tables 5-6 (the notes say so explicitly, and they
  differ — read them each time).

## 3. Tables re-typed

All values as printed, with the plus-or-minus signs restored. Retention rates are percentages.

### Table 1. "Retention rates of the volatile flavor substances in the meat analogs at different wheat gluten contents."

Note as printed: "Comparisons were carried out between values in the same row; values with different
letter(s) indicate a significant difference at p <= 0.05. **Samples 1-5 show the meat analogue
samples at 0, 10, 20, 30, and 40 % wheat gluten content, respectively.**"

| # | Volatile compound | 1 (0 % gluten) | 2 (10 %) | 3 (20 %) | 4 (30 %) | 5 (40 %) |
|---|---|---|---|---|---|---|
| 1 | β-Myrcene | 59.56 ± 1.15 b | 60.13 ± 1.34 b | 77.33 ± 1.29 c | 60.70 ± 1.49 b | 54.25 ± 1.57 a |
| 2 | 1-Methyl-4-isopropyl cyclohexadiene | 56.87 ± 1.54 ab | 59.06 ± 1.45 b | 64.08 ± 1.44 c | 58.96 ± 1.51 b | 54.25 ± 1.39 a |
| 3 | 3-Carene | 41.12 ± 1.58 a | 42.73 ± 1.49 a | 94.60 ± 1.34 b | 44.60 ± 1.62 a | 40.81 ± 1.48 a |
| 4 | Cinene | 36.79 ± 1.22 b | 39.44 ± 1.26 c | 46.48 ± 1.11 d | 40.80 ± 1.28 c | 29.17 ± 1.17 a |
| 5 | 1-methyl-4-acropropyl cyclohexane | 60.82 ± 1.29 b | 62.19 ± 1.16 b | 66.84 ± 1.12 c | 62.79 ± 1.07 b | 56.66 ± 1.15 a |
| 6 | 3,7-Dimethyl-1,3,6-octatriene | not detected | not detected | not detected | not detected | not detected |
| 7 | Caryophyllene | not detected | not detected | 9.79 ± 1.34 b | 4.31 ± 1.25 a | not detected |
| 8 | Estragole | 27.06 ± 1.43 c | 22.57 ± 1.32 b | 47.31 ± 1.47 d | 19.08 ± 1.11 a | 26.98 ± 1.26 c |
| 9 | Anethole | 95.42 ± 1.47 a | 91.78 ± 1.57 a | 91.46 ± 1.46 a | 91.30 ± 1.62 a | 95.50 ± 1.29 a |
| — | **Relative total retention of volatile flavors of total alkenes** | 39.39 ± 1.37 a | 41.96 ± 1.53 a | 41.99 ± 1.49 a | 55.32 ± 1.38 b | 42.50 ± 1.52 a |
| 10 | Maltol | 59.48 ± 1.58 a | 72.33 ± 1.61 b | 86.50 ± 1.59 c | 76.16 ± 1.43 b | 55.72 ± 1.29 a |
| 11 | Eugenol | 5.02 ± 1.06 a | 5.28 ± 1.52 a | 5.59 ± 1.37 a | 5.63 ± 1.26 a | 4.88 ± 1.29 a |
| — | **Relative total retention of volatile flavors of total phenols** | 40.72 ± 1.32 b | 32.25 ± 1.47 a | 38.81 ± 1.29 b | 46.05 ± 1.37 c | 40.90 ± 1.28 b |
| 12 | p-Anisaldehyde | 14.36 ± 0.84 b | 14.03 ± 0.91 b | 16.97 ± 0.92 c | 12.73 ± 0.79 a | 13.35 ± 0.88 ab |
| 13 | p-Isopropylbenzaldehyde | 11.38 ± 0.43 b | 11.54 ± 0.37 b | 9.69 ± 0.48 a | 11.61 ± 0.48 b | 10.29 ± 0.39 a |
| — | **Relative total retention of volatile flavors of total aldehydes** | 9.9 ± 1.21 a | 12.87 ± 1.17 a | 12.79 ± 1.31 a | 13.33 ± 1.22 a | 12.17 ± 1.57 a |
| 14 | Diethyl malonate | 73.41 ± 1.46 a | 87.46 ± 1.42 b | 92.04 ± 1.57 b | 90.63 ± 1.39 b | 70.81 ± 1.45 a |
| 15 | Cineole | not detected | not detected | not detected | not detected | not detected |
| 16 | 1-Methyl-4-(1-methylethyl)-Cyclohexane | 33.58 ± 1.46 a | 33.76 ± 1.52 a | 60.23 ± 1.49 b | 61.56 ± 1.62 b | 31.05 ± 1.52 a |
| — | **Relative total retention of volatile flavor compounds** | 35.46 ± 1.12 a | 42.04 ± 1.15 b | 44.07 ± 1.03 b | 43.78 ± 1.28 b | 33.39 ± 1.37 a |

### Table 2. "Retention rate of the volatile flavor substances in the meat analogs at different moisture contents of the raw materials."

Note as printed: "Comparisons were carried out between values of the same row ... **Samples 6-9 show
the meat analogue samples at 50, 60, 70, and 80 % moisture contents, respectively.**"

| # | Volatile compound | 6 (50 % moisture) | 7 (60 %) | 8 (70 %) | 9 (80 %) |
|---|---|---|---|---|---|
| 1 | β-Myrcene | 77.60 ± 1.37 c | 77.33 ± 1.42 c | 73.06 ± 1.19 b | 58.35 ± 1.24 a |
| 2 | 1-Methyl-4-isopropyl cyclohexadiene | 67.33 ± 1.19 c | 64.08 ± 1.22 b | 62.72 ± 1.47 b | 46.91 ± 1.36 a |
| 3 | 3-Carene | 95.27 ± 1.24 c | 94.60 ± 1.29 c | 58.50 ± 1.34 b | 34.63 ± 1.41 a |
| 4 | Cinene | 47.85 ± 1.40 c | 46.48 ± 1.51 c | 34.73 ± 1.22 b | 21.08 ± 1.33 a |
| 5 | 1-methyl-4-acropropyl cyclohexane | 68.04 ± 1.16 c | 66.84 ± 1.31 c | 23.76 ± 1.22 b | 16.92 ± 1.36 a |
| 6 | 3,7-Dimethyl-1,3,6-octatriene | not detected | not detected | not detected | not detected |
| 7 | Caryophyllene | 10.84 ± 0.89 b | 9.79 ± 0.79 b | 8.19 ± 0.97 ab | 6.82 ± 0.94 a |
| 8 | Estragole | 47.84 ± 1.07 c | 47.31 ± 1.14 c | 21.00 ± 1.21 b | 16.07 ± 1.19 a |
| 9 | Anethole | 92.07 ± 1.38 b | 91.46 ± 1.41 b | 11.73 ± 1.61 a | not detected |
| — | **Relative total retention of volatile flavors of total alkenes** | 56.32 ± 1.47 c | 55.32 ± 1.52 c | 32.63 ± 1.47 b | 22.31 ± 1.49 a |
| 10 | Maltol | 87.26 ± 1.02 c | 86.50 ± 1.01 c | 76.94 ± 1.20 b | 55.28 ± 1.37 a |
| 11 | Eugenol | 28.65 ± 0.98 d | 5.59 ± 0.55 c | 4.85 ± 0.72 b | 2.01 ± 0.68 a |
| — | **Relative total retention of volatile flavors of total phenols** | 57.96 ± 1.03 d | 46.05 ± 0.99 c | 40.90 ± 1.22 b | 28.65 ± 1.18 a |
| 12 | p-Anisaldehyde | 25.64 ± 1.12 d | 16.97 ± 0.57 c | 14.83 ± 0.69 b | 12.58 ± 0.81 a |
| 13 | p-Isopropylbenzaldehyde | 16.66 ± 0.96 d | 9.69 ± 0.84 c | 5.17 ± 0.76 b | 3.02 ± 0.68 a |
| — | **Relative total retention of volatile flavors of total aldehydes** | 21.15 ± 1.12 c | 13.33 ± 1.02 b | 10.00 ± 0.86 b | 7.80 ± 1.04 a |
| 14 | Diethyl malonate | 95.30 ± 1.49 c | 92.04 ± 1.51 c | 71.29 ± 1.24 b | 66.08 ± 1.17 a |
| 15 | Cineole | 11.64 ± 0.89 a | not detected | not detected | not detected |
| 16 | 1-Methyl-4-(1-methylethyl)-Cyclohexane | 67.41 ± 1.05 d | 60.23 ± 0.82 c | 49.02 ± 0.96 b | 28.83 ± 1.06 a |
| — | **Relative total retention of volatile flavor compounds** | 52.46 ± 1.12 d | 44.07 ± 1.06 c | 32.24 ± 1.22 b | 23.04 ± 1.03 a |

### Table 3. "Effects of wheat gluten content on the T2 relaxation time of the meat analogs."

| row | 0 % | 10 % | 20 % | 30 % | 40 % |
|---|---|---|---|---|---|
| T2b (ms) | 0.54 ± 0.14 a | 0.54 ± 0.09 a | 0.54 ± 0.10 a | 0.52 ± 0.13 a | 0.52 ± 0.23 a |
| T21 (ms) | 5.55 ± 0.56 c | 5.14 ± 0.40 c | 4.41 ± 0.35 b | 3.93 ± 0.23 a | 4.41 ± 0.36 b |
| T22 (ms) | 64.39 ± 2.72 c | 51.94 ± 2.10 b | 48.93 ± 1.17 b | 46.80 ± 1.16 a | 59.78 ± 1.18 c |
| Sample moisture content (%) | 49.67 ± 0.63 a | 53.24 ± 0.64 b | 56.19 ± 0.70 c | 56.81 ± 0.72 c | 57.05 ± 0.72 c |
| M2b (%) | 6.11 ± 0.14 a | 10.23 ± 0.16 b | 12.08 ± 0.13 c | 16.82 ± 0.19 d | 17.03 ± 0.19 d |
| M21 (%) | 91.29 ± 1.09 c | 85.70 ± 0.93 b | 85.98 ± 0.10 b | 82.25 ± 0.90 b | 80.35 ± 0.86 a |
| M22 (%) | 2.60 ± 0.04 c | 4.07 ± 0.05 d | 1.94 ± 0.02 b | 0.93 ± 0.01 a | 2.62 ± 0.03 c |

Footnote defines T2b / T21 / T22 as the relaxation of strongly / moderately / weakly bound protons
and M2b / M21 / M22 as their relaxation signal components.

### Table 4. "Effects of moisture content in the raw materials on the T2 relaxation time of the meat analogs."

| row | 50 % | 60 % | 70 % | 80 % |
|---|---|---|---|---|
| T2b (ms) | 0.55 ± 0.07 a | 0.54 ± 0.03 a | 0.55 ± 0.07 a | 0.56 ± 0.07 a |
| T21 (ms) | 3.08 ± 0.31 a | 4.41 ± 0.35 b | 6.11 ± 0.57 c | 7.39 ± 0.72 d |
| T22 (ms) | 59.22 ± 1.33 c | 48.93 ± 1.17 b | 47.09 ± 1.05 b | 36.04 ± 0.90 a |
| Sample moisture content (%) | 46.26 ± 0.52 a | 55.19 ± 0.70 b | 63.08 ± 0.72 c | 72.50 ± 0.80 d |
| M2b (%) | 13.89 ± 0.14 d | 12.08 ± 0.13 c | 10.01 ± 0.11 b | 5.21 ± 0.06 a |
| M21 (%) | 81.94 ± 0.79 a | 86.48 ± 0.88 b | 89.31 ± 0.90 c | 93.36 ± 0.95 d |
| M22 (%) | 4.17 ± 0.05 c | 1.44 ± 0.02 b | 0.68 ± 0.02 a | 1.43 ± 0.01 b |

### Table 5. "The relative percentage of protein secondary structure for the meat analogs at different wheat gluten contents."

Letters compare **down a column** here (the note says "between values of the same column").

| Wheat gluten (%) | antiparallel β-sheet (%) | Intramolecular β-sheet (%) | α-helix (%) | β-turn (%) | random coil (%) |
|---|---|---|---|---|---|
| 0 | 17.20 ± 0.16 ab | 24.24 ± 0.12 b | 15.35 ± 0.18 a | 28.36 ± 0.11 d | 14.85 ± 0.14 a |
| 10 | 16.85 ± 0.14 a | 24.95 ± 0.08 c | 16.05 ± 0.11 b | 26.30 ± 0.08 c | 15.85 ± 0.13 b |
| 20 | 17.45 ± 0.12 b | 25.41 ± 0.08 c | 14.90 ± 0.16 a | 26.23 ± 0.07 c | 16.01 ± 0.16 b |
| 30 | 19.48 ± 0.15 c | 23.06 ± 0.14 a | 17.25 ± 0.09 c | 23.19 ± 0.24 a | 17.02 ± 0.08 c |
| 40 | 17.51 ± 0.09 b | 25.41 ± 0.06 c | 15.98 ± 0.12 b | 24.07 ± 0.16 b | 17.03 ± 0.09 c |

### Table 6. "Relative percentage of protein secondary structure for the meat analogs with different moisture contents in the raw materials."

| Moisture (%) | antiparallel β-sheet (%) | Intramolecular β-sheet (%) | α-helix (%) | β-turn (%) | random coil (%) |
|---|---|---|---|---|---|
| 50 | 16.03 ± 1.02 a | 29.50 ± 0.08 c | 15.81 ± 0.13 b | 24.81 ± 0.12 a | 13.88 ± 0.11 a |
| 60 | 17.45 ± 1.06 b | 25.41 ± 0.11 b | 14.90 ± 0.08 b | 26.23 ± 0.11 b | 16.01 ± 0.08 c |
| 70 | 19.47 ± 0.89 c | 20.96 ± 0.09 a | 13.65 ± 0.11 a | 27.89 ± 0.08 c | 18.03 ± 0.13 d |
| 80 | 20.03 ± 0.95 d | 25.58 ± 0.14 b | 14.85 ± 0.09 b | 24.59 ± 0.09 a | 14.96 ± 0.10 b |

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| barrel zone set-points, zones 1-8 | 20 / 50 / 80 / **150** / 140 / 100 / 80 / 60 C | §2.2 |
| barrel L/D, screw compression ratio, die diameter | 26:1, 4.6:1, 1.00 cm | §2.2 |
| feed rate | **6 kg/h** in one sentence and **30 g/min** in another | §2.2 (Flags 2) |
| flavour dose | 1 % of the dry mass of the soy isolate | §2.2 |
| SPME head activation | 1 h at 270 C | §2.3 |
| sample for GC-MS | 2 g cut to 2 mm; 1 uL TMP internal standard; 5 mL headspace injector | §2.3 |
| library match threshold | above 80 % | §2.3 |
| FT-IR band windows | α-helix 1646-1662; antiparallel β-sheet 1608-1622 and 1682-1700; intramolecular β-sheet 1622-1639; β-turn 1662-1681; random coil 1637-1645 cm^-1 | §3.4 |
| T21 signal share | "accounted for > 80 % of the total signal" | §3.3 |
| relaxation-time interpretation | components below 1 ms are mainly non-exchangeable CH protons | §3.3 |
| the abstract's class ordering | "esters, alkanes, alkenes, phenols, aldehydes, and alcohols", largest to smallest | Abstract and Conclusions |

**Everything else in this paper is figure-only**: the two total-ion chromatograms (Figs. 1, 2), the
four SEM series (Figs. 3-6), the two T2 relaxation distributions (Figs. 7, 8) and the two amide-I
deconvolutions (Figs. 9, 10). No concentration in any unit of mass appears anywhere in the paper.

### Arithmetic on the printed numbers (all mine)

**1. The two tables share a sample, and it does not reconcile.** Table 1's column 3 (20 % gluten) and
Table 2's column 7 (60 % moisture) are the same extrudate — the centre point of the one-factor-at-a-
time design. **All fourteen detected individual-compound values agree exactly** (β-myrcene 77.33,
cyclohexadiene 64.08, 3-carene 94.60, cinene 46.48, the cyclohexane 66.84, caryophyllene 9.79,
estragole 47.31, anethole 91.46, maltol 86.50, eugenol 5.59, p-anisaldehyde 16.97,
p-isopropylbenzaldehyde 9.69, diethyl malonate 92.04, the methylethyl-cyclohexane 60.23), and so
does the grand total (44.07). **The standard deviations do not**: β-myrcene is ± 1.29 in Table 1 and
± 1.42 in Table 2, cinene ± 1.11 against ± 1.51, eugenol ± 1.37 against ± 0.55 (a **2.5x**
difference on the same sample). **And the three class sub-totals disagree outright**: Table 2's
column 7 prints alkenes 55.32, phenols 46.05 and aldehydes 13.33, which are Table 1's **column 4**
(30 % gluten) values, not column 3's (41.99, 38.81, 12.79). One of the two tables has its class
sub-total row shifted by one column (Flags 7).

**2. The class sub-totals are not means of their listed members.** For column 1, the seven detected
"alkene" rows average **53.95 %** against a printed total of **39.39 %**; the two "phenol" rows
average 32.25 against a printed 40.72; the two "aldehyde" rows average 12.87 against a printed 9.9.
**So "relative total retention" is computed over the whole chromatographic class, including
compounds that are not individually listed, and the sixteen named compounds are a subset.** The
paper never says how many compounds each class total covers. (The aldehyde total for column 2,
12.87, happens to equal the two listed aldehydes' mean for column 1 exactly — another sign of a
one-column shift.)

**3. The moisture effect, sized.** Grand total retention falls **52.46 -> 44.07 -> 32.24 -> 23.04 %**
across 50 -> 80 % material moisture: a **2.28x** fall overall, and it is monotone in every one of the
sixteen rows that is detected at more than one level. Per ten points of moisture the factor is
**1.32x (mine, geometric)**. The largest single-compound sensitivity is **eugenol, 28.65 -> 2.01 %,
a 14.3x fall**; the smallest among the well-detected is **β-myrcene, 77.60 -> 58.35 %, 1.33x**.
**Anethole falls off a cliff between 60 and 70 % moisture** (91.46 -> 11.73, a **7.8x** step in one
level) and is undetected at 80 %.

**4. The gluten effect, sized, and it is not monotone.** Grand total **35.46 -> 42.04 -> 44.07 ->
43.78 -> 33.39 %** across 0 -> 40 % gluten: a rise to a maximum near 20-30 % and then a fall.
0 -> 20 % gains **1.24x**; 20 -> 40 % loses **1.32x**. The whole span, best to worst, is **1.32x** —
**less than a fifth of the moisture effect on a log scale (mine: log(1.32)/log(2.28) = 0.34, so a
third)**. The paper's own explanation for the fall is mechanical: above 30 % gluten "a layered
structure and surface fractures gradually appeared", and "the surface rupture of the extrudate ...
results in a loss of some volatile substances due to volatilization."

**5. The aldehyde class is the worst-retained named class in the paper.** Across all nine samples the
aldehyde total runs **7.80 to 21.15 %**, against alkenes 22.31-56.32 %, phenols 28.65-57.96 % and a
grand total of 23.04-52.46 %. **Aromatic aldehydes lose four fifths of their load through this
process at the centre point.** The abstract's stated ordering puts aldehydes second-to-last, above
alcohols only — and **no alcohol is listed in either table** (Flags 8).

**6. Protein content of the extrudate (mine, and it is the only route to a loading).** At the centre
point (20 % gluten, 60 % material moisture), the dry blend is 1 part soy isolate at 94.2 % protein
plus 0.2 part gluten at 85.9 %, so the dry basis is **(94.2 + 0.2 x 85.9)/1.2 = 92.8 % protein**
(the 1 % flavour powder is ignored). The measured sample moisture at that point is **55.19 %**
(Table 4), so the extrudate is **0.4481 x 0.928 = 41.6 % protein by wet mass = 416 g protein per kg
of extrudate (mine)**. Across the whole design the figure runs from **~50 % protein at 50 %
material moisture** down to **~26 % at 80 %**. **This is not a g/L in the sense
`MatrixLoading.protein_g_per_l` means it** — the extrudate is a solid, not a solution, and no
density is reported — but it establishes that this matrix is **forty to fifty times more concentrated
in protein than any entry in `MATRIX_LOADING`** (the densest of which is `soy_paste_hong` at a
declared 100-200 g/L). Any comparison of this paper's retention with a dilute-solution binding
constant must carry that.

**7. The vapour-pressure ordering is not the retention ordering.** Anethole and 3-carene, the two
best-retained compounds at low moisture (91-95 %), are also two of the more volatile in the set;
eugenol, the worst (2-29 %), is among the least volatile. **Whatever governs retention here, it is
not simple volatility**, which is consistent with the paper's own mass-transfer-resistance argument
and inconsistent with a partition picture. The paper prints no physicochemical property for any
compound, so this cannot be quantified from the paper alone.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **None of the sixteen compounds is in the
registry.** The registry's 75-odd ids are Maillard and lipid-oxidation markers (pyrazines,
furans, thiols, aliphatic aldehydes, acrylamide, HMF, glycation markers); this paper's set is a
spice-extract profile — monoterpenes (β-myrcene, 3-carene, cinene, the p-menthane-type cyclohexanes
and cyclohexadiene), a sesquiterpene (caryophyllene), two phenylpropanoid ethers (estragole,
anethole), a pyranone (maltol), a phenol (eugenol), two aromatic aldehydes (p-anisaldehyde,
p-isopropylbenzaldehyde), a diester (diethyl malonate) and an ether (cineole). **`parameters_matrix.py`'s
`COMPOUND_STRUCTURE` carries none of them either**, and its `binding_class` vocabulary has no
terpene, no phenylpropanoid and no aromatic-aldehyde class. Registering any of these means extending
both registries, and `hdmf` (furaneol) is the only entry in `compounds.yml` even structurally near
maltol.

Every row below shares: **soy protein isolate (94.2 % protein) blended with wheat gluten (85.9 %
protein) and 1 % w/w of a β-cyclodextrin spice-flavour powder**, high-moisture twin-screw extrusion
with a **150 C peak zone set-point**, die 1.00 cm, retention measured against the unextruded raw
material by headspace GC-MS with a TMP internal standard, triplicate with duplicate samples.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **retention of the total flavour load vs material moisture** | **52.46 ± 1.12 / 44.07 ± 1.06 / 32.24 ± 1.22 / 23.04 ± 1.03** | % of the raw-material level | 50 / 60 / 70 / 80 % material moisture, 20 % wheat gluten, 150 C peak | Table 2, last row, p. 5 | **retention_percent** |
| **retention of the total flavour load vs wheat gluten** | **35.46 ± 1.12 / 42.04 ± 1.15 / 44.07 ± 1.03 / 43.78 ± 1.28 / 33.39 ± 1.37** | % | 0 / 10 / 20 / 30 / 40 % gluten, 60 % moisture | Table 1, last row, p. 5 | **retention_percent** |
| **retention of the ALDEHYDE class vs material moisture** | **21.15 ± 1.12 / 13.33 ± 1.02 / 10.00 ± 0.86 / 7.80 ± 1.04** | % | as above | Table 2, p. 5 | **retention_percent** — the class of interest, but **aromatic** aldehydes (Flags 6) |
| **retention of the ALDEHYDE class vs wheat gluten** | **9.9 ± 1.21 / 12.87 ± 1.17 / 12.79 ± 1.31 / 13.33 ± 1.22 / 12.17 ± 1.57** | % | as above | Table 1, p. 5 | **retention_percent** — no cell differs significantly from any other (all lettered "a") |
| retention of p-anisaldehyde vs moisture | **25.64 ± 1.12 / 16.97 ± 0.57 / 14.83 ± 0.69 / 12.58 ± 0.81** | % | as above | Table 2 | retention_percent |
| retention of p-isopropylbenzaldehyde vs moisture | **16.66 ± 0.96 / 9.69 ± 0.84 / 5.17 ± 0.76 / 3.02 ± 0.68** | % | as above | Table 2 | retention_percent |
| retention of the alkene class vs moisture / vs gluten | **56.32 / 55.32 / 32.63 / 22.31** and **39.39 / 41.96 / 41.99 / 55.32 / 42.50** | % | as above | Tables 2 and 1 | retention_percent (with the column-shift caveat, Flags 7) |
| retention of the phenol class vs moisture / vs gluten | **57.96 / 46.05 / 40.90 / 28.65** and **40.72 / 32.25 / 38.81 / 46.05 / 40.90** | % | as above | Tables 2 and 1 | retention_percent (same caveat) |
| retention of the sixteen individual compounds | see Tables 1 and 2 above, re-typed in full | % | as above | Tables 1 and 2, p. 5 | retention_percent |
| **moisture sensitivity of retention** | **2.28x** loss from 50 % to 80 % material moisture | — | 20 % gluten | 52.46/23.04 (mine) | **within_study_ratio** — the single most transferable statement in the paper |
| **gluten sensitivity of retention** | **1.32x** best-to-worst, non-monotone with a maximum at 20-30 % | — | 60 % moisture | 44.07/33.39 (mine) | **within_study_ratio** |
| eugenol moisture sensitivity | **14.3x** | — | 50 -> 80 % moisture | 28.65/2.01 (mine) | within_study_ratio — the extreme of the set |
| β-myrcene moisture sensitivity | **1.33x** | — | as above | 77.60/58.35 (mine) | within_study_ratio — the other extreme |
| anethole cliff | **7.8x** in one moisture step | — | 60 -> 70 % moisture | 91.46/11.73 (mine) | within_study_ratio |
| aldehyde class vs alkene class at the centre point | **13.33 against 55.32** (a **4.15x** disadvantage, mine) | — | 60 % moisture, 20 % gluten | Table 2 | within_study_ratio |
| measured sample moisture, gluten series | **49.67 ± 0.63 / 53.24 ± 0.64 / 56.19 ± 0.70 / 56.81 ± 0.72 / 57.05 ± 0.72** | % w/w of the extrudate | 0-40 % gluten | Table 3, p. 6 | **level_only** — note the sample moisture RISES with gluten even though material moisture was fixed |
| measured sample moisture, moisture series | **46.26 ± 0.52 / 55.19 ± 0.70 / 63.08 ± 0.72 / 72.50 ± 0.80** | % w/w | 50-80 % material moisture | Table 4, p. 6 | level_only — the extrudate always holds **less** water than the material was fed at |
| LF-NMR T2b / T21 / T22 and M2b / M21 / M22, both series | Tables 3 and 4 above | ms and % of signal | as above | Tables 3, 4, p. 6 | level_only (a water-mobility covariate, not a flavour measurement) |
| FT-IR secondary structure, both series | Tables 5 and 6 above | % of amide I | as above | Tables 5, 6, pp. 7-8 | level_only (a structural covariate) |
| protein content of the extrudate at the centre point | **~416** | g protein per kg of extrudate | 20 % gluten, 60 % material moisture | (mine) from Materials + Table 4 | **derived_assumption** — no density is printed, so this cannot be turned into g/L |
| soy isolate composition | **94.2 % protein, 2.8 % moisture, 1.7 % carbohydrate, 0.6 % fat** | % | as supplied | Materials, §2.1 | level_only — supplier's figures |
| wheat gluten composition | **85.9 % protein, 2.6 % moisture, 4.8 % carbohydrate, 2.6 % fat** | % | as supplied | Materials, §2.1 | level_only — supplier's figures |
| binding constant of any kind | — | — | — | **absent** | **the paper contains none** |

### Can these be put on the same basis as the shipped constants? No, and the reason is instructive.

**(a) There is no equilibrium here.** `REVERSIBLE_BINDING` multiplies a per-gram constant by a
loading to shift an air/matrix partition at a fixed temperature. This paper's retention is a
**survival fraction through an open, non-isothermal, steam-flashing process**. The two are not the
same quantity and no arithmetic converts one into the other. Attempting the registry's
K_g = (ratio − 1)/g would require a water baseline that does not exist here (the control is the
unextruded blend, not water) and a loading in g/L that the paper does not support.

**(b) The right home is a process-loss term or a benchmark, and neither exists yet.** The layer's
declared output is formulation-vs-formulation ratios; this paper supplies exactly that kind of
comparison (moisture 2.28x, gluten 1.32x) on nine matched samples. That is usable as a **within-study
ratio for a process term**, and it is the only such measurement in this five-paper cluster.

**(c) The matrix is forty to fifty times more concentrated than anything in `MATRIX_LOADING`.**
~416 g protein per kg of extrudate against `soy_paste_hong`'s declared 100-200 g/L and
`skim_milk`'s 33.9 g/L. Any binding-based reading of this paper's retention would sit far outside
the concentration range every shipped constant was measured in.

**(d) Nothing goes to `matrix_sites.py`.** No rate, no activation energy, no thiol, no amine, no
adduct. The paper mentions the covalent idea once, second-hand — "Proteins with higher contents of
lysine, arginine, and cysteine can potentially exhibit higher flavor binding capacities, as more
covalent bonds are involved (Wang & Arntfield, 2016)" — and measures nothing about it.

## 5. Flags

1. **This is not the "Guo 2020" the binding literature cites.** Bi 2022's reference list contains
   *Guo, J.; He, Z.; Wu, S.; Zeng, M.; Chen, J. (2020), Food Hydrocolloids 100, Article 105388,
   "Effects of concentration of flavor compounds on interaction between soy protein isolate and
   flavor compounds"* — a Klotz binding study on SPI, and the natural companion to this cluster.
   **That paper is not on disk.** The PDF here is Food Hydrocolloids **105**, article **105752**, a
   different study by a different group. **Request the 105388 paper**; if the matrix layer wants a
   soy non-covalent binding constant from 2020, that is where it is.
2. **The feed rate is stated twice and the two statements differ by 3.3x**: "the material blend was
   fed to the extruder at a rate of 6 kg/h" and, four sentences later, "The feed rate was kept
   constant at 30 g/min" (1.8 kg/h). Combined with the absence of any screw speed, residence time,
   melt temperature, die pressure or specific mechanical energy, **the thermal and shear history of
   the flavour is unreconstructable**. The only thermal number is the set-point ladder, whose peak
   is 150 C — a set-point, not a product temperature.
3. **The retention formula is never printed.** §2.3 describes the GC-MS and says the quantity is
   "the ratio of volatile flavor substances to the area of the TMP peak", and that both the extrudate
   and "the raw material for extrusion" were run. The retention percentage must be the ratio of those
   two normalised areas, but **the paper never says so**, never states the headspace equilibration
   temperature or time, and never reports an absolute concentration. Values above 90 % (anethole
   95.42, 3-carene 95.27) and the SPME head conditioning at **270 C** both deserve scrutiny against
   whatever formula was actually used.
4. **The flavour was dosed as a β-cyclodextrin inclusion powder and this is never controlled for.**
   β-cyclodextrin is a molecular host whose function is to trap volatiles; a "retention" measured on
   a cyclodextrin-delivered flavour is partly a measurement of the cyclodextrin's survival, not the
   protein matrix's. The powder also carries 11 % fat, another partition phase. **Neither appears
   again after the Materials section.** This is the single largest confounder in the paper: the
   quantity called "protein-matrix retention" has at least three unresolved phases in it.
5. **The design is one-factor-at-a-time with no factorial cell and no replication of the centre
   point across tables.** Nine extrudates: five gluten levels at 60 % moisture, four moisture levels
   at 20 % gluten, sharing one centre point. No gluten x moisture interaction can be estimated, and
   the paper's own results suggest one exists (gluten helps at 60 % moisture by raising structural
   tightness; at 80 % moisture the structure is already tight and the flavour is gone anyway).
6. **The two aldehydes are AROMATIC, not the aliphatic class the engine emits.** p-Anisaldehyde
   (4-methoxybenzaldehyde) and p-isopropylbenzaldehyde (cuminaldehyde) are conjugated aryl aldehydes;
   the repository's aldehyde species (`HEXANAL`, `NONANAL`, `DECADIENAL`, `FUR`) are aliphatic or
   furanic. `matrix_sites.py` treats furfural's aldehyde with the hexanal bracket "declared", and the
   same declaration would be needed here — **an aromatic aldehyde's Schiff-base reactivity and its
   volatility both differ from hexanal's**. Carry the 7.8-21.2 % survival as an analogy, never as a
   hexanal number.
7. **Table 1 and Table 2 share a sample and their class sub-totals disagree.** Column 3 of Table 1
   and column 7 of Table 2 are both "20 % gluten, 60 % moisture", and all fourteen individual
   compound values and the grand total match exactly — but the three class sub-totals printed in
   Table 2 column 7 (alkenes 55.32, phenols 46.05, aldehydes 13.33) are **Table 1 column 4's values**
   (30 % gluten), not column 3's (41.99, 38.81, 12.79). **One of the two tables has its sub-total
   rows shifted by one column.** Separately, the standard deviations on the shared sample differ
   between the two tables — up to **2.5x** on eugenol (± 1.37 vs ± 0.55) — which cannot happen if
   the same measurement is being reported twice. Until this is resolved, use the individual-compound
   rows (which agree) and treat the class sub-totals as suspect.
8. **The class labels do not match their members, and one class in the abstract does not exist in
   the tables.** The "total alkenes" block contains estragole and anethole, which are phenylpropanoid
   ethers, and caryophyllene, a sesquiterpene. The "total phenols" block contains maltol, a
   pyranone. The abstract and conclusions rank six classes "esters, alkanes, alkenes, phenols,
   aldehydes, and alcohols" — **no alcohol appears anywhere in Tables 1 or 2**, and the "ester"
   class rests on one compound, diethyl malonate. The class sub-totals are also **not** means of
   their listed members (section 3, item 2), so each class covers compounds the tables never name.
9. **The conclusions contradict the abstract on the paper's central claim.** The abstract says the
   gluten and moisture contents "affected the flavor characteristics of the meat analogs"; the
   conclusions say the class ordering "is not affected by the wheat gluten and moisture content".
   The two can be reconciled (the ordering is stable while the levels move), but as printed the
   second sentence reads as a denial of the first.
10. **No water activity, and moisture is not a_w.** The engine's Maillard side is indexed to water
    activity (`acrylamide_conditions.py`); this paper reports **moisture content** at 46-72 % w/w,
    which at those levels means a_w near 0.99 throughout but is not measured. Do not convert.
11. **"As moisture content increased from 60 % to 90 %"** (§3.2, describing the SEM series) — the
    experiment has no 90 % arm. A typographic error for 80 %.
12. **Two of the ten authors are employed by the supplier of the soy isolate** (Shandong Wan De Fu
    Industrial Group; Linyi Yuwang Plant Protein Industry, and the isolate came from Yuwang Group
    Ltd.), and the competing-interest statement reads "No conflict of interest." Recorded, not
    weighted.
13. **Replication is stated ambiguously**: "All statistical analyses were performed in triplicate
    with duplicate samples", while the NMR section says "Each measurement was performed in
    duplicate". Whether n = 3, 6 or 2 for any given table is not resolvable.
14. **What this paper does NOT contain**: any binding constant, partition coefficient or equilibrium
    quantity; any absolute concentration in mass units; any aliphatic aldehyde, ketone, pyrazine,
    thiol or other Maillard product; any protein loading in g/L; any water activity; any residence
    time or measured product temperature; any correlation between the structural covariates and the
    retention; any sensory measurement; any odour threshold; any supplementary material.
15. **What to request from the authors**: (i) the formula and the headspace conditions behind the
    retention percentages; (ii) the composition of the flavour powder, especially the cyclodextrin
    fraction and the absolute dose of each compound; (iii) which is right, 6 kg/h or 30 g/min, plus
    screw speed and residence time; (iv) how many compounds each class sub-total covers; (v) the
    resolution of the Table 1 / Table 2 sub-total column shift and of the differing standard
    deviations on the shared sample; (vi) a measured melt or die temperature.
16. **Registry gaps.** **None of the sixteen compounds is keyed in `data/keys/compounds.yml`**, and
    none has a `COMPOUND_STRUCTURE` entry in `parameters_matrix.py` — which also lacks any
    `binding_class` for terpenes, phenylpropanoids or aromatic aldehydes. `data/species/protein_matrices.yml`
    has `soy_isolate` (so the protein is known) but **there is no extruded-soy matrix in
    `MATRIX_LOADING`**, and the one that exists for soy (`soy_paste_hong`) is a paste at a declared
    100-200 g/L against this extrudate's ~416 g/kg. A companion source already in the corpus,
    `xiao2025_extraction.md`, covers low-moisture extrusion of soy blends and free-thiol survival;
    this paper is the high-moisture counterpart, and the two should be read together.
