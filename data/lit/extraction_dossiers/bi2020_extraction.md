# Bi 2020 — EXTRACTION (raw pea flour vs pea roasted whole at 160 C for 30 min then milled; HS-SPME GC-MS/GC-O, 26 odorants quantified against authentic standards in a deodorized pea-flour matrix, with odour thresholds and OAVs; aroma recombination and 33 omission models)

### The carried-in level and the cooked level of the same volatile in the same pea, printed side by side: hexanal 1260 ug/kg in raw pea flour and 324 ug/kg after 160 C / 30 min — real concentrations from matrix-matched five-point calibration curves, in a FLOUR (whole milled seed), not in an isolate.

**Source on disk:** `data/articles/bi2020.pdf` (10 pp., J. Agric. Food Chem. 2020, 68, 2718-2727).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/bi2020.txt`), whole file. Tables 1, 2, 3, 4 and 5 all came through clean and
are re-typed in full below; page 6 (journal p. 2723, the odour-threshold column of Table 4) was
additionally rendered at 200 dpi and read from the image to confirm two threshold values against
which the printed OAVs do not reconcile (Flags 3). Figure 1 (heat map of the 73 volatiles) and
Figure 2 (sensory spider diagrams, a and b) are images and are **figure-only**. Supporting
Information exists (effect of SPME extraction time, temperature and sample weight on total volatile
area and count) and is **not on disk**. Repo status before this dossier: `bi2026_extraction.md` is
on disk from the same first author; this 2020 paper had no dossier and is cited nowhere in `src/`,
`data/benchmarks/` or `results/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Characterization of Key Aroma Compounds in Raw and Roasted Peas (*Pisum sativum* L.) by Application of Instrumental and Sensory Techniques" |
| Authors | Shuang Bi, Xinxing Xu, Dongsheng Luo, Fei Lao, Xueli Pang, Qun Shen, Xiaosong Hu, Jihong Wu (corresponding, wjhcau@hotmail.com) — College of Food Science and Nutritional Engineering, China Agricultural University, Beijing; National Engineering Research Center for Fruit & Vegetable Processing |
| Venue | J. Agric. Food Chem. 2020, 68 (10), 2718-2727. Received 5 Dec 2019, revised 2 Feb 2020, accepted 4 Feb 2020, published 4 Feb 2020 |
| DOI | 10.1021/acs.jafc.9b07711 |
| Material | Peas, **Zhongwan variety**, grown at Datong, Shanxi Province, China (40 deg 08' N, 113 deg 24' E), harvested 2017, stored 4 C. Experiment run 2018-2019. |
| Naming | "raw peas" = whole dry seed milled to flour through a 0.30 mm screen. "roasted peas" = whole seed roasted 160 C / 30 min in a forced-air oven, then milled the same way. Both are analysed **as flour**. "AI" = aroma intensity, a 0-4 sniffing score. "OAV" = concentration / odour threshold in water. |
| Funding | National Key R&D Plan 2017YFD0401202; NSFC 31601483 |
| Data | no data-availability statement |

### How this relates to `bi2026_extraction.md`

Same first author (Shuang Bi), and both papers work on *Pisum sativum* of the Zhongwan family. They
are otherwise **different studies in different laboratories on different materials**, and they are
complementary rather than overlapping:

| | Bi 2020 (this paper) | Bi 2026 (`Bi2026.pdf`) |
|---|---|---|
| affiliation printed | China Agricultural University, Beijing (Jihong Wu group) | Beijing Technology and Business University, Flavor Science Laboratory (Ye Liu group) |
| material | dry pea seed milled to **flour**; and the same seed roasted 160 C / 30 min then milled | **pea milk**, 1:7 (w/w) seed:water slurry from a soymilk machine |
| cultivar as printed | "Zhongwan variety" | "Zhongwan No. 6 (ZW.6)" |
| the axis studied | **heat** (raw vs roasted, one condition) | **grinding time** (0 s to 8 min, eight quench points) |
| hexanal | **printed as a concentration**, 1260 ug/kg raw and 324 ug/kg roasted (Table 4) | **figure-only** (Fig. 1); the tables carry the lipid substrate side instead |
| lipid substrate | not measured | free linoleic acid 10.658 -> 0.309 ug/g over 8 min (Table 2 there) |
| what it gives Programme 7 | the **level** of the carried volatile, in a flour | the **shape in time** of the enzymatic step that makes it, in a milk |

So 2026 supplies the time course whose product level 2020 supplies — but in two different pea
matrices (flour vs milk), so they cannot be spliced into one series without an assumption. Neither
paper measures a protein isolate.

## 1. Why it matters

`tasks/roadmap_for_scientists.md` section 5d, Programme 7 part (ii), asks that "the isolate's own
volatiles enter as declared inputs with their measured levels and bands ... so `predict` reports
them as carried, with the matrix binding applied, rather than refusing". The two sources the
roadmap names for that levels table are Fischer 2021 and Zhang 2020b
(`fischer2021_extraction.md`, `zhang2020b_extraction.md`). Both of those measure **wet** pea
material: Zhang 2020b a 2 % protein pea milk (hexanal 164.18 ug/L), Fischer 2021 a pea protein
preparation in its own conditions. **This paper is the first on disk to print a carried-volatile
level in dry pea FLOUR, in ug/kg of flour**, which is the basis a dry recipe would use.

What it adds to that levels table, specifically:

1. **A raw, unheated, dry-basis hexanal level**: 1260 +/- 114 ug/kg in raw pea flour. This is the
   number the roadmap's argument needs — the model charges no hexanal at t = 0, and the panel's
   hexanal rows are under-predicted (`docs/assets/thiol_sink/24_fat_path_hexanal.png`, generated by
   `scripts/generators/build_story_figures.py`). A charge of order 1 mg/kg of flour is what this
   paper says a pea ingredient walks in with.
2. **A cooked counterpart at a stated cook**: 324 +/- 34.3 ug/kg after 160 C / 30 min of dry roast,
   i.e. hexanal **falls to 26 % of its raw level** across a real cook (mine, from Table 4). If the
   model charges hexanal as a declared input it must also be allowed to lose it; this is a measured
   loss factor at a stated temperature and time, from one laboratory, in one material.
3. **Nine other carried volatiles with levels in the same units and the same material**: nonanal,
   (E)-2-octenal, 1-hexanol, 1-pentanol, (Z)-2-penten-1-ol, benzyl alcohol, benzaldehyde,
   3-methylbutanoic acid, plus the raw/roasted pair for several of them.
4. **The Maillard side of the same cook, quantified**: nine pyrazines (243 to 16,100 ug/kg), maltol
   (103,000 ug/kg), ethyl maltol, furaneol, furfural, 3-methylbutanal, phenylacetaldehyde,
   guaiacol, 4-vinylguaiacol and dimethyl sulfide, all in the roasted flour. Twelve of these have
   registry ids. This is a **whole-seed roasting benchmark candidate** with a stated temperature
   and time, of the sort `results/validation/data_wishlist.md` asks for.
5. **Odour thresholds in water for 26 compounds, with their bibliographic source per compound**
   (van Gemert 2011; Jelen 2013; Preininger 2008). `data/species/off_flavour_targets.yml` carries
   a standing note that its thresholds "carry no verified citation and should be treated as
   uncited compilation values"; this table is a cited compilation for six of them and for twenty
   more.

What it does **not** give: no isolate, no concentrate, nothing with time resolution finer than
"raw vs 30 min", no rate constant, no moisture or water-activity statement, no lipid or free-fatty
-acid measurement, and no replicate cook at a second temperature.

## 2. Methods as they matter to a model

- **The material, exactly as described.** "Peas (Zhongwan variety) were cultivated at Datong in
  Shanxi Province, China ... and harvested in 2017 ... The peas were stored at 4 C after removing
  the impurities. Raw peas (500 g) were milled to flour using a coffee bean grinder (Joyoung) for
  30 s and then passing them through a 0.30 mm mesh screen." So the raw material is a **whole-seed
  flour**: cotyledon plus whatever passes 0.30 mm, at its native protein and lipid content. It is
  **not a protein isolate, not a protein concentrate, and not defatted**. No moisture content and
  no proximate composition are printed anywhere in the paper.
- **The cook.** "The roasted peas (500 g) were prepared using an Isotemp forced-air oven (Haier) at
  **160 C for 30 min**. These conditions were selected according to those commonly used in the
  processing of pea breads and were determined by a preliminary experiment to give a typical
  roasted flavor. After being roasted, the samples were milled as described for raw pea flour."
  The seeds are roasted **whole** and milled after; the raw sample is milled without roasting. No
  temperature log, no come-up-time correction, no mass-loss (moisture) figure across the roast.
  Because moisture is lost in a 160 C / 30 min roast and is not measured, the raw-to-roasted
  comparison is **not on a constant dry basis** (Flags 2).
- **Freshness.** "The raw and roasted peas were newly milled and roasted, respectively, before each
  experiment because the aroma changed over storage time." An explicit statement that the carried
  level is time-sensitive, with no number attached to it.
- **Extraction.** HS-SPME, 50/30 um DVB/CAR/PDMS fibre, 20 mm. **1.8 g NaCl dissolved in 5 mL
  distilled water in a 20 mL vial, then 1.5 g of raw or roasted pea flour, then 20 uL of internal
  standard.** Internal standard: 2,4,6-trimethylpyridine, 10 uL diluted 2 x 10^4-fold in methanol.
  Equilibration 50 C / 20 min; extraction 50 C / 50 min with intermittent shaking (stop 2 s after
  every 20 s); desorption at the 250 C splitless inlet, 5 min, fibre 20 mm in.
  **So every level in this paper is a headspace level from a rehydrated 1.5 g flour + 5 mL brine
  slurry at 50 C, back-calculated to the flour** — the matrix binding that Programme 7 wants to
  apply is already partly inside these numbers (Flags 5).
- **GC-MS.** Agilent 7890A + 5975 MSD; DB-5MS and HP-WAX, both 30 m x 0.25 mm x 0.25 um; 45 C
  (2 min) -> 240 C at 6 C/min, hold 5 min; He 1.0 mL/min. **Full scan (35-500 m/z, 5.2 scans/s) for
  the semi-quantitative survey of all 73 volatiles; selected-ion monitoring (GC-SIM) for the
  accurate quantification of the 26 odorants.** EI 70 eV, source 230 C, quadrupole 150 C.
- **GC-O.** Agilent 7890B + Gerstel ODP-3; same oven program; He 2 mL/min; sniff port 200 C; moist
  air 45 mL/min. **Four assessors**, each with >300 h of experience, each sniffing each sample
  twice = **eight sniffing runs**; a compound counts as a potent odorant at **frequency of
  detection >= 6 of 8**. Intensities scored 0-4.
- **Authentic standards and response factors — yes, and matrix-matched.** Analytical standards of
  22 named compounds from Sigma-Aldrich plus four pyrazines from YuanYe Biotechnology, all GC
  grade, plus C7-C30 alkanes for the retention indices. "**Five-point standard curves** were used
  to quantitate odorants in the pea samples. Standard stock solutions were prepared by dissolving
  the corresponding standard compound in methanol and then diluted to 1:10, 1:25, 1:50, 1:75, and
  1:100 strengths. Each standard solution was added to the matrix individually to make sure that
  the natural concentration of the studied compound was within the range ... **The calibration
  curves for each individual compound were obtained by plotting the area response ratio of the
  respective standard compound and 2,4,6-trimethylpyridine against their concentration ratio.**"
  The 26 calibration equations and their R^2 (0.940 to 0.999) are printed in Table 4. **This is a
  proper internal-standard calibration against authentic standards, in the sample's own matrix.**
- **The matrix blank.** "To obtain deodorized peas, raw pea samples were freeze-dried and milled as
  described in sample preparation. Then, the pea flours were extracted stepwise with methanol,
  dichloromethane, and pentane (each 200 mL) ... the residue was dried at room temperature." A
  separate deodorized matrix was made for the roasted sample by crushing the peas immediately after
  roasting and extracting the same way; "the absence of the abovementioned compounds was confirmed
  by GC-MS analysis". The calibration model solution is **1.8 g NaCl + 1.5 g deodorized pea flour +
  5 mL water**, i.e. the same loading as the sample.
- **The semi-quantitative survey (a different, weaker method).** "all the volatile components of raw
  and roasted peas were quantitated in terms of the internal standards by a **semiquantitative
  method**" in full-scan mode. Table 2's class sums come from this route. They are printed in
  ug/kg but they are **response-factor-1 numbers against 2,4,6-trimethylpyridine, not
  standard-calibrated** — house rule: a semi-quantitative area is not a concentration (Flags 1).
- **Basis of every concentration.** ug kg^-1 **of pea flour** (raw flour or roasted flour, each on
  its own as-milled mass). Table 4 footnote f: "The concentrations are the means of three repeated
  measurements +/- standard deviation (to three significant figures)." No dry-matter correction is
  stated anywhere.
- **OAV.** "OAV = C/T where C is the concentration of the volatile compound in the sample and T is
  the odor threshold of this compound, which was obtained from information available in the
  literature." Table 4 footnote c: "Odor threshold (OT) (ug kg^-1) **in water**". So the OAVs mix a
  flour concentration with a water threshold; they are potency indices, not activities in the
  matrix.
- **Sensory panel.** Ten trained non-smoking judges (4 M, 6 F, 23-38) from China Agricultural
  University; 10 g of flour in a covered 50 mL cup at room temperature; seven attributes agreed by
  discussion (grass-like, beans-like, fatty, popcorn-like, nutty, potato-like, smoky); 0-10 scale.
- **Recombination.** Two models: deodorized pea flour "with their natural content of water" plus
  the **9** (raw) and **20** (roasted) odorants with OAV >= 1, each at its original concentration;
  shaken 2 h; then rated on the same seven attributes.
- **Omission.** Triangle tests against the complete recombinate; **14 omission models for raw peas
  and 29 for roasted peas** (the model codes are listed verbatim in section 3 below).
- **Statistics.** Triplicate; ANOVA in SPSS v.18.0 with Tukey HSD post hoc; letters in the tables
  mark p < 0.05.

## 3. Tables re-typed

### Table 1. "Summary of Volatile Compounds Identified by Different Techniques Present in Pea Products"

A literature summary, not this paper's measurement — counts of compounds, no concentrations.
Footnote a: "Data not shown." Footnote b: "The numbers are consistent with those in references."

| no. | techniques | pea products | alcohols | aldehydes | ketones | esters | hydrocarbons | acids | benzene derivatives | phenols | others | total | ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | vacuum distillation; GC-MS | minced peas | 22 | a | | | | | | | | | 5 |
| 2 | vacuum sublimation; GC | unblanched green peas | 19 | | | | | | 4 | | | | 7 |
| 3 | vacuum sublimation; GC-MS; GC-O | unblanched green peas; pea shells | 29 | 17 | 14 | 12 | 16 | 0 | 5 | 0 | 9 | 102 | 6 |
| 4 | DHS-GC-MS; GC-O | blanched green peas | 12 | 9 | 3 | 5 | 8 | 0 | 0 | 0 | 10 | 47 | 38 |
| 5 | HS-SPME-GC-MS | pea flour | 8 | 7 | 2 | 3 | 7 | 0 | 6 | 0 | 5 | 38 | 9 |
| 6 | HS-SPME-GC-MS | pea flour | 7 | 4 | 2 | 1 | 8 | 0 | 6 | 0 | 5 | 32 | 11 |
| 7 | HS-SPME-GC-MS; SAFE-GC-MS; purge and trap-GC-MS; D-GC-O | pea flour | 20 | 9 | 21 | 3 | 16 | 2 | 6 | 1 | 9 | 87 | 8 |
| 8 | HS-SPME-GC-MS | raw/cooked peas | | | | | | | | | | | 10 |
| 9 | SAFE-GC-MS; GC-O | **pea/protein flour** | 19 | 5 | 16 | 4 | 4 | 7 | 3 | 0 | 1 | 59 | 25 |
| 10 | HS-SPME-GC-MS | cooked yellow peas | 11 | 4 | 2 | 1 | 10 | 0 | 4 | 0 | 4 | 36 | 12 |
| 11 | HS-SPME-GC-MS | cooked pea pastes | 9 | 3 | 3 | 4 | 6 | 5 | 1 | 1 | 4 | 36 | 35 |
| 12 | HS-SPME-GC-MS; GC-O | germinated pea flour | 22 | 21 | 13 | 6 | 10 | 5 | 17 | 1 | 13 | 108 | 39 |

Row 9's reference 25 is Murat et al. 2013, *Food Res. Int.* 53:31-41, "Characterisation of odour
active compounds along extraction process **from pea flour to pea protein extract**" — the only
entry in this table that reaches a protein extract, and it is **not on disk**. Flags 8.

### Table 2. "Quantitative Comparison of Different Chemical Classes in Raw and Roasted Peas"

Header as printed: `concentration (ug kg^-1)` and `percentage`. Footnote a: "Mean values in the same
row with different letters indicate that they are significantly different at p < 0.05; ND, not
detected." **These sums come from the full-scan semi-quantitative route** (section 2), not from the
calibrated Table 4 — see Flags 1.

| class | raw | roasted | raw (%) | roasted (%) |
|---|---|---|---|---|
| aldehydes | 2230 +/- 321 a | 2270 +/- 204 a | 17.08 | 1.33 |
| alcohols | 4480 +/- 727 a | 3500 +/- 484 a | 33.49 | 2.04 |
| pyrazines | ND | 39,700 +/- 2050 | ND | 23.15 |
| furanones | ND | 2780 +/- 125 | ND | 1.62 |
| pyranones | 45.2 +/- 7.96 b | 113,000 +/- 4040 a | 0.35 | 65.91 |
| other ketones | 232 +/- 42.1 b | 2940 +/- 145 a | 1.77 | 0.09 |
| esters | 84.8 +/- 18.4 a | 100 +/- 7.05 a | 0.65 | 0.06 |
| acids | 4810 +/- 520 a | 2450 +/- 43.1 b | 36.83 | 1.43 |
| hydrocarbons | 86.1 +/- 18.2 a | 73.9 +/- 8.04 a | 0.66 | 0.04 |
| phenols | 24.2 +/- 4.21 b | 5830 +/- 185 a | 0.19 | 3.40 |
| benzene derivatives | 1000 +/- 110 a | 1086 +/- 142 a | 7.68 | 0.63 |
| others | 170 +/- 32.3 b | 508 +/- 55.5 a | 1.30 | 0.30 |

The percentage columns do not sum to 100 in either column as printed (raw 100.00, roasted 98.57 —
mine); the raw column's own total, 13,162 ug/kg, is consistent with its percentages, so the roasted
percentages are the ones that do not close. Flags 4.

### Table 3. "Odorants Identified in Raw and Roasted Peas by GC-O and Means of AIs Calculation"

Footnotes as printed: a "Odor characteristics perceived by GC-O analysis"; b "MS: identification
based on the NIST 14 mass spectral database; LRIs: linear retention indices; odor: odor
descriptions; Std: confirmed by comparison with authentic standards"; c "Aroma intensities (AIs):
1 = low intensity; 2 = moderate; 3 = high; 4 = very high. Mean values in the same row with different
letters are significantly different at p < 0.05"; d "Retention indices on the HP-WAX column";
e "Retention indices on the DB-5MS column"; f "Not perceived".

| no. | compound | odorant description | identification | LRI HP-wax | LRI DB-5 | AI raw | AI roasted |
|---|---|---|---|---|---|---|---|
| 1 | dimethyl sulfide | cabbage, sulfur, sickly | MS/RI/odor/Std | 716 | 505 | f | 2.50 |
| 2 | 3-methylbutanal | malty | MS/RI/odor/Std | 914 | 649 | | 3.63 |
| 3 | **hexanal** | grass-like, green | MS/RI/odor/Std | 1066 | 803 | **3.00 a** | **2.50 b** |
| 4 | 1-pentanol | iodoform, rubber, phenolic | MS/RI/odor/Std | 1252 | 766 | 1.88 b | 2.50 a |
| 5 | 1-octen-3-one | mushroom-like | RI/odor/Std | 1299 | 975 | | 1.75 |
| 6 | (Z)-2-penten-1-ol | sour, pungent | MS/RI/odor/Std | 1318 | 891 | 1.50 | |
| 7 | 2,5-dimethylpyrazine | roasted, beans-like | MS/RI/odor/Std | 1323 | 917 | | 1.75 |
| 8 | 2,6-dimethylpyrazine | nutty | MS/RI/odor/Std | 1328 | 920 | | 2.50 |
| 9 | 2-ethylpyrazine | popcorn, beans-like | MS/RI/odor/Std | 1335 | 931 | | 2.63 |
| 10 | 2,3-dimethylpyrazine | malty | MS/RI/odor/Std | 1347 | 934 | | 2.75 |
| 11 | 1-hexanol | vegetable, herbaceous, beans-like | MS/RI/odor/Std | 1355 | 873 | 1.50 | |
| 12 | 2-ethyl-6-methylpyrazine | earthy | MS/RI/odor/Std | 1386 | 999 | | 2.50 |
| 13 | 2-ethyl-5-methylpyrazine | garlic | MS/RI/odor/Std | 1392 | 1002 | | 3.50 |
| 14 | nonanal | fatty, soapy | MS/RI/odor/Std | 1397 | 1106 | 1.75 | |
| 15 | 2,3,5-trimethylpyrazine | earthy, chocolate | MS/RI/odor/Std | 1403 | 990 | | 2.88 |
| 16 | (E)-2-octenal | fatty, green | MS/RI/odor/Std | 1433 | 1049 | 1.88 | |
| 17 | 2-ethyl-3,5-dimethylpyrazine | roasted, chocolate, cacao-like, nutty | MS/RI/odor/Std | 1445 | 1077 | | **4.00** |
| 18 | 3-ethyl-2,5-dimethylpyrazine | roasted | MS/RI/odor/Std | 1461 | — | | 3.50 |
| 19 | furfural | almond, roasted nut | MS/RI/odor/Std | 1472 | 836 | | 2.88 |
| 20 | benzaldehyde | almond, burnt sugar | MS/RI/odor/Std | 1530 | 969 | 2.75 a | 2.25 a |
| 21 | (E,Z)-2,6-nonadienal | cucumber-like | RI/odor/Std | 1588 | 1155 | 2.25 a | 2.50 a |
| 22 | benzeneacetaldehyde | flora, honey-like | MS/RI/odor/Std | 1651 | 1045 | | 3.25 |
| 23 | 3-methylbutanoic acid | sweaty | MS/RI/odor/Std | 1678 | 866 | **3.75 a** | 3.25 b |
| 24 | ethyl maltol | fruity | MS/RI/odor/Std | 1828 | — | | 3.50 |
| 25 | unknown | medicine-like | odor | 1843 | — | | 2.00 |
| 26 | 2-methoxyphenol | gammon-like, smoky | MS/RI/odor/Std | 1867 | 1092 | | 3.00 |
| 27 | benzyl alcohol | floral | MS/RI/odor/Std | 1882 | 1047 | 1.50 a | 1.88 a |
| 28 | maltol | caramel-like | MS/RI/odor/Std | 1968 | 1087 | | 3.50 |
| 29 | furaneol | caramel-like, maltol-like | MS/RI/odor/Std | 2050 | 1028 | | 2.75 |
| 30 | 2-methoxy-4-vinylphenol | smoky, clove-like | MS/RI/odor/Std | 2204 | 1194 | | 3.88 |

Ten odorants in raw peas (five aldehydes, four alcohols, one acid); 26 in roasted peas; six common
to both (hexanal, 1-pentanol, benzaldehyde, 3-methylbutanoic acid, benzyl alcohol,
(E,Z)-2,6-nonadienal).

### Table 4. "Odor Thresholds, Ions Used for Quantitation, Calibration Equations, Coefficients of Determination (R2), Concentrations, and OAVs of Odorants in Raw and Roasted Peas"

Rows are printed in descending order of the roasted-pea OAV, and that order is kept here.
Footnotes as printed: a "The numbers assigned to the compounds are consistent with those in
Table 3"; b "The compounds are listed according to the order of their OAVs in roasted pea samples";
c "**Odor threshold (OT) (ug kg^-1) in water** obtained from the literature"; d "Selected ions (m/z)
used in quantitative analysis"; e "Variables: x is the peak area relative to that of the internal
standard, 2,4,6-trimethylpyridine, and y is the concentration (ug kg^-1) in the pea sample relative
to that of the internal standard, 2,4,6-trimethylpyridine"; f "The concentrations are the means of
three repeated measurements +/- standard deviation (to three significant figures)"; g "Mean values
in the same row with different letters are significantly different at a level of p < 0.05"; h "OAV
(ratio of concentration to odor threshold)"; i "_: not detected"; j threshold from Jelen et al.
(ref 40); k threshold from van Gemert (ref 41); l threshold from Preininger et al. (ref 42).

| no. | compound | OT (ug/kg, water) | quant. ions | calibration equation | R2 | conc. raw (ug/kg) | conc. roasted (ug/kg) | OAV raw | OAV roasted |
|---|---|---|---|---|---|---|---|---|---|
| 17 | 2-ethyl-3,5-dimethylpyrazine | 0.16 j | 135, 136 | y = 0.001041x + 0.145414 | 0.992 | _ | 1250 +/- 35.1 | _ | 7822 |
| 18 | 3-ethyl-2,5-dimethylpyrazine | 8.6 k | 135, 136 | y = 0.000900x - 0.181500 | 0.997 | _ | 8980 +/- 651 | _ | 1044 |
| 24 | ethyl maltol | 100 l | 140, 139 | y = 0.000004x - 0.009267 | 0.985 | _ | 9930 +/- 848 | _ | 993 (Flags 3) |
| 30 | 2-methoxy-4-vinylphenol | 100 k | 135, 150 | y = 0.000040x - 0.124900 | 0.977 | _ | 5560 +/- 170 | _ | 556 (Flags 3) |
| 28 | maltol | 210 k | 126, 71 | y = 0.040200x - 0.167200 | 0.955 | _ | 103,000 +/- 3190 | _ | 490 |
| 2 | 3-methylbutanal | 2 k | 58, 44 | y = 0.001093x - 0.048500 | 0.940 | _ | 934 +/- 136 | _ | 467 |
| 23 | 3-methylbutanoic acid | 12 k | 60, 41 | y = 0.000005x - 0.033100 | 0.943 | **4580 +/- 471 a** | 2350 +/- 38.6 b | **382** | 196 |
| 13 | 2-ethyl-5-methylpyrazine | 16 k | 121, 122 | y = 0.000396x - 0.001607 | 0.998 | _ | 2760 +/- 82.5 | _ | 173 |
| 22 | benzeneacetaldehyde | 4 k | 91, 92 | y = 0.000413x - 0.002624 | 0.983 | _ | 421 +/- 7.62 | _ | 105 |
| 26 | 2-methoxyphenol | 3 k | 109, 124 | y = 0.000079x + 0.001071 | 0.995 | _ | 244 +/- 12.8 | _ | 81 |
| 3 | **hexanal** | **4.5 k** | 56, 44 | y = 0.003300x + 1.051800 | 0.997 | **1260 +/- 114 a** | **324 +/- 34.3 b** | **280** | **72** |
| 15 | 2,3,5-trimethylpyrazine | 23 k | 42, 122 | y = 0.000370x + 0.010085 | 0.984 | _ | 1640 +/- 35.5 | _ | 71 |
| 29 | furaneol | 60 k | 43, 57 | y = 0.000038x - 0.000325 | 0.999 | _ | 2780 +/- 125 | _ | 46 |
| 8 | 2,6-dimethylpyrazine | 1500 k | 108, 42 | y = 0.000272x + 0.008283 | 0.985 | _ | 16,100 +/- 953 | _ | 11 |
| 1 | dimethyl sulfide | 30 k | 62, 47 | y = 0.000915x + 0.007997 | 0.997 | _ | 212 +/- 11.1 | _ | 7 |
| 12 | 2-ethyl-6-methylpyrazine | 40 k | 121, 122 | y = 0.000593x + 0.031248 | 0.994 | _ | 243 +/- 7.48 | _ | 6 |
| 16 | (E)-2-octenal | 3 k | 70, 41 | y = 0.008600x + 0.048900 | 0.981 | **22.9 +/- 3.54** | _ | 8 | _ |
| 4 | 1-pentanol | 120 k | 55, 42 | y = 0.002400x - 0.026500 | 0.989 | **221 +/- 14.4 b** | 392 +/- 9.43 a | 2 | 3 |
| 27 | benzyl alcohol | 1000 k | 108, 79 | y = 0.000080x - 0.001300 | 0.998 | **1210 +/- 162 b** | 2590 +/- 382 a | 1 | 3 |
| 6 | (Z)-2-penten-1-ol | 720 k | 57, 27 | y = 0.000700x - 0.006500 | 0.986 | **901 +/- 119** | _ | 1 | _ |
| 7 | 2,5-dimethylpyrazine | 2600 k | 42, 108 | y = 0.000310x - 0.004711 | 0.985 | _ | 5960 +/- 77.6 | _ | 2 |
| 14 | nonanal | 40 k | 57, 41 | y = 0.008600x + 0.037900 | 0.993 | **69.8 +/- 7.39** | _ | 2 | _ |
| 11 | 1-hexanol | 500 k | 56, 43 | y = 0.005200x + 0.282600 | 0.988 | **595 +/- 58.5** | _ | 1 | _ |
| 19 | furfural | 282 k | 96, 95 | y = 0.000360x + 0.004434 | 0.991 | _ | 327 +/- 6.82 | _ | 1 |
| 20 | benzaldehyde | 24 k | 106, 77 | y = 0.001387x + 0.003463 | 0.999 | **54.4 +/- 7.20 a** | 15.1 +/- 0.32 b | 2 | <1 |
| 9 | 2-ethylpyrazine | 4000 k | 107, 108 | y = 0.002977x + 0.047103 | 0.995 | _ | 1810 +/- 98.8 | _ | <1 |
| 10 | 2,3-dimethylpyrazine | 400 k | 67, 108 | y = 0.000222x + 0.000797 | 0.999 | _ | 347.2 +/- 60.8 | _ | <1 |

Compounds 5, 21 and 25 (1-octen-3-one, (E,Z)-2,6-nonadienal, unknown) were sniffed but **not
quantified**: "except for three compounds that cannot be identified by MS, compound 5, 21, and 25".

### Table 5. "Results of Omission Experiments Performed on Aroma Reconstitutes of Peas"

Footnote a: "NS, no significant difference; _, not contained in this model. ***, 0.1 % significance
level. **, 1 % significance level. *, 5 % significance level."

| no. | compound(s) omitted | model I (raw peas) | model II (roasted peas) |
|---|---|---|---|
| 1 | all sulfur compound | _ | NS |
| 1-1 | dimethyl sulfide | _ | NS |
| 2 | all saturated aldehydes | *** | ** |
| 2-1 | 3-methylbutanal | _ | ** |
| 2-2 | **hexanal** | **\*\*\*** | **NS** |
| 2-3 | nonanal | NS | _ |
| 3 | all enal | *** | _ |
| 3-1 | (E)-2-octenal | *** | _ |
| 4 | all aromatic aldehydes | ** | ** |
| 4-1 | furfural | _ | NS |
| 4-2 | benzaldehyde | ** | _ |
| 4-3 | benzeneacetaldehyde | _ | ** |
| 5 | all furanone | _ | NS |
| 5-1 | furaneol | _ | NS |
| 6 | all alcohols | ** | ** |
| 6-1 | 1-pentanol | NS | * |
| 6-2 | (Z)-2-penten-1-ol | ** | _ |
| 6-3 | 1-hexanol | NS | _ |
| 6-4 | benzyl alcohol | * | ** |
| 7 | all pyrazines | _ | *** |
| 7-1 | 2,5-dimethylpyrazine | _ | NS |
| 7-2 | 2,6-dimethylpyrazine | _ | *** |
| 7-3 | 2-ethyl-6-methylpyrazine | _ | *** |
| 7-4 | 2-ethyl-5-methylpyrazine | _ | *** |
| 7-5 | 2,3,5-trimethylpyrazine | _ | *** |
| 7-6 | 2-ethyl-3,5-dimethylpyrazine | _ | *** |
| 7-7 | 3-ethyl-2,5-dimethylpyrazine | _ | ** |
| 8 | all acid | *** | ** |
| 8-1 | 3-methylbutanoic acid | *** | ** |
| 9 | all pyranones | _ | *** |
| 9-1 | ethyl maltol | _ | *** |
| 9-2 | maltol | _ | *** |
| 10 | all phenols | (blank as printed) | ** |
| 10-1 | 2-methoxy-4-vinylphenol | _ | * |
| 10-2 | 2-methoxyphenol | _ | * |

The raw-peas cell of row 10 ("all phenols") is **blank in the printed table** — neither a
significance mark nor the "_" used elsewhere. The running text lists the raw-pea omission models as
"2; 2-2; 2-3; 3; 3-1; 4; 4-2; 6; 6-1; 6-2; 6-3; 6-4; 8; and 8-1" (14 models), which does **not**
include row 10, so the blank should read "_". Recorded, not guessed.

### Numbers printed only in the running text

| quantity | value | where |
|---|---|---|
| total volatiles found | 73 compounds, >10 chemical classes | Results, first paragraph |
| odorants found by GC-O | 30 total; 10 in raw, 26 in roasted, 6 common | Results |
| odorants with OAV >= 1 | 9 in raw, 20 in roasted | Abstract and Results |
| aroma compounds significant in omission | 6 for raw peas, 15 for roasted peas (p < 0.05) | Abstract |
| fall in 3-methylbutanoic acid on roasting | "decreased by 48.65 %" | Results |
| OAV of 3-methylbutanoic acid | text says "decreased from **282** to 196"; **Table 4 and the abstract both say 382** | Results vs Table 4 (Flags 3) |
| raw-pea aroma leaders | 3-methylbutanoic acid OAV 382 and hexanal OAV 280 | Abstract |
| roasted-pea aroma leader | 2-ethyl-3,5-dimethylpyrazine, OAV 7822 | Results |
| recombination outcome | raw and roasted recombinates "similar" to the originals; the **fatty** attribute is noticeably higher in the raw recombinate than in raw pea flour | Results, Fig. 2a |

**Figure-only in this paper:** Figure 1 (heat map of all 73 volatiles, raw vs roasted) and Figure 2
(sensory profiles and recombination models, seven attributes, raw and roasted). Per house rule no
number is read off them.

### Arithmetic on the printed numbers (all mine)

1. **Hexanal survives the roast at 25.7 %** (324/1260). Equivalently the roast removes 74.3 % of
   the carried hexanal in 30 min at 160 C. On a first-order reading that is a pseudo-rate of
   ln(1260/324)/30 min = **0.0453 min^-1**, half-life **15.3 min** — but the seed is roasted whole
   and unmilled, moisture is lost and not measured, and hexanal is simultaneously being made and
   lost, so this is an apparent net disappearance, not a rate constant (see Flags 2 and 6).
2. **Benzaldehyde falls to 27.8 %** (15.1/54.4) and **3-methylbutanoic acid to 51.3 %**
   (2350/4580, i.e. a 48.7 % fall, matching the paper's "48.65 %").
3. **1-Pentanol rises 1.77x** (392/221) and **benzyl alcohol 2.14x** (2590/1210) across the roast:
   two carried volatiles that go **up**, which the paper reads as either survival or formation.
4. **Nonanal, (E)-2-octenal, (Z)-2-penten-1-ol and 1-hexanol are "not detected" in roasted flour**
   — all four are lipid-derived C6-C9 species and all four disappear below the SIM detection limit.
   The paper prints no detection limit, so "_" is an absence, not a zero (Flags 7).
5. **Total quantified carried (lipid-derived + pre-existing) volatiles in raw flour**, summing the
   nine raw entries of Table 4: 4580 + 1260 + 1210 + 901 + 595 + 221 + 69.8 + 54.4 + 22.9 =
   **8914 ug/kg**. Against Table 2's raw total of 13,162 ug/kg (mine, summing the class column),
   the calibrated set accounts for **68 %** of the semi-quantitative total.
6. **The pyranone class sum reconciles exactly with Table 4**: maltol 103,000 + ethyl maltol
   9930 = 112,930, printed as 113,000 in Table 2. The furanone class sum equals furaneol exactly
   (2780). The **aldehyde** class does not: Table 4's raw aldehydes sum to 1407 ug/kg against
   Table 2's 2230, the difference being aldehydes seen in the survey but never calibrated. So
   Table 2 mixes calibrated and semi-quantitative rows and cannot be treated uniformly (Flags 1).
7. **Every printed OAV reproduces as concentration / threshold to the rounding shown, except two.**
   I checked all 27 rows. Ethyl maltol: 9930/100 = 99.3, printed **993**. 2-Methoxy-4-vinylphenol:
   5560/100 = 55.6, printed **556**. Both are exactly 10x. The threshold column was re-read from
   the page image at 200 dpi and both thresholds are unambiguously **100**, so the error is in the
   OAVs, not in my reading of the text layer (Flags 3).
8. **Pyrazine total in roasted flour**, summing the nine Table 4 pyrazine rows: 1250 + 8980 + 2760
   + 1640 + 16,100 + 243 + 5960 + 1810 + 347.2 = **39,090 ug/kg**, against Table 2's 39,700 — a
   1.5 % difference, so the pyrazine class is essentially fully calibrated.

## 4. Numbers the repository can use

Every row shares the same material and cook unless stated: *Pisum sativum* L. cv. Zhongwan, Datong
(Shanxi) 2017 harvest, stored 4 C; **raw** = whole dry seed milled to flour through 0.30 mm;
**roasted** = whole seed 160 C / 30 min forced-air oven, then milled the same way; HS-SPME (1.5 g
flour + 1.8 g NaCl + 5 mL water, 50 C, 20 min equilibration + 50 min extraction) with GC-SIM
quantification against five-point matrix-matched calibration curves and 2,4,6-trimethylpyridine as
internal standard; triplicate; basis **ug kg^-1 of flour as milled**, moisture not stated.

| compound | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| hexanal | 1260 +/- 114 | ug/kg of raw pea flour | whole pea seed milled to flour, unheated | Table 4 p. 2723 | measured_level |
| hexanal | 324 +/- 34.3 | ug/kg of roasted pea flour | whole seed 160 C / 30 min, then milled | Table 4 | measured_level |
| hexanal | 0.257 | ratio roasted/raw, dimensionless | same cook | derived from Table 4 (mine) | within_study_ratio |
| nonanal | 69.8 +/- 7.39 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| nonanal | not detected | — | roasted flour | Table 4 | level_only (an absence, no LOD printed) |
| (E)-2-octenal | 22.9 +/- 3.54 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| (E)-2-octenal | not detected | — | roasted flour | Table 4 | level_only |
| 1-hexanol | 595 +/- 58.5 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| 1-hexanol | not detected | — | roasted flour | Table 4 | level_only |
| 1-pentanol | 221 +/- 14.4 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| 1-pentanol | 392 +/- 9.43 | ug/kg of roasted pea flour | 160 C / 30 min | Table 4 | measured_level |
| (Z)-2-penten-1-ol | 901 +/- 119 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| benzyl alcohol | 1210 +/- 162 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| benzyl alcohol | 2590 +/- 382 | ug/kg of roasted pea flour | 160 C / 30 min | Table 4 | measured_level |
| benzaldehyde | 54.4 +/- 7.20 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| benzaldehyde | 15.1 +/- 0.32 | ug/kg of roasted pea flour | 160 C / 30 min | Table 4 | measured_level |
| 3-methylbutanoic acid | 4580 +/- 471 | ug/kg of raw pea flour | as above | Table 4 | measured_level |
| 3-methylbutanoic acid | 2350 +/- 38.6 | ug/kg of roasted pea flour | 160 C / 30 min | Table 4 | measured_level |
| 3-methylbutanal | 934 +/- 136 | ug/kg of roasted pea flour | 160 C / 30 min; not detected raw | Table 4 | measured_level |
| benzeneacetaldehyde (phenylacetaldehyde) | 421 +/- 7.62 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| furfural | 327 +/- 6.82 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| furaneol (HDMF) | 2780 +/- 125 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| maltol | 103,000 +/- 3190 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| ethyl maltol | 9930 +/- 848 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-methoxyphenol (guaiacol) | 244 +/- 12.8 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-methoxy-4-vinylphenol (4-vinylguaiacol) | 5560 +/- 170 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| dimethyl sulfide | 212 +/- 11.1 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2,5-dimethylpyrazine | 5960 +/- 77.6 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2,6-dimethylpyrazine | 16,100 +/- 953 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2,3-dimethylpyrazine | 347.2 +/- 60.8 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-ethylpyrazine | 1810 +/- 98.8 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-ethyl-6-methylpyrazine | 243 +/- 7.48 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-ethyl-5-methylpyrazine | 2760 +/- 82.5 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2,3,5-trimethylpyrazine | 1640 +/- 35.5 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 2-ethyl-3,5-dimethylpyrazine | 1250 +/- 35.1 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| 3-ethyl-2,5-dimethylpyrazine | 8980 +/- 651 | ug/kg of roasted pea flour | as above | Table 4 | measured_level |
| total pyrazines | 39,700 +/- 2050 (survey) / 39,090 (sum of Table 4, mine) | ug/kg of roasted pea flour | as above; ND raw | Table 2 / derived | peak_area_only (Table 2) and derived (the sum) |
| class sums, all eleven classes, raw and roasted | see Table 2 above | ug/kg, semi-quantitative | as above | Table 2 p. 2721 | **peak_area_only** — response factor 1 against the internal standard, full scan, no authentic standard |
| odour threshold, hexanal | 4.5 | ug/kg in water | — | Table 4 footnote k, from van Gemert 2011 | threshold |
| odour threshold, nonanal | 40 | ug/kg in water | — | Table 4 footnote k | threshold |
| odour threshold, 1-hexanol | 500 | ug/kg in water | — | Table 4 footnote k | threshold |
| odour threshold, furfural | 282 | ug/kg in water | — | Table 4 footnote k | threshold |
| odour thresholds, the other 22 odorants | see the OT column of Table 4 above | ug/kg in water | — | Table 4, sources k / j / l per row | threshold |
| OAV, hexanal | 280 (raw), 72 (roasted) | dimensionless | flour concentration over water threshold | Table 4 | within_study_ratio |
| OAV, all 27 quantified odorants | see Table 4 above | dimensionless | as above | Table 4 | within_study_ratio (two rows do not reconcile — Flags 3) |
| aroma intensities, 30 odorants, raw and roasted | see Table 3 above | 0-4 sniff score, 4 assessors x 2 runs | GC-O, detection frequency >= 6/8 | Table 3 p. 2722 | level_only (a panel score, not a concentration) |
| omission significance, 34 models | see Table 5 above | p-level marks | triangle tests against the full recombinate | Table 5 p. 2724 | level_only |
| sensory profiles and recombination fit | — | seven attributes, 0-10 | raw and roasted, 10 judges | Fig. 2a, 2b | **figure_only** |
| heat map of all 73 volatiles | — | — | raw vs roasted | Fig. 1 | **figure_only** |
| apparent net hexanal disappearance | 0.0453 min^-1; half-life 15.3 min | min^-1, first-order reading | 160 C, 30 min, whole seed | derived from Table 4 (mine) | derived_assumption — **not a rate constant**; two points only, whole-seed geometry, no moisture correction, formation and loss confounded |

### What this can and cannot be used for

**Can:** supply the raw-flour side of Programme 7's levels table for hexanal, nonanal,
(E)-2-octenal, 1-hexanol, 1-pentanol, (Z)-2-penten-1-ol, benzyl alcohol, benzaldehyde and
3-methylbutanoic acid, in a **flour**, with a standard deviation from triplicates that can serve as
the band. Supply a measured survival fraction across one real dry cook. Supply a cited odour
threshold for six compounds `off_flavour_targets.yml` currently carries uncited.

**Cannot:** stand in for an isolate. The roadmap's sentence is about "a pea or soy isolate", and
the roadmap's own phrasing warns that an isolate is not a flour. Between flour and isolate lie
protein extraction, an alkaline or salt wash and a spray drying, each of which changes the volatile
load; the only paper in this article's own Table 1 that follows that path (Murat 2013, "from pea
flour to pea protein extract") is not on disk. Cannot supply a rate: two points, one temperature.
Cannot supply a benchmark row directly without a stated moisture, because the flour concentration
basis is as-milled mass.

## 5. Flags

1. **Table 2 is semi-quantitative and Table 4 is not, and the paper does not separate them
   typographically.** Both print `ug kg^-1`. Table 4's numbers come from five-point calibration
   curves against authentic standards in a deodorized pea-flour matrix and are concentrations.
   Table 2's class sums come from the full-scan survey "quantitated in terms of the internal
   standards by a semiquantitative method", i.e. every one of the 73 compounds is scaled by the
   response of 2,4,6-trimethylpyridine as if its response factor were 1. **Only Table 4's 27 rows
   may enter the repository as concentrations.** Table 2 is `peak_area_only` despite its unit. My
   own reconciliation (section 3, arithmetic 6 and 8) shows the pyranone, furanone and pyrazine
   class sums do match the calibrated rows, so those three classes are probably built from Table 4;
   the aldehyde class demonstrably is not.
2. **No moisture, no dry basis, and a cook that removes water.** The paper prints no moisture
   content for either flour and no mass loss across the roast. A 160 C / 30 min forced-air roast of
   whole seed removes water, so the roasted flour's kilogram is not the raw flour's kilogram. Every
   raw-to-roasted ratio in section 3 is therefore **an as-is ratio, not a dry-basis ratio**, and
   the true dry-basis fall in hexanal is larger than 74 % by whatever the moisture loss was. This
   is the single largest quantitative uncertainty in using the pair.
3. **Three internal inconsistencies in the printed numbers.** (i) The running text says the OAV of
   3-methylbutanoic acid "decreased from **282** to 196"; Table 4 and the abstract both say **382**,
   and 4580/12 = 381.7, so 282 is a typographical error and 382 is right. (ii) The OAV printed for
   ethyl maltol is 993 where 9930/100 = 99.3. (iii) The OAV printed for 2-methoxy-4-vinylphenol is
   556 where 5560/100 = 55.6. Both errors are exactly 10x; both thresholds were confirmed as 100
   from the rendered page image, so the OAVs are wrong, not the thresholds. All 24 other OAV rows
   reconcile. **Use the concentrations and recompute the OAVs; do not import the printed OAVs
   for these two compounds.**
4. **Table 2's roasted percentage column sums to 98.57 %, not 100** (mine). The raw column closes.
   Nothing depends on this, but it means the roasted percentages carry a rounding or a missing
   class.
5. **The concentrations are headspace-derived from a rehydrated slurry at 50 C, not from the dry
   flour.** 1.5 g of flour is suspended in 5 mL of saturated brine and held at 50 C for 70 min
   before the fibre is withdrawn. The calibration is matrix-matched (deodorized pea flour at the
   same loading), which corrects for the *average* binding but not for any difference in binding
   between raw and roasted flour — and roasting denatures protein and changes exactly that. **The
   raw-vs-roasted comparison assumes the matrix effect is unchanged by the cook; the paper does not
   test this.** For Programme 7, which explicitly wants "the matrix binding applied", note that the
   binding is already folded into these numbers once.
6. **Nothing here is resolved in time.** There is one roast: 160 C, 30 min. Raw versus roasted is a
   two-point contrast, not a series. No second temperature, no second duration, no intermediate
   sample, no come-up profile. Any rate read from this pair is a two-point apparent rate and is
   marked as such.
7. **"Not detected" has no limit attached.** Four lipid-derived species (nonanal, (E)-2-octenal,
   1-hexanol, (Z)-2-penten-1-ol) are marked "_" in the roasted column and pyrazines and furanones
   are "ND" in the raw column of Table 2. No detection limit or limit of quantification is printed
   anywhere in the paper. Treat these as "below an unstated limit", never as zero.
8. **What this paper does NOT contain**, and what to request: no protein isolate and no protein
   concentrate (the material is whole-seed flour throughout); no lipid, protein, moisture, free
   amino acid or sugar composition of either flour; no lipoxygenase or other enzyme activity; no
   2-pentylfuran and no methoxypyrazine — **two of the four volatiles the roadmap names are simply
   absent from this paper's 30 odorants**; no storage series; no second cultivar; no cultivar-level
   variation; no rate constant; no water activity. **To request:** (a) the numeric matrix behind
   Figure 1 (73 compounds x 2 samples), which would triple the compound coverage; (b) the moisture
   content of the two flours, without which the ratios cannot be put on a dry basis; (c) the
   Supporting Information (SPME optimisation), which is not on disk; (d) Murat et al. 2013,
   *Food Res. Int.* 53:31-41, this paper's own reference 25, the flour-to-protein-extract series
   that Programme 7 actually needs and that is **not on disk**.
9. **Registry gaps against `data/keys/compounds.yml`.** Present and directly keyable: `hexanal`,
   `nonanal`, `e_2_octenal`, `1_hexanol`, `benzaldehyde`, `furfural`, `hdmf` (furaneol),
   `phenylacetaldehyde` (benzeneacetaldehyde), `3_methylbutanal`, `4_vinylguaiacol`
   (2-methoxy-4-vinylphenol), `2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`,
   `2_3_dimethylpyrazine`, `2_ethylpyrazine`, `trimethylpyrazine` (2,3,5-trimethylpyrazine),
   `2_ethyl_3_5_dimethylpyrazine`, and the group id `pyrazines`. **Absent from the registry:**
   dimethyl sulfide (the registry has only `dimethyl_disulfide` and `dimethyl_trisulfide`), maltol,
   ethyl maltol, 2-methoxyphenol (guaiacol), 3-methylbutanoic acid, 1-pentanol,
   (Z)-2-penten-1-ol, benzyl alcohol, 1-octen-3-one, (E,Z)-2,6-nonadienal,
   2-ethyl-6-methylpyrazine, 2-ethyl-5-methylpyrazine and 3-ethyl-2,5-dimethylpyrazine.
   **The last three matter most**: they are three of the four highest-OAV pyrazines in roasted pea
   and the repository cannot name them.
10. **Registry gaps against `data/species/off_flavour_targets.yml`.** That file holds exactly six
    compounds: Hexanal, Nonanal, 1-Octen-3-ol, 2-Pentylfuran, 1-Hexanol, Furfural. This paper
    prints a raw-flour level for **four** of them (hexanal, nonanal, 1-hexanol; furfural only in
    the roasted flour) and prints **neither 1-octen-3-ol nor 2-pentylfuran at all** — so for two of
    the roadmap's four named carried volatiles this paper is silent, and the methoxypyrazines the
    roadmap also names are absent both from this paper and from that file (they appear in
    `compounds.yml` only, as the group `methoxypyrazines` and the single species
    `3_isobutyl_2_methoxypyrazine`). The hexanal threshold in that file, 4.5 ug/kg, is **the same
    value this paper uses**, attributed here to van Gemert 2011 rather than to the Belitz/Czerny
    pair the YAML cites; that is an independent corroboration of the number, not a second source
    for it.
