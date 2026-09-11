# Wang 2025 — EXTRACTION (a coffee-analogue solid food matrix built on chickpea flour: high-moisture twin-screw extrusion at 100 C, drying to 11 % moisture at 45 C, then fluidised-bed roasting at 170 C for 20 min; single-amino-acid, single-sugar and multiple-precursor batches with chlorogenic acid at 0-5 %; volatiles by HS-SPME GC-MS against one internal standard)

### A dry-roast Maillard model built out of a real legume flour, with the amino acid, the sugar and the chlorogenic acid charged at declared percentages and the pyrazine-versus-furan partition swinging with each — but its volatiles are quantified against a single internal standard with no authentic standards, so the ppb are semi-quantitative, and not one lipid-oxidation volatile appears anywhere in it.

**Source on disk:** `data/articles/Wang2025.pdf` (10 pp., Food Chemistry 492 (2025) 145552,
**open access, CC BY**). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Wang2025.txt`), whole file. Tables 1, 2, 3 and 4 came through clean and are
re-typed in full below; **note that Tables 2, 3 and 4 use the comma as the decimal separator
throughout** (e.g. "1446,91 ± 18,64" means 1446.91 ± 18.64), with a single inconsistent entry that
uses the point (Flags 5). Figures 1 (the process scheme), 2 (expansion ratio and bulk density),
3 (total furans and pyrazines by formulation) and 4 (the Pearson correlation heat map) are images
and are **figure-only**. **The Supplementary Material — Table S1 (CIE-LAB colour), Table S2 (the
oil-type volatile comparison), Table S3 (the 27 volatiles and 22 OAVs of the real coffee
reference), Table S4 (21 pyrazines across all sixteen amino acids), Figure S1 (the colour PCA) and
Figure S2 (the pH of every formulation) — is online-only and is NOT on disk.** That is a serious
loss for this paper: Table S3 holds every odour threshold and every odour-activity value, and
Figure S2 holds every pH. Repo status before this dossier: no dossier, and the paper is cited
nowhere in `src/`, `data/` or `results/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Generation of Maillard reaction-derived flavour compounds in coffee analogues by an extrusion-based solid model system" |
| Authors | Bei Wang, Marjanne J. Verhoeven, Vincenzo Fogliano (corresponding, vincenzo.fogliano@wur.nl), Teresa Oliviero — **Food Quality and Design Group, Wageningen University & Research, The Netherlands** |
| Venue | **Food Chemistry 492 (2025) 145552**. Received 29 Jan 2025, revised 11 Jul 2025, accepted 12 Jul 2025, online 14 Jul 2025 |
| DOI | **10.1016/j.foodchem.2025.145552** |
| Licence | Open access, CC BY |
| Matrix | **Chickpea flour** (Smaakt brand): **21 % protein, 6.7 % fat, 58 % carbohydrates** |
| Naming | **SFM** = solid food matrix (the model system); **EC** = extrudate cylinders; **CBA** = coffee bean analogue, i.e. the roasted extrudate; **AA** = a 1:1:1:1:1:1 mixture of Ile, Leu, Val, Asp, Gly and Ser; **S** = sucrose; **CA** = chlorogenic acid; **R** = rhamnose. Numbers after a code are the w/w % in the final extrudate, so AA5/S6/CA1 = 5 % amino acid mixture, 6 % sucrose, 1 % chlorogenic acid |
| Funding | China Scholarships Council No. 202208520011 |
| Data | "Data will be made available on request." |
| Lineage | The group's own dry-model line (Capuano 2010, Oliviero 2009); the "in-bean" biomimetic coffee work of Poisson 2020 and Cerny 2021; van Boekel's solid-versus-liquid model argument |

## 1. Why it matters

For the carried-volatile question this cluster was assembled around, the answer is short: **no
lipid-oxidation volatile is measured anywhere in this paper.** Hexanal, 2-pentylfuran,
1-octen-3-ol, nonanal, 1-hexanol and the methoxypyrazines do not appear in any of its four tables.
The chickpea flour brings 6.7 % fat into the extruder and the formulations add 10 % sunflower oil
on top, and the paper still reports no aldehyde other than the three Strecker aldehydes and
2-methylbutanal. Whether that is because they were not formed, not extracted by a 20-minute
headspace at 60 C, or not looked for, the paper does not say (Flags 2). **The Programme 7 levels
table gets no row from this paper.**

What it does give the repository is different and, for other parts of the roadmap, substantial:

1. **A blank that is a legume flour and nothing else.** The Blank formulation is **67 % chickpea
   flour + 33 % water**, extruded, dried and roasted at 170 C for 20 min, with no added amino acid,
   no added sugar, no added oil and no chlorogenic acid. It still produces methylpyrazine
   (293.45 ppb), trimethylpyrazine (101.74), 3-ethyl-2,5-dimethylpyrazine (129.06), furfural
   (163.97), 5-methylfurfural (26.89) and 2-furanmethanol (41.44). **That is the Maillard output of
   a chickpea flour's own free precursors under a stated dry roast**, and it is the closest thing in
   this cluster to "what the flour brings on its own" — for the Maillard side rather than the
   lipoxygenase side.
2. **A chlorogenic acid dose-response on the pyrazine-versus-furan partition.** Chlorogenic acid at
   0, 1 and 5 % w/w, at fixed amino acid and sugar. Going from 1 % to 5 % CA, **total pyrazines fall
   from 3416.36 ± 122.40 to 1753.36 ± 70.23 ppb** while **total furans rise nearly three-fold**.
   `data/keys/compounds.yml` carries `chlorogenic_acid` as a keyed compound, and the repository has
   a standing interest in phenolic suppression of pyrazine formation; this is a clean, dosed,
   solid-phase demonstration with a stated pH mechanism.
3. **A sucrose dose-response and an amino-acid dose-response on the same axis**, at 2 % vs 6 %
   sucrose and 1 % vs 5 % amino acid mixture, in the same matrix and cook.
4. **A full extrusion + roast process specification** — screw geometry, feed rate, water rate, screw
   speed, barrel temperature, die diameter, moisture at each stage, drying and roasting conditions.
   Together with `conti2025_extraction.md` and `conti2025b_extraction.md` this makes two independent
   descriptions of an extruded plant-protein cook in the corpus.
5. **A per-amino-acid pyrazine fingerprint.** Sixteen amino acids each at 10 % in chickpea flour,
   each extruded and roasted alone; Table 2 prints the six that were carried forward. Glycine gives
   trimethylpyrazine at 3626.89 ppb and tetramethylpyrazine at 1224.94; serine gives
   2-ethyl-6-methylpyrazine at 1673.74 and 2,6-diethylpyrazine at 1512.44; valine, isoleucine and
   leucine each give their own Strecker aldehyde and nothing else's. This bears directly on the
   roadmap's "amino-acid identity" programme (section 7 of the roadmap, "Strecker aldehydes within
   threefold on another laboratory"), though the quantification is too weak to test a threefold
   claim (Flags 1).

## 2. Methods as they matter to a model

- **The matrix, exactly as described.** "The **chickpea flour (Smaakt)** contained **21 % protein,
  6.7 % fat, and 58 % carbohydrates**." It is a **flour** — the whole milled seed, not a
  concentrate and not an isolate. No moisture, ash, free amino acid, free sugar or fatty-acid
  profile is printed, and no lipoxygenase status. The reference material is "**100 % Arabica coffee
  beans (Ekoplaza, strength 8)**".
- **The oils.** Conventional sunflower oil (Reddy): **10 % saturated, 27 % monounsaturated, 55 %
  polyunsaturated**. High-oleic sunflower oil (Reddy): **8.7 % saturated, 51 % monounsaturated,
  33 % polyunsaturated**, "containing at least 80 % oleic acid" per the text. Supplier figures.
- **The emulsions.** Each oil prepared as an oil-in-water emulsion at **9.5 % and 21 % oil-in-water**
  to reach **10 % and 15 % w/w oil in the final extrudate**. Aqueous phase: **3.6 % chickpea flour**
  as the emulsifier, dissolved in deionised water at 55 C, stirred at 600 rpm for 30 min.
  Homogenised on an IKA T25 ultra-turrax at **10,000 rpm for 10 min then 20,000 rpm for 5 min**.
- **Extrusion, in full.** **Thermo Fisher Process 11 parallel twin-screw extruder**, eight heating
  blocks, **barrel L/D = 40**, **barrel diameter 11 mm**, **die opening 6 mm**. Ingredients premixed
  outside and fed by a volumetric feeder at **1.41 kg/h**; tap water pumped by a peristaltic pump at
  **623 mL/h**; **screw speed 600 rpm**; **all eight heating blocks at 100 C**. Extrudate collected
  after steady state; **moisture of the extrudate 33 ± 1 % w/w**. Cut into cylinders about **7 mm
  long, 7-10 mm in diameter** depending on expansion. **No residence time, no die pressure, no melt
  temperature and no specific mechanical energy are printed** (Flags 3).
- **Drying and roasting.** Cylinders spread on a mesh and dried in an incubator at **45 C for 10 h**
  to **11 ± 1 % w/w moisture** (checked on a Kern DAB 200-2 moisture analyser), cooled in a
  desiccator, then roasted in a **Toper Optical fluidised-bed roaster at 170 C for 20 min**.
  Vacuum-packed and stored at **-20 C**. **One roast condition only: 170 C, 20 min. There is no
  temperature series and no time series** (Flags 4).
- **The three experimental batches.**
  1. *Single amino acid*: chickpea flour + one of sixteen amino acids at **10 % of the final
     extrudate**, **no emulsion**.
  2. *Single sugar*: the same with one of five sugars (glucose, sucrose, rhamnose, xylose,
     fructose) at **10 %**, **no emulsion**.
  3. *Multiple precursor*: chickpea flour + the six-amino-acid mixture + rhamnose + sucrose +
     chlorogenic acid, with **9.5 % sunflower-oil-in-water emulsion replacing the water** so the
     final oil is 10 %. Seven formulations, in Table 1.
- **pH.** "The samples (2 g) were extracted with 15 mL of boiling water in a 20-mL tube ... placed
  on a rotator overnight (room temperature). The extracts were then filtered with No. 4 Whatman
  filter paper before pH analysis." **Every pH value is in Figure S2, which is not on disk**; the
  running text quotes only two: **5.3** for CBA/AA1/S6/CA5 (the highest-chlorogenic-acid,
  highest-sucrose formulation) and **6.3** for the formulation without chlorogenic acid.
- **Volatile extraction.** HS-SPME, **50/30 um PDMS/DVB/CAR fibre**. The roasted analogue was ground
  in a coffee grinder and sieved through a **1.25 mm** test sieve; **0.25 g of ground sample in a
  10 mL vial**; **fibre inserted for 20 min at 60 C**; desorbed immediately in the injector.
  **Internal standard: 4-methyl-1-pentanol** (the amount added and its concentration are **not
  printed** — Flags 1).
- **GC-MS.** Thermo Trace GC Ultra + DSQ II with autosampler; **Stabilwax-DA (0.25 um x 0.25 mm x
  20 m)**; 1 uL injected splitless at 250 C; oven **40 C for 5 min, 12 C/min to 225 C, hold 5 min**;
  helium > 99.99 % at **1.2 mL/min**; ion source 250 C, 70 eV, transfer line 250 C; scan **m/z
  25-300**.
- **Identification.** "The volatile compounds were identified by using the mass spectra library
  **NIST014s** and further confirmed by the **retention index (RI)** calculated based on literature
  (van Den Dool & Kratz, 1963) ... using an n-alkane mixture (C7-C40) ... The values were then
  compared with previously published reports to verify their reliability."
- **Quantification — and this is the weak point.** Verbatim, and complete: "**The concentration of
  the identified compound was carried out based on the internal standard used.**" That is the whole
  of it. **No authentic standard is injected for any analyte, no calibration curve is built, and no
  response factor is determined**; the chemicals list in section 2.2 contains sugars, amino acids,
  chlorogenic acid, the C7-C40 alkanes and 3-heptanone (listed as "analytical standard") but **not
  one of the pyrazines, furans or aldehydes that Tables 2, 3 and 4 report**. So every ppb in this
  paper is a peak area ratioed to 4-methyl-1-pentanol with an implicit response factor of 1. Under
  the house rule that is **`peak_area_only`, not a concentration** (Flags 1). The reported precision
  is a repeatability statement only: "Mean concentration of aroma compounds expressed in part per
  billion of the sample powder with **RSD (%) below 15 % (n = 3)**".
- **Odour-activity values.** "OAV was calculated by the ratio of the concentration of each compound
  to its threshold value which was obtained from references (Caporaso 2018; Feng 2015; Liu 2021;
  van Gemert 2011). Compounds with OAV >= 1 were considered as potential contributors."
  **Every OAV and every threshold is in Table S3, which is not on disk**; the running text quotes
  only the coffee reference's top few.
- **Expansion ratio and bulk density.** Expansion ratio = extrudate diameter / die diameter
  (Vernier calliper). Bulk density by polenta displacement, with the equation printed.
- **Colour.** CIE L*a*b* on an Iris Visual Analyser (Alpha M.O.S.) with a 2592 x 1944 CCD, top
  lighting, X-rite Colour Checker Passport calibration, **25 extrudate cylinders per picture**,
  L*a*b* averaged over every pixel. Values in Table S1, **not on disk**.
- **Basis of every concentration.** "**parts per billion of the sample powder**", i.e. ug/kg of the
  ground roasted analogue at ~11 % moisture as roasted, n = 3.
- **Statistics.** Triplicate; one-way ANOVA at p < 0.05 in SPSS 28.0.1.1; Pearson correlation
  coefficients across the 12 volatiles common to the multiple-precursor formulations (Fig. 4).

## 3. Tables re-typed

### Table 1. "Composition and concentration (w/w) of formula used for multiple-precursor solid food matrix model systems"

Footnotes as printed: "In the formulation code the numbers following the ingredient abbreviations
indicate the concentration by weight of the ingredients in the final extrudates. So AA5/S6/CA1 means
5 % amino acid mixture, 6 % sucrose and 1 % Chlorogenic acid." and "a AA is a mixture of six amino
acids (**Ile, Leu, Val, Asp, Gly and Ser**) all at the same concentration."

| formulation | chickpea flour | water | oil | amino acid mixture (AA) | rhamnose | sucrose (A) | chlorogenic acid (CA) |
|---|---|---|---|---|---|---|---|
| Blank | 67.0 | 33.0 | 0 | 0 | 0 | 0 | 0 |
| AA1/S2/CA1 | 51.9 | 33.1 | 10.0 | 1 | 1 | 2 | 1 |
| AA1/S2 | 52.9 | 33.1 | 10.0 | 1 | 1 | 2 | 0 |
| AA1/S2/CA5 | 47.9 | 33.1 | 10.0 | 1 | 1 | 2 | 5 |
| AA1/S6/CA5 | 43.9 | 33.1 | 10.0 | 1 | 1 | 6 | 5 |
| AA1/S6/CA1 | 47.9 | 33.1 | 10.0 | 1 | 1 | 6 | 1 |
| AA5/S6/CA1 | 43.9 | 33.1 | 10.0 | 5 | 1 | 6 | 1 |

**Every row sums to exactly 100.0 % (mine, checked all seven).** Note that rhamnose is present at
1 % in all six precursor formulations and is never varied, and that "S" in the code refers to
sucrose only.

### Table 2. "Concentrations of volatile compounds identified from coffee bean analogue (CBA) containing six selected amino acids (Asp, Gly, Ser, Ile, Leu, and Val)"

Footnote as printed: "a n.d., not detectable. The abbreviation of amino acids represents the amino
acids contained in the coffee bean analogue. Mean concentration expressed in parts per billion with
RSD (%) below 15 % (n = 3)." **Each of these is chickpea flour + that one amino acid at 10 % of the
final extrudate, no emulsion, roasted 170 C / 20 min.** Decimal commas as printed are converted to
points here; the one entry that was printed with points is marked.

| compound | CBA-Asp | CBA-Gly | CBA-Ser | CBA-Ile | CBA-Leu | CBA-Val |
|---|---|---|---|---|---|---|
| isobutyraldehyde | n.d. | n.d. | n.d. | n.d. | n.d. | **1446.91 ± 18.64** |
| 2-methylbutanal | n.d. | n.d. | n.d. | **1156.53 ± 66.29** | n.d. | n.d. |
| 3-methyl-butanal | n.d. | n.d. | n.d. | n.d. | **1231.76 ± 22.8** | n.d. |
| methylpyrazine | 288.51 ± 11.94 | 34.35 ± 3.03 | 700.82 ± 63.72 | 254.48 ± 7.66 | n.d. | n.d. |
| 2,5-dimethylpyrazine | n.d. | 92.4 ± 4.53 | n.d. | n.d. | 23.7 ± 1.32 | 211.57 ± 7.86 |
| 2,6-dimethylpyrazine | n.d. | **556.68 ± 25.02** | n.d. | n.d. | n.d. | n.d. |
| 2,3-dimethylpyrazine | n.d. | 176.05 ± 4.03 | 349.84 ± 24.73 | n.d. | n.d. | n.d. |
| 2-ethyl-6-methylpyrazine | 92.87 ± 4.99 | 89.54 ± 2.46 | **1673.74 ± 42.95** | n.d. | n.d. | n.d. |
| 2-ethyl-5-methylpyrazine | 69.02 ± 1.78 | 100.2 ± 11.08 | 418.25 ± 6.87 | n.d. | 45.14 ± 1.41 | 123.34 ± 7.83 |
| 2-ethyl-3-methylpyrazine | 51.51 ± 1.87 | n.d. | 634.37 ± 7.13 | n.d. | n.d. | n.d. |
| trimethylpyrazine | n.d. | **3626.89 ± 49.22** | n.d. | n.d. | 54.72 ± 1.28 | 221.68 ± 1.19 |
| 2,6-diethylpyrazine | n.d. | n.d. | **1512.44 ± 2.19** | n.d. | n.d. | n.d. |
| 3-ethyl-2,5-dimethylpyrazine | 683.63 ± 14.41 | **1227.84 ± 27.23** | 453.35 ± 42.47 | 89.29 ± 8.08 | n.d. | 308.04 ± 8.52 |
| 2,3-dimethyl-5-ethylpyrazine | n.d. | n.d. | 727.51 ± 84.35 | n.d. | n.d. | n.d. |
| 2-ethyl-3,5-dimethylpyrazine | 27.38 ± 2.52 | n.d. | n.d. | n.d. | n.d. | 101.72 ± 2.32 |
| tetramethylpyrazine | n.d. | **1224.94 ± 87.75** | n.d. | n.d. | n.d. | n.d. |
| 2-methyl-3,5-diethylpyrazine | 222.2 ± 3.66 | n.d. | 544.43 ± 42.36 | n.d. | n.d. | n.d. |
| 3,5-dimethyl-2-propylpyrazine | n.d. | n.d. | 73.27 ± 6.1 | n.d. | n.d. | n.d. |
| 2,3-diethyl-5-methylpyrazine | 64.15 ± 2.81 | n.d. | 328.99 ± 23.09 | n.d. | n.d. | n.d. |
| 2,3,5-trimethyl-6-ethylpyrazine | n.d. | 484.04 ± 15.21 | n.d. | n.d. | n.d. | n.d. |
| 2-isoamyl-6-methylpyrazine | n.d. | n.d. | n.d. | 524.99 ± 15.16 | n.d. | n.d. |
| 2-isobutyl-3-methylpyrazine | n.d. | n.d. | n.d. | **153.60 ± 2.67** *(printed with decimal points, unlike every other entry)* | n.d. | n.d. |

### Table 3. "Concentrations of volatile compounds identified from coffee bean analogue (CBA) containing specific sugars (Glc, Suc, Rha, Xyl and Fru)"

Footnote as printed: "a n.d., not detectable. The abbreviation of sugar represents the sugar
contained in the coffee bean analogue. Mean concentration expressed in parts per billion with RSD
(%) below 15 % (n = 3)." **Chickpea flour + one sugar at 10 % of the final extrudate, no emulsion,
roasted 170 C / 20 min.**

| compound | CBA-Glc | CBA-Suc | CBA-Rha | CBA-Xyl | CBA-Fru |
|---|---|---|---|---|---|
| 3-methylfuran | 12.07 ± 1.44 | n.d. | 10.57 ± 1.21 | n.d. | n.d. |
| 2,5-dimethylfuran | n.d. | n.d. | 64.35 ± 3.4 | n.d. | n.d. |
| 2,3-butanedione | 39.21 ± 2.53 | n.d. | n.d. | n.d. | n.d. |
| 2,3-pentanedione | 61.76 ± 5.01 | n.d. | n.d. | n.d. | n.d. |
| methylpyrazine | n.d. | 708.77 ± 11.45 | n.d. | n.d. | n.d. |
| 2,5-dimethylpyrazine | n.d. | 412.19 ± 3.55 | n.d. | n.d. | n.d. |
| 2-ethyl-6-methylpyrazine | n.d. | 126.38 ± 3.74 | n.d. | n.d. | n.d. |
| trimethylpyrazine | n.d. | 483.97 ± 14.55 | n.d. | n.d. | n.d. |
| 5-methylfuran-2(3H)-one | n.d. | n.d. | 203.87 ± 1.93 | n.d. | n.d. |
| 3-ethyl-2,5-dimethylpyrazine | n.d. | 336.66 ± 15.68 | n.d. | n.d. | n.d. |
| **furfural** | 315.15 ± 3.75 | 2308.37 ± 35.29 | 949.28 ± 67.8 | **3748.97 ± 650** | 3200.33 ± 300.82 |
| furfuryl acetate | n.d. | n.d. | n.d. | 56.48 ± 9.8 | 238.97 ± 27.1 |
| 5-methylfurfural | 87.14 ± 5.57 | n.d. | **5505.51 ± 196.52** | 142.29 ± 8.64 | n.d. |
| 2-furanmethanol | 54.82 ± 4.6 | 689.14 ± 8.25 | 855.36 ± 67.61 | 634.93 ± 31.26 | 1245.43 ± 45.7 |
| 5-methyl-2(5H)-furanone | n.d. | n.d. | 629.98 ± 42.09 | n.d. | n.d. |
| 5-methylfurfuryl alcohol | n.d. | n.d. | 135.24 ± 6.08 | n.d. | n.d. |
| **furaneol** | n.d. | 76.8 ± 2.85 | **3583.32 ± 166.85** | 146.91 ± 38.33 | 117.16 ± 1.13 |

**Sucrose is the only sugar of the five that makes pyrazines at all** — because it is the only one
that is not a reducing sugar, so this is the opposite of what a naive reading would predict, and
the paper does not comment on it (Flags 6).

### Table 4. "Concentrations of volatile compounds obtained from coffee bean analogue (CBA) with different formulations"

Footnote as printed: "a n.d., not detectable; AA, a mixture of Ile, Leu, Val, Asp, Gly and Ser, and
the proportion of each amino acid was the same; S, sucrose; CA, chlorogenic acid. Numbers after the
ingredient abbreviations represent their content in the final extrudate (w/w). Mean concentration
expressed in parts per billion with RSD (%) below 15 % (n = 3)." **All six precursor formulations
carry 10 % sunflower oil and 1 % rhamnose; the Blank is chickpea flour and water only.**

| compound | Blank | AA1/S2/CA1 | AA1/S2 | AA1/S2/CA5 | AA1/S6/CA5 | AA1/S6/CA1 | AA5/S6/CA1 |
|---|---|---|---|---|---|---|---|
| 3-methylfuran | n.d. | n.d. | 10.85 ± 1.35 | 10.46 ± 0.48 | 10 ± 0.98 | 10.21 ± 0.36 | 4.52 ± 0.18 |
| 2-methylbutanal | n.d. | 123.84 ± 13.72 | 170.87 ± 19.75 | 119.38 ± 14.5 | 79.26 ± 0.85 | 69.93 ± 4.5 | **652.04 ± 38.18** |
| methylpyrazine | **293.45 ± 9.87** | 841.35 ± 14.35 | **1358.91 ± 167.44** | 410.44 ± 14.87 | 36.13 ± 5.18 | 468.84 ± 67.93 | 673.35 ± 32.59 |
| 2,5-dimethylpyrazine | n.d. | 223.98 ± 11.13 | 201.94 ± 25.45 | 32.96 ± 3.08 | 5.87 ± 0.52 | 42.7 ± 2.4 | 225.61 ± 24.43 |
| 2,3-dimethylpyrazine | n.d. | 129.5 ± 7.2 | 83.13 ± 12.04 | 59.14 ± 1.72 | 65.72 ± 0.29 | 142.82 ± 16.26 | 158.4 ± 8 |
| 2-ethyl-6-methylpyrazine | n.d. | 376 ± 18.12 | 419.15 ± 24.8 | 218.65 ± 7.3 | 102.84 ± 1.15 | 295.58 ± 12.32 | 534.46 ± 24.06 |
| 2-ethyl-5-methylpyrazine | n.d. | 421.67 ± 6.98 | 513.74 ± 22.54 | 340.78 ± 7.87 | 178.85 ± 10.54 | 322.58 ± 19.73 | 383.06 ± 18.85 |
| 2-ethyl-3-methylpyrazine | n.d. | 73.64 ± 7.27 | 80.09 ± 7.48 | 77.26 ± 4.55 | 45.5 ± 0.66 | 84.44 ± 4.27 | **158.4 ± 8** (Flags 5) |
| trimethylpyrazine | **101.74 ± 3.76** | 651.75 ± 35.15 | 885.59 ± 18.69 | 313 ± 15.81 | 119.25 ± 7.72 | 373.53 ± 29.58 | 686.65 ± 50.07 |
| 3-ethyl-2,5-dimethylpyrazine | **129.06 ± 8.34** | 698.46 ± 22.2 | 970.36 ± 27.94 | 301.14 ± 15.03 | 76.12 ± 4.57 | 385.35 ± 38.9 | **1407.86 ± 60.41** |
| furfural | **163.97 ± 12.75** | 715.98 ± 40.63 | 1081.92 ± 22.31 | 1827.56 ± 134.33 | **4158.53 ± 60.21** | 2007.71 ± 109.38 | 553.75 ± 55.95 |
| 3,5-dimethyl-2-isobutylpyrazine | n.d. | n.d. | n.d. | n.d. | n.d. | n.d. | 650 ± 29.47 |
| 5-methylfurfural (printed "5-methyifurfural") | **26.89 ± 1.76** | 203.99 ± 28.34 | 201.25 ± 6.45 | 1431.34 ± 67.13 | **3364.01 ± 58.67** | 985.72 ± 65.09 | n.d. |
| 2-isoamyl-6-methylpyrazine | n.d. | n.d. | n.d. | n.d. | n.d. | n.d. | 521.08 ± 20.45 |
| 2,5-dimethyl-3-(3-methylbutyl)pyrazine | n.d. | n.d. | n.d. | n.d. | n.d. | n.d. | 549.71 ± 12 |
| 2-furanmethanol | **41.44 ± 4.18** | 289.1 ± 22.4 | 410.56 ± 19.84 | 430.27 ± 24.81 | 527.3 ± 15.49 | **794.16 ± 57.24** | 268.91 ± 12.4 |
| furaneol | n.d. | 169.54 ± 1.82 | 173.45 ± 10.87 | 168.08 ± 9.14 | 125.25 ± 7.17 | 134.96 ± 23.02 | 45.09 ± 4.17 |

### Numbers printed only in the running text

| quantity | value | where |
|---|---|---|
| coffee reference | 27 volatiles identified, **22 with OAV >= 1** | section 3.2 (Table S3, **not on disk**) |
| coffee reference, furans and derivatives | **8607 ± 559 ppb** total | section 3.2 |
| coffee reference, pyrazines | **6875 ± 706 ppb** total | section 3.2 |
| coffee reference, heterocyclic N | **2525 ± 157 ppb** | section 3.2 |
| coffee reference, ketones | **2009 ± 95 ppb** | section 3.2 |
| coffee reference, top OAVs | **4-vinylguaiacol 666**; furfurylmethyl sulfide 352; 5-methylfurfural 344; furaneol 212; 3-ethyl-2,5-dimethylpyrazine 171; pyridine 110 | section 3.2 |
| coffee reference, OAV > 10 | guaiacol, phenylethyl alcohol, 2-methylpyrazine, 2-methylbutanal, 2-ethyl-6-methylpyrazine, 2-furanmethanol, acetate, butyrolactone, 2-ethyl-5-methylpyrazine | section 3.2 |
| single-amino-acid batch | **21 pyrazines** across all sixteen amino acids | section 3.3 (Table S4, **not on disk**) |
| amino acids giving pyrazines in most samples | methylpyrazine, 2,5-dimethylpyrazine, 2-ethyl-6-methylpyrazine, 2-ethyl-5-methylpyrazine, 3-ethyl-2,5-dimethylpyrazine | section 3.3 |
| single-sugar batch | **17 volatiles**: 10 furans and derivatives, 5 pyrazines, 2 diketones; rhamnose gave the most compounds, then sucrose | section 3.4 |
| multiple-precursor batch | **17 volatiles**: 11 pyrazines, 5 furanic compounds, 1 aldehyde; **12 detected in all precursor formulations**; methylpyrazine, trimethylpyrazine, 3-ethyl-2,5-dimethylpyrazine, furfural and 2-furanmethanol detected in every sample **including the Blank** | section 3.5 |
| pH, highest chlorogenic acid | **5.3** for CBA/AA1/S6/CA5 | section 3.5 (Fig. S2, **not on disk**) |
| pH, no chlorogenic acid | **6.3** | section 3.5.2 (Fig. S2, **not on disk**) |
| effect of removing chlorogenic acid | 3-ethyl-2,5-dimethylpyrazine **+38.9 %** and methylpyrazine **+61.5 %** in AA1/S2 vs AA1/S2/CA1 | section 3.5.1 |
| effect of 1 % -> 5 % chlorogenic acid | total pyrazines **3416.36 ± 122.40 -> 1753.36 ± 70.23 ppb** | section 3.5.1 |
| effect of 2 % -> 6 % sucrose on pyrazines | methylpyrazine **-44.3 %**, 2,5-dimethylpyrazine **-80.9 %** | section 3.5.1 (the text calls the low level "1 %"; Table 1 says 2 % — Flags 5) |
| effect of 1 % -> 5 % amino acid mixture | "the total content of pyrazine compounds **nearly doubled**" (my sum gives **2.81x** — Flags 5) | section 3.5.1 |
| effect of chlorogenic acid on furans | at 5 % CA, total furans up **nearly 3-fold** over 1 % CA; furfural **2.6-fold**, 5-methylfurfural **7.0-fold** | section 3.5.2 |
| effect of 2 % -> 6 % sucrose on furans | total furans **1378.61 ± 93.21 -> 3932.76 ± 255.1 ppb**; furfural and 2-furanmethanol each up about 3-fold | section 3.5.2 |
| effect of raising amino acid to 5 % on furans | total furans down by **about 1.5 times** | section 3.5.2 |
| Pearson correlations among the 12 common volatiles | 2,5-dimethylpyrazine with 3-ethyl-2,5-dimethylpyrazine **r = 0.88** and with 2-ethyl-6-methylpyrazine **r = 0.90**; 2-ethyl-5-methylpyrazine with 2-ethyl-3-methylpyrazine **r = 0.83** and with trimethylpyrazine **r = 0.95**; furaneol with 2-methylbutanal **r = -0.88**; furfural with methylpyrazine **r = -0.75**, with 2-ethyl-6-methylpyrazine **r = -0.90**, with 3-ethyl-2,5-dimethylpyrazine **r = -0.81**, with 2-ethyl-5-methylpyrazine **r = -0.86** | section 3.5, Fig. 4 |
| oil effect on volatiles | "**Neither the type nor the concentration of oil significantly influenced the number or diversity of volatile compounds produced**" | section 3.1.1 (Table S2, **not on disk**) |
| colour PCA | PC1 + PC2 explain **98.18 %**; oil-containing analogues cluster nearer the coffee reference | section 3.1.1, Fig. S1 (**not on disk**) |
| relative oxidation rates quoted | oleic acid **k = 10**, linoleic acid **k = 100** (from Erickson & List 1985 / Talbot 2016, not measured here) | section 3.1.1 |
| sunflower oil density quoted | 0.93 g/cm^3 | section 3.1.1 |
| green coffee composition quoted | sucrose **3.8-10.7 %** dry weight; chlorogenic acids **5-9 %** (Caporaso 2018, not measured here) | Introduction |

**Figure-only:** Figure 2 (expansion ratio and bulk density of every extrudate), Figure 3 (total
furans and total pyrazines by formulation) and Figure 4 (the correlation heat map, of which only the
nine coefficients quoted above are printed). Figure 1 is a process scheme with no numbers beyond
those already in the Methods.

### Arithmetic on the printed numbers (all mine)

1. **Table 1 closes exactly.** All seven formulations sum to 100.0 % w/w. A small point, and rare
   enough to be worth recording.
2. **The paper's own totals reproduce from Table 4 to the last decimal.** Total furans in
   AA1/S2/CA1 = 715.98 + 203.99 + 289.10 + 169.54 = **1378.61**, exactly the printed
   1378.61 ± 93.21. Total furans in AA1/S6/CA1 = 10.21 + 2007.71 + 985.72 + 794.16 + 134.96 =
   **3932.76**, exactly the printed 3932.76 ± 255.1. Total pyrazines in AA1/S2/CA1 = 841.35 +
   223.98 + 129.5 + 376 + 421.67 + 73.64 + 651.75 + 698.46 = **3416.35** against the printed
   3416.36; in AA1/S2/CA5 = 410.44 + 32.96 + 59.14 + 218.65 + 340.78 + 77.26 + 313 + 301.14 =
   **1753.37** against the printed 1753.36. **Four independent checks, all to within one unit in the
   second decimal.** Whatever else is true of these numbers, the table and the text describe the
   same data.
3. **The percentage claims also reproduce.** 970.36/698.46 = 1.389 (**+38.9 %**, as claimed);
   1358.91/841.35 = 1.615 (**+61.5 %**); 468.84/841.35 = 0.557 (**-44.3 %**); 42.7/223.98 = 0.191
   (**-80.9 %**). All four exact.
4. **One claim does not reproduce.** Total pyrazines in AA1/S6/CA1 = 2115.84 and in AA5/S6/CA1 =
   5948.58, a ratio of **2.81**, which the text calls "nearly doubled". It nearly tripled.
5. **The chickpea-flour blank's own Maillard output.** Summing the Blank column: **756.55 ppb** of
   volatiles in total, of which pyrazines 524.25 and furans 232.30 (mine). Against AA1/S2/CA1's
   5089.7 ppb (mine), the flour alone accounts for **15 %** of the total — so the added precursors
   do most, but not all, of the work, and a model that charges only the added precursors will
   under-predict by about that much in this system.
6. **The chlorogenic acid effect is a real dose-response and it goes both ways.** Pyrazines at
   0 / 1 / 5 % CA (holding AA at 1 % and sucrose at 2 %): **3812.72 / 3416.35 / 1753.37 ppb**
   (mine). Furans over the same series: **1866.42 / 1378.61 / 3897.71 ppb** (mine). So going from
   1 % to 5 % CA **halves the pyrazines and nearly triples the furans**, and the no-CA point sits
   above the 1 % point on pyrazines and above it on furans too — i.e. **the pyrazine response is
   monotone in CA but the furan response is not**, with a minimum at 1 %. The paper describes the
   furan side as a straightforward positive effect of CA and does not note the non-monotonicity.
7. **Furfural and the pyrazines really do trade off.** Across the six precursor formulations,
   furfural runs 715.98, 1081.92, 1827.56, 4158.53, 2007.71, 553.75 while total pyrazines run
   3416.35, 3812.72, 1753.37, 630.28, 2115.84, 5948.58 (mine). The rank correlation is strongly
   negative and the extremes are opposite: AA1/S6/CA5 has the most furfural and the fewest
   pyrazines; AA5/S6/CA1 has the most pyrazines and the least furfural. This is the paper's central
   observed pattern and it holds up in the raw table.
8. **The Strecker aldehydes are one-to-one with their amino acid, and nothing else makes them.**
   Table 2: isobutyraldehyde appears only with valine (1446.91 ppb), 2-methylbutanal only with
   isoleucine (1156.53), 3-methylbutanal only with leucine (1231.76). **No aldehyde is detected in
   the Asp, Gly or Ser samples at all.** That is a clean qualitative result, and it is the kind of
   thing the roadmap's amino-acid-identity programme wants — but at 10 % w/w of a single amino acid
   it is a very heavy charge, and the quantification cannot support a threefold comparison
   (Flags 1).
9. **A cross-link to `conti2025b_extraction.md`.** This paper prints
   **2,5-dimethyl-3-(3-methylbutyl)pyrazine** at 549.71 ± 12 ppb in the AA5/S6/CA1 analogue. That
   is the same compound Conti 2025b's running text names twice as significantly different between
   its two extrudates but never lists in its own Table 2. Two independent extruded plant-protein
   studies name it; only this one prints a number for it, and the registry has no id for it.

## 4. Numbers the repository can use

All rows: chickpea flour (21 % protein, 6.7 % fat, 58 % carbohydrate) as the matrix; twin-screw
extrusion at 100 C in all eight blocks, 600 rpm, 1.41 kg/h solids and 623 mL/h water, 6 mm die,
extrudate at 33 ± 1 % moisture; dried 45 C for 10 h to 11 ± 1 % moisture; **roasted in a fluidised
bed at 170 C for 20 min**; ground, sieved at 1.25 mm; HS-SPME 20 min at 60 C on 0.25 g;
GC-MS on a Stabilwax-DA; n = 3, RSD below 15 %; basis **ppb (ug/kg) of the ground roasted powder**.

| compound | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| **hexanal, 2-pentylfuran, 1-octen-3-ol, nonanal, 1-hexanol, any methoxypyrazine** | **not reported** | — | — | — | **absent — no lipid-oxidation volatile appears in any table of this paper** |
| **every value in Tables 2, 3 and 4** | see the tables above | printed as **ppb of the sample powder**; **a peak area ratioed to 4-methyl-1-pentanol with no authentic standard and no response factor** | as above | Tables 2, 3, 4 pp. 5-7 | **peak_area_only** (Flags 1) |
| the chickpea-flour Blank | methylpyrazine 293.45 ± 9.87; trimethylpyrazine 101.74 ± 3.76; 3-ethyl-2,5-dimethylpyrazine 129.06 ± 8.34; furfural 163.97 ± 12.75; 5-methylfurfural 26.89 ± 1.76; 2-furanmethanol 41.44 ± 4.18; **everything else n.d.** | ppb as above | 67 % chickpea flour + 33 % water, no oil, no added precursor, 170 C / 20 min | Table 4, Blank column | peak_area_only — **the most useful single column in the paper: a legume flour's own dry-roast Maillard output** |
| Strecker aldehyde from valine | isobutyraldehyde 1446.91 ± 18.64 | ppb as above | chickpea flour + 10 % L-valine, no oil, 170 C / 20 min | Table 2 | peak_area_only |
| Strecker aldehyde from isoleucine | 2-methylbutanal 1156.53 ± 66.29 | ppb as above | chickpea flour + 10 % L-isoleucine | Table 2 | peak_area_only |
| Strecker aldehyde from leucine | 3-methylbutanal 1231.76 ± 22.8 | ppb as above | chickpea flour + 10 % L-leucine | Table 2 | peak_area_only |
| the pyrazine fingerprint of six amino acids | 22 compounds x 6 amino acids, see Table 2 | ppb as above | chickpea flour + one amino acid at 10 %, no oil | Table 2 | peak_area_only |
| the volatile output of five sugars | 17 compounds x 5 sugars, see Table 3 | ppb as above | chickpea flour + one sugar at 10 %, no oil | Table 3 | peak_area_only |
| rhamnose as a furaneol precursor | furaneol 3583.32 ± 166.85 with rhamnose against 76.8 (sucrose), 146.91 (xylose), 117.16 (fructose), n.d. (glucose) | ppb as above | as above | Table 3 | peak_area_only — a **47-fold** within-study contrast (mine) |
| rhamnose as a 5-methylfurfural precursor | 5505.51 ± 196.52 against 142.29 (xylose), 87.14 (glucose), n.d. (sucrose, fructose) | ppb as above | as above | Table 3 | peak_area_only |
| chlorogenic acid dose-response on pyrazines | total pyrazines **3812.72 (0 % CA) / 3416.35 (1 %) / 1753.37 (5 %)** | ppb as above, my sums | AA 1 %, sucrose 2 %, rhamnose 1 %, oil 10 % | derived from Table 4 (mine); the 1 % and 5 % totals are printed in the text | within_study_ratio |
| chlorogenic acid dose-response on furans | total furans **1866.42 (0 %) / 1378.61 (1 %) / 3897.71 (5 %)** | ppb as above, my sums | as above | derived (mine) | within_study_ratio — **non-monotone** (Flags 6) |
| sucrose dose-response | 2 % -> 6 %: methylpyrazine -44.3 %, 2,5-dimethylpyrazine -80.9 %, total furans 1378.61 -> 3932.76 ppb | dimensionless and ppb | AA 1 %, CA 1 % | Table 4 and text | within_study_ratio |
| amino-acid dose-response | 1 % -> 5 %: total pyrazines **2115.84 -> 5948.58 ppb (2.81x, mine)**; total furans down about 1.5x | ppb | sucrose 6 %, CA 1 % | derived from Table 4 (mine) | within_study_ratio |
| pH of the analogues | **5.3** at 5 % CA / 6 % sucrose; **6.3** with no CA | pH of a boiling-water extract, 2 g in 15 mL, overnight | the multiple-precursor formulations | text; the full series is Fig. S2, **not on disk** | measured_level (two points only) |
| coffee reference class totals | furans 8607 ± 559; pyrazines 6875 ± 706; heterocyclic N 2525 ± 157; ketones 2009 ± 95 | ppb | 100 % Arabica, Ekoplaza strength 8, as bought | text, section 3.2 (Table S3, **not on disk**) | peak_area_only |
| coffee reference OAVs | 4-vinylguaiacol 666; furfurylmethyl sulfide 352; 5-methylfurfural 344; furaneol 212; 3-ethyl-2,5-dimethylpyrazine 171; pyridine 110 | dimensionless | as above | text (Table S3, **not on disk**) | within_study_ratio — **the thresholds behind them are not on disk** |
| Pearson correlations, nine printed | see the text table above | r | across the six precursor formulations, 12 common compounds | text, Fig. 4 | within_study_ratio |
| oil type and level, effect on volatiles | none significant, on number or diversity | — | sunflower vs high-oleic sunflower, 10 % vs 15 % | text (Table S2, **not on disk**) | level_only |
| the extrusion and roast specification | see section 2 | — | — | sections 2.3, 2.4 | measured_level (process set points; **no residence time, no melt temperature**) |
| chickpea flour composition | 21 % protein, 6.7 % fat, 58 % carbohydrate | % w/w | supplier declaration | section 2.1 | level_only (not measured here) |
| sunflower oil fatty-acid split | 10 / 27 / 55 % saturated / mono / poly; high-oleic 8.7 / 51 / 33 % | % of fatty acids | supplier declaration | section 2.1 | level_only |
| expansion ratio, bulk density, colour | — | — | all oil variants | Fig. 2, Table S1 (**not on disk**) | **figure_only** / not available |
| total furans and pyrazines by formulation | — | ppb | six formulations plus Blank | Fig. 3 | **figure_only** (but reconstructable by summing Table 4, as above) |

### What this can and cannot be used for

**Can:** supply a chickpea-flour dry-roast blank; supply a dosed chlorogenic-acid effect on the
pyrazine-versus-furan partition in a solid matrix; supply a per-amino-acid and per-sugar
qualitative fingerprint under one fixed cook; and supply a complete extrusion-plus-roast process
description that a future benchmark could reuse. All of it as **ratios and presences**, never as
levels.

**Cannot:** contribute to the carried-volatile levels table. Cannot supply an absolute
concentration, because nothing was calibrated against an authentic standard. Cannot supply a rate,
because there is one roast temperature and one roast time. Cannot supply an odour-activity value,
because the thresholds are in a supplementary table that is not on disk.

## 5. Flags

1. **Every "ppb" in this paper is a peak area, and the quantification method is described in a
   single sentence.** "The concentration of the identified compound was carried out based on the
   internal standard used." The internal standard is 4-methyl-1-pentanol; **its concentration and
   the amount added are never printed**. The chemicals list contains no pyrazine, no furan and no
   aldehyde standard. So a response factor of 1 is applied to every one of the ~40 reported
   compounds against a C6 alcohol. Under the house rule these are `peak_area_only`. **The
   within-study ratios are the usable quantity** — and the paper's own totals reproduce exactly
   from the tables (section 3, arithmetic 2), so the internal consistency is good even though the
   absolute scale is not established.
2. **No lipid-oxidation volatile is reported, in a system with 16.7 % total fat.** The chickpea
   flour is 6.7 % fat and the precursor formulations add 10 % sunflower oil, over half of which is
   linoleic acid, and the product is roasted at 170 C for 20 min. **Not one of hexanal, nonanal,
   2-pentylfuran, 1-octen-3-ol, the octadienones or the alkadienals appears in any table.** The
   authors go so far as to attribute the darker colour of the oiled samples partly to "lipid
   oxidation byproducts" reacting with amino acid residues, and to invoke the relative oxidation
   rates of oleic and linoleic acid — and then report no oxidation product. The paper does not say
   whether these compounds were absent, below detection, or excluded from the reported set. **This
   is the single most important thing to ask the authors.**
3. **The extrusion is set points, not thermal history.** Barrel blocks all at 100 C, 600 rpm,
   1.41 kg/h and 623 mL/h, die 6 mm, L/D 40. **No residence time, no melt temperature, no die
   pressure, no torque, no specific mechanical energy.** As with the Conti pair, the extrusion step
   cannot be reconstructed as a cook.
4. **One roast, and no time or temperature series.** 170 C for 20 min in a fluidised bed, on every
   sample. There is no second temperature, no second duration, no intermediate sample and no
   temperature log. **Nothing in this paper is resolved in time.** The only quantity that varies is
   composition.
5. **Five printed inconsistencies.** (i) **The decimal separator is a comma in Tables 2, 3 and 4**
   (e.g. "1446,91 ± 18,64"), except for one entry, 2-isobutyl-3-methylpyrazine in CBA-Ile, which is
   printed as "153.60 ± 2.67" with points. A reader who mixes the two conventions will be wrong by
   a factor of 100. (ii) In Table 4 the AA5/S6/CA1 column gives **2,3-dimethylpyrazine and
   2-ethyl-3-methylpyrazine both as exactly 158.4 ± 8** — two different compounds with byte-identical
   values, which is a duplication signal. (iii) Section 3.5.1 says "Two different concentrations
   (**1 % and 6 %**) of sucrose were used"; **Table 1 says 2 % and 6 %**, and the percentage changes
   quoted reproduce from the 2 % column, so the table is right. (iv) Section 3.5.2 says "In the
   absence of CA, which caused the pH to increase to 6.3, 3-methylfuran was detected at a
   concentration of 10.85 ± 1.35 ppb **in CBA/AA1/S2/CA5**" — but 10.85 is the value in the
   **CBA/AA1/S2** column (no chlorogenic acid), and CBA/AA1/S2/CA5 reads 10.46. The same paragraph
   says "the concentration of furanone in **CBA/AA1/S6/CA5** peaked at **173.45 ± 10.87 ppb**"; 173.45
   is furaneol in the **CBA/AA1/S2** column, and AA1/S6/CA5 reads 125.25. **Two sample labels are
   swapped in one paragraph**, and both swaps point at the no-chlorogenic-acid column. (v) The
   amino-acid effect on pyrazines is called "nearly doubled"; my sum of Table 4 gives 2.81x. The
   compound name "5-methyifurfural" in Table 4 is a typographical error for 5-methylfurfural.
6. **Two results the paper reports but does not explain, and one it does not notice.** (a) In
   Table 3, **sucrose is the only one of the five sugars that produces any pyrazine at all** — and
   sucrose is the only non-reducing sugar of the five. Glucose, rhamnose, xylose and fructose
   produce furans and diketones and no pyrazine whatever. The paper does not comment. (b) In
   Table 3, **rhamnose gives furaneol at 3583.32 ppb and 5-methylfurfural at 5505.51 ppb**, 30 to
   60 times any other sugar, which the paper does note and attributes to Illmann 2009. (c) My
   reconstruction of the chlorogenic-acid series (section 3, arithmetic 6) shows the **furan
   response is non-monotone in CA** — higher at 0 % than at 1 %, then much higher at 5 % — while the
   paper describes it as a straightforward positive effect.
7. **Most of the paper's evidence is in supplementary material that is not on disk.** Table S1
   (colour), **Table S2 (the oil-type volatile comparison, which is the sole basis for the claim
   that oil type does not matter)**, **Table S3 (the coffee reference's 27 compounds, their
   thresholds and their 22 OAVs)**, **Table S4 (21 pyrazines across all sixteen amino acids —
   ten amino acids' worth of data that never reach the main text)**, Figure S1 (colour PCA) and
   **Figure S2 (the pH of every formulation, which the mechanism in section 3.5 depends on)**.
   **The paper is open access (CC BY), so all of this is freely obtainable.**
8. **What this paper does NOT contain, and what to request.** No lipid-oxidation volatile; no
   measurement of the chickpea flour's own free amino acids, free sugars, moisture or fatty acids;
   no measurement of what survives the extruder before roasting (only the roasted analogue is
   analysed); no acrylamide, no HMF, no furan (the toxicant) despite the matrix and cook being
   ideal for all three; no residence time; no second roast condition; no odour threshold; no
   detection limit; no recovery. **To request:** (a) **the complete Supplementary Material**, which
   is open access and holds Tables S1-S4 and Figures S1-S2 — this is the highest-value, lowest-cost
   acquisition in the whole cluster; (b) whether hexanal and the other lipid-oxidation volatiles
   were looked for and, if so, their values; (c) the 4-methyl-1-pentanol loading, so the ppb could
   at least be put on a stated basis; (d) the extruder residence time; (e) the underlying values
   behind Figures 2, 3 and 4.
9. **Registry gaps against `data/keys/compounds.yml`.** Present and keyable: `methylpyrazine`,
   `2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`, `2_3_dimethylpyrazine`, `trimethylpyrazine`,
   `tetramethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine`, `furfural`, `hdmf` (furaneol),
   `2_3_butanedione`, `2_methylbutanal`, `3_methylbutanal`, `2_methylpropanal` (the paper's
   isobutyraldehyde), `4_vinylguaiacol` (the coffee reference's top odorant),
   `2_furfurylthiol` (named in the Introduction, not measured), **`chlorogenic_acid`** (the paper's
   dosed variable — a direct registry hit), and the class id `pyrazines`. **Absent:**
   2-ethyl-6-methylpyrazine, 2-ethyl-5-methylpyrazine, 2-ethyl-3-methylpyrazine, 2,6-diethylpyrazine,
   2,3-dimethyl-5-ethylpyrazine, 2-methyl-3,5-diethylpyrazine, 3,5-dimethyl-2-propylpyrazine,
   2,3-diethyl-5-methylpyrazine, 2,3,5-trimethyl-6-ethylpyrazine, 2-isoamyl-6-methylpyrazine,
   2-isobutyl-3-methylpyrazine, 3,5-dimethyl-2-isobutylpyrazine,
   **2,5-dimethyl-3-(3-methylbutyl)pyrazine** (also named by Conti 2025b), 3-methylfuran,
   2,5-dimethylfuran, 2,3-pentanedione, 5-methylfurfural, 2-furanmethanol (furfuryl alcohol),
   5-methylfuran-2(3H)-one, 5-methyl-2(5H)-furanone, 5-methylfurfuryl alcohol, furfuryl acetate,
   guaiacol, pyridine and furfurylmethyl sulfide. **5-Methylfurfural and 2-furanmethanol are the
   two most conspicuous gaps**: both appear in every one of this paper's tables, 5-methylfurfural
   carries the third-highest OAV in the coffee reference (344), and neither has a registry id.
10. **Registry gaps against `data/species/off_flavour_targets.yml`.** That file holds hexanal,
    nonanal, 1-octen-3-ol, 2-pentylfuran, 1-hexanol and furfural. **Only furfural appears in this
    paper**, and it appears as a *desirable* coffee-like compound rather than an off-flavour — the
    YAML's own note calls it "Almond, bread, bready, sweet" with a 3000 ug/kg threshold, and this
    paper measures it at up to 4158.53 ppb in a deliberately roasted product. That is a useful
    reminder that the file's framing (minimise and trap) is matrix-specific: in a coffee analogue
    the same molecule is the target, not the defect. No change to the file is supported by this
    paper, but the tension is worth recording.
