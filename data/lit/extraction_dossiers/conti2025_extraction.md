# Conti 2025 — EXTRACTION (soy protein CONCENTRATE with 1.5 % w/w thiamine, single-screw extruded at three moisture/temperature pairs — 30 %/180 C, 34 %/160 C, 38 %/140 C — then hydrated, salted and optionally oiled; texture profile analysis plus an eight-panellist in-vivo mastication study of bolus moisture and particle geometry)

### A texture and oral-processing paper, not a chemistry one: it prints NO concentration of any volatile, and everything it resolves in time is mastication time in the mouth, not reaction time — but it names the extrusion conditions and the thiamine loading of a soy-protein-concentrate meat analogue precisely, and it is the anchor point for the companion Part II on tastant release.

**Source on disk:** `data/articles/Conti2025.pdf` (11 pp., Food Research International 208 (2025)
116169). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/Conti2025.txt`), whole file. Tables 1, 2, 3 and 4 came through clean and are
re-typed in full below. Figures 1 (texture profile), 2 (particle-area histogram), 3A-D and 4A-H
(bolus parameters), 5 (PCA) and 6 (multiple factor analysis) are images and are **figure-only**; the
running text quotes a number of individual values off them and those are transcribed, marked as
quoted-from-text. **Supplementary Table 1** (significance of factors per particle-area range) and
the rest of Appendix A are online-only and are **not on disk**. Repo status before this dossier:
neither this paper nor its companion `Conti2025b.pdf` had a dossier, and neither is cited anywhere
in `src/`, `data/` or `results/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Oral processing of meat-flavour textured soy proteins — **Part I: Bolus properties and relationships with the texture profile analysis of the products**" |
| Authors | Ana Carolina Conti (corresponding, ac.conti@unesp.br), Chantal Septier, Emmanuel Denimal, Helene Laboure, Christian Salles |
| Affiliations | (a) Sao Paulo State University (Unesp), Ibilce, Department of Food Engineering and Technology, Sao Jose do Rio Preto, Brazil; (b) Centre des Sciences du Gout et de l'Alimentation, CNRS / INRAE / Institut Agro / Universite Bourgogne Europe, Dijon, France; (c) Institut Agro, Dijon |
| Venue | **Food Research International 208 (2025) 116169**. Received 18 Sep 2024, revised 17 Jan 2025, accepted 9 Mar 2025, online 11 Mar 2025 |
| DOI | **10.1016/j.foodres.2025.116169** |
| Companion | "In Part II, the focus will be on the release of sodium and glutamic acid from the products to the saliva and on their sensory aspects." That is `Conti2025b.pdf` (Food Research International 218:116938) — see `conti2025b_extraction.md` |
| Naming | **SPC** = soy protein concentrate; **TSP** = textured soy protein, the extrudate as it leaves the die; **FTSP** = *flavoured* textured soy protein, the TSP after hydration, salting and (for half the samples) oiling. "M" = moisture of the SPC on a dry basis before extrusion. The quoted temperature is that of **Zone 5** of the extruder barrel |
| Ethics | Inserm Ethics Evaluation Committee No. 20-754bis, approved March 2021 |
| Funding | FAPESP 2019/20911-4; CNPq 303602/2022-8; Institut Agro Dijon |
| Data | "Data will be made available on request." |

**Both Conti papers identified.** This file is **Part I**, FRI 208:116169, bolus properties and
texture. `Conti2025b.pdf` is **Part II**, FRI 218:116938, the release of sodium and glutamic acid
into saliva. They share one experimental design, one panel and one set of six products; Part I
measures what the mouth does to the food, Part II measures what comes out of the food into the
saliva. **Neither measures a volatile compound.**

## 1. Why it matters — and the honest answer is: not much, for the levels table

The task that brings this paper into the cluster is `tasks/roadmap_for_scientists.md` section 5d,
Programme 7 part (ii): charge the isolate's own volatiles as declared inputs with their measured
levels and bands. **This paper contains no volatile measurement of any kind.** There is no GC, no
mass spectrometer, no sniffing port, no headspace and no odour threshold anywhere in it. The words
"aroma" and "flavour" appear only in the product's name and in the description of thiamine as a
"meat aroma precursor"; the measurement is a texture analyser and a scanner.

So the entry it makes in the Programme 7 levels table is a **blank with a reason**, and that is
worth recording explicitly so the next reader does not re-open the PDF hoping for numbers.

What it does give the repository, and it is not nothing:

1. **A named soy protein CONCENTRATE with its protein specification**: Arcon SM (ADM Foods &
   Wellness, Decatur, Illinois), minimum 70 g protein per 100 g on a dry basis. The roadmap's
   levels table needs the material named exactly, and "concentrate at >= 70 % protein" is a
   different row from "isolate at 80 % protein" (`trikusuma2020_extraction.md`) and from "whole
   seed flour" (`bi2020_extraction.md`).
2. **A thiamine loading with a stated purpose**: "1.5 % thiamine (w/w) was added to each soy
   protein concentrate two hours before extrusion", as thiamine hydrochloride, purity > 99 %,
   explicitly as "a meat-like aroma precursor". The repository carries `thiamine_availability` as a
   modifier id in `data/keys/compounds.yml`, and thiamine degradation is the origin of the sulfur
   ladder the k6a synthesis works on. **A 1.5 % w/w charge of thiamine into an extruder is a real,
   industrially plausible dose from a peer-reviewed source**, and the repository has no such
   figure otherwise.
3. **Three extrusion conditions in full**, each a moisture/temperature pair, with the whole screw
   geometry, feed rate, screw speed and zone temperatures printed. A wave that wants to model
   extrusion as a cook needs exactly this level of description, and it is rarely printed.
4. **A hydration protocol and the resulting mass gain**: the extrudate's weight increases by
   **435 %, 411 % and 406 %** for the three conditions, from which a rehydrated moisture can be
   reasoned about.
5. **It points to the paper that actually did the chemistry.** Milani, Menis-Henrique & Conti
   (2022), "Thiamine as a new ingredient for obtaining textured soy protein with meat odour",
   *Journal of Food ...* — and Milani & Conti (2024), "Textured soy protein with meat odor as an
   ..." — are cited as the source of the flavour work and of the "intermediate condition
   (34 % M / 160 C) previously determined as optimal". **Neither is on disk.** That is the
   acquisition this dossier most clearly justifies.

## 2. Methods as they matter to a model

- **The material, exactly as described.** "**Soy protein concentrate (SPC) from Arcon SM, with a
  minimum of 70 g/100 g protein (dry basis), was supplied by ADM Foods & Wellness (Decatur,
  Illinois, USA).**" It is a **protein concentrate** — not an isolate, not a flour, not a whole
  seed. No lipid, ash, moisture-as-received, sugar or free amino acid content is printed. No
  cultivar, no batch, no lot.
- **The added precursor.** "Thiamine (vitamin B1), used as a meat aroma precursor, was supplied by
  Sigma-Aldrich as **thiamine hydrochloride (purity > 99 %)**." Loading: "the moisture of the SPC
  was adjusted to the desired value, i.e., 30 %, 34 % and 38 % (**dry basis**), and then **1.5 %
  thiamine (w/w)** was added to each soy protein concentrate **two hours before extrusion** ... The
  portions were kept at room temperature until extrusion."
- **The extruder, in full.** RXPQ Labor 24 **single-screw** extruder (INBRAMAQ, Ribeirao Preto,
  Brazil), five independent heating zones; helicoidally grooved barrel; screw with a large step,
  one exit, **compression ratio 3.3:1**, **length-to-diameter ratio 15.5:1**; pre-die with holes of
  **5.8 mm**; die **3.6 mm** round hole; **feed rate 170 g/min**; **screw speed 216 rpm**; zone 1
  off (approximately 40 C), zone 2 **60 C**, zone 3 **80 C**; **zone 4 always 15 C lower than zone
  5**; the quoted condition temperature is zone 5. Assays were run from lowest to highest
  temperature "to facilitate temperature changes".
- **The three conditions.** (i) most severe: **30 % moisture / 180 C**; (ii) intermediate:
  **34 % moisture / 160 C**; (iii) least severe: **38 % moisture / 140 C**. The intermediate
  condition was taken from Milani et al. 2022 as optimal for thiamine-flavoured SPC; the other two
  bracket it. **No residence time, no die pressure, no specific mechanical energy and no melt
  temperature is printed** — only the barrel-zone set points (Flags 2).
- **Hydration, salting and oiling.** "Samples approximately **2 cm in length** were cut from the
  TSPs using a mold and hydrated. The water was heated at **100 C, the heating was turned off**, and
  the TSPs were immersed in the hot water (**solid-liquid ratio of 1:4**) for **15 min**. Then,
  **100 g of hydrated TSPs** were mixed manually with **salt (1.0 g/100 g)** and **monosodium
  glutamate (0.4 g/100 g)** ... TSPs from the same extrusion condition were separated into two
  portions: one was added with **7.0 g of vegetable oil/100 g of hydrated TSPs**, while the other
  was not." Six FTSPs in total; all prepared at the same time and **frozen at -18 C**. The
  vegetable oil is not identified beyond "vegetable oil" (Flags 3).
- **Serving state.** Samples of 2.0 +/- 0.1 g were taken from the freezer one hour before each
  session and "heated in an oven to reach **50 C in the middle**, which is an adequate temperature
  for consuming this product". The texture analysis was done at the same 50 C so that the
  instrument and the mouth saw the same product.
- **Texture profile analysis.** TA.XT/Plus/50 (Stable Micro Systems) with Texture Exponent 32;
  eight fully hydrated ~2 cm samples per FTSP; 50 mm plastic round probe; test speed **1 mm/s**;
  compression to **50 % of sample height**; **5 s between the two compressions**. Hardness,
  cohesiveness, springiness and chewiness reported.
- **Panel.** Eight healthy panellists (6 female, 2 male, 23-60 years). Stimulated salivary flux
  measured by chewing a 5 x 5 cm piece of Parafilm for 5 min in triplicate, saliva expelled and
  weighed, density taken as 1.0 g/mL. Normality of mastication assessed on a carrot (2 cm diameter
  x 1 cm) chewed ten times, in triplicate.
- **Bolus collection with a time design.** Stage 1: chew normally to just before swallowing, spit,
  and **record the time** — this per-panellist, per-sample time becomes the reference. Stage 2: chew
  the same sample for **1/3, 2/3 and 3/3 of that reference time** and spit each bolus. Mouth rinsed
  once with a little water into the same pot. Order randomised.
- **Bolus measurement.** Moisture by evaporation at **103 C for 24 h**. Image analysis: boluses made
  up to 160 mL with distilled water (plus **1.5 mL of 20 % sodium dodecyl sulphate** for the oiled
  samples, to stop the oil forming a halo), shaken five times, spread on a glass plate on an Epson
  Perfection V850 Pro scanner, **8-bit greyscale at 400 dpi**, analysed in MATLAB 2019b with
  adaptive thresholding (Bradley & Roth 2007, preferred over Otsu so low-contrast particles are not
  suppressed). Particles with a side below 0.4 mm (area below 0.16 mm^2) and particles touching the
  image edge were discarded. Per particle: convex area, maximum and minimum Feret diameter,
  circularity (0 to 1, 1 = a circle), and mean intensity (0 = white to 250 = black, read as a proxy
  for thickness). **One number of particles per bolus; every other parameter is reported as the
  median over the particles of that bolus.**
- **Extraction and quantification of chemical species: none.** There is no chemical analysis in this
  paper beyond bolus moisture by oven drying. **No authentic standards, no response factors and no
  concentration of anything appears.**
- **Basis of the numbers that do appear.** Mastication time in **seconds**; bolus moisture in
  **% (w/w) of the wet bolus**, by loss on drying at 103 C for 24 h; particle counts as **counts per
  bolus**; areas in **mm^2**; Feret diameters in **mm**; circularity and mean intensity
  dimensionless; salivary flux in **mL/min**; carrot d50 in **mm**.
- **Statistics.** General linear model with panellist, extrusion condition, oil addition and
  mastication time as factors plus the three two-level interactions of interest, followed by Tukey;
  Statistica 7.0; alpha = 0.05. PCA on the bolus particles and multiple factor analysis over four
  tables (texture profile, mastication time, bolus moisture, bolus particles at 3/3) in XLSTAT.

## 3. Tables re-typed

### Table 1. "Stimulated salivary flux and particle size of the carrot bolus (mean ± SD; n = 3)"

| panellist | salivary flux (mL/min) | particle size d50 (mm) |
|---|---|---|
| V1 | 2.15 +/- 0.19 | 3.58 +/- 0.23 |
| V2 | 1.55 +/- 0.33 | 6.30 +/- 0.75 |
| V3 | 0.94 +/- 0.13 | 6.87 +/- 0.28 |
| V4 | 1.26 +/- 0.25 | 6.94 +/- 0.77 |
| V5 | 0.72 +/- 0.38 | 5.78 +/- 0.69 |
| V6 | 3.15 +/- 0.30 | 5.77 +/- 0.50 |
| V7 | 1.16 +/- 0.05 | 7.43 +/- 0.48 |
| V8 | 1.45 +/- 0.22 | 5.98 +/- 0.36 |

### Table 2. "Mastication times just before swallowing (mean ± SD in s) of the flavoured textured soy proteins (n = 3 for each panellist and n = 24 for the group)"

Footnotes as printed: "M = moisture of the soy protein concentrate (dry basis) before extrusion.
Temperature in Zone 5 of the extruder barrel. Different letters in the same line indicate different
statistical means according to the Tukey test (p <= 0.05)."

| panellist | 30 % M / 180 C, no oil | 34 % M / 160 C, no oil | 38 % M / 140 C, no oil | 30 % M / 180 C, oil | 34 % M / 160 C, oil | 38 % M / 140 C, oil |
|---|---|---|---|---|---|---|
| V1 | 12.7 +/- 1.2 | 19.7 +/- 0.6 | 21.0 +/- 3.6 | 18.3 +/- 4.0 | 20.0 +/- 1.7 | 20.0 +/- 2.0 |
| V2 | 13.7 +/- 2.1 | 26.0 +/- 5.0 | 28.3 +/- 2.9 | 16.7 +/- 4.5 | 27.0 +/- 3.6 | 26.0 +/- 2.0 |
| V3 | 11.0 +/- 4.4 | 16.7 +/- 1.5 | 16.3 +/- 1.5 | 12.7 +/- 3.5 | 16.3 +/- 2.5 | 15.0 +/- 1.7 |
| V4 | 30.3 +/- 6.0 | 28.0 +/- 5.0 | 31.0 +/- 2.6 | 29.3 +/- 0.6 | 28.7 +/- 5.7 | 28.7 +/- 7.5 |
| V5 | 19.0 +/- 2.0 | 19.3 +/- 3.5 | 18.0 +/- 4.6 | 17.0 +/- 2.6 | 21.7 +/- 2.1 | 20.3 +/- 4.0 |
| V6 | 13.3 +/- 1.2 | 16.7 +/- 1.5 | 18.3 +/- 1.2 | 14.7 +/- 2.5 | 17.7 +/- 1.5 | 16.7 +/- 1.2 |
| V7 | 8.0 +/- 2.0 | 11.7 +/- 3.8 | 13.3 +/- 0.6 | 9.0 +/- 1.7 | 12.3 +/- 2.3 | 13.7 +/- 3.5 |
| V8 | 8.0 +/- 1.0 | 9.7 +/- 1.5 | 11.7 +/- 3.1 | 7.0 +/- 1.0 | 8.3 +/- 0.6 | 11.7 +/- 5.0 |
| **group** | **14.5 b +/- 7.4** | **18.5 a +/- 6.6** | **19.8 a +/- 6.9** | **15.6 b +/- 7.0** | **19.0 a +/- 7.0** | **19.0 a +/- 6.6** |

### Table 3. "Moisture (%) of the boluses (mean ± SD; n = 24) from the flavoured textured soy proteins"

Footnotes as printed: "M = moisture of the soy protein concentrate (dry basis) before extrusion.
Temperature in Zone 5 of the extruder barrel. Different lower-case letters in the same column
indicate different statistical means according to the Tukey test (p <= 0.05). Different capital
letters in the same line indicate different statistical means according to the General Linear Model
(p <= 0.05)."

| extrusion condition | without vegetable oil | with vegetable oil |
|---|---|---|
| 30 % M / 180 C | 84.4 bA +/- 2.9 | 82.3 bB +/- 3.1 |
| 34 % M / 160 C | 86.2 aA +/- 2.7 | 83.2 aB +/- 2.8 |
| 38 % M / 140 C | 86.0 aA +/- 3.0 | 84.0 aB +/- 3.8 |

### Table 4. "Significance of factors from the general linear model for each bolus parameter from the flavoured texturized soy proteins"

Footnote as printed: "*p <= 0.05; **p <= 0.01; ***p <= 0.001. ns = not significant."

| factor | number of particles | convex area | maximum Feret diameter | minimum Feret diameter | circularity | mean intensity |
|---|---|---|---|---|---|---|
| panellist | *** | *** | *** | *** | *** | *** |
| extrusion condition (EC) | *** | *** | *** | *** | *** | *** |
| vegetable oil (VO) | *** | *** | *** | *** | ** | ns |
| time (T) | *** | *** | *** | *** | * | *** |
| EC x VO | ns | *** | *** | * | ns | *** |
| EC x T | ns | ns | ns | ns | ns | ns |
| VO x T | ns | ns | ns | ns | ns | ns |

### Numbers quoted in the running text off the figures

These are printed in the prose but their table is a figure, so they are recorded here as
text-quoted values with their figure of origin.

| quantity | value | condition | figure quoted from |
|---|---|---|---|
| hardness | 2.4 N (no oil), 2.3 N (oil) | 30 % M / 180 C | Fig. 1A |
| hardness | "much lower (~2 N)" at 180 C vs "> 6 N" at 160 or 140 C | all | Fig. 1A |
| hardness, highest | **10.2 N** | 34 % M / 160 C **with** oil | Fig. 1A |
| chewiness | 1.5 N (no oil), 1.4 N (oil) | 30 % M / 180 C | Fig. 1A |
| cohesiveness, lowest | **0.7** | 30 % M / 180 C | Fig. 1B |
| springiness | not significantly different across any factor | all | Fig. 1B |
| number of particles | **2180** | 30 % M / 180 C | Fig. 4A |
| circularity, lowest | **0.49** | 30 % M / 180 C | Fig. 4B |
| number of particles, with oil | **1826** | oil averaged over conditions | Fig. 4C |
| circularity, with oil | **0.53** | oil averaged over conditions | text ("figure not shown") |
| number of particles vs chewing time | **1303 at 1/3 -> 2438 at 3/3** | averaged | Fig. 4D |
| convex area at 1/3 | **0.56 mm^2** | averaged | Fig. 4E |
| maximum Feret diameter at 1/3 | **1.22 mm** | averaged | Fig. 4E family |
| minimum Feret diameter at 1/3 | **0.68 mm** | averaged | Fig. 4E family |
| mean intensity at 1/3 | **182** | averaged | text |
| circularity at 1/3 | **0.55** | averaged | text ("graph not shown") |
| convex area with oil | 0.64 mm^2 (30 % M / 180 C) and 0.63 mm^2 (34 % M / 160 C) | with oil | Fig. 4F |
| maximum Feret diameter, effect of oil | +6 % (30 % M / 180 C) and +3 % (34 % M / 160 C) | with oil | Fig. 4G |
| particles 0.13-0.26 mm^2 | about **32 % more** under 30 % M / 180 C | — | Fig. 3A |
| particles 4.37-7.59 mm^2 | **61** (30 % M / 180 C) falling to **50** (38 % M / 140 C) | — | text |
| effect of oil on small particles | **-12 %** (0.13-0.26 mm^2) and **-9 %** (0.26-0.48 mm^2); **+5 particles** in 22.91-39.81 mm^2 | with oil | Fig. 3B, text |
| effect of chewing time on particle count | **+70 %** across all area ranges from 1/3 to 3/3; **+90 %** in the 4.37-7.59 mm^2 range with oil | — | Fig. 3C, 3D |
| **hydration mass gain of the TSPs** | **+435 % (30 % M / 180 C), +411 % (34 % M / 160 C), +406 % (38 % M / 140 C)** | 15 min in water taken off the boil at 1:4 | text, section 3.4 |
| PCA variance explained | 83 % in the first two components | — | Fig. 5 |
| MFA variance explained | 89.4 % in two factors | — | Fig. 6 |
| atypical carrot boluses | 87.5 % (7 of 8 panellists) had d50 above the 4.0 mm normality cut-off | — | text |
| F statistics | mastication time: F_panellist(7;131) = 67.10 p < 0.001, F_extrusion(2;131) = 23.85 p < 0.001; bolus moisture: F_panellist(7;131) = 53.70, F_extrusion(2;131) = 13.30, F_oil(1;131) = 75.70, all p < 0.001; hardness and chewiness: F_extrusion*oil(2;42) = 8.56 p < 0.001 and 7.43 p = 0.002; cohesiveness: F_extrusion(2;42) = 16.40 p < 0.001 | — | text |

**Figure-only in this paper:** Figures 1 (texture profile, four variables x six samples), 2
(particle-count histograms by area range), 3A-D, 4A-H, 5A-B and 6A-B. Everything not quoted in the
prose above is not transcribed.

### Arithmetic on the printed numbers (all mine)

1. **The mastication-time effect is small next to the panellist effect.** The group means span
   14.5 to 19.8 s, a range of **5.3 s**, while individual panellist means at a fixed sample span
   8.0 to 30.3 s, a range of **22.3 s**. That is why panellist is a factor in the model, and it is
   why any release measurement built on this design (Part II) has to be read within-panellist.
2. **Vegetable oil costs the bolus about 2 moisture points.** Averaging the three conditions:
   85.5 % without oil against 83.2 % with, a **2.3-point** fall (mine). The paper quotes the
   extremes, "from 82.3 % to 84.0 %", which is the with-oil column's own range, not the effect
   size.
3. **The severe condition costs about 1.7 moisture points.** 84.4 vs the mean of 86.2 and 86.0
   without oil = **-1.75 points**; 82.3 vs the mean of 83.2 and 84.0 with oil = **-1.3 points**.
4. **A rough moisture of the hydrated FTSP before chewing (mine, and it is only a bound).** A mass
   gain of 435 % means 1 g of extrudate becomes 5.35 g. If the extrudate leaving the die were dry,
   the hydrated product would be (5.35 - 1)/5.35 = **81.3 % water**; at 411 % and 406 % the same
   arithmetic gives **80.4 %** and **80.2 %**. The extrudate is **not** dry — its moisture is never
   printed — so these are lower bounds on the hydrated moisture. They sit just below the measured
   bolus moistures (82.3 to 86.2 %), which is the expected direction once saliva is added, so the
   two sets are consistent. **This is arithmetic on an assumption and is marked as such.**
5. **The three extrusion conditions confound moisture with temperature by design.** Moisture rises
   30 -> 34 -> 38 % exactly as temperature falls 180 -> 160 -> 140 C. **No single-variable
   comparison exists in this paper**, so no effect can be attributed to temperature alone or to
   moisture alone. The authors say as much implicitly by naming the conditions "most severe",
   "intermediate" and "less severe". For a model this matters: the design gives one severity axis,
   not two.

## 4. Numbers the repository can use

All rows share: soy protein concentrate Arcon SM (ADM), >= 70 g protein/100 g dry basis, plus
1.5 % w/w thiamine hydrochloride added two hours before extrusion; single-screw RXPQ Labor 24;
hydrated 15 min at 1:4 in water taken off the boil; salted at 1.0 g NaCl and 0.4 g MSG per 100 g of
hydrated product; half the samples oiled at 7.0 g/100 g; served at 50 C.

| quantity | value | unit and basis | material and conditions | source location | evidence class |
|---|---|---|---|---|---|
| **any carried volatile (hexanal, 2-pentylfuran, 1-octen-3-ol, a methoxypyrazine, or any other)** | **not measured** | — | — | — | **absent — this paper contains no volatile analysis at all** |
| soy protein concentrate protein content | >= 70 | g protein / 100 g, dry basis | Arcon SM, ADM Foods & Wellness | Methods 2.1 | level_only (a supplier specification, not a measurement made here) |
| thiamine loading before extrusion | 1.5 | % w/w on the SPC | as thiamine hydrochloride, > 99 % pure, added 2 h before extrusion | Methods 2.2 | measured_level (a formulation, exactly stated) |
| extrusion moisture / zone-5 temperature pairs | 30 % / 180 C; 34 % / 160 C; 38 % / 140 C | % moisture dry basis; C | single-screw, die 3.6 mm, feed 170 g/min, 216 rpm, L/D 15.5:1, compression 3.3:1, zones 1-3 at ~40 / 60 / 80 C, zone 4 = zone 5 minus 15 C | Methods 2.2 | measured_level (process set points; **no residence time or melt temperature**) |
| salt and MSG addition | 1.0 g NaCl and 0.4 g monosodium glutamate per 100 g hydrated TSP | g / 100 g | after hydration, mixed by hand | Methods 2.2 | measured_level |
| vegetable oil addition | 7.0 | g / 100 g hydrated TSP | oil unidentified; added after hydration, so "free" in the product | Methods 2.2 | measured_level |
| hydration mass gain | +435 % / +411 % / +406 % | % of extrudate mass | 30 %/180 C, 34 %/160 C, 38 %/140 C; 15 min at 1:4 in off-boil water | text, section 3.4 | measured_level |
| implied hydrated moisture (lower bound) | 81.3 % / 80.4 % / 80.2 % | % w/w water | assumes the extrudate leaving the die is dry, which it is not | derived (mine) | **derived_assumption** |
| bolus moisture, six products | 84.4 / 86.2 / 86.0 (no oil) and 82.3 / 83.2 / 84.0 (oil) +/- 2.7 to 3.8 | % w/w of wet bolus, loss on drying 103 C / 24 h | n = 24 boluses each, eight panellists | Table 3 | measured_level |
| mastication time to swallow, six products | 14.5 / 18.5 / 19.8 (no oil) and 15.6 / 19.0 / 19.0 (oil) +/- 6.6 to 7.4 | s | group means, n = 24 | Table 2 | measured_level |
| mastication time per panellist | 7.0 to 31.0 | s | 8 panellists x 6 products, n = 3 each | Table 2 | measured_level |
| stimulated salivary flux | 0.72 to 3.15 | mL/min | Parafilm chewing, 5 min, triplicate, saliva density taken as 1.0 g/mL | Table 1 | measured_level |
| carrot bolus d50 | 3.58 to 7.43 | mm | 10 chews, triplicate | Table 1 | measured_level |
| hardness | ~2 N at 180 C; > 6 N at 160 and 140 C; **10.2 N** at 34 % M / 160 C with oil; 2.4 and 2.3 N at 30 % M / 180 C | N | TPA at 50 C, 50 mm probe, 1 mm/s, 50 % compression, 5 s between cycles, n = 8 | text quoting Fig. 1A | **figure_only, values quoted in the text** |
| chewiness | 1.5 N (no oil) and 1.4 N (oil) at 30 % M / 180 C | N | as above | text quoting Fig. 1A | figure_only, quoted |
| cohesiveness | 0.7 at 30 % M / 180 C (the lowest) | dimensionless | as above | text quoting Fig. 1B | figure_only, quoted |
| springiness | no significant difference across any factor | — | as above | Fig. 1B | figure_only |
| number of bolus particles | 2180 (30 % M / 180 C); 1826 (with oil); 1303 at 1/3 chewing time rising to 2438 at 3/3 | counts per bolus | image analysis at 400 dpi, particles below 0.16 mm^2 discarded | text quoting Figs. 4A, 4C, 4D | figure_only, quoted |
| particle circularity | 0.49 (30 % M / 180 C); 0.53 (with oil); 0.55 at 1/3 chewing time | dimensionless, 0 to 1 | as above, median per bolus | text quoting Fig. 4B and unshown graphs | figure_only, quoted |
| particle convex area, Feret diameters, mean intensity at 1/3 | 0.56 mm^2; 1.22 mm max Feret; 0.68 mm min Feret; 182 mean intensity | mm^2, mm, greyscale 0-250 | as above | text quoting Fig. 4E | figure_only, quoted |
| significance of every factor on every bolus parameter | see Table 4 above | p-level marks | GLM, alpha = 0.05 | Table 4 | level_only |
| complete particle-area distributions, PCA, MFA | — | — | — | Figs. 2, 3, 4, 5, 6 | **figure_only** |
| significance per particle-area range | — | — | — | Supplementary Table 1 | **not on disk** |

### What is resolved in time, and what is not

There **is** a time axis in this paper, and it must not be mistaken for a reaction time. It is
**mastication time in the human mouth**, expressed as fractions 1/3, 2/3 and 3/3 of each panellist's
own time-to-swallow for that product, which ranges from 7 to 31 s. Over that axis the paper resolves
particle count (1303 -> 2438), convex area, Feret diameters, circularity and mean intensity. It
resolves **no chemical quantity** in time: nothing is measured at more than one point of any
process, and the extrusion itself is a single pass at each condition with **no residence time
printed**. For the roadmap's question (iii) the answer is: **nothing here is a rate, and nothing
here is a concentration.**

## 5. Flags

1. **This paper does not answer the question the cluster was assembled to answer.** It prints no
   concentration of any volatile, in any material, in any unit. It has no chromatography and no
   sensory attribute related to aroma; the eight-person panel judges nothing about smell. **The
   Programme 7 levels table gets no row from Part I.** Recorded here so that this PDF is not
   re-opened for that purpose.
2. **The extrusion is described by set points, not by what the material experienced.** Feed rate,
   screw speed, geometry and five zone temperatures are given; **residence time, die pressure, melt
   temperature, torque and specific mechanical energy are not**. A thermal history cannot be
   reconstructed, so this cannot become a benchmark cook even for the thiamine chemistry.
3. **Two ingredients are unidentified.** The "vegetable oil" is never named — not soy, not
   sunflower, not rapeseed — and it is added *after* hydration and therefore never sees the
   extruder; and the water used for hydration is not characterised. For a lipid-oxidation model the
   oil's identity is the whole question, and it is missing.
4. **The design confounds moisture with temperature.** 30 %/180 C, 34 %/160 C and 38 %/140 C move
   both variables together in opposite directions. Nothing in this paper can separate a moisture
   effect from a temperature effect; the authors treat the three as a single "severity" ladder and
   so must anyone reading them.
5. **Most of the quantitative results live in figures.** Tables 1-4 carry salivary flux, mastication
   time, bolus moisture and a significance matrix. **Every texture value and every particle-geometry
   value is in a figure**, and reaches this dossier only through the numbers the authors chose to
   quote in the prose. Several are explicitly "figure not shown" or "graph not shown". Supplementary
   Table 1 is not on disk.
6. **The soy protein concentrate is characterised by one number.** ">= 70 g/100 g protein (dry
   basis)" is a supplier minimum, not a measurement. **No lipid content is printed anywhere** — and
   the lipid is what a carried-volatile budget would be built on. No moisture as received, no ash,
   no sugar, no free amino acid profile, no lipoxygenase status.
7. **What to request.** (a) **Milani, Menis-Henrique & Conti (2022), "Thiamine as a new ingredient
   for obtaining textured soy protein with meat odour"** — the paper that actually measured the
   volatiles from this exact system and fixed the 34 % M / 160 C condition as optimal; **not on
   disk**, and it is the acquisition this cluster most needs from the Conti group. (b) **Milani &
   Conti (2024), "Textured soy protein with meat odor as an ..."** — the follow-up on protein and
   emulsion stability, also cited here and also **not on disk**. (c) Supplementary Table 1 and the
   underlying data behind Figures 1 and 4 ("data will be made available on request"). (d) A
   composition sheet for Arcon SM, in particular its lipid content. (e) The identity of the
   vegetable oil.
8. **Registry gaps against `data/keys/compounds.yml`.** The only species this paper handles that
   the registry touches at all is **thiamine**, and the registry has no molecule for it — it carries
   `thiamine_availability`, a modifier id, not a compound. Sodium chloride and monosodium glutamate
   are not registry species either (the registry has `gmp` and `imp` as nucleotides but no
   glutamate). **Nothing this paper measures maps to a compound id**, because it measures no
   compounds.
9. **Registry gaps against `data/species/off_flavour_targets.yml`.** None of that file's six
   compounds (hexanal, nonanal, 1-octen-3-ol, 2-pentylfuran, 1-hexanol, furfural) appears in this
   paper. The file is silent on textured vegetable protein as a matrix, and this paper gives no
   basis to change that.
