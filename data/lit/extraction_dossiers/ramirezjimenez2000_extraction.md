# Ramírez-Jiménez, Guerra-Hernández & García-Villanova 2000 — EXTRACTION (browning indicators in 29 Spanish breads plus one controlled laboratory bake at 190 °C for 0-30 min; HMF by HPLC-UV, furosine by ion-pair HPLC, colour as 100 − L\*; six tables, all clean)

### THE FIRST PAIRED HMF-AND-BROWNING TIME SERIES IN A REAL FOOD ON DISK: Table 6 bakes one commercial dough at 190 °C and reports HMF and the colour index at six times, and HMF rises exponentially at 0.191 per minute (doubling every 3.6 min) while the colour index rises linearly — an ordering the engine can be asked about, but not scored on, because the paper never measures a sugar, an amino acid, a temperature inside the loaf or a water activity.

**Source on disk:** `data/articles/ramirez-jimenez2000.pdf` (6 pp., J. Agric. Food Chem. 2000, 48
(9), 4176-4181).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/ramirez-jimenez2000.txt`, 411 lines). **The text layer is clean throughout:
all six tables came through complete and are re-typed in full below**, together with the running
text. Two typographic artefacts of the ACS two-column scan are noted where they occur: the ACS house
style prints "r² = x" as "r ) x" in this text layer (the closing parenthesis is the equals sign),
and the furosine calibration equations arrive with spaces inserted inside the coefficients
(section 3). Figure 1 (chromatograms of crust and crumb) is an image and carries no number this
dossier needs. There is no supplementary material.

## 0. Identity

| field | value |
|---|---|
| Title | "Browning Indicators in Bread" |
| Authors | A. Ramírez-Jiménez, **E. Guerra-Hernández** (corresponding, ejguerra@platon.ugr.es), B. García-Villanova — Departamento de Nutrición y Bromatología, Facultad de Farmacia, Universidad de Granada, Campus Universitario de Cartuja, 18012 Granada, Spain |
| Venue | **J. Agric. Food Chem. 2000, 48 (9), 4176-4181**; published on the Web 2 August 2000 |
| DOI / article ID | `10.1021/jf9907687` (printed as `JF9907687`) |
| Naming | "browning indicators" = **furosine**, **HMF** and **colour**; colour index = **100 − L\*** (CIE L\*a\*b\*, reflectance, on **lyophilised** samples); furosine = ε-N-(furoylmethyl)-L-lysine, "formed during acid hydrolysis of the Amadori compounds fructosyl-lysine, lactulosyl-lysine, and maltulosyl-lysine" |
| Lineage | HMF method from García-Villanova et al. 1993; furosine method from Guerra & Corzo 1996 with the Resmini 1990 clean-up and the Delgado 1992 mobile phase; the first author's own 1998 pharmacy-degree memoir on toasted sliced bread is cited repeatedly for unpublished results |
| Companions on disk | `hamzalioglu2018_extraction.md` (the only measured HMF sink the corpus accepts), `kocadagli2016foodchem_extraction.md` (glucose/wheat flour, the closest matrix with constants), `goncuoglu2016_extraction.md` (HMF and a measured HMF sink in roasted hazelnut), `nguyen2016_extraction.md` (biscuit baking; the source of wave B20's glycation constants), `berk2021_extraction.md` (CML barriers) |

## 1. Why it matters

**Browning is the trunk's one out-of-sample success and it has never been tested in a food.** The
B1 hold-out is Martins' own melanoidin response: `parameters.py`'s `k_tdg_mel` note says so in
terms — "Martins fitted this k to the very melanoidin response that is held out, so a hold-out
evaluation that USES this k is a reproducibility check, not an out-of-sample test" — and
`kinetic_core_b1_fit_report.json` lists the browning response among
`holdout_species_never_read`, with every fitted term estimated on "the nine non-browning responses".
The browning readout itself is `species.melanoidin_repeat_units`, an A470 divided by an epsilon,
in an aqueous glucose/glycine pot at 80-120 °C.

**This paper is browning in bread.** It is the real-food browning paper of the cluster, and what it
actually offers is narrower and more useful than "browning in a food":

1. **Table 6 is a controlled, single-variable bake.** Six identical commercial doughs from the same
   brand and lot, baked in the authors' own laboratory at **190 °C** for **0, 10, 15, 20, 25 and
   30 min**, with **HMF (mg/kg dry matter) and 100 − L\* at each time**. That is the only isothermal
   time series in the paper and the only one where the composition is held fixed. It is a **paired**
   series: the same samples give the furanic marker and the colour marker, so the *relative* shapes
   of the two are measured on one matrix — exponential for HMF, linear for colour. The engine
   carries both quantities (HMF from wave B7's furanic channel; browning from `k_tdg_mel` and the
   repeat-unit readout), and their relative shapes is a question it can be asked.
2. **It measures furosine, which is an Amadori marker, alongside HMF and colour on the same breads.**
   `furosine` is one of the 75 ids in `data/keys/compounds.yml`, and wave B20's glycation arm
   (`r_glc_lysp` → `FLP`, then `r_flp_cml`, `r_flp_cel`, `r_flp_decay`) is built on exactly the
   compound furosine reports — protein-bound fructosyl-lysine. **This paper's central finding is
   that furosine goes UP and then DOWN while HMF and colour go up monotonically**, which is a direct
   structural statement about the Amadori pool's fate at high thermal load. Verbatim from the
   conclusions: "the indicators HMF and color increase with the length and temperature of baking
   process **whereas furosine decreases when higher intensity is reached**." And in the crust:
   "the higher the baking time the lower was the furosine content. In the crust, the high time
   produces a greater extension of the Maillard reaction and therefore a degradation of furosine."
   B20's `r_flp_decay` ("the dominant loss") is precisely that step, and this paper is an
   independent, real-food, second-laboratory observation of the phenomenon it encodes.
3. **It quantifies the crust/crumb split, which is the reason a bread cannot be modelled as one
   pot.** Table 4: HMF in the crumb is 0.6-2.2 mg/kg and in the crust 18.3-176.1 mg/kg, a factor of
   **24 to 104 (mine)** in the same loaf. Any attempt to use this paper as a benchmark must model
   a spatial gradient the engine does not have.

What this paper does **NOT** give the repository, and it is a long list that decides its usability:
**no sugar measurement of any kind** (not glucose, not fructose, not sucrose, not maltose — the
sugars that make the HMF are never quantified in any bread); **no amino acid or lysine measurement**
(protein is measured by Kjeldahl only, to normalise furosine); **no measured temperature inside any
loaf**; **no water activity**; **no pH**; **no rate constant**; **no activation energy**; **no
reference temperature**; and **no melanoidin measurement** — the colour is a reflectance index on a
freeze-dried crumb-plus-crust mixture, not an absorbance of a solution.

## 2. Methods as they matter to a model

- **Sample sets, five of them.**
  1. **Six common (white) breads A-F**, made by one commercial bakery to one formula: **wheat flour
     50 kg, water 27 kg, baker's yeast 2 kg, NaCl 1 kg, previously fermented dough 5 kg, and
     additives**. Sizes 30 to 1000 g. Fermentation and baking conditions per loaf in Table 1.
  2. **Nine special breads** (white with fruits; two whole white; bran, mixed-flour, oat, soy,
     whole mixed-flour and whole rye), same bakery, "similar" formula except for the cereal
     products.
  3. **Six commercial part-baked doughs**, same brand and lot, vacuum-packed, stored at room
     temperature; label ingredients wheat flour, water, baker's yeast, salt, enzymes, emulsifiers,
     dough conditioner; **weight 112-128 g, moisture 32 %**; label baking instruction 220 °C for
     12-15 min. **These are the Table 6 experiment**, baked instead at **190 °C for 10, 15, 20, 25
     and 30 min** in the authors' laboratory, with the unbaked dough as t = 0.
  4. **Nine commercial toasted sliced breads** (3-10 g each, **moisture 6 %**) and one case of small
     toasted sliced breads (2.5-2.8 g, **moisture 3 %**).
  5. **Four commercial snack breads** (1-3 g, **moisture 5 %**).
- **Moisture and water activity.** Moisture by AOAC 925.10 (gravimetric) and printed **per bread in
  Table 1** (18.6 to 34.3 % for the ovenbaked breads; 3-6 % for the toasted and snack products;
  flour "around 12 %"). **Water activity is never measured or estimated anywhere in this paper.**
- **pH.** **Not measured.** The introduction quotes Kroh 1994 that the Maillard reaction is favoured
  "at a pH of 4-7" and that caramelization needs "temperatures > 120 °C, pH < 3 or pH > 9, and low
  A_w", but no pH of any dough or bread is reported.
- **Temperature.** Only the **oven set-point** is known: 200 to 235 °C for the commercial breads
  (Table 1), 190 °C for the laboratory bake, 220 °C on the dough label. **No temperature was
  measured inside any loaf.** The introduction states, citing Hui 1991, that starch gelatinises and
  proteins denature at an **internal temperature of 60-80 °C**, and that during baking "the water
  content on the surface of the loaf becomes lower than in the middle and this, combined with the
  high temperature, is one of the factors that makes the crust different from the crumb". So the
  paper itself says the sample is not isothermal and not iso-moisture.
- **Colour.** CIE L\*a\*b\* by reflectance spectrophotometer **Elrepho 2000** (Datacolor S.A.,
  Spain), illuminant **D65**, calibrated against a **BaSO₄** standard. **"The samples were
  lyophilized prior to the analysis."** Duplicate samples. Reproducibility tested on commercial snack
  bread D (n = 7): **CV of L\* = 0.30 %**. Only **100 − L\*** is used as the index; a\* and b\* are
  measured but never reported.
- **HMF.** Method of García-Villanova 1993. **0.4 g** ground sample into a 10 mL centrifuge tube +
  **7 mL deionised water**, shaken 1 min, centrifuged 10 min at 5000 rpm; repeated twice more;
  supernatants clarified with **0.5 mL each of Carrez I (15 % w/v potassium ferrocyanide) and
  Carrez II (30 % w/v zinc acetate)**; centrifuged again; made to **25 mL**; 2 mL filtered at
  0.2 µm. Chromatography: Konic 500A, 20 µL loop, **Spherisorb S5 ODS2 (250 mm × 40 mm i.d.)**,
  mobile phase **water-acetonitrile 95:5**, **1 mL/min**, **UV at 284 nm**, HMF eluted in 8 min,
  15 min run. External standard, working range **0.02-0.5 mg/L**, calibration
  **Y = 292.17X − 0.27**, r = 0.9999 (n = 7). Duplicates.
  **Performance, printed:** CV **1.57 %** at low HMF (white bread E) and **2.60 %** at high HMF
  (crust of white bread D), both n = 7; recovery by standard addition (12.6-123.9 mg/kg added to
  the 3.4 mg/kg sample) **92.5-100 %, mean 96.2 %**, with "the highest accuracy ... at values
  < 72 mg/kg" — **and every sample except the crust of white bread D is below 72 mg/kg**, so the
  176.1 mg/kg value is the one outside the validated range.
- **Furosine.** Method of Guerra & Corzo 1996. **150 mg** sample hydrolysed with **4.5 mL of 7.95 M
  HCl at 110 °C for 24 h** in a Pyrex screw-cap vial with PTFE-faced septa, **N₂ bubbled 2 min**
  before sealing; filtered; **0.5 mL** onto a Sep-Pak C18 prewetted with 5 mL methanol and 10 mL
  water, eluted with **3 mL of 3 M HCl**, evaporated under vacuum; taken up in 3 mL of
  water:acetonitrile:formic acid **95:5:0.2**. Chromatography: Perkin-Elmer 250 with Waters 717
  autosampler and diode-array detector 235; **50 µL** injected on reversed-phase C18; mobile phase
  **5 mM sodium heptanesulfonate with 20 % acetonitrile and 0.2 % formic acid**, isocratic,
  **1.2 mL/min**, **UV at 280 nm**. **Calibration was by standard addition into a previously
  hydrolysed wheat flour**, two curves (n = 8 each), both r² = 0.9999 — the printed coefficients
  arrive from the text layer with spaces inside the numbers and are transcribed exactly as printed
  in section 3. Duplicates.
- **Protein.** Kjeldahl, AOAC 920.87 — used only to express furosine per 100 g of protein.
- **Statistics.** SPSS 7.5, correlations only. **No kinetic model, no regression of a rate, no
  activation energy.**
- **Reference temperature of any fitted constant.** **There is none: nothing in this paper is
  fitted as a rate.** The only regressions are the six correlations reported in section 3.

## 3. Tables re-typed

### Table 1. "Description, Flours, and Breadmaking Characteristics of Breads"

Column headings as printed: samples / descripn wt (g) / form / flours / fermentation Tª (°C)/time
(min) / baking Tª (°C)/time (min) / moisture %.

| sample | wt (g) | form | flours | fermentation °C / min | baking °C / min | moisture % |
|---|---:|---|---|---|---|---:|
| **white bread A** | 200 | stick | baking flour | 30-35 / 50 | 210 / 30 | 28.6 |
| white bread B | 30 | roll | baking flour | 30-35 / 30 | 235 / 16 | 27.2 |
| white bread C | 250 | (blank) | baking flour | **double ferment.** 30-35 / 15-20 then 30-35 / 30 | 200 / 50 | 30.4 |
| white bread D | 1000 | large round | baking flour | 30-35 / 50 | 200 / 60 | 33.2 |
| white bread E | 200 | stick | baking flour | 30-35 / 50 | 210 / 30 | 30.8 |
| white bread F | 700 | stick | baking flour | 30-35 / 60 | 200 / 50 | 31.8 |
| **white bread with fruits A** | 200 | stick | baking flour with fruits; orange peel; **glucose syrup; sucrose** | 30-35 / 45-60 | 220 / 15-18 | 18.6 |
| **whole white bread A** | 200 | stick | whole white flour-baking flour 1:2 | 30-35 / 50 | 210 / 35 | 30.2 |
| whole white bread B | 500 | large round | whole white flour-baking flour 1:2 | 30-35 / 50 | 210 / 50 | 32.9 |
| **bran bread** | 200 | stick | bran flour (wheat bran; wheat, soy and malt flours; germ wheat; whey)-baking flour 1:1 | 30-35 / 20-25 | 200 / 30 | 31.1 |
| mixed-flour bread | 200 | stick | mixed flour (wheat, corn, sesame, flax, oat, barley, millet whole soy and whole rye flours)-baking flour 1:1 | 30-35 / 20-25 | 200 / 30 | 34.3 |
| oat bread | 200 | stick | oat flour (oat fiber, wheat and oat flakes, wheat flour)-baking flour 1:1 | 30-35 / 20-25 | 200 / 30 | 24.5 |
| soy bread | 200 | stick | soy flour (soy granulated whole, wheat and rye flours)-baking flour 1:1 | 30-35 / 20-25 | 200 / 30 | 27.4 |
| whole mixed-flour bread | 200 | stick | whole mixed flour (wheat, corn, sesame, flax, oat, barley, millet, whole soy and whole rye flours, wheat and oat flakes, oat fibers and granulated soy)-baking flour 1:1 | 30-35 / 20-25 | 200 / 30 | 23.2 |
| whole rye bread | 200 | stick | whole rye flour (whole rye and malt flours, soy lecithin, citric acid, dairy solids and garrafin and guar gums)-baking flour 1:4 | 30-35 / 50 | 210 / 50 | 32.9 |

### Table 2. "Browning Indicators in Different Common Breads"

Footnotes as printed: ᵃ **mg/kg of dry matter**; ᵇ **mg/100 g of protein**.

| sample | HMF ᵃ | 100 − L\* | furosine ᵇ |
|---|---:|---:|---:|
| white bread A | 15.7 | 17.9 | 146.3 |
| white bread B | 21.8 | 17.0 | 125.4 |
| white bread C | **68.8** | **22.8** | 141.4 |
| white bread D | 40.1 | 18.4 | 165.4 |
| white bread E | **3.4** | 18.1 | 177.8 |
| white bread F | 11.8 | **15.9** | **208.1** |

### Table 3. "HMF and 100 − L\* Values in Special Breads"

Footnote ᵃ: **mg/kg dry matter**.

| sample | HMF ᵃ | 100 − L\* |
|---|---:|---:|
| white bread with fruits | **51.3** | **38.2** |
| whole white bread A | 7.4 | 20.1 |
| whole white bread B | 23.2 | 21.5 |
| bran bread | 23.4 | 27.8 |
| mixed-flour bread | 21.1 | 26.1 |
| oat bread | **4.8** | 21.4 |
| soy bread | 18.3 | 25.7 |
| whole mixed-flour bread | 8.7 | 23.8 |
| whole rye bread | 23.4 | 27.7 |

### Table 4. "HMF and Furosine in Crumb and Crust of Breads"

Footnotes: ᵃ **mg/kg of dry matter**; ᵇ **mg/100 g of protein**.

| sample | HMF ᵃ | furosine ᵇ |
|---|---:|---:|
| **white bread A** — crumb | 0.9 | 55.4 |
| white bread A — crust | 21.4 | 125.0 |
| **white bread D** — crumb | 1.7 | 42.8 |
| white bread D — crust | **176.1** | **75.5** |
| **whole white bread A** — crumb | 0.6 | 78.6 |
| whole white bread A — crust | 18.3 | **220.8** |
| **whole white bread B** — crumb | 2.2 | 95.0 |
| whole white bread B — crust | 73.3 | 170.3 |

### Table 5. "HMF and Color in Commercial Breads"

Footnote ᵃ: **mg/kg of dry matter**.

| sample | HMF ᵃ | 100 − L\* |
|---|---:|---:|
| sliced toasted bread A | 11.8 | 24.1 |
| sliced toasted bread B | 13.0 | 27.9 |
| sliced toasted bread C | 17.3 | 23.5 |
| sliced toasted bread D | **87.7** | 27.8 |
| sliced toasted bread E | 38.0 | 27.2 |
| sliced toasted bread F | 16.2 | 25.8 |
| sliced toasted bread G | 21.1 | **28.5** |
| sliced toasted bread H | 47.2 | 24.9 |
| sliced toasted bread I | 24.6 | 22.1 |
| snack A | 4.4 | 18.5 |
| snack B | **2.2** | 24.5 |
| snack C | 10.0 | **16.8** |
| snack D | 6.8 | 17.9 |

### ★ Table 6. "Behavior of Browning Indicators during Baking Time" — THE ONE CONTROLLED EXPERIMENT

Footnote ᵃ: **mg/kg of dry matter**. Conditions from the Methods: six commercial part-baked doughs
from the same brand and lot (weight 112-128 g, **moisture 32 %**), baked in the authors' laboratory
at **190 °C**; t = 0 is the unbaked dough.

| time (min) | HMF ᵃ (mg/kg dm) | 100 − L\* |
|---:|---:|---:|
| 0 | 0.06 | 13.7 |
| 10 | 0.47 | 14.1 |
| 15 | 1.27 | 16.3 |
| 20 | 3.07 | 17.2 |
| 25 | 7.38 | 19.6 |
| 30 | 19.58 | 20.6 |

### Every correlation printed in the paper

The ACS text layer renders "r² =" as "r )". The values as printed:

| pair | sample set | form | r² as printed |
|---|---|---|---|
| HMF vs 100 − L\* | common breads (Table 2) | linear | **0.7331** |
| HMF vs 100 − L\* | special breads (Table 3) | linear | **0.8347** |
| HMF vs 100 − L\* | commercial toasted breads (Table 5) | linear | **0.08** ("not significant") |
| HMF vs 100 − L\* | snack breads (Table 5) | linear | **0.7142** |
| HMF vs 100 − L\* | the laboratory bake (Table 6) | linear | **0.7391** |
| HMF vs 100 − L\* | the laboratory bake | exponential | **0.9080** |
| 100 − L\* vs baking time | the laboratory bake | linear | **0.9357** |
| 100 − L\* vs baking time | the laboratory bake | exponential | **0.9343** |
| HMF vs baking time | the laboratory bake | linear | **0.6523** |
| **HMF vs baking time** | **the laboratory bake** | **exponential** | **0.9988** |
| furosine vs HMF | common breads | linear, **r** (not r²) | **−0.4298** |
| furosine vs baking temperature | common breads | linear, **r** | **−0.6016** |

### Other numbers printed in the running text

| quantity | value | where |
|---|---|---|
| overall HMF range across all breads | **2.2 to 68.8 mg/kg** | Abstract |
| overall colour index range | **17.0 to 38.2** | Abstract |
| furosine range, common bread | **125 to 208 mg/100 g protein** | Abstract and Table 2 |
| crumb/crust HMF range | "0.9-1.76 mg/kg" — **an abstract typo; Table 4 shows 0.6-176.1 mg/kg** (Flags 2) | Abstract |
| crumb/crust furosine range | **43 to 221 mg/100 g protein** | Abstract |
| colour index of the **flours** | **5.4 to 11.4**; whole flours mean **8.8**; baking flour lowest, whole mixed and bran flours highest; whole rye flour **8**; flour moisture **around 12 %** | Results, Color |
| **HMF in the flours** | **"No HMF was detected in the flours."** | Results, HMF |
| mean 100 − L\*, snack breads | **19.4** | Results, Color |
| mean 100 − L\*, sliced toasted breads | **25.7** | Results, Color |
| mean HMF, snack breads | **5.8 mg/kg** | Results, HMF |
| toasting increment, sample H | **15 mg/kg before toasting → 47.2 mg/kg after** | Results, HMF |
| five toasted slices, same lot and same colour | **29.3, 36.0, 36.2, 41.9, 45.1 mg/kg HMF** | Results, HMF |
| crumb furosine, common breads | ~**50 mg/100 g protein** | Results, Furosine |
| crumb furosine, special breads | **78.6 and 95 mg/100 g protein** | Results, Furosine |
| crust furosine range | **75-221 mg/100 g protein** | Results, Furosine |
| literature furosine span quoted | **1 mg/100 g protein in wheat grain to 3000 mg/100 g protein in baby cereals** | Discussion |
| furosine in toasting (author's own 1998 memoir) | begins to fall **after 10 min** of toasting | Discussion |
| HMF calibration | Y = 292.17X − 0.27, r = 0.9999, range 0.02-0.5 mg/L, n = 7 | Methods |
| furosine calibration curves (**as printed, with the text layer's internal spaces**) | "Y ) 9 756 87 8.61X − 36 197.418 (range 0.0383-0.383 µg) r² ) 0.9999" and "Y ) 9 58 4 64 3.42 X + 15 901.2 (range 0.0193-0.0958 µg) r² ) 0.9999" — i.e. slopes of about **9.76 × 10⁶** and **9.58 × 10⁶** area units per µg (**my reading of the spacing; do not use the coefficients**) | Methods, Flags 3 |

### Arithmetic on the printed tables (all mine)

**1. HMF accumulates exponentially at 0.191 min⁻¹ at 190 °C, and the paper's own r² reproduces
exactly.** A least-squares line through ln(HMF) against time on all six points of Table 6 gives
**slope 0.1913 min⁻¹, intercept −2.727, R² = 0.9988** — the identical R² the paper prints for its
exponential fit, so my regression is the same one. **Doubling time 3.62 min.** A linear fit of HMF
against time gives **R² = 0.6523**, again exactly the paper's printed value. **This is the single
most transportable number in the paper: an apparent first-order accumulation constant for HMF in a
baking bread at an oven temperature of 190 °C.** It is an *accumulation* constant, not the rate
constant of any step — it is the net of every source and every sink, in a matrix whose water content
and internal temperature are both changing — and it is classed accordingly in section 4.

**2. The colour index rises linearly, and here my arithmetic and the paper's differ slightly.**
A least-squares line of 100 − L\* against time gives **slope 0.250 per minute, intercept 12.75,
R² = 0.9228 (mine)**, against the paper's printed **0.9357**. The difference, 0.013, is small but
real, and I cannot reproduce their value from the six printed points by ordinary least squares.
Recorded, not smoothed (Flags 4). What is not in doubt is the **shape contrast**: HMF is
exponential (R² 0.9988 exponential vs 0.6523 linear) and the colour index is linear (R² 0.92-0.94
linear, and the paper's own exponential fit is no better, 0.9343). Over 30 min at 190 °C **HMF rises
326-fold while the colour index rises 1.50-fold (mine)**.

**3. The crust/crumb factors.** HMF crust ÷ crumb: white bread A **23.8**, white bread D **103.6**,
whole white A **30.5**, whole white B **33.3** (all mine). Furosine crust ÷ crumb: **2.26, 1.76,
2.81, 1.79** (mine). **The two markers scale completely differently across the same gradient** —
HMF by one to two orders of magnitude, furosine by less than threefold — and in the most severely
baked loaf (white bread D, 1000 g, 200 °C for 60 min) the pattern inverts: it has **the highest
crust HMF (176.1) and the lowest crust furosine (75.5)** of the four. That is the furosine burn-off
the conclusions describe, visible in a single table.

**4. The furosine-vs-HMF anticorrelation.** The paper prints **r = −0.4298** for furosine against
HMF and **r = −0.6016** for furosine against baking temperature across the six common breads. As r²
those are **0.185 and 0.362 (mine)** — weak, and on six points, so neither would reach significance;
but both signs are negative, and the crust/crumb table and the conclusions say the same thing more
strongly. **Treat the direction as the finding and the coefficients as decoration.**

**5. Two loaves that isolate the water variable.** White breads A and E are the same flour, weight,
shape, fermentation (30-35 °C, 50 min) and bake (210 °C, 30 min); the only stated difference is that
"the dough of bread E had a smaller water content than that of bread A" — yet the *bread* E is the
wetter one (moisture 30.8 % against A's 28.6 %). Their markers separate hard: **HMF 3.4 vs 15.7
mg/kg (a factor of 4.6, mine)** and **furosine 177.8 vs 146.3 mg/100 g protein**, with colour almost
identical (18.1 vs 17.9). **The less-advanced loaf has more furosine and less HMF**, which is the
authors' own reading and is the cleanest single-pair demonstration of the Amadori-then-HMF ordering
in the paper.

**6. The correlation collapses across brands and holds within one.** HMF against colour: r² = 0.73
(one bakery's common breads), 0.83 (its special breads), 0.71 (snacks), 0.74 (one lot of dough,
laboratory-baked) — but **0.08 across nine commercial toasted breads from nine different brands**.
And five slices from the *same* lot at the *same* colour span **29.3 to 45.1 mg/kg HMF, a factor of
1.54 (mine)**. The paper's conclusion: "HMF is an indicator more sensible than color, being useful
in samples without distinct color." **For this repository the lesson is the reverse of the usual
one: colour and HMF are not interchangeable readouts of Maillard extent, and a model validated on
one of them is not thereby validated on the other.**

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Both of this paper's chemical markers are keyed:
**`hmf`** and **`furosine`** (and the family id `furosine_cml`). The colour index has no registry id
and no counterpart in the engine's units. Glucose, fructose, sucrose, maltose, lysine and every
amino acid remain absent — which is moot here, because the paper measures none of them.

Every row below shares: **wheat-flour breads or commercial part-baked doughs, oven-baked in air at a
set-point of 190 to 235 °C; moisture 18.6-34.3 % (ovenbaked) or 3-6 % (toasted/snack); no measured
internal temperature; no water activity; no pH; no sugar or amino acid measurement; samples
lyophilised before colour reading; HMF and furosine in duplicate.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| **net HMF accumulation** | **apparent first-order accumulation constant, ln(HMF) vs t** | **0.1913** (R² 0.9988; the paper prints the same R² for its own exponential fit) | **min⁻¹** | commercial dough, 32 % initial moisture, **oven 190 °C**, 0-30 min, one lot | **first order in time, not in any reactant** — it is a net accumulation, not a step | derived from Table 6 (**mine**); the paper prints the R², not the constant | **derived_assumption** — an accumulation rate, NOT a rate constant of any reaction; carries no reactant order and no reference temperature |
| " | doubling time of HMF | **3.62** | min | as above | — | derived (**mine**) | derived_assumption |
| " | HMF at 0 / 10 / 15 / 20 / 25 / 30 min | **0.06 / 0.47 / 1.27 / 3.07 / 7.38 / 19.58** | **mg/kg dry matter** | as above | — | **Table 6** | **level_only** — the only isothermal HMF time series in a real food on disk |
| **net browning** | 100 − L\* at 0 / 10 / 15 / 20 / 25 / 30 min | **13.7 / 14.1 / 16.3 / 17.2 / 19.6 / 20.6** | dimensionless colour index (CIE L\*, D65, on lyophilised sample) | as above | — | **Table 6** | **level_only** — and see the unit warning below |
| " | linear rise of the colour index | **0.250 per min**, intercept 12.75, R² 0.9228 (**mine**; the paper prints R² 0.9357) | index units / min | as above | — | derived from Table 6 (**mine**) | derived_assumption (and Flags 4) |
| shape contrast | HMF exponential vs colour linear | HMF: R² 0.9988 exponential against 0.6523 linear. Colour: 0.9357 linear against 0.9343 exponential (paper's own values) | — | as above | — | Results, laboratory bake | **within_study_ratio — the single most useful structural claim in the paper** |
| " | fold change over 30 min | HMF **×326**, colour index **×1.50** (both mine) | — | as above | — | derived (**mine**) | within_study_ratio |
| HMF vs colour | linear r² by sample set | 0.7331 (common) / 0.8347 (special) / **0.08** (nine brands of toasted) / 0.7142 (snacks) / 0.7391 (one lot, laboratory-baked) | — | see each table | — | Results | within_study_ratio |
| HMF | levels in 29 commercial breads | **2.2 to 87.7** | mg/kg dry matter | see Tables 2, 3, 5 | — | Tables 2, 3, 5 | **level_only** |
| HMF | **crumb** | **0.6 / 0.9 / 1.7 / 2.2** | mg/kg dry matter | four loaves, Table 4 | — | Table 4 | **level_only** |
| HMF | **crust** | **18.3 / 21.4 / 73.3 / 176.1** | mg/kg dry matter | the same four loaves | — | Table 4 | **level_only** (the 176.1 is outside the validated accuracy range, Flags 5) |
| " | crust ÷ crumb HMF | **23.8 / 103.6 / 30.5 / 33.3** | — | " | — | derived (**mine**) | within_study_ratio |
| **furosine (the Amadori marker)** | common breads | **125.4 to 208.1** | **mg/100 g of protein** | six white breads, 200-235 °C | — | Table 2 | **level_only** |
| " | crumb | **42.8 / 55.4 / 78.6 / 95.0** | mg/100 g protein | four loaves | — | Table 4 | **level_only** |
| " | crust | **75.5 / 125.0 / 170.3 / 220.8** | mg/100 g protein | the same four | — | Table 4 | **level_only** |
| " | crust ÷ crumb furosine | **1.76 / 2.26 / 1.79 / 2.81** | — | " | — | derived (**mine**) | within_study_ratio |
| **the Amadori maximum** | furosine RISES then FALLS with thermal load while HMF and colour rise monotonically | direction only | — | all sample sets, plus the author's 1998 memoir (falls after 10 min of toasting) | — | Conclusions and Results, Furosine | **within_study_ratio — the structural claim that bears on wave B20's `r_flp_decay`** |
| " | furosine vs HMF, and vs baking temperature | **r = −0.4298** and **r = −0.6016** (r, not r²; r² = 0.185 and 0.362, mine) | — | six common breads | — | Results, Furosine | within_study_ratio (weak, n = 6) |
| flours | **HMF not detected** | — | — | all breadmaking flours | — | Results, HMF | **measured null** — the model may set HMF(t=0) = 0 for a bread |
| flours | colour index | 5.4 to 11.4; whole flours mean 8.8; moisture ~12 % | — | — | — | Results, Color | level_only (the browning baseline of the matrix) |
| method performance | HMF CV, recovery | CV 1.57 % (low) and 2.60 % (high), n = 7 each; recovery 92.5-100 %, mean **96.2 %** | % | — | — | Results, HMF | measured (analytical) |
| method performance | L\* CV | **0.30 %**, n = 7 | % | — | — | Results, Color | measured (analytical) |
| — | **any rate constant, any barrier, any reference temperature** | **none** | — | — | — | — | **absent** |
| — | any sugar, amino acid, lysine, dicarbonyl, melanoidin, pH, water activity or internal temperature | **none measured** | — | — | — | — | **absent** |
| — | chromatograms of crust and crumb | — | — | — | — | Figure 1 | **figure_only** |

### Can this be set against the trunk's browning, and on what basis?

**The honest answer is: as an ordering, yes; as a number, no — and the obstacle is the readout, not
the chemistry.**

**(a) The browning unit does not convert.** The engine's browning is
`species.melanoidin_repeat_units` — **mmol of melanoidin repeat units per litre**, obtained as
A470 ÷ ε, where a repeat unit is *one 3-deoxyglucosone plus one glycine* (`MELANOIDIN_REPEAT_UNIT_CARBON = 8`,
`..._NITROGEN = 1`), which is the only definition Martins' step 9 supplies. This paper's browning is
**100 − L\***, a **reflectance** index measured on a **lyophilised solid** with a BaSO₄ white
reference under D65. Those are not the same physical quantity, they are not proportional to each
other in general, and no calibration between them exists in this paper or anywhere in the corpus.
**A 100 − L\* of 20.6 cannot be compared with a mmol/L of melanoidin.** What *can* be compared is
the shape against time and the ordering against treatment, and that is what section 4's table
offers.

**(b) The HMF unit converts, but the pot does not.** HMF in mg/kg of dry matter is a mass fraction,
and the engine answers in µg/L or mmol/L of a liquid pot. The conversion would need the bread's
water content at each time — which changes during the bake and is measured only for the finished
loaf — and, worse, it would need the **initial sugar and amino acid loading**, which this paper
never measures for any sample. **Without a sugar, there is no charge to start the model from.** So
Table 6 cannot become a benchmark row, however good it is.

**(c) What it CAN be used for, and this is not small.** Three claims, all of them structural, all of
them from a real food, and all of them testable against the engine as it stands:

1. **HMF rises exponentially and browning rises linearly over the same 30 minutes in the same
   sample.** Any model run on a bread-like charge that produces a linear HMF and a saturating
   browning has the two shapes the wrong way round. This is a paired, single-lot, single-oven
   observation, and it is the strongest thing this paper has.
2. **The Amadori marker turns over.** Furosine rises with thermal load and then falls, most sharply
   in the crust of the most severely baked loaf. Wave B20's `r_flp_decay` — "bound fructosyl-lysine
   → 3-deoxyglucosone + bound lysine ... the dominant loss" — predicts exactly that turnover, and
   this is an independent, second-laboratory, real-food observation of it. It is the first
   corroboration of B20's topology from outside the Nguyen/Berk lane.
3. **Colour and HMF decouple across matrices.** Within one bakery's product line they correlate at
   r² 0.71-0.83; across nine brands they correlate at r² 0.08; and five slices of the same lot at
   the same colour differ 1.5-fold in HMF. **A browning hold-out that passes says nothing about
   whether the HMF answer is right, and vice versa** — which is a caveat the trunk's one
   out-of-sample success should carry, since browning is the response it succeeds on and HMF is a
   response the engine also emits.

**(d) What CANNOT be transported.** No rate, no barrier, no reference temperature, no pot. The oven
set-points (190-235 °C) are not the reaction temperature; the paper says so itself when it puts the
internal gelatinisation temperature at 60-80 °C and describes the crust/crumb gradient. Reading the
190 °C of Table 6 as a reaction temperature would be reading an oven wall as a chemistry.

## 5. Flags

1. **This is a survey with one experiment inside it.** Twenty-nine of the thirty samples are
   commercial products differing in flour, size, shape, ingredients, fermentation and bake all at
   once; only Table 6 varies one thing (time) on one material. **Only Table 6 should ever be used
   quantitatively.** Everything in Tables 2, 3 and 5 is a level in an uncontrolled matrix.
2. **The abstract's crumb/crust HMF range is a typo.** It reads "Levels of HMF had a wide range
   (0.9-1.76 mg/kg)". Table 4 shows crumb 0.6-2.2 and crust 18.3-**176.1** mg/kg. The abstract has
   dropped a factor of 100 on the top end and taken 0.9 rather than 0.6 at the bottom. **Cite
   Table 4, never the abstract.**
3. **The furosine calibration coefficients are unusable as printed.** The text layer gives
   "Y ) 9 756 87 8.61X − 36 197.418" and "Y ) 9 58 4 64 3.42 X + 15 901.2", with spaces inside the
   numbers. My reading is 9 756 878.61 and 9 584 643.42 area units per µg, which is self-consistent
   (the two curves agree to 1.8 %) — **but it is a reading of broken spacing and must not be used as
   a number.** It does not matter for anything in this dossier, since only the resulting furosine
   values are used, and those are printed cleanly.
4. **I cannot reproduce the paper's colour-vs-time r².** Ordinary least squares on the six printed
   points of Table 6 gives R² = 0.9228; the paper prints 0.9357. Every other correlation I checked
   reproduces exactly (HMF exponential 0.9988, HMF linear 0.6523, HMF-vs-colour 0.7391 against the
   printed 0.7391). The one discrepancy is recorded and not resolved.
5. **The largest HMF value in the paper is outside the authors' own validated range.** They state
   "The highest accuracy was obtained at values < 72 mg/kg. All the samples presented values
   < 72 mg/kg except for the crust of white bread D" — which is the **176.1 mg/kg** value, the
   headline of Table 4 and the source of the 103-fold crust/crumb ratio. Carry it with that caveat.
   (The 87.7 mg/kg of toasted bread D in Table 5 is also above 72 and is not mentioned by the
   authors.)
6. **Colour was measured on LYOPHILISED samples.** "The samples were lyophilized prior to the
   analysis." So the index is not the colour of bread as eaten and not the colour of a crust *in
   situ*; it is the colour of a freeze-dried powder of the whole sample, crumb and crust together
   except where Table 4 separates them (and Table 4 reports no colour at all). Any comparison with
   a surface-colour measurement elsewhere in the corpus — the hazelnut paper's computer-vision L\*,
   for one — is comparing two different preparations.
7. **No temperature was measured inside anything.** Every temperature in this paper is an oven
   set-point or a label instruction. The paper itself describes the sample as non-isothermal and
   drying from the surface inward. **The 190 °C of Table 6 is an oven, not a reaction temperature**,
   and the apparent 0.191 min⁻¹ folds the whole heat-up and the whole moisture history into one
   number.
8. **The furosine turnover is inferred from a cross-section, not followed in time.** Table 6, the
   one time series, reports HMF and colour but **not furosine**. The turnover claim rests on
   comparing different breads and on the crust/crumb split, plus the first author's unpublished 1998
   memoir ("furosine levels began to descend after 10 min of the toasting process"). It is a strong
   and consistent direction; it is not a measured maximum on one sample.
9. **Furosine is an acid-hydrolysis proxy, not the Amadori compound.** It is formed *during the
   24 h in 7.95 M HCl at 110 °C*, at a yield that depends on which Amadori compound it came from —
   the paper names three (fructosyl-lysine from glucose, lactulosyl-lysine from lactose,
   maltulosyl-lysine from maltose) — and no conversion factor is given or applied. Wave B20's `FLP`
   is a molar pool of bound fructosyl-lysine; furosine in mg per 100 g of protein is that pool times
   an unknown, matrix-dependent recovery. **Comparable in direction, not in magnitude.**
10. **Bread B was baked at 235 °C, an ingredient-list bread contains added glucose syrup and
    sucrose, and one bread is a 1 kg loaf baked for an hour.** The sample set spans conditions no
    single model charge covers, and the highest-browning sample (white bread with fruits, 100 − L\*
    = 38.2, HMF 51.3 mg/kg) is the one whose sugars were added by the baker and never measured.
11. **What this paper does not contain**: any rate constant; any activation energy; any reference
    temperature; any sugar, amino acid or lysine measurement; any dicarbonyl; any melanoidin; any
    pH; any water activity; any internal temperature; any replicate count beyond "duplicate samples"
    (the n = 7 figures are method-precision runs, not sample replicates); any error bar on any
    entry of any of the six tables; and any supplementary material.
12. **What to request**: (i) the **sugar and free amino acid composition of the commercial dough**
    used for Table 6 — with it, Table 6 becomes a candidate benchmark row and without it it cannot
    be one; (ii) a **furosine time course on the same Table 6 doughs**, which would turn the
    Amadori-turnover claim from a cross-sectional inference into a measured maximum and would be a
    direct test of wave B20's `r_flp_decay`; (iii) the **a\* and b\*** values, measured but never
    reported, and any calibration between 100 − L\* and an absorbance; (iv) an internal temperature
    trace for the laboratory bake.
13. **Registry gaps against `data/keys/compounds.yml`**: `hmf` and `furosine` are both keyed, so
    this paper's two chemical markers are addressable — a rarity in this corpus. **What has no id
    is the colour index itself**: there is no browning or colour key in the registry, and the
    engine's own browning quantity (melanoidin repeat units, mmol/L) has no registry id either. If a
    real-food browning row is ever to be scored, that gap has to be closed first, along with a
    stated, sourced mapping between an absorbance-based melanoidin concentration and a reflectance
    index — which no paper on disk supplies.
