# Thong et al. 2024 — EXTRACTION (commercial plant-based burgers vs beef: descriptive sensory plus TCATA GC-O)

**Source on disk:** `data/articles/Thong2024.pdf` (the publisher's PDF, 12 pages, article 113848). Read and
checked by eye on 2026-10-09: every number below was read from the page image and cross-checked against
`pdftotext -layout`. Figures 2, 3 and 5 were re-rendered at 200-220 dpi to read labels. The supplementary
tables (S1 standards list, S2 ingredients, S3 sensory means, S4 tabulated volatiles) are **not on disk**.
Written for the plant-based meat-aroma question, not for the core fit.

| field | value |
|---|---|
| Title | "Comparison of differences in sensory, volatile odour-activity and volatile profile of commercial plant-based meats" |
| Authors | Aaron Thong, Vicki Wei Kee Tan, Geraldine Chan, Michelle Jie Ying Choy (SIFBI, A*STAR, Singapore), Ciarán G. Forde (Wageningen University) |
| Venue | Food Research International 177 (2024) 113848 |
| DOI | 10.1016/j.foodres.2023.113848 |

## 1. Methods

- **Products (Table 1, p. 2).** 11 burger products bought in Singapore, all in the sensory test; **6 minced
  products in GC-MS/O** (marked ×). Protein / fat (saturated) in g/100 g as printed:

  | code | main protein, key ingredients | protein | fat (sat.) | GC-MS/O |
  |---|---|---|---|---|
  | AM | beef (mince) | 18 | 13 (4) | × |
  | AP | beef, water, egg white, seasoning (patty) | 15 | 16 (6) | |
  | SM-1 | soy protein concentrate, coconut oil, sunflower oil, flavourings | 16.8 | 11.5 (5.3) | × |
  | SM-2 | soy protein, vegetable oil, pea protein, natural flavouring | 12.2 | 13.9 (1.2) | × |
  | SP | soy protein, pea protein, vegetable oil, seasoning | 14.1 | 18.8 (3.8) | |
  | PM | pea protein, pressed canola oil, refined coconut oil, rice protein, natural flavouring | 17.7 | 12.3 (4) | × |
  | PP | as PM (patty) | 17.7 | 12.4 (4.4) | |
  | MM | mycoprotein, egg white, wheat flour, vegetable oil, maize flour, wheat starch, textured wheat protein, natural flavouring | 14.5 | 2.0 (0.5) | × |
  | MP | mycoprotein, egg white, textured wheat protein, vegetable oil, flavouring | 16 | 8.1 (3.4) | |
  | GM | green spelt, oat flakes, spelt flakes, sunflower seeds, seasoning | 15 | 9 (1) | × |
  | GP | mushrooms, bulgur wheat, wheat gluten, sunflower oil, seasoning | 8.5 | 7.0 (0.7) | |

- **Cooking (§2.1.2, p. 3).** Mince shaped into 20 g (±0.5 g) patties, 1 cm (±0.1 cm) thick, cooked at
  180-200 °C for 3 min per side; frozen ready-to-cook patties at 180-200 °C for 4 min per side. Internal
  temperature at the end of cooking was 70-75 °C, and samples were served warm at 58 ± 2 °C. The cooking
  surface or appliance is not named.
- **Descriptive sensory (§2.1, Table 2, p. 3).** n = 21 (10 males), mean age 27.3 (±5.4), semi-trained in
  two 1-hour sessions. Attributes were scored on a 0-100 VAS, in triplicate across 3 sessions:
  - odour, before consuming: meaty, legume, off
  - taste/flavour: meaty, legume, salty, savoury, off
  - texture/mouthfeel: juiciness, chewiness, oily mouthfeel, flavour aftertaste
  - Statistics: linear mixed model (sample fixed, subject random), ANOVA with post-hoc Bonferroni at
    p = 0.05, and PCA (§2.3.1, p. 4).
- **Volatiles (§2.2.4-2.2.5, p. 4).**
  - Extraction: 100 g of cooked sample with 800 mL water, hydro-distilled at a simmer for 1 h; 100 mL of
    distillate collected. A PDMS stir bar (SBSE) was stirred in 2 mL of distillate for 30 min at 300 rpm.
  - GC: thermal desorption at 250 °C onto a DB-WAX UI column. The effluent was split equally between a
    Q-TOF MS (EI 70 eV, scan m/z 55-400) and an ODP 3 sniff port.
  - Annotation (§2.3.3, p. 5): ADAP-GC deconvolution, then matching against "an in-house library based on
    mass spectral and retention index similarity". Unknowns were "putatively identified" against NIST17.
    HCA used pareto-scaled peak areas. **Concentrations were not measured**: there is no internal-standard
    quantitation and no OAVs.
- **TCATA GC-O (§2.2.2, 2.2.6, 2.3.2, pp. 3-5).**
  - Panel: 27 recruited; **12 completed (6 males)**, mean age 31.4 (±5.2).
  - Descriptors and reference aromas (Table 3, p. 4): meaty (beef extract), legume (soybean extract), fatty
    ((E,E)-2,4-decadienal), nutty (2,5-dimethylpyrazine), sulfurous (furfuryl thiol), other.
  - Procedure: 6 s fading time; one sample per panellist per session, 72 sessions in all (12 × 6),
    complete block design.
  - Data treatment: recognition events binned in 0.05 min (3 s) segments; exponential smoothing with
    α = 0.8; noise filter s_t > 2, meaning at least 3 panellists (25 %) citing the same descriptor at once.
- **MFA (§2.3.4).** Sensory attributes and GC-O citation proportions combined (FactoMineR).

## 2. Findings that matter

### 2a. Panel training check (Table 4, p. 6)

| standard | descriptor | threshold as printed | recognised before | after |
|---|---|---|---|---|
| 2-methylpyrazine | nutty | 60,000 ppb | 27.3 % | 9.1 % |
| 2,5-dimethylpyrazine | nutty | 800 ppb | 27.3 % | 63.6 % |
| dimethyl trisulfide | sulfurous | 0.01 ppb | 54.5 % | 63.6 % |
| furfuryl thiol | sulfurous, meaty | 0.005 ppb | 72.7 % | 81.8 % |
| methional | sulfurous, legume | 0.2 ppb | 45.5 % | 54.5 % |
| (E,E)-2,4-decadienal | fatty | 0.07 ppb | 54.5 % | 63.6 % |
| 2,5-dimethyl-4-hydroxy-3(2H)-furanone | others | 4 ppb | 81.8 % | 100 % |
| bis(2-methyl-3-furyl) disulfide | meaty, sulfurous | **0.0007 ppt** (as printed) | 63.6 % | 81.8 % |
| 3-methyl indole | meaty, others | 0.05 ppb | 81.8 % | 81.8 % |
| total correct recognitions | | | 223 | 203 |
| total noise detections | | | 2910 | 1801 |

- Mean panellist accuracy rose from 56.6 ± 20.2 % to 66.7 ± 14.1 %; paired t-test, one-tailed P < 0.05
  (§3.2.1).
- The percentages are multiples of 1/11 (e.g. 7/11 = 63.6 %, derived here), so this check appears to rest
  on 11 panellists, not the 12 of the final panel. The paper does not say.
- The disulfide's "0.0007 ppt" is probably a unit slip. Christlbauer 2011 gives 0.0008 µg/L, i.e.
  0.8 ppt, in water. The threshold medium is not stated for any Table 4 row.
- "Total compound recognitions" fell (223 to 203) while the text reports improved recognition. The paper
  does not reconcile the two.

### 2b. Descriptive sensory (§3.1, Fig. 1, p. 5)

PCA: PC1 70.41 %, PC2 20.15 %. Beef (AM, AP) loads on meaty odour and meaty flavour; GM and GP (grain and
mushroom) load on legume and off attributes. GM and GP had significantly higher legume odour and flavour
than the other plant samples (p < 0.001), and off odour and off flavour tracked legume. SM-1 had "comparable
meaty odour and flavour intensity to the animal protein samples" (p. 5-6). The attribute means are in
Supplementary Table 3, which is not on disk. Texture (juiciness, chewiness, oily mouthfeel) was comparable
between several plant products and beef (p. 7).

### 2c. TCATA GC-O (§3.2.2, Fig. 2, p. 6-7)

- **Beef patty (AM): 21 odorants cited**, 9 of them sulfurous and 6 meaty. Each PBMA had "over 30 cited
  odour recognition peaks". Compounds eluting in the first 10 min were cited as sulfurous (SM-1, GM, PM)
  and meaty (SM-2).
- More sulfurous and nutty odorants in SM-1, GM and PM than in beef. More fatty and legume odorants in the
  PBMAs, especially GM, PM and SM-2.
- Citation proportions (CP) printed in the text:

  | odorant (putative) | descriptor | RT | CP |
  |---|---|---|---|
  | 2-ethyl-3,5-dimethylpyrazine | nutty | 9.8 min | SM-1 0.25, GM 0.6, PM 0.26; not perceived in AM |
  | 2,5-dimethyl-3-isoamylpyrazine | nutty | 12.9 min | SM-1 0.25; not perceived in AM |
  | dimethyl trisulfide | sulfurous | — | AM 0.3, SM-1 0.25, PM 0.75 |
  | 5,6-dihydro-2,4,6-trimethyl-4H-1,3,5-dithiazine | sulfurous | — | AM 0.25 |
  | (E,E)-2,4-decadienal | fatty (strongest fatty odour in PBMAs) | — | 0.25 to 0.75 |
  | (E)-2-nonenal, 2-undecanone, (E)-2-decenal | fatty | — | 0.25 to 0.47 |

  The fatty compounds were cited less in AM.
- An unidentified odour-active compound at RI 1280 was cited as meaty plus legume in AM and SM-1, with an
  HRMS fragment of C7H11NO. The authors propose an oxazole (p. 9). Not identified.

### 2d. Volatile fingerprint (Fig. 3 HCA heat map, p. 8; row z-scores, qualitative only)

There are two clusters: {AM, SM-1, PM} and {GM, MM, SM-2} (§3.3). Read from colour:
- **SM-1** is darkest for pyrazine, 2-(methoxymethyl)furan, 3-mercapto-2-butanone, tetrahydroquinoxaline,
  furfuryl mercaptan, 2-acetylthiazole, 2-methylpyrazine, 2,5-dimethylpyrazine, indole, the dithiazine and
  2-methyl-3-furanthiol.
- **SM-2** is darkest for the isoamyl methylpyrazines.
- **MM** is darkest for methional, acetophenone, furfural and 2-undecanol.
- **PM** is darkest for 2-(1-methylvinyl)thiophene, heptanal, 2-pentylfuran, hexanal, octanal, nonanal,
  DMTS and several alkylpyrazines.

The text adds that AM and SM-1 "had the lowest levels of 2,4-decadienal" (p. 10). The z-scores are relative
to the row, so no amounts can be read.

### 2e. Sensory-GC-O links (MFA, Figs. 4-5, pp. 9-10)

MFA Dim 1 explains 38.15 % and Dim 2 22.08 %; Dim 1 splits {AM, SM-1, PM} from {GM, MM, SM-2}. Meaty flavour
and odour intensity correlated positively with citations of sulfurous odorants, "identified as" (abstract:
"putatively identified as"):
- 2-methyl-3-furanthiol, RI 1338
- dimethyl trisulfide, RI 1394
- furfuryl mercaptan, RI 1441
- 2-(1-methylvinyl)thiophene, RI 1447
- methional, RI 1464
- 2-(methoxymethyl)furan, RI 1243

Off flavour and off odour correlated with fatty and legume citations, identified as (E,E)-3,5-octadien-2-one
(RI 1578), 2-undecanol (RI 1717) and (E,E)-2,4-decadienal (RI 1811) (p. 7).

**No correlation coefficients or p-values are printed for these links.** They are read from the
correlation circle. The text's RIs do not appear verbatim among the Fig. 5 labels, and the nearest labels
lie within about 5-15 RI units. Examples: Sulfurous@1395 against DMTS at 1394; Sulfurous@1446 against
FFT/thiophene at 1441/1447; Sulfurous@1249 against 1243; Fatty@1570 against octadienone at 1578. Some of these
vectors point along Dim 2 rather than towards the meaty sensory vectors, e.g. Sulfurous@1395 (read from
graph). The paper does not explain the mapping.

**Identification level.**
- Nine compounds had authentic standards in the panel training: 2-methylpyrazine, 2,5-dimethylpyrazine,
  DMTS, furfuryl thiol, methional, (E,E)-2,4-decadienal, HDMF, bis(2-methyl-3-furyl) disulfide and
  3-methylindole (§2.2.6). The full standards list (Table S1) is not on disk.
- In the main text, every odorant linked to sensory is called "putatively identified" (abstract; p. 9). Per
  compound, the main text never says which were confirmed against a standard by RI and spectrum. That
  information would be in Supplementary Table 4.
- 2-methyl-3-furanthiol, 2-(1-methylvinyl)thiophene, 2-(methoxymethyl)furan, the dithiazine,
  (E,E)-3,5-octadien-2-one and 2-undecanol were not among the training standards. On the evidence on disk
  they are **putative** (MS/RI library or NIST17).

**Authors' reading** (p. 9): odour activity was "independent of protein source", so most Maillard products in
the PBMAs come "from exogenous flavouring material (e.g., thiamine, yeast extracts)" rather than the base
protein. Even the best PBMAs (SM-1, PM) "fell short" of beef in meaty character. Balancing sulfurous against
fatty odorants is "needed".

## 3. What it means for the model and for formulation

**Target profile.** In this study, beef is a sparse profile (21 cited odorants) carried by sulfurous and meaty
notes. The PBMAs are denser (more than 30 peaks each) and over-supplied with fatty and legume notes. The
implied plant-based target has two parts:
- **Raise the sulfurous/meaty thiol set**: 2-methyl-3-furanthiol, furfuryl mercaptan (2-furfurylthiol),
  dimethyl trisulfide, methional and possibly 2-(1-methylvinyl)thiophene. Their citation proportions track
  rated meatiness.
- **Lower the lipid-derived carbonyls**: (E,E)-2,4-decadienal, (E,E)-3,5-octadien-2-one, 2-undecanol,
  (E)-2-nonenal, (E)-2-decenal and 2-undecanone. These track off and legume ratings, and beef had the least
  decadienal.

That the best-scoring plant products (SM-1 soy, PM pea) were richest in the thiols and pyrazines, and that
flavouring rather than protein separated them, points to precursor dosing (cysteine, thiamine, reducing
sugar, yeast extract) as the lever, not protein choice. All of this is correlational (MFA), unquantified,
and partly on putative identities.

**What the kinetic model emits** (grep of `src/kinetic_core/species*.py`, 2026-10-09):
- *Emitted as species:* 2-methyl-3-furanthiol (`MFT`), 2-furfurylthiol (`FFT`), their disulfides (`MFTD`,
  `FFTD`), methional (`MTAL`), methanethiol and dimethyl disulfide, 2-acetylthiazole (`ACTZ`), HDMF (`DMHF`),
  pyrazine / 2-methylpyrazine / 2,5-dimethylpyrazine (`PZ`, `MPZ`, `DMP`), furfural (`FUR`), and in the lipid
  lane hexanal, trans,trans-2,4-decadienal, nonanal and 2-pentylfuran.
- *Not emitted:*
  - dimethyl trisulfide (explicitly refused: no H2S on the methionine chain, `parameters_methionine.py`)
  - 2-(1-methylvinyl)thiophene, 2-(methoxymethyl)furan, 3-mercapto-2-butanone
  - 5,6-dihydro-2,4,6-trimethyl-4H-1,3,5-dithiazine (the paper cites it from 2,4-decadienal + cysteine
    model reactions, p. 9; no lane couples the lipid aldehydes to cysteine)
  - 2-ethyl-3,5-dimethylpyrazine and the isoamyl pyrazines; 2-ethylpyrazine
  - (E,E)-3,5-octadien-2-one, 2-undecanol, 2-undecanone, (E)-2-decenal
  - (E)-2-nonenal (in the live `engine.unrepresented_compounds` of
    `results/validation/core_prediction_uncertainty.json`)
  - 3-methylindole and indole

So the model covers the core of the "raise" side: MFT, FFT, methional and their disulfides. On the
"lower" side it covers only decadienal and hexanal. The paper gives no concentrations to score the model
against: citation proportions are not amounts.

## What it does not give

- Concentrations or OAVs of any compound. GC-O citation proportions are detection-frequency data at an
  undiluted extract.
- Per-compound identification level in the main text (in Supplementary Table 4, not on disk).
- Correlation coefficients or significance for the sensory-odorant links.
- Sensory attribute means (Supplementary Table 3, not on disk).
- Ingredient lists in full, flavouring identity or precursor content of the products (only the Table 1
  headline ingredients). Which "natural flavouring" carried thiamine or yeast extract is not stated.
- Cooking surface. The hydro-distillation (1 h simmer in water) can itself generate or lose volatiles; the
  paper does not test for artefacts.
