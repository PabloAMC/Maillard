# Thong et al. 2024 — EXTRACTION (commercial plant-based burgers vs beef: descriptive sensory plus TCATA GC-O)

**Source on disk:** `data/articles/Thong2024.pdf` (the publisher's PDF, 12 pages, article 113848). Read and
checked by eye on 2026-10-09: every number below was read from the page image and cross-checked against
`pdftotext -layout`. Figures 2, 3 and 5 were re-rendered at 200-220 dpi to read labels. The supplementary
information is now on disk as `data/articles/Thong2024_Supplementary.docx` and was read in full on
2026-10-09 (section 2x). Its numbering is S1 ingredients, S2 standards, S3 sensory means, S4 tabulated
volatiles. Written for the plant-based meat-aroma question, not for the core fit.

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
Supplementary Table 3, transcribed in section 2x.4 (SM-1 meaty odour 66.49 vs beef mince 72.48, same
letter group). Texture (juiciness, chewiness, oily mouthfeel) was comparable
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

**Identification level (revised 2026-10-09 from the SI, section 2x.2-2x.3).**
- Nine compounds were in the panel-training sniffing mixture: 2-methylpyrazine, 2,5-dimethylpyrazine,
  DMTS, furfuryl thiol, methional, (E,E)-2,4-decadienal, HDMF, bis(2-methyl-3-furyl) disulfide and
  3-methylindole (§2.2.6). The full standards list is Supplementary Table 2 (39 Sigma-Aldrich standards);
  it does not include HDMF or 3-methylindole, and lists the disulfide only as an impurity of the furfuryl
  mercaptan reagent.
- In the main text, every odorant linked to sensory is called "putatively identified" (abstract; p. 9).
  **Supplementary Table 4 confirms this per compound: its legend defines an "AS" (authentic standard) code,
  but no row uses it.** DMTS, furfuryl mercaptan and methional, whose standards were on hand, are coded
  "RI, MS"; so are 2-methyl-3-furanthiol and 2-(1-methylvinyl)thiophene, which had no standard.
- All six sensory-linked sulfurous odorants are therefore **putative (RI + MS)**, none standard-confirmed.
  2-methyl-3-furanthiol has the weakest RI match of them (Expt 1338 vs Lit 1318) and is **not detected in
  the beef**. 2-(methoxymethyl)furan carries odour "ND" in Table S4; the odour-active peak at its RI (1243)
  is 3,4-dimethylthiophene ("Sulfurous, Meaty", GM only).
- (E,E)-3,5-octadien-2-one and 2-undecanol are likewise "RI, MS" only.

**Authors' reading** (p. 9): odour activity was "independent of protein source", so most Maillard products in
the PBMAs come "from exogenous flavouring material (e.g., thiamine, yeast extracts)" rather than the base
protein. Even the best PBMAs (SM-1, PM) "fell short" of beef in meaty character. Balancing sulfurous against
fatty odorants is "needed".

## 2x. Supplementary information (read 2026-10-09)

**Source:** `data/articles/Thong2024_Supplementary.docx` (one Word file: Supplementary Tables 1-4 and
Supplementary Figure 1). Read in full on 2026-10-09: the tables from the document's own text (XML), the one
embedded image (Supp. Fig. 1, `word/media/image1.png`) by eye. The SI numbers its tables differently from
the guess in the Source paragraph above: **S1 = ingredients, S2 = standards, S3 = sensory means, S4 =
volatiles.** Brand names are printed in the SI (Table S1): AM/AP Master Grocer, SM-1 Impossible, SM-2/SP
vEEF, PM/PP Beyond, MM/MP Quorn, GM Seitenbacher, GP Amy's.

### 2x.1 Ingredients (Supplementary Table 1)

Presence of the precursor-relevant ingredients, read from the printed lists (a dash = not listed; nothing
here is a measured content). The right-hand column quotes the listed ingredient words.

| code | thiamine | yeast / yeast extract | heme / leghemoglobin | cysteine | flavouring | other Maillard-relevant items as printed |
|---|---|---|---|---|---|---|
| AM | – | – | – | – | – | "Beef (100%)" |
| AP | – | – | – | – | – (MSG, dextrose, sugar listed) | "Monosodium Glutamate", "Dextrose", "Sugar", "Onion powder", "Garlic Powder" |
| **SM-1** | **"Thiamine Hydrochloride (Vitamin B1)"** | **"Yeast Extract"** | **"Soy Leghemoglobin"** | – | "Natural Flavors" | "Cultured Dextrose", "Zinc Gluconate", B2/B3/B6/B12 |
| SM-2 | – | – | – | – | "Natural Flavours" | "Malt Extract (Barley)", "Beetroot Powder" |
| SP | "Vitamins (B3, B6, B2, B1, B12)" | "Yeast Extracts" | – | – | – | "Malt Powder (Barley)", "Mineral (Iron)", garlic, herbs, spices |
| **PM** | – (B3, B6, B12, pantothenate only) | **"Dried Yeast"** | – | – | "Natural Flavors" | "Beet Powder Colour", "Apple Extract", "Pomegranate Concentrate", "Cocoa Butter", "Zinc Sulfate" |
| PP | – | – | – | – | "Natural Flavors" | "Beet Juice Extract (for color)", "Apple Extract", "Pomegranate Extract", "Mung Bean Protein" |
| MM | "Wheat Flour (... Iron, Niacin, Thiamine)" (flour fortification) | – | – | – | "Natural Flavouring" | "Milk Proteins", "Dextrose", "Tetrasodium Diphosphate, Sodium Carbonate" (raising agents) |
| MP | – | – | – | – | "Flavouring (contains Smoke Flavourings)" | "Roasted Barley Malt Extract", "Colour: Plain Caramel", "Milk Proteins", onions |
| GM | – | "Food Yeast" | – | – | – | vegetables (onions, tomatoes, garlic, leeks, carrots, celery ...) |
| GP | – | – | – | – | – | "Organic Mushrooms", onions, walnuts, garlic |

- **No product lists cysteine** (or any free amino acid other than MSG in AP), and none lists a reaction
  flavour by name. Whatever sulfur precursor is inside "Natural Flavors" is not disclosed.
- **SM-1 (Impossible) is the only GC-analysed product with all three of thiamine HCl, yeast extract and a
  heme protein.** PM (Beyond) has dried yeast but no thiamine; MM (Quorn mince) has thiamine only as wheat
  flour fortification; SM-2 has neither (malt extract instead).

### 2x.2 Standards (Supplementary Table 2)

39 standards, all Sigma-Aldrich, are listed (columns: Standard | Company | CAS No. | Purity). Those that are
Maillard/meat-relevant, with purity as printed:

| Standard | CAS No. | Purity |
|---|---|---|
| 2-Methylpyrazine | 109-08-0 | ≥99% (FCC) (FG) |
| 2,5-Dimethylpyrazine | 123-32-0 | ≥98% (FG) |
| 2-Ethylpyrazine | 13925-00-3 | ≥98% (FG) |
| 2,3,5-Trimethylpyrazine | 14667-55-1 | ≥99% (FCC) (FG) |
| 2,3-Diethyl-5-methylpyrazine | 18138-04 (as printed; truncated) | ≥99% (FG) |
| Pyrazine | 290-37-9 | ≥99% (FG) |
| 5,6,7,8-Tetrahydroquinoxaline | 34413-35-9 | ≥97% (FG) |
| 2-Acetylthiazole | 24295-03-2 | ≥99% (FG) |
| 2-Methylthiophene | 554-14-3 | 98% |
| 3-mercapto-2-butanone | 40789-98-8 | ≥95% (FG) |
| Dimethyl trisulfide | 3658-80-8 | ≥98% (FG) |
| Furfuryl mercaptan | 98-02-2 | 98% (FG) |
| Furfuryl disulfide | 4437-20-1 | ≥95% (FG) |
| Furan, 3,3´-dithiobis[2-methyl- | 28588-75-2 | "Impurity in furfuryl mercaptan" |
| Methional | 3268-49-3 | ≥98% (FG) |
| Phenylacetaldehyde | 122-78-1 | ≥95% (FCC) (FG) |
| Furfuryl acetate | 623-17-6 | ≥98% (FG) |
| Pyrrole | 109-97-7 | ≥98% (FCC) (FG) |
| Indole | 120-72-9 | ≥97% (FG) |
| Hexanal; Heptanal; Octanal; Nonanal; Decanal | 66-25-1; 111-71-7; 124-13-0; 124-19-6; 112-31-2 | ≥95% (FG); ≥95% (FCC) (FG); ≥95% (FCC) (FG); ≥98% (FG); ≥97% (FG) |
| 2-Nonenal, (E)-; 2-Decenal, (E)-; 2,4-Decadienal, (E,E)- | 18829-56-6; 3913-81-3; 25152-84-5 | 97%; ≥95% (FCC) (FG); ≥90% (FG) |
| 2-Pentylfuran | 3777-69-3 | ≥98% (FG) |

The other 11 entries are alcohols (1-hexanol, 1-octanol, 1-octen-3-ol, 1-pentanol), ketones (3-octanone,
acetophenone), benzaldehyde, acids (nonanoic, octanoic), γ-caprolactone and γ-heptalactone.

- **Not in Table S2:** 2-methyl-3-furanthiol, 2-(1-methylvinyl)thiophene, 2-(methoxymethyl)furan, the
  dithiazine, (E,E)-3,5-octadien-2-one, 2-undecanol, 2-undecanone, 2-ethyl-3,5-dimethylpyrazine, the
  isoamyl pyrazines. Also absent: **HDMF and 3-methylindole**, although both are in the main text's
  9-standard sniffing mixture (Table 4). The bis(2-methyl-3-furyl) disulfide "standard" is an impurity in
  the furfuryl mercaptan reagent, not a purchased standard.

### 2x.3 Per-odorant identification (Supplementary Table 4)

Caption: "Volatile compounds **putatively identified** from the mass spectra and their respective peak areas
(total ion counts)". Codes: RI retention-index match, MS mass-spectral match, AM accurate mass of fragment /
molecular ion, **AS authentic standards**. Odour = GC-O descriptor, "ND" if none perceived. Peak areas are
TIC counts (not concentrations), columns AM, SM-1, GM, MM, PM, SM-2.

**Not one row of Table S4 carries the code "AS"**, including the compounds whose standards are listed in
Table S2 (DMTS, furfuryl mercaptan, methional, 2-acetylthiazole, the pyrazines, decadienal, the
aldehydes). Their identification column reads "RI, MS". So by the SI's own coding **no odorant is
confirmed against an authentic standard**; the standards appear to have served the training mixture and
perhaps the in-house RI library, but the SI does not say so. Every identity below is putative (level:
RI + MS, or accurate mass only for "Unk").

Sulfur-containing block (transcribed in full):

| Name | RT / min | Expt RI | Lit RI | Odour | AM | SM-1 | GM | MM | PM | SM-2 | Identification |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2-Methylthiophene | 4.837 | 1107 | 1108 | ND | ND | 67,584 | ND | ND | ND | ND | RI, MS |
| 2,4-Dimethylthiophene | 6.343 | 1205 | 1214 | Sulfurous | ND | ND | 17,991 | ND | ND | ND | AM, RI, MS |
| 3,4-Dimethylthiophene | 6.925 | 1243 | 1252 | Sulfurous, Meaty | ND | ND | 67,343 | ND | ND | ND | AM, RI, MS |
| 3-mercapto-2-butanone | 7.649 | 1310 | 1314 | ND | 7,265 | 70,660 | 1,561 | ND | 1,535 | 3,019 | RI, MS |
| 2-Methyl-3-furanthiol | 8.046 | 1338 | 1318 | Sulfurous, Meaty | ND | 2,148,091 | ND | ND | 91,730 | 12,709 | RI, MS |
| Dimethyl trisulfide | 8.944 | 1394 | 1396 | Sulfurous | 594,299 | 238,871 | 110,258 | 4,335 | 603,082 | 119,892 | RI, MS |
| Furfuryl mercaptan | 9.661 | 1441 | 1442 | Sulfurous | 2,598 | 128,129 | 5,230 | ND | 63,190 | 3,883 | RI, MS |
| 2-(1-Methylvinyl)thiophene | 9.767 | 1447 | 1433 | Sulfurous | ND | 1,916 | ND | ND | 57,298 | 9,001 | RI, MS |
| Methional | 10.025 | 1464 | 1464 | Sulfurous | 3,834 | 4,556 | 11,978 | 125,480 | 32,679 | 39,869 | RI, MS |
| Allyl (Z)-1-Propenyl disulfide | 10.398 | 1471 | 1464 | Sulfurous | ND | 2,788 | 9,878 | ND | ND | ND | RI, MS, AM |
| Unk C7H13NS | 10.474 | 1476 | NA | Sulfurous | 19,339 | 92,123 | ND | ND | ND | ND | AM |
| Diallyl disulfide | 10.483 | 1476 | 1475 | Sulfurous, Meaty, Legume | ND | ND | 21,120 | ND | ND | ND | RI, MS, AM |
| Thiazole, 2,4-dimethyl-5-propyl- | 10.988 | 1526 | 1515 | Sulfurous | ND | 3,919 | ND | 11,382 | ND | ND | RI, MS, AM |
| 2-Acetylthiazole | 13.061 | 1666 | 1656 | Sulfurous | 28,565 | 134,776 | 3,734 | 4,304 | 111,026 | 7,112 | RI, MS |
| 5,6-Dihydro-2,4,6-trimethyl-4H-1,3,5-dithiazine | 14.490 | 1766 | 1745 | Sulfurous | 439,174 | 767,273 | ND | ND | ND | ND | RI, MS |
| Unk C9H19NS2 | 16.782 | 1936 | NA | Meaty | 124,980 | 476,882 | ND | ND | 71,178 | 14,376 | AM |
| Unk C5H10S4 | 17.028 | 1955 | NA | ND | 57,586 | 58,178 | ND | ND | 211,296 | 2,601 | AM |
| Unk C9H19NS2 | 17.876 | 2025 | NA | Legume | 57,586 | 58,178 | ND | ND | 211,296 | 2,601 | AM |
| Unk C4H8S4 | 18.412 | 2074 | NA | Sulfurous | 88,392 | 442,465 | ND | ND | 7,342 | 3,609 | AM |
| Unk C10H21NS2 | 18.960 | 2131 | NA | Sulfurous | 93,582 | 145,442 | 5,113 | ND | 506,859 | 14,093 | AM |
| Furan, 3,3'-dithiobis[2-methyl- | 19.252 | 2165 | 2167 | Meaty | ND | 804,666 | ND | ND | ND | ND | RI, MS, AM |
| Unk C9H19NS | 19.443 | 2187 | NA | Meaty | 1,067,850 | 1,033,559 | ND | ND | 940,738 | ND | AM |
| Unk C2H4S4 | 21.964 | 2549 | NA | Sulfurous | 239,679 | 681,182 | ND | ND | 3,008 | 5,504 | AM |
| Total |  |  |  |  | 2,824,729 | 7,361,238 | 254,206 | 145,501 | 2,912,257 | 238,269 |  |

Two rows (Unk C5H10S4, RI 1955; Unk C9H19NS2, RI 2025) carry identical peak areas in all six products, so
one of them is probably a copy error in the SI.

Pyrazines block (in full):

| Name | RT / min | Expt RI | Lit RI | Odour | AM | SM-1 | GM | MM | PM | SM-2 | Identification |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pyrazine | 6.380 | 1224 | 1223 | ND | 2,845 | 20,924 | 894 | 2,276 | 4,327 | 6,594 | RI, MS |
| 2-Methylpyrazine | 7.149 | 1277 | 1277 | ND | 32,356 | 362,701 | 27,482 | 20,821 | 182,376 | 74,904 | RI, MS |
| 2,5-Dimethylpyrazine | 7.979 | 1332 | 1332 | Nutty | 24,751 | 253,443 | 11,385 | 12,571 | 201,214 | 85,149 | RI, MS |
| 2-Ethylpyrazine | 8.152 | 1343 | 1343 | ND | 4,217 | 52,630 | 6,013 | 5,060 | 42,701 | 26,386 | RI, MS |
| 2,3,5-Trimethylpyrazine | 9.068 | 1402 | 1412 | ND | 17,463 | 206,606 | 288,456 | 38,074 | 756,219 | 55,433 | RI, MS |
| 2-Ethyl-3,5-dimethylpyrazine | 9.824 | 1451 | 1455 | Nutty | 150,488 | 738,429 | 92,081 | 110,572 | 1,541,168 | 1,101,104 | RI, MS |
| 2,3-Dimethyl-5-ethylpyrazine | 10.079 | 1468 | 1460 | Nutty | 20,774 | 283,984 | 21,060 | 18,864 | 502,884 | 125,434 | RI, MS |
| 2,3-Diethyl-5-methylpyrazine | 10.569 | 1499 | 1497 | Nutty | 20,138 | 264,961 | 40,647 | 20,228 | 394,788 | 364,701 | RI, MS |
| Isoamyl methylpyrazine isomer 1 | 12.497 | 1627 | NA | ND | 35,329 | 196,849 | 82,733 | 23,665 | 223,939 | 472,937 | MS, AM |
| 2,5-Dimethyl-3-isoamylpyrazine | 12.999 | 1661 | 1666 | Nutty | 177,335 | 558,559 | ND | 58,951 | 480,573 | 1,142,558 | RI, MS |
| Isoamyl methylpyrazine isomer 2 | 13.448 | 1692 | NA | ND | ND | 62,507 | ND | 1,954 | 275,761 | 639,321 | MS, AM |
| 5,6,7,8-Tetrahydroquinoxaline | 14.197 | 1745 | 1752 | Nutty | 4,228 | 134,289 | ND | 2,759 | 26,968 | 67,451 | RI, MS |
| Total |  |  |  |  | 489,924 | 3,135,882 | 570,751 | 315,795 | 4,632,918 | 4,161,972 |  |

Carbonyls block (in full):

| Name | RT / min | Expt RI | Lit RI | Odour | AM | SM-1 | GM | MM | PM | SM-2 | Identification |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Hexanal | 4.703 | 1095 | 1096 | ND | 50,191 | 57,153 | 746,449 | 98,181 | 372,095 | 118,356 | RI, MS |
| 2-Heptanone | 5.913 | 1191 | 1182 | ND | 23,341 | 82,767 | 36,088 | 26,963 | 262,557 | 204,289 | RI, MS |
| Heptanal | 5.950 | 1193 | 1194 | ND | 15,762 | 44,893 | 19,506 | 20,243 | 83,071 | 84,979 | RI, MS |
| 3-Octanone | 6.915 | 1261 | 1262 | ND | 7,374 | 3,502 | ND | 31,997 | 26,017 | 13,117 | RI, MS |
| 2-Octanone | 7.378 | 1293 | 1287 | ND | 2,972 | 3,374 | 5,570 | 1,644 | 27,376 | 14,907 | RI, MS |
| Octanal | 7.440 | 1297 | 1297 | ND | 95,709 | 13,059 | 579,664 | 44,365 | 651,617 | 78,868 | RI, MS |
| Nonanal | 9.007 | 1398 | 1400 | ND | 618,435 | 97,031 | 2,355,486 | 155,903 | 1,927,209 | 283,168 | RI, MS |
| 3-Octen-2-one | 9.263 | 1415 | 1411 | ND | ND | ND | 1,256,971 | 54,906 | ND | 281,328 | RI, MS |
| 2-Octenal | 9.626 | 1438 | 1429 | ND | 30,516 | 8,771 | 475,255 | 14,092 | 78,929 | 39,479 | RI, MS |
| Furfural | 10.136 | 1471 | 1462 | ND | 46,820 | 107,479 | 775,819 | 1,839,564 | 723,670 | 317,046 | RI, MS |
| Benzaldehyde | 11.116 | 1535 | 1536 | ND | 944,100 | 611,468 | 742,246 | 960,758 | 773,874 | 571,900 | RI, MS |
| 2-Nonenal, (E)- | 11.202 | 1541 | 1544 | Fatty, Legume | 155,583 | 4,428 | 1,824,493 | 76,594 | 144,382 | 106,958 | RI, MS |
| 3,5-Octadien-2-one, (E,E)- | 11.767 | 1578 | 1570 | Fatty, Legume | ND | ND | 500,696 | 133,144 | ND | 515,811 | RI, MS |
| 2-Undecanone | 12.125 | 1602 | 1598 | Fatty | 24,952 | 40,388 | 89,589 | 396,228 | 318,057 | 60,214 | RI, MS |
| o-Tolualdehyde | 12.629 | 1636 | 1632 | ND | 3,700 | 1,947 | 57,696 | 40,630 | 29,403 | 13,593 | RI, MS |
| 2-Decenal, (E)- | 12.830 | 1650 | 1651 | Fatty | 303,048 | ND | 2,747,239 | 129,953 | 1,244,352 | 2,061,622 | RI, MS |
| Acetophenone | 13.069 | 1666 | 1662 | ND | 40,939 | 64,060 | 142,207 | 169,631 | 82,506 | 153,374 | RI, MS |
| 2-Undecenal | 14.384 | 1758 | 1751 | ND | 13,330 | 1,094 | 393,465 | 3,795 | 16,308 | 25,413 | RI, MS |
| 2,4-Decadienal, (E,E)- | 15.194 | 1816 | 1819 | Fatty | 111,328 | 41,247 | 2,459,319 | 568,793 | 857,815 | 1,485,642 | RI, MS |
| Benzeneacetaldehyde, a-ethylidene- | 16.845 | 1941 | 1929 | ND | 118,547 | 174,051 | 289,189 | 14,036 | 18,255 | 85,974 | RI, MS |
| 4-Hydroxy-3-methylacetophenone | 19.578 | 2203 | 2210 | ND | 3,429 | 181,097 | 1,362,350 | 4,583 | 1,330,153 | 458,888 | RI, MS |
| Total |  |  |  |  | 2,610,076 | 1,537,809 | 16,859,297 | 4,786,003 | 8,967,646 | 6,974,926 |  |

Odour-active rows of the "Alcohols" and "Others" blocks, plus indole (the remaining "Others" rows, the other
seven alcohols and the acids without odour are in the SI):

| Name | RT / min | Expt RI | Lit RI | Odour | AM | SM-1 | GM | MM | PM | SM-2 | Identification |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2-Undecanol (Alcohols block) | 13.811 | 1717 | 1717 | Fatty | ND | ND | 7,557 | 11,717 | ND | ND | RI, MS |
| 2-Pentylfuran | 6.550 | 1236 | 1237 | Meaty | 6,156 | 4,918 | 113,761 | 41,430 | 77,276 | 78,373 | RI, MS |
| Unk C7H11NO | 7.075 | 1271 | NA | Legume, Meaty | 7,340 | 20,891 | ND | ND | ND | 8,935 | AM |
| Unk C4H6NO | 8.472 | 1364 | NA | Meaty | ND | 71,685 | ND | ND | ND | 154,904 | AM |
| 2,5-Dimethylpyridine | 8.578 | 1371 | 1370 | Nutty | 3,474 | 11,969 | 38,756 | ND | 3,867 | 42,675 | RI, MS |
| Pyrrole | 10.807 | 1515 | 1520 | Nutty | 7,148 | 4,367 | ND | 3,723 | 8,872 | 1,859 | RI, MS |
| Unk C6H12NO | 14.689 | 1780 | NA | Legume | 2,117,068 | 2,041,100 | ND | ND | 1,794,662 | 1,192,313 | AM |
| Unk C10H17O3 | 20.175 | 2285 | NA | Meaty | 255,225 | 140,443 | 224,441 | 398,021 | 155,572 | 235,069 | AM |
| Unk C18H36O2 | 20.610 | 2349 | NA | Fatty, Meaty | 294,130 | 129,981 | 279,037 | 292,865 | 261,975 | 316,630 | AM |
| Indole | 21.383 | 2466 | 2466 | ND | 25,297 | 631,592 | 663,490 | 8,632 | 288,376 | 193,307 | RI, MS |
| Dodecanoic acid | 21.455 | 2477 | 2497 | Fatty | 101,996 | 287,333 | 901,803 | 143,008 | 105,528 | 85,435 | RI, MS |
| Total ("Others" block) |  |  |  |  | 4,120,126 | 5,501,526 | 9,439,738 | 1,488,889 | 5,531,204 | 7,387,802 |  |

What Table S4 says about the meaty sulfur set (within-compound comparisons across products only; TIC areas
of different compounds are not comparable, and none of this is a concentration):

- **2-methyl-3-furanthiol is not detected in beef (AM: ND).** It is the largest named sulfur peak in SM-1
  (2,148,091), small in PM (91,730) and SM-2 (12,709), absent in GM and MM. SM-1/PM = 23 (derived here,
  2,148,091 / 91,730). Its RI fit is the loosest of the named sulfur compounds: Expt 1338 against Lit 1318
  (20 units; derived here), and 1338 is 6 units from 2,5-dimethylpyrazine (1332).
- **Bis(2-methyl-3-furyl) disulfide** ("Furan, 3,3'-dithiobis[2-methyl-", RI, MS, AM; odour "Meaty") is
  found **only in SM-1** (804,666).
- **Furfuryl mercaptan**: SM-1 128,129, PM 63,190, beef 2,598 (SM-1/AM = 49, derived here).
- **DMTS**: beef 594,299 and PM 603,082 are similar; SM-1 238,871.
- **The dithiazine** occurs only in beef (439,174) and SM-1 (767,273).
- **2-acetylthiazole**: SM-1 134,776, PM 111,026, beef 28,565.
- **Methional** is largest in MM (125,480; beef 3,834). Furfural is also largest in MM (1,839,564; beef
  46,820).
- **2-(methoxymethyl)furan** (RI 1243) has odour **"ND"** in Table S4. The only odour-active peak at RI 1243
  is **3,4-dimethylthiophene** ("Sulfurous, Meaty", GM only). The main text's link of the sulfurous note at
  RI ~1243-1249 to 2-(methoxymethyl)furan is therefore not supported by the SI's own odour column.
- **2-pentylfuran** is described as "Meaty" by GC-O (not fatty), with the largest area in GM (113,761).
- Lipid-derived carbonyls: (E,E)-2,4-decadienal is lowest in SM-1 (41,247) and beef (111,328) and highest
  in GM (2,459,319), confirming the main text's p. 10 statement. (E,E)-3,5-octadien-2-one is absent in beef,
  SM-1 and PM and present in GM, MM and SM-2.

### 2x.4 Sensory means (Supplementary Table 3)

Caption: "Mean (±SEM) sensory ratings of 11 burger patty samples on a 0 – 100 Visual Analogue Scale (VAS)";
superscript 1 = p-value of the sample main effect; letters differ at p ≤ 0.05 (Bonferroni), "'a' always
represents the lowest value". Transcribed in full (row order as printed):

| Samples | Meaty Odour Intensity | Legume Odour Intensity | Off-odour Intensity | Meaty Flavour Intensity | Legume Flavour Intensity | Salty Taste Intensity | Savoury Taste Intensity | Off-flavour Intensity | Juiciness | Chewiness | Oily Mouthfeel | Flavour Aftertaste |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AM | 72.48 (4.79) e,f | 7.44 (5.07) a,b | 30.60 (6.04) a, b, c | 77.24 (5.03) e, f | 3.24 (5.31) b | 14.97 (4.57) b | 25.41 (4.74) b | 22.35 (5.73) a, b | 26.49 (3.90) c | 73.31 (4.55) e | 23.84 (4.76) a, b, c | 54.63 (5.41) a, b |
| AP | 89.77 (4.79) f | 4.25 (5.07) a | 13.38 (6.04) a | 93.23 (5.03) f | 2.32 (5.31) a, b | 52.83 (4.57) d | 67.98 (4.74) d | 5.20 (5.73) a | 67.32 (3.90) f | 72.03 (4.55) e | 52.71 (4.76) e | 68.78 (5.41) b, c |
| PM | 35.57 (4.79) c,d | 45.57 (5.07) c, d, e | 35.10 (6.04) a, b, c, d | 45.63 (5.03) b, c, d | 43.08 (5.31) c | 41.25 (4.57) c, d | 52.80 (4.74) c | 21.68 (5.73) a, b | 47.63 (3.90) d, e | 63.07 (4.55) d, e | 46.12 (4.76) d, e | 55.72 (5.41) a, b |
| PP | 40.41 (4.79) c, d | 46.21 (5.07) d, e | 28.32 (6.04) a, b, c | 54.96 (5.03) c, d | 39.49 (5.31) c | 48.82 (4.57) c, d | 65.18 (4.74) d | 17.46 (5.73) a, b | 81.71 (3.90) g | 52.85 (4.55) b, c, d | 59.00 (4.76) e | 60.90 (5.41) a, b, c |
| GM | 4.47 (4.79) a,b | 62.41 (5.07) e, f | 58.30 (6.04)  e | 4.67 (5.03) a | 68.78 (5.31) d | 51.59 (4.57) c, d | 52.21 (4.74) c, d | 54.03 (5.73) c, d | 11.27 (3.90) b | 44.30 (4.55) b, c | 23.91 (4.76) a, b, c | 74.25 (5.41) c |
| GP | 5.65 (4.79) b | 67.13 (5.07) f | 42.59 (6.04) b, c, d, e | 5.50 (5.03) a | 68.59 (5.31) d | 49.74 (4.57) c, d | 52.04 (4.74) c, d | 35.28 (5.73) b, c | 13.47 (3.90) b, c | 57.89 (4.55) c, d, e | 20.84 (4.76) a, b | 60.00 (5.41) a, b, c |
| MM | 9.66 (4.79)b | 54.57 (5.07) e, f | 54.74 (6.04) d, e | 6.13 (5.03) a | 48.18 (5.31) c, d | 9.07 (4.57) a, b | 4.13 (4.74) a | 58.98 (5.73) d | 4.15 (3.90) a, b | 20.23 (4.55) a | 6.28 (4.76) a | 47.51 (5.41) a |
| MP | 54.48 (4.79) d, e | 30.33 (5.07) c, d | 24.97 (6.04) a, b, c | 45.36 (5.03) b, c, d | 35.88 (5.31) c | 45.64 (4.57) c, d | 44.88 (4.74) c | 33.58 (5.73) b, c | 48.79 (3.90) d, e | 38.37 (4.55) b | 28.95 (4.76) b, c, d | 51.84 (5.41) a, b |
| SM-2 | 31.83 (4.79) c | 43.22 (5.07) c, d, e | 33.33 (6.04) a, b, c, d | 35.16 (5.03) b | 49.90 (5.31) c, d | 37.56 (4.57) c | 42.32 (4.74) b, c | 32.10 (5.73) b, c | 40.76 (3.90) d | 45.74 (4.55) b, c | 41.07 (4.76) c, d, e | 51.59 (5.41) a, b |
| SM-1 | 66.49 (4.79) e | 25.11 (5.07) b, c | 20.46 (6.04) a, b | 62.57 (5.03) d, e, f | 37.77 (5.31) c | 38.39 (4.57) c, d | 44.80 (4.74) c | 20.97 (5.73) a, b | 43.17 (3.90) d, e | 59.28 (4.55) c, d, e | 41.47 (4.76) c, d, e | 57.40 (5.41) a, b, c |
| SP | 32.31 (4.79) c | 47.07 (5.07) d, e, f | 44.36 (6.04) c, d, e | 41.28 (5.03) b, c | 47.64 (5.31) c | 45.03 (4.57) c, d | 55.66 (4.74) c, d | 33.56 (5.73) b, c | 55.50 (3.90) e, f | 50.45 (4.55) b, c, d | 53.33 (4.76) e | 63.61 (5.41) a, b, c |
| P-value1 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| DF | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 | 10,210 |
| F Ratio | 46.00 | 22.69 | 8.137 | 52.51 | 24.76 | 20.85 | 21.04 | 10.81 | 70.06 | 20.11 | 18.71 | 4.62 |

- Meaty odour: SM-1 66.49 (e) is not separable from beef mince AM 72.48 (e,f); the only other plant product
  that reaches "d,e" is MP 54.48 (Quorn patty, with smoke flavouring and roasted barley malt extract), which
  was not in the GC-MS/O set. PM 35.57 and SM-2 31.83 are in the "c" group.
- Meaty flavour: SM-1 62.57 (d,e,f) overlaps AM 77.24 (e,f); PM 45.63 and SM-2 35.16 are lower.
- Inside the GC-MS/O set, the meaty-odour order SM-1 > PM ≈ SM-2 > MM ≈ GM matches the order of the MFT
  areas (SM-1 > PM > SM-2; MM and GM not detected), but beef (the meatiest) has no detected MFT
  and little FFT; its sulfur areas are DMTS, the dithiazine and unknown sulfur compounds (Unk C9H19NS
  1,067,850). This is an observation from the two tables, not an analysis the paper reports.

### 2x.5 Supplementary Figure 1 (training aromagrams)

Three stacked traces against retention time (0-25+ min): (A) untrained and (B) trained panel aromagrams
(citation proportion 0-1, coloured by descriptor: Meaty, Legume, Fatty, Sulfurous, Nutty, Others) and (C)
the chromatogram of the 9-standard sniffing mixture, peaks numbered as in main-text Table 4. Peak positions,
read from graph, approx.: 1 at 7.4, 2 at 8.3, 3 at 9.3, 4 at 10.0, 5 at 10.4, 6 at 15.6, 7 at 18.4, 8 at
19.5, 9 at 22.1 min. Training visibly removes many scattered noise citations (B is sparser than A), consistent
with the main text's noise counts (2910 to 1801). These standard RTs sit about 0.3 min later than the Table
S4 RTs of the same compounds (e.g. DMTS 8.944, furfuryl mercaptan 9.661, decadienal 15.194; derived here by
comparison), so the training mixture ran under slightly different conditions. No intensities or amounts.

### 2x.6 What the SI changes in this dossier

1. **Identification level.** The main-text uncertainty is resolved in the weak direction: by Table S4's
   own codes no odorant was confirmed with an authentic standard ("AS" is never used). DMTS, furfuryl
   mercaptan and methional had standards available (Table S2) but are coded "RI, MS". 2-methyl-3-furanthiol,
   2-(1-methylvinyl)thiophene and 2-(methoxymethyl)furan had no standard at all.
2. **MFT is a plant-product marker here, not a beef one.** It was not detected in the beef mince.
3. **The "flavouring, not protein" reading** is consistent with the ingredient lists: the product richest in
   MFT/FFT/the MFT disulfide (SM-1) is the one with thiamine HCl + yeast extract + leghemoglobin. No product
   lists cysteine.
4. **Sensory means** are now available (Table S3).
5. Still no concentrations: Table S4 is TIC peak areas without internal standard or response factors.

## 3. What it means for the model and for formulation

**Target profile.** In this study, beef is a sparse profile (21 cited odorants) carried by sulfurous and meaty
notes. The PBMAs are denser (more than 30 peaks each) and over-supplied with fatty and legume notes. The
implied plant-based target has two parts:
- **Raise the sulfurous/meaty thiol set**: 2-methyl-3-furanthiol, furfuryl mercaptan (2-furfurylthiol),
  dimethyl trisulfide, methional and possibly 2-(1-methylvinyl)thiophene. Their citation proportions track
  rated meatiness. Caveat from the SI (section 2x.3): all are putative (RI + MS), and MFT was not detected
  in the beef at all, so for MFT the "target" is the best plant product (SM-1), not beef. Beef's sulfur
  signal is DMTS, the dithiazine and unidentified sulfur compounds.
- **Lower the lipid-derived carbonyls**: (E,E)-2,4-decadienal, (E,E)-3,5-octadien-2-one, 2-undecanol,
  (E)-2-nonenal, (E)-2-decenal and 2-undecanone. These track off and legume ratings, and beef had the least
  decadienal.

That the best-scoring plant products (SM-1 soy, PM pea) were richest in the thiols and pyrazines, and that
flavouring rather than protein separated them, points to precursor dosing as the lever, not protein choice.
The ingredient lists (Supplementary Table 1, section 2x.1) narrow this: SM-1, by far the richest in MFT, FFT
and the MFT disulfide, is the only analysed product listing thiamine HCl, yeast extract and leghemoglobin
together; PM lists dried yeast but no thiamine; **no product lists cysteine**, so any cysteine would have to
sit inside the undisclosed "Natural Flavors". All of this is correlational (MFA, ingredient presence),
unquantified, and on putative identities only.

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
against: citation proportions are not amounts, and the SI's Table S4 gives TIC peak areas without internal
standard or response factors. At most these support within-compound, between-product orderings (e.g. MFT
SM-1 >> PM > SM-2, not detected in beef), and the product precursor contents are not known.

## What it does not give

- Concentrations or OAVs of any compound. GC-O citation proportions are detection-frequency data at an
  undiluted extract.
- Confirmation of any odorant with an authentic standard: Supplementary Table 4 codes every named
  compound "RI, MS" (or accurate mass only); the "AS" code is defined but never used.
- Correlation coefficients or significance for the sensory-odorant links.
- Precursor amounts. Supplementary Table 1 gives the full ingredient lists (section 2x.1), but no
  quantities beyond a few printed percentages, and the composition of "Natural Flavors" / "Flavouring" is
  not disclosed. No product lists cysteine.
- Cooking surface. The hydro-distillation (1 h simmer in water) can itself generate or lose volatiles; the
  paper does not test for artefacts.
