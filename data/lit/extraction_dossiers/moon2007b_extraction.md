# Moon & Li-Chan 2007 (Food Res Int 40:1239) — EXTRACTION (file stem `moon2007b`)
DOI: **10.1016/j.foodres.2007.08.003** (printed on p. 1239 and in the PDF metadata title). Food Research International
**40**(10):**1239-1248**, December 2007; received 23 May 2007, accepted 14 Aug 2007. Publisher page
https://www.sciencedirect.com/science/article/abs/pii/S0963996907001378.
**Crossref check (2026-11-02, `api.crossref.org/works/10.1016/j.foodres.2007.08.003`)**: title "Changes in aroma
characteristics of simulated beef flavour by soy protein isolate assessed by descriptive sensory analysis and gas
chromatography", authors Soo-Yeun Moon, Eunice C.Y. Li-Chan, Food Res Int 40(10):1239-1248, issued 2007-12. **Matches
the PDF in every field.** (Crossref JSON parsed with a one-line `python3 -I`; no repo code run.)
**Stem check: `moon2007b` = .003 = pp. 1239-1248; `moon2007a` = .002 = pp. 1227-1238. Not swapped.** [!] Logically this
is the *first* paper of the pair: it defines the five sensory attributes, the panel and the 15 GC "indicator peaks" that
`moon2007a_extraction.md` (the follow-up, which cites this one as "accepted for publication") then uses. **The
headspace-suppression numbers the repo's literature search attributes to this DOI are in a, not here** (sec. 6).

### Commercial vegetarian simulated beef flavour (SBF, Mastertaste F96·49) at 150 or 500 mg × commercial soy protein isolate (Solae) at 0 / 250 / 500 mg in 5 g distilled water × unheated or heated 98 °C / 30 min (a 2 × 3 × 2 full factorial, 12 samples); descriptive sensory analysis (8 trained panellists, 6 retained, 2 replications, 5 attributes) and HS-SPME GC-FID (DVB/CAR/PDMS, 60 min at 60 °C). Output: **sensory means per sample (Table 7)**, correlations of sensory scores with SBF and SPI dose (Table 6), and 15 GC indicator peaks chosen by correlation with "beefy" and GC-O detection frequency (Table 8). **No GC peak area, no % headspace loss, no binding constant is printed**; the per-sample GC data exist only behind a correlation matrix in the online supplement (not on disk). No sulfur compound named.

**Read 2026-11-02 from `data/articles/moon2007b.pdf`** (10 pp., Elsevier, born-digital, text layer good). Prose from
`pdftotext -layout`. **Tables 2-8 checked on 150 dpi page renders by eye** (scratchpad `ingest/moon2007b_work/p-06.png`
Tables 4-6 + Fig. 1, `p-07.png` Table 7, `p-08.png` Table 8 + Fig. 2); minus signs in Tables 4 and 6 are visible only on
the render (the text layer drops them). Table 1 (design) and Tables 2-3 from the text layer, which is clean for them.
Fig. 1 (cobweb) plots the Table 7 means and was not digitised; Fig. 2 (PCA biplots) has no numeric axis and was not
digitised (PC variance shares are printed in the text).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Food, Nutrition and Health Program, Faculty of Land and Food Systems, **The University of British
Columbia**, Vancouver; corresponding (and senior) author **Eunice C.Y. Li-Chan**; first author Soo-Yeun Moon. GC run on
an instrument at the British Columbia Institute of Technology (acknowledgements); NSERC-funded. **Not in
`data/keys/laboratories.yml`: NEW**, proposed id **`ubc_li_chan`**, lab string **"UBC Vancouver (Li-Chan)"** — the same
proposal as `moon2007a_extraction.md`. No other dossier on disk is from this group. Same-lab companion:
**`moon2007a_extraction.md`** (same SBF, same SPI lot 02060631-532, same fibre / vial / water bath, five of its seven
panellists from this panel). Moon & Li-Chan 2004 (Food Chem 88:141) and Moon, Cliff & Li-Chan 2006 (Food Res Int 39:294)
are cited for the SPME method and the GC-MS / GC-O identities and are **not on disk**.

## 1. What this paper answers for the model

| model need | what Moon 2007b gives | tag |
|---|---|---|
| headspace loss of meat-flavour volatiles on adding plant protein | **not printed**: GC peak areas enter only correlations; per-sample areas sit (if anywhere) in the SI correlation matrix, not on disk | [!] |
| binding constants / % bound | **none** | [!] |
| protein-concentration dependence | **sensory only**: 0 / 250 / 500 mg SPI (≈ 0 / 4.2 / 8.1 % w/w protein in water + SPI [Z]) at two SBF doses, Table 7 | [M] |
| effect of heating protein with the flavour | 98 °C / 30 min of the SBF-SPI mixture, sensory only (Table 7); GC heated arm only in correlations | [M] |
| identities of beef-correlated volatiles in a commercial meat flavour | Table 8: 15 indicator peaks, 10 tentatively named (furans, terpenes, 3 alkylpyrazines, decanal), 5 unknown | [M] / [C] |
| sulfur odorants (MFT, FFT, methional, sulfides) | **none named** in the paper | [!] |
| perceptual suppression of beefy notes by soy protein | Table 6 r(beefy, SPI) = −0.696; Table 7 dose series | [M] |
| temperature / pH | single SPME T (60 °C); **pH never stated or measured** (distilled water) | [!] |

**Answers.**

1. **Compounds, protein, concentration, heating, pH, temperature, time, design** (confidence ~95 %). *Flavour*: one
   commercial **vegetarian beef-type reaction flavour** (Mastertaste F96·49: maltodextrin, autolysed yeast extract,
   natural flavours, onion powder, SiO2), ≥70 volatiles by GC-MS and 49 odour-active by GC-O in the authors' earlier work
   [C]; **no individual odorant is dosed or quantified**. *Protein*: **one commercial SPI** (Solae, 89.3 ± 0.1 % protein,
   N × 6.25); **no glycinin / β-conglycinin, no other plant protein**. *Doses*: SBF 150 or 500 mg ("low" / "high") and
   SPI 0 / 250 / 500 mg ("no" / "medium" / "high"), each in **5 g distilled deionised water** → SPI protein ≈ 0 / 4.2 /
   8.1 % w/w of water + SPI, ≈ 0 / 45 / 89 g/L of water [Z]. *Heating*: half the samples heated **98 °C / 30 min in a water bath, SBF and
   SPI together** (Table 1 footnote c). *pH*: **not stated, not measured, no buffer.** *Temperature / time*: sensory
   vials held **60 °C / 20 min** before serving; GC by HS-SPME **60 min at 60 °C**, stirred, no separate equilibration.
   *Design*: **a 2 (SBF) × 3 (SPI) × 2 (heat) full factorial, 12 samples** (Table 1), each in triplicate for GC and in
   two replication blocks for sensory. **Not a 3² factorial** (~95 %); the "3 2 full factorial" quote in the repo's
   literature-search table is not in this paper (sec. 6).
2. **Quantities** (confidence ~95 %). **Sensory**: mean intensity in cm on a 15 cm unstructured line (Table 7, n = 12:
   6 panellists × 2 replications), with Duncan letters; relative to a reference (150 mg SBF + 250 mg SPI) anchored
   "moderate". **GC**: FID peak areas from HS-SPME, used only as Pearson r with sensory scores (Table 8: r with "beefy"
   0.590-0.931 for the 15 indicator peaks; text: 75 of 112 peaks significantly correlated with some attribute, 64
   positively with roasted / beefy / yeasty and negatively with soymilk / cereal, 11 the reverse). **No absolute or
   relative headspace values, no % reduction, no binding constant, static or dynamic.** The headspace-suppression
   numbers of this system are in `moon2007a` Tables 1-2 (SPIF column; single SPI level, ~43 g/L).
3. **Sensory panel** (confidence ~90 %). Roasted, beefy and yeasty rise with SBF dose (r = 0.662 / 0.626 / 0.668) and
   fall with SPI dose (r = −0.692 / −0.696 / −0.659); soymilk-like and cereal do the reverse (r with SPI 0.714 / 0.711).
   **Beefy, unheated**: low SBF 8.47 → 6.71 → 2.54 cm at 0 / 250 / 500 mg SPI (−21 %, −70 % [Z]); high SBF 10.88 → 9.93 →
   7.72 (−9 %, −29 % [Z]). Roasted follows (low SBF 9.63 → 5.69 → 2.66; high 10.82 → 10.23 → 7.39). Soymilk-like rises
   3.52 → 6.73 → 9.98 (low SBF) and 2.53 → 3.25 → 6.25 (high SBF). **Tripling SBF (S1 → S4) is cancelled by 250 mg SPI**
   (beefy 8.47 vs 9.93, same Duncan group). **Heat (98 °C / 30 min) changes little** except at low SBF / medium SPI, where
   beefy falls 6.71 → 3.73 and roasted 5.69 → 3.46 (different Duncan letters); the authors call the heated trends "very
   similar" [Q]. PC1 carries 97.8 % of variance and splits SBF notes from SPI notes. [!] Because soy notes rise as beef
   notes fall, the beef suppression is **not separable into binding and perceptual masking** with these data
   (confidence ~70 % that masking contributes a substantial part; the paper does not test it).
4. **Usable data sets; T- and protein-dependence** (confidence ~80 %). One quantitative set the matrix ledger can hold:
   **Table 7 sensory means** (with Table 6 correlations and Table 2 ANOVA as its statistics) — the only
   protein-dose series (0 / 4.2 / 8.1 % protein) of a meat flavour on a plant protein in the corpus, **but sensory, not
   headspace**. The heating arm (98 °C) is a protein-plus-flavour cook, sensory only. **No temperature dependence of
   binding** (one GC temperature, numbers not printed). The **SI correlation matrix** is listed as `supplement_missing`
   (P(it prints per-sample peak areas, which would give a headspace dose series) ~25 %). Table 8 (indicator-peak
   identities and r) is not a matrix quantity and gets no entry.

**What the matrix code assumes for meat odorants on soy, and what this paper says** (registry as in
`moon2007a_extraction.md` sec. 1):

| odorant / class | code (soy) | Moon 2007b | verdict |
|---|---|---|---|
| MFT, FFT (thiol) | no constant; `SULFUR_BINDING_GAP`; thiol covalent class, ceiling zero | not named; GC-FID | **silent** (~95 %) |
| methional | not registered | not named | **silent** |
| alkylpyrazines | no binding class row (`no_binding_constant_for_class`) | three alkylpyrazines among the beef indicator peaks (r 0.590-0.601); no headspace numbers | **silent on magnitude**; a gives the numbers |
| n-aldehydes (decanal) | `kg_nonanal_soy` (dialysis) + 2.81×/CH2 | decanal an indicator peak (r 0.715), no numbers | **silent** |
| overall direction "soy protein lowers perceived meat flavour" | the layer ships no protein correction for sulfur or pyrazines, so a soy matrix predicts no loss for them | beefy −21 to −70 % (low SBF), −9 to −29 % (high SBF) with 4-8 % SPI | **contradicts a zero-correction default in direction** (~65 %, since part may be masking); gives no constant |
| heat × binding | "no scalar, three signs" | 98 °C / 30 min cook of SBF + SPI: sensory mostly unchanged, one cell down | **consistent** with "no scalar"; weak |

## 2. Conditions

- **SBF** [M]: "vegetarian beef type flavour F96·49", Mastertaste (Arlington Heights IL): maltodextrin, autolysed yeast
  extract, natural flavours, onion powder, silicon dioxide. Same lot as a. Doses 150 mg (L-SBF) or 500 mg (H-SBF) in 5 g.
- **SPI** [M]: Solae (St. Louis MO), lot 02060631-532, **89.3 ± 0.1 %** protein (LECO FP-428, N × 6.25). In `moon2007a`
  the same lot shows **no DSC transition 25-180 °C** (already denatured) and free SH 0.6 µmol/g protein. Doses 0 /
  250 / 500 mg (N-, M-, H-SPI).
- **Medium** [M]: 5 g distilled deionised water; **no buffer; pH not reported** [!]. Vessel: 15 mL Supelco vial,
  polypropylene hole cap, PTFE/silicone septum, magnetic stir bar (GC); foil-covered 15 mL vials (sensory).
- **Heat treatment** [M]: half of the 12 samples (code suffix H) "heated in water bath at 98 °C for 30 min" with SBF and
  SPI together (Table 1 footnote c). Cooling not described.
- **Sensory preparation** [M]: each of the 12 samples made at batch scale, then aliquoted into vials, held 60 °C / 20 min
  before serving [!] (a notes this aliquoting may lose volatiles and changes it to per-vial addition). Reference R:
  150 mg SBF + 250 mg SPI (the S3 recipe), anchored "moderate" (7.5 cm). Coffee beans offered between samples.
- **GC preparation** [M]: SBF + SPI mixed in 5 g water in the GC vial, HS-SPME 50/30 µm DVB/CAR/PDMS **60 min at 60 °C**
  (water bath ±2 °C), stirred; **all treatments in triplicate, mean of three** reported (only as correlations).
- **Panel** [M]: 10 trained (8 women, 2 men; UBC graduate students and staff of a soy-based meat-alternative company),
  8 h training in four 2 h sessions; 8 took part in the main test; **panellists 4 and 5 excluded** for low concordance /
  repeatability, leaving 6 (Table 7 n = 12). 6 sessions × 4 samples + reference; 3-digit codes, randomised within
  replication blocks.
- **Water activity / moisture**: dilute aqueous (~88-97 % water [Z]).

## 3. Analytical methods

- **Descriptive analysis** [M]: adapted from Zook & Pearce 1988; unstructured 15 cm line anchored 1.3 cm from each end
  ("Slight", "Moderate", "Intense"); score = distance from origin (cm). Attributes and definitions agreed in training:
  roasted (roast-beef crust top note), beefy (beef vs pork / chicken), yeasty (off-flavour, "a roasted smell gone bad"),
  soymilk-like (green, cooked soybean), cereal (cooked oatmeal / pasta). Training references: Knorr OXO beef bouillon,
  oven-roasted Canada AAA steak (crust vs interior), raw ground beef, Werther's caramel, wet cardboard, yeast extract,
  unsweetened soy beverage, cooked oatmeal, hay.
- **Statistics** [M]: ANOVA (sample, panellist, replication and interactions; Minitab 13.30), Duncan's test; adjusted F
  with MS(sample × replication) as the denominator where that interaction was significant (O'Mahony 1986); Pearson r for
  panel concordance and for sensory vs dose and vs GC peak area; PCA on the 12 × 5 mean matrix.
- **GC** [M]: HS-SPME / GC-FID per Moon & Li-Chan 2004 and Moon 2006 [C] (column and programme not restated); peak
  identities **tentative**, by retention time against Moon 2006 GC-MS; DF = GC-O detection frequency of 8 sniffers in
  Moon 2006 [C]. Indicator-peak rule: r with beefy significant (p ≤ 0.05) **and** DF > 50 % (≥ 4 of 8). **No
  quantification basis is reported because no quantity is reported.**

## 4. Data

### 4a. Table 1 (p. 1242): design [M]

| code | SBF (mg) | SPI (mg) | water (g) | heat (98 °C, 30 min) |
|---|---|---|---|---|
| S1 / S1H | 150 | 0 | 5 | no / yes |
| S2 / S2H | 500 | 0 | 5 | no / yes |
| S3 / S3H | 150 | 250 | 5 | no / yes |
| S4 / S4H | 500 | 250 | 5 | no / yes |
| S5 / S5H | 150 | 500 | 5 | no / yes |
| S6 / S6H | 500 | 500 | 5 | no / yes |

### 4b. Table 2 (p. 1243): ANOVA F-values, 8 panellists [M]

| attribute | sample S (df 11) | panellist P (7) | replication R (1) | S×P (77) | P×R (7) | S×R (11) | adjusted F, sample (11) |
|---|---|---|---|---|---|---|---|
| roasted | 21.06*** | 2.33* | 3.35 | 1.19 | 0.70 | 0.86 | 24.49*** |
| beefy | 15.74*** | 1.54 | 4.34* | 1.24 | 1.05 | 1.26 | 12.48*** |
| yeasty | 8.53*** | 11.66*** | 2.02 | 1.44 | 0.89 | 0.90 | 9.44*** |
| soymilk-like | 17.31*** | 6.94*** | 0.09 | 1.18 | 0.67 | 3.16** | 5.48** |
| cereal | 14.18*** | 7.42*** | 2.49 | 0.96 | 0.86 | 2.32* | 6.12** |

Tables 3-5 (panellist means with Duncan letters; panellist-vs-rest r; per-panellist replication F) are panel-performance
diagnostics; one line: panellists 4 and 5 have r −0.04 to 0.56 with the rest (others mostly > 0.6, yeasty weaker) and were
dropped. Not transcribed further (no model use).

### 4c. Table 6 (p. 1244): Pearson r of mean sensory score with dose (n = 12 sample means) [M, signs from render]

| attribute | r with SBF (P) | r with SPI (P) |
|---|---|---|
| roasted | 0.662* (0.019) | −0.692* (0.013) |
| beefy | 0.626* (0.030) | −0.696* (0.012) |
| yeasty | 0.668* (0.017) | −0.659* (0.020) |
| soymilk-like | −0.634* (0.027) | 0.714** (0.009) |
| cereal | −0.638* (0.026) | 0.711** (0.010) |

### 4d. Table 7 (p. 1245): mean intensity (cm, 15 cm line), 6 panellists × 2 reps, n = 12 [M]

Rearranged from the paper's per-attribute ascending lists into a sample × attribute grid; Duncan letters as printed
(different letters = p ≤ 0.05 within an attribute).

| sample (SBF / SPI / heat) | roasted | beefy | yeasty | soymilk-like | cereal |
|---|---|---|---|---|---|
| S1 (150 / 0 / –) | 9.63 ef | 8.47 bcd | 7.50 cd | 3.52 ab | 3.65 ab |
| S1H (150 / 0 / 98 °C) | 8.53 def | 10.13 cd | 8.43 cd | 4.02 ab | 3.73 ab |
| S2 (500 / 0 / –) | 10.82 fg | 10.88 d | 9.29 d | 2.53 ab | 2.32 ab |
| S2H (500 / 0 / 98 °C) | 12.14 g | 10.78 d | 9.41 d | 1.85 a | 1.73 a |
| S3 (150 / 250 / –) | 5.69 bc | 6.71 b | 4.17 ab | 6.73 cd | 6.84 cde |
| S3H (150 / 250 / 98 °C) | 3.46 ab | 3.73 a | 3.20 a | 7.87 de | 7.91 def |
| S4 (500 / 250 / –) | 10.23 fg | 9.93 cd | 7.97 cd | 3.25 ab | 3.57 ab |
| S4H (500 / 250 / 98 °C) | 10.63 fg | 9.92 cd | 8.43 cd | 2.88 ab | 2.33 ab |
| S5 (150 / 500 / –) | 2.66 a | 2.54 a | 3.11 a | 9.98 e | 9.99 f |
| S5H (150 / 500 / 98 °C) | 1.93 a | 2.69 a | 2.80 a | 8.98 e | 9.25 ef |
| S6 (500 / 500 / –) | 7.39 cde | 7.72 bc | 6.32 bc | 6.25 cd | 6.19 cd |
| S6H (500 / 500 / 98 °C) | 6.83 cd | 7.19 b | 6.97 cd | 4.73 bc | 4.69 bc |

Derived [Z]: beefy relative to the no-SPI sample at the same SBF dose and heat: unheated low SBF 1.00 / 0.79 / 0.30,
high SBF 1.00 / 0.91 / 0.71; heated low SBF 1.00 / 0.37 / 0.27, high SBF 1.00 / 0.92 / 0.67 (0 / 250 / 500 mg SPI).
The suppression is proportionally larger at the lower flavour dose, as a saturable binder or a masking background
would both produce.

### 4e. PCA (text, p. 1245) [M]

PC1 / PC2 / PC3 explain 97.8 / 1.0 / 0.7 % of variance (eigenvalues 4.8892 / 0.0478 / 0.0342). PC1: cereal and
soymilk-like positive, roasted / yeasty / beefy negative; sample order along PC1 (beef side first) S2, S4, S1, S6, S3, S5.
PC3 is mainly yeasty; heated samples except S2H sit higher on PC3 [Q].

### 4f. Table 8 (p. 1246): the 15 indicator peaks for "beefy" [M]

| IP | RT (min) | identification (Moon 2006 peak no.) | DF (%) | r with beefy |
|---|---|---|---|---|
| IP1 | 1.717 | 3-methylfuran (P1) | 50 | 0.710* |
| IP2 | 7.832 | 2-acetylfuran (P7) | 75 | 0.658* |
| IP3 | 12.046 | δ-3-carene (P17) | 88 | 0.872*** |
| IP4 | 15.772 | 2-ethyl-3,6-dimethylpyrazine (P24) | 88 | 0.601* |
| IP5 | 16.133 | unknown | 100 | 0.852*** |
| IP6 | 16.864 | unknown | 100 | 0.827** |
| IP7 | 18.199 | unknown | 75 | 0.605* |
| IP8 | 19.744 | 2,3-diethyl-5-methylpyrazine (P30) | 50 | 0.590* |
| IP9 | 22.062 | decanal (P35) | 75 | 0.715** |
| IP10 | 24.103 | 2-isoamyl-6-methylpyrazine (P38) | 50 | 0.599* |
| IP11 | 25.374 | unknown | 75 | 0.598* |
| IP12 | 25.703 | unknown | 75 | 0.794** |
| IP13 | 28.204 | δ-elemene (P44) | 50 | 0.931*** |
| IP14 | 30.291 | β-cubebene (P49) | 63 | 0.850*** |
| IP15 | 35.680 | calamenene (P71) | 75 | 0.883*** |

[!] The text (p. 1246) says DF "greater than 50 %" (then "4 or more ... out of 8"); five peaks sit at exactly 50 % (4 of 8), so "greater than" should read "at least". Retention
times are 0.02-0.13 min earlier than the same IPs in `moon2007a` (a separate run).

## 5. Flags and internal consistency

1. [!] **No headspace quantity is printed.** The abstract and conclusions speak of SPI-induced suppression of beef
   volatiles, but the GC results appear only as correlation coefficients; the numbers are in the follow-up (a) at one SPI
   level, or in this paper's SI (not on disk).
2. [!] **Binding and perceptual masking are confounded**: SPI adds soymilk / cereal odour while beef notes fall; the
   authors' "interaction between SBF and SPI considerably hindered the beefy notes from releasing" [Q] is a hypothesis.
3. [!] **pH not stated; unbuffered water**; the SBF (yeast extract) and SPI set it.
4. [!] **No sulfur compound named.** 5 of 15 indicator peaks are unknown and the detector is an FID; a sulfur odorant
   could hide among the unknowns, but nothing here shows one.
5. [!] Text p. 1245 says the three beef attributes decreased "significantly, in the order of S2, S4, and S6"; Table 7
   letters show S2 vs S4 **not** significantly different for roasted (fg / fg), beefy (d / cd) or yeasty (d / cd) — only
   S6 separates.
6. [!] Table 4 footnote: bold = r > 0.600 and p ≤ 0.05; panellist 1 beefy r = 0.597 (P 0.002) is not bold although
   significant. Table 5 footnote: bold = p ≤ 0.05, yet 0.098 and 0.095 are bold (the text uses p < 0.1). Table 5
   panellist 8 yeasty P prints "(0.2810" (typo).
7. [!] Sensory samples were batch-made and aliquoted (volatile loss possible, as a itself says); GC samples were made
   per vial. Sensory at 60 °C / 20 min headspace, GC at 60 °C / 60 min on a fibre: not the same headspace.
8. [!] Correlation-based peak selection on n = 12 sample means with 112 peaks tested and no multiplicity correction:
   some of the 75 "significant" peaks are expected false positives (~5-6 at α = 0.05 [Z]).
9. [!] The heat arm cooks flavour and protein together at 98 °C, so Maillard / degradation chemistry of the SBF itself
   (S1H vs S1: beefy 8.47 → 10.13) is confounded with any heat effect on binding.
10. [!] Cross-paper check (`moon2007a` Fig. 2 SPIF arm = this paper's S3 recipe): beefy 6.7 vs 6.71, soymilk 6.9 vs
    6.73, cereal 6.4 vs 6.84 agree; roasted 6.8 vs 5.69 and yeasty 7.0 vs 4.17 do not (different panel subset, per-vial
    preparation in a).

## 6. Existing records in the repo that cite this paper

`grep -rIn` over `data/`, `docs/`, `src/`, `scripts/`, `tasks/` for `10.1016/j.foodres.2007.08.003`,
`S0963996907001378`, `moon2007`, `Moon 2007`, `Li-Chan`, `simulated beef` (2026-11-02). Not in `data/keys/papers.yml`
(no key for either Moon paper; not edited). No ledger, registry or code record cites it. Two mentions, both in
`docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md`:

**l. 133 (V4 table row)**

| field | repo value | paper | verdict |
|---|---|---|---|
| paper / title | "Moon 2007", "Changes in aroma characteristics of simulated beef flavour by soy protein isolate assessed by descriptive sensory analysis and GC" | same (title ends "gas chromatography") | **matches** (abbreviated) |
| venue / DOI | Food Res Int, 10.1016/j.foodres.2007.08.003 | Food Res Int 40:1239, same DOI | **matches** |
| what it gives: "headspace loss of a simulated beef flavour ... on adding SPI" | — | this paper prints **no headspace value**; it prints sensory loss (Table 7) and correlations. The headspace loss (SBF + SPI as % of SBF alone) is in the **companion, 10.1016/j.foodres.2007.08.002**, Tables 1-2 | **mismatch: wrong DOI for the headspace data** (should cite .002, or both) |
| "(likely includes sulfur volatiles)" | — | **no sulfur compound named in either paper**; 12 named peaks are furans, terpenes, alkylpyrazines, decanal, limonene, 4-terpineol | **not supported**: speculation; the sulfur question rests on Moon 2006 (not on disk) |
| "the only meat-flavour-on-plant-protein series found" | — | true of the search at the time; the corpus now also holds Xu 2022 (pea vs 2-methyl-3-furanthiol, 2-methylpyrazine, 2,4-decadienal, a thiazole; `xu2022_extraction.md`) | **outdated** (search-scope claim, not a paper claim) |
| context quote "Primarily the 3 2 full factorial design was used" | — | **not in this paper or in a**; this paper's design is 2 × 3 × 2 (SBF × SPI × heat), a's is 6 treatments × 2 × 2 | **mismatch** (likely a citing paper's sentence attached to the wrong row; ~80 %) |
| class / verdict | B / MAYBE | — | repo judgement; given sec. 1, "B" (background) is right for this DOI and a is the data paper |

**l. 226** ("Next in line (MAYBE ...)": "..., Bredie 1998, Moon 2007, Reynolds 1959") — a queue entry, no claim to check;
both Moon papers are now on disk and read.

## 7. Proposed ledger entries (`data/lit/measurement_ledger/matrix.yml`; NOT written)

```yaml
  - id: moon2007b_tab7_sensory_spi_dose_beef_flavour
    dossier: moon2007b_extraction.md
    paper: "Moon & Li-Chan 2007 Food Res Int 40:1239 (10.1016/j.foodres.2007.08.003)"
    lab: "UBC Vancouver (Li-Chan)"
    system: "commercial vegetarian simulated beef flavour (Mastertaste F96.49) 150 or 500 mg x commercial soy protein isolate (Solae, 89.3 % protein) 0 / 250 / 500 mg in 5 g water (~0 / 45 / 89 g/L protein) x unheated or 98 C / 30 min heated together; 6 panellists x 2 replications; roasted, beefy, yeasty, soymilk-like, cereal"
    medium: "distilled deionised water, unbuffered; 15 mL vials held 60 C / 20 min before serving; 15 cm unstructured line scale, reference 150 mg SBF + 250 mg SPI = moderate"
    temperature_c: [60, 98]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20, 30]
    species: [roasted, beefy, yeasty, soymilk-like, cereal]
    kind: ranking
    location: "Table 7 (12 samples x 5 attributes, means with Duncan letters; dossier sec. 4d), Table 6 (r with SBF and SPI dose), Table 2 (ANOVA); Fig. 1 plots Table 7"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "the beefy related notes were substantially suppressed by increasing SPI content"
    flags: [sensory_only_no_headspace_numbers, binding_vs_masking_confounded, ph_not_stated, unbuffered, flavour_composition_unknown, no_sulfur_compound_identified, heat_arm_cooks_flavour_with_protein, two_of_eight_panellists_excluded, batch_aliquoted_samples, "also_relevant: aroma_other"]
    note: "The corpus's only protein-dose series (0 / ~4.2 / ~8.1 % SPI protein) of a meat flavour on a plant protein; sensory, not headspace. Beefy relative to no SPI: low flavour 1.00 / 0.79 / 0.30, high flavour 1.00 / 0.91 / 0.71 (unheated). Directional evidence against a zero protein correction for meat odorants on soy; no constant can be drawn. Headspace numbers of the same system (one SPI level): moon2007a_tab1_tab2_spi_headspace_suppression_beef_flavour. Unused."
  - id: moon2007b_si_gc_peak_area_sensory_correlation_matrix
    dossier: moon2007b_extraction.md
    paper: "Moon & Li-Chan 2007 Food Res Int 40:1239 (10.1016/j.foodres.2007.08.003)"
    lab: "UBC Vancouver (Li-Chan)"
    system: "same 12 samples (2 SBF x 3 SPI x 2 heat), 112 GC peaks"
    medium: "HS-SPME 50/30 um DVB/CAR/PDMS 60 min at 60 C, GC-FID, triplicate means"
    temperature_c: [60]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [60, 60]
    species: [unknown_beef_correlated_peaks, 3-methylfuran, 2-acetylfuran, delta-3-carene, 2-ethyl-3,6-dimethylpyrazine, 2,3-diethyl-5-methylpyrazine, decanal, 2-isoamyl-6-methylpyrazine, delta-elemene, beta-cubebene, calamenene]
    kind: other
    location: "Appendix A supplementary data (correlation matrix of mean GC peak areas vs the five sensory means), online at the DOI; not on disk"
    status: supplement_missing
    quant_class: semi_quantitative
    simulable: false
    condition_quote: ""
    flags: [si_not_on_disk, si_content_unknown_may_hold_only_coefficients, peak_area_only, "condition_not_stated: ph"]
    note: "Only route to a headspace protein-dose series (0 / 250 / 500 mg SPI) for this flavour; worth fetching only if the SI prints per-sample peak areas (P ~25 %). Indicator-peak identities and r (Table 8) are in the dossier and get no entry (not a matrix quantity)."
```

No entry for Table 8 (indicator-peak selection: identities from Moon 2006 and correlations, not a binding, partition or
threshold quantity) or for panel-performance Tables 3-5.

## 8. Supplementary material

**Yes: Appendix A** (p. 1247) — "The correlation matrix of the mean value of the GC peak areas and the mean scores of the
five sensory attributes", online with the article at doi:10.1016/j.foodres.2007.08.003 (ScienceDirect, Appendix A,
mmc file; not on disk). The text cites it as "supporting information available" for the 112-peak correlation analysis
(p. 1245). **Model-relevant only if it carries the per-sample mean peak areas**: then it would give the headspace
dose-response (0 / 250 / 500 mg SPI, two SBF doses, heated / unheated) that neither paper prints. If it holds only r
values, it adds nothing beyond Table 8. Recommendation: fetch once, low priority (P(useful) ~25 %). The more useful
missing paper for the sulfur question is Moon, Cliff & Li-Chan 2006 (Food Res Int 39:294; the GC-MS / GC-O
identification of this SBF), see `moon2007a_extraction.md` sec. 8.
