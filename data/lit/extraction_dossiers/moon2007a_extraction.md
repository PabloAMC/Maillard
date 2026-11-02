# Moon & Li-Chan 2007 (Food Res Int 40:1227) — EXTRACTION (file stem `moon2007a`)
DOI: **10.1016/j.foodres.2007.08.002** (printed on p. 1227 and in the PDF metadata title). Food Research International
**40**(10):**1227-1238**, December 2007; received 23 May 2007, accepted 14 Aug 2007. Publisher page
https://www.sciencedirect.com/science/article/abs/pii/S0963996907001366.
**Crossref check (2026-11-02, `api.crossref.org/works/10.1016/j.foodres.2007.08.002`)**: title "Assessment of added
ingredient effect on interaction of simulated beef flavour and soy protein isolate by gas chromatography, spectroscopy
and descriptive sensory analysis", authors Soo-Yeun Moon, Eunice C.Y. Li-Chan, Food Res Int 40(10):1227-1238, issued
2007-12. **Matches the PDF in every field.** (Crossref JSON was parsed with a one-line `python3 -I`, not the repo env; no
repo code was run.)
**Stem check: `moon2007a` = .002 = pp. 1227-1238; the companion `moon2007b` = .003 = pp. 1239-1248. The files are NOT
swapped.** [!] Logical order is the reverse of the stems: this paper (a) is the *follow-up* — it cites b as "Moon &
Li-Chan, accepted for publication" and takes b's 15 indicator peaks and five sensory attributes as given. Read
`moon2007b_extraction.md` first for the panel, the SBF and the indicator-peak definitions.

### Commercial soy protein isolate (Solae, 89.3 % protein, fully denatured by DSC) ± glucosamine / sucrose / PEG 8000 / ascorbic acid / ascorbic acid + PEG (100 mg each), 250 mg in 5 g 0.05 M Tris-HCl pH 7.4, pre-incubated 20 min at 23 or 60 °C, then 500 mg of a commercial vegetarian simulated beef flavour (SBF) added, HS-SPME (DVB/CAR/PDMS) 60 min at 23 or 60 °C, GC-FID. Output: **relative % peak area vs SBF alone (= 100 %)** for 17 (RT adsorption) and 53 (60 °C adsorption) beef-correlated peaks, 12 of them tentatively named; **no sulfur compound among them**. Plus SPI SH / S-S / surface hydrophobicity, FT-Raman structure, and a 7-panellist descriptive test. **No binding constant, no protein-dose series, no absolute concentration, no pH measured after ascorbic acid.**

**Read 2026-11-02 from `data/articles/moon2007a.pdf`** (12 pp., Elsevier, born-digital, text layer good). Prose from
`pdftotext -layout`. **Tables 1 and 2 (rotated, pp. 1231-1233) read from 200 dpi page renders by eye**
(scratchpad `ingest/moon2007a_work/p-05.png`, `p-06.png`, `p-07rot.png`; p. 1233 rotated 90° for reading); every cell
agrees with the text layer. Bold (significant vs SPIF) cells read from the renders. Tables 3-5 (pp. 1233-1236) from the
text layer, Table 3 also checked on the p. 1233 render. **Fig. 2 (cobweb, p. 1237) digitised from a 400 dpi crop**
(`fig2hi-11.png`), read error ±0.3 cm on the 15 cm line scale. Fig. 1 (Raman spectra) not digitised (Table 4 prints the
band intensities).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Food, Nutrition and Health Program, Faculty of Land and Food Systems, **The University of British
Columbia**, Vancouver; corresponding (and senior) author **Eunice C.Y. Li-Chan**; first author Soo-Yeun Moon (UBC
graduate fellowship). GC was run on an instrument at the British Columbia Institute of Technology (acknowledgements); the
FT-Raman at the UBC Laboratory for Molecular Biophysics. **Not in `data/keys/laboratories.yml`: NEW**, proposed id
**`ubc_li_chan`**, lab string **"UBC Vancouver (Li-Chan)"** (same proposal in `moon2007b_extraction.md`). No other
dossier in `data/lit/extraction_dossiers/` is from this group (grep for Li-Chan / British Columbia / UBC: none). Same-lab
companion: **`moon2007b_extraction.md`** (same SPI lot, same SBF, same SPME fibre and vial, same panel core). The two
earlier papers of the series (Moon & Li-Chan 2004 Food Chem 88:141, SPME method; Moon, Cliff & Li-Chan 2006 Food Res Int
39:294, GC-MS/GC-O of the SBF) are **not on disk**.

## 1. What this paper answers for the model

| model need | what Moon 2007a gives | tag |
|---|---|---|
| headspace suppression of meat-flavour volatiles by a plant protein | Tables 1-2, SPIF column: SBF+SPI peak area as % of SBF alone, 4 conditions, ~4.3 % w/w SPI | [M] |
| reversible binding constant K_g (L/g) per compound | **none printed**; apparent K_g from the SPIF ratio under a single-site assumption (sec. 4f) | [Z] |
| sulfur odorants (MFT, FFT, methional, sulfides) on plant protein | **none identified** among the reported peaks (FID, not sulfur-selective; 41 of 53 peaks "Unknown") | [!] |
| alkylpyrazines on plant protein | 2-ethyl-3,6-dimethylpyrazine, 2,3-diethyl-5-methylpyrazine, 2-isoamyl-6-methylpyrazine (60 °C adsorption only) | [M] |
| aldehydes on plant protein | decanal (60 °C adsorption only); no hexanal or nonanal in the tables | [M] |
| T-dependence of the protein effect | two contact (SPME) temperatures, 23 and 60 °C, and two protein pre-incubation temperatures, 23 and 60 °C | [M] |
| protein-concentration dependence | **none** (one SPI level, 250 mg); the dose series is sensory-only, in `moon2007b` | [!] |
| absolute concentrations / partition coefficients | **none**: relative peak areas only, no standards, no IS | [!] |
| soy protein site densities (free SH, S-S) | Table 3: SH 0.6, S-S 36.7 µmol/g protein (Ellman; turbidity caveat by the authors) | [M] |
| state of the commercial SPI | DSC: no transition 25-180 °C, i.e. already denatured [Q]; Raman: 14.4 % α-helix, 55.2 % β-sheet | [M] |
| sensory consequence | Fig. 2 cobweb, 3 samples (SPIF, AF, APF), 7 panellists × 2 reps | [F] |
| pH during the ascorbic-acid arms | **not measured**; my estimate is pH ~4 (sec. 5, flag 1) | [Z] / [!] |

**Answers (paper-specific questions; see `moon2007b_extraction.md` sec. 1 for the same questions on the companion).**

1. **Compounds, protein, concentrations, conditions, design** (confidence ~95 % on the facts as printed).
   *Flavour*: one commercial **vegetarian beef-type reaction flavour, Mastertaste F96x49** (maltodextrin, autolysed yeast
   extract, natural flavours, onion powder, SiO2), 500 mg per vial for GC; individual odorant contents unknown. The
   tables report 17 (RT adsorption) and 53 (60 °C adsorption) peaks that b found positively correlated with the "beefy"
   score; **12 are tentatively named** (3-methylfuran, 2-acetylfuran, limonene, δ-3-carene, 2-ethyl-3,6-dimethylpyrazine,
   2,3-diethyl-5-methylpyrazine, 4-terpineol, decanal, 2-isoamyl-6-methylpyrazine, δ-elemene, β-cubebene, calamenene),
   identifications carried over from Moon 2006 by retention time [C]. *Protein*: **one commercial SPI only** (Solae, lot
   02060631-532, 89.3 ± 0.1 % protein, N × 6.25); **no glycinin / β-conglycinin fractions, no other plant protein**.
   *Concentration*: 250 mg SPI (223 mg protein [Z]) in 5 g buffer + 500 mg SBF ≈ **3.9 % w/w protein of the vial
   contents, ~43 g/L on the aqueous phase [Z]**. *Heating*: the SPI is **already denatured as supplied** (DSC flat
   [Q]); the only heat step is a 20 min pre-incubation of SPI (± additive) at 23 or 60 °C **before** SBF is added.
   *pH*: 0.05 M Tris-HCl **pH 7.4 as prepared; never measured after the additives** (100 mg ascorbic acid in 5 mL of
   this buffer acidifies it, flag 1). *T / time*: flavour-protein contact = HS-SPME exposure, **60 min at 23 or 60 °C**
   (water bath ±2 °C), stirred; no separate equilibration before the fibre goes in. *Design*: **not a 3² factorial**.
   It is 6 treatments (SPIF control + 5 additives) × 2 pre-incubation T × 2 adsorption T, plus SBF alone per adsorption
   T; triplicate independent samples. The "3 2 full factorial" phrase in `LITERATURE_SEARCH_2026-10_VARIABLES.md` is not
   in either Moon paper (sec. 6).
2. **Quantities** (confidence ~95 %). **Relative % peak area, HS-SPME / GC-FID, each peak vs the same peak in SBF
   alone (100 %)**; SBF-alone absolute areas (FID counts, no units) printed. Values are **relative, not absolute; not
   static headspace** (SPME on a carboxen-containing fibre, 60 min, stirred; competitive and possibly non-equilibrium
   [!]). No % bound, no binding constant, no K_aw. Triplicate CVs 10.1 % (SBF alone) and 10.6 % (SBF + SPI) on peak areas
   [M]; no SD per cell; **no significance test of SPIF vs SBF** (bold marks only additive vs SPIF). **Zeros are
   non-detections** (no LOD stated) [!]. The SPI effect at ~43 g/L: most beef-correlated peaks fall to **~40-60 % of
   SBF alone**; limonene to 17-24 %; 3-methylfuran and 2-acetylfuran stay at 71-122 % (sec. 4b).
3. **Sensory** (confidence ~85 %; numbers [F] ±0.3 cm). 7 panellists (5 from b + 2 new) × 2 replications, 150 mg SBF +
   250 mg SPI in 5 g water, held 60 °C / 20 min. **APF (SPI pre-treated with ascorbic acid + PEG) vs SPIF**: roasted
   **9.4 vs 6.8** (significant), soymilk-like **3.7 vs 6.9** and cereal **4.2 vs 6.4** (both significant); beefy 8.8 vs
   6.7 and yeasty 6.0 vs 7.0 (not significant, p > 0.05). AF (ascorbic acid alone) moved the same way and was not
   significantly different from SPIF [Q]. The sensory test used water, not the Tris buffer of the GC, and 150 not 500 mg
   SBF [!].
4. **Usable data sets and the dependences the matrix layer lacks** (confidence ~80 %). Three: (i) Tables 1-2 relative
   areas (the SPI ratio and the additive effects), (ii) Table 3 SH / S-S / S0, (iii) Fig. 2 sensory. Sec. 7 gives
   YAML for each. **Protein-concentration dependence: none here** (one level). **Temperature dependence: weak, two
   contact temperatures**: for the large unknown peaks and δ-elemene the SPI ratio is essentially the same at 23 and
   60 °C (e.g. 31.881 min: 49 → 50 %; δ-elemene 42 → 48 %), while 3-methylfuran's suppression vanishes at 60 °C (71 → 106 %)
   (RT pre-incubation; sec. 4c). **60 °C is the boundary of `BINDING_AT_PROCESS_TEMPERATURE`, not beyond it**: the
   "no aqueous constant above 60 °C" limit stands (confidence ~90 %).

**What the matrix code assumes for meat odorants on soy, and what this paper says** (registry =
`src/kinetic_core/parameters_matrix.py` `REVERSIBLE_BINDING`, `SULFUR_BINDING_GAP`, `NO_ADDUCT_GATE_*`,
`ADDUCT_POSITIVE_CLASSES`, `BINDING_AT_PROCESS_TEMPERATURE`; `data/species/protein_matrices.yml` `soy_isolate`):

| odorant / class | what the code assumes (soy) | Moon 2007a | verdict |
|---|---|---|---|
| 2-methyl-3-furanthiol, 2-furfurylthiol (thiol) | no reversible constant; `SULFUR_BINDING_GAP` "measured once, not shippable" (pea DMDS/DMTS); thiol in `ADDUCT_POSITIVE_CLASSES`, covalent ceiling contributes zero | no thiol identified; 41 of 53 peaks unknown on FID | **silent** (~95 %) |
| methional | not in `COMPOUND_STRUCTURE`; no constant | not identified | **silent** |
| alkylpyrazines | no alkylpyrazine binding row (Hong 2020 2,3-dimethylpyrazine → `no_binding_constant_for_class`); class in `NO_ADDUCT_GATE_CLASSES`; only pyrazine K is Xu 2022 pea 2-methylpyrazine (~4.3e-2 L/g native, 25 °C, directional only) | 2-ethyl-3,6-dimethylpyrazine at 54 % (RT pre-inc.) / 64 % (60 °C pre-inc.) of SBF alone, 60 °C contact; apparent K_g 0.013-0.020 L/g [Z]; 2,3-diethyl-5-methylpyrazine 0 / 46 % | **contradicts the implicit zero correction** (a pyrazine is suppressed ~1.6-1.9× at 4 % SPI; ~70 %); same order as, 2-3× below, Xu's pea value; supports the no-adduct (reversible) classing only by default, since nothing tests reversibility |
| n-aldehydes (soy) | `kg_nonanal_soy` 4.38e-2 L/g (dialysis, pH 8, 25 °C); `kg_hexanal_soy_denatured` 1.47e-3 (gel filtration); chain-length slope 2.81×/CH2 (decanal by extrapolation 0.123 L/g [Z]); aldehyde method families never pooled | decanal 54-57 % at 60 °C → apparent K_g 0.018-0.020 L/g [Z], **~6-7× below** the dialysis-family extrapolation | **neither** (different method family by the code's own rule); the first headspace-family decanal-on-soy reading in the corpus, at the weak end (~60 % that a direct comparison would hold up after method / T / denaturation are matched) |
| binding vs temperature | "no scalar; three signs"; nothing aqueous above 60 °C | 23 → 60 °C contact: ratio flat for most peaks, falls in suppression for 3-methylfuran and calamenene detection | **supports "no scalar"** weakly; does not lift the 60 °C limit |
| soy_isolate free SH / S-S (native band 4.7-8.0 / 38-53 µmol/g) | sites charged once at the native density; heated_reference 0.75 µmol/g SH after 100 °C | **SH 0.6, S-S 36.7 µmol/g protein** on a commercial, DSC-denatured SPI | **contradicts the native band for a commercial SPI** (8-13× lower SH), matches the heated_reference; S-S just below band (~75 %) |

## 2. Conditions

- **Flavour** [M]: SBF = "vegetarian beef type flavour F96x49" (Mastertaste, Arlington Heights IL): maltodextrin,
  autolysed yeast extract, natural flavours, onion powder, silicon dioxide. **500 mg per GC vial; 150 mg per sensory
  sample.** Volatile composition from Moon 2006 [C] (≥70 compounds by GC-MS, 49 odour-active by GC-O, per `moon2007b`).
  [!] The "SBF alone" control already carries yeast-extract protein / peptides and maltodextrin: the reference is not a
  protein-free water leg.
- **Protein** [M]: soy protein isolate, Solae (St. Louis MO), lot 02060631-532, **89.3 ± 0.1 % protein** (LECO FP-428
  N combustion, N × 6.25). Commercial; DSC thermograms at 0.2 / 5 / 10 % SPI in pH 7.5 phosphate show **no endo- or
  exotherm 25-180 °C** ("data not shown") → "already considerably denatured" [Q]. No 7S / 11S fractions.
- **Additives** [M]: glucosamine (≥99 %), sucrose, ascorbic acid (≥99 %), polyethylene glycol (Mw ~8000), all Sigma;
  **100 mg each** per GC vial (AP = 100 mg AA for 20 min, then 100 mg PEG for another 20 min).
- **GC vial recipe** [M]: 100 mg additive + 250 mg SPI in a 15 mL Supelco vial with stir bar → + **5 g 0.05 M Tris-HCl
  pH 7.4** → incubate **20 min at RT (23 °C) or 60 °C** → + 500 mg SBF → cap (PTFE/silicone septum) → SPME. Controls: SPI
  + SBF without additive (SPIF), and SBF alone (no SPI). Headspace ~9-10 mL [Z: 15 mL vial minus ~5.5 mL].
  [!] Whether the 60 °C-incubated vials were cooled before SBF addition and RT adsorption is not stated.
- **Protein concentration** [Z]: 0.25 × 0.893 = 0.223 g protein in ~5.0-5.3 mL aqueous phase → **42-45 g/L** (I use
  43 g/L); 3.9 % w/w of the 5.75 g vial contents.
- **pH** [M]: 7.4 as prepared. **Not measured after additives or SBF.** [Z] 100 mg ascorbic acid = 0.57 mmol in 5 mL
  (0.11 M) against 0.05 M Tris at pH 7.4 (≈ 8 mM free base): the buffer is overwhelmed; with ~0.2 mmol of acid-binding
  by 223 mg protein (rough), the pH lands near **4** (the ascorbic acid pKa1 is 4.1), i.e. at or below the soy pI (~4.5)
  (confidence ~80 % that pH < 5). The authors treat ionic effects as negligible at "pH 7.4" (p. 1230) [Q].
- **Contact T / time** [M]: SPME exposure **60 min at 23 °C or 60 °C** (thermostatted water bath ±2 °C), stirred
  throughout. SBF-protein contact before the fibre is inserted: not stated, presumably minutes.
- **Sensory samples** [M]: 250 mg SPI (± 100 mg AA ± 100 mg PEG pre-treatment) in 5 g dd-**water** + **150 mg SBF**
  added per vial (not batch-aliquoted, to limit volatile loss, unlike b), foil-covered, **60 °C / 20 min** before
  serving. Reference R = 250 mg SPI + 150 mg SBF in 5 g dd-water, anchored "moderate".
- **Replicates** [M]: GC, three independently prepared samples per cell; SH / S-S duplicate; S0 triplicate; Raman
  duplicate (CV ≤1 %, significance from an assumed 5 % CV); sensory 7 panellists × 2 replications (n = 14).
- **Water activity / moisture**: dilute aqueous throughout (~87 % water [Z]); Raman on freeze-dried samples.

## 3. Analytical methods

- **Headspace volatiles** [M]: HS-SPME, 50/30 µm StableFlex **DVB/CAR/PDMS**, 60 min at 23 or 60 °C, then GC-FID "as
  previously described" (Moon & Li-Chan 2004 [C]; column, oven programme and desorption **not restated here** [!]).
  **Quantification: relative peak area** (each peak / same peak in SBF alone × 100). **No external standard, no internal
  standard, no LOD.** Identities are **tentative, by retention time match to Moon 2006 GC-MS** [C]; the FID cannot see
  sulfur selectively. DF = GC-O detection frequency from Moon 2006; "Beefy" = Pearson r with the beefy score from b.
  Peak-area CV 10.1 / 10.6 % (SBF / SBF + SPI triplicates). The repo's same-sample HS-SPME dispersion
  (`HS_SPME_SAME_SAMPLE_DISPERSION` 10-23 %) brackets these at its low end.
- **Statistics (GC)** [M]: GLM ANOVA, Fisher's LSD, 95 %; bold = additive differs from SPIF at p ≤ 0.05. **The ± on
  relative areas is not printed; SPIF vs SBF is not tested.**
- **SH / S-S** [M]: Ellman's DTNB (Beveridge 1974 as modified by Chung 2005) in pH 8.0 Tris-glycine-EDTA-SDS (10 mg/mL);
  total SH + SS by Thannhauser 1984 (NTSB, dark incubation); SS = (SH + SS) − SH, by difference [M];
  units µmol per g protein; ingredient blanks subtracted. Authors warn of **high protein-blank turbidity** even after
  14 950 g / 30 min with SDS [Q].
- **Surface hydrophobicity S0** [M]: PRODAN probe (Alizadeh & Li-Chan 2000), 1 % SPI stock in 0.05 M Tris pH 7.4 with
  40 mg additive per 10 mL, 60 °C / 20 min, diluted to 0.0025-0.010 % SPI; S0 = initial slope of net RFI vs % SPI
  (units %⁻¹); ± = SD, n = 3.
- **FT-Raman** [M]: Nexus 670, 0.5 W, 512 scans, 4 cm⁻¹, freeze-dried samples; normalised to Phe 1003 cm⁻¹; secondary
  structure by RSAP least squares on amide I.
- **Sensory** [M]: unstructured 15 cm line, anchors 1.3 cm from each end ("Slight" / "Moderate" / "Intense"), score in
  cm from origin.

## 4. Data

### 4a. Table 1 (p. 1231): relative % area, **adsorption at RT (23 °C)** [M]

SBF = absolute area of SBF alone (FID counts); all other columns % of it. Bold in the paper = significant vs SPIF
(p ≤ 0.05), marked here with **. DF = GC-O detection frequency (%), Beefy = Pearson r with beefy score (from b).

| IP | RT (min) | SBF area | RT-inc: SPIF | GF | SF | PF | AF | APF | 60 °C-inc: SPIF | GF | SF | PF | AF | APF | identification | DF | beefy r |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IP1 | 1.734 | 11512 | 71 | 94 | 68 | 60 | 88 | 48 | 42 | 42 | 40 | 32 | 55 | 56 | 3-methylfuran | 50 | 0.710* |
| | 13.271 | 2514 | 19 | 15 | 0 | 0 | 0 | 0 | 17 | 0 | 0 | 0 | 0 | 0 | limonene | | 0.713** |
| IP13 | 28.288 | 13837 | 42 | **72** | 50 | 50 | 63 | 64 | 50 | 45 | 42 | 34 | 57 | 57 | δ-elemene | 50 | 0.931*** |
| | 29.122 | 933 | 129 | **387** | 159 | 104 | 58 | 133 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.713** |
| | 29.971 | 15016 | 48 | **73** | 57 | 58 | **71** | **74** | 49 | 45 | 42 | 35 | 59 | 59 | unknown | | 0.928*** |
| | 30.729 | 1897 | 0 | **57** | **56** | 42 | **75** | **77** | 22 | 39 | 19 | 21 | **71** | **67** | unknown | | 0.929*** |
| | 31.271 | 4513 | 57 | 89 | 65 | 74 | **288** | **184** | 66 | 64 | 55 | 43 | **110** | 78 | unknown | | 0.638* |
| | 31.530 | 2735 | 24 | 83 | 29 | 50 | 70 | 66 | 49 | 53 | 0 | 26 | 87 | 78 | unknown | | 0.750** |
| | 31.881 | 113227 | 49 | **78** | 58 | 61 | 63 | **76** | 54 | 52 | 46 | 41 | 67 | 68 | unknown | | 0.924*** |
| | 32.231 | 3743 | 53 | 79 | 59 | 61 | 70 | 77 | 61 | 52 | 53 | 38 | 73 | 75 | unknown | | 0.879*** |
| | 33.166 | 3793 | 50 | **75** | 54 | 61 | 68 | 72 | 58 | 46 | 50 | 37 | 67 | 70 | unknown | | 0.826** |
| | 33.306 | 6149 | 50 | **79** | 58 | 63 | **73** | **83** | 56 | 55 | 49 | 45 | 71 | 72 | unknown | | 0.894*** |
| | 34.691 | 828 | 0 | 52 | 0 | 0 | 50 | 106 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.881*** |
| | 35.365 | 1501 | 57 | 80 | 0 | 0 | 26 | 57 | 0 | 0 | 0 | 0 | 0 | 29 | unknown | | 0.823** |
| IP15 | 35.759 | 417 | 0 | 105 | 0 | 0 | 98 | 220 | 0 | 0 | 0 | 0 | 85 | 90 | calamenene | 75 | 0.883*** |
| | 35.887 | 1624 | 0 | 40 | 32 | 215 | 58 | **393** | 81 | 181 | 77 | 147 | 82 | **296** | unknown | | 0.685* |
| | 36.322 | 5318 | 46 | 81 | 55 | 60 | 68 | 82 | 51 | 48 | 44 | 43 | 66 | 67 | unknown | | 0.822** |

SPIF = SBF + SPI; GF / SF / PF / AF / APF = SBF + SPI pre-incubated with glucosamine / sucrose / PEG / ascorbic acid /
ascorbic acid + PEG. Peak areas = mean of three replicate GC analyses.

### 4b. Table 2 (pp. 1232-1233): relative % area, **adsorption at 60 °C** [M]

| IP | RT (min) | SBF area | RT-inc: SPIF | GF | SF | PF | AF | APF | 60 °C-inc: SPIF | GF | SF | PF | AF | APF | identification | DF | beefy r |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IP1 | 1.734 | 6953 | 106 | 91 | 113 | 94 | **163** | **158** | 80 | 83 | 85 | 78 | **131** | **146** | 3-methylfuran | 50 | 0.710* |
| IP2 | 7.893 | 5537 | 122 | 83 | 141 | 89 | 83 | 74 | 77 | 62 | 75 | 65 | 59 | 63 | 2-acetylfuran | 75 | 0.658* |
| IP3 | 12.132 | 3360 | 66 | 58 | 67 | 54 | 62 | 49 | 20 | 31 | 32 | 17 | **43** | **45** | δ-3-carene | 88 | 0.872*** |
| | 13.271 | 1929 | 24 | 0 | 30 | 42 | 0 | 0 | 18 | 0 | 0 | 0 | 0 | 20 | limonene | | 0.713** |
| | 13.915 | 1352 | 33 | 34 | 39 | 31 | **240** | **174** | 0 | **82** | 25 | 0 | **133** | **148** | unknown | | 0.648* |
| | 15.014 | 6302 | 84 | 63 | 62 | 27 | 89 | 73 | 67 | 73 | 68 | 67 | 60 | 67 | unknown | | 0.773** |
| IP4 | 15.877 | 1761 | 54 | 51 | 88 | 80 | 87 | 47 | 64 | 50 | 53 | 52 | 48 | 51 | 2-ethyl-3,6-dimethylpyrazine | 88 | 0.601* |
| IP5 | 16.247 | 2543 | 49 | 23 | 68 | 75 | 82 | 49 | 45 | 35 | 36 | 32 | 43 | 49 | unknown | 100 | 0.852*** |
| IP6 | 16.953 | 1935 | 20 | 0 | 49 | 39 | 38 | 0 | 20 | 18 | 0 | 0 | 18 | 18 | unknown | 100 | 0.827** |
| IP8 | 19.917 | 2881 | 0 | 23 | 35 | 28 | 62 | 0 | 46 | 0 | 53 | 0 | 0 | 37 | 2,3-diethyl-5-methylpyrazine | 50 | 0.590* |
| | 20.608 | 3328 | 47 | 44 | 61 | 59 | **99** | **81** | 49 | 46 | 45 | 48 | **75** | 71 | 4-terpineol ("4-Terpeneol") | | 0.678* |
| IP9 | 22.151 | 3951 | 57 | 52 | 59 | 58 | 77 | 64 | 54 | 55 | 50 | 47 | 64 | **75** | decanal | 75 | 0.715** |
| IP10 | 24.232 | 468 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2-isoamyl-6-methylpyrazine | 50 | 0.599* |
| IP12 | 25.794 | 8033 | 41 | 47 | 46 | **51** | **50** | 40 | 38 | 42 | 45 | 50 | **71** | 43 | unknown | 75 | 0.794** |
| | 26.319 | 3792 | 19 | 18 | 40 | 41 | 38 | 34 | 0 | 0 | 17 | 0 | 16 | 19 | unknown | | 0.684* |
| | 26.954 | 2156 | 0 | 19 | 21 | **99** | **63** | 19 | 24 | 24 | 0 | 0 | 0 | 43 | unknown | | 0.752** |
| | 27.841 | 3710 | 16 | 35 | 38 | 56 | 22 | 95 | 40 | 22 | 35 | 21 | 0 | 0 | unknown | | 0.887*** |
| IP13 | 28.288 | 19545 | 48 | 36 | 55 | 55 | 62 | 42 | 50 | 42 | 47 | 42 | 46 | 46 | δ-elemene | 50 | 0.931*** |
| | 28.830 | 8259 | 44 | 41 | 53 | 61 | 35 | 28 | 43 | 41 | 47 | 50 | 25 | 26 | unknown | | 0.678* |
| | 29.122 | 7556 | 131 | 130 | 130 | 128 | 158 | 109 | 116 | 149 | 117 | 151 | 115 | 123 | unknown | | 0.713** |
| | 29.971 | 25200 | 46 | 39 | 52 | 53 | 62 | 54 | 51 | 45 | 49 | 44 | 50 | 53 | unknown | | 0.928*** |
| | 30.220 | 1260 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.617* |
| IP14 | 30.377 | 3031 | 16 | 14 | 49 | **50** | **62** | **60** | 40 | 17 | 15 | 19 | 36 | 36 | β-cubebene | 63 | 0.850*** |
| | 30.615 | 1639 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.645* |
| | 30.729 | 9478 | 34 | 29 | 42 | 39 | 46 | 42 | 40 | 37 | 40 | 35 | 37 | 43 | unknown | | 0.929*** |
| | 30.950 | 721 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.868*** |
| | 31.271 | 1849 | 88 | 0 | 106 | 102 | 0 | 0 | 34 | 0 | 46 | 41 | 0 | 0 | unknown | | 0.638* |
| | 31.317 | 11730 | 54 | 71 | 60 | 66 | **158** | **131** | 69 | 77 | 67 | 57 | 94 | **104** | unknown | | 0.914*** |
| | 31.530 | 4843 | 44 | 50 | 39 | 56 | 98 | 76 | 50 | 74 | 72 | 82 | 74 | 73 | unknown | | 0.750** |
| | 31.881 | 262026 | 50 | 43 | 57 | 57 | 55 | 50 | 55 | 47 | 54 | 49 | 48 | 53 | unknown | | 0.924*** |
| | 32.231 | 7117 | 41 | 25 | 43 | 44 | 53 | 47 | 42 | 38 | 40 | 33 | 54 | 48 | unknown | | 0.879*** |
| | 32.679 | 4244 | 42 | 29 | 48 | 44 | 30 | 30 | 42 | 39 | 40 | 36 | 40 | 48 | unknown | | 0.863*** |
| | 33.166 | 6537 | 43 | 38 | 44 | 45 | 56 | 53 | 51 | 45 | 45 | 38 | 56 | 57 | unknown | | 0.826** |
| | 33.306 | 24572 | 50 | 45 | 57 | 57 | 55 | 51 | 54 | 48 | 53 | 49 | 48 | 52 | unknown | | 0.894*** |
| | 33.620 | 4026 | 9 | 13 | 12 | 25 | 25 | 14 | 29 | 9 | 11 | 0 | 24 | 34 | unknown | | 0.779** |
| | 33.813 | 4478 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.689* |
| | 34.211 | 7156 | 31 | 17 | 31 | 49 | 26 | 30 | 32 | 26 | 33 | 37 | 15 | 24 | unknown | | 0.775** |
| | 34.346 | 11082 | 30 | 26 | 31 | 35 | 31 | 28 | 36 | 32 | 35 | 32 | 22 | 26 | unknown | | 0.787** |
| | 34.481 | 2104 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 62 | 47 | 0 | 0 | unknown | | 0.600* |
| | 34.691 | 9569 | 40 | 34 | 44 | 46 | 45 | 43 | 43 | 38 | 29 | 26 | 39 | 41 | unknown | | 0.881*** |
| | 34.871 | 4778 | 33 | 9 | 11 | 10 | 11 | 9 | 20 | 18 | 18 | 19 | 8 | 0 | unknown | | 0.781** |
| | 35.097 | 9752 | 27 | 35 | 43 | 45 | 40 | 38 | 38 | 36 | 38 | 36 | 38 | 42 | unknown | | 0.872*** |
| | 35.365 | 15997 | 34 | 31 | 34 | 37 | 37 | 36 | 36 | 33 | 33 | 31 | 36 | 40 | unknown | | 0.823** |
| IP15 | 35.759 | 16917 | 39 | 35 | 41 | 38 | 47 | 39 | 38 | 37 | 39 | 32 | 41 | 41 | calamenene | 75 | 0.883*** |
| | 35.887 | 20361 | 39 | 42 | 40 | **404** | 55 | **407** | 43 | 48 | 40 | **140** | 46 | **379** | unknown | | 0.685* |
| | 36.322 | 49502 | 36 | 33 | 40 | 40 | 46 | 43 | 38 | 35 | 37 | 34 | 41 | 43 | unknown | | 0.822** |
| | 36.646 | 11216 | 37 | 26 | 40 | 40 | 30 | 28 | 41 | 38 | 41 | 37 | 29 | 32 | unknown | | 0.780** |
| | 37.033 | 3176 | 0 | 61 | 26 | 54 | 64 | 54 | 59 | 27 | 69 | 23 | 41 | 43 | unknown | | 0.800** |
| | 37.347 | 12581 | 32 | 33 | 38 | 40 | 38 | 30 | 38 | 33 | 39 | 37 | 27 | 31 | unknown | | 0.778** |
| | 37.674 | 2119 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | unknown | | 0.815** |
| | 38.001 | 12714 | 45 | 39 | 47 | 47 | 31 | 26 | 38 | 31 | 51 | 53 | 35 | 37 | unknown | | 0.746** |
| | 38.425 | 6531 | 44 | 36 | 43 | 48 | 42 | 41 | 40 | 43 | 46 | 42 | 37 | 33 | unknown | | 0.634* |
| | 38.658 | 13153 | 40 | 36 | 43 | 52 | 39 | 38 | 39 | 36 | 44 | 44 | 34 | 30 | unknown | | 0.767** |

Text p. 1233 [M]: among the 15 indicator peaks, **2 / 4 / 4** were significantly increased by PEG / AA / AA + PEG
(Table 2). **Checked** [Z]: PF bold on IP12, IP14 (2); AF bold on IP1, IP3, IP12, IP14 (4); APF bold on IP1, IP3, IP9,
IP14 (4). Matches.

### 4c. The SPI effect alone (SPIF column) for the named compounds and the T contrasts [M, rearranged]

| compound | RT-inc / RT-ads | 60-inc / RT-ads | RT-inc / 60-ads | 60-inc / 60-ads |
|---|---|---|---|---|
| 3-methylfuran | 71 | 42 | 106 | 80 |
| 2-acetylfuran | – | – | 122 | 77 |
| limonene | 19 | 17 | 24 | 18 |
| δ-3-carene | – | – | 66 | 20 |
| 2-ethyl-3,6-dimethylpyrazine | – | – | 54 | 64 |
| 2,3-diethyl-5-methylpyrazine | – | – | 0 | 46 |
| 4-terpineol | – | – | 47 | 49 |
| decanal | – | – | 57 | 54 |
| 2-isoamyl-6-methylpyrazine | – | – | 0 | 0 (SBF area only 468) |
| δ-elemene | 42 | 50 | 48 | 50 |
| β-cubebene | – | – | 16 | 40 |
| calamenene | 0 (SBF area 417) | 0 | 39 | 38 |
| largest peak, 31.881 min (unknown) | 49 | 54 | 50 | 55 |

"–" = not in Table 1 (RT adsorption lists only 17 peaks). [Z] Median SPIF over the 53 Table 2 rows: ~40 % (RT-inc) and
~40 % (60 °C-inc); 6 rows (24.232, 30.220, 30.615, 30.950, 33.813, 37.674 min) are 0 in every column (below detection with SPI). **Contact temperature 23 → 60 °C** (RT
pre-incubation): flat within the ~15 % ratio noise for 31.881 / 29.971 / 33.306 / δ-elemene / limonene; up for
3-methylfuran (71 → 106). **Protein pre-incubation 23 → 60 °C** (60 °C contact): δ-3-carene 66 → 20, β-cubebene 16 → 40,
2-acetylfuran 122 → 77, 3-methylfuran 106 → 80, pyrazine IP4 54 → 64: mixed signs, as `BINDING_AT_PROCESS_TEMPERATURE`
records for other systems.

### 4d. Table 3 (p. 1233): SH, S-S, surface hydrophobicity [M]

| sample | total SH (µmol SH/g protein) | disulfide (µmol SS/g protein) | S0 (%⁻¹), mean ± SD |
|---|---|---|---|
| SPI | 0.6 b | 36.7 b | 145 ± 4 a |
| G (glucosamine) | 0.6 b | 36.8 b | 146 ± 2 a |
| S (sucrose) | 0.8 c | 38.6 b | 150 ± 2 a |
| P (PEG) | 0.5 a | 36.8 b | 170 ± 9 b |
| A (ascorbic acid) | 2.0 d | 31.0 a | 171 ± 6 b |
| AP (AA + PEG) | 2.4 e | 30.8 a | 198 ± 24 c |

SH / SS n = 2 (no ±); S0 n = 3. Text compares Boatright & Hettiarachchy 1995 (SH 8.3, SH + SS 52.9 µmol/g) and
Puppo 2004 (SH 0.29 at pH 8, 5.4 at pH 3) [C]. [Z] AA raises SH by 1.4 and lowers SS by 5.7 µmol/g: 2 × 5.7 = 11.4 µmol
of thiol should appear per gram and only 1.4 does; the authors note the shortfall [Q].

### 4e. Tables 4-5 (pp. 1235-1236): FT-Raman, one line each [M]

Table 4: S-S stretch shifts 522 → 537-541 cm⁻¹ with every additive (g-g-t → t-g-t); S-S intensity 0.55 (SPI) → 0.48 (A,
AP, significant); Trp 760 cm⁻¹ 0.74 → 0.58 (AP only); Tyr 850/830 ratio 0.97 → 0.60 (A) / 0.85 (AP); C-H bend 1448-1450
1.62 → 1.80 (P) / 1.81 (AP). Table 5 (% secondary structure, NOFD / SPI / G / S / P / A / AP): α-helix 14.4 / 13.6 / 13.1 /
15.0 / 16.5 / 18.9 / 17.5; β-sheet 55.2 / 57.2 / 59.6 / 57.2 / 55.1 / 56.8 / 56.9; β-turn 19.0 / 18.6 / 20.0 / 19.8 / 17.0 /
**1.4** / 7.6; unordered 11.5 / 10.6 / 7.3 / 8.0 / 11.4 / 22.9 / 18.0. Protein structure only; no model consumer.

### 4f. Apparent single-site constants from the SPIF ratio [Z]

If the SPI effect were reversible single-site binding with the fibre sampling the free fraction, R = 1/(1 + K_g·c_p), so
K_g = (1/R − 1)/c_p with c_p ≈ 43 g/L (band 42-45). **Assumption-laden** (SPME is competitive and may not be at
equilibrium; the SPI brings its own volatiles onto the carboxen; the SBF-alone leg is not protein-free; ±~15 % on each
ratio from the 10 % CVs).

| compound (60 °C contact) | R (RT-inc / 60-inc) | apparent K_g, L/g |
|---|---|---|
| 2-ethyl-3,6-dimethylpyrazine | 0.54 / 0.64 | 0.020 / 0.013 |
| 2,3-diethyl-5-methylpyrazine | 0 / 0.46 | not computable / 0.027 |
| decanal | 0.57 / 0.54 | 0.018 / 0.020 |
| 4-terpineol | 0.47 / 0.49 | 0.026 / 0.024 |
| δ-elemene | 0.48 / 0.50 | 0.025 / 0.023 |
| calamenene | 0.39 / 0.38 | 0.036 / 0.038 |
| limonene | 0.24 / 0.18 | 0.073 / 0.106 |
| δ-3-carene | 0.66 / 0.20 | 0.012 / 0.093 |
| 3-methylfuran, 2-acetylfuran | 1.06 / 0.80; 1.22 / 0.77 | ≤0 / 0.006; ≤0 / 0.007 |

Comparisons: registry `kg_nonanal_soy` 4.38e-2 × 2.81 (one CH2) → decanal **0.123 L/g** (dialysis family, 25 °C, pH 8),
predicting R = 0.16 at 43 g/L vs the 0.54-0.57 observed: **6-7× weaker binding here**. Xu 2022 2-methylpyrazine on native
pea, ~4.3e-2 L/g at 25 °C (`xu2022_extraction.md`): the Moon pyrazine values are **2-3× lower**, on a more alkylated
pyrazine, at 60 °C, on a denatured protein.

### 4g. Fig. 2 (p. 1237): cobweb, sensory means (cm on 15 cm line), n = 14 [F, ±0.3]

| attribute | SPIF | AF | APF | significance (text) |
|---|---|---|---|---|
| roasted | 6.8 | 7.9 | 9.4 | APF > SPIF, p ≤ 0.05 |
| beefy | 6.7 | 8.0 | 8.8 | n.s. |
| yeasty | 7.0 | 6.9 | 6.0 | n.s. |
| soymilk-like | 6.9 | 5.4 | 3.7 | APF < SPIF |
| cereal | 6.4 | 4.9 | 4.2 | APF < SPIF |

AF vs SPIF: no significant difference on any attribute [Q]. Cross-check with b (Table 7, S3 = the same 150 mg SBF +
250 mg SPI recipe): beefy 6.71, roasted 5.69, soymilk 6.73, cereal 6.84, yeasty 4.17 (`moon2007b_extraction.md` sec. 4)
— beefy, soymilk and cereal agree within 0.4 cm; yeasty differs by ~3 cm.

## 5. Flags and internal consistency

1. [!] **Ascorbic-acid arms are probably acidified to pH ~4, unmeasured** (sec. 2, [Z], ~80 %). 100 mg AA in 5 g of
   0.05 M Tris pH 7.4 (GC) and 40 mg in 10 mL (S0 assay) exceed the buffer. Near the soy pI the isolate aggregates and
   its binding changes; the paper attributes the AF / APF recovery of beef peaks to disulfide reduction only. The
   headline mechanism is confounded with pH.
2. [!] **Relative peak areas only; no standards, no IS, no LOD, no SD per cell, SPIF vs SBF never tested.** Zeros are
   censored non-detections (some SBF-alone areas are tiny: calamenene 417, IP10 468 at RT).
3. [!] **HS-SPME on DVB/CAR/PDMS for 60 min is competitive and not shown to be at equilibrium**; adding SPI adds soy
   volatiles to the fibre, which can lower SBF peaks without any binding. The SPIF ratio is "headspace availability",
   not partition.
4. [!] **The "SBF alone" reference is not a protein-free water leg**: it contains autolysed yeast extract and
   maltodextrin. It was also analysed once per adsorption T (no incubation arm).
5. [!] **Ratios above 100 %** (3-methylfuran 106 / 163, 2-acetylfuran 122 / 141, 29.122 min 116-387 %, 35.887 min up to
   407 %): salting-out, formation in the vial, or co-elution of an SPI / additive volatile; no blank of SPI + additive
   without SBF is reported. The 35.887 peak jumps ×4-10 with PEG, suggesting a PEG-derived co-eluter.
6. [!] **Identifications are tentative** (retention-time transfer from Moon 2006 GC-MS); 41 of 53 listed peaks (all 17 Table 1 peaks recur in Table 2) are
   unknown. **No sulfur compound is named**; an FID cannot flag one.
7. [!] **Two of the 15 indicator peaks of b (IP7, 18.199 min; IP11, 25.374 min) are absent from Table 2**, although every
   indicator peak is by definition beef-correlated.
8. [!] Text p. 1233 says "beta-cubebene and an unidentified peak (IP12 and IP14, respectively)"; Table 2 has IP12 =
   unknown and IP14 = β-cubebene — "respectively" is inverted.
9. [!] Text says the glucosamine effect "disappeared" under 60 °C incubation or 60 °C adsorption; Table 2 still has one
   bold GF cell (13.915 min, 60/60, 82 vs 0).
10. [!] **GC and sensory used different media and doses**: GC in Tris pH 7.4 with 500 mg SBF; sensory in water with
    150 mg SBF (sensory samples also not pH-controlled).
11. [!] SH result flagged by the authors themselves (blank turbidity); the SH rise does not balance the S-S fall
    (sec. 4d [Z]).
12. [!] Raman significance from an assumed 5 % CV, not the measured ≤1 % (p. 1230).
13. [!] Citation of b in the reference list ends its title with "gas chromatography–olfactometry"; b's printed title ends
    "gas chromatography" (cosmetic).
14. [!] GC column and oven programme are only referenced (Moon & Li-Chan 2004, not on disk), so retention times cannot
    be checked against an index.

## 6. Existing records in the repo that cite this paper

`grep -rIn` over `data/`, `docs/`, `src/`, `scripts/`, `tasks/` for `10.1016/j.foodres.2007.08.002`, `S0963996907001366`,
`moon2007`, `Moon 2007`, `Li-Chan`, `simulated beef` (2026-11-02): **no record cites this DOI.** Not in
`data/keys/papers.yml` (no key exists for either Moon paper; not edited). The only repo mention of the series is
`docs/guides/LITERATURE_SEARCH_2026-10_VARIABLES.md` l. 133 (and l. 226 "Moon 2007" in the MAYBE queue), which cites the
**companion** DOI .003; it is checked field by field in `moon2007b_extraction.md` sec. 6. Relevant here: that row's
"headspace loss of a simulated beef flavour ... on adding SPI" is **in this paper (.002), Tables 1-2, not in .003**, which
prints no headspace values.

## 7. Proposed ledger entries (`data/lit/measurement_ledger/matrix.yml`; NOT written)

```yaml
  - id: moon2007a_tab1_tab2_spi_headspace_suppression_beef_flavour
    dossier: moon2007a_extraction.md
    paper: "Moon & Li-Chan 2007 Food Res Int 40:1227 (10.1016/j.foodres.2007.08.002)"
    lab: "UBC Vancouver (Li-Chan)"
    system: "commercial vegetarian simulated beef flavour (Mastertaste F96x49, 500 mg) + commercial soy protein isolate (Solae, 89.3 % protein, DSC-denatured; 250 mg = ~43 g/L protein) +- 100 mg glucosamine / sucrose / PEG 8000 / ascorbic acid / ascorbic acid + PEG; SPI pre-incubated 20 min at 23 or 60 C before flavour addition; reference = flavour alone"
    medium: "5 g 0.05 M Tris-HCl pH 7.4 in 15 mL vial, stirred; HS-SPME 50/30 um DVB/CAR/PDMS 60 min at 23 or 60 C; GC-FID; relative % peak area vs flavour alone"
    temperature_c: [23, 60]
    ph: {initial: [7.4], measured_during: false, final: null}
    time_min: [20, 60]
    species: [3-methylfuran, 2-acetylfuran, limonene, delta-3-carene, 2-ethyl-3,6-dimethylpyrazine, 2,3-diethyl-5-methylpyrazine, 4-terpineol, decanal, 2-isoamyl-6-methylpyrazine, delta-elemene, beta-cubebene, calamenene, unknown_beef_correlated_peaks]
    kind: ratio
    location: "Table 1 (RT adsorption, 17 peaks) and Table 2 (60 C adsorption, 53 peaks), pp. 1231-1233; SPIF column = protein effect; dossier sec. 4a-4c, apparent K_g sec. 4f"
    status: digitised_in_dossier
    quant_class: semi_quantitative
    simulable: false
    condition_quote: "expressed as relative peak area (%) compared to the peak area of SBF"
    flags: [peak_area_only, hs_spme_competitive_not_equilibrium, reference_leg_contains_yeast_extract, single_protein_level, no_sulfur_compound_identified, tentative_rt_identification, zeros_are_non_detections, spif_vs_sbf_not_tested, ascorbic_arms_ph_unmeasured_likely_acidified, ratios_above_100_pct_unexplained, protein_denatured_as_supplied, "also_relevant: aroma_other"]
    note: "The corpus's only headspace series of a commercial meat (beef-type) reaction flavour against a plant protein. SPI at ~43 g/L leaves most beef-correlated peaks at ~40-60 % of flavour alone. Apparent single-site K_g (mine): 2-ethyl-3,6-dimethylpyrazine 0.013-0.020 L/g, decanal 0.018-0.020 L/g (6-7x below the dialysis-family extrapolation from kg_nonanal_soy; aldehyde method families are not pooled). Contact 23 vs 60 C: ratio flat for most peaks. No sulfur odorant identified (FID). Additive columns: ascorbic acid +- PEG raise several peaks, likely confounded by acidification (pH not measured). Simulable false: multicomponent flavour of unknown composition, no protein-free water leg."
  - id: moon2007a_tab3_spi_sh_ss_surface_hydrophobicity
    dossier: moon2007a_extraction.md
    paper: "Moon & Li-Chan 2007 Food Res Int 40:1227 (10.1016/j.foodres.2007.08.002)"
    lab: "UBC Vancouver (Li-Chan)"
    system: "commercial soy protein isolate (Solae, 89.3 % protein, no DSC transition 25-180 C) +- glucosamine, sucrose, PEG 8000, ascorbic acid, ascorbic acid + PEG"
    medium: "SH: Ellman DTNB in pH 8.0 Tris-glycine-EDTA, 10 mg/mL SDS; SH+SS: NTSB (Thannhauser); S0: PRODAN, 0.05 M Tris pH 7.4 after 60 C / 20 min"
    temperature_c: [60]
    ph: {initial: [7.4, 8.0], measured_during: false, final: null}
    time_min: [20, 40]
    species: [protein free thiol, protein disulfide, surface_hydrophobicity_S0]
    kind: endpoint
    location: "Table 3 (SH, SS n = 2 with letters; S0 mean +- SD n = 3), p. 1233; Raman Tables 4-5 (structure only) not part of this entry"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "should be considered with caution due to high values of the protein blanks"
    flags: [protein_site_density_candidate, commercial_isolate_denatured_as_supplied, free_sh_8_to_13x_below_soy_isolate_band, blank_turbidity_caveat_by_authors, ascorbic_arms_ph_unmeasured, "also_relevant: sulfur"]
    note: "Commercial SPI: free SH 0.6 and S-S 36.7 umol/g protein, i.e. at protein_matrices.yml soy_isolate heated_reference (Ruan 2014: 0.75 after 100 C / 30 min), far below the native free-SH band 4.7-8.0; S-S just below the 38-53 band. Evidence that a commercial isolate should be charged with heated-state thiol sites. Unused."
  - id: moon2007a_fig2_sensory_ascorbic_peg_treated_spi
    dossier: moon2007a_extraction.md
    paper: "Moon & Li-Chan 2007 Food Res Int 40:1227 (10.1016/j.foodres.2007.08.002)"
    lab: "UBC Vancouver (Li-Chan)"
    system: "150 mg simulated beef flavour + 250 mg SPI (untreated, ascorbic-acid treated, ascorbic acid + PEG treated) in 5 g water; 7 panellists x 2 replications; 5 attributes"
    medium: "15 mL vials, held 60 C / 20 min before serving; 15 cm unstructured line scale"
    temperature_c: [60]
    ph: {initial: null, measured_during: false, final: null}
    time_min: [20, 20]
    species: [roasted, beefy, yeasty, soymilk-like, cereal]
    kind: ranking
    location: "Fig. 2 cobweb (digitised +-0.3 cm, dossier sec. 4g); significance in text sec. 3.6"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: ""
    flags: [figure_only_digitised, sensory_medium_water_vs_gc_tris, sensory_dose_150mg_vs_gc_500mg, beefy_difference_not_significant, ph_not_stated]
    note: "Ascorbic acid + PEG pre-treatment of SPI raises roasted (6.8 -> 9.4) and lowers soymilk (6.9 -> 3.7) and cereal (6.4 -> 4.2); beefy 6.7 -> 8.8 n.s. SPIF arm agrees with moon2007b S3 on beefy, soymilk and cereal within 0.4 cm. Unused."
```

No entry for Tables 4-5 (Raman structure, no consumer) or the DSC (data not shown).

## 8. Supplementary material

**None referenced** (no Appendix A, no "supporting information"; checked by text search). The DSC thermograms are "data
not shown". The GC method details live in Moon & Li-Chan 2004 (Food Chem 88:141) and the compound identities in Moon,
Cliff & Li-Chan 2006 (Food Res Int 39:294), neither on disk; the 2006 GC-MS / GC-O paper is the one place a sulfur
odorant of this SBF could be named, so it is the paper to fetch if the sulfur question matters (P(it lists ≥1 sulfur
odorant in this SBF) ~65 %, a guess from the composition of yeast-extract reaction flavours). The companion's SI is
discussed in `moon2007b_extraction.md` sec. 8.
