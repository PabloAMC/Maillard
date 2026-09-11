# Bornhorst 2017b — EXTRACTION (mashed-potato model gel with 1-2 g/100 g D-ribose + 0.5-2 g/100 g L-lysine at pH 8.4 / 9.2 / 9.5, heated at 80, 90 and 100 C for up to 240 min; Arrhenius fits giving the corpus's only norfuraneol activation energies, 104.9-122.3 kJ/mol with z = 20.6-24.0 C, plus L* and a* barriers and a three-process pasteurization validation)

### PART II of the Washington State pair, and the source of the three activation energies `src/kinetic_core/parameters_sulfur.py` already carries as UNVALIDATED, NON-OPERATIVE alkaline priors: Ea(M-2) = 121.1 ± 8.1 / 122.3 ± 19.5 / 104.9 ± 8.9 kJ/mol — confirmed here from Table 1, together with the per-formula pH (8.4 / 9.2 / 9.5) that is the exact reason the sulfur lane strands them.

**Source on disk:** `data/articles/bornhorst2017b.pdf` (30+ pp. accepted manuscript, LWT — Food
Science and Technology; version of record
`https://www.sciencedirect.com/science/article/pii/S0023643817302384`, manuscript id
`3d0ea2393ab489cb7e1fd66acefdddbd`, Elsevier user licence, © 2017).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/bornhorst2017b.txt`). **This is an
accepted manuscript, double-spaced with line numbers, not the typeset article**: no volume, issue,
page range, DOI or publication date appears on its face, the figures are supplied separately with
captions, and **Tables 1-5 sit at the end after the reference list**. **Tables 1, 2, 3 and 4 came
through the text layer clean** and are re-typed in full below; their multi-row temperature blocks
extracted with the row labels offset, and each was checked cell by cell against the caption's
ordering and against the running text's stated ranges. **Table 5 is a grid of photographs and colour
maps and contains no number at all.** Figures 1 (M-2 vs time at 80 and 100 C), 2 (colour vs time),
3 (M-2, L* and a* against cook value) and 4 (pixel histograms) are images: **every concentration-time
datum and every colour-time datum in this paper is figure-only.** There is no supplementary
material.

**Repo status before this dossier.** `bornhorst2017b_extraction.md` is **quoted by name and
verbatim** in `data/lit/extraction_dossiers/k3_final_parameter_inventory.md` (§A.3.6(i), items B10.7
and B10.9), and this paper's three activation energies are carried operationally-inert in
`src/kinetic_core/parameters_sulfur.py` under `ALKALINE_PRIORS` with the source anchor "Bornhorst et
al. 2017b, LWT, Table 2 (1_R0.5_L / 1_R1_L / 2_R2_L)". **The dossier file itself was missing from
`data/lit/extraction_dossiers/`.** This document restores it from the primary source. Two
corrections to the registry entry fall out of that reading and are recorded in Flags 2 and 3.

## 0. Identity

| field | value |
|---|---|
| Title | "Thermal pasteurization process evaluation using mashed potato model food with Maillard reaction products" |
| Authors | Ellen R. Bornhorst, Juming Tang (corresponding, jtang@wsu.edu), Shyam S. Sablani, Gustavo V. Barbosa-Cánovas — Department of Biological Systems Engineering, Washington State University, Pullman, WA 99164-6120, USA. **Identical author list and address to Part I** |
| Venue | LWT — Food Science and Technology. Accepted manuscript, © 2017 published by Elsevier under the Elsevier user licence. **Article id S0023643817302384**; the manuscript carries no volume, page or DOI |
| Which of the pair | **PART II.** Part I is `bornhorst2017.pdf`, article id **S0023643816305825**, "Development of model food systems for thermal pasteurization applications based on Maillard reaction products", © 2016 |
| Keywords | quality; kinetics; chemical marker; microwave-assisted pasteurization |
| The marker | **M-2 = 4-hydroxy-5-methyl-3(2H)-furanone = norfuraneol**, registry id `norfuraneol`, engine species **`NF`** on the sulfur lane |
| Mechanism stated | "M-2 is formed as a result of **2,3 enolization of the Amadori compound** during the Maillard reaction between **D-ribose and an amine**" (from Kim et al. 1996) |
| The matrix | mashed potato model gel only (Part I's egg white and gellan are dropped, on Part I's own recommendation) |
| The process | **MAPS**, the Microwave Assisted Pasteurization System at Washington State University: patent-pending, **915 MHz** microwaves, single-mode cavities, **no overpressure**, four sections (preheating, microwave heating, holding, cooling), pilot scale, **8-20 oz (226.8-567.0 g) trays** |
| Naming | as Part I, plus D-value = time for a one-log (90 %) change, computed from k; z-value = the temperature change for a one-log change in D, computed from Ea |
| Companions on disk | `bornhorst2017_extraction.md` (Part I), `k3_final_parameter_inventory.md` |

## 1. Why it matters

**As with Part I, a correction to the premise first.** This is **not** a non-covalent binding study,
it measures no flavour compound and no aroma retention, and **it contains no whey protein**. Whey
appears only as a literature comparison — Lau et al. 2003's whey-protein-gel M-2 activation energies
at sterilization temperature, quoted as 64.0-122.3 kJ/mol. The matrix here is a **mashed potato and
gellan gel**. Nothing in this paper touches `src/kinetic_core/parameters_matrix.py` or
`src/kinetic_core/matrix_sites.py`: there is **no binding constant, no partition coefficient, no
protein loading in g/L and no aroma compound** anywhere in it.

**What it is, and exactly where it already sits in the repository.** It is the temperature-dependence
half of a two-part Maillard-marker kinetics study on **norfuraneol**, the engine's `NF` species and
the branch point `src/kinetic_core/sulfur.py` runs into 2-methyl-3-furanthiol (`r_nf_mft`) and
2-mercapto-3-pentanone (`r_nf_mp3p`). It supplies **the corpus's only norfuraneol activation
energies**. `src/kinetic_core/parameters_sulfur.py` carries them:

```
"quantity": "Ea, norfuraneol (M-2) net accumulation",
"value_set_kj_mol": [121.1, 122.3, 104.9],
"source_anchor": "Bornhorst et al. 2017b, LWT, Table 2 (1_R0.5_L / 1_R1_L / 2_R2_L)",
"rate_transfer": "not_licensed",   "operative": False,
```

with four mandatory qualifications travelling with them. **All four are confirmed from the primary
text**, and the confirmations are the main value of this dossier:

1. *"an APPARENT, LUMPED approach-to-plateau rate of NET accumulation with no destruction term at
   all"* — the paper says this itself, in words: **"The kinetic model utilized in this study only
   considered the formation of M-2 and did not take into account any elimination reactions occurring
   simultaneously for this intermediate compound."** (§2.3). That sentence is the licence for the
   module's forbidden-derivation entry "Bornhorst's norfuraneol k read as a DEGRADATION rate".
2. *"pH 8.4-9.5 against this module's pH 4.5-7"* — confirmed, and **this paper prints the pH per
   formula** where Part I gave only a range: **no precursors 5.2; 1_R,0.5_L 8.4; 1_R,1_L 9.2;
   2_R,2_L 9.5** (§3.1). So each Ea can now be tied to a specific pH: **121.1 at pH 8.4, 122.3 at
   pH 9.2, 104.9 at pH 9.5.**
3. *"three temperatures per fit, one of them imported from another paper"* — confirmed: the 80 and
   100 C rows are measured here, and **every 90 C row is copied from Part I**, footnoted on all three
   tables. So each Arrhenius line has **three points, one of which is a different manuscript's**.
4. *"CaCl2 in every row and the gellan matrix data withheld"* — confirmed: **0.13 g CaCl2·2H2O per
   100 g in every formula**, held constant deliberately because "Calcium chloride also affects the
   rate of Maillard browning (Kocadağli & Gӧkmen, 2016)". Gellan is absent from this paper entirely.

**What the paper adds beyond the three numbers.**

- **A z-value, and the transform that produced it.** z = 20.6-24.0 C for M-2, 20.8-28.8 C for L*,
  10.3-25.6 C for a*. I confirm below that every printed z reproduces from its Ea as
  **z = 2.303 R T²/Ea at T = 363.15 K (90 C)**, to within 0.05 C in all nine cases (mine). That fixes
  the reference temperature of the whole table, which the paper never states.
- **A dose-response on the barrier.** Ea(M-2) falls **121.1 -> 122.3 -> 104.9 kJ/mol** as precursor
  loading rises fourfold. `k3_final_parameter_inventory.md` item B10.7 already records that the
  norfuraneol Ea "is not a constant" and spans 64-122 kJ/mol across three studies as a function of
  precursor loading and matrix; **this paper is the within-study half of that claim**, and it shows
  the effect is smaller within one laboratory (a 1.17x spread) than across three (1.9x).
- **The plateau's temperature behaviour, which is the item the inventory quotes verbatim.** M-2∞ is
  essentially flat from 80 to 100 C. The inventory's sentence — *"the 1_R,0.5_L trend is
  0.57 ± 0.17 -> 0.43 ± 0.02, so the fall is inside one standard error at 80 C; the other two
  formulas are flat within noise. **This is suggestive corroboration, not a measurement of a sink.
  State it that way or not at all.**"* — is confirmed exactly against Table 1 below.
- **A process-validation block.** Three pasteurization processes matched to the same lethality
  (F90 = 10.9 / 11.0 / 10.9 min at the cold spot) but with different colour outcomes. That is a
  **quality-at-constant-safety comparison** and there is nothing like it elsewhere in the corpus.
- **A negative result the repository should carry.** Over 80-100 C the marker and the colour
  correlate **strongly with cook value (C100, z = 33 C) and only moderately with lethality (F90,
  z = 10 C)** — Pearson 0.89-0.96 against 0.39-0.58 for M-2. Part I, at one temperature, concluded
  the models were good for both. **Part II corrects Part I**, and the reason is printed: when several
  time-temperature histories are compared, "the model food z-value becomes important", and the
  markers' z-values (20-29 C) are far closer to the cook value's 33 C than to lethality's 10 C.

What this paper does NOT give the repository: any pH below 8.4 with precursors present; any water
activity (mentioned as a possible cause, never measured); any measured concentration of ribose or
lysine; any downstream product of M-2; any sulfur species; any aroma compound or binding constant;
any tabulated concentration-time datum; any density.

## 2. Methods as they matter to a model

- **The matrix, exactly as formulated (per 100 g).** 15 g instant mashed potato flakes (Oregon Potato
  Co., Boardman OR), 0.5 g low acyl gellan gum (Kelcogel F, CP Kelco), **0.13 g calcium chloride
  (CaCl2·2H2O, J.T. Baker)**, **1-2 g D-ribose**, **0.5-2 g L-lysine** (both Sigma-Aldrich), and
  80.37-84.37 g deionised distilled water. Identical to Part I's mashed-potato formula. Preparation:
  gellan and flakes mixed into water, heated to **90 C**, CaCl2 added, held **1 min at 90 C**, cooled
  to **60 C**, ribose and lysine mixed in, then cooled to 22 C to set. **The precursors never see
  more than 60 C before the kinetic run.**
- **Why CaCl2 is held constant, in the authors' words.** It sets the gel ("a strong, brittle, and
  heat stable gel"), **and** "Calcium chloride also affects the rate of Maillard browning (Kocadağli
  & Gӧkmen, 2016); for this reason, the amount of calcium chloride added to each formula was kept
  constant to maintain consistency." So calcium is a known accelerant present at a fixed level in
  every row and absent from every comparison.
- **The three formulas and their measured pH (22 C).** **1_R, 0.5_L -> pH 8.4**; **1_R, 1_L ->
  pH 9.2**; **2_R, 2_L -> pH 9.5**; **no added precursors -> pH 5.2**. This per-formula breakdown is
  new in Part II (Part I gave only "8.4-9.5" for this matrix) and it means **precursor loading and pH
  are still perfectly confounded** — they rise together by construction, and the paper names both as
  possible causes of the rate increase, plus a third: "lower water activity in sample formulas with
  more precursors could also have contribute[d]". **None of the three is separated and water activity
  is never measured** (Flags 1).
- **Thermal treatment (what "80 C" and "100 C" mean).** **1 mL cylindrical aluminium test cells**
  (Chung, Birla & Tang 2008), heated in a **water bath at 80 C** and an **ethylene glycol bath at
  100 C** (Haake DC 30), cooled in ice water at 0 C. **Come-up time 1.75 min**, defined as the time
  for the coldest spot to reach within **0.5 C** of set point, by calibrated type-T thermocouples.
  **80 C from 5 to 240 min and 100 C from 5 to 150 min, come-up time EXCLUDED**, in triplicate, with
  **the three replicates of each time point coming from three separate experimental batches** (so the
  replication is batch-level, not analytical). For colour, the slow 1_R, 0.5_L formula at 80 C was
  **extended to 360 min**. **The 90 C data are not measured here; they are taken from Part I.**
- **Method 1 — M-2 by HPLC.** As Part I: homogenised in **10 mmol/L H2SO4**, centrifuged, filtered;
  Agilent 1100 with diode-array detector, **100 x 7.8 mm fast acid analysis column** (Bio-Rad),
  **10 mmol/L H2SO4 at 1 mL/min**, **285 nm**, **25 uL injection**. **All three experimental
  replicates analysed twice (two analytical replicates).** External standard curves from commercial
  M-2. **Limit of detection 0.02 mg M-2 per LITRE**, from 3 x (standard deviation of the response) /
  (slope of the calibration curve), after Shrivastava & Gupta 2011 — note the LOD is per litre while
  every reported concentration is per gram of sample (Flags 6).
- **Method 2 — colour by computer vision.** Hardware from Pandit et al. 2007a (light pod, compact
  fluorescent bulbs, digital camera, acquisition software), settings and analysis from Part I:
  **15 frames per second, ISO 200, F11**, QPcard 203 reference card, quadratic RGB -> L*a*b*
  transform after Leon et al. 2006, MATLAB R2013a, a circle of **37 695 pixels** per sample. **b* was
  discarded** on the basis of Part I and preliminary data here (Pearson r with time below 0.6).
- **Method 3 — the two-step kinetic fit.** *Step one*, non-linear regression (Newton, SAS 9.2)
  fitting zero, first and second order; the first-order form is Equation 1,
  `C = C∞ − (C∞ − C0) exp(−k·t)`, multiplied by −1 for L*; **M-2_0 fixed at zero**; **all three
  replicates fed to the model rather than their means**, "to better reflect system variation".
  *Step two*, linear regression of the Arrhenius form, Equation 2, `k = A0 exp(−Ea/RT)`, with
  **R = 0.008314 kJ/K/mol**. D-value and z-value from Toledo 2007. **Explicitly: "The kinetic model
  utilized in this study only considered the formation of M-2 and did not take into account any
  elimination reactions occurring simultaneously for this intermediate compound."**
- **Method 4 — lethality and cook value.** `F90 = ∫10^((T−90)/z) dt` with **z = 10 C**;
  `C100 = ∫10^((T−100)/z) dt` with **z = 33 C**. Temperature measured at the **geometric centre of
  the test cell** (the cold spot) throughout, by calibrated type-T thermocouples. Correlation windows:
  the primary analysis used **F90 up to 100 min and C100 up to 60 min**, drawn from 0-240 min at 80 C,
  **0-90 min at 90 C (Part I's data)**, 0-5 min at 100 C for F90 and 0-45 min at 100 C for C100. A
  **secondary analysis restricted C100 to 20 min**, chosen because pilot-scale tests showed a maximum
  20 min cook value at the hot spot.
- **Method 5 — the validation.** **Duplicate, 280 g trays** of the **1_R, 1_L** formula (chosen for
  its correlations and its lower precursor cost), the two replicates from two separate batches.
  Rigid, high-barrier polypropylene trays (Printpack), **161 x 116 x 32 mm**, sealed with flexible
  lid-stock under **15 MPa of vacuum** (as printed — Flags 7). Cold spot located by initial tests and
  logged **every 2 s** with mobile metallic sensors (TMI-USA) after Luan et al. 2013; for the hot
  water processes the cold spot is the geometric centre of the tray. The three processes, all
  designed to 90 C for 10 min at the cold spot:

  | process | schedule | F90 at the cold spot | total time in 93 C water |
  |---|---|---:|---:|
  | **MAPS** | 30 min preheat in 61 C water; 2.3 min microwave heating with trays in 93 C water; 9 min hold in 93 C water; 5 min cooling in 23 C water | **10.9 min** | **11.3 min** |
  | **hot water, preheated** | 30 min preheat in 61 C water; 32.2 min in 93 C water; 10 min cooling in 5 C water | **11.0 min** | **32.2 min** |
  | **hot water, not preheated** | 38 min in 93 C water; 10 min cooling in 5 C water | **10.9 min** | **38 min** |

  After cooling to 22 C the trays were **sliced horizontally at ¼ and ½ the sample thickness**
  (quarter and middle layers) and the cut faces photographed. Colour maps used a **jet scale over
  L* 20-68 and a* 1-22**, ranges taken from the kinetic study's initial and saturation values.
  Histograms of the pixel values were normalised by the total pixel count.

## 3. Tables re-typed

### Table 1. "Predicted chemical marker M-2 concentration at saturation (M-2∞), reaction rate constant (k), decimal reduction time (D-value), activation energy (Ea), and thermal resistance constant (z-value) with estimated standard error (3 replicates) for mashed potato model food samples ... heated at 80, 90, and 100 °C. 90 °C reaction rates were repeated from Bornhorst et al. (2017) for comparison to other temperatures."

Footnote as printed: "*Data at 90 °C from **Bornhorst et al. (2016)**" — the caption says 2017 and the
footnote says 2016, for the same companion paper (Flags 8).

| Precursor amount | Temperature (°C) | M-2∞ (mg M-2/g sample) | k (10^-3 1/min) | D-value (min) | Ea (kJ/mol) | z-value (°C) |
|---|---:|---|---|---|---|---|
| 1_R, 0.5_L | 80 | 0.57 ± 0.17 | 1.5 ± 0.5 | 1588 ± 558 | | |
| 1_R, 0.5_L | 90 * | 0.47 ± 0.06 | 5.1 ± 0.9 | 449 ± 81 | **121.1 ± 8.1** | **20.8 ± 1.4** |
| 1_R, 0.5_L | 100 | 0.43 ± 0.02 | 13.2 ± 1.1 | 174 ± 14 | | |
| 1_R, 1_L | 80 | 0.18 ± 0.08 | 3.2 ± 1.8 | 720 ± 404 | | |
| 1_R, 1_L | 90 * | 0.24 ± 0.03 | 7.4 ± 1.8 | 310 ± 74 | **122.3 ± 19.5** | **20.6 ± 3.4** |
| 1_R, 1_L | 100 | 0.16 ± 0.01 | 30.0 ± 3.7 | 77 ± 9 | | |
| 2_R, 2_L | 80 | 0.15 ± 0.03 | 6.0 ± 1.8 | 387 ± 116 | | |
| 2_R, 2_L | 90 * | 0.10 ± 0.01 | 18.3 ± 4.1 | 126 ± 28 | **104.9 ± 8.9** | **24.0 ± 2.0** |
| 2_R, 2_L | 100 | 0.14 ± 0.01 | 40.3 ± 3.7 | 57 ± 5 | | |

### Table 2. "Predicted initial L* value (L*0), L* value at saturation (L*∞), reaction rate constant (k), decimal reduction time (D-value), activation energy (Ea), and thermal resistance constant (z-value) with estimated standard error (3 replicates) ... heated at 80, 90, and 100 °C."

| Model | Temp. (°C) | L*0 | L*∞ | k (10^-3 1/min) | D-value (min) | Ea (kJ/mol) | z-value (°C) |
|---|---:|---|---|---|---|---|---|
| 1_R, 0.5_L | 80 | 67.7 ± 0.6 | 46.0 ± 6.5 | 2.9 ± 1.3 | 805 ± 377 | | |
| 1_R, 0.5_L | 90 * | 68.2 ± 1.3 | 39.5 ± 13.1 | 6.0 ± 4.6 | 381 ± 274 | **87.7 ± 4.8** | **28.8 ± 1.6** |
| 1_R, 0.5_L | 100 | 69.1 ± 0.8 | 37.7 ± 2.3 | 14.2 ± 2.3 | 162 ± 26 | | |
| 1_R, 1_L | 80 | 63.6 ± 0.7 | 27.0 ± 6.9 | 4.6 ± 1.4 | 501 ± 151 | | |
| 1_R, 1_L | 90 * | 66.1 ± 1.0 | 35.4 ± 0.9 | 25.4 ± 2.9 | 91 ± 10 | **121.1 ± 36.6** | **20.8 ± 6.9** |
| 1_R, 1_L | 100 | 66.0 ± 1.6 | 33.1 ± 1.2 | 41.5 ± 6.3 | 55 ± 8 | | |
| 2_R, 2_L | 80 | 67.2 ± 1.5 | 25.6 ± 0.9 | 30.1 ± 2.8 | 76 ± 7 | | |
| 2_R, 2_L | 90 * | 62.5 ± 1.2 | 22.7 ± 0.5 | 69.3 ± 5.2 | 33 ± 2 | **107.2 ± 10.9** | **23.5 ± 2.4** |
| 2_R, 2_L | 100 | 64.5 ± 1.5 | 24.3 ± 1.1 | 213.8 ± 25.5 | 11 ± 1 | | |

### Table 3. "Predicted initial a* value (a*0), a* value at saturation (a*∞), reaction rate constant (k), decimal reduction time (D-value), activation energy (Ea), and thermal resistance constant (z-value) with estimated standard error (3 replicates) ... heated at 80, 90, and 100 °C."

| Model | Temp. (°C) | a*0 | a*∞ | k (10^-3 1/min) | D-value (min) | Ea (kJ/mol) | z-value (°C) |
|---|---:|---|---|---|---|---|---|
| 1_R, 0.5_L | 80 | 0.0 ± 0.4 | 37.0 ± 19.4 | 1.4 ± 0.9 | 1657 ± 1087 | | |
| 1_R, 0.5_L | 90 * | 1.7 ± 0.6 | 28.3 ± 5.3 | 6.5 ± 2.1 | 355 ± 115 | **158.1 ± 3.6** | **16.0 ± 0.4** |
| 1_R, 0.5_L | 100 | 0.0 ± 0.7 | 20.8 ± 0.8 | 24.9 ± 2.8 | 92 ± 10 | | |
| 1_R, 1_L | 80 | 1.6 ± 0.6 | 17.3 ± 0.9 | 12.9 ± 2.2 | 178 ± 30 | | |
| 1_R, 1_L | 90 * | 2.4 ± 0.6 | 20.7 ± 0.5 | 28.8 ± 3.0 | 80 ± 8 | **98.5 ± 7.7** | **25.6 ± 2.0** |
| 1_R, 1_L | 100 | 2.4 ± 0.9 | 20.9 ± 0.8 | 78.1 ± 12.4 | 29 ± 5 | | |
| 2_R, 2_L | 80 | 0.5 ± 0.9 | 27.9 ± 3.5 | 52.5 ± 13.4 | 44 ± 11 | | |
| 2_R, 2_L | 90 * | 2.7 ± 0.7 | 19.9 ± 1.0 | 169.9 ± 28.7 | 14 ± 2 | **245.6 ± 72.1** | **10.3 ± 3.3** |
| 2_R, 2_L | 100 | 3.0 ± 1.0 | 18.8 ± 0.6 | **4747 ± 1661** | 0.5 ± 0.2 | | |

### Table 4. "Pearson correlation coefficients for thermal lethality (F90) and cook value (C100) with experimental chemical marker M-2 concentration, L* value, and a* value (3 replicates) during heating at 80-100 °C. Chemical marker and color data were included for F90 up to 100 min and C100 up to 60 min."

Footnote: "* shows significant p-value of < 0.001."

| Thermal severity | Model | M-2 | L* | a* |
|---|---|---:|---:|---:|
| F90 | 1_R, 0.5_L | 0.58 * | −0.46 * | 0.58 * |
| F90 | 1_R, 1_L | 0.46 * | −0.54 * | 0.60 * |
| F90 | 2_R, 2_L | 0.39 * | −0.53 * | 0.40 (not significant) |
| C100 | 1_R, 0.5_L | 0.96 * | −0.61 * | 0.85 * |
| C100 | 1_R, 1_L | 0.89 * | −0.90 * | 0.85 * |
| C100 | 2_R, 2_L | 0.90 * | −0.71 * | 0.55 (not significant) |

### Table 5

"Example images of 1 g/100 g ribose, 1 g/100 g lysine mashed potato model food trays after a
conventional (hot water) pasteurization with and without preheating, Microwave Assisted
Pasteurization System (MAPS), and an unheated control. For each tray, the middle and quarter layers
are depicted by the original colored picture, L* value color map, and a* value color map."
**This "table" is a 3 x 6 grid of photographs and colour maps. It contains no number.**

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| pH at 22 C: no precursors / 1_R 0.5_L / 1_R 1_L / 2_R 2_L | **5.2 / 8.4 / 9.2 / 9.5** | §3.1 |
| M-2 first-order fit R², by formula | **> 0.99 / 0.97 / 0.94** | §3.1 |
| M-2 rate span across all formulas and temperatures | **1.5 to 40.3 x 10^-3 1/min** (D-values **57-1588 min**) | §3.1 |
| Arrhenius fit R² for M-2 | **above 0.98** for all three formulas | §3.1 |
| Ea and z for M-2, stated as ranges | **104.9-122.3 kJ/mol**, z **20.6-24.0 C** | §3.1 and Abstract |
| comparison: Lau et al. 2003, **whey protein** model food, 116-131 C | Ea **64.0-122.3 kJ/mol** | §3.1 |
| comparison: Pandit et al. 2006, mashed potato, formula **1.5 g/100 g ribose + 0 g/100 g lysine** | Ea **81.4-96.1 kJ/mol** | §3.1 |
| comparison: Kawai et al. 2004, glucose-lysine glassy matrices, 50-100 C | Ea **130-156 kJ/mol** for browning | §3.2 |
| L* first-order fit R², by formula | **above 0.84 / 0.92 / 0.95** | §3.2 |
| a* first-order fit R², by formula | **above 0.91 / 0.93 / 0.96** | §3.2 |
| L* rate span (D-values) | **2.9 to 213.8 x 10^-3 1/min** (**11-805 min**) | §3.2 |
| a* rate span (D-values) | **1.4 to 4747 x 10^-3 1/min** (**0.5-1657 min**) | §3.2 |
| Arrhenius fit R² for L* and a* | **above 0.92** | §3.2 |
| Ea and z ranges, L* / a* | **87.7-121.1 kJ/mol, z 20.8-28.8 C** and **98.5-245.6 kJ/mol, z 10.3-25.6 C** | §3.2 and Abstract |
| colour Ea range, whole study | **87.7-245.6 kJ/mol** | §3.2 |
| data restrictions applied for the fits | 1_R 1_L a* to **60 min at 100 C**; 2_R 2_L L* to **30 min at 100 C**; 2_R 2_L a* to **30 min at 80 C and 20 min at 100 C** | §3.2 |
| secondary correlation with C100 restricted to 20 min: M-2 | **0.79-0.96**, all three formulas significant | §3.3 |
| secondary correlation, C100 to 20 min: L* / a*, formulas 1_R 1_L and 2_R 2_L | **0.84-0.88** and **0.80-0.94** | §3.3 |
| secondary correlation, C100 to 20 min: L* and a*, formula 1_R 0.5_L | **0.23-0.42** (weak to moderate) | §3.3 |
| F90 at the cold spot: MAPS / hot water preheated / hot water not preheated | **10.9 / 11.0 / 10.9 min** | §2.4 |
| total time in 93 C water: MAPS / preheated / not preheated | **11.3 / 32.2 / 38 min** | §3.5 |
| median L* of the middle layer: control / MAPS / hot water preheated / hot water not preheated | **75.7 / 63.5 / 59.5 / 54.2** | §3.5 |
| median a* of the middle layer: control / MAPS / hot water preheated / hot water not preheated | **1.7 / 10.2 / 14.2 / 15.5** | §3.5 |
| interquartile range of the pixel values | **3-3.7** for L*, **1.4-2.4** for a*, similar across all treatments | §3.5 |
| colour-map ranges | L* **20-68**, a* **1-22** | §2.4 |
| M-2 limit of detection | **0.02 mg M-2/L** | §2.2 |
| safety target | **90 C for 10 min** at the cold spot, 6 log non-proteolytic *C. botulinum* (types B and E); the alternative is **70 C for 2 min**, 6 log *Listeria monocytogenes* | §1 |

**Everything time-resolved is figure-only**: Figure 1 (M-2 against time at 80 and 100 C), Figure 2
(colour against time), Figure 3 (M-2, L* and a* against cook value) and Figure 4 (normalised pixel
histograms for the middle layer). Table 5's images likewise carry no number.

### Arithmetic on the printed constants (all mine)

**1. Every z-value is computed from its Ea at 90 C.** The paper cites Toledo 2007 and does not state
the reference temperature. Using `z = 2.303 R T² / Ea` with R = 0.008314 kJ/K/mol and
**T = 363.15 K (90 C)**, so z = 2525.1/Ea:

| parameter | Ea printed | z from Ea at 90 C (mine) | z printed |
|---|---:|---:|---:|
| M-2, 1_R 0.5_L | 121.1 | 20.85 | 20.8 |
| M-2, 1_R 1_L | 122.3 | 20.65 | 20.6 |
| M-2, 2_R 2_L | 104.9 | 24.07 | 24.0 |
| L*, 1_R 0.5_L | 87.7 | 28.79 | 28.8 |
| L*, 1_R 1_L | 121.1 | 20.85 | 20.8 |
| L*, 2_R 2_L | 107.2 | 23.55 | 23.5 |
| a*, 1_R 0.5_L | 158.1 | 15.97 | 16.0 |
| a*, 1_R 1_L | 98.5 | 25.63 | 25.6 |
| a*, 2_R 2_L | 245.6 | 10.28 | 10.3 |

**All nine reproduce to within 0.05 C.** The reference temperature of every z in this paper is
therefore **90 C**, which is the midpoint of the 80-100 C window and the temperature Part I worked
at. A z-value read at any other reference is a different number.

**2. Every D-value is 2.303/k.** Checked on all nine M-2 rows and all eighteen colour rows: e.g.
2.303/0.0132 = 174.5 against a printed 174; 2.303/0.0300 = 76.8 against 77; 2.303/0.0403 = 57.1
against 57; 2.303/0.2138 = 10.8 against 11; 2.303/4.747 = 0.485 against 0.5. **This is a decimal
reduction time of the remaining distance to the plateau, computed from an approach-to-plateau
constant — it is NOT a 90 % loss time of the marker.** The same trap as in Part I.

**3. The Arrhenius fits reproduce from the two measured endpoints.** Using only the 80 and 100 C
rows (the two temperatures actually measured in this paper), Ea = R ln(k100/k80)/(1/353.15 −
1/373.15) = 54.79 x ln(k100/k80) kJ/mol:

| formula | k80 | k100 | ratio | Ea from two points (mine) | Ea printed (three points) |
|---|---:|---:|---:|---:|---:|
| 1_R, 0.5_L | 1.5 | 13.2 | 8.80 | **119.2** | 121.1 ± 8.1 |
| 1_R, 1_L | 3.2 | 30.0 | 9.38 | **122.6** | 122.3 ± 19.5 |
| 2_R, 2_L | 6.0 | 40.3 | 6.72 | **104.4** | 104.9 ± 8.9 |

**All three agree within 2 kJ/mol**, i.e. the imported 90 C point contributes almost nothing to the
barrier and the fit rests on the two temperatures measured here. That is reassuring about the import
and sobering about the leverage: **each Ea is essentially a two-point slope with a third point along
for the ride.**

**4. The barrier falls as loading rises, but only by 1.17x.** Ea(M-2) 121.1 -> 122.3 -> 104.9 kJ/mol
across a fourfold precursor increase. The first two are indistinguishable (their intervals overlap
completely: 113.0-129.2 against 102.8-141.8); the third is lower than the first by
**16.2 kJ/mol**, which is **just outside** the sum of their standard errors (8.1 + 8.9 = 17.0 —
so it is inside, marginally). **The correct statement is that this paper does NOT resolve a
loading-dependence of the barrier**, and item B10.7's cross-study 64-122 kJ/mol span is carried by
the differences *between* laboratories (Lau, Pandit, Bornhorst), not within this one.

**5. The plateau is flat with temperature, and the inventory's caution is exactly right.**
1_R, 0.5_L: 0.57 ± 0.17 -> 0.47 ± 0.06 -> 0.43 ± 0.02 mg/g at 80 / 90 / 100 C. The total fall is
0.14 mg/g against an 80 C standard error of **0.17** — **inside one standard error**. 1_R, 1_L:
0.18 -> 0.24 -> 0.16, **non-monotone**. 2_R, 2_L: 0.15 -> 0.10 -> 0.14, **non-monotone**. **There is
no temperature trend in M-2∞ in this paper.** The loading trend, by contrast, is clear at every
temperature: at 100 C, 0.43 -> 0.16 -> 0.14 mg/g (a **3.1x fall, mine**) while k rises 13.2 -> 30.0
-> 40.3 (a **3.05x rise, mine**).

**6. The a* barrier at the highest loading is not a measurement.** Ea(a*) for 2_R, 2_L is
**245.6 ± 72.1 kJ/mol**, a **29 % relative standard error**, fitted through a 100 C rate constant of
**4747 ± 1661 x 10^-3 1/min** — that is **k = 4.75 min^-1, a D-value of 0.5 min**, on a series whose
first sampling point is at 5 min and whose data were restricted to 20 min. **The reaction was over
before the first sample in every practical sense.** The same row's 80 C constant came from data
restricted to 30 min. Treat 245.6 kJ/mol as an artefact of an unresolved rate, not as a barrier; the
paper reports it without comment.

**7. The three processes were matched on safety and separated on quality — the useful result.**
F90 at the cold spot 10.9 / 11.0 / 10.9 min, i.e. **matched to within 1 %**. Median middle-layer L*
75.7 (unheated) -> 63.5 (MAPS) -> 59.5 (hot water, preheated) -> 54.2 (hot water, no preheat), so
**MAPS costs 12.2 L* units of the 21.5 that the harshest process costs (mine, 57 %)**. Median a*
1.7 -> 10.2 -> 14.2 -> 15.5: **MAPS delivers 8.5 of the 13.8 a* units (62 %, mine)**. The mechanism
is printed and it is time, not power: **11.3 min in 93 C water for MAPS against 32.2 and 38 min** —
a **2.8x and 3.4x** difference in hot-water contact (mine) at the same lethality.

**8. Cross-check against Part I.** Part I's Table 1 mashed-potato rows at 90 C are reproduced here
exactly for M-2 (0.47 ± 0.06 / 5.1 ± 0.9; 0.24 ± 0.03 / 7.4 ± 1.8; 0.10 ± 0.01 / 18.3 ± 4.1), for
every a* value, and for L* except one cell: **Part I prints the 1_R, 0.5_L L* rate as
6.0 ± 4.3 x 10^-3 1/min and Part II prints 6.0 ± 4.6** (Flags 9). Everything else matches to the
last digit.

**9. The markers are cook-value instruments, not lethality instruments, and the z-values say why.**
Mean Pearson r against C100 across the three formulas: M-2 **0.917**, |L*| **0.740**, a* **0.750**;
against F90: M-2 **0.477**, |L*| **0.510**, a* **0.527** (all mine, from Table 4). **M-2 correlates
with cook value 1.9x better than with lethality.** The stated reason checks out arithmetically: the
markers' z-values are 10.3-28.8 C, the cook value's z is 33 C and lethality's is 10 C, so **only the
2_R, 2_L a* row (z = 10.3 C) is anywhere near the lethality z — and it is the one row whose F90
correlation is not significant**, because its rate is unresolved (item 6). The one marker whose
z-value should have made it a lethality instrument is the one that was measured worst.

### Are these two manuscripts two parts of one study?

**Yes, and this one says so in its own introduction.** Part II cites "Bornhorst, Tang, Sablani, &
Barbosa-Cánovas, 2017" for all three model foods, states that "Bornhorst et al. (2017) compared the
performance and reaction kinetics at 90 °C for all three model foods ... They concluded that the
optimal model food for future research was mashed potato, **which is why this study utilized mashed
potato model food**", and names the gap: "Previous work by Bornhorst et al. (2017) ... determined
the chemical marker (M-2) and brown color formation kinetics at only one temperature (90 °C). The
temperature sensitivity ... is unknown." **Every 90 C row in Tables 1, 2 and 3 of this paper is
imported from Part I and footnoted as such**, so the two documents share data and neither Arrhenius
fit exists without the other. Part II also **revises** Part I's conclusion: Part I found the models
"equally good for both safety and quality evaluation at 90 °C"; Part II finds them "more relevant
for food quality evaluation than food safety", and explains the difference as a consequence of using
one temperature versus three. **Cite them together; do not treat Part II's Ea as independent of Part
I's 90 C k.** Cross-reference: `bornhorst2017_extraction.md`.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **M-2 IS keyed**: id `norfuraneol`, with
"4-hydroxy-5-methyl-3(2H)-furanone" already in its alias list; engine species **`NF`**, SULFUR lane.
**D-ribose, L-lysine, the Amadori compound and calcium chloride are all absent** — the registry
carries no Maillard reactants. There is no matrix entry for mashed potato anywhere in the
repository, and none is required: this paper supplies no quantity a matrix table consumes.

Every row below shares: **mashed potato + gellan model gel with 0.13 g/100 g CaCl2·2H2O, 1 mL
aluminium test cell, come-up time 1.75 min excluded, ice-water cooling, three replicates from three
separate batches, first-order approach-to-plateau fit with M-2_0 fixed at zero and NO elimination
term, Arrhenius reference temperature 90 C**. The 90 C rows are **Part I's data**.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **Ea, norfuraneol (M-2) net accumulation** | **121.1 ± 8.1 / 122.3 ± 19.5 / 104.9 ± 8.9** | kJ/mol | 80-100 C; **pH 8.4 / 9.2 / 9.5** respectively; 1_R0.5_L / 1_R1_L / 2_R2_L | **Table 1** | **measured_barrier** — the corpus's only norfuraneol Ea. Carried in `parameters_sulfur.ALKALINE_PRIORS`, `operative: False`, `rate_transfer: not_licensed` |
| **z-value, M-2** | **20.8 ± 1.4 / 20.6 ± 3.4 / 24.0 ± 2.0** | °C | as above, **reference 90 C** (mine, §3 item 1) | Table 1 | **measured_barrier** (a restatement of Ea, not independent information) |
| **k, M-2, 1_R 0.5_L at 80 / 90 / 100 C** | **1.5 ± 0.5 / 5.1 ± 0.9 / 13.2 ± 1.1** | 10^-3 min^-1 | pH 8.4 | Table 1 | **measured_rate** — approach-to-plateau FORMATION. **Reading as degradation is forbidden** (`parameters_sulfur.py`) |
| **k, M-2, 1_R 1_L at 80 / 90 / 100 C** | **3.2 ± 1.8 / 7.4 ± 1.8 / 30.0 ± 3.7** | 10^-3 min^-1 | pH 9.2 | Table 1 | measured_rate, same caveat |
| **k, M-2, 2_R 2_L at 80 / 90 / 100 C** | **6.0 ± 1.8 / 18.3 ± 4.1 / 40.3 ± 3.7** | 10^-3 min^-1 | pH 9.5 | Table 1 | measured_rate, same caveat |
| **M-2∞ (plateau), all three formulas at 80 / 90 / 100 C** | **0.57 ± 0.17, 0.47 ± 0.06, 0.43 ± 0.02**; **0.18 ± 0.08, 0.24 ± 0.03, 0.16 ± 0.01**; **0.15 ± 0.03, 0.10 ± 0.01, 0.14 ± 0.01** | mg M-2 per g of sample | as above | Table 1 | **level_only** — a fitted saturation. **No temperature trend survives the errors** (§3 item 5) |
| Ea and z, L* value | **87.7 ± 4.8 / 121.1 ± 36.6 / 107.2 ± 10.9** kJ/mol; z **28.8 ± 1.6 / 20.8 ± 6.9 / 23.5 ± 2.4** °C | kJ/mol, °C | 80-100 C, ref 90 C | Table 2 | **measured_barrier** (browning surrogate). The 1_R 1_L Ea carries a 30 % relative error |
| k, L*, nine rows | **2.9 / 6.0 / 14.2; 4.6 / 25.4 / 41.5; 30.1 / 69.3 / 213.8** (± in Table 2) | 10^-3 min^-1 | 80 / 90 / 100 C by formula | Table 2 | measured_rate (browning surrogate) |
| L*0 and L*∞, nine rows | L*0 **62.5-69.1**; L*∞ **22.7-46.0** | CIELAB L* | as above | Table 2 | level_only |
| Ea and z, a* value | **158.1 ± 3.6 / 98.5 ± 7.7 / 245.6 ± 72.1** kJ/mol; z **16.0 ± 0.4 / 25.6 ± 2.0 / 10.3 ± 3.3** °C | kJ/mol, °C | 80-100 C, ref 90 C | Table 3 | **measured_barrier**, except **2_R 2_L (245.6 ± 72.1) which should be refused** — its 100 C rate is unresolved (§3 item 6) |
| k, a*, nine rows | **1.4 / 6.5 / 24.9; 12.9 / 28.8 / 78.1; 52.5 / 169.9 / 4747** (± in Table 3) | 10^-3 min^-1 | 80 / 90 / 100 C by formula | Table 3 | measured_rate; the 4747 ± 1661 cell is **not resolved by the sampling schedule** |
| a*0 and a*∞, nine rows | a*0 **0.0-3.0**; a*∞ **17.3-37.0** | CIELAB a* | as above | Table 3 | level_only |
| **pH per formula** | **5.2 (no precursors) / 8.4 / 9.2 / 9.5** | pH at 22 C | mashed potato model | §3.1 | **level_only** — and the reason `sulfur.py` strands every rate above |
| D-values, M-2 / L* / a* | **57-1588 / 11-805 / 0.5-1657** | min | 80-100 C | Tables 1-3, §3.1-3.2 | **derived_assumption** — D = 2.303/k on an approach-to-plateau constant, **not a 90 % loss time of the marker** |
| Pearson r, M-2 vs C100 / vs F90 | **0.96 / 0.89 / 0.90** and **0.58 / 0.46 / 0.39** | — | 80-100 C, F90 to 100 min, C100 to 60 min | Table 4 | **within_study_ratio** (a statistic). All significant at p < 0.001 |
| Pearson r, L* vs C100 / vs F90 | **−0.61 / −0.90 / −0.71** and **−0.46 / −0.54 / −0.53** | — | as above | Table 4 | within_study_ratio |
| Pearson r, a* vs C100 / vs F90 | **0.85 / 0.85 / 0.55** and **0.58 / 0.60 / 0.40** | — | as above | Table 4 | within_study_ratio; the two 2_R 2_L cells are **not significant** |
| **cook-value versus lethality preference** | M-2 correlates **1.9x** better with C100 than with F90 (mean 0.917 against 0.477) | — | 80-100 C | (mine) from Table 4 | **within_study_ratio** — and it **overturns Part I's conclusion** |
| F90 at the cold spot: MAPS / hot water preheated / not preheated | **10.9 / 11.0 / 10.9** | min at 90 C, z = 10 C | 280 g trays, duplicate | §2.4 | **measured_rate**-adjacent process quantity; treat as level_only |
| time in 93 C water: MAPS / preheated / not preheated | **11.3 / 32.2 / 38** | min | as above | §3.5 | level_only |
| median middle-layer L*: control / MAPS / preheated / not preheated | **75.7 / 63.5 / 59.5 / 54.2** | CIELAB L* | 1_R 1_L formula, 280 g tray | §3.5 | **level_only** — a quality-at-constant-safety comparison, the paper's most transferable process result |
| median middle-layer a*: control / MAPS / preheated / not preheated | **1.7 / 10.2 / 14.2 / 15.5** | CIELAB a* | as above | §3.5 | level_only |
| interquartile range of pixel values | **3-3.7** (L*) and **1.4-2.4** (a*), similar across all treatments | CIELAB units | as above | §3.5 | level_only — the inherent colour variability of the matrix |
| M-2 limit of detection | **0.02** | mg M-2 per **litre** (unit mismatch with the mg/g results) | HPLC-DAD at 285 nm | §2.2 | **measured_bound** |
| Ea, M-2, **whey protein** model food at 116-131 C (Lau 2003) | **64.0-122.3** | kJ/mol | **NOT measured here** | §3.1, cited | **level_only, borrowed** |
| Ea, M-2, mashed potato at sterilization T, formula 1.5 g/100 g ribose + **0** lysine (Pandit 2006) | **81.4-96.1** | kJ/mol | **NOT measured here** | §3.1, cited | level_only, borrowed |
| Ea, Maillard browning, glucose-lysine glassy matrices 50-100 C (Kawai 2004) | **130-156** | kJ/mol | **NOT measured here** | §3.2, cited | level_only, borrowed |
| M-2 and colour time courses; pixel histograms | — | — | — | Figs. 1-4, Table 5 | **figure_only** |

### Can these be put on the same basis as the shipped constants?

**(a) The activation energies cannot be made operative, and the repository has already ruled.**
pH 8.4-9.5 against the sulfur module's pH 4.5-7. `sulfur.py` records the "alkaline block" as
stranded; `parameters_sulfur.ALKALINE_PRIORS` carries the three barriers with `operative: False`.
This dossier adds one thing the registry entry does not have: **which Ea goes with which pH**
(121.1 at 8.4, 122.3 at 9.2, 104.9 at 9.5), so a later wave that acquires a pH correction knows
where each number sits.

**(b) Two of the nine barriers should be refused outright even as priors.** Ea(a*) for 2_R, 2_L,
**245.6 ± 72.1 kJ/mol**, rests on a 100 C rate constant of 4.75 min^-1 (D = 0.5 min) fitted to data
whose first point is at 5 min — the reaction is complete before the first sample. Ea(L*) for
1_R, 1_L, **121.1 ± 36.6**, carries a 30 % relative error. Neither is a determination.

**(c) The Arrhenius fits are two-point slopes in practice.** Recomputing from only the 80 and 100 C
rows measured in this paper reproduces all three M-2 barriers to within 2 kJ/mol (§3 item 3). The
imported 90 C point adds precision, not independence. **Three temperatures, two of them
independent.**

**(d) z-values must carry their reference temperature.** Every z here is `2.303 R T²/Ea` at
**90 C** (verified on all nine, §3 item 1). A z quoted without that reference is not usable.

**(e) The cook-value finding is a structural claim and it transports.** Markers whose z-value is
20-29 C track a z = 33 C quality integral far better than a z = 10 C safety integral. That is a
statement about how time-temperature integrators work, not about mashed potato, and it applies to
any Maillard marker the engine might be asked to use as a process indicator.

**(f) Nothing here belongs to either matrix layer.** No binding constant, no partition coefficient,
no protein loading in g/L, no aroma compound, no site density. No entry should be created in
`parameters_matrix.py` or `matrix_sites.py` from this paper.

## 5. Flags

1. **Precursor loading, pH and water activity are all confounded, and the paper names all three
   without separating any.** "Model food systems with greater levels of added precursors also had
   higher pH ... In this pH range, higher pH could have contributed to the increased reaction rate
   ... Additionally, lower water activity in sample formulas with more precursors could also have
   contribute[d]." pH goes 5.2 -> 8.4 -> 9.2 -> 9.5 as loading rises; **water activity is never
   measured**. So the three Ea values are barriers of "the whole formula", not of a reaction at a
   defined pH and a_w.
2. **The registry's source anchor points at the wrong table.** `parameters_sulfur.ALKALINE_PRIORS`
   cites "Bornhorst et al. 2017b, LWT, **Table 2** (1_R0.5_L / 1_R1_L / 2_R2_L)" for the M-2
   activation energies. **The M-2 activation energies are in Table 1**; Table 2 is the L* table
   (whose barriers are 87.7 / 121.1 / 107.2, a different set that happens to contain a 121.1).
   The values carried, 121.1 / 122.3 / 104.9, are Table 1's and are correct; **only the table number
   is wrong**, and it is worth fixing because a reader following the anchor lands on a table that
   contains one of the same numbers for a different quantity.
3. **The registry labels the uncertainties `ci95_kj_mol`; the paper calls them standard errors.**
   Every Bornhorst table caption says "with estimated **standard error** (3 replicates)".
   `ALKALINE_PRIORS` stores 8.1 / 19.5 / 8.9 under the key `ci95_kj_mol`. **A standard error on three
   replicates is not a 95 % confidence interval** — a t-based 95 % interval on n = 3 is about
   **4.3x** the standard error. Anything that reads that field as a 95 % bound is under-stating the
   uncertainty by roughly a factor of four.
4. **Every k is an approach-to-plateau formation constant, and the paper states the omission
   explicitly**: "The kinetic model utilized in this study only considered the formation of M-2 and
   did not take into account any elimination reactions occurring simultaneously for this intermediate
   compound." The D-values compound the hazard: D = 2.303/k here is a decimal reduction time of the
   *remaining distance to the plateau*, not a 90 % loss time of the marker.
5. **One-third of every Arrhenius fit is imported from another manuscript.** All nine 90 C rows in
   Tables 1-3 are Part I's, footnoted. The fits are not independent of Part I and the two papers must
   be cited together.
6. **The limit of detection is in the wrong units for the results.** LOD "0.02 mg M-2/**L**" against
   concentrations reported in "mg M-2/**g sample**". Without the extraction volume-to-mass ratio the
   LOD cannot be compared to any reported value. Part I's method grinds 0.8 g into 8 mL, which if
   carried over would make the LOD **0.2 ug M-2 per g of sample (mine, 0.02 mg/L x 8 mL / 0.8 g)** —
   about 500x below the smallest reported plateau — but the volume ratio is not restated in this
   paper and that conversion is an inference.
7. **"sealed with flexible, plastic lid-stock under 15 MPa of vacuum".** 15 MPa is 150 atmospheres;
   a vacuum cannot exceed about 0.1 MPa. Presumably kPa or a differently defined setting. A
   typographic error, flagged so it is not propagated.
8. **The table footnotes and captions disagree on the companion paper's year.** All three captions
   say "repeated from Bornhorst et al. (**2017**)"; all three footnotes say "Data at 90 °C from
   Bornhorst et al. (**2016**)". Part I's copyright line reads © 2016 and the reference list of this
   paper cites it as 2017. The same document, two years.
9. **One cell disagrees between the two parts.** Part I's Table 1 gives the 1_R, 0.5_L mashed-potato
   L* rate constant at 90 C as **6.0 ± 4.3** x 10^-3 1/min; Part II's Table 2 gives **6.0 ± 4.6** for
   the same imported cell. Every other imported cell matches to the last digit. Minor, but it means
   the import was not a pure copy.
10. **Several fits use windows chosen after seeing the data.** 1_R, 1_L a* restricted to 60 min at
    100 C; 2_R, 2_L L* to 30 min at 100 C; 2_R, 2_L a* to 30 min at 80 C and 20 min at 100 C; and the
    1_R, 0.5_L colour series at 80 C extended to 360 min. Each restriction is disclosed and each is
    post-hoc, and they fall hardest on the fastest formula — which is the one whose barrier
    (245.6 ± 72.1 kJ/mol) is least believable.
11. **b* was excluded before this paper started**, on the strength of Part I plus "preliminary color
    data from this study" (Pearson r with time below 0.6). No b* number is reported anywhere.
12. **Section numbering skips.** The Results run 3.1, 3.2, 3.3, then **3.5**; there is no §3.4. No
    content appears to be missing, but a section number is.
13. **The validation is n = 2.** "Pasteurization validation studies were conducted in duplicate",
    two trays per process from two separate batches, and the colour comparison between MAPS and the
    two hot-water processes rests on those. The medians in §3.5 are medians of pixels within an
    image, not of replicates.
14. **This is an accepted manuscript, not the version of record.** No volume, pages, DOI or dates on
    the face of the document; figures separate; tables after the references. Any citation built from
    this file must fetch the typeset article for the bibliographic fields.
15. **What this paper does NOT contain**: any pH between 5.2 and 8.4 with precursors; any water
    activity measurement; any temperature below 80 C or above 100 C; any whey protein, egg white or
    gellan measurement (all three are Part I's); any aroma compound, binding constant or partition
    coefficient; any sulfur species or downstream product of M-2; any measurement of M-1 or M-3; any
    tabulated concentration-time datum; any density; any elimination term in the kinetic model; any
    supplementary material.
16. **What to request from the authors**: (i) the numeric M-2 and colour time courses behind
    Figures 1-3; (ii) the extraction volume-to-mass ratio so the LOD can be put on the reported
    basis; (iii) a pH-matched or a_w-matched control that separates the three confounded variables;
    (iv) water-activity measurements for the three formulas; (v) the covariance of each Arrhenius
    fit, given that one of three points is imported; (vi) whether the 2_R, 2_L a* series at 100 C
    was resolved at all, and if not a corrected Ea or its withdrawal; (vii) confirmation of the
    "15 MPa of vacuum" figure.
17. **Registry gaps against `data/keys/compounds.yml`**: **`norfuraneol` is present** and its alias
    list already covers this paper's name for M-2. **D-ribose, L-lysine and the Amadori compound are
    absent**, as is any matrix key for mashed potato — and none is needed, because this paper
    supplies no quantity a matrix table consumes. Two registry-side corrections belong to
    `src/kinetic_core/parameters_sulfur.py` rather than to `compounds.yml`: the source anchor's table
    number (Flags 2) and the `ci95_kj_mol` label (Flags 3).
