# Bornhorst 2017 — EXTRACTION (egg white, mashed potato and gellan model gels with 0-2 g/100 g D-ribose + 0-2 g/100 g L-lysine, heated at 90 C for 5-180 min in 18 x 4 mm aluminium test cells; first-order formation kinetics for the chemical marker M-2 / norfuraneol and for CIELAB L* and a*, with thermal-lethality and cook-value correlations)

### PART I of the Washington State pasteurization pair, and the corpus's only MATRIX-PAIRED norfuraneol measurement: at a single temperature (90 C) and a single formulation, egg white forms M-2 3.9x faster than mashed potato while reaching almost the same plateau — the matrix-transfer test `docs/reference/FIT_HOLDOUT_DECLARATION.md` already declares a STAR HOLD-OUT, written up here from the primary table for the first time.

**Source on disk:** `data/articles/bornhorst2017.pdf` (24 pp. accepted manuscript, LWT — Food Science
and Technology; version of record `https://www.sciencedirect.com/science/article/pii/S0023643816305825`,
manuscript id `11d43ba8dbee22c2c43a2631a0d5cf25`, Elsevier user licence, © 2016).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/bornhorst2017.txt`). **This is an
accepted manuscript, double-spaced with line numbers, not the typeset article**: it carries no
journal volume, issue, page range, DOI or publication date on its face, the figures are supplied
separately with captions at the end, and **Table 1 sits on the last page after the reference list**.
Table 1 came through the text layer clean and is re-typed in full below. Figures 1 (M-2 vs time),
2 (colour photographs of the model foods), 3 (M-2 against F90 and C100), 4 (L* against F90 and C100)
and 5 (a* against F90 and C100) are images: **every concentration-time datum, every colour datum and
every lethality datum in this paper is figure-only.** There is no supplementary material.

**Repo status before this dossier — and it is unusual.** `bornhorst2017_extraction.md` is **quoted by
name and verbatim** in `data/lit/extraction_dossiers/k3_final_parameter_inventory.md` (§A.3.6(i) and
the pH-wall discussion), its numbers are carried in
`src/kinetic_core/parameters_sulfur.py` under `ALKALINE_PRIORS` and in the module's
forbidden-derivation list, and `src/kinetic_core/sulfur.py` names "the Bornhorst norfuraneol Ea and
the whole alkaline block" as stranded. **The dossier file itself was missing from
`data/lit/extraction_dossiers/`.** This document restores it from the primary source; every number
the inventory attributes to this paper is checked against Table 1 below and they agree.

## 0. Identity

| field | value |
|---|---|
| Title | "Development of model food systems for thermal pasteurization applications based on Maillard reaction products" |
| Authors | Ellen R. Bornhorst, Juming Tang (corresponding, jtang@wsu.edu), Shyam S. Sablani, Gustavo V. Barbosa-Cánovas — Department of Biological Systems Engineering, Washington State University, Pullman, WA 99164-6120, USA |
| Venue | LWT — Food Science and Technology. Accepted manuscript, © 2016 published by Elsevier under the Elsevier user licence. **Article id S0023643816305825**; the manuscript carries no volume, page or DOI |
| Which of the pair | **PART I.** Its companion `bornhorst2017b.pdf` (article id **S0023643817302384**, "Thermal pasteurization process evaluation using mashed potato model food with Maillard reaction products", © 2017, same four authors, same address) is **PART II**. See "Are these two parts of one study?" below |
| Funding | USDA NIFA agreements 2011-68003-20096 and 2016-68003-24840; USDA National Needs Fellowship 2012-38420-19287 supporting E. Bornhorst's PhD |
| The marker | **M-2 = 4-hydroxy-5-methyl-3(2H)-furanone**, which is **norfuraneol**, registry id `norfuraneol` in `data/keys/compounds.yml` (the alias "4-hydroxy-5-methyl-3(2H)-furanone" is already listed there) and engine species **`NF`** on the sulfur lane |
| The other two markers, named but not measured | **M-1** = 2,3-dihydro-3,5-dihydroxy-6-methyl-4(H)-pyran-4-one (DDMP); **M-3** = 5-hydroxymethylfurfural, registry id `hmf`. Both are cited from Kim & Taub 1993 and Kim 1996 and neither is quantified here |
| Naming | M-2∞ = the fitted plateau (saturation) concentration; k = the fitted approach-to-plateau rate constant; F90 = accumulated thermal lethality in equivalent minutes at 90 C with z = 10 K; C100 = accumulated cook value in equivalent minutes at 100 C with z = 33 K; "0_R, 0_L" etc. = grams of D-ribose and L-lysine per 100 g of model food |
| Lineage | model foods from Zhang et al. 2014 (egg white), Zhang 2014/2015 (gellan) and Pandit et al. 2006 (mashed potato with xanthan, here re-made with gellan); regression method from Lau et al. 2003 (whey protein gels); test cell from Chung, Birla & Tang 2008; colour transform from Leon et al. 2006; F and C definitions from Toledo 2007 |
| Companions on disk | `bornhorst2017b_extraction.md` (Part II), `k3_final_parameter_inventory.md` (which already carries this paper's numbers as items B10.7, B10.8, B10.9 and B10.18) |

## 1. Why it matters

**First, a correction to the premise this dossier was commissioned under.** This paper is **not** a
non-covalent binding study, it measures no flavour compound, and **it contains no whey protein at
all**. Whey enters only through two citations — Lau et al. 2003, whose whey-protein-gel regression
method is borrowed, and Gupta et al. 2011, whey protein gels under pressure — and the paper's own
reason for abandoning whey is printed: the whey model is "not optimal for pasteurization due to
slower Maillard reaction kinetics and **high gelation temperatures (80 °C) for whey proteins**". The
three matrices actually measured are **egg white, mashed potato and gellan gum**. Nothing here
belongs to the matrix-retention layer (`src/kinetic_core/parameters_matrix.py`) or to the covalent
binding layer (`src/kinetic_core/matrix_sites.py`): **there is no binding constant, no partition
coefficient, no protein loading in g/L and no aroma compound in this document.**

**What it actually is, and where it already lives in the repository.** It is a **Maillard-product
formation-kinetics** paper on the sulfur lane's own species. M-2 is norfuraneol, engine species
`NF` in `src/kinetic_core/species_sulfur.py`, the branch point that
`src/kinetic_core/sulfur.py` runs into 2-methyl-3-furanthiol (`r_nf_mft`) and
2-mercapto-3-pentanone (`r_nf_mp3p`). So this paper's numbers are **rate constants and plateau
concentrations for a species the engine carries**, and the repository has already ruled on them:

- `src/kinetic_core/parameters_sulfur.py` carries the **Part II** activation energies
  (121.1 / 122.3 / 104.9 kJ/mol) in `ALKALINE_PRIORS` with `operative: False` and
  `rate_transfer: "not_licensed"`, because the pH here is 8.4-9.5 against the sulfur module's
  pH 4.5-7.
- The same module's forbidden-derivation list names, verbatim: **"Bornhorst's norfuraneol k read as
  a DEGRADATION rate"**, with the reason "THIS IS AN APPROACH-TO-PLATEAU FORMATION LAW, NOT A
  DEGRADATION LAW ... Anyone reading these k values as norfuraneol degradation rates would be
  inverting the paper." **This dossier confirms that reading from the primary text**: Equation 1 is
  `C = C∞ − (C∞ − C0) exp(−k·t)`, an approach-to-plateau, with M-2_0 fixed at zero. Every k in
  Table 1 is the rate at which M-2 *rises toward* its plateau.
- `docs/reference/FIT_HOLDOUT_DECLARATION.md` and the k3 inventory declare **"Bornhorst 2017 90 °C
  matrix pair (egg white vs mashed potato, M-2∞ and k)"** a **star HOLD-OUT** — "a matrix-transfer
  test at fixed T and fixed formulation — the model should predict a 3.9x k difference between two
  food gels it was not fitted to". **That 3.9x is confirmed here from Table 1: 19.9 / 5.1 = 3.90
  (mine)**, egg white over mashed potato at 1 g/100 g ribose + 0.5 g/100 g lysine.
- The **structural zero** — M-2 = 0 mg/g in all three matrices with no added precursors — is declared
  FIT ("a free, unambiguous zero; fitting it costs nothing and catches sign errors") and is
  confirmed here from §3.1.

**Why the paper is worth a dossier even though almost nothing in it is operative.** Three things:

1. **The matrix pair is a genuine controlled comparison.** Same marker, same analytical method, same
   temperature, same precursor loading, same test cell, same come-up time — and two food matrices
   that differ by an **18.6x lysine content (1.3 against 0.07 g/100 g, mine, from the USDA figures
   the paper quotes)**. Almost nothing else in the corpus isolates matrix that cleanly.
2. **The plateau moves the wrong way.** M-2∞ **falls** as precursor loading **rises** (0.54 -> 0.28
   -> 0.14 mg/g in egg white as ribose+lysine go from 1+0.5 to 2+2), while k rises. That is the
   signature of a marker that is an *intermediate* being consumed by its own downstream chemistry —
   which is exactly what `sulfur.py` models with `r_nf_decay` and the two thiol channels. **The
   paper says so** ("M-2 is an intermediate Maillard reaction product; as the reaction proceeds,
   furanones, such as M-2, may degrade into other smaller molecular weight color and flavor
   compounds") and **does not fit a destruction term**, which is precisely why its k must not be
   read as a degradation rate.
3. **A browning response measured on the same samples.** L* and a* first-order rates are printed for
   all nine treatments including gellan, so the paper offers a matched (marker, browning) pair.
   The engine's browning surrogate has no such matched anchor at pasteurization temperature.

What this paper does NOT give the repository: any activation energy (one temperature only — the Ea
is in Part II); any pH below 7.8 with precursors present; any measured concentration of ribose,
lysine or any sugar; any downstream product of M-2; any sulfur species; any aroma compound; any
tabulated concentration-time datum; any water activity; any whey.

## 2. Methods as they matter to a model

- **The three matrices, exactly as formulated (per 100 g of model food).**
  - **Egg white**: 25 g stabilised, **glucose-reduced** powdered egg white (JustWhites, Deb-El Food,
    Elizabeth NJ), 0-2 g D-ribose, 0-2 g L-lysine (both Sigma-Aldrich), remainder double-deionised
    water (71-75 g). Prepared by mixing the powder with **35 C** water for 10 min, heating at 35 C
    for 20 min to rehydrate, adding the precursors and mixing 30 min, filling the test cells, then
    **gelling at 70 C for 30 min** followed by ice water. **The gel is set by a thermal step that is
    itself a Maillard exposure** — 30 min at 70 C before the kinetic run begins (Flags 3).
  - **Mashed potato**: 15 g instant mashed potato flakes (Oregon Potato Co., Boardman OR), 0.5 g low
    acyl gellan gum, 0.13 g calcium chloride, 0-2 g D-ribose, 0-2 g L-lysine, remainder DDI water
    (80.37-84.37 g). Gellan mixed into 22 C water 5 min, flakes added, heated to **90 C**, CaCl2
    added, held at 90 C for **1 min**, cooled to **60 C** before the precursors are added and mixed
    5 min, then cooled to 22 C to set. **The precursors see 60 C at most before the run.**
  - **Gellan**: 1 g low acyl gellan gum (Kelcogel F, CP Kelco), 0.5 g titanium dioxide dispersed in
    glycerin and water (Wilton white-white icing colour, added to make the gel opaque), 0.26 g
    CaCl2·2H2O, 0-2 g D-ribose, 0-2 g L-lysine, remainder DDI water (printed as **84.24-98.24 g**,
    which does not close — Flags 5). TiO2 into 22 C water 5 min, gellan 5 min, heated to 90 C, CaCl2
    added, held 1 min at 90 C, cooled to **65 C** before the precursors, mixed 3 min, poured, cooled
    to 22 C.
- **The four formulas.** 0 g/100 g ribose + 0 g/100 g lysine (**0_R, 0_L**); 1 + 0.5 (**1_R, 0.5_L**,
  taken from Zhang et al. 2014); 1 + 1 (**1_R, 1_L**); 2 + 2 (**2_R, 2_L**). Twelve model foods in
  all; the 0_R, 0_L arm produces no marker and no colour and is not in Table 1.
- **pH, and it is the number that strands this paper.** Measured at 22 C. **Without precursors:
  egg white 6.0, mashed potato 5.2, gellan 6.1. With precursors: egg white 7.8-9.5, mashed potato
  8.4-9.5, gellan 9.8-9.9.** Adding free lysine raises the pH by up to 3.5 units, and the authors
  state this is part of the mechanism: "The higher pH of formulas with added precursors may have
  contributed to the increased rate of the Maillard reaction and M-2 formation ... This showed the
  importance of adding the precursors to adjust the pH". **So precursor loading and pH are
  confounded by construction and cannot be separated in this design** (Flags 1). The sulfur module's
  pH wall (`sulfur.py`, "pH 8.4-9.5 against this module's pH 4.5-7") rests on exactly these numbers.
- **Thermal treatment (what "90 C" means).** Cylindrical **aluminium test cells, 18 mm diameter x
  4 mm high** (Chung, Birla & Tang 2008), heated in an **ethylene glycol bath** (Haake DC 30). The
  **come-up time is 1.75 min**, defined as the time for the coldest spot to reach within **0.5 K** of
  target, measured with a calibrated type-T thermocouple. **All three model foods heated at 90 C for
  5 to 180 min, come-up time EXCLUDED**, then cooled in ice water at 0 C. **Triplicate.** So the
  reported times are isothermal holds with the ramp already subtracted — a cleaner thermal
  definition than most sources in the corpus.
- **Method 1 — M-2 by HPLC (what it measures).** Adapted from Zhang et al. 2014. **0.8 g of sample
  ground in 8 mL of 10 mmol/L H2SO4** extraction buffer, centrifuged, supernatant collected,
  filtered, sealed in a glass vial. Agilent 1100 HPLC with **diode-array detection**, a
  **100 x 7.8 mm fast acid analysis column** (Bio-Rad), **10 mmol/L H2SO4 mobile phase at
  1 mL/min**, **detection at 285 nm**. **25 uL injection, each sample analysed twice.** External
  standard curve from **commercial M-2** (Sigma-Aldrich). Result in **mg M-2 per g of sample**.
  **This is a UV-absorbance HPLC quantification against an authentic external standard** — no mass
  spectrometry, no internal standard, no isotope dilution.
- **Method 2 — colour by computer vision (what it measures).** CIELAB via a camera rig from Pandit
  et al. 2007a with new settings: white balance preset mode S on smooth white cardstock, aperture
  **F11**, **15 frames per second**, **ISO 200**, chosen against a **35-patch QPcard 203 reference
  card**. The reference card was photographed each session and used to correct colour and transform
  RGB to L*a*b* by the **quadratic model of Leon et al. 2006**. Analysis in **MATLAB R2013a** over a
  circle of **37 695 pixels** per sample. **This is a surface reflectance measurement of a gel disc,
  not an absorbance of a solution** — it cannot be compared to an A420 or A470 browning reading.
- **Method 3 — the kinetic fit.** Non-linear regression (Newton algorithm) in SAS 9.2, after Lau et
  al. 2003, fitting **zero, first and second order** rate equations. The first-order form is printed
  as Equation 1:

  > `C = C∞ − (C∞ − C0) exp(−k·t)`

  with C the parameter (M-2, L* or a*), C∞ its saturation value, C0 its initial value, k the rate
  constant in **1/min** and t the time in min. **For L* the equation is multiplied by −1** because
  L* falls. **M-2_0 was fixed at zero as a constant** because the measured initial M-2 was zero in
  every model. Order chosen by R². **Every k in Table 1 is therefore an approach-to-plateau
  constant of NET accumulation, with no destruction term in the model at all.**
- **Method 4 — lethality and cook value.** From Toledo 2007:
  `F90 = ∫ 10^((T−90)/z) dt` with **z = 10 K** for the target pathogen (non-proteolytic
  *C. botulinum*), and `C100 = ∫ 10^((T−100)/z) dt` with **z = 33 K** for overall food quality. The
  targets quoted from ECFF 2006 and FDA 2011 are **70 C for 2 min** (6 log *Listeria monocytogenes*)
  or **90 C for 10 min** (6 log non-proteolytic *C. botulinum* types B and E); *C. botulinum* at 90 C
  was chosen, which is why 90 C is the kinetic temperature.
- **Statistics.** SAS 9.2; Pearson correlation coefficients between the colour parameters and time,
  and between M-2 / L* / a* and F90 / C100; significance at **p < 0.05**. The lethality correlations
  were run **only over the first 60 min** at 90 C (an F90 range of 0 to about 60 min), because the
  safety minimum is F90 = 10 min and the hottest spot in the microwave system may reach about 50 min.

## 3. Tables re-typed

### Table 1 (last page). "Predicted M-2∞, L*0, L*∞, a*0, a*∞ and k with estimated standard error (3 replicates) for egg white, mashed potato, and gellan model food samples with added precursor amounts of 1 g/100g ribose and 0.5 g/100g lysine (1_R, 0.5_L), 1 g/100g ribose and 1 g/100g lysine (1_R, 1_L), and 2 g/100g ribose and 2 g/100g lysine (2_R, 2_L) during heating at 90 °C. The M-2 data were excluded for the gellan model due to low M-2 concentrations. M-2₀ concentration was assumed zero for all models."

Units exactly as printed: M-2∞ in **mg M-2 / g sample**; every k in **10^-3 1/min**; L* and a* are
dimensionless CIELAB coordinates.

| Model food | formula | M-2∞ (mg M-2/g sample) | k, M-2 (10^-3 1/min) | R² | L*0 | L*∞ | k, L* (10^-3 1/min) | R² | a*0 | a*∞ | k, a* (10^-3 1/min) | R² |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Egg white | 1_R, 0.5_L | 0.54 ± 0.01 | 19.9 ± 0.7 | 0.99 | 89.4 ± 0.8 | 59.0 ± 5.8 | 6.9 ± 2.2 | 0.91 | −2.0 ± 0.4 | 13.2 ± 0.6 | 17.5 ± 2.2 | 0.95 |
| Egg white | 1_R, 1_L | 0.28 ± 0.01 | 21.7 ± 1.8 | 0.99 | 85.3 ± 1.0 | 54.7 ± 2.6 | 12.5 ± 2.5 | 0.93 | −0.4 ± 0.7 | 14.6 ± 0.7 | 24.3 ± 3.9 | 0.89 |
| Egg white | 2_R, 2_L | 0.14 ± 0.01 | 22.0 ± 2.9 | 0.98 | 70.3 ± 0.9 | 38.7 ± 0.6 | 37.2 ± 3.3 | 0.97 | 6.6 ± 1.0 | 17.0 ± 0.4 | 67.8 ± 17.6 | 0.79 |
| Mashed potato | 1_R, 0.5_L | 0.47 ± 0.06 | 5.1 ± 0.9 | 0.99 | 68.2 ± 1.3 | 39.5 ± 13.1 | 6.0 ± 4.3 | 0.76 | 1.7 ± 0.6 | 28.3 ± 5.3 | 6.5 ± 2.1 | 0.93 |
| Mashed potato | 1_R, 1_L | 0.24 ± 0.03 | 7.4 ± 1.8 | 0.97 | 66.1 ± 1.0 | 35.4 ± 0.9 | 25.4 ± 2.9 | 0.95 | 2.4 ± 0.6 | 20.7 ± 0.5 | 28.8 ± 3.0 | 0.95 |
| Mashed potato | 2_R, 2_L | 0.10 ± 0.01 | 18.3 ± 4.1 | 0.93 | 62.5 ± 1.2 | 22.7 ± 0.5 | 69.3 ± 5.2 | 0.97 | 2.7 ± 0.7 | 19.9 ± 1.0 | 169.9 ± 28.7 | 0.97 |
| Gellan | 1_R, 0.5_L | — | — | — | 90.9 ± 2.1 | 66.7 ± 2.5 | 18.7 ± 6.3 | 0.68 | −0.1 ± 0.9 | 12.0 ± 1.3 | 17.0 ± 5.7 | 0.74 |
| Gellan | 1_R, 1_L | — | — | — | 91.5 ± 1.0 | 59.3 ± 1.0 | 23.2 ± 2.5 | 0.96 | 0.2 ± 0.6 | 12.8 ± 0.7 | 21.1 ± 4.0 | 0.90 |
| Gellan | 2_R, 2_L | — | — | — | 86.4 ± 1.5 | 43.3 ± 0.8 | 48.7 ± 4.8 | 0.96 | 0.2 ± 1.0 | 13.8 ± 0.4 | 99.3 ± 20.2 | 0.85 |

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| M-2 with no added precursors, all three matrices | **0 mg M-2/g sample** ("no significant M-2 formation") | §3.1 |
| pH at 22 C, 0_R 0_L: egg white / mashed potato / gellan | **6.0 / 5.2 / 6.1** | §3.1 |
| pH at 22 C with added precursors: egg white / mashed potato / gellan | **7.8-9.5 / 8.4-9.5 / 9.8-9.9** | §3.1 |
| gellan M-2 with precursors | **< 0.01 mg M-2/g sample**, "too low to measure accurately", excluded from the analysis | §3.1 |
| M-2 model fit, average R² across all treatments: zero / first / second order | **0.91 / 0.98 / 0.76** | §3.1 |
| egg white mean M-2 rate across the three formulas | **21.2 x 10^-3 1/min**, D-value average **108.8 min** | §3.1 |
| mashed potato M-2 rate range | **5.1 to 18.3 x 10^-3 1/min**, D-values **125.8-448.8 min** | §3.1 |
| egg white model composition at 25 g/100 g solids (USDA 2015) | ~**0** g/100 g total sugars, **1.3** lysine, **1.2** arginine, **0.8** methionine, **0.5** histidine | §3.1 |
| mashed potato model composition at 15 g/100 g solids (USDA 2015) | **0.5** g/100 g total sugars, **0.07** lysine, **0.06** arginine, **0.02** methionine and histidine | §3.1 |
| the four amino acids that lead to M-2 (from Pandit 2006) | lysine, arginine, methionine, histidine | §3.1 |
| a* model fit, average R²: zero / first / second | **0.71 / 0.89 / 0.69** | §3.2 |
| L* model fit, average R²: zero / first / second | **0.71 / 0.89 / 0.89** (six of nine treatments fit first and second order equally; two better to first; one better to second by 0.01 in R²) | §3.2 |
| L* vs time, Pearson r | all significant, **−0.72 to −0.94** | §3.2 |
| a* vs time, Pearson r | **0.71 to 0.95**, except mashed potato 2_R 2_L at **−0.17** (not significant), which improves to **0.88** when the data are cut to 20 min | §3.2 |
| b* vs time, Pearson r | **−0.87 to 0.90**; five significant positive, two significant negative, two not significant; four of nine between 0.02 and 0.57. **b* excluded from all kinetic analysis** | §3.2 |
| L* rate range across the three formulas: mashed potato / egg white / gellan | **63.3 / 30.3 / 30.0 x 10^-3 1/min** | §3.2 |
| a* rate range across the three formulas: mashed potato / gellan / egg white | **163.4 / 82.3 / 50.3 x 10^-3 1/min** | §3.2 |
| M-2 vs F90 and C100, Pearson r | all significant, **above 0.92**, average **0.97** over the six egg-white and mashed-potato treatments | §3.3 |
| L* vs F90 and C100, Pearson r | all significant, **above −0.73**, average **−0.87**, excluding mashed potato 1_R 0.5_L (worse than −0.7), which improves to **−0.93** when extended to 90 min | §3.3 |
| a* vs F90 and C100, Pearson r | all significant, **above 0.76**, average **0.86** | §3.3 |
| safety target | **F90 = 10 min** minimum at the cold spot for 6 log *C. botulinum*; the hot spot may reach **~50 min** | §3.3 |
| z-values used | **10 K** for lethality, **33 K** for cook value | §2.5 |
| come-up time | **1.75 min** to within 0.5 K of target | §2.2 |
| useful-temperature floor of each model | mashed potato **60 C and above**; egg white **70 C and above** (set by the temperature the precursors saw during preparation) | §3.4 |

**Everything time-resolved is figure-only.** Figure 1 (M-2 against time at 90 C, three replicates
with 95 % confidence intervals, egg white in A and mashed potato in B, with the fitted first-order
curves), Figure 2 (photographs of all twelve model foods during heating), and Figures 3, 4 and 5
(M-2, L* and a* against F90 and C100 over the first 60 min). Per house rule none is typed as a
number.

### Arithmetic on the printed constants (all mine)

**1. The matrix pair — the star hold-out, sized.** At fixed formulation and fixed temperature:

| formula | egg white k | mashed potato k | **k ratio (mine)** | egg white M-2∞ | potato M-2∞ | **M-2∞ ratio (mine)** |
|---|---:|---:|---:|---:|---:|---:|
| 1_R, 0.5_L | 19.9 | 5.1 | **3.90x** | 0.54 | 0.47 | **1.15x** |
| 1_R, 1_L | 21.7 | 7.4 | **2.93x** | 0.28 | 0.24 | **1.17x** |
| 2_R, 2_L | 22.0 | 18.3 | **1.20x** | 0.14 | 0.10 | **1.40x** |

**The matrix effect lives almost entirely in the RATE, not in the plateau** — up to 3.90x on k
against 1.15-1.40x on M-2∞ — **and it collapses as precursors are added**, from 3.90x at the lowest
loading to 1.20x at the highest. That collapse is the paper's own explanation working: egg white
carries **18.6x the lysine of mashed potato (1.3/0.07, mine)** and 1.2/0.06 = 20x the arginine, so at
low added precursor the matrix's native amino acids dominate and egg white runs fast; at 2 g/100 g
added lysine the added precursor swamps both matrices and they converge. **The hold-out's declared
3.9x is the LOW-loading cell only**, and a model that predicts 3.9x at every loading would be wrong
at the other two.

**2. The D-values check out, and they identify the transform.** The paper reports D-values without
defining them. D = ln(10)/k = 2.303/k reproduces every printed figure: egg white mean
2.303/0.0212 = **108.6 min** against a printed 108.8; mashed potato 2.303/0.0183 = **125.8 min**
against a printed 125.8 (exact) and 2.303/0.0051 = **451.6 min** against a printed 448.8 (the
difference is the rounding of k to 5.1). **So the "D-value" here is the first-order decimal
reduction time computed from the approach-to-plateau constant** — which is a decimal reduction time
of the *remaining distance to the plateau*, not of the marker. That is a second way the same trap
opens: a D-value invites reading as a destruction time, and it is not one.

**3. The plateau falls as the loading rises, and by almost the same factor in both matrices.** Egg
white M-2∞ 0.54 -> 0.28 -> 0.14 mg/g: a **3.86x fall (mine)** as ribose goes 1 -> 1 -> 2 and lysine
0.5 -> 1 -> 2. Mashed potato 0.47 -> 0.24 -> 0.10: a **4.70x fall (mine)**. Meanwhile k **rises**
1.11x (egg white) and **3.59x** (mashed potato). **More precursor gives faster formation of less
marker.** With no destruction term in the model, the only place that can come from is a
concentration-dependent sink on M-2 that the fit absorbs into C∞. `sulfur.py` carries exactly such
sinks (`r_nf_decay`, `r_nf_mft`, `r_nf_mp3p`, and the reductone-oxygen arm `ch_red_ox_nf`) — **so
this table is a qualitative check on the existence of an NF sink, and not a measurement of one.**

**4. The browning rate is the quantity that tracks precursor loading cleanly.** a* rate constants
rise monotonically with loading in all three matrices: egg white 17.5 -> 24.3 -> 67.8 (**3.87x**),
mashed potato 6.5 -> 28.8 -> 169.9 (**26.1x**), gellan 17.0 -> 21.1 -> 99.3 (**5.84x**), all mine.
L* rates likewise: 6.9 -> 12.5 -> 37.2 (**5.39x**), 6.0 -> 25.4 -> 69.3 (**11.6x**), 18.7 -> 23.2 ->
48.7 (**2.60x**). **The mashed-potato a* rate spans 26x across a 4x change in precursor, the steepest
dose response in the paper**, and it is why the authors chose that model food.

**5. Gellan browns without making the marker.** Gellan's L* and a* rates (18.7-48.7 and 17.0-99.3
x 10^-3 1/min) sit inside the range of the two food matrices, yet its M-2 is below 0.01 mg/g —
**more than 14x below the lowest measurable food value (0.14, mine)** and effectively zero. So in a
matrix with no protein and no native carbohydrate beyond the gellan itself, ribose + lysine still
brown at a normal rate and still produce essentially no norfuraneol. **Browning and norfuraneol are
decoupled**, which is a structural statement the engine can use and which no other corpus source
makes this cleanly.

**6. L*0 records how much browning happened before the run started.** Egg white L*0 falls
89.4 -> 85.3 -> 70.3 as precursors rise; mashed potato 68.2 -> 66.1 -> 62.5; gellan 90.9 -> 91.5 ->
86.4. **Egg white loses 19.1 L* units at t = 0 between its lowest and highest formula (mine)** —
after a 30 min gelation step at 70 C with the precursors already present. So the egg-white 2_R, 2_L
sample enters its "zero" already substantially browned (Flags 3).

**7. What is NOT derivable.** With one temperature there is **no activation energy in this paper**;
the Ea's are in Part II. There is also no way to separate the pH effect from the precursor effect
(they move together by construction), and no way to convert mg M-2 per gram of wet gel into mol/L
without a density, which is not printed.

### Are these two manuscripts two parts of one study?

**Yes — the same four authors, the same laboratory, the same marker, the same test cell, the same
analytical methods, and Part II's introduction picks up exactly where Part I stops.** Part I
(`bornhorst2017.pdf`, article id **S0023643816305825**, © 2016) *develops* three model foods and
characterises them at **one temperature, 90 C**, concluding that mashed potato is optimal. Part II
(`bornhorst2017b.pdf`, article id **S0023643817302384**, © 2017, "Thermal pasteurization process
evaluation using mashed potato model food with Maillard reaction products") takes **mashed potato
only** forward, adds the other temperatures to obtain **z-values (20.6-24.0 C for M-2, 20.8-28.8 C
for L*, 10.3-25.6 C for a*)** and the activation energies the repository carries, and validates
against a microwave-assisted pasteurization system and hot-water processes. Part II's own abstract
states the gap Part I left: "Model foods for pasteurization applications have been developed, but
studies on temperature sensitivity and validation are limited." **Read Table 1 here and Part II's
Table 2 together: this paper's 90 C mashed-potato column is one point of Part II's Arrhenius fit.**
Cross-reference: `bornhorst2017b_extraction.md`.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** **M-2 IS keyed**: id `norfuraneol`, display
"Norfuraneol (NF)", SMILES `CC1=C(O)C(=O)CO1`, InChIKey `DLVYTANECMRFGX-UHFFFAOYSA-N`, and the alias
list already contains "4-hydroxy-5-methyl-3(2H)-furanone" — the paper's exact name for M-2. The
engine species is **`NF`** (`src/kinetic_core/species_sulfur.py`), assigned to the SULFUR lane in
`engine.py`. **M-3 is also keyed** as `hmf`, though it is not measured here. **M-1 (DDMP) is not
keyed**, and neither is D-ribose, L-lysine, gellan gum or any matrix component — the registry
carries no Maillard reactants at all.

Every row below shares: **90 C, isothermal, come-up time of 1.75 min excluded, 5-180 min, 18 x 4 mm
aluminium test cell in an ethylene glycol bath, ice-water cooling, three replicates, first-order
approach-to-plateau fit with M-2_0 fixed at zero**. The pH is **7.8-9.9 depending on matrix and
formula** and is confounded with precursor loading.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **k, norfuraneol (M-2) net accumulation, egg white** | **19.9 ± 0.7 / 21.7 ± 1.8 / 22.0 ± 2.9** | 10^-3 min^-1 | 90 C, pH 7.8-9.5, 1_R0.5_L / 1_R1_L / 2_R2_L | Table 1 | **measured_rate** — an APPROACH-TO-PLATEAU FORMATION constant. **Reading it as a degradation rate is forbidden** (`parameters_sulfur.py`) |
| **k, norfuraneol (M-2) net accumulation, mashed potato** | **5.1 ± 0.9 / 7.4 ± 1.8 / 18.3 ± 4.1** | 10^-3 min^-1 | 90 C, pH 8.4-9.5, same three formulas | Table 1 | **measured_rate**, same caveat; this is the 90 C point of Part II's Arrhenius fit |
| **M-2∞ (plateau), egg white** | **0.54 ± 0.01 / 0.28 ± 0.01 / 0.14 ± 0.01** | mg M-2 per g of sample | as above | Table 1 | **level_only** — a fitted saturation, not a measured concentration |
| **M-2∞ (plateau), mashed potato** | **0.47 ± 0.06 / 0.24 ± 0.03 / 0.10 ± 0.01** | mg M-2 per g of sample | as above | Table 1 | **level_only** |
| M-2, gellan model | **< 0.01**, "too low to measure accurately" | mg M-2 per g of sample | 90 C, pH 9.8-9.9 | §3.1 | **measured_bound** — a one-sided ceiling, and a structural statement (browning without the marker) |
| **M-2 with no added precursors, all three matrices** | **0** | mg M-2 per g of sample | 90 C, 5-180 min, pH 5.2-6.1 | §3.1 | **structural_gate** — the free, unambiguous zero the FIT declaration already accepts |
| k, L* value change, egg white / mashed potato / gellan | **6.9 ± 2.2 / 12.5 ± 2.5 / 37.2 ± 3.3**; **6.0 ± 4.3 / 25.4 ± 2.9 / 69.3 ± 5.2**; **18.7 ± 6.3 / 23.2 ± 2.5 / 48.7 ± 4.8** | 10^-3 min^-1 | 90 C, computer-vision CIELAB on the gel surface | Table 1 | **measured_rate** (browning surrogate). Note L* fits first and second order equally well and was modelled first-order for comparability only |
| k, a* value change, egg white / mashed potato / gellan | **17.5 ± 2.2 / 24.3 ± 3.9 / 67.8 ± 17.6**; **6.5 ± 2.1 / 28.8 ± 3.0 / 169.9 ± 28.7**; **17.0 ± 5.7 / 21.1 ± 4.0 / 99.3 ± 20.2** | 10^-3 min^-1 | as above | Table 1 | **measured_rate** (browning surrogate) |
| L*0 and L*∞, all nine treatments | see Table 1 above (L*0 62.5-91.5; L*∞ 22.7-66.7) | CIELAB L* | as above | Table 1 | level_only — **L*0 is not a clean zero** (Flags 3) |
| a*0 and a*∞, all nine treatments | see Table 1 above (a*0 −2.0 to 6.6; a*∞ 12.0-28.3) | CIELAB a* | as above | Table 1 | level_only |
| **matrix-transfer ratio, k(egg white)/k(mashed potato)** | **3.90 / 2.93 / 1.20** at 1_R0.5_L / 1_R1_L / 2_R2_L | — | 90 C, fixed formulation | Table 1 (mine) | **within_study_ratio** — the ★ HOLD-OUT quantity; **the declared 3.9x is the low-loading cell only** |
| matrix-transfer ratio, M-2∞(egg white)/M-2∞(mashed potato) | **1.15 / 1.17 / 1.40** | — | as above | Table 1 (mine) | within_study_ratio |
| **plateau falls as loading rises** | **3.86x** (egg white) and **4.70x** (mashed potato) across a 4x precursor increase, while k RISES 1.11x and 3.59x | — | 90 C | Table 1 (mine) | **within_study_ratio** — qualitative evidence for an NF sink, **not a measurement of one** |
| D-value, M-2, egg white mean / mashed potato range | **108.8** and **125.8-448.8** | min | 90 C | §3.1 | **derived_assumption** — the paper's own arithmetic, D = 2.303/k on an approach-to-plateau constant; it is NOT a decimal reduction time of the marker |
| lysine content ratio, egg white / mashed potato model | **18.6x** (1.3 against 0.07 g per 100 g) | — | reconstituted models | §3.1 (mine, from the USDA figures the paper quotes) | **level_only, borrowed** — USDA 2015 composition tables, **not measured here** |
| total sugars, egg white / mashed potato model | **~0 / 0.5** | g per 100 g | reconstituted models | §3.1 | level_only, borrowed from USDA 2015 |
| pH with precursors, egg white / mashed potato / gellan | **7.8-9.5 / 8.4-9.5 / 9.8-9.9** | pH at 22 C | with 1-2 g/100 g ribose and 0.5-2 g/100 g lysine | §3.1 | **level_only** — and the reason the sulfur module strands every rate here |
| pH without precursors, egg white / mashed potato / gellan | **6.0 / 5.2 / 6.1** | pH at 22 C | 0_R, 0_L | §3.1 | level_only |
| correlation, M-2 vs F90 and C100 | **> 0.92**, mean **0.97** | Pearson r | first 60 min at 90 C, six treatments | §3.3 | **within_study_ratio** (a statistic, not a constant) |
| correlation, L* vs F90/C100 and a* vs F90/C100 | **mean −0.87** and **mean 0.86** | Pearson r | as above | §3.3 | within_study_ratio |
| activation energy for anything | — | kJ/mol | **one temperature only** | **absent — see Part II** | — |
| M-2 concentration-time data, colour-time data, lethality plots | — | — | — | Figs. 1-5 | **figure_only** |

### Can these be put on the same basis as the shipped constants?

**(a) The rate constants cannot be made operative, and the repository has already said so.** pH
7.8-9.9 against the sulfur module's declared pH 4.5-7 window; `sulfur.py` records the whole
"alkaline block" as stranded and `parameters_sulfur.ALKALINE_PRIORS` carries the Part II Ea's with
`rate_transfer: "not_licensed"` and `operative: False`. Nothing in Part I changes that: it is the
same matrix, the same precursors and the same pH range, at one temperature instead of three.

**(b) The matrix pair CAN be scored as a hold-out, and it should be scored per formula.** It is a
ratio, so the alkaline absolute scale cancels in the same way `parameters_matrix.py` lets Meynier's
partition ratios be used while refusing his absolute K. **But the ratio is not a single number**: it
is 3.90 / 2.93 / 1.20 across the three loadings, and a model scored only against 3.90 would be
scored against the easiest cell.

**(c) The structural zero costs nothing and should be fitted.** M-2 = 0 with no added precursors in
three different matrices at 90 C for up to 180 min. It is unambiguous, it needs no unit conversion,
and it catches sign errors.

**(d) The plateau values are fitted saturations, not concentrations.** M-2∞ comes out of a
non-linear regression whose model has no sink term. It is the asymptote the data implied over
180 min, and the 2_R, 2_L runs reached it fastest. Treat it as `level_only`, never as an equilibrium.

**(e) Nothing here belongs to either matrix layer.** No binding constant, no partition coefficient,
no protein loading in g/L, no aroma compound, no site density. `parameters_matrix.py` and
`matrix_sites.py` gain nothing from this paper, and no entry should be created in either from it.

## 5. Flags

1. **Precursor loading and pH are confounded by construction and cannot be separated.** Adding
   0.5-2 g of free L-lysine per 100 g raises the pH from 5.2-6.1 to 7.8-9.9 — up to 3.5 units — and
   the authors treat that as part of the mechanism, not as a nuisance: "This showed the importance of
   adding the precursors to adjust the pH and facilitate the Maillard reaction". **Every rate,
   plateau and colour trend attributed to "more precursor" in this paper is equally attributable to
   "higher pH".** No formula holds pH constant while moving the loading, and no pH is reported at
   90 C (all are at 22 C, and pH falls with temperature and with Maillard progress).
2. **Every k is an approach-to-plateau formation constant with NO destruction term in the model.**
   Equation 1 is `C = C∞ − (C∞ − C0)exp(−kt)` with M-2_0 fixed at zero. **Reading these as
   norfuraneol degradation rates inverts the paper**, and `parameters_sulfur.py` already carries that
   as a named forbidden derivation. The D-values compound the hazard: D = 2.303/k here is a decimal
   reduction time of the *remaining distance to the plateau*, not of the marker.
3. **The egg-white "zero" is not zero.** The egg-white gel is set by **30 min at 70 C with the
   precursors already mixed in**, before the kinetic clock starts. L*0 falls from 89.4 to 70.3 across
   the three formulas, so the 2_R, 2_L sample has already lost **19.1 L* units (mine)** at t = 0.
   Mashed potato and gellan see only 60-65 C during preparation and their L*0 spreads are 5.7 and
   5.1 units. **The egg-white kinetics therefore start from a different chemical state than the other
   two matrices, and the matrix-pair ratio inherits that.** The paper acknowledges the temperature
   floors (egg white "useful at 70 C and above", mashed potato "at 60 C and above") but does not
   connect them to the zero.
4. **M-2 is quantified by UV absorbance at 285 nm with an external standard, in a browning gel.** No
   internal standard, no mass spectrometry, no isotope dilution, and the sample matrix darkens by
   design over the run. The extraction is a 10 mmol/L H2SO4 grind — acid, which is a condition under
   which furanones interconvert. Systematic level error is plausible; the *shape* the kinetics is
   fitted to is more robust than the absolute mg/g.
5. **The gellan water range does not close.** The formula lists 1 g gellan + 0.5 g TiO2 + 0.26 g
   CaCl2 = 1.76 g plus 0-4 g of precursors, so the water should run **94.24-98.24 g**; the paper
   prints **84.24-98.24 g**. The egg-white (71-75) and mashed-potato (80.37-84.37) ranges both close
   exactly, so this is a typographic error in one digit, not a different formulation. Recorded so it
   is not propagated.
6. **The L* kinetics are order-ambiguous and were forced to first order for convenience.** Average
   R² is **0.89 for both first and second order**; six of nine treatments fit the two equally well.
   The paper says plainly why it chose first: "For ease of comparison with a* reaction rates in this
   study and previous literature ... L* value change was modeled with first order kinetics." **The L*
   rate constants are therefore a modelling choice, not a determination**, and their numerical values
   would change under a second-order fit.
7. **b* was dropped after the fact.** Nine treatments, correlations with time from −0.87 to +0.90,
   five significantly positive, two significantly negative, two not significant. b* was excluded from
   all kinetic analysis. That is a defensible decision, but it means the colour response is
   characterised on two of three CIELAB axes and the third disagrees with itself across matrices.
8. **Two correlation windows were adjusted after inspection.** The mashed potato 2_R, 2_L a*-vs-time
   correlation was **−0.17 over 180 min and 0.88 over 20 min**, and the 20 min window was adopted
   "for the remainder of the analysis". The mashed potato 1_R, 0.5_L L*-vs-lethality correlation was
   worse than −0.7 over 60 min and **−0.93 over 90 min**, and the paper reports both. Both changes
   are disclosed; both are post-hoc.
9. **Composition is borrowed, not measured.** The lysine, arginine, methionine, histidine and sugar
   contents of the egg white and mashed potato models come from **USDA 2015 food-composition
   tables**, not from an assay of the actual ingredients. The 18.6x lysine ratio that carries the
   paper's whole matrix explanation is therefore a generic-composition number. The egg white is also
   **glucose-reduced** by the supplier, which is why its "total sugars" is quoted as ~0 — so in the
   egg-white model **all the sugar is added ribose**, while mashed potato brings 0.5 g/100 g of its
   own.
10. **This is an accepted manuscript, not the version of record.** No volume, no pages, no DOI, no
    dates on the face of the document; the figures are separate; Table 1 follows the references. Any
    citation built from this file must fetch the typeset article for the bibliographic fields.
11. **What this paper does NOT contain**: any activation energy or z-value (one temperature — they
    are in Part II); any measurement at any temperature other than 90 C; any pH between 6.1 and 7.8;
    any whey protein; any aroma compound, binding constant or partition coefficient; any sulfur
    species or any downstream product of M-2; any measurement of M-1 or M-3; any tabulated
    concentration-time datum; any water activity; any density (so mg/g cannot become mol/L); any
    supplementary material.
12. **What to request from the authors**: (i) the numeric M-2 and colour time courses behind
    Figures 1-5; (ii) a pH-matched control that separates pH from precursor loading; (iii) an assay
    of the actual ingredients rather than USDA figures; (iv) the L* second-order fits, since they fit
    as well as the first-order ones; (v) whether any M-2 destruction was measured or attempted;
    (vi) the density of each gel, which would let mg/g become mol/L.
13. **Registry gaps against `data/keys/compounds.yml`**: **`norfuraneol` is present and its alias
    list already covers this paper's name for M-2**; `hmf` is present (M-3, unmeasured here). **M-1
    (2,3-dihydro-3,5-dihydroxy-6-methyl-4(H)-pyran-4-one, DDMP) is absent**, and so are D-ribose and
    L-lysine — the registry carries no Maillard reactants. There is no matrix entry for egg white,
    mashed potato or gellan anywhere in the repository, and none is needed: this paper supplies no
    quantity that a matrix table consumes.
