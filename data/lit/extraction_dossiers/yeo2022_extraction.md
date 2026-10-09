# Yeo 2022 (PhD thesis) — EXTRACTION (chicken triglyceride vs phospholipid in a defatted chicken matrix at 100 °C, a 30-150 min kinetics series, sunflower and soy lipids as substitutes, boiled beef and chicken GC-O)

**Source on disk:** `data/articles/Yeo2022.pdf` (the Reading thesis, 174 pp.; PDF page = printed page + 13).
Table of contents read first (`pdftotext -l 12`); then **Chapter 2** (pp. 38-65, beef and chicken broth
GC-O), **Chapter 3** (pp. 66-120, TG vs PL, LN-SDE and the kinetics experiment), **Chapter 4** (pp. 121-148,
sunflower and soy lipids) and the first pages of **Chapter 5** (pp. 153-156) read by eye on 2026-10-09.
Every table value below was read on the page image (Read tool) and cross-checked against
`pdftotext -layout`; they agree (pdftotext drops the µ in two units, restored from the image). Chapter 1
(literature review) and the appendices (RRFs, TLC) were not read. Written for the lipid lane and for
plant-lipid formulation, not for the core fit.

Published papers from this thesis on file: **none**. No chapter carries a publication note. `papers.yml`
holds two rows cited as "Yeo & Mottram (2023)", neither with a dossier:
`doi_10_1016_j_foodchem_2023_136009` (feeding `process_state_calibrations.json`
`yeo_mottram_2023_soy_lecithin_thiophene_uplift`: pea/soy protein isolate cysteine-ribose model, soy
lecithin 0.5-1 % w/w, alkylthiophene uplift 2.4-fold, sensory confirmation) and
`doi_10_48683_1926_00111736` (`deep_research_backlog.json`). Mottram is not a supervisor of this thesis and
**nothing in it matches that calibration record**: no protein isolate, no cysteine-ribose model, no lecithin
dose series, no sensory test of lipid samples beyond a free-choice profile, and no alkylthiophene detected
in any chapter (Ch. 4 p. 132: "Only 3-thiazolines and pyrazines were found"). If that record was built from
this work it is misattributed; its DOI and authorship need checking before it is used. The PDF prints no DOI.

| field | value |
|---|---|
| Title | "Elucidation of the role of different sources of phospholipids in meat aroma formation" |
| Author | HuiQi Yeo (as printed on the title page and in the PDF metadata) |
| Degree | PhD thesis, Department of Food and Nutritional Sciences, University of Reading, October 2022 (declaration dated 6 October 2022) |
| Supervisors | Jane K. Parker, Dimitrios P. Balagiannis, Jean H. Koek (Unilever) |
| Funding | Unilever and the Graduate School (p. iv) |
| DOI | none printed |

## 1. Methods

**Chapter 2, broths (pp. 41-43).** Retail beef silverside and chicken breast (Ross 308), trimmed, minced
(4.5 mm). 500 g meat + 500 g water boiled at 100 °C for 30 min in the flask, then Likens-Nickerson SDE for 2 h
with 30 mL pentane:diethyl ether 9:1, concentrated to 0.5 mL. **GC-O**: 3 trained assessors (Sniffin'
Sticks score ≥ 38), duplicate runs on Rxi-5Sil MS and ZB-Wax, intensity 1-10 (3 weak, 5 medium, 7 strong);
the printed "Intensity" is the **sum over duplicates and assessors, maximum 60**, "Detection frequency"
maximum 6. Identification by GC-MS, LRI and odour against authentic standards; two 3-thiazolines
synthesised to confirm identity. No quantification in Ch. 2.

**Chapter 3, the matrix (pp. 70-73, 81).** Chicken breast minced and freeze-dried, then two-stage Soxhlet:
petroleum ether (→ chicken triglyceride TG^C, Florisil-purified) then chloroform:methanol 2:1 (→ chicken
phospholipid PL^C, Folch-washed), 50-60 °C. The defatted meat lost on average **80 % of its sugars and sugar
phosphates and 40 % of its amino acids** to the extraction (p. 81); a preliminary SPME and free-choice
profiling test found adding them back unnecessary (pp. 81-83). Fatty acids by FAME GC-FID (Table 3-4).

**Chapter 3, LN-SDE (pp. 75-76; Table 3-3, p. 102).** 25 g sample = 5.8 g defatted meat + 0.5 g lipid
(TG^C, PL^C, 1:1, or MCT in the "DF" control) + 18.7 g water, i.e. **2 % lipid w/w** (derived here,
0.5/25); boiled 100 °C for 30 min, LN-SDE 2 h; GC-MS (HP-5MS). Internal standards 2-isopropylpyrazine (S/N
compounds) and 4-nonanone (rest); **"approx. ng" = semi-quantified with relative response factors from
authentic standards** (Appendix D, not read), n = 4, one-way ANOVA with Fisher LSD.

**Chapter 3, kinetics (pp. 73-74, 77; Table 3-1, p. 102).** 1 g sample in a 20 mL vial = 0.235 g defatted
meat + **0.020 g lipid** (2 % w/w) + 0.745 g water; CTG (TG^C), CTGPL (1:1) and CPL (PL^C). Vials stirred,
**left open 1 min to admit oxygen**, sealed, heated in a **100 °C water bath for 30, 60, 90, 120 and
150 min**. Immediately after heating, 5 µL methanol carrying 100 ng/µL 4-nonanone (printed "4-nonane",
p. 77) and 500 ng/µL (E,E)-2,4-decadienal-d4 injected through the septum, then HS-SPME (DVB/CAR/PDMS,
60 °C, 5 min incubation + 20 min extraction), GC-MS (ZB-5MSi) with SIM at m/z 81, 152, 156.
**(E,E)-2,4-Decadienal is quantified by stable-isotope dilution ("ng")**; every other compound is "approx. ng"
against 4-nonanone with RRFs (Appendix E, not read). Because the d4 standard goes in after heating, it
measures what is left, not formation and loss separately (Ch. 5 p. 155 says so). n = 4. The tables do not
say "per vial"; the amounts are per 1 g sample by the design (inferred here). **No rate constant is fitted
anywhere in the thesis.**

**Chapter 4 (pp. 124-128; Table 4-1, p. 144).** Same LN-SDE design and 2 % lipid, nine samples: TG, TG:PL
1:1 and PL from chicken (extracted as above), sunflower (TG^SF, KTC; PL^SF, Thew Arnott LEC5636) and soy
(TG^SY, Haepyo; PL^SY, Thew Arnott LeciTAs 5348). Two-way ANOVA (class × source). **The chicken columns of
Table 4-3 are the Chapter 3 LN-SDE data re-used** (identical values), not a new run.

## 2. Findings that matter

### 2a. The kinetics series, 100 °C, chicken lipids (Table 3-9, pp. 109-110)

| compound (ng or approx. ng) | sample | 30 min | 60 | 90 | 120 | 150 |
|---|---|---|---|---|---|---|
| (E,E)-2,4-decadienal (SIDA, ng) | CTG | 14.4 | 19.4 | 25.8 | 32.0 | 30.9 |
| | CTGPL | 239 | 301 | 149 | 107 | 96.7 |
| | CPL | 567 | 289 | 142 | 104 | 73.1 |
| hexanal | CTG | 147 | 184 | 214 | 367 | 439 |
| | CTGPL | 462 | 1280 | 1260 | 1370 | 1350 |
| | CPL | 606 | 1000 | 1250 | 1510 | 1200 |
| 2-pentylfuran | CTG | 23.5 | 58.0 | 92.9 | 229 | 288 |
| | CTGPL | 195 | 1180 | 1670 | 1920 | 2580 |
| | CPL | 467 | 1220 | 2390 | 3380 | 3630 |
| (E)-2-octenal | CTG | 6.71 | 8.22 | 12.0 | 19.9 | 28.1 |
| | CTGPL | 23.4 | 71.8 | 69.5 | 62.5 | 75.8 |
| | CPL | 38.8 | 54.3 | 48.4 | 45.0 | 39.8 |
| (E,Z)-2,4-nonadienal | CTG | n.d. | n.d. | 2.45 | 5.00 | 7.96 |
| | CTGPL | 3.68 | 27.7 | 41.7 | 45.8 | 63.8 |
| | CPL | 10.5 | 35.8 | 83.5 | 115 | 135 |
| octanal | CTG | 33.7 | 80.5 | 120 | 247 | 288 |
| | CPL | 252 | 378 | 494 | 509 | 434 |
| nonanal | CTG | 161 | 300 | 354 | 586 | 720 |
| | CPL | 882 | 1250 | 1410 | 1320 | 1180 |

SEM (pooled) hexanal 67.7, 2-pentylfuran 163, decadienal 19.8; p < 0.001 for every row. Also printed:
1-hexanol, heptanal, 1-octen-3-ol, (E)-2-nonenal, 1-octanol, decanal, 2-heptanone, 2,3-octanedione and
C11-C18 aldehydes.

Derived here from the table:

- **Decadienal in PL systems falls after 30 min.** Apparent first-order net loss in CPL:
  ln(567/289)/0.5 h = **1.35 h⁻¹** (30-60 min), ln(567/142)/1 h = 1.38 h⁻¹ (30-90), ln(567/73.1)/2 h =
  1.02 h⁻¹ (30-150); CTGPL 60-150 min 0.76 h⁻¹. Because formation continues, these are **lower bounds on a
  first-order decadienal loss constant at 100 °C** in this matrix. In CTG decadienal never falls.
- **Decadienal cannot be the main hexanal source.** CPL, 30→60 min: 278 ng decadienal lost (printed, p. 94)
  = 1.83 nmol (M 152.23); 394 ng hexanal gained (printed) = 3.93 nmol (M 100.16); octenal +15.5 ng = 0.12 nmol.
  Hexanal gained is **2.15×** the decadienal lost on a mole basis, so with at most one hexanal per
  decadienal, decadienal breakdown covers **≤ 46 %** of the hexanal rise. The thesis calls the 65 % hexanal
  rise "within reason" for decadienal degradation (p. 94); the arithmetic says most of it has another source.
  (Hexanal is approx. ng, decadienal is SIDA ng.)
- **The thesis's "4 and 40 times" is reversed.** It says decadienal and hexanal "were approximately 4 and 40
  times higher respectively in CPL ... after 30 min" (p. 95). The table gives decadienal 567/14.4 = **39×**
  and hexanal 606/147 = **4.1×**.
- **PL catalyses TG.** CTGPL carries half the PL of CPL yet matches or exceeds it from 60 min on (hexanal
  1280 v 1000; octenal 71.8 v 54.3; nonenal 34.4 v 17.7, printed p. 95). CTG hexanal is still rising at
  150 min (3.0× from 30 min); PL-containing samples plateau by 60-90 min.
- **2-Pentylfuran keeps accumulating after hexanal stops.** Molar 2-pentylfuran/hexanal (M 138.21 and
  100.16): CTG 0.116 → 0.475 and CPL 0.558 → 2.19 from 30 to 150 min. In headspace SPME at 60 °C the
  absolute ratio is extraction-weighted, but the **4-fold rise within a series** is not.
- Decadienal in CPL at 30 min, per lipid: 567 ng / 20 mg = **28 µg/g PL**.

### 2b. Lipid class and the Maillard side, LN-SDE, 30 min at 100 °C (Table 3-8, pp. 107-108; Table 4-3, pp. 146-147)

| compound (approx. ng) | CTG | CTGPL | CPL | PL/TG, chicken (derived) |
|---|---|---|---|---|
| 5-ethyl-2,4-dimethyl-3-thiazoline (I) | 28.2 | 45.2 | 63.6 | 2.26 |
| 5-ethyl-2,4-dimethyl-3-thiazoline (II) | 54.3 | 95.7 | 107 | 1.97 |
| tetramethylpyrazine | 70.9 | 99.3 | 144 | 2.03 |
| 2-ethyl-3,5-dimethylpyrazine | 275 | 342 | 346 | 1.26 |
| 1,3-benzothiazole | 113 | 307 | 480 | 4.25 |
| 3-(methylthio)propanal | 252 | 239 | 272 | 1.08 |
| dimethyl disulfide | 774 | 1030 | 655 | 0.85 |
| dimethyl trisulfide | 1450 | 2120 | 1430 | 0.99 |
| hexanal | 3870 | 3570 | 4760 | 1.23 |
| (E,E)-2,4-decadienal | 379 | 353 | 399 | 1.05 |

These rows are identical in Tables 3-8 and 4-3. **2-Methyl-3-furanthiol, 2-furylmethanethiol and
2-mercapto-3-pentanone were not detected** in any reconstituted sample, "likely ... present in extremely low
quantities due to the small sample size" (p. 86); no alkylthiophene, 2-alkylthiazole with a lipid-length
chain, or 2-pentylpyridine appears in any table. So the thesis **cannot** test MFT/FFT suppression by lipid.
What it shows is that phospholipid **raises** a lipid–Maillard crossover product: about twice the
5-ethyl-2,4-dimethyl-3-thiazoline (proposed from 3-mercapto-2-pentanone + acetaldehyde, Fig. 2-1, p. 64),
twice tetramethylpyrazine and four times benzothiazole, while the Strecker and methionine-derived products
(methional, DMDS, DMTS, 2-phenylacetaldehyde) are flat across lipid class (methional and DMDS: class effect
n.s., Table 4-3).

**A transcription conflict between the two tables.** For many lipid-derived rows the CTGPL and CPL columns
are swapped between Table 3-8 and Table 4-3 (same data): 2-pentylfuran 1110/733 (3-8) against 733/1110 (4-3);
likewise 1-pentanol, (E)-2-hexenal, (E)-2-(2-pentenyl)furan, (E)-2-heptenal, 1-octen-3-ol, (E)-2-octenal,
(E)-2-octen-1-ol, (E)-3-nonen-2-one, (E)-2-decenal, 1-octanol, 1-nonanol, all six ketones and C14-C16
aldehydes. Hexanal, heptanal, octanal, nonanal, decanal, 1-octen-3-one, (E)-2-nonenal, decadienal and every
Maillard row agree. Table 4-3's order matches the Ch. 3 text ("1-pentanol and E-2-(2-pentenyl)furan were
present in significantly higher quantities in CPL samples", p. 88), so **Table 3-8 is the likelier error**;
which is right is not decidable from the thesis. Values below use Table 4-3.

### 2c. Sunflower and soy lipids as substitutes (Table 4-2, pp. 144-145; Table 4-3, pp. 146-147; Table 4-4, p. 148)

| approx. ng | CTG | CPL | SFTG | SFPL | SYTG | SYPL |
|---|---|---|---|---|---|---|
| 5-ethyl-2,4-dimethyl-3-thiazoline (I) | 28.2 | 63.6 | 35.5 | 60.7 | 29.1 | 52.7 |
| 5-ethyl-2,4-dimethyl-3-thiazoline (II) | 54.3 | 107 | 57.2 | 119 | 48.9 | 119 |
| tetramethylpyrazine | 70.9 | 144 | 61.6 | 126 | 61.4 | 69.8 |
| 2-ethyl-3,5-dimethylpyrazine | 275 | 346 | 142 | 162 | 153 | 161 |
| 1,3-benzothiazole | 113 | 480 | 18.3 | 32.6 | n.d. | n.d. |
| dimethyl trisulfide | 1450 | 1430 | 1300 | 860 | 1230 | 1870 |
| 3-(methylthio)propanal | 252 | 272 | 186 | 174 | 178 | 158 |
| hexanal | 3870 | 4760 | 5120 | 6750 | 4800 | 4490 |
| (E,E)-2,4-decadienal | 379 | 399 | 1960 | 374 | 574 | 399 |
| 2-pentylfuran | 425 | 1110 | 522 | 887 | 554 | 604 |
| (E)-2-(2-pentenyl)furan | 9.91 | 47.8 | 14.0 | 13.0 | 14.4 | 14.8 |
| 1-octen-3-ol | 709 | 993 | 692 | 443 | 708 | 495 |
| hexadecanal | 47800 | 233000 | 16700 | 34300 | 19900 | 23000 |

Fatty acids (Table 4-2, % of measured FA): C18:2 c9,c12 is 23.0 / 21.8 (TG^C / PL^C), 57.3 / 58.6
(sunflower), 51.8 / 51.5 (soy); C18:3n-3 3.07 / 1.69, n.d. / 0.18, 5.53 / 4.75; C20:4n-6 0.47 / 7.65 in
chicken and n.d. in every plant lipid; total PUFA 27.4 / 36.3, 57.3 / 58.9, 57.4 / 56.3. Table 4-4 is a
literature compilation of PL classes (not measured): PC 46-61 / 31-45 / 22-35 %, PE 15-28 / 14-26 / 17-26 %,
PI 6-11 / 14-32 / 16-18 %, PA n.r. / 2-6 / 6-11 %, SM 3-8 % in chicken only.

Readings (ratios derived here):

- **The phospholipid thiazoline uplift transfers to plant PL.** PL/TG for thiazoline (I)/(II): chicken
  2.26/1.97, sunflower 1.71/2.08, soy 1.81/2.43. Plant PL give 0.83-1.11× the chicken-PL amount.
- **The pyrazine uplift transfers to sunflower, not soy.** Tetramethylpyrazine PL/TG: chicken 2.03,
  sunflower 2.05, soy 1.14. The thesis says "about twice the amount of tetramethylpyrazine in (TG)PL samples
  of chicken or soy origin" (p. 134); the table supports chicken or **sunflower**. 2-Ethyl-3,5-dimethylpyrazine
  is about half of chicken with any plant lipid (0.47×).
- **Benzothiazole is the largest source gap**: 480 (CPL) against 32.6 (SFPL) and n.d. (soy).
- **n-6 aldehydes do not scale with linoleate.** Plant lipids carry 2.2-2.7× chicken's C18:2, yet SFPL
  hexanal is 1.42× CPL and SYPL 0.94×; 1-octen-3-ol, (E)-2-nonenal and the C20-derived products are higher
  with chicken PL, which the thesis puts on PL^C's C20:4 and possible plant antioxidants (pp. 135-136).
  Sunflower TG is the decadienal outlier (1960, 5.2× CTG).
- The thesis's choice: **sunflower PL is "a more befitting alternative"** to chicken PL by PCA (p. 138), but
  SF(TG)PL still carry more hexanal, heptanal and decadienal and less DMDS, DMTS and
  2-ethyl-3,5-dimethylpyrazine than the chicken target; it suggests tocopherol to slow lipid oxidation.

### 2d. GC-O, boiled beef and chicken broths (Tables 2-2 and 2-3, pp. 57-62; Table 2-4, p. 64)

| odorant | beef: intensity (/60), frequency (/6) | chicken: intensity, frequency |
|---|---|---|
| 2-methyl-3-furanthiol | 41, 6 ("roasted, porridge oats, sulfur, faecal") | 32, 5 ("cooked meat, roasted, porridge oats") |
| 2-furylmethanethiol + 2-mercapto-3-pentanone | 39, 5 | 34, 5 |
| methional | 45.5, 6 | 33, 5 |
| (Z)-4-heptenal | 48.5, 6 | 38, 6 |
| 2-acetyl-1-pyrroline | 42, 6 | 39, 6 |
| hexanal | 36, 6 | 32, 6 |
| 5-ethyl-2,4-dimethyl-3-thiazoline (II) | 33, 6 ("fatty, beef fat") | 20, 3 ("fatty, grilled meat, savoury") |
| (E,E)-2,4-decadienal | 5, 1 | 9, 2 |
| bis(2-methyl-3-furyl) disulfide | 5, 1 (beef only) | not detected |
| 2,4-dimethyl-3-thiazoline | not detected | 27, 4 ("meaty, brothy") |
| benzothiazole | not detected | 18, 4 ("savoury, chicken fat, nutty") |

58 odour-active volatiles in beef and 65 in chicken (p. 43). 2-Pentylfuran is **not** among the odorants
of either broth; nor are any 2-alkylthiophenes. Table 2-4 prints literature odour thresholds in **µg/kg in
water**: (E)-2-heptenal 13, (E,E)-2,4-heptadienal 0.032, (E,E)-2,4-nonadienal 0.06, (E,Z)-2,6-nonadienal
0.0045-0.02, (E,E)-2,4-decadienal 0.027-0.07 (secondary: Flaig 2020, Milo & Grosch 1993, Buttery 1988).

## 3. What it means for the model and for formulation

**2-pentylfuran (`PENTYLFURAN_PER_HEXANAL`).** The lane forms 2-pentylfuran as hexanal × **0.160**
(`src/kinetic_core/parameters_lipid_b28.py`, `"linoleate_autoxidised": 2.4 / 15.0`, Frankel 1981 Table III,
a peak-area ratio). Yeo's LN-SDE molar ratios (derived, ×100.16/138.21): CTG 0.080, CTGPL 0.149, CPL 0.169
(Table 4-3), SFTG 0.074, SFPL 0.095, SYTG 0.084, SYPL 0.097, i.e. **0.46-1.06× the shipped ratio**, with the
plant lipids at 0.46-0.61×. That **tests** the constant in a cooked meat matrix and passes within a factor
of ~2. It does not pin it: both numbers are RRF semi-quantities, and 2a shows the ratio is **not a constant
in time** (4-fold rise over 30-150 min at 100 °C), which a fixed per-flux ratio cannot reproduce on long
cooks.

**Decadienal → hexanal (`DECADIENAL_TO_HEXANAL_GAP`).** The code says "No rate for it exists in the corpus".
Yeo gives the first observation bearing on it: a SIDA-quantified net loss of decadienal at 100 °C of
**≥ 1.0-1.4 h⁻¹** apparent first order in a PL-rich chicken matrix, and a stoichiometric ceiling showing
decadienal breakdown supplies **≤ 46 %** of the concurrent hexanal rise (2a). Neither is a rate constant for
the edge (loss of decadienal is not attributed to hexanal), so the gap stays declared; the bound is usable as
a check if the edge is ever enabled.

**Hydroperoxide decomposition and formation.** The live lane anchor is `k_LOOH_decomp` = 6.0e-3 h⁻¹ at 25 °C
(`parameters_lipid.py`, Schroen & Berton-Carabin 2022 Table 1) with `lipid.q10` a `declared_band` [2.0, 3.0],
centre 2.449 (`results/validation/core_prediction_uncertainty.json`); at 100 °C that is 1.09 / 4.97 /
22.7 h⁻¹ at Q10 2 / 2.449 / 3 (derived, factor Q10^7.5). Yeo prints no peroxide value and no
hydroperoxide, so nothing here pins that coordinate. Two qualitative tests: (i) hexanal in CTG still rises
3× between 30 and 150 min while PL samples plateau by 60-90 min, so the effective rate at 100 °C depends on
**lipid class and dispersion** at fixed temperature by an order of magnitude (decadienal 39× at 30 min), a
dimension the lane's carriers (PV, lipid fraction, FA profile) do not represent; (ii) continued accumulation
over 150 min after an oxygen top-up is formation during heating, which `LOOH_FORMATION_GAP` declares
unmodelled.

**Lane coupling.** `lipid.lane_coupling_verdict` sums the lipid and Maillard lanes with no cross term.
Methional and DMDS, the two engine species measured here, show no lipid-class effect (2b), which is
consistent with the direct sum for those species. The crossover products that do move (2× 3-thiazoline,
2× tetramethylpyrazine, 4× benzothiazole) are not engine species. Nothing here tests MFT or FFT.

**Formulation.** At 2 % lipid in a boiled matrix, swapping triglyceride for phospholipid is the larger lever
than swapping the source: it doubles the fatty/grilled-meat 3-thiazoline and tetramethylpyrazine with
chicken, sunflower or (for the thiazoline) soy PL. Sunflower lecithin reproduces chicken PL most closely on
this panel; soy lecithin gives the thiazoline but not the pyrazine. Neither plant PL gives benzothiazole or
the C20:4-derived odorants, and both push hexanal-family notes, so an antioxidant (the thesis suggests
tocopherol) or a lower dose is the open variable. None of this was sensory-tested in a plant matrix.

## What it does not give

- No rate constant, activation energy or temperature series: one temperature (100 °C), five times, one
  matrix; amounts are headspace-SPME semi-quantities except SIDA decadienal.
- No MFT, FFT, alkylthiophene, long-chain alkylthiazole or pentylpyridine data in any lipid experiment; no
  cysteine or ribose model system.
- No plant protein matrix: every substitution was made in **defatted chicken**, which had lost 80 % of its
  sugars and 40 % of its amino acids to the extraction.
- No peroxide value, hydroperoxide or lipid-oxidation state of any lipid; no measured PL class composition
  (Table 4-4 is literature).
- No quantitative GC-O (no AEDA/FD factors, no OAVs); thresholds in Table 2-4 are second-hand.
- No support for the `yeo_mottram_2023_soy_lecithin_thiophene_uplift` numbers (2.4-fold thiophene uplift,
  0.5-1 % lecithin window).
