# Wright 2006 — EXTRACTION (four commercial Cheddar-whey protein isolates rehydrated at 10 % solids, two with and two without a cabbage off-flavour; SAFE solvent extraction + GCO/AEDA on 2 sniffers; DMTS quantified by HS-SPME-GC-MS on two fibre chemistries; **orthonasal detection thresholds for dimethyl trisulfide measured in deodorised water AND in 10 % WPI by ASTM ascending forced choice on 80 panelists**; model-system confirmation)

### THE PAIRED THRESHOLD THE MATRIX LAYER HAS BEEN WAITING FOR, ON A SULFUR COMPOUND — AND ITS SIZE IS THE PROBLEM: **DMTS 0.07 ± 1.28 ppt in water against 0.80 ± 0.45 ppb in 10 % WPI, both confirmed, both from the same 80-panelist orthonasal 3-AFC study. That is a shift of 11 429x (mine)** — an order of magnitude beyond anything `REVERSIBLE_BINDING` can express (a `K_g` proxy of **114 L/g (mine)**, ~300x the largest constant in the table), and the water leg is itself **143x below the published water threshold the paper cites (mine)**. The matrix correction is real and enormous: applied to this paper's own samples it moves DMTS from OAV ~6 100-46 400 on the water threshold to **2.4 and 4.1 for the two cabbage-flavoured isolates and 0.55 and 0.54 for the two clean ones (mine)** — the only reading under which the sensory panel's verdict and the chemistry agree.

**Source on disk:** `data/articles/wright2006.pdf` (5 pp., J. Food Science 71(2) 2006, pp. C86-C90,
section "JFS C: Food Chemistry and Toxicology"). Read from the `pdftotext -layout` text layer
(`scratchpad/wright2006.txt`), **whole file**. **Tables 1, 2, 3 and 4 came through clean and are
re-typed in full below — that is every table in the paper; there is no supplementary material.**
**Figure 1** (WPI model-system cabbage flavour and aroma intensity and similarity scores) is an image
and is **figure-only**: its caption is typed below, but **not one intensity or similarity score from
the model-system experiment is printed as a number anywhere in this article**. Repo status before
this dossier: Wright 2006 is cited nowhere in `src/kinetic_core/` and has no extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Characterization of a Cabbage Off-flavor in Whey Protein Isolate" |
| Authors | Joy M. Wright, Mary E. Carunchia Whetstine, R. Evan Miracle, MaryAnne Drake (corresponding, maryanne_drake@ncsu.edu) — Dept. Food Science, **Southeast Dairy Foods Research Center, North Carolina State Univ., Raleigh, NC 27695** |
| Venue | **Journal of Food Science, Vol. 71, Nr. 2, 2006, pp. C86-C90.** MS 20050644, **submitted 25 October 2005, revised 21 November 2005, accepted 29 November 2005**; "Published on Web 2/27/2006"; "© 2006 Institute of Food Technologists" |
| DOI | **NO DOI IS PRINTED ANYWHERE IN THIS PDF.** There is no DOI on the title page, in the footer of any of the five pages, in the manuscript-number block, or in the acknowledgments. The only identifiers printed are the manuscript number **MS 20050644** and the departmental manuscript number **FSR 06-05**. Cite by volume/issue/pages |
| Funding | **California Dairy Research Foundation.** "Manuscript FSR 06-05 of the Dept. of Food Science, North Carolina State Univ." |
| Ethics | **Not stated.** No IRB number appears, despite an 80-person human sensory study |
| The compound | **Dimethyl trisulfide (DMTS)**, Sigma-Aldrich (Allentown, Pa.), **purity as purchased reported 98.5 %, re-evaluated on receipt by HS-SPME at 95 %**; stored at −20 C, all experiments **within 2 wk of purchase** |
| The matrix | **Four commercial whey protein isolates**, all **< 3 months old**, all manufactured from **Cheddar cheese whey**, from **5 domestic (US) sources** (five collected, four selected). Two with cabbage flavour, two without, chosen by a 7-panelist screen. **Rehydrated at 10 % solids (w/v) in deodorised deionised water** |
| Naming | **BET** = best estimate threshold; **group threshold** = geometric mean of individual BETs; **log3 FD** = flavour dilution factor from AEDA; **NB / AC** = neutral-basic / acidic extract fraction; **ND** = not detected |
| Companions on disk | `bornhorst2017_extraction.md` / `bornhorst2017b_extraction.md` (whey), `leksrisompong2010_extraction.md` (the same Drake laboratory's caseinate threshold/partition work — the corpus's other paired-threshold source), `anantharamkrishnan2020b_extraction.md` (the DMDS/DMTS adduct table), `k4b_paired_thresholds_and_browning.md` (the repository's paired-threshold register) |

## 1. Why it matters

**Confirming the numbers the task asked about, one at a time.**

| item asked | printed value | where |
|---|---|---|
| DMTS threshold in water | **0.07 ± 1.28** | Table 4, "Experimental threshold / Water (ppt)"; Abstract; Results p. C90 |
| its unit | **parts per trillion (ppt)** — spelled out in the Abstract as "0.07 ± 1.28 parts per trillion (ppt)" and abbreviated in the Table 4 column head | Abstract, Table 4 |
| DMTS threshold in WPI | **0.80 ± 0.45** | Table 4, "Experimental threshold / WPI (ppb)"; Abstract; Results p. C90 |
| its unit | **parts per billion (ppb)** — Abstract: "0.80 ± 0.45 ppb"; Table 4 column head "WPI (ppb)" | Abstract, Table 4 |
| panel size | **80** — Table 4 footnote (b): *"Best estimate threshold from 80 panelists."* Methods: *"Panelists (n = 80) were given these concentrations in a series..."* | Table 4 fn (b), Methods p. C88 |
| orthonasal or retronasal | **ORTHONASAL.** Abstract: *"Orthonasal thresholds of DMTS in deodorized water and WPI were determined by ascending forced choice analysis"*. Methods: *"Subjects were told to open the soufflé cups and to briefly sniff the headspace of each cup in the series."* Results: *"we determined orthonasal thresholds of DMTS in water and WPI"*. Table 4's own column header block is titled "Quantification and sensory **orthonasal** threshold values" | Abstract, Methods, Results, Table 4 title |
| detection or recognition | **DETECTION.** The task is *"to choose the 1 different sample from the 3 they were presented"*, against **2 blanks** (2 deodorised-water blanks for the water series, 2 WPI blanks for the WPI series). That is a difference-from-blank discrimination, not an identification of the odour quality, so it is a **detection threshold**. The paper calls it a "best estimate threshold" throughout and never uses the word "recognition" | Methods p. C88 |

**All four numbers confirmed as printed. Both legs measured in this study, by the same protocol, the
same panel of 80, on different days (*"Threshold testing was conducted for each compound on
different days"*).**

**The shift, computed here.** 0.80 ppb = 800 ppt. **800 / 0.07 = 11 428.6, i.e. the WPI threshold is
11 429x the water threshold (mine).** The paper never computes this ratio and never states it; it
says only *"The difference between water and WPI is expected, due to the components of the food
matrix (WPI) interacting with flavor compounds and impacting flavor release into the headspace."*

**Why this is the most important and the most dangerous paired threshold in the corpus.**

*Important*, for three reasons the repository can name precisely:

1. **It is on a SULFUR compound.** The matrix layer's `REVERSIBLE_BINDING` table holds 21 constants
   — esters, ketones, aldehydes, one lactone, one alcohol, furaneol — and **not one thiol, disulfide
   or trisulfide**. The repository's sulfur lane produces 2-methyl-3-furanthiol, 2-furfurylthiol and
   their disulfides, and every matrix shift it currently applies to them is 1.0 by absence of
   evidence. Wright is the first corpus source with a **measured, paired, same-panel water-vs-matrix
   threshold on a polysulfide**, and `trisulfide` is a class the module already carries in
   `ADDUCT_POSITIVE_CLASSES`.
2. **The correction changes the sign of the answer, and this paper proves it on its own samples.**
   Applying the water threshold to the four measured DMTS levels gives OAVs of **27 714 / 46 429 /
   6286 / 6143 (mine)** — all four isolates hugely over threshold, including the two the trained
   panel scored as having **no cabbage note at all (ND, Table 2)**. Applying the WPI threshold gives
   **2.43 / 4.06 / 0.55 / 0.54 (mine)** — the two cabbage isolates above 1, the two clean isolates
   below 1, exactly matching the panel. **This is a clean, within-paper demonstration that a
   water-threshold OAV can be wrong by four orders of magnitude and can flip a qualitative verdict.**
   It is direct evidence for the matrix layer's existence and against ever shipping a water-threshold
   OAV as if it described a protein matrix.
3. **It is the same laboratory as Leksrisompong 2010** (Drake, NC State), whose caseinate threshold
   and partition data already sit in `REVERSIBLE_BINDING` as three FIT rows, so the panel culture,
   the ASTM protocol and the "our experimental thresholds differ from published ones" habit are
   shared — which is both a consistency argument and a common-mode-error warning.

*Dangerous*, for four reasons that are equally precise:

1. **The magnitude is outside anything the layer can represent.** In the registry's form
   `K_g = (K_water/K_matrix − 1)/protein_g_per_L`, substituting the threshold ratio for the partition
   ratio at 100 g/L of WPI solids gives **114.3 L/g (mine)**; at a nominal 90 % protein it is
   **127.0 L/g (mine)**. The largest constant currently shipped is `kg_t_2_octenal_pea` at
   **3.834e-1 L/g**, and that one is quarantined. **This is roughly 300x the top of the shipped range
   (mine).** Either DMTS interacts with whey protein by a mechanism categorically unlike everything
   else in the table, or one of the two legs is wrong. Amendment 6 ruling 2 caps reversible binding
   at **~25 % of an observed log-shift**; a log-shift of 4.06 decades is 16x the entire
   log-shift budget the layer allows itself. **This number cannot be entered as a binding constant.**
2. **The water leg is 143x below the published water threshold the paper itself cites (mine).**
   Rychlik, Schieberle & Grosch 1998 give **0.01 ppb = 10 ppt**; Wright measures **0.07 ppt**.
   The paper acknowledges this head-on: *"Our results indicated that the orthonasal threshold for
   DMTS in water was much lower than previously reported—0.07 ppt"*, and pre-defends it with
   *"Thresholds can vary widely depending on testing procedure, number of panelists, and matrix
   used"* and *"Previous research in our laboratory has indicated large differences between
   previously published thresholds and experimental thresholds"*. **A 143x self-acknowledged
   disagreement on the WATER leg means the shift of 11 429x is not a matrix effect of known size.**
   Against the cited literature water value the WPI shift is only **80x (mine)** — still large, and
   143x smaller than the paper's own internal ratio. **Which baseline is used changes the answer by
   two orders of magnitude, and the repository must not silently pick one.**
3. **The ± quantities in Table 4 are not defined and cannot both be arithmetic standard
   deviations.** Water is **0.07 ± 1.28 ppt** — the ± is **18x the mean (mine)** and the lower bound
   is deeply negative, which is impossible for a concentration. WPI is **0.80 ± 0.45 ppb**, a
   perfectly ordinary 56 % arithmetic SD. A group BET is a **geometric** mean, so a natural companion
   statistic is a **geometric standard deviation** — a multiplicative factor, for which 1.28 is
   entirely sensible and 0.45 (being less than 1) is not. **So the two cells appear to carry two
   different statistics under one column heading, and the paper never says what either is**
   (Flags 1). Until that is resolved, no interval can be put on the shift.
4. **The two legs differ in more than the matrix.** The water blanks are deodorised water; the WPI
   blanks are 10 % WPI *"determined to be free of cabbage aroma and flavor by descriptive sensory
   analysis"* — i.e. free of the *target* note, but carrying the full whey background that Table 2
   scores at overall aroma intensity 2.27-2.83, cardboard 1.67-2.50, sweet aromatic ~2.0 and
   astringent 2.50. **A detection threshold measured against a blank that already has aroma is
   elevated by masking as well as by binding, and the two cannot be separated here** (Flags 3).

**What the repository should actually do with this.** Enter it as a **`threshold` pair with both legs
and the full provenance, and as a `measured_bound` on the sulfur matrix effect — never as a
`binding_constant` and never as a `K_g`.** The defensible statement is: *"on the one paired
measurement that exists for a polysulfide in a 10 % animal-protein matrix, the orthonasal detection
threshold was 11 429x higher in the matrix than in water by the same panel, and 80x higher than the
literature water value; both comparisons are far outside the range any reversible-binding term in
this model can produce, and the discrepancy is unexplained."* That is exactly the "UNEXPLAINED
RESIDUAL, quantified per compound" output the module's docstring says is the honest one.

## 2. Methods as they matter to a model

- **The pot (thresholds).** **15 mL** of each dilution poured into **clean, labelled 56-mL plastic
  cups**, sniffed as headspace. Cups opened by the subject.
- **The pot (descriptive sensory).** **30 mL** of rehydrated WPI in **56-mL plastic Sweetheart cups**
  with lids and 3-digit codes, stored at **5 C for less than 18 h**, **removed from refrigeration
  1 h before tasting**.
- **The pot (HS-SPME).** **20 mL** of rehydrated WPI + a stir bar + **1 g NaCl** in a **40-mL amber
  glass SPME vial**, sealed with a PTFE/silicone septum. So beta ≈ 1 (mine).
- **The pot (solvent extraction).** **40 g** of WPI rehydrated at 10 % solids, spiked with 20 uL of
  internal standards, **120 g NaCl** added to break the emulsion, extracted **3 x 50 mL ethyl
  ether** in 250-mL Teflon bottles, 30 min on a Roto Mix at top speed, centrifuged **735 x g for
  10 min**.
- **Protein loading — and this is the number every per-gram construction divides by.**
  **"reconstituted at 10% solids (w/v) in deodorized deionized water"** = **100 g of WPI powder per
  litre**. **The protein content of the isolate is NEVER stated** — no N x 6.38, no supplier
  specification, no assay. Commercial WPI is conventionally 90 % protein or better, so the protein
  loading is **~90 g/L (mine, assumed)**, but that is an assumption this paper does not support
  (Flags 5). **At 100 g/L this is by far the highest protein loading in the corpus's threshold work**
  — three times skim milk's 33.9 g/L, ten times the caseinate 10 g/L and pea 10 g/L rows.
- **pH.** **Never stated**, for any preparation.
- **Temperature.** Sensory: samples out of 5 C refrigeration **1 h before tasting**, so ambient.
  Threshold cups: no temperature stated at all. SPME: **48 C, stirred 30 min** equilibration, then
  fibre exposed **30 min** at a depth of 3.8 cm. SAFE: **50 C** circulating bath, 2 h under
  ~1e-5 to 1e-3 Torr vacuum.
- **Light control.** *"To avoid light oxidation, proteins were prepared for sensory and instrumental
  analysis with overhead lights off."* — worth noting for a paper about sulfur volatiles.
- **Threshold method and its family — the details that decide transferability.**
  - **Protocol: ASTM ascending forced choice method of limits.** Methods cites it as **"procedure
    E679-79 (ASTM 1992)"** while the reference list gives **"E-679-91. In: Annual book of standards.
    15.07. Philadelphia, Pa.: ASTM. p 35-9"** with year 1992 (Flags 8).
  - **With the correction factor of Lawless, Harono & Hernandez 2000**, J Sens Stud 15:437-47.
  - **Route: ORTHONASAL** (sniff the headspace of an opened cup).
  - **Criterion: DETECTION** (pick the odd sample out of three, two of which are blanks).
  - **Design: 3-AFC ascending**, **serial dilution factor of 3**, **five ascending series tested each
    time**, **each series presented in randomised order**, **n = 80 panelists**. Subjects **briefly
    instructed before testing** (i.e. **not a trained panel** — this is a naive-consumer threshold,
    unlike the 7-panelist trained descriptive panel used elsewhere in the paper; Flags 4).
  - **Palate/nose clearing: 1 min rest between each set of 3, and subjects instructed to sniff their
    own sleeve** between cups.
  - **Certainty judgment collected ("sure/not sure").**
  - **Individual BET rule, verbatim**: *"the geometric mean of the last concentration with an
    incorrect response and the 1st concentration with a correct response except for the following
    sequence: if the subject indicated a 'not sure' response for the correct choice, that
    concentration was increased by a factor of 1.41, to adjust for the possibility of a chance
    correct response (Lawless and others 2000)."*
  - **Group threshold = geometric mean of the individual BETs.**
  - **The two legs were run on different days.**
  - **The carrier: DMTS stock in METHANOL**, dosed into either water or WPI. *"Preliminary threshold
    analysis of blank solutions with added methanol indicated that methanol provided no discernable
    aroma."* **The volume of methanol per cup is never stated** (Flags 6).
- **Quantitation method: HS-SPME-GC-MS, and it is the weak leg of the paper.** WPI rehydrated at
  10 % solids, 100 mL spiked with **2.5 uL of 1-pentanol internal standard (stock 40.7 ppb)**, then
  20 mL taken into the vial. **Two fibre chemistries run in parallel**: 2-phase **75-um
  Carboxen-PDMS** and 3-phase **2 cm 50/30 um DVB/Carboxen/PDMS StableFlex** (Supelco), because
  Supelco recommended the 2-phase fibre as more responsive to sulfur. **Each sample run in triplicate
  on both fibres = 6 reps.** Desorbed **250 C for 5 min**, 7.6 cm depth, SPME inlet liner.
  GC-MS: Varian CP-3380 GC / **Saturn 2000 ion trap**, Rtx-5 30 m x 0.25 mm x 0.25 um, helium
  1 mL/min, oven 40 -> 250 C at 8 C/min with 5 min holds, transfer line 120 C, manifold 80 C, ion
  trap 150 C, automatic EI, m/z 35-350, EM 2135 V.
  **Calibration: external, by spiking DEODORISED DEIONISED WATER — not WPI — with DMTS over
  500 to 2000 ppt**, extracted and analysed the same way, standard curve of concentration against
  DMTS/IS peak-area ratio. **R^2 = 0.995 (3-phase) and 0.998 (2-phase).** The two fibres gave
  concentrations that *"were not different (P > 0.05) (data not shown)"* — though the text also says
  *"Our results did show a statistically significant difference between the 2 fibers"*, which
  contradicts it (Flags 7).
  **The calibration range (500-2000 ppt = 0.5-2.0 ppb) does not cover two of the four measured
  samples** (0.44 and 0.43 ppb are below it) and barely covers a third (1.94 ppb is inside;
  3.25 ppb is **above** it) — **so three of the four reported concentrations are extrapolations
  (mine)** (Flags 2).
- **GCO/AEDA.** HP 5890 series II with FID + sniffing port, splitless, 2 uL injected, **both a polar
  Rtx-Wax and a non-polar Rtx-5** (30 m x 0.25 mm x 0.25 um), effluent split **1:1** between FID and
  sniffing port through 1 m x 0.25 mm deactivated fused silica, oven 40 -> 200 C at 10 C/min with a
  3 min initial and 20 min final hold, FID and sniffing port at 250 C, humidified air 30 mL/min.
  **AEDA: stepwise 1:3 dilution with diethyl ether until no odorant detected; highest dilution
  reported as log3 FD.** Neutral/basic fractions sniffed on the **non-polar** column, acidic
  fractions on the **polar** column. **Two experienced sniffers**, each with >50 h training on dairy
  extracts, **each extract in duplicate**. Only **two of the four** WPI were extracted (one cabbage,
  one non-cabbage).
- **Descriptive sensory.** **7 experienced panelists (6 female, 1 male)**, each with **>150 h**
  general descriptive experience and **40 h specifically on dried whey proteins**, **15-point
  universal Spectrum intensity scale**, sensory language of Drake 2003 / Carunchia Whetstine 2005.
  **Each product evaluated by each panelist in duplicate**, randomised balanced block.
- **Model system.** DMTS at **1.05, 1.94 and 3.25 ppb** added to a commercial WPI without cabbage
  flavour, *"all of which were within the concentration range identified in WPI exhibiting cabbage
  flavor"*. **Methanol as the carrier for aroma evaluation** (nearly odourless), **95 % ethanol for
  flavour evaluation** (safe to consume). Blanks with the carrier and no DMTS. Evaluated by the same
  descriptive procedure, plus a **10-point similarity scale (1 = very different, 10 = identical to
  reference)** against a WPI that naturally had the cabbage flavour.
- **Statistics.** ANOVA by PROC GLM (SAS 8.2), **Fisher's LSD** post hoc.

## 3. Tables re-typed

### Table 1 (p. C87). "Whey protein sensory language and references^a"

| Term | Definition | References |
|---|---|---|
| **Overall aroma intensity** | The overall orthonasal aroma impact of the rehydrated sample | — |
| *Flavors (evaluated in the mouth)* | | |
| Sweet aromatic | The sweet aromatics associated with dairy products (diacetyl is 1 example) | Diacetyl (2,3-butanedione), mild cheddar, or Colby-jack shreds |
| Cardboard | Aromatics associated with wet cardboard | Pentanal, cardboard in water |
| Brothy | Aromatics associated with vegetable stock and boiled potatoes | Methional, broth from canned potatoes |
| Soapy | Aromatics associated with medium-chain fatty acids and soaps | Decanoic acid, unscented plain bar soap in water |
| **Cabbage** | **Sulfurous aromatic associated with cooked cruciferous vegetables** | **Dimethyl trisulfide, boiled fresh cut cabbage** |
| Bitter | Fundamental taste sensation elicited by caffeine, quinine | Caffeine (0.08 % in water)^b |
| Astringency | Drying sensation on the tongue and oral cavity surfaces | Black tea, alum |

Footnotes: *"(a) Adapted from Drake and others (2003) and Carunchia Whetstine and others (2005).
(b) Universal scale reference as described in Meilgaard and others (1999)."* **Note that "Overall
aroma intensity" is orthonasal and every other term is evaluated in the mouth**, so the "Cabbage"
scores in Table 2 are **in-mouth (retronasal) flavour** while the threshold in Table 4 is
**orthonasal** (Flags 9).

### Table 2 (p. C89). "Descriptive sensory analysis of whey protein isolates^a"

| Attribute | Non-cabbage WPI 1 | Non-cabbage WPI 2 | Cabbage WPI 1 | Cabbage WPI 2 |
|---|---:|---:|---:|---:|
| Overall aroma intensity | 2.27 c `[M]` | 2.83 c `[M]` | 3.33 b `[M]` | **4.67 a** `[M]` |
| Sweet aromatic | 2.05 a `[M]` | 2.00 a `[M]` | 1.83 a `[M]` | 1.83 a `[M]` |
| Cardboard | 1.67 b `[M]` | 2.50 a `[M]` | 3.00 a `[M]` | 2.67 a `[M]` |
| Brothy | ND | ND | ND | 2.33 `[M]` |
| **Cabbage** | **ND** | **ND** | **3.17 b** `[M]` | **4.00 a** `[M]` |
| Soapy | 1.17 a `[M]` | 1.17 a `[M]` | ND | ND |
| Bitter | 1.50 `[M]` | ND | ND | ND |
| Astringent | 2.50 a `[M]` | 2.50 a `[M]` | 1.67 a `[M]` | 2.33 a `[M]` |

Footnote: *"(a) Intensities are scored on a 15-point universal Spectrum™ scale where 0 = none and
15 = very high. Most dairy product flavors fall between 0 and 10. Means in a row followed by
different letters denote differences among samples (P < 0.05). ND = not detected."*

**The key cells: cabbage is ND in both clean isolates and 3.17 / 4.00 in the two off-flavoured ones,
and the higher score belongs to the isolate with the higher DMTS (3.25 vs 1.94 ppb).** Overall aroma
intensity orders the same way. **Note "Brothy 2.33" appears only in Cabbage WPI 2 and carries no
significance letter**, and its Table 1 reference compound is **methional** — which Table 3 shows at
log3 FD 2 in the cabbage extract (Flags 10).

### Table 3 (pp. C89). "Aroma extract dilution analysis of whey protein isolate (WPI) with and without the cabbage off-flavor"

Log3 FD^b (post-peak intensity^c):

| Nr | Compound | Fraction | Odor^a | Non-cabbage WPI 1 | Cabbage WPI 1 | RI DB-5 | RI DB-Wax | Method of Identification^e |
|---:|---|---|---|---|---|---:|---:|---|
| 1 | 2,3 Butanedione | NB | Buttery | <1 (1.38) | ND | 680 | 955 | RI, odor |
| 2 | Phenylacetaldehyde | NB | Rosy | ND | <1 (1.25) | 1044 | 1619 | RI, odor, MS |
| **3** | **Dimethyl disulfide** | NB | **Garlic** | **1 (2.5)** | **ND** | 777 | 1071 | RI, odor, MS |
| 4 | Butanoic acid | AC | Cheesey/rancid | 3 (3) | **6 (1.5)** | 840 | 1650 | RI, odor, MS |
| 5 | Methional | NB | Potato/brothy | 4 (1.67) | 2 (3) | 923 | 1441 | RI, odor, MS |
| 6 | 2-Acetyl-1-pyrroline^g | NB | Popcorn | <1 (1.5) | 2 (2.5) | 940 | 1291 | RI, odor |
| **7** | **Dimethyl trisulfide** | NB | **Cabbage** | **2 (2)** | **<1 (2.5)** | 981 | 1369 | RI, odor, MS |
| 8 | 1-Octen-3-one | NB | Mushroom | <1 (1.25) | 2 (2.7) | 984 | 1312 | RI, odor, MS |
| 9 | 1,5 Octadienone | NB | Earthy/musty | 2 (1.5) | ND | 988 | 1341 | RI, odor |
| 10 | Octanal | NB | Citrus/green | 3 (1.5) | ND | 1005 | 1051 | RI, odor |
| 11 | Hexanoic acid | AC | Sweaty | 4 (3) | ND | 1045 | 1275 | RI, odor, MS |
| 12 | 2,5-Dimethyl-4-hydroxy-3-(2H) furanone (Furaneol) | AC | Burnt sugar | <1 (1.5) | ND | 1072 | 1819 | RI, odor |
| 13 | 2-Isobutyl-3-methoxypyrazine | NB | Bell pepper/burnt | <1 (2) | ND | 1082 | 1403 | RI, odor |
| 14 | Nonanal | NB | Fatty/citrus | 2 (3) | <1 (2.7) | 1108 | 1385 | RI, odor |
| 15 | 3-Hydroxy-4,5-dimethyl-2-(5H)-furanone (sotolon) | AC | Maple/spicy | <1 (3) | ND | 1120 | 2047 | RI, odor |
| 16 | (Z)-2-Nonenal | NB | Fatty/green | <1 (1.5) | <1 (2.5) | 1126 | — | RI, odor^h |
| 17 | 2-Phenethanol | NB | Rosy | 1 (2.5) | <1 (1.5) | 1158 | 1873 | RI, odor |
| 18 | (E,Z)-2,6-nondienal | NB | Rosy/cucumber | 1 (2) | <1 (2.5) | 1160 | 1552 | RI, odor |
| 19 | E-2-nonenal | NB | Cucumber/old books | ND | <1 (2.5) | 1178 | 1584 | RI, odor |
| 20 | (E,E)-2,4-nonadienal | NB | Fatty | 3 (2) | <1 (2.5) | 1217 | 1609 | RI, odor |
| 21 | Delta-decalactone | NB | Coconut | <1 (1.5) | <1 (1.0) | 1491 | 1996 | RI, odor |
| 22 | Acetic acid | AC | Vinegar | 2 (1.5) | **6 (2.5)** | — | 1450 | RI, odor, MS |

Footnotes, verbatim: *"(a) Odor description at the GC-sniffing port. (b) The fractions were diluted
stepwise with diethyl ether at a ratio of 1/3 (v/v) and sniffed in duplicate. The dilution procedure
was followed until sniffers detected no odorants. The highest dilution was reported as the log3
flavor dilution (FD) factor (Grosch 1993). Flavor dilution factors were determined on a DB-5MS
column for NB compounds, and on a DB-Wax column for AC compounds. (c) Post peak intensities as
determined at the GC-sniffing port (van Ruth 2001). (d) Retention indices were calculated from gas
chromatography–olfactometry (GCO) data. (e) Compounds were identified by comparison with authentic
standards on the following criteria: retention index (RI) on DB-Wax and DB-5MS columns, odor
property at the GC-sniffing port, and mass spectra in the electron impact mode. Positive
identifications indicate that mass spectral data was compared with authentic standards.
(f) ND = not detected; RI = retention index; MS = mass spectra. (g) Compound identified by comparing
RI and aroma with literature (Avsar and others 2004). (h) Compound identified by comparing RI and
aroma with literature (Carunchia Whetstine and others 2005b)."*

**Note the column headers name the columns "Non-cabbage WPI 1" and "Cabbage WPI 1"** — only two of
the four isolates were solvent-extracted. **The instrument text says the columns are Rtx-Wax and
Rtx-5; the table and footnotes say DB-Wax and DB-5MS.** These are equivalent phases from different
manufacturers; the paper uses both naming conventions (Flags 11).

**THE ROW THAT MATTERS MOST, AND IT POINTS THE WRONG WAY. DMTS (row 7) has log3 FD = 2 in the
NON-cabbage WPI and <1 in the CABBAGE WPI** — the flavour dilution factor is **higher in the sample
with LESS DMTS** (0.44 vs 1.94 ppb by SPME). The paper says so plainly: *"DMTS was found in WPI with
and without cabbage flavor, and GCO results were inconclusive"*, and then gives the reasoning that
should be quoted in full because the repository's output layer says the same thing:

> *"GCO data alone cannot always be used to identify and confirm key flavor compounds. The aroma of
> compounds as they elute from the GC is not always indicative of the flavor they contribute;
> results are not quantitative, nor do they take into account the role of the food matrix and the
> collective impact of other volatile components. Therefore, threshold and model system analyses are
> necessary for confirmation of crucial flavor compounds."*

**Dimethyl disulfide (row 3) is present in the NON-cabbage WPI (log3 FD 1) and NOT DETECTED in the
cabbage WPI**, and is described at the port as **garlic**, not cabbage.

### Table 4 (p. C90). "Quantification and sensory orthonasal threshold values of dimethyl trisulfide"

| Compound | RI on DB-5 column | Cabbage WPI 1, conc. (ppb = ng/g)^a | Cabbage WPI 2 | Non-cabbage WPI 1 | Non-cabbage WPI 2 | Experimental threshold^b, **Water (ppt)** | Experimental threshold^b, **WPI (ppb)** | Reported threshold in water (ppb)^c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Dimethyl trisulfide** | 975 | **1.94 ± 0.26** `[M]` | **3.25 ± 0.61** `[M]` | **0.44 ± 0.25** `[M]` | **0.43 ± 0.18** `[M]` | **0.07 ± 1.28** `[M]` | **0.80 ± 0.45** `[M]` | **0.01** `[C]` |

Footnotes, verbatim: *"(a) Determined via solid-phase microextraction (SPME) analysis (mean values
are from 3 reps on each of 2 fibers = total of 6 reps). (b) Best estimate threshold from 80
panelists. (c) Rychlik and others (1998). Concentrations of dimethyl trisulfide are higher in
cabbage whey protein isolate (WPI) compared with non–cabbage WPI (P < 0.05). The sensory threshold
for dimethyl trisulfide in WPI is higher than the threshold in water (P < 0.05)."*

**Note the DB-5 retention index printed here (975) differs from Table 3's (981) for the same
compound in the same paper** (Flags 11). **Note also that the water column is in ppt and the WPI
column in ppb — a 1000x unit change within one table, both spelled out in the column headers and in
the Abstract, so it is deliberate and not a typo.**

### Figure 1 caption (p. C90), typed in full because its content is otherwise unrecoverable

> *"Figure 1—Whey protein isolate (WPI) model system cabbage flavor and aroma intensity and
> similarity to WPI with cabbage flavor. Dimethyl trisulfide (DMTS) was added to a WPI without
> cabbage flavor as described in the text. Different letters denote differences among samples
> (P < 0.05). Cabbage flavor and aroma intensities were scaled using a 15-point universal Spectrum™
> intensity scale in which most dairy flavors fall between 0 and 10. Similarity was scored on a
> 10-point similarity scale: 1 = very different, 10 = identical to reference, in which the reference
> was a WPI that naturally exhibited cabbage flavor."*

**Every value in this figure is figure-only.** The only textual result is: *"The addition of DMTS in
the concentration range naturally found in cabbage flavored WPI to a negative control (WPI without
cabbage flavor) resulted in a cabbage aroma and flavor very similar to the original WPI with cabbage
flavor (Figure 1)."* **No similarity score, no intensity, no p value is printed.**

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| the paper's own summary of the thresholds | *"The orthonasal thresholds for DMTS in water and WPI were 0.07 ± 1.28 parts per trillion (ppt) and 0.80 ± 0.45 ppb, respectively."* | Abstract |
| the paper's own summary of the levels | *"DMTS levels were 1.94 ± 0.26 and 3.25 ± 0.61 parts per billion (ppb) in WPI with cabbage flavor, and 0.44 ± 0.25 and 0.43 ± 0.18 ppb in those without cabbage flavor."* | Abstract |
| **published water threshold, cited** | **0.01 ppb** `[C]` — Rychlik, Schieberle & Grosch 1998, *Compilation of thresholds, odor qualities, and retention indices of key food odorants*, Deutsche Forschungsanstalt für Lebensmittelchemie, Garching | Results p. C89, Table 4 fn (c) |
| the authors' framing of their own disagreement with it | *"Our results indicated that the orthonasal threshold for DMTS in water was much lower than previously reported—0.07 ppt (Table 4), and in WPI the threshold was 0.80 ppb."* | Results p. C90 |
| the authors' explanation of the shift | *"The difference between water and WPI is expected, due to the components of the food matrix (WPI) interacting with flavor compounds and impacting flavor release into the headspace."* | Results p. C90 |
| the pre-defence of the low water value | *"Thresholds can vary widely depending on testing procedure, number of panelists, and matrix used (Meilgaard and others 1999). Previous research in our laboratory has indicated large differences between previously published thresholds and experimental thresholds."* | Results p. C89 |
| the purity caveat, raised by the authors themselves | *"The purity of DMTS used could also impact sensory threshold."* Purity as purchased **98.5 %**, re-measured on receipt by HS-SPME at **95 %** | Results p. C89, Methods p. C88 |
| number of aroma-active compounds found | **22** total; **8 positively identified, 14 tentatively**; all but compounds **9 and 16** previously reported in WPI | Results p. C88 |
| SPME standard-curve R^2 | **0.995** (3-phase fibre) and **0.998** (2-phase fibre) | Results p. C90 |
| SPME calibration range | **500 to 2000 ppt** (= 0.5-2.0 ppb), standards prepared by spiking **deodorised deionised water** | Methods p. C88 |
| the two fibres | *"Our results did show a statistically significant difference between the 2 fibers"* **and** *"Calculated DMTS concentrations from the 2 fibers were not different (P > 0.05) (data not shown)"* | Results p. C90 — **the two sentences contradict each other** (Flags 7) |
| model-system DMTS doses | **1.05, 1.94 and 3.25 ppb** | Methods p. C88 |
| the conclusion, in full | *"DMTS is likely present in some amount in most WPI. However, it is only when the product contains DMTS in concentrations above the threshold level that an off-flavor is apparent."* | Conclusions p. C90 |
| origin of DMDS and DMTS | *"DMDS and DMTS are degradation products of sulfur-containing amino acids."* Both previously found in WPC80 and WPI (Carunchia Whetstine 2005b) **without** documented cabbage flavour | Results p. C89 |
| every model-system intensity and similarity score | **figure-only (Fig. 1)** | — |

### Arithmetic on the printed values (all mine)

**1. THE SHIFT, which the paper never computes.** Putting both legs in ppt: water **0.07**, WPI
**0.80 ppb = 800 ppt**.

> **WPI threshold / water threshold = 800 / 0.07 = 11 428.6x, i.e. ~1.14e4 (mine); log10 shift =
> 4.06 decades (mine).**

**2. The same shift against the cited literature baseline (mine).** Rychlik's water value is
**0.01 ppb = 10 ppt**.

| baseline | WPI / baseline (mine) |
|---|---:|
| **Wright's own measured water leg (0.07 ppt)** | **11 429x** |
| **Rychlik 1998 cited water value (10 ppt)** | **80x** |

**and Wright's water leg is 10 / 0.07 = 143x BELOW the cited value (mine).** The choice of baseline
moves the answer by **143x**. **k2 sec. D.1 records that cross-study matrix/water threshold ratios
span 2000x with a 1-sigma band of 27-41x; the 80x reading sits inside that population and the
11 429x reading sits far outside it.** That is the single most useful sentence for deciding what to
do with this paper.

**3. The `K_g` proxy, computed so nobody has to and so nobody ships it (mine).** `K_g =
(K_water/K_matrix − 1)/protein_g_per_L`, substituting the threshold ratio for the partition ratio:

| protein basis | K_g (mine) | vs the largest shipped constant (`kg_t_2_octenal_pea`, 3.834e-1 L/g) |
|---|---:|---:|
| 100 g/L (WPI solids as stated) | **114.3 L/g** | **298x** |
| 90 g/L (assuming 90 % protein) | **127.0 L/g** | **331x** |
| on the Rychlik baseline, 100 g/L | **0.79 L/g** | **2.1x** |

**Do not ship any of these.** Three reasons: (i) a threshold ratio is not a partition ratio, and
k4b sec. B refutes partition-derived thresholds three independent ways on matched samples; (ii) even
the mildest reading (0.79 L/g on the cited baseline) exceeds every constant in the table, and the
paper's own reading exceeds them by ~300x; (iii) Amendment 6 ruling 2 caps reversible binding at
~25 % of an observed log-shift, and 25 % of 4.06 decades is 1.02 decades = 10.5x — so **even the
capped reversible-binding contribution the layer would allow itself falls 1000x short of explaining
this shift (mine)**. The honest output is an **UNEXPLAINED RESIDUAL of ~1000x on DMTS in a 10 %
whey-protein matrix (mine)**, which is exactly the module's stated design.

**4. OAVs computed both ways, on this paper's own four samples — this is the demonstration
(mine).**

| sample | DMTS (ppb) | panel's cabbage score | OAV on the WPI threshold 0.80 ppb (mine) | OAV on the water threshold 0.00007 ppb (mine) |
|---|---:|---:|---:|---:|
| Cabbage WPI 1 | 1.94 | **3.17** | **2.43** | 27 714 |
| Cabbage WPI 2 | 3.25 | **4.00** | **4.06** | 46 429 |
| Non-cabbage WPI 1 | 0.44 | **ND** | **0.55** | 6286 |
| Non-cabbage WPI 2 | 0.43 | **ND** | **0.54** | 6143 |

**The matrix threshold separates the four samples exactly as the trained panel did: both cabbage
isolates above OAV 1, both clean isolates below it. The water threshold puts all four at OAV
6000-46 000 and cannot distinguish them at all.** This is a **within-paper, four-sample validation
of the matrix-threshold concept against an independent sensory verdict** and it is the strongest
argument in the corpus for the matrix layer existing. The paper states the qualitative version —
*"concentrations of DMTS in WPI with cabbage flavor were above the orthonasal sensory threshold ...
Concentrations of DMTS in WPI without cabbage flavor were below the orthonasal sensory threshold"* —
but **prints none of these eight OAVs**.

**5. Concentration contrasts (mine).** Cabbage vs non-cabbage: **1.94/0.44 = 4.41x**;
**3.25/0.43 = 7.56x**; on the pair means, **2.595/0.435 = 5.97x**. Within the cabbage pair,
**3.25/1.94 = 1.68x**, and the cabbage sensory scores move **4.00/3.17 = 1.26x** in the same
direction — **a 1.68x chemical change producing a 1.26x scale change, which is a compressed
psychophysical response consistent with a power-law exponent near 0.5 (mine, and this paper fits no
such exponent).**

**6. Where the SPME numbers sit relative to their own calibration (mine).** The standards were
**500-2000 ppt = 0.50-2.00 ppb**. Of the four reported concentrations, **only 1.94 ppb is inside the
calibrated range**; 0.44 and 0.43 are **below** it and 3.25 is **above** it. **Three of the four are
extrapolations**, and the two low ones are the samples whose sub-threshold status carries the
paper's central conclusion.

**7. What this says about the repository's DMDS record (mine, and it is indirect).**
`SOURCE_CONTRADICTIONS["dimethyl_disulfide_adduct"]` carries the Anantharamkrishnan Table-2-vs-text
contradiction on DMDS, and `ADDUCT_POSITIVE_CLASSES` already contains `"trisulfide"` while the
`disulfide` class sits on the negative side. **Wright is consistent with that asymmetry but does not
test it**: DMDS appears at log3 FD 1 in the clean WPI and is **not detected** in the off-flavoured
one, is described as **garlic** not cabbage, and is **never quantified** — no concentration, no
threshold, no adduct search. **Wright adds nothing to the DMDS contradiction except a reminder that
DMDS and DMTS behave differently in the same whey matrix.**

## 4. Numbers the repository can use

**Registry mapping.** `dimethyl_trisulfide` is keyed in `data/keys/compounds.yml` (line 688) and has
**no `COMPOUND_STRUCTURE` entry** in `parameters_matrix.py`, though its class name `trisulfide`
already appears in `ADDUCT_POSITIVE_CLASSES`. `dimethyl_disulfide` is keyed (line 675) and does have
a structure record. **`MATRIX_LOADING` has no whey or WPI entry** — the existing keys are `water`,
`skim_milk` (33.9 g/L), `caseinate_1pct` (10 g/L), `gelatin_3pct` (30 g/L), `pea_protein_1pct`
(10 g/L) and `soy_paste_hong`. **A `wpi_10pct` loading at 100 g/L solids would be by far the highest
in the table** and would need the protein-vs-solids distinction recorded, because this paper never
measures protein content.

Every row below shares: **commercial Cheddar-whey WPI, < 3 months old, rehydrated at 10 % solids
(w/v) = 100 g/L in deodorised deionised water, prepared with overhead lights off; pH not stated;
threshold work by ASTM ascending forced-choice method of limits with the Lawless 2000 correction,
ORTHONASAL, DETECTION, 3-AFC, dilution factor 3, five ascending series, n = 80, group threshold =
geometric mean of individual BETs, DMTS delivered from a methanol stock.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **orthonasal DETECTION threshold, dimethyl trisulfide, in DEODORISED WATER** | **0.07 ± 1.28** | **ppt** (parts per trillion; = 7e-5 ug/L) | deodorised deionised water (boiled to 2/3 volume), **80 panelists**, 3-AFC ascending, ASTM E679, methanol carrier | **Table 4, p. C90**; Abstract; Results p. C90 | **threshold `[M]`** — **the ± is undefined and 18x the mean (mine)** (Flags 1) |
| **orthonasal DETECTION threshold, dimethyl trisulfide, in 10 % WPI** | **0.80 ± 0.45** | **ppb** (= 800 ppt = 0.80 ug/L) | 10 % solids WPI (= 100 g/L) in deodorised water, verified cabbage-free by descriptive panel, **80 panelists**, same protocol, **different day** | **Table 4, p. C90**; Abstract; Results p. C90 | **threshold `[M]`** — **the paired matrix leg** |
| **THE PAIRED SHIFT, WPI / water** | **11 429** (log10 = 4.06 decades) | — | same compound, same panel, same protocol, same laboratory | 800/0.07 **(mine)**; **the paper never computes it** | **measured_ratio** — **the corpus's only paired water/matrix threshold on a polysulfide**, and far outside the k2 sec. D.1 cross-study band of 27-41x (1 sigma) |
| **the same shift against the CITED water baseline** | **80** | — | WPI 0.80 ppb over Rychlik 1998's 0.01 ppb | 0.80/0.01 **(mine)** | **derived_assumption** — cross-study; inside the k2 band. **The two readings differ by 143x and the repository must carry both** |
| **published water threshold, DMTS** | **0.01** | ppb (= 10 ppt) | not stated by this paper | Table 4 fn (c), **CITED** from Rychlik, Schieberle & Grosch 1998 | **threshold `[C]`** |
| **Wright's water leg vs the cited water value** | **143x lower** | — | same compound, two laboratories | 0.01/0.00007 **(mine)**; the paper says "much lower than previously reported" | **within_study_ratio** *(strictly cross-study)* — **the single biggest reason not to ship the 11 429x** |
| **`K_g` proxy from the threshold shift** | **114.3** (at 100 g/L) / **127.0** (at 90 g/L) | L/g | — | (ratio − 1)/loading **(mine)** | **derived_assumption — DO NOT SHIP.** ~300x the largest shipped constant; a threshold ratio is not a partition ratio |
| **unexplained residual on DMTS in 10 % WPI** | **~1000x** after the maximum reversible-binding contribution Amendment 6 ruling 2 permits (25 % of a 4.06-decade log-shift = 10.5x) | — | as above | **(mine)** | **measured_bound** — the honest output form, and what the module's docstring prescribes |
| **DMTS concentration, Cabbage WPI 1** | **1.94 ± 0.26** | ppb (ng/g) | 10 % solids, HS-SPME 48 C / 30 min equilibration + 30 min extraction, 6 reps (3 x 2 fibres) | Table 4, p. C90 | **level_only** — **inside the 0.5-2.0 ppb calibration range** |
| **DMTS concentration, Cabbage WPI 2** | **3.25 ± 0.61** | ppb | as above | Table 4 | **level_only** — **ABOVE the calibration range (mine)** |
| **DMTS concentration, Non-cabbage WPI 1** | **0.44 ± 0.25** | ppb | as above | Table 4 | **level_only** — **BELOW the calibration range (mine)**; SD is 57 % of the mean |
| **DMTS concentration, Non-cabbage WPI 2** | **0.43 ± 0.18** | ppb | as above | Table 4 | **level_only** — **BELOW the calibration range (mine)** |
| **OAV in WPI, the four isolates** | **2.43 / 4.06 / 0.55 / 0.54** | — | concentration over the **WPI** threshold | **(mine)**; the paper states only "above"/"below" | **measured_ratio** — **matches the trained panel's cabbage verdict on all four samples** |
| **OAV on the water threshold, the four isolates** | **27 714 / 46 429 / 6286 / 6143** | — | concentration over the **water** threshold | **(mine)** | **derived_assumption** — recorded to show that the water-threshold OAV **cannot** distinguish the four samples and contradicts the panel on two of them |
| cabbage / non-cabbage DMTS concentration contrast | **4.41x**, **7.56x**; **5.97x** on pair means | — | within-study | **(mine)**; the paper says P < 0.05 | **within_study_ratio** |
| descriptive **cabbage** intensity, the four isolates | **ND / ND / 3.17 / 4.00** | 15-point Spectrum scale, **evaluated in the mouth** | 7 trained panelists (>150 h + 40 h whey), duplicate, randomised balanced block, 30 mL at ambient | Table 2, p. C89 | **measured_ratio** — the independent sensory verdict the OAVs are checked against |
| descriptive **overall aroma intensity** (ORTHONASAL), the four isolates | **2.27 / 2.83 / 3.33 / 4.67** | 15-point Spectrum scale | as above | Table 2, p. C89 | **level_only** — the only orthonasal attribute in Table 2 |
| the rest of Table 2 (sweet aromatic, cardboard, brothy, soapy, bitter, astringent) | see Table 2 above | 15-point Spectrum | as above | Table 2, p. C89 | **level_only** |
| **AEDA log3 FD, DMTS** | **2** in non-cabbage WPI 1; **<1** in cabbage WPI 1 | log3 flavour dilution | SAFE extract, non-polar column, 2 sniffers x duplicate | Table 3 row 7, p. C89 | **measured_bound — INVERTED against the quantitation.** The sample with 4.4x less DMTS has the higher FD factor (mine). **Direct evidence that an FD factor is not a concentration** |
| **AEDA log3 FD, dimethyl disulfide** | **1** in non-cabbage WPI 1; **ND** in cabbage WPI 1; odour at port = **"Garlic"** | log3 FD | as above | Table 3 row 3, p. C89 | **level_only** — DMDS is never quantified, never thresholded, and is **absent** from the off-flavoured sample |
| AEDA log3 FD, the other 20 compounds | see Table 3 above, with post-peak intensities and RIs on both columns | log3 FD (post-peak intensity) | as above | Table 3, p. C89 | **level_only** |
| **the paper's own warning about GCO** | *"GCO data alone cannot always be used to identify and confirm key flavor compounds... results are not quantitative, nor do they take into account the role of the food matrix"* | — | — | Results p. C89 | **structural_gate** — quotable justification for the repository's refusal to treat FD factors as quantitative |
| model-system confirmation | DMTS at **1.05, 1.94, 3.25 ppb** added to a clean WPI *"resulted in a cabbage aroma and flavor very similar to the original WPI with cabbage flavor"* | — | methanol carrier (aroma) or 95 % ethanol (flavour); 15-point intensity + 10-point similarity scales | Results p. C90, **Fig. 1** | **measured_bound** — **the outcome is stated but EVERY NUMBER is figure-only** |
| DMTS purity | **98.5 %** as purchased; **95 %** re-measured on receipt by HS-SPME | % | stored −20 C, used within 2 wk | Methods p. C88 | **level_only** — the authors flag it as a possible threshold confound |
| SPME calibration | R^2 **0.995** (3-phase) / **0.998** (2-phase); range **500-2000 ppt**; standards in **water**, not WPI | — | external calibration | Methods p. C88, Results p. C90 | **derived_assumption** — a water-calibrated curve applied to a 10 % protein matrix (Flags 2) |
| protein content of the WPI | — | — | — | **NEVER STATED** | **absent** — only "10 % solids (w/v)" |
| every model-system intensity, similarity and p value | — | — | — | Fig. 1 | **figure_only** |

### Can these be put on the same basis as what the repository already carries?

**(a) The threshold pair: YES as a `threshold` record with both legs, NO as a matrix-correction
factor.** File both values with `route = orthonasal`, `criterion = detection`, `protocol = ASTM
ascending forced choice method of limits + Lawless 2000 correction`, `n_panelists = 80`,
`panel_training = briefly instructed (naive)`, `medium = deodorized_water` and
`medium = wpi_10pct_solids`, `matrix_loading = 100 g/L solids, protein content UNSTATED`,
`carrier = methanol, volume unstated`, and **both** derived shifts (11 429x internal, 80x against
Rychlik) with the 143x baseline disagreement recorded on the record itself. **The lookup table's
`no_measured_threshold` state should NOT be cleared for DMTS silently** — it should become a state
that says "measured, paired, and the two available water baselines disagree by 143x".

**(b) The `K_g`: NO, at any loading, from either baseline.** §3 item 3 shows why: 114 L/g on the
paper's own numbers is ~300x the largest shipped constant, and even the mild 0.79 L/g reading
exceeds every one of the 21 existing rows. The layer's honest output here is the residual, not a
constant.

**(c) The OAV demonstration: YES, and it is the highest-value thing in the paper for the reporting
layer.** Four samples, an independent trained-panel verdict, and two candidate thresholds of which
only one reproduces the verdict. This belongs in the module's justification for the matrix layer
existing at all, alongside k2 sec. D.1.

**(d) The AEDA inversion: YES, as a `structural_gate` against treating FD factors quantitatively.**
DMTS has FD 2 where the concentration is 0.44 ppb and FD <1 where it is 1.94 ppb — a **4.4x
concentration difference read backwards by the dilution assay (mine)**, with the paper's own
explanation printed.

**(e) Nothing here is a rate, a binding constant, an activation energy or a partition coefficient.**
No headspace partition was measured; the SPME is a quantitation, not a partition determination; no
protein-side measurement of any kind was made.

**(f) The temperature question.** The thresholds were sniffed at an unstated ambient temperature
after 1 h out of 5 C refrigeration; the SPME was at 48 C. **Nothing here licenses a DMTS matrix
statement at process temperature**, and a whey isolate is in any case an animal protein — the
repository's plant-protein lane gets no transfer from it without a stated cross-protein assumption.

## 5. Flags

1. **The ± in Table 4 is undefined and the two threshold cells cannot both be arithmetic standard
   deviations.** Water: **0.07 ± 1.28 ppt** — the ± is **18.3x the mean (mine)** and the lower bound
   is −1.21 ppt, which is not a concentration. WPI: **0.80 ± 0.45 ppb** — an ordinary 56 % SD. A
   group BET is a **geometric** mean, whose natural dispersion statistic is a **geometric SD**, a
   multiplicative factor greater than 1: **1.28 works as a geometric SD and 0.45 does not** (it would
   place the geometric mean below its own lower bound). **So the two cells appear to carry different
   statistics under one heading, and the paper defines neither.** Consequence: **no interval can be
   placed on the 11 429x shift**, and the water leg's point value cannot be given an uncertainty at
   all. **This must be resolved before the pair is used quantitatively.**
2. **Three of the four DMTS concentrations lie outside the SPME calibration range (mine).** The
   standards span **500-2000 ppt (0.50-2.00 ppb)**. Measured: 0.44 (below), 0.43 (below), 1.94
   (inside), 3.25 (above). **The two below-range values are exactly the ones whose sub-threshold
   status carries the paper's conclusion.** Worse, **the standards were prepared in deodorised
   deionised WATER, not in WPI**, so the calibration ignores whatever matrix effect on SPME uptake
   the 10 % protein produces — and the paper's own thesis is that that matrix effect is enormous.
   **The reported WPI concentrations are therefore water-equivalent SPME responses, not
   matrix-matched quantitations**, and if protein suppresses headspace DMTS (which the threshold
   result says it does, by four decades), **the true total DMTS in these isolates is higher than
   reported by an unknown factor.** This is the deepest methodological problem in the paper and it is
   never addressed.
3. **The WPI blank is not odourless and the water blank is.** The WPI threshold's two blanks are
   10 % WPI *"determined to be free of cabbage aroma and flavor"* — free of the target note only.
   Table 2 scores clean WPI at overall aroma 2.27-2.83, cardboard 1.67-2.50, sweet aromatic ~2.0.
   **A detection threshold measured against an aromatic blank is elevated by masking and by
   attention/adaptation effects as well as by any physical binding**, and this design cannot separate
   them. **Part of the 11 429x is masking, and the paper attributes all of it to the matrix
   "impacting flavor release into the headspace".**
4. **The 80 threshold panelists were naive; the 7 descriptive panelists were highly trained.**
   *"Subjects were briefly instructed before testing."* Threshold values from a briefly instructed
   80-person panel and flavour intensities from a 7-person panel with >150 h experience are being
   combined into one argument. **They are not the same instrument and the OAV comparison of §3 item 4
   crosses that boundary** (though the fact that it comes out right on all four samples is a point in
   its favour).
5. **The protein content of the WPI is never stated.** "10 % solids (w/v)" is the only composition
   figure in the paper. WPI is conventionally >=90 % protein, but this paper does not say so, does not
   name the suppliers, and does not report ash, lactose or fat. **Every per-gram construction from
   this paper divides by a number the paper never measured**, and the difference between 100 g/L of
   solids and ~90 g/L of protein moves the `K_g` proxy by 11 % (mine) — trivial against the 300x
   problem, but it must be recorded as an assumption, not a measurement.
6. **The methanol carrier volume is never stated.** DMTS stocks were made in methanol and *"aliquots
   of the stock solution were placed into either water or WPI"*, then serially diluted 1:3. The
   authors checked that methanol blanks had no discernible aroma, which addresses odour but not
   solvent competition for hydrophobic protein sites — and methanol at even 1 % changes the
   air/water partition of a small sulfide. **Since the dilution series is 1:3 over five or six steps,
   the methanol concentration also VARIES down the series in step with the DMTS**, which is a
   confound at the low-concentration end where the BET is determined.
7. **The two SPME fibres are described as both different and not different, in adjacent sentences.**
   *"Our results did show a statistically significant difference between the 2 fibers."* immediately
   followed by *"Calculated DMTS concentrations from the 2 fibers were not different (P > 0.05) (data
   not shown)."* The most charitable reading is that the raw responses differed but the
   back-calculated concentrations did not; as printed the two sentences contradict each other, and
   **the supporting data is explicitly withheld ("data not shown")** while the reported means pool
   all six reps across both fibres.
8. **The ASTM standard is cited with two different numbers.** Methods: *"procedure E679-79 (ASTM
   1992)"*. Reference list: *"E-679-91. In: Annual book of standards. 15.07."* with year 1992.
   **E679-79 and E679-91 are different revisions**, and one of the two citations is wrong. The
   procedure as described (ascending forced choice, 3-AFC, BET as geometric mean of the last
   incorrect and first correct concentration) is standard to both.
9. **The threshold is ORTHONASAL and the sensory attribute it is validated against is RETRONASAL.**
   Table 1 states that "Overall aroma intensity" is *"the overall **orthonasal** aroma impact"* and
   that every other term, **cabbage included**, is a *"Flavor (evaluated in the mouth)"*. So the OAV
   check of §3 item 4 compares an orthonasal threshold against an in-mouth intensity. The correlation
   is nonetheless perfect across four samples, and Table 2's orthonasal overall-aroma column orders
   the same way, so the conclusion survives — **but the two measurements are not the same route and
   the paper does not say so.**
10. **The "Brothy 2.33" cell in Table 2 carries no significance letter and appears in only one
    sample.** Its reference compound is methional, which Table 3 shows at log3 FD 2 in the cabbage
    extract and 4 in the clean one — again inverted. Minor, but it means **Cabbage WPI 2 differs from
    the other three in more than DMTS**, and it is the sample carrying the highest DMTS, the highest
    cabbage score and the highest overall aroma intensity. **The two "cabbage" isolates are not
    replicates of one condition.**
11. **Two internal inconsistencies in identifiers.** (i) The DB-5 retention index for DMTS is
    **981 in Table 3** and **975 in Table 4** — same compound, same paper, same column type. (ii) The
    Methods name the GCO columns **Rtx-Wax** and **Rtx-5** (Restek) while Table 3 and its footnotes
    call them **DB-Wax** and **DB-5MS** (Agilent/J&W). Equivalent phases, different manufacturers;
    harmless for chemistry, but it means the RI values cannot be attributed to a specific column
    product.
12. **The AEDA result points the wrong way and the paper says so.** DMTS: log3 FD **2** in the clean
    WPI (0.44 ppb) and **<1** in the cabbage WPI (1.94 ppb). **A 4.4x concentration difference read
    backwards.** The authors handle this correctly and explicitly (*"GCO results were
    inconclusive"*), and the passage is worth quoting in the repository's own documentation — but it
    also means **the whole identification chain in this paper rests on the threshold and
    model-system work, not on the GCO**, and the GCO half of the study contributes no usable number.
13. **The model-system confirmation has no printed numbers.** Its outcome sentence is qualitative
    (*"very similar"*) and **every intensity, similarity score and p value lives in Figure 1**. The
    similarity scale (1-10) and the intensity scale (15-point Spectrum) are defined only in the
    caption. **The paper's third and final line of evidence is therefore unrecoverable from the text
    layer.**
14. **No IRB or ethics approval is reported** for a study with 80 human threshold subjects plus a
    7-member descriptive panel plus 2 sniffers. (Contrast Suppavorasatit 2012, which prints IRB
    Protocol 12038.) Not a scientific flaw, but it is an omission worth noting when the panel data is
    the load-bearing evidence.
15. **No DOI is printed in this article** (see §0). Any citation must go by volume, issue and pages —
    **J Food Sci 71(2):C86-C90 (2006)** — plus manuscript numbers MS 20050644 / FSR 06-05.
16. **This is an ANIMAL protein at 100 g/L.** The repository's active lane is plant protein at
    10-18 g/L. Transferring an 11 429x (or 80x) DMTS shift from 10 % whey isolate to 1 % pea isolate
    crosses **a 5.6-10x loading difference, a protein-family difference, and an unstated pH
    difference** simultaneously. `REVERSIBLE_BINDING`'s own precedent is that a per-gram constant
    scales with loading, but nothing in this paper measures a loading series, so **the loading
    dependence of the threshold shift is entirely untested here.**
17. **What this paper does NOT contain**: any protein content; any pH; any temperature for the
    threshold cups; any loading series; any partition coefficient; any headspace measurement other
    than the SPME quantitation; any binding constant; any rate; any adduct measurement; any
    quantitation of DMDS; any threshold for any compound other than DMTS; any matrix-matched
    calibration; any model-system number.
18. **What to request from the authors**: (i) **what the ± in Table 4 is** — arithmetic SD, geometric
    SD, or standard error — and the individual BETs behind both group thresholds, which would settle
    whether the 0.07 ppt water leg is even a stable estimate; (ii) why the water leg came out 143x
    below Rychlik's, and whether the water series' lowest step was ever below the point at which
    panelists were guessing; (iii) the methanol volume per cup and whether it varied down the
    dilution series; (iv) matrix-matched SPME calibration in 10 % WPI, or at minimum a spike-recovery
    figure, since three of four reported concentrations are extrapolations from a water-based curve;
    (v) the protein content, ash and lactose of the four isolates and their suppliers; (vi) the
    Figure 1 intensity and similarity values as numbers; (vii) whether DMDS was ever quantified;
    (viii) the temperature at which threshold cups were sniffed.
