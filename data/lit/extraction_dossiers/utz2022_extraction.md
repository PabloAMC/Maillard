# Utz 2022 — EXTRACTION (ten commercial pea protein ISOLATES, protein content 42-72 % by Dumas at N x 5.4; 27 odorants above threshold in isolate C quantified by 3-NPH derivatisation + UHPLC-MS/MS stable-isotope dilution; orthonasal aroma profile, complete and minimal recombination, and fifteen 3-AFC omission tests on a 16-member trained panel; nine key food odorants named)

### THE ISOLATE ITSELF, WHERE THE CORPUS HAD ONLY A BEVERAGE: Trikusuma 2020 decodes one pea **beverage** and shows ten of its 21 odorants above threshold before any heat; Utz decodes the **powder** and finds **27 compounds with OAV >= 1 and eight of them above OAV 1000**, with **3-methylbutanal at OAV 10 186 and hexanal at 6202** — but **every one of the 27 odour thresholds is a WATER threshold, and 24 of the 27 are CITED from the Leibniz-LSB@TUM database rather than measured here**, so this paper adds a very large levels table and **not one matrix-corrected threshold**. Two thresholds (2,3-octanedione, (E,E)-3,5-octadien-2-one) are ambiguous between measured-here and cited (§2.6 and the Table 3 footnote contradict each other) and are flagged rather than resolved.

**Source on disk:** `data/articles/Utz2022.pdf` (15 pp., Foods 2022, 11, 412).
Read from the `pdftotext -layout` text layer
(`scratchpad/Utz2022.txt`); **Tables 1, 2, 3 and 4 came through clean** and are re-typed in full
below. **Tables S1-S5 are supplementary and are NOT on disk** — S1 (MRM transitions of the 3-NPH
tagged odorants), S2 (MRM transitions of the pyrazines), S3 (internal standards, calibration curves
and R^2), **S4 (concentrations of every quantified odorant in ALL TEN pea proteins A-J)** and
**S5 (the OAVs behind the heatmap, for all ten proteins)**. **S4 and S5 are the two that matter and
neither is on disk**; everything below on proteins A, B and D-J is therefore restricted to what the
running text says about them. Figures 1 (derivatisation scheme), 2 (matrix calibration curves),
3 (MRM chromatogram), **4 (the aroma profile spider plot for C, its complete recombinant and its
minimal recombinant)**, **5 (the OAV heatmap over all ten proteins)** and 6 (KFO structures) are
images: **every intensity rating in Figure 4 and every cell of the heatmap in Figure 5 is
figure-only.** Repo status before this dossier: Utz 2022 is cited nowhere in
`src/kinetic_core/`, has no extraction dossier, and no benchmark bundle refers to it.

## 0. Identity

| field | value |
|---|---|
| Title | "Sensomics-Assisted Aroma Decoding of Pea Protein Isolates (*Pisum sativum* L.)" |
| Authors | Florian Utz, Andrea Spaccasassi, Johanna Kreissl, Timo D. Stark, Caren Tanger, Ulrich Kulozik, Thomas Hofmann, Corinna Dawid (corresponding, corinna.dawid@tum.de, Tel. +49-81-6171-2901) |
| Affiliations | 1 Chair of Food Chemistry and Molecular Sensory Science, TUM School of Life Sciences, Technical University of Munich, Lise-Meitner-Str. 34, 85354 Freising; 2 Leibniz-Institute for Food Systems Biology at TUM, same street; 3 Chair of Food and Bioprocess Engineering, TUM School of Life Sciences, Weihenstephaner Berg 1, 85354 Freising |
| Venue | **Foods 2022, 11, 412**. Received 10 January 2022; accepted 27 January 2022; published 30 January 2022. Academic Editors: Antonella Verzera and Federico Marini. Open access, CC BY |
| DOI | **`https://doi.org/10.3390/foods11030412`** — printed exactly so in the Citation block on p. 1 and again in the running footer of p. 1 |
| Funding | IGF Project of the FEI, supported via AiF within the Industrial Collective Research programme of the German Ministry of Economic Affairs and Energy (BMWi). **Project AiF 20197 N** |
| Data availability | "All data used in this publication are saved at the Chair for Food Chemistry and Molecular Sensory Science, Freising, Germany." — i.e. **not deposited anywhere public** |
| Conflicts | None declared |
| The decoded sample | **Pea protein C = Nutralys S85F, Roquette (Lestrem, France), batch W317M, protein content 68 %** — "a widely-used pea protein (*Pisum sativum* L.) within Europe". All the sensory work is on C alone |
| Companions on disk | `trikusuma2020_extraction.md` (the pea UHT beverage; **cited by this paper as ref. [11]**), `piornos2025_extraction.md`, `sagesser2023_extraction.md` / `sagesser2024_extraction.md`, `bi2022_extraction.md` (the same group's Bi 2020 pea paper is this paper's ref. [9]) |

## 1. Why it matters

**The gap it fills, stated against what the corpus already has.** `trikusuma2020_extraction.md`
decodes **one pea protein UHT beverage** — a 3 % w/w isolate suspension with carrageenan at pH 7.1 —
and its control column shows **10 of 21 odorants already above threshold before any heat is
applied**. That is the repository's only pea-isolate levels record, and it is a levels record *of a
formulated drink*. **Utz decodes the isolate powder itself**, in ten commercial lots from eight
producers in four countries, with a **16-member panel trained weekly for at least two years**, a
complete aroma recombination that succeeded, and **fifteen 3-AFC omission tests** — which is a
higher standard of sensory evidence than anything else in the pea corpus. It names **nine key food
odorants** and shows by heatmap clustering that they generalise across all ten isolates. If the
model wants to know what a pea protein isolate smells of before any Maillard chemistry starts, this
is the paper.

**The full odorant list with concentrations, thresholds and OAVs is in §3 below, all 27 rows.**

**The threshold question, which is the one that decides whether the repository can use this.**
The task asks, for each threshold, whether it was **measured by this paper** or **cited from a
database**. The paper answers in two places and **the two places disagree**:

> §2.6, "Determination of Odor Thresholds", in full: *"Odor thresholds were taken from the
> Leibniz-LSB@TUM odorant database or determined in water for 2,3-octanedione,
> (E,E)-3,5-octadien-2-one, and (E)-2-dodecenal according to literature [28,29]."*

> Table 3, footnote (c), in full: *"OT = Odor Threshold in water, taken from the Leibniz-LSB@TUM
> odorant database [29]."*

**Reading.** §2.6 says three compounds got their thresholds determined in water by this work
(following the method of Czerny et al. 2008 [28] and/or the database [29]); the footnote to the very
table that carries the numbers says the **whole OT column** came from the database. Two of the three
named compounds do appear in Table 3 — **2,3-octanedione at 29 ug/kg** and **(E,E)-3,5-octadien-2-one
at 27 ug/kg** — and the third, (E)-2-dodecenal, has no Table 3 row (its OAV in protein C is below 1;
§3.4 says it is above 1 only in the cluster-2 proteins). So the honest classification is:

- **24 of the 27 Table 3 thresholds: `[C]` CITED**, from the Leibniz-LSB@TUM odorant database
  (Kreissl, Mall, Steinhaus & Steinhaus, ref. [29], accessed 29 July 2021). Not measured here, and
  the database entry's own panel, criterion and date are not reproduced in this paper.
- **2 of the 27 (2,3-octanedione 29 ug/kg, (E,E)-3,5-octadien-2-one 27 ug/kg): `[C?]` AMBIGUOUS.**
  §2.6 says determined in water for this work; Table 3's footnote says taken from the database.
  **The paper contradicts itself and this dossier does not resolve it.**
- **1 named as determined ((E)-2-dodecenal) does not appear in Table 3 at all** and its value is
  therefore not printed anywhere in the article — it is in Table S5, off disk.

**Consequence for the repository, and it is a hard one.** Every threshold here is **in water**. The
paper says so three times (§2.6, Table 3 footnote (c), Figure 5 caption: "the corresponding odor
thresholds (OT) determined in water"). **There is no threshold in a pea matrix, no paired
water/matrix pair, and no measurement of any threshold in the presence of protein.** The matrix
layer's standing refusal of matrix-corrected thresholds is untouched by this paper. What Utz
supplies is the **numerator** of an OAV at high quality and the **denominator** at second hand.

**Where the OAVs actually live, and a 10x trap.** The concentrations in Table 3 are **per kilogram
of pea protein powder** (40 mg of powder suspended in 1.0 mL of acetonitrile/water). The sensory
work — aroma profile, recombination and every omission test — was done on **10 % pea protein
suspended in water/triacetin (97.5:2.5, w/w)**. So **the concentration the panel actually smelled is
about a tenth of the concentration tabulated (mine)**, and **every OAV in Table 3 is roughly 10x the
OAV of the stimulus the panel judged (mine)**. This does not damage the paper — the recombinant was
built "in native pea protein concentrations" so the ratios are preserved and the reconstitution
succeeded — but it means **an Utz OAV and a Trikusuma OAV are not the same kind of number**
(Trikusuma's are per litre of a 3 % beverage). Any cross-paper OAV comparison must be put on one
basis first (§3, "Arithmetic", item 3).

**What it says about the sulfur lane: essentially nothing, and that absence is itself informative.**
The additional method covered **2-acetyl-2-thiazoline** and **sotolon** among 26 further
sensometabolites, and neither reaches Table 3. **Methional is the only sulfur compound in the whole
27-row list, at 37 ug/kg and OAV 85, and its omission was not individually tested** (it fell in the
non-significant OAV < 100 block, test O2). There is no 2-methyl-3-furanthiol, no 2-furfurylthiol, no
dimethyl disulfide, no dimethyl trisulfide anywhere in this paper. **The repository's sulfur lane
gets no starting level from Utz.**

**And it says something sharp about pyrazines.** §3.1: eight pyrazines were searched for with
pre-tuned MRM transitions and **none was detectable in any pea protein**, even at higher sample
loading, while the identical method found them in coffee and highly roasted cocoa. The stated reason
is that **"the pea proteins investigated were hardly heat-treated during processing."** That is a
clean **measured negative** on the unheated starting state — worth having, because it says the
Maillard nitrogen-heterocycle channel is at zero before the model's cook begins. Note the tension
with the corpus: Trikusuma 2020 quantifies **2,5-dimethylpyrazine at 2.46 ug/L** in its unheated
control beverage (never above threshold), and Bi 2020 (this paper's ref. [9]) reports pyrazines in
*roasted* peas. Utz's negative is a limit-of-quantitation statement, not a zero (Flags 6).

## 2. Methods as they matter to a model

- **The pot — and there are two, which must not be conflated.**
  1. **Analytical pot**: **40 mg of dry protein isolate** suspended in **960 uL acetonitrile/water
     (50:50 v/v)** + **20 uL internal-standard mix**, equilibrated **overnight (>= 20 h) at room
     temperature under continuous shaking**, then derivatised. So the analytical loading is
     **40 mg/mL = 40 g/L of powder** in half-organic solvent. **This is an extraction, not a
     headspace or a partition measurement**: the tabulated concentrations are total content of the
     powder, not headspace-relevant concentrations.
  2. **Sensory pot**: **10 % pea protein** (w/w) suspended in **water/triacetin (97.5:2.5, w/w)**,
     **20 mL or 45 mL in closed sensory vials**, orthonasal, cabin temperature **20-25 C**.
- **Protein loading of the isolates.** Table 1: **42 % to 72 %**, by the **Dumas method** (Vario MAX
  cube, Elementar) with a **nitrogen conversion factor of 5.4**, not 6.25, "as proposed for pea
  proteins" (ref. [23], Mariotti 2008). **Protein C = 68 %.** These are `[M]` and they are the
  cleanest protein contents in the pea corpus — note that neither Bi 2022 nor Sun 2025 measured one
  (Bi's is never stated; Sun's is a supplier claim).
- **pH.** **Never stated.** Not for the analytical suspension, not for the sensory suspension.
- **Temperature and time.** Analytical: >= 20 h at room temperature (equilibration), then
  derivatisation **30 min at 40 C**. Sensory: 20-25 C. Storage of all isolates: **dark, 4 C**.
  **There is no heating step anywhere in this paper.** It is a decode of an unprocessed starting
  material.
- **Derivatisation chemistry (this is what makes the method work and what limits it).**
  **3-nitrophenylhydrazine hydrochloride (3-NPH)**, 20 uL of 200 mmol/L in acetonitrile/water, with
  **EDC** (N-(3-(dimethylamino)propyl)-N'-ethylcarbodiimide hydrochloride), 20 uL of 120 mmol/L in
  acetonitrile/water containing **6 % pyridine**, 30 min at 40 C. 3-NPH forms hydrazones with
  **carbonyls** and, with EDC activation, amides with **carboxylic acids**. **This is why the list is
  27 aldehydes, ketones and acids and nothing else**: a compound with no carbonyl and no carboxyl
  is invisible to this method by construction (Flags 5). Filtered through 0.45 um (Minisart RC 15),
  1 uL injected.
- **The measurement family: UHPLC-MS/MS with stable-isotope dilution assay (SIDA).** This is
  **NOT** GC-olfactometry and **NOT** headspace GC. Exion LC UHPLC + **QTRAP 6500+** (AB Sciex),
  **ESI+ at +5500 V**, Kinetex 1.7 um XB-C18, 100 x 2.1 mm, 100 Å; gradient of 0.1 % formic acid in
  water / in acetonitrile at 0.4 mL/min: 0 min 27 % B, 0.5 min 27 %, 1 min 50 %, 6 min 100 %,
  7 min 100 %, 7.5 min 27 %, 9 min 27 %. Nebuliser 55 psi, turbo gas 450 C, drying 65 psi, curtain
  35 psi, collision gas 1.5e-5 torr, unit resolution, scheduled MRM windows of **+/-30 s**. Analyst
  1.6.3 / MultiQuant.
- **Eight internal standards, all isotope-labelled** (in acetonitrile/water 50:50): 3-methylbutanal-
  d2 (22.0 ug/mL), hexanal-d12 (1.3 ug/mL), decanal-d2 (23.3 ug/mL), diacetyl-d6 (19.9 ug/mL),
  hexanoic acid-d3 (5.5 ug/mL), phenylacetic acid-13C2 (2.7 ug/mL), vanillin-d3 (1.8 ug/mL),
  gamma-nonalactone-d2 (15.5 ug/mL). Decanal-d2, vanillin-d3 and gamma-nonalactone-d2 were
  **synthesised at the Leibniz-LSB@TUM**. **Each of the 27 analytes is quantified against one of
  these eight**, and Table 3 names which — several are structurally distant from their standard
  (Flags 3).
- **Matrix-effect validation, which is unusually well done.** No analyte-free pea matrix exists, so
  calibration curves were run **with** 40 mg/mL of protein C present (standard addition) and
  **without** it. The curves showed "either the same slope, just shifted by a certain amount for the
  analytes present in pea protein, e.g., for hexanal and hexanoic acid, or congruent curves for no
  or low abundance, such as for (E,Z)-2,6-nonadienal and 2,3-octanedione (Figure 2)". Conclusion, as
  printed: **"matrix effects during ionization were fully compensated by the selected internal
  standards."** Recoveries then measured in triplicate in solvent: **80.5 % to 106.9 %** (Table 2).
- **LOD/LOQ.** Signal-to-noise 3 and 10 respectively. **LOD <0.1 to 5.8 nmol/L; LOQ <0.1 to
  19.2 nmol/L** (Table 2). The paper claims "All LOQs showed higher sensitivity than specific aroma
  thresholds or could be counterbalanced by increasing sample loading."
- **Sensory panel.** **Sixteen panelists (nine women, seven men, age 22-58)** from the Chair of Food
  Chemistry and Molecular Sensory Science and the Leibniz-LSB@TUM, **each trained weekly for a
  minimum of two years**, no history of known anosmia. **Orthonasal throughout** — the paper says so
  explicitly in §2.5 and again in Eq. (1)'s denominator ("orthonasal odor threshold"). QDA training
  used aqueous reference solutions at **10-fold odor thresholds**: hexanal 25.0 ug/L (grassy),
  3-isopropyl-2-methoxypyrazine 0.1 ug/L (beans-like), acetic acid 60 mg/L (sour),
  (E,E)-2,4-decadienal 0.32 ug/L (fatty), phenylacetic acid 0.68 mg/L (honey-like),
  2-ethyl-5-methylpyrazine 1.0 mg/L (nutty), 3-methylbutanal 5.0 ug/L (malty),
  2,3,5-trimethylpyrazine 0.12 mg/L (earthy). **Eight attributes rated 0 (not detectable) to 5 (very
  intense).**
- **The deodorised base for recombination.** Protein C (500 g) stirred in freshly distilled
  **n-pentane (2 x 1.5 L)** overnight, then **dichloromethane (2 x 1.5 L)** overnight at room
  temperature, dried under nitrogen. "The obtained pea protein powder **could not be sensorially
  related to pea by the panelists**." Recombinants: 10 % deodorised powder in water + the
  recombination solution in triacetin, 97.5:2.5 w/w, at **native pea protein concentrations**.
- **Omission tests.** **3-alternative forced choice (3-AFC)**, incomplete recombinant against
  complete recombinant, **13-14 panelists per test**, p by binomial distribution (ref. [27]).
  Fifteen tests: **O1 and O2 are BLOCKS** (all OAV < 10, and all OAV < 100 respectively);
  **O3-O15 are the thirteen individual odorants with OAV > 100**.
- **OAV definition, Eq. (1).** `OAV = concentration / orthonasal odor threshold`. **The concentration
  is the Table 3 per-kilogram-of-powder value and the threshold is a water threshold** — see §1's
  10x trap and Flags 1.
- **Statistics.** Concentrations are the mean of **n = 3 independent sample workups**, +/- SD, with
  RSD printed. Heatmap in R 4.0.4 with `ComplexHeatmap`, **log-transformed OAV**.

## 3. Tables re-typed

### Table 1 (p. 3). "List of analyzed pea proteins within FEI project AiF 20197 N."

| Code | Commercial Name | PC (a) | Batch | Producer |
|---|---|---:|---|---|
| A | Nutralys F85F | 67 % `[M]` | W084M | Roquette, Lestrem, France |
| B | Prestige | 42 % `[M]` | KO67X | Parrheim Foods, Saskatoon, Canada |
| **C** | **Nutralys S85F** | **68 %** `[M]` | **W317M** | **Roquette, Lestrem, France** |
| D | Bio Erbsen Protein | 70 % `[M]` | 15720-30118 | Golden Peanut, Garstedt, Germany |
| E | Bio Erbsen Protein | 68 % `[M]` | 170180323 | Piowald, Mühbrook, Germany |
| F | 1501018 | 68 % `[M]` | 83526272 | Döhler, Darmstadt, Germany |
| G | Pea Pro | 68 % `[M]` | 06102016 | LSP Sports Nutrition, Bonn, Germany |
| H | Empro E86HV | 72 % `[M]` | 41266 | Emsland-Stärke, Emlichheim, Germany |
| I | Empro E86 | 71 % `[M]` | 41266 | Emsland-Stärke, Emlichheim, Germany |
| J | Pisane C9 | 69 % `[M]` | 817021 | Cosucra Group, Warcoing, Belgium |

Footnote (a): "PC = Protein Content, determined using the Dumas method and a conversion factor of
5.4 (see Section 2.3)." **Note H and I share batch number 41266 while being sold as different
grades** (Flags 9).

### Table 2 (p. 8). "Performed validation experiments of important odorants in pea protein isolates."

| Validated Analyte | Add. (a) (umol/L) | Add. Found ± SD (b) (umol/L) | RSD (c) (%) | Recovery (%) | LOD (d) (nmol/L) | LOQ (e) (nmol/L) |
|---|---:|---:|---:|---:|---:|---:|
| 2-/3-methylbutanal, sum | 0.668 `[M]` | 0.667 ± 0.038 `[M]` | 5.7 | 99.8 | <0.1 | <0.1 |
| 2-methylbutanal | 0.346 `[M]` | 0.344 ± 0.009 `[M]` | 2.7 | 99.5 | <0.1 | <0.1 |
| hexanal | 0.333 `[M]` | 0.277 ± 0.004 `[M]` | 1.5 | **83.0** | 2.7 | 9.0 |
| heptanal | 0.301 `[M]` | 0.301 ± 0.006 `[M]` | 1.8 | 99.8 | <0.1 | <0.1 |
| methional | 0.370 `[M]` | 0.336 ± 0.008 `[M]` | 2.4 | 90.8 | 0.2 | 0.5 |
| (E)-2-octenal | 0.292 `[M]` | 0.259 ± 0.003 `[M]` | 1.0 | 88.8 | 0.5 | 1.7 |
| (E,E)-2,4-nonadienal | 0.274 `[M]` | 0.272 ± 0.004 `[M]` | 1.5 | 99.3 | 0.2 | 0.5 |
| (E,Z)-2,6-nonadienal | 0.294 `[M]` | 0.294 ± 0.008 `[M]` | 2.8 | 100.0 | 0.2 | 0.8 |
| (E,E)-2,4-decadienal | 0.246 `[M]` | 0.241 ± 0.003 `[M]` | 1.3 | 98.0 | 0.2 | 0.5 |
| (E)-2-undecenal | 0.278 `[M]` | 0.267 ± 0.020 `[M]` | 7.4 | 96.2 | 0.3 | 0.9 |
| (E)-2-dodecenal | 0.275 `[M]` | 0.258 ± 0.015 `[M]` | 5.7 | 93.7 | 1.4 | 4.7 |
| 2,3-octanedione | 0.293 `[M]` | 0.290 ± 0.003 `[M]` | 0.9 | 99.0 | 0.4 | 1.2 |
| (E,E)-3,5-octadien-2-one | 0.256 `[M]` | 0.264 ± 0.006 `[M]` | 2.1 | 103.2 | <0.1 | 0.3 |
| 2-undecanone | 0.291 `[M]` | 0.288 ± 0.011 `[M]` | 3.8 | 99.1 | <0.1 | <0.1 |
| hexanoic acid | 0.305 `[M]` | 0.288 ± 0.012 `[M]` | 4.3 | 94.4 | <0.1 | <0.1 |
| heptanoic acid | 0.280 `[M]` | 0.285 ± 0.013 `[M]` | 4.7 | 101.7 | **5.8** | **19.2** |
| phenylacetaldehyde | 0.286 `[M]` | 0.306 ± 0.009 `[M]` | 2.9 | **106.9** | <0.1 | 0.3 |
| 4-ethyl benzaldehyde | 0.400 `[M]` | 0.378 ± 0.023 `[M]` | 6.1 | 94.5 | 4.7 | 15.8 |
| vanillin | 0.285 `[M]` | 0.258 ± 0.012 `[M]` | 4.7 | 90.6 | 0.2 | 0.6 |
| γ-octalactone | 0.302 `[M]` | 0.243 ± 0.010 `[M]` | 4.0 | **80.5** | 1.0 | 3.2 |

Footnotes: "(a) Add. = Addition. (b) Add. found = Addition found, SD = Standard Deviation,
determined based on replicate sample workup and analysis (n = 3). (c) RSD = Relative Standard
Deviation. (d) LOD = Limit of Detection, determined based on a signal-to-noise ratio of 3.
(e) LOQ = Limit of Quantitation, determined based on a signal-to-noise ratio of 10."

**Note what is NOT validated here**: acetaldehyde, diacetyl, acetoin, benzaldehyde, acetic acid,
3-methylbutanoic acid, 2-methylbutanoic acid, nonanoic acid, octanoic acid, decanoic acid and
phenylacetic acid — **eleven of the 27 Table 3 rows**, all of them the ones quantified by the
"additional" dairy method (footnote (f) of Table 3). **No recovery, LOD or LOQ is printed for any of
them in this paper** (Flags 4).

### Table 3 (p. 9). "Concentrations of aroma-active compounds in pea protein C in descending OAV order."

| No. | Aroma-Active Analyte | Used IS | Odor Quality | Mean ± SD (a) (ug/kg) | RSD (b) (%) | OT (c) (ug/kg) | OAV (d) |
|---:|---|---|---|---:|---:|---:|---:|
| 1 | 3-methylbutanal (e) | 3-methylbutanal-d2 | malty | **5093** `[M]` (5360 ± 547) (e) | (10.2) | **0.5** `[C]` | **10186** |
| 2 | hexanal | hexanal-d12 | green, grassy | **14,886 ± 1904** `[M]` | 12.8 | **2.4** `[C]` | **6202** |
| 3 | acetaldehyde (f) | acetaldehyde-d3 | fresh, green | **72,197 ± 4231** `[M]` | 5.9 | **16** `[C]` | **4512** |
| 4 | (E,E)-2,4-decadienal | decanal-d2 | fatty, deep-fried | **101 ± 10** `[M]` | 9.8 | **0.027** `[C]` | **3736** |
| 5 | phenylacetaldehyde | phenylacetic acid-13C2 | flowery, honey-like | **6097 ± 303** `[M]` | 5.0 | **5.2** `[C]` | **1173** |
| 6 | (E,E)-2,4-nonadienal | decanal-d2 | fatty, green | **53 ± 7** `[M]` | 12.5 | **0.046** `[C]` | **1156** |
| 7 | (E)-2-octenal | hexanal-d12 | fatty, nutty | **907 ± 8** `[M]` | 0.8 | **1.7** `[C]` | **533** |
| 8 | diacetyl (f) | diacetyl-d6 | butter-like | **316 ± 25** `[M]` | 8.0 | **0.96** `[C]` | **329** |
| 9 | benzaldehyde (f) | phenylacetic acid-13C2 | bitter almond-like, marzipan-like | **37,201 ± 7850** `[M]` | 21.1 | **150** `[C]` | **248** |
| 10 | heptanal | hexanal-d12 | citrus-like, fatty | **1326 ± 130** `[M]` | 9.8 | **6.1** `[C]` | **217** |
| 11 | 2-methylbutanal | 3-methylbutanal-d2 | malty | **267 ± 8** `[M]` | 2.9 | **1.5** `[C]` | **178** |
| 12 | (E)-2-undecenal | decanal-d2 | soapy, metallic | **123 ± 8** `[M]` | 6.5 | **0.78** `[C]` | **157** |
| 13 | nonanoic acid (f) | octanoic acid-d15 | moldy, pungent | **2776 ± 167** `[M]` | 6.0 | **26** `[C]` | **107** |
| 14 | methional | hexanal-d12 | cooked potato-like | **37 ± 5** `[M]` | 4.7 | **0.43** `[C]` | **85** |
| 15 | acetic acid (f) | acetic acid-13C2 | vinegar-like | **262,037 ± 10,440** `[M]` | 4.0 | **5600** `[C]` | **47** |
| 16 | 3-methylbutanoic acid (f) | butyric acid-13C4 | sweaty | **23,274 ± 1425** `[M]` | 6.1 | **490** `[C]` | **47** |
| 17 | decanoic acid (f) | octanoic acid-d15 | soapy, musty | **51 ± 3** `[M]` | 5.5 | **3.5** `[C]` | **15** |
| 18 | vanillin | vanillin-d3 | vanilla-like, sweet | **561 ± 18** `[M]` | 3.2 | **53** `[C]` | **11** |
| 19 | (E,E)-3,5-octadien-2-one | diacetyl-d6 | woody, mushroom-like, green | **231 ± 18** `[M]` | 7.9 | **27** `[C?]` | **9** |
| 20 | hexanoic acid | hexanoic acid-d3 | sweaty | **36,802 ± 1978** `[M]` | 5.4 | **4800** `[C]` | **8** |
| 21 | octanoic acid (f) | octanoic acid-d15 | carrot-like, musty | **1536 ± 48** `[M]` | 3.1 | **190** `[C]` | **8** |
| 22 | phenylacetic acid (f) | phenylacetic acid-13C2 | honey-like, beeswax-like | **516 ± 61** `[M]` | 11.8 | **68** `[C]` | **8** |
| 23 | γ-octalactone | γ-nonalactone-d2 | coconut-like | **47 ± 5** `[M]` | 10.5 | **6.5** `[C]` | **7** |
| 24 | 2-methylbutanoic acid (f) | butyric acid-13C4 | malty, fruity, sweaty | **21,419 ± 1014** `[M]` | 4.7 | **3100** `[C]` | **7** |
| 25 | 2,3-octanedione | diacetyl-d6 | mushroom-like, dill-like, broccoli-like | **141 ± 13** `[M]` | 8.9 | **29** `[C?]` | **4.8** |
| 26 | 2-undecanone | decanal-d2 | soapy, green | **68 ± 4** `[M]` | 6.5 | **24** `[C]` | **2.8** |
| 27 | acetoin (f) | diacetyl-d6 | butter-like, carrot-like | **978 ± 7** `[M]` | 0.7 | **590** `[C]` | **1.7** |

Footnotes, verbatim: "(a) Mean = arithmetic mean, SD = Standard Deviation, determined based on
replicate sample workup and analysis (n = 3). (b) RSD = Relative Standard Deviation. (c) OT = Odor
Threshold in water, taken from the Leibniz-LSB@TUM odorant database [29]. (d) Concentration divided
by the odor threshold and expressed as OAV (Odor Activity Value). (e) The concentration of
3-methylbutanal was determined by subtraction the sum value (in brackets) from the 2-methylbutanal
concentration. (f) Quantified by an additional 3-NPH-UHPLC-MS/MS method published for dairy analysis
[16]."

**`[C?]` marks the two rows the paper contradicts itself about** (see §1): §2.6 says 2,3-octanedione
and (E,E)-3,5-octadien-2-one thresholds were "determined in water" for this work, footnote (c) says
the whole OT column came from the database. Both readings are printed; neither is retracted.

**Also NOT identified in ANY of the ten proteins** (§3.2, from Table S4, off disk):
**(E,Z)-2,6-nonadienal** and **4-ethylbenzaldehyde**. These are `[M]` measured negatives at the LOQs
of Table 2 (0.8 and 15.8 nmol/L).

### Table 4 (p. 11). "Omission experiments applied to the aroma model of pea protein C."

| Test | Odorant(s) Omitted (a) | OAV | p Value (%) | Significance (b) |
|---|---|---:|---:|---|
| O1 | 19-27 *(block)* | <10 | 23.8 | NS |
| O2 | 14-18 *(block)* | <100 | 19.1 | NS |
| **O3** | **1, 3-methylbutanal** | **10186** | **<0.1** | **\*\*\*** |
| **O4** | **2, hexanal** | **6202** | **2.6** | **\*** |
| **O5** | **3, acetaldehyde** | **4512** | **4.8** | **\*** |
| O6 | 4, (E,E)-2,4-decadienal | 3736 | 9.2 | NS |
| O7 | 5, phenylacetaldehyde | 1173 | 21.4 | NS |
| **O8** | **6, (E,E)-2,4-nonadienal** | **1156** | **4.8** | **\*** |
| **O9** | **7, (E)-2-octenal** | **533** | **1.3** | **\*** |
| O10 | 8, diacetyl | 329 | 15.6 | NS |
| **O11** | **9, benzaldehyde** | **248** | **4.0** | **\*** |
| **O12** | **10, heptanal** | **217** | **<0.1** | **\*\*\*** |
| **O13** | **11, 2-methylbutanal** | **178** | **0.7** | **\*\*** |
| O14 | 12, (E)-2-undecenal | 157 | 23.0 | NS |
| **O15** | **13, nonanoic acid** | **107** | **4.0** | **\*** |

Footnotes: "(a) Odorant numbers refer to Table 3. (b) NS, no significance (p > 5%); \*, significance
(5% >= p > 1%); \*\*, highly significance (1% >= p > 0.1%); \*\*\*, very highly significance
(p <= 0.1%)." **13-14 panelists per 3-AFC test.** Note the footnote's own inconsistency: the
double-asterisk band is defined as `1% >= p > 0.1%` in the footnote and described in the text as
`1% >= p < 0.1%`; the printed value for O13 (0.7 %) falls in the footnote's band either way.

**The nine key food odorants, as named in the Abstract and §4**, in descending OAV order:
**3-methylbutanal, hexanal, acetaldehyde, (E,E)-2,4-nonadienal, (E)-2-octenal, benzaldehyde,
heptanal, 2-methylbutanal, nonanoic acid.** **Eight of the nine are aldehydes**; the ninth is an
acid. **Not one is a Maillard nitrogen or sulfur heterocycle.**

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| number of analytes with OAV >= 1 in protein C | **27** | §3.2, p. 8 |
| the six with OAV > 1000, as quoted in the text | 3-methylbutanal (5.1 mg/kg; **OAV 10186**), hexanal (14.9 mg/kg; **OAV 6202**), acetaldehyde (72.2 mg/kg; **OAV 4512**), (E,E)-2,4-decadienal (101 ug/kg; **OAV 3741**), phenylacetaldehyde (6.1 mg/kg; **OAV 1173**), (E,E)-2,4-nonadienal (53 ug/kg; **OAV 1157**) | §3.2, p. 8 — **note 3741 and 1157 against Table 3's 3736 and 1156** (Flags 8) |
| recovery range across all validated analytes | **80.5 % to 106.9 %** | §3.1, p. 7 |
| LOD range / LOQ range | **<0.1 to 5.8 nmol/L / <0.1 to 19.2 nmol/L** | §3.1, p. 7 |
| **pyrazines** — 2,3-dimethyl-, 2,5-dimethyl-, 2,6-dimethyl-, 2,3,5-trimethyl-, 2-ethyl-, 2-ethyl-5(6)-methyl-, 3-isopropyl-2-methoxy-(5/6)-methyl-, 2-isobutyl-3-methoxypyrazine | **"could not be detected in pea protein, even when higher sample amounts were used"**; the same method found them in coffee and highly roasted cocoa; assumed below LOQ; reason given: **"the pea proteins investigated were hardly heat-treated during processing"** | §3.1, p. 7 |
| **the corpus's only prior pea-protein quantitation, as this paper reports it** | Murat 2013 [10], by GC-MS with hexanal-d12 as a relative internal standard: **hexanal 83 mg/kg**, **phenylacetaldehyde 2.1 mg/kg**; (E,E)-2,4-decadienal could not be quantitated (coeluting peaks) | §3.2, p. 9 `[C]` |
| Murat 2013 values that are **higher** than Utz's | heptanal **16.2 mg/kg**, nonanoic acid **4.3 mg/kg**, methional **97 ug/kg**, 2-undecanone **557 ug/kg**, 2,3-octanedione **4.1 mg/kg**, (E,E)-3,5-octadien-2-one **24.4 mg/kg** | §3.2, p. 10 `[C]` |
| Murat 2013 values that are **lower** than Utz's | (E)-2-octenal **312 ug/kg**, benzaldehyde **6.4 mg/kg**, vanillin **268 ug/kg** | §3.2, p. 10 `[C]` |
| acids Murat could not determine | 2- and 3-methylbutanoic, pentanoic, hexanoic, heptanoic, octanoic, dodecanoic | §3.2, p. 10 |
| the complete recombinant | 27 analytes at native concentrations in deodorised protein C; **"the results proved the quantified aroma compounds in their correct ratio and the successful aroma reconstitution"** | §3.3, p. 10, Fig. 4 |
| the minimal recombinant | the **nine KFO** only; **"highlighted very high similarity compared to the complete recombinant"** | §3.3, p. 11, Fig. 4 |
| **heatmap clusters over the ten isolates** | **cluster 1: A, B, J, C, G**; **cluster 2: D, E, F, H, I** | §3.4, p. 12, Fig. 5 |
| what separates the clusters | **(E)-2-dodecenal was highly present in cluster 2** and did not exceed OAV 1 in cluster 1; "the individual contribution of KFO aldehydes slightly differed among the examined pea proteins" | §3.4, p. 12 |
| which KFO clustered with the low-threshold aldehydes | **eight of the nine**: 3-methylbutanal, hexanal, acetaldehyde, (E,E)-2,4-nonadienal, 2-methylbutanal, benzaldehyde, (E)-2-octenal, heptanal — i.e. **all but nonanoic acid** | §3.4, p. 12 |
| **all ten proteins' concentrations and OAVs** | **Tables S4 and S5 — NOT on disk** | §3.4 |
| every intensity in the aroma profile | **figure-only (Fig. 4)** — eight attributes x three samples, no number printed | §3.3 |
| every heatmap cell | **figure-only (Fig. 5)**, log-transformed | §3.4 |

### Arithmetic on the printed values (all mine)

**1. Every OAV in Table 3 reproduces from its own two columns to within rounding (mine).** Spot
checks: 5093/0.5 = 10186 ✓; 14 886/2.4 = 6202.5 ✓; 72 197/16 = 4512.3 ✓; 6097/5.2 = 1172.5 ✓;
907/1.7 = 533.5 ✓; 316/0.96 = 329.2 ✓; 37 201/150 = 248.0 ✓; 1326/6.1 = 217.4 ✓; 267/1.5 = 178.0 ✓;
2776/26 = 106.8 ✓; 262 037/5600 = 46.8 ✓; 23 274/490 = 47.5 ✓; 561/53 = 10.6 ✓; 36 802/4800 = 7.67 ✓;
141/29 = 4.86 ✓; 978/590 = 1.66 ✓. **Two do not:** (E,E)-2,4-decadienal 101/0.027 = **3740.7**
against a tabulated **3736** (the text's 3741 is the correct one); (E,E)-2,4-nonadienal
53/0.046 = **1152.2** against a tabulated **1156** (the text's 1157 needs a concentration of 53.2).
**Both discrepancies are display rounding of the concentration column, not errors in the numbers**,
and both are under 0.4 % (Flags 8).

**2. The 3-methylbutanal row is a difference of two measurements, and its error is not propagated.**
Footnote (e): the sum of 2- and 3-methylbutanal was measured (5360 ± 547, RSD 10.2 %) and
2-methylbutanal separately (267 ± 8). **5360 − 267 = 5093 ✓ (mine).** But the printed SD column for
row 1 gives the SD **of the sum** in brackets, and no SD is printed for the difference. Propagating
in quadrature: sqrt(547^2 + 8^2) = **547.1, i.e. ±10.7 % on 5093 (mine)**. **The paper's single
highest-OAV compound therefore carries a ±11 % uncertainty that the table does not display.**

**3. Putting Utz on the same basis as Trikusuma 2020 (mine, and this is the comparison the corpus
needs).** Trikusuma measures a 3 % w/w beverage; `trikusuma2020_extraction.md` §3 item 4 already
converts its control column to a per-kilogram-of-isolate basis under a stated assumption. Utz's
Table 3 is already per kilogram of powder. **The two are directly comparable on that basis:**

| compound | Utz 2022, protein C (ug/kg powder) | Trikusuma 2020, control, per kg isolate (derived) | Utz / Trikusuma (mine) |
|---|---:|---:|---:|
| hexanal | 14 886 | 11 000 | **1.35x** |
| heptanal | 1326 | 232 | **5.7x** |
| methional | 37 | 18.3 | **2.0x** |
| (E,E)-2,4-nonadienal | 53 | 23.7 | **2.2x** |
| (E)-2-octenal | 907 | 10.3 | **88x** |
| (E,E)-2,4-decadienal | 101 | 2.0 | **51x** |

**Hexanal agrees to 1.35x across two laboratories, two continents, two commercial isolates and two
completely different analytical methods (LC-MS/MS on a solvent extract vs dynamic-headspace GC/MS on
a beverage). That is the strongest cross-validation of a pea-isolate level anywhere in the corpus.**
The alkenals do not agree at all — 88x and 51x — which is exactly the direction a headspace method
would err against an extraction method for reactive alpha,beta-unsaturated aldehydes that are
partly bound or adducted in the beverage (and Trikusuma's numbers pass through a 40-minute hydration
plus a carrageenan matrix that Utz's solvent extraction bypasses). **Neither ladder should be
transferred to the other; the hexanal agreement should be cited as a levels anchor.**

**4. The 10 % suspension correction (mine).** The panel smelled **10 % powder** in water/triacetin.
Taking the suspension density as 1.00 kg/L, a Table 3 value of C ug/kg of powder corresponds to
**0.1 x C ug/L in the cup, if every molecule partitioned out of the powder and none was retained by
the protein**. On that reading the in-cup OAVs are **one tenth of Table 3's**: 3-methylbutanal
~1019, hexanal ~620, acetaldehyde ~451, and the OAV = 1 line moves up to a **Table 3 OAV of 10**,
which is — strikingly — **exactly where the paper's own omission block O1 stops mattering** (all
OAV < 10 omitted, p = 23.8 %, not significant). **That coincidence is worth recording**: the
sensory result is consistent with the in-cup OAV threshold rather than the tabulated one. It is not
proof — protein binding pushes the true in-cup free concentration lower still, and triacetin is a
2.5 % co-solvent — but it means **the tabulated OAVs are an upper bound on what the panel
experienced, by about a factor of ten (mine)**.

**5. Cross-check on the recovery of hexanal (mine).** Table 2 gives hexanal recovery **83.0 %**, the
second-lowest in the table after gamma-octalactone's 80.5 %. If Table 3's 14 886 ug/kg is not
recovery-corrected (the paper does not say it is), the true content would be **14 886/0.830 =
17 935 ug/kg (mine)** and the OAV **7473 (mine)**. **The paper does not state whether the reported
concentrations are recovery-corrected.** Since the quantitation is SIDA against hexanal-d12 added
*before* workup, the labelled standard should already correct for extraction losses and the
recovery figure is a process check rather than a correction factor — but that is my inference, not
the paper's statement (Flags 7).

**6. What the omission tests actually establish about OAV as a ranking device.** Reading Table 4:
of the thirteen individually tested odorants, **nine were significant and four were not**, and the
four non-significant ones include **the #4 and #5 highest OAVs in the entire paper**
((E,E)-2,4-decadienal at 3736 and phenylacetaldehyde at 1173) while **the #10 OAV, heptanal at 217,
gave the joint-strongest result in the study (p < 0.1 %)**. **OAV rank and sensory necessity are
badly decorrelated here**: the Spearman relationship between OAV rank and significance is visibly
weak, and the paper does not remark on it. **This is a direct, same-panel, same-day measurement of
how far an OAV ranking can be trusted, and it says: not far.** For a repository whose output layer
is "ratios and rankings first, OAVs second", this is the most useful methodological finding in the
paper.

## 4. Numbers the repository can use

**Registry mapping.** Of the 27 compounds, `hexanal`, `3_methylbutanal`, `2_methylbutanal`,
`benzaldehyde`, `t_2_octenal` (= (E)-2-octenal), `methional`, `diacetyl`, `vanillin`, `furaneol`-
adjacent none, `butyric_acid` are the kinds of key already present in `data/keys/compounds.yml` /
`COMPOUND_STRUCTURE`; **acetaldehyde, (E,E)-2,4-decadienal, (E,E)-2,4-nonadienal, (E)-2-undecenal,
(E,E)-3,5-octadien-2-one, 2,3-octanedione, 2-undecanone, acetoin, gamma-octalactone,
phenylacetaldehyde, phenylacetic acid, nonanoic, octanoic, decanoic, hexanoic, heptanoic,
2-methylbutanoic and 3-methylbutanoic acids** should be checked individually before use. **`maltol`
is keyed in this repository and does NOT appear in Utz at all** — worth noting against Trikusuma,
which reports 48 000 ug/kg of maltol on a per-isolate basis in its control beverage. A compound that
large in one pea isolate and absent from a 27-row decode of another is a discrepancy the levels
layer should carry (Flags 10).

Every row below shares: **commercial pea protein isolate C = Nutralys S85F (Roquette), batch W317M,
68 % protein by Dumas at N x 5.4, stored dark at 4 C; 40 mg powder in 960 uL acetonitrile/water
50:50 + 20 uL IS, >= 20 h at room temperature, 3-NPH/EDC derivatisation 30 min at 40 C,
UHPLC-MS/MS SIDA, n = 3.** **Concentrations are per kilogram of POWDER, not per litre of any
beverage.** **Thresholds are ORTHONASAL, in WATER, and are cited unless marked otherwise.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **acetaldehyde** in pea isolate C | **72 197 ± 4231** | ug/kg powder | as above; additional dairy method | Table 3 row 3, p. 9 | **level_only** — the single largest odorant by mass in the paper |
| **hexanal** in pea isolate C | **14 886 ± 1904** | ug/kg powder | as above | Table 3 row 2 | **level_only** — and the corpus's best-corroborated pea level (1.35x against Trikusuma, mine) |
| **3-methylbutanal** in pea isolate C | **5093** (from a sum of 5360 ± 547 minus 267 ± 8) | ug/kg powder | as above | Table 3 row 1 + footnote (e) | **level_only** — **a difference of two measurements; propagated SD ±547, i.e. ±10.7 % (mine), NOT printed** |
| **phenylacetaldehyde** | **6097 ± 303** | ug/kg powder | as above | Table 3 row 5 | **level_only** |
| **benzaldehyde** | **37 201 ± 7850** | ug/kg powder | additional dairy method | Table 3 row 9 | **level_only** — RSD 21.1 %, the worst in the table |
| **(E)-2-octenal** | **907 ± 8** | ug/kg powder | as above | Table 3 row 7 | **level_only** — 88x Trikusuma's derived value (mine); **do not pool** |
| **heptanal** | **1326 ± 130** | ug/kg powder | as above | Table 3 row 10 | **level_only** |
| **2-methylbutanal** | **267 ± 8** | ug/kg powder | as above | Table 3 row 11 | **level_only** |
| **(E,E)-2,4-decadienal** | **101 ± 10** | ug/kg powder | as above | Table 3 row 4 | **level_only** |
| **(E,E)-2,4-nonadienal** | **53 ± 7** | ug/kg powder | as above | Table 3 row 6 | **level_only** |
| **methional** | **37 ± 5** | ug/kg powder | as above | Table 3 row 14 | **level_only** — **the ONLY sulfur compound in the entire paper** |
| **diacetyl** | **316 ± 25** | ug/kg powder | additional dairy method | Table 3 row 8 | **level_only** |
| **acetoin** | **978 ± 7** | ug/kg powder | additional dairy method | Table 3 row 27 | **level_only** |
| **vanillin** | **561 ± 18** | ug/kg powder | as above | Table 3 row 18 | **level_only** |
| **acetic acid** | **262 037 ± 10 440** | ug/kg powder | additional dairy method | Table 3 row 15 | **level_only** |
| the remaining 12 Table 3 rows | see Table 3 above, all with ± SD and RSD | ug/kg powder | as above | Table 3, p. 9 | **level_only** |
| **odour threshold in water, 3-methylbutanal** | **0.5** | ug/kg | orthonasal, water | Table 3 col. OT | **threshold `[C]`** — Leibniz-LSB@TUM database, ref. [29] |
| **odour threshold in water, hexanal** | **2.4** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, acetaldehyde** | **16** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, (E,E)-2,4-decadienal** | **0.027** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** — the lowest in the paper |
| **odour threshold in water, (E,E)-2,4-nonadienal** | **0.046** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, methional** | **0.43** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, (E)-2-octenal** | **1.7** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, diacetyl** | **0.96** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, benzaldehyde** | **150** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, heptanal / 2-methylbutanal / (E)-2-undecenal / phenylacetaldehyde** | **6.1 / 1.5 / 0.78 / 5.2** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, vanillin / nonanoic acid / acetic acid / 3-methylbutanoic acid / decanoic acid** | **53 / 26 / 5600 / 490 / 3.5** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, hexanoic / octanoic / phenylacetic / 2-methylbutanoic acid** | **4800 / 190 / 68 / 3100** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, γ-octalactone / 2-undecanone / acetoin** | **6.5 / 24 / 590** | ug/kg | orthonasal, water | Table 3 | **threshold `[C]`** |
| **odour threshold in water, 2,3-octanedione** | **29** | ug/kg | orthonasal, water | Table 3 + §2.6 | **threshold `[C?]` — AMBIGUOUS.** §2.6 says determined in water for this work; Table 3 footnote (c) says database. **Do not record as measured** |
| **odour threshold in water, (E,E)-3,5-octadien-2-one** | **27** | ug/kg | orthonasal, water | Table 3 + §2.6 | **threshold `[C?]` — AMBIGUOUS**, same contradiction |
| odour threshold, (E)-2-dodecenal | — | — | §2.6 names it as determined in water for this work | **never printed in the article** (Table S5, off disk) | **absent** |
| **OAV, all 27 compounds** | see Table 3, **10 186 down to 1.7** | — | concentration per kg of powder over a water threshold | Table 3 col. OAV | **derived_assumption** — a ratio of a `[M]` numerator to a `[C]` denominator on **mismatched bases** (powder vs cup). **~10x above the in-cup OAV (mine)** |
| training reference solutions, presented at **10-fold** their thresholds | hexanal 25.0 ug/L; 3-isopropyl-2-methoxypyrazine 0.1 ug/L; acetic acid 60 mg/L; (E,E)-2,4-decadienal 0.32 ug/L; phenylacetic acid 0.68 mg/L; 2-ethyl-5-methylpyrazine 1.0 mg/L; 3-methylbutanal 5.0 ug/L; 2,3,5-trimethylpyrazine 0.12 mg/L | in water | orthonasal, 20 mL | §2.5, p. 5 | **threshold `[C]`, implied** — dividing by 10 gives implied water thresholds of **hexanal 2.5 ug/L** (against Table 3's 2.4), **3-methylbutanal 0.5 ug/L** ✓, **(E,E)-2,4-decadienal 0.032 ug/L** (against Table 3's 0.027), **3-isopropyl-2-methoxypyrazine 0.01 ug/L**, **2-ethyl-5-methylpyrazine 100 ug/L**, **2,3,5-trimethylpyrazine 12 ug/L**, **acetic acid 6000 ug/L** (against Table 3's 5600), **phenylacetic acid 68 ug/L** ✓ (all mine). **Two of the four checkable ones do NOT match Table 3** (Flags 2) |
| **omission-test significance, all 15 tests** | see Table 4: **9 significant, 4 not, 2 blocks not** | p in %, 3-AFC, 13-14 panelists | orthonasal, 10 % deodorised protein + recombinant in water/triacetin 97.5:2.5 | Table 4, p. 11 | **measured_ratio** — the sensory evidence, and the strongest in the pea corpus |
| **the nine key food odorants** | 3-methylbutanal, hexanal, acetaldehyde, (E,E)-2,4-nonadienal, (E)-2-octenal, benzaldehyde, heptanal, 2-methylbutanal, nonanoic acid | — | as above, generalised across all ten isolates by heatmap clustering | Abstract, §3.3, §4, Fig. 6 | **structural_gate** — a validated membership list for what an unprocessed pea isolate smells of |
| **OAV rank does NOT predict sensory necessity** | the #4 and #5 OAVs ((E,E)-2,4-decadienal 3736, phenylacetaldehyde 1173) are **not significant** on omission; the #10 (heptanal 217) is **p < 0.1 %** | — | same panel, same recombinant | Table 4 vs Table 3 **(mine)** | **measured_bound** — a within-study measurement of how far an OAV ranking can be trusted |
| **protein content of ten commercial pea isolates** | **42, 67, 68, 68, 68, 68, 69, 70, 71, 72 %** (B, A, C, E, F, G, J, D, I, H) | % | **Dumas, N x 5.4** | Table 1, p. 3 | **measured_ratio** — directly usable as a `MATRIX_LOADING` basis, and the only measured pea protein contents in the corpus |
| **pyrazines are BELOW LOQ in all ten pea isolates** | eight pyrazines, none detected even at raised loading; LOQs from Table S2's method | — | unheated commercial isolates | §3.1, p. 7 | **measured_bound** — a **measured negative** on the nitrogen-heterocycle channel in the unheated starting state |
| **(E,Z)-2,6-nonadienal and 4-ethylbenzaldehyde not identified in ANY of the ten proteins** | below LOQ (0.8 and 15.8 nmol/L respectively) | — | as above | §3.2 p. 8 + Table 2 | **measured_bound** |
| **cross-study hexanal agreement, Utz vs Trikusuma, per kg of isolate** | **1.35x** (14 886 vs 11 000 ug/kg) | — | two labs, two methods, two commercial isolates, both unheated | **(mine)**, using `trikusuma2020_extraction.md` §3 item 4 | **measured_ratio** — cross-study, and the strongest levels anchor for pea hexanal in the corpus |
| **cross-study alkenal DISagreement, Utz vs Trikusuma** | **(E)-2-octenal 88x; (E,E)-2,4-decadienal 51x** | — | as above | **(mine)** | **derived_assumption** — cross-study, cross-method; the disagreement, not the values, is the finding |
| the ten isolates' full concentration and OAV matrices | — | — | — | **Tables S4 and S5, NOT on disk** | **absent** |
| every aroma-profile intensity, every heatmap cell | — | — | — | Figs. 4 and 5 | **figure_only** |

### Can these be put on the same basis as what the repository already carries?

**(a) The thresholds: no.** They are water thresholds and 24 of 27 are second-hand. The matrix layer
wants a **paired water/matrix threshold** and this paper has none. **Utz does not unlock any
matrix-corrected threshold.** It could, however, replace or corroborate individual water-threshold
entries — but only by importing the Leibniz-LSB@TUM database as a source in its own right, with its
own panel and criterion, which is a decision this dossier does not make.

**(b) The levels: yes, on a per-kilogram-of-isolate basis, and that basis must be stated.**
Trikusuma's derived per-isolate column and Utz's Table 3 are on the same basis and hexanal agrees to
1.35x. Any Utz level entered into a levels table must carry **`basis: per_kg_dry_isolate`**, not
per litre, and must not be compared to a beverage concentration without the 3 % (Trikusuma) or 10 %
(Utz sensory) conversion.

**(c) The protein contents: yes, and they are better than what the corpus has.** Ten commercial
isolates measured by Dumas at N x 5.4 with lot numbers. `MATRIX_LOADING` currently derives its pea
entry (`pea_protein_1pct`, 10 g/L) from Bi 2022's *isolate mass*, with **no protein content
anywhere**; Sun 2025's is a supplier claim of 90 %. **Utz's 68 % for Nutralys S85F is the kind of
number that should sit behind any pea loading**, and note that 68 % is a long way from the 90 %
Sun's supplier claims and from the implicit 100 % that "10 g/L of isolate" assumes.

**(d) The omission results are sensory evidence, not a parameter.** They belong in the reporting
layer as a check on the OAV ranking, and specifically as the answer to "how far can an OAV ranking
be trusted?" — which this paper measures at **four inversions in thirteen tests**.

**(e) Nothing here is a rate, a binding constant, a partition coefficient or an activation energy.**
There is no heating step anywhere in this paper. It is a starting-state decode.

## 5. Flags

1. **Every threshold is a WATER threshold and the OAVs are on a mismatched basis.** Concentrations
   are per kilogram of dry powder; thresholds are per kilogram of water; the panel smelled a 10 %
   suspension. **The tabulated OAVs are therefore about 10x the OAVs of the stimulus actually judged
   (mine)**, and they take no account of protein binding, which the corpus measures at
   1.8e-2 to 2.5e-1 L/g for pea protein on aldehydes (Bi 2022, Sun 2025). **An Utz OAV is an upper
   bound on perceived intensity, not an estimate of it.**
2. **The QDA training solutions imply water thresholds that disagree with Table 3 for two of the
   four checkable compounds.** §2.5 presents references at "10-fold odor thresholds": hexanal at
   25.0 ug/L implies **2.5** against Table 3's **2.4**; 3-methylbutanal at 5.0 ug/L implies **0.5**
   ✓; (E,E)-2,4-decadienal at 0.32 ug/L implies **0.032** against Table 3's **0.027** (a **19 %
   disagreement, mine**); acetic acid at 60 mg/L implies **6000** against Table 3's **5600** (a
   **7 % disagreement, mine**). Small, but it means **the paper is using at least two threshold
   values for the same compound in the same study**.
3. **Several analytes are quantified against a structurally distant internal standard.** Table 3's
   "Used IS" column: **methional** (a sulfur aldehyde) and **(E)-2-octenal** (an alkenal) are both
   quantified against **hexanal-d12** (an n-alkanal); **benzaldehyde** and **phenylacetaldehyde**
   against **phenylacetic acid-13C2** (an acid); **(E,E)-3,5-octadien-2-one** and **2,3-octanedione**
   against **diacetyl-d6** (a C4 diketone); **2-undecanone**, **(E)-2-undecenal**,
   **(E,E)-2,4-decadienal** and **(E,E)-2,4-nonadienal** all against **decanal-d2**. SIDA's whole
   virtue is that the standard co-elutes and co-ionises with the analyte; **eight labelled standards
   cannot do that for 27 analytes**. The matrix-effect validation of §3.1 was checked on four
   representatives only (hexanal, hexanoic acid, (E,Z)-2,6-nonadienal, 2,3-octanedione) and the
   conclusion "all analytes indicated the same behavior" is asserted, not shown.
4. **Eleven of the 27 rows have no validation data at all.** Everything marked footnote (f) —
   acetaldehyde, diacetyl, benzaldehyde, nonanoic acid, acetic acid, 3-methylbutanoic acid, decanoic
   acid, octanoic acid, phenylacetic acid, 2-methylbutanoic acid, acetoin — came from the
   "additional 3-NPH-UHPLC-MS/MS method published for dairy analysis [16]" and **appears nowhere in
   Table 2**. **No recovery, no LOD, no LOQ is printed for any of them in this paper**, including
   **acetaldehyde (OAV 4512, the #3 odorant and a named KFO)** and **benzaldehyde (OAV 248, a named
   KFO)**. Their standards (acetaldehyde-d3, acetic acid-13C2, butyric acid-13C4, octanoic acid-d15)
   also do not appear in the §2.4 IS list, which names only eight. **Chase ref. [16] before using
   any footnote-(f) row.**
5. **The derivatisation dictates the compound list, and the list has a hole shaped like the
   repository's main interest.** 3-NPH + EDC tags **carbonyls and carboxylic acids**. A thiol, a
   sulfide, a disulfide, a furan, a thiophene, a thiazole, an alcohol or a pyrazine without a
   carbonyl is **invisible by construction**. So the finding "eight of nine KFO are aldehydes" is
   partly a statement about the assay. **2-methyl-3-furanthiol, 2-furfurylthiol, 2-pentylfuran,
   1-octen-3-ol and every alcohol Trikusuma reports are absent from this paper because they could
   not have been seen**, not because they are not there. The pyrazine screen was a separate,
   underivatised extraction and is not subject to this.
6. **The pyrazine negative is a limit-of-quantitation statement, not a zero.** §3.1 says the
   pyrazines "were below the LOQ and thus they were excluded from further analysis". **No LOQ for
   any pyrazine is printed** (they are in Table S2, off disk), so the negative has no number
   attached. Trikusuma quantifies 2,5-dimethylpyrazine at 2.46 ug/L in an unheated pea beverage,
   which on a per-isolate basis is ~82 ug/kg — **not obviously below a plausible LC-MS/MS LOQ**.
   The two results may be reconcilable (different isolates, different methods) or may not.
7. **It is not stated whether the reported concentrations are recovery-corrected.** Recoveries run
   80.5-106.9 %. With SIDA the labelled standard added before workup should absorb extraction
   losses, making the recovery figure a process check rather than a correction — **but that is my
   inference and the paper never says it**. If the values are uncorrected, hexanal's true content
   is 20 % higher than tabulated (mine).
8. **Table 3 and the running text print different OAVs for two compounds.** (E,E)-2,4-decadienal:
   **3736** in Table 3, **3741** in §3.2, **3736** in Table 4. (E,E)-2,4-nonadienal: **1156** in
   Table 3, **1157** in §3.2 and in §3.3. (E)-2-undecenal: **157** in Tables 3 and 4, **158** in
   §3.3. All three are display-rounding artefacts under 0.4 % (mine, §3 item 1) and none changes any
   conclusion, but a reader quoting a single OAV should say which line they took it from.
9. **Proteins H and I share batch number 41266** while being sold as "Empro E86HV" and "Empro E86"
   with protein contents of 72 % and 71 %. Either the batch numbering is per-production-campaign
   rather than per-grade, or one entry is a transcription error. **H and I also fall in the same
   heatmap cluster (cluster 2)**, so nothing downstream depends on it; but they should not be
   treated as two independent samples.
10. **Maltol is absent from this 27-row decode and Trikusuma reports it at 48 000 ug/kg on a
    per-isolate basis.** Maltol has a compound-registry id in this repository. Maltol has a carbonyl
    and should have been derivatised and seen by this method; it is not in Table 2, not in Table 3,
    and not in the §3.1 analyte list. **Either it was not searched for, or it was below OAV 1 in
    protein C** (its water threshold is high, so a large concentration can still give a small OAV).
    The paper does not say which. **This is a real discrepancy between the corpus's two pea-isolate
    levels sources and it should be carried, not smoothed.**
11. **pH is never stated, anywhere.** Not for the extraction (which is half acetonitrile with
    pyridine and formic acid downstream), not for the 10 % sensory suspension. The corpus's binding
    and adduct layers are pH-gated (`PH_ADDUCT_GATE`; Leksrisompong's diacetyl binds at pH 7 and not
    at pH 5.5), so a pea suspension of unstated pH cannot be matched to any of them.
12. **Only ONE of the ten isolates was decoded sensorially.** All the recombination, all fifteen
    omission tests and the entire aroma profile are on protein C. The generalisation to the other
    nine rests on **heatmap clustering of OAVs alone (Fig. 5)** — no panel ever smelled A, B or
    D-J. The paper says "there were strong indications that the examined KFO of pea protein C were
    transferable"; that hedge should be preserved.
13. **The omission tests had 13-14 panelists, not 16.** §2.5 describes a 16-member panel; Table 4's
    footnote and §2.5 both say the 3-AFC tests were "based on the number of correctly identified
    samples and attended panelists (13-14)". At n = 13 in a 3-AFC, the significance boundary is
    coarse — the p values of 4.0 %, 4.8 % and 2.6 % that carry five of the nine KFO designations sit
    within one or two correct answers of non-significance. **The five single-asterisk KFO are the
    fragile ones**; only 3-methylbutanal, heptanal (both p < 0.1 %) and 2-methylbutanal (p = 0.7 %)
    are robust to a single panelist.
14. **The deodorised base is not characterised.** Pentane and dichloromethane extraction of 500 g of
    protein removes lipid as well as volatiles and may change the protein's binding surface, its
    surface hydrophobicity and hence its retention of the re-added odorants. The only check reported
    is that panelists "could not sensorially relate it to pea". **The recombinant's matrix is
    therefore not the same matrix as the reference sample C**, and the recombination's success is
    measured against that difference, not controlled for it.
15. **Triacetin at 2.5 % is a solvent in every sensory sample**, including the reference C. It
    cancels between reference and recombinant. It does not cancel against water thresholds taken
    from a database measured in plain water.
16. **The OAV ranking is measurably unreliable, by this paper's own data.** Four of thirteen
    individually tested odorants invert (§3 item 6). **Anyone using Table 3's OAV column as a
    priority ranking is using a quantity this paper itself falsified in four of thirteen cases.**
17. **What this paper does NOT contain**: any heating step; any threshold in a protein matrix; any
    paired water/matrix threshold; any binding constant, partition coefficient or retention
    measurement; any rate; any pH; any sulfur compound other than methional; any thiol; any pyrazine
    above LOQ; any furan; any alcohol; any 2-pentylfuran or 1-octen-3-ol; any retronasal
    measurement; any per-litre concentration; any of the ten-protein data as printed numbers.
18. **What to request from the authors or fetch online**: (i) **Tables S4 and S5** — the
    concentrations and OAVs for all ten isolates, which would turn this from one decode into a
    ten-lot distribution and is the single highest-value missing item; (ii) whether the two ambiguous
    thresholds (2,3-octanedione 29, (E,E)-3,5-octadien-2-one 27 ug/kg) were measured here or taken
    from the database, and by what panel and criterion; (iii) the (E)-2-dodecenal threshold, which
    §2.6 says was determined for this work and is printed nowhere; (iv) whether Table 3's
    concentrations are recovery-corrected; (v) the validation data (recovery, LOD, LOQ) for the
    eleven footnote-(f) analytes; (vi) the Figure 4 intensity ratings and the Figure 5 heatmap cells
    as numbers; (vii) whether maltol was searched for; (viii) the pH of the 10 % sensory suspension.
    Note the paper is **open access CC BY** and the supplementary is at
    `https://www.mdpi.com/article/10.3390/foods11030412/s1` — **S4 and S5 are freely downloadable
    and should be fetched before this dossier is used for anything quantitative about the ten-lot
    spread.**
