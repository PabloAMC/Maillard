# Suppavorasatit 2012 — EXTRACTION (soymilk from IA 3027 soybeans, ~3 % protein, pH 6.44-6.74; enzymatic deamidation by protein-glutaminase, 44 C / 40 U per g protein / 2 h, giving 66.4 % DD and 4.25 % DH; orthonasal odour DETECTION thresholds for vanillin and maltol measured **in soymilk** by ASTM E697-04 3-AFC on 21 and 27 panelists; Fechner and Stevens dose-response curves on a separate 12-member panel; protein solubility at pH 3.0 / 5.0 / 7.0)

### THE ANSWER TO THE QUESTION THIS PAPER WAS PULLED FOR: **YES, the odour detection thresholds WERE measured in soymilk — both of them, and by a proper ASTM forced-choice protocol — but there is NO WATER LEG.** The pairing is **control soymilk against DEAMIDATED soymilk** (9.61 vs 1.80 ug/mL vanillin; 23.8 vs 7.01 ug/mL maltol), i.e. protein-versus-modified-protein, **not protein-versus-water**. Both legs contain ~3 % soy protein. **This is therefore NOT the paired water/matrix threshold the matrix layer has been refusing every matrix-corrected threshold for want of**, and no `K_g`, no matrix-shift factor and no matrix-corrected threshold can be built from it. What it IS: a measured, same-panel, same-matrix demonstration that changing the protein's chemistry alone moves an odour threshold **5.3x for vanillin and 3.4x for maltol (mine)**.

**Source on disk:** `data/articles/suppavorasatit2012.pdf` (7 pp., J. Food Science 78(1) 2013,
pp. C1-C7). Read from the `pdftotext -layout` text layer
(`scratchpad/suppavorasatit2012.txt`), **whole file**. **Tables 1, 2, 3, 4 and 5 came through clean
and are re-typed in full below — that is every table in the paper; there is no supplementary
material.** Figures 1 (solubility bar chart), 2A/2B (vanillin Fechner and Stevens plots) and 3A/3B
(maltol Fechner and Stevens plots) are images — but **the four Stevens regression equations and
their R^2 are printed as text inside the Fig. 2B and 3B panels and came through the text layer, so
they ARE typed here**; every solubility bar and every point on every dose-response curve is
**figure-only** except the two Fechner intensity readings quoted in the running text. Repo status
before this dossier: Suppavorasatit 2012 is cited nowhere in `src/kinetic_core/` and has no
extraction dossier. **`maltol` has a compound-registry id in `data/keys/compounds.yml` (line 949) as
of today**; `vanillin` is also keyed.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Enzymatic Protein Deamidation on Protein Solubility and Flavor Binding Properties of Soymilk" |
| Authors | Inthawoot Suppavorasatit (Dept. of Food Technology, Faculty of Science, Chulalongkorn Univ., Phayatai Rd., Wangmai, Pathumwan, Bangkok 10330, Thailand); Soo-Yeun Lee and Keith R. Cadwallader (corresponding, cadwlldr@illinois.edu) — Dept. of Food Science and Human Nutrition, Univ. of Illinois at Urbana-Champaign, 1302 West Pennsylvania Ave., Urbana, IL 61801 |
| Venue | **Journal of Food Science, Vol. 78, Nr. 1, 2013, pp. C1-C7**, section "C: Food Chemistry". MS 20120601, **submitted 28 April 2012, accepted 21 October 2012**; copyright line "© 2012 Institute of Food Technologists". **The file is named `suppavorasatit2012.pdf` and the copyright is 2012, but the issue is 78(1) 2013** (Flags 12) |
| DOI | **`doi: 10.1111/j.1750-3841.2012.03012.x`** — printed exactly so at the foot of p. C1 |
| Funding | USDA National Inst. of Food and Agriculture, **Hatch project ILLU-698-366**; Royal Thai Government Scholarship (I. Suppavorasatit) |
| Ethics | **Univ. of Illinois Institutional Review Board, IRB Protocol Number 12038** |
| The two ligands | **vanillin** (4-hydroxy-3-methoxybenzaldehyde) and **maltol**, printed in this paper as "3-hydroxy-2-methoxyl-4H-pyran-4-one" — **that name is wrong**; maltol is 3-hydroxy-2-**methyl**-4H-pyran-4-one (Flags 11). Both **food grade, >98 % purity**, from FONA Intl. Inc. (Geneva, Ill.) |
| The matrix | Soymilk made from **soybean variety IA 3027** (Clarkson Grain Co., Cerro Gordo, Ill.), chosen "because of its high protein content and higher yield", by the method of Lozano et al. 2007; hot-filled into polyethylene bottles, ice-water cooled, capped, stored at **4 ± 1 C** |
| The enzyme | **Protein-glutaminase (PG) "Amano" 500, 500 U/g**, Amano Enzyme Inc. (Elgin, Ill.) |
| Naming | **FSM** = finished (untreated) soymilk; **CSM** = control soymilk (same 44 C / 2 h treatment, **no enzyme**); **DSM** = deamidated soymilk; **DD** = degree of deamidation; **DH** = degree of hydrolysis; **BET** = best estimate threshold; **group BET** = geometric mean of the individual BETs; **SPI/DSPI** = soy protein isolate / deamidated SPI, from the companion paper |
| Companions on disk | `damodaran1981_extraction.md` (soy binding, the registry's `soy_protein` rows), `guo2020_extraction.md` (soy isolate headspace), `leksrisompong2010_extraction.md` (the caseinate threshold/partition pair), `k4b_paired_thresholds_and_browning.md` (the repository's standing register of paired thresholds) |

## 1. Why it matters

**The question that sent us here, and its answer.** The matrix layer of
`src/kinetic_core/parameters_matrix.py` refuses every matrix-corrected threshold, and it says why:
k2 sec. D.1 measured that matrix-to-water threshold ratios span 2000x with a 1-sigma band of 27-41x
cross-study, so a general correction factor is refuted (a uniform 33x misplaces the two extreme
compounds by 10x and 28x **in opposite directions**), and there is a lookup table with an explicit
`no_measured_threshold` state instead. What would unlock a row in that table is **a threshold
measured in water and the same threshold measured in a plant-protein matrix, by the same panel with
the same criterion.** The task's hypothesis was that this paper supplies exactly that.

**It does not, and the reason is worth stating precisely.** The thresholds here are **real, measured,
in soymilk, by a named ASTM protocol**, and there are two of them per compound. But the two legs are:

- **CSM** — soymilk held 2 h at 44 C **without enzyme**. ~3 % soy protein. Threshold measured.
- **DSM** — the same soymilk held 2 h at 44 C **with protein-glutaminase**. ~3 % soy protein,
  66.4 % of its glutamine residues converted to glutamate. Threshold measured.

**Both legs are soymilk. Neither is water.** The words "water" and "odorless distilled water" appear
in this paper only in three places, and none of them is a threshold determination:

1. §"Flavor compounds": *"Vanillin and maltol solutions were prepared separately in odorless
   distilled water"* — that is the **stock solvent**, 1 mL of which is spiked into 14 mL of soymilk.
2. §"Sensory evaluation": *"Panelists were instructed to sniff odorless distilled water between
   aroma evaluations"* — a **palate/nose cleanser**.
3. §"Aroma intensity scaling": *"references (50 μg/mL vanillin and 250 μg/mL maltol in odorless
   distilled water) were developed"* — **anchor standards for the 15-point intensity scale**, not
   thresholds, and both are far above any plausible water threshold.

**The sentence that settles it**, from §"Determination of odor detection thresholds", quoted in
full:

> *"The test samples (flavored soymilk) were prepared by adding the same volume (1 mL) of various
> flavor compound solutions into 14 mL of either DSM or CSM (Table 1). Test samples were
> refrigerated (4 ± 1 °C) overnight (24 h) to allow for equilibration of flavor-matrix interactions
> before testing."*

And the confirming sentence from the Results, §"Effect of deamidation on odor detection threshold of
vanillin and maltol in soymilk":

> *"In this study, odor detection thresholds were measured and used to demonstrate the effect of
> deamidation on protein-flavor binding in soymilk."*

**So: measured in soymilk, yes. Against a water reference, no.** The only non-soymilk comparison in
the paper is a **citation** to Karagül-Yüceer et al. 2004 for the BET of vanillin (7.41 ppm) and
maltol (16.6 ppm) **in skim milk** — another protein matrix, from another laboratory, with another
panel. There is no water threshold, measured or cited, anywhere in this article.

**What the repository can therefore take, and what it cannot.**

- **CANNOT**: a `K_g`, a matrix shift, a water/matrix threshold ratio, or any entry in the
  `no_measured_threshold` lookup that would let the matrix layer emit a matrix-corrected threshold
  for vanillin or maltol. Every construction of that kind needs a protein-free leg and there is
  none. **The matrix layer's refusal stands, and this paper does not lift it.**
- **CAN**: two **absolute measured thresholds in a real 3 %-protein plant-protein beverage**, which
  is more than the corpus has for most compounds and is directly usable as a `threshold` record with
  `medium = soymilk`. Plus a **within-study, same-panel, same-matrix measurement of how far a
  protein-chemistry change alone moves a threshold**: **5.3x for vanillin, 3.4x for maltol (mine)** —
  a `within_study_ratio` of the cleanest possible construction, since the panel, the protocol, the
  base soymilk, the storage, the temperature and the serving vessel are all identical between legs
  and only the glutamine chemistry differs.
- **CAN, with care**: the paper's own mechanism claim. It attributes the CSM/DSM threshold gap to
  **loss of Schiff-base (covalent) capacity**: deamidation *"reduces the number of available amide
  groups by converting most of glutamine residues in soy protein to glutamic acids. Therefore, the
  potential to form covalent bonds (Schiff bases) with the carbonyl group of vanillin was
  decreased."* **That mechanism is chemically confused** — a Schiff base forms with a primary amine
  (lysine epsilon-amino), not with a glutamine side-chain amide, and protein-glutaminase does not
  touch lysine (Flags 3). The *observation* is sound; the *explanation* should not be imported.

**Why maltol specifically matters today.** `maltol` has a compound-registry id in this repository.
Its only other appearance in the pea/soy corpus is Trikusuma 2020, where maltol is quantified in a
pea UHT beverage (48 000 ug/kg on a derived per-isolate basis) and is a Maillard/caramelisation
product the model's own chemistry produces. **This paper gives maltol a measured orthonasal
detection threshold in a plant-protein beverage — 23.8 ug/mL = 23 800 ug/L in control soymilk** —
which is the first such number for maltol in the corpus. Note how high it is: three to four orders
of magnitude above the aldehyde thresholds in Utz 2022's water column. **Maltol is a weak odorant
per unit mass, and this paper measures how weak, in a matrix that resembles the model's own.**

## 2. Methods as they matter to a model

- **The pot (thresholds).** **1 mL of aqueous flavour stock + 14 mL of soymilk = 15 mL**, in a
  **125 mL FEP (Teflon) squeeze sniff bottle** (Nalge Nunc), covered with aluminium foil, 3-digit
  random code. So the matrix is **93.3 % soymilk, 6.7 % water (mine)**.
- **The pot (dose-response).** Same construction: 1 mL of flavour stock into 14 mL of CSM or DSM in
  sniff bottles.
- **The pot (aroma/taste difference tests).** 15 mL of soymilk in sniff bottles (aroma); 20 mL in
  1 oz (30 mL) translucent plastic cups with lids (taste).
- **Protein loading.** **"approximately 3% protein"** — stated in the Discussion as the reason
  soymilk and skim-milk BETs are similar (*"both fluids contain about the same amount of protein
  (approximately 3% protein)"*). That is **~30 g/L (mine)**. **It is an approximation offered in
  passing, not a measured composition of this batch**, and no total-protein figure for FSM, CSM or
  DSM is printed anywhere (the Lowry/DC assay was used for the *solubility* denominator only).
- **pH.** **FSM 6.74, CSM 6.44** — both printed in the solubility discussion. **The pH of DSM is
  never stated**, which is a real omission because deamidation converts amides to carboxylates and
  necessarily acidifies (Flags 5). Solubility was measured separately at pH 3.0, 5.0 and 7.0 in
  0.1 M acetate-phosphate buffer.
- **Temperature and time.**
  - Deamidation: **44 C for 2 h**, E/S = **40 U per g protein**. CSM: identical 2 h at 44 C, no
    enzyme.
  - Flavour-matrix equilibration: **refrigerated at 4 ± 1 C overnight (24 h)** — "to allow for
    equilibration of flavor-matrix interactions". For the dose-response samples, "at least 24 h".
  - Serving: **removed from refrigeration and conditioned at room temperature (25 C) for
    approximately 1 h** before evaluation.
  - Storage of soymilk: **4 ± 1 C**.
- **Threshold method and its family — read this before transferring the numbers.**
  - **Protocol: ASTM E697-04**, "Standard practice for determination of odor and taste thresholds by
    a forced-choice ascending concentration series method of limits", *"with some modifications
    (Watcharananun and others 2009)"*.
  - **Route: ORTHONASAL.** Printed explicitly: *"Odor detection thresholds were determined
    orthonasally using ASTM E697–04 protocol"*. Sniff bottles, foil-covered.
  - **Criterion: DETECTION, not recognition.** The task is *"identify the sample with the strongest
    odor"* / *"select the odd sample"* in a triad where two are unflavoured soymilk — that is a
    difference-from-blank detection task, not an identification of the odour quality.
  - **Design: 3-alternative forced choice (3-AFC)**, **6 sets of 3 samples**, each set = 2 blanks +
    1 flavoured, **served in ascending flavour concentration**, randomised within each set, guessing
    required when uncertain. Panelists were **told** the sets ascend.
  - **Statistic: group BET = the geometric mean of the individual BETs** (ASTM 2004).
  - **Panel sizes: 21 panelists for vanillin (5 male, 16 female, 22-48 y); 27 panelists for maltol
    (5 male, 22 female, 19-40 y).** Trained the day before in **2 practice sessions of 6 sets of
    3-AFC** each.
  - **The blanks are the corresponding soymilk** (CSM blanks for the CSM series, DSM blanks for the
    DSM series) — *"2 of them were plain soymilks (no flavor added)"*. So the panel is discriminating
    flavoured soymilk from unflavoured soymilk, and **the soymilk's own beany background is present
    in every sample including the blanks** (Table 3 lists "Beany" as a panel-generated aroma term).
- **The concentration ladders — and they are NOT the same for the two soymilks on vanillin.**
  Table 1 gives six sets at x/27, x/9, x/3, x, 3x, 9x. **Vanillin: CSM anchored at x = 45.00 ug/mL
  (ladder 1.667 to 405.0); DSM anchored at x = 15.00 ug/mL (ladder 0.556 to 135.0)** — a 3-fold
  shift of the whole ladder. **Maltol: CSM and DSM ladders are IDENTICAL** (x = 35.00, ladder 1.297
  to 315.0). The footnote says *"Concentrations were based on the literature and preliminary
  testing."* (Flags 2).
- **Dose-response method (a separate panel and a separate quantity).** **12 panelists** (2 male,
  10 female, 22-45 y) with prior aroma-intensity experience, plus **5 additional hours of training
  in five 1-h sessions**. Attributes **"vanilla"** and **"cotton candy"** by consensus. **15-point
  intensity scale.** Reference anchors: **50 ug/mL vanillin and 250 ug/mL maltol in odorless
  distilled water**. Five ascending concentrations in **5-fold increments**: **CSM 1000, 200, 40, 8,
  1.6 ug/mL; DSM 500, 100, 20, 4, 0.8 ug/mL** — again **different ladders for the two soymilks**,
  this time for both compounds.
  - **Fechner's law**: log(concentration) vs perceived intensity, expected sigmoidal.
  - **Stevens's power law**: `R = k*C^n`; `n` is the slope of log(intensity) vs log(concentration),
    fitted on **"the linear section of the plots"** only.
- **Difference testing.** **2-AFC with the warm-up method** of Thieme & O'Mahony 1990, **17 panelists
  (4 male, 13 female, 23-48 y)**, evaluations **in triplicate by each panelist**, analysed by
  **beta-binomial statistics** (IFPrograms 7.3) with a null probability of 0.5.
- **DD and DH.** DD by ammonia release (Sigma-Aldrich ammonia assay kit), expressed as the ratio of
  ammonia released by PG to the total glutamine, the latter measured by ammonia released on
  **2 N sulfuric acid at 100 C for 4 h**. DH as the percentage of protein remaining dissolved after
  **0.2 N TCA** precipitation over the total dissolved protein after the same complete acid
  hydrolysis.
- **Solubility.** 100 uL of sample in 1 mL of **0.1 M acetate-phosphate buffer at pH 3.0, 5.0 or
  7.0**, 1.5 mL microcentrifuge tubes, **25 C overnight**, vortexed, **3000 rpm (1000 x g) at 10 C
  for 10 min**, supernatant assayed by **Lowry after detergent solubilisation (Bio-Rad DC)**.
  Triplicate.
- **Statistics.** ANOVA and LSD in SAS 9.2. **The significance convention is printed backwards
  throughout**: *"to determine significant differences among treatments (P > 0.05)"*, and Figure 1's
  caption says *"same upper case letters within the same sample across different pH values indicate a
  significant difference (P > 0.05)"*. (Flags 8.)

## 3. Tables re-typed

### Table 1 (p. C3). "Final concentration of flavor compounds in test samples for threshold evaluation."

Concentration (ug/mL):

| Flavor compound | Soymilk | Set 1 (x/27) | Set 2 (x/9) | Set 3 (x/3) | Set 4 (x^a) | Set 5 (3x) | Set 6 (9x) |
|---|---|---:|---:|---:|---:|---:|---:|
| Vanillin | Control (CSM) | 1.667 `[M]` | 5.000 `[M]` | 15.00 `[M]` | **45.00** `[M]` | 135.0 `[M]` | 405.0 `[M]` |
| Vanillin | Deamidated (DSM) | 0.556 `[M]` | 1.667 `[M]` | 5.000 `[M]` | **15.00** `[M]` | 45.00 `[M]` | 135.0 `[M]` |
| Maltol | Control (CSM) | 1.297 `[M]` | 3.890 `[M]` | 11.67 `[M]` | **35.00** `[M]` | 105.0 `[M]` | 315.0 `[M]` |
| Maltol | Deamidated (DSM) | 1.297 `[M]` | 3.890 `[M]` | 11.67 `[M]` | **35.00** `[M]` | 105.0 `[M]` | 315.0 `[M]` |

Footnote (a): *"Concentrations were based on the literature and preliminary testing."* (The `x`
column header carries the footnote marker.) These are **stimulus concentrations**, i.e. design
values, not measurements of anything, but they are printed and they bound the thresholds:
**a BET can only lie inside its own ladder.** DSM vanillin's ladder starts at 0.556 ug/mL and its
BET came out at 1.80; CSM vanillin's starts at 1.667 and its BET came out at 9.61.

### Table 2 (p. C3). "Degree of deamidation (DD) and degree of hydrolysis (DH) of soymilk after deamidation^a by protein-glutaminase."

| | Average^b ± Standard deviation |
|---|---:|
| DD (%) | **66.4 ± 2.6** `[M]` |
| DH (%) | **4.25 ± 0.42** `[M]` |

Footnotes: *"(a) Deamidation under optimal conditions: reaction temperature of 44 °C,
enzyme:substrate (E/S) ratio of 40 U/g protein for 2h. (b) n = 3."*

### Table 3 (p. C4). "Terms generated by panelists to differentiate aroma and taste of soymilks."

| Attribute | Terms generated | Attribute | Terms generated |
|---|---|---|---|
| Aroma | Beany | Taste | Bitter |
|  | Burnt sugar |  | Salty |
|  | Caramel |  | Sweet |
|  | Chickeny |  |  |
|  | Creamy |  |  |
|  | Vanilla |  |  |

Qualitative; no number. Recorded because **"Beany" and "Chickeny" are the panel's own descriptors
for the unflavoured soymilk background present in every threshold blank**, and "Burnt sugar" and
"Caramel" overlap the maltol descriptor "cotton candy" (Flags 6).

### Table 4 (p. C4). "Group best estimate thresholds (BET) of vanillin and maltol in control (treated without enzyme; CSM) and deamidated soymilk (DSM)."

| Flavor compound | BET (ug/mL) ± SD, **CSM** | BET (ug/mL) ± SD, **DSM** |
|---|---:|---:|
| **Vanillin** | **9.61 ± 1.39** `[M]` | **1.80 ± 2.44** `[M]` |
| **Maltol** | **23.8 ± 2.61** `[M]` | **7.01 ± 3.71** `[M]` |

**This is the table the paper was pulled for. Both columns are soymilk. There is no water column and
no water row anywhere in this article.**

**Note the standard deviations.** DSM vanillin is **1.80 ± 2.44** — **the SD is 136 % of the mean
(mine)**. DSM maltol is **7.01 ± 3.71**, SD **53 % of the mean (mine)**. Against CSM's **14 %** and
**11 %** respectively. **The two deamidated legs are far noisier than the two control legs**, and
the vanillin/DSM cell is the least reliable number in the paper (Flags 1). The SD of a group BET
that is itself a geometric mean is also an odd statistic — an arithmetic SD on a log-scale
statistic; the paper does not say whether it is the SD of the individual BETs or something else.

### Table 5 (p. C6). "Stevens's power law exponent (n) and coefficients of determination (r^2) from the plots of log of flavor compound concentration versus log of perceived aroma intensity of control soymilk (CSM) and deamidated soymilk (DSM)."

| Protein | Vanillin, n | Vanillin, r^2 | Maltol, n | Maltol, r^2 |
|---|---:|---:|---:|---:|
| **CSM** | **0.4777** `[F]` | 0.99 `[F]` | **0.5151** `[F]` | 0.99 `[F]` |
| **DSM** | **0.5525** `[F]` | 0.99 `[F]` | **0.5702** `[F]` | 0.99 `[F]` |

**All four exponents are below 1**, i.e. perceived intensity grows more slowly than concentration in
both soymilks for both compounds. **Both DSM exponents exceed their CSM counterparts.**

### Regression equations printed inside Figures 2B and 3B (text layer, so typed)

| plot | sample | equation | R^2 |
|---|---|---|---:|
| Fig. 2B, vanillin | **DSM** | `y = 0.5525x − 0.2076` | **0.9973** `[F]` |
| Fig. 2B, vanillin | **CSM** | `y = 0.4777x − 0.2119` | **0.9990** `[F]` |
| Fig. 3B, maltol | **DSM** | `y = 0.5702x − 0.2361` | **0.9983** `[F]` |
| Fig. 3B, maltol | **CSM** | `y = 0.5151x − 0.2476` | **0.9997** `[F]` |

The slopes match Table 5 exactly; **the intercepts appear nowhere else in the paper** and the R^2
values are printed to four decimals here against Table 5's two.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| **the paper's own summary of the threshold shift** | *"Odor detection thresholds for the flavor compounds vanillin and maltol were approximately **5 and 3 fold lower**, respectively, in DSM than in CSM."* | Abstract |
| **vanillin BET, CSM** | **9.61 ppm** = 9.61 ug/mL | Results, p. C4 |
| **vanillin BET, DSM** | **about 1.80 ppm**, *"substantially lower (approximately 5-fold)"* | Results, p. C5 |
| **maltol BET, CSM** | **23.8 ppm** | Results, p. C5 |
| **maltol BET, DSM** | **7.02 ppm** *(the text)*, *"approximately 3 times lower"* — **Table 4 prints 7.01** | Results, p. C5 vs Table 4 (Flags 9) |
| **vanillin BET in SKIM MILK** | **7.41 ppm** `[C]` — Karagül-Yüceer et al. 2004, *"Evaluation of the character impact odorants in skim milk powder by sensory studies on model mixtures"*, J Sens Stud 19:1-3 | Results, p. C4 |
| **maltol BET in SKIM MILK** | **16.6 ppm** `[C]` — same source | Results, p. C5 |
| why soymilk and skim milk should be comparable | *"The BET values of soymilk and skim milk may be similar because both fluids contain about the same amount of protein (**approximately 3% protein**)."* | Results, p. C4 — **this is the only protein loading anywhere in the paper** |
| protein solubility at **pH 7.0** | *"no difference in protein solubility among soymilks"*, **nearly 100 %** for FSM, CSM and DSM | Results, p. C4, Fig. 1 |
| protein solubility at **pH 5.0** | **DSM approximately 35 %**; FSM and CSM **approximately 22 to 23 %** | Results, p. C4, Fig. 1 |
| protein solubility at **pH 3.0** | **approximately 24 %** for all three, no difference | Results, p. C4, Fig. 1 |
| pH of FSM / CSM | **6.74 / 6.44** | Results, p. C4 |
| pI of soy protein | **approximately 4.5** `[C]` (Hamada & Marshall 1989) | Results, p. C4 |
| **aroma difference CSM vs DSM** | **NO significant difference**: *"P < 0.4443; estimated probability of the data, 0.4902; and power of the test, 3.7%"*; no over-dispersion (gamma 0.0000) | Results, p. C4 |
| **taste difference CSM vs DSM** | **SIGNIFICANT**: *"P < 0.0116; estimated probability of the data, 0.7059; and power of the test, 75.4%"*; significant between-subject dispersion (gamma 0.3389) | Results, p. C4 |
| the authors' own caveat on the taste result | *"might be caused by factors other than taste per se (for example, viscosity, mouthfeel, and aroma by mouth)"* | Results, p. C4 |
| **Fechner reading at 50 ug/mL vanillin** | **DSM intensity 4.4**; **CSM intensity 3.0** (15-point scale) | Results, p. C6, Fig. 2a |
| DD of deamidated **SPI** under the same conditions | **43.7 %** `[C]` (Suppavorasatit et al. 2011) — against **66.4 %** for soymilk here | Results, p. C3 |
| DH of deamidated SPI | **4.81 %** `[C]` — against **4.25 %** here | Results, p. C3 |
| **n·K for VANILLIN binding, SPI vs DSPI at 25 C** | **88.8 x 10^4 M^-1** vs **9.69 x 10^4 M^-1**, described as *"9-times greater (P > 0.05)"* | Results, p. C5 `[C]` — **from the companion paper Suppavorasatit & Cadwallader 2012 (JAFC 60:7817-23), NOT measured here.** Note the printed statistic reads "greater (P > 0.05)", i.e. **not significant** (Flags 8) |
| **n·K for MALTOL binding, SPI vs DSPI at 25 C** | **303 x 10^4 M^-1** vs **79.6 x 10^4 M^-1**, *"about 4 times higher"* | Results, p. C5-C6 `[C]` — same companion paper |
| dose-response concentration ladders | **CSM: 1000, 200, 40, 8, 1.6 ug/mL**; **DSM: 500, 100, 20, 4, 0.8 ug/mL** (for both vanillin and maltol) | Methods, p. C6 |
| every solubility bar, every dose-response point | **figure-only** (Figs. 1, 2A, 3A) | — |

### Arithmetic on the printed values (all mine)

**1. The threshold shifts, exactly.**

| compound | BET CSM | BET DSM | CSM/DSM (mine) | the paper says |
|---|---:|---:|---:|---|
| vanillin | 9.61 | 1.80 | **5.34x** | "approximately 5-fold" ✓ |
| maltol | 23.8 | 7.01 (Table 4) | **3.40x** | "approximately 3 times lower" ✓ |
| maltol | 23.8 | 7.02 (text) | **3.39x** | — |

**Both reproduce the paper's own words.** Note that the vanillin shift **5.34x** is well inside its
own noise: with a DSM SD of ±2.44 on a mean of 1.80, the one-SD interval on the DSM leg is
**−0.64 to 4.24 ug/mL**, which includes zero and gives a CSM/DSM ratio anywhere from **2.3x to
infinity (mine)**. **The maltol shift is the more defensible of the two**: its DSM interval is
3.30 to 10.72, giving **2.2x to 7.2x (mine)**.

**2. Against the cited skim-milk values (mine).** vanillin: 9.61 / 7.41 = **1.30x**; maltol:
23.8 / 16.6 = **1.43x**. **Soymilk's control thresholds sit 1.3-1.4x above skim milk's** for both
compounds, in the same direction, on two independent panels 8 years apart — a modest and consistent
cross-matrix agreement, which is what the paper claims when it says the values are "close" and
"similar". **These are matrix-to-matrix, not matrix-to-water**, and cannot be turned into a water
baseline.

**3. The `K_g` that CANNOT be built, written out so nobody builds it.** The registry stores
`K_g = (K_water/K_matrix − 1) / protein_g_per_L`. Substituting a **threshold** ratio for a partition
ratio and the CSM/DSM pair for a water/matrix pair, at ~30 g/L, would give:

| compound | (BET_CSM/BET_DSM − 1)/30 |
|---|---:|
| vanillin | 1.45e-1 L/g |
| maltol | 7.98e-2 L/g |

**Do not ship either number, and do not record them anywhere but here as a warning.** Three separate
things are wrong with the construction: (i) **DSM is not water** — it contains the same ~30 g/L of
protein as CSM, merely deamidated, so the denominator protein loading is not the difference between
the legs and the division by 30 is meaningless; (ii) **a threshold ratio is not a partition ratio** —
k4b sec. B refutes partition-derived thresholds three independent ways on matched samples, and this
would be the same error run backwards; (iii) the numbers happen to land in the same decade as the
shipped pea and dairy constants, which makes them **plausible-looking and therefore dangerous**.
Recorded so a later reader who rediscovers the arithmetic knows it was already considered and
rejected.

**4. What the Stevens exponents actually say (mine).** Ratios: vanillin **0.5525/0.4777 = 1.157x**;
maltol **0.5702/0.5151 = 1.107x**. Both DSM slopes are **10-16 % steeper** than their CSM
counterparts. **This is a small effect on a fitted slope with no error bars anywhere**: Table 5 gives
no SD, no confidence interval and no n for the regression, and the fits are on **five points**, of
which only "the linear section" was used — so possibly three or four. **A 10 % difference between two
slopes fitted to four points each, with no stated uncertainty, is not a measurement of anything.**
The r^2 values of 0.997-0.9997 look impressive and are a property of fitting a line to a
five-point log ladder, not evidence about the difference between the two lines.

**5. The Fechner reading, converted.** At 50 ug/mL vanillin, DSM rated **4.4** against CSM's **3.0**
on a 15-point scale: a ratio of **1.47x (mine)** and an absolute gap of **1.4 scale points**. Both
figures are single readings pulled off a curve by the authors; **no SD is given for either**, and
50 ug/mL is not one of the five concentrations actually presented to the panel (CSM 1000/200/40/8/1.6
and DSM 500/100/20/4/0.8), so **both numbers are interpolations off Fig. 2a and neither is a
measurement** (Flags 7).

**6. The companion paper's binding ratios, for completeness (mine, on cited values).**
vanillin n·K SPI/DSPI = 88.8/9.69 = **9.16x**; maltol n·K SPI/DSPI = 303/79.6 = **3.81x**. Set
against this paper's threshold ratios of **5.34x** and **3.40x**: **maltol agrees to 1.12x between an
instrumental binding measurement on an isolate and a sensory threshold measurement on a soymilk;
vanillin disagrees by 1.72x.** That maltol agreement is a genuinely interesting cross-method
corroboration and it is the strongest thing this paper contributes to the corpus's method-boundary
question. **But it is cross-paper, cross-substrate (SPI suspension vs soymilk), and cross-DD (43.7 %
vs 66.4 %), so it is a coincidence of magnitude, not a validation.**

**7. Where the thresholds sit relative to the concentrations the model would predict.** Maltol's
control-soymilk BET is **23.8 ug/mL = 23 800 ug/L**. Trikusuma 2020's derived per-isolate maltol
level in an unheated pea beverage is 48 000 ug/kg of isolate, which in its own 3 % beverage is
**1440 ug/L** — **16.5x below this threshold (mine)**. **On these two numbers maltol would be below
threshold in that beverage**, cross-matrix and cross-legume caveats fully applying. Worth recording
because maltol is registry-keyed and the model can produce it.

## 4. Numbers the repository can use

**Registry mapping.** `maltol` is keyed (`data/keys/compounds.yml` line 949). `vanillin` is keyed.
Neither has a `COMPOUND_STRUCTURE` entry in `parameters_matrix.py` — vanillin would be an aromatic
aldehyde (a **phenolic** one, which is a class the structure table does not carry; the nearest
existing member is `4_ethylphenol`, class `phenol`) and maltol a **hydroxypyranone** (the nearest
relative in the module is `furaneol`, which appears as a `REVERSIBLE_BINDING` compound but likewise
has no structure entry). **`MATRIX_LOADING` has no `soymilk` entry**; the nearest is
`soy_paste_hong`. A soymilk medium at ~30 g/L would need creating, with the loading marked as an
approximation offered in the Discussion rather than a measured composition.

Every row below shares: **soymilk from IA 3027 soybeans, ~3 % protein, held 2 h at 44 C (with or
without protein-glutaminase at 40 U per g protein); 1 mL aqueous flavour stock into 14 mL soymilk in
a 125 mL Teflon sniff bottle; equilibrated 24 h at 4 C; conditioned 1 h at 25 C before evaluation;
ASTM E697-04, ORTHONASAL, DETECTION, 3-AFC ascending, group BET = geometric mean of individual
BETs.**

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **orthonasal odour DETECTION threshold, vanillin, in CONTROL soymilk** | **9.61 ± 1.39** | ug/mL (= ppm = 9610 ug/L) | ~3 % soy protein, pH ~6.44, 25 C at evaluation, **21 panelists**, 3-AFC, ASTM E697-04 | Table 4, p. C4 | **threshold `[M]`** — medium = **soymilk**, NOT water |
| **orthonasal odour DETECTION threshold, vanillin, in DEAMIDATED soymilk** | **1.80 ± 2.44** | ug/mL | as above, 66.4 % DD, **21 panelists** | Table 4, p. C4 | **threshold `[M]`** — **SD is 136 % of the mean (mine); the least reliable cell in the paper** |
| **orthonasal odour DETECTION threshold, maltol, in CONTROL soymilk** | **23.8 ± 2.61** | ug/mL (= 23 800 ug/L) | ~3 % soy protein, **27 panelists**, 3-AFC | Table 4, p. C4 | **threshold `[M]`** — **the corpus's first measured maltol threshold in a plant-protein beverage** |
| **orthonasal odour DETECTION threshold, maltol, in DEAMIDATED soymilk** | **7.01 ± 3.71** (Table 4; the text says 7.02) | ug/mL | as above, 66.4 % DD, **27 panelists** | Table 4, p. C4 | **threshold `[M]`** — SD 53 % of mean (mine) |
| **threshold ratio, control soymilk / deamidated soymilk, vanillin** | **5.34** | — | same panel, same base soymilk, same protocol; **only the glutamine chemistry differs** | 9.61/1.80 **(mine)**; the paper says "approximately 5-fold" | **within_study_ratio** — but the one-SD interval on the ratio is **2.3x to unbounded (mine)** |
| **threshold ratio, control / deamidated, maltol** | **3.40** | — | as above | 23.8/7.01 **(mine)**; the paper says "approximately 3 times lower" | **within_study_ratio** — one-SD interval **2.2x to 7.2x (mine)**; the more defensible of the two |
| **paired WATER threshold for either compound** | — | — | — | **DOES NOT EXIST IN THIS PAPER** | **absent** — see §1. No `K_g`, no matrix-shift factor and no matrix-corrected threshold can be built |
| vanillin BET in skim milk | **7.41** | ppm | skim-milk model mixtures, another panel, 2004 | Results p. C4, **CITED** from Karagül-Yüceer et al. 2004, J Sens Stud 19:1-3 | **threshold `[C]`** — matrix-to-matrix comparator only |
| maltol BET in skim milk | **16.6** | ppm | as above | Results p. C5, **CITED** | **threshold `[C]`** |
| soymilk vs skim milk threshold ratio, vanillin / maltol | **1.30 / 1.43** | — | cross-study, cross-matrix, cross-panel | **(mine)** | **derived_assumption** — two protein matrices at nominally the same 3 % protein, 1.3-1.4x apart in the same direction |
| **degree of deamidation** | **66.4 ± 2.6** | % of glutamine residues | 44 C, 40 U/g protein, 2 h, PG "Amano" 500 | Table 2, p. C3, n = 3 | **measured_ratio** — the treatment dose |
| **degree of hydrolysis** | **4.25 ± 0.42** | % | as above | Table 2, p. C3, n = 3 | **measured_ratio** — the **confound**: 4.25 % of peptide bonds were also cleaved (Flags 4) |
| protein solubility, DSM / CSM / FSM at pH 5.0 | **~35 / ~22-23 / ~22-23** | % | 0.1 M acetate-phosphate, 25 C overnight, Lowry | Results p. C4, Fig. 1 | **level_only** — approximations read off a bar chart by the authors |
| protein solubility at pH 7.0 and pH 3.0 | **~100 %** (all three) and **~24 %** (all three) | % | as above | Results p. C4, Fig. 1 | **level_only** |
| pH of untreated / control soymilk | **6.74 / 6.44** | — | as-made | Results p. C4 | **level_only** — **DSM's pH is never stated** |
| soymilk protein loading | **~3** | % w/v (≈ 30 g/L) | stated in passing as the reason soymilk and skim milk agree | Results p. C4 | **derived_assumption** — an approximation in the Discussion, **not a measured composition of this batch** |
| **aroma of CSM and DSM does NOT differ** | p < 0.4443, estimated probability 0.4902, **power 3.7 %** | — | 2-AFC warm-up, 17 panelists, triplicate, beta-binomial | Results p. C4 | **measured_bound** — but at **3.7 % power this is not evidence of no difference** (Flags 10) |
| **taste of CSM and DSM DOES differ** | p < 0.0116, estimated probability 0.7059, **power 75.4 %**, gamma 0.3389 | — | as above | Results p. C4 | **measured_ratio** — the authors themselves attribute it possibly to viscosity/mouthfeel rather than taste |
| Stevens exponent n, vanillin, CSM / DSM | **0.4777 / 0.5525** | — | 12 panelists, 15-point scale, 5-point ladder, **linear section only** | Table 5, p. C6 + Fig. 2B | **derived_assumption** — a fitted slope with **no stated uncertainty and no n for the regression** |
| Stevens exponent n, maltol, CSM / DSM | **0.5151 / 0.5702** | — | as above | Table 5, p. C6 + Fig. 3B | **derived_assumption** |
| Stevens intercepts (printed only inside Figs. 2B/3B) | vanillin **−0.2119** (CSM) / **−0.2076** (DSM); maltol **−0.2476** (CSM) / **−0.2361** (DSM) | log10 units | as above | Figs. 2B, 3B | **derived_assumption** |
| Fechner intensity at 50 ug/mL vanillin, DSM / CSM | **4.4 / 3.0** (ratio **1.47x, mine**) | 15-point scale | 12 panelists | Results p. C6, Fig. 2a | **derived_assumption** — **interpolated off a figure at a concentration never presented** (Flags 7) |
| n·K, vanillin, SPI vs DSPI at 25 C | **88.8 x 10^4** vs **9.69 x 10^4** (**9.16x, mine**) | M^-1 | aqueous model, isolate not soymilk | Results p. C5, **CITED** from Suppavorasatit & Cadwallader 2012, JAFC 60:7817-23 | **binding_constant `[C]`** — **not measured here; a different substrate, a different DD, and an instrumental method, so it must NOT be pooled with the sensory thresholds** |
| n·K, maltol, SPI vs DSPI at 25 C | **303 x 10^4** vs **79.6 x 10^4** (**3.81x, mine**) | M^-1 | as above | Results p. C5-C6, **CITED** | **binding_constant `[C]`** |
| **cross-method agreement on maltol** | binding ratio **3.81x** (instrumental, SPI, cited) against threshold ratio **3.40x** (sensory, soymilk, measured here) = **1.12x apart (mine)**; vanillin's equivalents are **9.16x** vs **5.34x** = **1.72x apart (mine)** | — | cross-paper, cross-substrate, cross-DD | **(mine)** | **derived_assumption** — a magnitude coincidence, recorded not shipped |
| stimulus ladders | Table 1 above, six sets per compound per soymilk | ug/mL | design values | Table 1, p. C3 | **level_only** — bounds on where each BET could have landed |
| every solubility bar, every dose-response point | — | — | — | Figs. 1, 2A, 3A | **figure_only** |

### Can these be put on the same basis as what the repository already carries?

**(a) The four thresholds: yes, as absolute `threshold` records with `medium = soymilk`.** They are
orthonasal detection thresholds by a named ASTM protocol on 21 and 27 panelists, which is a better
provenance than most threshold entries in the corpus. They must carry `medium`, `route =
orthonasal`, `criterion = detection`, `protocol = ASTM E697-04`, `n_panelists`, and the SD.

**(b) The CSM/DSM ratio: yes, as a `within_study_ratio`, and it is a genuinely novel kind of
observation for this repository** — a threshold shift produced by changing the protein's covalent
chemistry with everything else held fixed. But note what it does and does not license. It says
**"a 66.4 % conversion of glutamine to glutamate in a 3 % soy beverage lowers the vanillin threshold
5.3x and the maltol threshold 3.4x"**. It does **not** say what fraction of the *absolute* threshold
elevation in ordinary soymilk is due to protein at all, because there is no protein-free leg.

**(c) The absolute threshold elevation caused by soy protein: NOT MEASURABLE from this paper.** That
is the quantity the matrix layer wants and it is exactly the quantity the design omits. **If a water
leg existed the paper would unlock a `MATRIX_THRESHOLD` row for two compounds; it does not exist.**

**(d) Nothing here is a binding constant measured by these authors.** The two n·K pairs are cited
from the companion JAFC paper on **soy protein ISOLATE in an aqueous model**, not on soymilk. If the
repository wants a soy vanillin/maltol binding constant, **the source to fetch is Suppavorasatit &
Cadwallader 2012, JAFC 60:7817-7823**, not this article.

**(e) Nothing here is a rate, an activation energy or a partition coefficient.** The only
temperature that touches the flavour is 4 C (equilibration) and 25 C (serving). The 44 C is the
enzyme's, applied before flavouring.

## 5. Flags

1. **The two deamidated-soymilk thresholds are extremely noisy and one of them is unusable as a
   point.** Vanillin/DSM is **1.80 ± 2.44 ug/mL** — the standard deviation is **136 % of the mean
   (mine)**, so the one-SD interval crosses zero. Maltol/DSM is **7.01 ± 3.71**, SD **53 % of the
   mean**. The two control legs are at 14 % and 11 %. **The paper reports "approximately 5-fold" for
   vanillin without acknowledging that its own SD makes anything from ~2.3x upward consistent with
   the data (mine).** If only one threshold shift is carried forward, carry **maltol's**.
2. **The vanillin thresholds in CSM and DSM were measured on DIFFERENT concentration ladders.**
   Table 1: CSM runs 1.667-405.0 ug/mL (anchor x = 45.00), DSM runs 0.556-135.0 (anchor x = 15.00) —
   **the whole DSM ladder is 3x lower**. Maltol's two ladders are identical. A 3-AFC
   ascending-series BET depends on where the ladder starts and how many steps precede the detection
   point; two ladders offset by 3x are not the same instrument. **The maltol comparison is
   design-matched and the vanillin comparison is not**, which is a second, independent reason to
   prefer the maltol shift. The paper does not remark on this.
3. **The stated mechanism is chemically wrong and must not be imported.** The paper attributes the
   effect to lost **Schiff-base** capacity: *"deamidation by PG reduces the number of available
   amide groups by converting most of glutamine residues in soy protein to glutamic acids.
   Therefore, the potential to form covalent bonds (Schiff bases) with the carbonyl group of
   vanillin was decreased."* **A Schiff base forms between an aldehyde and a primary AMINE — the
   lysine epsilon-amino group — not with a glutamine side-chain AMIDE, which is not nucleophilic
   enough and is not what protein-glutaminase converts anything into.** `matrix_sites.py`'s
   `BINDING_CLASSES` is explicitly an aldehyde-to-lysine-amine channel and lysine is untouched by
   PG. Whatever moved these thresholds, **it was not the loss of a Schiff-base partner**; the more
   plausible routes are the charge/pI change, the conformational change, the 4.25 % hydrolysis, or
   the accompanying pH drop. **The observation stands; the explanation should be recorded as the
   paper's claim and flagged, not adopted.**
4. **Deamidation is confounded with 4.25 % hydrolysis.** Table 2: DH = 4.25 ± 0.42 %. The paper
   concedes the commercial enzyme *"could contain some residual protease activity"* and that peptide
   fragments may be released. So the DSM leg differs from CSM in **at least four ways**: 66.4 % less
   glutamine, correspondingly more glutamate, 4.25 % of peptide bonds cleaved, and (necessarily) a
   lower pH. **The threshold shift cannot be assigned to any one of them.**
5. **The pH of DSM is never stated.** FSM 6.74 and CSM 6.44 are printed; DSM's is not. Converting
   66 % of a protein's glutamine to glutamate releases ammonia and creates carboxylates and **must**
   change the pH of an unbuffered beverage. The corpus's adduct chemistry is pH-gated
   (`PH_ADDUCT_GATE`: carbonyl-lysine adduct formation is abolished at pH 3; Leksrisompong's
   diacetyl binds caseinate at pH 7 and not at pH 5.5), so **an unstated pH difference between the
   two legs is a live alternative explanation for the entire result.**
6. **The blanks are flavoured soymilk's own background, and the panel says that background overlaps
   the target odours.** Table 3's panel-generated aroma terms for plain soymilk include **"Burnt
   sugar"** and **"Caramel"** — which is precisely the maltol percept ("cotton candy"). A detection
   threshold measured against a blank that already smells faintly of the target is elevated by an
   unknown amount. The vanillin case is cleaner ("Vanilla" appears in Table 3, but as a term
   distinguishing the two soymilks, and the two soymilks did **not** differ in aroma).
7. **The Fechner comparison (4.4 vs 3.0 at 50 ug/mL) is an interpolation, not a measurement.**
   50 ug/mL is not one of the five concentrations presented to either panel — CSM got 1000, 200, 40,
   8, 1.6 and DSM got 500, 100, 20, 4, 0.8. **Both intensity values were read off Fig. 2a between
   points**, and no SD is given for either. It is quoted in the Results as though it were a result.
8. **The paper's significance notation is printed backwards in at least three places.** Methods:
   *"to determine significant differences among treatments (P > 0.05)"*. Figure 1 caption: *"same
   upper case letters within the same sample across different pH values indicate a significant
   difference (P > 0.05)"* — while the same caption's first half correctly uses *"not significantly
   different (P > 0.05)"*. Results, p. C5: the SPI/DSPI n·K difference is called *"9-times greater
   (P > 0.05)"*, which as printed says the 9-fold difference was **not significant**. **Every
   `P > 0.05` in this paper should be read as probably meaning `P < 0.05`, and none of them can be
   relied on as printed.**
9. **Table 4 and the running text disagree on maltol/DSM: 7.01 vs 7.02 ug/mL.** Trivial in
   magnitude (0.14 %) but it means a quoted maltol shift is 3.40x or 3.39x depending on the line
   taken. Table 4 is the table and should win.
10. **The "no aroma difference" result has 3.7 % power and is therefore not evidence of no
    difference.** The paper reports it as one (*"This result indicated that deamidation did not
    affect the overall aroma of the soymilk"*). At 3.7 % power with 17 panelists, a real and sizeable
    aroma difference would have been missed 96 % of the time. **This matters because the threshold
    experiment's blanks depend on CSM and DSM being aromatically equivalent**, and that equivalence
    is asserted on an essentially powerless test.
11. **Maltol's chemical name is printed wrongly.** *"maltol (3-hydroxy-2-methoxyl-4H-pyran-4-one)"*
    — maltol is 3-hydroxy-2-**methyl**-4H-pyran-4-one. The 2-methoxy compound is a different
    substance. The supplier (FONA), the purity (>98 %) and the descriptor ("cotton candy") all
    identify maltol correctly, so this is a slip in the Materials section, but anyone keying the
    compound from this paper's text would key the wrong molecule.
12. **Year ambiguity.** Copyright 2012, submitted and accepted 2012, DOI minted with a 2012 stem
    (`10.1111/j.1750-3841.2012.03012.x`), **issue and pagination 78(1) 2013, pp. C1-C7**. The file on
    disk is `suppavorasatit2012.pdf`. Cite as **J Food Sci 78(1):C1-C7 (2013)** with the 2012 DOI, or
    the citation will not resolve against the volume.
13. **The protein loading is an approximation from the Discussion, not a measurement.**
    "approximately 3% protein" is offered as an explanation for the soymilk/skim-milk agreement, not
    as a composition. **No total-protein figure for FSM, CSM or DSM is printed** — the Lowry assay
    appears only as the solubility denominator. Any per-gram construction from this paper divides by
    a number the paper never measured.
14. **Two different panels, two different sizes, for the two compounds.** Vanillin: 21 panelists
    (5 M, 16 F, 22-48 y). Maltol: 27 panelists (5 M, 22 F, 19-40 y). Aroma intensity: a third panel
    of 12. Difference testing: a fourth panel of 17. **The vanillin and maltol thresholds are not
    from the same panel** and should not be compared to each other as a within-panel ratio — only
    each compound's own CSM/DSM pair is within-panel.
15. **The 24 h equilibration was at 4 C and the evaluation at 25 C.** Flavour-protein equilibria are
    temperature-dependent, and the samples spent 1 h warming before being sniffed. Whatever
    equilibrium the 24 h at 4 C established, **the panel smelled a system 1 h into a 21 K
    re-equilibration**. The paper's own justification for the cold hold is "to allow for
    equilibration of flavor-matrix interactions", which the warming step then partly undoes.
16. **What this paper does NOT contain**: any threshold in water; any protein-free leg of any kind;
    any binding constant measured by these authors; any partition coefficient; any headspace
    measurement; any rate; any activation energy; any temperature above 44 C; any compound other than
    vanillin and maltol; any pH for DSM; any measured protein content; any error bar on a Stevens
    exponent; any of the dose-response points as numbers.
17. **What to request from the authors**: (i) **the water thresholds for vanillin and maltol from
    the same panel** — this single addition would convert the paper into the paired water/matrix
    record the matrix layer needs, and the panel and protocol already exist; (ii) the pH of DSM;
    (iii) the measured total protein of FSM, CSM and DSM; (iv) the individual BETs behind the four
    group BETs, so the geometric-mean statistic and its dispersion can be handled properly;
    (v) whether the ±SD in Table 4 is the SD of the individual BETs on a linear or a log scale;
    (vi) why the vanillin ladders differ between CSM and DSM while the maltol ladders do not;
    (vii) the raw intensity ratings behind Figs. 2 and 3, with SDs and the number of points used in
    each "linear section" fit.
