# Lee 2024 — EXTRACTION (the VOLATILE half of the same model-cake experiment as Lee 2022: the same glucose and glucose + leucine cakes, the same oven at 140/170/200 C and 25/50 Hz, with ten volatile markers trapped ON-LINE from the baking vapour on sorbent tubes and quantified by TD-GC-MS against two deuterated internal standards, plus a dynamic matrix-to-vapour transport model for furfural and HMF)

### THE PAPER THAT SHOWS NITROGEN LEAVING THE CAKE: **3-methylbutanal and all four pyrazines are found ONLY in the glucose + leucine cake and never in the glucose-only cake**, and 3-methylbutanal is the most abundant marker of the ten — so leucine's Strecker aldehyde, five of its six carbons, walks out of the matrix into the oven air. Read against Fang 2009, which finds no pyrazine nitrogen inside the polymer, this locates a nitrogen sink that competes with `MEL_N` and is invisible to any measurement made on the cake.

**Source on disk:** `data/articles/lee2024.pdf` (10 pp. including the HAL cover; the article is
Food Research International 2024, **183**, 114183, **open access CC BY-NC**).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/lee2024.txt`), which is the **typeset Elsevier PDF** and came through
cleanly except for subscripted units, which reflow badly; the unit strings in section 3.2.1 were
therefore confirmed by **rendering page 4 at 200 dpi** and are transcribed exactly as printed.
**Tables 1 and 2 came through the text layer intact** and are re-typed in full below. Figures 1,
2, 4 and the supplementary figures carry every concentration-time datum and are **figure_only**.
**Figure 4 panel C is a small parameter grid printed inside the figure**, holding the fitted
transport constants; per the house rule it is treated as figure_only and its values are **not**
typed — see section 3 and Flags 2 for what it contains and why that is the right call here.
**Supplementary material (Figs. S1-S4) exists at the article DOI and is NOT on disk.** Repo status
before this dossier: `lee2024.pdf` has **no extraction dossier** and is not cited anywhere in
`src/`.

**On the file names.** `lee2022.pdf` and `lee2024.pdf` share the first two lines of their title and
were checked against each other before either dossier was written. **They are two different
papers**: same authors, same model cake, same oven, same three set-points, but different journals,
different years and disjoint measured species. `lee2022.pdf` (Food Chemistry 376:131917) measures
the **non-volatile** precursors, alpha-dicarbonyls, furanic compounds and browning **in the
matrix**; this one measures **ten volatiles in the oven vapour**. This paper cites that one as
"Lee et al. (2022)" throughout and depends on its matrix concentrations for its transport model.
Both dossiers exist and should be read together — see `lee2022_extraction.md`.

## 0. Identity

| field | value |
|---|---|
| Title | "Unravelling caramelization and Maillard reactions in glucose and glucose + leucine model cakes: Formation and degradation kinetics of **volatile markers extracted during baking**" |
| Authors | J. (Jeehyun) Lee (Université Paris-Saclay, INRAE, AgroParisTech, UMR SayFood, Palaiseau; and INRAE, Institut Agro, STLO, Rennes); S. (Stéphanie) Roux; D. (listed as D., e-mail nicolas.) Descharles; B. (Barbara) Rega; **C. (Catherine) Bonazzi — corresponding** |
| Venue | Food Research International 2024, 183, 114183. Received 5 December 2023, revised 26 February 2024, accepted 28 February 2024, online 2 March 2024. **Open access, CC BY-NC.** Preprint: HAL hal-05318480 |
| DOI | 10.1016/j.foodres.2024.114183 |
| Naming | **G** = glucose-only cake; **G + L** = glucose + leucine. Markers: **3-MB** (3-methylbutanal), **P** (pyrazine), **2-MP** (2-methylpyrazine), **2,5-DMP**, **2,6-DMP**, **AA** (acetic acid), **F** (furfural), **5-MF** (5-methylfurfural — but abbreviated "5-MP" once in the reagent list), **FA** (furfuryl alcohol), **5-HMF**. **dIS** = deuterated internal standards |
| Lineage | the SayFood model-cake programme continued: **Lee 2020** (the cake, and the 18 volatiles identified qualitatively), **Lee 2021** (the TD-GC-MS quantitative method with isotope standard addition — the method paper this one depends on entirely), **Lee 2022** (the matrix half), Fehaili 2010 (the instrumented oven), Srivastava 2018, Rega 2009, Cepeda-Vázquez 2018. Chemistry cited to Knol 2010 and Martins & van Boekel 2005 (acetic acid from 1-deoxyosone), **Kocadagli & Gokmen 2016** (near-zero 1-deoxyosone formation direct from fructose), Nikolov & Yaylayan 2011 and Sanders 2003 (5-MF from HMF), Wnorowski & Yaylayan 2000, Yaylayan & Keyhani 2000 and Delatour 2020 (three competing routes to furfuryl alcohol) |
| Companion on disk | **`lee2022.pdf`** — dossier `lee2022_extraction.md` |
| Companions on disk | `fang2009_extraction.md` (**the paper that says pyrazines are NOT in the melanoidin** — read against this one), `knol2010_extraction.md` (cited here for acetic acid), `martins2005_extraction.md`, `kocadagli2016foodchem_extraction.md`, `mundt2004_extraction.md` (the glycine C/N this paper's leucine arithmetic must be contrasted with) |

## 1. Why it matters

**(a) It measures a nitrogen sink that leaves the food, and the cluster needs it.** The trunk
books unmeasured amine nitrogen into `MEL_N` (`src/kinetic_core/species.py`). This paper's
central, cleanest result is that **3-MB, P, 2-MP, 2,5-DMP and 2,6-DMP were "only found during the
baking of model cake G + L"** — they are absent from the glucose-only cake by construction,
because there is no nitrogen in it. All five are **nitrogen-derived volatiles trapped from the
oven air**, i.e. they have left the matrix. Four of them (the pyrazines) carry **two nitrogens
each**. So the amine nitrogen in a real baked matrix has at least three fates: the polymer, the
remaining free amino acid, and **the oven**. The trunk models the first two and not the third.

**(b) Read against Fang 2009 the point sharpens into a mechanism.** Fang & Schmidt-Rohr looked
for pyrazine nitrogen inside a glucose-glycine melanoidin by 15N NMR and found **none** — "no
significant 15N NMR signals of pyridines or pyrazines (between 250 and 350 ppm) ... are observed
in our or any published spectra". Lee finds pyrazines in the vapour. **Together: pyrazine
nitrogen is a volatile loss, not polymer nitrogen.** That is a genuine `MEL_N` leak, it is
quantified here (in the vapour, in the paper's own units), and nothing in the repository accounts
for it.

**(c) It corrects the arithmetic in the `lee2022_extraction.md` dossier on what a leucine
melanoidin's C/N would be — and the correction points the other way.** That dossier notes that
leucine is C6N1 against glycine's C2N1, so under the trunk's step-9 rule
(`MELANOIDIN_REPEAT_UNIT_CARBON = 8` = 6 from 3-deoxyglucosone + 2 from an intact glycine) a
leucine repeat unit would be **12 C per N**. That is the **upper** bound and it assumes the amine
arrives intact. This paper measures the other end. In Strecker degradation the amino acid gives up
its carboxyl as CO2 **and** its side chain plus alpha carbon as the Strecker aldehyde, while its
**nitrogen transfers to the dicarbonyl** as an alpha-aminocarbonyl. For leucine the aldehyde is
**3-methylbutanal, C5** — and this paper finds it is "**the most abundant of all the reaction
markers**", formed "very rapidly in large quantities, even under the mildest baking conditions",
and at 200 C above the calibration range. **Five of leucine's six carbons walk out of the cake as
vapour.** A fully Strecker-routed leucine therefore contributes **zero carbon and one nitrogen**
to the polymer, and the repeat unit becomes 6 + 0 = **6 C per N**. So the honest bracket for a
leucine melanoidin is **6 to 12**, not 12, and the measurement in this paper argues for the low
end. See section 3 arithmetic 4 and Flags 8 for the contrast with glycine, whose Strecker aldehyde
is **formaldehyde** — which Mundt & Wedzicha's radiochemistry shows is *retained* in the polymer,
not lost.

**(d) It is the same matrix, oven and formulations as Lee 2022, so the two are one dataset.**
Every conditioning caveat in `lee2022_extraction.md` applies here unchanged: oven set-point is not
product temperature, the leucine concentration is never printed, the product is a 1.3 cm cake with
a measured internal gradient. What is added here is the **vapour compartment**, and a dynamic
model linking the two.

**(e) It replicates and strengthens a null the repository can use.** Lee 2022 found no effect of
fan frequency on the temperature or water profiles; this paper tests the same variable on **ten
chemical markers by ANOVA** and finds p > 0.05 for the convection level on **eight of ten**
(the exceptions, furfural and furfuryl alcohol, are argued away as low-concentration noise near
the LOQ). That is a properly tested null on an oven variable, and it is what licenses pooling the
two fan settings in both papers.

**(f) It supplies a topology constraint on acetic acid.** The trunk carries `AA` (acetic acid) as
a measured B1 product. This paper places it "directly formed via the degradation of
1-deoxyosone", cites Knol 2010 and Martins & van Boekel 2005 for that route, and adds that it
"made it possible to highlight **another route of 1-deoxyosone degradation** other than those
leading to the formation of glyoxal, diacetyl and methylglyoxal". It also measures acetic acid
**roughly 10-fold higher in G + L than in G** — the largest formulation effect of any furanic or
acid marker in the paper.

What this paper does NOT give the repository: any chemical rate constant (the only fitted
constants are transport constants, Flags 2), any elemental analysis, any C/N, any measurement in
the solid matrix (that is the 2022 companion), any browning measurement, and any tabulated
concentration.

## 2. Methods as they matter to a model

Everything about the matrix, the formulations and the oven is identical to Lee 2022 and is set
out in `lee2022_extraction.md` section 2. What follows is what is new or restated here.

- **The cake and the oven, restated.** Ultrapure water, native corn starch (**12.4 % w/w water**,
  measured per **NF norm V05 707** — this paper cites the norm, the 2022 one does not),
  methylcellulose SGA7C and HPMC K250M (Dow). Glucose (Roquette) or glucose + L-leucine (Sigma,
  >= 99 %, food grade). **20 g of batter per aluminium mould, diameter 6.6 cm, mould height 4 cm,
  cake height 1.3 cm.** Instrumented pilot oven (Bongard), **volume printed here for the first
  time: V_oven = 96 L**. Set-points **140, 170, 200 C**; fan **25 and 50 Hz**.
- **Replication, which differs from the 2022 paper.** "The baking trials at **140 C/25 Hz,
  170 C/25 Hz, 170 C/50 Hz and 200 C/25 Hz** were conducted in triplicate, while the other
  combinations were performed once only." Figure captions give **n = 4 at 140 C, n = 6 at 170 C,
  n = 4 at 200 C** after pooling the two fan settings.
- **On-line vapour sampling — the method that defines this paper.** Vapour was drawn during
  baking through an **Air Toxics sorbent tube** (PerkinElmer) at an extraction flow rate of
  **50 mL/min**, per Lee 2021. Trapping intervals: **2-6, 8-12, 20-24, 33-37, 52-56 and
  86-90 min**, with the last replaced at 140 C by **91-95 and 116-120 min**. Each interval is
  assigned its **median time: 4, 10, 22, 35, 54 and 88 min** (93 and 118 at 140 C). **So every
  vapour datum is a 4-minute integral reported at its midpoint, not an instantaneous
  concentration** (Flags 3).
- **Internal standards, added after sampling.** Each tube was removed and **spiked with 1 uL of
  the two dIS (d4-pyrazine and d4-furfural, both at 8.00e-2 g/L)**, capped with Teflon, stored at
  room temperature and **analysed the same day**. Calibration solutions held a constant dIS
  concentration against increasing analyte concentrations, per Lee 2021. Note the standards
  correct for desorption and for the GC-MS, **not** for the trapping efficiency of the sorbent
  tube during baking.
- **Statistics.** The area ratios (marker / dIS) for the ten markers at 170 C, both fan settings,
  both formulations, gave a matrix of **10 variables x 72 individuals** (6 time intervals per
  trial). Linear model by `Linearmodel` (FactoMineR), one-way ANOVA (LSD), significance by
  "Turkey's test" (i.e. **Tukey**, printed as Turkey) at **p <= 0.05**, in R.
- **The transport model, which is the paper's second contribution.** The oven is treated as an
  **open, perfectly stirred reactor** of volume V_oven through which an air flow V̇_air
  circulates, with no volatile present at t = 0 and none entering with the air. The mass balance
  reduces to
  `(dC_vj/dt) = −C_vj·(V̇_air/V_oven) + S_j(t)`, with the source term split into an appearance and
  a disappearance term,
  `S_j(t) = k_a(T)·(dC_mw/dt)·C_mj(t) − k_d(T)·C_vj(t)`,
  where `C_mw` is the **water concentration in the matrix** — so the appearance term is driven by
  the **evaporation rate**, not by a partition coefficient. `k_a` is stated to lump three things:
  a matrix/vapour partition, the water evaporation flow rate, and a first-order kinetic term.
  Both constants are given a **reparameterised Arrhenius form**,
  `k(T) = k_ref·exp(−(Ea/R)(1/T − 1/T_ref))`. Water content was fitted as an exponential decay
  `C_mw = C_mw0·exp(−alpha·t)` and each marker's matrix concentration as a **sigmoid**
  `C_mj = a / (1 + b·exp(−c·t))`, both purely as smooth interpolants so the ODE solver had
  continuous inputs. Solved in MATLAB with `ode15s` and `fminsearch`, minimising a weighted sum of
  squares "designed to give the same weight to each condition regardless of the number of
  experimental points and the levels of concentration". **Run with V̇_air = 0.05 L/min and
  V_oven = 96 L**, on **furfural and 5-HMF only**, with `k_a` allowed to differ between G and
  G + L while `k_d` was held common.
- **Units, and a discrepancy.** Section 3.2.1 reports vapour concentrations "expressed in **mmol
  per L in the oven and per g of dry matter (DM) of product**", written `mmol·L_air^-1·g_DM^-1`
  (confirmed from the page render). Section 3.4 says that "the values in vapor previously shown as
  **nmol·L_air^-1·g_DM^-1** were multiplied by 10^-3 and V_oven (L)" to give the
  `umol·g_DM^-1` used in Figure 4. **Those two statements differ by a factor of 10^6** and cannot
  both be right (Flags 1).
- **What is not measured here.** Nothing in the solid matrix — no glucose, no leucine, no
  dicarbonyl, no browning (all of those are Lee 2022, and this paper takes its matrix
  concentrations from that dataset). No temperature inside the cake is reported here either; the
  paper says "the temperatures at different positions of model cakes ... were however measured and
  presented in a previous study (Lee et al., 2022)".

## 3. Tables re-typed

### Table 1. "Characteristics of the 10 volatile markers and 2 deuterated internal standards."

| compound | CAS No. | log Kow ᵃ | Pv (mm Hg at 25 C) ᵇ | selected ions (m/z) ᵈ |
|---|---|---:|---:|---|
| **Analytes (A)** | | | | |
| 3-methylbutanal | 590-86-3 | 1.23 | 50 | 58, 71 |
| Pyrazine | 290-37-9 | −0.2 | 10.81 | 53, 80 |
| 2-methylpyrazine | 109-08-0 | 0.21 | 8.06 | 67, 94 |
| 2,5-dimethylpyrazine | 123-32-0 | 0.63 | 3.18 | 42, 108 |
| 2,6-dimethylpyrazine (printed "2,6–imethylpyrazine") | 108-50-9 | 0.54 | 4.57 | 42, 108 |
| Acetic acid | 64-19-7 | −0.17 | 15.7 | 60 |
| Furfural | 98-01-1 | 0.41 | 2.21 | 95, 96 |
| 5-methylfurfural | 620-02-0 | 0.67 | 0.610 | 81, 109, 110 |
| Furfuryl alcohol | 98-00-0 | 0.28 | 0.609 | 81, 97, 98 |
| 5-hydroxymethylfurfural | 67-47-0 | −0.6 | **5.28 x 10^-3** | 97, 126 |
| **Deuterated internal standards (dIS)** | | | | |
| d4-pyrazine | 1758-62-9 | −0.2 | nd ᶜ | 56, 84 |
| d4-furfural | 1219803-80-1 | 0.4 | nd ᶜ | 99, 100 |

Footnotes as printed: **ᵃ** octanol/water partition coefficients from PubChem or EPI Suite v4.0
databases; **ᵇ** Pv: vapour pressure from the PubChem database; **ᶜ** not determined; **ᵈ**
quantifier ions are in bold. **The bold formatting that marks the quantifier ion did not survive
the text extraction**, so for the multi-ion rows it is not recoverable from this file which ion is
the quantifier (Flags 6).

Note the vapour-pressure span this table encodes: **3-methylbutanal at 50 mm Hg down to 5-HMF at
0.00528 mm Hg — a factor of ~9500**. That range is the whole reason the paper needs a transport
model rather than a partition coefficient.

### Table 2. "Results of the F-test according to the ANOVA model for 10 volatile markers."

Significant values are in bold in the original; bold is preserved below where the p-value is
<= 0.05.

| marker | R² | Time | Formula | Convection level | Time × Formula | Time × Convection level | Formula × Convection level |
|---|---:|---|---|---|---|---|---|
| 3-methylbutanal | 0.921 | **6.80·10⁻¹⁰** | **< 2.20·10⁻¹⁶** | 0.093 | **6.80·10⁻¹⁰** | **0.002** | 0.093 |
| Pyrazine | 0.490 | 0.186 | **4.58·10⁻⁰⁶** | 0.125 | 0.186 | 0.513 | 0.125 |
| 2-methylpyrazine | 0.994 | **< 2.00·10⁻¹⁶** | **< 2.00·10⁻¹⁶** | 0.070 | **< 2.00·10⁻¹⁶** | **0.039** | 0.070 |
| 2,5-dimethylpyrazine | 0.987 | **< 2.00·10⁻¹⁶** | **< 2.00·10⁻¹⁶** | 0.547 | **< 2.00·10⁻¹⁶** | 0.823 | 0.547 |
| 2,6-dimethylpyrazine (printed "2,6–imethylpyrazine") | 0.992 | **< 2.00·10⁻¹⁶** | **< 2.00·10⁻¹⁶** | 0.323 | **< 2.00·10⁻¹⁶** | 0.462 | 0.323 |
| Acetic acid | 0.894 | **< 2.20·10⁻¹⁶** | **2.39·10⁻¹²** | 0.101 | **3.78·10⁻¹⁴** | 0.205 | 0.060 |
| Furfural | 0.964 | **< 2.20·10⁻¹⁶** | **< 2.20·10⁻¹⁶** | **3.42·10⁻⁰⁷** | **3.16·10⁻¹³** | **3.34·10⁻⁰⁵** | 0.094 |
| 5-methylfurfural | 0.796 | **2.80·10⁻⁰⁹** | **1.02·10⁻⁰⁹** | 0.160 | **8.89·10⁻⁰⁹** | 0.606 | 0.943 |
| Furfuryl alcohol | 0.962 | **< 2.20·10⁻¹⁶** | **< 2.20·10⁻¹⁶** | **0.003** | **< 2.20·10⁻¹⁶** | 0.854 | 0.468 |
| 5-hydroxymethylfurfural | 0.363 | **0.029** | 0.326 | 0.280 | 0.4187 | 0.132 | 0.532 |

Read across: **the Formula effect (G vs G + L) is significant for nine of the ten markers** — the
sole exception being 5-HMF (p = 0.326), which also has by far the worst fit (**R² = 0.363**). The
**Convection level is non-significant for eight of ten**, the exceptions being furfural
(3.42·10⁻⁷) and furfuryl alcohol (0.003), which the text argues away as low-level noise near the
LOQ. Pyrazine itself has a poor model (R² = 0.490) and **no significant Time effect** (p = 0.186)
while its Formula effect is strongly significant — consistent with a compound that is either
present (G + L) or absent (G) but whose level does not move cleanly with time.

### The PubChem identifiers printed in the "chemical compounds studied" list

| compound | PubChem CID as printed |
|---|---|
| Glucose | 107526 |
| Leucine | 6106 |
| 3-Methylbutanal | 11552 |
| Pyrazine | 9261 |
| 2-Methylpyrazine | 7976 |
| 2,5-Dimethylpyrazine | 31252 |
| 2,6-Dimethylpyrazine | 7938 |
| Acetic acid | **176** |
| Furfural (printed "Furfual") | 7362 |
| 5-Methylfurfural | **176** |
| Furfuryl alcohol | 7361 |
| 5-Hydroxymethylfurfural | 237332 |
| d4-Pyrazine | 164233030 |
| d4-Furfural | 101759994 |

**CID 176 is printed twice**, for acetic acid and for 5-methylfurfural. One is wrong; Table 1's CAS
numbers (64-19-7 and 620-02-0) are unambiguous and should be used instead. The same "Furfual"
typo appears in the 2022 companion.

### Figure 4 panel C — what it contains, and why its values are not typed here

Figure 4's caption reads: "Modeling of furfural (A) and 5-hydroxymethylfurfural (B)
concentrations in vapor as a function of matrix concentrations, water evaporation and time
(Eq. (4)) for the 3 baking temperatures ... and the two formulas (G and G + L), and **identified
parameters (C)**." Panel C is a small typeset grid holding, for **furfural** and for **5-HMF**:
a **k_a,ref and an Ea_a for each formulation separately (G and G + L)**, and a single **k_d,ref
and Ea_d shared across formulations**. The units are those defined in section 2 — `k_a` in
(umol j in vapour)·g_DM·(umol j in matrix)⁻¹·(umol water in matrix)⁻¹ and `k_d` in min⁻¹, with Ea
in J/mol.

**Per the house rule these values are not typed**, and in this case that is also the substantively
right call for three reasons: (i) they are printed inside a figure panel, not in a table; (ii)
**they are not chemistry** — `k_a` is explicitly a lump of a partition coefficient, an evaporation
flow rate and a first-order term, and `k_d` has no identified mechanism at all (the paper offers
"deposition of molecules on the oven walls, some further reactions occurring in the vapor or even
gas leaks"); and (iii) **the reference temperature T_ref is never printed**, so the Arrhenius
reparameterisation cannot be evaluated at any temperature (Flags 2). Nothing in panel C can enter
the kinetic core. It is nonetheless the only place fitted constants appear in this paper and is
request #2 in Flags 12.

### Every quantity printed in the running text

**Vapour levels — units transcribed exactly as printed, see Flags 1.**

| quantity | value | conditions | where |
|---|---|---|---|
| **3-MB, P, 2-MP, 2,5-DMP, 2,6-DMP** | **found ONLY in G + L; absent from G** | all three set-points | Results 3.2.2; Fig. 3 |
| 3-MB abundance | "**the most abundant of all the reaction markers**"; formed "very rapidly in large quantities, even under the mildest baking conditions"; at 200 C "the quantities were **above the linearity range**" | all conditions | Results 3.2.2; Fig. 2 (hatched points) |
| acetic acid, G + L vs G | **roughly 10-fold higher** in the presence of leucine | 200 C | Results 3.2.1 |
| furfural among the furanics | "**F was found in the highest quantities**"; levels "**comparable in both formulae**" apart from the 86 min point | 200 C | Results 3.2.1 |
| 5-MF and FA maxima | "up to **2.0** and **1.5 mmol·L_air⁻¹·g_DM⁻¹**, respectively" | 200 C | Results 3.2.1 |
| 5-MF and FA, G + L vs G | both **5-fold higher** in G + L | 200 C | Results 3.2.1 |
| FA in G alone | the **lowest** of all furanic compounds, "up to **0.15 mmol·L_air⁻¹·g_DM⁻¹**" | G, 200 C | Results 3.2.1 |
| FA increase with leucine | "almost a **5-fold increase** between G and G + L at 54 min, 200 C" | 200 C | Results 3.3 |
| **5-HMF, the one marker LOWER with leucine** | up to **1.0** (G + L) against **1.5 mmol·L_air⁻¹·g_DM⁻¹** (G) | 200 C | Results 3.2.1 |
| time shape of AA, F, 5-MF, FA, 5-HMF | increased throughout baking at 200 C **until 56 min**; the 86 min point was **not significantly different in G** but **lower in G + L** | 200 C | Results 3.2.1 |
| pyrazine relative levels | P similar to 2,5-DMP; **2,6-DMP 2-fold lower than 2,5-DMP**; **2-MP the lowest of all pyrazines** | G + L | Results 3.2.2 |
| pyrazine time shape | P increased over time; the other three were **lower at 86 min than at 56 min** | G + L | Results 3.2.2 |
| convection level | **no significant influence** on 8 of 10 markers (p > 0.05); the two exceptions dismissed as near-LOQ noise | 170 C | Results 3.1; Table 2 |
| temperature effect | "a **nonlinear** accelerating impact", **bell-shaped kinetics for most markers at 200 C**; "a typical trend of the Arrhenius law, regardless of the nature of markers and precursors" | all | Abstract; Results 3.2 |
| transport model inputs | **V̇_air = 0.05 L·min⁻¹**, **V_oven = 96 L** | — | Results 3.4; Methods 2.2 |
| k_a dependence on formulation | "**considerable dependence on the formulae** for both F and 5-HMF", attributed to "the probable difference in the **retention capacity of the matrix because of the higher molecular weight polymers (melanoidins) formed in the G + L formula which could trap the volatile compounds**" | — | Results 3.4 |
| k_a and k_d vs volatility | both **lower for 5-HMF than for furfural**, "consistent with the lower volatility of the molecule" | — | Results 3.4 |
| k_d magnitude | "**far from negligible** as this was the only way that the model could obtain vapor concentrations which were lower at the end of the kinetic, the impact of air renewal being too weak" | — | Results 3.4 |
| number of volatiles identified in this cake previously | **at least 18** (Lee 2020); **10** quantified here | — | Introduction |
| volatiles formed in baked cereal goods generally | **more than 540** (Cho & Peterson 2010) | — | Introduction (**not measured here**) |

**Mechanistic statements the paper prints:**

- **Acetic acid comes from 1-deoxyosone**, "a new route of 1-deoxyosone degradation other than
  those leading to the formation of glyoxal, diacetyl and methylglyoxal" (citing Knol 2010,
  Martins & van Boekel 2005).
- **1-Deoxyosone forms faster via the Amadori product than via caramelisation.** "The very marked
  increase in quantity in the G + L model suggested that Leu accelerated the formation of
  1-deoxyosone", in agreement with Kocadagli & Gokmen 2016, "who estimated at almost zero the rate
  of 1-deoxyosone formation directly from fructose".
- **5-Methylfurfural is a direct thermal degradation product of 5-HMF**, previously reported only
  at 250 and 400 C (Nikolov & Yaylayan 2011; Sanders 2003), "indicating that this degradation is
  even active at lower temperatures".
- **Furfuryl alcohol has three candidate routes and this paper favours the third**: via Amadori
  degradation (Wnorowski & Yaylayan 2000), from gluconic/glucuronic acids and glyceraldehyde
  (Yaylayan & Keyhani 2000), or **by reduction of furfural** (Delatour 2020) — the last "thus
  validating this third hypothesis", because FA tracks furfural across formulations.
- **The pyrazine stoichiometry, stated explicitly.** "For one mole of 3-MB, one mole of
  alpha-aminocarbonyl compound may form. This latter intermediate may undergo further condensation
  with **another alpha-aminocarbonyl compound to form dihydropyrazine**, which then leads to
  pyrazine." So **one pyrazine consumes two nitrogens**. Further: "**All alpha-dicarbonyl
  intermediates may give rise to 3-MB as their substituted groups are not structure-determinant.
  However, this is not the case for pyrazines**: the functional groups of the alpha-dicarbonyl
  intermediates were conserved in the corresponding alpha-aminocarbonyl compound and thus also in
  the corresponding pyrazine." That is why 3-MB is abundant and pyrazines are not: **3-MB is the
  common product of every Strecker event, whereas each pyrazine requires a specific pair.**
- **The glyoxal explanation.** Lee 2022 could not quantify glyoxal in G + L; here the authors argue
  "it was consumed as soon as it was formed to participate in Strecker degradation to form P, 2-MP
  and 3-MB", supported by P, 2-MP and 3-MB being "detected promptly in the vapor".
- **The melanoidin definition, repeated verbatim from the 2022 paper**: these compounds "then
  undergo further reactions to form **melanoidins, which differ from the brown polymers in that
  they contain nitrogen**".

**All concentration-time data are in Figures 1, 2, 4 and Figs. S1-S4 and are figure_only.**

### Arithmetic on the printed numbers (all mine)

**1. The furanic ratios, G + L : G (mine, from the printed levels at 200 C).** 5-MF and FA both
**5x**; acetic acid **~10x**; furfural **~1x** ("comparable"); 5-HMF **0.67x** (1.0 / 1.5), the
only marker that goes **down** with leucine. Compare the *matrix* ratios from the 2022 companion
at the same condition: 5-HMF **2.5x up**, furfural **4x up**. **The two compartments disagree in
sign for 5-HMF and in magnitude for furfural**, which is exactly the transport effect this paper
exists to explain — the paper's own reading is that the G + L matrix, being richer in
high-molecular-weight melanoidins, **retains** volatiles better. Do not compare a vapour ratio
with a matrix ratio without that correction.

**2. The vapour-pressure span versus the measured behaviour (mine).** From Table 1, furfural's
Pv is 2.21 mm Hg and 5-HMF's is 5.28e-3, a factor of **419**. Yet in the matrix (Lee 2022) 5-HMF
is 10x more abundant than furfural, while in the vapour furfural "was found in the highest
quantities". **The rank order inverts between the two compartments**, and the ratio of ratios is
roughly the ratio of vapour pressures. That is the single most useful sanity check in this paper:
**a volatile marker measured in the headspace is not a measure of what is in the food**, and the
error is of order the vapour-pressure ratio.

**3. Nitrogen accounting on the pyrazines (mine, order of magnitude only).** Each pyrazine carries
**two** nitrogens and each is built from two alpha-aminocarbonyls, i.e. two Strecker events. Each
Strecker event also releases one 3-MB. So the ratio 3-MB : pyrazine-nitrogen tells how much of the
Strecker flux ends as volatile nitrogen versus as retained aminocarbonyl. **The paper gives no
number for that ratio** — all pyrazine and 3-MB levels are figure-only, and 3-MB is off-scale at
200 C — so the accounting cannot be closed. **This is the highest-value missing number in the
paper** and it is why the supplementary figures and the raw data are requested in Flags 12.

**4. The bracket on a leucine melanoidin's C/N, both ends (mine).** Under the trunk's step-9
rule the sugar side is fixed at 6 carbons per nitrogen. The amine side then runs:

| route for the amine | carbons the amine contributes | resulting C/N of the repeat unit |
|---|---|---|
| leucine incorporated **intact** (C6 N1) | 6 | **12.0** |
| leucine **decarboxylated only** (C5 N1) | 5 | **11.0** |
| leucine **fully Strecker-degraded** — CO2 lost, C5 lost as 3-methylbutanal, only the nitrogen transferred | **0** | **6.0** |

**This paper measures 3-methylbutanal as the most abundant of all ten markers**, escaping into the
oven air, which pushes the real system towards the bottom row. **So the bracket is 6 to 12 and the
evidence here favours the low end.** For contrast, the same three rows for **glycine** are 8.0,
7.0 and 6.0 — and Mundt & Wedzicha measure **7.64 ± 0.21**, i.e. between the decarboxylated and
the intact case, because glycine's Strecker aldehyde is **formaldehyde**, which Mundt's
radiochemistry shows is *retained* in the polymer ("it is necessary to include in the melanoidins
carbon and nitrogen atoms from the amino carbonyl ... **as well as the Strecker aldehyde
itself**"). **Leucine's Strecker aldehyde is not retained — it is measured leaving.** That is the
cleanest statement this cluster can make about why melanoidin C/N is amine-specific:
**it depends on the volatility of the amine's Strecker aldehyde.** All of this is arithmetic on
the trunk's own constant plus two measured facts; **no leucine melanoidin C/N is measured
anywhere.**

**5. The transport model has more parameters than the data can carry (mine).** Panel C fits, per
compound: k_a,ref and Ea_a **for each of two formulations**, plus a shared k_d,ref and Ea_d —
**six parameters per compound**. The data are six vapour time points x three temperatures x two
formulations = **36 points per compound**, but the matrix input C_mj(t) is itself a
three-parameter sigmoid fitted separately per condition, and the water input C_mw(t) a
two-parameter exponential. Counting the interpolants, the model is heavily parameterised relative
to its evidence, and **k_d has no independent measurement at all** — the paper says outright it
was needed only because "this was the only way that the model could obtain vapor concentrations
which were lower at the end of the kinetic". Treat panel C as a description, not an estimate.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** This paper has the **best registry
coverage of any paper in the cluster**: of the ten markers, **eight are keyed** —
`3_methylbutanal`, `pyrazines` (the class id), `methylpyrazine`, `2_5_dimethylpyrazine`,
`2_6_dimethylpyrazine`, `furfural`, `hmf`, and — for acetic acid — nothing (see below). Not
keyed: **acetic acid** (though the trunk carries it as the B1 species `AA`), **5-methylfurfural**,
**furfuryl alcohol**, **pyrazine itself** (the registry keys the class `pyrazines` and several
substituted members but not unsubstituted pyrazine), and **glucose, fructose and leucine**.

**Governing conditions on every row: the same model sponge cake as Lee 2022 (corn starch + MC +
HPMC, 20 g batter, 1.3 cm thick), glucose with or without leucine, baked in a 96 L instrumented
pilot oven at an OVEN SET-POINT of 140, 170 or 200 C with a fan at 25 or 50 Hz; vapour drawn at
50 mL/min onto an Air Toxics sorbent tube over six 4-minute intervals and reported at each
interval's midpoint; TD-GC-MS against d4-pyrazine and d4-furfural. The PRODUCT temperature is a
non-isothermal profile reported only in the 2022 companion.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **3-MB, pyrazine, 2-MP, 2,5-DMP, 2,6-DMP present in G + L and ABSENT in G** | — | — | all three set-points | Results 3.2.2; Fig. 3 | **level_only (a measured presence/absence; the paper's central result)** |
| 3-MB rank | **most abundant of all ten markers**; above the calibration range at 200 C | — | all conditions | Results 3.2.2 | level_only |
| **acetic acid, G + L : G** | **~10 : 1** | — | 200 C, vapour | Results 3.2.1 | **within_study_ratio** |
| **5-methylfurfural, G + L : G** | **5 : 1** | — | 200 C, vapour | Results 3.2.1 | **within_study_ratio** |
| **furfuryl alcohol, G + L : G** | **5 : 1** (and "almost a 5-fold increase ... at 54 min, 200 C") | — | 200 C, vapour | Results 3.2.1, 3.3 | **within_study_ratio** |
| furfural, G + L : G | **~1 : 1** ("comparable in both formulae") | — | 200 C, vapour, except the 86 min point | Results 3.2.1 | within_study_ratio |
| **5-HMF, G + L : G** | **0.67 : 1 — the only marker LOWER with leucine** | — | 200 C, vapour | Results 3.2.1 | **within_study_ratio** |
| 5-MF maximum | 2.0 | mmol·L_air⁻¹·g_DM⁻¹ **as printed** (Flags 1) | 200 C, G + L | Results 3.2.1 | level_only |
| FA maximum | 1.5 | same | 200 C, G + L | Results 3.2.1 | level_only |
| FA maximum in G | 0.15 | same | 200 C, G | Results 3.2.1 | level_only |
| 5-HMF maximum | 1.0 (G + L) / 1.5 (G) | same | 200 C | Results 3.2.1 | level_only |
| 2,6-DMP : 2,5-DMP | **1 : 2** | — | G + L | Results 3.2.2 | within_study_ratio |
| pyrazine : 2,5-DMP | **~1 : 1** (apart from the 86 min point) | — | G + L | Results 3.2.2 | within_study_ratio |
| 2-MP rank | **lowest of all pyrazines** | — | G + L | Results 3.2.2 | level_only |
| **fan frequency has no significant effect on 8 of 10 markers** | p > 0.05 | — | 170 C, ANOVA on 72 observations | Table 2 | **level_only (a properly tested null)** |
| formulation (G vs G + L) effect | significant for **9 of 10** markers; not for 5-HMF (p = 0.326) | — | 170 C | Table 2 | within_study_ratio (statistical) |
| model fit quality by marker | R² from **0.363 (5-HMF)** to **0.994 (2-MP)** | — | 170 C | Table 2 | measured level |
| time shape at 200 C | rise to 56 min, then flat in G and **falling in G + L** for AA, F, 5-MF, FA, 5-HMF; pyrazines except P also fall after 56 min | — | 200 C | Results 3.2.1, 3.2.2 | level_only (temporal) |
| log Kow and vapour pressure for all 10 markers | see Table 1 | — | literature values (PubChem / EPI Suite) | Table 1 | **level_only (borrowed, not measured here)** |
| furfural : 5-HMF vapour pressure ratio | 419 | — | 25 C | derived from Table 1 (mine) | derived_assumption |
| rank inversion of furfural and 5-HMF between matrix and vapour | furfural highest in vapour; 5-HMF ~10x furfural in the matrix (2022) | — | 200 C | Results 3.2.1 + Lee 2022 (mine) | within_study_ratio (**cross-paper**) |
| oven volume and air flow | V_oven = 96 L; V̇_air = 0.05 L·min⁻¹ | L, L·min⁻¹ | the transport model | Methods 2.2; Results 3.4 | measured level |
| topology: acetic acid from 1-deoxyosone | — | — | — | Results 3.3 | level_only (structural) |
| topology: 1-deoxyosone forms faster via the Amadori route than via caramelisation | — | — | — | Results 3.3 | level_only (structural) |
| topology: 5-methylfurfural is a direct degradation product of 5-HMF, active below 250 C | — | — | — | Results 3.3 | level_only (structural) |
| topology: furfuryl alcohol arises by reduction of furfural | — | — | — | Results 3.3 | level_only (structural, the paper's preferred hypothesis of three) |
| stoichiometry: 1 Strecker event → 1 3-MB + 1 alpha-aminocarbonyl; 2 aminocarbonyls → 1 pyrazine (**2 N per pyrazine**) | — | — | — | Results 3.3 | level_only (structural) |
| leucine repeat-unit C/N bracket under the trunk's step-9 rule | **6.0 (full Strecker) to 12.0 (intact)** | mol C per mol N | — | arithmetic on `species.py` plus this paper's 3-MB result (mine) | **derived_assumption** |
| all concentration-time courses, the model fits and the identified parameters | — | — | — | Figs. 1, 2, 4 (incl. panel C), S1-S4 | **figure_only** |

### How this bears on the C/N diagnostic and the trunk

**(a) It names a `MEL_N` leak the trunk does not have.** Nitrogen leaves the matrix as pyrazines
(two nitrogens each) and is measured doing so. Fang 2009 independently shows that nitrogen is
**not** in the polymer. The trunk routes all unmeasured amine nitrogen to `MEL_N`; in a real baked
matrix some fraction of it is in the oven air. **This paper cannot size the leak** (all levels are
figure-only and 3-MB is off-scale at the top condition), but it establishes that it exists and is
measurable, and it gives the stoichiometry (2 N per pyrazine, 1 Strecker event per 3-MB) that
would let it be sized from the deposited data.

**(b) It fixes the direction of the leucine C/N argument.** The `lee2022_extraction.md` dossier's
note that a leucine repeat unit would be 12 C per N is the **upper** bound only. With leucine's
Strecker aldehyde measured leaving the cake in the largest quantity of any marker, the realistic
range is **6 to 12**, and the mechanism that decides where in that range a real system sits is
**the volatility of the amine's Strecker aldehyde** — glycine's formaldehyde is retained (Mundt),
leucine's 3-methylbutanal is not (this paper). **Mundt & Wedzicha's 7.64 ± 0.21 is a glycine
number and must not be generalised to another amine.**

**(c) It gives the strongest available warning against reading headspace data as composition.**
Furfural and 5-HMF swap rank between the matrix and the vapour, by roughly their vapour-pressure
ratio. Any repository lane that scores a volatile marker measured in a headspace against a model
concentration in a matrix is comparing two different quantities, and this paper quantifies how
different.

**(d) Nothing here is a chemical rate.** The only fitted constants are the transport constants of
Figure 4C, which lump a partition coefficient, an evaporation rate and a first-order term, carry
an unstated T_ref, and include a disappearance constant with no identified mechanism. **They must
not enter `parameters.py` or any lane.**

**(e) What could become benchmark rows.** The **G + L : G vapour ratios** (acetic acid 10x,
5-MF 5x, FA 5x, furfural 1x, 5-HMF 0.67x) are within-study, matrix-matched, condition-matched and
internally-standardised, and a formulation ratio cancels any trapping-efficiency error that is
common to the two formulations. The **presence/absence result** for the five nitrogen-bearing
volatiles is a categorical row and is the cleanest thing in the paper.

## 5. Flags

1. **The vapour-concentration unit is printed two incompatible ways, differing by 10⁶.** Section
   3.2.1 says the concentrations are "expressed in **mmol** per L in the oven and per g of dry
   matter", written `mmol·L_air⁻¹·g_DM⁻¹` (confirmed from the page render, not an OCR artefact).
   Section 3.4 says "the values in vapor previously shown as **nmol**·L_air⁻¹·g_DM⁻¹ were
   multiplied by 10⁻³ and V_oven (L)" to reach the `umol·g_DM⁻¹` of Figure 4. Given that
   Figure 4's y-axes top out at a few umol·g_DM⁻¹ and V_oven = 96 L, **the nmol reading is the
   self-consistent one and the "mmol" of section 3.2.1 is almost certainly an error** — but that
   is an inference, so section 4 transcribes the printed "mmol" and flags it here. **Do not use
   any absolute vapour level from this paper until the unit is confirmed.** The ratios are
   unaffected.
2. **The only fitted constants in the paper are figure-trapped, mechanism-free and missing their
   reference temperature.** Figure 4C holds k_a,ref and Ea_a per formulation and a shared k_d,ref
   and Ea_d, for furfural and 5-HMF. **T_ref is never printed anywhere** — Eq. 5 defines it only
   as "chosen most of the time as the center for the experimental temperature domain" — so the
   Arrhenius form cannot be evaluated. `k_a` is by the authors' own statement a lump of three
   distinct physical quantities, and `k_d` is fitted with no independent evidence, its mechanism
   listed as "deposition of molecules on the oven walls, some further reactions occurring in the
   vapor or even gas leaks". **These are transport-model descriptors, not kinetics.**
3. **Every vapour datum is a 4-minute integral reported at its midpoint.** The tubes trap over
   2-6, 8-12, 20-24, 33-37, 52-56 and 86-90 min and each interval is plotted at 4, 10, 22, 35, 54
   and 88 min. In a system whose product temperature is still climbing steeply through the first
   30 min, a 4-minute window is not a point. **Six points per curve, each an average over a
   changing temperature.**
4. **The internal standards are added AFTER sampling and correct for the wrong step.** d4-pyrazine
   and d4-furfural are spiked onto the tube once the tube has been removed from the oven. They
   correct for thermal desorption, chromatography and ionisation. **They do not correct for the
   trapping efficiency of the sorbent tube under baking conditions**, which is where a
   temperature- and humidity-dependent bias would live — and the tube sees very different water
   loads at 4 min and at 88 min.
5. **Replication is uneven and is described differently from the 2022 companion.** Here:
   140 C/25 Hz, 170 C/25 Hz, 170 C/50 Hz and 200 C/25 Hz in triplicate; **140 C/50 Hz and
   200 C/50 Hz once only**. The figure captions then report n = 4, 6 and 4 at the three
   temperatures after pooling fan settings — so at 140 C and 200 C the pooled n is 3 + 1.
6. **Table 1's quantifier ions cannot be recovered from this file.** The footnote says "quantifier
   ions are in bold", and the bold formatting is lost in text extraction. Six of the twelve rows
   list two or three ions with no way to tell which is the quantifier. Read them from the typeset
   PDF or from Lee 2021 before transferring the method.
7. **Two printed identifier errors.** **PubChem CID 176 is given for both acetic acid and
   5-methylfurfural**; and "Furfual" appears for furfural, as it does in the 2022 companion. The
   CAS numbers in Table 1 are unambiguous and should be preferred. Also "2,6–imethylpyrazine" for
   2,6-dimethylpyrazine in both Table 1 and Table 2, and "Turkey's test" for Tukey's.
8. **The leucine C/N bracket in section 3 arithmetic 4 is mine, and no leucine melanoidin
   composition is measured anywhere in this cluster.** The bracket rests on the trunk's own
   step-9 stoichiometry, on standard Strecker chemistry, and on two measured facts (3-MB leaves
   here; formaldehyde is retained in Mundt). **It is a scoping calculation, not an extraction
   result**, and it is recorded because the 2022 dossier's single figure of 12 is only the upper
   bound and would otherwise stand alone.
9. **This paper measures nothing in the food.** Every number is a vapour concentration or a
   statistic on one. Its matrix concentrations are taken wholesale from Lee 2022, and its
   temperature history is not reported here at all. **The two papers must be read as one study**,
   and the conditioning caveats of `lee2022_extraction.md` (oven set-point is not product
   temperature; the leucine concentration is never printed; no moisture or a_w value; the product
   has a 30-50 min surface-to-centre lag) all apply.
10. **The pyrazine ANOVA is weak where it matters most.** Unsubstituted pyrazine has
    **R² = 0.490** and **no significant Time effect (p = 0.186)**, and 5-HMF has **R² = 0.363**
    with no significant Formula effect. Two of the ten markers are therefore poorly described by
    the model that the paper uses to justify pooling the fan settings, and one of those two is the
    parent pyrazine — the compound whose presence/absence carries the paper's nitrogen argument.
    The presence/absence result does not depend on the ANOVA; the level trends do.
11. **What this paper does not contain**: any chemical rate constant or activation energy; any
    tabulated concentration; any measurement in the solid matrix; any temperature measurement;
    any browning or colour measurement; any elemental analysis, C/N or polymer composition; any
    moisture or water-activity value; any leucine concentration; any nitrogen mass balance; any
    quantitative 3-MB value at 200 C (off-scale); any amine other than leucine; any sugar other
    than glucose; any isothermal condition. **The paper's data-availability statement reads
    "Data will be made available on request"** — unlike the 2022 companion, whose data are
    deposited under DOIs.
12. **What to request from the authors**, in priority order: (i) **the numeric vapour
    concentration series behind Figures 1, 2, S1 and S2**, with the unit confirmed — this would
    settle Flag 1 and would let the pyrazine/3-MB nitrogen accounting of section 3 arithmetic 3
    actually be closed; (ii) **Figure 4C's parameters with their T_ref**, in a typeset form;
    (iii) an extended calibration for 3-MB so the 200 C points are on scale; (iv) the
    supplementary figures S1-S4; (v) the quantifier-ion assignments of Table 1; (vi) whether the
    sorbent-tube trapping efficiency was ever checked against the changing water load during
    baking.
13. **Registry gaps against `data/keys/compounds.yml`**: keyed and covered — `3_methylbutanal`,
    `methylpyrazine`, `2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`, `pyrazines` (class),
    `furfural`, `hmf`. **Not keyed: acetic acid** (which the trunk carries as the B1 species `AA`
    — the same registry-versus-trunk gap the whole cluster records), **5-methylfurfural**,
    **furfuryl alcohol**, **unsubstituted pyrazine**, **glucose**, **fructose** and **leucine**.
    Unsubstituted pyrazine is the notable omission: it is the parent of a keyed class, it is one of
    the two compounds whose deuterated analogue is an internal standard here, and it is one of the
    five markers that carry this paper's nitrogen result.
