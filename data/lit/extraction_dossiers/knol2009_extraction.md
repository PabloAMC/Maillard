# Knol 2009 — EXTRACTION (potato crisps, three cold-sweetened genotypes, 1.5 mm slices deep-fried in rapeseed oil from 180 C down to a held 160 C, 1-12 min; three EMPIRICAL acrylamide models — Logistic-Fermi, Logistic-Exponential, Empirical — fitted with SD; precursor sugars and asparagine in mg/g dm)
### The real-food half of the scorecard's gap, and a warning label: this paper contains NO rate constant and NO activation energy — its "k1" and "k2" are logistic steepness parameters of a curve fitted without any mechanism — but it does carry four genotypes' acrylamide scale parameters (9.3e3 to 2.6e4 ug/kg dm), their precursor concentrations, and a regression that predicts one from the other.

**Source on disk:** `data/articles/knol2009.pdf` (7 pp., Food Chemistry 113 (2009) 103-109).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/knol2009.txt`). **Table 2 is
rotated 90 degrees on the page and the `-layout` extraction shredded it into fourteen
column-fragments in reverse order.** It was re-extracted with `pdftotext -f 5 -l 5 -raw`, which
returns the table row by row and cleanly; the two readings agree cell for cell, and the `-raw`
version is what is transcribed below. Tables 1 and 3 came through clean in `-layout`. Figures 1
(temperature profiles of oil and crisp surface), 2 (acrylamide vs frying time, three genotypes),
3 (water content vs time), 4 (the three model fits), 5 (parameter a vs sugars), 6 (predicted vs
experimental acrylamide) and 7 (L\*a\*b\* vs acrylamide) are images: **every acrylamide
concentration measured in this study is figure-only.** No supplementary material. Repo status
before this dossier: Knol 2009 appears in `k3_final_parameter_inventory.md` (§A.2 line 207 as one
of the three papers whose degradation parameters are REFUSED; §B8.5 for the authors' own refusal to
transfer; D.5 where the "real-food band 9.3 x 10^3-2.6 x 10^4 ug/kg dm" is declared **HOLD-OUT**
and the degradation parameters "neither"). It has **no extraction dossier**.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetic modelling: A tool to predict the formation of acrylamide in potato crisps" |
| Authors | Jeroen J. Knol (a), Gunilla Å. I. Viklund (b), Jozef P. H. Linssen (a, corresponding), Ingegerd M. Sjöholm (b), Kerstin I. Skog (b), Martinus A. J. S. van Boekel (a) — (a) Product Design and Quality Management Group, Wageningen University, NL; (b) Department of Food Technology, Engineering and Nutrition, Lund University, SE |
| Venue | Food Chemistry 113 (2009) 103-109. Received 18 April 2008, revised 9 June 2008, accepted 14 July 2008 |
| DOI | 10.1016/j.foodchem.2008.07.032 |
| Naming | a(T) = "scale factor" for the acrylamide concentration, in ug/kg; k1(T), k2(T) = **steepness parameters around the inflection points**, NOT rate constants; tc1(T), tc2(T) = **time characteristics for the inflection points** in formation and degradation; tau(T) = "characteristic time" of the exponential decay term; a, b, m, n = the four constants of the purely empirical model; dm = dry matter |
| Lineage | the three model forms are Corradini & Peleg 2006 (Crit. Rev. Food Sci. Nutr. 46:489); the experimental system is Viklund et al. 2007 (LWT 40:1066) and Viklund et al. 2008 (JSFA 88:305); the immediate predecessor is **Knol et al. 2008, Mol. Nutr. Food Res. 52:313-321** (the Bintje row of Tables 2 and 3), which is **NOT on disk** |
| Companions on disk | `knol2005_extraction.md` (the mechanistic aqueous model this paper explicitly contrasts itself with), `knol2010_extraction.md`, `claeys2005_extraction.md`, `devleeschouwer2009_extraction.md` |

## 1. Why it matters

The scorecard names the acrylamide lane's gap as "a second laboratory's constants; **real-food
matrices**" (`results/validation/data_wishlist.md`). The shipped lane
(`src/kinetic_core/parameters_acrylamide.py`, `src/kinetic_core/acrylamide.py`) is a named
mass-action network: asparagine + a carbonyl -> a Schiff-base / decarboxylated-Amadori route ->
acrylamide, with an elimination step and measured cysteine channels. **This paper is a real food,
and it is the one paper in this five-paper set that is.** It is also the one that supplies nothing
the network can be charged with.

That is not a defect of the extraction; it is the paper's own design and its own conclusion.
"Using models that only give a mathematical description of the formation or degradation of
acrylamide in food bypasses the problem of considering all the mechanisms that occur during
processing of foods." The parameters "were estimated by the software in such a way that the
residual sum of squares was minimised **without any mechanistic restraints, such as knowledge of
the underlying reaction networks**." The authors then state the transfer refusal that
`k3_final_parameter_inventory.md` §B8.5 already quotes, and it is confirmed here verbatim from two
places: "these mechanistic studies in model systems are not easily applied to real food systems"
(Introduction) and "in the attempt to model the formation of acrylamide one should realise that
**the parameters are only applicable for specific experimental conditions such as time-temperature
profile of frying, potato genotype, slice thickness and initial concentration of precursors**"
(Results 3.2).

So what the repository gets from this paper is four things, none of them a rate constant:

1. **A real-food acrylamide band with its precursor composition.** Four genotypes, four scale
   parameters a from 9.3 x 10^3 to 2.6 x 10^4 ug/kg dm, each paired with a measured reducing-sugar
   and asparagine content in mg/g dm (Table 1). That is exactly the shape of a hold-out benchmark:
   a matrix the lane can be charged with and a target it must not have been fitted to. The
   inventory's HOLD-OUT declaration is correct and this dossier supplies the composition needed to
   charge it (section 4).
2. **A quantitative precursor-to-product relation in a real matrix**: a = 1.06 x C_red.sugars -
   2.04, R^2 = 0.98 against reducing sugars (0.97 glucose, 0.91 fructose), across four genotypes.
   The shipped lane's whole claim to be mass-action rather than yield-based is that it "responds to
   precursor CONCENTRATION"; **this is a measured, real-food test of exactly that response**, and
   it says the response is close to linear in reducing sugar and essentially blind to asparagine
   (which varied only 9.93-11.60 mg/g dm here).
3. **A negative result on the acrylamide DEGRADATION that matters for `k_acr_dp`.** Every
   degradation parameter in this paper is unidentifiable: k2 = 3.5 ± 8 and tc2 = 0.60 ± 52 in the
   Logistic-Fermi model, tau = 22 ± 15 in the Logistic-Exponential. The authors ran Monte Carlo
   simulations, found the estimates approximately normal and **not** strongly correlated, and
   concluded "the low precision of the parameter estimates must be mainly due to the fact that the
   data do not contain enough information to extract precise model parameters for the degradation."
   They had *lengthened* the experiment to 12 min specifically to fix this, and it did not. The
   registry's decision to FIT `k_acr_dp` against three sources rather than transcribe one is
   supported.
4. **The one shape fact the panel's inverted-time-shape problem needs.** In a real crisp at a held
   160 C, acrylamide **rises steeply between 1.5 and 4 min and then decreases slowly** — the
   decrease is real but shallow, and the Logistic-Exponential model was preferred over the
   Logistic-Fermi precisely because the Fermi form drives acrylamide to zero at long times and "the
   Logistic-Exponential function with its exponential term for the degradation predicts a residual
   acrylamide concentration at prolonged heating times", which is what the data show. That is a
   qualitative constraint on `k_acr_dp` in a food matrix: **not zero, but not fast enough to
   consume the pool.**

What this paper is NOT: it is not a second laboratory's *constants*, because it has no constants in
the kinetic sense. Its "k1" has units of min^-1 and a value near 3, but multiplying an acrylamide
concentration by it would be meaningless — it is the steepness of a logistic curve at its
inflection point. **Nothing in Table 2 may enter `parameters_acrylamide.py`.**

## 2. Methods as they matter to a model

- **Matrix.** Potato tubers (*Solanum tuberosum* L.), genotypes **Hulda, Lady Rosetta and
  SW 91 102** (Viklund et al. 2008), grown in southern Sweden, harvested September 2006, wound-
  healed at 15 C for 2 weeks, temperature slowly reduced, then **stored at 4 C at 90-95 % RH for
  4 weeks** through November 2006 — i.e. deliberately **cold-sweetened**, which is why the reducing
  sugars are high (9.3-15.8 mg/g dm). The fourth genotype in Tables 2 and 3, **Bintje**, is not
  new work: it is Knol et al. 2008, included "as the experimental set-up was the same".
- **Preparation.** Tubers 6-10 cm diameter, 100-250 g; washed; **cut into 1.5 mm slices**;
  randomised into three batches per genotype. One batch fried; the other two analysed for
  precursors at the start (morning) and end (evening) of the experiment day.
- **Frying.** Rapeseed oil; **slices in three net cages stacked on top of each other**; oil
  preheated to **180 C**; on immersion the oil temperature fell rapidly; "After about 1.5-2.0 min,
  the oil reached a temperature of 160 C, and the thermostat setting was then changed to keep the
  temperature constant at 160 C for the rest of the frying session." Frying times **1, 1.5, 2, 2.5,
  3, 3.5, 4, 4.5, 5, 6, 8 and 12 min**; **duplicate experiments at 2, 2.5, 3 and 4 min**. All in one
  day. Crisps cooled to room temperature and stored in closed plastic bags at **-18 C**.
  **This is a strongly non-isothermal, non-isobaric, drying process and the models absorb all of it
  into empirical constants** — there is no separation of heat transfer, mass transfer and
  chemistry.
- **Temperature.** Oil and the **outer cell layer** of a slice logged by thermocouples (Knol et al.
  2008 method); Fig. 1 shows both plus their difference, against the earlier study. The text notes
  "The rapid increase of acrylamide in the crisps took place after the rapid decrease of the water
  content and when the temperature of the outer cell layer of the potato crisps had reached 160 C."
  **Fig. 1 is the only record of the temperature history and it is figure-only.**
- **Water activity / moisture.** Initial water content **78 % (Lady Rosetta and SW 91 102), 80 %
  (Hulda)**; after 1 min of frying **41 % (Lady Rosetta), 34 % (SW 91 102), 30 % (Hulda)**; "almost
  the same water content of the crisps after 3 min frying"; at longer times SW 91 102 slightly
  higher at **5 %**; industrial crisps are cited at ~2 %. Gravimetric, AOAC-984.25.
  **The a_w sweeps from ~1.0 to well below 0.3 during a single run.** The lane's declared
  water-activity window in `acrylamide_conditions.py` is 0.34-0.99 and it refuses to extrapolate
  outside it; a crisp finishes outside it. This is the single hardest obstacle to using this paper
  as a benchmark and it is a property of the food, not of the paper.
- **Acrylamide.** **LC-MS/MS** (Viklund et al. 2007). One extraction per sample, extract analysed
  in duplicate. **LOD 50 ug/kg dry matter** ("In crisps fried for 1 min, the acrylamide
  concentrations were below the detection limit of 50 ug/kg dry matter (dm)"). **All results are
  expressed per kg DRY MATTER.**
- **Precursors.** Fructose, glucose and sucrose by **GC-FID**; asparagine by **HPLC with
  fluorescence detection** (Olsson, Svensson & Roslund 2004). Duplicate analyses. Morning vs
  evening differences not significant (p > 0.05) except Hulda sucrose (p = 0.0072); "the average
  concentrations of all the measurements were used for the empirical modelling."
- **Colour.** Non-destructive computer vision system (Viklund et al. 2007), L\*a\*b\*.
- **Fitting.** Mathcad v13.1, `Minerr`, quasi-Newton, sum-of-squares minimisation. **SD of the
  parameters by linear approximation from the variance-covariance matrix and mean squares** — not
  HPD intervals (contrast Knol 2005 and 2010, same first author, which use Athena Visual Studio and
  the determinant criterion). Model discrimination by **corrected Akaike information criterion**,
  AIC = n ln(SS/n) + 2(p+1) and AICc = AIC + 2(p+1)(p+2)/(n-p), because n/p < 40 in every case.
  Precursor-content statistics in Minitab 13 (general linear model + Tukey).
- **The three models, as printed.**
  - **Logistic-Fermi (eq 1)**, 5 parameters:
    C(t) = a/(1+exp{k1[tc1 - t]}) - a/(1+exp[k1 tc1]) , all multiplied by 1/(1+exp{k2[t - tc2]}).
  - **Logistic-Exponential (eq 2)**, 4 parameters: the same logistic formation term multiplied by
    exp(-t/tau).
  - **Empirical (eq 3)**, 4 parameters: C(t) = a t^n / (b + t^m).
  **None of these is a rate law.** There is no concentration of any reactant in any of them; t is
  frying time in minutes and every other symbol is a fitted shape parameter.

## 3. Tables re-typed

### Table 1. "Acrylamide precursor concentrations (mg/g dm) with their variation (n = 2) of the raw potato slices from the three different genotypes used at the beginning of the experiments (morning) and at the end of all experiments (evening)"

| compound | genotype | morning | evening |
|---|---|---|---|
| Fructose | Hulcla *(printed; = Hulda, Flags 5)* | 3.33 ± 0.06 | 3.97 ± 0.13 |
| | Lady Rosetta | 7.28 ± 0.08 | 7.78 ± 0.36 |
| | SW 91 102 | 4.88 ± 0.15 | 4.72 ± 0.05 |
| Glucose | Hulda | 5.94 ± 0.11 | 6.21 ± 0.23 |
| | Lady Rosetta | 8.09 ± 0.03 | 8.03 ± 0.06 |
| | SW 91 102 | 7.34 ± 0.30 | 7.00 ± 0.04 |
| Sucrose | Hulda | 6.89 ± 0.16 | 5.26 ± 0.12 (a) |
| | Lady Rosetta | 7.39 ± 0.32 | 7.34 ± 0.15 |
| | SW 91 102 | 10.45 ± 0.36 | 11.51 ± 0.39 |
| Reducing sugars | Hulda | 9.27 ± 0.05 | 10.17 ± 0.10 |
| | Lady Rosetta | 15.37 ± 0.11 | 15.80 ± 0.42 |
| | SW 91 102 | 12.21 ± 0.45 | 11.71 ± 0.01 |
| Asparagine | Hulda | 10.53 ± 0.42 | 10.47 ± 0.25 |
| | Lady Rosetta | 11.60 ± 0.04 | 11.42 ± 0.04 |
| | SW 91 102 | 9.93 ± 0.07 | 10.05 ± 0.01 |

Footnote a: "Significant difference between morning and evening samples, p < 0.05."
Units **mg per g dry matter** throughout. Bintje's composition is not given here (it is in Knol
et al. 2008, not on disk).

**Conversions to mmol per g dm (mine).** Fructose and glucose 180.16 g/mol, sucrose 342.30,
asparagine 132.12. Morning/evening averages: **Hulda** reducing sugars 9.72 mg/g = **0.054 mmol/g
dm**, asparagine 10.50 mg/g = **0.0795 mmol/g dm**; **Lady Rosetta** 15.59 mg/g = **0.0865
mmol/g dm**, asparagine 11.51 mg/g = **0.0871 mmol/g dm**; **SW 91 102** 11.96 mg/g = **0.0664
mmol/g dm**, asparagine 9.99 mg/g = **0.0756 mmol/g dm**. So the **molar reducing-sugar : asparagine
ratio is 0.62 : 1 (Hulda), 0.99 : 1 (Lady Rosetta), 0.88 : 1 (SW 91 102)** — near-equimolar, which
is the ratio all three of the model-system papers in this set were built at, and sugar-limited in
Hulda.

### Table 2. "Estimates of parameters with their approximate SD (obtained by linear approximation) for the Logistic-Fermi, Logistic-Exponential and empirical models"

Rotated on the page; transcribed from the `-raw` extraction, which returns it row-wise.
Footnote a: "Data from Knol et al. (2008)." Units as printed: a in **ug/kg**, k1 and k2 in
**min^-1**, tc1, tc2 and tau in **min**; n, b and m carry no unit.

**Logistic-Fermi (5 parameters)**

| genotype | a (ug/kg) | k1 (min^-1) | tc1 (min) | k2 (min^-1) | tc2 (min) |
|---|---|---|---|---|---|
| Bintje (a) | 1.9 x 10^4 ± 8 x 10^2 | 3.0 ± 0.5 | 2.7 ± 0.1 | 3.5 ± 8 | 6.6 ± 1 |
| Lady Rosetta | 2.7 x 10^4 ± 4 x 10^4 | 3.0 ± 0.4 | 2.8 ± 0.1 | 0.053 ± 0.05 | 0.60 ± 52 |
| Hulda | 8.5 x 10^3 ± 2 x 10^2 | 2.8 ± 0.3 | 2.2 ± 0.04 | 1.7 ± 4 | 13 ± 2 |
| SW 91 102 | 7.1 x 10^3 ± 9 x 10^2 | 1.6 ± 0.4 | 3.0 ± 0.2 | 0.63 ± 1 | 14 ± 3 |

**Logistic-Exponential (4 parameters — the preferred model)**

| genotype | a (ug/kg) | k1 (min^-1) | tc1 (min) | tau (min) |
|---|---|---|---|---|
| Bintje (a) | 2.6 x 10^4 ± 6 x 10^3 | 2.8 ± 0.5 | 2.8 ± 0.1 | 15 ± 10 |
| Lady Rosetta | 1.4 x 10^4 ± 9 x 10^2 | 3.0 ± 0.4 | 2.8 ± 0.1 | 31 ± 9 |
| Hulda | 9.6 x 10^3 ± 5 x 10^2 | 2.5 ± 0.3 | 2.3 ± 0.1 | 41 ± 14 |
| SW 91 102 | 9.3 x 10^3 ± 2 x 10^3 | 1.5 ± 0.4 | 3.1 ± 0.3 | 22 ± 15 |

**Empirical (4 parameters)**

| genotype | a | n | b | m |
|---|---|---|---|---|
| Bintje (a) | 4.1 x 10^4 ± 2 x 10^4 | 6.2 ± 1 | 1.1 x 10^3 ± 1 x 10^3 | 6.7 ± 1 |
| Lady Rosetta | 2.5 x 10^4 ± 2 x 10^4 | 4.9 ± 2 | 3.5 x 10^2 ± 4 x 10^2 | 5.3 ± 2 |
| Hulda | 1.3 x 10^4 ± 7 x 10^3 | 4.9 ± 3 | 9.2 x 10^1 ± 2 x 10^2 | 5.1 ± 2 |
| SW 91 102 | 2.6 x 10^4 ± 3 x 10^5 | 3.2 ± 13 | 1.8 x 10^2 ± 2 x 10^3 | 3.8 ± 9 |

**The whole empirical block is unidentified**: every SD is at or above its estimate except two,
and SW 91 102's a carries an SD ten times the estimate. The authors say so.

### Table 3. "Model discrimination results for the Logistic-Fermi, Logistic-Exponential and empirical models"

Footnotes: A "Data from Knol et al. (2008)."; a "Residual mean squares."; b "Residual sum of
squares."; c "Corrected Akaike information criterion." The AICc values are **positive as printed**
and the ΔAIC column is internally consistent with them: for each genotype ΔAIC = AICc - min(AICc)
across the three models, checked row by row (Bintje 311.3 - 308.3 = 3.0; Lady Rosetta
447.4 - 444.4 = 3.0 and 450.6 - 444.4 = 6.2; Hulda 415.6 - 410.5 = 5.1 and 412.4 - 410.5 = 1.9;
SW 91 102 465.9 - 463.1 = 2.8 ~ 2.7 and 463.5 - 463.1 = 0.4). Lowest wins.

| genotype | | Logistic-Fermi | | | | | Logistic-Exponential | | | | | Empirical | | | |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| | p | n | MS | SS | AICc | ΔAIC | p | n | MS | SS | AICc | ΔAIC | p | n | MS | SS / AICc / ΔAIC |
| Bintje (A) | 5 | 20 | 3.3 x 10^6 | 5.0 x 10^7 | 311.3 | 3.0 | 4 | 20 | 3.3 x 10^6 | 5.3 x 10^7 | 308.5 | 0.2 | 4 | 20 | 3.3 x 10^6 | 5.2 x 10^7 / 308.3 / 0.0 |
| Lady Rosetta | 5 | 32 | 9.0 x 10^5 | 2.4 x 10^7 | 447.4 | 3.0 | 4 | 32 | 8.7 x 10^5 | 2.4 x 10^7 | 444.4 | 0.0 | 4 | 32 | 1.5 x 10^7 | 2.9 x 10^7 / 450.6 / 6.2 |
| Hulda | 5 | 32 | 2.8 x 10^5 | 7.7 x 10^6 | 410.5 | 0.0 | 4 | 32 | 3.5 x 10^5 | 9.9 x 10^6 | 415.6 | 5.1 | 4 | 32 | 4.7 x 10^6 | 8.9 x 10^6 / 412.4 / 1.9 |
| SW 91 102 | 5 | 32 | 1.6 x 10^6 | 4.3 x 10^7 | 465.9 | 2.7 | 4 | 32 | 1.6 x 10^6 | 4.4 x 10^7 | 463.1 | 0.0 | 4 | 32 | 2.2 x 10^8 | 4.4 x 10^7 / 463.5 / 0.4 |

Verdict in the text: "The model with the highest number of parameters (Logistic-Fermi) performs
less, though not substantially ... we prefer the Logistic-Exponential and empirical model. These
perform equally well ... If we take the results from the discrimination test and the precision of
the parameter estimation the Logistic-Exponential model is favoured." **Note that Hulda's own
ΔAIC prefers the Logistic-Fermi (0.0 against 5.1), against the overall verdict.**

### Equation 6 and the numbers in the running text

**a = 1.06 x C_red.sugars - 2.04** (eq 6), "where C_red.sugars is the total concentration (mg/g dm)
of the reducing sugars in the raw potato slices". **The unit of a is not restated in eq 6.** For
the equation to reproduce Table 2 it must return a in **10^3 ug/kg dm** (mine: Lady Rosetta,
1.06 x 15.59 - 2.04 = 14.5, against the printed a = 1.4 x 10^4 ug/kg — see Flags 2).

| quantity | value | where |
|---|---|---|
| acrylamide LOD | 50 ug/kg dm | Results 3.2 |
| acrylamide at 1 min | below LOD in all three genotypes | Results 3.2 |
| acrylamide time shape | "Frying between 1.5 and 4 min made the acrylamide concentrations increase rapidly for all three genotypes. **After 4 min, the acrylamide concentration decreased slowly**, with the largest decrease found in the Lady Rosetta crisps." | Results 3.2 |
| ranking | Lady Rosetta > Hulda > SW 91 102 in acrylamide | Results 3.2 |
| non-monotonicity vs sugar | SW 91 102 had higher fructose and glucose than Hulda yet **less** acrylamide | Results 3.2 |
| a vs precursors | R^2 = **0.91** (fructose), **0.97** (glucose), **0.98** (reducing sugars); "no clear relationship could be found" for any other parameter | Results 3.3, Fig. 5 |
| oil profile | preheated 180 C; reaches 160 C after 1.5-2.0 min; held at 160 C thereafter; industry cited at 180-190 C start, 150-175 C end | Methods 2.1, Results 3.1 |
| water content | 80 % (Hulda) / 78 % (LR, SW) raw; 30 / 41 / 34 % at 1 min; equal by 3 min; ~5 % for SW at long times; industrial crisps ~2 % | Results 3.2, Fig. 3 |
| colour vs acrylamide | R^2 = 0.28 (L\*), 0.82 (a\*), 0.85 (b\*), Hulda only | Results 3.3, Fig. 7 |
| prediction test | Logistic-Exponential with Hulda's k1, tc1, tau plus a from eq 6, tested against Viklund et al. 2008 crisps (Hulda, stored 6-24 weeks, fried 4 min): "close to the experimental values with some exceptions where the difference between the values was more than 15 %" | Results 3.3, Fig. 6 |
| transfer refusal | "the parameters are only applicable for specific experimental conditions such as time-temperature profile of frying, potato genotype, slice thickness and initial concentration of precursors"; and the relation "is based on four different genotypes. Ideally, this relationship should be established by using only one genotype with a large range of different sugar concentrations" | Results 3.2, 3.3 |

**Acrylamide concentration-versus-time data: FIGURE-ONLY** (Fig. 2, with SD bars, n = 4 at 2, 2.5,
3 and 4 min and n = 2 elsewhere; Fig. 4 with the three fitted curves). Water content: figure-only
except the values quoted above. Temperature profiles: figure-only. Colour: figure-only.

### Arithmetic on the printed parameters (all mine)

**1. What a means, and therefore what the real-food band is.** In both logistic models a is the
plateau of the formation term before degradation is applied; the second bracket in eqs 1 and 2
subtracts the t = 0 offset. So **a is very nearly the maximum acrylamide the crisp would reach with
no degradation**, in ug/kg dm. The Logistic-Exponential values are **2.6 x 10^4 (Bintje),
1.4 x 10^4 (Lady Rosetta), 9.6 x 10^3 (Hulda), 9.3 x 10^3 (SW 91 102) ug/kg dm** — this is exactly
the "9.3 x 10^3-2.6 x 10^4 ug/kg dm" band `k3_final_parameter_inventory.md` D.5 declares HOLD-OUT,
and it is confirmed here from the primary table.

**2. Equation 6 does not reproduce Table 2 well, and its unit is unstated.** With the
morning/evening average reducing sugars from Table 1: Lady Rosetta 15.59 -> 14.5 against a printed
14 (ratio 1.04); Hulda 9.72 -> 8.26 against 9.6 (0.86); SW 91 102 11.96 -> 10.64 against 9.3
(1.14), all in 10^3 ug/kg dm. **The R^2 = 0.98 is carried by Bintje**, whose composition is not in
this paper; inverting eq 6 on Bintje's a = 26 gives an implied **C_red.sugars = 26.5 mg/g dm
(mine)**, roughly 1.7x the highest genotype measured here. A four-point regression whose leverage
sits on the one point whose x-value is not printed is a weak instrument, and the authors say as
much ("Ideally, this relationship should be established by using only one genotype with a large
range of different sugar concentrations").

**3. Degradation is unidentified in every model.** Relative SD on the degradation parameters:
Logistic-Fermi k2 = 229 % (Bintje), 94 % (LR), 235 % (Hulda), 159 % (SW); tc2 = 15 %, **8667 %**,
15 %, 21 %; Logistic-Exponential tau = 67 %, 29 %, 34 %, **68 %**. By contrast the formation
parameters are tight: tc1 to 1-10 % and k1 to 12-27 % everywhere. **The data determine when
acrylamide appears and how steeply, and say almost nothing about how it goes away** — over a
12-minute fry.

**4. tau against the model-system elimination constants (mine, and the comparison is illegitimate
as a rate).** If one *forced* the Logistic-Exponential's exp(-t/tau) to be read as a first-order
decay, 1/tau would be 0.067 / 0.032 / 0.024 / 0.045 min^-1 for Bintje / LR / Hulda / SW, i.e.
**half-lives of 10 / 21 / 28 / 15 min**. The model-system elimination constants at 160 C are
0.1111 min^-1 (Claeys, aqueous), 0.10 min^-1 (De Vleeschouwer, a_w 0.92 powder) and 0.0881 min^-1
(Knol 2005, aqueous) — half-lives of 6.2, 6.9 and 7.9 min. **The crisp's apparent decay is 1.5x to
4.5x slower than any of them.** This is a suggestive coincidence of order and nothing more: tau
multiplies the whole formation term rather than acting on an acrylamide pool, the crisp's interior
is far below 160 C for much of the run, its a_w sweeps from 1.0 to below 0.05, and the authors
themselves say the degradation parameters carry no information. **It is recorded because it is the
only real-food number in the corpus that bears on the elimination at all, and because it points the
same way as the panel's complaint: in a food, acrylamide falls more slowly than the model-system
constants predict.**

**5. Precursor stability over the working day.** Only one of twelve morning/evening comparisons is
significant (Hulda sucrose, -24 %). Treat Table 1's averages as the charge composition.

## 4. Kinetic numbers the repository can use

**There are none.** No rate constant, no reaction order, no activation energy, no rate law. Table 2
contains shape parameters of curves fitted "without any mechanistic restraints". They are listed
below only so that nothing in them is ever mistaken for a rate.

**Registry mapping (`data/keys/compounds.yml`).** `acrylamide` is keyed. **Fructose, glucose,
sucrose and asparagine are NOT in the registry**, and neither is any potato matrix descriptor. A
charge built from Table 1 would need all four keyed plus a dry-matter basis.

| quantity | value | unit | conditions | order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| **acrylamide scale parameter a**, Logistic-Exponential | 2.6e4 / 1.4e4 / 9.6e3 / 9.3e3 (± 6e3 / 9e2 / 5e2 / 2e3) | **ug/kg dry matter** | Bintje / Lady Rosetta / Hulda / SW 91 102; 1.5 mm slices, oil 180 -> held 160 C, 1-12 min | **none — a logistic plateau, not a rate** | Table 2 p. 107 | **level_only** (a fitted asymptote, not a measurement at a time point) |
| acrylamide scale parameter a, Logistic-Fermi | 1.9e4 / 2.7e4 / 8.5e3 / 7.1e3 | ug/kg dm | same | none | Table 2 | level_only (Lady Rosetta's SD 4e4 exceeds its estimate) |
| formation inflection time tc1 | 2.8 / 2.8 / 2.3 / 3.1 (± 0.1 / 0.1 / 0.1 / 0.3) | min | same, Logistic-Exponential | none | Table 2 | level_only — **the tightest quantity in the paper**; "almost the same for all four potato genotypes" |
| formation steepness k1 | 2.8 / 3.0 / 2.5 / 1.5 (± 0.5 / 0.4 / 0.3 / 0.4) | min^-1 **as printed** | same | **NOT a rate constant** — logistic steepness | Table 2 | level_only — **must never enter a rate registry** |
| degradation characteristic time tau | 15 ± 10 / 31 ± 9 / 41 ± 14 / 22 ± 15 | min | same | none | Table 2 | level_only — SD 29-68 % of estimate; the authors call the degradation uninformative |
| Logistic-Fermi degradation k2, tc2 | k2 3.5 ± 8 / 0.053 ± 0.05 / 1.7 ± 4 / 0.63 ± 1; tc2 6.6 ± 1 / 0.60 ± 52 / 13 ± 2 / 14 ± 3 | min^-1 / min | same | none | Table 2 | **REFUSE** — SD >= estimate on most cells (this is the row `k3_final_parameter_inventory.md` §A.2 line 207 cites) |
| Empirical model a, n, b, m | see section 3 | — | same | none | Table 2 | **REFUSE** — every SD at or above its estimate |
| precursor composition, three genotypes | fructose 3.33-7.78, glucose 5.94-8.09, sucrose 5.26-11.51, **reducing sugars 9.27-15.80**, **asparagine 9.93-11.60** | **mg/g dry matter** | raw slices, cold-sweetened, morning and evening | — | Table 1 p. 105 | **measured level** — the charge composition for a hold-out benchmark |
| the same in molar terms (mine) | reducing sugars 0.054 / 0.0865 / 0.0664; asparagine 0.0795 / 0.0871 / 0.0756 | mmol/g dm | Hulda / Lady Rosetta / SW 91 102 | — | derived from Table 1 | derived_assumption (unit arithmetic only) |
| sugar : asparagine molar ratio (mine) | 0.62 / 0.99 / 0.88 | mol/mol | Hulda / Lady Rosetta / SW 91 102 | — | derived | derived_assumption |
| precursor -> acrylamide relation | a = 1.06 x C_red.sugars - 2.04; R^2 0.98 (reducing sugars), 0.97 (glucose), 0.91 (fructose) | a in 10^3 ug/kg dm (unit inferred, Flags 2); C in mg/g dm | four genotypes, one frying protocol | linear regression, not a rate law | eq 6, Results 3.3 | within_study_ratio |
| moisture trajectory | 78-80 % raw; 30-41 % at 1 min; converged by 3 min; ~5 % (SW) at long times | % w/w | during frying | — | Results 3.2 | level_only (the rest of Fig. 3 is figure-only) |
| acrylamide time shape | below 50 ug/kg dm at 1 min; rises steeply 1.5-4 min; **decreases slowly after 4 min**; a residual concentration persists at 12 min (this is why the Fermi form was rejected) | — | held 160 C | — | Results 3.2, 3.3 | **structural / level_only** |
| colour vs acrylamide | R^2 0.28 / 0.82 / 0.85 for L\* / a\* / b\* | — | Hulda only | — | Results 3.3 | within_study_ratio |
| acrylamide concentration vs time, all genotypes; temperature profiles; water content curves; colour | — | — | — | — | Figs. 1, 2, 3, 4, 5, 6, 7 | **figure_only** |

### What the repository can actually do with this paper

**(a) As a hold-out benchmark charge, with one blocking obstacle.** Table 1 gives a complete
precursor charge in mg/g dm for three genotypes, and Table 2 gives an acrylamide scale in ug/kg dm
for each. That is enough to state a benchmark. **The obstacle is water activity.** A crisp goes
from ~80 % water to ~5 % in three minutes; `acrylamide_conditions.py` declares its measured a_w
window as 0.34-0.99, refuses to extrapolate outside it, and has no moisture *trajectory* at all —
the lane takes one a_w for a whole run. A crisp is not a run at one a_w. Any benchmark built here
would have to declare a nominal a_w and carry the refusal, or the lane would need a moisture
trajectory it does not have. `results/validation/data_wishlist.md` already records that the
acrylamide lane's `moisture_aw` axis is "do-not-use ... blocked", and this paper is a good
illustration of why.

**(b) As the target of the precursor-response test.** The lane's justification for being
mass-action is that acrylamide should respond to precursor concentration. Eq 6 says that in a real
crisp the response to reducing sugar is close to proportional over 9-16 mg/g dm, while asparagine
barely varies and explains nothing. A lane charged with these three compositions should reproduce
the *ordering* Lady Rosetta > SW 91 102 > Hulda in a and should not require an asparagine term to
do it. **Note the paper's own counter-example**: SW 91 102 has more fructose and glucose than Hulda
but makes less acrylamide, which the authors attribute to unexplained "genotype specific factors".
Any lane that reproduces the regression will also reproduce that failure.

**(c) NOT as a source of constants, and not as a second laboratory's constants.** It has none.
The second-laboratory constants for this lane are in `knol2005_extraction.md` and
`knol2010_extraction.md`.

**(d) As corroboration that `k_acr_dp` should be fitted, not transcribed.** Three independent
statements now say the acrylamide degradation is unidentifiable in its own paper's fit: Knol 2005
("the model was not restrained by experimental data for the products formed in the degradation
reaction"), Knol 2010 (the barrier went negative and the step was deleted) and **this paper (every
degradation parameter's SD at or above its estimate, after the experiment was deliberately
lengthened to fix exactly that)**. The registry's GROUP 2 decision is correct and now rests on the
primary sources.

## 5. Flags

1. **This paper contains no rate constant, no reaction order and no activation energy.** Its k1 and
   k2 carry the unit min^-1 and look like rate constants; they are logistic steepness parameters.
   `k3_final_parameter_inventory.md` §A.2 line 205 already records that "Knol 2009 publishes no Ea
   at all"; confirmed. **Nothing in Table 2 may be transcribed into
   `src/kinetic_core/parameters_acrylamide.py`.**
2. **Equation 6's unit for a is not printed** and must be inferred as 10^3 ug/kg dm for the
   equation to reproduce Table 2 (mine). Its regression is over four genotypes, one of which
   (Bintje) has the highest leverage and whose sugar composition is **not in this paper** — it is
   in Knol et al. 2008, which is **not on disk**. Reconstructing it from eq 6 gives 26.5 mg/g dm
   (mine) and is circular.
3. **Table 2 is rotated on the page.** `pdftotext -layout` shreds it into fourteen fragments in
   reverse column order; the transcription above is from `pdftotext -f 5 -l 5 -raw`, which returns
   it row-wise, and the two readings were checked against each other cell by cell. Anyone
   re-extracting this table must use `-raw` or read the page image.
4. **Table 3's AICc column checks out arithmetically** (ΔAIC = AICc - min AICc, verified for all
   four genotypes, section 3), so nothing is lost there. **Table 2 is where the text layer drops
   characters**: exponents come through as bare digits ("1.9 x 10 4" for 1.9 x 10^4) and minus
   signs vanish in places. Every cell above was checked against the `-raw` extraction, which
   preserves the exponents; the units line of Table 2 reads "a (l g/kg)" in the text layer, i.e.
   **ug/kg** with the mu lost.
5. **"Hulcla" in Table 1's fructose row** is a typographical error for Hulda (every other row and
   the whole text say Hulda).
6. **Hulda's own model discrimination contradicts the paper's conclusion**: for Hulda the
   Logistic-Fermi has ΔAIC = 0.0 and the Logistic-Exponential 5.1. The overall preference for the
   Logistic-Exponential rests on the other three genotypes plus the parameter precision argument.
7. **The Bintje row is not this study's data** (Knol et al. 2008), and that paper is not on disk.
   It is included in Tables 2 and 3 without its precursor composition. Any use of the four-genotype
   band or of eq 6 inherits a source the repository does not hold.
8. **The whole run is non-isothermal, non-isobaric and drying, and the models absorb all of it.**
   Oil at 180 C falling to a held 160 C; the crisp's outer cell layer follows with a lag (Fig. 1);
   the interior is at wet-bulb temperature until the water has gone. The parameters are properties
   of *this frying protocol*, and the authors say so twice.
9. **Water activity sweeps outside the lane's declared window** (from ~1.0 to a moisture content of
   2-5 %). This blocks a straightforward benchmark charge; see section 4(a).
10. **One extraction per sample, analysed in duplicate**, and n = 2 at eight of the twelve frying
    times. The SDs in Fig. 2 are over four points only at 2, 2.5, 3 and 4 min.
11. **What this paper does NOT contain**: any activation energy; any rate constant; any reaction
    order; any mechanism; any model-system data; any absolute tabulated acrylamide concentration
    (only the fitted a); Bintje's precursor composition; any supplementary material.
12. **What to request from the authors**: (i) the acrylamide concentration-time data behind Fig. 2
    (three genotypes x twelve times, with replicates) — this is what would turn the paper from a
    parameter table into a benchmark; (ii) the crisp-surface and oil temperature profiles of Fig. 1
    as numbers; (iii) the water-content curves of Fig. 3 as numbers, which together with (ii) is
    what a moisture-trajectory lane would need; (iv) Knol et al. 2008 (Mol. Nutr. Food Res.
    52:313-321) for Bintje's composition.
13. **Registry gaps against `data/keys/compounds.yml`**: `acrylamide` present. **Fructose, glucose,
    sucrose and asparagine absent.** There is also no matrix or dry-matter-basis concept in the
    registry, and every number in this paper is per kg or per g of dry matter.
