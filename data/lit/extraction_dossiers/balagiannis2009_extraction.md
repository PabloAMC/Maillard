# Balagiannis 2009 — EXTRACTION (aqueous extract of defatted ox liver, 1 : 1 with water, no buffer, pH not stated, 20 mL in sealed 30-mL ampoules, 120 / 130 / 140 C, 5-240 min; glucose, twenty free amino acids and 2- and 3-methylbutanal followed; a six-equation multiresponse model with a pseudo-first-order glucose step, printed at 130 C with 95 % HPD intervals and three barriers)
### The only paper on disk that fits leucine- and isoleucine-derived aldehydes at 120-140 C with printed constants and intervals — and its central finding is that in a food matrix the Strecker aldehyde rate does not depend on the free amino-acid concentration at all, so no per-amino-acid second-order constant exists in it by construction.

**Source on disk:** `data/articles/balagiannis2009.pdf` (7 pp., 970,803 bytes, owner's download
2026-09-09). Born-digital ACS PDF, clean text layer read from `scratchpad/articles/balagiannis2009.txt`;
pages 9918-9920 rendered at 110 dpi to read Figure 1 (axes and units), Figure 2 (regression equations
printed on the plot), Figure 3 (scheme), Figure 4 (the six rate equations, an image) and Table 1 (in the
text layer, checked against the raster cell by cell). Table 1 re-typed in full below. Figures 1 and 5
(time courses and fits) are figure-only; Figure 2's two regression equations are printed text on the
figure and are recorded as such. No supplementary information exists for this paper. The companion
chapter on beef muscle + ribose is on disk as `balagiannis2010_extraction.md`; the earlier, superseded
version of this model (one liver batch, second-order first step) is a Weurman 2008 proceedings paper
(ref 23), not on disk.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetic Modeling of the Generation of 2- and 3-Methylbutanal in a Heated Extract of Beef Liver" |
| Authors | Dimitrios P. Balagiannis, Jane K. Parker*, D. Leo Pyle (deceased 25 March 2008), Neil Desforges, Bronislaw L. Wedzicha, Donald S. Mottram — University of Reading; Waltham Centre for Pet Nutrition; University of Leeds |
| Journal | J. Agric. Food Chem. 57 (21), 9916-9922 (2009); received 30 April 2009, revised 23 July 2009, accepted 10 September 2009, web 9 October 2009 |
| DOI | 10.1021/jf901443m (printed p. 9916) |
| PDF file name | `data/articles/balagiannis2009.pdf` |
| On disk vs SI | full paper; no SI |
| Naming | Int1, Int2 = the two lumped intermediate pools; M = other Maillard products formed with amino acids; M' = products formed from Int1 without amino acids (protein binding, organic acids); R_leu, R_ile = ratio of leucine / isoleucine to total free amino acids (fixed from Figure 2); F_leu, F_ile = fitted fraction of reacted leucine / isoleucine that becomes the aldehyde; "3MeBut", "2MeBut" in the equations; HPD = highest posterior density |
| Funding | BBSRC studentship with Waltham (Mars Petcare) |

## 1. Why it matters

The amino-acid-identity wave drafted in `results/validation/kinetic_core_b19_prereg_draft.md` names this
paper for "leucine and isoleucine at 120 to 140 C". It is the only source on disk with printed rate
constants and 95 % intervals for the methylbutanals in that window. But the constants are not what the
draft's row asks for. The model's first step is **pseudo-first order in glucose** with the amine
concentration deliberately removed, because across three liver batches whose free amino acids differed
three- to four-fold the aldehyde yield tracked glucose and not the amino acids; the Strecker step itself
(Int2 + leucine) is declared fast and never identified; leucine and isoleucine enter only as fixed
ratios to the total amino acid pool (0.111, 0.034) times fitted conversion fractions (2.3 %, 3.8 %).
What the paper gives the trunk is therefore: (a) a glucose-consumption constant and barrier in a
protein-rich food matrix at 120-140 C — k1 = 1.36e-2 min-1 at 130 C, Ea1 137 +/- 15 kJ/mol — to set
beside `MARTINS_M4`'s second-order `k_schiff` (1.6e-5 L/(mmol min) at 100 C, Ea 96.8 +/- 2.8, glucose +
glycine 200 + 200 mmol/L, pH 6.8); (b) the leucine : isoleucine selectivity per unit amino acid in one
pot (F_ile / F_leu = 1.6), a within-study ratio the wave can use as an identity factor; (c) a measured
first-order sink for 3-methylbutanal (k3 = 2.7e-3 min-1 at 130 C, Ea3 78 +/- 39 kJ/mol), which the trunk
has nowhere else; and (d) a warning: the aldehyde yield per glucose in a food matrix is 0.13 % (3-MB)
and 0.065 % (2-MB) (mine, section 4), two orders below what a clean pot would give, because half of the
intermediate goes to protein and most of the rest to melanoidin-type products. For the pyrazine step's
Strecker constants (`FROZEN_B18`, second order in dicarbonyl and glycine, Ea 103.1 and 114.9 kJ/mol) the
paper offers no test: its Int2 -> aldehyde step has no constant.

## 2. Methods as they matter to a model

- **Matrix.** Sliced ox liver (four retail batches, different dates) mixed with an equal mass of
  deionised water, homogenised 1 min, centrifuged 20 min at 29,800 g at 4 C, supernatant filtered
  (Whatman no. 3). **No buffer; pH never stated or measured in the text** (Flag 1). Concentrations are
  reported in **mmol/kg** of extract. Glucose (Fig. 1a, figure-only) is of order 1e2 mmol/kg at t = 0 and
  differs by batch ("batch 3 having 50 % more glucose than batch 1"); total free amino acids (Fig. 1b)
  differ 3-4-fold across batches ("Batch 2 had levels of amino acids 3-4 times higher than the other
  two"); "the ratio of free amino acids to sugars was of the order of 1 : 1"; protein and peptides
  present at unquantified levels ("an abundance of other amino groups, effectively present in excess").
- **Heating.** 20 mL aliquots sealed in 30 mL glass ampoules (about 10 mL headspace, air), oil bath,
  5-240 min. Batches 1, 2, 3 at 130 C; batch 4 at 120 C; batch 2 also at 140 C. Heat-up at 130 C: "2 min
  to reach 100 C and 4 min to reach 120 C". At least two replicate ampoules per time. Quench in dry
  ice / methanol at -50 C. Time zero convention (immersion vs temperature reached) not stated.
- **Volatiles.** 5 g heated extract + 10 mL water in a 250 mL flask with a Dreschel head, 60 C water
  bath, N2 40 mL/min for 1 h onto Tenax, purge 100 mL/min 10 min; ATD-GC-MS (Perkin-Elmer Clarus 500,
  DB-5 60 m x 0.32 mm x 1 µm), LRI on C6-C25 alkanes, identification by spectra and LRI against
  authentic compounds. **Quantification of 2- and 3-methylbutanal: "calibration curves using the
  standard addition method and headspace collected as described above"** (spiking levels not printed;
  1,2-dichlorobenzene is in the materials list and is presumably the internal standard, but the text
  never says so). Units mmol/kg (Fig. 1e,f ordinates).
- **Free amino acids.** 0.5 g homogenate + 10 mL 0.01 M HCl, 15 min, settle, centrifuge 7200 g 30 min;
  EZ-Faast derivatisation, GC-MS (Agilent 5975); L-norvaline (materials list) is the EZ-Faast internal
  standard. Twenty amino acids detected.
- **Sugars.** Dionex ion chromatography, CarboPac PA10, 96 : 4 water : 400 mM NaOH isocratic 30 min,
  pulsed amperometry; internal standard "trehalose solution (20 µg L-1)" — an implausibly low
  concentration for an IC internal standard, probably a misprint for mg L-1 (Flag 7); standards glucose,
  fructose, sucrose, ribose, maltose, mannose. Fructose and mannose rise to maxima of **8.3 and 3.7
  mmol/kg** and fall; omitted from the model.
- **Model (Figure 3, equations of Figure 4, re-typed):**
  (1) d[Glu]/dt = -k1 [Glu]
  (2) d[Int1]/dt = k1 [Glu] - k2 [Int1] R_leu F_leu - k2 [Int1] R_ile F_ile - k2 [Int1] (1 - R_leu F_leu - R_ile F_ile) - k2 [Int1]
  (3) d[3MeBut]/dt = k2 [Int1] R_leu F_leu - k3 [3MeBut]
  (4) d[2MeBut]/dt = k2 [Int1] R_ile F_ile - k4 [2MeBut]
  (5) d[M]/dt = k2 [Int1] (1 - R_leu F_leu - R_ile F_ile)
  (6) d[M']/dt = k2 [Int1]
  Note (mine): the four consumption terms in (2) sum to 2 k2 [Int1]; the M' path (glucose-derived
  intermediates binding protein, without amino acids) takes half of Int1 and the amino-acid paths the
  other half. R_leu = 0.111 and R_ile = 0.034 are fixed from Figure 2, not fitted. Arrhenius
  reparameterised per Brands & van Boekel 2002 (reference temperature not stated; constants printed at
  130 C). Ten parameters fitted on five data sets (three batches at 130 C, batch 2 at 140 C, batch 4 at
  120 C) in Athena Visual Studio. Alternative mechanisms were tried; the second-order first step of the
  2008 version "ceased to give an accurate prediction" across batches.
- **Figure 2 (printed regression text):** leucine vs total free amino acids: y = 0.1109 x - 0.9323, R2 =
  0.9588; isoleucine: y = 0.0342 x - 0.0415, R2 = 0.9655 (raw and heated extracts pooled; abscissa 0-200
  mmol/kg, figure-only apart from the equations).

## 3. Tables re-typed

### Table 1. "Optimal Estimates for the Parameters in the Revised Model (Figure 3)" (p. 9920)

| parameter | optimal estimate (a) | percentage CI as printed |
|---|---|---:|
| k1 (b) (min-1) | 1.36 x 10-2 +/- 9.22 x 10-4 | 7 % |
| k2 (min-1) | 5.77 x 10-2 +/- 2.88 x 10-2 | 50 % |
| k3 (min-1) | 2.72 x 10-3 +/- 1.01 x 10-3 | 37 % |
| k4 (min-1) | 2.56 x 10-4 +/- 8.75 x 10-4 | 342 % |
| F_leu | 2.33 x 10-2 +/- 2.45 x 10-3 | 11 % |
| F_ile | 3.83 x 10-2 +/- 4.28 x 10-3 | 11 % |
| Ea1 (kJ mol-1) | 1.37 x 10^2 +/- 1.52 x 10^1 | 11 % |
| Ea2 (kJ mol-1) | 4.87 x 10^1 +/- 8.40 x 10^1 | 173 % |
| Ea3 (kJ mol-1) | 7.82 x 10^1 +/- 3.90 x 10^1 | 50 % |
| Ea4 (kJ mol-1) | not determined | — |

(a) Optimal values and their 95 % higher posterior density intervals; values in parentheses are the
percentage confidence intervals. (b) Rate constants are at 130 C.

### Numbers in the running text

- Fructose and mannose maxima 8.3 and 3.7 mmol/kg (130 C, batch not stated).
- R_leu = 0.111, R_ile = 0.034 (R2 0.959, 0.966).
- "< 5 % of the reacted leucine and isoleucine was converted into methylbutanals".
- "about 80 % of the glucose remains unaccounted for" by amino-acid loss (the reason for the M' step).
- Amino acids: "rapid consumption" in the first 20 min, then level (regeneration from Amadori
  breakdown or proteolysis suggested).
- 3-Methylbutanal levels "started to decrease" after the maximum; 2-methylbutanal's decrease "not
  obvious".
- Comparators quoted: Jousse 2002 Ea 66.5 kJ/mol for the analogous second step; Cremer & Eichner 2000
  "120 and 124 kJ mol-1" for 3- and 2-methylbutanal (swapped against the source, see
  `cremer2000_extraction.md` Flag 4); Chan & Reineccius 1994, 3-MB 80.4 kJ/mol.

### Figures 1 and 5 — FIGURE-ONLY

Fig. 1 (130 C, batches 1-3): glucose (0-160 mmol/kg axis) falling to near zero by 240 min; total amino
acids (0-200); leucine (0-18); isoleucine (0-6); 3-methylbutanal (0-0.14 mmol/kg) and 2-methylbutanal
(0-0.10 mmol/kg) with a 15-20 min lag, a rise to a plateau near 90-180 min, 3-MB then falling. Fig. 5:
overlay of data and model for glucose, 3-MB, 2-MB (five panels). No value typed.

## 4. Kinetic numbers the repository can use

Registry mapping: 3-methylbutanal -> `3_methylbutanal`; 2-methylbutanal -> `2_methylbutanal`;
**glucose, fructose, mannose, leucine, isoleucine, trehalose, norvaline, 1,2-dichlorobenzene -> not in
`data/keys/compounds.yml`.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| glucose -> Int1 (all amine sources incl. protein; Schiff/Amadori lumped) | k1 at 130 C | 1.36e-2 +/- 9.2e-4 | min-1 | ox-liver extract 1 : 1 water, no buffer, pH unstated, glucose ~1e2 mmol/kg (figure-only), free amino acids ~1 : 1 with glucose plus protein, sealed ampoules with air, 120-140 C | pseudo-first order in glucose (amine term removed) | Table 1 | measured_rate (multiresponse fit; hides the whole amine pool) |
| same | Ea1 | 137 +/- 15 | kJ/mol | 120 / 130 / 140 C, five data sets | Arrhenius on k1 | Table 1 | measured_barrier |
| Int1 -> Int2 (Amadori-like -> dicarbonyl pool) | k2 at 130 C | 5.77e-2 +/- 2.9e-2 | min-1 | same | first order in Int1 (Int1 leaves at 2 k2 in total, eq 2) | Table 1 | measured_rate (model-identified; 50 % CI) |
| same | Ea2 | 48.7 +/- 84.0 | kJ/mol | same | — | Table 1 | measured_barrier in name only: the interval includes zero; do not use |
| 3-methylbutanal -> other products | k3 at 130 C | 2.72e-3 +/- 1.0e-3 (half-life 255 min, mine) | min-1 | same | first order | Table 1 | measured_rate (sink) |
| same | Ea3 | 78.2 +/- 39.0 | kJ/mol | same | — | Table 1 | measured_barrier (50 % CI) |
| 2-methylbutanal -> other products | k4 at 130 C | 2.56e-4 +/- 8.75e-4 | min-1 | same | first order | Table 1 | measured_rate, not distinguishable from zero |
| Int2 + leucine -> 3-methylbutanal | k | **none: declared fast, diffusion-controlled, not fitted** | — | — | — | text p. 9919 | — |
| leucine share of the amine pool | R_leu | 0.111 (regression 0.1109, R2 0.959) | mol/mol total free amino acids | raw and heated extracts, four batches | fixed input | Fig. 2 (printed equation), text | level_only |
| isoleucine share | R_ile | 0.034 (0.0342, R2 0.966) | mol/mol | same | fixed input | Fig. 2, text | level_only |
| fraction of reacted leucine -> 3-MB | F_leu | 2.33e-2 +/- 2.5e-3 | — | same | fitted constant | Table 1 | measured_rate-class fit parameter (yield fraction) |
| fraction of reacted isoleucine -> 2-MB | F_ile | 3.83e-2 +/- 4.3e-3 | — | same | fitted | Table 1 | as above |
| amino-acid identity, per unit amino acid | F_ile / F_leu | 1.64 (+/- ~16 %, mine from the two 11 % CIs) | — | one pot, 120-140 C | — | derived from Table 1 | within_study_ratio |
| product ratio in the pot | (R_leu F_leu) / (R_ile F_ile) = 3-MB : 2-MB formation rate | 2.586e-3 / 1.302e-3 = 1.99 | — | same | — | derived (mine) | within_study_ratio |
| molar yield per glucose consumed | 3-MB: R_leu F_leu / 2 = 1.29e-3 (0.13 %); 2-MB: 6.5e-4 (0.065 %) | mol/mol glucose | same; the /2 is the M' half of eq 2 | — | derived from eq 2 + Table 1 (mine) | derived_assumption |
| maximum 3-MB formation rate per glucose (quasi-steady Int1) | k1 R_leu F_leu / 2 = 1.76e-5 x [Glu] | min-1 x mmol/kg | 130 C | — | derived (mine) | derived_assumption |
| k1 extrapolated with Ea1 (mine) | 4.8e-3 (120 C); 5.1e-4 (100 C); 3.7e-2 (140 C) | min-1 | outside 120-140 C an extrapolation | — | derived from Table 1 | derived_assumption |
| fructose, mannose maxima | 8.3, 3.7 | mmol/kg | 130 C | — | text | level_only |
| glucose, amino acids, aldehydes vs time; model overlays | — | mmol/kg | 120-140 C | — | Figs. 1, 5 | figure_only |

**Can a per-amino-acid second-order constant in water be derived?** No, and the paper's argument is
the reason the wave should hear: in a matrix where the amine pool is in excess (free amino acids ~1 : 1
with glucose, plus protein), the aldehyde rate is set by glucose -> intermediate -> dicarbonyl and the
amino acid only partitions the dicarbonyl. The fitted quantities that carry amino-acid identity are
F_leu and F_ile, yield fractions with the fast step's constant divided out. The one identity number
that transports is the ratio F_ile / F_leu = 1.6: per mole of amino acid present, isoleucine is
converted to its aldehyde 1.6 times more efficiently than leucine in this pot (Huang 2017's aqueous
pots find leucine degraded about twice as fast as isoleucine, a different quantity — loss, not
aldehyde yield — and a different sugar ratio). To turn k1 into a second-order constant one would divide
by an amine molarity; the free amino acid total is figure-only (Fig. 1b, tens to ~180 mmol/kg by batch)
and the protein amine is unquantified, so the division is not licensed. For orientation only (mine,
derived_assumption): Martins' `k_schiff` 1.6e-5 L/(mmol min) at 100 C over 200 mmol/L glycine is a
pseudo-first-order 3.2e-3 min-1 in glucose; Balagiannis' k1 at 100 C by the printed barrier is
5.1e-4 min-1, a factor 6 lower, i.e. what Martins' constant would give over about 30 mmol/L of
glycine-equivalent amine at pH 6.8 — the same order as a liver extract's free amino acids. The two are
compatible in magnitude; their barriers are not: **Ea1 137 +/- 15 vs Martins' 96.8 +/- 2.8 kJ/mol, no
overlap** (Flag 3).

**Comparison with the trunk's Strecker constants.** `FROZEN_B18` measures the dicarbonyl + amino acid
step at 103.1 / 114.9 kJ/mol on fed glyoxal / methylglyoxal at 20 + 20 mmol/L; this paper's
corresponding step has no constant and its Ea2 (Int1 -> Int2, 48.7 +/- 84) is uninformative. The
paper's own reading — k1 and k2 rate-limiting, the Strecker step "diffusion controlled" — is the same
structure as the trunk's (supply-limited Strecker), with the difference that in the trunk at 200 + 200
mmol/L the Strecker step is second order in the amine and here the amine is saturating.

## 5. Flags

1. **pH is never stated.** No buffer, no measurement, no mention. A liver extract is typically mildly
   acidic, but nothing here says so; no row from this paper can carry a pH, and the trunk's
   `check_ph_homogeneity` would refuse it. Ask the authors for the initial and final pH of the extracts.
2. **The amine concentration is deliberately absent from k1.** k1 is a rate in a matrix, not a
   constant; it hides glucose ~1e2 mmol/kg (figure-only), free amino acids ~1 : 1 with glucose (text),
   and protein. It transports to other liver extracts of similar composition and nowhere else, as the
   authors say ("system dependent ... not yet been challenged").
3. **Ea1 137 +/- 15 kJ/mol against Martins' 96.8 +/- 2.8 for the same nominal step in water.** The
   intervals do not overlap. Possible reasons: (i) k1 here lumps Schiff, Amadori and the glucose ->
   protein path; (ii) three temperatures and a 120-140 C window vs Martins' 80-120; (iii) the heat-up (4
   min to 120 C at a 130 C bath) is inside the 5-min first point and biases the 130 and 140 C early
   points differently. Record the conflict; do not average.
4. **Ea2 is not a number.** 48.7 +/- 84.0 kJ/mol; the authors compare it to Jousse's 66.5 anyway. k2 is
   50 % uncertain. Neither should enter a fit as a prior.
5. **F_leu and F_ile are small (2.3 %, 3.8 %) and the M' path takes half of Int1** — the aldehyde yield
   per glucose is 0.13 % / 0.065 % (mine). A clean glucose + leucine pot will not reproduce these; they
   are a food-matrix statement. Cremer & Eichner's dry glucose + leucine pot conserves leucine as Leu +
   Fru-Leu + 3-MB for two hours — the opposite regime.
6. **Time-zero and heat-up.** 20 mL ampoules take 4 min to reach 120 C in a 130 C bath; the time axis
   starts at immersion (presumably). The 5-20 min "lag" in the aldehyde curves is partly heat-up, partly
   the Int1 build-up the model attributes it to.
7. **"Trehalose solution (20 µg L-1)"** as IC internal standard is almost certainly a misprint (a
   20 µg/L spike into 100 µL of sample is invisible on PAD); read as mg/L. Does not affect any kinetic
   number here.
8. **Standard-addition levels for the aldehydes are not printed**; no LOD; the methylbutanal ordinates
   run to 0.14 and 0.10 mmol/kg, so the early points (< 0.01 mmol/kg) are near the floor.
9. **Headspace contains air** (sealed ampoules, no purge). Irrelevant for the methylbutanals, relevant
   if anyone later uses this paper for sulfur volatiles.
10. **Concentrations are per kg of extract (1 : 1 liver : water)**, not per litre and not per kg of
    liver; a comparison with a mmol/L pot needs the extract density (~1.02) and is otherwise direct.
11. **What to request from the authors:** the batch-wise initial glucose and total / individual free
    amino acid concentrations (Fig. 1 data), the extract pH, the standard-addition levels, the
    reference temperature of the reparameterised Arrhenius form, and the raw time courses at 120 and
    140 C (Fig. 5 only).
12. **Registry gaps** (`data/keys/compounds.yml`): both aldehydes are keyed; glucose, leucine,
    isoleucine, fructose, mannose and the internal standards are not.
13. **Quotation errors in the paper:** Cremer & Eichner's 120 / 124 kJ/mol are assigned to the wrong
    aldehydes (source: 2-MB 120, 3-MB 124); ref 27 prints "2-acetyl-1-pyrrolone" for pyrroline and
    omits Reineccius from the RSC editors.
