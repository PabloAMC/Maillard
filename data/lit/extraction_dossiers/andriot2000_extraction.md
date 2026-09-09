# Andriot 2000 — EXTRACTION (three methyl ketones at 50 µL/L against 0-4 % beta-lactoglobulin in 50 mM NaCl at pH 3; static headspace at 30 C to 2700 s, a 16-panellist odour-intensity test at 21 C, and the Harrison-Hills release model; three global binding coefficients and nothing else numeric)

### THE BETA-LACTOGLOBULIN HALF OF THE BINDING BATCH, AND A PAPER WITH EXACTLY THREE USABLE NUMBERS IN IT: Kb = 330, 950 and 2 440 M^-1 for 2-heptanone, 2-octanone and 2-nonanone — no molar mass, no site count, no error bar, no table, and every other quantity in the paper locked inside a figure.

**Source on disk:** `data/articles/andriot2000.pdf` (6 pp., J. Agric. Food Chem. **2000**, 48 (9),
4246-4251).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/andriot2000.txt`). The text layer is clean modern typesetting and the
whole article — abstract, theory, methods, results, both tables and the reference list — came
through legibly; the only damage is to the mathematical display (subscripts and superscripts in
equations 1-6 are reflowed, and `K_b )` stands for `K_b =` throughout, an artefact of the journal's
"=" glyph). **The paper contains only two tables and neither holds a measurement**: Table 1 is a
symbol glossary with the model's default volumes, and Table 2 is a sensory ANOVA. **Figures 1-7 hold
every measured quantity in the study** — release curves (Fig. 1), partition coefficients (Fig. 2),
the normalised effective partition coefficient against protein concentration with the fitted binding
curves (Fig. 3), the binding affinities themselves as a bar chart (Fig. 4), a modelled release curve
(Fig. 5), mass transfer coefficients (Fig. 6) and the odour intensities (Fig. 7). All are images:
**everything in this paper except the three Kb values, the two retention ranges and the ANOVA is
figure-only.** There is no supplementary material. Repo status before this dossier: Andriot 2000 is
cited by `src/kinetic_core/parameters_matrix.py` in three `REVERSIBLE_BINDING` rows
(`kg_2_heptanone_blg`, `kg_2_octanone_blg`, `kg_2_nonanone_blg`, lines 420-437) as "andriot2000 via
k2 sec. (b) and sec. B.1", each carrying `"molar_basis": "recovered_by_arithmetic (36 800 g/mol
dimer)"`, and is named in `k2_matrix_and_thresholds.md`'s opening paragraph as one of the ten papers
whose figures were re-read at 200-900 dpi — but has **no extraction dossier**.

## 0. Identity

| field | value |
|---|---|
| Title | "Interactions between Methyl Ketones and β-Lactoglobulin: Sensory Analysis, Headspace Analysis, and Mathematical Modeling" |
| Authors | Isabelle Andriot (corresponding), Marcus Harrison, Nicole Fournier, Elisabeth Guichard — INRA Laboratoire de Recherches sur les Arômes, 17 Rue Sully, 21065 Dijon Cedex, France, **and** (Harrison) Institute of Food Research, Colney Lane, Norwich NR4 7AU, United Kingdom |
| Venue | J. Agric. Food Chem. **2000**, **48** (9), 4246-4251. Received 22 November 1999; revised and accepted 7 May 2000; published on the Web 5 August 2000. Part-financed by the French Ministry of Agriculture and Fisheries |
| DOI / article ID | 10.1021/jf991261z (printed as `JF991261Z`) |
| Protein | **commercial β-lactoglobulin**, Besnier Bridel Aliments (Châteaulin, France), **purity > 90 %**, dispersed in 50 mM NaCl and **adjusted to pH 3 with 1 N HCl** |
| Ligands | 2-heptanone, 2-octanone, 2-nonanone, gift of International Flavors and Fragrances (Longvic, France), **purity > 95 % by GC-MS**. Chosen because they are "flavor impact compounds in yogurt" (Ott 1997) |
| Quantity measured | **Kb, the "volatile-protein global binding coefficient" in M^-1**, obtained as the **fitting parameter** of eq 3 to the measured effective/actual partition-coefficient ratio against protein concentration (Fig. 3). It is a lumped constant, not an intrinsic per-site constant |
| Model | Harrison & Hills 1997 (J. Agric. Food Chem. 45, 1883-1890): penetration theory of interfacial mass transfer, binding assumed **reversible** and **not rate-limiting** |
| Relation to the repository's BLG matrix | `data/species/protein_matrices.yml` keys `blg` on the **monomer at 18 362 Da** with per-monomer counts (15 Lys, 1 free Cys, 2 S-S, 3 Arg). This paper's constants reach the code on a **36 800 g/mol** basis, i.e. the **dimer** — twice 18 400, so the two are consistent, but they are different bases and the arithmetic must not be mixed (Flags 2) |
| Second-hand content | Charles et al. 1996 (one 2-nonanone binding site per β-lactoglobulin dimer at pH 3); Sostmann & Guichard 1998 and Kinsella 1989 (bracketing binding curves drawn on Fig. 3, values not printed); Jouenne & Crouzet 1996 (agreement on retention) |

## 1. Why it matters

**(a) It is the beta-lactoglobulin leg of the reversible-binding table, and beta-lactoglobulin is a
matrix the engine actually carries.** `data/species/protein_matrices.yml` holds three matrices with
site densities: `blg`, `soy_isolate` and `pea_isolate`. Of those three, **only `blg` is built from a
sequence** — 15 lysines, 1 free cysteine, 2 disulfides and 3 arginines per 18 362 Da monomer — and
it is therefore the one matrix in the repository whose chemistry is known exactly. This paper is the
source of its three reversible constants:

| repository row | value in the code | this paper's Kb | check (mine) |
|---|---|---|---|
| `kg_2_heptanone_blg` | 8.97e-3 L/g | 330 M^-1 | 330 / 36 800 = **8.967e-3** ✓ |
| `kg_2_octanone_blg` | 2.58e-2 L/g | 950 M^-1 | 950 / 36 800 = **2.582e-2** ✓ |
| `kg_2_nonanone_blg` | 6.63e-2 L/g | 2 440 M^-1 | 2 440 / 36 800 = **6.630e-2** ✓ |

**All three reproduce with n = 1 site per 36 800 g/mol.** That confirms what "molar basis recovered
by arithmetic" means in the code: the paper prints **no** molar mass and **no** site count, and the
36 800 was reverse-engineered. Where it comes from is legible in the paper's own Introduction —
"Charles et al. (1996) have found **one binding site of 2-nonanone per β-lactoglobulin dimer at
pH 3**" — but Charles 1996 is a Weurman proceedings abstract, not this paper, and the number 36 800
appears nowhere in this article (Flags 2).

**(b) It supplies the corpus's ceiling on how much reversible binding can matter, and I reproduce
that ceiling here from first principles.** `k2_matrix_and_thresholds.md`'s four-line answer states
that "the measured ceiling on reversible hydrophobic binding is **1.3-3.7x at 4 % protein**
(Andriot) and **≤7.6x extrapolated to 10 %**", and that this is what refutes reversible binding as
the mechanism behind 100-6 700x threshold shifts. That ceiling is eq 3's factor **(1 + Kb·cb)**:

- at **4 % w/w** BLG, cb = 40 g/L ÷ 36 800 g/mol = **1.087e-3 M**, so the retention factor is
  1 + 330 x 1.087e-3 = **1.36x** (2-heptanone), **2.03x** (2-octanone), **3.65x** (2-nonanone)
  — **the printed 1.3-3.7x, reproduced (mine)**;
- at **10 %**, cb = 2.717e-3 M and the strongest binder gives 1 + 2 440 x 2.717e-3 = **7.63x** —
  **the printed ≤7.6x, reproduced (mine)**.

This is the load-bearing calculation behind `matrix_oav.py`'s `REVERSIBLE_LOG_SHIFT_CEILING` and
behind the module docstring's insistence that reversible binding is capped at about a quarter of a
log-shift. **It rests entirely on three numbers from this paper and on a molar mass this paper does
not print.**

**(c) What it does NOT do: it contains no threshold.** The odour-activity layer's refusal on protein
pots is untouched. This paper measures **odour intensity on a matched linear scale**, not detection
thresholds, and the intensities are in Fig. 7. There is no BET, no forced choice, no ASTM procedure
and no ppb value anywhere in it. It cannot add a row to `MATRIX_THRESHOLDS` for any matrix.

**(d) It is, however, the corpus's cleanest demonstration of the sensory consequence of binding.**
16 panellists, three aroma concentrations, three protein concentrations, a full three-way ANOVA
(Table 2), and a clear result: adding 1 % BLG **significantly reduces perceived odour intensity** for
all three ketones, but the effect is **weakest for 2-nonanone** — the compound that binds hardest.
The authors say so plainly: "this effect is **not well correlated** with the retention of the aromas
for the protein". That is a direct, same-study, same-panel refutation of the idea that headspace
retention translates into perception in proportion, and it is the same finding
`matrix_oav.py`'s docstring records from Baek 1999. **This paper is a second, independent witness for
that design decision, and the module does not currently cite it for that.**

## 2. Methods as they matter to a model

- **The solution, and its pH.** Commercial BLG powder dispersed in **50 mM NaCl**, **pH adjusted to 3
  with 1 N HCl**. Aroma solutions prepared **daily** in the same 50 mM NaCl at pH 3. **pH 3 is the
  single most consequential condition in this paper.** It is below every pH the engine's food pots
  run at, it is below `PH_ADDUCT_GATE_BELOW` (so no adduct chemistry would be expected), and it is
  the pH at which BLG's association state is least like its state in a neutral protein pot. **The
  paper neither measures nor discusses the protein's association state.**
- **Headspace analysis.** **40 mL amber flasks** with Supelco mininert valves. **5 mL of aqueous
  aroma solution + 5 mL of either NaCl solution or protein solution** = 10 mL liquid in a 40 mL
  flask, so **30 mL gas over 10 mL liquid** (matching Table 1's `vg = 3e-5 m3`, `va = 1e-5 m3`).
  **Final ketone concentration 50 µL L^-1** for every headspace run. **Protein at 0, 0.5, 1, 2, 3 and
  4 %** (the Methods sentence says "four final concentrations were studied" and then lists **five** —
  Flags 4). **Stirred and equilibrated at 30 C for times from 15 to 2 700 s.** "**Only one sample per
  flask was made**", and the analyses were "done in triplicate".
- **GC.** 1 mL vapour drawn with a gastight SGE syringe onto a Carlo Erba 8000 with a **DB-Wax**
  column (0.32 mm i.d., 30 m, 0.5 µm film); injector 250 C, detector 260 C; **H2 carrier at
  1.9 mL/min**; FID sampled every 50 ms by an in-house acquisition board. **No internal standard and
  no calibration procedure are described**, which is why every partition coefficient in this paper
  lives in a figure and none is printed.
- **Sensory analysis — a matched-intensity test, not a threshold test.** **16 panellists (8 male, 8
  female)**, trained in three sessions (one ketone per session) to rank **eight reference
  concentrations from 0.78 to 100 µL L^-1 in steps of 2**, then familiarised with judging intensity in
  the presence of protein by the matching test of Rousseau 1996. In three tasting sessions they
  judged **nine samples**: three aroma concentrations (**12.5, 50 and 100 µL L^-1**) crossed with
  three protein concentrations (**0, 0.5 and 1 %**), coded with three-digit numbers in a **Latin
  square** order, each ketone **in duplicate**. **20 mL samples in 60 mL brown screw-capped flasks,
  equilibrated 1 h at 21 C.** Panellists marked "their position on a continuous linear scale
  overlapping the eight reference points, reaching from <1 to >8", **judging the ketone's intensity
  "with abstraction of the protein odor"**. Three-way ANOVA on aroma concentration, BLG
  concentration and judge (random), with all interactions; means separated by Newman-Keuls.
- **Two temperatures, and they are different.** **Headspace at 30 C; sensory at 21 C.** The binding
  constants belong to the 30 C experiment. Nothing measures whether Kb differs between the two, and
  the sensory and instrumental arms are therefore not at the same temperature.
- **The model, and what it assumes.** Harrison & Hills 1997. Two assumptions stated in the
  Introduction: transport across the gas-liquid interface follows **penetration theory**, and the
  **exchange between bound and free states is fast enough not to be rate-limiting**. Binding is
  first-order (eq 2, `c_bf = Kb·c_b·c_ff`) and **reversible by assumption** — irreversible or covalent
  binding is acknowledged in the Introduction to exist (citing Hansen & Heinis 1991, 1992) and is
  then excluded from the model without test. Eq 3 gives the effective partition coefficient
  `K_ga^eff = K_ga / (1 + Kb·cb)`, with **cb taken as the total protein concentration**, which is
  what makes Kb a *global* (lumped n·K-like) constant rather than an intrinsic one.
- **How Kb was actually obtained.** Not from a Scatchard plot and not from dialysis: **the ratio
  K_ga^eff / K_ga was plotted against percentage BLG (Fig. 3) and Kb was the fitting parameter.**
  Fig. 3 also carries two comparison curves drawn from Sostmann & Guichard 1998 and Kinsella 1989,
  between which this paper's fits are said to fall. **No numerical value from either comparison
  source is printed.**
- **What is never measured.** No molar mass of the protein. No number of binding sites. No
  association state. No temperature series. No pH series (one pH, 3). No protein other than BLG. No
  aldehyde — all three ligands are methyl ketones. No error bar on Kb. No dissociation constant,
  despite the Conclusions claiming that "partition coefficients, **dissociation constants**, and mass
  transfer coefficients were quantified" (Flags 3).

## 3. Tables re-typed

### Table 1. "Description of Symbols with Corresponding Default Values Used in the Calculations"

**This is a glossary of the model, not a measurement.** Only four rows carry a value.

| symbol | description | value |
|---|---|---|
| A_ga | gas-liquid surface area (m2) | **5 x 10^-4** |
| c_a(t) | volatile concn in aq phase (mg/cm3) | (none) |
| c_g(t) | volatile concn in gaseous phase (mg/cm3) | (none) |
| c_tf(0) | initial volatile concn in aq phase (mg/cm3) | **0.041** |
| h_D | gas-liquid mass transfer coefficient (m/s) | (none) |
| K_ga | gas-liquid partition coefficient | (none) |
| K_ga^eff | effective gas-liquid partition coefficient | (none) |
| K_b | volatile-protein global binding coefficient (M^-1) | (none) |
| c_b | concn of β-lactoglobulin in aq phase (M) | (none) |
| v_g | vol of gas phase (m3) | **3 x 10^-5** |
| v_a | vol of aq phase (m3) | **1 x 10^-5** |
| t | time (s) | (none) |

### Table 2. "Results of Analysis of Variance Obtained in Sensory Analysis"

Footnote a: "*, significant at p ≤ 5 %; **, significant at p ≤ 1 %; ***, significant at p ≤ 0.1 %."

| factor | 2-heptanone F | 2-heptanone p | 2-octanone F | 2-octanone p | 2-nonanone F | 2-nonanone p |
|---|---:|---|---:|---|---:|---|
| aroma concn (A) | 102.92 | <0.0001*** | 107.36 | <0.0001*** | 113.17 | <0.0001*** |
| β-lactoglobulin concn (BLG) | 19.38 | <0.0001*** | 11.48 | <0.0001*** | **4.96** | **0.0138*** |
| judge (J) | 3.79 | <0.0001*** | 4.98 | 0.0002*** | 4.92 | <0.0001*** |
| A x BLG | 0.16 | 0.9576 | 1.03 | 0.4002 | 0.71 | 0.5852 |
| A x J | 2.08 | 0.0023** | 1.33 | 0.1348 | 1.31 | 0.1500 |
| BLG x J | 0.97 | 0.5196 | 1.21 | 0.2273 | **1.85** | **0.0091**** |
| A x BLG x J | 1.33 | 0.0887 | 0.54 | 0.9964 | 1.13 | 0.2822 |

**Read the BLG row across.** The protein effect on perceived intensity is overwhelming for
2-heptanone (F = 19.4), strong for 2-octanone (F = 11.5) and **weak for 2-nonanone (F = 4.96,
p = 0.014)** — which is the exact reverse of the headspace-retention ranking. The 2-nonanone
BLG x judge interaction is the only significant two-way interaction involving protein, and the
authors attribute it to panellists confusing the ketone's odour with the protein's own.

### Numbers printed in the running text (everything else in this paper is figure-only)

| quantity | value | where |
|---|---|---|
| **Kb, 2-heptanone** | **330** (M^-1, from Table 1's symbol definition) | Abstract and Results |
| **Kb, 2-octanone** | **950** | Abstract and Results |
| **Kb, 2-nonanone** | **2 440** | Abstract and Results |
| equilibrium time | "the equilibrium is reached after **~900 s** for the three methyl ketones" | Results |
| study duration | 2 700 s; sampling from **15 to 2 700 s** | Results, Methods |
| 2-octanone retention vs protein | "increasing the protein concentration increases the retention from **10 to 60 %**" | Results (describing Fig. 1) |
| 2-nonanone retention vs protein | "the retention is greater (**40-75 %**), reaching a **maximum at a protein concentration of 3 %**" | Results (describing Fig. 1) |
| sensory significance, 2-heptanone | 1 % protein reduces intensity significantly (p < 5 %) at **50 and 100 µL L^-1** (i.e. **not** at 12.5) | Results |
| sensory significance, 2-octanone | significant (p < 5 %) between 0 and 1 % protein at **every** aroma concentration | Results |
| sensory significance, 2-nonanone | significant **only at the lowest aroma concentration, 12.5 µL L^-1** | Results |
| binding site count | "one binding site of 2-nonanone per β-lactoglobulin **dimer** at pH 3" | Introduction — **Charles et al. 1996**, quoted |
| bracketing sources on Fig. 3 | this paper's Kb "fall between the values obtained by **Sostmann and Guichard (1998)** and **Kinsella (1989)**" | Results — **no numbers printed** |

**Partition coefficients, retention curves, binding-affinity bars, mass transfer coefficients and
odour intensities: FIGURE-ONLY.** Figures 1 (release of three ketones against 0-4 % BLG), 2
(experimentally determined partition coefficients), 3 (K_ga^eff/K_ga against % BLG with the fitted
and comparison curves), 4 (the three Kb values as a bar chart — the same three numbers as the text),
5 (modelled 2-octanone release at 2 % BLG), 6 (mass transfer coefficients) and 7 (odour intensities
with Newman-Keuls letters) are images. **Per house rule none is typed as a number.** Note that
Fig. 2's partition coefficients and Fig. 6's mass transfer coefficients exist nowhere else in the
paper: **this study measures gas-water partition coefficients for three ketones and prints not one of
them.**

### Arithmetic on the three printed constants (all mine)

**1. The retention ceiling, which is the number the repository leans on.** Retention factor
= 1 + Kb·cb, with cb = (% w/w x 10 g/L) / 36 800 g/mol:

| % BLG | cb (M) | 2-heptanone | 2-octanone | 2-nonanone |
|---:|---|---:|---:|---:|
| 0.5 | 1.359e-4 | 1.045x | 1.129x | 1.332x |
| 1 | 2.717e-4 | 1.090x | 1.258x | 1.663x |
| 2 | 5.435e-4 | 1.179x | 1.516x | 2.326x |
| 3 | 8.152e-4 | 1.269x | 1.774x | 2.989x |
| **4** | **1.087e-3** | **1.359x** | **2.033x** | **3.652x** |
| 10 (extrapolated) | 2.717e-3 | 1.897x | 3.582x | **7.631x** |

**This reproduces `k2_matrix_and_thresholds.md`'s "1.3-3.7x at 4 %" and "≤7.6x at 10 %" exactly.**
In log terms the largest effect in the table is **0.56 decades at 4 %** and **0.88 decades at 10 %**
— against threshold shifts of 100-6 700x (2 to 3.8 decades) reported for real matrices. **Reversible
hydrophobic binding, measured at its most favourable, accounts for a quarter of the observed
shift at most.** That is the whole argument, and it is arithmetic on three numbers.

**2. The chain-length rule, and how well it agrees with the soy paper.** Kb: 330 -> 950 -> 2 440,
i.e. **x2.88 then x2.57 per CH2**, geometric mean **x2.72**. Damodaran & Kinsella 1981 on soy gives
110 -> 310 -> 930 for the identical three ligands, **x2.82 then x3.00, geometric mean x2.91**.
**Two proteins, two methods (equilibrium dialysis vs headspace depletion), two laboratories, nineteen
years apart, and the per-CH2 slope agrees to 7 %.** In free energy at each paper's own temperature:
-605 cal/mol per CH2 here (30 C) against -632 cal/mol per CH2 on soy (25 C) — **a 4 % difference**.
This is the strongest cross-validation available in the binding batch and it says the *slope* is a
property of the methylene group, not of the protein.

**3. The absolute constants, protein against protein.** Per gram: 2-heptanone 8.97e-3 (BLG) vs
4.40e-3 (soy) = **2.04x**; 2-octanone 2.58e-2 vs 1.24e-2 = **2.08x**; 2-nonanone 6.63e-2 vs
3.72e-2 = **1.78x**. **Beta-lactoglobulin binds methyl ketones about twice as hard per gram as native
soy protein does, consistently across all three chain lengths.** The 1.78x on 2-nonanone is the
"agrees ... to 1.8x" that the code's note on `kg_2_nonanone_soy` claims, and it is now checkable
against both primary papers.

**4. What the sensory arm says against the headspace arm.** At 1 % BLG the headspace model predicts
retention factors of 1.09x, 1.26x and 1.66x for the three ketones — i.e. **2-nonanone should show
the largest perceptual drop**. The ANOVA shows the **smallest** protein effect for 2-nonanone
(F = 4.96 against 19.38 and 11.48) and significance at only one of three aroma concentrations.
**The instrumental and perceptual rankings are inverted.** The authors' explanation is confounding
with the protein's own odour; whatever the cause, this is a same-study, same-panel demonstration that
a headspace-calibrated free-fraction does not predict a perceptual shift — the same conclusion
`matrix_oav.py` draws from Baek 1999 and from Leksrisompong, from a third direction.

**5. A non-monotone point worth keeping.** 2-nonanone's retention "reach[es] a maximum at a protein
concentration of 3 %", i.e. it does **not** keep rising from 3 % to 4 %, whereas eq 3 is monotone in
cb by construction. **The measurement disagrees with the model at the top of the protein range**, in
the direction of less binding than predicted. The paper does not remark on it.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** **None of this paper's three ligands is
keyed**: `2_heptanone`, `2_octanone` and `2_nonanone` are all absent, though all three are used as
`compound` keys in `REVERSIBLE_BINDING` (and `2_heptanone` also appears in `matrix_oav.py`'s
`WATER_THRESHOLDS`). Every row below shares: **commercial β-lactoglobulin (>90 % pure) in 50 mM NaCl
adjusted to pH 3 with HCl; ketone at 50 µL L^-1; 5 mL aroma solution + 5 mL protein or NaCl solution
in a 40 mL amber mininert flask (10 mL liquid, 30 mL headspace, A_ga = 5e-4 m2); stirred at 30 C;
1 mL vapour to GC-FID on DB-Wax; Kb fitted to the K_ga^eff/K_ga curve of Fig. 3 through eq 3, with
c_b taken as total protein.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| 2-heptanone / β-lactoglobulin global binding coefficient | **330** | M^-1 | 30 C, pH 3, 50 mM NaCl, 0-4 % protein | Abstract, Results p. 4248; bar chart Fig. 4 | **binding_constant** (a fitted parameter, **no uncertainty given**) |
| 2-octanone / β-lactoglobulin | **950** | M^-1 | as above | Abstract, Results | **binding_constant** |
| 2-nonanone / β-lactoglobulin | **2 440** | M^-1 | as above | Abstract, Results | **binding_constant** |
| per-gram constant, 2-heptanone | 8.97e-3 | L/g | 30 C, pH 3; Kb / 36 800 g/mol | derived (mine); **already in `REVERSIBLE_BINDING`** | binding_constant — **the basis is recovered, not printed** (Flags 2) |
| per-gram constant, 2-octanone | 2.58e-2 | L/g | as above | derived (mine); already in the code | binding_constant, basis recovered |
| per-gram constant, 2-nonanone | 6.63e-2 | L/g | as above | derived (mine); already in the code | binding_constant, basis recovered |
| retention factor (1 + Kb·cb) at 4 % BLG | 1.36 / 2.03 / 3.65 | x | heptanone / octanone / nonanone, 30 C pH 3 | derived (mine) from the three Kb | within_study_ratio — **the ceiling `k2` sec. B.4 quotes** |
| retention factor at 10 % BLG | 1.90 / 3.58 / **7.63** | x | extrapolated beyond the measured range | derived (mine) | derived_assumption (**extrapolation**: 4 % is the highest protein measured) |
| chain-length effect on Kb | 2.72 | x per CH2 (geometric mean of 2.88 and 2.57) | 30 C, pH 3, C7-C9 methyl ketones | derived (mine) | within_study_ratio |
| chain-length effect in free energy | -605 | cal/mol per CH2 | 30 C | derived (mine), = -RT ln(Kb ratio) | within_study_ratio |
| BLG vs native soy, per gram | 2.04 / 2.08 / 1.78 | x | same three ligands; soy from Damodaran 1981 at 25 C by dialysis | derived (mine) | within_study_ratio — **cross-study, cross-method** |
| 2-octanone retention range, 0 -> 4 % protein | 10 to 60 | % | 30 C, pH 3 | Results p. 4248, describing Fig. 1 | level_only (**a prose reading of a figure by the authors**) |
| 2-nonanone retention range, 0 -> 4 % protein | 40 to 75 | %, **maximum at 3 %** | as above | Results, describing Fig. 1 | level_only (same caveat; and the maximum contradicts the monotone model, §3.5) |
| time to headspace equilibrium | ~900 | s | 30 C, stirred, 10 mL liquid / 30 mL gas | Results | level_only |
| sensory: effect of 1 % BLG on perceived odour intensity | F = 19.38 / 11.48 / **4.96**; p < 0.0001 / < 0.0001 / **0.0138** | — | 21 C, 16 panellists, 12.5-100 µL L^-1, 0-1 % BLG | Table 2 p. 4249 | level_only — **significance, not an effect size**; the intensities themselves are Fig. 7 |
| judge effect | F = 3.79 / 4.98 / 4.92, all p ≤ 0.0002 | — | as above | Table 2 | level_only |
| model default geometry | A_ga 5e-4 m2; v_g 3e-5 m3; v_a 1e-5 m3; c_tf(0) 0.041 mg/cm3 | — | the headspace flask | Table 1 p. 4247 | level_only |
| gas-water partition coefficients of the three ketones | — | — | pure water, 30 C | Fig. 2 | **figure_only** — measured and never printed |
| mass transfer coefficients h_D | — | m/s | 0-4 % BLG | Fig. 6 | **figure_only** |
| odour intensities with Newman-Keuls letters | — | 1-8 matched scale | 21 C, 9 sample combinations | Fig. 7 | **figure_only** |
| release curves, all ketones and all protein levels | — | headspace concn vs time | 15-2 700 s, 30 C | Figs. 1, 3, 5 | **figure_only** |
| binding site count, 1 per BLG dimer at pH 3 | 1 | site/dimer | pH 3 | Introduction — **Charles et al. 1996** | level_only, **second-hand**; the load-bearing assumption behind the 36 800 basis |

### What can and cannot be put next to these

**(a) The three constants are usable and are correctly transcribed into the code — but their per-gram
form carries an assumption the paper does not license.** Kb = 330 / 950 / 2 440 M^-1 is printed and
unambiguous. Turning it into L/g requires dividing by a molar mass, and **this paper prints none**.
The 36 800 g/mol used in `parameters_matrix.py` is the BLG dimer, justified by Charles 1996's "one
site per dimer at pH 3" quoted in this paper's Introduction. If the correct basis were the **monomer**
(18 362 Da, which is what `protein_matrices.yml` uses for `blg`), **every per-gram constant would
double** and so would the retention ceiling. **A factor-of-two risk sits under three FIT rows and
under the corpus's headline ceiling number.**

**(b) They are pH 3 constants, and nothing in this corpus tells you how they move.** The engine cooks
protein pots at pH 6-7. Beta-lactoglobulin's conformation and association state are strongly
pH-dependent; the binding site invoked here (the central calyx, per Wu 1999 on palmitate, cited in
the Introduction) is gated by the Tanford transition around pH 7. **No pH series exists in this
paper or, for this protein and these ligands, anywhere else in the batch.** Carrying these constants
to a neutral pot is an extrapolation and should be marked as one.

**(c) They are 30 C constants, and Damodaran gives the only licence to move them.** Damodaran &
Kinsella 1981 measured that a soy-ketone binding constant is unchanged between 25 and 45 C
(ΔH ≈ 0, entropy-driven hydrophobic association). If that carries to BLG — untested — then 30 C is
transportable across the same window. Above 45 C, nothing.

**(d) They are methyl ketones, and the aldehyde question is untouched.** All three ligands are
ketones. The engine's aldehyde species (`HEXANAL`, `NONANAL`, `DECADIENAL`, `FUR`,
`ME_9_OXONONANOATE`, `ME_13_OXO_TRIDECADIENOATE`) bind through `matrix_sites.py`'s covalent classes,
which is a rate and not an equilibrium constant, and this paper says nothing about them. Note that
the corpus's aldehyde-on-BLG rate bracket
(`saturated_aldehyde_amine`, 1e-6 to 1e-5 M^-1 s^-1 from anantharamkrishnan2020b) is a different
kind of quantity entirely and is not comparable to a Kb.

**(e) The sensory arm is a qualitative result and should be cited as one.** The three F values in
Table 2's BLG row rank the protein effect **opposite** to the retention ranking. That is usable as
evidence for a design decision; it is not usable as a number, because the intensities are in Fig. 7.

## 5. Flags

1. **Three numbers, no uncertainty, and they are fitting parameters.** Kb = 330, 950 and 2 440 M^-1
   are the fitted values of one parameter in eq 3 against the Fig. 3 curve of K_ga^eff/K_ga versus
   percentage protein. **There is no confidence interval, no standard error, no goodness-of-fit
   statistic and no replicate spread on any of them.** The headspace runs are described as triplicate
   with "only one sample per flask", so replication exists but is never propagated into Kb. Compare
   the constants' precision as quoted (three significant figures on 2 440) with the fact that the
   underlying plot has six protein levels and no error bars.
2. **The molar basis is recovered, not printed, and it is a factor-of-two risk.** The paper gives no
   molar mass for β-lactoglobulin anywhere, yet eq 3 needs c_b **in M** and the Methods state protein
   only in **% w/w**. So the authors themselves used a molar mass to produce Fig. 3, and did not say
   which. The repository's 36 800 g/mol is the dimer and reproduces the code's three L/g values
   exactly, and it is consistent with the Charles 1996 sentence in the Introduction ("one binding
   site per **dimer** at pH 3"). **But `protein_matrices.yml` keys the same protein on the 18 362 Da
   monomer**, and if the paper's `c_b` were monomer-molar the per-gram constants would all double.
   **Requesting the authors' `c_b` conversion is the single highest-value clarification in this
   dossier.** Note also that β-lactoglobulin's monomer-dimer equilibrium is itself pH-dependent, and
   pH 3 is at the edge of the range where the dimer is the dominant species; the paper measures
   nothing about this and neither does anything else in the corpus.
3. **The Conclusions claim a quantity the paper never reports.** "Using this model, partition
   coefficients, **dissociation constants**, and mass transfer coefficients were quantified." **No
   dissociation constant appears anywhere in the article**, in any table, in any figure caption or in
   the text. (A dissociation constant would be 1/Kb = 3.03, 1.05 and 0.41 mM — derivable, but not
   what the paper says it reports.) Read "dissociation constants" as a loose reference to Kb.
4. **The protein-concentration list does not match its own count.** Methods: "For protein, **four**
   final concentrations were studied (**0.5, 1, 2, 3, and 4 %**)" — five values. The abstract says
   "different concentrations of β-lactoglobulin (**0, 0.5, 1, 2, 3, and 4 %**)" — six, including the
   zero control. **Read as five protein levels plus a zero control**; the "four" is an error.
5. **Every measured quantity except three constants is in a figure.** This paper measures gas-water
   partition coefficients for three ketones (Fig. 2), release curves at six protein levels over
   2 700 s (Fig. 1), mass transfer coefficients (Fig. 6) and odour intensities for nine sample
   combinations (Fig. 7), and **prints none of them**. The two retention ranges quoted in the Results
   ("10 to 60 %", "40-75 %") are the authors' own prose readings of Fig. 1 and are the only numeric
   trace of the release data. `k2_matrix_and_thresholds.md` records that Andriot's Figures 1-7 were
   re-read at 200-900 dpi; **any value in the repository sourced to those figures is a figure-read and
   must stay classified `figure_only`**, not promoted to a measured number by having been read
   carefully.
6. **The sensory and instrumental arms disagree, at different temperatures, and the paper concedes
   it.** Retention ranks nonanone > octanone > heptanone; the perceptual protein effect ranks
   heptanone > octanone > nonanone. The authors write that the effect "is **not well correlated** with
   the retention of the aromas" and blame the protein's own odour, supported by the significant
   BLG x judge interaction for 2-nonanone. Add that the two arms were run at **30 C and 21 C**
   respectively, and in **different vessels** (40 mL mininert flask with 10 mL liquid vs 60 mL brown
   flask with 20 mL). **The two arms are not a matched comparison** and the disagreement cannot be
   fully attributed to perception.
7. **The model's assumption of reversibility is untested here.** Eq 2 assumes first-order reversible
   binding; the Introduction acknowledges that proteins bind volatiles irreversibly as well (Hansen &
   Heinis 1991, 1992) and then does not test for it. For methyl ketones at pH 3 the assumption is
   probably safe — ketones are poor Schiff-base formers and pH 3 disfavours amine chemistry — but it
   means **Kb is an upper bound on reversible binding only if none of the loss was irreversible**,
   and the experiment cannot distinguish the two. Contrast `parameters_matrix.py`'s treatment of
   Meynier's t-2-hexenal row, which is quarantined precisely because 22-33 % of the apparent
   "partition" is irreversible chemistry.
8. **The fit is admitted to be imperfect and the reason given is experimental, not model error.**
   "On first inspection of Figure 5 the correlation between theory and experiment is not ideal ...
   Deviations from the desired theoretical curves are a consequence of **inefficient stirring** in
   the aqueous phase." So the release-curve modelling — the third leg of the paper — rests on a fit
   the authors describe as not ideal, with the discrepancy attributed to a stirring artefact that is
   not characterised. **This affects the mass transfer coefficients (Fig. 6), not the Kb values,
   which come from equilibrium data.**
9. **A non-monotone retention point at the top of the protein range.** 2-nonanone's retention peaks
   at **3 %** protein and does not increase to 4 %, while eq 3 is monotone in protein concentration.
   Unremarked by the authors. It is a single point read from a figure, so it is weak evidence — but
   it is evidence in the direction of the model over-predicting binding at high protein, which is the
   regime a real protein pot sits in.
10. **What this paper does not contain**: any odour **threshold** of any kind; any pH other than 3;
    any temperature series; any protein other than β-lactoglobulin; any aldehyde; any molar mass;
    any binding-site count of its own; any uncertainty on any constant; any tabulated concentration,
    partition coefficient or intensity; any supplementary material.
11. **What to request from the authors**: (i) the molar mass used to convert % w/w to c_b — this is
    the factor-of-two question and it is one sentence; (ii) the numeric data behind Figures 1, 2, 6
    and 7, in particular the three gas-water partition coefficients, which are measured here and
    printed nowhere; (iii) an uncertainty on each Kb; (iv) whether the same experiment has been done
    at a food-relevant pH.
12. **Registry gaps against `data/keys/compounds.yml`**: **`2_heptanone`, `2_octanone` and
    `2_nonanone` are all absent from the registry**, and all three are used as `compound` keys in
    `REVERSIBLE_BINDING` for both the BLG and the soy rows — six FIT rows keyed on three names the
    compound registry does not know. This is the same finding as in `damodaran1981_extraction.md`
    Flags 12, from the other side of the cross-validation, and it is now the most repeated registry
    gap in the binding batch.
