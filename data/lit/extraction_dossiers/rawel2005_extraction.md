# Rawel, Meidtner & Kroll 2005 — EXTRACTION (six plant phenolics against BSA, HSA, soy glycinin and lysozyme by Hummel-Dreyer size-exclusion chromatography and quercetin-fluorescence enhancement, 0.05 M Hepes or Na-acetate, pH 3.5-7.4, 25 C; thirteen dissociation constants, seven site counts and six free energies, with pH, ionic-strength and thermal series on BSA)

### A DIFFERENT LIGAND CLASS ON THE SAME PROTEIN SITES: phenolics bind soy glycinin and the serum albumins two to three orders of magnitude harder than the methyl ketones of the binding batch do, and every one of the paper's own knobs — lower pH, more salt, more heat — makes the binding WEAKER, which is the opposite of what Damodaran measured for a carbonyl on soy.

**Source on disk:** `data/articles/jf0480290.pdf` (8 pp., J. Agric. Food Chem. **2005**, 53 (10),
4228-4235). **The filename is the ACS article ID, not an author-year key**; the dossier is filed as
`rawel2005_extraction.md`.
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/jf0480290.txt`), which is clean modern typesetting; all four tables came
through complete. **Page 3 (journal p. 4230) was additionally rendered at 200 dpi and read from the
image**, to settle the unit convention on the `K_D x 10^-6 (M)` column and the apparent
factor-of-1000 conflict between Table 1 and the Figure 1 caption — see Flags 1, which is resolved
from the figure's own axis. Figures 1-5 are images: **the chlorogenic-acid binding isotherm (Fig. 1),
the ferulic-acid binding capacities across eight proteins (Fig. 2), the quercetin-BSA fluorescence
fit (Fig. 3), the urea reversibility demonstration (Fig. 4) and the near-UV CD spectra of soy
glycinin (Fig. 5) are figure-only.** **Figure 2 is the paper's only measurement of ferulic-acid
binding to soy glycinin, gelatin and whey, and it is a bar chart — those numbers are not printed
anywhere.** There is a substantial **Supporting Information (parts I-IV)** covering the Hummel-Dreyer
principle, the chromatographic conditions, the derivation of the free-energy equation, the
fluorescence spectra, the urea experiment, the tryptophan quenching and the CD assignments;
**it is NOT on disk** and several statements in the paper are unverifiable without it (Flags 3).
Repo status before this dossier: **Rawel 2005 is not cited anywhere in `src/`** and has no
extraction dossier.

## 0. Identity

| field | value |
|---|---|
| Title | "Binding of Selected Phenolic Compounds to Proteins" |
| Authors | Harshadrai M. Rawel (corresponding), Karina Meidtner, Jürgen Kroll — University of Potsdam, Institute of Nutritional Science, A.-Scheunert Allee 114-116, D-14558 Nuthetal, Germany |
| Venue | J. Agric. Food Chem. **2005**, **53** (10), 4228-4235. Received 24 November 2004; revised 23 March 2005; accepted 24 March 2005; web 14 April 2005 |
| DOI / article ID | 10.1021/jf0480290 (printed as `JF0480290`) |
| Ligands | **chlorogenic acid** (MW 336 Da), **ferulic acid** (194 Da), **gallic acid** (188 Da), **quercetin**, **rutin** (quercetin-3-O-rhamnoglucoside), **isoquercetin** (quercetin-3-O-glucoside) |
| Proteins | **bovine serum albumin (BSA)** — the main model; **human serum albumin (HSA)**; **soy glycinin (SG)**; **lysozyme**; **trypsin**; and, for ferulic acid only, **gelatin, milk whey proteins, α-amylase and whole human serum** |
| Quantity measured | **K_D, a DISSOCIATION constant in M** (so a *small* number means *strong* binding — the opposite convention from the rest of the binding batch, which reports association constants K_b or K in M^-1), plus **n**, the number of sites per protein molecule, and **ΔG** in kJ/mol |
| Methods | **(i) Hummel-Dreyer size-exclusion chromatography (direct)**; **(ii) enhancement of quercetin's own fluorescence on binding (indirect)**; **(iii) far- and near-UV circular dichroism** for structure |
| Explicit scope | **noncovalent binding only.** The authors' earlier work was on covalent phenol-protein adducts (ref. 2), and they state that "covalent binding is reported to take place at conditions above pH 7" — the whole design here is below or at pH 7, and reversibility is demonstrated with 8 M urea (Fig. 4) |
| Soy glycinin preparation | from **defatted unheated soy flour** (~52 % protein, Sigma type 1) by the method of Thanh & Shibasaki 1976 — **the same isolation paper Damodaran & Kinsella 1981 cite for their absorptivity**. Kjeldahl protein **99.8 %**; SDS-PAGE purity **~95 %** by densitometry |

## 1. Why it matters

**(a) It is a third ligand class on the same protein sites, and it is far stronger than the other
two.** The matrix layer prices two kinds of small-molecule capture on protein: **reversible
hydrophobic binding** of ketones and aldehydes (`REVERSIBLE_BINDING` in
`src/kinetic_core/parameters_matrix.py`, from Damodaran on soy and Andriot on beta-lactoglobulin),
and **covalent adduction** of aldehydes and HMF to amine and thiol pools
(`BINDING_CLASSES` in `src/kinetic_core/matrix_sites.py`). This paper measures a third: **phenolics,
noncovalently, on soy glycinin and the albumins.** Inverting its dissociation constants to the
batch's usual convention:

| ligand / protein | K_D printed | K_A = 1/K_D (mine) | for comparison |
|---|---:|---:|---|
| quercetin / **soy glycinin** | 4.1 µM | **2.4 x 10^5 M^-1** | 2-nonanone / native soy = **930 M^-1** (Damodaran 1981) |
| quercetin / HSA | 1.9 µM | **5.3 x 10^5 M^-1** | 2-nonanone / β-lactoglobulin = **2 440 M^-1** (Andriot 2000) |
| quercetin / BSA | 9.1 µM | **1.1 x 10^5 M^-1** | |
| chlorogenic acid / BSA, pH 4.8 | 48 µM | **2.1 x 10^4 M^-1** | |
| ferulic acid / BSA, pH 4.8 | 58 µM | **1.7 x 10^4 M^-1** | |
| quercetin / lysozyme | 2 000 µM | **5 x 10^2 M^-1** | the one weak case, and the authors call it non-specific |

**Quercetin binds soy glycinin about 260 times harder than 2-nonanone binds native soy protein
(mine).** That is the finding a Maillard model should take from this paper: in a plant-protein pot
that also contains phenolics — which soy, pea, coffee and almost every real plant matrix does — the
protein's hydrophobic capacity is occupied first and hardest by the phenolics, not by the aroma
volatiles the odour-activity layer is tracking. **The repository has no term for this at all.**

**(b) Chlorogenic acid is a keyed compound in this repository and this is a binding constant for
it.** `data/keys/compounds.yml` carries `chlorogenic_acid` as an id (line 607). This paper gives
**four** numbers for it: K_D = 48 ± 8 µM with n = 1.0 ± 0.1 at pH 4.8, and K_D = 210 ± 30 µM with
n = 2.0 ± 0.1 at pH 7.0, both against BSA. There is no soy glycinin value for chlorogenic acid here
(Fig. 2 covers ferulic acid only). Still, this is the only protein-binding constant for a
registry-keyed phenolic anywhere in the corpus, and it is directly relevant to the coffee lane,
where chlorogenic acid is a major precursor.

**(c) It contradicts Damodaran on the direction of the heat effect, and the conflict is worth
carrying.** Damodaran & Kinsella 1981 measured that **heating soy protein to 90 C for 1 h INCREASED**
2-nonanone binding by 33 % (K 930 -> 1 240 M^-1, sites unchanged). Rawel measures that heating BSA to
90 C for 10 min **DECREASED** quercetin binding by 39 % (K_D 9.1 -> 15 µM, i.e. K_A 1.10e5 ->
6.7e4 M^-1, mine). Different protein, different ligand, different heat load — so the two are not
strictly in conflict — but they are the corpus's only two direct measurements of what a cook does to
a protein's binding capacity, **and they point in opposite directions**. Rawel's stated mechanism is
polymerisation of the heated protein and loss of exposed surface; Damodaran's is reorganisation of
subunits enhancing the existing hydrophobic sites. **`matrix_sites.py`'s note that "the sites are
charged once, at the start of the cook" and that heating's effect on site densities "is not
modelled" is exactly the gap these two papers straddle**, and neither of them resolves it.

**(d) It is a soy paper whose soy protein is a defined fraction.** `data/species/protein_matrices.yml`
builds `soy_isolate` from commercial isolates. Rawel's SG is **glycinin (11S) at 99.8 % Kjeldahl
protein and ~95 % purity, from unheated flour** — the cleanest soy preparation in the batch, and the
same Thanh & Shibasaki fractionation Damodaran used. Its molar mass is recoverable from the Methods
(Flags 2), which few of these papers allow.

**(e) What it does NOT do.** It contains **no odour threshold**, no volatile, no aroma compound, no
kinetics, no temperature-dependence of any binding constant in the Arrhenius sense (its "temperature"
series is a *pre-treatment*, not a measurement at temperature), and nothing about the covalent
chemistry the odour layer's aldehydes undergo. **It cannot add a row to `MATRIX_THRESHOLDS` and it
cannot change a `BINDING_CLASSES` bracket.**

## 2. Methods as they matter to a model

- **Two independent methods, and they measure different things.**
  - **Hummel-Dreyer (HD)**: isocratic size-exclusion chromatography on an Econo-Pac P-6 cartridge
    (exclusion > 6 kDa, Bio-Rad), diode-array detection, **column at 25 C**. Eluent A is buffer
    alone; eluent B is the phenolic at **0.5 mM** in the same buffer; a gradient mixer sets the
    working phenol concentration at **0.03-0.5 mM**. Flow **0.8-1 mL/min**, injection **30-50 µL**.
    Bound phenol quantified by **internal calibration** — a series of injections at fixed protein
    (**45-181 µM**) and rising phenol. The free concentration is taken as **the phenol concentration
    in the eluent buffer**, which is the method's defining approximation.
  - **Quercetin fluorescence (indirect)**: quercetin's own emission is **enhanced** when it binds,
    by resonance energy transfer. Excitation **370 nm** (slit 18 nm), emission scanned **380-900 nm**,
    peak height at the maximum, **526 ± 2 nm**. Method adapted from Guharay 2001 (ref. 13). Crucially,
    **the protein alone also fluoresces at 370 nm excitation, and its emission was subtracted** —
    the authors say this correction is theirs and not in the source method. Fitted with eq 10,
    ΔF = ΔF_max·[protein_total] / (K_D + [protein_total]), by non-linear least squares in Microcal
    Origin 6.0.
- **Buffers and pH, which is the paper's main axis.** **0.05 M Hepes** at pH 7, 6, 4.8 and 3.5, or
  **0.05 M Na-acetate** at the same pH values; also **0.067 M Na-phosphate at pH 7.4** for the
  "physiological" ferulic-acid comparison, and **0.1 M sodium phosphate** for all CD work (Hepes
  interferes with CD). **Gallic acid had to be moved out of Hepes entirely** because it gave
  irreproducible double peaks — "most likely the result of binding of gallic acid to the Hepes".
- **Concentrations, exactly as stated.** Protein bulk solutions: **BSA and HSA 10 g/L = 150 µM**;
  trypsin "**3534 g/L (150 µM)**" (an obvious typographical error, Flags 4); **SG 10 g/L = 29.1 µM**;
  **lysozyme 40 g/L = 2.74 mM**. Titration ranges: BSA and HSA **0-90 µM**, trypsin **0-90 µM**,
  lysozyme **0-1.64 mM**, SG **0-17.46 µM**, with "**10-20 concentrations for each protein
  (n = 4 for each concentration)**". Phenol bulk solutions in **DMSO**: quercetin 1 g/L = 3.31 mM,
  rutin 1 g/L = 1.50 mM, isoquercetin 1 g/L = 2.15 mM; working concentrations **15 µM quercetin,
  15 µM rutin, 1 / 5 / 15 µM isoquercetin**, with **final DMSO below 1 %**.
- **The three perturbation series, all on BSA, all with quercetin, all by fluorescence.**
  - **pH**: 0.05 M Hepes at pH **7 (standard), 6 and 5**.
  - **Ionic strength**: 0.05 M Hepes pH 7 with **50, 250 and 500 mM NaCl**.
  - **Temperature**: **"BSA solutions were incubated at 40, 60 and 90 C for 10 min and cooled to room
    temperature before addition of quercetin and measurement."** **This is a pre-treatment, not a
    binding temperature.** Every K_D in Table 2 is measured at room temperature (25 C); the
    temperature rows report how a *previously heated* protein binds when it is cold again. This
    distinction is easy to lose and it matters (Flags 5).
- **"Standard conditions" defined**: 0.05 M Hepes buffer, pH 7, **room temperature (25 C)**, 0 mM
  NaCl (Table 2 footnote a).
- **Reversibility test.** Repeated under standard conditions **with 8 M urea**, at 30 µM BSA and HSA,
  6.98 µM SG, 0.66 mM lysozyme, quercetin held at 15 µM. **The fluorescence enhancement is abolished
  by urea (Fig. 4), which is the paper's proof that the binding is noncovalent.** Confirmed for
  lysozyme, SG and HSA in the Supporting Information.
- **Circular dichroism.** Far-UV **178-260 nm** at 3 µM protein (0.2 g/L) in 0.1 M sodium phosphate,
  1 mm quartz cell, 1 nm steps, 50 nm/min, 4 s response; near-UV **250-320 nm** at 3 µM (0.2 g/L BSA,
  **1.029 g/L SG**), 5 mm path. Mean residue MW from the SWISS-PROT sequence. Deconvolution by
  **CDPro** with **CONTIN/LL, SELCON and CDSSTR** against a 48-protein reference set. Phenol
  concentrations for the structural work: **15, 30, 45, 60, 75 µM**.
- **Statistics.** "The analysis was repeated **at least four times**"; **the averaged data** were then
  fitted. Non-linear least squares in Microcal Origin 6.0. **The authors explicitly warn that the
  quoted standard deviations are fit uncertainties, not experimental ones**: "The K_D standard
  deviation generally underestimates the real uncertainties of the experiment, and every nonlinear
  least squares regression therefore actually represents the extent of the goodness of the curve
  fitting rather than the real error (16)." **That sentence should travel with every number in this
  paper.**
- **What could not be measured, and why.** Quercetin, rutin and isoquercetin **could not be run by
  HD at all** — low solubility, and worse, they bind the column material; adding 1 % DMSO or 20 %
  ethanol did not rescue it. **Trypsin gave no fluorescence enhancement with quercetin**, so no K_D.
  **Rutin and isoquercetin gave no fluorescence enhancement with any protein**, so the only glucoside
  results are the CD spectra. **α-Amylase was insoluble at pH 7.4** at the required concentration.

## 3. Tables re-typed

### Table 1. "Dissociation Constants (K_D) and Number of Sites (n) for Binding of Selected Phenolic Compounds as Determined by the HD (Model = One Site Binding)"

Header exactly as printed: `K_D x 10^-6 (M)`. **Read as K_D = <value> x 10^-6 M, i.e. the values are
in micromolar** — confirmed from Figure 1's own x-axis, which runs 0 to 1.2 **mM** with the
half-saturation near 0.2 mM = 200 µM against the table's 210 (Flags 1). All rows are BSA.

| protein | phenol | buffer | K_D x 10^-6 (M) | n |
|---|---|---|---:|---:|
| BSA | ferulic acid | 0.05 M Hepes, pH 4.8 | 58 ± 5 | 2.2 ± 0.1 |
| BSA | ferulic acid | 0.05 M Hepes, pH 7.0 | 249 ± 14 | 3.6 ± 0.1 |
| BSA | chlorogenic acid | 0.05 M Hepes, pH 4.8 | 48 ± 8 | 1.0 ± 0.1 |
| BSA | chlorogenic acid | 0.05 M Hepes, pH 7.0 | 210 ± 30 | 2.0 ± 0.1 |
| BSA | gallic acid | 0.05 M Na-Ac, pH 3.5 | 2300 ± 500 | 14 ± 3 |
| BSA | gallic acid | 0.05 M Na-Ac, pH 4.8 | 1000 ± 100 | 5.3 ± 0.4 |
| BSA | gallic acid | 0.05 M Na-Ac, pH 6.0 | 380 ± 80 | 3.7 ± 0.4 |

### Table 2. "Dissociation Constants (K_D) for Binding of Quercetin to Different Proteins Depending on the Reaction Conditions as Determined by the Enhancement of Its Fluorescence Intensity (Model = One Site Binding)"

Footnote a: "Means standard conditions: **0.05 M Hepes buffer, pH 7, room temperature (25 °C), and
0 mM NaCl**." Same `x 10^-6 (M)` convention. **No n column** — the fluorescence method as applied
here returns K_D only.

| protein | conditions | K_D x 10^-6 (M) |
|---|---|---:|
| BSA | standard (a) | 9.1 ± 0.5 |
| BSA | pH 6 | 14.9 ± 0.4 |
| BSA | pH 5 | 24 ± 2 |
| BSA | T = 40 °C | 10.2 ± 0.7 |
| BSA | T = 60 °C | 13 ± 1 |
| BSA | T = 90 °C | 15 ± 2 |
| BSA | 50 mM NaCl | 17.3 ± 0.4 |
| BSA | 250 mM NaCl | 22 ± 1 |
| BSA | 500 mM NaCl | 24.0 ± 0.9 |
| HSA | standard | 1.9 ± 0.2 |
| lysozym [sic] | standard | (2000 ± 800) |
| **sojaglycinin** [sic] | standard | **4.1 ± 0.6** |
| trypsin | standard | *(cell empty — no enhancement was observed)* |

The lysozyme value is printed **in parentheses** in the original, and the text explains why: it could
only be obtained "after setting ΔF_max = 0.05", and "the K_D value is rather high, meaning
nonspecific binding". **The parentheses are the authors' own flag and must be preserved.**
"T = 40/60/90 °C" means the BSA was **heated for 10 min and cooled before measurement** (Flags 5).

### Table 3. "Effect of Reaction Conditions on the Allocation of the Secondary Structure Elements According to Refs 15 and 16 by Applying 48 Protein Reference Set (190-240 nm) and CONTIN Method ([BSA] = 3 µM)"

All values are percentages of the secondary structure of **BSA**.

| condition | α-helix | β-strand | β-turn | unordered |
|---|---:|---:|---:|---:|
| **pH value** | | | | |
| pH 5 | 63.2 | 3.2 | 13.4 | 20.1 |
| pH 6 | 67.1 | 3.1 | 12.6 | 17.3 |
| pH 7 | 66.9 | 2.5 | 12.0 | 18.6 |
| **ionic strength** | | | | |
| 0 mM NaCl | 66.9 | 2.5 | 12.0 | 18.6 |
| 50 mM NaCl | 61.5 | 3.2 | 13.8 | 21.6 |
| 250 mM NaCl | 61.6 | 2.9 | 14.8 | 20.8 |
| 500 mM NaCl | 63.1 | 2.8 | 12.0 | 22.1 |
| **temperature** | | | | |
| 25 °C | 66.9 | 2.5 | 12.0 | 18.6 |
| 40 °C | 64.1 | 3.3 | 13.3 | 19.3 |
| 60 °C | 60.7 | 3.2 | 14.4 | 21.7 |
| **90 °C** | **50.6** | **6.5** | **16.3** | **26.4** |

(The pH 7 / 0 mM NaCl / 25 °C row is the same control repeated three times, as it must be.)

### Table 4. "Standard Gibbs Free Energy (ΔG) for Binding of Selected Phenolic Compounds as Determined Using Different Methods"

Two method columns; a blank cell means that method was not applied to that pair.

| protein | phenol | pH | ΔG (kJ/mol), quercetin fluorescence | ΔG (kJ/mol), HD method |
|---|---|---:|---:|---:|
| BSA | quercetin | 7.0 | -30.1 | |
| lysozym | quercetin | 7.0 | -16.8 | |
| **sojaglycinin** | **quercetin** | **7.0** | **-32.2** | |
| HSA | quercetin | 7.0 | -34.0 | |
| BSA | ferulic acid | 4.8 | | -25.6 |
| BSA | ferulic acid | 7.0 | | -22.0 |

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| ferulic-acid **binding capacity** across eight proteins | "a BC of **0-22 mg/g protein** was determined" | Results (describing **Fig. 2**, which is the only place the per-protein values exist) |
| HSA vs whole human serum, ferulic acid, same 7.7 g/L HSA | **10.3 ± 0.1** and **8.9 ± 0.1** mg ferulic acid / g protein | Results |
| ionic-strength asymptote | K_D "reaching an asymptotic value of **25 µM**" | Results |
| K_D at the isoelectric point (pH 5) | **24 µM** — "a similar K_D value" to the high-salt asymptote | Results |
| BSA isoelectric point | pH **5** (given as **5.1** three paragraphs later — Flags 6) | Results |
| BSA domain charges at pH 7 | domain I **-11**, domain II **-7**, domain III **+1** | Results (from Peters 1996, ref. 23) |
| lysozyme isoelectric point | **10.8** | Results |
| gallic acid pKa | pKa1 = **6.89**, pKa2 = **10.16** | Results (Tam & Takács-Novák 2001) |
| quercetin pKa | pKa1 = **7.3**, pKa2 = **8.4** | Results (Kopacz 2003) |
| molecular masses quoted | ferulic acid **194 Da**, chlorogenic acid **336 Da**, gallic acid **188 Da**; phenols generally 180-700 Da; proteins 14 000-350 000 Da | Results |
| BSA N->F transition | below pH **4.3** (4.3-3.5), unfolding of domain III | Results (Friedli 1996) |
| BSA thermal transitions | conformational change at **58.1 °C** (DSC, ref. 30); denaturation above **62 °C** (ref. 31) | Results, quoted |
| BSA control secondary structure by CONTIN | **65 % α-helix, 3 % β-strand, 12 % β-turn, 20 % unordered** | Results (the Table 3 pH 7 row prints 66.9 / 2.5 / 12.0 / 18.6 — Flags 6) |
| far-UV BSA spectrum | minima at **209 and 222 nm**, maximum at **190-195 nm** | Results |
| quercetin solubility | **< 1 g/L** in aqueous solution | Results |
| quercetin loss to membranes | "**60-80 %** of quercetin may bind to separation membranes applied in equilibrium dialysis or ultrafiltration" | Results, quoting Boulton 1998 — **the paper's reason for avoiding dialysis** |
| soy glycinin purity | Kjeldahl protein **99.8 %**; SDS-PAGE purity **~95 %** | Methods |

**Ferulic acid across eight proteins, the isotherms, the urea experiment and the CD spectra:
FIGURE-ONLY.** Figures 1 (chlorogenic acid / BSA isotherm), 2 (**ferulic-acid binding capacity for
gelatin, milk whey, soy glycinin, α-amylase, lysozyme, BSA, HSA and human serum, at pH 4.8 and
pH 7.4** — eluent 0.24 mM ferulic acid, protein 12 g/L, BC in mg/g), 3 (quercetin-BSA fluorescence
fit), 4 (urea abolishing the enhancement) and 5 (near-UV CD of soy glycinin with quercetin, rutin
and isoquercetin at 90 µM) are images. **Figure 2 is the paper's most repository-relevant
measurement — a phenolic binding capacity for soy glycinin, gelatin and whey side by side — and not
one of its values is printed.** Per house rule none is typed as a number.

### Fit quality printed in figure captions (the only goodness-of-fit numbers in the paper)

- Fig. 1 (chlorogenic acid / BSA, pH 7, HD): model = one site binding, **χ² = 0.01, R² = 0.96**,
  **n = 2.0 ± 0.1**, and the K_D quoted as **"K_D x 10^-6 = 0.21 ± 0.03 M"** — which conflicts with
  Table 1's 210 ± 30 by a factor of 1 000 and is resolved in Flags 1.
- Fig. 3 (quercetin / BSA, fluorescence): **χ² = 0.00001, R² = 0.99**, **ΔF_max = 0.237 ± 0.004**,
  and **K_D x 10^-6 = 9.1 ± 0.5** — which **agrees exactly** with Table 2's first row.

### Arithmetic on the printed constants (all mine)

**1. Inverting to association constants, so the batch can be compared on one scale.**
K_A = 1/K_D: quercetin **HSA 5.26e5**, **soy glycinin 2.44e5**, **BSA 1.10e5**, **lysozyme
5.0e2 M^-1**; chlorogenic acid on BSA **2.08e4** (pH 4.8) and **4.76e3** (pH 7.0); ferulic acid on
BSA **1.72e4** (pH 4.8) and **4.02e3** (pH 7.0); gallic acid on BSA **4.35e2** (pH 3.5), **1.00e3**
(pH 4.8), **2.63e3 M^-1** (pH 6.0).

**2. Phenolic against carbonyl, on soy.** Quercetin / soy glycinin K_A = **2.44e5 M^-1**;
2-nonanone / native whole soy protein (Damodaran 1981, dialysis, pH 8.0, 25 C) = **930 M^-1**.
**Ratio 262x.** Even against the strongest carbonyl in the batch (nonanal, 1 094 M^-1) the ratio is
**223x**. Different proteins (glycinin fraction vs whole soy), different methods, different pH — so
this is a cross-study, cross-method comparison and nothing more — but the size of the gap is not
something method differences plausibly explain. **The plant phenolics in a soy pot outcompete the
aroma carbonyls for the protein's hydrophobic capacity by two to three orders of magnitude in
affinity.**

**3. A per-gram constant for soy glycinin, conditional on a site count the paper does not give.**
The Methods let the molar mass be recovered: **SG at 10 g/L = 29.1 µM implies 343 600 g/mol (mine)**,
consistent with the glycinin hexamer. With **n = 1** (the model fitted, but n is not reported for the
fluorescence rows), k_g = n·K_A/MW = 2.44e5 / 343 600 = **0.71 L/g (mine)**. Against the
`REVERSIBLE_BINDING` soy rows (4.40e-3 to 4.38e-2 L/g) that is **16 to 160 times larger**. If n were
larger than 1, larger still. **This is a bound, not a value**, and it is marked as such because the
site count is missing.

**4. The pH effect, which is the paper's clearest quantitative result, and it has two parts.**
For the **HD** phenolic acids, moving pH 4.8 -> 7.0 **weakens** binding: ferulic K_D 58 -> 249 µM
(**4.3x weaker**), chlorogenic 48 -> 210 µM (**4.4x weaker**) — while **n rises**, 2.2 -> 3.6 and
1.0 -> 2.0. Gallic acid runs the **other way**: K_D 2 300 -> 1 000 -> 380 µM as pH goes 3.5 -> 4.8
-> 6.0 (**6.1x stronger**), with n falling 14 -> 5.3 -> 3.7. For **quercetin by fluorescence**,
lowering pH 7 -> 5 weakens binding: K_D 9.1 -> 24 µM (**2.6x weaker**). **So the pH dependence is
ligand-specific in both magnitude and sign**, which is the single most transferable lesson here for
a model that wants one pH factor.

**5. The salt effect.** K_D 9.1 -> 17.3 -> 22 -> 24.0 µM for 0 -> 50 -> 250 -> 500 mM NaCl:
**2.6x weaker binding across the range**, with most of it (1.9x) spent on the **first 50 mM**. The
authors note it plateaus at ~25 µM. **A pot with 0.5 % salt sits at the plateau**, so the salt
correction is a step, not a slope.

**6. The heat effect, and its comparison with Damodaran.** K_D 9.1 (unheated) -> 10.2 (40 C) ->
13 (60 C) -> 15 µM (90 C): **1.65x weaker binding after the strongest heat treatment (mine)**, and
the CD in Table 3 shows why — α-helix falls **66.9 % -> 50.6 %** and unordered structure rises
**18.6 % -> 26.4 %** over the same series. **Damodaran's soy/2-nonanone result under a comparable
treatment (90 C, but 1 h rather than 10 min) is 1.33x STRONGER.** Two papers, opposite signs
(§1(c)).

**7. Consistency of the free energies with the constants, which does NOT close.** ΔG = -RT ln K_A
at the stated 25 C (RT = 2.479 kJ/mol) gives: BSA/quercetin **-28.8** against a printed **-30.1**;
soy glycinin **-30.8** against **-32.2**; HSA **-32.7** against **-34.0**; lysozyme **-15.4** against
**-16.8**; ferulic acid pH 4.8 **-24.2** against **-25.6**; ferulic acid pH 7.0 **-20.6** against
**-22.0**. **Every printed value is 1.3 to 1.5 kJ/mol more negative than the recomputation, in the
same direction, on all six rows.** The implied temperature is **310-319 K**, not 298 K, and it is not
even constant across the six rows. The expression the authors used is in Supporting Information III
(eq 9), which is not on disk. **Prefer the K_D column; treat Table 4 as unreproducible until the SI
is retrieved** (Flags 3).

**8. Ferulic acid as a "metabolite sponge" check.** Whole human serum at 7.7 g/L HSA-equivalent binds
**8.9 mg/g**; pure HSA at the same 7.7 g/L binds **10.3 mg/g**. **HSA alone accounts for 86 % of the
serum's ferulic-acid binding capacity (mine)**, which is the paper's claim, quantified.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** **`chlorogenic_acid` is keyed** (id
`chlorogenic_acid`) — the only one of the six ligands that is. **Ferulic acid, gallic acid,
quercetin, rutin and isoquercetin are all absent.** (`4_vinylguaiacol`, ferulic acid's
decarboxylation product, is keyed, but that is a different compound.) On the protein side, soy is a
matrix in `data/species/protein_matrices.yml` but as a commercial **isolate**, not as the **glycinin
fraction** this paper purifies. Every row below shares: **25 C measurement temperature; 0.05 M Hepes
or Na-acetate buffer as stated; non-linear least squares fit of a ONE-SITE binding model in Microcal
Origin 6.0 to data averaged over at least four repetitions; and the authors' own warning that the ±
values are curve-fit uncertainties and "underestimate the real uncertainties of the experiment".**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| ferulic acid / BSA, dissociation constant | 58 ± 5 | µM | 0.05 M Hepes pH 4.8, 25 C, HD | Table 1 p. 4230 | **binding_constant** |
| ferulic acid / BSA, sites | 2.2 ± 0.1 | per protein molecule | as above | Table 1 | **binding_constant** |
| ferulic acid / BSA | 249 ± 14 | µM | 0.05 M Hepes pH 7.0, 25 C, HD | Table 1 | **binding_constant** |
| ferulic acid / BSA, sites | 3.6 ± 0.1 | — | as above | Table 1 | **binding_constant** |
| **chlorogenic acid / BSA** | **48 ± 8** | µM | 0.05 M Hepes **pH 4.8**, 25 C, HD | Table 1 | **binding_constant** — the only registry-keyed ligand in the paper |
| chlorogenic acid / BSA, sites | 1.0 ± 0.1 | — | as above | Table 1 | **binding_constant** |
| **chlorogenic acid / BSA** | **210 ± 30** | µM | 0.05 M Hepes **pH 7.0**, 25 C, HD | Table 1; isotherm Fig. 1 (χ² 0.01, R² 0.96) | **binding_constant** |
| chlorogenic acid / BSA, sites | 2.0 ± 0.1 | — | as above | Table 1, Fig. 1 caption | **binding_constant** |
| gallic acid / BSA | 2300 ± 500 | µM | 0.05 M Na-acetate pH 3.5, 25 C, HD | Table 1 | **binding_constant** — the authors call this **unspecific binding** and the fit unsaturated |
| gallic acid / BSA, sites | 14 ± 3 | — | as above | Table 1 | **binding_constant** — "indicative of an unspecific binding" |
| gallic acid / BSA | 1000 ± 100 | µM | 0.05 M Na-acetate pH 4.8, 25 C, HD | Table 1 | **binding_constant** |
| gallic acid / BSA, sites | 5.3 ± 0.4 | — | as above | Table 1 | **binding_constant** |
| gallic acid / BSA | 380 ± 80 | µM | 0.05 M Na-acetate pH 6.0, 25 C, HD | Table 1 | **binding_constant** |
| gallic acid / BSA, sites | 3.7 ± 0.4 | — | as above | Table 1 | **binding_constant** |
| quercetin / BSA | 9.1 ± 0.5 | µM | 0.05 M Hepes pH 7, 25 C, 0 mM NaCl (standard), fluorescence | Table 2 p. 4232; fit Fig. 3 (χ² 1e-5, R² 0.99) | **binding_constant** |
| quercetin / BSA | 14.9 ± 0.4 | µM | as standard but **pH 6** | Table 2 | **binding_constant** |
| quercetin / BSA | 24 ± 2 | µM | as standard but **pH 5** (= BSA's isoelectric point) | Table 2 | **binding_constant** |
| quercetin / BSA | 10.2 ± 0.7 | µM | BSA **pre-heated 40 C / 10 min**, then cooled; measured at 25 C | Table 2 | **binding_constant** (a pre-treatment, Flags 5) |
| quercetin / BSA | 13 ± 1 | µM | BSA **pre-heated 60 C / 10 min** | Table 2 | **binding_constant** |
| quercetin / BSA | 15 ± 2 | µM | BSA **pre-heated 90 C / 10 min** | Table 2 | **binding_constant** |
| quercetin / BSA | 17.3 ± 0.4 | µM | standard + **50 mM NaCl** | Table 2 | **binding_constant** |
| quercetin / BSA | 22 ± 1 | µM | standard + **250 mM NaCl** | Table 2 | **binding_constant** |
| quercetin / BSA | 24.0 ± 0.9 | µM | standard + **500 mM NaCl** | Table 2 | **binding_constant** |
| quercetin / HSA | 1.9 ± 0.2 | µM | standard | Table 2 | **binding_constant** — the strongest in the paper |
| **quercetin / soy glycinin** | **4.1 ± 0.6** | µM | standard (0.05 M Hepes pH 7, 25 C) | Table 2 | **binding_constant** — the repository-relevant row |
| quercetin / lysozyme | (2000 ± 800) | µM | standard; **only obtainable after fixing ΔF_max = 0.05** | Table 2 | **binding_constant** — **printed in parentheses by the authors**; non-specific by their own reading |
| quercetin / trypsin | **no value** | — | standard; **no fluorescence enhancement observed** | Table 2 | level_only (a measured null) |
| rutin and isoquercetin / any protein | **no value** | — | standard; **no fluorescence enhancement with any protein** | Results | level_only (a measured null) |
| association constants (inverted) | 1.72e4 / 4.02e3 (ferulic, pH 4.8 / 7.0); 2.08e4 / 4.76e3 (chlorogenic); 4.35e2 / 1.00e3 / 2.63e3 (gallic, pH 3.5 / 4.8 / 6.0); 1.10e5 (quercetin/BSA); 5.26e5 (HSA); **2.44e5 (soy glycinin)**; 5.0e2 (lysozyme) | M^-1 | as each row above | derived (mine), K_A = 1/K_D | binding_constant (arithmetic on printed values) |
| quercetin/soy glycinin vs 2-nonanone/native soy | **262** | x stronger in K_A | cross-study: Rawel pH 7 / 25 C fluorescence vs Damodaran pH 8 / 25 C dialysis | derived (mine) | within_study_ratio — **cross-study, cross-method, cross-protein-fraction** |
| soy glycinin molar mass, recovered | **343 600** | g/mol | from "10 g/L (29.1 µM)" | derived (mine) from Methods | derived_assumption (arithmetic on a printed pair) |
| per-gram binding constant, quercetin / soy glycinin, **if n = 1** | **0.71** | L/g | pH 7, 25 C | derived (mine) | derived_assumption — **n is NOT reported for any fluorescence row**; a lower bound if n > 1 |
| pH effect, ferulic and chlorogenic on BSA, 4.8 -> 7.0 | 4.3 and 4.4 | x weaker | HD | derived (mine) | within_study_ratio |
| pH effect, gallic on BSA, 3.5 -> 6.0 | 6.1 | x **stronger** | HD, Na-acetate | derived (mine) | within_study_ratio — **opposite sign to the other two acids** |
| pH effect, quercetin on BSA, 7 -> 5 | 2.6 | x weaker | fluorescence | derived (mine) | within_study_ratio |
| ionic-strength effect, 0 -> 500 mM NaCl | 2.6 | x weaker, plateauing at ~25 µM | quercetin / BSA, pH 7 | derived (mine); the asymptote is printed | within_study_ratio |
| thermal-pretreatment effect, 25 -> 90 C | 1.65 | x weaker | quercetin / BSA; BSA heated 10 min then cooled | derived (mine) | within_study_ratio — **contradicts Damodaran's +1.33x on soy/2-nonanone** |
| BSA secondary structure vs condition | see Table 3 | % α-helix / β-strand / β-turn / unordered | CD, CONTIN, 48-protein reference set, 3 µM BSA | Table 3 p. 4234 | level_only |
| BSA α-helix loss on heating to 90 C | 66.9 -> 50.6 | % | 10 min, 0.1 M Na-phosphate | Table 3 | level_only |
| standard Gibbs free energies | -30.1 / -16.8 / **-32.2** / -34.0 (quercetin on BSA / lysozyme / soy glycinin / HSA) and -25.6 / -22.0 (ferulic on BSA, pH 4.8 / 7.0) | kJ/mol | pH as stated | Table 4 p. 4234 | **REFUSE until the SI is retrieved** — none of the six reproduces from the printed K_D at 25 C (Flags 3) |
| ferulic acid binding capacity, eight proteins, two pH values | 0 to 22 (the range only) | mg/g protein | 0.24 mM ferulic acid, 12 g/L protein, pH 4.8 Hepes or pH 7.4 Na-phosphate | Results; **the per-protein values are Fig. 2** | **figure_only** except the printed range |
| ferulic acid on HSA vs whole human serum | 10.3 ± 0.1 and 8.9 ± 0.1 | mg/g protein | 7.7 g/L HSA in both | Results | level_only (printed) |
| binding isotherms; ferulic acid per protein; fluorescence titrations; urea reversibility; soy glycinin near-UV CD | — | — | — | Figs. 1-5 | **figure_only** |

### What can and cannot be put next to these

**(a) Nothing here goes into `REVERSIBLE_BINDING` as it stands, and the reason is the missing site
count.** Every row of that table is a **per-gram** constant, k_g = n·K/MW. This paper reports n for
the seven HD rows (all BSA — a protein the repository does not carry) and **reports no n at all for
the fluorescence rows, which are the ones on soy glycinin, HSA and lysozyme**. The soy glycinin
molar mass is recoverable (343 600 g/mol) but the site count is not, so the per-gram constant is a
bound, not a value.

**(b) The pH lesson is the most portable thing in the paper and it is a negative one.** Three
phenolic acids on one protein by one method give pH dependences of 4.3x, 4.4x and 6.1x — **and the
gallic acid one runs in the opposite direction**. `parameters_matrix.py` carries a
`PH_ADDUCT_GATE_BELOW` and a `PH_ADDUCT_GATE_UNCERTAIN_BELOW` as single thresholds. **This paper is
evidence that a single pH gate cannot be right across ligand classes**, at least for noncovalent
binding, because the sign itself is ligand-dependent (it tracks the ligand's own pKa and the
protein's charge state, both of which the paper works through).

**(c) The heat result is a genuine conflict and should be recorded as one, not averaged.** See
§1(c) and §3.6.

**(d) What cannot be transported.** Every constant is a **noncovalent** one, measured deliberately
below the pH at which the authors' own earlier work says covalent phenol-protein bonds form ("above
pH 7"). A Maillard pot at 100-180 C is exactly the regime where that covalent chemistry runs, and
this paper says nothing about it. Nothing here is a rate. Nothing here is at any temperature but 25 C
(the "temperature" rows are pre-treatments). Nothing here involves an aroma volatile.

**(e) A structural warning worth keeping.** Table 3 shows BSA's **secondary** structure surviving
everything except 90 C, while the near-UV CD (Fig. 5, soy glycinin) shows the **tertiary** structure
changing on phenol binding at every condition. The paper's conclusion — "the noncovalent binding of
phenolic compounds ... has no effect on the secondary structure of the proteins studied but causes
significant changes in the tertiary structure" — is the mechanism by which a phenolic-rich matrix
could change a protein's site accessibility without denaturing it. `protein_matrices.yml`'s
`amine_available_band` (0.4-1.0, "declared, not measured") is the parameter that would carry such an
effect, and this paper is a qualitative reason to think the band is real.

## 5. Flags

1. **The `K_D x 10^-6 (M)` header is strictly backwards, and one figure caption conflicts with its
   own table by a factor of 1 000 — but the reading is settled from the figure's axis.** Taken
   literally, "K_D x 10^-6 = 58" means K_D = 5.8 x 10^7 M, which is absurd; the intended reading is
   **K_D = 58 x 10^-6 M = 58 µM**, and this convention is used in both Table 1 and Table 2. The
   conflict: **Figure 1's caption prints "K_D x 10^-6 = 0.21 ± 0.03 M" for chlorogenic acid on BSA at
   pH 7, while Table 1 prints 210 ± 30 for the same pair.** Resolved by reading the figure itself at
   200 dpi: **its x-axis is `[CA_free] mM`, runs 0 to 1.2, and the curve is at half of its plateau
   (B ≈ 1.0 of n = 2.0) near 0.2 mM = 200 µM** — which matches Table 1's 210 µM. **So the figure
   caption's "0.21" is in mM and its "x 10^-6" is a misprint; Table 1 is correct.** A second,
   independent confirmation: Figure 3's caption prints "K_D x 10^-6 = 9.1 ± 0.5" and Table 2's first
   row prints 9.1 ± 0.5 — the two agree exactly, so the table convention is the µM one. **Use the
   tables; do not use the Figure 1 caption's number.**
2. **The site count is missing for exactly the rows the repository would want.** Table 1 gives n for
   seven BSA rows. **Table 2 gives no n at all** — and Table 2 is where soy glycinin, HSA and
   lysozyme live. The fluorescence model (eq 10) fits ΔF_max and K_D, not n, so the site count is not
   recoverable from the printed values. **Without n, no per-gram constant can be computed for soy
   glycinin**, which is the single most useful thing this paper could have given the matrix layer.
3. **Table 4's free energies do not reproduce from Table 1 and Table 2, and the derivation is in
   Supporting Information that is not on disk.** All six printed ΔG values are **1.3 to 1.5 kJ/mol
   more negative** than -RT ln(1/K_D) at the stated 25 C, and the offset is not a constant
   temperature (the implied T ranges over 310-319 K across the six rows). Possible explanations
   include a different reference state, an n-weighted constant, or a units slip; **none of them can
   be tested without SI part III**. **Refuse Table 4 and quote the constants.**
4. **A clear typographical error in the Methods concentrations.** "3534 g/L (150 µM) trypsin" — at
   trypsin's ~23 300 g/mol, 150 µM is **3.5 g/L**, so "3534" is almost certainly "3.534". It does not
   affect any result (trypsin gave no signal), but it is a warning about the care taken with the
   concentration list. The other three conversions all check out: BSA/HSA 10 g/L at 150 µM implies
   **66 700 g/mol** ✓; lysozyme 40 g/L at 2.74 mM implies **14 600 g/mol** ✓; SG 10 g/L at 29.1 µM
   implies **343 600 g/mol** ✓ (the glycinin hexamer).
5. **The "temperature" rows of Table 2 are NOT binding constants at temperature.** The protein was
   "incubated at 40, 60 and 90 °C for **10 min** and **cooled to room temperature before addition of
   quercetin and measurement**". So all three rows report the binding of **cold, previously heated**
   BSA. **There is no measurement anywhere in this paper of how a binding constant depends on the
   temperature at which the binding happens** — which is the quantity a thermal model needs. Anyone
   reading "T = 90 °C, K_D = 15 µM" as a hot-solution constant would be wrong.
6. **Two small internal inconsistencies.** (i) BSA's isoelectric point is given as **pH 5** in one
   paragraph and **pH 5.1** two paragraphs later. (ii) The CONTIN allocation for control BSA is
   quoted in the text as "**65 % α-helix, 3 % β-strand, 12 % β-turn, 20 % unordered**" while Table 3's
   pH 7 / 0 mM / 25 C row prints **66.9 / 2.5 / 12.0 / 18.6**. The text is evidently a rounded
   restatement, but they are not the same numbers. **Use Table 3.**
7. **The authors themselves say the error bars are too small, and it is the most important sentence
   in the paper.** "The K_D standard deviation generally underestimates the real uncertainties of
   the experiment, and every nonlinear least squares regression therefore actually represents the
   extent of the goodness of the curve fitting rather than the real error (16)." **Every ± in
   Tables 1 and 2 is a fit statistic on data that were already averaged over four repetitions before
   fitting** — so the replicate-to-replicate spread is never propagated at all. Treat all quoted
   precisions as lower bounds on the true uncertainty.
8. **Three of the seven Table 1 rows are compromised by the paper's own account.** Gallic acid
   "**was found to oxidize visibly**" at pH 7 and above; it had to be moved out of Hepes because it
   **binds the buffer**; a saturation curve could only be "tempted" by dropping the protein to 90 µM
   and "**was achieved only partly**"; at 30 µM protein and pH 3.5 the regression gave "a linear
   curve" rather than a saturation curve; and the authors write plainly that for gallic acid "the
   values of K_D are not optimal". The n = 14 ± 3 at pH 3.5 is called "indicative of an **unspecific**
   binding". **The three gallic-acid rows should be carried as a trend, not as three constants.**
9. **Three measured nulls, and they are useful.** (i) **Trypsin**: no enhancement of quercetin
   fluorescence, so no constant. (ii) **Rutin and isoquercetin**: no enhancement with **any** protein
   — so glycosylating quercetin at the 3-O position abolishes the fluorescence signal (not
   necessarily the binding: the near-UV CD in Fig. 5 shows rutin *does* perturb soy glycinin's
   tertiary structure, and isoquercetin perturbs SG but not BSA). **The method's blindness must not
   be read as an absence of binding.** (iii) Quercetin, rutin and isoquercetin **could not be run by
   HD at all** because they bind the column material — a direct demonstration that the two methods
   are not interchangeable.
10. **A covalent-binding artefact the authors found and reported honestly.** They tested whether
    *covalently* bound quercetin also enhances fluorescence at 526 nm, and found that "**the
    covalent-bound form is capable of producing a stronger effect on the emission as compared to the
    noncovalent-bound one**". So the fluorescence method **does not distinguish covalent from
    noncovalent binding by signal alone** — it distinguishes them by the urea control (Fig. 4). Any
    application of this method at a pH or temperature where covalent chemistry runs would conflate
    the two.
11. **Supporting Information parts I-IV are not on disk and several claims depend on them**: the
    Hummel-Dreyer principle and its calibration, all the chromatographic conditions per compound,
    the free-energy equation (Flags 3), the fluorescence spectra, the ΔF_max = 0.05 assumption for
    lysozyme, the urea controls for HSA/SG/lysozyme, the tryptophan-quenching evidence, the CDSSTR
    solvent artefact, and the near-UV CD of BSA. **Retrieve it before using Table 4 or the lysozyme
    row.**
12. **What this paper does not contain**: any odour threshold; any volatile or aroma compound; any
    rate constant; any binding constant measured at a temperature other than 25 C; any soy protein
    other than the glycinin fraction; any pea, lupin or other plant isolate; any covalent binding
    constant; any n for the fluorescence rows; any printed value for the ferulic-acid binding
    capacities of Figure 2.
13. **What to request**: (i) the Supporting Information, in full; (ii) the numeric values behind
    Figure 2 — eight proteins x two pH values of ferulic-acid binding capacity in mg/g, **including
    soy glycinin and gelatin**, which would be directly usable and are currently locked in a bar
    chart; (iii) site counts for the fluorescence rows, especially soy glycinin; (iv) the temperature
    at which Table 4's ΔG was computed.
14. **Registry gaps against `data/keys/compounds.yml`**: `chlorogenic_acid` is keyed and has four
    binding numbers here. **`ferulic_acid`, `gallic_acid`, `quercetin`, `rutin` and `isoquercetin`
    are all absent.** Ferulic acid is the more consequential omission for this repository: it is
    bound by **gelatin, milk whey and soy glycinin** in Figure 2, it is a ubiquitous cereal
    phenolic, and `4_vinylguaiacol` — its decarboxylation product and a keyed Maillard marker —
    is already in the registry, so the precursor is missing while the product is present.
