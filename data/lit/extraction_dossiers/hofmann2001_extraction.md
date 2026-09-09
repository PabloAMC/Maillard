# Hofmann & Schieberle 2002 — EXTRACTION (coffee brew and five model electrophile mixtures; 2-furfurylthiol 500 µg in 10 mL 0.1 mol/L phosphate pH 6.0 at 30 C for 30 min, and a real brew held at 80 C in a thermos for 210 min; SIDA quantification, [2H]-NMR and LC/MS structure work; the file on disk is named for the 2001 web date)

### THE PAPER `k_thioether` CITES PRINTS NO RATE CONSTANT AT ALL — and it does print, once, in one sentence, the one number this cluster was sent to look for: **400 µg of 2-furfurylthiol bound to 125 mg of coffee melanoidin**, which is 80 % of a 500 µg charge and, per gram of melanoidin, a measured thiol-binding **capacity of 0.028 mmol/g (mine)** — about **320 times smaller** than the 8-10 mmol/g titrated site density on which the `MELE` pool is sized. That ratio is a number the sink objective has never seen, and it points the opposite way from the B17 variant (a) fit, which was free to make sites and found the data indifferent.

**Source on disk:** `data/articles/hofmann2001.pdf` (8 pp., J. Agric. Food Chem. 2002, 50 (2), 319-326;
published on the web 7 December 2001, hence the file name). Read from the `pdftotext -layout` text
layer (`scratchpad/articles/hofmann2001.txt`, 397 lines). **Tables 1 and 2 came through clean** and
are re-typed in full below. Figures 1, 3 and 6 are the paper's three quantitative plots and are
**figure_only**: Figure 1 (thiol decay in the 80 C thermos brew), Figure 3 (FFT degraded by four
melanoidin fractions against their EPR radical activity) and Figure 6 (the 30 C time course of FFT
binding and of disulfide formation). Figures 2, 4, 5, 7, 8, 9 and 10 are NMR/MS spectra and reaction
schemes. There is no supplementary material. Repo status before this dossier: this paper is the
`source_anchor` of **`k_thioether`** in `src/kinetic_core/parameters_sulfur.py`, the authority for
the `MELE` species note in `src/kinetic_core/species_sulfur.py`, the source of the
`hofmann2002_brew_80C_FFT` hold-out row in `generate_kinetic_core_b2_holdout.py`, and a **FIT** and a
**★ HOLD-OUT** line in `docs/reference/FIT_HOLDOUT_DECLARATION.md` D.4. A **claim-by-claim
verification** of the code against this paper already exists on disk as
`data/lit/extraction_dossiers/hofmann2002_extraction.md` (written 2026-09-04, a different document
in a different format); **this dossier is the first full extraction of the paper itself** and it
re-types both tables, which that document did not.

## 0. Identity

| field | value |
|---|---|
| Title | "Chemical Interactions between Odor-Active Thiols and Melanoidins Involved in the Aroma Staling of Coffee Beverages" |
| Authors | Thomas Hofmann (corresponding) and Peter Schieberle — Deutsche Forschungsanstalt für Lebensmittelchemie, Lichtenbergstrasse 4, D-85748 Garching, Germany |
| Venue | J. Agric. Food Chem. 2002, 50 (2), 319-326. Received 25 June 2001, revised 11 October 2001, accepted 16 October 2001, **published on Web 07 December 2001** |
| DOI / article ID | 10.1021/jf010823n (printed as `JF010823N`) |
| The two names | The repository calls it "Hofmann 2002" (the print year, which is correct for citation) and the file is `hofmann2001.pdf` (the web year). They are one paper. |
| Naming | CROSSPY = 1,4-bis-(5-amino-5-carboxy-1-pentyl)pyrazinium radical cation; FFT = 2-furfurylthiol; MMBF = 3-mercapto-3-methylbutyl formate; MFT = 2-methyl-3-furanthiol; rFD = relative flavour dilution factor |
| Companions on disk | `hofmann2002_extraction.md` (the 2026-09-04 code-verification note on this same paper), `charlesbernard2005_extraction.md` (the titrated site density and the 25 C ladder), `stack2018` (the reversible quinone equilibrium the lane lumps into the same `MELE`), `kumazawa2003_extraction.md` (the FFT pH survival grid) |

## 1. Why it matters

**What it contributes to the thiol-sink question — and here there IS a number.** The paper prints,
on p. 325, one sentence of stoichiometry: "Although about 400 or 330 µg of 2-furfurylthiol was
bound to the coffee melanoidins or the pyrazinium-derived intermediates, respectively, less than
6 µg of the corresponding bis(2-furfuryl) disulfide was generated." The melanoidin charge for that
experiment is given in Methods as **125 mg in 10 mL** and the thiol charge as **500 µg**. Three
things follow, and only the first is already in the repository:

1. **80 % of the charged thiol is bound** (400 of 500 µg). This is the plateau the brief asked to be
   verified, and it verifies — but the arithmetic runs the other way round from the way it is
   usually stated. It is *400 µg bound, which is 80 % of the 500 µg charged*, not "80 % of 400 µg".
2. **The measured capacity is 0.028 mmol of thiol per gram of melanoidin (mine)**: 400 µg / 114.16
   g/mol = 3.50 µmol on 125 mg. The `MELE` pool is sized on Charles-Bernard's hydroxylamine
   titration, **8-10 mmol per g dry solids**, so the *carbonyl* density and the *thiol-consuming*
   density differ by a factor of about **320** in the same kind of matrix. That is not a
   contradiction — a hydroxylamine titration counts every carbonyl and only a small, radical-derived
   subset of them is a pyrazinium electrophile — but it is a **hard lower bound on a quantity the
   B17 variant (a) fit was allowed to invent**, and the fit's optimum (log10 yield 0.048, i.e. 1.1
   sites per osone) was never checked against it. A pot of 0.1 M sugar making 1.1 sites per osone
   would carry sites in the tens of millimolar; this paper's melanoidin at its own working
   concentration of 12.5 g/L carries **0.35 mmol/L** of thiol-consuming capacity.
3. **The disulfide branch is dead here**: < 6 µg of bis(2-furfuryl) disulfide, which consumes
   2 x 6/226.28 = 0.053 µmol = **6.1 µg of thiol (mine), i.e. 1.5 % of the 400 µg bound**. The
   authors' own reading: "neither the solution containing the diquaternary pyrazinium ions nor the
   coffee melanoidins are able to oxidize the thiol into its disulfide." This is the measured
   negative that B17 variant (b) reached by fitting.

**A thiol loss rate at a temperature other than 145 C: yes, two, and both are already in the
repository — but only one of them is a printed number.** At **80 C** in a real brew the text prints
"After 60 min the 2-furfurylthiol concentration decreased by a factor of more than four compared to
that of the fresh coffee brew" from an initial "about 16.0 ... µg", and complete loss by 210 min.
That is **k >= ln 4 / 60 min = 0.0231 /min = 3.85e-4 /s at 80 C (mine)** — the
`hofmann2002_brew_80C_FFT` hold-out row, and note that the printed wording makes it a **lower
bound**, not a point value. At **30 C**, Table 2's "17 % (15-19)" remaining after 30 min gives
**9.8e-4 /s (range 9.2e-4 to 1.05e-3, mine)**, which reproduces `k_thioether`'s quoted anchor
exactly. Both are already scored. What the sink objective has *not* seen is the capacity of point 2
and the fact in the next paragraph.

**What the code claims and what the paper actually prints — a correction that matters.** The
`k_thioether` entry is a **second-order** constant (`order=2`, k_ref 5.01e-4 L/(mmol·min)) whose
`source_anchor` reads "Hofmann & Schieberle 2002 Fig. 6 (9.4e-4 /s, 30 C, SIDA) and Table 2
(9.8e-4 /s)". **This paper prints no rate constant, no half-life, no reaction order and no
second-order constant anywhere.** Table 2 prints *relative amounts remaining after a fixed 30 min*;
Figure 6 is a plot. Every "/s" figure attributed to this paper is a first-order constant *derived*
from those by whoever read it, and the second order is supplied entirely by Charles-Bernard's site
density. The registry's `evidence_class="measured_rate"` on that row is therefore generous: what is
measured here is a **fractional loss over a fixed interval**, which is `within_study_ratio` in the
house scheme, and the rate is `derived_assumption` on top of an assumed first order. The number is
not in doubt — the pseudo-first-order arithmetic is elementary and reproduces — but the class is.

**What else the paper establishes, structurally.** (a) The binding is **covalent**: [2H]-NMR of
melanoidin pre-incubated with [2H2]-FFT shows a strong, line-broadened resonance at 3.0-4.2 ppm on
the macromolecule (Figure 2C), absent from both controls. (b) It is **not** reducible: earlier work
(ref 6) found dithioerythritol could not regenerate major amounts of free thiol, which excludes
disulfide linkage to cysteinyl residues. (c) The electrophile is a **pyrazinium** derived from the
CROSSPY radical, identified by LC/MS as thioether adducts at m/z 251, 267, 363 (and their [2H]
shifts to 253, 269, 367), and confirmed on a lysine-based model at m/z 537. (d) **Thiol binding
tracks EPR radical activity across four melanoidin size fractions** (Figure 3) — fraction IV,
highest radical activity, gave complete loss of FFT in 30 min at 30 C; fraction III, lowest radical
activity, the least binding. (e) The reaction is **general across thiols**: the same thioether
adducts form with MMBF (m/z 285, bis 431) and MFT.

What this paper does NOT give the repository: any rate constant; any activation energy; any reaction
order; any temperature series (30 C and 80 C are two different experiments in two different
matrices, and the paper never puts them on one axis); any pH other than 6.0 for the models; any
melanoidin concentration series; any measurement of how much melanoidin capacity remains after
binding; and any statement of the basis (per litre? per kg?) of the brew concentrations in Figure 1.

## 2. Methods as they matter to a model

- **The brew.** *Coffea arabica* var. caturra, Ecuador, medium roasted, colour value 12, ground in
  liquid nitrogen; **hot water at 95 C poured over 50 g powder per litre** in a paper filter. The
  volatile fraction was isolated from 100 mL of fresh brew by SAFE (ref 10).
- **The melanoidins.** 50 g of freshly ground powder extracted with **1 L of hot tap water at
  80-90 C**; defatted with dichloromethane; freeze-dried, **yield 12.5 g**. 1.25 g aliquots were
  redissolved in 20 mL water and either (i) ultrafiltered through a **3000 Da** cut-off (Diaflo YM3)
  to give **0.44 g of melanoidins** after freeze-drying — so the >3000 Da fraction is **35 % of the
  water extract (mine)** — or (ii) separated on Sephadex G-25 fine (75 x 5 cm) into four fractions,
  **I 258 mg, II 221 mg, III 570 mg, IV 141 mg** (1190 mg recovered of 1250 mg charged, 95 %, mine).
- **The headspace binding experiment (Table 2, and the 400 µg sentence).** **2-Furfurylthiol 500 µg
  in 10 mL of 0.1 mol/L phosphate at pH 6.0**, in a **septum-sealed 240 mL vessel**, **stored at
  30 C**, alone or with one of: chlorogenic acid 20 mg; chlorogenic acid pre-heated 5 min at 230 C;
  albumin + glucose 10 mg each, dry-heated 5 min at 230 C; albumin + glycolaldehyde 10 mg each,
  same; Nα-acetyl-L-lysine + glycolaldehyde 10 mg each, same; **or coffee melanoidins 125 mg**.
  **Triplicate**, and Table 2 prints the range. The charge is therefore **438 µmol/L FFT and
  12.5 g/L melanoidin** (mine, and the same pair the registry's `_THIOETHER_CONDITIONS` string
  carries).
- **The comparative AEDA (Table 1).** The total volatile fraction from 10 mL of brew, alone or
  remixed with **125 mg of melanoidins (MW > 3000 Da)** — "their 'natural' concentrations" — in
  10 mL, incubated 30 min, then static headspace: stepwise reduced injection volumes from **25 mL
  down to 0.1 mL**, giving relative flavour dilution factors **1 to 256**, onto a 60 m x 0.32 mm
  RTX-5 (3 µm film), 30 -> 230 C at 6 C/min, to a sniffing port or an Incos XL MS. **The incubation
  temperature is printed three different ways in this paper** — Table 1 footnote a says "incubated
  for 30 min at 30 C", the Results text says "a stored model (30 min; 40 C)", and the Static
  Headspace Analysis method says the vessel "were equilibrated for 30 min at 45 C". See Flags 3.
- **Quantification.** **Stable isotope dilution assay** with [2H2]-FFT, [2H8]-3-methyl-2-butenethiol,
  [2H6]-MMBF and [2H4]-bis(2-furfuryl)disulfide as internal standards (ref 11); ion-trap MS (ITD 800)
  in **chemical ionisation with methanol** as reactant gas. This is the strongest quantification
  method in this cluster and it is what makes the 400 µg and < 6 µg figures usable.
- **The thermos experiment (Figure 1).** An original coffee brew **kept at 80 C in a thermos flask**,
  sampled at **0, 30, 60, 90 and 210 min**, FFT and MMBF by SIDA.
- **The labelling experiment.** Total melanoidins 100 mg in 10 mL tap water + **[2H2]-FFT 1.0 mg**,
  **90 min at 30 C** in a closed vessel, then ultrafiltered (>3000 Da), retentate taken up in 1 mL
  water, **[2H] NMR at 500 MHz**. Controls: melanoidin alone; [2H2]-FFT alone.
- **The model electrophile.** 1,4-Diethyl pyrazinium diquaternary salt **1.0 mg in 2 mL tap water**
  with **500 µg** of FFT, MMBF or MFT added individually, **30 min at 30 C**, then LC/MS (LCQ, ESI,
  direct injection). Separately, glycolaldehyde 20 mg + Nα-acetyl-L-lysine 20 mg + 100 µL tap water
  **heated in an open beaker 5 min at 230 C** to generate CROSSPY (checked by EPR and LC/MS on a
  10 mg aliquot in 3 mL water), then FFT 500 µg, 30 min at 30 C, LC/MS and LC/MS².

## 3. Tables re-typed

### Table 1. "Results of Comparative Aroma Dilution Analyses of the Headspaces of Isolated Coffee Brew Volatiles Incubated Either in the Absence (I) or Presence of Coffee Melanoidins (II)^a"

Footnote a exactly as printed: "Aqueous solutions (10 mL) containing the total volatile fraction
isolated from a freshly prepared coffee brew (10 mL) were incubated for 30 min at 30 °C. Model II
contained 125 mg of melanoidins (molecular weight > 3000 Da)."

| odorant | aroma quality | rFD factor I | rFD factor II |
|---|---|---:|---:|
| butane-2,3-dione | buttery | 256 | 128 |
| pentane-2,3-dione | buttery | 128 | 128 |
| 3-methylbutanal | malty | 64 | 64 |
| 2-methylbutanal | malty | 32 | 64 |
| acetaldehyde | fruity | 32 | 32 |
| methional | potato-like | 32 | 16 |
| **2-furfurylthiol** | roasty, sulfury | **32** | **2** |
| 2-ethyl-3,5-dimethylpyrazine | earthy | 32 | 32 |
| 2,3-diethyl-5-methylpyrazine | earthy | 32 | 32 |
| 2-methoxyphenol | phenolic | 16 | 32 |
| dimethyl trisulfide | cabbage-like | 16 | 32 |
| 2-isobutyl-3-methoxypyrazine | green, earthy | 16 | 16 |
| **3-methyl-2-butenthiol** | foxy, skunky | **8** | **1** |
| **3-mercapto-3-methylbutyl formate** | catty | **8** | **2** |
| **2-methyl-3-furanthiol** | meatlike | **4** | **2** |
| **methanethiol** | cabbage-like | **2** | **<1** |

**rFD drops (mine, I/II):** FFT **16x**, 3-methyl-2-butenthiol **8x**, MMBF **4x**, MFT **2x**,
methanethiol **> 2x** (below detection in II), butane-2,3-dione 2x, methional 2x. Two non-thiols go
**up** 2x (2-methoxyphenol, dimethyl trisulfide) and one (2-methylbutanal) goes up 2x; the pyrazines
and the C4/C5 diones are otherwise unchanged. The text's summary: "the aroma impacts of odorants
belonging to other chemical classes, such as the 2,3-diones, the phenols, or the pyrazines, were not
significantly changed." **An rFD factor is a dilution step on a factor-of-two ladder; it is not a
concentration and a 16x rFD drop is not a 16x concentration drop** (Flags 5).

### Table 2. "Relative Amounts of 'Free' 2-Furfurylthiol Present in the Headspaces of Aqueous Solutions of 2-Furfurylthiol Stored in the Presence of Different Model Mixtures"

Footnote a: "Compounds were intimately mixed and dry-heated for 5 min at 230 °C." Footnote b:
"Relative amount of 'free' thiol is given as the mean of triplicates. Ranges of data measured are
given in parentheses."

| 2-furfurylthiol stored in the presence of | rel. amount of 2-furfurylthiol^b |
|---|---|
| no additive (control) | 100 |
| chlorogenic acid (20 mg) | 92 (88−94) |
| thermally pretreated chlorogenic acid (20 mg)^a | 86 (82−90) |
| albumin/glucose (10 mg each)^a | 58 (55−61) |
| albumin/glycolaldehyde (10 mg each)^a | 31 (28−34) |
| Nα-acetyl-L-lysine/glycolaldehyde (10 mg each)^a | 17 (15−19) |

The conditions belong to the Methods block quoted in section 2: **500 µg FFT in 10 mL 0.1 mol/L
phosphate pH 6.0, 30 min at 30 C, 240 mL septum-sealed vessel.** Note that **coffee melanoidin
itself is not a row of Table 2** — its 30 min point lives only in Figure 6 and in the 400 µg
sentence.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| FFT in the fresh brew | "about **16.0** ... µg" | p. 322 (basis not printed — see Flags 4) |
| MMBF in the fresh brew | "about ... **8.2** µg" | p. 322 |
| FFT after 60 min at 80 C | "decreased by a **factor of more than four**" | p. 322 |
| FFT after 210 min at 80 C | "**complete loss** of 2-furfurylthiol"; "only small amounts of 3-mercapto-3-methylbutyl formate were detectable" | p. 322 |
| **FFT bound to coffee melanoidins** | "**about 400** µg" | p. 325 |
| **FFT bound to the pyrazinium intermediates** | "**about 330** µg" | p. 325 |
| **bis(2-furfuryl) disulfide generated** | "**less than 6 µg**" | p. 325 |
| chlorogenic acid, both forms | "the losses in both models were **below 20 %**" | p. 323 |
| albumin/glucose | thiol "decreased by a factor of nearly two" | p. 323 |
| albumin/glycolaldehyde | thiol "reduced to below 30 %" | p. 323 |
| Nα-acetyl-lysine/glycolaldehyde | "only **17 %** of the 'free' 2-furfurylthiol was left" | p. 323 |
| melanoidin GPC fraction IV | "**complete loss** of 2-furfurylthiol after 30 min at 30 C"; highest EPR radical activity | p. 323 and Figure 3 |
| melanoidin GPC fraction III | "least effective in 2-furfurylthiol binding"; lowest radical activity | p. 323 |
| melanoidin and pyrazinium kinetics | "showed **similar kinetics** of 2-furfurylthiol degradation" | p. 325 |
| [2H2]-FFT alone, [2H] NMR | 3.67 ppm (deuterated methylene) and 4.70 ppm (natural [2H] in tap water) | p. 322, Fig. 2A |
| melanoidin + [2H2]-FFT, [2H] NMR | "additional strong resonance between **3.0 and 4.2 ppm**", strongly line-broadened | p. 322, Fig. 2C |
| pyrazinium hydrolysis products (ESI m/z) | CROSSPY-type radical cation **138**; 2-hydroxy-1,4-diethyl-1,4-dihydropyrazine **155**; dihydroxy **171**; bis-hydroxy dimer **309** | p. 325, Fig. 5 |
| FFT + pyrazinium adducts (ESI m/z) | **251** (base, 100 %, mono-thioether), **267** (hydroxy-thioether), **363** (bis-thioether); [2H2] shifts to **253**, **269**, **367** | p. 325, Fig. 7 |
| MMBF + pyrazinium adducts | **285** (base, 100 %), **431** (bis) | p. 325-326, Fig. 8A |
| MFT + pyrazinium | "the expected thioether derivatives" (no m/z printed in the text) | p. 326, Fig. 8B |
| lysine/glycolaldehyde + FFT | quasi-molecular ion **537**; MS² base ions **424** (loss of 113) and **455** (loss of 82) | p. 326, Fig. 9 |
| melanoidin freeze-dry yield | 12.5 g from 50 g coffee powder | Methods |
| >3000 Da fraction | 0.44 g from a 1.25 g aliquot | Methods |
| GPC fractions | I 258 mg, II 221 mg, III 570 mg, IV 141 mg | Methods |

**Figure-only in this paper:** every point of Figure 1 (FFT and MMBF against 0/30/60/90/210 min at
80 C) except the printed 16.0, 8.2, "factor of more than four" and "complete loss"; every point of
Figure 3 (µg FFT degraded by fractions I-IV against relative radical activity); every point of
Figure 6 (the 30 C time courses of FFT loss and of disulfide formation, melanoidin and pyrazinium)
except the printed 400 / 330 / < 6 µg end-points. Per house rule they are not typed as numbers here.

### Arithmetic on the printed numbers (all mine)

**1. Pseudo-first-order constants from Table 2**, assuming first order over the fixed 1800 s and
that the headspace tracks the solution:

| additive | fraction left | k (1/s) | k (1/min) |
|---|---:|---:|---:|
| chlorogenic acid | 0.92 | 4.6e-5 | 2.8e-3 |
| heated chlorogenic acid | 0.86 | 8.4e-5 | 5.0e-3 |
| albumin/glucose | 0.58 | 3.0e-4 | 1.8e-2 |
| albumin/glycolaldehyde | 0.31 | 6.5e-4 | 3.9e-2 |
| **Nα-acetyl-lysine/glycolaldehyde** | **0.17 (0.15-0.19)** | **9.8e-4 (9.2e-4 to 1.05e-3)** | **5.9e-2** |

The last row reproduces `k_thioether`'s quoted 9.8e-4 /s to two figures. The ladder spans **21x**
from chlorogenic acid to the lysine model — one condition, one temperature, one interval, and
therefore a clean **within_study_ratio** on *which electrophile precursor matters*, which is the
paper's real contribution: **a lysine + glycolaldehyde melanoidin is 21x more thiol-consuming than
chlorogenic acid and 3.3x more than an albumin + glucose melanoidin, at equal mass.**

**2. The capacity, three ways.** 400 µg FFT bound / 114.16 g/mol = **3.50 µmol** bound.
- per gram of melanoidin: 3.50 µmol / 0.125 g = **28.0 µmol/g = 0.0280 mmol/g**
- per litre of the model solution: 3.50 µmol / 0.010 L = **0.350 mmol/L**
- as a fraction of charge: 400/500 = **80 %**, leaving 100 µg = **87.6 µmol/L** free.
Against Charles-Bernard's hydroxylamine-titrated **8-10 mmol/g**, the thiol-consuming density is
**0.28-0.35 % of the titrated carbonyl density**, i.e. a factor of **290-360 (mine)**. Note this is
a **lower bound on capacity**: 20 % of the thiol was still free at the plateau, so the pool may not
have been exhausted, and the experiment cannot distinguish an exhausted pool from an equilibrium.

**3. The disulfide share.** < 6 µg of bis(2-furfuryl) disulfide, MW 226.28, = **< 0.0265 µmol**,
consuming **< 0.053 µmol = < 6.1 µg** of thiol. Against the 400 µg bound that is **< 1.5 %**; against
the 500 µg charged, **< 1.2 %**. Confirms `parameters_sulfur.py`'s note verbatim.

**4. The 80 C brew constant, and why it is a bound.** "A factor of more than four" in 60 min gives
**k >= ln(4)/60 min = 0.0231 /min = 3.85e-4 /s**; "complete loss" by 210 min is consistent with
anything above about 0.02 /min but places no upper bound because "complete" is a detection limit,
not a zero. The row `hofmann2002_brew_80C_FFT` carries 0.023 /min; **it should carry a `>=`.**

**5. The cross-temperature comparison the paper invites and the repository forbids.** The 30 C model
gives 5.9e-2 /min and the 80 C brew **at least** 2.3e-2 /min. Taken naively that is a **negative**
apparent activation energy. The registry's `PROHIBITED_DERIVATIONS` already names this derivation
and forbids it, on the ground (from `k1_kinetic_parameters.md` §1d) that the real brew's electrophile
pool was partly consumed during extraction. **This paper supplies the direct evidence for that
explanation and not merely the anomaly**: Figure 3 shows binding tracks *radical activity*, and the
radical is generated by roasting and consumed by its own redox cycle (Figure 4), so a brew that has
stood is a matrix with fewer sites than a freshly reconstituted melanoidin. The prohibition stands
and this dossier strengthens its stated reason.

**6. What the 330 µg pyrazinium figure cannot be turned into.** The mass of diquaternary salt used
in the Figure 6 experiment is **not printed** (Methods gives 1.0 mg in 2 mL for the LC/MS work,
which is a different vessel and a different volume). So there is no per-mole site capacity for the
pyrazinium arm, only the observation that it binds nearly as much thiol as 125 mg of melanoidin.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** Keyed: `2_furfurylthiol`, `2_methyl_3_furanthiol`,
`methanethiol`, `methional`, `dimethyl_trisulfide`, `chlorogenic_acid`, `2_3_butanedione`,
`3_methylbutanal`, `2_methylbutanal`, `acetaldehyde`, `2_ethyl_3_5_dimethylpyrazine`,
`3_isobutyl_2_methoxypyrazine`, `bis_2_methyl_3_furyl_disulfide`. **Not keyed:** 3-mercapto-3-
methylbutyl formate, 3-methyl-2-butenthiol, bis(2-furfuryl) disulfide, pentane-2,3-dione,
2,3-diethyl-5-methylpyrazine, 2-methoxyphenol, and the melanoidin/CROSSPY/pyrazinium species (which
the model carries as the atom-free `MELE` site pool). See Flags 8.

Rows below marked "30 C model" share: **FFT 500 µg in 10 mL of 0.1 mol/L phosphate at pH 6.0, 30 C,
30 min, 240 mL sealed vessel, triplicate, static headspace with SIDA.** Rows marked "80 C brew"
share: **filter brew from 50 g powder/L, thermos flask, 80 C, SIDA.**

| step / quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|
| FFT remaining, no additive | 100 | % of control | 30 C model | none printed | Table 2 | within_study_ratio |
| FFT remaining, chlorogenic acid 20 mg | 92 (88-94) | % | 30 C model | none | Table 2 | within_study_ratio |
| FFT remaining, heated chlorogenic acid 20 mg | 86 (82-90) | % | 30 C model | none | Table 2 | within_study_ratio |
| FFT remaining, albumin/glucose 10 mg each | 58 (55-61) | % | 30 C model | none | Table 2 | within_study_ratio |
| FFT remaining, albumin/glycolaldehyde 10 mg each | 31 (28-34) | % | 30 C model | none | Table 2 | within_study_ratio |
| **FFT remaining, Nα-acetyl-lysine/glycolaldehyde 10 mg each** | **17 (15-19)** | % | 30 C model | none | Table 2 | **within_study_ratio** — the anchor `k_thioether` is derived from |
| pseudo-first-order constant, lysine/glycolaldehyde model | **9.8e-4 (9.2e-4 to 1.05e-3)** | 1/s | 30 C, pH 6.0 | **first order assumed by me**; the paper states no order | derived from Table 2 (mine) | **derived_assumption** (NOT measured_rate — see Flags 1) |
| electrophile-precursor ladder | 21x from chlorogenic acid to lysine/glycolaldehyde; 3.3x albumin/glucose to lysine/glycolaldehyde | — | 30 C model, equal mass | — | derived from Table 2 (mine) | within_study_ratio |
| **FFT bound to coffee melanoidin** | **~400** | µg, from a 500 µg charge on 125 mg melanoidin in 10 mL | 30 C, pH 6.0, 12.5 g/L melanoidin | mass balance, not a rate | p. 325 | **fed_intermediate_yield** (a bound amount from a fed thiol) |
| **melanoidin thiol-binding capacity** | **>= 0.0280** | mmol thiol per g melanoidin (>3000 Da) | as above | — | derived from p. 325 + Methods (mine) | **derived_assumption** — a LOWER bound; 20 % of the thiol was still free |
| same, per volume | >= 0.350 | mmol/L at 12.5 g/L | as above | — | derived (mine) | derived_assumption |
| FFT bound to pyrazinium intermediates | ~330 | µg | 30 C, mass of salt **not printed** | — | p. 325 | fed_intermediate_yield (unusable per gram — Flags 6) |
| **bis(2-furfuryl) disulfide formed** | **< 6** | µg | 30 C, alongside 400 µg bound | — | p. 325 | **threshold** (an upper bound) |
| disulfide share of the thiol flux | **< 1.5** | % | 30 C | — | derived (mine) | within_study_ratio |
| FFT in fresh brew | ~16.0 | µg (basis not printed) | 80 C brew, t = 0 | — | p. 322 | level_only |
| MMBF in fresh brew | ~8.2 | µg (basis not printed) | 80 C brew, t = 0 | — | p. 322 | level_only |
| FFT loss at 80 C | **>= 0.0231** | 1/min | 80 C brew, 0-60 min | first order assumed by me | derived from p. 322 (mine) | derived_assumption (a **lower** bound) |
| FFT at 210 min, 80 C | complete loss | — | 80 C brew | — | p. 322 | threshold |
| rFD drop on adding melanoidin | FFT 16x, 3-methyl-2-butenthiol 8x, MMBF 4x, MFT 2x, methanethiol >2x; pyrazines and diones unchanged | dilution steps | 10 mL brew volatiles ± 125 mg melanoidin, 30 min | — | Table 1 | **threshold** (rFD is an odour-dilution factor, not a concentration — Flags 5) |
| melanoidin fraction IV | complete FFT loss in 30 min at 30 C; highest EPR radical activity | — | GPC fraction, 30 C | — | p. 323 + Fig. 3 | threshold / figure_only for the plot |
| melanoidin fraction ordering | binding tracks radical activity: IV > I ≈ II > III | — | 30 C | — | Fig. 3 + p. 323 | **figure_only** |
| covalent binding, [2H] NMR | resonance 3.0-4.2 ppm on the >3000 Da retentate after 90 min at 30 C | ppm | 100 mg melanoidin + 1.0 mg [2H2]-FFT in 10 mL | — | Fig. 2C | level_only (a structural observation, no amount) |
| time courses at 30 C (FFT and disulfide, both matrices) | — | µg vs min | 30 C model | — | Figure 6 | **figure_only** |
| time course at 80 C (FFT and MMBF) | — | µg vs min | 80 C brew | — | Figure 1 | **figure_only** |

### Can these be put on the same basis as the shipped constant? Step by step.

- **The 30 C fractional losses can, and they already are.** The pseudo-first-order arithmetic is
  sound *given* first order, and the SIDA quantification is the best in this cluster. What must
  travel with the number is that (i) the paper prints no order, so first order is the reader's
  assumption; (ii) the rate silently contains the electrophile concentration of one 12.5 g/L coffee
  melanoidin at pH 6.0, which is exactly why the module carries `MELE` explicitly.
- **The capacity is new to the objective and it can be used as a bound, not as a fit target.**
  0.0280 mmol/g at 12.5 g/L is a *lower* bound on the site density of a real, fully browned coffee
  melanoidin. Any wave that gives the lane a site pool should be required to state what site
  concentration its yield implies in the reference pot and to compare it with **0.35 mmol/L in a
  matrix that is already 12.5 g/L of finished melanoidin**. B17 variant (a) implied roughly 1.1
  sites per osone with no such comparison made.
- **The 80 C row cannot be put on an Arrhenius line with the 30 C row and the registry already says
  so.** This dossier adds the mechanism for why (radical-derived sites are consumed by their own
  redox cycle, Figure 4) and notes that the 80 C figure is a lower bound, which makes the negative
  apparent barrier *worse*, not better.
- **Nothing here transports to 145 C.** The highest temperature at which a thiol loss is measured in
  this paper is 80 C. The 230 C in the Methods is a dry-heating step that *makes* the electrophile;
  no thiol is present at that temperature.

## 5. Flags

1. **The paper prints no rate constant, no order and no barrier.** Every "/s" figure the repository
   attributes to it is derived by assuming first order over a fixed interval, and the second order in
   `k_thioether` comes entirely from Charles-Bernard's site density, not from here. The registry's
   `evidence_class="measured_rate"` on that row overstates what this source supports; the honest
   classes are `within_study_ratio` for Table 2 and `derived_assumption` for the constant. The
   *value* is not in dispute.
2. **The 400 µg is a single sentence with no error bar, no replicate count and no plateau
   demonstration.** It is the paper's only stoichiometric statement about binding and the whole
   capacity argument in section 3 rests on it. The three figures behind it (Figure 6's melanoidin
   curve, its pyrazinium curve, its disulfide curve) are plots. **Whether the melanoidin was
   saturated is not established**: 20 % of the thiol remained free, which is equally consistent with
   an exhausted pool and with a stalled or equilibrium-limited reaction.
3. **The incubation temperature of the Table 1 experiment is printed three different ways**: 30 C
   (Table 1 footnote a), 40 C (Results, p. 321) and 45 C (Methods, Static Headspace Analysis).
   These may be two different steps — incubation and headspace equilibration — described loosely,
   but as printed they conflict. Table 1's own footnote is the most specific and is what this
   dossier's table reproduces. Do not build a temperature axis on Table 1.
4. **The brew concentrations have no basis printed in the text layer**: "the freshly prepared coffee
   beverage contained about 16.0 or 8.2 µg of FFT or MMBF". Per litre, per kilogram, or in the
   sampled aliquot is not said in the running text; the Figure 1 axis presumably carries it and the
   figure is an image. The existing `hofmann2002_extraction.md` records µg/kg. **Verify from the
   figure before any absolute level from this paper enters a benchmark**; the *ratio* and the decay
   shape are unaffected.
5. **rFD factors are not concentrations.** Table 1's values are steps on a factor-of-two dilution
   ladder read by a human nose; a drop from 32 to 2 means the odour survived four fewer halvings,
   which bounds a concentration drop only loosely and is confounded by the odour threshold of the
   compound in that matrix. The paper's own abstract says FFT's *concentration* "decreased by a
   factor of 16", which is an over-reading of its own Table 1. Use Table 1 for the **ordering**
   (thiols affected, pyrazines and diones not) and never as a level.
6. **The pyrazinium arm of the Figure 6 experiment has no stated charge**, so its 330 µg cannot be
   turned into a capacity per mole of electrophile. This is the number that would let the lane size
   `MELE` in molecular rather than mass units, and it is one question to the authors away.
7. **What to request from the authors**: (i) the numeric data behind Figures 1, 3 and 6, especially
   the melanoidin curve's t = 0 and its late points, which would settle whether the plateau is a
   capacity or an equilibrium; (ii) the mass of diquaternary salt in the Figure 6 experiment;
   (iii) the basis (per L or per kg) of the Figure 1 concentrations; (iv) a melanoidin concentration
   series at 30 C, which would give the reaction order in sites directly and is the single
   experiment that would convert this paper's anchor from a pseudo-first-order constant into a real
   bimolecular one; (v) whether the residual 20 % free thiol is released on further melanoidin
   addition.
8. **Registry gaps against `data/keys/compounds.yml`**: **3-mercapto-3-methylbutyl formate** and
   **bis(2-furfuryl) disulfide** are the two that matter — MMBF is the second thiol in the 80 C
   hold-out row and the FFT disulfide is the species the "< 6 µg" bound is stated on, and the lane
   models both. Also absent: 3-methyl-2-butenthiol, pentane-2,3-dione, 2,3-diethyl-5-methylpyrazine,
   2-methoxyphenol.
9. **What this paper does not contain**: any melanoidin concentration series; any pH series (the
   models are all pH 6.0 and the brew is unbuffered and unmeasured); any temperature series on a
   single matrix; any activation energy; any measurement of the melanoidin's remaining capacity after
   binding; any regeneration experiment of its own (the dithioerythritol negative is cited to ref 6,
   not measured here); any statement of whether the 20 % of thiol left free at the plateau is
   releasable; and any supplementary material.
