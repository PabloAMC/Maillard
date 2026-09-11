# Coukos 2021 — EXTRACTION (methylglyoxal + glutathione or N-acetylcysteine + arginine or metformin, **1 mM each in PBS**, time courses at **25 °C** and dose–response at **37 °C for 24 h**, by LC–MS as **relative peak areas only**; plus a model peptide, BSA proteomics, HeLa metabolomics and an MRP1 knockdown)

### THE HEADLINE: this paper prints **NO rate constant, no reaction order and no activation energy for anything** — its "kinetic characterization" is a set of relative-abundance time courses and one **logistic fit whose parameters are never reported**. What it does supply is a **within-study kinetic ordering** on the adduct branch (thiol consumed faster than guanidine; hemithioacetal fastest, then declining after **2–6 h at 25 °C** as the thermodynamically favoured imidazolone and mercaptomethylimidazole crosslink accumulate) and a third irreversible **carbon**-sink for the thiol. **The MICA crosslink is a pure double condensation, not a redox step:** its printed search mass, **+210.1116 on cysteine**, equals **methylglyoxal + arginine − 2 H₂O = 210.1117 (mine)** to 1 × 10⁻⁴ Da. The sulfur ends up bonded to a ring carbon. **No disulfide is formed, measured or discussed anywhere in this paper.**

**Source on disk:** `data/articles/coukos2021.pdf` (9 pp., ACS Chem. Biol. 2021, 16, 2453−2461; open access, © 2021 The Authors, published by ACS).
Read from the `pdftotext -layout` text layer, which came through clean apart from the rotated ACS download banner in the left margin.
**This paper prints NO TABLES in the main text.** Its two tables — **extended data Table S1** (the MRM parameters for the metabolomics) and **extended data Table S2** (cloning primers) — are in the **Supporting Information, which is NOT on disk** (`data/articles/` holds `coukos2021.pdf` only). **Figures S1–S4 are also off disk**: S1 (structures of MGO-derived modifications and the model compounds), S2 (the 24 h screen showing all four reaction pairs form imidazolone, hemithioacetal and MICA), S3 (targeted fragmentation structures / synthetic standards), S4 (the shABCC1 metformin data). Figures 1–4 are in the PDF and are all plots, chromatograms, blots or structural depictions; per house rule their axis values are recorded as figure-only.
**Repo status before this dossier:** Coukos 2021 is cited nowhere in `src/kinetic_core/`, nowhere in `data/lit/reaction_rules.yml`, and has no extraction dossier. The term "mercaptomethylimidazole" appears nowhere in the repository.

## 0. Identity

| field | value |
|---|---|
| Title | "Methylglyoxal Forms Diverse Mercaptomethylimidazole Crosslinks with Thiol and Guanidine Pairs in Endogenous Metabolites and Proteins" |
| Authors | John S. Coukos, Raymond E. Moellering (corresponding, rmoellering@uchicago.edu, ORCID 0000-0002-2043-7838) — Department of Chemistry, The University of Chicago, Chicago, Illinois 60637, USA |
| Venue | **ACS Chemical Biology 2021, 16, 2453−2461.** Received 15 July 2021; accepted 16 September 2021; published 28 September 2021 |
| DOI | **10.1021/acschembio.1c00553** — printed in the page footer of every page as `https://doi.org/10.1021/acschembio.1c00553`, and again in the "Complete contact information" block as `https://pubs.acs.org/10.1021/acschembio.1c00553` |
| Funding / interest | NIH MSTP training grant T32GM007281 (J.S.C.); NSF-CAREER CHE-1945442 (R.E.M.); Alfred P. Sloan Foundation FG-2020-12839 (R.E.M.). "The authors declare no competing financial interest." |
| Contributions | "J.S.C. performed **all** experiments. R.E.M. supervised research." |
| Computational content | **NONE — and this is worth stating explicitly because the brief anticipated otherwise.** There is no DFT, no quantum chemistry, no molecular dynamics, no docking and no calculated energy anywhere in this paper. The only structure-derived object is Figure 2D, a **depiction** of the observed BSA modification sites drawn on the crystal structure **PDB 4F5S** — a rendering, not a calculation. **Nothing was excluded under the standing no-DFT policy, because there was nothing to exclude.** |
| The α-dicarbonyl | **methylglyoxal (MGO)**, synthesised in-house by acid hydrolysis of MG-1,1-dimethylacetal and purified by fractional distillation; concentration assigned by the aminoguanidine/320 nm colorimetric assay against a 3-amino-1,2,4-triazine calibration curve |
| The thiols | **glutathione (GSH, reduced)** and **N-acetylcysteine (NAC)**; also **free cysteine** (in cells), the **CRV2 model peptide** (synthesised in-house), and **BSA Cys58**, "the only nondisulfide cysteine in BSA" |
| The guanidines | **arginine** and the biguanide drug **metformin** |
| Naming | **MICA** = mercapto**m**ethyl**i**midazole **c**rosslink between **c**ysteine and **a**rginine (the paper's own coinage, from Bollong et al., ref 20); **MG-H1** = the hydroimidazolone arginine adduct; **CEA** / **CEL** = 1-carboxyethyl-arginine / -lysine; **Klac** = D-lactoyllysine |

## 1. Why it matters

**For `results/validation/kinetic_core_b27_prereg.md`, this paper is a third adduct-branch paper, and
it is the one that gives the repository the least arithmetic and the most structure.**

**(1) It supplies no constant of any kind, and the brief's question about it can be answered flatly.**
The word "kinetic" appears throughout, but every kinetic statement in this paper is either an
**ordering** ("more rapid consumption of the thiol-containing reactants than the guanidine-containing
reactants") or a **time to a turning point** ("after 2−6 h … the level of hemithioacetals began to
decrease"). The only fit anywhere is the **"logistic fit of MICA formation"** in Fig. 2B, and
**no fitted parameter — no rate, no half-time, no plateau — is printed for it.** The Methods say
outright what the *y* axis is: "Relative metabolite abundance was quantified by the integrated peak
area for each extracted ion chromatogram … **normalized to the most abundant peak in a given time- or
dose-response series**." **There is not one absolute concentration in the entire in vitro dataset,
and therefore not one number from which a rate constant could be recovered.** The paper says so of
itself in the Discussion: "relative quantification of these metabolites was used in this study
primarily to gauge kinetics and dose−response profiles … it is important that the **absolute levels**
… be determined in future work."

**(2) It adds a THIRD irreversible carbon-sink for the thiol, and that matters to a wave arguing
about sinks.** The repository has refused three sulfur-sink structures (B17(b) reversible disulfide,
B17(a) saturable thioether, B25 irreversible addition to deoxypentosones). Coukos describes a fourth
mechanism, and it is neither of the two the repository has tried: the reversible hemithioacetal, once
formed, is **captured by a guanidine and dehydrated twice** into a stable heteroaromatic
mercaptomethylimidazole. Sulfur is then bonded to an imidazole ring carbon, permanently. **A thiol
that ends up in a MICA crosslink can never make a disulfide.** In a pot that contains both a thiol
and an α-dicarbonyl and a guanidine, this is a competitor to *both* branches of the B27 argument.
**But it has a hard prerequisite: a guanidine.** The pot the wave must reproduce
(`fed_nf_cys_MFT` / `whitfield_nf_cys_MFT`) is charged with **cysteine** and norfuraneol and has no
arginine and no biguanide. **So the MICA route is structurally unavailable in the pot that must show
the 35 % disulfide share** — which is a genuine, if narrow, point *in favour* of leaving the sulfur
lane's competitor set as it is. It is available in any pot carrying protein (arginine side chains),
which most of the repository's fed systems do.

**(3) It sharpens the branch question the pre-registration has to keep straight, by mass arithmetic.**
The proteomics search parameters print exact mass shifts, and they settle the redox question for
every product this paper observes without any need to trust a mechanism drawing:

- **MICA on cysteine: +210.1116.** MGO (72.02113) + arginine (174.11168) − 2 H₂O (36.02113) =
  **210.1117 (mine)**. A **double condensation**. No electrons move.
- **MG-H1 on arginine: +54.0106.** MGO − H₂O = **54.0106 (mine)**. A **single condensation**.
- **CEL on lysine and CEA on arginine: +72.0211**, which is MGO **with no water lost** — a net
  addition, consistent with the intramolecular hydride shift these adducts are known to require, but
  the paper neither says so nor measures it, so that reading is mine and is not load-bearing here.

**None of these is a thiol oxidation.** The repository's instruction to keep the adduct and redox
branches apart is served by this paper better than by any other in the batch, because the separation
here is arithmetic rather than mechanistic.

**(4) It corroborates, at a third laboratory, the shape that Zheng 2022 and Zheng 2023 report.**
Thiol capture is fast and reversible; the stable products are slow. Coukos's Fig. 1C shows the
hemithioacetal rising first and then **falling after 2–6 h at 25 °C** while MICA and imidazolone
rise; Zheng 2023 shows the GSH-MGO adduct falling after **8 h at 37 °C** as kaempferol takes the MGO;
Zheng 2022 shows a plateau at 6 h. Three independent groups, three chemistries, same shape. **The
thiol wins the race and loses the war** is now a corpus-wide observation and not a single paper's
claim.

**What this paper does NOT do for the wave.** No rate constant, order or barrier — for either branch.
No disulfide, anywhere. Nothing above 37 °C. Nothing outside PBS at physiological pH. No α-diketone
(MGO again has the aldehyde carbon that both the hemithioacetal and the MICA route need). No
norfuraneol and no route from any substrate to an α-dicarbonyl, so pre-reg §3(a) is untouched. And
because everything is normalised within a series, **not even a within-study ratio between two
different species can be taken from it** (Flags 2).

## 2. Methods as they matter to a model

**The α-dicarbonyl was made and assayed in-house (Methods, p. 2458).** 6 mL of MG-1,1-dimethylacetal
into 100 mL of **2.5 % (v/v) sulfuric acid**, refluxed **1 h**; purified by **fractional distillation
under reduced pressure**, first fraction discarded for methanol. Concentration assigned by diluting
into **50 mM sodium phosphate pH 7.4** to below an estimated 2 mM, reacting with an **equal volume of
40 mM aminoguanidine** in phosphate buffer for **5–6 h at 37 °C**, and reading **absorbance at
320 nm** against a calibration curve of serial dilutions of **3-amino-1,2,4-triazine**. Fractions
then diluted to **50 mM stocks** in phosphate buffer, pH confirmed 7.4, stored at **−80 °C**. This is
the only absolute quantification in the paper and it is of the reagent, not of any product.

**Pot A — in vitro dose–response (Methods, p. 2458; Fig. 1D).**
**1 mM arginine or metformin-HCl + 1 mM N-acetylcysteine or glutathione (reduced) + MGO at 0, 0.1,
0.2, 0.5, 1, 2, 5 or 10 mM**, in **PBS**, **37 °C**, **24 h**. Reactions diluted **1:2 with 0.1 %
trifluoroacetic acid in H₂O** and frozen at **−20 °C** for later analysis. The CRV2 peptide arm used
**1 mM peptide** under the same conditions.

**Pot B — in vitro time course (Methods, p. 2458; Fig. 1C).**
**1 mM arginine or metformin-HCl + 1 mM NAC or GSH + 1 mM MGO** in **PBS**, at **25 °C** — note the
temperature difference from Pot A. **Sampled once an hour by autosampler for LC−MS, starting at
0 min.** The CRV2 peptide arm again at 1 mM.

- **Buffer and pH.** PBS throughout, i.e. ≈10 mM phosphate + 137 mM NaCl + 2.7 mM KCl at
  **pH ≈ 7.4**. The paper never states the pH of its PBS or its exact composition.
- **Atmosphere.** **Not controlled and not stated.** Bench incubations in an autosampler tray and a
  37 °C incubator; no degassing, no inert gas, no oxygen statement anywhere in the paper.
- **Replication.** **n = 4 independent biological replicates**, mean with S.E.M. (Figs. 1C, 1D, 2A–C).
- **How everything was quantified — the load-bearing methods statement.** Agilent **6540 Q-TOF**
  MS/MS with 1290 UHPLC and 1260 nanoLC-Chip, **positive ion**, mass window 50−1000 *m/z*, capillary
  3.5 kV, drying gas 300 °C at 8 L/min, nebuliser 35 psi, fragmentor 150 V. Phenomenex Gemini C18
  **50 × 4.6 mm, 5 µm** at **0.4 mL/min**; A = 0.1 % TFA in H₂O, B = 0.1 % TFA in CH₃CN; gradient
  0 % B (0–2 min), 0→30 % B (2–5), 30→100 % B (5–6), 100 % B (6–7), 100→0 % B (7–8), 0 % B (8–11).
  (Peptide runs: 0→60 % B over 0–5 min, then as above.)
  **"Relative metabolite abundance was quantified by the integrated peak area for each extracted ion
  chromatogram with a mass window of ±0.1 and normalized to the most abundant peak in a given time−
  or dose−response series."**
  **There is no calibration curve, no internal standard and no response factor for any in vitro
  product. Every in vitro *y* axis is a normalised peak area on a 0–1 scale.**
- **What was NOT measured.** MGO itself in any reaction. Any disulfide (cystine, GSSG,
  NAC-disulfide). Any absolute concentration of any product. Any mass balance. Any second pH. Any
  temperature above 37 °C. Any dissolved oxygen.

**Pot C — the route test (Fig. 2C).** MGO **pre-equilibrated for 24 h** with **either NAC or
arginine**, then the opposing partner added and MICA formation monitored by LC−MS. The result is
qualitative: NAC-first gave MICA "much faster and in higher yield" than arginine-first. **No number
is printed for either.**

**Pot D — BSA proteomics (Methods, p. 2458).** **BSA 0.5 mg/mL in PBS + 0.5 mM MGO + 0.5 mM arginine,
37 °C, 24 h.** Then 3× dialysis into fresh PBS on 30 kDa Amicon filters, 1 mM MgCl₂, trypsin at
**1:100 trypsin/protein overnight at 37 °C**, C18 desalting, lyophilisation. LC−MS/MS on an
Easy-nLC 1000 + **Q Exactive HF** orbitrap, PepMap RSLC C18 75 µm × 15 cm, 45 °C, 0.3 µL/min; full MS
at 120 000 resolution, 375−1500 *m/z*, AGC 1e6, max IT 60 ms; top-10 data-dependent HCD at 30 000
resolution, AGC 1e5, NCE 27, isolation window 2.0 *m/z*, dynamic exclusion 20 s. Searched with
**ProLuCID / IP2** against a concatenated target/decoy UniProt BSA database; up to **two differential
modification sites per peptide** from **oxidised methionine +15.9949 (M), MICA +210.1116 (C),
MG-H1 +54.0106 (R), CEL +72.0211 (K), CEA +72.0211 (R)**; delta mass cutoff 10 ppm; **FDR 1 % at the
peptide level**; reported peptides required in **at least two different experiments**.
**Note that no disulfide or oxidised-cysteine modification was in the search list** other than
methionine oxidation — so this experiment could not have found a cystine even if one had formed
(Flags 3).

**Pot E — cells.** HeLa (and HEK293T for virus), RPMI 1640 + 10 % FBS + 1 % pen/strep. Two million
cells per 10 cm plate, 24 h, then **8 h** with MGO (up to **0.5 mM**), metformin, or both, in 5 mL
medium. Metabolome extracted into **300 µL of cold 80:20 MeOH/H₂O with 1 µL of 10 mM d3-serine** as
internal standard; extracellular from 200 µL medium into 800 µL cold MeOH with 3 µL of 10 mM
d3-serine. Targeted MRM on an **Agilent 6460 QQQ**, positive ion, capillary 4.0 kV, drying gas 300 °C
at 5 L/min, nebuliser 45 psi, delta EMV(+) 200; **MRM parameters in extended data Table S1 (off
disk)**; same Gemini C18 column; **peak areas normalised to the internal standard** — so the *cell*
data are internal-standard-normalised while the *in vitro* data are series-normalised.
**MRM standards** were made by incubating **10 mM arginine or metformin-HCl + 10 mM cysteine or
glutathione + 10 mM MGO in PBS at 37 °C for 24 h** and optimising transitions with Agilent MassHunter
Optimizer.

**Pot F — the transporter.** ABCC1 knocked down in HeLa by shRNA (pLKO.1 puro; scramble = SHC002);
treated with **0.5 mM MGO for 8 h**; MRP1 and PGK1 by Western (anti-MRP1 1:1000 CST #72202, anti-PGK1
1:3000 sc-130335). Statistics: one-way ANOVA (Figs. 3E–H, 4A) or unpaired Student's *t* (Figs. 4D–G).

## 3. Tables re-typed

**THE MAIN TEXT PRINTS NO TABLES.** The only two tables are **extended data Table S1** (MRM
parameters) and **extended data Table S2** (cloning primers), both in the Supporting Information and
**both off disk**. Nothing quantitative is lost from the chemistry by their absence — S1 is
instrument settings and S2 is oligonucleotides — but it means **no MRM transition used in this paper
is recoverable from what is on disk**.

### Every number printed in the main text

| quantity | value | class | where |
|---|---|---|---|
| **any rate constant, reaction order or activation energy** | **— none printed, for any reaction, anywhere in the paper** | — | — |
| **MICA modification mass shift, on cysteine** | **+210.1116** | `[M]` (search parameter) | Methods, BSA proteomics, p. 2458 |
| **MG-H1 mass shift, on arginine** | **+54.0106** | `[M]` | Methods, p. 2458 |
| **CEL mass shift, on lysine** | **+72.0211** | `[M]` | Methods, p. 2458 |
| **CEA mass shift, on arginine** | **+72.0211** | `[M]` | Methods, p. 2458 |
| oxidised methionine mass shift | +15.9949 | `[M]` | Methods, p. 2458 |
| time at which the hemithioacetal begins to fall, 1 mM each, **25 °C** | **2–6 h** | `[M]` | Results, p. 2455 |
| thiol vs guanidine consumption | thiols consumed **more rapidly** than guanidines; "thiols are more reactive toward MGO" | `[M]`, no number | Results, p. 2455 |
| product order of appearance | **hemithioacetal first**, then MICA and imidazolone | `[M]`, no number | Results, p. 2455 |
| arginine vs metformin | **arginine more reactive toward MGO than metformin** | `[M]`, no number | Results, p. 2456 |
| equilibrium behaviour, 24 h at 37 °C | imidazolone and MICA reach maximal formation at **lower MGO equivalents** than hemithioacetal; hemithioacetal accumulates significantly **only at multiple equivalents** of MGO | `[M]`, no number | Results, p. 2456 |
| intramolecular vs intermolecular MICA | in the CRV2 peptide, MICA forms **more rapidly and at lower MGO equivalents** than in any intermolecular pair | `[M]`, no number | Results, p. 2456; Fig. 2A,B |
| route test | NAC-preincubated MGO gave MICA "**much faster and in higher yield**" than arginine-preincubated | `[M]`, no number | Results, p. 2456; Fig. 2C |
| BSA MICA site | **a single MICA site, at Cys58** — "the only nondisulfide cysteine in BSA" | `[M]` | Results, p. 2456 |
| BSA, MICA on surface arginines | **none identified** | `[M]` | Results, p. 2456 |
| MGO dose range in cells | up to **0.5 mM**, 8 h | `[M]` | Results, p. 2456 |
| effect of MGO in HeLa | significantly reduced free **arginine** but **not glutathione**; dose-dependent rise in MG-H1-arginine and GSH-Arg-MICA | `[M]`, figure-only values | Results, p. 2456; Fig. 3E–H |
| MICA metabolites detected in cells | **GSH-Arg-MICA** (± metformin) and **GSH-Met-MICA** (metformin only). **Cys-Arg-MICA and Cys-Met-MICA: not detected** | `[M]` | Results, p. 2456; Fig. 3A–D |
| free glutathione in the medium | **not detectable** in control or MGO-treated cells | `[M]` | Results, p. 2457 |
| ABCC1 knockdown | significantly **higher intracellular** and **lower extracellular** GSH-Arg-MICA; intracellular GSH and arginine **unchanged** | `[M]`, figure-only values | Results, p. 2457; Fig. 4D–G |
| MGO concentrations in cells, literature | "low-to-mid **micromolar**" | `[C]` (refs 26,27) | Results, p. 2456 |
| significance thresholds | * p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001 | — | Figs. 3, 4 captions |

**Figure-only** (per house rule, not typed): every point of Fig. 1C, Fig. 1D, Fig. 2A, **Fig. 2B (the
logistic fits — the parameters are not printed either)**, Fig. 2C, Figs. 3E–H, Fig. 4A and Figs.
4D–G; the chromatograms of Figs. 3A–D and 4B; the Western blot of Fig. 4C; the mass spectrum of
Fig. 2E; and all of Figs. S1–S4, which are off disk in any case.

### Arithmetic on the printed numbers (all mine)

**1. The MICA mass shift proves the branch — no redox step.** Using monoisotopic masses:
MGO C₃H₄O₂ = 72.02113; arginine C₆H₁₄N₄O₂ = 174.11168; H₂O = 18.010565.

| product | printed shift | my composition | my value | agreement |
|---|---:|---|---:|---|
| **MICA on Cys** | **+210.1116** | MGO + Arg − **2** H₂O | **210.1117** | **1e-4 Da** |
| **MG-H1 on Arg** | **+54.0106** | MGO − **1** H₂O | **54.0106** | exact |
| CEL on Lys / CEA on Arg | +72.0211 | MGO, **no water lost** | 72.0211 | exact |

**MICA and MG-H1 are condensations: two waters and one water respectively, and no change in
oxidation state of anything.** The thiol sulfur in a MICA crosslink is bonded to an sp² imidazole
carbon. **This is adduct-branch chemistry by mass balance, not by assertion.** (The CEL/CEA row's
net MGO addition with no water lost implies an internal hydride shift; that reading is mine, the
paper does not discuss it, and nothing in this dossier rests on it.)

**2. The two in vitro pots run at different temperatures and this is easy to miss.** The **time
courses are at 25 °C** (Methods, "In Vitro Time-Course MICA Formation Experiments … in PBS at
25 °C"; Fig. 1C caption confirms) and the **dose–responses are at 37 °C** (Fig. 1D caption and
Methods). Fig. 2B fits the peptide time course **together with** "time-course experiments in
Figure 1D" — but Fig. 1D is the **dose–response** at 37 °C, not a time course. **Either the Fig. 2B
caption mislabels its own source or the two temperatures are pooled in one fit** (Flags 4). Since no
fitted parameter is printed, nothing downstream depends on the answer, but a reader should not treat
Fig. 2B as a single-temperature result.

**3. There is no ratio to be taken.** Because every in vitro series is normalised to *its own* most
abundant peak, **the hemithioacetal curve and the MICA curve in the same panel are on different
implicit scales**, and no MS response factor is given for any species. So the corpus's usual
`within_study_ratio` extraction is **unavailable for this paper**: one cannot say "MICA reached x %
of the hemithioacetal", only "MICA rose while the hemithioacetal fell". That is the single largest
practical limitation of this dataset.

## 4. Numbers the repository can use

| quantity | value | unit | conditions | anchor | evidence class |
|---|---|---|---|---|---|
| **rate constant for thiol + MGO (either branch)** | **— none printed** | — | — | — | **absent** |
| **reaction order for anything** | **— none printed** | — | — | — | **absent** |
| **activation energy for anything** | **— none printed** (single-temperature per experiment; two temperatures across experiments but never the same reaction at both) | — | — | — | **absent** |
| **rate, order or barrier for thiol → disulfide** | **— none, and no disulfide is measured, sought or discussed anywhere** | — | — | — | **absent** |
| time at which the hemithioacetal turns over | **2–6 h** | h | 1 mM GSH or NAC + 1 mM Arg or metformin + 1 mM MGO, PBS ≈pH 7.4, **25 °C** | Results, p. 2455; Fig. 1C | **`measured_bound`** — a turning point, not a rate; the only time-like number in the paper |
| kinetic ordering, reactant consumption | **thiol > guanidine**, and **arginine > metformin** | — | as above | Results, pp. 2455–2456 | **`level_only`** (an ordering; no magnitude is recoverable — see §3 item 3) |
| kinetic ordering, product appearance | **hemithioacetal → then MICA and imidazolone** | — | as above | Results, p. 2455 | `level_only` |
| thermodynamic ordering | imidazolone and MICA maximal at **lower MGO equivalents** than hemithioacetal | — | 1 mM each + 0.1–10 mM MGO, PBS, **37 °C, 24 h** | Results, p. 2456; Fig. 1D | `level_only` |
| MICA is a double condensation, not a redox step | **+210.1116 = MGO + Arg − 2 H₂O** | Da | proteomic search parameter | Methods p. 2458; my arithmetic reproduces it to 1e-4 Da | **`measured_ratio`** (a mass identity, and the paper's most transferable single fact) |
| MG-H1 is a single condensation | **+54.0106 = MGO − H₂O** | Da | as above | Methods, p. 2458 (mine, exact) | `measured_ratio` |
| MICA site occupancy on a real protein | **exactly one site, Cys58**, the only non-disulfide cysteine in BSA; **zero** MICA on surface arginines | sites | BSA 0.5 mg/mL + 0.5 mM MGO + 0.5 mM Arg, PBS, 37 °C, 24 h, 1 % FDR, ≥2 experiments | Results, p. 2456; Fig. 2D,E | **`measured_bound`** — a bound on how promiscuous the crosslink is on a protein surface |
| MICA requires a guanidine partner | structurally, yes — cysteine + MGO alone gives only the reversible hemithioacetal | — | throughout | Results, p. 2456; Fig. 1A | **`derived_assumption`** (mine, from the mechanism and the mass balance) — **relevant because `whitfield_nf_cys_MFT` has no arginine** |
| proximity effect | intramolecular MICA (CRV2 peptide, Cys and Arg in one chain) forms **faster and at lower MGO equivalents** than any intermolecular pair | — | 1 mM peptide + 1 mM MGO, PBS, 25 °C (time course) / 37 °C 24 h (dose) | Results, p. 2456; Fig. 2A,B | `level_only` |
| route to MICA | via the **hemithioacetal first**, not via a preformed imidazolone | — | 24 h pre-equilibration test, PBS | Results, p. 2456; Fig. 2C | `level_only` — a mechanistic ordering with no number |
| any absolute concentration of any in vitro product | **— none** (every in vitro axis is a peak area normalised to the most abundant peak in its own series) | — | — | Methods, p. 2458 | **absent** — see Flags 2 |

### Adduct branch or redox branch?

**Every reaction in this paper is ADDUCT branch, and for once the assignment is arithmetic rather
than inferential.**

| species | branch | how it was decided |
|---|---|---|
| **hemithioacetal (thiol + MGO)** | **ADDUCT, reversible** | The paper's own framing throughout ("reversible hemithioacetal formation on cysteine thiolates"; "the hemithioacetal modification is reversible"). Nucleophilic thiolate addition to the aldehyde carbon. C–S bond. |
| **MICA (mercaptomethylimidazole crosslink)** | **ADDUCT, irreversible** | **Mass identity: +210.1116 on cysteine = MGO + arginine − 2 H₂O (mine, 210.1117).** Two dehydrations and no change of oxidation state. The sulfur ends bonded to an imidazole ring carbon. Described as "stable" and "thermodynamically favored" throughout. |
| **MG-H1 imidazolone** | **ADDUCT** | +54.0106 = MGO − H₂O (mine, exact). A guanidine adduct; **no sulfur involved at all**. |
| CEL / CEA | **ADDUCT** | +72.0211, a net MGO addition on a lysine or arginine. No sulfur. |
| **any disulfide** | **— NOT PRESENT** | The word "disulfide" occurs in this paper exactly twice, both times as the phrase "the only **nondisulfide** cysteine in BSA" — i.e. as a way of saying which cysteine was free, not as a product. **No cystine, no GSSG, no NAC-disulfide is measured, searched for, or mentioned.** The proteomic search list contains oxidised methionine and four MGO adducts, and **no cysteine oxidation state at all** (Flags 3). |

**The one thing this paper does say that bears on the redox branch, indirectly.** In HeLa cells at up
to 0.5 mM MGO for 8 h, **free glutathione was not significantly reduced** while free arginine was
(Results, p. 2456; Fig. 3E,F). If MGO were an efficient thiol oxidant at 37 °C, a 0.5 mM dose against
a millimolar intracellular GSH pool over 8 hours would be expected to show *something* on the GSH
axis. It did not. **That is a cell-based, confounded, figure-only observation** — cells regenerate
GSH enzymatically, so the null is not clean — but it points the same way as Zheng 2022's "limited"
GSSG, at a third laboratory.

## 5. Flags

1. **TEMPERATURE TRANSFER, AND HERE IT IS ALMOST ACADEMIC BECAUSE THERE IS NOTHING TO TRANSFER.**
   The two in vitro pots are at **25 °C** (time courses) and **37 °C** (dose–responses); the wave
   concerns **140 °C**, gaps of **115 K** and **103 K**. But the deeper problem is that **this paper
   has no constant to carry**, so the usual Arrhenius accounting does not even apply. What *is*
   transferable — the branch identities, the ordering, the mass balances — is structural and
   temperature-robust in kind, though not in degree: at 140 °C a hemithioacetal with a sub-second
   lifetime may never accumulate at all, and the "2–6 h turnover" observation would collapse to
   something unobservable. **Nothing in this paper licenses any statement about a rate at Maillard
   temperature, and the paper contains two temperatures that are never applied to the same reaction,
   so not even a crude two-point barrier can be extracted** (Flags 4).
2. **EVERY IN VITRO NUMBER IS A RELATIVE PEAK AREA NORMALISED WITHIN ITS OWN SERIES.** No calibration
   curve, no internal standard, no response factor for any in vitro product. This has two hard
   consequences: (a) **no absolute concentration exists**, so no rate constant could ever be
   recovered from these data even by re-fitting them; (b) **no cross-species ratio exists**, because
   two curves in the same panel are each normalised to their own series maximum and the MS response
   of a hemithioacetal and of a mercaptomethylimidazole are not the same. The corpus's
   `within_study_ratio` class is simply unavailable here. The authors concede the point in the
   Discussion.
3. **The BSA proteomic search could not have detected a disulfide.** The differential modification
   list is oxidised methionine, MICA, MG-H1, CEL and CEA — **no cysteine oxidation, no cystine, no
   sulfinic/sulfonic acid, and no disulfide crosslink search**. Combined with **"two total
   differential modification sites per peptide"** and **trypsin with three missed cleavages**, this
   is a search designed to find MGO adducts and structurally blind to thiol redox chemistry.
   **The absence of a disulfide finding in this paper is not evidence of a disulfide's absence.**
4. **Figure 2B's caption does not match the Methods.** It reads "Logistic fit of MICA formation for
   peptide time-course experiment **and time-course experiments in Figure 1D**" — but Figure 1D is
   the **24 h dose–response at 37 °C**, not a time course, and the actual time courses (Figure 1C)
   are at **25 °C**. Either the reference is to the wrong panel or two temperatures are pooled into
   one fit. **No fitted parameter is printed for Fig. 2B**, so nothing quantitative rests on it, but
   the panel should not be cited as a single-condition kinetic result.
5. **Only one α-dicarbonyl, and it is again an α-oxoaldehyde.** MGO's aldehyde carbon is what both
   the hemithioacetal and (therefore) MICA require. **2,3-pentanedione, 2,4-pentanedione and
   3,4-hexanedione have no aldehyde carbon**, so neither product class in this paper is available to
   the α-diketones the pre-registration's §3(a) source step would make. Same structural caveat as in
   `zheng2022_extraction.md` Flags 4 and `zheng2023_extraction.md` Flags 8; it applies to all three.
6. **MICA needs a guanidine, and the pot that must show the disulfide share has none.**
   `fed_nf_cys_MFT` and `whitfield_nf_cys_MFT` carry cysteine, norfuraneol and (in the H₂S pots)
   hydrogen sulfide — no arginine, no protein, no biguanide. **This mechanism is unavailable there.**
   It *is* available in any fed system carrying protein, which is most of them, so if the MICA route
   is ever added it must be gated on a guanidine pool and not charged globally.
7. **The thiols are again aliphatic (GSH, NAC, cysteine, a peptide cysteine, BSA Cys58).** No
   heteroaromatic thiol anywhere. See `zheng2022_extraction.md` Flags 9.
8. **PBS is not characterised.** The paper says "PBS" and never states its composition, ionic
   strength or pH. Standard PBS is ≈10 mM phosphate at pH 7.4 with 137 mM NaCl, but that is an
   assumption, and the chloride is a difference from every other pot in this batch.
9. **Atmosphere is uncontrolled and unmentioned throughout** — no degassing, no inert gas, no
   oxygen statement, in a paper about thiol chemistry over 24 h. Unlike the Zheng papers, this one
   does not even remark on autoxidation.
10. **The MGO is made in-house and quantified by a colorimetric proxy.** Acid hydrolysis of the
    dimethylacetal in 2.5 % H₂SO₄, fractional distillation, then the aminoguanidine/320 nm assay
    against a **3-amino-1,2,4-triazine** calibrant — i.e. the standard is the *product* of the
    derivatisation, not MGO. Sound and conventional, but it means the stated 1 mM is an assay-derived
    figure with unstated uncertainty, and no purity or water content is reported for the distillate.
11. **All Supporting Information is off disk.** Figures S1 (structures), **S2 (the 24 h screen
    showing that all four thiol/guanidine pairs form imidazolone, hemithioacetal and MICA — the
    paper's breadth claim)**, S3 (synthetic standards and fragmentation), S4 (shABCC1 metformin data),
    and extended data Tables S1 (MRM parameters) and S2 (primers). **Nothing in the main text
    recovers any of them.**
12. **What to request from the authors:** (i) **absolute concentrations, or at minimum MS response
    factors, for the hemithioacetal, MICA and MG-H1 species**, without which no constant can ever be
    fitted to Fig. 1C; (ii) the **logistic fit parameters** behind Fig. 2B and confirmation of which
    experiments it pools; (iii) whether **any disulfide (GSSG, NAC-disulfide, BSA Cys58 oxidation)**
    was ever detected in the LC−MS traces and simply not searched or not reported; (iv) any run at a
    temperature above 37 °C; (v) the time course repeated at 37 °C so that a two-point barrier for
    hemithioacetal turnover could be estimated.
13. **What this paper does NOT contain:** any main-text table; any rate constant; any reaction order;
    any activation energy; any absolute in vitro concentration; any calibration curve for a product;
    any disulfide; any cysteine-oxidation search; any temperature above 37 °C; any pH other than
    PBS's; any α-diketone; any heteroaromatic thiol; any food matrix; any volatile measurement; and
    **any DFT, quantum-chemical, molecular-dynamics or docking calculation — the PDB 4F5S in
    Figure 2D is a rendering of observed modification sites, not a computation, so nothing was
    excluded under the no-DFT policy.**
