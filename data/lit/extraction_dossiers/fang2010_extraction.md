# Fang 2010 — EXTRACTION (the glucose half of the same experiment as Fang 2009: equimolar glucose + glycine co-precipitated and heated dry at 125 C for 2 h, with glucose labelled one carbon at a time — 13C1, 13C2, 13C3, 13C6 and U-13C6 — against 15N-glycine, and every carbon's fate resolved by quantitative solid-state NMR with spectral editing)

### THE SUGAR-SIDE COMPANION: it says what happens to the six carbons the trunk books into `MEL_C` — **about 40 % of melanoidin carbon is aromatic or alkene and ~51 % is still alkyl**, C6 barely changes while C1 makes new C-C bonds, **between one-quarter and one-half of all glucose carbon backbones are fragmented**, and no single structure accounts for more than 15 % of the polymer.

**Source on disk:** `data/articles/fang2010.pdf` (10 pp., J. Agric. Food Chem. **2011**, 59 (2),
481-490 — the file name says 2010 because the paper was published on the web 28 December 2010;
see Flags 1). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/fang2010.txt`). **The paper has exactly one table, Table 1, and it contains
no numbers at all** — it is a grid of nine drawn chemical structures lettered A, B, C, D, Ea, Eb,
F, G, H. Its caption came through the text layer; the structures did not, so page 488 was
**rendered at 190 dpi (`pdftoppm -r 190 -f 8 -l 8`)** and the structures are described in words in
section 3. Every percentage in this paper is printed in the running text, not in a table.
Figures 1-10 are pulse sequences and NMR spectra and are **figure_only**; **Figure 4 is a bar
chart of per-carbon functional-group percentages with numeric labels and is treated as
figure_only per house rule** — the values that also appear in the running text are carried, the
rest are not. **Supporting Information exists and is NOT on disk**, and in this case it holds the
paper's main quantitative table — see Flags 2.
Repo status before this dossier: `fang2010.pdf` has **no extraction dossier** and is not cited
anywhere in `src/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Alkyl and Other Major Structures in 13C-Labeled Glucose-Glycine Melanoidins Identified by Solid-State Nuclear Magnetic Resonance" |
| Authors | Xiaowen Fang and Klaus Schmidt-Rohr (corresponding), Department of Chemistry, Iowa State University, Ames, Iowa 50011 |
| Venue | J. Agric. Food Chem. **2011**, 59 (2), 481-490. © 2010 American Chemical Society, **published on Web 28 December 2010** |
| DOI | 10.1021/jf102917v |
| Naming | "melanoidin" = the **soluble HMW dialysis retentate**, MWCO 6000-8000, unless the insoluble fraction is named explicitly. "C1"-"C6" are **glucose** carbon positions throughout (in the companion paper the same labels mean glycine carbons — see Flags 3). "Monomer unit" = a structural fragment identified from 2D correlations, lettered A-H in Table 1 |
| Lineage | same COST Action 919 dry-reaction protocol as the companion; explicitly the sugar-side sequel — "We have previously studied the fate of glycine in this material (20) and are now focusing on the more complex transformations of glucose", where ref. 20 **is Fang & Schmidt-Rohr 2009**. Argues against three named melanoidin models: **Cämmerer & Kroh's** (ref. 12, no pyrroles, no C2 in alkyl sites), the **enamine/imine models** (refs. 13, 14), and **Tressl's** polycondensed-pyrrole/furan model (ref. 4) |
| Companion on disk | **`fang2009_extraction.md` is the same experiment with the labels on the amino acid.** The two papers share one sample preparation and must be read as one study |
| Companions on disk | `mundt2004_extraction.md` (the C/N measurement on the same chemistry), `cammerer1994_extraction.md` (this paper's ref. 12, whose model it partly refutes), `hofmann2000_extraction.md` and relatives (Hofmann is ref. 9 and, with Ames, ref. 11 — the source of the "insoluble fraction is the majority" claim), `adams2008_extraction.md` (refs. 22-23, the Ghent thermal-degradation work on melanoidins from this same protocol) |

## 1. Why it matters

The trunk routes carbon into `MEL_C` through one step: Martins' step 9, "3-DG + Gly ->
melanoidins", which contributes **six intact glucose-derived carbons** per event
(`MELANOIDIN_REPEAT_UNIT_CARBON = 8` in `src/kinetic_core/species.py`, of which 6 are the sugar's
and 2 the amine's). Fang 2009 gave the amine side of that unit. **This paper gives the sugar
side**, and it makes three statements that bear on the trunk directly:

**(a) The six carbons do not travel together.** "Overall, we estimate that between 1/4 and 1/2 of
all glucose carbon backbones are fragmented." The intact-C6 assumption in the repeat unit holds
for at most three-quarters of events and possibly only half. Every carbon that leaves the
backbone either ends up in the polymer anyway (as a separate fragment) or does not — and the
trunk's `FRAG_C` pool exists for exactly the second case. `species.py` calls `FRAG_C` "carbon
leaving a measured step in an unmeasured co-product" and says its size should be "visible rather
than hidden". This paper is the corpus's best evidence that the pool is chemically real and not
small.

**(b) The polymer is not aromatic, which contradicts how melanoidins are usually pictured.**
"Aromatic and alkene carbons make up only 40 % of the total in the melanoidin"; "alkyl carbons
account for ~51 % of sugar carbons, and more than 85 % of them are protonated". The trunk does
not model structure and does not need to, but this is the number that explains why an elemental
C/N near 7.6 and an H/C near 1.5 (Mundt 2004, arithmetic 3 in that dossier) are consistent with
each other: the polymer is half unreacted-looking sugar chain. Any future attempt to give
`MEL_C` a molecular interpretation has to start here.

**(c) It contradicts a structure the B7 furanic block might be assumed to feed.**
`src/kinetic_core/species.py` carries HMF on the trunk (B7 block) and `data/keys/compounds.yml`
keys `hmf` and `furfural`. Hodge's scheme has 5-hydroxymethyl-2-furaldehyde as the melanoidin
precursor with C1 as the aldehyde and C6 intact. Fang looks for it and does not find it: "**no
aldehyde bonded to a furan ring is observed in our sample**" (while an aldehyde bonded to a
*pyrrole* is clearly detected). His conclusion is that C1 is split off — C1 + C5 fragmentation —
before the furan enters the polymer. **So HMF, if it is a melanoidin precursor here, does not
enter the polymer as HMF.** Furan rings and their associated C6H2 are ~11 % of sugar carbons;
pyrroles ~8 %; imidazolium plus oxazolium ~6 %.

Two further findings worth carrying:

- **The dry reaction needs the amine.** "Without glycine, dry heating of glucose at 125 C for
  2 h does not lead to major structural changes" — the neat-glucose spectrum keeps its three
  original peaks, "only broadened". At this temperature and time, **caramelisation alone produces
  nothing that looks like a melanoidin**. That is a clean null for any lane that competes
  caramelisation against the Maillard route in a low-moisture matrix.
- **The insoluble fraction is finally characterised.** Fang 2009 covered the ~50 % insoluble
  product with the single word "similar". Here it is measured: "similar spectral features and
  percentages as in the HMW fraction ... with **larger aromatic and smaller O-alkyl** components
  in the insoluble fraction as the only significant differences", the two differing "mostly by a
  higher degree of polymerization or cross-linking, not by further major chemical
  transformations". That substantially softens Flag 2 of the 2009 dossier — though the
  supporting figure is in the missing SI.

What this paper does NOT give the repository: any rate constant, any activation energy, any time
course, any temperature series, any elemental analysis, any C/N, and **any measurement on the
amino acid** (that is the 2009 paper).

## 2. Methods as they matter to a model

- **One pot, one condition.** "An equimolar **coprecipitated** mixture of glucose and glycine was
  heated for **2 h at 125 C**", following the COST Action 919 protocol. Note the word
  *coprecipitated* — Fang 2009 describes the same step as dissolving in water and freeze-drying,
  which is what produces the co-precipitate. **This is the same preparation as the 2009 dry
  reaction**, and the paper says so ("More details are given in ref 20"). Low moisture; no water
  added.
- **Reactants.** D-glucose anhydrous 99+ % (MW **180.16** printed), glycine 98 % (FW **75.06**
  printed), Acros Organics. Labelled glucose at 13C1, 13C2, 13C3, 13C6 and U-13C6 (all 99 %) from
  Cambridge Isotope Laboratories; glucose-1,2-13C2 (99 %) from Isotec; glycine-15N (98 %).
- **Fractionation.** Whatman 41 ashless filter paper separates the insoluble from the soluble
  fraction; the soluble part is dialysed (Fisherbrand regenerated cellulose, **MWCO 6000-8000**,
  5.10 mL/cm) and the retentate is the HMW melanoidin. **Both** the soluble HMW and the insoluble
  fraction were measured here; the paper then follows COST Action 919 in reporting the soluble
  HMW fraction "to enable comparison of our results with those of other research groups".
- **Ratio series.** Melanoidins were also made at **9:1 and 1:0 glucose:glycine** (Figure 2), not
  only 1:1. The 1:0 sample is the neat-glucose control described above. The 9:1 comparison is
  qualitative only (Flags 6).
- **The quantitative measurement.** 13C **direct-polarisation** with a Hahn echo at 14 kHz MAS on
  a Bruker DSX400 (100 MHz 13C, 40 MHz 15N), 4 mm triple-resonance probe; recycle delays set from
  measured T1 so residual signal is < 5 %, i.e. **all carbon sites fully relaxed**. Recycle delays
  **100-220 s**, 32-128 scans. Only these DP spectra are quantitative; the CP/TOSS spectra (6.5
  kHz, 1 ms contact, 3 s recycle) are for routine characterisation, and the 13C{15N} REDOR
  results are called **semiquantitative** by the paper itself.
- **The natural-abundance correction, which is larger here than in the companion.** 13C at
  natural abundance is 1.1 % of all carbon, but in a sample where only **one** of six glucose
  carbons is labelled that background is **8 % of the spectral intensity**. It was subtracted,
  using a synthetic background built as the sum of the U-13C6-glucose and U-13C-glycine DP spectra
  scaled by 0.011 — because "long T1C relaxation times made it impractical" to measure it directly
  and CP under-represents nonprotonated aromatic carbons. **Every single-label percentage in this
  paper depends on that constructed background.**
- **The C4 + C5 spectrum is a difference, not a measurement.** No 13C4 or 13C5 glucose was bought.
  "The spectrum of 13C4 plus 13C5 was obtained by **subtracting** the sum of the spectra of
  glucose-13C1, -13C2, -13C3 and -13C6 from the spectrum of glucose-13C6." So every C4/C5 number
  carries the accumulated error of five spectra (Flags 4).
- **Spectral editing techniques, and what each one separates.** CH-only by dipolar DEPT; CH2-only
  by three-spin coherence selection at 5.787 kHz; **alkyl (sp3) selection by a 13C CSA filter**
  with a 38 us filter time, which is what resolves the 115-100 ppm overlap between alkyl and
  aromatic carbons; nonprotonated + CH3 by 40 or 68 us of gated decoupling; carbons near nitrogen
  by 13C{15N} REDOR (5 kHz, Ntr = 1.6 ms) and by a **MELODI** sequence (7 kHz, Ntr = 5.6 ms) that
  separates C-N from C--C-N. Connectivities by 2D 15N-13C HSQC (7 kHz, ~10 h per spectrum) and 2D
  13C-13C spin exchange (7 kHz, 50 ms and 10 ms mixing, ~8 h per spectrum). Spatial mixing by
  1H-13C HetCor with Lee-Goldburg CP and 0.45 ms of 1H spin diffusion.
- **What is not measured.** No time course, no temperature series, no water content, no molecular
  weight distribution, no yield, no elemental analysis, no absorbance, no small-molecule
  quantitation, and no pH (there is no solvent).

## 3. Tables re-typed

### Table 1. "Alkyl and Other Major Structural Fragments Identified in Melanoidins Made from Glucose Reacted with Glycine in a Molar 1:1 Ratio in a Dry Reaction"

**This table contains no numbers.** It is a 3 x 3 grid of drawn structures, lettered to match the
cross-peak labels in Figure 9. Read from the page render and described here; the abundances,
which are printed in the running text and not in the table, are attached in brackets.

| label | structure as drawn (glucose carbon positions numbered 1-6 where the paper numbers them) | abundance printed in the text |
|---|---|---|
| **A** | R−C(=O)−N(R)−**C1**H2− joined to a saturated ring: **C2** as a quaternary O−C−O (ketal) bearing HO, then **C3**(OH)−**C4**(OH)−**C5**−O closing the ring, with **C6**H2−R exocyclic | **~10 % of all sugar C** |
| **B** | R−O/N−**C1**(=O)−**C2**H(−N(R/H)−CH2−COO−)−**C3**H2−**C4**H(OH)−**C5**H(OH)−**C6**H2−R — an open chain with the glycine still attached as N−CH2−COO− | **~8 % of all sugar C** |
| **C** | H**C1**(=O)− bonded to a **pyrrole** ring carrying **C6**H2−R, the ring N substituted with −CH2−COO− (i.e. the glycine) | **~1.5 % of all C** |
| **D** | H3**C1**−**C2**(=O)−O−R — an acetyl ester | **~1 % of all C** |
| **Ea** | R−C(=O)−N(R)−**C1**(R)(O−R)−**C2**(=O)− bonded to a **pyrrole** ring numbered 3,4,5,6, ring N substituted with R | **Ea + Eb together ~4 % of all sugar C** |
| **Eb** | the open-chain analogue: R−C(=O)−N(R)−**C1**(R)(O−R)−**C2**(=O)−**C3**H2−**C4**H(OH)−**C5**H(OH)−**C6**H2−R | (included in the 4 % above) |
| **F** | R−**C1**H(R)− bonded to a **furan** ring (positions 3,4) with **C6**H2−O−R on the far side; the C1 label is drawn parenthesised "(1)" | **furans and their associated C6H2 ~11 % of total sugar C** |
| **G** | R−**C1**H(R)− bonded to a **pyrrole** ring with **C6**H2−R/H, ring N substituted with −CH2−COO−; C1 again drawn "(1)" | **pyrroles including associated alkyl-C6 ~8 % of all glucose C** |
| **H** | R−**C2**H(OH)−**C3**H=**C4**H−**C5**H((O)H)−**C6**H2−R — a non-aromatic **alkene** | **~6 % of all glucose C** |
| (not lettered) | imidazolium and oxazolium rings, with associated alkyl groups | **~6 %** |

Nine lettered structures plus imidazolium and oxazolium make the **11 "monomer units"** the
abstract counts. **My sum of the printed abundances: 10 + 8 + 1.5 + 1 + 4 + 11 + 8 + 6 + 6 =
55.5 %**, against the paper's own "about half of all the glucose carbon" — consistent.

### Every quantity printed in the running text and the abstract

Grouped by what it is about. All refer to the **soluble HMW fraction, glucose:glycine 1:1, dry,
125 C, 2 h** unless a row says otherwise.

**Bulk composition of the melanoidin (the numbers the repository is most likely to want):**

| quantity | value | where |
|---|---|---|
| **aromatic + alkene carbons** | **40 % of all C** (the paper writes "≤ 40 %" once and "about 40 %" twice) | Abstract; Results, "sp2 Hybridized Carbons"; Figure 4 discussion |
| **all sp2-hybridised carbons** | **49 % of total sugar carbons** | Results, "sp2 Hybridized Carbons" |
| **alkyl (sp3) carbons** | **~51 % of sugar carbons**, of which **> 85 % are protonated** | Results, "Alkyl Fractions" |
| carbonyl (C=O) carbons | **~9 %** overall; **4-17 %** depending on which glucose carbon | Results; Figure 4 |
| short alkyl segments associated with the aromatic fraction | 10-20 % | Results, "Alkyl Fractions" |
| longer alkyl segments and their associated C=O | 40-50 % | Results, "Alkyl Fractions" |
| carbon resonating 140-115 ppm | **19 % of all C**, of which **2/3** is aromatic C bonded to N → **~12 % of all carbons are aromatic and N-bonded** | Results |
| aromatic C two bonds from N or O (116.5-104 ppm), after di-O-alkyl correction | **~19 %**, = total pyrrole + furan including one associated alkyl C | Results, "Pyrroles" |

**Per-glucose-carbon fates:**

| quantity | value | where |
|---|---|---|
| alkyl share of C3 and of C4+C5 | **~52 %** each | Results, "Alkyl Fractions" (from Figure 4) |
| alkyl share of C6 | **82 %** | Results, "Alkyl Fractions" |
| C6 remaining as OCH2 | **50 %** | Results, "Synopsis" |
| C4 and C5 remaining as OCH | **25 %** each | Results, "Synopsis" |
| nonprotonated C1 between 150 and 20 ppm | **29 % of total C1** | Results, "New C-C Bond Formation" |
| alkyl-CH C1 | **20 % of total C1** | same |
| C−CH2−C C1 | **5 % of total C1** | same |
| **C1 carbons forming NEW C−C bonds** | **"about half"** (abstract: "more than half") | Results; Abstract |
| C6 in C−CH2−C units | **16 %** | Results |
| C6 in aromatic rings (nonprotonated) | **5 %** | Results |
| C1 and C2 directly bonded to N, within the 140-115 ppm band | **~85 %** each | Results |
| C3 directly bonded to N, same band | **9 %** | Results |
| C6 directly bonded to N | "very little"; in the dominant OCH2 form not even within two bonds | Results |
| carbons in the 116.5-104 ppm range that are two bonds from N | ~50 % | Results |
| protonated sp2 carbons near 132 ppm, C3 and C4 | **~7 %** each (the HC3=C4H alkene) | Results, "Alkene Structure" |
| C2's share of the total ketone signal near 200 ppm | **half** | Results, "Synopsis" |
| COO / NC=O intensity, the fragmentation marker | **most intense for C2, ~7 % of C2** | Results, "Fragmentation" |

**Structural inventory (the section-3 table above carries these too):**
furans ~11 %; pyrroles ~8 %; imidazolium + oxazolium ~6 %; alkene ~6 %; structure A ~10 %;
B ~8 %; Ea+Eb ~4 %; C ~1.5 %; D ~1 %.

**Fragmentation and structural-diversity claims:**

| quantity | value | where |
|---|---|---|
| **glucose carbon backbones fragmented** | **between 1/4 and 1/2** | Results, "Fragmentation" |
| dominant fragmentation site | **C2**, i.e. C1 + C5 or C2 + C4 splitting | Results, "Fragmentation" |
| isolated NCH3 and OCH3 glucose-derived carbons | **completely absent** | Results, "Fragmentation" |
| aldehyde bonded to a **furan** ring | **not observed** (an aldehyde bonded to a **pyrrole** is clearly detected) | Results, "Furans" |
| pyridines | **insignificant**; 13CH or 13C{15N} signals near 160 ppm "observed very weakly only" | Results, "Synopsis" |
| largest single peak in the 13C1 spectrum | **< 15 % of total intensity** | Results, "Synopsis"; Abstract |
| number of distinct structures implied | **at least 10** of significant concentration (10 editing peaks of similar intensity) | Results, "Synopsis" |
| the 11 identified units together | **> half of all glucose carbon** | Abstract; Results |
| structural diversity by carbon | C2 one peak ~20 %; C3 two peaks ~20 %; C4/C5 25 % OCH; C6 50 % OCH2 | Results, "Synopsis" |

**Sample and technique facts:**

| quantity | value | where |
|---|---|---|
| natural-abundance 13C | 1.1 % of all C, contributing **8 % of spectral intensity** in single-label samples | Experimental |
| **insoluble vs soluble HMW fraction** | "similar spectral features and percentages", differing only by **larger aromatic and smaller O-alkyl** in the insoluble; interpreted as a higher degree of polymerisation or cross-linking, **not** further chemical transformation | Results, "Insoluble vs HMW Fractions" |
| insoluble fraction's share of the product | "the majority of this material", cited to Hofmann, Ames, Krome & Faist 2001 | Results, ref. 11 |
| small molecules (< 1 kDa) as a share of Maillard products | **< 25 %**, cited to the same ref. 11 | Introduction |
| **neat glucose, 125 C, 2 h, no glycine** | spectrum "basically unchanged", three original peaks "only broadened"; amorphous or oligomerised **without significant structural change** | Results, "Overall Spectra"; Figure 2c |
| 9:1 glucose:glycine vs 1:1 | more ketones (220-185 ppm), **more furans** (150 and 110 ppm), **fewer N-heterocycles** such as pyrrole (135 and 110 ppm), more O−CH remaining | Results, "Overall Spectra"; Figure 2b,e — **qualitative, no numbers given** |
| unfractionated 1:1 reaction mixture | "completely dominated by −O−CH signals"; HMW melanoidin is "apparently only a small fraction" of it, so it was not studied further | Results |
| neat glucose 13C shifts before reaction | 93 (anomeric O−CH−O, C1), 72 (the four CH−OH), 62 ppm (CH2−OH, C6) | Results |
| 1H spin-diffusion length | 0.45 ms with D ≈ 0.5 nm²/ms equilibrates a **~2 nm** sphere; aromatic and O-alkyl segments are mixed on that scale, refuting a "core-shell" model | Results, "Mixing of Alkyl and Aromatic Units" |

**All ten figures are FIGURE-ONLY.** Figure 4's bar-chart labels are not typed here; the values
that also appear in the running text are in the tables above and are sourced to the text.

### Arithmetic on the printed numbers (all mine)

**1. The composition closes.** Alkyl 51 % + sp2 49 % = 100 % of sugar carbons. Within the sp2
half, aromatic + alkene is 40 % and carbonyl ~9 %, summing to 49 %. **The paper's four
categories are internally consistent**, which is what one expects when they come from integrating
one quantitative spectrum into four windows.

**2. The structural inventory sums to 55.5 %** (section 3 table), against the claimed "about
half". Given that the individual figures are quoted to one or two significant digits with "ca."
in front of each, the agreement is as good as it can be. **No component reaches 15 %**, which is
the paper's own diversity criterion, and the largest single one is the 11 % furan block.

**3. What this does to the trunk's repeat unit (mine).** `MELANOIDIN_REPEAT_UNIT_CARBON = 8`
assumes six glucose-derived carbons arrive together from one 3-deoxyglucosone. Fang's
fragmentation estimate is 25-50 % of backbones broken. If a broken backbone contributes on
average **half** its carbon to the polymer and the other half elsewhere, the sugar term in the
repeat unit falls from 6.00 to between **5.25 and 4.50** carbons — which would push the model's
C/N **down**, in the same direction as the amine-side correction from Fang 2009 and Mundt 2004,
and by a much larger amount. **This is a scoping calculation on a range the paper states as
"between 1/4 and 1/2", not a proposed constant.** It is recorded because it shows the sugar side
of the repeat unit is at least as uncertain as the amine side, and the repository currently
treats the sugar side as exact.

**4. Cross-check against Mundt 2004's H/C (mine).** Mundt's glucose-glycine melanoidin at
MW > 12500 gives H/C = **1.50** (that dossier, arithmetic 3). A fully aromatic polymer would sit
near 0.5-1.0; a saturated sugar chain near 2.0. Fang's finding that **51 % of the carbon is still
alkyl and > 85 % of that is protonated** is exactly what an H/C of 1.5 requires. **Two
laboratories, two techniques, two decades apart, and the same picture of a half-aliphatic
polymer.** Neither paper cites the other on this point; the agreement is mine.

**5. The 9:1 comparison is the only nitrogen-content lever in the paper, and it is
unquantified (mine).** Going from 1:1 to 9:1 glucose:glycine, Fang reports more furans and fewer
pyrroles — i.e. **less nitrogen in the polymer at higher sugar:amine ratio**, which is the same
direction as Cämmerer & Kroh's C/N-versus-ratio result that this paper's own introduction says
found the elemental composition "little affected" by a 4:1 to 1:1 change. **The paper gives no
number for the 9:1 melanoidin's composition**, so nothing here can be turned into a C/N versus
ratio curve. See Flags 6.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Glucose, glycine and any melanoidin
pool are **absent**. The registry does key `hmf`, `furfural`, `furan`, `2_acetylfuran`,
`2_pentylfuran`, `2_acetyl_1_pyrroline`, `hdmf`, `hemf`, `norfuraneol` and eight pyrazines.
This paper bears on `hmf` in the negative (no furan-bound aldehyde in the polymer) and on the
pyrazines in the negative (pyridines insignificant; pyrazines were ruled out in the companion).
Nothing here has a positive registry surface.

Conditions shared by every row: **equimolar glucose + glycine, co-precipitated from water by
freeze-drying, heated dry at 125 C for 2 h in a closed vessel; soluble HMW fraction after
filtration and dialysis at MWCO 6000-8000; percentages measured by quantitative 13C
direct-polarisation MAS NMR with the natural-abundance background subtracted.** No time, no
temperature series, no solvent.

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **aromatic + alkene carbon in the melanoidin** | **40** | % of total C | as above | Abstract; Results p. 484 | **elemental_analysis** (composition by quantitative NMR, not by combustion) |
| all sp2-hybridised carbon | 49 | % of total sugar C | as above | Results p. 485 | elemental_analysis |
| **alkyl (sp3) carbon** | **~51** | % of sugar C | as above | Results p. 485 | elemental_analysis |
| protonated share of that alkyl carbon | > 85 | % | as above | Results p. 485 | within_study_ratio |
| carbonyl carbon | ~9 | % of total C | as above | Results p. 485 | elemental_analysis |
| carbon in the 140-115 ppm band | 19 | % of all C | as above | Results p. 486 | elemental_analysis |
| aromatic carbon bonded to N | ~12 | % of all C | as above (19 % × 2/3) | Results p. 486 (the paper's own arithmetic) | within_study_ratio |
| pyrrole + furan carbon including one associated alkyl C | ~19 | % of glucose C | as above | Results p. 486 | within_study_ratio |
| **furan rings + associated C6H2** | **~11** | % of total sugar C | as above | Results p. 486; structure F | within_study_ratio |
| **pyrrole rings + associated alkyl-C6** | **~8** | % of all glucose C | as above | Results p. 486; structure G | within_study_ratio |
| imidazolium + oxazolium + associated alkyl | ~6 | % of all glucose C | as above | Results p. 486 | within_study_ratio |
| HC3=C4H alkene unit | ~6 | % of all glucose C | as above | Results p. 486; structure H | within_study_ratio |
| structure A (N-acyl + ketal ring) | ~10 | % of all sugar C | as above | Results p. 487; Table 1 | within_study_ratio |
| structure B (open chain with glycine attached) | ~8 | % of all sugar C | as above | Results p. 487; Table 1 | within_study_ratio |
| structures Ea + Eb (ketone C2 with NCq C1) | ~4 | % of all sugar C | as above | Results p. 488; Table 1 | within_study_ratio |
| structure C (pyrrole-2-carbaldehyde with N-glycine) | ~1.5 | % of all C | as above | Results p. 487; Table 1 | within_study_ratio |
| structure D (acetyl ester from C2 + C4 fragmentation) | ~1 | % of all C | as above | Results p. 488; Table 1 | within_study_ratio |
| **glucose carbon backbones fragmented** | **between 25 and 50** | % | as above | Results p. 486, "Fragmentation" | **within_study_ratio** (author's own estimate, stated as a range) |
| dominant fragmentation site | C2 (C1+C5 or C2+C4) | — | as above | Results p. 486 | level_only |
| COO/NC=O marker intensity on C2 | ~7 | % of C2 | as above | Results p. 486 | within_study_ratio |
| alkyl share of C3 | ~52 | % of C3 | as above | Results p. 485 | within_study_ratio |
| alkyl share of C4+C5 | ~52 | % of C4+C5 | as above, **difference spectrum** (Flags 4) | Results p. 485 | within_study_ratio |
| alkyl share of C6 | 82 | % of C6 | as above | Results p. 485 | within_study_ratio |
| C6 remaining as OCH2 | 50 | % of C6 | as above | Results p. 489 | within_study_ratio |
| C4, C5 remaining as OCH | 25 | % each | as above, difference spectrum | Results p. 489 | within_study_ratio |
| nonprotonated C1 (150-20 ppm) | 29 | % of C1 | as above | Results p. 485 | within_study_ratio |
| alkyl-CH C1 | 20 | % of C1 | as above | Results p. 486 | within_study_ratio |
| C−CH2−C C1 | 5 | % of C1 | as above | Results p. 486 | within_study_ratio |
| **C1 forming new C−C bonds** | **~50** (abstract: "more than half") | % of C1 | as above | Results p. 486; Abstract | within_study_ratio |
| C6 in C−CH2−C | 16 | % of C6 | as above | Results p. 486 | within_study_ratio |
| C6 in aromatic rings, nonprotonated | 5 | % of C6 | as above | Results p. 486 | within_study_ratio |
| C1 and C2 directly N-bonded (140-115 ppm band) | ~85 | % each | as above | Results p. 486 | within_study_ratio (semiquantitative REDOR — Flags 5) |
| C3 directly N-bonded, same band | 9 | % | as above | Results p. 486 | within_study_ratio (same caveat) |
| protonated sp2 near 132 ppm, C3 and C4 | ~7 | % each | as above | Results p. 486 | within_study_ratio |
| C2's share of the ketone signal | 50 | % of ketone C | as above | Results p. 489 | within_study_ratio |
| largest single structure | < 15 | % of total intensity | as above, 13C1 sample | Results p. 489; Abstract | within_study_ratio |
| the 11 units together | > 50 | % of glucose C | as above | Abstract; Results p. 489 | within_study_ratio |
| **aldehyde bonded to a furan ring** | **not observed** | — | as above | Results p. 488 | level_only (**a measured null; bears on `hmf`**) |
| aldehyde bonded to a pyrrole ring | clearly detected | — | as above | Results p. 488; structure C | level_only |
| pyridines | insignificant, "very weakly only" | — | as above | Results p. 489 | level_only (a measured null) |
| isolated NCH3 / OCH3 glucose carbons | completely absent | — | as above | Results p. 486 | level_only (a measured null) |
| **neat glucose heated 125 C, 2 h, no amine** | no major structural change; three original peaks, "only broadened" | — | dry, 125 C, 2 h, glucose alone | Results p. 484; Figure 2c | level_only (**a measured null on caramelisation at this severity**) |
| insoluble vs soluble HMW fraction | similar percentages; insoluble has larger aromatic and smaller O-alkyl | — | same reaction | Results p. 484 | level_only (supporting figure is in the missing SI) |
| 9:1 vs 1:1 glucose:glycine | more ketones and furans, fewer pyrroles at 9:1 | — | same dry reaction, ratio varied | Results p. 483; Figure 2 | level_only (**no numbers given**) |
| implied sugar-side carbon per repeat unit if fragments split evenly | 4.5 to 5.25 | mol C | as above | derived from the 25-50 % range (mine) | derived_assumption (scoping only) |
| all NMR spectra, pulse sequences and the Figure 4 bar chart | — | — | — | Figures 1-10 | **figure_only** |

### How this bears on the trunk

**(a) Nothing here is a rate and nothing here is a benchmark.** One endpoint, one temperature,
one time. No constant transports into `parameters.py` or any lane.

**(b) It is the strongest evidence in the corpus that `FRAG_C` is chemically real.** 25-50 % of
glucose backbones broken, with C2 the dominant break point, means the trunk's practice of routing
unreported carbon into a visible accounting pool is closer to the chemistry than an assumption
that every step-9 event carries an intact C6.

**(c) It is a soft argument against the sugar half of `MELANOIDIN_REPEAT_UNIT_CARBON = 8`, in
the same direction as the amine half.** Both corrections lower the model's C/N, which is
currently floored at 8.0 while the measurements (Mundt 2004: 7.64 ± 0.21; Fang 2009's implied
~7.4-7.6) sit below it. **This paper cannot size the sugar-side correction** — "between 1/4 and
1/2" is too wide, and it does not say what the broken fragments do next.

**(d) The `hmf` null is worth recording against the B7 furanic block.** `species.py`'s B7 block
carries HMF on the trunk and `furanic.py` gives it two sources. Fang finds furan rings in the
polymer at ~11 % of sugar carbon but **no furan-bound aldehyde**, and argues C1 is split off
before incorporation. So HMF is a product, not a monomer, in this matrix — which is consistent
with how the repository already treats it (a measured species with its own sinks) and is a reason
not to add an HMF-to-melanoidin edge without evidence.

**(e) The neat-glucose control is a usable null for a caramelisation lane.** 125 C for 2 h of dry
glucose with no amine produced no structural change beyond amorphisation. Any model that
generates browning polymer from sugar alone at that severity is contradicted by a direct
measurement.

## 5. Flags

1. **The file is named `fang2010.pdf` but the paper is J. Agric. Food Chem. 2011, 59 (2),
   481-490.** The ACS copyright line reads 2010 and web publication was 28 December 2010. **Cite
   it as 2011** in any bibliography; the dossier keeps the file-name year in its own name to match
   the corpus convention.
2. **The Supporting Information is missing and it holds the paper's main quantitative table.**
   The SI note reads: "Contributions of glucose carbons to melanoidins, main features of the alkyl
   carbons in the glucose:glycine 1:1 melanoidins, and figure of quantitative 13C NMR spectra of
   insoluble melanoidins." The first item is precisely the per-carbon contribution table that
   would let this paper's percentages be used as anything more than text quotations, and the third
   is the only evidence behind the insoluble-vs-soluble similarity claim. **Fetch it before any of
   these numbers enters a registry.**
3. **"C1" means a different atom in this paper than in its companion.** Here C1-C6 are **glucose**
   carbons; in Fang 2009 C1 and C2 are the **glycine** carboxyl and methylene. The two papers share
   a sample, an author pair and a journal, and the collision is easy to make. Always say which
   molecule.
4. **The C4 and C5 numbers are a five-spectrum difference.** No 13C4 or 13C5 glucose was used;
   the C4+C5 spectrum is `13C6(uniform) − (13C1 + 13C2 + 13C3 + 13C6(single))`. Its error is the
   quadrature sum of five spectra each already carrying an 8 % constructed natural-abundance
   background subtraction. **The C4/C5 rows in section 4 are the least reliable in the paper**,
   and C4 and C5 are never resolved from each other at all.
5. **The nitrogen-proximity percentages are semiquantitative by the paper's own word.** "To obtain
   **semiquantitative** information about the distribution of carbons that are within two bonds
   from 15N, the 13C{15N} REDOR pulse sequences ... were applied." The ~85 % N-bonded figure for
   C1 and C2 and the 9 % for C3 carry that label; only the DP/echo integrals are quantitative.
6. **The 9:1 ratio series is described and never quantified.** More ketones, more furans, fewer
   pyrroles, more O−CH — every one of those is a direction with no number, read off Figure 2. It
   is the only handle in this paper on how melanoidin composition varies with the sugar:amine
   ratio, and it cannot be used. Cämmerer & Kroh's opposing claim (composition "little affected"
   from 4:1 to 1:1) is quoted in the introduction from ref. 12, not tested here.
7. **Several assignments are forward-referenced to a paper that is not in this corpus.** "Specific
   information on the nitrogen-containing rings can be obtained from 2D 15N-13C HSQC and 3D
   15N-13C-13C spectra **to be discussed in a future publication**", and again "confirmed by the
   analysis of 15N-13C-13C spectra, **detailed in a forthcoming publication**". The imidazolium
   and oxazolium assignments — and therefore the ~6 % figure attached to them — rest on that
   unpublished work. The same forward reference appears in Fang 2009. **Search for the third paper
   in the series.**
8. **Structure H's abundance rests on a background correction the paper does not detail.** "Taking
   background signals into account, this component represents about 6 % of all glucose C." No
   error and no method for that correction is given.
9. **Every abundance in this paper is quoted with "ca." or "about" and none carries an error
   bar.** The companion paper (Fang 2009) gives ± 2 % on its Table 2 and ± 3-6 % on Tables 3 and 4.
   This one gives nothing. Treat every percentage here as one significant figure.
10. **What this paper does not contain**: any rate constant; any activation energy; any time
    course or intermediate time; any temperature series; any water content or water-activity
    point; any combustion elemental analysis, C/N or H/C; any absorbance or extinction
    coefficient; any molecular-weight distribution beyond the dialysis cut-off; any yield; any
    quantitative measurement on the amino acid (that is the companion); any amine other than
    glycine; any sugar other than glucose; any real food.
11. **What to request from the authors**: (i) the Supporting Information, especially the
    "Contributions of glucose carbons to melanoidins" table; (ii) the forthcoming 15N-13C-13C
    paper that carries the imidazolium and oxazolium assignments; (iii) quantitative percentages
    for the **9:1** melanoidin, which would give a composition-versus-ratio point the corpus does
    not have; (iv) a combustion CHN analysis of these same samples — the same request as for the
    companion, and the one that would let both papers speak to the C/N question directly.
12. **Registry gaps against `data/keys/compounds.yml`**: glucose, glycine and any melanoidin pool
    are unkeyed. The furan and pyrrole *rings* this paper finds inside the polymer are not the
    volatile furans (`furan`, `furfural`, `hmf`, `2_acetylfuran`, `2_pentylfuran`) the registry
    keys — they are polymer-bound and unextractable — and the two must never be totalled together.
    The trunk's `MEL_C`, `MEL_N` and `FRAG_C` remain network-local names with no registry surface.
