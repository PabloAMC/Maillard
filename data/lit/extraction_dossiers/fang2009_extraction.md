# Fang 2009 — EXTRACTION (equimolar glucose + glycine melanoidins made two ways — freeze-dried solid at 125 C for 2 h, and 0.056 M aqueous at pH 8, 100 C for 1 week — with 13C- and 15N-labelled glycine, dissected by quantitative solid-state NMR: how much glycine carbon is in the polymer, how much left as CO2, and what fraction of the C1-C2 and C2-N bonds survive)

### THE ACCOUNTING PAPER FOR THE AMINE: it does not measure a melanoidin C/N, but it prints the two numbers a C/N is made of — **glycine carbon is 24 % (dry) / 22 % (solution) of all melanoidin carbon against 25 % in the reactants**, and **27 ± 4 % (dry) / 33 ± 4 % (solution) of glycine reactant loses its carboxyl as CO2** — and its NMR says the nitrogen almost always stays attached to the amino acid's own C2, which is the assumption the trunk's `MEL_N` repeat-unit count rests on.

**Source on disk:** `data/articles/fang2009.pdf` (11 pp., J. Agric. Food Chem. 2009, 57 (22), 10701-10711).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/fang2009.txt`). **Table 1 came through the text layer clean.** Tables 2, 3
and 4 did **not** — they are graphics (Table 2's column headers are drawn chemical structures)
and the text layer returned only their captions and footnotes. Page 10709 was therefore
**rendered at 200 dpi (`pdftoppm -r 200 -f 9 -l 9`) and read from the image**; all three are
re-typed in full below, with Table 2's structural column headers written out in words and
formulae. Figures 1-10 are NMR spectra and are **figure_only** throughout; the paper's own
practice is to print every integrated percentage in the text or in a table, so almost nothing
useful is trapped in a figure here. **Supporting Information (Figures S1-S5) exists and is NOT
on disk** — see Flags 9. Repo status before this dossier: `fang2009.pdf` has **no extraction
dossier** and is not cited anywhere in `src/`.

## 0. Identity

| field | value |
|---|---|
| Title | "Fate of the Amino Acid in Glucose-Glycine Melanoidins Investigated by Solid-State Nuclear Magnetic Resonance (NMR)" |
| Authors | Xiaowen Fang and Klaus Schmidt-Rohr (corresponding), Department of Chemistry, Iowa State University, Ames, Iowa 50011 |
| Venue | J. Agric. Food Chem. 2009, 57 (22), 10701-10711. Received 16 June 2009, revised 30 September 2009, accepted 3 October 2009, web 29 October 2009. NSF grant CHE-0138117 |
| DOI | 10.1021/jf9020587 |
| Naming | "melanoidin" = the **soluble high-molecular-weight (HMW) dialysis retentate** at MWCO 6000-8000. "C1" = the carboxyl carbon of glycine; "C2" = the methylene carbon; "N" = the amine nitrogen. "Strecker degradation" is used **narrowly and deliberately**: only fragmentation that loses C1 as CO2. Fragmentation that keeps C1 is called "other degradation" and deamination is called separately — see Flags 3, because this is the definitional point on which the paper disagrees with its predecessors |
| Lineage | melanoidins made to the **European COST Action 919** standard protocol for the solvent-free glucose-glycine reaction (ref. 23, Ames 2002), following Adams 2003 and Tehrani 2002 for the dry route and Hedges 1978 / Ikan 1986 for the solution route. Argues against Cämmerer & Kroh 1995 (the enamine-C claim), against Kato & Tsuchida 1981 and Yaylayan & Kaminsky 1998 (imine/enamine models), and against Benzing-Purdie & Ripmeester 1987 (which it says over-counted Strecker degradation) |
| Companion on disk | **`fang2010.pdf` is the direct isotope companion** (J. Agric. Food Chem. 2011, 59:481, the labelled-glucose half of the same programme) — dossier `fang2010_extraction.md` |
| Companions on disk | `mundt2004_extraction.md` (the Leeds radiochemical + microanalysis measurement of the same system's C/N — the paper this one must be read against), `cammerer1994_extraction.md` (this paper's ref. 17, whose enamine claim it disproves), `adams2008_extraction.md` (the same Ghent thermal-degradation school as its refs. 19-20) |

## 1. Why it matters

The trunk books unmeasured amine nitrogen into `MEL_N` and unmeasured carbon into `FRAG_C`
(`src/kinetic_core/species.py`, the two pool entries and the note above them). Neither pool is
measured by any experiment in the fit corpus. **This paper is a direct measurement of where the
amino acid actually goes**, on the same chemistry (glucose + glycine, 1:1) and by a method that
sees bonds rather than bulk composition. Four things it settles or constrains:

**(a) The nitrogen-per-repeat-unit assumption.** `melanoidin_repeat_units` in `species.py`
returns `MEL_N`, not `MEL_C/8`, and the docstring's argument is that "every step-9 event
contributes exactly one nitrogen, and every carbon-only addition to the polymer ... grows an
existing unit rather than creating a new one". That argument needs the amino acid's nitrogen to
stay with the polymer and stay countable. Fang measures it: **75 % (dry) / 63 % (solution) of
melanoidin nitrogen is still directly bonded to the glycine's own C2 carbon**, and the 2D HSQC
shows that nitrogen sits in a dozen different chemical environments "usually without breaking
the bond to glycine C2". The nitrogen is not lost and it is not scrambled off its own carbon.
The counting assumption survives.

**(b) The carbon side of the same repeat unit, which is where the trunk is wrong.** The trunk's
`MELANOIDIN_REPEAT_UNIT_CARBON = 8` books an **intact** glycine (2 C) onto a 3-deoxyglucosone
(6 C). Fang measures the depletion directly: glycine carbon is **25 % of reactant carbon** by
construction (glucose C6 + glycine C2) but only **24 % (dry) / 22 % (solution) of the carbon in
the melanoidin**, and the missing ~2-3 percentage points are C1 leaving as CO2. In molecule
terms, **27 ± 4 % (dry) and 33 ± 4 % (solution) of the glycine reactant loses C1 as CO2**. So
some, but not most, of the amine arrives decarboxylated. This is the **same direction** as
Mundt & Wedzicha's radiochemistry but a **much smaller magnitude** (Mundt: ~2/3 of incorporated
glycine decarboxylated), and section 4 does the arithmetic that reconciles them into a C/N.

**(c) `FRAG_C` has a competitor the trunk does not model.** `species.py` says the unassigned
fragment carbon pool exists so that "the total carbon balance closes exactly". Fang's headline
is that **carbon leaves the liquid as CO2** — one carbon per Strecker event — and CO2 is not a
species in the trunk and is not `FRAG_C` (which is dissolved fragment carbon, nominally still in
the pot). Any future rewrite of step 9 that decarboxylates has to decide where that carbon goes.

**(d) A structural veto list.** The paper rules several structures out of the polymer by direct
measurement: **no enamines, no imines, no pyrazines, no pyridines** in any significant amount,
in either preparation. That matters because `data/keys/compounds.yml` carries `pyrazines`,
`methylpyrazine`, `2_3_dimethylpyrazine`, `2_5_`, `2_6_`, `trimethylpyrazine`,
`tetramethylpyrazine` and `2_ethyl_3_5_dimethylpyrazine` as volatile *products*. Fang's finding
is not that pyrazines are absent from a Maillard pot — they plainly are not — but that they are
**not part of the high-polymer melanoidin**, so nitrogen routed to `MEL_N` and nitrogen routed
to a pyrazine are separate sinks that must not be conflated. See Flags 5.

What this paper does NOT give the repository: any rate constant, any activation energy, any time
course, any temperature series, any elemental analysis, any C/N, any absorbance, and any
molecular weight beyond the dialysis cut-off. It is a bond-accounting paper at two fixed
endpoints.

## 2. Methods as they matter to a model

Two melanoidins are made, and **they differ in almost every dimension at once** — phase,
temperature, time, pH and yield. Nothing here is a controlled comparison of one variable
(Flags 1).

- **Dry reaction (the COST Action 919 protocol).** Equimolar glucose and glycine, **5.6 mmol
  each**, dissolved in 15 mL E-pure water in a scintillation vial, **freeze-dried for 1 day**
  (Labconco Freezon 4.5), then the dried solid heated in a **closed, preheated oven at 125 C for
  2 h** (VMR Scientific Products model 1430). "After 2 h of heating, the color of the sample
  turned to dark brown and its volume increased." Cooled to room temperature in a desiccator.
  **This is a low-moisture matrix**, which is the axis on which most of the corpus is silent.
- **Solution reaction.** A 500 mL flask containing **100 mL of D-glucose and glycine, 0.056 M
  each**, in a **pH 8 buffer**, sealed and heated in a preheated oven at **100 C for 1 week**.
  "After the first 48 h, the color of the solution turned brown and **the pH dropped to ~5**.
  Na2HPO4·7H2O was added to adjust the pH back to **8.5** at room temperature. Another 5 days of
  reaction caused **the pH of the solution to drop to < 6**." So the pH history is 8 → 5 → 8.5 →
  < 6 over 7 days, with one manual intervention. Read "pH 8" as a starting condition only.
- **Work-up, both routes.** The dry product was redissolved in 100 mL E-pure water with 1 h of
  magnetic stirring "to dissolve as much sample as possible", filtered twice through Whatman 41
  ashless paper, and the filtrate transferred to 30 cm of Fisherbrand regenerated-cellulose
  dialysis tubing, **MWCO 6000-8000**, 5.10 mL/cm. Submerged in **2 L of E-pure water at 4 C**,
  water changed every 12 h, **four changes, 48 h total**. Retentate freeze-dried.
- **Yields, and which fraction was analysed.** Dry route: the soluble HMW fraction studied here
  "accounts for **~7 % of reactant mass**"; the **insoluble fraction is the main product
  (~50 % yield)** and "shows 13C and 15N NMR spectra **similar** to those of the HMW fraction
  studied here". Solution route: HMW yield **~20 %**, and **no insoluble fraction was obtained**.
  So the measured object is a 7 % (dry) or 20 % (solution) minority cut of the product mass, with
  a stated but unquantified similarity claim covering the missing 50 %.
- **Labels used.** D-glucose U-13C6 (99 %); glycine 1-13C, 2-13C, 1,2-13C2, 15N, and
  glycine-2-13C-15N. Each label answers a different question: `1,2-13C2` + J-modulation gives
  the survival of the **C1-C2 bond**; `2-13C-15N` + REDOR in both detection directions gives the
  survival of the **C2-N bond**; `13C6` glucose against `15N` glycine separates glucose carbon
  from glycine carbon.
- **The quantitative measurement.** 13C **direct-polarisation** (DP) with a Hahn echo at 14 kHz
  MAS, Bruker DSX 400 (100 MHz 13C, 40 MHz 15N), 4 mm triple-resonance MAS probe. Recycle delays
  were set from measured T1 so that "all carbons are fully relaxed" — this is what makes the
  percentages quantitative, and it is why the DP spectra, not the cross-polarisation ones, carry
  the numbers. Cross-polarisation (CP/TOSS) spectra are used for routine analysis and for the
  minor signals, and are explicitly **not** quantitative for 15N (nonprotonated N is
  under-represented; see the Nt:NH discussion).
- **Corrections applied by the authors.** (i) Natural-abundance 13C background, "~5 % of the
  total intensity", subtracted. (ii) A scaling factor of **1/0.75** for 15N{1H} dephasing,
  calibrated on 15N-t-BOC-L-proline. (iii) A **1.12 (dry) / 1.2 (solution)** conversion factor
  from "fraction of glycine carbon in the product" to "fraction of glycine reactant carbon",
  derived from the C2 excess. Every number in Tables 3 and 4 has been through that conversion;
  the numbers in Table 2 and in Figure 2 have not.
- **What is not measured.** No elemental analysis, no molecular weight distribution, no
  absorbance, no yield of any small molecule, no time course, no temperature series, and **no
  measurement of nitrogen loss** — the paper tracks where nitrogen *is*, never how much of it
  left as ammonia.

## 3. Tables re-typed

### Table 1. "Approximate Amounts of the Forms of Nitrogen in the Two Melanoidins Studied, as a Percentage of Total Glycine Nitrogen"

| sample | aromatic N (pyrrolic) (%) | amide (NC=O) (%) | amine (%) | N−H (%) | nonprotein Nt (%) |
|---|---|---|---|---|---|
| dry | 39 | 33 | 28 | 21 ± 2 | 79 ± 2 |
| solution | 31 | 53 | 16 | 15 ± 5 | 85 ± 5 |

Header exactly as printed, including the column label "**nonprotein Nt**", which from the text
plainly means **nonprotonated** tertiary nitrogen (Flags 6). The first three columns sum to 100
in both rows (mine: 39 + 33 + 28 = 100; 31 + 53 + 16 = 100), and the last two are a separate
protonated/nonprotonated partition of the same total.

### Table 2. "Structural Units into which Glycine C2 (in Bold; with C1 in Italics and Underlined) Has Been Transformed during the Maillard Reaction and Their Percentages of Total Glycine Carbon in the Melanoidins"

Read from the 200 dpi page render. The column headers are drawn structures; I give them as
formulae, preserving which atom the percentage is counting. Footnote as printed: "Percentages of
N−C moieties from glycine are in bold; those of C1−C2 moieties are in italics. Thus, values for
N−C1−C2 moieties are in bold italics. **Error margins: ±2 %.** † Uncertain assignment."

**Upper block:**

| | N−**CH2**−*COO* (intact glycine) | Other N−**CH2** | O=C−N**CH3** | C**CH2**C | C**CH3** | *C***CH**† |
|---|---|---|---|---|---|---|
| Dry | **33 %** | 8 % | 3 % | 3 % | 2 % | 1.2 % |
| Sol. | **25 %** | 9 % | 3 % | 5 % | 2 % | 4 % |

**Lower block:**

| | pyrrolic ring **C**−*COO* | imidazolium N⟨ring⟩N−**CH2** | Other ⟨ring⟩**N**− | Other aromatic **C** | **C**OO | Total **C2** |
|---|---|---|---|---|---|---|
| Dry | 2 % | 1.4 + 1.4 % | 0 % | 0 % | 1 % | **56 %** |
| Sol. | 3 % | 0.7 + 0.7 % | 3.5 % | 1 % | 3 % | **60 %** |

The "1.4 + 1.4 %" and "0.7 + 0.7 %" are printed exactly so — the imidazolium column carries two
distinguishable carbon positions in the ring and the paper adds them without combining them.
The "Total C2" column is the sum of all C2 assignments and reproduces the 56 % / 60 % figures
quoted throughout the text.

### Table 3. "Fractions of C in Glycine Reactants Incorporated into Melanoidin after Dry Reaction or Lost as CO2"

Two overlapping partitions of the same 100 % of glycine reactant, printed as two banded rows.

*Band 1 — by which bonds survive:*

| C1−C2 (spanning) | | C2−N (spanning) | | | | Total |
|---|---|---|---|---|---|---|
| **70 ± 5 %** | | **75 ± 4 %** | | | | |
| Isolated C1 | Other C1-C2 | C1-C2-N | Other C2-N | Isolated C2 | C1 lost | |
| 3 % | 6 % | **59 ± 4 %** | 11 % | 7 % | 15 % | **101 %** |

*Band 2 — by which degradation happened:*

| Other degradation | No degr. | Other | Strecker degradation | Total |
|---|---|---|---|---|
| 9 % | **59 ± 4 %** | 3 % | **27 ± 4 %** | **99 %** |

Footnote as printed: "'Other C1−C2' refers to C1−C2 pairs that are not bonded to N (deamination
products), 'other C2−N' to C2−N fragments not bonded to C1 (decarboxylation products), 'isolated
C2' to C2 bonded to neither C1 nor N, and 'other' or 'other degradation' to degradation which,
unlike Strecker degradation, does not involve loss of C1. The sums of the percentages within the
two full rows are compatible with 100 % within the error margins, while the sum of the first row
exceeds 100 % because the same C2 carbon (in C1−C2−N moieties) can contribute to both
categories."

### Table 4. "Fractions of C in Glycine Reactants Incorporated into Melanoidin after Solution Reaction or Lost as CO2"

*Band 1:*

| C1−C2 (spanning) | | C2−N (spanning) | | | | Total |
|---|---|---|---|---|---|---|
| **54 ± 4 %** | | **62 ± 3 %** | | | | |
| Isolated C1 | Other C1-C2 | C1-C2-N | Other C2-N | Isolated C2 | C1 lost | |
| 6 % | 12 % | **42 ± 4 %** | 14 % | 10 % | 17 % | **101 %** |

*Band 2:*

| Other degradation | No degr. | Other | Strecker degradation | Total |
|---|---|---|---|---|
| 18 % | **42 ± 4 %** | 6 % | **33 ± 4 %** | **99 %** |

Footnote: "Percentages within the two full rows add up to 100 % within the error margins. The
terminology is the same as in Table 3."

### Numbers printed in the running text and figure captions

| quantity | dry reaction | solution reaction | where |
|---|---|---|---|
| glycine C as % of total **reactant** carbon | 25 % | 25 % | Results (structural, 1:1 molar, C6 + C2) |
| **glycine C as % of total carbon in the melanoidin** | **24 %** | **22 %** | Results, "Extent of Glycine Loss"; Figure 2 caption |
| of that glycine carbon, the share that is C1 | 44 % | 40 % | Results; Figure 2 caption |
| of that glycine carbon, the share that is C2 | 56 % | 60 % | Results; Figure 2 caption |
| C1 as % of **all** melanoidin carbon | ~10 % (0.42 × 23 %) | — | Results (the paper works the ~23 % average) |
| C2 as % of **all** melanoidin carbon | ~13 % (0.58 × 23 %) | — | Results |
| carbon lost to CO2, as % of all melanoidin carbon | ~3 % | — | Results |
| glycine C1 chemical shift | single peak at **172 ppm** | **174 ppm** | Results, "Fate of Glycine C1"; Figure 2b,e |
| C1 CSA powder pattern | resembles **esters** significantly | resembles typical **COOH** | Results, from SUPER (Figure S5 — not on disk) |
| C1 remaining in COO moieties | ≥ 90 % | ≥ 90 % | Abstract |
| amides among C1 forms | ~5 % | ~5 % | Abstract; Results |
| NCH2 resonance near 50 ppm, as % of all glycine C | 44 % | 39 % | Results, "Fate of Glycine C2" |
| **C1−C2 bond intact, as % of glycine C in the melanoidin** | **78 ± 6 %** | **65 ± 5 %** | Results, "Fate of the C1-C2 Bond" |
| C1−C2 bond intact, as % of glycine **reactant** | **70 %** (= 78/1.12) | **54 %** (= 65/1.2) | Results; Tables 3, 4 |
| C1 carbons still bonded to C2 (J-dephasing) | ~95 % | ~87 % | Results |
| C2 carbons still bonded to C1 (J-dephasing) | 64 % | 50 % | Results |
| the same, cross-checked via quantitative 13C | 42 % vs 36 % of all glycine C (**the paper calls this discrepancy unexplained**) | 35 % vs 30 % | Results |
| **glycine C2 bonded to N** (13C{15N} REDOR) | **82 %** | **62 %** | Results, "Fate of the C-N Bond" |
| **N bonded to glycine C2** (15N{13C} REDOR) | **75 %** | **63 %** | Results; Figure 6 caption |
| the C2−N / N−C2 asymmetry | 82 vs 75, a 7 % gap attributed to one N bonding two C2 (imidazolium) | 62 vs 63, no gap | Results |
| imidazolium N bonded to three glycine carbons | **1.4 % of all N** | "significantly smaller" | Results (from Figure 8) |
| amide N peak at 110 ppm | **6.2 % of all N**, two-thirds protonated | — | Results |
| glycine carbons forming peptide bonds (H−N−C=O) or O=C−NH2 end groups | ~6 % | — | Results |
| C2=O carbons | ~1.5 % | — | Results (from Figure 2) |
| inferred NC1=O amide contribution | ~5 % | — | Results |
| pyrrolic N as % of all N | < 40 % | < 40 % | Results (both melanoidins) |
| N in **regular** pyrrole rings (~150 ppm) | < 20 % | < 20 % | Results |
| amide + amine N together | > 50 % | > 50 % | Results |
| **Nt : NH ratio, 13C-detected** | **79 : 21** | **92 : 8** | Results, "(Non)protonation of Nitrogen" |
| Nt : NH ratio, 15N-detected (authors call it under-estimated) | 75 : 25 | 75 : 25 | Results |
| nonprotonated N overall | > 75 %, and > 90 % of pyrrole N | > 75 % | Results; Abstract says > 78 % |
| **Strecker degradation, % of glycine reactant** | **27 ± 4 %** (a first route gave 21 %, revised upward) | **33 ± 4 %** | Results; Tables 3, 4 |
| fragmentation **without** C1 loss, % of glycine reactant | 5 % → 9 % as "other degradation" | 12 % → 18 % | Results; Tables 3, 4 |
| glycine remaining fully intact (C1-C2-N) | **59 %** | **42 %** | Results; Tables 3, 4; Conclusion |
| "C1 lost" cross-check | 18 % − 2 % = 16 % one way, 21/2 = 11 % the other; **27 ± 4 % adopted** | 23 % / 17 % | Results |
| C2−N fragments not bonded to C1 | (8 + 3 + 1.4)/1.12 = 11 % | 14 % | Results |
| C1−C2 fragments not bonded to N | 2(1.2 + 2)/1.12 = 5.7 % | 12 % | Results |
| isolated C2 (bonded to neither C1 nor the original N) | (3 + 2 + 1.4 + 1)/1.12 = 6.6 % | 10 % | Results |
| **structures ruled OUT of the polymer** | enamines, imines, pyrazines, pyridines — no significant 15N signal between 250-350 ppm or 300-370 ppm, no 13C at 146 or ~150 ppm | same | Results, three separate sections; Abstract |
| glucose HMW yield | ~7 % of reactant mass (soluble HMW); ~50 % insoluble, spectra "similar" | ~20 % HMW, **no** insoluble fraction | Experimental |

**Every NMR spectrum (Figures 1-10) is FIGURE-ONLY.** The paper prints its integrals in text and
tables, so no number below is read off a figure.

### Arithmetic on the printed numbers (all mine)

**1. The two 13C shares and the CO2 loss reconcile.** Dry: 25 % reactant → 24 % product is a
1 percentage-point drop; solution 25 % → 22 % is 3 points. The paper's own reconciliation uses
the ~23 % average and gets "~3 %". **My check on the dry case**: if 27 % of glycine reactant
loses one of its two carbons, glycine carbon in the surviving pool falls by 27/2 = **13.5 %
relative**, i.e. 25 % → 21.6 % of a *fixed* total carbon. The observed dry value is 24 %, not
21.6 %. **The two do not close on their own**, and the reason is that glucose carbon is also
being lost from the polymer at the same time (as CO2, water and volatiles) so the denominator is
not fixed. This is not an error in the paper — it is a warning that these percentages are
*shares of the retained polymer*, never absolute recoveries. Nothing here supports a mass
balance on the pot.

**2. Strecker's share among the glycine that made it into the polymer (mine).** Table 3's
"C1 lost" is **15 %** of glycine reactant carbon (dry) / **17 %** (solution). Since one lost C1
is half of one glycine's carbon, 15 % of carbon corresponds to 30 % of glycine molecules, which
brackets the adopted 27 ± 4 %. Consistent.

**3. An implied melanoidin C/N from this paper alone (mine, assumption-laden — this is the row
that connects to Mundt 2004).** Take the dry reaction. The conversion factor 1.12 says the
glycine reactant carried 1.12 times the carbon that its representatives carry in the melanoidin,
so **each glycine-derived unit in the polymer carries 2/1.12 = 1.786 carbons**. That glycine
carbon is 24 % of all melanoidin carbon, so **total melanoidin carbon per glycine-derived unit =
1.786 / 0.24 = 7.44**. If each such unit brings exactly one nitrogen, **C/N ≈ 7.4 (dry)**. The
same arithmetic on the solution numbers: 2/1.2 = 1.667, / 0.22 = **C/N ≈ 7.6 (solution)**.

  Assumptions, stated plainly: (i) every glycine-derived unit in the polymer retains exactly one
  nitrogen — Table 3 says 7 % of glycine reactant carbon is "isolated C2" with no N, which pushes
  the true C/N **up**, while any nitrogen that entered without carbon pushes it **down**;
  (ii) all melanoidin nitrogen is glycine nitrogen, which is true by construction here (glycine
  is the only nitrogen source); (iii) the retained HMW fraction is representative. **This is my
  arithmetic, not the paper's, and the paper prints no C/N anywhere.** Its value is that it
  lands on **7.4-7.6 by a completely independent technique**, against Mundt & Wedzicha's
  microanalytical **7.64 ± 0.21** and radiochemical **7.61** on the same chemistry. Three methods
  in two laboratories converging near 7.5 is the strongest statement the corpus can currently
  make about melanoidin C/N.

**4. The trunk's floor, restated against this paper (mine).** `MELANOIDIN_REPEAT_UNIT_CARBON = 8`
gives a model floor of C/N = 8.0. Fang's dry-reaction unit carries 1.786 glycine carbons instead
of the trunk's 2.000. Substituting only that term into the trunk's own repeat unit — 6 from
3-deoxyglucosone plus 1.786 from the amine — gives **7.79**, still above Mundt's 7.64 but below
the trunk's floor. **The gap between the trunk and the measurements is at least partly the
decarboxylation this paper quantifies**, and the size of Fang's effect (0.21 C per N) is about a
third of the size of Mundt's (0.70 C per N). The two papers agree on the sign and disagree on
the magnitude (Flags 3).

**5. The two preparations differ by more than the phase (mine).** Every degradation index is
worse in solution: Strecker 27 → 33 %, other degradation 9 → 18 %, intact glycine 59 → 42 %,
C1−C2 survival 70 → 54 %, C2−N survival 75 → 62 %. But the solution run is **1 week at 100 C**
and the dry run is **2 h at 125 C** — a ~84-fold difference in time against a 25 C difference in
temperature. **Nothing in this paper separates "wet vs dry" from "long vs short".** Do not read
these as a water-activity effect.

**6. Nitrogen environments are not conserved between the two routes (mine, from Table 1).** Amide
nitrogen goes 33 → 53 % and amine nitrogen 28 → 16 % from dry to solution, while pyrrolic goes
39 → 31 %. The nitrogen is present in both cases but in materially different chemistry. A model
that carries nitrogen as a single `MEL_N` number cannot see this, and does not need to — but any
future claim that `MEL_N` "is" a particular functional group is unsupported.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Glucose, glycine, CO2, ammonia and
melanoidins are all **absent**; the registry carries no Maillard reactants and no polymer pool.
The registry does carry eight pyrazine ids (`pyrazines`, `methylpyrazine`,
`2_3_dimethylpyrazine`, `2_5_dimethylpyrazine`, `2_6_dimethylpyrazine`, `trimethylpyrazine`,
`tetramethylpyrazine`, `2_ethyl_3_5_dimethylpyrazine`) plus `2_acetyl_1_pyrroline` — this paper
bears on all of them only in the negative sense of Flags 5.

Conditions shared by every "dry" row: **freeze-dried equimolar glucose + glycine, 5.6 mmol each,
closed preheated oven, 125 C, 2 h, low moisture; soluble HMW fraction (~7 % of reactant mass)
after dialysis at MWCO 6000-8000.** Conditions shared by every "solution" row: **0.056 M glucose
+ 0.056 M glycine, initial pH 8 (falling to < 6, readjusted once to 8.5 at 48 h), sealed flask,
100 C, 7 days; HMW fraction ~20 % yield, same dialysis.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **glycine carbon as a share of melanoidin carbon** | **24** | % of total C | dry | Results p. 10703; Figure 2 caption | **elemental_analysis** (by quantitative 13C DP NMR, not by combustion) |
| " | **22** | % of total C | solution | same | elemental_analysis |
| glycine carbon in the reactants | 25 | % of total C | both, structural | Results | derived_assumption (1:1 molar, C6 + C2) |
| C1 share of glycine carbon in the melanoidin | 44 / 40 | % | dry / solution | Results; Figure 2 caption | within_study_ratio |
| C2 share of glycine carbon in the melanoidin | 56 / 60 | % | dry / solution | Results; Figure 2 caption | within_study_ratio |
| **glycine reactant undergoing Strecker degradation (C1 lost as CO2)** | **27 ± 4 / 33 ± 4** | % of glycine reactant | dry / solution | Tables 3, 4 band 2 | **within_study_ratio** |
| glycine reactant fragmenting without C1 loss ("other degradation") | 9 / 18 | % of glycine reactant | dry / solution | Tables 3, 4 band 2 | within_study_ratio |
| **glycine reactant remaining fully intact (C1-C2-N)** | **59 ± 4 / 42 ± 4** | % of glycine reactant | dry / solution | Tables 3, 4 | within_study_ratio |
| C1−C2 bond surviving | 78 ± 6 / 65 ± 5 | % of glycine C **in the melanoidin** | dry / solution | Results p. 10705 | within_study_ratio |
| the same, on a reactant basis | 70 ± 5 / 54 ± 4 | % of glycine **reactant** | dry / solution | Tables 3, 4 | within_study_ratio |
| **C2−N bond surviving** | **75 ± 4 / 62 ± 3** | % of glycine reactant | dry / solution | Tables 3, 4 | within_study_ratio |
| glycine C2 still bonded to N | 82 / 62 | % of glycine C2 | dry / solution | Results p. 10706, 13C{15N} REDOR | within_study_ratio |
| **melanoidin N still bonded to glycine C2** | **75 / 63** | % of total N | dry / solution | Results p. 10706, 15N{13C} REDOR; Figure 6 caption | within_study_ratio |
| C1 lost, as carbon | 15 / 17 | % of glycine reactant carbon | dry / solution | Tables 3, 4 band 1 | within_study_ratio |
| isolated C1 | 3 / 6 | % of glycine reactant C | dry / solution | Tables 3, 4 | within_study_ratio |
| "other C1-C2" (deamination products) | 6 / 12 | % | dry / solution | Tables 3, 4 | within_study_ratio |
| "other C2-N" (decarboxylation products) | 11 / 14 | % | dry / solution | Tables 3, 4 | within_study_ratio |
| isolated C2 (no C1, no original N) | 7 / 10 | % | dry / solution | Tables 3, 4 | within_study_ratio |
| intact N−CH2−COO unit | 33 / 25 | % of total glycine C in the melanoidin | dry / solution | Table 2 upper | within_study_ratio (± 2 %) |
| other N−CH2 | 8 / 9 | % of total glycine C | dry / solution | Table 2 | within_study_ratio (± 2 %) |
| O=C−NCH3 | 3 / 3 | % of total glycine C | dry / solution | Table 2 | within_study_ratio (± 2 %) |
| C−CH2−C | 3 / 5 | % of total glycine C | dry / solution | Table 2 | within_study_ratio (± 2 %) |
| C−CH3 | 2 / 2 | % of total glycine C | dry / solution | Table 2 | within_study_ratio (± 2 %) |
| C−CH (uncertain assignment) | 1.2 / 4 | % of total glycine C | dry / solution | Table 2, footnote † | within_study_ratio (author-flagged uncertain) |
| pyrrolic C−COO | 2 / 3 | % of total glycine C | dry / solution | Table 2 lower | within_study_ratio (± 2 %) |
| imidazolium N−CH2 (two positions) | 1.4 + 1.4 / 0.7 + 0.7 | % of total glycine C | dry / solution | Table 2 lower | within_study_ratio (± 2 %) |
| other ring N− | 0 / 3.5 | % of total glycine C | dry / solution | Table 2 lower | within_study_ratio (± 2 %) |
| other aromatic C | 0 / 1 | % of total glycine C | dry / solution | Table 2 lower | within_study_ratio (± 2 %) |
| COO (from C2) | 1 / 3 | % of total glycine C | dry / solution | Table 2 lower | within_study_ratio (± 2 %) |
| **aromatic (pyrrolic) N** | **39 / 31** | % of total glycine N | dry / solution | Table 1 | within_study_ratio |
| **amide (NC=O) N** | **33 / 53** | % of total glycine N | dry / solution | Table 1 | within_study_ratio |
| **amine N** | **28 / 16** | % of total glycine N | dry / solution | Table 1 | within_study_ratio |
| protonated N (N−H) | 21 ± 2 / 15 ± 5 | % of total glycine N | dry / solution | Table 1 | within_study_ratio |
| **nonprotonated (tertiary) N** | **79 ± 2 / 85 ± 5** | % of total glycine N | dry / solution | Table 1 | within_study_ratio |
| Nt : NH by 13C detection | 79 : 21 / 92 : 8 | ratio | dry / solution | Results p. 10708 | within_study_ratio |
| imidazolium N bonded to three glycine C | 1.4 | % of all N | dry | Results p. 10707 | within_study_ratio |
| amide N at 110 ppm | 6.2 | % of all N | dry | Results p. 10707 | within_study_ratio |
| glycine C forming peptide bonds / O=C−NH2 ends | ~6 | % | dry | Results p. 10707 | within_study_ratio |
| C1 remaining in COO | ≥ 90 | % of incorporated C1 | both | Abstract | level_only |
| enamines, imines, pyrazines, pyridines in the polymer | **not detected** | — | both | Results, three sections; Abstract | level_only (a measured null) |
| soluble HMW yield | ~7 / ~20 | % of reactant mass | dry / solution | Experimental | level_only |
| insoluble fraction | ~50 % (spectra "similar") / **none formed** | % yield | dry / solution | Experimental | level_only |
| **implied melanoidin C/N** | **~7.4 (dry), ~7.6 (solution)** | mol C per mol N | as above | derived from the 24 %/22 % shares and the 1.12/1.2 factors (mine) | **derived_assumption** (assumes 1 N per glycine-derived unit; the paper prints no C/N) |
| glycine C per glycine-derived polymer unit | 1.786 / 1.667 | mol C | dry / solution | derived from the 1.12 / 1.2 factors (mine) | derived_assumption |
| all NMR spectra | — | — | — | Figures 1-10 | **figure_only** |

### How this bears on the trunk's melanoidin pools

**(a) It supports keeping `MEL_N` as the repeat-unit counter.** 75 % (dry) of melanoidin nitrogen
is still bonded to the glycine C2 it arrived with, and the HSQC shows that nitrogen is "an active
site for reactions incorporating glycine into melanoidins, but usually without breaking the bond
to glycine C2". The one nitrogen per amine event that `melanoidin_repeat_units` counts is a real,
persistent object. The one structural exception the paper finds — an imidazolium nitrogen bonded
to **two** glycine C2 carbons — is quantified at **1.4 % of all N** in the dry reaction, i.e.
negligible for a counting argument.

**(b) It is evidence against the carbon side of `MELANOIDIN_REPEAT_UNIT_CARBON = 8`, but weaker
evidence than Mundt 2004.** Fang's effect size is 0.21 carbon per nitrogen (the amine brings
1.786 C, not 2.000); Mundt's is 0.70. Both point the same way — the model's floor of 8.0 is too
high — and neither is a direct measurement of a bulk C/N in a form the trunk can consume without
the assumptions in section 3, arithmetic 3.

**(c) It does not touch a single rate.** No constant here transports into
`parameters.py` or any lane. The paper's endpoints are 2 h at 125 C and 7 days at 100 C, with no
intermediate times. Nothing here can be a fit target or a hold-out.

**(d) The CO2 sink is the honest gap.** 27-33 % of glycine molecules lose a carbon to a gas the
trunk has no species for. Whatever the repository does about this, it should be recorded as a
declared gap rather than absorbed into `FRAG_C`, which `species.py` defines as carbon "leaving a
measured step in an unmeasured co-product" and not as a gas leaving the system.

## 5. Flags

1. **The two melanoidins differ in phase, temperature, time and pH simultaneously.** Dry: 125 C,
   2 h, no added water. Solution: 100 C, 168 h, pH 8 → < 6 with one manual readjustment. The
   dossier reports the two columns side by side because the paper does; **they are not a
   controlled water-activity comparison** and must never be quoted as one.
2. **The analysed material is a minority fraction, and the majority is only claimed to be
   similar.** In the dry reaction the soluble HMW retentate is ~7 % of reactant mass while the
   insoluble fraction is ~50 %, and the only evidence offered that they match is that the
   insoluble one "shows 13C and 15N NMR spectra **similar** to those of the HMW fraction studied
   here" — no spectra of it are shown and no similarity metric is given. An elemental ratio or a
   bond census measured on a 7 % dialysis retentate is not the same object as a model's lumped
   melanoidin pool, and this is the sharpest instance of that rule in the corpus.
3. **This paper and Mundt 2004 disagree on how much glycine is decarboxylated, and both are on
   disk.** Fang: **27 ± 4 % (dry) / 33 ± 4 % (solution) of glycine reactant** undergoes Strecker
   decarboxylation. Mundt & Wedzicha: **~2/3 of the glycine incorporated into the polymer** is
   decarboxylated (0.662 of 0.951 mol per mol glucose). The two are not stated on the same base
   — Fang's denominator is all glycine reactant, Mundt's is incorporated glycine — but the gap
   survives any reasonable re-basing. Fang's own introduction says the literature spread on this
   quantity runs "**between 0 and 65 %**" (his refs. 21 = Wedzicha & Kaputo 1992, the Leeds group
   again, and 24 = Feather & Huang 1985), and his discussion argues that earlier NMR work
   "significantly overestimated" Strecker degradation by attributing all fragmentation,
   deamination included, to it. **Do not average them and do not pick one silently.** The
   conditions differ enormously (Mundt: 70 C, pH 5.5, days, aqueous; Fang dry: 125 C, 2 h,
   solid), and Mundt's own cited direction — more amino-acid incorporation at lower temperature —
   is consistent with more decarboxylated glycine at 70 C than at 125 C.
4. **The paper's own internal cross-check on Strecker degradation does not close, and the authors
   say so.** Two routes give 16 % and 11 % for "lost C1" in the dry reaction; a third gives 21 %;
   the adopted value is **27 ± 4 %**, chosen on the argument that "21 % underestimates the
   extent". The separate C2-bonding cross-check disagrees by 42 vs 36 % (dry) and 35 vs 30 %
   (solution), and the paper states outright that "the origin of the moderate discrepancy between
   the two calculations is not clear". Carry the ± 4 % and carry this note with it.
5. **The pyrazine null is about the polymer, not about the pot.** Fang finds no significant 15N
   signal for pyrazines or pyridines (250-350 ppm) or imines (300-370 ppm), and no 13C at 146 or
   ~150 ppm, in either melanoidin. That is a strong statement about the **high-molecular-weight
   fraction**, which is all he measured. `data/keys/compounds.yml` keys eight pyrazines as
   volatile products, and they are found in real systems by every headspace paper in the corpus.
   **The correct reading is that pyrazine nitrogen and `MEL_N` are different sinks**, not that
   pyrazines do not form.
6. **Table 1's column label "nonprotein Nt" is almost certainly a typesetting slip for
   "nonprotonated Nt".** Every discussion of that column in the text is about protonation, the
   column pairs with "N−H (%)", and 21 + 79 = 100 and 15 + 85 = 100. Transcribed as printed
   above; read as nonprotonated.
7. **The nitrogen percentages are "approximate" by the table's own title** and rest on peak
   assignments the paper defers: "these assignments will be justified in a future publication,
   based on 15N-13C-13C NMR experiments". The imidazolium (~170 ppm) and oxazolium (~200 ppm)
   assignments in particular are forward-referenced, not established here.
8. **No nitrogen mass balance anywhere.** Every nitrogen number is a *share of the nitrogen that
   is in the polymer*. The paper never measures how much glycine nitrogen left as ammonia, and
   ammonia release is explicitly named in the introduction as a Strecker product that "may be
   liberated". So this paper cannot tell the trunk how much nitrogen `MEL_N` should hold, only
   what the nitrogen that is there looks like.
9. **The Supporting Information is not on disk.** Figures S1-S5 carry the JCC pulse sequence, the
   L-leucine-13C1,2 calibration, the 13C{15N{1H}}-HSQC-REDOR sequence, the two model-compound
   calibrations (75 % and 8 % residual signals) and the SUPER powder patterns for glycine C1.
   **Fetch it** — the 1/0.75 scaling factor that every nonprotonated-N number depends on is
   calibrated there.
10. **Two spectrometer-side caveats the paper raises itself.** (i) 15N cross-polarisation
    under-represents nonprotonated nitrogen, which is why the 15N-detected Nt:NH of 75:25 is
    called "(underestimated)" and the 13C-detected 79:21 / 92:8 pair is preferred. (ii) The
    13C{15N} REDOR experiments "did not use direct polarization of 13C", so "the relative peak
    intensities are not fully quantitative" and had to be renormalised against the quantitative
    DP integrals of Figure 2c,f. Only the DP numbers are quantitative on their own.
11. **What this paper does not contain**: any rate constant; any activation energy; any time
    course or intermediate time point; any temperature series; any pH series at constant
    everything-else; any combustion elemental analysis or C/N; any absorbance or extinction
    coefficient; any molecular-weight distribution; any nitrogen mass balance; any water-activity
    measurement; any amine other than glycine; any sugar other than glucose; any real food.
12. **What to request from the authors**: (i) the Supporting Information PDF; (ii) 13C and 15N
    spectra and integrals of the **insoluble** dry-reaction fraction, which is 50 % of the product
    and is currently covered only by the word "similar"; (iii) a combustion CHN analysis of the
    same two samples, which would convert the implied C/N of section 3 arithmetic 3 from a
    derived number into a measurement and would settle the disagreement with Mundt 2004
    directly; (iv) whether any glycine nitrogen was recovered as ammonia.
13. **Registry gaps against `data/keys/compounds.yml`**: glucose, glycine, carbon dioxide,
    ammonia and any melanoidin pool are all unkeyed among the 75 ids. The trunk's `Glc`, `Gly`,
    `MEL_C`, `MEL_N` and `FRAG_C` are network-local names, not registry ids. A bond-census
    benchmark built from this paper has no registry surface at all and would need one invented —
    which is a decision, not an omission, given that `species.py` already declares why a polymer
    pool has no molecular weight.
