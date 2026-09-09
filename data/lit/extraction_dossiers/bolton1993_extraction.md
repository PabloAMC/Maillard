# Bolton, Reineccius, Liardon & Huynh-Ba 1993/1994 — EXTRACTION (a bouillon-type solids blend — MSG 45.19, NaCl 45.19, IMP 5.73, D-glucose 2.16, thiamine-HCl 1.08, cysteine-HCl 0.43, D-xylose 0.22 g/100 g — at 30 % solids, pH 5.5-5.8, a_w 0.83, sealed 125 mL vials, 120 C for 1 h; MFT by dichloromethane extraction and GC/MS-SIM against dodecane, with 34S-labelled cysteine to split the sulfur's origin)

### THE ONE ISOTOPE EXPERIMENT IN THE CORPUS THAT SPLITS MFT's SULFUR BETWEEN THIAMINE AND CYSTEINE IN A REAL SAVOURY MATRIX, AND IT COMES OUT 92 : 8 FOR THIAMINE — with the harder finding alongside it that removing thiamine gave **no detectable MFT at all**, so in this pot the whole pentose-plus-sulfide lane the engine spends most of its constants on contributes under about two percent.

**Source on disk:** `data/articles/bolton1993.pdf` (9 pp., ACS Symposium Series 543, *Thermally
Generated Flavors: Maillard, Microwave, and Extrusion Processes*, ed. Parliment, Morello & McGorrin,
chapter 22, pp. 270-278). Read from the `pdftotext -layout` text layer
(`scratchpad/articles/bolton1993.txt`, 640 lines). The running text, Table I, the Materials and
Methods, the Conclusions and all 25 references came through legible; the scanned page carries a
"Downloaded by KTH ROYAL INST OF TECHNOLOGY" watermark down the left margin of every page, which the
extractor interleaves into the text but which does not obscure it. **Table II's individual replicate
values were partly garbled** by the extractor ("4J&", "H6S", "JSQ"), so journal page 275 was rendered
at 300 dpi with `pdftoppm` and Table II is re-typed in section 3 from that render, with every mean and
standard deviation re-computed from the replicates as a check (they all reproduce). Figures 1 and 2
(the two proposed MFT formation schemes) are chemical-structure drawings and are figure-only; their
structures are additionally mangled in the text layer and were not relied on. There is no
supplementary material — this is a 1993 book chapter.

## 0. Identity

| field | value |
|---|---|
| Title | "Role of Cysteine in the Formation of 2-Methyl-3-furanthiol in a Thiamine-Cysteine Model System" |
| Authors | T. A. Bolton and G. A. Reineccius — Department of Food Science and Nutrition, University of Minnesota, St. Paul, MN 55108; **R. Liardon and T. Huynh-Ba** — Nestlé Research Center, Lausanne, Switzerland |
| Venue | Chapter 22 in *Thermally Generated Flavors*, ACS Symposium Series **543**, American Chemical Society, Washington DC, pp. **270-278** |
| **DOI** | **10.1021/bk-1994-0543.ch022** (printed in the left margin of every page) |
| **Is this the same work as the repository's "Bolton 1994"?** | ★ **YES — it is one and the same chapter, not a companion.** Four independent confirmations on the page itself: (i) the printed DOI `10.1021/bk-1994-0543.ch022` is **exactly** the `source_doi` in `data/benchmarks/thiamine_cys_glucose_120C_Bolton1994.json`; (ii) the page range 270-278 matches that file's citation, "ACS Symp. Ser. 1994, 543 (Thermally Generated Flavors), 270-278"; (iii) the author list matches to the fourth author; (iv) the numbers match — Table I's solids composition and Table II's 389 ± 57 / 489 ± 90 / not detected / 425 ± 103 are precisely what the benchmark's `content_verification.quoted_values` records, and that block already states it was "read from data/articles/bolton1993.pdf on 2026-09-04". **There is no second Bolton paper.** |
| Why the two years | the chapter's own margin reads **"Publication Date: November 30, 1993"** and the article was **"RECEIVED February 9, 1993"**, while the copyright line reads **"0097-6156/94/0543-0270$06.00/0 © 1994 American Chemical Society"** and the DOI carries `bk-1994`. The file stem uses the **publication date** (1993); the benchmark and the repository's hold-out label use the **volume year** (1994). Both are correct citations of the same object; **the repository should keep one convention and say which.** |
| Paper type | **Isotope-labelling / mechanism attribution**, with quantification of one compound. No kinetics: one temperature, one time, one pH, no time course. |
| Predecessor | ref 21 = Bolton, Reineccius & Liardon, University of Minnesota, **unpublished data, 1989-1990** — the source of the claim that "MFT and the characteristic meat-like aromas could not be produced without thiamine" in this model system, and of the earlier BMFD observation. **Not available.** |
| Label provenance | the 34S-cysteine-HCl monohydrate was **synthesised at Nestlé Research Center, Lausanne**, at **100 % 34S**, per ref 22 = Huynh-Ba & Fay, *J. Labelled Compd. Radiopharm.* 1990, 28, 1185 |
| Companions on disk | `cerny2008_extraction.md` (the HMP identification — this chapter's Figure 1 is built entirely on HMP as the thiamine intermediate), `cerny2007_extraction.md` (the other thiamine-vs-sugar split in the corpus), `mottram1995_extraction.md` (same HMP claim, cited to the same 1979/1984/1993 sources), `hofmann1998b_extraction.md`, `schutte1972_extraction.md`, `whitfield1988_extraction.md` |

## 1. Why it matters

This chapter bears on **three** named parts of the sulfur lane at once, and on the benchmark panel.

**(a) It is the direct measurement of the branch `PROHIBITED_DERIVATIONS` forbids anyone to fix.**
`parameters_sulfur.py` carries, as a prohibited derivation, "any fixed thiamine:sugar MFT branch
fraction", enforced by requiring every split to be "a ratio of mass-action fluxes computed at run
time". The evidence cited there is Cerny 2007 Table 5, where a 2× change in precursor loading moves
the xylose share of MFT from 15 % to 46 %. **This chapter is a second, independent observation of the
same branch**, by isotope rather than by inference, in a different matrix: with 100 % 34S-cysteine
present alongside thiamine, **only 7.5-8.0 % of the MFT carried cysteine's sulfur** (12.50 ± 0.76 %
and 12.11 ± 1.36 % measured against a 4.58 ± 0.09 % natural-abundance baseline). The two studies
disagree by a factor of two to six on where the branch sits, which is exactly the prohibition's point:
the split moves with the pot, and no constant may encode it. **This dossier adds the second data
point that makes that prohibition an observation rather than a policy.**

**(b) It puts a hard ceiling on the entire pentose-plus-sulfide lane — in this matrix.** Model system
**III** contained 34S-cysteine, D-xylose, D-glucose and IMP (a ribose source) and **no thiamine**, and
produced **no detectable MFT** in any of three replicates, against ~440 ng in the thiamine-containing
systems, with a stated detection limit of **ca. 10 ng**. So in this pot the sugar-plus-sulfide route to
MFT is **≤ ~2 % of the thiamine route (mine, 10/440)**. The engine's sugar-side MFT edges —
`r_ddp_mft_hs`, `r_nf_mft`, `r_ha_mp_mft`, together with `k_mgo_mp` and `k_glc_ha` feeding them — are
the most heavily parameterised part of the network, and here they produce nothing measurable. The
conditional must travel with the finding (Flags 5): a_w 0.83, a near-saturated brine, 120 °C, one
hour, pH 5.5-5.8. It is a bound on the lane **in a low-water savoury matrix**, not a refutation of the
lane in dilute water where Hofmann measures it.

**(c) It says the thiamine route in the engine is missing an edge.** The authors' own conclusion is
that MFT here forms "primarily ... by direct cyclization of 5-hydroxy-3-mercaptopentan-2-one resulting
from thiamine degradation" — which is exactly `r_thi_hmp` followed by `r_hmp_mft`, the engine's
thiamine lane, whose two constants `k_thi_hmp` and `k_hmp_mft` are the only two steps in the
`parameters_sulfur.py` step table carrying an **empty note string**. But the labelled 7.5-8.0 % is
attributed by the authors to the *second* thiamine pathway of their Figure 1, in which HMP loses its
own -SH, becomes 3,5-dihydroxypentan-2-one, cyclises to an oxofuran, and then **takes up exogenous
H2S**: "the direct contribution of cysteine to MFT must occur by the addition of 34S-H2S to the
oxofuran intermediate." **The engine has no such edge.** `r_hmp_mft` is one first-order step in which
HMP keeps its own sulfur; there is no HMP -> oxofuran -> +H2S -> MFT branch and no oxofuran species.
So about 8 % of the MFT in this pot travels a route with thiamine's carbon and cysteine's sulfur that
the network cannot express, and the engine would have to attribute all of it to the sugar lane.

**(d) It is a live benchmark, and this dossier checks it.** `thiamine_cys_glucose_120C_Bolton1994.json`
targets **11.7 ppb MFT** (Table II system I, 389 ng / 33.3 g) with a factor-3 contract, and
`docs/guides/INTRODUCTION.md` lists it among the five laboratories of the panel. **Every printed value
in that bundle reproduces against the chapter**, including the conditions block. Two small molar-mass
discrepancies in the derived concentrations are recorded in Flags 3, and one claim in the bundle's
`concentration_note` — that MSG, NaCl and IMP "are not reactive precursors on this path" — is too
strong and is discussed in Flags 4.

What this chapter does NOT give: any rate, any barrier, any time course, any second temperature, any
pH variation, any measurement of thiamine or cysteine loss, any other volatile in numbers, any yield
of MFT on thiamine, and any usable measurement of the disulfide.

## 2. Methods as they matter to a model

- **The four model systems.** All four share MSG, NaCl, IMP, D-glucose; systems I, II and IV carry
  thiamine-HCl; system I carries **unlabelled** cysteine-HCl and systems II, III, IV carry **100 %
  34S-labelled** cysteine-HCl; systems III and IV carry D-xylose. In one line: **I** = thiamine +
  unlabelled Cys, no xylose (the control that fixes the natural-abundance baseline in situ); **II** =
  thiamine + 34S-Cys, no xylose; **III** = **no thiamine**, 34S-Cys + xylose; **IV** = thiamine +
  34S-Cys + xylose. Exact composition in section 3, Table I.
- **Reagents.** Monosodium glutamate (food grade, Takeda), NaCl (food grade, Morton), inosine-5'-
  monophosphate (food grade, Takeda), **D-glucose monohydrate** (USP, Sigma), thiamine-HCl (USP,
  Sigma), **cysteine-HCl monohydrate** (food grade, Takeda) **or 34S-labelled cysteine-HCl
  monohydrate** (100 % 34S, Nestlé Lausanne), D-xylose (food grade, Takeda). **Both glucose and
  cysteine are charged as monohydrates** — this matters for the molarities (Flags 3).
- **The pot.** "Samples were prepared **in triplicate** using **10 g of the model system and 23.3 g of
  distilled water** (at 30 % solids level complete solubility was observed) and were thermally
  processed at **120 °C for one hour** in **125 mL glass vials sealed with Teflon-lined septa and
  aluminum crimp caps**." So **33.3 g of fill in a 125 mL vial**, i.e. roughly **92 mL of air
  headspace**, sealed, **no inert-gas purge is described** — the oxygen inventory is whatever the air
  holds. Stirring is not mentioned.
- **pH and water activity — both printed, and both unusual.** "The pH of the model system ranged from
  **5.5 to 5.8 before and after processing** and the initial a_w of the system was **0.83** at ambient
  temperature." The pH stability across a 120 °C hour is itself a result: this pot is buffered by its
  own gram-scale MSG and IMP, unlike almost every dilute model system in the corpus. **a_w 0.83 is
  the lowest water activity of any MFT measurement in this dossier set**, and the engine has no a_w
  term (the benchmark bundle says so explicitly).
- **Work-up.** "Post processing, samples were **immediately diluted to 10 % solids** with distilled
  water", then extracted **5 × 20 mL with dichloromethane** in a 250 mL separatory funnel;
  **dodecane** added as internal standard **at 1 ppm to the 100 mL solvent extract**; dried over
  anhydrous MgSO4; concentrated under N2 to **0.10 mL**. **Extraction, not headspace** — which is why
  the benchmark classes it `internal_standard_gcms_sim`.
- **GC/MS-SIM.** Splitless injection; HP 5890 GC with an HP 5970 low-resolution mass selective
  detector in **selected ion monitoring**; oven 35 °C for 2 min, **3 °C/min** to 250 °C, held 20 min;
  helium at 15 psi head pressure; **DB-5**, 30 m × 0.32 mm i.d. × 1.0 µm film. MFT identified against
  a **pure standard** by mass spectrum and retention time. Quantified by internal standard: "The
  relative response factor of MFT/C12 was calculated from the area sum of **three integrated ions**
  from the SIM spectra of calibration solutions of dodecane (**m/z 58, 85, 170**) and MFT (**m/z 71,
  113, 114**)."
- **The isotope measurement.** The 34S share is reported as **m/z(116/114) × 100** (footnote a to
  Table II) — the ratio of the 34S-bearing molecular ion to the 32S one. **The natural-abundance
  baseline was measured on a pure MFT standard and printed: 4.5 %**, and system I's in-situ value is
  4.58 ± 0.09 %, so the baseline is confirmed twice. The "net" incorporation is the labelled system's
  value minus this baseline.
- **Replication.** **Genuine triplicate, with every individual replicate printed** — the only paper in
  this five-dossier set that prints its raw replicates. Means and SDs re-computed in section 3 and all
  six reproduce exactly.
- **Detection limit.** Given in the Results for the no-thiamine case as **ca. 10 ng** (per 33.3 g
  sample), i.e. **~0.3 ppb**.
- **What was looked for and not reported.** Bis-(2-methyl-3-furyl) disulfide (BMFD — the engine's
  `MFTD`): "only **trace amounts** of BMFD were detected in model systems I, II and IV ... Only partial
  spectra of BMFD were obtained by GC/MS at the appropriate retention time and, as a result, **data for
  BMFD are inconclusive and not included in this study**." (Flags 6.)

## 3. Tables re-typed

Both tables are re-typed in full. Table I came through the text layer cleanly; Table II was read from
a 300 dpi render of journal page 275 because the extractor garbled three of its twelve replicate
values.

### Table I. "Solids Composition of Model Systems (g/100 g)"

| Component | I | II | III | IV |
|---|---:|---:|---:|---:|
| MSG | 45.19 | 45.19 | 45.19 | 45.19 |
| NaCl | 45.19 | 45.19 | 45.19 | 45.19 |
| IMP | 5.73 | 5.73 | 5.73 | 5.73 |
| D-Glucose | 2.16 | 2.16 | 2.16 | 2.16 |
| Thiamine-HCl | 1.08 | 1.08 | **–** | 1.08 |
| Cysteine-HCl | 0.43 | – | – | – |
| **34S** Cysteine-HCl | – | 0.43 | 0.43 | 0.43 |
| D-Xylose | – | – | 0.22 | 0.22 |

(The column heads are printed as "/ II Ill rv" in the scan's italic face; they are systems I, II, III,
IV. Columns I and II sum to 99.78 g/100 g and columns III and IV to 98.92 and 100.00 respectively —
see Flags 2.)

### Table II. "Percent 34S Incorporated into MFT Formed from Cysteine and 34S Labeled Cysteine"

Footnote a: `m/z(116/114) × 100`. Footnote b: "ng MFT/g **based on 33.3 g of thermally-processed
sample**" — the column head reads "ng MFT/g" but the quantity is **ng per 33.3 g sample**, not per
gram (Flags 1).

| Sample | Replicate | Percent 34S ᵃ | ng MFT ᵇ |
|---|---|---:|---:|
| **I.** unlabeled cysteine, no xylose | 1 | 4.54 | 330 |
| | 2 | 4.52 | 444 |
| | 3 | 4.68 | 393 |
| | **mean** | **4.58 ± 0.09** | **389 ± 57** |
| **II.** with 34S cysteine, no xylose | 1 | 11.66 | 405 |
| | 2 | 13.11 | 477 |
| | 3 | 12.75 | 585 |
| | **mean** | **12.50 ± 0.76** | **489 ± 90** |
| **III.** with 34S cysteine, **no thiamine** | 1 | **not detected** | – |
| | 2 | – | – |
| | 3 | – | – |
| **IV.** with 34S cysteine (+ xylose) | 1 | 11.48 | 306 |
| | 2 | 11.19 | 489 |
| | 3 | 13.68 | 480 |
| | **mean** | **12.11 ± 1.36** | **425 ± 103** |

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| natural 34S/32S abundance (m/z 116/114) in a **pure MFT standard** | **4.5 %** | Results |
| net labelled MFT | "A net **7.5 to 8.0 %** of the MFT formed contained labeled sulfur (model systems II and IV vs I)" | Results and Abstract |
| the same expressed as a mass | "approximately **29 to 39 ng** per 33.3 g of processed model system" | Results (see the arithmetic below — my recomputation gives 32-39 ng) |
| MFT detection limit | "below detectable limits (**ca. 10 ng**)" | Results |
| the comparison the no-thiamine result is set against | "**ca. 440 ng** formed in model systems I, II, and IV" | Results |
| pH | **5.5 to 5.8 before and after processing** | Process |
| water activity | **0.83** at ambient temperature, initial | Process |
| fill and vessel | 10 g solids + 23.3 g distilled water in a **125 mL** vial | Process |
| temperature and time | **120 °C, one hour** | Process |
| extraction | 5 × 20 mL CH2Cl2; dodecane at **1 ppm** in the 100 mL extract; concentrated to **0.10 mL** | Analysis |
| SIM ions | dodecane m/z **58, 85, 170**; MFT m/z **71, 113, 114** | Analysis |
| BMFD | "only **trace amounts** ... data for BMFD are **inconclusive and not included**" | Results |
| effect of added xylose | "**added xylose did not influence** the amount of MFT formed in model systems III and IV" | Results |
| Hincelin's contrasting result (**cited, ref 24, not measured here**) | thiamine + added xylose at pH 7.0, 140 °C increased MFT **4 to 5 times** over thiamine alone | Results |
| Grosch & Zeiler-Hilgart's method effect (**cited, ref 13**) | a **4-5 fold** increase in recovered MFT by Likens-Nickerson distillation-extraction over direct extraction | Results |
| the pH trend in the older literature (**cited, refs 15, 16, 18, 25**) | MFT "decreased in concentration with increasing pH in the range of **2.3 to 9.5**" | Results |

### Arithmetic on the printed numbers (all mine)

**1. Table II is internally consistent — every mean and SD reproduces.** System I: (330+444+393)/3 =
**389.0**, sample SD **57.1**. System II: (405+477+585)/3 = **489.0**, SD **90.6**. System IV:
(306+489+480)/3 = **425.0**, SD **103.2**. Percent 34S: I (4.54+4.52+4.68)/3 = **4.580**, SD 0.087;
II **12.507**, SD 0.755; IV **12.117**, SD 1.362. All six agree with the printed values to the last
digit shown. **This is the best-documented table in the five-paper set.**

**2. The net incorporation, and a small discrepancy in the printed mass.** Net = labelled − baseline:
system II 12.50 − 4.58 = **7.92 %**; system IV 12.11 − 4.58 = **7.53 %**. That reproduces the
abstract's "7.5 to 8.0 %" exactly. Converting to a mass **within each system**: 7.92 % × 489 =
**38.7 ng** (II) and 7.53 % × 425 = **32.0 ng** (IV), i.e. **32-39 ng**, against the text's
"approximately 29 to 39 ng". The lower bound 29 is reproduced only by applying system IV's percentage
to system I's MFT (7.53 % × 389 = 29.3), which crosses systems. The difference is immaterial to any
conclusion; it is recorded because the house rule is that printed arithmetic gets checked.

**3. The concentration basis, and the ppb conversion.** 389 ng in 33.3 g = **11.7 ng/g = 11.7 ppb**,
which is the benchmark's target. The other two systems are 489/33.3 = **14.7 ppb** and 425/33.3 =
**12.8 ppb**. The detection limit, 10 ng, is **0.30 ppb**.

**4. Molar charges per vial (mine), and the basis question.** 10 g of solids gives, per vial:
MSG 4.519 g, NaCl 4.519 g, IMP 0.573 g, glucose monohydrate 0.216 g, thiamine-HCl 0.108 g,
cysteine-HCl monohydrate 0.043 g, xylose 0.022 g. Using the **stated hydrate forms**: thiamine-HCl
(M 337.27) **0.320 mmol**; cysteine-HCl monohydrate (M 175.63) **0.245 mmol**; glucose monohydrate
(M 198.17) **1.090 mmol**; xylose (M 150.13) **0.147 mmol**; NaCl **77.3 mmol**; MSG monohydrate
(M 187.13) **24.1 mmol**. Expressed **per litre of added water** (23.3 mL), the basis the benchmark
uses: thiamine **13.7 mM**, cysteine **10.5 mM**, glucose **46.8 mM**, xylose **6.3 mM**, and — the
numbers nobody has written down — **NaCl 3.3 mol/L and MSG 1.03 mol/L**. Expressed per 33.3 g of
total fill instead: thiamine 9.6 mM, cysteine 7.4 mM, glucose 32.7 mM. **The two bases differ by
1.43×**, and the benchmark's choice (added water) is the higher one; whichever is used should be
stated, because it is a 43 % swing on every precursor.

**5. The salt explains the water activity (mine).** ~3.3 mol/kg NaCl plus ~1.0 mol/kg MSG in the
added water is a near-saturated brine, and a_w 0.83 is what that predicts. **The internal consistency
is a check that Table I, the 30 % solids figure and the printed a_w all describe the same pot** — and
it makes plain that this is a bouillon, not a model solution.

**6. Yield of MFT on thiamine (mine).** 389 ng / 114.17 g mol⁻¹ = **3.41 nmol** against 0.320 mmol of
thiamine charged = **0.0011 mol %**, i.e. about **11 parts per million of the thiamine** ends as MFT.
Even at the highest system (489 ng) it is 0.0013 mol %. **The MFT channel is a trace branch off
thiamine**, which is the same order of magnitude as the trace-branch finding in the TTCA pot (see
`zhai2023b_extraction.md` section 3) and is worth knowing before anyone expects a thiamine lane to
carry appreciable flux.

**7. The ceiling on the sugar route (mine).** ≤ 10 ng (system III's limit) against ~440 ng (the
authors' own round figure for I, II, IV) is **≤ 2.3 %**. Against system I specifically, ≤ 10/389 =
**≤ 2.6 %**. Either way: **in this matrix the sugar-plus-cysteine lane to MFT is under about 2.5 % of
the thiamine lane.**

**8. Added xylose does nothing here, and the null is well-powered enough to say so (mine).** System II
(34S-Cys, no xylose) 489 ± 90 against system IV (34S-Cys, + xylose) 425 ± 103. The difference is
64 ng with a pooled SD of about 97 on n = 3 each; a two-sample t gives t ≈ 0.81, nowhere near
significance. The paper's own wording — "added xylose did not influence the amount of MFT formed" — is
supported. **The 34S share is likewise unmoved: 12.50 ± 0.76 against 12.11 ± 1.36.** Xylose at 6.3 mM
adds neither MFT nor labelled sulfur. This is a **direct contradiction of Hincelin et al. 1992** as
the paper cites it (thiamine + xylose, pH 7.0, 140 °C, 4-5× more MFT), and the paper names the likely
reasons: different pH, different temperature, different matrix.

## 4. Kinetic numbers the repository can use

**There is no rate constant and no activation energy in this chapter** — one temperature, one time,
one pH. What it supplies is a **level**, a **branch fraction by isotope**, and two **negatives**, all
in a matrix the engine has never been fitted on.

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Keyed: `2_methyl_3_furanthiol`,
`bis_2_methyl_3_furyl_disulfide`, `hydrogen_sulfide`, `thiamine_availability`. **Not keyed:**
thiamine as a reactant, cysteine, D-glucose, D-xylose, MSG/glutamate, IMP, ribose,
5-hydroxy-3-mercaptopentan-2-one (the engine's `HMP`), 3,5-dihydroxypentan-2-one,
2-methyl-3-oxotetrahydrofuran, 2-methyl-4,5-dihydro-3-furanthiol, or dodecane.

Every row below shares: the Table I solids blend at **30 % solids** (10 g solids + 23.3 g distilled
water), **pH 5.5-5.8 unbuffered but self-buffered, a_w 0.83**, sealed **125 mL** vial with ~92 mL of
**air** headspace, **120 °C for 60 min**, **triplicate**, diluted to 10 % solids post-process,
dichloromethane extraction, **GC/MS-SIM against dodecane at 1 ppm** with a measured relative response
factor.

| step / observable | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| **MFT level, thiamine + unlabelled Cys (the benchmark row)** | MFT | **389 ± 57** ng per 33.3 g = **11.7 ± 1.7 ppb** | ng / ppb (w/w) | system I | **nothing is fitted in this chapter** | Table II p. 275 | **measured_rate**-class *level*: strictly a **level_only** with a real internal standard and a measured response factor — the strongest quantification class in this five-paper set |
| MFT level, thiamine + 34S-Cys | MFT | **489 ± 90** ng = **14.7 ppb** | as above | system II | — | Table II | **level_only** |
| MFT level, thiamine + 34S-Cys + xylose | MFT | **425 ± 103** ng = **12.8 ppb** | as above | system IV | — | Table II | **level_only** |
| **MFT without thiamine** | MFT | **not detected**, all three replicates; limit **ca. 10 ng = 0.30 ppb** | ng | system III (34S-Cys + xylose + glucose + IMP, **no thiamine**) | — | Table II + Results | **threshold** — a one-sided bound, and the most consequential single result in the chapter |
| **the sugar-plus-sulfide lane's ceiling in this matrix** | MFT(no thiamine) / MFT(thiamine) | **≤ 2.3-2.6 %** | — | as above | — | derived from the printed ≤10 ng and 389-489 ng (mine) | **derived_assumption** built on a `threshold` and a `level_only` |
| natural-abundance baseline, measured on a pure standard | m/z 116/114 | **4.5** | % | pure MFT standard | — | Results | **level_only** |
| the same baseline measured **in situ** | m/z 116/114, system I | **4.58 ± 0.09** | % | system I | — | Table II | **level_only** — an in-pot confirmation of the baseline, which is what makes the net incorporation trustworthy |
| **34S share of MFT, labelled systems** | m/z 116/114 | **12.50 ± 0.76** (II) and **12.11 ± 1.36** (IV) | % | systems II, IV | — | Table II | **measured** isotope ratio; class **within_study_ratio** |
| **★ the cysteine-sulfur share of MFT (the branch fraction)** | net 34S | **7.92 %** (II) and **7.53 %** (IV); the paper's "**7.5 to 8.0 %**" | % of the MFT formed | 120 °C, 60 min, pH 5.5-5.8, a_w 0.83, thiamine present | — | Table II minus system I's baseline (mine; the range is printed) | **within_study_ratio** — ★ **a direct isotopic measurement of the thiamine : exogenous-sulfur split, 92 : 8** |
| the same as a mass | labelled MFT | 32-39 ng (mine); the text prints "approximately 29 to 39 ng" | ng per 33.3 g | as above | — | Results (mine, section 3 arithmetic 2) | **derived_assumption** |
| **effect of added xylose** | MFT(+xylose) / MFT(no xylose), both 34S | **425/489 = 0.87**, not significant (t ≈ 0.81, n = 3, mine) | — | systems IV vs II | — | Table II (mine) | **within_study_ratio** — a **null**, and a well-documented one |
| effect of added xylose on the sulfur split | 12.11 ± 1.36 vs 12.50 ± 0.76 | unchanged | % | as above | — | Table II | **within_study_ratio** (a second null) |
| MFT yield on charged thiamine | MFT / thiamine | ~**0.0011** | mol % | system I | — | derived (mine, section 3 arithmetic 6) | **derived_assumption** |
| **BMFD (the engine's `MFTD`)** | disulfide | "trace amounts" only; **spectra partial, data declared inconclusive and excluded** | — | systems I, II, IV | — | Results | **peak_area_only**, and by the authors' own statement **not usable** (Flags 6) |
| pH stability across the cook | 5.5-5.8 before **and after** | — | — | as above | — | Process | **level_only** — a genuinely useful condition datum |
| water activity | **0.83** | — | initial, ambient | as above | — | Process | **level_only** |
| Figures 1 and 2 (the two proposed MFT schemes) | — | — | — | — | Figs. 1, 2 | **figure_only**, and additionally **adapted from refs 9, 12, 18 and 20** — not results of this work |

### Can these be put on the same basis as the sulfur lane's constants?

**(a) The level — already is, as a benchmark, and this dossier verifies it.** 11.7 ppb at 120 °C /
60 min with the conditions block as printed. Everything the bundle quotes is confirmed. The two
things a user of that bundle should know are in Flags 3 (two precursor molarities are computed from
anhydrous rather than the stated monohydrate masses, making them 7-11 % high) and Flags 4 (the pot's
dominant amine is glutamate at ~1 mol/L, twenty times the glucose, and it is not in the network).

**(b) The 8 % branch fraction — as a scored ratio, never as a parameter.** This is precisely the
quantity `PROHIBITED_DERIVATIONS` refuses to let anyone hard-code, and precisely the quantity a
mass-action network should be able to **predict**. It is an excellent panel row: dimensionless,
isotopically measured, replicated three times, with the natural-abundance baseline measured twice.
**Its condition set must travel with it** — at a_w 0.83 and 3.3 M NaCl, and with the caveat in (d).

**(c) The two negatives — as one-sided rows.** "No MFT without thiamine at a 10 ng limit" and "added
xylose changes nothing" are the kind of falsifying observation the lane most needs, because the
engine's sugar-side MFT edges are its most heavily parameterised and its least directly constrained.
A model that predicts appreciable MFT for system III is wrong **in this matrix**, and a model whose
MFT jumps when 6.3 mM xylose is added is wrong here too.

**(d) What blocks a clean reading of (b) and (c): the engine cannot express the labelled route.** The
authors attribute the 7.5-8.0 % to thiamine's *second* pathway — HMP losing its sulfur, cyclising to
an oxofuran, and picking up cysteine's H2S. The engine has only `r_hmp_mft`, a single step that keeps
HMP's own sulfur. So an engine scored against the 8 % would have to produce it from the **sugar**
lane, which system III says produces nothing. **The row is scoreable only if the scoring is on the
isotope share, and it will then be measuring an edge the network does not have.** That is worth
saying before the row is installed, not after.

**(e) What cannot be transported at all.** Nothing here is a rate. The single temperature makes a
barrier impossible, and pairing this 120 °C level with any other laboratory's level to extract one
would repeat exactly the two-point, two-lab, two-matrix error `PROHIBITED_DERIVATIONS` already forbids
twice.

## 5. Flags

1. **Table II's column header is wrong, by a factor of 33.3.** The head reads "**ng MFT/g**"; footnote
   b reads "ng MFT/g **based on 33.3 g of thermally-processed sample**"; the Abstract reads
   "**389 to 489 ng/33.3 g model system**". The quantity is **ng per sample**, and dividing by 33.3
   gives the 11.7-14.7 ppb the benchmark carries. **A reader taking the header literally would be
   33.3× high** — i.e. would read 389 ppb instead of 11.7 ppb. The benchmark bundle already handles
   this correctly; the trap is in the chapter, and it should be recorded.
2. **Table I's columns do not sum to 100.** Systems I and II sum to **99.78** g/100 g, system III to
   **98.92**, system IV to **100.00**. The 100.00 for system IV means the composition is a normalised
   recipe rather than a weighed one, and that the missing 0.22-1.08 g in the other columns is simply
   the omitted component with nothing put back in its place. **The consequence is that systems I-III
   are not at the same total solids as IV**, by up to 1.1 %, which is negligible against a ±13-24 %
   MFT spread but means the four systems are not exactly matched.
3. **Two of the benchmark's precursor molarities are computed from the wrong salt form.** The Materials
   section states **"D-glucose monohydrate (USP, Sigma)"** and **"cysteine-HCl monohydrate"**. Using
   the monohydrate molar masses (198.17 and 175.63) gives, per litre of added water, **glucose
   46.8 mM** and **cysteine 10.5 mM**; the bundle carries **51.5 mM** and **11.7 mM**, which are the
   anhydrous values (180.16 and 157.62). The bundle is therefore **10 % high on glucose and 11 % high
   on cysteine**. Thiamine-HCl at 13.7 mM is correct. This is well inside the bundle's factor-3
   contract and changes no conclusion, but it should be corrected rather than carried.
4. **The bundle's claim that MSG, NaCl and IMP "are not reactive precursors on this path" is too
   strong, and the omission is large.** By my arithmetic the pot contains **~1.03 mol/L
   monosodium glutamate** and **~3.3 mol/L NaCl** in the added water. Glutamate is a free amino acid
   at **twenty times the glucose concentration** and **eighty times the cysteine**; it is a Maillard
   reactant (glycosylamine formation, Strecker degradation) and it will compete for every carbonyl in
   the pot. IMP at ~0.57 g per vial is a ribose source, and the paper's own introduction cites Zhang &
   Ho's IMP-plus-cysteine route to MFT as one of the mechanisms being tested. **The narrow claim that
   survives is the one the isotope data actually support** — that in *this* pot the MFT sulfur is 92 %
   thiamine's and the no-thiamine system made none — **not** the general claim that these three
   components are inert. A benchmark that omits a 1 mol/L amine from a Maillard pot is modelling a
   different pot; the difference should be declared rather than dismissed.
5. **Every finding here is conditional on a low-water, high-salt matrix, and the paper says so
   itself.** a_w **0.83**, near-saturated brine, pH 5.5-5.8, 120 °C, 60 min, direct solvent
   extraction. The authors list the reasons their numbers differ from others': "differences in
   precursor concentrations, pH, processing times and temperatures, and analytical methodologies",
   and specifically note that the cysteine-ribose systems that *do* make MFT used higher precursor
   concentrations, phosphate buffer and 130-140 °C. **The "sugar route contributes nothing"
   conclusion must never be quoted without a_w 0.83 attached**, and the engine has no a_w term with
   which to reproduce the condition.
6. **The disulfide data are declared unusable by the authors and must not be resurrected.** "Only
   partial spectra of BMFD were obtained by GC/MS at the appropriate retention time and, as a result,
   data for BMFD are inconclusive and not included in this study." What survives is a **direction**:
   BMFD was present only in traces at 120 °C in a sealed vial with ~92 mL of air headspace. Set beside
   `THIOL_CHANNELS`'s `oxidative_dimerisation` entry, which records the dimer carrying up to 49 % of
   the MFT pool at 115 °C (Zhang 2024's redox series), that is a **qualitative tension worth noting
   and no more**: Zhang's system carries a deliberate oxidant (cystine) and this one does not, which
   is exactly the distinction the channel's oxidant gate is built to express. This chapter is
   therefore weak **support** for the gate, not evidence against it — and it supplies no number in
   either direction.
7. **A 4-5× method effect on MFT recovery is cited in this very paper.** Grosch & Zeiler-Hilgart
   (ref 13) found Likens-Nickerson distillation-extraction recovers **4-5 times** more MFT than direct
   extraction, and this study uses **direct extraction**. **The absolute level, 11.7 ppb, may be a
   fourfold to fivefold under-estimate on that basis alone**, and the benchmark's factor-3 contract
   does not cover a fourfold method bias. The *ratios* — the 8 % isotope share, the xylose null, the
   no-thiamine null — are immune to it, since all are within-method.
8. **The chapter is dated two ways and the repository uses both.** File stem `bolton1993`, hold-out
   label "Bolton 1994", DOI `bk-1994-0543.ch022`, printed publication date 30 November 1993, received
   9 February 1993, copyright 1994. **It is one work.** Pick one convention — the DOI's 1994 is the
   least ambiguous — and note the other in the citation so that nobody goes looking for a second
   paper. This dossier is filed under the file stem, as the brief requires.
9. **No thiamine or cysteine was measured.** Neither precursor's disappearance is quantified, so
   nothing here constrains `k_thi_hmp`, `k_thi_mesh`, or any cysteine sink — only the MFT that
   emerges at the end of one hour. The engine's thiamine lane gains a target, not a rate.
10. **One temperature, one time, one pH.** No kinetics of any kind can be extracted, and the single
    120 °C point must not be paired with another laboratory's point to make a barrier.
11. **What this chapter does not contain**: any rate or barrier; any time course; any second
    temperature or pH; any measurement of thiamine, cysteine, glucose or xylose consumption; any other
    volatile in numbers (only MFT is quantified, and BMFD is excluded); any usable disulfide datum; any
    mass balance; any figure with data in it (both figures are borrowed reaction schemes); any
    supplementary material; and any access to ref 21, the unpublished 1989-1990 data on which the
    "no MFT without thiamine" claim originally rested.
12. **What to request from the authors** (recognising that this is a 1993 chapter and the request is
    largely archival): (i) the unpublished 1989-1990 data of ref 21, which is where the
    temperature- and solids-dependence of MFT in this model system lives — **a second temperature in
    this exact matrix would be the single most valuable addition**; (ii) the BMFD spectra and areas,
    even as trace data, since the corpus has no disulfide observation in a low-a_w matrix; (iii)
    whether MFT was ever measured in a system with thiamine but **no** cysteine, which would close the
    branch measurement from the other side; (iv) the relative response factor itself and the
    calibration range; (v) whether the 34S-cysteine was checked for label scrambling into the pot's
    other sulfur pools before use.
13. **Registry gaps against `data/keys/compounds.yml`**: `2_methyl_3_furanthiol` and
    `bis_2_methyl_3_furyl_disulfide` are present, and `thiamine_availability` exists as a
    non-reactant id. **Absent: thiamine as a reactant, cysteine, D-glucose, D-xylose, glutamate/MSG,
    IMP and ribose** — i.e. every precursor in this benchmark's own pot. Also absent, and relevant to
    the mechanism the chapter argues for: **5-hydroxy-3-mercaptopentan-2-one** (the engine's `HMP`;
    see `cerny2008_extraction.md`), **3,5-dihydroxypentan-2-one**, **2-methyl-3-oxotetrahydrofuran**
    (the "oxofuran" that takes up the labelled H2S — the engine has no species for it) and
    **2-methyl-4,5-dihydro-3-furanthiol**.
