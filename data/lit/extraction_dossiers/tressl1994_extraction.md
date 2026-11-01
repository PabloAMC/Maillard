# Tressl, Kersten, Nittka & Rewicki 1994 — EXTRACTION
DOI: 10.1021/bk-1994-0564.ch018 (printed in the page margin; Crossref 2026-11-01: same title, authors Tressl / Kersten / Nittka / Rewicki, ACS Symposium Series 564, pp. 224-235, 1994)
### [1-13C]- and [6-13C]-D-glucose (and [1-13C]-D-arabinose, [1-13C]-glucose + methionine) with L-cysteine, equimolar, aqueous, initial pH 5, 160 °C, 1.5 h; GC/MS isotopomer patterns and 1H-NMR of sulfur heterocycles. **NO 2-methyl-3-furanthiol, NO norfuraneol, NO 2-furfurylthiol from glucose in any table.** ONE temperature. Position-specific single labels, not the 1:1 [13C6] design. What it does show: every C5 (and the C4) sulfur heterocycle from glucose + Cys is glucose C2-C6 (C-1 lost, via 3-deoxyglucosone), and every 1-deoxyglucosone product keeps all six carbons.

**Read 2026-11-01 from `data/articles/tressl1994.pdf`** (12 pp., pp. 224-235; ABBYY OCR of a scan, ACS
"Downloaded by KTH" stamp). The text layer garbles every table and scheme; **Tables I-V and Schemes 2-6 were read
from 130 dpi page renders, Tables IV and V again from 300 dpi crops** (scratchpad `ing1101/tr-*.png`,
`t4crop.png`, `t5crop.png`). All label shares below are printed numbers [M]; no figure is digitised (the paper has
no data figures).

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a
figure · **[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived
by me, never printed · **[!]** flag.

**Laboratory.** Tressl, Institut für Biotechnologie, Technische Universität Berlin (with Rewicki, Organische
Chemie, FU Berlin). **Not in `docs/reference/contrast_benchmark/lab_ids.yml`.** It is NOT `tu_berlin_kroh` (a
different institute and PI; no shared author). Proposed new id: `tu_berlin_tressl`. Independent of BTBU, CAAS,
Nestle / Firmenich and DFA.

## 1. What this paper answers for the model

The brief's questions (engine error E1; E6 / E11): does hexose + cysteine make MFT through an intact C5, with
which label pattern; what fraction of MFT comes from free norfuraneol vs the intact sugar chain, at which T?

| model need | what Tressl 1994 gives | tag |
|---|---|---|
| MFT label pattern from glucose + Cys | **none: MFT (2-methyl-3-furanthiol) is not reported** in any table or the text (grep of the text layer and every table render: no furanthiol, no "methyl-3-fur") | [M] / [!] |
| MFT share via free norfuraneol | **none**: norfuraneol is not reported, not fed, not measured | [!] |
| temperature dependence | **none**: one condition, 160 °C, 1.5 h, pH 5 | [M] |
| which glucose carbon a hexose-derived C5 loses | **C-1.** Furfuryl alcohol (18), 2-hydroxymethylthiophene (22), 5-hydroxymethyltetrahydrothiophene-3-one (23) and 3-thiolanone (24, a C4) are UNLABELLED from [1-13C]-glucose and 100 % singly labelled from [6-13C]-glucose (Table IV): glucose C2-C6 | [M] |
| route to that C5 | 3-deoxyglucosone, α-dicarbonyl cleavage to "2-deoxypentose" + [13C]-formic acid (C-1) (Scheme 4); "All C5-compounds" 100 % singly labelled from [6-13C], unlabelled from [1-13C] | [Q] / [P] (scheme) |
| 1-deoxyglucosone products: intact C6 or C5? | **intact C6.** All 1-DG products (28, 32, 39, 40, 41, 42, 43) are 100 % singly labelled from BOTH [1-13C]- and [6-13C]-glucose (Tables II, V); thioacetylformoin (39) "unmistakably demonstrates" the intact hexose chain | [M] / [Q] |
| pentose control | from [1-13C]-arabinose, 18, 22, 23, furfurylmercaptan (2-furfurylthiol) have intact C5 skeletons, 100 % singly labelled | [M] / [Q] |
| amounts | Tables I and II in "ppm" (basis not stated: per sugar? per extract?) | [M] / [!] |

**Answers.**

1. **Hexose + Cys -> MFT via an intact C5: NOT ANSWERED** (confidence the paper holds no MFT datum ~97 %; the
   remaining doubt is an MFT entry under a name I did not recognise, but no listed structure is a 2-methyl-3-thio
   furan). Tressl 1994 does not report MFT at all, so it gives no label pattern, no free-NF share and no
   temperature. It cannot move E1, E6 or E11 directly.
2. **What it does settle, for the hexose -> C5 entry in general [M]:** at 160 °C / pH 5, every C5 sulfur heterocycle
   identified from glucose + Cys is glucose C2-C6, i.e. C-1 is lost, through 3-deoxyglucosone (Scheme 4). The C5
   so formed has no methyl group; it gives the furfuryl skeleton (furfuryl alcohol, 2-hydroxymethylthiophene, the
   thiolanones), not the 2-methylfuran skeleton of MFT [Z].
3. **Against a 1-DG -> C5 entry, weakly [Z]:** the engine's off-by-default `r_odg_dpo` (1-deoxyglucosone ->
   1-deoxypentosone + C1, "the C-6 loss") would leave an unlabelled C5 from [6-13C]-glucose and a labelled methyl
   C5 from [1-13C]-glucose. No such product is reported: every 1-DG product here keeps all six carbons. This is
   absence of a reported C5 product, not a measured zero (no detection limit; MFT and norfuraneol were not
   targeted). Weight: low (P ~0.25 that it bears on the C-6-loss entry at all).
4. **Implication for hexose MFT [Z, a lead]:** the only hexose -> C5 cleavage this lab documents (C-1 loss from
   3-DG) cannot by itself give MFT's 2-methyl group. A hexose MFT that is an intact C5 (BTBU, CAAS, 1:1 [13C6]) must
   then come from a C-6 (or C-1 from the other end) loss on a 1-deoxy chain, or from a reduced / dehydrated 3-DG C5
   (a 1,4-dideoxypentosone-type unit) that this paper did not see. The [13C6] 1:1 design of Zhao 2019b / Du 2023
   cannot tell C-1 from C-6 loss; **a [1-13C]- vs [6-13C]-glucose + Cys pot with MFT measured is the decisive
   experiment, and Tressl's design is that experiment minus the MFT measurement.**

## 2. Conditions

| item | as printed |
|---|---|
| sugar | [1-13C]-D-glucose, [6-13C]-D-glucose, [1-13C]-D-arabinose (single-position labels; enrichment not printed) |
| amine | L-cysteine (or L-methionine); Table II also "AS / Pro" and "AS / Cys" pots, AS = Val, Ile, Leu, Phe (mixtures of primary amino acids) |
| ratio / concentration | "Equimolar amounts" [V]; absolute concentrations not printed (method "as described (7)" = Tressl et al. 1993 JAFC 41:547, not on disk) |
| medium / buffer | "aqueous solution at pH 5" [V]; buffer not named; pH initial only, drift not reported |
| temperature / time | 160 °C, 1.5 h (one condition) |
| vessel / headspace / O2 | not stated |
| work-up / analysis | pentane / ether extraction; capillary GC/MS; 1H-NMR (CDCl3) for 21, 23, 40; label share from fragment-ion intensities corrected for natural 13C (e.g. 39: (M-COCH3)/(M-CO13CH3) and m/z 43/44) |

## 3. Quantitative results

**Table I** (Strecker / β-elimination / transamination products, [1-13C]-glucose, "formation [ppm]") [M]:
pyruvate (1) 50 (Cys); 2-mercaptopropionic acid (2) 600 (Cys); 3-mercaptopropionic acid (3) 180 (Cys);
2-oxo-5-thiahexanoic acid (4) 10 (Met); 5-thiahexanoic acid (5) 1050 (Met); 6: 410 (Met); 3,5-dimethyl-1,2,4-
trithiolane (7) 20 (Cys); 8: 15 (Cys); methional (9) 1300 (Met); 10: 10 (Met); methylpyrazine-type 11: 8 (Cys) /
80 (Met); 12: 90 (Met); 2-formyl-5-methylpyrrole (13) 70 / 20; 2-acetylpyrrole (14) 110 / 100. Products 1-5 are
unlabelled in both [1-] and [6-13C] systems (intact amino-acid skeletons). In the Met system 2-methylthiopropionic
acid is a 20 : 80 mixture of unlabelled : singly labelled [M].

**Table II** ([1-13C]-glucose, ppm, AS/Pro vs AS/Cys) [M]: 15 (pyrrole lactones) 880 / 650; 16 40 / 1100;
HMF (17) 230 / 3420; furfuryl alcohol (18) 170 / 330; 28 (dihydro-dihydroxy-methylpyranone) 120 / 125; furaneol
(32) 310 / 300; cyclotene (33) 400 / 230; 2-acetylfuran (34) 35 / 90. "All C6-compounds are 100% singly labeled"
[V]. Cysteine raises the 3-DG products (HMF x15, 18 x1.9) but not the 1-DG products (28, 32 flat) [Z from M].

**Table IV** (100 % singly labelled vs unlabelled, Cys) [M]:

| compound | [1-13C]-glucose | [6-13C]-glucose |
|---|---|---|
| 5-hydroxymethylfurfurylthiol (19) | label on CH2SH (100 % M+1) | label on CH2OH |
| 5-methylfurfurylthiol (20) | label on CH2SH | label on CH3 |
| 5-acetyltetrahydrothiophene-3-one (21) | label at C-2; M 0 % : M+1 100 % : M+2 0 % | label on acetyl CH3 |
| 2-acetylpyrrole (14) | singly labelled (ring C-5 in the drawing) | singly labelled (position not resolved in the render) |
| furfuryl alcohol (18) | **unlabelled** | label on CH2OH ([1-13C]-arabinose the same) |
| 2-hydroxymethylthiophene (22) | **unlabelled** | label on CH2OH |
| 5-hydroxymethyltetrahydrothiophene-3-one (23) | **unlabelled** | label at C-2 (ring CH2-S) |
| 3-thiolanone (24, C4) | **unlabelled** | label at C-2 |

**Table V** (1-DG products, glucose + Cys; isotopomer shares as printed) [M]:

| compound | [1-13C]-glucose | [6-13C]-glucose |
|---|---|---|
| 28 | label on the 6-methyl | label on ring C-2 (CH2-O) |
| thioacetylformoin (39) | [5-13CH3] 55 % : [2-13CH3] 45 % | 45 % : 55 % |
| thiofuraneol (40) | [5-13CH3] 54 % : [2-13CH3] 46 % | 46 % : 54 % |
| furaneol (32) | [5-13CH3] 48 % : [2-13CH3] 52 % | 52 % : 48 % |
| 2,5-dimethyltetrahydrofuran-3-one (41) | 54 % : 46 % | not given |
| 42 / 2-acetylthiophene (43) | [13CH3]-42; 43 acetyl CH3 | [13CH2OH]-42; 43 ring C-5 |

Text: the 50 : 50 expected for a symmetric diacetylformoin is not met (55 : 45, 54 : 46, 48 : 52); the authors
propose a cysteine-specific route through intermediate D (Scheme 6) [P].

## 4. Mechanism evidence and its kind

- **Labelling (kind 1), position-specific single 13C**, one laboratory (Tressl, TU Berlin), one condition. Shown:
  C-1 loss for the glucose C5 heterocycles [M]; intact C6 for the 1-DG products [M]; transamination (3, 4, 5) and
  β-elimination (pyruvate, 2-mercaptopropionic acid) from unlabelled amino-acid-skeleton products [M/Q].
- **Proposed, not shown [P]:** Schemes 3-6 (3,4-dideoxyosone + Strecker H2S for 19 / 20 / 21; α-dicarbonyl cleavage
  of 3-DG to 2-deoxypentose + formic acid; thiofuraneol via D). Formic acid itself was not measured.
- No mass balance, no yields per sugar, no replicates stated, no detection limits.

## 5. Bearing on open questions

| item | verdict | why |
|---|---|---|
| E1 (hexose intact-C5 MFT; ENGINE_VS_LITERATURE secs. 3, 5) | **LEAVES UNCHANGED** | no MFT; weakly against a 1-DG -> C5 entry (no 1-DG-derived C5 product among those reported), and shows the hexose C5 that does form loses C-1 via 3-DG [Z] |
| E6 (MFT via free NF) / E11 (per-route barriers) | **LEAVES UNCHANGED** | no MFT, no NF, one T |
| F25-c (hexose FFT ~94 % an intact glucose C5, Zhao 2019b, BTBU) | **SUPPORTS the skeleton, from a second lab** | furfuryl alcohol and 2-hydroxymethylthiophene from glucose are an intact C2-C6 unit; FFT itself is listed only for arabinose ("All C5-compounds" may include it, ambiguous) [Q] |
| F23 / F24 (MFT intact C5 vs C2 + C3) | LEAVES UNCHANGED | no MFT |
| evidence map 2.12 / triage (hexose sulfur supply) | adds | cysteine raises 3-DG products (HMF x15) and leaves 1-DG products flat (Table II) at 160 °C |

## 6. Disagreements found

1. **With the brief / plan expectation (tasks/overnight_log_2026-11-01.md: "Tressl 1994 the top lead" for E1, E6 /
   E11)**: the paper has no MFT, no NF, one temperature. Cause: the lead was taken from the title and abstract
   ("[13C]-labeled sugars, cysteine"), which name furans and thiophenones but not MFT. Not a PDF problem: the file
   is the right paper (title, authors, pages, DOI match).
2. **With the engine's `r_odg_dpo` design note** ("1-deoxyglucosone -> 1-deoxypentosone + C1 (the C-6 loss)",
   parameters_sulfur.HexoseC5Entry): Tressl finds only C-1 loss (3-DG) for glucose C5 products and intact C6 for all
   1-DG products. Quantity: the carbon a hexose C5 loses; values C-6 (engine, assumed) vs C-1 (Tressl, for
   non-MFT C5s). Candidate causes: different products (MFT not measured here), 160 °C only, chemistry (a minor
   1-DG C-6 loss below the reporting level). Not a contradiction of the engine's MFT, a gap.

## 7. Ledger-ready measurements

1. **Table IV** isotopomer assignments (18, 22, 23, 24 unlabelled from [1-13C], labelled from [6-13C]): the C-1-loss
   signature of hexose C5 heterocycles at 160 °C / pH 5 (a qualitative label row; two readers trivial).
2. **Table V** isotopomer shares (39 55 / 45, 40 54 / 46, 32 48 / 52, 41 54 / 46; reversed for [6-13C]): intact C6
   for the 1-DG products, the asymmetry the authors attribute to a Cys-specific route.
3. **Table II** AS/Pro vs AS/Cys ppm (HMF 230 vs 3420): a within-study cysteine effect on 3-DG products, units
   unstated, so a ratio row only.
