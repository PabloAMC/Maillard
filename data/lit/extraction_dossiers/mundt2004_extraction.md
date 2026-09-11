# Mundt 2004 — EXTRACTION (glucose 0.25 M or maltose 0.25 M + glycine 0.25 M, 0.2 M acetate buffer pH 5.5, aqueous, 70.0 ± 0.1 C, heated to A470 = 2; melanoidins isolated by dialysis at MW > 3500 and MW > 12500 and characterised by CHN microanalysis AND by 14C labelling of sugar, whole glycine and C-1 of glycine)

### THE PAPER THE C/N DIAGNOSTIC HAS BEEN WAITING FOR: it prints the elemental carbon-to-nitrogen ratio of a glucose-glycine melanoidin with NO protein backbone — **C/N = 7.64 ± 0.21 by microanalysis and 7.61 by an independent radiochemical route** — and its radiochemical half shows WHY that number is below the trunk's structural floor of 8: about two-thirds of the glycine enters the polymer DECARBOXYLATED, carrying one carbon per nitrogen instead of two.

**Source on disk:** `data/articles/mundt2004.pdf` (5 pp., J. Agric. Food Chem. 2004, 52 (13), 4256-4260).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/mundt2004.txt`). **Table 1 — the only table in the paper — came through
clean** and is re-typed in full below; its two-tier header ("Melanoidins MW >12500" /
"Melanoidins MW >3500" as sub-headings inside the body) is preserved. The four displayed
equations in the Radiochemical Investigation section came through with their fractions broken
across lines but every numeral is legible, and they are re-typed below with the fractions
restored. Figures 1 (dialysis time course), 2 (A470 vs incorporated 14C, glucose-glycine),
3 (the same, maltose-glycine) and 4 (glucose release / maltose loss / melanoidin formation vs
time in the maltose reaction) are images: **the only slopes that reach the reader are the three
the authors themselves type into the text for Figure 2** (0.575, 0.927, 1.993). Figure 3's slopes
are NOT printed — only the two ratios derived from them. There is no supplementary material.
Repo status before this dossier: `mundt2004.pdf` has **no extraction dossier** and no citation
anywhere in `src/`; it is a new source for the melanoidin-composition question.

## 0. Identity

| field | value |
|---|---|
| Title | "Comparative Study of the Composition of Melanoidins from Glucose and Maltose" |
| Authors | Sandra Mundt and Bronislaw L. Wedzicha (corresponding), Procter Department of Food Science, University of Leeds, Leeds LS2 9JT, United Kingdom |
| Venue | J. Agric. Food Chem. 2004, 52 (13), 4256-4260. Received 24 November 2003, revised 6 April 2004, accepted 20 April 2004, web 3 June 2004 |
| DOI / article ID | 10.1021/jf035381p (printed as `JF035381P`) |
| Naming | "melanoidins" = the **non-dialysable retentate** at a stated molecular-weight cut-off, nothing else; "residue" or "subunit" = one sugar or amino-acid molecule's worth of atoms found in that retentate; "decarboxylated glycine" = the Strecker amino ketone plus the Strecker aldehyde route, i.e. glycine minus CO2 (C1 N1) |
| Lineage | the Leeds melanoidin school: Wedzicha & Kaputo 1992, Wedzicha & Vasiliauskaite 1997 (1:0.92 glucose:glycine), Wedzicha & Leong 2000 (melanoidin concentration expressed as incorporated sugar residues via U-14C). Compared throughout against Cämmerer & Kroh 1995 (C/N 7.22 under similar conditions) and Cämmerer, Jalyschko & Kroh 2002 (intact carbohydrate side chains), and against Martins & van Boekel 2003 for the temperature/pH direction |
| Companions on disk | `cammerer1994_extraction.md` (the same laboratory family's elemental-composition work), `martins2005_extraction.md` and `martins2003*_extraction.md` (the glucose-glycine network the trunk is fitted on; Martins & van Boekel 2003 is this paper's ref. 14), `brands2001_extraction.md` / `brands2002b_extraction.md` (the glucose-**casein** melanoidin the trunk's diagnostic is currently compared against), `knol2005_extraction.md` (whose eps = 282 L/mol/cm comes from Leong 1999, the same Leeds group) |

## 1. Why it matters

The trunk carries the brown polymer as two elemental pools, `MEL_C` and `MEL_N`, in mmol of
element per litre (`src/kinetic_core/species.py`, the two `Species("MEL_C", ...)` /
`Species("MEL_N", ...)` entries and the melanoidin note above them). Their quotient is computed
by `melanoidin_c_over_n` in that file and surfaced by `CoreRun.melanoidin_c_over_n` in
`src/kinetic_core/integrate.py`. `results/validation/kinetic_core_b1_fit_report.json` reports it
under `melanoidin_c_over_n_directional_diagnostic`, where it is scored **only for sign**:

```
measured_brands_glucose_casein_120C: {10: 4.01, 30: 4.10, 60: 4.22,
                                      unheated_casein_reference: 3.97}
predicted_model_c_over_n:            {10: 8.42, 30: 9.13, 60: 9.94}
commensurable_in_level: false
why_not: "Brands' melanoidin is protein-bound (casein backbone, unheated reference C/N 3.97),
          so its level is set by the protein, not by the Maillard carbon.
          Only the SIGN of the change transfers."
```

**This paper removes the reason for `commensurable_in_level: false`.** Mundt & Wedzicha's
melanoidin has no protein in it at all: the pot is glucose plus glycine plus acetate buffer, so
every carbon and every nitrogen in the isolated polymer came from the Maillard reaction itself.
There is no backbone setting the level. The measured C/N is therefore an object of the **same
kind** as `MEL_C`/`MEL_N` — a Maillard-carbon-to-Maillard-nitrogen ratio — in a way Brands'
casein number never was. It is the first level check the diagnostic has ever had.

It also supplies the second half of what a level check needs, which is a mechanism for the
mismatch. The trunk's melanoidin repeat unit is declared in `species.py` as
`MELANOIDIN_REPEAT_UNIT_CARBON = 8`, `MELANOIDIN_REPEAT_UNIT_NITROGEN = 1`, sourced to Martins &
van Boekel 2005 Table 2 step 9 ("3-DG + Gly -> melanoidins"): six carbons from 3-deoxyglucosone
plus **two** from an intact glycine. So the model's C/N has a hard structural **floor of exactly
8.0**, reached when no carbon-only addition has occurred, and every carbon-only addition (a
trapped methylglyoxal, say) drives it up. The measured value, 7.64 ± 0.21, is **below that
floor**. This paper's radiochemical half says why in one sentence: about two-thirds of the
glycine that enters the polymer has already lost its carboxyl carbon as CO2, so it contributes
**one** carbon per nitrogen, not two. Section 4 works the consequence through.

Two further things this paper bears on:

- **`FRAG_C`, the unassigned-fragment-carbon pool** (`species.py`). Mundt's summary sentence is
  that "melanoidins from glucose-glycine are formed with the conservation of all carbon and
  nitrogen atoms from glucose and glycine **except for carbon dioxide**". In the trunk's
  bookkeeping the only sinks are the polymer and `FRAG_C`; there is no CO2 sink. Whatever the
  trunk books as fragment carbon, this paper says at least one carbon per decarboxylated glycine
  residue leaves the liquid entirely.
- **The melanoidin extinction coefficient.** The `kinetic_core_b1_holdout_report.json` browning
  hold-out is run against Martins' eps = 0.64 L/(mmol*cm) at 470 nm for glycine melanoidins, and
  the report's second hold-out is Knol's eps = 282 L/(mol*cm) for asparagine. This paper's
  Figure 2 is exactly the calibration those coefficients come from — A470 plotted against
  incorporated 14C-glucose — and it prints the slope: **A470/[U-14C glucose] = 0.575 per mM**
  (section 3). That is a third, independent A470-per-sugar-residue number, at 70 C and pH 5.5.
  See Flags 5 for why it is not a drop-in replacement.

What this paper does NOT give the repository: any rate constant, any activation energy, any
temperature series, any concentration-time table, any pH series, and any melanoidin below
3500 Da. It is a composition paper, not a kinetics paper.

## 2. Methods as they matter to a model

- **Pot.** Glucose (0.25 M) **or** maltose (0.25 M) with glycine (0.25 M), dissolved in water;
  before make-up an aliquot equal to 10 % of the final volume of a **2.0 M sodium acetate /
  glacial acetic acid pH 5.5 buffer** was added, giving **0.2 M acetate** in the final mixture,
  and the pH was then adjusted to 5.5 with NaOH. Dilute aqueous, a_w ~ 1.0. Chemicals AnalaR
  grade (Sigma / Aldrich).
- **Heating.** Water bath at **70.0 ± 0.1 C**. This is a *low* temperature by the standards of
  the corpus and the reaction is run for **days**: the maltose time course in Figure 4 runs to
  120 h. There is no temperature series — one temperature, one pH, two sugars.
- **Endpoint for the microanalysis samples.** Heated "until an absorbance at 470 nm of 2 units
  was reached." So the two melanoidins compared in Table 1 are matched on **colour**, not on
  time and not on conversion. Nothing in the paper says how long that took for either sugar.
- **Isolation — this is the step that defines the object measured.** The brown solution was
  either dialysed against water (**10 x 10 L**) in Visking tubing (**MW > 12500**) or
  continuously against water for **10 days** in dialysis cassettes (**MW > 3500**). Retentate
  vacuum-dried on a rotary evaporator at 35 C, repeated with 2 x 3 mL ethanol, then kept in a
  desiccator at room temperature for up to 2 days before microanalysis (Chemistry Department,
  University of Leeds). **The number reported is the composition of the retained high-polymer
  fraction only.** Everything that passes the membrane — every coloured or colourless Maillard
  product below the cut-off, and all unreacted glucose and glycine — is discarded before the
  analysis.
- **Adequacy of dialysis.** Established by counting the 14C activity of a retentate over 8 days
  (Figure 1); the activity "appears to have reached its minimum value after the sixth day", and
  all mixtures were then dialysed for **10 days** to be safe.
- **Microanalysis.** CHN only. **C %, H % and N % are printed; oxygen is not measured and not
  printed** (see section 3 for my by-difference arithmetic and Flags 3 for why it is soft).
  The MW > 12500 entries are means of **n = 4** with standard deviations; the MW > 3500 entries
  are **single measurements with no error**, stated as such in the table footnote.
- **Radiochemistry.** 25 mL mixtures at the same concentrations, spiked individually with
  ~0.925 MBq of D-[U-14C]glucose, [U-14C]glycine or [1-14C]glycine (glucose runs), or 0.75 MBq
  D-[U-14C]maltose or 0.925 MBq of either glycine label (maltose runs). 2 mL aliquots withdrawn
  at intervals, injected into MW > 3500 cassettes, dialysed 10 days, made up to 10 mL, 1 mL
  counted in 10 mL scintillant for 100 min on a Packard TR1500C. Counting efficiency checked
  with an internal standard added to every vial after counting; the difference was "negligibly
  small", so **no quench correction was applied**.
- **The three labels and what each sees.** [U-14C]glucose counts **sugar-derived carbon**;
  [U-14C]glycine counts **all glycine-derived carbon** (2 C per intact residue, but only 1 C per
  decarboxylated residue, and — the key point — the specific activity of a decarboxylated
  uniformly-labelled glycine is **half** that of the whole molecule, because the label is spread
  evenly over both carbons); [1-14C]glycine labels **the carboxyl carbon only**, so it counts
  **intact, non-decarboxylated glycine residues alone**. The difference between the two glycine
  labels is the decarboxylated fraction, and the factor of 1/2 in the authors' equation is that
  specific-activity halving.
- **Linearity.** "It is striking that all data sets are excellent straight-line graphs,
  indicating that the relationship between the absorbance and concentration is constant
  throughout the browning reaction." So the composition is **constant over the whole course of
  browning** at this temperature and pH — not just at the endpoint. That is a stronger statement
  than a single endpoint analysis and it is what licenses comparing a ratio at all.
- **Sugar assay (maltose runs only).** R-Biopharm enzyme kit (alpha-glucosidase, hexokinase,
  glucose-6-phosphate dehydrogenase) for maltose and D-glucose, chosen for specificity without
  separation. 1 mL aliquots, cooled to room temperature, diluted, assayed per the supplier.
- **What is not controlled.** No inert atmosphere is described; no stirring is described; the pH
  is set at t = 0 by a 0.2 M acetate buffer against 0.25 M glycine and 0.25 M sugar and is
  **never re-measured**. Read 5.5 as an initial pH.

## 3. Tables re-typed

### Table 1. "Elemental Composition of Nondialyzable Melanoidins Prepared from Sugar (0.25 M)-Glycine (0.25 M) Reaction Mixtures in Acetate Buffer (0.2 M), pH 5.5, Heated at 70.0 ± 0.1 °C"

Footnote `a` exactly as printed: "Errors are standard deviations (n = 4). Data with no errors
shown were single measurements."

| reaction mixture | C% | H% | N% | C/N |
|---|---|---|---|---|
| **Melanoidins MW >12500** | | | | |
| glucose−glycine | 42.9 ± 1.4 | 5.4 ± 0.2 | 6.6 ± 0.4 | **7.64 ± 0.21** |
| maltose−glycine | 43.6 ± 0.8 | 5.7 ± 0.1 | 4.84 ± 0.1 | 10.55 ± 0.3 |
| **Melanoidins MW >3500** | | | | |
| glucose−glycine | 42.25 | 5.2 | 6.25 | **7.88** |
| maltose−glycine | 43.55 | 5.65 | 4.75 | 10.7 |

That is the whole of Table 1 and it is the only table in the paper.

### The four displayed equations of the Radiochemical Investigation, re-typed

**Glucose-glycine, slopes of Figure 2** (A470 against incorporated label, all concentrations
millimolar; these three numbers are *printed in the text*, not read off the figure):

| slope | value | unit as the authors define it |
|---|---|---|
| A470 / [U-14C glc] | **0.575** | per mM of incorporated glucose residues |
| A470 / [U-14C gly] | **0.927** | per mM of incorporated glycine-derived carbon (see the 1/2 caveat) |
| A470 / [1-14C gly] | **1.993** | per mM of incorporated **intact** glycine residues |

**Glucose-glycine, the composition derived from them** (printed):

- `[1-14C gly] / [U-14C glc]` = **0.289 mol/mol of glucose** — intact, non-decarboxylated glycine.
- `[U-14C gly] / [U-14C glc]` = **0.62** — the total glycine-carbon count in specific-activity units.
- `[decarbox. gly] / [U-14C glc]` = 2 x (0.62 − 0.289) = **0.662 mol/mol of glucose**.
- Ratio of glucose- to glycine-derived subunits = **1 : 0.95** (i.e. 1 : 0.289 + 0.662), "in
  remarkably good agreement with the previously reported value of 1:0.92" (Wedzicha &
  Vasiliauskaite 1997).

**Maltose-glycine, the same quantities** (printed; the underlying Figure 3 slopes are NOT printed):

- `[1-14C gly] / [U-14C mal]` = **0.538 mol/mol of maltose** — intact glycine.
- (the total-glycine ratio appears only inside the next equation, as **0.991**)
- `[decarbox. gly] / [U-14C mal]` = 2 x (0.991 − 0.538) = **0.906 mol/mol of maltose**.
- Total maltose- to glycine-derived subunits = **1 : 1.44**.

**The radiochemical C/N reconstruction** (printed in the Comparison section, worked in full by
the authors for glucose):

> whole glycine : decarboxylated glycine : glucose = **0.289 : 0.662 : 1**; multiplying by the
> carbons each contributes (6 for glucose, 2 for whole glycine, 1 for decarboxylated glycine) and
> summing gives **7.24 C atoms** to **0.289 + 0.662 = 0.951 N atoms**, hence **C/N = 7.61** for
> the glucose-glycine melanoidins — against **7.64 ± 0.21** by microanalysis, which the authors
> call "identical".
>
> The same calculation for maltose (12 C per maltose residue) gives **C/N = 9.69**, against
> **10.53 ± 0.22** by microanalysis — "slightly lower". (On the 10.53 ± 0.22 / 10.55 ± 0.3
> discrepancy between the text and Table 1, see Flags 2.)

The authors' own conclusion from that agreement, quoted because it is the structural claim:
to reproduce the microanalysis it is "necessary to include in the melanoidins carbon and
nitrogen atoms from the **amino carbonyl** formed in the Strecker degradation reaction as well
as the Strecker aldehyde itself", and "the amino acid is incorporated into melanoidins
completely except for carbon dioxide, when only the decarboxylated portion of the amino acid is
incorporated."

### Numbers printed in the running text (everything else is figure-only)

| quantity | value | where |
|---|---|---|
| C/N, glucose-glycine, abstract rounding | 7.6 ± 0.2 | Abstract |
| C/N, maltose-glycine, abstract rounding | 10.5 ± 0.2 | Abstract |
| C/N, glucose-glycine, Cämmerer & Kroh 1995 under "similar conditions" | 7.22 | Results, citing ref. 7 (**not measured here**) |
| C/N invariance with molecular weight | the MW > 12500 and MW > 3500 values "were found to be the same", so C/N "does not change with molecular weight in this range" | Results |
| direction of C/N with temperature and pH | amino-acid incorporation **increases** with **decreasing** temperature, giving a **lower** C/N; same direction for decreasing pH | Results, citing refs. 7 and 14 (Martins & van Boekel 2003) |
| decarboxylated fraction of incorporated glycine | "approximately two-thirds", in **both** the glucose and the maltose reaction | Results and Summary |
| dialysis equilibration | 14C activity at its minimum after the **6th** day; 10 days used | Methods, Figure 1 |
| glucose released in the maltose-glycine reaction | rises to **~13 mM** over the 120 h observation period, with an induction phase | Results, Figure 4 |
| maltose lost over the same period | "falls by a similar amount" (no number) | Results, Figure 4 |
| melanoidin formed in the maltose reaction | **~2 mM** based on maltose subunits, over the observation period | Results |
| whole-maltose incorporation needed to explain the high maltose C/N | **1-2 mM**, "below our limit of detection of concentration changes of maltose (present at an initial concentration of 250 mM)" | Results |
| Figure 4 trend-line exponents | power functions with exponents **1.44** (glucose release) and **2.21** (melanoidin formation) | Figure 4 caption |
| Figure 4 melanoidin scaling | melanoidin values "multiplied by five" for shape comparison only | Figure 4 caption |

**Everything in Figures 1, 2, 3 and 4 is FIGURE-ONLY** and is not typed as a number here, except
the three Figure 2 slopes and the two Figure 4 exponents, which the authors print in text and
caption respectively.

### Arithmetic on the printed numbers (all mine)

**1. The C/N values recomputed from the printed C % and N % (mine).** Using C = 12.011,
N = 14.007 g/mol:

| fraction | (C%/12.011) / (N%/14.007) | printed C/N | agreement |
|---|---|---|---|
| glucose, MW > 12500 | 3.572 / 0.4712 = **7.58** | 7.64 ± 0.21 | 0.8 % low |
| maltose, MW > 12500 | 3.630 / 0.3455 = **10.51** | 10.55 ± 0.3 | 0.4 % low |
| glucose, MW > 3500 | 3.518 / 0.4462 = **7.88** | 7.88 | exact |
| maltose, MW > 3500 | 3.626 / 0.3391 = **10.69** | 10.7 | exact |

The two single-measurement rows reproduce exactly; the two n = 4 rows are ~0.5 % off, which is
what you expect when a mean C/N is computed **per replicate and then averaged** rather than from
the averaged percentages. The table is internally consistent. **The printed C/N is the primary
number and should be used, not my recomputation.**

**2. Oxygen by difference (mine, soft — see Flags 3).** Assuming the polymer is C, H, N, O only
and the ash is zero: glucose MW > 12500 gives **O ≈ 45.1 %**; glucose MW > 3500, **46.3 %**;
maltose MW > 12500, **45.9 %**; maltose MW > 3500, **46.05 %**. These are not printed and the
paper never claims the residue is only CHNO.

**3. H/C atomic ratio (mine).** glucose MW > 12500: (5.4/1.008)/(42.9/12.011) = **1.50**;
maltose MW > 12500: **1.57**. Both are far above the "conjugated carbon double bonds" skeleton
the introduction invokes, i.e. most of the polymer mass is not chromophore.

**4. Carbon per nitrogen from each side, glucose-glycine (mine, from the printed 0.289 : 0.662 : 1).**
Sugar side: 6 C for 0.951 N = **6.31 C/N**. Amine side: (2 x 0.289 + 1 x 0.662) = 1.240 C for
0.951 N = **1.30 C/N**. Total **7.61**. Compare the trunk's step-9 repeat unit: sugar side
**6.00**, amine side **2.00**, total **8.00**. The whole gap between the trunk's floor and the
measurement is in the amine term, and it is 0.70 C per N.

**5. What the trunk's floor would be if step 9 decarboxylated (mine).** With the same
two-thirds/one-third split and exactly one amine per sugar (which the trunk assumes), the amine
contributes 2/3 x 1 + 1/3 x 2 = **1.333 C per N** and the floor becomes **7.33**. Mundt's 7.61
is above that because his polymer takes only 0.951 amine residues per sugar residue, not 1.000.
So the measured value sits **between** a fully-intact-glycine model (8.00) and a
two-thirds-decarboxylated model with perfect 1:1 stoichiometry (7.33), which is exactly where the
measured 1 : 0.95 stoichiometry puts it. **This is a consistency check on the paper, not a
proposed model change.**

**6. Glycine consumed per sugar and the nitrogen budget (mine).** At 1 : 0.95, forming the
~2 mM of melanoidin the maltose run reaches consumes ~2.9 mM of glycine out of 250 mM — about
**1.2 %** of the initial amine. The polymer is a small sink in this pot; the C/N ratio is a
composition of what *did* polymerise and says nothing about how much did.

**7. Maltose: the two methods disagree and the direction is informative (mine).** Radiochemical
9.69 vs microanalytical 10.53. The radiochemical route counts only atoms that arrived as
*labelled maltose or labelled glycine*; the microanalysis counts every carbon in the retentate.
The microanalytical value being **higher** means there is carbon in the maltose polymer that the
labelled-maltose count does not see per unit nitrogen — consistent with the authors' own
proposal of intact maltose side chains at the 1-2 mM level. For **glucose the two agree to
0.4 %** (7.61 vs 7.64), so no such unaccounted carbon is needed there. **The glucose number is
the clean one, and it is the one the trunk needs.**

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Neither glucose, nor maltose, nor
glycine, nor melanoidins is keyed; the registry is a product/marker list and carries no Maillard
reactants and no polymer pool. Nothing in this paper has a registry id.

Every row below shares these conditions unless stated otherwise: **0.25 mol/L sugar + 0.25 mol/L
glycine, 0.2 mol/L sodium acetate/acetic acid, initial pH 5.5, dilute aqueous (a_w ~ 1.0), water
bath at 70.0 ± 0.1 C, heated to A470 = 2, polymer isolated by exhaustive dialysis at the stated
molecular-weight cut-off, CHN microanalysis of the dried retentate.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **C/N of glucose-glycine melanoidin, MW > 12500** | **7.64 ± 0.21** | mol C per mol N (atomic) | as above, n = 4 | Table 1 p. 4257 | **elemental_analysis** |
| C/N of glucose-glycine melanoidin, MW > 3500 | 7.88 | mol C per mol N | as above, **single measurement** | Table 1 | elemental_analysis |
| C/N of maltose-glycine melanoidin, MW > 12500 | 10.55 ± 0.3 (text prints 10.53 ± 0.22; Flags 2) | mol C per mol N | as above, n = 4 | Table 1 | elemental_analysis |
| C/N of maltose-glycine melanoidin, MW > 3500 | 10.7 | mol C per mol N | as above, single measurement | Table 1 | elemental_analysis |
| C, H, N of glucose-glycine melanoidin, MW > 12500 | C 42.9 ± 1.4; H 5.4 ± 0.2; N 6.6 ± 0.4 | mass % | as above, n = 4 | Table 1 | elemental_analysis |
| C, H, N of glucose-glycine melanoidin, MW > 3500 | C 42.25; H 5.2; N 6.25 | mass % | as above, single | Table 1 | elemental_analysis |
| C, H, N of maltose-glycine melanoidin, MW > 12500 | C 43.6 ± 0.8; H 5.7 ± 0.1; N 4.84 ± 0.1 | mass % | as above, n = 4 | Table 1 | elemental_analysis |
| C, H, N of maltose-glycine melanoidin, MW > 3500 | C 43.55; H 5.65; N 4.75 | mass % | as above, single | Table 1 | elemental_analysis |
| **C/N of glucose-glycine melanoidin, radiochemical route** | **7.61** | mol C per mol N | MW > 3500 cassettes, 14C, same pot | Results p. 4259 | elemental_analysis (independent method, authors' own arithmetic) |
| C/N of maltose-glycine melanoidin, radiochemical route | 9.69 | mol C per mol N | as above | Results p. 4259 | elemental_analysis (same) |
| **intact glycine incorporated per glucose residue** | **0.289** | mol/mol | as above, from Figure 2 slopes | Results p. 4258 | within_study_ratio |
| **decarboxylated glycine incorporated per glucose residue** | **0.662** | mol/mol | as above | Results p. 4258 | within_study_ratio |
| total glycine-derived subunits per glucose residue | 0.951 (printed as the ratio 1 : 0.95) | mol/mol | as above | Results p. 4258 | within_study_ratio |
| **decarboxylated share of incorporated glycine** | **0.662 / 0.951 = 0.696, i.e. ~2/3** (mine; the paper says "approximately two-thirds") | — | as above | derived from the two rows above (mine) | within_study_ratio |
| intact glycine incorporated per maltose residue | 0.538 | mol/mol | as above, maltose pot | Results p. 4258 | within_study_ratio |
| decarboxylated glycine incorporated per maltose residue | 0.906 | mol/mol | as above | Results p. 4258 | within_study_ratio |
| total glycine-derived subunits per maltose residue | 1.44 | mol/mol | as above | Results p. 4258 | within_study_ratio |
| A470 per mM incorporated glucose residues | 0.575 | absorbance units per mM | 1 cm assumed but **not stated**; 70 C, pH 5.5, MW > 3500 | Results p. 4258 (Figure 2 slope, printed in text) | within_study_ratio (see Flags 5) |
| A470 per mM incorporated glycine-derived carbon | 0.927 | absorbance units per mM | as above | Results p. 4258 | within_study_ratio |
| A470 per mM incorporated intact glycine | 1.993 | absorbance units per mM | as above | Results p. 4258 | within_study_ratio |
| C/N is independent of molecular weight above 3500 | 7.64 vs 7.88 (glucose); 10.55 vs 10.7 (maltose) | — | as above | Table 1 + Results | within_study_ratio (structural) |
| oxygen by difference | 45.1 / 46.3 / 45.9 / 46.05 | mass % | the four Table 1 rows in order | derived (mine) | derived_assumption (assumes CHNO only, zero ash) |
| H/C atomic ratio | 1.50 (glucose, MW > 12500); 1.57 (maltose) | mol/mol | as above | derived (mine) | derived_assumption |
| glucose released during maltose-glycine browning | ~13 | mM at 120 h | maltose pot only | Results (text figure quoted from Figure 4) | figure_only (level quoted in text) |
| melanoidin formed, maltose pot | ~2 | mM in maltose subunits, over 120 h | maltose pot only | Results | figure_only (level quoted in text) |
| C/N, glucose-glycine, from Cämmerer & Kroh 1995 | 7.22 | mol C per mol N | "similar conditions", **not measured here** | Results, ref. 7 | level_only (borrowed) |
| glucose release / maltose loss / melanoidin time courses; A470 vs 14C plots; dialysis time course | — | — | — | Figures 1-4 | **figure_only** |

### What the trunk's diagnostic should be compared with, and on what basis

**The comparison object is 7.64 ± 0.21 (Table 1, glucose-glycine, MW > 12500), with 7.88 as the
MW > 3500 companion and 7.61 as the independent radiochemical confirmation.** All three are the
same pot at 70 C and pH 5.5. The trunk's `melanoidin_c_over_n` at 120 C runs 8.42 → 9.13 → 9.94
over 10 → 30 → 60 min (`kinetic_core_b1_fit_report.json`). Five things have to be said before
those two numbers are put on one line.

**(a) The kind of object now matches, and that is the change this paper makes.** The fit
report's `why_not` string turns entirely on the casein backbone: Brands' unheated reference is
already C/N 3.97 before any Maillard chemistry, so his 4.01 → 4.22 is a small Maillard
perturbation on a protein. **Mundt's polymer contains no protein.** Its carbon and nitrogen are
Maillard carbon and Maillard nitrogen and nothing else. The `commensurable_in_level: false`
verdict was correct about Brands and is **not** correct about this paper. The diagnostic can be
compared in level here.

**(b) But an isolated fraction is not the lumped pool, and this is the load-bearing caveat.**
`MEL_C`/`MEL_N` count every carbon and every nitrogen the model routes into the polymer, of any
size, soluble or not, coloured or not, from the first step-9 event onwards. Mundt's number is the
composition of what **stayed inside a dialysis membrane** at 3500 or 12500 Da after ten days
against 100 L of water. Every low-molecular-weight brown product, every oligomer below the
cut-off and every unreacted reactant is thrown away before the analysis. The two objects can
have different C/N for reasons that have nothing to do with the model being wrong. **The paper
does supply the one piece of evidence that softens this**: the C/N of the > 12500 and the > 3500
fractions agree (7.64 vs 7.88, and 10.55 vs 10.7), so composition is at least flat across the
high-polymer range, and Figure 2's straight lines say it is flat in time as well. It is still an
extrapolation to say it is flat all the way down to a dimer, and the paper does not say that.

**(c) The temperature and pH are wrong for the fit report's row, and the paper says which way
that pushes.** The diagnostic's numbers are at **120 C**; Mundt's are at **70 C**. The paper
states the direction explicitly — amino-acid incorporation increases as temperature falls, so
C/N falls with falling temperature, and the same for pH — and cites Cämmerer & Kroh and Martins
& van Boekel for it. So a 120 C melanoidin should have a C/N **above** 7.64. **The trunk's
8.42-9.94 is therefore in the right direction and of the right order**, and the honest statement
is a **bracket, not a match**: 7.64 at 70 C is a lower bound on what a 120 C glucose-glycine
melanoidin should show, and the trunk clears it. The paper carries no temperature series of its
own, so nothing here converts 7.64 to a 120 C expectation.

**(d) The trunk's structural floor is above the measurement, and that is a real finding.**
`MELANOIDIN_REPEAT_UNIT_CARBON = 8` with `..._NITROGEN = 1` means the model's C/N **cannot go
below 8.0** by construction. The measured value at 70 C is 7.64 ± 0.21 — about 1.7 analytical
standard deviations below the floor — and the radiochemical route independently gives 7.61. So
the disagreement is not in the fitted rate constants; it is in the **stoichiometry of step 9**,
which books an intact C2 N1 glycine into the polymer while the measurement says roughly
two-thirds of the glycine arrives as C1 N1 having lost CO2. The size of the effect is 0.70
carbon per nitrogen (section 3, arithmetic 4), or **8.8 % of the model's floor value**. Two
consequences worth stating plainly:

  - The trunk's `melanoidin_repeat_units` helper returns `MEL_N` on the argument that "every
    step-9 event contributes exactly one nitrogen". **That argument survives this paper
    untouched** — a decarboxylated glycine still carries exactly one nitrogen. Only the carbon
    count changes.
  - If step 9 were rewritten to emit CO2 for two-thirds of events, one carbon per event would
    have to go somewhere. It cannot go to `MEL_C` and it is not a fragment in solution, so the
    trunk would need a CO2 sink it does not have, or would have to route it to `FRAG_C` and
    accept that `FRAG_C` then contains gas. **This is a design question, not an extraction
    finding, and it is recorded as Flags 8, not as a recommendation.**

**(e) The comparison is a level check on a diagnostic, and must stay one.** The fit report
already marks the C/N row `"role": "directional diagnostic, not a fit target and not a scored
hold-out"`. Nothing in this paper makes it a fit target: there is no time course, no temperature
series, and one pot. What it can become is a **second, level-commensurable directional row
alongside the Brands row**, with the honest annotation that the level agreement is a bracket
(7.64 at 70 C < 8.42 at 120 C, direction correct) and that the model's floor of 8.0 is a
falsifiable structural claim which this measurement is the first evidence against.

**(f) What cannot be transported at all.** No rate, no barrier, no time, no temperature
dependence, no water-activity point, no amine other than glycine, no real food. The three A470
slopes are the one quantity that touches a fitted object (the extinction coefficient), and they
are not directly usable — see Flags 5.

## 5. Flags

1. **The measured object is a dialysis retentate, not a model pool.** Say this every time the
   7.64 is quoted. An elemental ratio measured on an isolated MW > 12500 fraction and the C/N of
   a lumped `MEL_C`/`MEL_N` sink are different objects; they agree only if composition is flat
   in molecular weight all the way down, which the paper demonstrates over 3500-12500 and
   nowhere else.
2. **The maltose C/N is printed twice with two different values and two different errors.**
   Table 1 says **10.55 ± 0.3**; the Comparison section, referring to the same microanalysis,
   says **10.53 ± 0.22**; the Abstract rounds to **10.5 ± 0.2**. The glucose value is stable
   (7.64 ± 0.21 in the table, 7.64 ± 0.21 in the text, 7.6 ± 0.2 in the abstract). Use the table
   for both and note the maltose discrepancy. It does not affect anything the trunk needs.
3. **Oxygen is never measured.** Only C, H and N are printed. No ash, no sulfur (there is none in
   this system), no residual moisture figure after the 35 C rotary evaporation and 2 days in a
   desiccator. My by-difference oxygen in section 3 assumes CHNO and zero ash; **the C/N ratio
   itself does not depend on any of this**, which is precisely the authors' own point that C/N is
   "the most reliable parameter obtained from microanalysis data ... as it was calculated with no
   model assumption". Do not carry the by-difference oxygen into anything.
4. **One temperature, one pH, one amine, one endpoint.** 70.0 C, pH 5.5, glycine, A470 = 2.
   There is no series of any kind in this paper. The temperature and pH directions it asserts
   are borrowed from refs. 7 and 14, not measured here.
5. **The three A470 slopes are not a drop-in extinction coefficient.** They are absorbance per
   mM of *incorporated 14C label in the MW > 3500 retentate*, at an unstated (presumably 1 cm)
   path length, at 70 C and pH 5.5. Martins' eps = 0.64 L/(mmol*cm) and Knol's 282 L/(mol*cm) are
   absorbance per mM of *melanoidin* on their own definitions of a melanoidin molecule. The
   0.575 figure is per **glucose residue incorporated**, which is a third definition again. They
   can be compared only after the definitions are reconciled, and the path length must be
   confirmed from the paper's original figure axes.
6. **Figure 3's slopes are not printed.** The maltose composition numbers (0.538, 0.991, 0.906,
   1.44) are the authors' arithmetic on slopes the reader never sees. The glucose numbers are
   traceable to three printed slopes; the maltose ones are not. Weight them accordingly.
7. **The two glycine labels do not measure the same thing, and one of them needs a
   specific-activity correction the authors apply by hand.** The factor of 1/2 in
   `[decarbox.gly] = 2 x (0.62 − 0.289)` is because a uniformly-labelled glycine that loses one
   of its two carbons retains half its specific activity. This is correct as written, but it
   means the decarboxylated count is a **difference of two numbers of similar size** (0.62 and
   0.289) multiplied by two, so its relative error is roughly twice the relative error of either
   slope, and **no error bar is given for it anywhere**. The 0.662 and the "two-thirds" that
   follows from it are point estimates with no stated uncertainty.
8. **The trunk's step-9 stoichiometry is the thing this paper contradicts, and the fix is not
   free.** `MELANOIDIN_REPEAT_UNIT_CARBON = 8` asserts an intact glycine. This paper says ~2/3 of
   incorporated glycine is decarboxylated. Changing the constant would change the browning
   readout's conversion between `MEL_N` and absorbance, would require a CO2 sink the trunk has
   no species for, and would touch the conservation invariant `network.validate_balance()`
   enforces at import. **Recorded here as an observation with a mechanism attached; it is a
   design decision for a later wave, not an extraction result.** Note also that the trunk's
   source for the constant, Martins & van Boekel 2005 step 9, is the same laboratory family
   Mundt cites (ref. 14) for the temperature and pH direction.
9. **"Approximately two-thirds" is the authors' own rounding of 0.696.** They use the phrase for
   both sugars: glucose gives 0.662/0.951 = 0.696 and maltose 0.906/1.444 = 0.627 (both mine).
   The two are not equal and neither is exactly 2/3. Carry the underlying pairs, not the phrase.
10. **What this paper does not contain**: any rate constant; any activation energy; any
    concentration-time table; any temperature or pH series; any melanoidin below 3500 Da; any
    amine other than glycine; any water-activity point; any real-food matrix; any oxygen or ash
    determination; any molecular-weight distribution finer than the two dialysis cut-offs; any
    error bar on the radiochemical composition; any supplementary material.
11. **What to request from the authors**: (i) the numeric slopes behind Figure 3 with their
    standard errors, so the maltose composition becomes traceable the way the glucose one is;
    (ii) the heating times to A470 = 2 for the two sugars, which the paper never states and
    without which the composition cannot be placed on any time axis; (iii) the path length for
    the A470 measurements; (iv) whether the microanalysis was corrected for residual moisture or
    ash; (v) the per-replicate C, H and N values behind the n = 4 means, which would settle the
    0.8 % gap between the printed C/N and my recomputation from the printed percentages.
12. **Registry gaps against `data/keys/compounds.yml`**: none of `glucose`, `maltose`, `glycine`
    or any melanoidin pool is keyed among the 75 ids. The trunk's internal names `Glc`, `Gly`,
    `MEL_C`, `MEL_N` and `FRAG_C` are network-local and are not registry ids. A composition
    benchmark built from this paper would need at least glucose and glycine keyed, plus a
    convention for keying a polymer pool that has no molecular weight — which
    `species.py` already declares as the reason `MEL_C` and `MEL_N` have no entry in
    `MOLECULAR_WEIGHT_G_PER_MOL`.
