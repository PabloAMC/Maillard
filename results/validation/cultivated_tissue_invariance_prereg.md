# Pre-registration: does the kinetic layer change what to add to cultivated tissue? (2026-09-13)

## 1. The question, and why it comes before the module

Cultivated muscle and fat are grown. They are not exercised, not aged, not bled. Their free amino
acid, sugar, nucleotide and lipid pools therefore differ from slaughtered tissue, and those pools are
this engine's inputs. A cultivated tissue enters the front door as a **spec**: precursors in
millimolar, a temperature programme, pH, water activity. No species, no reaction and no parameter is
added. Rule 1 of the roadmap ("one engine") is untouched, and `core_panel_scores.json` cannot move,
because the scorecard scores the engine against benchmarks and a spec is an input.

That is also the problem. If the engine adds nothing to a question that is entirely about inputs,
then the work is a literature corpus with a model bolted on for decoration.

**The claim under test.** A flavour chemist handed a composition table for cultivated muscle beside
beef can already rank which precursor to restore: the one most missing, among those the meaty lanes
consume. The engine is worth running only if it disagrees, and only if the disagreement survives both
the uncertainty in the composition and the uncertainty in the rates.

This matters more than usual here because the engine's disagreement, if any, would come from the
sulfur lane's competition and saturation terms, and the sulfur lane is the lane the core panel scores
worst (0 of 37 strict-ready; the thiol sink unpinned, `docs/validation/thiol_sink_candidates.md`).
The place where the model would add value is the place it is least able to.

## 2. The design

**No composition is asserted.** Three specs are built: a beef reference, cultivated muscle,
cultivated fat. Every precursor concentration is a declared **interval**, not a value, and the run
sweeps the box. A range with a source carries it; a range without one is labelled `stub` and is
widened, not narrowed. Nothing written by this run enters `data/`; the box and the artifact live
under `results/cultivated_tissue_invariance/`.

Sweeping the box rather than guessing point values is the whole design. If the ordering is stable
across a deliberately wide box, the stub values did not matter. If it flips inside the box, the
answer is that composition must be measured before anything is predicted, which is itself a result.

**Declared before the run, not chosen after:**

- **Cooking programme.** 100 °C, 20 min, pH 6.0, a_w 0.98, matrix `water` with no protein loading.
  100 °C is the lower edge of the fit corpus (100–145 °C), so the run sits at the corpus boundary
  and not outside it. A second programme at 140 °C, 5 min is run as a robustness arm.
- **Target set.** 2-methyl-3-furanthiol, 2-furfurylthiol, methional, 2,5-dimethylpyrazine,
  hexanal, and the Strecker aldehydes the panel already scores. Ranked on summed odour-activity
  ratio over this set.
- **Candidate precursors.** Those the trunk, sulfur and Strecker lanes consume and that plausibly
  differ between grown and slaughtered tissue: ribose, ribose-5-phosphate, glucose, cysteine,
  methionine, thiamine, IMP/inosine, leucine, the lipid pool.

**Per draw of the composition box:**

1. **Naive ranking N.** Candidate precursors ordered by absolute molar deficit against that draw's
   beef reference, restricted to the candidate set. This is the ranking that needs no model.
2. **Engine ranking E.** Each precursor restored separately to its beef level, the programme
   integrated, precursors ranked by the resulting change in summed OAV over the target set.
3. Record whether `top(E) == top(N)`, and Kendall tau between E and N.

**Per candidate reversal**, the parameter envelope is then drawn on top of the composition draw, so
a reversal is only counted if it survives the rates as well as the composition.

## 3. What counts as success, declared before the run

- **T1 — the engine earns the module.** In at least 20 % of composition draws `top(E) != top(N)`;
  the disagreement is *concentrated*, meaning one specific reversal (precursor A over precursor B)
  accounts for at least half of those draws rather than scattering across pairs; and that reversal
  holds in at least 80 % of the parameter-envelope draws taken on top of it.
- **T2 — the engine is decoration.** Mean Kendall tau between E and N is at least 0.9 **and**
  `top(E) == top(N)` in at least 90 % of draws. The composition gap map ships; no module is built.
- **T3 — the engine refuses.** If the declared targets are named refusals under the cultivated specs
  in at least half of the draws, the sweep answers nothing. The outcome is T2, and the refusal
  reasons are the deliverable: they name exactly which measurement is missing.
- **Anything between T1 and T2 is indeterminate**, is reported as indeterminate, and is not a
  licence to build the module. An indeterminate result is resolved by measuring composition, not by
  rerunning with a narrower box.

## 4. What this run will not do

- It will not add a species, a reaction or a parameter to `src/kinetic_core`.
- It will not write into `data/`.
- It will not move `results/validation/core_panel_scores.json`. A test asserts it is unchanged with
  this branch's code present.
- It will not claim a composition for cultivated tissue. The box is a sensitivity device. No number
  in it is a measurement, and the artifact labels each range `sourced` or `stub`.
- It will not be cited as evidence about cultivated meat. It is evidence about this engine.

## 5. What follows from the outcome

**On T1**, the module is built, inside this repository, with the separation carried by provenance
rather than by a second project:

- `data/tissues/cultivated.yml` — composition rows, each with `class: declared_assumption`, a band,
  and an extraction dossier under `data/lit/extraction_dossiers/`. Never merged into
  `protein_matrices.yml`, whose entries are measured site densities.
- `results/validation/cultivated_panel_scores.md` — a second scorecard with its own benchmark count.
  It reads zero benchmarks, which is correct and is the headline.
- A flag on every output row and report header of any run whose spec names a cultivated tissue,
  carried by the existing refusal machinery, so the module cannot emit an unflagged number.
- Two tests: the flag cannot be stripped; the core scorecard is unchanged.
- Its own guide page under `docs/guides/`, cross-linked once from `INTRODUCTION.md`, not a new
  section inside it.

**On T2 or T3**, the deliverable is the composition corpus and its gap map: what has been measured in
cultured muscle and fat, what has not, and the one experiment that would close the largest gap. That
corpus is the asset either way, and it is citable without this engine.

## 6. Amendments before the run (2026-09-13, same day, no draw taken yet)

The engine was asked, before any sweep, which of section 2's declared names it can take. The
answers below change the design. Every change is recorded here, dated, before the first draw, and
none of section 3's logic is relaxed; where a threshold is re-declared, the reason is the size of
the candidate set, which was not known when section 3 was written.

**A1 — the lanes do not compose, so the engine cannot rank across the meaty set at once.** The
three Maillard lanes refuse each other (`engine.resolve_lanes`). Cysteine, ribose and thiamine
force the sulfur lane; methional and the pyrazines run on the trunk only and are carried inert in
the sulfur network; hexanal needs a lipid carrier and no tissue-fat carrier exists. So:

- The **sulfur arm is the decision arm.** Precursors ribose, cysteine, thiamine, glucose. Targets
  2-methyl-3-furanthiol, 2-furfurylthiol, furfural, hydrogen sulfide, methanethiol, furaneol.
- The **trunk arm is refused by construction.** Asked for methional, 2-methylpyrazine,
  2,5-dimethylpyrazine and furaneol on glucose + glycine, the engine refuses methional (the
  methionine chain did not ship, wave B22), returns the pyrazines at ~1e-13 µg/L from a sugar +
  amine pot, and has no water threshold for any of the four. No decision metric exists on that
  arm. Recorded as a refusal; not swept.
- **Cultivated fat is refused by construction.** No route from tissue lipid to any lane. Recorded.
- **Hexanal is dropped from the target set** for the same reason. Recorded.

**A2 — three of the eight declared candidate precursors are not species in any lane.** Leucine,
IMP/inosine and ribose-5-phosphate are refused (`UNMAPPED PRECURSORS`). The engine cannot rank
them. The naive ranking is therefore computed over the four the engine takes, so that E and N are
compared on the same set; the naive ranking over all eight is also reported, for the reader, with
the three unrankable ones marked. This is the first structural refusal of the run and is reported
as such in the artifact; it does not trigger T3, which is about targets.

**A3 — the decision metric is summed odour-activity over targets with a measured water threshold.**
Of the sulfur arm's six targets, three carry a threshold in the corpus (MFT 0.005 µg/L, FFT
0.006 µg/L, furfural 3000 µg/L; Zhou 2023 SI Table S2). Hydrogen sulfide, methanethiol and
furaneol have `no_measured_threshold_for_this_matrix` and are reported as concentration ratios,
outside the metric. In practice the metric is MFT + FFT.

**A4 — the naive ranking is fractional deficit, not absolute.** Section 2 said absolute molar
deficit. Glucose sits at millimolar and thiamine at micromolar, so absolute deficit would rank
glucose first in nearly every draw, which no chemist would do; it would be a straw man that hands
T1 to the engine. N is now `1 − cultivated/beef` per precursor, ties broken by absolute deficit. A
precursor whose cultivated draw is at or above its beef draw has nothing to restore and is dropped
from that draw's ranking in both N and E.

**A5 — Kendall tau re-declared for four candidates.** With n = 4 one adjacent swap gives τ = 0.67
and τ ≥ 0.9 means identical order, which makes T2 unreachable by construction. T2 is now
`top(E) == top(N)` in at least 90 % of draws **and** mean τ ≥ 0.75. T1 is unchanged.

**A6 — declared, fixed, not swept:** phosphate buffer 0.03 mol/L in every spec (the muscle
phosphate pool, order of magnitude; identical across arms, so it cancels in the ranking);
`matrix: water`, no protein loading; a_w 0.98; pH 6.0. Two programmes as declared: 100 °C for
20 min and 140 °C for 5 min.

**A7 — the composition box.** 200 draws per programme, seed 0, each precursor log-uniform within
its range because the ranges span decades. Every range in `composition_box.yml` is labelled
`stub`; none is sourced. That is the honest state today and is what the module, if built, would
replace. Ranges are deliberately wide and overlapping.

**A8 — the parameter envelope.** For each composition draw in which `top(E) != top(N)`, 50 joint
draws of the core's priors (`kinetic_core.uncertainty.sample_draws`; the sulfur lane's identified
coordinates from the Laplace covariance, the thiol sink's flat direction uniform across its band)
re-run the two arms. The observable multipliers (K_aw, HS-SPME) are the same in both arms and
cancel. The reversal "survives" a draw if the same precursor still wins.

## 7. Outcome (2026-09-13, the declared run: 200 draws per programme, seed 0, 50 envelope draws per reversal)

Artifact: `results/cultivated_tissue_invariance/cultivated_tissue_invariance.{json,md}`; the box that
was swept is echoed beside it. Wall time 33 min at 100 °C, 80 min at 140 °C.

**Verdict on both programmes: indeterminate.** Read off section 3 mechanically:

| clause | threshold | 100 °C / 20 min | 140 °C / 5 min |
|---|---|---|---|
| T3: draws the engine refused | ≥ 50 % | 0 % | 0 % |
| T1: draws where top(E) ≠ top(N) | ≥ 20 % | **40 %** | **38 %** |
| T1: share of disagreements carried by the dominant reversal | ≥ 50 % | 27 % (cysteine over glucose) | 32 % (cysteine over glucose) |
| T1: envelope survival of that reversal | ≥ 80 % | **100 %** | **100 %** |
| T2: top agreement | ≥ 90 % | 60 % | 62 % |
| T2: mean Kendall τ | ≥ 0.75 | 0.35 | 0.33 |

T3 did not fire: the engine answered every draw on the sulfur arm. T2 did not fire: the engine's
ranking is far from the naive one. T1 failed on exactly one clause, concentration. So by the rule
declared in section 3 the result is indeterminate, and section 3 says what that means: it is
resolved by measuring composition, not by rerunning with a narrower box, and it is not a licence to
build the module.

**What the disagreement actually is.** It is not one reversal. It is one direction, spread across
pairs. Over the draws the engine answered, restoring each precursor to its beef level changed the
metric (MFT + FFT + furfural odour-activity) by these medians, relative to the cultivated baseline:

| precursor restored | 100 °C / 20 min | 140 °C / 5 min |
|---|---|---|
| ribose | +124 % | +306 % |
| cysteine | +69 % | +193 % |
| thiamine | +66 % | +9 % |
| glucose | +4 % | +6 % |

The engine puts glucose last in 87 of the 104 draws where glucose and cysteine were both
restorable; the naive ranking puts it last in 25. Glucose is millimolar in beef, so it often carries
a large fractional deficit, and the fractional-deficit rule promotes it. The engine knows that a
hexose reaches the thiols only through the unidentified entry, at a fraction of the pentose's yield,
and demotes it every time. That is the direction, and because which *pair* it shows up in depends on
which precursors happened to be restorable in the draw, no single pair reached half.

**What this does and does not say about the engine.** The engine's disagreement with the naive
rule is robust to the rates (survival 100 % for the dominant reversal on both programmes; 95 % and
99 % averaged over every reversal). It is also, mostly, textbook: pentose over hexose for the meaty
thiols, cysteine as the sulfur donor. A flavour chemist would demote glucose by hand. So the sweep
shows the engine beats *no chemistry*; it was not sharp enough to show it beats *a chemist*, and no
model-free baseline can encode what a chemist knows without becoming a model. That limit was named
in section 1 and it held.

The one row a chemist would not produce by hand is thiamine: worth as much as cysteine at 100 °C,
nearly nothing at 140 °C. That is a kinetic claim, it comes from the thiamine route's own
temperature dependence, and it sits in the lane the core panel scores worst. It is the sweep's only
candidate for a prediction the model adds, and it is exactly the kind that needs the reference-pot
experiment before it is believed.

**What follows, per section 5.** T1 did not fire, so the module is not built. The deliverable is
the composition corpus and its gap map: every range in the box is a stub, and the ranking's
sensitivity to the box is the evidence that composition must be measured first. When a measured
range exists for a precursor, its `stub` label flips to `sourced` with the dossier cited, the test
in `tests/unit/test_cultivated_tissue_invariance.py` that asserts "nothing is sourced today" is
edited to say which, and the sweep is rerun on the narrower box as a *prediction*, not as a
resolution of this verdict.

**A lesson for the next pre-registration of a ranking comparison.** Section 3's concentration
clause assumed a disagreement would look like one reversal. A consistent direction that shows up
in different pairs reads as scattered under that clause and lands as indeterminate. The clause was
right to exist and wrong in shape: the next such pre-registration should ask whether one
*precursor* moves consistently in one direction, not whether one *pair* recurs. Recorded in
`tasks/lessons.md`; not applied retroactively here.

## 8. The composition read (2026-09-14): what exists, what does not

Section 7 said the indeterminate verdict is resolved by measuring composition. The first step is
to find out what has already been measured. One day of reading, through automated fetches of
open-access texts and abstracts (the dossiers say which, and say that none was checked against a
PDF by eye). Seven dossiers under `data/lit/extraction_dossiers/`: `muroya2019`, `bischof2023`,
`kim2024b`, `joo2022`, `koutsidis2008b`, `hwang2026`, `lombardiboccia2005`.

**Beef: every range now has a source.** Cysteine, leucine, IMP and ribose 5-phosphate from one
CE-TOFMS table (Muroya 2019, three steers, 0 to 14 days); glucose, leucine and IMP from one NMR
table (Bischof 2023, fourteen bulls, 0 to 28 days); ribose and thiamine only at second hand (a
cited point plus an abstract's fold-change; two reviews' ranges). Two beef ranges moved against
the stub: cysteine is 0.002 to 0.14 mM, not 0.05 to 0.5 (the stub's centre was above the measured
top), and thiamine is 0.0004 to 0.005 mM, about twofold below the stub.

**Cultured muscle: none of the four rankable precursors has been measured.** Free ribose, free
cysteine, thiamine and free glucose in cultured muscle of any species: no measurement found. What
exists is IMP (two primaries, four decades apart: 0.11 mg/kg in a pig gelatin construct, Kim
2024b; 1.98 mmol/kg in bovine 2D tissue, Joo 2022) and one confounded free-leucine point (Kim
2024b; the scaffold alone carried more leucine than the construct). The one untargeted
metabolomics comparison (Park 2025, chicken) reports fold-changes only. Joo 2022 prints cysteine
as a percent of total amino acids with no absolute total, which does not convert.

So the read closed the beef side and left the side that decides the ranking open. The engine's
ranking depends on the cultivated-side values of ribose, cysteine, thiamine and glucose, and
all four remain stubs. That is the gap map, and it is short.

**The one experiment.** Muroya 2019's CE-TOFMS panel quantified cysteine, ribose 5-phosphate,
IMP and leucine in beef in one run. The same panel on washed cultured bovine myotubes, with free
ribose and free glucose by GC-MS or enzymatic assay and thiamine by thiochrome HPLC on the same
extract, beside a beef sample handled identically, turns every cultivated stub into a sourced
range in one campaign. Three biological replicates; two harvest washes (none, PBS) because the
medium is the obvious confound; one 24-hour, 2 °C ageing arm, because whether the IMP-to-ribose
route runs at all in a construct is the single fact that most changes the ribose stub. Perhaps
thirty samples. No cooking, no GC-O, no aroma work: composition only.

**The prediction run.** The sweep is rerun on the current box (7 sourced, 2 secondary, 5 stub
ranges) to `results/cultivated_tissue_invariance/measured_box_prediction.{json,md}`, with the
same design and the same statistics against the same thresholds, labelled a prediction and not a
verdict. Its outcome is recorded below when it finishes. What it can show: how the beef-side
narrowing moves the rankings. What it cannot show: anything the four cultivated stubs decide.

## 9. The dossiers checked against the PDFs (2026-09-14, later the same day)

Section 8's read was automated and unchecked; its dossiers said so. Four of the PDFs are now on disk
(`koutsidis2008a`, `koutsidis2008b`, `Bischof2023`, `Kim2024`, plus `Dashmaa2026`, which is the review
filed as `hwang2026` after its corresponding author). Every number in those dossiers' section 2 was read
against the paper's own table by eye. Three of the seven dossiers could not be checked (`muroya2019`,
`joo2022`, `lombardiboccia2005`: no PDF on disk); their "Source on disk" lines still say unchecked.

**What the check found.**

- `kim2024b`: every number matched. One misattribution (the construct's 90 % moisture is in the text
  and Fig. 3A, not Table 3).
- `koutsidis2008b`: the abstract-only dossier is replaced by the tables. Beef free ribose 0.25 to
  1.67 mmol/kg over 21 days (n = 16), free glucose 7.33 to 10.3, free cysteine 0.05 to 0.16, IMP 6.27 to
  2.63, ribose 5-phosphate 0.04 flat. The first primary beef ribose values on file.
- `koutsidis2008a` (new): the companion paper, 30 steers at 10 days, two breeds, two diets, with the range
  across individual animals printed: ribose 0.57 to 1.08 mmol/kg, glucose 6.94 to 10.6, cysteine 0.05 to
  0.17, leucine 0.78 to 2.40.
- `bischof2023`: **every number in the automated dossier was wrong.** Table 1 prints glucose as two
  anomer rows (α and β); the fetch returned the α row for one breed and the β row for the other under
  "glucose", and returned numbers for IMP, inosine, hypoxanthine, leucine and methionine that appear
  nowhere in the table. The dossier's own suspicion (repeated columns) was a symptom of a worse fault.
  Corrected: free glucose (α + β) 4.43 ± 2.61 to 10.01 ± 1.42 µmol/g; IMP 1.00 to 3.19; leucine 0.29 to
  1.53.
- `hwang2026`: the one sentence taken was read correctly, but the paper it cites, Aliani et al. 2013
  (Meat Science 94, 55-62), is titled "Post-slaughter changes in ATP metabolites, reducing and
  phosphorylated sugars in **chicken** meat". The review says "beef". The 2.3 mM ribose and 11 mM glucose
  that section 8 used as corners were chicken values, and this review is now cited on no range.

**What moved in the box** (`BOX` in `src/cultivated_tissue_invariance.py`; `STUB_BOX` untouched):

| beef range | section 8 | now | why |
|---|---|---|---|
| ribose | 0.4 to 2.5 mM, secondary | **0.33 to 2.2 mM, sourced** | Koutsidis 2008b Table 2 over 21 d; 2008a inside it |
| glucose | 1.2 to 11 mM, sourced + secondary corner | **2.4 to 15 mM, sourced** | Bischof anomers summed, mean ∓ SD; both Koutsidis inside |
| cysteine | 0.002 to 0.14 mM | **0.002 to 0.23 mM** | Koutsidis 2008a's spread across 30 animals sets the top |
| leucine | 0.35 to 2.4 mM | 0.29 to 3.2 mM | corrected Bischof rows; Koutsidis 2008a's spread |
| IMP, ribose 5-phosphate | unchanged | unchanged | the new primaries sit inside Muroya 2019's spans |
| thiamine | 0.00044 to 0.0049 mM, secondary | unchanged | Lombardi-Boccia 2005 has no PDF on disk |

The box now has 6 beef ranges sourced from a measuring paper's table read by eye, 1 secondary
(thiamine), and the same 5 stubs on the cultivated side. Section 8's line "7 sourced, 2 secondary, 5 stub"
described the box before this check; a prediction run was started on that box and was not committed, so
no artifact of it exists and none is cited.

**Cultured muscle: still nothing.** The four PDFs were searched for free ribose, free cysteine, thiamine
and free glucose in cultured muscle cells or constructs. Kim 2024b mentions glucose only in a cited
C2C12 fasting experiment; the two Koutsidis papers and Bischof 2023 are slaughtered beef; the review has
no cultured-meat content. The four rankable cultivated stubs stay stubs and the test that asserts it is
unchanged. The gap map of section 8 stands.

**What was not done, and why.** No dossier for Aliani 2013: it is a chicken paper and cannot source a
beef range, so the plan's step "point the glucose and ribose notes at aliani2013" was wrong in premise
and was not carried out. The sweep was not rerun on the corrected box; that is a separate decision,
and section 7 says what such a run would and would not show. The pattern behind the Bischof failure is
recorded in `tasks/lessons.md`.

**Addendum, same evening: the three remaining PDFs arrived** (`Muroya2019`, `Joo2022`,
`lombardi-boccia2005`) and were checked the same way.

- `muroya2019`: every number in Table 1 matched. The check added one caveat: the Wagyu muscle blocks
  were 36.5 % fat and 47.7 % moisture, the fat was trimmed by hand before analysis, and the trimmed
  lean's moisture is not printed, so the 75 % used in the conversion is an assumption that could
  understate the mM by a factor under two. Thiamine in this paper is relative content only (Table 3).
- `joo2022`: every number matched, and the methods settled the open question: the amino acids were
  acid-hydrolysed (6 M HCl, 110 °C, 24 h), so Table 1's cysteine is protein cysteine as a share of the
  total. It was rightly not taken. IMP 1.98 mmol/kg in bovine cultured tissue confirmed.
- `lombardiboccia2005`: the excerpt's two numbers were right; the per-cut table is now on file, with the
  DOI confirmed (10.1016/j.jfca.2003.10.007). Beef total thiamine 0.01 to 0.08 ± 0.01 mg/100 g across
  five raw cuts, and **not detected in any cut after cooking**, the behaviour the sulfur lane's thiamine
  route presupposes. The box's beef thiamine range is now sourced from Table 2 alone, 0.00044 to
  0.0040 mM (the unread review that supplied the 0.0049 corner is no longer cited).

Every range the box cites now rests on a table read by eye from a PDF on disk: 7 beef ranges sourced,
0 secondary, 5 cultivated stubs. The cultivated side is unchanged: none of the three papers measures
free ribose, cysteine, thiamine or glucose in cultured muscle (Joo 2022 measured only IMP and
hydrolysed amino acids in its cultured tissue).

## 10. The thiamine row against the literature (2026-10-09)

Section 7 named one row a chemist would not write by hand: restoring thiamine to its beef level is worth
about as much as restoring cysteine at 100 °C / 20 min (+66 % vs +69 %) and almost nothing at 140 °C /
5 min (+9 % vs +193 %). This section asks what the literature already says about it, before any
experiment. Seven PDFs were read by eye (dossiers under `data/lit/extraction_dossiers/`): `brehm2020`
(with its SI), `ramaswamy1990`, `mauri1992`, `mulley1975`, `jhoo2002`, `thomas2014`, `madruga1997`. Two
more were read and set aside: `yang2011` holds thiamine fixed in every run, so it cannot separate
thiamine from cysteine; `schieberle2000` is already in the B16 fit and would be circular.

**What the engine actually runs.** `k_thi_hmp` (thiamine to 5-hydroxy-3-mercapto-2-pentanone) is a rate
fitted at 145 °C, frozen without a recorded uncertainty, and carried to other temperatures by the frozen
lumped formation barrier, 64.08 kJ/mol (`kinetic_core_b9_fit_report.json`; B10's route split was
re-merged). Neither is sampled by the envelope, so section 7's 100 % envelope survival says nothing about
this row. At 64 kJ/mol the engine converts about 1.8 times more thiamine at 140 °C / 5 min than at
100 °C / 20 min; the row's 100 °C advantage therefore comes from the competing sugar route gaining more
from heat, not from thiamine converting faster at 100 °C.

**No paper measures thiamine to MFT, or to its precursor, at two temperatures.** What exists:

| source | conditions | what it measures | apparent Ea of total thiamine loss |
|---|---|---|---|
| Mauri 1992, Table 1 | 0.07 M phosphate, pH 5.5 and 4.0, aw 0.95, 3 µM, 80-100 °C | total loss | 114.6-117.6 kJ/mol (pH 5.5), 121.3 (pH 4.0) |
| Ramaswamy 1990, Table 2 | water, 110-150 °C | total loss | 102.6-118.0 kJ/mol |
| Ramaswamy 1990, Table 2 | water + glucose, glycine, ascorbate | total loss | 71.1-73.4 kJ/mol |
| Brehm 2020, SI Table S4 | 0.1 M phosphate pH 6.5, 100 mM, 80/120/180 °C, one point at 120 min | residual thiamine; HMP, MFT and 3-mercapto-2-pentanone trapped as thiamine thioethers | 22.6-45.2 kJ/mol, derived here; weak (batch-to-batch factor ~2, heat-up unreported, 312.0 vs 412.0 printed for one cell) |
| Mulley 1975, Table 1 | 0.1 M phosphate pH 4.5-6.5, 129.4 °C only | D values | none; k triples from pH 5.5 to 6.5 |

Total loss is a ceiling on the HMP branch, not its value. The engine's 64 kJ/mol sits inside the
measured spread. The case nearest a meat pot (pH 5.5, micromolar thiamine: Mauri) is about 51 kJ/mol
higher; if it applies, thiamine at 100 °C runs about 6 times slower than the engine assumes and most of
the 100 °C advantage goes. If the reactive-mixture value applies (~72), the difference is about 1.3 times
and the row roughly stands. The literature does not choose between them.

**Direction in real meat, at one temperature each.** Thomas 2014 (model hams, 69 °C, peak areas, no
calibration): thiamine is the only additive that raised MFT and the furyl disulfides; cysteine with
xylose or fructose stayed at the control, at every dose, though a clear MFT rise needed thiamine doses
far above native. Madruga 1997 (beef psoas, 140 °C / 30 min, approximate TIC quantification, one portion
per arm): free MFT 13 ng/100 g with IMP, 3 with cysteine, trace with thiamine and in the blank. Both
agree in direction with the row: thiamine matters at low temperature and not at 140 °C, where the ribose
source leads. Neither runs at 100 °C, and the row's claim is about magnitude there.

**A sink the engine lacks.** Jhoo 2002 shows thiamine's pyrimidine fragment trapping MFT as an
almost odourless thioether (MAMP; 110 °C, pH 6.5, molar concentrations). The engine has no such step. At
beef thiamine (0.0004-0.004 mM) it is negligible; in thiamine-dosed reference pots it is not, and it
would bias any test that adds thiamine at millimolar levels.

**Where the row stands.** Direction: weakly supported at both ends. Magnitude at 100 °C: unresolved,
with the closest-to-meat kinetics leaning against it. Estimated probability that the row survives the
reference-pot test: about 0.3 (it was about 0.3 before this read; the kinetics moved it down and the two
meat studies moved it back). Nothing in the box, the engine or the fit changes. A measured barrier for
`k_thi_hmp` would change a core parameter and moves `core_panel_scores.json`, so it belongs to a
pre-registered refit wave, not to this branch.

**The test that decides it.** Thiamine at beef level (micromolar, so Jhoo's sink stays negligible)
against cysteine at beef level, each in a ribose-containing pH 5.6 pot, at 100 °C / 20 min and
140 °C / 5 min, with MFT by stable isotope dilution. Three arms (thiamine, cysteine, both) plus a blank
at each temperature, in triplicate: 24 samples, and it can share the extraction bench of section 8's composition campaign.

## 11. The prediction run on the corrected box (2026-10-09)

Section 8 declared it and section 9 explained why its first attempt was not committed. Run natively, with the
declared design (200 draws per programme, seed 0, 50 envelope draws per reversal), on the box after the
by-eye checks: 7 beef ranges sourced, 5 cultivated stubs. Artifact:
`results/cultivated_tissue_invariance/measured_box_prediction.{json,md}`, box echoed beside it. Wall time
20 min at 100 °C, 56 min at 140 °C. **It is labelled a prediction, and it resolves nothing in section 7.**

| clause | threshold | 100 °C / 20 min | 140 °C / 5 min |
|---|---|---|---|
| T3: draws refused | ≥ 50 % | 0 % | 0 % |
| T1: top(E) ≠ top(N) | ≥ 20 % | 51 % | 46 % |
| T1: share carried by the dominant reversal | ≥ 50 % | **64 %** (ribose over glucose) | **74 %** (ribose over glucose) |
| T1: envelope survival of that reversal | ≥ 80 % | 93 % | 97 % |
| T2: top agreement / mean τ | ≥ 90 % / ≥ 0.75 | 49 % / 0.15 | 54 % / 0.24 |

Read mechanically, both programmes meet T1. **That does not license section 5's module**, for two reasons
declared before this run (sections 7 and 8). A prediction on a box whose cultivated side is still all
stubs cannot resolve a verdict that section 7 said only a composition measurement resolves. And the
reversal is the one section 7 already described: glucose demoted. The beef glucose range rose (2.4 to 15 mM,
sourced) and the beef cysteine range fell (top 0.23 mM, against a cultivated stub reaching 0.5). So
glucose is now usually the largest fractional deficit, cysteine is restorable in only 40 of 200 draws,
and the scattered direction of section 7 collapses onto one pair. The engine is doing what a flavour
chemist would do by hand: pentose over hexose for the meaty thiols. Section 7's limit stands. This shows
the engine beats no chemistry, not that it beats a chemist.

**Median change in the metric on restoring each precursor to its beef level** (over the draws where it was
restorable; section 7's run reproduced exactly from its own artifact):

| precursor | 100 °C, section 7 → now | 140 °C, section 7 → now | draws restorable now |
|---|---|---|---|
| ribose | +124 % → +111 % | +306 % → +287 % | 162 |
| cysteine | +69 % → +31 % | +193 % → +96 % | 40 |
| thiamine | +66 % → +24 % | +9 % → +4 % | 85 |
| glucose | +4 % → +9 % | +6 % → +12 % | 162 |

The thiamine row keeps its shape: comparable to cysteine at 100 °C (+24 % vs +31 %) and close to nothing at
140 °C. Both shrank because the sourced beef cysteine and thiamine ranges are lower than the stubs, so
less is restored. Section 10's caveats on that row apply unchanged.
