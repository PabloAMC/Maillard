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
