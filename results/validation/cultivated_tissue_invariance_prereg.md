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
