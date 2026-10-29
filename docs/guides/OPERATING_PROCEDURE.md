# Operating procedure: plan, execute, revise; and what to do when data do not fit

*Written 2026-10-29 (strategy M1; owner: "how to set up a plan, execute, and know when to revise; how to notice
when the data does not fit a model and figure out the wisest possible directions"). Every rule here is the answer
to a failure in `tasks/lessons.md`, cited in brackets. It is short on purpose; change it when it stops helping.*

## 1. The plan

- **One living plan:** `docs/guides/STRATEGY_2026-10-29.md` for the order and the reasons, mirrored by the task list
  the owner sees. One step is active at a time. A finding that belongs to a later step is written into that step,
  not acted on now [plan clarity, 2026-10-29].
- **Every step states**, before it starts: the question; the decision its result will change (no decision, no
  step); the deliverable; the test for "done"; a budget (hours of work, core-hours); and its expected result.

## 2. The step cycle

1. **Frame.** Question, decision, done-test, budget, expectation (one paragraph).
2. **Look up before building.** The measurement index for every direct measurement of the species and steps
   involved; the conflicts register; the record for the same idea tried before [B59b's audit missed Zhai 2019
   Fig. 3A; "Schieberle as fit rows" was B16].
3. **Cheap probe.** The change evaluated with nothing refitted: rows, flux budget, hold-outs, on one core in
   minutes [15 refits launched before a five-minute probe].
4. **Design and short pre-registration** (about a page): data and their conditions checked against the paper's own
   text [Zhai 2021 rows were pH 5.5, labelled pH 7]; what is free and why it is identified; arms; checks; the rule
   for what goes forward.
5. **Design review** by an independent agent, given the code and the data, not my conclusions, before any compute
   [the stage-1 review found the missing pH laws after the compute].
6. **Run** through `src/kinetic_core/fitting.py` (two starts; identification; leave-one-laboratory-out with the
   predictive check) [false convergence; a relative finite-difference step], niced, within the core cap. Agents
   write their findings to a file as they go [two agents lost].
7. **Results.** Numbers printed by the generator from its JSON; prose only interprets [a sigma shift written as a
   fold]. Results review by an independent agent.
8. **Decide and record.** Forward / carry as alternatives / redesign / stop, with the reason; commit at the gate.

**Data guards (2026-10-29, owner: "ensure we do not make the same mistake again").**
- Fit rows come only from keyed tables (`data/lit/fit_rows/*.yml`), one row per measurement, each pointing at its
  measurement-ledger id. Conditions, quantification class, series, laboratory and duplicate links come from the
  ledger, not from code. No generator types a measured value.
- A row enters a fit only when its ledger entry is verified against the PAPER (not the dossier) by two independent
  readers, with conditions, units, quantification, duplicates and a sample of values checked
  (`data/lit/source_verification.yml`; `scripts/ci/source_verification_gate.py`).
- Engine switches are set inside a function or context manager, never at import
  (`scripts/ci/import_time_switch_gate.py`).
- Every constant cites its evidence by a key that resolves to a dossier and ledger entry (M4 builds the reaction
  evidence table from those keys).

**A sub-model is done only when (mechanism checks; a fit to its own rows proves little):**
1. it reproduces its OWN species, intermediates included, across laboratories, under the measurement model;
2. the diagnostics show no structural flag (one-sign groups, trends with T / pH / time / concentration), or each
   flag is explained by the misfit protocol;
3. no fitted constant sits on a bound or far from an independent measurement of the same step, and leave-one-
   laboratory-out moves no constant beyond its sigma without being reported;
4. element balances (C, N, S) close where the source measures totals; isotope-labelling shares are reproduced
   where they exist;
5. one prediction OUTSIDE the fitted conditions (another laboratory, temperature, pH, concentration or time scale)
   was frozen before it was scored, and is reported whatever it says;
6. where the literature offers a rival mechanism, it is fitted as an arm, and the conditions where the two diverge
   are ranked by the experiment-value tool.

## 3. When to revise the plan (triggers)

Stop the step, write what changed in a few lines, update the strategy and the task list, then continue, when:
- a probe or a review contradicts the step's premise;
- a result falls outside the expectation written in step 1 (an unexpected result is information about the plan,
  not only about the model);
- a data error or duplicate is found in inputs the plan relies on;
- the done-test fails twice, or the budget is exceeded twofold;
- the owner asks.
At each phase gate, re-read the strategy against the north-star metrics (leave-one-laboratory-out error on
contrasts and levels, interval coverage, evidence grounding, blind tests) and ask whether the next phase is still
the most valuable one.

## 4. When the data do not fit (the misfit protocol)

**Notice it.** The diagnostics report (M7) shows the residuals by series, laboratory, species, temperature, pH,
time and concentration. A misfit is structural, not noise, when: it has one sign across laboratories or conditions;
it correlates with a condition variable; a fitted constant sits on a bound or far from its measured analogue;
leave-one-group-out finds predictive tension; or a held-out series has the right level but the wrong shape.

**Then, in this order** (cheapest and most frequent cause first):
1. **Integrity.** Re-read the source: conditions, units, the analyte's definition and standard, duplicates,
   digitisation. In this repository this has been the most common cause.
2. **Conditions.** An unrecorded or mis-specified condition (heat-up, pH path, buffer, oxygen, headspace,
   concentration) becomes an input or a latent condition with a prior from facts, never a laboratory term.
3. **Measurement.** Method class, calibration, internal-standard artefacts: the measurement model's job.
4. **Mechanism.** Candidate mechanisms from the hypothesis layer, the literature and the conflicts register; each
   must predict something beyond the misfit and is tested on data outside it; prefer the one that explains most
   with measured constants.
5. **Unresolved.** Record it (conflicts or choices register), carry the alternatives, rank the experiment that
   would settle it.

**Never:** tune an unidentified constant to absorb a misfit; drop the inconvenient data; average conflicting
laboratories by row count; judge a change by whether a score improved.

## 5. Choosing the wisest direction

Rank candidate next steps by: (what it changes for the purpose: flux, predictions, decisions) × (probability it is
right) ÷ (cost). In practice: reading and digitising before fitting; measured constants before fitted ones; a
change that removes a fitted constant before one that adds one; one change at a time; and the cheapest test that
could prove a direction wrong before investing in it.
