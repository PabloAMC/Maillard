# Contrast benchmark, scoring step 1: the shipped engine on the held-out contrasts it can answer (prereg v1)

*Plan step 2 ("the contrast benchmark, built on the split": `docs/guides/STATUS_AND_PLAN.md`, plan from 2026-11-01;
handoff step 3). Written 2026-11-01, BEFORE any engine prediction on any row named here. Generator to be written:
`scripts/generators/score_contrast_benchmark.py`; outputs `results/validation/kinetic_core_contrast_benchmark_score_v1.{json,md}`.*

## 0. What the author had seen when writing this

- The split files and their power tables (`docs/reference/contrast_benchmark/README.md` secs. 5, 8 - 11).
- A coverage inventory (2026-11-01, read-only): every held-out row passed through `engine.declare_envelope`, which
  decides answerability WITHOUT integrating. It gave rows answerable / refused and the refusal reasons, per table.
  No engine number for any row was computed or read.
- The measured values of the scored tables were seen while reading their headers: Meynier 1995's header states that
  furfural, norfuraneol and furfuryl alcohol fall with pH in every table; Grosch 1992 and Xu 2007 values were not read.
  This is a forking-path risk for the Meynier furan expectation in sec. 7, declared here.

## 1. Questions

- **Q1 (coverage).** Of the held-out contrasts, how many can the shipped engine answer at all, and why not the rest?
  This is the headline: a contrast on a factor the engine does not take as an input is answered "no change" by
  construction, which is the no-change null.
- **Q2 (score).** On the held-out contrasts it can answer, does the shipped engine beat both nulls in log-ratio
  error, and does it get the measured direction more often than chance?

## 2. Item set (frozen by this document)

**A. Coverage (Q1).** Every item below is classified by the generator, deterministically, from `declare_envelope`
and the table / cell fields, into exactly one class:
`answerable`, `factor_not_an_input` (the arms differ only in a factor the engine ignores: buffer species or
concentration beyond pH, lipid, metals, O2, atmosphere, reheating), `species_absent`, `operator_missing` (label
shares, rates, percent remaining, isotopomer %), `inputs_missing` (no printed charge, time or temperature),
`refused_other` (the envelope's own reason, recorded verbatim).
- the 65 contrasts in any v2 TEST(E) set (`split_v2_per_change.yml`), each read from its atlas cell;
- every one-factor pot pair of the 15 keyed held-out tables (`data/lit/fit_rows/heldout_*.yml`).
The classification of the atlas cells is done once by hand (they are free text), recorded in the generator as a
reviewed table with one reason per cell, as `build_contrast_split.py` does for its factor tags.

**B. Scored set (Q2).** The one-factor pot pairs, inside one `ledger_id`, of the three keyed tables with answerable
rows:

| table | lab | varied factor | species |
|---|---|---|---|
| `heldout_grosch1992_tab2` | DFA Garching (Grosch) | composition (thiamine vs TPP, +/- Cys, +/- ribose, +/- GSH) | MFT, FFT, MFTD |
| `heldout_meynier1995_tab1_2_4` | Reading (Mottram) | pH 4.5 - 6.5 | FUR, NF |
| `heldout_xu2007_tab1` | China Agricultural Univ. (Gao Y.) | pH 5.6 / 7.0 / 8.0 | MFT, FFT, MFTD, FUR, ACTZ |

- F27#15@grosch1992 is the same data as the Grosch table; it is counted once.
- **Split status, stated with every number:** these rows are "keyed held-out (not in a frozen split)" (README
  secs. 10 - 11). None of the three papers is in v2's `calibration_input_papers` (checked 2026-11-01; `meynier2004`
  is, and is a different paper). The score is a score of the INCUMBENT; it tests no change, so rule (i)'s
  per-change DEV / TEST does not apply, and the calibration-input exclusion does.
- A pair enters only if both rows are answerable, both are finite and > 0 (a `not_detected` / `trace` cell is not a
  row), and the two pots differ in exactly one factor of the table's design. Pairs are all such combinations, as in
  `build_contrast_scores.py`.

## 3. The engine

- The shipped configuration (`Config` defaults at the generator's commit), through the front door
  `engine.predict(FormulationSpec, targets)`, so every refusal is the shipped one. One call per pot.
- Each pot's spec: the table's `initial_mmol_L`, `t_c`, `minutes`, `label_ph`; the buffer passed as the engine
  takes it (sulfur lane: through pH only). A pot the envelope refuses is `refused_other` in A and leaves B.
- Units cancel within a pair (both arms in the same table unit; the engine in ug/L).

## 4. Models

- **engine:** p = log10(engine_a / engine_b).
- **no_change:** p = 0.
- **typical_yield:** p = log10(L_a / L_b), with L the pot's limiting precursor as `fit_thiol_lolo.limiting_precursor`
  defines it, unchanged. It is 0 on every pH pair.

## 5. Metrics

Per pair, measured m = log10(a / b).
- **E (log-ratio error):** median |p - m| over pairs, per model; and the share within x3.
- **S (direction):** on resolved pairs (|m| > log10 1.4, the panel's MEASUREMENT_NOISE_FOLD), the share where
  sign(p) = sign(m) and p != 0. The nulls have no sign on pH pairs, so S is read against chance 0.5.
- **Unit of resampling: the ladder** = (ledger_id, species). Pairs within a ladder share rows, so the 90 % intervals
  come from a ladder-cluster bootstrap (2,000 resamples, seed 20261101). A lab-cluster bootstrap is reported too; with
  3 labs it is degenerate, and it is shown only to make that visible.
- Every number is reported per table, per species and pooled, with the pair, ladder and lab counts beside it.

## 6. Decision rules (fixed)

- **D1 (log-ratio):** the engine beats a null if its pooled median |p - m| is lower AND the 90 % ladder-bootstrap
  interval of the paired difference (null error - engine error, per pair, median) excludes 0.
- **D2 (direction):** the engine's S is above chance if the 90 % ladder-bootstrap interval of S lies above 0.5.
- **D3 (scope):** 3 laboratories. Every verdict is PROVISIONAL. Nothing ships, no constant moves, no split changes
  on this result. Q1's coverage table is not provisional: it is a property of the code.
- A per-table result is reported but decides nothing.

## 7. Expectations (judgement, stated before the run)

- Q1: answerable share of the 65 v2 TEST contrasts <= 10 % (the inventory already showed this; it is not a
  prediction, it is recorded so the reader knows it was known).
- Meynier FUR / NF pH ladders: S >= 0.7, P ~0.6 (forking-path risk, sec. 0).
- Xu 2007 pH ladders: S >= 0.7, P ~0.4 (MFT's pH direction is one of the few signature groups that pass beyond
  chance on the joint fits; the shipped engine is not those fits).
- Grosch composition: S >= 0.7, P ~0.5; the engine beats no_change on E, P ~0.45.
- D1 against no_change, pooled: P ~0.4. D2: P ~0.45.

## 8. What this does not do

No fit, no engine change, no split change, no edit of a table. It does not score DEVELOPMENT (the existing
`contrast_scores.json` covers the panel's in-sample pairs). The keyed tables enter a frozen split only through a new,
dated version (README sec. 3).

## 9. Revision log

- v1, 2026-11-01: first version.
