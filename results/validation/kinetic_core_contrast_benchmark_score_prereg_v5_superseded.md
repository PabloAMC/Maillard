# Contrast benchmark, scoring step 1: coverage of the held-out contrasts, and the shipped engine on the clean ones (prereg v5)

*Plan W1.1, the yardstick (`docs/guides/STATUS_AND_PLAN.md`, the 2026-11-02 plan; PR #46 lands it on main). v1 and v2 written 2026-11-01 and kept
as `kinetic_core_contrast_benchmark_score_prereg_v1_superseded.md` / `_v2_superseded.md`; their reviews are
`kinetic_core_contrast_benchmark_score_prereg_v1_review.md` (REVISE) and `_v2_review.md` (REVISE). v3 (2026-11-03,
`_v3_superseded.md`) was reviewed in `_v3_review.md` (REVISE); v4 (`_v4_superseded.md`) in `_v4_review.md` (REVISE).
v5, the same day, applies the v4 review and the owner's rule of 2026-11-03 (plan sec. 7: comparisons with tolerances
from the measurements' own noise, no round-number thresholds, no pass marks), BEFORE any engine integration of any
row named here. Generator:
`scripts/generators/score_contrast_benchmark.py`. `--items` (envelope only, no integration) wrote
`kinetic_core_contrast_benchmark_items.{json,md}`, whose numbers are frozen in secs. 2 and 5 below; `--score` writes
`kinetic_core_contrast_benchmark_score_v1.{json,md}`. A diagnostic: it reports comparisons with their intervals,
and nothing ships on it.*

## 0. What the authors had seen

- The split files and their power tables (`docs/reference/contrast_benchmark/README.md` secs. 5, 8 - 11).
- The envelope verdicts (`engine.declare_envelope`, no integration) of every keyed pot x species and of the five atlas
  specs in sec. 2A. No engine number for any row named here was computed or read by the author or either reviewer of
  v2 - v4 (the v2 reviewer saw the two exposed Xu rows, below). To test the v4 review's B8 fix, the v5 author ran
  `engine_value` on one synthetic pot that is in no table (ribose + Cys 20 + 20 mM, 120 C, 30 min, pH 6.2).
- Measured values: every header and row of the three scored tables, and the atlas cells. **The sec. 7 expectations
  are not blind.**
- **Repository exposure:** Xu 2007 MFT / FFT vs pH are directional claims XU07-PH-01 / 02
  (`docs/validation/directional_claims_panel.yml`, status "agree"); the v2 reviewer saw their two score rows. Meynier
  1995's furfural pH fold is a hold-out claim scored since B2 (`kinetic_core_b2_1_holdout_report.md`); Meynier's pH
  drift is a named pH-model control. Hence the `dev_exposed` tag (sec. 2D). Meynier's NF is NOT exposed: the B2 claims
  are MFT / FFT / FUR (Table 3 and the furfural fold), and no claim or fit row uses Meynier's NF. Grosch 1992 Table 2 is
  ledger `digitised_in_dossier`, not an in-model row.

## 1. Questions

- **Q1 (coverage; the headline).** Of the held-out contrasts, how many can the shipped engine answer at all, and why
  not the rest? Reported as two families, never pooled: the 65 v2 TEST atlas cells, and the one-factor pairs of the 15
  keyed held-out tables (which sit outside every frozen split; README secs. 10 - 11).
- **Q2 (score; a small diagnostic).** On the CLEAN held-out pairs it can answer, does the engine get the measured
  direction more often than chance, and is its log-ratio error lower than "no change"?

## 2. Item sets (frozen here)

**Classes (Q1), first that applies:**
1. `inputs_missing`: no printed charge (or an amount printed as null), numeric time, temperature or label pH;
2. `operator_missing`: the observable is a label share, rate, percent remaining, isotopomer %, an activation energy
   or a pH slope (units in `OPERATOR_UNITS`). The unit is the ITEM's: its own rows in arms a and b, operator_missing
   if any of them is in `OPERATOR_UNITS`; with no row in either arm, the units of that species' rows in the pots'
   ledger entry (v3 review B7);
3. `species_absent`: the observed species has no engine target;
4. `precursor_unmapped`: a charged precursor has no engine name;
5. `factor_not_translated`: the arms differ only in a factor this scorer does not translate into the engine's inputs.
   Three have no engine input (lipid and phospholipid on the sulfur lane, where a lipid carrier is a lane conflict;
   a lipid aldehyde, no engine name; reheating); the rest exist on `ProcessSpec` but the keyed tables give them as
   free text or in a form the engine does not take (buffer species / concentration, atmosphere incl. N2, protein,
   matrix) (v3 review S4);
6. `refused_other`: the envelope's own reason, verbatim (first reason of the first refused arm);
7. `answerable`.

**A. Atlas family (65 cells).** Classes 1 - 5 by hand in `ATLAS_COVERAGE` (one reason per cell, read from the cell,
first that applies; `inputs_missing` only where the cell says a charge / time / temperature / pH is not printed,
including `ph_basis: not_measured`, or the system is a food whose charges it does not give). **A cell is answerable if
ANY of its legs is;** such cells carry a spec naming the leg (`ATLAS_SPECS`), and the envelope assigns class 6 or 7. No
hand-classed cell has an answerable leg (checked by the v3 reviewer).
Both tables are hashed (`atlas_table_sha256`) and committed before the score run, at a commit named in the score
artifact. Frozen result:

| class | cells |
|---|---|
| operator_missing | 23 |
| inputs_missing | 17 (v3: 13; + F19#17, F20#13 Shibamoto pH not printed, F22#5 Mizutani, F22#6 Obata pH not measured) |
| species_absent | 12 |
| precursor_unmapped | 7 |
| answerable | 4 (F22#7 Li 2022 H2S vs time; F26#11 Meynier **Table 3** 2-acetylthiazole vs pH; F27#15 Grosch thiamine / ribose / Cys; F29#11 Quan 2020 acrylamide 100 vs 130 C) |
| factor_not_translated | 1 |
| refused_other | 1 (F16#17: the methionine chain did not ship) |

No atlas cell is de-duplicated against a keyed table (v2 review S3): F26#11 is Meynier Table 3, which the keyed file
excludes; F27#15 spans Tables II, III and VIII and the keyed file is Table II. Keyed tables that reproduce atlas cells
are TAGGED (`REPRODUCES_ATLAS`; Grosch carries F27#15 and F27#16), not merged. One paper may sit in different classes in
the two families (Zamora 2015: F30#15 `species_absent`; its keyed items `inputs_missing`, no `initial_mmol_L`): the
families are never pooled.

**B. Keyed family (Q1).** An item is (pot pair, species with a row, an upper limit, an nd / tr entry, or a
`present_not_quantified` / `not_quantified` entry in either arm; that order of precedence when one pot lists a species
twice).
A pair is two pots of one `ledger_id` that differ in EXACTLY ONE factor, where a factor is one composition key
(`initial_mmol_L` or `amounts_umol`, present / absent or its amount) or one of `label_ph`, `t_c`, `minutes`, `buffer`,
`phosphate_mM`, `lipid`, `phospholipid_pct`, `aldehyde`, `post_heat`, `atmosphere`, `protein_mg_per_mL`, `matrix`. The
descriptive fields `aldehyde`, `lipid`, `phospholipid_pct` count only when no composition key differs (Zhao 2017b's
`aldehyde` restates its `V_` aldehyde key). A swap (one key out, another in) is two factors, so not a pair. Envelope:
one `declare_envelope` call per (pot, species), so one target's refusal cannot refuse another's. Frozen result: 236 pot
pairs, 1,424 items:

| class | items |
|---|---|
| species_absent | 739 |
| inputs_missing | 397 |
| refused_other | 96 |
| answerable | 73 |
| factor_not_translated | 54 |
| operator_missing | 53 |
| precursor_unmapped | 12 |

Per table in the items artifact. The v3 review's B7 moved 78 Zhao 2019b items out of `operator_missing` (69 to
`species_absent`, 9 to `factor_not_translated`); the reviewer's count was 68 / 9 with a table-wide unit for the one
item with no row in either arm, which the per-ledger rule sends to `species_absent`. All 73 answerable items are in
the three scored tables (Grosch 18, Meynier 40, Xu 15). A species refused in one arm only is Q1 class 6 and still a Q2
item (sec. 2D), so Q2 has 75 items (Grosch 20).

**C. The frozen name map (v2 review B4).** A YAML key K maps to the FIRST front-door alias (dict order of
`PRECURSOR_ALIASES` / `TARGET_ALIASES`) whose value equals K case-insensitively; no such alias = unmapped. Among
others it maps GSH (glutathione), MET, THI, PENT, Glc, Cys, Gly, LYS (lysine; the envelope then refuses LYS + PENT as
a lane conflict, class 6), MFTD, ACTZ, MPZ, DMP, MTAL, DMDS, FUR, NF, MFT, FFT, and leaves TPP, LEU, CAMA and the `V_`
keys unmapped. **The resolved map printed in the items artifact is authoritative.**

**D. Scored pairs (Q2).** Tables `heldout_grosch1992_tab2`, `heldout_meynier1995_tab1_2_4`, `heldout_xu2007_tab1`.
- **Pair rule (narrow):** two pots of one `ledger_id` differing in exactly one `initial_mmol_L` key or `label_ph`
  (exact for these three tables: none of their pairs differs in an unseen field). Precursors all mapped.
- **Answerable:** the species is in-envelope in at least one arm. A species refused in BOTH arms is not a Q2 item
  (it is Q1 class 6: Meynier MPZ / DMP with glycine, MTAL / DMDS with methionine, every LYS pair, Xu MPZ). A species
  refused in ONE arm stays an item and scores as an engine zero (sec. 4).
- **Exposure:** `dev_exposed` = Meynier FUR (every amine), Xu MFT, Xu FFT. **The readings (sec. 6) are of clean pairs
  only**; exposed pairs are reported beside them as DEVELOPMENT.
- **Methionine:** the Meynier Table 4 pairs (methionine recorded, not charged on the sulfur lane) are Q1 `answerable`,
  flagged; all of them (clean and exposed) reported in their own row, outside the readings.
- **Measurement noise (v5; plan sec. 7):** each printed value carries the log10 scatter its own paper states:
  sigma = hypot(rel / ln 10, read_error), where rel is the relative SD of the REPORTED MEAN. Grosch prints SD and N per
  row: rel = SD / value / sqrt(N). Meynier states "triplicate analyses, CV < 20 %": rel = 0.20 / sqrt(3). Xu states
  "mean CV < 20 %" with n not given: rel = 0.20. `read_error` is the row's (half the last printed digit, in dex).
- **Detected:** both arms printed > 0; m = log10(a / b); sigma_m = hypot(sigma_a, sigma_b). It enters S with the
  weight w = 2 Phi(|m| / sigma_m) - 1: how surely the measured sign is real given the papers' own scatter (0 = a coin,
  1 = certain). No resolution cut-off (v4's |m| > log10 1.4 is withdrawn).
- **Censored:** one arm printed (value v, scatter sigma_v), the other nd or tr with bound L: the LOD for nd (Meynier
  0.01 ug / 10 mg ribose; Xu 0.07 ng/mL; Grosch the row's `upper_limits_ug`), the tr ceiling for tr (Meynier 0.1; Xu
  0.3). The measured sign is the printed arm's; w = max(0, 2 Phi(log10(v / L) / sigma_v) - 1). A printed value at or
  below the other arm's limit has w = 0 and is excluded (v4's v / L > 1.4 is withdrawn). Excluded as before: nd with no
  limit printed, `present_not_quantified`, both arms censored. Censored pairs enter S only.
- **Ladder** = (ledger_id, species, varied factor).

**Frozen Q2 counts.** The first three columns count all Q2 items (exposed and methionine overlap: Meynier Met FUR); the
rest are the decision set (clean, methionine excluded):

| table | all Q2 items | of which exposed | of which Met | E: detected / ladders | S: detected | S: censored | S weight | S ladders |
|---|---|---|---|---|---|---|---|---|
| Grosch | 20 | 0 | 0 | 8 / 4 | 8 | 4 | 11.00 | 8 |
| Meynier | 40 | 20 | 20 | 10 / 1 (NF, glycine) | 10 | 0 | 9.24 | 1 |
| Xu | 15 | 6 | 0 | 5 / 3 | 5 | 4 | 8.72 | 3 |
| **total** | 75 | 26 | 20 | **23 / 8** | **23** | **8** | **28.96** | **12** |

One S item has a refused arm (Grosch FFT no. 1 vs no. 3: FFT is unreachable from thiamine alone; it scores as a
miss, sec. 4). The lowest weights are the near-ties: Meynier NF pH 5.0 vs 5.5 (0.27), Grosch MFT no. 1 vs no. 3
(0.48), Grosch FFT no. 6 vs no. 8 (0.74). Per item in the items artifact.

## 3. The engine (frozen)

- The shipped configuration: the run asserts `Config.from_globals().digest() == Config.shipped().digest()` and records it.
- One `engine.predict(spec, [target])` per (pot, species), target = the sec. 2C alias; the value is
  `require(<alias>)` in ug/L (the prediction is keyed by the requested alias; v4 review B8). Ratios only are used.
- **Process:** `ThermalProgram.isothermal(t_c, minutes)`, `ph = label_ph`,
  `BufferSpec(kind="pyrophosphate", pyrophosphate_mol_l=0.2, declared=True, source=...)` for every pot of the three
  tables (all three print 0.2 M pyrophosphate), `vessel=None`, no `ph_final` (Xu prints one: reported, not used), no
  atmosphere (air). No `ph_drift` is passed, which means the frozen `core_ph_drift()` applies (`engine.py`,
  `_integrate_program`; the `ProcessSpec` docstring saying otherwise is stale). **Declared:** Grosch is an open reflux, and `vessel=None` charges the sealed-pot
  default reservoir.
- **Sensitivity (reported beside the primary run):** the same pots with `BufferSpec(kind="clamped")`, because the
  dynamic pH drifts Meynier's pots while the paper reports <= 0.2 unit. Refusals it meets that the frozen envelope
  verdicts did not are counted and reported.

## 4. Models and per-item values

- **engine:** p = log10(engine_a / engine_b). **no_change:** p = 0.
- **Only one null is scored.** `typical_yield` reduces to log10(L_a / L_b) with L = min(charge), the same in both arms
  of every scored pair, so it equals `no_change`; reported as identical, not as a second test.
- **Engine zeros and refusals:** if the engine returns 0, or refuses (`not pred.answered`), in exactly one arm, the
  item is a sign MISS in S and has error +inf in E. In BOTH arms: p = 0 (a miss in S; error |m| in E). Counts reported.
  In the primary run a refusal is allowed only in an item's frozen `refused_arm`; any other refusal aborts the run.
  **Any other exception aborts the run** (v4 review N2: catching it would turn a bug into a "zero").
- **Numerical resolution:** the engine integrates at rtol 1e-8; |p| < 100 x 1e-8 / ln 10 (4.3e-7 dex) is numerical
  noise and counts as p = 0 (a miss in S). This replaces v4's log10(1.01), a round number.

## 5. Metrics

- **S (direction):** over S items, the WEIGHTED share sum(w x hit) / sum(w), where hit = 1 when sign(p) = the
  measured sign and p != 0; a coin scores 0.5. The unweighted share is reported beside it.
- **E (log-ratio):** over E items, the LOWER median (`statistics.median_low`) of |p - m| for engine and no_change,
  and the median of the per-item difference |m| - |p - m| (positive = the engine's error is smaller). Rank-based, so
  +inf errors never interpolate to NaN. The share within x3 is a descriptive column, never read.
- **Intervals:** ladder-cluster bootstrap, 2,000 resamples, seed 20261101; 90 % intervals = the 5th and 95th
  percentiles with `method="inverted_cdf"`. Reported for S and for the E difference: pooled, per table, the median of
  per-ladder values, and a lab-cluster bootstrap (3 labs, shown to make its degeneracy visible, never read).
- **What the data can detect (context, not a gate; v3 review S1).** Model: in a presence / absence ladder each pair's
  sign agrees with the measured sign with probability q; in a pH ladder the engine's latent value walks over the pots
  in pH order, each adjacent step agreeing with the measured adjacent step with probability q and magnitude Exp(1), and
  the pair signs follow; a refused arm is always a miss. 2,000 simulations per q, each computing the weighted S and its
  interval over the same 2,000 ladder resamples (inverted_cdf, v4 review N1), seed 20261101. The bound variant draws
  once per presence ladder. On the 31 S items in 12 ladders, P(the 90 % interval of the weighted S lies above 0.5):

  | q | 0.5 (false-pass rate) | 0.6 | 0.7 | 0.8 | 0.9 |
  |---|---|---|---|---|---|
  | per pair | 0.047 | 0.169 | 0.427 | 0.714 | 0.918 |
  | one draw per presence ladder | 0.056 | 0.177 | 0.391 | 0.661 | 0.903 |

## 6. How the result is read (frozen sentences; no pass mark)

The result is a comparison with its interval, read by these sentences, written before any integration:
- **Direction** (the 90 % interval of the weighted S on the decision set): above 0.5 -> "better than a coin on these
  held-out directions"; containing 0.5 -> "not distinguishable from a coin on these data (see the detection curve for
  what they could detect)"; below 0.5 -> "worse than a coin on these held-out directions".
- **Log-ratio error** (the 90 % interval of the median difference |m| - |p - m|): above 0 -> "smaller log-ratio error
  than predicting no change"; containing 0 -> "not distinguishable from predicting no change on these data"; below 0 ->
  "larger log-ratio error than predicting no change".
- The detection curve (sec. 5) is printed beside both, so a "not distinguishable" is read against what these 31 items
  could have shown: with the engine right on 70 % of steps, the interval clears 0.5 only ~43 % of the time.
- **Status:** PROVISIONAL (3 laboratories). It is a yardstick reading: nothing ships, no constant moves, no split
  changes. Under plan sec. 7 a change ships on its own written case, and this benchmark is that case's item 7 (side
  effects), never its gate. Q1's coverage tables are not provisional: they are a property of the code.

## 7. Expectations (judgement; not blind, sec. 0)

Recorded earlier and kept (v2, v3): answerable share of the 65 atlas cells <= 10 % (now 4 / 65); exposed pairs agree in
sign, P ~0.75; the Meynier NF ladder above 0.5, P ~0.5; Grosch above 0.5, P ~0.55.

Restated in v5 for the new readings, before any integration: the direction interval lies above 0.5, P ~0.3; it
contains 0.5, P ~0.6; it lies below 0.5, P ~0.1. The log-ratio-error interval lies above 0, P ~0.25; contains 0, P ~0.6;
below 0, P ~0.15. (The engine's pH and amine laws are the open errors E4 / E8, so the per-step agreement is likely
0.6 - 0.7, where the detection curve is 0.17 - 0.43.)

## 8. What this does not do, and the departure from plan step 2

No fit, no engine change, no split change, no table edit. **Departure (v2 review S5):** the 11-01 plan's step 2 and
the 11-02 plan's W1.1 ask for the engine and two nulls on DEVELOPMENT and TEST separately. The
DEVELOPMENT side is the in-panel within-study contrast score that already exists (`contrast_scores.json`,
`scripts/generators/build_contrast_scores.py`); it is cited, not re-run. The keyed held-out rows are not in a frozen
split version, so they are reported as "keyed held-out (not in a frozen split)", not as v2 TEST.

## 9. Revision log

- v1, 2026-11-01: first version.
- v2, 2026-11-01, after the v1 review (all blocking and should-fix items).
- v5, 2026-11-03, after the v4 review and the owner's plan sec. 7 (no round-number thresholds, no pass marks):
  every cut-off replaced by the papers' own noise (S items weighted by 2 Phi(|m| / sigma_m) - 1; censored items by
  their margin over the limit; the inert tolerance tied to the solver's rtol); D1 - D4 replaced by readings of the
  comparisons with their intervals and the detection curve; B8 the prediction read by its alias, a refusal the only
  None, other exceptions abort, an unforeseen refusal aborts the primary run; N1 inverted_cdf in the simulation; N3 E
  per table and by labs, shares within x3; N4 the methionine block holds all Met pairs; N5 new refusals in the clamped
  run reported; N6 the stale 0.722 gone with the gates.
- v4, 2026-11-03, after the v3 review: B7 the item's own unit for class 2 (re-frozen keyed table); S1 the power curve,
  the false-pass rate and D3 read in both directions (a fail inconclusive below 0.6 power at q = 0.7); S2 pH in the
  hand rule (four atlas cells to `inputs_missing`), the any-leg rule, the Zamora note; S3 a raise in both arms is p = 0,
  lower medians and `inverted_cdf` percentiles; S4 the restating clause on composition keys only, class 5 renamed
  `factor_not_translated` with true reasons, the item sources and their precedence; minor: the authoritative map,
  F27#16, 73 vs 75, `core_ph_drift()`, relabelled count headers, noise-limited items reported, and `--score` re-derives
  the items and refuses unless they equal the committed record and the frozen files are committed.
- v3, 2026-11-03, after the v2 review: B4 the name map by rule (GSH, LYS, MPZ / DMP / MTAL / DMDS mapped; one
  envelope call per pot x species); B5 the broad factor rule and the item unit for Q1, the narrow rule kept for Q2;
  B6 zeros / refusals as misses with error +inf, the inert tolerance, the censored inclusion rule; S1 the power model,
  computed and frozen (0.722); S2 D1's interval above 0; S3 no de-duplication, two Q1 families, keyed tables tagged; S4
  the frozen count table; S5 the departure stated; minor: the Config assertion, the BufferSpec keywords, atlas class 6
  / 7 from specs, the Met pairs' class, the hand table's commit, the Grosch vessel.
