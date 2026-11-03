# Contrast benchmark, scoring step 1: coverage of the held-out contrasts, and the shipped engine on the clean ones (prereg v2)

*Plan step 2 (`docs/guides/STATUS_AND_PLAN.md`, plan from 2026-11-01; handoff step 3). v1 written 2026-11-01 and kept
as `kinetic_core_contrast_benchmark_score_prereg_v1_superseded.md`; v2 the same day after the independent review
(`kinetic_core_contrast_benchmark_score_prereg_review.md`, verdict REVISE). Written BEFORE any engine prediction on any
row named here. Generator: `scripts/generators/score_contrast_benchmark.py`; outputs
`results/validation/kinetic_core_contrast_benchmark_score_v1.{json,md}`.*

## 0. What the author had seen

- The split files and their power tables (`docs/reference/contrast_benchmark/README.md` secs. 5, 8 - 11).
- A read-only coverage inventory: held-out rows passed through `engine.declare_envelope` (no integration), giving
  answerable / refused rows and reasons per table. No engine number for a row was computed or read in this session.
- Measured values: the headers of all three scored tables print measured values (Meynier: the furans fall with pH;
  Grosch lines 4 - 6, 19 - 20; Xu lines 4 - 5), and atlas cells F27#15 / F27#16 print theirs. **The sec. 7
  expectations are therefore not blind.**
- **Repository exposure (review B1):** the engine has been checked against parts of these data before: Xu 2007 MFT /
  FFT vs pH are directional claims XU07-PH-01 / 02 (`docs/validation/directional_claims_panel.yml:3154-3200`, status
  "agree"); Meynier 1995's furfural pH fold is a hold-out claim scored since B2
  (`kinetic_core_b2_1_holdout_report.md:36`); Meynier's pH drift is a named pH-model control. Those values were not
  re-read for this document. Grosch 1992 Table 2 is ledger `digitised_in_dossier`, not an in-model row.

## 1. Questions

- **Q1 (coverage; the headline).** Of the held-out contrasts, how many can the shipped engine answer at all, and why
  not the rest? A contrast on a factor the engine does not take as an input is answered "no change" by construction.
- **Q2 (score; a small diagnostic).** On the CLEAN held-out pairs it can answer, does the engine get the measured
  direction more often than chance, and is its log-ratio error lower than "no change"?

## 2. Item set (frozen here)

**A. Coverage (Q1).** Each item gets exactly one class, by this precedence (first that applies):
1. `inputs_missing`: no printed charge, time or temperature;
2. `operator_missing`: label share, rate, percent remaining, isotopomer %;
3. `species_absent`: the observed species has no engine target;
4. `precursor_unmapped`: a charged precursor has no engine name (e.g. TPP, GSH, LYS beside a pentose);
5. `factor_not_an_input`: the arms differ only in a factor the engine ignores (buffer species / concentration beyond
   pH, lipid, metals, O2 / atmosphere, reheating);
6. `refused_other`: the envelope's own reason, verbatim;
7. `answerable`.

Items: the 65 contrasts in any v2 TEST(E) set, and every one-factor pot pair (sec. 2B's pair rule) of the 15 keyed
held-out tables. The 65 atlas cells are classified BY HAND into a reviewed table inside the generator
(`ATLAS_COVERAGE`, one reason per cell). That table is committed, with its sha256 printed in the artifact, before the
generator's first scoring run; the score run asserts the hash.

**B. Scored pairs (Q2).** Tables `heldout_grosch1992_tab2`, `heldout_meynier1995_tab1_2_4`, `heldout_xu2007_tab1`.
- **Pair rule.** Two pots of one `ledger_id` that differ in exactly one factor, where a factor is one
  `initial_mmol_L` key (present / absent or its amount) or `label_ph`. A swap (one key out, another in: Grosch no. 6 vs
  no. 7, Cys for GSH) is not a pair. All such combinations, per species.
- **Exposure tag (review B1).** `dev_exposed` = Meynier FUR (every amine), Xu MFT and Xu FFT. Everything else is
  `clean`. **D1 and D2 are decided on clean pairs only**; exposed pairs are reported beside them as DEVELOPMENT.
- **Detected pairs:** both rows finite and > 0. **Censored pairs (review S2, secondary):** one arm not detected; the
  measured direction is known and |m| >= log10(value / LOD), with the table's per-pot LOD where printed; scored for
  sign only.
- **Engine zeros (review S4):** a pair where the engine returns 0 (or a refusal) for an arm stays in the item count
  and is reported as `engine_zero` / `engine_refused`; it is never dropped silently.
- **De-duplication:** F27#15 / F27#16@grosch1992 and F26#11@meynier1995 are the same data as the tables; counted once.
- **Split status:** "keyed held-out (not in a frozen split)" (README secs. 10 - 11); none of the three papers is in v2's
  `calibration_input_papers`. DFA Garching and Reading are DEVELOPMENT laboratories under the lab rule (v1); the score
  is of the incumbent and tests no change, so this is reported, not used to exclude.

## 3. The engine (frozen)

- The shipped `Config` at the generator's commit; its `Config.digest()` is recorded. One `engine.predict` call per pot.
- **Name map (review B2), frozen in the generator as `PRECURSOR_MAP` / `TARGET_MAP`:** PENT -> ribose; THI ->
  thiamine; Cys -> cysteine; Gly -> glycine; MET -> methionine; FUR -> furfural; NF -> norfuraneol; MFT, FFT, MFTD,
  ACTZ -> the engine's own target names. TPP, GSH and LYS stay unmapped (class 4): a lysine + ribose pot would be a
  lane conflict, and the generator does not choose a lane for it. Every mapping is checked against
  `PRECURSOR_ALIASES` / `TARGET_ALIASES` at start; a missing name fails the run instead of being guessed.
- **Methionine (review B2):** the sulfur lane records methionine but does not charge it (`engine.py:2549-2557`), so
  Meynier Table 4 runs as ribose without an amine. Its pairs are reported in their own row and excluded from D1 / D2.
- **Process:** `t_c`, `minutes`, `label_ph`; buffer `BufferSpec` pyrophosphate 0.2 M for every pot of the three
  tables (all three print 0.2 M pyrophosphate); `vessel=None`; no `ph_final` (Xu prints one: reported, not used).
- **Sensitivity (decides nothing):** the same pots with the pH clamped at `label_ph` (`BufferSpec(kind="clamped")`),
  because the dynamic pH drifts Meynier's pots while the paper reports <= 0.2 unit.

## 4. Models

- **engine:** p = log10(engine_a / engine_b).
- **no_change:** p = 0.
- **Only one null is scored (review B3).** `typical_yield` reduces to log10(L_a / L_b) with L = min(charge), which is
  the same in both arms of every scored pair, so it equals `no_change` here. It is reported as identical, not as a
  second test. A null with a sign needs DEVELOPMENT contrasts of the same factor step; none is built in step 1.

## 5. Metrics

- **S (direction):** on resolved detected pairs (|m| > log10 1.4) plus censored pairs, the share with sign(p) = sign(m)
  and p != 0; chance 0.5.
- **E (log-ratio):** median |p - m| over detected pairs, engine and no_change; share within x3.
- **Ladder** = (ledger_id, species, varied factor). 90 % intervals by ladder-cluster bootstrap (2,000 resamples, seed
  20261101). Reported: the pooled value, the median of per-ladder values, and each table separately with its
  pair / ladder / lab counts. A lab-cluster bootstrap is shown to make its degeneracy (3 labs) visible.
- **Power line (review S1), computed before scoring from the item counts only:** the probability that D2 passes when
  the true per-ladder sign agreement is 0.8, under the ladder bootstrap, by simulation with the actual ladder sizes.

## 6. Decision rules (fixed)

- **D1:** the engine beats no_change on E if, on clean detected pairs, its pooled median |p - m| is lower AND the 90 %
  ladder-bootstrap interval of the per-pair difference (no_change error - engine error, median) excludes 0.
- **D2:** S is above chance if, on clean pairs, the 90 % ladder-bootstrap interval of S lies above 0.5 AND at least
  two of the three tables have S > 0.5 (so one table cannot carry the pool).
- **D3:** if the power line is < 0.6, D1 / D2 are reported as UNDERPOWERED and decide nothing.
- **D4:** whatever D1 - D3 say, the verdict is PROVISIONAL (3 laboratories). Nothing ships, no constant moves, no split
  changes. Q1's coverage table is not provisional: it is a property of the code.

## 7. Expectations (judgement; not blind, sec. 0)

- Q1: answerable share of the 65 v2 TEST contrasts <= 10 % (known from the inventory; recorded, not predicted).
- Power line < 0.6 (UNDERPOWERED), P ~0.6: ~10 clean ladders.
- Clean S > 0.5 in >= 2 tables, P ~0.45. D1 passes, P ~0.3.
- Exposed pairs (Meynier FUR, Xu MFT / FFT) agree in sign, P ~0.75 (the panel has scored them before).

## 8. What this does not do

No fit, no engine change, no split change, no table edit. It does not score DEVELOPMENT beyond the exposed pairs.

## 9. Revision log

- v1, 2026-11-01: first version.
- v2, 2026-11-01, after the review (all blocking and should-fix items):
  B1 the exposure tag, D1 / D2 on clean pairs; B2 the frozen name map, the pyrophosphate BufferSpec, vessel and
  ph_final, methionine excluded, the clamped-pH sensitivity; B3 one null, said plainly; S1 the power line and the
  two-of-three-tables rule; S2 censored pairs for sign; S3 sec. 0 rewritten; S4 the pair rule, the swap exclusion, the
  de-duplication, engine zeros kept; minor: the class precedence, the hashed hand table, per-ladder medians.
