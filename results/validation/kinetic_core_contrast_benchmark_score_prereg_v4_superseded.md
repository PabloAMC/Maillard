# Contrast benchmark, scoring step 1: coverage of the held-out contrasts, and the shipped engine on the clean ones (prereg v4)

*Active order item 1 (`docs/guides/STATUS_AND_PLAN.md`, "Active order from here"). v1 and v2 written 2026-11-01 and kept
as `kinetic_core_contrast_benchmark_score_prereg_v1_superseded.md` / `_v2_superseded.md`; their reviews are
`kinetic_core_contrast_benchmark_score_prereg_v1_review.md` (REVISE) and `_v2_review.md` (REVISE). v3 (2026-11-03,
`_v3_superseded.md`) was reviewed in `_v3_review.md` (REVISE); v4 written the same day after that review, BEFORE any
engine integration of any row named here. Generator:
`scripts/generators/score_contrast_benchmark.py`. `--items` (envelope only, no integration) wrote
`kinetic_core_contrast_benchmark_items.{json,md}`, whose numbers are frozen in secs. 2 and 5 below; `--score` writes
`kinetic_core_contrast_benchmark_score_v1.{json,md}`. A diagnostic under STRATEGY_2026-10-29.md sec. 9: scores
diagnose, nothing ships on them.*

## 0. What the authors had seen

- The split files and their power tables (`docs/reference/contrast_benchmark/README.md` secs. 5, 8 - 11).
- The envelope verdicts (`engine.declare_envelope`, no integration) of every keyed pot x species and of the five atlas
  specs in sec. 2A. No engine number for any row named here was computed or read by the author or either reviewer of
  v2 / v3 (the v2 reviewer saw the two exposed Xu rows, below).
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
- **Exposure:** `dev_exposed` = Meynier FUR (every amine), Xu MFT, Xu FFT. **D1 and D2 are decided on clean pairs
  only**; exposed pairs are reported beside them as DEVELOPMENT.
- **Methionine:** the Meynier Table 4 pairs (methionine recorded, not charged on the sulfur lane) are Q1 `answerable`,
  flagged; reported in their own row and excluded from D1 / D2.
- **Detected:** both arms printed > 0; m = log10(a / b). **Resolved:** |m| > log10 1.4.
- **Censored (v2 review B6):** one arm printed, the other nd or tr. The bound is the LOD for nd (Meynier 0.01 ug /
  10 mg ribose; Xu 0.07 ng/mL; Grosch the row's `upper_limits_ug`) and the tr ceiling for tr (Meynier 0.1; Xu 0.3).
  Included iff printed value / bound > 1.4; the measured sign is the printed arm's. Excluded: nd with no limit printed,
  `present_not_quantified`, both arms censored. Censored pairs enter S only.
- **Ladder** = (ledger_id, species, varied factor).

- **Noise-limited (reported only, v3 review m6):** a detected item where both arms print an SD and value +/- SD
  overlap (Grosch FFT no. 6 vs no. 8). The rule is not changed.

**Frozen Q2 counts.** The first three columns count all Q2 items (exposed and methionine overlap: Meynier Met FUR); the
rest are the decision set (clean, methionine excluded):

| table | all Q2 items | of which exposed | of which Met | E: detected / ladders | S: resolved detected | S: censored | S ladders |
|---|---|---|---|---|---|---|---|
| Grosch | 20 | 0 | 0 | 8 / 4 | 7 | 2 | 6 |
| Meynier | 40 | 20 | 20 | 10 / 1 (NF, glycine) | 9 | 0 | 1 |
| Xu | 15 | 6 | 0 | 5 / 3 | 5 | 4 | 3 |
| **total** | 75 | 26 | 20 | **23 / 8** | **21** | **6** | **10** |

No S item in the decision set has a refused arm. Meynier contributes one ladder (9 of 27 S items).

## 3. The engine (frozen)

- The shipped configuration: the run asserts `Config.from_globals().digest() == Config.shipped().digest()` and records it.
- One `engine.predict(spec, [target])` per (pot, species), target = the sec. 2C alias; the value is
  `require(<engine key>)` in ug/L (ratios only are used).
- **Process:** `ThermalProgram.isothermal(t_c, minutes)`, `ph = label_ph`,
  `BufferSpec(kind="pyrophosphate", pyrophosphate_mol_l=0.2, declared=True, source=...)` for every pot of the three
  tables (all three print 0.2 M pyrophosphate), `vessel=None`, no `ph_final` (Xu prints one: reported, not used), no
  atmosphere (air). No `ph_drift` is passed, which means the frozen `core_ph_drift()` applies (`engine.py`,
  `_integrate_program`; the `ProcessSpec` docstring saying otherwise is stale). **Declared:** Grosch is an open reflux, and `vessel=None` charges the sealed-pot
  default reservoir.
- **Sensitivity (decides nothing):** the same pots with `BufferSpec(kind="clamped")`, because the dynamic pH drifts
  Meynier's pots while the paper reports <= 0.2 unit.

## 4. Models and per-item values

- **engine:** p = log10(engine_a / engine_b). **no_change:** p = 0.
- **Only one null is scored.** `typical_yield` reduces to log10(L_a / L_b) with L = min(charge), the same in both arms
  of every scored pair, so it equals `no_change`; reported as identical, not as a second test.
- **Engine zeros and refusals (v2 review B6, v3 review S3):** if the engine returns 0, or raises, in exactly one arm,
  the item is a sign MISS in S and has error +inf in E. Zero or raise in BOTH arms: p = 0 (a miss in S; error |m| in E).
  Counts reported.
- **Inert tolerance:** |p| < log10(1.01) is p = 0 (a miss in S).

## 5. Metrics

- **S (direction):** over S items, the share with sign(p) = measured sign and p != 0; chance 0.5.
- **E (log-ratio):** over E items, the LOWER median (`statistics.median_low`) of |p - m| for engine and no_change;
  share within x3. Rank-based, so +inf errors never interpolate to NaN (v3 review S3).
- **Bootstrap:** ladder-cluster, 2,000 resamples, seed 20261101; 90 % intervals = the 5th and 95th percentiles with
  `method="inverted_cdf"`. Reported: the
  pooled value, the median of per-ladder values, each table separately, and a lab-cluster bootstrap (3 labs, shown to
  make its degeneracy visible).
- **Power (v2 review S1, v3 review S1), computed from the items alone, frozen.** Model: in a presence / absence
  ladder each pair's sign agrees with the measured sign with probability q; in a pH ladder the engine's latent value
  walks over the pots in pH order, each adjacent step agreeing with the measured adjacent step (from the printed
  values, nd = below every printed value) with probability q and magnitude Exp(1), and the pair signs follow; a
  refused arm is always a miss. 2,000 simulations per q, each with D2 applied in full over the same 2,000 ladder
  resamples, seed 20261101. The bound variant draws once per presence ladder (its pairs share one engine law). On the
  27 S items in 10 ladders, P(D2 passes):

  | q | 0.5 (false-pass rate) | 0.6 | 0.7 (likely, sec. 7) | 0.8 | 0.9 |
  |---|---|---|---|---|---|
  | per pair | 0.070 | 0.215 | **0.456** | **0.726** | 0.925 |
  | one draw per presence ladder | 0.066 | 0.185 | 0.390 | 0.641 | 0.880 |

  The v3 reviewer's independent implementation gave 0.06 / 0.21 / 0.45 / 0.73 / 0.93.

## 6. Decision rules (fixed)

- **D1:** the engine beats no_change on E if, on the E items, its pooled median |p - m| is lower AND the 90 %
  ladder-bootstrap interval of the median per-item difference (|m| - |p - m|) lies ABOVE 0.
- **D2:** S is above chance if the 90 % ladder-bootstrap interval of S lies above 0.5 AND at least two of the three
  tables have S > 0.5.
- **D3 (restated, v3 review S1):** a PASS of D1 or D2 decides when the power at q = 0.8 is >= 0.6 (it is 0.726; the
  false-pass rate is 0.070). A FAIL decides only when the power at the likely q = 0.7 is >= 0.6; it is 0.456, so **a
  fail is reported as INCONCLUSIVE**, not as evidence that the engine's directions are at chance. D1 has no power line
  of its own and is read through D2's.
- **D4:** whatever D1 - D3 say, the verdict is PROVISIONAL (3 laboratories). Nothing ships, no constant moves, no split
  changes. Q1's coverage tables are not provisional: they are a property of the code.

## 7. Expectations (judgement; not blind, sec. 0)

v2's, kept as written: answerable share of the 65 atlas cells <= 10 % (now 4 / 65, frozen in sec. 2A); the power line
< 0.6, P ~0.6 (it came out 0.722 once GSH was mapped and censored pairs entered); clean S > 0.5 in >= 2 tables, P ~0.45;
D1 passes, P ~0.3; exposed pairs agree in sign, P ~0.75.

v3 additions, written before any integration (unchanged in v4): D2 passes, P ~0.35 (power 0.72 assumes 0.8 agreement; the engine's
pH and amine laws are the open errors E4 / E8, so the true agreement is likely lower); the Meynier NF ladder has S > 0.5,
P ~0.5; Grosch S > 0.5, P ~0.55 (thiamine and cysteine directions are what the sulfur lane was built on).

## 8. What this does not do, and the departure from plan step 2

No fit, no engine change, no split change, no table edit. **Departure (v2 review S5):** plan step 2
(`STATUS_AND_PLAN.md`, "Order" item 2) asks for the engine and both nulls on DEVELOPMENT and TEST separately. The
DEVELOPMENT side is the in-panel within-study contrast score that already exists (`contrast_scores.json`,
`scripts/generators/build_contrast_scores.py`); it is cited, not re-run. The keyed held-out rows are not in a frozen
split version, so they are reported as "keyed held-out (not in a frozen split)", not as v2 TEST.

## 9. Revision log

- v1, 2026-11-01: first version.
- v2, 2026-11-01, after the v1 review (all blocking and should-fix items).
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
