# Review: kinetic_core_contrast_benchmark_score_prereg.md (v1)

## Blocking

**B1. Two of the three tables are not held out from the incumbent.**
- **Xu 2007.** These exact pots (ribose + Cys 50/50, 0.2 M pyrophosphate, 140 C / 1 h, pH 5.6 / 7.0 / 8.0) are already the directional claims XU07-PH-01 (MFT) and XU07-PH-02 (FFT) (`docs/validation/directional_claims_panel.yml:3154-3200`). Engine values for them are committed in `core_directional_scores.json`, with status "agree".
- **Meynier 1995.** `meynier_FUR_pH_fold` ("furfural falls 15-49x") has been scored in every B-panel since B2 and passes (`generate_kinetic_core_b2_2_holdout.py:437-458`; `kinetic_core_b8_panel.md:72,150`). PH-01/02 use its protocol. Its pH drift is a named control on the pH model (`kinetic_core_thiol_lolo_prereg.md:176`). Its Table 3 is `in_model_check` (`sulfur.yml:3614`).
- **Labs.** Under the lab rule, `dfa_garching` and `reading_mottram_parker` are DEVELOPMENT (`in_v1_test: false`, `split_v2_per_change.yml:1086-1116`).

So sec. 0 ("no engine number ... computed") and sec. 2 line 49 are wrong for the repository.

*Fix:* tag each pair `clean` or `dev_exposed` (Meynier FUR/NF; Xu MFT/FFT). Decide D1 and D2 on clean pairs only, and report the exposed pairs as DEVELOPMENT, as plan step 2 asks.

**B2. Sec. 3 as written cannot be run, and the parts it leaves open change the answer.**
- **Names.** The YAML keys are not front-door names. `PENT`, `THI`, `MET`, `LYS`, `TPP`, `FUR`, `NF`, `MFTD` and `ACTZ` are not in `PRECURSOR_ALIASES` / `TARGET_ALIASES` (I checked), so `declare_envelope` refuses them as UNMAPPED (`engine.py:1141-1169`). The key-to-name map has to be frozen, and it decides lanes:
  - "lysine" becomes `Lys` (acrylamide lane) next to ribose (sulfur lane): a LANE CONFLICT.
  - `MET` is "recorded_not_charged" on the sulfur lane (`engine.py:2549-2557`), so Meynier Table 4 is scored as ribose with no amine. Class it `factor_not_an_input`, or report it separately.
- **Buffer.** "Sulfur lane: through pH only" is false. The sulfur lane runs the dynamic pH with `buffer_spec`, and the default is unbuffered `none` (`engine.py:2596-2618`).
- **Freeze:**
  - `BufferSpec("pyrophosphate", 0.2)` for every pot (as PH-01/02 and XU07 do);
  - `vessel=None`;
  - `ph_final` not passed (Xu prints one);
  - the engine commit or digest, not just "Config defaults at the generator's commit".
- **Sensitivity.** Declare a non-deciding clamped-pH sensitivity run. The engine drifts Meynier's pH by -0.4 to -0.8 units, while the paper measured <= 0.2.

**B3. `typical_yield` is the same as `no_change` on every scored pair.** `limiting_precursor` takes the minimum of the charge (`fit_thiol_lolo.py:289-291`). That is THI 0.06 in every Grosch pot, PENT 35 in every Meynier pot and 50 in every Xu pot, so p = 0 everywhere. "Beats both nulls" is therefore one test. The function also cannot be called "unchanged": it needs a `fit_joint` row key, and `fit_rows.load` refuses these tables.

*Fix:* add a null that can take a sign, for example the median measured log-ratio of the same factor step (+Cys, +ribose, +GSH, per pH unit) on DEVELOPMENT contrasts. Otherwise say plainly that there is one null.

## Should-fix

**S1. Power (rule ii).** My count, if the inventory's answerable rows hold, is about 50 pairs in about 10 ladders:

| table | pairs | ladders | detail |
|---|---|---|---|
| Meynier | 36 | 4 | Gly FUR 6, Gly NF 10, Met FUR 10, Met NF 10; Lys refused |
| Grosch | 8 | 2 | MFT 7, FFT 1; TPP has no rows |
| Xu | 6 | 4 | MFT, FUR, MFTD 1 each; ACTZ 3; FFT 0 |

The pooled numbers are about 72 % Meynier. *Fix:* add a simulated-sign power line before the run, and pool with each table weighted equally, or require the per-table directions to agree.

**S2. Outcome-selected drop.** The "both > 0" rule (line 52) removes the strongest contrasts: Xu's FFT 22.09 -> nd, MFT and FUR at pH 8 (E8's own consumer), and Grosch's TPP. *Fix:* add a secondary censored-sign item with m >= log10(value / LOD), using `upper_limits_ug`, Xu's ~0.07 ng/mL and Meynier's ~0.01 ug.

**S3. Sec. 0 disclosure.** Measured values were in front of the author. The Grosch header (lines 4-6, 19-20) and the Xu header (lines 4-5) print them, and so do atlas cells F27#15 and F27#16 (4.7 / 33 ug, x7), which the author hand-classifies for Q1. The sec. 7 expectations for these tables are therefore not blind.

**S4. Pair rule.**
- Define a "factor" as one `initial_mmol_L` key or `label_ph`.
- State that Grosch 6 vs 7 (Cys swapped for GSH, which the header lists as a target contrast) is excluded.
- Deduplicate F27#16 and F26#11@meynier1995 as well as F27#15.
- Fix what happens when the engine gives 0 for an arm. `build_contrast_scores.py:43` silently drops such pairs, which is a selection that depends on the engine.

## Minor

- The classes in sec. 2A need a precedence order. An unmapped precursor (TPP) fits no named class.
- Freeze or hash the hand classification of the 65 atlas cells before any run.
- The E (log-ratio) metric's pooled median hides ladder sizes of 1 to 10; report a median of ladder medians beside it.

**Verdict: REVISE.**
