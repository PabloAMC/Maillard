# Review: kinetic_core_contrast_benchmark_score_prereg.md (v3)

*Independent reviewer (a fresh agent without the author's context), 2026-11-03. Read-only: nothing edited, staged or
committed; two scratch scripts re-derived the counts through `docker_maillard.sh`; `--items` was not re-run.
**Exposure:** no engine number of any kind; the envelope verdicts and refusal reasons in the items JSON / MD, and the
statement text (not the scores) of `directional_claims_panel.yml`; `core_directional_scores.*` not opened. The v3
document is kept as `kinetic_core_contrast_benchmark_score_prereg_v3_superseded.md`; line numbers below are v3's and
the generator's at review time.*

## v2 items: resolved?

| v2 item | status | where |
|---|---|---|
| B4 name map | **Resolved.** Map by rule, case-insensitive, first alias. GSH, LYS, MPZ / DMP / MTAL / DMDS, MFTD and ACTZ map. LYS + PENT is class 6 ("LANE CONFLICT"), checked in the JSON. One `declare_envelope` per (pot, species). The prose list is incomplete (minor m1). | prereg 88-92; code 169-182, 389-394 |
| B5 broad pair rule / item unit | **Resolved in the text.** The code departs from it in two places: B7 (unit per item), and a latent bug in the restating clause (S4). | prereg 67-73, 94-96; code 187-243 |
| B6 zeros / tolerance / censoring | **Resolved.** `_measured` follows the rule. Two edge cases are still open: a raise in both arms, and a median with infinite values (S3). | prereg 104-108, 140-142; code 440-457 |
| S1 power model | **Resolved as a model** (specified, implemented, reproduced). **The way D3 uses it misleads** (S1 below). | prereg 151-156; code 508-551 |
| S2 D1 interval above 0 | **Resolved.** | prereg 160-161 |
| S3 no de-duplication, two families | **Resolved.** The Grosch tag misses F27#16 (m2). | prereg 63-65; code 59-67 |
| S4 frozen counts | **Resolved.** Re-derived from the YAML: Grosch 20 / 0 / 0 / 8-4 / 7 / 2 / 6; Meynier 40 / 20 / 20 / 10-1 / 9 / 0 / 1; Xu 15 / 6 / 0 / 5-3 / 5 / 4 / 3; 27 S items in 10 ladders. All match. | prereg 111-120 |
| S5 departure from plan step 2 | **Resolved.** | prereg 181-185 |
| minor: Config, BufferSpec, atlas specs, Met class, commit, vessel | **Resolved.** One stale point under m4: `ph_drift=None` uses `core_ph_drift()`. | prereg 47-51, 102-103, 124-131; code 273-276, 632-633 |

Also checked and holding: 236 pot pairs recomputed independently (same set); the narrow and broad rules agree on the
three scored tables; every censored bound and ratio (Grosch 1.27, 0.80, 1.35 and 8.0; Xu 17.6-56); Meynier NF has no
claim and no fit-row use; Xu 2008's duplicate pH 5.6 columns are only `digitised_in_dossier`; class precedence in
`classify_pair` follows sec. 2.

## Blocking

**B7. Class 2 is decided by a table-wide unit, not the item's unit. This mislabels 77 frozen Q1 items.**
- `table_unit` (code 265-267) takes `sorted(units)[0]` over all of a species' rows in the table, and `classify_pair`
  (359-361) applies it to every pair.
- Zhao 2019b prints MFT, FFT, ACTZ and 23 other species in both `percent_of_isotopomers` (Table 2) and
  `relative_peak_area` (Table 1). `percent_...` sorts first, so all become `operator_missing`, including pairs whose
  rows in both arms are relative peak areas. Of Zhao 2019b's 131 `operator_missing` items, 77 have only
  `relative_peak_area` rows: 9 (MFT / FFT / ACTZ, +/- fat) are class 5, `lipid`; 68 are `species_absent`.
- The frozen keyed table should read: operator_missing 54 (not 131), species_absent 738 (not 670),
  factor_not_an_input 54 (not 45). Q1 is the headline and "a property of the code", so text and code must agree
  before the freeze. Larralde and Sedehi also have mixed units; their class is unaffected.
- **Fix:** take the unit from the item's own rows in arms a and b; write the rule for disagreeing arms or no row
  (e.g. "operator_missing if any arm's row is in `OPERATOR_UNITS`"). Re-run `--items` and re-freeze sec. 2B.

## Should-fix

**S1. D3 is met only at an effect size the authors themselves doubt (prereg 151-156, 164-165, 175-176).** 0.722 is
P(D2 passes) at per-step agreement 0.8; sec. 7 expects lower. Independent implementation, 1,000 simulations each, same
bootstrap and seed:

| agreement | P(D2 passes), as specified | P, one draw per presence ladder |
|---|---|---|
| 0.5 | 0.06 (the false-pass rate) | 0.06 |
| 0.6 | 0.21 | 0.17 |
| 0.7 | 0.45 | 0.40 |
| 0.8 | 0.73 | 0.67 |
| 0.9 | 0.93 | 0.90 |

The test is correctly sized, but a D2 failure carries little information below ~0.75 agreement, where sec. 7 puts
its prior. Presence ladders draw each pair independently (code 518-521); Grosch's Cys ladder is three backgrounds of
one engine law, so this flatters power by ~0.06. **Fix:** print the curve and the false-pass rate; restate D3 so D2
decides only in the pass direction and a fail is "inconclusive" unless power >= 0.6 at an agreement the authors call
likely; optionally the one-draw variant as a bound.

**S2. Hand rule for atlas cells (prereg 47-49).** Class 1 includes "label pH" but the hand clause limits
`inputs_missing` to charge / time / temperature: F19#17 and F20#13 (Shibamoto 1977) carry `ph_note: not printed`, so
class 1 should precede `species_absent` (inputs_missing 15, species_absent 12), or say pH is excluded. Multi-leg cells
(F26#11, F29#11, F27#15) were made answerable by the answerable leg with no stated rule. **Fix:** "a cell is answerable
if any leg is; the spec names the leg"; check no hand-classed cell has an answerable leg (none found). F30#15 and the
keyed Zamora items give one paper different classes: acceptable (families not pooled); say so.

**S3. Sec. 4 leaves two scoring edge cases open.** An engine raise in BOTH arms is undefined: choose now (p = 0, or
dropped and reported). D1's median and bootstrap percentile with +/-inf can return NaN (interpolation between -inf and
+inf, 0 x inf). **Fix:** a rank-based estimator (`method="inverted_cdf"`, or the lower median), stated.

**S4. Code differs from the text with no effect today.**
- Restating clause (code 238-242): `not (f in RESTATING and diff)` checks `diff`, which by then holds scalar factors
  too; a `buffer` + `lipid` pair would count as one-factor. No effect today (same 236 pairs). **Fix:** composition
  differences computed separately.
- `NOT_AN_INPUT` (code 193-196; prereg 42-43) says the engine ignores what it does not: `ProcessSpec` takes
  `atmosphere` ("argon" / "air" / "air_cu"), `protein_g_per_l`, `matrix` and buffer concentration; lipid with the
  sulfur lane is a lane conflict, not "ignored". **Fix:** narrow class 5, or call it "factor not translated here".
- `species_arms` (code 246-262) reads `present_not_quantified` as an item source (the text lists row / nd / tr /
  upper limit only); ignores Grosch `not_quantified`; `setdefault` lets `not_detected` hide a later `upper_limits_ug`
  bound. **Fix:** align text or code.

## Minor

- **m1.** Prereg 90-92 omits mapped DMHF, PZ, GO, MGO, TDG; "every `V_` key" misses non-`V_` unmapped keys (DMP26,
  FFTD, MP2P / MP3P, MFT_FFT_mixed_disulfide, THIOPH_*, PENT as a target). Say the artifact's map is authoritative.
- **m2.** `REPRODUCES_ATLAS` for Grosch should add F27#16 (keyed no2-no5 is TPP +/- Cys).
- **m3.** One-arm refusals are Q1 class 6 and Q2 items, so 73 Q1-answerable != 75 Q2 items (Grosch 18 vs 20); say so.
- **m4.** Prereg 129: "no `ph_drift`" means the frozen `core_ph_drift()` (`engine.py:2616-2618`; the `ProcessSpec`
  docstring is stale). Say so, or the clamped sensitivity looks vacuous.
- **m5.** "Decision set" headers count all Q2 items in the Q2 / exposed / Met columns, and exposed and Met overlap
  (10 Met FUR). Relabel.
- **m6.** Grosch FFT no6 vs no8 (1.2 +/- 1.1, n = 3 vs 2.3 +/- 1.0, n = 4) is a resolved S item within the printed
  scatter. Report it as noise-limited; do not change the rule now.
- **m7.** The items artifact was generated on a dirty tree. Have `--score` re-derive the items and assert equality with
  the committed JSON, not only the hash.

**Verdict: REVISE.**
