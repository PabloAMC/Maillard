# Review: kinetic_core_contrast_benchmark_score_prereg.md (v4)

*Narrow confirmation review (a fresh agent without the author's context), 2026-11-03. Read-only. **Exposure:** no
engine prediction of any kind; the items JSON summary fields, the atlas `ph_basis` tags and the engine source for
`require` / `predict`; `core_directional_scores.*` not opened; `predict` and `--score` not run. The v4 document is
kept as `kinetic_core_contrast_benchmark_score_prereg_v4_superseded.md`; line numbers are v4's and the generator's at
review time.*

## v3 review items

| v3 item | status | where |
|---|---|---|
| B7 class 2 by the item's unit | **RESOLVED.** `item_units` takes the units of the item's own rows in arms a and b, falling back to the ledger. The JSON gives 739 / 397 / 96 / 73 / 54 / 53 / 12 = 1,424 over 236 pairs, the same as the prereg table. | prereg 40-43, 89-101; code 281-289, 381-383 |
| S1 power curve, D3 | **RESOLVED.** Both curves are in the JSON (per pair 0.0695 / 0.215 / 0.456 / 0.726 / 0.9255; one draw 0.0655 / ... / 0.88). D3 is read in both directions in `decisions`. One leftover under Non-blocking (percentile method). | prereg 176-189, 197-200; code 225-228, 533-590, 705-709 |
| S2 atlas hand rule | **RESOLVED.** pH in class 1, the any-leg rule, the Zamora note. No hand cell outside `inputs_missing` has `ph_basis: not_measured`; F22#5 and F22#6 are `not_measured`. Atlas counts match the JSON (sha matches the code). | prereg 54-58, 65, 74-76; code 73-77, 95, 99, 105-106 |
| S3 raise in both arms; NaN | **RESOLVED.** Both arms zero gives p = 0, `median_low`, and `inverted_cdf` (numpy 2.5.3 returns -inf, not NaN, where the default method gives NaN). | prereg 162-163, 170-171; code 612-614, 638, 651, 655-657 |
| S4 restating clause / class 5 / item sources | **RESOLVED.** The restating clause checks `composition_diff` only; `NOT_TRANSLATED` gives the true reasons; precedence row > upper limit > nd > tr > pnq / `not_quantified`. | prereg 46-50, 78-85; code 199-209, 254-257, 260-278 |
| m1 - m7 | **RESOLVED** (m4: `engine.py:2618` confirmed; m6: Grosch has 1 S item flagged; m7: `strip_volatile` drops `git` and `generated_on`, input sha256s are content hashes). | prereg 105-109, 74, 102-103, 151-152, 131-139, 128-129, 231-232 |

## New blocking

- **B8. `engine_value` asks for the wrong key, so every prediction silently becomes a "zero"** (code 603-606; prereg
  sec. 3 line 147). `pred.require(engine.TARGET_ALIASES[tg])` passes the engine key (`"MFT"`), but
  `CorePrediction.concentrations_ug_per_l` is keyed by the requested target string (`engine.py` ~2815). Confirmed for
  all six Q2 species (ACTZ, FFT, FUR, MFT, MFTD, NF): `require` raises `OutOfEnvelope` ("not among the requested
  targets"), which is caught and returned as None for every arm, so p = 0 everywhere, S = 0, engine E = no_change, and
  D1 / D2 read "FAIL: INCONCLUSIVE" with no error raised. **Fix:** `if not pred.answered: return None; return
  pred.require(tg)`; in `score_payload`, a value may be None only for the item's frozen `refused_arm` (primary run),
  otherwise abort; prereg sec. 3 to "`require(<alias>)`".

## Non-blocking

- **N1.** The power simulation (`d2_passes`, code 522) uses the default linear `nanpercentile`; the score path and
  sec. 5 use `inverted_cdf`. Pass `method="inverted_cdf"`, re-run `--items`, re-freeze sec. 5.
- **N2.** Prereg sec. 4 counts any "raise" as a miss; after the B8 fix only a refusal is None. State that any other
  exception aborts the run; do not catch it.
- **N3.** Sec. 5 asks for per-table and lab-cluster results; `score_block` gives both for S only. Add E `per_table`
  and `ci90_lab` (and a within-x3 share), or narrow the wording.
- **N4.** The methionine block requires `exposure == "clean"` and leaves out the exposed Met FUR items, while the
  prereg puts all Met pairs in their own row. Drop the filter, or say Met FUR is reported under exposed only.
- **N5.** The clamped sensitivity can meet refusals the frozen record never saw. Report them as a count and exempt
  them from the B8 guard.
- **N6.** Sec. 7 says "came out 0.722"; the frozen value is 0.726.

**Verdict: REVISE.**
