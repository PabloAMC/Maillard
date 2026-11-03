# Review: kinetic_core_contrast_benchmark_score_prereg.md (v5)

*Independent reviewer (a fresh agent without the author's context), 2026-11-03. Read-only. **Exposure:** the measured
headers and rows of the three scored tables, the items JSON (measured weights, envelope verdicts), the engine source
for `predict`, `require`, the integrator tolerances and the pH fixed point. Ran: a re-derivation of `items_payload()`
(envelope only; identical to the committed record after stripping volatile fields), `declare_envelope` on the clamped
specs (all 53 pot x species give the primary run's verdicts), the 11 unit tests (pass). Did not run `predict` or
`--score`; did not open `core_directional_scores.*` or `contrast_scores.json`. The v5 document is kept as
`kinetic_core_contrast_benchmark_score_prereg_v5_superseded.md`; line numbers are v5's.*

## v4 review items

| item | status | where |
|---|---|---|
| B8 wrong key | **RESOLVED.** `predict` read by the requested alias; a refusal is the only None; a None outside the frozen `refused_arm` aborts the primary run. Engine confirmed: `predict` returns an unanswered prediction rather than raising (`engine.py:2724-2726`); concentrations keyed by compound name (`engine.py:2170-2184`). Monkeypatched test. | code 609-619, 754, 763-765; prereg 155-156, 172-175; test 75-84 |
| N1 percentile method | **RESOLVED.** `inverted_cdf` in the simulation; curve re-frozen. | code 543; prereg 192-199 |
| N2 other exceptions | **RESOLVED.** Not caught: any other exception aborts. | prereg 175; code 609-619 |
| N3 E by table / lab | **RESOLVED in substance.** Residue: the per-table and per-ladder values of the E *difference* (non-blocking 8). | code 711-718; prereg 186-188 |
| N4 methionine block | **RESOLVED.** `methionine_all` holds every Met pair. | code 769, 772; prereg 122-123 |
| N5 clamped refusals | **RESOLVED.** Counted, exempt from the abort, rendered. | code 763-765, 773, 924-925; prereg 163-165 |
| N6 stale 0.722 | **RESOLVED.** Only in the v3 log entry. | prereg 255 |

## Noise model against the tables

- **Grosch** prints mean +/- SD over N trials per row; SD / value / sqrt(N) is the SD of the mean (MFT no. 1 vs no. 3:
  sigma_m 0.069, w 0.478).
- **Meynier** states "triplicate analyses, CV < 20 % except some low-concentration compounds", each pot a mean of 3;
  0.20 / sqrt(3) gives sigma_m 0.071. 20 % is an upper bound, so the weights are conservative; the exception clause is
  not applied, and "analyses" may be repeat injections rather than repeat reactions.
- **Xu** states "mean CV < 20 %, max 47 %", n not stated; 0.20 gives sigma_m 0.123, conservative as single experiments,
  but 0.20 is a mean over compounds, not each compound's own.
- **Detection limits** match the headers (Meynier 0.01 / 0.1, Xu 0.07 / 0.3, Grosch per row).
- **Frozen numbers** match the JSON: 20 / 40 / 15; E 8/4, 10/1, 5/3; S 23 + 8 = 31 items in 12 ladders; weight sum
  28.96; near-tie weights 0.268 / 0.478 / 0.743; both detection curves.
- **The weighting is sound.** Under a flat prior Phi(|m| / sigma) is P(measured sign = true sign); for an engine
  agreeing with the true sign with probability a, E[hit] = 0.5 + (a - 0.5) w, so weighting by w is the locally most
  powerful linear statistic. The censored weight is a lower bound on the same probability. Under a coin engine
  E[S_w] = 0.5 except for the one forced miss (Grosch FFT no. 1 vs no. 3), which lowers it to 0.484; the simulated
  false-pass rate already includes it, so reading against 0.5 is conservative.

**Round numbers that remain.** Deciding: the 90 % interval level (selects the reading sentence); `INERT_DEX = 100 x
rtol` (decides p = 0; wrong basis, item 2). Describing only: `BAND_DEX = log10 3`, `N_BOOT` / `N_SIM` = 2000 and the
seed, the agreement grid, the sec. 7 probabilities.

## Blocking

None.

## Non-blocking (1 - 3 before `--score`)

1. A frozen sentence differs: code says "(see the **power** curve ...)", prereg says "**detection** curve". Fix the
   code and the `render_score` label.
2. The inert tolerance is tied to the wrong solver. On the sulfur lane the floor is the pH fixed point (tolerance 1e-3
   pH; residual floors at 3e-5 - 6e-4 pH, `config.py:192-194`, `sulfur.py:1672`); the repo measured cross-platform
   drift of a log10 step at ~3.3e-6 dex (`provenance.py:132-139`); both well above 4.3e-7 dex. Practical risk small
   (every varied factor is a reactive precursor or pH). Fix: name the pH fixed point; take each item's floor as
   |p(shipped) - p(tighter `ph.fixed_point.tolerance`)|, or at least the measured drift.
3. `PhFixedPointNotConverged` only warns and the scorer neither records nor stops on it. Record a per-prediction flag.
4. Normal weights for SDs from N = 2 - 4 overstate certainty (MFT no. 1 vs no. 3: 0.48 -> ~0.42 under t). Use
   Student-t with Welch-Satterthwaite dof for rows with `sd` / `n`; re-freeze.
5. Meynier's exception and Xu's mean CV (max 47 %): say both; optionally a description with the maximum CV.
6. Say the coin's expectation on this set is 0.484 (the forced miss).
7. E ignores the measurement noise: say why (both rivals face the same m; the ladder bootstrap carries the scatter;
   within-study ratios cancel lab offsets, so H0 does not enter).
8. N3 residue: add `per_table_difference` and `difference_median_of_ladders`, or narrow the text.
9. Grosch's 8 S ladders share pots (no. 6, no. 8): the ladder bootstrap is somewhat narrow; note it (description).
10. Minor: `_committed_at` should check the items MD; the revision log is out of order.

**Verdict: APPROVE WITH CHANGES.**
