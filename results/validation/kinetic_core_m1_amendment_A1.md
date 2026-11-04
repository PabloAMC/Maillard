# M1 amendment A1: deterministic exploratory pot failures

- Date: 2026-11-04, 22:45. **Written before any decisive number** (no Laplace for C1, no decisive score): the
  only fields read were optimiser and engine fields (status, cost, evaluations, failed pot runs and their pots, the
  per-segment failure counts) and the optimum's own residual vector, to recompute its cost. No decisive datum was read.
- Scope: prereg v3.3 sec. 4.3's rule "a record with failed_pot_runs > 0 is re-run, same settings". The rule was
  written for failures caused by machine state (H0, amendment A2: wall-clock timeouts under SIGSTOP and memory
  thrash). It cannot clear a deterministic failure: a re-run reproduces it exactly, and the chain stops at the
  Laplace step forever.

## What happened

C1's starts 1 and 2 each logged one failed pot run: `b60:kang_freecys_100` exceeded the 120 s wall-clock alarm.

| record | failed (total / in segments) | cost | 0.5 r.r at the optimum | evaluations |
|---|---|---|---|---|
| C1 s1 (first run, kept as `_failed_run1`) | 1 / 1 | 154.32736219217875 | 154.3273621921789 | 226 |
| C1 s1 (same-settings re-run, quiet machine, load 4 - 10) | 1 / 1 | 154.32736219217875 | 154.3273621921789 | 226 |
| C1 s2 (first run, kept as `_failed_run1`) | 1 / 1 | 154.326490521694 | 154.3264905216941 | 202 |

The re-run reproduced start 1 to the last digit (the same optimiser path, the same failing pot), on a nearly idle
machine: the failure is a deterministic stiff point the optimiser explores, not machine state. Every failure is inside
the optimiser's segments; the re-evaluation at the optimum added none, and the optimum's residuals give its cost. A
failed exploratory point enters the objective as a penalty the optimiser steers away from; it cannot alter a clean
optimum's predictions.

## The rule, as amended

A fit record with failed pot runs is accepted only if all three hold:
1. every failure happened inside the optimiser's segments (the optimum's own evaluation was clean);
2. the cost recomputed from the optimum's residuals equals the record's cost (relative 1e-9);
3. a same-settings re-run reproduced it exactly: the same cost (relative 1e-12), evaluations, failure count and
   failed pots.

Otherwise the record is re-run, as before (a machine-state failure does not reproduce). Implemented as
`failure_accepted` in `fit_m1_supply_combinations.py`, used by `--fit`, `--laplace` and `--score`, unit-tested. The
acceptance and its reason are printed and recorded. The incumbent's starts have no failure; the rule applies to C1
only. No estimator, threshold or verdict rule changes.
