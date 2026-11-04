# H0 amendment A2: three pair fits re-run because of wall-clock pot timeouts

- Date: 2026-11-04. **Written before any H0 number was computed; no pair residual was read.** The only fields read
  from the pair records were optimiser and engine fields: `complete`, `converged_all`, `status`, `history` lengths,
  `evaluations`, `leak_prone` counts, `freeze_commit` and `failed_pot_runs`.
- Scope: three pair records, the scorer's record precondition, and this note. No estimator, threshold or verdict
  rule changes. The prereg, the runner and its settings are unchanged.

## What happened

The ten pair fits ran on a 16 GB machine whose memory was mostly held by other applications (a virtual machine
held 10 GB). Swap filled and the fits thrashed. To stop the thrashing, the session paused two fits with SIGSTOP for
part of the run and resumed them later.

Every pot integration runs under a 120 s wall-clock alarm (`build_submodel_probe.run_pot_limited`, via
`fit_joint_deoxyosone_thiol.POT_SECONDS`). Wall-clock time kept passing while a fit was paused or starved of memory,
so integrations "timed out" and counted as failed pot runs. A failed run enters the objective as a penalty residual
(`fitting.Problem.failed_residual`). The three affected records therefore depended on machine state, not on the
chemistry:

| pair | failed pot runs | machine state during the run |
|---|---|---|
| dalian_polytechnic_chen_d + dfa_hofmann_schieberle | 55 | ran through the worst thrashing (load average up to 180) |
| dalian_polytechnic_chen_d + ochanomizu_murata | 117 | paused twice by SIGSTOP |
| dfa_hofmann_schieberle + ochanomizu_murata | 50 | paused once by SIGSTOP |

The other seven records report 0 failed pot runs, as did every EX-056 fold.

## The rule, fixed before any residual was read

- A pair record with `failed_pot_runs > 0` is re-run from scratch, with the same settings, unpaused.
- The re-run comes from a worktree at the freeze commit 1fb63f5a, because the H0 branch now carries amendment A1 and
  the runner refuses a changed frozen file. Its record therefore carries the same freeze commit as the other seven.
- This is not "a non-converged fit re-run with other settings", which prereg sec. 7.4 forbids: all three had
  converged, the settings are identical, and the criterion is an engine failure count, not a residual or its effect
  on tau.
- The three contaminated records are kept, unread, outside the repository (the session scratchpad), as the record of
  what was replaced.
- The scorer's precondition now refuses any pair record with `failed_pot_runs > 0` (a unit test pins it), so a
  timeout-contaminated record can never be scored silently.

## Outcome of the re-runs

The three re-runs ran 17:31 - 20:00 on 2026-11-04 from a detached worktree at 1fb63f5a, three at a time, unpaused,
after the virtual machine was closed. Read from optimiser and engine fields only:

| pair | failed pot runs | optimiser | wall time |
|---|---|---|---|
| dalian_polytechnic_chen_d + dfa_hofmann_schieberle | 0 | warm 2 (203 ev), reset 2 | 8895 s |
| dalian_polytechnic_chen_d + ochanomizu_murata | 0 | warm 3 (228 ev) | 3159 s |
| dfa_hofmann_schieberle + ochanomizu_murata | 0 | warm 3 (262 ev), reset 3 | 8841 s |

All ten pair records now report 0 failed pot runs, every fit converged (status > 0), and every record carries freeze
commit 1fb63f5a. The timeouts were the pauses and the memory starvation, not the chemistry.

## Lesson

Never pause (SIGSTOP) a fit whose engine runs under a wall-clock limit, and never run fits on a machine short of
memory. Reduce concurrency at launch instead.
