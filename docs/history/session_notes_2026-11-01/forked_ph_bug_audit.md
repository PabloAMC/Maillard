# Forked-pool pH bug: blast radius audit (origin/main 19cf422f, 2026-11-01)

## 1. Mechanism
- `scripts/generators/fit_joint_deoxyosone_thiol.py:359-379` (`make_raw_many`): `parallel.map_ordered(run_task, ...)` runs each
  pot in a forked worker. `run_task` (`:293`) calls `install(job)` (`:497-501`, `job_config` `:465-494`). That sets
  `ph.trunk_solutes_in_balance = job["trunk_solutes"]`, which is True for every entry-arm trial (`trial`, `:267`).
  The worker returns only the states (`build_submodel_probe.run_pot`, `:130`). The PARENT then calls
  `probe.predict` (`:375`, `build_submodel_probe.py:191-209`) -> `observations.ph_cooled / ph_diluted / base_to_hold`
  -> `ph_state.titratable_inventory`, which reads `ph_state.TRUNK_SOLUTES_IN_BALANCE` at call time (`ph_state.py:650, 699`).
  The parent never installs a job, so that read sees the default False: Gly / AMA / SB / PAMA are left out of the cooled-sample
  balance. The integration itself had them (`sulfur.py:1807`).
- Run serially, the loop predicts in the last task's state (True, which is uniform across a trial), so serial is correct.
- Jacobian points evaluated inside workers run their pots serially, so they are correct too. Only parent-side calls are wrong:
  the hold-out and score-only reports.
- Only one other switch is read after integration, `ph.b51` (cysteine table, ionic strength). No arm varies it.
- Only observables that use the charge balance are affected (pH, diluted pH, NaOH uptake, titratable acidity). In the joint
  fit those are the `ph_acid_submodel` pH rows only, which are hold-out (Cao 2017 `cao2017_gly_x1`) or score-only (Zhou 2023).

## 2. Other parallel generators: not affected
- `fit_ph_acid_submodel`: the parent sets the flag True at `:238` and the tasks set it True at `:136`.
- `predict_ph_acid_holdout`: reports FA and AA only.
- `fit_deoxyosone_submodel`: its arms change no charge-balance switch.
- `build_submodel_probe`: installs no per-task state, so the fork inherits the parent's.
- B2.3 and B64 fits: no per-task state.
- `score_*_signatures` and `feasibility_deoxyosone`: they predict in the worker right after the install.

## 3. Affected tracked artifacts
Twenty entry-arm joint-fit records (.json and .md), all run with MAILLARD_WORKERS 8 or 9 (scratch run scripts):
- e, e_conv, e_fixed
- r3e(_fixed), r3en(_fixed), r3es(_fixed)
- r4f_ea64 / 100 / 133
- r3en_b7anull_w02 / w04 / w02_fixed
- b7d, b7s, b7t, b7st, b7as, b7ass (all _w02)

Affected fields: `holdout_raw_at_optimum`, `holdout_raw_at_start` and `holdout_after_scale_at_optimum`, Cao pH only (7 rows).

Fields that are not affected:
- Zhou score-only pH: identical, because those pots have no trunk solutes.
- theta, cost, Laplace, fitted-row residuals, noise, signatures, chance baselines, Cerny checks, entry checks and LOPO.

Recomputed serially at each recorded optimum (`phbug_recompute.py`, results in `phbug_all.json`). With the solutes off it
reproduces every recorded value to 1e-9. With the solutes on, Cao pH moves by up to 0.16 pH:

| records | recorded RMS | corrected RMS | verdict |
|---|---|---|---|
| R4f x3 | 0.34-0.36 FAIL | 0.21-0.22 | **flips to PASS** |
| E / R3E family | 0.064-0.077 | 0.061-0.081 | mean sign flips (-0.03 -> +0.07); still within |
| R3En / R3Es / B7a | 0.24-0.29 | 0.11-0.16 | within either way |

Text that cites the wrong number:
- `docs/guides/STATUS_AND_PLAN.md:33` ("Cao hold-out pH within 0.06")
- `docs/reference/experiment_register.yml:173` (arm E, "0.06"; corrected RMS 0.080).

EX-030 (R4f) stays rejected: its conclusion rests on other grounds.

## 4. Shipped predictions
Not affected:
- No shipped generator uses the pool with per-task state.
- The flag ships False.
- No joint-fit theta is shipped, and theta is not affected anyway.

## 5. Minimal regeneration
Do not refit. A refit would take about 20 h on 8 workers and would overwrite frozen evidence.
1. Add a re-scoring artifact for the 20 records (new generator; the scratch script is a prototype). About 10 min, serial,
   one process.
2. Correct the two text lines above and annotate EX-030.
3. Optional: regenerate R4f x3 and then e (about 1 h on the fixed pool).

## Separate latent issue (not caused by the fork; identical in serial)
`predict_ph_acid_holdout.jobs_for` computes the pH-stat run's held hot pH in the parent with the solutes OFF: 6.891, which is
what is frozen for arms A, B and C. The fit (`integrate_pot`) would use 6.760 with them ON. This affects the frozen
sub-model 1 pH-stat FA / AA predictions. Not quantified further.
