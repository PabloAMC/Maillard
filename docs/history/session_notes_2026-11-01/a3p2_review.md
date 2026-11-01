# Review: PR #36 (a3-config-phase2, 82a81359), independent

Verdict: **MERGE**, with three follow-ups (F1-F3). I found no change to default-mode numbers. Everything below is either latent or applies only with the flag on.

## What I ran (worktree at origin/a3-config-phase2, niced, one process at a time)
- tests/unit/test_engine_config.py and test_guard_gates.py with MAILLARD_WORKERS=2: 41 passed. test_cysteine_route_calibration, test_parallel and test_submodel_probe_time_limit: 20 passed.
- import_time_switch_gate passes. build_wave_manifest --check reports 96 frozen and none changed.
- The values the migrated generators set all construct and install without refusal: B21 (both branches), B40, B41, B59a's alternatives, CR-02's Q10s.
- applied() gives every binding back the very object it held.
- Fit coordinate bounds sit inside the field bands, with equality at the edges. fitting._forward_jacobian steps inward at an upper bound, so a fit never trips a band.
- No script or test sets a build key, apart from the one test that expects the refusal.
- Not run: the artifact freshness gate (I stopped it after 10 minutes because the machine was loaded) and the full suite.

## Findings, by severity

**F1 MEDIUM: records already on main are wrong, and nothing flags them.** This concerns parallel.py and the lessons entry. The diagnosis is correct. On main, a serial run predicted every row in the last task's state, and a forked run predicted in the parent's state. Only job_config sets ph.trunk_solutes_in_balance, so the parent never had it on.

The switches the predictions read are ph.b51 and ph.trunk_solutes_in_balance. Both are fixed per arm, not per trial, so the branch's per-pot run_under() equals main's serial result. Serial is the right truth: the observation operator should use the integration's own charge balance.

The fitted entry rows (ph_acid_davidek_acids: AA, PENT, FA amounts) do not read the charge balance, so the fits themselves are unaffected. But the hold-out and score-only pH rows in any record from a forked run are wrong, and records do not log the worker count (s042_b7st_w02 has no such field). Add a backlog item: re-score the hold-out pH of the entry-arm records serially.

**F2 MEDIUM (latent): config.update() now imports about 20 modules before a ship-rule child sets its fit report.** Config.from_globals() imports every binding's module, including uncertainty, panel, fit_targets and ea_class_prior. uncertainty.py:102 and :107 compute LAPLACE_PATH and PROFILE_PATH at import from engine._B2_FIT_REPORT.

In the `_child` of b51_1, b58, b61, b62 and b63 ship rules, update now runs before `engine._B2_FIT_REPORT = report`. Uncertainty therefore binds to the shipped report, not the candidate's. Today no `_child` reads uncertainty (score_panel, directional and declare_envelope do not), so the outputs are unchanged.

The claim "same place, same timing" is still not literally true, and the WAVES.md digest evidence cannot see state outside the Config. Fix: have from_globals read only modules already in sys.modules, or move the report swap ahead of update.

**F3 MEDIUM (flag off): the persistent pool carries and resets only the keyed state.** These are not carried: engine._B2_FIT_REPORT (12 assignments across 10 scripts, invisible to the coverage scan because of its leading underscore), the ea_class_prior lru_caches (filled from core_parameters under whichever task ran first, then kept across tasks), and engine._REPORT_CACHE.

A spawned worker gets the defaults where a fork inherits the parent's values. A picklable fn mapped after a fit-report swap would compute silently with the shipped report. No current call site does this: the only picklable fns are the four sub-model fits and predict_ph_acid_holdout, and the rest are closures that fall back to fork. Before the flag defaults on, either refuse persistent mode when _B2_FIT_REPORT is not the default, or carry it in the Config.

**L1 LOW: the allowlist test catches public names only.** A new public UPPERCASE call-time global fails until it is keyed or listed, so the test does catch a new switch. It is blind to:
- private globals (_B2_FIT_REPORT);
- callables (b60_verification monkeypatches engine.integrate_sulfur);
- getattr by string;
- reads from src/*.py outside kinetic_core.

The list was drafted from the scan, and its 'numeric' category is broad.

**L2 LOW: every update() replaces all mapping keys with fresh dicts.** I confirmed that the identities of THREE_DEOXY_EXIT_PH and TRUNK_PH_EXPONENT_OVERRIDE change. No in-place mutation or identity read exists, so this is harmless. A Config stores a mapping as a tuple of items in insertion order, so equal dicts in a different order give different digests. That only produces false negatives.

**L3 LOW (intended): CR-02 now changes lipid's q10_used and its decades.** I did not re-run choice_robustness.json, so the claim that the artifact is unchanged is still unverified.

## Frozen generators
The edit follows the manifest's own rule: the manifest was rebuilt and WAVES.md has a line in the same commit. The gate extension is correct and tested. I sampled 22 generators: B21, B39-41, B48, B51 (fit, ship, 51.1), B58 (fit, ship), B59 (holdout, 59a fit and diagnosis, 59b fit and anchor), B60 (fit, verification), B61-B63 (ships), B64 and route calibration. All set the same values at the same point. Import-time setters stay at import time, and restores keep the saved objects.
